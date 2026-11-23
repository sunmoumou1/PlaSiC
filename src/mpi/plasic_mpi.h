#ifndef PLASIC_MPI_H
#define PLASIC_MPI_H

#include <stdint.h>

/*
 * PlaSiC 并行层的统一门面：整个模式只通过这里的 mp_* 函数使用并行能力；
 * 其他模块（runtime/、physics/、math_tools/、tools/）只包含本头文件，
 * 绝不直接包含 <mpi.h>。
 *
 * 具体后端由 plasic_mpi.c 中的编译期开关 PLASIC_USE_MPI 选择：
 * PLASIC_USE_MPI=1 使用 mpi_native_backend（真正的 MPI 消息传递）；
 * PLASIC_USE_MPI=0 使用 mpi_serial_backend（单进程、语义等价的退化实现）。
 *
 * 本头文件只依赖 <stdint.h>，因此串行构建完全不需要 MPI。
 */
/*
 * Unified facade of the PlaSiC parallel layer: the whole model uses
 * parallelism only through the mp_* functions declared here; every other
 * module (runtime/, physics/, math_tools/, tools/) includes this header and
 * never includes <mpi.h> directly.
 *
 * The concrete backend is chosen at compile time in plasic_mpi.c through
 * PLASIC_USE_MPI: 1 selects mpi_native_backend (real MPI message passing),
 * 0 selects mpi_serial_backend (a single-process, semantically equivalent
 * degradation).
 *
 * This header depends only on <stdint.h>, so serial builds need no MPI at all.
 */

/*
 * mp_start 的启动状态码；集合通信成功码统一为 0。
 */
/*
 * Startup status codes for mp_start; the success code of every collective
 * operation is 0 as well.
 */
enum plasic_mp_start_status
{
    /* 启动成功。 */
    /* Startup succeeded. */
    PLASIC_MP_START_OK = 0,
    /* 实际进程数与编译期常量 PLASIC_NPRO 不一致。 */
    /* The actual process count differs from the compile-time PLASIC_NPRO. */
    PLASIC_MP_START_SIZE_MISMATCH = 1,
    /* MPI 环境本身初始化失败。 */
    /* The MPI environment itself failed to initialize. */
    PLASIC_MP_START_ERROR = 2
};

/*
 * 生命周期：初始化所选后端，确定进程总数与本进程编号。
 * expected_processes 是编译期常量 PLASIC_NPRO；实际进程数写入 *processes；
 * 本进程号缓存在内部并由 mp_rank() 读出；后端原始状态码写入 *mpi_error。
 * 进程数不匹配时返回 PLASIC_MP_START_SIZE_MISMATCH。
 */
/*
 * Lifecycle: initialize the selected backend and determine the total process
 * count and this process's rank. expected_processes is the compile-time
 * constant PLASIC_NPRO; the actual count is written to *processes, the rank
 * is cached internally for mp_rank(), and the raw backend status is written
 * to *mpi_error. A count mismatch yields PLASIC_MP_START_SIZE_MISMATCH.
 */
int mp_start(
    int32_t expected_processes, int32_t *processes, int32_t *mpi_error);

/*
 * 关闭所选后端并复位内部状态。
 */
/*
 * Shut down the selected backend and reset the internal state.
 */
int32_t mp_stop(void);

/*
 * 身份与布局：以下四个函数不通信，只读取启动时缓存的状态。
 */
/*
 * Identity and layout: the four functions below communicate nothing; they
 * only read state cached at startup.
 */

/* 本进程在通信域中的编号（从 0 开始）。 */
/* This process's rank inside the communicator (starting from 0). */
int32_t mp_rank(void);

/* 通信域中的进程总数。 */
/* Total number of processes in the communicator. */
int32_t mp_size(void);

/* 本进程是否为根进程（rank 0）。 */
/* Whether this process is the root process (rank 0). */
int32_t mp_is_root(void);

/*
 * 本进程的数据块在全局数组中的起始偏移，即 rank * local_rows。
 */
/*
 * Starting offset of this process's data block inside the global array, that
 * is, rank * local_rows.
 */
int32_t mp_local_offset(int32_t local_rows);

/*
 * 集合通信：成功时返回 0，失败时返回非零状态码。
 * 注意：count 始终是元素个数，不是字节数。
 */
/*
 * Collective operations: zero on success and a non-zero status on failure.
 * Note that count is always an element count, never a byte count.
 */

/*
 * 广播：根进程把 buffer 中的 count 个元素发送给通信域内所有进程。
 */
/*
 * Broadcast: the root process sends the count elements of buffer to every
 * process in the communicator.
 */
int32_t mp_broadcast_integer(int32_t *buffer, int32_t count);

/* 广播 float 数组的版本。 */
/* Broadcast version for a float array. */
int32_t mp_broadcast_real(float *buffer, int32_t count);

/*
 * 全局归约：对 buffer 中的 count 个元素逐元素求和，并把结果发回每个进程；
 * 调用返回后所有 rank 持有完全相同的全局和。
 */
/*
 * Global reduction: sum the count elements of buffer element by element and
 * send the result back to every process; after the call every rank holds
 * exactly the same sums.
 */
int32_t mp_allreduce_real(float *buffer, int32_t count);

/* Share the status and message of the lowest rank reporting a nonzero status.
 * Returns a backend communication error independently of the shared status. */
int32_t mp_sync_error(int32_t *status, char *message, int32_t bytes);

/*
 * 全收集格点场：每个进程贡献自己在每层中的 local_rows 个格点，MPI 按 rank
 * 顺序拼成 global_rows 长的全局层并发给所有进程。
 * global 缓冲区至少容纳 global_rows * levels 个 float。
 */
/*
 * Allgather of a grid field: every process contributes its local_rows grid
 * points within each level, and MPI concatenates them in rank order into the
 * global_rows-long global level sent to every process. The global buffer must
 * hold at least global_rows * levels floats.
 */
int32_t mp_allgather_grid(
    float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels);

/*
 * 谱收集：数据布局与 mp_allgather_grid 相同，但完整结果只写入根进程的
 * global 缓冲区。
 *
 * 与格点分块不同，谱分块允许最后一个 rank 的局部块越过 global_rows 尾部：
 * 当 PLASIC_NPRO 不整除 PLASIC_NRSP 时，每个 rank 仍预留 local_rows 个槽位，
 * 但最后一个 rank 只有 clamp(global_rows - rank * local_rows, 0, local_rows)
 * 个槽位有效。后端负责只通信这些有效槽位（native 后端使用 MPI_Gatherv）。
 */
/*
 * Spectral gather: same data layout as mp_allgather_grid, but the complete
 * result is written only into the root process's global buffer.
 *
 * Unlike the grid decomposition, the spectral block of the last rank may
 * extend past the end of global_rows: when PLASIC_NPRO does not divide
 * PLASIC_NRSP, every rank still reserves local_rows slots, but the last rank
 * owns only clamp(global_rows - rank * local_rows, 0, local_rows) valid ones.
 * The backend communicates only those valid slots (the native backend uses
 * MPI_Gatherv).
 */
int32_t mp_gather_spectral(
    float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels);

/*
 * 谱全收集：每个进程都得到完整的全局谱。局部块的语义与 mp_gather_spectral
 * 相同，native 后端使用 MPI_Allgatherv。
 */
/*
 * Spectral allgather: every process receives the complete global spectrum.
 * The local-block semantics match mp_gather_spectral, and the native backend
 * uses MPI_Allgatherv.
 */
int32_t mp_allgather_spectral(
    float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels);

/*
 * 与相邻 rank 交换一行（elements 个 float）纬度边界数据。
 * 全局纬度索引自北向南增大，因此 rank r 的北侧邻居是 r-1，南侧邻居是 r+1。
 *
 *   send_north：本进程北侧边缘行，发往 rank r-1
 *   send_south：本进程南侧边缘行，发往 rank r+1
 *   recv_north：接收来自 rank r-1 的行
 *   recv_south：接收来自 rank r+1 的行
 *
 * 边界进程在某一方向没有邻居；对应的接收缓冲区保持不变，以便调用者预填
 * 镜像（零梯度）边界值。
 */
/*
 * Exchange one latitude-boundary row of `elements` floats with the adjacent
 * ranks. Global latitude indices increase southward, so rank r's northern
 * neighbour is rank r-1 and its southern neighbour is rank r+1.
 *
 *   send_north: this rank's northern edge row, sent toward rank r-1
 *   send_south: this rank's southern edge row, sent toward rank r+1
 *   recv_north: receives the row from rank r-1
 *   recv_south: receives the row from rank r+1
 *
 * Boundary ranks have no neighbour in one direction; the corresponding
 * receive buffer is left untouched so the caller can pre-fill it with a
 * mirrored (zero-gradient) edge.
 */
int32_t mp_exchange_latitude_halo(
    float *recv_north, float *recv_south, const float *send_north,
    const float *send_south, int32_t elements);

#endif
