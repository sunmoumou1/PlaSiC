/*
 * PlaSiC 并行接口的分发层（门面的实现）：根据编译期开关 PLASIC_USE_MPI 把每个
 * mp_* 调用转发给 native 或 serial 后端，并在这里统一做参数校验。
 *
 * 注意：这不是运行期变量而是编译期常量。src/Makefile 通过
 * -DPLASIC_USE_MPI=$(MPI) 注入它，同一份源码编译出的两个可执行文件因此走两条
 * 完全不同的后端代码路径。
 */
/*
 * Dispatcher of the PlaSiC parallel interface (the implementation behind the
 * facade): forwards every mp_* call to the native or serial backend according
 * to the compile-time switch PLASIC_USE_MPI, and performs the shared argument
 * validation here.
 *
 * Note that this is a compile-time constant rather than a runtime variable.
 * src/Makefile injects it with -DPLASIC_USE_MPI=$(MPI), so the two executables
 * built from the same sources follow two completely different backend code
 * paths.
 */

#include "plasic_mpi.h"
#include "mpi_serial_backend.h"

#ifndef PLASIC_USE_MPI
/* 编译命令行未指定时默认选择串行后端。 */
/* Defaults to the serial backend when the compiler command line omits it. */
#define PLASIC_USE_MPI 0
#endif

#if PLASIC_USE_MPI
/* 只有 MPI 构建才包含 native 后端；这里是它唯一的包含点。 */
/* Only MPI builds include the native backend; this is its only inclusion
 * point. */
#include "mpi_native_backend.h"
#endif

#include <stddef.h>

/*
 * 启动时缓存的身份状态，mp_rank()/mp_size()/mp_is_root() 直接读取它们。
 */
/*
 * Identity state cached at startup and read directly by
 * mp_rank()/mp_size()/mp_is_root().
 */
static int32_t active_rank;
static int32_t active_processes = 1;

/*
 * 非零表示后端已经启动；mp_stop() 据此决定是否执行后端停止流程。
 */
/*
 * Non-zero once the backend has started; mp_stop() uses it to decide whether
 * the backend shutdown must be performed.
 */
static int backend_started;

/*
 * 集合操作的公共参数校验：指针非空、行数与层数非负。
 */
/*
 * Shared argument validation for collectives: non-NULL pointers and
 * non-negative row and level counts.
 */
static int valid_collective(
    const float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels)
{
    return global != NULL && local != NULL && global_rows >= 0 &&
           local_rows >= 0 && levels >= 0;
}

/*
 * 初始化所选后端；完整语义见 plasic_mpi.h 中 mp_start 的说明。
 */
/*
 * Initialize the selected backend; see mp_start in plasic_mpi.h for the full
 * semantics.
 */
int mp_start(
    int32_t expected_processes, int32_t *processes, int32_t *mpi_error)
{
    if (processes == NULL || mpi_error == NULL || expected_processes <= 0)
    {
        return PLASIC_MP_START_ERROR;
    }

    /* 先恢复单进程初值，保证后端启动失败时内部状态仍然自洽。 */
    /* Reset to the single-process defaults so the internal state stays
     * consistent even if the backend fails to start. */
    active_rank = 0;
    active_processes = 1;
    backend_started = 0;
    *mpi_error = 0;

#if PLASIC_USE_MPI
    {
        const int status = (int)mpi_native_start(
            expected_processes, &active_processes, &active_rank, mpi_error);

        /*
         * 进程数不匹配时 MPI 环境其实已经建立，因此仍标记为已启动，
         * 使随后的 mp_stop() 能正常调用 MPI_Finalize()。
         */
        /*
         * On a size mismatch the MPI environment is in fact already up, so it
         * is still marked as started and the later mp_stop() can call
         * MPI_Finalize() normally.
         */
        if (status == PLASIC_MP_START_OK ||
            status == PLASIC_MP_START_SIZE_MISMATCH)
        {
            backend_started = 1;
        }
        *processes = active_processes;
        return status;
    }
#else
    {
        const int status = (int)mpi_serial_start(
            expected_processes, &active_processes, &active_rank, mpi_error);

        *processes = active_processes;
        return status;
    }
#endif
}

/*
 * 关闭所选后端，然后把缓存的进程身份复位为单进程状态。
 */
/*
 * Shut down the selected backend and reset the cached process identity to the
 * single-process state.
 */
int32_t mp_stop(void)
{
    int32_t status = 0;

#if PLASIC_USE_MPI
    /* 只有真正启动过的 native 后端才需要关闭。 */
    /* Only a native backend that actually started needs to be shut down. */
    if (backend_started)
    {
        status = mpi_native_stop();
    }
#else
    status = mpi_serial_stop();
#endif
    backend_started = 0;
    active_rank = 0;
    active_processes = 1;
    return status;
}

/* 本进程编号。 */
/* This process's rank. */
int32_t mp_rank(void)
{
    return active_rank;
}

/* 进程总数。 */
/* Total number of processes. */
int32_t mp_size(void)
{
    return active_processes;
}

/* 是否为根进程（rank 0）。 */
/* Whether this process is the root (rank 0). */
int32_t mp_is_root(void)
{
    return active_rank == 0;
}

/* 本进程数据块在全局数组中的起始偏移。 */
/* Starting offset of this process's block in the global array. */
int32_t mp_local_offset(int32_t local_rows)
{
    return active_rank * local_rows;
}

/*
 * 广播 count 个 int32_t 元素；根进程是 rank 0。
 */
/*
 * Broadcast count int32_t elements; rank 0 is the root.
 */
int32_t mp_broadcast_integer(int32_t *buffer, int32_t count)
{
    if (buffer == NULL || count < 0)
    {
        return -1;
    }
#if PLASIC_USE_MPI
    return mpi_native_broadcast_integer(buffer, count);
#else
    return mpi_serial_broadcast_integer(buffer, count);
#endif
}

/* 广播 count 个 float 元素。 */
/* Broadcast count float elements. */
int32_t mp_broadcast_real(float *buffer, int32_t count)
{
    if (buffer == NULL || count < 0)
    {
        return -1;
    }
#if PLASIC_USE_MPI
    return mpi_native_broadcast_real(buffer, count);
#else
    return mpi_serial_broadcast_real(buffer, count);
#endif
}

/*
 * 对 count 个 float 逐元素求全局和（MPI_SUM），结果发回所有进程。
 */
/*
 * Element-wise global summation (MPI_SUM) over count floats, with the result
 * returned to every process.
 */
int32_t mp_allreduce_real(float *buffer, int32_t count)
{
    if (buffer == NULL || count < 0)
    {
        return -1;
    }
#if PLASIC_USE_MPI
    return mpi_native_allreduce_real(buffer, count);
#else
    return mpi_serial_allreduce_real(buffer, count);
#endif
}

int32_t mp_sync_error(int32_t *status, char *message, int32_t bytes)
{
    if (status == NULL || message == NULL || bytes < 1)
    {
        return -1;
    }
#if PLASIC_USE_MPI
    return mpi_native_sync_error(status, message, bytes, active_rank);
#else
    return 0;
#endif
}

/*
 * 全收集格点场：把各进程的局部纬度带按 rank 顺序拼成全局场并发给所有进程。
 */
/*
 * Allgather of a grid field: concatenate the local latitude bands of all ranks
 * in rank order and deliver the global field to every process.
 */
int32_t mp_allgather_grid(
    float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels)
{
    if (!valid_collective(
            global, local, global_rows, local_rows, levels))
    {
        return -1;
    }
#if PLASIC_USE_MPI
    return mpi_native_allgather_grid(
        global, local, global_rows, local_rows, levels);
#else
    return mpi_serial_allgather_grid(
        global, local, global_rows, local_rows, levels);
#endif
}

/*
 * 谱收集：完整结果只写入根进程的 global 缓冲区。
 */
/*
 * Spectral gather: the complete result is written only into the root's global
 * buffer.
 */
int32_t mp_gather_spectral(
    float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels)
{
    if (!valid_collective(
            global, local, global_rows, local_rows, levels))
    {
        return -1;
    }
#if PLASIC_USE_MPI
    return mpi_native_gather_spectral(
        global, local, global_rows, local_rows, levels);
#else
    return mpi_serial_gather_spectral(
        global, local, global_rows, local_rows, levels);
#endif
}

/*
 * 谱全收集：每个进程都得到完整的全局谱。
 */
/*
 * Spectral allgather: every process receives the complete global spectrum.
 */
int32_t mp_allgather_spectral(
    float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels)
{
    if (!valid_collective(
            global, local, global_rows, local_rows, levels))
    {
        return -1;
    }
#if PLASIC_USE_MPI
    return mpi_native_allgather_spectral(
        global, local, global_rows, local_rows, levels);
#else
    return mpi_serial_allgather_spectral(
        global, local, global_rows, local_rows, levels);
#endif
}

/*
 * 纬度 halo 交换：完整语义（邻居关系、边界处理）见 plasic_mpi.h 中的接口
 * 说明。
 */
/*
 * Latitude halo exchange: the full semantics (neighbour relations, boundary
 * handling) are documented with the interface declaration in plasic_mpi.h.
 */
int32_t mp_exchange_latitude_halo(
    float *recv_north, float *recv_south, const float *send_north,
    const float *send_south, int32_t elements)
{
    if (recv_north == NULL || recv_south == NULL ||
        send_north == NULL || send_south == NULL || elements < 0)
    {
        return -1;
    }
#if PLASIC_USE_MPI
    return mpi_native_exchange_latitude_halo(
        recv_north, recv_south, send_north, send_south, elements);
#else
    return mpi_serial_exchange_latitude_halo(
        recv_north, recv_south, send_north, send_south, elements);
#endif
}
