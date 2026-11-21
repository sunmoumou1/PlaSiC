/*
 * PlaSiC 的原生 MPI 后端：所有集合操作最终落到 MPI_Bcast、MPI_Allreduce、
 * MPI_Allgather、MPI_Gather 或 MPI_Sendrecv，是全工程唯一直接使用 <mpi.h>
 * 的实现文件。
 *
 * 完整的 MPI 背景与逐函数讲解见
 * docs-site/docs/7-appendix/7.7-mpi-parallel-basics.md。
 */
/*
 * PlaSiC's native MPI backend: every collective operation ends up in MPI_Bcast,
 * MPI_Allreduce, MPI_Allgather, MPI_Gather, or MPI_Sendrecv; this is the only
 * implementation file in the project that uses <mpi.h> directly.
 *
 * See docs-site/docs/7-appendix/7.7-mpi-parallel-basics.md for the full background
 * and a function-by-function walkthrough.
 */

#include "mpi_native_backend.h"

#include <mpi.h>

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/* 模型数据以 rank 0 为根；错误消息由报告该错误的进程广播。 */
/* Model data use rank 0 as root. Error messages come from the reporting rank. */
enum
{
    ROOT_RANK = 0
};

/*
 * 启动时保存的通信域（MPI_COMM_WORLD 的副本）；所有通信都在它上面进行，
 * 关闭后复位为 MPI_COMM_NULL。
 */
/*
 * Communicator saved at startup (a copy of MPI_COMM_WORLD); all communication
 * takes place on it, and it is reset to MPI_COMM_NULL on shutdown.
 */
static MPI_Comm active_world = MPI_COMM_NULL;

/* 元素个数必须非负。 */
/* The element count must be non-negative. */
static int valid_count(int32_t count)
{
    return count >= 0;
}

/*
 * 启动 MPI 环境并回答“我是谁、一共有多少人”。
 * 依次调用 MPI_Init、MPI_Comm_size、MPI_Comm_rank，把结果写入调用者的输出
 * 参数；实际进程数与 expected_processes 不一致时返回
 * PLASIC_MP_START_SIZE_MISMATCH（此时 MPI 环境已经建立，不匹配只是配置错误）；
 * MPI 自身失败时返回 PLASIC_MP_START_ERROR。
 */
/*
 * Start the MPI environment and answer "who am I and how many of us are
 * there". MPI_Init, MPI_Comm_size, and MPI_Comm_rank are called in turn and
 * their results are written to the caller's output parameters. A mismatch
 * between the actual process count and expected_processes yields
 * PLASIC_MP_START_SIZE_MISMATCH (the MPI environment is already up; only the
 * configuration is wrong); an MPI failure yields PLASIC_MP_START_ERROR.
 */
int32_t mpi_native_start(
    int32_t expected_processes, int32_t *processes, int32_t *rank,
    int32_t *mpi_error)
{
    int rank_value = 0;
    int size_value = 1;
    int status = MPI_SUCCESS;

    if (processes == NULL || rank == NULL || mpi_error == NULL)
    {
        return PLASIC_MP_START_ERROR;
    }
    status = MPI_Init(NULL, NULL);
    /* MPI_COMM_WORLD 在 MPI_Init 之后才有效，保存下来供后续所有通信用。 */
    /* MPI_COMM_WORLD becomes valid only after MPI_Init; save it for all later
     * communication. */
    active_world = MPI_COMM_WORLD;
    if (status == MPI_SUCCESS)
    {
        status = MPI_Comm_size(active_world, &size_value);
    }
    if (status == MPI_SUCCESS)
    {
        status = MPI_Comm_rank(active_world, &rank_value);
    }
    /* 无论成功与否都把已有信息交回调用者，*mpi_error 保留原始状态码。 */
    /* Hand back whatever information was obtained; *mpi_error keeps the raw
     * status code. */
    *mpi_error = (int32_t)status;
    *processes = (int32_t)size_value;
    *rank = (int32_t)rank_value;
    if (status != MPI_SUCCESS)
    {
        return PLASIC_MP_START_ERROR;
    }
    if (size_value != expected_processes)
    {
        return PLASIC_MP_START_SIZE_MISMATCH;
    }
    return PLASIC_MP_START_OK;
}

/*
 * 先做一次 Barrier 再 MPI_Finalize：Barrier 是集合同步点，保证所有进程都到达
 * 之后才一起关闭，避免某个进程提前 Finalize 时其他进程仍在通信。
 */
/*
 * Barrier first, then MPI_Finalize: the barrier is a collective synchronization
 * point that guarantees every process has arrived before the environment is
 * torn down, so no process finalizes while others are still communicating.
 */
int32_t mpi_native_stop(void)
{
    int status = MPI_Barrier(active_world);
    if (status == MPI_SUCCESS)
    {
        status = MPI_Finalize();
    }
    active_world = MPI_COMM_NULL;
    return (int32_t)status;
}

/*
 * 把根进程（rank 0）的 count 个 int32_t 元素广播给通信域内所有进程。
 * count 是元素个数（不是字节数）。
 */
/*
 * Broadcast the count int32_t elements of the root (rank 0) to every process
 * in the communicator. count is an element count, not a byte count.
 */
int32_t mpi_native_broadcast_integer(int32_t *buffer, int32_t count)
{
    if (buffer == NULL || !valid_count(count))
    {
        return (int32_t)MPI_ERR_ARG;
    }
    return (int32_t)MPI_Bcast(buffer, (int)count, MPI_INT32_T,
                              ROOT_RANK, active_world);
}

/* 与上面相同，广播的是 count 个 float 元素。 */
/* Same as above, broadcasting count float elements. */
int32_t mpi_native_broadcast_real(float *buffer, int32_t count)
{
    if (buffer == NULL || !valid_count(count))
    {
        return (int32_t)MPI_ERR_ARG;
    }
    return (int32_t)MPI_Bcast(buffer, (int)count, MPI_FLOAT,
                              ROOT_RANK, active_world);
}

int32_t mpi_native_sync_error(
    int32_t *status, char *message, int32_t bytes, int32_t rank)
{
    int32_t reporter = *status == 0 ? INT32_MAX : rank;
    int rc = MPI_Allreduce(MPI_IN_PLACE, &reporter, 1, MPI_INT32_T,
                           MPI_MIN, active_world);
    if (rc != MPI_SUCCESS || reporter == INT32_MAX)
    {
        return (int32_t)rc;
    }
    rc = MPI_Bcast(status, 1, MPI_INT32_T, reporter, active_world);
    if (rc == MPI_SUCCESS)
    {
        rc = MPI_Bcast(message, bytes, MPI_CHAR, reporter, active_world);
    }
    return (int32_t)rc;
}

/*
 * 按层收集格点场的公共实现。
 * 每一层中，每个进程贡献自己连续的 local_rows 个 float，MPI 按 rank 顺序把它们
 * 拼进 global_rows 长的全局层：gather_to_all 非零时用 MPI_Allgather（每个进程
 * 都得到完整结果），为零时用 MPI_Gather（只有根进程写 global）。
 *
 * 之所以能这样直接拼接，是因为格点分块约定保证每个进程的块大小相同、按 rank
 * 顺序连续排列（MPI 构建要求 PLASIC_NPRO 整除 PLASIC_NLAT，见
 * state/plasicmod_state.h），因此总元素数恰好是 size * local_rows。
 *
 * 注意：recvcount 传的是 local_rows，即“来自每个进程”的元素数，而不是接收总数。
 */
/*
 * Shared per-level gather implementation for grid fields. Within each level
 * every process contributes its local_rows consecutive floats, and MPI
 * concatenates them in rank order into the global_rows-long global level:
 * gather_to_all non-zero selects MPI_Allgather (everyone receives the full
 * result), zero selects MPI_Gather (only the root's global buffer is written).
 *
 * This concatenation works only because the grid decomposition convention
 * gives every process an equal-sized block laid out contiguously in rank order
 * (MPI builds require PLASIC_NPRO to divide PLASIC_NLAT; see
 * state/plasicmod_state.h), so the total is exactly size * local_rows.
 *
 * Note: recvcount is local_rows, the number of elements received from each
 * process, not the total number received.
 */
static int gather_columns(float *global, const float *local,
                          int32_t global_rows, int32_t local_rows,
                          int32_t levels, int gather_to_all)
{
    int32_t level;
    int status = MPI_SUCCESS;

    if (global == NULL || local == NULL || global_rows < 0 || local_rows < 0 ||
        levels < 0)
    {
        return MPI_ERR_ARG;
    }
    /* 逐层循环；一旦某层失败立即停止后续通信。 */
    /* Loop level by level; stop further communication as soon as one level
     * fails. */
    for (level = 0; level < levels && status == MPI_SUCCESS; ++level)
    {
        if (gather_to_all)
        {
            /* 发送与接收指针都按层偏移，接收侧每层跨 global_rows。 */
            /* Both the send and receive pointers are offset by level; the
             * receive side strides global_rows per level. */
            status = MPI_Allgather(local + (size_t)level * (size_t)local_rows,
                                   (int)local_rows, MPI_FLOAT,
                                   global + (size_t)level * (size_t)global_rows,
                                   (int)local_rows, MPI_FLOAT, active_world);
        }
        else
        {
            /* MPI_Gather 是 MPI_Allgather 的“只给根进程”版本。 */
            /* MPI_Gather is the root-only version of MPI_Allgather. */
            status = MPI_Gather(local + (size_t)level * (size_t)local_rows,
                                (int)local_rows, MPI_FLOAT,
                                global + (size_t)level * (size_t)global_rows,
                                (int)local_rows, MPI_FLOAT,
                                ROOT_RANK, active_world);
        }
    }
    return status;
}

/*
 * 谱分解中 rank 实际贡献/接收的有效元素个数。
 *
 * 谱块长度 local_rows = NSPP 是按 ceil(NRSP / NPRO) 预留的，因此当 NPRO
 * 不整除 NRSP 时，最后一个 rank 只剩一个部分块：
 *
 *     valid(rank) = clamp(global_rows - rank * local_rows, 0, local_rows)
 *
 * 前 NPRO-1 个 rank 仍得到完整的 local_rows；padding 槽位不参与通信。
 *
 * Number of valid elements contributed/received by a rank in the spectral
 * decomposition. The block length local_rows = NSPP is reserved as
 * ceil(NRSP / NPRO), so when NPRO does not divide NRSP the last rank owns
 * only a partial block:
 *
 *     valid(rank) = clamp(global_rows - rank * local_rows, 0, local_rows)
 *
 * The first NPRO-1 ranks still get a full local_rows block; padding slots
 * never take part in communication.
 */
static int32_t spectral_block_count(
    int32_t rank_value, int32_t global_rows, int32_t local_rows)
{
    const int64_t start = (int64_t)rank_value * (int64_t)local_rows;
    int64_t valid = (int64_t)global_rows - start;

    if (valid > (int64_t)local_rows)
    {
        valid = (int64_t)local_rows;
    }
    if (valid < 0)
    {
        valid = 0;
    }
    return (int32_t)valid;
}

/*
 * 按层收集谱场的公共实现（mp_gather_spectral / mp_allgather_spectral）。
 *
 * 与格点场的 gather_columns 不同，谱分块允许最后一个 rank 的表观块越过
 * global_rows 的尾部（padding 槽位），因此必须使用可变计数的
 * MPI_Allgatherv / MPI_Gatherv：
 *
 *     counts[r]        = min(local_rows, max(global_rows - r * local_rows, 0))
 *     displacements[r] = r * local_rows
 *
 * 接收侧每层仍从 global + level * global_rows 开始，按 rank 顺序恰好写满
 * global_rows 个有效元素，于是全局谱布局与旧实现完全一致。只要每个 rank
 * 的计数相等（NPRO 整除 NRSP 时 counts 全为 local_rows），Allgatherv 的
 * 结果与原来的 Allgather 逐位相同。
 *
 * Shared per-level gather implementation for spectral fields
 * (mp_gather_spectral / mp_allgather_spectral).
 *
 * Unlike the grid gather_columns, the spectral block of the last rank may
 * extend past global_rows (its padding slots), so variable-count
 * MPI_Allgatherv / MPI_Gatherv is required:
 *
 *     counts[r]        = min(local_rows, max(global_rows - r * local_rows, 0))
 *     displacements[r] = r * local_rows
 *
 * The receive side still starts each level at global + level * global_rows
 * and fills exactly global_rows valid elements in rank order, so the global
 * spectral layout is identical to the previous implementation. Whenever every
 * rank has the same count (all local_rows when NPRO divides NRSP), the
 * Allgatherv result is bit-for-bit identical to the old Allgather.
 */
static int gather_spectral_columns(
    float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels, int gather_to_all)
{
    int rank_value = 0;
    int size_value = 0;
    int *counts = NULL;
    int *displacements = NULL;
    int32_t level;
    int32_t processor;
    int status = MPI_SUCCESS;

    if (global == NULL || local == NULL || global_rows < 0 || local_rows < 0 ||
        levels < 0)
    {
        return MPI_ERR_ARG;
    }
    if (MPI_Comm_rank(active_world, &rank_value) != MPI_SUCCESS ||
        MPI_Comm_size(active_world, &size_value) != MPI_SUCCESS)
    {
        return MPI_ERR_OTHER;
    }
    if (size_value > 0)
    {
        counts = (int *)malloc((size_t)size_value * sizeof(*counts));
        displacements =
            (int *)malloc((size_t)size_value * sizeof(*displacements));
        if (counts == NULL || displacements == NULL)
        {
            free(counts);
            free(displacements);
            return MPI_ERR_NO_MEM;
        }
        for (processor = 0; processor < size_value; ++processor)
        {
            counts[processor] = (int)spectral_block_count(
                processor, global_rows, local_rows);
            displacements[processor] =
                (int)((int64_t)processor * (int64_t)local_rows);
        }
    }

    /* 逐层循环；一旦某层失败立即停止后续通信。 */
    /* Loop level by level; stop further communication as soon as one level
     * fails. */
    for (level = 0; level < levels && status == MPI_SUCCESS; ++level)
    {
        const float *send_buffer =
            local + (size_t)level * (size_t)local_rows;
        float *receive_buffer =
            global + (size_t)level * (size_t)global_rows;
        const int send_count = counts != NULL ? counts[rank_value] : 0;

        if (gather_to_all)
        {
            status = MPI_Allgatherv(
                send_buffer, send_count, MPI_FLOAT,
                receive_buffer, counts, displacements, MPI_FLOAT,
                active_world);
        }
        else
        {
            status = MPI_Gatherv(
                send_buffer, send_count, MPI_FLOAT,
                receive_buffer, counts, displacements, MPI_FLOAT,
                ROOT_RANK, active_world);
        }
    }
    free(counts);
    free(displacements);
    return status;
}

/*
 * 全收集格点场：把各进程的局部纬度带拼成完整全局场，并发给每个进程。
 */
/*
 * Allgather a grid field: assemble the local latitude bands into the complete
 * global field and deliver it to every process.
 */
int32_t mpi_native_allgather_grid(
    float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels)
{
    return (int32_t)gather_columns(global, local, global_rows, local_rows,
                                   levels, 1);
}

/*
 * 谱收集：与 allgather 相同，但完整结果只写入根进程的 global。
 */
/*
 * Spectral gather: same as the allgather variant, but the complete result is
 * written only into the root process's global buffer.
 */
int32_t mpi_native_gather_spectral(
    float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels)
{
    return (int32_t)gather_spectral_columns(
        global, local, global_rows, local_rows, levels, 0);
}

/* 谱全收集：每个进程都得到完整的全局谱。 */
/* Spectral allgather: every process receives the complete global spectrum. */
int32_t mpi_native_allgather_spectral(
    float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels)
{
    return (int32_t)gather_spectral_columns(
        global, local, global_rows, local_rows, levels, 1);
}

/*
 * 对 count 个 float 做 MPI_SUM 的逐元素全局归约，结果发回所有进程。
 * MPI_Allreduce 要求 sendbuf 与 recvbuf 不是同一块内存，因此先分配临时缓冲接收
 * 结果，成功后再拷回 buffer。
 */
/*
 * Element-wise MPI_SUM global reduction over count floats, with the result sent
 * back to every process. MPI_Allreduce requires sendbuf and recvbuf to be
 * distinct, so a temporary receive buffer is allocated and the result is copied
 * back into buffer on success.
 */
int32_t mpi_native_allreduce_real(float *buffer, int32_t count)
{
    float *temporary;
    int status;

    if (buffer == NULL || !valid_count(count))
    {
        return MPI_ERR_ARG;
    }
    /* count == 0 时仍需要一个合法指针，因此至少分配 1 个元素。 */
    /* count == 0 still needs a valid pointer, so allocate at least one
     * element. */
    temporary = (float *)malloc((size_t)(count == 0 ? 1 : count) *
                                sizeof(*temporary));
    if (temporary == NULL)
    {
        return MPI_ERR_NO_MEM;
    }
    status = MPI_Allreduce(buffer, temporary, (int)count, MPI_FLOAT,
                           MPI_SUM, active_world);
    /* 归约成功且确有元素时才把结果拷回。 */
    /* Copy the result back only when the reduction succeeded and there was
     * data. */
    if (status == MPI_SUCCESS && count != 0)
    {
        memcpy(buffer, temporary, (size_t)count * sizeof(*buffer));
    }
    free(temporary);
    return (int32_t)status;
}

/*
 * 纬度 halo 交换：src/mpi/ 中唯一使用点对点通信的函数，服务于三维海洋的南北
 * 差分。
 *
 * 全局纬度索引自北向南增大，所以北侧邻居是 rank-1、南侧邻居是 rank+1。边界进程
 * 在缺失的一侧使用 MPI_PROC_NULL：发送退化为空操作，接收立即返回且不修改接收
 * 缓冲区，从而保留调用者预填的镜像（零梯度）边界值。
 *
 * 两次 MPI_Sendrecv 把“发送”和“接收”合成一个原子操作，两个方向天然互补，
 * 因此不会出现互相等待的死锁；两个阶段使用不同 tag，避免消息混淆（MPI 以
 * “通信域 + 来源 + tag”唯一标识一条消息）。
 */
/*
 * Latitude halo exchange: the only function in src/mpi/ that uses point-to-point
 * communication; it serves the north-south differences of the 3-D ocean.
 *
 * Global latitude indices increase from north to south, so the northern
 * neighbour is rank-1 and the southern neighbour is rank+1. At a boundary the
 * missing side uses MPI_PROC_NULL: the send becomes a no-op and the receive
 * returns immediately without touching the receive buffer, preserving the
 * caller's pre-filled mirrored (zero-gradient) edge values.
 *
 * The two MPI_Sendrecv calls fuse send and receive into one atomic operation;
 * the directions are naturally complementary, so the mutual-wait deadlock
 * cannot occur. The two stages use different tags so their messages cannot be
 * confused (MPI identifies a message uniquely by communicator, source, and
 * tag).
 */
int32_t mpi_native_exchange_latitude_halo(
    float *recv_north, float *recv_south, const float *send_north,
    const float *send_south, int32_t elements)
{
    enum
    {
        HALO_NORTH_TAG = 41,
        HALO_SOUTH_TAG = 42
    };
    int rank_value = 0;
    int size_value = 1;
    int north_neighbor;
    int south_neighbor;
    MPI_Status status;

    if (recv_north == NULL || recv_south == NULL ||
        send_north == NULL || send_south == NULL || !valid_count(elements))
    {
        return (int32_t)MPI_ERR_ARG;
    }
    if (MPI_Comm_rank(active_world, &rank_value) != MPI_SUCCESS ||
        MPI_Comm_size(active_world, &size_value) != MPI_SUCCESS)
    {
        return (int32_t)MPI_ERR_OTHER;
    }
    /* 北侧邻居 rank-1，南侧邻居 rank+1；边界一侧没有邻居。 */
    /* Northern neighbour rank-1, southern neighbour rank+1; a boundary rank
     * has no neighbour on one side. */
    north_neighbor = rank_value > 0 ? rank_value - 1 : MPI_PROC_NULL;
    south_neighbor = rank_value < size_value - 1 ? rank_value + 1
                                                 : MPI_PROC_NULL;

    /*
     * 第一次交换（tag = HALO_NORTH_TAG）：向北发送，从南接收
     *   send_north → north_neighbor (rank-1)，成为对方的南侧 halo
     *   recv_south ← south_neighbor (rank+1)，接收对方的北侧边缘行
     */
    /*
     * First exchange (tag = HALO_NORTH_TAG): send north, receive from south
     *   send_north → north_neighbor (rank-1), becoming its southern halo
     *   recv_south ← south_neighbor (rank+1), receiving its northern edge row
     */
    if (MPI_Sendrecv(send_north, (int)elements, MPI_FLOAT, north_neighbor,
                     HALO_NORTH_TAG, recv_south, (int)elements, MPI_FLOAT,
                     south_neighbor, HALO_NORTH_TAG, active_world,
                     &status) != MPI_SUCCESS)
    {
        return (int32_t)MPI_ERR_OTHER;
    }
    /*
     * 第二次交换（tag = HALO_SOUTH_TAG）：向南发送，从北接收，与第一次互为镜像
     *   send_south → south_neighbor (rank+1)，成为对方的北侧 halo
     *   recv_north ← north_neighbor (rank-1)，接收对方的南侧边缘行
     */
    /*
     * Second exchange (tag = HALO_SOUTH_TAG): send south, receive from north,
     * the mirror image of the first exchange
     *   send_south → south_neighbor (rank+1), becoming its northern halo
     *   recv_north ← north_neighbor (rank-1), receiving its southern edge row
     */
    return (int32_t)MPI_Sendrecv(
        send_south, (int)elements, MPI_FLOAT, south_neighbor,
        HALO_SOUTH_TAG, recv_north, (int)elements, MPI_FLOAT, north_neighbor,
        HALO_SOUTH_TAG, active_world, &status);
}
