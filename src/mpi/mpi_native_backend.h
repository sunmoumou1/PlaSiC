#ifndef PLASIC_MPI_NATIVE_BACKEND_H
#define PLASIC_MPI_NATIVE_BACKEND_H

/*
 * 原生 MPI 后端的声明：下面的函数由 mpi_native_backend.c 用真正的 MPI 调用
 * 实现，这两个文件是全工程唯一直接使用 <mpi.h> 的地方。
 *
 * 状态码约定：成功返回 0（即 MPI_SUCCESS），失败返回原始 MPI 错误码；
 * 参数不合法时返回 MPI_ERR_ARG。
 *
 * 本头文件只应由 plasic_mpi.c 在 PLASIC_USE_MPI=1 时包含。
 */
/*
 * Declarations of the native MPI backend: the functions below are implemented
 * in mpi_native_backend.c with real MPI calls, and these two files are the
 * only place in the project that uses <mpi.h> directly.
 *
 * Status convention: zero (MPI_SUCCESS) on success, the raw MPI error code on
 * failure, and MPI_ERR_ARG for invalid arguments.
 *
 * This header should be included only by plasic_mpi.c when PLASIC_USE_MPI=1.
 */

#include "plasic_mpi.h"

#include <stdint.h>

/*
 * 启动 MPI 并返回进程总数与编号；期望进程数不匹配时返回
 * PLASIC_MP_START_SIZE_MISMATCH，MPI 自身出错时返回 PLASIC_MP_START_ERROR。
 */
/*
 * Start MPI and report the total process count and rank; a count mismatch
 * yields PLASIC_MP_START_SIZE_MISMATCH and an MPI failure yields
 * PLASIC_MP_START_ERROR.
 */
int32_t mpi_native_start(
    int32_t expected_processes, int32_t *processes, int32_t *rank,
    int32_t *mpi_error);

/* Barrier 同步后调用 MPI_Finalize 关闭 MPI 环境。 */
/* Barrier synchronization followed by MPI_Finalize to tear down MPI. */
int32_t mpi_native_stop(void);

/* 由根进程（rank 0）广播 count 个 int32_t 元素（count 为元素个数）。 */
/* Broadcast count int32_t elements from the root (rank 0); count is an element
 * count. */
int32_t mpi_native_broadcast_integer(int32_t *buffer, int32_t count);

/* 由根进程广播 count 个 float 元素。 */
/* Broadcast count float elements from the root. */
int32_t mpi_native_broadcast_real(float *buffer, int32_t count);

/* 全收集格点场到所有进程：MPI_Allgather 的封装。 */
/* Allgather a grid field to every process: a wrapper around MPI_Allgather. */
int32_t mpi_native_allgather_grid(
    float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels);

/* 收集全局谱到根进程：MPI_Gather 的封装。 */
/* Gather the global spectrum onto the root: a wrapper around MPI_Gather. */
int32_t mpi_native_gather_spectral(
    float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels);

/* 全收集全局谱到所有进程：MPI_Allgather 的封装。 */
/* Allgather the global spectrum to every process: a wrapper around
 * MPI_Allgather. */
int32_t mpi_native_allgather_spectral(
    float *global, const float *local, int32_t global_rows,
    int32_t local_rows, int32_t levels);

/* 对 count 个 float 做逐元素 MPI_SUM 归约，结果发给所有进程。 */
/* Element-wise MPI_SUM reduction over count floats, with the result sent to
 * every process. */
int32_t mpi_native_allreduce_real(float *buffer, int32_t count);
int32_t mpi_native_sync_error(
    int32_t *status, char *message, int32_t bytes, int32_t rank);

/*
 * 用两次 MPI_Sendrecv 与相邻 rank 交换纬度 halo 行；
 * 详见 plasic_mpi.h 中 mp_exchange_latitude_halo 的说明。
 */
/*
 * Exchange latitude halo rows with the adjacent ranks through two
 * MPI_Sendrecv calls; see mp_exchange_latitude_halo in plasic_mpi.h.
 */
int32_t mpi_native_exchange_latitude_halo(
    float *recv_north, float *recv_south, const float *send_north,
    const float *send_south, int32_t elements);

#endif
