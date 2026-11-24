#include "runtime_internal.h"

static int require_restart_integer(const char *name, int32_t *value)
{
    const int status = restart_read_integer(name, value);

    if (status != PLASIC_RESTART_OK)
    {
        runtime_set_error("restart integer '%s' could not be read (status %d)", name, status);
        return PLASIC_RUNTIME_RESTART_FAILED;
    }

    return PLASIC_RUNTIME_OK;
}

/*
 * 读取 restart 中必须存在的整数记录；记录缺失时给出明确错误。
 * 目前只有 nstep 使用该函数：nstep 是模型唯一时钟，缺失即无法解释时间。
 */
static int require_restart_integer_mandatory(const char *name, int32_t *value)
{
    const int status = restart_read_integer(name, value);

    if (status == PLASIC_RESTART_NOT_FOUND)
    {
        runtime_set_error(
            "restart file has no mandatory '%s' record", name);
        return PLASIC_RUNTIME_RESTART_FAILED;
    }
    if (status != PLASIC_RESTART_OK)
    {
        runtime_set_error(
            "restart integer '%s' could not be read (status %d)", name,
            status);
        return PLASIC_RUNTIME_RESTART_FAILED;
    }

    return PLASIC_RUNTIME_OK;
}

static int require_restart_array(const char *name, float *values,
                                 int32_t rows, int32_t columns)
{
    const int status = restart_read_array(name, values, rows, columns);

    if (status != PLASIC_RESTART_OK)
    {
        runtime_set_error("restart array '%s' could not be read (status %d)", name, status);
        return PLASIC_RUNTIME_RESTART_FAILED;
    }

    return PLASIC_RUNTIME_OK;
}

static int require_restart_array_slice(
    const char *name, float *values, int32_t global_rows,
    int32_t local_offset, int32_t local_rows, int32_t columns)
{
    const int status = restart_read_array_slice(
        name, values, global_rows, local_offset, local_rows, columns);

    if (status != PLASIC_RESTART_OK)
    {
        runtime_set_error("restart array slice '%s' could not be read (status %d)", name, status);
        return PLASIC_RUNTIME_RESTART_FAILED;
    }

    return PLASIC_RUNTIME_OK;
}

/*
 * 从完整的全局谱数组中提取当前 rank 负责的局部谱块。
 *
 * 当 PLASIC_NPRO 不整除 PLASIC_NRSP 时，最后一个 rank 只有部分槽位有效：
 * 有效部分从全局数组复制，剩余槽位清零，作为既不参与通信也不参与物理解算
 * 的 padding。调用者必须保证 local 的每个 rank 块都按 PLASIC_NSPP 个槽位
 * 分配（见 state/plasicmod_state.h 中的 PLASIC_NESP）。
 *
 * Extract the local spectral block owned by the current MPI rank from the
 * complete global spectral array.
 *
 * When PLASIC_NPRO does not divide PLASIC_NRSP, the last rank owns only a
 * partial block: the valid part is copied from the global array and the
 * remaining slots are zeroed as padding that neither participates in
 * communication nor enters any physical computation. The caller must ensure
 * every rank-local block holds PLASIC_NSPP slots (see PLASIC_NESP in
 * state/plasicmod_state.h).
 */
static void split_global_spectral(
    const float *global, float *local, int32_t levels)
{
    const int32_t rank = mp_rank();
    const size_t rank_offset = (size_t)rank * (size_t)PLASIC_NSPP;
    const int32_t valid_count = plasic_local_spectral_count(rank);
    const size_t padding_count = (size_t)(PLASIC_NSPP - valid_count);
    int32_t level;

    for (level = 0; level < levels; ++level)
    {
        float *destination =
            local + (size_t)level * (size_t)PLASIC_NSPP;

        if (valid_count > 0)
        {
            memcpy(
                destination,
                global + (size_t)level * (size_t)PLASIC_NRSP + rank_offset,
                (size_t)valid_count * sizeof(float));
        }
        if (padding_count > 0)
        {
            memset(
                destination + valid_count, 0,
                padding_count * sizeof(float));
        }
    }
}

static int initialize_cold_start_atmosphere(
    const char *surface_data_directory)
{
    int level;
    int status;

    status = surface_data_read_local(
        surface_data_directory, PLASIC_NTRU, PLASIC_NLAT, PLASIC_NLON,
        "orography", 1, mp_local_offset(PLASIC_NHOR), PLASIC_NHOR,
        runtime_fields->cold_start_orography_grid);
    if (status != SURFACE_DATA_OK)
    {
        runtime_set_error(
            "cold start cannot read NetCDF surface field 'orography' "
            "from '%s' (status %d)",
            surface_data_directory, status);
        return PLASIC_RUNTIME_IO_FAILED;
    }

    status = sht_grid_to_scalar(
        1, runtime_fields->cold_start_orography_grid, plasic_so);
    if (status != PLASIC_SHT_OK)
    {
        runtime_set_error("cold-start orography transform failed (status %d)", status);
        return PLASIC_RUNTIME_NUMERICAL_FAILURE;
    }

    /*
     * 对全球谱地形做无量纲化预处理
     *
     * Preprocess the global spectral orography into non-dimensional form
     */
    surface_prepare_spectral(PLASIC_NRSP, plasic_cv, plasic_so);
    split_global_spectral(plasic_so, plasic_sop, 1);

    /*
     * 谱对数地表气压与地形静力平衡
     * Spectral log surface pressure and hydrostatic balance with the orography
     */
    {
        const float balance =
            -plasic_cv * plasic_cv /
            (plasic_gascon * plasic_tgr);
        int32_t index;

        for (index = 2; index < PLASIC_NRSP; ++index)
        {
            plasic_sp[index] = balance * plasic_so[index];
        }

        for (level = 0; level < PLASIC_NLEV; ++level)
        {
            plasic_sz[(size_t)level * (size_t)PLASIC_NRSP + 2U] =
                plasic_plavor;
        }
    }

    /*
     * 冷启动没有历史层，故把当前谱状态复制为跃蛙的旧时间层。
     * Cold start has no history levels, so the current spectral state is copied to the leapfrog old time level.
     */
    {
        split_global_spectral(plasic_sz, plasic_szp, PLASIC_NLEV);
        split_global_spectral(plasic_sp, plasic_spp, 1);
        memcpy(plasic_szm, plasic_szp, sizeof(plasic_szm));
        memcpy(plasic_spm, plasic_spp, sizeof(plasic_spm));
    }

    for (level = 0; level < PLASIC_NLEV; ++level)
    {
        const size_t grid_offset = (size_t)level * (size_t)PLASIC_NHOR;
        int32_t point;

        for (point = 0; point < PLASIC_NHOR; ++point)
        {
            plasic_dt[grid_offset + (size_t)point] = plasic_t0[level];
        }
    }

    {
        const size_t surface_offset = (size_t)PLASIC_NHOR * PLASIC_NLEV;
        int32_t point;

        for (point = 0; point < PLASIC_NHOR; ++point)
        {
            plasic_dt[surface_offset + (size_t)point] =
                plasic_tgr;
        }
    }

    /*
     * 冷启动的 nstep 由调用方根据 --calendar 与 --start-* 参数先用
     * cal2step() 算好，本函数只负责重置输出累积计数。
     *
     * The cold-start nstep is computed by the caller from --calendar and the
     * --start-* options before this function runs; only the output
     * accumulation counter is reset here.
     */
    plasic_naccuout = 0;
    return PLASIC_RUNTIME_OK;
}

int runtime_load_atmosphere_restart(
    const char *path, const char *surface_data_directory)
{
    /*
     * 从 restart 文件恢复大气状态（热启动路径）。
     * Restore the atmospheric state from a restart file (warm-start path).
     */
    int32_t stored_nlat;
    int32_t stored_nlon;
    int32_t stored_nlev;
    int32_t stored_nrsp;
    int status;

    if (path == NULL)
    {
        return initialize_cold_start_atmosphere(surface_data_directory);
    }

    status = restart_open_read(path);
    if (status != PLASIC_RESTART_OK)
    {
        runtime_set_error("cannot open restart '%s' (status %d)", path, status);
        return PLASIC_RUNTIME_RESTART_FAILED;
    }

/*
 * 下面三个宏封装“读取 + 失败即跳转 failed 清理退出”的样板代码：
 *   REQUIRE_INTEGER(name, dest)          读一个标量整数；
 *   REQUIRE_ARRAY(name, dest, rows, columns)
 *                                        读整个数组（每个 rank 都读完整份），
 *                                        rows×columns 个元素按行维展平存放；
 *   REQUIRE_SLICE(name, dest, global_rows, local_offset, local_rows, columns)
 *                                        只读全局数组中属于本 rank 的
 *                                        [local_offset, local_offset+local_rows) 行，
 *                                        用于 *m 谱块和格点条带的并行切片读取；
 *
 * The following three macros encapsulate the boilerplate of "read and jump to
 * failed cleanup on error":
 *   REQUIRE_INTEGER(name, dest)          read one scalar integer;
 *   REQUIRE_ARRAY(name, dest, rows, columns)
 *                                        read an entire array (each rank reads
 *                                        the full copy), with rows×columns
 *                                        elements stored flattened along the
 *                                        row dimension;
 *   REQUIRE_SLICE(name, dest, global_rows, local_offset, local_rows, columns)
 *                                        read only the rows
 *                                        [local_offset, local_offset+local_rows)
 *                                        of the global array that belong to
 *                                        this rank, used for the *m spectral
 *                                        blocks and grid strips;
 */
#define REQUIRE_INTEGER(name, destination)                       \
    do                                                           \
    {                                                            \
        status = require_restart_integer((name), (destination)); \
        if (status != PLASIC_RUNTIME_OK)                         \
        {                                                        \
            goto failed;                                         \
        }                                                        \
    } while (0)

#define REQUIRE_ARRAY(name, destination, rows, columns)               \
    do                                                                \
    {                                                                 \
        status = require_restart_array((name), (destination), (rows), \
                                       (columns));                    \
        if (status != PLASIC_RUNTIME_OK)                              \
        {                                                             \
            goto failed;                                              \
        }                                                             \
    } while (0)
#define REQUIRE_SLICE(name, destination, global_rows, local_offset, \
                      local_rows, columns)                          \
    do                                                              \
    {                                                               \
        status = require_restart_array_slice(                       \
            (name), (destination), (global_rows), (local_offset),   \
            (local_rows), (columns));                               \
        if (status != PLASIC_RUNTIME_OK)                            \
        {                                                           \
            goto failed;                                            \
        }                                                           \
    } while (0)

    REQUIRE_INTEGER("nlat", &stored_nlat);
    REQUIRE_INTEGER("nlon", &stored_nlon);
    REQUIRE_INTEGER("nlev", &stored_nlev);
    REQUIRE_INTEGER("nrsp", &stored_nrsp);
    if (stored_nlat != PLASIC_NLAT || stored_nlon != PLASIC_NLON ||
        stored_nlev != PLASIC_NLEV || stored_nrsp != PLASIC_NRSP)
    {
        runtime_set_error("restart grid is N%d/%d levels/%d spectral values, "
                          "but this executable is N%d/%d/%d",
                          stored_nlat, stored_nlev, stored_nrsp,
                          PLASIC_NLAT, PLASIC_NLEV, PLASIC_NRSP);
        status = PLASIC_RUNTIME_RESTART_FAILED;
        goto failed;
    }

    /*
     * nstep 是模型唯一时钟，强制要求存在。
     * calendar 决定 nstep↔日期的映射：新文件总是写出，旧文件缺失时按
     * Gregorian（0）处理，因为旧版本本来就是 Gregorian 行为。
     */
    status = require_restart_integer_mandatory("nstep", &plasic_nstep);
    if (status != PLASIC_RUNTIME_OK)
    {
        goto failed;
    }
    if (plasic_nstep < 0)
    {
        runtime_set_error("restart nstep %d is negative", plasic_nstep);
        status = PLASIC_RUNTIME_RESTART_FAILED;
        goto failed;
    }
    {
        int32_t stored_calendar = PLASIC_CALENDAR_GREGORIAN;
        const int32_t calendar_status =
            restart_read_integer("calendar", &stored_calendar);

        if (calendar_status == PLASIC_RESTART_NOT_FOUND)
        {
            stored_calendar = PLASIC_CALENDAR_GREGORIAN;
        }
        else if (calendar_status != PLASIC_RESTART_OK)
        {
            runtime_set_error(
                "restart calendar record could not be read (status %d)",
                calendar_status);
            status = PLASIC_RUNTIME_RESTART_FAILED;
            goto failed;
        }
        else if (!calendar_is_valid(stored_calendar))
        {
            runtime_set_error(
                "restart declares unknown calendar kind %d",
                stored_calendar);
            status = PLASIC_RUNTIME_RESTART_FAILED;
            goto failed;
        }
        plasic_calendar = stored_calendar;
    }
    REQUIRE_INTEGER("naccuout", &plasic_naccuout);
    status = restart_read_seed(
        "seed", plasic_seed,
        (int32_t)(sizeof(plasic_seed) / sizeof(plasic_seed[0])));

    REQUIRE_ARRAY("sz", plasic_sz, PLASIC_NRSP, PLASIC_NLEV);
    REQUIRE_ARRAY("sd", plasic_sd, PLASIC_NRSP, PLASIC_NLEV);
    REQUIRE_ARRAY("st", plasic_st, PLASIC_NRSP, PLASIC_NLEV);
    REQUIRE_ARRAY("sq", plasic_sq, PLASIC_NRSP, PLASIC_NLEV);
    REQUIRE_ARRAY("sr", plasic_sr, PLASIC_NRSP, PLASIC_NLEV);
    REQUIRE_ARRAY("sp", plasic_sp, PLASIC_NRSP, 1);
    REQUIRE_ARRAY("so", plasic_so, PLASIC_NRSP, 1);

    REQUIRE_SLICE(
        "szm", plasic_szm, PLASIC_NRSP,
        mp_rank() * PLASIC_NSPP, PLASIC_NSPP,
        PLASIC_NLEV);
    REQUIRE_SLICE(
        "sdm", plasic_sdm, PLASIC_NRSP,
        mp_rank() * PLASIC_NSPP, PLASIC_NSPP,
        PLASIC_NLEV);
    REQUIRE_SLICE(
        "stm", plasic_stm, PLASIC_NRSP,
        mp_rank() * PLASIC_NSPP, PLASIC_NSPP,
        PLASIC_NLEV);
    REQUIRE_SLICE(
        "sqm", plasic_sqm, PLASIC_NRSP,
        mp_rank() * PLASIC_NSPP, PLASIC_NSPP,
        PLASIC_NLEV);
    REQUIRE_SLICE(
        "spm", plasic_spm, PLASIC_NRSP,
        mp_rank() * PLASIC_NSPP, PLASIC_NSPP, 1);

#define REQUIRE_GRID(name, destination, columns)   \
    REQUIRE_SLICE(                                 \
        (name), (destination), PLASIC_NUGP,        \
        mp_local_offset(PLASIC_NHOR), PLASIC_NHOR, \
        (columns))

    REQUIRE_GRID("dls", plasic_dls, 1);
    REQUIRE_GRID("dwetfac", plasic_dwetfac, 1);
    REQUIRE_GRID("dalb", plasic_dalb, 1);
    REQUIRE_GRID("dz0", plasic_dz0, 1);
    REQUIRE_GRID("dicec", plasic_dicec, 1);
    REQUIRE_GRID("diced", plasic_diced, 1);
    REQUIRE_GRID("dwatc", plasic_dwatc, 1);
    REQUIRE_GRID("dust3", plasic_dust3, 1);
    REQUIRE_GRID("dcc", plasic_dcc, PLASIC_NLEV);
    REQUIRE_GRID("dql", plasic_dql, PLASIC_NLEV);
    REQUIRE_GRID("dqsat", plasic_dqsat, PLASIC_NLEV);
    REQUIRE_GRID(
        "dt", plasic_dt + PLASIC_NHOR * PLASIC_NLEV, 1);
    REQUIRE_GRID(
        "dq", plasic_dq + PLASIC_NHOR * PLASIC_NLEV, 1);

#define LOAD_ACCUMULATOR(name, destination) \
    REQUIRE_GRID((name), (destination), 1)

    LOAD_ACCUMULATOR("aprl", plasic_aprl);
    LOAD_ACCUMULATOR("aprc", plasic_aprc);
    LOAD_ACCUMULATOR("aprs", plasic_aprs);
    LOAD_ACCUMULATOR("aevap", plasic_aevap);
    LOAD_ACCUMULATOR("ashfl", plasic_ashfl);
    LOAD_ACCUMULATOR("alhfl", plasic_alhfl);
    LOAD_ACCUMULATOR("asmelt", plasic_asmelt);
    LOAD_ACCUMULATOR("asndch", plasic_asndch);
    LOAD_ACCUMULATOR("acc", plasic_acc);
    LOAD_ACCUMULATOR("assol", plasic_assol);
    LOAD_ACCUMULATOR("asthr", plasic_asthr);
    LOAD_ACCUMULATOR("atsol", plasic_atsol);
    LOAD_ACCUMULATOR("atthr", plasic_atthr);
    LOAD_ACCUMULATOR("ataux", plasic_ataux);
    LOAD_ACCUMULATOR("atauy", plasic_atauy);
    LOAD_ACCUMULATOR("atsolu", plasic_atsolu);
    LOAD_ACCUMULATOR("assolu", plasic_assolu);
    LOAD_ACCUMULATOR("asthru", plasic_asthru);
    LOAD_ACCUMULATOR("aqvi", plasic_aqvi);
    LOAD_ACCUMULATOR("atsa", plasic_atsa);
    LOAD_ACCUMULATOR("ats0", plasic_ats0);
    LOAD_ACCUMULATOR("atsama", plasic_atsama);
    LOAD_ACCUMULATOR("atsami", plasic_atsami);
#undef LOAD_ACCUMULATOR
#undef REQUIRE_GRID

    status = restart_close();

    split_global_spectral(plasic_sz, plasic_szp, PLASIC_NLEV);
    split_global_spectral(plasic_sd, plasic_sdp, PLASIC_NLEV);
    split_global_spectral(plasic_st, plasic_stp, PLASIC_NLEV);
    split_global_spectral(plasic_sq, plasic_sqp, PLASIC_NLEV);
    split_global_spectral(plasic_sr, plasic_srp, PLASIC_NLEV);
    split_global_spectral(plasic_sp, plasic_spp, 1);
    split_global_spectral(plasic_so, plasic_sop, 1);
    return PLASIC_RUNTIME_OK;

failed:
    (void)restart_close();
    return status;

#undef REQUIRE_ARRAY
#undef REQUIRE_INTEGER
#undef REQUIRE_SLICE
}

int runtime_prepare_orography_grid(void)
{
    size_t horizontal;
    int status;

    status = sht_scalar_to_grid(
        1, plasic_so, runtime_fields->surface_geopotential);

    if (status != PLASIC_SHT_OK)
    {
        runtime_set_error("orography synthesis failed (status %d)", status);
        return PLASIC_RUNTIME_NUMERICAL_FAILURE;
    }

    for (horizontal = 0; horizontal < PLASIC_NHOR; ++horizontal)
    {
        runtime_fields->surface_geopotential[horizontal] =
            runtime_fields->surface_geopotential[horizontal] *
            plasic_cv * plasic_cv;
    }
    return PLASIC_RUNTIME_OK;
}

/*
 * 时间步长。
 *
 * 用户给出步长（--step-minutes，或 JSON 配置中的 step_minutes）时直接使用，
 * 不检验其稳定性：步长是否合理完全由用户负责。未给出时使用当前分辨率的
 * 默认步长 PLASIC_DEFAULT_STEP_SECONDS（见 state/plasicmod_state.h）：
 * T21 30、T31 18、T42 12、T85 7.5、T170 2.5 分钟，保证长时间气候积分稳定、
 * 不发生数值发散。
 *
 * 日历按“每天整数步”计时（calmod_init(ntspd)），所以步长必须把一个太阳日
 * 分成整数步。做不到的步长无法表示，这里直接报错，而不是悄悄换成另一个步长。
 *
 * Time step.
 *
 * A step given by the user (--step-minutes, or step_minutes in the JSON
 * config) is used as given, without any stability check: whether it is
 * reasonable is entirely the user's responsibility. Otherwise the default of
 * this resolution, PLASIC_DEFAULT_STEP_SECONDS in state/plasicmod_state.h, is
 * used: 30 minutes at T21, 18 at T31, 12 at T42, 7.5 at T85 and 2.5 at T170,
 * short enough to keep long climate integrations stable.
 *
 * The calendar counts a whole number of steps per day (calmod_init(ntspd)),
 * so the step must divide the solar day into whole steps. A step that cannot
 * is rejected rather than silently replaced by a different one.
 */
static int select_time_step(double step_minutes)
{
    // lroundf 将单精度浮点数四舍五入为长整型整数
    // lroundf rounds a single-precision float to a long integer
    const double solar_day = (double)lroundf(plasic_solar_day);
    double step_seconds = (double)PLASIC_DEFAULT_STEP_SECONDS;
    double steps_per_day;
    double whole_steps;

    if (!isfinite(step_minutes) || step_minutes < 0.0)
    {
        runtime_set_error("invalid time step of %g minutes", step_minutes);
        return PLASIC_RUNTIME_BAD_ARGUMENT;
    }
    if (step_minutes > 0.0)
    {
        step_seconds = step_minutes * 60.0;
    }

    steps_per_day = solar_day / step_seconds;
    whole_steps = round(steps_per_day);
    if (whole_steps < 1.0 || whole_steps > (double)INT32_MAX ||
        fabs(steps_per_day - whole_steps) > 1.0e-9 * steps_per_day)
    {
        runtime_set_error(
            "a time step of %g minutes does not divide the %g-minute day "
            "into a whole number of steps (%.6g steps per day)",
            step_seconds / 60.0, solar_day / 60.0, steps_per_day);
        return PLASIC_RUNTIME_BAD_ARGUMENT;
    }

    plasic_ntspd = (int32_t)whole_steps;
    plasic_mpstep = (float)(step_seconds / 60.0);
    return PLASIC_RUNTIME_OK;
}

int runtime_initialize_dynamics(double step_minutes)
{
    int status;
    int latitude;
    int longitude;

    planet_defaults(
        &plasic_sidereal_day, &plasic_solar_day, &plasic_kap,
        &plasic_gascon, &plasic_ra1, &plasic_ra2,
        &plasic_ra4, &plasic_pnu, &plasic_ga,
        &plasic_plarad);

    status = select_time_step(step_minutes);
    if (status != PLASIC_RUNTIME_OK)
    {
        return status;
    }
    plasic_nafter = plasic_ntspd;
    plasic_nkits = 0;

    /*
     * 日历算术只有一个实现（tools/calmod_kernels.c）。这里只把“每天多少步”
     * 写入日历模块，日历类型由冷启动参数或 restart 文件决定。
     *
     * Calendar arithmetic lives only in tools/calmod_kernels.c. Only the
     * number of steps per day is registered here; the calendar kind itself is
     * selected by the cold-start option or by the restart file.
     */
    if (calmod_init(plasic_ntspd) != CALMOD_OK)
    {
        runtime_set_error(
            "could not initialize the calendar with %d steps per day",
            plasic_ntspd);
        return PLASIC_RUNTIME_INITIALIZATION_FAILED;
    }

    plasic_ww = plasic_two_pi / plasic_sidereal_day;
    plasic_cpd = plasic_gascon / plasic_kap;
    plasic_cpv_cpd_minus1 = plasic_water_heat_capacity / plasic_cpd - 1.0f;
    plasic_cv = plasic_plarad * plasic_ww;
    plasic_ct = plasic_cv * plasic_cv / plasic_gascon;
    plasic_pnu21 = 1.0f - 2.0f * plasic_pnu;
    plasic_rdbrv = plasic_gascon / plasic_water_gas_constant;

    plasic_deltsec = plasic_solar_day / (float)plasic_ntspd;
    plasic_deltsec2 = 2.0f * plasic_deltsec;
    plasic_delt = plasic_two_pi / (float)plasic_ntspd;
    plasic_delt2 = 2.0f * plasic_delt;

    gauss_inigau(PLASIC_NLAT, plasic_sid, plasic_gwd);
    inilat(PLASIC_NLAT, plasic_sid, plasic_csq, plasic_rcs);

    /*
     * MPI 下每个 rank 只保存连续的一段纬度。初始化例程先生成全局纬向数组，
     * 再把当前 rank 对应的 NLPP 行移到数组开头；rank 0 已经从 offset=0 开始，
     * 因而不需要移动。
     *
     * Under MPI, each rank stores only one contiguous segment of latitudes.
     * The initialization routines first generate the global latitude arrays,
     * then move the NLPP rows belonging to the current rank to the beginning
     * of the arrays; rank 0 already starts at offset=0 and therefore needs no
     * move.
     */
    if (!mp_is_root())
    {
        const size_t latitude_offset =
            (size_t)mp_rank() * (size_t)PLASIC_NLPP;

        memmove(
            plasic_sid, plasic_sid + latitude_offset,
            (size_t)PLASIC_NLPP * sizeof(plasic_sid[0]));
        memmove(
            plasic_gwd, plasic_gwd + latitude_offset,
            (size_t)PLASIC_NLPP * sizeof(plasic_gwd[0]));
        memmove(
            plasic_csq, plasic_csq + latitude_offset,
            (size_t)PLASIC_NLPP * sizeof(plasic_csq[0]));
        memmove(
            plasic_rcs, plasic_rcs + latitude_offset,
            (size_t)PLASIC_NLPP * sizeof(plasic_rcs[0]));
    }

    for (latitude = 0; latitude < PLASIC_NLPP; ++latitude)
    {
        plasic_cola[latitude] = sqrtf(plasic_csq[latitude]);
        for (longitude = 0; longitude < PLASIC_NLON; ++longitude)
        {
            plasic_rcsq[latitude * PLASIC_NLON + longitude] =
                1.0f / plasic_csq[latitude];
        }
    }

    initpm_vertical(
        PLASIC_NLEV, plasic_neqsig, plasic_sigmah_configured,
        plasic_sigmah, plasic_dsigma, plasic_rdsig,
        plasic_sigma);

    /*
     * 高层海绵层（模式顶 Rayleigh 摩擦）。
     *
     * 参考 10 层配置只对模式顶两层施加线性摩擦：sigma 界面约 0.077 的
     * 顶层取 20 天，界面约 0.162 的第 2 层取 100 天，覆盖约 16% 的大气质量。
     * 原 C 实现用 `if (PLASIC_NLEV == 10)` 限定这组预设，其他层数（例如
     * T85/L25）的模式顶完全没有阻尼，高层半隐式离散的计算模态会从舍入误差
     * 放大为网格尺度数值爆炸（T85/L25 分别在约第 160 天和第 200 天出现了
     * 模式顶/平流层下部的两次数值爆炸；动力学单独积分也能复现）。
     *
     * 这里按“覆盖相同大气质量”的原则把该预设推广到所有层数：
     *   sigma 界面 <= 0.095 的层：20 天；
     *   sigma 界面 <= 0.17  的层：100 天；
     *   其余层保持 0。
     * 对 10 层配置，该规则恰好复现原来的 20/100 天两层预设；层数增加时
     * 海绵层自动向下延伸，保持与参考配置一致的高层阻尼厚度。数值单位为天，
     * initpm_spectral() 会把正的 τ 换算为无量纲阻尼率 1/(2πτ)，再作用于
     * 散度和涡度倾向。
     *
     * Upper sponge layer (Rayleigh friction at the model top). The reference
     * 10-level configuration damps only its top two levels: 20 days for the
     * top interface (sigma ~ 0.077) and 100 days for the second (sigma ~
     * 0.162), covering about 16% of the atmospheric mass. The previous C
     * implementation enabled this preset only for NLEV==10, leaving every
     * other level count (for example T85/L25) with a completely undamped model
     * top; the top-level computational mode of the semi-implicit
     * discretisation then amplified from round-off into a grid-scale blow-up
     * (T85/L25 produced such blow-ups near day 160 at the top level and near
     * day 200 in the lower stratosphere; both also appear in dynamics-only
     * integrations).
     *
     * The preset is generalised to every level count with the same damped
     * mass fraction:
     *   levels with sigma interface <= 0.095 : 20 days;
     *   levels with sigma interface <= 0.17  : 100 days;
     *   all deeper levels                     : 0 (no friction).
     * For the 10-level configuration this reproduces the original 20/100-day
     * two-level preset exactly; deeper grids extend the sponge downwards so
     * that the damped layer keeps the reference mass depth. Values are in
     * days; initpm_spectral() converts a positive tau into the
     * non-dimensional damping rate 1/(2*pi*tau) applied to the divergence and
     * vorticity tendencies.
     */
    {
        int32_t level;

        for (level = 0; level < PLASIC_NLEV; ++level)
        {
            const float sigma_interface = plasic_sigmah[level];

            if (sigma_interface <= 0.095f)
            {
                plasic_tfrc[level] = 20.0f;
            }
            else if (sigma_interface <= 0.17f)
            {
                plasic_tfrc[level] = 100.0f;
            }
        }
    }

    initpm_spectral(
        PLASIC_NLEV, PLASIC_NTRU, plasic_nhdiff, PLASIC_NRSP,
        plasic_two_pi, plasic_ct,
        plasic_ndel, plasic_restim, plasic_tfrc,
        plasic_tdissd, plasic_tdissz, plasic_tdisst, plasic_tdissq,
        plasic_t0, plasic_damp,
        runtime_fields->friction_rate,
        runtime_fields->divergence_diffusion_rate,
        runtime_fields->vorticity_diffusion_rate,
        runtime_fields->temperature_diffusion_rate,
        runtime_fields->humidity_diffusion_rate, plasic_Lnk,
        plasic_nindex, plasic_spnorm,
        runtime_fields->nondimensional_reference_temperature);

    initsi(
        PLASIC_NLEV, plasic_kap,
        runtime_fields->nondimensional_reference_temperature,
        plasic_sigmah,
        plasic_dsigma,
        plasic_rdsig, plasic_tkp, plasic_t01s2,
        plasic_g, plasic_c, plasic_tau);

    makebm(
        PLASIC_NLEV, PLASIC_NTRU, plasic_delt,
        runtime_fields->nondimensional_reference_temperature,
        plasic_dsigma,
        plasic_g, plasic_tau, plasic_bm1);

    status = sht_init(
        PLASIC_NTRU, PLASIC_NLAT, PLASIC_NLON, PLASIC_NLPP,
        mp_rank() * PLASIC_NLPP);
    if (status != PLASIC_SHT_OK)
    {
        runtime_set_error("SHTns initialization failed (status %d)", status);
        return PLASIC_RUNTIME_INITIALIZATION_FAILED;
    }

    split_global_spectral(plasic_Lnk, plasic_Lnkpp, PLASIC_NLEV);
    return PLASIC_RUNTIME_OK;
}
