/*
 * 月度 NetCDF 输出的运行时实现：把每个时间步的瞬时诊断量累加为月平均，
 * 并按 CF/CMIP 约定写成一组 NetCDF 文件。
 * Runtime implementation of monthly NetCDF output: it accumulates per-step
 * instantaneous diagnostics into monthly means and writes them as a set of
 * CF/CMIP-style NetCDF files.
 *
 * 工作流程 / Workflow:
 *   open()       校验输出目录、构建坐标与缓冲区、聚合陆地掩膜，并启用输出。
 *                Validate the output directory, build coordinates and buffers,
 *                gather the land mask and enable output.
 *   accumulate() 每个时间步把瞬时诊断量累加到当前月份的和；跨月时先写上月。
 *                Add each time step's instantaneous diagnostics to the current
 *                month's sums, flushing the previous month when the date rolls
 *                over.
 *   close()      写出最后一个完整月份，关闭全部文件并释放资源；可重复调用。
 *                Flush the final complete month, close every file and release
 *                resources; safe to call repeatedly.
 *   abort()      不写盘直接丢弃状态并释放资源，用于打开失败或异常退出。
 *                Discard state and release resources without writing; used on
 *                open failures and abnormal termination.
 *
 * 并行语义 / Parallel semantics: 每个进程累加自己的本地列；flush 时用
 * mp_allgather_grid() 拼成全局场，只有根进程持有文件句柄并执行 NetCDF I/O。
 * Every rank accumulates its local columns; flush gathers a global field with
 * mp_allgather_grid(), and only the root holds file handles and performs
 * NetCDF I/O.
 */

#include "runtime_internal.h"

#include "co2_forcing.h"
#include "ozone_forcing.h"
#include "monthly_netcdf.h"
#include "tools/monthly_netcdf_writer.h"

#include <sys/stat.h>

enum
{
    /* 输出文件路径的最大字节数。 */
    /* Maximum number of bytes in an output file path. */
    MONTHLY_PATH_BYTES = 4096,
    /* 掩膜策略：不掩膜 / 只保留海洋 / 只保留陆地。 */
    /* Masking policy: none / ocean only / land only. */
    MASK_NONE = 0,
    MASK_OCEAN = 1,
    MASK_LAND = 2
};

/*
 * 每个输出变量对应的模式内部诊断量来源：source_value() 根据这些枚举值从模式
 * 状态数组或物理模块中读取原始场。
 * Model-internal diagnostic source of each output variable: source_value()
 * uses these enumerators to read the raw field from the model state arrays or
 * the physics module.
 */
typedef enum monthly_source
{
    /* 大气三维场：纬向风、经向风、温度、比湿、云量。 */
    /* Atmospheric 3-D fields: zonal wind, meridional wind, temperature,
     * specific humidity and cloud fraction. */
    SOURCE_UA,
    SOURCE_VA,
    SOURCE_TA,
    SOURCE_HUS,
    SOURCE_CL,
    /* 地表大气场：地面气压、地表温度、近地面气温。 */
    /* Surface atmospheric fields: surface pressure, surface skin temperature
     * and near-surface air temperature. */
    SOURCE_PS,
    SOURCE_TS,
    SOURCE_TAS,
    /* 海洋与海冰表面场：海表温度、海冰密集度、海冰厚度。 */
    /* Ocean and sea-ice surface fields: sea-surface temperature, sea-ice
     * concentration and sea-ice thickness. */
    SOURCE_TOS,
    SOURCE_SICONC,
    SOURCE_SITHICK,
    /* 陆面场：积雪深度与土壤水储量。 */
    /* Land-surface fields: snow depth and soil water storage. */
    SOURCE_SND,
    SOURCE_MRSO,
    /* 地表通量：降水、蒸发、感热与潜热通量。 */
    /* Surface fluxes: precipitation, evaporation, sensible and latent heat. */
    SOURCE_PR,
    SOURCE_EVSPSBL,
    SOURCE_HFSS,
    SOURCE_HFLS,
    /* 辐射通量：大气顶与地表的短波、长波辐射。 */
    /* Radiation fluxes: shortwave and longwave at the top of atmosphere and
     * at the surface. */
    SOURCE_RSDT,
    SOURCE_RSUT,
    SOURCE_RLUT,
    SOURCE_RSDS,
    SOURCE_RSUS,
    SOURCE_RLDS,
    SOURCE_RLUS,
    /* 地表反照率、海洋混合层深度、深海温度与盐度、CO2 浓度。 */
    /* Surface albedo, ocean mixed-layer depth, deep-ocean temperature and
     * salinity, and CO2 concentration. */
    SOURCE_ALBEDO,
    SOURCE_MLD,
    SOURCE_CO2,
    SOURCE_THETAO,
    SOURCE_SO
} monthly_source;

/*
 * 单个输出变量的静态定义。
 * Static definition of one output variable
 */
typedef struct monthly_definition
{
    /* NetCDF 变量名，同时用作输出文件名 <name>.nc。 */
    /* NetCDF variable name, also used as the output file name <name>.nc. */
    const char *name;
    /* CF standard_name 属性。 */
    /* CF standard_name attribute. */
    const char *standard_name;
    /* 人类可读的 long_name 属性。 */
    /* Human-readable long_name attribute. */
    const char *long_name;
    /* 变量的物理单位。 */
    /* Physical units of the variable. */
    const char *units;
    /* 读取原始诊断值所用的来源枚举。 */
    /* Source enumerator used to read the raw diagnostic value. */
    monthly_source source;
    /* 垂直层数：1 为单层场，>1 为带垂直坐标的三维场。 */
    /* Number of vertical layers: 1 for a single-level field, >1 for a 3-D
     * field that carries a vertical coordinate. */
    int32_t layers;
    /* 掩膜策略：MASK_NONE / MASK_OCEAN / MASK_LAND。 */
    /* Masking policy: MASK_NONE / MASK_OCEAN / MASK_LAND. */
    int32_t mask;
} monthly_definition;

/*
 * 单个输出变量的运行时实例
 * Runtime instance of one output variable
 */
typedef struct monthly_variable
{
    /* 指向 definitions[] 中对应的静态定义。 */
    /* Points to the matching static definition in definitions[]. */
    const monthly_definition *definition;
    /* 根进程持有的 NetCDF 句柄；其他进程的 ncid 始终为 -1。 */
    /* NetCDF handle held by the root rank; ncid is always -1 elsewhere. */
    monthly_netcdf_file file;
    /* 本地列上的时间累加和，尺寸为 PLASIC_NHOR * layers。 */
    /* Time-accumulated sums over the local columns, sized PLASIC_NHOR * layers. */
    float *sum;
} monthly_variable;

/* 大气变量统一使用模式的垂直层数。 */
/* Atmospheric variables use the model's vertical layer count. */
#define ATMOS_LEVELS PLASIC_NLEV

/*
 * 全部月度输出变量的定义表，顺序即变量枚举顺序。
 * 每列依次为：name、standard_name、long_name、units、source、layers、mask。
 * Table of all monthly output variables, in the same order as the source
 * enum. The columns are: name, standard_name, long_name, units, source,
 * layers, mask.
 */
static const monthly_definition definitions[] = {
    {"ua", "eastward_wind", "Zonal wind", "m s-1", SOURCE_UA, ATMOS_LEVELS, MASK_NONE},
    {"va", "northward_wind", "Meridional wind", "m s-1", SOURCE_VA, ATMOS_LEVELS, MASK_NONE},
    {"ta", "air_temperature", "Air temperature", "K", SOURCE_TA, ATMOS_LEVELS, MASK_NONE},
    {"hus", "specific_humidity", "Specific humidity", "kg kg-1", SOURCE_HUS, ATMOS_LEVELS, MASK_NONE},
    {"cl", "cloud_area_fraction_in_atmosphere_layer", "Cloud fraction", "1", SOURCE_CL, ATMOS_LEVELS, MASK_NONE},
    {"ps", "surface_air_pressure", "Surface pressure", "Pa", SOURCE_PS, 1, MASK_NONE},
    {"ts", "surface_temperature", "Surface skin temperature", "K", SOURCE_TS, 1, MASK_NONE},
    {"tas", "air_temperature", "Near-surface air temperature", "K", SOURCE_TAS, 1, MASK_NONE},
    {"tos", "sea_surface_temperature", "Sea surface temperature", "K", SOURCE_TOS, 1, MASK_OCEAN},
    {"siconc", "sea_ice_area_fraction", "Sea-ice area fraction", "1", SOURCE_SICONC, 1, MASK_OCEAN},
    {"sithick", "sea_ice_thickness", "Sea-ice thickness", "m", SOURCE_SITHICK, 1, MASK_OCEAN},
    {"snd", "surface_snow_thickness", "Snow depth", "m", SOURCE_SND, 1, MASK_NONE},
    {"mrso", "lwe_thickness_of_soil_moisture_content", "Soil water storage", "m", SOURCE_MRSO, 1, MASK_LAND},
    {"pr", "precipitation_flux", "Total precipitation flux", "kg m-2 s-1", SOURCE_PR, 1, MASK_NONE},
    {"evspsbl", "water_evaporation_flux", "Surface evaporation flux", "kg m-2 s-1", SOURCE_EVSPSBL, 1, MASK_NONE},
    {"hfss", "surface_upward_sensible_heat_flux", "Surface sensible heat flux", "W m-2", SOURCE_HFSS, 1, MASK_NONE},
    {"hfls", "surface_upward_latent_heat_flux", "Surface latent heat flux", "W m-2", SOURCE_HFLS, 1, MASK_NONE},
    {"rsdt", "toa_incoming_shortwave_flux", "TOA incoming shortwave radiation", "W m-2", SOURCE_RSDT, 1, MASK_NONE},
    {"rsut", "toa_outgoing_shortwave_flux", "TOA outgoing shortwave radiation", "W m-2", SOURCE_RSUT, 1, MASK_NONE},
    {"rlut", "toa_outgoing_longwave_flux", "TOA outgoing longwave radiation", "W m-2", SOURCE_RLUT, 1, MASK_NONE},
    {"rsds", "surface_downwelling_shortwave_flux_in_air", "Surface downwelling shortwave radiation", "W m-2", SOURCE_RSDS, 1, MASK_NONE},
    {"rsus", "surface_upwelling_shortwave_flux_in_air", "Surface upwelling shortwave radiation", "W m-2", SOURCE_RSUS, 1, MASK_NONE},
    {"rlds", "surface_downwelling_longwave_flux_in_air", "Surface downwelling longwave radiation", "W m-2", SOURCE_RLDS, 1, MASK_NONE},
    {"rlus", "surface_upwelling_longwave_flux_in_air", "Surface upwelling longwave radiation", "W m-2", SOURCE_RLUS, 1, MASK_NONE},
    {"albedo", "surface_albedo", "Surface albedo", "1", SOURCE_ALBEDO, 1, MASK_NONE},
    {"mlotst", "ocean_mixed_layer_thickness", "Ocean mixed-layer depth", "m", SOURCE_MLD, 1, MASK_OCEAN},
    {"thetao", "sea_water_potential_temperature", "Sea water potential temperature", "K", SOURCE_THETAO, PLASIC_DEEP_NLEV, MASK_OCEAN},
    {"so", "sea_water_salinity", "Sea water salinity", "1e-3", SOURCE_SO, PLASIC_DEEP_NLEV, MASK_OCEAN},
    {"co2", "mole_fraction_of_carbon_dioxide_in_air", "Prescribed atmospheric carbon dioxide mole fraction", "1e-6", SOURCE_CO2, 1, MASK_NONE}};

/*
 * 三维海洋变量使用深海深度坐标轴，而不是大气 sigma 轴。
 * 该标志由 source 推导，使上面的变量表保持紧凑。
 * Three-dimensional ocean variables use the deep-ocean depth axis rather than
 * the atmospheric sigma axis.  The flag is derived from the source so the
 * variable table above stays compact.
 */
static int is_depth_source(monthly_source source)
{
    return source == SOURCE_THETAO || source == SOURCE_SO;
}

enum
{
    /* 输出变量总数，由定义表长度自动推导。 */
    /* Total number of output variables, derived from the definition table. */
    MONTHLY_VARIABLE_COUNT =
        (int)(sizeof(definitions) / sizeof(definitions[0]))
};

static monthly_variable variables[MONTHLY_VARIABLE_COUNT];
static float *global_buffer;
static float *global_land_mask;
static float monthly_latitudes[PLASIC_NLAT];
static float monthly_longitudes[PLASIC_NLON];
static char monthly_directory[MONTHLY_PATH_BYTES];
static int monthly_enabled;
static int files_created;
static int active_year;
static int active_month;
static int32_t sample_count;
static const float *ocean_depth_coordinate;
static int32_t ocean_depth_count;

/*
 * 把 writer 返回的 NetCDF 错误转成运行时错误信息，并统一返回 I/O 失败状态。
 * Translate a writer NetCDF error into a runtime error message and return the
 * common I/O-failure status.
 */
static int report_netcdf_error(const monthly_netcdf_error *error)
{
    runtime_set_error(
        "NetCDF %s failed for %s: %s", error->operation, error->variable,
        monthly_netcdf_nc_error(error->nc_status));
    return PLASIC_RUNTIME_IO_FAILED;
}

/*
 * 读取一个变量的原始诊断值，并完成必要的单位换算与符号翻转。
 * Read the raw diagnostic value of one variable, applying the necessary unit
 * conversion and sign flip.
 *
 * 深海温度/盐度可能不可用，此时返回 0 而不报错。
 * Deep temperature/salinity may be unavailable, in which case 0 is returned
 * instead of raising an error.
 *
 * 约定了 CF 标准的方向：水通量换算为 kg m-2 s-1，辐射和湍流通量统一为
 * “向上为正”。
 * CF conventions are applied: water fluxes are converted to kg m-2 s-1, and
 * radiation and turbulent fluxes are made positive upward.
 */
static float source_value(
    monthly_source source, size_t horizontal, int32_t level,
    const float *sea_surface_temperature, const float *deep_temperature,
    const float *deep_salinity)
{
    /* cell 是三维场索引；surface 指向存放地表值的附加层。 */
    /* cell indexes a 3-D field; surface points to the extra layer holding surface values. */
    const size_t cell = horizontal + (size_t)level * PLASIC_NHOR;
    const size_t surface = horizontal + (size_t)PLASIC_NLEV * PLASIC_NHOR;

    switch (source)
    {
    /* 大气三维场：cell 为 (level, horizontal) 对应的格点索引。 */
    /* Atmospheric 3-D fields: cell is the grid index for (level, horizontal). */
    case SOURCE_UA:
        return plasic_du[cell];
    case SOURCE_VA:
        return plasic_dv[cell];
    case SOURCE_TA:
        return plasic_dt[cell];
    case SOURCE_HUS:
        return plasic_dq[cell];
    case SOURCE_CL:
        return plasic_dcc[cell];
    /* 地表场：ps/tas 为二维诊断量，ts 取温度场附加层中的地表值。 */
    /* Surface fields: ps and tas are 2-D diagnostics, while ts reads the
     * surface value from the extra layer of the temperature field. */
    case SOURCE_PS:
        return plasic_dp[horizontal];
    case SOURCE_TS:
        return plasic_dt[surface];
    case SOURCE_TAS:
        return plasic_dtsa[horizontal];
    /* 海洋与海冰表面场。 */
    /* Ocean and sea-ice surface fields. */
    case SOURCE_TOS:
        return sea_surface_temperature[horizontal];
    case SOURCE_SICONC:
        return plasic_dicec[horizontal];
    case SOURCE_SITHICK:
        return plasic_diced[horizontal];
    /* 陆面场：积雪与土壤水储量。 */
    /* Land-surface fields: snow and soil water storage. */
    case SOURCE_SND:
        return plasic_dsnow[horizontal];
    case SOURCE_MRSO:
        return plasic_dwatc[horizontal];
    /* 地表通量：降水率（m/s 等效水柱）乘以 1000 得 kg m-2 s-1；
     * 蒸发和湍流热通量在模式中以向下为正，取负号后输出为向上为正。 */
    /* Surface fluxes: precipitation rate (m/s water equivalent) times 1000
     * gives kg m-2 s-1; evaporation and turbulent heat fluxes are positive
     * downward in the model and are negated to be positive upward. */
    case SOURCE_PR:
        return 1000.0f * (plasic_dprl[horizontal] + plasic_dprc[horizontal]);
    case SOURCE_EVSPSBL:
        return -1000.0f * plasic_devap[horizontal];
    case SOURCE_HFSS:
        return -plasic_dshfl[horizontal];
    case SOURCE_HFLS:
        return -plasic_dlhfl[horizontal];
    /* 辐射通量：模式以向下为正、向上的通量存为负值；这里按 CF 约定输出，
     * 向下分量保持向下为正，向上分量翻转为向上为正。 */
    /* Radiation fluxes: the model stores positive downward with upward fluxes
     * negative; following CF conventions the outputs keep downwelling
     * components positive downward and flip upwelling components to positive
     * upward. */
    case SOURCE_RSDT:
        return plasic_dfd[horizontal];
    case SOURCE_RSUT:
        return -plasic_dfu[horizontal];
    case SOURCE_RLUT:
        return -plasic_dftu[horizontal];
    case SOURCE_RSDS:
        return plasic_dfd[surface];
    case SOURCE_RSUS:
        return -plasic_dfu[surface];
    case SOURCE_RLDS:
        return plasic_dftd[surface];
    case SOURCE_RLUS:
        return -plasic_dftu[surface];
    case SOURCE_ALBEDO:
        return plasic_dalb[horizontal];
    case SOURCE_MLD:
        return plasic_dmld[horizontal];
    /* 深海温度/盐度来自深海模块，可能不可用；不可用时按 0 处理。 */
    /* Deep temperature/salinity come from the deep-ocean module and may be
     * unavailable; they contribute 0 in that case. */
    case SOURCE_THETAO:
        return deep_temperature == NULL ? 0.0f : deep_temperature[cell];
    case SOURCE_SO:
        return deep_salinity == NULL ? 0.0f : deep_salinity[cell];
    /* CO2 体积混合比由物理模块提供，单位为 ppmv（1e-6）。 */
    /* CO2 volume mixing ratio is provided by the physics module in ppmv
     * (1e-6). */
    case SOURCE_CO2:
        return plasic_physics_co2_field()[horizontal];
    }
    return 0.0f;
}

/*
 * 按需创建全部月度输出文件；NetCDF 创建只在根进程执行。
 * Create all monthly output files on demand; NetCDF creation runs on the root
 * rank only.
 *
 * 变量层数大于 1 时才选择垂直坐标；深海变量在深海层数与坐标长度一致时使用
 * 深度轴，否则回退到大气 sigma 层。
 * A vertical coordinate is needed only when the variable has more than one
 * layer; deep-ocean variables use the depth axis when its length matches,
 * and the atmospheric sigma levels otherwise.
 */
static int ensure_variable_files(void)
{
    int status = PLASIC_RUNTIME_OK;

    if (files_created)
    {
        /* 文件只需创建一次；后续月份直接复用已打开的句柄。 */
        /* Files are created only once; later months reuse the open handles. */
        return PLASIC_RUNTIME_OK;
    }

    if (mp_is_root())
    {
        /* 根进程串行创建全部变量文件；任何失败都会跳出循环并广播给其他进程。 */
        /* The root creates every variable file sequentially; any failure
         * breaks the loop and is broadcast to the other ranks. */
        for (int index = 0; index < MONTHLY_VARIABLE_COUNT; ++index)
        {
            monthly_variable *variable = &variables[index];
            const monthly_definition *definition = variable->definition;
            monthly_netcdf_spec spec;
            monthly_netcdf_error error;
            const float *levels = NULL;
            int writer_status;

            /* 用静态定义填充 writer 所需的形状与 CF 元数据。 */
            /* Fill the shape and CF metadata the writer needs from the static
             * definition. */
            spec.name = definition->name;
            spec.standard_name = definition->standard_name;
            spec.long_name = definition->long_name;
            spec.units = definition->units;
            spec.layers = definition->layers;
            spec.depth_axis = is_depth_source(definition->source);
            /* 三维场才需要垂向坐标：深海变量在深度轴可用且层数匹配时使用
             * 深度轴，否则回退到大气 sigma 层。 */
            /* A vertical coordinate is needed only for 3-D fields: deep-ocean
             * variables use the depth axis when it is available and its length
             * matches, and otherwise fall back to the atmospheric sigma
             * levels. */
            if (definition->layers > 1)
            {
                levels = is_depth_source(definition->source) &&
                                 ocean_depth_coordinate != NULL &&
                                 ocean_depth_count == definition->layers
                             ? ocean_depth_coordinate
                             : plasic_sigma;
            }
            /* 每个变量对应一个 CF/CMIP 风格文件；writer 只依赖传入的坐标与元数据。 */
            /* One CF/CMIP-style file per variable; the writer depends only on
             * the coordinates and metadata passed in. */
            writer_status = monthly_netcdf_create(
                &variable->file, &spec, monthly_directory,
                runtime_co2_forcing_source(), runtime_ozone_forcing_source(),
                calmod_cf_calendar_name(plasic_calendar),
                monthly_latitudes, PLASIC_NLAT,
                monthly_longitudes, PLASIC_NLON,
                levels, definition->layers, &error);
            // signature: int monthly_netcdf_create(
            //     monthly_netcdf_file *file, const monthly_netcdf_spec *spec,
            //     const char *directory, const char *co2_forcing_source,
            //     const char *ozone_forcing_source, const char *calendar,
            //     const float *latitudes, int32_t latitude_count,
            //     const float *longitudes, int32_t longitude_count,
            //     const float *levels, int32_t level_count,
            //     monthly_netcdf_error *error)

            if (writer_status == MONTHLY_NETCDF_PATH_TOO_LONG)
            {
                runtime_set_error("monthly NetCDF path is too long");
                status = PLASIC_RUNTIME_IO_FAILED;
                break;
            }
            if (writer_status != MONTHLY_NETCDF_OK)
            {
                status = report_netcdf_error(&error);
                break;
            }
        }
    }

    /* 所有进程必须就文件是否创建成功达成一致。 */
    /* Every rank must agree on whether file creation succeeded. */
    status = runtime_broadcast_status(status);
    if (status == PLASIC_RUNTIME_OK)
    {
        files_created = 1;
    }
    return status;
}

/*
 * 结束当前月份的累加并写出月平均记录。
 * Finish the current month's accumulation and write its monthly-mean record.
 *
 * 只有当累加步数恰好等于 当月天数 × 每天的步数 时才写文件；否则认为该月不完整
 * （例如运行在月中开始或结束），清空累加和并等待下一个月。
 * A record is written only when the sample count equals days-in-month times
 * steps-per-day; otherwise the month is considered incomplete (for example,
 * the run started or ended mid-month), the sums are discarded and the next
 * month begins.
 */
static int flush_month(void)
{
    /* 时间轴单位是“天”，起点 1850-01-01；base 把第 0 年起算的天数平移到该历元。 */
    /* The time axis is in days since 1850-01-01; base shifts year-0 day counts
     * onto that epoch. */
    const int64_t base = days_before_year(plasic_calendar, 1850);
    const int64_t begin =
        days_before_month(plasic_calendar, active_year, active_month) - base;
    const int64_t end =
        begin + days_in_month(plasic_calendar, active_year, active_month);
    /* 记录时刻取月区间中点，bounds 保存月首、月末，便于按 CF 约定描述月平均。 */
    /* The record time is the midpoint of the month, and bounds hold the month
     * start and end so the monthly mean is CF-compliant. */
    const double time_value = 0.5 * ((double)begin + (double)end);
    const double bounds[2] = {(double)begin, (double)end};
    int status = PLASIC_RUNTIME_OK;

    /* 没有任何样本（例如上一次 flush 刚清空）时无事可做。 */
    /* Nothing to do when no sample was accumulated (for example, right after
     * a previous flush). */
    if (sample_count <= 0)
    {
        return PLASIC_RUNTIME_OK;
    }
    if (sample_count !=
        days_in_month(plasic_calendar, active_year, active_month) *
            plasic_ntspd)
    {
        /* 月份不完整：丢弃累加和，不产生文件记录。 */
        /* Incomplete month: discard the sums and emit no record. */
        for (int index = 0; index < MONTHLY_VARIABLE_COUNT; ++index)
        {
            const size_t local_count =
                (size_t)PLASIC_NHOR *
                (size_t)variables[index].definition->layers;

            memset(variables[index].sum, 0,
                   local_count * sizeof(*variables[index].sum));
        }
        sample_count = 0;
        return PLASIC_RUNTIME_OK;
    }

    status = ensure_variable_files();
    if (status != PLASIC_RUNTIME_OK)
    {
        return status;
    }

    for (int index = 0; index < MONTHLY_VARIABLE_COUNT; ++index)
    {
        /* 逐变量独立处理：先本地求平均、再全局聚合、最后根进程落盘。 */
        /* Handle one variable at a time: local average, then global gather,
         * then root writes the record. */
        monthly_variable *variable = &variables[index];
        const size_t local_count =
            (size_t)PLASIC_NHOR * (size_t)variable->definition->layers;
        const size_t global_count =
            (size_t)PLASIC_NUGP * (size_t)variable->definition->layers;

        /* 本地累加和除以样本数得到该进程负责列的月平均。 */
        /* Divide the local sums by the sample count to get this rank's monthly
         * column means. */
        for (size_t cell = 0; cell < local_count; ++cell)
        {
            variable->sum[cell] /= (float)sample_count;
        }

        /* 把分布在各进程上的列拼成完整的全球经纬网格。 */
        /* Assemble the columns distributed over the ranks into a complete global latitude-longitude grid. */
        if (mp_allgather_grid(
                global_buffer, variable->sum, PLASIC_NUGP,
                PLASIC_NHOR, variable->definition->layers) != 0)
        {
            runtime_set_error("could not gather monthly %s", variable->definition->name);
            return PLASIC_RUNTIME_NUMERICAL_FAILURE;
        }

        if (mp_is_root())
        {
            /* 掩膜值必须与 writer 写入 _FillValue 属性的取值一致（1e20）。 */
            /* The fill value must match the _FillValue attribute defined by
             * the writer (1e20). */
            const float fill_value = 1.0e20f;
            monthly_netcdf_error error;
            int writer_status;

            /* 按变量策略把不属于该表面类型的格点填成 _FillValue；掩膜按水平
             * 列查找，因此同一列的所有层会被一起掩掉。 */
            /* Apply the variable's masking policy by filling cells that do not
             * belong to the requested surface type; the mask is looked up per
             * horizontal column, so all layers of a column are masked
             * together. */
            if (variable->definition->mask != MASK_NONE)
            {
                for (size_t cell = 0; cell < global_count; ++cell)
                {
                    const size_t horizontal = cell % PLASIC_NUGP;
                    const int land = global_land_mask[horizontal] >= 0.5f;
                    if ((variable->definition->mask == MASK_OCEAN && land) ||
                        (variable->definition->mask == MASK_LAND && !land))
                    {
                        global_buffer[cell] = fill_value;
                    }
                }
            }
            writer_status = monthly_netcdf_write_record(
                &variable->file, time_value, bounds, global_buffer, &error);
            if (writer_status != MONTHLY_NETCDF_OK)
            {
                status = report_netcdf_error(&error);
            }
            else
            {
                status = PLASIC_RUNTIME_OK;
            }
        }

        /* 写盘结果由根进程统一广播，保证所有进程看到同样的成败。 */
        /* The write result is broadcast from the root so every rank observes the same outcome. */
        status = runtime_broadcast_status(status);
        if (status != PLASIC_RUNTIME_OK)
        {
            return status;
        }
        /* 该变量已写出，清空累加和以便下一个月份使用。 */
        /* The variable has been written; clear its sums for the next month. */
        memset(variable->sum, 0, local_count * sizeof(*variable->sum));
    }
    sample_count = 0;
    return PLASIC_RUNTIME_OK;
}

/*
 * 启用月度 NetCDF 输出并准备全部运行时资源。
 * Enable monthly NetCDF output and prepare all runtime resources.
 */
int runtime_monthly_netcdf_open(const char *directory)
{
    float allocation_failure = 0.0f;
    int status = PLASIC_RUNTIME_OK;
    struct stat information;

    /* 先清掉上一次运行可能遗留的状态与资源。 */
    /* First clear any state and resources left over from a previous run. */
    runtime_monthly_netcdf_abort();
    if (directory == NULL)
    {
        return PLASIC_RUNTIME_OK;
    }

    /* 只在根进程检查路径长度和目录是否存在，再把结论广播出去。 */
    /* Only the root checks the path length and directory existence, then
     * broadcasts the verdict. */
    if (mp_is_root() &&
        (strlen(directory) >= sizeof(monthly_directory) ||
         stat(directory, &information) != 0 ||
         !S_ISDIR(information.st_mode)))
    {
        runtime_set_error("monthly output directory does not exist: %s", directory);
        status = PLASIC_RUNTIME_IO_FAILED;
    }

    status = runtime_broadcast_status(status);
    if (status != PLASIC_RUNTIME_OK)
    {
        return status;
    }
    /* 路径检查通过后复制目录名；后续 writer 用它拼出 <name>.nc。 */
    /* After the path check passes, copy the directory; the writer later builds
     * <name>.nc from it. */
    (void)snprintf(
        monthly_directory, sizeof(monthly_directory), "%s", directory);

    /* 高斯格点存的是纬度正弦值，这里还原为度数坐标。 */
    /* The Gaussian grid stores the sine of latitude; convert it back to degrees for the coordinate. */
    for (int32_t latitude = 0; latitude < PLASIC_NLAT; ++latitude)
    {
        monthly_latitudes[latitude] = (float)(asin(plasic_sid[latitude]) * 180.0 / 3.14159265358979323846);
    }

    /* 经度按 0..360 均匀分布。 */
    /* Longitudes are uniform over 0..360. */
    for (int32_t longitude = 0; longitude < PLASIC_NLON; ++longitude)
    {
        monthly_longitudes[longitude] =
            360.0f * (float)longitude / (float)PLASIC_NLON;
    }

    /* 深海变量需要的垂向坐标，由物理模块提供。 */
    /* The vertical coordinate needed by deep-ocean variables, provided by the physics module. */
    ocean_depth_coordinate = plasic_physics_ocean_depth(&ocean_depth_count);

    /* global_buffer 按最深的变量分配，供所有变量的 gather 复用；
     * global_land_mask 保存随后聚合得到的陆海掩膜。 */
    /* global_buffer is sized for the deepest variable and reused by every
     * gather; global_land_mask holds the land mask gathered below. */
    global_buffer = (float *)malloc(
        (size_t)PLASIC_NUGP *
        (size_t)(PLASIC_NLEV > PLASIC_DEEP_NLEV ? PLASIC_NLEV
                                                : PLASIC_DEEP_NLEV) *
        sizeof(*global_buffer));
    global_land_mask = (float *)malloc(
        (size_t)PLASIC_NUGP * sizeof(*global_land_mask));

    /* 初始化每个变量的定义、根进程专用句柄与本地累加和。 */
    /* Initialize each variable's definition, root-only handle and local sums. */
    for (int index = 0; index < MONTHLY_VARIABLE_COUNT; ++index)
    {
        variables[index].definition = &definitions[index];
        variables[index].file.ncid = -1;
        variables[index].sum = (float *)calloc(
            (size_t)PLASIC_NHOR * (size_t)definitions[index].layers,
            sizeof(*variables[index].sum));
        if (variables[index].sum == NULL)
        {
            allocation_failure = 1.0f;
        }
    }
    if (global_buffer == NULL || global_land_mask == NULL)
    {
        allocation_failure = 1.0f;
    }

    /* 用一次全局归约确认所有进程分配成功，避免部分进程继续运行。 */
    /* One global reduction confirms that every rank allocated successfully, preventing split execution. */
    if (mp_allreduce_real(&allocation_failure, 1) != 0 ||
        allocation_failure != 0.0f)
    {
        runtime_set_error("could not allocate monthly NetCDF buffers");
        runtime_monthly_netcdf_abort();
        return PLASIC_RUNTIME_INITIALIZATION_FAILED;
    }
    /* 提前把陆地掩膜聚合成全局场，供 flush 的掩膜步骤使用。 */
    /* Gather the land mask once so the masking step of flush can use it. */
    if (mp_allgather_grid(
            global_land_mask, plasic_dls, PLASIC_NUGP,
            PLASIC_NHOR, 1) != 0)
    {
        runtime_set_error("could not gather land mask for monthly output");
        runtime_monthly_netcdf_abort();
        return PLASIC_RUNTIME_NUMERICAL_FAILURE;
    }
    /* 一切就绪，标记输出已启用；active_year/month 保持 -1 直到第一个样本到来。 */
    /* Everything is ready, so mark output as enabled; active_year/month stay
     * -1 until the first sample arrives. */
    monthly_enabled = 1;
    files_created = 0;
    active_year = -1;
    active_month = -1;
    sample_count = 0;
    return PLASIC_RUNTIME_OK;
}

/*
 * 累加一个时间步的瞬时诊断量；模型日期跨月时先写出上一个月。
 * Accumulate the instantaneous diagnostics of one time step; the previous
 * month is flushed first when the model date rolls over.
 */
int runtime_monthly_netcdf_accumulate(int year, int month)
{
    const float *sea_surface_temperature = NULL;
    const float *deep_temperature = NULL;
    const float *deep_salinity = NULL;
    int32_t layers = 0;
    int32_t deep_layers = 0;
    int status;

    if (!monthly_enabled)
    {
        return PLASIC_RUNTIME_OK;
    }
    /* 日期发生变化：先把上一个月的平均写出去。 */
    /* The date changed: write out the previous month's mean first. */
    if (active_year >= 0 && (year != active_year || month != active_month))
    {
        status = flush_month();
        if (status != PLASIC_RUNTIME_OK)
        {
            return status;
        }
    }
    /* 该月的第一个样本决定累加窗口属于哪个月。 */
    /* The first sample of a month fixes which year-month the window covers. */
    if (sample_count == 0)
    {
        active_year = year;
        active_month = month;
    }
    status = plasic_physics_frame_field(
        PLASIC_STREAM_SEA_SURFACE_TEMPERATURE,
        &sea_surface_temperature, &layers);
    if (status != PLASIC_RUNTIME_OK || sea_surface_temperature == NULL)
    {
        runtime_set_error("sea surface temperature is unavailable for monthly output");
        return PLASIC_RUNTIME_INITIALIZATION_FAILED;
    }
    /* 深海场允许不可用；此时 source_value 会按 0 处理。 */
    /* Deep fields may be unavailable; source_value then contributes 0. */
    (void)plasic_physics_frame_field(
        PLASIC_STREAM_DEEP_OCEAN_TEMPERATURE, &deep_temperature,
        &deep_layers);
    (void)plasic_physics_frame_field(
        PLASIC_STREAM_DEEP_OCEAN_SALINITY, &deep_salinity, &deep_layers);

    /* 把本步所有变量累加到各自的 sum 数组。 */
    /* Accumulate every variable of this step into its own sum array. */
    for (int index = 0; index < MONTHLY_VARIABLE_COUNT; ++index)
    {
        monthly_variable *variable = &variables[index];

        /* 遍历该变量的全部垂直层，再遍历本进程负责的水平列。 */
        /* Walk every vertical layer, then the local horizontal columns owned by
         * this rank. */
        for (int32_t level = 0; level < variable->definition->layers; ++level)
        {
            for (size_t horizontal = 0; horizontal < PLASIC_NHOR; ++horizontal)
            {
                const size_t cell = horizontal + (size_t)level * PLASIC_NHOR;
                variable->sum[cell] += source_value(
                    variable->definition->source, horizontal, level,
                    sea_surface_temperature, deep_temperature,
                    deep_salinity);
            }
        }
    }
    /* 本步已计入样本数；flush 时用它判断月份是否完整。 */
    /* The step is now counted; flush uses the total to decide whether the
     * month is complete. */
    ++sample_count;
    return PLASIC_RUNTIME_OK;
}

/*
 * 结束月度输出：写出最后一个月、关闭全部文件并释放资源；可重复调用。
 * Finish monthly output: flush the final month, close every file and release
 * resources; safe to call repeatedly.
 */
int runtime_monthly_netcdf_close(void)
{
    int status = PLASIC_RUNTIME_OK;

    /* 若仍处于启用状态，先把最后一个已完整累加的月份写出去。 */
    /* While still enabled, flush the last month if it was fully accumulated. */
    if (monthly_enabled)
    {
        status = flush_month();
    }

    for (int index = 0; index < MONTHLY_VARIABLE_COUNT; ++index)
    {
        monthly_variable *variable = &variables[index];

        /* 只有根进程持有有效句柄；非根进程始终是 -1。 */
        /* Only the root holds a valid handle; it is always -1 elsewhere. */
        if (mp_is_root() && variable->file.ncid >= 0)
        {
            monthly_netcdf_error error;
            const int close_status =
                monthly_netcdf_close(&variable->file, &error);

            if (close_status != MONTHLY_NETCDF_OK && status == PLASIC_RUNTIME_OK)
            {
                status = report_netcdf_error(&error);
            }
        }
        variable->file.ncid = -1;
        free(variable->sum);
        variable->sum = NULL;
    }
    /* 句柄与缓冲区全部释放后复位全局状态，保证可重复调用。 */
    /* Reset global state once every handle and buffer has been released, so
     * the call is repeatable. */
    free(global_buffer);
    global_buffer = NULL;
    free(global_land_mask);
    global_land_mask = NULL;
    monthly_enabled = 0;
    files_created = 0;
    sample_count = 0;
    monthly_directory[0] = '\0';
    /* 关闭阶段的错误（若有）也广播给所有进程，保持返回码一致。 */
    /* Broadcast any close-stage error so every rank returns the same code. */
    return runtime_broadcast_status(status);
}

/*
 * 中止月度输出：不做任何写盘，直接关闭句柄并释放全部资源。
 * 可在 open() 失败、初始化失败或进程退出时重复调用；从未启用时是空操作。
 * Abort monthly output: close handles and release every resource without
 * writing anything. It may be called repeatedly after an open failure, an
 * initialization failure or at process teardown, and is a no-op when output
 * was never enabled.
 */
void runtime_monthly_netcdf_abort(void)
{
    /* 逐个变量释放句柄与累加和；非根进程的句柄本来就是 -1。 */
    /* Release each variable's handle and sums; handles are already -1 on
     * non-root ranks. */
    for (int index = 0; index < MONTHLY_VARIABLE_COUNT; ++index)
    {
        monthly_variable *variable = &variables[index];

        if (variable->file.ncid >= 0)
        {
            monthly_netcdf_abort(&variable->file);
        }
        variable->file.ncid = -1;
        free(variable->sum);
        variable->sum = NULL;
    }
    free(global_buffer);
    global_buffer = NULL;
    free(global_land_mask);
    global_land_mask = NULL;
    monthly_enabled = 0;
    files_created = 0;
    sample_count = 0;
    monthly_directory[0] = '\0';
}
