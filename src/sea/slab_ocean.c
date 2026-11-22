#include "slab_ocean.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

/*
 * 单层 slab 混合层海洋的数值核。
 *
 * 参数名与教程 4.3 / 4.4 示意图中的符号对应：
 *
 *   符号      含义                                   代码名
 *   -------   ------------------------------------   ---------------------------
 *   T_o       混合层（slab）温度，即 SST             sst
 *   T_o^*     本步试算温度（双精度工作量）           trial_sst
 *   T_f       海水冻结点 271.25 K                    freezing_temperature
 *   H         混合层厚度（默认 50 m）                mixed_layer_depth
 *   dt_mix    slab 时间步 = N 个海冰步               mixed_layer_timestep
 *   rho_o     海水密度 1030 kg m^-3                  seawater_density
 *   c_p,o     海水比热 4180 J kg^-1 K^-1             seawater_heat_capacity
 *   C_o       单位面积热容量 rho_o c_p,o H           （由上面三项组合）
 *   Q_i->o    海冰交给 slab 的通量                   ice_to_ocean_flux
 *   Q_deep    3-D 深海交给 slab 的通量               deep_ocean_flux
 *   Q_o->i    slab 返回给海冰的通量                  -ice_flux
 *
 * 通量单位 W m^-2，Q_i->o、Q_deep 为正时使 slab 升温；ice_flux 使用海洋侧
 * 符号（为正时表示海洋向冰索取热量），ocean_copy_outputs() 取负后作为
 * Q_o->i 交给海冰。
 */

/*
 * 初始化单层混合层海洋的参数和时间步长。
 *
 * 参数说明：
 *   nhor                 : 水平网格点数。
 *   steps_per_day        : 每个太阳日内的 slab 时间步数。
 *   solar_day            : 一个太阳日对应的模型时间长度（s）。
 *   restoring_timescale  : 输入/输出；输入是以“天”为单位的 Newtonian
 *                          restoring 时间常数，输出转换为秒。
 *   layer_depth          : 混合层厚度 H，对所有水平网格点相同。
 *   mixed_layer_depth    : 输出；将 layer_depth 扩展到所有水平网格点。
 *   mixed_layer_timestep : 输出；一个 slab 时间步的长度（s）。
 */
void ocean_initialize_parameters(
    int32_t nhor, int32_t steps_per_day, float solar_day,
    float *restoring_timescale, float layer_depth,
    float *mixed_layer_depth, float *mixed_layer_timestep)
{
    /*
     * 单层混合层：层厚 layer_depth 对所有水平点相同，
     * 因此直接复制成混合层厚度场 mixed_layer_depth(horizontal)。
     */
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        mixed_layer_depth[horizontal] = layer_depth;
    }

    /* 一个太阳日均分为 steps_per_day 个 slab 时间步。 */
    *mixed_layer_timestep = solar_day / (float)steps_per_day;

    /* 将以“太阳日”为单位给出的时间常数转换为模型时间单位。 */
    const float restoring_timescale_in_days = *restoring_timescale;
    const float restoring_timescale_in_seconds =
        solar_day * restoring_timescale_in_days;
    *restoring_timescale = restoring_timescale_in_seconds;
}

/*
 * 在两个相邻月的水平场之间做线性插值。
 *
 * 参数说明：
 *   month1/month2 : 两个时间切片的零基索引。
 *   weight        : 第二个时间切片的权重；第一个切片权重为 1 - weight。
 *   monthly       : 按 [month][horizontal] 连续存储的输入场。
 *   field         : 输出的当前时刻水平场。
 *
 * 当 weight = 0 时完全取 month1；weight = 1 时完全取 month2。
 * 函数本身不限制 weight 范围，因此调用端也可以用它进行线性外推。
 */
void ocean_interpolate_cycle(
    int32_t nhor, int32_t month1, int32_t month2, float weight,
    const float *monthly, float *field)
{
    /* 两个月份在连续数组中的起始偏移。 */
    const size_t first = (size_t)month1 * (size_t)nhor;
    const size_t second = (size_t)month2 * (size_t)nhor;

    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        /* 标准线性插值：(1-w) * x1 + w * x2。 */
        field[horizontal] =
            (1.0f - weight) * monthly[first + (size_t)horizontal] + weight * monthly[second + (size_t)horizontal];
    }
}

/*
 * 把 slab 的结果交给海冰模块（示意图第 7 步的出口）。
 *
 * 参数说明：
 *   sst/mixed_layer_depth : slab 的 T_o 和 H。
 *   ice_flux              : slab 内部的冰通量（海洋侧符号）。
 *   ice_sst/ice_mixed_layer_depth/ocean_to_ice_flux : 海冰模块中的
 *                           ice.sst、ice.mixed_layer_depth、ice.ocean_to_ice_flux。
 *
 * 注意：Q_o->i = -ice_flux。海冰在下一步的冰能量账中使用它
 * （ice_make_ice 把它并入 Q_net），因此 Q_o->i 比 Q_i->o 滞后一个海冰步。
 */
void ocean_copy_outputs(
    int32_t nhor, const float *sst, const float *mixed_layer_depth,
    const float *ice_flux, float *ice_sst, float *ice_mixed_layer_depth,
    float *ocean_to_ice_flux)
{
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        ice_sst[horizontal] = sst[horizontal];
        ice_mixed_layer_depth[horizontal] = mixed_layer_depth[horizontal];

        /* 海洋侧 ice_flux 取负，得到海冰侧 Q_o->i。 */
        ocean_to_ice_flux[horizontal] = -ice_flux[horizontal];
    }
}

/*
 * 将海冰模块在耦合窗口内平均好的强迫复制到 slab 的工作数组
 * （示意图第 7 步的入口），并构造冰上总雪/冰水当量。
 *
 * input_ice_to_ocean_flux 是 ice.ocean_heat_accumulator 的窗口平均，
 * 即 slab 看到的 Q_i->o。
 * ice_and_snow 把冰厚乘以冰密度后除以 1000，再与雪量相加；
 * 若 ice_density 以 kg/m^3、厚度以 m 表示，则该项对应以水密度
 * 1000 kg/m^3 归一化的水当量。
 */
void ocean_copy_inputs(
    int32_t nhor, float ice_density, const float *input_ice_cover,
    const float *input_ice_thickness, const float *input_ice_to_ocean_flux,
    const float *input_freshwater_flux, const float *input_stress_x,
    const float *input_stress_y, const float *input_friction_velocity_cubed,
    const float *input_snow, const float *input_climatological_ice_thickness,
    float *ice_cover, float *ice_thickness, float *ice_to_ocean_flux,
    float *precipitation_minus_evaporation, float *stress_x,
    float *stress_y, float *friction_velocity_cubed, float *ice_and_snow,
    float *climatological_ice_thickness)
{
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        /* 以下为逐水平点的一一复制，不做缩放或筛选。 */
        ice_cover[horizontal] = input_ice_cover[horizontal];
        ice_thickness[horizontal] = input_ice_thickness[horizontal];
        ice_to_ocean_flux[horizontal] = input_ice_to_ocean_flux[horizontal];
        precipitation_minus_evaporation[horizontal] =
            input_freshwater_flux[horizontal];
        stress_x[horizontal] = input_stress_x[horizontal];
        stress_y[horizontal] = input_stress_y[horizontal];
        friction_velocity_cubed[horizontal] =
            input_friction_velocity_cubed[horizontal];

        /* 冰的水当量与积雪量相加，得到冰雪总量相关内部场。 */
        ice_and_snow[horizontal] = input_ice_thickness[horizontal] * ice_density / 1000.0f + input_snow[horizontal];

        climatological_ice_thickness[horizontal] =
            input_climatological_ice_thickness[horizontal];
    }
}

/*
 * slab 时间步的前半部分：试算 T_o^*，并建立冰通量 ice_flux。
 *
 * 判断“有冰/无冰”使用的是本步海冰更新之后的 h_i（ice_thickness）。
 *
 *   无冰海点（h_i <= 0）：T_o^* = T_o^n + (dt_mix/C_o)(Q_i->o + Q_deep)；
 *   有冰海点（h_i > 0） ：T_o^* = T_o^n，ice_flux = -(Q_i->o + Q_deep)，
 *                         即这部分热量不改变 T_o，而是交给冰。
 *
 * 随后若步初温度 T_o^n < T_f（注意检查的是旧值 sst，而不是 T_o^*），
 * 把 T_o^* 拉回 T_f，所需热量 (C_o/dt_mix)(T_f - T_o^*) 记入 ice_flux。
 *
 * 参数说明：
 *   mixed_layer_timestep   : dt_mix。
 *   seawater_density       : rho_o。
 *   seawater_heat_capacity : c_p,o。
 *   freezing_temperature   : T_f。
 *   land_mask              : 海陆掩膜，land_mask < 1 表示海点。
 *   ice_thickness          : 本步更新后的 h_i；<= 0 视为无冰。
 *   ice_to_ocean_flux      : Q_i->o。
 *   deep_ocean_flux        : Q_deep（3-D 深海在 slab 之前算出；slab-only
 *                            构建中为 0）。
 *   mixed_layer_depth      : H。
 *   sst                    : 步初温度 T_o^n（只读）。
 *   trial_sst              : 输出；T_o^*（双精度）。
 *   ice_flux               : 输出；冰通量（海洋侧符号）。
 */
void ocean_mksst_begin(
    int32_t nhor, float mixed_layer_timestep, float seawater_density,
    float seawater_heat_capacity, float freezing_temperature,
    const float *land_mask, const float *ice_thickness,
    const float *ice_to_ocean_flux, const float *deep_ocean_flux,
    const float *mixed_layer_depth, const float *sst, double *trial_sst,
    float *ice_flux)
{
    /*
     * dt_mix / (rho_o c_p,o)。再除以 H 就是 dt_mix / C_o，
     * 即把通量换成本步温度增量的系数。
     */
    const float temperature_response_factor = mixed_layer_timestep / (seawater_density * seawater_heat_capacity);

    /* 处理每个水平点的表层热收支和初始冰通量。 */
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        /* T_o^* 从步初温度 T_o^n 开始，用双精度工作变量累加。 */
        trial_sst[horizontal] = (double)sst[horizontal];

        /* 每个时间步开始时先清零冰通量工作量。 */
        ice_flux[horizontal] = 0.0f;

        /*
         * 无冰海点：T_o^* = T_o^n + (Q_i->o + Q_deep) dt_mix / C_o。
         */
        if (land_mask[horizontal] < 1.0f && ice_thickness[horizontal] <= 0.0f)
        {
            const float increment =
                (ice_to_ocean_flux[horizontal] + deep_ocean_flux[horizontal]) *
                temperature_response_factor / mixed_layer_depth[horizontal];
            trial_sst[horizontal] += (double)increment;
        }

        /*
         * 有冰海点：T_o 不变，Q_i->o + Q_deep 以相反符号记入 ice_flux，
         * 经 Q_o->i = -ice_flux 交还给冰。
         */
        if (land_mask[horizontal] < 1.0f && ice_thickness[horizontal] > 0.0f)
        {
            ice_flux[horizontal] = -deep_ocean_flux[horizontal] - ice_to_ocean_flux[horizontal];
        }
    }

    /*
     * 保留原始 MKIFLUX 的条件语义：条件检查的是步初温度 sst（T_o^n），
     * 而不是已更新的 trial_sst。这一区别会影响本步刚因热通量越过冻结点
     * 的格点（它们要到下一步才被拉回 T_f），review 时不可擅自替换。
     */
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        if (sst[horizontal] < freezing_temperature && land_mask[horizontal] < 1.0f)
        {
            /*
             * 把 T_o^* 调整到 T_f 所需的热量换算为通量并累加到 ice_flux：
             *   (C_o/dt_mix)(T_f - T_o^*)
             *     = (T_f - T_o^*) * H / temperature_response_factor。
             * 随后把 T_o^* 设为 T_f。
             */
            ice_flux[horizontal] =
                (float)((double)ice_flux[horizontal] +
                        ((double)freezing_temperature - trial_sst[horizontal]) *
                            (double)mixed_layer_depth[horizontal] /
                            (double)temperature_response_factor);
            trial_sst[horizontal] = (double)freezing_temperature;
        }
    }
}

/*
 * slab 时间步的后半部分：冰下 T_o 固定在 T_f，并写回 T_o。
 *
 * 仅当预报海洋（prognostic_ocean > 0）且不使用 SST restoring
 * （flux_correction == 0，默认）时，对有冰海点补上
 * (C_o/dt_mix)(T_f - T_o^*)，使 T_o^* 回到 T_f，并把这部分热量记入
 * ice_flux。最后 T_o = (float)T_o^*。
 *
 * 参数说明：
 *   prognostic_ocean       : 大于 0 表示使用预报海洋。
 *   flux_correction        : 0 表示关闭 SST restoring；1 表示在线计算
 *                            Newtonian restoring 修正通量。
 *   mixed_layer_timestep   : dt_mix。
 *   seawater_density       : rho_o。
 *   seawater_heat_capacity : c_p,o。
 *   freezing_temperature   : T_f。
 *   land_mask              : 海陆掩膜；小于 1.0 表示海点。
 *   ice_thickness          : 本步更新后的 h_i；大于 0 表示有冰。
 *   mixed_layer_depth      : H。
 *   trial_sst              : 输入/输出；T_o^*（双精度）。
 *   ice_flux               : 输入/输出；冰通量，函数内进行累加。
 *   sst                    : 输出；新的 T_o（单精度）。
 */
void ocean_mksst_finish(
    int32_t nhor, int32_t prognostic_ocean, int32_t flux_correction,
    float mixed_layer_timestep, float seawater_density,
    float seawater_heat_capacity, float freezing_temperature,
    const float *land_mask, const float *ice_thickness,
    const float *mixed_layer_depth, double *trial_sst, float *ice_flux,
    float *sst)
{
    /*
     * dt_mix / (rho_o c_p,o)，再除以 H 即 dt_mix / C_o。
     */
    const float temperature_response_factor = mixed_layer_timestep / (seawater_density * seawater_heat_capacity);
    /*
     * 只有在预报海洋开启且不使用海洋热通量修正时，才在此处显式处理
     * 有海冰格点的冻结点约束。
     */
    if (prognostic_ocean > 0 && flux_correction == 0)
    {
        for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
        {
            /* 仅处理“存在海冰”的海洋格点，陆地格点保持不变。 */
            if (ice_thickness[horizontal] > 0.0f && land_mask[horizontal] < 1.0f)
            {
                /*
                 * 把 T_o^* 调整到 T_f 所需的通量 (C_o/dt_mix)(T_f - T_o^*)；
                 * H / temperature_response_factor 即 C_o / dt_mix。
                 */
                const float flux_to_freezing_temperature =
                    (float)(((double)freezing_temperature - trial_sst[horizontal]) *
                            (double)mixed_layer_depth[horizontal] /
                            (double)temperature_response_factor);
                /*
                 * 再把该通量换算成温度增量。该写法与上面的反推公式配对，
                 * 理想算术下会使 T_o^* 到达 T_f。
                 */
                const float increment =
                    flux_to_freezing_temperature *
                    temperature_response_factor / mixed_layer_depth[horizontal];
                trial_sst[horizontal] += (double)increment;
                /* 这部分热量由冰提供：并入 ice_flux（即 -Q_o->i）。 */
                ice_flux[horizontal] += flux_to_freezing_temperature;
            }
        }
    }
    /*
     * T_o^* 用双精度完成中间计算；时间步结束时写回单精度 T_o。
     */
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        sst[horizontal] = (float)trial_sst[horizontal];
    }
}

/*
 * 可选的 SST restoring（flux_correction == 1；默认关闭）：应用 Newtonian
 * restoring 通量，并在海洋显热与海冰相变热之间重新分配能量。
 *
 * 对每个水平格点，本函数使用同一个混合层热收支关系：
 *
 *   C_o = rho_o c_p,o H
 *   dT  = F dt_mix / C_o
 *
 * 其中 C_o 是单位面积混合层热容量，F 是进入混合层的热通量。
 * 因而，通量与温度变化之间可以互相换算：
 *
 *   temperature increment = F * temperature_response_factor / H
 *   flux needed for dT    = dT * H / temperature_response_factor
 *
 * current_climatological_sst 是插值到当前模型时刻的气候态 SST
 * （调用方传入 ocean.interpolated_sst），作为 restoring 目标 T_clim。
 *
 * 根据海冰状态，本函数有三条路径：
 *
 *   1. 无海冰：
 *      restoring 通量全部改变 SST。
 *
 *   2. 有气候态海冰（climatological_ice_thickness > 0）：
 *      将 T_clim 作为局部目标温度，只施加不会越过该温度的部分；
 *      未用于改变 SST 的通量通过 residual_ice_flux 转交给海冰模块。
 *
 *   3. 有实际海冰但没有气候态海冰：
 *      先将 restoring 通量全部施加到 SST；若 SST 高于海水冻结点，
 *      再计算融化现有海冰所需的潜热，并从海洋中扣除。
 *
 * 第二遍循环还会消耗已有的正 ice_flux，将高于约束温度的 SST 冷却下来；
 * 最后强制 SST 不低于 T_f，并把对应的冻结/结冰能量记入 residual_ice_flux。
 * ice_flux 使用海洋侧符号，之后传给海冰模块时会在 ocean_copy_outputs() 中反号。
 */
static void add_flux_correction(
    int32_t nhor, float mixed_layer_timestep, float restoring_timescale,
    float seawater_density, float ice_density, float seawater_heat_capacity,
    float latent_heat_ice, float freezing_temperature, const float *land_mask,
    const float *current_climatological_sst, const float *mixed_layer_depth,
    const float *ice_thickness, const float *climatological_ice_thickness,
    const float *ice_to_ocean_flux, float *sst, float *ice_flux,
    float *residual_ice_flux)
{
    /*
     * 热通量到温度变化的公共系数 dt_mix / (rho_o c_p,o)；
     * 对单位面积通量 F 还需再除以 H，才能得到 dT。
     */
    const float temperature_response_factor = mixed_layer_timestep / (seawater_density * seawater_heat_capacity);
    /* 体积热容量 rho_o c_p,o，用于从温差反推单位面积热通量。 */
    const float volumetric_heat_capacity = seawater_density * seawater_heat_capacity;
    /*
     * corrected_sst：SST 的双精度工作副本。一个时间步内可能连续执行
     *       restoring、融冰冷却和冻结点截断，因此使用双精度减少累计舍入误差。
     *
     * constraint_temperature：第二遍循环使用的局部温度上限/目标：
     *       - 默认是 T_f；
     *       - 气候态海冰格点改为 T_clim。
     */
    double *corrected_sst = malloc((size_t)nhor * sizeof(*corrected_sst));
    float *constraint_temperature =
        malloc((size_t)nhor * sizeof(*constraint_temperature));

    /* 内存分配失败时清理并终止，避免继续使用空指针。 */
    if (corrected_sst == NULL || constraint_temperature == NULL)
    {
        free(corrected_sst);
        free(constraint_temperature);
        abort();
    }

    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        /* residual_ice_flux 只记录本次修正产生的残余通量，不继承上一个时间步。 */
        residual_ice_flux[horizontal] = 0.0f;
        constraint_temperature[horizontal] = freezing_temperature;
        corrected_sst[horizontal] = (double)sst[horizontal];

        /* 陆地格点不参与海洋热通量修正。 */
        if (land_mask[horizontal] >= 1.0f)
        {
            continue;
        }

        /*
         * 根据当前 SST 偏差即时计算 restoring 通量。
         *
         * restoring_timescale (tau) > 0 时：
         *
         *   F_R = rho_o c_p,o H (T_clim - T_o) / tau
         *
         * 这是方程 dT/dt = (T_clim - T_o)/tau 乘以混合层单位面积
         * 热容量后的结果。正通量使海洋升温，负通量使海洋降温。
         * tau <= 0 时，F_R 取一步内恰好回到 T_clim 所需的通量再扣去 Q_i->o。
         */
        const float restoring_flux =
            restoring_timescale > 0.0f ?
                volumetric_heat_capacity *
                    (current_climatological_sst[horizontal] - sst[horizontal]) *
                    mixed_layer_depth[horizontal] / restoring_timescale :
                volumetric_heat_capacity / mixed_layer_timestep *
                    (current_climatological_sst[horizontal] - sst[horizontal]) *
                    mixed_layer_depth[horizontal] -
                    ice_to_ocean_flux[horizontal];

        /*
         * 无海冰：restoring 通量全部进入海洋显热。
         *
         *   T_o <- T_o + F_R dt_mix / C_o
         */
        if (ice_thickness[horizontal] <= 0.0f)
        {
            const float increment =
                restoring_flux * temperature_response_factor /
                mixed_layer_depth[horizontal];
            corrected_sst[horizontal] += (double)increment;
        }
        /*
         * 气候态海冰：T_clim 作为局部目标温度。只有“朝向目标温度”的
         * restoring 通量才允许改变 SST，并且不能越过目标温度；
         * 剩余通量转为海冰侧的残余热通量。
         */
        else if (climatological_ice_thickness[horizontal] > 0.0f)
        {
            constraint_temperature[horizontal] =
                current_climatological_sst[horizontal];
            /*
             * 当前 SST 高于 T_clim，且 restoring 为冷却方向。
             * 到达目标温度所需的通量为：
             *
             *   F_target = C_o (T_clim - T_o) / dt_mix < 0
             *
             * fmax(F_target, F_R) 使实际冷却量不超过 restoring 提供的量，
             * 从而不会把 SST 冷却到 T_clim 以下。
             */
            if (corrected_sst[horizontal] > (double)current_climatological_sst[horizontal] && restoring_flux < 0.0f)
            {
                const float flux_to_climatology =
                    (float)(((double)current_climatological_sst[horizontal] -
                             corrected_sst[horizontal]) *
                            (double)mixed_layer_depth[horizontal] /
                            (double)temperature_response_factor);
                /* 实际用于改变 SST 的通量：取两者中较大的负值。 */
                const float applied_flux =
                    fmaxf(flux_to_climatology, restoring_flux);
                const float increment =
                    applied_flux * temperature_response_factor /
                    mixed_layer_depth[horizontal];
                corrected_sst[horizontal] += (double)increment;
                /*
                 * 未用于 SST 的 restoring 部分转交海冰：
                 *
                 *   residual_ice_flux = -(F_R - F_applied)
                 *
                 * 这里的负号来自海洋侧 ice_flux 与海冰侧通量的符号约定。
                 */
                residual_ice_flux[horizontal] =
                    residual_ice_flux[horizontal] -
                    (restoring_flux - applied_flux);
            }
            /*
             * 当前 SST 低于 T_clim，且 restoring 为增暖方向。
             * 此时 F_target > 0，使用 fmin(F_target, F_R)，避免升温越过目标。
             */
            else if (corrected_sst[horizontal] < (double)current_climatological_sst[horizontal] && restoring_flux > 0.0f)
            {
                const float flux_to_climatology =
                    (float)(((double)current_climatological_sst[horizontal] -
                             corrected_sst[horizontal]) *
                            (double)mixed_layer_depth[horizontal] /
                            (double)temperature_response_factor);
                /* 实际只使用到达目标温度所需的那一部分正通量。 */
                const float applied_flux =
                    fminf(flux_to_climatology, restoring_flux);
                const float increment =
                    applied_flux * temperature_response_factor /
                    mixed_layer_depth[horizontal];
                corrected_sst[horizontal] += (double)increment;
                residual_ice_flux[horizontal] =
                    residual_ice_flux[horizontal] -
                    (restoring_flux - applied_flux);
            }
            /*
             * 其他方向组合会把 SST 推离目标温度：
             * 不让 restoring 改变 SST，而把全部通量交给海冰侧处理。
             */
            else
            {
                residual_ice_flux[horizontal] -= restoring_flux;
            }
        }
        /*
         * 实际海冰但无气候态海冰：先按普通混合层热收支施加 restoring。
         * 与气候态海冰分支不同，这里不把 T_clim 当作 SST 的硬目标。
         */
        else
        {
            const float increment =
                restoring_flux * temperature_response_factor /
                mixed_layer_depth[horizontal];
            corrected_sst[horizontal] += (double)increment;
            /*
             * 若施加后 SST 高于冻结点，则将当前冰量对应的融化潜热
             * 从海洋显热中扣除：
             *
             *   F_melt = rho_i h_i L_i / dt_mix
             *
             * 这是单位面积、单位时间的等效融冰通量。它转入 residual_ice_flux，
             * 后续传给海冰模块后会改变海冰潜热储量。
             */
            if (corrected_sst[horizontal] > (double)freezing_temperature)
            {
                const float melt_flux =
                    ice_thickness[horizontal] * ice_density * latent_heat_ice / mixed_layer_timestep;
                const float cooling =
                    melt_flux * temperature_response_factor /
                    mixed_layer_depth[horizontal];
                corrected_sst[horizontal] -= (double)cooling;
                residual_ice_flux[horizontal] -= melt_flux;
            }
        }
    }

    /*
     * 第二遍：处理之前已经存在的 ice_flux，并执行最终冻结点约束。
     * 第一遍主要决定 restoring 通量如何分配；第二遍则把已有海冰通量
     * 与 SST 状态重新协调，避免同一份能量同时被计入海洋和海冰。
     */
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        /*
         * 若 SST 高于局部约束温度，且 ice_flux > 0，尝试用已有海冰通量冷却海洋。
         * 到达约束温度所需的通量为：
         *
         *   F_target = C_o (T_constraint - T_o) / dt_mix < 0.
         *
         * 由于 ice_flux 是可用的正通量，实际抵消量不能超过 -ice_flux，
         * 所以使用 max(F_target, -ice_flux)。
         */
        if (corrected_sst[horizontal] > (double)constraint_temperature[horizontal] &&
            ice_flux[horizontal] > 0.0f && land_mask[horizontal] < 1.0f)
        {
            const float flux_to_constraint_temperature =
                (float)(((double)constraint_temperature[horizontal] -
                         corrected_sst[horizontal]) *
                        (double)mixed_layer_depth[horizontal] /
                        (double)temperature_response_factor);
            /* 限制实际冷却量，不能消耗超过已有的 ice_flux。 */
            const float applied_flux =
                fmaxf(flux_to_constraint_temperature,
                      -ice_flux[horizontal]);
            const float increment =
                applied_flux * temperature_response_factor /
                mixed_layer_depth[horizontal];
            corrected_sst[horizontal] += (double)increment;
            ice_flux[horizontal] += applied_flux;
        }
        /*
         * 最终硬性冻结点约束：SST 不能低于 T_f。
         * 若 T_o < T_f，需要补充：
         *
         *   F_freeze = C_o (T_f - T_o) / dt_mix > 0.
         *
         * 这部分能量使海水回到冻结点，并作为结冰/冻结残余通量转交海冰。
         */
        if (corrected_sst[horizontal] < (double)freezing_temperature && land_mask[horizontal] < 1.0f)
        {
            residual_ice_flux[horizontal] =
                (float)((double)residual_ice_flux[horizontal] +
                        ((double)freezing_temperature - corrected_sst[horizontal]) *
                            (double)mixed_layer_depth[horizontal] /
                            (double)temperature_response_factor);
            corrected_sst[horizontal] = (double)freezing_temperature;
        }
        /* 将双精度工作值写回单精度表层 SST。 */
        sst[horizontal] = (float)corrected_sst[horizontal];
    }

    /* 释放本函数申请的临时数组。 */
    free(corrected_sst);
    free(constraint_temperature);
}

/*
 * 完成一个 slab 时间步：生成/应用热通量修正、累积诊断量，并按输出周期
 * 计算时间平均。
 *
 * 返回值：
 *   0 - 当前步不是输出步，诊断量继续累积；
 *   1 - 当前步到达输出周期，诊断累积量已除以样本数得到平均值。
 *
 * 关键控制参数：
 *   prognostic_ocean = 0：使用规定的气候态 SST，不预报海洋温度；
 *   prognostic_ocean > 0：使用预报海洋（默认）。
 *
 *   flux_correction = 0：不使用海洋热通量修正（默认）；
 *   flux_correction = 1：根据气候态 SST 与当前 SST 的偏差在线计算
 *                        Newtonian restoring 修正通量。
 *
 * 诊断累积数组（restart 记录 yheata / yifluxa / yfldoa）：
 *   ice_to_ocean_flux_accumulator - Q_i->o 累计；
 *   ice_flux_accumulator          - ice_flux（-Q_o->i）累计；
 *   deep_ocean_flux_accumulator   - Q_deep 累计。
 *
 * accumulation_count 记录自上次重置以来累计的样本数。
 */
int32_t ocean_finish_step(
    int32_t nhor, int32_t prognostic_ocean, int32_t flux_correction,
    int32_t nstep, int32_t output_interval, float mixed_layer_timestep,
    float restoring_timescale, float seawater_density, float ice_density,
    float seawater_heat_capacity, float latent_heat_ice,
    float freezing_temperature, const float *land_mask,
    const float *current_climatological_sst, const float *mixed_layer_depth,
    const float *ice_thickness, const float *climatological_ice_thickness,
    const float *ice_to_ocean_flux, const float *deep_ocean_flux,
    float *sst, float *ice_flux, float *residual_ice_flux,
    float *ice_to_ocean_flux_accumulator, float *ice_flux_accumulator,
    float *deep_ocean_flux_accumulator, int32_t *accumulation_count)
{
    /*
     * 无预报海洋：直接把 SST 设为气候态值。
     */
    if (prognostic_ocean == 0)
    {
        for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
        {
            if (land_mask[horizontal] < 1.0f)
            {
                sst[horizontal] = current_climatological_sst[horizontal];
            }
        }
    }
    /* 预报海洋且启用方案 1：在线计算并应用恢复通量。 */
    else if (flux_correction == 1)
    {
        add_flux_correction(
            nhor, mixed_layer_timestep, restoring_timescale,
            seawater_density, ice_density, seawater_heat_capacity,
            latent_heat_ice, freezing_temperature, land_mask,
            current_climatological_sst, mixed_layer_depth, ice_thickness,
            climatological_ice_thickness, ice_to_ocean_flux, sst, ice_flux,
            residual_ice_flux);
    }

    /*
     * 启用在线修正方案时，将 add_flux_correction 产生的残余通量
     * 合并到主海冰热通量 ice_flux。
     */
    if (prognostic_ocean > 0 && flux_correction == 1)
    {
        for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
        {
            if (land_mask[horizontal] < 1.0f)
            {
                ice_flux[horizontal] += residual_ice_flux[horizontal];
            }
        }
    }

    /*
     * 累积本时间步的各类诊断量。这里不区分海陆，保持调用方提供数组的
     * 原始约定；陆地值应由上游逻辑保证为合适的值。
     */
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        ice_to_ocean_flux_accumulator[horizontal] += ice_to_ocean_flux[horizontal];
        deep_ocean_flux_accumulator[horizontal] += deep_ocean_flux[horizontal];
        ice_flux_accumulator[horizontal] += ice_flux[horizontal];
    }
    /* 本次时间步已贡献一个诊断样本。 */
    *accumulation_count += 1;

    /*
     * 只有当 nstep 能被 output_interval 整除时才生成输出平均值。
     * 调用方应保证 output_interval 非零。
     */
    if (nstep % output_interval != 0)
    {
        return 0;
    }
    /* 用累计样本数把各诊断累计和转换为时间平均值。 */
    const float count = (float)*accumulation_count;
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        ice_to_ocean_flux_accumulator[horizontal] /= count;
        deep_ocean_flux_accumulator[horizontal] /= count;
        ice_flux_accumulator[horizontal] /= count;
    }
    /* 返回 1 通知调用方：本步已形成可输出的平均诊断量。 */
    return 1;
}

/*
 * 重置 slab 诊断累计量，为下一个输出累计窗口做准备。
 *
 * 所有逐格点累计数组清零，同时把累计样本计数 accumulation_count 置为 0。
 * 通常应在 finish_step 返回 1、诊断量完成输出之后调用。
 */
void ocean_reset_diagnostics(
    int32_t nhor, float *ice_to_ocean_flux_accumulator,
    float *ice_flux_accumulator, float *deep_ocean_flux_accumulator,
    int32_t *accumulation_count)
{
    /* 逐水平格点清空全部诊断累计数组。 */
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        ice_to_ocean_flux_accumulator[horizontal] = 0.0f;
        ice_flux_accumulator[horizontal] = 0.0f;
        deep_ocean_flux_accumulator[horizontal] = 0.0f;
    }
    /* 清空累计样本数，使下一输出窗口从零开始计数。 */
    *accumulation_count = 0;
}
