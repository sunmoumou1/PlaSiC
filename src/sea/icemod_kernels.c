// ============================================================================
// 热力海冰数值核（零层热力学海冰）。
//
// 变量名与教程 4.3 的两张示意图（sea_ice_algorithm.svg、
// sea_ice_flux_structure.svg）中的符号一一对应：
//
//   符号        含义                                       代码名
//   ---------   ----------------------------------------   ---------------------------
//   T_i         冰雪表面（表皮）温度                         surface_temperature
//   T_i^n       步初表面温度，先截到不高于 T_m               previous_surface_temperature
//   T_i^*       半隐式试算表面温度                           trial_surface_temperature
//   T_o         slab 混合层温度（SST）                       sst
//   T_m, T_f    融点 273.16 K、海水冻结点 271.25 K           melting_temperature,
//                                                            freezing_temperature
//   h_i, h_s    冰厚、几何雪深 h_s = (rho_w/rho_s) S         ice_thickness, snow_depth
//   S           冰上积雪（m 水当量）                         snow
//   P_s         降雪率（m 水当量 / s）                       snowfall
//   k_bar       冰雪层厚度加权导热率                         mean_conductivity
//   K_z         冰雪柱等效导热 k_bar/(h_i+h_s)               effective_conductance
//   C_i         表皮热容量 rho_i c_p,i h_min                 skin_heat_capacity
//   Q_a^mod     修正大气净热通量 Q_a - Q_max                 modified_atmospheric_heat
//   Q_m         表面达到 T_m 后的融化能量                    melt_energy
//   Q_m^rem     融雪后剩下的 Q_m                             remaining_melt_energy
//   Q_c         冰雪柱导热；开阔水 Q_a^mod，薄冰 0           conductive_flux
//   Q_f         冻结点导热 K_z (T_f - T_o)                   freezing_flux
//   Q_s->i      积雪压沉转冰（flooding）通量                 snow_to_ice_flux
//   Q_o->i      slab 上一步返回给海冰的通量                  ocean_to_ice_flux
//   Q_net       到达冰的全部热量（冰能量账的输入）           net_heat_flux
//   Q_ice       冰能量账：Q_net + Q_s->i - Q_f               ice_energy_flux
//   Q_ex        超过现有冰潜热的融化能量                     excess_melt_flux
//   Q_i->o      交给 slab 的通量                             ice_to_ocean_flux
//   c_i^*       连续紧密度（预报量）                         compactness
//   c_i         0/1 海冰掩膜                                 ice_cover
//
// 通量单位均为 W m^-2。Q_a、Q_c、Q_net、Q_ice、Q_i->o 向下为正；
// Q_m、Q_ice、Q_o->i 为正时融冰，为负时长冰。
// rho_w = 1000 kg m^-3 是把雪水当量换成质量的淡水密度，代码中直接写作
// 1000.0f；seawater_density（rho_o = 1030）只用于浮力和 flooding。
// ============================================================================
#include "icemod_kernels.h"

#include <math.h>
#include <stddef.h>

// 月气候态的布局为 field[month][horizontal]；同一个月的 nhor 个格点连续。
// 用 size_t 完成乘法，避免大网格下 int32_t 索引溢出。
static size_t monthly_index(int32_t horizontal, int32_t month, int32_t nhor)
{
    return (size_t)horizontal + (size_t)month * (size_t)nhor;
}

// ----------------------------------------------------------------------------
// 冰能量账（示意图第 5 步）：按潜热守恒更新 h_i，并得到交给 slab 的 Q_i->o。
//
// 进入本函数时：
//   net_heat_flux     = Q_c - rho_w L_s P_s + Q_m^rem（由 ice_snow 累加）；
//   freezing_flux     = Q_f = K_z (T_f - T_o)；
//   ocean_to_ice_flux = Q_o->i，slab 上一步返回；
//   snow_to_ice_flux  = Q_s->i <= 0。
//
// 本函数先把 Q_o->i 并入 Q_net，然后：
//   - 若 T_o <= T_f 或已有冰：Q_ice = Q_net + Q_s->i - Q_f，
//     h_i <- h_i - Q_ice dt / (rho_i L_i)，Q_i->o = Q_f；
//   - 否则（无冰且 T_o > T_f 的开阔水）：不改变冰，Q_i->o = Q_net。
// 若 h_i 试算 <= 0，超额融化能量 Q_ex = -h_i rho_i L_i / dt 退回 Q_i->o，
// 而不是被截断丢失。正 Q_ice 使冰变薄，负 Q_ice 使冰增长。
// ----------------------------------------------------------------------------
void ice_make_ice(
    int32_t nhor, float freezing_temperature, float ice_density,
    float latent_heat_ice, float timestep, const float *land_mask,
    const float *ocean_to_ice_flux, const float *sst,
    const float *snow_to_ice_flux, const float *freezing_flux,
    float *net_heat_flux, float *ice_energy_flux, float *excess_melt_flux,
    float *ice_to_ocean_flux, float *ice_thickness)
{
    /*
     * 每 1 m 海冰在一个时间步内对应的潜热通量：
     *     ice_energy_per_timestep = rho_i L_i / dt.
     * 因此，后面用“相变通量 / ice_energy_per_timestep”就可以得到
     * 本步的冰厚变化（单位：m）。
     */
    const float ice_energy_per_timestep = ice_density * latent_heat_ice / timestep;

    /* 每个水平格点都是一个彼此独立的冰雪柱。 */
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        /*
         * Q_o->i 是 slab 的冻结点约束返回给冰的通量。它只存在于海洋
         * 格点，先并入 Q_net，之后和其它冰相变能量统一结算。
         */
        if (land_mask[horizontal] < 0.5f)
        {
            net_heat_flux[horizontal] =
                net_heat_flux[horizontal] + ocean_to_ice_flux[horizontal];
        }

        /*
         * 只有以下两种情况允许改变冰的潜热储量：
         *
         *   1. 海水已经达到冻结点（T_o <= T_f），可以生成新冰；
         *   2. 已经存在海冰（h_i > 0），可以融冰或继续增厚。
         *
         * 无冰且 SST 高于冻结点时，热量只能留给海洋显热，不能用于
         * 凭空生成海冰：此时 Q_net 全部交给 slab。
         */
        if (sst[horizontal] <= freezing_temperature ||
            ice_thickness[horizontal] > 0.0f)
        {
            /* 冰能量账：Q_ice = Q_net + Q_s->i - Q_f。 */
            ice_energy_flux[horizontal] = net_heat_flux[horizontal] +
                                          snow_to_ice_flux[horizontal] -
                                          freezing_flux[horizontal];

            /* 潜热守恒：h_i <- h_i - Q_ice dt / (rho_i L_i)。 */
            ice_thickness[horizontal] = ice_thickness[horizontal] -
                                        ice_energy_flux[horizontal] /
                                            ice_energy_per_timestep;

            /*
             * Q_c - Q_f = K_z (T_i - T_f) 留在冰里；冰底以下的 Q_f 是
             * 冰与 slab 之间的接口通量，交给 slab。
             */
            ice_to_ocean_flux[horizontal] = freezing_flux[horizontal];
        }
        else
        {
            /* 温暖的开阔水：不结算冰，Q_i->o = Q_net。 */
            ice_to_ocean_flux[horizontal] = net_heat_flux[horizontal];
        }

        /*
         * 如果试算冰厚不为正，说明融化能量超过了现有冰的全部潜热。
         * h_i 不能变成负数；超出的能量 Q_ex 必须退回海洋预算，而不能丢失。
         */
        if (ice_thickness[horizontal] <= 0.0f)
        {
            /*
             * 下面三行保持原有的浮点运算和符号关系（h_i <= 0）：
             *   - ice_to_ocean_flux：增加退回 slab 的 Q_ex；
             *   - ice_energy_flux：从 Q_ice 中扣除没能用上的 Q_ex；
             *   - excess_melt_flux：记录 Q_ex。
             */
            ice_to_ocean_flux[horizontal] = ice_to_ocean_flux[horizontal] -
                                            ice_thickness[horizontal] *
                                                ice_energy_per_timestep;
            ice_energy_flux[horizontal] = ice_energy_flux[horizontal] +
                                          ice_thickness[horizontal] *
                                              ice_energy_per_timestep;
            excess_melt_flux[horizontal] = excess_melt_flux[horizontal] -
                                           ice_thickness[horizontal] *
                                               ice_energy_per_timestep;
            ice_thickness[horizontal] = 0.0f;
        }

        /* 避免浮点舍入留下会误触发“已有冰”分支的幽灵冰。 */
        if (fabsf(ice_thickness[horizontal]) < 1.0e-8f)
        {
            ice_thickness[horizontal] = 0.0f;
        }
    }
}

// ----------------------------------------------------------------------------
// 将厚度变化映射为连续紧密度 c_i^*（示意图第 6 步）。
// 它用于低成本地表示冰缘扩张/退缩；最终传给 atmosphere/sea 的 ice_cover
// 会在调用方按阈值 c_thr 二值化，不能解释为格点内的真实部分覆盖面积。
// ----------------------------------------------------------------------------
void ice_compactness(
    int32_t nhor, const float *old_thickness, const float *new_thickness,
    float *compactness)
{
    // h_r = 0.5 m：经验厚度尺度，控制新冰向未覆盖面积扩张的速度。
    const float reference_thickness = 0.5f;

    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        // h_i 增加时，只让尚未覆盖的 1-c_i^* 部分扩张。
        if (new_thickness[horizontal] > old_thickness[horizontal])
        {
            compactness[horizontal] = compactness[horizontal] +
                                      (1.0f - compactness[horizontal]) *
                                          (new_thickness[horizontal] -
                                           old_thickness[horizontal]) /
                                          reference_thickness;
            compactness[horizontal] =
                fminf(compactness[horizontal], 1.0f);
        }

        // h_i 减少时，c_i^* 按相对厚度损失衰减；此分支只有 old_thickness > 0
        // 时才有物理意义，调用顺序保证该条件。
        if (new_thickness[horizontal] < old_thickness[horizontal])
        {
            compactness[horizontal] = compactness[horizontal] +
                                      compactness[horizontal] *
                                          (new_thickness[horizontal] -
                                           old_thickness[horizontal]) /
                                          (2.0f * old_thickness[horizontal]);
            compactness[horizontal] =
                fmaxf(compactness[horizontal], 0.0f);
        }
        // 冰的潜热储库耗尽时，连续紧密度也必须归零。
        if (new_thickness[horizontal] <= 0.0f)
        {
            compactness[horizontal] = 0.0f;
        }
    }
}

// ----------------------------------------------------------------------------
// 积雪（示意图第 4 步）：更新冰上积雪 S，执行“先融雪、再融冰”和
// 雪转冰两条能量路径，并把相应项累加进 net_heat_flux（Q_net）。
//
// 调用前 net_heat_flux 已等于 Q_c（薄冰为 0，其热量全部在 Q_m 中）。本函数依次：
//   1. c_i = 0 处的降雪落在开阔水面：Q_net -= rho_w L_s P_s；
//   2. 有冰（h_i > 0）时，Q_m 先融雪：
//        M_s = max[0, min(S/dt, Q_m/(rho_w L_s))]，S -= M_s dt，
//        Q_m^rem = Q_m - rho_w L_s M_s，Q_net += Q_m^rem；
//      无冰时残留的 S 全部融化：Q_net -= rho_w L_s S/dt；
//   3. 若 rho_w S > (rho_o - rho_i) h_i（雪把冰面压到海面以下），把多余
//      的雪转成冰，得到 Q_s->i <= 0；冰厚本身在 ice_make_ice 中统一结算。
// ----------------------------------------------------------------------------
void ice_snow(
    int32_t nhor, float minimum_ice_cover, float seawater_density,
    float ice_density, float latent_heat_snow, float timestep,
    const float *land_mask, const float *ice_cover, float *snowfall,
    const float *melt_energy, const float *ice_thickness, float *snow,
    float *net_heat_flux, float *snow_melt, float *snow_melt_flux,
    float *snow_to_ice_flux, float *snow_to_ice_rate)
{
    // 所有过程均为本地一维冰雪柱过程，格点之间没有雪或海冰输送。
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        // c_i 未达到保雪阈值时，把降雪视为落在开阔水面：不建立 S，
        // 而是立即从热量预算扣除其融化潜热 rho_w L_s P_s。
        if (ice_cover[horizontal] < minimum_ice_cover &&
            snowfall[horizontal] > 0.0f)
        {
            // 1000 kg m^-3 将雪水当量转换为面质量；负号表示这份能量用于融雪。
            net_heat_flux[horizontal] = net_heat_flux[horizontal] -
                                        snowfall[horizontal] * 1000.0f *
                                            latent_heat_snow;
            snow_melt[horizontal] = snowfall[horizontal] + snow_melt[horizontal];
            snow_melt_flux[horizontal] = snow_melt_flux[horizontal] +
                                         snowfall[horizontal] * 1000.0f *
                                             latent_heat_snow;
            snowfall[horizontal] = 0.0f;
        }

        // 存在冰时：Q_m 先改变 S，剩余部分 Q_m^rem 才能改变 h_i。
        if (ice_thickness[horizontal] > 0.0f)
        {
            // rho_w L_s S/dt：一步内融尽全部积雪所需的通量。
            const float melt_all_snow_flux = snow[horizontal] * 1000.0f *
                                             latent_heat_snow / timestep;
            const int32_t melts_all_snow =
                melt_energy[horizontal] > melt_all_snow_flux;
            const int32_t melts_some_snow =
                !melts_all_snow && melt_energy[horizontal] > 0.0f &&
                snow[horizontal] > 0.0f;

            // Q_m^rem = Q_m - rho_w L_s M_s：融尽积雪时为余量，部分融雪时为零，
            // 不融雪（Q_m <= 0，或无雪）时就是 Q_m 本身。
            const float remaining_melt_energy =
                melts_all_snow
                    ? melt_energy[horizontal] - melt_all_snow_flux
                    : (melts_some_snow ? 0.0f : melt_energy[horizontal]);

            if (melts_all_snow)
            {
                // M_s = S/dt
                snow_melt[horizontal] = snow_melt[horizontal] +
                                        snow[horizontal] / timestep;
                snow_melt_flux[horizontal] = snow_melt_flux[horizontal] +
                                             snow[horizontal] * 1000.0f *
                                                 latent_heat_snow / timestep;
                snow[horizontal] = 0.0f;
            }
            else if (melts_some_snow)
            {
                // M_s = Q_m/(rho_w L_s)
                snow_melt[horizontal] = snow_melt[horizontal] +
                                        melt_energy[horizontal] /
                                            latent_heat_snow / 1000.0f;
                snow[horizontal] = snow[horizontal] -
                                   melt_energy[horizontal] * timestep /
                                       latent_heat_snow / 1000.0f;
                snow_melt_flux[horizontal] = snow_melt_flux[horizontal] +
                                             melt_energy[horizontal];
            }
            // Q_net += Q_m^rem
            net_heat_flux[horizontal] = net_heat_flux[horizontal] +
                                        remaining_melt_energy;
        }
        // 无冰时，积雪没有可依附的冰；残余 S 全部融化，热量取自 Q_net。
        else
        {
            snow_melt[horizontal] = snow_melt[horizontal] +
                                    snow[horizontal] / timestep;
            net_heat_flux[horizontal] = net_heat_flux[horizontal] -
                                        snow[horizontal] * 1000.0f *
                                            latent_heat_snow / timestep;
            snow_melt_flux[horizontal] = snow_melt_flux[horizontal] +
                                         snow[horizontal] * 1000.0f *
                                             latent_heat_snow / timestep;
            snow[horizontal] = 0.0f;
        }

        // snow_to_ice_rate 是质量转换诊断；每步重新计算，不能沿用旧值。
        snow_to_ice_rate[horizontal] = 0.0f;
        // 由阿基米德浮力条件 rho_w S > (rho_o - rho_i) h_i
        // 判断积雪是否压低冰雪界面并触发雪转冰（flooding）。
        if (land_mask[horizontal] < 0.5f &&
            1000.0f * snow[horizontal] >
                (seawater_density - ice_density) * ice_thickness[horizontal])
        {
            // 目标厚度 h_target = (rho_w S + rho_i h_i)/rho_o，
            // 即转换后冰雪界面恰回到海平面。
            const float target_thickness =
                (1000.0f * snow[horizontal] +
                 ice_density * ice_thickness[horizontal]) /
                seawater_density;
            // 把同样的质量从 S 移到 h_i：dS = -(rho_i/rho_w) dh <= 0，dh >= 0。
            const float flooding_thickness_change =
                target_thickness - ice_thickness[horizontal];
            const float flooding_snow_change =
                -flooding_thickness_change * ice_density / 1000.0f;

            // 此处只更新雪；冰的潜热变化必须和其他 Q_ice 项一起在 ice_make_ice 中结算。
            snow[horizontal] = snow[horizontal] + flooding_snow_change;

            // Q_s->i = rho_w L_s dS / dt <= 0（负值使 Q_ice 变小，即长冰）。
            snow_to_ice_flux[horizontal] =
                flooding_snow_change * 1000.0f * latent_heat_snow / timestep;
            snow_to_ice_rate[horizontal] =
                flooding_thickness_change * ice_density / 1000.0f / timestep;
            // 与 PlaSim 相同，转成冰的雪也记入融雪诊断
            // （flooding_snow_change <= 0，下面两行实际是加）。
            snow_melt[horizontal] = snow_melt[horizontal] -
                                    flooding_snow_change / timestep;
            snow_melt_flux[horizontal] = snow_melt_flux[horizontal] -
                                         flooding_snow_change * 1000.0f *
                                             latent_heat_snow / timestep;
        }
    }
}

// ----------------------------------------------------------------------------
// 冰雪柱导热（示意图第 3 步）：诊断 Q_c 和冻结点导热 Q_f。
//
// 雪水当量先转换为几何雪深 h_s = rho_w S/rho_s，再使用
// 厚度加权导热率 k_bar 和 K_z = k_bar/(h_i+h_s)。对 h_i >= h_min：
//
//   conductive_flux = diagnosed_conductive_flux = Q_c = K_z (T_i - T_o)
//   freezing_flux                                = Q_f = K_z (T_f - T_o)
//
// 两者之差 Q_c - Q_f = K_z (T_i - T_f) 是冰底（T_f）与表面之间的导热，
// 留在冰能量账里；Q_f 是冰底与 slab 之间的部分，交给 slab。
// h_i < h_min 时不解析冰雪柱，Q_f = 0：开阔水（h_i = 0）Q_c = Q_a^mod；
// 薄冰（0 < h_i < h_min）Q_c = 0，因为 Q_a^mod 已在第 2 步作为 Q_m 进入 Q_net。
// PlaSim 的 mkcflux 在薄冰上同样令 xcflux = xheat，使 Q_a^mod 被算两次；
// 这里有意不沿用。
// h_i >= h_max 时 K_z = 0（关闭垂直导热，冰不再因导热增厚）。
// ----------------------------------------------------------------------------
void ice_conductive_flux(
    int32_t nhor, int32_t ocean_flux_correction, int32_t prognostic_ice,
    float minimum_thickness, float maximum_thickness,
    float freezing_temperature, float snow_density, float ice_conductivity,
    float snow_conductivity, const float *ice_thickness, const float *snow,
    const float *surface_temperature, const float *sst,
    const float *modified_atmospheric_heat,
    float *diagnosed_conductive_flux, float *conductive_flux,
    float *freezing_flux)
{
    // 逐格点计算等效导热系数及相应通量。
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        diagnosed_conductive_flux[horizontal] = 0.0f;
        // h_i < h_min 时没有独立、可解析的冰雪柱，Q_f = 0。
        //   开阔水（h_i = 0）：Q_c = Q_a^mod，大气热量经 Q_net 交给 slab，
        //     或在 T_o <= T_f 时直接结冰；
        //   薄冰（0 < h_i < h_min）：Q_a^mod 已在第 2 步作为 Q_m 进入，由
        //     ice_snow 先融雪、余量 Q_m^rem 进入 Q_net，故 Q_c = 0，
        //     否则 Q_a^mod 会在 Q_net 中算两次。
        if (ice_thickness[horizontal] < minimum_thickness)
        {
            conductive_flux[horizontal] =
                ice_thickness[horizontal] > 0.0f
                    ? 0.0f
                    : modified_atmospheric_heat[horizontal];
            freezing_flux[horizontal] = 0.0f;
        }
        // 厚冰：以串联冰雪层的等效导热参数化连接 T_i 和 T_o。
        else
        {
            // h_s = (rho_w/rho_s) S
            const float snow_depth = 1000.0f / snow_density * snow[horizontal];
            // k_bar = (k_s h_s + k_i h_i)/(h_s + h_i)
            const float mean_conductivity =
                (snow_conductivity * snow_depth +
                 ice_conductivity * ice_thickness[horizontal]) /
                (snow_depth + ice_thickness[horizontal]);
            // k_bar，或 h_i >= h_max 时为 0；K_z = effective_conductivity/(h_s+h_i)。
            const float effective_conductivity =
                maximum_thickness >= 0.0f &&
                        ice_thickness[horizontal] >= maximum_thickness
                    ? 0.0f
                    : mean_conductivity;
            // Q_c = K_z (T_i - T_o)
            diagnosed_conductive_flux[horizontal] =
                effective_conductivity *
                (surface_temperature[horizontal] - sst[horizontal]) /
                (snow_depth + ice_thickness[horizontal]);
            conductive_flux[horizontal] = diagnosed_conductive_flux[horizontal];
            // Q_f = K_z (T_f - T_o)
            freezing_flux[horizontal] =
                effective_conductivity * (freezing_temperature - sst[horizontal]) /
                (snow_depth + ice_thickness[horizontal]);
        }
    }

    // 气候态冰路径或海洋 SST 在线订正路径不使用冻结点导热，
    // 故显式关掉 Q_f。该开关与 OCEANMOD 共享。
    if (ocean_flux_correction == 1 || prognostic_ice == 0)
    {
        for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
        {
            freezing_flux[horizontal] = 0.0f;
        }
    }
}


// ----------------------------------------------------------------------------
// 表面温度（示意图第 2 步）：更新冰雪表面温度 T_i，并把到达融点后不能
// 继续升温的能量转为 Q_m。
//
// 对 h_i >= h_min 的厚冰，表层有效热容量取 C_i = rho_i c_p,i h_min，
// 并半隐式求解表面能量平衡
//
//   (C_i/dt)(T_i^* - T_i^n) = Q_a^mod - K_z (T_i^* - T_o)
//   => T_i^* = [(C_i/dt) T_i^n + Q_a^mod + K_z T_o] / [C_i/dt + K_z]
//
// 若 T_i^* > T_m，则 T_i^(n+1) = T_m，并输出融化能量
//   Q_m = Q_a^mod + K_z (T_o - T_m) - (C_i/dt)(T_m - T_i^n) > 0；
// 否则 T_i^(n+1) = T_i^*，Q_m = 0。
// 薄冰/无冰不建立独立表层储热：T_i^(n+1) = T_o；若仍有薄冰（0 < h_i < h_min），
// Q_m = Q_a^mod（正负均可）。这是薄冰上 Q_a^mod 进入冰能量账的唯一途径，
// ice_conductive_flux 相应令薄冰 Q_c = 0。
// ----------------------------------------------------------------------------
void ice_skin_temperature(
    int32_t nhor, float minimum_thickness, float maximum_thickness,
    float melting_temperature, float ice_density, float ice_heat_capacity,
    float snow_density, float ice_conductivity, float snow_conductivity,
    float timestep, const float *ice_thickness, const float *snow,
    const float *sst, const float *modified_atmospheric_heat,
    float *surface_temperature, float *surface_storage_flux,
    float *melt_energy)
{
    // C_i：零层模型只给最薄的有效表皮（h_min）配置热容量，而非求解冰内温度剖面。
    const float skin_heat_capacity =
        minimum_thickness * ice_density * ice_heat_capacity;

    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        surface_storage_flux[horizontal] = 0.0f;
        melt_energy[horizontal] = 0.0f;

        // T_i^n 不允许高于 T_m；超过融点的旧值没有可保留的冰雪显热意义。
        const float previous_surface_temperature =
            fminf(melting_temperature, surface_temperature[horizontal]);
        // h_s = (rho_w/rho_s) S
        const float snow_depth = 1000.0f / snow_density * snow[horizontal];

        // 厚冰
        if (ice_thickness[horizontal] >= minimum_thickness)
        {
            // C_i/dt
            const float skin_heat_capacity_per_timestep =
                skin_heat_capacity / timestep;
            // k_bar = (k_s h_s + k_i h_i)/(h_s + h_i)
            const float mean_conductivity =
                (snow_conductivity * snow_depth +
                 ice_conductivity * ice_thickness[horizontal]) /
                (snow_depth + ice_thickness[horizontal]);

            // k_bar/(h_i + h_s)，尚未考虑 h_max 截断
            const float unlimited_conductance =
                mean_conductivity /
                (ice_thickness[horizontal] + snow_depth);

            // K_z：h_i >= h_max 时为 0
            const float effective_conductance =
                ice_thickness[horizontal] >= maximum_thickness &&
                        maximum_thickness >= 0.0f
                    ? 0.0f
                    : unlimited_conductance;

            // 分子中的 Q_a^mod + K_z T_o。
            const float surface_heat_input =
                modified_atmospheric_heat[horizontal] +
                effective_conductance * sst[horizontal];

            // T_i^* = [(C_i/dt) T_i^n + Q_a^mod + K_z T_o] / [C_i/dt + K_z]
            const float trial_surface_temperature =
                (skin_heat_capacity_per_timestep * previous_surface_temperature +
                 surface_heat_input) /
                (skin_heat_capacity_per_timestep + effective_conductance);
            surface_temperature[horizontal] = trial_surface_temperature;

            // T_i 不能高于 T_m：剩余能量按先融雪、后融冰的次序处理。
            if (trial_surface_temperature > melting_temperature)
            {
                // Q_m = Q_a^mod + K_z(T_o - T_m) - (C_i/dt)(T_m - T_i^n)
                melt_energy[horizontal] =
                    modified_atmospheric_heat[horizontal] +
                    effective_conductance *
                        (sst[horizontal] - melting_temperature) -
                    skin_heat_capacity_per_timestep *
                        (melting_temperature - previous_surface_temperature);
                surface_temperature[horizontal] = melting_temperature;
            }
            // (C_i/dt)(T_i^(n+1) - T_i^n)：本步留在冰雪表皮的显热通量（诊断）。
            surface_storage_flux[horizontal] =
                skin_heat_capacity_per_timestep *
                (surface_temperature[horizontal] - previous_surface_temperature);
        }
        // 薄冰/无冰：海水就是大气看到的表面，不设独立冰表温度。
        else
        {
            surface_temperature[horizontal] = sst[horizontal];
            // 虽未达到可解析厚冰阈值，残余薄冰仍以全部 Q_a^mod 作为 Q_m。
            if (ice_thickness[horizontal] > 0.0f)
            {
                melt_energy[horizontal] = modified_atmospheric_heat[horizontal];
            }
        }
    }
}


// ----------------------------------------------------------------------------
// 在两个相邻月气候态间插值。
// ----------------------------------------------------------------------------
void ice_interpolate_climatology(
    int32_t nhor, int32_t previous_month1, int32_t previous_month2,
    float previous_weight1, float previous_weight2, int32_t month1,
    int32_t month2, float weight1, float weight2,
    const float *climatological_sst, const float *climatological_cover,
    const float *climatological_thickness, float *previous_sst,
    float *current_sst, float *current_cover, float *current_thickness)
{
    /*
     * 这个函数只负责“查表并插值”，不推进海冰状态：
     *
     *   - climatological_* 是按月份保存的固定气候态数据；
     *   - month1/month2 是当前时刻两侧的月份；
     *   - weight1/weight2 是两个月份的插值权重，通常满足
     *     weight1 + weight2 = 1；
     *   - previous_* 表示前一个模型时刻使用的月份和权重。
     *
     * 月份数组的内存布局是 field[month][horizontal]
     */
    for (int32_t horizontal = 0; horizontal < nhor; ++horizontal)
    {
        /*
         * 保存前一个模型时刻的气候态海表温度。
         * 它主要用于需要比较相邻模型时刻的季节性过程。
         */
        previous_sst[horizontal] =
            previous_weight1 *
                climatological_sst[monthly_index(horizontal, previous_month1, nhor)] +
            previous_weight2 *
                climatological_sst[monthly_index(horizontal, previous_month2, nhor)];

        /* 当前模型时刻的气候态海表温度。 */
        current_sst[horizontal] =
            weight1 * climatological_sst[monthly_index(horizontal, month1, nhor)] +
            weight2 * climatological_sst[monthly_index(horizontal, month2, nhor)];

        /* 当前模型时刻的气候态海冰厚度。 */
        current_thickness[horizontal] =
            weight1 *
                climatological_thickness[monthly_index(horizontal, month1, nhor)] +
            weight2 *
                climatological_thickness[monthly_index(horizontal, month2, nhor)];

        /* 当前模型时刻的气候态海冰覆盖率。 */
        current_cover[horizontal] =
            weight1 * climatological_cover[monthly_index(horizontal, month1, nhor)] +
            weight2 * climatological_cover[monthly_index(horizontal, month2, nhor)];

        // 数值防护：厚度和覆盖率是非负状态量。
        current_cover[horizontal] =
            fmaxf(current_cover[horizontal], 0.0f);
        current_thickness[horizontal] =
            fmaxf(current_thickness[horizontal], 0.0f);

        /* 没有冰厚时，覆盖率也必须为零，保持两个状态量的一致性。 */
        // 没有潜热储库就不能声明存在海冰覆盖。
        if (current_thickness[horizontal] <= 0.0f)
        {
            current_cover[horizontal] = 0.0f;
        }
    }
}


// ----------------------------------------------------------------------------
// 从月气候态覆盖率构造经验海冰厚度，供非预报冰或冷启动初始化使用。
//
// 这里的 thickness 不是通过海冰潜热方程预报出来的厚度，而是根据给定的
// 连续覆盖率 compactness 诊断出的“参考厚度”。北、南两个纬向半球使用不同
// 的经验映射；北侧还乘以一个逐月修正因子，以表示季节性厚度差异。
//
// 三个二维/三维量的布局都是：
//
//     field[month][latitude][longitude]
//
// northern_factor 有 14 个元素，包含月气候态使用的年界 halo；调用方必须
// 保证 month 索引位于该数组范围内。
// ----------------------------------------------------------------------------
void ice_make_thickness(
    int32_t nlon, int32_t nlat, int32_t nmonths,
    const float *compactness, float *thickness)
{
    /* 北侧海冰厚度的逐月经验修正；最后两个值对应年界 halo 月份。 */
    static const float northern_factor[14] = {
        0.912f, 0.942f, 1.0f, 1.058f, 1.124f, 1.161f, 1.175f,
        1.058f, 0.931f, 0.883f, 0.88f, 0.876f, 0.912f, 0.942f};
    /*
     * 北侧经验关系：
     *   cover < 0.1       ：只有少量冰，赋予最小厚度；
     *   0.1 <= cover < 0.9：在最小和最大厚度之间线性插值；
     *   cover >= 0.9      ：赋予最大厚度。
     */
    const float northern_minimum_cover = 0.1f;
    const float northern_maximum_cover = 0.9f;
    const float northern_minimum_thickness = 0.25f;
    const float northern_maximum_thickness = 3.0f;
    /*
     * 南侧经验关系：
     *   cover < 0.25       ：厚度为 0.25 + cover；
     *   cover >= 0.25      ：厚度固定为 0.50 m。
     * 南侧不使用月修正。
     */
    const float southern_minimum_cover = 0.25f;
    const float southern_minimum_thickness = 0.25f;
    const float southern_maximum_thickness = 0.50f;

    /* 一个“月份 × 纬度 × 经度”平面的格点数。 */
    const size_t plane = (size_t)nlon * (size_t)nlat;

    /*
     * 先清零全部输出。这样 cover <= 0 的格点无需显式赋值，仍保持无冰；
     * 同时也保证输出数组中不会残留调用前的旧厚度。
     */
    for (size_t index = 0; index < plane * (size_t)nmonths; ++index)
    {
        thickness[index] = 0.0f;
    }

    /* 按 month → latitude → longitude 遍历，最内层 longitude 与内存布局一致。 */
    for (int32_t month = 0; month < nmonths; ++month)
    {
        for (int32_t latitude = 0; latitude < nlat; ++latitude)
        {
            for (int32_t longitude = 0; longitude < nlon; ++longitude)
            {
                /* 将三维下标转换为 field[month][latitude][longitude] 的一维下标。 */
                const size_t index = (size_t)longitude +
                                     (size_t)latitude * (size_t)nlon +
                                     (size_t)month * plane;
                const float cover = compactness[index];
                
                /* 数组前半纬向带采用北侧经验关系；真实南北方向由网格纬度排序决定。 */
                if (latitude < nlat / 2)
                {
                    /* 覆盖率很高时直接采用北侧最大厚度，再乘季节因子。 */
                    if (cover >= northern_maximum_cover)
                    {
                        thickness[index] =
                            northern_maximum_thickness * northern_factor[month];
                    }

                    /* 中等覆盖率：在最小厚度和最大厚度之间线性插值。 */
                    else if (cover >= northern_minimum_cover)
                    {
                        thickness[index] =
                            (northern_minimum_thickness +
                             (northern_maximum_thickness -
                              northern_minimum_thickness) *
                                 (cover - northern_minimum_cover) /
                                 (northern_maximum_cover -
                                  northern_minimum_cover)) *
                            northern_factor[month];
                    }
                    
                    /* 只要存在少量海冰，就至少赋予北侧最小厚度。 */
                    else if (cover > 0.0f)
                    {
                        thickness[index] =
                            northern_minimum_thickness * northern_factor[month];
                    }
                }
                /* 后半纬向带采用较薄、无逐月因子的南侧经验关系。 */
                else
                {
                    /* 南侧覆盖率达到阈值后使用固定的最大经验厚度。 */
                    if (cover >= southern_minimum_cover)
                    {
                        thickness[index] = southern_maximum_thickness;
                    }
                    /* 南侧低覆盖率时，厚度按 0.25 + cover 增加。 */
                    else if (cover > 0.0f)
                    {
                        thickness[index] =
                            southern_minimum_thickness + cover;
                    }
                }
            }
        }
    }
}
