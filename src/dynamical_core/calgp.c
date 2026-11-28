#include "plasic_kernels.h"

#include <stddef.h>
#include <stdlib.h>

/*
 * calcgp 独立成文件后仍使用文件内静态辅助函数，避免向动力核公开实现细节。
 * 三维场统一按 [level][point] 展平，同一层的水平点连续存放。
 */
static size_t field_index(int32_t point, int32_t level, int32_t point_count)
{
    return (size_t)level * (size_t)point_count + (size_t)point;
}

static float max_zero(float value)
{
    return value > 0.0f ? value : 0.0f;
}

/*
 * 在每个高斯格点的一根 σ 气柱内计算绝热非线性动力。
 *
 * 本文件遵守一个可读性约束：一个变量在整个生命周期内只代表一个物理量。
 * 数组可以持续累加同一种物理量，但不能被当作另一种物理量的临时工作区。
 * 变量名直接说明数学含义，例如 zonal_robert_wind=U=u cosφ、
 * surface_pressure_ratio=exp(π)=ps/p0、mass_weighted_humidity=χ=(ps/p0)q。
 *
 * 本函数先诊断
 *   P=v_h·∇π，A=D+P，Abar=Σ_k Δσ_k A_k，
 * 再由 σdot_{k+1/2}=σ_{k+1/2} Abar-Σ_{j<=k}Δσ_j A_j
 * 构造温度、动量和水汽的垂直项，最后输出供谱分析使用的格点源项/通量。
 * 顶、底边界的 σdot 都为零，所以顶层、内部层和底层必须分别离散。
 *
 * 三维整层场按 [level][point] 展平；界面场只有 level_count-1 层。
 * column_pressure_advection 只输出显式非线性柱积分 Abar_exp=ΣΔσP；线性柱散度
 * ΣΔσDbar 在 spectrala_implicit 中补入（而 ∂π/∂τ=-Abar）。
 * scaled_moist_geopotential 仅为湿虚温引起的位势修正，并非完整位势。
 */
void calcgp(
    int32_t horizontal_count, int32_t level_count,
    float kappa, float water_vapor_heat_capacity_correction,
    float dry_to_vapor_gas_constant_ratio, float reference_surface_pressure,
    float rotation_rate, const float *reciprocal_cosine_squared,
    const float *zonal_robert_wind,
    const float *zonal_log_surface_pressure_gradient,
    const float *meridional_robert_wind,
    const float *meridional_log_surface_pressure_gradient,
    const float *mass_weighted_humidity,
    const float *surface_pressure_ratio,
    const float *temperature_perturbation,
    const float *reference_temperature, const float *divergence,
    const float *sigma_thickness, const float *lower_interface_sigma,
    const float *absolute_vorticity,
    const float *inverse_double_sigma_thickness,
    const float *reference_temperature_difference,
    const float *compression_matrix,
    const float *kappa_reference_temperature,
    const float *hydrostatic_integral,
    float *pressure_vertical_velocity, float *temperature_source,
    float *humidity_vertical_source, float *zonal_momentum_flux,
    float *meridional_momentum_flux, float *column_pressure_advection,
    float *scaled_moist_geopotential)
{
    // L 个整层之间只有 L-1 个内部界面；模式顶 σ=0 和地表 σ=1 不另存数组。
    const int32_t interface_count = level_count - 1;

    // field_values：整层三维场总元素数 = 水平点数 × 垂直层数。
    const size_t field_values =
        (size_t)horizontal_count * (size_t)level_count;

    // interface_values：界面场总元素数 = 水平点数 × (垂直层数 - 1)。
    const size_t interface_values =
        (size_t)horizontal_count * (size_t)interface_count;
    
    // P_k=v_h·∇π：气流沿地表对数气压梯度运动产生的局地变化率。
    float *surface_log_pressure_advection =
        (float *)malloc(field_values * sizeof(float));
    
    // Tv-T0
    float *virtual_temperature_perturbation =
        (float *)malloc(field_values * sizeof(float));

    // Tv-T
    float *virtual_temperature_moisture_increment =
        (float *)malloc(field_values * sizeof(float));

    // Tv/(1+cpv_cpd_minus1*q)-T0：含湿空气热容修正的压缩加热温度扰动；cpv_cpd_minus1 是对应常数。
    float *compression_temperature_perturbation =
        (float *)malloc(field_values * sizeof(float));

    // 自模式顶累计到各内部界面的 ΣΔσ(D+P)，用于诊断 σdot 和 ω。
    float *cumulative_mass_divergence_at_interface =
        (float *)malloc(interface_values * sizeof(float));

    // 定义在相邻整层界面上的 σ 垂直速度。
    float *sigma_dot = (float *)malloc(interface_values * sizeof(float));

    // 由 sigma_dot 与相邻层温度差构造的界面垂直输送项。
    float *temperature_vertical =
        (float *)malloc(interface_values * sizeof(float));

    // 由 sigma_dot 与相邻层湿度和构造的界面项。
    float *humidity_vertical =
        (float *)malloc(interface_values * sizeof(float));

    // 由 sigma_dot 与相邻层 Robert 纬向风 U 差构造的界面项。
    float *zonal_vertical =
        (float *)malloc(interface_values * sizeof(float));

    // 由 sigma_dot 与相邻层 Robert 经向风 V 差构造的界面项。
    float *meridional_vertical =
        (float *)malloc(interface_values * sizeof(float));

    // 整根气柱的 ΣΔσD，在构造所有界面 sigma_dot 时保持这一含义不变。
    float *column_integrated_divergence =
        (float *)malloc((size_t)horizontal_count * sizeof(float));

    // 自模式顶累计到当前整层的 ΣΔσP，用于参考温度的显式垂直输送。
    float *cumulative_pressure_advection =
        (float *)malloc((size_t)horizontal_count * sizeof(float));

    // C(P)：地表气压平流 P 的垂直矩阵积分，逐目标层重新计算。
    float *pressure_advection_vertical_integral =
        (float *)malloc((size_t)horizontal_count * sizeof(float));

    // C(P+D)：总质量散度 A=P+D 的垂直矩阵积分，逐目标层重新计算。
    float *mass_divergence_vertical_integral =
        (float *)malloc((size_t)horizontal_count * sizeof(float));

    int32_t level;
    int32_t point;

    // 任一 malloc 失败就直接 abort
    if (surface_log_pressure_advection == NULL ||
        virtual_temperature_perturbation == NULL ||
        virtual_temperature_moisture_increment == NULL ||
        compression_temperature_perturbation == NULL ||
        cumulative_mass_divergence_at_interface == NULL || sigma_dot == NULL ||
        temperature_vertical == NULL || humidity_vertical == NULL ||
        zonal_vertical == NULL || meridional_vertical == NULL ||
        column_integrated_divergence == NULL ||
        cumulative_pressure_advection == NULL ||
        pressure_advection_vertical_integral == NULL ||
        mass_divergence_vertical_integral == NULL)
    {
        abort();
    }

    // 阶段 1：由 (χ, expπ, T') 恢复 q、T、Tv，并计算 P=v_h·∇π。
    for (level = 0; level < level_count; ++level)
    {
        for (point = 0; point < horizontal_count; ++point)
        {
            const size_t index = field_index(point, level, horizontal_count);

            // physical_humidity = max(χ/(ps/p0), 0)。
            const float physical_humidity =
                max_zero(mass_weighted_humidity[index] /
                         surface_pressure_ratio[point]);

            // 由 T'=T-T0 恢复总温度 T。
            const float total_temperature =
                temperature_perturbation[index] +
                reference_temperature[level];

            // 虚温乘子：1 + (1/rdbrv - 1)q；rdbrv=Rd/Rv。
            const float virtual_factor =
                1.0f +
                (1.0f / dry_to_vapor_gas_constant_ratio - 1.0f) *
                    physical_humidity;

            // Robert 风和 Robert 梯度都含 cosφ 因子；乘 1/cos²φ 后得到 v_h·∇π。
            surface_log_pressure_advection[index] =
                reciprocal_cosine_squared[point] *
                (zonal_robert_wind[index] *
                     zonal_log_surface_pressure_gradient[point] +
                 meridional_robert_wind[index] *
                     meridional_log_surface_pressure_gradient[point]);

            // Tv-T0：相对于参考温度的虚温扰动，供动量通量使用。
            virtual_temperature_perturbation[index] =
                total_temperature * virtual_factor -
                reference_temperature[level];

            // 单独保存 Tv-T。
            virtual_temperature_moisture_increment[index] =
                virtual_temperature_perturbation[index] -
                temperature_perturbation[index];
            
            // Tv/(1+cpv_cpd_minus1*q)-T0：同时考虑虚温与湿空气热容的压缩加热温度扰动。
            compression_temperature_perturbation[index] =
                total_temperature * virtual_factor /
                    (1.0f + water_vapor_heat_capacity_correction *
                                physical_humidity) -
                reference_temperature[level];
        }
    }

    // 阶段 2：自模式顶向下累计 ΣΔσD 和 ΣΔσP；二者之和就是 ΣΔσA。
    for (point = 0; point < horizontal_count; ++point)
    {
        // level=0 时 field_index(point, 0, horizontal_count) 就等于 point，
        // 因此这里直接把 point 转为 size_t 作为顶层下标。
        const size_t top = (size_t)point;

        // 整柱积分从顶层贡献开始累计。
        column_integrated_divergence[point] =
            sigma_thickness[0] * divergence[top];

        // 输出量始终表示整柱 ΣΔσP；此处从顶层贡献开始累计。
        column_pressure_advection[point] =
            sigma_thickness[0] * surface_log_pressure_advection[top];

        // 独立数组保存截至顶层下界面的 ΣΔσ(D+P)。
        cumulative_mass_divergence_at_interface[top] =
            column_integrated_divergence[point] +
            column_pressure_advection[point];

        // 顶层上界为 σ=0，层中心的 ω/(ps Ω) = -(Δσ_0/2) D_0。
        pressure_vertical_velocity[top] =
            -0.5f * sigma_thickness[0] * divergence[top] *
            surface_pressure_ratio[point] * reference_surface_pressure *
            rotation_rate;
    }

    // 继续累计内部整层（不含顶层 0 和底层 level_count-1）。
    for (level = 1; level < level_count - 1; ++level)
    {
        for (point = 0; point < horizontal_count; ++point)
        {
            const size_t index = field_index(point, level, horizontal_count);

            column_integrated_divergence[point] =
                column_integrated_divergence[point] +
                sigma_thickness[level] * divergence[index];
            column_pressure_advection[point] =
                column_pressure_advection[point] +
                sigma_thickness[level] *
                    surface_log_pressure_advection[index];

            cumulative_mass_divergence_at_interface[index] =
                column_integrated_divergence[point] +
                column_pressure_advection[point];
        }
    }

    // 将最底层也并入整柱累计。
    // 最底层没有对应的内部下界面，因此只完成两个整柱积分。
    for (point = 0; point < horizontal_count; ++point)
    {
        const size_t bottom =
            field_index(point, level_count - 1, horizontal_count);
        column_integrated_divergence[point] =
            column_integrated_divergence[point] +
            sigma_thickness[level_count - 1] * divergence[bottom];
        column_pressure_advection[point] =
            column_pressure_advection[point] +
            sigma_thickness[level_count - 1] *
                surface_log_pressure_advection[bottom];
    }

    // 阶段 3：计算内部界面的 σdot 和界面下方整层中心的 ω。
    for (level = 0; level < interface_count; ++level)
    {
        const float half_thickness_below = 0.5f * sigma_thickness[level + 1];
        for (point = 0; point < horizontal_count; ++point)
        {
            // index 指向界面上方整层 level 的平面位置。
            const size_t index = field_index(point, level, horizontal_count);
            // below 指向下一整层 level+1。
            const size_t below = field_index(point, level + 1, horizontal_count);

            // σdot_{k+1/2}=σ_{k+1/2}Σ_jΔσ_jA_j-Σ_{j<=k}Δσ_jA_j。
            sigma_dot[index] =
                lower_interface_sigma[level] *
                    (column_integrated_divergence[point] +
                     column_pressure_advection[point]) -
                cumulative_mass_divergence_at_interface[index];

            // 连续方程给出 ω/(ps Ω) = σ_k P_k - ∫_0^σ_k (D+P) dσ。
            // σ_k P_k 中的半层 P_k 与积分中的对应项相消。
            pressure_vertical_velocity[below] =
                (lower_interface_sigma[level] *
                     surface_log_pressure_advection[below] -
                 cumulative_mass_divergence_at_interface[index] -
                 half_thickness_below * divergence[below]) *
                surface_pressure_ratio[point] * reference_surface_pressure *
                rotation_rate;

            // 温度垂直项使用相邻层温度扰动差值：(下层 - 上层)。
            temperature_vertical[index] =
                sigma_dot[index] *
                (temperature_perturbation[below] -
                 temperature_perturbation[index]);

            // χ 是守恒变量，界面值取相邻层平均；这里先存 σdot(χ_k+χ_{k+1})，
            // 回收到整层时再乘 1/(2Δσ_k)，所以“用和”不是漏写差分。
            humidity_vertical[index] =
                sigma_dot[index] *
                (mass_weighted_humidity[below] +
                 mass_weighted_humidity[index]);

            // U 的垂直项使用相邻层差值：(下层 - 上层)。
            zonal_vertical[index] =
                sigma_dot[index] *
                (zonal_robert_wind[below] - zonal_robert_wind[index]);

            // V 的垂直项同样使用相邻层差值：(下层 - 上层)。
            meridional_vertical[index] =
                sigma_dot[index] *
                (meridional_robert_wind[below] -
                 meridional_robert_wind[index]);
        }
    }

    // 阶段 4：顶层只有下界面；上边界 σdot=0，因此只保留下界面贡献。
    for (point = 0; point < horizontal_count; ++point)
    {
        const size_t index = (size_t)point;

        // 从顶层开始累计 ΣΔσP；该数组在后续各层始终保持这一含义。
        cumulative_pressure_advection[point] =
            surface_log_pressure_advection[index] * sigma_thickness[0];

        // 顶层 temperature_source 由三部分组成：
        //   1) T' D；
        //   2) -κ × compression_temperature_perturbation × D；
        //   3) 温度扰动和参考温度的下界面垂直输送。
        temperature_source[index] =
            temperature_perturbation[index] * divergence[index] -
            kappa * compression_temperature_perturbation[index] * divergence[index] -
            inverse_double_sigma_thickness[0] *
                (temperature_vertical[index] +
                reference_temperature_difference[0] *
                (lower_interface_sigma[0] * column_pressure_advection[point] - cumulative_pressure_advection[point]));

        // 顶层水汽垂直源项只使用其下方界面，并带负号。
        humidity_vertical_source[index] =
            -inverse_double_sigma_thickness[0] *
            humidity_vertical[index];

        // F_U：绝对涡度项、纬向 π 梯度项和 U 的下界面垂直输送。
        zonal_momentum_flux[index] =
            meridional_robert_wind[index] * absolute_vorticity[index] -
            zonal_log_surface_pressure_gradient[point] *
                virtual_temperature_perturbation[index] -
            inverse_double_sigma_thickness[0] * zonal_vertical[index];

        // F_V：绝对涡度项、经向 π 梯度项和 V 的下界面垂直输送。
        meridional_momentum_flux[index] =
            -zonal_robert_wind[index] * absolute_vorticity[index] -
            meridional_log_surface_pressure_gradient[point] *
                virtual_temperature_perturbation[index] -
            inverse_double_sigma_thickness[0] *
                meridional_vertical[index];
    }

    // 阶段 5：内部层同时接收上下两个界面贡献，并用 C 矩阵计算垂直积分。
    for (level = 1; level < level_count - 1; ++level)
    {
        for (point = 0; point < horizontal_count; ++point)
        {
            // 先用 inner=0（顶层）对应的 C 系数初始化两个积分。
            // 下标 level*level_count 即矩阵第 level 行第 0 列。
            const size_t top = (size_t)point;
            pressure_advection_vertical_integral[point] =
                compression_matrix[(size_t)level * (size_t)level_count] *
                surface_log_pressure_advection[top];

            // C(P+D)=C(A) 用于压缩加热的层平均 ω/p。
            mass_divergence_vertical_integral[point] =
                compression_matrix[(size_t)level * (size_t)level_count] *
                (surface_log_pressure_advection[top] + divergence[top]);
        }
        // 单独作用域仅用于限制 inner 的生命周期，不改变控制流。
        {
            int32_t inner;
            // 将 inner=1 ... level 的贡献逐层累加到两个 C 积分。
            // 因而对当前目标层 level，积分只使用 C 的第 0..level 列。
            for (inner = 1; inner <= level; ++inner)
            {
                // 取 C 矩阵的 [level][inner] 系数。
                const float coefficient =
                    compression_matrix[
                        (size_t)level * (size_t)level_count +
                        (size_t)inner];
                for (point = 0; point < horizontal_count; ++point)
                {
                    const size_t index =
                        field_index(point, inner, horizontal_count);
                    pressure_advection_vertical_integral[point] =
                        pressure_advection_vertical_integral[point] +
                        coefficient * surface_log_pressure_advection[index];
                    mass_divergence_vertical_integral[point] =
                        mass_divergence_vertical_integral[point] +
                        coefficient *
                            (surface_log_pressure_advection[index] + divergence[index]);
                }
            }
        }
        // 使用刚得到的 C(P) 与 C(P+D) 计算当前内部层的输出。
        for (point = 0; point < horizontal_count; ++point)
        {
            const size_t index = field_index(point, level, horizontal_count);
            // upper_interface 指向当前整层与其上方整层之间的界面（level-1）。
            // 当前层下方界面在线性存储中直接使用 index（对应 interface level）。
            const size_t upper_interface =
                field_index(point, level - 1, horizontal_count);
            // 把当前层 ΔσP 加入自顶向下累计量，供参考温度垂直输送使用。
            cumulative_pressure_advection[point] =
                cumulative_pressure_advection[point] +
                surface_log_pressure_advection[index] * sigma_thickness[level];
            
            // 内部层 temperature_source 同时包含：
            //   - T'D；
            //   - 湿压缩加热温度扰动 * (P-C(A))；
            //   - κT0 * (P-C(P))；
            //   - 上下两个界面的 temperature_vertical 与 sigma 相关修正。
            temperature_source[index] =
                temperature_perturbation[index] * divergence[index] +
                kappa * compression_temperature_perturbation[index] *
                    (surface_log_pressure_advection[index] -
                     mass_divergence_vertical_integral[point]) +
                kappa_reference_temperature[level] *
                    (surface_log_pressure_advection[index] -
                     pressure_advection_vertical_integral[point]) -
                inverse_double_sigma_thickness[level] *
                    (temperature_vertical[index] +
                     temperature_vertical[upper_interface] +
                     column_pressure_advection[point] *
                         (reference_temperature_difference[level] *
                              lower_interface_sigma[level] +
                          reference_temperature_difference[level - 1] *
                              lower_interface_sigma[level - 1]) -
                     cumulative_pressure_advection[point] *
                         (reference_temperature_difference[level - 1] +
                          reference_temperature_difference[level]) +
                     surface_log_pressure_advection[index] *
                         sigma_thickness[level] *
                         reference_temperature_difference[level - 1]);

            // 内部层水汽垂直源项使用“下方界面 - 上方界面”的通量差。
            humidity_vertical_source[index] =
                -inverse_double_sigma_thickness[level] *
                (humidity_vertical[index] -
                 humidity_vertical[upper_interface]);
            
            // F_U 同时接收当前层水平项和 U 在上下两个界面的垂直输送。
            zonal_momentum_flux[index] =
                meridional_robert_wind[index] * absolute_vorticity[index] -
                zonal_log_surface_pressure_gradient[point] *
                    virtual_temperature_perturbation[index] -
                inverse_double_sigma_thickness[level] *
                    (zonal_vertical[index] + zonal_vertical[upper_interface]);

            // F_V 对应使用 V 在上下两个界面的垂直输送。
            meridional_momentum_flux[index] =
                -zonal_robert_wind[index] * absolute_vorticity[index] -
                meridional_log_surface_pressure_gradient[point] *
                    virtual_temperature_perturbation[index] -
                inverse_double_sigma_thickness[level] *
                    (meridional_vertical[index] +
                     meridional_vertical[upper_interface]);
        }
    }

    // 阶段 6：底层下边界 σdot=0，只保留唯一的上界面贡献。
    level = level_count - 1;
    // 与内部层相同，先计算当前底层对应的 C 垂直积分。
    for (point = 0; point < horizontal_count; ++point)
    {
        const size_t top = (size_t)point;
        pressure_advection_vertical_integral[point] =
            compression_matrix[(size_t)level * (size_t)level_count] *
            surface_log_pressure_advection[top];
        mass_divergence_vertical_integral[point] =
            compression_matrix[(size_t)level * (size_t)level_count] *
            (surface_log_pressure_advection[top] + divergence[top]);
    }

    // 单独作用域用于局部声明 inner。
    {
        int32_t inner;
        // 底层积分覆盖 inner=1 ... level_count-1，
        // 与内部层的 inner<=level 一致，只是此时 level 已是最后一层。
        for (inner = 1; inner < level_count; ++inner)
        {
            const float coefficient =
                compression_matrix[
                    (size_t)level * (size_t)level_count + (size_t)inner];
            for (point = 0; point < horizontal_count; ++point)
            {
                const size_t index =
                    field_index(point, inner, horizontal_count);
                pressure_advection_vertical_integral[point] =
                    pressure_advection_vertical_integral[point] +
                    coefficient * surface_log_pressure_advection[index];
                mass_divergence_vertical_integral[point] =
                    mass_divergence_vertical_integral[point] +
                    coefficient *
                        (surface_log_pressure_advection[index] + divergence[index]);
            }
        }
    }
    // 使用底层积分结果和唯一的上方界面，计算底层各输出量。
    for (point = 0; point < horizontal_count; ++point)
    {
        const size_t index = field_index(point, level, horizontal_count);
        // interface = level-1，即底层与倒数第二层之间的唯一相邻界面。
        const size_t interface =
            field_index(point, level - 1, horizontal_count);
        
        // 底层温度源项不含下方界面，只保留唯一的上方界面贡献。
        temperature_source[index] =
            temperature_perturbation[index] * divergence[index] +
            kappa * compression_temperature_perturbation[index] *
                (surface_log_pressure_advection[index] -
                 mass_divergence_vertical_integral[point]) +
            kappa_reference_temperature[level] *
                (surface_log_pressure_advection[index] -
                 pressure_advection_vertical_integral[point]) -
            inverse_double_sigma_thickness[level] *
                (temperature_vertical[interface] +
                 reference_temperature_difference[level - 1] *
                     (lower_interface_sigma[level - 1] * column_pressure_advection[point] - cumulative_pressure_advection[point]));

        // 底层水汽垂直源项的符号与顶层对应边界项相反。
        humidity_vertical_source[index] =
            inverse_double_sigma_thickness[level] *
            humidity_vertical[interface];

        // 底层 F_U 只包含唯一上方界面的 U 垂直输送修正。
        zonal_momentum_flux[index] =
            meridional_robert_wind[index] * absolute_vorticity[index] -
            zonal_log_surface_pressure_gradient[point] *
                virtual_temperature_perturbation[index] -
            inverse_double_sigma_thickness[level] *
                zonal_vertical[interface];

        // 底层 F_V 只包含唯一上方界面的 V 垂直输送修正。
        meridional_momentum_flux[index] =
            -zonal_robert_wind[index] * absolute_vorticity[index] -
            meridional_log_surface_pressure_gradient[point] *
                virtual_temperature_perturbation[index] -
            inverse_double_sigma_thickness[level] *
                meridional_vertical[interface];
    }

    // 阶段 7：静力矩阵 G 沿气柱积分 δTv=Tv-T，得到湿非线性位势修正。
    // 乘 2/rcsq=2cos²φ 后，可与 U²+V²=2cos²φ K 共用后续谱散度算子。
    for (level = 0; level < level_count; ++level)
    {
        for (point = 0; point < horizontal_count; ++point)
        {
            // hydrostatic_sum 保存 G 当前行与该水平点整根垂直列的点积。
            float hydrostatic_sum = 0.0f;
            int32_t inner;
            const size_t index = field_index(point, level, horizontal_count);
            // 遍历全部垂直层 inner=0 ... level_count-1。
            for (inner = 0; inner < level_count; ++inner)
            {
                // 累加 G[level][inner] * (Tv-T)[inner][point]。
                hydrostatic_sum =
                    hydrostatic_sum +
                    hydrostatic_integral[(size_t)level *
                                             (size_t)level_count +
                                         (size_t)inner] *
                        virtual_temperature_moisture_increment[
                            field_index(point, inner, horizontal_count)];
            }
            scaled_moist_geopotential[index] =
                hydrostatic_sum * 2.0f /
                reciprocal_cosine_squared[point];
        }
    }

    // 阶段 8：输出数组归调用者所有；此处只释放柱诊断工作区。
    free(mass_divergence_vertical_integral);
    free(pressure_advection_vertical_integral);
    free(cumulative_pressure_advection);
    free(column_integrated_divergence);
    free(meridional_vertical);
    free(zonal_vertical);
    free(humidity_vertical);
    free(temperature_vertical);
    free(sigma_dot);
    free(cumulative_mass_divergence_at_interface);
    free(compression_temperature_perturbation);
    free(virtual_temperature_moisture_increment);
    free(virtual_temperature_perturbation);
    free(surface_log_pressure_advection);
}
