#ifndef PLASIC_ICEMOD_KERNELS_H
#define PLASIC_ICEMOD_KERNELS_H

#include <stdint.h>

// 参数名与教程 4.3 示意图符号的对应关系见 icemod_kernels.c 文件头。

// 冰能量账：Q_net += Q_o->i；Q_ice = Q_net + Q_s->i - Q_f；更新 h_i；
// 输出 Q_i->o 与超额融化 Q_ex。
void ice_make_ice(
    int32_t nhor, float freezing_temperature, float ice_density,
    float latent_heat_ice, float timestep, const float *land_mask,
    const float *ocean_to_ice_flux, const float *sst,
    const float *snow_to_ice_flux, const float *freezing_flux,
    float *net_heat_flux, float *ice_energy_flux, float *excess_melt_flux,
    float *ice_to_ocean_flux, float *ice_thickness);

// 由 h_i 的变化更新连续紧密度 c_i^*。
void ice_compactness(
    int32_t nhor, const float *old_thickness, const float *new_thickness,
    float *compactness);

// 积雪：开阔水降雪、先融雪后融冰（Q_m -> M_s, Q_m^rem）、雪转冰 Q_s->i；
// 把相应项累加到 net_heat_flux（Q_net）。
void ice_snow(
    int32_t nhor, float minimum_ice_cover, float seawater_density,
    float ice_density, float latent_heat_snow, float timestep,
    const float *land_mask, const float *ice_cover, float *snowfall,
    const float *melt_energy, const float *ice_thickness, float *snow,
    float *net_heat_flux, float *snow_melt, float *snow_melt_flux,
    float *snow_to_ice_flux, float *snow_to_ice_rate);

// 冰雪柱导热 Q_c = K_z (T_i - T_o) 与冻结点导热 Q_f = K_z (T_f - T_o)；
// h_i < h_min 时 Q_f = 0，开阔水 Q_c = Q_a^mod，薄冰 Q_c = 0（热量在 Q_m 中）。
void ice_conductive_flux(
    int32_t nhor, int32_t ocean_flux_correction, int32_t prognostic_ice,
    float minimum_thickness, float maximum_thickness,
    float freezing_temperature, float snow_density, float ice_conductivity,
    float snow_conductivity, const float *ice_thickness, const float *snow,
    const float *surface_temperature, const float *sst,
    const float *modified_atmospheric_heat,
    float *diagnosed_conductive_flux, float *conductive_flux,
    float *freezing_flux);

// 半隐式表面温度 T_i 与融化能量 Q_m。
void ice_skin_temperature(
    int32_t nhor, float minimum_thickness, float maximum_thickness,
    float melting_temperature, float ice_density, float ice_heat_capacity,
    float snow_density, float ice_conductivity, float snow_conductivity,
    float timestep, const float *ice_thickness, const float *snow,
    const float *sst, const float *modified_atmospheric_heat,
    float *surface_temperature, float *surface_storage_flux,
    float *melt_energy);

void ice_interpolate_climatology(
    int32_t nhor, int32_t previous_month1, int32_t previous_month2,
    float previous_weight1, float previous_weight2, int32_t month1,
    int32_t month2, float weight1, float weight2,
    const float *climatological_sst, const float *climatological_cover,
    const float *climatological_thickness, float *previous_sst,
    float *current_sst, float *current_cover, float *current_thickness);

void ice_make_thickness(
    int32_t nlon, int32_t nlat, int32_t nmonths,
    const float *compactness, float *thickness);

#endif
