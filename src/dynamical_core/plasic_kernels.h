#ifndef PLASIC_MAIN_KERNELS_H
#define PLASIC_MAIN_KERNELS_H

#include <stdint.h>

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
    float *scaled_moist_geopotential);

void exp_array(
    int32_t count, const float *log_values, float *exponential_values);

void gridpointa_products(
    int32_t count, const float *zonal_wind, const float *meridional_wind,
    const float *temperature, const float *humidity,
    const float *scaled_moist_geopotential, float *zonal_temperature_flux,
    float *meridional_temperature_flux,
    float *scaled_kinetic_plus_geopotential,
    float *zonal_humidity_flux, float *meridional_humidity_flux);

void gridpointa_snapshots(
    int32_t horizontal_count, int32_t level_count, float reference_pressure,
    float velocity_scale, const float *rcsq,
    const float *surface_pressure_ratio,
    const float *zonal_wind, const float *meridional_wind,
    float *saved_surface_pressure, float *saved_zonal_wind,
    float *saved_meridional_wind);

void gridpointd_zero(
    int32_t horizontal_count, int32_t level_count,
    float *zonal_tendency, float *meridional_tendency,
    float *temperature_tendency, float *humidity_tendency,
    float *physical_zonal_tendency, float *physical_meridional_tendency,
    float *physical_temperature_tendency, float *physical_humidity_tendency);

void gridpointd_physical(
    int32_t horizontal_count, int32_t level_count,
    float reference_pressure, float velocity_scale, float temperature_scale,
    const float *rcsq, const float *log_surface_pressure,
    const float *zonal_wind, const float *meridional_wind,
    const float *temperature, const float *reference_temperature,
    const float *humidity, float *surface_pressure,
    float *physical_zonal_wind, float *physical_meridional_wind,
    float *physical_temperature, float *physical_humidity);

void gridpointd_tendencies(
    int32_t horizontal_count, int32_t level_count,
    float reference_pressure, float velocity_scale, float temperature_scale,
    float rotation_rate, const float *rcsq, const float *surface_pressure,
    const float *physical_zonal_tendency,
    const float *physical_meridional_tendency,
    const float *physical_temperature_tendency,
    const float *physical_humidity_tendency,
    float *zonal_tendency, float *meridional_tendency,
    float *temperature_tendency, float *humidity_tendency);

void spectrala_implicit(
    int32_t spectral_count, int32_t level_count, int32_t advection_enabled,
    int32_t process_id, int32_t root_process, float time_step,
    const int32_t *wavenumber, const float *vertical_g,
    const float *reference_temperature, const float *inverse_matrix,
    const float *surface_geopotential, const float *sigma_thickness,
    const float *vertical_tau, const float *column_pressure_advection,
    const float *explicit_divergence_tendency,
    const float *explicit_temperature_tendency,
    float *column_mass_divergence, float *centered_divergence,
    float *complete_temperature_tendency,
    const float *old_surface_log_pressure, const float *old_divergence,
    const float *old_vorticity, const float *old_temperature,
    const float *old_humidity, float *saved_surface_log_pressure,
    float *saved_divergence, float *saved_vorticity,
    float *saved_temperature, float *saved_humidity);

void spectrala_finalize(
    int32_t spectral_count, int32_t level_count, int32_t advection_enabled,
    int32_t process_id, int32_t root_process, int32_t short_step_count,
    float double_time_step, float filter_weight, float filter_center_weight,
    const float *column_mass_divergence,
    const float *centered_divergence,
    const float *vorticity_tendency,
    const float *complete_temperature_tendency,
    const float *humidity_tendency,
    const float *saved_surface_log_pressure,
    const float *saved_divergence, const float *saved_vorticity,
    const float *saved_temperature, const float *saved_humidity,
    float *filtered_surface_log_pressure, float *filtered_divergence,
    float *filtered_vorticity, float *filtered_temperature,
    float *filtered_humidity, float *surface_log_pressure,
    float *divergence, float *vorticity, float *temperature, float *humidity);

void spectrald_add_physics(
    int32_t count, float double_time_step,
    const float *vorticity_tendency, const float *temperature_tendency,
    const float *divergence_tendency, const float *humidity_tendency,
    float *vorticity, float *temperature, float *divergence, float *humidity);

void spectrald_dissipation(
    int32_t spectral_count, int32_t level_count, int32_t process_id,
    int32_t root_process, int32_t sponge_enabled, float double_time_step,
    float planetary_vorticity, float sponge_rate,
    const float *temperature_drag_rate,
    const float *temperature_diffusion_rate,
    const float *humidity_diffusion_rate, const float *friction_rate,
    const float *divergence_diffusion_rate,
    const float *vorticity_diffusion_rate,
    const float *spectral_diffusion_scale,
    const float *restoration_temperature, const float *divergence,
    const float *vorticity, const float *temperature, const float *humidity,
    float *surface_log_pressure, float *temperature_tendency,
    float *humidity_tendency, float *divergence_tendency,
    float *vorticity_tendency, float *divergence_friction,
    float *divergence_diffusion_tendency, float *vorticity_friction,
    float *vorticity_diffusion_tendency);

void spectrald_finalize(
    int32_t spectral_count, int32_t level_count, int32_t short_step_count,
    float double_time_step, float filter_weight,
    const float *divergence_tendency, const float *vorticity_tendency,
    const float *temperature_tendency, const float *humidity_tendency,
    float *filtered_surface_log_pressure, float *filtered_divergence,
    float *filtered_vorticity, float *filtered_temperature,
    float *filtered_humidity, float *surface_log_pressure,
    float *divergence, float *vorticity, float *temperature, float *humidity);

void initpm_vertical(
    int32_t level_count, int32_t sigma_mode, const float *configured_sigma,
    float *sigma_half, float *sigma_thickness,
    float *reciprocal_double_sigma, float *sigma_full);

void initpm_spectral(
    int32_t level_count, int32_t truncation, int32_t diffusion_offset,
    int32_t spectral_count, float two_pi,
    float temperature_scale, const int32_t *diffusion_order,
    const float *restoration_time_days, const float *friction_time_days,
    const float *divergence_diffusion_time_days,
    const float *vorticity_diffusion_time_days,
    const float *temperature_diffusion_time_days,
    const float *humidity_diffusion_time_days,
    const float *reference_temperature_kelvin, float *temperature_drag_rate,
    float *friction_rate, float *divergence_diffusion_rate,
    float *vorticity_diffusion_rate, float *temperature_diffusion_rate,
    float *humidity_diffusion_rate, float *spectral_diffusion_scale,
    int32_t *wavenumber, float *spectral_norm,
    float *nondimensional_reference_temperature);

void makebm(
    int32_t level_count, int32_t truncation, float time_step,
    const float *reference_temperature, const float *sigma_thickness,
    const float *vertical_g, const float *vertical_tau,
    float *inverse_matrix);

void initsi(
    int32_t level_count, float kappa,
    const float *reference_temperature, const float *sigma_half,
    const float *sigma_thickness, const float *reciprocal_double_sigma,
    float *kappa_reference_temperature,
    float *reference_temperature_difference, float *vertical_g,
    float *vertical_c, float *vertical_tau);

void inilat(
    int32_t latitude_count, const double *sine_latitude,
    float *cosine_squared, float *reciprocal_cosine);

#endif
