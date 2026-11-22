#ifndef PLASIC_SLAB_OCEAN_H
#define PLASIC_SLAB_OCEAN_H

#include <stdint.h>

/*
 * Numerical kernels for oceanmod.f90. Arrays use the model's column-major
 * storage. Month indices passed to the interpolation kernel are the original
 * 0..13 indices used by OCEANMOD.
 *
 * The mixed-layer ocean is forced to a single layer and has no vertical or
 * horizontal diffusion: each ocean point only exchanges heat with the
 * atmosphere, the sea ice and the deep ocean.
 *
 * Parameter names follow the symbols of the tutorial figures (4.3, 4.4):
 *   sst                    T_o           mixed_layer_depth   H
 *   trial_sst              T_o^*         mixed_layer_timestep dt_mix
 *   seawater_density       rho_o         seawater_heat_capacity c_p,o
 *   ice_to_ocean_flux      Q_i->o        deep_ocean_flux     Q_deep
 *   ice_flux               -Q_o->i       freezing_temperature T_f
 * so that C_o = rho_o c_p,o H.
 */

void ocean_initialize_parameters(
    int32_t nhor, int32_t steps_per_day, float solar_day,
    float *restoring_timescale, float layer_depth,
    float *mixed_layer_depth, float *mixed_layer_timestep);

void ocean_interpolate_cycle(
    int32_t nhor, int32_t month1, int32_t month2, float weight,
    const float *monthly, float *field);

void ocean_copy_outputs(
    int32_t nhor, const float *sst, const float *mixed_layer_depth,
    const float *ice_flux, float *ice_sst, float *ice_mixed_layer_depth,
    float *ocean_to_ice_flux);

void ocean_copy_inputs(
    int32_t nhor, float ice_density, const float *input_ice_cover,
    const float *input_ice_thickness, const float *input_ice_to_ocean_flux,
    const float *input_freshwater_flux, const float *input_stress_x,
    const float *input_stress_y, const float *input_friction_velocity_cubed,
    const float *input_snow, const float *input_climatological_ice_thickness,
    float *ice_cover, float *ice_thickness, float *ice_to_ocean_flux,
    float *precipitation_minus_evaporation, float *stress_x,
    float *stress_y, float *friction_velocity_cubed, float *ice_and_snow,
    float *climatological_ice_thickness);

void ocean_mksst_begin(
    int32_t nhor, float mixed_layer_timestep, float seawater_density,
    float seawater_heat_capacity, float freezing_temperature,
    const float *land_mask, const float *ice_thickness,
    const float *ice_to_ocean_flux, const float *deep_ocean_flux,
    const float *mixed_layer_depth, const float *sst, double *trial_sst,
    float *ice_flux);

void ocean_mksst_finish(
    int32_t nhor, int32_t prognostic_ocean, int32_t flux_correction,
    float mixed_layer_timestep, float seawater_density,
    float seawater_heat_capacity, float freezing_temperature,
    const float *land_mask, const float *ice_thickness,
    const float *mixed_layer_depth, double *trial_sst, float *ice_flux,
    float *sst);

int32_t ocean_finish_step(
    int32_t nhor, int32_t prognostic_ocean, int32_t flux_correction,
    int32_t nstep, int32_t output_interval, float mixed_layer_timestep,
    float restoring_timescale, float seawater_density, float ice_density,
    float seawater_heat_capacity, float latent_heat_ice,
    float freezing_temperature, const float *land_mask,
    const float *climatological_sst, const float *mixed_layer_depth,
    const float *ice_thickness, const float *climatological_ice_thickness,
    const float *ice_to_ocean_flux, const float *deep_ocean_flux,
    float *sst, float *ice_flux, float *residual_ice_flux,
    float *ice_to_ocean_flux_accumulator, float *ice_flux_accumulator,
    float *deep_ocean_flux_accumulator, int32_t *accumulation_count);

void ocean_reset_diagnostics(
    int32_t nhor, float *ice_to_ocean_flux_accumulator,
    float *ice_flux_accumulator, float *deep_ocean_flux_accumulator,
    int32_t *accumulation_count);

#endif
