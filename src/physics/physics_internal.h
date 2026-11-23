#ifndef PLASIC_PHYSICS_INTERNAL_H
#define PLASIC_PHYSICS_INTERNAL_H

#include "plasic_physics.h"

#include "fluxmod_kernels.h"
#include "land/landmod_kernels.h"
#include "math_tools/sht_kernels.h"
#include "mpi/plasic_mpi.h"
#include "plasic_status.h"
#include "state/plasicmod_state.h"
#include "radmod_kernels.h"
#include "rainmod_kernels.h"
#include "sea/icemod_kernels.h"
#include "sea/slab_ocean.h"
#include "sea/3docean.h"
#include "sea/seamod_kernels.h"
#include "tools/calmod_kernels.h"
#include "tools/miscmod_kernels.h"
#include "tools/plasic_stream.h"
#include "tools/restartmod_io.h"
#include "tools/surface_data_io.h"
#include "vdiff_kernel.h"

#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Private physics state. Prognostic and accumulated fields remain allocated
 * for the complete run and keep their original flat array layout.
 */
typedef struct misc_state
{
    int32_t moisture_fixer;
    float nudging_time;
    float moisture_weight[PLASIC_NHOR];
} misc_state;

typedef struct flux_state
{
    int32_t vertical_diffusion;
    int32_t sensible_heat;
    int32_t evaporation;
    int32_t stress;
    int32_t surface_temperature_mode;
    float minimum_wind;
    float mixing_length;
    float stability_b;
    float stability_c;
    float stability_d;
    float heat_transfer[PLASIC_NHOR];
    float momentum_transfer[PLASIC_NHOR];
} flux_state;

typedef struct rain_state
{
    int32_t beta_mode;
    int32_t large_scale_rain;
    int32_t convective_rain;
    int32_t clouds;
    int32_t dry_adjustment;
    int32_t surface_convection;
    int32_t mix_momentum;
    int32_t shallow_convection;
    int32_t evaporate_precipitation;
    int32_t beta_exponent;
    float critical_relative_humidity;
    float prescribed_beta;
    float cloud_vertical_velocity_1;
    float cloud_vertical_velocity_2;
    float deep_pressure;
    float shallow_top_pressure;
    float shallow_diffusion;
    float evaporation_rate;
    float cloud_factor;
    float critical_humidity[PLASIC_NLEV];
    int32_t convective_layer[PLASIC_NHOR * PLASIC_NLEV];
    int32_t convective_layer_count[PLASIC_NHOR];
    float convective_rain_layer[PLASIC_NHOR * PLASIC_NLEV];
    float large_scale_rain_layer[PLASIC_NHOR * PLASIC_NLEV];
    float convective_snow_layer[PLASIC_NHOR * PLASIC_NLEV];
    float large_scale_snow_layer[PLASIC_NHOR * PLASIC_NLEV];
} rain_state;

typedef struct radiation_state
{
    int32_t use_synthetic_ozone;
    int32_t external_ozone_column_reported;
    int32_t shortwave;
    int32_t longwave;
    int32_t rayleigh_scattering;
    float solar_constant;
    float water_continuum;
    float cloud_tuning_1;
    float cloud_tuning_2;
    float cloud_tuning_3;
    float cloud_absorption;
    float ozone_a0;
    float ozone_a1;
    float ozone_seasonal_amplitude;
    float ozone_height;
    float ozone_thickness;
    float ozone_seasonal_offset;
    float orbital_eccentricity;
    float obliquity_radians;
    float mean_longitude_equinox;
    float perihelion_plus_pi;
    float sun_distance_factor;
    float cosine_zenith[PLASIC_NHOR];
    float ozone[PLASIC_NHOR * PLASIC_NLEV];
    float carbon_dioxide[PLASIC_NHOR * PLASIC_NLEV];
    float longwave_tendency[PLASIC_NHOR * PLASIC_NLEV];
    float shortwave_tendency[PLASIC_NHOR * PLASIC_NLEV];
} radiation_state;

enum
{
    PLASIC_LAND_SOIL_LEVELS = 5,
    PLASIC_CLIMATOLOGY_MONTHS = 14,
    PLASIC_LAND_CLIMATOLOGY_MONTHS = PLASIC_CLIMATOLOGY_MONTHS
};

typedef struct land_state
{
    int32_t prognostic_temperature;
    int32_t prognostic_water;
    float land_albedo;
    float snow_albedo_minimum;
    float snow_albedo_maximum;
    float glacier_albedo_minimum;
    float glacier_albedo_maximum;
    float land_roughness;
    float land_wetness;
    float full_wetness_fraction;
    float glacier_height;
    float top_layer_depth;
    float maximum_snow_depth;
    float maximum_soil_water;
    float snow_density;
    float soil_conductivity;
    float ice_conductivity;
    float snow_conductivity;
    float soil_heat_capacity;
    float ice_heat_capacity;
    float snow_heat_capacity;
    float soil_depth[PLASIC_LAND_SOIL_LEVELS];
    float surface_temperature[PLASIC_NHOR];
    float previous_surface_temperature[PLASIC_NHOR];
    float surface_humidity[PLASIC_NHOR];
    float soil_temperature[PLASIC_NHOR * PLASIC_LAND_SOIL_LEVELS];
    float snow_temperature[PLASIC_NHOR];
    float snow_depth[PLASIC_NHOR];
    float soil_water_flux[PLASIC_NHOR];
    float temperature_climatology[
        PLASIC_NHOR * PLASIC_LAND_CLIMATOLOGY_MONTHS];
    float wetness_climatology[
        PLASIC_NHOR * PLASIC_LAND_CLIMATOLOGY_MONTHS];
    float albedo_climatology[
        PLASIC_NHOR * PLASIC_LAND_CLIMATOLOGY_MONTHS];
    float interpolated_temperature[PLASIC_NHOR];
    float interpolated_wetness[PLASIC_NHOR];
    float roughness_climatology[PLASIC_NHOR];
    float interpolated_albedo[PLASIC_NHOR];
} land_state;

/*
 * The original Fortran model kept the sea coupler, thermodynamic sea ice,
 * and mixed-layer ocean in three modules.  Keep the same ownership here:
 * these are prognostic or accumulated states and therefore must survive a
 * restart; they are not temporary work arrays.
 */
typedef struct sea_state
{
    int32_t coupling_interval;
    int32_t accumulation_count;
    float sea_albedo;
    float ice_albedo;
    float sea_roughness;
    float ice_roughness;
    float sea_relative_humidity;
    float ice_relative_humidity;
    float charnock;
    float surface_temperature[PLASIC_NHOR];
    float surface_humidity[PLASIC_NHOR];
    float sst[PLASIC_NHOR];
    float coupled_sst[PLASIC_NHOR];
    float coupled_mixed_layer_depth[PLASIC_NHOR];
    float coupled_surface_temperature[PLASIC_NHOR];
    float coupled_ice_cover[PLASIC_NHOR];
    float coupled_ice_thickness[PLASIC_NHOR];
    float coupled_snow[PLASIC_NHOR];
    float coupled_snow_melt[PLASIC_NHOR];
    float coupled_snow_change[PLASIC_NHOR];
    float heat_accumulator[PLASIC_NHOR];
    float freshwater_accumulator[PLASIC_NHOR];
    float snowfall_accumulator[PLASIC_NHOR];
    float stress_x_accumulator[PLASIC_NHOR];
    float stress_y_accumulator[PLASIC_NHOR];
    float friction_accumulator[PLASIC_NHOR];
    float sensible_heat_accumulator[PLASIC_NHOR];
    float sensible_derivative_accumulator[PLASIC_NHOR];
    float latent_heat_accumulator[PLASIC_NHOR];
    float latent_derivative_accumulator[PLASIC_NHOR];
    float shortwave_accumulator[PLASIC_NHOR];
    float longwave_accumulator[PLASIC_NHOR];
} sea_state;

typedef struct ocean_state
{
    int32_t accumulation_count;
    int32_t output_interval;
    int32_t prognostic_ocean;
    int32_t flux_correction;
    /* Newtonian SST restoring time (flux_correction == 1): days as set,
     * seconds after ocean_initialize_parameters(). */
    float restoring_timescale;
    float timestep;
    float layer_depth;
    /*
     * Slab mixed layer.  Symbols are those of the tutorial figures in
     * 4.3 / 4.4: T_o = temperature, H = mixed_layer_depth,
     * Q_i->o = ice_to_ocean_flux, Q_deep = deep_ocean_flux and
     * Q_o->i = -ice_flux.
     */
    float land_mask[PLASIC_NHOR];
    float temperature[PLASIC_NHOR];
    float mixed_layer_depth[PLASIC_NHOR];
    float ice_cover[PLASIC_NHOR];
    float ice_thickness[PLASIC_NHOR];
    float climatological_ice_thickness[PLASIC_NHOR];
    /* Q_i->o: sea-ice flux handed to the slab, averaged over the window */
    float ice_to_ocean_flux[PLASIC_NHOR];
    /* Q_deep: deep-ocean flux into the slab (0 in slab-only builds) */
    float deep_ocean_flux[PLASIC_NHOR];
    float precipitation_minus_evaporation[PLASIC_NHOR];
    float stress_x[PLASIC_NHOR];
    float stress_y[PLASIC_NHOR];
    float friction_velocity_cubed[PLASIC_NHOR];
    /* -Q_o->i: heat the slab cannot take under ice or at T_f */
    float ice_flux[PLASIC_NHOR];
    float residual_ice_flux[PLASIC_NHOR];
    float climatological_sst[
        PLASIC_NHOR * PLASIC_CLIMATOLOGY_MONTHS];
    float interpolated_sst[PLASIC_NHOR];
    float ice_and_snow[PLASIC_NHOR];
    float ice_to_ocean_flux_accumulator[PLASIC_NHOR];
    float ice_flux_accumulator[PLASIC_NHOR];
    float deep_ocean_flux_accumulator[PLASIC_NHOR];
} ocean_state;

/*
 * Halo-extended ocean storage: one ghost latitude row on each side of the
 * process's interior rows, so meridional finite differences stay two-sided
 * in parallel runs.
 */
#define PLASIC_DEEP_EXT_NLAT (PLASIC_NLPP + 2)
#define PLASIC_DEEP_EXT_HOR (PLASIC_NLON * PLASIC_DEEP_EXT_NLAT)

/*
 * Three-dimensional deep-ocean state.  This sits below the single-layer
 * mixed-layer slab in `ocean_state`.  The slab remains the atmosphere-facing
 * surface layer; the deep ocean stores a full three-dimensional temperature,
 * salinity and velocity field on the same Gaussian grid, halo-extended in
 * latitude.
 */
typedef struct deep_ocean_state
{
    int32_t enabled;
    int32_t initialized;
    int32_t restart_loaded;
    int32_t nlev;
    float timestep;
    float radius;
    float solar_day;
    float sidereal_day;
    float gravity;
    float cp;
    float rho0;
    float alpha;
    float beta;
    float kappa_v;
    float kappa_h;
    float t_freeze;
    float layer_thickness[PLASIC_DEEP_NLEV];
    float layer_center_depth[PLASIC_DEEP_NLEV];
    float sine_latitude[PLASIC_DEEP_EXT_NLAT];
    float cosine_latitude[PLASIC_DEEP_EXT_NLAT];
    float land_mask[PLASIC_DEEP_EXT_HOR];
    float temperature[PLASIC_DEEP_EXT_HOR * PLASIC_DEEP_NLEV];
    float salinity[PLASIC_DEEP_EXT_HOR * PLASIC_DEEP_NLEV];
    float u_velocity[PLASIC_DEEP_EXT_HOR * PLASIC_DEEP_NLEV];
    float v_velocity[PLASIC_DEEP_EXT_HOR * PLASIC_DEEP_NLEV];
    float w_velocity[PLASIC_DEEP_EXT_HOR * PLASIC_DEEP_NLEV];
    float density[PLASIC_DEEP_EXT_HOR * PLASIC_DEEP_NLEV];
} deep_ocean_state;

typedef struct ice_state
{
    int32_t prognostic_ice;
    int32_t prognostic_snow;
    int32_t prognostic_skin_temperature;
    int32_t ocean_coupling_interval;
    int32_t output_interval;
    int32_t maximum_thickness_correction;
    int32_t output_accumulation_count;
    int32_t ocean_accumulation_count;
    int32_t diagnosed_climatological_thickness;
    float timestep;
    float minimum_thickness;
    float maximum_thickness;
    float compactness_threshold;
    float minimum_ice_cover;
    /*
     * Field names follow the symbols of the tutorial figures in 4.3
     * (sea_ice_algorithm.svg, sea_ice_flux_structure.svg); the full
     * symbol table is at the top of src/sea/icemod_kernels.c.  Fluxes are
     * W m^-2 and positive downward / towards melting.
     */
    float land_mask[PLASIC_NHOR];
    /* T_i */
    float surface_temperature[PLASIC_NHOR];
    /* T_o, returned by the slab */
    float sst[PLASIC_NHOR];
    float mixed_layer_depth[PLASIC_NHOR];
    /* h_i, c_i (0/1 mask of compactness) and S (m water equivalent) */
    float ice_thickness[PLASIC_NHOR];
    float ice_cover[PLASIC_NHOR];
    float snow[PLASIC_NHOR];
    /* M_s and rho_w L_s M_s (diagnostics) */
    float snow_melt[PLASIC_NHOR];
    float snow_melt_flux[PLASIC_NHOR];
    /* Q_ice = Q_net + Q_s->i - Q_f (+ Q_max), the ice energy ledger */
    float ice_energy_flux[PLASIC_NHOR];
    float snow_change[PLASIC_NHOR];
    /* Q_m (Q_a^mod on thin ice, 0 < h_i < h_min) */
    float melt_energy[PLASIC_NHOR];
    /* snow converted to ice by flooding (m water equivalent s^-1) */
    float snow_to_ice_rate[PLASIC_NHOR];
    /* prognostic compactness c_i^* (restart xicecc); ice_cover is its 0/1 mask */
    float compactness[PLASIC_NHOR];
    /* Q_a and Q_a^mod = Q_a - Q_max */
    float atmospheric_heat[PLASIC_NHOR];
    float modified_atmospheric_heat[PLASIC_NHOR];
    /* Q_c = K_z (T_i - T_o); Q_a^mod on open water, 0 on thin ice */
    float conductive_flux[PLASIC_NHOR];
    /* Q_f = K_z (T_f - T_o) */
    float freezing_flux[PLASIC_NHOR];
    /* Q_net = Q_c - rho_w L_s P_s|open + Q_m^rem + Q_o->i */
    float net_heat_flux[PLASIC_NHOR];
    /* Q_i->o: flux handed to the slab */
    float ice_to_ocean_flux[PLASIC_NHOR];
    /* Q_max: ice above h_max, fed back next step through Q_a^mod */
    float maximum_thickness_flux[PLASIC_NHOR];
    /* Q_ex: melt energy beyond the remaining ice, returned to Q_i->o */
    float excess_melt_flux[PLASIC_NHOR];
    /* P_s (m water equivalent s^-1) */
    float snowfall[PLASIC_NHOR];
    float freshwater_flux[PLASIC_NHOR];
    float stress_x[PLASIC_NHOR];
    float stress_y[PLASIC_NHOR];
    float friction_velocity_cubed[PLASIC_NHOR];
    /* Q_o->i: returned by the slab, used in the next ice step */
    float ocean_to_ice_flux[PLASIC_NHOR];
    float surface_storage_flux[PLASIC_NHOR];
    float diagnosed_conductive_flux[PLASIC_NHOR];
    /* Q_s->i <= 0: latent heat of flooded snow */
    float snow_to_ice_flux[PLASIC_NHOR];
    float climatological_sst[
        PLASIC_NHOR * PLASIC_CLIMATOLOGY_MONTHS];
    float climatological_ice_cover[
        PLASIC_NHOR * PLASIC_CLIMATOLOGY_MONTHS];
    float climatological_ice_thickness[
        PLASIC_NHOR * PLASIC_CLIMATOLOGY_MONTHS];
    float current_climatological_sst[PLASIC_NHOR];
    float previous_climatological_sst[PLASIC_NHOR];
    float current_climatological_cover[PLASIC_NHOR];
    float current_climatological_thickness[PLASIC_NHOR];
    float current_flux_correction[PLASIC_NHOR];
    /*
     * Ocean-coupling sums over the coupling window, averaged and handed to
     * the slab by ocean_step(); ocean_heat_accumulator is the sum of Q_i->o
     * (restart record cheat).
     */
    float ocean_heat_accumulator[PLASIC_NHOR];
    float ocean_freshwater_accumulator[PLASIC_NHOR];
    float ocean_stress_x_accumulator[PLASIC_NHOR];
    float ocean_stress_y_accumulator[PLASIC_NHOR];
    float ocean_friction_accumulator[PLASIC_NHOR];
    float ocean_snow_accumulator[PLASIC_NHOR];
    /* Output-interval diagnostic sums of the fields above. */
    float flux_correction_accumulator[PLASIC_NHOR];
    float atmospheric_heat_accumulator[PLASIC_NHOR];
    float ocean_to_ice_flux_accumulator[PLASIC_NHOR];
    float melt_energy_accumulator[PLASIC_NHOR];
    float ice_to_ocean_flux_accumulator[PLASIC_NHOR];
    float maximum_thickness_flux_accumulator[PLASIC_NHOR];
    float excess_melt_flux_accumulator[PLASIC_NHOR];
    float snow_melt_accumulator[PLASIC_NHOR];
    float ice_energy_flux_accumulator[PLASIC_NHOR];
    float surface_storage_accumulator[PLASIC_NHOR];
    float diagnosed_conductive_accumulator[PLASIC_NHOR];
    float snow_to_ice_flux_accumulator[PLASIC_NHOR];
    float freshwater_diagnostic_accumulator[PLASIC_NHOR];
    float snow_to_ice_rate_accumulator[PLASIC_NHOR];
} ice_state;

typedef struct physics_workspace
{
    float level_a[PLASIC_NHOR * PLASIC_NLEV];
    float horizontal_a[PLASIC_NHOR];
    float horizontal_b[PLASIC_NHOR];
    float horizontal_c[PLASIC_NHOR];
    float horizontal_d[PLASIC_NHOR];
    float horizontal_e[PLASIC_NHOR];
    float horizontal_f[PLASIC_NHOR];
    float horizontal_g[PLASIC_NHOR];
    float horizontal_h[PLASIC_NHOR];
} physics_workspace;

extern misc_state misc;
extern flux_state flux;
extern rain_state rain;
extern radiation_state radiation;
extern land_state land;
extern sea_state sea;
extern ocean_state ocean;
extern deep_ocean_state deep_ocean;
extern ice_state ice;
extern physics_workspace *work;
extern char physics_error_message[512];

void physics_set_error(const char *format, ...);
void physics_interpolate_land_surface(void);
void physics_interpolate_ice_climatology(void);
void physics_interpolate_ocean_climatology(void);

void physics_ocean3d_initialize(int from_restart);
void physics_ocean3d_advance(void);

/*
 * Copy between the halo-extended deep-ocean layout and a contiguous local
 * slice suitable for restart I/O.  `levels` is the number of vertical levels.
 */
void physics_ocean3d_pack_local(
    const float *extended, float *local, int32_t levels);
void physics_ocean3d_unpack_local(
    const float *local, float *extended, int32_t levels);

int physics_initialize(
    const char *restart_path, const char *surface_data_directory,
    int32_t sea_coupling_interval);
void physics_finalize(void);
int physics_step(void);
int physics_write_restart(void);
int physics_frame_field(
    int32_t variable_id, const float **values, int32_t *layers);

#endif
