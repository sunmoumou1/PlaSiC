#ifndef PLASIC_RUNTIME_INTERNAL_H
#define PLASIC_RUNTIME_INTERNAL_H

#include "plasic_runtime.h"

#include "dynamical_core/plasic_kernels.h"
#include "math_tools/gaussmod_kernels.h"
#include "math_tools/sht_kernels.h"
#include "mpi/plasic_mpi.h"
#include "physics/plasic_physics.h"
#include "state/plasicmod_state.h"
#include "tools/calmod_kernels.h"
#include "tools/p_earth_kernels.h"
#include "tools/plasic_stream.h"
#include "tools/restartmod_io.h"
#include "tools/surface_data_io.h"
#include "tools/surfmod_kernels.h"

#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*
 * Private runtime state shared by the small runtime implementation units.
 */
typedef struct runtime_workspace
{
    float cold_start_orography_grid[PLASIC_NHOR];

    /* 绝热格点动力的诊断量和通量。 */
    /* Diagnostics and fluxes of the gridpoint adiabatic dynamics. */
    float zonal_log_surface_pressure_gradient[PLASIC_NHOR];
    float surface_pressure_ratio[PLASIC_NHOR];
    float adiabatic_temperature_source[PLASIC_NHOR * PLASIC_NLEV];
    float humidity_vertical_source[PLASIC_NHOR * PLASIC_NLEV];
    float zonal_momentum_flux[PLASIC_NHOR * PLASIC_NLEV];
    float meridional_momentum_flux[PLASIC_NHOR * PLASIC_NLEV];
    float column_pressure_advection[PLASIC_NHOR];
    float scaled_moist_geopotential[PLASIC_NHOR * PLASIC_NLEV];
    float zonal_temperature_flux[PLASIC_NHOR * PLASIC_NLEV];
    float meridional_temperature_flux[PLASIC_NHOR * PLASIC_NLEV];
    float scaled_kinetic_plus_geopotential[PLASIC_NHOR * PLASIC_NLEV];
    float zonal_humidity_flux[PLASIC_NHOR * PLASIC_NLEV];
    float meridional_humidity_flux[PLASIC_NHOR * PLASIC_NLEV];

    /* 球谐分析得到的全局显式 tendency；不同物理量不共享谱缓存。 */
    /* Globally explicit tendencies obtained from spherical-harmonic analysis;
     * different physical quantities do not share spectral buffers. */
    float column_pressure_advection_spectrum[PLASIC_NRSP];
    float divergence_tendency_spectrum[PLASIC_NRSP * PLASIC_NLEV];
    float temperature_tendency_spectrum[PLASIC_NRSP * PLASIC_NLEV];
    float vorticity_tendency_spectrum[PLASIC_NRSP * PLASIC_NLEV];
    float humidity_tendency_spectrum[PLASIC_NRSP * PLASIC_NLEV];

    /* 半隐式与跃蛙时间推进工作量。 */
    /* Working arrays for the semi-implicit and leapfrog time integration. */
    float saved_surface_log_pressure[PLASIC_NSPP];
    float saved_divergence[PLASIC_NSPP * PLASIC_NLEV];
    float saved_vorticity[PLASIC_NSPP * PLASIC_NLEV];
    float saved_temperature[PLASIC_NSPP * PLASIC_NLEV];
    float saved_humidity[PLASIC_NSPP * PLASIC_NLEV];
    float column_mass_divergence[PLASIC_NSPP];
    float centered_divergence[PLASIC_NSPP * PLASIC_NLEV];
    float complete_temperature_tendency[PLASIC_NSPP * PLASIC_NLEV];
    float nondimensional_reference_temperature[PLASIC_NLEV];
    float friction_rate[PLASIC_NLEV];
    float divergence_diffusion_rate[PLASIC_NLEV];
    float vorticity_diffusion_rate[PLASIC_NLEV];
    float temperature_diffusion_rate[PLASIC_NLEV];
    float humidity_diffusion_rate[PLASIC_NLEV];
    float divergence_friction_tendency[PLASIC_NSPP * PLASIC_NLEV];
    float divergence_diffusion_tendency[PLASIC_NSPP * PLASIC_NLEV];
    float vorticity_friction_tendency[PLASIC_NSPP * PLASIC_NLEV];
    float vorticity_diffusion_tendency[PLASIC_NSPP * PLASIC_NLEV];

    /* 输出诊断各自拥有与名称一致的存储。 */
    /* Output diagnostics each have storage matching their names. */
    float surface_geopotential[PLASIC_NHOR];
    float wind_speed[PLASIC_NHOR * PLASIC_NLEV];
    float total_precipitation_rate[PLASIC_NHOR];
    float convective_precipitation_rate[PLASIC_NHOR];
    float large_scale_precipitation_rate[PLASIC_NHOR];
    float snowfall_rate[PLASIC_NHOR];
    float total_cloud_cover[PLASIC_NHOR];
    float relative_humidity[PLASIC_NHOR * PLASIC_NLEV];
    float relative_vorticity[PLASIC_NHOR * PLASIC_NLEV];
    float physical_divergence[PLASIC_NHOR * PLASIC_NLEV];
    float geopotential_height[PLASIC_NHOR * PLASIC_NLEV];
} runtime_workspace;

extern runtime_workspace *runtime_fields;
extern FILE *progress_stream;
extern int progress_enabled;
extern int frames_enabled;
extern plasic_stream_writer frame_writer;
extern double integration_start_time;
extern char runtime_error_message[512];

void runtime_set_error(const char *format, ...);
int runtime_broadcast_status(int status);

int runtime_load_atmosphere_restart(
    const char *path, const char *surface_data_directory);
int runtime_prepare_orography_grid(void);
int runtime_initialize_dynamics(double step_minutes);

int runtime_gridpoint_adiabatic(void);
int runtime_spectral_adiabatic(void);
int runtime_gridpoint_diabatic(int dynamics_only);
int runtime_spectral_diabatic(void);

int runtime_state_is_finite(void);
void runtime_accumulate_output_fields(void);
int runtime_write_step_outputs(
    const plasic_runtime_options *options, int32_t completed_steps);
int runtime_open_live_outputs(const plasic_runtime_options *options);
int runtime_close_step_outputs(const plasic_runtime_options *options);

#endif
