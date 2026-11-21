#ifndef PLASIC_PHYSICS_H
#define PLASIC_PHYSICS_H

#include <stdint.h>

/* sea_coupling_interval: atmospheric steps per sea-ice/ocean step (>= 1). */
int plasic_physics_initialize(
    const char *restart_path, const char *surface_data_directory,
    int32_t sea_coupling_interval);
const char *plasic_physics_error(void);
int plasic_physics_step(void);
/* Write component records while restartmod_io already owns an output file. */
int plasic_physics_write_restart(void);
/*
 * Return a local latitude-slice field for the live frame stream.
 * The pointer remains owned by the physics runtime.
 */
int plasic_physics_frame_field(
    int32_t variable_id, const float **values, int32_t *layers);
int plasic_physics_set_co2_bands(const float *bands, int32_t band_count);
const float *plasic_physics_co2_field(void);
/*
 * Set radiation ozone to current + weight * (neighbor - current). Both are
 * local [lev][lat][lon] slabs of ozone mass fractions (kg/kg).
 */
void plasic_physics_interpolate_ozone(const double *current, const double *neighbor, float weight);
/*
 * Return the centre depth (metres) of every three-dimensional deep-ocean
 * level and write the level count to *count.  The pointer remains owned by
 * the physics runtime.
 */
const float *plasic_physics_ocean_depth(int32_t *count);
void plasic_physics_finalize(void);

#endif
