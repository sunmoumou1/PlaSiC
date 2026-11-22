#ifndef PLASIC_3DOCEAN_H
#define PLASIC_3DOCEAN_H

#include <stdint.h>

/*
 * PlaSiC three-dimensional deep-ocean kernels.
 *
 * The 3-D ocean lives below the single-layer mixed-layer / slab ocean.  The
 * slab remains the fast atmospheric boundary (SST seen by the atmosphere);
 * these kernels advance a genuine three-dimensional temperature, salinity and
 * velocity field on the model's Gaussian grid.
 *
 * Parallel layout
 * ---------------
 * The Gaussian grid is distributed over latitude rows.  Every process owns
 * `nlat_local` interior rows; the kernels operate on a halo-extended array
 * with one ghost row on each side, so meridional finite differences are
 * two-sided everywhere.  The physics layer fills the ghost rows by exchanging
 * the boundary rows with the neighbouring processes (see
 * mp_exchange_latitude_halo), and deep_ocean_step() calls the same exchange
 * through `exchange_halo` whenever it changes interior rows that a
 * neighbour's stencil reads.  Only the global polar ghosts are mirrored, so
 * the result does not depend on the number of processes.
 *
 * Extended storage is column-major:
 *
 *   index(horizontal, level) = horizontal + level * nhor_ext
 *   horizontal               = longitude + extended_row * nlon
 *   nhor_ext                 = nlon * (nlat_local + DEEP_OCEAN_HALO_ROWS)
 *
 * Extended row 0 is the north ghost, rows 1..nlat_local are the interior
 * rows (increasing southward), and row nlat_local+1 is the south ghost.
 * All temperatures are in Kelvin, salinities in g/kg, velocities in m/s,
 * distances in metres and fluxes in W/m^2.
 */

#define DEEP_OCEAN_HALO_ROWS 2
#define DEEP_OCEAN_EXT_NLAT(nlat_local) \
    ((nlat_local) + DEEP_OCEAN_HALO_ROWS)
#define DEEP_OCEAN_EXT_HOR(nlon, nlat_local) \
    ((nlon) * DEEP_OCEAN_EXT_NLAT(nlat_local))

typedef struct deep_ocean_kernel_config
{
    int32_t nlev;            /* number of vertical levels                    */
    int32_t nlon;            /* longitudes per latitude row                  */
    int32_t nlat_local;      /* interior latitude rows owned by this process */
    float dt;                /* ocean time step (s)                          */
    float radius;            /* planetary radius (m)                         */
    float solar_day;         /* length of the model solar day (s)            */
    float sidereal_day;      /* length of the sidereal day (s)               */
    float gravity;           /* gravitational acceleration (m/s^2)           */
    float cp;                /* seawater heat capacity (J/kg/K)              */
    float rho0;              /* reference seawater density (kg/m^3)          */
    float alpha;             /* thermal expansion coefficient (1/K)          */
    float beta;              /* haline contraction coefficient (1/(g/kg))    */
    float kappa_v;           /* vertical diffusivity (m^2/s)                 */
    float kappa_h;           /* horizontal diffusivity (m^2/s)               */
    float t_freeze;          /* freezing point of seawater (K)               */
    /*
     * 刷新 halo 扩展场前 `levels` 层的 ghost 行（层间跨度为 nlon *
     * (nlat_local + DEEP_OCEAN_HALO_ROWS)）：进程边界处的 ghost 行取相邻进程
     * 的边缘行，全球极点处的 ghost 行取相邻内部行的镜像。成功返回 0。
     * 为 NULL 表示单个进程拥有全部纬度，所有 ghost 行都取镜像。
     *
     * Refresh the ghost rows of the first `levels` levels of a halo-extended
     * field (level stride nlon * (nlat_local + DEEP_OCEAN_HALO_ROWS)): a
     * ghost row at a process boundary receives the neighbouring process's
     * edge row, a ghost row at a global pole mirrors the adjacent interior
     * row.  Returns 0 on success.  NULL means one process owns every
     * latitude, so every ghost row is mirrored.
     */
    int (*exchange_halo)(float *field, int32_t levels);
} deep_ocean_kernel_config;

/*
 * Build the default vertical grid: the thickness (metres) of each layer and
 * the depth (metres) of its centre.  The top layer matches the standard 50 m
 * mixed layer so the slab and the 3-D ocean share a surface layer.  The mesh
 * is generated for any user-selected `nlev`.
 */
void deep_ocean_build_levels(
    int32_t nlev, float *layer_thickness, float *layer_center_depth);

/*
 * Seed the halo-extended 3-D ocean.  Interior rows start from the supplied
 * surface temperature with an exponentially decaying stratification towards a
 * cold abyss; ghost rows are mirrored from the adjacent interior row.
 */
void deep_ocean_seed_state(
    const deep_ocean_kernel_config *config, const float *land_mask,
    const float *surface_temperature, const float *layer_thickness,
    float *temperature, float *salinity, float *u_velocity,
    float *v_velocity, float *w_velocity);

/*
 * Advance the halo-extended 3-D ocean by one time step.
 *
 * The top layer is first reset to the slab temperature (the atmosphere-facing
 * mixed layer), exactly as the reference LSG surface forcing does.  The step
 * then computes density and diagnostic geostrophic plus wind-driven
 * velocities, advects and diffuses temperature and salinity, applies
 * convective adjustment and finally reports the heat flux that the deep ocean
 * exchanges with the slab in `deep_ocean_flux` (Q_deep in W m^-2, positive
 * into the slab; the caller passes ocean.deep_ocean_flux).  Land points get
 * Q_deep = 0 and are excluded from every horizontal stencil through masked,
 * zero-flux coastal gradients.
 *
 * `sine_latitude` and `cosine_latitude` are local, halo-extended latitude
 * coordinates of length nlat_local + DEEP_OCEAN_HALO_ROWS; `slab_temperature`,
 * `freshwater_flux`, `stress_x`, `stress_y` and `deep_ocean_flux` are
 * contiguous local slices of the process's interior rows.
 */
void deep_ocean_step(
    const deep_ocean_kernel_config *config, const float *sine_latitude,
    const float *cosine_latitude, const float *land_mask,
    const float *slab_temperature, const float *freshwater_flux,
    const float *stress_x, const float *stress_y,
    const float *layer_thickness, float *temperature, float *salinity,
    float *u_velocity, float *v_velocity, float *w_velocity,
    float *density, float *deep_ocean_flux);

#endif
