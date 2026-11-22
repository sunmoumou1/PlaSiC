#include "3docean.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

/*
 * PlaSiC three-dimensional deep-ocean kernels.
 *
 * The numerical design is deliberately conservative: an explicit, upwind,
 * flux-limited advection-diffusion scheme with convective adjustment and a
 * diagnostic geostrophic plus wind-driven velocity field.  It is stable at
 * the atmospheric coupling time step and produces a genuine baroclinic
 * three-dimensional ocean without requiring a barotropic solver.  The
 * reference PLASIM LSG model is used as the physical guide.
 *
 * Arrays are halo-extended in latitude so that every meridional difference is
 * two-sided and every process computes only its own latitude band.  Coastal
 * gradients use the land mask so that land points never leak into the ocean
 * stencils.
 */

enum
{
    DEEP_OCEAN_NSUB_MAX = 16
};

static size_t field_index(int32_t horizontal, int32_t level, int32_t nhor)
{
    return (size_t)horizontal + (size_t)level * (size_t)nhor;
}

static float clamp_float(float value, float lower, float upper)
{
    return fmaxf(lower, fminf(upper, value));
}

static float seawater_density(
    const deep_ocean_kernel_config *config, float temperature,
    float salinity)
{
    const float celsius = temperature - 273.15f;
    return config->rho0 *
           (1.0f - config->alpha * (celsius - 4.0f) +
            config->beta * (salinity - 35.0f));
}

static int32_t extended_nhor(const deep_ocean_kernel_config *config)
{
    return config->nlon * DEEP_OCEAN_EXT_NLAT(config->nlat_local);
}

static int is_land(const float *land_mask, int32_t horizontal)
{
    return land_mask[horizontal] >= 0.5f;
}

/*
 * 刷新 halo 扩展场前 `levels` 层的 ghost 行。并行时由物理层的交换钩子把
 * 相邻进程的边缘行写进 ghost 行（极点处保留镜像）；单进程或交换失败时，
 * 所有 ghost 行取相邻内部行的镜像（零梯度）。
 *
 * 以前这里对所有进程都只做镜像：进程边界处表层的 ghost 行被本进程的边缘行
 * 覆盖，经向梯度被人为置零，表层的经向平流和扩散在进程边界处被切断。
 * 4 进程的 T85 积分因此在 45°N/45°S 出现约 2 K 的海温断层。
 *
 * Refresh the ghost rows of the first `levels` levels of a halo-extended
 * field.  In parallel the physics layer's exchange hook writes the
 * neighbouring process's edge row into each ghost row (the global polar
 * ghosts stay mirrored); with one process, or if the exchange fails, every
 * ghost row mirrors its adjacent interior row (zero gradient).
 *
 * The step used to mirror every ghost row on every process.  At a process
 * boundary that overwrote the top-layer ghost row with the process's own edge
 * row, set the meridional gradient there to zero and cut the top layer's
 * meridional advection and diffusion at the boundary; a four-process T85
 * integration developed SST steps of about 2 K at 45N and 45S.
 */
static void refresh_ghost_rows(
    const deep_ocean_kernel_config *config, float *field, int32_t levels)
{
    const int32_t nlon = config->nlon;
    const int32_t nlat_local = config->nlat_local;
    const int32_t nhor = extended_nhor(config);
    int32_t level;
    int32_t i;

    if (config->exchange_halo != NULL &&
        config->exchange_halo(field, levels) == 0)
    {
        return;
    }
    for (level = 0; level < levels; ++level)
    {
        for (i = 0; i < nlon; ++i)
        {
            field[field_index(i, level, nhor)] =
                field[field_index(nlon + i, level, nhor)];
            field[field_index((nlat_local + 1) * nlon + i, level, nhor)] =
                field[field_index(nlat_local * nlon + i, level, nhor)];
        }
    }
}

void deep_ocean_build_levels(
    int32_t nlev, float *layer_thickness, float *layer_center_depth)
{
    float depth = 0.0f;
    int32_t level;

    for (level = 0; level < nlev; ++level)
    {
        /*
         * A 50 m surface layer, thickening linearly with depth.  For the
         * default 16 levels this sums to about 5000 m, close to the global
         * mean ocean depth.
         */
        const float thickness =
            level == 0 ? 50.0f : 50.0f + 40.0f * (float)(level - 1);
        layer_thickness[level] = thickness;
        layer_center_depth[level] = depth + 0.5f * thickness;
        depth += thickness;
    }
}

/*
 * Meridional grid spacing (metres, positive northwards) at an interior
 * extended row.  The latitude coordinates are local and halo-extended, so the
 * centred difference is always available; at the global poles the physics
 * layer supplies a linearly extrapolated ghost coordinate.
 */
static float latitude_dy(
    const deep_ocean_kernel_config *config, const float *sine_latitude,
    int32_t extended_row, float cosine_latitude)
{
    float delta_sine = 0.5f * fabsf(
        sine_latitude[extended_row - 1] - sine_latitude[extended_row + 1]);

    if (delta_sine < 1.0e-7f)
    {
        delta_sine = 1.0e-7f;
    }
    return config->radius * delta_sine / fmaxf(1.0e-4f, cosine_latitude);
}

static int32_t clamped_dy_row(
    const deep_ocean_kernel_config *config, int32_t extended_row)
{
    if (extended_row < 1)
    {
        return 1;
    }
    if (extended_row > config->nlat_local)
    {
        return config->nlat_local;
    }
    return extended_row;
}

void deep_ocean_seed_state(
    const deep_ocean_kernel_config *config, const float *land_mask,
    const float *surface_temperature, const float *layer_thickness,
    float *temperature, float *salinity, float *u_velocity,
    float *v_velocity, float *w_velocity)
{
    const int32_t nlon = config->nlon;
    const int32_t nlat_local = config->nlat_local;
    const int32_t nhor = extended_nhor(config);
    const float abyssal_temperature = 275.15f; /* 2 degrees Celsius */
    /*
     * 比 2 °C 冷的水柱整柱取海表温度（深层温度取 min(2 °C, SST)）。若深层仍是
     * 2 °C，冷而重的表层水压在暖而轻的水上；slab 每步又把表层重置到冰点，于是
     * 对流不断把深层热量翻到表层，以 2000 W/m² 的上限加热 slab，初始海冰在
     * 几天内融光，之后至少一年都无法结冰。等温的冷水柱靠参考盐度廓线保持稳定。
     *
     * Columns colder than 2 degC are seeded isothermal at the surface
     * temperature (deep temperature min(2 degC, SST)). With 2 degC water below,
     * cold dense surface water would sit on warm light water; the slab resets
     * the top layer to its own (freezing) temperature every step, so
     * convective adjustment would keep bringing deep heat up and feed the slab
     * at the 2000 W/m2 flux limit, melting the initial sea ice within days and
     * keeping polar seas ice-free for at least a year. An isothermal cold
     * column is held stable by the reference salinity profile.
     */
    float depth = 0.0f;
    int32_t level;

    for (level = 0; level < config->nlev; ++level)
    {
        const float thickness = layer_thickness[level];
        const float center_depth = depth + 0.5f * thickness;
        const float profile = expf(-center_depth / 500.0f);
        const float salt_profile =
            34.7f + 0.2f * (1.0f - expf(-center_depth / 1000.0f));
        int32_t row;
        int32_t i;

        for (row = 0; row < nlat_local; ++row)
        {
            const int32_t local_row = row * nlon;
            const int32_t extended_row = (row + 1) * nlon;

            for (i = 0; i < nlon; ++i)
            {
                const int32_t local = local_row + i;
                const int32_t extended = extended_row + i;
                const size_t index = field_index(extended, level, nhor);
                const float surface = surface_temperature[local];
                const float deep = fminf(abyssal_temperature, surface);

                temperature[index] = deep + (surface - deep) * profile;
                salinity[index] = salt_profile;
                u_velocity[index] = 0.0f;
                v_velocity[index] = 0.0f;
                w_velocity[index] = 0.0f;
            }
        }
        /*
         * Mirror the innermost interior rows into the vertical ghost rows.
         * The physics layer replaces these with true neighbour values during
         * the first halo exchange.
         */
        for (i = 0; i < nlon; ++i)
        {
            const size_t north_ghost = field_index(i, level, nhor);
            const size_t north_edge = field_index(nlon + i, level, nhor);
            const size_t south_edge =
                field_index(nlat_local * nlon + i, level, nhor);
            const size_t south_ghost =
                field_index((nlat_local + 1) * nlon + i, level, nhor);

            temperature[north_ghost] = temperature[north_edge];
            salinity[north_ghost] = salinity[north_edge];
            u_velocity[north_ghost] = u_velocity[north_edge];
            v_velocity[north_ghost] = v_velocity[north_edge];
            w_velocity[north_ghost] = w_velocity[north_edge];
            temperature[south_ghost] = temperature[south_edge];
            salinity[south_ghost] = salinity[south_edge];
            u_velocity[south_ghost] = u_velocity[south_edge];
            v_velocity[south_ghost] = v_velocity[south_edge];
            w_velocity[south_ghost] = w_velocity[south_edge];
        }
        depth += thickness;
    }
    (void)land_mask;
}

static float latitude_dx(
    const deep_ocean_kernel_config *config, float cosine_latitude)
{
    const float coslat = fmaxf(1.0e-4f, cosine_latitude);
    const float two_pi = 6.28318530717958647692f;
    return config->radius * coslat * two_pi / (float)config->nlon;
}

/*
 * Masked zonal gradient.  When a neighbour is land the difference is taken
 * one-sided against the remaining ocean neighbour; when both sides are land
 * the gradient is zero, so coastlines never inject land values.
 */
static float zonal_gradient(
    const float *field, const float *land_mask, int32_t horizontal,
    int32_t level, int32_t nhor, int32_t nlon, float dx)
{
    const int32_t i = horizontal % nlon;
    const int32_t row = horizontal - i;
    const int32_t left = row + (i == 0 ? nlon - 1 : i - 1);
    const int32_t right = row + (i == nlon - 1 ? 0 : i + 1);
    const size_t base = (size_t)level * (size_t)nhor;
    const float self = field[base + (size_t)horizontal];

    if (is_land(land_mask, left) && is_land(land_mask, right))
    {
        return 0.0f;
    }
    if (is_land(land_mask, left))
    {
        return (field[base + (size_t)right] - self) / dx;
    }
    if (is_land(land_mask, right))
    {
        return (self - field[base + (size_t)left]) / dx;
    }
    return (field[base + (size_t)right] - field[base + (size_t)left]) /
           (2.0f * dx);
}

/*
 * Masked meridional gradient.  Global latitude indices increase southward, so
 * `north` means the lower extended-row index and `dy` is positive northwards.
 */
static float meridional_gradient(
    const float *field, const float *land_mask, int32_t horizontal,
    int32_t level, int32_t nhor, int32_t nlon, float dy)
{
    const size_t base = (size_t)level * (size_t)nhor;
    const int32_t north = horizontal - nlon;
    const int32_t south = horizontal + nlon;
    const float self = field[base + (size_t)horizontal];

    if (is_land(land_mask, north) && is_land(land_mask, south))
    {
        return 0.0f;
    }
    if (is_land(land_mask, north))
    {
        return (self - field[base + (size_t)south]) / dy;
    }
    if (is_land(land_mask, south))
    {
        return (field[base + (size_t)north] - self) / dy;
    }
    return (field[base + (size_t)north] - field[base + (size_t)south]) /
           (2.0f * dy);
}

/*
 * Diagnostic velocity field: geostrophic balance from the pressure gradient
 * plus a surface Ekman component driven by the wind stress.  The Coriolis
 * parameter is regularised near the equator and all speeds are capped so the
 * explicit advection remains stable.
 */
static void diagnose_velocities(
    const deep_ocean_kernel_config *config, const float *sine_latitude,
    const float *cosine_latitude, const float *land_mask,
    const float *pressure, const float *stress_x, const float *stress_y,
    const float *layer_thickness, float *u_velocity, float *v_velocity)
{
    const int32_t nlon = config->nlon;
    const int32_t nhor = extended_nhor(config);
    const int32_t rows = DEEP_OCEAN_EXT_NLAT(config->nlat_local);
    const float omega = 6.28318530717958647692f / config->sidereal_day;
    const float minimum_sine = sinf(2.0f * 3.14159265359f / 180.0f);
    const float speed_cap = 0.5f;
    int32_t extended_row;

    for (extended_row = 0; extended_row < rows; ++extended_row)
    {
        const int32_t dy_row = clamped_dy_row(config, extended_row);
        const float coslat = cosine_latitude[extended_row];
        const float dx = latitude_dx(config, coslat);
        const float dy = latitude_dy(config, sine_latitude, dy_row, coslat);
        const int is_interior =
            extended_row >= 1 && extended_row <= config->nlat_local;
        const int32_t local_row = extended_row - 1;
        float coriolis =
            2.0f * omega * sine_latitude[extended_row];
        int32_t i;
        int32_t level;

        if (!is_interior)
        {
            continue;
        }
        if (fabsf(coriolis) < 2.0f * omega * minimum_sine)
        {
            coriolis = coriolis < 0.0f ?
                           -2.0f * omega * minimum_sine :
                           2.0f * omega * minimum_sine;
        }
        for (i = 0; i < nlon; ++i)
        {
            const int32_t horizontal = extended_row * nlon + i;
            const float h0 = fmaxf(1.0f, layer_thickness[0]);
            float ekman_u = 0.0f;
            float ekman_v = 0.0f;

            if (is_land(land_mask, horizontal))
            {
                for (level = 0; level < config->nlev; ++level)
                {
                    const size_t index = field_index(horizontal, level, nhor);
                    u_velocity[index] = 0.0f;
                    v_velocity[index] = 0.0f;
                }
                continue;
            }
            if (is_interior)
            {
                const int32_t local = local_row * nlon + i;
                ekman_u =
                    stress_y[local] / (config->rho0 * coriolis * h0);
                ekman_v =
                    -stress_x[local] / (config->rho0 * coriolis * h0);
            }
            for (level = 0; level < config->nlev; ++level)
            {
                const size_t index = field_index(horizontal, level, nhor);
                float u_value = -meridional_gradient(
                                    pressure, land_mask, horizontal, level,
                                    nhor, nlon, dy) /
                                (config->rho0 * coriolis);
                float v_value = zonal_gradient(
                                    pressure, land_mask, horizontal, level,
                                    nhor, nlon, dx) /
                                (config->rho0 * coriolis);

                if (level == 0 && is_interior)
                {
                    u_value += ekman_u;
                    v_value += ekman_v;
                }
                u_velocity[index] =
                    clamp_float(u_value, -speed_cap, speed_cap);
                v_velocity[index] =
                    clamp_float(v_value, -speed_cap, speed_cap);
            }
        }
    }
}

/*
 * Vertical velocity diagnosed from the horizontal divergence so that the
 * three-dimensional flow is non-divergent in the column-integrated sense.
 */
static void diagnose_vertical_velocity(
    const deep_ocean_kernel_config *config, const float *sine_latitude,
    const float *cosine_latitude, const float *land_mask,
    const float *layer_thickness, const float *u_velocity,
    const float *v_velocity, float *w_velocity)
{
    const int32_t nlon = config->nlon;
    const int32_t nhor = extended_nhor(config);
    const float speed_cap = 1.0e-6f;
    int32_t row;
    int32_t i;

    for (row = 0; row < config->nlat_local; ++row)
    {
        const int32_t extended_row = row + 1;
        const float coslat = cosine_latitude[extended_row];
        const float dx = latitude_dx(config, coslat);
        const float dy =
            latitude_dy(config, sine_latitude, extended_row, coslat);

        for (i = 0; i < nlon; ++i)
        {
            const int32_t horizontal = extended_row * nlon + i;
            float w_value = 0.0f;
            int32_t level;

            w_velocity[field_index(horizontal, 0, nhor)] = 0.0f;
            if (is_land(land_mask, horizontal))
            {
                for (level = 1; level < config->nlev; ++level)
                {
                    w_velocity[field_index(horizontal, level, nhor)] = 0.0f;
                }
                continue;
            }
            for (level = 1; level < config->nlev; ++level)
            {
                const float divergence =
                    zonal_gradient(
                        u_velocity, land_mask, horizontal, level - 1, nhor,
                        nlon, dx) +
                    meridional_gradient(
                        v_velocity, land_mask, horizontal, level - 1, nhor,
                        nlon, dy);
                w_value -= divergence * layer_thickness[level - 1];
                w_value = clamp_float(w_value, -speed_cap, speed_cap);
                w_velocity[field_index(horizontal, level, nhor)] = w_value;
            }
        }
    }
}

/*
 * Explicit advection plus diffusion tendency for temperature and salinity.
 * Advection uses first-order upwinding with zero flux across coastlines;
 * diffusion is a masked Laplacian with a zero-gradient condition at land.
 */
static void tracer_tendencies(
    const deep_ocean_kernel_config *config, const float *sine_latitude,
    const float *cosine_latitude, const float *land_mask,
    const float *layer_thickness, const float *temperature,
    const float *salinity, const float *u_velocity,
    const float *v_velocity, const float *w_velocity, float *tendency_t,
    float *tendency_s)
{
    const int32_t nlon = config->nlon;
    const int32_t nhor = extended_nhor(config);
    int32_t row;
    int32_t i;

    for (row = 0; row < config->nlat_local; ++row)
    {
        const int32_t extended_row = row + 1;
        const float coslat = cosine_latitude[extended_row];
        const float dx = latitude_dx(config, coslat);
        const float dy =
            latitude_dy(config, sine_latitude, extended_row, coslat);

        for (i = 0; i < nlon; ++i)
        {
            const int32_t horizontal = extended_row * nlon + i;
            int32_t level;

            if (is_land(land_mask, horizontal))
            {
                for (level = 0; level < config->nlev; ++level)
                {
                    const size_t index =
                        field_index(horizontal, level, nhor);
                    tendency_t[index] = 0.0f;
                    tendency_s[index] = 0.0f;
                }
                continue;
            }
            for (level = 0; level < config->nlev; ++level)
            {
                const size_t index = field_index(horizontal, level, nhor);
                const int32_t row_base = horizontal - i;
                const int32_t east =
                    row_base + (i == nlon - 1 ? 0 : i + 1);
                const int32_t west =
                    row_base + (i == 0 ? nlon - 1 : i - 1);
                const size_t north_index =
                    field_index(horizontal - nlon, level, nhor);
                const size_t south_index =
                    field_index(horizontal + nlon, level, nhor);
                const size_t east_index = field_index(east, level, nhor);
                const size_t west_index = field_index(west, level, nhor);
                const int east_land = is_land(land_mask, east);
                const int west_land = is_land(land_mask, west);
                const int north_land = is_land(land_mask, horizontal - nlon);
                const int south_land = is_land(land_mask, horizontal + nlon);
                const float u_value = u_velocity[index];
                const float v_value = v_velocity[index];
                const float w_value = w_velocity[index];
                const float spacing_above =
                    level > 0 ?
                        0.5f * (layer_thickness[level] +
                                layer_thickness[level - 1]) :
                        1.0f;
                const float spacing_below =
                    level < config->nlev - 1 ?
                        0.5f * (layer_thickness[level] +
                                layer_thickness[level + 1]) :
                        1.0f;
                float d_temperature_dx;
                float d_temperature_dy;
                float d_temperature_dz;
                float d_salinity_dx;
                float d_salinity_dy;
                float d_salinity_dz;
                float east_temperature;
                float west_temperature;
                float north_temperature;
                float south_temperature;
                float east_salinity;
                float west_salinity;
                float north_salinity;
                float south_salinity;
                float laplacian_t;
                float laplacian_s;

                /* Zonal upwind advection: no inflow across a coastline. */
                if (u_value >= 0.0f)
                {
                    d_temperature_dx =
                        west_land ?
                            0.0f :
                            (temperature[index] -
                             temperature[west_index]) /
                                dx;
                    d_salinity_dx =
                        west_land ?
                            0.0f :
                            (salinity[index] - salinity[west_index]) / dx;
                }
                else
                {
                    d_temperature_dx =
                        east_land ?
                            0.0f :
                            (temperature[east_index] -
                             temperature[index]) /
                                dx;
                    d_salinity_dx =
                        east_land ?
                            0.0f :
                            (salinity[east_index] - salinity[index]) / dx;
                }
                if (v_value >= 0.0f)
                {
                    d_temperature_dy =
                        south_land ?
                            0.0f :
                            (temperature[index] -
                             temperature[south_index]) /
                                dy;
                    d_salinity_dy =
                        south_land ?
                            0.0f :
                            (salinity[index] - salinity[south_index]) / dy;
                }
                else
                {
                    d_temperature_dy =
                        north_land ?
                            0.0f :
                            (temperature[north_index] -
                             temperature[index]) /
                                dy;
                    d_salinity_dy =
                        north_land ?
                            0.0f :
                            (salinity[north_index] - salinity[index]) / dy;
                }
                if (w_value >= 0.0f)
                {
                    d_temperature_dz =
                        level > 0 ?
                            (temperature[index] -
                             temperature[field_index(
                                 horizontal, level - 1, nhor)]) /
                                spacing_above :
                            0.0f;
                    d_salinity_dz =
                        level > 0 ?
                            (salinity[index] -
                             salinity[field_index(
                                 horizontal, level - 1, nhor)]) /
                                spacing_above :
                            0.0f;
                }
                else
                {
                    d_temperature_dz =
                        level < config->nlev - 1 ?
                            (temperature[field_index(
                                 horizontal, level + 1, nhor)] -
                             temperature[index]) /
                                spacing_below :
                            0.0f;
                    d_salinity_dz =
                        level < config->nlev - 1 ?
                            (salinity[field_index(
                                 horizontal, level + 1, nhor)] -
                             salinity[index]) /
                                spacing_below :
                            0.0f;
                }

                east_temperature =
                    east_land ? temperature[index]
                              : temperature[east_index];
                west_temperature =
                    west_land ? temperature[index]
                              : temperature[west_index];
                north_temperature =
                    north_land ? temperature[index]
                               : temperature[north_index];
                south_temperature =
                    south_land ? temperature[index]
                               : temperature[south_index];
                east_salinity =
                    east_land ? salinity[index] : salinity[east_index];
                west_salinity =
                    west_land ? salinity[index] : salinity[west_index];
                north_salinity =
                    north_land ? salinity[index] : salinity[north_index];
                south_salinity =
                    south_land ? salinity[index] : salinity[south_index];

                laplacian_t =
                    (east_temperature - 2.0f * temperature[index] +
                     west_temperature) /
                        (dx * dx) +
                    (south_temperature - 2.0f * temperature[index] +
                     north_temperature) /
                        (dy * dy);
                laplacian_s =
                    (east_salinity - 2.0f * salinity[index] +
                     west_salinity) /
                        (dx * dx) +
                    (south_salinity - 2.0f * salinity[index] +
                     north_salinity) /
                        (dy * dy);

                tendency_t[index] =
                    -(u_value * d_temperature_dx +
                      v_value * d_temperature_dy +
                      w_value * d_temperature_dz) +
                    config->kappa_h * laplacian_t;
                tendency_s[index] =
                    -(u_value * d_salinity_dx +
                      v_value * d_salinity_dy +
                      w_value * d_salinity_dz) +
                    config->kappa_h * laplacian_s;

                if (level > 0)
                {
                    tendency_t[index] +=
                        config->kappa_v *
                        (temperature[field_index(horizontal, level - 1, nhor)] -
                         temperature[index]) /
                        spacing_above / layer_thickness[level];
                    tendency_s[index] +=
                        config->kappa_v *
                        (salinity[field_index(horizontal, level - 1, nhor)] -
                         salinity[index]) /
                        spacing_above / layer_thickness[level];
                }
                if (level < config->nlev - 1)
                {
                    tendency_t[index] +=
                        config->kappa_v *
                        (temperature[field_index(horizontal, level + 1, nhor)] -
                         temperature[index]) /
                        spacing_below / layer_thickness[level];
                    tendency_s[index] +=
                        config->kappa_v *
                        (salinity[field_index(horizontal, level + 1, nhor)] -
                         salinity[index]) /
                        spacing_below / layer_thickness[level];
                }
            }
        }
    }
}

static void apply_convective_adjustment(
    const deep_ocean_kernel_config *config, const float *land_mask,
    const float *layer_thickness, float *temperature, float *salinity,
    float *density)
{
    const int32_t nlon = config->nlon;
    const int32_t nhor = extended_nhor(config);
    const int32_t nlev = config->nlev;
    int32_t row;
    int32_t i;

    for (row = 0; row < config->nlat_local; ++row)
    {
        const int32_t extended_row = row + 1;

        for (i = 0; i < nlon; ++i)
        {
            const int32_t horizontal = extended_row * nlon + i;
            int pass;

            if (is_land(land_mask, horizontal))
            {
                continue;
            }
            for (pass = 0; pass < nlev; ++pass)
            {
                int mixed = 0;
                int32_t level;

                for (level = 1; level < nlev; ++level)
                {
                    const size_t upper =
                        field_index(horizontal, level - 1, nhor);
                    const size_t lower =
                        field_index(horizontal, level, nhor);
                    float h_upper;
                    float h_lower;
                    float h_total;

                    if (salinity[upper] <= 0.0f || salinity[lower] <= 0.0f)
                    {
                        continue;
                    }
                    if (density[lower] >= density[upper] - 1.0e-4f)
                    {
                        continue;
                    }
                    h_upper = layer_thickness[level - 1];
                    h_lower = layer_thickness[level];
                    h_total = h_upper + h_lower;
                    temperature[upper] =
                        (temperature[upper] * h_upper +
                         temperature[lower] * h_lower) /
                        h_total;
                    temperature[lower] = temperature[upper];
                    salinity[upper] =
                        (salinity[upper] * h_upper +
                         salinity[lower] * h_lower) /
                        h_total;
                    salinity[lower] = salinity[upper];
                    density[upper] = seawater_density(
                        config, temperature[upper], salinity[upper]);
                    density[lower] = density[upper];
                    mixed = 1;
                }
                if (mixed == 0)
                {
                    break;
                }
            }
        }
    }
}

void deep_ocean_step(
    const deep_ocean_kernel_config *config, const float *sine_latitude,
    const float *cosine_latitude, const float *land_mask,
    const float *slab_temperature, const float *freshwater_flux,
    const float *stress_x, const float *stress_y,
    const float *layer_thickness, float *temperature, float *salinity,
    float *u_velocity, float *v_velocity, float *w_velocity,
    float *density, float *deep_ocean_flux)
{
    const int32_t nlon = config->nlon;
    const int32_t nhor = extended_nhor(config);
    const int32_t nlev = config->nlev;
    float *pressure = NULL;
    float *tendency_t = NULL;
    float *tendency_s = NULL;
    int32_t row;
    int32_t i;
    int32_t level;
    int step_count;
    int substep;

    pressure = (float *)malloc(
        (size_t)nhor * (size_t)nlev * sizeof(*pressure));
    tendency_t = (float *)malloc(
        (size_t)nhor * (size_t)nlev * sizeof(*tendency_t));
    tendency_s = (float *)malloc(
        (size_t)nhor * (size_t)nlev * sizeof(*tendency_s));
    if (pressure == NULL || tendency_t == NULL || tendency_s == NULL)
    {
        free(pressure);
        free(tendency_t);
        free(tendency_s);
        return;
    }

    /* 1. Surface boundary: the slab sets the interior top-layer temperature
     *    and the net freshwater flux dilutes or concentrates the surface
     *    salinity.  The new top-layer edge rows are then exchanged, so the
     *    ghost rows hold the neighbouring process's top layer. */
    for (row = 0; row < config->nlat_local; ++row)
    {
        const int32_t extended_row = row + 1;
        const int32_t local_row = row * nlon;
        const int32_t ext_row = extended_row * nlon;

        for (i = 0; i < nlon; ++i)
        {
            const int32_t horizontal = ext_row + i;
            const int32_t local = local_row + i;
            const float h0 = fmaxf(1.0f, layer_thickness[0]);
            const float freshwater = freshwater_flux[local] * config->dt;

            if (is_land(land_mask, horizontal))
            {
                continue;
            }
            temperature[field_index(horizontal, 0, nhor)] =
                slab_temperature[local];
            if (h0 + freshwater > 1.0f)
            {
                const float top_salinity =
                    salinity[field_index(horizontal, 0, nhor)];
                salinity[field_index(horizontal, 0, nhor)] =
                    top_salinity * h0 / (h0 + freshwater);
            }
            {
                const size_t top = field_index(horizontal, 0, nhor);
                salinity[top] = clamp_float(salinity[top], 0.0f, 60.0f);
                temperature[top] =
                    fmaxf(config->t_freeze, temperature[top]);
            }
        }
    }
    refresh_ghost_rows(config, temperature, 1);
    refresh_ghost_rows(config, salinity, 1);

    /* 2. Equation of state and hydrostatic pressure on all extended rows. */
    for (row = 0; row < DEEP_OCEAN_EXT_NLAT(config->nlat_local); ++row)
    {
        for (i = 0; i < nlon; ++i)
        {
            const int32_t horizontal = row * nlon + i;
            float column = 0.0f;

            for (level = 0; level < nlev; ++level)
            {
                const size_t index = field_index(horizontal, level, nhor);
                float rho;

                if (is_land(land_mask, horizontal))
                {
                    density[index] = config->rho0;
                    pressure[index] = 0.0f;
                    continue;
                }
                rho = seawater_density(
                    config, temperature[index], salinity[index]);
                density[index] = rho;
                column += rho * layer_thickness[level];
                pressure[index] =
                    config->gravity * column -
                    0.5f * config->gravity * rho *
                        layer_thickness[level];
            }
        }
    }

    /* 3. Diagnostic velocities and column continuity. */
    diagnose_velocities(
        config, sine_latitude, cosine_latitude, land_mask, pressure,
        stress_x, stress_y, layer_thickness, u_velocity, v_velocity);
    /* The divergence at an edge row needs the neighbour's new velocity. */
    refresh_ghost_rows(config, u_velocity, nlev);
    refresh_ghost_rows(config, v_velocity, nlev);
    diagnose_vertical_velocity(
        config, sine_latitude, cosine_latitude, land_mask, layer_thickness,
        u_velocity, v_velocity, w_velocity);

    /* 4. Sub-cycled advection-diffusion for temperature and salinity. */
    step_count = 1;
    {
        const float cap_speed = 0.5f;
        const float cap_vertical = 1.0e-6f;
        const float smallest_spacing = 50.0f;
        const float rate =
            cap_speed / 1.0e5f + cap_vertical / smallest_spacing;

        step_count = 1 + (int)(rate * config->dt * 4.0f);
        if (step_count > DEEP_OCEAN_NSUB_MAX)
        {
            step_count = DEEP_OCEAN_NSUB_MAX;
        }
        if (step_count < 1)
        {
            step_count = 1;
        }
    }
    for (substep = 0; substep < step_count; ++substep)
    {
        const float sub_dt = config->dt / (float)step_count;

        tracer_tendencies(
            config, sine_latitude, cosine_latitude, land_mask,
            layer_thickness, temperature, salinity, u_velocity,
            v_velocity, w_velocity, tendency_t, tendency_s);
        for (row = 0; row < config->nlat_local; ++row)
        {
            const int32_t extended_row = row + 1;

            for (i = 0; i < nlon; ++i)
            {
                const int32_t horizontal = extended_row * nlon + i;

                if (is_land(land_mask, horizontal))
                {
                    continue;
                }
                for (level = 0; level < nlev; ++level)
                {
                    const size_t index =
                        field_index(horizontal, level, nhor);
                    temperature[index] = clamp_float(
                        temperature[index] + sub_dt * tendency_t[index],
                        config->t_freeze, 320.0f);
                    salinity[index] = clamp_float(
                        salinity[index] + sub_dt * tendency_s[index], 0.0f,
                        60.0f);
                }
            }
        }
        /* Refresh every level's ghost tracers between sub-cycles. */
        refresh_ghost_rows(config, temperature, nlev);
        refresh_ghost_rows(config, salinity, nlev);
    }

    /* 5. Recompute density, then remove static instabilities. */
    for (row = 0; row < DEEP_OCEAN_EXT_NLAT(config->nlat_local); ++row)
    {
        for (i = 0; i < nlon; ++i)
        {
            const int32_t horizontal = row * nlon + i;

            for (level = 0; level < nlev; ++level)
            {
                const size_t index = field_index(horizontal, level, nhor);

                if (is_land(land_mask, horizontal))
                {
                    continue;
                }
                density[index] = seawater_density(
                    config, temperature[index], salinity[index]);
            }
        }
    }
    apply_convective_adjustment(
        config, land_mask, layer_thickness, temperature, salinity, density);

    /* 6. Report the heat flux exchanged with the slab mixed layer:
     *    Q_deep = rho_o c_p,o (T^(0)_new - T_o^n) h_0 / dt, clamped to
     *    +-2000 W m^-2, where T^(0)_new is the top layer after steps 2-5 and
     *    T_o^n = slab_temperature the slab value it started from.  Positive
     *    Q_deep warms the slab (tutorial figures 4.3 / 4.5). */
    for (row = 0; row < config->nlat_local; ++row)
    {
        const int32_t extended_row = row + 1;
        const int32_t local_row = row * nlon;
        const int32_t ext_row = extended_row * nlon;

        for (i = 0; i < nlon; ++i)
        {
            const int32_t horizontal = ext_row + i;
            const int32_t local = local_row + i;
            const float h0 = fmaxf(1.0f, layer_thickness[0]);
            float flux;

            if (is_land(land_mask, horizontal))
            {
                deep_ocean_flux[local] = 0.0f;
                continue;
            }
            flux = config->rho0 * config->cp *
                   (temperature[field_index(horizontal, 0, nhor)] -
                    slab_temperature[local]) *
                   h0 / config->dt;
            deep_ocean_flux[local] = clamp_float(flux, -2000.0f, 2000.0f);
        }
    }

    free(pressure);
    free(tendency_t);
    free(tendency_s);
}
