#include "physics_internal.h"

static int initialize_misc(void)
{
    misc.moisture_fixer = 1;
    misc.nudging_time = 10.0f;
    misc_initialize(
        PLASIC_NLON, PLASIC_NLPP, plasic_gwd,
        6.28318530718, plasic_ww, &misc.nudging_time,
        misc.moisture_weight);
    return PLASIC_RUNTIME_OK;
}

static void initialize_flux(void)
{
    flux.vertical_diffusion = 1;
    flux.sensible_heat = 1;
    flux.evaporation = 1;
    flux.stress = 1;
    flux.surface_temperature_mode = 2;
    flux.minimum_wind = 1.0f;
    flux.mixing_length = 160.0f;
    flux.stability_b = 5.0f;
    flux.stability_c = 5.0f;
    flux.stability_d = 5.0f;
}

static void initialize_rain(void)
{
    rain.beta_mode = 1;
    rain.large_scale_rain = 1;
    rain.convective_rain = 1;
    rain.clouds = 1;
    rain.dry_adjustment = 1;
    rain.surface_convection = 1;
    rain.mix_momentum = 0;
    rain.shallow_convection = 1;
    rain.evaporate_precipitation = 1;
    rain.beta_exponent = 3;
    rain.critical_relative_humidity = 0.0f;
    rain.prescribed_beta = 0.0f;
    rain.cloud_vertical_velocity_1 = -0.1f;
    rain.cloud_vertical_velocity_2 = 0.0f;
    rain.deep_pressure = 999999.0f;
    rain.shallow_top_pressure = 70000.0f;
    rain.shallow_diffusion = 10.0f;
    rain.evaporation_rate = 0.01f;
    rain_initialize_defaults(
        PLASIC_NTRU, PLASIC_NLEV, plasic_sigma,
        &rain.shallow_convection,
        rain.critical_humidity);
    rain.cloud_factor = rain_cloud_factor(
        rain.cloud_vertical_velocity_1, rain.cloud_vertical_velocity_2);
}

static void initialize_radiation(void)
{
    /*
     * The standard PlaSiC calendar begins in AD 1.  The original Berger
     * Earth-orbit series evaluates to these parameters for that year.
     * Keeping the evaluated Earth constants here removes the large Fortran
     * coefficient tables while preserving the standard experiment.
     */
    const float orbital_obliquity_degrees = 23.695250f;
    const float orbital_perihelion_degrees = 68.836624f;
    const float eccentricity_squared = 0.017465f * 0.017465f;
    const float eccentricity_cubed =
        eccentricity_squared * 0.017465f;
    const float beta = sqrtf(1.0f - eccentricity_squared);

    radiation.use_synthetic_ozone = 1;
    radiation.external_ozone_column_reported = 0;
    radiation.shortwave = 1;
    radiation.longwave = 1;
    radiation.rayleigh_scattering = 1;
    radiation.solar_constant = 1365.0f;
    radiation.water_continuum = 0.024f;
    radiation.cloud_tuning_1 = 0.077f;
    radiation.cloud_tuning_2 = 0.065f;
    radiation.cloud_tuning_3 = 0.0055f;
#if PLASIC_NTRU == 42 && PLASIC_NLEV == 10
    radiation.water_continuum = 0.0285f;
#endif
    radiation.cloud_absorption = 0.100f;
    radiation.ozone_a0 = 0.25f;
    radiation.ozone_a1 = 0.11f;
    radiation.ozone_seasonal_amplitude = 0.08f;
    radiation.ozone_height = 20000.0f;
    radiation.ozone_thickness = 5000.0f;
    radiation.ozone_seasonal_offset = 90.0f;

    radiation.orbital_eccentricity = 0.017465f;
    radiation.obliquity_radians =
        orbital_obliquity_degrees * 3.14159265359f / 180.0f;
    radiation.perihelion_plus_pi =
        (orbital_perihelion_degrees + 180.0f) *
        3.14159265359f / 180.0f;
    radiation.mean_longitude_equinox =
        2.0f *
        ((0.5f * radiation.orbital_eccentricity +
          0.125f * eccentricity_cubed) *
             (1.0f + beta) * sinf(radiation.perihelion_plus_pi) -
         0.25f * eccentricity_squared * (0.5f + beta) *
             sinf(2.0f * radiation.perihelion_plus_pi) +
         0.125f * eccentricity_cubed * (1.0f / 3.0f + beta) *
             sinf(3.0f * radiation.perihelion_plus_pi));
    memset(
        radiation.carbon_dioxide, 0,
        sizeof(radiation.carbon_dioxide));
}

/*
 * Read one deep-ocean restart record.  The file stores the contiguous local
 * latitude slice, while the in-memory state is halo-extended, so the record is
 * unpacked after a successful read.  A missing record reports failure and
 * leaves the cold-start state untouched.
 */
static int read_deep_restart_field(const char *name, float *extended)
{
    float *local = (float *)malloc(
        (size_t)PLASIC_NHOR * (size_t)PLASIC_DEEP_NLEV * sizeof(float));
    int ok = 0;

    if (local != NULL &&
        restart_read_array_slice(
            name, local, PLASIC_NUGP,
            mp_local_offset(PLASIC_NHOR), PLASIC_NHOR,
            PLASIC_DEEP_NLEV) == PLASIC_RESTART_OK)
    {
        physics_ocean3d_unpack_local(local, extended, PLASIC_DEEP_NLEV);
        ok = 1;
    }
    free(local);
    return ok;
}

static int load_physics_restart(const char *path)
{
    int status;

    if (path == NULL)
    {
        /*
         * Cold start: no restart state is loaded; the surface-derived
         * initialization below fills every prognostic field.
         */
        return PLASIC_RUNTIME_OK;
    }
    status = restart_open_read(path);
    if (status != PLASIC_RESTART_OK)
    {
        return PLASIC_RUNTIME_RESTART_FAILED;
    }

#define LOAD_INTEGER(name, destination)                                         \
    do                                                                          \
    {                                                                           \
        status = restart_read_integer((name), &(destination));                  \
        if (status != PLASIC_RESTART_OK)                                        \
        {                                                                       \
            (void)restart_close();                                     \
            return PLASIC_RUNTIME_RESTART_FAILED;                               \
        }                                                                       \
    } while (0)
#define LOAD_ARRAY(name, values, columns)                                       \
    do                                                                          \
    {                                                                           \
        status = restart_read_array_slice(                                      \
            (name), (values), PLASIC_NUGP,                                      \
            mp_local_offset(PLASIC_NHOR), PLASIC_NHOR, (columns));              \
        if (status != PLASIC_RESTART_OK)                                        \
        {                                                                       \
            (void)restart_close();                                     \
            return PLASIC_RUNTIME_RESTART_FAILED;                               \
        }                                                                       \
    } while (0)

    /* LANDMOD persistent state. */
    {
        int32_t restart_soil_levels;
        LOAD_INTEGER("nlsoil", restart_soil_levels);
        if (restart_soil_levels != PLASIC_LAND_SOIL_LEVELS)
        {
            (void)restart_close();
            return PLASIC_RUNTIME_RESTART_FAILED;
        }
    }
    LOAD_ARRAY("dtsl", land.surface_temperature, 1);
    LOAD_ARRAY("dtsm", land.previous_surface_temperature, 1);
    LOAD_ARRAY("dqs", land.surface_humidity, 1);
    LOAD_ARRAY("dwmax", plasic_dwmax, 1);
    LOAD_ARRAY("dtcl", land.temperature_climatology,
               PLASIC_LAND_CLIMATOLOGY_MONTHS);
    LOAD_ARRAY("dwcl", land.wetness_climatology,
               PLASIC_LAND_CLIMATOLOGY_MONTHS);
    LOAD_ARRAY("dsnowt", land.snow_temperature, 1);
    LOAD_ARRAY("dsnowz", land.snow_depth, 1);
    LOAD_ARRAY("dsoilt", land.soil_temperature,
               PLASIC_LAND_SOIL_LEVELS);
    LOAD_ARRAY("dglac", plasic_dglac, 1);
    LOAD_ARRAY("dz0clim", land.roughness_climatology, 1);
    LOAD_ARRAY("dalbcl", land.albedo_climatology,
               PLASIC_LAND_CLIMATOLOGY_MONTHS);

    /* SEAMOD accumulators.  dts is restored after sea initialization. */
    LOAD_INTEGER("naccua", sea.accumulation_count);
    LOAD_ARRAY("dts", sea.surface_temperature, 1);
    LOAD_ARRAY("cheata", sea.heat_accumulator, 1);
    LOAD_ARRAY("cpmea", sea.freshwater_accumulator, 1);
    LOAD_ARRAY("cprsa", sea.snowfall_accumulator, 1);
    LOAD_ARRAY("ctauxa", sea.stress_x_accumulator, 1);
    LOAD_ARRAY("ctauya", sea.stress_y_accumulator, 1);
    LOAD_ARRAY("cust3a", sea.friction_accumulator, 1);
    LOAD_ARRAY("cshfla", sea.sensible_heat_accumulator, 1);
    LOAD_ARRAY("cshdta", sea.sensible_derivative_accumulator, 1);
    LOAD_ARRAY("clhfla", sea.latent_heat_accumulator, 1);
    LOAD_ARRAY("clhdta", sea.latent_derivative_accumulator, 1);
    LOAD_ARRAY("cswfla", sea.shortwave_accumulator, 1);
    LOAD_ARRAY("clwfla", sea.longwave_accumulator, 1);

    /* ICEMOD prognostic, climatological, coupling and diagnostic state. */
    LOAD_INTEGER("naccuice", ice.output_accumulation_count);
    LOAD_INTEGER("naccuo", ice.ocean_accumulation_count);
    LOAD_INTEGER("nicec2d", ice.diagnosed_climatological_thickness);
    LOAD_ARRAY("xls", ice.land_mask, 1);
    LOAD_ARRAY("xts", ice.surface_temperature, 1);
    LOAD_ARRAY("xicec", ice.ice_cover, 1);
    LOAD_ARRAY("xiced", ice.ice_thickness, 1);
    LOAD_ARRAY("xsnow", ice.snow, 1);
    LOAD_ARRAY("xclicec", ice.climatological_ice_cover,
               PLASIC_CLIMATOLOGY_MONTHS);
    LOAD_ARRAY("xcliced", ice.climatological_ice_thickness,
               PLASIC_CLIMATOLOGY_MONTHS);
    LOAD_ARRAY("xclsst", ice.climatological_sst,
               PLASIC_CLIMATOLOGY_MONTHS);
    LOAD_ARRAY("cheat", ice.ocean_heat_accumulator, 1);
    LOAD_ARRAY("cpme", ice.ocean_freshwater_accumulator, 1);
    LOAD_ARRAY("ctaux", ice.ocean_stress_x_accumulator, 1);
    LOAD_ARRAY("ctauy", ice.ocean_stress_y_accumulator, 1);
    LOAD_ARRAY("cust3", ice.ocean_friction_accumulator, 1);
    LOAD_ARRAY("csnow", ice.ocean_snow_accumulator, 1);
    LOAD_ARRAY("xflxicea", ice.flux_correction_accumulator, 1);
    LOAD_ARRAY("xheata", ice.atmospheric_heat_accumulator, 1);
    LOAD_ARRAY("xcfluxr", ice.maximum_thickness_flux, 1);
    LOAD_ARRAY("xofluxa", ice.ocean_to_ice_flux_accumulator, 1);
    LOAD_ARRAY("xqmelta", ice.melt_energy_accumulator, 1);
    LOAD_ARRAY("xcfluxa", ice.ice_to_ocean_flux_accumulator, 1);
    LOAD_ARRAY("xcfluxra", ice.maximum_thickness_flux_accumulator, 1);
    LOAD_ARRAY("xcfluxna", ice.excess_melt_flux_accumulator, 1);
    LOAD_ARRAY("xsmelta", ice.snow_melt_accumulator, 1);
    LOAD_ARRAY("ximelta", ice.ice_energy_flux_accumulator, 1);
    LOAD_ARRAY("xtsfluxa", ice.surface_storage_accumulator, 1);
    LOAD_ARRAY("xfluxca", ice.diagnosed_conductive_accumulator, 1);
    LOAD_ARRAY("xscflxa", ice.snow_to_ice_flux_accumulator, 1);
    LOAD_ARRAY("xcpmea", ice.freshwater_diagnostic_accumulator, 1);
    LOAD_ARRAY("xstoia", ice.snow_to_ice_rate_accumulator, 1);
    LOAD_ARRAY("xicecc", ice.compactness, 1);

    /*
     * OCEANMOD state.  The mixed-layer ocean is forced to a single layer and
     * has no vertical or horizontal diffusion, so a restart that carries a
     * different layer count is rejected.
     */
    {
        int32_t restart_ocean_levels;
        LOAD_INTEGER("nlev_oce", restart_ocean_levels);
        if (restart_ocean_levels != 1)
        {
            (void)restart_close();
            return PLASIC_RUNTIME_RESTART_FAILED;
        }
    }
    LOAD_INTEGER("naccuoce", ocean.accumulation_count);
    LOAD_ARRAY("yls", ocean.land_mask, 1);
    LOAD_ARRAY("ysst", ocean.temperature, 1);
    LOAD_ARRAY("yiflux", ocean.ice_flux, 1);
    LOAD_ARRAY("yclsst", ocean.climatological_sst,
               PLASIC_CLIMATOLOGY_MONTHS);
    LOAD_ARRAY("yheata", ocean.ice_to_ocean_flux_accumulator, 1);
    LOAD_ARRAY("yifluxa", ocean.ice_flux_accumulator, 1);
    LOAD_ARRAY("yfldoa", ocean.deep_ocean_flux_accumulator, 1);
    LOAD_ARRAY("yfldo", ocean.deep_ocean_flux, 1);
    /*
     * yfldo 是三维深海洋最后一步交给 slab 的热通量。slab 构建没有深海洋来
     * 更新它；从三维海洋的重启文件启动时，旧值会冻结成一个永久热源（极地
     * 个别点可达 2000 W/m²），所以 slab 构建把它清零。
     *
     * yfldo is the heat flux the three-dimensional deep ocean handed to the
     * slab on its last step.  A slab build has no deep ocean to update it, so
     * a value restored from a three-dimensional-ocean restart would stay
     * frozen as a permanent heat source (up to 2000 W/m2 at single polar
     * points); a slab build therefore starts from zero.
     */
    if (deep_ocean.enabled == 0)
    {
        memset(ocean.deep_ocean_flux, 0, sizeof(ocean.deep_ocean_flux));
    }

    /*
     * Three-dimensional deep ocean.  Older restart files do not contain these
     * records, so a missing record is not an error: the deep ocean is then
     * re-seeded from the mixed-layer temperature during initialization.
     */
    if (deep_ocean.enabled != 0)
    {
        int32_t deep_levels = 0;
        int deep_ok = 1;

        if (restart_read_integer(
                "deep_nlev", &deep_levels) !=
                PLASIC_RESTART_OK ||
            deep_levels != PLASIC_DEEP_NLEV)
        {
            deep_ok = 0;
        }
        if (deep_ok)
        {
            deep_ok = read_deep_restart_field(
                "deept", deep_ocean.temperature);
        }
        if (deep_ok)
        {
            deep_ok = read_deep_restart_field(
                "deeps", deep_ocean.salinity);
        }
        if (deep_ok)
        {
            deep_ok = read_deep_restart_field(
                "deepu", deep_ocean.u_velocity);
        }
        if (deep_ok)
        {
            deep_ok = read_deep_restart_field(
                "deepv", deep_ocean.v_velocity);
        }
        if (deep_ok)
        {
            deep_ok = read_deep_restart_field(
                "deepw", deep_ocean.w_velocity);
        }
        deep_ocean.restart_loaded = deep_ok;
    }

#undef LOAD_ARRAY
#undef LOAD_INTEGER

    status = restart_close();
    return status == PLASIC_RESTART_OK ?
               PLASIC_RUNTIME_OK : PLASIC_RUNTIME_RESTART_FAILED;
}

static int read_surface_field(
    const char *directory, const char *field_name, int32_t months, float *values,
    int optional)
{
    const int status = surface_data_read_local(
        directory, PLASIC_NTRU, PLASIC_NLAT, PLASIC_NLON, field_name, months,
        mp_local_offset(PLASIC_NHOR), PLASIC_NHOR, values);

    if (status == SURFACE_DATA_OK ||
        (optional && status == SURFACE_DATA_NOT_FOUND))
    {
        return PLASIC_RUNTIME_OK;
    }
    physics_set_error(
        "cold start cannot read NetCDF surface field '%s' from '%s' "
        "(status %d)",
        field_name, directory, status);
    return PLASIC_RUNTIME_IO_FAILED;
}

static int load_surface_initialization(const char *directory)
{
    int status;
    int32_t point;

#define READ_SURFACE(field_name, months, values)                              \
    do                                                                        \
    {                                                                         \
        status = read_surface_field(                                          \
            directory, (field_name), (months), (values), 0);                  \
        if (status != PLASIC_RUNTIME_OK)                                      \
        {                                                                     \
            return status;                                                    \
        }                                                                     \
    } while (0)

    /* Orography is intentionally not loaded here: the atmospheric cold
     * start reads it exactly once, transforms it to spectral space, and
     * shares it through state. */
    READ_SURFACE("land_sea_mask", 1, plasic_dls);
    for (point = 0; point < PLASIC_NHOR; ++point)
    {
        /*
         * LANDINI, ICEINI and OCEANINI all convert the fractional source
         * mask to a binary model mask.  This is essential because restart
         * land temperatures use 1.e20 sentinels at sea points.
         */
        plasic_dls[point] =
            plasic_dls[point] > 0.5f ? 1.0f : 0.0f;
    }
    memcpy(ice.land_mask, plasic_dls, sizeof(ice.land_mask));
    memcpy(ocean.land_mask, plasic_dls, sizeof(ocean.land_mask));

    READ_SURFACE(
        "surface_temperature", PLASIC_CLIMATOLOGY_MONTHS,
        land.temperature_climatology);
    memcpy(
        ice.climatological_sst, land.temperature_climatology,
        sizeof(ice.climatological_sst));
    memcpy(
        ocean.climatological_sst, land.temperature_climatology,
        sizeof(ocean.climatological_sst));
    READ_SURFACE(
        "sea_ice_cover", PLASIC_CLIMATOLOGY_MONTHS,
        ice.climatological_ice_cover);
    {
        float maximum_ice_cover = 0.0f;
        size_t index;

        for (index = 0;
             index < sizeof(ice.climatological_ice_cover) /
                         sizeof(ice.climatological_ice_cover[0]);
             ++index)
        {
            maximum_ice_cover = fmaxf(
                maximum_ice_cover,
                ice.climatological_ice_cover[index]);
        }
        if (maximum_ice_cover > 5.0f)
        {
            for (index = 0;
                 index < sizeof(ice.climatological_ice_cover) /
                             sizeof(ice.climatological_ice_cover[0]);
                 ++index)
            {
                ice.climatological_ice_cover[index] *= 0.01f;
            }
        }
        for (index = 0;
             index < sizeof(ocean.climatological_sst) /
                         sizeof(ocean.climatological_sst[0]);
             ++index)
        {
            ocean.climatological_sst[index] =
                fmaxf(ocean.climatological_sst[index], 271.25f);
            /*
             * The source climatology carries small interpolation over-
             * and undershoots; a physical cover fraction stays in [0,1]
             * (negative covers would produce negative sea roughness).
             */
            ice.climatological_ice_cover[index] = fminf(
                1.0f, fmaxf(0.0f, ice.climatological_ice_cover[index]));
            if (plasic_dls[index % PLASIC_NHOR] >= 1.0f)
            {
                ice.climatological_ice_cover[index] = 0.0f;
            }
        }
    }
    READ_SURFACE("roughness_length", 1, land.roughness_climatology);
    READ_SURFACE(
        "background_albedo", PLASIC_CLIMATOLOGY_MONTHS,
        land.albedo_climatology);
    READ_SURFACE("glacier_fraction", 1, plasic_dglac);
    for (point = 0; point < PLASIC_NHOR; ++point)
    {
        plasic_dglac[point] =
            plasic_dglac[point] > 0.5f ? 1.0f : 0.0f;
    }

    READ_SURFACE("soil_water_capacity", 1, plasic_dwmax);

#undef READ_SURFACE
    return PLASIC_RUNTIME_OK;
}

static void initialize_land_defaults(void)
{
    land.prognostic_temperature = 1;
    land.prognostic_water = 1;
    land.land_albedo = 0.2f;
    land.snow_albedo_minimum = 0.4f;
    land.snow_albedo_maximum = 0.8f;
    land.glacier_albedo_minimum = 0.6f;
    land.glacier_albedo_maximum = 0.8f;
    land.land_roughness = 2.0f;
    land.land_wetness = 0.25f;
    land.full_wetness_fraction = 0.4f;
    land.glacier_height = -1.0f;
    land.top_layer_depth = 0.20f;
    land.maximum_snow_depth = 5.0f;
    land.maximum_soil_water = 0.5f;
    land.snow_density = 330.0f;
    land.soil_conductivity = 1.8f;
    land.ice_conductivity = 2.03f;
    land.snow_conductivity = 0.31f;
    land.soil_heat_capacity = 2.4e6f;
    land.ice_heat_capacity = 2.07e6f;
    land.snow_heat_capacity = 0.6897e6f;
    land.soil_depth[0] = 0.4f;
    land.soil_depth[1] = 0.8f;
    land.soil_depth[2] = 1.6f;
    land.soil_depth[3] = 3.2f;
    land.soil_depth[4] = 6.4f;
}

static void initialize_deep_ocean_defaults(void)
{
    deep_ocean.enabled = PLASIC_DEEP_OCEAN;
    deep_ocean.initialized = 0;
    deep_ocean.restart_loaded = 0;
    deep_ocean.nlev = PLASIC_DEEP_NLEV;
    deep_ocean.timestep = 0.0f;
    deep_ocean.radius = plasic_plarad;
    deep_ocean.solar_day = plasic_solar_day;
    deep_ocean.sidereal_day = plasic_sidereal_day;
    deep_ocean.gravity = plasic_ga;
    deep_ocean.cp = 4180.0f;
    deep_ocean.rho0 = 1030.0f;
    deep_ocean.alpha = 2.0e-4f;
    deep_ocean.beta = 7.5e-4f;
    deep_ocean.kappa_v = 2.0e-4f;
    deep_ocean.kappa_h = 500.0f;
    deep_ocean.t_freeze = 271.25f;
    deep_ocean_build_levels(
        PLASIC_DEEP_NLEV, deep_ocean.layer_thickness,
        deep_ocean.layer_center_depth);
}

static void initialize_sea_ice_ocean_defaults(void)
{
    sea.coupling_interval = 1;
    sea.sea_albedo = 0.069f;
    sea.ice_albedo = 0.7f;
    sea.sea_roughness = 1.5e-5f;
    sea.ice_roughness = 0.001f;
    sea.sea_relative_humidity = 1.0f;
    sea.ice_relative_humidity = 1.0f;
    sea.charnock = 0.018f;

    ice.prognostic_ice = 1;
    ice.prognostic_snow = 1;
    ice.prognostic_skin_temperature = 1;
    ice.ocean_coupling_interval = 1;
    ice.output_interval = 32;
    ice.maximum_thickness_correction = 1;
    ice.minimum_thickness = 0.1f;
    ice.maximum_thickness = 9.0f;
    ice.compactness_threshold = 0.5f;
    ice.minimum_ice_cover = 0.5f;
    /*
     * One atmospheric step.  sea_step() resets it to the length of the
     * coupling window before every ice_step().
     */
    ice.timestep = plasic_solar_day / (float)plasic_ntspd;

    ocean.output_interval = 32;
    ocean.prognostic_ocean = 1;
    ocean.flux_correction = 0;
    ocean.restoring_timescale = 0.0f;
    ocean.layer_depth = 50.0f;

    initialize_deep_ocean_defaults();
}

void physics_interpolate_ice_climatology(void)
{
    int32_t previous_month1;
    int32_t previous_month2;
    int32_t current_month1;
    int32_t current_month2;
    int32_t previous_step = plasic_nstep;
    int32_t current_step = plasic_nstep + 1;
    float previous_weight2;
    float current_weight2;

    momint(plasic_calendar, previous_step, &previous_month1,
            &previous_month2, &previous_weight2);
    momint(plasic_calendar, current_step, &current_month1,
            &current_month2, &current_weight2);
    ice_interpolate_climatology(
        PLASIC_NHOR, previous_month1, previous_month2,
        1.0f - previous_weight2, previous_weight2,
        current_month1, current_month2, 1.0f - current_weight2,
        current_weight2, ice.climatological_sst,
        ice.climatological_ice_cover,
        ice.climatological_ice_thickness,
        ice.previous_climatological_sst,
        ice.current_climatological_sst,
        ice.current_climatological_cover,
        ice.current_climatological_thickness);
}

void physics_interpolate_ocean_climatology(void)
{
    int32_t first_month;
    int32_t second_month;
    int32_t next_step = plasic_nstep + 1;
    float weight;

    momint(plasic_calendar, next_step, &first_month, &second_month,
            &weight);
    ocean_interpolate_cycle(
        PLASIC_NHOR, first_month, second_month, weight,
        ocean.climatological_sst, ocean.interpolated_sst);
}

/*
 * Cold start: seed every prognostic field from the packaged surface
 * climatologies for the authoritative model start month. Runs after
 * load_surface_initialization so the
 * climatological fields are available, and before initialize_sea_ice_ocean
 * so the coupler sees a consistent land/ice/ocean state.
 */
void physics_interpolate_land_surface(void);

static int initialize_cold_start_physics(void)
{
    const size_t plane = (size_t)PLASIC_NUGP;
    float *global_cover = NULL;
    float *global_thickness = NULL;
    int status = PLASIC_RUNTIME_OK;
    int32_t point;
    int32_t month;
    int32_t level;
    const int32_t initial_month = plasic_ndatim[1];

    if (initial_month < 1 || initial_month > 12)
    {
        physics_set_error(
            "cold-start calendar produced invalid month %d", initial_month);
        return PLASIC_RUNTIME_INITIALIZATION_FAILED;
    }

    for (point = 0; point < PLASIC_NHOR; ++point)
    {
        const int is_land = plasic_dls[point] >= 1.0f;

        /*
         * Restart land temperatures use 1.e20 sentinels at sea points;
         * reproduce that convention for the cold start.
         */
        land.surface_temperature[point] =
            is_land ? land.temperature_climatology[
                          (size_t)initial_month * PLASIC_NHOR +
                          (size_t)point] :
                      1.0e20f;
        land.previous_surface_temperature[point] =
            land.surface_temperature[point];
        land.snow_temperature[point] = land.surface_temperature[point];
        land.surface_humidity[point] = is_land ? 0.0f : 1.0e20f;
        plasic_dwatc[point] =
            is_land ? 0.25f * plasic_dwmax[point] : 0.0f;
        land.snow_depth[point] = 0.0f;
        for (level = 0; level < PLASIC_LAND_SOIL_LEVELS; ++level)
        {
            land.soil_temperature[
                (size_t)point + (size_t)level * PLASIC_NHOR] =
                land.surface_temperature[point];
        }

        /* The 14-slot layout is [Dec halo, Jan, ..., Dec, Jan halo]. */
        ocean.temperature[point] = ocean.climatological_sst[
            (size_t)initial_month * PLASIC_NHOR + (size_t)point];
        ice.surface_temperature[point] = ice.climatological_sst[
            (size_t)initial_month * PLASIC_NHOR + (size_t)point];
        ice.ice_cover[point] = ice.climatological_ice_cover[
            (size_t)initial_month * PLASIC_NHOR + (size_t)point];
        /*
         * Seed the prognostic compactness c_i^* once from the initial
         * cover; ice_step() carries it afterwards (restart record xicecc).
         */
        ice.compactness[point] = ice.ice_cover[point];
        land.interpolated_albedo[point] = land.albedo_climatology[
            (size_t)initial_month * PLASIC_NHOR + (size_t)point];
        ice.snow[point] = 0.0f;
    }
    for (month = 0; month < PLASIC_LAND_CLIMATOLOGY_MONTHS; ++month)
    {
        for (point = 0; point < PLASIC_NHOR; ++point)
        {
            land.wetness_climatology[
                (size_t)point + (size_t)month * PLASIC_NHOR] = 0.25f;
        }
    }

    /*
     * LANDINI derives the exchange fields (roughness, albedo, wetness,
     * snow cover and the surface-layer temperature/humidity seen by the
     * atmosphere) from the freshly seeded state.  The restart path reads
     * these grids from the restart file instead.  The gridpoint pressure
     * is not available yet, so the exchange uses the constant reference
     * pressure.  Sea points are overwritten by sea_initialize later.
     */
    land_surface_exchange(
        PLASIC_NHOR, PLASIC_LAND_SOIL_LEVELS, 1, 1,
        plasic_psurf, plasic_rdbrv, plasic_ra1,
        plasic_ra2, plasic_ra4, plasic_tmelt,
        land.snow_albedo_minimum, land.snow_albedo_maximum,
        land.full_wetness_fraction, plasic_dp, plasic_dls,
        land.surface_temperature, land.snow_depth, plasic_dwatc,
        plasic_dwmax, land.interpolated_albedo,
        land.soil_temperature, land.roughness_climatology,
        land.previous_surface_temperature, land.surface_humidity,
        plasic_dsnow, plasic_dalb, plasic_dwetfac,
        plasic_dz0, plasic_dtsoil, plasic_dtd2,
        plasic_dtd3, plasic_dtd4, plasic_dtd5,
        plasic_dt + (size_t)PLASIC_NHOR * PLASIC_NLEV,
        plasic_dq + (size_t)PLASIC_NHOR * PLASIC_NLEV);
    land_apply_glacier_surface(
        PLASIC_NHOR, 1, land.maximum_snow_depth,
        land.glacier_albedo_minimum, land.glacier_albedo_maximum,
        plasic_tmelt, plasic_dls, plasic_dglac,
        land.surface_temperature, land.snow_depth, plasic_dsnow,
        plasic_dalb, plasic_dwetfac);

    /*
     * Diagnose the climatological ice thickness from the 12-month cover.
     * ice_make_thickness expects the global (longitude, latitude, month)
     * grid, so gather every monthly field, compute on the full grid, and
     * copy the local slice back.  Accumulators start at zero for a fresh
     * integration.
     */
    global_cover =
        (float *)malloc(plane * PLASIC_CLIMATOLOGY_MONTHS * sizeof(float));
    global_thickness =
        (float *)malloc(plane * PLASIC_CLIMATOLOGY_MONTHS * sizeof(float));
    if (global_cover == NULL || global_thickness == NULL)
    {
        status = PLASIC_RUNTIME_INITIALIZATION_FAILED;
        goto done;
    }
    for (month = 0; month < PLASIC_CLIMATOLOGY_MONTHS; ++month)
    {
        status = mp_allgather_grid(
            global_cover + (size_t)month * plane,
            ice.climatological_ice_cover + (size_t)month * PLASIC_NHOR,
            PLASIC_NUGP, PLASIC_NHOR, 1) == 0
                     ? PLASIC_RUNTIME_OK
                     : PLASIC_RUNTIME_NUMERICAL_FAILURE;
        if (status != PLASIC_RUNTIME_OK)
        {
            goto done;
        }
    }
    ice_make_thickness(
        PLASIC_NLON, PLASIC_NLAT, PLASIC_CLIMATOLOGY_MONTHS,
        global_cover, global_thickness);
    for (month = 0; month < PLASIC_CLIMATOLOGY_MONTHS; ++month)
    {
        memcpy(
            ice.climatological_ice_thickness +
                (size_t)month * PLASIC_NHOR,
            global_thickness + (size_t)month * plane +
                mp_local_offset(PLASIC_NHOR),
            (size_t)PLASIC_NHOR * sizeof(float));
    }
    memcpy(
        ice.ice_thickness,
        ice.climatological_ice_thickness +
            (size_t)initial_month * PLASIC_NHOR,
        sizeof(ice.ice_thickness));
    ice.diagnosed_climatological_thickness = 1;

    memset(sea.heat_accumulator, 0, sizeof(sea.heat_accumulator));
    memset(
        sea.freshwater_accumulator, 0, sizeof(sea.freshwater_accumulator));
    memset(sea.snowfall_accumulator, 0, sizeof(sea.snowfall_accumulator));
    memset(sea.stress_x_accumulator, 0, sizeof(sea.stress_x_accumulator));
    memset(sea.stress_y_accumulator, 0, sizeof(sea.stress_y_accumulator));
    memset(sea.friction_accumulator, 0, sizeof(sea.friction_accumulator));
    memset(
        sea.sensible_heat_accumulator, 0,
        sizeof(sea.sensible_heat_accumulator));
    memset(
        sea.sensible_derivative_accumulator, 0,
        sizeof(sea.sensible_derivative_accumulator));
    memset(
        sea.latent_heat_accumulator, 0, sizeof(sea.latent_heat_accumulator));
    memset(
        sea.latent_derivative_accumulator, 0,
        sizeof(sea.latent_derivative_accumulator));
    memset(
        sea.shortwave_accumulator, 0, sizeof(sea.shortwave_accumulator));
    memset(
        sea.longwave_accumulator, 0, sizeof(sea.longwave_accumulator));
    sea.accumulation_count = 0;

    memset(
        ice.ocean_heat_accumulator, 0, sizeof(ice.ocean_heat_accumulator));
    memset(
        ice.ocean_freshwater_accumulator, 0,
        sizeof(ice.ocean_freshwater_accumulator));
    memset(
        ice.ocean_stress_x_accumulator, 0,
        sizeof(ice.ocean_stress_x_accumulator));
    memset(
        ice.ocean_stress_y_accumulator, 0,
        sizeof(ice.ocean_stress_y_accumulator));
    memset(
        ice.ocean_friction_accumulator, 0,
        sizeof(ice.ocean_friction_accumulator));
    memset(
        ice.ocean_snow_accumulator, 0, sizeof(ice.ocean_snow_accumulator));
    memset(
        ice.flux_correction_accumulator, 0,
        sizeof(ice.flux_correction_accumulator));
    memset(
        ice.atmospheric_heat_accumulator, 0,
        sizeof(ice.atmospheric_heat_accumulator));
    memset(
        ice.ocean_to_ice_flux_accumulator, 0,
        sizeof(ice.ocean_to_ice_flux_accumulator));
    memset(
        ice.melt_energy_accumulator, 0,
        sizeof(ice.melt_energy_accumulator));
    memset(
        ice.ice_to_ocean_flux_accumulator, 0,
        sizeof(ice.ice_to_ocean_flux_accumulator));
    memset(
        ice.maximum_thickness_flux_accumulator, 0,
        sizeof(ice.maximum_thickness_flux_accumulator));
    memset(
        ice.excess_melt_flux_accumulator, 0,
        sizeof(ice.excess_melt_flux_accumulator));
    memset(
        ice.snow_melt_accumulator, 0, sizeof(ice.snow_melt_accumulator));
    memset(
        ice.ice_energy_flux_accumulator, 0, sizeof(ice.ice_energy_flux_accumulator));
    memset(
        ice.surface_storage_accumulator, 0,
        sizeof(ice.surface_storage_accumulator));
    memset(
        ice.diagnosed_conductive_accumulator, 0,
        sizeof(ice.diagnosed_conductive_accumulator));
    memset(
        ice.snow_to_ice_flux_accumulator, 0,
        sizeof(ice.snow_to_ice_flux_accumulator));
    memset(
        ice.freshwater_diagnostic_accumulator, 0,
        sizeof(ice.freshwater_diagnostic_accumulator));
    memset(
        ice.snow_to_ice_rate_accumulator, 0, sizeof(ice.snow_to_ice_rate_accumulator));
    ice.output_accumulation_count = 0;
    ice.ocean_accumulation_count = 0;

    memset(ocean.ice_to_ocean_flux_accumulator, 0, sizeof(ocean.ice_to_ocean_flux_accumulator));
    memset(
        ocean.ice_flux_accumulator, 0, sizeof(ocean.ice_flux_accumulator));
    memset(
        ocean.deep_ocean_flux_accumulator, 0,
        sizeof(ocean.deep_ocean_flux_accumulator));
    ocean.accumulation_count = 0;

done:
    free(global_cover);
    free(global_thickness);
    return status;
}

static void initialize_sea_ice_ocean(int from_restart)
{
    const int32_t surface_level = PLASIC_NLEP - 1;

    ocean_initialize_parameters(
        PLASIC_NHOR, plasic_ntspd, plasic_solar_day,
        &ocean.restoring_timescale, ocean.layer_depth,
        ocean.mixed_layer_depth, &ocean.timestep);
    physics_interpolate_ocean_climatology();
    ocean_copy_outputs(
        PLASIC_NHOR, ocean.temperature, ocean.mixed_layer_depth,
        ocean.ice_flux, ice.sst, ice.mixed_layer_depth,
        ice.ocean_to_ice_flux);
    physics_interpolate_ice_climatology();

    memcpy(sea.coupled_surface_temperature, ice.surface_temperature,
           sizeof(sea.coupled_surface_temperature));
    memcpy(sea.coupled_sst, ice.sst, sizeof(sea.coupled_sst));
    memcpy(sea.coupled_mixed_layer_depth, ice.mixed_layer_depth,
           sizeof(sea.coupled_mixed_layer_depth));
    memcpy(sea.coupled_ice_cover, ice.ice_cover,
           sizeof(sea.coupled_ice_cover));
    memcpy(sea.coupled_ice_thickness, ice.ice_thickness,
           sizeof(sea.coupled_ice_thickness));
    memcpy(sea.coupled_snow, ice.snow, sizeof(sea.coupled_snow));

    /*
     * SEAMOD initializes dts from ICEMOD, then the original restart path
     * restores dts.  Preserve that ordering with a temporary copy.  On a
     * cold start sea_initialize must also derive the surface humidity,
     * wetness, albedo and roughness of the sea points (nrestart == 0);
     * a restart keeps the values read from the restart file.
     */
    memcpy(work->horizontal_h, sea.surface_temperature,
           sizeof(sea.surface_temperature));
    sea_initialize(
        PLASIC_NHOR, surface_level, from_restart, plasic_rdbrv,
        plasic_ra1, plasic_ra2, plasic_ra4,
        plasic_tmelt, plasic_psurf, sea.sea_albedo,
        sea.ice_albedo, sea.sea_roughness, sea.ice_roughness,
        sea.sea_relative_humidity, sea.ice_relative_humidity,
        plasic_dls, sea.coupled_surface_temperature,
        sea.coupled_sst, sea.coupled_mixed_layer_depth,
        sea.coupled_ice_cover, sea.coupled_ice_thickness,
        sea.coupled_snow, sea.surface_temperature,
        plasic_dicec, plasic_diced, plasic_dsnow,
        sea.sst, plasic_dmld, plasic_dt, plasic_dq,
        sea.surface_humidity, plasic_dwetfac, plasic_dalb,
        plasic_dz0);
    memcpy(sea.surface_temperature, work->horizontal_h,
           sizeof(sea.surface_temperature));

    /*
     * Seed or restore the three-dimensional deep ocean.  On cold start it is
     * profiled from the mixed-layer temperature; on restart it is read back
     * from the deep records when they are present.
     */
    physics_ocean3d_initialize(from_restart);
}

int physics_initialize(
    const char *restart_path, const char *surface_data_directory,
    int32_t sea_coupling_interval)
{
    int status;

    physics_error_message[0] = '\0';
    work = (physics_workspace *)calloc(1, sizeof(*work));
    if (work == NULL)
    {
        return PLASIC_RUNTIME_INITIALIZATION_FAILED;
    }
    (void)initialize_misc();
    initialize_flux();
    initialize_rain();
    initialize_radiation();
    initialize_land_defaults();
    initialize_sea_ice_ocean_defaults();
    /* A run option, not a restart record: each run may choose its own. */
    sea.coupling_interval = sea_coupling_interval;
    status = load_physics_restart(restart_path);
    if (status != PLASIC_RUNTIME_OK)
    {
        return status;
    }
    if (restart_path == NULL)
    {
        status = load_surface_initialization(surface_data_directory);
        if (status != PLASIC_RUNTIME_OK)
        {
            return status;
        }
        status = initialize_cold_start_physics();
        if (status != PLASIC_RUNTIME_OK)
        {
            return status;
        }
    }
    /* On restart, the records loaded above are authoritative. In particular,
     * do not overwrite dls, climatologies, or any prognostic/coupling state
     * from surface files; --data-dir is deliberately irrelevant. */
    initialize_sea_ice_ocean(restart_path != NULL);
    return PLASIC_RUNTIME_OK;
}

void physics_finalize(void)
{
    free(work);
    work = NULL;
}
