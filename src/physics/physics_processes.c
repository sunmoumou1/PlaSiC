#include "physics_internal.h"

static void orbital_declination(
    float calendar_day, float days_per_year, float *declination,
    float *distance_factor)
{
    static const float vernal_equinox_day = 80.5f;
    const float mean_longitude =
        radiation.mean_longitude_equinox +
        (calendar_day - vernal_equinox_day) * 6.28318530718f /
            days_per_year;
    const float anomaly =
        mean_longitude - radiation.perihelion_plus_pi;
    const float sin_anomaly = sinf(anomaly);
    const float true_longitude =
        mean_longitude +
        radiation.orbital_eccentricity *
            (2.0f * sin_anomaly +
             radiation.orbital_eccentricity *
                 (1.25f * sinf(2.0f * anomaly) +
                  radiation.orbital_eccentricity *
                      ((13.0f / 12.0f) * sinf(3.0f * anomaly) -
                       0.25f * sin_anomaly)));
    const float inverse_distance =
        (1.0f +
         radiation.orbital_eccentricity *
             cosf(true_longitude - radiation.perihelion_plus_pi)) /
        (1.0f - radiation.orbital_eccentricity *
                    radiation.orbital_eccentricity);

    *declination = asinf(
        sinf(radiation.obliquity_radians) * sinf(true_longitude));
    *distance_factor = inverse_distance * inverse_distance;
}

static int radiation_step(void)
{
    float calendar_day;
    float declination;
    float distance_factor;
    int32_t day;
    int32_t minute_of_day;
    int32_t days_per_year;
    int status;

    rad_zero_fluxes(
        PLASIC_NHOR, PLASIC_NLEP, plasic_dfu, plasic_dfd,
        plasic_dftu, plasic_dftd, plasic_dswfl,
        plasic_dlwfl);
    day = ndayofyear(plasic_calendar, plasic_nstep);
    days_per_year = days_in_year(plasic_calendar, plasic_ndatim[0]);
    minute_of_day = plasic_ndatim[3] * 60 + plasic_ndatim[4];
    calendar_day =
        (float)day + (float)minute_of_day / 1440.0f;
    orbital_declination(
        calendar_day, (float)days_per_year, &declination,
        &distance_factor);
    rad_solar_angles(
        PLASIC_NLON, PLASIC_NLPP, minute_of_day,
        3.14159265359f, 6.28318530718f, declination,
        distance_factor, plasic_sid, plasic_cola,
        radiation.cosine_zenith,
        &radiation.sun_distance_factor);
    if (radiation.use_synthetic_ozone == 1)
    {
        status = rad_synthetic_ozone(
            PLASIC_NLON, PLASIC_NLPP, PLASIC_NLEV, calendar_day,
            days_per_year,
            6.28318530718f, radiation.ozone_a0, radiation.ozone_a1,
            radiation.ozone_seasonal_amplitude, radiation.ozone_height,
            radiation.ozone_thickness, radiation.ozone_seasonal_offset,
            plasic_gascon, plasic_ga,
            plasic_sid, plasic_sigmah, plasic_dsigma,
            plasic_dp, plasic_dt, radiation.ozone);
        if (status == 0)
        {
            return PLASIC_RUNTIME_NUMERICAL_FAILURE;
        }
    }
    if (radiation.shortwave != 0)
    {
        status = rad_shortwave(
            PLASIC_NHOR, PLASIC_NLEV,
            radiation.rayleigh_scattering, plasic_ga,
            radiation.solar_constant, radiation.sun_distance_factor,
            radiation.cloud_tuning_1, radiation.cloud_tuning_2,
            radiation.cloud_tuning_3,
            plasic_sigma, plasic_dsigma, plasic_dp,
            plasic_dalb, plasic_dls, plasic_dicec,
            plasic_dq, plasic_dql, plasic_dcc,
            radiation.ozone, radiation.cosine_zenith, plasic_dt,
            plasic_dfu, plasic_dfd, plasic_dswfl);
        if (status == 0)
        {
            return PLASIC_RUNTIME_NUMERICAL_FAILURE;
        }
    }
    if (radiation.longwave != 0)
    {
        status = rad_longwave(
            PLASIC_NHOR, PLASIC_NLEV, plasic_ga,
            radiation.water_continuum, radiation.cloud_absorption,
            plasic_sigma,
            plasic_sigmah, plasic_dsigma, plasic_dp,
            plasic_dls, plasic_dq, plasic_dt,
            plasic_dcc, plasic_dql, radiation.ozone,
            radiation.carbon_dioxide, plasic_dftu,
            plasic_dftd, plasic_dlwfl);
        if (status == 0)
        {
            return PLASIC_RUNTIME_NUMERICAL_FAILURE;
        }
    }
    rad_compute_tendencies(
        PLASIC_NHOR, PLASIC_NLEV, plasic_ga, plasic_cpd,
        plasic_cpv_cpd_minus1, plasic_dsigma, plasic_dp,
        plasic_dq, plasic_dswfl, plasic_dlwfl,
        plasic_dflux, plasic_dtdt,
        radiation.shortwave_tendency, radiation.longwave_tendency,
        work->level_a);
    return PLASIC_RUNTIME_OK;
}

static int moisture_fixer(void)
{
    float totals[2];

    if (misc.moisture_fixer == 0)
    {
        return PLASIC_RUNTIME_OK;
    }
    misc_prepare_fixer(
        PLASIC_NHOR, PLASIC_NLEV, plasic_dq, plasic_dsigma,
        plasic_dp, misc.moisture_weight, work->level_a,
        work->horizontal_a, work->horizontal_b);
    if (misc_sum(PLASIC_NHOR, work->horizontal_a) > 0.0f)
    {
        misc_fix_columns(
            PLASIC_NHOR, PLASIC_NLEV, work->level_a,
            work->horizontal_a, work->horizontal_b);
    }
    if (misc_sum(PLASIC_NHOR, work->horizontal_a) > 0.0f)
    {
        misc_fix_latitudes(
            PLASIC_NLON, PLASIC_NLPP, PLASIC_NLEV, work->level_a,
            work->horizontal_a, work->horizontal_b);
    }
    misc_local_totals(
        PLASIC_NHOR, work->horizontal_a, work->horizontal_b, totals);
    if (mp_allreduce_real(totals, 2) != 0)
    {
        return PLASIC_RUNTIME_NUMERICAL_FAILURE;
    }
    misc_finish_fixer(
        PLASIC_NHOR, PLASIC_NLEV, totals[0], totals[1], plasic_dp,
        plasic_delt2, plasic_psurf, work->level_a,
        plasic_dq, plasic_gqdt);
    return PLASIC_RUNTIME_OK;
}

static void surface_fluxes(void)
{
    const int32_t low_level = PLASIC_NLEV - 1;
    const int32_t surface_level = PLASIC_NLEP - 1;

    flux_surface_exchange_coefficients(
        PLASIC_NHOR, low_level, surface_level,
        flux.surface_temperature_mode, 0.4f, plasic_kap,
        plasic_sigma[PLASIC_NLEV - 1], plasic_rdbrv,
        plasic_gascon, plasic_ga, flux.minimum_wind,
        flux.stability_b, flux.stability_c, flux.stability_d,
        plasic_dt, plasic_dq, plasic_du, plasic_dv,
        plasic_dz0, plasic_dls, plasic_dtsa,
        flux.heat_transfer, flux.momentum_transfer, work->horizontal_a,
        work->horizontal_b, work->horizontal_c, work->horizontal_d,
        work->horizontal_e);
    if (flux.stress != 0)
    {
        flux_mkstress(
            PLASIC_NHOR, low_level, surface_level, plasic_ndheat,
            plasic_ga, plasic_deltsec2, plasic_gascon,
            plasic_dsigma[PLASIC_NLEV - 1], plasic_cpd,
            plasic_cpv_cpd_minus1, flux.momentum_transfer, plasic_dt,
            plasic_dq, plasic_du, plasic_dv, plasic_dp,
            plasic_dudt, plasic_dvdt, plasic_dtdt,
            plasic_dtaux, plasic_dtauy, plasic_dust3,
            work->horizontal_f, work->horizontal_g);
    }
    if (flux.sensible_heat != 0)
    {
        flux_mkshfl(
            PLASIC_NHOR, low_level, surface_level, plasic_kap,
            plasic_sigma[PLASIC_NLEV - 1], plasic_ga,
            plasic_deltsec2, plasic_gascon,
            plasic_dsigma[PLASIC_NLEV - 1], plasic_cpd,
            plasic_cpv_cpd_minus1, flux.heat_transfer, plasic_dt,
            plasic_dq, plasic_dp, plasic_dtdt,
            plasic_dshfl, plasic_dshdt, plasic_dtsa,
            work->horizontal_f, work->horizontal_g);
    }
    if (flux.evaporation != 0)
    {
        flux_mkevap(
            PLASIC_NHOR, low_level, surface_level, plasic_ga,
            plasic_deltsec, plasic_deltsec2,
            plasic_gascon, plasic_dsigma[PLASIC_NLEV - 1],
            plasic_tmelt, plasic_lv, plasic_ls,
            plasic_ra2, plasic_ra4, flux.heat_transfer,
            plasic_dt, plasic_dq, plasic_dp,
            plasic_dwetfac, plasic_dls, plasic_dwatc,
            plasic_dqdt, plasic_devap, plasic_dlhfl,
            plasic_dlhdt, work->horizontal_f);
    }
}

static void vertical_diffusion_step(void)
{
    if (flux.vertical_diffusion == 0)
    {
        return;
    }
    vdiff(
        PLASIC_NHOR, PLASIC_NLEV, plasic_ndheat,
        flux.minimum_wind, flux.mixing_length, flux.stability_b,
        flux.stability_c, flux.stability_d, plasic_deltsec2,
        plasic_ga, plasic_gascon, plasic_kap,
        plasic_rdbrv, plasic_cpd, plasic_cpv_cpd_minus1,
        plasic_sigma, plasic_sigmah, plasic_dsigma,
        plasic_dt, plasic_dq, plasic_du,
        plasic_dv, plasic_dtdt, plasic_dqdt,
        plasic_dudt, plasic_dvdt);
}

static int prepare_convective_humidity_tendency(void)
{
    int status = sht_local_scalar_to_grid(
        PLASIC_NLEV, plasic_sqt, PLASIC_NSPP, plasic_dqt);

    if (status != PLASIC_SHT_OK)
    {
        return status == PLASIC_SHT_INVALID_ARGUMENT
                   ? PLASIC_RUNTIME_BAD_ARGUMENT
               : status == PLASIC_SHT_ALLOCATION_FAILED
                   ? PLASIC_RUNTIME_INITIALIZATION_FAILED
                   : PLASIC_RUNTIME_NUMERICAL_FAILURE;
    }
    rain_scale_dqt(
        PLASIC_NHOR, PLASIC_NLEV, plasic_psurf, plasic_ww,
        plasic_dp, plasic_dqt);
    return PLASIC_RUNTIME_OK;
}

static int rain_step(void)
{
    int status;
    int level;

    rain_zero_precipitation(
        PLASIC_NHOR, plasic_dprl, plasic_dprc,
        plasic_dprs);
    if (rain.convective_rain != 0)
    {
        status = prepare_convective_humidity_tendency();
        if (status != PLASIC_RUNTIME_OK)
        {
            return status;
        }
        status = rain_kuo(
            PLASIC_NHOR, PLASIC_NLEV, rain.beta_mode,
            rain.surface_convection, rain.shallow_convection,
            rain.mix_momentum, rain.beta_exponent, plasic_ndheat,
            rain.prescribed_beta, rain.critical_relative_humidity,
            rain.deep_pressure, rain.shallow_top_pressure,
            rain.shallow_diffusion, plasic_deltsec2,
            plasic_kap, plasic_rdbrv, plasic_ra1,
            plasic_ra2, plasic_tmelt, plasic_ra4,
            plasic_ls, plasic_lv, plasic_cpd,
            plasic_cpv_cpd_minus1, plasic_ga, plasic_gascon,
            plasic_sigma, plasic_sigmah, plasic_dsigma,
            plasic_dp, plasic_dt, plasic_dq,
            plasic_du, plasic_dv, plasic_dqt,
            plasic_dtdt, plasic_dqdt, plasic_dudt,
            plasic_dvdt, rain.convective_layer,
            rain.convective_layer_count, rain.convective_rain_layer,
            rain.convective_snow_layer);
        if (status != 0)
        {
            return PLASIC_RUNTIME_NUMERICAL_FAILURE;
        }
    }
    if (rain.dry_adjustment != 0)
    {
        status = rain_dry_adjustment(
            PLASIC_NHOR, PLASIC_NLEV, plasic_kap,
            plasic_deltsec2, plasic_sigma, plasic_dsigma,
            plasic_dt, plasic_dq, plasic_dtdt,
            plasic_dqdt);
        if (status != 0)
        {
            return PLASIC_RUNTIME_NUMERICAL_FAILURE;
        }
    }
    if (rain.large_scale_rain != 0)
    {
        memset(rain.large_scale_rain_layer, 0,
               sizeof(rain.large_scale_rain_layer));
        memset(rain.large_scale_snow_layer, 0,
               sizeof(rain.large_scale_snow_layer));
        for (level = 0; level < PLASIC_NLEV; ++level)
        {
            const size_t offset =
                (size_t)level * (size_t)PLASIC_NHOR;
            rain_lsp_prepare(
                PLASIC_NHOR, plasic_sigma[level],
                plasic_deltsec2, plasic_lv,
                plasic_cpd, plasic_cpv_cpd_minus1, plasic_rdbrv,
                plasic_ra1, plasic_ra2, plasic_tmelt,
                plasic_ra4, plasic_dq + offset,
                plasic_dqdt + offset, plasic_dt + offset,
                plasic_dtdt + offset, plasic_dp,
                work->horizontal_a, work->horizontal_b,
                work->horizontal_c, work->horizontal_d,
                work->horizontal_e, work->horizontal_f);
            rain_lsp_finish(
                PLASIC_NHOR, plasic_dsigma[level],
                plasic_sigma[level], plasic_deltsec2,
                plasic_rdbrv, plasic_ra1, plasic_ra2,
                plasic_tmelt, plasic_ra4, plasic_dp,
                work->horizontal_a, work->horizontal_b,
                work->horizontal_c, work->horizontal_d,
                work->horizontal_e, work->horizontal_f,
                plasic_dqdt + offset, plasic_dtdt + offset,
                rain.large_scale_rain_layer + offset);
        }
    }
    status = rain_fall(
        PLASIC_NHOR, PLASIC_NLEV, rain.evaporate_precipitation,
        rain.evaporation_rate, plasic_deltsec2, plasic_ls,
        plasic_lv, plasic_cpd, plasic_cpv_cpd_minus1,
        plasic_tmelt, plasic_rdbrv, plasic_ra1,
        plasic_ra2, plasic_ra4, plasic_ga,
        plasic_sigma, plasic_dsigma, plasic_dp,
        plasic_dt, plasic_dq, rain.convective_rain_layer,
        rain.large_scale_rain_layer, rain.convective_snow_layer,
        rain.large_scale_snow_layer, plasic_dtdt,
        plasic_dqdt, plasic_dprc, plasic_dprl,
        plasic_dprs);
    if (status != 0)
    {
        return PLASIC_RUNTIME_NUMERICAL_FAILURE;
    }
    if (rain.clouds != 0)
    {
        status = rain_clouds(
            PLASIC_NHOR, PLASIC_NLEV, plasic_solar_day,
            plasic_deltsec2, plasic_rdbrv, plasic_ra1,
            plasic_ra2, plasic_tmelt, plasic_ra4,
            rain.cloud_factor, rain.cloud_vertical_velocity_1,
            rain.cloud_vertical_velocity_2, plasic_gascon,
            plasic_ga, rain.convective_layer_count,
            rain.convective_layer, plasic_dprc, plasic_dp,
            plasic_sigma, plasic_sigmah, plasic_dsigma,
            rain.critical_humidity, plasic_dt, plasic_dtdt,
            plasic_dq, plasic_dqdt, plasic_dw,
            plasic_dcc, plasic_dql, plasic_dqvi);
        if (status != 0)
        {
            return PLASIC_RUNTIME_NUMERICAL_FAILURE;
        }
    }
    return PLASIC_RUNTIME_OK;
}

void physics_interpolate_land_surface(void)
{
    int32_t first_month;
    int32_t second_month;
    int32_t next_step = plasic_nstep + 1;
    float weight;

    momint(
        plasic_calendar, next_step, &first_month,
        &second_month, &weight);
    if (land.prognostic_temperature == 0 ||
        land.prognostic_water == 0)
    {
        land_interpolate_climatology(
            PLASIC_NHOR, first_month, second_month, weight,
            land.temperature_climatology, land.wetness_climatology,
            land.interpolated_temperature, land.interpolated_wetness);
    }
    land_interpolate_albedo(
        PLASIC_NHOR, first_month, second_month, weight,
        land.albedo_climatology, land.interpolated_albedo);
}

static int land_step(void)
{
    int status;

    physics_interpolate_land_surface();
    if (land.prognostic_temperature != 0)
    {
        status = land_tands(
            PLASIC_NHOR, PLASIC_LAND_SOIL_LEVELS, plasic_deltsec,
            land.top_layer_depth, land.maximum_snow_depth,
            plasic_tmelt, plasic_ls, plasic_lv,
            land.snow_density, land.soil_conductivity,
            land.ice_conductivity, land.snow_conductivity,
            land.soil_heat_capacity, land.ice_heat_capacity,
            land.snow_heat_capacity, land.soil_depth, plasic_dls,
            plasic_dglac, plasic_dshfl, plasic_dlhfl,
            plasic_dflux + (size_t)PLASIC_NHOR * PLASIC_NLEV,
            plasic_dprs, plasic_devap, plasic_dprl,
            plasic_dprc, land.previous_surface_temperature,
            land.soil_temperature, land.surface_temperature,
            land.snow_depth, land.snow_temperature, land.soil_water_flux,
            plasic_dsmelt, plasic_dsndch, work->horizontal_a,
            work->horizontal_b, work->horizontal_c);
        if (status != 0)
        {
            return PLASIC_RUNTIME_NUMERICAL_FAILURE;
        }
    }
    else
    {
        land_prescribed_temperature(
            PLASIC_NHOR, plasic_dls, plasic_devap,
            plasic_dprl, plasic_dprc,
            land.interpolated_temperature, land.soil_water_flux,
            land.surface_temperature);
    }
    if (land.prognostic_water != 0)
    {
        land_update_soil_water(
            PLASIC_NHOR, plasic_deltsec, plasic_dls,
            land.soil_water_flux, plasic_dwmax, plasic_dwatc);
    }
    else
    {
        land_prescribed_water(
            PLASIC_NHOR, plasic_dls, land.interpolated_wetness,
            plasic_dwatc);
    }

    land_surface_exchange(
        PLASIC_NHOR, PLASIC_LAND_SOIL_LEVELS, 0, 0,
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
        PLASIC_NHOR, 0, land.maximum_snow_depth,
        land.glacier_albedo_minimum, land.glacier_albedo_maximum,
        plasic_tmelt, plasic_dls, plasic_dglac,
        land.surface_temperature, land.snow_depth, plasic_dsnow,
        plasic_dalb, plasic_dwetfac);
    return PLASIC_RUNTIME_OK;
}

static int ocean3d_exchange_halo(float *field, int32_t levels);

static void ocean3d_fill_config(deep_ocean_kernel_config *config)
{
    config->exchange_halo = ocean3d_exchange_halo;
    config->nlev = deep_ocean.nlev;
    config->nlon = PLASIC_NLON;
    config->nlat_local = PLASIC_NLPP;
    config->dt = deep_ocean.timestep;
    config->radius = deep_ocean.radius;
    config->solar_day = deep_ocean.solar_day;
    config->sidereal_day = deep_ocean.sidereal_day;
    config->gravity = deep_ocean.gravity;
    config->cp = deep_ocean.cp;
    config->rho0 = deep_ocean.rho0;
    config->alpha = deep_ocean.alpha;
    config->beta = deep_ocean.beta;
    config->kappa_v = deep_ocean.kappa_v;
    config->kappa_h = deep_ocean.kappa_h;
    config->t_freeze = deep_ocean.t_freeze;
}

void physics_ocean3d_pack_local(
    const float *extended, float *local, int32_t levels)
{
    const int32_t nlon = PLASIC_NLON;
    const int32_t rows = PLASIC_NLPP;
    int32_t level;
    int32_t row;
    int32_t i;

    for (level = 0; level < levels; ++level)
    {
        const size_t ext_base = (size_t)level * PLASIC_DEEP_EXT_HOR;
        const size_t local_base = (size_t)level * PLASIC_NHOR;

        for (row = 0; row < rows; ++row)
        {
            for (i = 0; i < nlon; ++i)
            {
                local[local_base + (size_t)(i + row * nlon)] =
                    extended[ext_base + (size_t)(i + (row + 1) * nlon)];
            }
        }
    }
}

void physics_ocean3d_unpack_local(
    const float *local, float *extended, int32_t levels)
{
    const int32_t nlon = PLASIC_NLON;
    const int32_t rows = PLASIC_NLPP;
    int32_t level;
    int32_t row;
    int32_t i;

    for (level = 0; level < levels; ++level)
    {
        const size_t ext_base = (size_t)level * PLASIC_DEEP_EXT_HOR;
        const size_t local_base = (size_t)level * PLASIC_NHOR;

        for (row = 0; row < rows; ++row)
        {
            for (i = 0; i < nlon; ++i)
            {
                extended[ext_base + (size_t)(i + (row + 1) * nlon)] =
                    local[local_base + (size_t)(i + row * nlon)];
            }
        }
    }
}

/*
 * Exchange the latitude ghost rows of a halo-extended ocean field.  The ghost
 * buffers are pre-filled with the mirrored interior edge so that the global
 * polar ranks keep a zero-gradient condition.
 */
static int ocean3d_exchange_halo(float *field, int32_t levels)
{
    const int32_t nlon = PLASIC_NLON;
    const int32_t rows = PLASIC_NLPP;
    const int32_t elements = nlon * levels;
    const size_t bytes = (size_t)elements * sizeof(float);
    float *send_north = (float *)malloc(bytes);
    float *send_south = (float *)malloc(bytes);
    float *recv_north = (float *)malloc(bytes);
    float *recv_south = (float *)malloc(bytes);
    int32_t level;
    int32_t i;

    if (send_north == NULL || send_south == NULL || recv_north == NULL ||
        recv_south == NULL)
    {
        free(send_north);
        free(send_south);
        free(recv_north);
        free(recv_south);
        return -1;
    }
    for (level = 0; level < levels; ++level)
    {
        const size_t base = (size_t)level * PLASIC_DEEP_EXT_HOR;
        const size_t packed = (size_t)level * (size_t)nlon;

        for (i = 0; i < nlon; ++i)
        {
            send_north[packed + (size_t)i] =
                field[base + (size_t)(i + nlon)];
            send_south[packed + (size_t)i] =
                field[base + (size_t)(i + rows * nlon)];
        }
    }
    memcpy(recv_north, send_north, bytes);
    memcpy(recv_south, send_south, bytes);
    if (mp_exchange_latitude_halo(
            recv_north, recv_south, send_north, send_south, elements) != 0)
    {
        free(send_north);
        free(send_south);
        free(recv_north);
        free(recv_south);
        return -1;
    }
    for (level = 0; level < levels; ++level)
    {
        const size_t base = (size_t)level * PLASIC_DEEP_EXT_HOR;
        const size_t packed = (size_t)level * (size_t)nlon;

        for (i = 0; i < nlon; ++i)
        {
            field[base + (size_t)i] = recv_north[packed + (size_t)i];
            field[base + (size_t)(i + (rows + 1) * nlon)] =
                recv_south[packed + (size_t)i];
        }
    }
    free(send_north);
    free(send_south);
    free(recv_north);
    free(recv_south);
    return 0;
}

static void ocean3d_build_land_mask(void)
{
    const int32_t nlon = PLASIC_NLON;
    const int32_t rows = PLASIC_NLPP;
    int32_t row;
    int32_t i;

    for (row = 0; row < rows; ++row)
    {
        for (i = 0; i < nlon; ++i)
        {
            deep_ocean.land_mask[i + (row + 1) * nlon] =
                ocean.land_mask[i + row * nlon];
        }
    }
    for (i = 0; i < nlon; ++i)
    {
        deep_ocean.land_mask[i] = deep_ocean.land_mask[i + nlon];
        deep_ocean.land_mask[i + (rows + 1) * nlon] =
            deep_ocean.land_mask[i + rows * nlon];
    }
    (void)ocean3d_exchange_halo(deep_ocean.land_mask, 1);
}

/*
 * Exchange a single scalar per latitude row (used for the local latitude
 * coordinate).  This is deliberately separate from the horizontal field
 * exchange above, whose row width is the number of longitudes.
 */
static void ocean3d_exchange_scalar_halo(float *coordinate)
{
    const int32_t rows = PLASIC_NLPP;
    float send_north = coordinate[1];
    float send_south = coordinate[rows];
    float recv_north = send_north;
    float recv_south = send_south;

    if (mp_exchange_latitude_halo(
            &recv_north, &recv_south, &send_north, &send_south, 1) != 0)
    {
        return;
    }
    coordinate[0] = recv_north;
    coordinate[rows + 1] = recv_south;
}

/*
 * Build the local, halo-extended latitude coordinates.  PlaSiC moves each
 * rank's latitude block to the front of `plasic_sid`/`plasic_cola`, so only
 * indices 0..NLPP-1 are meaningful.  Interior ghosts are exchanged with the
 * neighbouring ranks; the global polar ghosts are linearly extrapolated so
 * the centred meridional difference stays finite and physical.
 */
static void ocean3d_build_latitudes(void)
{
    const int32_t rows = PLASIC_NLPP;
    int32_t i;

    for (i = 0; i < rows; ++i)
    {
        deep_ocean.sine_latitude[i + 1] = (float)plasic_sid[i];
        deep_ocean.cosine_latitude[i + 1] = plasic_cola[i];
    }
    deep_ocean.sine_latitude[0] = deep_ocean.sine_latitude[1];
    deep_ocean.sine_latitude[rows + 1] = deep_ocean.sine_latitude[rows];
    deep_ocean.cosine_latitude[0] = deep_ocean.cosine_latitude[1];
    deep_ocean.cosine_latitude[rows + 1] = deep_ocean.cosine_latitude[rows];
    ocean3d_exchange_scalar_halo(deep_ocean.sine_latitude);
    ocean3d_exchange_scalar_halo(deep_ocean.cosine_latitude);
    if (mp_rank() == 0)
    {
        deep_ocean.sine_latitude[0] =
            2.0f * deep_ocean.sine_latitude[1] -
            deep_ocean.sine_latitude[2];
        deep_ocean.cosine_latitude[0] = deep_ocean.cosine_latitude[1];
    }
    if (mp_rank() == mp_size() - 1)
    {
        deep_ocean.sine_latitude[rows + 1] =
            2.0f * deep_ocean.sine_latitude[rows] -
            deep_ocean.sine_latitude[rows - 1];
        deep_ocean.cosine_latitude[rows + 1] =
            deep_ocean.cosine_latitude[rows];
    }
}

void physics_ocean3d_initialize(int from_restart)
{
    deep_ocean_kernel_config config;

    if (deep_ocean.enabled == 0)
    {
        return;
    }
    deep_ocean.timestep = ocean.timestep;
    ocean3d_build_land_mask();
    ocean3d_build_latitudes();
    if (from_restart != 0 && deep_ocean.restart_loaded != 0)
    {
        /* Interior rows were restored from the restart; refresh ghosts. */
    }
    else
    {
        ocean3d_fill_config(&config);
        deep_ocean_seed_state(
            &config, deep_ocean.land_mask, ocean.temperature,
            deep_ocean.layer_thickness, deep_ocean.temperature,
            deep_ocean.salinity, deep_ocean.u_velocity,
            deep_ocean.v_velocity, deep_ocean.w_velocity);
    }
    (void)ocean3d_exchange_halo(deep_ocean.temperature, PLASIC_DEEP_NLEV);
    (void)ocean3d_exchange_halo(deep_ocean.salinity, PLASIC_DEEP_NLEV);
    (void)ocean3d_exchange_halo(deep_ocean.u_velocity, PLASIC_DEEP_NLEV);
    (void)ocean3d_exchange_halo(deep_ocean.v_velocity, PLASIC_DEEP_NLEV);
    (void)ocean3d_exchange_halo(deep_ocean.w_velocity, PLASIC_DEEP_NLEV);
    deep_ocean.initialized = 1;
}

void physics_ocean3d_advance(void)
{
    deep_ocean_kernel_config config;

    if (deep_ocean.enabled == 0 || deep_ocean.initialized == 0)
    {
        return;
    }
    (void)ocean3d_exchange_halo(deep_ocean.temperature, PLASIC_DEEP_NLEV);
    (void)ocean3d_exchange_halo(deep_ocean.salinity, PLASIC_DEEP_NLEV);
    (void)ocean3d_exchange_halo(deep_ocean.u_velocity, PLASIC_DEEP_NLEV);
    (void)ocean3d_exchange_halo(deep_ocean.v_velocity, PLASIC_DEEP_NLEV);
    (void)ocean3d_exchange_halo(deep_ocean.w_velocity, PLASIC_DEEP_NLEV);
    ocean3d_fill_config(&config);
    deep_ocean_step(
        &config, deep_ocean.sine_latitude, deep_ocean.cosine_latitude,
        deep_ocean.land_mask, ocean.temperature,
        ocean.precipitation_minus_evaporation, ocean.stress_x,
        ocean.stress_y, deep_ocean.layer_thickness,
        deep_ocean.temperature, deep_ocean.salinity,
        deep_ocean.u_velocity, deep_ocean.v_velocity,
        deep_ocean.w_velocity, deep_ocean.density, deep_ocean.qflux);
    memcpy(ocean.deep_ocean_flux, deep_ocean.qflux,
           sizeof(ocean.deep_ocean_flux));
}

static void ocean_step(void)
{
    double temperature[PLASIC_NHOR];
    int output_due;

    ocean_copy_inputs(
        PLASIC_NHOR, 920.0f, ice.ice_cover, ice.ice_thickness,
        ice.ocean_heat_accumulator, ice.ocean_freshwater_accumulator,
        ice.ocean_stress_x_accumulator, ice.ocean_stress_y_accumulator,
        ice.ocean_friction_accumulator, ice.ocean_snow_accumulator,
        ice.current_climatological_thickness, ocean.ice_cover,
        ocean.ice_thickness, ocean.heat,
        ocean.precipitation_minus_evaporation, ocean.stress_x,
        ocean.stress_y, ocean.friction_velocity_cubed,
        ocean.ice_and_snow, ocean.climatological_ice_thickness);
    physics_interpolate_ocean_climatology();

    /*
     * Run the three-dimensional deep ocean before the slab so that the slab
     * heat budget below can use the deep-ocean heat flux computed here.
     */
    physics_ocean3d_advance();

    ocean_mksst_begin(
        PLASIC_NHOR, ocean.timestep, 1030.0f,
        4180.0f, 271.25f, ocean.land_mask, ocean.ice_thickness,
        ocean.heat, ocean.deep_ocean_flux, ocean.mixed_layer_depth,
        ocean.temperature, temperature, ocean.ice_flux);
    /*
     * The mixed-layer ocean is forced to a single layer with no vertical or
     * horizontal diffusion, so the slab only exchanges heat through the
     * atmosphere, sea ice and deep-ocean fluxes above.
     */
    ocean_mksst_finish(
        PLASIC_NHOR, ocean.prognostic_ocean,
        ocean.flux_correction, ocean.timestep, 1030.0f, 4180.0f,
        271.25f, ocean.land_mask, ocean.ice_thickness,
        ocean.mixed_layer_depth, temperature, ocean.ice_flux,
        ocean.temperature);

    output_due = ocean_finish_step(
        PLASIC_NHOR, ocean.prognostic_ocean, ocean.flux_correction,
        plasic_nstep, ocean.output_interval, ocean.timestep,
        ocean.cooling_timescale, 1030.0f, 920.0f, 4180.0f, 3.28e5f,
        271.25f, ocean.land_mask, ocean.interpolated_sst,
        ocean.mixed_layer_depth, ocean.ice_thickness,
        ocean.climatological_ice_thickness, ocean.heat,
        ocean.deep_ocean_flux, ocean.temperature,
        ocean.ice_flux, ocean.residual_ice_flux, ocean.heat_accumulator,
        ocean.ice_flux_accumulator,
        ocean.deep_ocean_flux_accumulator,
        &ocean.accumulation_count);
    if (output_due != 0)
    {
        ocean_reset_diagnostics(
            PLASIC_NHOR, ocean.heat_accumulator,
            ocean.ice_flux_accumulator,
            ocean.deep_ocean_flux_accumulator,
            &ocean.accumulation_count);
    }
    ocean_copy_outputs(
        PLASIC_NHOR, ocean.temperature, ocean.mixed_layer_depth,
        ocean.ice_flux, ice.sst, ice.mixed_layer_depth, ice.ocean_heat);
}

static void reset_ice_output_accumulators(void)
{
    memset(ice.flux_correction_accumulator, 0,
           sizeof(ice.flux_correction_accumulator));
    memset(ice.atmospheric_heat_accumulator, 0,
           sizeof(ice.atmospheric_heat_accumulator));
    memset(ice.ocean_heat_diagnostic_accumulator, 0,
           sizeof(ice.ocean_heat_diagnostic_accumulator));
    memset(ice.residual_melt_accumulator, 0,
           sizeof(ice.residual_melt_accumulator));
    memset(ice.conductive_flux_accumulator, 0,
           sizeof(ice.conductive_flux_accumulator));
    memset(ice.maximum_thickness_flux_accumulator, 0,
           sizeof(ice.maximum_thickness_flux_accumulator));
    memset(ice.negative_ice_flux_accumulator, 0,
           sizeof(ice.negative_ice_flux_accumulator));
    memset(ice.snow_melt_accumulator, 0,
           sizeof(ice.snow_melt_accumulator));
    memset(ice.ice_melt_accumulator, 0,
           sizeof(ice.ice_melt_accumulator));
    memset(ice.surface_storage_accumulator, 0,
           sizeof(ice.surface_storage_accumulator));
    memset(ice.diagnosed_conductive_accumulator, 0,
           sizeof(ice.diagnosed_conductive_accumulator));
    memset(ice.snow_conversion_accumulator, 0,
           sizeof(ice.snow_conversion_accumulator));
    memset(ice.freshwater_diagnostic_accumulator, 0,
           sizeof(ice.freshwater_diagnostic_accumulator));
    memset(ice.snow_to_ice_accumulator, 0,
           sizeof(ice.snow_to_ice_accumulator));
    ice.output_accumulation_count = 0;
}

static void accumulate_ice_diagnostics(void)
{
    size_t horizontal;

    for (horizontal = 0; horizontal < PLASIC_NHOR; ++horizontal)
    {
        ice.flux_correction_accumulator[horizontal] +=
            ice.current_flux_correction[horizontal];
        ice.atmospheric_heat_accumulator[horizontal] +=
            ice.atmospheric_heat[horizontal];
        ice.ocean_heat_diagnostic_accumulator[horizontal] +=
            ice.ocean_heat[horizontal];
        ice.residual_melt_accumulator[horizontal] +=
            ice.residual_melt_energy[horizontal];
        ice.ice_melt_accumulator[horizontal] +=
            ice.ice_melt_flux[horizontal];
        ice.snow_melt_accumulator[horizontal] +=
            ice.snow_melt_flux[horizontal];
        ice.conductive_flux_accumulator[horizontal] +=
            ice.conductive_flux[horizontal];
        ice.maximum_thickness_flux_accumulator[horizontal] +=
            ice.maximum_thickness_flux[horizontal];
        ice.negative_ice_flux_accumulator[horizontal] +=
            ice.negative_ice_flux[horizontal];
        ice.surface_storage_accumulator[horizontal] +=
            ice.surface_storage_flux[horizontal];
        ice.diagnosed_conductive_accumulator[horizontal] +=
            ice.diagnosed_conductive_flux[horizontal];
        ice.snow_conversion_accumulator[horizontal] +=
            ice.snow_conversion_flux[horizontal];
        ice.freshwater_diagnostic_accumulator[horizontal] +=
            ice.freshwater_flux[horizontal];
        ice.snow_to_ice_accumulator[horizontal] +=
            ice.snow_to_ice[horizontal];
    }
    ++ice.output_accumulation_count;
    if (ice.output_interval > 0 &&
        plasic_nstep % ice.output_interval == 0)
    {
        /*
         * The Fortran module forms averages for model output and immediately
         * clears them.  This standalone executable currently writes restart
         * and comparison diagnostics instead of legacy unit-71 output.
         */
        reset_ice_output_accumulators();
    }
}

static void ice_step(void)
{
    const float freezing_temperature = 271.25f;
    const float seawater_density = 1030.0f;
    const float ice_density = 920.0f;
    const float snow_density = 330.0f;
    const float ice_heat_capacity = 2070.0f;
    const float ice_conductivity = 2.03f;
    const float snow_conductivity = 0.31f;
    const float latent_heat_ice = 3.28e5f;
    const float latent_heat_snow = 3.337e5f;
    const float ice_energy_per_timestep =
        ice_density * latent_heat_ice / ice.timestep;
    size_t horizontal;

    memset(ice.ice_melt_flux, 0, sizeof(ice.ice_melt_flux));
    memset(ice.snow_melt, 0, sizeof(ice.snow_melt));
    memset(ice.snow_melt_flux, 0, sizeof(ice.snow_melt_flux));
    memset(ice.surface_storage_flux, 0, sizeof(ice.surface_storage_flux));
    memset(ice.diagnosed_conductive_flux, 0,
           sizeof(ice.diagnosed_conductive_flux));
    memset(ice.snow_conversion_flux, 0,
           sizeof(ice.snow_conversion_flux));
    memset(ice.negative_ice_flux, 0, sizeof(ice.negative_ice_flux));
    memset(ice.snow_change, 0, sizeof(ice.snow_change));

    memcpy(ice.atmospheric_heat, sea.heat_accumulator,
           sizeof(ice.atmospheric_heat));
    memcpy(ice.modified_atmospheric_heat, ice.atmospheric_heat,
           sizeof(ice.modified_atmospheric_heat));
    memcpy(ice.freshwater_flux, sea.freshwater_accumulator,
           sizeof(ice.freshwater_flux));
    memcpy(ice.snowfall, sea.snowfall_accumulator, sizeof(ice.snowfall));
    memcpy(ice.stress_x, sea.stress_x_accumulator, sizeof(ice.stress_x));
    memcpy(ice.stress_y, sea.stress_y_accumulator, sizeof(ice.stress_y));
    memcpy(ice.friction_velocity_cubed, sea.friction_accumulator,
           sizeof(ice.friction_velocity_cubed));
    if (ice.maximum_thickness_correction > 0)
    {
        for (horizontal = 0; horizontal < PLASIC_NHOR; ++horizontal)
        {
            ice.modified_atmospheric_heat[horizontal] -=
                ice.maximum_thickness_flux[horizontal];
        }
    }
    physics_interpolate_ice_climatology();

    memcpy(work->horizontal_a, ice.snow, sizeof(ice.snow));
    memset(work->horizontal_b, 0, sizeof(ice.ice_thickness));
    for (horizontal = 0; horizontal < PLASIC_NHOR; ++horizontal)
    {
        if (ice.land_mask[horizontal] < 0.5f)
        {
            work->horizontal_b[horizontal] = ice.ice_thickness[horizontal];
        }
        if (ice.prognostic_snow == 1 &&
            ice.ice_cover[horizontal] >= ice.minimum_ice_cover)
        {
            ice.snow[horizontal] += ice.snowfall[horizontal] * ice.timestep;
            ice.snowfall[horizontal] = 0.0f;
        }
    }

    /*
     * ice.compactness is the prognostic compactness c_i^*: it is carried
     * from step to step (restart record xicecc) and only updated by
     * ice_compactness() below.  ice.ice_cover is the thresholded 0/1
     * mask derived from it and must not be copied back here, otherwise
     * sub-threshold growth is lost every step and the ice edge can only
     * retreat.
     */

    if (ice.prognostic_skin_temperature == 1)
    {
        ice_skin_temperature(
            PLASIC_NHOR, ice.minimum_thickness, ice.maximum_thickness,
            plasic_tmelt, ice_density, ice_heat_capacity,
            snow_density, ice_conductivity, snow_conductivity,
            ice.timestep, ice.ice_thickness, ice.snow, ice.sst,
            ice.modified_atmospheric_heat, ice.surface_temperature,
            ice.surface_storage_flux, ice.residual_melt_energy);
    }
    else
    {
        for (horizontal = 0; horizontal < PLASIC_NHOR; ++horizontal)
        {
            ice.surface_temperature[horizontal] =
                ice.ice_thickness[horizontal] +
                            1000.0f / snow_density * ice.snow[horizontal] >=
                        0.1f ?
                    ice.current_climatological_sst[horizontal] :
                    ice.sst[horizontal];
        }
    }
    ice_conductive_flux(
        PLASIC_NHOR, ocean.flux_correction, ice.prognostic_ice,
        ice.minimum_thickness, ice.maximum_thickness,
        freezing_temperature, snow_density, ice_conductivity,
        snow_conductivity, ice.ice_thickness, ice.snow,
        ice.surface_temperature, ice.sst, ice.modified_atmospheric_heat,
        ice.diagnosed_conductive_flux, ice.conductive_flux,
        ice.freezing_flux);
    if (ice.prognostic_snow == 1)
    {
        ice_snow(
            PLASIC_NHOR, ice.minimum_ice_cover, seawater_density,
            ice_density, latent_heat_snow, ice.timestep, ice.land_mask,
            ice.ice_cover, ice.snowfall, ice.residual_melt_energy,
            ice.ice_thickness, ice.snow, ice.conductive_flux,
            ice.snow_melt, ice.snow_melt_flux,
            ice.snow_conversion_flux, ice.snow_to_ice);
    }
    else
    {
        memset(ice.snow, 0, sizeof(ice.snow));
    }
    ice_make_ice(
        PLASIC_NHOR, freezing_temperature, ice_density,
        latent_heat_ice, ice.timestep, ice.land_mask,
        ice.modified_ocean_heat, ice.sst, ice.snow_conversion_flux,
        ice.freezing_flux, ice.conductive_flux, ice.ice_melt_flux,
        ice.negative_ice_flux, ice.ice_thickness);
    ice_compactness(
        PLASIC_NHOR, work->horizontal_b, ice.ice_thickness,
        ice.compactness);

    if (ice.prognostic_ice == 0)
    {
        for (horizontal = 0; horizontal < PLASIC_NHOR; ++horizontal)
        {
            if (ice.ice_thickness[horizontal] > 0.0f ||
                ice.current_climatological_thickness[horizontal] > 0.0f)
            {
                ice.current_flux_correction[horizontal] =
                    (ice.ice_thickness[horizontal] -
                     ice.current_climatological_thickness[horizontal]) *
                        ice_energy_per_timestep -
                    ice.conductive_flux[horizontal];
                ice.conductive_flux[horizontal] = 0.0f;
                ice.ice_thickness[horizontal] =
                    ice.current_climatological_thickness[horizontal];
                ice.ice_cover[horizontal] =
                    ice.current_climatological_cover[horizontal];
            }
            else
            {
                ice.current_flux_correction[horizontal] = 0.0f;
                ice.ice_cover[horizontal] = 0.0f;
            }
        }
    }
    else
    {
        memcpy(ice.ice_cover, ice.compactness,
               sizeof(ice.ice_cover));
        if (ocean.flux_correction == 0)
        {
            memset(ice.current_flux_correction, 0,
                   sizeof(ice.current_flux_correction));
        }
    }

    memset(work->horizontal_c, 0, sizeof(ice.ice_thickness));
    memset(ice.maximum_thickness_flux, 0,
           sizeof(ice.maximum_thickness_flux));
    for (horizontal = 0; horizontal < PLASIC_NHOR; ++horizontal)
    {
        if (ice.maximum_thickness >= 0.0f &&
            ice.ice_thickness[horizontal] > ice.maximum_thickness)
        {
            work->horizontal_c[horizontal] =
                (ice.ice_thickness[horizontal] - ice.maximum_thickness) *
                ice_energy_per_timestep;
            ice.ice_thickness[horizontal] = ice.maximum_thickness;
            ice.maximum_thickness_flux[horizontal] +=
                work->horizontal_c[horizontal];
            ice.ice_melt_flux[horizontal] +=
                work->horizontal_c[horizontal];
        }
        ice.ice_cover[horizontal] =
            ice.ice_cover[horizontal] >= ice.compactness_threshold ?
                1.0f : 0.0f;
        if (ice.ice_thickness[horizontal] <= 0.0f &&
            ice.snow[horizontal] > 0.0f)
        {
            const float flux_value =
                ice.snow[horizontal] * 1000.0f * latent_heat_snow /
                ice.timestep;
            ice.snow_melt[horizontal] +=
                ice.snow[horizontal] / ice.timestep;
            ice.snow_melt_flux[horizontal] += flux_value;
            ice.conductive_flux[horizontal] -= flux_value;
            ice.snow[horizontal] = 0.0f;
        }
    }
    accumulate_ice_diagnostics();

    for (horizontal = 0; horizontal < PLASIC_NHOR; ++horizontal)
    {
        if (ice.land_mask[horizontal] < 0.5f)
        {
            ice.ocean_heat_accumulator[horizontal] +=
                ice.conductive_flux[horizontal];
            ice.ocean_freshwater_accumulator[horizontal] +=
                ice.freshwater_flux[horizontal];
            ice.ocean_stress_x_accumulator[horizontal] +=
                ice.stress_x[horizontal];
            ice.ocean_stress_y_accumulator[horizontal] +=
                ice.stress_y[horizontal];
            ice.ocean_friction_accumulator[horizontal] +=
                ice.friction_velocity_cubed[horizontal];
            ice.ocean_snow_accumulator[horizontal] += ice.snow[horizontal];
        }
    }
    ++ice.ocean_accumulation_count;
    if (plasic_nstep % ice.ocean_coupling_interval == 0)
    {
        const float inverse_count =
            1.0f / (float)ice.ocean_accumulation_count;

        /*
         * The ocean, including the deep ocean, integrates all the sea-ice
         * steps it is about to average.  With the fixed
         * ocean_coupling_interval = 1 this is exactly the current ice step.
         */
        ocean.timestep = (float)ice.ocean_accumulation_count * ice.timestep;
        deep_ocean.timestep = ocean.timestep;
        for (horizontal = 0; horizontal < PLASIC_NHOR; ++horizontal)
        {
            if (ice.land_mask[horizontal] < 0.5f)
            {
                ice.ocean_heat_accumulator[horizontal] *= inverse_count;
                ice.ocean_freshwater_accumulator[horizontal] *= inverse_count;
                ice.ocean_stress_x_accumulator[horizontal] *= inverse_count;
                ice.ocean_stress_y_accumulator[horizontal] *= inverse_count;
                ice.ocean_friction_accumulator[horizontal] *= inverse_count;
                ice.ocean_snow_accumulator[horizontal] *= inverse_count;
            }
        }
        ocean_step();
        memcpy(ice.modified_ocean_heat, ice.ocean_heat,
               sizeof(ice.modified_ocean_heat));
        for (horizontal = 0; horizontal < PLASIC_NHOR; ++horizontal)
        {
            if (ice.land_mask[horizontal] < 0.5f)
            {
                ice.ocean_heat_accumulator[horizontal] = 0.0f;
                ice.ocean_freshwater_accumulator[horizontal] = 0.0f;
                ice.ocean_stress_x_accumulator[horizontal] = 0.0f;
                ice.ocean_stress_y_accumulator[horizontal] = 0.0f;
                ice.ocean_friction_accumulator[horizontal] = 0.0f;
                ice.ocean_snow_accumulator[horizontal] = 0.0f;
            }
        }
        ice.ocean_accumulation_count = 0;
    }

    for (horizontal = 0; horizontal < PLASIC_NHOR; ++horizontal)
    {
        if (ice.land_mask[horizontal] < 0.5f)
        {
            ice.snow_change[horizontal] +=
                (ice.snow[horizontal] - work->horizontal_a[horizontal]) /
                ice.timestep;
        }
    }
    memcpy(sea.coupled_surface_temperature, ice.surface_temperature,
           sizeof(sea.coupled_surface_temperature));
    memcpy(sea.coupled_ice_cover, ice.ice_cover,
           sizeof(sea.coupled_ice_cover));
    memcpy(sea.coupled_ice_thickness, ice.ice_thickness,
           sizeof(sea.coupled_ice_thickness));
    memcpy(sea.coupled_snow, ice.snow, sizeof(sea.coupled_snow));
    memcpy(sea.coupled_snow_melt, ice.snow_melt,
           sizeof(sea.coupled_snow_melt));
    memcpy(sea.coupled_snow_change, ice.snow_change,
           sizeof(sea.coupled_snow_change));
    memcpy(sea.coupled_sst, ice.sst, sizeof(sea.coupled_sst));
    memcpy(sea.coupled_mixed_layer_depth, ice.mixed_layer_depth,
           sizeof(sea.coupled_mixed_layer_depth));
}

static void sea_step(void)
{
    const int32_t surface_level = PLASIC_NLEP - 1;
    const int coupling_due = sea_prepare_step(
        PLASIC_NHOR, surface_level, plasic_nstep,
        sea.coupling_interval, plasic_nkits,
        &sea.accumulation_count, plasic_dls, plasic_dshfl,
        plasic_dswfl, plasic_dlwfl, plasic_dlhfl,
        plasic_dprl, plasic_dprc, plasic_devap,
        plasic_dprs, plasic_dtaux, plasic_dtauy,
        plasic_dust3, plasic_dshdt,
        plasic_dlhdt, sea.heat_accumulator,
        sea.freshwater_accumulator, sea.snowfall_accumulator,
        sea.stress_x_accumulator, sea.stress_y_accumulator,
        sea.friction_accumulator,
        sea.sensible_heat_accumulator,
        sea.sensible_derivative_accumulator,
        sea.latent_heat_accumulator,
        sea.latent_derivative_accumulator,
        sea.shortwave_accumulator, sea.longwave_accumulator);

    if (coupling_due != 0)
    {
        /*
         * Sea ice integrates the whole coupling window in one step.  The
         * window is the number of atmospheric steps whose fluxes were just
         * averaged: normally sea.coupling_interval, but the first window
         * after a start or restart can be shorter or longer.
         */
        ice.timestep = (float)sea.accumulation_count * plasic_deltsec;
        ice_step();
        sea_finish_coupling(
            PLASIC_NHOR, &sea.accumulation_count, plasic_dls,
            sea.coupled_surface_temperature, sea.coupled_sst,
            sea.coupled_mixed_layer_depth, sea.coupled_ice_cover,
            sea.coupled_ice_thickness, sea.coupled_snow,
            sea.coupled_snow_melt, sea.coupled_snow_change,
            sea.surface_temperature, sea.sst, plasic_dmld,
            plasic_dicec, plasic_diced, plasic_dsnow,
            plasic_dsmelt, plasic_dsndch,
            sea.heat_accumulator, sea.freshwater_accumulator,
            sea.snowfall_accumulator, sea.stress_x_accumulator,
            sea.stress_y_accumulator,
            sea.friction_accumulator, sea.sensible_heat_accumulator,
            sea.sensible_derivative_accumulator,
            sea.latent_heat_accumulator,
            sea.latent_derivative_accumulator,
            sea.shortwave_accumulator, sea.longwave_accumulator);
    }
    sea_update_surface(
        PLASIC_NHOR, surface_level, plasic_rdbrv,
        plasic_ra1, plasic_ra2, plasic_ra4,
        plasic_tmelt, sea.charnock, plasic_gascon,
        plasic_ga, sea.sea_albedo, sea.ice_albedo,
        sea.sea_roughness, sea.ice_roughness,
        sea.sea_relative_humidity, sea.ice_relative_humidity,
        plasic_dls, sea.surface_temperature, plasic_dicec,
        plasic_dt, plasic_dq, plasic_dp,
        plasic_dtaux, plasic_dtauy, sea.surface_humidity,
        plasic_dwetfac, plasic_dalb, plasic_dz0);
}

int physics_step(void)
{
    int status;

    status = moisture_fixer();
    if (status != PLASIC_RUNTIME_OK)
    {
        return status;
    }
    surface_fluxes();
    vertical_diffusion_step();
    if (plasic_nrad > 0)
    {
        status = radiation_step();
        if (status != PLASIC_RUNTIME_OK)
        {
            return status;
        }
    }
    status = rain_step();
    if (status != PLASIC_RUNTIME_OK)
    {
        return status;
    }
    status = land_step();
    if (status != PLASIC_RUNTIME_OK)
    {
        return status;
    }
    sea_step();
    return PLASIC_RUNTIME_OK;
}
