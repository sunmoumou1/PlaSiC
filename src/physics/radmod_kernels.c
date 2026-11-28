#include "radmod_kernels.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#define RAD_AT(array, horizontal, level, nhor) \
    ((array)[(size_t)(horizontal) + (size_t)(nhor) * (size_t)(level)])

static float rad_min(float left, float right)
{
    return left < right ? left : right;
}

static float rad_max(float left, float right)
{
    return left > right ? left : right;
}

static float rad_square(float value)
{
    return value * value;
}

static float rad_cube(float value)
{
    return (value * value) * value;
}

static float *rad_allocate(size_t count)
{
    return (float *)calloc(count, sizeof(float));
}

static void rad_free_all(float **arrays, size_t count)
{
    size_t index;

    for (index = 0; index < count; ++index)
    {
        free(arrays[index]);
    }
}

void rad_zero_fluxes(
    int32_t nhor, int32_t nlep, float *dfu, float *dfd, float *dftu,
    float *dftd, float *dswfl, float *dlwfl)
{
    const size_t count = (size_t)nhor * (size_t)nlep;
    size_t index;

    for (index = 0; index < count; ++index)
    {
        dfu[index] = 0.0f;
        dfd[index] = 0.0f;
        dftu[index] = 0.0f;
        dftd[index] = 0.0f;
        dswfl[index] = 0.0f;
        dlwfl[index] = 0.0f;
    }
}

void rad_compute_tendencies(
    int32_t nhor, int32_t nlev, float gravity, float dry_air_heat_capacity,
    float cpv_cpd_minus1, const float *dsigma, const float *surface_pressure,
    const float *humidity, const float *shortwave_flux,
    const float *longwave_flux, float *total_flux,
    float *temperature_tendency, float *shortwave_tendency,
    float *longwave_tendency, float *radiative_tendency)
{
    int32_t level;
    int32_t horizontal;
    const int32_t nlep = nlev + 1;

    for (level = 0; level < nlep; ++level)
    {
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t index =
                (size_t)horizontal + (size_t)nhor * (size_t)level;
            total_flux[index] = longwave_flux[index] + shortwave_flux[index];
        }
    }

    for (level = 0; level < nlev; ++level)
    {
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t current =
                (size_t)horizontal + (size_t)nhor * (size_t)level;
            const size_t below = current + (size_t)nhor;
            const float denominator =
                ((dsigma[level] * surface_pressure[horizontal]) *
                 dry_air_heat_capacity) *
                (1.0f + cpv_cpd_minus1 * humidity[current]);
            const float total =
                (-gravity * (total_flux[below] - total_flux[current])) /
                denominator;
            const float shortwave =
                (-gravity *
                 (shortwave_flux[below] - shortwave_flux[current])) /
                denominator;
            const float longwave =
                (-gravity *
                 (longwave_flux[below] - longwave_flux[current])) /
                denominator;

            radiative_tendency[current] = total;
            temperature_tendency[current] =
                temperature_tendency[current] + total;
            shortwave_tendency[current] = shortwave;
            longwave_tendency[current] = longwave;
        }
    }
}

void rad_solar_angles(
    int32_t nlon, int32_t nlat, int32_t minute_of_day,
    float pi, float two_pi, float declination, float distance_factor,
    const double *sin_latitude, const float *cos_latitude,
    float *cos_zenith, float *stored_distance_factor)
{
    const float longitude_scale = two_pi / (float)nlon;
    const float time_scale = two_pi / 1440.0f;
    const float sin_declination = sinf(declination);
    const float cos_declination = cosf(declination);
    int32_t latitude;
    int32_t longitude;

    for (latitude = 0; latitude < nlat; ++latitude)
    {
        for (longitude = 0; longitude < nlon; ++longitude)
        {
            const size_t horizontal =
                (size_t)longitude + (size_t)nlon * (size_t)latitude;
            const float hour_angle =
                (float)minute_of_day * time_scale +
                (float)longitude * longitude_scale - pi;
            const float second_term =
                (cos_latitude[latitude] * cos_declination) *
                cosf(hour_angle);
            const float zenith =
                (float)((double)sin_declination * sin_latitude[latitude] +
                        (double)second_term);

            cos_zenith[horizontal] = zenith > 0.0f ? zenith : 0.0f;
        }
    }
    *stored_distance_factor = distance_factor;
}

int32_t rad_synthetic_ozone(
    int32_t nlon, int32_t nlat, int32_t nlev, float calendar_day,
    int32_t days_per_year, float two_pi,
    float a0, float a1, float seasonal_amplitude, float profile_height,
    float profile_thickness, float seasonal_offset,
    float gas_constant, float gravity, const double *sin_latitude,
    const float *sigma_half, const float *dsigma,
    const float *surface_pressure, const float *temperature,
    float *ozone)
{
    const int32_t nhor = nlon * nlat;
    const float ozone_density = 2.14f;
    const float ozone_factor = 100.0f / ozone_density;
    const float seasonal_cosine = cosf(
        (two_pi * (calendar_day - seasonal_offset)) /
        (float)days_per_year);
    const float profile_constant =
        expf(-profile_height / profile_thickness);
    float *column_amount = rad_allocate((size_t)nhor);
    float *height = rad_allocate((size_t)nhor);
    float *remaining = rad_allocate((size_t)nhor);
    int32_t latitude;
    int32_t longitude;
    int32_t level;

    if (column_amount == NULL || height == NULL || remaining == NULL)
    {
        float *arrays[] = {column_amount, height, remaining};
        rad_free_all(arrays, sizeof(arrays) / sizeof(arrays[0]));
        return 0;
    }

    for (latitude = 0; latitude < nlat; ++latitude)
    {
        const double latitude_sine = sin_latitude[latitude];
        const float amount = (float)((double)a0 + (double)a1 * fabs(latitude_sine) +
                                     ((double)seasonal_amplitude * latitude_sine) *
                                         (double)seasonal_cosine);

        for (longitude = 0; longitude < nlon; ++longitude)
        {
            const size_t horizontal =
                (size_t)longitude + (size_t)nlon * (size_t)latitude;
            column_amount[horizontal] = amount;
            remaining[horizontal] = amount;
            height[horizontal] = 0.0f;
        }
    }

    for (level = nlev - 1; level >= 1; --level)
    {
        int32_t horizontal;
        const float sigma_ratio_log =
            logf(sigma_half[level - 1] / sigma_half[level]);

        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t index =
                (size_t)horizontal + (size_t)nhor * (size_t)level;
            float layer_amount;

            height[horizontal] = height[horizontal] -
                                 ((temperature[index] * gas_constant) / gravity) *
                                     sigma_ratio_log;
            layer_amount =
                -(column_amount[horizontal] +
                  column_amount[horizontal] * profile_constant) /
                    (1.0f + expf((height[horizontal] - profile_height) /
                                 profile_thickness)) +
                remaining[horizontal];
            ozone[index] =
                (layer_amount * gravity) /
                ((ozone_factor * dsigma[level]) *
                 surface_pressure[horizontal]);
            remaining[horizontal] = remaining[horizontal] - layer_amount;
        }
    }

    for (longitude = 0; longitude < nhor; ++longitude)
    {
        ozone[longitude] =
            (remaining[longitude] * gravity) /
            ((ozone_factor * dsigma[0]) * surface_pressure[longitude]);
    }

    {
        float *arrays[] = {column_amount, height, remaining};
        rad_free_all(arrays, sizeof(arrays) / sizeof(arrays[0]));
    }
    return 1;
}

static float rad_ozone_absorption(float amount)
{
    const float first =
        (0.02118f * amount) /
        (1.0f + 0.042f * amount + 0.000323f * rad_square(amount));
    const float second =
        (1.082f * amount) /
        powf(1.0f + 138.6f * amount, 0.805f);
    const float third =
        (0.0658f * amount) /
        (1.0f + rad_cube(103.6f * amount));

    return first + second + third;
}

static float rad_water_absorption(float amount)
{
    return (2.9f * amount) /
           (powf(1.0f + 141.5f * amount, 0.635f) + 5.925f * amount);
}

int32_t rad_shortwave(
    int32_t nhor, int32_t nlev, int32_t rayleigh_scattering,
    float gravity, float solar_constant, float distance_factor,
    float cloud_tuning_1, float cloud_tuning_2, float cloud_tuning_3,
    const float *sigma, const float *dsigma,
    const float *surface_pressure, float *surface_albedo,
    const float *land_mask, const float *sea_ice,
    const float *humidity, const float *cloud_liquid,
    const float *cloud_cover, const float *ozone,
    const float *cos_zenith, const float *temperature,
    float *upward_flux, float *downward_flux, float *net_flux)
{
    const float zero = 1.0e-6f;
    const float solar_band_1 = 0.517f;
    const float solar_band_2 = 0.483f;
    const float water_magnification = 1.66f;
    const float ozone_magnification = 1.9f;
    const float ozone_factor = 100.0f / 2.14f;
    const int32_t nlep = nlev + 1;
    const size_t interface_count = (size_t)nhor * (size_t)nlep;
    const size_t layer_count = (size_t)nhor * (size_t)nlev;
    const size_t vector_count = (size_t)nhor;
    const size_t workspace_count =
        8U * interface_count + 20U * layer_count + 42U * vector_count;
    float *workspace = rad_allocate(workspace_count);
    int32_t *has_sun =
        (int32_t *)calloc((size_t)nhor, sizeof(int32_t));
    float *cursor;
    float *zt1;
    float *zt2;
    float *zr1s;
    float *zr2s;
    float *zrl1;
    float *zrl2;
    float *zrl1s;
    float *zrl2s;
    float *ztb1;
    float *ztb2;
    float *ztb1u;
    float *ztb2u;
    float *zrb1;
    float *zrb2;
    float *zrb1s;
    float *zrb2s;
    float *zo3l;
    float *zxo3l;
    float *zwvl;
    float *zywvl;
    float *zrcs;
    float *zrcsu;
    float *zrcl1;
    float *zrcl2;
    float *zrcl1s;
    float *zrcl2s;
    float *ztcl2;
    float *ztcl2s;
    float *zftop1;
    float *zftop2;
    float *zcs;
    float *zm;
    float *zo3;
    float *zo3t;
    float *zxo3t;
    float *zto3;
    float *zto3u;
    float *zto3t;
    float *zto3tu;
    float *zwv;
    float *zwvt;
    float *zywvt;
    float *ztwv;
    float *ztwvu;
    float *ztwvt;
    float *ztwvtu;
    float *zra1;
    float *zra2;
    float *zra1s;
    float *zra2s;
    float *zta1;
    float *zta2;
    float *zta1s;
    float *zta2s;
    float *z1mrabr;
    float *zlwp;
    float *ztau;
    float *zlog;
    float *zb2;
    float *zom0;
    float *zuz;
    float *zun;
    float *zr;
    float *zexp;
    float *zu;
    float *zb1;
    float *zfu1;
    float *zfu2;
    float *zfd1;
    float *zfd2;
    int32_t horizontal;
    int32_t level;

    if (workspace == NULL || has_sun == NULL)
    {
        free(workspace);
        free(has_sun);
        return 0;
    }

    cursor = workspace;
#define RAD_TAKE_INTERFACE(name)   \
    do                             \
    {                              \
        name = cursor;             \
        cursor += interface_count; \
    } while (0)
#define RAD_TAKE_LAYER(name)   \
    do                         \
    {                          \
        name = cursor;         \
        cursor += layer_count; \
    } while (0)
#define RAD_TAKE_VECTOR(name)   \
    do                          \
    {                           \
        name = cursor;          \
        cursor += vector_count; \
    } while (0)
    RAD_TAKE_INTERFACE(zt1);
    RAD_TAKE_INTERFACE(zt2);
    RAD_TAKE_INTERFACE(zr1s);
    RAD_TAKE_INTERFACE(zr2s);
    RAD_TAKE_INTERFACE(zrl1);
    RAD_TAKE_INTERFACE(zrl2);
    RAD_TAKE_INTERFACE(zrl1s);
    RAD_TAKE_INTERFACE(zrl2s);
    RAD_TAKE_LAYER(ztb1);
    RAD_TAKE_LAYER(ztb2);
    RAD_TAKE_LAYER(ztb1u);
    RAD_TAKE_LAYER(ztb2u);
    RAD_TAKE_LAYER(zrb1);
    RAD_TAKE_LAYER(zrb2);
    RAD_TAKE_LAYER(zrb1s);
    RAD_TAKE_LAYER(zrb2s);
    RAD_TAKE_LAYER(zo3l);
    RAD_TAKE_LAYER(zxo3l);
    RAD_TAKE_LAYER(zwvl);
    RAD_TAKE_LAYER(zywvl);
    RAD_TAKE_LAYER(zrcs);
    RAD_TAKE_LAYER(zrcsu);
    RAD_TAKE_LAYER(zrcl1);
    RAD_TAKE_LAYER(zrcl2);
    RAD_TAKE_LAYER(zrcl1s);
    RAD_TAKE_LAYER(zrcl2s);
    RAD_TAKE_LAYER(ztcl2);
    RAD_TAKE_LAYER(ztcl2s);
    RAD_TAKE_VECTOR(zftop1);
    RAD_TAKE_VECTOR(zftop2);
    RAD_TAKE_VECTOR(zcs);
    RAD_TAKE_VECTOR(zm);
    RAD_TAKE_VECTOR(zo3);
    RAD_TAKE_VECTOR(zo3t);
    RAD_TAKE_VECTOR(zxo3t);
    RAD_TAKE_VECTOR(zto3);
    RAD_TAKE_VECTOR(zto3u);
    RAD_TAKE_VECTOR(zto3t);
    RAD_TAKE_VECTOR(zto3tu);
    RAD_TAKE_VECTOR(zwv);
    RAD_TAKE_VECTOR(zwvt);
    RAD_TAKE_VECTOR(zywvt);
    RAD_TAKE_VECTOR(ztwv);
    RAD_TAKE_VECTOR(ztwvu);
    RAD_TAKE_VECTOR(ztwvt);
    RAD_TAKE_VECTOR(ztwvtu);
    RAD_TAKE_VECTOR(zra1);
    RAD_TAKE_VECTOR(zra2);
    RAD_TAKE_VECTOR(zra1s);
    RAD_TAKE_VECTOR(zra2s);
    RAD_TAKE_VECTOR(zta1);
    RAD_TAKE_VECTOR(zta2);
    RAD_TAKE_VECTOR(zta1s);
    RAD_TAKE_VECTOR(zta2s);
    RAD_TAKE_VECTOR(z1mrabr);
    RAD_TAKE_VECTOR(zlwp);
    RAD_TAKE_VECTOR(ztau);
    RAD_TAKE_VECTOR(zlog);
    RAD_TAKE_VECTOR(zb2);
    RAD_TAKE_VECTOR(zom0);
    RAD_TAKE_VECTOR(zuz);
    RAD_TAKE_VECTOR(zun);
    RAD_TAKE_VECTOR(zr);
    RAD_TAKE_VECTOR(zexp);
    RAD_TAKE_VECTOR(zu);
    RAD_TAKE_VECTOR(zb1);
    RAD_TAKE_VECTOR(zfu1);
    RAD_TAKE_VECTOR(zfu2);
    RAD_TAKE_VECTOR(zfd1);
    RAD_TAKE_VECTOR(zfd2);
#undef RAD_TAKE_INTERFACE
#undef RAD_TAKE_LAYER
#undef RAD_TAKE_VECTOR

    for (horizontal = 0; horizontal < nhor; ++horizontal)
    {
        zftop1[horizontal] =
            ((solar_band_1 * solar_constant) * distance_factor) *
            cos_zenith[horizontal];
        zftop2[horizontal] =
            ((solar_band_2 * solar_constant) * distance_factor) *
            cos_zenith[horizontal];
        has_sun[horizontal] =
            zftop1[horizontal] + zftop2[horizontal] > zero;
        zcs[horizontal] = 1.0f;
    }

    {
        const float mean_cosine = 0.5f;
        const float zb3 =
            (cloud_tuning_1 * sqrtf(mean_cosine)) / mean_cosine;
        const float zb4 = cloud_tuning_2 * sqrtf(mean_cosine);
        const float zb5 =
            (cloud_tuning_3 * mean_cosine) * mean_cosine;
        size_t index;

        for (index = 0; index < layer_count; ++index)
        {
            zrcl1[index] = 0.0f;
            zrcl2[index] = 0.0f;
            ztcl2[index] = 1.0f;
            zrcl1s[index] = 0.0f;
            zrcl2s[index] = 0.0f;
            ztcl2s[index] = 1.0f;
        }
        for (level = 0; level < nlev; ++level)
        {
            for (horizontal = 0; horizontal < nhor; ++horizontal)
            {
                const size_t layer =
                    (size_t)horizontal + (size_t)nhor * (size_t)level;

                if (has_sun[horizontal] && cloud_cover[layer] > 0.0f)
                {
                    float exponent_argument;

                    zlwp[horizontal] = rad_min(
                        1000.0f,
                        (((1000.0f * cloud_liquid[layer]) *
                          surface_pressure[horizontal]) /
                         gravity) *
                            dsigma[level]);
                    ztau[horizontal] =
                        2.0f * powf(
                                   log10f(zlwp[horizontal] + 1.5f), 3.9f);
                    zlog[horizontal] =
                        logf(1000.0f / ztau[horizontal]);
                    zb2[horizontal] =
                        zb4 / logf(3.0f + 0.1f * ztau[horizontal]);
                    zom0[horizontal] = rad_min(
                        0.9999f,
                        1.0f - zb5 * zlog[horizontal]);
                    zun[horizontal] = 1.0f - zom0[horizontal];
                    zuz[horizontal] = zun[horizontal] +
                                      (2.0f * zb2[horizontal]) * zom0[horizontal];
                    zu[horizontal] =
                        sqrtf(zuz[horizontal] / zun[horizontal]);
                    exponent_argument =
                        (ztau[horizontal] *
                         sqrtf(zuz[horizontal] * zun[horizontal])) /
                        mean_cosine;
                    zexp[horizontal] =
                        expf(rad_min(25.0f, exponent_argument));
                    zr[horizontal] =
                        rad_square(zu[horizontal] + 1.0f) *
                            zexp[horizontal] -
                        rad_square(zu[horizontal] - 1.0f) /
                            zexp[horizontal];
                    zrcl1s[layer] = 1.0f -
                                    1.0f /
                                        (1.0f + zb3 * ztau[horizontal]);
                    ztcl2s[layer] =
                        (4.0f * zu[horizontal]) / zr[horizontal];
                    zrcl2s[layer] =
                        ((rad_square(zu[horizontal]) - 1.0f) /
                         zr[horizontal]) *
                        (zexp[horizontal] -
                         1.0f / zexp[horizontal]);

                    zb1[horizontal] =
                        cloud_tuning_1 * sqrtf(cos_zenith[horizontal]);
                    zb2[horizontal] =
                        (cloud_tuning_2 * sqrtf(cos_zenith[horizontal])) /
                        logf(3.0f + 0.1f * ztau[horizontal]);
                    zom0[horizontal] = rad_min(
                        0.9999f,
                        1.0f -
                            ((cloud_tuning_3 * cos_zenith[horizontal]) *
                             cos_zenith[horizontal]) *
                                zlog[horizontal]);
                    zun[horizontal] = 1.0f - zom0[horizontal];
                    zuz[horizontal] = zun[horizontal] +
                                      (2.0f * zb2[horizontal]) * zom0[horizontal];
                    zu[horizontal] =
                        sqrtf(zuz[horizontal] / zun[horizontal]);
                    exponent_argument =
                        (ztau[horizontal] *
                         sqrtf(zuz[horizontal] * zun[horizontal])) /
                        cos_zenith[horizontal];
                    zexp[horizontal] =
                        expf(rad_min(25.0f, exponent_argument));
                    zr[horizontal] =
                        rad_square(zu[horizontal] + 1.0f) *
                            zexp[horizontal] -
                        rad_square(zu[horizontal] - 1.0f) /
                            zexp[horizontal];
                    zrcl1[layer] = 1.0f -
                                   1.0f /
                                       (1.0f +
                                        (zb1[horizontal] * ztau[horizontal]) /
                                            cos_zenith[horizontal]);
                    ztcl2[layer] =
                        (4.0f * zu[horizontal]) / zr[horizontal];
                    zrcl2[layer] =
                        ((rad_square(zu[horizontal]) - 1.0f) /
                         zr[horizontal]) *
                        (zexp[horizontal] -
                         1.0f / zexp[horizontal]);
                    zrcl1[layer] =
                        zcs[horizontal] * zrcl1[layer] +
                        (1.0f - zcs[horizontal]) * zrcl1s[layer];
                    ztcl2[layer] =
                        zcs[horizontal] * ztcl2[layer] +
                        (1.0f - zcs[horizontal]) * ztcl2s[layer];
                    zrcl2[layer] =
                        zcs[horizontal] * zrcl2[layer] +
                        (1.0f - zcs[horizontal]) * zrcl2s[layer];
                }
                zcs[horizontal] = zcs[horizontal] *
                                  (1.0f - cloud_cover[layer]);
            }
        }
    }

    for (horizontal = 0; horizontal < nhor; ++horizontal)
    {
        if (has_sun[horizontal])
        {
            zm[horizontal] = 35.0f /
                             sqrtf(1.0f +
                                   (1224.0f * cos_zenith[horizontal]) *
                                       cos_zenith[horizontal]);
            zcs[horizontal] = 1.0f;
            zo3t[horizontal] = 0.0f;
            zxo3t[horizontal] = 0.0f;
            zwvt[horizontal] = 0.0f;
            zywvt[horizontal] = 0.0f;
        }
    }
    for (level = 0; level < nlev; ++level)
    {
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t layer =
                (size_t)horizontal + (size_t)nhor * (size_t)level;

            if (has_sun[horizontal])
            {
                float water_amount;

                zo3[horizontal] =
                    (((ozone_factor * dsigma[level]) *
                      surface_pressure[horizontal]) *
                     ozone[layer]) /
                    gravity;
                zo3t[horizontal] =
                    zo3t[horizontal] + zo3[horizontal];
                zxo3t[horizontal] =
                    zcs[horizontal] *
                        (zxo3t[horizontal] +
                         zm[horizontal] * zo3[horizontal]) +
                    (1.0f - zcs[horizontal]) *
                        (zxo3t[horizontal] +
                         ozone_magnification * zo3[horizontal]);
                zo3l[layer] = zo3t[horizontal];
                zxo3l[layer] = zxo3t[horizontal];
                water_amount = 0.1f * dsigma[level];
                water_amount = water_amount * humidity[layer];
                water_amount = water_amount * surface_pressure[horizontal];
                water_amount = water_amount / gravity;
                water_amount =
                    water_amount * sqrtf(273.0f / temperature[layer]);
                water_amount = water_amount * sigma[level];
                water_amount = water_amount * surface_pressure[horizontal];
                zwv[horizontal] = water_amount / 100000.0f;
                zwvt[horizontal] =
                    zwvt[horizontal] + zwv[horizontal];
                zywvt[horizontal] =
                    zcs[horizontal] *
                        (zywvt[horizontal] +
                         zm[horizontal] * zwv[horizontal]) +
                    (1.0f - zcs[horizontal]) *
                        (zywvt[horizontal] +
                         water_magnification * zwv[horizontal]);
                zwvl[layer] = zwvt[horizontal];
                zywvl[layer] = zywvt[horizontal];
                zcs[horizontal] = zcs[horizontal] *
                                  (1.0f - cloud_cover[layer]);
                zrcs[layer] = 0.0f;
                zrcsu[layer] = 0.0f;
            }
        }
    }

    for (horizontal = 0; horizontal < nhor; ++horizontal)
    {
        if (has_sun[horizontal])
        {
            const size_t lowest =
                (size_t)horizontal +
                (size_t)nhor * (size_t)(nlev - 1);

            zta1[horizontal] = 1.0f;
            zta1s[horizontal] = 1.0f;
            zra1[horizontal] = 0.0f;
            zra1s[horizontal] = 0.0f;
            zta2[horizontal] = 1.0f;
            zta2s[horizontal] = 1.0f;
            zra2[horizontal] = 0.0f;
            zra2s[horizontal] = 0.0f;
            zto3t[horizontal] = 1.0f;
            zo3[horizontal] = zxo3t[horizontal] +
                              ozone_magnification * zo3t[horizontal];
            zto3tu[horizontal] = 1.0f -
                                 rad_ozone_absorption(zo3[horizontal]) / solar_band_1;
            ztwvt[horizontal] = 1.0f;
            zwv[horizontal] = zywvt[horizontal] +
                              water_magnification * zwvt[horizontal];
            ztwvtu[horizontal] = 1.0f -
                                 rad_water_absorption(zwv[horizontal]) / solar_band_2;
            zrcs[lowest] =
                (0.219f / (1.0f + 0.816f * cos_zenith[horizontal]) *
                     zcs[horizontal] +
                 0.144f *
                     (1.0f - zcs[horizontal] - cloud_cover[lowest])) *
                (float)rayleigh_scattering;
            zrcsu[lowest] =
                0.144f * (1.0f - cloud_cover[lowest]) *
                (float)rayleigh_scattering;
        }
    }

    for (level = 0; level < nlev; ++level)
    {
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t layer =
                (size_t)horizontal + (size_t)nhor * (size_t)level;

            if (has_sun[horizontal])
            {
                zt1[layer] = zta1[horizontal];
                zt2[layer] = zta2[horizontal];
                zr1s[layer] = zra1s[horizontal];
                zr2s[layer] = zra2s[horizontal];
                zrb1[layer] = zrcs[layer] +
                              zrcl1[layer] * cloud_cover[layer];
                zrb1s[layer] = zrcsu[layer] +
                               zrcl1s[layer] * cloud_cover[layer];

                zo3[horizontal] = zxo3l[layer];
                zto3[horizontal] =
                    (1.0f -
                     rad_ozone_absorption(zo3[horizontal]) /
                         solar_band_1) /
                    zto3t[horizontal];
                zto3t[horizontal] =
                    zto3t[horizontal] * zto3[horizontal];
                zo3[horizontal] = zxo3t[horizontal] +
                                  ozone_magnification *
                                      (zo3t[horizontal] - zo3l[layer]);
                zto3u[horizontal] =
                    zto3tu[horizontal] /
                    (1.0f -
                     rad_ozone_absorption(zo3[horizontal]) /
                         solar_band_1);
                zto3tu[horizontal] =
                    zto3tu[horizontal] / zto3u[horizontal];
                ztb1[layer] =
                    1.0f -
                    (1.0f - zto3[horizontal]) *
                        (1.0f - cloud_cover[layer]) -
                    zrb1[layer];
                ztb1u[layer] =
                    1.0f -
                    (1.0f - zto3u[horizontal]) *
                        (1.0f - cloud_cover[layer]) -
                    zrb1s[layer];
                z1mrabr[horizontal] =
                    1.0f /
                    (1.0f - zra1s[horizontal] * zrb1s[layer]);
                zra1[horizontal] = zra1[horizontal] +
                                   ((zta1[horizontal] * zrb1[layer]) *
                                    zta1s[horizontal]) *
                                       z1mrabr[horizontal];
                zta1[horizontal] =
                    (zta1[horizontal] * ztb1[layer]) *
                    z1mrabr[horizontal];
                zra1s[horizontal] = zrb1s[layer] +
                                    ((ztb1u[layer] * zra1s[horizontal]) * ztb1[layer]) *
                                        z1mrabr[horizontal];
                zta1s[horizontal] =
                    (ztb1u[layer] * zta1s[horizontal]) *
                    z1mrabr[horizontal];

                zrb2[layer] = zrcl2[layer] * cloud_cover[layer];
                zrb2s[layer] = zrcl2s[layer] * cloud_cover[layer];
                zwv[horizontal] = zywvl[layer];
                ztwv[horizontal] =
                    (1.0f -
                     rad_water_absorption(zwv[horizontal]) /
                         solar_band_2) /
                    ztwvt[horizontal];
                ztwvt[horizontal] =
                    ztwvt[horizontal] * ztwv[horizontal];
                zwv[horizontal] = zywvt[horizontal] +
                                  water_magnification *
                                      (zwvt[horizontal] - zwvl[layer]);
                ztwvu[horizontal] =
                    ztwvtu[horizontal] /
                    (1.0f -
                     rad_water_absorption(zwv[horizontal]) /
                         solar_band_2);
                ztwvtu[horizontal] =
                    ztwvtu[horizontal] / ztwvu[horizontal];
                ztb2[layer] =
                    1.0f -
                    (1.0f - ztwv[horizontal]) *
                        (1.0f - cloud_cover[layer]) -
                    (1.0f - ztcl2[layer]) * cloud_cover[layer];
                ztb2u[layer] =
                    1.0f -
                    (1.0f - ztwvu[horizontal]) *
                        (1.0f - cloud_cover[layer]) -
                    (1.0f - ztcl2s[layer]) * cloud_cover[layer];
                z1mrabr[horizontal] =
                    1.0f /
                    (1.0f - zra2s[horizontal] * zrb2s[layer]);
                zra2[horizontal] = zra2[horizontal] +
                                   ((zta2[horizontal] * zrb2[layer]) *
                                    zta2s[horizontal]) *
                                       z1mrabr[horizontal];
                zta2[horizontal] =
                    (zta2[horizontal] * ztb2[layer]) *
                    z1mrabr[horizontal];
                zra2s[horizontal] = zrb2s[layer] +
                                    ((ztb2u[layer] * zra2s[horizontal]) * ztb2[layer]) *
                                        z1mrabr[horizontal];
                zta2s[horizontal] =
                    (ztb2u[layer] * zta2s[horizontal]) *
                    z1mrabr[horizontal];
            }
        }
    }

    for (horizontal = 0; horizontal < nhor; ++horizontal)
    {
        if (has_sun[horizontal])
        {
            const size_t surface =
                (size_t)horizontal +
                (size_t)nhor * (size_t)nlev;
            const float direct_ocean_albedo = rad_min(
                0.05f / (cos_zenith[horizontal] + 0.15f), 0.15f);

            zt1[surface] = zta1[horizontal];
            zt2[surface] = zta2[horizontal];
            zr1s[surface] = zra1s[horizontal];
            zr2s[surface] = zra2s[horizontal];
            zra1s[horizontal] = surface_albedo[horizontal];
            zra2s[horizontal] = surface_albedo[horizontal];
            surface_albedo[horizontal] =
                land_mask[horizontal] * surface_albedo[horizontal] +
                (1.0f - land_mask[horizontal]) * sea_ice[horizontal] *
                    surface_albedo[horizontal] +
                (1.0f - land_mask[horizontal]) *
                    (1.0f - sea_ice[horizontal]) * direct_ocean_albedo;
            zra1[horizontal] = surface_albedo[horizontal];
            zra2[horizontal] = surface_albedo[horizontal];
        }
    }

    for (level = nlev - 1; level >= 0; --level)
    {
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t layer =
                (size_t)horizontal + (size_t)nhor * (size_t)level;
            const size_t below = layer + (size_t)nhor;

            if (has_sun[horizontal])
            {
                const float denominator_1 =
                    1.0f - zra1s[horizontal] * zrb1s[layer];
                const float denominator_2 =
                    1.0f - zra2s[horizontal] * zrb2s[layer];

                zrl1[below] = zra1[horizontal];
                zrl2[below] = zra2[horizontal];
                zrl1s[below] = zra1s[horizontal];
                zrl2s[below] = zra2s[horizontal];
                zra1[horizontal] = zrb1[layer] +
                                   ((ztb1[layer] * zra1[horizontal]) * ztb1u[layer]) /
                                       denominator_1;
                zra1s[horizontal] = zrb1s[layer] +
                                    ((ztb1u[layer] * zra1s[horizontal]) * ztb1u[layer]) /
                                        denominator_1;
                zra2[horizontal] = zrb2[layer] +
                                   ((ztb2[layer] * zra2[horizontal]) * ztb2u[layer]) /
                                       denominator_2;
                zra2s[horizontal] = zrb2s[layer] +
                                    ((ztb2u[layer] * zra2s[horizontal]) * ztb2u[layer]) /
                                        denominator_2;
            }
        }
    }
    for (horizontal = 0; horizontal < nhor; ++horizontal)
    {
        if (has_sun[horizontal])
        {
            zrl1[horizontal] = zra1[horizontal];
            zrl2[horizontal] = zra2[horizontal];
            zrl1s[horizontal] = zra1s[horizontal];
            zrl2s[horizontal] = zra2s[horizontal];
        }
    }

    for (level = 0; level < nlep; ++level)
    {
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t interface =
                (size_t)horizontal + (size_t)nhor * (size_t)level;

            if (has_sun[horizontal])
            {
                z1mrabr[horizontal] =
                    1.0f /
                    (1.0f - zr1s[interface] * zrl1s[interface]);
                zfd1[horizontal] =
                    zt1[interface] * z1mrabr[horizontal];
                zfu1[horizontal] =
                    -zt1[interface] * zrl1[interface] *
                    z1mrabr[horizontal];
                z1mrabr[horizontal] =
                    1.0f /
                    (1.0f - zr2s[interface] * zrl2s[interface]);
                zfd2[horizontal] =
                    zt2[interface] * z1mrabr[horizontal];
                zfu2[horizontal] =
                    -zt2[interface] * zrl2[interface] *
                    z1mrabr[horizontal];
                upward_flux[interface] =
                    zfu1[horizontal] * zftop1[horizontal] +
                    zfu2[horizontal] * zftop2[horizontal];
                downward_flux[interface] =
                    zfd1[horizontal] * zftop1[horizontal] +
                    zfd2[horizontal] * zftop2[horizontal];
                net_flux[interface] =
                    upward_flux[interface] + downward_flux[interface];
            }
        }
    }

    free(workspace);
    free(has_sun);
    return 1;
}

int32_t rad_longwave(
    int32_t nhor, int32_t nlev, float gravity,
    float water_continuum, float cloud_absorption,
    const float *sigma, const float *sigma_half,
    const float *dsigma, const float *surface_pressure,
    const float *land_mask, const float *humidity,
    const float *temperature, const float *cloud_cover,
    const float *cloud_liquid, const float *ozone,
    const float *carbon_dioxide, float *upward_flux,
    float *downward_flux, float *net_flux)
{
    const float stefan_boltzmann = 5.67e-8f;
    const float molecular_weight_air = 0.0289644f;
    const float molecular_weight_co2 = 0.0440098f;
    const float volume_to_mass_co2 =
        molecular_weight_co2 / molecular_weight_air;
    const float co2_density = 1.9635f;
    const float ozone_density = 2.14f;
    const float water_factor = 0.1f;
    const float co2_factor = 100.0f / co2_density;
    const float ozone_factor = 100.0f / ozone_density;
    const float zero = 1.0e-6f;
    const int32_t nlep = nlev + 1;
    const size_t interface_count = (size_t)nhor * (size_t)nlep;
    const size_t boundary_count =
        (size_t)nhor * (size_t)(nlep + 1);
    const size_t layer_count = (size_t)nhor * (size_t)nlev;
    const size_t vector_count = (size_t)nhor;
    const size_t workspace_count =
        2U * boundary_count + 2U * interface_count +
        7U * layer_count + 13U * vector_count;
    float *workspace = rad_allocate(workspace_count);
    float *cursor;
    float *zbu;
    float *zbd;
    float *zst4h;
    float *zst4;
    float *ztau;
    float *zq;
    float *zqo3;
    float *zqco2;
    float *ztausf;
    float *ztaucs;
    float *ztaucc0;
    float *ztaucc;
    float *ztau0;
    float *zsumwv;
    float *zsumo3;
    float *zsumco2;
    float *zsfac;
    float *zah2o;
    float *zaco2;
    float *zao3;
    float *zth2o;
    float *zbdl;
    float *zeps;
    float *zps2;
    float zao30;
    float zco20;
    float zh2o0a;
    float zh2o0;
    float zao3c;
    float zaco2c;
    float zah2oc;
    float zth2oc;
    int32_t horizontal;
    int32_t level;

    if (workspace == NULL)
    {
        return 0;
    }

    cursor = workspace;
    zbu = cursor;
    cursor += boundary_count;
    zbd = cursor;
    cursor += boundary_count;
    zst4h = cursor;
    cursor += interface_count;
    zst4 = cursor;
    cursor += interface_count;
    ztau = cursor;
    cursor += layer_count;
    zq = cursor;
    cursor += layer_count;
    zqo3 = cursor;
    cursor += layer_count;
    zqco2 = cursor;
    cursor += layer_count;
    ztausf = cursor;
    cursor += layer_count;
    ztaucs = cursor;
    cursor += layer_count;
    ztaucc0 = cursor;
    cursor += layer_count;
    ztaucc = cursor;
    cursor += vector_count;
    ztau0 = cursor;
    cursor += vector_count;
    zsumwv = cursor;
    cursor += vector_count;
    zsumo3 = cursor;
    cursor += vector_count;
    zsumco2 = cursor;
    cursor += vector_count;
    zsfac = cursor;
    cursor += vector_count;
    zah2o = cursor;
    cursor += vector_count;
    zaco2 = cursor;
    cursor += vector_count;
    zao3 = cursor;
    cursor += vector_count;
    zth2o = cursor;
    cursor += vector_count;
    zbdl = cursor;
    cursor += vector_count;
    zeps = cursor;
    cursor += vector_count;
    zps2 = cursor;

    zao30 = 0.209f * powf(7.0e-5f, 0.436f);
    zco20 = 0.0676f * powf(0.01022f, 0.421f);
    zh2o0a = 0.846f * powf(3.59e-5f, 0.243f);
    zh2o0 = 0.832f * powf(0.0286f, 0.26f);
    zao3c =
        0.209f * powf(0.01f + 7.0e-5f, 0.436f) - zao30 -
        0.0212f * log10f(0.01f);
    zaco2c = 0.0676f * powf(1.01022f, 0.421f) - zco20;
    zah2oc =
        0.846f * powf(0.01f + 3.59e-5f, 0.243f) - zh2o0a -
        0.24f * log10f(0.02f);
    zth2oc =
        1.0f -
        (0.832f * powf(2.0f + 0.0286f, 0.26f) - zh2o0) +
        0.1196f * logf(2.0f - 0.6931f);

    for (horizontal = 0; horizontal < nhor; ++horizontal)
    {
        zps2[horizontal] =
            surface_pressure[horizontal] * surface_pressure[horizontal];
    }

    for (level = 0; level < nlep; ++level)
    {
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t index =
                (size_t)horizontal + (size_t)nhor * (size_t)level;
            const float temperature_squared =
                temperature[index] * temperature[index];

            upward_flux[index] = 0.0f;
            downward_flux[index] = 0.0f;
            zst4[index] = stefan_boltzmann *
                          (temperature_squared * temperature_squared);
        }
    }

    for (horizontal = 0; horizontal < nhor; ++horizontal)
    {
        zst4h[horizontal] = zst4[horizontal] -
                            ((zst4[horizontal] -
                              zst4[(size_t)horizontal + (size_t)nhor]) *
                             sigma[0]) /
                                (sigma[0] - sigma[1]);
    }
    for (level = 1; level < nlev; ++level)
    {
        const int32_t above_level = level - 1;

        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t current =
                (size_t)horizontal + (size_t)nhor * (size_t)level;
            const size_t above = current - (size_t)nhor;

            zst4h[current] =
                (zst4[current] *
                     (sigma[above_level] - sigma_half[above_level]) +
                 zst4[above] *
                     (sigma_half[above_level] - sigma[level])) /
                (sigma[above_level] - sigma[level]);
        }
    }
    for (horizontal = 0; horizontal < nhor; ++horizontal)
    {
        const size_t lowest =
            (size_t)horizontal + (size_t)nhor * (size_t)(nlev - 1);
        const size_t above = lowest - (size_t)nhor;
        const size_t surface = lowest + (size_t)nhor;

        if ((zst4[lowest] - zst4h[lowest]) *
                (zst4[lowest] - zst4[surface]) >
            0.0f)
        {
            zst4h[surface] = zst4[surface];
        }
        else
        {
            zst4h[surface] = zst4[lowest] +
                             ((zst4[lowest] - zst4[above]) *
                              (1.0f - sigma[nlev - 1])) /
                                 (sigma[nlev - 1] - sigma[nlev - 2]);
        }
        zbd[horizontal] = 0.0f;
        zeps[horizontal] =
            land_mask[horizontal] + 0.98f * (1.0f - land_mask[horizontal]);
        RAD_AT(zbu, horizontal, nlep, nhor) =
            zeps[horizontal] * zst4[surface];
    }

    for (level = 0; level < nlev; ++level)
    {
        float factor_1 = sigma[level] * dsigma[level];
        const float factor_2 = volume_to_mass_co2 * 1.0e-6f;
        float factor_3 = -1.66f * cloud_absorption;

        factor_1 = factor_1 / gravity;
        factor_1 = factor_1 / 100000.0f;
        factor_3 = factor_3 * 1000.0f;
        factor_3 = factor_3 * dsigma[level];
        factor_3 = factor_3 / gravity;
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t layer =
                (size_t)horizontal + (size_t)nhor * (size_t)level;

            zsfac[horizontal] = factor_1 * zps2[horizontal];
            zq[layer] =
                (water_factor * zsfac[horizontal]) * humidity[layer];
            zqo3[layer] =
                (ozone_factor * zsfac[horizontal]) * ozone[layer];
            zqco2[layer] =
                ((co2_factor * factor_2) * zsfac[horizontal]) *
                carbon_dioxide[layer];
            ztaucc0[layer] =
                1.0f - cloud_cover[layer] *
                           (1.0f -
                            expf((factor_3 * cloud_liquid[layer]) *
                                 surface_pressure[horizontal]));
        }
    }

    for (level = 0; level < nlev; ++level)
    {
        int32_t second_level;

        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            ztaucc[horizontal] = 1.0f;
            zsumwv[horizontal] = 0.0f;
            zsumo3[horizontal] = 0.0f;
            zsumco2[horizontal] = 0.0f;
        }

        for (second_level = level;
             second_level < nlev;
             ++second_level)
        {
            for (horizontal = 0; horizontal < nhor; ++horizontal)
            {
                const size_t layer =
                    (size_t)horizontal +
                    (size_t)nhor * (size_t)second_level;
                float clear_transmissivity;

                zsumwv[horizontal] =
                    zsumwv[horizontal] + zq[layer];
                zsumo3[horizontal] =
                    zsumo3[horizontal] + zqo3[layer];
                zsumco2[horizontal] =
                    zsumco2[horizontal] + zqco2[layer];

                if (zsumwv[horizontal] <= 0.01f)
                {
                    zah2o[horizontal] =
                        0.846f *
                            powf(zsumwv[horizontal] + 3.59e-5f, 0.243f) -
                        zh2o0a;
                }
                else
                {
                    zah2o[horizontal] =
                        0.24f * log10f(zsumwv[horizontal] + 0.01f) +
                        zah2oc;
                }
                if (water_continuum > 0.0f)
                {
                    zah2o[horizontal] = rad_min(
                        zah2o[horizontal] +
                            (1.0f -
                             expf(-water_continuum * zsumwv[horizontal])),
                        1.0f);
                }

                if (zsumco2[horizontal] <= 1.0f)
                {
                    zaco2[horizontal] =
                        0.0676f *
                            powf(zsumco2[horizontal] + 0.01022f, 0.421f) -
                        zco20;
                }
                else
                {
                    zaco2[horizontal] =
                        0.0546f * log10f(zsumco2[horizontal]) + zaco2c;
                }

                if (zsumwv[horizontal] <= 2.0f)
                {
                    zth2o[horizontal] = 1.0f -
                                        (0.832f *
                                             powf(zsumwv[horizontal] + 0.0286f, 0.26f) -
                                         zh2o0);
                }
                else
                {
                    zth2o[horizontal] = rad_max(
                        0.0f,
                        zth2oc -
                            0.1196f *
                                logf(zsumwv[horizontal] - 0.6931f));
                }

                if (zsumo3[horizontal] <= 0.01f)
                {
                    zao3[horizontal] =
                        0.209f *
                            powf(zsumo3[horizontal] + 7.0e-5f, 0.436f) -
                        zao30;
                }
                else
                {
                    zao3[horizontal] =
                        0.0212f * log10f(zsumo3[horizontal]) + zao3c;
                }

                clear_transmissivity =
                    1.0f - zah2o[horizontal] - zao3[horizontal] -
                    zaco2[horizontal] * zth2o[horizontal];
                ztaucs[layer] = rad_min(
                    1.0f - zero,
                    rad_max(zero, clear_transmissivity));
                ztaucc[horizontal] =
                    ztaucc[horizontal] * ztaucc0[layer];
                ztau[layer] =
                    ztaucs[layer] * ztaucc[horizontal];
            }
        }

        if (level == 0)
        {
            for (horizontal = 0; horizontal < nhor; ++horizontal)
            {
                ztau0[horizontal] = 1.0f - zero;
            }
            for (second_level = 0;
                 second_level < nlev;
                 ++second_level)
            {
                for (horizontal = 0; horizontal < nhor; ++horizontal)
                {
                    const size_t layer =
                        (size_t)horizontal +
                        (size_t)nhor * (size_t)second_level;
                    const size_t upper_half =
                        (size_t)horizontal +
                        (size_t)nhor * (size_t)second_level;
                    const size_t lower_half = upper_half + (size_t)nhor;
                    const size_t full = layer;

                    ztau0[horizontal] =
                        ((ztaucs[layer] / ztau0[horizontal]) *
                         ztaucc0[layer]);
                    ztau0[horizontal] = rad_min(
                        1.0f - zero,
                        rad_max(zero, ztau0[horizontal]));
                    if ((zst4[full] - zst4h[upper_half]) *
                            (zst4[full] - zst4h[lower_half]) >
                        0.0f)
                    {
                        RAD_AT(zbd, horizontal, second_level + 1, nhor) =
                            0.5f * zst4[full] +
                            0.25f *
                                (zst4h[upper_half] + zst4h[lower_half]);
                        RAD_AT(zbu, horizontal, second_level + 1, nhor) =
                            RAD_AT(
                                zbd, horizontal, second_level + 1, nhor);
                    }
                    else
                    {
                        RAD_AT(zbd, horizontal, second_level + 1, nhor) =
                            (zst4h[lower_half] -
                             ztau0[horizontal] * zst4h[upper_half]) /
                                (1.0f - ztau0[horizontal]) -
                            (zst4h[upper_half] - zst4h[lower_half]) /
                                logf(ztau0[horizontal]);
                        RAD_AT(zbu, horizontal, second_level + 1, nhor) =
                            zst4h[upper_half] + zst4h[lower_half] -
                            RAD_AT(
                                zbd, horizontal, second_level + 1, nhor);
                    }
                    ztau0[horizontal] = ztaucs[layer];
                }
            }
        }

        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t interface =
                (size_t)horizontal + (size_t)nhor * (size_t)level;

            upward_flux[interface] = upward_flux[interface] -
                                     RAD_AT(zbu, horizontal, level + 1, nhor);
            downward_flux[interface] = downward_flux[interface] +
                                       RAD_AT(zbd, horizontal, level, nhor);
            zbdl[horizontal] =
                RAD_AT(zbd, horizontal, level, nhor) -
                RAD_AT(zbd, horizontal, level + 1, nhor);
        }
        for (second_level = level;
             second_level < nlev;
             ++second_level)
        {
            for (horizontal = 0; horizontal < nhor; ++horizontal)
            {
                const size_t layer =
                    (size_t)horizontal +
                    (size_t)nhor * (size_t)second_level;
                const size_t target_interface =
                    (size_t)horizontal +
                    (size_t)nhor * (size_t)(second_level + 1);

                upward_flux[(size_t)horizontal +
                            (size_t)nhor * (size_t)level] =
                    upward_flux[(size_t)horizontal +
                                (size_t)nhor * (size_t)level] -
                    (RAD_AT(zbu, horizontal, second_level + 2, nhor) -
                     RAD_AT(zbu, horizontal, second_level + 1, nhor)) *
                        ztau[layer];
                downward_flux[target_interface] =
                    downward_flux[target_interface] +
                    zbdl[horizontal] * ztau[layer];
            }
        }
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            RAD_AT(ztausf, horizontal, level, nhor) =
                RAD_AT(ztau, horizontal, nlev - 1, nhor);
        }
    }

    for (horizontal = 0; horizontal < nhor; ++horizontal)
    {
        const size_t surface =
            (size_t)horizontal + (size_t)nhor * (size_t)nlev;

        upward_flux[surface] = upward_flux[surface] -
                               RAD_AT(zbu, horizontal, nlep, nhor);
        downward_flux[surface] = downward_flux[surface] +
                                 RAD_AT(zbd, horizontal, nlev, nhor);
        zeps[horizontal] =
            (1.0f - zeps[horizontal]) * downward_flux[surface];
        /* Include the unattenuated reflected flux at the surface. */
        upward_flux[surface] -= zeps[horizontal];
    }
    for (level = 0; level < nlev; ++level)
    {
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t interface =
                (size_t)horizontal + (size_t)nhor * (size_t)level;

            upward_flux[interface] = upward_flux[interface] -
                                     RAD_AT(ztausf, horizontal, level, nhor) *
                                         zeps[horizontal];
        }
    }
    for (level = 0; level < nlep; ++level)
    {
        for (horizontal = 0; horizontal < nhor; ++horizontal)
        {
            const size_t interface =
                (size_t)horizontal + (size_t)nhor * (size_t)level;
            net_flux[interface] =
                upward_flux[interface] + downward_flux[interface];
        }
    }

    free(workspace);
    return 1;
}
