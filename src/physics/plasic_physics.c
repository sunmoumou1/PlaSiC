#include "physics_internal.h"

char physics_error_message[512];
misc_state misc;
flux_state flux;
rain_state rain;
radiation_state radiation;
land_state land;
sea_state sea;
ocean_state ocean;
deep_ocean_state deep_ocean;
ice_state ice;
physics_workspace *work;

void physics_set_error(const char *format, ...)
{
    va_list arguments;

    va_start(arguments, format);
    (void)vsnprintf(
        physics_error_message, sizeof(physics_error_message), format,
        arguments);
    va_end(arguments);
}

const char *plasic_physics_error(void)
{
    return physics_error_message[0] == '\0' ?
               "atmospheric physics initialization failed" :
               physics_error_message;
}

int plasic_physics_initialize(
    const char *restart_path, const char *surface_data_directory,
    int32_t sea_coupling_interval)
{
    return physics_initialize(
        restart_path, surface_data_directory, sea_coupling_interval);
}

void plasic_physics_finalize(void)
{
    physics_finalize();
}

int plasic_physics_step(void)
{
    return physics_step();
}

int plasic_physics_write_restart(void)
{
    return physics_write_restart();
}

int plasic_physics_frame_field(
    int32_t variable_id, const float **values, int32_t *layers)
{
    return physics_frame_field(variable_id, values, layers);
}

int plasic_physics_set_co2_bands(const float *bands, int32_t band_count)
{
    int32_t latitude;

    if (bands == NULL || band_count != 12)
    {
        return PLASIC_RUNTIME_BAD_ARGUMENT;
    }
    for (latitude = 0; latitude < PLASIC_NLPP; ++latitude)
    {
        const double degrees = asin(plasic_sid[latitude]) *
                               180.0 / 3.14159265358979323846;
        int32_t band = (int32_t)floor((degrees + 90.0) / 15.0);
        int32_t level;
        int32_t longitude;

        band = band < 0 ? 0 : (band >= band_count ? band_count - 1 : band);
        for (level = 0; level < PLASIC_NLEV; ++level)
        {
            for (longitude = 0; longitude < PLASIC_NLON; ++longitude)
            {
                const size_t horizontal =
                    (size_t)longitude + (size_t)PLASIC_NLON *
                                            (size_t)latitude;
                radiation.carbon_dioxide[
                    horizontal + (size_t)level * PLASIC_NHOR] = bands[band];
            }
        }
    }
    return PLASIC_RUNTIME_OK;
}

void plasic_physics_interpolate_ozone(const double *current, const double *neighbor, float weight)
{
    for (size_t i = 0; i < (size_t)PLASIC_NHOR * PLASIC_NLEV; ++i)
        radiation.ozone[i] = (float)(current[i] + weight * (neighbor[i] - current[i]));
    radiation.use_synthetic_ozone = 0;
}

const float *plasic_physics_co2_field(void)
{
    return radiation.carbon_dioxide;
}

const float *plasic_physics_ocean_depth(int32_t *count)
{
    if (count != NULL)
    {
        *count = PLASIC_DEEP_NLEV;
    }
    return deep_ocean.layer_center_depth;
}
