#include "physics_internal.h"

/*
 * Write one deep-ocean restart record.  The in-memory field is halo-extended
 * but restart files store only the contiguous local latitude slice, so the
 * interior rows are packed first.
 */
static int write_deep_restart_field(const char *name, const float *extended)
{
    float *local = (float *)malloc(
        (size_t)PLASIC_NHOR * (size_t)PLASIC_DEEP_NLEV * sizeof(float));
    int status;

    if (local == NULL)
    {
        return PLASIC_RUNTIME_IO_FAILED;
    }
    physics_ocean3d_pack_local(extended, local, PLASIC_DEEP_NLEV);
    status = restart_write_distributed_grid(
        name, local, PLASIC_NUGP, PLASIC_NHOR, PLASIC_DEEP_NLEV);
    free(local);
    return status == PLASIC_RESTART_OK ? PLASIC_RUNTIME_OK
                                       : PLASIC_RUNTIME_IO_FAILED;
}

int physics_write_restart(void)
{
    int status;

#define WRITE_PHYSICS_INTEGER(name, value)                                     \
    do                                                                          \
    {                                                                           \
        status = restart_write_root_integer((name), (value));                   \
        if (status != PLASIC_RESTART_OK)                                        \
        {                                                                       \
            return PLASIC_RUNTIME_IO_FAILED;                                    \
        }                                                                       \
    } while (0)
#define WRITE_PHYSICS_ARRAY(name, values, columns)                              \
    do                                                                          \
    {                                                                           \
        status = restart_write_distributed_grid(                               \
            (name), (values), PLASIC_NUGP, PLASIC_NHOR, (columns));             \
        if (status != PLASIC_RESTART_OK)                                        \
        {                                                                       \
            return PLASIC_RUNTIME_IO_FAILED;                                    \
        }                                                                       \
    } while (0)

    /* LANDMOD */
    WRITE_PHYSICS_INTEGER("nlsoil", PLASIC_LAND_SOIL_LEVELS);
    WRITE_PHYSICS_ARRAY("dtsl", land.surface_temperature, 1);
    WRITE_PHYSICS_ARRAY("dtsm", land.previous_surface_temperature, 1);
    WRITE_PHYSICS_ARRAY("dqs", land.surface_humidity, 1);
    WRITE_PHYSICS_ARRAY("dwmax", plasic_dwmax, 1);
    WRITE_PHYSICS_ARRAY("dtcl", land.temperature_climatology,
                        PLASIC_LAND_CLIMATOLOGY_MONTHS);
    WRITE_PHYSICS_ARRAY("dwcl", land.wetness_climatology,
                        PLASIC_LAND_CLIMATOLOGY_MONTHS);
    WRITE_PHYSICS_ARRAY("dsnowt", land.snow_temperature, 1);
    WRITE_PHYSICS_ARRAY("dsnowz", land.snow_depth, 1);
    WRITE_PHYSICS_ARRAY("dsoilt", land.soil_temperature,
                        PLASIC_LAND_SOIL_LEVELS);
    WRITE_PHYSICS_ARRAY("dglac", plasic_dglac, 1);
    WRITE_PHYSICS_ARRAY("dz0clim", land.roughness_climatology, 1);
    WRITE_PHYSICS_ARRAY("dalbcl", land.albedo_climatology,
                        PLASIC_LAND_CLIMATOLOGY_MONTHS);

    /* SEAMOD */
    WRITE_PHYSICS_INTEGER("naccua", sea.accumulation_count);
    WRITE_PHYSICS_ARRAY("dts", sea.surface_temperature, 1);
    WRITE_PHYSICS_ARRAY("cheata", sea.heat_accumulator, 1);
    WRITE_PHYSICS_ARRAY("cpmea", sea.freshwater_accumulator, 1);
    WRITE_PHYSICS_ARRAY("cprsa", sea.snowfall_accumulator, 1);
    WRITE_PHYSICS_ARRAY("ctauxa", sea.stress_x_accumulator, 1);
    WRITE_PHYSICS_ARRAY("ctauya", sea.stress_y_accumulator, 1);
    WRITE_PHYSICS_ARRAY("cust3a", sea.friction_accumulator, 1);
    WRITE_PHYSICS_ARRAY("cshfla", sea.sensible_heat_accumulator, 1);
    WRITE_PHYSICS_ARRAY("cshdta", sea.sensible_derivative_accumulator, 1);
    WRITE_PHYSICS_ARRAY("clhfla", sea.latent_heat_accumulator, 1);
    WRITE_PHYSICS_ARRAY("clhdta", sea.latent_derivative_accumulator, 1);
    WRITE_PHYSICS_ARRAY("cswfla", sea.shortwave_accumulator, 1);
    WRITE_PHYSICS_ARRAY("clwfla", sea.longwave_accumulator, 1);

    /* ICEMOD */
    WRITE_PHYSICS_INTEGER("naccuice", ice.output_accumulation_count);
    WRITE_PHYSICS_INTEGER("naccuo", ice.ocean_accumulation_count);
    WRITE_PHYSICS_INTEGER("nicec2d",
                          ice.diagnosed_climatological_thickness);
    WRITE_PHYSICS_ARRAY("xls", ice.land_mask, 1);
    WRITE_PHYSICS_ARRAY("xts", ice.surface_temperature, 1);
    WRITE_PHYSICS_ARRAY("xicec", ice.ice_cover, 1);
    WRITE_PHYSICS_ARRAY("xiced", ice.ice_thickness, 1);
    WRITE_PHYSICS_ARRAY("xsnow", ice.snow, 1);
    WRITE_PHYSICS_ARRAY("xclicec", ice.climatological_ice_cover,
                        PLASIC_CLIMATOLOGY_MONTHS);
    WRITE_PHYSICS_ARRAY("xcliced", ice.climatological_ice_thickness,
                        PLASIC_CLIMATOLOGY_MONTHS);
    WRITE_PHYSICS_ARRAY("xclsst", ice.climatological_sst,
                        PLASIC_CLIMATOLOGY_MONTHS);
    WRITE_PHYSICS_ARRAY("cheat", ice.ocean_heat_accumulator, 1);
    WRITE_PHYSICS_ARRAY("cpme", ice.ocean_freshwater_accumulator, 1);
    WRITE_PHYSICS_ARRAY("ctaux", ice.ocean_stress_x_accumulator, 1);
    WRITE_PHYSICS_ARRAY("ctauy", ice.ocean_stress_y_accumulator, 1);
    WRITE_PHYSICS_ARRAY("cust3", ice.ocean_friction_accumulator, 1);
    WRITE_PHYSICS_ARRAY("csnow", ice.ocean_snow_accumulator, 1);
    WRITE_PHYSICS_ARRAY("xflxicea", ice.flux_correction_accumulator, 1);
    WRITE_PHYSICS_ARRAY("xheata", ice.atmospheric_heat_accumulator, 1);
    WRITE_PHYSICS_ARRAY("xcfluxr", ice.maximum_thickness_flux, 1);
    WRITE_PHYSICS_ARRAY("xofluxa",
                        ice.ocean_to_ice_flux_accumulator, 1);
    WRITE_PHYSICS_ARRAY("xqmelta", ice.melt_energy_accumulator, 1);
    WRITE_PHYSICS_ARRAY("xcfluxa", ice.ice_to_ocean_flux_accumulator, 1);
    WRITE_PHYSICS_ARRAY("xcfluxra",
                        ice.maximum_thickness_flux_accumulator, 1);
    WRITE_PHYSICS_ARRAY("xcfluxna",
                        ice.excess_melt_flux_accumulator, 1);
    WRITE_PHYSICS_ARRAY("xsmelta", ice.snow_melt_accumulator, 1);
    WRITE_PHYSICS_ARRAY("ximelta", ice.ice_energy_flux_accumulator, 1);
    WRITE_PHYSICS_ARRAY("xtsfluxa", ice.surface_storage_accumulator, 1);
    WRITE_PHYSICS_ARRAY("xfluxca",
                        ice.diagnosed_conductive_accumulator, 1);
    WRITE_PHYSICS_ARRAY("xscflxa", ice.snow_to_ice_flux_accumulator, 1);
    WRITE_PHYSICS_ARRAY("xcpmea",
                        ice.freshwater_diagnostic_accumulator, 1);
    WRITE_PHYSICS_ARRAY("xstoia", ice.snow_to_ice_rate_accumulator, 1);
    WRITE_PHYSICS_ARRAY("xicecc", ice.compactness, 1);

    /* OCEANMOD */
    WRITE_PHYSICS_INTEGER("nlev_oce", 1);
    WRITE_PHYSICS_INTEGER("naccuoce", ocean.accumulation_count);
    WRITE_PHYSICS_ARRAY("yls", ocean.land_mask, 1);
    WRITE_PHYSICS_ARRAY("ysst", ocean.temperature, 1);
    WRITE_PHYSICS_ARRAY("yiflux", ocean.ice_flux, 1);
    WRITE_PHYSICS_ARRAY("yclsst", ocean.climatological_sst,
                        PLASIC_CLIMATOLOGY_MONTHS);
    WRITE_PHYSICS_ARRAY("yheata", ocean.ice_to_ocean_flux_accumulator, 1);
    WRITE_PHYSICS_ARRAY("yifluxa", ocean.ice_flux_accumulator, 1);
    WRITE_PHYSICS_ARRAY("yfldoa",
                        ocean.deep_ocean_flux_accumulator, 1);
    WRITE_PHYSICS_ARRAY("yfldo", ocean.deep_ocean_flux, 1);

    /* Three-dimensional deep ocean. */
    if (deep_ocean.enabled != 0 && deep_ocean.initialized != 0)
    {
        WRITE_PHYSICS_INTEGER("deep_nlev", PLASIC_DEEP_NLEV);
        if (write_deep_restart_field("deept", deep_ocean.temperature) !=
                PLASIC_RUNTIME_OK ||
            write_deep_restart_field("deeps", deep_ocean.salinity) !=
                PLASIC_RUNTIME_OK ||
            write_deep_restart_field("deepu", deep_ocean.u_velocity) !=
                PLASIC_RUNTIME_OK ||
            write_deep_restart_field("deepv", deep_ocean.v_velocity) !=
                PLASIC_RUNTIME_OK ||
            write_deep_restart_field("deepw", deep_ocean.w_velocity) !=
                PLASIC_RUNTIME_OK)
        {
            return PLASIC_RUNTIME_IO_FAILED;
        }
    }

#undef WRITE_PHYSICS_ARRAY
#undef WRITE_PHYSICS_INTEGER
    return PLASIC_RUNTIME_OK;
}

int physics_frame_field(
    int32_t variable_id, const float **values, int32_t *layers)
{
    static float deep_temperature_local[
        PLASIC_NHOR * PLASIC_DEEP_NLEV];
    static float deep_salinity_local[
        PLASIC_NHOR * PLASIC_DEEP_NLEV];

    if (values == NULL || layers == NULL)
    {
        return PLASIC_RUNTIME_BAD_ARGUMENT;
    }
    *layers = 1;
    switch (variable_id)
    {
    case PLASIC_STREAM_SEA_SURFACE_TEMPERATURE:
        *values = ocean.temperature;
        return PLASIC_RUNTIME_OK;
    case PLASIC_STREAM_DEEP_OCEAN_TEMPERATURE:
        /*
         * Report the deep field as unavailable when the three-dimensional
         * ocean is disabled or not initialized yet, so callers can skip it.
         */
        if (deep_ocean.enabled == 0 || deep_ocean.initialized == 0)
        {
            *values = NULL;
            *layers = 0;
            return PLASIC_RUNTIME_BAD_ARGUMENT;
        }
        physics_ocean3d_pack_local(
            deep_ocean.temperature, deep_temperature_local,
            PLASIC_DEEP_NLEV);
        *values = deep_temperature_local;
        *layers = PLASIC_DEEP_NLEV;
        return PLASIC_RUNTIME_OK;
    case PLASIC_STREAM_DEEP_OCEAN_SALINITY:
        if (deep_ocean.enabled == 0 || deep_ocean.initialized == 0)
        {
            *values = NULL;
            *layers = 0;
            return PLASIC_RUNTIME_BAD_ARGUMENT;
        }
        physics_ocean3d_pack_local(
            deep_ocean.salinity, deep_salinity_local,
            PLASIC_DEEP_NLEV);
        *values = deep_salinity_local;
        *layers = PLASIC_DEEP_NLEV;
        return PLASIC_RUNTIME_OK;
    default:
        *values = NULL;
        *layers = 0;
        return PLASIC_RUNTIME_BAD_ARGUMENT;
    }
}
