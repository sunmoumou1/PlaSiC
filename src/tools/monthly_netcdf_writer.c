#include "monthly_netcdf_writer.h"

#include <netcdf.h>
#include <stdio.h>
#include <string.h>

enum
{
    MONTHLY_NETCDF_PATH_BYTES = 4096,
    MONTHLY_NETCDF_BOUNDS = 2
};

static void record_error(
    monthly_netcdf_error *error, const char *operation,
    const char *variable, int nc_status)
{
    if (error != NULL)
    {
        error->operation = operation;
        error->variable = variable;
        error->nc_status = nc_status;
    }
}

static int put_text_attribute(
    int ncid, int variable_id, const char *name, const char *value,
    monthly_netcdf_error *error)
{
    const int status = nc_put_att_text(
        ncid, variable_id, name, strlen(value), value);

    if (status != NC_NOERR)
    {
        record_error(error, "attribute write", name, status);
        return MONTHLY_NETCDF_DEFINE_FAILED;
    }
    return MONTHLY_NETCDF_OK;
}

int monthly_netcdf_create(
    monthly_netcdf_file *file, const monthly_netcdf_spec *spec,
    const char *directory, const char *co2_forcing_source,
    const char *ozone_forcing_source, const char *calendar,
    const float *latitudes, int32_t latitude_count,
    const float *longitudes, int32_t longitude_count,
    const float *levels, int32_t level_count,
    monthly_netcdf_error *error)
{
    char path[MONTHLY_NETCDF_PATH_BYTES];
    int has_levels;
    const char *vertical_name;
    int dimensions[4];
    int time_dimension;
    int level_dimension = -1;
    int latitude_dimension;
    int longitude_dimension;
    int bounds_dimension;
    int latitude_id;
    int longitude_id;
    int level_id = -1;
    int status;
    float fill_value = 1.0e20f;

    if (file == NULL || spec == NULL || directory == NULL ||
        co2_forcing_source == NULL || ozone_forcing_source == NULL ||
        calendar == NULL || calendar[0] == '\0' ||
        (latitude_count > 0 && latitudes == NULL) ||
        (longitude_count > 0 && longitudes == NULL) ||
        (spec->layers > 1 && (levels == NULL || level_count != spec->layers)))
    {
        record_error(
            error, "argument validation",
            spec != NULL ? spec->name : "monthly", NC_NOERR);
        return MONTHLY_NETCDF_BAD_ARGUMENT;
    }
    has_levels = spec->layers > 1;
    vertical_name = spec->depth_axis ? "depth" : "lev";

    file->ncid = -1;
    file->time_id = -1;
    file->bounds_id = -1;
    file->data_id = -1;
    file->records = 0;
    file->spec = *spec;
    file->latitude_count = latitude_count;
    file->longitude_count = longitude_count;

    if (snprintf(path, sizeof(path), "%s/%s.nc", directory, spec->name) >=
        (int)sizeof(path))
    {
        return MONTHLY_NETCDF_PATH_TOO_LONG;
    }
    status = nc_create(path, NC_NETCDF4 | NC_CLOBBER, &file->ncid);
    if (status != NC_NOERR)
    {
        record_error(error, "create", spec->name, status);
        return MONTHLY_NETCDF_CREATE_FAILED;
    }
    if ((status = nc_def_dim(file->ncid, "time", NC_UNLIMITED,
                             &time_dimension)) != NC_NOERR ||
        (status = nc_def_dim(file->ncid, "lat", (size_t)latitude_count,
                             &latitude_dimension)) != NC_NOERR ||
        (status = nc_def_dim(file->ncid, "lon", (size_t)longitude_count,
                             &longitude_dimension)) != NC_NOERR ||
        (status = nc_def_dim(file->ncid, "bnds", MONTHLY_NETCDF_BOUNDS,
                             &bounds_dimension)) != NC_NOERR)
    {
        record_error(error, "dimension definition", spec->name, status);
        return MONTHLY_NETCDF_DEFINE_FAILED;
    }
    if (has_levels &&
        (status = nc_def_dim(file->ncid, vertical_name,
                             (size_t)spec->layers,
                             &level_dimension)) != NC_NOERR)
    {
        record_error(error, "level definition", spec->name, status);
        return MONTHLY_NETCDF_DEFINE_FAILED;
    }
    if ((status = nc_def_var(file->ncid, "time", NC_DOUBLE, 1,
                             &time_dimension, &file->time_id)) != NC_NOERR ||
        (status = nc_def_var(file->ncid, "time_bnds", NC_DOUBLE, 2,
                             (int[]){time_dimension, bounds_dimension},
                             &file->bounds_id)) != NC_NOERR ||
        (status = nc_def_var(file->ncid, "lat", NC_FLOAT, 1,
                             &latitude_dimension, &latitude_id)) != NC_NOERR ||
        (status = nc_def_var(file->ncid, "lon", NC_FLOAT, 1,
                             &longitude_dimension, &longitude_id)) != NC_NOERR)
    {
        record_error(error, "coordinate definition", spec->name, status);
        return MONTHLY_NETCDF_DEFINE_FAILED;
    }
    if (has_levels &&
        (status = nc_def_var(file->ncid, vertical_name, NC_FLOAT, 1,
                             &level_dimension, &level_id)) != NC_NOERR)
    {
        record_error(error, "vertical coordinate definition", spec->name, status);
        return MONTHLY_NETCDF_DEFINE_FAILED;
    }
    if (has_levels)
    {
        dimensions[0] = time_dimension;
        dimensions[1] = level_dimension;
        dimensions[2] = latitude_dimension;
        dimensions[3] = longitude_dimension;
    }
    else
    {
        dimensions[0] = time_dimension;
        dimensions[1] = latitude_dimension;
        dimensions[2] = longitude_dimension;
    }
    status = nc_def_var(
        file->ncid, spec->name, NC_FLOAT, has_levels ? 4 : 3,
        dimensions, &file->data_id);
    if (status != NC_NOERR ||
        (status = nc_put_att_float(file->ncid, file->data_id,
                                   "_FillValue", NC_FLOAT, 1,
                                   &fill_value)) != NC_NOERR ||
        (status = nc_def_var_deflate(file->ncid, file->data_id,
                                     1, 1, 2)) != NC_NOERR)
    {
        record_error(error, "data variable definition", spec->name, status);
        return MONTHLY_NETCDF_DEFINE_FAILED;
    }
    if (put_text_attribute(
            file->ncid, file->time_id, "standard_name", "time", error) !=
            MONTHLY_NETCDF_OK ||
        put_text_attribute(
            file->ncid, file->time_id, "units",
            "days since 1850-01-01 00:00:00", error) != MONTHLY_NETCDF_OK ||
        put_text_attribute(
            file->ncid, file->time_id, "calendar",
            calendar, error) != MONTHLY_NETCDF_OK ||
        put_text_attribute(
            file->ncid, file->time_id, "bounds", "time_bnds", error) !=
            MONTHLY_NETCDF_OK ||
        put_text_attribute(
            file->ncid, latitude_id, "standard_name", "latitude", error) !=
            MONTHLY_NETCDF_OK ||
        put_text_attribute(
            file->ncid, latitude_id, "units", "degrees_north", error) !=
            MONTHLY_NETCDF_OK ||
        put_text_attribute(
            file->ncid, longitude_id, "standard_name", "longitude", error) !=
            MONTHLY_NETCDF_OK ||
        put_text_attribute(
            file->ncid, longitude_id, "units", "degrees_east", error) !=
            MONTHLY_NETCDF_OK ||
        put_text_attribute(
            file->ncid, file->data_id, "standard_name",
            spec->standard_name, error) != MONTHLY_NETCDF_OK ||
        put_text_attribute(
            file->ncid, file->data_id, "long_name",
            spec->long_name, error) != MONTHLY_NETCDF_OK ||
        put_text_attribute(
            file->ncid, file->data_id, "units", spec->units, error) !=
            MONTHLY_NETCDF_OK ||
        put_text_attribute(
            file->ncid, file->data_id, "cell_methods", "time: mean", error) !=
            MONTHLY_NETCDF_OK ||
        put_text_attribute(
            file->ncid, NC_GLOBAL, "Conventions", "CF-1.10", error) !=
            MONTHLY_NETCDF_OK ||
        put_text_attribute(
            file->ncid, NC_GLOBAL, "title",
            "PlaSiC monthly simulation output", error) !=
            MONTHLY_NETCDF_OK ||
        put_text_attribute(
            file->ncid, NC_GLOBAL, "frequency", "mon", error) !=
            MONTHLY_NETCDF_OK ||
        put_text_attribute(
            file->ncid, NC_GLOBAL, "co2_forcing_source",
            co2_forcing_source, error) != MONTHLY_NETCDF_OK ||
        put_text_attribute(
            file->ncid, NC_GLOBAL, "ozone_forcing_source",
            ozone_forcing_source, error) != MONTHLY_NETCDF_OK)
    {
        return MONTHLY_NETCDF_DEFINE_FAILED;
    }
    if (has_levels &&
        (put_text_attribute(
             file->ncid, level_id, "long_name",
             spec->depth_axis ? "ocean depth" : "sigma level", error) !=
             MONTHLY_NETCDF_OK ||
         put_text_attribute(
             file->ncid, level_id, "standard_name",
             spec->depth_axis ? "depth" : "atmosphere_sigma_coordinate",
             error) != MONTHLY_NETCDF_OK ||
         put_text_attribute(
             file->ncid, level_id, "units",
             spec->depth_axis ? "m" : "1", error) != MONTHLY_NETCDF_OK ||
         put_text_attribute(
             file->ncid, level_id, "axis", "Z", error) !=
             MONTHLY_NETCDF_OK ||
         put_text_attribute(
             file->ncid, level_id, "positive", "down", error) !=
             MONTHLY_NETCDF_OK))
    {
        return MONTHLY_NETCDF_DEFINE_FAILED;
    }
    status = nc_enddef(file->ncid);
    if (status != NC_NOERR)
    {
        record_error(error, "end definition", spec->name, status);
        return MONTHLY_NETCDF_DEFINE_FAILED;
    }
    status = nc_put_var_float(file->ncid, latitude_id, latitudes);
    if (status == NC_NOERR)
    {
        status = nc_put_var_float(file->ncid, longitude_id, longitudes);
    }
    if (status == NC_NOERR && has_levels)
    {
        status = nc_put_var_float(file->ncid, level_id, levels);
    }
    if (status != NC_NOERR)
    {
        record_error(error, "coordinate write", spec->name, status);
        return MONTHLY_NETCDF_WRITE_FAILED;
    }
    return MONTHLY_NETCDF_OK;
}

int monthly_netcdf_write_record(
    monthly_netcdf_file *file, double time_value, const double bounds[2],
    const float *data, monthly_netcdf_error *error)
{
    size_t time_start[1];
    size_t time_count[1];
    size_t bounds_start[2];
    size_t bounds_count[2];
    size_t data_start[4];
    size_t data_count[4];
    int status;

    if (file == NULL || data == NULL)
    {
        record_error(error, "argument validation", "monthly", NC_NOERR);
        return MONTHLY_NETCDF_BAD_ARGUMENT;
    }
    time_start[0] = file->records;
    time_count[0] = 1;
    bounds_start[0] = file->records;
    bounds_start[1] = 0;
    bounds_count[0] = 1;
    bounds_count[1] = 2;
    data_start[0] = file->records;
    data_start[1] = 0;
    data_start[2] = 0;
    data_start[3] = 0;
    data_count[0] = 1;
    data_count[1] = (size_t)file->spec.layers;
    data_count[2] = (size_t)file->latitude_count;
    data_count[3] = (size_t)file->longitude_count;
    if (file->spec.layers == 1)
    {
        data_count[1] = (size_t)file->latitude_count;
        data_count[2] = (size_t)file->longitude_count;
    }
    status = nc_put_vara_double(
        file->ncid, file->time_id, time_start, time_count, &time_value);
    if (status == NC_NOERR)
    {
        status = nc_put_vara_double(
            file->ncid, file->bounds_id, bounds_start, bounds_count, bounds);
    }
    if (status == NC_NOERR)
    {
        status = nc_put_vara_float(
            file->ncid, file->data_id, data_start, data_count, data);
    }
    if (status == NC_NOERR)
    {
        status = nc_sync(file->ncid);
    }
    if (status != NC_NOERR)
    {
        record_error(error, "monthly write", file->spec.name, status);
        return MONTHLY_NETCDF_WRITE_FAILED;
    }
    ++file->records;
    return MONTHLY_NETCDF_OK;
}

int monthly_netcdf_close(monthly_netcdf_file *file, monthly_netcdf_error *error)
{
    if (file == NULL)
    {
        record_error(error, "argument validation", "monthly", NC_NOERR);
        return MONTHLY_NETCDF_BAD_ARGUMENT;
    }
    if (file->ncid >= 0)
    {
        const int status = nc_close(file->ncid);

        file->ncid = -1;
        if (status != NC_NOERR)
        {
            record_error(error, "close", file->spec.name, status);
            return MONTHLY_NETCDF_CLOSE_FAILED;
        }
    }
    return MONTHLY_NETCDF_OK;
}

void monthly_netcdf_abort(monthly_netcdf_file *file)
{
    if (file != NULL)
    {
        if (file->ncid >= 0)
        {
            (void)nc_close(file->ncid);
        }
        file->ncid = -1;
        file->records = 0;
    }
}

const char *monthly_netcdf_nc_error(int nc_status)
{
    return nc_strerror(nc_status);
}
