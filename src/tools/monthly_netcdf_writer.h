#ifndef PLASIC_MONTHLY_NETCDF_WRITER_H
#define PLASIC_MONTHLY_NETCDF_WRITER_H

#include <stddef.h>
#include <stdint.h>

enum monthly_netcdf_status
{
    MONTHLY_NETCDF_OK = 0,
    MONTHLY_NETCDF_BAD_ARGUMENT = 1,
    MONTHLY_NETCDF_PATH_TOO_LONG = 2,
    MONTHLY_NETCDF_CREATE_FAILED = 3,
    MONTHLY_NETCDF_DEFINE_FAILED = 4,
    MONTHLY_NETCDF_WRITE_FAILED = 5,
    MONTHLY_NETCDF_CLOSE_FAILED = 6
};

/* NetCDF shape and CF metadata of one monthly output variable. */
typedef struct monthly_netcdf_spec
{
    const char *name;
    const char *standard_name;
    const char *long_name;
    const char *units;
    int32_t layers;
    int32_t depth_axis;
} monthly_netcdf_spec;

/* Diagnostic detail for a failed writer call. */
typedef struct monthly_netcdf_error
{
    const char *operation;
    const char *variable;
    int nc_status;
} monthly_netcdf_error;

/* Open NetCDF handle owned by one writer call sequence. */
typedef struct monthly_netcdf_file
{
    monthly_netcdf_spec spec;
    int ncid;
    int time_id;
    int bounds_id;
    int data_id;
    int32_t latitude_count;
    int32_t longitude_count;
    size_t records;
} monthly_netcdf_file;

/*
 * Create the CF-NetCDF file, dimensions, coordinates and the data variable for
 * one monthly output field. Coordinates are passed explicitly so the writer
 * holds no model or runtime state; a non-NULL levels array is required when
 * the spec has more than one layer. `calendar` is the CF calendar name of the
 * time axis ("proleptic_gregorian" or "360_day"). On any failure after
 * nc_create() the file stays open and must be released with
 * monthly_netcdf_abort() or monthly_netcdf_close().
 */
int monthly_netcdf_create(
    monthly_netcdf_file *file, const monthly_netcdf_spec *spec,
    const char *directory, const char *co2_forcing_source,
    const char *ozone_forcing_source, const char *calendar,
    const float *latitudes, int32_t latitude_count,
    const float *longitudes, int32_t longitude_count,
    const float *levels, int32_t level_count,
    monthly_netcdf_error *error);

/* Append one averaged record including its time, bounds and coordinates. */
int monthly_netcdf_write_record(
    monthly_netcdf_file *file, double time_value, const double bounds[2],
    const float *data, monthly_netcdf_error *error);

int monthly_netcdf_close(monthly_netcdf_file *file, monthly_netcdf_error *error);

/* Close silently and reset the handle; used on error paths. */
void monthly_netcdf_abort(monthly_netcdf_file *file);

/* Text of a raw NetCDF status code, for caller-side diagnostics. */
const char *monthly_netcdf_nc_error(int nc_status);

#endif
