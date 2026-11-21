/* Monthly ozone on the native [lev][lat][lon] grid. Only rank 0 uses
 * NetCDF; one validated month is shared when the model date changes.
 * With no file, no buffer is allocated and radiation keeps its own profile. */
#include "ozone_forcing.h"
#include "runtime_internal.h"
#include "plasic_config.h"

#include <netcdf.h>

static struct
{
    char path[PLASIC_CONFIG_PATH_BYTES];
    int ncid, o3;
    int64_t first_month, active_month;
    size_t records;
    double fill, missing;
    double *raw;     /* root only, one full horizontal level */
    float *plane;    /* broadcast one level at a time */
} forcing = {.ncid = -1, .active_month = -1};

static int file_error(const char *detail)
{
    runtime_set_error("ozone forcing '%s': %s", forcing.path, detail);
    return PLASIC_RUNTIME_IO_FAILED;
}

/* Shape is checked before NetCDF copies dimension IDs or coordinate values
 * into fixed-size buffers. Never interpret packed coordinates or ozone as
 * physical values: netCDF-C does not apply CF scale/offset automatically. */
static int variable(const char *name, int rank, const int *dims, int *id,
                    int integer)
{
    int ndims, actual[4];
    nc_type type;
    if (nc_inq_varid(forcing.ncid, name, id) != NC_NOERR ||
        nc_inq_varndims(forcing.ncid, *id, &ndims) != NC_NOERR ||
        ndims != rank ||
        nc_inq_vardimid(forcing.ncid, *id, actual) != NC_NOERR ||
        memcmp(actual, dims, (size_t)rank * sizeof(*dims)) != 0 ||
        nc_inq_vartype(forcing.ncid, *id, &type) != NC_NOERR ||
        (integer ? type != NC_INT : (type != NC_FLOAT && type != NC_DOUBLE)) ||
        nc_inq_att(forcing.ncid, *id, "scale_factor", NULL, NULL) != NC_ENOTATT ||
        nc_inq_att(forcing.ncid, *id, "add_offset", NULL, NULL) != NC_ENOTATT)
    {
        runtime_set_error("ozone forcing: invalid shape, type or packing for '%s'", name);
        return PLASIC_RUNTIME_IO_FAILED;
    }
    return PLASIC_RUNTIME_OK;
}

static int has_units(int id, const char *units)
{
    char value[32];
    nc_type type;
    size_t len;
    return nc_inq_att(forcing.ncid, id, "units", &type, &len) == NC_NOERR &&
           type == NC_CHAR && len == strlen(units) && len <= sizeof(value) &&
           nc_get_att_text(forcing.ncid, id, "units", value) == NC_NOERR &&
           memcmp(value, units, len) == 0;
}

static int read_marker(int id, const char *name, double *value)
{
    size_t count;
    int status = nc_inq_attlen(forcing.ncid, id, name, &count);
    *value = NAN;
    if (status == NC_ENOTATT)
        return PLASIC_RUNTIME_OK;
    if (status != NC_NOERR || count != 1 ||
        nc_get_att_double(forcing.ncid, id, name, value) != NC_NOERR)
        return file_error("expected a numeric scalar _FillValue or missing_value");
    return PLASIC_RUNTIME_OK;
}

static int check_header(void)
{
    const char *const names[] = {"time", "lev", "lat", "lon"};
    const char *const units[] = {"", "1", "degrees_north", "degrees_east"};
    const size_t sizes[] = {0, PLASIC_NLEV, PLASIC_NLAT, PLASIC_NLON};
    int dims[4], date;
    for (int axis = 0; axis < 4; ++axis)
    {
        size_t length;
        int id;
        if (nc_inq_dimid(forcing.ncid, names[axis], &dims[axis]) != NC_NOERR ||
            nc_inq_dimlen(forcing.ncid, dims[axis], &length) != NC_NOERR ||
            length == 0 || (axis != 0 && length != sizes[axis]))
            return file_error("expected time and model-sized lev, lat, lon dimensions");
        if (axis == 0)
        {
            forcing.records = length;
            continue;
        }
        if (variable(names[axis], 1, &dims[axis], &id, 0) != PLASIC_RUNTIME_OK)
            return PLASIC_RUNTIME_IO_FAILED;
        if (!has_units(id, units[axis]))
            return file_error("unexpected lev, lat or lon coordinate units");
        double fill, missing;
        if (read_marker(id, "_FillValue", &fill) != PLASIC_RUNTIME_OK ||
            read_marker(id, "missing_value", &missing) != PLASIC_RUNTIME_OK)
            return PLASIC_RUNTIME_IO_FAILED;
        for (size_t i = 0; i < length; ++i)
        {
            double value;
            const double expected = axis == 1 ? (double)plasic_sigma[i] :
                axis == 2 ? asin(plasic_sid[i]) * (180.0 / 3.14159265358979323846) :
                            (double)i * 360.0 / PLASIC_NLON;
            const double tolerance = axis == 1 ? 1e-6 : 1e-4;
            if (nc_get_var1_double(forcing.ncid, id, &i, &value) != NC_NOERR ||
                !isfinite(value) || value == fill || value == missing ||
                fabs(value - expected) > tolerance)
            {
                runtime_set_error("ozone forcing: '%s' coordinate at index %zu "
                                  "does not match the model grid", names[axis], i);
                return PLASIC_RUNTIME_IO_FAILED;
            }
        }
    }

    if (variable("yyyymm", 1, dims, &date, 1) != PLASIC_RUNTIME_OK ||
        variable("o3", 4, dims, &forcing.o3, 0) != PLASIC_RUNTIME_OK)
        return PLASIC_RUNTIME_IO_FAILED;
    if (!has_units(forcing.o3, "kg kg-1"))
        return file_error("o3 must be mass fraction with units 'kg kg-1'");
    if (read_marker(forcing.o3, "_FillValue", &forcing.fill) != PLASIC_RUNTIME_OK ||
        read_marker(forcing.o3, "missing_value", &forcing.missing) != PLASIC_RUNTIME_OK)
        return PLASIC_RUNTIME_IO_FAILED;

    double date_fill, date_missing;
    if (read_marker(date, "_FillValue", &date_fill) != PLASIC_RUNTIME_OK ||
        read_marker(date, "missing_value", &date_missing) != PLASIC_RUNTIME_OK)
        return PLASIC_RUNTIME_IO_FAILED;
    for (size_t i = 0; i < forcing.records; ++i)
    {
        int code;
        if (nc_get_var1_int(forcing.ncid, date, &i, &code) != NC_NOERR ||
            code == date_fill || code == date_missing ||
            code < 1 || code % 100 < 1 || code % 100 > 12)
            return file_error("invalid yyyymm month code");
        int64_t month = (int64_t)(code / 100) * 12 + code % 100 - 1;
        if (i == 0)
            forcing.first_month = month;
        if (month < forcing.first_month || (uint64_t)(month - forcing.first_month) != i)
            return file_error("yyyymm must list consecutive months in ascending order");
    }
    return PLASIC_RUNTIME_OK;
}

int runtime_ozone_forcing_open(const char *path)
{
    if (path == NULL)
        return PLASIC_RUNTIME_OK;
    if (path[0] == '\0' || strlen(path) >= sizeof(forcing.path))
    {
        runtime_set_error("expected a nonempty ozone forcing file path");
        return PLASIC_RUNTIME_BAD_ARGUMENT;
    }
    strcpy(forcing.path, path);
    int status = PLASIC_RUNTIME_OK;
    if (mp_is_root())
    {
        int ncid;
        int result = nc_open(path, NC_NOWRITE, &ncid);
        if (result != NC_NOERR)
            status = file_error(nc_strerror(result));
        else
        {
            forcing.ncid = ncid;
            status = check_header();
        }
    }
    status = runtime_broadcast_status(status);
    if (status != PLASIC_RUNTIME_OK)
    {
        if (!mp_is_root())
            (void)file_error("see rank 0 for the ozone input error");
        return status;
    }
    if (mp_is_root())
        forcing.raw = malloc((size_t)PLASIC_NUGP * sizeof(*forcing.raw));
    forcing.plane = malloc((size_t)PLASIC_NUGP * sizeof(*forcing.plane));
    float failed = forcing.plane == NULL || (mp_is_root() && forcing.raw == NULL);
    if (mp_allreduce_real(&failed, 1) != 0 || failed != 0.0f)
        return file_error("could not allocate the ozone level buffer");
    return PLASIC_RUNTIME_OK;
}

static int read_level(int year, int month, int64_t index, size_t k)
{
    /* Decode in double so an out-of-range or missing input is checked before
     * narrowing. Default NetCDF floating fill values are outside [0,1]. */
    const size_t count[] = {1, 1, PLASIC_NLAT, PLASIC_NLON};
    if (index < 0 || (uint64_t)index >= forcing.records)
    {
        runtime_set_error("ozone forcing does not cover model date %d-%02d", year, month);
        return PLASIC_RUNTIME_BAD_ARGUMENT;
    }
    const size_t start[] = {(size_t)index, k, 0, 0};
    int status = nc_get_vara_double(forcing.ncid, forcing.o3, start, count, forcing.raw);
    if (status != NC_NOERR)
        return file_error(nc_strerror(status));
    for (size_t i = 0; i < PLASIC_NUGP; ++i)
    {
        double q = forcing.raw[i];
        if (!isfinite(q) || q < 0.0 || q > 1.0 ||
            q == forcing.fill || q == forcing.missing)
        {
            runtime_set_error("ozone forcing: invalid o3 at %d-%02d "
                              "(lev=%zu, lat=%zu, lon=%zu; zero-based)",
                              year, month, k, i / PLASIC_NLON, i % PLASIC_NLON);
            return PLASIC_RUNTIME_IO_FAILED;
        }
        forcing.plane[i] = (float)q;
    }
    return PLASIC_RUNTIME_OK;
}

int runtime_ozone_forcing_apply(int year, int month)
{
    const int64_t key = (int64_t)year * 12 + month - 1;
    if (forcing.plane == NULL || forcing.active_month == key)
        return PLASIC_RUNTIME_OK;
    for (size_t k = 0; k < PLASIC_NLEV; ++k)
    {
        int status = PLASIC_RUNTIME_OK;
        if (mp_is_root())
            status = read_level(year, month, key - forcing.first_month, k);
        status = runtime_broadcast_status(status);
        if (status != PLASIC_RUNTIME_OK)
        {
            if (!mp_is_root())
                (void)file_error("see rank 0 for the ozone input error");
            return status;
        }
        if (mp_broadcast_real(forcing.plane, PLASIC_NUGP) != 0)
            return file_error("could not broadcast the ozone level");
        plasic_physics_set_ozone_level(k, forcing.plane);
    }
    forcing.active_month = key;
    return PLASIC_RUNTIME_OK;
}

const char *runtime_ozone_forcing_source(void)
{
    return forcing.plane != NULL ? forcing.path : NULL;
}

void runtime_ozone_forcing_close(void)
{
    if (forcing.ncid >= 0)
        (void)nc_close(forcing.ncid);
    free(forcing.raw);
    free(forcing.plane);
    forcing.raw = NULL;
    forcing.plane = NULL;
    forcing.ncid = -1;
    forcing.active_month = -1;
}
