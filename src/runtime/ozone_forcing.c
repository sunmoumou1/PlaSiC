#include "ozone_forcing.h"

#include "runtime_internal.h"

#include <netcdf.h>
#include <float.h>

/* Each rank reads its own latitude slab. The file stays open across months. */
static struct
{
    int file;
    int variable;
    size_t months;
    int64_t first_month;
    int64_t cached_record[2];
    double mass_factor;
    double fill_value;
    double missing_value;
    double *slabs;
    char source[4096];
} forcing = {.file = -1, .cached_record = {-1, -1}};

static const size_t local_values = (size_t)PLASIC_NHOR * PLASIC_NLEV;

void runtime_ozone_forcing_close(void)
{
    if (forcing.file >= 0)
        (void)nc_close(forcing.file);
    free(forcing.slabs);
    forcing.slabs = NULL;
    forcing.file = -1;
    forcing.cached_record[0] = forcing.cached_record[1] = -1;
    forcing.source[0] = '\0';
}

const char *runtime_ozone_forcing_source(void)
{
    return forcing.file < 0 ? "synthetic" : forcing.source;
}

/* Check rank before retrieving dimension IDs into the fixed-size buffer. */
static int find_variable(
    const char *name, int rank, const int *dimensions, int *id, nc_type *type)
{
    int actual_rank, actual_dims[4];

    return nc_inq_varid(forcing.file, name, id) == NC_NOERR &&
           nc_inq_varndims(forcing.file, *id, &actual_rank) == NC_NOERR &&
           actual_rank == rank &&
           nc_inq_vardimid(forcing.file, *id, actual_dims) == NC_NOERR &&
           memcmp(actual_dims, dimensions, (size_t)rank * sizeof(int)) == 0 &&
           nc_inq_vartype(forcing.file, *id, type) == NC_NOERR &&
           nc_inq_att(forcing.file, *id, "scale_factor", NULL, NULL) == NC_ENOTATT &&
           nc_inq_att(forcing.file, *id, "add_offset", NULL, NULL) == NC_ENOTATT;
}

static int read_month_index(int time_dimension)
{
    int year_id, month_id;
    nc_type year_type, month_type;
    int years[256], months[256];

    if (!find_variable("year", 1, &time_dimension, &year_id, &year_type) ||
        !find_variable("month", 1, &time_dimension, &month_id, &month_type) ||
        year_type != NC_INT || month_type != NC_INT)
        return 0;

    for (size_t start = 0; start < forcing.months;)
    {
        const size_t remaining = forcing.months - start;
        const size_t count = remaining < 256 ? remaining : 256;
        if (nc_get_vara_int(forcing.file, year_id, &start, &count, years) != NC_NOERR ||
            nc_get_vara_int(forcing.file, month_id, &start, &count, months) != NC_NOERR)
            return 0;
        for (size_t i = 0; i < count; ++i)
        {
            const int64_t date = (int64_t)years[i] * 12 + months[i] - 1;
            if (years[i] < 0 || months[i] < 1 || months[i] > 12)
                return 0;
            if (start + i == 0)
                forcing.first_month = date;
            if (date - forcing.first_month != (int64_t)(start + i))
                return 0;
        }
        start += count;
    }
    return 1;
}

static int check_coordinate(const char *name, int dimension, size_t count)
{
    int id;
    nc_type type;
    const int latitude = strcmp(name, "lat") == 0;
    const int sigma = strcmp(name, "lev") == 0;
    const size_t offset = latitude ? (size_t)mp_local_offset(PLASIC_NLPP) : 0;
    double coordinates[PLASIC_NLON > PLASIC_NLEV ? PLASIC_NLON : PLASIC_NLEV];

    if (!find_variable(name, 1, &dimension, &id, &type) ||
        (type != NC_FLOAT && type != NC_DOUBLE) ||
        nc_get_vara_double(forcing.file, id, &offset, &count, coordinates) != NC_NOERR)
        return 0;
    for (size_t i = 0; i < count; ++i)
    {
        const double expected = sigma ? plasic_sigma[i] :
            (latitude ? asin(plasic_sid[i]) * (180.0 / 3.14159265358979323846)
                      : 360.0 * (double)i / PLASIC_NLON);
        const double tolerance = sigma ? 1e-6 : 1e-4;
        if (!isfinite(coordinates[i]) || fabs(coordinates[i] - expected) > tolerance)
            return 0;
    }
    return 1;
}

static int read_missing_value(const char *name, nc_type stored_type, double *value)
{
    size_t length;
    const int status = nc_inq_attlen(forcing.file, forcing.variable, name, &length);
    *value = NAN;
    if (status == NC_ENOTATT)
        return 1;
    if (status != NC_NOERR || length != 1 ||
        nc_get_att_double(forcing.file, forcing.variable, name, value) != NC_NOERR)
        return 0;
    /* Compare a marker against values as stored, before any unit conversion. */
    if (stored_type == NC_FLOAT && fabs(*value) <= FLT_MAX)
        *value = (double)(float)*value;
    return 1;
}

static int read_units(void)
{
    char text[32];
    char *string = NULL;
    const char *units = NULL;
    nc_type type;
    size_t length;

    if (nc_inq_att(forcing.file, forcing.variable, "units", &type, &length) != NC_NOERR)
        return 0;
    if (type == NC_CHAR && length < sizeof(text) &&
        nc_get_att_text(forcing.file, forcing.variable, "units", text) == NC_NOERR)
    {
        text[length] = '\0';
        units = text;
    }
    else if (type == NC_STRING && length == 1 &&
             nc_get_att_string(forcing.file, forcing.variable, "units", &string) == NC_NOERR)
    {
        units = string;
    }

    forcing.mass_factor = 0.0;
    if (units != NULL && strcmp(units, "kg kg-1") == 0)
        forcing.mass_factor = 1.0;
    if (units != NULL && strcmp(units, "mol mol-1") == 0)
        forcing.mass_factor = 47.9982 / 28.9644; /* O3 / dry-air molar masses */
    if (string != NULL)
        (void)nc_free_string(1, &string);
    return forcing.mass_factor != 0.0;
}

static int read_header(void)
{
    const char *names[] = {"time", "lev", "lat", "lon"};
    const size_t sizes[] = {0, PLASIC_NLEV, PLASIC_NLAT, PLASIC_NLON};
    const size_t local_sizes[] = {0, PLASIC_NLEV, PLASIC_NLPP, PLASIC_NLON};
    int dimensions[4];
    nc_type type;
    size_t length;

    for (int axis = 0; axis < 4; ++axis)
    {
        if (nc_inq_dimid(forcing.file, names[axis], &dimensions[axis]) != NC_NOERR ||
            nc_inq_dimlen(forcing.file, dimensions[axis], &length) != NC_NOERR ||
            length == 0 || (axis != 0 && length != sizes[axis]))
        {
            runtime_set_error("invalid ozone forcing dimension '%s' in '%s'", names[axis], forcing.source);
            return 0;
        }
        if (axis == 0)
            forcing.months = length;
        else if (!check_coordinate(names[axis], dimensions[axis], local_sizes[axis]))
        {
            runtime_set_error("ozone forcing coordinate '%s' in '%s' must match the model grid",
                              names[axis], forcing.source);
            return 0;
        }
    }
    if (forcing.months > INT32_MAX || !read_month_index(dimensions[0]))
    {
        runtime_set_error("ozone forcing '%s' requires consecutive integer year(time), month(time) records",
                          forcing.source);
        return 0;
    }
    if (!find_variable("o3", 4, dimensions, &forcing.variable, &type) ||
        (type != NC_FLOAT && type != NC_DOUBLE))
    {
        runtime_set_error("ozone forcing '%s' requires unpacked floating point o3(time,lev,lat,lon)",
                          forcing.source);
        return 0;
    }
    if (!read_units())
    {
        runtime_set_error("ozone forcing o3 units must be 'kg kg-1' or 'mol mol-1' in '%s'",
                          forcing.source);
        return 0;
    }
    if (!read_missing_value("_FillValue", type, &forcing.fill_value) ||
        !read_missing_value("missing_value", type, &forcing.missing_value))
    {
        runtime_set_error("invalid ozone forcing missing-value attributes in '%s'", forcing.source);
        return 0;
    }
    return 1;
}

int runtime_ozone_forcing_open(const char *path)
{
    int status = PLASIC_RUNTIME_IO_FAILED;

    runtime_ozone_forcing_close();
    if (path == NULL || path[0] == '\0')
        return PLASIC_RUNTIME_OK;
    if (strlen(path) >= sizeof(forcing.source))
    {
        runtime_set_error("ozone forcing NetCDF path is too long");
        return PLASIC_RUNTIME_BAD_ARGUMENT;
    }

    (void)snprintf(forcing.source, sizeof(forcing.source), "%s", path);
    const int open_status = nc_open(path, NC_NOWRITE, &forcing.file);
    if (open_status != NC_NOERR)
    {
        forcing.file = -1;
        runtime_set_error("cannot open ozone forcing '%s': %s", path, nc_strerror(open_status));
    }
    else if (read_header())
    {
        forcing.slabs = malloc(2 * local_values * sizeof(*forcing.slabs));
        if (forcing.slabs == NULL)
        {
            runtime_set_error("could not allocate ozone forcing slab");
            status = PLASIC_RUNTIME_INITIALIZATION_FAILED;
        }
        else
            status = PLASIC_RUNTIME_OK;
    }
    status = runtime_sync_error(status);
    if (status != PLASIC_RUNTIME_OK)
        runtime_ozone_forcing_close();
    return status;
}

static int read_month(int64_t record, double *slab)
{
    const size_t start[] = {(size_t)record, 0, (size_t)mp_local_offset(PLASIC_NLPP), 0};
    const size_t count[] = {1, PLASIC_NLEV, PLASIC_NLPP, PLASIC_NLON};
    int status;

    const int64_t date = forcing.first_month + record;
    const int year = (int)(date / 12);
    const int month = (int)(date % 12) + 1;

    status = nc_get_vara_double(forcing.file, forcing.variable, start, count, slab);
    if (status != NC_NOERR)
    {
        runtime_set_error("could not read ozone forcing for %d-%02d (%s) in '%s'",
                          year, month, nc_strerror(status), forcing.source);
        return PLASIC_RUNTIME_IO_FAILED;
    }
    for (size_t i = 0; i < local_values; ++i)
    {
        const double value = slab[i];
        const double mass_fraction = value * forcing.mass_factor;
        if (!isfinite(value) || value == forcing.fill_value ||
            value == forcing.missing_value || mass_fraction < 0.0 || mass_fraction > 1.0)
        {
            runtime_set_error("invalid or missing ozone forcing value at %d-%02d in '%s'",
                              year, month, forcing.source);
            return PLASIC_RUNTIME_IO_FAILED;
        }
        slab[i] = mass_fraction;
    }
    return PLASIC_RUNTIME_OK;
}

/* Consecutive records occupy different slots, including across December. */
static int cached_month(int64_t record, const double **values)
{
    const int slot = (int)(record % 2);
    double *slab = forcing.slabs + (size_t)slot * local_values;
    if (forcing.cached_record[slot] != record)
    {
        forcing.cached_record[slot] = -1;
        const int status = runtime_sync_error(read_month(record, slab));
        if (status != PLASIC_RUNTIME_OK)
            return status;
        forcing.cached_record[slot] = record;
    }
    *values = slab;
    return PLASIC_RUNTIME_OK;
}

int runtime_ozone_forcing_apply(void)
{
    const int year = plasic_ndatim[0], month = plasic_ndatim[1];
    const int64_t record = (int64_t)year * 12 + month - 1 - forcing.first_month;
    const double *current, *neighbor;
    int32_t month_a, month_b;
    float weight;

    if (forcing.file < 0)
        return PLASIC_RUNTIME_OK;
    if (record < 0 || record >= (int64_t)forcing.months)
    {
        runtime_set_error("ozone forcing '%s' does not cover model date %d-%02d",
                          forcing.source, year, month);
        return PLASIC_RUNTIME_BAD_ARGUMENT;
    }

    int status = cached_month(record, &current);
    if (status != PLASIC_RUNTIME_OK)
        return status;
    momint(plasic_calendar, plasic_nstep, &month_a, &month_b, &weight);
    const int64_t adjacent = record + (month_b - month_a);
    neighbor = current;
    if (weight > 0 && adjacent >= 0 && adjacent < (int64_t)forcing.months)
    {
        status = cached_month(adjacent, &neighbor);
        if (status != PLASIC_RUNTIME_OK)
            return status;
    }
    plasic_physics_interpolate_ozone(current, neighbor, weight);
    return PLASIC_RUNTIME_OK;
}
