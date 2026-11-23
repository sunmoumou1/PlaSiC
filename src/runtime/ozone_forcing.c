/*
 * 外部臭氧按月提供在模式高斯网格与 sigma 层上。每个 MPI 进程仅读取本地纬度带。
 * Monthly ozone is supplied on the native Gaussian and sigma grid. Each MPI
 * rank reads its latitude slab and caches the two bracketing monthly fields.
 * 辐射使用质量混合比（kg/kg），在相邻月的月中值之间线性插值。
 * Radiation takes mass mixing ratios (kg/kg), linearly interpolated between
 * monthly anchors. The first and last anchors extend to their calendar edges.
 */
#include "ozone_forcing.h"

#include "runtime_internal.h"

#include <ctype.h>
#include <float.h>
#include <netcdf.h>

static struct
{
    int file;
    int variable;
    size_t months;
    int64_t first_month;
    int64_t cached_record[2];
    int64_t validated_low, validated_high;
    double mass_factor;
    double fill_value;
    double missing_value;
    double *slabs;
    char source[4096];
} forcing = {.file = -1, .cached_record = {-1, -1}};

/* Molar masses in g/mol. */
static const double ozone_molar_mass = 47.9982, dry_air_molar_mass = 28.9644;
static const size_t local_values = (size_t)PLASIC_NHOR * PLASIC_NLEV;

/* Mass fractions down to this far below zero are regridding noise and read as zero. */
static const double negative_tolerance = 1e-9;

void runtime_ozone_forcing_close(void)
{
    if (forcing.file >= 0)
        (void)nc_close(forcing.file);
    free(forcing.slabs);
    forcing.slabs = NULL;
    forcing.file = -1;
    forcing.cached_record[0] = forcing.cached_record[1] = -1;
    forcing.validated_high = -1;
    forcing.source[0] = '\0';
}

const char *runtime_ozone_forcing_source(void)
{
    return forcing.file >= 0 ? forcing.source : "synthetic";
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

/* Parse a ratio or a product with the inverse of the same unit. */
static double units_factor(const char *units)
{
    char compact[32];
    size_t length = 0;

    for (; *units != '\0'; ++units)
    {
        if (isspace((unsigned char)*units))
            continue;
        if (length + 1 >= sizeof(compact))
            return 0.0;
        compact[length++] = *units;
    }
    compact[length] = '\0';
    const char *unit;
    double factor;
    if (strncmp(compact, "kg", 2) == 0)
    {
        unit = "kg";
        factor = 1.0;
    }
    else if (strncmp(compact, "mol", 3) == 0)
    {
        unit = strncmp(compact, "mole", 4) == 0 ? "mole" : "mol";
        factor = ozone_molar_mass / dry_air_molar_mass;
    }
    else
        return 0.0;

    const size_t unit_length = strlen(unit);
    const char *rest = compact + unit_length;
    if (*rest == '/')
        return strcmp(rest + 1, unit) == 0 ? factor : 0.0;
    if (*rest == '*')
        ++rest;
    if (strncmp(rest, unit, unit_length) != 0)
        return 0.0;
    rest += unit_length;
    return strcmp(rest, "-1") == 0 || strcmp(rest, "^-1") == 0 ||
           strcmp(rest, "**-1") == 0 ? factor : 0.0;
}

static int read_units(void)
{
    char *text = NULL;
    nc_type type;
    size_t length;

    if (nc_inq_att(forcing.file, forcing.variable, "units", &type, &length) != NC_NOERR)
        return 0;
    forcing.mass_factor = 0.0;
    if (type == NC_CHAR && length < SIZE_MAX)
    {
        /* Fortran writers can space-pad a fixed-width CHAR attribute. */
        text = malloc(length + 1);
        if (text != NULL && nc_get_att_text(forcing.file, forcing.variable, "units", text) == NC_NOERR)
        {
            while (length && (text[length - 1] == '\0' || isspace((unsigned char)text[length - 1])))
                --length;
            text[length] = '\0';
            if (memchr(text, '\0', length) == NULL)
                forcing.mass_factor = units_factor(text);
        }
        free(text);
    }
    else if (type == NC_STRING && length == 1 &&
             nc_get_att_string(forcing.file, forcing.variable, "units", &text) == NC_NOERR)
    {
        if (text != NULL)
            forcing.mass_factor = units_factor(text);
        (void)nc_free_string(1, &text);
    }
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
        runtime_set_error("ozone forcing o3 units must be 'kg kg-1' or 'mol mol-1' "
                          "(also written 'kg/kg', 'kg kg**-1', 'mol/mol', 'mole mole-1', 'mole/mole' ...) in '%s'",
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

static int sync_ozone_error(int32_t status)
{
    if (mp_sync_error(&status, runtime_error_message, sizeof(runtime_error_message)) != 0)
    {
        runtime_set_error("could not synchronize ozone forcing error");
        return PLASIC_RUNTIME_NUMERICAL_FAILURE;
    }
    return status;
}

static int open_local_file(const char *path)
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
    return status;
}

/* All ranks must take the same branches when loading and caching records. */
static int check_shared_time_axis(void)
{
    int32_t local[3] = {0, 0, 0};
    if (forcing.file >= 0)
    {
        local[0] = (int32_t)(forcing.first_month / 12);
        local[1] = (int32_t)(forcing.first_month % 12) + 1;
        local[2] = (int32_t)forcing.months;
    }
    int32_t root[3];
    memcpy(root, local, sizeof(root));
    int status = PLASIC_RUNTIME_OK;
    if (mp_broadcast_integer(root, 3) != 0)
    {
        runtime_set_error("could not broadcast ozone forcing time axis");
        status = PLASIC_RUNTIME_IO_FAILED;
    }
    else if (memcmp(local, root, sizeof(root)) != 0)
    {
        runtime_set_error("ozone forcing time axis differs on rank %d: %d-%02d (%d records), "
                          "rank 0 has %d-%02d (%d records)",
                          mp_rank(), local[0], local[1], local[2], root[0], root[1], root[2]);
        status = PLASIC_RUNTIME_BAD_ARGUMENT;
    }
    return sync_ozone_error(status);
}

int runtime_ozone_forcing_open(const char *path)
{
    int status = sync_ozone_error(open_local_file(path));
    if (status == PLASIC_RUNTIME_OK)
        status = check_shared_time_axis();
    if (status != PLASIC_RUNTIME_OK)
        runtime_ozone_forcing_close();
    return status;
}

static int read_month(int64_t record, double *slab, float warnings[2])
{
    const size_t latitude_offset = (size_t)mp_local_offset(PLASIC_NLPP);
    const size_t start[] = {(size_t)record, 0, latitude_offset, 0};
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
            value == forcing.missing_value || mass_fraction < -negative_tolerance || mass_fraction > 1.0)
        {
            const size_t horizontal = i % PLASIC_NHOR;
            runtime_set_error("invalid or missing ozone forcing value at %d-%02d: "
                              "o3(time=%lld, lev=%zu, lat=%zu, lon=%zu) = %g "
                              "(0-based indices) in '%s'",
                              year, month, (long long)record, i / PLASIC_NHOR,
                              latitude_offset + horizontal / PLASIC_NLON,
                              horizontal % PLASIC_NLON, value, forcing.source);
            return PLASIC_RUNTIME_IO_FAILED;
        }
        if (warnings != NULL)
        {
            if (mass_fraction > 1e-4)
                warnings[0] = 1.0f;
            if (mass_fraction < -1e-12)
                warnings[1] = 1.0f;
        }
        slab[i] = mass_fraction < 0.0 ? 0.0 : mass_fraction;
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
        const int status = sync_ozone_error(read_month(record, slab, NULL));
        if (status != PLASIC_RUNTIME_OK)
            return status;
        forcing.cached_record[slot] = record;
    }
    *values = slab;
    return PLASIC_RUNTIME_OK;
}

/*
 * File record of the month containing step, and the record its mid-month
 * interpolation reads with weight. Records may lie outside the file.
 */
static int64_t step_records(int32_t step, const int32_t date[7], int64_t *neighbor, float *weight)
{
    int32_t month_a, month_b;

    momint(plasic_calendar, step, &month_a, &month_b, weight);
    const int64_t record = (int64_t)date[0] * 12 + date[1] - 1 - forcing.first_month;
    *neighbor = *weight > 0 ? record + (month_b - month_a) : record;
    return record;
}

static int covered(int64_t record)
{
    return record >= 0 && record < (int64_t)forcing.months;
}

/*
 * Check that the file covers every month from the current step through the
 * last of steps, and read every record the run will interpolate from, so a
 * gap or bad value stops the run before its first step.
 */
int runtime_ozone_forcing_validate(int32_t steps)
{
    const int64_t last_step = (int64_t)plasic_nstep + steps - 1;
    int64_t first, last, first_neighbor, last_neighbor;
    int32_t end_date[7];
    float weight;

    if (forcing.file < 0)
        return PLASIC_RUNTIME_OK;
    step2cal(plasic_calendar, (int32_t)last_step, end_date);
    first = step_records(plasic_nstep, plasic_ndatim, &first_neighbor, &weight);
    last = step_records((int32_t)last_step, end_date, &last_neighbor, &weight);
    if (!covered(first) || !covered(last))
    {
        const int64_t file_end = forcing.first_month + (int64_t)forcing.months - 1;
        runtime_set_error("ozone forcing covers %d-%02d to %d-%02d, but the run needs "
                          "%d-%02d to %d-%02d (file '%s')",
                          (int)(forcing.first_month / 12), (int)(forcing.first_month % 12) + 1,
                          (int)(file_end / 12), (int)(file_end % 12) + 1,
                          plasic_ndatim[0], plasic_ndatim[1], end_date[0], end_date[1],
                          forcing.source);
        return PLASIC_RUNTIME_BAD_ARGUMENT;
    }

    /* Interpolation reaches one record beyond the run only at its two ends. */
    const int64_t low = first_neighbor < first && covered(first_neighbor) ? first_neighbor : first;
    const int64_t high = last_neighbor > last && covered(last_neighbor) ? last_neighbor : last;
    if (low >= forcing.validated_low && high <= forcing.validated_high)
        return PLASIC_RUNTIME_OK;
    int64_t retained[2] = {-1, -1};
    float warnings[2] = {0.0f, 0.0f};
    /* Reading can overwrite both cache slots. Publish their new identities
     * only when all ranks have succeeded, including on a caller's retry. */
    forcing.cached_record[0] = forcing.cached_record[1] = -1;
    forcing.validated_high = -1;
    /* Validate local slabs before synchronizing. Reverse order retains
     * the first interpolation pair in the cache when integration starts. */
    int status = PLASIC_RUNTIME_OK;
    for (int64_t record = high; record >= low && status == PLASIC_RUNTIME_OK; --record)
    {
        const int slot = (int)(record % 2);
        status = read_month(record, forcing.slabs + (size_t)slot * local_values, warnings);
        retained[slot] = record;
    }
    status = sync_ozone_error(status);
    if (status != PLASIC_RUNTIME_OK)
        return status;
    if (mp_allreduce_real(warnings, 2) != 0)
    {
        runtime_set_error("could not summarize ozone forcing range across ranks");
        return PLASIC_RUNTIME_NUMERICAL_FAILURE;
    }
    memcpy(forcing.cached_record, retained, sizeof(retained));
    forcing.validated_low = low;
    forcing.validated_high = high;
    if (mp_is_root() && warnings[0] > 0)
        fprintf(stderr, "External ozone exceeds 1e-4 kg/kg in the requested records; check units for Earth ozone.\n");
    if (mp_is_root() && warnings[1] > 0)
        fprintf(stderr, "External ozone below -1e-12 kg/kg was set to zero; check the input regridding.\n");
    return PLASIC_RUNTIME_OK;
}

int runtime_ozone_forcing_apply(void)
{
    const double *current, *adjacent;
    int64_t record, neighbor;
    float weight;

    if (forcing.file < 0)
        return PLASIC_RUNTIME_OK;
    record = step_records(plasic_nstep, plasic_ndatim, &neighbor, &weight);
    if (!covered(record))
    {
        runtime_set_error("ozone forcing '%s' does not cover model date %d-%02d",
                          forcing.source, plasic_ndatim[0], plasic_ndatim[1]);
        return PLASIC_RUNTIME_BAD_ARGUMENT;
    }

    int status = cached_month(record, &current);
    if (status != PLASIC_RUNTIME_OK)
        return status;
    adjacent = current;
    if (neighbor != record && covered(neighbor))
    {
        status = cached_month(neighbor, &adjacent);
        if (status != PLASIC_RUNTIME_OK)
            return status;
    }
    plasic_physics_interpolate_ozone(current, adjacent, weight);
    return PLASIC_RUNTIME_OK;
}
