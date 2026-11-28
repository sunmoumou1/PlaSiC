/*
 * JSON configuration reader for the PlaSiC model.
 *
 * Values are read from a flat JSON object with the vendored cJSON library;
 * unknown keys are ignored and missing keys keep their current value.
 */

#include "plasic_config.h"

#include "tools/calmod_kernels.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../third_party/cJSON/cJSON.h"

/* Upper bound on the accepted configuration file size (2 MiB). */
#define PLASIC_CONFIG_MAX_BYTES (2 * 1024 * 1024)

/* Read the whole file into a NUL-terminated buffer; NULL on failure. */
static char *read_file(const char *path)
{
    FILE *stream = fopen(path, "rb");
    long size;
    char *text;

    if (stream == NULL)
        return NULL;
    if (fseek(stream, 0, SEEK_END) != 0 || (size = ftell(stream)) < 0 ||
        size > PLASIC_CONFIG_MAX_BYTES || fseek(stream, 0, SEEK_SET) != 0)
    {
        fclose(stream);
        return NULL;
    }
    text = malloc((size_t)size + 1);
    if (text != NULL && fread(text, 1, (size_t)size, stream) == (size_t)size)
        text[size] = '\0';
    else
    {
        free(text);
        text = NULL;
    }
    fclose(stream);
    return text;
}

/* String option; nullable options treat JSON null or "" as "not set". */
static int load_string(
    const cJSON *root, const char *key, char *storage, size_t bytes,
    const char **option, int nullable)
{
    /* Look up the requested key in the JSON object. */
    const cJSON *item = cJSON_GetObjectItemCaseSensitive(root, key);
    size_t length;

    /* Missing keys are allowed and leave the option unchanged. */
    if (item == NULL)
        return 1;

    /* JSON null is accepted only for nullable options. */
    if (cJSON_IsNull(item))
    {
        if (!nullable)
            return 0;

        *option = NULL;
        return 1;
    }

    /* The value must be a JSON string. */
    if (!cJSON_IsString(item))
        return 0;

    length = strlen(item->valuestring);

    /* Ensure the string, including '\0', fits in the destination buffer. */
    if (length >= bytes)
        return 0;

    memcpy(storage, item->valuestring, length + 1);

    /* For nullable options, an empty string is treated as "not set". */
    *option = (nullable && length == 0) ? NULL : storage;

    return 1;
}

/* Integer option; rejects non-integral and out-of-range numbers. */
static int load_integer(const cJSON *root, const char *key, int32_t *option)
{
    const cJSON *item = cJSON_GetObjectItemCaseSensitive(root, key);
    double value;

    if (item == NULL)
        return 1;
    if (!cJSON_IsNumber(item))
        return 0;
    value = item->valuedouble;
    if (value < (double)INT32_MIN || value > (double)INT32_MAX ||
        value != (double)(int32_t)value)
        return 0;
    *option = (int32_t)value;
    return 1;
}

/* Floating-point option. */
static int load_float(const cJSON *root, const char *key, float *option)
{
    const cJSON *item = cJSON_GetObjectItemCaseSensitive(root, key);

    if (item == NULL)
        return 1;
    if (!cJSON_IsNumber(item))
        return 0;
    *option = (float)item->valuedouble;
    return 1;
}

/*
 * Time step in minutes: a positive number (fractions allowed), or null for
 * the default of the resolution. The value is not checked for stability.
 */
static int load_step_minutes(const cJSON *root, const char *key, double *option)
{
    const cJSON *item = cJSON_GetObjectItemCaseSensitive(root, key);

    if (item == NULL)
        return 1;
    if (cJSON_IsNull(item))
    {
        *option = 0.0;
        return 1;
    }
    if (!cJSON_IsNumber(item) || !(item->valuedouble > 0.0))
        return 0;
    *option = item->valuedouble;
    return 1;
}

/* Boolean option (JSON true/false). */
static int load_boolean(const cJSON *root, const char *key, int32_t *option)
{
    const cJSON *item = cJSON_GetObjectItemCaseSensitive(root, key);

    if (item == NULL)
        return 1;
    if (!cJSON_IsBool(item))
        return 0;
    *option = cJSON_IsTrue(item) ? 1 : 0;
    return 1;
}

/* Calendar option: "gregorian" or "360day". */
static int load_calendar(const cJSON *root, const char *key, int32_t *option)
{
    const cJSON *item = cJSON_GetObjectItemCaseSensitive(root, key);
    int32_t calendar;

    if (item == NULL)
        return 1;
    if (!cJSON_IsString(item))
        return 0;
    calendar = calmod_calendar_from_name(item->valuestring);
    if (calendar < 0)
        return 0;
    *option = calendar;
    return 1;
}

/*
 * Load runtime options from the JSON file at `path`. String options are
 * copied into the caller-provided `storage` buffers, which must outlive
 * `options`. Returns 1 on success, 0 on failure with a message in `error`.
 */
int plasic_config_load(
    const char *path, plasic_runtime_options *options,
    plasic_config_storage *storage, char *error, size_t error_bytes)
{
    char *text = read_file(path);
    const char *where;
    cJSON *root;
    int valid;

    root = cJSON_Parse(text);
    if (root == NULL)
    {
        where = cJSON_GetErrorPtr();
        snprintf(
            error, error_bytes, "config is not valid JSON near byte %lld",
            where != NULL ? (long long)(where - text) : 0LL);
        free(text);
        return 0;
    }
    free(text);
    if (!cJSON_IsObject(root))
    {
        cJSON_Delete(root);
        snprintf(error, error_bytes, "config root must be a JSON object");
        return 0;
    }
    /*
     * Every load must succeed; short-circuit evaluation stops at the first
     * invalid value.
     */
#define LOAD_STRING(key, field, nullable)                                    \
    load_string(                                                             \
        root, key, storage->field, sizeof(storage->field), &options->field,  \
        nullable)
    valid =
        LOAD_STRING("restart", restart_input, 1) &&
        LOAD_STRING("data_dir", surface_data_directory, 0) &&
        LOAD_STRING("output", restart_output, 0) &&
        LOAD_STRING("progress", progress_output, 1) &&
        LOAD_STRING("frames", frame_output, 1) &&
        LOAD_STRING("co2_forcing", co2_forcing_input, 1) &&
        LOAD_STRING("ozone_forcing", ozone_forcing_input, 1) &&
        LOAD_STRING("monthly_output", monthly_output_directory, 1) &&
        load_float(root, "co2_ppm", &options->co2_ppm) &&
        load_integer(root, "steps", &options->steps) &&
        load_integer(root, "frame_interval", &options->frame_interval) &&
        load_integer(root, "progress_interval", &options->progress_interval) &&
        load_integer(
            root, "sea_coupling_interval", &options->sea_coupling_interval) &&
        load_integer(root, "start_year", &options->start_year) &&
        load_integer(root, "start_yday", &options->start_yday) &&
        load_integer(root, "start_hour", &options->start_hour) &&
        load_integer(root, "start_minute", &options->start_minute) &&
        load_integer(root, "reset_nstep", &options->reset_nstep) &&
        load_step_minutes(root, "step_minutes", &options->step_minutes) &&
        load_calendar(root, "calendar", &options->calendar) &&
        load_boolean(root, "dynamics_only", &options->dynamics_only);
#undef LOAD_STRING
    cJSON_Delete(root);
    if (!valid)
    {
        snprintf(
            error, error_bytes,
            "config contains an invalid or overlong value");
        return 0;
    }
    return 1;
}
