#ifndef PLASIC_CONFIG_H
#define PLASIC_CONFIG_H

#include "plasic_runtime.h"

#include <stddef.h>

enum
{
    PLASIC_CONFIG_PATH_BYTES = 4096
};

typedef struct plasic_config_storage
{
    char restart_input[PLASIC_CONFIG_PATH_BYTES];
    char surface_data_directory[PLASIC_CONFIG_PATH_BYTES];
    char restart_output[PLASIC_CONFIG_PATH_BYTES];
    char progress_output[PLASIC_CONFIG_PATH_BYTES];
    char frame_output[PLASIC_CONFIG_PATH_BYTES];
    char co2_forcing_input[PLASIC_CONFIG_PATH_BYTES];
    char ozone_forcing_input[PLASIC_CONFIG_PATH_BYTES];
    char monthly_output_directory[PLASIC_CONFIG_PATH_BYTES];
} plasic_config_storage;

int plasic_config_load(
    const char *path, plasic_runtime_options *options,
    plasic_config_storage *storage, char *error, size_t error_bytes);

#endif
