#ifndef PLASIC_RUNTIME_H
#define PLASIC_RUNTIME_H

#include "plasic_status.h"

#include <stdint.h>

typedef struct plasic_runtime_options
{
    const char *restart_input;
    const char *surface_data_directory;
    const char *restart_output;
    const char *progress_output;
    const char *frame_output;
    const char *co2_forcing_input;
    float co2_ppm;
    const char *monthly_output_directory;
    int32_t steps;
    int32_t frame_interval;
    int32_t progress_interval;
    /* Calendar kind: -1 = unspecified, 0 = Gregorian, 1 = 12 x 30-day (360day).
     * An unspecified calendar defaults to Gregorian on a cold start and to the
     * restart file's calendar on a warm start. */
    int32_t calendar;
    /* Cold-start date: year 0-based, yday 1-based, hour/minute 0-based. */
    int32_t start_year;
    int32_t start_yday;
    int32_t start_hour;
    int32_t start_minute;
    /* Warm-start clock rebase; negative means "not requested". */
    int32_t reset_nstep;
    int32_t dynamics_only;
    /*
     * Atmosphere-sea ice coupling interval in atmospheric steps (>= 1).
     * Sea-surface fluxes are accumulated every step; sea ice and the ocean
     * run once per interval and integrate the whole interval.
     */
    int32_t sea_coupling_interval;
    /*
     * Time step in minutes (--step-minutes, or step_minutes in the JSON
     * config); fractions such as 7.5 are allowed. Zero means "not set" and
     * selects the default of the resolution (PLASIC_DEFAULT_STEP_SECONDS).
     * A user step is used as given, without any stability check.
     */
    double step_minutes;
    /* Optional monthly ozone NetCDF on the model grid. */
    const char *ozone_forcing_input;
} plasic_runtime_options;

int plasic_runtime_initialize(const plasic_runtime_options *options);
/* Like initialize(), with a root-only error callback before MPI and I/O teardown. */
int plasic_runtime_initialize_with_diagnostic(
    const plasic_runtime_options *options, void (*on_error)(const char *message));
int plasic_runtime_run(const plasic_runtime_options *options);
int plasic_runtime_write_restart(const char *path);
void plasic_runtime_finalize(void);
const char *plasic_runtime_error(void);

#endif
