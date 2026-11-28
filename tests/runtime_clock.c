/* Check the closing restart step at initialization and for later run intervals. */
#include "plasic_runtime.h"
#include "state/plasicmod_state.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    plasic_runtime_options options = {
        .restart_input = "data/earth_t21_l10.restart",
        .surface_data_directory = "data/surface",
        .co2_ppm = 300.0f,
        .calendar = -1,
        .start_yday = 1,
        .steps = 1,
        .reset_nstep = INT32_MAX,
        .dynamics_only = 1,
        .sea_coupling_interval = 1
    };

    if (plasic_runtime_initialize(&options) != PLASIC_RUNTIME_BAD_ARGUMENT)
    {
        fprintf(stderr, "Initialization accepted a run beyond INT32_MAX\n");
        return EXIT_FAILURE;
    }
    options.reset_nstep = INT32_MAX - 1;
    if (plasic_runtime_initialize(&options) != PLASIC_RUNTIME_OK)
    {
        fprintf(stderr, "Initialization at INT32_MAX - 1 failed: %s\n", plasic_runtime_error());
        return EXIT_FAILURE;
    }

    int result = EXIT_FAILURE;
    options.steps = 2;
    if (plasic_runtime_run(&options) != PLASIC_RUNTIME_BAD_ARGUMENT || plasic_nstep != INT32_MAX - 1)
        goto finish;
    options.steps = 1;
    if (plasic_runtime_run(&options) == PLASIC_RUNTIME_OK && plasic_nstep == INT32_MAX)
        result = EXIT_SUCCESS;

finish:
    if (result != EXIT_SUCCESS)
        fprintf(stderr, "Clock check failed at nstep=%d: %s\n", plasic_nstep, plasic_runtime_error());
    plasic_runtime_finalize();
    return result;
}
