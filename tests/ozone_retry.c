/* Exercise successive calls on one initialized MPI model, including an I/O failure. */
#include "plasic_runtime.h"
#include "state/plasicmod_state.h"

int main(int argc, char **argv)
{
    if (argc != 3)
        return 1;
    plasic_runtime_options options = {
        .ozone_forcing_input = argv[1], .surface_data_directory = argv[2],
        .co2_ppm = 280.0f, .steps = 1, .step_minutes = 45, .calendar = -1, .start_year = 1850,
        .start_yday = 1, .reset_nstep = -1, .sea_coupling_interval = 1,
        .frame_interval = 1, .progress_interval = 1
    };
    if (plasic_runtime_initialize(&options) != PLASIC_RUNTIME_OK ||
        plasic_runtime_run(&options) != PLASIC_RUNTIME_OK)
        return 1;
    const int32_t step = plasic_nstep;
    options.steps = 32 * 72;
    if (plasic_runtime_run(&options) != PLASIC_RUNTIME_IO_FAILED || plasic_nstep != step)
        return 1;
    options.steps = 1;
    const int status = plasic_runtime_run(&options);
    const int passed = status == PLASIC_RUNTIME_OK && plasic_nstep == step + 1;
    plasic_runtime_finalize();
    return passed ? 0 : 1;
}
