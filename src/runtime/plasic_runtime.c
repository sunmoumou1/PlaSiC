#include "runtime_internal.h"
#include "co2_forcing.h"
#include "ozone_forcing.h"
#include "monthly_netcdf.h"

runtime_workspace *runtime_fields;
FILE *progress_stream;
int progress_enabled;
int frames_enabled;
plasic_stream_writer frame_writer;
double integration_start_time;
char runtime_error_message[512];


/*
 * Set the runtime error message.
 */
void runtime_set_error(const char *format, ...)
{
    /*
     * va_list is used to access the variadic argument list.
     */
    va_list arguments;

    /*
     * Initialize the variadic argument list.
     */
    va_start(arguments, format);

    /*
     * Format the error message according to format and arguments, and write it
     * into runtime_error_message.
     *
     * Parameter meanings:
     *   runtime_error_message
     *       Target character buffer used to hold the final error string.
     *
     *   sizeof(runtime_error_message)
     *       Total capacity of the target buffer, in bytes.
     *       Passing the actual buffer size to vsnprintf() limits the maximum
     *       write length and prevents the formatted result from overflowing
     *       into other memory.
     *
     *   format
     *       printf-style format string.
     *
     *   arguments
     *       Variadic argument list initialized by va_start().
     *
     * The explicit cast to (void) indicates that the return value of
     * vsnprintf() is intentionally ignored. The return value of vsnprintf()
     * usually indicates "the number of characters that would have been
     * generated if the buffer space were sufficient (not including the
     * terminating '\0')", so it could also be used to detect truncation;
     * however, the current implementation does not perform that check and only
     * stores the error message that fits.
     */
    (void)vsnprintf(runtime_error_message, sizeof(runtime_error_message), format, arguments);

    /*
     * End access to the variadic argument list arguments.
     * Used in pairs with va_start(). After calling va_end(), arguments must
     * not be used again unless va_start() is executed again or it is
     * reinitialized in a valid way.
     */
    va_end(arguments);
}

/* Return the most recent runtime error message. */
const char *plasic_runtime_error(void)
{
    return runtime_error_message;
}

void plasic_runtime_finalize(void)
{
    /*
     * This function can be safely called if any step of initialization fails.
     */
    if (progress_stream != NULL)
    {
        (void)fclose(progress_stream);
        progress_stream = NULL;
    }
    progress_enabled = 0;

    if (frame_writer.stream != NULL)
    {
        (void)plasic_stream_close(&frame_writer);
    }
    frames_enabled = 0;

    runtime_monthly_netcdf_abort();
    runtime_co2_forcing_close();
    runtime_ozone_forcing_close();
    plasic_physics_finalize();
    sht_finalize();
    free(runtime_fields);
    runtime_fields = NULL;
    (void)mp_stop();
}

int runtime_sync_error(int32_t status)
{
    if (mp_sync_error(&status, runtime_error_message, sizeof(runtime_error_message)) != 0)
    {
        runtime_set_error("could not synchronize runtime error");
        return PLASIC_RUNTIME_NUMERICAL_FAILURE;
    }
    return status;
}

static int apply_monthly_forcing(void)
{
    const int year = plasic_ndatim[0];
    const int month = plasic_ndatim[1];
    const int status = runtime_co2_forcing_apply(year, month);
    return status == PLASIC_RUNTIME_OK ? runtime_ozone_forcing_apply() : status;
}

/* The final restart clock includes the increment after the last step. */
static int validate_run_length(int32_t steps)
{
    int status = PLASIC_RUNTIME_OK;
    if (steps <= 0 || (int64_t)plasic_nstep + steps > INT32_MAX)
    {
        runtime_set_error("cannot advance %d steps from nstep=%d within the model clock range",
                          steps, plasic_nstep);
        status = PLASIC_RUNTIME_BAD_ARGUMENT;
    }
    return runtime_sync_error(status);
}

static int initialize(const plasic_runtime_options *options)
{
    /*
     * Single entry point for the overall initialization of a PlaSiC run. 
     */
    int status; /* Uniformly receives the runtime status code returned by each initialization substep. */

    /*
     * Support repeated initialization within the same process: first close any
     * leftover from the previous run, then start from a clean state. The error message for this
     * run is also cleared here.
     */
    plasic_runtime_finalize();
    runtime_error_message[0] = '\0';

    {
        const int has_forcing =
            options->co2_forcing_input != NULL &&
            options->co2_forcing_input[0] != '\0';
        const int has_fixed =
            isfinite(options->co2_ppm) && options->co2_ppm > 0.0f;

        if (!isfinite(options->co2_ppm) || options->co2_ppm < 0.0f ||
            has_forcing + has_fixed != 1)
        {
            runtime_set_error(
                "specify exactly one CO2 mode: --co2-forcing FILE or "
                "--co2-ppm PPMV");
            return PLASIC_RUNTIME_BAD_ARGUMENT;
        }
    }

    /*
     * 时间选项的防御性校验：CLI/JSON 已经检查过，这里保证直接调用运行时
     * 接口的程序也能得到同样的行为。
     *
     * Defensive validation of the time options: the CLI/JSON layer checks
     * these too, but direct runtime callers must get the same behaviour.
     */
    if (options->steps <= 0)
    {
        runtime_set_error("--steps must be a positive integer");
        return PLASIC_RUNTIME_BAD_ARGUMENT;
    }
    if (options->sea_coupling_interval < 1)
    {
        runtime_set_error(
            "--sea-coupling-interval must be a positive integer");
        return PLASIC_RUNTIME_BAD_ARGUMENT;
    }
    if (options->restart_input != NULL &&
        (options->start_year != 0 || options->start_yday != 1 ||
         options->start_hour != 0 || options->start_minute != 0))
    {
        runtime_set_error(
            "warm start reads its date from the restart file; start-date "
            "options only apply to cold starts");
        return PLASIC_RUNTIME_BAD_ARGUMENT;
    }
    if (options->restart_input == NULL && options->reset_nstep >= 0)
    {
        runtime_set_error("--reset-nstep requires a restart file");
        return PLASIC_RUNTIME_BAD_ARGUMENT;
    }

    {
        int32_t mpi_error;
        int32_t process_count;
        int mp_start_status;

        /*
         * plasic_mpi uniformly manages the serial and native MPI backends.
         * PLASIC_NPRO must match the actual number of processes, otherwise the
         * decomposition sizes of the grid and spectral arrays are no longer valid.
         */
        mp_start_status = mp_start(PLASIC_NPRO, &process_count, &mpi_error);

        if (mp_start_status != PLASIC_MP_START_OK)
        {
            runtime_set_error(
                "parallel initialization failed: expected %d processes, "
                "got %d (backend status %d)",
                PLASIC_NPRO, process_count, mpi_error);
            return PLASIC_RUNTIME_UNSUPPORTED_CONFIGURATION;
        }
    }

    /*
     * Initialize the global arrays and default scalars managed by
     * plasicmod_state
     */
    plasic_state_init();
    plasic_mypid = mp_rank();

    /*
     * Allocate the temporary workspace for this run.
     * calloc zeroing ensures that buffers not yet written do not carry leftover
     * values from the previous run.
     */
    runtime_fields = (runtime_workspace *)calloc(1, sizeof(*runtime_fields));
    status = PLASIC_RUNTIME_OK;
    if (runtime_fields == NULL)
    {
        runtime_set_error("could not allocate runtime fields");
        status = PLASIC_RUNTIME_INITIALIZATION_FAILED;
    }
    else
    {
        status = runtime_initialize_dynamics();
    }
    status = runtime_sync_error(status);
    if (status != PLASIC_RUNTIME_OK)
    {
        return status;
    }

    /*
     * 时间基准：
     *   冷启动：日历由 --calendar 选择，起始步由 --start-year/--start-yday/
     *           --start-hour/--start-minute 经 cal2step() 算出；
     *   热启动：日历与 nstep 都来自 restart 文件（见下方加载函数）。
     *
     * Time base:
     *   cold start: the calendar comes from --calendar and the starting step
     *               is computed from the --start-* options by cal2step();
     *   warm start: both the calendar and nstep come from the restart file.
     */
    if (options->restart_input == NULL)
    {
        /* An unspecified --calendar means Gregorian for a cold start. */
        const int32_t cold_calendar =
            options->calendar >= 0 ? options->calendar
                                   : PLASIC_CALENDAR_GREGORIAN;
        const int32_t calendar_status = cal2step(
            cold_calendar, options->start_year, options->start_yday,
            options->start_hour, options->start_minute, &plasic_nstep);

        plasic_calendar = cold_calendar;
        if (calendar_status != CALMOD_OK)
        {
            runtime_set_error(
                "invalid cold-start date (year %d, yday %d, %02d:%02d, %s "
                "calendar): %s",
                options->start_year, options->start_yday, options->start_hour,
                options->start_minute,
                calmod_calendar_name(cold_calendar),
                calmod_status_text(calendar_status));
            status = PLASIC_RUNTIME_BAD_ARGUMENT;
        }
        status = runtime_sync_error(status);
        if (status != PLASIC_RUNTIME_OK)
        {
            return status;
        }
    }

    status = runtime_load_atmosphere_restart(
        options->restart_input, options->surface_data_directory);
    status = runtime_sync_error(status);
    if (status != PLASIC_RUNTIME_OK)
    {
        return status;
    }

    /*
     * 热启动时日历以 restart 为准：--calendar 只是默认值，若与文件冲突
     * 给出提示并继续使用文件中的日历。
     *
     * On a warm start the restart file owns the calendar; a conflicting
     * --calendar option is reported and ignored.
     */
    if (options->restart_input != NULL && options->calendar >= 0 &&
        options->calendar != plasic_calendar && plasic_mypid == 0)
    {
        fprintf(
            stderr,
            "ignoring --calendar %s: restart file declares %s\n",
            calmod_calendar_name(options->calendar),
            calmod_calendar_name(plasic_calendar));
    }

    step2cal(plasic_calendar, plasic_nstep, plasic_ndatim);

    /*
     * --reset-nstep 只在热启动时生效：把 restart 中的 nstep 直接改成指定
     * 值（物理状态不变，只重贴日期标签），并在终端打印新的模型日期。
     *
     * --reset-nstep only applies to a warm start: it replaces the nstep read
     * from the restart file while the physical state stays as stored, and the
     * new model date is printed.
     */
    if (options->restart_input != NULL && options->reset_nstep >= 0)
    {
        const int32_t previous_nstep = plasic_nstep;

        plasic_nstep = options->reset_nstep;
        step2cal(plasic_calendar, plasic_nstep, plasic_ndatim);

        if (plasic_mypid == 0)
        {
            printf(
                "--reset-nstep: nstep %d -> %d, model date "
                "%d-%02d-%02d %02d:%02d (%s calendar)\n",
                previous_nstep, plasic_nstep, plasic_ndatim[0],
                plasic_ndatim[1], plasic_ndatim[2], plasic_ndatim[3],
                plasic_ndatim[4], calmod_calendar_name(plasic_calendar));
            fflush(stdout);
        }
    }

    status = validate_run_length(options->steps);
    if (status != PLASIC_RUNTIME_OK)
        return status;

    status = runtime_sync_error(runtime_prepare_orography_grid());
    if (status != PLASIC_RUNTIME_OK)
    {
        return status;
    }

    /*
     * Correction: the higher the mean orography, the lower the reference pressure.
     */
    plasic_psurf = expf(
        logf(plasic_psurf) -
        plasic_so[0] * plasic_cv * plasic_cv /
            sqrtf(2.0f) / (plasic_gascon * plasic_tgr));

    status = plasic_physics_initialize(
        options->restart_input, options->surface_data_directory,
        options->sea_coupling_interval);
    if (status != PLASIC_RUNTIME_OK)
    {
        runtime_set_error("%s", plasic_physics_error());
        return status;
    }

    status = runtime_co2_forcing_open(
        options->co2_forcing_input, options->co2_ppm);
    if (status == PLASIC_RUNTIME_OK)
    {
        status = runtime_ozone_forcing_open(options->ozone_forcing_input);
    }
    if (status == PLASIC_RUNTIME_OK)
    {
        status = runtime_ozone_forcing_validate(options->steps);
    }
    if (status == PLASIC_RUNTIME_OK)
    {
        status = apply_monthly_forcing();
    }
    if (status != PLASIC_RUNTIME_OK)
    {
        return status;
    }

    status = runtime_monthly_netcdf_open(
        options->monthly_output_directory);
    if (status != PLASIC_RUNTIME_OK)
    {
        return status;
    }

    return runtime_open_live_outputs(options);
}

int plasic_runtime_initialize_with_diagnostics(const plasic_runtime_options *options, FILE *stream)
{
    const int status = initialize(options);
    if (status != PLASIC_RUNTIME_OK)
    {
        /* The shutdown barrier keeps peers alive until root has reported the error. */
        if (stream != NULL && mp_is_root())
        {
            fprintf(stream, "PlaSiC initialization failed: %s\n", plasic_runtime_error());
            fflush(stream);
        }
        plasic_runtime_finalize();
    }
    return status;
}

int plasic_runtime_initialize(const plasic_runtime_options *options)
{
    return plasic_runtime_initialize_with_diagnostics(options, NULL);
}

// The most important function is the one below.
int plasic_runtime_run(const plasic_runtime_options *options)
{
    /*
     * PlaSiC's main time integration loop: each iteration advances one
     * dynamical time step Δτ.
     */
    int32_t step; /* Relative step number within this plasic_runtime_run() call, starting from 0. */
    int status;   /* Runtime status code uniformly returned by each stage. */

    /*
     * Advance step by step; options->steps is always >= 1 because the CLI and
     * the JSON loader reject a missing or zero step count.
     */
    status = validate_run_length(options->steps);
    if (status != PLASIC_RUNTIME_OK)
        return status;

    for (step = 0; step < options->steps; ++step)
    {
        step2cal(plasic_calendar, plasic_nstep, plasic_ndatim);

        status = apply_monthly_forcing();
        if (status != PLASIC_RUNTIME_OK)
        {
            return status;
        }

        /*
         * Stage 1
         */
        status = runtime_gridpoint_adiabatic();
        if (status != PLASIC_RUNTIME_OK)
        {
            return status;
        }

        /*
         * Stage 2
         */
        status = runtime_spectral_adiabatic();
        if (status != PLASIC_RUNTIME_OK)
        {
            return status;
        }

        /*
         * Stage 3
         */
        status = runtime_gridpoint_diabatic(options->dynamics_only);
        if (status != PLASIC_RUNTIME_OK)
        {
            return status;
        }

        /*
         * Stage 4
         */
        status = runtime_spectral_diabatic();
        if (status != PLASIC_RUNTIME_OK)
        {
            return status;
        }

        if (!runtime_state_is_finite())
        {
            runtime_set_error("non-finite atmospheric state after nstep=%d", plasic_nstep);
            return PLASIC_RUNTIME_NUMERICAL_FAILURE;
        }

        if (!options->dynamics_only)
        {
            runtime_accumulate_output_fields();
            status = runtime_monthly_netcdf_accumulate(
                plasic_ndatim[0], plasic_ndatim[1]);
            if (status != PLASIC_RUNTIME_OK)
            {
                return status;
            }
        }

        ++plasic_nstep;
        step2cal(plasic_calendar, plasic_nstep, plasic_ndatim);

        status = runtime_write_step_outputs(options, step + 1);
        if (status != PLASIC_RUNTIME_OK)
        {
            runtime_set_error("could not write live model output");
            return status;
        }
    }

    /*
     * After all requested time steps complete successfully
     */
    status = runtime_monthly_netcdf_close();
    if (status != PLASIC_RUNTIME_OK)
    {
        return status;
    }

    status = runtime_close_step_outputs(options);
    if (status != PLASIC_RUNTIME_OK)
    {
        runtime_set_error("could not close live model output");
        return status;
    }

    return PLASIC_RUNTIME_OK;
}

int runtime_broadcast_status(int status)
{
    /* In MPI runs, all ranks must make the same decision about success/failure to prevent some processes from hanging in later communication. */
    int32_t broadcast_status_value = (int32_t)status;

    if (mp_broadcast_integer(&broadcast_status_value, 1) != 0)
    {
        return PLASIC_RUNTIME_IO_FAILED;
    }
    return (int)broadcast_status_value;
}

int plasic_runtime_write_restart(const char *path)
{
    /*
     * Save the current model state as a restart file that can continue the
     * integration.
     */
    int status; /* Runtime status code of the current write stage. */

    /*
     * Only the root process opens the shared restart file; the result is then
     * broadcast so that all ranks exit together if file creation fails, rather
     * than some processes continuing into subsequent MPI communication.
     */
    status = PLASIC_RUNTIME_OK;
    if (mp_is_root() &&
        restart_open_write(path) != PLASIC_RESTART_OK)
    {
        runtime_set_error("cannot create restart '%s'", path);
        status = PLASIC_RUNTIME_IO_FAILED;
    }
    status = runtime_broadcast_status(status);
    if (status != PLASIC_RUNTIME_OK)
    {
        return status;
    }

    /*
     * The two macros uniformly handle the repeated "write record -> check
     * status -> jump on failure" flow:
     *   WRITE_INTEGER: write a single integer record;
     *   WRITE_ARRAY: write a complete array already held by the root process.
     * The underlying restart file writes column-contiguously, and the stored
     * row count describes the actual length of each column.
     */
#define WRITE_INTEGER(name, value)                         \
    do                                                     \
    {                                                      \
        status = restart_write_root_integer((name), (value)); \
        if (status != PLASIC_RESTART_OK)                   \
        {                                                  \
            goto failed;                                   \
        }                                                  \
    } while (0)
#define WRITE_ARRAY(name, values, rows, columns)             \
    do                                                       \
    {                                                        \
        status = restart_write_root_array(                   \
            (name), (values), (rows), (columns));            \
        if (status != PLASIC_RESTART_OK)                     \
        {                                                    \
            goto failed;                                     \
        }                                                    \
    } while (0)

    /* First write the scalar metadata needed for recovery; the reader uses it to check whether the grid and spectral configuration match. */
    WRITE_INTEGER("nstep", plasic_nstep);
    WRITE_INTEGER("calendar", plasic_calendar);
    WRITE_INTEGER("naccuout", plasic_naccuout);
    WRITE_INTEGER("nlat", PLASIC_NLAT);
    WRITE_INTEGER("nlon", PLASIC_NLON);
    WRITE_INTEGER("nlev", PLASIC_NLEV);
    WRITE_INTEGER("nrsp", PLASIC_NRSP);

    /* The random number seed is written by the root process, and the I/O result is synchronized to all ranks. */
    status = PLASIC_RUNTIME_OK;
    if (mp_is_root() &&
        restart_write_seed(
            "seed", plasic_seed,
            (int32_t)(sizeof(plasic_seed) /
                      sizeof(plasic_seed[0]))) != PLASIC_RESTART_OK)
    {
        runtime_set_error("restart seed could not be written");
        status = PLASIC_RUNTIME_IO_FAILED;
    }
    status = runtime_broadcast_status(status);
    if (status != PLASIC_RUNTIME_OK)
    {
        goto failed;
    }

    /*
     * Write the complete spectral arrays of the current time level. 
     */
    WRITE_ARRAY("sz", plasic_sz, PLASIC_NRSP, PLASIC_NLEV);
    WRITE_ARRAY("sd", plasic_sd, PLASIC_NRSP, PLASIC_NLEV);
    WRITE_ARRAY("st", plasic_st, PLASIC_NRSP, PLASIC_NLEV);
    WRITE_ARRAY("sq", plasic_sq, PLASIC_NRSP, PLASIC_NLEV);
    WRITE_ARRAY("sr", plasic_sr, PLASIC_NRSP, PLASIC_NLEV);
    WRITE_ARRAY("sp", plasic_sp, PLASIC_NRSP, 1);
    WRITE_ARRAY("so", plasic_so, PLASIC_NRSP, 1);

    /*
     * Write the history arrays of the leapfrog previous time level. The *m
     * arrays are stored in per-rank blocks, and the helper function first
     * gathers them into complete spectral arrays before the root process
     * writes them
     */
#define WRITE_PARTIAL(name, values, columns)                         \
    do                                                               \
    {                                                                \
        status = restart_write_distributed_spectral(                 \
            (name), (values), PLASIC_NRSP, PLASIC_NSPP, (columns));  \
        if (status != PLASIC_RESTART_OK)                             \
        {                                                            \
            goto failed;                                             \
        }                                                            \
    } while (0)
    WRITE_PARTIAL("szm", plasic_szm, PLASIC_NLEV);
    WRITE_PARTIAL("sdm", plasic_sdm, PLASIC_NLEV);
    WRITE_PARTIAL("stm", plasic_stm, PLASIC_NLEV);
    WRITE_PARTIAL("sqm", plasic_sqm, PLASIC_NLEV);
    WRITE_PARTIAL("spm", plasic_spm, 1);
#undef WRITE_PARTIAL

    /*
     * Write local grid-point physical quantities and accumulators. Each rank
     * only holds its own latitudinal band, and the helper function assembles
     * them into a complete global grid before writing; columns=1 denotes a
     * single-level field and columns=NLEV denotes a three-dimensional field.
     */
#define WRITE_GRID(name, values, columns)                         \
    do                                                            \
    {                                                             \
        status = restart_write_distributed_grid(                  \
            (name), (values), PLASIC_NUGP, PLASIC_NHOR, (columns)); \
        if (status != PLASIC_RESTART_OK)                          \
        {                                                         \
            goto failed;                                          \
        }                                                         \
    } while (0)
    WRITE_GRID("dls", plasic_dls, 1);
    WRITE_GRID("dwetfac", plasic_dwetfac, 1);
    WRITE_GRID("dalb", plasic_dalb, 1);
    WRITE_GRID("dz0", plasic_dz0, 1);
    WRITE_GRID("dicec", plasic_dicec, 1);
    WRITE_GRID("diced", plasic_diced, 1);
    WRITE_GRID("dwatc", plasic_dwatc, 1);
    WRITE_GRID("dust3", plasic_dust3, 1);
    WRITE_GRID("dcc", plasic_dcc, PLASIC_NLEV);
    WRITE_GRID("dql", plasic_dql, PLASIC_NLEV);
    WRITE_GRID("dqsat", plasic_dqsat, PLASIC_NLEV);
    WRITE_GRID("dt", plasic_dt + PLASIC_NHOR * PLASIC_NLEV, 1);
    WRITE_GRID("dq", plasic_dq + PLASIC_NHOR * PLASIC_NLEV, 1);
    WRITE_GRID("aprl", plasic_aprl, 1);
    WRITE_GRID("aprc", plasic_aprc, 1);
    WRITE_GRID("aprs", plasic_aprs, 1);
    WRITE_GRID("aevap", plasic_aevap, 1);
    WRITE_GRID("ashfl", plasic_ashfl, 1);
    WRITE_GRID("alhfl", plasic_alhfl, 1);
    WRITE_GRID("asmelt", plasic_asmelt, 1);
    WRITE_GRID("asndch", plasic_asndch, 1);
    WRITE_GRID("acc", plasic_acc, 1);
    WRITE_GRID("assol", plasic_assol, 1);
    WRITE_GRID("asthr", plasic_asthr, 1);
    WRITE_GRID("atsol", plasic_atsol, 1);
    WRITE_GRID("atthr", plasic_atthr, 1);
    WRITE_GRID("ataux", plasic_ataux, 1);
    WRITE_GRID("atauy", plasic_atauy, 1);
    WRITE_GRID("atsolu", plasic_atsolu, 1);
    WRITE_GRID("assolu", plasic_assolu, 1);
    WRITE_GRID("asthru", plasic_asthru, 1);
    WRITE_GRID("aqvi", plasic_aqvi, 1);
    WRITE_GRID("atsa", plasic_atsa, 1);
    WRITE_GRID("ats0", plasic_ats0, 1);
    WRITE_GRID("atsama", plasic_atsama, 1);
    WRITE_GRID("atsami", plasic_atsami, 1);
#undef WRITE_GRID

    /* The physics module appends its own land surface, ocean, sea ice and vegetation restart records. */
    status = plasic_physics_write_restart();
    if (status != PLASIC_RUNTIME_OK)
    {
        runtime_set_error("component restart records could not be written");
        goto failed;
    }

    /* After all records are written successfully, close the file and synchronize the close result to all ranks. */
    status = PLASIC_RUNTIME_OK;
    if (mp_is_root() &&
        restart_close() != PLASIC_RESTART_OK)
    {
        runtime_set_error("cannot close restart output '%s'", path);
        status = PLASIC_RUNTIME_IO_FAILED;
    }
    return runtime_broadcast_status(status);

failed:
    /* If any record write fails, close the already opened file to avoid leaving an unclosed file handle. */
    if (mp_is_root())
    {
        (void)restart_close();
    }
    return PLASIC_RUNTIME_IO_FAILED;

#undef WRITE_ARRAY
#undef WRITE_INTEGER
}
