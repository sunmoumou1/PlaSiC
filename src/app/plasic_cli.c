/*
 * PlaSiC standalone C program entry point.
 *
 * Note:
 * This file itself is mainly responsible for program startup and scheduling;
 * the actual PlaSiC initialization, integration, and restart I/O logic
 * are implemented by functions such as plasic_runtime_*.
 */

#include "plasic_runtime.h"

#include "plasic_cli.h"
#include "plasic_config.h"
#include "state/plasicmod_state.h"
#include "tools/calmod_kernels.h"

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Default surface data directory.
 */
#define PLASIC_DEFAULT_SURFACE_DATA "data/surface"

/*
 * Print command-line usage.
 */
static void print_usage(FILE *stream, const char *program)
{
    fprintf(stream,
            "Usage: %s [--config FILE] [--restart FILE] [--data-dir DIR]\n"
            "          --steps N [--output FILE]\n"
            "          [--progress FILE] [--frames FILE]\n"
            "          [--frame-interval N] [--progress-interval N]\n"
            "          [--monthly-output DIR]\n"
            "          [--co2-forcing FILE | --co2-ppm PPMV] [--ozone-forcing FILE]\n"
            "          [--calendar gregorian|360day]\n"
            "          [--start-year YEAR] [--start-yday DAY]\n"
            "          [--start-hour HOUR] [--start-minute MINUTE]\n"
            "          [--reset-nstep N]\n"
            "          [--sea-coupling-interval N] [--dynamics-only]\n\n"
            "Standalone C PlaSiC (T%d, %d latitude rows, %d levels).\n"
            "--steps is mandatory. A cold start (no --restart) begins at\n"
            "year --start-year, day-of-year --start-yday, --start-hour:--start-minute\n"
            "(defaults: year 0, day 1, 00:00). A warm start takes nstep and the\n"
            "calendar from the restart file; --reset-nstep rebases nstep.\n"
            "For every resolution, omitting --restart performs a cold start\n"
            "from NetCDF surface data under --data-dir. Restart files are\n"
            "used only when explicitly selected by CLI or JSON config.\n"
            "--sea-coupling-interval N runs sea ice and the ocean once every\n"
            "N atmospheric steps (default 1); sea-surface fluxes are still\n"
            "averaged over every step of the interval.\n",
            program, PLASIC_NTRU, PLASIC_NLAT, PLASIC_NLEV);
}

/*
 * Parse a string as a non-negative integer within the int32_t range.
 *
 * Parameters:
 *   text  - string to parse, e.g. "100"
 *   value - result written upon successful parsing
 *
 * Returns:
 *   1 - parsing succeeded
 *   0 - parsing failed
 *
 * Currently allowed numeric range:
 *
 *     0 <= value <= INT32_MAX
 */
static int parse_positive_integer(const char *text, int32_t *value)
{
    char *end = NULL;
    long parsed;

    /*
     * Basic parameter checks:
     */
    if (text == NULL || value == NULL || text[0] == '\0')
    {
        return 0;
    }

    /*
     * strtol() may set errno on range errors.
     *
     * errno must be cleared before the call;
     * otherwise it is impossible to distinguish an error produced by
     * the current strtol call from a leftover errno.
     */
    errno = 0;

    /*
     * Parse in base 10.【按十进制解析】
     *
     * end will ultimately point to the first character that was not parsed.
     */
    parsed = strtol(text, &end, 10);

    /*
     * Any of the following conditions makes the input invalid:
     */
    if (errno != 0 || end == text || *end != '\0' ||
        parsed < 0 || parsed > INT32_MAX)
    {
        return 0;
    }

    /*
     * At this point it is safe to convert to int32_t.
     */
    *value = (int32_t)parsed;

    return 1;
}

static int parse_positive_float(const char *text, float *value)
{
    char *end = NULL;
    float parsed;

    if (text == NULL || value == NULL || text[0] == '\0')
    {
        return 0;
    }
    errno = 0;
    parsed = strtof(text, &end);
    if (errno != 0 || end == text || *end != '\0' ||
        !isfinite(parsed) || parsed <= 0.0f)
    {
        return 0;
    }
    *value = parsed;
    return 1;
}

// The most important function in this document
int plasic_cli_run(int argc, char **argv)
{
    plasic_runtime_options options = {
        .restart_input = NULL,

        .surface_data_directory = PLASIC_DEFAULT_SURFACE_DATA,

        .restart_output = "status",

        .progress_output = NULL,

        .frame_output = NULL,

        .co2_forcing_input = NULL,

        .co2_ppm = 0.0f,

        .monthly_output_directory = NULL,

        /*
         * --steps has no default: an explicit positive value is mandatory.
         * Zero marks "not provided" until the argument scan finishes.
         */
        .steps = 0,

        .frame_interval = 4,

        .progress_interval = 1,

        /* Cold start defaults to year 0, day 1, 00:00 in the Gregorian calendar. */
        .calendar = -1,

        .start_year = 0,

        .start_yday = 1,

        .start_hour = 0,

        .start_minute = 0,

        /* Negative reset_nstep means "keep the nstep stored in the restart". */
        .reset_nstep = -1,

        .dynamics_only = 0,

        /* Sea ice and the ocean run every atmospheric step by default. */
        .sea_coupling_interval = 1};

    plasic_config_storage config_storage = {0};

    char config_error[256];

    const char *config_path = NULL;

    int argument;
    int status;

    /*
     * ================================================================
     * First CLI scan: look only for --config
     * ================================================================
     *
     * The reason for scanning --config separately first
     * is to implement the following priority:
     *
     *     defaults
     *       ↓
     *     config file
     *       ↓
     *     command line
     *
     * That is:
     *
     * First use the config file to modify options,
     * then fully parse the CLI,
     * and values explicitly given on the CLI override those in the config file.
     */
    for (argument = 1; argument < argc; ++argument)
    {
        /*
         * Find:
         *     --config FILE
         * and ensure that a FILE argument actually follows.
         */
        if (strcmp(argv[argument], "--config") == 0 &&
            argument + 1 < argc)
        {
            /*
             * Save the config file path.
             */
            config_path = argv[++argument];
        }
    }

    /*
     * If --config was specified on the CLI,
     * load the config file before formally parsing the other command-line arguments.
     *
     * plasic_config_load() updates options based on the config file.
     */
    if (config_path != NULL &&
        !plasic_config_load(
            config_path, &options, &config_storage, config_error,
            sizeof(config_error)))
    {
        /*
         * Config file loading failed:
         *
         * Output:
         *   - config file path
         *   - the specific error message returned by the loader
         *
         * Then exit directly.
         */
        fprintf(stderr, "Cannot load --config %s: %s\n",
                config_path, config_error);
        return EXIT_FAILURE;
    }

    /*
     * ================================================================
     * Second CLI scan: formally parse all command-line arguments
     * ================================================================
     *
     * The config file has now been fully loaded.
     *
     * Therefore any modification to options here
     * overrides the corresponding configuration in config,
     * implementing:
     *
     *     CLI > config > defaults
     */
    for (argument = 1; argument < argc; ++argument)
    {
        /*
         * --config has already been handled in the first scan.
         *
         * The second scan only needs to skip both
         * of the following argv items:
         *
         *     --config FILE
         */
        if (strcmp(argv[argument], "--config") == 0 &&
            argument + 1 < argc)
        {
            ++argument;
        }

        else if (strcmp(argv[argument], "--restart") == 0 &&
                 argument + 1 < argc)
        {
            options.restart_input = argv[++argument];
        }

        else if (strcmp(argv[argument], "--data-dir") == 0 &&
                 argument + 1 < argc)
        {
            options.surface_data_directory = argv[++argument];
        }

        else if (strcmp(argv[argument], "--output") == 0 &&
                 argument + 1 < argc)
        {
            options.restart_output = argv[++argument];
        }

        else if (strcmp(argv[argument], "--progress") == 0 &&
                 argument + 1 < argc)
        {
            options.progress_output = argv[++argument];
        }

        else if (strcmp(argv[argument], "--frames") == 0 &&
                 argument + 1 < argc)
        {
            options.frame_output = argv[++argument];
        }

        else if (strcmp(argv[argument], "--co2-forcing") == 0 &&
                 argument + 1 < argc)
        {
            options.co2_forcing_input = argv[++argument];
        }

        else if (strcmp(argv[argument], "--ozone-forcing") == 0 &&
                 argument + 1 < argc)
        {
            options.ozone_forcing_input = argv[++argument];
        }

        else if (strcmp(argv[argument], "--co2-ppm") == 0 &&
                 argument + 1 < argc)
        {
            if (!parse_positive_float(argv[++argument], &options.co2_ppm))
            {
                fprintf(stderr, "Invalid --co2-ppm value: %s\n", argv[argument]);
                return EXIT_FAILURE;
            }
        }

        else if (strcmp(argv[argument], "--monthly-output") == 0 &&
                 argument + 1 < argc)
        {
            options.monthly_output_directory = argv[++argument];
        }

        else if (strcmp(argv[argument], "--steps") == 0 &&
                 argument + 1 < argc)
        {
            if (!parse_positive_integer(argv[++argument], &options.steps) ||
                options.steps < 1)
            {
                fprintf(stderr, "Invalid --steps value: %s\n", argv[argument]);
                return EXIT_FAILURE;
            }
        }

        else if (strcmp(argv[argument], "--frame-interval") == 0 &&
                 argument + 1 < argc)
        {
            if (!parse_positive_integer(
                    argv[++argument], &options.frame_interval) ||
                options.frame_interval < 1)
            {
                fprintf(
                    stderr, "Invalid --frame-interval value: %s\n",
                    argv[argument]);
                return EXIT_FAILURE;
            }
        }

        else if (strcmp(argv[argument], "--progress-interval") == 0 &&
                 argument + 1 < argc)
        {
            if (!parse_positive_integer(
                    argv[++argument], &options.progress_interval) ||
                options.progress_interval < 1)
            {
                fprintf(
                    stderr, "Invalid --progress-interval value: %s\n",
                    argv[argument]);
                return EXIT_FAILURE;
            }
        }

        else if (strcmp(argv[argument], "--sea-coupling-interval") == 0 &&
                 argument + 1 < argc)
        {
            if (!parse_positive_integer(
                    argv[++argument], &options.sea_coupling_interval) ||
                options.sea_coupling_interval < 1)
            {
                fprintf(
                    stderr, "Invalid --sea-coupling-interval value: %s\n",
                    argv[argument]);
                return EXIT_FAILURE;
            }
        }

        else if (strcmp(argv[argument], "--calendar") == 0 &&
                 argument + 1 < argc)
        {
            const int32_t calendar =
                calmod_calendar_from_name(argv[++argument]);

            if (calendar < 0)
            {
                fprintf(stderr,
                        "Invalid --calendar value: %s "
                        "(expected gregorian or 360day)\n",
                        argv[argument]);
                return EXIT_FAILURE;
            }
            options.calendar = calendar;
        }

        else if (strcmp(argv[argument], "--start-year") == 0 &&
                 argument + 1 < argc)
        {
            if (!parse_positive_integer(
                    argv[++argument], &options.start_year))
            {
                fprintf(stderr, "Invalid --start-year value: %s\n",
                        argv[argument]);
                return EXIT_FAILURE;
            }
        }

        else if (strcmp(argv[argument], "--start-yday") == 0 &&
                 argument + 1 < argc)
        {
            if (!parse_positive_integer(
                    argv[++argument], &options.start_yday) ||
                options.start_yday < 1)
            {
                fprintf(
                    stderr, "Invalid --start-yday value: %s\n",
                    argv[argument]);
                return EXIT_FAILURE;
            }
        }

        else if (strcmp(argv[argument], "--start-hour") == 0 &&
                 argument + 1 < argc)
        {
            if (!parse_positive_integer(
                    argv[++argument], &options.start_hour) ||
                options.start_hour > 23)
            {
                fprintf(
                    stderr, "Invalid --start-hour value: %s\n",
                    argv[argument]);
                return EXIT_FAILURE;
            }
        }

        else if (strcmp(argv[argument], "--start-minute") == 0 &&
                 argument + 1 < argc)
        {
            if (!parse_positive_integer(
                    argv[++argument], &options.start_minute) ||
                options.start_minute > 59)
            {
                fprintf(
                    stderr, "Invalid --start-minute value: %s\n",
                    argv[argument]);
                return EXIT_FAILURE;
            }
        }

        else if ((strcmp(argv[argument], "--reset-nstep") == 0 ||
                  strcmp(argv[argument], "--reset_nstep") == 0) &&
                 argument + 1 < argc)
        {
            if (!parse_positive_integer(
                    argv[++argument], &options.reset_nstep))
            {
                fprintf(stderr, "Invalid --reset-nstep value: %s\n",
                        argv[argument]);
                return EXIT_FAILURE;
            }
        }

        else if (strcmp(argv[argument], "--dynamics-only") == 0)
        {
            options.dynamics_only = 1;
        }

        else if (strcmp(argv[argument], "--help") == 0 ||
                 strcmp(argv[argument], "-h") == 0)
        {
            print_usage(stdout, argv[0]);
            return EXIT_SUCCESS;
        }

        /*
         * Unknown argument / incomplete argument
         */
        else
        {
            fprintf(stderr, "Unknown or incomplete option: %s\n",
                    argv[argument]);

            /*
             * Also print the full usage to stderr
             */
            print_usage(stderr, argv[0]);

            return EXIT_FAILURE;
        }
    }

    if (options.steps <= 0)
    {
        fprintf(stderr,
                "--steps N is required and must be a positive integer\n");
        return EXIT_FAILURE;
    }

    if (options.restart_input != NULL)
    {
        /*
         * Warm start: nstep and the calendar come from the restart file.
         * --reset-nstep is the only supported clock rebase.
         */
        if (options.start_year != 0 || options.start_yday != 1 ||
            options.start_hour != 0 || options.start_minute != 0)
        {
            fprintf(stderr,
                    "warm start reads its date from the restart file; "
                    "--start-year/--start-yday/--start-hour/--start-minute "
                    "only apply to cold starts (use --reset-nstep to rebase "
                    "the clock)\n");
            return EXIT_FAILURE;
        }
    }
    else
    {
        if (options.reset_nstep >= 0)
        {
            fprintf(stderr, "--reset-nstep requires --restart\n");
            return EXIT_FAILURE;
        }
        {
            /* An unspecified --calendar means Gregorian for a cold start. */
            const int32_t cold_calendar =
                options.calendar >= 0 ? options.calendar
                                      : PLASIC_CALENDAR_GREGORIAN;

            if (options.start_yday < 1 ||
                options.start_yday >
                    days_in_year(cold_calendar, options.start_year))
            {
                fprintf(
                    stderr,
                    "Invalid --start-yday value: %d (year %d has %d days in "
                    "the %s calendar)\n",
                    options.start_yday, options.start_year,
                    days_in_year(cold_calendar, options.start_year),
                    calmod_calendar_name(cold_calendar));
                return EXIT_FAILURE;
            }
        }
    }

    /*
     * PlaSiC runtime initialization
     */
    status = plasic_runtime_initialize_with_diagnostics(&options, stderr);

    /*
     * Runtime initialization failed.
     */
    if (status != PLASIC_RUNTIME_OK)
    {
        return EXIT_FAILURE;
    }

    if (options.restart_input == NULL && plasic_mypid == 0)
    {
        fprintf(stderr, "No --restart given; performing a cold start from "
                        "surface data in '%s'.\n",
                options.surface_data_directory);
    }

    /*
     * ================================================================
     * Run PlaSiC time integration
     * ================================================================
     */
    status = plasic_runtime_run(&options);

    /*
     * Integration failed.
     */
    if (status != PLASIC_RUNTIME_OK)
    {
        /*
         * In a parallel environment, only rank 0 prints the error.
         */
        if (plasic_mypid == 0)
        {
            fprintf(stderr, "PlaSiC integration failed: %s\n",
                    plasic_runtime_error());
        }

        /*
         * Release runtime resources.
         */
        plasic_runtime_finalize();

        return EXIT_FAILURE;
    }

    /*
     * Write the final restart
     */
    status = plasic_runtime_write_restart(options.restart_output);

    /*
     * Restart write failed.
     */
    if (status != PLASIC_RUNTIME_OK)
    {
        if (plasic_mypid == 0)
        {
            fprintf(stderr, "PlaSiC restart write failed: %s\n",
                    plasic_runtime_error());
        }

        plasic_runtime_finalize();

        return EXIT_FAILURE;
    }

    /*
     * ================================================================
     * Normal completion notice
     * ================================================================
     */
    if (plasic_mypid == 0)
    {
        printf(
            "Standalone C PlaSiC completed %d step(s); nstep=%d; "
            "restart=%s\n",
            options.steps,
            plasic_nstep,
            options.restart_output);
    }

    plasic_runtime_finalize();

    return EXIT_SUCCESS;
}
