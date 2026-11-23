"""Integration checks for a serial T21 L10 plasic.x (make -C src test-ozone)."""

import argparse
import json
import os
import re
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from scipy.io import netcdf_file


MODEL = os.environ.get("PLASIC_TEST_MODEL")
MPI_MODEL = os.environ.get("PLASIC_TEST_MPI_MODEL")
MPIRUN = os.environ.get("MPIRUN", "mpirun")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from experiments.historical_3d_500y.read_restart import read_records

SURFACE = ROOT / "src/data/surface"
DIMS = ("time", "lev", "lat", "lon")


def write_forcing(path, values, units="kg kg-1", months=(1, 2), dtype="d",
                  marker_attribute="missing_value", years=1850):
    # Default PlaSiC stretched sigma grid, North-to-South Gaussian latitudes.
    s = np.arange(11) / 10
    interfaces = 0.75 * s + 1.75 * s**3 - 1.5 * s**4
    latitudes = np.degrees(np.arcsin(np.polynomial.legendre.leggauss(32)[0][::-1]))
    with netcdf_file(path, "w") as ds:
        for name, size in zip(DIMS, (len(months), 10, 32, 64)):
            ds.createDimension(name, size)
        ds.createVariable("lev", "f", ("lev",))[:] = (interfaces[:-1] + interfaces[1:]) / 2
        ds.createVariable("lat", "f", ("lat",))[:] = latitudes
        ds.createVariable("lon", "f", ("lon",))[:] = np.arange(64) * 360 / 64
        ds.createVariable("year", "i", ("time",))[:] = years
        ds.createVariable("month", "i", ("time",))[:] = months
        o3 = ds.createVariable("o3", dtype, DIMS)
        o3.units = units
        setattr(o3, marker_attribute, np.float64(0.123))
        o3[:] = values


class MonthlyOzone(unittest.TestCase):
    def setUp(self):
        if MODEL is None and not self._testMethodName.startswith("test_mpi_"):
            self.skipTest("set PLASIC_TEST_MODEL to a serial T21 L10 plasic.x")
        work = tempfile.TemporaryDirectory()
        self.addCleanup(work.cleanup)
        self.directory = Path(work.name)
        self.ozone = self.directory / "forcing.nc"
        self.config = self.directory / "run.json"
        # Start in the last 45 minutes of January and finish after the midnight
        # diagnostic reset, so the restart retains the February fluxes.
        self.config.write_text(json.dumps(dict(
            steps=3, data_dir=str(SURFACE), start_year=1850, start_yday=31,
            start_hour=23, start_minute=15, output=str(self.directory / "status"),
        )))
        co2 = self.directory / "co2.csv"
        co2.write_text("datetime,datenum,data,year,month,lat\n" + "".join(
            f"monthly,0,{280 + month},1850,{month},{-82.5 + 15 * band}\n"
            for month in range(1, 13) for band in range(12)
        ))
        self.co2 = ["--co2-forcing", co2]

    def run_model(self, *extra, error=None, mpi=False, wrapper=()):
        command = [*shlex.split(MPIRUN), "-np", "2", *wrapper, MPI_MODEL] if mpi else [MODEL]
        env = os.environ.copy()
        if mpi and os.geteuid() == 0:
            env.update(OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
        result = subprocess.run(
            [*command, "--config", self.config, "--step-minutes", "45", *extra],
            capture_output=True, text=True, timeout=120, env=env,
        )
        if error is None:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(error, result.stderr)
        return result

    def fluxes(self, *options):
        # Restart accumulators at the top of the atmosphere retain the sum of
        # the net shortwave and longwave fluxes even for a partial month.
        self.run_model(*self.co2, *options)
        fields = read_records(self.directory / "status")
        return np.stack([fields[name] for name in ("atsol", "atthr")])

    def test_known_300_du_column(self):
        # Obtain the actual first-step pressure through the model's independent
        # frame output, before prescribing a column that varies horizontally.
        sys.path.insert(0, str(ROOT / "plasic_app/src"))
        try:
            from plasic_app.protocol.pstream import PStreamReader
        finally:
            sys.path.pop(0)
        frames = self.directory / "initial.pstream"
        self.run_model(*self.co2, "--steps", "1", "--frames", frames, "--frame-interval", "1")
        pressure = next(frame.data.reshape(32, 64) for frame in PStreamReader(frames).read_new()
                        if frame.name == "surface_pressure").astype(float)
        mu = np.polynomial.legendre.leggauss(32)[0][::-1]
        # P2 has zero area mean. Inserting all ozone into the top layer also
        # detects treating a layer mean as a point value or omitting Δsigma.
        target_du = (300.0 + 100.0 * (3.0 * mu**2 - 1.0) / 2.0)[:, None]
        top_dsigma = 0.75 * 0.1 + 1.75 * 0.1**3 - 1.5 * 0.1**4
        q = np.zeros((2, 10, 32, 64))
        q[:, 0] = target_du * 2.14e-5 * 9.80665 / (pressure * top_dsigma)
        for units, factor in (("kg kg-1", 1.0), ("mole mole-1", 28.9644 / 47.9982)):
            with self.subTest(units=units):
                write_forcing(self.ozone, q * factor, units=units)
                result = self.run_model(*self.co2, "--ozone-forcing", self.ozone)
                columns = re.findall(r"global mean column = ([0-9.]+) DU", result.stderr)
                self.assertEqual(len(columns), 1, result.stderr)
                self.assertAlmostEqual(float(columns[0]), 300.0, delta=.01)

    def test_monthly_radiation_and_units(self):
        synthetic = self.fluxes()
        values = np.full((2, 10, 32, 64), 1e-8)
        write_forcing(self.ozone, values)
        config = json.loads(self.config.read_text())
        config["ozone_forcing"] = str(self.ozone)
        self.config.write_text(json.dumps(config))
        control = self.fluxes()

        # Change only February, in a subset of the top layer and northern hemisphere.
        values[1, 0, :16, 32:] = 1e-5
        write_forcing(self.ozone, values)
        forced = self.fluxes()
        mole_file = self.directory / "mole-fraction.nc"
        write_forcing(mole_file, values * (28.9644 / 47.9982), units="mol mol-1")
        converted = self.fluxes("--ozone-forcing", mole_file)
        for base, low, high, mol in zip(synthetic, control, forced, converted):
            self.assertGreater(np.max(np.abs(base - low)), 0.01)
            self.assertGreater(np.max(np.abs(low - high)), 0.01)
            np.testing.assert_allclose(mol, high, atol=1e-4)

    def test_interpolated_radiation(self):
        # Compare a single radiation step against the exact same model state
        # driven by an independently computed, spatially uniform ozone field.
        low, high = 1e-8, 1e-5
        for calendar, year, yday, months, years, weight in (
            ("gregorian", 1850, 31, (1, 2), 1850, 15.5 / 29.5),
            ("gregorian", 1850, 32, (1, 2), 1850, 16.5 / 29.5),
            ("gregorian", 1852, 60, (2, 3), 1852, 14.5 / 30),
            # Day 60 is March 1 in 1900 and February 29 in 2000.
            ("gregorian", 1900, 60, (2, 3), 1900, 15.0 / 29.5),
            ("gregorian", 2000, 60, (2, 3), 2000, 14.5 / 30),
            ("360_day", 1851, 1, (12, 1), (1850, 1851), 16 / 30),
        ):
            with self.subTest(calendar=calendar, year=year, yday=yday):
                self.config.write_text(json.dumps(dict(
                    steps=1, data_dir=str(SURFACE), start_year=year, start_yday=yday,
                    start_hour=12, calendar=calendar, output=str(self.directory / "status"),
                )))
                options = ("--co2-ppm", "300", "--ozone-forcing", self.ozone)
                def radiation(values):
                    write_forcing(self.ozone, values, months=months, years=years)
                    self.run_model(*options)
                    fields = read_records(self.directory / "status")
                    return np.stack([fields[k] for k in ("atsol", "atthr")])

                expected = radiation(low + weight * (high - low))
                actual = radiation(np.array([low, high])[:, None, None, None])
                self.assertGreater(np.max(np.abs(expected - radiation(low))), 0.01)
                np.testing.assert_allclose(actual, expected, rtol=1e-6, atol=1e-4)

    def test_clear_config_with_empty_cli(self):
        synthetic = self.fluxes()
        config = json.loads(self.config.read_text())
        config["ozone_forcing"] = str(self.directory / "does-not-exist.nc")
        self.config.write_text(json.dumps(config))
        np.testing.assert_array_equal(synthetic, self.fluxes("--ozone-forcing", ""))

    def test_float_storage(self):
        # A double attribute describing 0.123 differs from its float32 field
        # representation. An adjacent, valid float must still be accepted.
        options = ["--co2-ppm", "300", "--ozone-forcing", self.ozone]
        for marker in ("missing_value", "_FillValue"):
            with self.subTest(marker=marker):
                write_forcing(self.ozone, 0.123, dtype="f", marker_attribute=marker)
                self.run_model(*options, error="invalid or missing")
        values = np.full((2, 10, 32, 64), 1e-7, dtype=np.float32)
        values[0, 0, 0, 0] = np.nextafter(np.float32(0.123), np.float32(0))
        write_forcing(self.ozone, values, dtype="f")
        self.run_model(*options)

    def test_mpi_initialization_diagnostic_once(self):
        if MPI_MODEL is None:
            self.skipTest("requires a two-rank T21 L10 plasic.x")
        result = self.run_model(mpi=True, error="specify exactly one CO2 mode")
        self.assertEqual(result.stderr.count("PlaSiC initialization failed:"), 1)

    @unittest.skipUnless(shutil.which("ncdump") and shutil.which("ncgen"), "requires ncdump and ncgen")
    def test_string_units(self):
        write_forcing(self.ozone, 1e-7, dtype="f")
        expected = self.fluxes("--ozone-forcing", self.ozone)
        cdl = subprocess.check_output(["ncdump", self.ozone], text=True)
        cdl = cdl.replace('o3:units = "kg kg-1"', 'string o3:units = "kg kg-1"')
        self.assertIn("string o3:units", cdl)
        subprocess.run(["ncgen", "-4", "-o", self.ozone], input=cdl, text=True, check=True)
        np.testing.assert_array_equal(expected, self.fluxes("--ozone-forcing", self.ozone))

    def test_mpi_earth_scale_advisory(self):
        if MPI_MODEL is None:
            self.skipTest("requires a two-rank T21 L10 plasic.x")
        for mass_fraction, expected_warnings in ((0.995e-4, 0), (1.005e-4, 1)):
            with self.subTest(mass_fraction=mass_fraction):
                values = np.full((2, 10, 32, 64), 1e-7)
                # A single cell on rank one straddles the mass-fraction threshold
                # while the stored mole fraction remains below it in both runs.
                values[:, 0, -1, 0] = mass_fraction * (28.9644 / 47.9982)
                write_forcing(self.ozone, values, units="mole/mole")
                result = self.run_model(*self.co2, "--ozone-forcing", self.ozone, mpi=True)
                self.assertEqual(result.stderr.count("check units for Earth ozone"), expected_warnings)
                self.assertEqual(result.stderr.count("global mean column ="), 1)

    def test_mpi_late_month_unit_advisory(self):
        if MPI_MODEL is None:
            self.skipTest("requires a two-rank T21 L10 plasic.x")
        config = json.loads(self.config.read_text())
        config.update(start_yday=16, start_hour=0, start_minute=0, steps=2)
        self.config.write_text(json.dumps(config))
        values = np.full((3, 10, 32, 64), 1e-7)
        values[2, 0, -1, 0] = 0.2  # March is outside the two-step run.
        for february, expected in ((1e-7, 0), (0.2, 1)):
            values[1, 0, -1, 0] = february
            write_forcing(self.ozone, values, months=(1, 2, 3))
            result = self.run_model(*self.co2, "--ozone-forcing", self.ozone, mpi=True)
            self.assertEqual(result.stderr.count("check units for Earth ozone"), expected)

    def test_mpi_retry_after_failed_extension(self):
        driver = Path(MPI_MODEL).with_name("ozone_retry") if MPI_MODEL else None
        if driver is None or not driver.is_file():
            self.skipTest("requires the two-rank API driver built by make test-ozone-mpi")
        values = np.full((3, 10, 32, 64), 1e-7)
        values[2, 1, -1, 0] = np.nan  # The expanded run encounters bad March data on rank one.
        write_forcing(self.ozone, values, months=(1, 2, 3))
        env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
        result = subprocess.run([*shlex.split(MPIRUN), "-np", "2", driver, self.ozone, SURFACE],
                                capture_output=True, text=True, env=env, timeout=45)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_mpi_southern_slab_failure(self):
        if MPI_MODEL is None:
            self.skipTest("requires a two-rank T21 L10 plasic.x")
        options = ["--co2-ppm", "300", "--ozone-forcing", self.ozone]
        # Start at the January anchor. February is first needed on step two,
        # so only the startup scan of the whole run reads it before stepping.
        config = json.loads(self.config.read_text())
        config.update(start_yday=16, start_hour=0, start_minute=0)
        self.config.write_text(json.dumps(config))
        for record, month in enumerate(("01", "02")):
            with self.subTest(month=month):
                values = np.full((2, 10, 32, 64), 1e-7)
                values[record, 9, 31, 63] = -1.0
                write_forcing(self.ozone, values)
                error = f"invalid or missing ozone forcing value at 1850-{month}"
                result = self.run_model(*options, mpi=True, error=error)
                self.assertEqual(result.stderr.count(error), 1, result.stderr)
                # Rank 1 reports the global latitude index of its slab.
                self.assertIn(f"o3(time={record}, lev=9, lat=31, lon=63) = -1", result.stderr)
                self.assertIn("PlaSiC initialization failed", result.stderr)

    def rank_one_options(self, option, value):
        launcher = self.directory / "rank_options.py"
        launcher.write_text(
            "import os, sys\n"
            "rank = os.environ.get('OMPI_COMM_WORLD_RANK', os.environ.get('PMI_RANK', os.environ.get('PMIX_RANK')))\n"
            "args = sys.argv[3:]\n"
            "if rank == '1': args += sys.argv[1:3]\n"
            "os.execv(args[0], args)\n"
        )
        return (sys.executable, launcher, option, value)

    def test_mpi_local_ozone_time_axis(self):
        if MPI_MODEL is None:
            self.skipTest("requires a two-rank T21 L10 plasic.x")
        write_forcing(self.ozone, 1e-7, months=(1, 2))
        local = self.directory / "local.nc"
        options = ("--co2-ppm", "300", "--ozone-forcing", self.ozone)
        for months in ((1,), (2, 3), (1, 2), ()):
            with self.subTest(months=months):
                if months:
                    write_forcing(local, 1e-7, months=months)
                wrapper = self.rank_one_options("--ozone-forcing", local if months else "")
                error = None if months == (1, 2) else "ozone forcing time axis differs on rank 1"
                result = self.run_model(*options, mpi=True, wrapper=wrapper, error=error)
                if error:
                    self.assertEqual(result.stderr.count("PlaSiC initialization failed"), 1)

    def test_warm_start_clock_reset(self):
        # Save an initial state at the January/February boundary. Use the
        # same restart state both times, changing only its calendar label.
        self.run_model(*self.co2, "--steps", "1")
        initial = self.directory / "initial.restart"
        (self.directory / "status").rename(initial)
        nstep = read_records(initial)["nstep"].view("<i4").item()
        self.config.write_text(json.dumps(dict(
            steps=2, restart=str(initial), output=str(self.directory / "status"),
        )))
        options = ["--ozone-forcing", self.ozone]
        write_forcing(self.ozone, 1e-7, months=(2,))
        self.run_model(*self.co2, *options)
        # Shift back by one day while keeping the February state.
        reset = ["--reset-nstep", str(nstep - 32)]
        self.run_model(*self.co2, *options, *reset, error="1850-01")
        write_forcing(self.ozone, 1e-7, months=(1,))
        self.run_model(*self.co2, *options, *reset)

    def test_final_restart_clock_limit(self):
        self.run_model("--co2-ppm", "300", "--calendar", "360day", "--steps", "1")
        output = self.directory / "status"
        initial = output.rename(self.directory / "initial.restart")
        self.config.write_text(json.dumps(dict(steps=1, restart=str(initial), output=str(output))))
        limit = 2**31 - 1
        year, yday = divmod(limit // 32, 360)
        write_forcing(self.ozone, 1e-7, months=(1 + yday // 30,), years=year)
        for forcing in ((), ("--ozone-forcing", self.ozone)):
            with self.subTest(forcing=forcing):
                options = ("--co2-ppm", "300", "--dynamics-only", *forcing, "--reset-nstep")
                self.run_model(*options, str(limit - 1))
                records = read_records(output)
                self.assertEqual(records["nstep"].view("<i4").item(), limit)
                output.unlink()
                self.run_model(*options, str(limit), error="model clock range")
                self.assertFalse(output.exists())

    @unittest.skipUnless(shutil.which("ncdump"), "requires ncdump")
    def test_monthly_provenance(self):
        # The monthly writer emits only complete model months.
        self.config.write_text(json.dumps(dict(
            steps=28 * 32, data_dir=str(SURFACE), start_year=1850, start_yday=32,
            output=str(self.directory / "status"),
        )))
        write_forcing(self.ozone, 1e-7, months=(2,))
        for label, args, source in (("synthetic", [], "synthetic"),
                                    ("external", ["--ozone-forcing", self.ozone], str(self.ozone))):
            with self.subTest(label=label):
                output = self.directory / label
                output.mkdir()
                self.run_model(*self.co2, *args, "--monthly-output", output)
                header = subprocess.check_output(["ncdump", "-h", output / "ta.nc"], text=True)
                self.assertIn(f':ozone_forcing_source = "{source}"', header)
                self.assertIn(f':co2_forcing_source = "{self.co2[1]}"', header)
                from netCDF4 import Dataset
                with Dataset(output / "toz.nc") as column, Dataset(output / "ps.nc") as pressure:
                    self.assertEqual(column["toz"].units, "1e-5 m")
                    self.assertEqual(column["toz"].cell_methods, "time: mean")
                    self.assertEqual(column.ozone_forcing_source, source)
                    np.testing.assert_array_equal(column["time_bnds"][:], pressure["time_bnds"][:])
                    actual = column["toz"][:]
                    if label == "external":
                        # Constant mass fraction makes the column proportional to mean pressure.
                        expected = 1e-7 * pressure["ps"][:].astype(float) / (9.80665 * 2.14e-5)
                        np.testing.assert_allclose(actual, expected, rtol=1e-5)
                    else:
                        mu = np.sin(np.deg2rad(column["lat"][:].astype(float)))
                        days = 32 + np.arange(28 * 32) / 32
                        season = np.cos(2 * np.pi * (days - 90) / 365).mean()
                        expected = 1000 * (.25 + .11 * abs(mu) + .08 * mu * season)
                        np.testing.assert_allclose(actual, np.broadcast_to(expected[None, :, None], actual.shape),
                                                   rtol=0, atol=.002)

    def test_startup_validation(self):
        # The three-step run reaches February on its last step. A gap or a bad
        # February value must stop the run before the first step.
        options = ["--co2-ppm", "300", "--ozone-forcing", self.ozone]
        write_forcing(self.ozone, 1e-7, months=(1,))
        result = self.run_model(*options, error="but the run needs 1850-01 to 1850-02")
        self.assertIn("PlaSiC initialization failed", result.stderr)
        values = np.full((2, 10, 32, 64), 1e-7)
        values[1, 3, 20, 7] = np.nan
        write_forcing(self.ozone, values)
        result = self.run_model(*options, error="o3(time=1, lev=3, lat=20, lon=7) = nan")
        self.assertIn("PlaSiC initialization failed", result.stderr)

    def test_small_negative_values(self):
        # Regridding noise down to -1e-9 kg/kg reads as zero ozone.
        values = np.full((2, 10, 32, 64), 1e-7)
        values[:, 0, 5, 9] = 0.0
        write_forcing(self.ozone, values)
        zero = self.fluxes("--ozone-forcing", self.ozone)
        for mass, warnings in ((-1e-12, 0), (-5e-10, 1)):
            values[:, 0, 5, 9] = mass
            write_forcing(self.ozone, values)
            np.testing.assert_array_equal(zero, self.fluxes("--ozone-forcing", self.ozone))
            result = self.run_model(*self.co2, "--ozone-forcing", self.ozone)
            self.assertEqual(result.stderr.count("check the input regridding"), warnings)

    def test_unit_spellings(self):
        values = np.full((2, 10, 32, 64), 1e-7)
        for reference, spellings in (("kg kg-1", ("kg kg**-1", "kg kg^-1", "kg/kg", " kg  kg-1 ", "kg*kg^-1")),
                                     ("mol mol-1", ("mol mol**-1", "mol/mol", "mole mole-1", "mole/mole", "mol mol-1".ljust(256)))):
            write_forcing(self.ozone, values, units=reference)
            expected = self.fluxes("--ozone-forcing", self.ozone)
            for units in spellings:
                with self.subTest(units=units):
                    write_forcing(self.ozone, values, units=units)
                    np.testing.assert_array_equal(expected, self.fluxes("--ozone-forcing", self.ozone))

    def test_invalid_inputs(self):
        options = ["--co2-ppm", "300", "--ozone-forcing", self.ozone]
        self.run_model(*options, error="cannot open ozone forcing")
        write_forcing(self.ozone, 1e-6, months=(1,))
        self.run_model(*options, error="1850-02")
        write_forcing(self.ozone, 1e-6, months=(1, 3))
        self.run_model(*options, error="consecutive")
        for units in ("ppmv", "kg / kg**", "mol/mol^", "kg^kg-1", "kg kg***-1", "mol mol-1".ljust(80) + "x"):
            with self.subTest(units=units):
                write_forcing(self.ozone, 1e-6, units=units)
                self.run_model(*options, error="units")
        for invalid in (float("nan"), -2e-9, 0.123):
            with self.subTest(invalid=invalid):
                write_forcing(self.ozone, invalid)
                self.run_model(*options, error="invalid or missing")
        write_forcing(self.ozone, 1e-6)
        with netcdf_file(self.ozone, "a") as ds:
            lat = ds.variables["lat"]
            lat[:] = lat[::-1].copy()
        self.run_model(*options, error="coordinate 'lat'")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", nargs="?", type=Path)
    parser.add_argument("--mpi-model", type=Path)
    parser.add_argument("--mpirun", default=MPIRUN)
    args, unittest_args = parser.parse_known_args()
    MODEL = args.model.resolve() if args.model else MODEL
    MPI_MODEL = args.mpi_model.resolve() if args.mpi_model else MPI_MODEL
    MPIRUN = args.mpirun
    unittest.main(argv=[sys.argv[0], *unittest_args])
