"""Integration checks for a serial T21 L10 plasic.x (make -C src test-ozone)."""

import argparse
import json
import os
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
from experiments.forecast_era5.scripts.restart_io import decode_float, read_raw_records

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
            [*command, "--config", self.config, *extra],
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
        fields = dict(read_raw_records(self.directory / "status"))
        return np.stack([decode_float(fields[name], 32 * 64) for name in ("atsol", "atthr")])

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
                    fields = dict(read_raw_records(self.directory / "status"))
                    return np.stack([decode_float(fields[k], 32 * 64) for k in ("atsol", "atthr")])

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

    @unittest.skipUnless(shutil.which("ncdump") and shutil.which("ncgen"), "requires ncdump and ncgen")
    def test_string_units(self):
        write_forcing(self.ozone, 1e-7, dtype="f")
        expected = self.fluxes("--ozone-forcing", self.ozone)
        cdl = subprocess.check_output(["ncdump", self.ozone], text=True)
        cdl = cdl.replace('o3:units = "kg kg-1"', 'string o3:units = "kg kg-1"')
        self.assertIn("string o3:units", cdl)
        subprocess.run(["ncgen", "-4", "-o", self.ozone], input=cdl, text=True, check=True)
        np.testing.assert_array_equal(expected, self.fluxes("--ozone-forcing", self.ozone))

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

    def test_mpi_nonroot_surface_error(self):
        if MPI_MODEL is None:
            self.skipTest("requires a two-rank T21 L10 plasic.x")
        # Only rank 1 uses the incomplete data directory.
        surface = self.directory / "surface"
        (surface / "T21").mkdir(parents=True)
        orography = "N032_surf_orography.nc"
        launcher = self.directory / "local_surface.py"
        launcher.write_text(
            "import os, sys\n"
            "rank = os.environ.get('OMPI_COMM_WORLD_RANK', os.environ.get('PMI_RANK', os.environ.get('PMIX_RANK')))\n"
            "args = sys.argv[2:]\n"
            "if rank == '1': args += ['--data-dir', sys.argv[1]]\n"
            "os.execv(args[0], args)\n"
        )
        for field in ("orography", "land_sea_mask"):
            with self.subTest(field=field):
                if field == "land_sea_mask":
                    (surface / "T21" / orography).symlink_to(SURFACE / "T21" / orography)
                result = self.run_model("--co2-ppm", "300", mpi=True,
                                        wrapper=(sys.executable, launcher, surface),
                                        error=f"surface field '{field}'")
                self.assertEqual(result.stderr.count("PlaSiC initialization failed"), 1, result.stderr)
                self.assertIn(str(surface), result.stderr)

    def test_warm_start_clock_reset(self):
        # Save an initial state at the January/February boundary. Use the
        # same restart state both times, changing only its calendar label.
        self.run_model(*self.co2, "--steps", "1")
        initial = self.directory / "initial.restart"
        (self.directory / "status").rename(initial)
        nstep = int(np.frombuffer(dict(read_raw_records(initial))["nstep"], dtype="<i4")[0])
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
        values[:, 0, 5, 9] = -1e-12
        write_forcing(self.ozone, values)
        np.testing.assert_array_equal(zero, self.fluxes("--ozone-forcing", self.ozone))

    def test_unit_spellings(self):
        values = np.full((2, 10, 32, 64), 1e-7)
        for reference, spellings in (("kg kg-1", ("kg kg**-1", "kg kg^-1", "kg/kg", " kg  kg-1 ")),
                                     ("mol mol-1", ("mol mol**-1", "mol/mol"))):
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
        write_forcing(self.ozone, 1e-6, units="ppmv")
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
