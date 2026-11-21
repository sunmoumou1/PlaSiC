"""Integration checks for a serial T21 L10 plasic.x (make -C src test-ozone)."""

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from scipy.io import netcdf_file


MODEL = Path(sys.argv.pop(1)).resolve()
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from experiments.forecast_era5.scripts.restart_io import decode_float, read_raw_records

SURFACE = ROOT / "src/data/surface"
DIMS = ("time", "lev", "lat", "lon")


def write_forcing(path, values, units="kg kg-1", months=(1, 2)):
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
        ds.createVariable("year", "i", ("time",))[:] = 1850
        ds.createVariable("month", "i", ("time",))[:] = months
        o3 = ds.createVariable("o3", "d", DIMS)
        o3.units = units
        o3.missing_value = np.float64(0.123)
        o3[:] = values


class MonthlyOzone(unittest.TestCase):
    def setUp(self):
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

    def run_model(self, *extra, error=None):
        result = subprocess.run(
            [MODEL, "--config", self.config, *extra], capture_output=True, text=True,
        )
        if error is None:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(error, result.stderr)

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
        jan_control = self.fluxes("--steps", "1")
        control = self.fluxes()

        # Change only February, in a subset of the top layer and northern hemisphere.
        values[1, 0, :16, 32:] = 1e-5
        write_forcing(self.ozone, values)
        np.testing.assert_allclose(jan_control, self.fluxes("--steps", "1"), atol=1e-4)
        forced = self.fluxes()
        mole_file = self.directory / "mole-fraction.nc"
        write_forcing(mole_file, values * (28.9644 / 47.9982), units="mol mol-1")
        converted = self.fluxes("--ozone-forcing", mole_file)
        for base, low, high, mol in zip(synthetic, control, forced, converted):
            self.assertGreater(np.max(np.abs(base - low)), 0.01)
            self.assertGreater(np.max(np.abs(low - high)), 0.01)
            np.testing.assert_allclose(mol, high, atol=1e-4)

    def test_invalid_inputs(self):
        options = ["--co2-ppm", "300", "--ozone-forcing", self.ozone]
        self.run_model(*options, error="cannot open ozone forcing")
        write_forcing(self.ozone, 1e-6, months=(1,))
        self.run_model(*options, error="1850-02")
        write_forcing(self.ozone, 1e-6, months=(1, 3))
        self.run_model(*options, error="consecutive")
        write_forcing(self.ozone, 1e-6, units="ppmv")
        self.run_model(*options, error="units")
        for invalid in (float("nan"), -1e-50, 0.123):
            with self.subTest(invalid=invalid):
                write_forcing(self.ozone, invalid)
                self.run_model(*options, error="invalid or missing")
        write_forcing(self.ozone, 1e-6)
        with netcdf_file(self.ozone, "a") as ds:
            lat = ds.variables["lat"]
            lat[:] = lat[::-1].copy()
        self.run_model(*options, error="coordinate 'lat'")


if __name__ == "__main__":
    unittest.main()
