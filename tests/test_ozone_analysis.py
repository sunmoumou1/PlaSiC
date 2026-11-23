"""Test area/time weighting and the scientific contrast of experiment commands."""
import importlib.util
import json
from datetime import datetime
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import netCDF4 as nc
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


analysis = load('analysis', 'experiments/ozone_historical/fig1/compute_temperature.py')
runner = load('runner', 'experiments/ozone_historical/run_experiments.py')


verify = load('verify_restart', 'experiments/ozone_historical/verify_restart.py')


class HistoricalAnalysis(unittest.TestCase):
    def test_continuous_restart_checks_monthly_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            continuous, first, second = (Path(tmp) / name for name in ('full', 'first', 'second'))
            for directory, start, count in ((continuous, 0, 24), (first, 0, 12), (second, 12, 12)):
                (directory / 'monthly').mkdir(parents=True)
                (directory / 'final.restart').write_bytes(b'same final state')
                with nc.Dataset(directory / 'monthly/tas.nc', 'w') as ds:
                    ds.createDimension('time', count)
                    temp = ds.createVariable('tas', 'f4', ('time',))
                    temp.units = 'K'
                    temp.valid_range = np.array([150., 330.], dtype='f4')
                    temp[:] = np.arange(start, start + count) + 280.
            verify.compare(continuous, first, second)
            with nc.Dataset(second / 'monthly/tas.nc', 'a') as ds:
                ds['tas'].units = 'degree_Celsius'
            with self.assertRaisesRegex(AssertionError, 'attribute units'):
                verify.compare(continuous, first, second)

    def test_gaussian_and_calendar_weighting(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'tas.nc'
            nodes, _ = np.polynomial.legendre.leggauss(32)
            bounds_dates = [datetime(1850, m, 1) for m in range(1, 13)] + [datetime(1851, 1, 1)]
            days = np.array([(b-a).days for a, b in zip(bounds_dates[:-1], bounds_dates[1:])])
            with nc.Dataset(path, 'w') as ds:
                for name, size in [('time', 12), ('lat', 32), ('lon', 64), ('bnds', 2)]:
                    ds.createDimension(name, size)
                ds.createVariable('lat', 'f8', ('lat',))[:] = np.rad2deg(np.arcsin(nodes[::-1]))
                t = ds.createVariable('time', 'f8', ('time',))
                t.units, t.calendar, t.bounds = 'days since 1850-01-01', 'proleptic_gregorian', 'time_bnds'
                edges = nc.date2num(bounds_dates, t.units, t.calendar)
                t[:] = (edges[:-1]+edges[1:])/2
                ds.createVariable('time_bnds', 'f8', ('time', 'bnds'))[:] = np.stack([edges[:-1],edges[1:]],axis=1)
                v = ds.createVariable('tas', 'f8', ('time', 'lat', 'lon'))
                v.units = 'K'
                v[:] = 280 + np.arange(12)[:,None,None] + nodes[::-1][None,:,None]**2
            years, values = analysis.annual_temperature(path, 1850, 1850)
            np.testing.assert_array_equal(years, [1850])
            np.testing.assert_allclose(values, 280 + np.average(np.arange(12),weights=days) + 1/3, rtol=1e-14)
            with nc.Dataset(path, 'a') as ds:
                ds['time_bnds'][-1, 1] -= 1
            with self.assertRaisesRegex(ValueError, 'complete calendar month'):
                analysis.annual_temperature(path, 1850, 1850)

    def test_forcing_contrasts_and_total_steps(self):
        for case in runner.CASES:
            cmd = runner.command(Path('/model'), Path('/ic'), Path('/o3'), Path('/output'), case)
            self.assertEqual(cmd[cmd.index('--steps')+1], '1928480')
            self.assertEqual('--co2-ppm' in cmd, case == 'ozone_only')
            self.assertEqual('--co2-forcing' in cmd, case != 'ozone_only')
            self.assertIn('--ozone-forcing', cmd)
            self.assertEqual(cmd[cmd.index('--step-minutes') + 1], '45')
        cmd = runner.command(Path('/model'), Path('/ic'), Path('/o3'), Path('/out'),
                             'co2_only', mpirun='srun -n 2')
        self.assertEqual(cmd[:4], ['srun', '-n', '2', '/model'])

    def test_historical_clock_with_archived_and_current_executables(self):
        with patch.object(runner.subprocess, 'check_output') as help_text:
            for usage, expected in (('--ozone-forcing FILE', ()),
                                    ('--ozone-forcing FILE --step-minutes M', ('--step-minutes', '45'))):
                with self.subTest(usage=usage):
                    help_text.return_value = usage
                    self.assertEqual(runner.historical_clock('/model'), expected)
                    cmd = runner.command('/model', '/ic', '/o3', Path('/out'), 'ozone_only', clock=expected)
                    self.assertEqual('--step-minutes' in cmd, bool(expected))
                    self.assertEqual(cmd[cmd.index('--steps') + 1], '1928480')
            help_text.return_value = '--steps N'
            with self.assertRaisesRegex(ValueError, 'external ozone'):
                runner.historical_clock('/model')

    def test_gregorian_centuries(self):
        self.assertEqual(runner.year_start(1850), 21622368)
        for year, feb_days in ((1900, 28), (2000, 29)):
            with self.subTest(year=year):
                self.assertEqual(runner.year_start(year + 1) - runner.year_start(year),
                                 (337 + feb_days) * 32)
                # Only February is warm, so using a 28-day February in 2000
                # changes the annual mean by a measurable, known amount.
                with tempfile.TemporaryDirectory() as tmp:
                    path = Path(tmp) / 'tas.nc'
                    dates = [datetime(year, m, 1) for m in range(1, 13)] + [datetime(year + 1, 1, 1)]
                    with nc.Dataset(path, 'w') as ds:
                        for name, size in [('time', 12), ('lat', 1), ('lon', 1), ('bnds', 2)]:
                            ds.createDimension(name, size)
                        ds.createVariable('lat', 'f8', ('lat',))[:] = 0.
                        time = ds.createVariable('time', 'f8', ('time',))
                        time.units, time.calendar, time.bounds = 'days since 1850-01-01', 'proleptic_gregorian', 'time_bnds'
                        edges = nc.date2num(dates, time.units, time.calendar)
                        time[:] = (edges[:-1] + edges[1:]) / 2
                        ds.createVariable('time_bnds', 'f8', ('time', 'bnds'))[:] = np.column_stack((edges[:-1], edges[1:]))
                        temp = ds.createVariable('tas', 'f8', ('time', 'lat', 'lon'))
                        temp.units = 'K'
                        temp[:] = 280.
                        temp[1] = 290.
                    result = analysis.annual_temperature(path, year, year)[1][0]
                    self.assertAlmostEqual(result, 280 + 10 * feb_days / (337 + feb_days), places=11)

    def test_resume_requires_both_checkpoint_checksums(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case = 'co2_ozone'
            directory = root / case / '1850'
            (directory / 'monthly').mkdir(parents=True)
            initial = root / 'initial.restart'
            initial.write_bytes(b'initial state')
            args = SimpleNamespace(output=root, restart=initial, historical_ozone=root / 'ozone.nc')
            signature = dict(inputs={}, case=case, year=1850, restart=runner.digest(initial))
            hashes = {}
            for filename in ('final.restart', 'monthly/tas.nc'):
                path = directory / filename
                path.write_bytes(filename.encode())
                hashes[filename] = runner.digest(path)
            manifest = directory / 'run.json'
            for missing in hashes:
                incomplete = {key: value for key, value in hashes.items() if key != missing}
                manifest.write_text(json.dumps(dict(signature=signature, output_sha256=incomplete)))
                with self.assertRaisesRegex(ValueError, 'expected checksums'):
                    runner.segment(args, case, 1850, {})
            manifest.write_text(json.dumps(dict(signature=signature, output_sha256=hashes)))
            runner.segment(args, case, 1850, {})
            for name in ('ta', 'co2', 'thetao', 'rsdt', 'rsut', 'rlut'):
                rel = f'monthly/{name}.nc'
                (directory / rel).write_text(name)
                hashes[rel] = runner.digest(directory / rel)
            manifest.write_text(json.dumps(dict(signature=signature, checkpoint_version=2, output_sha256=hashes)))
            runner.segment(args, case, 1850, {})
            (directory / 'monthly/thetao.nc').write_text('damaged')
            with self.assertRaisesRegex(ValueError, 'thetao.nc: checksum changed'):
                runner.segment(args, case, 1850, {})


if __name__ == '__main__':
    unittest.main()
