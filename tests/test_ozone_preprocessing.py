"""Numerical and metadata tests for the historical ozone preparation workflow."""
import importlib.util
import os
import subprocess
from pathlib import Path
import tempfile
import unittest

import netCDF4 as nc
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('prepare', ROOT / 'experiments/ozone_historical/prepare_ozone.py')
prepare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)


class Preprocessing(unittest.TestCase):
    def test_pressure_overlap_and_latitude(self):
        # Source bounds at 0, 20, 60, 100 kPa, deliberately descending.
        p = np.array([80000., 40000., 10000.])
        bounds = [[100000., 60000.], [60000., 20000.], [20000., -1e-9]]
        lat = np.array([90., 0., -90.])
        lon = np.array([180., 0., 90., 270.])
        source = np.array([1., 3., 7.]) * 1e-6
        v = np.broadcast_to(source[:, None, None] + 1e-9*lat[None, :, None], (3, 3, 4))
        ps = np.array([[80000., 110000.], [90000., 95000.]])
        target_lat = np.array([30., -60.])
        interfaces = np.array([0., .25, 1.])
        actual = prepare.regrid(v, p, bounds, lat, lon, interfaces, target_lat, [359., 1.], ps)
        # Independent rectangle integration of the source, including below 100 kPa.
        column = (20000*7e-6 + 40000*3e-6 + (ps - 60000)*1e-6 + 1e-9*target_lat[:, None]*ps)
        dp = np.diff(interfaces)[:, None, None] * ps
        np.testing.assert_allclose((actual*dp).sum(axis=0), column, rtol=1e-14)
        self.assertAlmostEqual(actual[0, 0, 0], 7e-6 + 30e-9, places=15)

    def test_periodic_longitude_and_invalid_coordinates(self):
        v = np.tile(np.array([0., 1., 0., 1.])*1e-6, (2, 2, 1))
        bounds = [[0., 50000.], [50000., 100000.]]
        result = prepare.regrid(v, [1000., 100000.], bounds, [-80., 80.], [0., 90., 180., 270.],
                                [0., .5, 1.], [-89., 89.], [359., 1.], np.full((2, 2), 101325.))
        np.testing.assert_allclose(result, 1e-6/90, rtol=1e-13)
        with self.assertRaisesRegex(ValueError, 'unique'):
            prepare.regrid(v, [1000., 100000.], bounds, [-80., 80.], [0., 90., 180., 360.],
                           [0., 1.], [0.], [0.], np.array([[1e5]]))

    def test_conservative_cells_across_greenwich(self):
        # The source includes both poles and straddles the target's longitude seam.
        with tempfile.TemporaryDirectory() as tmp, nc.Dataset(Path(tmp) / 'grid.nc', 'w') as grid:
            grid.createDimension('bnds', 2)
            for axis, centres, edges in [('lat', [-90., 0., 90.], [-90., -45., 45., 90.]),
                                         ('lon', [45., 135., 225., 315.], [0., 90., 180., 270., 360.])]:
                grid.createDimension(axis, len(centres))
                coordinate = grid.createVariable(axis, 'f8', (axis,))
                coordinate[:] = centres
                coordinate.bounds = axis + '_bounds'
                grid.createVariable(coordinate.bounds, 'f8', (axis, 'bnds'))[:] = np.column_stack((edges[:-1], edges[1:]))
            _, lat, lon = prepare.sigma_grid(nlat=4)
            wy, wx = prepare.area_weights(grid, lat, lon)
            raw = np.zeros((3, 4))
            raw[-1, 0] = 1.  # One northern polar sector.
            remapped = wy @ raw @ wx.T
            band_areas = np.polynomial.legendre.leggauss(4)[1][::-1] / 2
            self.assertAlmostEqual(band_areas @ remapped.mean(axis=1), (1 - np.sin(np.pi/4)) / 8)
            self.assertAlmostEqual(remapped[0, 0], remapped[0, 1] / 2)
            np.testing.assert_array_equal(remapped[-1], 0.)

    def test_published_month_units_and_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            with nc.Dataset(Path(tmp)/'source.nc', 'w') as ds:
                ds.createDimension('time', 3)
                time = ds.createVariable('time', 'f8', ('time',))
                time.units = 'months since 1850-01-01 00:00'
                time.calendar = 'standard'
                time[:] = [0.5, 1.5, 12.5]
                np.testing.assert_array_equal(prepare.month_indices(time), [22200, 22201, 22212])
                time[0] = np.ma.masked
                with self.assertRaisesRegex(ValueError, 'missing'):
                    prepare.month_indices(time)

    def test_expected_files_required(self):
        with self.assertRaisesRegex(ValueError, 'four requested'):
            prepare.prepare([], Path('missing'), Path('missing'))


class FullConversion(unittest.TestCase):
    def test_four_file_conversion_and_fixed_seasonal_cycle(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sigma, lat, lon = prepare.sigma_grid()
            pressure_file = root/'ps.nc'
            with nc.Dataset(pressure_file, 'w') as ds:
                for name, size in [('time', 12), ('lat', 32), ('lon', 64)]:
                    ds.createDimension(name, size)
                ds.createVariable('lat', 'f8', ('lat',))[:] = lat
                ds.createVariable('lon', 'f8', ('lon',))[:] = lon
                t = ds.createVariable('time', 'f8', ('time',))
                t.units, t.calendar = 'months since 1850-01-01 00:00', 'standard'
                t[:] = np.arange(12) + .5
                ps = ds.createVariable('ps', 'f8', ('time', 'lat', 'lon'))
                ps.units = 'Pa'
                ps[:] = 100000.
            files = []
            start = 0
            for period, count in zip(prepare.PERIODS, [600, 600, 600, 180]):
                path = root/(prepare.PREFIX + period + '.nc')
                files.append(path)
                with nc.Dataset(path, 'w') as ds:
                    ds.source_id = 'UReading-CCMI-1-0'
                    for name, size in [('time', count), ('plev', 2), ('lat', 2), ('lon', 2)]:
                        ds.createDimension(name, size)
                    for name, values in [('plev', [1., 1000.]), ('lat', [-90., 90.]), ('lon', [0., 180.])]:
                        ds.createVariable(name, 'f8', (name,))[:] = values
                    ds['plev'].units = 'hPa'
                    ds['plev'].bounds = 'bounds_plev'
                    ds.createDimension('bnds', 2)
                    for axis, bounds in [('lat', [[-90., 0.], [0., 90.]]),
                                         ('lon', [[-90., 90.], [90., 270.]])]:
                        ds[axis].bounds = axis + '_bnds'
                        ds.createVariable(ds[axis].bounds, 'f8', (axis, 'bnds'))[:] = bounds
                    ds.createVariable('bounds_plev', 'f8', ('plev', 'bnds'))[:] = [[0., 500.], [500., 1050.]]
                    t = ds.createVariable('time', 'f8', ('time',))
                    t.units, t.calendar = 'months since 1850-01-01 00:00', 'standard'
                    t[:] = start + np.arange(count) + .5
                    o3 = ds.createVariable('vmro3', 'f4', ('time', 'plev', 'lat', 'lon'))
                    o3.units = 'mole mole-1'
                    o3[:] = (1e-6 + (start + np.arange(count)) * 1e-10)[:, None, None, None]
                start += count
            for fixed in (False, True):
                output = root/f'converted-{fixed}.nc'
                interrupted = output.with_suffix('.partial.nc')
                interrupted.write_bytes(b'unfinished earlier run')
                prepare.prepare(files, pressure_file, output, fixed)
                self.assertEqual(interrupted.read_bytes(), b'unfinished earlier run')
                with nc.Dataset(output) as ds:
                    self.assertEqual(ds.horizontal_remapping, 'conservative')
                    self.assertEqual(ds['o3'].chunking(), [1, 10, 16, 64])
                    self.assertEqual(ds['o3'].shape, (1980, 10, 32, 64))
                    np.testing.assert_allclose(ds['o3'][1979], 1e-6 + (11 if fixed else 1979)*1e-10, rtol=1e-6)
                    self.assertEqual(int(ds['year'][-1]), 2014)
                    self.assertEqual(int(ds['month'][-1]), 12)
                    self.assertIn('sha256', ds.source_files)
                    self.assertEqual(ds.reference_ps_sha256, prepare.sha256(pressure_file))
                    self.assertEqual(ds.converter_sha256, prepare.sha256(Path(prepare.__file__)))
                    self.assertEqual(ds['o3'].standard_name, 'mole_fraction_of_ozone_in_air')
                    self.assertEqual(ds['lev'].formula_terms, 'sigma: lev ps: ps ptop: ptop')
                    self.assertEqual(float(ds['ptop'][...]), 0.)
                    self.assertEqual(ds['ps'].dtype, np.dtype('f8'))
                    np.testing.assert_array_equal(ds['ps'][19], 100000.)
                    np.testing.assert_allclose(ds['lev_bnds'][:].mean(axis=1), ds['lev'][:])
                    t = ds['time']
                    b = ds[t.bounds][:]
                    dates = nc.num2date(t[:], t.units, t.calendar)
                    for year, length in ((1900, 28), (2000, 29)):
                        feb = (year - 1850) * 12 + 1
                        self.assertEqual(dates[feb].month, 2)
                        self.assertEqual(dates[feb].year, year)
                        self.assertEqual(b[feb, 1] - b[feb, 0], length)
                        self.assertEqual(t[feb], b[feb].mean())
                model = os.environ.get('PLASIC_TEST_MODEL')
                if model:
                    subprocess.run([str(Path(model).resolve()), '--co2-ppm', '284.316853',
                                    '--steps', '1', '--start-year', '1850',
                                    '--data-dir', str(ROOT/'src/data/surface'),
                                    '--ozone-forcing', str(output), '--output', str(root/'restart')],
                                   check=True, capture_output=True, timeout=120)
                    self.assertTrue((root/'restart').is_file())

if __name__ == '__main__':
    unittest.main()
