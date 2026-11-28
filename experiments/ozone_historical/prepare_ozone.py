"""Convert UReading-CCMI-1-0 monthly ozone to PlaSiC's native sigma grid.

The prescribed sigma field uses an explicitly supplied, fixed monthly surface
pressure climatology. It does not follow instantaneous surface-pressure changes.
"""
from __future__ import annotations

import argparse
from datetime import date
import json
import os
from pathlib import Path
import re
import sys
import subprocess
from tempfile import TemporaryDirectory

import netCDF4 as nc
import numpy as np
from scipy.interpolate import RegularGridInterpolator

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from experiments.ozone_historical.common import sha256


PERIODS = ('185001-189912', '190001-194912', '195001-199912', '200001-201412')
PREFIX = 'vmro3_input4MIPs_ozone_CMIP_UReading-CCMI-1-0_gn_'


def finite(variable, key=slice(None)):
    a = np.ma.asarray(variable[key])
    if np.any(np.ma.getmaskarray(a)) or not np.isfinite(a).all():
        raise ValueError(f'{variable.name}: missing or nonfinite values')
    return np.asarray(a, dtype=float)


def month_indices(time):
    """The published CMIP6 file uses CF nonstandard 'months since' units."""
    values = finite(time)
    match = re.fullmatch(r'months since (\d{4})-(\d{2})-01(?: 00:00(?::00)?)?', time.units.strip())
    if match:
        if getattr(time, 'calendar', 'standard') not in ('standard', 'gregorian'):
            raise ValueError('months-since input requires the published standard calendar')
        year, month = map(int, match.groups())
        return year * 12 + month - 1 + np.floor(values).astype('i8')
    dates = nc.num2date(values, time.units, getattr(time, 'calendar', 'standard'))
    return np.array([d.year * 12 + d.month - 1 for d in dates])


def sigma_edges(nlev=10):
    s = np.arange(nlev + 1) / nlev
    return .75 * s + 1.75 * s**3 - 1.5 * s**4


def sigma_grid(nlev=10, nlat=32):
    half = sigma_edges(nlev)
    return ((half[:-1] + half[1:]) / 2,
            np.degrees(np.arcsin(np.polynomial.legendre.leggauss(nlat)[0][::-1])),
            np.arange(nlat * 2) * 180 / nlat)


def area_weights(dataset, target_lat, target_lon):
    """Fraction of each Gaussian cell covered by each published source cell."""
    source_lat, source_lon = finite(dataset['lat']), finite(dataset['lon'])
    source_y = np.sin(np.deg2rad(finite(dataset[dataset['lat'].bounds])))
    source_x = finite(dataset[dataset['lon'].bounds])
    # Align each periodic longitude interval with its centre in [0, 360).
    source_x = (source_x - source_lon[:, None] + 180) % 360 - 180 + (source_lon % 360)[:, None]
    source_y = np.sort(source_y[np.argsort(source_lat)], axis=1)
    source_x = np.sort(source_x[np.argsort(source_lon % 360)], axis=1)
    mu, bands = np.polynomial.legendre.leggauss(len(target_lat))
    np.testing.assert_allclose(np.sin(np.deg2rad(target_lat)), mu[::-1], atol=1e-6, rtol=0)
    edges = np.r_[-1., -1. + np.cumsum(bands)]
    edges[-1] = 1.
    target_y = np.column_stack((edges[:-1], edges[1:]))[::-1]
    half_width = 180. / len(target_lon)
    target_x = np.asarray(target_lon)[:, None] + np.array([-half_width, half_width])

    def overlap(target, source):
        width = np.minimum(target[:, 1, None], source[:, 1]) - np.maximum(target[:, 0, None], source[:, 0])
        return np.maximum(width, 0.) / (target[:, 1] - target[:, 0])[:, None]

    wy = overlap(target_y, source_y)
    wx = sum(overlap(target_x, source_x + shift) for shift in (-360., 0., 360.))
    for weights in (wy, wx):
        np.testing.assert_allclose(weights.sum(axis=1), 1., rtol=0, atol=1e-10)
    return wy, wx


def regrid(values, pressure, bounds, lat, lon, interfaces, target_lat, target_lon, ps, weights=None):
    """Integrate horizontally remapped mixing ratios over target pressure layers.

    Production uses published cell overlap weights. Point interpolation is
    retained for comparison with previously prepared bilinear inputs.
    """
    if np.any(values < 0) or np.any(values > 1):
        raise ValueError('ozone mole fraction must lie in [0, 1]')
    pressure, lat, lon = map(np.asarray, (pressure, lat, lon))
    orders = [np.argsort(pressure), np.argsort(lat), np.argsort(lon % 360)]
    p, y, x = pressure[orders[0]], lat[orders[1]], (lon % 360)[orders[2]]
    for axis in (p, y, x):
        if not np.isfinite(axis).all() or np.any(np.diff(axis) <= 0):
            raise ValueError('source coordinates must be finite and unique')
    if np.any(ps <= 0) or not np.isfinite(ps).all():
        raise ValueError('reference surface pressure must be finite and positive')
    cells = np.sort(np.asarray(bounds)[orders[0]], axis=1)
    if (cells.shape != (len(p), 2) or not np.isfinite(cells).all()
            or np.any(p < cells[:, 0]) or np.any(p > cells[:, 1])
            or abs(cells[0, 0]) > 1e-5 or not np.allclose(cells[:-1, 1], cells[1:, 0],
                                                          rtol=1e-8, atol=1e-8)):
        raise ValueError('source pressure cells must join from zero and contain their pressure levels')
    edges = np.r_[0., cells[:, 1]]
    if np.any(np.diff(edges) <= 0):
        raise ValueError('source pressure cells must have positive thickness')
    # Put the pressure and horizontal axes in increasing order for interpolation.
    sorted_values = values[np.ix_(*orders)]
    if weights is not None:
        wy, wx = weights
        horizontal = (wy @ sorted_values) @ wx.T
    else:
        x = np.r_[x[-1] - 360, x, x[0] + 360]
        wrapped = np.concatenate((sorted_values[..., -1:], sorted_values, sorted_values[..., :1]), axis=2)
        query = np.stack(np.meshgrid(np.clip(target_lat, y[0], y[-1]), np.asarray(target_lon) % 360,
                                     indexing='ij'), axis=-1)
        horizontal = RegularGridInterpolator((y, x), wrapped.transpose(1, 2, 0))(query).transpose(2, 0, 1)
    totals = np.cumsum(horizontal * np.diff(edges)[:, None, None], axis=0)
    cumulative = np.concatenate((np.zeros_like(totals[:1]), totals), axis=0)

    # Extend the lowest mixing ratio down to the local surface where ps exceeds
    # the lowest published boundary. The model top is the 0 Pa interface.
    target = np.asarray(interfaces)[:, None, None] * ps
    index = np.clip(np.searchsorted(edges, target, side='right') - 1, 0, len(p) - 1)
    iy, ix = np.indices(ps.shape)
    integral = cumulative[index, iy, ix] + (target - edges[index]) * horizontal[index, iy, ix]
    return np.diff(integral, axis=0) / np.diff(target, axis=0)


def write_cf_grid(out, months, sigma, interfaces, lat, lon, reference_ps):
    """Describe the source calendar months and the pressure grid used by the converter."""
    out.Conventions = 'CF-1.10'
    out.createDimension('bnds', 2)
    for name, values, standard_name, units, axis in (
        ('lev', sigma, 'atmosphere_sigma_coordinate', '1', 'Z'),
        ('lat', lat, 'latitude', 'degrees_north', 'Y'),
        ('lon', lon, 'longitude', 'degrees_east', 'X'),
    ):
        coordinate = out.createVariable(name, 'f8', (name,))
        coordinate[:] = values
        coordinate.setncatts(dict(standard_name=standard_name, units=units, axis=axis))
    out['lev'].setncatts(dict(positive='down', bounds='lev_bnds', computed_standard_name='air_pressure',
                             formula_terms='sigma: lev ps: ps ptop: ptop'))
    bounds = out.createVariable('lev_bnds', 'f8', ('lev', 'bnds'))
    bounds[:] = np.column_stack((interfaces[:-1], interfaces[1:]))
    bounds.formula_terms = 'sigma: lev_bnds ps: ps ptop: ptop'
    ptop = out.createVariable('ptop', 'f8')
    ptop.setncatts(dict(units='Pa', standard_name='air_pressure_at_top_of_atmosphere_model',
                        long_name='Pressure at the model top'))
    ptop.assignValue(0.)
    pressure = out.createVariable('ps', 'f8', ('time', 'lat', 'lon'), zlib=True, complevel=1)
    pressure.setncatts(dict(standard_name='surface_air_pressure', units='Pa',
                            long_name='Monthly reference surface pressure used for vertical remapping'))
    for start in range(0, len(months), 12):
        pressure[start:start + 12] = reference_ps

    time = out.createVariable('time', 'f8', ('time',))
    time.setncatts(dict(standard_name='time', axis='T', calendar='proleptic_gregorian',
                       units='days since 1850-01-01', bounds='time_bnds',
                       comment='Centres of source calendar-month means. PlaSiC interpolates using year and month.'))
    base = date(1850, 1, 1)
    edges = [(date(i // 12, i % 12 + 1, 1) - base).days
             for i in range(int(months[0]), int(months[-1]) + 2)]
    out.createVariable('time_bnds', 'f8', ('time', 'bnds'))[:] = np.column_stack((edges[:-1], edges[1:]))
    time[:] = (np.asarray(edges[:-1]) + edges[1:]) / 2
    for name, values in [('year', months // 12), ('month', months % 12 + 1)]:
        coord = out.createVariable(name, 'i4', ('time',))
        coord.long_name = 'Calendar ' + name + ' read by PlaSiC'
        coord[:] = values
    out.converter_sha256 = sha256(Path(__file__))
    try:
        script = Path(__file__).resolve()
        subprocess.run(['git', 'diff', '--quiet', 'HEAD', '--', script.name], cwd=script.parent,
                       check=True, capture_output=True)
        revision = subprocess.run(['git', 'log', '-1', '--format=%H', '--', script.name],
                                  cwd=script.parent, capture_output=True, text=True, check=True)
    except (OSError, subprocess.CalledProcessError):
        return
    if revision.stdout.strip():
        out.converter_git_revision = revision.stdout.strip()


def prepare(files, pressure_file, output, fixed=False, checksums=None, nlev=10):
    expected = [PREFIX + p + '.nc' for p in PERIODS]
    files = sorted(map(Path, files))
    if [p.name for p in files] != expected:
        raise ValueError('supply exactly the four requested UReading-CCMI-1-0 historical files')
    if output.exists():
        raise FileExistsError(output)
    if nlev < 1:
        raise ValueError("the number of sigma levels must be positive")
    with nc.Dataset(pressure_file) as ds:
        ps = finite(ds['ps'])
        if ds['ps'].units != 'Pa' or ds['ps'].dimensions != ('time', 'lat', 'lon'):
            raise ValueError('reference ps must be ps(time,lat,lon) in Pa')
        nlat = len(ds.dimensions['lat'])
        if ps.shape != (12, nlat, 2 * nlat):
            raise ValueError('reference ps must contain twelve monthly fields on a PlaSiC Gaussian grid')
        sigma, lat, lon = sigma_grid(nlev=nlev, nlat=nlat)
        interfaces = sigma_edges(nlev)
        if not np.array_equal(month_indices(ds['time']) % 12, np.arange(12)):
            raise ValueError('reference ps must be ordered January through December')
        np.testing.assert_allclose(finite(ds['lat']), lat, atol=1e-4, rtol=0)
        np.testing.assert_allclose(finite(ds['lon']), lon, atol=1e-4, rtol=0)
    provenance = [{'file': p.name, 'sha256': sha256(p)} for p in files]
    if checksums is not None:
        observed = {entry['file']: entry['sha256'] for entry in provenance}
        if observed != checksums:
            raise ValueError('input SHA256 checksums differ from ESGF v20160711')
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent, prefix=f'.{output.name}.') as work:
        temporary = Path(work) / output.name
        with nc.Dataset(temporary, 'w', format='NETCDF4_CLASSIC') as out:
            for name, size in zip(('time', 'lev', 'lat', 'lon'), (1980, nlev, nlat, len(lon))):
                out.createDimension(name, size)
            months = np.arange(1850 * 12, 2015 * 12)
            write_cf_grid(out, months, sigma, interfaces, lat, lon, ps)
            o3 = out.createVariable('o3', 'f4', ('time', 'lev', 'lat', 'lon'),
                                    zlib=True, complevel=1, chunksizes=(1, nlev, min(nlat, 16), len(lon)))
            o3.setncatts(dict(standard_name='mole_fraction_of_ozone_in_air', units='mol mol-1',
                             long_name='Ozone mole fraction averaged over each model pressure layer',
                             cell_methods='time: mean'))
            out.source_files = json.dumps(provenance)
            out.reference_ps_sha256 = sha256(pressure_file)
            out.horizontal_remapping = 'conservative'
            out.interpolation = 'area overlap on native pressure levels and Gaussian latitude bands; pressure-overlap means on published plev bounds; deepest source value extended below source grid'
            out.vertical_mapping = 'p=sigma*fixed monthly reference ps; not instantaneous surface pressure'
            out.ozone_experiment = 'repeated 1850 seasonal cycle' if fixed else 'historical 1850-2014'
            position = 0
            for path in files:
                with nc.Dataset(path) as ds:
                    if ds.source_id != 'UReading-CCMI-1-0':
                        raise ValueError(f'{path}: incorrect source_id')
                    variable = ds['vmro3']
                    if variable.dimensions != ('time', 'plev', 'lat', 'lon'):
                        raise ValueError('vmro3 must have dimensions time,plev,lat,lon')
                    if variable.units.strip() not in ('mole mole-1', 'mol mol-1'):
                        raise ValueError('unsupported vmro3 units')
                    punit = ds['plev'].units.strip()
                    if punit not in ('Pa', 'hPa'):
                        raise ValueError('plev must use Pa or hPa')
                    pressure = finite(ds['plev']) * (100 if punit == 'hPa' else 1)
                    bounds = finite(ds[ds['plev'].bounds]) * (100 if punit == 'hPa' else 1)
                    dates = month_indices(ds['time'])
                    np.testing.assert_array_equal(dates, months[position:position + len(dates)])
                    source_lat, source_lon = finite(ds['lat']), finite(ds['lon'])
                    weights = area_weights(ds, lat, lon)
                    for i, date in enumerate(dates):
                        if not fixed or position < 12:
                            o3[position] = regrid(finite(variable, i), pressure, bounds,
                                                  source_lat, source_lon, interfaces, lat, lon, ps[date % 12], weights)
                        else:
                            o3[position] = o3[position % 12]
                        position += 1
            if position != 1980:
                raise ValueError(f'expected 1980 months, got {position}')
        # Publish the closed NetCDF atomically; a concurrent writer may have finished first.
        os.link(temporary, output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('files', nargs=4, type=Path)
    parser.add_argument('--reference-ps', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--levels', type=int, default=10, help='number of stretched model sigma levels (historical experiment: 10)')
    parser.add_argument('--fixed-1850', action='store_true')
    args = parser.parse_args()
    sums_file = Path(__file__).with_name('SHA256SUMS')
    expected = {name: digest for digest, name in (line.split() for line in sums_file.read_text().splitlines())}
    prepare(args.files, args.reference_ps, args.output, args.fixed_1850, expected, args.levels)
