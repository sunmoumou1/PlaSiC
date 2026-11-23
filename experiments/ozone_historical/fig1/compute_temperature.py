"""Validate complete experiments and write annual global surface-air temperatures."""
import csv
import json
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from experiments.ozone_historical.common import CASES, RUNS

OUTPUT = RUNS / 'data_output/fig1'


def annual_temperature(path, first=1850, last=2014):
    with nc.Dataset(path) as ds:
        var = ds['tas']
        if var.units != 'K' or var.dimensions != ('time', 'lat', 'lon'):
            raise ValueError('tas must be surface air temperature in K on (time,lat,lon)')
        t = ds['time']
        dates = nc.num2date(t[:], t.units, t.calendar)
        months = np.array([d.year * 12 + d.month - 1 for d in dates])
        np.testing.assert_array_equal(months, np.arange(first * 12, (last + 1) * 12))
        bounds = np.asarray(ds[t.bounds][:])
        duration = np.diff(bounds, axis=1).ravel()
        if np.any(duration <= 0) or not np.isfinite(duration).all():
            raise ValueError('invalid monthly time bounds')
        # Bounds must describe whole months, not a partially completed final month.
        bd = nc.num2date(bounds.ravel(), t.units, t.calendar).reshape(-1, 2)
        for index, (a, b) in enumerate(bd):
            if (a.year * 12 + a.month - 1 != months[index] or
                b.year * 12 + b.month - 1 != months[index] + 1 or
                any(d.day != 1 or d.hour or d.minute or d.second for d in (a, b))):
                raise ValueError('time bounds do not cover a complete calendar month')
        values = np.ma.asarray(var[:])
        if np.ma.getmaskarray(values).any() or not np.isfinite(values).all():
            raise ValueError('tas contains missing/nonfinite values')
        nodes, weights = np.polynomial.legendre.leggauss(len(ds['lat']))
        lat = np.asarray(ds['lat'][:])
        order = np.argsort(lat)
        np.testing.assert_allclose(np.sin(np.deg2rad(lat[order])), nodes, atol=1e-6, rtol=0)
        monthly = np.average(np.asarray(values)[:, order].mean(axis=-1, dtype=np.float64), weights=weights, axis=-1)
        years = np.arange(first, last + 1)
        annual = np.array([np.average(monthly[12*i:12*i+12], weights=duration[12*i:12*i+12])
                           for i in range(len(years))])
    return years, annual


def antarctic_spring(path, years):
    """Area mean south of 60°S in the top sigma layer during September–November."""
    means = []
    with nc.Dataset(path) as ds:
        field = ds['ta']
        if field.dimensions != ('time', 'lev', 'lat', 'lon') or field.units != 'K':
            raise ValueError(f'{path}: expected air temperature in K on native model layers')
        lev = np.asarray(ds['lev'][:])
        if not np.all(np.diff(lev) > 0) or not np.isclose(lev[0], .0383, atol=1e-5):
            raise ValueError(f'{path}: expected the historical L10 top layer')
        nodes, weights = np.polynomial.legendre.leggauss(len(ds['lat']))
        order = np.argsort(ds['lat'][:])
        np.testing.assert_allclose(np.sin(np.deg2rad(ds['lat'][:][order])), nodes, atol=1e-6, rtol=0)
        edges = np.r_[-1., -1 + np.cumsum(weights)]
        cap = np.maximum(0., np.minimum(edges[1:], np.sin(np.deg2rad(-60.))) - edges[:-1])
        time = ds['time']
        dates = nc.num2date(time[:], time.units, time.calendar)
        np.testing.assert_array_equal([d.year * 12 + d.month - 1 for d in dates],
                                      np.arange(years[0] * 12, (years[-1] + 1) * 12))
        bounds = np.asarray(ds[time.bounds][:])
        for i, year in enumerate(years):
            season = slice(i * 12 + 8, i * 12 + 11)
            edges = nc.num2date(bounds[season], time.units, time.calendar)
            for month, (start, end) in enumerate(edges, start=9):
                if ((start.year, start.month, end.year, end.month) != (year, month, year, month + 1)
                        or any(d.day != 1 or d.hour or d.minute or d.second for d in (start, end))):
                    raise ValueError(f'{path}: incomplete spring monthly bounds')
            temp = np.ma.asarray(field[season, 0])
            if np.any(np.ma.getmaskarray(temp)) or not np.isfinite(temp).all():
                raise ValueError(f'{path}: missing top-layer air temperature')
            zonal = np.asarray(temp).mean(axis=-1, dtype=np.float64)
            monthly = np.average(zonal[:, order], weights=cap, axis=-1)
            means.append(np.average(monthly, weights=np.diff(bounds[season], axis=1).ravel()))
    return means


def main():
    series = {}
    spring = {}
    for case in CASES:
        case_dir = RUNS / case
        if (case_dir / 'monthly/tas.nc').is_file():
            if (case_dir / '1850').exists():
                raise ValueError(f'{case_dir}: mixed continuous and checkpointed output')
            years, values = annual_temperature(case_dir / 'monthly/tas.nc')
            if (case_dir / 'monthly/ta.nc').exists():
                spring[case] = np.array(antarctic_spring(case_dir / 'monthly/ta.nc', years))
        else:
            years = np.arange(1850, 2015)
            values = np.array([annual_temperature(case_dir / str(year) / 'monthly/tas.nc',
                                                  year, year)[1][0] for year in years])
            if (case_dir / '1850/monthly/ta.nc').exists():
                spring[case] = np.array([antarctic_spring(case_dir / str(y) / 'monthly/ta.nc', [y])[0]
                                         for y in years])
        series[case] = values
        print(f'{case}: min={values.min():.6f}, max={values.max():.6f}, NaN fraction={np.isnan(values).mean():.3f}, units=K')
    OUTPUT.mkdir(parents=True, exist_ok=True)
    anomalies = {c: v - v[(years >= 1850) & (years <= 1900)].mean() for c, v in series.items()}
    summary = {
        'temperature_variable': 'tas',
        'units': 'K',
        'anomaly_baseline': 'Each series uses its own 1850-1900 mean.',
        'period': [int(years[0]), int(years[-1])],
        'experiments': {
            case: {
                'years': len(years),
                'baseline_1850_1900_K': float(values[years <= 1900].mean()),
                'mean_anomaly_2005_2014_K': float(anomalies[case][years >= 2005].mean()),
                'anomaly_2014_K': float(anomalies[case][-1]),
                'temperature_2014_K': float(values[-1]),
            } for case, values in series.items()
        },
    }
    contrasts = {f'co2_ozone_minus_{case}_K': series['co2_ozone'] - series[case]
                 for case in ('co2_only', 'ozone_only')}
    summary['paired_surface_temperature_differences_K'] = {
        name: {'mean_1850_1900': float(values[years <= 1900].mean()),
               'mean_2005_2014': float(values[years >= 2005].mean()), 'year_2014': float(values[-1])}
        for name, values in contrasts.items()
    }
    if spring:
        if set(spring) != set(CASES):
            raise ValueError('top-layer temperature is missing for a case')
        recent = years >= 1980
        polar = {}
        for case, values in spring.items():
            polar[case] = {f'mean_{a}_{b}_K': float(values[(years >= a) & (years <= b)].mean())
                           for a, b in ((1960, 1979), (1980, 1999), (2000, 2014))}
            polar[case]['slope_1980_2014_K_per_decade'] = float(10 * np.polyfit(years[recent], values[recent], 1)[0])
        summary['antarctic_spring_top_layer'] = {
            'field': 'ta', 'sigma_interfaces': [0., .0766],
            'domain': '60–90°S; September–November; area and month-length weighted',
            'experiments': polar,
        }
    (OUTPUT / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    with (OUTPUT / 'annual_temperature.csv').open('w', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(['year', *[c + '_K' for c in CASES], *[c + '_anomaly_K' for c in CASES],
                         *contrasts, *[c + '_antarctic_son_top_K' for c in spring]])
        for i, year in enumerate(years):
            writer.writerow([year, *[series[c][i] for c in CASES], *[anomalies[c][i] for c in CASES],
                             *[v[i] for v in contrasts.values()], *[v[i] for v in spring.values()]])


if __name__ == '__main__':
    main()
