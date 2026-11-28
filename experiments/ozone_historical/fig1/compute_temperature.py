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


def main():
    series = {}
    for case in CASES:
        case_dir = RUNS / case
        if (case_dir / 'monthly/tas.nc').is_file():
            if (case_dir / '1850').exists():
                raise ValueError(f'{case_dir}: mixed continuous and checkpointed output')
            years, values = annual_temperature(case_dir / 'monthly/tas.nc')
        else:
            years = np.arange(1850, 2015)
            values = np.array([annual_temperature(case_dir / str(year) / 'monthly/tas.nc',
                                                  year, year)[1][0] for year in years])
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
    (OUTPUT / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    with (OUTPUT / 'annual_temperature.csv').open('w', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(['year', *[c + '_K' for c in CASES], *[c + '_anomaly_K' for c in CASES],
                         *contrasts])
        for i, year in enumerate(years):
            writer.writerow([year, *[series[c][i] for c in CASES], *[anomalies[c][i] for c in CASES],
                             *[v[i] for v in contrasts.values()]])


if __name__ == '__main__':
    main()
