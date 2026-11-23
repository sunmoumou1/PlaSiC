"""Compare source and model ozone columns at every monthly input anchor."""
import argparse
import json
from pathlib import Path
import sys

import netCDF4 as nc
import numpy as np
from scipy.interpolate import RegularGridInterpolator

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from experiments.ozone_historical.common import RUNS, sha256
from experiments.ozone_historical.prepare_ozone import area_weights, finite, month_indices, sigma_edges

EXPERIMENT = Path(__file__).resolve().parents[1]
# Hydrostatic column for a mole fraction integrated over Pa. Use PlaSiC's
# ozone density at STP (2.14 kg/m³) and standard Earth gravity.
DU_PER_PA = (47.9982 / 28.9644) / (9.80665 * 2.14e-5)


def horizontal(field, lat, lon, target_lat, target_lon):
    """Bilinear sampling with a periodic longitude and constant polar end caps."""
    extended = np.concatenate((field, field[:, :1]), axis=1)
    interp = RegularGridInterpolator((lat, np.r_[lon, lon[0] + 360]), extended)
    return interp(tuple(np.meshgrid(np.clip(target_lat, lat.min(), lat.max()),
                                    (target_lon - lon[0]) % 360 + lon[0], indexing='ij')))


def pressure_thickness(dataset, surface):
    # Integrate each published cell directly down to local ground. Extend the
    # deepest ozone value below its lower boundary, as in the production map.
    bounds = finite(dataset[dataset['plev'].bounds]) * 100
    if dataset['plev'].units != 'hPa':
        raise ValueError('expected the published pressure bounds in hPa')
    tops = np.maximum(0., bounds.min(axis=1))
    bottoms = bounds.max(axis=1)
    bottoms[np.argmax(bottoms)] = np.inf
    lower = np.minimum(bottoms[:, None, None], surface)
    return np.maximum(0., lower - tops[:, None, None])


def compute(data, reference, output):
    checksums = {name: digest for digest, name in
                 (line.split() for line in (EXPERIMENT / 'SHA256SUMS').read_text().splitlines())}
    rows, zonal = [], []
    with nc.Dataset(reference) as ref, nc.Dataset(data / 'ozone_historical.nc') as model:
        if model.reference_ps_sha256 != sha256(reference):
            raise ValueError('reference pressure differs from conversion input')
        if {s['file']: s['sha256'] for s in json.loads(model.source_files)} != checksums:
            raise ValueError('model forcing identifies different source files')
        ps = finite(ref['ps'])
        lat, lon = finite(model['lat']), finite(model['lon'])
        for axis in ('lat', 'lon'):
            np.testing.assert_allclose(ref[axis][:], model[axis][:], atol=1e-4, rtol=0)
        weights = np.polynomial.legendre.leggauss(len(lat))[1][::-1] / 2
        dsigma = np.diff(sigma_edges(len(model['lev'])))
        for name, digest in sorted(checksums.items()):
            if sha256(data / name) != digest:
                raise ValueError(f'{name}: source checksum mismatch')
            with nc.Dataset(data / name) as source:
                sy, sx = finite(source['lat']), finite(source['lon'])
                lat_bounds = np.deg2rad(finite(source[source['lat'].bounds]))
                lon_bounds = finite(source[source['lon'].bounds])
                area = np.diff(np.sin(lat_bounds), axis=1) * np.diff(lon_bounds, axis=1).T / 720
                np.testing.assert_allclose(area.sum(), 1.)
                # Native columns use the reference pressure sampled on the
                # source grid. Matched columns use exactly the target ps.
                raw_dp = [pressure_thickness(source, horizontal(p[::-1], lat[::-1], lon, sy, sx)) for p in ps]
                matched_dp = [pressure_thickness(source, p) for p in ps]
                conservative = getattr(model, 'horizontal_remapping', 'bilinear') == 'conservative'
                if conservative:
                    wy, wx = area_weights(source, lat, lon)
                for i, date in enumerate(month_indices(source['time'])):
                    year, m = divmod(int(date), 12)
                    t = len(rows)
                    if date != 1850 * 12 + t or (model['year'][t], model['month'][t]) != (year, m + 1):
                        raise ValueError(f'{name}: mismatched month at record {i}')
                    raw = finite(source['vmro3'], i)
                    native = np.einsum('kij,kij->ij', raw, raw_dp[m]) * DU_PER_PA
                    if conservative:
                        colocated = (wy @ raw) @ wx.T
                        # Every pressure level must retain its native global area mean.
                        np.testing.assert_allclose(colocated.mean(axis=2) @ weights,
                                                   (raw * area).sum(axis=(1, 2)), rtol=1e-12, atol=1e-18)
                    else:
                        colocated = horizontal(raw.transpose(1, 2, 0), sy, sx, lat, lon).transpose(2, 0, 1)
                    matched = np.einsum('kij,kij->ij', colocated, matched_dp[m]) * DU_PER_PA
                    model_column = np.einsum('k,kij->ij', dsigma, finite(model['o3'], t)) * ps[m] * DU_PER_PA
                    largest = np.max(np.abs(model_column - matched))
                    if largest > 1e-4:
                        raise ValueError(f'{year}-{m+1:02}: column mismatch {largest:.6g} DU')
                    rows.append((year, m + 1, (area * native).sum(), weights @ matched.mean(axis=1),
                                 weights @ model_column.mean(axis=1), largest))
                    zonal.append(model_column.mean(axis=1))
    if len(rows) != 1980:
        raise ValueError('expected 1980 consecutive monthly anchors')
    output.mkdir(parents=True, exist_ok=True)
    np.savetxt(output / 'monthly_columns.csv', rows, delimiter=',', comments='',
               header='year,month,native_source_DU,colocated_source_DU,model_DU,max_cell_difference_DU',
               fmt=['%d', '%d', '%.9f', '%.9f', '%.9f', '%.9f'])
    climatology = np.reshape(zonal, (165, 12, len(lat))).mean(axis=0).T
    np.savetxt(output / 'zonal_seasonal_cycle.csv', np.column_stack((lat, climatology)), delimiter=',',
               header='latitude_deg,' + ','.join(f'month_{m:02}_DU' for m in range(1, 13)),
               comments='', fmt='%.7f')
    print(f'Global model column {np.min(rows, axis=0)[4]:.3f} to {np.max(rows, axis=0)[4]:.3f} DU; '
          f'zonal climatology {climatology.min():.2f} to {climatology.max():.2f} DU; '
          f'maximum cell difference {max(r[-1] for r in rows):.8f} DU; missing fraction 0.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=EXPERIMENT / 'data')
    parser.add_argument('--reference-ps', type=Path)
    parser.add_argument('--output', type=Path, default=RUNS / 'data_output/fig2')
    args = parser.parse_args()
    compute(args.data, args.reference_ps or args.data / 'reference/monthly/ps.nc', args.output)
