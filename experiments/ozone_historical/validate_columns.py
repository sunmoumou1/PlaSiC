"""Compare source and model ozone columns at every monthly input anchor."""
import argparse
import json
from pathlib import Path
import sys

import netCDF4 as nc
import numpy as np

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from experiments.ozone_historical.common import sha256
from experiments.ozone_historical.prepare_ozone import area_weights, finite, month_indices, sigma_edges

EXPERIMENT = Path(__file__).resolve().parent
# Hydrostatic column for a mole fraction integrated over Pa. Use PlaSiC's
# ozone density at STP (2.14 kg/m³) and standard Earth gravity.
DU_PER_PA = (47.9982 / 28.9644) / (9.80665 * 2.14e-5)


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


def validate(data, reference):
    checksums = {name: digest for digest, name in
                 (line.split() for line in (EXPERIMENT / 'SHA256SUMS').read_text().splitlines())}
    global_means = []
    max_column_error = 0.0
    with nc.Dataset(reference) as ref, nc.Dataset(data / 'ozone_historical.nc') as model:
        if model.reference_ps_sha256 != sha256(reference):
            raise ValueError('reference pressure differs from conversion input')
        if {s['file']: s['sha256'] for s in json.loads(model.source_files)} != checksums:
            raise ValueError('model forcing identifies different source files')
        if getattr(model, 'horizontal_remapping', None) != 'conservative':
            raise ValueError('expected area-conservative ozone input')
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
                lat_bounds = np.deg2rad(finite(source[source['lat'].bounds]))
                lon_bounds = finite(source[source['lon'].bounds])
                area = np.diff(np.sin(lat_bounds), axis=1) * np.diff(lon_bounds, axis=1).T / 720
                np.testing.assert_allclose(area.sum(), 1.)
                # Integrate the original pressure cells down to the same local surface
                # used by the converter to construct each target sigma layer.
                matched_dp = [pressure_thickness(source, p) for p in ps]
                wy, wx = area_weights(source, lat, lon)
                for i, date in enumerate(month_indices(source['time'])):
                    year, m = divmod(int(date), 12)
                    t = len(global_means)
                    if date != 1850 * 12 + t or (model['year'][t], model['month'][t]) != (year, m + 1):
                        raise ValueError(f'{name}: mismatched month at record {i}')
                    raw = finite(source['vmro3'], i)
                    colocated = (wy @ raw) @ wx.T
                    # Each pressure level must retain its native global area mean.
                    np.testing.assert_allclose(colocated.mean(axis=2) @ weights,
                                               (raw * area).sum(axis=(1, 2)), rtol=1e-12, atol=1e-18)
                    matched = np.einsum('kij,kij->ij', colocated, matched_dp[m]) * DU_PER_PA
                    model_column = np.einsum('k,kij->ij', dsigma, finite(model['o3'], t)) * ps[m] * DU_PER_PA
                    largest = np.max(np.abs(model_column - matched))
                    if largest > 1e-4:
                        raise ValueError(f'{year}-{m+1:02}: column mismatch {largest:.6g} DU')
                    global_means.append(weights @ model_column.mean(axis=1))
                    max_column_error = max(max_column_error, largest)
    if len(global_means) != 1980:
        raise ValueError('expected 1980 consecutive monthly anchors')
    print(f'Checked {len(global_means)} months. Global ozone column ranged from '
          f'{min(global_means):.3f} to {max(global_means):.3f} DU. '
          f'Maximum cell difference was {max_column_error:.8f} DU (tolerance 0.0001 DU).')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=EXPERIMENT / 'data')
    parser.add_argument('--reference-ps', type=Path)
    args = parser.parse_args()
    validate(args.data, args.reference_ps or args.data / 'reference/monthly/ps.nc')
