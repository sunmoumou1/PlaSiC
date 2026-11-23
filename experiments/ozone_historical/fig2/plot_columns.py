"""Plot the monthly climatology of prescribed, vertically integrated ozone."""
import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.colors import BoundaryNorm
import numpy as np

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from experiments.ozone_historical.common import RUNS


def plot(directory, final=False):
    table = np.loadtxt(directory / 'zonal_seasonal_cycle.csv', delimiter=',', skiprows=1)
    lat, cycle = table[:, 0], table[:, 1:]
    if cycle.shape[1] != 12 or not np.isfinite(table).all():
        raise ValueError('expected a complete latitude by calendar-month climatology')
    # Gaussian cell areas determine the latitude edges, including both poles.
    weights = np.polynomial.legendre.leggauss(len(lat))[1]
    edges = np.rad2deg(np.arcsin(np.clip(np.r_[-1., -1 + weights.cumsum()], -1, 1)))
    plt.rcParams.update({'font.family': 'sans-serif',
                         'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
                         'font.size': 8, 'axes.linewidth': .6, 'pdf.fonttype': 42})
    fig, ax = plt.subplots(figsize=(145/25.4, 72/25.4), layout='constrained', dpi=600)
    levels = np.arange(225, 476, 25)
    cmap = plt.get_cmap('viridis')
    field = ax.pcolormesh(np.arange(.5, 13), edges, cycle[np.argsort(lat)], cmap=cmap,
                          norm=BoundaryNorm(levels, cmap.N), rasterized=True)
    fig.colorbar(field, ax=ax, label='Ozone column (DU)', ticks=levels[::2], aspect=25, pad=.025)
    ax.set_xticks(np.arange(1, 13), 'Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec'.split())
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], ['90°S', '60°S', '30°S', '0°', '30°N', '60°N', '90°N'])
    ax.set(xlabel='Calendar month', ylabel='Latitude', ylim=(-90, 90))
    for extension in (('png', 'pdf') if final else ('png',)):
        fig.savefig(directory / f'zonal_seasonal_cycle.{extension}', dpi=600)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=RUNS / 'data_output/fig2')
    parser.add_argument('--final', action='store_true')
    args = parser.parse_args()
    plot(args.output, args.final)
