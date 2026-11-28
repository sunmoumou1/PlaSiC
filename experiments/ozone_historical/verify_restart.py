"""Compare a continuous two-year replay with the two completed joint-forcing years."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import netCDF4 as nc
import numpy as np

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from experiments.ozone_historical.common import RUNS, sha256
from experiments.ozone_historical.run_experiments import command, DEFAULT_LAUNCHER, ROOT


def compare_metadata(reference, other, label):
    if set(reference.ncattrs()) != set(other.ncattrs()):
        raise ValueError(f'{label}: attribute inventory differs')
    for key in reference.ncattrs():
        np.testing.assert_array_equal(reference.getncattr(key), other.getncattr(key),
                                      err_msg=f'{label}: attribute {key} differs')


def compare(continuous, first, second):
    if sha256(continuous / 'final.restart') != sha256(second / 'final.restart'):
        raise ValueError('continuous and segmented end restarts differ')
    filenames = {p.name for p in (first / 'monthly').glob('*.nc')}
    if not filenames or any({p.name for p in (d / 'monthly').glob('*.nc')} != filenames
                            for d in (second, continuous)):
        raise ValueError('monthly file inventories differ')
    for name in sorted(filenames):
        with nc.Dataset(continuous / 'monthly' / name) as full, nc.Dataset(first / 'monthly' / name) as a, \
                nc.Dataset(second / 'monthly' / name) as b:
            for segment in (a, b):
                compare_metadata(full, segment, name)
            if len(full.dimensions['time']) != 24 or len(a.dimensions['time']) != 12 or len(b.dimensions['time']) != 12:
                raise ValueError(f'{name}: expected two complete years')
            if set(full.variables) != set(a.variables) or set(full.variables) != set(b.variables):
                raise ValueError(f'{name}: variable inventories differ')
            for key, var in full.variables.items():
                for i, segment in enumerate((a, b)):
                    recorded = segment[key]
                    compare_metadata(var, recorded, f'{name}/{key}')
                    part = slice(12*i, 12*i+12) if var.dimensions[:1] == ('time',) else slice(None)
                    if var.dimensions != recorded.dimensions or var.dtype != recorded.dtype:
                        raise ValueError(f'{name}/{key}: inconsistent dimensions or type')
                    # Compare the stored representation, including missing values.
                    var.set_auto_maskandscale(False)
                    recorded.set_auto_maskandscale(False)
                    if np.asarray(var[part]).tobytes() != np.asarray(recorded[:]).tobytes():
                        raise ValueError(f'{name}/{key}: continuous and segmented values differ')
    return len(filenames)


def main(args):
    first, second = (args.output / 'co2_ozone' / str(year) for year in (1850, 1851))
    if args.continuous:
        return compare(args.continuous, first, second)
    manifest = json.loads((first / 'run.json').read_text())
    argv = manifest['command']
    model = Path(argv[argv.index('--restart') - 1])
    restart = Path(argv[argv.index('--restart') + 1])
    ozone = Path(argv[argv.index('--ozone-forcing') + 1])
    co2 = Path(argv[argv.index('--co2-forcing') + 1])
    for path in (model, restart, ozone, co2):
        if sha256(path) != manifest['signature']['inputs'][str(path)]:
            raise ValueError(f'{path}: input differs from the segmented run')
    with TemporaryDirectory() as tmp:
        directory = Path(tmp)
        (directory / 'monthly').mkdir()
        # Recover the exact clock arguments from the original segmented run.
        clock = ()
        if '--step-minutes' in argv:
            clock = ('--step-minutes', argv[argv.index('--step-minutes') + 1])
        cmd = command(model, restart, ozone, directory, 'co2_ozone', 1851, args.mpirun, clock=clock)
        cmd[cmd.index('--co2-forcing') + 1] = str(co2)
        subprocess.run(cmd, cwd=ROOT / 'src', check=True)
        return compare(directory, first, second)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=RUNS, help='archived annual case directories')
    parser.add_argument('--mpirun', default=DEFAULT_LAUNCHER)
    parser.add_argument('--continuous', type=Path, help='verify an existing two-year replay directory')
    count = main(parser.parse_args())
    print(f'End restart matches. Compared {count} monthly files across 24 months.')
