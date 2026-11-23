"""Download and verify the ozone inputs, prepare the grid, run three cases and plot the results."""
import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
from tempfile import TemporaryDirectory

import netCDF4 as nc

ROOT = Path(__file__).resolve().parents[2]
if not __package__:
    sys.path.insert(0, str(ROOT))
from experiments.ozone_historical.common import RUNS, sha256
from experiments.ozone_historical.run_experiments import DEFAULT_LAUNCHER, historical_clock, lock_run, year_start

HERE = Path(__file__).resolve().parent
ESGF = ('https://esgf-node.ornl.gov/thredds/fileServer/user_pub_work/input4MIPs/'
        'CMIP6/CMIP/UReading/UReading-CCMI-1-0/atmos/mon/vmro3/gn/v20160711')
INITIAL = ROOT / 'experiments/historical_3d_500y/initial_conditions/icA.restart'


def prepare_inputs(args, data, run):
    checksums = {name: digest for digest, name in
                 (line.split() for line in (HERE / 'SHA256SUMS').read_text().splitlines())}
    files = [data / name for name in checksums]
    for path in files:
        if path.exists():
            if sha256(path) != checksums[path.name]:
                raise ValueError(f'{path}: archived source checksum mismatch')
            continue
        partial = path.with_suffix('.nc.part')
        download = ['curl', '--fail', '--location', '--retry', '3', '--output', partial]
        url = f'{ESGF}/{path.name}'
        resumed = False
        if partial.is_file():
            try:
                run(*download, '--continue-at', '-', url)
                resumed = sha256(partial) == checksums[path.name]
            except subprocess.CalledProcessError:
                # Some replicas cannot serve a range from the requested offset.
                pass
        if not resumed:
            # A complete but corrupt partial cannot be repaired by appending.
            run(*download, url)
            if sha256(partial) != checksums[path.name]:
                raise ValueError(f'{partial}: download checksum mismatch')
        partial.replace(path)

    if args.model:
        model = args.model.resolve()
    else:
        run('make', 'all', 'MPI=1', 'NPRO=2', 'NLAT=32', 'NLEV=10')
        model = ROOT / 'build/T21_L10_MPI_P2_O16/plasic.x'
    reference = data / 'reference'
    ps = args.reference_ps.resolve() if args.reference_ps else reference / 'monthly/ps.nc'
    if not ps.exists() and not args.reference_ps:
        # Publish the reference only after the full pressure year finishes.
        with TemporaryDirectory(dir=data) as work:
            staging = Path(work) / 'reference'
            (staging / 'monthly').mkdir(parents=True)
            run(*shlex.split(args.mpirun), model, '--restart', INITIAL, *historical_clock(model),
                '--co2-ppm', '284.316853',
                '--reset-nstep', str(year_start(1850)), '--steps', str(year_start(1851) - year_start(1850)),
                '--monthly-output', staging / 'monthly', '--output', staging / 'final.restart')
            staging.rename(reference)
    reference_hash = sha256(ps)
    history, fixed = data / 'ozone_historical.nc', data / 'ozone_1850_cycle.nc'
    for path in (history, fixed):
        if not path.exists():
            options = ['--fixed-1850'] if path == fixed else []
            run(sys.executable, HERE / 'prepare_ozone.py', *files, '--reference-ps', ps, '--output', path, *options)
        with nc.Dataset(path) as prepared:
            sources = {s['file']: s['sha256'] for s in json.loads(prepared.source_files)}
            if getattr(prepared, 'horizontal_remapping', 'bilinear') != 'conservative':
                raise ValueError(f'{path}: this pipeline requires area-conservative inputs, use a new data directory')
            if sources != checksums or prepared.reference_ps_sha256 != reference_hash:
                raise ValueError(f'{path}: existing conversion used different inputs')

    return history, fixed, ps, model


def pipeline(args):
    data, runs = args.data.resolve(), args.output.resolve()
    data.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, PLASIC_OZONE_RUNS=str(runs))

    def run(*command):
        subprocess.run(command, cwd=ROOT / 'src', env=env, check=True)

    # Resume partial downloads and publish each conversion with one preparer at a time.
    with lock_run(data):
        history, fixed, ps, model = prepare_inputs(args, data, run)

    run(sys.executable, HERE / 'fig2/compute_columns.py', '--data', data, '--reference-ps', ps)
    run(sys.executable, HERE / 'fig2/plot_columns.py', '--final')
    run(sys.executable, HERE / 'run_experiments.py', '--historical-ozone', history, '--fixed-ozone', fixed,
        '--model', model, '--jobs', str(args.jobs), '--mpirun', args.mpirun)
    run(sys.executable, HERE / 'fig1/compute_temperature.py')
    run(sys.executable, HERE / 'fig1/plot_temperature.py', '--final')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=HERE / 'data')
    parser.add_argument('--output', type=Path, default=RUNS)
    parser.add_argument('--model', type=Path, help='use an existing T21/L10 two-rank executable')
    parser.add_argument('--reference-ps', type=Path, help='use an existing 12-month model pressure climatology')
    parser.add_argument('--jobs', type=int, choices=(1, 2, 3), default=1)
    parser.add_argument('--mpirun', default=DEFAULT_LAUNCHER, help='complete two-rank launch command')
    pipeline(parser.parse_args())
