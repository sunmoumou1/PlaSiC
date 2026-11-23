"""Run the three 1850–2014 forcing experiments from one common spin-up state."""
from __future__ import annotations
import argparse
from collections import deque
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from contextlib import contextmanager
import sys
from tempfile import mkdtemp
import netCDF4 as nc
import numpy as np
import json
import os
from pathlib import Path
import shlex
import subprocess
import time

STEP_MINUTES = 45
STEPS_PER_DAY = 1440 // STEP_MINUTES
HISTORICAL_CLOCK = ('--step-minutes', str(STEP_MINUTES))
DEFAULT_LAUNCHER = 'mpirun --bind-to none -np 2'
CHECKPOINT_FILES = ('final.restart', 'monthly/tas.nc')
ROOT = Path(__file__).resolve().parents[2]
CO2 = ROOT / 'src/data/forcing/CO2/mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_CMIP_UoM-CMIP-1-2-0_gn-15x360deg_000001-201412.csv'


sys.path.insert(0, str(ROOT))
from experiments.ozone_historical.common import CASES, RUNS, sha256 as digest
from experiments.historical_3d_500y.read_restart import read_records
from experiments.ozone_historical.fig1.compute_temperature import annual_temperature


def year_start(year):
    """Model step at 00:00 on January 1. PlaSiC counts Gregorian year zero."""
    leap_days = (year + 3) // 4 - (year + 99) // 100 + (year + 399) // 400
    return (365 * year + leap_days) * STEPS_PER_DAY


def historical_clock(model):
    """Pre-override T21 builds used 45 minutes; newer builds need an explicit step."""
    usage = subprocess.check_output([model, '--help'], text=True)
    if '--ozone-forcing' not in usage:
        raise ValueError(f'{model}: monthly external ozone support is required')
    return HISTORICAL_CLOCK if '--step-minutes' in usage else ()


def command(model, restart, ozone, directory, case, end_year=2014, mpirun=DEFAULT_LAUNCHER,
            first_year=1850, *, clock=HISTORICAL_CLOCK):
    steps = year_start(end_year + 1) - year_start(first_year)
    result = [*shlex.split(mpirun), str(model),
              '--restart', str(restart),
              '--steps', str(steps), *clock, '--ozone-forcing', str(ozone),
              '--monthly-output', str(directory / 'monthly'),
              '--progress', str(directory / 'progress.jsonl'), '--progress-interval', '3200',
              '--output', str(directory / 'final.restart')]
    if first_year == 1850:
        result += ['--reset-nstep', str(year_start(1850))]
    return result + (['--co2-ppm', '284.316853'] if case == 'ozone_only'
                     else ['--co2-forcing', str(CO2)])


@contextmanager
def lock_run(output):
    with (output / '.run.lock').open('a+b') as lock:
        if os.name == 'nt':
            import msvcrt
            if lock.tell() == 0:
                lock.write(b'\0')
                lock.flush()
            lock.seek(0)
            msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
        else:
            import fcntl
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield


def checkpoint_files(directory):
    monthly = sorted(p.relative_to(directory).as_posix() for p in (directory / 'monthly').glob('*.nc'))
    required = {'monthly/' + name + '.nc' for name in ('tas', 'ta', 'co2', 'thetao', 'rsdt', 'rsut', 'rlut')}
    if not required.issubset(monthly):
        raise ValueError(f'{directory}: missing monthly model diagnostics')
    return ('final.restart', *monthly)


def segment(args, case, year, inputs):
    """Resume a verified year or publish a newly completed one with its restart."""
    case_dir = args.output / case
    case_dir.mkdir(exist_ok=True)
    if (case_dir / 'monthly').exists():
        raise ValueError(f'{case_dir}: contains a continuous run; choose a separate output directory')
    destination = case_dir / str(year)
    restart = args.restart if year == 1850 else case_dir / str(year - 1) / 'final.restart'
    ozone = args.fixed_ozone if case == 'co2_only' else args.historical_ozone
    signature = dict(inputs=inputs, case=case, year=year, restart=digest(restart))
    if destination.exists():
        previous = json.loads((destination / 'run.json').read_text())
        if previous['signature'] != signature:
            raise ValueError(f'{destination}: inputs differ from completed segment')
        checksums = previous.get('output_sha256')
        version = previous.get('checkpoint_version', 1)
        if version not in (1, 2):
            raise ValueError(f'{destination}: unsupported checkpoint version {version}')
        files = CHECKPOINT_FILES if version == 1 else checkpoint_files(destination)
        if not isinstance(checksums, dict) or checksums.keys() != set(files):
            raise ValueError(f'{destination}: expected checksums for {", ".join(files)}')
        for name in files:
            if digest(destination / name) != checksums[name]:
                raise ValueError(f'{destination / name}: checksum changed')
        return
    # The next year may consume this restart only after all monthly files pass validation.
    directory = Path(mkdtemp(prefix=f'.{year}-', dir=case_dir))
    (directory / 'monthly').mkdir()
    cmd = command(args.model, restart, ozone, directory, case, year, args.mpirun, year, clock=args.clock)
    manifest = dict(command=cmd, signature=signature, cwd=str(ROOT / 'src'), start_unix=time.time())
    (directory / 'run.json').write_text(json.dumps(manifest, indent=2))
    with (directory / 'run.log').open('w') as log:
        subprocess.run(cmd, cwd=ROOT / 'src', stdout=log, stderr=subprocess.STDOUT, check=True)
    # The shared historical reader preserves each word as float32 bits.
    step = read_records(directory / 'final.restart')['nstep'].view('<i4').item()
    if step != year_start(year + 1):
        raise ValueError(f'{directory}: final restart clock differs from requested end')
    tas = directory / 'monthly/tas.nc'
    _, temperatures = annual_temperature(tas, year, year)
    with nc.Dataset(tas) as ds:
        if case == 'ozone_only':
            expected_co2 = f'fixed CO2: {np.float32(284.316853):.9g} ppmv'
        else:
            expected_co2 = str(CO2)
        if ds.co2_forcing_source != expected_co2:
            raise ValueError(f'{tas}: incorrect CO₂ forcing source {ds.co2_forcing_source}')
        if ds.ozone_forcing_source != str(ozone):
            raise ValueError(f'{tas}: incorrect ozone forcing source')
    manifest.update(checkpoint_version=2, finish_unix=time.time(), annual_tas_K=float(temperatures[0]),
                    output_sha256={f: digest(directory / f) for f in checkpoint_files(directory)})
    (directory / 'run.json').write_text(json.dumps(manifest, indent=2))
    directory.rename(destination)
    print(f'{case} {year}: complete ({temperatures[0]:.5f} K)', flush=True)


def run(args):
    for p in (args.model, args.restart, args.historical_ozone, args.fixed_ozone, CO2):
        if not p.is_file():
            raise FileNotFoundError(p)
    args.clock = historical_clock(args.model)
    for path, expected in ((args.historical_ozone, 'historical 1850-2014'),
                           (args.fixed_ozone, 'repeated 1850 seasonal cycle')):
        with nc.Dataset(path) as ds:
            if getattr(ds, 'ozone_experiment', None) != expected:
                raise ValueError(f'{path}: expected {expected}')
            np.testing.assert_array_equal(ds['year'][:], np.repeat(np.arange(1850, 2015), 12))
            np.testing.assert_array_equal(ds['month'][:], np.tile(np.arange(1, 13), 165))
    inputs = {str(p): digest(p) for p in (args.model, args.restart, args.historical_ozone, args.fixed_ozone, CO2)}
    args.output.mkdir(parents=True, exist_ok=True)
    with lock_run(args.output):
        # Use separate processes: netCDF-C readers must not share Python threads.
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            ready = deque((case, 1850) for case in CASES)
            active = {}
            while ready or active:
                while ready and len(active) < args.jobs:
                    case, year = ready.popleft()
                    task = pool.submit(segment, args, case, year, inputs)
                    active[task] = (case, year)
                done, _ = wait(active, return_when=FIRST_COMPLETED)
                for task in done:
                    task.result()
                    case, year = active.pop(task)
                    if year < args.end_year:
                        ready.append((case, year + 1))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', type=Path, default=ROOT / 'build/T21_L10_MPI_P2_O16/plasic.x')
    parser.add_argument('--restart', type=Path, default=ROOT / 'experiments/historical_3d_500y/initial_conditions/icA.restart')
    parser.add_argument('--historical-ozone', required=True, type=Path)
    parser.add_argument('--fixed-ozone', required=True, type=Path)
    parser.add_argument('--output', type=Path, default=RUNS)
    parser.add_argument('--end-year', type=int, default=2014)
    parser.add_argument('--jobs', type=int, choices=(1, 2, 3), default=1, help='simultaneous two-rank simulations')
    parser.add_argument('--mpirun', default=DEFAULT_LAUNCHER, help='complete two-rank MPI launch command (e.g. srun -n 2)')
    args = parser.parse_args()
    if not 1850 <= args.end_year <= 2014:
        parser.error('--end-year must be between 1850 and 2014')
    for key in ('model', 'restart', 'historical_ozone', 'fixed_ozone', 'output'):
        setattr(args, key, getattr(args, key).resolve())
    run(args)
