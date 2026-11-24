"""Verify shared T21/L10 initialization errors are printed once under MPI."""
import argparse
from pathlib import Path
import shlex
import subprocess
from tempfile import TemporaryDirectory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", type=Path)
    parser.add_argument("--mpirun", default="mpirun")
    args = parser.parse_args()
    source = Path(__file__).resolve().parents[1] / "src"
    with TemporaryDirectory() as tmp:
        missing = Path(tmp) / "missing.csv"
        output = Path(tmp) / "final.restart"
        command = [*shlex.split(args.mpirun), "-np", "2", str(args.model.resolve()),
                   "--restart", "data/earth_t21_l10.restart", "--steps", "1", "--output", str(output)]
        for flags, message in [([], "specify exactly one CO2 mode"),
                               (["--co2-forcing", str(missing)], f"cannot open CO2 forcing file '{missing}'")]:
            result = subprocess.run([*command, *flags], cwd=source, capture_output=True, text=True, timeout=60)
            prefix = "PlaSiC initialization failed: "
            if result.returncode == 0 or result.stderr.count(prefix) != 1 or prefix + message not in result.stderr:
                raise AssertionError(result.stderr)
            if output.exists():
                raise AssertionError("Failed initialization produced a restart")


if __name__ == "__main__":
    main()
