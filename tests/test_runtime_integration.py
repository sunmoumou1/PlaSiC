"""Exercise the serial grid/restart compatibility matrix from make test-runtime."""

import argparse
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("t21", "t21-l12", "t31", "t42", "surface", "earth-t21", "earth-t42"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        output = Path(tmp) / "final.restart"

        def run(model, *options, expected_error=None):
            output.unlink(missing_ok=True)
            result = subprocess.run(
                [model, "--co2-ppm", "280", "--steps", "1", "--output", output, *options],
                capture_output=True, text=True, timeout=120,
            )
            if expected_error is None:
                if result.returncode or not output.is_file() or not output.stat().st_size:
                    raise RuntimeError(f"{model} failed to write a restart: {result.stderr}")
            elif result.returncode == 0 or expected_error not in result.stderr:
                raise RuntimeError(f"{model} expected {expected_error!r}: {result.stderr}")

        for model in (args.t21, args.t21_l12, args.t31, args.t42):
            run(model, "--data-dir", args.surface)
            initial = output.replace(Path(tmp) / "initial.restart")
            run(model, "--restart", initial)
        run(args.t21, "--restart", args.earth_t21)
        run(args.t42, "--restart", args.earth_t42)
        for model in (args.t21_l12, args.t31, args.t42):
            run(model, "--restart", args.earth_t21, expected_error="restart grid")


if __name__ == "__main__":
    main()
