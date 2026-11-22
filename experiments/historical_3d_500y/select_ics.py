#!/usr/bin/env python3
"""Pick two initial conditions from the year-500 ensemble.

The two branches are chosen to sample internal variability rather than two
nearly identical states: A is the initial condition with the warmest final
near-surface air temperature, B the coldest.  The surface temperature is read
from the last ``step`` event of each one-year segment's progress log.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def last_surface_temperature(path: Path) -> float | None:
    value = None
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            line = line.strip()
            if not line:
                continue
            try:
                event = json.loads(line)
            except json.JSONDecodeError:
                continue
            if event.get("event") == "step":
                mean = event.get("global_mean", {})
                if "surface_air_temperature_K" in mean:
                    value = float(mean["surface_air_temperature_K"])
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    entries: list[tuple[str, float]] = []
    for restart in sorted(args.dir.glob("ic_*.restart")):
        if restart.name in {"icA.restart", "icB.restart"}:
            continue
        index = restart.stem.replace("ic_", "")
        progress = args.dir / f"ic_{index}.progress.jsonl"
        if not progress.exists():
            # ic_00 has no dedicated progress log; use the year-500 segment.
            progress = args.dir.parent / "spinup" / "progress_y500.jsonl"
        if not progress.exists():
            continue
        temperature = last_surface_temperature(progress)
        if temperature is None:
            continue
        entries.append((restart.name, temperature))

    if len(entries) < 2:
        raise SystemExit(f"need at least two initial conditions, found {entries}")

    entries.sort(key=lambda item: item[1])
    coldest = entries[0]
    warmest = entries[-1]
    print("initial conditions (name, final surface air temperature K):")
    for name, temperature in entries:
        print(f"  {name:20s} {temperature:.4f}")

    with args.out.open("w", encoding="utf-8") as stream:
        stream.write(f"A {warmest[0]} {warmest[1]:.6f}\n")
        stream.write(f"B {coldest[0]} {coldest[1]:.6f}\n")


if __name__ == "__main__":
    main()
