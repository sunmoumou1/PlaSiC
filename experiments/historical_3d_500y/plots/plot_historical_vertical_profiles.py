#!/usr/bin/env python3
"""Rebuild Figure 4 with large Arial text and explicit pressure coordinates.

Usage: python historical_3d_500y/plots/plot_historical_vertical_profiles.py

Figure contract
---------------
Claim: the historical simulations have different vertical temperature-change
structures in the atmosphere and ocean.
Evidence: (a) global atmospheric temperature changes at all 10 sigma levels;
(b) global ocean temperature changes at all 16 depth levels, for members A/B.
Archetype: quantitative grid with two parallel profile panels.
Both members are displayed individually; no confidence interval is inferred.
Changes retain the original equal-month 2005-2014 minus 1855-1859 control
calculation, Gaussian latitude weights and valid-ocean-cell normalization.
Atmospheric pressure is sigma times the control global-mean surface pressure.
The positive pressure axis remains logarithmic, with explicit hPa tick values.
No explanatory callouts, interpolation or smoothing are added.
Exports: PNG only, with large Arial text at 450 dpi.
Inputs: archived profiles when available, otherwise local monthly NetCDF runs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.ticker import FixedLocator, FuncFormatter, MultipleLocator, NullLocator


ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "figures/data/historical_vertical_profiles"
IMAGES = ROOT / "figures"
STEM = "fig4_vertical_ocean_profiles"
BRANCHES = ("A", "B")
CONTROL_PERIOD = (1855, 1859)
HISTORICAL_PERIOD = (2005, 2014)
PRESSURE_TICKS = [50, 70, 100, 200, 300, 400, 600, 800, 1000]
INK = "#222222"


def register_arial() -> None:
    for name in ("Arial.ttf", "Arial Bold.ttf"):
        path = Path("/System/Library/Fonts/Supplemental") / name
        if path.is_file():
            font_manager.fontManager.addfont(path)
    for weight in ("normal", "bold"):
        try:
            font_manager.findfont(
                font_manager.FontProperties(family="Arial", weight=weight),
                fallback_to_default=False,
            )
        except ValueError as exc:
            raise RuntimeError("Install Arial before generating this figure.") from exc


def validate_profiles(metrics: dict) -> dict:
    result = {}
    for branch in BRANCHES:
        item = {key: np.asarray(metrics[branch][key], dtype=np.float64) for key in (
            "pressure_hpa", "delta_ta", "ocean_depth_m", "delta_thetao"
        )}
        for coordinate, temperature in (
            ("pressure_hpa", "delta_ta"), ("ocean_depth_m", "delta_thetao")
        ):
            y, x = item[coordinate], item[temperature]
            if y.ndim != 1 or y.size < 2 or x.shape != y.shape:
                raise ValueError("Every level requires one temperature change.")
            if not np.isfinite(x).all() or not np.isfinite(y).all():
                raise ValueError("Missing profile values are not silently removed.")
            if not np.all(y > 0) or not np.all(np.diff(y) > 0):
                raise ValueError("Pressure and depth must be positive and increasing.")
        result[branch] = item
    return result


def read_archive() -> dict:
    metadata = json.loads((DATA / "provenance.json").read_text())
    metrics = {branch: {} for branch in BRANCHES}
    for archive in metadata["archives"]:
        path = DATA / archive["file"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != archive["sha256"]:
            raise ValueError(f"Archive checksum mismatch: {path.name}")
        with path.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        if set(row["branch"] for row in rows) != set(BRANCHES):
            raise ValueError("The archive must contain exactly members A and B.")
        for branch in BRANCHES:
            selected = [row for row in rows if row["branch"] == branch]
            metrics[branch][archive["coordinate_key"]] = [
                float(row["vertical_coordinate"]) for row in selected
            ]
            metrics[branch][archive["temperature_key"]] = [
                float(row["temperature_change_K"]) for row in selected
            ]
    return validate_profiles(metrics)


def period_profile(path: Path, variable: str, period: tuple[int, int]):
    """Reproduce the original monthly Gaussian-weighted global mean."""
    import netCDF4 as nc4

    start, end = period
    with nc4.Dataset(path) as ds:
        time = ds["time"]
        dates = nc4.num2date(time[:], time.units, time.calendar)
        months = np.array([d.year * 12 + d.month - 1 for d in dates])
        selected = (months >= start * 12) & (months < (end + 1) * 12)
        if not np.array_equal(months[selected], np.arange(start * 12, (end + 1) * 12)):
            raise ValueError(f"Incomplete monthly coverage in {path.name}.")
        values = np.ma.filled(np.ma.asarray(ds[variable][selected], dtype=np.float64), np.nan)
        lat = np.asarray(ds["lat"][:], dtype=np.float64)
        coordinate = "lev" if variable == "ta" else "depth"
        vertical = np.asarray(ds[coordinate][:], dtype=np.float64) if variable != "ps" else None
        expected_units = "Pa" if variable == "ps" else "K"
        if ds[variable].units != expected_units:
            raise ValueError(f"Expected {expected_units} for {variable}.")
    nodes, weights = np.polynomial.legendre.leggauss(lat.size)
    closest = np.argmin(np.abs(nodes[None, :] - np.sin(np.deg2rad(lat))[:, None]), axis=1)
    weights = weights[closest]
    weights /= weights.sum()
    valid = np.isfinite(values)
    if valid.all():
        monthly = np.tensordot(values.mean(axis=-1), weights, axes=([-1], [0]))
    else:
        if variable != "thetao":
            raise ValueError("Only the existing ocean land/depth mask is allowed.")
        weighted = weights[None, :, None] * valid
        denominator = weighted.sum(axis=(-2, -1))
        if np.any(denominator == 0):
            raise ValueError("A depth level contains no valid ocean cells.")
        monthly = np.nansum(np.where(valid, values, 0) * weighted, axis=(-2, -1)) / denominator
    return monthly.mean(axis=0), vertical


def read_raw_profiles(raw_root: Path) -> dict:
    metrics = {}
    for branch in BRANCHES:
        control = raw_root / f"control_{branch}" / "monthly"
        historical = raw_root / f"historical_{branch}" / "monthly"
        ta_control, sigma = period_profile(control / "ta.nc", "ta", CONTROL_PERIOD)
        ta_history, history_sigma = period_profile(historical / "ta.nc", "ta", HISTORICAL_PERIOD)
        ocean_control, depth = period_profile(control / "thetao.nc", "thetao", CONTROL_PERIOD)
        ocean_history, history_depth = period_profile(historical / "thetao.nc", "thetao", HISTORICAL_PERIOD)
        surface_pressure, _ = period_profile(control / "ps.nc", "ps", CONTROL_PERIOD)
        if not np.array_equal(sigma, history_sigma) or not np.array_equal(depth, history_depth):
            raise ValueError("Control and historical runs must share vertical levels.")
        metrics[branch] = {
            "pressure_hpa": sigma * surface_pressure / 100,
            "delta_ta": ta_history - ta_control,
            "ocean_depth_m": depth,
            "delta_thetao": ocean_history - ocean_control,
        }
    return validate_profiles(metrics)


def make_figure(metrics: dict, *, publication: bool = False):
    width, height = (183 / 25.4, 96 / 25.4) if publication else (12.4, 6.5)
    body, label, title, panel = (7, 7, 7, 8) if publication else (15, 16, 17, 19)
    line_width, marker_size = (1.1, 3.1) if publication else (2.2, 5.8)
    style = {
        "font.family": "Arial", "font.sans-serif": ["Arial"],
        "font.size": body, "text.color": INK,
        "axes.labelsize": label, "axes.labelcolor": INK,
        "axes.edgecolor": INK, "axes.linewidth": 0.7,
        "axes.spines.right": False, "axes.spines.top": False, "axes.grid": False,
        "xtick.labelsize": body, "ytick.labelsize": body,
        "xtick.color": INK, "ytick.color": INK,
        "xtick.direction": "out", "ytick.direction": "out",
        "legend.frameon": False, "legend.fontsize": body,
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "pdf.use14corefonts": False, "text.usetex": False,
        "figure.facecolor": "white", "savefig.facecolor": "white",
    }
    with plt.rc_context(style):
        fig = plt.figure(figsize=(width, height))
        atmosphere = fig.add_axes([0.095, 0.22, 0.365, 0.67])
        ocean = fig.add_axes([0.595, 0.22, 0.365, 0.67])
        for ax, temperature, coordinate in (
            (atmosphere, "delta_ta", "pressure_hpa"),
            (ocean, "delta_thetao", "ocean_depth_m"),
        ):
            ax.axvline(0, color="#888888", linestyle="--", linewidth=0.8, zorder=0)
            for branch, colour, marker, linestyle in (
                ("A", "#D55E00", "o", "-"), ("B", "#0072B2", "s", "--")
            ):
                item = metrics[branch]
                ax.plot(item[temperature], item[coordinate], color=colour,
                        linewidth=line_width, marker=marker, markersize=marker_size,
                        markeredgewidth=0.6, linestyle=linestyle, label=f"Branch {branch}")
            ax.tick_params(length=3.5, width=0.7, pad=4)
            ax.set_xlabel("Temperature change (K)", labelpad=7)

        atmosphere.set_yscale("log")
        pressures = np.concatenate([metrics[b]["pressure_hpa"] for b in BRANCHES])
        atmosphere.set_ylim(max(1050, pressures.max() * 1.04), min(30, pressures.min() / 1.08))
        atmosphere.yaxis.set_major_locator(FixedLocator(PRESSURE_TICKS))
        atmosphere.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
        atmosphere.yaxis.set_minor_locator(NullLocator())
        atmosphere.set_ylabel("Pressure (hPa)", labelpad=7)
        atmospheric_changes = np.concatenate([metrics[b]["delta_ta"] for b in BRANCHES])
        atmosphere.set_xlim(np.floor((atmospheric_changes.min() - 0.15) * 2) / 2,
                            np.ceil((atmospheric_changes.max() + 0.15) * 2) / 2)
        atmosphere.xaxis.set_major_locator(MultipleLocator(1))

        ocean_changes = np.concatenate([metrics[b]["delta_thetao"] for b in BRANCHES])
        depths = np.concatenate([metrics[b]["ocean_depth_m"] for b in BRANCHES])
        ocean.set_xlim(min(-0.03, ocean_changes.min() - 0.02),
                       np.ceil((ocean_changes.max() + 0.02) * 10) / 10 + 0.02)
        ocean.set_ylim(np.ceil(depths.max() / 1000) * 1000, 0)
        ocean.yaxis.set_major_locator(MultipleLocator(1000))
        ocean.xaxis.set_major_locator(MultipleLocator(0.1))
        ocean.set_ylabel("Depth (m)", labelpad=7)

        handles, labels = atmosphere.get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.53, 0.025),
                   ncol=2, handlelength=2.3, columnspacing=2.4,
                   fontsize=7 if publication else 14.5)
        for letter, letter_x, title_x, heading in (
            ("a", 0.024, 0.095, "Atmospheric vertical structure"),
            ("b", 0.524, 0.595, "Ocean temperature change"),
        ):
            fig.text(letter_x, 0.965, letter, fontsize=panel, weight="bold", va="top")
            fig.text(title_x, 0.961, heading, fontsize=title, va="top")
    return fig


def render_vertical_profiles(metrics: dict | None = None, output_dir: Path = IMAGES) -> None:
    register_arial()
    if metrics is None:
        metrics = read_archive() if (DATA / "provenance.json").is_file() else read_raw_profiles(ROOT)
    else:
        metrics = validate_profiles(metrics)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    web = make_figure(metrics)
    web.savefig(output_dir / f"{STEM}.png", dpi=450)
    plt.close(web)
    print(f"Wrote {STEM}.png (450 dpi, large Arial text) to {output_dir}.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=IMAGES)
    parser.add_argument("--raw-root", type=Path, help="Recompute from the original monthly NetCDF output.")
    args = parser.parse_args()
    metrics = read_raw_profiles(args.raw_root) if args.raw_root else None
    render_vertical_profiles(metrics, args.output_dir)


if __name__ == "__main__":
    main()
