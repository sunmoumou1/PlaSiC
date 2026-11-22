#!/usr/bin/env python3
"""Rebuild Figure 6: longwave, absorbed shortwave and net TOA changes.

Usage: python historical_3d_500y/plots/plot_historical_radiation_budget.py
       python historical_3d_500y/plots/plot_historical_radiation_budget.py --refresh-data

Figure contract
---------------
Question: how do PlaSiC and the six reference models partition their historical
TOA radiation changes between outgoing longwave and absorbed shortwave?
Evidence: (a) delta rlut; (b) delta (rsdt - rsut); (c) each model's 2005-2014
net TOA change = delta (rsdt - rsut) - delta rlut. Positive net is downward.
Archetype: quantitative grid; two parallel time series and a model comparison.
All 165 annual values and all eight simulations are retained. Every simulation,
including each PlaSiC branch, uses its own 1850-1859 historical mean as baseline.
Baseline subtraction precedes smoothing and ensemble summaries; no control-run
or shared multi-model baseline is used.
Months have equal weight; global means use Gaussian weights for PlaSiC and
latitude bounds for CMIP6, with uniform longitude weights on both grids.
Time-series display uses centred 9-year means with shortened endpoint windows.
CMIP6 shading is one sample SD across the six smoothed model series; the panel-c
error bar is one sample SD across their unsmoothed decadal means, not a CI.
Both PlaSiC branches are displayed individually without an inferred interval.
Exports: large Arial PNG at 450 dpi in historical_3d_500y/figures/.
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
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator


ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "figures/data/historical_radiation_budget"
IMAGES = ROOT / "figures"
STEM = "fig6_radiation_budget"
YEARS = np.arange(1850, 2015)
BASELINE = (1850, 1859)
BASELINE_TYPE = "Each simulation's own historical run"
LAST_DECADE = (2005, 2014)
CMIP6_MODELS = (
    "BCC-CSM2-MR", "CESM2", "CESM2-WACCM", "FGOALS-f3-L", "GFDL-CM4", "MRI-ESM2-0",
)
PLASIC_MODELS = ("PlaSiC branch A", "PlaSiC branch B")
MODELS = (*CMIP6_MODELS, *PLASIC_MODELS)
FLUXES = ("rlut", "rsdt", "rsut")
QUANTITIES = (*FLUXES, "delta_rlut", "delta_sw", "delta_net")
INK, REFERENCE, SPREAD = "#222222", "#636363", "#D3D7DB"
COLORS = {PLASIC_MODELS[0]: "#D55E00", PLASIC_MODELS[1]: "#0072B2"}


def register_arial() -> None:
    """Require actual Arial regular and bold instead of silently substituting."""
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


def smooth(values: np.ndarray, window: int = 9) -> np.ndarray:
    """Centred mean; normalize truncated endpoint windows by actual counts."""
    kernel = np.ones(window)
    return np.convolve(values, kernel, mode="same") / np.convolve(
        np.ones_like(values), kernel, mode="same"
    )


def validate_data(data: dict) -> dict:
    if set(data) != set(MODELS):
        raise ValueError("Expected all six CMIP6 models and both PlaSiC branches.")
    for model in MODELS:
        item = data[model]
        for key in QUANTITIES:
            values = np.asarray(item[key], dtype=np.float64)
            if values.shape != YEARS.shape or not np.isfinite(values).all():
                raise ValueError(f"Incomplete/non-finite annual {key} for {model}.")
            item[key] = values
        baseline_years = (YEARS >= BASELINE[0]) & (YEARS <= BASELINE[1])
        for variable in FLUXES:
            own_mean = item[variable][baseline_years].mean()
            if not np.isclose(item["baseline"][variable], own_mean, rtol=0, atol=1e-10):
                raise ValueError(
                    f"{model} {variable}: baseline must be its own "
                    f"{BASELINE[0]}-{BASELINE[1]} historical mean."
                )
        expected = {
            "delta_rlut": item["rlut"] - item["baseline"]["rlut"],
            "delta_sw": (item["rsdt"] - item["rsut"])
            - (item["baseline"]["rsdt"] - item["baseline"]["rsut"]),
            "delta_net": item["delta_sw"] - item["delta_rlut"],
        }
        for key, values in expected.items():
            if not np.allclose(item[key], values, rtol=0, atol=1e-10):
                raise ValueError(f"Radiation-budget identity failed: {model} {key}.")
    return data


def monthly_global(paths: list[Path], variable: str, *, gaussian: bool):
    """Read bounded time blocks, checking complete global grids and months."""
    import netCDF4 as nc4

    month_parts, value_parts = [], []
    first_grid = None
    for path in paths:
        with nc4.Dataset(path) as ds:
            field = ds[variable]
            if field.dimensions != ("time", "lat", "lon") or field.units != "W m-2":
                raise ValueError(f"Unexpected dimensions or units in {path.name}.")
            lat, lon = (np.asarray(ds[key][:], dtype=float) for key in ("lat", "lon"))
            if first_grid is None:
                first_grid = lat, lon
            elif not all(np.array_equal(a, b) for a, b in zip(first_grid, (lat, lon))):
                raise ValueError("Split time files must share an identical grid.")
            if not np.allclose(np.diff(lon), 360 / lon.size, rtol=0, atol=1e-4):
                raise ValueError("Expected a complete, uniform global longitude grid.")
            if gaussian:
                nodes, weights = np.polynomial.legendre.leggauss(lat.size)
                indices = np.argmin(abs(nodes[:, None] - np.sin(np.deg2rad(lat))), axis=0)
                if np.unique(indices).size != lat.size or not np.allclose(
                    nodes[indices], np.sin(np.deg2rad(lat)), atol=1e-6, rtol=0
                ):
                    raise ValueError("PlaSiC latitude centres must form a Gaussian grid.")
                weights = weights[indices]
            else:
                bounds = np.asarray(ds["lat_bnds"][:], dtype=float)
                if bounds.shape != (lat.size, 2):
                    raise ValueError("Expected latitude bounds for every CMIP6 row.")
                weights = abs(np.diff(np.sin(np.deg2rad(bounds)), axis=1)[:, 0])
                ordered = np.sort(bounds, axis=1)[np.argsort(lat)]
                contiguous = (np.isclose(ordered[0, 0], -90)
                              and np.isclose(ordered[-1, 1], 90)
                              and np.allclose(ordered[:-1, 1], ordered[1:, 0], atol=1e-5, rtol=0))
                if not contiguous:
                    # The FGOALS files leave 0.2-degree gaps between nominal
                    # cells. Recover full cells only for a complete regular
                    # global grid. No rows/cells are dropped or interpolated.
                    centres = np.sort(lat)
                    step = 180 / lat.size
                    if not np.allclose(centres, np.arange(lat.size) * step - 90 + step / 2,
                                       rtol=0, atol=1e-5):
                        raise ValueError("Cannot recover invalid bounds on an irregular latitude grid.")
                    edges = np.r_[-90, (centres[:-1] + centres[1:]) / 2, 90]
                    recovered = np.diff(np.sin(np.deg2rad(edges)))
                    weights = recovered[np.argsort(np.argsort(lat))]
            weights = weights / weights.sum()
            time = ds["time"]
            dates = nc4.num2date(time[:], time.units, getattr(time, "calendar", "standard"))
            month_parts.append(np.array([d.year * 12 + d.month - 1 for d in dates]))
            values = []
            # CESM2 stores missing_value as float64 on float32 fields. Cast
            # sentinels to the storage type before comparison to avoid the
            # library's unsafe-cast warning and still reject every missing cell.
            field.set_auto_mask(False)
            missing = []
            for attribute in ("_FillValue", "missing_value"):
                if attribute in field.ncattrs():
                    sentinel = np.asarray(field.getncattr(attribute), dtype=field.dtype)
                    missing.extend(np.atleast_1d(sentinel).astype(float))
            for start in range(0, len(dates), 120):
                block = np.asarray(field[start:start + 120], dtype=np.float64)
                if not np.isfinite(block).all() or any(np.any(block == value) for value in missing):
                    raise ValueError("Missing radiation cells are not silently excluded.")
                values.extend(np.asarray(block).mean(axis=-1) @ weights)
            value_parts.append(np.array(values))
    months, values = np.concatenate(month_parts), np.concatenate(value_parts)
    order = np.argsort(months)
    months, values = months[order], values[order]
    if not np.array_equal(months, np.arange(months[0], months[-1] + 1)):
        raise ValueError("Monthly radiation series contains gaps or duplicate months.")
    return months, values, first_grid


def read_raw_data(raw_root: Path, cmip6_root: Path) -> tuple[dict, list[dict]]:
    """Calculate each flux on its own native grid before forming paired changes."""
    data, sources = {}, []
    for model in MODELS:
        annual, baseline, model_grid = {}, {}, None
        for variable in FLUXES:
            gaussian = model in PLASIC_MODELS
            if gaussian:
                branch = model[-1]
                paths = [raw_root / f"historical_{branch}/monthly/{variable}.nc"]
            else:
                paths = sorted(cmip6_root.glob(f"{variable}_Amon_{model}_historical_*.nc"))
                # Mixing realizations would overweight a model and duplicate months.
                if len({path.name.split("_historical_")[1].split("_")[0] for path in paths}) != 1:
                    raise ValueError(f"Expected one realization of {model} {variable}.")
            if not paths or not all(path.is_file() for path in paths):
                raise FileNotFoundError(f"Missing {variable} inputs for {model}.")
            months, values, grid = monthly_global(paths, variable, gaussian=gaussian)
            if model_grid is not None and not all(
                np.array_equal(a, b) for a, b in zip(model_grid, grid)
            ):
                raise ValueError(f"The three radiation fluxes must share a grid: {model}.")
            model_grid = grid
            if not np.array_equal(months, np.arange(1850 * 12, 2015 * 12)):
                raise ValueError(f"Expected every month in 1850-2014 for {model}.")
            annual[variable] = values.reshape(len(YEARS), 12).mean(axis=1)
            selected = (months >= BASELINE[0] * 12) & (months < (BASELINE[1] + 1) * 12)
            expected_months = np.arange(BASELINE[0] * 12, (BASELINE[1] + 1) * 12)
            if not np.array_equal(months[selected], expected_months):
                raise ValueError(f"Incomplete {BASELINE[0]}-{BASELINE[1]} historical baseline.")
            baseline[variable] = float(values[selected].mean())
            for path in paths:
                with path.open("rb") as stream:
                    digest = hashlib.file_digest(stream, "sha256").hexdigest()
                sources.append({"file": path.name, "model": model,
                                "run": "historical",
                                "sha256": digest, "bytes": path.stat().st_size})
        annual["baseline"] = baseline
        annual["delta_rlut"] = annual["rlut"] - baseline["rlut"]
        annual["delta_sw"] = annual["rsdt"] - annual["rsut"] - (baseline["rsdt"] - baseline["rsut"])
        annual["delta_net"] = annual["delta_sw"] - annual["delta_rlut"]
        data[model] = annual
        print(f"Computed all three radiation fluxes: {model}", flush=True)
    return validate_data(data), sources


def write_archive(data: dict, sources: list[dict]) -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    path = DATA / "annual_radiation.csv"
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["model", "year", *QUANTITIES])
        for model in MODELS:
            for index, year in enumerate(YEARS):
                writer.writerow([model, year, *(repr(float(data[model][q][index])) for q in QUANTITIES)])
    metadata = {
        "file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "units": "W m-2", "years": [1850, 2014], "baseline_years": list(BASELINE),
        "baseline_type": BASELINE_TYPE,
        "baselines": {model: data[model]["baseline"] for model in MODELS},
        "annual_mean": "Arithmetic mean of 12 monthly global means, not day weighted.",
        "spatial_mean": "PlaSiC Gaussian quadrature; CMIP6 sin(latitude bounds); equal longitudes.",
        "latitude_bound_repair": "FGOALS-f3-L: nominal 0.8-degree cells on a complete 1-degree grid leave gaps; use centre midpoints and poles to recover contiguous cells. No data removed.",
        "net_identity": "delta_net = delta_sw - delta_rlut; delta_sw = delta(rsdt - rsut)",
        "sources": sources,
    }
    (DATA / "provenance.json").write_text(json.dumps(metadata, indent=2) + "\n")


def read_archive() -> dict:
    metadata = json.loads((DATA / "provenance.json").read_text())
    if metadata.get("baseline_years") != list(BASELINE) or metadata.get("baseline_type") != BASELINE_TYPE:
        raise ValueError("Outdated radiation baseline; rebuild with --refresh-data.")
    path = DATA / metadata["file"]
    if hashlib.sha256(path.read_bytes()).hexdigest() != metadata["sha256"]:
        raise ValueError("Radiation source-data checksum mismatch.")
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != len(MODELS) * len(YEARS) or {r["model"] for r in rows} != set(MODELS):
        raise ValueError("Expected exactly 1,320 model-year rows.")
    data = {}
    for model in MODELS:
        selected = [row for row in rows if row["model"] == model]
        if not np.array_equal([int(r["year"]) for r in selected], YEARS):
            raise ValueError(f"Annual coverage/order mismatch for {model}.")
        data[model] = {q: np.array([float(r[q]) for r in selected]) for q in QUANTITIES}
        data[model]["baseline"] = metadata["baselines"][model]
    return validate_data(data)


def make_figure(data: dict, *, publication: bool = False):
    width, height = (183 / 25.4, 131 / 25.4) if publication else (13, 9.3)
    body, title, panel = (7, 7, 8) if publication else (15, 16, 19)
    member_width, mean_width, reference_width = (0.9, 1, 0.4) if publication else (2, 2.3, 0.85)
    style = {
        "font.family": "Arial", "font.sans-serif": ["Arial"],
        "font.size": body, "text.color": INK, "axes.labelsize": body,
        "axes.labelcolor": INK, "axes.edgecolor": INK,
        "xtick.labelsize": body, "ytick.labelsize": body,
        "xtick.color": INK, "ytick.color": INK, "axes.unicode_minus": True,
        "axes.linewidth": 0.65 if publication else 0.9,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
        "xtick.direction": "out", "ytick.direction": "out",
        "legend.frameon": False, "svg.fonttype": "none", "pdf.fonttype": 42,
        "pdf.use14corefonts": False, "text.usetex": False,
        "figure.facecolor": "white", "savefig.facecolor": "white",
    }
    with plt.rc_context(style):
        fig = plt.figure(figsize=(width, height))
        for index, (quantity, heading, ylabel) in enumerate((
            ("delta_rlut", "Outgoing longwave radiation", "Δrlut (W/m²)"),
            ("delta_sw", "Net downward shortwave radiation", "Δ(rsdt − rsut) (W/m²)"),
        )):
            left = 0.085 if index == 0 else 0.585
            ax = fig.add_axes([left, 0.555, 0.38, 0.29])
            stack = np.stack([smooth(data[m][quantity]) for m in CMIP6_MODELS])
            mean, sd = stack.mean(axis=0), stack.std(axis=0, ddof=1)
            ax.fill_between(YEARS, mean - sd, mean + sd, color=SPREAD, alpha=0.55,
                            linewidth=0, zorder=0)
            for values in stack:
                ax.plot(YEARS, values, color=REFERENCE, alpha=0.5, linewidth=reference_width)
            ax.plot(YEARS, mean, color=INK, linewidth=mean_width, zorder=3)
            for model, linestyle in zip(PLASIC_MODELS, ("-", "--")):
                ax.plot(YEARS, smooth(data[model][quantity]), color=COLORS[model],
                        linestyle=linestyle, linewidth=member_width, zorder=4)
            ax.axhline(0, color="#888888", linestyle="--", linewidth=0.6, zorder=1)
            ax.set(xlim=(1850, 2014), xlabel="Year", ylabel=ylabel)
            ax.set_xticks([1850, 1900, 1950, 2000])
            ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
            ax.margins(y=0.10)
            fig.text(left - 0.045, 0.883, "ab"[index], fontsize=panel, weight="bold", va="top")
            fig.text(left, 0.88, heading, fontsize=title, va="top")

        handles = [
            Line2D([], [], color=REFERENCE, alpha=0.65, linewidth=reference_width),
            (Patch(facecolor=SPREAD, alpha=0.55), Line2D([], [], color=INK, linewidth=mean_width)),
            Line2D([], [], color=COLORS[PLASIC_MODELS[0]], linewidth=member_width),
            Line2D([], [], color=COLORS[PLASIC_MODELS[1]], linewidth=member_width, linestyle="--"),
        ]
        fig.legend(handles, ["CMIP6 models", "CMIP6 mean ± 1 s.d.", *PLASIC_MODELS],
                   handler_map={tuple: HandlerTuple(ndivide=1)},
                   loc="upper center", bbox_to_anchor=(0.52, 0.949), ncol=4,
                   fontsize=7 if publication else 14, handlelength=1.8,
                   columnspacing=1.5, borderaxespad=0)
        fig.text(0.525, 0.985,
                 f"Relative to each run's {BASELINE[0]}–{BASELINE[1]} mean · Centred 9-year time series",
                 ha="center", va="top", fontsize=body)

        # Keep each model visible; form net from paired SW and LW before averaging.
        ax = fig.add_axes([0.245, 0.095, 0.665, 0.32])
        selected = (YEARS >= LAST_DECADE[0]) & (YEARS <= LAST_DECADE[1])
        model_values = np.array([data[m]["delta_net"][selected].mean() for m in CMIP6_MODELS])
        labels = [*CMIP6_MODELS, "CMIP6 mean", *PLASIC_MODELS]
        values = np.array([*model_values, model_values.mean(),
                           *(data[m]["delta_net"][selected].mean() for m in PLASIC_MODELS)])
        sd = model_values.std(ddof=1)
        positions = np.array([0, 1, 2, 3, 4, 5, 6.4, 7.8, 8.8])
        colors = [REFERENCE] * 6 + [INK, *(COLORS[m] for m in PLASIC_MODELS)]
        ax.barh(positions, values, height=0.56, color=colors, edgecolor="none", zorder=2)
        ax.errorbar(values[6], positions[6], xerr=sd, fmt="none", ecolor=INK,
                    elinewidth=0.7 if publication else 1.3,
                    capsize=2 if publication else 4, zorder=3)
        ax.axvline(0, color="#888888", linestyle="--", linewidth=0.6, zorder=0)
        ax.set_yticks(positions, labels)
        ax.set_ylim(positions[-1] + 0.65, -0.65)
        # Numeric annotations sit beyond bar/error extents on either side of zero.
        lower = min(0, float(values.min()), float(values[6] - sd))
        upper = max(0, float(values.max()), float(values[6] + sd))
        span = upper - lower
        for index, (position, value) in enumerate(zip(positions, values)):
            direction = 1 if value >= 0 else -1
            extent = value + direction * (sd if index == 6 else 0)
            ax.text(extent + direction * span * 0.025, position, f"{value:+.2f}".replace("-", "−"),
                    ha="left" if direction > 0 else "right", va="center", fontsize=body)
        ax.set_xlim(lower - span * 0.17, upper + span * 0.20)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
        ax.set_xlabel("Δ(rsdt − rsut) − Δrlut (W/m²)")
        ax.tick_params(axis="y", length=0, pad=5)
        ax.spines["left"].set_visible(False)
        fig.text(0.04, 0.474, "c", fontsize=panel, weight="bold", va="top")
        fig.text(0.085, 0.471, "Net TOA radiation budget, 2005–2014", fontsize=title, va="top")
        fig.text(0.965, 0.471, "Positive = net downward", fontsize=body, va="top", ha="right")
    return fig


def render_radiation_budget(data=None, output_dir: Path = IMAGES) -> None:
    register_arial()
    data = read_archive() if data is None else validate_data(data)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    web = make_figure(data)
    web.savefig(output_dir / f"{STEM}.png", dpi=450)
    plt.close(web)
    print(f"Wrote {STEM}.png (450 dpi, large Arial text) to {output_dir}.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=IMAGES)
    parser.add_argument("--refresh-data", action="store_true", help="Rebuild the small source-data archive from NetCDF.")
    parser.add_argument("--raw-root", type=Path, default=ROOT)
    parser.add_argument("--cmip6-root", type=Path, default=ROOT / "CMIP_data")
    args = parser.parse_args()
    if args.refresh_data or not (DATA / "provenance.json").is_file():
        data, sources = read_raw_data(args.raw_root, args.cmip6_root)
        write_archive(data, sources)
    render_radiation_budget(output_dir=args.output_dir)


if __name__ == "__main__":
    main()
