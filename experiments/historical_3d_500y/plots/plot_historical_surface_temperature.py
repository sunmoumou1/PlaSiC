#!/usr/bin/env python3
"""Rebuild the historical-experiment surface-temperature figure offline.

Usage: python historical_3d_500y/plots/plot_historical_surface_temperature.py

Observations: historical_3d_500y/observation_data/
Outputs: historical_3d_500y/figures/

Figure contract
---------------
Claim: compare the evolution of two PlaSiC historical branches with four
observational GMST products, without implying identical temperature measures.
Evidence: the main time series compares warming; the CO2 inset documents the
forcing. Archetype: quantitative comparison with a subordinate forcing panel.
Temperature changes use each series' own 1850-1900 mean. GISTEMP starts in 1880,
so its 1850-1900 reference is estimated with NASA's overlap-period method using
HadCRUT5, Berkeley Earth and NOAA; no pre-1880 GISTEMP values are invented.
Both model members are drawn individually, without a shaded range.
Observations show product best estimates, without uncertainty bands.
Exports: PNG only, with large Arial text at 450 dpi.
Model inputs use the experiment's annual means and monthly CO2 output.
The derived source-data CSV and QA report are exported alongside the figures.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from matplotlib import font_manager
from matplotlib.ticker import FixedLocator, MultipleLocator


ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "observation_data"
IMAGES = ROOT / "figures"
MODEL_DATA = ROOT / "historical_annual_means.csv"
CO2_DATA = ROOT / "historical_A" / "monthly" / "co2.nc"
STEM = "fig3_surface_temperature"
FIRST_YEAR, LAST_YEAR = 1850, 2014
BASELINE_START, BASELINE_END = 1850, 1900
GISTEMP_OVERLAP_START, GISTEMP_OVERLAP_END = 1880, 1899
INK = "#222222"

# Wong's accessible palette; solid model curves and dashed observations.
OBSERVATIONS = {
    "HadCRUT5": {"color": "#0072B2", "linestyle": "--"},
    "GISTEMP": {"color": "#D55E00", "linestyle": "--"},
    "Berkeley Earth": {"color": "#009E73", "linestyle": "--"},
    "NOAAGlobalTemp": {"color": "#CC79A7", "linestyle": "--"},
}
MEMBERS = {
    "A": {"color": "#555555", "linestyle": "-"},
    "B": {"color": "#8A6F40", "linestyle": "-"},
}


def register_arial() -> None:
    """Require real Arial rather than silently substituting another typeface."""
    for name in ("Arial.ttf", "Arial Bold.ttf", "Arial Italic.ttf"):
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


def read_annual_co2(years: np.ndarray) -> np.ndarray:
    """Read the experiment's monthly, spatially uniform CO2 forcing in ppm."""
    with nc4.Dataset(CO2_DATA) as dataset:
        time = dataset.variables["time"]
        dates = nc4.num2date(time[:], units=time.units, calendar=time.calendar)
        monthly_years = np.array([date.year for date in dates])
        co2 = np.asarray(dataset.variables["co2"][:], dtype=np.float64)
    if not np.array_equal(np.unique(monthly_years), years):
        raise ValueError("CO2 forcing must cover the same years as model temperatures.")
    if not np.isfinite(co2).all():
        raise ValueError("CO2 forcing must contain finite monthly values.")
    return np.array([co2[monthly_years == year].mean() for year in years])


def read_observations() -> dict[str, tuple[np.ndarray, np.ndarray]]:
    with (DATA / "hadcrut5_annual.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    hadcrut = np.array([
        (int(row["Time"]), float(row["Anomaly (deg C)"])) for row in rows
    ])

    with (DATA / "gistemp_v4.csv").open(newline="") as stream:
        next(stream)  # Product title; the second line is the CSV header.
        rows = list(csv.DictReader(stream))
    gistemp = np.array([
        (int(row["Year"]),
         np.nan if row["J-D"].strip() == "***" else float(row["J-D"]))
        for row in rows
    ])

    # Berkeley column 2 uses air temperature over sea ice (recommended product).
    # The source uncertainty and alternative sea-ice fields remain in the input.
    berkeley = np.loadtxt(DATA / "berkeley_earth_annual.txt", comments="%",
                          usecols=(0, 1), encoding="latin-1")
    # NOAA column 3 is error variance, not an anomaly or an uncertainty half-width.
    noaa = np.loadtxt(DATA / "noaa_global_temp_v6.1.0.asc", usecols=(0, 1))
    return {
        "HadCRUT5": (hadcrut[:, 0].astype(int), hadcrut[:, 1]),
        "GISTEMP": (gistemp[:, 0].astype(int), gistemp[:, 1]),
        "Berkeley Earth": (berkeley[:, 0].astype(int), berkeley[:, 1]),
        "NOAAGlobalTemp": (noaa[:, 0].astype(int), noaa[:, 1]),
    }


def reference_mean(years: np.ndarray, values: np.ndarray,
                   start: int = BASELINE_START, end: int = BASELINE_END) -> float:
    """Compute a reference mean only when every requested year is available."""
    if years.shape != values.shape or not np.isfinite(values).all():
        raise ValueError("Every plotted year must have one finite annual value.")
    if not np.all(np.diff(years) == 1):
        raise ValueError("Annual years must be unique, ordered and consecutive.")
    baseline = (years >= start) & (years <= end)
    if not np.array_equal(years[baseline], np.arange(start, end + 1)):
        raise ValueError(f"The complete {start}-{end} reference period is required.")
    return float(values[baseline].mean())


def rebase(years: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Subtract the series' own complete 51-year, 1850-1900 mean."""
    return values - reference_mean(years, values)


def rebase_gistemp(years: np.ndarray, values: np.ndarray,
                  references: dict) -> tuple[np.ndarray, dict]:
    """Apply the three-record overlap correction described in NASA's FAQ.

    G(t) - mean(G, 1880:1899) + mean_j[mean(R_j, 1880:1899)
                                           - mean(R_j, 1850:1900)].
    The correction estimates a reference offset, not missing annual values.
    Compute it from the input observations rather than hard-coding a rounded offset.
    """
    corrections = {}
    for name in ("HadCRUT5", "Berkeley Earth", "NOAAGlobalTemp"):
        ref_years, ref_values = references[name]
        corrections[name] = (
            reference_mean(ref_years, ref_values,
                           GISTEMP_OVERLAP_START, GISTEMP_OVERLAP_END)
            - reference_mean(ref_years, ref_values)
        )
    adjustment = float(np.mean(list(corrections.values())))
    overlap_mean = reference_mean(years, values,
                                  GISTEMP_OVERLAP_START, GISTEMP_OVERLAP_END)
    return values - overlap_mean + adjustment, {
        "method": "NASA FAQ overlap-period conversion; estimated 1850-1900 baseline",
        "source": "https://data.giss.nasa.gov/gistemp/faq/",
        "overlap_period": [GISTEMP_OVERLAP_START, GISTEMP_OVERLAP_END],
        "reference_product_offsets_K": corrections,
        "mean_reference_offset_K": adjustment,
        "original_GISTEMP_overlap_mean_K": overlap_mean,
        "total_offset_added_to_original_GISTEMP_K": adjustment - overlap_mean,
        "caveat": "GISTEMP has no observed 1850-1879 values; its baseline is estimated, not a direct 51-year GISTEMP mean.",
    }


def centred_mean(values: np.ndarray, window: int = 9) -> np.ndarray:
    """Equal-year centred mean, using shortened windows at both endpoints."""
    kernel = np.ones(window)
    return np.convolve(values, kernel, mode="same") / np.convolve(
        np.ones_like(values), kernel, mode="same"
    )


def assemble_data(metrics: dict | None = None) -> tuple[dict, dict]:
    """Read local observations and either annual CSV or live legacy model metrics."""
    model = np.genfromtxt(MODEL_DATA, delimiter=",", names=True)
    years = model["year"].astype(int)
    if not np.array_equal(years, np.arange(FIRST_YEAR, LAST_YEAR + 1)):
        raise ValueError("Model annual means must cover all years from 1850 to 2014.")
    series = {}
    for branch in ("A", "B"):
        values = model[f"tas_{branch}_K"]
        if metrics is not None:
            if not np.array_equal(metrics[branch]["tas_years"], years):
                raise ValueError("Legacy model metrics do not cover 1850-2014.")
            # The control offset cancels when rebasing the anomaly series.
            values = np.asarray(metrics[branch]["tas_anomaly"])
        series[f"PlaSiC {branch}"] = (years, rebase(years, values))
    ensemble = (series["PlaSiC A"][1] + series["PlaSiC B"][1]) / 2
    series["PlaSiC ensemble"] = (years, ensemble)

    audit = {}
    selected_observations = {}
    for name, (obs_years, values) in read_observations().items():
        selected = (obs_years >= FIRST_YEAR) & (obs_years <= LAST_YEAR)
        expected_start = 1880 if name == "GISTEMP" else FIRST_YEAR
        if not np.array_equal(obs_years[selected], np.arange(expected_start, LAST_YEAR + 1)):
            raise ValueError(f"Incomplete historical observations: {name}")
        selected_observations[name] = (obs_years[selected], values[selected])
        audit[name] = {
            "source_rows": len(obs_years), "plotted_rows": int(selected.sum()),
            "excluded_rows": int((~selected).sum()),
            "selection": "1850 <= year <= 2014; no within-period rows removed",
        }
    for name, (obs_years, values) in selected_observations.items():
        if name == "GISTEMP":
            anomaly, conversion = rebase_gistemp(obs_years, values, selected_observations)
            audit[name]["baseline_conversion"] = conversion
        else:
            anomaly = rebase(obs_years, values)
        series[name] = (obs_years, anomaly)
    series["CO2"] = (years, read_annual_co2(years))
    return series, audit


def make_figure(series: dict):
    """Draw the surface-temperature comparison for PNG export."""
    width, height = 10.8, 5.9
    body, label, legend, inset = 15, 16, 14.5, 13
    line_width, model_width = 1.35, 2.6

    style = {
        "font.family": "Arial", "font.sans-serif": ["Arial"],
        "font.size": body, "text.color": INK,
        "axes.labelcolor": INK, "axes.edgecolor": INK,
        "xtick.color": INK, "ytick.color": INK,
        "axes.labelsize": label, "xtick.labelsize": body,
        "ytick.labelsize": body, "axes.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False, "axes.unicode_minus": True,
        "xtick.direction": "out", "ytick.direction": "out",
        "legend.frameon": False, "legend.fontsize": legend,
        "text.usetex": False,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white",
    }
    with plt.rc_context(style):
        fig, ax = plt.subplots(figsize=(width, height))
        fig.subplots_adjust(left=0.105, right=0.98, bottom=0.13, top=0.83)
        years, branch_a = series["PlaSiC A"]
        branch_b = series["PlaSiC B"][1]
        member_handles = []
        for branch, values in (("A", branch_a), ("B", branch_b)):
            line, = ax.plot(years, values, **MEMBERS[branch],
                            lw=line_width, label=f"PlaSiC member {branch}", zorder=2)
            member_handles.append(line)
        model_line, = ax.plot(
            years, centred_mean(series["PlaSiC ensemble"][1]), color=INK,
            linestyle="-", lw=model_width,
            label="PlaSiC ensemble (9-year mean)", zorder=4,
        )
        observation_handles = []
        for name, appearance in OBSERVATIONS.items():
            obs_years, values = series[name]
            line, = ax.plot(obs_years, values, **appearance, lw=line_width,
                            label=name, zorder=3)
            observation_handles.append(line)

        ax.axhline(0, color="#666666", lw=0.65, zorder=0)
        ax.set_xlim(FIRST_YEAR, LAST_YEAR + 2)
        ax.set_ylim(-0.5, 1.3)
        ax.xaxis.set_major_locator(FixedLocator([1850, 1880, 1910, 1940, 1970, 2000, 2014]))
        ax.yaxis.set_major_locator(MultipleLocator(0.3))
        ax.tick_params(length=3.5, width=0.8, pad=4)
        ax.set_xlabel("Year", labelpad=8)
        ax.set_ylabel("Temperature change, ΔT (K)", labelpad=8)

        # This inset occupies the empty upper-left region, above all data curves.
        forcing = ax.inset_axes([0.055, 0.615, 0.385, 0.32])
        forcing.plot(series["CO2"][0], series["CO2"][1], color="#6C604F",
                     lw=line_width * 1.15)
        forcing.set_xlim(FIRST_YEAR, LAST_YEAR)
        forcing.set_ylim(275, 410)
        forcing.set_xticks([1850, 1950, 2014])
        forcing.set_yticks([280, 340, 400])
        forcing.tick_params(labelsize=inset, length=2.5, width=0.7, pad=3)
        forcing.set_title("Atmospheric CO2 (ppm)", loc="left", fontsize=inset, pad=6)
        fig.legend(handles=[*member_handles, model_line], loc="upper center",
                   bbox_to_anchor=(0.54, 0.995), ncol=3, handlelength=2.6,
                   columnspacing=2.2, handletextpad=0.6)
        fig.legend(handles=observation_handles, loc="upper center",
                   bbox_to_anchor=(0.54, 0.925), ncol=4, handlelength=2.6,
                   columnspacing=2.0, handletextpad=0.6)
        return fig


def export_source_data(series: dict, audit: dict, output_dir: Path = IMAGES) -> None:
    """Export plotted annual values; retain GISTEMP's absent early years as blanks."""
    names = ["PlaSiC A", "PlaSiC B", "PlaSiC ensemble", *OBSERVATIONS]
    lookup = {name: dict(zip(*series[name])) for name in names}
    smooth = centred_mean(series["PlaSiC ensemble"][1])
    with (output_dir / "fig3_source_data.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["year", *[f"{name} anomaly (K)" for name in names],
                         "PlaSiC ensemble 9-year mean (K)", "CO2 (ppm)"])
        for i, year in enumerate(series["PlaSiC ensemble"][0]):
            writer.writerow([int(year), *[
                f"{lookup[name][year]:.8f}" if year in lookup[name] else ""
                for name in names
            ], f"{smooth[i]:.8f}", f"{series['CO2'][1][i]:.6f}"])
    audit["baseline"] = [BASELINE_START, BASELINE_END]
    audit["smoothing"] = "Model ensemble only: centred 9-year mean, shortened endpoint windows."
    audit["model_members"] = "Both members shown as separate annual curves; no shaded range."
    audit["observations"] = "Annual product best estimates; uncertainty bands are not displayed."
    audit["comparison_caveat"] = "PlaSiC TAS is global surface air temperature; observed GMST blends land air temperature with ocean SST."
    audit["summary_K"] = {}
    for name in names:
        years, values = series[name]
        audit["summary_K"][name] = {
            "2005-2014": float(values[(years >= 2005) & (years <= 2014)].mean()),
            "2014": float(values[years == 2014][0]),
        }
    (output_dir / "qa.json").write_text(json.dumps(audit, indent=2) + "\n")


def render_surface_temperature(metrics: dict | None = None, output_dir: Path = IMAGES) -> None:
    register_arial()
    series, audit = assemble_data(metrics)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    web = make_figure(series)
    web.savefig(output_dir / f"{STEM}.png", dpi=450)
    plt.close(web)
    export_source_data(series, audit, output_dir)
    print(f"Wrote {STEM}.png (450 dpi, large Arial text) to {output_dir}.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=IMAGES)
    args = parser.parse_args()
    render_surface_temperature(output_dir=args.output_dir)


if __name__ == "__main__":
    main()
