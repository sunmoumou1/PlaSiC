#!/usr/bin/env python3
"""Rebuild Figure 5 from historical surface-temperature changes.

Usage: python historical_3d_500y/plots/plot_historical_warming_pattern.py

Figure contract
---------------
Claim: the two-member historical simulation has spatially uneven warming,
with stronger northern high-latitude warming.
Evidence: (a) the ensemble surface-temperature change on a Robinson map;
(b) the zonal means of both members and their arithmetic ensemble mean.
Archetype: quantitative comparison, with the map as the dominant panel.
The change is the equal-month 2005-2014 mean minus each member's equal-month
1855-1859 control mean. All 32 x 64 cells and both members are retained.
Individual members are shown instead of an inferred confidence interval.
Python/Cartopy draws the native grid without smoothing. Coastlines use the
bundled Natural Earth coordinates, so rendering requires no downloads.
Exports: PNG only, with large Arial text at 450 dpi.
Inputs: archived fields when available, otherwise local monthly NetCDF runs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import FuncFormatter, MultipleLocator


ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "figures" / "data" / "historical_warming_pattern"
IMAGES = ROOT / "figures"
COASTLINES = ROOT.parent / "plasic_app/src/plasic_app/resources/coastlines_110m.json"
STEM = "fig5_warming_pattern"
BRANCHES = ("A", "B")
CONTROL_PERIOD = (1855, 1859)
HISTORICAL_PERIOD = (2005, 2014)
INK = "#222222"


def register_arial() -> None:
    """Require actual Arial regular and bold, including on macOS."""
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


def validate_fields(lat, lon, deltas) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lat, lon, deltas = (np.asarray(v, dtype=np.float64) for v in (lat, lon, deltas))
    if lat.ndim != 1 or lon.ndim != 1 or deltas.shape != (2, lat.size, lon.size):
        raise ValueError("Expected two complete latitude x longitude fields.")
    if not all(np.isfinite(v).all() for v in (lat, lon, deltas)):
        raise ValueError("Coordinates and both members must contain only finite values.")
    if not (np.all(np.diff(lat) > 0) or np.all(np.diff(lat) < 0)):
        raise ValueError("Latitudes must be strictly monotone.")
    if not np.all(np.abs(lat) < 90):
        raise ValueError("Expected Gaussian latitude centres strictly between the poles.")
    if lon.size < 2 or lon.size % 2 or not np.allclose(np.diff(lon), 360 / lon.size):
        raise ValueError("Expected one complete, uniformly spaced longitude cycle.")
    if not np.isclose(lon[0], 0):
        raise ValueError("Expected native longitudes starting at 0 degrees.")
    return lat, lon, deltas


def fields_from_runs(runs: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Keep the legacy period means, checking grid and monthly coverage."""
    lat = np.asarray(runs["A"]["historical"]["tas"]["lat"])
    lon = np.asarray(runs["A"]["historical"]["tas"]["lon"])
    deltas = []
    for branch in BRANCHES:
        means = []
        for run_name, (start, end) in (
            ("historical", HISTORICAL_PERIOD), ("control", CONTROL_PERIOD)
        ):
            item = runs[branch][run_name]["tas"]
            if not (np.array_equal(item["lat"], lat) and np.array_equal(item["lon"], lon)):
                raise ValueError("Both members and controls must share the same grid.")
            month_ids = np.rint(np.asarray(item["years"]) * 12).astype(int)
            selected = (month_ids >= start * 12) & (month_ids < (end + 1) * 12)
            if not np.array_equal(month_ids[selected], np.arange(start * 12, (end + 1) * 12)):
                raise ValueError(f"Incomplete monthly coverage: {branch} {run_name}.")
            values = np.ma.asarray(item["data"])[selected].astype(np.float64)
            if np.ma.getmaskarray(values).any() or not np.isfinite(values).all():
                raise ValueError("Missing temperatures are not silently excluded.")
            means.append(np.asarray(values).mean(axis=0))
        deltas.append(means[0] - means[1])
    return validate_fields(lat, lon, deltas)


def read_raw_fields(raw_root: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read only TAS when explicitly rebuilding from the original NetCDF runs."""
    import netCDF4 as nc4

    runs = {}
    for branch in BRANCHES:
        runs[branch] = {}
        for run_name in ("historical", "control"):
            with nc4.Dataset(raw_root / f"{run_name}_{branch}" / "monthly/tas.nc") as ds:
                time = ds["time"]
                dates = nc4.num2date(time[:], time.units, time.calendar)
                runs[branch][run_name] = {"tas": {
                    "lat": np.asarray(ds["lat"][:]),
                    "lon": np.asarray(ds["lon"][:]),
                    "years": np.array([d.year + (d.month - 1) / 12 for d in dates]),
                    "data": ds["tas"][:],
                }}
    return fields_from_runs(runs)


def read_archive() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    metadata = json.loads((DATA / "provenance.json").read_text())
    path = DATA / metadata["archive"]["file"]
    if hashlib.sha256(path.read_bytes()).hexdigest() != metadata["archive"]["sha256"]:
        raise ValueError("The archived warming fields have a checksum mismatch.")
    with np.load(path, allow_pickle=False) as archive:
        return validate_fields(archive["lat"], archive["lon"], archive["deltas_K"])


def map_cells(lat, lon, field) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reorder cells with longitude, close the seam and cover both polar caps.

    Longitude edges begin at -180 degrees. Splitting the native cell centred
    on the dateline into two half-cells preserves the value and its total area.
    Latitude edges bisect the Gaussian quadrature weights in sin(latitude),
    so the displayed cell areas cover the sphere without adding polar data.
    """
    longitude = (lon + 180) % 360 - 180
    order = np.argsort(longitude)
    longitude, field = longitude[order], field[:, order]
    spacing = 360 / lon.size
    lon_edges = np.r_[-180, longitude + spacing / 2, 180]
    field = np.column_stack([field, field[:, 0]])
    lat_order = np.argsort(lat)
    lat = lat[lat_order]
    field = field[lat_order]
    nodes, weights = np.polynomial.legendre.leggauss(lat.size)
    if not np.allclose(np.sin(np.deg2rad(lat)), nodes, atol=2e-7, rtol=0):
        raise ValueError("Expected the native Gaussian latitude grid.")
    sin_edges = np.r_[-1, -1 + np.cumsum(weights)]
    lat_edges = np.rad2deg(np.arcsin(np.clip(sin_edges, -1, 1)))
    lat_edges[[0, -1]] = (-90, 90)
    return lon_edges, lat_edges, field


def coastline_segments() -> list[np.ndarray]:
    """Keep coastlines in their native -180..180 degree coordinate system."""
    payload = json.loads(COASTLINES.read_text(encoding="utf-8"))
    segments = []
    for line in payload["lines"]:
        points = np.asarray(line, dtype=float)
        if len(points) < 2:
            continue
        # Break at dateline jumps, never discard western-hemisphere coordinates.
        breaks = np.flatnonzero(np.abs(np.diff(points[:, 0])) > 180) + 1
        segments.extend(part for part in np.split(points, breaks) if len(part) >= 2)
    return segments


def make_figure(lat, lon, deltas, *, publication: bool = False):
    width, height = (183 / 25.4, 86 / 25.4) if publication else (13, 6.1)
    body, title, panel = (7, 7, 8) if publication else (15, 17, 19)
    mean_width, member_width = (1.2, 0.65) if publication else (2.6, 1.35)
    style = {
        "font.family": "Arial", "font.sans-serif": ["Arial"],
        "font.size": body, "text.color": INK, "axes.labelsize": body,
        "axes.labelcolor": INK, "axes.edgecolor": INK,
        "xtick.labelsize": body, "ytick.labelsize": body,
        "xtick.color": INK, "ytick.color": INK,
        "axes.linewidth": 0.7, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": False,
        "xtick.direction": "out", "ytick.direction": "out",
        "legend.frameon": False, "legend.fontsize": body,
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "pdf.use14corefonts": False, "text.usetex": False,
        "figure.facecolor": "white", "savefig.facecolor": "white",
    }
    with plt.rc_context(style):
        fig = plt.figure(figsize=(width, height))
        geographic = ccrs.PlateCarree()
        ax_map = fig.add_axes([0.055, 0.23, 0.61, 0.64],
                              projection=ccrs.Robinson(central_longitude=0))
        ax_map.set_global()
        field = deltas.mean(axis=0)
        lon_edges, lat_edges, cells = map_cells(lat, lon, field)
        # Zero is near-white, negative changes are blue, warming is yellow-red.
        cmap = LinearSegmentedColormap.from_list("warming", [
            (0, "#5A90B6"), (0.5 / 3.5, "#FFF9EB"),
            (1 / 3.5, "#FDE0A2"), (1.75 / 3.5, "#F6AD56"),
            (2.5 / 3.5, "#E36A39"), (1, "#9D2638"),
        ])
        cmap.set_over("#701A30")
        mesh = ax_map.pcolormesh(
            lon_edges, lat_edges, cells, transform=geographic,
            shading="flat", cmap=cmap, norm=Normalize(-0.5, 3),
            edgecolors="none", antialiased=False, rasterized=True, zorder=1,
        )
        for points in coastline_segments():
            ax_map.plot(points[:, 0], points[:, 1], transform=geographic,
                        color=INK, linewidth=0.45 if publication else 0.8, zorder=3)
        ax_map.spines["geo"].set_visible(True)
        ax_map.spines["geo"].set_linewidth(0.6 if publication else 0.8)
        labels = ax_map.gridlines(
            crs=geographic, draw_labels={"bottom": "x"},
            x_inline=False, y_inline=False,
            xlocs=[-120, 0, 120], ylocs=[-60, -30, 0, 30, 60],
            rotate_labels=False,
        )
        labels.xlines = labels.ylines = False
        labels.top_labels = labels.right_labels = labels.left_labels = False
        labels.geo_labels = False
        # Gridliner creates these labels lazily when saving, after this style
        # context ends, so the family must be explicit here too.
        labels.xlabel_style = labels.ylabel_style = {
            "size": body, "color": INK, "family": "Arial",
        }
        labels.xpadding = labels.ypadding = 4 if publication else 7
        # Place latitude labels at the curved left boundary; Gridliner can
        # otherwise suppress high-latitude labels or repeat them on the right.
        for latitude in (-60, -30, 0, 30, 60):
            point = ax_map.projection.transform_point(-180, latitude, geographic)
            label = "0°" if latitude == 0 else f"{abs(latitude)}°{'N' if latitude > 0 else 'S'}"
            ax_map.annotate(label, xy=point, xycoords="data",
                            xytext=(-4 if publication else -7, 0),
                            textcoords="offset points", ha="right", va="center",
                            fontsize=body, annotation_clip=False)
        cax = fig.add_axes([0.17, 0.105, 0.38, 0.027])
        extend = "both" if field.min() < -0.5 else "max"
        cbar = fig.colorbar(mesh, cax=cax, orientation="horizontal", extend=extend)
        cbar.set_ticks([-0.5, 0, 1, 2, 3])
        cbar.set_label("Surface air temperature change (K)", labelpad=5)
        cbar.outline.set_linewidth(0.6)
        cbar.ax.tick_params(length=3, width=0.6, pad=3)

        ax_zonal = fig.add_axes([0.765, 0.23, 0.21, 0.64])
        ax_zonal.axvline(0, color="#888888", linewidth=0.7, linestyle="--", zorder=0)
        for branch, values, color, linestyle in zip(
            BRANCHES, deltas, ("#D55E00", "#0072B2"), ("-", "--")
        ):
            ax_zonal.plot(values.mean(axis=-1), lat, color=color,
                          linewidth=member_width, linestyle=linestyle,
                          label=f"Branch {branch}")
        ax_zonal.plot(field.mean(axis=-1), lat, color=INK,
                      linewidth=mean_width, label="Ensemble mean")
        zonal = deltas.mean(axis=-1)
        ax_zonal.set_xlim(min(-0.3, float(zonal.min()) - 0.15),
                          np.ceil(float(zonal.max()) + 0.15))
        ax_zonal.set_ylim(-90, 90)
        ax_zonal.set_yticks([-60, -30, 0, 30, 60])
        ax_zonal.yaxis.set_major_formatter(FuncFormatter(
            lambda value, _: "0°" if value == 0 else f"{abs(value):g}°{'N' if value > 0 else 'S'}"
        ))
        ax_zonal.xaxis.set_major_locator(MultipleLocator(1))
        ax_zonal.set_xlabel("Temperature change (K)", labelpad=6)
        ax_zonal.set_ylabel("Latitude", labelpad=4)
        handles, names = ax_zonal.get_legend_handles_labels()
        ax_zonal.legend([handles[i] for i in (2, 0, 1)], [names[i] for i in (2, 0, 1)],
                        loc="center right", bbox_to_anchor=(1, 0.40),
                        fontsize=6.5 if publication else 13,
                        handlelength=1.7, labelspacing=0.55, borderaxespad=0.1)
        for letter, x, title_x, heading in (
            ("a", 0.022, 0.063, "Surface warming (2005–2014)"),
            ("b", 0.721, 0.765, "Zonal mean"),
        ):
            fig.text(x, 0.965, letter, fontsize=panel, weight="bold", va="top")
            fig.text(title_x, 0.961, heading, fontsize=title, va="top")
        fig.text(0.063, 0.905, "Relative to the 1855–1859 control", fontsize=body, va="top")
    return fig


def render_warming_pattern(fields=None, output_dir: Path = IMAGES) -> None:
    register_arial()
    if fields is None:
        lat, lon, deltas = read_archive() if (DATA / "provenance.json").is_file() else read_raw_fields(ROOT)
    else:
        lat, lon, deltas = validate_fields(*fields)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    web = make_figure(lat, lon, deltas)
    web.savefig(output_dir / f"{STEM}.png", dpi=450)
    plt.close(web)
    print(f"Wrote {STEM}.png (450 dpi, large Arial text) to {output_dir}.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=IMAGES)
    parser.add_argument("--raw-root", type=Path, help="Recompute from the original NetCDF runs.")
    args = parser.parse_args()
    fields = read_raw_fields(args.raw_root) if args.raw_root else None
    render_warming_pattern(fields, args.output_dir)


if __name__ == "__main__":
    main()
