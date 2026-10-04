#!/usr/bin/env python3
"""Field maps and SEEPS climatology figures for the ERA5-initialized forecasts.

中文说明：为 ERA5 初始化的预报绘制场图与 SEEPS 气候态图。

Only the figures used by the textbook page are produced:

中文说明：只生成教材页面用到的三类图：

* ``state_<case>``   ERA5/PlaSiC maps of T700, u250 and q700 for the direct
  analysis initialization;
  中文说明：``state_<case>``——直接分析初值方案下 T700、u250 和 q700 的
  ERA5/PlaSiC 对比图；
* ``zoom_<case>``    regional evolution for a case with a focus region
  (nudged initialization);
  中文说明：``zoom_<case>``——带关注区域的个例（松弛同化初值）的区域演变图；
* ``seeps_climatology``  SEEPS dry fraction and 24 h wet thresholds.
  中文说明：``seeps_climatology``——SEEPS 干比例与 24 小时湿阈值图。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_ic  # noqa: E402
import common  # noqa: E402
import frame_io  # noqa: E402

# ---------------------------------------------------------------------------
# Publication style
# 出版风格
# ---------------------------------------------------------------------------
# 中文说明：统一的论文图件风格（无衬线字体、细轴线、无上/右边框、无图例边框）。
# English: unified publication style (sans-serif font, thin axes, no top/right
# spines, frameless legend).
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 8.0,
    "axes.linewidth": 0.6,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "legend.frameon": False,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.0,
    "ytick.major.size": 2.0,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

LAT = build_ic.MODEL_LAT                   # 模式纬度 / model latitudes
LON = build_ic.MODEL_LON                   # 模式经度 / model longitudes
DATA_CRS = ccrs.PlateCarree()              # 数据坐标参考系 / data CRS
MAP_CRS = ccrs.Robinson(central_longitude=0)  # 绘图投影（Robinson）/ map projection
# 把 [0, 360) 经度平移到 [-180, 180)，使地图接缝落在太平洋以外。
# Shift [0, 360) longitudes to [-180, 180) so the map seam is placed sensibly.
LON_SHIFTED = ((np.roll(LON, LON.size // 2) + 180.0) % 360.0) - 180.0
LON_EDGES = np.linspace(-180.0, 180.0, LON.size + 1)  # pcolormesh 经度边界 / pcolormesh longitude edges
GLOBAL_REGION = (-180.0, 180.0, -90.0, 90.0)          # 全球区域 / global region

VERIF_ROOT = common.DATA_ROOT / "verification"  # 检验归档目录 / verification archive directory
FIELD_LABELS = {
    "t": "Temperature",
    "u": "Zonal wind",
    "q": "Specific humidity",
}  # 中文说明：状态图的英文标题标签 / English title labels for the state maps.


def load(case_key: str, method: str) -> dict:
    """Load one case/method verification archive.

    中文说明：载入某个例/初值方法的检验归档文件。
    """
    path = VERIF_ROOT / f"{case_key}_{method}.npz"
    if not path.exists():
        raise FileNotFoundError(path)
    return dict(np.load(path, allow_pickle=True))


def land_mask() -> np.ndarray:
    """Model land mask from the winter direct frame stream.

    中文说明：从冬季个例的 direct 帧流中取得模式陆地掩膜。
    """
    path = common.FRAME_ROOT / "winter2020" / "direct_seg1.frames"
    grouped = frame_io.load_stream(path)
    return grouped["land_mask"][0].values[0]


def add_coastlines(ax, land, color="0.25", linewidth=0.3, alpha=0.9):
    """Draw the model land mask as a coastline contour.

    中文说明：把模式陆地掩膜的 0.5 等值线画成海岸线。
    """
    shifted = np.roll(land, LON.size // 2, axis=-1)
    ax.contour(LON_SHIFTED, LAT, shifted, levels=[0.5], colors=color,
               linewidths=linewidth, alpha=alpha, transform=DATA_CRS)


def plot_global(ax, field, cmap, vmin, vmax, land=None):
    """Draw a model-grid field on a Robinson map with a seam-safe longitude roll.

    中文说明：在 Robinson 投影地图上绘制模式网格场，并通过经度滚动保证接缝连续。
    """
    shifted = np.roll(np.asarray(field), LON.size // 2, axis=-1)
    # 首列再补一次，形成闭合的 pcolormesh 网格 / append the first column to close the mesh
    mesh = ax.pcolormesh(LON_EDGES, LAT, np.concatenate([shifted, shifted[:, :1]], axis=-1),
                         cmap=cmap, vmin=vmin, vmax=vmax, shading="auto",
                         rasterized=True, transform=DATA_CRS)
    ax.set_global()
    if land is not None:
        add_coastlines(ax, land)
    return mesh


def field_color_scale(name: str) -> tuple[float, float, str, str]:
    """Return (vmin, vmax, cmap, unit) for a state-map field.

    中文说明：返回状态图某个变量对应的 (最小值, 最大值, 色标, 单位)。
    """
    if name == "t":
        return 210.0, 285.0, "RdYlBu_r", "K"
    if name == "u":
        return -45.0, 45.0, "RdBu_r", "m s⁻¹"
    if name == "q":
        return 0.0, 8.0, "YlGnBu", "g kg⁻¹"
    raise KeyError(name)


# ---------------------------------------------------------------------------
# Figure 1: state evolution (ERA5 vs PlaSiC)
# 图 1：状态演变（ERA5 对 PlaSiC）
# ---------------------------------------------------------------------------
def figure_state(case: common.ForecastCase, method: str = "direct") -> plt.Figure:
    """Compare ERA5 and PlaSiC fields in a large, web-readable map plate.

    中文说明：用大幅面、适合网页阅读的图版对比 ERA5 与 PlaSiC 场。

    Each variable occupies one column.  For each saved lead time, ERA5 is the
    upper map and PlaSiC is the lower map, so a reader can compare one pair
    without scanning across a wide sequence of tiny panels.  The 700-hPa
    temperature map is used in place of the former T500 map.

    中文说明：每个变量占一列。对每个保存的时效，上图为 ERA5、下图为 PlaSiC，
    读者无需扫视一长串小面板就能逐对比较。温度图使用 700 hPa 取代早先的
    T500。
    """
    data = load(case.key, method)
    land = land_mask()
    lead_hours = data["map_times_h"].astype(int)

    variables = [("t", 700), ("u", 250), ("q", 700)]
    nrow, ncol = 2 * len(lead_hours), len(variables)
    # Three wide map columns and fourteen source/lead rows.  At the report's
    # single-column width each map remains roughly twice as wide as in the
    # former seven-column layout.
    # 中文说明：三列宽地图、十四个“数据源/时效”行。在报告的单栏宽度下，
    # 每幅图仍约为早先七列布局时的两倍宽。
    fig = plt.figure(figsize=(11.0, 17.0))
    grid = fig.add_gridspec(nrow, ncol, hspace=0.035, wspace=0.045,
                            left=0.085, right=0.925, top=0.94, bottom=0.085)

    for vi, (name, level) in enumerate(variables):
        key = f"{name}{level}"
        vmin, vmax, cmap, unit = field_color_scale(name)
        model = data[f"map_model_{key}"].astype(float)
        truth = data[f"map_truth_{key}"].astype(float)
        if name == "q":
            # 比湿 kg/kg → g/kg / specific humidity kg/kg -> g/kg
            model = model * 1000.0
            truth = truth * 1000.0
        for lead_index, lead_hour in enumerate(lead_hours):
            # 每个时效两行：奇数行为 ERA5，偶数行为 PlaSiC。
            # Two rows per lead: ERA5 on top, PlaSiC below.
            for source, field, row_offset in (("ERA5", truth, 0), ("PlaSiC", model, 1)):
                row = 2 * lead_index + row_offset
                ax = fig.add_subplot(grid[row, vi], projection=MAP_CRS)
                plot_global(ax, field[lead_index], cmap, vmin, vmax, land)
                ax.xaxis.set_visible(False)
                ax.yaxis.set_visible(False)
                ax.tick_params(pad=0.8, length=1.2)
                if lead_index == 0 and row_offset == 0:
                    ax.set_title(f"{FIELD_LABELS[name]} {level} hPa",
                                 fontsize=7.2, pad=3.0, fontweight="bold")
                if vi == 0:
                    # 左侧竖排标注数据源与时效 / vertical labels for source and lead
                    ax.text(-0.055, 0.5, source, transform=ax.transAxes,
                            fontsize=6.2, color="0.2", ha="right", va="center",
                            rotation=90, rotation_mode="anchor")
                    if row_offset == 0:
                        label = "Analysis" if lead_hour == 0 else f"+{lead_hour} h"
                        ax.text(-0.14, 0.5, label, transform=ax.transAxes,
                                fontsize=6.2, color="0.25", ha="right", va="center",
                                rotation=90, rotation_mode="anchor")
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        # A horizontal bar under each variable column avoids shrinking any
        # map to make room for a vertical colorbar.
        # 中文说明：每列下方放置水平色标，避免为竖直色标而压缩地图。
        cax = fig.add_axes([0.105 + vi * 0.285, 0.035, 0.225, 0.012])
        cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                          orientation="horizontal")
        cb.set_label(f"[{unit}]", fontsize=7.0, labelpad=2.0)
        cb.ax.tick_params(labelsize=8.5, length=1.5, pad=1.0)
    fig.suptitle(f"{case.key}: ERA5 analysis (top) and PlaSiC forecast (bottom)",
                 fontsize=9.0, y=0.968)
    return fig


# ---------------------------------------------------------------------------
# Figure 2: regional zoom (typhoon case)
# 图 2：区域放大（台风个例）
# ---------------------------------------------------------------------------
# 中文说明：ZOOM_FIELDS 每项为
# (键名, 图例标签, 单位, 色标下限, 色标上限, 色标, 数值缩放系数)。
# English: each item is (key, label, unit, vmin, vmax, cmap, scale).
# (key, label, unit, vmin, vmax, cmap, scale)
ZOOM_FIELDS = [
    ("mslp", "Mean sea-level pressure", "hPa", 960.0, 1035.0, "RdBu_r", 0.01),
    ("wind10", "10 m wind speed", "m s⁻¹", 0.0, 32.0, "YlOrRd", 1.0),
    ("vort850", "850 hPa vorticity", "10⁻⁵ s⁻¹", -10.0, 45.0, "PuOr", 1.0e5),
    ("z500", "500 hPa height", "m", 5300.0, 6000.0, "viridis", 1.0),
    ("t500", "500 hPa temperature", "K", 240.0, 280.0, "Spectral_r", 1.0),
    ("wind250", "250 hPa wind speed", "m s⁻¹", 0.0, 80.0, "magma", 1.0),
    ("q700", "700 hPa specific humidity", "g kg⁻¹", 0.0, 14.0, "YlGnBu", 1000.0),
    ("tcc", "Total cloud cover", "%", 0.0, 100.0, "Greys", 100.0),
    ("t2m", "2 m temperature", "K", 265.0, 310.0, "RdYlBu_r", 1.0),
]


def zoom_norm(key: str, vmin: float, vmax: float) -> mpl.colors.Normalize:
    """Colour normalisation for one zoom field.

    中文说明：某个放大区域变量的颜色归一化方式。

    Relative vorticity is dominated by a two-cell-wide tropical-cyclone core
    over a weak large-scale background, so a power norm keeps the background
    structure visible without washing out the storm.

    中文说明：相对涡度由热带气旋核心（约两个格点宽）主导，背景大尺度信号很弱；
    幂律归一化既能保持背景结构可见，又不会让风暴核心过曝。
    """
    if key == "vort850":
        return mpl.colors.PowerNorm(gamma=0.6, vmin=vmin, vmax=vmax)
    return mpl.colors.Normalize(vmin=vmin, vmax=vmax)


def relative_vorticity(u: np.ndarray, v: np.ndarray,
                       lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """Relative vorticity (s^-1) on the model lon / Gaussian-lat grid.

    中文说明：在模式经度/高斯纬度网格上计算相对涡度（s⁻¹）。

    Points whose five-point stencil touches a missing (below-ground) wind
    value are masked so the plateau edge does not leak into the map.

    中文说明：五点差分模板触及缺测（地下）风值的格点会被掩膜，避免高原边缘的
    伪信号泄漏到图中。
    """
    radius = 6371220.0
    lat_rad = np.deg2rad(np.asarray(lat, dtype=float))
    lon_rad = np.deg2rad(np.asarray(lon, dtype=float))
    cos_lat = np.cos(lat_rad)[:, None]
    valid = np.isfinite(u) & np.isfinite(v)
    u_filled = np.where(valid, u, 0.0)
    v_filled = np.where(valid, v, 0.0)
    # 球面坐标下的相对涡度：ζ = (∂v/∂λ − ∂(u·cosφ)/∂φ) / (a·cosφ)
    # Relative vorticity in spherical coordinates.
    dv_dlon = np.gradient(v_filled, lon_rad, axis=1)
    du_dmu = np.gradient(u_filled * cos_lat, lat_rad, axis=0)
    vorticity = (dv_dlon - du_dmu) / (radius * cos_lat)
    # 逐方向滚动掩膜，剔除模板触及缺测值的点 / erode the mask around missing points
    stencil = valid.copy()
    for axis, step in ((0, 1), (0, -1), (1, 1), (1, -1)):
        stencil &= np.roll(valid, step, axis=axis)
    vorticity[~stencil] = np.nan
    return vorticity


def zoom_series(data: dict, key: str) -> tuple[np.ndarray, np.ndarray]:
    """(ERA5, PlaSiC) stacks for one zoom field, before region slicing.

    中文说明：返回某个放大变量在区域切片之前的 (ERA5, PlaSiC) 堆叠数组。
    """
    if key == "vort850":
        # 涡度由 u850/v850 直接诊断 / vorticity is diagnosed from u850/v850
        truth = np.stack([relative_vorticity(np.asarray(u), np.asarray(v), LAT, LON)
                          for u, v in zip(data["map_truth_u850"], data["map_truth_v850"])])
        model = np.stack([relative_vorticity(np.asarray(u), np.asarray(v), LAT, LON)
                          for u, v in zip(data["map_model_u850"], data["map_model_v850"])])
        return truth, model
    if key == "wind250":
        # 250 hPa 风速由 u250、v250 合成 / wind speed from u250 and v250
        truth = np.hypot(data["map_truth_u250"], data["map_truth_v250"])
        model = np.hypot(data["map_model_u250"], data["map_model_v250"])
        return truth, model
    return data[f"map_truth_{key}"], data[f"map_model_{key}"]


def figure_zoom(case: common.ForecastCase, method: str = "nudged") -> plt.Figure:
    """ERA5 and PlaSiC regional zoom over the western North Pacific.

    中文说明：西北太平洋区域的 ERA5 与 PlaSiC 放大对比图。

    Nine diagnostics are shown for the typhoon case; each variable occupies
    two rows so that every ERA5 map sits directly above its PlaSiC partner.

    中文说明：台风个例展示九个诊断量；每个变量占两行，使每一幅 ERA5 图都
    正好位于其 PlaSiC 对应图的上方。
    """
    data = load(case.key, method)
    land = land_mask()
    lead_hours = data["map_times_h"].astype(int)
    west, east, south, north = case.focus_region
    # 按关注区域裁剪经纬度 / slice longitudes and latitudes to the focus region
    lon_mask = (LON >= west) & (LON <= east)
    lat_mask = (LAT >= south) & (LAT <= north)
    lon_zoom = LON[lon_mask]
    lat_zoom = LAT[lat_mask]

    nvar = len(ZOOM_FIELDS)
    nrow, ncol = 2 * nvar, len(lead_hours)
    figure_width, figure_height = 8.0, 13.3
    fig = plt.figure(figsize=(figure_width, figure_height))
    left, right, top, bottom = 0.14, 0.885, 0.955, 0.05
    grid = fig.add_gridspec(nrow, ncol, hspace=0.06, wspace=0.06,
                            left=left, right=right, top=top, bottom=bottom)
    for fi, (key, label, unit, vmin, vmax, cmap, scale) in enumerate(ZOOM_FIELDS):
        truth_series, model_series = zoom_series(data, key)
        norm = zoom_norm(key, vmin, vmax)
        for col in range(ncol):
            # 上行为 ERA5，下行为 PlaSiC / ERA5 on top, PlaSiC below
            for row, field in enumerate([truth_series, model_series]):
                values = np.asarray(field[col], dtype=float)[lat_mask][:, lon_mask] * scale
                ax = fig.add_subplot(grid[2 * fi + row, col], projection=DATA_CRS)
                ax.pcolormesh(lon_zoom, lat_zoom, values, cmap=cmap, norm=norm,
                              shading="auto", rasterized=True, transform=DATA_CRS)
                add_coastlines(ax, land)
                ax.set_xlim(west, east)
                ax.set_ylim(south, north)
                ax.set_xticks([110, 130, 150, 170])
                ax.set_yticks([10, 25, 40, 55])
                if 2 * fi + row == nrow - 1:
                    # 仅最底行显示横轴刻度标签 / only the bottom row shows x labels
                    ax.set_xticklabels(["110°E", "130°E", "150°E", "170°E"],
                                       fontsize=7.0, rotation=45, ha="right",
                                       rotation_mode="anchor")
                else:
                    ax.set_xticklabels([])
                if col == 0:
                    ax.set_yticklabels(["10°N", "25°N", "40°N", "55°N"], fontsize=7.0)
                else:
                    ax.set_yticklabels([])
                ax.tick_params(pad=0.8, length=1.4)
                if row == 0 and fi == 0:
                    # 第一行变量标注各列时效 / the first variable row labels each lead
                    head = "Analysis" if lead_hours[col] == 0 else f"+{lead_hours[col]} h"
                    ax.set_title(head, fontsize=7.0, pad=1.5)
                if key == "tcc" and row == 1 and col == 0:
                    # 模式没有初始云量场，给出文字说明 / model has no initial cloud field
                    ax.text(0.5, 0.5, "no initial cloud field", transform=ax.transAxes,
                            fontsize=7.0, color="0.45", ha="center", va="center")
        # 每行左侧竖排标注数据源 ERA5/PlaSiC / vertical source labels per row
        for row, tag in enumerate(["ERA5", "PlaSiC"]):
            fig.text(0.085, top - (2 * fi + row + 0.5) / nrow * (top - bottom),
                     tag, fontsize=7.0, ha="right", va="center", rotation=90,
                     rotation_mode="anchor", color="0.25")
        # 每个变量的分组标签 / group label for each variable block
        fig.text(0.028, top - (2 * fi + 1.0) / nrow * (top - bottom),
                 label, fontsize=8.0, ha="center", va="center", rotation=90,
                 rotation_mode="anchor", fontweight="bold")
        block_height = (top - bottom) / nvar
        cax = fig.add_axes([0.895, top - (fi + 0.8) * block_height, 0.009, 0.6 * block_height])
        cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
        cb.set_label(f"[{unit}]", fontsize=7.0, labelpad=1.0)
        cb.ax.tick_params(labelsize=8.5, length=1.2, pad=0.8)
    fig.suptitle(
        f"{case.key}: western North Pacific evolution "
        f"({west:.0f}–{east:.0f}°E, {south:.0f}–{north:.0f}°N)",
        fontsize=8.5, y=0.982,
    )
    return fig


# ---------------------------------------------------------------------------
# Figure 3: SEEPS climatological weights (dry fraction and wet threshold)
# 图 3：SEEPS 气候权重（干比例与湿阈值）
# ---------------------------------------------------------------------------
SEEPS_TRUTH_ROOT = common.DATA_ROOT / "wb2_truth"  # WB2 真值缓存 / cached WB2 truth
SEEPS_LEAD_HOURS = 24                               # 湿阈值展示时效 / lead shown for the threshold


def figure_seeps() -> plt.Figure:
    """Dry fraction and 24 h wet threshold behind the SEEPS score.

    中文说明：SEEPS 评分背后的干比例与 24 小时湿阈值。

    Panel (a) is the time-mean dry fraction p1 used by the score.  Panels
    (b)-(d) show the 24 h wet threshold of the three cases at +24 h, read from
    the cached WeatherBench 2 products written by ``fetch_wb2_truth.py``.
    Points without a cached climatological threshold are left white.

    中文说明：面板 (a) 为评分使用的时间平均干比例 p1；面板 (b)–(d) 展示三个
    个例在 +24 h 的 24 小时湿阈值，数据来自 ``fetch_wb2_truth.py`` 写出的
    WeatherBench 2 缓存产品。没有气候阈值的格点保持白色。
    """
    for case in common.CASES:
        path = SEEPS_TRUTH_ROOT / f"{case.key}_precip24.npz"
        if not path.exists():
            raise FileNotFoundError(path)
    p1 = np.load(SEEPS_TRUTH_ROOT / "seeps_p1.npz")["p1"].astype(float)
    land = land_mask()

    fig = plt.figure(figsize=(7.4, 4.9))
    grid = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.05,
                            left=0.015, right=0.985, top=0.95, bottom=0.035)

    def strip_ticks(ax) -> None:
        """Hide axis ticks on a map panel.

        中文说明：隐藏地图面板的坐标刻度。
        """
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

    def horizontal_bar(ax, mappable, label: str):
        """Attach a compact horizontal colour bar to a map panel.

        中文说明：为地图面板附上紧凑的水平色标。
        """
        bar = fig.colorbar(mappable, ax=ax, orientation="horizontal",
                           fraction=0.05, pad=0.03, aspect=26)
        bar.ax.tick_params(labelsize=8.5, length=1.5, pad=1.0)
        bar.set_label(label, fontsize=7.5, labelpad=1.0)
        return bar

    # (a) climatological dry fraction
    # (a) 气候干比例
    ax = fig.add_subplot(grid[0, 0], projection=MAP_CRS)
    mesh = plot_global(ax, p1, "BrBG", 0.0, 1.0)
    add_coastlines(ax, land)
    shifted = np.roll(p1, LON.size // 2, axis=-1)
    # 用虚线标出 SEEPS 评分中 p1 的有效区间 [0.1, 0.85]。
    # Dashed lines mark the valid p1 range [0.1, 0.85] used by SEEPS.
    ax.contour(LON_SHIFTED, LAT, shifted, levels=[0.1, 0.85], colors="0.1",
               linewidths=0.45, linestyles="--", transform=DATA_CRS)
    strip_ticks(ax)
    ax.set_title("(a) dry fraction p1 (1990-2019 mean)", fontsize=8.0, pad=2.5)
    horizontal_bar(ax, mesh, "p1 [-]")

    # (b)-(d) wet threshold of each case at +24 h
    # (b)–(d) 各个例 +24 h 的湿阈值
    for index, case in enumerate(common.CASES):
        data = np.load(SEEPS_TRUTH_ROOT / f"{case.key}_precip24.npz")
        leads = np.asarray(data["leads_hours"], dtype=float)
        lead_index = int(np.argmin(np.abs(leads - SEEPS_LEAD_HOURS)))
        threshold = data["precip24_threshold_mm"][lead_index].astype(float)
        ax = fig.add_subplot(grid[divmod(1 + index, 2)],
                             projection=MAP_CRS)
        mesh = plot_global(ax, threshold, "YlGnBu", None, None)
        # 阈值跨多个量级，使用对数色标 / thresholds span orders of magnitude: log colour scale
        mesh.set_norm(mpl.colors.LogNorm(vmin=0.25, vmax=35.0))
        add_coastlines(ax, land)
        strip_ticks(ax)
        letter = "bcd"[index]
        ax.set_title(f"({letter}) {case.key}: wet threshold at +{SEEPS_LEAD_HOURS} h",
                     fontsize=8.0, pad=2.5)
        horizontal_bar(ax, mesh, "24 h threshold [mm]")

    return fig


# ---------------------------------------------------------------------------
# Saving
# 保存
# ---------------------------------------------------------------------------
def save(fig: plt.Figure, name: str) -> None:
    """Save editable vector files and a 600-dpi raster preview.

    中文说明：把图件以 400 dpi 的 PNG 保存并关闭。
    """
    common.FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    base = common.FIGURE_ROOT / name
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"[figure] {name}")


def main_for_cases(cases: list[common.ForecastCase]) -> None:
    """Generate the per-case maps used by the report without parsing argv.

    中文说明：不解析命令行参数，直接生成报告所需的逐个例图件。

    Also used by the experiment driver (``run_workflow.py``).
    中文说明：试验驱动脚本（``run_workflow.py``）也会调用本函数。
    """
    for case in cases:
        try:
            save(figure_state(case, "direct"), f"state_{case.key}")
            if case.focus_region != GLOBAL_REGION:
                # 只有设置了关注区域的个例才画放大图 / only cases with a focus region get a zoom figure
                save(figure_zoom(case, "nudged"), f"zoom_{case.key}")
        except FileNotFoundError as exc:
            print(f"[skip] {exc}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="*", default=None)
    args = parser.parse_args()
    cases = [case for case in common.CASES if not args.cases or case.key in args.cases]
    main_for_cases(cases)
    try:
        save(figure_seeps(), "seeps_climatology")
    except FileNotFoundError as exc:
        print(f"[skip] {exc}")


if __name__ == "__main__":
    main()
