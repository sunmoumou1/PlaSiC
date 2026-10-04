#!/usr/bin/env python3
"""Build PlaSiC initial conditions and nudging targets from ERA5.

The grid and atmospheric level count follow :mod:`common`, so the same code
works with T85/L25 by default and with another compiled model configuration.

Pipeline per case:
1. load ERA5 pressure-level and single-level NetCDF files;
2. regrid every field bilinearly from the 1.5 deg regular grid to the model's
   Gaussian grid (north -> south, longitude fastest);
3. correct the ERA5 surface pressure hydrostatically to the model orography;
4. interpolate T, q, u, v from ERA5 pressure levels to the configured sigma
   levels;
5. analyse the grid fields into PlaSiC spectra with the model's own SHTns
   kernels (``tools/plasic_era5_tools``);
6. assemble a restart file that keeps the template's land/ocean static fields
   but replaces the atmosphere and the key surface prognostic fields.

The script also writes a compact binary nudging file holding the hourly target
spectra used by the Newtonian-relaxation initialization driver.
"""
from __future__ import annotations

import argparse
import datetime as dt
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
import restart_io  # noqa: E402

TOOL = common.TOOLS_ROOT / "plasic_era5_tools"                # 谱分析辅助工具 / spectral helper tool
TEMPLATE_RESTART = common.TEMPLATE_RESTART  # 与模式网格/层数匹配的模板重启文件 / matching template restart

# ERA5 ``swvl1..4`` are volumetric soil-water fractions (m³/m³).  These are
# the thicknesses of the four ERA5 soil layers, in metres; multiplying a
# fraction by a layer thickness therefore gives water-equivalent depth in m.
# 中文说明：ERA5 ``swvl1..4`` 是体积含水量（m³/m³），这里的数值是四层
# 土壤层厚度（m）；体积含水量乘以层厚后得到水当量厚度（m）。
ERA5_SOIL_DEPTHS = np.array([0.07, 0.21, 0.72, 1.89])


# ---------------------------------------------------------------------------
# Gaussian grid helpers
# 高斯网格辅助函数
# ---------------------------------------------------------------------------
def gaussian_latitudes(nlat: int) -> np.ndarray:
    """Gaussian latitudes in degrees, north -> south (matches the model).
    """
    nodes, _ = np.polynomial.legendre.leggauss(nlat)
    return np.degrees(np.arcsin(nodes[::-1]))


def gaussian_quadrature_weights(nlat: int) -> np.ndarray:
    """Normalised Gaussian quadrature weights, ordered north to south.
    """
    _, weights = np.polynomial.legendre.leggauss(nlat)
    weights = weights[::-1].astype(np.float64)
    return weights / weights.sum()


MODEL_LAT = gaussian_latitudes(common.NLAT)                 # 模式纬度 / model latitudes
MODEL_WEIGHTS = gaussian_quadrature_weights(common.NLAT)    # 面积权重 / area weights
MODEL_LON = np.arange(common.NLON) * (360.0 / common.NLON)  # 模式经度 / model longitudes
MODEL_COS_LAT = np.cos(np.deg2rad(MODEL_LAT))               # cos(纬度) / cosine of latitude


def model_orography() -> tuple[np.ndarray, float]:
    """Synthesize the model orography grid (m) and reference pressure (Pa).

    The first call invokes the helper tool once; the result is cached under
    ``data/`` for subsequent runs.
    中文说明：首次调用时会运行一次辅助工具，结果缓存到 ``data/`` 目录供后续运行复用。
    """
    # Include the horizontal grid in the cache name.  A cache generated for
    # another T resolution cannot be reshaped safely into this configuration.
    # 中文说明：缓存文件名包含水平网格尺寸；其他 T 分辨率生成的缓存不能直接
    # reshape 到当前配置，避免切换网格后误用旧地形。
    cache = common.DATA_ROOT / (
        f"model_orography_t{common.NTRU}_nlat{common.NLAT}.bin"
    )
    meta_cache = common.DATA_ROOT / (
        f"model_orography_t{common.NTRU}_nlat{common.NLAT}.txt"
    )
    expected_bytes = common.NLAT * common.NLON * np.dtype("<f4").itemsize
    if common.NTRU == 85 and common.NLAT == 128:
        # Reuse the original untagged T85 cache when it has the expected size,
        # preserving existing archives while keeping custom grids isolated.
        legacy_cache = common.DATA_ROOT / "model_orography.bin"
        legacy_meta = common.DATA_ROOT / "model_orography.txt"
        if (legacy_cache.exists() and legacy_meta.exists()
                and legacy_cache.stat().st_size == expected_bytes):
            cache, meta_cache = legacy_cache, legacy_meta
    if (not cache.exists() or not meta_cache.exists()
            or cache.stat().st_size != expected_bytes):
        subprocess.run(
            [str(TOOL), "orography", str(TEMPLATE_RESTART), str(cache), str(meta_cache)],
            check=True,
        )
    orography = np.fromfile(cache, dtype="<f4").reshape(common.NLAT, common.NLON)
    psurf = float(
        [line for line in meta_cache.read_text().splitlines() if line.startswith("psurf")][0].split()[1]
    )

    # the helper writes the dimensionless geopotential scaled by (a*Omega)^2;
    # divide by g to obtain metres
    # 中文说明：辅助工具写出的是按 (a·Ω)² 无量纲化的位势，除以 g 得到米。
    return (orography.astype(np.float64) / common.GRAVITY), psurf


TOP_BLEND_FULL_LEVELS = 3   # levels 0..2 take the model background
                            # 中文说明：第 0–2 层完全采用模式背景
TOP_BLEND_TAPER_LEVELS = 2  # levels 3..4 taper back to the ERA5 analysis
                            # 中文说明：第 3–4 层平滑过渡回 ERA5 分析场


def top_level_weights(nlev: int = common.NLEV) -> np.ndarray:
    """Weight of the model background in the top-of-model blend.

    中文说明：模式顶层混合中“模式背景场”所占的权重。
    """
    weights = np.zeros(nlev)
    for level in range(nlev):
        if level < TOP_BLEND_FULL_LEVELS:
            weights[level] = 1.0
        elif level < TOP_BLEND_FULL_LEVELS + TOP_BLEND_TAPER_LEVELS:
            # 线性锥形过渡：两层内从 1 降到 0 / linear taper from 1 to 0 over two levels
            weights[level] = (
                TOP_BLEND_FULL_LEVELS + TOP_BLEND_TAPER_LEVELS - level
            ) / TOP_BLEND_TAPER_LEVELS
    return weights


def blend_top_levels(spectra: dict, background: dict, weights: np.ndarray) -> dict:
    """Blend the top-of-model spectral coefficients toward the model background.

    中文说明：把模式顶部的谱系数向模式背景场混合。

    The upper sigma levels are not suitable to accept a full analysis increment; inserting the full ERA5 truncation there can excite a grid-scale computational mode. Above the blend depth the model background from the configured template restart is therefore kept, with a short taper back to the ERA5 analysis.
    """
    blended = {}
    for name, values in spectra.items():
        if name == "sp":
            blended[name] = values  # 地面气压谱不参与混合 / surface-pressure spectra are kept as-is
            continue
        field = np.asarray(values, dtype=np.float64).reshape(common.NLEV, common.NRSP)
        background_field = np.asarray(background[name], dtype=np.float64).reshape(
            common.NLEV, common.NRSP
        )
        # 逐层加权混合：field + (background − field)·w / per-level weighted blend
        blended[name] = (
            field + (background_field - field) * weights[:, None]
        ).reshape(-1)
    return blended


def background_spectra() -> dict:
    """Read the configured template restart for the top-level background.
    """
    records = dict(restart_io.read_raw_records(TEMPLATE_RESTART))
    return {
        name: restart_io.decode_float(records[name], common.NRSP * common.NLEV)
        for name in ("st", "sd", "sz", "sq")
    }


# ---------------------------------------------------------------------------
# ERA5 loading / regridding
# ERA5 读取与重网格化
# ---------------------------------------------------------------------------
def open_pl(paths: "list[Path]") -> xr.Dataset:
    """Open and concatenate pressure-level files along the time axis.
    """
    datasets = [xr.open_dataset(path) for path in paths]
    ds = xr.concat(datasets, dim="valid_time") if len(datasets) > 1 else datasets[0]
    # drop duplicate times introduced by the concat boundary, keep first
    # 中文说明：去掉拼接边界引入的重复时刻，保留第一次出现的那个。
    _, index = np.unique(ds["valid_time"].values, return_index=True)
    ds = ds.isel(valid_time=np.sort(index))
    return ds


def open_time_datasets(paths: "list[Path]") -> xr.Dataset:
    """Open exact-time files produced by the daily CDS downloader.

    中文说明：打开按天拆分的精确时刻文件（由 CDS 下载器生成）。
    """
    paths = sorted(paths)
    if not paths:
        raise FileNotFoundError("no ERA5 files found")
    datasets = [xr.open_dataset(path) for path in paths]
    if len(datasets) == 1:
        return datasets[0]
    
    # 兼容 CDS 的 ``valid_time`` 与 ``time`` 两种时间维命名。
    # Support both the ``valid_time`` and ``time`` dimension names from CDS.
    dim = "valid_time" if "valid_time" in datasets[0].dims else "time"
    ds = xr.concat(datasets, dim=dim)
    if dim in ds:
        _, index = np.unique(ds[dim].values, return_index=True)
        ds = ds.isel({dim: np.sort(index)})
    return ds


def to_model_grid(da: xr.DataArray) -> np.ndarray:
    """Bilinear regrid of a (..., latitude, longitude) field to the model grid.
    """
    lon = da["longitude"].values
    if lon.max() > 180.0:
        # 0..360 经度转 -180..180 / convert 0..360 longitudes to -180..180
        da = da.assign_coords(longitude=(lon + 180.0) % 360.0 - 180.0)
    da = da.sortby("longitude")
    lon360 = da["longitude"].values % 360.0
    order = np.argsort(lon360)
    da = da.isel(longitude=order).assign_coords(longitude=lon360[order])
    
    # 在经度 360° 处补一列，保证跨 0° 经线的插值连续。
    # Add a wrap column at 360 deg so interpolation across the 0-deg seam is continuous.
    wrap = da.isel(longitude=0).assign_coords(longitude=360.0)
    da = xr.concat([da, wrap], dim="longitude")
    out = da.interp(longitude=MODEL_LON, latitude=MODEL_LAT)
    return out.transpose(*da.dims).values.astype(np.float64)


def hydrostatic_ps_correction(ps_era5, orog_era5, temperature, orog_model, psurf):
    """Reduce ERA5 surface pressure to the model's (smoother) orography.
    """
    # 局地标高 H = R·T / g（温度取下限 200 K，避免极端值）
    # local scale height H = R*T/g, with temperature floored at 200 K
    scale_height = common.GASCON * np.maximum(temperature, 200.0) / common.GRAVITY
    return ps_era5 * np.exp(-(orog_model - orog_era5) / scale_height)


def sigma_interp(levels, fields, ps, p_target):
    """Linear-in-log-p interpolation of several fields to target pressures.

    ``levels`` are the ERA5 pressure levels (Pa, increasing downward),
    ``fields`` a list of (..., nlev, nlat, nlon) arrays, ``ps`` the surface
    pressure (Pa) and ``p_target`` (nlev_model, nlat, nlon) the target
    pressures.  Below-ground targets take the deepest ERA5 level value.
    """
    out = []
    for field in fields:
        log_target = np.log(p_target)
        
        # vectorised per level-pair search
        # 中文说明：向量化的“所在层对”查找：idx 为下边界层序号。
        idx = np.searchsorted(levels, p_target, side="right") - 1

        idx = np.clip(idx, 0, len(levels) - 2)
        p_lo = levels[idx]
        p_hi = levels[idx + 1]
        
        # 对数气压线性权重 w ∈ [0, 1] / linear-in-log-p weight, clipped to [0, 1]
        w = (log_target - np.log(p_lo)) / (np.log(p_hi) - np.log(p_lo))
        w = np.clip(w, 0.0, 1.0)
        low = np.take_along_axis(field, idx, axis=0)
        high = np.take_along_axis(field, idx + 1, axis=0)
        interp = low * (1.0 - w) + high * w
        out.append(interp)
    return out


@dataclass
class ModelState:
    """One analysis mapped onto the model grid and sigma levels.
    """

    when: dt.datetime
    t: np.ndarray       # (nlev, nlat, nlon) K
    q: np.ndarray       # (nlev, nlat, nlon) kg/kg
    u: np.ndarray       # (nlev, nlat, nlon) m/s
    v: np.ndarray       # (nlev, nlat, nlon) m/s
    ps: np.ndarray      # (nlat, nlon) Pa, corrected to the model orography / 地面气压（Pa，已订正到模式地形）
    surface: dict       # 地表辅助字段字典 / dictionary of auxiliary surface fields


class CaseBuilder:
    """Build all initial-condition products for one forecast case.
    """

    def __init__(self, case: common.ForecastCase):
        self.case = case
        self.t0 = common.init_datetime(case)
        self.directory = common.ERA5_ROOT / case.key
        self.orography, self.psurf = model_orography()
        self._cache: dict = {}  # 已打开的 NetCDF 数据集缓存

    # -- loading ------------------------------------------------------------
    # -- 数据读取 -----------------------------------------------------------
    def _pl_dataset(self, kind: str) -> xr.Dataset:
        """Open the pressure-level archive of ``nudge`` or ``verif``.
        """
        key = f"pl_{kind}"
        if key not in self._cache:
            if kind == "verif":
                paths = sorted(self.directory.glob("verif6_pl*.nc"))
            else:
                paths = sorted(self.directory.glob("nudge_pl_*.nc"))
            self._cache[key] = open_pl(paths)
        return self._cache[key]

    def _sl_dataset(self, kind: str) -> xr.Dataset:
        """Open a single-level archive, e.g. ``surface_ic``.
        """
        key = f"sl_{kind}"
        if key not in self._cache:
            paths = sorted(self.directory.glob(f"{kind}.nc"))
            paths += sorted(self.directory.glob(f"{kind}_*.nc"))
            self._cache[key] = open_time_datasets(paths)
        return self._cache[key]

    # -- state construction -------------------------------------------------
    # -- 状态构建 -----------------------------------------------------------
    def build_state(self, when: dt.datetime) -> ModelState:
        """Map the ERA5 analysis at ``when`` onto model grid and sigma levels.
        """
        # t0 及以前取松弛窗口数据，之后取检验数据。
        # Use the nudging-window archive up to t0, the verification archive after t0.
        ds = self._pl_dataset("nudge" if when <= self.t0 else "verif")
        time_index = int(np.argmin(np.abs(ds["valid_time"].values - np.datetime64(when))))
        selected = ds.isel(valid_time=time_index)

        # 气压层单位由 hPa 转 Pa，并严格按升序排列（向下递增）。
        # Convert pressure levels from hPa to Pa and sort ascending (increasing downward).
        levels = np.asarray(selected["pressure_level"].values, dtype=np.float64) * 100.0
        order = np.argsort(levels)
        levels = levels[order]

        t = to_model_grid(selected["t"])[order]
        q = to_model_grid(selected["q"])[order]
        u = to_model_grid(selected["u"])[order]
        v = to_model_grid(selected["v"])[order]

        sl = self._sl_dataset("surface_ic")
        sl_index = int(np.argmin(np.abs(sl["valid_time"].values - np.datetime64(when))))
        sl_sel = sl.isel(valid_time=sl_index)
        ps_era5 = to_model_grid(sl_sel["sp"]).astype(np.float64)
        orog_era5 = to_model_grid(sl_sel["z"]).astype(np.float64) / common.GRAVITY
        skt = to_model_grid(sl_sel["skt"]).astype(np.float64)
        t2m = to_model_grid(sl_sel["t2m"]).astype(np.float64)
        sst = to_model_grid(sl_sel["sst"]).astype(np.float64)
        ice = to_model_grid(sl_sel["siconc"]).astype(np.float64)
        # ERA5 SST/sea-ice are undefined over land; fall back to skin temperature / zero ice there.
        sst = np.where(np.isfinite(sst), sst, skt)
        ice = np.where(np.isfinite(ice), ice, 0.0)

        ps = hydrostatic_ps_correction(ps_era5, orog_era5, t2m, self.orography, self.psurf)
        # 地面气压限制在 400–1080 hPa 的物理合理范围。
        # Clip surface pressure to a physically plausible 400-1080 hPa range.
        ps = np.clip(ps, 40000.0, 108000.0)

        # 每个 σ 层的目标气压 = σ_full × ps / target pressure of each sigma level
        p_target = common.SIGMA_FULL[:, None, None] * ps[None, :, :]
        t_lev, q_lev, u_lev, v_lev = sigma_interp(
            levels, [t, q, u, v], ps, p_target
        )
        # 比湿上限 45 g/kg，避免插值产生过饱和的极端值。
        # Cap specific humidity at 45 g/kg to avoid supersaturated extremes.
        q_lev = np.clip(q_lev, 0.0, 0.045)

        surface = {
            "skt": skt,
            "t2m": t2m,
            "sst": sst,
            "ice": ice,
            "ps_era5": ps_era5,
            "orog_era5": orog_era5,
            # 四层土壤温度 (stl1..4) 与体积含水量 (swvl1..4, m³/m³)。
            # Four soil-temperature layers (stl1..4) and volumetric soil-water
            # fractions (swvl1..4, m³/m³).
            "soil_temperature": [
                to_model_grid(sl_sel[f"stl{level}"]).astype(np.float64) for level in range(1, 5)
            ],
            "soil_water": [
                to_model_grid(sl_sel[f"swvl{level}"]).astype(np.float64) for level in range(1, 5)
            ],
            "snow_depth": to_model_grid(sl_sel["sd"]).astype(np.float64),
        }
        return ModelState(when=when, t=t_lev, q=q_lev, u=u_lev, v=v_lev, ps=ps, surface=surface)

    # -- spectral analysis --------------------------------------------------
    # -- 谱分析 -------------------------------------------------------------
    def analyse(self, state: ModelState, workdir: Path) -> dict[str, np.ndarray]:
        """Analyse gridded fields into PlaSiC spectral coefficients.
        """
        workdir.mkdir(parents=True, exist_ok=True)
        state.t.astype("<f4").tofile(workdir / "t.bin")
        state.q.astype("<f4").tofile(workdir / "q.bin")
        # the model's winds are Robert-form, nondimensional: U = u cos(lat)/cv
        (state.u * MODEL_COS_LAT[None, :, None] / common.CV).astype("<f4").tofile(workdir / "u.bin")
        (state.v * MODEL_COS_LAT[None, :, None] / common.CV).astype("<f4").tofile(workdir / "v.bin")
        state.ps.astype("<f4").tofile(workdir / "ps.bin")
        outdir = workdir / "spectra"
        outdir.mkdir(exist_ok=True)
        env = __import__("os").environ.copy()
        # 把参考气压通过环境变量传给辅助工具 / pass the reference pressure via the environment
        env["PLASIC_PSURF"] = repr(self.psurf)
        subprocess.run(
            [
                str(TOOL), "analyse", str(common.NLEV),
                str(workdir / "t.bin"), str(workdir / "u.bin"),
                str(workdir / "v.bin"), str(workdir / "q.bin"),
                str(workdir / "ps.bin"), str(outdir),
            ],
            check=True, env=env,
        )
        # Return temperature, divergence, vorticity, humidity and surface-pressure spectra.
        return {
            name: np.fromfile(outdir / f"{name}.bin", dtype="<f4")
            for name in ("st", "sd", "sz", "sq", "sp")
        }

    # -- restart assembly ---------------------------------------------------
    # -- 重启文件组装 -------------------------------------------------------
    def surface_replacement_arrays(self, state: ModelState, template_records) -> dict[str, np.ndarray]:
        """Return physically consistent surface records in model layout.

        Keeping this conversion separate lets the nudging driver relax the
        slowly varying land/ocean state as well as the atmospheric spectra.
        The previous workflow updated only atmosphere records, leaving the
        nudged restart with the surface state from t0−window.
        """
        arrays: dict[str, np.ndarray] = {}

        # --- surface state -------------------------------------------------
        # dls >= 0.5 marks land (same mask convention as the model).
        dls = restart_io.decode_float(dict(template_records)["dls"]).reshape(common.NLAT, common.NLON)
        land = dls >= 0.5

        def land_or_sentinel(field, sentinel):
            """Copy a field but write ``sentinel`` over ocean points.
            """
            out = np.array(field, dtype="<f4")
            out[~land] = sentinel
            return out

        skt = state.surface["skt"]
        sst = state.surface["sst"]
        t2m = state.surface["t2m"]
        ice = np.clip(state.surface["ice"], 0.0, 1.0)
        ice[land] = 0.0  # 陆地上不存在海冰 / no sea ice over land

        # The model's slab ocean cannot be colder than the freezing point, and
        # its sea-ice thermodynamics needs a consistent cover/thickness pair.
        # ERA5 SST is undefined under ice, so fall back to the freezing point
        # there; the ice surface temperature comes from the ERA5 skin
        # temperature.
        freezing_point = 271.25
        ocean_temperature = np.where(np.isfinite(sst), np.maximum(sst, freezing_point), freezing_point)
        ice_mask = ice >= 0.5
        ice_surface_temperature = np.where(ice_mask, np.minimum(skt, freezing_point), ocean_temperature)
        sea_surface_temperature = np.where(ice_mask, ice_surface_temperature, ocean_temperature)

        arrays["dts"] = sea_surface_temperature   # 海表温度 / sea surface temperature
        arrays["ysst"] = ocean_temperature        # 海洋混合层温度 / ocean mixed-layer temperature
        arrays["xts"] = (
            np.where(ice_mask, np.minimum(skt, 272.0), freezing_point)
        )                                         # 海冰表面温度 / sea-ice surface temperature
        arrays["xicec"] = ice                     # 海冰密集度 / sea-ice concentration

        # Ice thickness from the template's own climatology, so that a cover
        # of one is never paired with a zero thickness.
        # 中文说明：海冰厚度取自模板自带的气候态，避免出现覆盖度为 1 但厚度为 0
        # 的不自洽组合。
        climatology = restart_io.decode_float(dict(template_records)["xcliced"]).reshape(
            14, common.NLAT, common.NLON
        )
        month = state.when.month
        thickness = np.where(ice_mask, np.maximum(climatology[month], 0.5), 0.0)
        arrays["xiced"] = thickness  # 海冰厚度 / sea-ice thickness

        # 陆面温度类记录：陆地上使用皮温，海洋上置哨兵值。
        # Land temperature records use skin temperature on land and sentinels over ocean.
        dtsl = land_or_sentinel(skt, 1.0e20)
        arrays["dtsl"] = dtsl   
        arrays["dtsm"] = dtsl 
        arrays["dsnowt"] = land_or_sentinel(skt, 1.0e20)

        soil_layers = []
        for level in range(5):
            # ERA5 只有 4 层土壤温度；模式需要 5 层，最深两层共用第 4 层值。
            # ERA5 has 4 soil layers, the model 5; the deepest two share layer 4.
            source = state.surface["soil_temperature"][min(level, 3)]
            soil_layers.append(land_or_sentinel(source, 1.0e20))
        arrays["dsoilt"] = np.concatenate([x.ravel() for x in soil_layers])

        # ERA5 swvl1..4 are m³/m³.  Integrating each fraction over its layer
        # thickness (m) gives the column water-equivalent depth (m).  This is
        # the unit used by PlaSiC's dwatc/dwmax records, so no factor of 1000
        # is applied here.
        # 中文说明：ERA5 swvl1..4 的单位是 m³/m³，乘以各层厚度（m）后得到
        # 柱状水当量厚度（m）。PlaSiC 的 dwatc/dwmax 也使用 m，因此不乘 1000。
        soil_water = np.zeros((common.NLAT, common.NLON))
        for layer, depth in zip(state.surface["soil_water"], ERA5_SOIL_DEPTHS):
            soil_water += layer * depth
        dwmax = restart_io.decode_float(dict(template_records)["dwmax"]).reshape(common.NLAT, common.NLON)
        soil_water = np.clip(soil_water, 0.0, np.maximum(dwmax, 1.0))
        soil_water[~land] = 0.0
        # 柱状土壤水当量厚度（m）/ column soil-water equivalent depth (m)
        arrays["dwatc"] = soil_water

        snow = np.clip(state.surface["snow_depth"], 0.0, 5.0)
        snow[~land] = 0.0
        arrays["dsnowz"] = snow  # 积雪深度（水当量）/ snow depth (water equivalent)

        # atmosphere-facing surface temperature/humidity of the lowest level
        arrays["dt"] = skt
        # 饱和混合比（Tetens 公式）→ 饱和比湿 / saturation mixing ratio (Tetens) -> specific humidity
        sat_mixing = common.GASCON / 461.51 * 610.78 * np.exp(
            17.2693882 * (skt - 273.16) / (skt - 35.86)
        ) / state.ps
        sat_specific = sat_mixing / (1.0 - (1.0 / (common.GASCON / 461.51) - 1.0) * sat_mixing)
        arrays["dq"] = sat_specific
        return arrays

    def build_restart(self, spectra: dict[str, np.ndarray], state: ModelState,
                      template_records, nstep: int) -> "list[tuple[str, bytes]]":
        """Assemble a restart record list from spectra and surface fields.
        """
        replacements: dict[str, bytes] = {}
        for name in ("st", "sd", "sz", "sq"):
            values = np.asarray(spectra[name], dtype="<f4")
            replacements[name] = restart_io.encode_float(values)
            # the leapfrog old time level starts from the same state
            # 中文说明：蛙跳格式的旧时间层从同一状态出发（冷启动，无时间倾斜）。
            replacements[f"{name}m"] = restart_io.encode_float(values)
        replacements["sp"] = restart_io.encode_float(np.asarray(spectra["sp"], dtype="<f4"))
        replacements["spm"] = replacements["sp"]

        for name, values in self.surface_replacement_arrays(state, template_records).items():
            replacements[name] = restart_io.encode_float(values)

        # clock
        # 中文说明：写入模式时钟（时间步计数器）。
        replacements["nstep"] = restart_io.encode_int(nstep)
        return restart_io.replace_records(template_records, replacements)


# ---------------------------------------------------------------------------
# Nudging target file
# 松弛目标文件
# ---------------------------------------------------------------------------
NUDGE_MAGIC = b"PLASICND"  # file magic


def write_nudge_file(path: Path, entries: "list[tuple[int, dict[str, np.ndarray]]]") -> None:
    """Binary sequence of spectral target states used by the nudging driver.

    中文说明：写出松弛驱动使用的谱目标状态二进制序列。

    Layout: magic ``PLASICND``, int32 header [version, nt, nlev, nrsp], then
    for each time: int64 step followed by sp, st, sq, sz, sd spectra.
    """
    with path.open("wb") as stream:
        stream.write(NUDGE_MAGIC)
        stream.write(np.asarray([1, len(entries), common.NLEV, common.NRSP], dtype="<i4").tobytes())
        for step, spectra in entries:
            stream.write(np.asarray([step], dtype="<i8").tobytes())
            for name in ("sp", "st", "sq", "sz", "sd"):
                stream.write(np.asarray(spectra[name], dtype="<f4").tobytes())


def read_nudge_file(path: Path) -> "list[tuple[int, dict[str, np.ndarray]]]":
    """Read a nudging target file written by :func:`write_nudge_file`.
    """
    with path.open("rb") as stream:
        if stream.read(8) != NUDGE_MAGIC:
            raise ValueError("not a PlaSiC nudging file")
        version, ntimes, nlev, nrsp = np.frombuffer(stream.read(16), dtype="<i4")
        if int(nlev) != common.NLEV or int(nrsp) != common.NRSP:
            raise ValueError(
                "nudging archive grid mismatch: "
                f"archive has nlev={int(nlev)}, nrsp={int(nrsp)}; "
                f"configured model has nlev={common.NLEV}, nrsp={common.NRSP}"
            )
        entries = []
        for _ in range(int(ntimes)):
            step = int(np.frombuffer(stream.read(8), dtype="<i8")[0])
            spectra = {}
            for name in ("sp", "st", "sq", "sz", "sd"):
                # sp 为单层谱，其余变量为逐层谱
                # sp is a single-level spectrum; others are level stacks
                count = nrsp if name == "sp" else nrsp * nlev
                spectra[name] = np.frombuffer(stream.read(4 * count), dtype="<f4").copy()
            entries.append((step, spectra))
    return entries


def write_surface_targets(path: Path, entries: "list[tuple[int, dict[str, np.ndarray]]]") -> None:
    """Store gridded surface targets for the slow part of initialization.

    中文说明：将初始化过程中变化较慢的地表目标场保存为压缩 ``.npz`` 文件。
    ``entries`` 中的每一项是 ``(step, fields)``：``step`` 是模式时间步，
    ``fields`` 是“字段名 → 格点数组”的映射。写出的文件包含一个一维的
    ``step`` 数组，以及每个地表字段对应的三维数组，数组维度统一为
    ``(time, nlat, nlon)``，这样读取程序可以按时间索引直接取得完整格点场。

    English: Save the slowly varying surface target fields used during
    initialization to a compressed ``.npz`` archive. Each item in ``entries``
    is ``(step, fields)``, where ``step`` is the model time step and ``fields``
    maps a field name to its gridded array. The archive stores one one-dimensional
    ``step`` array plus one three-dimensional array per field. The field arrays
    use the common layout ``(time, nlat, nlon)``, allowing the reader to select
    one complete grid at a time with a single time index.
    """
    if not entries:
        raise ValueError("surface target list is empty")

    # 所有时刻应包含同一组字段；字段名按字典序固定，保证输出顺序稳定，
    # 也便于比较不同个例生成的文件。
    # Every time entry is expected to contain the same fields. Sorting the names
    # fixes a deterministic output order and makes files from different cases
    # easier to inspect and compare.
    names = sorted(entries[0][1])

    # 时间步使用 64 位整数保存，避免较长积分或较大时间步发生溢出。
    # Store steps as int64 so long integrations and large step numbers are safe.
    payload = {"step": np.asarray([step for step, _ in entries], dtype=np.int64)}

    for name in names:
        payload[name] = np.stack([np.asarray(values[name], dtype="<f4") for _, values in entries])

    np.savez_compressed(path, **payload)


def read_surface_targets(path: Path) -> list[tuple[int, dict[str, np.ndarray]]]:
    """Read the ``.npz`` surface targets into (step, fields) entries.
    """
    archive = np.load(path)
    steps = np.asarray(archive["step"], dtype=np.int64)
    names = [name for name in archive.files if name != "step"]
    return [(int(step), {name: np.asarray(archive[name][i]).copy() for name in names})
            for i, step in enumerate(steps)]


# ---------------------------------------------------------------------------
# Driver
# 主流程
# ---------------------------------------------------------------------------
def build_case(case: common.ForecastCase, nudge: bool = True) -> None:
    """Build every initial-condition product for one case.

    Writes ``direct.restart`` and, when ``nudge`` is true,
    ``nudge_start.restart`` plus the hourly nudging targets.
    """
    builder = CaseBuilder(case)
    out_dir = common.IC_ROOT / case.key
    out_dir.mkdir(parents=True, exist_ok=True)
    template = restart_io.read_raw_records(TEMPLATE_RESTART)

    weights = top_level_weights()
    background = background_spectra()

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)

        # direct analysis at t0
        # 中文说明：t0 时刻的直接分析初值。
        state_t0 = builder.build_state(builder.t0)
        era5_spectra_t0 = builder.analyse(state_t0, tmpdir / "t0")
        spectra_t0 = blend_top_levels(era5_spectra_t0, background, weights)
        nstep_t0 = common.nstep_for_datetime(builder.t0)
        direct = builder.build_restart(spectra_t0, state_t0, template, nstep_t0)
        restart_io.write_raw_records(out_dir / "direct.restart", direct)
        print(f"[{case.key}] direct IC written (nstep={nstep_t0})")

        if not nudge:
            return

        # analysis 12 h earlier: start of the nudging window
        # 中文说明：窗口起点（t0 前 12 小时）的分析场，作为松弛积分的起始状态。
        start = builder.t0 - dt.timedelta(hours=common.NUDGE_WINDOW_HOURS)
        state_start = builder.build_state(start)
        spectra_start = blend_top_levels(
            builder.analyse(state_start, tmpdir / "start"), background, weights
        )
        nstep_start = common.nstep_for_datetime(start)
        start_restart = builder.build_restart(spectra_start, state_start, template, nstep_start)
        restart_io.write_raw_records(out_dir / "nudge_start.restart", start_restart)
        print(f"[{case.key}] nudging start IC written (nstep={nstep_start})")

        # Hourly targets are pure mapped ERA5 spectra.  The background blend
        # above is used only to assemble the starting restart; applying it to
        # these targets would incorrectly pull the fifth layer toward the
        # fixed template during every cycle.  run_forecasts skips the top four
        # layers when inserting the atmospheric increment.
        # 中文说明：逐小时目标只使用 ERA5 谱系数；模板混合仅用于构造起始
        # restart，不能应用到目标，否则第 5 层仍会逐小时向固定模板松弛。
        entries = []
        surface_entries = []
        for hour in range(common.NUDGE_WINDOW_HOURS + 1):
            when = start + dt.timedelta(hours=hour)
            # 最后一个时刻直接复用 t0 的分析结果，避免重复计算。
            # Reuse the t0 analysis at the final hour to avoid duplicate work.
            state = state_t0 if hour == common.NUDGE_WINDOW_HOURS else builder.build_state(when)
            spectra = (
                era5_spectra_t0
                if hour == common.NUDGE_WINDOW_HOURS
                else builder.analyse(state, tmpdir / f"h{hour:02d}")
            )
            entries.append((common.nstep_for_datetime(when), spectra))
            surface_entries.append((common.nstep_for_datetime(when),
                                    builder.surface_replacement_arrays(state, template)))
        write_nudge_file(out_dir / "nudge_targets.bin", entries)
        write_surface_targets(out_dir / "surface_targets.npz", surface_entries)
        print(f"[{case.key}] nudging targets written ({len(entries)} times)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="*", default=None)
    parser.add_argument("--no-nudge", action="store_true")
    args = parser.parse_args()
    for case in common.CASES:
        if args.cases and case.key not in args.cases:
            continue
        build_case(case, nudge=not args.no_nudge)


if __name__ == "__main__":
    main()
