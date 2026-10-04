#!/usr/bin/env python3
"""Shared configuration for the ERA5-initialized PlaSiC forecast experiment.

中文说明：ERA5 初始化的 PlaSiC 预报试验的共享配置模块。

The module centralises the experiment directory layout, model-grid and state
constants, ERA5 request configuration and the three forecast cases.  Grid,
atmospheric-level and ocean settings are read from environment variables so
the workflow follows the model build it is run with.  Defaults are T85/L25
with the slab (non-three-dimensional) ocean.

中文说明：本模块集中管理试验目录结构、模式网格与状态常数、ERA5 下载配置
以及三个预报个例。网格、大气层数和海洋模式都从环境变量读取，使工作流可以
跟随实际编译的模式配置；默认配置为 T85/L25 和单层混合层海洋（不启用三维海洋）。
"""
from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path

import numpy as np

# ----------------------------------------------------------------------------
# Paths
# 路径配置
# ----------------------------------------------------------------------------
# 中文说明：所有目录都从本文件位置推导，保证脚本无论从哪里执行都能找到数据。
# English: all directories are derived from this file's location, so the
# scripts work no matter which working directory they are launched from.
EXPERIMENT_ROOT = Path(__file__).resolve().parents[1]     # experiments/forecast_era5 试验根目录 / experiment root
PROJECT_ROOT = EXPERIMENT_ROOT.parents[1]                 # 仓库根目录 / repository root
SRC_ROOT = PROJECT_ROOT / "src"                           # 模式源码目录 / model source tree
DATA_ROOT = EXPERIMENT_ROOT / "data"                      # 试验数据目录 / experiment data
ERA5_ROOT = DATA_ROOT / "era5"                            # ERA5 原始下载目录 / raw ERA5 downloads
IC_ROOT = DATA_ROOT / "initial_conditions"                # 初始条件目录 / initial conditions
RUN_ROOT = EXPERIMENT_ROOT / "runs"                       # 积分输出目录 / integration outputs
FRAME_ROOT = EXPERIMENT_ROOT / "frames"                   # 模式帧流目录 / model frame streams
FIGURE_ROOT = EXPERIMENT_ROOT / "figures"                 # 论文图件目录 / publication figures
TOOLS_ROOT = EXPERIMENT_ROOT / "tools"                    # 谱分析辅助工具目录 / spectral helper tools

# ----------------------------------------------------------------------------
# Model grid / build configuration
# 模式网格与编译配置
# ----------------------------------------------------------------------------
# The experiment follows the model executable selected below.  Set these
# variables before importing/running the workflow when using another build:
#   PLASIC_NLAT, PLASIC_NLEV, PLASIC_NPRO, PLASIC_MPI
#   PLASIC_OCEAN_MODEL=slab|3d, PLASIC_DEEP_OCEAN_LEVELS
#   PLASIC_BUILD_DIR, PLASIC_TEMPLATE_RESTART, PLASIC_EXECUTABLE
# 中文说明：运行其他编译配置时，在启动脚本前设置上述环境变量即可。三维海洋
# 只有在 ``PLASIC_OCEAN_MODEL=3d`` 时启用，层数由 ``PLASIC_DEEP_OCEAN_LEVELS``
# 单独决定，不再强制为 16 层。


def _env_int(name: str, default: int, minimum: int = 1) -> int:
    """Read a positive integer environment override with a useful error."""
    raw = os.environ.get(name)
    if raw is None or not raw.strip():
        return default
    try:
        value = int(raw)
    except ValueError as exc:
        raise ValueError(f"{name} must be an integer, got {raw!r}") from exc
    if value < minimum:
        raise ValueError(f"{name} must be >= {minimum}, got {value}")
    return value


def _env_bool(name: str, default: bool) -> bool:
    """Read a conventional boolean environment override."""
    raw = os.environ.get(name)
    if raw is None or not raw.strip():
        return default
    value = raw.strip().lower()
    if value in {"1", "true", "yes", "on"}:
        return True
    if value in {"0", "false", "no", "off"}:
        return False
    raise ValueError(f"{name} must be a boolean, got {raw!r}")


# T85 is the default (128 Gaussian latitudes, 256 longitudes).  NLAT is the
# primary grid setting, matching src/Makefile; PLASIC_NTRU is also accepted as
# a convenience when NLAT is omitted.
# 中文说明：默认是 T85（128 个高斯纬圈、256 个经度点）。NLAT 与 src/Makefile
# 保持一致；未设置 NLAT 时也可以用 PLASIC_NTRU 直接指定谱截断。
if "PLASIC_NLAT" in os.environ:
    NLAT = _env_int("PLASIC_NLAT", 128)
    NTRU = _env_int("PLASIC_NTRU", (2 * NLAT - 1) // 3)
else:
    NTRU = _env_int("PLASIC_NTRU", 85)
    # Gaussian grids use an even number of latitudes; choose the smallest
    # matching even NLAT for a truncation supplied on its own.
    NLAT = (3 * NTRU + 2) // 2
    NLAT += NLAT % 2
if NTRU != (2 * NLAT - 1) // 3:
    raise ValueError(
        "PLASIC_NLAT and PLASIC_NTRU must satisfy the model truncation formula"
    )
NLON = 2 * NLAT
NLEV = _env_int("PLASIC_NLEV", 25)
NRSP = (NTRU + 1) * (NTRU + 2)

# The model chooses its atmospheric time step from NLAT.  Keep the same rule
# here so nudge windows, frame intervals and date-to-step conversion stay in
# sync for resolutions whose model day contains whole-hour intervals.
# 中文说明：模式根据 NLAT 选择基本时间步；这里复现同一规则，并要求每天的
# 步数可以整除 24，使松弛窗口、帧间隔和日期到步数的换算保持一致。
if NLAT <= 32:
    MINUTES_PER_STEP = 45
elif NLAT <= 48:
    MINUTES_PER_STEP = 36
else:
    MINUTES_PER_STEP = (30 * 64) // NLAT
NTSPD = _env_int("PLASIC_STEPS_PER_DAY", int(round(86400 / (MINUTES_PER_STEP * 60))))
if NTSPD % 2:
    NTSPD += 1
STEPS_PER_6H = NTSPD // 4
STEPS_PER_HOUR = NTSPD // 24
if NTSPD % 24 != 0:
    raise ValueError(
        f"PLASIC_STEPS_PER_DAY={NTSPD} is not divisible by 24; "
        "the forecast workflow requires whole-hour frame intervals"
    )
FORECAST_DAYS = 10                          # 预报天数 / forecast length in days
FORECAST_STEPS = FORECAST_DAYS * NTSPD      # 预报总步数 / total forecast steps

# Ocean/build selection.  The default deliberately matches a T85/L25 slab
# build; enabling 3-D ocean does not imply any particular layer count.
# 中文说明：默认使用 T85/L25 单层混合层海洋；开启三维海洋后，层数仍由用户单独设置。
OCEAN_MODEL = os.environ.get("PLASIC_OCEAN_MODEL", "slab").strip().lower()
if OCEAN_MODEL not in {"slab", "3d"}:
    raise ValueError("PLASIC_OCEAN_MODEL must be 'slab' or '3d'")
USE_DEEP_OCEAN = OCEAN_MODEL == "3d"
DEEP_OCEAN_LEVELS = _env_int("PLASIC_DEEP_OCEAN_LEVELS", 16)
MODEL_NPRO = _env_int("PLASIC_NPRO", 4)
MODEL_MPI = _env_bool("PLASIC_MPI", MODEL_NPRO > 1)
if not MODEL_MPI and MODEL_NPRO > 1:
    raise ValueError("PLASIC_NPRO must be 1 when PLASIC_MPI is disabled")

# Match src/Makefile's build-directory convention.  Explicit paths always win,
# which is useful when a build uses a project-specific directory name.
_MPI_SUFFIX = f"_MPI_P{MODEL_NPRO}" if MODEL_MPI else ""
_OCEAN_SUFFIX = f"_O{DEEP_OCEAN_LEVELS}" if USE_DEEP_OCEAN else ""
_MODEL_BUILD_DEFAULT = PROJECT_ROOT / "build" / (
    f"T{NTRU}_L{NLEV}{_MPI_SUFFIX}{_OCEAN_SUFFIX}"
)
MODEL_BUILD_DIR = Path(os.environ.get("PLASIC_BUILD_DIR", _MODEL_BUILD_DEFAULT))
# Slab-ocean forecasts use a copy of the packaged restart with the optional
# three-dimensional ocean records removed.  Keeping that file beside the
# experiment makes the model/data contract explicit and prevents a slab run
# from accidentally inheriting a deep-ocean state.  Other grids still fall
# back to their matching source restart unless an explicit path is supplied.
_SLAB_TEMPLATE_DEFAULT = EXPERIMENT_ROOT / f"earth_t{NTRU}_l{NLEV}_slab.restart"
_SOURCE_TEMPLATE_DEFAULT = SRC_ROOT / f"data/earth_t{NTRU}_l{NLEV}.restart"
TEMPLATE_RESTART = Path(os.environ.get(
    "PLASIC_TEMPLATE_RESTART",
    _SLAB_TEMPLATE_DEFAULT
    if not USE_DEEP_OCEAN and _SLAB_TEMPLATE_DEFAULT.exists()
    else _SOURCE_TEMPLATE_DEFAULT,
))

PLARAD = 6371220.0                          # 地球半径（m）/ planetary radius (m)
SIDEREAL_DAY = 86164.0916                   # 恒星日长度（s）/ sidereal day (s)
GASCON = 287.0                              # 干空气气体常数（J kg-1 K-1）/ gas constant of dry air
GRAVITY = 9.80665                           # 重力加速度（m s-2）/ gravitational acceleration
TWO_PI = 2.0 * np.pi                        # 2π
CV = PLARAD * (TWO_PI / SIDEREAL_DAY)       # 模式 Robert 形风速的无量纲化系数 a·Ω / scaling a*Omega for the model's dimensionless Robert-form winds


def sigma_levels(nlev: int = NLEV) -> tuple[np.ndarray, np.ndarray]:
    """Reproduce initpm_vertical() (stretched grid, sigma_mode neither 1 nor -1).

    Returns (sigma_half, sigma_full); sigma_half[k] is the lower interface of
    layer k (the model top sigma=0 is implicit) and sigma_full the layer centre.

    中文说明：复现模式初始化例程 initpm_vertical() 的 σ 层设置（拉伸网格，
    sigma_mode 既不是 1 也不是 -1 的情形）。

    返回 (sigma_half, sigma_full)：sigma_half[k] 是第 k 层的下界面（模式
    顶部的 sigma=0 是隐含的，不在数组里），sigma_full 是层中心。
    """
    level = np.arange(1, nlev + 1, dtype=np.float64)   # 层序号 1..nlev / level indices 1..nlev
    s = level / nlev                                   # 归一化层序号 / normalized level number
    # 拉伸函数：低层较薄、高层较厚，与模式 initpm_vertical() 中的多项式一致。
    # Stretching polynomial matching the model's initpm_vertical().
    sigma_half = 0.75 * s + 1.75 * s**3 - 1.5 * s**4
    sigma_thickness = np.empty(nlev)
    sigma_thickness[0] = sigma_half[0]                 # 最底层厚度从模式顶起算 / bottom layer thickness measured from the top
    sigma_thickness[1:] = np.diff(sigma_half)          # 其余各层厚度 / thickness of the remaining layers
    sigma_full = np.empty(nlev)
    sigma_full[0] = 0.5 * sigma_half[0]                # 第 0 层中心 / centre of the first layer
    sigma_full[1:] = 0.5 * (sigma_half[:-1] + sigma_half[1:])  # 相邻界面中点 / midpoint between adjacent interfaces
    return sigma_half, sigma_full


SIGMA_HALF, SIGMA_FULL = sigma_levels()

# ----------------------------------------------------------------------------
# ERA5 request configuration
# ERA5 下载请求配置
# ----------------------------------------------------------------------------
# All pressure levels from 1000 to 10 hPa that ERA5 offers, used for building
# the initial condition and the nudging targets.
# 中文说明：ERA5 提供的 1000–10 hPa 全部气压层，用于构建初始条件和松弛同化目标。
ERA5_IC_LEVELS = [
    1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700, 650, 600,
    550, 500, 450, 400, 350, 300, 250, 225, 200, 175, 150, 125, 100, 70, 50,
    30, 20, 10,
]
# Reduced set for verification, covering the troposphere and lower stratosphere
# at standard levels.
# 中文说明：用于检验的精简层次集合，以标准气压层覆盖对流层和平流层下部。
ERA5_VERIF_LEVELS = [
    1000, 975, 950, 925, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450,
    400, 350, 300, 250, 200, 150, 100, 70, 50, 30, 10,
]

ERA5_GRID = 1.5  # 下载网格分辨率（度），接近默认目标网格 / download grid resolution in degrees


@dataclass(frozen=True)
class ForecastCase:
    """Definition of one forecast case (immutable).

    中文说明：单个预报个例的定义（不可变数据类）。
    """

    key: str                # 个例标识 / case identifier
    init_time: str          # ISO "YYYY-MM-DDTHH" 初始时刻 / ISO initialization time
    label_zh: str           # 中文标签 / Chinese label
    description_zh: str     # 中文描述 / Chinese description
    # 关注区域 (west, east, south, north)，默认全球
    # focus region (west, east, south, north); global by default
    focus_region: tuple[float, float, float, float] = (-180.0, 180.0, -90.0, 90.0)


# 三个检验个例分别代表隆冬、盛夏和台风活跃期，覆盖不同的环流形势。
# The three cases sample midwinter, midsummer and an active typhoon period,
# covering different circulation regimes.
CASES: list[ForecastCase] = [
    ForecastCase(
        key="winter2020",
        init_time="2020-01-15T00",
        label_zh="个例 A：2020-01-15 00 UTC（北半球隆冬）",
        description_zh="北半球冬季，极涡强盛、中纬度急流活跃，东亚一次寒潮过程前后。",
    ),
    ForecastCase(
        key="summer2020",
        init_time="2020-07-15T00",
        label_zh="个例 B：2020-07-15 00 UTC（北半球盛夏）",
        description_zh="北半球夏季，副热带高压控制、对流活跃，检验模式对暖季环流与湿过程的预报。",
    ),
    ForecastCase(
        key="typhoon2020",
        init_time="2020-09-01T00",
        label_zh="个例 C：2020-09-01 00 UTC（台风 Maysak）",
        description_zh="西北太平洋台风 Maysak（2020 年第 9 号台风）发展盛期，检验模式对强天气系统的预报。",
        focus_region=(100.0, 180.0, 0.0, 60.0),
    ),
]

# Relaxation time scale of the Newtonian nudging initialization (hours).
# 中文说明：牛顿松弛（nudging）初始化的松弛时间尺度（小时）。
NUDGE_WINDOW_HOURS = 12     # 同化窗口长度 / assimilation window length
NUDGE_TAU_HOURS = 1.0       # 大气谱系数的松弛时间尺度 / relaxation time scale for atmospheric spectra
# Surface reservoirs evolve on a much slower time scale than the free
# atmosphere.  Relaxing them over 24 h removes the old bug where a nudged
# restart retained the t0−window land/ocean state indefinitely.
# 中文说明：地表状态量的演变比自由大气慢得多，因此用 24 小时的时间尺度对其
# 松弛。这样修复了旧版本的一个缺陷：松弛后的重启文件中，陆地/海洋状态一直
# 停留在 t0−窗口 时刻而不再更新。
NUDGE_SURFACE_TAU_HOURS = 24.0  # 地表状态的松弛时间尺度 / relaxation time scale for surface state


def init_datetime(case: ForecastCase):
    """Return the case initialization time as a ``datetime.datetime``.

    中文说明：返回个例的初始时刻，类型为 ``datetime.datetime``。
    """
    import datetime as dt

    return dt.datetime.fromisoformat(case.init_time)


def nstep_for_datetime(when, calendar: str = "gregorian") -> int:
    """Model step counter for a Gregorian datetime (epoch = 0000-01-01 00:00).

    中文说明：把公历日期时间换算为模式步数计数器（起算时刻为 0000-01-01 00:00）。

    The counter advances with the configured model time step and applies the
    proleptic Gregorian leap-year rule.
    中文说明：计数器按当前配置的模式时间步递增，并采用推测公历
    （proleptic Gregorian calendar）的闰年规则。
    """
    import datetime as dt

    if calendar != "gregorian":
        raise ValueError(calendar)
    # 该年之前的天数（含闰年修正）/ days before this year, including leap-year correction
    year_days = 365 * when.year + (when.year + 3) // 4 - (when.year + 99) // 100 + (when.year + 399) // 400
    # day of year, 1-based / 一年中的第几天（从 1 开始计数）
    day_of_year = (when - dt.datetime(when.year, 1, 1)).days + 1
    days = year_days + (day_of_year - 1)
    # 总步数 = 天数 × 每天步数 + 小时步数 + 分钟折算步数
    # total = days * steps/day + hour * steps/hour + minute-converted steps
    return days * NTSPD + when.hour * STEPS_PER_HOUR + when.minute * NTSPD // 60
