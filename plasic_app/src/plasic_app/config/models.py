from __future__ import annotations

# asdict:
#   将 dataclass 实例递归转换成普通 Python dict。
# dataclass:
#   自动生成 __init__、__repr__、__eq__ 等常用方法。
# field:
#   用于对 dataclass 字段进行更细粒度配置。
#   后面 ExperimentConfig 中会使用 default_factory，
#   防止多个 ExperimentConfig 实例共享同一个可变对象。
from dataclasses import asdict, dataclass, field

# SHA-256 哈希算法。
from hashlib import sha256

import json

from pathlib import Path

from typing import Any

# ---------------------------------------------------------------------------
# 不同谱分辨率对应的网格尺寸与谱截断参数
# ---------------------------------------------------------------------------
#
# 每个 value 都是一个三元组：
#
#     (nlon, nlat, truncation)
#
# 例如：
#
#     "T21": (64, 32, 21)
#
# 表示：
# - 经度方向网格点数 nlon = 64
# - 纬度方向网格点数 nlat = 32
# - 谱截断阶数 truncation = 21
#
# T21 / T31 / T42 / T85 是谱模式中常见的三角截断分辨率表示方式。
RESOLUTION_GRIDS = {
    "T21": (64, 32, 21),
    "T31": (96, 48, 31),
    "T42": (128, 64, 42),
    "T85": (256, 128, 85),
}


# ---------------------------------------------------------------------------
# 每个模拟日需要执行的动力学时间步数，与 C 模型的默认步长
# PLASIC_DEFAULT_STEP_SECONDS（src/state/plasicmod_state.h）一致：
# T21 30、T31 18、T42 12、T85 7.5 分钟。
# Dynamical steps per simulated day, matching the default time step of the C
# model, PLASIC_DEFAULT_STEP_SECONDS in src/state/plasicmod_state.h:
# 30, 18, 12 and 7.5 minutes at T21, T31, T42 and T85.
# ---------------------------------------------------------------------------
STEPS_PER_DAY = {
    "T21": 48,
    "T31": 80,
    "T42": 120,
    "T85": 192,
}


# ---------------------------------------------------------------------------
# 不同分辨率默认均执行 cold start。已有 restart 仍可由用户显式选择。
# ---------------------------------------------------------------------------
# 空字符串意味着不向 C runtime 提供 --restart；T21/T31/T42/T85 语义一致。
RESTART_FILES = {
    "T21": "",
    "T31": "",
    "T42": "",
    "T85": "",
}


# ---------------------------------------------------------------------------
# 与默认 restart 文件匹配的大气垂直层数
# ---------------------------------------------------------------------------
ATMOSPHERE_LEVELS = {
    "T21": 10,
    "T31": 12,
    "T42": 16,
    "T85": 25,
}


# ===========================================================================
# BuildConfig
# ===========================================================================
#
# BuildConfig 描述“编译模型”时需要确定的配置。
#
# 可以把整个程序配置粗略区分为：
#
#     BuildConfig
#         编译期设置
#
#     RunConfig
#         运行期设置
#
#     VisualizationConfig
#         可视化设置
#
# BuildConfig 中的设置如果发生变化，
# 很可能意味着需要重新编译 C 模型。
@dataclass
class BuildConfig:
    resolution: str = "T21"
    atmosphere_levels: int = 10
    mpi_enabled: bool = False
    mpi_processes: int = 2
    c_compiler: str = "cc"
    mpi_compiler: str = "mpicc"
    mpi_launcher: str = "mpirun"

    @property
    def nlon(self) -> int:

        return RESOLUTION_GRIDS[self.resolution][0]

    @property
    def nlat(self) -> int:

        return RESOLUTION_GRIDS[self.resolution][1]

    @property
    def truncation(self) -> int:

        return RESOLUTION_GRIDS[self.resolution][2]

    @property
    def steps_per_day(self) -> int:

        return STEPS_PER_DAY[self.resolution]

    @property
    def spectral_rows(self) -> int:
        """
        计算谱空间中的 row 数量。

        当前公式为：

            (truncation + 1) * (truncation + 2)

        例如 T21：

            truncation = 21

            spectral_rows
            = (21 + 1) * (21 + 2)
            = 22 * 23
            = 506

        该值只作诊断展示；C runtime 允许进程数不整除谱行数，最后一个
        rank 的部分块由 padding 补齐。
        """

        return (self.truncation + 1) * (self.truncation + 2)

    @property
    def process_count(self) -> int:

        return self.mpi_processes if self.mpi_enabled else 1

    @property
    def default_restart(self) -> str:

        return RESTART_FILES.get(self.resolution, "")

    @property
    def default_atmosphere_levels(self) -> int:

        return ATMOSPHERE_LEVELS.get(self.resolution, 10)

    def validate(self) -> list[str]:

        # 初始化错误列表。
        errors: list[str] = []

        # ---------------------------------------------------------------
        # 1. 检查 resolution 是否支持
        # ---------------------------------------------------------------
        if self.resolution not in RESOLUTION_GRIDS:
            errors.append(f"Unsupported resolution: {self.resolution}")
            return errors

        # ---------------------------------------------------------------
        # 2. 检查大气层数
        # ---------------------------------------------------------------
        if not 5 <= self.atmosphere_levels <= 50:
            errors.append("Atmospheric levels must be between 5 and 50")

        # ---------------------------------------------------------------
        # 3. MPI 相关校验
        # ---------------------------------------------------------------
        # 谱分块允许最后一个 rank 持有带 padding 的部分块（C runtime 用
        # MPI_Allgatherv / MPI_Gatherv 只通信有效槽位），因此谱行数不再约束
        # 进程数，只有纬度行数必须能被进程数整除。
        if self.mpi_enabled:
            if self.mpi_processes < 2:
                errors.append("MPI requires at least 2 processes")
            elif self.nlat % self.mpi_processes != 0:
                errors.append(
                    "MPI process count must divide latitude rows "
                    f"{self.nlat}"
                )

        # 返回所有收集到的错误。
        return errors

    def build_key(self) -> str:
        """
        根据完整 BuildConfig 内容生成稳定的 build key。

        这个 key 可以用于区分不同编译配置。

        例如以下配置变化：

            T21 -> T42
            atmosphere_levels 10 -> 16
            MPI False -> True
            compiler gcc -> clang

        都会导致 build_key 发生变化。

        这样可以避免不同编译配置互相覆盖。
        """
        encoded = json.dumps(asdict(self), sort_keys=True, ensure_ascii=False).encode(
            "utf-8"
        )
        return sha256(encoded).hexdigest()[:12]

    def build_directory_name(self) -> str:
        """
        根据编译配置生成 build directory 名称。

        非 MPI 示例可能类似：

            T21_L10_ab12cd34ef56

        MPI 示例可能类似：

            T21_L10_MPI_P4_ab12cd34ef56
        """

        suffix = f"_MPI_P{self.mpi_processes}" if self.mpi_enabled else ""

        return (
            f"{self.resolution}_L{self.atmosphere_levels}"
            f"{suffix}_{self.build_key()}"
        )


# ===========================================================================
# RunConfig
# ===========================================================================
#
# RunConfig 描述的是：
#
#     模型已经编译完成以后，
#     一次具体 simulation run 所需要的运行时配置。
#
# 与 BuildConfig 不同，
# RunConfig 的变化通常不意味着一定要重新编译 C 模型。
@dataclass
class RunConfig:

    steps: int = 32

    # 冷启动日历与起始日期：日历二选一，起始日期用"年内第几天 + 时:分"。
    # Cold-start calendar and date: one of two calendars, date given as
    # day-of-year plus hour:minute.
    calendar: str = "gregorian"

    start_year: int = 0

    start_yday: int = 1

    start_hour: int = 0

    start_minute: int = 0

    # 热启动时钟重贴标签：-1 表示保留 restart 中的 nstep。
    # Warm-start clock rebase: -1 keeps the nstep stored in the restart.
    reset_nstep: int = -1

    co2_forcing: str = (
        "data/forcing/CO2/"
        "mole-fraction-of-carbon-dioxide-in-air_input4MIPs_"
        "GHGConcentrations_CMIP_UoM-CMIP-1-2-0_"
        "gn-15x360deg_000001-201412.csv"
    )

    restart: str = ""

    data_dir: str = "data/surface"

    dynamics_only: bool = False

    # 每隔多少 integration steps 输出一次 live frame。
    # 默认每 4 步输出一次。
    frame_interval: int = 4

    progress_interval: int = 1

    # 大气–海冰耦合间隔：每隔多少个大气步运行一次海冰和海洋；默认 1 即每步耦合。
    # Atmosphere-sea ice coupling interval: sea ice and the ocean run once
    # every this many atmospheric steps; the default 1 couples every step.
    sea_coupling_interval: int = 1

    ozone_forcing: str = ""

    def validate(self) -> list[str]:
        """
        校验单次运行配置。
        """

        errors: list[str] = []

        if self.steps < 1:
            errors.append("Integration steps must be at least 1")

        if self.calendar not in ("gregorian", "360day"):
            errors.append("Calendar must be 'gregorian' or '360day'")

        if self.start_year < 0:
            errors.append("Calendar start year cannot be negative")

        if self.start_yday < 1:
            errors.append("Cold-start day of year must be at least 1")

        if not 0 <= self.start_hour <= 23:
            errors.append("Cold-start hour must be between 0 and 23")

        if not 0 <= self.start_minute <= 59:
            errors.append("Cold-start minute must be between 0 and 59")

        if self.reset_nstep < -1:
            errors.append(
                "reset_nstep must be -1 (disabled) or a non-negative step"
            )

        if not self.co2_forcing:
            errors.append("A monthly CO2 forcing CSV is required")

        if self.frame_interval < 1:
            errors.append("Live frame interval must be at least 1")

        if self.progress_interval < 1:
            errors.append("Progress interval must be at least 1")

        if self.sea_coupling_interval < 1:
            errors.append("Sea ice coupling interval must be at least 1")

        # 只有 cold start 需要 surface data；restart 分支允许完全省略。
        if not self.data_dir and not self.restart:
            errors.append("A surface data directory is required")

        return errors


# ===========================================================================
# VisualizationConfig
# ===========================================================================
#
# 这一部分只描述 UI / 可视化层面的配置，
# 不直接控制底层物理模式。
@dataclass
class VisualizationConfig:

    buffer_frames: int = 320

    # 是否自动确定绘图 levels。
    # 例如对于 contour / color scale，
    # True 时可能根据实际数据范围自动设置。
    auto_levels: bool = True

    # 地球显示方式。
    # 当前仅允许：
    #     "map"
    #     "globe"
    view_mode: str = "map"

    # 默认选择显示的变量名称。
    selected_variable: str = "surface_temperature"

    # 默认选择的垂直层。
    # 0 的具体含义由后续可视化/变量协议实现决定。
    selected_level: int = 0

    def validate(self) -> list[str]:
        """
        校验可视化相关配置。
        """

        errors: list[str] = []

        if not 2 <= self.buffer_frames <= 512:
            errors.append("Buffered frame groups must be between 2 and 512")

        if self.view_mode not in {"map", "globe"}:
            errors.append("Earth view must be map or globe")

        return errors


# ===========================================================================
# ExperimentConfig
# ===========================================================================
#
# ExperimentConfig 是整个实验配置的顶层对象。
#
# 它把：
#
#     BuildConfig
#     RunConfig
#     VisualizationConfig
#
# 三部分组合在一起。
#
# 可以理解为：
#
# ExperimentConfig
# ├── build
# ├── run
# └── visualization
@dataclass
class ExperimentConfig:

    # field(default_factory=BuildConfig)
    #
    # 这里故意没有写：
    #
    #     build: BuildConfig = BuildConfig()
    #
    # 而是使用 default_factory。
    #
    # 这样每创建一个 ExperimentConfig，
    # 都会单独创建一个新的 BuildConfig 实例，
    # 避免多个 ExperimentConfig 共享同一个对象。
    build: BuildConfig = field(default_factory=BuildConfig)

    run: RunConfig = field(default_factory=RunConfig)

    visualization: VisualizationConfig = field(default_factory=VisualizationConfig)

    def validate(self) -> list[str]:

        return (
            self.build.validate() + self.run.validate() + self.visualization.validate()
        )

    def to_dict(self) -> dict[str, Any]:

        return asdict(self)

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "ExperimentConfig":
        """
        从普通 dict 恢复 ExperimentConfig。

        这是 to_dict() 的对应反向构造入口。

        这里允许 payload 缺少某些配置块。

        如果缺少：
            "build"

        则：
            payload.get("build", {})
            -> {}

        然后：
            BuildConfig(**{})
            -> 使用 BuildConfig 自身默认值。

        因此它天然支持部分配置。
        """

        build_payload = dict(payload.get("build", {}))
        # 旧版预设可能仍包含已删除的海洋层数；忽略该字段以保持可读。
        build_payload.pop("ocean_levels", None)

        run_payload = dict(payload.get("run", {}))
        # 旧版预设可能仍包含已删除的植被模式；忽略该字段以保持可读。
        run_payload.pop("vegetation", None)
        # 旧版日历字段已被 start_yday/start_hour/start_minute 与 reset_nstep 取代。
        run_payload.pop("start_month", None)
        run_payload.pop("reset_calendar", None)

        return cls(
            build=BuildConfig(**build_payload),
            run=RunConfig(**run_payload),
            visualization=VisualizationConfig(**payload.get("visualization", {})),
        )

    @staticmethod
    def _resolve_model_path(value: str, c_root: Path) -> Path:
        """
        将模型配置中的路径解析成实际文件系统路径。
        """

        path = Path(value).expanduser()

        return path if path.is_absolute() else c_root / path

    def c_runtime_config(self, c_root: Path, run_directory: Path) -> dict[str, Any]:

        # 构造 C runtime 接收的基础配置。
        config: dict[str, Any] = {
            "output": str(run_directory / "final.restart"),
            "progress": str(run_directory / "progress.txt"),
            "frames": str(run_directory / "frames.pstream"),
            "steps": self.run.steps,
            "calendar": self.run.calendar,
            "start_year": self.run.start_year,
            "start_yday": self.run.start_yday,
            "start_hour": self.run.start_hour,
            "start_minute": self.run.start_minute,
            "reset_nstep": self.run.reset_nstep,
            "co2_forcing": str(
                self._resolve_model_path(self.run.co2_forcing, c_root)
            ),
            "monthly_output": str(run_directory / "monthly"),
            "frame_interval": self.run.frame_interval,
            "progress_interval": self.run.progress_interval,
            "dynamics_only": self.run.dynamics_only,
            "sea_coupling_interval": self.run.sea_coupling_interval,
        }

        if self.run.ozone_forcing:
            config["ozone_forcing"] = str(
                self._resolve_model_path(self.run.ozone_forcing, c_root)
            )

        if self.run.data_dir:
            config["data_dir"] = str(
                self._resolve_model_path(self.run.data_dir, c_root)
            )

        # ------------------------------------------------------------------
        # Restart / Cold Start 处理
        # ------------------------------------------------------------------
        if self.run.restart:

            config["restart"] = str(self._resolve_model_path(self.run.restart, c_root))

        return config
