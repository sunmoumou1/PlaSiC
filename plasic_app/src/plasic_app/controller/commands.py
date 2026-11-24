from __future__ import annotations

# dataclass 用于简化“纯数据结构”类的定义。
# 这里 BuildInvocation 和 RunInvocation 都主要用于封装一次命令调用所需的信息。
from dataclasses import dataclass

import os

from pathlib import Path

import shutil

import socket

import sys

# BuildConfig / ExperimentConfig 是项目内部定义的配置模型。
from plasic_app.config.models import BuildConfig, ExperimentConfig

from plasic_app.paths import app_build_root


@dataclass(frozen=True)
class BuildInvocation:
    """
    描述一次“模型构建命令”的完整调用信息。

    这个类本身并不执行构建，只负责保存执行构建时所需的全部参数。

    frozen=True 表示实例创建完成后不可修改，可以把它视为一个
    immutable value object。
    """

    # 实际需要执行的程序。
    #
    # 在 build_invocation() 中固定为：
    #
    #     "make"
    #
    # 也就是说，真正进行编译的是 Makefile。
    program: str

    # 传递给 program 的命令行参数。
    #
    # 这里使用 tuple 而不是 list，与 frozen dataclass 的不可变设计一致。
    #
    # 例如最终可能类似：
    #
    #     (
    #         "BUILD_DIR=/path/to/build/...",
    #         "NLAT=64",
    #         "NLEV=10",
    #         "NPRO=4",
    #         "MPI=1",
    #         "MPI_CC=mpicc",
    #     )
    arguments: tuple[str, ...]

    # 执行构建命令时使用的 working directory。
    #
    # 对于 make 来说，这一点尤其重要，因为 make 默认会在当前工作目录中寻找 Makefile。
    working_directory: Path

    # 构建命令运行时使用的环境变量。
    #
    # 这里保存的是一个独立的 dict[str, str]，而不是直接持有 os.environ。
    environment: dict[str, str]

    # 本次构建预计生成的 PlaSiC 可执行文件路径。
    # 该字段只是描述“构建完成后 executable 应该在哪里”，
    # BuildInvocation 自身并不检查文件是否实际存在。
    executable: Path

    def as_list(self) -> list[str]:
        """
        将 program 和 arguments 合并成标准的 argv 风格列表。

        例如：

            program = "make"
            arguments = ("NLAT=64", "MPI=0")

        返回：

            ["make", "NLAT=64", "MPI=0"]

        这种形式通常可以直接交给 subprocess 一类的进程执行接口。
        """

        # *self.arguments 使用 iterable unpacking，
        # 把参数 tuple 中的每一项依次展开到 list 中。
        return [self.program, *self.arguments]


@dataclass(frozen=True)
class RunInvocation:
    """
    描述一次“运行 PlaSiC 模型”的完整进程调用信息。

    与 BuildInvocation 类似：
    这个对象本身不启动进程，只负责封装启动进程所需的参数。

    BuildInvocation 与 RunInvocation 分开的原因是：
    构建阶段和运行阶段需要保存的信息略有不同。

    例如：
    - 构建阶段需要记录最终生成的 executable；
    - 运行阶段则只需要描述即将执行的 program / arguments / cwd / env。
    """

    # 最外层实际启动的程序。
    program: str

    # 传递给当前 Python 解释器的参数。
    #
    # run_invocation() 会通过：
    #
    #     python -m plasic_app.controller.process_wrapper ...
    #
    # 先启动项目自己的 process_wrapper，
    # 再由 wrapper 使用 exec 替换为真正的模型命令。
    arguments: tuple[str, ...]

    # 运行进程时的工作目录。
    working_directory: Path

    # 运行进程时使用的环境变量。
    # run_invocation() 会基于当前 os.environ 创建副本，
    # 并额外调整 PYTHONPATH。
    environment: dict[str, str]

    def as_list(self) -> list[str]:
        """
        将 program 和 arguments 合并成一个完整命令列表。

        例如可能返回：

            [
                "/usr/bin/python",
                "-m",
                "plasic_app.controller.process_wrapper",
                "/path/to/plasic.x",
                "--config",
                "/path/to/model.json",
            ]

        如果开启 MPI，中间还会插入 MPI launcher 及进程数参数。
        """

        return [self.program, *self.arguments]


def _conda_bin_directories() -> tuple[Path, ...]:
    """Return nearby Conda ``bin`` directories in a stable order.

    The desktop app is commonly run from a lightweight Python environment,
    while MPI is installed in another environment (for example xesmf_env).
    Conda environments under the same installation are used as a fallback
    after PATH for compiler and launcher wrappers.
    """
    prefixes = [Path(sys.prefix)]
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        prefixes.append(Path(conda_prefix))

    roots: list[Path] = []
    for prefix in prefixes:
        root = (
            prefix.parent
            if prefix.parent.name == "envs"
            else prefix / "envs"
        )
        if root not in roots:
            roots.append(root)

    directories: list[Path] = []
    for root in roots:
        try:
            environments = sorted(root.iterdir())
        except OSError:
            continue
        for environment in environments:
            bin_directory = environment / "bin"
            if bin_directory.is_dir() and bin_directory not in directories:
                directories.append(bin_directory)
    return tuple(directories)


def _resolve_tool(command: str, environment: dict[str, str]) -> str:
    """Resolve a configured tool without requiring its Conda env to be active."""
    expanded = Path(command).expanduser()
    if expanded.is_absolute() or os.sep in command:
        return str(expanded)

    resolved = shutil.which(command, path=environment.get("PATH"))
    if resolved is not None:
        return resolved

    for bin_directory in _conda_bin_directories():
        candidate = bin_directory / command
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    return command


def _macos_default_mpi_interface() -> str | None:
    """Prefer the physical primary Mac interface over VPN ``utun`` devices."""
    if sys.platform != "darwin":
        return None
    try:
        names = {name for _index, name in socket.if_nameindex()}
    except OSError:
        return None
    return "en0" if "en0" in names else None


def _configure_mpi_runtime_environment(environment: dict[str, str]) -> None:
    """Avoid MPICH/libfabric selecting a macOS VPN tunnel as its NIC."""
    interface = environment.get("PLASIC_MPI_INTERFACE")
    if not interface:
        interface = _macos_default_mpi_interface()
    if interface:
        environment.setdefault("FI_SOCKETS_IFACE", interface)


def build_invocation(config: BuildConfig, c_root: Path) -> BuildInvocation:
    """
    根据 BuildConfig 生成一次 PlaSiC 构建调用。

    参数：
        config:
            当前构建配置，决定模型分辨率、垂直层数、MPI 设置、
            编译器等 Makefile 参数。

        c_root:
            PlaSiC C 源代码所在目录。
            后续执行 make 时会把它作为 working directory。

    返回：
        BuildInvocation

    注意：
        这个函数只负责“构造 invocation”，不会实际调用 make。
    """
    # 最终目录结构大致类似：
    #     <project>/build/app/<configuration-specific-directory>/
    build_directory = app_build_root() / config.build_directory_name()

    # 构造传递给 make 的变量参数。
    arguments = [
        # 指定此次编译产物的输出目录。
        f"BUILD_DIR={build_directory}",
        f"NLAT={config.nlat}",
        f"NLEV={config.atmosphere_levels}",
        f"NPRO={config.process_count}",
    ]

    # 将当前 Python 进程的环境变量复制成普通 dict。
    # 使用 dict(os.environ) 的意义是：
    # 创建独立副本，而不是之后直接修改全局的 os.environ。
    # 因此调用方可以把 environment 交给子进程，而不会因为这里的操作
    # 改变当前 Python 主进程的环境。
    environment = dict(os.environ)

    # 根据配置决定构建 MPI 版本还是非 MPI 版本。
    if config.mpi_enabled:
        # Conda MPICH wrappers may remember a build-time compiler name that
        # is unavailable at runtime. Open MPI ignores MPICH_CC, so this is
        # safe for both common implementations and preserves user overrides.
        environment.setdefault("MPICH_CC", config.c_compiler)
        mpi_compiler = _resolve_tool(config.mpi_compiler, environment)
        # MPI 模式：
        # 告诉 Makefile 开启 MPI，并指定 MPI C compiler wrapper。
        arguments.extend(
            [
                # 开启 Makefile 中的 MPI 构建逻辑。
                "MPI=1",
                # 指定 MPI 模式下的 C 编译器。
                # 常见值可能是 mpicc，但具体值完全由配置决定。
                f"MPI_CC={mpi_compiler}",
            ]
        )
    else:
        arguments.extend(
            [
                # 禁用 MPI 编译。
                "MPI=0",
                # 指定普通 C 编译器。
                f"CC={config.c_compiler}",
            ]
        )

    # 约定 make 构建完成后，PlaSiC 主可执行文件位于：
    #     <build_directory>/plasic.x
    # 这里只计算目标路径，不检查该文件现在是否存在。
    executable = build_directory / "plasic.x"

    # 将上面生成的全部信息封装为不可变的 BuildInvocation。
    return BuildInvocation(
        # 最终实际使用 make 执行编译。
        program="make",
        # list 转成 tuple，使 Invocation 中的 arguments 不可变。
        arguments=tuple(arguments),
        # make 从 C 源代码根目录执行。
        working_directory=c_root,
        # 继承当前进程环境变量。
        environment=environment,
        # 保存构建完成后的目标 executable 路径。
        executable=executable,
    )


def run_invocation(
    config: ExperimentConfig,
    c_root: Path,
    executable: Path,
    run_directory: Path,
    initial_restart: Path | None,
) -> RunInvocation:
    """
    根据配置生成一次 PlaSiC 模型运行命令。

    运行参数不落盘为 model_config.json，
    而是直接以命令行选项传给 C runtime：

        plasic.x --steps ... --progress ... --frames ...

    最终进程结构不是简单地直接执行 plasic.x，
    而是：

        当前 Python
            ↓
        python -m plasic_app.controller.process_wrapper
            ↓
        [MPI launcher，可选]
            ↓
        plasic.x [运行参数...]

    process_wrapper 内部会使用 execvpe()，
    因此 wrapper 成功后会被真正的模型进程替换。
    """

    environment = dict(os.environ)

    build = config.build

    if build.mpi_enabled:
        _configure_mpi_runtime_environment(environment)

    source_root = Path(__file__).resolve().parents[2]

    # 检查调用当前程序的环境中是否已经设置 PYTHONPATH。
    # environment.get() 在不存在 PYTHONPATH 时返回 None。
    existing_pythonpath = environment.get("PYTHONPATH")

    # 构造运行 process_wrapper 时使用的 PYTHONPATH。
    environment["PYTHONPATH"] = (
        str(source_root)
        if not existing_pythonpath
        # 如果已经存在 PYTHONPATH：
        # 把 source_root 添加到原有 PYTHONPATH 的最前面。
        else os.pathsep.join((str(source_root), existing_pythonpath))
    )

    # 将 ExperimentConfig 转成 C runtime 认识的一组命令行选项。
    runtime = config.c_runtime_config(c_root, run_directory)

    # 构造真正传递给 PlaSiC executable 的参数。
    # 运行目录中的每一个输出路径都显式给出，
    # 与旧版 model_config.json 中的字段一一对应。
    model_arguments = [
        str(executable),
        "--steps",
        str(runtime["steps"]),
        "--output",
        runtime["output"],
        "--progress",
        runtime["progress"],
        "--frames",
        runtime["frames"],
        "--frame-interval",
        str(runtime["frame_interval"]),
        "--progress-interval",
        str(runtime["progress_interval"]),
        "--co2-forcing",
        runtime["co2_forcing"],
        "--monthly-output",
        runtime["monthly_output"],
        "--calendar",
        runtime["calendar"],
        "--start-year",
        str(runtime["start_year"]),
        "--start-yday",
        str(runtime["start_yday"]),
        "--start-hour",
        str(runtime["start_hour"]),
        "--start-minute",
        str(runtime["start_minute"]),
    ]

    if runtime.get("ozone_forcing"):
        model_arguments.extend(["--ozone-forcing", runtime["ozone_forcing"]])

    if runtime.get("data_dir"):
        model_arguments.extend(["--data-dir", runtime["data_dir"]])

    # 热启动时用 --reset-nstep 重贴时钟；-1 表示保留 restart 中的 nstep。
    if runtime.get("reset_nstep", -1) >= 0:
        model_arguments.extend(["--reset-nstep", str(runtime["reset_nstep"])])

    if runtime.get("dynamics_only"):
        model_arguments.append("--dynamics-only")

    # 大气–海冰耦合间隔只在不等于默认值 1 时传入。C 端缺省同样是 1，
    # 这样尚不认识该选项的旧 plasic.x 仍能运行默认配置。
    if runtime.get("sea_coupling_interval", 1) != 1:
        model_arguments.extend(
            ["--sea-coupling-interval", str(runtime["sea_coupling_interval"])]
        )

    # prepare() 已经把用户选择的 restart 文件链接/复制到运行目录中，
    # 这里优先使用本次运行目录里的 initial.restart。
    if initial_restart is not None:
        model_arguments.extend(["--restart", str(initial_restart)])

    # 如果启用了 MPI，则不能直接运行 plasic.x，
    # 而需要通过 MPI launcher 启动多个模型进程。
    if build.mpi_enabled:
        mpi_launcher = _resolve_tool(build.mpi_launcher, environment)
        command = [
            mpi_launcher,
            "-np",
            str(build.mpi_processes),
            *model_arguments,
        ]
    else:
        command = model_arguments

    # 最终并不直接让外部进程执行 `command`，
    # 而是先启动当前 Python 解释器，并运行项目自己的 process_wrapper。
    #
    # 最终命令大致为：
    #
    # 非 MPI：
    #
    #     <python> \
    #       -m plasic_app.controller.process_wrapper \
    #       <plasic.x> \
    #       [运行参数...]
    #
    # MPI：
    #
    #     <python> \
    #       -m plasic_app.controller.process_wrapper \
    #       <mpi_launcher> \
    #       -np <N> \
    #       <plasic.x> \
    #       [运行参数...]
    #
    # process_wrapper 随后会使用 os.execvpe()，
    # 用 command 中的程序替换 wrapper 自身。
    return RunInvocation(
        program=sys.executable,
        arguments=(
            "-m",
            "plasic_app.controller.process_wrapper",
            *command,
        ),
        working_directory=c_root,
        environment=environment,
    )
