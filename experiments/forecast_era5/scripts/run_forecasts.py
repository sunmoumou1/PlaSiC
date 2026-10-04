#!/usr/bin/env python3
"""Run the PlaSiC 10-day forecasts (direct and nudged initial conditions).

Every forecast is integrated in two segments so that the frame stream carries
both an hourly view of the first day (for the initial-shock analysis) and the
6-hourly view of the remaining nine days:

Two initial conditions are compared for every case:

* ``direct``  the ERA5 analysis at t0, mapped to the configured model grid
  (old time level equal to the current one);
* ``nudged``  a 12-h Newtonian-relaxation initialization (tau = 1 h) toward the
  hourly ERA5 analyses, implemented by 1-h restart recycling, ending at t0.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_ic  # noqa: E402
import common  # noqa: E402
import frame_io  # noqa: E402
import restart_io  # noqa: E402

# 模式可执行文件和并行进程数均由 common 的配置统一管理；仍可用
# PLASIC_EXECUTABLE / PLASIC_NPRO 在环境变量中覆盖。默认是 T85/L25、P4、
# 不启用三维海洋的构建。
# The executable and process count come from common so they stay aligned with
# the selected grid and ocean build.  PLASIC_EXECUTABLE / PLASIC_NPRO remain
# available as environment overrides.
PROGRAM = Path(os.environ.get("PLASIC_EXECUTABLE", common.MODEL_BUILD_DIR / "plasic.x"))
NPRO = common.MODEL_NPRO  # MPI 进程数 / number of MPI ranks
SEGMENT1_STEPS = common.NTSPD                      # 第一段 24 小时 / segment 1: 24 h
SEGMENT1_FRAME_INTERVAL = common.STEPS_PER_HOUR    # 第一段帧间隔 1 小时 / segment 1 frame interval: 1 h
SEGMENT2_STEPS = common.FORECAST_STEPS - SEGMENT1_STEPS  # 第二段 216 小时 / segment 2: remaining hours
SEGMENT2_FRAME_INTERVAL = common.STEPS_PER_6H      # 第二段帧间隔 6 小时 / segment 2 frame interval: 6 h

# Keep the four uppermost sigma layers freely evolving during nudging.  ERA5
# relaxation starts at the fifth layer (zero-based index 4).
# 中文说明：松弛循环中最上面的四个 sigma 层保持模式自由演化；从第 5 层
# （下标 4）开始才应用 ERA5 松弛。
NUDGE_FREE_TOP_LEVELS = 4


def run_model(restart: Path, steps: int, output: Path, frames: Path | None,
              frame_interval: int, log: Path, co2_ppm: float = 410.0) -> None:
    """Run one model segment from a restart file.

    ``frames`` may be ``None`` for the short nudging cycles, which do not need
    a frame stream; stdout/stderr are redirected to ``*.stdout`` next to the
    progress log.
    """
    log.parent.mkdir(parents=True, exist_ok=True)

    if common.MODEL_MPI:
        launcher = ["mpirun", "-np", str(NPRO), str(PROGRAM)]
    else:
        launcher = [str(PROGRAM)]

    command = launcher + [
        "--restart", str(restart),
        "--co2-ppm", repr(co2_ppm),
        "--steps", str(steps),
        "--output", str(output),
        "--progress", str(log),
        "--progress-interval", str(max(1, steps // 8)),
    ]
    if frames is not None:
        command += ["--frames", str(frames), "--frame-interval", str(frame_interval)]
    started = time.time()
    with log.with_suffix(".stdout").open("w") as stream:
        subprocess.run(command, check=True, stdout=stream, stderr=subprocess.STDOUT)

    print(f"    {steps} steps in {time.time() - started:.0f} s -> {output.name}")


def blend_restart(model_restart: Path, target: dict[str, np.ndarray],
                  gamma: float, output: Path) -> None:
    """Insert the analysis increment into the current time level only.

    The model's own previous (filtered) leapfrog level is kept, so the
    integration stays continuous and the Robert-Asselin filter damps the
    increment smoothly instead of being restarted from a cold state every hour.
    """
    records = restart_io.read_raw_records(model_restart)
    current = dict(records)
    replacements: dict[str, bytes] = {}
    for name in ("st", "sd", "sz", "sq"):
        model = restart_io.decode_float(current[name]).reshape(common.NLEV, common.NRSP)
        target_values = np.asarray(target[name], dtype=np.float64).reshape(
            common.NLEV, common.NRSP
        )
        # Leave levels 0..3 untouched.  From level 4 downward:
        # new = target + gamma * (model - target), i.e. 63.2% target and
        # 36.8% post-step model for the default one-hour relaxation.
        blended = model.copy()
        blended[NUDGE_FREE_TOP_LEVELS:] = target_values[NUDGE_FREE_TOP_LEVELS:] + gamma * (
            model[NUDGE_FREE_TOP_LEVELS:] - target_values[NUDGE_FREE_TOP_LEVELS:]
        )
        replacements[name] = restart_io.encode_float(blended.reshape(-1))

    # Surface pressure has no vertical levels and keeps the existing update.
    model = restart_io.decode_float(current["sp"])
    target_values = np.asarray(target["sp"], dtype=np.float64)
    replacements["sp"] = restart_io.encode_float(
        target_values + gamma * (model - target_values)
    )
    restart_io.write_raw_records(
        output, restart_io.replace_records(records, replacements)
    )


def blend_surface_restart(model_restart: Path, target: dict[str, np.ndarray],
                          gamma: float, output: Path) -> None:
    """Relax slowly varying surface records toward the matching ERA5 hour.

    Missing/land sentinel values are copied from the model, so a finite
    analysis never contaminates the model's land/ocean bookkeeping.
    """
    records = restart_io.read_raw_records(model_restart)
    current = dict(records)
    replacements: dict[str, bytes] = {}
    for name, target_values in target.items():
        if name not in current:
            continue
        model = restart_io.decode_float(current[name]).astype(np.float64)
        target_array = np.asarray(target_values, dtype=np.float64).reshape(model.shape)
        # 哨兵值（≥1e10 或非有限值）不参与松弛 / sentinels (>=1e10 or non-finite) are not relaxed
        valid = np.isfinite(target_array) & (np.abs(target_array) < 1.0e10)
        blended = np.where(valid, target_array + gamma * (model - target_array), model)
        replacements[name] = restart_io.encode_float(blended)
    restart_io.write_raw_records(output, restart_io.replace_records(records, replacements))


def frames_match(path: Path, expected_count: int, expected_interval: int) -> bool:
    """Return true only for a frame stream with the requested cadence.

    This prevents a pre-6-hour archive from being silently reused after the
    output schedule changes.
    中文说明：这样可以防止旧的（非 6 小时节奏的）归档在输出计划变更后被悄悄复用。
    """
    if not path.exists():
        return False
    try:
        grouped = frame_io.load_stream(path)
        frames = next(iter(grouped.values()))
        steps = np.asarray([frame.step for frame in frames], dtype=np.int64)
    except (OSError, ValueError, StopIteration):
        return False
    return (steps.size == expected_count
            and (steps.size < 2 or np.all(np.diff(steps) == expected_interval)))


def forecast(ic_path: Path, run_dir: Path, frame_dir: Path, tag: str,
             force: bool) -> None:
    """Run the two-segment 10-day forecast for one initial condition.
    """
    frames1 = frame_dir / f"{tag}_seg1.frames"
    frames2 = frame_dir / f"{tag}_seg2.frames"

    expected1 = SEGMENT1_STEPS // SEGMENT1_FRAME_INTERVAL
    expected2 = SEGMENT2_STEPS // SEGMENT2_FRAME_INTERVAL

    if (frames_match(frames1, expected1, SEGMENT1_FRAME_INTERVAL)
            and frames_match(frames2, expected2, SEGMENT2_FRAME_INTERVAL)
            and not force):
        print(f"    [skip] {tag} frames exist")
        return

    print(f"  [{tag}] segment 1: 0-24 h, hourly frames")
    run_model(ic_path, SEGMENT1_STEPS, run_dir / f"{tag}_seg1.restart",
              frames1, SEGMENT1_FRAME_INTERVAL, run_dir / f"{tag}_seg1.jsonl")

    print(f"  [{tag}] segment 2: 24-240 h, 6-hourly frames")
    # 第二段从第一段结束时的重启文件继续。/ segment 2 continues from segment 1's restart.
    run_model(run_dir / f"{tag}_seg1.restart", SEGMENT2_STEPS,
              run_dir / f"{tag}_final.restart", frames2, SEGMENT2_FRAME_INTERVAL,
              run_dir / f"{tag}_seg2.jsonl")


def nudged_ic(case: common.ForecastCase, run_dir: Path, ic_dir: Path, force: bool = False) -> Path:
    """Build the nudged initial condition by 1-hour restart recycling.

    Each hour, the model is first integrated for one hour from the current
    restart. The atmospheric spectral records are then relaxed toward the
    matching hourly ERA5 target with weight ``gamma``. If surface targets are
    available, the slowly varying surface records are relaxed in a second step
    with their longer time scale. The final state is written as
    ``nudged.restart``.

    The first target entry corresponds to ``nudge_start.restart`` itself, so
    the loop starts at entry 1 and advances one hour before applying each later
    target. The step check at the end of every cycle keeps the restart clock and
    target archive on the same time axis.

    中文说明：目标序列的第 0 项与 ``nudge_start.restart`` 属于同一时刻，因而
    循环从第 1 项开始；每次先前进 1 小时，再应用新的目标场。每个循环末尾的
    步数检查用于确保重启文件时钟与目标文件时间轴没有错位。
    """
    nudged_ic_path = ic_dir / "nudged.restart"

    if nudged_ic_path.exists() and not force:
        print(f"    [skip] {nudged_ic_path.name} exists")
        return nudged_ic_path

    # nudge_targets.bin contains atmospheric spectral targets at the hourly
    # model steps, including both the window start and the final t0 target.
    print(f"  nudging initialization ({common.NUDGE_WINDOW_HOURS} h, "
          f"tau = {common.NUDGE_TAU_HOURS} h)")
    entries = build_ic.read_nudge_file(ic_dir / "nudge_targets.bin")

    # Surface targets were added after the atmospheric target format. Keep the
    # optional read for compatibility with older IC directories: those runs
    # still perform atmospheric nudging and simply keep their surface state.
    # 中文说明：地表目标文件是在大气目标格式之后加入的，因此允许文件缺失，
    # 以兼容旧的初值目录；旧目录仍会执行大气松弛，但保留原有地表状态。
    surface_path = ic_dir / "surface_targets.npz"
    surface_entries = build_ic.read_surface_targets(surface_path) if surface_path.exists() else []

    gamma = float(np.exp(-1.0 / common.NUDGE_TAU_HOURS))

    # The restart produced by build_ic is the state at the beginning of the
    # window (entry 0). Each loop replaces state with the output of the latest
    # cycle, so the next one continues from the already relaxed state.
    # 中文说明：build_ic 生成的重启文件位于松弛窗口起点（目标序列第 0 项）。
    # 每轮循环都把 state 更新为最新松弛后的文件，下一轮从该连续状态继续。
    state = ic_dir / "nudge_start.restart"
    for hour in range(1, len(entries)):
        step, target = entries[hour]
        cycle_out = run_dir / f"nudge_cycle_{hour:02d}.restart"

        # First advance the physical model by one hour without writing frames;
        # these short cycles only provide the forecast state to be nudged.
        # 中文说明：先向前积分一个小时且不输出帧；这些短循环只需要生成
        # 待松弛的模式状态，不需要保存预报帧流。
        run_model(state, common.STEPS_PER_HOUR, cycle_out, None, 1,
                  run_dir / f"nudge_cycle_{hour:02d}.jsonl")

        blended = run_dir / f"nudge_blend_{hour:02d}.restart"

        # Relax only the current atmospheric time level. blend_restart preserves
        # the previous filtered leapfrog level so the next cycle remains
        # dynamically continuous.
        # 中文说明：这里只松弛当前的大气时间层；blend_restart 会保留上一层
        # 已滤波的蛙跳状态，使下一轮积分保持动力学连续。
        blend_restart(cycle_out, target, gamma, blended)

        if surface_entries:
            surface_target = surface_entries[hour][1]
            surface_blended = run_dir / f"nudge_surface_blend_{hour:02d}.restart"
            surface_gamma = float(np.exp(-1.0 / common.NUDGE_SURFACE_TAU_HOURS))
            blend_surface_restart(blended, surface_target, surface_gamma, surface_blended)
            state = surface_blended
        else:
            state = blended

        nstep = restart_io.decode_int(
            dict(restart_io.read_raw_records(state))["nstep"]
        )
        assert nstep == step, (nstep, step)

    # Copy the final recycled state to the stable output name. Re-reading and
    # writing raw records preserves the complete restart record layout.
    # 中文说明：将最后一轮循环状态写入固定输出文件名；通过原始记录读写，
    # 保留重启文件的完整记录布局。
    restart_io.write_raw_records(nudged_ic_path, restart_io.read_raw_records(state))
    print(f"    nudged IC written -> {nudged_ic_path.name}")
    return nudged_ic_path


def run_case(case: common.ForecastCase, methods: "list[str]", force: bool) -> None:
    """Run the requested initialization methods for one case.
    """
    run_dir = common.RUN_ROOT / case.key
    frame_dir = common.FRAME_ROOT / case.key
    ic_dir = common.IC_ROOT / case.key
    run_dir.mkdir(parents=True, exist_ok=True)
    frame_dir.mkdir(parents=True, exist_ok=True)

    if "direct" in methods:
        print(f"[{case.key}] direct forecast")
        forecast(ic_dir / "direct.restart", run_dir, frame_dir, "direct", force)
    if "nudged" in methods:
        print(f"[{case.key}] nudged forecast")
        ic = nudged_ic(case, run_dir, ic_dir, force=force)
        forecast(ic, run_dir, frame_dir, "nudged", force)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="*", default=None)
    parser.add_argument("--methods", nargs="*", default=["direct", "nudged"])
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    for case in common.CASES:
        if args.cases and case.key not in args.cases:
            continue
        run_case(case, args.methods, args.force)


if __name__ == "__main__":
    main()
