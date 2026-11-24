#!/bin/bash

# =============================================================================
# PlaSiC T21/L10 + 三维海洋（16 层）历史气候实验流水线
# PlaSiC T21/L10 + 3-D ocean (16 levels) historical experiment pipeline
#
# 整个实验分为三个主要阶段：
# The experiment consists of three main stages:
#
#   Phase 1a:
#     从冷启动海洋状态开始进行 500 年 spin-up，
#     分为 5 个连续的 100 年模拟段。
#     Perform a 500-year ocean spin-up from a cold-start state,
#     split into five consecutive 100-year segments.
#
#   Phase 1b:
#     在 spin-up 第 500 年附近继续积分，并保存 10 个相邻年份的 restart，
#     作为候选初始条件。
#     Continue integration around spin-up year 500 and save 10 restart states
#     from consecutive years as candidate initial conditions.
#
#   Phase 1c:
#     根据近地表气温，从 10 个候选状态中选择：
#       A = 表面最暖的状态
#       B = 表面最冷的状态
#     Select two states from the 10 candidates based on surface air temperature:
#       A = warmest surface state
#       B = coldest surface state
#
#   Phase 2:
#     分别从 A 和 B 两个初始条件启动：
#       1. 1850-1859 的 10 年固定 CO2 控制实验
#       2. 1850-2014 的历史 CO2 强迫实验
#     From both initial conditions A and B, launch:
#       1. a 10-year fixed-CO2 control experiment for 1850-1859
#       2. a historical CO2-forced experiment for 1850-2014
#
#   Phase 3:
#     如果 analyze.py 存在，则自动执行后处理、统计分析和绘图。
#     If analyze.py exists, automatically perform post-processing,
#     statistical analysis, and figure generation.
#
# 使用方法：
# Usage:
#
#   bash run_pipeline.sh spinup
#       只运行 Phase 1：spin-up + 候选初始条件生成 + 初始条件选择
#       Run Phase 1 only: spin-up + candidate IC generation + IC selection
#
#   bash run_pipeline.sh historical
#       只运行 Phase 2 和 Phase 3；
#       要求 Phase 1 已经完成并生成 selection.txt。
#       Run Phase 2 and Phase 3 only;
#       Phase 1 must already have produced selection.txt.
#
#   bash run_pipeline.sh all
#       依次运行完整流水线：Phase 1 + Phase 2 + Phase 3。
#       Run the complete pipeline: Phase 1 + Phase 2 + Phase 3.
#
# 可恢复性：
# Resumability:
#
#   脚本设计为可断点续跑。已经存在的 spin-up restart 或候选初始条件
#   将被自动跳过，因此中途终止后可以重新执行脚本。
#   The pipeline is designed to be resumable. Existing spin-up restart files
#   and candidate initial-condition files are automatically skipped, allowing
#   the script to be rerun after an interruption.
# =============================================================================


# -----------------------------------------------------------------------------
# Shell 基本设置
# Basic shell settings
# -----------------------------------------------------------------------------

# 如果使用了未定义变量，则立即报错。
# Exit with an error when an undefined variable is referenced.
#
# 注意：这里没有启用 "set -e"，因此普通命令失败不会自动退出；
# 关键步骤通过显式的 "|| { ...; exit 1; }" 进行错误处理。
# Note: "set -e" is intentionally not enabled here, so ordinary command
# failures do not automatically terminate the script; critical steps are
# handled explicitly using "|| { ...; exit 1; }".
set -u


# -----------------------------------------------------------------------------
# 路径与全局配置
# Paths and global configuration
# -----------------------------------------------------------------------------

# 获取当前脚本所在目录的绝对路径。
# Get the absolute path of the directory containing this script.
#
# "$0" 是当前脚本路径；
# dirname "$0" 得到脚本所在目录；
# cd + pwd 将其转换为绝对路径。
# "$0" is the current script path;
# dirname "$0" gives its parent directory;
# cd + pwd converts it to an absolute path.
HERE="$(cd "$(dirname "$0")" && pwd)"


# PlaSiC 源代码目录。
# PlaSiC source-code directory.
#
# 假设目录结构类似：
# Assumed directory layout:
#
#   project/
#   ├── src/
#   └── experiments/
#       └── run_pipeline.sh
#
SRC="$(cd "$HERE/../src" && pwd)"


# 编译后 PlaSiC 可执行文件的位置。
# Location of the compiled PlaSiC executable.
#
# 配置含义：
# Configuration encoded in the directory name:
#
#   T21      : T21 水平谱截断
#              T21 horizontal spectral truncation
#
#   L10      : 10 个大气垂直层
#              10 atmospheric vertical levels
#
#   MPI_P2   : MPI 并行版本，使用 2 个进程
#              MPI version configured for 2 processes
#
#   O16      : 三维海洋模型具有 16 个深海层
#              3-D ocean model with 16 deep-ocean levels
#
BUILD="$HERE/../build/T21_L10_MPI_P2_O16/plasic.x"


# Python 解释器路径。
# Python interpreter path.
#
# 如果环境变量 PLASIC_PYTHON 已定义，则优先使用它；
# 否则使用下面给出的默认 Conda Python。
# If the environment variable PLASIC_PYTHON is defined, use it;
# otherwise fall back to the default Conda Python below.
PYTHON="${PLASIC_PYTHON:-/opt/homebrew/Caskroom/miniconda/base/envs/mainpython/bin/python}"


# MPI 启动命令。
# MPI launcher command.
#
# --bind-to none：
#   不将 MPI rank 强制绑定到特定 CPU core。
#   Do not bind individual MPI ranks to specific CPU cores.
#
# -np 2：
#   每个 PlaSiC 实例使用 2 个 MPI rank。
#   Run each PlaSiC instance with 2 MPI ranks.
MPIRUN="mpirun --bind-to none -np 2"


# 工业化前固定 CO2 浓度，单位 ppm。
# Fixed preindustrial CO2 concentration in ppm.
#
# 该值用于 spin-up 和 control 实验。
# This value is used for the spin-up and control experiments.
CO2_FIXED="284.316853"


# 历史 CO2 强迫文件。
# Historical CO2 forcing file.
#
# 该 CSV 文件包含随时间变化的 CO2 浓度，
# 用于 1850-2014 historical 实验。
# This CSV file contains time-varying CO2 concentrations
# used by the 1850-2014 historical experiment.
CO2_CSV="data/forcing/CO2/mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_CMIP_UoM-CMIP-1-2-0_gn-15x360deg_000001-201412.csv"


# -----------------------------------------------------------------------------
# 模型日历换算辅助函数
# Model-calendar helper functions
# -----------------------------------------------------------------------------
#
# 当前流水线固定使用公历（Gregorian）和 T21/L10 分辨率：
# The pipeline uses the Gregorian calendar and the T21/L10 resolution:
#
#   NTSPD = 48：每个太阳日 48 个动力步（T21 默认步长 30 分钟/步）。
#   NTSPD = 48: 48 dynamical steps per solar day (the T21 default step of
#   30 minutes).
#
# 发布的 500 年结果是用旧的 45 分钟步长（NTSPD = 32）算出的；现在的流水线
# 使用新的默认步长，因此会多花约一半的时间。
# The published 500-year results were computed with the earlier 45-minute
# step (NTSPD = 32); the pipeline now uses the new default step and takes
# about half as long again.
#
# 纪元固定为 nstep = 0 ≡ 第 0 年 1 月 1 日 00:00，与模型内部一致。
# The epoch is nstep = 0 ≡ year 0, January 1, 00:00, as in the model.
NTSPD=48


# 返回第 0 年 1 月 1 日到 year 年 1 月 1 日的公历天数。
# Return the number of Gregorian days from year 0 to January 1 of `year`.
days_before_year() {

    local y="$1"

    echo $(( 365 * y + (y + 3) / 4 - (y + 99) / 100 + (y + 399) / 400 ))
}


# 返回 year 年 1 月 1 日 00:00 对应的模型步数（nstep）。
# Return the model step (nstep) of January 1, 00:00 of `year`.
year_start_step() {

    echo $(( $(days_before_year "$1") * NTSPD ))
}


# 返回从 sy 年 1 月 1 日到 ey 年 12 月 31 日 24:00 的模型步数。
# Return the number of model steps from January 1 of `sy` through the end of `ey`.
segment_steps() {

    local sy="$1"
    local ey="$2"

    echo $(( ( $(days_before_year $((ey + 1))) - $(days_before_year "$sy") ) * NTSPD ))
}


# 从第一个命令行参数读取要运行的阶段。
# Read the requested pipeline phase from the first command-line argument.
#
# 如果没有给出参数，则默认执行 "all"。
# If no argument is supplied, default to "all".
PHASE="${1:-all}"


# -----------------------------------------------------------------------------
# 日志函数
# Logging function
# -----------------------------------------------------------------------------

log() {

    # 在终端打印带时间戳的日志，同时追加到 pipeline.log。
    # Print a timestamped message to the terminal and append it to pipeline.log.
    #
    # "$*" 表示传递给 log() 的所有参数。
    # "$*" represents all arguments passed to log().
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" \
        | tee -a "$HERE/pipeline.log"
}


# -----------------------------------------------------------------------------
# 检查并构建 PlaSiC
# Check and build PlaSiC
# -----------------------------------------------------------------------------

require_build() {

    # "-x" 检查 BUILD 指向的文件是否存在且可执行。
    # "-x" tests whether BUILD exists and is executable.
    if [ ! -x "$BUILD" ]; then

        log "building $BUILD"

        # 进入源代码目录并编译 PlaSiC。
        # Enter the source directory and compile PlaSiC.
        #
        # MPI=1
        #   启用 MPI。
        #   Enable MPI support.
        #
        # NPRO=2
        #   编译配置使用 2 个 MPI 进程。
        #   Configure the build for 2 MPI processes.
        #
        # NLAT=32
        #   T21 对应 32 个高斯纬度。
        #   T21 uses 32 Gaussian latitudes.
        #
        # NLEV=10
        #   大气具有 10 个垂直层。
        #   Use 10 atmospheric vertical levels.
        #
        # OCEAN_MODEL=3d
        #   使用三维海洋模型。
        #   Enable the 3-D ocean model.
        #
        # DEEP_OCEAN_LEVELS=16
        #   海洋垂直方向使用 16 层。
        #   Use 16 deep-ocean vertical levels.
        #
        # 编译输出全部追加到 pipeline.log。
        # Append all build output to pipeline.log.
        (
            cd "$SRC" &&
            make all \
                MPI=1 \
                NPRO=2 \
                NLAT=32 \
                NLEV=10 \
                OCEAN_MODEL=3d \
                DEEP_OCEAN_LEVELS=16
        ) >>"$HERE/pipeline.log" 2>&1 \
            || {
                log "BUILD FAILED"
                exit 1
            }
    fi
}


# -----------------------------------------------------------------------------
# PlaSiC 模型统一启动函数
# Common PlaSiC model launcher
# -----------------------------------------------------------------------------

run_model() {

    # 第一个参数固定作为日志文件路径。
    # The first argument is always interpreted as the log-file path.
    local logfile="$1"

    # 删除第一个参数，使剩余参数全部作为 PlaSiC 命令行参数。
    # Remove the first argument so all remaining arguments can be passed
    # directly to PlaSiC.
    shift

    # 在源代码目录运行模型。
    # Run the model from the source directory.
    #
    # 这样相对路径形式的 forcing/data 文件可以从 SRC 正确解析。
    # This allows relative forcing/data paths to be resolved from SRC.
    #
    # stdout 和 stderr 都追加到指定 logfile。
    # Both stdout and stderr are appended to the requested logfile.
    (
        cd "$SRC" &&
        $MPIRUN "$BUILD" "$@"
    ) >>"$logfile" 2>&1
}


# =============================================================================
# Phase 1a
# 500 年海洋 spin-up：5 × 100 年
# 500-year ocean spin-up: 5 × 100-year segments
# =============================================================================
#
# 模型日历年份使用：
# Model calendar years:
#
#   Segment 1: 1850-1949
#   Segment 2: 1950-2049
#   Segment 3: 2050-2149
#   Segment 4: 2150-2249
#   Segment 5: 2250-2349
#
# 这里的 1850-2349 主要作为模型内部日历标签。
# 由于 CO2 始终固定为工业化前值，因此并不代表真实历史强迫。
# Here, 1850-2349 primarily serves as the model's internal calendar labeling.
# Since CO2 remains fixed at the preindustrial value, these years do not
# represent actual historical forcing.
# =============================================================================

phase_spinup() {

    # 创建 spin-up 输出目录；如果已经存在则不报错。
    # Create the spin-up output directory; do nothing if it already exists.
    mkdir -p "$HERE/spinup"

    log "phase 1a: 500-year ocean spin-up (5 x 100 years)"


    # i   : 当前 segment 编号，范围 0-4
    # sy  : segment 起始年份
    # ey  : segment 结束年份
    # inp : 输入 restart
    # out : 输出 restart
    #
    # i   : current segment index, ranging from 0 to 4
    # sy  : segment start year
    # ey  : segment end year
    # inp : input restart file
    # out : output restart file
    local i sy ey inp out


    # 总共运行 5 个连续的 100 年 segment。
    # Run five consecutive 100-year segments.
    for i in 0 1 2 3 4; do

        # 计算当前 segment 起始年份。
        # Compute the start year of the current segment.
        #
        # i=0 -> 1850
        # i=1 -> 1950
        # ...
        sy=$((1850 + 100 * i))


        # 计算当前 segment 结束年份。
        # Compute the end year of the current segment.
        #
        # 每个 segment 包含完整的 100 个年份。
        # Each segment contains exactly 100 complete calendar years.
        ey=$((1949 + 100 * i))


        # 定义当前 segment 的最终 restart 文件名。
        # Define the final restart filename for this segment.
        #
        # Segment 1 -> spinup_y100.restart
        # Segment 2 -> spinup_y200.restart
        # ...
        # Segment 5 -> spinup_y500.restart
        out="$HERE/spinup/spinup_y$((100 * (i + 1))).restart"


        # 如果目标 restart 已存在，认为该 segment 已经完成，
        # 因此直接跳过，实现断点续跑。
        # If the target restart already exists, assume this segment has
        # completed successfully and skip it, enabling resumability.
        if [ -f "$out" ]; then

            log "  skip segment $((i + 1))/5 ($sy-$ey): $out exists"
            continue
        fi


        # 第一个 segment 从 PlaSiC 默认 Earth restart 启动。
        # The first segment starts from the standard PlaSiC Earth restart.
        if [ "$i" -eq 0 ]; then

            inp="$SRC/data/earth_t21_l10.restart"

        else

            # 后续 segment 从上一个 100 年 segment 的输出继续运行。
            # Every later segment continues from the output restart produced
            # by the preceding 100-year segment.
            inp="$HERE/spinup/spinup_y$((100 * i)).restart"
        fi


        log "  segment $((i + 1))/5: $sy-$ey  (in: $inp)"


        # 启动当前 100 年 spin-up segment。
        # Launch the current 100-year spin-up segment.
        run_model "$HERE/spinup/run.log" \
            \
            --restart "$inp" \
            \
            --co2-ppm "$CO2_FIXED" \
            \
            --reset-nstep "$(year_start_step "$sy")" \
            --steps "$(segment_steps "$sy" "$ey")" \
            \
            --progress "$HERE/spinup/progress_y$((100 * (i + 1))).jsonl" \
            \
            --progress-interval 17532 \
            \
            --output "$out" \
            \
            || {
                # 当前 segment 失败时立即终止整个流水线。
                # Abort the entire pipeline immediately if this segment fails.
                log "  SPIN-UP SEGMENT $((i + 1)) FAILED"
                exit 1
            }

        # 参数说明：
        # Parameter description:
        #
        # --restart "$inp"
        #   从指定 restart 状态启动模型。
        #   Initialize the model from the specified restart state.
        #
        # --co2-ppm "$CO2_FIXED"
        #   整个 spin-up 使用固定 CO2 浓度。
        #   Keep CO2 fixed throughout the spin-up.
        #
        # --reset-nstep "$(year_start_step "$sy")"
        #   把 restart 中的 nstep 直接改到 sy 年 1 月 1 日 00:00，
        #   只重贴日期标签，不改变物理状态。
        #   Rebase the nstep read from the restart to January 1, 00:00 of
        #   year sy without touching the physical state.
        #
        # --steps "$(segment_steps "$sy" "$ey")"
        #   从 sy 年 1 月 1 日积分到 ey 年 12 月 31 日 24:00，
        #   步数按公历天数 × NTSPD 算出。
        #   Integrate from January 1 of sy through the end of ey, with the
        #   step count derived from Gregorian day counts × NTSPD.
        #
        # --progress ...
        #   将运行过程中的诊断信息写入 JSON Lines 文件。
        #   Write runtime progress diagnostics to a JSON Lines file.
        #
        # --progress-interval 17532
        #   每隔指定模型步数/内部时间间隔写一次 progress。
        #   Write progress information at the specified model-step/internal
        #   interval.
        #
        # --output "$out"
        #   模拟结束时保存最终 restart。
        #   Save the final restart state at the end of the segment.
    done


    # 所有 5 个 segment 完成后，spinup_y500.restart 即为
    # 500 年 spin-up 后的平衡候选状态。
    # After all five segments complete, spinup_y500.restart represents
    # the candidate equilibrated state after 500 years of spin-up.
    log "phase 1a complete: $HERE/spinup/spinup_y500.restart"
}


# =============================================================================
# Phase 1b
# 在 spin-up 第 500 年附近保存多个平衡初始条件
# Save multiple equilibrated initial conditions around spin-up year 500
# =============================================================================

phase_initial_conditions() {

    # 创建候选初始条件目录。
    # Create the candidate initial-condition directory.
    mkdir -p "$HERE/initial_conditions"


    # Phase 1a 最终的 500 年 restart。
    # Final 500-year restart produced by Phase 1a.
    local src500="$HERE/spinup/spinup_y500.restart"


    # 如果 spinup_y500.restart 不存在，则无法继续。
    # The pipeline cannot continue if spinup_y500.restart is missing.
    [ -f "$src500" ] || {
        log "missing $src500"
        exit 1
    }


    # 将第 500 年 spin-up 状态保存为第一个候选初始条件 ic_00。
    # Save the year-500 spin-up state as the first candidate IC, ic_00.
    cp -f "$src500" "$HERE/initial_conditions/ic_00.restart"


    log "phase 1b: saving initial conditions ic_00 (year 500) .. ic_09 (year 509)"


    # prev：
    #   前一个候选状态，用作下一年的 restart 输入。
    #
    # j：
    #   候选 IC 编号 1-9。
    #
    # yr：
    #   实际用于 PlaSiC 日历的年份。
    #
    # out：
    #   当前候选 restart 输出文件。
    #
    # prev:
    #   Previous candidate state, used as the restart for the next year.
    #
    # j:
    #   Candidate IC index from 1 to 9.
    #
    # yr:
    #   Calendar year used by PlaSiC.
    #
    # out:
    #   Current candidate restart output file.
    local prev="$src500" j yr out


    # 在第 500 年状态基础上继续逐年积分 9 年，
    # 从而得到总共 10 个相邻年份的候选状态。
    # Continue integrating one year at a time for nine additional years,
    # producing 10 candidate states from consecutive years in total.
    for j in 1 2 3 4 5 6 7 8 9; do

        # spinup_y500.restart 对应模型日历 2349 年末。
        # 因此 j=1 使用 2350 年，j=9 使用 2358 年。
        # spinup_y500.restart corresponds to the end of model year 2349.
        # Therefore j=1 uses year 2350 and j=9 uses year 2358.
        yr=$((2349 + j))


        # 生成两位数字编号：
        # Generate a zero-padded two-digit index:
        #
        #   1 -> ic_01.restart
        #   ...
        #   9 -> ic_09.restart
        out="$HERE/initial_conditions/ic_$(printf '%02d' "$j").restart"


        # 如果该候选状态已经存在，则跳过重新计算。
        # If this candidate state already exists, skip recomputation.
        if [ -f "$out" ]; then

            log "  skip ic_$(printf '%02d' "$j") ($yr): exists"


            # 即使跳过，也必须把 prev 更新为该文件，
            # 这样下一年仍然可以从正确的 restart 继续。
            # Even when skipping, update prev to this file so the next
            # year continues from the correct restart.
            prev="$out"

            continue
        fi


        log "  ic_$(printf '%02d' "$j"): year $yr"


        # 从前一年候选状态继续运行一年。
        # Integrate one additional year from the previous candidate state.
        run_model "$HERE/initial_conditions/run.log" \
            \
            --restart "$prev" \
            \
            --co2-ppm "$CO2_FIXED" \
            \
            --reset-nstep "$(year_start_step "$yr")" \
            --steps "$(segment_steps "$yr" "$yr")" \
            \
            --progress "$HERE/initial_conditions/ic_$(printf '%02d' "$j").progress.jsonl" \
            \
            --progress-interval 17532 \
            \
            --output "$out" \
            \
            || {
                # 任意候选 IC 生成失败时立即停止。
                # Abort immediately if generation of any candidate IC fails.
                log "  IC $j FAILED"
                exit 1
            }


        # 当前输出成为下一次循环的输入。
        # The current output becomes the restart input for the next iteration.
        prev="$out"
    done


    # 统计实际生成的候选 restart 数量。
    # Count the number of candidate restart files actually present.
    #
    # wc -l 计算文件数量；
    # tr -d ' ' 删除 wc 可能输出的前导空格。
    # wc -l counts the files;
    # tr -d ' ' removes any leading spaces produced by wc.
    log "phase 1b complete: $(ls "$HERE/initial_conditions"/ic_*.restart | wc -l | tr -d ' ') initial conditions"
}


# =============================================================================
# Phase 1c
# 从候选状态中选择暖、冷两个初始条件
# Select warm and cold initial conditions from the candidate states
# =============================================================================

phase_select() {

    # selection.txt 保存最终选择结果。
    # selection.txt stores the final IC-selection result.
    local sel="$HERE/initial_conditions/selection.txt"


    # 如果 selection.txt 已存在，则认为此前已经完成选择。
    # If selection.txt already exists, assume the selection has already
    # been completed.
    if [ -f "$sel" ]; then

        log "phase 1c: selection already exists"
        return
    fi


    log "phase 1c: selecting two initial conditions"


    # 调用 Python 脚本 select_ics.py 对候选 IC 进行分析。
    # Invoke select_ics.py to analyze the candidate initial conditions.
    #
    # --dir
    #   候选 restart 所在目录。
    #   Directory containing candidate restart files.
    #
    # --out
    #   最终选择结果文件。
    #   Output file containing the selection.
    #
    # Python 的 stdout 和 stderr 都写入 pipeline.log。
    # Both stdout and stderr from Python are appended to pipeline.log.
    "$PYTHON" "$HERE/select_ics.py" \
        --dir "$HERE/initial_conditions" \
        --out "$sel" \
        >>"$HERE/pipeline.log" 2>&1 \
        || {
            log "SELECTION FAILED"
            exit 1
        }


    # -------------------------------------------------------------------------
    # 将被选中的两个状态固定复制为 icA.restart 和 icB.restart，
    # 从而使 Phase 2 不需要知道原始 ic_XX.restart 的编号。
    #
    # Materialize the two selected states as icA.restart and icB.restart,
    # so Phase 2 does not need to know their original ic_XX.restart indices.
    # -------------------------------------------------------------------------

    local a b


    # 假设 selection.txt 中包含类似：
    # Assume selection.txt contains entries similar to:
    #
    #   A ic_03.restart ...
    #   B ic_07.restart ...
    #
    # awk 找到第一列为 A/B 的行，并取第二列文件名。
    # awk finds rows whose first column is A/B and extracts the filename
    # from the second column.
    a="$(awk '$1=="A"{print $2}' "$sel")"
    b="$(awk '$1=="B"{print $2}' "$sel")"


    # A = 表面最暖候选状态。
    # A = candidate with the warmest surface state.
    cp -f \
        "$HERE/initial_conditions/$a" \
        "$HERE/initial_conditions/icA.restart"


    # B = 表面最冷候选状态。
    # B = candidate with the coldest surface state.
    cp -f \
        "$HERE/initial_conditions/$b" \
        "$HERE/initial_conditions/icB.restart"


    log "  A = $a  -> icA.restart"
    log "  B = $b  -> icB.restart"
}


# =============================================================================
# Phase 2
# 从 A/B 两个初始条件分别运行 control 和 historical
# Run control and historical experiments from both A/B initial conditions
# =============================================================================
#
# 总共启动四个 PlaSiC 实验：
# Four PlaSiC experiments are launched in total:
#
#   control_A
#   historical_A
#   control_B
#   historical_B
#
# 每个模型实例自身使用：
# Each model instance uses:
#
#   mpirun -np 2
#
# 并且四个任务通过 Bash 后台运行机制并发执行。
# All four jobs are launched concurrently using Bash background processes.
# =============================================================================

phase_historical() {

    # Phase 2 依赖 Phase 1c 的选择结果。
    # Phase 2 depends on the selection generated by Phase 1c.
    local sel="$HERE/initial_conditions/selection.txt"


    [ -f "$sel" ] || {
        log "missing selection; run 'spinup' phase first"
        exit 1
    }


    # 确保 PlaSiC 已经编译。
    # Ensure that the PlaSiC executable has been built.
    require_build


    # tag：
    #   当前 ensemble/initial-condition 标签 A 或 B。
    #
    # icfile：
    #   当前使用的 restart 初始条件。
    #
    # cdir：
    #   control 实验目录。
    #
    # hdir：
    #   historical 实验目录。
    #
    # tag:
    #   Current ensemble/initial-condition label, A or B.
    #
    # icfile:
    #   Restart file used as the current initial condition.
    #
    # cdir:
    #   Control-experiment directory.
    #
    # hdir:
    #   Historical-experiment directory.
    local tag icfile cdir hdir


    # 分别处理 A 和 B 两套初始条件。
    # Process initial conditions A and B separately.
    for tag in A B; do


        # 根据 tag 选择对应的 restart 文件。
        # Select the appropriate restart file based on tag.
        if [ "$tag" = "A" ]; then

            icfile="$HERE/initial_conditions/icA.restart"

        else

            icfile="$HERE/initial_conditions/icB.restart"
        fi


        # 当前 control 实验目录。
        # Current control-experiment directory.
        cdir="$HERE/control_$tag"


        # 当前 historical 实验目录。
        # Current historical-experiment directory.
        hdir="$HERE/historical_$tag"


        # 为 control 和 historical 创建月平均输出目录。
        # Create monthly-output directories for both experiments.
        mkdir -p "$cdir/monthly" "$hdir/monthly"


        log "phase 2: launching control_$tag (1850-1859) and historical_$tag (1850-2014)"


        # ---------------------------------------------------------------------
        # Control experiment
        # 控制实验
        #
        # 使用固定工业化前 CO2：
        # Use fixed preindustrial CO2:
        #
        #   1850-1859，共 10 年
        #   1850-1859, for a total of 10 years
        # ---------------------------------------------------------------------
        (
            run_model "$cdir/run.log" \
                \
                --restart "$icfile" \
                \
                --co2-ppm "$CO2_FIXED" \
                \
                --reset-nstep "$(year_start_step 1850)" \
                --steps "$(segment_steps 1850 1859)" \
                \
                --monthly-output "$cdir/monthly" \
                \
                --progress "$cdir/progress.jsonl" \
                \
                --progress-interval 480 \
                \
                --output "$cdir/final.restart"


            # 保存模型进程退出状态。
            # Save the model process exit status.
            #
            # 0 通常表示成功；
            # 非 0 表示运行过程中发生错误。
            # 0 normally means success;
            # nonzero indicates an error.
            echo $? > "$cdir/exit_code"

        ) &

        # 最后的 "&" 表示该 control 实验在后台运行，
        # Bash 不等待它完成，而是立即继续启动其他实验。
        # The trailing "&" launches the control experiment in the background,
        # so Bash immediately proceeds to launch other experiments.


        # ---------------------------------------------------------------------
        # Historical experiment
        # 历史强迫实验
        #
        # 使用随时间变化的历史 CO2：
        # Use time-varying historical CO2 forcing:
        #
        #   1850-2014，共 165 年
        #   1850-2014, for a total of 165 years
        # ---------------------------------------------------------------------
        (
            run_model "$hdir/run.log" \
                \
                --restart "$icfile" \
                \
                --co2-forcing "$CO2_CSV" \
                \
                --reset-nstep "$(year_start_step 1850)" \
                --steps "$(segment_steps 1850 2014)" \
                \
                --monthly-output "$hdir/monthly" \
                \
                --progress "$hdir/progress.jsonl" \
                \
                --progress-interval 480 \
                \
                --output "$hdir/final.restart"


            # 保存 historical 模型进程退出状态。
            # Save the exit status of the historical run.
            echo $? > "$hdir/exit_code"

        ) &

    done


    # 等待当前 shell 启动的所有后台实验完成。
    # Wait until all background experiments launched by this shell finish.
    #
    # 因为 A/B 各启动 control + historical，
    # 因此这里最多同时等待四个 PlaSiC 实例。
    # Since A/B each launch one control and one historical experiment,
    # this waits for up to four PlaSiC instances.
    wait


    log "phase 2 complete"


    # 汇总四个实验的 exit code。
    # Summarize the exit codes of all four experiments.
    for tag in A B; do

        # "2>/dev/null" 防止 exit_code 文件缺失时把 cat 的错误信息
        # 输出到终端。
        # "2>/dev/null" suppresses cat error messages if an exit_code file
        # is unexpectedly missing.
        log "  control_$tag exit=$(cat "$HERE/control_$tag/exit_code" 2>/dev/null)"

        log "  historical_$tag exit=$(cat "$HERE/historical_$tag/exit_code" 2>/dev/null)"
    done
}


# =============================================================================
# Phase 3
# 后处理、统计分析与绘图
# Post-processing, statistical analysis, and figure generation
# =============================================================================

phase_analysis() {

    # 只有 analyze.py 存在时才执行分析。
    # Run the analysis only if analyze.py exists.
    if [ -f "$HERE/analyze.py" ]; then

        log "phase 3: analysis and figures"


        # 切换到实验目录执行 analyze.py，
        # 使 Python 脚本可以方便地使用相对路径访问：
        #
        #   control_A/
        #   control_B/
        #   historical_A/
        #   historical_B/
        #
        # Change to the experiment directory before executing analyze.py,
        # allowing the Python script to conveniently access relative paths:
        #
        #   control_A/
        #   control_B/
        #   historical_A/
        #   historical_B/
        (
            cd "$HERE" &&
            "$PYTHON" analyze.py
        ) >>"$HERE/analysis.log" 2>&1 \
            || {
                log "ANALYSIS FAILED (see analysis.log)"
                exit 1
            }

        # 后处理失败时返回非零状态，避免把不完整的图表误报为成功。
        # A post-processing failure is fatal so incomplete figures are not
        # reported as a successful experiment.
    fi
}


# =============================================================================
# 主程序入口
# Main program entry point
# =============================================================================

# 记录流水线启动时间和请求执行的阶段。
# Record pipeline startup and the requested phase.
log "===== pipeline start (phase=$PHASE) ====="


# 无论执行 spinup、historical 还是 all，
# 首先检查模型 executable 是否存在；不存在则自动编译。
# Regardless of whether spinup, historical, or all is requested,
# first ensure that the model executable exists, building it if necessary.
require_build


# 根据第一个命令行参数决定执行哪一部分。
# Select which parts of the pipeline to execute based on the first argument.
case "$PHASE" in

    spinup)

        # Phase 1a:
        #   500 年海洋 spin-up
        #   500-year ocean spin-up
        #
        # Phase 1b:
        #   生成 10 个候选初始条件
        #   Generate 10 candidate initial conditions
        #
        # Phase 1c:
        #   选择暖状态 A 和冷状态 B
        #   Select warm state A and cold state B
        phase_spinup
        phase_initial_conditions
        phase_select
        ;;


    historical)

        # Phase 2:
        #   从 A/B 分别启动 control 和 historical 实验
        #   Run control and historical experiments from A/B
        #
        # Phase 3:
        #   自动执行分析
        #   Run analysis automatically
        phase_historical
        phase_analysis
        ;;


    all)

        # 完整执行 Phase 1。
        # Run the complete Phase 1.
        phase_spinup
        phase_initial_conditions
        phase_select

        # 随后执行 Phase 2 和 Phase 3。
        # Then run Phase 2 and Phase 3.
        phase_historical
        phase_analysis
        ;;


    *)

        # 无法识别的 phase 参数。
        # Unrecognized phase argument.
        #
        # 使用 exit code 2 表示命令行参数/用法错误。
        # Exit code 2 conventionally indicates invalid command-line usage.
        echo "unknown phase: $PHASE" >&2
        exit 2
        ;;
esac


# 所有请求的阶段执行结束。
# All requested phases have finished.
log "===== pipeline done (phase=$PHASE) ====="
