# =============================================================================
# 用户可配置变量
# =============================================================================
#
# 本区域集中放置日常构建和运行时最常修改的变量。
# 大多数变量使用 “?=” 设置默认值，因此可以在命令行中直接覆盖，例如：
#
#   make NLAT=96 NLEV=20
#   make MPI=1 NPRO=8
#   make run STEPS=64
#   make run RESTART=path/to/restart
#
# 通常情况下，使用者只需要阅读和修改这一部分；
# 后面的“构建实现”区域主要描述 Makefile 内部的派生变量和编译规则。


# -----------------------------------------------------------------------------
# 编译工具配置
# -----------------------------------------------------------------------------
#
# “?=” 表示仅在变量尚未被定义时赋值。
# 因此，用户可以在命令行或环境变量中覆盖默认编译器，例如：
#   make CC=gcc
#   make CC=clang
# 如果外部没有指定 CC，则默认使用系统的 cc 命令。
CC ?= cc

# 静态库归档工具。
# AR 通常指向 ar，用于将多个目标文件打包成 .a 静态库。
# 同样可以通过命令行覆盖，例如：
#   make AR=llvm-ar
AR ?= ar

# Python is used only by conversion/integration tests, not by the C runtime.
PYTHON ?= python3


# -----------------------------------------------------------------------------
# 模型网格配置
# -----------------------------------------------------------------------------

# 纬向网格参数，默认值为 64。
# 该值会通过 CPPFLAGS 传递给 C 预处理器，
# 最终在 C 代码中以 PLASIC_NLAT 宏的形式使用。
NLAT ?= 48

# 垂直层数，默认值为 16。
# 最终会定义预处理宏 PLASIC_NLEV。
NLEV ?= 13

# 处理器或并行分区数量，默认值为 1。
# 最终会定义预处理宏 PLASIC_NPRO。
NPRO ?= 2

# 海洋模式后端。当前支持：
#   slab - 原始单层混合层海洋（混合层深度海洋）；
#   3d   - 三维深海洋，与混合层耦合，默认启用。
# 三维海洋的编译开关为 PLASIC_DEEP_OCEAN。
OCEAN_MODEL ?= 3d

# 三维深海洋的垂向层数，默认 16 层，最终定义 PLASIC_DEEP_NLEV。
DEEP_OCEAN_LEVELS ?= 16


# -----------------------------------------------------------------------------
# MPI 配置
# -----------------------------------------------------------------------------
#
# MPI=1 selects the real MPI implementation and an MPI compiler wrapper.
# The default serial build uses the explicitly named mpi_serial_backend.
#
# 例如：
#   make MPI=1 NPRO=8
#   make MPI=1 MPI_CC=mpicc MPIRUN=mpirun NPRO=8
MPI ?= 1
MPI_CC ?= mpicc
MPIRUN ?= mpirun


# -----------------------------------------------------------------------------
# 构建产物目录配置
# -----------------------------------------------------------------------------

# 构建产物根目录。
# 默认指向项目根目录下的 build/（Makefile 位于 src/，因此使用 ../build），
# 避免构建产物污染源代码目录。
BUILD_ROOT ?= ../build

# 供从项目根目录执行的测试命令使用；将 BUILD_ROOT 规范化为绝对路径，
# 这样覆盖 BUILD_ROOT 时，测试程序仍会读取本次构建得到的可执行文件。
ABS_BUILD_ROOT := $(abspath $(BUILD_ROOT))

# 构建产物目录。
# 默认目录名会根据 NLAT、NLEV 和 MPI 配置自动生成。
# 例如默认 NLAT=64、NLEV=16、MPI=0 时：
#   ../build/T42_L16
#
# MPI=1、NPRO=4 时会附加 MPI 后缀，例如：
#   ../build/T42_L16_MPI_P4
#
# 三维海洋启用时再附加海洋层数后缀，例如 ../build/T42_L16_O16，这样修改
# DEEP_OCEAN_LEVELS 会使用独立目录，避免旧的 .o 与新的结构体布局混用。
#
# 也可以直接覆盖完整构建目录：
#   make BUILD_DIR=/tmp/plasic-build
OCEAN_BUILD_SUFFIX = $(if $(filter 3d,$(OCEAN_MODEL)),_O$(DEEP_OCEAN_LEVELS),)
BUILD_DIR ?= $(BUILD_ROOT)/T$(TRUNCATION)_L$(NLEV)$(MPI_BUILD_SUFFIX)$(OCEAN_BUILD_SUFFIX)


# -----------------------------------------------------------------------------
# 运行产物目录配置
# -----------------------------------------------------------------------------

# 运行产物根目录。
# 默认指向项目根目录下的 run/（Makefile 位于 src/，因此使用 ../run），
# 避免运行产物污染源代码目录。
RUN_ROOT ?= ../run

# 运行产物目录。
# 默认目录名会根据 NLAT、NLEV 和 MPI 配置自动生成
# 例如默认 NLAT=64、NLEV=16、MPI=0 时：
#   ../run/T42_L16
#
# 也可以直接覆盖完整运行目录：
#   make run RUN_DIR=/tmp/plasic-run
RUN_DIR ?= $(RUN_ROOT)/T$(TRUNCATION)_L$(NLEV)$(MPI_BUILD_SUFFIX)$(OCEAN_BUILD_SUFFIX)


# -----------------------------------------------------------------------------
# 标准运行配置
# -----------------------------------------------------------------------------

# 模型运行步数。
# 原先 run 规则固定使用 --steps 32；现在提取为变量，默认行为保持不变。
# 例如：
#   make run STEPS=128
STEPS ?= 100

# 每隔多少步向 $(RUN_DIR)/progress.txt 追加一条 JSON 进度记录。
# 默认每一步都写；长积分可调大该值以减少文件体积，例如：
#   make run STEPS=20000 PROGRESS_INTERVAL=100
PROGRESS_INTERVAL ?= 1

# RESTART 指定初始 restart 文件路径，可在命令行覆盖：
#   make run RESTART=path/to/restart
#
# 当前默认值保持为 ../status，与原 Makefile 行为一致。
# 如果希望让 Makefile 自动检查 $(RUN_DIR)/status，并在不存在时不传
# --restart，可以显式将 RESTART 设为空：
#   make run RESTART=
# 当 RESTART 为空时，按如下优先级确定输入 restart：
#   1. 运行目录 $(RUN_DIR)/status 已存在 → 自动续跑；
#   2. 否则不传 --restart，所有分辨率都从 NetCDF surface data 冷启动。

# RESTART ?= ../status
RESTART ?=

# Historical experiment defaults. `make historical` deliberately uses a
# separate T21/L10 MPI build and requires an explicit initial condition:
#
#   make historical INITIAL_CONDITION=data/earth_t21_l10.restart
#
# The model CLI no longer has --end-year or --reset-calendar: a warm start
# rebases the clock with --reset-nstep and takes an explicit --steps. The two
# convenience values below are derived here from the start/end years with
# proleptic-Gregorian day counts and the T21 step length (ntspd = 32).
HISTORICAL_NLAT ?= 32
HISTORICAL_NLEV ?= 10
HISTORICAL_NPRO ?= 2
HISTORICAL_START_YEAR ?= 1850
HISTORICAL_END_YEAR ?= 2014
HISTORICAL_NTSPD ?= 32
HISTORICAL_RESET_NSTEP ?= $(shell $(PYTHON) -c "y=$(HISTORICAL_START_YEAR); print((365*y+(y+3)//4-(y+99)//100+(y+399)//400)*$(HISTORICAL_NTSPD))")
HISTORICAL_STEPS ?= $(shell $(PYTHON) -c "d=lambda y:365*y+(y+3)//4-(y+99)//100+(y+399)//400; print((d($(HISTORICAL_END_YEAR)+1)-d($(HISTORICAL_START_YEAR)))*$(HISTORICAL_NTSPD))")
HISTORICAL_PROGRESS_INTERVAL ?= 320
INITIAL_CONDITION ?=
CO2_FORCING ?= data/forcing/CO2/mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_CMIP_UoM-CMIP-1-2-0_gn-15x360deg_000001-201412.csv
CO2_PPM ?=
HISTORICAL_OUTPUT_DIR ?= $(RUN_ROOT)/historical_T21_L10_MPI_P$(HISTORICAL_NPRO)_$(HISTORICAL_START_YEAR)-$(HISTORICAL_END_YEAR)


# -----------------------------------------------------------------------------
# 编译和链接参数
# -----------------------------------------------------------------------------

# CPPFLAGS 通常用于保存预处理宏定义（-D）和头文件搜索路径（-I）。
# 项目所需的宏定义和 -I. 会在后面的“构建实现”区域自动追加。
CPPFLAGS ?=

# 默认优化等级为 -O2。
# 如果外部已经定义 CFLAGS，则这里不会自动添加 -O2。
CFLAGS ?= -O2

# WARNFLAGS 独立于 CFLAGS，便于调用者整体替换警告策略。
#
# -Wall：
#   启用常用警告。
#
# -Wextra：
#   启用更多额外警告。
#
# -Wpedantic：
#   检查不符合所选 C 标准的扩展或写法。
#
# -Wconversion：
#   检查可能改变数值或精度的隐式类型转换。
#
# -Wshadow：
#   检查局部变量遮蔽外层变量、全局变量或参数的情况。
#
# -Wstrict-prototypes：
#   要求 C 函数声明提供明确的参数列表，
#   避免使用类似 “int func();” 的非严格原型。
WARNFLAGS ?= -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes

# 额外链接器参数。
LDFLAGS ?=

# 额外链接库。
# 项目所需的数学库 -lm 会在后面的“构建实现”区域自动追加。
LDLIBS ?=

# SHTns provides the complete scalar and vector spherical-harmonic transform.
# Use a prefix for a local installation, or let pkg-config/system paths find it:
#   make SHTNS_PREFIX=/path/to/shtns
SHTNS_PREFIX ?=
SHTNS_CFLAGS ?= $(if $(strip $(SHTNS_PREFIX)),-I$(SHTNS_PREFIX)/include,$(shell pkg-config --cflags shtns 2>/dev/null))
SHTNS_LIBS ?= $(if $(strip $(SHTNS_PREFIX)),-L$(SHTNS_PREFIX)/lib -lshtns,$(shell pkg-config --libs shtns 2>/dev/null))
FFTW_CFLAGS ?= $(shell pkg-config --cflags fftw3 2>/dev/null)
FFTW_LIBS ?= $(shell pkg-config --libs fftw3 2>/dev/null)

# NetCDF-C is the only supported runtime surface-data reader. Override these
# variables for a non-system installation; no text/SRA fallback is compiled.
NETCDF_CFLAGS ?= $(shell pkg-config --cflags netcdf 2>/dev/null)
NETCDF_LIBS ?= $(shell pkg-config --libs netcdf 2>/dev/null)

# LAPACK provides the matrix-inversion implementation used by
# math_tools/linear_algebra.c. macOS ships LAPACK in Accelerate; other Unix
# platforms can use pkg-config or fall back to the conventional -llapack.
HOST_OS := $(shell uname -s)
ifeq ($(HOST_OS),Darwin)
LAPACK_LIBS ?= -framework Accelerate
else
LAPACK_LIBS ?= $(shell pkg-config --libs lapack 2>/dev/null)
ifeq ($(strip $(LAPACK_LIBS)),)
LAPACK_LIBS = -llapack
endif
endif

ifeq ($(strip $(SHTNS_LIBS)),)
SHTNS_LIBS = -lshtns
endif
ifeq ($(strip $(FFTW_LIBS)),)
FFTW_LIBS = -lfftw3
endif
ifeq ($(strip $(NETCDF_LIBS)),)
NETCDF_LIBS = -lnetcdf
endif

# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# 构建实现
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# 从这里开始主要是 Makefile 的内部实现：
#   1. 根据用户配置计算派生变量；
#   2. 定义源文件和目标文件；
#   3. 定义编译、链接、运行和清理规则。
#
# 一般使用者通常不需要修改下面的内容。


# -----------------------------------------------------------------------------
# MPI 派生配置
# -----------------------------------------------------------------------------

ifeq ($(MPI),1)
CC := $(MPI_CC)
MPI_SOURCE = mpi/mpi_native_backend.c
MPI_BUILD_SUFFIX = _MPI_P$(NPRO)
RUNNER = $(MPIRUN) -np $(NPRO)
else
MPI_SOURCE =
MPI_BUILD_SUFFIX =
RUNNER =
endif


# -----------------------------------------------------------------------------
# 网格派生参数
# -----------------------------------------------------------------------------
#
# 谱截断分辨率的计算公式为：
#   (NLAT * 2 - 1) / 3
#
# shell 的 expr 执行整数运算，因此结果会进行整数除法。
# 例如默认 NLAT=64 时：
#   (64 * 2 - 1) / 3 = 127 / 3 = 42
#
# 因此默认 NLAT=64、NLEV=16、MPI=0 时：
#   BUILD_DIR = ../build/T42_L16
#   RUN_DIR   = ../run/T42_L16
#
# 反斜杠用于转义圆括号，避免圆括号被 shell 解释为子 shell 语法。
TRUNCATION = $(shell expr \( $(NLAT) \* 2 - 1 \) / 3)


# -----------------------------------------------------------------------------
# Restart 派生配置
# -----------------------------------------------------------------------------

# 根据 RESTART 的值确定实际传入模型的 restart 文件。
#
# $(strip ...) 去除首尾空白后判断是否为空。
# $(wildcard ...) 在 Makefile 解析阶段检查文件是否存在；
# 若运行目录已有 status 则自动作为续跑输入。
ifeq ($(strip $(RESTART)),)
  RESTART_INPUT = $(wildcard $(RUN_DIR)/status)
# wildcard 按照通配符模式查找实际存在的文件；存在就返回文件名，不存在就返回空字符串。
else
  RESTART_INPUT = $(RESTART)
endif

# 将非空的 RESTART_INPUT 转换为 --restart 标志；为空则不传该参数，
# 模型据此执行冷启动。
ifneq ($(strip $(RESTART_INPUT)),)
  RESTART_FLAG = --restart $(RESTART_INPUT)
else
  RESTART_FLAG =
endif

# Normal runs use the external forcing by default. Setting CO2_PPM switches
# only `make run` to a spatially and temporally fixed concentration.
ifneq ($(strip $(CO2_PPM)),)
  RUN_CO2_FLAG = --co2-ppm $(CO2_PPM)
else
  RUN_CO2_FLAG = --co2-forcing "$(CO2_FORCING)"
endif


# -----------------------------------------------------------------------------
# C 预处理器参数
# -----------------------------------------------------------------------------
#
# CPPFLAGS 通常用于保存：
#
#   1. 预处理宏定义（-D）
#   2. 头文件搜索路径（-I）
#
# 使用 “+=” 而不是 “=” 可以在默认情况下追加项目所需参数，
# 而不是在这里重新定义整组 CPPFLAGS。

# 将 Makefile 中的网格和运行配置转换为 C 预处理宏。
#
# 编译后，C 源代码中可以直接使用：
#
#   PLASIC_NLAT
#   PLASIC_NLEV
#   PLASIC_NPRO
CPPFLAGS += -DPLASIC_NLAT=$(NLAT) -DPLASIC_NLEV=$(NLEV) -DPLASIC_NPRO=$(NPRO)

# 定义三维深海洋的垂向层数和启用开关。
CPPFLAGS += -DPLASIC_DEEP_NLEV=$(DEEP_OCEAN_LEVELS)
ifeq ($(OCEAN_MODEL),3d)
CPPFLAGS += -DPLASIC_DEEP_OCEAN=1
else
CPPFLAGS += -DPLASIC_DEEP_OCEAN=0
endif

# 指示当前构建是否启用真实 MPI 实现。
CPPFLAGS += -DPLASIC_USE_MPI=$(MPI)

# 将当前目录加入头文件搜索路径。
# 这样源文件中的：
#   #include "xxx.h"
#
# 可以从当前工程目录中查找相应头文件。
# SHTns/FFTW 只提供数学变换实现，不属于 PlaSiC 模式核心源码。
CPPFLAGS += -I.
CPPFLAGS += $(SHTNS_CFLAGS) $(FFTW_CFLAGS) $(NETCDF_CFLAGS)


# -----------------------------------------------------------------------------
# C 编译参数
# -----------------------------------------------------------------------------

# 要求编译器按照 C11 标准编译。
# -fno-fast-math
# -f：用于控制某个编译功能/优化选项
# no-：表示禁用
# 禁用可能破坏严格 IEEE 浮点语义的激进优化。
# 这对于需要可重复数值结果的气候或物理模型尤其重要。
CFLAGS += -std=c11 -fno-fast-math


# -----------------------------------------------------------------------------
# 链接参数
# -----------------------------------------------------------------------------

# 链接数学库 libm。
# 许多 Unix 类系统要求显式使用 -lm，才能解析：
#   sin
#   cos
#   sqrt
# 等数学函数。
LDLIBS += $(SHTNS_LIBS) $(FFTW_LIBS) $(NETCDF_LIBS) $(LAPACK_LIBS) -lm


# -----------------------------------------------------------------------------
# 内核源文件列表
# -----------------------------------------------------------------------------
#
# 第三方 vendored 源码（cJSON）位于仓库根目录，即本 Makefile 所在 src/ 的上一级。
# 通过 VPATH 让 KERNEL_SOURCES 中的 third_party/... 路径能被 make 解析，
# 同时目标文件仍然生成在 $(BUILD_DIR)/third_party/... 下。
VPATH = ..
#
# KERNEL_SOURCES 保存模型核心计算、物理过程、状态管理、输入输出和辅助模块对应的 C 源文件。
#
# 每行末尾的反斜杠表示变量定义在下一行继续。
#
# 后续这些源文件会被：
#
#   1. 分别编译为目标文件 .o
#   2. 打包为静态库 libplasic.a
#
KERNEL_SOURCES = \
	tools/calmod_kernels.c \
	tools/co2_forcing_table.c \
	physics/fluxmod_kernels.c \
	math_tools/gaussmod_kernels.c \
	math_tools/linear_algebra.c \
	sea/icemod_kernels.c \
	land/landmod_kernels.c \
	math_tools/sht_kernels.c \
	tools/miscmod_kernels.c \
	tools/monthly_netcdf_writer.c \
	mpi/plasic_mpi.c \
	mpi/mpi_serial_backend.c \
	$(MPI_SOURCE) \
	sea/slab_ocean.c \
	sea/3docean.c \
	tools/p_earth_kernels.c \
	third_party/cJSON/cJSON.c \
	plasic_config.c \
	dynamical_core/calgp.c \
	dynamical_core/plasic_kernels.c \
	tools/plasic_stream.c \
	state/plasicmod_state.c \
	physics/radmod_kernels.c \
	physics/rainmod_kernels.c \
	tools/restartmod_io.c \
	sea/seamod_kernels.c \
	tools/surface_data_io.c \
	tools/surfmod_kernels.c \
	physics/vdiff_kernel.c


# -----------------------------------------------------------------------------
# 运行时和主程序源文件
# -----------------------------------------------------------------------------
#
# 这些源文件不会被打包进核心静态库，而是在最终链接阶段
# 直接与 libplasic.a 一起链接生成可执行程序。
#
# 从命名上看，它们分别负责：
#
#   runtime/：
#     初始化、时间积分与在线输出；plasic_runtime.c 只保留公共入口。
#
#   physics/physics_*.c：
#     分别管理物理状态初始化、过程调度和输入输出。
#
#   app/：
#     plasic_standalone.c 只包含 main，plasic_cli.c 负责参数策略。
RUNTIME_SOURCES = \
	runtime/plasic_runtime.c \
	runtime/co2_forcing.c \
	runtime/monthly_netcdf.c \
	runtime/runtime_initialization.c \
	runtime/runtime_dynamics.c \
	runtime/runtime_output.c \
	physics/plasic_physics.c \
	physics/physics_initialization.c \
	physics/physics_processes.c \
	physics/physics_io.c \
	app/plasic_cli.c \
	app/plasic_standalone.c


# -----------------------------------------------------------------------------
# 根据源文件列表生成目标文件路径
# -----------------------------------------------------------------------------
#
# 模式替换语法：
#
#   $(变量:旧模式=新模式)
#
# 会将每个匹配 %.c 的文件名转换为：
#
#   $(BUILD_DIR)/%.o
#
# 例如默认配置下：
#
#   calmod_kernels.c
#
# 会被转换为：
#
#   ../build/T42_L16/calmod_kernels.o
KERNEL_OBJECTS = $(KERNEL_SOURCES:%.c=$(BUILD_DIR)/%.o)
RUNTIME_OBJECTS = $(RUNTIME_SOURCES:%.c=$(BUILD_DIR)/%.o)


# -----------------------------------------------------------------------------
# 最终构建产物
# -----------------------------------------------------------------------------

# 核心静态库路径。
LIBRARY = $(BUILD_DIR)/libplasic.a

# 最终可执行程序路径。
PROGRAM = $(BUILD_DIR)/plasic.x


# -----------------------------------------------------------------------------
# 伪目标声明
# -----------------------------------------------------------------------------
#
# .PHONY 表示这些目标名称不对应真实文件。
#
# 即使源码目录中恰好存在名为 all、clean、run 等文件，
# make 仍然会执行对应规则，而不会错误地认为目标已经是最新状态。
.PHONY: all clean help strict run historical historical-run test test-data test-runtime


# ----------------------------------------------------------------------------
# Command overview
# ----------------------------------------------------------------------------
#
# Keep the most common entry points discoverable without requiring users to
# read the full Makefile first.
help:
	@echo "PlaSiC build targets:"
	@echo "  make                  Build the model"
	@echo "  make run              Build and run a short integration"
	@echo "  make test             Run data and runtime tests"
	@echo "  make strict           Build with warnings treated as errors"
	@echo "  make clean            Remove build and run products"
	@echo ""
	@echo "Common options:"
	@echo "  MPI=0                 Build without MPI"
	@echo "  NLAT=64 NLEV=16      Set horizontal and vertical resolution"
	@echo "  NPRO=4                Set the number of MPI processes"
	@echo "  STEPS=100             Set the number of run steps"
	@echo "  SHTNS_PREFIX=/path    Use a custom SHTns installation"


# -----------------------------------------------------------------------------
# 默认构建目标
# -----------------------------------------------------------------------------
#
# Makefile 中第一个普通目标通常是默认目标。
# 执行：
#   make
# 等价于：
#   make all
# all 本身没有命令，只要求最终程序 $(PROGRAM) 构建完成。
all: $(PROGRAM)


# -----------------------------------------------------------------------------
# 链接最终可执行程序
# -----------------------------------------------------------------------------
#
# 最终程序依赖：
#   1. 运行时目标文件 $(RUNTIME_OBJECTS)
#   2. 核心静态库 $(LIBRARY)
#
# 只要其中任意依赖不存在或比最终程序更新，
# make 就会重新执行下面的链接命令。

# 当目标程序 $(PROGRAM) 需要生成时，把若干 .o 目标文件和库链接起来，生成最终的可执行文件。
# 目标: 依赖1 依赖2 ...
# 	生成目标时要执行的命令
$(PROGRAM): $(RUNTIME_OBJECTS) $(LIBRARY)
	# $@ 是 Makefile 的一个自动变量，表示：当前规则的目标文件名。
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(RUNTIME_OBJECTS) $(LIBRARY) $(LDLIBS)


# -----------------------------------------------------------------------------
# 构建核心静态库
# -----------------------------------------------------------------------------
#
# libplasic.a 依赖所有核心目标文件。
# 自动变量含义：
#   $@：
#     当前规则的目标，即 $(LIBRARY)。
#   $^：
#     当前规则的全部依赖文件，
#     即 $(KERNEL_OBJECTS)，并自动去除重复项。
#
# ar 参数含义：
#   r：
#     将目标文件插入静态库；同名成员存在时进行替换。
#   c：
#     如果静态库不存在，则创建它。
#   s：
#     为静态库创建或更新符号索引，便于链接器查找符号。
$(LIBRARY): $(KERNEL_OBJECTS)
	$(AR) rcs $@ $^


# -----------------------------------------------------------------------------
# 通用 C 源文件编译规则
# -----------------------------------------------------------------------------
#
# 这是一个模式规则。
# 它描述了如何将当前目录中的任意：
#   %.c
# 编译为构建目录中的：
#   $(BUILD_DIR)/%.o
#
# 自动变量含义：
#   $@：
#     当前目标文件路径，例如 ../build/T42_L16/math_tools/sht_kernels.o。
#   $<：
#     第一个依赖文件，即对应的 C 源文件。
#
# $(dir $@)：
#   提取目标文件所在目录，并保留末尾斜杠。
#
# @mkdir：
#   命令前的 @ 会阻止 make 输出该命令本身，
#   但 mkdir 的错误信息仍然会正常显示。
#
# mkdir -p：
#   在编译前递归创建构建目录；
#   如果目录已经存在，也不会报错。
#
# -MMD（Makefile make dependency）：
#   编译时自动生成用户头文件依赖文件 .d。
#   系统头文件通常不会写入该依赖文件。
#
# -MP (Maekfile phony)：
#   为依赖文件中涉及的头文件生成空目标，
#   避免某个头文件被删除后 make 因“没有规则可生成该文件”而直接失败。

# -MMD
# │
# └─ “帮我记录这个 .o 到底依赖哪些用户头文件”
#    → 生成 .d

# -MP
# │
# └─ “顺便给这些头文件放一个空规则”
#    → 防止旧 .d 引用已经删除的头文件时
#       make 直接报 No rule to make target

# -c：
#   仅编译，不执行链接。
#
# -o $@：
#   将输出写入目标路径。
$(BUILD_DIR)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARNFLAGS) -MMD -MP -c $< -o $@

# -----------------------------------------------------------------------------
# 严格编译检查
# -----------------------------------------------------------------------------
# 该目标专门为 strict 及其依赖追加 -Werror。
strict: CFLAGS += -Werror
strict: clean all


# -----------------------------------------------------------------------------
# 标准运行目标
# -----------------------------------------------------------------------------
#
# run 依赖 all，因此程序不存在或源码发生变化时会先自动构建。
#
# 运行产物统一写入 $(RUN_DIR)，
# 运行参数：
#   $(RESTART_FLAG)
#     指定初始 restart 文件。由 RESTART 变量控制：
#       make run RESTART=path/to/restart
#     当 RESTART 为空时，若 $(RUN_DIR)/status 已存在则自动续跑，
#     否则模型执行冷启动。
#
#   --steps $(STEPS)
#     执行 $(STEPS) 个模型步，默认值为 32。
#
#   --progress $(RUN_DIR)/progress.txt
#     持续追加 JSON 进度记录，便于在运行过程中监控模型状态。
#
#   --output $(RUN_DIR)/status
#     将最终状态或重启状态写入 status 文件。
#
# 命令被拆成多行，仅用于提高可读性；
# 行末反斜杠使 shell 将多行视为一条完整命令。
run: all
	@mkdir -p $(RUN_DIR)
	$(RUNNER) $(PROGRAM) $(RESTART_FLAG) --steps $(STEPS) \
		$(RUN_CO2_FLAG) \
		--progress $(RUN_DIR)/progress.txt \
		--progress-interval $(PROGRESS_INTERVAL) \
		--output $(RUN_DIR)/status

# Build and run the standard transient historical experiment. A recursive make
# is used so the target has independent T21/L10 MPI defaults without changing
# the normal `make` and `make run` configurations.
historical:
	@test -n "$(INITIAL_CONDITION)" || { \
		echo "INITIAL_CONDITION is required." >&2; \
		echo "Example: make historical INITIAL_CONDITION=data/earth_t21_l10.restart" >&2; \
		exit 2; \
	}
	@test -f "$(INITIAL_CONDITION)" || { \
		echo "Initial condition not found: $(INITIAL_CONDITION)" >&2; \
		exit 2; \
	}
	@test -f "$(CO2_FORCING)" || { \
		echo "CO2 forcing not found: $(CO2_FORCING)" >&2; \
		exit 2; \
	}
	$(MAKE) historical-run MPI=1 \
		NLAT=$(HISTORICAL_NLAT) NLEV=$(HISTORICAL_NLEV) \
		NPRO=$(HISTORICAL_NPRO) \
		INITIAL_CONDITION="$(INITIAL_CONDITION)" \
		CO2_FORCING="$(CO2_FORCING)" \
		HISTORICAL_START_YEAR=$(HISTORICAL_START_YEAR) \
		HISTORICAL_END_YEAR=$(HISTORICAL_END_YEAR) \
		HISTORICAL_PROGRESS_INTERVAL=$(HISTORICAL_PROGRESS_INTERVAL) \
		HISTORICAL_OUTPUT_DIR="$(HISTORICAL_OUTPUT_DIR)"

historical-run: all
	@mkdir -p "$(HISTORICAL_OUTPUT_DIR)/monthly"
	$(RUNNER) $(PROGRAM) \
		--restart "$(INITIAL_CONDITION)" \
		--co2-forcing "$(CO2_FORCING)" \
		--reset-nstep $(HISTORICAL_RESET_NSTEP) \
		--steps $(HISTORICAL_STEPS) \
		--monthly-output "$(HISTORICAL_OUTPUT_DIR)/monthly" \
		--progress "$(HISTORICAL_OUTPUT_DIR)/progress.txt" \
		--progress-interval $(HISTORICAL_PROGRESS_INTERVAL) \
		--output "$(HISTORICAL_OUTPUT_DIR)/final.restart"

# Surface-data integrity: validate the complete 45-file archive, grid-specific
# filenames, dimensions, time axes, and finite field values.
test-data:
	cd .. && $(PYTHON) -m unittest discover -s tests -p 'test_*.py'

# End-to-end acceptance matrix. Separate BUILD_DIR values avoid collisions
# between the grid/layer configurations used for positive and negative tests.
test-runtime:
	$(MAKE) all NLAT=32 NLEV=10 BUILD_DIR=$(BUILD_ROOT)/tests/T21_L10
	$(MAKE) all NLAT=32 NLEV=12 BUILD_DIR=$(BUILD_ROOT)/tests/T21_L12
	$(MAKE) all NLAT=48 NLEV=12 BUILD_DIR=$(BUILD_ROOT)/tests/T31_L12
	$(MAKE) all NLAT=64 NLEV=16 BUILD_DIR=$(BUILD_ROOT)/tests/T42_L16
	cd .. && $(PYTHON) tests/test_runtime_integration.py \
		--t21 "$(ABS_BUILD_ROOT)/tests/T21_L10/plasic.x" \
		--t21-l12 "$(ABS_BUILD_ROOT)/tests/T21_L12/plasic.x" \
		--t31 "$(ABS_BUILD_ROOT)/tests/T31_L12/plasic.x" \
		--t42 "$(ABS_BUILD_ROOT)/tests/T42_L16/plasic.x" \
		--surface src/data/surface \
		--earth-t21 src/data/earth_t21_l10.restart \
		--earth-t42 src/data/earth_t42_l16.restart

test: test-data test-runtime


# -----------------------------------------------------------------------------
# 清理构建产物
# -----------------------------------------------------------------------------
# 删除整个构建产物根目录（BUILD_ROOT）和运行产物根目录（RUN_ROOT），
# 而不只删除当前 BUILD_DIR / RUN_DIR。
#
# 这意味着不同 NLAT、NLEV 配置生成的全部构建与运行产物都会被删除。
#
# rm -rf：
#
#   -r：
#     递归删除目录及其内容。
#
#   -f：
#     忽略不存在的文件，并且不进行交互式确认。
clean:
	rm -rf $(BUILD_ROOT) $(RUN_ROOT)


# -----------------------------------------------------------------------------
# 自动包含头文件依赖
# -----------------------------------------------------------------------------
# 编译规则中的 -MMD -MP 会为每个目标文件生成对应的 .d 文件。
#
# 例如默认配置下：
#   ../build/T42_L16/plasic_runtime.o
#
# 会对应：
#   ../build/T42_L16/plasic_runtime.d
#
# 这些 .d 文件记录源文件所包含的项目头文件。
# 当头文件发生变化时，make 可以据此重新编译受影响的源文件。
#
# “-include” 前面的减号表示：
#   如果 .d 文件尚不存在，不将其视为错误。
# 这对于第一次构建非常重要，因为第一次编译之前还没有生成 .d 文件。
# 两组目标文件的后缀替换：
#   $(KERNEL_OBJECTS:.o=.d)
#   $(RUNTIME_OBJECTS:.o=.d)
-include $(KERNEL_OBJECTS:.o=.d) $(RUNTIME_OBJECTS:.o=.d)
