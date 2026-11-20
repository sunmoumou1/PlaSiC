#                     ┌─────────────────┐
#                     │    用户操作 UI   │
#                     │ ParameterForm   │
#                     │ buttons/combo   │
#                     └────────┬────────┘
#                              │
#                              ▼
#                     ┌─────────────────┐
#                     │   MainWindow    │
#                     │  orchestration  │
#                     └────────┬────────┘
#                              │
#             ┌────────────────┼────────────────┐
#             │                │                │
#             ▼                ▼                ▼
#    ┌────────────────┐ ┌──────────────┐ ┌───────────────┐
#    │ Config objects │ │ProcessController│ │ RunManager    │
#    │ BuildConfig    │ │ build / run   │ │prepare run dir │
#    │ RunConfig      │ │ subprocess    │ │config/log/etc  │
#    └────────┬───────┘ └───────┬──────┘ └───────────────┘
#             │                  │
#             └──────────┬───────┘
#                        ▼
#                ┌───────────────┐
#                │  C model core │
#                │ PlaSiC binary │
#                └───────┬───────┘
#                        │
#           ┌────────────┴─────────────┐
#           ▼                          ▼
#     progress file              frames.pstream
#           │                          │
#           ▼                          ▼
#    ProgressTail                PStreamReader
#           │                          │
#           └────────────┬─────────────┘
#                        ▼
#                   FrameBuffer
#                        │
#                        ▼
#                  LiveFieldView



# =============================================================================
# PlaSiC Studio 主窗口模块
#
# Review 导读：
# 1. 本文件负责桌面 GUI 的主窗口与交互编排，不直接实现数值模式本身。
# 2. 主要职责包括：参数收集、模型构建/启动、运行状态管理、输出轮询、实时可视化、预设读写。
# 3. 本注释版仅新增中文 `#` 注释；所有原始代码行均保持不变。
# =============================================================================
from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
from typing import Any

# --- Qt / PySide6 -------------------------------------------------------------
# QTimer 驱动周期轮询；Qt 提供枚举；QIcon 用于窗口图标；
# QTextCursor 用于把日志光标移动到文本末尾；其余控件组成主界面。
from PySide6.QtCore import QTimer, Qt
from PySide6.QtGui import QIcon, QTextCursor
from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QSplitter,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

# --- PlaSiC 应用内部模块 ------------------------------------------------------
# config：构建、运行、可视化配置模型及参数/变量目录。
from plasic_app.config import (
    BuildConfig,
    ExperimentConfig,
    ParameterCatalog,
    RunConfig,
    VariableCatalog,
    VisualizationConfig,
)

# controller.commands：把配置转换为“构建命令”和“运行命令”的调用描述。
from plasic_app.controller.commands import (
    BuildInvocation,
    build_invocation,
    run_invocation,
)

# ProcessController：负责真正启动/停止外部进程，并通过 Qt signal 回传状态/日志。
# protocol：增量读取进度文件和二进制 frame stream，并维护帧缓冲。
# RunManager：创建一次实验所需的运行目录、配置文件和输出路径。
# UI/visualization 模块：参数表单、主题、自定义分段控件和实时场视图。
from plasic_app.controller.process_controller import ProcessController
from plasic_app.paths import discover_c_root
from plasic_app.protocol import FrameBuffer, ProgressTail, PStreamReader
from plasic_app.runs import PreparedRun, RunManager
from plasic_app.ui.branding import logo_pixmap
from plasic_app.ui.fonts import fixed_width_font
from plasic_app.ui.parameter_form import ParameterForm
from plasic_app.ui.theme import STUDIO_STYLESHEET
from plasic_app.ui.widgets import SegmentedControl
from plasic_app.visualization.live_view import LiveFieldView


# =============================================================================
# MainWindow：整个 PlaSiC Studio 的顶层窗口与流程协调器。
#
# 可以把该类理解为“UI + orchestration”层：
# - UI 事件从这里进入；
# - 配置对象在这里组装和校验；
# - 构建/运行请求交给 ProcessController；
# - 模型输出再由这里轮询并推送到 LiveFieldView。
# =============================================================================
class MainWindow(QMainWindow):
    # 初始化窗口、运行期状态对象，并在最后启动固定周期的输出轮询定时器。
    def __init__(self) -> None:
        super().__init__()
        # 基础窗口属性：标题、默认尺寸以及允许缩放到的最小尺寸。
        self.setWindowTitle("PlaSiC Studio")
        self.resize(1520, 940)
        self.setMinimumSize(1180, 560)

        # 发现 C 数值核心目录，并加载参数/变量元数据目录。
        # 这些 catalog 后续既驱动表单生成，也用于变量显示名称、色标等 UI 配置。
        self.c_root = discover_c_root()
        self.parameter_catalog = ParameterCatalog.load()
        self.variable_catalog = VariableCatalog.load()
        # 运行管理器负责创建实验目录；进程控制器负责外部 build/run 子进程。
        self.run_manager = RunManager()
        self.controller = ProcessController(self)

        # current_build：最近一次构建命令；current_run：当前已准备好的实验。
        # 两者使用 Optional 状态是因为应用启动时还没有任何构建或运行。
        self.current_build: BuildInvocation | None = None
        self.current_run: PreparedRun | None = None

        # progress_tail / frame_reader 分别用于增量读取进度记录和二进制场数据。
        # frame_buffer 是内存中的有限历史缓存，供变量/层级切换时快速重绘。
        self.progress_tail: ProgressTail | None = None
        self.frame_reader: PStreamReader | None = None
        self.frame_buffer = FrameBuffer()
        self._updating_time_controls = False
        self._follow_latest_time = True
        # 当前帧流的能力信息：文件头的海洋层数，以及已经出现过的变量名。
        # 目录中列出但当前运行无法产生的变量（例如未启用三维海洋时的
        # deep_ocean_temperature）会据此禁用并在状态栏说明原因。
        self._stream_ocean_levels: int | None = None
        self._stream_variables: set[str] = set()
        # 数据来源：True 表示通过 Replay Stream 打开的已保存帧流。
        self._replay_active = False
        # 数据来源是否已经读尽（回放读完或本次运行结束），用于区分
        # “还没有输出”和“确实没有该变量”。
        self._stream_complete = False
        # 内部运行状态（idle/building/running/stopping），与状态标签文案解耦。
        self._state = "idle"
        # pending_run=True 表示用户点了 Start Run，但可执行文件不存在，
        # 因而先异步 build；build 成功后需要自动续接真正的 run。
        self.pending_run = False

        # UI 构造与 signal-slot 连接分开，便于 review 控件布局和业务事件绑定。
        self._build_ui()
        self._connect_signals()

        # 每 150 ms 轮询一次增量输出。QTimer 在 GUI 线程中触发，
        # 因此 _poll_outputs 必须避免一次读取过多数据阻塞界面。
        self.poll_timer = QTimer(self)
        self.poll_timer.setInterval(150)
        self.poll_timer.timeout.connect(self._poll_outputs)
        self.poll_timer.start()
        self.playback_timer = QTimer(self)
        self.playback_timer.setInterval(500)
        self.playback_timer.timeout.connect(self._advance_playback)
        self._set_state("idle")

    # -------------------------------------------------------------------------
    # 构建整个窗口的可视层级。这里只创建/布局控件，不启动模型业务流程。
    # -------------------------------------------------------------------------
    def _build_ui(self) -> None:
        # 根 QWidget + 垂直布局：顶部 header，下面是左右分栏 splitter。
        root = QWidget()
        root.setObjectName("Root")
        root.setMinimumHeight(760)
        root_layout = QVBoxLayout(root)
        root_layout.setContentsMargins(16, 14, 16, 16)
        root_layout.setSpacing(10)
        self.app_scroll = QScrollArea()
        self.app_scroll.setObjectName("AppScroll")
        self.app_scroll.setFrameShape(QFrame.Shape.NoFrame)
        self.app_scroll.setWidgetResizable(True)
        self.app_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.app_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.app_scroll.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.app_scroll.setWidget(root)
        self.app_scroll.verticalScrollBar().setSingleStep(48)
        self.scroll_content = root
        self.setCentralWidget(self.app_scroll)

        # ===== 顶部 Header：Logo、标题、Help 按钮、运行状态 pill =====
        header = QFrame()
        header.setObjectName("HeaderBar")
        header_layout = QHBoxLayout(header)
        header_layout.setContentsMargins(14, 10, 14, 10)
        header_layout.setSpacing(11)

        # 从包资源加载 PlaSiC Logo。
        # 窗口图标与 Header 左侧的小图标共用这一张位图。
        logo = logo_pixmap()
        self.setWindowIcon(QIcon(logo))
        app_logo = QLabel()
        app_logo.setObjectName("AppLogo")
        app_logo.setAlignment(Qt.AlignmentFlag.AlignCenter)
        # 统一缩放 Logo 到 40×40（逻辑像素），并保持原始宽高比与平滑插值。
        # 按屏幕的 devicePixelRatio 生成更高分辨率的位图并标注像素比，
        # 否则 Retina 屏会把 40 px 位图再放大 2 倍，Logo 显得发虚。
        dpr = self.devicePixelRatioF()
        header_logo = logo.scaled(
            round(40 * dpr),
            round(40 * dpr),
            Qt.AspectRatioMode.KeepAspectRatio,
            Qt.TransformationMode.SmoothTransformation,
        )
        header_logo.setDevicePixelRatio(dpr)
        app_logo.setPixmap(header_logo)
        header_layout.addWidget(app_logo)
        # 标题区用一个内嵌的垂直布局承载主标题和副标题。
        title_layout = QVBoxLayout()
        title_layout.setSpacing(0)
        title = QLabel("PlaSiC Studio")
        title.setObjectName("AppTitle")
        subtitle = QLabel("CLIMATE SYSTEMS LAB")
        subtitle.setObjectName("AppSubtitle")
        title_layout.addWidget(title)
        title_layout.addWidget(subtitle)
        header_layout.addLayout(title_layout)
        # stretch 把 Help 和状态标签推到 Header 的右侧。
        header_layout.addStretch(1)
        self.configuration_toggle = QPushButton("Hide Configuration")
        self.configuration_toggle.setCheckable(True)
        self.configuration_toggle.setChecked(True)
        self.configuration_toggle.setToolTip(
            "Show or hide the experiment configuration sidebar"
        )
        header_layout.addWidget(self.configuration_toggle)
        self.help_button = QPushButton("Help")
        self.help_button.setObjectName("HelpButton")
        header_layout.addWidget(self.help_button)
        self.state_label = QLabel()
        self.state_label.setObjectName("StatePill")
        header_layout.addWidget(self.state_label)
        root_layout.addWidget(header)

        # 主体使用水平 QSplitter，用户可以调整左右面板宽度；
        # setChildrenCollapsible(False) 防止任一面板被拖到完全折叠。
        self.main_splitter = QSplitter(Qt.Orientation.Horizontal)
        self.main_splitter.setChildrenCollapsible(False)
        root_layout.addWidget(self.main_splitter, 1)

        # ===== 左侧 Card：实验参数、快捷时长、Build/Run 控制、进度和日志 =====
        self.configuration_panel = QFrame()
        self.configuration_panel.setObjectName("Card")
        left_layout = QVBoxLayout(self.configuration_panel)
        left_layout.setContentsMargins(15, 14, 15, 15)
        left_layout.setSpacing(9)
        setup_eyebrow = QLabel("Configuration")
        setup_eyebrow.setProperty("role", "eyebrow")
        left_layout.addWidget(setup_eyebrow)
        settings_title = QLabel("Experiment Setup")
        settings_title.setProperty("role", "section")
        left_layout.addWidget(settings_title)
        # ParameterForm 根据 parameter_catalog 动态生成所有配置项。
        self.form = ParameterForm(self.parameter_catalog)
        left_layout.addWidget(self.form, 1)

        # 预设操作：恢复默认值 / 从 JSON 加载 / 保存为 JSON。
        preset_buttons = QHBoxLayout()
        preset_buttons.setSpacing(8)
        self.defaults_button = QPushButton("Defaults")
        self.load_button = QPushButton("Load Preset")
        self.save_button = QPushButton("Save Preset")
        preset_buttons.addWidget(self.defaults_button)
        preset_buttons.addWidget(self.load_button)
        preset_buttons.addWidget(self.save_button)
        left_layout.addLayout(preset_buttons)

        # Quick Start 用“模拟天数”这种更直观的概念设置 integration steps。
        quick_title = QLabel("Quick Start (length in days)")
        quick_title.setProperty("role", "section")
        left_layout.addWidget(quick_title)
        quick_buttons = QHBoxLayout()
        quick_buttons.setSpacing(6)
        self.quick_1day = QPushButton("1 Day")
        self.quick_10days = QPushButton("10 Days")
        self.quick_30days = QPushButton("30 Days")
        self.quick_1year = QPushButton("1 Year")
        # 四个快捷按钮统一打上 role=quick，样式由 stylesheet 根据属性选择器控制。
        for button in (
            self.quick_1day,
            self.quick_10days,
            self.quick_30days,
            self.quick_1year,
        ):
            button.setProperty("role", "quick")
            quick_buttons.addWidget(button)
        left_layout.addLayout(quick_buttons)
        # 该提示会根据 resolution 和 steps 动态更新“步数 ≈ 模拟天数”。
        self.run_hint = QLabel("Current run length: 32 steps \u2248 1.0 day")
        self.run_hint.setProperty("role", "secondary")
        self.run_hint.setWordWrap(True)
        left_layout.addWidget(self.run_hint)

        # ===== 运行控制区：Build、Start Run、Stop =====
        run_title = QLabel("Run Control")
        run_title.setProperty("role", "section")
        left_layout.addWidget(run_title)
        control_buttons = QHBoxLayout()
        control_buttons.setSpacing(8)
        self.build_button = QPushButton("Build Model")
        self.run_button = QPushButton("Start Run")
        # role 属性只影响样式语义：primary 强调主操作，danger 强调停止操作。
        self.run_button.setProperty("role", "primary")
        self.stop_button = QPushButton("Stop")
        self.stop_button.setProperty("role", "danger")
        control_buttons.addWidget(self.build_button)
        control_buttons.addWidget(self.run_button)
        control_buttons.addWidget(self.stop_button)
        left_layout.addLayout(control_buttons)

        # 进度条内部范围设为 0..1000，而不是 0..100；
        # 这样 completed/total 映射时可以保留 0.1% 级别的显示分辨率。
        self.progress = QProgressBar()
        self.progress.setRange(0, 1000)
        self.progress.setValue(0)
        self.progress_detail = QLabel("No run started")
        self.progress_detail.setProperty("role", "secondary")
        self.progress_detail.setWordWrap(True)
        left_layout.addWidget(self.progress)
        left_layout.addWidget(self.progress_detail)

        # 只读日志框显示编译输出、模型标准输出/错误信息等。
        self.log = QTextEdit()
        self.log.setReadOnly(True)
        self.log.setFont(fixed_width_font())
        self.log.setPlaceholderText("Build and model output will appear here.")
        self.log.setMinimumHeight(135)
        self.log.setMaximumHeight(190)
        left_layout.addWidget(self.log)
        self.main_splitter.addWidget(self.configuration_panel)

        # ===== 右侧 Card：实时诊断与场可视化 =====
        right = QFrame()
        right.setObjectName("Card")
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(16, 14, 16, 14)
        right_layout.setSpacing(9)

        # 可视化 Header：标题说明 + 右侧地图/球面等 view mode 分段控件。
        visual_header = QHBoxLayout()
        visual_title_layout = QVBoxLayout()
        visual_title_layout.setSpacing(1)
        visual_eyebrow = QLabel("Visualization")
        visual_eyebrow.setProperty("role", "eyebrow")
        visual_title = QLabel("Live Climate Field")
        visual_title.setProperty("role", "section")
        visual_title_layout.addWidget(visual_eyebrow)
        visual_title_layout.addWidget(visual_title)
        visual_header.addLayout(visual_title_layout)
        visual_header.addStretch(1)
        # view_mode 的候选值和显示标签来自参数目录，避免在 UI 中重复硬编码。
        view_definition = self.parameter_catalog.by_name("view_mode")
        view_labels = view_definition["choice_labels"]
        self.view_mode = SegmentedControl(
            [(view_labels[value], value) for value in view_definition["choices"]]
        )
        visual_header.addWidget(self.view_mode)
        right_layout.addLayout(visual_header)

        # 可视化工具条：变量选择、垂直层选择、历史 frame stream 回放。
        visual_controls = QHBoxLayout()
        visual_controls.setSpacing(8)
        self.variable_combo = QComboBox()
        # ComboBox 的显示文本是“变量名 (单位)”，item data 保存稳定的内部变量名。
        for variable in self.variable_catalog.variables:
            self.variable_combo.addItem(
                f"{variable['label']} ({variable['unit']})",
                variable["name"],
            )
        # 默认优先选择 surface_temperature；若目录中不存在，则退回第 0 项。
        default_index = self.variable_combo.findData("surface_temperature")
        self.variable_combo.setCurrentIndex(max(default_index, 0))
        visual_controls.addWidget(self.variable_combo, 1)
        # 垂直层下拉框会在收到具体 frame 后，根据 frame.layers 动态重建。
        self.level_label = QLabel("Sigma model level")
        self.level_label.setProperty("role", "secondary")
        visual_controls.addWidget(self.level_label)
        self.level_combo = QComboBox()
        self.level_combo.addItem("0", 0)
        self.level_combo.setMinimumWidth(72)
        visual_controls.addWidget(self.level_combo)
        self.time_label = QLabel("Time")
        self.time_label.setProperty("role", "secondary")
        visual_controls.addWidget(self.time_label)
        self.time_combo = QComboBox()
        self.time_combo.setMinimumWidth(205)
        self.time_combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        visual_controls.addWidget(self.time_combo)
        self.play_button = QPushButton("Play")
        self.play_button.setCheckable(True)
        visual_controls.addWidget(self.play_button)
        self.playback_speed = QComboBox()
        for label, interval in (("1×", 800), ("2×", 400), ("4×", 200)):
            self.playback_speed.addItem(label, interval)
        self.playback_speed.setCurrentIndex(1)
        self.playback_speed.setMaximumWidth(66)
        visual_controls.addWidget(self.playback_speed)
        self.open_frames_button = QPushButton("Replay Stream")
        visual_controls.addWidget(self.open_frames_button)
        right_layout.addLayout(visual_controls)
        # LiveFieldView 是实际绘制 2D/Globe/时间序列等内容的自定义可视化组件。
        self.live_view = LiveFieldView()
        right_layout.addWidget(self.live_view, 1)
        self.main_splitter.addWidget(right)
        # 初始 splitter 宽度比例偏向右侧可视化区域；最后统一应用 Studio 主题样式。
        self._configuration_width = 430
        self.main_splitter.setSizes([self._configuration_width, 1050])
        self.setStyleSheet(STUDIO_STYLESHEET)

    # -------------------------------------------------------------------------
    # 统一建立 Qt signal-slot 连接。
    # Review 时可把本函数看成“用户操作/进程事件 → 对应处理函数”的路由表。
    # -------------------------------------------------------------------------
    def _connect_signals(self) -> None:
        # 参数变化：进入统一的 _parameter_changed，以便处理 build 失效和联动更新。
        self.form.valueChanged.connect(self._parameter_changed)
        # Preset 与 build/run 控件绑定。
        self.defaults_button.clicked.connect(self.form.reset_defaults)
        self.save_button.clicked.connect(self._save_preset)
        self.load_button.clicked.connect(self._load_preset)
        self.build_button.clicked.connect(self._build)
        self.run_button.clicked.connect(self._run)
        # Stop 直接委托给 ProcessController；Controller 决定如何终止正在运行的进程。
        self.stop_button.clicked.connect(self.controller.stop_run)
        self.open_frames_button.clicked.connect(self._open_frames)
        self.help_button.clicked.connect(self._show_help)
        self.configuration_toggle.toggled.connect(self._toggle_configuration)
        # 快捷时长按钮通过 lambda 固定传入模拟天数。
        self.quick_1day.clicked.connect(lambda: self._apply_day_preset(1))
        self.quick_10days.clicked.connect(lambda: self._apply_day_preset(10))
        self.quick_30days.clicked.connect(lambda: self._apply_day_preset(30))
        self.quick_1year.clicked.connect(lambda: self._apply_day_preset(365))
        # view mode 改变无需经过 MainWindow 再包装，直接推给 LiveFieldView。
        self.view_mode.valueChanged.connect(self.live_view.set_view_mode)
        # 变量或垂直层改变时，从现有 frame buffer 立即重绘最新帧。
        self.variable_combo.currentIndexChanged.connect(self._variable_changed)
        self.level_combo.currentIndexChanged.connect(self._show_latest_frame)
        self.time_combo.currentIndexChanged.connect(self._time_changed)
        self.play_button.toggled.connect(self._toggle_playback)
        self.playback_speed.currentIndexChanged.connect(self._playback_speed_changed)
        self.live_view.primary_control.valueChanged.connect(self._visualization_mode_changed)
        # 外部进程控制器发出的事件：日志、状态、build 完成、run 完成。
        self.controller.logReceived.connect(self._append_log)
        self.controller.stateChanged.connect(self._set_state)
        self.controller.buildFinished.connect(self._build_finished)
        self.controller.runFinished.connect(self._run_finished)
        # signal 连接后主动做一次初始化同步，确保色标和时长提示与当前默认值一致。
        self._variable_changed()
        self._update_run_hint()

    def _toggle_configuration(self, visible: bool) -> None:
        """Collapse the configuration card while keeping its header toggle available."""
        sizes = self.main_splitter.sizes()
        if not visible and sizes and sizes[0] > 0:
            self._configuration_width = sizes[0]
        self.configuration_panel.setVisible(visible)
        self.configuration_toggle.setText(
            "Hide Configuration" if visible else "Show Configuration"
        )
        if visible:
            total = max(sum(sizes), self.width())
            width = min(self._configuration_width, max(280, total - 500))
            self.main_splitter.setSizes([width, max(1, total - width)])

    # 用户切换变量后：更新该变量对应的 colormap，并用缓存中的最新帧重绘。
    def _variable_changed(self) -> None:
        # currentData() 取的是 addItem 时存入的内部变量名，而不是界面显示文本。
        name = str(self.variable_combo.currentData())
        definition = self.variable_catalog.by_name(name)
        self.live_view.set_variable_definition(definition)
        self.live_view.set_colormap(str(definition["colormap"]))
        # 海洋层变量（如深海温度/盐度）的垂直轴是深度层，而非大气 sigma 层。
        dimensions = tuple(definition.get("dimensions", ()))
        self.level_label.setText(
            "Ocean layer"
            if definition.get("component") == "ocean" and "level" in dimensions
            else "Sigma model level"
        )
        self._follow_latest_time = True
        self._refresh_time_controls()
        self._show_latest_frame()

    def _visualization_mode_changed(self, mode: str) -> None:
        self.view_mode.setVisible(mode == "global")

    def _refresh_time_controls(self) -> None:
        name = str(self.variable_combo.currentData())
        frames = self.frame_buffer.history(name)
        previous_index = self.time_combo.currentIndex()
        previous_text = self.time_combo.currentText()
        self._updating_time_controls = True
        self.time_combo.blockSignals(True)
        self.time_combo.clear()
        for index, frame in enumerate(frames):
            self.time_combo.addItem(
                f"{frame.year:04d}-{frame.month:02d}-{frame.day:02d} "
                f"{frame.hour:02d}:{frame.minute:02d} · s{frame.step}",
                index,
            )
        if frames:
            preserved = self.time_combo.findText(previous_text)
            target = len(frames) - 1 if self._follow_latest_time else (
                preserved if preserved >= 0 else min(max(previous_index, 0), len(frames) - 1)
            )
            self.time_combo.setCurrentIndex(target)
        self.time_combo.blockSignals(False)
        self._updating_time_controls = False
        has_time = len(frames) > 1
        self.time_combo.setEnabled(bool(frames))
        self.play_button.setEnabled(has_time)
        self.playback_speed.setEnabled(has_time)
        if not has_time and self.play_button.isChecked():
            self.play_button.setChecked(False)

    def _time_changed(self) -> None:
        if not self._updating_time_controls:
            self._follow_latest_time = self.time_combo.currentIndex() == self.time_combo.count() - 1
        self._show_latest_frame()

    def _toggle_playback(self, playing: bool) -> None:
        self.play_button.setText("Pause" if playing else "Play")
        if playing and self.time_combo.count() > 1:
            self._follow_latest_time = False
            self.playback_timer.start(int(self.playback_speed.currentData()))
        else:
            self.playback_timer.stop()

    def _playback_speed_changed(self) -> None:
        if self.playback_timer.isActive():
            self.playback_timer.setInterval(int(self.playback_speed.currentData()))

    def _advance_playback(self) -> None:
        count = self.time_combo.count()
        if count < 2:
            self.play_button.setChecked(False)
            return
        self.time_combo.setCurrentIndex((self.time_combo.currentIndex() + 1) % count)

    # 参数表单任意值变化时触发。这里处理会影响构建产物或其他字段的联动逻辑。
    def _parameter_changed(self, name: str, _value: object) -> None:
        definition = self.parameter_catalog.by_name(name)
        # compile-phase 参数变化意味着现有 executable 可能已与配置不一致，
        # 因此把 current_build 清空并在状态标签上明确提示“构建配置已变化”。
        if definition["phase"] == "compile":
            self.current_build = None
            self.state_label.setText("Build configuration changed")
        # resolution 改变时，restart 和大气垂直层数也必须同步到该分辨率的默认值。
        if name == "resolution":
            self._sync_resolution_defaults(str(_value))
        # resolution 或 steps 改变都会影响“模拟天数”换算提示。
        if name in {"resolution", "steps"}:
            self._update_run_hint()
        if name == "auto_levels":
            self.live_view.set_color_mode("dynamic" if bool(_value) else "fixed")

    # 当 resolution 改变时同步与该网格绑定的 restart 文件和大气层数。
    def _sync_resolution_defaults(self, resolution: str) -> None:
        # 临时构造 BuildConfig，只为查询该分辨率下定义好的默认 restart/levels。
        build = BuildConfig(resolution=resolution)
        updates: dict[str, object] = {
            "restart": build.default_restart,
            "atmosphere_levels": build.default_atmosphere_levels,
        }
        # 通过 form.set_values 一次性更新关联字段；是否继续发 valueChanged 由表单实现决定。
        self.form.set_values(updates)

    # 根据当前表单中的 resolution 查询该分辨率标定的“每模拟日步数”。
    def _steps_per_day(self) -> int:
        resolution = str(self.form.value("resolution"))
        return BuildConfig(resolution=resolution).steps_per_day

    # 把用户选择的模拟天数转换为 integration steps。
    def _apply_day_preset(self, days: int) -> None:
        # 若 resolution 不存在于 BuildConfig 的已知映射中，直接忽略本次快捷设置。
        try:
            steps_per_day = self._steps_per_day()
        except KeyError:
            return
        # steps = 天数 × steps_per_day；真正的数据校验仍由后续 config.validate() 完成。
        self.form.set_values({"steps": days * steps_per_day})

    # 根据当前 resolution 和 steps 刷新 UI 上的运行时长说明。
    def _update_run_hint(self) -> None:
        # 表单值可能在编辑过程中暂时不可转换，因此这里只做防御式读取，失败则不更新提示。
        try:
            resolution = str(self.form.value("resolution"))
            steps = int(self.form.value("steps"))
            steps_per_day = BuildConfig(resolution=resolution).steps_per_day
        except (KeyError, ValueError, TypeError):
            return
        # 防御 steps_per_day=0，避免除零；正常配置下该值应为正数。
        days = steps / steps_per_day if steps_per_day else 0.0
        self.run_hint.setText(
            f"Current run length: {steps} steps \u2248 {days:.1f} days "
            f"({steps_per_day} steps/day at {resolution})"
        )

    # 弹出静态 Quick Guide。内容使用 Qt 支持的 HTML 子集进行排版。
    def _show_help(self) -> None:
        QMessageBox.information(
            self,
            "PlaSiC \u2014 Quick Guide",
            # 以下字符串在 Python 中会自动拼接为一个 HTML 文本参数。
            "<h3>PlaSiC Studio</h3>"
            "<p>PlaSiC is a simplified general circulation model (GCM) for "
            "Earth's atmosphere, ocean, sea ice and land surface. "
            "This app lets you configure, run and visualize experiments "
            "without leaving the desktop.</p>"
            "<h4>Running your first experiment</h4>"
            "<ol>"
            "<li><b>Build Model</b> compiles the C core for the selected "
            "resolution.</li>"
            "<li>Pick a length with <b>Quick Start</b> (in days) or set "
            "<i>Integration steps</i> directly.</li>"
            "<li><b>Start Run</b> launches the model; live fields and the "
            "global-mean time series update as it runs.</li>"
            "<li><b>Stop</b> terminates the run early if needed.</li>"
            "</ol>"
            "<h4>Exploring results</h4>"
            "<p>Switch between the <b>2D Map</b> and the interactive "
            "<b>Globe</b> (drag to rotate, scroll to zoom, double-click to "
            "reset). Choose any of the live variables; 3D atmosphere "
            "variables also offer a vertical-level selector.</p>"
            "<h4>Replay Stream \u2014 replaying a saved run</h4>"
            "<p>Every run writes a binary <code>frames.pstream</code> file "
            "inside its run directory. <b>Replay Stream</b> lets you load "
            "any of these files to revisit the simulation after it has "
            "finished:</p>"
            "<ol>"
            "<li>Click <b>Replay Stream</b> in the visualization toolbar."
            "</li>"
            "<li>In the file dialog, navigate to the <code>runs/</code> "
            "folder (inside the project root) and pick the "
            "<code>frames.pstream</code> from the run you want. Each run "
            "lives in a timestamped sub-folder such as "
            "<code>runs/20250818T120000Z-a1b2c3d4/</code>.</li>"
            "<li>The map, globe and global-mean time series populate from "
            "the stream. Frames load in batches of 128 per refresh, so "
            "very long streams appear progressively rather than all at "
            "once.</li>"
            "<li>While replaying, switch variables or vertical levels "
            "with the dropdowns \u2014 the visualisation updates from the "
            "buffered history without re-reading the file.</li>"
            "<li>The replay does not re-execute the model. To continue a "
            "run from its final state, start a new run and point "
            "<i>Initial restart</i> at the <code>final.restart</code> file "
            "inside the original run directory.</li>"
            "</ol>"
            "<h4>Reproducibility</h4>"
            "<p>Every run is saved under <code>runs/&lt;run-id&gt;</code> with "
            "its config, build info, log and <code>frames.pstream</code>. "
            "Use <b>Save Preset</b> / <b>Load Preset</b> to share exact "
            "configurations, and <b>Replay Stream</b> to revisit any saved "
            "frame stream.</p>",
        )

    # -------------------------------------------------------------------------
    # 从当前 UI 状态组装完整 ExperimentConfig。
    # 这是“界面状态 → 领域配置对象”的集中转换点。
    # -------------------------------------------------------------------------
    def _configuration(self) -> ExperimentConfig:
        # values() 返回 ParameterForm 当前所有参数的统一字典。
        values = self.form.values()
        # BuildConfig/RunConfig 只取各自 dataclass 声明的字段，
        # 避免把 UI 中属于其他阶段的参数误传给构造函数。
        build = BuildConfig(**{key: values[key] for key in asdict(BuildConfig())})
        run = RunConfig(**{key: values[key] for key in asdict(RunConfig())})
        # 可视化配置并非全部来自 ParameterForm；view mode、变量和层级来自右侧控件。
        visualization = VisualizationConfig(
            buffer_frames=int(values["buffer_frames"]),
            auto_levels=self.live_view.auto_levels,
            view_mode=self.view_mode.value(),
            selected_variable=str(self.variable_combo.currentData()),
            selected_level=int(self.level_combo.currentData() or 0),
        )
        # 最后将三个子配置组合成一次实验的不可分割配置对象。
        return ExperimentConfig(build=build, run=run, visualization=visualization)

    # 创建配置并执行领域级校验；失败时集中弹窗并返回 None。
    def _validated_configuration(
        self,
    ) -> ExperimentConfig | None:
        # validate() 返回错误字符串列表，MainWindow 只负责展示，不复制校验规则。
        config = self._configuration()
        errors = config.validate()
        if errors:
            QMessageBox.critical(self, "Invalid Parameters", "\n".join(errors))
            return None
        return config

    # 用户显式点击 Build Model 的入口。
    def _build(self) -> None:
        # Build 前始终重新读取并校验当前 UI，避免使用陈旧配置。
        config = self._validated_configuration()
        if config is None:
            return
        # 显式 Build 不应该在完成后自动启动模型，因此清除 pending_run。
        self.pending_run = False
        # build_invocation 把 BuildConfig + C 根目录解析成具体编译命令/目标产物描述。
        self.current_build = build_invocation(config.build, self.c_root)
        # 真正的异步/外部进程执行交给 ProcessController。
        self.controller.start_build(self.current_build)

    # 用户点击 Start Run 的入口：必要时先 build，再自动续接运行。
    def _run(self) -> None:
        config = self._validated_configuration()
        if config is None:
            return
        # 始终基于当前 build 配置计算期望的 executable 路径和构建命令。
        invocation = build_invocation(config.build, self.c_root)
        self.current_build = invocation
        # 如果目标可执行文件不存在，则把本次 run 标记为“等待 build 完成”。
        # _build_finished() 会检查 pending_run，并在成功后继续 _start_model()。
        if not invocation.executable.is_file():
            self.pending_run = True
            self.controller.start_build(invocation)
            return
        # 已有 executable 时跳过构建，直接准备并启动模型。
        self._start_model(config, invocation)

    # -------------------------------------------------------------------------
    # 真正启动一次模型实验
    # -------------------------------------------------------------------------
    def _start_model(
        self,
        config: ExperimentConfig,
        build: BuildInvocation,
    ) -> None:
        # RunManager.prepare 会创建 run-id 目录、写入 config.json，
        # 并准备好日志/进度/frame 路径等。
        # 任何文件系统或配置准备错误都在此转成用户可见弹窗。
        try:
            prepared = self.run_manager.prepare(config, self.c_root)
        except (OSError, ValueError) as error:
            QMessageBox.critical(self, "Cannot Create Experiment", str(error))
            return
        
        # 保存本次 PreparedRun，并建立两个“增量读取器”：progress 与 frame stream。
        self.current_run = prepared
        self.progress_tail = ProgressTail(prepared.progress)
        self.frame_reader = PStreamReader(prepared.frames)

        # 为这次实验重建帧缓存，容量来自 visualization.buffer_frames。
        self.frame_buffer = FrameBuffer(config.visualization.buffer_frames)
        self._stream_ocean_levels = None
        self._stream_variables.clear()
        self._stream_complete = False
        self._replay_active = False
        self.live_view.reset_data_context()
        self._follow_latest_time = True
        self._refresh_time_controls()

        # LiveFieldView 的时间历史也使用同一上限，保持 UI 与 FrameBuffer 的留存策略一致。
        self.live_view.set_history_limit(config.visualization.buffer_frames)

        # 新实验启动前重置进度，并在详情标签中显示可追踪的 run_id。
        self.progress.setValue(0)
        self.progress_detail.setText(f"Experiment: {prepared.run_id}")

        # run_invocation 生成真正启动数值核心的命令；
        # C runtime 参数通过命令行直接传递，不再依赖 model_config.json。
        invocation = run_invocation(
            config,
            self.c_root,
            build.executable,
            prepared.directory,
            prepared.initial_restart,
        )

        # controller.start_run 负责子进程生命周期，并把输出写入 prepared.log。
        self.controller.start_run(invocation, prepared.log)

    # Build 子进程结束后的回调。success/message 由 ProcessController 通过 signal 提供。
    def _build_finished(self, success: bool, message: str) -> None:
        # 构建失败时必须取消 pending_run，避免后续错误地自动启动。
        if not success:
            self.pending_run = False
            QMessageBox.critical(self, "Build Failed", message)
            return
        # 成功时把一条明确日志追加到 UI。
        self._append_log("Build completed.\n")
        # 如果这次 build 是 Start Run 自动触发的，则在成功后续接模型启动。
        if self.pending_run:
            self.pending_run = False
            # Build 过程中用户仍可能改动表单，所以这里重新校验当前配置。
            # 只有 config 有效且 current_build 仍存在时才启动。
            config = self._validated_configuration()
            if config is not None and self.current_build is not None:
                self._start_model(config, self.current_build)

    # 模型进程结束后的回调。
    def _run_finished(self, success: bool, exit_code: int) -> None:
        # 进程刚退出时文件尾可能还有尚未读入的数据，因此先做最后一次输出轮询。
        self._poll_outputs()
        # 正常结束只记日志；非零退出码弹出警告。
        if success:
            self._append_log("Model run completed.\n")
        elif exit_code != 0:
            QMessageBox.warning(self, "Model Stopped", f"Exit code: {exit_code}")

    # 根据 controller 状态统一刷新状态文本、样式属性以及 Build/Run/Stop 可用性。
    def _set_state(self, state: str) -> None:
        # 内部状态值与用户可读文案分离，便于状态机和 UI 文案独立维护。
        self._state = state
        labels = {
            "idle": "Ready",
            "building": "Building",
            "running": "Model Running",
            "stopping": "Stopping",
        }
        self.state_label.setText(labels.get(state, state))
        # Qt 动态 property 改变后，主动 unpolish/polish 以立即重新应用 stylesheet 选择器。
        self.state_label.setProperty("state", state)
        self.state_label.style().unpolish(self.state_label)
        self.state_label.style().polish(self.state_label)
        # building/running/stopping 都视为“活动中”，防止重复发起 Build/Run。
        # Stop 只在真正 running 时启用。
        active = state in {"building", "running", "stopping"}
        self.build_button.setEnabled(not active)
        self.run_button.setEnabled(not active)
        self.stop_button.setEnabled(state == "running")

    # 把新日志追加到 QTextEdit 末尾，并自动滚动到最新内容。
    def _append_log(self, text: str) -> None:
        self.log.moveCursor(QTextCursor.MoveOperation.End)
        self.log.insertPlainText(text)
        self.log.ensureCursorVisible()

    # -------------------------------------------------------------------------
    # 定时器每 150 ms 调用：增量读取 progress 与 frame stream。
    # 这是实时 UI 数据流的核心函数，必须控制单次工作量以避免卡住 GUI 线程。
    # -------------------------------------------------------------------------
    def _poll_outputs(self) -> None:
        # ===== 1) 读取新的进度记录 =====
        if self.progress_tail is not None:
            # read_new() 只返回自上次读取位置之后新增的进度记录。
            try:
                for event in self.progress_tail.read_new():
                    self._handle_progress(event)
            # 读取异常不会让整个 GUI 崩溃；记录错误后禁用 progress reader。
            except (OSError, ValueError) as error:
                self._append_log(f"Progress read error: {error}\n")
                self.progress_tail = None

        # ===== 2) 读取新的场数据帧 =====
        if self.frame_reader is not None:
            try:
                # 原注释强调两层内存约束：
                # - 每次 GUI poll 最多读取 128 帧；
                # - FrameBuffer/LiveFieldView 只保留配置指定长度的历史。
                # max_frames=128 是关键的背压/响应性保护，长历史文件不会一次性全部读入内存。
                frames = self.frame_reader.read_new(max_frames=128)

            # frame stream 损坏或读取失败时记录错误、停用 reader，并结束本次 poll。
            except (OSError, ValueError) as error:
                self._append_log(f"Frame stream read error: {error}\n")
                self.frame_reader = None
                return

            # 有新帧时先全部进入统一缓存，再只筛选当前用户正在查看的变量用于即时绘制。
            capabilities_changed = False
            # 文件头只在第一次成功读取时解析；一拿到就更新流能力。
            header = self.frame_reader.header
            if header is not None and header.ocean_levels != self._stream_ocean_levels:
                self._stream_ocean_levels = header.ocean_levels
                capabilities_changed = True

            if frames:
                self.frame_buffer.extend(frames)
                # 记录流中出现过的变量，旧帧流即使文件头未声明深海层数，
                # 也能通过实际出现的三维海洋帧识别能力。
                names = {frame.name for frame in frames}
                if not names <= self._stream_variables:
                    self._stream_variables.update(names)
                    capabilities_changed = True

                # selected 是变量内部 name；不同变量的帧可能交错出现在同一个 stream 中。
                selected = str(self.variable_combo.currentData())
                selected_frames = [frame for frame in frames if frame.name == selected]
                # 只有本批次包含当前变量时才更新右侧视图。
                if selected_frames:
                    self._refresh_time_controls()
                    self._show_latest_frame()

            # 回放/运行结束且文件已读尽时，才能断定某变量确实不存在。
            # 放在渲染之后，避免切换变量时的提示被本次渲染覆盖。
            if (
                not self._stream_complete
                and self._state != "running"
                and self.frame_reader.exhausted
            ):
                self._stream_complete = True
                capabilities_changed = True

            if capabilities_changed:
                self._apply_stream_capabilities()

    # 解析一条 progress.txt 进度记录，并映射到进度条/详情文本。
    def _handle_progress(self, event: dict[str, Any]) -> None:
        # total 至少取 1，防止异常记录中的 total_steps=0 导致除零。
        completed = int(event["completed_steps"])
        total = max(int(event["total_steps"]), 1)
        # 进度映射到 0..1000，与 _build_ui 中 QProgressBar 的范围一致。
        self.progress.setValue(round(1000 * completed / total))
        # date 是模型模拟时间，而不是本机 wall-clock 时间。
        date = event["date"]
        self.progress_detail.setText(
            f"{completed}/{total} steps · "
            f"{date['year']:04d}-{date['month']:02d}-"
            f"{date['day']:02d} "
            f"{date['hour']:02d}:{date['minute']:02d} · "
            f"{event['steps_per_second']:.2f} step/s"
        )

    # 当变量或层级变化时，从 FrameBuffer 中找到该变量最新帧并立即重绘。
    def _show_latest_frame(self) -> None:
        name = str(self.variable_combo.currentData())
        frames = self.frame_buffer.history(name)
        self.live_view.set_frames(frames)
        index = self.time_combo.currentIndex()
        frame = frames[index] if frames and 0 <= index < len(frames) else (frames[-1] if frames else None)
        # 尚未收到该变量帧时清空旧场并说明原因，避免把上一个变量的图像
        # 继续当作当前变量显示。
        if frame is None:
            self.live_view.show_unavailable(
                self._missing_variable_message(name),
            )
            return
        # _level_for_frame 同时负责确保 level_combo 与 frame.layers 一致。
        level = self._level_for_frame(frame)
        self.live_view.show_frame(
            frame,
            level,
            self.live_view.auto_levels,
        )

    def _requires_deep_ocean(self, definition: dict[str, Any]) -> bool:
        """目录中以三维海洋层为垂直轴的变量（深海温度、盐度）。"""
        return (
            definition.get("component") == "ocean"
            and "level" in tuple(definition.get("dimensions", ()))
        )

    def _stream_has_deep_ocean_frames(self) -> bool:
        for name in self._stream_variables:
            try:
                definition = self.variable_catalog.by_name(name)
            except KeyError:
                continue
            if self._requires_deep_ocean(definition):
                return True
        return False

    def _deep_ocean_available(self) -> bool:
        """三维海洋变量在当前帧流中是否可用。

        新写入的帧流会在文件头声明深海层数；旧帧流可能仍写 1，但流中
        实际出现了三维海洋帧，因此已见变量也作为可用性依据。只有在数据
        来源读尽且仍未见到三维海洋帧时才判定为不可用。
        """
        if self._stream_has_deep_ocean_frames():
            return True
        if self._stream_ocean_levels is None:
            return True
        if self._stream_ocean_levels > 1:
            return True
        return not self._stream_complete

    def _deep_ocean_unavailable_reason(self) -> str:
        return (
            "the current frame stream contains no 3-D ocean frames; "
            "if this run should have them, rebuild with OCEAN_MODEL=3d"
        )

    def _missing_variable_message(self, name: str) -> str:
        try:
            definition = self.variable_catalog.by_name(name)
        except KeyError:
            return f"No {name} frames in the current source."
        label = str(definition.get("label", name))
        if self._requires_deep_ocean(definition) and not self._deep_ocean_available():
            return f"{label} is unavailable: {self._deep_ocean_unavailable_reason()}."
        if self._replay_active:
            return f"The loaded frame stream contains no {label} frames."
        if self._state == "running":
            return f"Waiting for the model to output {label} frames…"
        if self.frame_reader is None:
            return "Waiting for model output…"
        return f"No {label} frames in this run."

    def _apply_stream_capabilities(self) -> None:
        """根据帧流能力启用/禁用目录中无法产生的变量。"""
        available = self._deep_ocean_available()
        tooltip = "" if available else (
            "Unavailable: " + self._deep_ocean_unavailable_reason() + "."
        )
        model = self.variable_combo.model()
        for definition in self.variable_catalog.variables:
            if not self._requires_deep_ocean(definition):
                continue
            index = self.variable_combo.findData(str(definition["name"]))
            if index < 0:
                continue
            item = model.item(index)
            item.setEnabled(available)
            item.setToolTip(tooltip)
        current = str(self.variable_combo.currentData())
        try:
            current_definition = self.variable_catalog.by_name(current)
        except KeyError:
            return
        if available or not self._requires_deep_ocean(current_definition):
            return
        fallback = self.variable_combo.findData("sea_surface_temperature")
        if fallback < 0 or not model.item(fallback).isEnabled():
            fallback = next(
                (
                    index
                    for index in range(self.variable_combo.count())
                    if model.item(index).isEnabled()
                ),
                self.variable_combo.currentIndex(),
            )
        self.variable_combo.setCurrentIndex(fallback)
        note = (
            f"{current_definition.get('label', current)} is unavailable: "
            f"{self._deep_ocean_unavailable_reason()}."
        )
        self.live_view.status_label.setText(note)

    # 保证垂直层选择器与当前 frame 的层数一致，并返回最终合法的当前层索引。
    def _level_for_frame(self, frame) -> int:
        # currentData 为空时退回 0，兼容 combo 初始化或清空重建的瞬间状态。
        current_level = int(self.level_combo.currentData() or 0)
        # 仅当层数变化时重建下拉框，避免每帧都重复清空/填充 UI。
        if self.level_combo.count() != frame.layers:
            # 重建过程中临时阻断 signal，防止 clear/add/setCurrentIndex 递归触发重绘。
            self.level_combo.blockSignals(True)
            self.level_combo.clear()
            # 层索引从 0 到 layers-1，显示文本与 item data 都使用该整数索引。
            for level in range(frame.layers):
                self.level_combo.addItem(str(level), level)
            # 尽可能保留用户之前选择的层；若新 frame 层数更少则夹到最大合法层。
            self.level_combo.setCurrentIndex(min(current_level, frame.layers - 1))
            # 重建结束恢复 signal。
            self.level_combo.blockSignals(False)
        has_vertical = frame.layers > 1
        self.level_label.setVisible(has_vertical)
        self.level_combo.setVisible(has_vertical)
        # 返回 combo 当前 data；若仍为空则安全回退到第 0 层。
        return int(self.level_combo.currentData() or 0)

    # 将当前完整 ExperimentConfig 保存为格式化 JSON 预设。
    def _save_preset(self) -> None:
        # Qt 文件对话框返回 (path, selected_filter)；这里不需要 selected_filter。
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Experiment Preset",
            "plasic-preset.json",
            "JSON (*.json)",
        )
        # 用户取消保存时 path 为空，直接返回。
        if not path:
            return
        # 保存的是完整 configuration，而不是只保存 ParameterForm，因而包含可视化配置。
        config = self._configuration()
        # ensure_ascii=False 保留非 ASCII 字符；indent/sort_keys 让文件稳定且适合 diff/review。
        # 末尾显式补 '\n'，符合常见文本文件约定。
        Path(path).write_text(
            json.dumps(
                config.to_dict(),
                indent=2,
                sort_keys=True,
                ensure_ascii=False,
            )
            + "\n",
            encoding="utf-8",
        )

    # 从 JSON 预设恢复实验配置，并把不同子配置回填到对应 UI 控件。
    def _load_preset(self) -> None:
        # 只允许选择 JSON 文件；用户取消则直接返回。
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Experiment Preset", "", "JSON (*.json)"
        )
        if not path:
            return
        # 文件读取、JSON 解析和领域对象反序列化统一放在 try 中。
        try:
            # ExperimentConfig.from_dict 是 JSON 字典 → 强类型配置对象的入口。
            config = ExperimentConfig.from_dict(
                json.loads(Path(path).read_text(encoding="utf-8"))
            )
            # BuildConfig 与 RunConfig 展平成 ParameterForm 可一次性消费的 values 字典；
            # visualization 中属于表单的两个字段也一并回填。
            values = {
                **asdict(config.build),
                **asdict(config.run),
                "buffer_frames": config.visualization.buffer_frames,
                "auto_levels": config.visualization.auto_levels,
            }
            # set_values 负责批量刷新表单控件。
            self.form.set_values(values)
            self.live_view.set_color_mode(
                "dynamic" if config.visualization.auto_levels else "fixed"
            )
            # selected_variable 用稳定的 item data 查找，而不是依赖可能变化的显示 label。
            index = self.variable_combo.findData(config.visualization.selected_variable)
            # 预设中的变量若当前目录仍存在才切换，避免 findData=-1 导致非法索引。
            if index >= 0:
                self.variable_combo.setCurrentIndex(index)
            # view mode 由 SegmentedControl 自己恢复。
            self.view_mode.set_value(config.visualization.view_mode)
        # 无论是磁盘错误、字段类型错误、领域校验错误还是 JSON 语法错误，都转成统一弹窗。
        except (OSError, ValueError, TypeError, json.JSONDecodeError) as error:
            QMessageBox.critical(self, "Invalid Preset", str(error))

    # 打开历史 frames.pstream，用与实时运行相同的 reader/buffer/view 管线进行回放。
    def _open_frames(self) -> None:
        # 文件对话框优先过滤 PlaSiC stream，同时允许用户选择其他扩展名文件。
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Open PlaSiC Frame Stream",
            "",
            "PlaSiC stream (*.pstream);;All files (*)",
        )
        # 用户取消选择时不改变当前 reader/buffer。
        if not path:
            return
        
        # 新建 reader 并重置 FrameBuffer，避免历史回放与之前运行的数据混在一起。
        self.frame_reader = PStreamReader(Path(path))
        self.frame_buffer = FrameBuffer(int(self.form.value("buffer_frames")))
        self._stream_ocean_levels = None
        self._stream_variables.clear()
        self._stream_complete = False
        self._replay_active = True
        self.live_view.reset_data_context()
        self._follow_latest_time = True
        self._refresh_time_controls()

        # LiveFieldView 的历史容量同步为当前表单中的 buffer_frames。
        self.live_view.set_history_limit(int(self.form.value("buffer_frames")))
        
        # 立即主动 poll 一次，让用户选择文件后无需等待下一个 150 ms 定时器 tick。
        self._poll_outputs()
