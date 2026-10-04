#!/usr/bin/env python3
"""Error-growth figures built from the verification archives.

中文说明：基于检验归档数据绘制误差增长图。

English: the script deliberately uses large panels and reads every score from
``data/verification`` so figures cannot silently diverge from the tables in
the report.

中文说明：脚本刻意使用大尺寸面板，并且所有评分都从 ``data/verification``
读取，保证图件不会悄悄偏离报告中的表格数值。
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

OUT = common.FIGURE_ROOT  # 图件输出目录 / figure output directory
CASES = [c.key for c in common.CASES]
CASE_LABELS = {"winter2020": "Winter", "summer2020": "Summer", "typhoon2020": "Typhoon"}
# 中文说明：图例中的个例英文简称 / English short labels for the cases.
METHOD_LABELS = {"direct": "Direct analysis", "nudged": "12-h nudging"}
# 中文说明：两种初值方法的图例标签 / legend labels of the two initialization methods.
COLORS = {"direct": "#1B4F72", "nudged": "#C45A3C"}
# 中文说明：两种初值方法的曲线颜色 / curve colours of the two initialization methods.

# 中文说明：全局出版风格参数（字体、线宽、刻度等）。
# English: global publication-style settings (fonts, line widths, ticks, ...).
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 10,
    "axes.labelsize": 10,
    "axes.titlesize": 11,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "legend.frameon": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})


def load(case: str, method: str) -> dict[str, np.ndarray]:
    """Load one verified case/method archive.

    中文说明：载入某个“个例 + 初值方法”的检验归档文件。
    """
    return dict(np.load(common.DATA_ROOT / "verification" / f"{case}_{method}.npz"))


def save(fig: plt.Figure, stem: str) -> None:
    """Save editable vector files and a 600-dpi raster preview.

    中文说明：把图件以 300 dpi 的 PNG 保存并关闭。
    """
    OUT.mkdir(parents=True, exist_ok=True)
    base = OUT / stem
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    plt.close(fig)


def pick(data: dict, key: str, level: int | None, scale: float = 1.0) -> np.ndarray:
    """Select an optional level and apply a unit-conversion scale.

    中文说明：可选地选取某个层次，并施加单位换算系数。
    """
    values = data[key]
    if level is not None:
        values = values[level]
    return np.asarray(values, dtype=float) * scale


def overview() -> None:
    """Plot area-weighted RMSE separately for each archived case.

    中文说明：为每个已归档个例单独绘制面积加权 RMSE 增长曲线。

    The report uses these panels as case diagnostics.  Earlier versions
    averaged all three cases and shaded their standard deviation, which made
    the band look like an uncertainty interval and hid case-specific
    behaviour.  Keeping the two initialization methods in each panel while
    plotting one case at a time makes the comparison explicit.
    中文说明：报告用这些面板做单一个例诊断。早期版本把三个个例平均并绘制
    标准差阴影，容易让人误以为是置信区间，而且掩盖了个例之间的差异。现在
    每个面板只画一个例的两种初值方法，比较更加清晰。
    """
    specs = [
        ("rmse_t", 0, 1.0, "T850 RMSE (K)"),
        ("rmse_t", 1, 1.0, "T500 RMSE (K)"),
        ("rmse_z", 0, 1.0, "Z500 RMSE (m)"),
        ("rmse_u", 1, 1.0, "U250 RMSE (m s$^{-1}$)"),
        ("rmse_q", 1, 1000.0, "Q700 RMSE (g kg$^{-1}$)"),
        ("rmse_t2m", None, 1.0, "2-m temperature RMSE (K)"),
    ]
    # 中文说明：specs 中每项为 (归档键, 层次索引, 单位换算, y 轴标签)。
    # English: each spec is (archive key, level index, unit scale, y label).
    for case in CASES:
        fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.6), sharex=True)
        for index, (ax, (key, level, scale, ylabel)) in enumerate(zip(axes.flat, specs)):
            for method in ("direct", "nudged"):
                data = load(case, method)
                lead = data["all_lead_hours"]
                values = pick(data, key, level, scale)
                ax.plot(lead, values, lw=2.2, color=COLORS[method],
                        label=METHOD_LABELS[method])
            ax.set_ylabel(ylabel)
            ax.set_xlim(0, 240)
            ax.set_xticks([0, 48, 96, 144, 192, 240])
            ax.grid(axis="y", color="#D9DEE2", lw=0.7)
            # 面板字母标签 a–f / panel letters a-f
            ax.text(0.02, 0.94, "abcdef"[index], transform=ax.transAxes,
                    fontweight="bold", fontsize=12, va="top")
        for ax in axes[1]:
            ax.set_xlabel("Lead time (h)")
        axes[0, 0].legend(loc="upper left", ncol=2)
        label = CASE_LABELS.get(case, case)
        fig.suptitle(f"Area-weighted RMSE growth: {label} 2020", fontsize=14, y=1.01)
        fig.tight_layout()
        save(fig, f"overview_error_growth_{case}")


def main() -> None:
    overview()


if __name__ == "__main__":
    main()
