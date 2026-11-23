#!/usr/bin/env python3
"""PlaSiC physics time step: call order of the dynamical core and the
physical parameterizations (circular schematic, English labels).

Figure contract
---------------
Core conclusion : In every physics step plasic_physics_step() calls
                  1. moisture fixer -> 2. surface fluxes -> 3. vertical
                  diffusion -> 4. radiation -> 5. moist processes (six
                  sub-processes) -> 6. land surface -> 7. sea ice & slab ocean
                  (optional 3-D ocean), closing the loop
                  "dynamical core -> physics -> surface -> dynamical core".
Evidence chain  : single hero schematic; charcoal arcs = same-step call order;
                  amber dashes = this-step flux forcing (including
                  precipitation from moist processes to the surface); teal
                  dashes = clouds -> radiation (next step); blue dashes =
                  surface boundary for the next step; red arc = tendency
                  feedback to the dynamical core.
Archetype       : schematic-led composite (single schematic panel).
Backend         : Python / matplotlib (saved preference).
Export          : 600 dpi PNG written next to the tutorial images (default);
                  editable-text SVG / PDF on request (--formats).
Typography      : Arial only, sized for a ~700 px wide documentation column.
Palette         : neutral slate cards with a navy anchor; saturated colour is
                  reserved for the five arrow meanings.

Source of truth : src/physics/physics_processes.c (physics_step order),
                  docs-site/docs/3-parameterizations/*,
                  docs-site/docs/4-model-and-coupling/4.1.

Usage
-----
    python plot_plasic_parameterizations.py
    python plot_plasic_parameterizations.py --formats png svg pdf --outdir out
"""

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch

STEM = "plasic_atmospheric_parameterizations"
IMAGES_DIR = Path(__file__).resolve().parent.parent

# ----------------------------------------------------------------------------
# typography: Arial for every glyph
# ----------------------------------------------------------------------------
for _name in ("Arial.ttf", "Arial Bold.ttf", "Arial Italic.ttf",
              "Arial Bold Italic.ttf"):
    _path = Path("/System/Library/Fonts/Supplemental") / _name
    if _path.exists():
        font_manager.fontManager.addfont(str(_path))

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)

FS_CORE = 24.0      # "Dynamical Core"
FS_TITLE = 21.5     # station titles
FS_BADGE = 16.5     # step numbers
FS_TEXT = 18.0      # sub-processes, variable labels, caption, legend
CAP = 0.716         # Arial cap height / em (used to centre text optically)

# ----------------------------------------------------------------------------
# palette: neutral structure, colour only where it carries meaning
# ----------------------------------------------------------------------------
INK = "#1F2A37"          # primary text
MUTED = "#4F5B6A"        # secondary text
NAVY = "#1F3A5F"         # dynamical core and step badges

ATM_FILL = "#EDF1F6"     # atmospheric-physics cards
ATM_EDGE = "#C5D0DD"
PILL_EDGE = "#D3DBE5"

SURF_FILL = "#F3EEE5"    # land / sea-ice / ocean cards (warm stone)
SURF_EDGE = "#D7C9B1"
SURF_GROUP = "#FBF9F5"   # surface container
SURF_DASH = "#AE9C7E"
SURF_TEXT = "#7A6648"

ORDER = "#353E4A"                       # call order (solid)
FORCING = ("#DE8B2C", "#B4661A")        # (line, text) flux forcing now
CLOUD = ("#2A9579", "#1E7A62")          # clouds -> radiation, next step
BOUNDARY = ("#2F6DB3", "#275D9A")       # surface boundary, next step
TENDENCY = ("#C93A30", "#AE2F27")       # tendency feedback (solid)

LW_RING = 2.6
LW_LINK = 2.3
DASH = (0, (3.2, 2.0))                  # scaled by the line width
HEAD = "-|>,head_length=0.45,head_width=0.22"
HEAD_SCALE = 19.0

# ----------------------------------------------------------------------------
# canvas (1 data unit = 0.9 inch, equal aspect)
# ----------------------------------------------------------------------------
W, H = 16.0, 14.0
UNIT_IN = 0.9
PT = 72.0 * UNIT_IN                     # points per data unit

fig = plt.figure(figsize=(W * UNIT_IN, H * UNIT_IN))
ax = fig.add_axes((0.0, 0.0, 1.0, 1.0))
ax.set_xlim(0.0, W)
ax.set_ylim(0.0, H)
ax.set_aspect("equal")
ax.axis("off")
RENDERER = fig.canvas.get_renderer()


@dataclass(frozen=True)
class Box:
    cx: float
    cy: float
    w: float
    h: float

    @property
    def left(self):
        return self.cx - self.w / 2.0

    @property
    def right(self):
        return self.cx + self.w / 2.0

    @property
    def top(self):
        return self.cy + self.h / 2.0

    @property
    def bottom(self):
        return self.cy - self.h / 2.0

    def pt(self, fx, fy, pad=0.07):
        """Point relative to the box: fx, fy in [-1, 1]; pad pushes outwards."""
        return (self.cx + fx * (self.w / 2.0 + pad),
                self.cy + fy * (self.h / 2.0 + pad))


def text_width(text, size, weight="normal"):
    """Rendered width of a single line of text in data units."""
    t = ax.text(0.0, 0.0, text, fontsize=size, fontweight=weight)
    width = t.get_window_extent(renderer=RENDERER).width
    t.remove()
    return width / fig.dpi / UNIT_IN


def baseline(cy, size):
    """Baseline that centres capital letters optically on cy."""
    return cy - CAP * size / PT / 2.0


def put(x, y, text, size, color=INK, weight="normal", ha="center", z=6):
    ax.text(x, baseline(y, size), text, ha=ha, va="baseline", fontsize=size,
            fontweight=weight, color=color, zorder=z)


def rounded(box, fc, ec="none", lw=1.1, rounding=0.17, ls="-", z=3):
    ax.add_patch(FancyBboxPatch(
        (box.left, box.bottom), box.w, box.h,
        boxstyle=f"round,pad=0,rounding_size={rounding}",
        facecolor=fc, edgecolor=ec, linewidth=lw, linestyle=ls, zorder=z))


# ----------------------------------------------------------------------------
# stations
# ----------------------------------------------------------------------------
PAD_X = 0.32
BADGE_R = 0.27
BADGE_GAP = 0.19
H_ONE = 0.96
H_TWO = 1.50


def badge(x, y, number, z=5):
    ax.add_patch(Circle((x, y), BADGE_R, facecolor=NAVY, edgecolor="none",
                        zorder=z))
    put(x, y, str(number), FS_BADGE, color="white", weight="bold", z=z + 1)


def station(cx, cy, title, number=None, sub=None, fill=ATM_FILL,
            edge=ATM_EDGE, color=INK, size=FS_TITLE, pad_x=PAD_X, h=None):
    title_w = text_width(title, size, "bold")
    sub_w = text_width(sub, FS_TEXT) if sub else 0.0
    lead = 2.0 * BADGE_R + BADGE_GAP if number is not None else 0.0
    w = 2.0 * pad_x + lead + max(title_w, sub_w)
    box = Box(cx, cy, w, h or (H_TWO if sub else H_ONE))
    rounded(box, fill, edge)

    if number is not None:
        badge(box.left + pad_x + BADGE_R, cy, number)
        ha, x_text = "left", box.left + pad_x + lead
    else:
        ha, x_text = "center", cx
    if sub is None:
        put(x_text, cy, title, size, color=color, weight="bold", ha=ha)
    else:
        put(x_text, cy + 0.25, title, size, color=color, weight="bold", ha=ha)
        put(x_text, cy - 0.29, sub, FS_TEXT, color=MUTED, ha=ha)
    return box


def process_group(left, top, number, title, items):
    """Card holding the ordered sub-processes of the moist physics."""
    pill_pad, pill_h, pitch = 0.30, 0.52, 0.61
    pad, head_h = 0.22, 0.66
    pill_w = max(text_width(s, FS_TEXT) for s in items) + 2.0 * pill_pad
    head_w = 2.0 * BADGE_R + BADGE_GAP + text_width(title, FS_TITLE, "bold")
    w = max(pill_w, head_w) + 2.0 * pad
    h = pad + head_h + 0.10 + (len(items) - 1) * pitch + pill_h + pad
    box = Box(left + w / 2.0, top - h / 2.0, w, h)
    rounded(box, ATM_FILL, ATM_EDGE, rounding=0.22)

    head_y = top - pad - head_h / 2.0
    head_x = box.cx - head_w / 2.0
    badge(head_x + BADGE_R, head_y, number)
    put(head_x + 2.0 * BADGE_R + BADGE_GAP, head_y, title, FS_TITLE,
        weight="bold", ha="left")

    first = top - pad - head_h - 0.10 - pill_h / 2.0
    for index, text in enumerate(items):
        y = first - index * pitch
        rounded(Box(box.cx, y, pill_w, pill_h), "white", PILL_EDGE, lw=0.9,
                rounding=0.11, z=4)
        put(box.cx, y, text, FS_TEXT, z=5)
    return box


# ----------------------------------------------------------------------------
# arrows
# ----------------------------------------------------------------------------
def arc_control(p1, p2, rad):
    """Control point of matplotlib's arc3 connection (equal-aspect axes)."""
    (x1, y1), (x2, y2) = p1, p2
    return ((x1 + x2) / 2.0 + rad * (y2 - y1),
            (y1 + y2) / 2.0 - rad * (x2 - x1))


def arrow(p1, p2, color, lw=LW_LINK, rad=0.0, dashed=False, z=2):
    """Arc arrow; dashed shafts get a solid head so the tip stays crisp."""
    style = f"arc3,rad={rad}"
    if not dashed:
        ax.add_patch(FancyArrowPatch(
            p1, p2, arrowstyle=HEAD, mutation_scale=HEAD_SCALE, linewidth=lw,
            color=color, connectionstyle=style, shrinkA=0.0, shrinkB=0.0,
            capstyle="round", zorder=z))
        return
    ax.add_patch(FancyArrowPatch(
        p1, p2, arrowstyle="-", linewidth=lw, color=color, linestyle=DASH,
        connectionstyle=style, shrinkA=0.0, shrinkB=0.30 * HEAD_SCALE,
        zorder=z))
    cx, cy = arc_control(p1, p2, rad)
    dx, dy = p2[0] - cx, p2[1] - cy
    norm = math.hypot(dx, dy)
    tail = (p2[0] - 0.02 * dx / norm, p2[1] - 0.02 * dy / norm)
    ax.add_patch(FancyArrowPatch(
        tail, p2, arrowstyle=HEAD, mutation_scale=HEAD_SCALE, linewidth=lw,
        color=color, shrinkA=0.0, shrinkB=0.0, zorder=z))


def rect_edge(box, target, pad=0.07):
    """Where the line from the box centre towards target leaves the box."""
    dx, dy = target[0] - box.cx, target[1] - box.cy
    sx = (box.w / 2.0 + pad) / abs(dx) if dx else float("inf")
    sy = (box.h / 2.0 + pad) / abs(dy) if dy else float("inf")
    s = min(sx, sy)
    return box.cx + dx * s, box.cy + dy * s


def ring(a, b, color=ORDER, lw=LW_RING, rad=-0.20, start=None, end=None):
    p1 = start or rect_edge(a, (b.cx, b.cy))
    p2 = end or rect_edge(b, (a.cx, a.cy))
    arrow(p1, p2, color, lw=lw, rad=rad)


def note(x, y, text, color, ha="center", z=7):
    ax.text(x, y, text, ha=ha, va="center", fontsize=FS_TEXT, color=color,
            linespacing=1.3, zorder=z,
            bbox=dict(boxstyle="round,pad=0.16", facecolor="white",
                      edgecolor="none", alpha=0.95))


def legend(y_top, pitch=0.50, sample=0.62, gap=0.18, col_gap=0.75):
    """Two-row key: solid meanings on the left, dashed meanings to the right."""
    columns = [
        [(ORDER, False, "Call order"),
         (TENDENCY[0], False, "Tendency feedback")],
        [(FORCING[0], True, "Flux forcing (now)"),
         (BOUNDARY[0], True, "Surface boundary (next step)")],
        [(CLOUD[0], True, "Clouds → radiation (next step)")],
    ]
    widths = [sample + gap + max(text_width(t, FS_TEXT) for *_, t in col)
              for col in columns]
    x = (W - sum(widths) - col_gap * (len(columns) - 1)) / 2.0
    for col, width in zip(columns, widths):
        for row, (color, dashed, text) in enumerate(col):
            y = y_top - row * pitch
            arrow((x, y), (x + sample, y), color, lw=2.3, dashed=dashed, z=5)
            put(x + sample + gap, y, text, FS_TEXT, color=MUTED, ha="left")
        x += width + col_gap


# ----------------------------------------------------------------------------
# layout
# ----------------------------------------------------------------------------
def draw():
    dyn = station(8.40, 12.95, "Dynamical Core", fill=NAVY, edge=NAVY,
                  color="white", size=FS_CORE, pad_x=0.50, h=1.08)
    s1 = station(12.25, 11.25, "Moisture fixer", 1)
    s2 = station(13.45, 7.95, "Surface fluxes", 2)
    s3 = station(12.45, 4.70, "Vertical diffusion (PBL)", 3)
    s4 = station(9.00, 2.85, "Radiation (SW + LW)", 4)
    s6 = station(3.30, 7.95, "Land surface", 6, fill=SURF_FILL,
                 edge=SURF_EDGE)
    s7 = station(4.45, 10.95, "Sea ice & slab ocean", 7,
                 sub="(optional 3D ocean)", fill=SURF_FILL, edge=SURF_EDGE)
    s5 = process_group(1.20, 6.22, 5, "Moist processes", [
        "Cumulus convection",
        "Shallow-convection diffusion",
        "Dry convective adjustment",
        "Large-scale condensation",
        "Precipitation fallout",
        "Diagnostic clouds",
    ])

    # surface group: land (6) and sea (7) share one container
    surf = Box(4.30, 9.46, 6.70, 5.08)
    rounded(surf, SURF_GROUP, SURF_DASH, lw=1.4, rounding=0.28,
            ls=(0, (5, 3.5)), z=1)
    ax.text(6.10, 9.45, "Land / ocean split\nby land sea mask",
            ha="center", va="center", fontsize=FS_TEXT, color=SURF_TEXT,
            linespacing=1.3, zorder=5)

    # call order around the ring; the closing arc feeds the dynamical core
    ring(dyn, s1)
    ring(s1, s2)
    ring(s2, s3)
    ring(s3, s4)
    ring(s4, s5, start=s4.pt(-1, 0), end=s5.pt(1, -0.62), rad=-0.15)
    ring(s5, s6, start=s5.pt(0.55, 1), end=s6.pt(0.42, -1), rad=-0.10)
    ring(s6, s7)
    ring(s7, dyn, color=TENDENCY[0], lw=2.9)

    # data exchanged between stations
    edge_x = surf.right + 0.06
    floor_y = surf.bottom - 0.06
    arrow(s2.pt(-1, -0.30), (edge_x, 8.75), FORCING[0], rad=-0.04,
          dashed=True)
    arrow((edge_x, 9.60), s2.pt(-1, 0.45), BOUNDARY[0], rad=0.04,
          dashed=True)
    arrow(s4.pt(0.30, 1), (7.35, floor_y), FORCING[0], rad=0.04, dashed=True)
    arrow((7.00, floor_y), s4.pt(0.12, 1), BOUNDARY[0], rad=0.04,
          dashed=True)
    arrow(s5.pt(-0.75, 1), (1.55, floor_y), FORCING[0], rad=-0.25,
          dashed=True)
    arrow(s5.pt(1, 0.10), s4.pt(-0.62, 1), CLOUD[0], rad=-0.35, dashed=True)

    note(9.45, 7.45, "dshfl, dlhfl, devap, dshdt,\ndlhdt, dtaux, dtauy, dust3",
         FORCING[1])
    note(10.05, 9.95, "Ts, qs, dz0, dwetfac\n(next-step surface fluxes)",
         BOUNDARY[1])
    note(10.25, 6.30, "dswfl, dlwfl, dflux\n(this step)", FORCING[1])
    note(10.25, 5.58, "dalb, dicec (next step)", BOUNDARY[1])
    note(6.80, 4.75, "dcc, dql\n(next radiation)", CLOUD[1])
    note(3.00, 6.57, "dprl, dprc, dprs", FORCING[1])
    ax.text(0.40, 13.20,
            "This-step tendencies dtdt, dqdt, dudt, dvdt\n"
            "→ spectral integration (applied now)",
            ha="left", va="center", fontsize=FS_TEXT, color=TENDENCY[1],
            linespacing=1.3, zorder=7)

    legend(y_top=0.92)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--outdir", type=Path, default=IMAGES_DIR,
                        help="output directory (default: tutorial images)")
    parser.add_argument("--formats", nargs="+", default=["png"],
                        choices=["png", "svg", "pdf"])
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()

    draw()
    args.outdir.mkdir(parents=True, exist_ok=True)
    for fmt in args.formats:
        out = args.outdir / f"{STEM}.{fmt}"
        fig.savefig(out, dpi=args.dpi, facecolor="white")
        print("saved", out)


if __name__ == "__main__":
    main()
