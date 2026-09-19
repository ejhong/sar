"""Figure styling shared by all experiments (palette from the dataviz reference instance)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
BLUE_STEPS = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
              "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
CMAP_SEQ = LinearSegmentedColormap.from_list("seqblue", BLUE_STEPS)
CMAP_SEQ_R = LinearSegmentedColormap.from_list("seqblue_r", BLUE_STEPS[::-1])
CMAP_GRAY = LinearSegmentedColormap.from_list("paper_gray", ["#1a1a19", "#5a5750", "#a29d92", "#e6e1d6", "#fcfcfb"])
CMAP_DIV = LinearSegmentedColormap.from_list("div", ["#0d366b", "#2a78d6", "#f0efec", "#e34948", "#7a1414"])


def style():
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "font.family": "sans-serif", "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 9.5, "axes.titlesize": 10.5, "axes.labelsize": 9.5, "axes.titleweight": "medium",
        "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "axes.grid": False, "grid.color": GRID, "grid.linewidth": 0.8, "grid.linestyle": "-",
        "lines.linewidth": 2.0, "lines.solid_capstyle": "round", "lines.solid_joinstyle": "round",
        "legend.frameon": False, "legend.fontsize": 8.5,
        "savefig.dpi": 170, "savefig.bbox": "tight", "savefig.pad_inches": 0.08,
        "text.color": INK,
    })


def db(a, floor=-40):
    a = np.abs(a).astype(float)
    a = 20 * np.log10(a / (a.max() + 1e-30) + 1e-30)
    return np.clip(a, floor, 0)


def imshow_slc(ax, slc, x, r, floor=-40, vmin=None, vmax=None, **kw):
    """Amplitude in dB, azimuth on the horizontal axis, slant range increasing downward."""
    im = ax.imshow(db(slc, floor).T, cmap=CMAP_GRAY, aspect="equal", origin="upper",
                   extent=[x[0], x[-1], r[-1], r[0]], vmin=floor if vmin is None else vmin,
                   vmax=0 if vmax is None else vmax, interpolation="nearest", **kw)
    ax.set_xlabel("azimuth (m)")
    ax.set_ylabel("slant range (m)  to  toward radar is up")
    return im


def finish(fig, path):
    fig.savefig(path)
    plt.close(fig)
    return path
