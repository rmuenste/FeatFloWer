"""Dark-background variant of the compare_tencate.py plot (same data, same tool, restyled).
Usage: python3 tencate_plot_dark.py RUN_LOG CASE TCREF OUT.png
"""
import sys, runpy, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cycler import cycler
BG, FG, GRID = "#2a2d31", "#f2f2f2", "#5a5f66"
LINE, MARK = "#4fd8ff", "#ffb84a"          # DNS line cyan, PIV markers amber
plt.rcParams.update({
    "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
    "axes.edgecolor": FG, "axes.labelcolor": FG, "xtick.color": FG, "ytick.color": FG, "text.color": FG,
    "axes.titlecolor": FG, "axes.prop_cycle": cycler(color=[LINE, MARK]),
    "lines.linewidth": 2.4, "lines.markersize": 5.5, "lines.markeredgewidth": 0.0,
    "axes.grid": True, "grid.color": GRID, "grid.alpha": 0.6, "grid.linewidth": 0.8,
    "legend.facecolor": BG, "legend.edgecolor": FG, "legend.framealpha": 0.9,
    "font.size": 12, "axes.titlesize": 14, "axes.labelsize": 13,
})
import matplotlib.axes as ma
_plot = ma.Axes.plot
def plot(self, *a, **k):
    if "ms" in k: k["ms"] = 5.5          # the tool hard-codes ms=3 for the PIV markers
    return _plot(self, *a, **k)
ma.Axes.plot = plot
import matplotlib.figure as mf
_sf = mf.Figure.savefig
def sf(self, path, *a, **k):
    k["dpi"] = 300; k["facecolor"] = BG
    self.tight_layout()
    return _sf(self, path, *a, **k)
mf.Figure.savefig = sf
log, case, tcref, out = sys.argv[1:5]
sys.argv = ["compare_tencate.py", log, "--case", case, "--tcref", tcref, "--plot", out]
runpy.run_path("/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tools/compare_tencate.py", run_name="__main__")
