"""Dark-background variant of the d62_jeffery_analysis.py plot (same data, same tool, restyled).
Usage: python3 jeffery_plot_dark.py AXIS_LOG OUT.png   (gammadot 0.2, tmin 0.5, seams 40.01,80.01 as in the datasheet row)
"""
import sys, runpy, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
BG, FG, GRID = "#2a2d31", "#f2f2f2", "#5a5f66"
REMAP = {"tab:blue": "#4fd8ff", "tab:red": "#ffb84a", "tab:green": "#8dff6a", "k": FG}
plt.rcParams.update({
    "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
    "axes.edgecolor": FG, "axes.labelcolor": FG, "xtick.color": FG, "ytick.color": FG, "text.color": FG,
    "axes.titlecolor": FG, "grid.color": GRID, "grid.linewidth": 0.8,
    "legend.facecolor": BG, "legend.edgecolor": FG, "legend.framealpha": 0.9,
    "font.size": 12, "axes.labelsize": 13,
})
import matplotlib.axes as ma
_plot = ma.Axes.plot
def plot(self, *a, **k):
    if k.get("color") in REMAP: k["color"] = REMAP[k["color"]]
    if "lw" in k: k["lw"] = 2.2
    if "ms" in k: k["ms"] = 5.0; k["mew"] = 1.3
    return _plot(self, *a, **k)
ma.Axes.plot = plot
_grid = ma.Axes.grid
ma.Axes.grid = lambda self, *a, **k: _grid(self, True, alpha=0.6)
_legend = ma.Axes.legend
def legend(self, *a, **k):
    k = {**k, "fontsize": 11}
    if k.get("ncol") == 4:            # the axis-component panel: put the wide legend above the axes
        k.update(loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False)
    return _legend(self, *a, **k)
ma.Axes.legend = legend
import matplotlib.figure as mf
_sf = mf.Figure.savefig
def sf(self, path, *a, **k):
    k["dpi"] = 300; k["facecolor"] = BG
    return _sf(self, path, *a, **k)
mf.Figure.savefig = sf
_st = mf.Figure.suptitle
mf.Figure.suptitle = lambda self, t, **k: _st(self, t, **{**k, "fontsize": 13})
log, out = sys.argv[1:3]
sys.argv = ["d62_jeffery_analysis.py", log, "--gammadot", "0.2", "--tmin", "0.5", "--seams", "40.01,80.01", "--plot", out]
runpy.run_path("/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tools/d62_jeffery_analysis.py", run_name="__main__")
