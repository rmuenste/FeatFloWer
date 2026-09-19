"""Hindered settling (D3.2) geometry-only hand-off: cloud snapshot, trails, reference sphere, RZ plot, preview.

Usage: python3 hindered_export.py CLOUD_RUNDIR UT_RUNDIR OUTDIR T_SNAP
Units: sphere diameter d = 1; walled column x,y in [-3,3], z in [0,24]; floor z = 0; u_t = 0.4061 (row d32_ut_ref).
"""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

cloud, utdir, out, T = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
os.makedirs(out, exist_ok=True)
UT = 0.4061
d = np.loadtxt(os.path.join(cloud, "particle_force.log"), comments="#")
t, ip = d[:, 0], d[:, 1].astype(int)
N = ip.max()
def snap(T0):
    m = np.abs(t - T0) < 0.0026
    s = d[m]; s = s[np.argsort(s[:, 1])]
    return s
s = snap(T)
x, y, z, vx, vy, vz = s[:, 8], s[:, 9], s[:, 10], s[:, 11], s[:, 12], s[:, 13]
ratio = -vz / UT
print("N=%d t=%.3f mean vz=%.4f (U/u_t=%.3f) ratio range %.2f..%.2f z %.2f..%.2f" % (N, s[0, 0], vz.mean(), -vz.mean() / UT, ratio.min(), ratio.max(), z.min(), z.max()))
with open(os.path.join(out, "cloud_t%02d.csv" % T), "w") as f:
    f.write("# %s: %d spheres d = 1 at t = %.3f; ratio = -vz/u_t with u_t = %.4f (isolated sphere in the same column)\n" % (os.path.basename(cloud), N, s[0, 0], UT))
    f.write("id,x,y,z,vx,vy,vz,settling_ratio\n")
    for i in range(N):
        f.write("%d,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.4f\n" % (s[i, 1], x[i], y[i], z[i], vx[i], vy[i], vz[i], ratio[i]))
# trails: positions over the 3 t.u. before the snapshot, every 0.25
with open(os.path.join(out, "trails_t%02d.csv" % T), "w") as f:
    f.write("# positions of every sphere at 13 instants from t-3 to t (0.25 apart); draw as short polylines behind each sphere\n")
    f.write("time,id,x,y,z\n")
    for T0 in np.arange(T - 3.0, T + 0.001, 0.25):
        q = snap(T0)
        for i in range(N):
            f.write("%.2f,%d,%.5f,%.5f,%.5f\n" % (q[0, 0], q[i, 1], q[i, 8], q[i, 9], q[i, 10]))
# reference single sphere
u = np.loadtxt(os.path.join(utdir, "particle_force.log"), comments="#")
i = int(np.argmin(np.abs(u[:, 0] - T)))
ref = dict(time=float(u[i, 0]), x=float(u[i, 8]), y=float(u[i, 9]), z=float(u[i, 10]), vz=float(u[i, 13]))
print("reference sphere at t=%.2f: z=%.3f vz=%.4f" % (ref["time"], ref["z"], ref["vz"]))
json.dump(dict(units="sphere diameter d = 1; column x,y in [-3,3], z in [0,24]; floor z = 0",
               cloud_run=os.path.basename(cloud), N=int(N), snapshot_time=float(s[0, 0]),
               mean_settling_ratio=float(-vz.mean() / UT), u_t=UT, reference_sphere=ref,
               cloud_extent=dict(x=[float(x.min()), float(x.max())], y=[float(y.min()), float(y.max())], z=[float(z.min()), float(z.max())])),
          open(os.path.join(out, "snapshot.json"), "w"), indent=2)

# ---- RZ plot (dark): the ladder numbers from tools/d32_ladder_analysis.py --window 15 25 --ut 0.4061
ladder = {20: [(0.0193, 0.895), (0.0188, 0.901), (0.0205, 0.879)], 40: [(0.0453, 0.802), (0.0378, 0.842), (0.0382, 0.834)],
          80: [(0.0790, 0.684), (0.0756, 0.706), (0.0791, 0.678)], 120: [(0.1054, 0.597), (0.1126, 0.561), (0.1147, 0.605)]}
n_fit = 4.58
BG, FG, GRID = "#2a2d31", "#f2f2f2", "#5a5f66"
plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG, "axes.edgecolor": FG,
                     "axes.labelcolor": FG, "xtick.color": FG, "ytick.color": FG, "text.color": FG, "axes.titlecolor": FG,
                     "grid.color": GRID, "grid.alpha": 0.6, "axes.grid": True, "font.size": 12, "axes.labelsize": 13,
                     "legend.facecolor": BG, "legend.edgecolor": FG, "legend.framealpha": 0.9})
cols = {20: "#4fd8ff", 40: "#8dff6a", 80: "#ffb84a", 120: "#ff7ab6"}
fig, ax = plt.subplots(figsize=(9, 6.5))
phi = np.linspace(0.005, 0.14, 200)
ax.fill_between(phi, (1 - phi) ** 2.7, (1 - phi) ** 3.0, color=FG, alpha=0.18, label="unbounded Rowe / Richardson-Zaki band, n = 2.7 to 3.0")
ax.plot(phi, (1 - phi) ** n_fit, color=FG, lw=2.2, ls="--", label="fit through the origin, n = %.2f" % n_fit)
for Nn, pts in ladder.items():
    ax.plot([p[0] for p in pts], [p[1] for p in pts], "o", ms=9, mfc=cols[Nn], mec=BG, mew=1.0, label="N = %d, three seeds" % Nn)
ax.set_xlabel("cloud volume fraction phi"); ax.set_ylabel("U / u_t")
ax.set_xlim(0, 0.14); ax.set_ylim(0.5, 1.0)
ax.set_title("Hindered settling in the 6 d walled column: confined exponent 4.6 vs unbounded 2.7 to 3.0", fontsize=13)
ax.legend(loc="lower left")
fig.tight_layout(); fig.savefig(os.path.join(out, "hindered_rz_plot_dark_300dpi.png"), dpi=300)

# ---- preview: side view of the cloud coloured by settling ratio, trails, reference sphere box
fig, ax = plt.subplots(figsize=(7, 12))
ax.set_facecolor(BG); ax.grid(False)
ax.add_patch(plt.Rectangle((-3, 0), 6, 24, fill=False, ec=FG, lw=1.2, alpha=0.6))
ax.add_patch(plt.Rectangle((4, 0), 6, 24, fill=False, ec=FG, lw=1.2, alpha=0.35))
order = np.argsort(y)      # far spheres first
sc = ax.scatter(x[order], z[order], c=ratio[order], cmap="viridis", vmin=0.0, vmax=1.2, s=170, edgecolors="none")
ax.scatter([4 + 3 + ref["x"]], [ref["z"]], c=[1.0], cmap="viridis", vmin=0.0, vmax=1.2, s=170, edgecolors=FG, linewidths=0.8)
ax.text(7, ref["z"] + 1.2, "isolated sphere\nu_t reference", ha="center", fontsize=9.5, color=FG)
ax.text(0, 23.0, "N = %d cloud, t = %.0f\nmean U / u_t = %.2f" % (N, T, -vz.mean() / UT), ha="center", fontsize=10, color=FG)
cb = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02); cb.set_label("settling speed / u_t"); cb.ax.yaxis.set_tick_params(color=FG)
ax.set_xlim(-3.5, 10.5); ax.set_ylim(-0.5, 24.5); ax.set_aspect("equal"); ax.set_xlabel("x / d"); ax.set_ylabel("z / d")
ax.set_title("Hindered settling snapshot (side view), colour = own settling speed / u_t", fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(out, "cloud_preview_dark.png"), dpi=200)
print("wrote", sorted(os.listdir(out)))
