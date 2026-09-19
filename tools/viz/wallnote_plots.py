"""Dark wall-note plots for the viscometer (eta_r vs phi) and Oberbeck (anisotropy ratio ladder, absolutes).
Numbers are the datasheet's (rows d52_v21_einstein, d52_v22_phi10, d52_v23_phi20, d61_v123_l3, d61_v4_resolution, d61_v5_halfsize, d61_review_corrections)."""
import sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
BG, FG, GRID = "#2a2d31", "#f2f2f2", "#5a5f66"
plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG, "axes.edgecolor": FG,
                     "axes.labelcolor": FG, "xtick.color": FG, "ytick.color": FG, "text.color": FG, "axes.titlecolor": FG,
                     "grid.color": GRID, "grid.alpha": 0.6, "axes.grid": True, "font.size": 12, "axes.labelsize": 13,
                     "legend.facecolor": BG, "legend.edgecolor": FG, "legend.framealpha": 0.9})
outv, outo = sys.argv[1], sys.argv[2]
# ---- viscometer ------------------------------------------------------------
phi = np.linspace(0, 0.25, 300)
ein = 1 + 2.5 * phi
bat = 1 + 2.5 * phi + 6.2 * phi ** 2
phim = 0.64
kd = (1 - phi / phim) ** (-2.5 * phim)
meas = [(0.05, 1.1062), (0.10, 1.2454), (0.20, 1.7143)]
comp = [(0.05, 1.0996, "Einstein composite"), (0.10, 1.2552, "Batchelor composite"), (0.20, 1.7305, "Krieger-Dougherty composite")]
fig, ax = plt.subplots(figsize=(9, 6.5))
ax.plot(phi, ein, color="#8dff6a", lw=2.0, ls="--", label="Einstein, 1 + 2.5 phi")
ax.plot(phi, bat, color="#ffb84a", lw=2.0, ls="--", label="Batchelor, 1 + 2.5 phi + 6.2 phi^2")
ax.plot(phi, kd, color="#ff7ab6", lw=2.0, ls="--", label="Krieger-Dougherty, phi_m = 0.64")
ax.plot([c[0] for c in comp], [c[1] for c in comp], "s", ms=9, mfc="none", mec=FG, mew=1.6, label="closure evaluated on the measured phi(r,z)")
ax.plot([m[0] for m in meas], [m[1] for m in meas], "o", ms=11, mfc="#4fd8ff", mec=BG, mew=1.0, label="DNS, T(phi) / T(0) at the plateau")
for (p, e), (_, c, _) in zip(meas, comp):
    ax.annotate("%+.2f %%" % ((e / c - 1) * 100), (p, e), xytext=(10, -18), textcoords="offset points", fontsize=11, color="#4fd8ff")
ax.set_xlim(0, 0.25); ax.set_ylim(0.95, 2.1); ax.set_xlabel("particle volume fraction phi"); ax.set_ylabel("relative viscosity eta_r")
ax.set_title("Numerical viscometer: suspension viscosity from the torque on the rotating cylinder", fontsize=13)
ax.legend(loc="upper left"); fig.tight_layout(); fig.savefig(outv, dpi=300)
# ---- oberbeck ---------------------------------------------------------------
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5.6))
rungs = ["2b/h = 9.5\nfull size", "2b/h = 19\nfull size", "2b/h = 9.5\nhalf size"]
ratio = [1.16993, 1.17000, 1.14148]; ref = 1.14532
xs = np.arange(3)
a1.axhspan(ref * 0.98, ref * 1.02, color=FG, alpha=0.15, label="gate band, +/- 2 %")
a1.axhline(ref, color=FG, lw=2.0, ls="--", label="Oberbeck Y^A / X^A = 1.14532")
a1.plot(xs, ratio, "o", ms=12, mfc="#4fd8ff", mec=BG, label="DNS, F_perp / F_par")
for x, r in zip(xs, ratio):
    a1.annotate("%+.2f %%" % ((r / ref - 1) * 100), (x, r), xytext=(0, 12), textcoords="offset points", ha="center", fontsize=11, color="#4fd8ff")
a1.set_xticks(xs); a1.set_xticklabels(rungs); a1.set_ylim(1.10, 1.20); a1.set_ylabel("anisotropy ratio"); a1.set_title("Anisotropy ratio: lattice systematic collapses at half size", fontsize=12)
a1.legend(loc="lower left", fontsize=10)
par = [-1.55, -1.82, 1.36]; perp = [0.29, 0.14, 0.70]
w = 0.32
a2.axhspan(-3, 3, color=FG, alpha=0.12, label="gate band, +/- 3 %")
a2.bar(xs - w / 2, par, w, color="#4fd8ff", label="axis parallel to the force")
a2.bar(xs + w / 2, perp, w, color="#ffb84a", label="axis perpendicular to the force")
a2.axhline(0, color=FG, lw=1.0)
a2.set_xticks(xs); a2.set_xticklabels(rungs); a2.set_ylim(-4, 4); a2.set_ylabel("drag deviation from Oberbeck [%] (a_eff-corrected)")
a2.set_title("Absolute drags", fontsize=12); a2.legend(loc="upper left", fontsize=10)
fig.suptitle("Oberbeck anisotropic drag on a prolate spheroid, r_e = 2: closing ratio -0.34 %", fontsize=13)
fig.tight_layout(); fig.savefig(outo, dpi=300); print("wrote plots")
