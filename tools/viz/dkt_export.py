"""DKT stroboscopic hand-off: stage positions, ghost sequence, wall-note plot, preview.

Usage: python3 dkt_export.py RUNDIR OUTDIR
Reads RUNDIR/particle_force.log (time ip fx fy fz tx ty tz px py pz vx vy vz; ip=1 leader, ip=2 trailer).
Units: sphere diameter d = 1, box 6 x 6 x 24, floor z = 0.
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

rundir, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
d = np.loadtxt(os.path.join(rundir, "particle_force.log"), comments="#")
ip = d[:, 1].astype(int)
a, b = d[ip == 1], d[ip == 2]
n = min(len(a), len(b)); a, b = a[:n], b[:n]
t = a[:, 0]
x1, z1, x2, z2 = a[:, 8], a[:, 10], b[:, 8], b[:, 10]
vz1, vz2 = a[:, 13], b[:, 13]
dx, dz = x2 - x1, z2 - z1
sep = np.hypot(dx, dz)
tilt = np.degrees(np.arctan2(np.abs(dx), dz))          # 0 = trailer straight above the leader
tilt_signed = np.degrees(np.arctan2(dx, dz))            # >0: trailer to +x of the leader

def at(T):
    i = int(np.argmin(np.abs(t - T))); return i

# derived event times from the log itself
i_kiss = int(np.argmax(sep <= 1.005))
i_sepon = int(np.argmax((sep > 1.02) & (t > t[i_kiss])))
i_horiz = int(np.argmax((tilt >= 90.0) & (t > t[i_kiss])))
i_vmax_tr = int(np.argmin(vz2[: i_kiss]))              # trailer's fastest drafting speed before the kiss
print("t_kiss=%.2f  separation onset=%.2f  horizontal=%.2f  trailer v_max=%.3f at t=%.2f" %
      (t[i_kiss], t[i_sepon], t[i_horiz], vz2[i_vmax_tr], t[i_vmax_tr]))

stages = [
    ("start",        at(0.0),    "Release",             "Two equal spheres released from rest, the trailer 1.5 d above the leader and offset 0.05 d sideways."),
    ("drafting",     at(10.0),   "Drafting",            "The trailer sits in the leader's wake, feels less drag and falls faster: the gap closes."),
    ("kissing",      i_kiss,     "Kissing",             "Contact. Centres 1 d apart, the pair nearly vertical (tilt %.1f deg)."),
    ("rolling",      at(25.0),   "Rolling contact",     "Still touching, the doublet rolls over: tilt %.1f deg and growing, frictionless contact through the liquid film."),
    ("tumbling",     at(28.0),   "Tumbling",            "The tumble accelerates, tilt %.1f deg; the wake asymmetry that started it now drives it."),
    ("separation",   i_sepon,    "Separation",          "The gap opens (centres %.2f d apart) while the pair keeps rotating."),
    ("horizontal",   i_horiz,    "Side by side",        "The pair passes through horizontal (tilt %.0f deg), %.2f d apart, falling together."),
    ("exchange",     at(40.0),   "Role exchange",       "The former trailer is now below and falling faster (%.3f vs %.3f); the sequence is complete."),
]
rows = []
for key, i, label, desc in stages:
    fill = {
        "kissing": (tilt[i],), "rolling": (tilt[i],), "tumbling": (tilt[i],), "separation": (sep[i],),
        "horizontal": (tilt[i], sep[i]), "exchange": (vz2[i], vz1[i]),
    }.get(key, ())
    rows.append(dict(stage=key, label=label, time=round(float(t[i]), 3),
                     leader=dict(x=round(float(x1[i]), 5), y=0.0, z=round(float(z1[i]), 5), vz=round(float(vz1[i]), 5)),
                     trailer=dict(x=round(float(x2[i]), 5), y=0.0, z=round(float(z2[i]), 5), vz=round(float(vz2[i]), 5)),
                     separation_d=round(float(sep[i]), 4), tilt_deg=round(float(tilt_signed[i]), 2),
                     description=desc % fill if fill else desc))
json.dump(dict(units="sphere diameter d = 1; box x,y in [-3,3], z in [0,24]; floor z = 0; time in box units",
               source=os.path.basename(rundir), stages=rows), open(os.path.join(out, "stages.json"), "w"), indent=2)
with open(os.path.join(out, "stages.csv"), "w") as f:
    f.write("stage,label,time,leader_x,leader_z,trailer_x,trailer_z,separation_d,tilt_deg,leader_vz,trailer_vz\n")
    for r in rows:
        f.write("%s,%s,%.3f,%.5f,%.5f,%.5f,%.5f,%.4f,%.2f,%.5f,%.5f\n" % (r["stage"], r["label"], r["time"], r["leader"]["x"], r["leader"]["z"],
                r["trailer"]["x"], r["trailer"]["z"], r["separation_d"], r["tilt_deg"], r["leader"]["vz"], r["trailer"]["vz"]))

# ghost sequence every 2 t.u. and the full trajectory
with open(os.path.join(out, "ghosts_every_2tu.csv"), "w") as f:
    f.write("time,leader_x,leader_z,trailer_x,trailer_z,separation_d,tilt_deg\n")
    for T in np.arange(0.0, 40.01, 2.0):
        i = at(T); f.write("%.1f,%.5f,%.5f,%.5f,%.5f,%.4f,%.2f\n" % (t[i], x1[i], z1[i], x2[i], z2[i], sep[i], tilt_signed[i]))
with open(os.path.join(out, "trajectory_full.csv"), "w") as f:
    f.write("time,leader_x,leader_z,leader_vz,trailer_x,trailer_z,trailer_vz,separation_d,tilt_deg\n")
    for i in range(0, n):
        f.write("%.3f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.4f,%.2f\n" % (t[i], x1[i], z1[i], vz1[i], x2[i], z2[i], vz2[i], sep[i], tilt_signed[i]))

# ---- dark wall-note plot: descent, separation and tilt, stage markers -------
BG, FG, GRID = "#2a2d31", "#f2f2f2", "#5a5f66"
C_LEAD, C_TRAIL, C_SEP, C_TILT = "#4fd8ff", "#ffb84a", "#8dff6a", "#ff7ab6"
plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG, "axes.edgecolor": FG,
                     "axes.labelcolor": FG, "xtick.color": FG, "ytick.color": FG, "text.color": FG, "axes.titlecolor": FG,
                     "grid.color": GRID, "grid.alpha": 0.6, "axes.grid": True, "font.size": 12, "axes.labelsize": 13,
                     "legend.facecolor": BG, "legend.edgecolor": FG, "legend.framealpha": 0.9})
fig, (p1, p2) = plt.subplots(2, 1, figsize=(11, 8.5), sharex=True)
p1.plot(t, z1, color=C_LEAD, lw=2.4, label="leader z")
p1.plot(t, z2, color=C_TRAIL, lw=2.4, label="trailer z")
p1.set_ylabel("height z / d"); p1.legend(loc="lower left")
p2.plot(t, sep, color=C_SEP, lw=2.4, label="centre distance / d")
p2b = p2.twinx()
p2b.plot(t, tilt, color=C_TILT, lw=2.4, label="tilt from vertical [deg]")
p2b.set_ylabel("tilt [deg]", color=C_TILT); p2b.tick_params(axis="y", colors=C_TILT); p2b.grid(False)
p2.set_ylabel("distance / d", color=C_SEP); p2.tick_params(axis="y", colors=C_SEP); p2.set_xlabel("t")
p2.set_ylim(0.9, 2.6); p2b.set_ylim(0, 120)
for r in rows[2:]:
    for p in (p1, p2):
        p.axvline(r["time"], color=FG, alpha=0.35, lw=1.0, ls="--")
    p1.text(r["time"], p1.get_ylim()[1] * 0.98, " " + r["label"], rotation=90, va="top", ha="right", fontsize=10, color=FG, alpha=0.9)
h1, l1 = p2.get_legend_handles_labels(); h2, l2 = p2b.get_legend_handles_labels()
p2.legend(h1 + h2, l1 + l2, loc="upper left")
fig.suptitle("Drafting, kissing, tumbling: frictionless contact, D/h = 8 (t_kiss = %.2f, separation %.2f, horizontal %.2f)" %
             (t[i_kiss], t[i_sepon], t[i_horiz]), fontsize=13)
fig.tight_layout()
fig.savefig(os.path.join(out, "dkt_phases_plot_dark_300dpi.png"), dpi=300)

# ---- stroboscopic preview (x-z side view, ghosts every 2 t.u., stages labelled) ----
fig, ax = plt.subplots(figsize=(6.5, 13))
ax.set_facecolor(BG); ax.grid(False)
ax.add_patch(plt.Rectangle((-3, 0), 6, 24, fill=False, ec=FG, lw=1.2, alpha=0.6))
ghosts = np.arange(0.0, 40.01, 2.0)
for k, T in enumerate(ghosts):
    i = at(T); al = 0.12 + 0.5 * k / len(ghosts)
    ax.add_patch(Circle((x1[i], z1[i]), 0.5, fc=C_LEAD, ec="none", alpha=al))
    ax.add_patch(Circle((x2[i], z2[i]), 0.5, fc=C_TRAIL, ec="none", alpha=al))
for r in rows:
    zl = 0.5 * (r["leader"]["z"] + r["trailer"]["z"]); xl = max(r["leader"]["x"], r["trailer"]["x"]) + 0.8
    ax.annotate("%s  (t = %.1f)" % (r["label"], r["time"]), xy=(xl - 0.25, zl), xytext=(3.4, zl), color=FG, fontsize=9.5,
                va="center", arrowprops=dict(arrowstyle="-", color=FG, alpha=0.5, lw=0.8))
ax.set_xlim(-3.5, 7.5); ax.set_ylim(-0.5, 24.5); ax.set_aspect("equal")
ax.set_xlabel("x / d"); ax.set_ylabel("z / d")
ax.set_title("Stroboscopic sequence, ghosts every 2 t.u. (cyan leader, amber trailer)", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(out, "strobe_preview_dark.png"), dpi=200)
print("wrote", sorted(os.listdir(out)))
