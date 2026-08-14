import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

R="/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/"
CASES={"fric8":R+"q2p1_dns_rundir_dkt_offset/particle_force.log",
       "fric16":R+"q2p1_dns_rundir_dkt16_offset/particle_force.log",
       "nofric":R+"q2p1_dns_rundir_dkt_nofric_long/particle_force.log"}
def load(p):
    d=np.loadtxt(p,comments="#"); ip=d[:,1].astype(int)
    a,b=d[ip==1],d[ip==2]; n=min(len(a),len(b)); a,b=a[:n],b[:n]
    dx,dz=b[:,8]-a[:,8], b[:,10]-a[:,10]
    return dict(t=a[:,0],x1=a[:,8],z1=a[:,10],x2=b[:,8],z2=b[:,10],
                tilt=np.degrees(np.arctan2(np.abs(dx),dz)),sep=np.hypot(dx,dz),dx=dx,dz=dz)
D={k:load(v) for k,v in CASES.items()}

C8,C16,CNF="#b5473f","#3d6fa8","#0b7d8c"
GRID="#d7dee0"; INK="#1d2b30"; MUT="#5d747c"
plt.rcParams.update({"font.size":10,"axes.edgecolor":"#9fb0b5","axes.labelcolor":INK,
    "text.color":INK,"xtick.color":MUT,"ytick.color":MUT,"axes.titlesize":11})

fig,(axA,axB)=plt.subplots(2,1,figsize=(7.2,8.6),gridspec_kw=dict(height_ratios=[1,1.22],hspace=0.28))

# ---------------- Panel A : tilt vs time
axA.axhline(90,color=MUT,ls=":",lw=1)
axA.text(3.0,91.5,"pair horizontal",color=MUT,fontsize=8.5)
axA.plot(D["fric8"]["t"],D["fric8"]["tilt"],color=C8,lw=1.8,label="friction on, D/h = 8")
axA.plot(D["fric16"]["t"],D["fric16"]["tilt"],color=C16,lw=1.8,label="friction on, D/h = 16")
axA.plot(D["nofric"]["t"],D["nofric"]["tilt"],color=CNF,lw=2.4,label="frictionless, D/h = 8")

for tk,c,lab in ((18.04,C8,None),(18.34,C16,None),(18.04,CNF,"kiss")):
    axA.axvline(tk,color=c,ls="--",lw=0.9,alpha=.55)
axA.annotate("kissing\nt = 18.0 – 18.3",xy=(18.2,7),xytext=(11.0,34),fontsize=8.5,color=MUT,
    ha="center",arrowprops=dict(arrowstyle="->",color=MUT,lw=.9,shrinkB=2))
axA.axvline(29.64,color=CNF,ls="-.",lw=1.1,alpha=.75)
axA.annotate("separation onset\nt = 29.64",xy=(29.64,20),xytext=(34.6,13),fontsize=8.5,color=CNF,
    ha="center",arrowprops=dict(arrowstyle="->",color=CNF,lw=.9,shrinkB=2))
axA.annotate("locked at ~4°",xy=(24.9,4.1),xytext=(24.2,24),fontsize=9,color=C8,ha="right",
    arrowprops=dict(arrowstyle="->",color=C8,lw=1.0,shrinkB=2))
axA.annotate("109°",xy=(40,109.1),xytext=(37.4,97),fontsize=9,color=CNF,ha="center",
    arrowprops=dict(arrowstyle="->",color=CNF,lw=1.0,shrinkB=2))
axA.set_xlim(0,41); axA.set_ylim(0,120)
axA.set_xlabel("time  $t$"); axA.set_ylabel("tilt of the pair axis  [deg]")
axA.set_title("A   Contact friction decides whether the doublet tumbles",loc="left",fontweight="bold")
axA.grid(True,color=GRID,lw=.7); axA.set_axisbelow(True)
axA.legend(loc="upper left",frameon=True,fontsize=9,framealpha=.95)

# ---------------- Panel B : trajectories, frictionless
n=D["nofric"]
axB.plot(n["x1"],n["z1"],color="#b5473f",lw=1.9,label="sphere 1  (initial leader)",zorder=3)
axB.plot(n["x2"],n["z2"],color="#3d6fa8",lw=1.9,label="sphere 2  (initial trailer)",zorder=3)
step=int(len(n["t"])/26)
for i in range(0,len(n["t"]),step):                      # dumbbell connectors
    axB.plot([n["x1"][i],n["x2"][i]],[n["z1"][i],n["z2"][i]],color="#9aa9ad",lw=.8,alpha=.75,zorder=2)
ann=[(18.04,"kiss  t = 18.0",(1.05,13.9),"center"),
     (32.92,"pair horizontal\nt = 32.9",(1.30,10.4),"center"),
     (40.00,"t = 40:  2.37 d apart,  109°",(1.45,7.45),"center")]
for tq,txt,xyt,ha in ann:
    j=np.argmin(abs(n["t"]-tq))
    axB.plot([n["x1"][j],n["x2"][j]],[n["z1"][j],n["z2"][j]],color=INK,lw=2.0,zorder=4)
    axB.annotate(txt,xy=((n["x1"][j]+n["x2"][j])/2,(n["z1"][j]+n["z2"][j])/2),xytext=xyt,
        fontsize=8.5,color=INK,ha=ha,arrowprops=dict(arrowstyle="->",color=INK,lw=.85,shrinkB=4))
axB.scatter([n["x1"][0],n["x2"][0]],[n["z1"][0],n["z2"][0]],s=46,facecolor="white",
            edgecolor=["#b5473f","#3d6fa8"],zorder=5,lw=1.6)
axB.scatter([n["x1"][-1],n["x2"][-1]],[n["z1"][-1],n["z2"][-1]],s=52,
            color=["#b5473f","#3d6fa8"],zorder=5)
axB.annotate("roles exchanged: sphere 2\nis now the lower one,\nand falling faster",
    xy=(n["x2"][-1],n["z2"][-1]),xytext=(1.10,2.35),fontsize=8.8,color=INK,ha="center",
    arrowprops=dict(arrowstyle="->",color=INK,lw=1.0,shrinkB=5))
axB.set_xlabel("$x$  [sphere diameters]"); axB.set_ylabel("$z$  [sphere diameters]")
axB.set_title("B   The frictionless run, drawn as it falls",loc="left",fontweight="bold")
axB.grid(True,color=GRID,lw=.7); axB.set_axisbelow(True)
axB.legend(loc="upper right",frameon=True,fontsize=9,framealpha=.95)
axB.set_xlim(-2.45,2.65); axB.set_ylim(0.9,22.6)

# equal-aspect inset: pair orientation vector
ins=axB.inset_axes([0.035,0.60,0.275,0.345])
ins.plot(n["dx"],n["dz"],color=CNF,lw=1.6)
ins.scatter([n["dx"][0]],[n["dz"][0]],s=18,facecolor="white",edgecolor=CNF,lw=1.3,zorder=4)
ins.scatter([n["dx"][-1]],[n["dz"][-1]],s=22,color=CNF,zorder=4)
ins.axhline(0,color=MUT,lw=.7,ls=":"); ins.axvline(0,color=MUT,lw=.7,ls=":")
ins.set_aspect("equal"); ins.set_title("pair axis, true aspect",fontsize=7.8,color=MUT,pad=3)
ins.set_xlabel("$\\Delta x$",fontsize=7.5,labelpad=1); ins.set_ylabel("$\\Delta z$",fontsize=7.5,labelpad=1)
ins.tick_params(labelsize=6.5,length=2); ins.grid(True,color=GRID,lw=.5); ins.set_axisbelow(True)
ins.set_facecolor("#f7fafa")

fig.suptitle("Drafting–kissing–tumbling: the tumble is unlocked by removing contact friction",
             fontsize=11.5,fontweight="bold",y=0.977)
fig.savefig("/tmp/artsrc/dkt_tumble_fixed.png",dpi=140,bbox_inches="tight",facecolor="white")
print("written")
