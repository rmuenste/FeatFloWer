"""G4 slab profile from the FeatFloWer VTK output (ASCII vtu pieces; ELEpsilonF/ELAlphaP are CELL
data = the element-wise constant fields the source uses, Velocity is Q1 vertex data):
usage: vft_slab_profile.py <rundir> <index> [<index> ...]
Per z-slab of 0.1 m (10 slabs of the 0.1 x 0.1 x 1 m box): slab mean of eps_e * <u_z>_e (cell eps
times the 8-vertex mean of u_z, volume weighted), slab mean of alpha_p, and the ratio eps u_z / (|w| alpha)
with |w| = 0.3 m/s; plus the box integral / MOM (the crude vertex-average G5)."""
import sys, os, re, glob
run, idx = sys.argv[1], sys.argv[2:]
W, MOM, NS = 0.3, 4.24115e-6, 10
def read_piece(path):
    txt = open(path).read()
    def arr(name, ncomp):
        m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name, txt, re.S)
        vals = [float(x) for x in m.group(1).split()]
        return [vals[i:i+ncomp] for i in range(0, len(vals), ncomp)] if ncomp > 1 else vals
    pts = arr("Points", 3); eps = arr("ELEpsilonF", 1); alp = arr("ELAlphaP", 1); vel = arr("Velocity", 3)
    m = re.search(r'<DataArray[^>]*Name="connectivity"[^>]*>(.*?)</DataArray>', txt, re.S)
    con = [int(x) for x in m.group(1).split()]
    return pts, eps, alp, vel, [con[i:i+8] for i in range(0, len(con), 8)]
for ix in idx:
    pv = os.path.join(run, "_vtk", "main.%s.pvtu" % ix)
    pieces = re.findall(r'Source="([^"]+)"', open(pv).read())
    flux = [0.0]*NS; alpha = [0.0]*NS; vol = [0.0]*NS
    for pc in pieces:
        pts, eps, alp, vel, cells = read_piece(os.path.join(run, "_vtk", pc))
        for ic, c in enumerate(cells):
            xs = [pts[i][0] for i in c]; ys = [pts[i][1] for i in c]; zs = [pts[i][2] for i in c]
            v = (max(xs)-min(xs))*(max(ys)-min(ys))*(max(zs)-min(zs))
            zc = sum(zs)/8.0; k = min(NS-1, max(0, int(zc/0.1)))
            flux[k] += v*eps[ic]*sum(vel[i][2] for i in c)/8.0
            alpha[k] += v*alp[ic]
            vol[k] += v
    tot = sum(flux)
    print(f"{os.path.basename(os.path.normpath(run))} index {ix}: box integral eps u_z / MOM = {tot/MOM:.4f} (cell eps x vertex-mean u_z); V = {sum(vol):.4f}")
    print("  slab  z-range     <eps u_z> [m/s]   <alpha>     eps u_z/(|w| alpha)   share of integral")
    for k in range(NS):
        fm = flux[k]/vol[k]; am = alpha[k]/vol[k]
        r = fm/(W*am) if am > 1e-12 else float("nan")
        print(f"   {k+1:2d}  {0.1*k:.1f}-{0.1*(k+1):.1f}   {fm:+.4e}   {am:.4e}   {r:+8.4f}            {flux[k]/tot:+.4f}")
