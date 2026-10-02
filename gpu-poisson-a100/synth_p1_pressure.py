#!/usr/bin/env python3
"""Synthetic pressure-Poisson system with the structure, block layout and conditioning of
FeatFloWer's Q2/P1 pressure operator (the "C matrix" of the discrete projection), for GPU
solver tests.

Construction = the real thing in miniature: C = B^T Ml^-1 B, with
  * pressure space  P1-DISCONTINUOUS on a uniform hexahedral grid (4 unknowns per cell:
    cell mean and three slopes; the same layout as the code, Create_LinMatStruct);
  * velocity space  triquadratic Q2 on the 27-node lattice, as in the code (Q2/P1-disc is
    inf-sup stable; a Q1 velocity leaves spurious pressure modes and was rejected in testing);
    C then has the 27-cell vertex-neighbourhood pattern, 4 x 27 = 108 entries/row interior,
    the measured 103/row average with boundaries;
  * B = discrete divergence (gradient transpose) assembled exactly with a 2x2x2 Gauss rule;
  * Ml = lumped velocity mass matrix (the code's Ml), so C is a genuine SPD (up to the constant
    pressure mode) discrete Laplacian on the pressure space: lambda_min = O(1) (~pi^2 on the
    unit box with Neumann walls), lambda_max = O(h^-2).
The constant mode is pinned by a small shift (--shift). Values are those of the construction,
not of the code (Q1 vs Q2, no FBM), the spectrum and sparsity are what a solver sees.

Output: CSR with 64-bit row pointers, 32-bit column indices, float64 values (.npz) + the RHS
b = C x_true for a smooth x_true, so a GPU solve can be verified against x_true.

Usage: synth_p1_pressure.py NX NY NZ [--out F.npz] [--shift 1e-6] [--dry-run]
D6.4 tier-S fine level = 160 x 160 x 480 cells (12.288 M cells, 49.15 M pressure unknowns).
Memory while building (scipy, this script): ~45 bytes per nnz of B plus C -> the full size
needs ~0.5 TB of host RAM; build it in slabs or on the GPU host with the --slab option
(TODO) - the reduced sizes (<= 128x128x384) build on a 1 TB node.
"""
import argparse, numpy as np

def budget(ncell, per_row=103.0, idx=4, val=8):
    nu = 4*ncell; nnz = per_row*nu
    return nu, nnz, nnz*(val+idx) + (nu+1)*8

def assemble(nx, ny, nz, shift=1e-6):
    import scipy.sparse as sp
    h = 1.0/nx; ncell = nx*ny*nz
    # Q2 velocity: nodes on the (2nx+1) x (2ny+1) x (2nz+1) lattice (vertices, edge/face mids, centres)
    mx, my, mz = 2*nx+1, 2*ny+1, 2*nz+1; nvert = mx*my*mz
    cid = np.arange(ncell).reshape(nx, ny, nz)
    vid = np.arange(nvert).reshape(mx, my, mz)
    # 3x3x3 Gauss points (exact for the Q2 x P1 integrands)
    g = np.sqrt(3/5); gp1 = np.array([-g, 0.0, g]); gw1 = np.array([5/9, 8/9, 5/9])
    gp = np.array([(a,b,c) for a in gp1 for b in gp1 for c in gp1]); gw = np.array([wa*wb*wc for wa in gw1 for wb in gw1 for wc in gw1])
    # 1-D quadratic Lagrange basis on [-1,1] at nodes -1,0,1 and its derivative
    def l1(i, t):  return [t*(t-1)/2, 1-t*t, t*(t+1)/2][i]
    def dl1(i, t): return [t-0.5, -2*t, t+0.5][i]
    nodes = [(i,j,k) for i in (0,1,2) for j in (0,1,2) for k in (0,1,2)]      # 27 local nodes
    def dN(v, p):
        i,j,k = v
        return np.array([dl1(i,p[0])*l1(j,p[1])*l1(k,p[2]), l1(i,p[0])*dl1(j,p[1])*l1(k,p[2]), l1(i,p[0])*l1(j,p[1])*dl1(k,p[2])]) * (2/h)
    P = np.column_stack([np.ones(len(gp)), gp[:,0], gp[:,1], gp[:,2]])     # P1-disc basis at gp
    wq = gw*(h/2)**3
    Bloc = np.zeros((3, 4, 27))
    for d in range(3):
        for q in range(4):
            for iv, v in enumerate(nodes):
                Bloc[d, q, iv] = sum(wq[gpi]*P[gpi, q]*dN(v, gp[gpi])[d] for gpi in range(len(gp)))
    I, J, K = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    cells = cid.ravel()
    vglob = np.stack([vid[2*I+i, 2*J+j, 2*K+k].ravel() for (i,j,k) in nodes], axis=1)   # ncell x 27
    Bs = []
    for d in range(3):
        rows = np.repeat(4*cells[:,None] + np.arange(4)[None,:], 27, axis=1).ravel()
        cols = np.tile(vglob[:,None,:], (1,4,1)).ravel()
        vals = np.tile(Bloc[d].ravel(), ncell)
        Bs.append(sp.coo_matrix((vals,(rows,cols)), shape=(4*ncell, nvert)).tocsr())
    # lumped Q2 mass: row sums of the consistent mass = integral of each basis function
    # (1-D: nodes -1,0,1 integrate to 1/3, 4/3, 1/3 on [-1,1]); product form, times (h/2)^3
    w1 = np.array([1/3, 4/3, 1/3])
    Ml = np.zeros(nvert)
    for iv,(i,j,k) in enumerate(nodes):
        np.add.at(Ml, vglob[:,iv], w1[i]*w1[j]*w1[k]*(h/2)**3)
    Mli = sp.diags(1.0/Ml)
    C = sum(B @ Mli @ B.T for B in Bs)
    C = (C + C.T)*0.5
    C = C + sp.diags(np.full(4*ncell, shift))
    C = C.tocsr(); C.sum_duplicates(); C.sort_indices()
    return C, h

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("nx",type=int); ap.add_argument("ny",type=int); ap.add_argument("nz",type=int)
    ap.add_argument("--out"); ap.add_argument("--shift",type=float,default=1e-6); ap.add_argument("--dry-run",action="store_true")
    a=ap.parse_args(); ncell=a.nx*a.ny*a.nz
    nu,nnz,by = budget(ncell)
    print(f"cells {ncell:,} unknowns {nu:,} nnz~{nnz:,.0f} CSR(f64,i32) {by/1e9:.1f} GB  vector {nu*8/1e9:.2f} GB")
    if a.dry_run: return
    A,h = assemble(a.nx,a.ny,a.nz, shift=a.shift)
    X,Y,Z=np.meshgrid((np.arange(a.nx)+0.5)/a.nx,(np.arange(a.ny)+0.5)/a.ny,(np.arange(a.nz)+0.5)/a.nz,indexing="ij")
    nu=A.shape[0]; xt=np.zeros(nu); xt[0::4]=(np.cos(np.pi*X)*np.cos(np.pi*Y)*np.cos(np.pi*Z)).ravel(); b=A@xt
    print(f"built: nu {nu:,} nnz {A.nnz:,} ({A.nnz/nu:.1f}/row) max|A-A^T| {abs(A-A.T).max():.1e} CSR(f64,i32) {(A.data.nbytes+4*A.nnz+8*(nu+1))/1e9:.2f} GB")
    out=a.out or f"p1_pressure_{a.nx}x{a.ny}x{a.nz}.npz"
    np.savez(out, indptr=A.indptr.astype(np.int64), indices=A.indices.astype(np.int32), data=A.data.astype(np.float64), b=b, x_true=xt, nx=a.nx, ny=a.ny, nz=a.nz, h=h)
    print("wrote", out)

if __name__=="__main__": main()
