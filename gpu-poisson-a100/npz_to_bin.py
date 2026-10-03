#!/usr/bin/env python3
"""Dump a gpu-poisson-a100 instance (.npz from synth_p1_pressure.py / gpu_amg_cg.py --build) to raw
binary files for the C driver hypre_boomer_pcg.c (route 2 of README.md §5: CUDA hypre BoomerAMG).

Written next to PREFIX (keep PREFIX on node-local scratch, the files are as large as the .npz):
  PREFIX.meta.txt     key=value lines: nx ny nz h nu nnz
  PREFIX.indptr.i64   CSR row pointer, int64, nu+1
  PREFIX.indices.i32  CSR column indices, int32, nnz
  PREFIX.data.f64     CSR values, float64, nnz
  PREFIX.b.f64        stored smooth RHS b = C x_true
  PREFIX.xtrue.f64    x_true
  PREFIX.rough.f64    the cpu_reference.py rough RHS (rng(1) standard normal, ones-mean removed)
  PREFIX.cand.f64     4 candidate vectors [e, x, y, z] in the P1-disc basis (gpu_amg_cg.candidates_linear),
                      vector c contiguous at offset c*nu (shape (4, nu) row-major)
Usage: python3 npz_to_bin.py --npz F.npz --out PREFIX [--force]
"""
import argparse, os, sys, time
import numpy as np


def candidates_linear(nx, ny, nz, h):
    """Same as gpu_amg_cg.candidates_linear: mean dof = value at the cell centre, slope dof = h/2 (xi in [-1,1])."""
    nu = 4 * nx * ny * nz
    X, Y, Z = np.meshgrid((np.arange(nx) + 0.5) / nx, (np.arange(ny) + 0.5) / ny, (np.arange(nz) + 0.5) / nz, indexing="ij")
    B = np.zeros((4, nu)); B[0, 0::4] = 1.0
    for c, F in enumerate((X, Y, Z), start=1):
        B[c, 0::4] = F.ravel(); B[c, c::4] = h / 2
    return B


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--npz", required=True); ap.add_argument("--out", required=True, help="output PREFIX (scratch!)")
    ap.add_argument("--force", action="store_true", help="rewrite even if PREFIX.meta.txt exists")
    a = ap.parse_args()
    meta = a.out + ".meta.txt"
    if os.path.exists(meta) and not a.force:
        print(f"{meta} exists, nothing to do (use --force)"); return
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    T0 = time.time()
    d = np.load(a.npz)
    nx, ny, nz, h = int(d["nx"]), int(d["ny"]), int(d["nz"]), float(d["h"])
    indptr = d["indptr"].astype(np.int64, copy=False)
    nu = len(indptr) - 1; nnz = int(indptr[-1])
    print(f"[{time.time()-T0:6.1f}s] {a.npz}: {nx}x{ny}x{nz} nu {nu:,} nnz {nnz:,}", flush=True)

    def dump(name, arr, dtype):
        arr = np.ascontiguousarray(arr, dtype=dtype)
        t = time.time(); p = f"{a.out}.{name}"; arr.tofile(p)
        print(f"[{time.time()-T0:6.1f}s] wrote {p} ({arr.nbytes/1e9:.2f} GB, {time.time()-t:.0f} s)", flush=True)
        del arr

    dump("indptr.i64", indptr, np.int64); del indptr
    dump("indices.i32", d["indices"], np.int32)
    dump("data.f64", d["data"], np.float64)
    dump("b.f64", d["b"], np.float64)
    dump("xtrue.f64", d["x_true"], np.float64)
    one = np.ones(nu); rng = np.random.default_rng(1); b = rng.standard_normal(nu); b -= one * (one @ b) / nu
    dump("rough.f64", b, np.float64); del b, one
    dump("cand.f64", candidates_linear(nx, ny, nz, h), np.float64)
    with open(meta, "w") as f:
        f.write(f"nx={nx}\nny={ny}\nnz={nz}\nh={h!r}\nnu={nu}\nnnz={nnz}\nsource={os.path.abspath(a.npz)}\n")
    print(f"[{time.time()-T0:6.1f}s] wrote {meta}; done", flush=True)


if __name__ == "__main__":
    main()
