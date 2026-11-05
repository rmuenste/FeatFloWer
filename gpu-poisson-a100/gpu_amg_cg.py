#!/usr/bin/env python3
"""GPU (CuPy) Jacobi-CG and AMG-CG harness for the synthetic P1-disc pressure-Poisson
operator of gpu-poisson-a100/README.md (route 1: CuPy + PyAMG).

What it does
  * loads an instance (.npz from synth_p1_pressure.py) or builds one slab-wise with the
    same construction (C = B^T Ml^-1 B, Q2/P1-disc, 2x2x2... 3x3x3 Gauss, lumped Q2 mass,
    symmetrised, constant mode pinned by --shift) -- the slab builder is verified against
    the shipped instances with --verify-build (bitwise-identical structure, data to 1e-15);
  * puts the fine operator on the device in FULL storage (CSR, f64 values, int32 column
    indices, int64 row pointers -> no 2^31 nnz limit; matvec with cuSPARSE on row slabs of
    < 2^31 nnz, or a warp-per-row CUDA kernel) or in SYMMETRIC storage (upper triangle only,
    y = U x + (U - D)^T x with a warp-per-row kernel and f64 atomics);
  * builds the AMG hierarchy with PyAMG on the host (several variants, see VARIANTS), copies
    the coarse levels + transfer operators to the device in f32, and runs CG (f64) with
    one V-cycle (damped Jacobi, 2 pre/2 post, dense coarse inverse on the device) as the
    preconditioner; plain Jacobi-CG (and 4x4 block-Jacobi-CG) as the baselines;
  * writes one JSON with everything the brief asks for (instance, nu, nnz, storage, library
    versions, device, peak device memory, AMG complexity/levels, iterations and wall time for
    Jacobi-CG and AMG-CG, final relative residual, relative error vs x_true).

Usage examples
  python3 gpu_amg_cg.py --instance instances/p1_40x40x120.npz --storage full --amg sa_csr_spec,geo_lin --out results/gpu_p1_40x40x120.json
  python3 gpu_amg_cg.py --build 128 128 384 --save /scratch/.../p1_128x128x384.npz --storage full --amg geo_lin --out results/gpu_p1_128x128x384.json
  python3 gpu_amg_cg.py --instance F.npz --storage sym --validate-sym ...
  --cpu-only runs the same code with numpy/scipy (login-node logic test on the small instance).
"""
import argparse, json, os, sys, time, threading, platform, resource, math
import numpy as np
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
I32MAX = 2**31 - 1
SLAB_NNZ = 2**30          # cuSPARSE SpMV returned CUSPARSE_STATUS_INTERNAL_ERROR on a 2.1e9-nnz int32 slab; 2^30 per slab is safe

# ----------------------------------------------------------------------------------------
# instance construction (same construction as synth_p1_pressure.assemble, slab-wise in x)
# ----------------------------------------------------------------------------------------
def _element_data(h):
    g = np.sqrt(3 / 5); gp1 = np.array([-g, 0.0, g]); gw1 = np.array([5 / 9, 8 / 9, 5 / 9])
    gp = np.array([(a, b, c) for a in gp1 for b in gp1 for c in gp1])
    gw = np.array([wa * wb * wc for wa in gw1 for wb in gw1 for wc in gw1])
    def l1(i, t):  return [t * (t - 1) / 2, 1 - t * t, t * (t + 1) / 2][i]
    def dl1(i, t): return [t - 0.5, -2 * t, t + 0.5][i]
    nodes = [(i, j, k) for i in (0, 1, 2) for j in (0, 1, 2) for k in (0, 1, 2)]
    def dN(v, p):
        i, j, k = v
        return np.array([dl1(i, p[0]) * l1(j, p[1]) * l1(k, p[2]),
                         l1(i, p[0]) * dl1(j, p[1]) * l1(k, p[2]),
                         l1(i, p[0]) * l1(j, p[1]) * dl1(k, p[2])]) * (2 / h)
    P = np.column_stack([np.ones(len(gp)), gp[:, 0], gp[:, 1], gp[:, 2]])
    wq = gw * (h / 2) ** 3
    Bloc = np.zeros((3, 4, 27))
    for d in range(3):
        for q in range(4):
            for iv, v in enumerate(nodes):
                Bloc[d, q, iv] = sum(wq[gpi] * P[gpi, q] * dN(v, gp[gpi])[d] for gpi in range(len(gp)))
    w1 = np.array([1 / 3, 4 / 3, 1 / 3])
    wm = np.array([w1[i] * w1[j] * w1[k] * (h / 2) ** 3 for (i, j, k) in nodes])
    return nodes, Bloc, wm


def assemble_box(nxl, ny, nz, h):
    """C (symmetrised, no shift) of an nxl x ny x nz box of cells of size h; local ids,
    cell id = (i*ny + j)*nz + k, velocity on the (2nxl+1)(2ny+1)(2nz+1) Q2 lattice."""
    nodes, Bloc, wm = _element_data(h)
    ncell = nxl * ny * nz
    mx, my, mz = 2 * nxl + 1, 2 * ny + 1, 2 * nz + 1; nvert = mx * my * mz
    vid = np.arange(nvert).reshape(mx, my, mz)
    I, J, K = np.meshgrid(np.arange(nxl), np.arange(ny), np.arange(nz), indexing="ij")
    cells = np.arange(ncell)
    vglob = np.stack([vid[2 * I + i, 2 * J + j, 2 * K + k].ravel() for (i, j, k) in nodes], axis=1)
    Ml = np.zeros(nvert)
    for iv in range(27):
        np.add.at(Ml, vglob[:, iv], wm[iv])
    Mli = sp.diags(1.0 / Ml)
    rows = np.repeat(4 * cells[:, None] + np.arange(4)[None, :], 27, axis=1).ravel()
    cols = np.tile(vglob[:, None, :], (1, 4, 1)).ravel()
    C = None
    for d in range(3):
        B = sp.coo_matrix((np.tile(Bloc[d].ravel(), ncell), (rows, cols)), shape=(4 * ncell, nvert)).tocsr()
        Cd = B @ Mli @ B.T
        C = Cd if C is None else C + Cd
        del B, Cd
    C = (C + C.T) * 0.5
    return C.tocsr()


def x_true_of(nx, ny, nz):
    X, Y, Z = np.meshgrid((np.arange(nx) + 0.5) / nx, (np.arange(ny) + 0.5) / ny, (np.arange(nz) + 0.5) / nz, indexing="ij")
    xt = np.zeros(4 * nx * ny * nz); xt[0::4] = (np.cos(np.pi * X) * np.cos(np.pi * Y) * np.cos(np.pi * Z)).ravel()
    return xt


def build_slab(args):
    """Rows of the GLOBAL matrix for cells i in [i0, i1): assemble [i0-1, i1+1) (one-cell halo
    gives the exact lumped mass at every node of the slab's cells), take the slab rows,
    add the shift, offset the columns to global ids. Returns (indptr i64, indices i32, data, b_rows)."""
    nx, ny, nz, i0, i1, shift, outfile = args
    h = 1.0 / nx
    a0, a1 = max(i0 - 1, 0), min(i1 + 1, nx)
    C = assemble_box(a1 - a0, ny, nz, h)
    pc = ny * nz                       # cells per x-plane
    r0, r1 = 4 * (i0 - a0) * pc, 4 * (i1 - a0) * pc
    Cs = C[r0:r1]                       # slab rows, local columns
    Cs = (Cs + sp.diags(np.full(Cs.shape[0], shift), offsets=r0, shape=Cs.shape)).tocsr()
    Cs.sum_duplicates(); Cs.sort_indices()
    col_off = 4 * a0 * pc
    indices = (Cs.indices.astype(np.int64) + col_off).astype(np.int32)
    xt = x_true_of(nx, ny, nz)
    Cg = sp.csr_matrix((Cs.data, indices, Cs.indptr), shape=(Cs.shape[0], 4 * nx * pc))
    b = Cg @ xt
    out = (Cs.indptr.astype(np.int64), indices, Cs.data.astype(np.float64), b)
    if outfile:
        np.savez(outfile, indptr=out[0], indices=out[1], data=out[2], b=out[3])
        return outfile
    return out


def build_instance(nx, ny, nz, shift=1e-6, workers=6, slab_unknowns=3.0e6, tmpdir=None, log=print):
    """Slab-wise build of the whole instance (multiprocessing). Returns dict like np.load of an npz."""
    import multiprocessing as mp
    pc = ny * nz
    w = max(1, int(slab_unknowns / (4 * pc)) - 2)
    bounds = [(i, min(i + w, nx)) for i in range(0, nx, w)]
    log(f"[build] {nx}x{ny}x{nz}: {len(bounds)} slabs of <= {w} planes (+halo), {workers} workers")
    t0 = time.time()
    tasks = [(nx, ny, nz, i0, i1, shift, (os.path.join(tmpdir, f"slab_{s:03d}.npz") if tmpdir else None)) for s, (i0, i1) in enumerate(bounds)]
    if workers > 1:
        with mp.get_context("fork").Pool(workers) as pool:
            parts = pool.map(build_slab, tasks, chunksize=1)
    else:
        parts = [build_slab(t) for t in tasks]
    log(f"[build] slabs done in {time.time()-t0:.0f} s, concatenating")
    if tmpdir:
        loaded = [np.load(p) for p in parts]
        parts = [(d["indptr"], d["indices"], d["data"], d["b"]) for d in loaded]
    nu = 4 * nx * pc
    indptr = np.zeros(nu + 1, dtype=np.int64); off = 0; r = 0
    for (ip, ix, da, bb) in parts:
        n = len(ip) - 1; indptr[r + 1:r + n + 1] = ip[1:] + off; off += ip[-1]; r += n
    indices = np.concatenate([p[1] for p in parts]); data = np.concatenate([p[2] for p in parts]); b = np.concatenate([p[3] for p in parts])
    if tmpdir:
        for p in tasks:
            try: os.remove(p[-1])
            except OSError: pass
    log(f"[build] done: nu {nu:,} nnz {len(indices):,} ({len(indices)/nu:.1f}/row) in {time.time()-t0:.0f} s, maxrss {maxrss_gb():.1f} GB")
    return dict(indptr=indptr, indices=indices, data=data, b=b, x_true=x_true_of(nx, ny, nz), nx=nx, ny=ny, nz=nz, h=1.0 / nx)


def maxrss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6   # kB -> GB (decimal)


# ----------------------------------------------------------------------------------------
# backend
# ----------------------------------------------------------------------------------------
class Backend:
    def __init__(self, gpu):
        self.gpu = gpu
        if gpu:
            import cupy as cp, cupyx.scipy.sparse as csp, cupyx.cusparse as cusparse
            self.cp, self.csp, self.cusparse = cp, csp, cusparse
            self.xp = cp
            self.pool = cp.get_default_memory_pool()
        else:
            self.xp = np; self.cp = None
    def asarray(self, a, dtype=None):
        a = np.asarray(a) if dtype is None else np.asarray(a, dtype=dtype)
        return self.cp.asarray(a) if self.gpu else a
    def asnumpy(self, a):
        return self.cp.asnumpy(a) if self.gpu else np.asarray(a)
    def sync(self):
        if self.gpu: self.cp.cuda.Device().synchronize()
    def csr(self, A, dtype):
        """scipy csr -> device csr of given value dtype (canonical, int32 indices)."""
        A = sp.csr_matrix(A); A.sum_duplicates(); A.sort_indices()
        if not self.gpu:
            return A.astype(dtype)
        cp = self.cp
        M = self.csp.csr_matrix((cp.asarray(A.data.astype(dtype)), cp.asarray(A.indices.astype(np.int32)), cp.asarray(A.indptr.astype(np.int32))), shape=A.shape)
        M._has_canonical_format = True
        return M
    def dev_used_gb(self):
        if not self.gpu: return 0.0
        free, total = self.cp.cuda.runtime.memGetInfo(); return (total - free) / 1e9
    def pool_gb(self):
        return self.pool.total_bytes() / 1e9 if self.gpu else 0.0


class DevMemMonitor(threading.Thread):
    """Samples the device's used memory (driver view: total - free) every `dt` s."""
    def __init__(self, be, dt=0.25):
        super().__init__(daemon=True); self.be = be; self.dt = dt; self.peak = 0.0; self.stop_flag = False
    def run(self):
        if not self.be.gpu: return
        self.be.cp.cuda.Device(0).use()
        while not self.stop_flag:
            self.peak = max(self.peak, self.be.dev_used_gb()); time.sleep(self.dt)
    def stop(self):
        self.stop_flag = True; self.join(timeout=5); return self.peak


# ----------------------------------------------------------------------------------------
# fine-level operators
# ----------------------------------------------------------------------------------------
KERNEL_SRC = r"""
extern "C" {
__global__ void csr_spmv_warp(const long long* __restrict__ indptr, const int* __restrict__ indices,
                              const double* __restrict__ data, const double* __restrict__ x,
                              double* __restrict__ y, const int nrows)
{
    const int warp = (blockIdx.x * blockDim.x + threadIdx.x) >> 5;
    const int lane = threadIdx.x & 31;
    if (warp >= nrows) return;
    const long long p0 = indptr[warp], p1 = indptr[warp + 1];
    double s = 0.0;
    for (long long p = p0 + lane; p < p1; p += 32) s += data[p] * x[indices[p]];
    for (int off = 16; off > 0; off >>= 1) s += __shfl_down_sync(0xffffffffu, s, off);
    if (lane == 0) y[warp] = s;
}
/* y += U x + (U - D)^T x  for the upper triangle U (incl. diagonal) in CSR; y must be zeroed. */
__global__ void sym_spmv_warp(const long long* __restrict__ indptr, const int* __restrict__ indices,
                              const double* __restrict__ data, const double* __restrict__ x,
                              double* __restrict__ y, const int nrows)
{
    const int warp = (blockIdx.x * blockDim.x + threadIdx.x) >> 5;
    const int lane = threadIdx.x & 31;
    if (warp >= nrows) return;
    const long long p0 = indptr[warp], p1 = indptr[warp + 1];
    const double xr = x[warp];
    double s = 0.0;
    for (long long p = p0 + lane; p < p1; p += 32) {
        const int c = indices[p]; const double a = data[p];
        s += a * x[c];
        if (c != warp) atomicAdd(&y[c], a * xr);
    }
    for (int off = 16; off > 0; off >>= 1) s += __shfl_down_sync(0xffffffffu, s, off);
    if (lane == 0) atomicAdd(&y[warp], s);
}
}
"""


def upper_triangle(indptr, indices, data, chunk_rows=2_000_000):
    """Upper triangle (col >= row) of a CSR given as arrays; chunked to bound host memory."""
    nu = len(indptr) - 1
    parts_i, parts_d, counts = [], [], np.zeros(nu, dtype=np.int64)
    for r0 in range(0, nu, chunk_rows):
        r1 = min(r0 + chunk_rows, nu)
        p0, p1 = indptr[r0], indptr[r1]
        sub_i = indices[p0:p1]; sub_d = data[p0:p1]
        rows = np.repeat(np.arange(r0, r1, dtype=np.int32), np.diff(indptr[r0:r1 + 1]))
        m = sub_i >= rows
        parts_i.append(sub_i[m]); parts_d.append(sub_d[m])
        counts[r0:r1] = np.bincount(rows[m] - r0, minlength=r1 - r0)
    u_indptr = np.zeros(nu + 1, dtype=np.int64); np.cumsum(counts, out=u_indptr[1:])
    return u_indptr, np.concatenate(parts_i), np.concatenate(parts_d)


class FineOperator:
    """SPD fine-level operator on the device. storage 'full' (matvec 'cusparse' on row slabs
    with <= 2^30 nnz each, or 'kernel') or 'sym' (upper triangle, kernel with atomics)."""
    def __init__(self, be, indptr, indices, data, storage="full", matvec="cusparse", log=print):
        self.be, self.xp, self.storage, self.mode = be, be.xp, storage, matvec
        self.n = len(indptr) - 1; self.nnz_full = int(indptr[-1])
        t0 = time.time()
        # diagonal (host side, cheap in chunks)
        diag = np.empty(self.n)
        for r0 in range(0, self.n, 2_000_000):
            r1 = min(r0 + 2_000_000, self.n); p0, p1 = indptr[r0], indptr[r1]
            rows = np.repeat(np.arange(r0, r1, dtype=np.int32), np.diff(indptr[r0:r1 + 1]))
            m = indices[p0:p1] == rows; diag[r0:r1] = data[p0:p1][m]      # one diagonal entry per row
        self.diag_host = diag
        self.diag = be.asarray(diag)
        if storage == "sym":
            ui, uj, ud = upper_triangle(indptr, indices, data)
            self.nnz_stored = len(uj)
            self.indptr = be.asarray(ui); self.indices = be.asarray(uj); self.data = be.asarray(ud)
            del ui, uj, ud
            if be.gpu:
                self.kern = be.cp.RawModule(code=KERNEL_SRC).get_function("sym_spmv_warp")
            else:
                self.U = sp.csr_matrix((self.data, self.indices, self.indptr), shape=(self.n, self.n))
        else:
            self.nnz_stored = self.nnz_full
            if be.gpu and matvec == "cusparse":
                self.slabs = []
                r0 = 0
                while r0 < self.n:
                    # largest r1 with indptr[r1]-indptr[r0] <= I32MAX
                    r1 = int(np.searchsorted(indptr, indptr[r0] + SLAB_NNZ, side="right")) - 1
                    r1 = min(max(r1, r0 + 1), self.n)
                    p0, p1 = indptr[r0], indptr[r1]
                    cp = be.cp
                    M = be.csp.csr_matrix((cp.asarray(data[p0:p1]), cp.asarray(indices[p0:p1]), cp.asarray((indptr[r0:r1 + 1] - p0).astype(np.int32))), shape=(r1 - r0, self.n))
                    M._has_canonical_format = True
                    self.slabs.append((r0, r1, M))
                    r0 = r1
                log(f"[op] full storage, cuSPARSE, {len(self.slabs)} row slab(s)")
            else:
                self.indptr = be.asarray(indptr); self.indices = be.asarray(indices); self.data = be.asarray(data)
                if be.gpu:
                    self.kern = be.cp.RawModule(code=KERNEL_SRC).get_function("csr_spmv_warp")
                else:
                    self.A = sp.csr_matrix((self.data, self.indices, self.indptr), shape=(self.n, self.n))
        be.sync()
        self.bytes = self.nnz_stored * 12 + (self.n + 1) * 8
        log(f"[op] storage={storage} matvec={matvec} nnz_stored={self.nnz_stored:,} ({self.bytes/1e9:.1f} GB), upload {time.time()-t0:.1f} s")
        self._tmp = None

    def matvec(self, x, out=None):
        xp = self.xp
        if out is None: out = xp.empty_like(x)
        if self.storage == "sym":
            if self.be.gpu:
                out.fill(0.0)
                nthreads = 256; nblocks = (self.n * 32 + nthreads - 1) // nthreads
                self.kern((nblocks,), (nthreads,), (self.indptr, self.indices, self.data, x, out, np.int32(self.n)))
            else:
                out[:] = self.U @ x + self.U.T @ x - self.diag * x
        elif self.be.gpu and self.mode == "cusparse":
            for (r0, r1, M) in self.slabs:
                y = self.be.cusparse.spmv(M, x)
                out[r0:r1] = y
        elif self.be.gpu:
            nthreads = 256; nblocks = (self.n * 32 + nthreads - 1) // nthreads
            self.kern((nblocks,), (nthreads,), (self.indptr, self.indices, self.data, x, out, np.int32(self.n)))
        else:
            out[:] = self.A @ x
        return out

    def time_matvec(self, reps=20):
        x = self.xp.ones(self.n); y = self.xp.empty_like(x)
        self.matvec(x, y); self.be.sync(); t0 = time.time()
        for _ in range(reps): self.matvec(x, y)
        self.be.sync(); return (time.time() - t0) / reps

    def free(self):
        for a in ("slabs", "indptr", "indices", "data", "diag", "U", "A"):
            if hasattr(self, a): delattr(self, a)
        if self.be.gpu: self.be.pool.free_all_blocks()


# ----------------------------------------------------------------------------------------
# AMG hierarchy (PyAMG on the host) -> device V-cycle
# ----------------------------------------------------------------------------------------
VARIANTS = {
    "sa_csr_spec":   "README spec: SA on the dof-level CSR, B=ones (pyamg default), strength=symmetric, aggregate=standard, smooth=jacobi, pyamg default improve_candidates",
    "sa_csr_mean":   "SA on the dof-level CSR with the TRUE near-nullspace B=e (1 on cell means, 0 on slopes), no candidate improvement",
    "sa_bsr_mean":   "SA on the 4x4 cell blocks (BSR), B=e",
    "sa_bsr_lin":    "SA on the 4x4 cell blocks, B=[1,x,y,z] in the P1-disc basis (coarse space reproduces linears)",
    "geo_lin":       "cell-level 2x2x2 geometric aggregation (predefined), B=[1,x,y,z], unsmoothed P (= nested P1-disc coarse spaces, Galerkin)",
    "geo_lin_jac":   "as geo_lin with Jacobi-smoothed prolongation",
}


def candidates_linear(nx, ny, nz, h):
    """[1, x, y, z] in the P1-disc basis: mean dof = value at the cell centre, slope dof = h/2 (xi in [-1,1])."""
    nu = 4 * nx * ny * nz
    X, Y, Z = np.meshgrid((np.arange(nx) + 0.5) / nx, (np.arange(ny) + 0.5) / ny, (np.arange(nz) + 0.5) / nz, indexing="ij")
    B = np.zeros((nu, 4)); B[0::4, 0] = 1.0
    for c, F in enumerate((X, Y, Z), start=1):
        B[0::4, c] = F.ravel(); B[c::4, c] = h / 2
    return B


def geometric_aggops(nx, ny, nz, max_coarse_cells):
    """List of 2x2x2 cell aggregation operators (CSR, cells x aggregates) down to <= max_coarse_cells cells."""
    ops, dims = [], (nx, ny, nz)
    while dims[0] * dims[1] * dims[2] > max_coarse_cells and max(dims) > 1:
        cx, cy, cz = [(d + 1) // 2 for d in dims]
        I, J, K = np.meshgrid(np.arange(dims[0]), np.arange(dims[1]), np.arange(dims[2]), indexing="ij")
        agg = ((I // 2) * cy + (J // 2)) * cz + (K // 2)
        n = dims[0] * dims[1] * dims[2]
        ops.append(sp.csr_matrix((np.ones(n), agg.ravel(), np.arange(n + 1)), shape=(n, cx * cy * cz)))
        dims = (cx, cy, cz)
    return ops, dims


def build_pyamg(A_csr_host, variant, grid, h, max_coarse=1000, log=print):
    import pyamg
    nu = A_csr_host.shape[0]; nx, ny, nz = grid
    e = np.zeros(nu); e[0::4] = 1.0
    t0 = time.time()
    kw = dict(strength="symmetric", aggregate="standard", smooth="jacobi", max_coarse=max_coarse, max_levels=20)
    if variant == "sa_csr_spec":
        ml = pyamg.smoothed_aggregation_solver(A_csr_host, **kw)
    elif variant == "sa_csr_mean":
        ml = pyamg.smoothed_aggregation_solver(A_csr_host, B=e[:, None], improve_candidates=None, **kw)
    else:
        tb = time.time(); A_bsr = A_csr_host.tobsr(blocksize=(4, 4)); log(f"[amg] tobsr {time.time()-tb:.1f} s, blocks {A_bsr.nnz//16:,}")
        if variant == "sa_bsr_mean":
            ml = pyamg.smoothed_aggregation_solver(A_bsr, B=e[:, None], improve_candidates=None, **kw)
        elif variant == "sa_bsr_lin":
            ml = pyamg.smoothed_aggregation_solver(A_bsr, B=candidates_linear(nx, ny, nz, h), improve_candidates=None, **kw)
        elif variant in ("geo_lin", "geo_lin_jac"):
            ops, dims = geometric_aggops(nx, ny, nz, max_coarse_cells=max_coarse // 4)
            log(f"[amg] geometric levels: {grid} -> ... -> {dims} ({len(ops)} coarsenings)")
            ml = pyamg.smoothed_aggregation_solver(A_bsr, B=candidates_linear(nx, ny, nz, h), improve_candidates=None,
                                                   strength="symmetric", aggregate=[("predefined", {"AggOp": op}) for op in ops],
                                                   smooth=("jacobi" if variant == "geo_lin_jac" else None), max_coarse=0, max_levels=len(ops) + 1)
        else:
            raise ValueError(variant)
        del A_bsr
    setup = time.time() - t0
    nnz = [int(l.A.nnz) for l in ml.levels]; n = [int(l.A.shape[0]) for l in ml.levels]
    info = dict(variant=variant, description=VARIANTS[variant], levels=len(ml.levels), sizes=n, nnz_per_level=nnz,
                operator_complexity=float(sum(nnz) / nnz[0]), grid_complexity=float(sum(n) / n[0]), host_setup_s=setup)
    log(f"[amg] {variant}: {len(ml.levels)} levels, n={n}, OC={info['operator_complexity']:.3f}, GC={info['grid_complexity']:.3f}, setup {setup:.1f} s, maxrss {maxrss_gb():.1f} GB")
    return ml, info


class DeviceVCycle:
    """V-cycle preconditioner: level 0 = the f64 FineOperator (smoother in f64), coarse levels
    in f32 (CSR), transfer operators in f32, damped Jacobi (2 pre / 2 post), dense coarse inverse."""
    def __init__(self, be, fine_op, ml, omega_factor=4.0 / 3.0, nu1=2, nu2=2, log=print):
        self.be, self.xp, self.A0, self.nu1, self.nu2 = be, be.xp, fine_op, nu1, nu2
        xp = self.xp
        t0 = time.time()
        self.L = len(ml.levels)
        self.A, self.Dinv, self.P, self.R, self.omega = [None], [None], [], [], []
        for l, lev in enumerate(ml.levels[:-1]):
            self.P.append(be.csr(lev.P, np.float32)); self.R.append(be.csr(lev.R, np.float32))
        for lev in ml.levels[1:-1]:
            A = be.csr(lev.A, np.float32); self.A.append(A); self.Dinv.append(1.0 / A.diagonal())
        Ac = sp.csr_matrix(ml.levels[-1].A).toarray()
        self.Ainv = be.asarray(np.linalg.inv(Ac))
        self.nc = Ac.shape[0]
        # fine level: D^-1 in f64; damping omega_factor / rho(D^-1 A), rho by power iteration
        self.Dinv[0] = 1.0 / fine_op.diag
        self.omega.append(omega_factor / self._rho(lambda v, out: fine_op.matvec(v, out), self.Dinv[0], fine_op.n, np.float64))
        for l in range(1, self.L - 1):
            A = self.A[l]
            self.omega.append(omega_factor / self._rho(lambda v, out, A=A: self._set(out, A @ v), self.Dinv[l], A.shape[0], np.float32))
        be.sync()
        self.bytes = sum(int(M.data.nbytes + M.indices.nbytes + M.indptr.nbytes) for M in self.P + self.R + self.A[1:]) + self.Ainv.nbytes + sum(int(d.nbytes) for d in self.Dinv)
        self.setup_s = time.time() - t0
        log(f"[vcycle] {self.L} levels on device, {self.bytes/1e9:.2f} GB, omega={[round(float(w),3) for w in self.omega]}, coarse n={self.nc}, setup {self.setup_s:.1f} s")

    @staticmethod
    def _set(out, v):
        out[:] = v; return out

    def _rho(self, mv, Dinv, n, dtype, iters=25):
        xp = self.xp
        v = xp.ones(n, dtype=dtype) + 0.01 * xp.asarray(np.random.default_rng(3).standard_normal(n).astype(dtype))
        v /= xp.linalg.norm(v); w = xp.empty_like(v); rho = 0.0
        for _ in range(iters):
            mv(v, w); w *= Dinv
            rho = float(xp.dot(v, w)); nrm = float(xp.linalg.norm(w)); v = w / nrm; w = xp.empty_like(v)
        return 1.02 * max(rho, nrm)

    def _smooth(self, l, A_mv, x, b, sweeps, r):
        Dinv, om = self.Dinv[l], self.omega[l]
        for _ in range(sweeps):
            A_mv(x, r); r *= -1.0; r += b          # r = b - A x
            r *= Dinv; x += om * r
        return x

    def _level_mv(self, l):
        if l == 0: return lambda v, out: self.A0.matvec(v, out)
        A = self.A[l]; return lambda v, out: self._set(out, A @ v)

    def solve_level(self, l, b):
        xp = self.xp
        if l == self.L - 1:
            return self.Ainv @ b.astype(np.float64) if b.dtype != np.float64 else self.Ainv @ b
        mv = self._level_mv(l)
        r = xp.empty_like(b)
        x = self.omega[l] * (self.Dinv[l] * b)                       # first sweep from x=0
        x = self._smooth(l, mv, x, b, self.nu1 - 1, r)
        mv(x, r); r *= -1.0; r += b
        rc = self.R[l] @ (r if r.dtype == np.float32 else r.astype(np.float32))
        ec = self.solve_level(l + 1, rc)
        corr = self.P[l] @ ec.astype(np.float32)
        x += corr if x.dtype == np.float32 else corr.astype(np.float64)
        x = self._smooth(l, mv, x, b, self.nu2, r)
        return x

    def __call__(self, r):
        return self.solve_level(0, r)


# ----------------------------------------------------------------------------------------
# Krylov
# ----------------------------------------------------------------------------------------
def pcg(be, A_mv, b, M, tol=1e-8, maxiter=1000, flexible=True, log=print, name="cg", log_every=0):
    """Preconditioned CG, f64, relative residual ||r||/||b|| <= tol. Returns (x, info)."""
    xp = be.xp
    n = b.shape[0]; x = xp.zeros(n); r = b.copy(); q = xp.empty_like(b)
    bn = float(xp.linalg.norm(b)); z = M(r); p = z.copy(); rz = float(xp.dot(r, z))
    hist = []; it = 0; conv = False
    be.sync(); t0 = time.time()
    for it in range(1, maxiter + 1):
        A_mv(p, q); alpha = rz / float(xp.dot(p, q))
        x += alpha * p; r -= alpha * q
        rn = float(xp.linalg.norm(r)); hist.append(rn / bn)
        if log_every and it % log_every == 0: log(f"  [{name}] it {it} relres {rn/bn:.3e} ({time.time()-t0:.1f} s)")
        if rn <= tol * bn: conv = True; break
        if flexible: r_old = r.copy()
        z = M(r)
        rz_new = float(xp.dot(r, z))
        if flexible:
            beta = float(xp.dot(z, r - r_old)) / rz if rz != 0 else 0.0
            beta = max(beta, 0.0)
        else:
            beta = rz_new / rz
        rz = rz_new
        p = z + beta * p
    be.sync(); dt = time.time() - t0
    A_mv(x, q); true_res = float(xp.linalg.norm(b - q)) / bn
    return x, dict(iters=it, converged=conv, seconds=dt, seconds_per_iter=dt / max(it, 1), final_relres=hist[-1] if hist else None,
                   true_relres=true_res, hist_every_10=[float(h) for h in hist[::10]][:400])


# ----------------------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--instance", help=".npz instance (load; if --build is also given and the file exists, it is loaded instead of rebuilt)")
    ap.add_argument("--build", nargs=3, type=int, metavar=("NX", "NY", "NZ"), help="build the instance slab-wise (same construction as synth_p1_pressure.py)")
    ap.add_argument("--save", help="write the built instance to this .npz (scratch!)")
    ap.add_argument("--shift", type=float, default=1e-6)
    ap.add_argument("--workers", type=int, default=6, help="build workers (each ~16 GB host RAM)")
    ap.add_argument("--tmpdir", help="scratch dir for per-slab files during the build (default: in-memory)")
    ap.add_argument("--verify-build", action="store_true", help="rebuild the loaded instance slab-wise and compare (small instances)")
    ap.add_argument("--storage", choices=["full", "sym"], default="full")
    ap.add_argument("--matvec", choices=["cusparse", "kernel"], default="cusparse", help="full-storage matvec")
    ap.add_argument("--validate-sym", action="store_true", help="compare the symmetric matvec with the full one on this instance")
    ap.add_argument("--amg", default="geo_lin", help="comma list of AMG variants (first = primary); 'none' to skip. Known: " + ", ".join(VARIANTS))
    ap.add_argument("--max-coarse", type=int, default=1000)
    ap.add_argument("--omega-factor", type=float, default=4.0 / 3.0, help="Jacobi damping = factor / rho(D^-1 A)")
    ap.add_argument("--sweeps", type=int, nargs=2, default=(2, 2), metavar=("PRE", "POST"))
    ap.add_argument("--tol", type=float, default=1e-8)
    ap.add_argument("--jacobi-maxiter", type=int, default=20000)
    ap.add_argument("--amg-maxiter", type=int, default=500)
    ap.add_argument("--no-jacobi", action="store_true"); ap.add_argument("--block-jacobi", action="store_true", help="also run 4x4 block-Jacobi-CG")
    ap.add_argument("--rough-rhs", action="store_true", help="also solve the cpu_reference.py rough RHS (rng(1) normal, ones-mean removed) with Jacobi-CG")
    ap.add_argument("--pyamg-cpu-solve", action="store_true", help="also run pyamg's own CPU CG-accelerated solve for the primary variant (small instances)")
    ap.add_argument("--cpu-only", action="store_true"); ap.add_argument("--no-flexible", action="store_true")
    ap.add_argument("--out", help="result JSON"); ap.add_argument("--tag", default="")
    args = ap.parse_args()

    T0 = time.time()
    def log(*a):
        print(f"[{time.time()-T0:8.1f}s]", *a, flush=True)

    res = dict(tag=args.tag, timestamp=time.strftime("%Y-%m-%dT%H:%M:%S"), host=platform.node(), slurm_job=os.environ.get("SLURM_JOB_ID"),
               python=sys.version.split()[0], numpy=np.__version__, scipy=__import__("scipy").__version__)
    try:
        import pyamg; res["pyamg"] = pyamg.__version__
    except ImportError: res["pyamg"] = None
    gpu = not args.cpu_only
    be = Backend(gpu)
    if gpu:
        cp = be.cp; d = cp.cuda.Device(0); props = cp.cuda.runtime.getDeviceProperties(0)
        res["cupy"] = cp.__version__; res["device"] = props["name"].decode() if isinstance(props["name"], bytes) else str(props["name"])
        res["device_memory_total_gb"] = d.mem_info[1] / 1e9; res["cuda_runtime"] = cp.cuda.runtime.runtimeGetVersion(); res["cuda_driver"] = cp.cuda.runtime.driverGetVersion()
        res["compute_capability"] = d.compute_capability
        log(f"device {res['device']} ({res['device_memory_total_gb']:.1f} GB), cupy {cp.__version__}")
    else:
        res["device"] = "cpu-only (numpy/scipy)"
    mon = DevMemMonitor(be); mon.start()

    # ---- instance ----
    if args.instance and os.path.exists(args.instance):
        t0 = time.time(); d = np.load(args.instance)
        inst = {k: d[k] for k in ("indptr", "indices", "data", "b", "x_true")}; inst.update(nx=int(d["nx"]), ny=int(d["ny"]), nz=int(d["nz"]), h=float(d["h"]))
        log(f"loaded {args.instance} in {time.time()-t0:.1f} s"); res["instance"] = args.instance
    elif args.build:
        nx, ny, nz = args.build
        if args.tmpdir: os.makedirs(args.tmpdir, exist_ok=True)
        tb = time.time(); inst = build_instance(nx, ny, nz, shift=args.shift, workers=args.workers, tmpdir=args.tmpdir, log=log)
        res["build_seconds"] = time.time() - tb; res["build_maxrss_gb"] = maxrss_gb()
        res["instance"] = args.save or f"built in memory {nx}x{ny}x{nz}"
        if args.save:
            ts = time.time(); np.savez(args.save, **inst); log(f"saved {args.save} ({os.path.getsize(args.save)/1e9:.1f} GB) in {time.time()-ts:.0f} s")
            res["instance"] = args.save
    else:
        sys.exit("need --instance or --build")
    indptr, indices, data = inst["indptr"].astype(np.int64, copy=False), inst["indices"].astype(np.int32, copy=False), inst["data"]
    nu = len(indptr) - 1; nnz = int(indptr[-1]); grid = (inst["nx"], inst["ny"], inst["nz"]); h = inst["h"]
    res.update(grid=list(grid), h=h, nu=nu, nnz=nnz, nnz_per_row=nnz / nu, shift=args.shift)
    log(f"instance {grid}: nu {nu:,} nnz {nnz:,} ({nnz/nu:.1f}/row), CSR f64/i32 {(nnz*12+8*(nu+1))/1e9:.2f} GB")
    cpu_ref = os.path.join(HERE, "results", f"p1_{grid[0]}x{grid[1]}x{grid[2]}.cpu.json")
    if os.path.exists(cpu_ref):
        try:
            cr = json.load(open(cpu_ref))
            res["cpu_reference"] = dict(file=os.path.relpath(cpu_ref, HERE), nnz=cr.get("nnz"), condition=cr.get("condition"), lambda_min=cr.get("lambda_min"), lambda_max=cr.get("lambda_max"),
                                        cg_jacobi_iters_rough_rhs=(cr.get("cg_jacobi") or {}).get("iters"), cg_jacobi_seconds_cpu=(cr.get("cg_jacobi") or {}).get("seconds"),
                                        cg_block_jacobi_iters_rough_rhs=(cr.get("cg_block_jacobi_4x4") or {}).get("iters"), smooth_rhs_rel_err_vs_x_true=cr.get("smooth_rhs_rel_err_vs_x_true"),
                                        nnz_matches_instance=(cr.get("nnz") == nnz))
            log(f"cpu reference {cpu_ref}: Jacobi-CG {res['cpu_reference']['cg_jacobi_iters_rough_rhs']} it (rough RHS), cond {cr.get('condition'):.3g}, nnz match {res['cpu_reference']['nnz_matches_instance']}")
        except Exception as ex:
            log(f"cpu reference unreadable: {ex!r}")

    if args.verify_build:
        t0 = time.time(); ref = build_instance(*grid, shift=args.shift, workers=min(args.workers, 4), log=log)
        same_struct = bool(np.array_equal(ref["indptr"], indptr) and np.array_equal(ref["indices"], indices))
        dd = float(np.max(np.abs(ref["data"] - data)) / np.max(np.abs(data))) if same_struct else None
        db = float(np.max(np.abs(ref["b"] - inst["b"])) / np.max(np.abs(inst["b"])))
        res["verify_build"] = dict(same_structure=same_struct, max_rel_data_diff=dd, max_rel_b_diff=db, seconds=time.time() - t0)
        log(f"[verify-build] same structure {same_struct}, max rel data diff {dd}, b diff {db}")
        del ref

    # ---- fine operator on the device ----
    op = FineOperator(be, indptr, indices, data, storage=args.storage, matvec=args.matvec, log=log)
    res.update(storage=args.storage, matvec=("kernel" if args.storage == "sym" else args.matvec), nnz_stored=op.nnz_stored, fine_operator_gb=op.bytes / 1e9)
    res["matvec_seconds"] = op.time_matvec(); res["matvec_effective_GBps"] = op.bytes / res["matvec_seconds"] / 1e9
    log(f"matvec {res['matvec_seconds']*1e3:.2f} ms ({res['matvec_effective_GBps']:.0f} GB/s effective on stored bytes); device used {be.dev_used_gb():.1f} GB")
    if args.validate_sym:
        other = FineOperator(be, indptr, indices, data, storage=("full" if args.storage == "sym" else "sym"), matvec="kernel", log=log)
        rng = np.random.default_rng(7); errs = []
        for trial in range(3):
            x = be.asarray(rng.standard_normal(nu)); y1 = op.matvec(x); y2 = other.matvec(x)
            errs.append(float(be.xp.linalg.norm(y1 - y2) / be.xp.linalg.norm(y1)))
        res["validate_sym"] = dict(max_rel_diff=max(errs), trials=errs, other_storage=other.storage, other_matvec_seconds=other.time_matvec())
        log(f"[validate-sym] max rel diff full vs sym matvec: {max(errs):.2e} (other matvec {res['validate_sym']['other_matvec_seconds']*1e3:.2f} ms)")
        other.free(); del other

    b = be.asarray(inst["b"]); xt = be.asarray(inst["x_true"]); xp = be.xp
    xt_norm = float(xp.linalg.norm(xt))
    A_mv = lambda v, out: op.matvec(v, out)

    # ---- Jacobi-CG baseline ----
    if not args.no_jacobi:
        Dinv = 1.0 / op.diag
        log(f"Jacobi-CG on the stored b (tol {args.tol}, cap {args.jacobi_maxiter})")
        x, info = pcg(be, A_mv, b, lambda r: Dinv * r, tol=args.tol, maxiter=args.jacobi_maxiter, flexible=False, log=log, name="jacobi-cg", log_every=500)
        info["rel_err_vs_x_true"] = float(xp.linalg.norm(x - xt)) / xt_norm
        res["jacobi_cg"] = info; log(f"Jacobi-CG: {info['iters']} it, {info['seconds']:.2f} s, relres {info['final_relres']:.2e}, err vs x_true {info['rel_err_vs_x_true']:.2e}")
        if args.rough_rhs:
            rng = np.random.default_rng(1); br = rng.standard_normal(nu); br -= br.mean()
            brd = be.asarray(br)
            _, info = pcg(be, A_mv, brd, lambda r: Dinv * r, tol=args.tol, maxiter=args.jacobi_maxiter, flexible=False, log=log, name="jacobi-cg-rough", log_every=1000)
            res["jacobi_cg_rough_rhs"] = info; log(f"Jacobi-CG rough RHS (cpu_reference.py RHS): {info['iters']} it, {info['seconds']:.2f} s")
            del brd
        if args.block_jacobi:
            # 4x4 block diagonal inverse (host), applied with einsum on the device
            nc = nu // 4; blk = np.zeros((nc, 4, 4))
            for r0 in range(0, nu, 2_000_000):
                r1 = min(r0 + 2_000_000, nu); p0, p1 = indptr[r0], indptr[r1]
                rows = np.repeat(np.arange(r0, r1), np.diff(indptr[r0:r1 + 1])); cols = indices[p0:p1].astype(np.int64)
                m = (cols // 4) == (rows // 4)
                blk[rows[m] // 4, rows[m] % 4, cols[m] % 4] = data[p0:p1][m]
            Binv = be.asarray(np.linalg.inv(blk)); del blk
            M = lambda r: xp.einsum("cij,cj->ci", Binv, r.reshape(-1, 4)).ravel()
            _, info = pcg(be, A_mv, b, M, tol=args.tol, maxiter=args.jacobi_maxiter, flexible=False, log=log, name="bjacobi-cg", log_every=500)
            res["block_jacobi_cg"] = info; log(f"block-Jacobi(4x4)-CG: {info['iters']} it, {info['seconds']:.2f} s")
            del Binv
        del Dinv, x
    res["device_used_gb_after_baseline"] = be.dev_used_gb()

    # ---- AMG-CG ----
    variants = [v for v in args.amg.split(",") if v and v != "none"]
    if variants:
        t0 = time.time()
        A_host = sp.csr_matrix((data, indices, indptr), shape=(nu, nu))        # scipy may upcast indices to int64 when nnz > 2^31
        A_host.has_canonical_format = True
        log(f"host scipy CSR for pyamg: index dtype {A_host.indices.dtype}, {time.time()-t0:.1f} s, maxrss {maxrss_gb():.1f} GB")
        res["amg_variants"] = {}
        for vi, variant in enumerate(variants):
            try:
                ml, info = build_pyamg(A_host, variant, grid, h, max_coarse=args.max_coarse, log=log)
                info["host_maxrss_gb_after_setup"] = maxrss_gb()
                if args.pyamg_cpu_solve and vi == 0:
                    tc = time.time(); resid = []
                    ml.solve(inst["b"], tol=args.tol, maxiter=args.amg_maxiter, accel="cg", residuals=resid)
                    info["pyamg_cpu_cg_iters"] = len(resid) - 1; info["pyamg_cpu_cg_seconds"] = time.time() - tc
                    log(f"[amg] pyamg CPU cg solve (its own smoothers): {len(resid)-1} iterations, {time.time()-tc:.1f} s")
                vc = DeviceVCycle(be, op, ml, omega_factor=args.omega_factor, nu1=args.sweeps[0], nu2=args.sweeps[1], log=log)
                del ml
                info.update(device_setup_s=vc.setup_s, hierarchy_device_gb=vc.bytes / 1e9, omega=[float(w) for w in vc.omega], sweeps=list(args.sweeps),
                            omega_factor=args.omega_factor, device_used_gb_after_setup=be.dev_used_gb())
                x, cg = pcg(be, A_mv, b, vc, tol=args.tol, maxiter=args.amg_maxiter, flexible=not args.no_flexible, log=log, name=f"amg-cg/{variant}", log_every=10)
                cg["rel_err_vs_x_true"] = float(xp.linalg.norm(x - xt)) / xt_norm
                # one V-cycle timing
                be.sync(); tv = time.time(); vc(b); be.sync(); cg["vcycle_seconds"] = time.time() - tv
                info["cg"] = cg
                log(f"AMG-CG [{variant}]: {cg['iters']} it, {cg['seconds']:.2f} s ({cg['seconds_per_iter']*1e3:.0f} ms/it, V-cycle {cg['vcycle_seconds']*1e3:.0f} ms), relres {cg['final_relres']:.2e}, err vs x_true {cg['rel_err_vs_x_true']:.2e}")
                del vc, x
                if be.gpu: be.pool.free_all_blocks()
            except Exception as ex:                       # record, do not hide
                import traceback; info = dict(variant=variant, error=repr(ex), traceback=traceback.format_exc()[-2000:])
                log(f"AMG variant {variant} FAILED: {ex!r}")
            # release device memory of this variant only AFTER the except block: the live exception's traceback
            # keeps the half-built hierarchy referenced (job 147436's geo_lin attempt was confounded by this)
            import gc; ml = vc = x = None; gc.collect()
            if be.gpu:
                be.pool.free_all_blocks(); info["device_used_gb_after_release"] = be.dev_used_gb()
            res["amg_variants"][variant] = info
        prim = res["amg_variants"].get(variants[0], {})
        res["amg_cg"] = dict(variant=variants[0], **{k: prim.get(k) for k in ("levels", "operator_complexity", "grid_complexity", "host_setup_s", "device_setup_s", "hierarchy_device_gb")}, **(prim.get("cg") or {}))
        del A_host

    res["peak_device_memory_gb"] = mon.stop(); res["device_mempool_total_gb"] = be.pool_gb(); res["host_maxrss_gb"] = maxrss_gb(); res["total_seconds"] = time.time() - T0
    log(f"peak device memory {res['peak_device_memory_gb']:.1f} GB (mempool {res['device_mempool_total_gb']:.1f} GB), host maxrss {res['host_maxrss_gb']:.1f} GB, total {res['total_seconds']:.0f} s")
    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w") as f: json.dump(res, f, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
        log(f"wrote {args.out}")


if __name__ == "__main__":
    main()
