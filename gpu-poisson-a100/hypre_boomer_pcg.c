/* hypre_boomer_pcg.c -- BoomerAMG-PCG (and Jacobi-PCG baseline) on a synthetic P1-disc pressure-Poisson
 * instance of gpu-poisson-a100 (README.md §5 route 2), for a CUDA build of hypre 3.1.0 (mixed-int:
 * HYPRE_BigInt 64-bit global ids, HYPRE_Int 32-bit local counts).
 *
 * Input: the raw files written by npz_to_bin.py (PREFIX.meta.txt, .indptr.i64, .indices.i32, .data.f64,
 * .b.f64, .xtrue.f64, .rough.f64, .cand.f64). Rows are split across the MPI ranks on x-plane boundaries
 * (multiples of 4*ny*nz rows, so the 4 unknowns of a cell never straddle ranks); every rank needs
 * < 2^31 local non-zeros (HYPRE_Int row pointers), hence the half/full-size instances need 3/6 ranks on
 * the single GPU (ranks time-share the A100).
 *
 * Assembly paths: --direct (default) fills a hypre_ParCSRMatrix (diag/offd split done here on the host,
 * 12 B/nnz) and migrates it to the device: the only path whose device footprint is the CSR itself.
 * --ij uses HYPRE_IJMatrixSetValues2 with device arrays (the brief's request): hypre keeps a COO stash
 * (24 B/nnz) plus sort temporaries until Assemble, so it is used for the correctness check on the small
 * instance and documented as not fitting the large ones.
 *
 * Output: one JSON per run (--out) with the fields of the CuPy JSONs (library, settings, levels, operator
 * complexity, setup s, solve s, iterations, final relres, true relres, rel err vs x_true, device memory
 * samples), plus hypre's own printouts in the job log. A preliminary JSON is written as soon as the matrix
 * is on the device so an OOM in the AMG setup still leaves a record.
 */
#define _FILE_OFFSET_BITS 64
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdarg.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <mpi.h>
#include "HYPRE.h"
#include "HYPRE_utilities.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_krylov.h"
#include "_hypre_utilities.h"
#include "_hypre_parcsr_mv.h"
#include "_hypre_parcsr_ls.h"
#include "_hypre_krylov.h"
#if defined(HYPRE_USING_CUDA)
#include <cuda_runtime.h>
#endif

#define MAXLEV 64
#define NSAMP 32

typedef struct {
   const char *bin, *out, *label;
   int ij, host, jacobi, amg, rough_rhs;
   int num_functions, nodal, nodal_diag, interp_vectors, interp_vec_variant, interp_vec_qmax;
   int coarsen, interp, relax, agg_levels, agg_interp, pmax, max_levels, max_coarse, num_sweeps, keep_transpose, mod_rap2;
   double theta, tol; int maxit, spmv_vendor, spgemm_vendor, amg_print, cheby_order, cheby_variant; double cheby_fraction;
} opts_t;

static int rank, nprocs;
static double T0;
static double now(void) { struct timeval tv; gettimeofday(&tv, NULL); return tv.tv_sec + 1e-6 * tv.tv_usec; }
static void logp(const char *fmt, ...) __attribute__((format(printf, 1, 2)));
static void logp(const char *fmt, ...) {
   if (rank) return;
   va_list ap; va_start(ap, fmt); printf("[%8.1fs] ", now() - T0); vprintf(fmt, ap); printf("\n"); fflush(stdout); va_end(ap);
}
static void die(const char *msg) { fprintf(stderr, "rank %d: FATAL %s\n", rank, msg); fflush(stderr); MPI_Abort(MPI_COMM_WORLD, 1); exit(1); }
static void sync_all(void) {
#if defined(HYPRE_USING_CUDA)
   cudaDeviceSynchronize();
#endif
   MPI_Barrier(MPI_COMM_WORLD);
}
static double maxrss_gb(void) { struct rusage ru; getrusage(RUSAGE_SELF, &ru); return ru.ru_maxrss / 1e6; }

/* device memory samples (cudaMemGetInfo is per device, i.e. the aggregate over all ranks' contexts) */
static struct { const char *label; double used_gb; } samp[NSAMP]; static int nsamp = 0; static double dev_total_gb = 0, dev_peak_gb = 0;
static double devmem_used_gb(void) {
#if defined(HYPRE_USING_CUDA)
   size_t f = 0, t = 0; cudaDeviceSynchronize(); if (cudaMemGetInfo(&f, &t) != cudaSuccess) return -1; dev_total_gb = t / 1e9; return (t - f) / 1e9;
#else
   return 0;
#endif
}
static void sample(const char *label) {
   sync_all(); double u = devmem_used_gb(); double umax; MPI_Allreduce(&u, &umax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   if (umax > dev_peak_gb) dev_peak_gb = umax;
   if (nsamp < NSAMP) { samp[nsamp].label = label; samp[nsamp].used_gb = umax; nsamp++; }
   logp("[devmem] %-28s used %.2f GB of %.1f GB (peak %.2f), host maxrss %.1f GB", label, umax, dev_total_gb, dev_peak_gb, maxrss_gb());
}

/* ---- raw file helpers ---- */
static void read_at(const char *path, void *dst, int64_t offset_bytes, int64_t nbytes) {
   FILE *f = fopen(path, "rb"); if (!f) { fprintf(stderr, "cannot open %s\n", path); die("open"); }
   if (fseeko(f, (off_t) offset_bytes, SEEK_SET)) die("fseeko");
   char *p = (char *) dst; int64_t left = nbytes;
   while (left > 0) { size_t chunk = left > (1 << 30) ? (1 << 30) : (size_t) left; size_t got = fread(p, 1, chunk, f); if (got != chunk) { fprintf(stderr, "short read %s (%zu of %zu)\n", path, got, chunk); die("read"); } p += got; left -= got; }
   fclose(f);
}
static char *pathf(const char *prefix, const char *suffix) { char *s = (char *) malloc(strlen(prefix) + strlen(suffix) + 2); sprintf(s, "%s.%s", prefix, suffix); return s; }
static int cmp_ll(const void *a, const void *b) { long long x = *(const long long *) a, y = *(const long long *) b; return (x > y) - (x < y); }
static HYPRE_Int bsearch_ll(const long long *arr, HYPRE_Int n, long long key) {
   HYPRE_Int lo = 0, hi = n - 1; while (lo <= hi) { HYPRE_Int mid = lo + (hi - lo) / 2; if (arr[mid] == key) return mid; if (arr[mid] < key) lo = mid + 1; else hi = mid - 1; } return -1;
}

/* ---- instance ---- */
typedef struct {
   int nx, ny, nz; double h; long long nu, nnz;
   long long r0, r1; HYPRE_Int nloc; long long nnz_loc;
   int64_t *ip; int32_t *jx; double *ax;          /* local CSR rows (global columns) */
   double *b, *xt, *rough, *cand;                  /* local slices; cand: 4 x nloc */
} inst_t;

static void read_meta(const char *prefix, inst_t *I) {
   char *m = pathf(prefix, "meta.txt"); FILE *f = fopen(m, "r"); if (!f) { fprintf(stderr, "cannot open %s\n", m); die("meta"); }
   char line[512]; I->nx = I->ny = I->nz = 0; I->nu = I->nnz = 0; I->h = 0;
   while (fgets(line, sizeof line, f)) {
      long long v; double dv;
      if (sscanf(line, "nx=%d", &I->nx) == 1) continue;
      if (sscanf(line, "ny=%d", &I->ny) == 1) continue;
      if (sscanf(line, "nz=%d", &I->nz) == 1) continue;
      if (sscanf(line, "h=%lf", &dv) == 1) { I->h = dv; continue; }
      if (sscanf(line, "nu=%lld", &v) == 1) { I->nu = v; continue; }
      if (sscanf(line, "nnz=%lld", &v) == 1) { I->nnz = v; continue; }
   }
   fclose(f); free(m);
   if (!I->nx || !I->nu || !I->nnz) die("bad meta");
}

static void load_local(const char *prefix, inst_t *I) {
   long long plane = 4LL * I->ny * I->nz;
   int p0 = (int) (((long long) I->nx * rank) / nprocs), p1 = (int) (((long long) I->nx * (rank + 1)) / nprocs);
   I->r0 = plane * p0; I->r1 = plane * p1; I->nloc = (HYPRE_Int) (I->r1 - I->r0);
   I->ip = (int64_t *) malloc(sizeof(int64_t) * (I->nloc + 1));
   char *p = pathf(prefix, "indptr.i64"); read_at(p, I->ip, 8 * I->r0, 8LL * (I->nloc + 1)); free(p);
   I->nnz_loc = I->ip[I->nloc] - I->ip[0];
   if (sizeof(HYPRE_Int) == 4 && I->nnz_loc >= 2147483647LL) { fprintf(stderr, "rank %d: local nnz %lld >= 2^31: HYPRE_Int is 32-bit (mixed-int build); use more ranks\n", rank, I->nnz_loc); die("nnz_loc"); }
   I->jx = (int32_t *) malloc(4 * (size_t) I->nnz_loc); I->ax = (double *) malloc(8 * (size_t) I->nnz_loc);
   p = pathf(prefix, "indices.i32"); read_at(p, I->jx, 4 * I->ip[0], 4 * I->nnz_loc); free(p);
   p = pathf(prefix, "data.f64");    read_at(p, I->ax, 8 * I->ip[0], 8 * I->nnz_loc); free(p);
   I->b = (double *) malloc(8 * (size_t) I->nloc); I->xt = (double *) malloc(8 * (size_t) I->nloc); I->rough = (double *) malloc(8 * (size_t) I->nloc); I->cand = (double *) malloc(8 * 4 * (size_t) I->nloc);
   p = pathf(prefix, "b.f64");     read_at(p, I->b, 8 * I->r0, 8LL * I->nloc); free(p);
   p = pathf(prefix, "xtrue.f64"); read_at(p, I->xt, 8 * I->r0, 8LL * I->nloc); free(p);
   p = pathf(prefix, "rough.f64"); read_at(p, I->rough, 8 * I->r0, 8LL * I->nloc); free(p);
   p = pathf(prefix, "cand.f64");  for (int c = 0; c < 4; c++) read_at(p, I->cand + (size_t) c * I->nloc, 8 * (c * I->nu + I->r0), 8LL * I->nloc); free(p);
   for (HYPRE_Int i = I->nloc; i >= 0; i--) I->ip[i] -= I->ip[0];   /* local row pointer */
   long long mx; MPI_Allreduce(&I->nnz_loc, &mx, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
   logp("instance %dx%dx%d: nu %lld nnz %lld (%.1f/row); %d ranks, planes/rank ~%d, max local nnz %lld (%.2f G, %s 2^31; HYPRE_Int %d B)",
        I->nx, I->ny, I->nz, I->nu, I->nnz, (double) I->nnz / I->nu, nprocs, I->nx / nprocs, mx, mx / 1e9, mx < 2147483647LL ? "<" : ">=", (int) sizeof(HYPRE_Int));
}

/* ---- direct ParCSR assembly: diag/offd split on the host, then migrate ---- */
static hypre_ParCSRMatrix *build_direct(inst_t *I, HYPRE_MemoryLocation memloc) {
   HYPRE_Int nloc = I->nloc; long long r0 = I->r0, r1 = I->r1;
   HYPRE_Int nnz_d = 0, nnz_o = 0;
   for (long long k = 0; k < I->nnz_loc; k++) { long long c = I->jx[k]; if (c >= r0 && c < r1) nnz_d++; else nnz_o++; }
   long long *ocols = (long long *) malloc(sizeof(long long) * (nnz_o > 0 ? nnz_o : 1)); HYPRE_Int m = 0;
   for (long long k = 0; k < I->nnz_loc; k++) { long long c = I->jx[k]; if (!(c >= r0 && c < r1)) ocols[m++] = c; }
   qsort(ocols, (size_t) nnz_o, sizeof(long long), cmp_ll);
   HYPRE_Int ncol_o = 0; for (HYPRE_Int k = 0; k < nnz_o; k++) if (k == 0 || ocols[k] != ocols[k - 1]) ocols[ncol_o++] = ocols[k];
   HYPRE_BigInt row_starts[2] = { (HYPRE_BigInt) r0, (HYPRE_BigInt) r1 };
   hypre_ParCSRMatrix *A = hypre_ParCSRMatrixCreate(MPI_COMM_WORLD, (HYPRE_BigInt) I->nu, (HYPRE_BigInt) I->nu, row_starts, row_starts, ncol_o, nnz_d, nnz_o);
   hypre_ParCSRMatrixInitialize_v2(A, HYPRE_MEMORY_HOST);
   hypre_CSRMatrix *D = hypre_ParCSRMatrixDiag(A), *O = hypre_ParCSRMatrixOffd(A);
   HYPRE_Int *Di = hypre_CSRMatrixI(D), *Dj = hypre_CSRMatrixJ(D), *Oi = hypre_CSRMatrixI(O), *Oj = hypre_CSRMatrixJ(O);
   HYPRE_Real *Da = hypre_CSRMatrixData(D), *Oa = hypre_CSRMatrixData(O); HYPRE_BigInt *cmap = hypre_ParCSRMatrixColMapOffd(A);
   for (HYPRE_Int k = 0; k < ncol_o; k++) cmap[k] = (HYPRE_BigInt) ocols[k];
   HYPRE_Int kd = 0, ko = 0;
   /* hypre convention (relied on by DiagScale, the relaxations, strength-of-connection, ...): the DIAGONAL entry is
      the FIRST entry of each row of the diag block -- the IJ interface guarantees it, a direct fill must do it. */
   HYPRE_Int missing_diag = 0;
   for (HYPRE_Int i = 0; i < nloc; i++) {
      Di[i] = kd; Oi[i] = ko; long long gi = r0 + i; int found = 0;
      for (int64_t k = I->ip[i]; k < I->ip[i + 1]; k++) if (I->jx[k] == gi) { Dj[kd] = i; Da[kd++] = I->ax[k]; found = 1; break; }
      if (!found) missing_diag++;
      for (int64_t k = I->ip[i]; k < I->ip[i + 1]; k++) {
         long long c = I->jx[k]; if (c == gi) continue;
         if (c >= r0 && c < r1) { Dj[kd] = (HYPRE_Int) (c - r0); Da[kd++] = I->ax[k]; }
         else { HYPRE_Int j = bsearch_ll(ocols, ncol_o, c); if (j < 0) die("offd map"); Oj[ko] = j; Oa[ko++] = I->ax[k]; }
      }
   }
   if (missing_diag) { fprintf(stderr, "rank %d: %d rows without a stored diagonal entry\n", rank, (int) missing_diag); die("diagonal"); }
   Di[nloc] = kd; Oi[nloc] = ko; free(ocols);
   free(I->jx); free(I->ax); free(I->ip); I->jx = NULL; I->ax = NULL; I->ip = NULL;
   hypre_ParCSRMatrixSetNumNonzeros(A); hypre_ParCSRMatrixSetDNumNonzeros(A);
   long long mo, mo_l = ncol_o; MPI_Allreduce(&mo_l, &mo, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
   logp("[direct] ParCSR on host (diagonal first in every row): global nnz %.0f, max offd cols/rank %d", hypre_ParCSRMatrixDNumNonzeros(A), (int) mo);
   if (memloc == HYPRE_MEMORY_DEVICE) { double t = now(); hypre_ParCSRMatrixMigrate(A, HYPRE_MEMORY_DEVICE); sync_all(); logp("[direct] migrated to device in %.1f s", now() - t); }
   return A;
}

/* ---- IJ assembly from device arrays (HYPRE_IJMatrixSetValues2), chunked by rows ---- */
static hypre_ParCSRMatrix *build_ij(inst_t *I, HYPRE_MemoryLocation memloc) {
   HYPRE_IJMatrix ij; HYPRE_Int nloc = I->nloc;
   HYPRE_IJMatrixCreate(MPI_COMM_WORLD, (HYPRE_BigInt) I->r0, (HYPRE_BigInt) (I->r1 - 1), (HYPRE_BigInt) I->r0, (HYPRE_BigInt) (I->r1 - 1), &ij);
   HYPRE_IJMatrixSetObjectType(ij, HYPRE_PARCSR);
   HYPRE_IJMatrixInitialize_v2(ij, memloc);
   const int64_t chunk_nnz = 1LL << 28; HYPRE_Int i0 = 0; int nchunks = 0;
   while (i0 < nloc) {
      HYPRE_Int i1 = i0; while (i1 < nloc && I->ip[i1 + 1] - I->ip[i0] <= chunk_nnz) i1++; if (i1 == i0) i1 = i0 + 1;
      HYPRE_Int nr = i1 - i0; int64_t k0 = I->ip[i0], k1 = I->ip[i1]; HYPRE_Int nk = (HYPRE_Int) (k1 - k0);
      HYPRE_Int *ncols_h = hypre_TAlloc(HYPRE_Int, nr, HYPRE_MEMORY_HOST), *ridx_h = hypre_TAlloc(HYPRE_Int, nr, HYPRE_MEMORY_HOST);
      HYPRE_BigInt *rows_h = hypre_TAlloc(HYPRE_BigInt, nr, HYPRE_MEMORY_HOST), *cols_h = hypre_TAlloc(HYPRE_BigInt, nk, HYPRE_MEMORY_HOST);
      for (HYPRE_Int i = 0; i < nr; i++) { ncols_h[i] = (HYPRE_Int) (I->ip[i0 + i + 1] - I->ip[i0 + i]); ridx_h[i] = (HYPRE_Int) (I->ip[i0 + i] - k0); rows_h[i] = (HYPRE_BigInt) (I->r0 + i0 + i); }
      for (HYPRE_Int k = 0; k < nk; k++) cols_h[k] = (HYPRE_BigInt) I->jx[k0 + k];
      HYPRE_Int *ncols = hypre_TAlloc(HYPRE_Int, nr, memloc), *ridx = hypre_TAlloc(HYPRE_Int, nr, memloc);
      HYPRE_BigInt *rows = hypre_TAlloc(HYPRE_BigInt, nr, memloc), *cols = hypre_TAlloc(HYPRE_BigInt, nk, memloc); HYPRE_Real *vals = hypre_TAlloc(HYPRE_Real, nk, memloc);
      hypre_TMemcpy(ncols, ncols_h, HYPRE_Int, nr, memloc, HYPRE_MEMORY_HOST); hypre_TMemcpy(ridx, ridx_h, HYPRE_Int, nr, memloc, HYPRE_MEMORY_HOST);
      hypre_TMemcpy(rows, rows_h, HYPRE_BigInt, nr, memloc, HYPRE_MEMORY_HOST); hypre_TMemcpy(cols, cols_h, HYPRE_BigInt, nk, memloc, HYPRE_MEMORY_HOST);
      hypre_TMemcpy(vals, I->ax + k0, HYPRE_Real, nk, memloc, HYPRE_MEMORY_HOST);
      HYPRE_IJMatrixSetValues2(ij, nr, ncols, rows, ridx, cols, vals);
      hypre_TFree(ncols, memloc); hypre_TFree(ridx, memloc); hypre_TFree(rows, memloc); hypre_TFree(cols, memloc); hypre_TFree(vals, memloc);
      hypre_TFree(ncols_h, HYPRE_MEMORY_HOST); hypre_TFree(ridx_h, HYPRE_MEMORY_HOST); hypre_TFree(rows_h, HYPRE_MEMORY_HOST); hypre_TFree(cols_h, HYPRE_MEMORY_HOST);
      i0 = i1; nchunks++;
   }
   sample("after IJ SetValues (stash)");
   double t = now(); HYPRE_IJMatrixAssemble(ij); sync_all(); logp("[ij] %d SetValues2 chunks, Assemble %.1f s", nchunks, now() - t);
   free(I->jx); free(I->ax); free(I->ip); I->jx = NULL; I->ax = NULL; I->ip = NULL;
   hypre_ParCSRMatrix *A; HYPRE_IJMatrixGetObject(ij, (void **) &A);   /* the IJ wrapper is kept alive (owns the ParCSR) */
   hypre_ParCSRMatrixSetNumNonzeros(A); hypre_ParCSRMatrixSetDNumNonzeros(A);
   logp("[ij] ParCSR global nnz %.0f", hypre_ParCSRMatrixDNumNonzeros(A));
   return A;
}

static hypre_ParVector *make_vec(inst_t *I, const double *host_or_null, HYPRE_MemoryLocation memloc) {
   HYPRE_BigInt part[2] = { (HYPRE_BigInt) I->r0, (HYPRE_BigInt) I->r1 };
   hypre_ParVector *v = hypre_ParVectorCreate(MPI_COMM_WORLD, (HYPRE_BigInt) I->nu, part);
   hypre_ParVectorInitialize_v2(v, memloc);
   if (host_or_null) hypre_TMemcpy(hypre_VectorData(hypre_ParVectorLocalVector(v)), host_or_null, HYPRE_Real, I->nloc, memloc, HYPRE_MEMORY_HOST);
   else hypre_ParVectorSetConstantValues(v, 0.0);
   return v;
}
static void vec_to_host(hypre_ParVector *v, double *dst, HYPRE_Int n, HYPRE_MemoryLocation memloc) { hypre_TMemcpy(dst, hypre_VectorData(hypre_ParVectorLocalVector(v)), HYPRE_Real, n, HYPRE_MEMORY_HOST, memloc); }

/* ---- results ---- */
typedef struct {
   int ran, converged, iters, levels; double setup_s, solve_s, final_relres, true_relres, rel_err, oc, gc, hier_nnz_A, hier_nnz_P, hier_nnz_R, hier_gb_est;
   double sizes[MAXLEV], nnzs[MAXLEV], hist[64]; int nhist; double dev_after_setup_gb, dev_after_solve_gb;
} res_t;

static double norm2_global(const double *a, HYPRE_Int n) { double s = 0, g; for (HYPRE_Int i = 0; i < n; i++) s += a[i] * a[i]; MPI_Allreduce(&s, &g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); return sqrt(g); }

static void run_pcg(const opts_t *o, inst_t *I, hypre_ParCSRMatrix *A, hypre_ParVector *b, hypre_ParVector *x, const double *xt_host, int use_amg, HYPRE_MemoryLocation memloc, res_t *R, const char *tag) {
   HYPRE_Solver pcg, amg = NULL; memset(R, 0, sizeof *R); R->ran = 1;
   HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &pcg);
   HYPRE_PCGSetMaxIter(pcg, o->maxit); HYPRE_PCGSetTol(pcg, o->tol); HYPRE_PCGSetTwoNorm(pcg, 1); HYPRE_PCGSetRelChange(pcg, 0);
   HYPRE_PCGSetPrintLevel(pcg, 2); HYPRE_PCGSetLogging(pcg, 1);
   if (use_amg) {
      HYPRE_BoomerAMGCreate(&amg);
      HYPRE_BoomerAMGSetPrintLevel(amg, o->amg_print); HYPRE_BoomerAMGSetMaxIter(amg, 1); HYPRE_BoomerAMGSetTol(amg, 0.0);
      HYPRE_BoomerAMGSetCoarsenType(amg, o->coarsen); HYPRE_BoomerAMGSetInterpType(amg, o->interp); HYPRE_BoomerAMGSetRelaxType(amg, o->relax);
      if (o->relax == 16) { HYPRE_BoomerAMGSetChebyOrder(amg, o->cheby_order); HYPRE_BoomerAMGSetChebyFraction(amg, o->cheby_fraction); HYPRE_BoomerAMGSetChebyVariant(amg, o->cheby_variant); HYPRE_BoomerAMGSetChebyEigEst(amg, 10); }
      HYPRE_BoomerAMGSetNumSweeps(amg, o->num_sweeps); HYPRE_BoomerAMGSetAggNumLevels(amg, o->agg_levels); if (o->agg_interp > 0) HYPRE_BoomerAMGSetAggInterpType(amg, o->agg_interp);
      HYPRE_BoomerAMGSetPMaxElmts(amg, o->pmax); HYPRE_BoomerAMGSetStrongThreshold(amg, o->theta); HYPRE_BoomerAMGSetMaxLevels(amg, o->max_levels); HYPRE_BoomerAMGSetMaxCoarseSize(amg, o->max_coarse);
      HYPRE_BoomerAMGSetKeepTranspose(amg, o->keep_transpose); HYPRE_BoomerAMGSetModuleRAP2(amg, o->mod_rap2); HYPRE_BoomerAMGSetRAP2(amg, 0);
      if (o->num_functions > 1) HYPRE_BoomerAMGSetNumFunctions(amg, o->num_functions);
      if (o->nodal) { HYPRE_BoomerAMGSetNodal(amg, o->nodal); HYPRE_BoomerAMGSetNodalDiag(amg, o->nodal_diag); }
      if (o->interp_vectors > 0) {
         /* candidates [e, x, y, z] in the P1-disc basis; 3 -> [x, y, z] (the per-function constants are interpolated anyway), 4 -> all */
         int nv = o->interp_vectors; int first = (nv == 3) ? 1 : 0; if (nv > 4) nv = 4;
         HYPRE_ParVector *vecs = (HYPRE_ParVector *) malloc(sizeof(HYPRE_ParVector) * nv);
         for (int c = 0; c < nv; c++) vecs[c] = (HYPRE_ParVector) make_vec(I, I->cand + (size_t) (first + c) * I->nloc, memloc);
         HYPRE_BoomerAMGSetInterpVectors(amg, nv, vecs); HYPRE_BoomerAMGSetInterpVecVariant(amg, o->interp_vec_variant);
         if (o->interp_vec_qmax > 0) HYPRE_BoomerAMGSetInterpVecQMax(amg, o->interp_vec_qmax);
         logp("[amg] %d interpolation vectors (%s), GM/LN variant %d, qmax %d", nv, first ? "x,y,z" : "e,x,y,z", o->interp_vec_variant, o->interp_vec_qmax);
      }
      HYPRE_PCGSetPrecond(pcg, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, amg);
   } else {
      HYPRE_PCGSetPrecond(pcg, (HYPRE_PtrToSolverFcn) HYPRE_ParCSRDiagScale, (HYPRE_PtrToSolverFcn) HYPRE_ParCSRDiagScaleSetup, NULL);
   }
   hypre_ParVectorSetConstantValues(x, 0.0);
   sync_all(); double t = now();
   HYPRE_ParCSRPCGSetup(pcg, (HYPRE_ParCSRMatrix) A, (HYPRE_ParVector) b, (HYPRE_ParVector) x);
   sync_all(); R->setup_s = now() - t;
   logp("[%s] setup %.3f s", tag, R->setup_s);
   if (use_amg) {
      hypre_ParAMGData *ad = (hypre_ParAMGData *) amg; int nl = hypre_ParAMGDataNumLevels(ad); R->levels = nl;
      hypre_ParCSRMatrix **Aa = hypre_ParAMGDataAArray(ad), **Pa = hypre_ParAMGDataPArray(ad), **Ra = hypre_ParAMGDataRArray(ad);
      double sn = 0, sz = 0;
      for (int l = 0; l < nl && l < MAXLEV; l++) { hypre_ParCSRMatrixSetDNumNonzeros(Aa[l]); R->nnzs[l] = hypre_ParCSRMatrixDNumNonzeros(Aa[l]); R->sizes[l] = (double) hypre_ParCSRMatrixGlobalNumRows(Aa[l]); sn += R->nnzs[l]; sz += R->sizes[l]; }
      for (int l = 0; l < nl - 1; l++) { hypre_ParCSRMatrixSetDNumNonzeros(Pa[l]); R->hier_nnz_P += hypre_ParCSRMatrixDNumNonzeros(Pa[l]); if (Ra && Ra[l] && Ra[l] != Pa[l]) { hypre_ParCSRMatrixSetDNumNonzeros(Ra[l]); R->hier_nnz_R += hypre_ParCSRMatrixDNumNonzeros(Ra[l]); } }
      R->hier_nnz_A = sn; R->oc = sn / R->nnzs[0]; R->gc = sz / R->sizes[0];
      R->hier_gb_est = 12e-9 * (sn - R->nnzs[0] + R->hier_nnz_P + R->hier_nnz_R);
      if (!rank) { printf("[%s] %d levels, OC %.3f, GC %.3f, nnz(P) %.3g, nnz(R) %.3g, hierarchy est. %.1f GB; sizes:", tag, nl, R->oc, R->gc, R->hier_nnz_P, R->hier_nnz_R, R->hier_gb_est); for (int l = 0; l < nl; l++) printf(" %.0f", R->sizes[l]); printf("\n"); fflush(stdout); }
   }
   sample(use_amg ? "after AMG setup" : "after DiagScale setup"); R->dev_after_setup_gb = samp[nsamp - 1].used_gb;
   sync_all(); t = now();
   HYPRE_ParCSRPCGSolve(pcg, (HYPRE_ParCSRMatrix) A, (HYPRE_ParVector) b, (HYPRE_ParVector) x);
   sync_all(); R->solve_s = now() - t;
   HYPRE_Int it; HYPRE_Real rr; HYPRE_Int conv; HYPRE_PCGGetNumIterations(pcg, &it); HYPRE_PCGGetFinalRelativeResidualNorm(pcg, &rr); HYPRE_PCGGetConverged(pcg, &conv);
   R->iters = it; R->final_relres = rr; R->converged = conv;
   hypre_PCGData *pd = (hypre_PCGData *) pcg;
   if (pd->norms && pd->norms[0] > 0) { R->nhist = 0; for (int i = 0; i <= it && R->nhist < 64; i += 10) R->hist[R->nhist++] = pd->norms[i] / pd->norms[0]; }
   /* true residual and error vs x_true */
   hypre_ParVector *r = make_vec(I, NULL, memloc); hypre_ParVectorCopy(b, r); hypre_ParCSRMatrixMatvec(-1.0, A, x, 1.0, r);
   R->true_relres = sqrt(hypre_ParVectorInnerProd(r, r)) / sqrt(hypre_ParVectorInnerProd(b, b)); hypre_ParVectorDestroy(r);
   if (xt_host) { double *xh = (double *) malloc(8 * (size_t) I->nloc); vec_to_host(x, xh, I->nloc, memloc); double *d = (double *) malloc(8 * (size_t) I->nloc); for (HYPRE_Int i = 0; i < I->nloc; i++) d[i] = xh[i] - xt_host[i]; R->rel_err = norm2_global(d, I->nloc) / norm2_global(xt_host, I->nloc); free(xh); free(d); }
   logp("[%s] %d it (%s), solve %.3f s (%.1f ms/it), final relres %.3e, true relres %.3e, rel err vs x_true %.3e", tag, (int) it, conv ? "converged" : "NOT converged", R->solve_s, 1e3 * R->solve_s / (it > 0 ? it : 1), (double) rr, R->true_relres, R->rel_err);
   sample("after solve"); R->dev_after_solve_gb = samp[nsamp - 1].used_gb;
   if (amg) HYPRE_BoomerAMGDestroy(amg);
   HYPRE_ParCSRPCGDestroy(pcg);
   sample(use_amg ? "after AMG destroy" : "after PCG destroy");
}

static void json_res(FILE *f, const char *key, const res_t *R, int amg) {
   if (!R->ran) { fprintf(f, "  \"%s\": null", key); return; }
   fprintf(f, "  \"%s\": {\"iters\": %d, \"converged\": %s, \"setup_seconds\": %.6f, \"seconds\": %.6f, \"seconds_per_iter\": %.6g, \"final_relres\": %.6e, \"true_relres\": %.6e, \"rel_err_vs_x_true\": %.6e, \"device_used_gb_after_setup\": %.3f, \"device_used_gb_after_solve\": %.3f",
           key, R->iters, R->converged ? "true" : "false", R->setup_s, R->solve_s, R->solve_s / (R->iters > 0 ? R->iters : 1), R->final_relres, R->true_relres, R->rel_err, R->dev_after_setup_gb, R->dev_after_solve_gb);
   if (amg) {
      fprintf(f, ", \"levels\": %d, \"operator_complexity\": %.6f, \"grid_complexity\": %.6f, \"nnz_P_total\": %.0f, \"nnz_R_total\": %.0f, \"hierarchy_gb_est_12B_per_nnz\": %.3f, \"sizes\": [", R->levels, R->oc, R->gc, R->hier_nnz_P, R->hier_nnz_R, R->hier_gb_est);
      for (int l = 0; l < R->levels; l++) fprintf(f, "%s%.0f", l ? ", " : "", R->sizes[l]);
      fprintf(f, "], \"nnz_per_level\": [");
      for (int l = 0; l < R->levels; l++) fprintf(f, "%s%.0f", l ? ", " : "", R->nnzs[l]);
      fprintf(f, "]");
   }
   fprintf(f, ", \"hist_every_10\": ["); for (int i = 0; i < R->nhist; i++) fprintf(f, "%s%.6e", i ? ", " : "", R->hist[i]); fprintf(f, "]}");
}

static void write_json(const opts_t *o, const inst_t *I, const char *status, const res_t *Rj, const res_t *Rr, const res_t *Ra, HYPRE_MemoryLocation memloc, double total_s) {
   if (rank) return;
   FILE *f = fopen(o->out, "w"); if (!f) { perror(o->out); return; }
   char tbuf[64]; time_t tt = time(NULL); strftime(tbuf, sizeof tbuf, "%Y-%m-%dT%H:%M:%S", localtime(&tt));
   const char *job = getenv("SLURM_JOB_ID"); char host[256]; gethostname(host, sizeof host);
   fprintf(f, "{\n  \"library\": \"hypre\", \"hypre_version\": \"%s\", \"hypre_release_date\": \"%s\", \"status\": \"%s\", \"label\": \"%s\", \"timestamp\": \"%s\", \"host\": \"%s\", \"slurm_job\": %s%s%s,\n",
           HYPRE_RELEASE_VERSION, HYPRE_RELEASE_DATE, status, o->label, tbuf, host, job ? "\"" : "", job ? job : "null", job ? "\"" : "");
   fprintf(f, "  \"build\": {\"cuda\": %s, \"mixedint\": %s, \"bigint\": %s, \"unified_memory\": %s, \"sizeof_HYPRE_Int\": %d, \"sizeof_HYPRE_BigInt\": %d},\n",
#if defined(HYPRE_USING_CUDA)
           "true",
#else
           "false",
#endif
#if defined(HYPRE_MIXEDINT)
           "true",
#else
           "false",
#endif
#if defined(HYPRE_BIGINT)
           "true",
#else
           "false",
#endif
#if defined(HYPRE_USING_UNIFIED_MEMORY)
           "true",
#else
           "false",
#endif
           (int) sizeof(HYPRE_Int), (int) sizeof(HYPRE_BigInt));
   fprintf(f, "  \"instance\": \"%s\", \"grid\": [%d, %d, %d], \"h\": %.10g, \"nu\": %lld, \"nnz\": %lld, \"nnz_per_row\": %.4f,\n", o->bin, I->nx, I->ny, I->nz, I->h, I->nu, I->nnz, (double) I->nnz / I->nu);
   fprintf(f, "  \"execution\": \"%s\", \"assembly\": \"%s\", \"mpi_ranks\": %d, \"storage\": \"full (ParCSR diag+offd, f64 values, int32 local indices)\", \"fine_operator_gb\": %.3f,\n",
           memloc == HYPRE_MEMORY_DEVICE ? "device" : "host", o->ij ? "HYPRE_IJMatrixSetValues2 (device arrays)" : "direct hypre_ParCSRMatrix fill + hypre_ParCSRMatrixMigrate", nprocs, (12.0 * I->nnz + 4.0 * I->nu) / 1e9);
   fprintf(f, "  \"settings\": {\"solver\": \"PCG (two-norm, rel. residual)\", \"tol\": %.3g, \"max_iter\": %d, \"num_functions\": %d, \"nodal\": %d, \"nodal_diag\": %d, \"interp_vectors\": %d, \"interp_vec_variant\": %d, \"interp_vec_qmax\": %d, "
           "\"coarsen_type\": %d, \"interp_type\": %d, \"relax_type\": %d, \"num_sweeps\": %d, \"agg_num_levels\": %d, \"agg_interp_type\": %d, \"p_max_elmts\": %d, \"strong_threshold\": %.3f, \"max_levels\": %d, \"max_coarse_size\": %d, \"keep_transpose\": %d, \"mod_rap2\": %d, \"spmv_use_vendor\": %d, \"spgemm_use_vendor\": %d, \"cheby_order\": %d, \"cheby_fraction\": %.3f, \"cheby_variant\": %d},\n",
           o->tol, o->maxit, o->num_functions, o->nodal, o->nodal_diag, o->interp_vectors, o->interp_vec_variant, o->interp_vec_qmax, o->coarsen, o->interp, o->relax, o->num_sweeps, o->agg_levels, o->agg_interp, o->pmax, o->theta, o->max_levels, o->max_coarse, o->keep_transpose, o->mod_rap2, o->spmv_vendor, o->spgemm_vendor, o->cheby_order, o->cheby_fraction, o->cheby_variant);
   json_res(f, "jacobi_cg", Rj, 0); fprintf(f, ",\n"); json_res(f, "jacobi_cg_rough_rhs", Rr, 0); fprintf(f, ",\n"); json_res(f, "amg_cg", Ra, 1); fprintf(f, ",\n");
   fprintf(f, "  \"device_memory_total_gb\": %.3f, \"peak_device_memory_gb\": %.3f, \"device_memory_samples\": [", dev_total_gb, dev_peak_gb);
   for (int i = 0; i < nsamp; i++) fprintf(f, "%s[\"%s\", %.3f]", i ? ", " : "", samp[i].label, samp[i].used_gb);
   fprintf(f, "],\n  \"host_maxrss_gb_rank0\": %.3f, \"total_seconds\": %.3f\n}\n", maxrss_gb(), total_s);
   fclose(f);
   logp("wrote %s (%s)", o->out, status);
}

static void usage(void) {
   if (!rank) fprintf(stderr, "usage: hypre_boomer_pcg --bin PREFIX --out F.json [--label S] [--ij] [--host] [--no-amg] [--jacobi] [--rough-rhs]\n"
      "  [--num-functions 4] [--nodal 0] [--nodal-diag 0] [--interp-vectors 0|3|4] [--interp-vec-variant 2] [--interp-vec-qmax 0]\n"
      "  [--coarsen 8] [--interp 6] [--relax 18] [--sweeps 1] [--agg-levels 1] [--agg-interp 0(default)] [--pmax 4] [--theta 0.25] [--max-levels 25] [--max-coarse 9]\n"
      "  [--keep-transpose 1] [--mod-rap2 1] [--spmv-vendor 1] [--spgemm-vendor 0] [--tol 1e-8] [--maxit 500] [--amg-print 1] [--cheby-order 2 --cheby-fraction 0.3 --cheby-variant 0 (with --relax 16)]\n");
   MPI_Finalize(); exit(2);
}

int main(int argc, char **argv) {
   MPI_Init(&argc, &argv); MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &nprocs); T0 = now();
   opts_t o = { .bin = NULL, .out = NULL, .label = "", .ij = 0, .host = 0, .jacobi = 0, .amg = 1, .rough_rhs = 0,
                .num_functions = 4, .nodal = 0, .nodal_diag = 0, .interp_vectors = 0, .interp_vec_variant = 2, .interp_vec_qmax = 0,
                .coarsen = 8, .interp = 6, .relax = 18, .agg_levels = 1, .agg_interp = 0, .pmax = 4, .max_levels = 25, .max_coarse = 9, .num_sweeps = 1, .keep_transpose = 1, .mod_rap2 = 1,
                .theta = 0.25, .tol = 1e-8, .maxit = 500, .spmv_vendor = 1, .spgemm_vendor = 0, .amg_print = 1, .cheby_order = 2, .cheby_variant = 0, .cheby_fraction = 0.3 };
   for (int i = 1; i < argc; i++) {
#define OPT_S(name, field) if (!strcmp(argv[i], name) && i + 1 < argc) { o.field = argv[++i]; continue; }
#define OPT_I(name, field) if (!strcmp(argv[i], name) && i + 1 < argc) { o.field = atoi(argv[++i]); continue; }
#define OPT_D(name, field) if (!strcmp(argv[i], name) && i + 1 < argc) { o.field = atof(argv[++i]); continue; }
      OPT_S("--bin", bin) OPT_S("--out", out) OPT_S("--label", label)
      if (!strcmp(argv[i], "--ij")) { o.ij = 1; continue; } if (!strcmp(argv[i], "--host")) { o.host = 1; continue; }
      if (!strcmp(argv[i], "--no-amg")) { o.amg = 0; continue; } if (!strcmp(argv[i], "--jacobi")) { o.jacobi = 1; continue; } if (!strcmp(argv[i], "--rough-rhs")) { o.rough_rhs = 1; continue; }
      OPT_I("--num-functions", num_functions) OPT_I("--nodal", nodal) OPT_I("--nodal-diag", nodal_diag) OPT_I("--interp-vectors", interp_vectors) OPT_I("--interp-vec-variant", interp_vec_variant) OPT_I("--interp-vec-qmax", interp_vec_qmax)
      OPT_I("--coarsen", coarsen) OPT_I("--interp", interp) OPT_I("--relax", relax) OPT_I("--sweeps", num_sweeps) OPT_I("--agg-levels", agg_levels) OPT_I("--agg-interp", agg_interp) OPT_I("--pmax", pmax)
      OPT_I("--max-levels", max_levels) OPT_I("--max-coarse", max_coarse) OPT_I("--keep-transpose", keep_transpose) OPT_I("--mod-rap2", mod_rap2) OPT_I("--spmv-vendor", spmv_vendor) OPT_I("--spgemm-vendor", spgemm_vendor)
      OPT_I("--maxit", maxit) OPT_I("--amg-print", amg_print) OPT_D("--theta", theta) OPT_D("--tol", tol) OPT_I("--cheby-order", cheby_order) OPT_I("--cheby-variant", cheby_variant) OPT_D("--cheby-fraction", cheby_fraction)
      if (!rank) fprintf(stderr, "unknown option %s\n", argv[i]);
      usage();
   }
   if (!o.bin || !o.out) usage();
   HYPRE_Initialize();
   HYPRE_MemoryLocation memloc = o.host ? HYPRE_MEMORY_HOST : HYPRE_MEMORY_DEVICE;
#if !defined(HYPRE_USING_CUDA)
   memloc = HYPRE_MEMORY_HOST; o.host = 1;
#endif
   HYPRE_SetMemoryLocation(memloc); HYPRE_SetExecutionPolicy(o.host ? HYPRE_EXEC_HOST : HYPRE_EXEC_DEVICE);
#if defined(HYPRE_USING_CUDA)
   HYPRE_SetSpMVUseVendor(o.spmv_vendor); HYPRE_SetSpGemmUseVendor(o.spgemm_vendor);
#endif
   logp("hypre %s (%s), %s, %d MPI ranks, label '%s', bin %s", HYPRE_RELEASE_VERSION, HYPRE_RELEASE_DATE, o.host ? "HOST execution" : "DEVICE execution (CUDA)", nprocs, o.label, o.bin);
   sample("after HYPRE_Initialize");
   inst_t I; read_meta(o.bin, &I); double t = now(); load_local(o.bin, &I); logp("local rows read in %.1f s (host maxrss %.1f GB)", now() - t, maxrss_gb());
   t = now(); hypre_ParCSRMatrix *A = o.ij ? build_ij(&I, memloc) : build_direct(&I, memloc); logp("matrix assembled in %.1f s", now() - t);
   sample("matrix on device");
   hypre_ParVector *b = make_vec(&I, I.b, memloc), *x = make_vec(&I, NULL, memloc);
   res_t Rj = { 0 }, Rr = { 0 }, Ra = { 0 };
   write_json(&o, &I, "matrix_on_device", &Rj, &Rr, &Ra, memloc, now() - T0);
   /* matvec timing (fine operator) */
   { hypre_ParVector *y = make_vec(&I, NULL, memloc); hypre_ParCSRMatrixMatvec(1.0, A, b, 0.0, y); sync_all(); double tm = now(); for (int k = 0; k < 10; k++) hypre_ParCSRMatrixMatvec(1.0, A, b, 0.0, y); sync_all(); tm = (now() - tm) / 10;
     logp("matvec %.2f ms (%.0f GB/s effective on %.1f GB)", 1e3 * tm, (12.0 * I.nnz + 4.0 * I.nu) / 1e9 / tm, (12.0 * I.nnz + 4.0 * I.nu) / 1e9); hypre_ParVectorDestroy(y); }
   if (o.jacobi) { run_pcg(&o, &I, A, b, x, I.xt, 0, memloc, &Rj, "jacobi-cg"); write_json(&o, &I, "jacobi_done", &Rj, &Rr, &Ra, memloc, now() - T0); }
   if (o.rough_rhs) { hypre_ParVector *br = make_vec(&I, I.rough, memloc); run_pcg(&o, &I, A, br, x, NULL, 0, memloc, &Rr, "jacobi-cg/rough"); hypre_ParVectorDestroy(br); write_json(&o, &I, "rough_done", &Rj, &Rr, &Ra, memloc, now() - T0); }
   if (o.amg) { run_pcg(&o, &I, A, b, x, I.xt, 1, memloc, &Ra, "amg-cg"); }
   write_json(&o, &I, "done", &Rj, &Rr, &Ra, memloc, now() - T0);
   hypre_ParVectorDestroy(b); hypre_ParVectorDestroy(x); hypre_ParCSRMatrixDestroy(A);
   free(I.b); free(I.xt); free(I.rough); free(I.cand);
   HYPRE_Finalize(); MPI_Finalize(); return 0;
}
