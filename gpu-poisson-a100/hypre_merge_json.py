#!/usr/bin/env python3
"""Merge the per-variant JSONs of hypre_boomer_pcg.c into one results/hypre_<instance>.json with the layout of the
CuPy JSONs (jacobi_cg, jacobi_cg_rough_rhs, amg_variants{label: ...}, amg_cg = primary variant, peak_device_memory_gb).
Usage: python3 hypre_merge_json.py --out results/hypre_p1_40x40x120.json [--primary LABEL] results/hypre_p1_40x40x120__*.json
"""
import argparse, json, sys

ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); ap.add_argument("--primary", default=None); ap.add_argument("files", nargs="+")
a = ap.parse_args()
runs = []
for f in a.files:
    try:
        d = json.load(open(f)); d["_file"] = f; runs.append(d)
    except Exception as ex:
        print(f"skip {f}: {ex!r}", file=sys.stderr)
if not runs: sys.exit("no readable JSONs")
base_keys = ["library", "hypre_version", "hypre_release_date", "build", "instance", "grid", "h", "nu", "nnz", "nnz_per_row", "storage", "fine_operator_gb", "device_memory_total_gb", "host"]
out = {k: runs[0].get(k) for k in base_keys}
out["slurm_jobs"] = sorted({str(r.get("slurm_job")) for r in runs})
out["timestamp"] = max(r.get("timestamp", "") for r in runs)
out["jacobi_cg"] = None; out["jacobi_cg_rough_rhs"] = None; out["amg_variants"] = {}
for r in runs:
    lab = r.get("label") or r["_file"]
    if r.get("jacobi_cg") and out["jacobi_cg"] is None:
        out["jacobi_cg"] = dict(r["jacobi_cg"], execution=r.get("execution"), mpi_ranks=r.get("mpi_ranks"), label=lab)
    if r.get("jacobi_cg_rough_rhs") and out["jacobi_cg_rough_rhs"] is None:
        out["jacobi_cg_rough_rhs"] = dict(r["jacobi_cg_rough_rhs"], label=lab)
    if r.get("settings") and (r.get("amg_cg") or r.get("status") != "done"):
        out["amg_variants"][lab] = dict(status=r.get("status"), execution=r.get("execution"), assembly=r.get("assembly"), mpi_ranks=r.get("mpi_ranks"),
                                        settings=r.get("settings"), cg=r.get("amg_cg"), peak_device_memory_gb=r.get("peak_device_memory_gb"),
                                        device_memory_samples=r.get("device_memory_samples"), host_maxrss_gb_rank0=r.get("host_maxrss_gb_rank0"), total_seconds=r.get("total_seconds"), file=r["_file"])
prim = a.primary or next((l for l, v in out["amg_variants"].items() if v.get("cg") and v.get("execution") == "device"), None)
out["primary_variant"] = prim
out["amg_cg"] = (out["amg_variants"].get(prim) or {}).get("cg") if prim else None
out["peak_device_memory_gb"] = max([r.get("peak_device_memory_gb") or 0 for r in runs] + [0])
json.dump(out, open(a.out, "w"), indent=1)
print(f"wrote {a.out}: {len(out['amg_variants'])} AMG variants, primary {prim}, jacobi {'yes' if out['jacobi_cg'] else 'no'}, peak device {out['peak_device_memory_gb']:.2f} GB")
