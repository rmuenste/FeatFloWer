# Reproducing the DKT mesh set — agent instructions

Audience: an agent (or person) who needs to rebuild `dkt_mesh_v1/` from nothing,
or adapt it to a different diameter, box, or resolution ladder.

Everything here was executed and verified on 2026-08-02. Every expected number
below is a real observed value, not an estimate — if your run disagrees with one
of them, stop and find out why rather than continuing.

---

## 0. What you are building

A closed box **6d × 6d × 24d** (d = sphere diameter) for the drafting–kissing–
tumbling benchmark, as a FeatFloWer coarse mesh plus Cartesian partitionings that
match PE's `decomposeDomain`.

| artefact | what | size |
|---|---|---|
| `coarseB/` | 12×12×48 = 6,912 el, coarse cell = d/2 | 62 MB |
| `coarseB/_mesh/DKT_B_108` | `x3-y3-z12`, 108 subdomains × 64 coarse el | |
| `coarseB/_mesh/DKT_B_864` | `x6-y6-z24`, 864 subdomains × 8 coarse el | |
| `coarseC/` | 24×24×96 = 55,296 el, coarse cell = d/4 | 82 MB |
| `coarseC/_mesh/DKT_C_864` | `x6-y6-z24`, 864 subdomains × 64 coarse el | |

Total 25,747 files. `coarseC` is a fallback for the case where 8 coarse elements
per subdomain proves too thin; `coarseB` is the primary path.

**No custom script is needed.** An earlier plan called for a `make_box_tri.py`;
`pe_partpy/unit_cube_vtk.py` already is exactly that.

---

## 1. Environment

Two things bite here. Both are silent-ish failures if you skip them.

```bash
module load python/3.13.5
```

**Required for every shell that touches the venv.** The venv's `bin/python`
dynamically links `libpython3.13.so.1.0`; without the module you get
`error while loading shared libraries`. The system `python3` is 3.9.25, too old
for the pinned `numpy==2.3.1`.

If the venv does not exist yet:

```bash
GL2=<repo>/meshing_component/gmsh-learning2
uv venv --python /sfw/python/3.13.5/bin/python $GL2/.venv
uv pip install --python $GL2/.venv/bin/python \
    -r $GL2/gmsh-learning/requirements.txt vtk pytest
```

`vtk` and `pytest` are **not** in `requirements.txt`, but `unit_cube_vtk.py`
(stage 1) and the test suite need them.

`gmsh` is **not** needed — this is the pure-`pe_partpy` box pipeline, no `.geo`
involved. (`module load gmsh/4.13.1` is only for the `.geo` generators.)

Verify: `$GL2/workflow/check_env.sh` should print `=== OK ===`.

### pe_partpy version requirement

You need `pe_partpy` at or after the `axis_uniform` boundary-precision fix
(rmuenste/pe_partpy#1). Before that fix, `axis_uniform` used a span-relative
tolerance that silently mis-partitions meshes whose `.tri` was written with few
decimals. The meshes built here are full precision and would be unaffected, but
a regenerated or hand-edited mesh could trip it, and the failure is silent —
it still prints `The partitioning was successful!`.

Quick check that you have the fix:

```bash
grep -n "_snap_to_levels" $GL2/pe_partpy/partitioner/part.py    # present => fixed
```

---

## 2. Decide the parameters before running anything

### 2.1 Diameter `d`

The committed meshes are **nondimensional with d = 1**, so bounds are
x, y ∈ [−3, 3] and z ∈ [0, 24].

To use a physical `d`, **regenerate with scaled bounds**. Do *not* scale the
`.tri` by hand: the `.par` files carry quoted type-4 plane descriptors
(`'4 1d0 0d0 0d0 3d0'`) that encode the same bounds, and they would silently
disagree with the mesh.

### 2.2 Coarse resolution

`D/h` counts **Q2 elements** per diameter, `h = dvol^(1/3)` — the same
convention as `DNS_RESOLUTION` and `ComputeCFL`. (Velocity nodal spacing is
`h/2`; LBM/IBM papers quote nodal, so halve their numbers before comparing.)

`{8, 16, 32}` is a doubling sequence and regular refinement halves `h`, so this
is **one coarse mesh at three consecutive levels**, not three meshes. Pick the
coarse cell so the target `D/h` lands on the levels you intend to run:

| coarse cell | coarse mesh | D/h=8 | 16 | 32 |
|---|---|---|---|---|
| d/2 (**B**) | 12×12×48 | L3 | L4 | L5 |
| d/4 (**C**) | 24×24×96 | L2 | L3 | L4 |

### 2.3 Partition split

Two hard constraints and one soft one:

1. **Each per-axis count must divide the corresponding cell count** — on the
   coarse mesh *and*, implicitly, on every refined level (refinement preserves
   divisibility, so checking the coarse mesh is enough).
2. **The `x?-y?-z?` spec must equal the `(px, py, pz)` passed to PE's
   `decomposeDomain`.** Subdomain numbering loops x outermost, z innermost:
   `subId = ((ix-1)*partY + (iy-1))*partZ + iz` — **z varies fastest**.
3. Soft: keep subdomains close to cubic. For this 6:6:24 box, exact AR = 1.0
   requires `px = py = k, pz = 4k`, i.e. **ranks = 4k³** → 4, 32, 108, 256, 864,
   2048. Nothing exists between 256 and 864.

> **`x5-y5-z20` looks valid and is not.** 6/5 = 24/20 = 1.2 so the aspect is 1.0,
> but the fine element counts are `6·D/h = 2^(m+1)·3`, which never has a factor 5
> — at D/h = 32 the mesh is 192×192×768 and 192/5 = 38.4. Any split must divide
> the *element counts*, not just look right geometrically.

Helper for picking:

```bash
$GL2/.venv/bin/python $GL2/.claude/skills/box-mesh-pipeline/scripts/pick_partition_dims.py \
    --nx 12 --ny 12 --nz 48 --total 108
```

Note it optimises for cubic subdomains, which is *not* always what you want —
if PE is configured for pure z-slabs you must match that instead.

---

## 3. Build

`$GL2` = `<repo>/meshing_component/gmsh-learning2`, `$D` = target folder.

### Coarse B

```bash
module load python/3.13.5
mkdir -p $D/coarseB && cd $D

$GL2/.venv/bin/python $GL2/pe_partpy/unit_cube_vtk.py --out /tmp/dkt_B.vtk \
  --nx 12 --ny 12 --nz 48 --x0 -3 --x1 3 --y0 -3 --y1 3 --z0 0 --z1 24

$GL2/.venv/bin/python $GL2/pe_partpy/tri2vtk_converter.py /tmp/dkt_B.vtk /tmp/dkt_B.tri

$GL2/.venv/bin/python $GL2/pe_partpy/gen_par_from_tri.py /tmp/dkt_B.tri --outdir $D/coarseB

$GL2/.venv/bin/python $GL2/workflow/tag_box_boundaries.py --case-dir $D/coarseB
```

Stage notes:

- `unit_cube_vtk.py` writes **only** `.vtk`; it cannot emit `.tri` directly.
- `tri2vtk_converter.py` goes **both** directions, chosen from the file
  extensions of its two positional arguments, despite the name.
- `gen_par_from_tri.py` is the **axis-aligned** detector — correct for a box.
  Do not use `gen_par_from_tri_regions.py` (region growing, for curved geometry)
  or `..._by_normals.py` here. It copies the `.tri` into `--outdir` and writes
  `file.prj`, so the source `.tri` can live in `/tmp`.
- `tag_box_boundaries.py` rewrites the six `.par` files in place with quoted
  type-4 plane descriptors and btype `Wall`. **The quoting is load-bearing**: the
  CFD reads that line with a Fortran list-directed READ, which stops at the first
  blank in an unquoted string.

### Partition

```bash
cd $D/coarseB
$GL2/.venv/bin/python $GL2/pe_partpy/PyPartitioner.py 108 axis_uniform x3-y3-z12 DKT_B_108 file.prj -f v2
$GL2/.venv/bin/python $GL2/pe_partpy/PyPartitioner.py 864 axis_uniform x6-y6-z24 DKT_B_864 file.prj -f v2
```

- **`-f v2` is required.** It produces `sub%04d/GRID%04d.tri`, which is what
  `source/include/PartitionReader.f90` builds. The default `v1` gives 3-digit
  names the reader will not find.
- **Use the strategy name.** Numeric codes (`-4`) are rejected; you get a message
  naming the replacement.
- **Run `PyPartitioner.py` directly, from the case folder.** It writes `_mesh/`
  relative to the current directory, which is what you want. Do *not* use
  `workflow/run_partition_to_vtu.py` for this: it invokes
  `pe_partpy/PyPartitioner.py` by *relative* path so it only works with cwd =
  `$GL2`, and it always writes to `$GL2/_mesh/`, not next to your case.

### Coarse C

Identical, with `--nx 24 --ny 24 --nz 96`, output name `dkt_C`, folder
`coarseC`, and only the 864 split (`x6-y6-z24`, mesh name `DKT_C_864`).

### Optional: ParaView file

```bash
$GL2/.venv/bin/python $GL2/pe_partpy/tri2vtk_converter.py $D/coarseB/file.prj -proj $D/coarseB
# writes $D/coarseB/main.vtu
```

---

## 4. Verify — expected values

### Coarse meshes

```bash
sed -n 3p $D/coarseB/dkt_B.tri
sed -n 3p $D/coarseC/dkt_C.tri
```

```
 6912  8281 0 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE
55296 60625 0 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE
```

8281 = 13·13·49, 60625 = 25·25·97. Boundary node counts:

| | xmin/xmax/ymin/ymax | zmin/zmax |
|---|---|---|
| B | 637 (= 13·49) | 169 (= 13·13) |
| C | 2425 (= 25·97) | 625 (= 25·25) |

Bounds must be x,y ∈ [−3,3], z ∈ [0,24]. Check with:

```bash
awk '/DCORVG/{f=1;next} /KVERT/{f=0} f{print}' $D/coarseB/dkt_B.tri | \
awk '{if(NR==1){x1=x2=$1;y1=y2=$2;z1=z2=$3}
      if($1<x1)x1=$1; if($1>x2)x2=$1; if($2<y1)y1=$2
      if($2>y2)y2=$2; if($3<z1)z1=$3; if($3>z2)z2=$3}
     END{printf "x[%g,%g] y[%g,%g] z[%g,%g]\n",x1,x2,y1,y2,z1,z2}'
```

> Keep this as **two** awk stages. Collapsing it into one breaks `NR==1`, which
> then matches the file header instead of the first coordinate line, leaving the
> min/max uninitialised and printing plausible-looking garbage. (This bit me.)

### Partitions

```bash
for m in DKT_B_108:108:64 DKT_B_864:864:8; do
  IFS=: read name n exp <<< "$m"
  tot=0; cnt=0; bad=0
  for d in $D/coarseB/_mesh/$name/sub*/; do
    c=$(sed -n 3p ${d}GRID0001.tri | awk '{print $1}')
    tot=$((tot+c)); cnt=$((cnt+1)); [ "$c" != "$exp" ] && bad=$((bad+1))
  done
  echo "$name: $cnt subs, $tot el, expected $exp each, uneven=$bad"
done
```

Expected — `108 subs, 6912 el, 64 each, uneven=0` and `864 subs, 6912 el, 8 each,
uneven=0`; `DKT_C_864` gives `864 subs, 55296 el, 64 each, uneven=0`.

### Subdomain geometry and ordering

`DKT_B_864` subdomains must be 1.0³ cubes, `DKT_B_108` 2.0³:

```
sub0001  x[-3.00,-2.00] y[-3.00,-2.00] z[0.00,1.00]
sub0002  x[-3.00,-2.00] y[-3.00,-2.00] z[1.00,2.00]   <- z stepped: z is fastest
sub0025  x[-3.00,-2.00] y[-2.00,-1.00] z[0.00,1.00]   <- y stepped after 24 z-slabs
sub0864  x[ 2.00, 3.00] y[ 2.00, 3.00] z[23.00,24.00]
```

That `sub0002` steps z and `sub0025` steps y is the check that ordering matches
PE's `decomposeDomain(6,6,24)`. If your build steps x first, the CFD and PE
decompositions will disagree and particles will be attributed to the wrong ranks.

### Wall times (this machine, for sanity)

| step | wall |
|---|---|
| coarse B partition, 108 | 4.4 s |
| coarse B partition, 864 | 32.3 s |
| coarse C partition, 864 | 39.4 s |

---

## 5. Things that will trip you

**`NBCT = 0` and all-zero `KNPR`.** These meshes differ from the DeViSoR
`benchSym/mesh.tri` (which has `NBCT = 1` and 652 nonzero KNPR). Both are
**benign** and were traced:

- `%nbct` has no consumer — read at `mesh_refine.f90:245`, echoed back at `:37`,
  hard-set to 1 in the JSON path at `:1007`. FeatFloWer writes `NBCT = 0` itself
  at `:101`.
- `mesh_refine.f90:1492` does `mesh1%knpr = 0` on every refined level, so coarse
  KNPR is discarded at the first refinement. Walls come from the `.par` files via
  `Parametrization`.

Do not "fix" these to match benchSym.

**KNPR record layout does matter.** `readTriCoarse` reads it with one `read` per
vertex (`do i=1,nvt; read(cunit,*) knpr(i)`), so values must be **one per
record**. A single long line would satisfy the first read and then hit EOF. The
toolkit writes one per line; verify with:

```bash
f=$D/coarseB/dkt_B.tri
ln=$(grep -n KNPR $f | head -1 | cut -d: -f1)
echo "$(( $(wc -l < $f) - ln )) records vs NVT $(sed -n 3p $f | awk '{print $2}')"
```

**Partition counts that do not divide.** `axis_uniform` now refuses two cases
that used to produce structurally-valid but broken output: an element matching no
slab, and two boundaries collapsing onto one mesh level (an empty subdomain —
which still consumes a rank, so the CFD decomposition stops matching PE's). If
you see `ValueError: ... at least one subdomain would be empty`, your rank count
exceeds what the mesh resolution supports on that axis.

**Disk and file count.** 864 subdomains × ~14 files each. The full set is 25,747
files / 144 MB. Partitioning 864 ways takes ~30–40 s and creates a lot of inodes;
do not run it somewhere with a tight quota.

**Mesh folders are untracked.** `dkt_mesh_v1/` is deliberately not in git, same as
`benchSym/`, `ten_cate_mesh_v1/`, `pipemesh_v1/`. Do not add it.

---

## 6. Adapting

| change | what to do |
|---|---|
| different `d` | rescale all six `--x0…--z1` bounds; **regenerate**, never edit the `.tri` |
| different box aspect | recompute the AR = 1.0 family: `px=py=k, pz=(Lz/Lx)·k` |
| different `D/h` ladder | pick the coarse cell so the targets land on levels you will run; keep the ladder a doubling sequence or you need more than one mesh |
| pure z-slabs for PE | use `x1-y1-zN`; N must divide the z cell count |
| periodic instead of closed | use `workflow/tag_pipe_boundaries.py` as the pattern (`Periodic` btype, coupling by coordinate matching with `dPeriodicity`) instead of `tag_box_boundaries.py` |

For the reasoning behind the specific choices here — why 864 ranks, why coarse C
exists, the weak-scaling property, and the decomposition-invariance check worth
running at D/h = 16 — see `README.md` in this folder.
