# Mesh fixture (untracked, like the DNS DKT fixture)

`mesh/coarseB/*` and `mesh/partitions/DKT_B_108` are symlinks into
`<main checkout>/dkt_mesh_v1/coarseB` — the DNS campaign's 6d × 6d × 24d closed
column (12 × 12 × 48 coarse elements, coarse cell d/2, d = 1, x,y ∈ [−3,3],
z ∈ [0,24], all six faces `Wall` with quoted type-4 plane descriptors) and its
axis-uniform x3-y3-z12 partition (108 subdomains). Record and reproduction
recipe: `dkt_mesh_v1/README.md`, `dkt_mesh_v1/REPRODUCE.md`. `.gitignore`
excludes `*.tri`/`*.par`, so the fixture is not committed; recreate the links
with the commands in RUNBOOK.md section 2 on a fresh clone.
