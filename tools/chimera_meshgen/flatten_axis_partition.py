#!/usr/bin/env python3
"""Flatten an axis-aligned subgrid partition into the single-host GRID000k layout.

Periodic runs need an axis-uniform Cartesian partition (subdomain faces
must coincide with their periodic images; parentcomm.f90
CheckFaceClaimDecode).  tools/PyPartitioner.py produces such a split only
as SUBGRIDS (method -123 with NSubPart = 8 gives a 2x2x2 median split):

    python3 tools/PyPartitioner.py 1 -123 8 CHIBOX6 _data/CHIMERA_BOX6/box6.prj

which writes  _mesh/CHIBOX6/sub000k/GRID.tri + <face>.par  for k = 1..8.
The solver maps ranks to subgrids by HOST NAME (get_pid.f90), so on a
single host every worker reads sub0001/GRID000<rank>.tri.  This script
rewrites the subgrid layout into that flat form:

    _mesh/<folder>/sub0001/GRID000k.tri        := sub000k/GRID.tri
    _mesh/<folder>/sub0001/<face>_000k.par     := sub000k/<face>.par
    _mesh/<folder>/sub0001/GRID.tri, <face>.par, GRID.prj := top-level copies

and removes sub0002..sub000N.  Run it once after the partitioner; the
deck then uses SimPar@SubMeshNumber = 1 and -np N+1.

Usage:
  flatten_axis_partition.py _mesh/CHIBOX6
"""

import os
import shutil
import sys


def main(argv):
    if len(argv) != 1:
        sys.exit(__doc__)
    root = argv[0]
    subs = sorted(d for d in os.listdir(root)
                  if d.startswith("sub") and os.path.isdir(os.path.join(root, d)))
    if not subs:
        sys.exit("no sub* directories under %s" % root)
    pars = sorted(f for f in os.listdir(root) if f.endswith(".par"))
    if not os.path.exists(os.path.join(root, "GRID.tri")):
        sys.exit("missing %s/GRID.tri" % root)

    # collect the per-subgrid coarse meshes first (sub0001 is rewritten)
    staged = []
    for k, sub in enumerate(subs, start=1):
        src = os.path.join(root, sub)
        tri = os.path.join(src, "GRID.tri")
        if not os.path.exists(tri):
            sys.exit("missing %s" % tri)
        with open(tri, "rb") as f:
            tri_bytes = f.read()
        par_bytes = {}
        for p in pars:
            with open(os.path.join(src, p), "rb") as f:
                par_bytes[p] = f.read()
        staged.append((k, tri_bytes, par_bytes))

    flat = os.path.join(root, "sub0001")
    for sub in subs:
        shutil.rmtree(os.path.join(root, sub))
    os.makedirs(flat)
    for k, tri_bytes, par_bytes in staged:
        with open(os.path.join(flat, "GRID%04d.tri" % k), "wb") as f:
            f.write(tri_bytes)
        for p, b in par_bytes.items():
            with open(os.path.join(flat, "%s_%04d.par" % (p[:-4], k)), "wb") as f:
                f.write(b)
    shutil.copy(os.path.join(root, "GRID.tri"), os.path.join(flat, "GRID.tri"))
    shutil.copy(os.path.join(root, "GRID.prj"), os.path.join(flat, "GRID.prj"))
    for p in pars:
        shutil.copy(os.path.join(root, p), os.path.join(flat, p))
    print("flattened %d subgrids of %s into sub0001/GRID0001..%04d" % (len(staged), root, len(staged)))


if __name__ == "__main__":
    main(sys.argv[1:])
