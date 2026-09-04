#!/usr/bin/env python3
"""One-shot axis-aligned partition of a periodic box for a single-host run.

Runs tools/PyPartitioner.py with the axis-split method (-123: median split
along x, y and z -> 2x2x2 = 8 subgrids, NPart = 1 each) and flattens the
result into the single-host layout (sub0001/GRID0001..0008) with
flatten_axis_partition.py.  Meant for featflower_test partition stages,
which execute exactly one command without a shell.

Usage (from the run directory, like PyPartitioner.py):
  partition_periodic_box.py <MeshName> <ProjectFile> [nsub]
e.g.
  partition_periodic_box.py CHIBOX6 _data/CHIMERA_BOX6/box6.prj
"""

import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.dirname(HERE)


def main(argv):
    if len(argv) not in (2, 3):
        sys.exit(__doc__)
    mesh, prj = argv[0], argv[1]
    nsub = int(argv[2]) if len(argv) == 3 else 8
    method = {2: "-1", 4: "-12", 8: "-123"}.get(nsub)
    if method is None:
        sys.exit("nsub must be 2, 4 or 8 (axis split along 1, 2 or 3 axes)")
    cmd = [sys.executable, os.path.join(TOOLS, "PyPartitioner.py"), "1", method,
           str(nsub), mesh, prj]
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)
    subprocess.run([sys.executable, os.path.join(HERE, "flatten_axis_partition.py"),
                    os.path.join("_mesh", mesh)], check=True)


if __name__ == "__main__":
    main(sys.argv[1:])
