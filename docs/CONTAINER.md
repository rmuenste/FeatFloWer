# Containerization Guide (Docker/Podman)

This image builds FeatFloWer with GCC/OpenMP/Fortran, OpenMPI, and CGAL (via CMake ExternalProject). It uses low-noise CMake config and logs build output to a file.

## Build Image

- Default repo/ref used inside the image:
  - `FF_REPO=https://github.com/rmuenste/FeatFloWer.git`
  - `FF_REF=master`
- Build:
  - `docker build -t featflower:latest .`
  - Or pin a ref: `docker build --build-arg FF_REF=<tag-or-commit> -t featflower:latest .`

## Partition a Mesh

- Mount your working dir and meshes; set `Q2P1_MESH_DIR`:
  - `docker run --rm -v "$PWD":/work -w /work -v "$MESHDIR":/meshes -e Q2P1_MESH_DIR=/meshes featflower:latest \
    python /opt/featflower/tools/PyPartitioner.py 3 1 1 NEWFAC _adc/2D_FAC/2Dbench.prj`

## Run an Application (quiet)

- Example (q2p1_fc_ext):
  - `docker run --rm -v "$PWD":/work -w /work featflower:latest \
     bash -lc "mpirun -np 4 /opt/featflower/applications/q2p1_fc_ext/q2p1_fc_ext > run.log 2>&1 || tail -n 300 run.log"`

## Validate a Mesh with CGAL Tool

- `docker run --rm -v "$PWD":/work -w /work featflower:latest \
   bash -lc "validate_triangle_mesh /opt/featflower/tools/check_manifold/cube.off > validate.log 2>&1 && tail -n 200 validate.log"`

## Podman Notes

- Works with the same Dockerfile: `podman build -t featflower:latest .`
- Use `--userns=keep-id` and, on SELinux systems, `--security-opt label=disable` if volume permissions are restrictive.

