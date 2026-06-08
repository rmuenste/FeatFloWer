# Containerization Guide

The repository Dockerfile builds FeatFloWer and its applications in an
Ubuntu 24.04 image. The same Dockerfile can be used with Docker or Podman.

The image contains:

- GCC, G++, and GFortran
- CMake and Ninja
- OpenMPI
- BLAS and LAPACK development libraries
- Boost and `program_options`
- Python 3
- the repository source at `/workspace`
- the configured build tree at `/workspace/build`

The default image does not enable optional PE, CGAL, HYPRE, or other external
features unless the Dockerfile is extended with additional CMake arguments and
dependencies.

## Current Image Model

The Dockerfile copies the current checkout into the image, initializes
submodules, configures the project, and builds it during `docker build`.

This differs from older container workflows that:

- cloned `FF_REPO` at `FF_REF`,
- installed into `/opt/featflower`,
- used a separate runtime stage, or
- enabled CGAL and mesh-validation tools by default.

Those paths and build arguments are not part of the current Dockerfile.

## Build the Image

Run from the repository root:

```bash
docker build -t featflower:latest .
```

The supported build arguments are:

- `CMAKE_BUILD_TYPE` (default: `Release`)
- `BUILD_APPLICATIONS` (default: `ON`)
- `Q2P1_BUILD_ID` (default: `generic-linux-gcc-release`)
- `USERNAME` (default: `featflower`)
- `USER_UID` (default: `1000`)
- `USER_GID` (default: `1000`)

To match the container user to the current host user:

```bash
docker build -t featflower:latest \
  --build-arg USER_UID="$(id -u)" \
  --build-arg USER_GID="$(id -g)" .
```

To build only the libraries:

```bash
docker build -t featflower:libs \
  --build-arg BUILD_APPLICATIONS=OFF .
```

The image is built from the current working tree. To build a particular tag or
commit, check out that revision on the host before running `docker build`.

## Important Build-Context Note

The repository currently has no `.dockerignore`. Docker therefore sends the
entire working tree as build context, including local build directories and
other untracked files.

For smaller and more reproducible image builds, keep large generated
directories outside the checkout or add an appropriate `.dockerignore` in a
separate change.

## Start an Interactive Container

Use the source and binaries already embedded in the image:

```bash
docker run --rm -it featflower:latest
```

Inside the container:

```text
/workspace                  source tree
/workspace/build            build tree
/workspace/build/applications
```

Do not mount the host repository over `/workspace` when you intend to use the
prebuilt binaries. Such a mount hides both the source and build tree embedded
in the image.

Mount input or output data at a separate path instead:

```bash
mkdir -p container-work

docker run --rm -it \
  -v "$PWD/container-work":/work \
  -w /work \
  featflower:latest
```

## Rebuild Host Source Inside the Container

To compile the current host checkout rather than use the embedded build, mount
it at a different path:

```bash
docker run --rm -it \
  -v "$PWD":/src \
  -w /src \
  featflower:latest \
  bash
```

Then configure a separate build directory:

```bash
cmake -S . -B build-container \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_APPLICATIONS=ON \
  -DQ2P1_BUILD_ID=generic-linux-gcc-release

cmake --build build-container -- -j"$(nproc)"
```

Because the image runs as the configured non-root user, generated host files
use that user's UID and GID.

## Mesh Repository and Partitioning

Mount the external mesh repository separately and set `Q2P1_MESH_DIR`:

```bash
docker run --rm -it \
  -v "$MESHDIR":/meshes:ro \
  -v "$PWD/container-work":/work \
  -e Q2P1_MESH_DIR=/meshes \
  featflower:latest \
  bash
```

The image build itself does not know the runtime mesh directory, so application
`_adc` links may need to be created or refreshed in the build application
directory:

```bash
ln -sfn /meshes /workspace/build/applications/q2p1_fc_ext/_adc
```

The maintained partitioner wrapper is documented in
[featflower_partitioner_usage_guide.md](featflower_partitioner_usage_guide.md).
The legacy script remains available at `/workspace/tools/PyPartitioner.py`.

Example from the application build directory:

```bash
cd /workspace/build/applications/q2p1_fc_ext
python3 /workspace/tools/PyPartitioner.py \
  3 1 1 NEWFAC _adc/2D_FAC/2Dbench.prj
```

The partitioner needs access to the built METIS library. Use the application's
staging target where one exists, or ensure `libmetis.so` is available in the
application directory or library search path.

## Run an Application

Applications are built under `/workspace/build/applications`.

Example:

```bash
cd /workspace/build/applications/q2p1_fc_ext
mpirun -np 4 ./q2p1_fc_ext
```

For a quiet run with a mounted output directory:

```bash
docker run --rm \
  -v "$PWD/container-work":/work \
  featflower:latest \
  bash -lc '
    cd /workspace/build/applications/q2p1_fc_ext
    mpirun -np 4 ./q2p1_fc_ext > /work/run.log 2>&1 ||
      { tail -n 300 /work/run.log; exit 1; }
  '
```

Small single-node MPI runs can execute inside the container. Multi-node runs
require the cluster's container-aware MPI integration and are outside the scope
of this image.

## Run Tests

CTest is enabled by the top-level CMake configuration:

```bash
docker run --rm featflower:latest \
  ctest --test-dir /workspace/build --output-on-failure
```

Some tests require external meshes, application runtime data, or MPI settings
that are not self-contained in the image.

## Optional CGAL Tools

The current Dockerfile does not configure `USE_CGAL=ON`, so
`validate_triangle_mesh` is not built in the default image.

To provide CGAL tooling, update the image dependencies and CMake invocation,
then build the `validate_triangle_mesh` target. See
`tools/check_manifold/README.md` for tool-specific usage.

## Podman

Build with:

```bash
podman build -t featflower:latest .
```

Run while preserving the host user mapping:

```bash
podman run --rm -it \
  --userns=keep-id \
  -v "$PWD/container-work":/work \
  -w /work \
  featflower:latest
```

On SELinux systems, add a suitable volume label such as `:Z`, or use
`--security-opt label=disable` when required by local policy.

## Cleanup

List and remove images with:

```bash
docker images featflower
docker image rm featflower:<tag>
```

Containers started with `--rm` are removed automatically when they exit.
