# Containerizing FeatFloWer

This guide summarizes how to build and use the bundled Docker image for FeatFloWer. Images are based on Ubuntu 24.04, install GCC/GFortran, OpenMPI, BLAS/LAPACK, and Boost (including `program_options`); submodules are fetched over HTTPS as configured in `.gitmodules`, so no SSH keys are needed during the build. The image defaults to the portable build id `generic-linux-gcc-release` so CMake does not rely on CPU autodetection inside the container.

## Build the Image

From the repository root, run:

```bash
docker build -t featflower:latest .
```

Pass additional CMake configuration via build arguments, for example:

```bash
docker build -t featflower:debug \
  --build-arg CMAKE_BUILD_TYPE=Debug \
  --build-arg BUILD_APPLICATIONS=OFF \
  --build-arg Q2P1_BUILD_ID=generic-linux-gcc-debug .
```

For a production-ready release image with all applications prebuilt:

```bash
docker build -t featflower:latest \
  --build-arg CMAKE_BUILD_TYPE=Release \
  --build-arg BUILD_APPLICATIONS=ON \
  --build-arg Q2P1_BUILD_ID=generic-linux-gcc-release .
```

## Run Interactive Sessions

Launch a shell inside the container with your project directory mounted:

```bash
docker run --rm -it \
  -v "$(pwd)":/workspace \
  featflower:latest
```

Mounted volumes inherit host permissions; adjust `USERNAME/UID/GID` build arguments if needed to match your user. The Dockerfile adapts to existing identities: if the chosen GID or UID are already in use (for example the `ubuntu` user shipped with some base images), it renames that account/group to the requested `USERNAME` instead of failing, so overriding `USER_UID`/`USER_GID` to match your host remains safe.

Use the following run workflow:
- Start: `docker run --rm -it -v "$(pwd)":/workspace featflower:latest`
- Inside the container, applications live in `/workspace/build/applications`. Launch them directly or copy binaries out as needed.
- Shutdown: type `exit` or press `Ctrl+D`; the container stops immediately because runs use `--rm`.

## Executing MPI Jobs

For small-scale runs you can invoke MPI inside the container directly:

```bash
mpirun -np 4 ./applications/q2p1_devel/q2p1_devel
```

For multi-node deployments, consider using container-aware MPI launchers provided by your cluster environment.

## Running Tests

Images are built with CTest enabled. From within the container:

```bash
cmake --build build --target test
ctest --test-dir build --output-on-failure
```

If you need to rebuild after modifying the source on the host, run the configure and build commands inside the container again.

## Image Cleanup

List images with `docker images featflower` and remove unused ones using `docker image rm featflower:<tag>`. Because containers are started with `--rm`, they disappear automatically once you exit.
