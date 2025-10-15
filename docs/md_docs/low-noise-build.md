# Low‑Noise Builds and Logging

This guide documents the recommended CMake options, parameters, and commands to configure and build FeatFloWer with minimal console noise while capturing full logs to files for later inspection.

## Configure (out‑of‑source, low noise)

Use Unix Makefiles and enable JSON diagnostics to keep output machine‑parsable and compact:

```bash
cmake -S . -B build -G "Unix Makefiles" \
  -DCMAKE_RULE_MESSAGES=OFF \
  -DCMAKE_MESSAGE_LOG_LEVEL=WARNING \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DCMAKE_C_FLAGS='-fdiagnostics-format=json' \
  -DCMAKE_CXX_FLAGS='-fdiagnostics-format=json' \
  -DCMAKE_Fortran_FLAGS='-fdiagnostics-format=json' \
  -DUSE_CGAL=ON -DUSE_CGAL_LOCAL=OFF
```

Notes
- Out‑of‑source only: use `-S . -B build` (in‑source is rejected).
- Generator: prefer `-G "Unix Makefiles"` to avoid Ninja’s Fortran module collisions.
- Message filtering: `CMAKE_RULE_MESSAGES=OFF` + `CMAKE_MESSAGE_LOG_LEVEL=WARNING` reduces CMake chatter.
- JSON diagnostics: `-fdiagnostics-format=json` on all languages keeps compiler output structured and compact.
- CGAL: `-DUSE_CGAL=ON -DUSE_CGAL_LOCAL=OFF` uses the ExternalProject flow (preferred).

Optional flags you may add as needed (examples):
- `-DUSE_PE=ON` to enable particle engine features if required by your target.
- `-DSHOW_BUILD_IDS=ON` to list supported toolchains.

## Build quietly and capture logs

Build with parallel jobs, redirect all output to a log, and only print the tail on failure:

```bash
# Build to a file inside the build tree
cmake --build build -j"$(nproc)" > build/build.log 2>&1 \
  || (tail -n 200 build/build.log && exit 1)
```

Alternatives and tips
- Pass `-- -s` to make for even less rule echoing: `cmake --build build -- -j"$(nproc)" -s`.
- You can also call make directly: `make -C build -j"$(nproc)" > build/build.log 2>&1`.
- Keep logs under `build/` so they don’t clutter the repo and remain per‑configuration.

## Running tests quietly

Run CTest with parallelism, capture to a log, and only show context on failure:

```bash
ctest --test-dir build -j8 > build/ctest.log 2>&1 \
  || (tail -n 200 build/ctest.log && exit 1)

# Filter to a subset
ctest --test-dir build -j8 -R "q2p1_.*" > build/ctest_q2p1.log 2>&1 \
  || (tail -n 200 build/ctest_q2p1.log && exit 1)
```

## Summarizing warnings/errors from logs

Because compiler diagnostics are JSON, you can quickly extract warnings/errors from the log:

```bash
# Using ripgrep if available (preferred)
rg '"(kind|level)":\s*"(warning|error)"' build/build.log || true

# Using grep as a fallback
grep -E '"(kind|level)":\s*"(warning|error)"' build/build.log || true
```

To inspect only the last entries after a long build:

```bash
tail -n 200 build/build.log
```

## Example: end‑to‑end low‑noise session

```bash
# 1) Configure
cmake -S . -B build -G "Unix Makefiles" \
  -DCMAKE_RULE_MESSAGES=OFF -DCMAKE_MESSAGE_LOG_LEVEL=WARNING \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DCMAKE_C_FLAGS='-fdiagnostics-format=json' \
  -DCMAKE_CXX_FLAGS='-fdiagnostics-format=json' \
  -DCMAKE_Fortran_FLAGS='-fdiagnostics-format=json' \
  -DUSE_CGAL=ON -DUSE_CGAL_LOCAL=OFF

# 2) Build quietly and log
cmake --build build -j"$(nproc)" > build/build.log 2>&1 \
  || (tail -n 200 build/build.log && exit 1)

# 3) (Optional) Run tests quietly
ctest --test-dir build -j8 > build/ctest.log 2>&1 \
  || (tail -n 200 build/ctest.log && exit 1)

# 4) Summarize issues (if any)
rg '"(kind|level)":\s*"(warning|error)"' build/build.log || true
```

## Ninja caveat (Fortran modules)

If you must use Ninja, ensure each Fortran module source is compiled by exactly one target and set per‑target module directories to avoid duplicate `.mod` outputs, e.g. in CMake:

```cmake
set_target_properties(<tgt> PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/<tgt>")
```

Until that is arranged, prefer `-G "Unix Makefiles"` for low‑noise builds involving multiple applications.

