# FeatFloWer Guide 05: Agent Workflow to Create a New Ready-to-Use Application

This guide is a practical checklist for agents to create a new FeatFloWer application under `applications/`, including optional PE coupling in serial and parallel modes.

It covers:
- Application scaffolding (Fortran + app-local CMake)
- PE bridge wiring (`commf2c_*`)
- PE setup function implementation (serial and optional parallel)
- Integration into top-level CMake
- Runtime staging target (including METIS library copy)
- Validation and commit workflow (including submodule handling)

## 1. Naming and Scope

Choose one canonical application name and use it consistently:
- Folder: `applications/q2p1_<name>`
- Executable target: `q2p1_<name>`
- Main program file: `q2p1_<name>.f90`
- PE bridge symbol: `commf2c_<name>` (Fortran call) / `commf2c_<name>_` (C symbol)
- PE serial setup function: `setup<Name>Serial(...)`
- Optional parallel setup function: `setup<Name>(MPI_Comm ex0)`

## 2. Scaffold the Application Folder

Use a similar app as template (for Q2/P1 + FBM/PE, `q2p1_creep`, `q2p1_fc_ext`, `q2p1_drill` are common references).

Typical files to create:
- `applications/q2p1_<name>/CMakeLists.txt`
- `applications/q2p1_<name>/app_init.f90`
- `applications/q2p1_<name>/q2p1_<name>.f90`
- `applications/q2p1_<name>/_data/q2p1_param.dat`
- Optional run/test helpers (`tests/`, start scripts)
- If PE serial setup loads JSON: `applications/q2p1_<name>/example.json`

## 3. Update Fortran Sources

## 3.1 Main Program

In `q2p1_<name>.f90`:
- Update `PROGRAM ...` and `END PROGRAM ...` names.
- Keep standard loop and calls unless the physics workflow differs.

## 3.2 app_init.f90

In `General_init_ext(...)`:
- Keep MPI group split excluding rank 0 (master) if app follows standard PE coupling pattern.
- In `#ifdef HAVE_PE`, call your bridge function on nonzero ranks:

```fortran
#ifdef HAVE_PE
 if (myid .ne. 0) then
   call commf2c_<name>(MPI_COMM_WORLD, MPI_Comm_Ex0, myid)
 end if
#endif
```

- If your app defines `myGDATNEW`, call it explicitly (avoid unresolved `GDATNEW`):

```fortran
CALL myGDATNEW(CSimPar,0)
```

## 4. Integrate Application into CMake

## 4.1 App-local CMakeLists

In `applications/q2p1_<name>/CMakeLists.txt`:
- Define source list with `app_init.f90` + `q2p1_<name>.f90`.
- Create executable `q2p1_<name>`.
- Link `${FF_APPLICATION_LIBS}`.
- Include `${FF_APPLICATION_INCLUDE_PATH}`.
- Set Fortran compile flags and linker language.
- Use `createDefaultDirectories(...)` for out-of-source builds.
- Copy only required runtime files (for example `example.json`, `_data` assets, mesh/config helpers) to binary dir.

Legacy items to skip for new applications:
- Do **not** add `add_test(...)` entries in the app CMake (CTest path is legacy).
- Do **not** copy legacy test/starter scripts by default:
  - `tests/`
  - `q2p1_fc_ext_start.py`
  - `q2p1_ctest_start.py`

## 4.2 Top-level applications registration

Add in `applications/CMakeLists.txt`:

```cmake
add_subdirectory(q2p1_<name>)
```

## 4.3 Add staging target (recommended)

For runtime convenience, add a stage target that copies METIS into the app binary dir:

```cmake
add_custom_target(q2p1_<name>_stage
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          ${CMAKE_BINARY_DIR}/extern/libraries/metis-4.0.3/Lib/libmetis.so
          ${CMAKE_CURRENT_BINARY_DIR}/libmetis.so
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS q2p1_<name> metis
  COMMENT "Staging q2p1_<name> runtime files in build directory"
)
```

Extend this target with additional runtime copies if needed (`_data/MG.dat`, `example.json`, partitioner tools, helper scripts).

## 5. Add PE Bridge Functions

Files involved:
- `libs/pe/pe/interface/c2f_interface.h`
- `libs/pe/src/interface/c2f_interface.cpp`

## 5.1 Header declaration

Add declaration (MPI signature):

```cpp
extern "C" void commf2c_<name>_(MPI_Fint *Fcomm, MPI_Fint *FcommEx0, int *remoteRank);
```

## 5.2 c2f_interface.cpp wiring

Add both branches:

1. Parallel PE branch (`#if HAVE_MPI`):
- Either call a real parallel setup function (see section 7), or
- Provide explicit runtime error stub if parallel mode is intentionally unsupported.

2. Serial PE branch (`#ifdef PE_SERIAL_MODE`):
- Add function that calls your serial setup:

```cpp
extern "C" void commf2c_<name>_(int *Fcomm, int *FcommEx0, int *remoteRank) {
  pe::setup<Name>Serial(*remoteRank);
}
```

## 6. Implement PE Serial Setup

Primary file:
- `libs/pe/pe/interface/sim_setup_serial.h`

Add inline function:

```cpp
inline void setup<Name>Serial(int cfd_rank) {
  pe::logging::Logger::setCustomRank(cfd_rank);
  SimulationConfig::loadFromFile("example.json");
  auto &config = SimulationConfig::getInstance();
  config.setCfdRank(cfd_rank);

  WorldID world = theWorld();
  CollisionSystemID cs = theCollisionSystem();
  applyOptionalLubricationParams(*cs, config);

  world->setGravity(config.getGravity());
  world->setLiquidSolid(true);
  world->setLiquidDensity(config.getFluidDensity());
  world->setViscosity(config.getFluidViscosity());
  world->setDamping(1.0);
  TimeStep::stepsize(config.getStepsize());

  // Create your geometry/particles here
}
```

Serial setup recommendations:
- Always set logger custom rank for per-domain PE log separation.
- Use `example.json` for adjustable setup parameters.
- Validate critical config values and throw clear `std::runtime_error` messages.
- Print concise setup summary on representative rank (`cfd_rank == 1`).

## 7. Optional: Implement True Parallel PE Setup

If parallel PE support is required (not just a stub), add/setup:

1. Declare function in `libs/pe/pe/interface/sim_setup.h`:

```cpp
void setup<Name>(MPI_Comm ex0);
```

2. Implement function in a new interface header under `libs/pe/pe/interface/`, usually:
- `setup_<name>.h` with `void setup<Name>(MPI_Comm ex0)`

3. Include that header from `libs/pe/src/interface/sim_setup.cpp`.

4. In `libs/pe/src/interface/c2f_interface.cpp` parallel branch, call `setup<Name>(CcommEx0)`.

Parallel setup should:
- Set PE MPI communicator (`mpisystem->setComm(ex0)`)
- Validate process grid dimensions from config
- Perform domain decomposition
- Create bodies with ownership/domain logic

## 8. Runtime File Integration

Ensure runtime assets are present in app binary directory:
- `_data/q2p1_param.dat`
- `_data/MG.dat` (if required)
- `example.json` (if serial PE setup loads it)
- Mesh/project assets (`_adc` links, project files)
- `libmetis.so` (via stage target)

If using out-of-source builds, verify copy rules under `IF(${OUT_OF_SOURCE_BUILD})`.
Avoid adding legacy CTest helper copies in new app templates.

## 9. Build and Verify

Recommended check sequence from repo root:

```bash
cmake -S . -B build-check \
  -DBUILD_APPLICATIONS=ON \
  -DBUILD_TESTING=OFF \
  -DUSE_PE=ON \
  -DUSE_PE_SERIAL_MODE=ON

cmake --build build-check --target q2p1_<name> -j8
cmake --build build-check --target q2p1_<name>_stage -j8
```

Validation criteria:
- `q2p1_<name>` links successfully.
- `q2p1_<name>_stage` exists and copies `libmetis.so` to app binary dir.
- No missing-symbol errors for `commf2c_<name>_`.

## 10. Commit Workflow (Important with `libs/pe` Submodule)

`libs/pe` is a git submodule. If you changed PE interface files:

1. Commit inside `libs/pe` first.
2. Commit superproject second, including updated submodule pointer plus application files.

This keeps application and PE interface changes reproducible.

## 11. Agent Checklist

Use this final checklist before handing over:
- [ ] App folder exists with renamed target/program/file names.
- [ ] `applications/CMakeLists.txt` includes `add_subdirectory(q2p1_<name>)`.
- [ ] `app_init.f90` calls `commf2c_<name>` in `#ifdef HAVE_PE` block.
- [ ] `c2f_interface.cpp` has both parallel and serial `commf2c_<name>_` branches.
- [ ] Serial PE setup exists in `sim_setup_serial.h` and compiles.
- [ ] `example.json` is present and copied (if serial setup loads it).
- [ ] Stage target `q2p1_<name>_stage` copies `libmetis.so`.
- [ ] Target builds and stage target runs.
- [ ] Submodule commit order handled correctly (if PE files changed).
