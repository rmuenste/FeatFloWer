# q2p1_sse Geometry Handling

Primary files:

- `source/src_util/ReadExtrud3DParameters.f90`
- `source/src_quadLS/QuadSc_Sigma_User.f90`
- `applications/q2p1_sse/app_init.f90`
- `applications/s3d_mesher/s3d_mesher.f90`

## Geometry Pipeline

```text
setup.e3d / Extrud3D.dat
  -> ReadS3Dfile
       fills mySigma machine, segment, material, and object data
  -> s3d_mesher
       creates _data/meshDir/file.prj when automatic meshing is enabled
  -> partitioner
       creates partition data for compute ranks
  -> General_init_ext
       reads/refines mesh, parametrizes boundary, smooths/deforms mesh
  -> CGAL/FBM setup
       registers STL/OFF geometry and particle-like rigid bodies
```

## Machine Types

The same executable is configured by `mySigma%cType`:

- `SSE`: single-screw geometry.
- `TSE`: twin-screw geometry.
- `DIE`: die geometry/gendie-style setup.
- `XSE`: broader release mode that accepts multiple extruder types.
- `NETZSCH`: accepted parser type; exact runtime coverage should be checked
  before refactoring.

## Segment And Object Handling

`ReadS3Dfile` parses the machine segment list from E3D sections. Segment data is
stored in `mySigma%mySegment`. The parser handles analytic screw elements and
file-backed geometry forms.

`Setup_STL_Segments` registers file-backed segment geometry. It assigns CGAL
indices and writes `mesh_names.offs` for segments with `ART` values:

- `STL`
- `STL_L`
- `STL_R`
- `STL_LR`

For left/right split geometry, separate left and right OFF lists receive
separate CGAL index arrays.

## Runtime Geometry Routines

`QuadSc_Sigma_User.f90` contains global extrusion geometry state and many
distance/shape routines. Examples include shell distance handling and element
shape routines for different screw/segment types.

Important implication: this module is both a data contract and an implementation
container. Any restructuring should first separate the stable state API from the
geometry algorithms, at least in documentation.

## Mesh And Boundary Steps

`General_init_ext` performs the mesh-side geometry work:

- `readTriCoarse`
- `refineMesh`
- coarse and refined communication setup
- `E013_CreateComm_coarse` and `E013_CreateComm`
- `E011_CreateComm`
- `InitParametrization_STRCT`
- `ProlongateParametrization_STRCT`
- `DeterminePointParametrization_STRCT`
- `ParametrizeBndryPoints_STRCT`
- `ProlongateCoordinates_SSE`
- `UmbrellaSmoother_STRCT`

These calls form the boundary between abstract E3D geometry and the actual FE
mesh coordinates used by the solver.

