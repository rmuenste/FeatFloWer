# FeatFloWer Communication Structures (MPI)

This document explains the communication-related data structures built in
`source/src_mpi/parentcomm.f90` by `PARENTCOMM`, with explicit intent and
semantics for each variable when it is introduced. The goal is to help new
users connect mesh entities (vertices, elements, faces) to the parallel
communication layout.

---

## 1. Mesh Context (Entities and Connectivity)

The communication structures are derived from the core mesh arrays:

- `NVT`: number of vertices in the current mesh partition.
- `NEL`: number of elements in the current mesh partition.
- `NAT`: number of faces in the current mesh partition.
- `DCORVG(1:3, ivt)`: vertex coordinates for vertex `ivt`.
- `DCORAG(1:3, iface)`: face-center coordinates for face `iface`.
- `KVERT(1:8, iel)`: element-to-vertex connectivity (8 vertices per hex).
- `KAREA(1:6, iel)`: element-to-face connectivity (6 faces per hex).

`PARENTCOMM` maps local subdomain entities to a master (coarse) reference mesh
using coordinate matching via octrees.

---

## 2. High-Level Output of `PARENTCOMM`

`PARENTCOMM` builds three categories of mappings:

1. **Vertex mapping**: local vertex index ↔ master vertex index
2. **Element mapping**: local element index ↔ master element index
3. **Face-based neighbor structure**: for each neighbor rank, the interface
   faces and their owning elements and local face sides.

It also constructs a **vertex communication scheme** that counts shared
vertices between subdomains to size exchange buffers.

---

## 3. Core Communication Structures and Semantics

### 3.1 Global container `coarse` (master mapping storage)

`coarse` is a module-level structure (defined in `PP3D_MPI`) used as a shared
container for all mapping tables.

- `coarse%pElem`  
  **Intent:** total number of elements in the master mesh.  
  **Semantics:** value equals `NEL` on the master rank.

- `coarse%pFace`  
  **Intent:** total number of faces in the master mesh.  
  **Semantics:** value equals `NAT` on the master rank.

- `coarse%pVert`  
  **Intent:** total number of vertices in the master mesh.  
  **Semantics:** value equals `NVT` on the master rank.

- `coarse%pNEL(pID)`  
  **Intent:** number of elements owned by subdomain `pID`.  
  **Semantics:** filled on master from `MPI_Allgather` of `NEL`.

- `coarse%pNVT(pID)`  
  **Intent:** number of vertices owned by subdomain `pID`.  
  **Semantics:** filled on master from `MPI_Allgather` of `NVT`.

- `coarse%pDX(1:coarse%pFace)`  
  **Intent:** per-face scratch array for communication (sizes based on master
  face count).  
  **Semantics:** allocated but not filled in `PARENTCOMM` itself.

- `coarse%DX(1:coarse%pFace)`  
  **Intent:** local scratch array for non-master ranks.  
  **Semantics:** allocated but not filled in `PARENTCOMM` itself.

#### Vertex mapping tables

- `coarse%myVERTLINK(1:NVT)`  
  **Intent:** local-to-master vertex map for the current rank.  
  **Semantics:** `coarse%myVERTLINK(ivt_local) = ivt_master` or `0` if unmatched.

- `coarse%pVERTLINK(pID, 1:coarse%pVert)`  
  **Intent:** master-side table storing all subdomains' vertex mappings.  
  **Semantics:** for subdomain `pID`, the entry is the master vertex index that
  matches local vertex `ivt_local` at that position.

#### Element mapping tables

- `coarse%myELEMLINK(1:NEL)`  
  **Intent:** local-to-master element map for the current rank.  
  **Semantics:** `coarse%myELEMLINK(iel_local) = iel_master` or `0` if unmatched.

- `coarse%pELEMLINK(pID, 1:coarse%pElem)`  
  **Intent:** master-side table storing all subdomains' element mappings.  
  **Semantics:** for subdomain `pID`, the entry is the master element index
  that matches the local element at that position.

---

### 3.2 Temporary coordinate buffers (matching by position)

These arrays exist only inside `PARENTCOMM` to match entities by coordinates.

- `dCOORDINATES(1:3, 1:count)`  
  **Intent:** master-side reference coordinates to build the octree.  
  **Semantics:** holds master `DCORVG` (vertices), element centroids, or
  `DCORAG` (face centers) depending on the stage.

- `pCOORDINATES(1:3, 1:NEL)`  
  **Intent:** non-master local element centroids.  
  **Semantics:** centroid of each local element computed from `KVERT`.

---

### 3.3 Vertex communication scheme

`PARENTCOMM` builds a vertex sharing graph to support later MPI exchanges.

- `xComm(1:coarse%pVert, 1:subnodes)`  
  **Intent:** master-side boolean map of which global vertex is present in
  which subdomain.  
  **Semantics:** `xComm(ivt_master, pID) = .TRUE.` if subdomain `pID` owns a
  local vertex mapped to master vertex `ivt_master`.

- `VerticeCommunicationScheme(1:subnodes)`  
  **Intent:** per-rank counts of shared vertices with every other rank.  
  **Semantics:** for the current rank, entry `pJD` is the number of master
  vertices that are shared with rank `pJD`. This is sent to each subdomain.

---

### 3.4 Face mapping and neighbor interface structure

The face mapping is built using face centers and a periodic-aware octree.

#### Local face mapping arrays

- `myFACELINK(1:2, 1:pNAT)`  
  **Intent:** map master face indices to local faces and owner rank.  
  **Semantics:** for master face index `J`, `myFACELINK(1,J)` is the local face
  index (or 0), and `myFACELINK(2,J)` is the owning rank id.

- `myFACEPINK(1:2, 1:pNAT)`  
  **Intent:** collective sum across subdomains to identify shared faces.  
  **Semantics:** after `MPI_Allreduce`, `myFACEPINK` contains the sum of
  `myFACELINK` across ranks; mismatches indicate shared faces.

- `nFACELISTS(1:subnodes)`  
  **Intent:** per-rank count of interface faces with each neighbor.  
  **Semantics:** `nFACELISTS(pJD)` is the number of local faces shared with
  neighbor rank `pJD`.

#### Main neighbor structure: `mg_mpi(1)%parST`

`mg_mpi` is the MPI communication structure used by the multigrid framework.
`PARENTCOMM` fills its neighbor list for the finest level (`mg_mpi(1)`).

- `mg_mpi(1)%NeighNum`  
  **Intent:** number of neighboring ranks sharing a face interface.  
  **Semantics:** count of nonzero entries in `nFACELISTS`.

- `mg_mpi(1)%UE(1:NAT)`  
  **Intent:** per-face scratch/work array for later MPI exchanges.  
  **Semantics:** allocated here; contents filled elsewhere.

Each neighbor has a `parST(pID)` record:

- `mg_mpi(1)%parST(pID)%Neigh`  
  **Intent:** neighbor rank id (relative index within subnodes).  
  **Semantics:** the rank to exchange data with.

- `mg_mpi(1)%parST(pID)%Num`  
  **Intent:** number of interface faces with this neighbor.  
  **Semantics:** equals `nFACELISTS(Neigh)` at allocation time.

- `mg_mpi(1)%parST(pID)%FaceLink(1:2, 1:Num)`  
  **Intent:** face indices associated with the neighbor interface.  
  **Semantics:** `FaceLink(1,k)` is the local face index; `FaceLink(2,k)` holds
  a second face index (currently set equal to the local face).

- `mg_mpi(1)%parST(pID)%ElemLink(1:2, 1:Num)`  
  **Intent:** owning elements for the interface faces.  
  **Semantics:** `ElemLink(1,k)` is the local element containing
  `FaceLink(1,k)`, and `ElemLink(2,k)` is the local element containing
  `FaceLink(2,k)`.

- `mg_mpi(1)%parST(pID)%SideLink(1:Num)`  
  **Intent:** local face side index (1..6) inside the owning element.  
  **Semantics:** `SideLink(k)` is the `KAREA` face index for `ElemLink(1,k)`.

- `mg_mpi(1)%parST(pID)%SDVect(1:Num)`  
  **Intent:** send buffer for scalar data on interface faces.  
  **Semantics:** allocated here; filled during communication steps.

- `mg_mpi(1)%parST(pID)%RDVect(1:Num)`  
  **Intent:** receive buffer for scalar data on interface faces.  
  **Semantics:** allocated here; filled during communication steps.

- `mg_mpi(1)%parST(pID)%SVVect(1:Num)`  
  **Intent:** send buffer for vector data on interface faces.  
  **Semantics:** allocated here; filled during communication steps.

- `mg_mpi(1)%parST(pID)%RVVect(1:Num)`  
  **Intent:** receive buffer for vector data on interface faces.  
  **Semantics:** allocated here; filled during communication steps.

- `mg_mpi(1)%parST(pID)%PE(1:Num)`  
  **Intent:** per-interface-face scratch array (often for pressure/flux data).  
  **Semantics:** allocated here; populated by later MPI exchange routines.

---

## 4. Execution Flow Summary (What happens when `PARENTCOMM` runs)

1. **Vertex stage**  
   - Master broadcasts all vertex coordinates (`DCORVG`).  
   - Each rank matches its local vertices to master vertices.  
   - Master gathers `coarse%pVERTLINK` and builds vertex sharing counts.

2. **Element stage**  
   - Master computes element centroids and broadcasts them.  
   - Each rank matches local elements to master elements.  
   - Master gathers `coarse%pELEMLINK` and element counts.

3. **Face stage**  
   - Master broadcasts face centers (`DCORAG`).  
   - Each rank matches local faces to master faces (with periodic matching).  
   - Shared faces are detected; neighbor lists are created.  
   - For each neighbor, interface faces are linked to local elements and sides
     via `KAREA`.

---

## 5. Notes on Matching Strategy

- **Octree matching** is used for fast coordinate lookups.
- **`DEpsPrec`** is the geometric tolerance for a match.
- **Periodic faces** are included via `FindInPeriodicOctTree` with
  `dPeriodicity`.

---

## 6. Where This Is Used Next

The structures created in `PARENTCOMM` are later used to:

- size and allocate send/receive buffers,
- exchange per-face, per-vertex, and per-element data between neighbors,
- ensure consistent boundary conditions across subdomain interfaces.

For implementation details, see:
`source/src_mpi/parentcomm.f90`.

---

## 7. `CREATECOMM` (Multigrid-Level Communication)

`CREATECOMM` (in `source/src_mpi/pp3d_mpi.f90`) builds the communication
structures for a *refined multigrid level* `ILEV` by expanding the neighbor
interface from `ILEV-1`. It is called on non-master ranks and performs three
main tasks:

1. **Refine interface faces**: each coarse interface face becomes 4 faces.
2. **Rebuild local face↔element/side mapping** for the refined level.
3. **Align face ordering across neighbors** using face-center coordinates.
4. **Optionally build vertex exchange maps** (when `BLIN = .TRUE.`).

### 7.1 Inputs and base structures

- `ILEV`  
  **Intent:** current multigrid level being constructed.  
  **Semantics:** level `ILEV` is derived from `ILEV-1`.

- `NAT`, `NEL`, `NVT`  
  **Intent:** number of faces/elements/vertices at this level.  
  **Semantics:** define local array extents for coordinates and connectivity.

- `DCORAG(1:3, iface)`  
  **Intent:** face center coordinates used for interface matching.  
  **Semantics:** used to identify equivalent faces across neighbors.

- `DCORVG(1:3, ivt)`  
  **Intent:** vertex coordinates, used for optional vertex exchange mapping.  
  **Semantics:** used in `BLIN` path to match shared vertices across ranks.

- `KAREA(1:6, iel)`  
  **Intent:** element-to-face connectivity.  
  **Semantics:** maps each element face to a global face index.

- `KADJ(1:6, iel)`  
  **Intent:** element adjacency.  
  **Semantics:** used to walk across elements and build 2x2 refined patches.

- `KVERT(1:8, iel)`  
  **Intent:** element-to-vertex connectivity.  
  **Semantics:** used to extract vertices of interface faces in `BLIN` path.

- `BLIN`  
  **Intent:** toggle for vertex-based communication setup.  
  **Semantics:** if `.TRUE.`, builds E011 vertex maps; if `.FALSE.`, returns
  after face-based setup.

### 7.2 Refining the neighbor interface (`mg_mpi(ILEV)`)

`CREATECOMM` starts from `mg_mpi(ILEV-1)` and builds `mg_mpi(ILEV)`:

- `mg_mpi(ILEV)%NeighNum`  
  **Intent:** number of neighbor ranks at this level.  
  **Semantics:** copied from `mg_mpi(ILEV-1)%NeighNum`.

- `mg_mpi(ILEV)%UE(1:NAT)`  
  **Intent:** per-face scratch array for interface exchanges.  
  **Semantics:** allocated here; filled in later exchange routines.

For each neighbor `pID`, the interface list is refined:

- `nSIZE = mg_mpi(ILEV-1)%parST(pID)%Num * 4`  
  **Intent:** refined face count (each coarse face splits into 4).  
  **Semantics:** new `parST(pID)` arrays are sized to `nSIZE`.

Allocated arrays for `mg_mpi(ILEV)%parST(pID)`:

- `ElemLink(1:2, 1:nSIZE)`  
  **Intent:** local element indices for each refined interface face.  
  **Semantics:** `ElemLink(1,:)` is local ordering; `ElemLink(2,:)` is aligned
  to neighbor ordering after matching.

- `FaceLink(1:2, 1:nSIZE)`  
  **Intent:** face indices for each refined interface face.  
  **Semantics:** `FaceLink(1,:)` local face indices; `FaceLink(2,:)` neighbor-
  aligned ordering.

- `SideLink(1:nSIZE)`  
  **Intent:** local face index (1..6) inside the owning element.  
  **Semantics:** used to reconstruct face vertices and adjacency.

- `CoragLinkX/Y/Z(1:2, 1:nSIZE)`  
  **Intent:** face center coordinates for local (index 1) and neighbor (index 2)
  orderings.  
  **Semantics:** populated before/after exchange to match face ordering.

- `PE`, `SDVect`, `RDVect`, `SVVect`, `RVVect`  
  **Intent:** communication buffers for scalar/vector exchange.  
  **Semantics:** allocated here; populated by later routines.

- `ElemLin1/2`, `FaceLin1/2`  
  **Intent:** sorted versions of element/face lists.  
  **Semantics:** allocated for later ordering/sorting steps.

### 7.3 Refinement pattern (coarse face → 4 fine faces)

For each coarse interface face (from `mg_mpi(ILEV-1)`):

- `IEL = ElemLink(1,I)`  
  **Intent:** local element containing the coarse interface face.  

- `IAT = SideLink(I)`  
  **Intent:** local face number (1..6) on `IEL`.

Using `KADJ`, the code constructs a 2x2 patch of elements across faces 3 and 6:
`IEL1 .. IEL8`. Based on `IAT`, it selects four refined faces:

- `(IE1,IA1,IM1)`, `(IE2,IA2,IM2)`, `(IE3,IA3,IM3)`, `(IE4,IA4,IM4)`  
  **Intent:** each triple defines a refined interface face:
  `IEk` = owning element, `IAk` = local face side, `IMk` = global face id.

These are written into `ElemLink(1,:)`, `FaceLink(1,:)`, and `SideLink(:)`.

### 7.4 Sorting and neighbor alignment by face centers

After local construction:

- `SORTALL` is called on `FaceLink(1,:)` with associated arrays  
  **Intent:** stable local ordering of interface faces.

Then each rank exchanges `CoragLinkX/Y/Z(1,:)` with neighbors and fills
`CoragLinkX/Y/Z(2,:)` from received data. Matching is performed by coordinate
comparison with periodic offsets:

- **Match rule:**  
  `(abs(local - remote) < DEpsPrec)` or  
  `abs(abs(local - remote) - dPeriodicity(axis)) < DEpsPrec`.

When a match is found:

- `FaceLink(2,J) = FaceLink(1,I)`  
  **Intent:** neighbor-aligned face ordering.

- `ElemLink(2,J) = ElemLink(1,I)`  
  **Intent:** neighbor-aligned element ordering.

### 7.5 Optional vertex exchange map (`BLIN = .TRUE.`)

If `BLIN` is enabled, a vertex-based communication structure is built for
linear-element operations (E011-style exchanges).

Per neighbor:

- `CoorSP(pJD)%dCoor(1:3, 1:Num)`  
  **Intent:** coordinates of local interface vertices shared with neighbor.  

- `CoorSP(pJD)%iCoor(1:Num)`  
  **Intent:** corresponding local vertex indices.

These lists are exchanged, and matching vertices are detected by coordinate
equality within `DEpsPrec`. The result is stored in:

- `E011ST(pID)%VertLink(1:2, 1:Num)`  
  **Intent:** local vertex mapping for E011 exchange.  
  **Semantics:** both columns store the matching local vertex index; used as a
  stable ordering for send/receive buffers.

Buffers for E011 exchange are also allocated:

- `E011ST(pID)%SVVect/RVVect`  
  **Intent:** vector data exchange on shared vertices.

- `E011ST(pID)%SDVect/RDVect`  
  **Intent:** scalar data exchange on shared vertices.

- `E011ST(pID)%SBVect/RBVect`  
  **Intent:** logical (boolean) data exchange on shared vertices.

These buffers are consumed in routines like `E011Sum`, `E011_ABSMIN`, and
`E011Knpr`.

---

## 8. Related Documentation

- Mesh data structures and connectivity:
  `docs/md_docs/mesh_structure.md`
