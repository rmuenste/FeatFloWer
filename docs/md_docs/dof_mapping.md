## Subroutine NDFGL 
This subroutine is located in extern/feat3d/src/ndfgl.f
This is the final piece of the core puzzle! `NDFGL` is the glue between the mesh data (`KVERT`, `KEDGE`, `KAREA`) and the basis function definitions (`E012`, `E013`). This code tells us exactly how the degrees of freedom (DOFs) are numbered and associated with the mesh entities.

Let's break down what `NDFGL` tells us:

**1. General DOF Numbering Scheme:**
The code reveals a standard hierarchical numbering scheme for degrees of freedom in a complex FEM code:
*   **Vertex DOFs:** Global DOF numbers for vertices are simply their global vertex numbers. There are `NVT` of these, numbered 1 to `NVT`.
*   **Edge DOFs:** Global DOF numbers for edges are their global edge numbers offset by `NVT`. There are `NET` of these, numbered `NVT+1` to `NVT+NET`.
*   **Face (Area) DOFs:** Global DOF numbers for faces are their global face numbers offset by `NVT+NET`. There are `NAT` of these, numbered `NVT+NET+1` to `NVT+NET+NAT`.
*   **Element (Volume/Cell) DOFs:** Some elements (like Q2) have DOFs associated with the element interior (the "center of gravity"). These are numbered by their element number `IEL` offset by `NVT+NET+NAT`.

This is confirmed by the block starting `DO 1 IVE=1,NVE ...`:
```fortran
! --- Collects DOFs associated with Vertices, Edges, and Faces ---
DO 1 IVE=1,NVE
  JVG(IVE)=KVERT(IVE,IEL) ! Global vertex number
  JVL(IVE)=IVE            ! Local vertex number (1 to 8)
ENDDO
...
DO 2 IEE=1,NEE
  JVG(NKE+IEE)=KEDGE(IEE,IEL)+NVT ! Global edge number offset by NVT
  JVL(NKE+IEE)=NKE+IEE          ! Local edge number (9 to 20)
ENDDO
...
DO 3 IAE=1,NAE
  JVG(NKE+IAE)=KAREA(IAE,IEL)+NVT+NET ! Global face number offset by NVT+NET
  JVL(NKE+IAE)=NKE+IAE              ! Local face number (21 to 26)
ENDDO
```
*   `JVG` collects the global DOF numbers.
*   `JVL` collects the corresponding local DOF numbers.

**2. Handling of Specific Element Types (`IELTYP`):**
The `GOTO` statement directs the flow based on `IELTYP`. This confirms our previous suspicion about how different element types are handled.

*   **`IELTYP = 12` (P1 Discontinuous Pressure):**
    ```fortran
    KDFG(1)= 4*(IEL-1)+1
    KDFG(2)= 4*(IEL-1)+2
    ...
    ```
    This is extremely telling. The global DOFs for a P1 discontinuous element are not associated with any shared mesh entity (vertex, edge, face). They are numbered `4*IEL-3, 4*IEL-2, 4*IEL-1, 4*IEL`, meaning each element gets its own private block of 4 DOFs. This is the definition of a discontinuous element space. The `P` array in `GetForces` is indexed this way. `KDFL` is just `1,2,3,4`.

*   **`IELTYP = 13` (Q2 Continuous Velocity):**
    The `GOTO` directs this to `label 130`.
    ```fortran
    130   CONTINUE
          DO 131 IKE=1,NKE
    131   KDFG(IKE)=JVG(IKE)
    C *** Dof corresponding to the center of gravity
          KDFG(NKE+1)=NVT+NET+NAT+IEL
    ```
    This assembles the DOFs for the 27-node Q2 element:
    *   `JVG` contains the 8 vertex DOFs, 12 edge DOFs, and 6 face DOFs (total `NKE=26`). These are copied to `KDFG`.
    *   An additional DOF is added for the element's interior (the 27th node): its global number is `NVT+NET+NAT+IEL`.
    *   This perfectly constructs the list of 27 global DOF numbers for the Q2 element, which are then used to index into the global `U1,U2,U3` arrays.
    *   `KDFL` gets the corresponding local numbers `JVL` (1 to 26) plus the 27th for the center.

**3. The `NGLS` Subroutine and Sorting (`IPAR`):**
*   `NGLS` is a simple insertion sort ("Piksort").
*   It sorts the `JVG` array (global DOFs) and applies the same swaps to the `JVL` array (local DOFs) to maintain the correct correspondence.
*   `IPAR` controls whether this sorting happens. In our routines (`GetForces`, etc.), `IPAR=1`, which means sorting is *not* done because `IPAR.GE.0` is false.
    *Correction*: The calls in `GetForces` are `NDFGL(IEL,1,IELTYP,...)`. `IPAR=1`, so `IPAR.GE.0` is TRUE. **Sorting IS performed**.
*   **Why sort?** Sorting `KDFG` can be very important for cache efficiency and for constructing sparse matrices in other parts of the FEM solver (e.g., the assembly of the stiffness matrix). When looping through `KDFG` to add contributions, accessing memory in a monotonically increasing way is much faster. For our post-processing `GetForces` routine, the impact is less critical but it's good practice.

**4. The `KDFL` array:**
*   This array maps the sorted `KDFG` back to the original local node numbering scheme used by the element routines (`E013A`).
*   For example, before sorting, you might have:
    `JVG = [100, 50, 200]` (Global DOFs)
    `JVL = [  1,  2,   3]` (Local DOFs)
*   After sorting (done by `NGLS`):
    `JVG` becomes `[50, 100, 200]`
    `JVL` becomes `[ 2,   1,   3]`
*   Now, when the main loop processes the first entry:
    *   `IG = KDFG(1)` will be `50`.
    *   `KDFL(1)` will be `2`.
    *   The `DBAS` array is accessed using `KDFL(1)`, so it gets the basis function for local node 2 (`DBAS(..., 2, ...)`).
    *   The global data arrays `U1, U2, U3` are accessed using `IG`, so `U1(50)`.
    This correctly associates the data for global node 50 with the basis function for local node 2.

**Summary of What `NDFGL` Tells Us:**

1.  **DOF Numbering is Hierarchical:** The global DOF indices are structured based on entity type: vertices first, then edges, then faces, then element interiors.
2.  **Discontinuous vs. Continuous Elements:** The routine explicitly handles different numbering for discontinuous elements (like P1 pressure, `IELTYP=12`) where DOFs are local to the element, and continuous elements (like Q2 velocity, `IELTYP=13`) where DOFs are shared on common vertices, edges, and faces.
3.  **Q2 Element Assembly is Confirmed:** The logic for `IELTYP=13` confirms our understanding of how the 27 DOFs (8 vertex + 12 edge + 6 face + 1 cell) are assembled from the mesh connectivity data.
4.  **`KDFL` is the "Remapping Key":** The `KDFL` array is crucial. It ensures that after the global DOF list `KDFG` is sorted for efficiency, we can still find the correct corresponding local basis function shape from `DBAS`.

This routine provides the final, low-level link between the abstract element definition (`E013`) and the concrete mesh data (`KVERT`, etc.), showing exactly how they are connected to form the global solution space. We now have a complete top-to-bottom picture of the core FEM evaluation loop.