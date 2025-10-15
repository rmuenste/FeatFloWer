module evaluations
  USE var_QuadScalar,ONLY:knvt,knet,knat,knel,tMultiMesh
  !-------------------------------------------------------------------------------------------------
  ! A module for computing the force acting on a cylinder     
  ! 
  ! 
  !-------------------------------------------------------------------------------------------------

contains
!************************************************************************
      SUBROUTINE GetVelocityAtElementCentroid(IEL_TARGET, U1,U2,U3, &
     &                                    KVERT,DCORVG,ELE, VelCentroid)
!
! Purpose: Calculates the velocity (U1, U2, U3) at the centroid of
!          a specified hexahedral element (IEL_TARGET).
!
! ELE should point to the routine for evaluating Q2 basis functions (e.g., E013).
! VelCentroid is an array of 3 REAL*8 values for the output velocity.
!************************************************************************
      USE PP3D_MPI, ONLY:myid ! Only for potential debug messages if needed
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
!
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,&
                 NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
!
      INTEGER IEL_TARGET
      REAL*8  U1(*),U2(*),U3(*)
      REAL*8  DCORVG(NNDIM,*)
      INTEGER KVERT(NNVE,*) ! KAREA, KEDGE not needed for centroid velocity
      INTEGER KDFG(NNBAS),KDFL(NNBAS)
!
      REAL*8  VelCentroid(NNDIM) ! Output: Velocity (U,V,W) at centroid
!
      REAL*8  XI1_centroid, XI2_centroid, XI3_centroid
!
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,&
                      DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),&
                      IEL,NDIM ! IEL will be set to IEL_TARGET
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,&
                      NVAR,NEAR,NBCT,NVBD,NEBD,NABD
!     COMMON /CUB/ not needed as we are not doing cubature integration here
      COMMON /COAUX1/ KDFG,KDFL,IDFL
!
! *** user COMMON blocks (copied for completeness, may not all be needed)
      INTEGER  VIPARM
      DIMENSION VIPARM(100)
      EQUIVALENCE (IAUSAV,VIPARM)
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,&
                    IMASS,IMASSL,IUPW,IPRECA,IPRECB,&
                    ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,&
                    INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,&
                    ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,&
                    IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA
 
      SAVE
C
      IF (ICHECK.GE.997) CALL OTRC('GetVelCentroid','NEW')
      IER = 0
      VelCentroid(1:NNDIM) = 0.0D0
!
! --- Define centroid coordinates in reference element [-1,1]^3 ---
      XI1_centroid = 0.0D0
      XI2_centroid = 0.0D0
      XI3_centroid = 0.0D0
!
! --- Initialize BDER to get only basis function VALUES ---
      DO I= 1,NNDER
        BDER(I)=.FALSE.
      END DO
      BDER(1)=.TRUE. ! We only need N_i (value), not derivatives
!
! --- Setup Element Type (assuming Q2 for velocity like E013) ---
      IELTYP = -1 ! Prepare for ELE to set this
      CALL ELE(0D0,0D0,0D0,IELTYP) ! ELE is E013, IELTYP becomes 13
      IDFL=NDFL(IELTYP) ! Should be 27 for Q2 elements
      IF (IDFL .NE. NNBAS .AND. IELTYP .EQ. 13) THEN
        CALL WERR(-1, 'GetVelCentroid: IDFL mismatch for E013')
        GOTO 99998
      END IF
!
! --- Set the current element number for COMMON /ELEM/ ---
      IEL = IEL_TARGET
      IF (IEL_TARGET .LT. 1 .OR. IEL_TARGET .GT. NEL) THEN
        CALL WERR(-1, 'GetVelCentroid: Invalid IEL_TARGET')
        GOTO 99998
      END IF
!
! --- Get DOF mapping for the target element ---
      CALL NDFGL(IEL_TARGET,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99998 ! Error in NDFGL
!
! --- Evaluate basis functions at the reference centroid (0,0,0) ---
! --- For IPAR=0, ELE (E013) will calculate DBAS directly.
! --- We don't need to pre-calculate Jacobian here because E013 with IPAR=0
! --- if it needed physical coordinates for some reason would calculate them.
! --- However, for just values, E013A output DHELP(...,1,...) is directly used.
! --- E013 with IPAR=0, when only BDER(1)=.TRUE., simply copies DHELP(IDFL,1,1)
! --- (calculated by E013A(xi1,xi2,xi3, DHELP,1)) to DBAS(1,IDFL,1).
! --- No Jacobian transformation is applied to *values*.
!
      CALL ELE(XI1_centroid, XI2_centroid, XI3_centroid, 0) ! IPAR = 0
      IF (IER.LT.0) GOTO 99998 ! Error in ELE
!
! --- Interpolate velocity components at the centroid ---
      DO I=1,IDFL
        IG = KDFG(I)
        ! DBAS(1,KDFL(I),1) contains the value of the KDFL(I)-th basis function
        ! evaluated at the reference centroid.
        ! Note: KDFL(I) is the local basis function number corresponding to
        ! the I-th DOF of the element.
        VelCentroid(1) = VelCentroid(1) + U1(IG) * DBAS(1,KDFL(I),1)
        VelCentroid(2) = VelCentroid(2) + U2(IG) * DBAS(1,KDFL(I),1)
        VelCentroid(3) = VelCentroid(3) + U3(IG) * DBAS(1,KDFL(I),1)
      END DO
!
      GOTO 99999 ! Normal exit
!
99998 CONTINUE ! Error exit
      VelCentroid(1:NNDIM) = -1.0D+30 ! Indicate error with a large negative value
      IF (myid .EQ. 0) WRITE(*,*) 'Error in GetVelCentroid, IER=', IER
C
99999 CONTINUE
      END SUBROUTINE GetVelocityAtElementCentroid

!How to Use:
!You would call this subroutine from your main program after you have the mesh information and the nodal velocity solution:
!
!REAL*8 MyVelAtCentroid(3)
!      INTEGER TargetElementID
!      ! ... set TargetElementID ...
!      ! ... U1, U2, U3, KVERT, DCORVG are available ...
!      CALL GetVelocityAtElementCentroid(TargetElementID, U1,U2,U3, &
!     &                                KVERT,DCORVG,E013, MyVelAtCentroid)
!      ! Now MyVelAtCentroid(1), (2), (3) contain ux, uy, uz at centroid
!      ! of element TargetElementID
!
!Make sure that the ELE argument in the call actually points to E013 (or whatever the Q2 element routine is named if it's different in the calling scope).
end module evaluations

!
!**The Challenge: Mapping Arbitrary Physical Point to Reference Coordinates**
!
!1.  **Input:** You'll provide:
!    *   The target element `IEL_TARGET`.
!    *   The **physical coordinates $(x_p, y_p, z_p)$** of the arbitrary point where you want to evaluate the velocity.
!    *   Nodal velocities `U1, U2, U3`, mesh data `KVERT, DCORVG`, element routine `ELE`.
!
!2.  **The Crucial Step:** Before we can use `ELE` (e.g., `E013`) with `IPAR=0` to evaluate basis functions, we need the coordinates of your arbitrary physical point $(x_p, y_p, z_p)$ **transformed into the reference element coordinates $(\xi_{1p}, \xi_{2p}, \xi_{3p})$ for that specific element `IEL_TARGET`**.
!
!    The mapping we have discussed so far (using `DJmn` coefficients and the Jacobian in `GetForces`) is for:
!    *   Reference $(\xi_1, \xi_2, \xi_3)$  $\rightarrow$ Physical $(x,y,z)$
!
!    We need the **inverse mapping**:
!    *   Physical $(x_p, y_p, z_p)$ (within `IEL_TARGET`) $\rightarrow$ Reference $(\xi_{1p}, \xi_{2p}, \xi_{3p})$
!
!    This inverse mapping for a general isoparametric element (especially non-linear like Q2 geometry, though our geometry seems trilinear) is non-trivial. It usually involves solving a system of non-linear equations iteratively, for example, using a Newton-Raphson method.
!
!    The equations to solve are:
!    $x_p = \sum_{k=1}^{N_{geom}} N_k^{geom}(\xi_1, \xi_2, \xi_3) x_k$
!    $y_p = \sum_{k=1}^{N_{geom}} N_k^{geom}(\xi_1, \xi_2, \xi_3) y_k$
!    $z_p = \sum_{k=1}^{N_{geom}} N_k^{geom}(\xi_1, \xi_2, \xi_3) z_k$
!    where $N_k^{geom}$ are the shape functions used for the geometry mapping (trilinear, 8-node in our case based on `DJmn` calculations). We need to find $(\xi_1, \xi_2, \xi_3)$ given $(x_p, y_p, z_p)$ and the nodal coordinates $(x_k, y_k, z_k)$ of the element `IEL_TARGET`.
!
!**Do we need additional information? Potentially, yes:**
!
!*   **Inverse Mapping Routine:** Does the existing codebase have a utility subroutine that performs this inverse mapping? (e.g., `FindReferenceCoords(IEL, Xphys, Xi_ref)`).
!    *   If **yes**, then the task is straightforward: call this routine first, then proceed as with the centroid calculation, but using the found $(\xi_{1p}, \xi_{2p}, \xi_{3p})$.
!    *   If **no**, we would need to implement one. For a trilinear hexahedron, the system is non-linear but manageable with Newton-Raphson.
!
!*   **Which Element Contains the Point?** The current request is "velocity at an arbitrary point *in a hexahedral element*". This implies you already know which element `IEL_TARGET` contains your physical point $(x_p, y_p, z_p)$. If you only have $(x_p, y_p, z_p)$ and need to *first find which element it's in*, that's an additional "point location" or "search" problem.
!
!**Let's assume for now:**
!
!1.  You **know** which element `IEL_TARGET` contains the physical point $(x_p, y_p, z_p)$.
!2.  We **need to implement or assume availability of an inverse mapping** function.
!
!**High-Level Plan for `GetVelocityAtArbitraryPoint`:**
!
!```fortran
!SUBROUTINE GetVelocityAtArbitraryPoint(IEL_TARGET, X_phys, U1,U2,U3, &
!                                       KVERT,DCORVG,ELE, VelAtPoint, IER_map)
!    ! IEL_TARGET: The element assumed to contain X_phys
!    ! X_phys(3): The physical coordinates (xp, yp, zp) of the point
!    ! U1,U2,U3: Nodal velocities
!    ! KVERT, DCORVG: Mesh data
!    ! ELE: Element routine (e.g., E013)
!    ! VelAtPoint(3): Output velocity
!    ! IER_map: Error code from inverse mapping (0 if successful)
!
!    IMPLICIT DOUBLE PRECISION ...
!    ! ... Declarations ...
!    REAL*8 Xi_ref(3) ! To store reference coordinates
!
!    ! 1. Perform Inverse Mapping:
!    !    Convert X_phys within IEL_TARGET to reference coordinates Xi_ref.
!    !    This is the crucial step. Let's assume a hypothetical function for now:
!    CALL PhysicalToRefCoords(IEL_TARGET, X_phys, KVERT, DCORVG, Xi_ref, IER_map)
!    IF (IER_map .NE. 0) THEN
!        ! Point might be outside the element, or convergence failed
!        VelAtPoint = HugeNegativeValue_or_NaN 
!        RETURN
!    END IF
!
!    ! Optional: Check if Xi_ref is within [-1,1]^3 (or [0,1]^3 depending on E013A)
!    ! If not, the point X_phys was not actually in IEL_TARGET, or mapping error.
!
!    ! 2. Initialize BDER for values only (same as centroid calculation)
!    DO I= 1,NNDER
!        BDER(I)=.FALSE.
!    END DO
!    BDER(1)=.TRUE.
!
!    ! 3. Setup Element Type (same as centroid)
!    IELTYP = -1
!    CALL ELE(0D0,0D0,0D0,IELTYP) 
!    IDFL=NDFL(IELTYP)
!    IEL = IEL_TARGET ! Set for common block
!
!    ! 4. Get DOF mapping (same as centroid)
!    CALL NDFGL(IEL_TARGET,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
!    IF (IER.LT.0) THEN /* Handle error */ VelAtPoint = ...; RETURN; END IF
!
!    ! 5. Evaluate Basis Functions at Xi_ref (same as centroid, but with Xi_ref)
!    CALL ELE(Xi_ref(1), Xi_ref(2), Xi_ref(3), 0) ! IPAR = 0
!    IF (IER.LT.0) THEN /* Handle error */ VelAtPoint = ...; RETURN; END IF
!
!    ! 6. Interpolate Velocity (same as centroid)
!    VelAtPoint(1:NNDIM) = 0.0D0
!    DO I=1,IDFL
!        IG = KDFG(I)
!        VelAtPoint(1) = VelAtPoint(1) + U1(IG) * DBAS(1,KDFL(I),1)
!        VelAtPoint(2) = VelAtPoint(2) + U2(IG) * DBAS(1,KDFL(I),1)
!        VelAtPoint(3) = VelAtPoint(3) + U3(IG) * DBAS(1,KDFL(I),1)
!    END DO
!
!END SUBROUTINE
!```
!
!**Implementation of `PhysicalToRefCoords` (Newton-Raphson for Trilinear Hex):**
!
!This is the non-trivial part if it doesn't exist. For a trilinear mapping:
!$x(\xi) = \sum N_k(\xi) x_k$
!$y(\xi) = \sum N_k(\xi) y_k$
!$z(\xi) = \sum N_k(\xi) z_k$
!
!We want to find $\xi = (\xi_1, \xi_2, \xi_3)$ such that:
!$f_1(\xi) = \sum N_k(\xi) x_k - x_p = 0$
!$f_2(\xi) = \sum N_k(\xi) y_k - y_p = 0$
!$f_3(\xi) = \sum N_k(\xi) z_k - z_p = 0$
!
!The Newton-Raphson iteration is:
!$\mathbf{J}_{\text{map}}(\xi^{(n)}) \Delta\xi^{(n)} = -\mathbf{f}(\xi^{(n)})$
!$\xi^{(n+1)} = \xi^{(n)} + \Delta\xi^{(n)}$
!
!Where $\mathbf{J}_{\text{map}}$ is the Jacobian of the *geometric mapping functions* $f_1, f_2, f_3$ with respect to $\xi_1, \xi_2, \xi_3$. This is precisely the `DJAC` matrix we calculate in `GetForces` (but derived from the 8-node trilinear geometry shape functions, not the 27-node Q2 field shape functions).
!
!The `DJmn` coefficients pre-calculated in `GetForces` are exactly what's needed to form this `DJAC` for the geometry mapping.
!
!```fortran
!SUBROUTINE PhysicalToRefCoords(IEL_TARGET, X_phys, KVERT_map, DCORVG_map, &
!                               Xi_ref_out, IER_map)
!    ! KVERT_map, DCORVG_map are KVERT, DCORVG
!    ! Uses Newton-Raphson to find reference coordinates Xi_ref_out
!    ! corresponding to physical point X_phys in element IEL_TARGET.
!    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!    PARAMETER (NNVE_geom=8, NNDIM_geom=3, MAX_ITER_MAP=20, TOL_MAP=1.D-8)
!    INTEGER IEL_TARGET, IER_map, ITER
!    REAL*8 X_phys(NNDIM_geom), Xi_ref_out(NNDIM_geom)
!    REAL*8 Xi_current(NNDIM_geom), F_vec(NNDIM_geom), Delta_Xi(NNDIM_geom)
!    REAL*8 Geom_Jacobian(NNDIM_geom,NNDIM_geom), Geom_Jacobian_inv(NNDIM_geom,NNDIM_geom)
!    REAL*8 X_calc(NNDIM_geom)
!    REAL*8 DX_geom(NNVE_geom), DY_geom(NNVE_geom), DZ_geom(NNVE_geom) ! Element nodal coords
!    REAL*8 DJ11, DJ12, ..., DJ83 ! Geometry mapping coefficients
!
!    INTEGER KVERT_map(NNVE_geom,*)
!    REAL*8  DCORVG_map(NNDIM_geom,*)
!    EXTERNAL MatInv3x3 ! Hypothetical 3x3 matrix inversion
!
!    IER_map = 0
!!
!!   1. Get nodal coordinates for the element IEL_TARGET
!    DO IVE=1,NNVE_geom
!        JP=KVERT_map(IVE,IEL_TARGET)
!        DX_geom(IVE)=DCORVG_map(1,JP)
!        DY_geom(IVE)=DCORVG_map(2,JP)
!        DZ_geom(IVE)=DCORVG_map(3,JP)
!    END DO
!!
!!   2. Calculate DJmn coefficients for geometry mapping (trilinear)
!!      These are based on DX_geom, DY_geom, DZ_geom (similar to GetForces)
!!      Q8_geom = 0.125D0
!    DJ11 = (DX_geom(1)+...+DX_geom(8))*Q8_geom
!    ! ... all other DJmn based on DX_geom, DY_geom, DZ_geom ...
!    DJ83 = (-DZ_geom(1)+...+ -DZ_geom(8))*Q8_geom 
!!
!!   3. Newton-Raphson Iteration
!    Xi_current = 0.0D0 ! Initial guess: center of reference element
!
!    DO ITER = 1, MAX_ITER_MAP
!!       a. Calculate current physical position X_calc based on Xi_current
!!          and the geometric mapping DJmn terms (like XX,YY,ZZ in GetForces)
!!          Geom_Jacobian here is J_map(xi_current)
!        Geom_Jacobian(1,1)=DJ21+DJ51*Xi_current(2)+DJ61*Xi_current(3)+DJ81*Xi_current(2)*Xi_current(3)
!        ! ... all components of Geom_Jacobian ...
!        Geom_Jacobian(3,3)=DJ43+DJ63*Xi_current(1)+DJ73*Xi_current(2)+DJ83*Xi_current(1)*Xi_current(2)
!
!        X_calc(1)=DJ11+Geom_Jacobian(1,1)*Xi_current(1)+DJ31*Xi_current(2)+DJ41*Xi_current(3)+DJ71*Xi_current(2)*Xi_current(3) 
!        X_calc(2)=DJ12+DJ22*Xi_current(1)+Geom_Jacobian(2,2)*Xi_current(2)+DJ42*Xi_current(3)+DJ62*Xi_current(1)*Xi_current(3)
!        X_calc(3)=DJ13+DJ23*Xi_current(1)+DJ33*Xi_current(2)+Geom_Jacobian(3,3)*Xi_current(3)+DJ53*Xi_current(1)*Xi_current(2)
!
!!       b. Calculate residual vector F_vec
!        F_vec(1) = X_calc(1) - X_phys(1)
!        F_vec(2) = X_calc(2) - X_phys(2)
!        F_vec(3) = X_calc(3) - X_phys(3)
!!
!!       c. Check for convergence
!        IF (MAXVAL(ABS(F_vec)) .LT. TOL_MAP) THEN
!            Xi_ref_out = Xi_current
!            RETURN ! Converged
!        END IF
!!
!!       d. Calculate update Delta_Xi = - Geom_Jacobian_inv * F_vec
!        CALL MatInv3x3(Geom_Jacobian, Geom_Jacobian_inv, DET_geom_jac) ! Invert Geom_Jacobian
!        IF (ABS(DET_geom_jac) .LT. 1.D-12) THEN
!            IER_map = -2 ! Singular Jacobian
!            RETURN
!        END IF
!        Delta_Xi(1) = -(Geom_Jacobian_inv(1,1)*F_vec(1) + Geom_Jacobian_inv(1,2)*F_vec(2) + Geom_Jacobian_inv(1,3)*F_vec(3))
!        Delta_Xi(2) = -(Geom_Jacobian_inv(2,1)*F_vec(1) + Geom_Jacobian_inv(2,2)*F_vec(2) + Geom_Jacobian_inv(2,3)*F_vec(3))
!        Delta_Xi(3) = -(Geom_Jacobian_inv(3,1)*F_vec(1) + Geom_Jacobian_inv(3,2)*F_vec(2) + Geom_Jacobian_inv(3,3)*F_vec(3))
!!
!!       e. Update Xi_current
!        Xi_current = Xi_current + Delta_Xi
!    END DO
!!
!    IER_map = -1 ! Failed to converge
!    Xi_ref_out = Xi_current ! Return last estimate
!END SUBROUTINE
!```
!
!**Summary of what's needed vs. what we have:**
!
!*   **We can definitely do it.**
!*   The main new piece of logic is the `PhysicalToRefCoords` subroutine. This routine itself reuses the geometric mapping coefficients (`DJmn`) that are already calculated in `GetForces` (though for the 8-node trilinear geometry, not the 27-node field).
!*   If such an inverse mapping utility already exists in your legacy code, then the `GetVelocityAtArbitraryPoint` becomes very similar to `GetVelocityAtElementCentroid`, just using the mapped reference coordinates instead of $(0,0,0)$.
!
!So, the answer is: **Yes, we can handle this task, but it requires either finding an existing inverse mapping routine or implementing one (like the Newton-Raphson sketch above).** The rest of the machinery is in place.






!!************************************************************************
!      SUBROUTINE PhysicalToRefCoords(IEL_TARGET, X_phys, KVERT_map, DCORVG_map, &
!                                    and the `PhysicalToRefCoords` by incorporating this. I'll make the Newton-Raphson a bit more completeXi_ref_out, IER_map)
!!
!! Purpose: Converts physical coordinates X_phys to and embed it, or at least detail it clearly.
!
!!          Xi_ref_out for the given element IEL_TARGET, assuming a
!!          trilinear its helper `Local_PhysicalToRefCoords`:**
!
!!************************************************************************
!      SUBROUTINE GetVelocityAtArbitraryPoint(IEL_TARGET, X_phys, &
!     &                                isoparametric mapping for geometry.
!!
!! Inputs:
!!   IEL_TARGET : Global index of the element
!     U1,U2,U3, KVERT,DCORVG,ELE, &
!     &                                    Vel!   X_phys(3)  : Physical coordinates (xp, yp, zp) of the point
!!AtPoint, IER_overall)
!!
!! Purpose: Calculates the velocity (U1, U2, U   KVERT_map  : Element vertex connectivity (NNVE_geom rows, NEL cols)
!!   DCORVG_3) at an arbitrary
!!          physical point X_phys known to be within element IEL_TARGET.
!!
!map : Global vertex coordinates (NNDIM_geom rows, NVT cols)
!!
!! Outputs:
!! IEL_TARGET: The global index of the element containing X_phys.
!! X_phys(3!   Xi_ref_out(3): Reference coordinates (xi1, xi2, xi3)
!!   IER): The physical coordinates (xp, yp, zp) of the point.
!! U1,U2,U3: Nodal velocity arrays.
!! KVERT, DCORVG: Mesh data.
!! ELE: Function pointer to_map      : Error/Status code for the mapping
!!                  0 = Success
!!                 -1 = Failed to converge within MAX_ITER_MAP
!!                 -2 = Singular or near-singular Jacobian encountered
!!                  element routine (e.g., E013 for Q2 velocity).
!! VelAtPoint(3): Output array for velocity (Ux, Uy, Uz) at X_phys.
!! IER_overall: Overall error indicator (0-3 = Converged point is outside reference element bounds
!!
!! Parameters for geometry mapping (trilinear, 8- for success).
!!************************************************************************
!      USE PP3D_MPI, ONLY:myid ! Fornode)
!      PARAMETER (NNVE_geom=8, NNDIM_geom=3)
!      PARAMETER debug
!      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL (Q8_geom=0.125D0) ! 1.0/8.0
!! Parameters(B)
!      CHARACTER SUB*6,FMT*15,CPARAM*120
! for Newton-Raphson
!      PARAMETER (MAX_ITER_MAP=20, TOL_MAP_RES=1.!
!      PARAMETER (NNBAS=27,NNDER=10,NNVE=8,NNDIM=D-8, TOL_MAP_DXI=1.D-7)
!
!      IMPLICIT DOUBLE PREC3) ! NNCUBP not needed
!      PARAMETER (Q8=0.125D0)
!!
!ISION (A-H,O-U,W-Z)
!
!      INTEGER IEL_TARGET, IER_map,      INTEGER IEL_TARGET, IER_overall, IER_map
!      REAL*8  X_phys(NNDIM)
!      REAL*8  U1(*),U2(*),U3(*)
!      REAL ITER, IVE, JP, I, J
!      REAL*8 X_phys(NNDIM_geom), Xi_ref_out(NNDIM_geom)
!      REAL*8 Xi_current(NNDIM_*8  DCORVG(NNDIM,*)
!      INTEGER KVERT(NNVE,*)
!      geom), F_vec(NNDIM_geom), Delta_Xi(NNDIM_geom)
!      REALINTEGER KDFG(NNBAS),KDFL(NNBAS)
!!
!      REAL*8  Vel*8 Geom_Jacobian(NNDIM_geom,NNDIM_geom)
!      REAL*8 X_calcAtPoint(NNDIM)
!      REAL*8  Xi_ref(NNDIM) ! Reference coordinates of(NNDIM_geom)
!      REAL*8 DX_elem_nodes(NNVE_geom), DY_ X_phys
!!
!      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
!      COMMON /ERRCTL/ IER,ICHECKelem_nodes(NNVE_geom), DZ_elem_nodes(NNVE_geom)
!      REAL*8 DJ11, DJ12, DJ13, DJ21, DJ22, DJ23 ! IER is used by ELE
!      COMMON /CHAR/   SUB,FMT(3),CPARAM
!      , DJ31, DJ32, DJ33, &
!               DJ41, DJ42,COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC( DJ43, DJ51, DJ52, DJ53, DJ61, DJ62,NNDIM,NNDIM),DETJ,&
!                      DBAS(NNDIM,NNBAS,NND DJ63, &
!               DJ71, DJ72, DJ73, DJ81, DJER),BDER(NNDER),KVE(NNVE),&
!                      IEL,NDIM_82, DJ83
!      REAL*8 DET_geom_jac
!
!      INTEGER KVERT_map(NNVE_geom,*)
!      REAL*8  DCORVG_map(NNDIM_geom,*)
!common ! Renamed NDIM in COMMON to avoid conflict
!      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE_common,NEE,NAE,NVEL,NEEL,NVED,&
!                      NVAR      LOGICAL BDER_geom(4) ! Not strictly needed here, but for consistency if using ELE for geom,NEAR,NBCT,NVBD,NEBD,NABD
!      COMMON /COAUX1/ KDFG,
!
!!     Intrinsic functions
!      INTRINSIC ABS, MAXVAL
!
!      IER_map = 0
!!
!!KDFL,IDFL
!!
!      SAVE
!C
!      IF (ICHECK.GE.99 --- 1. Get nodal coordinates for the target element's vertices ---
!      DO IVE=1,NNVE_geom7) CALL OTRC('GetVelArbPt','NEW')
!      IER_overall = 0
!      Vel
!        JP=KVERT_map(IVE,IEL_TARGET)
!        DX_elem_nodes(IVE)=DCAtPoint(1:NNDIM) = 0.0D0
!!
!! --- Set the currentORVG_map(1,JP)
!        DY_elem_nodes(IVE)=DCORVG_map element number for COMMON /ELEM/ ---
!      IEL = IEL_TARGET
!      IF (IEL_TARGET .(2,JP)
!        DZ_elem_nodes(IVE)=DCORVG_map(3,JPLT. 1 .OR. IEL_TARGET .GT. NEL) THEN
!        IF (myid .)
!      END DO
!!
!! --- 2. Calculate DJmn coefficients for the trilinear geometry mapping ---
!!EQ. 0) WRITE(*,*) 'GetVelArbPt: Invalid IEL_TARGET:', IEL_TARGET     (Copied and adapted from GetForces, using DX_elem_nodes etc.)
!      DJ11=( DX
!        IER_overall = -10
!        GOTO 99998
!      END IF
!_elem_nodes(1)+DX_elem_nodes(2)+DX_elem_nodes(3)+DX!
!! --- 1. Perform Inverse Mapping: Physical to Reference Coordinates ---
!      CALL Local_PhysicalToRefCoords(_elem_nodes(4)+DX_elem_nodes(5)+DX_elem_nodes(6)+DXIEL_TARGET, X_phys, KVERT, DCORVG, &
!     &                               Xi_ref, I_elem_nodes(7)+DX_elem_nodes(8))*Q8_geom
!      DJ12=( DY_ER_map)
!      IF (IER_map .NE. 0) THEN
!        IF (myidelem_nodes(1)+DY_elem_nodes(2)+DY_elem_nodes(3)+DY_ .EQ. 0) WRITE(*,*) 'GetVelArbPt: Inverse map failed for IEL',elem_nodes(4)+DY_elem_nodes(5)+DY_elem_nodes(6)+DY_ IEL_TARGET, &
!     &                              ' err:', IER_map
!        IER_overall = Ielem_nodes(7)+DY_elem_nodes(8))*Q8_geom
!      DJ13=(ER_map
!        GOTO 99998
!      END IF
!
!! --- Optional: Check if mapped DZ_elem_nodes(1)+DZ_elem_nodes(2)+DZ_elem_nodes(3)+DZ_elem_nodes(4)+DZ_elem_nodes(5)+DZ_elem_nodes(6)+ Xi_ref is within bounds [-1,1] (sanity check) ---
!!     DO I=1,NDZ_elem_nodes(7)+DZ_elem_nodes(8))*Q8_geom
!      DJ2NDIM
!!        IF (Xi_ref(I) .LT. -1.0001D0 .1=(-DX_elem_nodes(1)+DX_elem_nodes(2)+DX_elem_nodes(3)-OR. Xi_ref(I) .GT. 1.0001D0) THEN
!!DX_elem_nodes(4)-DX_elem_nodes(5)+DX_elem_nodes(6)+           IF (myid .EQ. 0) WRITE(*,*) 'GetVelArbPt: Xi_refDX_elem_nodes(7)-DX_elem_nodes(8))*Q8_geom
!      DJ2 out of bounds for IEL', IEL_TARGET, &
!!    &                                  ' Xi_ref:', Xi_ref2=(-DY_elem_nodes(1)+DY_elem_nodes(2)+DY_elem_nodes(
!!           IER_overall = -11
!!           GOTO 99998
!!        END IF3)-DY_elem_nodes(4)-DY_elem_nodes(5)+DY_elem_nodes(6)+DY_elem_nodes(7)-DY_elem_nodes(8))*Q8_geom
!      
!!     END DO
!!
!! --- 2. Initialize BDER for values only ---
!      DO I= DJ23=(-DZ_elem_nodes(1)+DZ_elem_nodes(2)+DZ_elem_1,NNDER
!        BDER(I)=.FALSE.
!      END DO
!      BDER(nodes(3)-DZ_elem_nodes(4)-DZ_elem_nodes(5)+DZ_elem_1)=.TRUE.
!!
!! --- 3. Setup Element Type (assuming Q2 for velocity like E0nodes(6)+DZ_elem_nodes(7)-DZ_elem_nodes(8))*Q8_geom
!      DJ31=(-DX_elem_nodes(1)-DX_elem_nodes(2)+DX_13) ---
!      IELTYP = -1
!      CALL ELE(0D0,0D0,elem_nodes(3)+DX_elem_nodes(4)-DX_elem_nodes(5)-DX_0D0,IELTYP) ! ELE is E013, IELTYP becomes 13
!      elem_nodes(6)+DX_elem_nodes(7)+DX_elem_nodes(8))*Q8IDFL=NDFL(IELTYP)
!      IF (IDFL .NE. NNBAS .AND._geom
!      DJ32=(-DY_elem_nodes(1)-DY_elem_nodes(2)+ IELTYP .EQ. 13) THEN
!         IF (myid .EQ. 0) WRITE(*DY_elem_nodes(3)+DY_elem_nodes(4)-DY_elem_nodes(5)-,*) 'GetVelArbPt: IDFL mismatch for E013'
!         IER_overall =DY_elem_nodes(6)+DY_elem_nodes(7)+DY_elem_nodes(8))* -12
!         GOTO 99998
!      END IF
!!
!! --- 4Q8_geom
!      DJ33=(-DZ_elem_nodes(1)-DZ_elem_nodes(. Get DOF mapping for the target element ---
!      CALL NDFGL(IEL_TARGET,1,IELTY2)+DZ_elem_nodes(3)+DZ_elem_nodes(4)-DZ_elem_nodes(P,KVERT,KEDGE,KAREA,KDFG,KDFL) ! KEDGE, KAREA5)-DZ_elem_nodes(6)+DZ_elem_nodes(7)+DZ_elem_nodes(8))*Q8_geom
!      DJ41=(-DX_elem_nodes(1)-DX_elem_ dummy
!      IF (IER.LT.0 .OR. IER_global.NE.0) THEN ! Inodes(2)-DX_elem_nodes(3)-DX_elem_nodes(4)+DX_elem_ER_global is COMMON /ERRCTL/ IER
!         IF (myid .EQ. 0) WRITE(*nodes(5)+DX_elem_nodes(6)+DX_elem_nodes(7)+DX_elem_nodes(8))*Q8_geom
!      DJ42=(-DY_elem_nodes(1)-DY_,*) 'GetVelArbPt: NDFGL error, IER=', IER_global
!         IER_elem_nodes(2)-DY_elem_nodes(3)-DY_elem_nodes(4)+DY_overall = -13
!         GOTO 99998
!      END IF
!      IER_global =elem_nodes(5)+DY_elem_nodes(6)+DY_elem_nodes(7)+DY_ 0 ! Reset global IER before calling ELE
!!
!! --- 5. Evaluate Basis Functions at the mappedelem_nodes(8))*Q8_geom
!      DJ43=(-DZ_elem_nodes(1)- reference point Xi_ref ---
!      CALL ELE(Xi_ref(1), Xi_ref(2), XiDZ_elem_nodes(2)-DZ_elem_nodes(3)-DZ_elem_nodes(4)+DZ_elem_nodes(5)+DZ_elem_nodes(6)+DZ_elem_nodes(7)+_ref(3), 0) ! IPAR = 0
!      IF (IER_global.LT.0) THEN
!         IF (myid .EQ. 0) WRITE(*,*) 'GetVelArbDZ_elem_nodes(8))*Q8_geom
!      DJ51=( DX_elem_nodes(1)-Pt: ELE error, IER=', IER_global
!         IER_overall = -14
!         GDX_elem_nodes(2)+DX_elem_nodes(3)-DX_elem_nodes(4)+OTO 99998
!      END IF
!!
!! --- 6. Interpolate Velocity Components ---
!      DODX_elem_nodes(5)-DX_elem_nodes(6)+DX_elem_nodes(7)-DX_elem_nodes(8))*Q8_geom
!      DJ52=( DY_elem_nodes( I=1,IDFL
!        IG = KDFG(I)
!        ! DBAS(1,1)-DY_elem_nodes(2)+DY_elem_nodes(3)-DY_elem_nodes(4)+KDFL(I),1) contains value of basis func KDFL(I) at Xi_ref
!        DY_elem_nodes(5)-DY_elem_nodes(6)+DY_elem_nodes(7)-VelAtPoint(1) = VelAtPoint(1) + U1(IG) * DBAS(1,KDY_elem_nodes(8))*Q8_geom
!      DJ53=( DZ_elem_nodes(1)-DFL(I),1)
!        VelAtPoint(2) = VelAtPoint(2) + UDZ_elem_nodes(2)+DZ_elem_nodes(3)-DZ_elem_nodes(4)+2(IG) * DBAS(1,KDFL(I),1)
!        VelAtPoint(DZ_elem_nodes(5)-DZ_elem_nodes(6)+DZ_elem_nodes(7)-3) = VelAtPoint(3) + U3(IG) * DBAS(1,KDFLDZ_elem_nodes(8))*Q8_geom
!      DJ61=( DX_elem_nodes(1)-DX(I),1)
!      END DO
!!
!      GOTO 99999 ! Normal exit_elem_nodes(2)-DX_elem_nodes(3)+DX_elem_nodes(4)-DX_elem_
!!
!99998 CONTINUE ! Error exit
!      VelAtPoint(1:NNDIM) = -1nodes(5)+DX_elem_nodes(6)+DX_elem_nodes(7)-DX_elem_.0D+30 ! Indicate error
!C
!99999 CONTINUE
!      END SUBROUTINE GetVelocitynodes(8))*Q8_geom
!      DJ62=( DY_elem_nodes(1)-DY_elem_nodesAtArbitraryPoint
!!************************************************************************
!
!!************************************************************************
!      SUBROUTINE Local(2)-DY_elem_nodes(3)+DY_elem_nodes(4)-DY_elem_nodes(5)+_PhysicalToRefCoords(IEL_map, X_phys_map, &
!     &                             KVERT_lDY_elem_nodes(6)+DY_elem_nodes(7)-DY_elem_nodes(8))*, DCORVG_l, Xi_ref_map, IER_map_l)
!!
!! PurposeQ8_geom
!      DJ63=( DZ_elem_nodes(1)-DZ_elem_nodes(2)-DZ: Finds reference coordinates Xi_ref_map for a physical point
!!          X_phys_map known to be within element_elem_nodes(3)+DZ_elem_nodes(4)-DZ_elem_nodes(5)+DZ_elem_ IEL_map.
!!          Uses Newton-Raphson for an 8-node trilinear geometry mappingnodes(6)+DZ_elem_nodes(7)-DZ_elem_nodes(8))*Q8_geom.
!!
!      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z
!      DJ71=( DX_elem_nodes(1)+DX_elem_nodes(2)-DX_elem_)
!      PARAMETER (NNVE_geom=8, NNDIM_geom=3, MAX_ITER_NRnodes(3)-DX_elem_nodes(4)-DX_elem_nodes(5)-DX_elem_=25, TOL_NR=1.D-9)
!      PARAMETER (Q8_geom=0nodes(6)+DX_elem_nodes(7)+DX_elem_nodes(8))*Q8_geom.125D0)
!!
!      INTEGER IEL_map, IER_map_l,
!      DJ72=( DY_elem_nodes(1)+DY_elem_nodes(2)-DY_elem_nodes(3)-DY_elem_nodes(4)-DY_elem_nodes(5)-DY_ ITER_NR, IVE, JP
!      REAL*8  X_phys_map(NNDIM_geom),elem_nodes(6)+DY_elem_nodes(7)+DY_elem_nodes(8))*Q8_geom
!      DJ73=( DZ_elem_nodes(1)+DZ_elem_nodes(2)- Xi_ref_map(NNDIM_geom)
!      REAL*8  Xi_curr(NNDDZ_elem_nodes(3)-DZ_elem_nodes(4)-DZ_elem_nodes(5)-IM_geom), F_resid(NNDIM_geom), dXi(NNDIM_geom)
!      DZ_elem_nodes(6)+DZ_elem_nodes(7)+DZ_elem_nodes(8))*REAL*8  GeomJ(NNDIM_geom,NNDIM_geom), GeomJ_invQ8_geom
!      DJ81=(-DX_elem_nodes(1)+DX_elem_nodes((NNDIM_geom,NNDIM_geom)
!      REAL*8  X_calc_curr2)-DX_elem_nodes(3)+DX_elem_nodes(4)+DX_elem_nodes(5)-DX_elem_nodes(6)+DX_elem_nodes(7)-DX_elem_nodes(8))*(NNDIM_geom)
!      REAL*8  Nodes_X(NNVE_geom), Nodes_Y(NNVE_geom), Nodes_Z(NNVE_geom)
!      REAL*8  DJ11,Q8_geom
!      DJ82=(-DY_elem_nodes(1)+DY_elem_nodes( DJ12, DJ13, DJ21, DJ22, DJ23, DJ31,2)-DY_elem_nodes(3)+DY_elem_nodes(4)+DY_elem_nodes(5)-DY_elem_nodes(6)+DY_elem_nodes(7)-DY_elem_nodes(8))* DJ32, DJ33
!      REAL*8  DJ41, DJ42, DJ43, DJ51, DJ52, DJ53, DJ61, DJ62, DJ63Q8_geom
!      DJ83=(-DZ_elem_nodes(1)+DZ_elem_nodes(
!      REAL*8  DJ71, DJ72, DJ73, DJ81, DJ82)-DZ_elem_nodes(3)+DZ_elem_nodes(4)+DZ_elem_nodes(5)-DZ_elem_nodes(6)+DZ_elem_nodes(7)-DZ_elem_nodes(2, DJ83
!      REAL*8  DET_GeomJ
!!
!      INTEGER KVERT_8))*Q8_geom
!!
!! --- 3. Newton-Raphson Iteration ---
!      Xi_currentl(NNVE_geom,*)
!      REAL*8  DCORVG_l(NNDIM_(1) = 0.0D0 ! Initial guess: center of reference element [-1,1]^3
!      geom,*)
!!
!      IER_map_l = 0
!!
!! --- 1. Get nodal coordinates ofXi_current(2) = 0.0D0
!      Xi_current(3) = 0 the element's vertices ---
!      DO IVE=1,NNVE_geom
!        JP=KVERT_l.0D0
!
!      DO ITER = 1, MAX_ITER_MAP
!!       a. Calculate current(IVE,IEL_map)
!        Nodes_X(IVE)=DCORVG_l(1,JP)
!        Nodes_Y(IVE)=DCORVG_l(2,JP)
!        Nodes_Z(IVE physical position X_calc(Xi_current)
!!          and the Jacobian of the geometric map J_map(Xi_)=DCORVG_l(3,JP)
!      END DO
!!
!! --- 2. Calculate DJmn coefficients for trilinear geometry mapping ---
!      DJ11=(Nodes_X(1)+Nodes_X(2)+current)
!!          (This is the same J_map = DJAC as in GetForces, evaluated at Xi_Nodes_X(3)+Nodes_X(4)+Nodes_X(5)+Nodes_X(6)+current)
!        Geom_Jacobian(1,1)=DJ21+DJ51*Xi_current(2Nodes_X(7)+Nodes_X(8))*Q8_geom
!      DJ12=(Nodes_)+DJ61*Xi_current(3)+DJ81*Xi_current(2)*Xi_currentY(1)+Nodes_Y(2)+Nodes_Y(3)+Nodes_Y(4)+Nodes_(3)
!        Geom_Jacobian(1,2)=DJ31+DJ51*XiY(5)+Nodes_Y(6)+Nodes_Y(7)+Nodes_Y(8))*Q8_current(1)+DJ71*Xi_current(3)+DJ81*Xi_current(1)*Xi__geom
!      DJ13=(Nodes_Z(1)+Nodes_Z(2)+Nodes_Z(current(3)
!        Geom_Jacobian(1,3)=DJ41+DJ61*3)+Nodes_Z(4)+Nodes_Z(5)+Nodes_Z(6)+Nodes_Z(Xi_current(1)+DJ71*Xi_current(2)+DJ81*Xi_current(7)+Nodes_Z(8))*Q8_geom
!      DJ21=(-Nodes_X(1)+Nodes_1)*Xi_current(2)
!        Geom_Jacobian(2,1)=DJ22+X(2)+Nodes_X(3)-Nodes_X(4)-Nodes_X(5)+Nodes_DJ52*Xi_current(2)+DJ62*Xi_current(3)+DJ82*Xi_current(2)*Xi_current(3)
!        Geom_Jacobian(2,2)=X(6)+Nodes_X(7)-Nodes_X(8))*Q8_geom
!      DJ22=(-DJ32+DJ52*Xi_current(1)+DJ72*Xi_current(3)+Nodes_Y(1)+Nodes_Y(2)+Nodes_Y(3)-Nodes_Y(4)-Nodes_DJ82*Xi_current(1)*Xi_current(3)
!        Geom_Jacobian(Y(5)+Nodes_Y(6)+Nodes_Y(7)-Nodes_Y(8))*Q8_geom
!      DJ23=(-Nodes_Z(1)+Nodes_Z(2)+Nodes_Z(2,3)=DJ42+DJ62*Xi_current(1)+DJ72*Xi_3)-Nodes_Z(4)-Nodes_Z(5)+Nodes_Z(6)+Nodes_Z(current(2)+DJ82*Xi_current(1)*Xi_current(2)
!        Geom7)-Nodes_Z(8))*Q8_geom
!      DJ31=(-Nodes_X(1)-_Jacobian(3,1)=DJ23+DJ53*Xi_current(2)+DJ6Nodes_X(2)+Nodes_X(3)+Nodes_X(4)-Nodes_X(5)-3*Xi_current(3)+DJ83*Xi_current(2)*Xi_current(3)Nodes_X(6)+Nodes_X(7)+Nodes_X(8))*Q8_geom
!      
!        Geom_Jacobian(3,2)=DJ33+DJ53*Xi_current(DJ32=(-Nodes_Y(1)-Nodes_Y(2)+Nodes_Y(3)+Nodes_1)+DJ73*Xi_current(3)+DJ83*Xi_current(1)*Xi_current(3)
!        Geom_Jacobian(3,3)=DJ43+DJ63*Y(4)-Nodes_Y(5)-Nodes_Y(6)+Nodes_Y(7)+Nodes_Xi_current(1)+DJ73*Xi_current(2)+DJ83*Xi_current(Y(8))*Q8_geom
!      DJ33=(-Nodes_Z(1)-Nodes_Z(2)+Nodes_Z(3)+Nodes_Z(4)-Nodes_Z(5)-Nodes_Z(1)*Xi_current(2)
!
!        X_calc(1)=DJ11+Geom_Jacob6)+Nodes_Z(7)+Nodes_Z(8))*Q8_geom
!      DJ41=(-ian(1,1)*Xi_current(1)+Geom_Jacobian(1,2)*Xi_currentNodes_X(1)-Nodes_X(2)-Nodes_X(3)-Nodes_X(4)+(2)+Geom_Jacobian(1,3)*Xi_current(3)
!        X_calcNodes_X(5)+Nodes_X(6)+Nodes_X(7)+Nodes_X(8))*(2)=DJ12+Geom_Jacobian(2,1)*Xi_current(1)+GeQ8_geom
!      DJ42=(-Nodes_Y(1)-Nodes_Y(2)-Nodes_om_Jacobian(2,2)*Xi_current(2)+Geom_Jacobian(2,3Y(3)-Nodes_Y(4)+Nodes_Y(5)+Nodes_Y(6)+Nodes_)*Xi_current(3)
!        X_calc(3)=DJ13+Geom_JacobianY(7)+Nodes_Y(8))*Q8_geom
!      DJ43=(-Nodes_Z((3,1)*Xi_current(1)+Geom_Jacobian(3,2)*Xi_current1)-Nodes_Z(2)-Nodes_Z(3)-Nodes_Z(4)+Nodes_Z((2)+Geom_Jacobian(3,3)*Xi_current(3)
!!
!!       5)+Nodes_Z(6)+Nodes_Z(7)+Nodes_Z(8))*Q8_geomb. Calculate residual vector F_vec = X_calc - X_phys
!        F_vec(1) = X_calc(1) - X_phys(1)
!        F_vec(2) = X_calc
!      DJ51=(Nodes_X(1)-Nodes_X(2)+Nodes_X(3)-Nodes_X(4)+Nodes_X(5)-Nodes_X(6)+Nodes_X(7)-(2) - X_phys(2)
!        F_vec(3) = X_calc(3Nodes_X(8))*Q8_geom
!      DJ52=(Nodes_Y(1)-Nodes_) - X_phys(3)
!!
!!       c. Check for convergence based on residual
!        IF (MAXY(2)+Nodes_Y(3)-Nodes_Y(4)+Nodes_Y(5)-Nodes_VAL(ABS(F_vec)) .LT. TOL_MAP_RES) THEN
!            Xi_ref_Y(6)+Nodes_Y(7)-Nodes_Y(8))*Q8_geom
!      DJ5out = Xi_current
!            GOTO 100 ! Converged based on residual
!        END IF
!!3=(Nodes_Z(1)-Nodes_Z(2)+Nodes_Z(3)-Nodes_Z(4)+
!!       d. Solve J_map * Delta_Xi = -F_vec for Delta_Xi
!!          (Nodes_Z(5)-Nodes_Z(6)+Nodes_Z(7)-Nodes_Z(8))*Q8_geom
!      DJ61=(Nodes_X(1)-Nodes_X(2)-Nodes_Requires inverting Geom_Jacobian or solving linear system)
!        DET_geom_jac = Geom_Jacobian(X(3)+Nodes_X(4)-Nodes_X(5)+Nodes_X(6)+Nodes_X(1,1)*(Geom_Jacobian(2,2)*Geom_Jacobian(3,3)-Geom_7)-Nodes_X(8))*Q8_geom
!      DJ62=(Nodes_Y(1)-Nodes_YJacobian(3,2)*Geom_Jacobian(2,3)) &
!                     - Geom_Jacobian((2)-Nodes_Y(3)+Nodes_Y(4)-Nodes_Y(5)+Nodes_Y(1,2)*(Geom_Jacobian(2,1)*Geom_Jacobian(3,3)-6)+Nodes_Y(7)-Nodes_Y(8))*Q8_geom
!      DJ63=(Geom_Jacobian(3,1)*Geom_Jacobian(2,3)) &
!                     + Geom_JacobNodes_Z(1)-Nodes_Z(2)-Nodes_Z(3)+Nodes_Z(4)-Nodes_ian(1,3)*(Geom_Jacobian(2,1)*Geom_Jacobian(3,Z(5)+Nodes_Z(6)+Nodes_Z(7)-Nodes_Z(8))*Q82)-Geom_Jacobian(3,1)*Geom_Jacobian(2,2))
!
!        IF (ABS_geom
!      DJ71=(Nodes_X(1)+Nodes_X(2)-Nodes_X((DET_geom_jac) .LT. 1.D-12) THEN
!            IER_map = -3)-Nodes_X(4)-Nodes_X(5)-Nodes_X(6)+Nodes_X(7)+Nodes_X(8))*Q8_geom
!      DJ72=(Nodes_Y(1)+2 ! Singular Jacobian
!            Xi_ref_out = Xi_current ! Return last estimate
!            RETURN
!        END IFNodes_Y(2)-Nodes_Y(3)-Nodes_Y(4)-Nodes_Y(5)-
!        
!        INV_DET = 1.0D0 / DET_geom_jac
!        ! ExplicitNodes_Y(6)+Nodes_Y(7)-Nodes_Y(8))*Q8_geom
!       inverse for 3x3 matrix (adjugate / determinant)
!        Delta_Xi(1) = -DJ73=(Nodes_Z(1)+Nodes_Z(2)-Nodes_Z(3)-Nodes_Z(INV_DET * ( &
!             (Geom_Jacobian(2,2)*Geom_Jacobian(34)-Nodes_Z(5)-Nodes_Z(6)+Nodes_Z(7)+Nodes_Z(,3) - Geom_Jacobian(2,3)*Geom_Jacobian(3,2))*F_vec(8))*Q8_geom
!      DJ81=(-Nodes_X(1)+Nodes_X(2)-Nodes_1) + &
!             (Geom_Jacobian(1,3)*Geom_Jacobian(3X(3)+Nodes_X(4)+Nodes_X(5)-Nodes_X(6)+Nodes_X(,2) - Geom_Jacobian(1,2)*Geom_Jacobian(3,3))*F_vec(2) + &
!             (Geom_Jacobian(1,2)*Geom_Jacob7)-Nodes_X(8))*Q8_geom
!      DJ82=(-Nodes_Y(1)+Nodes_Y(2)-Nodes_Y(3)+Nodes_Y(4)+Nodes_Y(5)-Nodes_ian(2,3) - Geom_Jacobian(1,3)*Geom_Jacobian(2,Y(6)+Nodes_Y(7)-Nodes_Y(8))*Q8_geom
!      DJ82))*F_vec(3) )
!        Delta_Xi(2) = -INV_DET * (3=(-Nodes_Z(1)+Nodes_Z(2)-Nodes_Z(3)+Nodes_Z( &
!             (Geom_Jacobian(2,3)*Geom_Jacobian(3,1)4)+Nodes_Z(5)-Nodes_Z(6)+Nodes_Z(7)-Nodes_Z( - Geom_Jacobian(2,1)*Geom_Jacobian(3,3))*F_vec(8))*Q8_geom
!!
!! --- 3. Newton-Raphson Iteration ---
!      Xi_curr1) + &
!             (Geom_Jacobian(1,1)*Geom_Jacobian(3 = 0.0D0 ! Initial guess: center of reference element [-1,1]^3
!!
!      DO,3) - Geom_Jacobian(1,3)*Geom_Jacobian(3,1))*F ITER_NR = 1, MAX_ITER_NR
!!       a. Calculate current physical position X_calc_vec(2) + &
!             (Geom_Jacobian(1,3)*Geom_Jacobian(2,1) - Geom_Jacobian(1,1)*Geom_Jacobian(2,_curr using Xi_curr
!!          and the geometric mapping (forward mapping)
!        GeomJ(1,13))*F_vec(3) )
!        Delta_Xi(3) = -INV_DET * ()=DJ21+DJ51*Xi_curr(2)+DJ61*Xi_curr(3)+DJ8 &
!             (Geom_Jacobian(2,1)*Geom_Jacobian(3,2)1*Xi_curr(2)*Xi_curr(3)
!        GeomJ(1,2)= - Geom_Jacobian(2,2)*Geom_Jacobian(3,1))*F_vec(DJ31+DJ51*Xi_curr(1)+DJ71*Xi_curr(3)+1) + &
!             (Geom_Jacobian(1,2)*Geom_Jacobian(3DJ81*Xi_curr(1)*Xi_curr(3)
!        GeomJ(1,,1) - Geom_Jacobian(1,1)*Geom_Jacobian(3,2))*F3)=DJ41+DJ61*Xi_curr(1)+DJ71*Xi_curr(2)+DJ81*Xi_curr(1)*Xi_curr(2)
!        GeomJ(_vec(2) + &
!             (Geom_Jacobian(1,1)*Geom_Jacob2,1)=DJ22+DJ52*Xi_curr(2)+DJ62*Xi_ian(2,2) - Geom_Jacobian(1,2)*Geom_Jacobian(2,curr(3)+DJ82*Xi_curr(2)*Xi_curr(3)
!        Geom1))*F_vec(3) )
!!
!!       e. Update Xi_current
!        Xi_current = XiJ(2,2)=DJ32+DJ52*Xi_curr(1)+DJ72*_current + Delta_Xi
!!
!!       f. Check for convergence based on step size
!        IF (MAXVAL(ABS(Delta_Xi)) .LT. TOL_MAP_DXI) THEN
!             Xi_ref_out = Xi_Xi_curr(3)+DJ82*Xi_curr(1)*Xi_curr(3)
!        GeomJ(2,3)=DJ42+DJ62*Xi_curr(1)+DJ7current
!             GOTO 100 ! Converged based on step size
!        END IF
!      END DO
!!
!2*Xi_curr(2)+DJ82*Xi_curr(1)*Xi_curr(2)      IER_map = -1 ! Failed to converge within MAX_ITER_MAP
!      Xi_ref_out = Xi
!        GeomJ(3,1)=DJ23+DJ53*Xi_curr(2)+DJ63*Xi_curr(3)+DJ83*Xi_curr(2)*Xi_curr(_current ! Return last estimate
!      RETURN
!
!100   CONTINUE ! Converged
!!     Optional: Check if Xi3)
!        GeomJ(3,2)=DJ33+DJ53*Xi_curr(_ref_out is within reference bounds [-1,1]
!      DO I = 1, NNDIM_geom1)+DJ73*Xi_curr(3)+DJ83*Xi_curr(1)*Xi_
!          IF (Xi_ref_out(I) < -1.0D0 - TOL_MAP_REScurr(3)
!        GeomJ(3,3)=DJ43+DJ63*Xi_ .OR. Xi_ref_out(I) > 1.0D0 + TOL_MAP_RES)curr(1)+DJ73*Xi_curr(2)+DJ83*Xi_curr(1)* THEN
!              IER_map = -3 ! Converged point outside expected bounds
!              ! This might happen if XXi_curr(2)
!
!        X_calc_curr(1)=DJ11+GeomJ(1,1)*Xi_curr(1)+GeomJ(1,2)*Xi_curr(2)+Ge_phys was actually not in IEL_TARGET
!              EXIT
!          END IF
!      END DO
!
!END SUBROUTINE PhysicalToRefCoords

!************************************************************************
      SUBROUTINE EvaluateVelocityAtPoint(TargetPointPhys, Target_IEL, &
     &                                 U1_global, U2_global, U3_global, &
     &                                 KVERT, DCORVG, ELE, &
     &                                 VelocityAtPoint, RefCoords, Found)
! INPUTS:
!   TargetPointPhys(3)  : R*8, Desired physical point (xp, yp, zp)
!   Target_IEL          : I*4, The element IEL we assume the point is in
!   U1_global, U2_global, U3_global: R*8, Global nodal velocity arrays
!   KVERT               : I*4, Element vertex connectivity
!   DCORVG              : R*8, Global vertex coordinates
!   ELE                 : External subroutine (E013)
! OUTPUTS:
!   VelocityAtPoint(3)  : R*8, Interpolated velocity (ux, uy, uz) at TargetPointPhys
!   RefCoords(3)        : R*8, Calculated reference coords (xi1, xi2, xi3)
!   Found               : LOGICAL, .TRUE. if inverse mapping converged within bounds
!************************************************************************
      USE PP3D_MPI, ONLY:myid ! Maybe for error messages
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
!
      PARAMETER (NNBAS=27,NNDER=10,NNVE=8,NNDIM=3) ! NNBAS for Q2 velocity
      PARAMETER (Q8=0.125D0)
!
      REAL*8  TargetPointPhys(NNDIM)
      REAL*8  U1_global(*), U2_global(*), U3_global(*)
      REAL*8  DCORVG(NNDIM,*)
      INTEGER KVERT(NNVE,*)
      INTEGER KDFG_vel(NNBAS), KDFL_vel(NNBAS) ! For velocity DOFs
!
      REAL*8  VelocityAtPoint(NNDIM), RefCoords(NNDIM)
      LOGICAL Found
!
! Local variables for Newton-Raphson
      REAL*8  CurrentRefCoords(NNDIM), PrevRefCoords(NNDIM)
      REAL*8  CurrentPhysCoords(NNDIM)
      REAL*8  Residual(NNDIM)
      REAL*8  JacobianInvMap(NNDIM,NNDIM) ! This is DJAC from GetForces
      REAL*8  DetJacobianInvMap
      REAL*8  DeltaRefCoords(NNDIM)
      INTEGER Iter, MaxIter
      REAL*8  Tolerance, ErrorNorm
!
! Common blocks (copied for ELE, etc.)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(NNDIM,NNDIM),DETJ,&
                      DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),&
                      IEL_common,NDIM_common ! Renamed to avoid conflict with Target_IEL
      COMMON /COAUX1/ KDFG_vel, KDFL_vel, IDFL_vel ! Renamed for clarity

      SAVE
C
      IF (ICHECK.GE.997) CALL OTRC('EvalVelAtPt','NEW')
      IER = 0
      Found = .FALSE.
      VelocityAtPoint = 0.0D0
      RefCoords = -999.0D0 ! Initialize to an unlikely value

! --- Parameters for Newton-Raphson ---
      MaxIter = 20
      Tolerance = 1.0D-8
!
! --- Set IEL_common for common block /ELEM/ ---
      IEL_common = Target_IEL
!
! --- Initialize BDER for Q2 basis function values (no derivatives needed for interpolation) ---
!     However, E013 might require derivative flags for its internal workings,
!     and the Jacobian calculation (for Newton-Raphson) *does* need derivatives.
      DO I= 1,NNDER
        BDER(I)=.FALSE.
      END DO
      DO I=1,4 ! BDER(1) for value, BDER(2,3,4) for d/dx, d/dy, d/dz needed by Jacobian_inv_map
        BDER(I)=.TRUE.
      END DO
!
! --- Setup Element Type for Q2 Velocity (E013) ---
      IELTYP_vel = -1
      CALL ELE(0D0,0D0,0D0,IELTYP_vel) ! ELE is E013, IELTYP_vel becomes 13
      IDFL_vel = NDFL(IELTYP_vel)       ! Should be NNBAS = 27
      IF (IDFL_vel .NE. NNBAS .AND. IELTYP_vel .EQ. 13) THEN
          IF (myid .EQ. 0) WRITE(*,*) 'EvalVel: IDFL_vel mismatch E013'
          RETURN
      END IF
!
! --- Get element vertex coordinates for Target_IEL ---
      DO IVE=1,NNVE
        JP=KVERT(IVE,Target_IEL)
        KVE(IVE)=JP ! For common block /ELEM/ if E013 uses it directly
        DX(IVE)=DCORVG(1,JP)
        DY(IVE)=DCORVG(2,JP)
        DZ(IVE)=DCORVG(3,JP)
      END DO
!
! --- Pre-calculate DJmn coefficients for geometry mapping (trilinear) ---
!     These are used to calculate the forward map x(xi) and its Jacobian J_xi
      DJ11=( DX(1)+DX(2)+DX(3)+DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
      DJ12=( DY(1)+DY(2)+DY(3)+DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
!     ... (all DJmn terms as in GetForces, up to DJ83) ...
      DJ83=(-DZ(1)+DZ(2)-DZ(3)+DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8

! === Step 1: Inverse Isoparametric Mapping (Newton-Raphson) ===
      CurrentRefCoords = 0.0D0 ! Initial guess: center of reference element [-1,1]^3

      DO Iter = 1, MaxIter
        PrevRefCoords = CurrentRefCoords
        XI1_iter = CurrentRefCoords(1)
        XI2_iter = CurrentRefCoords(2)
        XI3_iter = CurrentRefCoords(3)

        ! Calculate Jacobian of forward map (J_xi) at CurrentRefCoords
        ! This is the DJAC(3,3) matrix calculation from GetForces
        JacobianInvMap(1,1)=DJ21+DJ51*XI2_iter+DJ61*XI3_iter+DJ81*XI2_iter*XI3_iter
        JacobianInvMap(1,2)=DJ31+DJ51*XI1_iter+DJ71*XI3_iter+DJ81*XI1_iter*XI3_iter
        ! ... (all 9 components of JacobianInvMap = J_xi) ...
        JacobianInvMap(3,3)=DJ43+DJ63*XI1_iter+DJ73*XI2_iter+DJ83*XI1_iter*XI3_iter

        DetJacobianInvMap = JacobianInvMap(1,1)*(JacobianInvMap(2,2)*JacobianInvMap(3,3)-JacobianInvMap(3,2)*JacobianInvMap(2,3))&
                         -JacobianInvMap(2,1)*(JacobianInvMap(1,2)*JacobianInvMap(3,3)-JacobianInvMap(3,2)*JacobianInvMap(1,3))&
                         +JacobianInvMap(3,1)*(JacobianInvMap(1,2)*JacobianInvMap(2,3)-JacobianInvMap(2,2)*JacobianInvMap(1,3))

        IF (ABS(DetJacobianInvMap) < 1.0D-12) THEN
            IF (myid .EQ. 0) WRITE(*,*) 'EvalVel: Singular Jacobian in Newton Iter', Iter
            Found = .FALSE.
            RETURN
        END IF

        ! Calculate current physical coords x(xi_iter), y(xi_iter), z(xi_iter)
        ! This is the XX, YY, ZZ calculation from GetForces
        CurrentPhysCoords(1)=DJ11+JacobianInvMap(1,1)*XI1_iter+DJ31*XI2_iter+DJ41*XI3_iter+DJ71*XI2_iter*XI3_iter ! Check this carefully for correct DJmn terms
        CurrentPhysCoords(2)=DJ12+DJ22*XI1_iter+JacobianInvMap(2,2)*XI2_iter+DJ42*XI3_iter+DJ62*XI1_iter*XI3_iter ! Check this carefully
        CurrentPhysCoords(3)=DJ13+DJ23*XI1_iter+DJ33*XI2_iter+JacobianInvMap(3,3)*XI3_iter+DJ53*XI1_iter*XI2_iter ! Check this carefully

        ! Calculate residual f = x(xi_iter) - x_p
        Residual(1) = CurrentPhysCoords(1) - TargetPointPhys(1)
        Residual(2) = CurrentPhysCoords(2) - TargetPointPhys(2)
        Residual(3) = CurrentPhysCoords(3) - TargetPointPhys(3)

        ! Check for convergence of residual
        ErrorNorm = SQRT(Residual(1)**2 + Residual(2)**2 + Residual(3)**2)
        IF (ErrorNorm < Tolerance) THEN
          Found = .TRUE.
          EXIT ! Newton iteration converged
        END IF

        ! Solve J_xi * DeltaRefCoords = -Residual for DeltaRefCoords
        ! (i.e., DeltaRefCoords = - (J_xi)^-1 * Residual)
        ! Need to compute inverse of JacobianInvMap (J_xi)
        ! For 3x3, can do it explicitly (adjugate / determinant)
        DetInv = 1.0D0 / DetJacobianInvMap
        DeltaRefCoords(1) = -DetInv * ( &
            (JacobianInvMap(2,2)*JacobianInvMap(3,3) - JacobianInvMap(2,3)*JacobianInvMap(3,2)) * Residual(1) + &
            (JacobianInvMap(1,3)*JacobianInvMap(3,2) - JacobianInvMap(1,2)*JacobianInvMap(3,3)) * Residual(2) + &
            (JacobianInvMap(1,2)*JacobianInvMap(2,3) - JacobianInvMap(1,3)*JacobianInvMap(2,2)) * Residual(3)  )
        DeltaRefCoords(2) = -DetInv * ( &
            (JacobianInvMap(2,3)*JacobianInvMap(3,1) - JacobianInvMap(2,1)*JacobianInvMap(3,3)) * Residual(1) + &
            (JacobianInvMap(1,1)*JacobianInvMap(3,3) - JacobianInvMap(1,3)*JacobianInvMap(3,1)) * Residual(2) + &
            (JacobianInvMap(1,3)*JacobianInvMap(2,1) - JacobianInvMap(1,1)*JacobianInvMap(2,3)) * Residual(3)  )
        DeltaRefCoords(3) = -DetInv * ( &
            (JacobianInvMap(2,1)*JacobianInvMap(3,2) - JacobianInvMap(2,2)*JacobianInvMap(3,1)) * Residual(1) + &
            (JacobianInvMap(1,2)*JacobianInvMap(3,1) - JacobianInvMap(1,1)*JacobianInvMap(3,2)) * Residual(2) + &
            (JacobianInvMap(1,1)*JacobianInvMap(2,2) - JacobianInvMap(1,2)*JacobianInvMap(2,1)) * Residual(3)  )

        CurrentRefCoords = CurrentRefCoords + DeltaRefCoords

        ! Optional: Check if CurrentRefCoords are still within reference element bounds [-1,1]
        ! If it goes too far out, the point might not be in this element.
        IF (ANY(ABS(CurrentRefCoords) > 1.0D0 + 10.0D0*Tolerance)) THEN ! Allow small overshoot
            ! IF (myid .EQ. 0) WRITE(*,*) 'EvalVel: Ref coords out of bounds', CurrentRefCoords
            ! Found = .FALSE. ! Could stop here or let it try to converge back
        END IF
      END DO ! Newton Iteration loop

      IF (.NOT. Found) THEN
        IF (myid .EQ. 0) WRITE(*,*) 'EvalVel: Newton did not converge. ErrorNorm:', ErrorNorm
        RETURN
      END IF

! === Step 2: Interpolate Velocity using Q2 Basis Functions ===
      RefCoords = CurrentRefCoords ! Store the converged reference coordinates

      ! Get global DOFs for velocity for Target_IEL
      CALL NDFGL(Target_IEL,1,IELTYP_vel,KVERT,KEDGE,KAREA,KDFG_vel,KDFL_vel) ! KEDGE, KAREA dummy if not used by NDFGL for Q2
      IF (IER.LT.0) THEN
          IF (myid .EQ. 0) WRITE(*,*) 'EvalVel: NDFGL error'
          Found = .FALSE.
          RETURN
      END IF

      ! Set BDER to calculate only values for interpolation
      DO I= 1,NNDER; BDER(I)=.FALSE.; END DO
      BDER(1)=.TRUE. ! Only need basis function values

      ! Evaluate Q2 basis functions at found RefCoords
      ! E013 with IPAR=0 uses XI1,XI2,XI3 directly, stores values in DBAS(1,:,1)
      ! It does not use the pre-computation DHELP for IPAR=0.
      ! DJAC and DETJ in COMMON /ELEM/ are needed if E013 also transforms derivatives
      ! even for IPAR=0, but we only set BDER(1)=.TRUE.
      ! For safety, ensure DJAC and DETJ are set for the *current physical element*
      ! although E013 for values only might not use them.
      XI1_ref = RefCoords(1)
      XI2_ref = RefCoords(2)
      XI3_ref = RefCoords(3)
      
      ! Recalculate DJAC and DETJ for the element at RefCoords just in case E013 uses them
      ! even when only values are requested (IPAR=0, BDER(1)=.TRUE.).
      ! This is the same DJAC used in Newton iteration, J_xi.
      DJAC(1,1)=DJ21+DJ51*XI2_ref+DJ61*XI3_ref+DJ81*XI2_ref*XI3_ref
      ! ... (all 9 components of DJAC) ...
      DJAC(3,3)=DJ43+DJ63*XI1_ref+DJ73*XI2_ref+DJ83*XI1_ref*XI3_ref
      DETJ = DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))&
                -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))&
                +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))

      CALL ELE(RefCoords(1), RefCoords(2), RefCoords(3), 0) ! IPAR = 0 mode for E013
      IF (IER.LT.0) THEN
          IF (myid .EQ. 0) WRITE(*,*) 'EvalVel: ELE(IPAR=0) error'
          Found = .FALSE.
          RETURN
      END IF

      VelocityAtPoint = 0.0D0
      DO I=1,IDFL_vel ! Should be NNBAS = 27
        IG_vel = KDFG_vel(I)
        BasisFuncVal = DBAS(1, KDFL_vel(I), 1) ! KDFL_vel(I) is the local basis func index

        VelocityAtPoint(1) = VelocityAtPoint(1) + U1_global(IG_vel) * BasisFuncVal
        VelocityAtPoint(2) = VelocityAtPoint(2) + U2_global(IG_vel) * BasisFuncVal
        VelocityAtPoint(3) = VelocityAtPoint(3) + U3_global(IG_vel) * BasisFuncVal
      END DO

99999 CONTINUE
      END SUBROUTINE EvaluateVelocityAtPoint

FUNCTION TrilinearShapeFunc(NodeK, xi1, xi2, xi3) RESULT(N_val)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    INTEGER NodeK
    REAL*8 xi1, xi2, xi3, N_val
    REAL*8 NodeRefCoords(3,8) ! Predefined for nodes 1-8 of ref cube [-1,1]^3
    ! Example: NodeRefCoords(:,1) = (/-1.0D0, -1.0D0, -1.0D0/)
    !          NodeRefCoords(:,2) = (/ 1.0D0, -1.0D0, -1.0D0/) etc.
    ! ... (fill NodeRefCoords) ...

    N_val = 0.125D0 * (1.0D0 + NodeRefCoords(1,NodeK)*xi1) * &
                      (1.0D0 + NodeRefCoords(2,NodeK)*xi2) * &
                      (1.0D0 + NodeRefCoords(3,NodeK)*xi3)
END FUNCTION TrilinearShapeFunc