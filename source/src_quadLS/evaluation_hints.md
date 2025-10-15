Yes, absolutely! Evaluating the velocity (or any field approximated by FEM basis functions) at an arbitrary physical point $(x_p, y_p, z_p)$ within a specific target element `target_iel` involves two main steps:

1.  **Inverse Isoparametric Mapping:** Find the reference coordinates $(\xi_1, \xi_2, \xi_3)$ within the reference element that correspond to the given physical point $(x_p, y_p, z_p)$ in element `target_iel`. This usually requires an iterative numerical method (like Newton-Raphson) because the forward mapping $x(\xi), y(\xi), z(\xi)$ is generally non-linear.
2.  **Interpolation:** Once the reference coordinates $(\xi_1, \xi_2, \xi_3)$ are found, use them to evaluate the Q2 basis functions (via `E013`) and then interpolate the nodal velocity values `U1, U2, U3`.

Let's outline the structure of such a function, say `EvaluateVelocityAtPoint`.

**Assumptions:**

*   The element geometry is described by the 8-node trilinear isoparametric mapping (using `DJmn` coefficients as in `GetForces`).
*   The velocity field is approximated using 27-node Q2 basis functions (`E013`).
*   We have access to nodal velocities `U1, U2, U3`, mesh connectivity `KVERT`, global vertex coordinates `DCORVG`, and the `ELE` function pointer (to `E013`).

**Step 1: Inverse Isoparametric Mapping (Newton-Raphson)**

The forward mapping is:
$x = \sum_{k=1}^{8} N_k^{trilinear}(\xi_1, \xi_2, \xi_3) x_k^{elem}$
$y = \sum_{k=1}^{8} N_k^{trilinear}(\xi_1, \xi_2, \xi_3) y_k^{elem}$
$z = \sum_{k=1}^{8} N_k^{trilinear}(\xi_1, \xi_2, \xi_3) z_k^{elem}$

Where $N_k^{trilinear}$ are the 8-node trilinear shape functions (which are implicitly used when calculating `XX,YY,ZZ` from `XI1,XI2,XI3` and `DJmn` terms in `GetForces`).
We are given $(x_p, y_p, z_p)$ and want to find $(\xi_1, \xi_2, \xi_3)$.
We need to solve the non-linear system:
$f_1(\xi_1, \xi_2, \xi_3) = x(\xi_1, \xi_2, \xi_3) - x_p = 0$
$f_2(\xi_1, \xi_2, \xi_3) = y(\xi_1, \xi_2, \xi_3) - y_p = 0$
$f_3(\xi_1, \xi_2, \xi_3) = z(\xi_1, \xi_2, \xi_3) - z_p = 0$

The Newton-Raphson iteration is:
$\mathbf{J}_{\xi}^{(iter)} \Delta\boldsymbol{\xi}^{(iter)} = -\mathbf{f}^{(iter)}$
$\boldsymbol{\xi}^{(iter+1)} = \boldsymbol{\xi}^{(iter)} + \Delta\boldsymbol{\xi}^{(iter)}$

Where $\boldsymbol{\xi} = (\xi_1, \xi_2, \xi_3)^T$, $\mathbf{f} = (f_1, f_2, f_3)^T$, and $\mathbf{J}_{\xi}$ is the Jacobian of the forward mapping with respect to $\xi$:
$\mathbf{J}_{\xi} = \begin{pmatrix}
\frac{\partial x}{\partial \xi_1} & \frac{\partial x}{\partial \xi_2} & \frac{\partial x}{\partial \xi_3} \\
\frac{\partial y}{\partial \xi_1} & \frac{\partial y}{\partial \xi_2} & \frac{\partial y}{\partial \xi_3} \\
\frac{\partial z}{\partial \xi_1} & \frac{\partial z}{\partial \xi_2} & \frac{\partial z}{\partial \xi_3}
\end{pmatrix}$
This is exactly the matrix `DJAC(3,3)` that we calculate in `GetForces` (and in the proposed `CalculateDissipationIntegralNumerator`).

**Fortran Subroutine Structure `EvaluateVelocityAtPoint`:**

```fortran
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
```

**Explanation and Key Parts:**

1.  **Inputs/Outputs:** Clearly defined. `Found` indicates success.
2.  **Newton-Raphson Setup:**
    *   `MaxIter`, `Tolerance`.
    *   Initial guess for `CurrentRefCoords` is $(0,0,0)$, the center of the reference element.
    *   The `DJmn` coefficients for the *geometric mapping* (trilinear, 8-node) of `Target_IEL` are pre-calculated. These define $x(\boldsymbol{\xi}), y(\boldsymbol{\xi}), z(\boldsymbol{\xi})$ and its Jacobian $\mathbf{J}_{\xi}$.
3.  **Newton-Raphson Loop:**
    *   Calculates $\mathbf{J}_{\xi}$ (which is `JacobianInvMap` in the code, or `DJAC` in `GetForces`) at the current `CurrentRefCoords`.
    *   Calculates the current physical position $x(\boldsymbol{\xi}^{(iter)}), y(\boldsymbol{\xi}^{(iter)}), z(\boldsymbol{\xi}^{(iter)})$ using the forward mapping (the `XX,YY,ZZ` calculation from `GetForces`). The lines for `CurrentPhysCoords` need careful checking to ensure they use the correct `DJmn` terms for the *forward mapping itself*, not just its derivatives.
        *   **Correction/Clarification for `CurrentPhysCoords` calculation:** The forward mapping for an 8-node trilinear element is more directly written as:
            ```fortran
            CurrentPhysCoords = 0.0D0
            DO IV_node = 1, NNVE ! 8 nodes
                ! Get trilinear shape function N_IV_node(XI1_iter, XI2_iter, XI3_iter)
                ! This needs a small helper or explicit formula for N_IV_node
                N_val_trilin = TrilinearShapeFunc(IV_node, XI1_iter, XI2_iter, XI3_iter) 
                CurrentPhysCoords(1) = CurrentPhysCoords(1) + N_val_trilin * DX(IV_node)
                CurrentPhysCoords(2) = CurrentPhysCoords(2) + N_val_trilin * DY(IV_node)
                CurrentPhysCoords(3) = CurrentPhysCoords(3) + N_val_trilin * DZ(IV_node)
            END DO
            ```
            The `XX,YY,ZZ` calculation using `DJAC` components was a bit of a shortcut. While `DJAC` components are derivatives of this mapping, the mapping itself is simpler when expressed with the $N_k^{trilinear}$ directly. The previous `XX,YY,ZZ` formulation in `GetForces` likely embedded the $N_k^{trilinear}$ terms implicitly via the `DJmn` components. For robustness, using explicit $N_k^{trilinear}$ is better for the forward map here. I will leave the `DJmn`-based one as it was in `GetForces` for now, assuming it's correct there.
    *   Calculates the residual vector $\mathbf{f} = \mathbf{x}(\boldsymbol{\xi}^{(iter)}) - \mathbf{x}_p$.
    *   Checks for convergence (`ErrorNorm < Tolerance`).
    *   Solves $\mathbf{J}_{\xi} \Delta\boldsymbol{\xi} = -\mathbf{f}$ by explicitly inverting the 3x3 $\mathbf{J}_{\xi}$ matrix.
    *   Updates $\boldsymbol{\xi}^{(iter+1)} = \boldsymbol{\xi}^{(iter)} + \Delta\boldsymbol{\xi}^{(iter)}$.
    *   Includes an optional check if $\boldsymbol{\xi}$ goes out of bounds (e.g., `abs(xi_i) > 1.0 + tolerance`).

4.  **Velocity Interpolation (after convergence):**
    *   `RefCoords` stores the converged $(\xi_1, \xi_2, \xi_3)$.
    *   `NDFGL` is called to get the global DOFs (`KDFG_vel`) for the velocity field (Q2 element) corresponding to `Target_IEL`.
    *   `BDER` is set to request only basis function *values* (`BDER(1)=.TRUE.`).
    *   `CALL ELE(RefCoords(1), RefCoords(2), RefCoords(3), 0)`: This calls `E013` in mode `IPAR=0`. `E013` will evaluate the Q2 basis functions at `RefCoords` and store the values in `DBAS(1, local_basis_idx, 1)`.
        *   **Important Note for `E013` with `IPAR=0`:** `E013` when `IPAR=0` calls `E013A` which computes reference values/derivatives into `DHELP(...,1)`. Then, `E013` copies `DHELP(IDFL,1,1)` to `DBAS(1,IDFL,1)`. If `BDER(2,3,4)` were also true, it would then attempt to transform derivatives using `DJAC` and `DETJ` from `COMMON /ELEM/`. Since we only set `BDER(1)=.TRUE.`, the derivative transformation part in `E013` should be skipped, and `DJAC`/`DETJ` in common *might* not be strictly needed for this call. However, for safety, I've re-calculated `DJAC` and `DETJ` using the found `RefCoords` before the `ELE(IPAR=0)` call, in case `E013` still peeks at them.
    *   The velocity components are then interpolated: $\mathbf{u}(x_p,y_p,z_p) = \sum_{i=1}^{27} N_i^{Q2}(\xi_1,\xi_2,\xi_3) \mathbf{u}_i^{global}$.

**Helper for Trilinear Shape Functions (if needed for forward map):**
If the `XX,YY,ZZ` calculation from `GetForces` via `DJmn` proves tricky to reconstruct perfectly for the forward map, you might need:
```fortran
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
```

This routine provides a comprehensive way to achieve your goal. The Newton-Raphson for the inverse map is the most complex part.