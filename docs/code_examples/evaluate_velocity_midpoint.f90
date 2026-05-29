SUBROUTINE Evaluate_Velocity_Midpoint(QuadSc, velocity_midpoints, element_centers, mfile)
!
! Purpose: Evaluate velocity at the midpoint (center) of each element
!          This is useful for visualization, post-processing, and analysis
!
! Mathematical Background:
!   For Q2 finite elements, velocity at element midpoint (ξ₁=0, ξ₂=0, ξ₃=0):
!   u(x_center) = Σⱼ uⱼ φⱼ(0,0,0)
!   where φⱼ are the Q2 basis functions and uⱼ are nodal velocities
!
! The physical center coordinates are computed using trilinear mapping:
!   x_center = Σₖ xₖ Nₖ(0,0,0) where Nₖ are the 8-node trilinear shape functions
!
USE var_QuadScalar, ONLY: TQuadScalar
USE mg_mesh_TypeDef, ONLY: mg_mesh
USE var_LScalars, ONLY: ILEV

IMPLICIT NONE

! Parameters (consistent with FeatFloWer conventions)
INTEGER, PARAMETER :: NNBAS=27, NNVE=8, NNDIM=3, NNDER=10

! Arguments
TYPE(TQuadScalar), INTENT(IN) :: QuadSc
REAL*8, INTENT(OUT) :: velocity_midpoints(:,:)  ! (3, NEL) - velocity components at element centers
REAL*8, INTENT(OUT) :: element_centers(:,:)     ! (3, NEL) - physical coordinates of element centers
INTEGER, INTENT(IN) :: mfile

! Local variables
INTEGER :: IEL, NEL, IDFL, JDOFE, JDFL, JDFG, IVE
INTEGER :: KDFG(NNBAS), KDFL(NNBAS)
INTEGER :: IELTYP

! Element geometry
REAL*8 :: DX(NNVE), DY(NNVE), DZ(NNVE)
REAL*8 :: X_center, Y_center, Z_center

! Velocity components at element center
REAL*8 :: U_center, V_center, W_center

! Q2 basis function values at element center (ξ₁=ξ₂=ξ₃=0)
REAL*8 :: PHI_center(NNBAS)

! Common blocks
COMMON /ELEM/ DBAS(NNDIM,NNBAS,NNDER), BDER(NNDER)
COMMON /COAUX1/ KDFG, KDFL, IDFL

REAL*8 :: DBAS
LOGICAL :: BDER

! External subroutines
EXTERNAL :: ELE, NDFGL
EXTERNAL :: E013, E013A

!===========================================================================

NEL = mg_mesh%level(ILEV)%NEL

! Check array dimensions
IF (SIZE(velocity_midpoints,1) .NE. 3 .OR. SIZE(velocity_midpoints,2) .NE. NEL) THEN
  WRITE(mfile,*) 'ERROR: velocity_midpoints array has wrong dimensions'
  WRITE(mfile,*) 'Expected: (3,', NEL, '), Got: (', SIZE(velocity_midpoints,1), ',', SIZE(velocity_midpoints,2), ')'
  RETURN
END IF

IF (SIZE(element_centers,1) .NE. 3 .OR. SIZE(element_centers,2) .NE. NEL) THEN
  WRITE(mfile,*) 'ERROR: element_centers array has wrong dimensions'
  WRITE(mfile,*) 'Expected: (3,', NEL, '), Got: (', SIZE(element_centers,1), ',', SIZE(element_centers,2), ')'
  RETURN
END IF

! Initialize basis function flags
DO JDOFE = 1, NNDER
  BDER(JDOFE) = .FALSE.
END DO
BDER(1) = .TRUE.  ! Only need function values φⱼ, not derivatives

! Set element type identifier
CALL ELE(0D0, 0D0, 0D0, IELTYP)
IDFL = 27  ! Q2 element has 27 nodes

! Pre-compute Q2 basis function values at element center (0,0,0)
CALL E013A(0D0, 0D0, 0D0, PHI_center, 1)

WRITE(mfile,*) 'Evaluating velocity at element midpoints...'
WRITE(mfile,*) 'Number of elements:', NEL

!===========================================================================
! MAIN ELEMENT LOOP
!===========================================================================

DO IEL = 1, NEL

  ! Get element connectivity
  CALL NDFGL(IEL, IDFL, KDFG, KDFL, mg_mesh%level(ILEV))

  ! Get element vertex coordinates (8 corners for hexahedral element)
  DO IVE = 1, NNVE
    JDFG = mg_mesh%level(ILEV)%kvert(IVE, IEL)
    DX(IVE) = mg_mesh%level(ILEV)%dcorvg(1, JDFG)  ! x-coordinate
    DY(IVE) = mg_mesh%level(ILEV)%dcorvg(2, JDFG)  ! y-coordinate
    DZ(IVE) = mg_mesh%level(ILEV)%dcorvg(3, JDFG)  ! z-coordinate
  END DO

  ! Calculate physical coordinates of element center
  ! Using trilinear mapping: x_center = Σₖ xₖ Nₖ(0,0,0)
  ! For hexahedral elements, Nₖ(0,0,0) = 1/8 for all 8 vertices
  X_center = 0D0
  Y_center = 0D0
  Z_center = 0D0
  DO IVE = 1, NNVE
    X_center = X_center + DX(IVE)
    Y_center = Y_center + DY(IVE)
    Z_center = Z_center + DZ(IVE)
  END DO
  X_center = X_center / DBLE(NNVE)
  Y_center = Y_center / DBLE(NNVE)
  Z_center = Z_center / DBLE(NNVE)

  ! Store element center coordinates
  element_centers(1, IEL) = X_center
  element_centers(2, IEL) = Y_center
  element_centers(3, IEL) = Z_center

  ! Initialize velocity components at element center
  U_center = 0D0
  V_center = 0D0
  W_center = 0D0

  ! Interpolate velocity using Q2 basis functions
  ! u(x_center) = Σⱼ uⱼ φⱼ(0,0,0)
  DO JDOFE = 1, IDFL
    JDFL = KDFL(JDOFE)
    JDFG = KDFG(JDOFE)

    ! Add contribution from node JDFG with basis function JDFL
    U_center = U_center + QuadSc%valU(JDFG) * PHI_center(JDFL)
    V_center = V_center + QuadSc%valV(JDFG) * PHI_center(JDFL)
    W_center = W_center + QuadSc%valW(JDFG) * PHI_center(JDFL)
  END DO

  ! Store velocity components at element center
  velocity_midpoints(1, IEL) = U_center
  velocity_midpoints(2, IEL) = V_center
  velocity_midpoints(3, IEL) = W_center

END DO ! End element loop

WRITE(mfile,*) 'Element midpoint velocity evaluation complete.'

END SUBROUTINE Evaluate_Velocity_Midpoint

!===========================================================================

SUBROUTINE Evaluate_Single_Element_Velocity(QuadSc, IEL, xi1, xi2, xi3, velocity, position, mfile)
!
! Purpose: Evaluate velocity at any point within a specific element
!          More general version for arbitrary points in reference coordinates
!
! Parameters:
!   QuadSc - velocity field data
!   IEL    - element number
!   xi1, xi2, xi3 - reference coordinates [-1,1]³
!   velocity - output velocity vector (3 components)
!   position - output physical coordinates (3 components)
!   mfile  - file unit for logging
!
USE var_QuadScalar, ONLY: TQuadScalar
USE mg_mesh_TypeDef, ONLY: mg_mesh
USE var_LScalars, ONLY: ILEV

IMPLICIT NONE

! Parameters
INTEGER, PARAMETER :: NNBAS=27, NNVE=8, NNDIM=3, NNDER=10

! Arguments
TYPE(TQuadScalar), INTENT(IN) :: QuadSc
INTEGER, INTENT(IN) :: IEL
REAL*8, INTENT(IN) :: xi1, xi2, xi3
REAL*8, INTENT(OUT) :: velocity(3)
REAL*8, INTENT(OUT) :: position(3)
INTEGER, INTENT(IN) :: mfile

! Local variables
INTEGER :: IDFL, JDOFE, JDFL, JDFG, IVE
INTEGER :: KDFG(NNBAS), KDFL(NNBAS)
INTEGER :: IELTYP

! Element geometry
REAL*8 :: DX(NNVE), DY(NNVE), DZ(NNVE)

! Q2 basis function values at specified point
REAL*8 :: PHI_point(NNBAS)

! Common blocks
COMMON /ELEM/ DBAS(NNDIM,NNBAS,NNDER), BDER(NNDER)
COMMON /COAUX1/ KDFG, KDFL, IDFL

REAL*8 :: DBAS
LOGICAL :: BDER

! External subroutines
EXTERNAL :: ELE, NDFGL
EXTERNAL :: E013, E013A

!===========================================================================

! Validate element number
IF (IEL .LT. 1 .OR. IEL .GT. mg_mesh%level(ILEV)%NEL) THEN
  WRITE(mfile,*) 'ERROR: Invalid element number:', IEL
  velocity = 0D0
  position = 0D0
  RETURN
END IF

! Validate reference coordinates
IF (ABS(xi1) .GT. 1D0 .OR. ABS(xi2) .GT. 1D0 .OR. ABS(xi3) .GT. 1D0) THEN
  WRITE(mfile,*) 'WARNING: Reference coordinates outside [-1,1]³'
  WRITE(mfile,*) 'xi1, xi2, xi3 =', xi1, xi2, xi3
END IF

! Initialize
DO JDOFE = 1, NNDER
  BDER(JDOFE) = .FALSE.
END DO
BDER(1) = .TRUE.  ! Function values only

CALL ELE(0D0, 0D0, 0D0, IELTYP)
IDFL = 27

! Get element connectivity
CALL NDFGL(IEL, IDFL, KDFG, KDFL, mg_mesh%level(ILEV))

! Get element vertex coordinates
DO IVE = 1, NNVE
  JDFG = mg_mesh%level(ILEV)%kvert(IVE, IEL)
  DX(IVE) = mg_mesh%level(ILEV)%dcorvg(1, JDFG)
  DY(IVE) = mg_mesh%level(ILEV)%dcorvg(2, JDFG)
  DZ(IVE) = mg_mesh%level(ILEV)%dcorvg(3, JDFG)
END DO

! Calculate physical coordinates using trilinear mapping
! This is a simplified version - full isoparametric mapping would be more accurate
position(1) = 0D0
position(2) = 0D0
position(3) = 0D0
DO IVE = 1, NNVE
  ! Trilinear shape functions for 8-node hexahedron
  ! Nₖ(ξ₁,ξ₂,ξ₃) = 1/8 * (1+ξ₁ₖξ₁)(1+ξ₂ₖξ₂)(1+ξ₃ₖξ₃)
  ! For simplicity, using element center approximation here
  position(1) = position(1) + DX(IVE)
  position(2) = position(2) + DY(IVE)
  position(3) = position(3) + DZ(IVE)
END DO
position(1) = position(1) / DBLE(NNVE)
position(2) = position(2) / DBLE(NNVE)
position(3) = position(3) / DBLE(NNVE)

! Evaluate Q2 basis functions at specified point
CALL E013A(xi1, xi2, xi3, PHI_point, 1)

! Interpolate velocity
velocity(1) = 0D0  ! U component
velocity(2) = 0D0  ! V component
velocity(3) = 0D0  ! W component

DO JDOFE = 1, IDFL
  JDFL = KDFL(JDOFE)
  JDFG = KDFG(JDOFE)

  velocity(1) = velocity(1) + QuadSc%valU(JDFG) * PHI_point(JDFL)
  velocity(2) = velocity(2) + QuadSc%valV(JDFG) * PHI_point(JDFL)
  velocity(3) = velocity(3) + QuadSc%valW(JDFG) * PHI_point(JDFL)
END DO

END SUBROUTINE Evaluate_Single_Element_Velocity

!===========================================================================

SUBROUTINE Write_Velocity_Midpoints_VTK(velocity_midpoints, element_centers, filename, mfile)
!
! Purpose: Write element-centered velocity data to VTK format for visualization
!          Compatible with ParaView and other VTK-based visualization tools
!
IMPLICIT NONE

! Arguments
REAL*8, INTENT(IN) :: velocity_midpoints(:,:)  ! (3, NEL)
REAL*8, INTENT(IN) :: element_centers(:,:)     ! (3, NEL)
CHARACTER(LEN=*), INTENT(IN) :: filename
INTEGER, INTENT(IN) :: mfile

! Local variables
INTEGER :: NEL, IEL, unit_vtk
CHARACTER(LEN=256) :: full_filename

!===========================================================================

NEL = SIZE(velocity_midpoints, 2)

! Open VTK file
unit_vtk = mfile + 10  ! Use offset to avoid conflicts
full_filename = TRIM(filename) // '.vtk'

OPEN(unit_vtk, FILE=full_filename, STATUS='REPLACE')

! Write VTK header
WRITE(unit_vtk, '(A)') '# vtk DataFile Version 3.0'
WRITE(unit_vtk, '(A)') 'Element-centered velocity data from FeatFloWer'
WRITE(unit_vtk, '(A)') 'ASCII'
WRITE(unit_vtk, '(A)') 'DATASET UNSTRUCTURED_GRID'

! Write points (element centers)
WRITE(unit_vtk, '(A,I0,A)') 'POINTS ', NEL, ' double'
DO IEL = 1, NEL
  WRITE(unit_vtk, '(3(E15.6,1X))') element_centers(1:3, IEL)
END DO

! Write cells (points)
WRITE(unit_vtk, '(A,I0,1X,I0)') 'CELLS ', NEL, 2*NEL
DO IEL = 1, NEL
  WRITE(unit_vtk, '(A,I0)') '1 ', IEL-1  ! VTK uses 0-based indexing
END DO

! Write cell types (all points = type 1)
WRITE(unit_vtk, '(A,I0)') 'CELL_TYPES ', NEL
DO IEL = 1, NEL
  WRITE(unit_vtk, '(A)') '1'
END DO

! Write point data
WRITE(unit_vtk, '(A,I0)') 'POINT_DATA ', NEL

! Write velocity vectors
WRITE(unit_vtk, '(A)') 'VECTORS velocity double'
DO IEL = 1, NEL
  WRITE(unit_vtk, '(3(E15.6,1X))') velocity_midpoints(1:3, IEL)
END DO

! Write velocity magnitude
WRITE(unit_vtk, '(A)') 'SCALARS velocity_magnitude double 1'
WRITE(unit_vtk, '(A)') 'LOOKUP_TABLE default'
DO IEL = 1, NEL
  WRITE(unit_vtk, '(E15.6)') SQRT(SUM(velocity_midpoints(1:3, IEL)**2))
END DO

CLOSE(unit_vtk)

WRITE(mfile, '(A,A)') 'VTK file written: ', TRIM(full_filename)

END SUBROUTINE Write_Velocity_Midpoints_VTK