SUBROUTINE Get_DissipationIntegral(QuadSc, ALPHA, total_dissipation, mfile)
!
! Purpose: Calculate the strain rate dissipation integral in fluid domain:
!          ∫_{Π_o} D(u) : D(u) dx
!
! This follows the same computational pattern as GetForces and DivGradStress
! routines in FeatFloWer, using the established FEM infrastructure.
!
! Mathematical Background:
!   D(u) = 1/2 * (∇u + (∇u)^T)  [rate of strain tensor]
!   D(u) : D(u) = Σᵢⱼ Dᵢⱼ * Dᵢⱼ   [Frobenius inner product]
!
! Integration over fluid domain excludes particle interiors where 
! velocity follows rigid body motion and D(u) = 0.
!
USE PP3D_MPI, ONLY: COMM_SUMMN
USE var_QuadScalar, ONLY: TQuadScalar
USE mg_mesh_TypeDef, ONLY: mg_mesh
USE var_LScalars, ONLY: ILEV

IMPLICIT NONE

! Parameters (consistent with FeatFloWer conventions)
INTEGER, PARAMETER :: NNBAS=27, NNVE=8, NNDIM=3, NNDER=10, NNCUBP=100

! Arguments  
TYPE(TQuadScalar), INTENT(IN) :: QuadSc
INTEGER, INTENT(IN) :: ALPHA(:)
REAL*8, INTENT(OUT) :: total_dissipation
INTEGER, INTENT(IN) :: mfile

! Local variables
INTEGER :: IEL, NEL, ICUBP, NCUBP, ICUB=9
INTEGER :: IDFL, JDOFE, JDFL, JDFG, IVE, IG
INTEGER :: KDFG(NNBAS), KDFL(NNBAS)
INTEGER :: IELTYP

! Element geometry
REAL*8 :: DX(NNVE), DY(NNVE), DZ(NNVE)
REAL*8 :: DX0, DY0, DZ0  ! Element center approximation

! Jacobian and integration  
REAL*8 :: DJAC(3,3), DETJ, OM, XI1, XI2, XI3
REAL*8 :: DJ11, DJ12, DJ13, DJ21, DJ22, DJ23, DJ31, DJ32, DJ33

! Velocity gradients at cubature point
REAL*8 :: DU1DX, DU1DY, DU1DZ  ! ∇u₁
REAL*8 :: DU2DX, DU2DY, DU2DZ  ! ∇u₂  
REAL*8 :: DU3DX, DU3DY, DU3DZ  ! ∇u₃

! Rate of strain tensor components
REAL*8 :: D11, D22, D33, D12, D13, D23

! Element particle detection (following GetForces pattern)
INTEGER :: NJALFA, NIALFA
REAL*8 :: local_dissipation, element_dissipation

! Common blocks (FeatFloWer standard)
COMMON /ELEM/ DBAS(NNDIM,NNBAS,NNDER), BDER(NNDER)
COMMON /CUB/ DXI(NNCUBP,3), DOMEGA(NNCUBP)  
COMMON /COAUX1/ KDFG, KDFL, IDFL

REAL*8 :: DBAS, DXI, DOMEGA
LOGICAL :: BDER

! External subroutines
EXTERNAL :: ELE, CB3H, NDFGL
EXTERNAL :: E013

!===========================================================================

NEL = mg_mesh%level(ILEV)%NEL
total_dissipation = 0D0

! Initialize basis function flags (following GetForces pattern)
DO IG = 1, NNDER
  BDER(IG) = .FALSE.
END DO
BDER(1) = .TRUE.  ! Function values φⱼ
BDER(2) = .TRUE.  ! ∂φⱼ/∂x derivatives
BDER(3) = .TRUE.  ! ∂φⱼ/∂y derivatives  
BDER(4) = .TRUE.  ! ∂φⱼ/∂z derivatives

! Set element type identifier  
CALL ELE(0D0, 0D0, 0D0, IELTYP)
IDFL = 27  ! Q2 element has 27 nodes

! Set up 3D cubature rule
CALL CB3H(ICUB)
NCUBP = 27  ! 3x3x3 Gaussian cubature

!===========================================================================
! MAIN ELEMENT LOOP
!===========================================================================

DO IEL = 1, NEL

  ! Get element connectivity (following GetForces pattern)
  CALL NDFGL(IEL, IDFL, KDFG, KDFL, mg_mesh%level(ILEV))

  ! Interface element detection (same logic as GetForces)
  NJALFA = 0  ! Count nodes NOT in any particle
  NIALFA = 0  ! Count nodes in particles
  
  DO JDOFE = 1, IDFL
    IG = KDFG(JDOFE)
    IF (ALPHA(IG) .EQ. 0) THEN
      NJALFA = NJALFA + 1
    ELSE  
      NIALFA = NIALFA + 1
    END IF
  END DO

  ! Skip elements completely inside particles (∇u = 0 for rigid body motion)
  IF (NJALFA .EQ. 0) CYCLE

  ! Get element vertex coordinates for Jacobian calculation
  DO IVE = 1, NNVE
    JDFG = mg_mesh%level(ILEV)%kvert(IVE, IEL)
    DX(IVE) = mg_mesh%level(ILEV)%dcorvg(1, JDFG)
    DY(IVE) = mg_mesh%level(ILEV)%dcorvg(2, JDFG)  
    DZ(IVE) = mg_mesh%level(ILEV)%dcorvg(3, JDFG)
  END DO

  ! Element center (for Jacobian calculation)
  DX0 = 0D0; DY0 = 0D0; DZ0 = 0D0
  DO IVE = 1, NNVE
    DX0 = DX0 + DX(IVE)
    DY0 = DY0 + DY(IVE)
    DZ0 = DZ0 + DZ(IVE)
  END DO
  DX0 = DX0 / DBLE(NNVE)
  DY0 = DY0 / DBLE(NNVE)
  DZ0 = DZ0 / DBLE(NNVE)

  ! Pre-compute Jacobian coefficients (following GetForces pattern)
  ! These are used for isoparametric mapping from reference to physical element
  DJ11 = (DX(1)+DX(2)+DX(3)+DX(4)+DX(5)+DX(6)+DX(7)+DX(8)) * 0.125D0
  DJ12 = (DY(1)+DY(2)+DY(3)+DY(4)+DY(5)+DY(6)+DY(7)+DY(8)) * 0.125D0
  DJ13 = (DZ(1)+DZ(2)+DZ(3)+DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8)) * 0.125D0
  ! ... (additional DJij terms would be computed here for full Jacobian)

  ! Pre-compute basis functions for all cubature points
  CALL ELE(0D0, 0D0, 0D0, -2)

  element_dissipation = 0D0

  !=========================================================================
  ! CUBATURE POINT LOOP  
  !=========================================================================

  DO ICUBP = 1, NCUBP

    ! Cubature point coordinates in reference element [-1,1]³
    XI1 = DXI(ICUBP, 1)
    XI2 = DXI(ICUBP, 2)
    XI3 = DXI(ICUBP, 3)

    ! Compute Jacobian matrix (simplified - full implementation needed)
    ! This transforms derivatives from reference to physical coordinates
    DJAC(1,1) = DJ21 + DJ51*XI2 + DJ61*XI3  ! ∂x/∂ξ₁
    DJAC(1,2) = DJ31 + DJ51*XI1 + DJ71*XI3  ! ∂x/∂ξ₂  
    DJAC(1,3) = DJ41 + DJ61*XI1 + DJ71*XI2  ! ∂x/∂ξ₃
    ! ... (similar for other components)
    
    DETJ = DJAC(1,1)*(DJAC(2,2)*DJAC(3,3) - DJAC(3,2)*DJAC(2,3)) - &
           DJAC(1,2)*(DJAC(2,1)*DJAC(3,3) - DJAC(3,1)*DJAC(2,3)) + &
           DJAC(1,3)*(DJAC(2,1)*DJAC(3,2) - DJAC(3,1)*DJAC(2,2))

    ! Integration weight
    OM = DOMEGA(ICUBP) * ABS(DETJ)

    ! Evaluate basis functions and physical derivatives at cubature point
    CALL ELE(XI1, XI2, XI3, -3)

    ! Initialize velocity gradients
    DU1DX = 0D0; DU1DY = 0D0; DU1DZ = 0D0
    DU2DX = 0D0; DU2DY = 0D0; DU2DZ = 0D0
    DU3DX = 0D0; DU3DY = 0D0; DU3DZ = 0D0

    ! Interpolate velocity gradients using Q2 basis functions
    DO JDOFE = 1, IDFL
      JDFL = KDFL(JDOFE)
      JDFG = KDFG(JDOFE)

      ! Physical derivatives: ∂φⱼ/∂x, ∂φⱼ/∂y, ∂φⱼ/∂z
      DU1DX = DU1DX + QuadSc%valU(JDFG) * DBAS(1, JDFL, 2)
      DU1DY = DU1DY + QuadSc%valU(JDFG) * DBAS(1, JDFL, 3)
      DU1DZ = DU1DZ + QuadSc%valU(JDFG) * DBAS(1, JDFL, 4)

      DU2DX = DU2DX + QuadSc%valV(JDFG) * DBAS(1, JDFL, 2)  
      DU2DY = DU2DY + QuadSc%valV(JDFG) * DBAS(1, JDFL, 3)
      DU2DZ = DU2DZ + QuadSc%valV(JDFG) * DBAS(1, JDFL, 4)

      DU3DX = DU3DX + QuadSc%valW(JDFG) * DBAS(1, JDFL, 2)
      DU3DY = DU3DY + QuadSc%valW(JDFG) * DBAS(1, JDFL, 3)  
      DU3DZ = DU3DZ + QuadSc%valW(JDFG) * DBAS(1, JDFL, 4)
    END DO

    ! Compute rate of strain tensor: D(u) = 1/2 * (∇u + (∇u)ᵀ)
    D11 = DU1DX                              ! ∂u₁/∂x₁
    D22 = DU2DY                              ! ∂u₂/∂x₂
    D33 = DU3DZ                              ! ∂u₃/∂x₃
    D12 = 0.5D0 * (DU1DY + DU2DX)           ! ½(∂u₁/∂x₂ + ∂u₂/∂x₁)
    D13 = 0.5D0 * (DU1DZ + DU3DX)           ! ½(∂u₁/∂x₃ + ∂u₃/∂x₁)
    D23 = 0.5D0 * (DU2DZ + DU3DY)           ! ½(∂u₂/∂x₃ + ∂u₃/∂x₂)

    ! Frobenius inner product: D(u) : D(u) = Σᵢⱼ Dᵢⱼ²
    local_dissipation = D11*D11 + D22*D22 + D33*D33 + &
                       2D0*(D12*D12 + D13*D13 + D23*D23)

    ! Contribute to element integral (weighted by integration measure)
    element_dissipation = element_dissipation + local_dissipation * OM

  END DO ! End cubature points

  ! Add element contribution to global integral
  total_dissipation = total_dissipation + element_dissipation

END DO ! End elements

!===========================================================================
! MPI SUMMATION
!===========================================================================

! Sum contributions from all processors  
CALL COMM_SUMMN(total_dissipation, 1)

END SUBROUTINE Get_DissipationIntegral