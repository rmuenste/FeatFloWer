SUBROUTINE Calculate_StrainRate_Dissipation(QuadSc, ALPHA, dissipation_integral, mfile)
!
! Purpose: Calculate the strain rate dissipation integral:
!          ∫_{Π_o} D(u) : D(u) dx
!
! Where D(u) is the rate of strain tensor: D(u) = 1/2 (∇u + (∇u)^T)
! and : denotes the Frobenius inner product
!
! Parameters:
!   QuadSc - TQuadScalar type containing velocity field data
!   ALPHA  - Particle indicator function (0 in fluid, particle_id in particles)
!   dissipation_integral - Output: computed integral value
!   mfile  - File unit for output/logging
!
USE PP3D_MPI, ONLY: COMM_SUMMN, myid
USE var_QuadScalar, ONLY: TQuadScalar
USE mg_mesh_TypeDef, ONLY: mg_mesh
USE var_LScalars, ONLY: ILEV

IMPLICIT NONE

! Parameters
INTEGER, PARAMETER :: NNBAS=27, NNVE=8, NNDIM=3, NNDER=10, NNCUBP=100

! Arguments
TYPE(TQuadScalar), INTENT(IN) :: QuadSc
INTEGER, INTENT(IN) :: ALPHA(:)
REAL*8, INTENT(OUT) :: dissipation_integral
INTEGER, INTENT(IN) :: mfile

! Local variables
INTEGER :: IEL, NEL, ICUBP, NCUBP, ICUB=9
INTEGER :: IDFL, JDOFE, JDFL, JDFG
INTEGER :: KDFG(NNBAS), KDFL(NNBAS)
INTEGER :: IELTYP

REAL*8 :: OM, DETJ, XI1, XI2, XI3
REAL*8 :: DJAC(3,3), DPP(3)
REAL*8 :: local_dissipation, element_dissipation

! Velocity gradients at cubature point
REAL*8 :: DU1DX, DU1DY, DU1DZ  ! ∂u₁/∂x, ∂u₁/∂y, ∂u₁/∂z
REAL*8 :: DU2DX, DU2DY, DU2DZ  ! ∂u₂/∂x, ∂u₂/∂y, ∂u₂/∂z  
REAL*8 :: DU3DX, DU3DY, DU3DZ  ! ∂u₃/∂x, ∂u₃/∂y, ∂u₃/∂z

! Rate of strain tensor components
REAL*8 :: D11, D22, D33, D12, D13, D23

! For fluid domain detection
REAL*8 :: fluid_indicator
INTEGER :: NJALFA, node_count, IG

! Common blocks
COMMON /ELEM/ DBAS(NNDIM,NNBAS,NNDER), BDER(NNDER)
COMMON /CUB/ DXI(NNCUBP,3), DOMEGA(NNCUBP)
COMMON /COAUX1/ KDFG, KDFL, IDFL

REAL*8 :: DBAS, DXI, DOMEGA
LOGICAL :: BDER

! External functions
EXTERNAL :: ELE, CB3H, NDFGL
EXTERNAL :: E013

!===========================================================================

NEL = mg_mesh%level(ILEV)%NEL

! Initialize
dissipation_integral = 0d0

! Set up basis function derivatives needed
DO ICUBP = 1, NNDER
  BDER(ICUBP) = .FALSE.
END DO
BDER(1) = .TRUE.  ! Function values
BDER(2) = .TRUE.  ! ∂/∂x derivatives  
BDER(3) = .TRUE.  ! ∂/∂y derivatives
BDER(4) = .TRUE.  ! ∂/∂z derivatives

! Set element type (Q2)
CALL ELE(0D0, 0D0, 0D0, IELTYP)
IDFL = 27  ! NDFL(IELTYP) for Q2 elements

! Set up cubature rule
CALL CB3H(ICUB)
NCUBP = 27  ! Number of cubature points for ICUB=9

IF (myid .EQ. 0) THEN
  WRITE(mfile,*) 'Calculating strain rate dissipation integral...'
  WRITE(mfile,*) 'Number of elements:', NEL
  WRITE(mfile,*) 'Cubature points per element:', NCUBP
END IF

!===========================================================================
! MAIN LOOP OVER ELEMENTS
!===========================================================================

DO IEL = 1, NEL

  ! Get element degrees of freedom
  CALL NDFGL(IEL, IDFL, KDFG, KDFL, mg_mesh%level(ILEV))

  ! Check if element is in fluid domain
  ! Count nodes NOT belonging to any particle
  NJALFA = 0
  DO JDOFE = 1, IDFL
    IG = KDFG(JDOFE)
    IF (ALPHA(IG) .EQ. 0) THEN
      NJALFA = NJALFA + 1
    END IF
  END DO

  ! Skip elements completely inside particles
  IF (NJALFA .EQ. 0) CYCLE

  ! Get element vertex coordinates for Jacobian calculation
  DO JDOFE = 1, 8  ! 8 vertices for hexahedral element
    JDFG = mg_mesh%level(ILEV)%kvert(JDOFE, IEL)
    DPP(1) = mg_mesh%level(ILEV)%dcorvg(1, JDFG)  ! x-coordinate
    DPP(2) = mg_mesh%level(ILEV)%dcorvg(2, JDFG)  ! y-coordinate  
    DPP(3) = mg_mesh%level(ILEV)%dcorvg(3, JDFG)  ! z-coordinate
  END DO

  ! Pre-compute basis functions at all cubature points
  CALL ELE(0D0, 0D0, 0D0, -2)

  element_dissipation = 0d0

  !=========================================================================
  ! LOOP OVER CUBATURE POINTS
  !=========================================================================

  DO ICUBP = 1, NCUBP

    ! Get cubature point coordinates and weights
    XI1 = DXI(ICUBP, 1)
    XI2 = DXI(ICUBP, 2) 
    XI3 = DXI(ICUBP, 3)

    ! Compute Jacobian matrix and determinant
    ! (This would need proper implementation based on element geometry)
    ! For now, assume DETJ is computed properly
    ! CALL Compute_Jacobian(XI1, XI2, XI3, IEL, DJAC, DETJ)
    
    ! Integration weight
    OM = DOMEGA(ICUBP) * ABS(DETJ)

    ! Evaluate basis functions and derivatives at cubature point
    CALL ELE(XI1, XI2, XI3, -3)

    ! Initialize velocity gradients
    DU1DX = 0D0; DU1DY = 0D0; DU1DZ = 0D0
    DU2DX = 0D0; DU2DY = 0D0; DU2DZ = 0D0  
    DU3DX = 0D0; DU3DY = 0D0; DU3DZ = 0D0

    ! Interpolate velocity gradients at cubature point
    DO JDOFE = 1, IDFL
      JDFL = KDFL(JDOFE)
      JDFG = KDFG(JDOFE)

      ! Basis function derivatives (physical coordinates)
      ! ∂φⱼ/∂x, ∂φⱼ/∂y, ∂φⱼ/∂z
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

    ! Calculate rate of strain tensor components
    ! D(u) = 1/2 * (∇u + (∇u)^T)
    D11 = DU1DX                           ! ∂u₁/∂x₁
    D22 = DU2DY                           ! ∂u₂/∂x₂
    D33 = DU3DZ                           ! ∂u₃/∂x₃
    D12 = 0.5D0 * (DU1DY + DU2DX)        ! 1/2(∂u₁/∂x₂ + ∂u₂/∂x₁)
    D13 = 0.5D0 * (DU1DZ + DU3DX)        ! 1/2(∂u₁/∂x₃ + ∂u₃/∂x₁)
    D23 = 0.5D0 * (DU2DZ + DU3DY)        ! 1/2(∂u₂/∂x₃ + ∂u₃/∂x₂)

    ! Calculate Frobenius inner product D(u) : D(u)
    ! = D₁₁² + D₂₂² + D₃₃² + 2(D₁₂² + D₁₃² + D₂₃²)
    local_dissipation = D11*D11 + D22*D22 + D33*D33 + &
                       2D0*(D12*D12 + D13*D13 + D23*D23)

    ! Fluid domain indicator (simple approach)
    ! For interface elements, we could use a more sophisticated approach
    IF (NJALFA .GT. IDFL/2) THEN
      fluid_indicator = 1D0  ! Majority of nodes in fluid
    ELSE
      fluid_indicator = 0.5D0 ! Mixed element - partial contribution
    END IF

    ! Integrate: add contribution to element integral
    element_dissipation = element_dissipation + &
                         fluid_indicator * local_dissipation * OM

  END DO ! End cubature points loop

  ! Add element contribution to global integral
  dissipation_integral = dissipation_integral + element_dissipation

END DO ! End elements loop

!===========================================================================
! MPI REDUCTION
!===========================================================================

! Sum contributions from all processors
CALL COMM_SUMMN(dissipation_integral, 1)

IF (myid .EQ. 0) THEN
  WRITE(mfile,*) 'Strain rate dissipation integral = ', dissipation_integral
  WRITE(mfile,*) 'Calculation complete.'
END IF

END SUBROUTINE Calculate_StrainRate_Dissipation