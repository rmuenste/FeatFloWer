!=========================================================================
! QuadSc_corrections.f90
!
! Solution correction and RHS assembly functions for QuadScalar module
! Extracted from QuadSc_main.f90 for better code organization
!=========================================================================
!
!=========================================================================
SUBROUTINE Velocity_Correction()
  INTEGER i

  ! *** Update of U = U~ - k M^-1 B P

  ILEV = NLMAX
  qlMat     => mg_qlMat(ILEV)
  BXMat     => mg_BXMat(ILEV)%a
  BYMat     => mg_BYMat(ILEV)%a
  BZMat     => mg_BZMat(ILEV)%a
  MlRhoPmat => mg_MlRhoPmat(ILEV)%a

  CALL B_Mul_U(qlMat%ColA,qlMAt%LdA,BXMat,BYMat,BZMat,LinSc%ValP(NLMAX)%x,&
    QuadSc%defU,QuadSc%defV,QuadSc%defW,QuadSc%ndof,-TSTEP,0d0)

  CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)

  CALL Boundary_QuadScalar_Def()

  DO I=1,QuadSc%ndof
  QuadSc%valU(i) = QuadSc%valU(i) - QuadSc%defU(i)/MlRhoPmat(i)
  QuadSc%valV(i) = QuadSc%valV(i) - QuadSc%defV(i)/MlRhoPmat(i)
  QuadSc%valW(i) = QuadSc%valW(i) - QuadSc%defW(i)/MlRhoPmat(i)
  END DO

END SUBROUTINE Velocity_Correction
!=========================================================================
!
!=========================================================================
SUBROUTINE Pressure_Correction()
  INTEGER I
  real*8 dR
  REAL*8 daux

  if (GAMMA.gt.0d0) THEN
   CALL SETLEV(2)
   P1iMMat   => mg_P1iMMat(NLMAX)%a

   daux = THSTEP*(GAMMA + Properties%Viscosity(1))
   CALL AddDiffPrec(P1iMMat,LinSc%rhsP(NLMAX)%x,LinSc%ValP(NLMAX)%x,knel(nlmax),daux)
  END IF

  DO I=1,lMat%nu
   LinSc%valP(NLMAX)%x(i) = LinSc%valP(NLMAX)%x(i) + LinSc%valP_old(i)
   LinSc%P_new(i) = 1.5d0*LinSc%valP(NLMAX)%x(i) - 0.5d0*LinSc%valP_old(i)
  END DO

END SUBROUTINE Pressure_Correction
!=========================================================================
!
!=========================================================================
SUBROUTINE AddPressureGradient()
  INTEGER I,J,IEL
  REAL*8 ddx,ddy,ddz,ddp

  ILEV = NLMAX
  qlMat    => mg_qlMat(ILEV)
  BXMat    => mg_BXMat(ILEV)%a
  BYMat    => mg_BYMat(ILEV)%a
  BZMat    => mg_BZMat(ILEV)%a
  CALL B_Mul_U(qlMat%ColA,qlMAt%LdA,BXMat,BYMat,BZMat,LinSc%valP(NLMAX)%x,&
    QuadSc%defU,QuadSc%defV,QuadSc%defW,QuadSc%ndof,TSTEP,1d0)

END SUBROUTINE AddPressureGradient
!=========================================================================
!
!=========================================================================
SUBROUTINE AddPressureGradientWithJump()
  INTEGER I,J,IEL,ivt
  REAL*8 ddx,ddy,ddz,ddp,P(3)

  IF (.not.allocated(dPeriodicVector)) ALLOCATE(dPeriodicVector(QuadSc%ndof))
  dPeriodicVector = 0d0

  DO i=1,SIZE(LinSc%AuxP(NLMAX)%x)
   IF (MOD(i,4).EQ.1) then
    LinSc%AuxP(NLMAX)%x(i) = 1d1
   ELSE
    LinSc%AuxP(NLMAX)%x(i) = 0d0
   END IF
  END DO

  CALL B_Mul_U(qlMat%ColA,qlMAt%LdA,BXMat,BYMat,BZMat,LinSc%AuxP(NLMAX)%x,&
       QuadSc%AuxU,QuadSc%AuxV,QuadSc%AuxW,QuadSc%ndof,+1d0,0d0)

  DO I=1,QuadSc%ndof
   IF ((abs(myQ2Coor(1,i)).LT.+1e-1).and.(myQ2Coor(2,i).LT.+0d0).and.myid.eq.8) THEN
    dPeriodicVector(i) = QuadSc%auxU(i)
   ELSE
    dPeriodicVector(i) = 0d0
   END IF
  END DO

  CALL B_Mul_U(qlMat%ColA,qlMAt%LdA,BXMat,BYMat,BZMat,LinSc%valP(NLMAX)%x,&
  QuadSc%defU,QuadSc%defV,QuadSc%defW,QuadSc%ndof,TSTEP,1d0)

  DO I=1,QuadSc%ndof
   QuadSc%defU(i) = QuadSc%defU(i) + TSTEP*dPeriodicVector(i)
  END DO

END SUBROUTINE AddPressureGradientWithJump
!=========================================================================
!
!=========================================================================
SUBROUTINE AddPeriodicPressureGradient()
  INTEGER I,J,IEL
  REAL*8 ddx,ddy,ddz,ddp

  CALL B_Mul_U(qlMat%ColA,qlMAt%LdA,BXMat,BYMat,BZMat,LinSc%valP(NLMAX)%x,&
  QuadSc%defU,QuadSc%defV,QuadSc%defW,QuadSc%ndof,TSTEP,1d0)

  DO I=1,QuadSc%ndof
   QuadSc%defW(i) = QuadSc%defW(i) + TSTEP*dPeriodicVector(i)
  END DO

END SUBROUTINE AddPeriodicPressureGradient
!=========================================================================
!
!=========================================================================
SUBROUTINE AddGravForce
  INTEGER I,LINT
  REAL*8 DAUX
  EXTERNAL E013

  ILEV=NLMAX
  CALL SETLEV(2)

  DAUX=SQRT(Properties%Gravity(1)**2d0+Properties%Gravity(2)**2d0+Properties%Gravity(3)**2d0)

  IF (DAUX.GT.0d0) THEN
    CALL Grav_QuadSc(QuadSc%defU,QuadSc%defV,QuadSc%defW,mgDensity(ILEV)%x,&
      Properties%Gravity,QuadSc%ndof,&
      mg_mesh%level(ilev)%kvert,&
      mg_mesh%level(ilev)%karea,&
      mg_mesh%level(ilev)%kedge,&
      KWORK(L(KLINT(NLMAX))),&
      mg_mesh%level(ilev)%dcorvg,&
      tstep,E013)

  END IF

END SUBROUTINE AddGravForce
!=========================================================================
!
!=========================================================================
SUBROUTINE OperatorRegenaration(iType)
  INTEGER iType
  LOGICAL bHit

  bHit = .FALSE.

  IF (iType.EQ.myMatrixRenewal%D) THEN
    CALL Create_DiffMat(QuadSc)
    bHit = .TRUE.
  END IF

  IF (iType.EQ.myMatrixRenewal%K) THEN
    CALL Create_KMat(QuadSc)
    bHit = .TRUE.
  END IF

  IF (iType.EQ.myMatrixRenewal%M) THEN
    CALL Create_MRhoMat()
    bHit = .TRUE.
  END IF

  IF (iType.EQ.myMatrixRenewal%S) THEN
    CALL Create_SMat(QuadSc)
    bHit = .TRUE.
  END IF

  IF (iType.EQ.myMatrixRenewal%C) THEN

    CALL Create_BMat()

    IF (myid.ne.master) THEN
      CALL Fill_QuadLinParMat()
    END IF

    CALL SetSlipOnBandBT()

    CALL Create_CMat(QuadSc%knprU,QuadSc%knprV,QuadSc%knprW,LinSc%knprP,LinSc%prm%MGprmIn%MinLev,LinSc%prm%MGprmIn%CrsSolverType)
    IF (myid.ne.master) THEN
      CALL Create_ParCMat(QuadSc%knprU,QuadSc%knprV,QuadSc%knprW)
    END IF
    bHit = .TRUE.
  END IF

  IF (iType.EQ.1.and.bNS_Stabilization) then
   CALL Create_hDiffMat()
   bHit = .TRUE.
  END IF

  IF (iType.EQ.3.and.(NewtonForBurgers.ne.0d0)) then
   CALL Create_barMMat_iso(QuadSc)
   bHit = .TRUE.
  END IF

  IF (myid.EQ.ShowID.AND.bHit) WRITE(MTERM,'(A)', advance='yes') " "

END SUBROUTINE OperatorRegenaration
!=========================================================================
!
!=========================================================================
SUBROUTINE  AddViscoStress()
  EXTERNAL E013

  ILEV=NLMAX
  CALL SETLEV(2)

  CALL AssembleViscoStress(QuadSc%defU,QuadSc%defV,QuadSc%defW,BndrForce,&
    ViscoSc%Val11,ViscoSc%Val22,ViscoSc%Val33,ViscoSc%Val12,ViscoSc%Val13,ViscoSc%Val23,&
    mg_mesh%level(ilev)%kvert,&
    mg_mesh%level(ilev)%karea,&
    mg_mesh%level(ilev)%kedge,&
    mg_mesh%level(ilev)%dcorvg,&
    Properties%Viscosity(2),tstep,Properties%ViscoLambda,ViscoElasticForce,E013)

  CALL Comm_SummN(ViscoElasticForce,3)

END SUBROUTINE  AddViscoStress
