! TODO: new name Navier-Stokes + Q2 or similar
MODULE Transport_Q2P1  

USE def_QuadScalar
! USE PP3D_MPI
USE PP3D_MPI, ONLY:myid,master,E011Sum,COMM_Maximum,COMM_MaximumX,COMM_Minimum,&
                   COMM_NLComplete,Comm_Summ,Comm_SummN,&
                   myMPI_Barrier,coarse
USE Parametrization,ONLY : InitBoundaryStructure,ReviseWallBC,myParBndr,&
ParametrizeQ2Nodes

USE Sigma_User, ONLY: mySigma,myThermodyn,myProcess,mySetup,myMultiMat,BKTPRELEASE
! USE PP3D_MPI, ONLY:E011Sum,E011True_False,Comm_NLComplete,&
!               Comm_Maximum,Comm_Summ,knprmpi,myid,master
! USE LinScalar, ONLY: AddSurfaceTension
use fbm
use fbm_particle_reynolds, only: fbm_compute_particle_reynolds, fbm_compute_particle_reynolds_interface, &
                                 fbm_compute_particle_reynolds_interface_extended, &
                                 fbm_compute_particle_reynolds_farfield, &
                                 fbm_compute_particle_reynolds_reference_shell

use var_QuadScalar, only: QuadSc, LinSc, ViscoSc, PLinSc, Viscosity, &
                          bPrintParticleReynolds

use EL_CONFIG, only: el_apply_fluid_feedback, el_prescribed_field, el_shear_rate, &
                     el_prescribed_umax, el_cylinder_center, el_cylinder_radius, &
                     el_cylinder_axis, el_initial_field
use EL_CONFIG, only: el_write_diagnostics, el_fluid_gravity, el_meso_filter
use EL_GEOMETRY, only: EL_DOMAIN_Z_CENTER
use EL_DIAGNOSTICS, only: EL_WRITE_MOMENTUM_DIAGNOSTICS, &
                          EL_CAPTURE_MOMENTUM_REFERENCE, el_momentum_reference_set, &
                          EL_WRITE_SS_RADIUS, EL_FLUID_PAIR_STAGE
use EL_FIELDS, only: EL_APPLY_FLUID_FEEDBACK_SOURCE, el_field_data

use, intrinsic :: ieee_arithmetic

IMPLICIT NONE

REAL*8, ALLOCATABLE :: ST_force(:)
REAL*8 :: Density_Secondary=1d0,Density_Primary=1d0
REAL*8 :: myPowerLawFluid(3),ViscoElasticForce(3)
REAL*8 :: Sigma=0.034D0,DiracEps=0.00625d0
INTEGER, ALLOCATABLE :: QuadScBoundary(:)
INTEGER PressureSample(2)
REAL tttt0,tttt1

!!!!!!!!!!   Artificial - TimeStepControl !!!!!!!!!!!
REAL*8  :: xTimeStepMultiplier=1d0,old_TSTEP
INTEGER :: TimeStepIncrease = 0,TimeStepIncreaseCrit = 3, MaxSmootheningSteps = 32
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! interfaces for the fbm_update and fbm_geom function
! handlers that process the dynamics update and
! the geometric computations for fbm objects
include 'fbm_geom_include.h'
include 'fbm_up_include.h'
include 'fbm_vel_bc_include.h'

! The handler function for the dynamics update
procedure(update_fbm_handler), pointer :: fbm_up_handler_ptr => null()
procedure(fbm_force_wrapper), pointer :: fbm_force_handler_ptr => null()

! The handler function for the geometry update
procedure(fbm_geom_handler), pointer :: fbm_geom_handler_ptr => null()

! The handler function for the velocity boundary condition update
! for the fictitious boundary object
procedure(fbm_velBC_handler), pointer :: fbm_vel_bc_handler_ptr => null()

CONTAINS
include 'QuadSc_utilities.f90'
include 'QuadSc_corrections.f90'
include 'QuadSc_mesh_operations.f90'
include 'QuadSc_geometry_utilities.f90'
include 'QuadSc_material_properties.f90'
include 'QuadSc_integration.f90'
include 'QuadSc_force_torque_calc.f90'
#include "QuadSc_boundary.f90"
include 'QuadSc_initialization.f90'
!=========================================================================
! 
!=========================================================================
SUBROUTINE Transport_Q2P1_UxyzP(mfile,inl_u,itns)

use cinterface, only: calculateDynamics,calculateFBM
use fbm, only: fbm_updateFBM
INTEGER mfile,INL,inl_u,itns
REAL*8  ResU,ResV,ResW,DefUVW,RhsUVW,DefUVWCrit
REAL*8  ResP,DefP,RhsPG,defPG,defDivU,DefPCrit, global_lubrication
INTEGER INLComplete,I,J,IERR,iOuter,iITER
integer :: error_indicator
external E013

CALL updateFBMGeometry()
CALL report_and_reset_hashgrid_stats()

thstep = tstep*(1d0-theta)

CALL OperatorRegenaration(2)

CALL OperatorRegenaration(3)

! -------------------------------------------------
! Compute the momentum equations
! -------------------------------------------------

! GOTO 1
IF (myid.ne.master) THEN

 CALL ZTIME(tttt0)

 ! Assemble the right hand side
 CALL Matdef_General_QuadScalar(QuadSc,1)

 ! Add the pressure gradient to the rhs
 CALL AddPressureGradient()
END IF

 ! Add the viscoelastic stress to the rhs
 IF(bViscoElastic)THEN
   CALL AddViscoStress()
 END IF

IF (myid.ne.master) THEN
 ! Add the gravity force to the rhs
 CALL AddGravForce()
 CALL AddConstantForce()

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 ! Store the constant right hand side
 QuadSc%rhsU = QuadSc%defU
 QuadSc%rhsV = QuadSc%defV
 QuadSc%rhsW = QuadSc%defW

! Set dirichlet boundary conditions on the solution
 CALL Boundary_QuadScalar_Val()

END IF

thstep = tstep*theta

IF (myid.ne.master) THEN

 ! Assemble the defect vector and fine level matrix
 CALL Matdef_General_QuadScalar(QuadSc,-1)

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 QuadSc%auxU = QuadSc%defU
 QuadSc%auxV = QuadSc%defV
 QuadSc%auxW = QuadSc%defW
 CALL E013Sum3(QuadSc%auxU,QuadSc%auxV,QuadSc%auxW)

 ! Save the old solution
 CALL LCP1(QuadSc%valU,QuadSc%valU_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valV,QuadSc%valV_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valW,QuadSc%valW_old,QuadSc%ndof)

 ! Compute the norm of the defect
 CALL Resdfk_General_QuadScalar(QuadSc,ResU,ResV,ResW,DefUVW,RhsUVW)

END IF

CALL COMM_Maximum(RhsUVW)
DefUVWCrit=MAX(RhsUVW*QuadSc%prm%defCrit,QuadSc%prm%MinDef)

CALL Protocol_QuadScalar(mfile,QuadSc,0,&
     ResU,ResV,ResW,DefUVW,DefUVWCrit," Momentum equation ")

CALL ZTIME(tttt1)
myStat%tDefUVW = myStat%tDefUVW + (tttt1-tttt0)

DO INL=1,QuadSc%prm%NLmax
INLComplete = 0

! Calling the solver
CALL Solve_General_QuadScalar(QuadSc,Boundary_QuadScalar_Val,&
Boundary_QuadScalar_Mat,Boundary_QuadScalar_Mat_9,mfile)

!!!!          Checking the quality of the result           !!!!
!!!! ----------------------------------------------------- !!!!

CALL OperatorRegenaration(3)

IF (myid.ne.master) THEN
! Restore the constant right hand side
 CALL ZTIME(tttt0)
 QuadSc%defU = QuadSc%rhsU
 QuadSc%defV = QuadSc%rhsV
 QuadSc%defW = QuadSc%rhsW
END IF

IF (myid.ne.master) THEN

 ! Assemble the defect vector and fine level matrix
 CALL Matdef_General_QuadScalar(QuadSc,-1)

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 QuadSc%auxU = QuadSc%defU
 QuadSc%auxV = QuadSc%defV
 QuadSc%auxW = QuadSc%defW
 CALL E013Sum3(QuadSc%auxU,QuadSc%auxV,QuadSc%auxW)

 ! Save the old solution
 CALL LCP1(QuadSc%valU,QuadSc%valU_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valV,QuadSc%valV_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valW,QuadSc%valW_old,QuadSc%ndof)

 ! Compute the defect
 CALL Resdfk_General_QuadScalar(QuadSc,ResU,ResV,ResW,DefUVW,RhsUVW)

END IF

! Checking convergence rates against criterions
RhsUVW=DefUVW
CALL COMM_Maximum(RhsUVW)
CALL Protocol_QuadScalar(mfile,QuadSc,INL,&
     ResU,ResV,ResW,DefUVW,RhsUVW)
IF (ISNAN(RhsUVW)) stop

IF ((DefUVW.LE.DefUVWCrit).AND.&
    (INL.GE.QuadSc%prm%NLmin)) INLComplete = 1

CALL COMM_NLComplete(INLComplete)
CALL ZTIME(tttt1)
myStat%tDefUVW = myStat%tDefUVW + (tttt1-tttt0)

IF (INLComplete.eq.1) GOTO 1
!IF (timens.lt.tstep+1d-8) GOTO 1

END DO

1 CONTINUE

! return
myStat%iNonLin = myStat%iNonLin + INL
inl_u = INL

! -------------------------------------------------
! Compute the pressure correction
! -------------------------------------------------
IF (myid.ne.0) THEN

 CALL ZTIME(tttt0)
 ! Save the old solution
 LinSc%valP_old = LinSc%valP(NLMAX)%x
 LinSc%valP(NLMAX)%x = 0d0

 ! Assemble the right hand side (RHS=1/k B^T U~)
 CALL Matdef_General_LinScalar(LinSc,QuadSc,PLinSc,1)

!  ! Assemble the right hand side (RHS:=RHS-C*Q)
!  CALL Matdef_General_LinScalar(LinSc,QuadSc,PLinSc,2)

 ! Save the right hand side
 LinSc%rhsP(NLMAX)%x = LinSc%defP(NLMAX)%x

 CALL ZTIME(tttt1)
 myStat%tDefP = myStat%tDefP + (tttt1-tttt0)
END IF

! Calling the solver
CALL Solve_General_LinScalar(LinSc,PLinSc,QuadSc,Boundary_LinScalar_Mat,Boundary_LinScalar_Def,mfile)

CALL Protocol_LinScalar(mfile,LinSc," Pressure-Poisson equation")

2 CONTINUE

IF (myid.ne.0) THEN
 CALL ZTIME(tttt0)
 CALL Velocity_Correction()
 CALL Pressure_Correction()
 CALL ZTIME(tttt1)
 myStat%tCorrUVWP = myStat%tCorrUVWP + (tttt1-tttt0)
END IF

CALL QuadScP1toQ2(LinSc,QuadSc)

CALL FAC_GetForces(mfile)

CALL GetNonNewtViscosity()

#ifdef HAVE_PE
if (myid.eq. 1) then
#ifdef PE_SERIAL_MODE
  write(*,*)'fbm force (SERIAL PE mode)'
#else
  write(*,*)'fbm force (PARALLEL PE mode)'
#endif
end if
#endif 
! Calculate the forces
call fbm_updateForces(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                      LinSc%valP(NLMAX)%x,&
                      fbm_force_handler_ptr)

if (bPrintParticleReynolds) then
  call fbm_compute_particle_reynolds(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                                     Viscosity,Properties%Density(1),mfile, E013)
end if

! Step the particle simulation
call fbm_updateFBM(Properties%Density(1),tstep,timens,&
                   Properties%Gravity,mfile,myid,&
                   QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                   LinSc%valP(NLMAX)%x,&
                   fbm_up_handler_ptr) 

IF (myid.ne.0) THEN
 CALL STORE_OLD_MESH(mg_mesh%level(NLMAX+1)%dcorvg)
END IF
 
 CALL UmbrellaSmoother_ext(0d0,nUmbrellaSteps)
 
IF (myid.ne.0) THEN
 CALL STORE_NEW_MESH(mg_mesh%level(NLMAX+1)%dcorvg)
END IF
 
 CALL GET_MESH_VELO()
 
 ILEV=NLMAX
 CALL SETLEV(2)
 CALL SetUp_myQ2Coor( mg_mesh%level(ILEV)%dcorvg,&
                      mg_mesh%level(ILEV)%dcorag,&
                      mg_mesh%level(ILEV)%kvert,&
                      mg_mesh%level(ILEV)%karea,&
                      mg_mesh%level(ILEV)%kedge)

 CALL updateFBMGeometry()
 CALL report_and_reset_hashgrid_stats()

RETURN

END SUBROUTINE Transport_q2p1_UxyzP
!========================================================================================
!
!========================================================================================
SUBROUTINE Transport_q2p1_UxyzP_el(mfile,inl_u,itns)

  USE EL_TRANSFER, ONLY: EL_PARTICLE_MESH_PASS, EL_ADVANCE_PARTICLES, &
                         el_pair_expected_impulse
  USE EL_DIAGNOSTICS, ONLY: EL_NEWTON_PAIR_BEGIN, EL_NEWTON_PAIR_END, &
                            EL_FLUID_PAIR_BEGIN, EL_FLUID_PAIR_END, &
                            EL_MESO_FILTER_APPLY
  USE var_QuadScalar, ONLY: bConstForce, ConstForce, myQ2Coor

  INTEGER, INTENT(IN) :: mfile, itns
  INTEGER, INTENT(OUT) :: inl_u
  REAL*8 :: dummy_velocity(1), el_ext_force(3)
  LOGICAL :: prescribed_active

  dummy_velocity = 0.0d0
  inl_u = 0
  prescribed_active = TRIM(el_prescribed_field).NE.'none'

  IF (ALLOCATED(FictKNPR)) FictKNPR = 0
  IF (myid.NE.master .AND. prescribed_active) CALL EL_IMPOSE_PRESCRIBED_FIELD()

  ! Viscosity convention: Properties%Viscosity is KINEMATIC (param_parser
  ! stores Prop@Viscosity as-is and converts Prop@DynVisc via /Density);
  ! the EL closures take the kinematic value and form mu = rho*nu
  ! internally (el_forces.f90 dynamic_viscosity). Do NOT pre-multiply by
  ! density here -- that double-counts rho for any fluid with rho /= 1
  ! (caught by the ten Cate SI validation case).
  ! Evaluate hydrodynamic particle forces from the current fluid state.
  IF (myid.EQ.master) THEN
    CALL EL_PARTICLE_MESH_PASS(mg_mesh, NLMAX, dummy_velocity, &
      dummy_velocity, dummy_velocity, dummy_velocity, Properties%Density(1), &
      Properties%Viscosity(1), Properties%Gravity, tstep, mfile, -itns)
  ELSE
    CALL EL_PARTICLE_MESH_PASS(mg_mesh, NLMAX, QuadSc%valU, QuadSc%valV, &
      QuadSc%valW, LinSc%valP(NLMAX)%x, Properties%Density(1), Properties%Viscosity(1), &
      Properties%Gravity, tstep, mfile, -itns)
  END IF

  ! Capture the element-integrated total-momentum reference at the true initial
  ! state, before the first particle or fluid update. The pre-advance EL pass
  ! above initializes epsilon_f but does not advance particle or fluid state.
  IF (myid.NE.master .AND. .NOT.prescribed_active .AND. &
      .NOT.el_momentum_reference_set) THEN
    CALL EL_CAPTURE_MOMENTUM_REFERENCE(QuadSc%valU, QuadSc%valV, &
      QuadSc%valW, mg_mesh, NLMAX, Properties%Density(1), &
      el_field_data%epsilon_f)
  END IF

  IF (myid.NE.master) CALL EL_NEWTON_PAIR_BEGIN()
  CALL EL_ADVANCE_PARTICLES()
  IF (myid.NE.master) CALL EL_NEWTON_PAIR_END(el_pair_expected_impulse, &
    timens, mfile, itns)

  ! Refresh volume fraction and diagnostic feedback after particle motion.
  IF (myid.EQ.master) THEN
    CALL EL_PARTICLE_MESH_PASS(mg_mesh, NLMAX, dummy_velocity, &
      dummy_velocity, dummy_velocity, dummy_velocity, Properties%Density(1), &
      Properties%Viscosity(1), Properties%Gravity, tstep, mfile, itns)
  ELSE
    CALL EL_PARTICLE_MESH_PASS(mg_mesh, NLMAX, QuadSc%valU, QuadSc%valV, &
      QuadSc%valW, LinSc%valP(NLMAX)%x, Properties%Density(1), Properties%Viscosity(1), &
      Properties%Gravity, tstep, mfile, itns)
  END IF

  ! Radial-migration monitor; internally gated on cylinder domain + audit
  ! frequency, and (unlike the momentum audit) also active in frozen runs.
  IF (myid.NE.master) CALL EL_WRITE_SS_RADIUS(timens, mfile, itns)

  IF (ALLOCATED(FictKNPR)) FictKNPR = 0
  IF (.NOT.prescribed_active) THEN
    ! Fluid-side momentum audit around the whole fluid solve: body force
    ! per unit volume = rho*(ConstForce + Gravity if fluid gravity is on).
    ! EL_INTEGRATE_FLUID_MOMENTUM relies on the GLOBAL multigrid level
    ! (/MGPAR/ ILEV) for its NDFGL structures; the preceding EL mesh pass
    ! restores whatever coarse level the previous solver phase used, so pin
    ! ILEV to the coupling level before each audit call (the post-solve END
    ! call is pinned too, defensively).
    el_ext_force = 0.0d0
    IF (bConstForce) el_ext_force = el_ext_force + &
      Properties%Density(1)*ConstForce
    IF (el_fluid_gravity) el_ext_force = el_ext_force + &
      Properties%Density(1)*Properties%Gravity
    IF (myid.NE.master .AND. el_apply_fluid_feedback .AND. &
        ALLOCATED(el_field_data%fluid_feedback_source)) THEN
      ILEV = NLMAX
      CALL SETLEV(2)
      CALL EL_FLUID_PAIR_BEGIN(QuadSc%valU, QuadSc%valV, &
        QuadSc%valW, mg_mesh, NLMAX, Properties%Density(1), &
        el_field_data%epsilon_f)
    END IF
    CALL Transport_q2p1_UxyzP_fluid_core(mfile,inl_u,itns,.FALSE.)
    ! Mesoscale lane-mode filter (periodic settling boxes): applied inside
    ! the audited window so its (near-zero) momentum residue lands in the
    ! EL_FLUID_PAIR mismatch and is nulled by the ELMomentumFix deadbeat.
    IF (myid.NE.master .AND. el_meso_filter) THEN
      ILEV = NLMAX
      CALL SETLEV(2)
      CALL SetUp_myQ2Coor(mg_mesh%level(ILEV)%dcorvg, &
                          mg_mesh%level(ILEV)%dcorag, &
                          mg_mesh%level(ILEV)%kvert, &
                          mg_mesh%level(ILEV)%karea, &
                          mg_mesh%level(ILEV)%kedge)
      CALL EL_MESO_FILTER_APPLY(QuadSc%valW, mg_mesh, NLMAX, myQ2Coor, &
        QuadSc%ndof)
    END IF
    IF (myid.NE.master .AND. el_apply_fluid_feedback .AND. &
        ALLOCATED(el_field_data%fluid_feedback_source)) THEN
      ILEV = NLMAX
      CALL SETLEV(2)
      CALL EL_FLUID_PAIR_END(QuadSc%valU, QuadSc%valV, QuadSc%valW, &
        mg_mesh, NLMAX, Properties%Density(1), el_field_data%epsilon_f, &
        el_field_data%fluid_feedback_source, el_ext_force, tstep, timens, &
        mfile, itns)
    END IF
  END IF

  IF (myid.NE.master) THEN
    CALL EL_WRITE_MOMENTUM_DIAGNOSTICS(QuadSc%valU, QuadSc%valV, &
      QuadSc%valW, timens, mfile, itns, mg_mesh, NLMAX, &
      Properties%Density(1), Properties%Viscosity(1), el_field_data%epsilon_f)
  END IF

END SUBROUTINE Transport_q2p1_UxyzP_el

SUBROUTINE EL_IMPOSE_PRESCRIBED_FIELD()

  INTEGER :: i, ndof
  REAL*8 :: zc, r2, radius2

  IF (TRIM(el_prescribed_field).EQ.'none') RETURN

  ILEV = NLMAX
  CALL SETLEV(2)
  CALL SetUp_myQ2Coor(mg_mesh%level(ILEV)%dcorvg, &
                      mg_mesh%level(ILEV)%dcorag, &
                      mg_mesh%level(ILEV)%kvert, &
                      mg_mesh%level(ILEV)%karea, &
                      mg_mesh%level(ILEV)%kedge)

  ndof = mg_mesh%level(ILEV)%nvt + mg_mesh%level(ILEV)%net + &
         mg_mesh%level(ILEV)%nat + mg_mesh%level(ILEV)%nel

  SELECT CASE (TRIM(el_prescribed_field))
  CASE ('linear_shear')
    zc = EL_DOMAIN_Z_CENTER()
    DO i=1,ndof
      QuadSc%valU(i) = el_shear_rate*(myQ2Coor(3,i)-zc)
      QuadSc%valV(i) = 0.0d0
      QuadSc%valW(i) = 0.0d0
    END DO
  CASE ('poiseuille')
    radius2 = el_cylinder_radius*el_cylinder_radius
    IF (radius2.LE.0.0d0) RETURN
    DO i=1,ndof
      QuadSc%valU(i) = 0.0d0
      QuadSc%valV(i) = 0.0d0
      QuadSc%valW(i) = 0.0d0
      SELECT CASE (TRIM(el_cylinder_axis))
      CASE ('x')
        r2 = (myQ2Coor(2,i)-el_cylinder_center(2))**2 + &
             (myQ2Coor(3,i)-el_cylinder_center(3))**2
        QuadSc%valU(i) = el_prescribed_umax*MAX(0.0d0,1.0d0-r2/radius2)
      CASE ('y')
        r2 = (myQ2Coor(1,i)-el_cylinder_center(1))**2 + &
             (myQ2Coor(3,i)-el_cylinder_center(3))**2
        QuadSc%valV(i) = el_prescribed_umax*MAX(0.0d0,1.0d0-r2/radius2)
      CASE DEFAULT
        r2 = (myQ2Coor(1,i)-el_cylinder_center(1))**2 + &
             (myQ2Coor(2,i)-el_cylinder_center(2))**2
        QuadSc%valW(i) = el_prescribed_umax*MAX(0.0d0,1.0d0-r2/radius2)
      END SELECT
    END DO
  END SELECT

END SUBROUTINE EL_IMPOSE_PRESCRIBED_FIELD
!========================================================================================
!
!========================================================================================
! One-shot initial condition for COUPLED runs (ELInitialField=linear_shear):
! seeds the exact single-phase Couette solution u_x = G*(z-z_c) so the run
! skips the viscous startup (L^2/nu time units). Called once from the app
! driver after init (istart=0 only); afterwards the fluid solve owns the
! field, unlike EL_IMPOSE_PRESCRIBED_FIELD which re-freezes it every step.
SUBROUTINE EL_IMPOSE_INITIAL_FIELD()

  INTEGER :: i, ndof
  REAL*8 :: zc

  IF (TRIM(el_initial_field).NE.'linear_shear') RETURN
  IF (myid.EQ.master) RETURN

  ILEV = NLMAX
  CALL SETLEV(2)
  CALL SetUp_myQ2Coor(mg_mesh%level(ILEV)%dcorvg, &
                      mg_mesh%level(ILEV)%dcorag, &
                      mg_mesh%level(ILEV)%kvert, &
                      mg_mesh%level(ILEV)%karea, &
                      mg_mesh%level(ILEV)%kedge)

  ndof = mg_mesh%level(ILEV)%nvt + mg_mesh%level(ILEV)%net + &
         mg_mesh%level(ILEV)%nat + mg_mesh%level(ILEV)%nel

  zc = EL_DOMAIN_Z_CENTER()
  DO i=1,ndof
    QuadSc%valU(i) = el_shear_rate*(myQ2Coor(3,i)-zc)
    QuadSc%valV(i) = 0.0d0
    QuadSc%valW(i) = 0.0d0
  END DO

  ! Re-assert Dirichlet values so wall DOFs carry the BC, not the raw IC.
  CALL Boundary_QuadScalar_Val()

END SUBROUTINE EL_IMPOSE_INITIAL_FIELD
!========================================================================================
!
!========================================================================================
SUBROUTINE Transport_q2p1_UxyzP_fc_ext(mfile,inl_u,itns)

INTEGER, INTENT(IN) :: mfile, itns
INTEGER, INTENT(OUT) :: inl_u

CALL Transport_q2p1_UxyzP_fluid_core(mfile,inl_u,itns,.TRUE.)

END SUBROUTINE Transport_q2p1_UxyzP_fc_ext
!========================================================================================
!
! Shared low-level Q2/P1 solve. enable_fbm is true only for the validated FBM transport.
!========================================================================================
SUBROUTINE Transport_q2p1_UxyzP_fluid_core(mfile,inl_u,itns,enable_fbm)
use cinterface, only: calculateDynamics,calculateFBM
use fbm, only: fbm_updateFBM, fbm_velBCTest,fbm_testFBMGeom
use PP3D_MPI, only: Barrier_myMPI, Sum_myMPI
#ifdef HAVE_PE
use dem_query, only: numLocalParticles, numTotalParticles, &
                     get_lubrication_enabled, get_lubrication_stage_diag, &
                     isSphere, getParticleOrientation
#endif
external E013

INTEGER, INTENT(IN) :: mfile,itns
INTEGER, INTENT(OUT) :: inl_u
LOGICAL, INTENT(IN) :: enable_fbm
INTEGER INL
REAL*8  ResU,ResV,ResW,DefUVW,RhsUVW,DefUVWCrit
REAL*8  ResP,DefP,RhsPG,defPG,defDivU,DefPCrit, global_lubrication
INTEGER INLComplete,I,J,IERR,iITER
INTEGER iaux_lev
real*8 px, py, pz
real*8 dres_buf(2), dofs_pp, d_over_h
real*8 dLubF, dLubDiss, dLubJmax
integer nLubPairs, nLubSat
real*8 part_axis(3)
integer ipart_ax
integer k
k=1

IF (enable_fbm) THEN
  CALL updateFBMGeometry()
  CALL report_and_reset_hashgrid_stats()
END IF

thstep = tstep*(1d0-theta)

CALL OperatorRegenaration(2)

CALL OperatorRegenaration(3)

! -------------------------------------------------
! Compute the momentum equations
! -------------------------------------------------
! GOTO 1
IF (myid.ne.master) THEN

 CALL ZTIME(tttt0)

 ! Assemble the right hand side
 CALL Matdef_General_QuadScalar(QuadSc,1)

! Add the pressure gradient
  CALL AddPressureGradient()

! Add the pressure gradient with the jump term to the rhs
!   CALL AddPressureGradientWithJump()
END IF

 ! Add the viscoelastic stress to the rhs
 IF(bViscoElastic)THEN
   CALL AddViscoStress()
 END IF

IF (myid.ne.master) THEN

 ! Add the gravity force to the rhs
 CALL AddGravForce()
 CALL AddConstantForce()
 CALL EL_APPLY_MOMENTUM_FIX()
 IF (el_apply_fluid_feedback) THEN
   CALL EL_APPLY_FLUID_FEEDBACK_SOURCE(QuadSc%defU, QuadSc%defV, &
     QuadSc%defW, tstep)
 END IF

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 ! Store the constant right hand side
 QuadSc%rhsU = QuadSc%defU
 QuadSc%rhsV = QuadSc%defV
 QuadSc%rhsW = QuadSc%defW

 ! Set dirichlet boundary conditions on the solution
 CALL Boundary_QuadScalar_Val()

END IF

thstep = tstep*theta

IF (myid.ne.master) THEN

 ! Assemble the defect vector and fine level matrix
 CALL Matdef_General_QuadScalar(QuadSc,-1)

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 QuadSc%auxU = QuadSc%defU
 QuadSc%auxV = QuadSc%defV
 QuadSc%auxW = QuadSc%defW
 CALL E013Sum3(QuadSc%auxU,QuadSc%auxV,QuadSc%auxW)
!  CALL E013Sum(QuadSc%auxU)
!  CALL E013Sum(QuadSc%auxV)
!  CALL E013Sum(QuadSc%auxW)

 ! Save the old solution
 CALL LCP1(QuadSc%valU,QuadSc%valU_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valV,QuadSc%valV_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valW,QuadSc%valW_old,QuadSc%ndof)

 ! Compute the norm of the defect
 CALL Resdfk_General_QuadScalar(QuadSc,ResU,ResV,ResW,DefUVW,RhsUVW)

END IF

CALL COMM_MaximumX(RhsUVW)
DefUVWCrit=MAX(RhsUVW*QuadSc%prm%defCrit,QuadSc%prm%MinDef)

CALL Protocol_QuadScalar(mfile,QuadSc,0,&
     ResU,ResV,ResW,DefUVW,DefUVWCrit," Momentum equation ")

CALL ZTIME(tttt1)
myStat%tDefUVW = myStat%tDefUVW + (tttt1-tttt0)

! CALL ALStructExtractor()

DO INL=1,QuadSc%prm%NLmax
INLComplete = 0

! Calling the solver
CALL Solve_General_QuadScalar(QuadSc,Boundary_QuadScalar_Val,&
Boundary_QuadScalar_Mat,Boundary_QuadScalar_Mat_9,mfile)

!!!!          Checking the quality of the result           !!!!
!!!! ----------------------------------------------------- !!!!

CALL OperatorRegenaration(3)

IF (myid.ne.master) THEN
! Restore the constant right hand side
 CALL ZTIME(tttt0)
 QuadSc%defU = QuadSc%rhsU
 QuadSc%defV = QuadSc%rhsV
 QuadSc%defW = QuadSc%rhsW
END IF

IF (myid.ne.master) THEN

 ! Assemble the defect vector and fine level matrix
 CALL Matdef_General_QuadScalar(QuadSc,-1)

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_QuadScalar_Def()

 QuadSc%auxU = QuadSc%defU
 QuadSc%auxV = QuadSc%defV
 QuadSc%auxW = QuadSc%defW
 CALL E013Sum3(QuadSc%auxU,QuadSc%auxV,QuadSc%auxW)
! CALL E013Sum(QuadSc%auxU)
! CALL E013Sum(QuadSc%auxV)
! CALL E013Sum(QuadSc%auxW)

 ! Save the old solution
 CALL LCP1(QuadSc%valU,QuadSc%valU_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valV,QuadSc%valV_old,QuadSc%ndof)
 CALL LCP1(QuadSc%valW,QuadSc%valW_old,QuadSc%ndof)

 ! Compute the defect
 CALL Resdfk_General_QuadScalar(QuadSc,ResU,ResV,ResW,DefUVW,RhsUVW)

END IF

! Checking convergence rates against criterions
RhsUVW=DefUVW
CALL COMM_MaximumX(RhsUVW)
CALL Protocol_QuadScalar(mfile,QuadSc,INL,&
     ResU,ResV,ResW,DefUVW,RhsUVW)
IF (ISNAN(RhsUVW)) stop

IF ((DefUVW.LE.DefUVWCrit).AND.&
    (INL.GE.QuadSc%prm%NLmin)) INLComplete = 1

CALL COMM_NLComplete(INLComplete)
CALL ZTIME(tttt1)
myStat%tDefUVW = myStat%tDefUVW + (tttt1-tttt0)

IF (INLComplete.eq.1) GOTO 1
!IF (timens.lt.tstep+1d-8) GOTO 1

END DO

1 CONTINUE

! return
myStat%iNonLin = myStat%iNonLin + INL
inl_u = INL

! Mid-solve momentum snapshot for the EL fluid audit: after the momentum
! (Burgers) solve, before the projection. EL_INTEGRATE_FLUID_MOMENTUM
! needs the global ILEV pinned to the coupling level for its NDFGL
! lookups; restore the solver's level afterwards.
IF (myid.NE.master .AND. el_apply_fluid_feedback .AND. &
    ALLOCATED(el_field_data%fluid_feedback_source)) THEN
  iaux_lev = ILEV
  ILEV = NLMAX
  CALL SETLEV(2)
  CALL EL_FLUID_PAIR_STAGE(QuadSc%valU, QuadSc%valV, QuadSc%valW, &
    mg_mesh, NLMAX, Properties%Density(1), el_field_data%epsilon_f)
  ILEV = iaux_lev
  CALL SETLEV(2)
END IF

! -------------------------------------------------
! Compute the pressure correction
! -------------------------------------------------
CALL MonitorVeloMag(QuadSc)

IF (myid.ne.0) THEN

 CALL ZTIME(tttt0)
 ! Save the old solution
 LinSc%valP_old = LinSc%valP(NLMAX)%x
 LinSc%valP(NLMAX)%x = 0d0

 ! Assemble the right hand side (RHS=1/k B^T U~)
 CALL Matdef_General_LinScalar(LinSc,QuadSc,PLinSc,1)

!  ! Assemble the right hand side (RHS:=RHS-C*Q)
!  CALL Matdef_General_LinScalar(LinSc,QuadSc,PLinSc,2)

 ! Save the right hand side
 LinSc%rhsP(NLMAX)%x = LinSc%defP(NLMAX)%x

 CALL ZTIME(tttt1)
 myStat%tDefP = myStat%tDefP + (tttt1-tttt0)
END IF

! Calling the solver
CALL Solve_General_LinScalar(LinSc,PLinSc,QuadSc,Boundary_LinScalar_Mat,Boundary_LinScalar_Def,mfile)

CALL Protocol_LinScalar(mfile,LinSc," Pressure-Poisson equation")

2 CONTINUE

IF (myid.ne.0) THEN
 CALL ZTIME(tttt0)
 !if (myid.eq.1) write(*,*) 'no correction ... '
 CALL Velocity_Correction()
 CALL Pressure_Correction()
 CALL ZTIME(tttt1)
 myStat%tCorrUVWP = myStat%tCorrUVWP + (tttt1-tttt0)
END IF

CALL QuadScP1toQ2(LinSc,QuadSc)

CALL FAC_GetForces(mfile)
CALL FAC_GetSurfForces(mfile)

!CALL VISC_GetTorqueVolume(mfile)
!CALL VISC_GetTorqueVolume(mfile)

CALL GetNonNewtViscosity()

IF (bNS_Stabilization) THEN
 CALL ExtractVeloGradients()
END IF

#ifdef HAVE_PE
IF (enable_fbm) THEN
  if (myid.eq. 1) write(*,*)'fbm force'
END IF
#endif

total_lubrication = 0.0d0

! Calculate the forces
if(enable_fbm .and. .not. skipFBMForce) then
  call fbm_updateForces(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                        LinSc%valP(NLMAX)%x,&
                        fbm_force_handler_ptr)
end if


if (enable_fbm .and. bPrintParticleReynolds) then
  call fbm_compute_particle_reynolds_interface_extended(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                                                        FictKNPR,Viscosity,Properties%Density(1),mfile, E013, 2)

  call fbm_compute_particle_reynolds_farfield(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                                              FictKNPR,Viscosity,Properties%Density(1),mfile,E013)

  ! Keep the robust shell-based reference velocity as the default ParticleRe value.
  call fbm_compute_particle_reynolds_reference_shell(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                                                     FictKNPR,Viscosity,Properties%Density(1),mfile,E013)
end if

!call Sum_myMPI(total_lubrication, global_lubrication)
!call Get_SoosWallShear(mfile)
!call Get_DissipationIntegral(mfile)


#ifdef HAVE_PE
IF (enable_fbm) THEN
  if (myid.eq. 1) write(*,*)'fbm update'
END IF
#endif

if(enable_fbm .and. .not. skipFBMDynamics) then
  ! Step the particle simulation
  call fbm_updateFBM(Properties%Density(1),tstep,timens,&
                     Properties%Gravity,mfile,myid,&
                     QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                     LinSc%valP(NLMAX)%x,&
                     fbm_up_handler_ptr) 
end if

!IF (myid.ne.0) THEN
! CALL STORE_OLD_MESH(mg_mesh%level(NLMAX+1)%dcorvg)
!END IF
! 
!!if (myid.eq. 1) write(*,*)'umbrella smoother'
!!CALL UmbrellaSmoother_ext(0d0,nUmbrellaSteps)
! 
!IF (myid.ne.0) THEN
! CALL STORE_NEW_MESH(mg_mesh%level(NLMAX+1)%dcorvg)
!END IF
! 
! CALL GET_MESH_VELO()
! 
! ILEV=NLMAX
! CALL SETLEV(2)
! CALL SetUp_myQ2Coor( mg_mesh%level(ILEV)%dcorvg,&
!                      mg_mesh%level(ILEV)%dcorag,&
!                      mg_mesh%level(ILEV)%kvert,&
!                      mg_mesh%level(ILEV)%karea,&
!                      mg_mesh%level(ILEV)%kedge)

!CALL updateFBMGeometry()

CALL MonitorVeloMag(QuadSc)
CALL ComputeCFL(QuadSc)

IF (myid .EQ. 1) THEN
  IF (bPrintCFL) THEN
    WRITE(*,'(A,F10.4)') ' CFL (fluid): ', cfl_global
  END IF
  IF (bPrintParticleCFL) THEN
    WRITE(*,'(A,F10.4)') ' CFL (particle): ', cfl_particle_global
    WRITE(*,'(A,ES12.4,A,ES12.4,A,ES12.4)') &
      '   CFL components: vp_max=', vp_max_global, &
      '  h_min=', h_min_global, '  dt=', tstep
  END IF
END IF

#ifdef HAVE_PE
IF (bPrintParticleState .AND. enable_fbm) THEN
  ! Rank-symmetric reduction of the per-rank inside-DOF counts stored by
  ! QuadScalar_FictKnpr (the classification block itself runs under
  ! rank-asymmetric control flow and must stay collective-free). The flag
  ! is parsed identically on every rank, so gating the collective is safe.
  dres_buf(1) = DBLE(dns_inside_dofs_local)
  dres_buf(2) = 0d0
#ifdef PE_SERIAL_MODE
  ! Serial PE: every worker holds the full particle world - count it once.
  IF (myid .EQ. 1) dres_buf(2) = DBLE(numTotalParticles())
#else
  ! Parallel PE: each particle is local to exactly one rank - sum is global.
  IF (myid .NE. 0) dres_buf(2) = DBLE(numLocalParticles())
#endif
  CALL Comm_SummN(dres_buf, 2)
  IF (myid .EQ. 1) THEN
    dofs_pp = dres_buf(1) / MAX(dres_buf(2), 1d0)
    d_over_h = 0d0
    IF (h_min_global .GT. 0d0) d_over_h = 2d0*particle_rad_global/h_min_global
    WRITE(*,'(A,ES16.8,A,F12.1,A,I9,A,ES16.8,A,F10.3)') &
      'DNS_RESOLUTION time= ', timens, ' dofs_per_particle= ', dofs_pp, &
      ' nparticles= ', NINT(dres_buf(2)), ' h_min= ', h_min_global, &
      ' D_over_h= ', d_over_h
    WRITE(*,'(A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8)') &
      'DNS_CFL time= ', timens, ' cfl_fluid= ', cfl_global, &
      ' cfl_particle= ', cfl_particle_global, ' vp_max= ', vp_max_global, &
      ' dt= ', tstep
#ifdef PE_SERIAL_MODE
    ! D6.1: orientation record for non-spherical bodies (world-frame
    ! direction of the body a-axis). Sphere-only runs print nothing here.
    ! Serial PE only: rank 1 holds the full particle world.
    DO ipart_ax = 0, NINT(dres_buf(2))-1
      IF (.NOT. isSphere(ipart_ax)) THEN
        CALL getParticleOrientation(ipart_ax, part_axis)
        WRITE(*,'(A,ES16.8,A,I6,A,3ES17.8)') &
          'DNS_PART_AXIS time= ', timens, ' ip= ', ipart_ax, &
          ' axis= ', part_axis
      END IF
    END DO
#endif
    ! Lubrication-stage record (D2.2 G2 instrumentation): rank-1 PE instance;
    ! in PE serial mode every rank holds the full world, so this is the global
    ! record. For a single sphere-wall pair F_total IS the applied lubrication
    ! force and J_max/dt its normal component. Printed only while the add-on
    ! is enabled, so lubrication-off runs keep byte-identical stdout.
    IF (get_lubrication_enabled() .NE. 0) THEN
      CALL get_lubrication_stage_diag(dLubF, dLubDiss, dLubJmax, nLubPairs, nLubSat)
      WRITE(*,'(A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A,I6,A,I6)') &
        'DNS_LUB time= ', timens, ' F_total= ', dLubF, &
        ' J_max= ', dLubJmax, ' dissipation= ', dLubDiss, &
        ' n_pairs= ', nLubPairs, ' n_saturated= ', nLubSat
    END IF
  END IF
END IF
#endif

RETURN

END SUBROUTINE Transport_q2p1_UxyzP_fluid_core
!=========================================================================
! 
!=========================================================================
! Include custom implementations of the Q2 transport equation
include 'QuadSc_transport_extensions.f90'
!=========================================================================
! 
!=========================================================================
include 'QuadSc_handlers.f90'
!=========================================================================
! 
!=========================================================================
subroutine InitMeshDeform(mfile, mgMesh)
USE var_QuadScalar, ONLY : nMainUmbrellaSteps,tMultiMesh
USE var_QuadScalar, ONLY : nUmbrellaStepsLvl
use umbrella_smoother, only : us_UmbrellaSmoother
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
USE def_FEAT

implicit none

integer, intent(in) :: mfile

type(tMultiMesh), intent(inout) :: mgMesh

! local variables
real :: dttt0,dttt1
integer :: i

CALL myMPI_Barrier()
CALL ztime(dttt0)

ilev=nlmax
call setlev(2)

if (nMainUmbrellaSteps.ne.0) then
 nMainUmbrellaSteps =nMainUmbrellaSteps + (mySetup%MeshResolution-1)
end if

do i=1,nMainUmbrellaSteps

 call us_UmbrellaSmoother(0d0, nUmbrellaStepsLvl, mgMesh, QuadSc)

 ilev=nlmax
 call setlev(2)

 call SetUp_myQ2Coor( mgMesh%level(ilev)%dcorvg,&
                      mgMesh%level(ilev)%dcorag,&
                      mgMesh%level(ilev)%kvert,&
                      mgMesh%level(ilev)%karea,&
                      mgMesh%level(ilev)%kedge)

END DO

call myMPI_Barrier()
call ztime(dttt1)
if (myid.eq.1) write(mfile,"(A,F6.1,A)") "Time used for mesh smoothening was: ", dttt1-dttt0, "s!"
if (myid.eq.1) write(*,"(A,F6.1,A)") "Time used for mesh smoothening was: ", dttt1-dttt0, "s!"

end subroutine InitMeshDeform
!=========================================================================
! 
!=========================================================================
subroutine InitOperators(mfile, mgMesh,bCreate)
use PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
use var_QuadScalar, only : tMultiMesh,bUseDumpedMixerGeometry,bCGALGeometryInitialized
implicit none

integer, intent(in) :: mfile
logical :: bCreate
type(tMultiMesh), intent(inout) :: mgMesh

! local variables
integer :: i

ilev = mgMesh%nlmax
call setlev(2)

call SetUp_myQ2Coor( mgMesh%level(ilev)%dcorvg,&
                     mgMesh%level(ilev)%dcorag,&
                     mgMesh%level(ilev)%kvert,&
                     mgMesh%level(ilev)%karea,&
                     mgMesh%level(ilev)%kedge)

call StoreOrigCoor(mgMesh%level(mgMesh%nlmax)%dcorvg)
call store_old_mesh(mgMesh%level(mgMesh%nlmax)%dcorvg)

ilev = mgMesh%nlmax
call setlev(2)

if (myid.ne.0 .and. .not.bUseDumpedMixerGeometry) call updateMixerGeometry(mfile)

! call FilterColdElements(mfile)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       PRESS BC        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (mySetup%bPressureFBM) THEN
 ilev=nlmin
 CALL SETLEV(2)
 CALL SetPressBC_NewGen(mgMesh)
 ! send them to the master
 ilev=nlmin
 CALL SETLEV(2)
 CALL SendPressBCElemsToCoarse(LinSc%knprP(ilev)%x,nel)
 
 ! propagate the structure consistently to the upper level
 do ilev=nlmin+1,nlmax
  CALL SETLEV(2)
  CALL GetMG_KNPRP(mgMesh)
 end do
! 
 ! Set up the boundary condition types (knpr)
 DO ILEV=NLMIN,NLMAX
  CALL SETLEV(2)
  CALL IncludeFBM_BCs(mgMesh)
  CALL QuadScalar_Knpr()
 END DO
END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       PRESS BC        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (bCreate) THEN
 CALL InitializeProlRest(QuadSc,LinSc)
 
 CALL MemoryPrint(1,'w','CGALOUT0')
#ifdef USE_CGAL
      if (bCGALGeometryInitialized) then
       CALL Release_cgal_structures()
       bCGALGeometryInitialized = .false.
      end if
#endif
 CALL MemoryPrint(1,'w','CGALOUT1')
 
 !!! for the SSE app it is assumed to have a constant density distribution, which depends !!!!
 !!! only on the local tempertaure and material distribution                              !!!!
 CALL UpdateDensityDistribution_XSE(mfile)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 call OperatorRegenaration(1)
 call OperatorRegenaration(2)
 call OperatorRegenaration(3)

 CALL SetUp_HYPRE_Solver(LinSc,PLinSc,mfile)
 
 CALL Create_MMat()
END IF

end subroutine InitOperators
!=========================================================================
! 
!=========================================================================
END MODULE Transport_Q2P1
