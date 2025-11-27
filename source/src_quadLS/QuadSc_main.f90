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
                                 fbm_compute_particle_reynolds_farfield

use var_QuadScalar, only: QuadSc, LinSc, ViscoSc, PLinSc, Viscosity

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

call fbm_compute_particle_reynolds(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                                   Viscosity,Properties%Density(1),mfile, E013)

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

RETURN

END SUBROUTINE Transport_q2p1_UxyzP
!========================================================================================
!
!========================================================================================
SUBROUTINE Transport_q2p1_UxyzP_fc_ext(mfile,inl_u,itns)
use cinterface, only: calculateDynamics,calculateFBM
use fbm, only: fbm_updateFBM, fbm_velBCTest,fbm_testFBMGeom
use PP3D_MPI, only: Barrier_myMPI, Sum_myMPI
external E013

INTEGER mfile,INL,inl_u,itns
REAL*8  ResU,ResV,ResW,DefUVW,RhsUVW,DefUVWCrit
REAL*8  ResP,DefP,RhsPG,defPG,defDivU,DefPCrit, global_lubrication
INTEGER INLComplete,I,J,IERR,iITER
real*8 px, py, pz
integer k
k=1

CALL updateFBMGeometry()

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

!CALL DNA_GetTorques(mfile)
!CALL DNA_GetTorques(mfile)

CALL GetNonNewtViscosity()

IF (bNS_Stabilization) THEN
 CALL ExtractVeloGradients()
END IF


!call fbm_testBasicFBM(0.0d0, 0.0d0, 0.1275d0, FictKNPR_uint64(1))
#ifdef HAVE_PE 
if (myid.eq. 1) write(*,*)'fbm force'
#endif 

total_lubrication = 0.0d0
! Calculate the forces
call fbm_updateForces(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                      LinSc%valP(NLMAX)%x,&
                      fbm_force_handler_ptr)


!call fbm_compute_particle_reynolds_interface(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
!                                  FictKNPR,Viscosity,Properties%Density(1),mfile, E013)
!call fbm_compute_particle_reynolds_interface_extended(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
!                                                      FictKNPR,Viscosity,Properties%Density(1),mfile, E013, 2)
!
!call fbm_compute_particle_reynolds_farfield(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
!                                            Viscosity,Properties%Density(1),mfile)

!call Sum_myMPI(total_lubrication, global_lubrication)
!call DNA_GetSoosForce(mfile)
!call Get_DissipationIntegral(mfile)


#ifdef HAVE_PE 
if (myid.eq. 1) write(*,*)'fbm update'
#endif 

! Step the particle simulation
call fbm_updateFBM(Properties%Density(1),tstep,timens,&
                   Properties%Gravity,mfile,myid,&
                   QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                   LinSc%valP(NLMAX)%x,&
                   fbm_up_handler_ptr) 

!call fbm_velBCTest()

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

RETURN

END SUBROUTINE Transport_q2p1_UxyzP_fc_ext 
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
use var_QuadScalar, only : tMultiMesh
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

if (myid.ne.0) call updateMixerGeometry(mfile)

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
 CALL Release_cgal_structures()
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
