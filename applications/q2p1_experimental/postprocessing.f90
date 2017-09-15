SUBROUTINE mySolToFile(iOutput)
USE def_FEAT
USE Transport_Q2P1,ONLY:QuadSc,LinSc,bViscoElastic
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
USE Transport_Q1,ONLY:Tracer
USE PP3D_MPI, ONLY:myid

IMPLICIT NONE
INTEGER iOutput

! -------------- workspace -------------------
INTEGER  NNWORK
PARAMETER (NNWORK=1)
INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)

INTEGER            :: KWORK(1)
REAL               :: VWORK(1)
DOUBLE PRECISION   :: DWORK(NNWORK)

COMMON       NWORK,IWORK,IWMAX,L,DWORK
EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
! -------------- workspace -------------------
INTEGER ifilen,iOut,nn
DATA ifilen/0/

IF (iOutput.LT.0) THEN
 ifilen=ifilen+1
 iOut=MOD(ifilen+insavn-1,insavn)+1
ELSE
 iOut = iOutput
END IF

CALL WriteSol_Velo(iOut,0,QuadSc%ValU,QuadSc%ValV,QuadSc%ValW)
nn = knel(nlmax)
CALL WriteSol_Pres(iOut,0,LinSc%ValP(NLMAX)%x,LinSc%AuxP(NLMAX)%x(1),&
     LinSc%AuxP(NLMAX)%x(nn+1),LinSc%AuxP(NLMAX)%x(2*nn+1),LinSc%AuxP(NLMAX)%x(3*nn+1),nn)
! CALL WriteSol_Coor(iOut,0,DWORK(L(KLCVG(NLMAX))),QuadSc%AuxU,QuadSc%AuxV,QuadSc%AuxW,QuadSc%ndof)
CALL WriteSol_Time(iOut)

if(bViscoElastic)then
  CALL WriteSol_Visco(iOut,0)
end if

if(myid.eq.1 .and. myFBM%nParticles.gt.0)then
  call writeparticles(iOut)
end if

END SUBROUTINE mySolToFile
!
! ----------------------------------------------
!
SUBROUTINE myOutput_Profiles(iOutput)
USE def_FEAT
USE Transport_Q2P1,ONLY:QuadSc,LinSc,PressureToGMV,&
    Viscosity,Distance,Distamce,mgNormShearStress
USE Transport_Q1, ONLY:Tracer
USE PP3D_MPI, ONLY:myid,showid,Comm_Summ
USE var_QuadScalar,ONLY:myExport,myFBM,mg_mesh
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
! USE PLinScalar,ONLY:PLinScP1toQ1,OutputInterphase,PLinLS,&
!                dNorm,IntPhaseElem,FracFieldQ1
IMPLICIT NONE
INTEGER iOutput,mfile

! -------------- workspace -------------------
INTEGER  NNWORK
PARAMETER (NNWORK=1)
INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)

INTEGER            :: KWORK(1)
REAL               :: VWORK(1)
DOUBLE PRECISION   :: DWORK(NNWORK)

COMMON       NWORK,IWORK,IWMAX,L,DWORK
EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
! -------------- workspace -------------------

IF     (myExport%Format.EQ."GMV") THEN

 IF (myid.NE.0) THEN
  NLMAX = NLMAX + 1
  ILEV = myExport%Level
  CALL SETLEV(2)
  CALL Output_GMV_fields(iOutput,&
    mg_mesh%level(ILEV)%dcorvg,&
    mg_mesh%level(ILEV)%kvert)
  NLMAX = NLMAX - 1
 END IF

ELSEIF (myExport%Format.EQ."VTK") THEN

 IF (myid.NE.0) THEN
  NLMAX = NLMAX + 1
  ILEV = myExport%Level
  CALL SETLEV(2)
  CALL Output_VTK_piece(iOutput,&
    mg_mesh%level(ILEV)%dcorvg,&
    mg_mesh%level(ILEV)%kvert)
  NLMAX = NLMAX - 1
 ELSE
  CALL Output_VTK_main(iOutput)
 END IF

END IF

END SUBROUTINE 
!
! ----------------------------------------------
!
subroutine sim_finalize(dttt0, filehandle)
include 'defs_include.h'
use Mesh_Structures, only: release_mesh
USE PP3D_MPI, ONLY : myid,master,showid,Barrier_myMPI
use var_QuadScalar, only: istep_ns,mg_mesh
use sol_out, only: write_sol_to_file

real, intent(inout) :: dttt0
integer, intent(in) :: filehandle

integer :: ierr
integer :: terminal = 6
real :: time,time_passed

CALL ZTIME(time)

time_passed = time - dttt0
CALL StatOut(time_passed,0)

CALL StatOut(time_passed,terminal)

! Save the final solution vector in unformatted form
istep_ns = istep_ns - 1
call write_sol_to_file(insavn, timens)

IF (myid.eq.showid) THEN
  WRITE(*,*) "PP3D_LES has successfully finished. "
  WRITE(filehandle,*) "PP3D_LES has successfully finished. "
END IF

call release_mesh(mg_mesh)
CALL Barrier_myMPI()
CALL MPI_Finalize(ierr)

end subroutine sim_finalize
