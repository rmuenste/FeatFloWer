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

!CALL WriteSol_Time(iOut)

if(bViscoElastic)then
  CALL WriteSol_Visco(iOut,0)
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
!
! IF (myid.NE.0) THEN
!  NLMAX = NLMAX + 1
!  ILEV = myExport%Level
!  CALL SETLEV(2)
!  CALL Output_VTK_piece(iOutput,&
!    mg_mesh%level(ILEV)%dcorvg,&
!    mg_mesh%level(ILEV)%kvert)
!  NLMAX = NLMAX - 1
! ELSE
!  CALL Output_VTK_main(iOutput)
! END IF

END IF

!if(myid.eq.1 .and. myFBM%nParticles.gt.0)then
!  call writeparticles(iOutput)
!end if

END SUBROUTINE 
!
!----------------------------------------------
!
SUBROUTINE TimeStepCtrl(dt,inlU,inlT, filehandle)

  USE PP3D_MPI,only :myid,ShowID

  INTEGER IADTIM

  REAL*8  TIMEMX,DTMIN,DTMAX

  COMMON /NSADAT/ TIMEMX,DTMIN,DTMAX,IADTIM

  integer, intent(in) :: filehandle

  INTEGER :: inlU,inlT
  INTEGER :: iMem,nMEm=2
  REAL*8  :: dt, dt_old
  CHARACTER(len=9) :: char_dt
  DATA iMem/0/

  IF (IADTIM.EQ.0) RETURN

  iMem = iMem + 1
  dt_old = dt
  IF (((inlU.GT.3).OR. (inlT.GT.5)).AND.iMem.GE.nMem) THEN
    dt=MAX(dt/1.1d0,DTMIN)
    WRITE(char_dt,'(D9.2)') dt
    READ (char_dt,'(D9.2)') dt
  END IF
  IF (((inlU.LT.3).AND.(inlT.LT.4)).AND.iMem.GE.nMem) THEN
    dt=MIN(1.1d0*dt,DTMAX)
    WRITE(char_dt,'(D9.2)') dt
    READ (char_dt,'(D9.2)') dt
  END IF

  IF (dt.NE.dt_old.AND.myid.eq.ShowID) THEN
    WRITE(MTERM,1) dt_old,dt
    WRITE(filehandle,1) dt_old,dt
  END IF

  IF (dt.NE.dt_old) iMem = 0

  1  FORMAT('Time step change from ',D9.2,' to ',D9.2)

END SUBROUTINE TimeStepCtrl
!
! ----------------------------------------------
!
subroutine postprocessing_fc_ext(dout, iogmv, inlU,inlT,filehandle)

include 'defs_include.h'

implicit none

integer, intent(in) :: filehandle

integer, intent(inout) :: iogmv
real, intent(inout) :: dout

INTEGER :: inlU,inlT,MFILE

! Output the solution in GMV or GiD format
IF (itns.eq.1) THEN
  CALL ZTIME(myStat%t0)
  CALL Output_Profiles(0)
  CALL ZTIME(myStat%t1)
  myStat%tGMVOut = myStat%tGMVOut + (myStat%t1-myStat%t0)
END IF

IF(dout.LE.(timens+1e-10)) THEN

  iOGMV = NINT(timens/dtgmv)
  IF (itns.ne.1) THEN
    CALL ZTIME(myStat%t0)
    CALL Output_Profiles(iOGMV)
    CALL ZTIME(myStat%t1)
    myStat%tGMVOut = myStat%tGMVOut + (myStat%t1-myStat%t0)
  END IF
  dout=dout+dtgmv

  ! Save intermediate solution to a dump file
  IF (insav.NE.0.AND.itns.NE.1) THEN
    IF (MOD(iOGMV,insav).EQ.0) THEN
      CALL ZTIME(myStat%t0)
      CALL SolToFile(-1)
      CALL ZTIME(myStat%t1)
      myStat%tDumpOut = myStat%tDumpOut + (myStat%t1-myStat%t0)
    END IF
  END IF

END IF

! Timestep control
CALL TimeStepCtrl(tstep,inlU,inlT,filehandle)

! Interaction from user
CALL ProcessControl(filehandle,mterm)

end subroutine postprocessing_fc_ext
!
! ----------------------------------------------
!
subroutine handle_statistics(dttt0, istepns)

include 'defs_include.h'

implicit none

real, intent(inout) :: dttt0
integer, intent(inout) :: istepns
real :: dttx = 0.0

! Statistics reset
IF (istepns.eq.1) THEN
  CALL ResetTimer()
  CALL ZTIME(dttt0)
END IF

IF (MOD(istepns,10).EQ.5) THEN
  CALL ZTIME(dttx)
  CALL StatOut(dttx-dttt0,0)
END IF

end subroutine handle_statistics
!
! ----------------------------------------------
!
subroutine print_time(dtimens, dtimemx, dt, istepns, istepmaxns, ufile,uterm)

include 'defs_include.h'
USE PP3D_MPI,only :myid,ShowID

implicit none

real*8, intent(inout) :: dtimens, dtimemx, dt
integer, intent(inout) :: istepns , istepmaxns
integer, intent(inout) :: ufile, uterm

IF (myid.eq.showid) THEN
  write(uTERM,5)
  write(uFILE,5)
  write(uTERM,'(A5,D12.4,A1,D12.4,A8,I8,A1,I8,A7,D12.4)')&
    "time:", dtimens,"/",dtimemx," | itns:",istepns,"/",istepmaxns,&
    " | dt:",dt
  write(uFILE,'(A5,D12.4,A1,D12.4,A8,I8,A1,I8,A7,D12.4)')&
    "time:", dtimens,"/",dtimemx," | itns:",istepns,"/",istepmaxns,&
    " | dt:",dt
  write(uTERM,5)
  write(uFILE,5)
END IF

5 FORMAT(104('='))

end subroutine print_time
!
! ----------------------------------------------
!
subroutine release_mesh() 
USE PP3D_MPI, ONLY : myid,master,showid
USE var_QuadScalar, only : mg_mesh
implicit none

integer :: i
integer :: maxlevel

maxlevel = mg_mesh%maxlevel

if(associated(mg_mesh%level(maxlevel)%dcorvg))then
  deallocate(mg_mesh%level(maxlevel)%dcorvg)
  mg_mesh%level(maxlevel)%dcorvg => null()
end if

end subroutine release_mesh
!
! ----------------------------------------------
!
      SUBROUTINE StatOut_mod(time_passed,myOutFile)
      USE def_FEAT
      USE PP3D_MPI, ONLY : myid,master,showid,subnodes
      USE var_QuadScalar, ONLY : myStat,bNonNewtonian,myMatrixRenewal
      IMPLICIT NONE

      Real, intent(in) :: time_passed

      REAL*8 daux,daux1,ds
      INTEGER myFile,myOutFile,itms,istat
      LOGICAL bExist

      itms = min(itns-1,nitns-1)
      ds = DBLE(subnodes)

      IF (myid.eq.showid) THEN

      IF (myOutFile.eq.0) THEN
       myFile = 669
       OPEN (UNIT=myFile, FILE='_data/Statistics.txt',action='write',iostat=istat)
       if(istat .ne. 0)then
         write(*,*)"Could not open file for writing in StatOut(). "
       stop          
       end if
      ELSE
       myFile = myOutFile
      END IF

      daux = myStat%tKMat+myStat%tDMat+myStat%tMMat+myStat%tCMat+myStat%tSMat
      daux1 = myStat%tGMVOut +myStat%tDumpOut
      WRITE(myFile,*) 
      WRITE(myFile,8) " Overall time            ",time_passed
      WRITE(myFile,*)  
      WRITE(myFile,8) " Solving time            ",time_passed-daux-daux1
      WRITE(myFile,*) 
      WRITE(myFile,8) " Operator assembly time  ",daux
      WRITE(myFile,8) "  Convection matrix      ",myStat%tKMat
      WRITE(myFile,8) "  Deformation matrix     ",myStat%tSMat
      WRITE(myFile,8) "  Diffusion matrix       ",myStat%tDMat
      WRITE(myFile,8) "  Mass matrix            ",myStat%tMMat
      WRITE(myFile,8) "  Reactive term          ",myStat%tCMat
      WRITE(myFile,*) 
      WRITE(myFile,8) " Output time             ",daux1
      WRITE(myFile,8) "  GMV output             ",myStat%tGMVOut
      WRITE(myFile,8) "  Dump file output       ",myStat%tDumpOut

      IF (myOutFile.eq.0) THEN
       CLOSE (myFile)
      END IF

      END IF

8     FORMAT(A24,' : ',F12.4,'s')

      END SUBROUTINE StatOut_mod
!
! ----------------------------------------------
!
subroutine sim_finalize(dttt0, filehandle)

USE PP3D_MPI, ONLY : myid,master,showid,Barrier_myMPI

real, intent(inout) :: dttt0
integer, intent(in) :: filehandle

integer :: ierr
integer :: terminal = 6
real :: time,time_passed


CALL ZTIME(time)

time_passed = time - dttt0
CALL StatOut_mod(time_passed,filehandle)

CALL StatOut_mod(time_passed,terminal)

! Save the final solution vector in unformatted form
!CALL mySolToFile(-1)
CALL Output_Profiles(0)

IF (myid.eq.showid) THEN
  WRITE(MTERM,*) "CC3D_iso_adaptive has successfully finished. "
  WRITE(filehandle,*) "CC3D_iso_adaptive has successfully finished. "
END IF

call release_mesh()

CALL Barrier_myMPI()
CALL MPI_Finalize(ierr)

end subroutine sim_finalize
