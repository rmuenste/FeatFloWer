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
USE Transport_Q1,ONLY:Tracer
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
  CALL myOutput_VTK_piece(iOutput,&
    mg_mesh%level(ILEV)%dcorvg,&
    mg_mesh%level(ILEV)%kvert)
  NLMAX = NLMAX - 1
 ELSE
  CALL myOutput_VTK_main(iOutput)
 END IF

END IF

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
subroutine postprocessing_laplace(dout, iogmv, inlU,inlT,filehandle)

include 'defs_include.h'

implicit none

integer, intent(in) :: filehandle

integer, intent(inout) :: iogmv
real, intent(inout) :: dout

INTEGER :: inlU,inlT,MFILE

! Output the solution in GMV or GiD format
IF (itns.eq.1) THEN
  CALL ZTIME(myStat%t0)
  CALL myOutput_Profiles(0)
  CALL ZTIME(myStat%t1)
  myStat%tGMVOut = myStat%tGMVOut + (myStat%t1-myStat%t0)
END IF

IF(dout.LE.(timens+1e-10)) THEN

  iOGMV = NINT(timens/dtgmv)
  IF (itns.ne.1) THEN
    CALL ZTIME(myStat%t0)
    CALL myOutput_Profiles(iOGMV)
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

end subroutine postprocessing_laplace
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
subroutine sim_finalize(dttt0, filehandle)

USE PP3D_MPI, ONLY : myid,master,showid,Barrier_myMPI

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
!CALL SolToFile(-1)

IF (myid.eq.showid) THEN
  WRITE(MTERM,*) "PP3D_LES has successfully finished. "
  WRITE(filehandle,*) "PP3D_LES has successfully finished. "
END IF

call release_mesh()

CALL Barrier_myMPI()
CALL MPI_Finalize(ierr)

end subroutine sim_finalize

SUBROUTINE myOutput_VTK_piece(iO,dcoor,kvert)
USE def_FEAT
USE  PP3D_MPI, ONLY:myid,showid,subnodes
USE Transport_Q2P1,ONLY: QuadSc,LinSc,Viscosity,Distance,Distamce,mgNormShearStress,myALE
USE Transport_Q2P1,ONLY: MixerKnpr,FictKNPR,ViscoSc, myBoundary,myQ2Coor
USE Transport_Q1,ONLY:Tracer
USE var_QuadScalar,ONLY:myExport, Properties, bViscoElastic,myFBM
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
USE var_QuadScalar,ONLY:mg_mesh

IMPLICIT NONE
REAL*8 dcoor(3,*)
INTEGER kvert(8,*),iO,ioffset,ive,ivt,iField,i,istat
CHARACTER fileid*(5),filename*(27),procid*(3)
INTEGER NoOfElem,NoOfVert
REAL*8,ALLOCATABLE ::  tau(:,:)
REAL*8 psi(6)
integer :: iunit = 908070
integer :: iPhase
real*8 :: dp,daux,daux2

NoOfElem = KNEL(ILEV)
NoOfVert = KNVT(ILEV)

filename=" "
WRITE(filename(1:),'(A,I5.5,A4)') '_vtk/res_node_***.',iO,".vtu"

IF(myid.eq.showid) WRITE(*,'(104("="))') 
IF(myid.eq.showid) WRITE(*,*) "Outputting vtk file into ",filename
WRITE(filename(15:17),'(I3.3)') myid

OPEN (UNIT=iunit,FILE=filename,action='write',iostat=istat)
if(istat .ne. 0)then
  write(*,*)"Could not open file for writing. "
  stop          
end if

write(iunit, *)"<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">"
write(iunit, *)"  <UnstructuredGrid>"
write(iunit, *)"    <Piece NumberOfPoints=""",KNVT(ILEV),""" NumberOfCells=""",NoOfElem,""">"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the node field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <PointData>"

DO iField=1,SIZE(myExport%Fields)

 SELECT CASE(ADJUSTL(TRIM(myExport%Fields(iField))))
 CASE('Velocity')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Velocity",""" NumberOfComponents=""3"" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,3E16.7)')"        ",1.0,&
     0.0,&
     0.0
  end do
  write(iunit, *)"        </DataArray>"

   
 CASE('Stress')
  if(bViscoElastic)then

  ALLOCATE(tau(6,NoOfVert))
  DO i=1,NoOfVert
   psi = [ViscoSc%Val11(i),ViscoSc%Val22(i),ViscoSc%Val33(i),&
          ViscoSc%Val12(i),ViscoSc%Val13(i),ViscoSc%Val23(i)]
   CALL ConvertPsiToTau(psi,tau(:,i))   
     tau(1,i) = (tau(1,i) - 1d0)/Properties%ViscoLambda
     tau(2,i) = (tau(2,i) - 1d0)/Properties%ViscoLambda
     tau(3,i) = (tau(3,i) - 1d0)/Properties%ViscoLambda
     tau(4,i) = (tau(4,i) - 0d0)/Properties%ViscoLambda
     tau(5,i) = (tau(5,i) - 0d0)/Properties%ViscoLambda
     tau(6,i) = (tau(6,i) - 0d0)/Properties%ViscoLambda
  END DO

  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Stress",""" NumberOfComponents=""6"" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,6E16.7)')"        ",REAL(tau(1,ivt)),REAL(tau(2,ivt)),REAL(tau(3,ivt)),REAL(tau(4,ivt)),REAL(tau(5,ivt)),REAL(tau(6,ivt))
  end do

  DEALLOCATE(tau)
  write(iunit, *)"        </DataArray>"
  end if
 
 CASE('MeshVelo')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","MeshVelocity",""" NumberOfComponents=""3"" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,3E16.7)')"        ",REAL(myALE%MeshVelo(1,ivt)),REAL(myALE%MeshVelo(2,ivt)),REAL(myALE%MeshVelo(3,ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('Pressure_V')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Pressure_V",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(LinSc%Q2(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('Temperature')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Temperature",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(Tracer%Val(NLMAX)%x(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('Mixer')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Mixer",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(FictKNPR(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('Viscosity')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Viscosity",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(Viscosity(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('Monitor')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Monitor",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(myALE%monitor(ivt))
  end do
  write(iunit, *)"        </DataArray>"
 CASE('PhaseProp_V')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","PhaseProp_V",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(myBoundary%iPhase(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 END SELECT 

END DO

write(iunit, '(A)')"    </PointData>"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the cell field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <CellData>"

DO iField=1,SIZE(myExport%Fields)

 SELECT CASE(ADJUSTL(TRIM(myExport%Fields(iField))))
 CASE('Pressure_E')
  IF (ILEV.EQ.NLMAX-1) THEN
   write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Pressure_E",""" format=""ascii"">"
   do ivt=1,NoOfElem
    ive = 4*(ivt-1)+1
    write(iunit, '(A,E16.7)')"        ",REAL(LinSc%ValP(NLMAX-1)%x(ive))
   end do
   write(iunit, *)"        </DataArray>"
  END IF
  IF (ILEV.EQ.NLMAX) THEN
   write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Pressure_E",""" format=""ascii"">"
   do ivt=1,NoOfElem
    ive = 4*(ivt-1)+1
    IF (ivt.le.NoOfElem/8d0) THEN
     dP = LinSc%ValP(NLMAX-1)%x(ive)
    ELSE
     ive = ivt - NoOfElem/8d0 - 1
     ive = (ive - mod(ive,7))/7 + 1
     ive = 4*(ive-1)+1
     dP = LinSc%ValP(NLMAX-1)%x(ive)
    END IF
    write(iunit, '(A,E16.7)')"        ",REAL(dP)
   end do
   write(iunit, *)"        </DataArray>"
  END IF

 CASE('PhaseProp_E')
  write(iunit, '(A,A,A)')"        <DataArray type=""UInt8"" Name=""","PhaseProp_E",""" format=""ascii"">"
  do ivt=1,NoOfElem

  iPhase = 0
!   IF (ivt.gt.KNVT(ILEV-1)+KNET(ILEV-1)+KNAT(ILEV-1)) THEN
!    iPhase = myBoundary%iPhase(KNVT(ILEV-1)+KNET(ILEV-1)+KNAT(ILEV-1)+ivt)
!  end if
!
!   IF (ivt.le.NoOfElem/8d0) THEN
!    iPhase = myBoundary%iPhase(KNVT(ILEV-1)+KNET(ILEV-1)+KNAT(ILEV-1)+ivt)
!    write(*,*)"phase if: ",iPhase
!   ELSE
!    ive = ivt - NoOfElem/8d0 - 1
!    ive = (ive - mod(ive,7))/7 + 1
!    iPhase = myBoundary%iPhase(KNVT(ILEV-1)+KNET(ILEV-1)+KNAT(ILEV-1)+ive)
!    write(*,*)"phase else: ",iPhase,daux2
!   END IF
   write(iunit, '(A,I3)')"        ",iPhase
   
!    IF (bPhase) THEN
!    ELSE
!     write(iunit, '(A,I3)')"        ",0
!    END IF
  end do
  write(iunit, *)"        </DataArray>"

 CASE('Viscosity_E')
  IF (ILEV.EQ.NLMAX-1) THEN
   write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Viscosity_E",""" format=""ascii"">"
   do ivt=1,NoOfElem
    write(iunit, '(A,E16.7)')"        ",REAL(Viscosity((nvt+net+nat+ivt)))
   end do
   write(iunit, *)"        </DataArray>"
  END IF

 END SELECT

END DO

write(iunit, '(A)')"    </CellData>"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the mesh data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"      <Points>"
write(iunit, '(A)')"        <DataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3"" format=""ascii"" RangeMin=""0"" RangeMax=""1.0"">"
do ivt=1,NoOfVert
 write(iunit,'(A10,3E16.7)')"          ",REAL(myQ2Coor(1,ivt)),REAL(myQ2Coor(2,ivt)),REAL(myQ2Coor(3,ivt))
end do
write(iunit, *)"        </DataArray>"
write(iunit, *)"      </Points>"

write(iunit, *)"      <Cells>"
write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""connectivity"" format=""ascii"" RangeMin=""0"" RangeMax=""",NoOfElem-1,""">"
do ive=1,NoOfElem   
 write(iunit, '(8I10)')kvert(1,ive)-1,kvert(2,ive)-1,kvert(3,ive)-1,kvert(4,ive)-1,&
                       kvert(5,ive)-1,kvert(6,ive)-1,kvert(7,ive)-1,kvert(8,ive)-1
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""offsets"" format=""ascii"" RangeMin=""8"" RangeMax=""",8*NoOfElem,""">"
ioffset=NoOfElem/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')ive*8,(ive+1)*8,(ive+2)*8,(ive+3)*8,(ive+4)*8,(ive+5)*8,(ive+6)*8,(ive+7)*8
end do

do ive=ioffset+1,NoOfElem
 write(iunit, '(I10)')ive*8
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A)')"        <DataArray type=""UInt8"" Name=""types"" format=""ascii"" RangeMin=""12"" RangeMax=""12"">"
ioffset=NoOfElem/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')12,12,12,12,12,12,12,12
end do
do ive=ioffset+1,NoOfElem
 write(iunit, '(I10)')12
end do
write(iunit, '(A)')"        </DataArray>"
 
write(iunit, *)"      </Cells>"
write(iunit, *)"    </Piece>"
   
write(iunit, *)"  </UnstructuredGrid>"
write(iunit, *)"</VTKFile>"
close(iunit)

END SUBROUTINE myOutput_VTK_piece



SUBROUTINE myOutput_VTK_main(iO)
USE  PP3D_MPI, ONLY:myid,showid,subnodes
USE var_QuadScalar,ONLY:myExport,bViscoElastic
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
USE def_FEAT

IMPLICIT NONE
INTEGER iO,iproc,iField
INTEGER :: iMainUnit=555
CHARACTER mainname*(20) 
CHARACTER filename*(26)

integer :: istat

! generate the file name
mainname=' '
WRITE(mainname(1:),'(A,I5.5,A5)') '_vtk/main.',iO,'.pvtu'

OPEN (UNIT=imainunit,FILE=mainname,action='write',iostat=istat)
if(istat .ne. 0)then
  write(*,*)"Could not open file for writing. "
  stop          
end if

write(imainunit, *)"<VTKFile type=""PUnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">"
write(imainunit, *)"  <PUnstructuredGrid GhostLevel=""0"">"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the node field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(imainunit, '(A)')"    <PPointData>"

DO iField=1,SIZE(myExport%Fields)

 SELECT CASE(ADJUSTL(TRIM(myExport%Fields(iField))))
 CASE('Velocity')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Velocity",""" NumberOfComponents=""3""/>"
 CASE('Stress')
  if(bViscoElastic)then
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Stress",""" NumberOfComponents=""6""/>"
  end if
 CASE('MeshVelo')
  write(imainunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","MeshVelocity",""" NumberOfComponents=""3""/>"
 CASE('Pressure_V')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Pressure_V","""/>"
 CASE('Temperature')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Temperature","""/>"
 CASE('Mixer')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Mixer","""/>"
 CASE('Viscosity')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Viscosity","""/>"
 CASE('Monitor')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Monitor","""/>"
 CASE('PhaseProp_V')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","PhaseProp_V","""/>"
 END SELECT
END DO

write(imainunit, '(A)')"    </PPointData>"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the cell field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(imainunit, '(A)')"    <PCellData>"
DO iField=1,SIZE(myExport%Fields)

 SELECT CASE(ADJUSTL(TRIM(myExport%Fields(iField))))
 CASE('Pressure_E')
   write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Pressure_E","""/>"
 CASE('Viscosity_E')
!  WRITE(*,*) myExport%Level,myExport%LevelMax,myExport%Level.EQ.myExport%LevelMax
  IF (myExport%Level.EQ.myExport%LevelMax) THEN
   write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Viscosity_E","""/>"
  END IF
 CASE('PhaseProp_E')
   write(imainunit, '(A,A,A)')"       <PDataArray type=""UInt8"" Name=""","PhaseProp_E","""/>"
 END SELECT
END DO
write(imainunit, '(A)')"    </PCellData>"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the mesh data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(imainunit, *)"    <PPoints>"
write(imainunit, *)"      <PDataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3""/>"
write(imainunit, *)"    </PPoints>"

do iproc=1,subnodes
 filename=" "
 WRITE(filename(1:),'(A9,I3.3,A1,I5.5,A4)') 'res_node_',iproc,'.',iO,".vtu"
 write(imainunit, '(A,A,A)')"      <Piece Source=""",trim(adjustl(filename)),"""/>"  
end do
write(imainunit, *)"  </PUnstructuredGrid>"
write(imainunit, *)"  </VTKFile>"
close(imainunit)

END SUBROUTINE myOutput_VTK_main
