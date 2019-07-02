module visualization_out

use cinterface

use var_QuadScalar, only:knvt,knet,knat,knel,tMultiMesh,tQuadScalar,tLinScalar, &
                         t1DOutput,MlRhoMat, mg_MlRhomat, MixerKNPR, temperature, mg_mesh

use Sigma_User, only: tOutput, tSigma,dMinOutputPressure

use  PP3D_MPI, only:myid,showid,subnodes, comm_summn
use  PP3D_MPI, only:master,COMM_Maximum,COMM_Maximumn,&
                   COMM_Minimumn,COMM_NLComplete,Comm_Summ,myMPI_Barrier
!------------------------------------------------------------------------------------------------
! A module for output routines for vtk, gmv or other output formats.
! For an application that needs strongly different visualization output
! than the standard, these output routines should be implemented here.
!------------------------------------------------------------------------------------------------

contains
!
!-------------------------------------------------------------------------------------------------
! A routine for outputting fields for an sse application
!-------------------------------------------------------------------------------------------------
! @param sExport An export structure
! @param iOutput The output idx of the visulation file
! @param sQuadSc The velocity solution structure of the mesh
! @param sLinSc The pressure solution structure of the mesh
! @param sTracer The solution of the scalar tracer equation
! @param visc The viscosity solution
! @param dist The distance solution
! @param shear The shear rate solution
! @param mgMesh The mesh that will be written out
!subroutine viz_output_fields(sExport, iOutput, sQuadSc, sLinSc, sTracer, visc, dist, shear, mgMesh)
subroutine viz_output_fields(sExport, iOutput, sQuadSc, sLinSc, visc, screw, shell, shear, mgMesh)

use var_QuadScalar, only:tExport

USE PP3D_MPI, ONLY:myid
USE def_FEAT

USE Transport_Q1,ONLY:Tracer

implicit none

type(tExport), intent(in) :: sExport

integer, intent(in) :: iOutput

type(tQuadScalar), intent(in) :: sQuadSc

type(tLinScalar), intent(in) :: sLinSc

real*8, dimension(:), intent(in) :: visc

real*8, dimension(:), intent(in) :: screw

real*8, dimension(:), intent(in) :: shell

real*8, dimension(:), intent(in) :: shear

type(tMultiMesh), intent(in) :: mgMesh

! locals
integer :: ioutput_lvl

!type(fieldPtr), dimension(3) :: packed

if (sExport%Format .eq. "VTK") then

 call viz_OutputHistogram(iOutput, sQuadSc, mgMesh%nlmax)

 call viz_OutPut_1D(iOutput, sQuadSc, sLinSc, Tracer, mgMesh%nlmax)
 
 if (myid.ne.0) then

  ioutput_lvl = sExport%Level

  call viz_write_vtu_process(iOutput,&
    mgMesh%level(ioutput_lvl)%dcorvg,&
    mgMesh%level(ioutput_lvl)%kvert, &
    sQuadSc, sLinSc, visc, screw, shell, shear, ioutput_lvl,&
    mgMesh)
 else

  call viz_write_pvtu_main(iOutput)
 end if

end if

end subroutine viz_output_fields
!
!------------------------------------------------------------------------------------------------
! Routine for writing out a VTK visualization of a process
!------------------------------------------------------------------------------------------------
! @param iO Output file idx
! @param dcoor Coordinate array of the mesh
! @param kvert Vertices at an element connectivity array of the mesh
! @param sQuadSc The velocity solution structure of the mesh
! @param sLinSc The pressure solution structure of the mesh
! @param visc The viscosity solution
! @param dist The distance solution
! @param shear The shearrate solution
! @param ioutput_lvl Output level of the mesh
! @param mgMesh The mesh that will be written out
subroutine viz_write_vtu_process(iO,dcoor,kvert, sQuadSc, sLinSc, visc, screw, shell, shear,&
                                 ioutput_lvl, mgMesh)

use var_QuadScalar,only:myExport,MixerKnpr

implicit none

integer, intent(in) :: iO

real*8, dimension(:,:), intent(in) :: dcoor

integer, dimension(:,:), intent(in) :: kvert

type(tQuadScalar), intent(in) :: sQuadSc

type(tLinScalar), intent(in) :: sLinSc

real*8, dimension(:), intent(in) :: visc

real*8, dimension(:), intent(in) :: screw

real*8, dimension(:), intent(in) :: shell

real*8, dimension(:), intent(in) :: shear

integer, intent(in) :: ioutput_lvl

type(tMultiMesh), intent(in) :: mgMesh

! locals
real*8 :: DXX,DYY,DZZ

real*8 :: myPI = dATAN(1d0)*4d0

integer :: NoOfElem,NoOfVert

integer :: istat, ioffset,ive,ivt,iField

! local
integer :: iunit = 908070

integer :: inlmax

character fileid*(5),filename*(27),procid*(3)

NoOfElem = KNEL(ioutput_lvl)
NoOfVert = KNVT(ioutput_lvl)

inlmax = mgMesh%maxlevel

filename=" "

write(filename(1:),'(A,I5.5,A4)') '_vtk/res_node_***.',iO,".vtu"

if(myid.eq.showid) write(*,'(104("="))')
if(myid.eq.showid) write(*,*) "Outputting vtk file into ",filename
write(filename(15:17),'(I3.3)') myid

open (unit=iunit,file=filename,action='write',iostat=istat)
if(istat .ne. 0)then
  write(*,*)"Could not open file for writing. "
  stop
end if

write(iunit, *)"<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">"
write(iunit, *)"  <UnstructuredGrid>"
write(iunit, *)"    <Piece NumberOfPoints=""",KNVT(ioutput_lvl),""" NumberOfCells=""",NoOfElem,""">"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the node field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <PointData>"

do iField=1,size(myExport%Fields)

 select case(adjustl(trim(myExport%Fields(iField))))
 case('Velocity')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Velocity [m/s]",""" NumberOfComponents=""3"" format=""ascii"">"
  do ivt=1,NoOfVert
   DXX = 1d-2*sQuadSc%ValU(ivt)
   DYY = 1d-2*sQuadSc%ValV(ivt)
   DZZ = 1d-2*sQuadSc%ValW(ivt) !- myProcess%FrameVelocity
   write(iunit, '(A,3E16.7)')"        ",REAL(DXX),REAL(DYY),REAL(DZZ)
  end do
  write(iunit, *)"        </DataArray>"

! CASE('MeshVelo')
!  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","MeshVelocity",""" NumberOfComponents=""3"" format=""ascii"">"
!  do ivt=1,NoOfVert
!   write(iunit, '(A,3E16.7)')"        ",REAL(myALE%MeshVelo(1,ivt)),REAL(myALE%MeshVelo(2,ivt)),REAL(myALE%MeshVelo(3,ivt))
!  end do
!  write(iunit, *)"        </DataArray>"

 case('Pressure_V')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Pressure [bar]",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(1e-5*0.1d0*sLinSc%Q2(ivt)-dMinOutputPressure)
  end do
  write(iunit, *)"        </DataArray>"

! CASE('Temperature')
!  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Temperature",""" format=""ascii"">"
!  do ivt=1,NoOfVert
!   write(iunit, '(A,E16.7)')"        ",REAL(Tracer%Val(NLMAX)%x(ivt))
!  end do
!  write(iunit, *)"        </DataArray>"
!
 case('Shell')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Shell",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(1d1*Shell(ivt))
  end do
  write(iunit, *)"        </DataArray>"

  case('FBM')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","FBM",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(MixerKNPR(ivt))
  end do
  write(iunit, *)"        </DataArray>"

  case('Screw')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Screw",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(1d1*Screw(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 case('Viscosity')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Viscosity [Pa s]",""" format=""ascii"">"

  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(0.1d0*visc(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 case('Temperature')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Temperature [C]",""" format=""ascii"">"

  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(temperature(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 case('NormShearRate')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","NormShearRate [1/s]",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(shear(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 case('Partitioning')
  write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","Partition",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,I10)')"        ", myid
  end do
  write(iunit, *)"        </DataArray>"

 end select

end do


write(iunit, '(A)')"    </PointData>"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the cell field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <CellData>"

do iField=1,size(myExport%Fields)

 select case(adjustl(trim(myExport%Fields(iField))))

 case('Pressure_E')

  if (ioutput_lvl.EQ.inlmax-1)then
   write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Pressure_E",""" format=""ascii"">"
   do ivt=1,NoOfElem
    ive = 4*(ivt-1)+1
    write(iunit, '(A,E16.7)')"        ",REAL(sLinSc%ValP(inlmax-1)%x(ive))
   end do
   write(iunit, *)"        </DataArray>"
  end if

 end select

end do

write(iunit, '(A)')"    </CellData>"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the mesh data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"      <Points>"
write(iunit, '(A)')"        <DataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3"" format=""ascii"" RangeMin=""0"" RangeMax=""1.0"">"

do ivt=1,NoOfVert

 DXX = 1d1*dcoor(1,ivt)
 DYY = 1d1*dcoor(2,ivt)
 DZZ = 1d1*dcoor(3,ivt)

 write(iunit,'(A10,3E16.7)')"          ",REAL(DXX),REAL(DYY),REAL(DZZ)
!  write(iunit,'(A10,3E16.7)')"          ",REAL(dcoor(1,ivt)),REAL(dcoor(2,ivt)),REAL(dcoor(3,ivt))
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

end subroutine viz_write_vtu_process
!
!-------------------------------------------------------------------------------------------------
! A routine for outputting fields for an sse application
!-------------------------------------------------------------------------------------------------
! @param iO Output idx of the visualization file
!
subroutine viz_write_pvtu_main(iO)
USE  PP3D_MPI, ONLY:myid,showid,subnodes
USE var_QuadScalar,ONLY:myExport
USE def_FEAT

IMPLICIT NONE
INTEGER iO,iproc,iField, istat
INTEGER :: iMainUnit=555
CHARACTER mainname*(20)
CHARACTER filename*(26)

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
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Velocity [m/s]",""" NumberOfComponents=""3""/>"
! CASE('MeshVelo')
!  write(imainunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","MeshVelocity",""" NumberOfComponents=""3""/>"
 CASE('Pressure_V')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Pressure [bar]","""/>"
! CASE('Temperature')
!  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Temperature","""/>"
 CASE('Shell')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Shell","""/>"
 CASE('FBM')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","FBM","""/>"
 CASE('Screw')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Screw","""/>"
 CASE('Viscosity')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Viscosity [Pa s]","""/>"
 CASE('Temperature')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Temperature [C]","""/>"
 CASE('NormShearRate')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","NormShearRate [1/s]","""/>"
!  CASE('Monitor')
!  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Monitor","""/>"
!
!  CASE('FBM')
!  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","FBM","""/>"

 case('Partitioning')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Int32"" Name=""","Partition","""/>"

 END SELECT
END DO


write(imainunit, '(A)')"    </PPointData>"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the cell field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(imainunit, '(A)')"    <PCellData>"
DO iField=1,SIZE(myExport%Fields)

 SELECT CASE(ADJUSTL(TRIM(myExport%Fields(iField))))
 CASE('Pressure_E')
!  WRITE(*,*) myExport%Level,myExport%LevelMax,myExport%Level.EQ.myExport%LevelMax
  IF (myExport%Level.EQ.myExport%LevelMax) THEN
   write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Pressure_E","""/>"
  END IF

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

end subroutine viz_write_pvtu_main
!
!-------------------------------------------------------------------------------------------------
! A routine for outputting fields for an sse application
!-------------------------------------------------------------------------------------------------
!
! @param iO Output file idx
! @param sQuadSc The velocity solution structure of the mesh
! @param maxlevel The maximum grid level used in computation (former NLMAX)
SUBROUTINE viz_OutputHistogram(iOut, sQuadSc, maxlevel)
use var_QuadScalar, only : mg_MlRhoMat, MixerKNPR, ShearRate, Viscosity
USE PP3D_MPI, ONLY:myid,showid,Comm_Summ,Comm_Summn
use Sigma_User, only : myOutput
implicit none

interface
  subroutine c_write_histogram() bind(C, name="c_write_histogram")
    use cinterface, only: c1dOutput
    use iso_c_binding
!    type(c1dOutput) :: thetype
!    character(kind=c_char) :: dataName(*)
  end subroutine c_write_histogram
end interface

interface
  subroutine c_add_histogram(histogram) bind(C, name="c_add_histogram")
    use cinterface, only: histoData
    use iso_c_binding
    type(histoData) :: histogram
  end subroutine c_add_histogram
end interface

interface
  subroutine c_init_histogram_section() bind(C, name="c_init_histogram_section")
    use cinterface, only: c1dOutput
    use iso_c_binding
!    type(c1dOutput) :: thetype
!    character(kind=c_char) :: dataName(*)
  end subroutine c_init_histogram_section
end interface



integer :: iOut

type(tQuadScalar), intent(in) :: sQuadSc

integer, intent(in) :: maxlevel

! local variables
type(histoData) :: histData
integer :: i,j,nBin
real*8 logShear,logVisco
real*8,allocatable ::BinEta(:,:),BinVis(:,:)

real*8,allocatable, dimension(:), target :: HistoEta
real*8,allocatable, dimension(:), target :: HistoVis
real*8,allocatable, dimension(:), target :: JsonBin

real*8 Ml_i,dauxVis,dauxEta
real*8 minBinEta,mAXBinEta,minBinVis,mAXBinVis,dBinEta,dBinVis,meanEta,meanVis
character*100 cHistoFile

nBin = myOutput%nOfHistogramBins
if (.not.allocated(HistoEta)) allocate(HistoEta(nBin))
if (.not.allocated(HistoVis)) allocate(HistoVis(nBin))
if (.not.allocated(JsonBin)) allocate(JsonBin(nBin))
if (.not.allocated(BinEta))   allocate(BinEta(3,nBin))
if (.not.allocated(BinVis))   allocate(BinVis(3,nBin))

IF (myid.ne.0) THEN

HistoVis = 0d0
HistoEta = 0d0

BinEta(1,   1)=-1d1
BinVis(1,   1)=-1d1
BinEta(2,nBin)=+1d1
BinVis(2,nBin)=+1d1

minBinEta = log10(myOutput%HistogramShearMin)
maxBinEta = log10(myOutput%HistogramShearMax)
minBinVis = log10(myOutput%HistogramViscoMin)
maxBinVis = log10(myOutput%HistogramViscoMax)

dBinEta = (maxBinEta-minBinEta)/DBLE(nBin-2)
dBinVis = (maxBinVis-minBinVis)/DBLE(nBin-2)

DO i=2,nBin

 BinEta(2,i-1)= minBinEta + DBLE(i-2)*dBinEta
 BinVis(2,i-1)= minBinVis + DBLE(i-2)*dBinVis
 BinEta(3,i-1)= minBinEta + (DBLE(i)-2.5d0)*dBinEta
 BinVis(3,i-1)= minBinVis + (DBLE(i)-2.5d0)*dBinVis
 BinEta(1,i)  = BinEta(2,i-1)
 BinVis(1,i)  = BinVis(2,i-1)

END DO
 BinEta(3,nBin)= minBinEta + (DBLE(nBin)-1.5d0)*dBinEta
 BinVis(3,nBin)= minBinVis + (DBLE(nBin)-1.5d0)*dBinVis

! IF (myid.eq.1) then
!  DO i=1,nBin
!   write(*,'(4ES14.4)') BinEta(:,i), BinVis(:,i)
!  end do
! end if
! pause

DO i=1,sQuadSc%ndof

 IF (MixerKNPR(i).eq.0) THEN
  logShear = LOG10(Shearrate(i))
  logVisco = LOG10(0.1d0*Viscosity(i))
  Ml_i = mg_MlRhoMat(maxlevel)%a(i)

  DO j=1,nBin
   IF (logShear.ge.BinEta(1,j).and.logShear.lt.BinEta(2,j)) HistoEta(j)  = HistoEta(j)  + Ml_i
   IF (logVisco.ge.BinVis(1,j).and.logVisco.lt.BinVis(2,j)) HistoVis(j)  = HistoVis(j)  + Ml_i
  END DO

 END IF

END DO
END IF

CALL COMM_SUMMn(HistoVis,nBin)
CALL COMM_SUMMn(HistoEta,nBin)

IF (myid.eq.1) THEN
 dauxVis = 0d0; dauxEta = 0d0
 DO i=1,nBin
  dauxVis = dauxVis + HistoVis(i)
  dauxEta = dauxEta + HistoEta(i)
 END DO

 meanVis = 0d0
 meanEta = 0d0

 DO i=1,nBin
  meanVis = meanVis + HistoVis(i)*BinVis(3,i)/dauxVis
  meanEta = meanEta + HistoEta(i)*BinEta(3,i)/dauxEta
 END DO


 WRITE(cHistoFile,'(A,I5.5,A)') "_hist/h_",iOut, ".txt"
 WRITE(*,'(3A)')"'",ADJUSTL(TRIM(cHistoFile)),"'"
 OPEN(FILE=ADJUSTL(TRIM(cHistoFile)),UNIT=652)
 DO i=1,nBin
  WRITE(652,'(2(2(ES14.4),A))') BinEta(3,i),HistoEta(i)/dauxEta, " | ",BinVis(3,i),HistoVis(i)/dauxVis
 END DO
  WRITE(652,'(A)') "------------------------------------------------------"
  WRITE(652,'(2(1(ES14.4),A))') 10d0**meanEta, " | ",10d0**meanVis

 CLOSE(652)

 DO i=1,nBin
  HistoEta(i) = HistoEta(i)/dauxEta
  HistoVis(i) = HistoVis(i)/dauxVis
  JsonBin(i) = BinEta(3,i)
 END DO

 call c_init_histogram_section()

 histData%binPos = c_loc(JsonBin)
 histData%binHeight = c_loc(HistoEta)
 histData%mean = 10d0**meanEta 
 histData%length = nBin 
 call cf_make_c_char(C_CHAR_"Eta"//C_NULL_CHAR, histData%name )

 call c_add_histogram(histData)

 DO i=1,nBin
  JsonBin(i) = BinVis(3,i)
 END DO

 histData%binPos = c_loc(JsonBin)
 histData%binHeight = c_loc(HistoVis)
 histData%mean = 10d0**meanVis 
 call cf_make_c_char(C_CHAR_"Viscosity"//C_NULL_CHAR, histData%name )

 call c_add_histogram(histData)

END IF

end subroutine viz_OutputHistogram
!
!-------------------------------------------------------------------------------------------------
! A routine for outputting fields for an sse application
!-------------------------------------------------------------------------------------------------
!
! @param iO Output file idx
! @param sQuadSc The velocity solution structure of the mesh
! @param maxlevel The maximum grid level used in computation (former NLMAX)
SUBROUTINE viz_CreateHistogram(field, mass,nData,minV,maxV,dCrit,bLog)
USE PP3D_MPI, ONLY:myid,showid,Comm_Summn
use Sigma_User, only : myOutput

implicit none

integer, intent(in) :: nData

logical, intent(in) :: bLog

real*8, intent(in) :: field(*), mass(*)

real*8, intent(inout) :: minV,maxV

real*8, intent(in) :: dCrit

real*8 :: d_low,d_high

! local variables
integer :: i,j,nBin
real*8 logShear,dMinCrit
real*8,allocatable :: HistoEta(:),BinEta(:,:)
real*8 Ml_i,dauxEta
real*8 minBinEta,mAXBinEta,dBinEta,meanEta

nBin = 100

if (.not.allocated(HistoEta)) allocate(HistoEta(nBin))
if (.not.allocated(BinEta))   allocate(BinEta(3,nBin))

IF (myid.ne.0) THEN

HistoEta = 0d0

if (bLog) then
 minBinEta = log10(minV)
 maxBinEta = log10(maxV)
else
 minBinEta = minV
 maxBinEta = maxV
end if

BinEta(1,   1)=minBinEta
BinEta(2,nBin)=maxBinEta

dBinEta = (maxBinEta-minBinEta)/DBLE(nBin-2)

DO i=2,nBin

 BinEta(2,i-1)= minBinEta + DBLE(i-2)*dBinEta
 BinEta(3,i-1)= minBinEta + (DBLE(i)-2.5d0)*dBinEta
 BinEta(1,i)  = BinEta(2,i-1)

END DO
BinEta(3,nBin)= minBinEta + (DBLE(nBin)-1.5d0)*dBinEta


DO i=1,nData

if (bLog) then
  logShear = LOG10(field(i))
else
  logShear = field(i)
end if

  Ml_i = mass(i)

  DO j=1,nBin
   IF (logShear.ge.BinEta(1,j).and.logShear.lt.BinEta(2,j)) HistoEta(j)  = HistoEta(j)  + Ml_i
  END DO

END DO
END IF

CALL COMM_SUMMn(HistoEta,nBin)

dauxEta = 0d0
DO i=1,nBin
 dauxEta = dauxEta + HistoEta(i)
END DO

dMinCrit = dauxEta*dCrit

d_low  = BinEta(1,1)
d_high = BinEta(2,nBin)

dauxEta = 0d0
DO i=1,nBin
 dauxEta = dauxEta + HistoEta(i)
 if (dauxEta.gt.dMinCrit) then
  if (bLog) then
   minV = 10d0**BinEta(1,i)
  else
   minV = BinEta(1,i)
  end if
  exit
 end if
END DO

dauxEta = 0d0
DO i=nBin,1,-1
 dauxEta = dauxEta + HistoEta(i)
 if (dauxEta.gt.dMinCrit) then
  if (bLog) then
   maxV = 10d0**BinEta(2,i)
  else
   maxV = BinEta(2,i)
  end if
  exit
 end if
END DO

end subroutine viz_CreateHistogram
!
!-------------------------------------------------------------------------------------------------
! A wrapper routine for outputting 1D fields for an sse application
!-------------------------------------------------------------------------------------------------
! @param iO Output file idx
! @param sQuadSc The velocity solution structure of the mesh
! @param sLinSc The pressure solution structure of the mesh
! @param sTracer Scalar solution structure
! @param maxlevel The maximum grid level used in computation (former NLMAX)
! @param ptempSim Optional Parameter: Wether it is a temperature-simulation or not
subroutine viz_OutPut_1D(iOut, sQuadSc, sLinSc, sTracer, maxlevel,btempSim)

USE PP3D_MPI, ONLY:myid
USE def_FEAT
use def_LinScalar, only: lScalar
USE var_QuadScalar, ONLY:ShearRate,Viscosity, my1DOut, my1Dout_nol,ScrewDist
use Sigma_User, only: myOutput, mySigma
use iniparser
use iso_c_binding, only: C_CHAR, C_NULL_CHAR

implicit none

interface
  subroutine c_write_json(thetype, dataName) bind(C, name="c_write_json")
    use cinterface, only: c1dOutput
    use iso_c_binding
    type(c1dOutput) :: thetype
    character(kind=c_char) :: dataName(*)
  end subroutine c_write_json
end interface

interface
  subroutine c_add_json_array(thetype, dataName) bind(C, name="c_add_json_array")
    use cinterface, only: c1dOutput
    use iso_c_binding
    type(c1dOutput) :: thetype
    character(kind=c_char) :: dataName(*)
  end subroutine c_add_json_array 
end interface

interface
  subroutine c_init_json_output(thetype) bind(C, name="c_init_json_output")
    use cinterface, only: c1dOutput
    use iso_c_binding
    type(c1dOutput) :: thetype
  end subroutine c_init_json_output
end interface

type(tQuadScalar), intent(in) :: sQuadSc

type(tLinScalar), intent(in) :: sLinSc

type(lScalar), intent(in) :: sTracer

integer, intent(in) :: maxlevel

logical, intent(in), optional :: btempSim

integer :: iOut

integer i,j
character cf*17,cf2*30
character(8)  :: cdate
character(10) :: ctime
character(5)  :: czone
integer,dimension(8) :: values
character(100) :: command
! logical :: bTemperatureSimulation
integer :: iVerlaufMax
integer :: ifile
logical :: bfileExists

type(c1dOutput) :: thestruct

! For dumping the e3dfile
type(t_parlist) :: parameterlistModifiedKTP1D

! The layout of the t1DOutput structure
type t1dname
 character(kind = c_char, len=20) :: name
 character(kind = c_char, len=20) :: unit
end type t1dname

type(t1dname) :: names1d(11)
 call cf_make_c_char(C_CHAR_"VelocityZ"//C_NULL_CHAR, names1d(1)%name )
 call cf_make_c_char(C_CHAR_"Pressure"//C_NULL_CHAR, names1d(2)%name )
 call cf_make_c_char(C_CHAR_"Temperature"//C_NULL_CHAR, names1d(3)%name )
 call cf_make_c_char(C_CHAR_"Viscosity"//C_NULL_CHAR, names1d(4)%name )
 call cf_make_c_char(C_CHAR_"ShearRate"//C_NULL_CHAR, names1d(5)%name )
 call cf_make_c_char(C_CHAR_"VelocityX"//C_NULL_CHAR, names1d(6)%name )
 call cf_make_c_char(C_CHAR_"VelocityY"//C_NULL_CHAR, names1d(7)%name )
 call cf_make_c_char(C_CHAR_"VelocityM"//C_NULL_CHAR, names1d(8)%name )
 call cf_make_c_char(C_CHAR_"ShearSG"//C_NULL_CHAR, names1d(9)%name )
 call cf_make_c_char(C_CHAR_"ShearSS"//C_NULL_CHAR, names1d(10)%name )
 call cf_make_c_char(C_CHAR_"VelocityWSS"//C_NULL_CHAR, names1d(11)%name )

 call cf_make_c_char(C_CHAR_'m/s'//C_NULL_CHAR, names1d(1)%unit )
 call cf_make_c_char(C_CHAR_'bar'//C_NULL_CHAR, names1d(2)%unit )
 call cf_make_c_char(C_CHAR_'C'//C_NULL_CHAR, names1d(3)%unit )
 call cf_make_c_char(C_CHAR_'kg/m/s'//C_NULL_CHAR, names1d(4)%unit )
 call cf_make_c_char(C_CHAR_'1/s'//C_NULL_CHAR, names1d(5)%unit )
 call cf_make_c_char(C_CHAR_'m/s'//C_NULL_CHAR, names1d(6)%unit )
 call cf_make_c_char(C_CHAR_'m/s'//C_NULL_CHAR, names1d(7)%unit )
 call cf_make_c_char(C_CHAR_'m/s'//C_NULL_CHAR, names1d(8)%unit )
 call cf_make_c_char(C_CHAR_'1/s'//C_NULL_CHAR, names1d(9)%unit )
 call cf_make_c_char(C_CHAR_'1/s'//C_NULL_CHAR, names1d(10)%unit )
 call cf_make_c_char(C_CHAR_'m/s'//C_NULL_CHAR, names1d(11)%unit )

! The layout of the t1DOutput structure
!TYPE t1DOutput
! REAL*8, ALLOCATABLE :: dMean(:),dMin(:),dMax(:),dLoc(:)
! CHARACTER cName*20
!END TYPE t1DOutput
!TYPE(t1DOutput) :: my1DOut(11)
!REAL*8, ALLOCATABLE :: my1DIntervals(:,:),my1DWeight(:)
!INTEGER my1DOut_nol

!  if (present(btempSim)) then
!    bTemperatureSimulation = btempSim
!  else
!    bTemperatureSimulation = .FALSE.
!  end if
!
!  ! The number of Verlauf-Outputs in the 1D-Files depends on wether it is a temperature-simulation
!  if (bTemperatureSimulation) then
!   ! If it is a Temperature-Simulation we have 8 Fields:
!   ! Pressure, Velocity_Mag, AxialVelocity, RotX-Velocity, RotY-Velocity, Shearrate, Viscosity, Temperature
!   ! For each field we output 3 Quantites: Min, Max, Med.
!   ! Therefore, we then have 3*8=24 Verlauf-Sections
!   iVerlaufMax = 24
!  else
!   ! If it is not a Temperature-Simulation we have only 7 Fields:
!   ! Pressure, Velocity_Mag, AxialVelocity, RotX-Velocity, RotY-Velocity, Shearrate, Viscosity
!   ! For each field we output 3 Quantites: Min, Max, Med.
!   ! Therefore, we then have 3*7=21 Verlauf-Sections
!   iVerlaufMax = 21
!  end if


 iVerlaufMax = 33

 my1DOut(1)%cName = 'VelocityZ_[m/s]'
 CALL viz_OutPut_1D_sub(sQuadSc%valW,sQuadSc%valV,sQuadSc%valU,1, my1Dout_nol, my1DOut, myOutput, mySigma, sQuadSc, maxlevel)

 my1DOut(2)%cName = 'Pressure_[bar]'
 CALL viz_OutPut_1D_sub(sLinSc%Q2,sLinSc%Q2,sLinSc%Q2,2, my1Dout_nol, my1DOut, myOutput, mySigma, sQuadSc, maxlevel)

 my1DOut(3)%cName = 'Temperature_[C]'
 CALL viz_OutPut_1D_sub(Temperature,sLinSc%Q2,sLinSc%Q2,3, my1Dout_nol, my1DOut, myOutput, mySigma, sQuadSc, maxlevel)

 my1DOut(4)%cName = 'Viscosity_[kg/m/s]'
 CALL viz_OutPut_1D_sub(Viscosity,sLinSc%Q2,sLinSc%Q2,4, my1Dout_nol, my1DOut, myOutput, mySigma, sQuadSc, maxlevel)

 my1DOut(5)%cName = 'ShearRate_[1/s]'
 CALL viz_OutPut_1D_sub(ShearRate,sLinSc%Q2,sLinSc%Q2,5, my1Dout_nol, my1DOut, myOutput, mySigma, sQuadSc, maxlevel)

 my1DOut(6)%cName = 'VelocityX_[m/s]'
 CALL viz_OutPut_1D_sub(sQuadSc%valW,sQuadSc%valV,sQuadSc%valU,6, my1Dout_nol, my1DOut, myOutput, mySigma, sQuadSc, maxlevel)

 my1DOut(7)%cName = 'VelocityY_[m/s]'
 CALL viz_OutPut_1D_sub(sQuadSc%valW,sQuadSc%valV,sQuadSc%valU,7, my1Dout_nol, my1DOut, myOutput, mySigma, sQuadSc, maxlevel)

 my1DOut(8)%cName = 'VelocityM_[m/s]'
 CALL viz_OutPut_1D_sub(sQuadSc%valW,sQuadSc%valV,sQuadSc%valU,8, my1Dout_nol, my1DOut, myOutput, mySigma, sQuadSc, maxlevel)

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE") THEN
 my1DOut(9)%cName = 'ShearSG_[1/s]'
 CALL OutPut_1D_subExtra(ShearRate,ScrewDist,9, my1Dout_nol, my1DOut, myOutput, mySigma, sQuadSc, maxlevel)

 my1DOut(10)%cName = 'ShearSS_[1/s]'
 CALL OutPut_1D_subExtra(ShearRate,ScrewDist,10, my1Dout_nol, my1DOut, myOutput, mySigma, sQuadSc, maxlevel)

 my1DOut(11)%cName = 'VelocityWSS_[m/s]'
 CALL OutPut_1D_subExtra(sQuadSc%valW,ScrewDist,11, my1Dout_nol, my1DOut, myOutput, mySigma, sQuadSc, maxlevel)
END IF

 thestruct%length = my1DOut_nol
 thestruct%dmean = c_loc(my1DOut(1)%dMean)
 thestruct%dmin = c_loc(my1DOut(1)%dMin)
 thestruct%dmax = c_loc(my1DOut(1)%dMax)
 thestruct%dLoc = c_loc(my1DOut(1)%dLoc) 

 call cf_make_c_char(C_CHAR_"cm"//C_NULL_CHAR, thestruct%unit_name)

 call c_init_json_output(thestruct) 

 IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE") THEN
  do i=1,11
    thestruct%dmean = c_loc(my1DOut(i)%dMean)
    thestruct%dmin = c_loc(my1DOut(i)%dMin)
    thestruct%dmax = c_loc(my1DOut(i)%dMax)
    thestruct%dLoc = c_loc(my1DOut(i)%dLoc) 
    call cf_make_c_char(names1d(i)%unit, thestruct%unit_name)
    call c_add_json_array(thestruct, names1d(i)%name)
  end do
 ELSE
  do i=1,8
    thestruct%dmean = c_loc(my1DOut(i)%dMean)
    thestruct%dmin = c_loc(my1DOut(i)%dMin)
    thestruct%dmax = c_loc(my1DOut(i)%dMax)
    thestruct%dLoc = c_loc(my1DOut(i)%dLoc) 
    call cf_make_c_char(names1d(i)%unit, thestruct%unit_name)
    call c_add_json_array(thestruct, names1d(i)%name)
  end do
 END IF

 call inip_init(parameterlistModifiedKTP1D)
 ! write(*,*)'------------------This should crash-----------------------------'
IF (myid.eq.1) THEN

 ! We need to dump the e3dfile to the 1d-files, so lets read it in
 call inip_readfromfile(parameterlistModifiedKTP1D,trim(adjustl("_data/Extrud3D.dat")))
 ! KPT needs them modified, all sections should be indeted by one level
 IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE") call INIP_indentAllSections(parameterlistModifiedKTP1D,"InputSigmaFile")
 IF (ADJUSTL(TRIM(mySigma%cType)).EQ."SSE") call INIP_indentAllSections(parameterlistModifiedKTP1D,"InputREXFile")
 IF (ADJUSTL(TRIM(mySigma%cType)).EQ."DIE") call INIP_indentAllSections(parameterlistModifiedKTP1D,"InputDIEFile")
 

 WRITE(cf2,'(A,I4.4,A)') '_1D/extrud3d_',iOut,'.res'

 call inip_openFileForWriting(cf2, ifile, INIP_REPLACE, bfileExists, .TRUE.)
! OPEN(UNIT=120,FILE=TRIM(ADJUSTL(cf2)))
 WRITE(ifile,'(A)')"[SigmaFileInfo]"
 WRITE(ifile,'(A)')"FileType=ResultsExtrud3d"
 call date_and_time(cdate,ctime,czone,values)
 WRITE(ifile,'(8A)')"Date=",cdate(7:8),"/",cdate(5:6),"/",cdate(3:4)
 WRITE(ifile,'(A)')"Extrud3dVersion=Extrud3d 2.0"
 WRITE(ifile,'(A,I2.2)')"counter_pos=",my1DOut_nol
 WRITE(ifile,'(A,I2.2)')"counter_verl=",iVerlaufMax
!  WRITE(120,'(A,E12.4)')"TimeLevel=",timens
 WRITE(ifile,'(A)') "[InputSigmaFile]"
 WRITE(ifile,'(("#")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-"))')
! CLOSE(120)
 call inip_dumpToUnit(parameterlistModifiedKTP1D,ifile)

! command = ' '
! command = "cat _data/Extrud3D.dat >> "//TRIM(ADJUSTL(cf2))
! CALL system(TRIM(ADJUSTL(command)))

 ! Write the 1D-Header.
 ! It depends on iVerlaufMax as this parameter tells us how many sections we
 ! are going to find.
!  call write_1d_header(iOut,my1DOut_nol,iVerlaufMax)
!  OPEN(UNIT=120,FILE=TRIM(ADJUSTL(cf2)),ACCESS='APPEND')

 WRITE(ifile,'(("#")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-")("-COPY-"))')

 WRITE(ifile,'(A)')"[Positions]"
 WRITE(ifile,'(A)')"ID=ELEMENT_LENGTH"
 WRITE(ifile,'(A)')"Unit=0"
 WRITE(ifile,'(A)')"SIUnit=mm"
 DO i=0,my1DOut_nol-1
!  WRITE(ifile,'(A,I2.2,A,E14.6,1X)') "POS",i,"=",1d1*my1DOut(1)%dLoc(i+1)
  CALL outputLine('POS',1d1*my1DOut(1)%dLoc(i+1))
 END DO
 WRITE(ifile,'(A)')"[Ergebnisse]"
! Pressure
 WRITE(ifile,'(A,A,(100("/")))')"#///////////","PRESSURE"
 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf0]"
 WRITE(ifile,'(A)')"ID=PRESSURE_MIN"
 WRITE(ifile,'(A)')"Unit=38"
 WRITE(ifile,'(A)')"SIUnit=bar"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(2)%dMin(i+1)-dMinOutputPressure)
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf1]"
 WRITE(ifile,'(A)')"ID=PRESSURE_MAX"
 WRITE(ifile,'(A)')"Unit=38"
 WRITE(ifile,'(A)')"SIUnit=bar"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(2)%dMax(i+1)-dMinOutputPressure)
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf2]"
 WRITE(ifile,'(A)')"ID=PRESSURE_MED"
 WRITE(ifile,'(A)')"Unit=38"
 WRITE(ifile,'(A)')"SIUnit=bar"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(2)%dMean(i+1)-dMinOutputPressure)
 END DO

! Shearrate
 WRITE(ifile,'(A,A,(100("/")))')"#///////////","SHEARRATE"
 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf3]"
 WRITE(ifile,'(A)')"ID=GAMMA_P_MIN"
 WRITE(ifile,'(A)')"Unit=17"
 WRITE(ifile,'(A)')"SIUnit=1/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(5)%dMin(i+1))
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf4]"
 WRITE(ifile,'(A)')"ID=GAMMA_P_MAX"
 WRITE(ifile,'(A)')"Unit=17"
 WRITE(ifile,'(A)')"SIUnit=1/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(5)%dMax(i+1))
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf5]"
 WRITE(ifile,'(A)')"ID=GAMMA_P_MED"
 WRITE(ifile,'(A)')"Unit=17"
 WRITE(ifile,'(A)')"SIUnit=1/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(5)%dMean(i+1))
 END DO

! Viscosity
 WRITE(ifile,'(A,A,(100("/")))')"#///////////","VISCOSITY"
 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf6]"
 WRITE(ifile,'(A)')"ID=ETA_MIN"
 WRITE(ifile,'(A)')"Unit=35"
 WRITE(ifile,'(A)')"SIUnit=Pa s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(4)%dMin(i+1))
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf7]"
 WRITE(ifile,'(A)')"ID=ETA_MAX"
 WRITE(ifile,'(A)')"Unit=35"
 WRITE(ifile,'(A)')"SIUnit=Pa s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(4)%dMax(i+1))
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf8]"
 WRITE(ifile,'(A)')"ID=ETA_MED"
 WRITE(ifile,'(A)')"Unit=35"
 WRITE(ifile,'(A)')"SIUnit=Pa s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(4)%dMean(i+1))
 END DO

! AxialVelocity
 WRITE(ifile,'(A,A,(100("/")))')"#///////////","AXIALVELOCITY"
 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf9]"
 WRITE(ifile,'(A)')"ID=AX_V_MIN"
 WRITE(ifile,'(A)')"Unit=20"
 WRITE(ifile,'(A)')"SIUnit=m/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(1)%dMin(i+1))
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf10]"
 WRITE(ifile,'(A)')"ID=AX_V_MAX"
 WRITE(ifile,'(A)')"Unit=20"
 WRITE(ifile,'(A)')"SIUnit=m/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(1)%dMax(i+1))
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf11]"
 WRITE(ifile,'(A)')"ID=AX_V_MED"
 WRITE(ifile,'(A)')"Unit=20"
 WRITE(ifile,'(A)')"SIUnit=m/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(1)%dMean(i+1))
 END DO

! Velocity-X
 WRITE(ifile,'(A,A,(100("/")))')"#///////////","VELOCITY_XCOMP"
 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf12]"
 WRITE(ifile,'(A)')"ID=ROTX_V_MIN"
 WRITE(ifile,'(A)')"Unit=20"
 WRITE(ifile,'(A)')"SIUnit=m/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(6)%dMin(i+1))
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf13]"
 WRITE(ifile,'(A)')"ID=ROTX_V_MAX"
 WRITE(ifile,'(A)')"Unit=20"
 WRITE(ifile,'(A)')"SIUnit=m/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(6)%dMax(i+1))
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf14]"
 WRITE(ifile,'(A)')"ID=ROTX_V_MED"
 WRITE(ifile,'(A)')"Unit=20"
 WRITE(ifile,'(A)')"SIUnit=m/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(6)%dMean(i+1))
 END DO

! Velocity-Y
 WRITE(ifile,'(A,A,(100("/")))')"#///////////","VELOCITY_YCOMP"
 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf15]"
 WRITE(ifile,'(A)')"ID=ROTY_V_MIN"
 WRITE(ifile,'(A)')"Unit=20"
 WRITE(ifile,'(A)')"SIUnit=m/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(7)%dMin(i+1))
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf16]"
 WRITE(ifile,'(A)')"ID=ROTY_V_MAX"
 WRITE(ifile,'(A)')"Unit=20"
 WRITE(ifile,'(A)')"SIUnit=m/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(7)%dMax(i+1))
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf17]"
 WRITE(ifile,'(A)')"ID=ROTY_V_MED"
 WRITE(ifile,'(A)')"Unit=20"
 WRITE(ifile,'(A)')"SIUnit=m/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(7)%dMean(i+1))
 END DO

! Velocity-Mag
 WRITE(ifile,'(A,A,(100("/")))')"#///////////","VELOCITY_MAG"
 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf18]"
 WRITE(ifile,'(A)')"ID=MAG_MIN"
 WRITE(ifile,'(A)')"Unit=20"
 WRITE(ifile,'(A)')"SIUnit=m/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(8)%dMin(i+1))
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf19]"
 WRITE(ifile,'(A)')"ID=MAG_MAX"
 WRITE(ifile,'(A)')"Unit=20"
 WRITE(ifile,'(A)')"SIUnit=m/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(8)%dMax(i+1))
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf20]"
 WRITE(ifile,'(A)')"ID=MAG_MED"
 WRITE(ifile,'(A)')"Unit=20"
 WRITE(ifile,'(A)')"SIUnit=m/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(8)%dMean(i+1))
 END DO

! Temperature
 WRITE(ifile,'(A,A,(100("/")))')"#///////////","TEMPERATURE"
 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf21]"
 WRITE(ifile,'(A)')"ID=TEMPERATURE_MIN"
 WRITE(ifile,'(A)')"Unit=23"
 WRITE(ifile,'(A)')"SIUnit=C"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(3)%dMin(i+1))
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf22]"
 WRITE(ifile,'(A)')"ID=TEMPERATURE_MAX"
 WRITE(ifile,'(A)')"Unit=23"
 WRITE(ifile,'(A)')"SIUnit=C"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(3)%dMax(i+1))
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf23]"
 WRITE(ifile,'(A)')"ID=TEMPERATURE_MED"
 WRITE(ifile,'(A)')"Unit=23"
 WRITE(ifile,'(A)')"SIUnit=C"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(3)%dMean(i+1))
 END DO

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE") THEN
! ShearSG
 WRITE(ifile,'(A,A,(100("/")))')"#///////////","SHEAR_SG"
 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf24]"
 WRITE(ifile,'(A)')"ID=GAMMA_SG_MIN"
 WRITE(ifile,'(A)')"Unit=17"
 WRITE(ifile,'(A)')"SIUnit=1/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(9)%dMin(i+1))
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf25]"
 WRITE(ifile,'(A)')"ID=GAMMA_SG_MAX"
 WRITE(ifile,'(A)')"Unit=17"
 WRITE(ifile,'(A)')"SIUnit=1/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(9)%dMax(i+1))
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf26]"
 WRITE(ifile,'(A)')"ID=GAMMA_SG_MED"
 WRITE(ifile,'(A)')"Unit=17"
 WRITE(ifile,'(A)')"SIUnit=1/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(9)%dMean(i+1))
 END DO

! ShearSS
 WRITE(ifile,'(A,A,(100("/")))')"#///////////","SHEAR_SS"
 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf27]"
 WRITE(ifile,'(A)')"ID=GAMMA_SS_MIN"
 WRITE(ifile,'(A)')"Unit=17"
 WRITE(ifile,'(A)')"SIUnit=1/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(10)%dMin(i+1))
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf28]"
 WRITE(ifile,'(A)')"ID=GAMMA_SS_MAX"
 WRITE(ifile,'(A)')"Unit=17"
 WRITE(ifile,'(A)')"SIUnit=1/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(10)%dMax(i+1))
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf29]"
 WRITE(ifile,'(A)')"ID=GAMMA_SS_MED"
 WRITE(ifile,'(A)')"Unit=17"
 WRITE(ifile,'(A)')"SIUnit=1/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(10)%dMean(i+1))
 END DO

! Velocity-W-SS
 WRITE(ifile,'(A,A,(100("/")))')"#///////////","VELO_Z_SS"
 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf30]"
 WRITE(ifile,'(A)')"ID=VELO_Z_SS_MIN"
 WRITE(ifile,'(A)')"Unit=20"
 WRITE(ifile,'(A)')"SIUnit=m/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(11)%dMin(i+1))
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf31]"
 WRITE(ifile,'(A)')"ID=VELO_Z_SS_MAX"
 WRITE(ifile,'(A)')"Unit=20"
 WRITE(ifile,'(A)')"SIUnit=m/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(11)%dMax(i+1))
 END DO

 WRITE(ifile,'(A)')"[Ergebnisse/Verlauf32]"
 WRITE(ifile,'(A)')"ID=VELO_Z_SS_MED"
 WRITE(ifile,'(A)')"Unit=20"
 WRITE(ifile,'(A)')"SIUnit=m/s"
 DO i=0,my1DOut_nol-1
  CALL outputLine('ST',my1DOut(11)%dMean(i+1))
 END DO

END IF

 CLOSE(ifile)


end if
call inip_done(parameterlistModifiedKTP1D)

contains
SUBROUTINE outputLine(c,dV)
character*(*) c
REAL*8 dV

IF (abs(dV).ge.1d29) THEN
 IF     (i.ge.0     .and.i.lt.    10) THEN
  WRITE(ifile,'(A,I1.1,A,1X)') c,i,"= _INVAILD_"
 ELSEIF (i.ge.10    .and.i.lt.   100) THEN
  WRITE(ifile,'(A,I2.2,A,1X)') c,i,"= _INVALID_"
 ELSEIF (i.ge.100   .and.i.lt.  1000) THEN
  WRITE(ifile,'(A,I3.3,A,1X)') c,i,"= _INVALID_"
 ELSEIF (i.ge.1000  .and.i.lt. 10000) THEN
  WRITE(ifile,'(A,I4.4,A,1X)') c,i,"= _INVALID_"
 ELSE

 ENDIF 
!  WRITE(ifile,'(A,I2.2,A,1X)') c,i,"= _INVALID_"
ELSE
 IF     (i.ge.0     .and.i.lt.    10) THEN
  WRITE(ifile,'(A,I1.1,A,E14.6,1X)') c,i,"=",dV
 ELSEIF (i.ge.10    .and.i.lt.   100) THEN
  WRITE(ifile,'(A,I2.2,A,E14.6,1X)') c,i,"=",dV
 ELSEIF (i.ge.100   .and.i.lt.  1000) THEN
  WRITE(ifile,'(A,I3.3,A,E14.6,1X)') c,i,"=",dV
 ELSEIF (i.ge.1000  .and.i.lt. 10000) THEN
  WRITE(ifile,'(A,I4.4,A,E14.6,1X)') c,i,"=",dV
 ELSE

 ENDIF 
!  WRITE(ifile,'(A,I2.2,A,E14.6,1X)') c,i,"=",dV
END IF

END SUBROUTINE outputLine

end subroutine


!
!-------------------------------------------------------------------------------------------------
! The particular routine for outputting 1D fields for an sse application
!-------------------------------------------------------------------------------------------------
! @param iO Output file idx
! @param sQuadSc The velocity solution structure of the mesh
! @param sLinSc The pressure solution structure of the mesh
! @param sTracer Scalar solution structure
! @param sTracer Scalar solution structure
subroutine viz_OutPut_1D_sub(dField1,dField2,dField3,i1D, my1DOut_nol, my1DOutput, myOutput, mySigma, sQuadSc, maxlevel)
real*8 :: dField1(*),dField2(*),dField3(*)
integer :: i1D
integer :: my1DOut_nol
type(t1DOutput), dimension(:) :: my1DOutput
type(tOutput) :: myOutput
type(tSigma) :: mySigma
type(tQuadScalar), intent(in) :: sQuadSc
integer, intent(in) :: maxlevel

type tHist
 real*8, allocatable :: x(:),m(:)
 integer n
end type tHist

type (tHist), allocatable :: myHist(:)

! local variables
real*8  :: dMinSample, dMaxSample
integer :: i,j,jj
real*8  :: dZ,dWidth,daux,dScale

real*8, dimension(:,:), allocatable :: my1DIntervals
real*8, dimension(:), allocatable :: my1DWeight

my1DOut_nol = myOutput%nOf1DLayers

dMinSample = 0d0
dMaxSample = mySigma%L

if (.not.allocated(my1DWeight)) ALLOCATE(my1DWeight(my1DOut_nol))
if (.not.allocated(my1DIntervals)) ALLOCATE(my1DIntervals(my1DOut_nol,2))
if (.not.allocated(my1DOutput(i1D)%dMean)) ALLOCATE(my1DOutput(i1D)%dMean(my1DOut_nol))
if (.not.allocated(my1DOutput(i1D)%dMin))  ALLOCATE(my1DOutput(i1D)%dMin(my1DOut_nol))
if (.not.allocated(my1DOutput(i1D)%dMax))  ALLOCATE(my1DOutput(i1D)%dMax(my1DOut_nol))
if (.not.allocated(my1DOutput(i1D)%dLoc))  ALLOCATE(my1DOutput(i1D)%dLoc(my1DOut_nol))

allocate(myHist(my1DOut_nol))
do i=1,my1DOut_nol
 allocate(myHist(i)%x(SIZE(sQuadSc%ValU)))
 allocate(myHist(i)%m(SIZE(sQuadSc%ValU)))
 myHist(i)%x = 0d0
 myHist(i)%m = 0d0
 myHist(i)%n = 0
end do

if (myid.ne.0) THEN

 if (i1D.EQ.1) dScale = 1d-2    !z-velo
 if (i1D.EQ.2) dScale = 1d-6    !Pressure
 if (i1D.EQ.3) dScale = 1d0     !Temperature
 if (i1D.EQ.4) dScale = 1d-1    !viscosity
 if (i1D.EQ.5) dScale = 1d0     !gamma
 if (i1D.EQ.6) dScale = 1d-2    !x-velo
 if (i1D.EQ.7) dScale = 1d-2    !y-velo
 if (i1D.EQ.8) dScale = 1d-2    !mag-velo

!  IF (i1D.EQ.1) dScale = 1d-2   !z-velo
!  IF (i1D.EQ.2) dScale = 1d-4   !Pressure
!  IF (i1D.EQ.3) dScale = 1d0    !Temperature
!  IF (i1D.EQ.4) dScale = 1d0    !gamma
!  IF (i1D.EQ.5) dScale = 1d-1   !viscosity
!  IF (i1D.EQ.6) dScale = 1d-2   !x-velo
!  IF (i1D.EQ.7) dScale = 1d-2   !y-velo
!  IF (i1D.EQ.8) dScale = 1d-2   !mag-velo
!
 dWidth = (dMaxSample-dMinSample)/DBLE(my1DOut_nol)

 DO i=1,my1DOut_nol
  my1DIntervals(i,1) = dMinSample + DBLE(i-1)*dWidth
  my1DIntervals(i,2) = dMinSample + DBLE(i)*dWidth
 END DO

 my1DOutput(i1D)%dMean   = 0d0
 my1DWeight           = 0d0
 my1DOutput(i1D)%dMin    = 1d30
 my1DOutput(i1D)%dMax    =-1d30
 MlRhoMat => mg_MlRhoMat(maxlevel)%a

 DO i=1,SIZE(sQuadSc%ValU)

  IF (MixerKNPR(i).eq.0) THEN
   jj=0
   dZ = mg_mesh%level(maxlevel)%dcorvg(3,i)
   DO j=1,my1DOut_nol
    IF (dZ.GE.my1DIntervals(j,1).AND.dZ.LE.my1DIntervals(j,2)) THEN
     jj = j
     EXIT
    END IF
   END DO

   IF (jj.NE.0) THEN
    myHist(jj)%n = myHist(jj)%n + 1
    daux = dScale*dField1(i)
    IF (i1D.EQ.6) daux = dScale*dField2(i)
    IF (i1D.EQ.7) daux = dScale*dField3(i)
    IF (i1D.EQ.8) daux = dScale*SQRT(dField1(i)**2d0+dField2(i)**2d0+dField3(i)**2d0)
    my1DOutput(i1D)%dMean(jj)   = my1DOutput(i1D)%dMean(jj)    + daux*MlRhoMat(i)
    my1DWeight(jj)              = my1DWeight(jj) + MlRhoMat(i)
    my1DOutput(i1D)%dMin(jj)    = MIN(my1DOutput(i1D)%dMin(jj),daux)
    my1DOutput(i1D)%dMax(jj)    = MAX(my1DOutput(i1D)%dMax(jj),daux)
    myHist(jj)%x(myHist(jj)%n) = daux
    myHist(jj)%m(myHist(jj)%n) = MlRhoMat(i)
   END IF
  END IF
 END DO
END IF

CALL COMM_SUMMN(my1DOutput(i1D)%dMean,my1DOut_nol)
CALL COMM_SUMMN(my1DWeight,my1DOut_nol)
CALL COMM_Maximumn(my1DOutput(i1D)%dMax,my1DOut_nol)
CALL COMM_Minimumn(my1DOutput(i1D)%dMin,my1DOut_nol)

DO j=1,my1DOut_nol
 if (i1D.eq.5) then
  CALL viz_CreateHistogram(myHist(j)%x, myHist(j)%m, myHist(j)%n,my1DOutput(i1D)%dMin(j),my1DOutput(i1D)%dMax(j),myOutput%CutDtata_1D,.true.)
 else
  CALL viz_CreateHistogram(myHist(j)%x, myHist(j)%m, myHist(j)%n,my1DOutput(i1D)%dMin(j),my1DOutput(i1D)%dMax(j),myOutput%CutDtata_1D,.false.)
 end if
END DO

IF (myid.ne.0) THEN
 DO i=1,my1DOut_nol
  if (my1DWeight(i).eq.0d0) then
   my1DOutput(i1D)%dMean(i)=1d30
  else
   my1DOutput(i1D)%dMean(i) = my1DOutput(i1D)%dMean(i)/my1DWeight(i)
  end if
  my1DOutput(i1D)%dLoc(i)  = 0.5d0*(my1DIntervals(i,1)+my1DIntervals(i,2))
 END DO
END IF

deallocate(myHist)

IF (myid.ne.0) THEN
 if (i1D.eq.2) then
  dMinOutputPressure = my1DOutput(i1D)%dMean(1)
  DO j=1,my1DOut_nol
   dMinOutputPressure = MIN(dMinOutputPressure,my1DOutput(i1D)%dMean(j))
  END DO
  IF (myid.eq.1) THEN
   write(*,*) 'dMinOutputPressure: ',dMinOutputPressure,myid
!   write(*,*) 'dMinOutputPressureValues: ',my1DOutput(i1D)%dMean(:)
  end if
 end if
end if


end subroutine viz_OutPut_1D_sub
!
! ----------------------------------------------
!
SUBROUTINE  OutPut_1D_subEXTRA(dField1,ScrewDist,i1D, my1DOut_nol, my1DOutput, myOutput, mySigma, sQuadSc, maxlevel)
use, intrinsic :: ieee_arithmetic
implicit none
real*8 :: dField1(*),ScrewDist(2,*)
integer :: i1D
integer :: my1DOut_nol
type(t1DOutput), dimension(:) :: my1DOutput
type(tOutput) :: myOutput
type(tSigma) :: mySigma
type(tQuadScalar), intent(in) :: sQuadSc
integer, intent(in) :: maxlevel

type tHist
 real*8, allocatable :: x(:),m(:)
 integer n
end type tHist

type (tHist), allocatable :: myHist(:)

! local variables
real*8  :: dMinSample, dMaxSample
integer :: i,j,jj
real*8  :: dX,dY,dZ,dWidth,daux,dScale,dR,dDist,dRadius

real*8, dimension(:,:), allocatable :: my1DIntervals
real*8, dimension(:), allocatable :: my1DWeight

LOGICAL bValid
INTEGER iSeg

my1DOut_nol = myOutput%nOf1DLayers
dMinSample = 0d0
dMaxSample = mySigma%L

if (.not.allocated(my1DWeight)) ALLOCATE(my1DWeight(my1DOut_nol))
if (.not.allocated(my1DIntervals)) ALLOCATE(my1DIntervals(my1DOut_nol,2))
if (.not.allocated(my1DOutput(i1D)%dMean)) ALLOCATE(my1DOutput(i1D)%dMean(my1DOut_nol))
if (.not.allocated(my1DOutput(i1D)%dMin))  ALLOCATE(my1DOutput(i1D)%dMin(my1DOut_nol))
if (.not.allocated(my1DOutput(i1D)%dMax))  ALLOCATE(my1DOutput(i1D)%dMax(my1DOut_nol))
if (.not.allocated(my1DOutput(i1D)%dLoc))  ALLOCATE(my1DOutput(i1D)%dLoc(my1DOut_nol))

allocate(myHist(my1DOut_nol))
do i=1,my1DOut_nol
 allocate(myHist(i)%x(SIZE(sQuadSc%ValU)))
 allocate(myHist(i)%m(SIZE(sQuadSc%ValU)))
 myHist(i)%x = 0d0
 myHist(i)%m = 0d0
 myHist(i)%n = 0
end do

IF (myid.ne.0) THEN

 IF (i1D.EQ.9)  dScale = 1d0     !Shear SG
 IF (i1D.EQ.10) dScale = 1d0     !Shear SS
 IF (i1D.EQ.11) dScale = 1d-3    !VelozityZ SS

 dWidth = (dMaxSample-dMinSample)/DBLE(my1DOut_nol)

 DO i=1,my1DOut_nol
  my1DIntervals(i,1) = dMinSample + DBLE(i-1)*dWidth
  my1DIntervals(i,2) = dMinSample + DBLE(i)*dWidth
 END DO

 my1DOutput(i1D)%dMean   = 0d0
 my1DWeight              = 0d0
 my1DOutput(i1D)%dMin    = 1d30
 my1DOutput(i1D)%dMax    =-1d30
 MlRhoMat => mg_MlRhoMat(maxlevel)%a


 DO iSeg=1,mySigma%NumberOfSeg

 DO i=1,SIZE(sQuadSc%ValU)

  IF (MixerKNPR(i).eq.0) THEN
   jj=0
   dZ = mg_mesh%level(maxlevel)%dcorvg(3,i)
   DO j=1,my1DOut_nol
    IF (dZ.GE.my1DIntervals(j,1).AND.dZ.LE.my1DIntervals(j,2)) THEN
     jj = j
     EXIT
    END IF
   END DO

   IF (jj.NE.0) THEN
    bValid = .FALSE.
    IF (i1D.EQ.9) THEN
     dX = mg_mesh%level(maxlevel)%dcorvg(1,i)
     dY = mg_mesh%level(maxlevel)%dcorvg(2,i)
     dDist = ABS(min(ScrewDist(1,i),ScrewDist(2,i)))
!     dDist = ABS(Distamce(i))
     IF (dY.gt.0d0)  dR = SQRT(dX**2d0+(dY-mySigma%a/2d0)**2d0)
     IF (dY.le.0d0)  dR = SQRT(dX**2d0+(dY+mySigma%a/2d0)**2d0)
!      write(*,*) mySigma%mySegment(iSeg)%delta,mySigma%mySegment(iSeg)%Ds,dDist,dR
!      pause
     IF (dR.gt.0.5d0*mySigma%mySegment(iSeg)%Ds.and.dDist.le.mySigma%mySegment(iSeg)%delta) THEN
      daux = dScale*dField1(i)
      bValid = .TRUE.
     END IF
    END IF
    IF (i1D.EQ.10) THEN
     dX = mg_mesh%level(maxlevel)%dcorvg(1,i)
     dY = mg_mesh%level(maxlevel)%dcorvg(2,i)
     IF (ABS(ScrewDist(1,i)).lt.mySigma%mySegment(iSeg)%s.and.ABS(ScrewDist(2,i)).lt.mySigma%mySegment(iSeg)%s)  THEN
!      IF (dY.gt.-0.5d0*mySigma%a.and.dY.lt.+0.5d0*mySigma%a.and.&
!          dX.gt.-0.5d0*mySigma%s.and.dX.lt.+0.5d0*mySigma%s)  THEN
      daux = dScale*dField1(i)
      bValid = .TRUE.
     END IF
    END IF
    IF (i1D.EQ.11) THEN

     IF (ieee_is_finite(mySigma%DZz)) THEN
      dRadius = 0.5d0*mySigma%DZz
     ELSE 
      dRadius = SQRT((0.5*mySigma%Dz_out)**2d0 - (0.5*mySigma%a)**2d0) + mySigma%W
     END IF 
     
!      dRadius = ((0.5d0*mySigma%Dz)**2d0 - (0.5d0*mySigma%a)**2d0)**0.5d0 + mySigma%V
     dX = mg_mesh%level(maxlevel)%dcorvg(1,i)
     dY = mg_mesh%level(maxlevel)%dcorvg(2,i)
     dR = SQRT(dX**2d0+dY**2d0)
     
     IF (dR.lt.dRadius)  THEN
      daux = dScale*dField1(i)
      bValid = .TRUE.
     END IF
    END IF

    IF (bValid) THEN
     myHist(jj)%n = myHist(jj)%n + 1
     my1DOutput(i1D)%dMean(jj)   = my1DOutput(i1D)%dMean(jj)    + daux*MlRhoMat(i)
     my1DWeight(jj)              = my1DWeight(jj) + MlRhoMat(i)
     my1DOutput(i1D)%dMin(jj)    = MIN(my1DOutput(i1D)%dMin(jj),daux)
     my1DOutput(i1D)%dMax(jj)    = MAX(my1DOutput(i1D)%dMax(jj),daux)
     myHist(jj)%x(myHist(jj)%n) = daux
     myHist(jj)%m(myHist(jj)%n) = MlRhoMat(i)
    END IF
   END IF
  END IF
 END DO
 END DO
END IF

CALL COMM_SUMMN(my1DOutput(i1D)%dMean,my1DOut_nol)
CALL COMM_SUMMN(my1DWeight,my1DOut_nol)
CALL COMM_Maximumn(my1DOutput(i1D)%dMax,my1DOut_nol)
CALL COMM_Minimumn(my1DOutput(i1D)%dMin,my1DOut_nol)

DO j=1,my1DOut_nol
 if (i1D.eq.9.or.i1D.eq.10) then
  CALL viz_CreateHistogram(myHist(j)%x, myHist(j)%m, myHist(j)%n,my1DOutput(i1D)%dMin(j),my1DOutput(i1D)%dMax(j),myOutput%CutDtata_1D,.true.)
 else
  CALL viz_CreateHistogram(myHist(j)%x, myHist(j)%m, myHist(j)%n,my1DOutput(i1D)%dMin(j),my1DOutput(i1D)%dMax(j),myOutput%CutDtata_1D,.false.)
 end if
END DO

IF (myid.ne.0) THEN
 DO i=1,my1DOut_nol
  if (my1DWeight(i).eq.0d0) then
   my1DOutput(i1D)%dMean(i)=1d30
  else
   my1DOutput(i1D)%dMean(i) = my1DOutput(i1D)%dMean(i)/my1DWeight(i)
  end if
  my1DOutput(i1D)%dLoc(i)  = 0.5d0*(my1DIntervals(i,1)+my1DIntervals(i,2))
 END DO
END IF

deallocate(myHist)

END SUBROUTINE  OutPut_1D_subEXTRA
!
!-------------------------------------------------------------------------------------------------
! A wrapper routine for outputting 1d fields for an sse application
!-------------------------------------------------------------------------------------------------
! @param sExport An export structure
! @param iOutput The output idx of the visulation file
! @param sQuadSc The velocity solution structure of the mesh
! @param sLinSc The pressure solution structure of the mesh
! @param sTracer The solution of the scalar tracer equation
! @param visc The viscosity solution
! @param dist The distance solution
! @param shear The shear rate solution
! @param mgMesh The mesh that will be written out
!subroutine viz_output_fields(sExport, iOutput, sQuadSc, sLinSc, sTracer, visc, dist, shear, mgMesh)
subroutine viz_output_1D_fields(sExport, iOutput, sQuadSc, sLinSc, visc, dist, shear, mgMesh)

use var_QuadScalar, only:tExport

USE PP3D_MPI, ONLY:myid
USE def_FEAT

USE Transport_Q1,ONLY:Tracer

implicit none

interface
  subroutine c_write_json_output() bind(C, name="c_write_json_output")
    use cinterface, only: c1dOutput
    use iso_c_binding
  end subroutine c_write_json_output
end interface

type(tExport), intent(in) :: sExport

integer, intent(in) :: iOutput

type(tQuadScalar), intent(in) :: sQuadSc

type(tLinScalar), intent(in) :: sLinSc

real*8, dimension(:), intent(in) :: visc

real*8, dimension(:), intent(in) :: dist

real*8, dimension(:), intent(in) :: shear

type(tMultiMesh), intent(in) :: mgMesh

! locals
integer :: ioutput_lvl


 call viz_OutputHistogram(iOutput, sQuadSc, mgMesh%nlmax)

! call viz_OutputRegionHistogram(iOutput)

 call viz_OutPut_1D(iOutput, sQuadSc, sLinSc, Tracer, mgMesh%nlmax)

 call c_write_json_output() 

end subroutine viz_output_1D_fields

end module visualization_out
