module vizsualization_out

use var_QuadScalar, only:knvt,knet,knat,knel,tMultiMesh,tQuadScalar,tLinScalar 
use  PP3D_MPI, only:myid,showid,subnodes
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
! @param visc The viscosity solution  
! @param dist The distance solution  
! @param shear The shear rate solution  
! @param mgMesh The mesh that will be written out
subroutine viz_output_fields(sExport, iOutput, sQuadSc, sLinSc, visc, dist, shear, mgMesh)

use var_QuadScalar, only:tExport

implicit none

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

!type(fieldPtr), dimension(3) :: packed

if (sExport%Format .eq. "VTK") then

 if (myid.ne.0) then

  ioutput_lvl = sExport%Level

  call viz_write_vtu_process(iOutput,&
    mgMesh%level(ioutput_lvl)%dcorvg,&
    mgMesh%level(ioutput_lvl)%kvert, &
    sQuadSc, sLinSc, visc, dist, shear, ioutput_lvl,&
    mgMesh)
 else

  call viz_write_pvtu_main(ioutput_lvl)
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
subroutine viz_write_vtu_process(iO,dcoor,kvert, sQuadSc, sLinSc, visc, dist, shear,&
                                 ioutput_lvl, mgMesh)

use var_QuadScalar,only:myExport,MixerKnpr

implicit none

integer, intent(in) :: iO

real*8, dimension(:,:), intent(in) :: dcoor

integer, dimension(:,:), intent(in) :: kvert

type(tQuadScalar), intent(in) :: sQuadSc 

type(tLinScalar), intent(in) :: sLinSc 

real*8, dimension(:), intent(in) :: visc

real*8, dimension(:), intent(in) :: dist

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
   write(iunit, '(A,E16.7)')"        ",REAL(1e-5*0.1d0*sLinSc%Q2(ivt))
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
!   write(iunit, '(A,E16.7)')"        ",REAL(1d1*DistanceToSurface(1,ivt))
   write(iunit, '(A,E16.7)')"        ",REAL(1.0)
  end do
  write(iunit, *)"        </DataArray>"

 case('Screw')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Screw",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(1d1*dist(ivt))
  end do
  write(iunit, *)"        </DataArray>"

  case('Viscosity')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Viscosity [Pa s]",""" format=""ascii"">"

  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(0.1d0*visc(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('NormShearRate')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","NormShearRate [1/s]",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(shear(ivt))
  end do
  write(iunit, *)"        </DataArray>"

!  CASE('Monitor')
!  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Monitor",""" format=""ascii"">"
!  do ivt=1,NoOfVert
!   write(iunit, '(A,E16.7)')"        ",REAL(myALE%Monitor(ivt))
!  end do
!  write(iunit, *)"        </DataArray>"
!
!  CASE('FBM')
!  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","FBM",""" format=""ascii"">"
!  do ivt=1,NoOfVert
!   write(iunit, '(A,E16.7)')"        ",REAL(MixerKNPR(ivt))
!  end do
!  write(iunit, *)"        </DataArray>"

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
 CASE('Screw')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Screw","""/>"
 CASE('Viscosity')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Viscosity [Pa s]","""/>"
 CASE('NormShearRate')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","NormShearRate [1/s]","""/>"
!  CASE('Monitor')
!  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Monitor","""/>"
!
!  CASE('FBM')
!  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","FBM","""/>"


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
end module vizsualization_out
