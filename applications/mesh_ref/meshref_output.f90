Module MeshRefOutput

use MeshRefDef

 contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Output_VTK

IMPLICIT NONE
INTEGER nel,nvt
INTEGER ive,ivt,ioffset
INTEGER :: iunit=123
CHARACTER*(100) filename

nel = mg_mesh%level(ilev)%nel
nvt = mg_mesh%level(ilev)%nvt

write(*,*) 'nel=',nel

filename=" "
WRITE(filename(1:),'(A)') adjustl(trim(cOutputFolder))//"/mesh.vtu"

WRITE(*,'(104("="))') 
WRITE(*,*) "Outputting vtk file into ",filename

OPEN (UNIT=iunit,FILE=filename)

write(iunit, *)"<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">"
write(iunit, *)"  <UnstructuredGrid>"
write(iunit, *)"    <Piece NumberOfPoints=""",nvt,""" NumberOfCells=""",nel,""">"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the node field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <PointData>"

 write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","MarkerV",""" format=""ascii"">"
 do ivt=1,nvt
  write(iunit, '(A,I10)')"        ",markerV(ivt)
 end do
 write(iunit, *)"        </DataArray>"

 if (allocated(level)) then
  write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","level",""" format=""ascii"">"
  do ivt=1,nvt
   write(iunit, '(A,I10)')"        ",level(ivt)
  end do
  write(iunit, *)"        </DataArray>"
 end if

 write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","ID",""" format=""ascii"">"
 do ivt=1,nvt
  write(iunit, '(A,E16.7)')"        ",REAL(ivt)
 end do
 write(iunit, *)"        </DataArray>"

 write(iunit, '(A)')"    </PointData>"

write(iunit, '(A)')"      <Points>"
write(iunit, '(A)')"        <DataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3"" format=""ascii"" RangeMin=""0"" RangeMax=""1.0"">"
do ivt=1,nvt
 write(iunit,'(A10,3E16.7)')"          ",REAL(mg_mesh%level(ilev)%dcorvg(1,ivt)),REAL(mg_mesh%level(ilev)%dcorvg(2,ivt)),REAL(mg_mesh%level(ilev)%dcorvg(3,ivt))
end do
write(iunit, *)"        </DataArray>"
write(iunit, *)"      </Points>"


write(iunit, '(A)')"    <CellData>"

write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","MarkerE",""" format=""ascii"">"
do ivt=1,nel
 write(iunit, '(A,I10)')"        ",MarkerE(ivt)
end do
write(iunit, *)"        </DataArray>"

if (allocated(AreaIntensity)) then
 write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","AreaIntensity",""" format=""ascii"">"
 do ivt=1,nel
  write(iunit, '(A,3E16.7)')"        ",AreaIntensity(2,ivt)
 end do
 write(iunit, *)"        </DataArray>"
end if

write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","OrigID",""" format=""ascii"">"
do ivt=1,nel
 write(iunit, '(A,I10)')"        ",ivt
end do
write(iunit, *)"        </DataArray>"

write(iunit, '(A)')"    </CellData>"

write(iunit, *)"      <Cells>"
write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""connectivity"" format=""ascii"" RangeMin=""0"" RangeMax=""",nel-1,""">"
do ive=1,nel   
 write(iunit, '(8I10)')mg_mesh%level(ilev)%kvert(1,ive)-1,mg_mesh%level(ilev)%kvert(2,ive)-1,mg_mesh%level(ilev)%kvert(3,ive)-1,mg_mesh%level(ilev)%kvert(4,ive)-1,&
                       mg_mesh%level(ilev)%kvert(5,ive)-1,mg_mesh%level(ilev)%kvert(6,ive)-1,mg_mesh%level(ilev)%kvert(7,ive)-1,mg_mesh%level(ilev)%kvert(8,ive)-1
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""offsets"" format=""ascii"" RangeMin=""8"" RangeMax=""",8*nel,""">"
ioffset=nel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')ive*8,(ive+1)*8,(ive+2)*8,(ive+3)*8,(ive+4)*8,(ive+5)*8,(ive+6)*8,(ive+7)*8
end do

do ive=ioffset+1,nel
 write(iunit, '(I10)')ive*8
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A)')"        <DataArray type=""UInt8"" Name=""types"" format=""ascii"" RangeMin=""12"" RangeMax=""12"">"
ioffset=nel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')12,12,12,12,12,12,12,12
end do
do ive=ioffset+1,nel
 write(iunit, '(I10)')12
end do
write(iunit, '(A)')"        </DataArray>"
 
write(iunit, *)"      </Cells>"
write(iunit, *)"    </Piece>"
   
write(iunit, *)"  </UnstructuredGrid>"
write(iunit, *)"</VTKFile>"
close(iunit)

END SUBROUTINE Output_VTK
!----------------------------------------------------------
SUBROUTINE Output_RefVTK

IMPLICIT NONE
INTEGER nel,nvt
INTEGER i,j, iel,jel,ive,ivt,ioffset,n
integer nnel,nnvt
INTEGER :: iunit=123
CHARACTER*(100) filename


nel = mg_mesh%level(ilev)%nel
nvt = mg_mesh%level(ilev)%nvt

nnvt = 0
nnel = 0

DO iel=1,nel
 if (.not.allocated(myRF(iel)%myRF)) then
  nnel = nnel + myRF(iel)%nOfElem
  nnvt = nnvt + myRF(iel)%nOfVert
 else
  do jel =1,myRF(iel)%nOfElem
    nnel = nnel + myRF(iel)%myRF(jel)%nOfElem
    nnvt = nnvt + myRF(iel)%myRF(jel)%nOfVert
  end do
 end if
end do


filename=" "
WRITE(filename(1:),'(A)') adjustl(trim(cOutputFolder))//"/RefMesh.vtu"

WRITE(*,'(104("="))') 
WRITE(*,*) "Outputting vtk file into ",filename

OPEN (UNIT=iunit,FILE=filename)

write(iunit, *)"<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">"
write(iunit, *)"  <UnstructuredGrid>"
write(iunit, *)"    <Piece NumberOfPoints=""",nnvt,""" NumberOfCells=""",nnel,""">"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the node field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <PointData>"

!  write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","MarkerV",""" format=""ascii"">"
!  do ivt=1,nvt
!   write(iunit, '(A,I10)')"        ",markerV(ivt)
!  end do
!  write(iunit, *)"        </DataArray>"
! 
!  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","ID",""" format=""ascii"">"
!  do ivt=1,nvt
!   write(iunit, '(A,E16.7)')"        ",REAL(ivt)
!  end do
!  write(iunit, *)"        </DataArray>"

write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","knpr",""" format=""ascii"">"
do iel=1,nel
 if (.not.allocated(myRF(iel)%myRF)) then
  do i=1,myRF(iel)%nOfVert
    write(iunit,'(A10,I0)')"          ",myRF(iel)%knpr(i)
  end do
 else
  do jel=1,myRF(iel)%nOfElem
   do i=1,myRF(iel)%myRF(jel)%nOfVert
     write(iunit,'(A10,I0)')"          ",myRF(iel)%myRF(jel)%knpr(i)
   end do
  end do
 end if  
end do
write(iunit, *)"        </DataArray>"

write(iunit, '(A)')"    </PointData>"

write(iunit, '(A)')"      <Points>"
write(iunit, '(A)')"        <DataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3"" format=""ascii"" RangeMin=""0"" RangeMax=""1.0"">"
do iel=1,nel
 if (.not.allocated(myRF(iel)%myRF)) then
  do i=1,myRF(iel)%nOfVert
   write(iunit,'(A10,3E16.7)')"          ",REAL(myRF(iel)%dcoor(1,i)),REAL(myRF(iel)%dcoor(2,i)),REAL(myRF(iel)%dcoor(3,i))
  end do
 else
  do jel=1,myRF(iel)%nOfElem
   do i=1,myRF(iel)%myRF(jel)%nOfVert
    write(iunit,'(A10,3E16.7)')"          ",REAL(myRF(iel)%myRF(jel)%dcoor(1,i)),REAL(myRF(iel)%myRF(jel)%dcoor(2,i)),REAL(myRF(iel)%myRF(jel)%dcoor(3,i))
   end do
  end do
 end if
end do
write(iunit, *)"        </DataArray>"
write(iunit, *)"      </Points>"

write(iunit, '(A)')"    <CellData>"

write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","NREFEL",""" format=""ascii"">"
do iel=1,nel
 n=myRF(iel)%nOfElem
 if (.not.allocated(myRF(iel)%myRF)) then
  do i=1,myRF(iel)%nOfElem
   write(iunit, '(A,I10)')"        ",n !MarkerE(ivt)
  end do
 else
  do jel=1,myRF(iel)%nOfElem
   do i=1,myRF(iel)%myRF(jel)%nOfElem
    write(iunit, '(A,I10)')"        ",n !MarkerE(ivt)
   end do
  end do
 end if
end do
write(iunit, *)"        </DataArray>"

write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","PATCH_ID",""" format=""ascii"">"
do iel=1,nel
 n=myRF(iel)%nOfElem
 if (.not.allocated(myRF(iel)%myRF)) then
  do i=1,myRF(iel)%nOfElem
   write(iunit, '(A,I10)')"        ",myRF(iel)%patchID
  end do
 else
  do jel=1,myRF(iel)%nOfElem
   do i=1,myRF(iel)%myRF(jel)%nOfElem
    write(iunit, '(A,I10)')"        ",myRF(iel)%myRF(jel)%patchID
   end do
  end do
 end if
end do
write(iunit, *)"        </DataArray>"

write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","Selection",""" format=""ascii"">"
do iel=1,nel
 n=myRF(iel)%nOfElem
 if (.not.allocated(myRF(iel)%myRF)) then
  do i=1,myRF(iel)%nOfElem
   write(iunit, '(A,I10)')"        ",myRF(iel)%selection
  end do
 else
  do jel=1,myRF(iel)%nOfElem
   do i=1,myRF(iel)%myRF(jel)%nOfElem
    write(iunit, '(A,I10)')"        ",myRF(iel)%myRF(jel)%selection
   end do
  end do
 end if
end do
write(iunit, *)"        </DataArray>"

write(iunit, '(A)')"    </CellData>"

write(iunit, *)"      <Cells>"
write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""connectivity"" format=""ascii"" RangeMin=""0"" RangeMax=""",nel-1,""">"
nnvt = 0
do iel=1,nel
 if (.not.allocated(myRF(iel)%myRF)) then
  do i=1,myRF(iel)%nOfElem
   write(iunit, '(8I10)') nnvt + myRF(iel)%kvert(1,i) - 1, nnvt + myRF(iel)%kvert(2,i) - 1, nnvt + myRF(iel)%kvert(3,i) - 1, nnvt + myRF(iel)%kvert(4,i) - 1,&
                          nnvt + myRF(iel)%kvert(5,i) - 1, nnvt + myRF(iel)%kvert(6,i) - 1, nnvt + myRF(iel)%kvert(7,i) - 1, nnvt + myRF(iel)%kvert(8,i) - 1 
  end do
  nnvt = nnvt + myRF(iel)%nOfVert
 else
  do jel=1,myRF(iel)%nOfElem
   do i=1,myRF(iel)%myRF(jel)%nOfElem
    write(iunit, '(8I10)') nnvt + myRF(iel)%myRF(jel)%kvert(1,i) - 1, nnvt + myRF(iel)%myRF(jel)%kvert(2,i) - 1, nnvt + myRF(iel)%myRF(jel)%kvert(3,i) - 1, nnvt + myRF(iel)%myRF(jel)%kvert(4,i) - 1,&
                           nnvt + myRF(iel)%myRF(jel)%kvert(5,i) - 1, nnvt + myRF(iel)%myRF(jel)%kvert(6,i) - 1, nnvt + myRF(iel)%myRF(jel)%kvert(7,i) - 1, nnvt + myRF(iel)%myRF(jel)%kvert(8,i) - 1 
   end do
   nnvt = nnvt + myRF(iel)%myRF(jel)%nOfVert
  end do
 end if
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""offsets"" format=""ascii"" RangeMin=""8"" RangeMax=""",8*nel,""">"
ioffset=nnel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')ive*8,(ive+1)*8,(ive+2)*8,(ive+3)*8,(ive+4)*8,(ive+5)*8,(ive+6)*8,(ive+7)*8
end do

do ive=ioffset+1,nnel
 write(iunit, '(I10)')ive*8
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A)')"        <DataArray type=""UInt8"" Name=""types"" format=""ascii"" RangeMin=""12"" RangeMax=""12"">"
ioffset=nnel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')12,12,12,12,12,12,12,12
end do
do ive=ioffset+1,nnel
 write(iunit, '(I10)')12
end do
write(iunit, '(A)')"        </DataArray>"
 
write(iunit, *)"      </Cells>"
write(iunit, *)"    </Piece>"
   
write(iunit, *)"  </UnstructuredGrid>"
write(iunit, *)"</VTKFile>"
close(iunit)

END SUBROUTINE Output_RefVTK
!----------------------------------------------------------
SUBROUTINE Output_UniqueRefVTK

IMPLICIT NONE
INTEGER nel,nvt
INTEGER i,j, iel,jel,ive,ivt,ioffset,n
integer nnel,nnvt
INTEGER :: iunit=123
CHARACTER*(100) filename


nel = mg_mesh%level(ilev)%nel
nvt = mg_mesh%level(ilev)%nvt

nnvt = 0
nnel = 0

DO iel=1,nel
  nnel = nnel + myRF(iel)%nUniqueElems
  nnvt = nnvt + myRF(iel)%nUniquePoints
end do


filename=" "
WRITE(filename(1:),'(A)') adjustl(trim(cOutputFolder))//"/UniqueRefMesh.vtu"

WRITE(*,'(104("="))') 
WRITE(*,*) "Outputting vtk file into ",filename

OPEN (UNIT=iunit,FILE=filename)

write(iunit, *)"<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">"
write(iunit, *)"  <UnstructuredGrid>"
write(iunit, *)"    <Piece NumberOfPoints=""",nnvt,""" NumberOfCells=""",nnel,""">"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the node field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <PointData>"

!  write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","MarkerV",""" format=""ascii"">"
!  do ivt=1,nvt
!   write(iunit, '(A,I10)')"        ",markerV(ivt)
!  end do
!  write(iunit, *)"        </DataArray>"
! 
!  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","ID",""" format=""ascii"">"
!  do ivt=1,nvt
!   write(iunit, '(A,E16.7)')"        ",REAL(ivt)
!  end do
!  write(iunit, *)"        </DataArray>"

! write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","knpr",""" format=""ascii"">"
! do iel=1,nel
!  if (.not.allocated(myRF(iel)%myRF)) then
!   do i=1,myRF(iel)%nOfVert
!     write(iunit,'(A10,I0)')"          ",myRF(iel)%knpr(i)
!   end do
!  else
!   do jel=1,myRF(iel)%nOfElem
!    do i=1,myRF(iel)%myRF(jel)%nOfVert
!      write(iunit,'(A10,I0)')"          ",myRF(iel)%myRF(jel)%knpr(i)
!    end do
!   end do
!  end if  
! end do
! write(iunit, *)"        </DataArray>"

write(iunit, '(A)')"    </PointData>"

write(iunit, '(A)')"      <Points>"
write(iunit, '(A)')"        <DataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3"" format=""ascii"" RangeMin=""0"" RangeMax=""1.0"">"
do iel=1,nel
  do i=1,myRF(iel)%nUniquePoints
   write(iunit,'(A10,3E16.7)')"          ",REAL(myRF(iel)%dUniquedCoor(1,i)),REAL(myRF(iel)%dUniquedCoor(2,i)),REAL(myRF(iel)%dUniquedCoor(3,i))
  end do
end do
write(iunit, *)"        </DataArray>"
write(iunit, *)"      </Points>"

write(iunit, '(A)')"    <CellData>"

! write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","NREFEL",""" format=""ascii"">"
! do iel=1,nel
!  n=myRF(iel)%nOfElem
!  if (.not.allocated(myRF(iel)%myRF)) then
!   do i=1,myRF(iel)%nOfElem
!    write(iunit, '(A,I10)')"        ",n !MarkerE(ivt)
!   end do
!  else
!   do jel=1,myRF(iel)%nOfElem
!    do i=1,myRF(iel)%myRF(jel)%nOfElem
!     write(iunit, '(A,I10)')"        ",n !MarkerE(ivt)
!    end do
!   end do
!  end if
! end do
! write(iunit, *)"        </DataArray>"

write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","PATCH_ID",""" format=""ascii"">"
do iel=1,nel
  do i=1,myRF(iel)%nUniqueElems
   write(iunit, '(A,I10)')"        ",myRF(iel)%patchID
  end do
end do
write(iunit, *)"        </DataArray>"

! write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","Selection",""" format=""ascii"">"
! do iel=1,nel
!  n=myRF(iel)%nOfElem
!  if (.not.allocated(myRF(iel)%myRF)) then
!   do i=1,myRF(iel)%nOfElem
!    write(iunit, '(A,I10)')"        ",myRF(iel)%selection
!   end do
!  else
!   do jel=1,myRF(iel)%nOfElem
!    do i=1,myRF(iel)%myRF(jel)%nOfElem
!     write(iunit, '(A,I10)')"        ",myRF(iel)%myRF(jel)%selection
!    end do
!   end do
!  end if
! end do
! write(iunit, *)"        </DataArray>"

write(iunit, '(A)')"    </CellData>"

write(iunit, *)"      <Cells>"
write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""connectivity"" format=""ascii"" RangeMin=""0"" RangeMax=""",nel-1,""">"
nnvt = 0
do iel=1,nel
  do i=1,myRF(iel)%nUniqueElems
   write(iunit, '(8I10)') nnvt + myRF(iel)%kUniqueElem(1,i) - 1, nnvt + myRF(iel)%kUniqueElem(2,i) - 1, nnvt + myRF(iel)%kUniqueElem(3,i) - 1, nnvt + myRF(iel)%kUniqueElem(4,i) - 1,&
                          nnvt + myRF(iel)%kUniqueElem(5,i) - 1, nnvt + myRF(iel)%kUniqueElem(6,i) - 1, nnvt + myRF(iel)%kUniqueElem(7,i) - 1, nnvt + myRF(iel)%kUniqueElem(8,i) - 1 
  end do
  nnvt = nnvt + myRF(iel)%nUniquePoints
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""offsets"" format=""ascii"" RangeMin=""8"" RangeMax=""",8*nel,""">"
ioffset=nnel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')ive*8,(ive+1)*8,(ive+2)*8,(ive+3)*8,(ive+4)*8,(ive+5)*8,(ive+6)*8,(ive+7)*8
end do

do ive=ioffset+1,nnel
 write(iunit, '(I10)')ive*8
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A)')"        <DataArray type=""UInt8"" Name=""types"" format=""ascii"" RangeMin=""12"" RangeMax=""12"">"
ioffset=nnel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')12,12,12,12,12,12,12,12
end do
do ive=ioffset+1,nnel
 write(iunit, '(I10)')12
end do
write(iunit, '(A)')"        </DataArray>"
 
write(iunit, *)"      </Cells>"
write(iunit, *)"    </Piece>"
   
write(iunit, *)"  </UnstructuredGrid>"
write(iunit, *)"</VTKFile>"
close(iunit)

END SUBROUTINE Output_UniqueRefVTK
!----------------------------------------------------------
SUBROUTINE Output_ReducedRefVTK

IMPLICIT NONE
INTEGER nel,nvt
INTEGER i,j, iel,jel,ive,ivt,ioffset,n
integer nnel,nnvt
INTEGER :: iunit=123
CHARACTER*(100) filename

nnel = nReducedElems
nnvt = nReducedPoints

filename=" "
WRITE(filename(1:),'(A)') adjustl(trim(cOutputFolder))//"/ReducedRefMesh.vtu"

WRITE(*,'(104("="))') 
WRITE(*,*) "Outputting vtk file into ",filename

OPEN (UNIT=iunit,FILE=filename)

write(iunit, *)"<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">"
write(iunit, *)"  <UnstructuredGrid>"
write(iunit, *)"    <Piece NumberOfPoints=""",nnvt,""" NumberOfCells=""",nnel,""">"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the node field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <PointData>"

!  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","KNPR",""" format=""ascii"">"
!  do ivt=1,nnvt
!   write(iunit, '(A,E16.7)')"        ",REAL(MergedMeshKnpr(ivt))
!  end do
!  write(iunit, *)"        </DataArray>"

write(iunit, '(A)')"    </PointData>"
write(iunit, '(A)')"      <Points>"
write(iunit, '(A)')"        <DataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3"" format=""ascii"" RangeMin=""0"" RangeMax=""1.0"">"
do i=1,nnvt
   write(iunit,'(A10,3E16.7)')"          ",REAL(ReducedMeshCoor(1,i)),REAL(ReducedMeshCoor(2,i)),REAL(ReducedMeshCoor(3,i))
end do
write(iunit, *)"        </DataArray>"
write(iunit, *)"      </Points>"

write(iunit, '(A)')"    <CellData>"

write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","K",""" format=""ascii"">"
do iel=1,nnel
 write(iunit, '(A,E16.7)')"        ",REAL(iel)
end do
write(iunit, *)"        </DataArray>"

write(iunit, '(A)')"    </CellData>"

write(iunit, *)"      <Cells>"
write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""connectivity"" format=""ascii"" RangeMin=""0"" RangeMax=""",nel-1,""">"
do iel=1,nnel
   write(iunit, '(8I10)') ReducedMeshElem(1,iel) - 1,ReducedMeshElem(2,iel) - 1,ReducedMeshElem(3,iel) - 1,ReducedMeshElem(4,iel) - 1,&
                          ReducedMeshElem(5,iel) - 1,ReducedMeshElem(6,iel) - 1,ReducedMeshElem(7,iel) - 1,ReducedMeshElem(8,iel) - 1
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""offsets"" format=""ascii"" RangeMin=""8"" RangeMax=""",8*nel,""">"
ioffset=nnel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')ive*8,(ive+1)*8,(ive+2)*8,(ive+3)*8,(ive+4)*8,(ive+5)*8,(ive+6)*8,(ive+7)*8
end do

do ive=ioffset+1,nnel
 write(iunit, '(I10)')ive*8
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A)')"        <DataArray type=""UInt8"" Name=""types"" format=""ascii"" RangeMin=""12"" RangeMax=""12"">"
ioffset=nnel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')12,12,12,12,12,12,12,12
end do
do ive=ioffset+1,nnel
 write(iunit, '(I10)')12
end do
write(iunit, '(A)')"        </DataArray>"
 
write(iunit, *)"      </Cells>"
write(iunit, *)"    </Piece>"
   
write(iunit, *)"  </UnstructuredGrid>"
write(iunit, *)"</VTKFile>"
close(iunit)

END SUBROUTINE Output_ReducedRefVTK
! ----------------------------------------------
SUBROUTINE Output_ReducedCleanRefVTK
use Sigma_User, only : myProcess

IMPLICIT NONE
INTEGER nel,nvt
INTEGER i,j, iel,jel,ive,ivt,ioffset,n,iInflow
integer nnel,nnvt
INTEGER :: iunit=123
CHARACTER*(100) filename

nnel = nReducedCleanElems
nnvt = nReducedCleanPoints

filename=" "
WRITE(filename(1:),'(A)') adjustl(trim(cOutputFolder))//"/ReducedCleanRefMesh.vtu"

WRITE(*,'(104("="))') 
WRITE(*,*) "Outputting vtk file into ",filename

OPEN (UNIT=iunit,FILE=filename)

write(iunit, *)"<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">"
write(iunit, *)"  <UnstructuredGrid>"
write(iunit, *)"    <Piece NumberOfPoints=""",nnvt,""" NumberOfCells=""",nnel,""">"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the node field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <PointData>"

 write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Outflow",""" format=""ascii"">"
 do ivt=1,nnvt
  IF (ParCleanList%Outflow(ivt)) THEN
   write(iunit, '(A,E16.7)')"        ",REAL(1.0)
  else
   write(iunit, '(A,E16.7)')"        ",REAL(0.0)
  end if
 end do
 write(iunit, *)"        </DataArray>"

do iInflow=1,myProcess%nOfInflows
 write(iunit, '(A,A,I0,A)')"        <DataArray type=""Float32"" Name=""","Inflow",iInflow,""" format=""ascii"">"
 do ivt=1,nnvt
  IF (ParCleanList%Inflow(iInflow,ivt)) THEN
   write(iunit, '(A,E16.7)')"        ",REAL(1.0)
  else
   write(iunit, '(A,E16.7)')"        ",REAL(0.0)
  end if
 end do
 write(iunit, *)"        </DataArray>"
end do

write(iunit, '(A)')"    </PointData>"
write(iunit, '(A)')"      <Points>"
write(iunit, '(A)')"        <DataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3"" format=""ascii"" RangeMin=""0"" RangeMax=""1.0"">"
do i=1,nnvt
   write(iunit,'(A10,3E16.7)')"          ",REAL(ReducedCleanMeshCoor(1,i)),REAL(ReducedCleanMeshCoor(2,i)),REAL(ReducedCleanMeshCoor(3,i))
end do
write(iunit, *)"        </DataArray>"
write(iunit, *)"      </Points>"

write(iunit, '(A)')"    <CellData>"

write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","K",""" format=""ascii"">"
do iel=1,nnel
 write(iunit, '(A,E16.7)')"        ",REAL(iel)
end do
write(iunit, *)"        </DataArray>"

write(iunit, '(A)')"    </CellData>"

write(iunit, *)"      <Cells>"
write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""connectivity"" format=""ascii"" RangeMin=""0"" RangeMax=""",nel-1,""">"
do iel=1,nnel
   write(iunit, '(8I10)') ReducedCleanMeshElem(1,iel) - 1,ReducedCleanMeshElem(2,iel) - 1,ReducedCleanMeshElem(3,iel) - 1,ReducedCleanMeshElem(4,iel) - 1,&
                          ReducedCleanMeshElem(5,iel) - 1,ReducedCleanMeshElem(6,iel) - 1,ReducedCleanMeshElem(7,iel) - 1,ReducedCleanMeshElem(8,iel) - 1
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""offsets"" format=""ascii"" RangeMin=""8"" RangeMax=""",8*nel,""">"
ioffset=nnel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')ive*8,(ive+1)*8,(ive+2)*8,(ive+3)*8,(ive+4)*8,(ive+5)*8,(ive+6)*8,(ive+7)*8
end do

do ive=ioffset+1,nnel
 write(iunit, '(I10)')ive*8
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A)')"        <DataArray type=""UInt8"" Name=""types"" format=""ascii"" RangeMin=""12"" RangeMax=""12"">"
ioffset=nnel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')12,12,12,12,12,12,12,12
end do
do ive=ioffset+1,nnel
 write(iunit, '(I10)')12
end do
write(iunit, '(A)')"        </DataArray>"
 
write(iunit, *)"      </Cells>"
write(iunit, *)"    </Piece>"
   
write(iunit, *)"  </UnstructuredGrid>"
write(iunit, *)"</VTKFile>"
close(iunit)

END SUBROUTINE Output_ReducedCleanRefVTK
! ----------------------------------------------
SUBROUTINE Output_ReducedRefTRIandPAR
use Sigma_User, only : myProcess
IMPLICIT NONE
INTEGER nel,nvt
INTEGER i,j, iel,jel,ive,ivt,ioffset,n,iInflow,nBC,jat
integer nnel,nnvt
INTEGER :: iunit=123
CHARACTER*(256) cf

WRITE(cf,'(A)') ADJUSTL(TRIM(cOutputFolder))//'/'//ADJUSTL(TRIM(cReducedMeshdFolder))//'/'//adjustl(trim(cReducedGridFile))
WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cf))//"'"
OPEN(UNIT=iunit,FILE=ADJUSTL(TRIM(cf)))
WRITE(iunit,*) 'Coarse mesh exported by DeViSoR TRI3D exporter'
WRITE(iunit,*) 'Parametrisierung PARXC, PARYC, TMAXC'
WRITE(iunit,'(2I8,A)') nReducedCleanElems,nReducedCleanPoints, " 1 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE"

WRITE(iunit,'(A)') 'DCORVG'
DO i = 1,nReducedCleanPoints
 WRITE(iunit,'(3ES13.5)') MeshOutputScaleFactor*ReducedCleanMeshCoor(:,i)
END DO

WRITE(iunit,'(A)') 'KVERT'
DO i = 1,nReducedCleanElems
 WRITE(iunit,'(8I8)') ReducedCleanMeshElem(:,i)
END DO

WRITE(iunit,'(A)') 'KNPR'
DO i = 1,nReducedCleanPoints
  WRITE(iunit,'(I8)') 0
END DO

CLOSE(iunit)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
WRITE(cf,'(A)') ADJUSTL(TRIM(cOutputFolder))//'/'//ADJUSTL(TRIM(cReducedMeshdFolder))//'/Wall.par'
WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cf))//"'"
OPEN(UNIT=iunit,FILE=ADJUSTL(TRIM(cf)))
nBC = 0
DO i=1,nReducedCleanPoints
 if (ParCleanList%Wall(i)) nBC = nBC + 1
END DO
write(iunit,'(I0,A)') nBC,' '//'Wall'
write(iunit,'(A,I0," ",4F10.4,A)') "' '"
DO i=1,nReducedCleanPoints
 if (ParCleanList%Wall(i)) write(iunit,'(I0)') i
END DO
CLOSE(iunit)

WRITE(cf,'(A)') ADJUSTL(TRIM(cOutputFolder))//'/'//ADJUSTL(TRIM(cReducedMeshdFolder))//'/Wall.pls'
WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cf))//"'"
OPEN(UNIT=iunit,FILE=ADJUSTL(TRIM(cf)))
write(iunit,'(I0,A)') ParCleanList%nWA,' '//'Wall'
write(iunit,'(A,I0," ",4F10.4,A)') "' '"
DO i=1,ParCleanList%nWA
 write(iunit,'(4(I0,A))') ParCleanList%tWA%i(1,i),',',ParCleanList%tWA%i(2,i),',',ParCleanList%tWA%i(3,i),',',ParCleanList%tWA%i(4,i)
END DO
CLOSE(iunit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
DO jat=1,5
 WRITE(cf,'(A,I0,A)') ADJUSTL(TRIM(cOutputFolder))//'/'//ADJUSTL(TRIM(cReducedMeshdFolder))//'/SideWall',jat,'.par'
 WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cf))//"'"
 OPEN(UNIT=iunit,FILE=ADJUSTL(TRIM(cf)))
 nBC = 0
 DO i=1,nReducedCleanPoints
  if (ParCleanList%SideWall(jat,i)) nBC = nBC + 1
 END DO
 write(iunit,'(I0,A)') nBC,' '//'Wall'
 write(iunit,'(A,I0," ",4F10.4,A)') "' '"
 DO i=1,nReducedCleanPoints
  if (ParCleanList%SideWall(jat,i)) write(iunit,'(I0)') i
 END DO
 CLOSE(iunit)

 WRITE(cf,'(A,I0,A)') ADJUSTL(TRIM(cOutputFolder))//'/'//ADJUSTL(TRIM(cReducedMeshdFolder))//'/SideWall',jat,'.pls'
 WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cf))//"'"
 OPEN(UNIT=iunit,FILE=ADJUSTL(TRIM(cf)))
 write(iunit,'(I0,A)') ParCleanList%nSA(jat),' '//'Wall'
 write(iunit,'(A,I0," ",4F10.4,A)') "' '"
 DO i=1,ParCleanList%nSA(jat)
  write(iunit,'(4(I0,A))') ParCleanList%tSA(jat)%i(1,i),',',ParCleanList%tSA(jat)%i(2,i),',',ParCleanList%tSA(jat)%i(3,i),',',ParCleanList%tSA(jat)%i(4,i)
 END DO
 CLOSE(iunit)
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
WRITE(cf,'(A)') ADJUSTL(TRIM(cOutputFolder))//'/'//ADJUSTL(TRIM(cReducedMeshdFolder))//'/Outflow.par'
WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cf))//"'"
OPEN(UNIT=iunit,FILE=ADJUSTL(TRIM(cf)))
nBC = 0
DO i=1,nReducedCleanPoints
 if (ParCleanList%Outflow(i)) nBC = nBC + 1
END DO
write(iunit,'(I0,A)') nBC,' '//'Outflow'
write(iunit,'(A,I0," ",4F10.4,A)') "' '"
DO i=1,nReducedCleanPoints
 if (ParCleanList%Outflow(i))  write(iunit,'(I0)') i
END DO
CLOSE(iunit)

WRITE(cf,'(A)') ADJUSTL(TRIM(cOutputFolder))//'/'//ADJUSTL(TRIM(cReducedMeshdFolder))//'/Outflow.pls'
WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cf))//"'"
OPEN(UNIT=iunit,FILE=ADJUSTL(TRIM(cf)))
write(iunit,'(I0,A)') ParCleanList%nOA,' '//'Outflow'
write(iunit,'(A,I0," ",4F10.4,A)') "' '"
DO i=1,ParCleanList%nOA
 write(iunit,'(4(I0,A))') ParCleanList%tOA%i(1,i),',',ParCleanList%tOA%i(2,i),',',ParCleanList%tOA%i(3,i),',',ParCleanList%tOA%i(4,i)
END DO
CLOSE(iunit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
do iInflow=1,myProcess%nOfInflows
 WRITE(cf,'(A,I0,A)') ADJUSTL(TRIM(cOutputFolder))//'/'//ADJUSTL(TRIM(cReducedMeshdFolder))//'/Inflow',iInflow,'.par'
 WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cf))//"'"
 OPEN(UNIT=iunit,FILE=ADJUSTL(TRIM(cf)))
 nBC = 0
 DO i=1,nReducedCleanPoints
  if (ParCleanList%Inflow(iInflow,i)) nBC = nBC + 1
 END DO
 IF (nBC.eq.0) THEN
  write(*,*) 'Either the inflow radius too small or the meshresolution too coarse, or the inflow position is not matching with the bounding box!'
  write(*,*) 'Setup is not suitable for simulation //  program forced to exit!'
  STOP
 END IF
 write(iunit,'(I0,A,I0)') nBC,' '//'Inflow',-iInflow
 write(iunit,'(A,I0," ",4F10.4,A)') "' '"
 DO i=1,nReducedCleanPoints
  if (ParCleanList%Inflow(iInflow,i))  write(iunit,'(I0)') i
 END DO
 CLOSE(iunit)
END DO

do iInflow=1,myProcess%nOfInflows
 WRITE(cf,'(A,I0,A)') ADJUSTL(TRIM(cOutputFolder))//'/'//ADJUSTL(TRIM(cReducedMeshdFolder))//'/Inflow',iInflow,'.pls'
 WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cf))//"'"
 OPEN(UNIT=iunit,FILE=ADJUSTL(TRIM(cf)))
 IF (ParCleanList%nIA(iInflow).eq.0) THEN
  write(*,*) 'Either the inflow radius too small or the meshresolution too coarse, or the inflow position is not matching with the bounding box!'
  write(*,*) 'Setup is not suitable for simulation //  program forced to exit!'
  STOP
 END IF
 write(iunit,'(I0,A,I0)') ParCleanList%nIA(iInflow),' '//'Inflow',-iInflow
 write(iunit,'(A,I0," ",4F10.4,A)') "' '"
 DO i=1,ParCleanList%nIA(iInflow)
  write(iunit,'(4(I0,A))') ParCleanList%tIA(iInflow)%i(1,i),',',ParCleanList%tIA(iInflow)%i(2,i),',',ParCleanList%tIA(iInflow)%i(3,i),',',ParCleanList%tIA(iInflow)%i(4,i)
 END DO
 CLOSE(iunit)
END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11



WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cOutputFolder))//'/'//ADJUSTL(TRIM(cReducedMeshdFolder))//'/file.prj'//"'"
open(file=ADJUSTL(TRIM(cOutputFolder))//'/'//ADJUSTL(TRIM(cReducedMeshdFolder))//'/file.prj',unit=iunit)
write(iunit,'(a)') adjustl(trim(cReducedGridFile))
write(iunit,'(a)') 'Wall.par'
write(iunit,'(a)') 'Outflow.par'
do iInflow=1,myProcess%nOfInflows
 write(iunit,'(A,I0,A)') 'Inflow',iInflow,'.par'
END DO
do jat=1,5
 write(iunit,'(A,I0,A)') 'SideWall',jat,'.par'
END DO
close(iunit)

WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cOutputFolder))//'/'//ADJUSTL(TRIM(cReducedMeshdFolder))//'/file_pls.prj'//"'"
open(file=ADJUSTL(TRIM(cOutputFolder))//'/'//ADJUSTL(TRIM(cReducedMeshdFolder))//'/file_pls.prj',unit=iunit)
write(iunit,'(a)') adjustl(trim(cReducedGridFile))
write(iunit,'(a)') 'Wall.pls'
write(iunit,'(a)') 'Outflow.pls'
do iInflow=1,myProcess%nOfInflows
 write(iunit,'(A,I0,A)') 'Inflow',iInflow,'.pls'
END DO
do jat=1,5
 write(iunit,'(A,I0,A)') 'SideWall',jat,'.pls'
END DO
close(iunit)

END SUBROUTINE Output_ReducedRefTRIandPAR
! ----------------------------------------------

SUBROUTINE Output_ReducedRefTRI()
IMPLICIT NONE
INTEGER nel,nvt
INTEGER i,j, iel,jel,ive,ivt,ioffset,n
integer nnel,nnvt
INTEGER :: iunit=123
CHARACTER*(256) cf

WRITE(cf,'(A)') ADJUSTL(TRIM(cOutputFolder))//'/'//adjustl(trim(cReducedGridFile))
WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cf))//"'"
OPEN(UNIT=iunit,FILE=ADJUSTL(TRIM(cf)))
WRITE(iunit,*) 'Coarse mesh exported by DeViSoR TRI3D exporter'
WRITE(iunit,*) 'Parametrisierung PARXC, PARYC, TMAXC'
WRITE(iunit,'(2I8,A)') nReducedElems,nReducedPoints, " 1 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE"

WRITE(iunit,'(A)') 'DCORVG'
DO i = 1,nReducedPoints
 WRITE(iunit,'(3ES13.5)') ReducedMeshCoor(:,i)
END DO

WRITE(iunit,'(A)') 'KVERT'
DO i = 1,nReducedElems
 WRITE(iunit,'(8I8)') ReducedMeshElem(:,i)
END DO

WRITE(iunit,'(A)') 'KNPR'
DO i = 1,nReducedPoints
  WRITE(iunit,'(I8)') 0
END DO

CLOSE(iunit)

END SUBROUTINE Output_ReducedRefTRI
! ----------------------------------------------
SUBROUTINE Output_MergedRefVTK

IMPLICIT NONE
INTEGER nel,nvt
INTEGER i,j, iel,jel,ive,ivt,ioffset,n
integer nnel,nnvt
INTEGER :: iunit=123
CHARACTER*(100) filename

nel = mg_mesh%level(ilev)%nel
nvt = mg_mesh%level(ilev)%nvt

nnel = nUniqueElems
nnvt = nUniquePoints

filename=" "
WRITE(filename(1:),'(A)') adjustl(trim(cOutputFolder))//"/MergedRefMesh.vtu"

WRITE(*,'(104("="))') 
WRITE(*,*) "Outputting vtk file into ",filename

OPEN (UNIT=iunit,FILE=filename)

write(iunit, *)"<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">"
write(iunit, *)"  <UnstructuredGrid>"
write(iunit, *)"    <Piece NumberOfPoints=""",nnvt,""" NumberOfCells=""",nnel,""">"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the node field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <PointData>"

!  write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","MarkerV",""" format=""ascii"">"
!  do ivt=1,nvt
!   write(iunit, '(A,I10)')"        ",markerV(ivt)
!  end do
!  write(iunit, *)"        </DataArray>"
! 
!  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","ID",""" format=""ascii"">"
!  do ivt=1,nvt
!   write(iunit, '(A,E16.7)')"        ",REAL(ivt)
!  end do
!  write(iunit, *)"        </DataArray>"

 write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","KNPR",""" format=""ascii"">"
 do ivt=1,nnvt
  write(iunit, '(A,E16.7)')"        ",REAL(MergedMeshKnpr(ivt))
 end do
 write(iunit, *)"        </DataArray>"

 ! write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","knpr",""" format=""ascii"">"
! do iel=1,nel
!  if (.not.allocated(myRF(iel)%myRF)) then
!   do i=1,myRF(iel)%nOfVert
!     write(iunit,'(A10,I0)')"          ",myRF(iel)%knpr(i)
!   end do
!  else
!   do jel=1,myRF(iel)%nOfElem
!    do i=1,myRF(iel)%myRF(jel)%nOfVert
!      write(iunit,'(A10,I0)')"          ",myRF(iel)%myRF(jel)%knpr(i)
!    end do
!   end do
!  end if  
! end do
! write(iunit, *)"        </DataArray>"

write(iunit, '(A)')"    </PointData>"
write(iunit, '(A)')"      <Points>"
write(iunit, '(A)')"        <DataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3"" format=""ascii"" RangeMin=""0"" RangeMax=""1.0"">"
do i=1,nnvt
   write(iunit,'(A10,3E16.7)')"          ",REAL(MergedMeshCoor(1,i)),REAL(MergedMeshCoor(2,i)),REAL(MergedMeshCoor(3,i))
end do
write(iunit, *)"        </DataArray>"
write(iunit, *)"      </Points>"

write(iunit, '(A)')"    <CellData>"


! write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","NREFEL",""" format=""ascii"">"
! do iel=1,nel
!  n=myRF(iel)%nOfElem
!  if (.not.allocated(myRF(iel)%myRF)) then
!   do i=1,myRF(iel)%nOfElem
!    write(iunit, '(A,I10)')"        ",n !MarkerE(ivt)
!   end do
!  else
!   do jel=1,myRF(iel)%nOfElem
!    do i=1,myRF(iel)%myRF(jel)%nOfElem
!     write(iunit, '(A,I10)')"        ",n !MarkerE(ivt)
!    end do
!   end do
!  end if
! end do
! write(iunit, *)"        </DataArray>"

write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","PATCH_ID",""" format=""ascii"">"
do iel=1,nel
  do i=1,myRF(iel)%nUniqueElems
   write(iunit, '(A,I10)')"        ",myRF(iel)%patchID 
  end do
end do
write(iunit, *)"        </DataArray>"

write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","DIST",""" format=""ascii"">"
do iel=1,nUniqueElems
 write(iunit, '(A,E16.7)')"        ",REAL(MergedMeshDist(iel))
end do
write(iunit, *)"        </DataArray>"

write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","EL_ID",""" format=""ascii"">"
do iel=1,nUniqueElems
 write(iunit, '(A,E16.7)')"        ",REAL(iel)
end do
write(iunit, *)"        </DataArray>"

! write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","Selection",""" format=""ascii"">"
! do iel=1,nel
!  n=myRF(iel)%nOfElem
!  if (.not.allocated(myRF(iel)%myRF)) then
!   do i=1,myRF(iel)%nOfElem
!    write(iunit, '(A,I10)')"        ",myRF(iel)%selection
!   end do
!  else
!   do jel=1,myRF(iel)%nOfElem
!    do i=1,myRF(iel)%myRF(jel)%nOfElem
!     write(iunit, '(A,I10)')"        ",myRF(iel)%myRF(jel)%selection
!    end do
!   end do
!  end if
! end do
! write(iunit, *)"        </DataArray>"

write(iunit, '(A)')"    </CellData>"

write(iunit, *)"      <Cells>"
write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""connectivity"" format=""ascii"" RangeMin=""0"" RangeMax=""",nel-1,""">"
do iel=1,nnel
   write(iunit, '(8I10)') MergedMeshElem(1,iel) - 1,MergedMeshElem(2,iel) - 1,MergedMeshElem(3,iel) - 1,MergedMeshElem(4,iel) - 1,&
                          MergedMeshElem(5,iel) - 1,MergedMeshElem(6,iel) - 1,MergedMeshElem(7,iel) - 1,MergedMeshElem(8,iel) - 1
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""offsets"" format=""ascii"" RangeMin=""8"" RangeMax=""",8*nel,""">"
ioffset=nnel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')ive*8,(ive+1)*8,(ive+2)*8,(ive+3)*8,(ive+4)*8,(ive+5)*8,(ive+6)*8,(ive+7)*8
end do

do ive=ioffset+1,nnel
 write(iunit, '(I10)')ive*8
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A)')"        <DataArray type=""UInt8"" Name=""types"" format=""ascii"" RangeMin=""12"" RangeMax=""12"">"
ioffset=nnel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')12,12,12,12,12,12,12,12
end do
do ive=ioffset+1,nnel
 write(iunit, '(I10)')12
end do
write(iunit, '(A)')"        </DataArray>"
 
write(iunit, *)"      </Cells>"
write(iunit, *)"    </Piece>"
   
write(iunit, *)"  </UnstructuredGrid>"
write(iunit, *)"</VTKFile>"
close(iunit)

END SUBROUTINE Output_MergedRefVTK
! ----------------------------------------------
SUBROUTINE Output_RefTriMesh()
USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:mg_mesh,myBoundary
USE Parametrization, ONLY : myParBndr,nBnds
IMPLICIT NONE
INTEGER i,j,iBnds,nnvt,nnel
CHARACTER cf*(256)

nnvt = 0
nnel = 0

DO i=1,nel
 nnel = nnel + myRF(i)%nOfElem
 nnvt = nnvt + myRF(i)%nOfVert
end do

WRITE(cf,'(A)') ADJUSTL(TRIM(cOutputFolder))//'/Ref'//adjustl(trim(cProjectGridFile))
WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cf))//"'"
OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cf)))
WRITE(1,*) 'Coarse mesh exported by DeViSoR TRI3D exporter'
WRITE(1,*) 'Parametrisierung PARXC, PARYC, TMAXC'
WRITE(1,'(2I8,A)') nNEL,nNVT, " 1 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE"

WRITE(1,'(A)') 'DCORVG'
do iel=1,nel
 do i=1,myRF(iel)%nOfVert
  write(1,'(A10,3E16.7)')"          ",REAL(myRF(iel)%dcoor(1,i)),REAL(myRF(iel)%dcoor(2,i)),REAL(myRF(iel)%dcoor(3,i))
 end do
end do

WRITE(1,'(A)') 'KVERT'
nnvt = 0
do iel=1,nel
 do i=1,myRF(iel)%nOfElem
  write(1, '(8I10)') nnvt + myRF(iel)%kvert(1,i), nnvt + myRF(iel)%kvert(2,i), nnvt + myRF(iel)%kvert(3,i), nnvt + myRF(iel)%kvert(4,i),&
                         nnvt + myRF(iel)%kvert(5,i), nnvt + myRF(iel)%kvert(6,i), nnvt + myRF(iel)%kvert(7,i), nnvt + myRF(iel)%kvert(8,i) 
 end do
 nnvt = nnvt + myRF(iel)%nOfVert
end do

WRITE(1,'(A)') 'KNPR'
do iel=1,nel
 do i=1,myRF(iel)%nOfVert
  WRITE(1,'(I8)') 0
 end do
end do

CLOSE(1)

END SUBROUTINE Output_RefTriMesh
! ----------------------------------------------
SUBROUTINE Output_TriMesh()
USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:mg_mesh,myBoundary
USE Parametrization, ONLY : myParBndr,nBnds
IMPLICIT NONE
INTEGER i,j,iBnds
CHARACTER cf*(256)

WRITE(cf,'(A)') ADJUSTL(TRIM(cOutputFolder))//'/'//adjustl(trim(cProjectGridFile))
WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cf))//"'"
OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cf)))
WRITE(1,*) 'Coarse mesh exported by DeViSoR TRI3D exporter'
WRITE(1,*) 'Parametrisierung PARXC, PARYC, TMAXC'
WRITE(1,'(2I8,A)') mg_mesh%level(ilev)%NEL,mg_mesh%level(ilev)%NVT, " 1 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE"

WRITE(1,'(A)') 'DCORVG'
DO i = 1,mg_mesh%level(ilev)%nvt
 WRITE(1,'(3ES13.5)') mg_mesh%level(ilev)%dcorvg(:,i)
END DO

WRITE(1,'(A)') 'KVERT'
DO i = 1,mg_mesh%level(ilev)%nel
 WRITE(1,'(8I8)') mg_mesh%level(ILEV)%kvert(:,i)
END DO

WRITE(1,'(A)') 'KNPR'
IF (allocated(myBoundary%bWall)) then
 DO i = 1,mg_mesh%level(ilev)%nvt
  IF (myBoundary%bWall(i)) THEN
   WRITE(1,'(I8)') 1
  ELSE
   WRITE(1,'(I8)') 0
  END IF
 END DO
ELSE
 DO i = 1,mg_mesh%level(ilev)%nvt
  WRITE(1,'(I8)') 0
 END DO
END IF

CLOSE(1)

! OPEN(UNIT=2,FILE=ADJUSTL(TRIM(cOutputFolder))//'/'//adjustl(trim(cShortProjectFile)))
! !WRITE(cf,'(A11,I3.3,A4)') 'cMESH_',iO, '.tri'
! WRITE(2,'(A)') adjustl(trim(cProjectGridFile))
!  
! DO iBnds = 1, nBnds
!  cf = ' '
!  WRITE(cf,'(A)') ADJUSTL(TRIM(cOutputFolder))//"/"//ADJUSTL(TRIM(myParBndr(iBnds)%Names))//".par"
!  WRITE(2,'(A)') ADJUSTL(TRIM(myParBndr(iBnds)%Names))//".par"
!  WRITE(*,*) "Outputting actual parametrization into: '"//ADJUSTL(TRIM(cf))//"'"
!  OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cf)))
!  j=0
!  DO i=1,mg_mesh%level(ilev)%nvt
!   IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
!    j = j + 1
!   END IF
!  END DO
!  WRITE(1,'(I8,A)') j," "//myParBndr(iBnds)%Types
!  WRITE(1,'(A)')    "'"//ADJUSTL(TRIM(myParBndr(iBnds)%Parameters))//"'"
!  DO i=1,mg_mesh%level(ilev)%nvt
!   IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
!    WRITE(1,'(I8,A)') i
!   END IF
!  END DO
!  CLOSE(1)
! END DO
! CLOSE(2)

END SUBROUTINE Output_TriMesh
! ----------------------------------------------
SUBROUTINE Output_MergedRefTriMeshPar()

USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:mg_mesh,myBoundary
USE Parametrization, ONLY : myParBndr,nBnds
use Sigma_User, only : myProcess,mySigma

IMPLICIT NONE
INTEGER i,j,iloc,i1,i2,i3,i4,i5,i6
integer iat, iel, ivt(4), ilev,iX,jX
CHARACTER cf*(256),ctxt*(256),cInputFile*(256),cVal*(256),cKey*(256)
real*8 P0(3),box(3),bound(3,2)
real*8 :: dEps=1d-4
logical :: bExist
type(tMultiMesh),save :: mg_NewMesh
INTEGER NeighA(4,6),nbox(3),iInflow,jInflow
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
REAL*8 P(3),Q(3),dist,mindist,dSize,minDistP!,dVolume
logical bToMarkFace
logical, allocatable :: bInflowMarker(:,:)
logical, allocatable, dimension(:,:) :: HexSide
integer, allocatable :: InflowToSideMapper(:)
real*8 OneSideCloseD,SideClose
integer OneSideCloseI,iSide
real*8 dn(3)
integer iInflowSide(2)

real*8 :: DA(3), DB(3), dCenter(3), dAux1, dAux2 
real*8 :: dAdC, dBdC, dPdA, dPdB, dR

 !------------------------------------------------------------------
 cInputFile = ADJUSTL(TRIM(cIntputFolder))//'/'//'setup.e3d'
 inquire(file=cInputFile,Exist=bExist)
 if (bExist) then
  call ReadS3Dfile(cInputFile)
 else
  write(*,*) 'file: "',adjustl(trim(cInputFile)),'" does not exist!'
  STOP
 end if
 !------------------------------------------------------------------

P0  = 10d0*mySigma%DIE_Start  
box = 10d0*mySigma%DIE_Length

! nbox = mySigma%DIE_Voxels
!  cInputFile = ADJUSTL(TRIM(cIntputFolder))//'/'//'param.txt'
! 
!  cKey='geometryStart'
!  CALL GetValueFromFile(cInputFile,cVal,cKey)
!  read(cVal,*) P0
! 
!  cKey='geometryLength'
!  CALL GetValueFromFile(cInputFile,cVal,cKey)
!  read(cVal,*) box
! 
!  cKey='voxelAmount'
!  CALL GetValueFromFile(cInputFile,cVal,cKey)
!  read(cVal,*) nBox
 

!  dVolume = box(1)*box(2)*box(3)/(nbox(1)*nbox(2)*nbox(3)) ! volume of one voxel
!  dSize   = dVolume**(1d0/3d0)
!  WRITE(*,*) "dSize  =",dSize
!  
 dSize = (box(1)*box(2)*box(3))**(1d0/3d0)
 dEps  = dSize/1e5

 bound(:,1) = P0
 bound(:,2) = P0 + box
 OverallBoundingBox = bound
 write(*,*) bound(:,1)
 write(*,*) bound(:,2)

 i1=0
 i2=0
 i3=0
 i4=0
 i5=0
 i6=0
 open(file=ADJUSTL(TRIM(cOutputFolder))//'/meshDir_BU/x-.par',unit=11)
 open(file=ADJUSTL(TRIM(cOutputFolder))//'/meshDir_BU/x+.par',unit=12)
 open(file=ADJUSTL(TRIM(cOutputFolder))//'/meshDir_BU/y-.par',unit=13)
 open(file=ADJUSTL(TRIM(cOutputFolder))//'/meshDir_BU/y+.par',unit=14)
 open(file=ADJUSTL(TRIM(cOutputFolder))//'/meshDir_BU/z-.par',unit=15)
 open(file=ADJUSTL(TRIM(cOutputFolder))//'/meshDir_BU/z+.par',unit=16)
 allocate(HexSide(6,nUniquePoints))
 HexSide = .false.

 do i=1,nUniquePoints
  P0 = MergedMeshCoor(:,i)
  if (abs(P0(1)-bound(1,1)).lt.dEps) then
   i1 = i1 + 1  
   HexSide(1,i) = .true.
  end if
  if (abs(P0(1)-bound(1,2)).lt.dEps) then
   i2 = i2 + 1  
   HexSide(2,i) = .true.
  end if
  if (abs(P0(2)-bound(2,1)).lt.dEps) then
   i3 = i3 + 1  
   HexSide(3,i) = .true.
  end if
  if (abs(P0(2)-bound(2,2)).lt.dEps) then
   i4 = i4 + 1  
   HexSide(4,i) = .true.
  end if
  if (abs(P0(3)-bound(3,1)).lt.dEps) then
   i5 = i5 + 1  
   HexSide(5,i) = .true.
  end if
  if (abs(P0(3)-bound(3,2)).lt.dEps) then
   i6 = i6 + 1  
   HexSide(6,i) = .true.
  end if
 end do

 CALL writeBC(11,i1,1d0,0d0,0d0,-MeshOutputScaleFactor*bound(1,1),mySigma%DIE_SymmetryBC(1),'Wall','100')
 CALL writeBC(12,i2,1d0,0d0,0d0,-MeshOutputScaleFactor*bound(1,2),mySigma%DIE_SymmetryBC(4),'Wall','100')
 CALL writeBC(13,i3,0d0,1d0,0d0,-MeshOutputScaleFactor*bound(2,1),mySigma%DIE_SymmetryBC(2),'Wall','010')
 CALL writeBC(14,i4,0d0,1d0,0d0,-MeshOutputScaleFactor*bound(2,2),mySigma%DIE_SymmetryBC(5),'Wall','010')
 CALL writeBC(15,i5,0d0,0d0,1d0,-MeshOutputScaleFactor*bound(3,1),.FALSE.,'Wall','001')
 CALL writeBC(16,i6,0d0,0d0,1d0,-MeshOutputScaleFactor*bound(3,2),.FALSE.,'Outflow','001')
 
!  write(11,'(I0,A)') i1,' '//'Wall'
!  write(11,'(A,I0," ",4F10.4,A)') "'",4,1.0, 0.0, 0.0, -MeshOutputScaleFactor*bound(1,1), "'"
!  write(12,'(I0,A)') i2,' '//'Wall'
!  write(12,'(A,I0," ",4F10.4,A)') "'",4,1.0, 0.0, 0.0, -MeshOutputScaleFactor*bound(1,2), "'"
!  write(13,'(I0,A)') i3,' '//'Wall'
!  write(13,'(A,I0," ",4F10.4,A)') "'",4,0.0, 1.0, 0.0, -MeshOutputScaleFactor*bound(2,1), "'"
!  write(14,'(I0,A)') i4,' '//'Wall'
!  write(14,'(A,I0," ",4F10.4,A)') "'",4,0.0, 1.0, 0.0, -MeshOutputScaleFactor*bound(2,2), "'"
!  write(15,'(I0,A)') i5,' '//'Wall'
!  write(15,'(A,I0," ",4F10.4,A)') "'",4,0.0, 0.0, 1.0, -MeshOutputScaleFactor*bound(3,1), "'"
!  write(16,'(I0,A)') i6,' '//'Outflow'
!  write(16,'(A,I0," ",4F10.4,A)') "'",4,0.0, 0.0, 1.0, -MeshOutputScaleFactor*bound(3,2), "'"

 write(*,*) i1,i2,i3,i4,i5,i6

 do i=1,nUniquePoints
  P0 = MergedMeshCoor(:,i)
  if (abs(P0(1)-bound(1,1)).lt.dEps) then
   write(11,*) i
  end if
  if (abs(P0(1)-bound(1,2)).lt.dEps) then
   write(12,*) i
  end if
  if (abs(P0(2)-bound(2,1)).lt.dEps) then
   write(13,*) i
  end if
  if (abs(P0(2)-bound(2,2)).lt.dEps) then
   write(14,*) i
  end if
  if (abs(P0(3)-bound(3,1)).lt.dEps) then
   write(15,*) i
  end if
  if (abs(P0(3)-bound(3,2)).lt.dEps) then
   write(16,*) i
  end if
 end do

 close(11)
 close(12)
 close(13)
 close(14)
 close(15)
 close(16)

 allocate(bInflowMarker(myProcess%nOfInflows,nUniquePoints))
 
 bInflowMarker = .false.
 
 !  deallocate(mg_mesh%level)
!  !!! building up the new mesh structures
 mg_NewMesh%maxlevel = mg_Mesh%maxlevel
 allocate(mg_NewMesh%level(mg_NewMesh%maxlevel))
 allocate(mg_NewMesh%level(1)%dcorvg(3,nUniquePoints))
 allocate(mg_NewMesh%level(1)%knpr(nUniquePoints))
 allocate(mg_NewMesh%level(1)%kvert(8,nUniqueElems))
 mg_NewMesh%level(1)%nvt = nUniquePoints
 mg_NewMesh%level(1)%nel = nUniqueElems
 mg_NewMesh%level(1)%knpr = 0
 mg_NewMesh%level(1)%dcorvg = MergedMeshCoor
 mg_NewMesh%level(1)%kvert  = MergedMeshElem
 call refineMesh(mg_NewMesh, mg_NewMesh%maxlevel)  
 
 DO ilev = 1,mg_NewMesh%maxlevel
  write(*,*) mg_NewMesh%level(ilev)%nvt,mg_NewMesh%level(ilev)%nel,mg_NewMesh%level(ilev)%net,mg_NewMesh%level(ilev)%nat
 END DO

 allocate(InflowToSideMapper(myProcess%nOfInflows))
 do iInflow=1,myProcess%nOfInflows
  Q = myProcess%myInflow(iInflow)%Center
  dn = myProcess%myInflow(iInflow)%normal
  
  iInflowSide(:) = 0
  if (abs(dn(1)).gt.0.9d0) iInflowSide(:)= [1,2]
  if (abs(dn(2)).gt.0.9d0) iInflowSide(:)= [3,4]
  if (abs(dn(3)).gt.0.9d0) iInflowSide(:)= [5,6]
  if (iInflowSide(1).eq.0.or.iInflowSide(1).eq.0) THEN
   WRITE(*,*) 'Inflow sides were not identified. program stops'
   STOP 60
  END IF
  
  OneSideCloseD = 1d8
  OneSideCloseI = 0
  
  do iSide=1,6
   SideClose =1d8
   do i=1,nUniquePoints
    if (HexSide(iSide,i).and.(iSide.eq.iInflowSide(1).or.iSide.eq.iInflowSide(2))) then
     P = MeshOutputScaleFactor*mg_NewMesh%level(1)%dcorvg(:,i)
     dist = sqrt((P(1)-Q(1))**2d0 + (P(2)-Q(2))**2d0 + (P(3)-Q(3))**2d0) 
     if (dist.lt.SideClose) SideClose = dist
    end if
   end do
   if (OneSideCloseD.gt.SideClose) then
    OneSideCloseD = SideClose
    OneSideCloseI = iSide
   end if
  end do
  InflowToSideMapper(iInflow) = OneSideCloseI
 end do

 write(*,*) 'InflowToSideMapper = ',InflowToSideMapper

 DO iel = 1, mg_NewMesh%level(1)%nel
 
  do iat = 1,6
   
   if (mg_NewMesh%level(1)%kadj(iat,iel).eq.0) then
    ivt(1) = mg_NewMesh%level(1)%kvert(NeighA(1,iat),iel)
    ivt(2) = mg_NewMesh%level(1)%kvert(NeighA(2,iat),iel)
    ivt(3) = mg_NewMesh%level(1)%kvert(NeighA(3,iat),iel)
    ivt(4) = mg_NewMesh%level(1)%kvert(NeighA(4,iat),iel)
    
    P = 0d0
    DO i=1,4
     P = P + 0.25d0*mg_NewMesh%level(1)%dcorvg(:,ivt(i))
    end do
    
    minDistP = 1d8
    DO i=1,4
     Q = mg_NewMesh%level(1)%dcorvg(:,ivt(i))
     dist = sqrt((P(1)-Q(1))**2d0 + (P(2)-Q(2))**2d0 + (P(3)-Q(3))**2d0) 
     IF (dist.lt.minDistP) minDistP = dist
    end do
    
    iSide = 0
    do i=1,6
     if (HexSide(i,ivt(1)).and.HexSide(i,ivt(2)).and.HexSide(i,ivt(3)).and.HexSide(i,ivt(4))) then
      iSide = i
     end if
    end do
    
!     ! lets check if the face is really a boundary face (only for thew brick!)
!     iX =0
!     jX =0 
!     DO i=1,3
!      DO j=1,2
!       if (abs(bound(i,j) - P(i)).lt.0.1d0*dSize) then
!        iX = i
!        jX = j
!       end if
!      end do
!     END DO
!     
!     IF (iX.eq.0.or.jX.eq.0) THEN
!      WRITE(*,*) 'not identified face has been detected... '
!      stop
!     END IF

      mindist = 1d8
      jInflow = 0
      
      P = MeshOutputScaleFactor*P ! scaling to cm
      DO iInflow=1,myProcess%nOfInflows
       if (iSide.eq.InflowToSideMapper(iInflow)) then
        Q = myProcess%myInflow(iInflow)%Center
        dist = sqrt((P(1)-Q(1))**2d0 + (P(2)-Q(2))**2d0 + (P(3)-Q(3))**2d0)
        if (dist.lt.mindist) then
         jInflow = iInflow 
         mindist = dist
        end if
       end if
      END DO
      
      bToMarkFace = .false.
      if (jInflow.ge.1.and.jInflow.le.myProcess%nOfInflows) then
       if (iSide.eq.InflowToSideMapper(jInflow)) then
        do i=1,4
         if (myProcess%myInflow(jInflow)%iBCtype.LT.5) then
           P = mg_NewMesh%level(1)%dcorvg(:,ivt(i))
           Q = myProcess%myInflow(jInflow)%Center/MeshOutputScaleFactor
           dist = sqrt((P(1)-Q(1))**2d0 + (P(2)-Q(2))**2d0 + (P(3)-Q(3))**2d0)
           if (dist.lt.myProcess%myInflow(jInflow)%outerradius/MeshOutputScaleFactor+minDistP) then
  !          if (dist.lt.myProcess%myInflow(jInflow)%outerradius+minDistP) then
            bToMarkFace = .true.
           end if
         elseif (myProcess%myInflow(jInflow)%iBCType.EQ.5) then
           P = mg_NewMesh%level(1)%dcorvg(:,ivt(i))
           
           dCenter = myProcess%myInflow(jInflow)%center/MeshOutputScaleFactor
           DA = myProcess%myInflow(jInflow)%midpointA/MeshOutputScaleFactor
           DB = myProcess%myInflow(jInflow)%midpointB/MeshOutputScaleFactor

           dAux1 = DOT_PRODUCT(DA-dCenter, DA-dCenter)**2&
                  -DOT_PRODUCT(P-dCenter, DA-dCenter)**2+0.5*minDistP
           dAux2 = DOT_PRODUCT(DB-dCenter, DB-dCenter)**2&
                  -DOT_PRODUCT(P-dCenter, DB-dCenter)**2+0.5*minDistP
           IF ( (dAux1.GE.0D0).and.(dAux2.GE.0D0) ) THEN
             bToMarkFace = .true.
           END IF
         elseif (myProcess%myInflow(jInflow)%iBCType.EQ.6) then
           P = mg_NewMesh%level(1)%dcorvg(:,ivt(i))
           
           dCenter = myProcess%myInflow(jInflow)%center/MeshOutputScaleFactor
           DA = myProcess%myInflow(jInflow)%midpointA/MeshOutputScaleFactor
           DB = myProcess%myInflow(jInflow)%midpointB/MeshOutputScaleFactor

           dAdC = DOT_PRODUCT(dA-dCenter, DA-dCenter)
           dBdC = DOT_PRODUCT(dB-dCenter, DB-dCenter)

           dPdA = DOT_PRODUCT(P-dCenter, DA-dCenter)
           dPdB = DOT_PRODUCT(P-dCenter, DB-dCenter)

           dR = NORM2(DB-dCenter)

           if ( (dPdA**2.LE.dAdC**2+0.5*minDistP).and.(dPdB**2.LE.dBdC**2+0.5*minDistP) ) THEN
               bToMarkFace = .true.
           elseif ( (dPdA.GE.dAdC).and.(NORM2(P-DA).LE.dR+0.5*minDistP) ) THEN
               bToMarkFace = .true.
           elseif ( (dPdA**2.GE.dAdC**2).and.(NORM2(P-2*dCenter+DA).LE.dR+0.5*minDistP) ) THEN
               bToMarkFace = .true.
           end if
         end if
        end do
       end if
       
       if (bToMarkFace) then
 !       WRITE(*,*) iel,iat
        bInflowMarker(jInflow,ivt(:)) = .TRUE.
       end if
      end if

!     WRITE(*,'(I0," ",I0,3ES12.4,":",6ES12.4)') iX,jX,P,bound
   end if
   
  end do
 
 END DO

 DO iInflow=1,myProcess%nOfInflows
  jInflow = 0
  DO i=1,nUniquePoints
   if (bInflowMarker(iInflow,i)) jInflow = jInflow + 1
  END DO
  WRITE(cInputFile,'(A,I0,A)') ADJUSTL(TRIM(cOutputFolder))//'/meshDir_BU/inflow_',iInflow,'.par'
  open(file=ADJUSTL(TRIM(cInputFile)),unit=5)
  write(5,'(I0,A,I0)') jInflow, ' Inflow-',iInflow
  write(5,'(A)') "' '"
  DO i=1,nUniquePoints
   if (bInflowMarker(iInflow,i)) THEN
    write(5,'(I0)') i
   END IF
  END DO
  close(5)
 END DO

 open(file=ADJUSTL(TRIM(cOutputFolder))//'/meshDir_BU/file.prj',unit=5)
 write(5,'(a)') 'Merged_'//adjustl(trim(cProjectGridFile))

 write(5,'(a)') 'x+.par'
 write(5,'(a)') 'x-.par'
 write(5,'(a)') 'y+.par'
 write(5,'(a)') 'y-.par'
 write(5,'(a)') 'z+.par'
 write(5,'(a)') 'z-.par'
 DO iInflow=1,myProcess%nOfInflows
  write(5,'(A,I0,A)') 'inflow_',iInflow,'.par'
 END DO

 close(5)
 
 contains
 
 subroutine writeBC(iunit,iN,p1,p2,p3,p4,bBC,cBC1,cBC2)
 integer iunit,iN
 real*8 p1,p2,p3,p4
 logical bBC
 character cBC1*(*),cBC2*(*)
 
 if (bBC) then
  write(iunit,'(I0,A)') iN,' '//'Symmetry'//ADJUsTL(TRIM(cBC2))
 else
  write(iunit,'(I0,A)') iN,' '//ADJUsTL(TRIM(cBC1))
 end if
 write(iunit,'(A,I0," ",4F10.4,A)') "'",4,p1, p2, p3, p4, "'"

 end subroutine writeBC
 
END SUBROUTINE Output_MergedRefTriMeshPar
! ----------------------------------------------
SUBROUTINE Output_MergedRefTriMeshParOLD()

USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:mg_mesh,myBoundary
USE Parametrization, ONLY : myParBndr,nBnds
use Sigma_User, only : myProcess

IMPLICIT NONE
INTEGER i,j,iloc,i1,i2,i3,i4,i5,i6
CHARACTER cf*(256),ctxt*(256),cInputFile*(256),cVal*(256),cKey*(256)
CHARACTER cBC(6)*(256)
real*8 P0(3),box(3),bound(3,2)
real*8 :: dEps=1d-4
logical :: bExist

 !------------------------------------------------------------------
 cInputFile = ADJUSTL(TRIM(cIntputFolder))//'/'//'setup.e3d'
 inquire(file=cInputFile,Exist=bExist)
 if (bExist) then
  call ReadS3Dfile(cInputFile)
 else
  write(*,*) 'file: "',adjustl(trim(cInputFile)),'" does not exist!'
  STOP
 end if
 !------------------------------------------------------------------

 cInputFile = ADJUSTL(TRIM(cIntputFolder))//'/'//'param.txt'

 cKey='geometryStart'
 CALL GetValueFromFile(cInputFile,cVal,cKey)
 read(cVal,*) P0

 cKey='geometryLength'
 CALL GetValueFromFile(cInputFile,cVal,cKey)
 read(cVal,*) box

 cKey='XMinBC'
 CALL GetValueFromFile(cInputFile,cVal,cKey)
 read(cVal,*) cBC(1)

 cKey='XMaxBC'
 CALL GetValueFromFile(cInputFile,cVal,cKey)
 read(cVal,*) cBC(2)

 cKey='YMinBC'
 CALL GetValueFromFile(cInputFile,cVal,cKey)
 read(cVal,*) cBC(3)

 cKey='YMaxBC'
 CALL GetValueFromFile(cInputFile,cVal,cKey)
 read(cVal,*) cBC(4)

 cKey='ZMinBC'
 CALL GetValueFromFile(cInputFile,cVal,cKey)
 read(cVal,*) cBC(5)

 cKey='ZMaxBC'
 CALL GetValueFromFile(cInputFile,cVal,cKey)
 read(cVal,*) cBC(6)

bound(:,1) = P0
bound(:,2) = P0 + box
write(*,*) bound(:,1)
write(*,*) bound(:,2)

i1=0
i2=0
i3=0
i4=0
i5=0
i6=0
open(file=ADJUSTL(TRIM(cOutputFolder))//'/meshDir/x-.par',unit=11)
open(file=ADJUSTL(TRIM(cOutputFolder))//'/meshDir/x+.par',unit=12)
open(file=ADJUSTL(TRIM(cOutputFolder))//'/meshDir/y-.par',unit=13)
open(file=ADJUSTL(TRIM(cOutputFolder))//'/meshDir/y+.par',unit=14)
open(file=ADJUSTL(TRIM(cOutputFolder))//'/meshDir/z-.par',unit=15)
open(file=ADJUSTL(TRIM(cOutputFolder))//'/meshDir/z+.par',unit=16)

do i=1,nUniquePoints
 P0 = MergedMeshCoor(:,i)
 if (abs(P0(1)-bound(1,1)).lt.dEps) then
  i1 = i1 + 1  
 end if
 if (abs(P0(1)-bound(1,2)).lt.dEps) then
  i2 = i2 + 1  
 end if
 if (abs(P0(2)-bound(2,1)).lt.dEps) then
  i3 = i3 + 1  
 end if
 if (abs(P0(2)-bound(2,2)).lt.dEps) then
  i4 = i4 + 1  
 end if
 if (abs(P0(3)-bound(3,1)).lt.dEps) then
  i5 = i5 + 1  
 end if
 if (abs(P0(3)-bound(3,2)).lt.dEps) then
  i6 = i6 + 1  
 end if
end do

write(11,'(I0,A)') i1,' '//ADJUSTL(TRIM(cBC(1)))
write(11,'(A,I0," ",4F10.4,A)') "'",4,1.0, 0.0, 0.0, -MeshOutputScaleFactor*bound(1,1), "'"
write(12,'(I0,A)') i2,' '//ADJUSTL(TRIM(cBC(2)))
write(12,'(A,I0," ",4F10.4,A)') "'",4,1.0, 0.0, 0.0, -MeshOutputScaleFactor*bound(1,2), "'"
write(13,'(I0,A)') i3,' '//ADJUSTL(TRIM(cBC(3)))
write(13,'(A,I0," ",4F10.4,A)') "'",4,0.0, 1.0, 0.0, -MeshOutputScaleFactor*bound(2,1), "'"
write(14,'(I0,A)') i4,' '//ADJUSTL(TRIM(cBC(4)))
write(14,'(A,I0," ",4F10.4,A)') "'",4,0.0, 1.0, 0.0, -MeshOutputScaleFactor*bound(2,2), "'"
write(15,'(I0,A)') i5,' '//ADJUSTL(TRIM(cBC(5)))
write(15,'(A,I0," ",4F10.4,A)') "'",4,0.0, 0.0, 1.0, -MeshOutputScaleFactor*bound(3,1), "'"
write(16,'(I0,A)') i6,' '//ADJUSTL(TRIM(cBC(6)))
write(16,'(A,I0," ",4F10.4,A)') "'",4,0.0, 0.0, 1.0, -MeshOutputScaleFactor*bound(3,2), "'"

write(*,*) i1,i2,i3,i4,i5,i6

do i=1,nUniquePoints
 P0 = MergedMeshCoor(:,i)
 if (abs(P0(1)-bound(1,1)).lt.dEps) then
  write(11,*) i
 end if
 if (abs(P0(1)-bound(1,2)).lt.dEps) then
  write(12,*) i
 end if
 if (abs(P0(2)-bound(2,1)).lt.dEps) then
  write(13,*) i
 end if
 if (abs(P0(2)-bound(2,2)).lt.dEps) then
  write(14,*) i
 end if
 if (abs(P0(3)-bound(3,1)).lt.dEps) then
  write(15,*) i
 end if
 if (abs(P0(3)-bound(3,2)).lt.dEps) then
  write(16,*) i
 end if
end do

close(11)
close(12)
close(13)
close(14)
close(15)
close(16)

open(file=ADJUSTL(TRIM(cOutputFolder))//'/meshDir/file.prj',unit=5)
write(5,'(a)') 'Merged_'//adjustl(trim(cProjectGridFile))

write(5,'(a)') 'x+.par'
write(5,'(a)') 'x-.par'
write(5,'(a)') 'y+.par'
write(5,'(a)') 'y-.par'
write(5,'(a)') 'z+.par'
write(5,'(a)') 'z-.par'

close(5)


END SUBROUTINE Output_MergedRefTriMeshParOLD
! ----------------------------------------------
SUBROUTINE Output_MergedRefTriMesh()

USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:mg_mesh,myBoundary
USE Parametrization, ONLY : myParBndr,nBnds
IMPLICIT NONE
INTEGER i,j
CHARACTER cf*(256)

WRITE(cf,'(A)') ADJUSTL(TRIM(cOutputFolder))//'/meshDir_BU/Merged_'//adjustl(trim(cProjectGridFile))
WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cf))//"'"
OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cf)))
WRITE(1,*) 'Coarse mesh exported by DeViSoR TRI3D exporter'
WRITE(1,*) 'Parametrisierung PARXC, PARYC, TMAXC'
WRITE(1,'(2I8,A)') nUniqueElems,nUniquePoints, " 1 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE"

WRITE(1,'(A)') 'DCORVG'
DO i = 1,nUniquePoints
 WRITE(1,'(3ES13.5)') MeshOutputScaleFactor*MergedMeshCoor(:,i)
END DO

WRITE(1,'(A)') 'KVERT'
DO i = 1,nUniqueElems
 WRITE(1,'(8I8)') MergedMeshElem(:,i)
END DO

WRITE(1,'(A)') 'KNPR'
DO i = 1,nUniquePoints
  WRITE(1,'(I8)') 0
END DO

CLOSE(1)

END SUBROUTINE Output_MergedRefTriMesh
! ----------------------------------------------
SUBROUTINE Output_ParticleMergedRefTriMesh()

USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:mg_mesh,myBoundary
USE Parametrization, ONLY : myParBndr,nBnds
IMPLICIT NONE
INTEGER i,j
CHARACTER cf*(256)

WRITE(cf,'(A)') ADJUSTL(TRIM(cOutputFolder))//'/meshDir/Merged_'//adjustl(trim(cProjectGridFile))
WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cf))//"'"
OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cf)))
WRITE(1,*) 'Coarse mesh exported by DeViSoR TRI3D exporter'
WRITE(1,*) 'Parametrisierung PARXC, PARYC, TMAXC'
WRITE(1,'(2I8,A)') nUniqueElems,nUniquePoints, " 1 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE"

WRITE(1,'(A)') 'DCORVG'
DO i = 1,nUniquePoints
 WRITE(1,'(3ES13.5)') MergedMeshCoor(:,i)
END DO

WRITE(1,'(A)') 'KVERT'
DO i = 1,nUniqueElems
 WRITE(1,'(8I8)') MergedMeshElem(:,i)
END DO

WRITE(1,'(A)') 'KNPR'
DO i = 1,nUniquePoints
  WRITE(1,'(I8)') MergedMeshKnpr(i)
END DO

CLOSE(1)

! OPEN(UNIT=2,FILE=ADJUSTL(TRIM(cOutputFolder))//'/'//adjustl(trim(cShortProjectFile)))
! !WRITE(cf,'(A11,I3.3,A4)') 'cMESH_',iO, '.tri'
! WRITE(2,'(A)') adjustl(trim(cProjectGridFile))
!  
! DO iBnds = 1, nBnds
!  cf = ' '
!  WRITE(cf,'(A)') ADJUSTL(TRIM(cOutputFolder))//"/"//ADJUSTL(TRIM(myParBndr(iBnds)%Names))//".par"
!  WRITE(2,'(A)') ADJUSTL(TRIM(myParBndr(iBnds)%Names))//".par"
!  WRITE(*,*) "Outputting actual parametrization into: '"//ADJUSTL(TRIM(cf))//"'"
!  OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cf)))
!  j=0
!  DO i=1,mg_mesh%level(ilev)%nvt
!   IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
!    j = j + 1
!   END IF
!  END DO
!  WRITE(1,'(I8,A)') j," "//myParBndr(iBnds)%Types
!  WRITE(1,'(A)')    "'"//ADJUSTL(TRIM(myParBndr(iBnds)%Parameters))//"'"
!  DO i=1,mg_mesh%level(ilev)%nvt
!   IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
!    WRITE(1,'(I8,A)') i
!   END IF
!  END DO
!  CLOSE(1)
! END DO
! CLOSE(2)

END SUBROUTINE Output_ParticleMergedRefTriMesh

SUBROUTINE Output_ParticleMergedRefTriMeshPar()
USE MESH_Structures
use Sigma_User, only : mySigma
integer i,j,nKNPR,nParam,ParIndex
character cF*(256)
real*8 dc(3),P(3),minCoor(3),maxCoor(3),myRand(3)
integer iminCoor(3),imaxCoor(3)
REAL*8 :: dr=0.22d0
type(tMultiMesh) :: mg_NewMesh
CHARACTER cG*(256)
logical bFound


NLMAX = 2
NLMIN = 1
mg_NewMesh%maxlevel = 2
allocate(mg_NewMesh%level(mg_NewMesh%maxlevel))

WRITE(cG,'(A)') ADJUSTL(TRIM(cOutputFolder))//'/meshDir/Merged_'//adjustl(trim(cProjectGridFile))
call readTriCoarse(adjustl(trim(cG)), mg_NewMesh)
write(*,*) 'LoadedFile: ', adjustl(trim(cG))

!!! building up the mesh structures
call refineMesh(mg_NewMesh, mg_NewMesh%maxlevel)  

! do iel=1,mg_NewMesh%level(1)%nel
!  write(*,*) iel,':',mg_NewMesh%level(1)%kadj(:,iel)
!  pause
! end do
! pause

ParIndex = 0
do iel=1,mg_NewMesh%level(1)%nel
 bFound=.false.
 do j=1,8
  ivt = mg_NewMesh%level(1)%kvert(j,iel)
  if (mg_NewMesh%level(1)%knpr(ivt).gt.0) then
   bFound=.true.
  end IF
 end do
 
 if (bFound) THEN
  ParIndex = ParIndex + 1
  nParam = 0
  do j=1,8
   ivt = mg_NewMesh%level(1)%kvert(j,iel)
   if (mg_NewMesh%level(1)%knpr(ivt).gt.0) then
    nParam = nParam + 1
!    write(*,*) 'ivt ',ivt
    mg_NewMesh%level(1)%knpr(ivt) = -ParIndex
   end if
  end do
  CALL SetRecKNPR(iel)
!   write(*,*) ParIndex,nParam
 end if
! 
end do

write(cF,'(A)') ADJUSTL(TRIM(cOutputFolder))//'/meshDir/file.prj'
open(file=ADJUSTL(TRIM(cF)),unit=12)
write(cF,'(A)') ADJUSTL(TRIM(cOutputFolder))//'/meshDir/file.prj'
write(12,'(A)')  'Merged_'//adjustl(trim(cProjectGridFile))

minCoor=+1d30
maxCoor=-1d30
do ivt=1,mg_NewMesh%level(1)%nvt
 P = mg_NewMesh%level(1)%dcorvg(:,ivt)
 if (minCoor(1).gt.P(1)) minCoor(1) = P(1)
 if (minCoor(2).gt.P(2)) minCoor(2) = P(2)
 if (minCoor(3).gt.P(3)) minCoor(3) = P(3)
 if (maxCoor(1).lt.P(1)) maxCoor(1) = P(1)
 if (maxCoor(2).lt.P(2)) maxCoor(2) = P(2)
 if (maxCoor(3).lt.P(3)) maxCoor(3) = P(3)
end do

iminCoor =0
imaxCoor =0
DO i1=1,3
 do ivt=1,mg_NewMesh%level(1)%nvt
  P = mg_NewMesh%level(1)%dcorvg(:,ivt)
  if (abs(minCoor(i1)-P(i1)).lt.1e-3) iminCoor(i1) = iminCoor(i1) + 1
  if (abs(maxCoor(i1)-P(i1)).lt.1e-3) imaxCoor(i1) = imaxCoor(i1) + 1
 end do
END DO
! write(*,*) minCoor,iminCoor
! write(*,*) maxCoor,imaxCoor

write(cF,'(A,I0,A)') 'X-.par'
write(12,'(A)')  ADJUSTL(TRIM(cF))
write(cF,'(A,I0,A)') ADJUSTL(TRIM(cOutputFolder))//'/meshDir/X-.par'
open(file=ADJUSTL(TRIM(cF)),unit=11)
if (mySigma%DIE_SymmetryBC(1)) then
 write(iunit,'(I0,A)') iminCoor(1),' '//'Symmetry100'
else
 write(iunit,'(I0,A)') iminCoor(1),' '//'Wall'
end if
! write(11,*) iminCoor(1), 'Wall'
write(11,'(A,I0,4ES12.4,A)') "'",4,1.0,0.0,0.0,-minCoor(1),"'"
do ivt=1,mg_NewMesh%level(1)%nvt
 P = mg_NewMesh%level(1)%dcorvg(:,ivt)
 if (abs(minCoor(1)-P(1)).lt.1e-3) then
  write(11,*) ivt
 end if
end do
close(11)

write(cF,'(A,I0,A)') 'X+.par'
write(12,'(A)')  ADJUSTL(TRIM(cF))
write(cF,'(A,I0,A)') ADJUSTL(TRIM(cOutputFolder))//'/meshDir/X+.par'
open(file=ADJUSTL(TRIM(cF)),unit=11)
if (mySigma%DIE_SymmetryBC(4)) then
 write(iunit,'(I0,A)') imaxCoor(1),' '//'Symmetry100'
else
 write(iunit,'(I0,A)') imaxCoor(1),' '//'Wall'
end if
! write(11,*) imaxCoor(1), 'Wall'
write(11,'(A,I0,4ES12.4,A)') "'",4,1.0,0.0,0.0,-maxCoor(1),"'"
do ivt=1,mg_NewMesh%level(1)%nvt
 P = mg_NewMesh%level(1)%dcorvg(:,ivt)
 if (abs(maxCoor(1)-P(1)).lt.1e-3) then
  write(11,*) ivt
 end if
end do
close(11)

write(cF,'(A,I0,A)') 'Y-.par'
write(12,'(A)')  ADJUSTL(TRIM(cF))
write(cF,'(A,I0,A)') ADJUSTL(TRIM(cOutputFolder))//'/meshDir/Y-.par'
open(file=ADJUSTL(TRIM(cF)),unit=11)
if (mySigma%DIE_SymmetryBC(2)) then
 write(iunit,'(I0,A)') iminCoor(2),' '//'Symmetry010'
else
 write(iunit,'(I0,A)') iminCoor(2),' '//'Wall'
end if
! write(11,*) iminCoor(2), 'Wall'
write(11,'(A,I0,4ES12.4,A)') "'",4,0.0,1.0,0.0,-minCoor(2),"'"
do ivt=1,mg_NewMesh%level(1)%nvt
 P = mg_NewMesh%level(1)%dcorvg(:,ivt)
 if (abs(minCoor(2)-P(2)).lt.1e-3) then
  write(11,*) ivt
 end if
end do
close(11)

write(cF,'(A,I0,A)') 'Y+.par'
write(12,'(A)')  ADJUSTL(TRIM(cF))
write(cF,'(A,I0,A)') ADJUSTL(TRIM(cOutputFolder))//'/meshDir/Y+.par'
open(file=ADJUSTL(TRIM(cF)),unit=11)
if (mySigma%DIE_SymmetryBC(2)) then
 write(iunit,'(I0,A)') imaxCoor(2),' '//'Symmetry010'
else
 write(iunit,'(I0,A)') imaxCoor(2),' '//'Wall'
end if
! write(11,*) imaxCoor(2), 'Wall'
write(11,'(A,I0,4ES12.4,A)') "'",4,0.0,1.0,0.0,-maxCoor(2),"'"
do ivt=1,mg_NewMesh%level(1)%nvt
 P = mg_NewMesh%level(1)%dcorvg(:,ivt)
 if (abs(maxCoor(2)-P(2)).lt.1e-3) then
  write(11,*) ivt
 end if
end do
close(11)

write(cF,'(A,I0,A)') 'Z-.par'
write(12,'(A)')  ADJUSTL(TRIM(cF))
write(cF,'(A,I0,A)') ADJUSTL(TRIM(cOutputFolder))//'/meshDir/Z-.par'
open(file=ADJUSTL(TRIM(cF)),unit=11)
write(11,*) iminCoor(3), 'Wall'
write(11,'(A,I0,4ES12.4,A)') "'",4,0.0,0.0,1.0,-minCoor(3),"'"
do ivt=1,mg_NewMesh%level(1)%nvt
 P = mg_NewMesh%level(1)%dcorvg(:,ivt)
 if (abs(minCoor(3)-P(3)).lt.1e-3) then
  write(11,*) ivt
 end if
end do
close(11)

write(cF,'(A,I0,A)') 'Z+.par'
write(12,'(A)')  ADJUSTL(TRIM(cF))
write(cF,'(A,I0,A)') ADJUSTL(TRIM(cOutputFolder))//'/meshDir/Z+.par'
open(file=ADJUSTL(TRIM(cF)),unit=11)
write(11,*) imaxCoor(3), 'Wall'
write(11,'(A,I0,4ES12.4,A)') "'",4,0.0,0.0,1.0,-maxCoor(3),"'"
do ivt=1,mg_NewMesh%level(1)%nvt
 P = mg_NewMesh%level(1)%dcorvg(:,ivt)
 if (abs(maxCoor(3)-P(3)).lt.1e-3) then
  write(11,*) ivt
 end if
end do
close(11)

do iParam=1,ParIndex
 n = 0
 dc = 0d0
 do ivt=1,mg_NewMesh%level(1)%nvt
  if (mg_NewMesh%level(1)%knpr(ivt).eq.-iParam) then
   n = n + 1
   dc = dc + mg_NewMesh%level(1)%dcorvg(:,ivt)
!   if (iParam.eq.80) write(*,'(3ES12.4)') mg_NewMesh%level(1)%dcorvg(:,ivt)
  end if
 end do
 
 CALL RANDOM_NUMBER(myRand)
 myRand = 2d0*myRand-1d0
 myRand = dr*(-0.05d0)*myRand
 
 if (n.eq.26) then
  dc = dc/dble(n)
  write(cF,'(A,I0,A)') 'SPHERE_',iParam,'.par'
  write(12,'(A)')  ADJUSTL(TRIM(cF))
  write(cF,'(A,I0,A)') ADJUSTL(TRIM(cOutputFolder))//'/meshDir/SPHERE_',iParam,'.par'
  open(file=ADJUSTL(TRIM(cF)),unit=11)
  write(11,*) n, 'BALL'
  write(11,'(A,I0,7ES12.4,A)') "'",7,dc,dr,1.0,1.0,1.0,"'"
!  write(11,'(A,I0,7ES12.4,A)') "'",7,dc+myRand,dr,1.0,1.0,1.0,"'"
  do ivt=1,mg_NewMesh%level(1)%nvt
   if (mg_NewMesh%level(1)%knpr(ivt).eq.-iParam) then
    write(11,*) ivt
   end if
  end do
  close(11)
 end if
end do

close(12)

write(cF,'(A,I0,A)') ADJUSTL(TRIM(cOutputFolder))//'/meshDir/PATCHID.ele'
open(file=ADJUSTL(TRIM(cF)),unit=11)
write(11,*) mg_NewMesh%level(1)%nel, 'PATCHID'
write(11,'(A)') "' '"
do iel=1,mg_mesh%level(ilev)%nel
  do i=1,myRF(iel)%nUniqueElems
   write(11, *)  myRF(iel)%patchID 
  end do
end do
close(11)

 CONTAINS

RECURSIVE SUBROUTINE SetRecKNPR(IEL8)
implicit none
INTEGER IEL8
integer i8,ivt8,iat8,jel8
logical bContinue1,bContinue2

DO iat8=1,6
 jel8 = mg_NewMesh%level(1)%kadj(iat8,iel8)
 bContinue1 = .false.
 bContinue2 = .false.
 if (JEL8.gt.0) THEN
  do i8=1,8
   ivt8 = mg_NewMesh%level(1)%kvert(i8,jel8)
   if (mg_NewMesh%level(1)%knpr(ivt8).eq.-ParIndex) THEN
    bContinue1 = .true.
   end if
   if (mg_NewMesh%level(1)%knpr(ivt8).gt.0) THEN
    bContinue2 = .true.
   end if
  end do
!  write(*,*) iat8,bContinue1, bContinue1, JEL8
  if (bContinue1.and.bContinue2) then
   do i8=1,8
    ivt8 = mg_NewMesh%level(1)%kvert(i8,jel8)
    if (mg_NewMesh%level(1)%knpr(ivt8).gt.0) THEN
     nParam = nParam + 1
     mg_NewMesh%level(1)%knpr(ivt8) = -ParIndex
    end if
   end do
   CALL SetRecKNPR(jel8)
  end if
 end if
end do

END SUBROUTINE SetRecKNPR
 
END SUBROUTINE Output_ParticleMergedRefTriMeshPar
! ----------------------------------------------
SUBROUTINE Output_ParticleFineVTK

IMPLICIT NONE
INTEGER nel,nvt
INTEGER ive,ivt,ioffset,ibnds
INTEGER :: iunit=123
CHARACTER*(100) filename
real*8 dBCValue 

nel = mg_ParticleMesh%level(ilev)%nel
nvt = mg_ParticleMesh%level(ilev)%nvt

write(*,*) 'nel=',nel

filename=" "
WRITE(filename(1:),'(A)') adjustl(trim(cOutputFolder))//"/FineParticleMesh.vtu"

WRITE(*,'(104("="))') 
WRITE(*,*) "Outputting vtk file into ",filename

OPEN (UNIT=iunit,FILE=filename)

write(iunit, *)"<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">"
write(iunit, *)"  <UnstructuredGrid>"
write(iunit, *)"    <Piece NumberOfPoints=""",nvt,""" NumberOfCells=""",nel,""">"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the node field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <PointData>"

 write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","ID",""" format=""ascii"">"
 do ivt=1,nvt
  write(iunit, '(A,E16.7)')"        ",REAL(ivt)
 end do
 write(iunit, *)"        </DataArray>"
 
 write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","BC",""" format=""ascii"">"
 do ivt=1,nvt
  dBCValue = 0d0
  DO iBnds=1,nBnds
   if (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).eq.'BALL') THEN
    if (myParBndr(ibnds)%Bndr(ilev)%Vert(ivt)) then
     dBCValue = 1d0
    end if
   end if
  end do
  write(iunit, '(A,E16.7)')"        ",REAL(dBCValue)

 end do
 write(iunit, *)"        </DataArray>"

 write(iunit, '(A)')"    </PointData>"

write(iunit, '(A)')"      <Points>"
write(iunit, '(A)')"        <DataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3"" format=""ascii"" RangeMin=""0"" RangeMax=""1.0"">"
do ivt=1,nvt
 write(iunit,'(A10,3E16.7)')"          ",REAL(mg_ParticleMesh%level(ilev)%dcorvg(1,ivt)),REAL(mg_ParticleMesh%level(ilev)%dcorvg(2,ivt)),REAL(mg_ParticleMesh%level(ilev)%dcorvg(3,ivt))
end do
write(iunit, *)"        </DataArray>"
write(iunit, *)"      </Points>"

write(iunit, '(A)')"    <CellData>"
write(iunit, '(A,A,A)')"        <DataArray type=""Int32"" Name=""","OrigID",""" format=""ascii"">"
do ivt=1,nel
 write(iunit, '(A,I10)')"        ",ivt
end do
write(iunit, *)"        </DataArray>"

write(iunit, '(A)')"    </CellData>"

write(iunit, *)"      <Cells>"
write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""connectivity"" format=""ascii"" RangeMin=""0"" RangeMax=""",nel-1,""">"
do ive=1,nel   
 write(iunit, '(8I10)')mg_ParticleMesh%level(ilev)%kvert(1,ive)-1,mg_ParticleMesh%level(ilev)%kvert(2,ive)-1,mg_ParticleMesh%level(ilev)%kvert(3,ive)-1,mg_ParticleMesh%level(ilev)%kvert(4,ive)-1,&
                       mg_ParticleMesh%level(ilev)%kvert(5,ive)-1,mg_ParticleMesh%level(ilev)%kvert(6,ive)-1,mg_ParticleMesh%level(ilev)%kvert(7,ive)-1,mg_ParticleMesh%level(ilev)%kvert(8,ive)-1
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""offsets"" format=""ascii"" RangeMin=""8"" RangeMax=""",8*nel,""">"
ioffset=nel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')ive*8,(ive+1)*8,(ive+2)*8,(ive+3)*8,(ive+4)*8,(ive+5)*8,(ive+6)*8,(ive+7)*8
end do

do ive=ioffset+1,nel
 write(iunit, '(I10)')ive*8
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A)')"        <DataArray type=""UInt8"" Name=""types"" format=""ascii"" RangeMin=""12"" RangeMax=""12"">"
ioffset=nel/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')12,12,12,12,12,12,12,12
end do
do ive=ioffset+1,nel
 write(iunit, '(I10)')12
end do
write(iunit, '(A)')"        </DataArray>"
 
write(iunit, *)"      </Cells>"
write(iunit, *)"    </Piece>"
   
write(iunit, *)"  </UnstructuredGrid>"
write(iunit, *)"</VTKFile>"
close(iunit)

END SUBROUTINE Output_ParticleFineVTK

END Module MeshRefOutput
