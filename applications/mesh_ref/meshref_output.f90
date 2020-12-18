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
do iel=1,nnvt
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
WRITE(*,*) "Outputting actual Coarse mesh into: Ref'"//ADJUSTL(TRIM(cf))//"'"
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

OPEN(UNIT=2,FILE=ADJUSTL(TRIM(cOutputFolder))//'/'//adjustl(trim(cShortProjectFile)))
!WRITE(cf,'(A11,I3.3,A4)') 'cMESH_',iO, '.tri'
WRITE(2,'(A)') adjustl(trim(cProjectGridFile))
 
DO iBnds = 1, nBnds
 cf = ' '
 WRITE(cf,'(A)') ADJUSTL(TRIM(cOutputFolder))//"/"//ADJUSTL(TRIM(myParBndr(iBnds)%Names))//".par"
 WRITE(2,'(A)') ADJUSTL(TRIM(myParBndr(iBnds)%Names))//".par"
 WRITE(*,*) "Outputting actual parametrization into: '"//ADJUSTL(TRIM(cf))//"'"
 OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cf)))
 j=0
 DO i=1,mg_mesh%level(ilev)%nvt
  IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
   j = j + 1
  END IF
 END DO
 WRITE(1,'(I8,A)') j," "//myParBndr(iBnds)%Types
 WRITE(1,'(A)')    "'"//ADJUSTL(TRIM(myParBndr(iBnds)%Parameters))//"'"
 DO i=1,mg_mesh%level(ilev)%nvt
  IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
   WRITE(1,'(I8,A)') i
  END IF
 END DO
 CLOSE(1)
END DO
CLOSE(2)

END SUBROUTINE Output_TriMesh
! ----------------------------------------------
SUBROUTINE Output_MergedRefTriMeshPar()

USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:mg_mesh,myBoundary
USE Parametrization, ONLY : myParBndr,nBnds
IMPLICIT NONE
INTEGER i,j,iloc,i1,i2,i3,i4,i5,i6
CHARACTER cf*(256),ctxt*(256)
real*8 P0(3),box(3),bound(3,2)
real*8 :: dEps=1d-4

OPEN(file='Out/param.txt',unit=5)
read(5,'(a)') ctxt
write(*,*) "'"//adjustl(trim(ctxt))//"'"
iloc = INDEX(ctxt,"=")
write(*,*) iloc
read(ctxt(iloc+1:),*) P0

read(5,'(a)') ctxt
iloc = INDEX(ctxt,"=")
write(*,*) iloc
read(ctxt(iloc+1:),*) box

close(5)

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

write(11,'(I0,A)') i1,' Wall'
write(11,'(A,I0," ",4F10.4,A)') "'",4,1.0, 0.0, 0.0, -MeshOutputScaleFactor*bound(1,1), "'"
write(12,'(I0,A)') i2,' Wall'
write(12,'(A,I0," ",4F10.4,A)') "'",4,1.0, 0.0, 0.0, -MeshOutputScaleFactor*bound(1,2), "'"
write(13,'(I0,A)') i3,' Wall'
write(13,'(A,I0," ",4F10.4,A)') "'",4,0.0, 1.0, 0.0, -MeshOutputScaleFactor*bound(2,1), "'"
write(14,'(I0,A)') i4,' Inflow-2'
write(14,'(A,I0," ",4F10.4,A)') "'",4,0.0, 1.0, 0.0, -MeshOutputScaleFactor*bound(2,2), "'"
write(15,'(I0,A)') i5,' Inflow-1'
write(15,'(A,I0," ",4F10.4,A)') "'",4,0.0, 0.0, 1.0, -MeshOutputScaleFactor*bound(3,1), "'"
write(16,'(I0,A)') i6,' Outflow'
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


END SUBROUTINE Output_MergedRefTriMeshPar
! ----------------------------------------------
SUBROUTINE Output_MergedRefTriMesh()

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

END SUBROUTINE Output_MergedRefTriMesh

END Module MeshRefOutput