MODULE MeshProcDef
USE types, ONLY: tMultiMesh,tMesh
USE Parametrization, ONLY: InitParametrization_STRCT,ParametrizeBndryPoints_STRCT,&
    DeterminePointParametrization_STRCT,ProlongateParametrization_STRCT,myParBndr,nBnds
USE var_QuadScalar
USE PP3D_MPI, ONLY:myid,master

implicit none

!type(tMultiMesh) mg_Mesh
!CHARACTER*(200) :: cProjectFolder,cProjectFile,cMeshFile

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GetParameters

 OPEN(1,file='param.cfg')

 READ(1,*) cProjectFolder
 WRITE(*,*) adjustl(trim(cProjectFolder))
 READ(1,*) cProjectFile
 cProjectFile = adjustl(trim(cProjectFolder))//adjustl(trim(cProjectFile))
 WRITE(*,*) adjustl(trim(cProjectFile))
 READ(1,*) cProjectGridFile
 cProjectGridFile = adjustl(trim(cProjectFolder))//adjustl(trim(cProjectGridFile))
 WRITE(*,*) adjustl(trim(cProjectGridFile))

 CLOSE(1)

END SUBROUTINE GetParameters
!----------------------------------------------------------
SUBROUTINE Output_VTK

IMPLICIT NONE
INTEGER nel,nvt
INTEGER ive,ivt,ioffset
INTEGER :: iunit=123
CHARACTER*(100) filename

nel = mg_mesh%level(ilev)%nel
nvt = mg_mesh%level(ilev)%nvt

filename=" "
WRITE(filename(1:),'(A)') "mesh.vtu"

WRITE(*,'(104("="))') 
WRITE(*,*) "Outputting vtk file into ",filename

OPEN (UNIT=iunit,FILE=filename,buffered="yes")

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

 write(iunit, '(A)')"    </PointData>"

write(iunit, '(A)')"      <Points>"
write(iunit, '(A)')"        <DataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3"" format=""ascii"" RangeMin=""0"" RangeMax=""1.0"">"
do ivt=1,nvt
 write(iunit,'(A10,3E16.7)')"          ",REAL(mg_mesh%level(ilev)%dcorvg(1,ivt)),REAL(mg_mesh%level(ilev)%dcorvg(2,ivt)),REAL(mg_mesh%level(ilev)%dcorvg(3,ivt))
end do
write(iunit, *)"        </DataArray>"
write(iunit, *)"      </Points>"

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
SUBROUTINE GetFileList()
 CHARACTER(LEN=200) :: string,cFile
 CHARACTER cWD*200,cSub*200
 INTEGER lenCommand,i,iPos,LenStr,iEnd,iBnds
 
 nBnds = 0
 OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cProjectFolder))//'/'//ADJUSTL(TRIM(cProjectFile)))
 DO
  READ(1,FMT='(200A)',IOSTAT=iEnd) string
  IF (iEnd.EQ.-1) EXIT
  LenStr = LEN(ADJUSTL(TRIM(string)))
  IF (LenStr.gt.4) THEN
   cFile = ADJUSTL(TRIM(string))
   IF (cFile(LenStr-3:LenStr).EQ.".par") nBnds = nBnds + 1
  END IF
 END DO

 ALLOCATE (myParBndr(nBnds))

 Write(*,*)'Number of boundary parametrizations: ',nBnds
 DO iBnds = 1, nBnds
  ALLOCATE (myParBndr(iBnds)%Bndr(mg_Mesh%NLMIN:mg_Mesh%NLMAX+1))
 END DO

 nBnds = 0
 REWIND(1)

 DO
  READ(1,FMT='(200A)',IOSTAT=iEnd) string
  IF (iEnd.EQ.-1) EXIT
  LenStr = LEN(ADJUSTL(TRIM(string)))
  IF (LenStr.gt.4) THEN
   cFile = ADJUSTL(TRIM(string))
   IF (cFile(LenStr-3:LenStr).EQ.".par") THEN
    Write(*,*)'Boundary parametrization: "'//ADJUSTL(TRIM(cFile))//'"'
    nBnds = nBnds + 1
    READ(cFile(1:LenStr-4),"(A)") myParBndr(nBnds)%Names
   END IF
   IF (cFile(LenStr-3:LenStr).EQ.".tri") THEN
    READ(cFile(1:LenStr),"(A)") cGridFileName
    Write(*,*)'Mesh file: "'//ADJUSTL(TRIM(cFile))//'"'
   END IF
  END IF
 END DO
 CLOSE(1)

END SUBROUTINE GetFileList
!----------------------------------------------------------
END MODULE MeshProcDef
