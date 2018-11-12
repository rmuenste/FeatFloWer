MODULE MeshProcDef
USE types, ONLY: tMultiMesh,tMesh
USE Parametrization, ONLY: InitParametrization_STRCT,ParametrizeBndryPoints_STRCT,&
    DeterminePointParametrization_STRCT,ProlongateParametrization_STRCT,myParBndr,nBnds
USE var_QuadScalar
USE PP3D_MPI, ONLY:myid,master
USE MESH_Structures
use iniparser

implicit none

INTEGER lTriOutputLevel,lVTUOutputLevel

CHARACTER*(200) :: cOutputFolder
LOGICAL :: bA_MD=.false.
LOGICAL :: bPDE_MD=.false.
Logical :: bDefTensor = .true.

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GetParameters
CHARACTER*(10) :: cParam

!  CHARACTER*(200) :: cmd
 LOGICAL bExist

 OPEN(1,file='param.cfg')

 READ(1,*) cProjectFolder
 WRITE(*,*) adjustl(trim(cProjectFolder))
 READ(1,*) cProjectFile
 cProjectFile = adjustl(trim(cProjectFolder))//adjustl(trim(cProjectFile))
 WRITE(*,*) adjustl(trim(cProjectFile))
 READ(1,*) cProjectGridFile
 cProjectGridFile = adjustl(trim(cProjectFolder))//adjustl(trim(cProjectGridFile))
 WRITE(*,*) adjustl(trim(cProjectGridFile))
 READ(1,*) mg_Mesh%nlmax
 mg_Mesh%nlmin = 1
 WRITE(*,*) 'Min and Max levels: ', mg_Mesh%nlmin,mg_Mesh%nlmax
 mg_Mesh%nlmin = 1
 READ(1,*) cParam
 call inip_toupper_replace(cParam)
 IF (trim(adjustl(cParam)).eq.'ON') THEN
  bA_MD=.TRUE.
  WRITE(*,*) 'Alghebraic Mesh Smoother is ON!'
 ELSE
  bA_MD=.FALSE.
  WRITE(*,*) 'Alghebraic Mesh Smoother is OFF!'
 END IF
 READ(1,*) cParam
 call inip_toupper_replace(cParam)
 bPDE_MD=.FALSE.
 IF (trim(adjustl(cParam)).eq.'LAPLACE'.OR.trim(adjustl(cParam)).eq.'ON') THEN
  bPDE_MD=.TRUE.
  bDefTensor = .false.
  WRITE(*,*) 'PDE based Mesh Smoother with LAPLACE!'
 END IF
 IF (trim(adjustl(cParam)).eq.'DEFTENSOR') THEN
  bPDE_MD=.TRUE.
  bDefTensor = .true.
  WRITE(*,*) 'PDE based Mesh Smoother with DEFTENSOR!'
 END IF
 IF (.NOT.bPDE_MD) THEN
  WRITE(*,*) 'PDE based Mesh Smoother is OFF!'
 END IF

 READ(1,*) nUmbrellaSteps
 WRITE(*,*) '# of Umbrella Steps: ', nUmbrellaSteps
 READ(1,*) dCGALtoRealFactor
 WRITE(*,'(A,ES12.4)') ' CGAL Scaling factor: ', dCGALtoRealFactor
 READ(1,*) cOutputFolder
 WRITE(*,*) 'Output Folder: "'//adjustl(trim(cOutputFolder))//'"'
 READ(1,*) lTriOutputLevel
 WRITE(*,*) 'Outputlevel for the ".tri" file: ', lTriOutputLevel
 READ(1,*) lVTUOutputLevel
 WRITE(*,*) 'Outputlevel for the ".vtu" file: ', lVTUOutputLevel

 CLOSE(1)

 INQUIRE(DIRECTORY=adjustl(trim(cOutputFolder)),EXIST=bExist)
 if (.not.bExist) then
!   cmd='mkdir '//adjustl(trim(cOutputFolder))
  CALL system('mkdir '//adjustl(trim(cOutputFolder)))
 end if

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
WRITE(filename(1:),'(A)') adjustl(trim(cOutputFolder))//"/mesh.vtu"

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
SUBROUTINE SeqUmbrella(ilev,nProjStep)
INTEGER ilev,nProjStep
!!
INTEGER ndof
REAL*8 , ALLOCATABLE :: a1(:),a2(:),a3(:),a4(:),a5(:),a6(:)

ndof = mg_mesh%level(ilev)%nvt + mg_mesh%level(ilev)%net + mg_mesh%level(ilev)%nat + mg_mesh%level(ilev)%nel

ALLOCATE(a1(ndof))
ALLOCATE(a2(ndof))
ALLOCATE(a3(ndof))
ALLOCATE(a4(ndof))
ALLOCATE(a5(ndof))
ALLOCATE(a6(ndof))

CALL SeqEdgeRunner(a1,a2,a3,a4,a5,a6,&
    mg_mesh,ilev,&
    mg_mesh%level(ilev)%dcorvg,&
    mg_mesh%level(ilev)%kvert,&
    mg_mesh%level(ilev)%kedge,&
    mg_mesh%level(ilev)%nel,&
    mg_mesh%level(ilev)%nvt,&
    mg_mesh%level(ilev)%net,&
    nProjStep)
                   
DEALLOCATE(a1,a2,a3,a4,a5,a6)
  
END SUBROUTINE SeqUmbrella
!----------------------------------------------------------
SUBROUTINE SeqEdgeRunner(f,x,y,z,w,v,mgMesh,ilevel,&
  dcorvg,kvert,kedge,nel,nvt,net,nProjStep)
 
IMPLICIT NONE

REAL*8 f(*),x(*),y(*),z(*),w(*),v(*)
INTEGER nel,nvt,net,nProjStep
REAL*8, intent(inout), dimension(:,:) :: dcorvg
integer, intent(in), dimension(:,:) :: kvert,kedge
        
integer :: ilevel
type(tMultiMesh) :: mgMesh

INTEGER NeighE(2,12)
DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
INTEGER i,j,k,ivt1,ivt2,iProjStep,iaux,iel
REAL*8 WeightE,P1(3),P2(3),daux2,daux1,PX,PY,PZ
REAL*8 :: dOmega = 0.25d0
REAL*8 :: DIST
REAL*8, ALLOCATABLE :: myVol(:)

DO k=nvt+1,nvt+net
v(k) = 1d0
END DO

ALLOCATE(myVol(nel+1))

DO iProjStep=1,nProjStep

myVol = 0e0
CALL  SETARE(myVol,nel,kvert,dcorvg)

f(1:nvt) = 0d0
w(1:nvt) = 0d0

DO iel=1,nel
  DO i=1,8
    j = kvert(i,iel)
    f(j) = f(j) + abs(myVol(iel))
    w(j) = w(j) + 1d0
  END DO
END DO

DO i=1,nvt
  f(i) = f(i)/w(i)
END DO

! DO i=1,nvt
!   PX = dcorvg(1,i)
!   PY = dcorvg(2,i)
!   PZ = dcorvg(3,i)
!   f(i) = 1d0
! END DO

x(1:nvt) = 0d0
y(1:nvt) = 0d0
z(1:nvt) = 0d0
w(1:nvt) = 0d0

k=1
DO i=1,nel
  DO j=1,12
  IF (k.eq.kedge(j,i)) THEN
    ivt1 = kvert(NeighE(1,j),i)
    ivt2 = kvert(NeighE(2,j),i)
    P1(:) = dcorvg(:,ivt1)
    P2(:) = dcorvg(:,ivt2)

    daux1 = ABS(f(ivt1))
    daux2 = ABS(f(ivt2))
    WeightE = 1d0/(v(nvt + k))

    x(ivt1) = x(ivt1) + WeightE*P2(1)*daux2
    y(ivt1) = y(ivt1) + WeightE*P2(2)*daux2
    z(ivt1) = z(ivt1) + WeightE*P2(3)*daux2
    w(ivt1) = w(ivt1) + WeightE*daux2

    x(ivt2) = x(ivt2) + WeightE*P1(1)*daux1
    y(ivt2) = y(ivt2) + WeightE*P1(2)*daux1
    z(ivt2) = z(ivt2) + WeightE*P1(3)*daux1
    w(ivt2) = w(ivt2) + WeightE*daux1

    k = k + 1
  END IF
  END DO
END DO

DO i=1,nvt
 PX = x(i)/w(i)
 PY = y(i)/w(i)
 PZ = z(i)/w(i)
 dcorvg(1,i) = MAX(0d0,(1d0-dOmega))*dcorvg(1,i) + dOmega*PX
 dcorvg(2,i) = MAX(0d0,(1d0-dOmega))*dcorvg(2,i) + dOmega*PY
 dcorvg(3,i) = MAX(0d0,(1d0-dOmega))*dcorvg(3,i) + dOmega*PZ
 END DO

 CALL ParametrizeBndryPoints_STRCT(mgMesh,ilevel)

END DO

END SUBROUTINE SeqEdgeRunner
!----------------------------------------------------------
END MODULE MeshProcDef
