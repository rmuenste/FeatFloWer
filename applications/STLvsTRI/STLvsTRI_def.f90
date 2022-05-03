MODULE mSTLvsTRI

implicit none

type tOFFmesh
 REAL*8,allocatable :: coor(:,:)
 integer, allocatable :: kvert(:,:)
 integer nel,nvt
end type tOFFmesh
type (tOFFmesh) :: OFFmesh 
CHARACTER*(200) :: cProjectFolder,cShortProjectFile,cOFFMeshFile,cProjectGridFile,cAreaIntensityFile
integer iOK,iTriang,iRecLevel
TYPE tBoxMesh
 REAL*8 extent(3,2)
 INTEGER division(3)
END TYPE tBoxMesh
TYPE(tBoxMesh) BoxMesh

TYPE tTriMesh
 REAL*8, allocatable :: x(:),y(:),z(:)
 REAL*8, allocatable :: d(:,:,:) 
END TYPE tTriMesh
TYPE(tTriMesh) TriMesh

Type triplet
 REAl*8 Q(3)
end type triplet

 CONTAINS
 
SUBROUTINE readOFFMesh(cF)
CHARACTER*(*) cF
integer i,iaux

OPEN(File=trim(cF),unit=1)
read(1,*)
read(1,*) OFFmesh%nvt,OFFmesh%nel

allocate(OFFmesh%coor(OFFmesh%nvt,3))
allocate(OFFmesh%kvert(OFFmesh%nel,3))

DO i=1,OFFmesh%nvt
 read(1,*) OFFmesh%coor(i,:)
END DO

DO i=1,OFFmesh%nel
 read(1,*) iaux,OFFmesh%kvert(i,:)
END DO
OFFmesh%kvert =  OFFmesh%kvert + 1

close(1)

END SUBROUTINE readOFFMesh
!
! -------------------------------------------------------------------------------
!
SUBROUTINE CheckForIntersection
REAL*8  P1(3),P2(3),P3(3)
integer i,iDir

allocate(TriMesh%x(BoxMesh%Division(1)))
allocate(TriMesh%y(BoxMesh%Division(2)))
allocate(TriMesh%z(BoxMesh%Division(3)))
allocate(TriMesh%d(BoxMesh%Division(1),BoxMesh%Division(2),BoxMesh%Division(3)))
TriMesh%x = 0d0
TriMesh%y = 0d0
TriMesh%z = 0d0
TriMesh%d = 0d0
DO i=1,BoxMesh%division(1)
 TriMesh%x(i) = BoxMesh%Extent(1,1) + dble(i-1)*(BoxMesh%Extent(1,2)-BoxMesh%Extent(1,1))/(BoxMesh%division(1)-1)
END DO
DO i=1,BoxMesh%division(2)
 TriMesh%y(i) = BoxMesh%Extent(2,1) + dble(i-1)*(BoxMesh%Extent(2,2)-BoxMesh%Extent(2,1))/(BoxMesh%division(2)-1)
END DO
DO i=1,BoxMesh%division(3)
 TriMesh%z(i) = BoxMesh%Extent(3,1) + dble(i-1)*(BoxMesh%Extent(3,2)-BoxMesh%Extent(3,1))/(BoxMesh%division(3)-1)
END DO


OPEN(FILE='EXP.stl',unit=2)
WRITE(2,'(A)') "solid Created by STLtoTRI"

DO iTriang=1,OFFmesh%nel
 P1 = OFFmesh%coor(OFFmesh%kvert(iTriang,1),:)
 P2 = OFFmesh%coor(OFFmesh%kvert(iTriang,2),:)
 P3 = OFFmesh%coor(OFFmesh%kvert(iTriang,3),:)
 iDir = 1
 CALL DecomposeTriangle(P1,P2,P3,iDir)
END DO

WRITE(2,'(A)')  "endsolid Created by STLtoTRI"
close(2)

END SUBROUTINE CheckForIntersection
!
! -------------------------------------------------------------------------------
!
RECURSIVE SUBROUTINE DecomposeTriangle(P1,P2,P3,iDir)
REAL*8, INTENT(in) :: P1(3),P2(3),P3(3)
integer iDir
REAL*8 X(3),dmin,dmax,dLowerBorder,dIncrement
integer iCase,jCase,i,iMin,newDir,iLevel
REAL*8 :: Q0(3),Q1(3),Q2(3),Q3(3),R2(3),R3(3),dFact1,dFact2,dLevel,dLevel0
Type (triplet) t(3)
logical :: bDo

if (iDir.gt.3) THEN
 call SaveTriangle(P1,P2,P3)
 RETURN
end if

t(1)%Q = P1
t(2)%Q = P2
t(3)%Q = P3

iMin = -1
dmin = +1d8

DO i=1,3
 if (t(i)%Q(iDir).lt.dmin) THEN
  iMin = i
  dmin = t(i)%Q(iDir)
 end if
end do

IF (iMin.ne.1) THEN
 X = t(1)%Q
 t(1)%Q = t(iMin)%Q
 t(iMin)%Q = X
END IF

IF (t(2)%Q(iDir).gt.t(3)%Q(iDir)) then
 X = t(3)%Q
 t(3)%Q = t(2)%Q
 t(2)%Q = X
END IF


iCase = 1
IF (t(1)%Q(iDir).eq.t(2)%Q(iDir)) iCase = 2
IF (t(2)%Q(iDir).eq.t(3)%Q(iDir)) iCase = 3
jCase = iCase

call FindLowerBorder(t(1)%Q(iDir),dLowerBorder,iDir)

dIncrement = (BoxMesh%extent(iDir,2)-BoxMesh%extent(iDir,1))/DBLE(BoxMesh%division(iDir)-1)

! write(*,*) "--------------------------"
! write(*,*) t(1)%Q
! write(*,*) t(2)%Q
! write(*,*) t(3)%Q
! write(*,*) X
! write(*,*) iDir,iMin
! write(*,*) "--------------------------"

! write(*,'(4ES12.4,I5)') t(1)%Q(iDir),t(2)%Q(iDir),t(3)%Q(iDir),dLowerBorder,iDir

if (iCase.eq.1) then
 if (dLowerBorder+dIncrement.ge.t(3)%Q(iDir)) THEN
  newDir = iDir + 1
  Q1 = t(1)%Q
  Q2 = t(2)%Q
  Q3 = t(3)%Q
  CALL DecomposeTriangle(Q1,Q2,Q3,newDir)
 ELSE

  newDir = iDir + 1
  
  dLevel0 = dLowerBorder
  iLevel  = 0
  
  dLevel = MIN(dLevel0 + DBLE(iLevel+1)*dIncrement,t(2)%Q(iDir))
  IF (dLevel.ge.dLevel0 + DBLE(iLevel+1)*dIncrement) iLevel = iLevel + 1
  
  dFact1 = (dLevel - t(1)%Q(iDir)) / (t(2)%Q(iDir) - t(1)%Q(iDir))
  dFact2 = (dLevel - t(1)%Q(iDir)) / (t(3)%Q(iDir) - t(1)%Q(iDir))
  Q1 = t(1)%Q
  Q2 = t(1)%Q + dFact1*(t(2)%Q-t(1)%Q)
  Q3 = t(1)%Q + dFact2*(t(3)%Q-t(1)%Q)
  CALL DecomposeTriangle(Q1,Q2,Q3,newDir)
!   write(*,*) t(1)%Q(iDir),t(2)%Q(iDir),dLevel,dIncrement
  
  bDo = .true.
  DO while(bDo)
   R2 = Q2
   R3 = Q3
   if (dLevel.lt.t(2)%Q(iDir)) then
    dLevel = MIN(dLevel0 + DBLE(iLevel+1)*dIncrement,t(2)%Q(iDir))
    IF (dLevel.ge.dLevel0 + DBLE(iLevel+1)*dIncrement) iLevel = iLevel + 1
    
    dFact1 = (dLevel - t(1)%Q(iDir)) / (t(2)%Q(iDir) - t(1)%Q(iDir))
    dFact2 = (dLevel - t(1)%Q(iDir)) / (t(3)%Q(iDir) - t(1)%Q(iDir))
    Q2 = t(1)%Q + dFact1*(t(2)%Q-t(1)%Q)
    Q3 = t(1)%Q + dFact2*(t(3)%Q-t(1)%Q)
!     write(*,*) "A",iLevel
   else
    dLevel = MIN(dLevel0 + DBLE(iLevel+1)*dIncrement,t(3)%Q(iDir))
    IF (dLevel.ge.dLevel0 + DBLE(iLevel+1)*dIncrement) iLevel = iLevel + 1
    
    dFact1 = (dLevel - t(2)%Q(iDir)) / (t(3)%Q(iDir) - t(2)%Q(iDir))
    dFact2 = (dLevel - t(1)%Q(iDir)) / (t(3)%Q(iDir) - t(1)%Q(iDir))
    Q2 = t(2)%Q + dFact1*(t(3)%Q-t(2)%Q)
    Q3 = t(1)%Q + dFact2*(t(3)%Q-t(1)%Q)
!     write(*,*) "B",iLEvel
   end if
!   pause
   
   if (dLevel.eq.t(3)%Q(iDir)) THEN
    CALL DecomposeTriangle(R2,R3,Q2,newDir)
    bDo = .false.
   ELSE
!     write(*,*) "++++++++++++++++++++++++++++"
!     write(*,*) r2
!     write(*,*) R3
!     write(*,*) q2
!     write(*,*) "++++++++++++++++++++++++++++"
    CALL DecomposeTriangle(R2,R3,Q2,newDir)
    CALL DecomposeTriangle(R3,Q2,Q3,newDir)
   END IF
  
  END DO
 
 END IF
end if

if (iCase.eq.2) then

 if (dLowerBorder+dIncrement.ge.t(3)%Q(iDir)) THEN
  newDir = iDir + 1
  Q1 = t(1)%Q
  Q2 = t(2)%Q
  Q3 = t(3)%Q
  CALL DecomposeTriangle(Q1,Q2,Q3,newDir)
 ELSE

  newDir = iDir + 1
  
  dLevel0 = dLowerBorder
  iLevel  = 0
  
  Q2 = t(2)%Q
  Q3 = t(1)%Q
  
  bDo = .true.
  DO while(bDo)
   R2 = Q2
   R3 = Q3
   dLevel = MIN(dLevel0 + DBLE(iLevel+1)*dIncrement,t(3)%Q(iDir))
   IF (dLevel.ge.dLevel0 + DBLE(iLevel+1)*dIncrement) iLevel = iLevel + 1
   
   dFact1 = (dLevel - t(2)%Q(iDir)) / (t(3)%Q(iDir) - t(2)%Q(iDir))
   dFact2 = (dLevel - t(1)%Q(iDir)) / (t(3)%Q(iDir) - t(1)%Q(iDir))
   Q2 = t(2)%Q + dFact1*(t(3)%Q-t(2)%Q)
   Q3 = t(1)%Q + dFact2*(t(3)%Q-t(1)%Q)
!     write(*,*) "B",iLEvel
!   pause
   
   if (dLevel.eq.t(3)%Q(iDir)) THEN
    CALL DecomposeTriangle(R2,R3,Q2,newDir)
    bDo = .false.
   ELSE
    CALL DecomposeTriangle(R2,R3,Q2,newDir)
    CALL DecomposeTriangle(R3,Q2,Q3,newDir)
   END IF
  
  END DO
 
 END IF
end if

if (iCase.eq.3) then
 if (dLowerBorder+dIncrement.ge.t(3)%Q(iDir)) THEN
  newDir = iDir + 1
  Q1 = t(1)%Q
  Q2 = t(2)%Q
  Q3 = t(3)%Q
  CALL DecomposeTriangle(Q1,Q2,Q3,newDir)
 ELSE

  newDir = iDir + 1
  
  dLevel0 = dLowerBorder
  iLevel  = 0
  
  dLevel = MIN(dLevel0 + DBLE(iLevel+1)*dIncrement,t(2)%Q(iDir))
  IF (dLevel.ge.dLevel0 + DBLE(iLevel+1)*dIncrement) iLevel = iLevel + 1
  
  dFact1 = (dLevel - t(1)%Q(iDir)) / (t(2)%Q(iDir) - t(1)%Q(iDir))
  dFact2 = (dLevel - t(1)%Q(iDir)) / (t(3)%Q(iDir) - t(1)%Q(iDir))
  Q1 = t(1)%Q
  Q2 = t(1)%Q + dFact1*(t(2)%Q-t(1)%Q)
  Q3 = t(1)%Q + dFact2*(t(3)%Q-t(1)%Q)
  CALL DecomposeTriangle(Q1,Q2,Q3,newDir)
  
  bDo = .true.
  DO while(bDo)
   R2 = Q2
   R3 = Q3
   
   dLevel = MIN(dLevel0 + DBLE(iLevel+1)*dIncrement,t(3)%Q(iDir))
   IF (dLevel.ge.dLevel0 + DBLE(iLevel+1)*dIncrement) iLevel = iLevel + 1
   
   dFact1 = (dLevel - t(1)%Q(iDir)) / (t(2)%Q(iDir) - t(1)%Q(iDir))
   dFact2 = (dLevel - t(1)%Q(iDir)) / (t(3)%Q(iDir) - t(1)%Q(iDir))
   Q2 = t(1)%Q + dFact1*(t(2)%Q-t(1)%Q)
   Q3 = t(1)%Q + dFact2*(t(3)%Q-t(1)%Q)
  
   CALL DecomposeTriangle(R2,R3,Q2,newDir)
   CALL DecomposeTriangle(R3,Q2,Q3,newDir)
   
   if (dLevel.eq.t(3)%Q(iDir)) THEN
    bDo = .false.
   END IF
  
  END DO
 
 END IF
end if

END SUBROUTINE DecomposeTriangle
!
! -------------------------------------------------------------------------------
!
SUBROUTINE SaveTriangle(P1,P2,P3)
REAL*8, INTENT(in) :: P1(3),P2(3),P3(3)
integer j,dir(3)
REAL*8 :: P(3),dnx,dny,dnz,darea,n1(3),n2(3)
REAL*8 :: deps = 1d-8

J = 1
IF (.not.(P1(J).ge.BoxMesh%extent(J,1)-deps.and.P2(J).ge.BoxMesh%extent(J,1)-deps.and.P3(J).ge.BoxMesh%extent(J,1)-deps)) return
IF (.not.(P1(J).le.BoxMesh%extent(J,2)+deps.and.P2(J).le.BoxMesh%extent(J,2)+deps.and.P3(J).le.BoxMesh%extent(J,2)+deps)) return
J = 2
IF (.not.(P1(J).ge.BoxMesh%extent(J,1)-deps.and.P2(J).ge.BoxMesh%extent(J,1)-deps.and.P3(J).ge.BoxMesh%extent(J,1)-deps)) return
IF (.not.(P1(J).le.BoxMesh%extent(J,2)+deps.and.P2(J).le.BoxMesh%extent(J,2)+deps.and.P3(J).le.BoxMesh%extent(J,2)+deps)) return
J = 3
IF (.not.(P1(J).ge.BoxMesh%extent(J,1)-deps.and.P2(J).ge.BoxMesh%extent(J,1)-deps.and.P3(J).ge.BoxMesh%extent(J,1)-deps)) return
IF (.not.(P1(J).le.BoxMesh%extent(J,2)+deps.and.P2(J).le.BoxMesh%extent(J,2)+deps.and.P3(J).le.BoxMesh%extent(J,2)+deps)) return

P = (P1+P2+P3)/3d0

DO J=1,3
 Dir(J) =  INT(DBLE(BoxMesh%division(J)-1)*(P(J)-BoxMesh%extent(J,1))/(BoxMesh%extent(J,2)-BoxMesh%extent(J,1)))+1
 IF (Dir(J).lt.1.or.Dir(J).gt.BoxMesh%division(J)-1) THEN
  WRITE(*,*) 'au',J,dir(j)
  stop
 END IF
END DO

n1=P3-P1
n2=P2-P1

DNX = +n2(2)*n1(1) - n2(3)*n1(2)
DNY = -n2(1)*n1(3) + n2(3)*n1(1)
DNZ = +n2(1)*n1(2) - n2(2)*n1(1)

dArea = 0.5d0*SQRT(DNX*DNX + DNY*DNY + DNZ*DNZ)


TriMesh%d(Dir(1),Dir(2),Dir(3)) = TriMesh%d(Dir(1),Dir(2),Dir(3)) + dArea

! WRITE(*,*) 'saving triangle'

WRITE(2,'(A)') "facet normal 1.0 0.0 0.0"
WRITE(2,'(A)') " outer loop"
WRITE(2,'(A,3ES12.4)') "   vertex",P1
WRITE(2,'(A,3ES12.4)') "   vertex",P2
WRITE(2,'(A,3ES12.4)') "   vertex",P3
WRITE(2,'(A)') " endloop"
WRITE(2,'(A)')  "endfacet"


END SUBROUTINE SaveTriangle
!
! -------------------------------------------------------------------------------
!
SUBROUTINE FindLowerBorder(X,V,iDir)
REAL*8 X,V
integer iDir
integer iIncrement
real*8 dIncrement

iIncrement = INT((BoxMesh%division(iDir)-1)*(X-BoxMesh%extent(iDir,1))/(BoxMesh%extent(iDir,2)-BoxMesh%extent(iDir,1)))
V = BoxMesh%extent(iDir,1) + dble(iIncrement)*(BoxMesh%extent(iDir,2)-BoxMesh%extent(iDir,1))/DBLE(BoxMesh%division(iDir)-1)

if (V.gt.X) V = V - (BoxMesh%extent(iDir,2)-BoxMesh%extent(iDir,1))/DBLE(BoxMesh%division(iDir)-1)

END SUBROUTINE FindLowerBorder
!
! -------------------------------------------------------------------------------
!
SUBROUTINE Output_TriMesh()
implicit none
INTEGER i,j,k
CHARACTER cf*(256)
integer nvt,nel
INTEGER :: iunit=123

nvt = BoxMesh%division(1)*BoxMesh%division(2)*BoxMesh%division(3)
nel = (BoxMesh%division(1)-1)*(BoxMesh%division(2)-1)*(BoxMesh%division(3)-1)

WRITE(cf,'(A)') ADJUSTL(TRIM(cProjectFolder))//'/Coarse_meshDir/'//adjustl(trim(cProjectGridFile))
WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cf))//"'"
OPEN(UNIT=iunit,FILE=ADJUSTL(TRIM(cf)))
WRITE(iunit,*) 'Coarse mesh exported by DeViSoR TRI3D exporter'
WRITE(iunit,*) 'Parametrisierung PARXC, PARYC, TMAXC'
WRITE(iunit,'(2I8,A)') NEL,NVT, " 1 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE"

WRITE(iunit,'(A)') 'DCORVG'
 do i=1,BoxMesh%division(1)
  do j=1,BoxMesh%division(2)
   do k=1,BoxMesh%division(3)
    WRITE(iunit,'(3ES13.5)') TriMesh%x(i),TriMesh%y(j),TriMesh%z(k)
   end do
  end do
 end do

WRITE(iunit,'(A)') 'KVERT'
 do i=1,BoxMesh%division(1)-1
  do j=1,BoxMesh%division(2)-1
   do k=1,BoxMesh%division(3)-1
    write(iunit, '(8I8)') &
     k + (j-1)*BoxMesh%division(3)     + (i-1)*BoxMesh%division(2)*BoxMesh%division(3),&
     k + (j-1)*BoxMesh%division(3) + 1 + (i-1)*BoxMesh%division(2)*BoxMesh%division(3),&
     k + (j)*BoxMesh%division(3)+1     + (i-1)*BoxMesh%division(2)*BoxMesh%division(3),&
     k + (j)*BoxMesh%division(3)       + (i-1)*BoxMesh%division(2)*BoxMesh%division(3),&
     
     k + (j-1)*BoxMesh%division(3)     + (i)*BoxMesh%division(2)*BoxMesh%division(3) ,&
     k + (j-1)*BoxMesh%division(3) + 1 + (i)*BoxMesh%division(2)*BoxMesh%division(3) ,&
     k + (j)*BoxMesh%division(3)+1     + (i)*BoxMesh%division(2)*BoxMesh%division(3) ,&
     k + (j)*BoxMesh%division(3)       + (i)*BoxMesh%division(2)*BoxMesh%division(3) 
   end do
  end do
 end do

WRITE(iunit,'(A)') 'KNPR'
DO i = 1,nvt
  WRITE(iunit,'(I8)') 0
END DO

CLOSE(iunit)

WRITE(cf,'(A)') ADJUSTL(TRIM(cProjectFolder))//'/Coarse_meshDir/'//adjustl(trim(cShortProjectFile))
OPEN(UNIT=iunit,FILE=ADJUSTL(TRIM(cf)))
  WRITE(iunit,'(A)') cProjectGridFile
CLOSE(iunit)

END SUBROUTINE Output_TriMesh
!
! -------------------------------------------------------------------------------
!
SUBROUTINE Output_AreaIntenisty()
INTEGER i,j,k
CHARACTER cf*(256)
integer nvt,nel
INTEGER :: iunit=123

nvt = BoxMesh%division(1)*BoxMesh%division(2)*BoxMesh%division(3)
nel = (BoxMesh%division(1)-1)*(BoxMesh%division(2)-1)*(BoxMesh%division(3)-1)

WRITE(cf,'(A)') ADJUSTL(TRIM(cProjectFolder))//'/'//adjustl(trim(cAreaIntensityFile))
WRITE(*,*) "Outputting Area Intensity into: '"//ADJUSTL(TRIM(cf))//"'"
OPEN(UNIT=iunit,FILE=ADJUSTL(TRIM(cf)))

do i=1,BoxMesh%division(1)-1
 do j=1,BoxMesh%division(2)-1
  do k=1,BoxMesh%division(3)-1
   write(iunit, '(A,E16.7)')"        ",REAL(TriMesh%d(i,j,k))
  end do
 end do
end do

close(iunit)

END SUBROUTINE Output_AreaIntenisty
!
! -------------------------------------------------------------------------------
!
SUBROUTINE Output_VTK()

IMPLICIT NONE
INTEGER ive,ivt,ioffset
INTEGER i,j,k
INTEGER :: iunit=123
CHARACTER*(100) filename
integer nvt,nel


nvt = BoxMesh%division(1)*BoxMesh%division(2)*BoxMesh%division(3)
nel = (BoxMesh%division(1)-1)*(BoxMesh%division(2)-1)*(BoxMesh%division(3)-1)

filename=" "
WRITE(filename(1:),'(A)') "mesh.vtu"

WRITE(*,'(104("="))') 
WRITE(*,*) "Outputting vtk file into ",filename

OPEN (UNIT=iunit,FILE=filename)

write(iunit, *)"<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">"
write(iunit, *)"  <UnstructuredGrid>"
write(iunit, *)"    <Piece NumberOfPoints=""",nvt,""" NumberOfCells=""",nel,""">"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the node field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <PointData>"

 write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","ID",""" format=""ascii"">"
 ivt = 0
 do i=1,BoxMesh%division(1)
  do j=1,BoxMesh%division(2)
   do k=1,BoxMesh%division(3)
    ivt = ivt + 1
    write(iunit, '(A,E16.7)')"        ",REAL(ivt)
   end do
  end do
 end do
 write(iunit, *)"        </DataArray>"

 write(iunit, '(A)')"    </PointData>"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the cell field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <CellData>"

 write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","SurfIntensity",""" format=""ascii"">"
 do i=1,BoxMesh%division(1)-1
  do j=1,BoxMesh%division(2)-1
   do k=1,BoxMesh%division(3)-1
    write(iunit, '(A,E16.7)')"        ",REAL(TriMesh%d(i,j,k))
   end do
  end do
 end do
 write(iunit, *)"        </DataArray>"

write(iunit, '(A)')"    </CellData>"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the mesh data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(iunit, '(A)')"      <Points>"
write(iunit, '(A)')"        <DataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3"" format=""ascii"" RangeMin=""0"" RangeMax=""1.0"">"
 do i=1,BoxMesh%division(1)
  do j=1,BoxMesh%division(2)
   do k=1,BoxMesh%division(3)
    write(iunit, '(A,3E16.7)')"        ",TriMesh%x(i),TriMesh%y(j),TriMesh%z(k)
   end do
  end do
 end do
write(iunit, *)"        </DataArray>"
write(iunit, *)"      </Points>"

write(iunit, *)"      <Cells>"
write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""connectivity"" format=""ascii"" RangeMin=""0"" RangeMax=""",nel-1,""">"
 do i=1,BoxMesh%division(1)-1
  do j=1,BoxMesh%division(2)-1
   do k=1,BoxMesh%division(3)-1
    write(iunit, '(A,8I10)')"        ",&
     -1 + k + (j-1)*BoxMesh%division(3)     + (i-1)*BoxMesh%division(2)*BoxMesh%division(3),&
     -1 + k + (j-1)*BoxMesh%division(3) + 1 + (i-1)*BoxMesh%division(2)*BoxMesh%division(3),&
     -1 + k + (j)*BoxMesh%division(3)+1     + (i-1)*BoxMesh%division(2)*BoxMesh%division(3),&
     -1 + k + (j)*BoxMesh%division(3)       + (i-1)*BoxMesh%division(2)*BoxMesh%division(3),&
      
     -1 + k + (j-1)*BoxMesh%division(3)     + (i)*BoxMesh%division(2)*BoxMesh%division(3) ,&
     -1 + k + (j-1)*BoxMesh%division(3) + 1 + (i)*BoxMesh%division(2)*BoxMesh%division(3) ,&
     -1 + k + (j)*BoxMesh%division(3)+1     + (i)*BoxMesh%division(2)*BoxMesh%division(3) ,&
     -1 + k + (j)*BoxMesh%division(3)       + (i)*BoxMesh%division(2)*BoxMesh%division(3) 
   end do
  end do
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
!------------------------------------------------------------
! ----------------------------------------------
SUBROUTINE GetValueFromFile(cFx,cVx,cKx)
use iniparser

CHARACTER*(256) cFx,cVx,cKx
CHARACTER*(256) ctxt,string
integer iEnd,iloc

OPEN(file=ADJUSTL(TRIM(cFx)),unit=698,action='read')
!write(*,*) ':'//ADJUSTL(TRIM(cFx))//':'

call inip_toupper_replace(cKx)

do 

 read(698,FMT='(A256)',IOSTAT=iEnd) string
 IF (iEnd.EQ.-1) EXIT
 iloc = INDEX(string,"=")
! write(*,*) iloc
 read(string(:iloc-1),*) ctxt
 call inip_toupper_replace(ctxt)
 
!  write(*,*) 'keywords: |'//ADJUSTL(TRIM(ctxt))//'|,|'//ADJUSTL(TRIM(cKx))//"|"
 
 if (ADJUSTL(TRIM(ctxt)).eq.ADJUSTL(TRIM(cKx))) THEN
!  write(*,*) 'keyword found!: ',ADJUSTL(TRIM(ctxt))
  read(string(iloc+1:),'(A256)') cVx
  write(*,'(A)') 'keyword found!: '//ADJUSTL(TRIM(ctxt))//' : '//ADJUSTL(TRIM(cVx))
  GOTO 1 
 END IF

end do

WRITE(*,*) 'Keyword '//ADJUSTL(TRIM(cKx))//' was not found!'
STOP

1 close(698)

END SUBROUTINE GetValueFromFile

END 