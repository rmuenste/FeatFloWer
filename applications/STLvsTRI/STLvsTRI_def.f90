MODULE mSTLvsTRI
USE PP3D_MPI, ONLY:myid

implicit none

type tOFFmesh
 REAL*8,allocatable :: coor(:,:)
 integer, allocatable :: kvert(:,:)
 integer nel,nvt
end type tOFFmesh
type (tOFFmesh) :: OFFmesh 
CHARACTER*(200) :: cProjectFolder,cShortProjectFile,cProjectGridFile,cAreaIntensityFile
CHARACTER*(256) :: cOFFMeshFile
integer iOK,iTriang,iRecLevel
TYPE tBoxMesh
 REAL*8 extent(3,2)
 REAL*8 dsize(3)
 INTEGER division(3)
END TYPE tBoxMesh
TYPE(tBoxMesh) BoxMesh

TYPE tTriMesh
 integer nx,ny,nz
 REAL*8, allocatable :: x(:),y(:),z(:)
 REAL*8, allocatable :: d(:,:,:) 
 REAL*8, allocatable :: I(:,:,:) 
END TYPE tTriMesh
TYPE(tTriMesh) TriMesh,CoarseTriMesh

Type triplet
 REAl*8 Q(3)
end type triplet

REAL*8 DomainVolume,VoxelVolume,UnityArea
INTEGER NumberOfElements
INTEGER :: NumberOfElementsInit=1500
REAL*8 :: dSurfIntCrit = 2.5d0
!REAL*8 :: dSurfIntCrit = 6.0d0

 CONTAINS
 
SUBROUTINE readOFFMesh(cFF,cF)
CHARACTER*(*) cF,cFF
integer i,iaux
integer num,k,lt,i0,i1,invt,inel,iNum
character*256 t
character*256, allocatable ::  cOFFfiles(:)

t = '.off'

lt = len(adjustl(trim(t))) - 1
k = 1
num = 0
do
!   print *,cF(k:)
   i = index(cF(k:),adjustl(trim(t)))
   if (i==0) exit
   num = num + 1
   k = k + i + lt + 1
end do
 
write(*,*) "'"//trim(adjustl(cF))//"'", num

allocate(cOFFfiles(num))

k = 1
num = 0
do
   i = index(cF(k:),adjustl(trim(t)))
   if (i==0) exit
   i0 = k
   i1 = (k-1) + i + lt
   num = num + 1
   write(*,*) i0,i1
   read(cF(i0:i1),'(A)') cOFFfiles(num)
   write(*,*) "'"//trim(adjustl(cOFFfiles(num)))//"'"
   k = k + i + lt + 1
end do

OFFmesh%nvt = 0
OFFmesh%nel = 0

DO iNum = 1 , Num
 OPEN(File=trim(adjustl(cFF))//'/'//trim(adjustl(cOFFfiles(num))),unit=1)
 read(1,*)
 read(1,*) invt,inel
 
 OFFmesh%nvt = OFFmesh%nvt + invt
 OFFmesh%nel = OFFmesh%nel + inel
 
 close(1)

END DO

allocate(OFFmesh%coor(OFFmesh%nvt,3))
allocate(OFFmesh%kvert(OFFmesh%nel,3))

OFFmesh%nvt = 0
OFFmesh%nel = 0

DO iNum = 1 , Num

 OPEN(File=trim(adjustl(cFF))//'/'//trim(adjustl(cOFFfiles(num))),unit=1)
 read(1,*)
 read(1,*) invt,inel

 DO i=1,invt
  read(1,*) OFFmesh%coor(OFFmesh%nvt+i,:)
 END DO

 DO i=1,inel
  read(1,*) iaux,OFFmesh%kvert(OFFmesh%nel+i,:)
 END DO

 OFFmesh%nvt = OFFmesh%nvt + invt
 OFFmesh%nel = OFFmesh%nel + inel
 
 close(1)
end do

OFFmesh%kvert =  OFFmesh%kvert + 1

END SUBROUTINE readOFFMesh
! -------------------------------------------------------------------------------
SUBROUTINE DeallocateStructures

deallocate(TriMesh%x)
deallocate(TriMesh%y)
deallocate(TriMesh%z)
deallocate(TriMesh%d)
deallocate(TriMesh%I)

if (allocated(CoarseTriMesh%x)) deallocate(CoarseTriMesh%x)
if (allocated(CoarseTriMesh%y)) deallocate(CoarseTriMesh%y)
if (allocated(CoarseTriMesh%z)) deallocate(CoarseTriMesh%z)
if (allocated(CoarseTriMesh%d)) deallocate(CoarseTriMesh%d)
if (allocated(CoarseTriMesh%I)) deallocate(CoarseTriMesh%I)

END SUBROUTINE DeallocateStructures
!
! -------------------------------------------------------------------------------
!
SUBROUTINE CheckForIntersection
REAL*8  P1(3),P2(3),P3(3)
integer i,iDir

allocate(TriMesh%x(BoxMesh%Division(1)))
allocate(TriMesh%y(BoxMesh%Division(2)))
allocate(TriMesh%z(BoxMesh%Division(3)))
allocate(TriMesh%d(BoxMesh%Division(1)-1,BoxMesh%Division(2)-1,BoxMesh%Division(3)-1))
allocate(TriMesh%I(BoxMesh%Division(1)-1,BoxMesh%Division(2)-1,BoxMesh%Division(3)-1))
TriMesh%x = 0d0
TriMesh%y = 0d0
TriMesh%z = 0d0
TriMesh%d = 0d0
TriMesh%I = 0d0
DO i=1,BoxMesh%division(1)
 TriMesh%x(i) = BoxMesh%Extent(1,1) + dble(i-1)*(BoxMesh%Extent(1,2)-BoxMesh%Extent(1,1))/(BoxMesh%division(1)-1)
END DO
DO i=1,BoxMesh%division(2)
 TriMesh%y(i) = BoxMesh%Extent(2,1) + dble(i-1)*(BoxMesh%Extent(2,2)-BoxMesh%Extent(2,1))/(BoxMesh%division(2)-1)
END DO
DO i=1,BoxMesh%division(3)
 TriMesh%z(i) = BoxMesh%Extent(3,1) + dble(i-1)*(BoxMesh%Extent(3,2)-BoxMesh%Extent(3,1))/(BoxMesh%division(3)-1)
END DO


OPEN(FILE=ADJUSTL(TRIM(cProjectFolder))//'/EXP.stl',unit=2)
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
  RETURN
 END IF
END DO

n1=P3-P1
n2=P2-P1

DNX = +n2(2)*n1(1) - n2(1)*n1(2)
DNY = -n2(1)*n1(3) + n2(3)*n1(1)
DNZ = +n2(3)*n1(2) - n2(2)*n1(3)

dArea = 0.5d0*SQRT(DNX*DNX + DNY*DNY + DNZ*DNZ)


TriMesh%d(Dir(1),Dir(2),Dir(3)) = TriMesh%d(Dir(1),Dir(2),Dir(3)) + ABS(dArea)

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

nvt = CoarseTriMesh%nx*CoarseTriMesh%ny*CoarseTriMesh%nz
nel = (CoarseTriMesh%nx-1)*(CoarseTriMesh%ny-1)*(CoarseTriMesh%nz-1)
! nvt = BoxMesh%division(1)*BoxMesh%division(2)*BoxMesh%division(3)
! nel = (BoxMesh%division(1)-1)*(BoxMesh%division(2)-1)*(BoxMesh%division(3)-1)

WRITE(cf,'(A)') ADJUSTL(TRIM(cProjectFolder))//'/Coarse_meshDir'
myid = 0
CALL CheckIfFolderIsThereCreateIfNot(cf,-1)
myid = 0


WRITE(cf,'(A)') ADJUSTL(TRIM(cProjectFolder))//'/Coarse_meshDir/'//adjustl(trim(cProjectGridFile))
WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cf))//"'"
OPEN(UNIT=iunit,FILE=ADJUSTL(TRIM(cf)))
WRITE(iunit,*) 'Coarse mesh exported by DeViSoR TRI3D exporter'
WRITE(iunit,*) 'Parametrisierung PARXC, PARYC, TMAXC'
WRITE(iunit,'(2I8,A)') NEL,NVT, " 1 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE"

WRITE(iunit,'(A)') 'DCORVG'
 do i=1,CoarseTriMesh%nx
  do j=1,CoarseTriMesh%ny
   do k=1,CoarseTriMesh%nz
    WRITE(iunit,'(3ES13.5)') CoarseTriMesh%x(i),CoarseTriMesh%y(j),CoarseTriMesh%z(k)
   end do
  end do
 end do

WRITE(iunit,'(A)') 'KVERT'
 do i=1,CoarseTriMesh%nx-1
  do j=1,CoarseTriMesh%ny-1
   do k=1,CoarseTriMesh%nz-1
    write(iunit, '(8I8)') &
     k + (j-1)*CoarseTriMesh%nz    + (i-1)*CoarseTriMesh%ny*CoarseTriMesh%nz,&
     k + (j-1)*CoarseTriMesh%nz+ 1 + (i-1)*CoarseTriMesh%ny*CoarseTriMesh%nz,&
     k + (j)*CoarseTriMesh%nz+1     + (i-1)*CoarseTriMesh%ny*CoarseTriMesh%nz,&
     k + (j)*CoarseTriMesh%nz      + (i-1)*CoarseTriMesh%ny*CoarseTriMesh%nz,&
     
     k + (j-1)*CoarseTriMesh%nz    + (i)*CoarseTriMesh%ny*CoarseTriMesh%nz,&
     k + (j-1)*CoarseTriMesh%nz+ 1 + (i)*CoarseTriMesh%ny*CoarseTriMesh%nz,&
     k + (j)*CoarseTriMesh%nz+1     + (i)*CoarseTriMesh%ny*CoarseTriMesh%nz,&
     k + (j)*CoarseTriMesh%nz      + (i)*CoarseTriMesh%ny*CoarseTriMesh%nz
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

do i=1,CoarseTriMesh%nx-1
 do j=1,CoarseTriMesh%ny-1
  do k=1,CoarseTriMesh%nz-1
   
! do i=1,BoxMesh%division(1)-1
!  do j=1,BoxMesh%division(2)-1
!   do k=1,BoxMesh%division(3)-1
   write(iunit, '(A,E16.7)')"        ",REAL(CoarseTriMesh%I(i,j,k))
!    write(iunit, '(A,E16.7)')"        ",REAL(TriMesh%I(i,j,k))
  end do
 end do
end do

close(iunit)

END SUBROUTINE Output_AreaIntenisty
!
! -------------------------------------------------------------------------------
!
SUBROUTINE Output_CoarseMeshVTK()

IMPLICIT NONE
INTEGER ive,ivt,ioffset
INTEGER i,j,k
INTEGER :: iunit=123
CHARACTER*(100) filename
integer nvt,nel

nvt = CoarseTriMesh%nx*CoarseTriMesh%ny*CoarseTriMesh%nz
nel = (CoarseTriMesh%nx-1)*(CoarseTriMesh%ny-1)*(CoarseTriMesh%nz-1)

filename=" "
WRITE(filename(1:),'(A)') ADJUSTL(TRIM(cProjectFolder))//"/CoarseTriMesh.vtu"

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
 do i=1,CoarseTriMesh%nx
  do j=1,CoarseTriMesh%ny
   do k=1,CoarseTriMesh%nz
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
 
 do i=1,CoarseTriMesh%nx-1
  do j=1,CoarseTriMesh%ny-1
   do k=1,CoarseTriMesh%nz-1
    write(iunit, '(A,E16.7)')"        ",REAL(CoarseTriMesh%I(i,j,k))
   end do
  end do
 end do
 write(iunit, *)"        </DataArray>"
 
write(iunit, '(A)')"    </CellData>"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the mesh data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(iunit, '(A)')"      <Points>"
write(iunit, '(A)')"        <DataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3"" format=""ascii"" RangeMin=""0"" RangeMax=""1.0"">"
 do i=1,CoarseTriMesh%nx
  do j=1,CoarseTriMesh%ny
   do k=1,CoarseTriMesh%nz
    write(iunit, '(A,3E16.7)')"        ",CoarseTriMesh%x(i),CoarseTriMesh%y(j),CoarseTriMesh%z(k)
   end do
  end do
 end do
write(iunit, *)"        </DataArray>"
write(iunit, *)"      </Points>"

write(iunit, *)"      <Cells>"
write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""connectivity"" format=""ascii"" RangeMin=""0"" RangeMax=""",nel-1,""">"
 do i=1,CoarseTriMesh%nx-1
  do j=1,CoarseTriMesh%ny-1
   do k=1,CoarseTriMesh%nz-1
    write(iunit, '(A,8I10)')"        ",&
     -1 + k + (j-1)*CoarseTriMesh%nz     + (i-1)*CoarseTriMesh%ny*CoarseTriMesh%nz,&
     -1 + k + (j-1)*CoarseTriMesh%nz + 1 + (i-1)*CoarseTriMesh%ny*CoarseTriMesh%nz,&
     -1 + k + (j)*CoarseTriMesh%nz+1     + (i-1)*CoarseTriMesh%ny*CoarseTriMesh%nz,&
     -1 + k + (j)*CoarseTriMesh%nz       + (i-1)*CoarseTriMesh%ny*CoarseTriMesh%nz,&
      
     -1 + k + (j-1)*CoarseTriMesh%nz     + (i)*CoarseTriMesh%ny*CoarseTriMesh%nz ,&
     -1 + k + (j-1)*CoarseTriMesh%nz + 1 + (i)*CoarseTriMesh%ny*CoarseTriMesh%nz ,&
     -1 + k + (j)*CoarseTriMesh%nz+1     + (i)*CoarseTriMesh%ny*CoarseTriMesh%nz ,&
     -1 + k + (j)*CoarseTriMesh%nz       + (i)*CoarseTriMesh%ny*CoarseTriMesh%nz 
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

END SUBROUTINE Output_CoarseMeshVTK
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
WRITE(filename(1:),'(A)') ADJUSTL(TRIM(cProjectFolder))//"/mesh.vtu"

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

 write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Area",""" format=""ascii"">"
 do i=1,BoxMesh%division(1)-1
  do j=1,BoxMesh%division(2)-1
   do k=1,BoxMesh%division(3)-1
    write(iunit, '(A,E16.7)')"        ",REAL(TriMesh%d(i,j,k))
   end do
  end do
 end do
 write(iunit, *)"        </DataArray>"

 write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","SurfIntensity",""" format=""ascii"">"
 
 do i=1,BoxMesh%division(1)-1
  do j=1,BoxMesh%division(2)-1
   do k=1,BoxMesh%division(3)-1
    write(iunit, '(A,E16.7)')"        ",REAL(TriMesh%d(i,j,k)/UnityArea)
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
SUBROUTINE CoarseSurfIntensity()
integer i,j,k
integer i0,j0,k0,i1,j1,k1,i2,j2,k2
integer ii,jj,kk
REAL*8 daux
REAL*8 dCrit(3)

 CoarseTriMesh%nx = INT(0.5d0*(BoxMesh%Division(1)-1)+1)
 CoarseTriMesh%ny = INT(0.5d0*(BoxMesh%Division(2)-1)+1)
 CoarseTriMesh%nz = INT(0.5d0*(BoxMesh%Division(3)-1)+1)
 
WRITE(*,*) BoxMesh%Division(1:3)
WRITE(*,*) CoarseTriMesh%nx,CoarseTriMesh%nY,CoarseTriMesh%nZ

allocate(CoarseTriMesh%x(CoarseTriMesh%nx))
allocate(CoarseTriMesh%y(CoarseTriMesh%ny))
allocate(CoarseTriMesh%z(CoarseTriMesh%nz))
allocate(CoarseTriMesh%I(CoarseTriMesh%nx,CoarseTriMesh%ny-1,CoarseTriMesh%nz-1))

 CoarseTriMesh%x = 0d0
 CoarseTriMesh%y = 0d0
 CoarseTriMesh%z = 0d0
 CoarseTriMesh%I = 0d0
DO i=1,CoarseTriMesh%nx 
 CoarseTriMesh%x(i) = BoxMesh%Extent(1,1) + dble(i-1)*(BoxMesh%Extent(1,2)-BoxMesh%Extent(1,1))/(CoarseTriMesh%nx -1)
END DO
DO i=1,CoarseTriMesh%ny
 CoarseTriMesh%y(i) = BoxMesh%Extent(2,1) + dble(i-1)*(BoxMesh%Extent(2,2)-BoxMesh%Extent(2,1))/(CoarseTriMesh%ny -1)
END DO
DO i=1,CoarseTriMesh%nz
 CoarseTriMesh%z(i) = BoxMesh%Extent(3,1) + dble(i-1)*(BoxMesh%Extent(3,2)-BoxMesh%Extent(3,1))/(CoarseTriMesh%nz -1)
END DO


do i=1,BoxMesh%division(1)-2,2
 do j=1,BoxMesh%division(2)-2,2
  do k=1,BoxMesh%division(3)-2,2
  
  i0 = (i+1)/2
  j0 = (j+1)/2
  k0 = (k+1)/2
  
!    daux  = max(TriMesh%I(i,j,k  ),TriMesh%I(i+1,j,k  ),TriMesh%I(i,j+1,k  ),TriMesh%I(i+1,j+1,k  ),&
!                TriMesh%I(i,j,k+1),TriMesh%I(i+1,j,k+1),TriMesh%I(i,j+1,k+1),TriMesh%I(i+1,j+1,k+1))
!    daux = 0.125d0*daux

   daux = 0d0
   
   DO i2=-1,1
    DO j2=-1,1
     DO k2=-1,1
     
       dCrit = 0d0
       
       DO i1=0,1
        DO j1=0,1
         DO k1=0,1
          II =  i + i1 + i2
          JJ =  j + j1 + j2
          KK =  k + k1 + k2
          IF ((II.ge.1.and.II.le.BoxMesh%division(1)-1).and.&
              (JJ.ge.1.and.JJ.le.BoxMesh%division(2)-1).and.&
              (KK.ge.1.and.KK.le.BoxMesh%division(3)-1)) THEN
               
              dCrit(1) = dCrit(1) + TriMesh%I(II,JJ,KK)
              dCrit(2) = dCrit(2) + 1d0
              dCrit(3) = max(dCrit(3),TriMesh%I(II,JJ,KK))
          END IF
         END DO
        END DO
       END DO
       
!        if (i0.eq.9.and.j0.eq.5.and.k0.eq.18) then
!         write(*,*) daux,dCrit(1)/dCrit(2),dCrit(3)
!        end if
       if (dCrit(2).eq.8d0) then
        daux = max(daux,2d0*dCrit(1)/dCrit(2))
       end if
   
     END DO
    END DO
   END DO
    
   CoarseTriMesh%I(i0,j0,k0)  = daux
  end do
 end do
end do


! do i=1,BoxMesh%division(1)-2,2
!  do j=1,BoxMesh%division(2)-2,2
!   do k=1,BoxMesh%division(3)-2,2
! 
!    i1 = (i+1)/2
!    j1 = (j+1)/2
!    k1 = (k+1)/2
!    
!    DO ii=i-1,i+1,2
!     IF (ii.ge.1.and.ii.le.BoxMesh%division(1).and.CoarseTriMesh%I(i1,j1,k1).gt.0d0) then
!      daux  = (TriMesh%I(i,j,k) + TriMesh%I(ii,j,k))/2d0
!      CoarseTriMesh%I(i1,j1,k1)  = MAX(CoarseTriMesh%I(i1,j1,k1),daux)
!     END IF
!    END DO
!    
!    DO jj=j-1,j+1,2
!     IF (jj.ge.1.and.j.le.BoxMesh%division(2).and.CoarseTriMesh%I(i1,j1,k1).gt.0d0) then
!      daux  = (TriMesh%I(i,j,k) + TriMesh%I(i,jj,k))/2d0
!      CoarseTriMesh%I(i1,j1,k1)  = MAX(CoarseTriMesh%I(i1,j1,k1),daux)
!     END IF
!    END DO
! ! 
!    DO kk=k-1,k+1,2
!     IF (kk.ge.1.and.k.le.BoxMesh%division(3).and.CoarseTriMesh%I(i1,j1,k1).gt.0d0) then
!      daux  = (TriMesh%I(i,j,k) + TriMesh%I(i,j,kk))/2d0
!      CoarseTriMesh%I(i1,j1,k1)  = MAX(CoarseTriMesh%I(i1,j1,k1),daux)
!     END IF
!    END DO
! 
!   end do
!  end do
! end do
 
END SUBROUTINE CoarseSurfIntensity
!------------------------------------------------------------
SUBROUTINE ConstructSurfIntensity()
integer i,j,k

do i=1,BoxMesh%division(1)-1
 do j=1,BoxMesh%division(2)-1
  do k=1,BoxMesh%division(3)-1
   TriMesh%I(i,j,k) = TriMesh%d(i,j,k)/UnityArea
  end do
 end do
end do
 
END SUBROUTINE ConstructSurfIntensity
!------------------------------------------------------------
SUBROUTINE QualityCheck(dMaxAreaIntensity,nCrit)
integer nCrit
real*8 dMaxAreaIntensity
integer i,j,k

nCrit = 0
dMaxAreaIntensity = 0d0

do i=1,BoxMesh%division(1)-1
 do j=1,BoxMesh%division(2)-1
  do k=1,BoxMesh%division(3)-1
   dMaxAreaIntensity = MAX(dMaxAreaIntensity,TriMesh%I(i,j,k))
   IF (TriMesh%d(i,j,k)/UnityArea.gt.dSurfIntCrit) nCrit = nCrit + 1
  end do
 end do
end do
 
END SUBROUTINE QualityCheck
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