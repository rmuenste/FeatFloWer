MODULE MeshProcDef
implicit none

INTEGER :: NumberOfSurfaces=0,nParFiles=0,nStitchFiles=0
CHARACTER, ALLOCATABLE :: cParFile(:)*(200),cStitchFile(:)*(200)
INTEGER,ALLOCATABLE  :: kClassicEdge(:,:)
INTEGER,ALLOCATABLE  :: kvert(:,:),karea(:,:),kedge(:,:),kfaces(:,:),kNeighE(:,:),kFaceNeigh(:,:)
INTEGER,ALLOCATABLE  :: knpr(:),kE1(:,:),kE2(:),kE3(:,:)
REAL*8 ,ALLOCATABLE  :: dcorvg(:,:)
INTEGER nArea,ninArea
integer iix,iiy,iiz,iii
CHARACTER*(100) :: cProjectFolder,cProjectFile
character cfile*(200)

INTEGER,ALLOCATABLE  :: khelp(:,:)

INTEGER nel,nvt,iArea,inArea,net

INTEGER i,j,k,l,m

TYPE tParamList
 integer :: nvt,nel
 REAL*8 , allocatable :: dcoor(:,:)
 INTEGER, allocatable :: kvert(:,:),locInd(:),ID(:),JD(:)
END TYPE tParamList
TYPE(tParamList), allocatable :: myPar(:)

TYPE tCube
 INTEGER :: n=0
 INTEGER, ALLOCATABLE :: L(:)
END TYPE tCube

TYPE tElement
 REAL*8  :: dC(3)
 INTEGER :: iC(3)
END TYPE tElement

TYPE tOctTree
 integer :: nn = 30,nx,ny,nz
 real*8 xmax,xmin,ymax,ymin,zmax,zmin,dx,dy,dz,dL
 TYPE(tCube), ALLOCATABLE :: E(:,:,:)
 TYPE(tElement), allocatable :: M(:)
END TYPE tOctTree
TYPE (tOctTree) OctTree 


CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ReadMEsh()
CHARACTER sCommand*(200)

sCommand = 'rm -fr '//ADJUSTL(TRIM(cProjectFolder))//'/*.par'
CALL system(ADJUSTL(TRIM(sCommand)))
sCommand = 'rm -fr '//ADJUSTL(TRIM(cProjectFolder))//'/*.vtp'
CALL system(ADJUSTL(TRIM(sCommand)))
sCommand = 'rm -fr '//ADJUSTL(TRIM(cProjectFolder))//'/*.prj'
CALL system(ADJUSTL(TRIM(sCommand)))

OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cProjectFolder))//"/mesh.tri")
READ(1,*)
READ(1,*)
READ(1,*) NEL,NVT

ALLOCATE(dcorvg(3,nvt))
ALLOCATE(kvert(8,nel))

READ(1,*)
DO i=1,nvt
 READ(1,*) dcorvg(:,i)
END DO

READ(1,*)
DO i=1,nel
 READ(1,*) kvert(:,i)
END DO

CLOSE(1)
END SUBROUTINE ReadMEsh
!----------------------------------------------------------
SUBROUTINE BuildKedge()
INTEGER ind(2),jnd(2)
INTEGER NeighE(2,12)
REAL*8 dfrac,frac
LOGICAL bFound
DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/

dfrac = 0d0
ALLOCATE(kE1(100,nvt))
ALLOCATE(kE2(nvt))
ALLOCATE(kE3(100,nvt))
ALLOCATE(kClassicEdge(12,nel))
kE1 = 0
kE2 = 0
net = 0

DO i=1,nel
  
 frac = dble(i)/dble(nel)
 if (frac.gt.dfrac) then
  dfrac = dfrac + 0.1d0
  write(*,"(I3.3,A)") INT(1d2*frac), "% is done"
 END IF

 DO j=1,12
  bFound = .FALSE.
  ind = kvert(NeighE(:,j),i)
  CALL SORT(ind,2)
  DO k=1,kE2(ind(1))
   IF (kE1(k,ind(1)).eq.ind(2)) THEN
    bFound = .TRUE.
    GOTO 1
   END IF
  END DO

1 CONTINUE

  IF (.NOT.bFound) THEN
   net = net + 1 
   kedge(:,net) = kvert(NeighE(:,j),i)
   kE2(ind(1)) = kE2(ind(1)) + 1
   kE1(kE2(ind(1)),ind(1)) = ind(2)
  END IF

 END DO
END DO

kE1 = 0
kE2 = 0
kE3 = 0
dfrac = 0d0
ALLOCATE(kedge(2,net))
net = 0

DO i=1,nel
  
 frac = dble(i)/dble(nel)
 if (frac.gt.dfrac) then
  dfrac = dfrac + 0.1d0
  write(*,"(I3.3,A)") INT(1d2*frac), "% is done"
 END IF

 DO j=1,12
  bFound = .FALSE.
  ind = kvert(NeighE(:,j),i)
  CALL SORT(ind,2)
  DO k=1,kE2(ind(1))
   IF (kE1(k,ind(1)).eq.ind(2)) THEN
    bFound = .TRUE.
    kClassicEdge(j,i) = KE3(k,ind(1))
    GOTO 2
   END IF
  END DO

2 CONTINUE

  IF (.NOT.bFound) THEN
   net = net + 1 
   kedge(1,net) = kvert(NeighE(1,j),i)
   kedge(2,net) = kvert(NeighE(2,j),i)
!    write(*,*) net,kedge(:,net)
   kE2(ind(1)) = kE2(ind(1)) + 1
   kE1(kE2(ind(1)),ind(1)) = ind(2)
   kE3(kE2(ind(1)),ind(1)) = net
   kClassicEdge(j,i) = net
  END IF

 END DO
END DO

DEALLOCATE(kE1,kE2,kE3)


DO i=1,nel
!  write(*,'(12I10)') kClassicEdge(:,i)
end do

WRITE(*,*) "number of elments is: ", nel
WRITE(*,*) "number of all unique edges is: ", net

END SUBROUTINE BuildKedge
!----------------------------------------------------------
SUBROUTINE BuildKareaOld()
INTEGER ind(4),jnd(4)
INTEGER NeighE(2,12),NeighA(4,6)
REAL*8 dfrac,frac
LOGICAL bFound
DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
INTEGER jArea
LOGICAL THERE1,THERE2,THERE3


INQUIRE( FILE=ADJUSTL(TRIM(cProjectFolder))//"/kfaces.txt", EXIST=THERE1) 
INQUIRE( FILE=ADJUSTL(TRIM(cProjectFolder))//"/karea.txt",  EXIST=THERE2) 
INQUIRE( FILE=ADJUSTL(TRIM(cProjectFolder))//"/kNeighE.txt",EXIST=THERE3) 

IF (THERE1.and.THERE2.and.THERE3) THEN
 OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cProjectFolder))//"/kfaces.txt")
 READ(1,*) nArea
 ALLOCATE(kfaces(4,nArea))
 DO i=1,nArea
  READ(1,*) kfaces(:,i)
 END DO
 CLOSE(1)

 OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cProjectFolder))//"/karea.txt")
 READ(1,*) 
 ALLOCATE(karea(6,nel))
 DO i=1,nel
  READ(1,*) karea(:,i)
 END DO
 CLOSE(1)

 OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cProjectFolder))//"/kNeighE.txt")
 READ(1,*) 
 ALLOCATE(kNeighE(2,nArea))
 DO i=1,nArea
  READ(1,*) kNeighE(:,i)
 END DO
 
ELSE

 iArea = 0
 inArea = 0
 dfrac = 0d0
 ALLOCATE(khelp(4,6*nel))

 DO i=1,nel
  
  frac = dble(i)/dble(nel)
  if (frac.gt.dfrac) then
   dfrac = dfrac + 0.1d0
   write(*,"(I3.3,A)") INT(1d2*frac), "% is done"
  END IF
  
  DO j=1,6
   bFound = .false.
   ind = kvert(NeighA(:,j),i)
   CALL SORT(ind,4)

   DO k=1,i-1
    DO l=1,6
     jnd = kvert(NeighA(:,l),k)
     CALL SORT(jnd,4)

     IF (ind(1).eq.jnd(1).and.ind(2).eq.jnd(2).and.&
         ind(3).eq.jnd(3).and.ind(4).eq.jnd(4)) then
      inArea = inArea + 1 
      bFound = .true.
      GOTO 1
     END IF
    END DO
   END DO

 1 if (.not.bFound) THEN
    iArea = iArea + 1
    khelp(:,iArea) = kvert(NeighA(:,j),i)
   END IF

  END DO
 END DO

 WRITE(*,*) "number of elments is: ", nel
 WRITE(*,*) "number of all unique faces is: ", iArea
 WRITE(*,*) "number of double faces is: ", inArea
 WRITE(*,*) "number of single faces is: ", iArea-inArea

 nArea = iArea
 ninArea = inArea


 ALLOCATE(kfaces(4,nArea))
 DO iArea=1,nArea
  kfaces(:,iArea) = khelp(:,iArea)
 END DO
 deallocate(khelp)

 ALLOCATE(karea(6,nel))
 ALLOCATE(kNeighE(2,nArea))
 karea = 0
 kNeighE = 0

 iArea = 0
 inArea = 0
 dfrac = 0d0


 DO i=1,nel
  
  frac = dble(i)/dble(nel)
  if (frac.gt.dfrac) then
   dfrac = dfrac + 0.1d0
   write(*,"(I3.3,A)") INT(1d2*frac), "% is done"
  END IF
  
  DO j=1,6
   bFound = .false.
   ind = kvert(NeighA(:,j),i)
   CALL SORT(ind,4)

   DO jArea=1,iArea
    jnd = kfaces(:,jArea)
    CALL SORT(jnd,4)
    IF (ind(1).eq.jnd(1).and.ind(2).eq.jnd(2).and.&
        ind(3).eq.jnd(3).and.ind(4).eq.jnd(4)) then
     bFound = .true.
 !     write(*,*) 'jArea:',jArea
     karea(j,i) = jArea
     kNeighE(2,jArea) = i
     GOTO 2
    END IF
   END DO

 2 if (.not.bFound) THEN
    iArea = iArea + 1
    karea(j,i) = iArea
    kNeighE(1,iArea) = i
   END IF

  END DO
 END DO

 OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cProjectFolder))//"/kfaces.txt")
 WRITE(1,*) nArea
 DO i=1,nArea
  WRITE(1,*) kfaces(:,i)
 END DO
 CLOSE(1)

 OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cProjectFolder))//"/karea.txt")
 WRITE(1,*) nel
 DO i=1,nel
  WRITE(1,*) karea(:,i)
 END DO
 CLOSE(1)

 OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cProjectFolder))//"/kNeighE.txt")
 WRITE(1,*) nArea
 DO i=1,nArea
  WRITE(1,*) kNeighE(:,i)
 END DO
 CLOSE(1)
END IF


END SUBROUTINE BuildKareaOld
!----------------------------------------------------------
SUBROUTINE BuildFaceNeighborhood()
integer iel,iiel,jel,iat,jat,kat,lat,jarea,larea,kk
integer ind(4),jnd(4),knd(4)
integer NeighI(4,6),iCount
integer edgei1(4),edgei2(4),iEE1,iEE2
integer edgej1(4),edgej2(4),iFF1,iFF2
logical :: bShow=.false.

data NeighI/2,3,4,5, 1,3,5,6, 1,2,4,6, 1,3,5,6, 1,2,4,6, 2,3,4,5/


allocate(kFaceNeigh(4,nArea))
kFaceNeigh=0

! DO iel=1,nel
!  write(*,'(I8,A,6I8,A,2I8)') iel,' = ',karea(:,iel)
! END DO
! 
! DO iArea=1,nArea
!  write(*,'(I8,A,2I8)') iArea,' = ',kneighE(:,iArea)
! END DO

iarea = 1
kk = 0
DO iel=1,nel
do iat = 1,6
 if (karea(iat,iel).eq.iarea) then
  if (kNeighE(2,iarea).eq.0) then
   if (bShow) write(*,'(2I0,A$)') kk,iarea,' :'
   ind = kfaces(:,iarea)
! ! !    write(*,*) iel,'ind= ',ind
   DO jat=1,4
    kat = NeighI(jat,iat)
    jarea = karea(kat,iel)
    iiel = iel
99  continue
    jel = kNeighE(2,jarea)
    if (jel.eq.iiel) jel = kNeighE(1,jarea)
    if (jel.ne.0) then
     jnd = kfaces(:,jarea)
     CALL CheckCommon(iEE1)
     iCount = 0
     do lat = 1,6
      larea = karea(lat,jel)
      if (larea.ne.jarea) then
       jnd = kfaces(:,larea)
       CALL CheckCommon(iFF1)
       if (iFF1.eq.2) then
        iCount = iCount + 1
! ! !         write(*,*) jel,'jnd= ',jnd,iCount
        if (kNeighE(2,larea).ne.0) then
!         write(*,'(A,I1,A$)') ' X',iCount,' '
         jarea = larea
         iiel = jel
         goto 99
        else
         kFaceNeigh(jat,iarea) = larea
         if (bShow) write(*,'(A,I1,A$)') ' f',iCount,' '
        end if
       end if
      end if
     end do
     if (iCount.eq.0) then
      if (bShow) write(*,'(A,I1,A$)') '  ',iCount,' '
     end if
    else
     kFaceNeigh(jat,iarea) = jarea
     if (bShow) write(*,'(A,I1,A$)') ' F',0,' '
    end if 
   END DO
   kk = kk + 1 
   if (bShow) write(*,'(A)') ' |'
!    pause
  end if
  iarea=iarea+1
 end if
end do
END DO
iarea = iarea-1

write(*,*) 'k = ',iarea,'kk = ',kk

CONTAINS

SUBROUTINE CheckCommon(if1)
integer ik,il,if1,ix
logical bF

if1 = 0

do ik = 1,4
 bF = .false.
 do il = 1,4
  if (ind(ik).eq.jnd(il)) bF = .true.
 end do
 if (bF) then
  if1 = if1 + 1
 end if
end do

END SUBROUTINE CheckCommon

END SUBROUTINE BuildFaceNeighborhood
!----------------------------------------------------------
SUBROUTINE ReadParametrization()
integer, allocatable :: iVerts(:)
integer in,nn,iFile,iiArea
integer nnn,iiiArea,LenStr,iReason

allocate(iVerts(nvt))

iVerts = 0
DO iFile=1,nStitchFiles
 OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cProjectFolder))//"/STITCH/"//ADJUSTL(TRIM(cStitchFile(iFile))))
 READ(1,*) nn
 READ(1,*) 
 DO i=1,nn
  READ(1,*) in
  iVerts(in) = 1
 END DO
 CLOSE(1) 
END DO

nn = 0
do i=1,nvt
 if (iVerts(i).eq.1) nn = nn + 1
end do

iiArea = 0
do iArea=1,nArea
   IF (iVerts(kfaces(1,iArea)).eq.1.and.iVerts(kfaces(1,iArea)).eq.1.and.&
       iVerts(kfaces(3,iArea)).eq.1.and.iVerts(kfaces(4,iArea)).eq.1) THEN
      iiArea = iiArea + 1
      kNeighE(2,iArea) = -(NumberOfSurfaces + 1)
   END IF
end do
IF (iiArea.gt.0) then
 NumberOfSurfaces = NumberOfSurfaces + 1
 WRITE(*,*) 'Stitched parametrization file:',NumberOfSurfaces, 'with nFaces, nVerts:', iiArea,nn
end if

DO iFile=1,nParFiles

 LenStr = LEN(ADJUSTL(TRIM(cParFile(iFile))))
 IF (LenStr.gt.4) THEN
  cFile = ADJUSTL(TRIM(cParFile(iFile)))
  IF     (cFile(LenStr-3:LenStr).EQ.".csv") THEN
   iVerts = 0
   OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cProjectFolder))//"/PAR/"//ADJUSTL(TRIM(cParFile(iFile))))
   READ(1,*) 
   nn = 0
   DO 
    READ(1,*,IOSTAT=iReason) in
    IF (iReason > 0)  THEN
     WRITE(*,*)'  ... something wrong by reading of the parametrization file: "'//ADJUSTL(TRIM(cParFile(iFile)))//'"'
    ELSE IF (iReason < 0) THEN
     EXIT
!     WRITE(*,*)'  ... end of file reached ... '
    ELSE
     nn = nn + 1
!     WRITE(*,*)'  ... do normal stuff ... '
    END IF
   END DO
   REWIND(1)
   READ(1,*) 
   DO i=1,nn
    READ(1,*) in
    iVerts(in) = 1
   END DO
   CLOSE(1) 
   
   iiArea = 0
   do iArea=1,nArea
      IF (iVerts(kfaces(1,iArea)).eq.1.and.iVerts(kfaces(1,iArea)).eq.1.and.&
          iVerts(kfaces(3,iArea)).eq.1.and.iVerts(kfaces(4,iArea)).eq.1) THEN
         iiArea = iiArea + 1
         kNeighE(2,iArea) = -(NumberOfSurfaces + 1)
      END IF
   end do
   NumberOfSurfaces = NumberOfSurfaces + 1
   WRITE(*,*) 'Read parametrization file:',NumberOfSurfaces, 'with nFaces, nVerts:', iiArea,nn
  
  ELSEIF (cFile(LenStr-3:LenStr).EQ.".par") THEN
   iVerts = 0
   OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cProjectFolder))//"/PAR/"//ADJUSTL(TRIM(cParFile(iFile))))
   READ(1,*) nn
   READ(1,*) 
   DO i=1,nn
    READ(1,*) in
    iVerts(in) = 1
   END DO
   CLOSE(1) 

   iiArea = 0
   do iArea=1,nArea
      IF (iVerts(kfaces(1,iArea)).eq.1.and.iVerts(kfaces(1,iArea)).eq.1.and.&
          iVerts(kfaces(3,iArea)).eq.1.and.iVerts(kfaces(4,iArea)).eq.1) THEN
         iiArea = iiArea + 1
         kNeighE(2,iArea) = -(NumberOfSurfaces + 1)
      END IF
   end do
   NumberOfSurfaces = NumberOfSurfaces + 1
   WRITE(*,*) 'Read parametrization file:',NumberOfSurfaces, 'with nFaces, nVerts:', iiArea,nn
  ELSE
   WRITE(*,*) 'Unknown format of the parametrization file: "'//ADJUSTL(TRIM(cParFile(iFile)))//'"'
  END IF
 
 END IF
   
END DO !iFile

WRITE(*,*) NumberOfSurfaces, ' parametrizations from files have been applied'

END SUBROUTINE ReadParametrization
!----------------------------------------------------------
SUBROUTINE BuildFaceLists(iMinFaces,dCrit)
integer iMinFaces
real*8 dCrit

real*8 dn(3)
integer kk,kkk,jArea

do iArea=1,nArea
 if (kNeighE(2,iArea).eq.0) then
  dn = getNorm(iArea)
  kk = 0
  CALL FillList(iArea,dn)
  if (kk.lt.iMinFaces) then
   do jArea=1,nArea 
    if (kNeighE(2,jArea).eq.-(NumberOfSurfaces+1)) kNeighE(2,jArea) = 0
   end do
  else
   NumberOfSurfaces = NumberOfSurfaces + 1
  end if
 end if
end do

!write(*,*) 'NumberOfSurfaces=',NumberOfSurfaces

! open(1,file='faces.txt')
! do iArea=1,nArea
! ! if (kNeighE(2,iArea).eq.0)  write(1,*) 'problem'
!  if (kNeighE(2,iArea).lt.0) write(1,*) iArea,kNeighE(2,iArea)
! end do
! close(1)

contains
RECURSIVE SUBROUTINE FillList(iat,dni)
real*8 dni(3)
integer ia,iat,jat
real*8 dnj(3),daux

 kNeighE(2,iat) = -(NumberOfSurfaces+1)
 do ia=1,4
  jat = kFaceNeigh(ia,iat)
  if (kNeighE(2,jat).eq.0) then

   dnj = getNorm(jat)
   daux = (dni(1)*dnj(1) + dni(2)*dnj(2) + dni(3)*dnj(3))
   if ((1d0-daux).lt.1d-8) THEN
    daux = 0d0
   else
    daux = 45d0*dabs(dacos(dabs(daux))/DATAN(1d0))
   end if
   if (daux.lt.dCrit) then
!    if (1d0-abs(daux).lt.dCrit) then
    kk = kk + 1
    CALL FillList(jat,dnj)
   end if
  end if
 end do

END  SUBROUTINE FillList

Function getNorm(iA)
integer :: iA
real*8, dimension(3) :: P1,P2,P3,P4,A2,A3,A4,DN1,DN2
real*8 :: DNAR1,DNAR2
real*8, dimension(3) :: getNorm

 P1(:) = dcorvg(:,kfaces(1,IA))
 P2(:) = dcorvg(:,kfaces(2,IA))
 P3(:) = dcorvg(:,kfaces(3,IA))
 P4(:) = dcorvg(:,kfaces(4,IA))
 
 A2 = P2 - P1
 A3 = P3 - P1
 A4 = P4 - P1

 DN1(1) = (A3(2)*A2(3))-(A3(3)*A2(2))
 DN1(2) = (A3(3)*A2(1))-(A3(1)*A2(3))
 DN1(3) = (A3(1)*A2(2))-(A3(2)*A2(1))
 DNAR1=SQRT(DN1(1)**2d0 + DN1(2)**2d0 + DN1(3)**2d0)
 DN1(1) = DN1(1)/DNAR1
 DN1(2) = DN1(2)/DNAR1
 DN1(3) = DN1(3)/DNAR1

 DN2(1) = (A4(2)*A3(3))-(A4(3)*A3(2))
 DN2(2) = (A4(3)*A3(1))-(A4(1)*A3(3))
 DN2(3) = (A4(1)*A3(2))-(A4(2)*A3(1))
 DNAR2=SQRT(DN2(1)**2d0 + DN2(2)**2d0 + DN2(3)**2d0)
 DN2(1) = DN2(1)/DNAR2
 DN2(2) = DN2(2)/DNAR2
 DN2(3) = DN2(3)/DNAR2

!  write(*,*) 'jat',dn1,dn2
 getNorm=0.5d0*(DN1+DN2)

End Function getNorm

END SUBROUTINE BuildFaceLists
!----------------------------------------------------------
SUBROUTINE ExportParFiles()
character cfile*(200)
integer, allocatable :: iVerts(:)
integer nn,nm,ivt,iel,iS

allocate(iVerts(nvt))
allocate(myPar(0:NumberOfSurfaces))

DO iS = 1,NumberOfSurfaces
 iVerts = 0
 nn = 0
 allocate(myPar(iS)%locInd(nvt))
 myPar(iS)%locInd = 0
 do iArea=1,nArea
  if (kNeighE(2,iArea).eq.-iS) then
   nn = nn + 1
   iVerts(kfaces(:,iArea)) = 1
  end if
 end do
 nm = 0
 do ivt=1,nvt
  if (iVerts(ivt).eq.1) nm = nm + 1
 end do

 allocate(myPar(iS)%dcoor(3,nm))
 myPar(iS)%dcoor = 0d0
 allocate(myPar(iS)%kvert(4,nn))
 allocate(myPar(iS)%ID(nm))
 
 myPar(iS)%kvert = 0
 myPar(iS)%ID = 0
 myPar(iS)%nel = nn
 myPar(iS)%nvt = nm
 write(*,*) iS,' : ', myPar(iS)%nel,myPar(iS)%nvt
 
 write(cFile,'(A,I2.2,A)') trim(cProjectFolder)//'/',iS,'.par'
 open(1,file=trim(cFile))
 write(1,*) myPar(iS)%nvt, 'Wall'
 write(1,*) '" "'
 nm = 0
 do ivt=1,nvt
  if (iVerts(ivt).eq.1) THEN
   write(1,*) ivt
   nm = nm + 1
   myPar(iS)%dcoor(:,nm) = dcorvg(:,ivt)
   myPar(iS)%locInd(ivt) = nm
   myPar(iS)%ID(nm) = ivt
  end if
 end do
 
 nn = 0
 do iArea=1,nArea
  if (kNeighE(2,iArea).eq.-iS) then
   nn = nn + 1
   myPar(iS)%kvert(:,nn) = myPar(iS)%locInd(kfaces(:,iArea))
  end if
 end do
 
!  write(1,*) 
!  do ivt=1,myPar(iS)%nvt
!   write(1,'(I8,3ES14.6)') ivt,myPar(iS)%dcoor(:,ivt)
!  end do
! 
!  write(1,*) 
!  do iel=1,myPar(iS)%nel
!   write(1,'(5I8)') iel,myPar(iS)%kvert(:,iel)
!  end do
! 0
!  write(1,*) 
!  do ivt=1,myPar(iS)%nvt
!   write(1,'(I8)') myPar(iS)%ID(ivt)
!  end do
 
 close(1)
 
END DO

 iVerts = 0
 nn = 0
 nm = 0
 DO iS=1,NumberOfSurfaces
  DO ivt =1,myPar(iS)%nvt
   iVerts(myPar(iS)%ID(ivt)) = 1
  END DO
  do ivt=1,nvt
   if (iVerts(ivt).eq.1) nm = nm + 1
  end do
  nn = nn + myPar(iS)%nel
 END DO

 myPar(0)%nel = nn
 myPar(0)%nvt = nm
 allocate(myPar(0)%locInd(nvt))
 myPar(0)%locInd = 0
 allocate(myPar(0)%dcoor(3,nm))
 myPar(0)%dcoor = 0d0
 allocate(myPar(0)%kvert(4,nn))
 allocate(myPar(0)%ID(nm))
 allocate(myPar(0)%JD(nn))

 nm = 0
 do ivt=1,nvt
  if (iVerts(ivt).eq.1) then
   nm = nm + 1
   myPar(0)%dcoor(:,nm) = dcorvg(:,ivt)
   myPar(0)%locInd(ivt) = nm
   myPar(0)%ID(nm) = ivt
  end if
 end do
 
 nn = 0
 do iArea=1,nArea
  if (kNeighE(2,iArea).lt.0) then
   nn = nn + 1
   myPar(0)%kvert(:,nn) = myPar(0)%locInd(kfaces(:,iArea))
   myPar(0)%JD(nn) = -kNeighE(2,iArea)
  end if
 end do

 END SUBROUTINE ExportParFiles
!----------------------------------------------------------
SUBROUTINE SORT(LW,N)
INTEGER LW(N),LWA
INTEGER I,J,N

DO I=2,N
 DO J=N,I,-1
  IF (LW(J).LT.LW(J-1)) THEN
   LWA     = LW(J)
   LW(J)   = LW(J-1)
   LW(J-1) = LWA
  END IF
 END DO
END DO

END SUBROUTINE SORT
!----------------------------------------------------------
SUBROUTINE Output_SingleSurfToVTK
integer iS,ivt,iel,nn
character cfile*(200)
logical bExist

 iS = 0
 write(cFile,'(A,I2.2,A)') trim(cProjectFolder)//'/',iS,'.vtp'
 open(1,file=trim(cFile))
 write(1,'(A)')'<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian" header_type="UInt64">'
 write(1,'(A)')'  <PolyData>'
 write(1,'(A,I0,A,I0,A)')'    <Piece NumberOfPoints="',myPar(iS)%nvt,'" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="',myPar(iS)%nel,'">'
 write(1,'(A)')'      <PointData>'
 write(1,'(A)')'        <DataArray type="Float32" Name="ID" format="ascii" RangeMin="17" RangeMax="2074.25">'
 do ivt=1,myPar(iS)%nvt
  write(1,'(A10,E16.7)')"          ",REAL(myPar(iS)%ID(ivt))
 end do
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'      </PointData>'
 write(1,'(A)')'      <CellData>'
 write(1,'(A)')'        <DataArray type="Float32" Name="JD" format="ascii" RangeMin="17" RangeMax="2074.25">'
 do ivt=1,myPar(iS)%nel
  write(1,'(A10,E16.7)')"          ",REAL(myPar(iS)%JD(ivt))
 end do
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'      </CellData>'
 write(1,'(A)')'      <Points>'
 write(1,'(A)')'        <DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii" RangeMin="18.299993964" RangeMax="72.450446961">'
 do ivt=1,myPar(iS)%nvt
  write(1,'(A10,3E16.7)')"          ",REAL(myPar(iS)%dcoor(1,ivt)),REAL(myPar(iS)%dcoor(2,ivt)),REAL(myPar(iS)%dcoor(3,ivt))
 end do
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'      </Points>'
 write(1,'(A)')'      <Verts>'
 write(1,'(A)')'        <DataArray type="Int64" Name="connectivity" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">'
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'        <DataArray type="Int64" Name="offsets" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">'
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'      </Verts>'
 write(1,'(A)')'      <Lines>'
 write(1,'(A)')'        <DataArray type="Int64" Name="connectivity" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">'
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'        <DataArray type="Int64" Name="offsets" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">'
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'      </Lines>'
 write(1,'(A)')'      <Strips>'
 write(1,'(A)')'        <DataArray type="Int64" Name="connectivity" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">'
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'        <DataArray type="Int64" Name="offsets" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">'
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'      </Strips>'
 write(1,'(A)')'      <Polys>'
 write(1,'(A)')'        <DataArray type="Int64" Name="connectivity" format="ascii" RangeMin="0" RangeMax="1371">'
 do ivt=1,myPar(iS)%nel
  write(1,'(A10,4I8)')"          ",myPar(iS)%kvert(1,ivt)-1,myPar(iS)%kvert(2,ivt)-1,myPar(iS)%kvert(3,ivt)-1,myPar(iS)%kvert(4,ivt)-1
 end do
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'        <DataArray type="Int64" Name="offsets" format="ascii" RangeMin="4" RangeMax="5488">'                                                                     
 nn = 0
 do ivt=1,myPar(iS)%nel
   nn = nn + 4
  write(1,'(A10,I8)')"          ",nn
 end do
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'      </Polys>'
 write(1,'(A)')'    </Piece>'
 write(1,'(A)')'  </PolyData>'
 write(1,'(A)')'</VTKFile>'

#ifdef __GFORTRAN__
 INQUIRE(file=trim(cProjectFolder)//'/PAR',EXIST=bExist)
#elif defined __INTEL_COMPILER
 INQUIRE(directory=trim(cProjectFolder)//'/PAR',EXIST=bExist)
#else
 INQUIRE(file=trim(cProjectFolder)//'/PAR',EXIST=bExist)
#endif

 if (.not.bExist) then
  call system('mkdir '//trim(cProjectFolder)//'/PAR')
 end if
 write(cFile,'(A,I2.2,A)') trim(cProjectFolder)//'/PAR/mesh.tri'
 open(1,file=trim(cFile))
 WRITE(1,*) 'Coarse mesh exported by ParQ2P1 TRI3D exporter'
 WRITE(1,*) 'Parametrisierung PARXC, PARYC, TMAXC'
 WRITE(1,'(2I8,A)') NEL,NVT, " 1 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE"

 write(1,*) 'DCORVG'
 DO i=1,nvt
 WRITE(1,'(3ES16.8)') 1d-1*dcorvg(:,i)
 END DO

 write(1,*) 'KVERT'
 DO i=1,nel
  write(1,'(8I8)') kvert(:,i)
 END DO

 write(1,*) 'KNPR'
 DO i=1,nvt
  write(1,'(I8)') 0
 END DO

 CLOSE(1)
 
END SUBROUTINE Output_SingleSurfToVTK
!----------------------------------------------------------
SUBROUTINE Output_SurfToVTK
integer iS,ivt,iel,nn

DO iS = 1,NumberOfSurfaces
 write(cFile,'(A,I2.2,A)') trim(cProjectFolder)//'/',iS,'.vtp'
 open(1,file=trim(cFile))
 write(1,'(A)')'<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian" header_type="UInt64">'
 write(1,'(A)')'  <PolyData>'
 write(1,'(A,I0,A,I0,A)')'    <Piece NumberOfPoints="',myPar(iS)%nvt,'" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="',myPar(iS)%nel,'">'
 write(1,'(A)')'      <PointData>'
 write(1,'(A)')'        <DataArray type="Float32" Name="ID" format="ascii" RangeMin="17" RangeMax="2074.25">'
 do ivt=1,myPar(iS)%nvt
  write(1,'(A10,E16.7)')"          ",REAL(myPar(iS)%ID(ivt))
 end do
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'      </PointData>'
!  write(1,'(A)')'      <CellData>'
!  write(1,'(A)')'        <DataArray type="Float32" Name="ID" format="ascii" RangeMin="17" RangeMax="2074.25">'
!  write(1,'(A)')'        </DataArray>'
!  write(1,'(A)')'      </CellData>'
 write(1,'(A)')'      <Points>'
 write(1,'(A)')'        <DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii" RangeMin="18.299993964" RangeMax="72.450446961">'
 do ivt=1,myPar(iS)%nvt
  write(1,'(A10,3E16.7)')"          ",REAL(myPar(iS)%dcoor(1,ivt)),REAL(myPar(iS)%dcoor(2,ivt)),REAL(myPar(iS)%dcoor(3,ivt))
 end do
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'      </Points>'
 write(1,'(A)')'      <Verts>'
 write(1,'(A)')'        <DataArray type="Int64" Name="connectivity" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">'
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'        <DataArray type="Int64" Name="offsets" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">'
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'      </Verts>'
 write(1,'(A)')'      <Lines>'
 write(1,'(A)')'        <DataArray type="Int64" Name="connectivity" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">'
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'        <DataArray type="Int64" Name="offsets" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">'
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'      </Lines>'
 write(1,'(A)')'      <Strips>'
 write(1,'(A)')'        <DataArray type="Int64" Name="connectivity" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">'
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'        <DataArray type="Int64" Name="offsets" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">'
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'      </Strips>'
 write(1,'(A)')'      <Polys>'
 write(1,'(A)')'        <DataArray type="Int64" Name="connectivity" format="ascii" RangeMin="0" RangeMax="1371">'
 do ivt=1,myPar(iS)%nel
  write(1,'(A10,4I8)')"          ",myPar(iS)%kvert(1,ivt)-1,myPar(iS)%kvert(2,ivt)-1,myPar(iS)%kvert(3,ivt)-1,myPar(iS)%kvert(4,ivt)-1
 end do
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'        <DataArray type="Int64" Name="offsets" format="ascii" RangeMin="4" RangeMax="5488">'                                                                     
 nn = 0
 do ivt=1,myPar(iS)%nel
   nn = nn + 4
  write(1,'(A10,I8)')"          ",nn
 end do
 write(1,'(A)')'        </DataArray>'
 write(1,'(A)')'      </Polys>'
 write(1,'(A)')'    </Piece>'
 write(1,'(A)')'  </PolyData>'
 write(1,'(A)')'</VTKFile>'
END DO


END SUBROUTINE Output_SurfToVTK
!----------------------------------------------------------
SUBROUTINE Output_VTK

IMPLICIT NONE
INTEGER ive,ivt,ioffset
INTEGER :: iunit=123
CHARACTER*(100) filename

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

!   write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Velocity",""" NumberOfComponents=""3"" format=""ascii"">"
!   do ivt=1,NoOfVert
!    write(iunit, '(A,3E16.7)')"        ",REAL(QuadSc%ValU(ivt)),REAL(QuadSc%ValV(ivt)),REAL(QuadSc%ValW(ivt))
!   end do
!   write(iunit, *)"        </DataArray>"
 

 write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","ID",""" format=""ascii"">"
 do ivt=1,nvt
  write(iunit, '(A,E16.7)')"        ",REAL(ivt)
 end do
 write(iunit, *)"        </DataArray>"

!  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","distance",""" format=""ascii"">"
!  do ivt=1,nvt
!   write(iunit, '(A,E16.7)')"        ",REAL(distance(ivt))
!  end do
!  write(iunit, *)"        </DataArray>"
! 
!  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","monitor",""" format=""ascii"">"
!  do ivt=1,nvt
!   write(iunit, '(A,E16.7)')"        ",REAL(monitor(ivt))
!  end do
!  write(iunit, *)"        </DataArray>"

 write(iunit, '(A)')"    </PointData>"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the cell field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write(iunit, '(A)')"    <CellData>"
! 
!    write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Pressure_E",""" format=""ascii"">"
!    do ivt=1,NoOfElem
!     ive = 4*(ivt-1)+1
!     write(iunit, '(A,E16.7)')"        ",REAL(LinSc%ValP(NLMAX-1)%x(ive))
!    end do
!    write(iunit, *)"        </DataArray>"
! 
! write(iunit, '(A)')"    </CellData>"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the mesh data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(iunit, '(A)')"      <Points>"
write(iunit, '(A)')"        <DataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3"" format=""ascii"" RangeMin=""0"" RangeMax=""1.0"">"
do ivt=1,nvt
 write(iunit,'(A10,3E16.7)')"          ",REAL(dcorvg(1,ivt)),REAL(dcorvg(2,ivt)),REAL(dcorvg(3,ivt))
end do
write(iunit, *)"        </DataArray>"
write(iunit, *)"      </Points>"

write(iunit, *)"      <Cells>"
write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""connectivity"" format=""ascii"" RangeMin=""0"" RangeMax=""",nel-1,""">"
do ive=1,nel   
 write(iunit, '(8I10)')kvert(1,ive)-1,kvert(2,ive)-1,kvert(3,ive)-1,kvert(4,ive)-1,&
                       kvert(5,ive)-1,kvert(6,ive)-1,kvert(7,ive)-1,kvert(8,ive)-1
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
SUBROUTINE BuildOctTree
integer i1
real*8 dmax,dmin,p(3),dSize

ALLOCATE(OctTree%M(nel))

dSize = 0d0
DO i=1,nel
 p = 0d0
 do j=1,8
  p = p + 0.125d0*dcorvg(:,kvert(j,i))
 end do
 OctTree%M(i)%dc = p
 
 dmax = -1e30
 dmin = +1e30
 do j=1,8
  if (dmax.lt.dcorvg(1,kvert(j,i))) dmax = dcorvg(1,kvert(j,i))
  if (dmin.gt.dcorvg(1,kvert(j,i))) dmin = dcorvg(1,kvert(j,i))
 end do
 dSize = max(dSize,dmax-dmin)

 dmax = -1e30
 dmin = +1e30
 do j=1,8
  if (dmax.lt.dcorvg(2,kvert(j,i))) dmax = dcorvg(2,kvert(j,i))
  if (dmin.gt.dcorvg(2,kvert(j,i))) dmin = dcorvg(2,kvert(j,i))
 end do
 dSize = max(dSize,dmax-dmin)
 
 dmax = -1e30
 dmin = +1e30
 do j=1,8
  if (dmax.lt.dcorvg(3,kvert(j,i))) dmax = dcorvg(3,kvert(j,i))
  if (dmin.gt.dcorvg(3,kvert(j,i))) dmin = dcorvg(3,kvert(j,i))
 end do
 dSize = max(dSize,dmax-dmin)
end do

OctTree%xmin = +1d30
OctTree%xmax = -1d30

OctTree%ymin = +1d30
OctTree%ymax = -1d30

OctTree%zmin = +1d30
OctTree%zmax = -1d30

DO i=1,nvt
 if (OctTree%xmin.gt.dcorvg(1,i)) OctTree%xmin=dcorvg(1,i)
 if (OctTree%ymin.gt.dcorvg(2,i)) OctTree%ymin=dcorvg(2,i)
 if (OctTree%zmin.gt.dcorvg(3,i)) OctTree%zmin=dcorvg(3,i)
 
 if (OctTree%xmax.lt.dcorvg(1,i)) OctTree%xmax=dcorvg(1,i)
 if (OctTree%ymax.lt.dcorvg(2,i)) OctTree%ymax=dcorvg(2,i)
 if (OctTree%zmax.lt.dcorvg(3,i)) OctTree%zmax=dcorvg(3,i)
end do

OctTree%dx = OctTree%xmax - OctTree%xmin
OctTree%dy = OctTree%ymax - OctTree%ymin
OctTree%dz = OctTree%zmax - OctTree%zmin

write(*,*) 'OctTree%xmin,OctTree%xmax: ',OctTree%xmin,OctTree%xmax,OctTree%dx
write(*,*) 'OctTree%ymin,OctTree%ymax: ',OctTree%ymin,OctTree%ymax,OctTree%dy
write(*,*) 'OctTree%zmin,OctTree%zmax: ',OctTree%zmin,OctTree%zmax,OctTree%dz

if (OctTree%dx.ge.OctTree%dy.and.OctTree%dx.ge.OctTree%dz) OctTree%dL = OctTree%dx
if (OctTree%dy.ge.OctTree%dx.and.OctTree%dy.ge.OctTree%dz) OctTree%dL = OctTree%dy
if (OctTree%dz.ge.OctTree%dx.and.OctTree%dz.ge.OctTree%dy) OctTree%dL = OctTree%dz

OctTree%nn = MAX(1,NINT(OctTree%dL / (1.1d0*dSize)))
WRITE(*,*) 'max cube size:', dSize, OctTree%nn

OctTree%nx = max(1,nint(DBLE(OctTree%nn)*OctTree%dx/OctTree%dL))
OctTree%ny = max(1,nint(DBLE(OctTree%nn)*OctTree%dy/OctTree%dL))
OctTree%nz = max(1,nint(DBLE(OctTree%nn)*OctTree%dz/OctTree%dL))
write(*,*) 'Limiting size: ',OctTree%dL,OctTree%nx,OctTree%ny,OctTree%nz

ALLOCATE(OctTree%E(0:OctTree%nx+1,0:OctTree%ny+1,0:OctTree%nz+1))

DO i=1,nel

 p = 0d0
 do j=1,8
  p = p + 0.125d0*dcorvg(:,kvert(j,i))
 end do
 OctTree%M(i)%dc = p
 
 do i1=1,OctTree%nx
  dmin = OctTree%xmin + dble(i1-1)*OctTree%dx/dble(OctTree%nx)
  dmax = OctTree%xmin + dble(i1-0)*OctTree%dx/dble(OctTree%nx)
  if (p(1).ge.dmin.and.p(1).le.dmax) then
   iiX = i1
  end if  
 end do

 do i1=1,OctTree%ny
  dmin = OctTree%ymin + dble(i1-1)*OctTree%dy/dble(OctTree%ny)
  dmax = OctTree%ymin + dble(i1-0)*OctTree%dy/dble(OctTree%ny)
  if (p(2).ge.dmin.and.p(2).le.dmax) then
   iiY = i1
  end if  
 end do

 do i1=1,OctTree%nz
  dmin = OctTree%zmin + dble(i1-1)*OctTree%dz/dble(OctTree%nz)
  dmax = OctTree%zmin + dble(i1-0)*OctTree%dz/dble(OctTree%nz)
  if (p(3).ge.dmin.and.p(3).le.dmax) then
   iiZ = i1
  end if  
 end do

 OctTree%E(iix,iiy,iiz)%n = OctTree%E(iix,iiy,iiz)%n + 1
!  write(*,'(3I,3ES12.4)'), iix,iiy,iiz, p 
end do

! pause
do iix=1,OctTree%nx
 do iiy=1,OctTree%ny
  do iiz=1,OctTree%nz
   ALLOCATE(OctTree%E(iix,iiy,iiz)%L(OctTree%E(iix,iiy,iiz)%n))
   OctTree%E(iix,iiy,iiz)%n = 0
  end do
 end do
end do

DO i=1,nel

 p = 0d0
 do j=1,8
  p = p + 0.125d0*dcorvg(:,kvert(j,i))
 end do
 
 do i1=1,OctTree%nx
  dmin = OctTree%xmin + dble(i1-1)*OctTree%dx/dble(OctTree%nx)
  dmax = OctTree%xmin + dble(i1-0)*OctTree%dx/dble(OctTree%nx)
  if (p(1).ge.dmin.and.p(1).le.dmax) then
   iiX = i1
!    exit
  end if  
 end do

 do i1=1,OctTree%ny
  dmin = OctTree%ymin + dble(i1-1)*OctTree%dy/dble(OctTree%ny)
  dmax = OctTree%ymin + dble(i1-0)*OctTree%dy/dble(OctTree%ny)
  if (p(2).ge.dmin.and.p(2).le.dmax) then
   iiY = i1
!    exit
  end if  
 end do

 do i1=1,OctTree%nz
  dmin = OctTree%zmin + dble(i1-1)*OctTree%dz/dble(OctTree%nz)
  dmax = OctTree%zmin + dble(i1-0)*OctTree%dz/dble(OctTree%nz)
  if (p(3).ge.dmin.and.p(3).le.dmax) then
   iiZ = i1
!    exit
  end if  
 end do

 OctTree%E(iix,iiy,iiz)%n = OctTree%E(iix,iiy,iiz)%n + 1
 OctTree%E(iix,iiy,iiz)%L(OctTree%E(iix,iiy,iiz)%n) = i
 
 OctTree%M(i)%ic = [iix,iiy,iiz]
end do

iii = 0
do iix=1,OctTree%nx
 do iiy=1,OctTree%ny
  do iiz=1,OctTree%nz
   iii = iii + OctTree%E(iix,iiy,iiz)%n
!    if (OctTree%E(iix,iiy,iiz)%n.ne.0) write(*,'(4I,A,999I)') iix,iiy,iiz,OctTree%E(iix,iiy,iiz)%n,' : ',OctTree%E(iix,iiy,iiz)%L(:)
  end do
 end do
end do

write(*,*) 'nn = ', iii, nel
END SUBROUTINE BuildOctTree
!----------------------------------------------------------
SUBROUTINE BuildKarea()
INTEGER ind(4),jnd(4)
INTEGER NeighE(2,12),NeighA(4,6)
REAL*8 dfrac,frac
LOGICAL bFound
DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
INTEGER jArea
LOGICAL THERE1,THERE2,THERE3

INQUIRE( FILE=ADJUSTL(TRIM(cProjectFolder))//"/kfaces.txt", EXIST=THERE1) 
INQUIRE( FILE=ADJUSTL(TRIM(cProjectFolder))//"/karea.txt",  EXIST=THERE2) 
INQUIRE( FILE=ADJUSTL(TRIM(cProjectFolder))//"/kNeighE.txt",EXIST=THERE3) 

IF (THERE1.and.THERE2.and.THERE3) THEN
 OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cProjectFolder))//"/kfaces.txt")
 READ(1,*) nArea
 ALLOCATE(kfaces(4,nArea))
 DO i=1,nArea
  READ(1,*) kfaces(:,i)
 END DO
 CLOSE(1)

 OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cProjectFolder))//"/karea.txt")
 READ(1,*) 
 ALLOCATE(karea(6,nel))
 DO i=1,nel
  READ(1,*) karea(:,i)
 END DO
 CLOSE(1)

 OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cProjectFolder))//"/kNeighE.txt")
 READ(1,*) 
 ALLOCATE(kNeighE(2,nArea))
 DO i=1,nArea
  READ(1,*) kNeighE(:,i)
 END DO
 
ELSE

 iArea = 0
 inArea = 0
 dfrac = 0d0
 ALLOCATE(khelp(4,6*nel))

 DO i=1,nel
  
  frac = dble(i)/dble(nel)
  if (frac.gt.dfrac) then
   dfrac = dfrac + 0.1d0
   write(*,"(I3.3,A)") INT(1d2*frac), "% is done"
  END IF
  
  DO j=1,6
   bFound = .false.
   ind = kvert(NeighA(:,j),i)
   CALL SORT(ind,4)

   !!!  Find the pair !!!
   do iix=OctTree%M(i)%iC(1)-1,OctTree%M(i)%iC(1)+1
    do iiy=OctTree%M(i)%iC(2)-1,OctTree%M(i)%iC(2)+1
     do iiz=OctTree%M(i)%iC(3)-1,OctTree%M(i)%iC(3)+1
      do iii = 1,OctTree%E(iix,iiy,iiz)%N
       k = OctTree%E(iix,iiy,iiz)%L(iii)
       if (k.lt.i) then
	DO l=1,6
	 jnd = kvert(NeighA(:,l),k)
	 CALL SORT(jnd,4)

	 IF (ind(1).eq.jnd(1).and.ind(2).eq.jnd(2).and.&
	     ind(3).eq.jnd(3).and.ind(4).eq.jnd(4)) then
	  inArea = inArea + 1 
	  bFound = .true.
	  GOTO 1
	 END IF
	END DO
       end if
      end do
     end do
    end do
   end do
   
 1 if (.not.bFound) THEN
    iArea = iArea + 1
    khelp(:,iArea) = kvert(NeighA(:,j),i)
   END IF

  END DO
 END DO

 WRITE(*,*) "number of elments is: ", nel
 WRITE(*,*) "number of all unique faces is: ", iArea
 WRITE(*,*) "number of double faces is: ", inArea
 WRITE(*,*) "number of single faces is: ", iArea-inArea

 nArea = iArea
 ninArea = inArea

 ALLOCATE(kfaces(4,nArea))
 DO iArea=1,nArea
  kfaces(:,iArea) = khelp(:,iArea)
 END DO
 deallocate(khelp)

 ALLOCATE(karea(6,nel))
 ALLOCATE(kNeighE(2,nArea))
 karea = 0
 kNeighE = 0

 iArea = 0
 inArea = 0
 dfrac = 0d0


 DO i=1,nel
  
  frac = dble(i)/dble(nel)
  if (frac.gt.dfrac) then
   dfrac = dfrac + 0.1d0
   write(*,"(I3.3,A)") INT(1d2*frac), "% is done"
  END IF
  
  DO j=1,6
   bFound = .false.
   ind = kvert(NeighA(:,j),i)
   CALL SORT(ind,4)

   !!!  Find the pair !!!
   do iix=OctTree%M(i)%iC(1)-1,OctTree%M(i)%iC(1)+1
    do iiy=OctTree%M(i)%iC(2)-1,OctTree%M(i)%iC(2)+1
     do iiz=OctTree%M(i)%iC(3)-1,OctTree%M(i)%iC(3)+1
      do iii = 1,OctTree%E(iix,iiy,iiz)%N
       k = OctTree%E(iix,iiy,iiz)%L(iii)
       if (k.lt.i) then
	DO l=1,6
	 jnd = kvert(NeighA(:,l),k)
	 CALL SORT(jnd,4)

	 IF (ind(1).eq.jnd(1).and.ind(2).eq.jnd(2).and.&
	     ind(3).eq.jnd(3).and.ind(4).eq.jnd(4)) then
	  bFound = .true.
	  jArea = karea(l,k)
          karea(j,i) = jArea
          kNeighE(2,jArea) = i
	  GOTO 2
	 END IF
	END DO
       end if
      end do
     end do
    end do
   end do
   
 2 if (.not.bFound) THEN
    iArea = iArea + 1
    karea(j,i) = iArea
    kNeighE(1,iArea) = i
   END IF

  END DO
 END DO

 OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cProjectFolder))//"/kfaces.txt")
 WRITE(1,*) nArea
 DO i=1,nArea
  WRITE(1,*) kfaces(:,i)
 END DO
 CLOSE(1)

 OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cProjectFolder))//"/karea.txt")
 WRITE(1,*) nel
 DO i=1,nel
  WRITE(1,*) karea(:,i)
 END DO
 CLOSE(1)

 OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cProjectFolder))//"/kNeighE.txt")
 WRITE(1,*) nArea
 DO i=1,nArea
  WRITE(1,*) kNeighE(:,i)
 END DO
 CLOSE(1)
END IF
 
END SUBROUTINE BuildKarea
!----------------------------------------------------------
END MODULE MeshProcDef
