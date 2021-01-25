MODULE MeshRefDef
USE types, ONLY: tMultiMesh,tMesh
USE Parametrization, ONLY: InitParametrization_STRCT,ParametrizeBndryPoints_STRCT,&
    DeterminePointParametrization_STRCT,ProlongateParametrization_STRCT,myParBndr,nBnds
USE var_QuadScalar
USE PP3D_MPI, ONLY:myid,master
USE MESH_Structures
use iniparser
use meshrefvar

implicit none


CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CutMesh()
integer iel,i
real*8 dc(3),p(3),dist1,dS,dRad
integer :: isin, ipc = 0

allocate(MergedMeshDist(nUniqueElems))

do iel=1,nUniqueElems

 dc = 0d0
 do i=1,8
  p = MergedMeshCoor(:,MergedMeshElem(i,iel))
  dc = dc + 0.125d0*p
 end do
 
 dRad = 1d8
 do i=1,8
  dist1 = sqrt((dc(1)-p(1))**2d0 + (dc(2)-p(2))**2d0 + (dc(3)-p(3))**2d0)
  if (dist1.lt.dRad) dRad = dist1
 end do

!  isin = 0
!  call isinelementid(dc(1),dc(2),dc(3),0,isin)
! 
!  if(isin .gt. 0)then
!    dS = -1d0
!  else
!    dS = +1d0
!  end if
!  
!  call getdistanceid(dc(1),dc(2),dc(3),dist1,ipc)
!  dist1 = dist1*dS
 dist1 = 1d0
 
 if (dist1.lt.dRad) THEN
  dist1=+1d0
 else
  dist1=-1d0
 end if
 
 MergedMeshDist(iel)=dist1
 
end do


END SUBROUTINE CutMesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE InitMarker()
logical bM
integer i,j,k,l,m,ivt
real*8 myRand 

ilev = mg_Mesh%nlmax

nel = mg_mesh%level(ilev)%nel
nvt = mg_mesh%level(ilev)%nvt

allocate(markerE(nel))
allocate(markerV(nvt))

markerE = 0
markerV = 0

! j=0
! DO i=1,nel
!  CALL RANDOM_NUMBER(myRand)
!  if (myRand.gt.0.95d0) then
!   j = j +1
!   markerE(i) = 1
!  end if
! end do
IF (initfield.eq.0) then
 CALL Initfield0(markerE,mg_mesh%level(ilev)%kvert,mg_mesh%level(ilev)%dcorvg,nel)
end if

IF (initfield.eq.1) then
 CALL Initfield1(markerE,mg_mesh%level(ilev)%kvert,mg_mesh%level(ilev)%dcorvg,nel,RefinementThickness)
end if

IF (initfield.eq.2) then
 CALL Initfield2(markerE,mg_mesh%level(ilev)%kvert,mg_mesh%level(ilev)%dcorvg,nel,RefinementThickness)
end if

END SUBROUTINE InitMarker
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SetUpMarker()
logical bM
integer i,j,k,l,m,ivt
integer istep

ilev = mg_Mesh%nlmax

nel = mg_mesh%level(ilev)%nel
nvt = mg_mesh%level(ilev)%nvt

!!! First the Marker Corner poionts of the 1st level have to be created
DO i=1,nel
 if (markerE(i).Eq.1) then
  do k=1,8
   j = mg_mesh%level(ilev)%kvert(k,i)
   markerV(j) = 1
  end do
 end if
end do

!!! Construction of the 2nd level Markers
DO i=1,nel
 if (markerE(i).Eq.1) then
  do j=1,8
   ivt = mg_mesh%level(ilev)%kvert(j,i)
   do k=1,mg_mesh%level(ilev)%nvel
    l = mg_mesh%level(ilev)%kvel(k,ivt)
    if (l.ne.0) then
     if (MarkerE(l).ne.1) MarkerE(l) = 3
    end if
   end do
  end do
 end if
end do


END SUBROUTINE SetUpMarker
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GetParameters()
CHARACTER*(10) :: cParam
integer ilong

!  CHARACTER*(200) :: cmd
 LOGICAL bExist

 OPEN(1,file='param.cfg')

 READ(1,*) cIntputFolder
 
 READ(1,*) cProjectFolder
 cProjectFolder = ADJUSTL(TRIM(cIntputFolder))//'/'//ADJUSTL(TRIM(cProjectFolder))
 iLong = LEN(ADJUSTL(TRIM(cProjectFolder)))+1
 WRITE(cProjectFolder(iLong:),"(A)") "/"
 WRITE(*,*) adjustl(trim(cProjectFolder))
 
 READ(1,*) cShortProjectFile
 cProjectFile = adjustl(trim(cProjectFolder))//adjustl(trim(cShortProjectFile))
 WRITE(*,*) adjustl(trim(cProjectFile))
 
 CALL ExtractMeshfile()
 
! READ(1,*) mg_Mesh%nlmax
 mg_Mesh%nlmax = 1
 mg_Mesh%nlmin = 1
 WRITE(*,*) 'Min and Max levels: ', mg_Mesh%nlmin,mg_Mesh%nlmax
 mg_Mesh%nlmin = 1
 READ(1,*) dCGALtoRealFactor
 WRITE(*,'(A,ES12.4)') ' CGAL Scaling factor: ', dCGALtoRealFactor
 READ(1,*) cOutputFolder
 WRITE(*,*) 'Output Folder: "'//adjustl(trim(cOutputFolder))//'"'
 READ(1,*) lTriOutputLevel
 WRITE(*,*) 'Outputlevel for the ".tri" file: ', lTriOutputLevel
 READ(1,*) lVTUOutputLevel
 WRITE(*,*) 'Outputlevel for the ".vtu" file: ', lVTUOutputLevel
 READ(1,*) RefinementThickness
 WRITE(*,*) 'Refinement thickness is set to be: ', RefinementThickness
 READ(1,*) initfield
 WRITE(*,*) 'Refinement field is set to be: ', initfield
 READ(1,*) MeshOutputScaleFactor
 WRITE(*,*) 'MeshOutputScaleFactor is set to be: ', MeshOutputScaleFactor
 CLOSE(1)

#ifdef __GFORTRAN__
 INQUIRE(file=adjustl(trim(cOutputFolder)),EXIST=bExist)
#elif defined __INTEL_COMPILER
 INQUIRE(directory=adjustl(trim(cOutputFolder)),EXIST=bExist)
#else
 INQUIRE(file=adjustl(trim(cOutputFolder)),EXIST=bExist)
#endif

 if (.not.bExist) then
!   cmd='mkdir '//adjustl(trim(cOutputFolder))
  CALL system('mkdir '//adjustl(trim(cOutputFolder)))
 end if

END SUBROUTINE GetParameters
!----------------------------------------------------------
SUBROUTINE ExtractMeshfile()
IMPLICIT NONE
INTEGER LenStr,iEnd
CHARACTER(LEN=200) :: string,cFile
logical :: bFound=.false.

OPEN(unit=2,file=adjustl(trim(cProjectFile)))

DO
 READ(2,FMT='(200A)',IOSTAT=iEnd) string
 IF (iEnd.EQ.-1) EXIT
 LenStr = LEN(ADJUSTL(TRIM(string)))
 IF (LenStr.gt.4) THEN
  cFile = ADJUSTL(TRIM(string))
  IF (cFile(LenStr-3:LenStr).EQ.".tri") THEN
   cProjectGridFile = adjustl(trim(cFile))
!    cProjectGridFile = adjustl(trim(cProjectFolder))//'/'//adjustl(trim(cFile))
   Write(*,*) 'Mesh file: "'//ADJUSTL(TRIM(cProjectGridFile))//'"'
   bFound=.true.
   exit
  END IF
 END IF
END DO

CLOSE(2)
if (.not.bFound) then
 Write(*,*) 'Mesh file was NOT found in Project file! '
end if

END SUBROUTINE ExtractMeshfile
! ----------------------------------------------
SUBROUTINE RotatePatch(C,ID)
integer ID(8),IDR(8),IDR1(8),IDR2(8),IDR3(8)
logical :: C(8),R1(8),R2(8),R3(8),P(8)
integer i,iSum,jSum,j,jSum_init

jSum_init = 64*8
jSum = jSum_init
 
 DO i=1,8
  IF (C(i)) THEN
   
   R1 = C(RROT(:,i))
   IDR1 = ID(RROT(:,i))
   iSum = 0
   do j=1,8
    if (R1(j)) iSum = iSum + j*j
   end do
   if (iSum.lt.jSum) then
    P = R1
    IDR = IDR1
    jSum = iSum
   end if
   
   R2 = R1(R2ROT)
   IDR2 = IDR1(R2ROT)
   iSum = 0
   do j=1,8
    if (R2(j)) iSum = iSum + j*j
   end do
   if (iSum.lt.jSum) then
    P = R2
    IDR = IDR2
    jSum = iSum
   end if
   
   R3 = R1(R3ROT)
   IDR3 = IDR1(R3ROT)
   iSum = 0
   do j=1,8
    if (R3(j)) iSum = iSum + j*j
   end do
   if (iSum.lt.jSum) then
    P = R3
    IDR = IDR3
    jSum = iSum
   end if
   
  END IF
 END DO

 if (jSum.lt.jSum_init) then
  ID = IDR
  C=P
 end if
!   WRITE(*,'(8L,A,8L,A,8L,A,8L)') C , " :: ", R1,  " | ", R2,  " | ", R3

END SUBROUTINE RotatePatch
!
!------------------------------------------------------------------
!
SUBROUTINE DetermineTemplate(C,iT)
logical C(8),R(8)
integer iT
integer i
logical bSame

iT = 0

do i=1,22
 R = templates(:,i)
 bSame=.true.
 call ComparePatch(C,R,bSame)
!  write(*,*) C,":",R,":",bSame
 if (.not.bSame) then
  iT = i
  RETURN
 end if
end do

END SUBROUTINE DetermineTemplate
!
!------------------------------------------------------------------
!
SUBROUTINE ComparePatch(C,R,b)
logical :: C(8),R(8)
logical b
integer i,j

j = 0
DO i=1,8
 if ((C(i).and.R(i)).or.((.not.C(i)).and.(.not.R(i)))) j = j + 1
END DO
 
if (j.eq.8) b=.false.
!  write(*,*) C,":",R,":",j
 
END SUBROUTINE ComparePatch
!
!------------------------------------------------------------------
!
SUBROUTINE CleanUpPatches()

integer iel,jel,i,j,jj,iat,nUniquePoints,nAllPoints,nAllElements,il
integer nTotalPoints,nTotalElements
real*8 Pi(3),Pj(3),dist
logical bFound

integer I8(8)
real*8 P8(3,8),dVol,dTol

ilev = mg_Mesh%nlmax
nel = mg_mesh%level(ilev)%nel
nvt = mg_mesh%level(ilev)%nvt

nTotalPoints = 0
nTotalElements = 0
nAllPoints = 0
nAllElements = 0

do iel=1,nel

 nAllPoints = 0
 nAllElements = 0

 if (allocated(myRF(iel)%myRF)) then
  CALL CountSubPatch(myRF(iel))
 else
  nAllPoints = myRF(iel)%nOfVert
  nAllElements = myRF(iel)%nOfElem
 end if
 
 allocate(myRF(iel)%dUniquedCoor(3,nAllPoints))
 allocate(myRF(iel)%kUniqueElem(8,nAllElements))
 
 nUniquePoints = 0
 nAllElements  = 0
 if (allocated(myRF(iel)%myRF)) then
  il = 0
  CALL CleanUpSubPatch(myRF(iel))
 else
  nUniquePoints = myRF(iel)%nOfVert
  nAllElements = myRF(iel)%nOfElem
  myRF(iel)%dUniquedCoor(1:3,1:nUniquePoints) = myRF(iel)%dcoor(1:3,1:nUniquePoints)
  myRF(iel)%kUniqueElem(1:8,1:nAllElements) = myRF(iel)%kvert(1:8,1:nAllElements)
  
  I8 = mg_mesh%level(ilev)%kvert(:,iel)
  do jj=1,8
   P8(:,jj) = mg_mesh%level(ilev)%dcorvg(:,I8(jj))
  end do
  CALL getVol(P8,dVol)
  dTol = 0.01d0*dVol**(0.333333d0)
  myRF(iel)%dTol = min(myRF(iel)%dTol,dTol)
 end if

 myRF(iel)%nUniquePoints = nUniquePoints
 myRF(iel)%nUniqueElems  = nAllElements
 nTotalPoints   = nTotalPoints + nUniquePoints
 nTotalElements = nTotalElements + nAllElements
!  write(*,*) myRF(iel)%dTol

end do

write(*,*) "nTotalElements,nTotalPoints",nTotalElements,nTotalPoints

! nTotalPoints = 0
! nTotalElements = 0
! do iel=1,nel
!  nTotalPoints   = nTotalPoints + myRF(iel)%nUniquePoints
!  nTotalElements = nTotalElements + myRF(iel)%nUniqueElems
! enddo
! write(*,*) "nTotalElements,nTotalPoints",nTotalElements,nTotalPoints

 CONTAINS
 
 RECURSIVE SUBROUTINE CountSubPatch(M)
 type(RefinerMesh) :: M
 integer jj

 if (allocated(M%myRF)) then
  DO jj=1,M%nOfElem
   call CountSubPatch(M%myRF(jj))
  END DO
 else
  nAllPoints = nAllPoints + M%nOfVert
  nAllElements = nAllElements + M%nOfElem
 end if
 
 END SUBROUTINE CountSubPatch
! 
!---------------------------------------------------------------
!
 RECURSIVE SUBROUTINE CleanUpSubPatch(M)
 type(RefinerMesh) :: M
 integer ii,jj,kk,I8(8)
 real*8 Pi(3),Pj(3),dist,P8(3,8),dVol
 real*8 :: dTol=1d-4
 integer iFound
 integer, allocatable :: iBuff(:)
 
 allocate(iBuff(M%nOfVert))
 
 if (allocated(M%myRF)) then
  il = il + 1
  DO jj=1,M%nOfElem
   call CleanUpSubPatch(M%myRF(jj))
  END DO
 else
 
  I8 = mg_mesh%level(ilev)%kvert(:,iel)
  do jj=1,8
   P8(:,jj) = mg_mesh%level(ilev)%dcorvg(:,I8(jj))
  end do
  CALL getVol(P8,dVol)
  dTol = 0.01d0*dVol**(0.333333d0)
  myRF(iel)%dTol = min(myRF(iel)%dTol,dTol)
  
  DO jj = 1,M%nOfVert
   Pi = M%dcoor(:,jj)
   
   iFound = 0
   DO kk = 1, nUniquePoints
    Pj = myRF(iel)%dUniquedCoor(:,kk)
    dist =  sqrt((Pi(1)-Pj(1))**2d0 + (Pi(2)-Pj(2))**2d0 + (Pi(3)-Pj(3))**2d0)
    if (dist.lt.dTol) then
     iFound = kk
     exit
    END IF
   END DO
   
   if (iFound.gt.0) then
     iBuff(jj) = iFound
   else
     nUniquePoints = nUniquePoints + 1
     myRF(iel)%dUniquedCoor(:,nUniquePoints) = Pi
     iBuff(jj) = nUniquePoints
   end if
  END DO
  
  DO jj=1,M%nOfElem
   do kk=1,8
    myRF(iel)%kUniqueElem(kk,nAllElements+jj) =  iBuff(M%kvert(kk,jj))
   end do
  END DO
  
  nAllElements = nAllElements + M%nOfElem
  
 end if
 
 deallocate(iBuff)
 
!   WRITE(*,*) "<zsdsaddasdsa sadsa sa "
!   do jj=1,nUniquePoints
!     write(*,'(3ES12.3)') myRF(iel)%dUniquedCoor(:,jj)
!   end do
!   pause
!   do jj=1,nAllElements
!     write(*,'(8I8)') myRF(iel)%kUniqueElem(:,jj)
!   end do

 
 END SUBROUTINE CleanUpSubPatch
 
END SUBROUTINE CleanUpPatches
!
!------------------------------------------------------------------
!
SUBROUTINE CleanUpMesh()

logical, allocatable :: bDone(:)
integer iel,jel,i,j,iat,ivt,ivel
real*8 Pi(3),Pj(3),dist
integer iFound

ilev = mg_Mesh%nlmax

nel = mg_mesh%level(ilev)%nel
nvt = mg_mesh%level(ilev)%nvt

allocate(bDone(nel))

bDone = .false.
nUniquePoints = 0
nUniqueElems = 0

do iel=1,nel

 do i=1,myRF(iel)%nUniquePoints
 
  Pi = myRF(iel)%dUniquedCoor(:,i)
  
  iFound = 0
  
  do iat=1,8
   ivt = mg_mesh%level(ilev)%kvert(iat,iel)
   do ivel=1,mg_mesh%level(ilev)%nvel
    jel =  mg_mesh%level(ilev)%kvel(ivel,ivt)
    if (jel.gt.0) then
     if (bDone(jel)) then
      do j=1,myRF(jel)%nUniquePoints
       Pj = myRF(jel)%dUniquedCoor(:,j)
       dist =  sqrt((Pi(1)-Pj(1))**2d0 + (Pi(2)-Pj(2))**2d0 + (Pi(3)-Pj(3))**2d0)
       if (dist.lt.myRF(jel)%dTol) then
        iFound = 1
        GOTO 1
       end if
      end do
     end if
    end if
   end do
  end do
1 continue
  if (iFound.ge.1) then
  else
   nUniquePoints = nUniquePoints + 1 
  end if
 end do
 bDone(iel) = .true.
 nUniqueElems = nUniqueElems + myRF(iel)%nUniqueElems
end do

write(*,*) "nUniquePoints=",nUniquePoints,"nUniqueElems=",nUniqueElems

allocate(MergedMeshCoor(3,nUniquePoints))
allocate(MergedMeshElem(8,nUniqueElems))

bDone = .false.
nUniquePoints = 0
nUniqueElems = 0

do iel=1,nel

 allocate(myRF(iel)%PointerToMerged(myRF(iel)%nUniquePoints))

 do i=1,myRF(iel)%nUniquePoints
 
  Pi = myRF(iel)%dUniquedCoor(:,i)
  
  iFound = 0
  
  do iat=1,8
   
   ivt = mg_mesh%level(ilev)%kvert(iat,iel)
   
   do ivel=1,mg_mesh%level(ilev)%nvel
   
      jel =  mg_mesh%level(ilev)%kvel(ivel,ivt)
      if (jel.gt.0) then
       if (bDone(jel)) then
        do j=1,myRF(jel)%nUniquePoints
         Pj = myRF(jel)%dUniquedCoor(:,j)
         dist =  sqrt((Pi(1)-Pj(1))**2d0 + (Pi(2)-Pj(2))**2d0 + (Pi(3)-Pj(3))**2d0)
         if (dist.lt.myRF(jel)%dTol) then
          iFound = myRF(jel)%PointerToMerged(j)
          if (iFound.lt.1) then
           write(*,*) 'o uuu ...'
          end if
          GOTO 2
         end if
        end do
       end if
      end if
   
   end do
  end do

2 continue

  if (iFound.ge.1) then
!    write(*,*) "iFound:", iFound
   myRF(iel)%PointerToMerged(i) = iFound
  else
   nUniquePoints = nUniquePoints + 1 
   MergedMeshCoor(:,nUniquePoints) = Pi
   myRF(iel)%PointerToMerged(i) = nUniquePoints
  end if
 end do
 
 do i=1,myRF(iel)%nUniqueElems
  do j=1,8
   ivt = myRF(iel)%kUniqueElem(j,i)
   MergedMeshElem(j,nUniqueElems+i) = myRF(iel)%PointerToMerged(ivt)
  end do
 end do
 
 bDone(iel) = .true.
 nUniqueElems = nUniqueElems + myRF(iel)%nUniqueElems
 
end do


END SUBROUTINE CleanUpMesh
!
!------------------------------------------------------------------
!
SUBROUTINE getVol(P,dVol)

REAL*8 P(3,8),dVol
REAL*8 :: A1 = 1d0/6d0

 dVol=A1*((DABS((P(1,4)-P(1,1))*(P(2,4)-P(2,3))*(P(3,4)-P(3,8))+(P(2,4)-P(2,1))*  &
      (P(3,4)-P(3,3))*(P(1,4)-P(1,8))+(P(3,4)-P(3,1))*(P(1,4)-P(1,3))*(P(2,4)-P(2,8))- &
      (P(1,4)-P(1,8))*(P(2,4)-P(2,3))*(P(3,4)-P(3,1))-(P(2,4)-P(2,8))*(P(3,4)-P(3,3))* &
      (P(1,4)-P(1,1))-(P(3,4)-P(3,8))*(P(1,4)-P(1,3))*(P(2,4)-P(2,1))))+       &
      (DABS((P(1,2)-P(1,3))*(P(2,2)-P(2,1))*(P(3,2)-P(3,6))+(P(2,2)-P(2,3))*   &
      (P(3,2)-P(3,1))*(P(1,2)-P(1,6))+(P(3,2)-P(3,3))*(P(1,2)-P(1,1))*(P(2,2)-P(2,6))- &
      (P(1,2)-P(1,6))*(P(2,2)-P(2,1))*(P(3,2)-P(3,3))-(P(2,2)-P(2,6))*(P(3,2)-P(3,1))* &
      (P(1,2)-P(1,3))-(P(3,2)-P(3,6))*(P(1,2)-P(1,1))*(P(2,2)-P(2,3))))+       &
      (DABS((P(1,5)-P(1,8))*(P(2,5)-P(2,6))*(P(3,5)-P(3,1))+(P(2,5)-P(2,8))*   &
      (P(3,5)-P(3,6))*(P(1,5)-P(1,1))+(P(3,5)-P(3,8))*(P(1,5)-P(1,6))*(P(2,5)-P(2,1))- &
      (P(1,5)-P(1,1))*(P(2,5)-P(2,6))*(P(3,5)-P(3,8))-(P(2,5)-P(2,1))*(P(3,5)-P(3,6))* &
      (P(1,5)-P(1,8))-(P(3,5)-P(3,1))*(P(1,5)-P(1,6))*(P(2,5)-P(2,8))))+       &
      (DABS((P(1,7)-P(1,6))*(P(2,7)-P(2,8))*(P(3,7)-P(3,3))+(P(2,7)-P(2,6))*   &
      (P(3,7)-P(3,8))*(P(1,7)-P(1,3))+(P(3,7)-P(3,6))*(P(1,7)-P(1,8))*(P(2,7)-P(2,3))- &
      (P(1,7)-P(1,3))*(P(2,7)-P(2,8))*(P(3,7)-P(3,6))-(P(2,7)-P(2,3))*(P(3,7)-P(3,8))* &
      (P(1,7)-P(1,6))-(P(3,7)-P(3,3))*(P(1,7)-P(1,8))*(P(2,7)-P(2,6))))+       &
      (DABS((P(1,1)-P(1,3))*(P(2,1)-P(2,8))*(P(3,1)-P(3,6))+(P(2,1)-P(2,3))*   &
      (P(3,1)-P(3,8))*(P(1,1)-P(1,6))+(P(3,1)-P(3,3))*(P(1,1)-P(1,8))*(P(2,1)-P(2,6))- &
      (P(1,1)-P(1,6))*(P(2,1)-P(2,8))*(P(3,1)-P(3,3))-(P(2,1)-P(2,6))*(P(3,1)-P(3,8))* &
      (P(1,1)-P(1,3))-(P(3,1)-P(3,6))*(P(1,1)-P(1,8))*(P(2,1)-P(2,3)))))
   dVol = abs(dVol)

END SUBROUTINE getVol
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

END MODULE MeshRefDef
