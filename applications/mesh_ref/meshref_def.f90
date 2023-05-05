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
type(tMultiMesh),save :: mg_ReducedMesh
type(tMultiMesh),save :: mg_ParticleMesh


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
 CALL Initfield0(markerE,mg_mesh%level(ilev)%kvert,mg_mesh%level(ilev)%dcorvg,nel,RefinementThickness)
end if

IF (initfield.eq.1) then
 CALL Initfield1(markerE,mg_mesh%level(ilev)%kvert,mg_mesh%level(ilev)%dcorvg,nel,RefinementThickness)
end if

IF (initfield.eq.2) then
 CALL Initfield2(markerE,mg_mesh%level(ilev)%kvert,mg_mesh%level(ilev)%dcorvg,nel,RefinementThickness)
end if

IF (initfield.eq.3) then
 CALL Initfield3(markerE,mg_mesh%level(ilev)%kvert,mg_mesh%level(ilev)%dcorvg,nel,RefinementThickness)
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

 OPEN(1,file='param_meshref.cfg')

!  READ(1,*) cIntputFolder
 
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
!  READ(1,*) cOutputFolder
 cOutputFolder = adjustl(trim(cIntputFolder))
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
SUBROUTINE CleanUpSmallPatch(RFin,RFOut,E)
implicit none
type(RefinerMesh) RFin,RFOut
integer i,j,iP
real*8 P(3),Q(3),dVol,dTol,Diff
real*8 E(3,8)
integer, allocatable :: ihelp(:)

CALL getVol(E,dVol)
dTol = 0.01d0*dVol**(0.333333d0)

allocate(iHelp(RFin%nOfVert))
iHelp = 0

iP = 0

do i=1,RFin%nOfVert
 P = RFin%dcoor(:,i)
 if (iHelp(i).eq.0) then
  iP = iP + 1
  iHelp(i) = iP
  
  DO j=i+1,RFin%nOfVert
   Q = RFin%dcoor(:,j)
   DIFF = SQRT((P(1)-Q(1))**2d0 + (P(2)-Q(2))**2d0 + (P(3)-Q(3))**2d0)
   if (DIFF.lt.dTol) THEN
    iHelp(j) = iP
   end if
  END DO
  
 end if
end do

!write(*,*) iHelp
RFOut%nOfVert = iP
RFOut%nOfElem = RFin%nOfElem
allocate(RFOut%dcoor(3,RFOut%nOfVert))
allocate(RFOut%kvert(8,RFOut%nOfElem))
allocate(RFOut%knpr(RFOut%nOfVert))

do i=1,RFin%nOfVert
 RFOut%dcoor(:,iHelp(i)) = RFIn%dcoor(:,i)
end do

do i=1,RFin%nOfElem
 DO j=1,8
  RFOut%kvert(j,i) = iHelp(RFIn%kvert(j,i))
 END DO
end do

do i=1,RFin%nOfVert
 RFOut%knpr(iHelp(i)) = RFIn%knpr(i)
end do

! do i=1,RFIn%nOfVert
!  write(*,*) RFIn%dcoor(:,i) 
! end do

!write(*,*)  'RFOut%nOfVert, RFOut%nOfElem'
!write(*,*)  RFOut%nOfVert, RFOut%nOfElem

!do i=1,RFOut%nOfVert
! write(*,'(I8,3ES12.4)') i,RFOut%dcoor(:,i) 
!end do

!do i=1,RFOut%nOfElem
!  write(*,'(9I8)') i,RFOut%kvert(:,i) 
!end do

!do i=1,RFOut%nOfVert
! write(*,*) i,RFOut%knpr(i)
!end do

!pause

END SUBROUTINE CleanUpSmallPatch
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
 allocate(myRF(iel)%kUniqueKnpr(nAllPoints))
 
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
!   write(*,*) nAllElements
!   write(*,*) myRF(iel)%knpr(1:nAllPoints)
  myRF(iel)%kUniqueKnpr(1:nUniquePoints) = myRF(iel)%knpr(1:nUniquePoints)
    
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
allocate(MergedMeshknpr(nUniquePoints))

bDone = .false.
nUniquePoints = 0
nUniqueElems = 0

do iel=1,nel

 allocate(myRF(iel)%PointerToMerged(myRF(iel)%nUniquePoints))

 do i=1,myRF(iel)%nUniquePoints
 
  Pi = myRF(iel)%dUniquedCoor(:,i)
  !!!!!! CHECKPOINT !!!!!
!   if (myRF(iel)%kUniqueKnpr(i).ne.0) write(*,*) myRF(iel)%kUniqueKnpr(i)
  !!!!!! CHECKPOINT !!!!!
  
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
   MergedMeshKnpr(nUniquePoints) = myRF(iel)%kUniqueKnpr(i)
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
! ----------------------------------------------
SUBROUTINE CreateReducedMesh
use geometry_processing, only : GetDistToSTL
integer, allocatable :: iVertKeep(:),iElemKeep(:)
real*8, allocatable :: Distance(:)
integer i,j,k,jElem,jVert,ivt
real*8 :: pC(3),dist,distC,P(3),dMaxDist,dMinDist
real*8 :: p8(3,8),D8(8)
logical :: bKeep,bPrint

allocate(iVertKeep(nUniquePoints))
allocate(iElemKeep(nUniqueElems))
allocate(Distance(nUniquePoints))
Distance = 1d8

iVertKeep=0
iElemKeep=0

jElem = 0
jVert = 0

DO i=1,nUniqueElems
 
 pC = 0d0
 do j=1,8
  ivt = MergedMeshElem(j,i)
  P = MergedMeshCoor(:,ivt)
  pC = pC + 0.125d0*P
 end do
 
 dMaxDist = -1d8
 dMinDist = +1d8
 do j=1,8
  ivt = MergedMeshElem(j,i)
  P = MergedMeshCoor(:,ivt)
  distC = SQRT((PC(1)-P(1))**2d0 + (PC(2)-P(2))**2d0 + (PC(3)-P(3))**2d0)
  if (distC.gt.dMaxDist) dMaxDist = distC
  if (distC.lt.dMaxDist) dMinDist = distC
 end do

 CALL GetDistToSTL(pC(1),pC(2),pC(3),1,distC,.TRUE.)
 
 if (distC.lt.dMaxDist) then
  ! Clear removal --> element is too far
  IF (distC.lt.-dMaxDist) THEN
   jElem = jElem + 1
   iElemKeep(i) = jElem
   iVertKeep(MergedMeshElem(:,i)) = 1
  else
  
   bKeep = .FALSE.
   
   IF (distC.lt.1d0*dMinDist) then
   
    bKeep = .TRUE.
    jElem = jElem + 1
    iElemKeep(i) = jElem
    iVertKeep(MergedMeshElem(:,i)) = 1

   ELSE
   
    do j=1,8
     ivt = MergedMeshElem(j,i)
     P = MergedMeshCoor(:,ivt)
     IF (abs(Distance(ivt)).gt.1d7) THEN
      CALL GetDistToSTL(p(1),p(2),p(3),1,dist,.TRUE.)
      Distance(ivt) = dist
     else
      dist = Distance(ivt)
     END IF
     
     IF (dist.lt.0d0) then
      bKeep = .TRUE.
      jElem = jElem + 1
      iElemKeep(i) = jElem
      iVertKeep(MergedMeshElem(:,i)) = 1
      exit
     END IF
     
    end do
    
   END IF
   
   IF (.not.bKeep) THEN
    P8(1,:) = MergedMeshCoor(1,MergedMeshElem(:,i)) 
    P8(2,:) = MergedMeshCoor(2,MergedMeshElem(:,i)) 
    P8(3,:) = MergedMeshCoor(3,MergedMeshElem(:,i)) 
    D8(:)   = Distance(MergedMeshElem(:,i))
    bPrint = .false.
!     IF (i.eq.93265) THEN
!      bPrint = .true.
!      WRITE(*,*) "bKeep",bKeep
!      DO k=1,8
!       WRITE(*,'(3ES12.4,A,ES12.4)') P8(:,k)," : ",D8(k)
!      END DO
!     END IF
    ! we should figure out if the hex intersects the triangulation
    CALL TriangulationIntersectsnElement(P8,D8,PC,distC,bKeep,bPrint)
!     IF (i.eq.62402) THEN
!      pause
!     end if
    IF (bKeep) THEN
     write(*,*) "Element ",i," has been determined to be added to list"
     bKeep = .TRUE.
     jElem = jElem + 1
     iElemKeep(i) = jElem
     iVertKeep(MergedMeshElem(:,i)) = 1
    END IF
   END IF
   
  END IF
 end if
 
end do

DO i=1,nUniquePoints
 if (iVertKeep(i).gt.0) then
  jVert = jVert + 1
  iVertKeep(i) = jVert
 end if
END DO

allocate(ReducedMeshElem(8,jElem))
allocate(ReducedMeshCoor(3,jVert))

DO i=1,nUniqueElems
 IF (iElemKeep(i).gt.0) then
  do j=1,8
   ReducedMeshElem(j,iElemKeep(i)) = iVertKeep(MergedMeshElem(j,i))
  end do
 END IF 
end do

DO i=1,nUniquePoints
 IF (iVertKeep(i).gt.0) then
  do j=1,3
   ReducedMeshCoor(j,iVertKeep(i)) = MergedMeshCoor(j,i)
  end do
 END IF 
end do

nReducedElems = jElem
nReducedPoints = jVert
deallocate(iVertKeep)
deallocate(iElemKeep)
deallocate(Distance)

END SUBROUTINE CreateReducedMesh
! ----------------------------------------------
SUBROUTINE CreateCleanReducedMesh
use geometry_processing, only : GetDistToSTL
use Sigma_User, only : myProcess
implicit none
integer, allocatable :: iVertKeep(:),iElemKeep(:)
integer i,j,k,ia,jat,jElem,iarea,jVert,ivt,iStartElem,nBCFace,iInflow,nInFlow,nOutFlow
real*8 :: dN(3),pC(3),dist,P(3),Q(3),dMinDist,minDistP,maxDistP
logical bInflowArea,bQuad
INTEGER NeighA(4,6),iQuad(4)
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

real*8 :: DA(3), DB(3), dCenter(3), dAux1, dAux2 
real*8 :: dAdC, dBdC, dPdA, dPdB, dR

nel = mg_ReducedMesh%level(1)%nel
nvt = mg_ReducedMesh%level(1)%nvt
nat = mg_ReducedMesh%level(1)%nat

allocate(ReducedMeshBC(nat))
ReducedMeshBC = -2

nBCFace = 0
DO i=1,nel
 DO ia= 1,6
  j = mg_ReducedMesh%level(1)%kadj(ia,i)
  IF (j.eq.0) then
   ReducedMeshBC(mg_ReducedMesh%level(1)%karea(ia,i)) = -1
   nBCFace = nBCFace + 1
  end if
 end do
end do

nOutFlow = 0
nInFlow = 0

DO i=1,nel
 DO ia= 1,6
  iarea = mg_ReducedMesh%level(1)%karea(ia,i)
  j = mg_ReducedMesh%level(1)%kadj(ia,i)
  IF (j.eq.0) then
  
   P = 0d0
   DO ivt = 1,4
    k = mg_ReducedMesh%level(1)%kvert(neighA(ivt,ia),i)
    P  = P + 0.25d0*mg_ReducedMesh%level(1)%dcorvg(:,k)
   END DO
   
   minDistP = +1d8
   maxDistP = -1d8
   DO ivt = 1,4
    k = mg_ReducedMesh%level(1)%kvert(neighA(ivt,ia),i)
    Q = mg_ReducedMesh%level(1)%dcorvg(:,k)
    dist = sqrt((P(1)-Q(1))**2d0 + (P(2)-Q(2))**2d0 + (P(3)-Q(3))**2d0) 
    IF (dist.lt.minDistP) minDistP = dist
    IF (dist.gt.maxDistP) maxDistP = dist 
   END DO

   IF (abs(P(3)-OverallBoundingBox(3,2)).lt.0.5d0*minDistP) THEN
     nOutFlow = nOutFlow + 1
     ReducedMeshBC(iarea) = 0
   end if
   DO iInflow=1,myProcess%nOfInflows
    pC = myProcess%myInflow(iInflow)%center/MeshOutputScaleFactor
    
    bInflowArea = .false.
    IF ((abs(P(1)-OverallBoundingBox(1,1)).lt.0.5d0*minDistP).or.&
        (abs(P(1)-OverallBoundingBox(1,2)).lt.0.5d0*minDistP).or.&
        (abs(P(2)-OverallBoundingBox(2,1)).lt.0.5d0*minDistP).or.&
        (abs(P(2)-OverallBoundingBox(2,2)).lt.0.5d0*minDistP).or.&
        (abs(P(3)-OverallBoundingBox(3,1)).lt.0.5d0*minDistP)) THEN

     DO ivt = 1,4
      k = mg_ReducedMesh%level(1)%kvert(neighA(ivt,ia),i)
      Q = mg_ReducedMesh%level(1)%dcorvg(:,k)
      dist = sqrt((Q(1)-pC(1))**2d0 + (Q(2)-pC(2))**2d0 + (Q(3)-pC(3))**2d0) 
      
      IF (myProcess%myInflow(iInflow)%iBCtype.eq.2) then
       IF (dist.le.myProcess%myInflow(iInflow)%outerradius/MeshOutputScaleFactor+0.5d0*minDistP.and.&
           dist.gt.myProcess%myInflow(iInflow)%innerradius/MeshOutputScaleFactor-0.5d0*minDistP) THEN
        dN = myProcess%myInflow(iInflow)%normal
        jat =0
        IF (abs(P(1)-OverallBoundingBox(1,1)).lt.0.5d0*minDistP) jat = 1
        IF (abs(P(1)-OverallBoundingBox(1,2)).lt.0.5d0*minDistP) jat = 2
        IF (abs(P(2)-OverallBoundingBox(2,1)).lt.0.5d0*minDistP) jat = 3
        IF (abs(P(2)-OverallBoundingBox(2,2)).lt.0.5d0*minDistP) jat = 4
        IF (abs(P(3)-OverallBoundingBox(3,1)).lt.0.5d0*minDistP) jat = 5
        
        IF (abs(dN(1)).gt.0.99d0.and.(jat.eq.1.or.jat.eq.2)) bInflowArea = .true.
        IF (abs(dN(2)).gt.0.99d0.and.(jat.eq.3.or.jat.eq.4)) bInflowArea = .true.
        IF (abs(dN(3)).gt.0.99d0.and.(jat.eq.5)) bInflowArea = .true.
       END IF
      END IF
      
      IF (myProcess%myInflow(iInflow)%iBCtype.eq.1.or.&
          myProcess%myInflow(iInflow)%iBCtype.eq.3.or.&
          myProcess%myInflow(iInflow)%iBCtype.eq.4) then
       IF (dist.le.myProcess%myInflow(iInflow)%outerradius/MeshOutputScaleFactor+0.5d0*minDistP) THEN
        dN = myProcess%myInflow(iInflow)%normal
        jat =0
        IF (abs(P(1)-OverallBoundingBox(1,1)).lt.0.5d0*minDistP) jat = 1
        IF (abs(P(1)-OverallBoundingBox(1,2)).lt.0.5d0*minDistP) jat = 2
        IF (abs(P(2)-OverallBoundingBox(2,1)).lt.0.5d0*minDistP) jat = 3
        IF (abs(P(2)-OverallBoundingBox(2,2)).lt.0.5d0*minDistP) jat = 4
        IF (abs(P(3)-OverallBoundingBox(3,1)).lt.0.5d0*minDistP) jat = 5
        
        IF (abs(dN(1)).gt.0.99d0.and.(jat.eq.1.or.jat.eq.2)) bInflowArea = .true.
        IF (abs(dN(2)).gt.0.99d0.and.(jat.eq.3.or.jat.eq.4)) bInflowArea = .true.
        IF (abs(dN(3)).gt.0.99d0.and.(jat.eq.5)) bInflowArea = .true.
       END IF
      END IF

      IF (myProcess%myInflow(iInflow)%iBCType.EQ.5) THEN
         k = mg_ReducedMesh%level(1)%kvert(neighA(ivt,ia),i)
         Q = mg_ReducedMesh%level(1)%dcorvg(:,k)
         
         dCenter = myProcess%myInflow(iInflow)%center/MeshOutputScaleFactor
         DA = myProcess%myInflow(iInflow)%midpointA/MeshOutputScaleFactor
         DB = myProcess%myInflow(iInflow)%midpointB/MeshOutputScaleFactor

         dAux1 = DOT_PRODUCT(DA-dCenter, DA-dCenter)**2&
                -DOT_PRODUCT(Q-dCenter, DA-dCenter)**2+0.5*minDistP
         dAux2 = DOT_PRODUCT(DB-dCenter, DB-dCenter)**2&
                -DOT_PRODUCT(Q-dCenter, DB-dCenter)**2+0.5*minDistP
         IF ( (dAux1.GE.0D0).and.(dAux2.GE.0D0) ) THEN
           bInflowArea = .true.
         END IF
       END IF
       IF (myProcess%myInflow(iInflow)%iBCType.EQ.6) THEN
           k = mg_ReducedMesh%level(1)%kvert(neighA(ivt,ia),i)
           Q = mg_ReducedMesh%level(1)%dcorvg(:,k)
           
           dCenter = myProcess%myInflow(iInflow)%center/MeshOutputScaleFactor
           DA = myProcess%myInflow(iInflow)%midpointA/MeshOutputScaleFactor
           DB = myProcess%myInflow(iInflow)%midpointB/MeshOutputScaleFactor

           dAdC = DOT_PRODUCT(dA-dCenter, DA-dCenter)
           dBdC = DOT_PRODUCT(dB-dCenter, DB-dCenter)

           dPdA = DOT_PRODUCT(Q-dCenter, DA-dCenter)
           dPdB = DOT_PRODUCT(Q-dCenter, DB-dCenter)

           dR = NORM2(DB-dCenter)

           if ( (dPdA**2.LE.dAdC**2+0.5*minDistP).and.(dPdB**2.LE.dBdC**2+0.5*minDistP) ) THEN
               bInflowArea = .true.
           elseif ( (dPdA.GE.dAdC).and.(NORM2(Q-DA).LE.dR+0.5*minDistP) ) THEN
               bInflowArea = .true. 
           elseif ( (dPdA**2.GE.dAdC**2).and.(NORM2(Q-2*dCenter+DA).LE.dR+0.5*minDistP) ) THEN
               bInflowArea = .true.
           end if
        END IF ! Check for boundary Type
      
     END DO
    END IF
     

    if (bInflowArea) THEN
     nInFlow = nInFlow + 1
     ReducedMeshBC(iarea) = iInflow
    END If 
   END DO
   
  end if
 end do
end do

DO i=1,nel
 DO ia= 1,6
  iarea = mg_ReducedMesh%level(1)%karea(ia,i)
  j = mg_ReducedMesh%level(1)%kadj(ia,i)
  IF (j.eq.0) then
  
   if (ReducedMeshBC(iarea).eq.-1) THEN
    P = 0d0
    DO ivt = 1,4
     k = mg_ReducedMesh%level(1)%kvert(neighA(ivt,ia),i)
     P  = P + 0.25d0*mg_ReducedMesh%level(1)%dcorvg(:,k)
    END DO
    
    jat =0
    IF (abs(P(1)-OverallBoundingBox(1,1)).lt.0.5d0*minDistP) jat = 1
    IF (abs(P(1)-OverallBoundingBox(1,2)).lt.0.5d0*minDistP) jat = 2
    IF (abs(P(2)-OverallBoundingBox(2,1)).lt.0.5d0*minDistP) jat = 3
    IF (abs(P(2)-OverallBoundingBox(2,2)).lt.0.5d0*minDistP) jat = 4
    IF (abs(P(3)-OverallBoundingBox(3,1)).lt.0.5d0*minDistP) jat = 5
    
    if (jat.gt.0) ReducedMeshBC(iarea) = 1000 + jat 
        
   end if
   
  end if
 end do
end do

WRITE(*,'(A,8I8)') "number of Boundary faces // nvt:",nBCFace,nat,nInFlow,nOutFlow

ALLOCATE(ParList%Wall(nvt))
ALLOCATE(ParList%SideWall(5,nvt))
ALLOCATE(ParList%Outflow(nvt))
ALLOCATE(ParList%Inflow(myProcess%nOfInflows,nvt))
ParList%SideWall = .false.
ParList%Wall = .false.
ParList%Outflow = .false.
ParList%Inflow = .false.

ALLOCATE(ParList%WallA(nat))
ALLOCATE(ParList%SideWallA(5,nat))
ALLOCATE(ParList%OutflowA(nat))
ALLOCATE(ParList%InflowA(myProcess%nOfInflows,nat))
ParList%SideWallA = .false.
ParList%WallA = .false.
ParList%OutflowA = .false.
ParList%InflowA = .false.

DO i=1,nel
 DO ia= 1,6
  iarea = mg_ReducedMesh%level(1)%karea(ia,i)
  IF (ReducedMeshBC(iarea).eq.-1) THEN
   ParList%Wall(mg_ReducedMesh%level(1)%kvert(neighA(:,ia),i)) = .true.
   ParList%WallA(iarea) = .true.
  END IF
  IF (ReducedMeshBC(iarea).gt.1000) THEN
   jat = ReducedMeshBC(iarea) - 1000
   ParList%SideWall(jat,mg_ReducedMesh%level(1)%kvert(neighA(:,ia),i)) = .true.
   ParList%SideWallA(jat,iarea) = .true.
  END IF
  IF (ReducedMeshBC(iarea).eq.0) THEN
   ParList%Outflow(mg_ReducedMesh%level(1)%kvert(neighA(:,ia),i)) = .true.
   ParList%OutflowA(iarea) = .true.
  END IF
  DO iInflow=1,myProcess%nOfInflows
   IF (ReducedMeshBC(iarea).eq.iInflow) THEN
!     WRITE(*,*) iInflow,mg_ReducedMesh%level(1)%kvert(neighA(:,ia),i)
    ParList%Inflow(iInflow,mg_ReducedMesh%level(1)%kvert(neighA(:,ia),i)) = .true.
    ParList%InflowA(iInflow,iarea) = .true.
   END IF
  END DO
 END DO
END DO

allocate(iVertKeep(nvt))
allocate(iElemKeep(nel))
iVertKeep=0
iElemKeep=0

jElem = 0
jVert = 0

dMinDist = +1d8
DO i=1,nel
 
 pC = 0d0
 do j=1,8
  ivt = mg_ReducedMesh%level(1)%kvert(j,i)
  P = mg_ReducedMesh%level(1)%dcorvg(:,ivt)
  pC = pC + 0.125d0*P
 end do
 
 CALL GetDistToSTL(pC(1),pC(2),pC(3),1,dist,.TRUE.)
 
 if (dist.lt.dMinDist) then
  iStartElem = i
  dMinDist=dist
 end if
 
end do

CALL MarkElems(iStartElem)

jElem = 0
DO i=1,nel
 IF (iElemKeep(i).gt.0) then
  jElem = jElem + 1
  iElemKeep(i) = jElem
  iVertKeep(mg_ReducedMesh%level(1)%kvert(:,i)) = 1
 END IF
END DO

DO i=1,nvt
 if (iVertKeep(i).gt.0) then
  jVert = jVert + 1
  iVertKeep(i) = jVert
 end if
END DO

WRITE(*,*)"Number of Marked Elems/Verts:",jElem,jVert
WRITE(*,*)"Number of un-Marked Elems/Verts:",nel-jElem,nvt-jVert

nReducedCleanPoints = jVert
nReducedCleanElems = jElem

allocate(ReducedCleanMeshElem(8,jElem))
allocate(ReducedCleanMeshCoor(3,jVert))

DO i=1,nel
 IF (iElemKeep(i).gt.0) then
  do j=1,8
   ReducedCleanMeshElem(j,iElemKeep(i)) = iVertKeep(mg_ReducedMesh%level(1)%kvert(j,i))
  end do
 END IF 
end do

DO i=1,nvt
 IF (iVertKeep(i).gt.0) then
  do j=1,3
   ReducedCleanMeshCoor(j,iVertKeep(i)) = mg_ReducedMesh%level(1)%dcorvg(j,i)
  end do
 END IF 
end do

ALLOCATE(ParCleanList%Wall(nReducedCleanPoints))
ALLOCATE(ParCleanList%SideWall(5,nReducedCleanPoints))
ALLOCATE(ParCleanList%Outflow(nReducedCleanPoints))
ALLOCATE(ParCleanList%Inflow(myProcess%nOfInflows,nReducedCleanPoints))
ParCleanList%Wall = .false.
ParCleanList%SideWall = .false.
ParCleanList%Outflow = .false.
ParCleanList%Inflow = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

ALLOCATE(ParCleanList%nIA(myProcess%nOfInflows))
ALLOCATE(ParCleanList%tIA(myProcess%nOfInflows))
ParCleanList%nIA = 0
ParCleanList%nSA = 0
ParCleanList%nOA = 0
ParCleanList%nWA = 0

DO i=1,nel
 DO ia= 1,6
  iarea = mg_ReducedMesh%level(1)%karea(ia,i)
  
  iQuad = [iVertKeep(mg_ReducedMesh%level(1)%kvert(neighA(1,ia),i)),&
           iVertKeep(mg_ReducedMesh%level(1)%kvert(neighA(2,ia),i)),&
           iVertKeep(mg_ReducedMesh%level(1)%kvert(neighA(3,ia),i)),&
           iVertKeep(mg_ReducedMesh%level(1)%kvert(neighA(4,ia),i))]
  bQuad = .false.
  if (iQuad(1).gt.0.and.iQuad(2).gt.0.and.iQuad(3).gt.0.and.iQuad(4).gt.0) bQuad = .true.

  if (bQuad) then
   if (ParList%WallA(iarea)) then 
    ParCleanList%nWA = ParCleanList%nWA + 1 
   end if
   
   if (ParList%OutflowA(iarea)) then 
    ParCleanList%nOA = ParCleanList%nOA + 1 
   end if
   
   DO jat=1,5
    if (ParList%SideWallA(jat,iarea)) then 
     ParCleanList%nSA(jat) = ParCleanList%nSA(jat)+ 1 
    end if
   END DO
   
   DO iInflow=1,myProcess%nOfInflows
    if (ParList%InflowA(iInflow,iarea)) then 
     ParCleanList%nIA(iInflow) = ParCleanList%nIA(iInflow)+ 1 
    end if
   END DO
   
  end if  
 end do  
end do

write(*,*) ParCleanList%nWA,ParCleanList%nOA,ParCleanList%nIA,ParCleanList%nSA

allocate(ParCleanList%tOA%i(4,ParCleanList%nOA))
allocate(ParCleanList%tWA%i(4,ParCleanList%nWA))
DO iInflow=1,myProcess%nOfInflows
 allocate(ParCleanList%tIA(iInflow)%i(4,ParCleanList%nIA(iInflow)))
END DO
DO jat=1,5
 allocate(ParCleanList%tSA(jat)%i(4,ParCleanList%nSA(jat)))
END DO
ParCleanList%nIA = 0
ParCleanList%nSA = 0
ParCleanList%nOA = 0
ParCleanList%nWA = 0

DO i=1,nel
 DO ia= 1,6
  iarea = mg_ReducedMesh%level(1)%karea(ia,i)
  
  iQuad = [iVertKeep(mg_ReducedMesh%level(1)%kvert(neighA(1,ia),i)),&
           iVertKeep(mg_ReducedMesh%level(1)%kvert(neighA(2,ia),i)),&
           iVertKeep(mg_ReducedMesh%level(1)%kvert(neighA(3,ia),i)),&
           iVertKeep(mg_ReducedMesh%level(1)%kvert(neighA(4,ia),i))]
  bQuad = .false.
  if (iQuad(1).gt.0.and.iQuad(2).gt.0.and.iQuad(3).gt.0.and.iQuad(4).gt.0) bQuad = .true.

  if (bQuad) then
   if (ParList%WallA(iarea)) then 
    ParCleanList%nWA = ParCleanList%nWA + 1 
    ParCleanList%tWA%i(:,ParCleanList%nWA) = iQuad
   end if
   
   if (ParList%OutflowA(iarea)) then 
    ParCleanList%nOA = ParCleanList%nOA + 1 
    ParCleanList%tOA%i(:,ParCleanList%nOA) = iQuad
   end if
   
   DO jat=1,5
    if (ParList%SideWallA(jat,iarea)) then 
     ParCleanList%nSA(jat) = ParCleanList%nSA(jat)+ 1 
     ParCleanList%tSA(jat)%i(:,ParCleanList%nSA(jat)) = iQuad
    end if
   END DO
   
   DO iInflow=1,myProcess%nOfInflows
    if (ParList%InflowA(iInflow,iarea)) then 
     ParCleanList%nIA(iInflow) = ParCleanList%nIA(iInflow)+ 1 
     ParCleanList%tIA(iInflow)%i(:,ParCleanList%nIA(iInflow)) = iQuad
    end if
   END DO
   
  end if  
 end do  
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1




DO i=1,nvt
 IF (iVertKeep(i).gt.0) then
  if (ParList%Wall(i))    THEN
   ParCleanList%Wall(iVertKeep(i)) = .TRUE.
  end if
  do jat=1,5
   if (ParList%SideWall(jat,i))    THEN
    ParCleanList%SideWall(jat,iVertKeep(i)) = .TRUE.
   end if
  end do
  if (ParList%Outflow(i)) THEN
   ParCleanList%Outflow(iVertKeep(i)) = .TRUE.
  end if
  DO iInflow=1,myProcess%nOfInflows
   IF (ParList%Inflow(iInflow,i)) THEN
!     write(*,*) iInflow,i
    ParCleanList%Inflow(iInflow,iVertKeep(i)) = .TRUE.
   END IF
  END DO
 END IF
END DO

deallocate(iVertKeep)
deallocate(iElemKeep)

 CONTAINS
 
RECURSIVE SUBROUTINE MarkElems(IE)
implicit none
INTEGER IE
INTEGER IAT,JEL

jElem = jElem + 1
iElemKeep(IE)=1
! write(*,*) 'marking', IE
! pause

DO iat= 1,6
 jel = mg_ReducedMesh%level(1)%kadj(iat,IE)
 IF (jel.ne.0) THEN
  IF (iElemKeep(jel).eq.0) then 
   CALL MarkElems(jel)
  END IF
 END IF
END DO

END SUBROUTINE MarkElems

END SUBROUTINE CreateCleanReducedMesh
! ----------------------------------------------
SUBROUTINE TriangulationIntersectsnElement(P8,D8,PC,DistC,bKeep,bP)
REAL*8 P8(3,8),D8(8),PC(3),DistC
LOGICAL bKeep,bP
LOGICAL bEdgeF,bEdge,bEdge1,bEdge2,bEdge3,bEdge4
REAL*8 dMinDist
INTEGER i,iMinDist
INTEGER neighV(3,8),Neigh1,Neigh2,Neigh3,NeighR1,NeighR2
REAL*8  PMIN1(3),PMIN2(3),DMIN1,DMIN2
DATA neighV /2,4,5, 1,3,6, 2,4,7, 1,3,8, 1,6,8, 5,7,2, 6,8,3, 5,7,4/

dMinDist = 1d8
do i=1,8
 if (D8(i).lt.dMinDist) THEN
  dMinDist = D8(i)
  iMinDist = i
 end if
end DO 

Neigh1 = neighV(1,iMinDist)
Neigh2 = neighV(2,iMinDist)
Neigh3 = neighV(3,iMinDist)

NeighR1 = Neigh1
NeighR2 = Neigh2

!kick out the 1st edge
IF (D8(Neigh1).gt.D8(Neigh2).and.D8(Neigh1).gt.D8(Neigh3)) THEN
 NeighR1 = Neigh2
 NeighR2 = Neigh3
END IF

!kick out the 2nd edge
IF (D8(Neigh2).gt.D8(Neigh1).and.D8(Neigh2).gt.D8(Neigh3)) THEN
 NeighR1 = Neigh1
 NeighR2 = Neigh3
END IF

!kick out the 3rd edge
IF (D8(Neigh3).gt.D8(Neigh1).and.D8(Neigh3).gt.D8(Neigh2)) THEN
 NeighR1 = Neigh1
 NeighR2 = Neigh2
END IF

PMIN1 = P8(:,iMinDist)
DMIN1 = D8(iMinDist)
IF (D8(NeighR1).lt.D8(NeighR2)) THEN
 PMIN2 = P8(:,NeighR1)
 DMIN2 = D8(NeighR1)
ELSE
 PMIN2 = P8(:,NeighR2)
 DMIN2 = D8(NeighR2)
END IF

CALL updatePoint(PC,DistC)

if (bP) then
 write(*,*) P8(:,iMinDist)
 write(*,*) P8(:,NeighR1)
 write(*,*) P8(:,NeighR2)
end if

CALL GetEdgeLengthAndMid(P8(:,iMinDist),P8(:,NeighR1),bEdge1,.true.)

CALL GetEdgeLengthAndMid(P8(:,iMinDist),P8(:,NeighR2),bEdge2,.true.)

CALL GetEdgeLengthAndMid(P8(:,NeighR2),P8(:,NeighR1),bEdge3,.true.)

CALL GetEdgeLengthAndMid(P8(:,iMinDist),PC,bEdge4,.true.)

if (bEdge1.or.bEdge2.or.bEdge3.or.bEdge4) THEN
!  bKeep = .TRUE.
 CALL GetEdgeLengthAndMid(PMIN1,PMIN2,bEdge,.true.)
!  if (bEdge) bKeep = .TRUE.
 if (bEdge) THEN
  CALL GetEdgeLengthAndMid(PMIN1,PMIN2,bEdgeF,.false.)
  if (bEdgeF) bKeep = .TRUE.
 END IF
 
END IF

if (bP) then
 WRITE(*,*) "XX1",PMIN1,DMIN1
 WRITE(*,*) "XX2",PMIN2,DMIN2
 pause
end if
CONTAINS
 
SUBROUTINE GetEdgeLengthAndMid(P1,P2,bX,bUpdate)
use geometry_processing, only : GetDistToSTL
REAL*8 P1(3),P2(3)
LOGICAL bX,bUpdate
REAL*8 daux,P12(3),dL

dL = DSQRT((P1(1)-P2(1))**2d0 + (P1(2)-P2(2))**2d0 + (P1(3)-P2(3))**2d0 )
P12 = 0.5d0*(P1 + P2)

CALL GetDistToSTL(p12(1),p12(2),p12(3),1,daux,.TRUE.)

IF (daux.lt.0.25d0*dL) THEN
 bX = .TRUE.
ELSE
 bX = .FALSE.
END IF

IF (bP) THEN
 write(*,*) 0.25d0*dL, daux, bX
end if

if (bUpdate) CALL updatePoint(P12,daux)

END SUBROUTINE GetEdgeLengthAndMid

SUBROUTINE  updatePoint(myP,myD)
REAL*8 myP(3),myD
real*8 dFactor

dFactor = 1d0
if (myD.gt.0d0) dFactor = 0.95d0
if (myD.lt.0d0) dFactor = 1.05d0

IF      (dFactor*myD.le.DMIN1) THEN
 PMIN2 = PMIN1
 DMIN2 = DMIN1
 PMIN1 = myP
 DMIN1 = myD
ELSE IF (dFactor*myD.le.DMIN2) THEN 
 PMIN2 = myP
 DMIN2 = myD
END IF

END SUBROUTINE  updatePoint

END SUBROUTINE TriangulationIntersectsnElement

END MODULE MeshRefDef
! ----------------------------------------------
