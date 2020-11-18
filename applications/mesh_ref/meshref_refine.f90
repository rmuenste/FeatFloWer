Module MeshRefRefine

use MeshRefDef

 contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CreateRefinedMesh
integer i,j,k,l,ii,ivt,jj
real*8 E(3,8)
logical bVert(8),bVertCopy(8)
integer iVert(8)
logical bFound,bUnique
logical, allocatable :: mytemp(:,:)
integer, allocatable :: mytempCount(:)

integer iCombination, nOfCombinations,nAdditionalElems

ilev = mg_Mesh%nlmax

nel = mg_mesh%level(ilev)%nel
nvt = mg_mesh%level(ilev)%nvt

nOfCombinations = 0
allocate(mytemp(8,22),mytempCount(22))
mytempCount = 0

11 continue

nAdditionalElems = 0
do i=1,nel
 if (markerE(i).eq.3) then
 
   bFound = .false.
   bVert = .FALSE.
   DO j=1,8
    ivt = mg_mesh%level(ilev)%kvert(j,i)
    do k=1,mg_mesh%level(ilev)%nvel
     l = mg_mesh%level(ilev)%kvel(k,ivt)
     if (l.ne.0) then
      if (MarkerE(l).eq.1.or.MarkerE(l).eq.2) then
       bVert(j) = .TRUE.
      end if
     end if
    end do
   END DO
  
   CALL RotatePatch(bVert,iVert)
   CALL DetermineTemplate(bVert,myTemplate)
   
   if (myTemplate.eq.1.or.&
       myTemplate.eq.2.or.&
       myTemplate.eq.5.or.&
       myTemplate.eq.8.or.&
       myTemplate.eq.10.or.&
       myTemplate.eq.15.or.&
       myTemplate.eq.18.or.&
       myTemplate.eq.21.or.&
       myTemplate.eq.22) then
    else
     markerE(i)= 1
     
     do j=1,8
      ivt = mg_mesh%level(ilev)%kvert(j,i)
      do k=1,mg_mesh%level(ilev)%nvel
       l = mg_mesh%level(ilev)%kvel(k,ivt)
       if (l.ne.0) then
        if (MarkerE(l).ne.1) MarkerE(l) = 3
       end if
      end do
     end do
     
     nAdditionalElems = nAdditionalElems + 1
     
    end if
    
   
 end if 
end do

write(*,*) nAdditionalElems, " additional element has been set to full refinement"

if (nAdditionalElems.gt.0) goto 11

allocate(myRF(nel))

do i=1,nel

 if (markerE(i).eq.0) then
 
   myRF(i)%nOfElem = nNonRefScheme**3
   myRF(i)%nOfVert = (nNonRefScheme+1)**3
   allocate(myRF(i)%kvert(8,myRF(i)%nOfElem))
   allocate(myRF(i)%dcoor(3,myRF(i)%nOfVert))
   allocate(myRF(i)%knpr(myRF(i)%nOfVert))
   myRF(i)%kvert=0
   do j=1,8
    E(:,j) = mg_mesh%level(ilev)%dcorvg(:,mg_mesh%level(ilev)%kvert(j,i))
   end do
   
   Call FillUpRefined(myRF(i),E,nNonRefScheme+1)
   myRF(i)%patchID = -1
 end if

 if (markerE(i).eq.1.or.markerE(i).eq.2) then
 
   myRF(i)%nOfElem = nRefScheme**3
   myRF(i)%nOfVert = (nRefScheme+1)**3
   allocate(myRF(i)%kvert(8,myRF(i)%nOfElem))
   allocate(myRF(i)%dcoor(3,myRF(i)%nOfVert))
   allocate(myRF(i)%knpr(myRF(i)%nOfVert))
   myRF(i)%knpr=0
   do j=1,8
    E(:,j) = mg_mesh%level(ilev)%dcorvg(:,mg_mesh%level(ilev)%kvert(j,i))
   end do
   
   Call FillUpRefined(myRF(i),E,nRefScheme+1)
   myRF(i)%patchID = 0

 end if

 if (markerE(i).eq.3) then
 
   bFound = .false.
   bVert = .FALSE.
   DO j=1,8
    ivt = mg_mesh%level(ilev)%kvert(j,i)
    do k=1,mg_mesh%level(ilev)%nvel
     l = mg_mesh%level(ilev)%kvel(k,ivt)
     if (l.ne.0) then
      if (MarkerE(l).eq.1.or.MarkerE(l).eq.2) then
       bVert(j) = .TRUE.
      end if
     end if
    end do
   END DO

   ii = 0
   do j=1,8
    if (bVert(j))  ii = ii + 1
   end  do
   
   bVertCopy = bVert
   CALL RotatePatch(bVert,iVert)
   CALL DetermineTemplate(bVert,myTemplate)
   bVert = bVertCopy
   
   if (myTemplate.eq.1.or.myTemplate.eq.2.or.myTemplate.eq.8) then
   
    iVert = mg_mesh%level(ilev)%kvert(:,i)

    CALL RotatePatch(bVert,iVert)
   
    do j=1,8
     E(:,j) = mg_mesh%level(ilev)%dcorvg(:,iVert(j))
    end do
    
    CALL FillUpRefinedElement(myRF(i),E,ii,bvert,.false.)
    myRF(i)%patchID = ii
    bFound=.true.
!     write(*,*) myTemplate,ii
    
   end if

   if (ii.eq.8) then
   
    myRF(i)%nOfElem = nRefScheme**3
    myRF(i)%nOfVert = (nRefScheme+1)**3
    allocate(myRF(i)%kvert(8,myRF(i)%nOfElem))
    allocate(myRF(i)%dcoor(3,myRF(i)%nOfVert))
    allocate(myRF(i)%knpr(myRF(i)%nOfVert))
    myRF(i)%knpr=0
    do j=1,8
     E(:,j) = mg_mesh%level(ilev)%dcorvg(:,mg_mesh%level(ilev)%kvert(j,i))
    end do
    
    Call FillUpRefined(myRF(i),E,nRefScheme+1)
    myRF(i)%patchID = ii
    bFound=.true.
    
   end if

   
   if (myTemplate.eq.5.or.myTemplate.eq.10.or.myTemplate.eq.15.or.myTemplate.eq.18.or.myTemplate.eq.21) then
!    if ((ii.eq.3.or.ii.eq.5.or.ii.eq.6.or.ii.eq.7)) then
   
    iVert = mg_mesh%level(ilev)%kvert(:,i)

!      write(*,'(8L,8((", "),I0))') bVert,iVert
     CALL RotatePatch(bVert,iVert)
!      write(*,'(8L,8((", "),I0))') bVert,iVert
   
    do j=1,8
     E(:,j) = mg_mesh%level(ilev)%dcorvg(:,iVert(j))
    end do
    
    CALL FillUpRefinedElement(myRF(i),E,ii,bvert,.true.)
    myRF(i)%patchID = ii
    bFound=.true.
    
    !!!!!!!!!!! RECURSIVE REFINEMENT !!!!!!!!!!!!!!!
    allocate(myRF(i)%myRF(myRF(i)%nOfElem))
    
    DO iel=1,myRF(i)%nOfElem
     
     DO j=1,8
      ivt = myRF(i)%kvert(j,iel)
      iVert(j) = ivt
!       write(*,*) ivt
      if (myRF(i)%knpr(ivt).eq.1) then
       bVert(j)=.true.
      else
       bVert(j)=.false.
      endif
     END DO
     
     jj = 0
     do j=1,8
      if (bVert(j))  jj = jj + 1
     end  do
     
     CALL RotatePatch(bVert,iVert)
!      write(*,'(4(I0,(": ")),8L,8((", "),I0))') jj,iel,myRF(i)%nOfElem,ii,bVert,iVert
     
     do j=1,8
      E(:,j) = myRF(i)%dcoor(:,iVert(j))
     end do
     
     if (jj.gt.0) then
      myTemplate=22+jj
      CALL FillUpRefinedElement(myRF(i)%myRF(iel),E,jj+8,bvert,.false.)
      myRF(i)%myRF(iel)%selection = jj
      myRF(i)%myRF(iel)%patchID = ii
!       write(*,*) 'jj=',jj, bVert,cPatches(jj+8),cSPatches(jj+8)
      
     else
!       write(*,*) 'jj=',jj
      myRF(i)%myRF(iel)%nOfElem = 1
      myRF(i)%myRF(iel)%nOfVert = 8
      allocate(myRF(i)%myRF(iel)%kvert(8,myRF(i)%myRF(iel)%nOfElem))
      allocate(myRF(i)%myRF(iel)%dcoor(3,myRF(i)%myRF(iel)%nOfVert))
      allocate(myRF(i)%myRF(iel)%knpr(myRF(i)%myRF(iel)%nOfVert))
      myRF(i)%myRF(iel)%kvert=0
      
      Call FillUpRefined(myRF(i)%myRF(iel),E,1+1)
      myRF(i)%myRF(iel)%patchID = ii
      myRF(i)%myRF(iel)%selection = jj
     
     end if
    end do
    
   end if
   
   if (.not.bFound) then
     
    iVert = mg_mesh%level(ilev)%kvert(:,i)
    CALL RotatePatch(bVert,iVert)
    
   
    do j=1,8
     E(:,j) = mg_mesh%level(ilev)%dcorvg(:,iVert(j))
    end do
    CALL FillUpRefinedElement(myRF(i),E,0,bVert,.false.)
    
    myRF(i)%patchID = -2

   end if

   bVert = bVertCopy
   iVert = 0
   CALL RotatePatch(bVert,iVert)

   bUnique = .true.
   DO iCombination = 1,nOfCombinations
    CALL ComparePatch(mytemp(:,iCombination),bVert,bUnique)
    if (.not.bUnique) exit
   END DO
   
   if (bUnique) then
    nOfCombinations = nOfCombinations+1
    mytemp(:,nOfCombinations) = bVert
    mytempCount(iCombination) = 1
   else
    mytempCount(iCombination) = mytempCount(iCombination) + 1
   end if
   
 end if

end do

 DO iCombination = 1,nOfCombinations
   j = 0
   do i=1,8
    if (mytemp(i,iCombination)) j = j + 1
   end do
   
   write(*,*) j,' :: ', mytemp(:,iCombination),mytempCount(iCombination)
   
 END DO

 END SUBROUTINE CreateRefinedMesh
! ----------------------------------------------
SUBROUTINE FillUpRefined(RF,E,nCBP)
integer nCBP
type(RefinerMesh) RF
real*8 E(3,8)
real*8 :: Q8=0.125d0
real*8 DJ(8,3),XX,YY,ZZ
real*8 :: xi1,xi2,xi3,DJAC(3,3),DETJ
integer i,j,k,ii

DJ(1,1)=( E(1,1)+E(1,2)+E(1,3)+E(1,4)+E(1,5)+E(1,6)+E(1,7)+E(1,8))*Q8
DJ(1,2)=( E(2,1)+E(2,2)+E(2,3)+E(2,4)+E(2,5)+E(2,6)+E(2,7)+E(2,8))*Q8
DJ(1,3)=( E(3,1)+E(3,2)+E(3,3)+E(3,4)+E(3,5)+E(3,6)+E(3,7)+E(3,8))*Q8
DJ(2,1)=(-E(1,1)+E(1,2)+E(1,3)-E(1,4)-E(1,5)+E(1,6)+E(1,7)-E(1,8))*Q8
DJ(2,2)=(-E(2,1)+E(2,2)+E(2,3)-E(2,4)-E(2,5)+E(2,6)+E(2,7)-E(2,8))*Q8
DJ(2,3)=(-E(3,1)+E(3,2)+E(3,3)-E(3,4)-E(3,5)+E(3,6)+E(3,7)-E(3,8))*Q8
DJ(3,1)=(-E(1,1)-E(1,2)+E(1,3)+E(1,4)-E(1,5)-E(1,6)+E(1,7)+E(1,8))*Q8
DJ(3,2)=(-E(2,1)-E(2,2)+E(2,3)+E(2,4)-E(2,5)-E(2,6)+E(2,7)+E(2,8))*Q8
DJ(3,3)=(-E(3,1)-E(3,2)+E(3,3)+E(3,4)-E(3,5)-E(3,6)+E(3,7)+E(3,8))*Q8
DJ(4,1)=(-E(1,1)-E(1,2)-E(1,3)-E(1,4)+E(1,5)+E(1,6)+E(1,7)+E(1,8))*Q8
DJ(4,2)=(-E(2,1)-E(2,2)-E(2,3)-E(2,4)+E(2,5)+E(2,6)+E(2,7)+E(2,8))*Q8
DJ(4,3)=(-E(3,1)-E(3,2)-E(3,3)-E(3,4)+E(3,5)+E(3,6)+E(3,7)+E(3,8))*Q8
DJ(5,1)=( E(1,1)-E(1,2)+E(1,3)-E(1,4)+E(1,5)-E(1,6)+E(1,7)-E(1,8))*Q8
DJ(5,2)=( E(2,1)-E(2,2)+E(2,3)-E(2,4)+E(2,5)-E(2,6)+E(2,7)-E(2,8))*Q8
DJ(5,3)=( E(3,1)-E(3,2)+E(3,3)-E(3,4)+E(3,5)-E(3,6)+E(3,7)-E(3,8))*Q8
DJ(6,1)=( E(1,1)-E(1,2)-E(1,3)+E(1,4)-E(1,5)+E(1,6)+E(1,7)-E(1,8))*Q8
DJ(6,2)=( E(2,1)-E(2,2)-E(2,3)+E(2,4)-E(2,5)+E(2,6)+E(2,7)-E(2,8))*Q8
DJ(6,3)=( E(3,1)-E(3,2)-E(3,3)+E(3,4)-E(3,5)+E(3,6)+E(3,7)-E(3,8))*Q8
DJ(7,1)=( E(1,1)+E(1,2)-E(1,3)-E(1,4)-E(1,5)-E(1,6)+E(1,7)+E(1,8))*Q8
DJ(7,2)=( E(2,1)+E(2,2)-E(2,3)-E(2,4)-E(2,5)-E(2,6)+E(2,7)+E(2,8))*Q8
DJ(7,3)=( E(3,1)+E(3,2)-E(3,3)-E(3,4)-E(3,5)-E(3,6)+E(3,7)+E(3,8))*Q8
DJ(8,1)=(-E(1,1)+E(1,2)-E(1,3)+E(1,4)+E(1,5)-E(1,6)+E(1,7)-E(1,8))*Q8
DJ(8,2)=(-E(2,1)+E(2,2)-E(2,3)+E(2,4)+E(2,5)-E(2,6)+E(2,7)-E(2,8))*Q8
DJ(8,3)=(-E(3,1)+E(3,2)-E(3,3)+E(3,4)+E(3,5)-E(3,6)+E(3,7)-E(3,8))*Q8

ii = 0
do i=1,nCBP
 do j=1,nCBP
  do k=1,nCBP
   
   xi1 = -1d0 + 2d0*DBLE(i-1)/DBLE(nCBP-1)
   xi2 = -1d0 + 2d0*DBLE(j-1)/DBLE(nCBP-1)
   xi3 = -1d0 + 2d0*DBLE(k-1)/DBLE(nCBP-1)
   ii = ii + 1

   DJAC(1,1)=DJ(2,1) + DJ(5,1)*XI2 + DJ(6,1)*XI3 + DJ(8,1)*XI2*XI3
   DJAC(1,2)=DJ(3,1) + DJ(5,1)*XI1 + DJ(7,1)*XI3 + DJ(8,1)*XI1*XI3
   DJAC(1,3)=DJ(4,1) + DJ(6,1)*XI1 + DJ(7,1)*XI2 + DJ(8,1)*XI1*XI2
   DJAC(2,1)=DJ(2,2) + DJ(5,2)*XI2 + DJ(6,2)*XI3 + DJ(8,2)*XI2*XI3
   DJAC(2,2)=DJ(3,2) + DJ(5,2)*XI1 + DJ(7,2)*XI3 + DJ(8,2)*XI1*XI3
   DJAC(2,3)=DJ(4,2) + DJ(6,2)*XI1 + DJ(7,2)*XI2 + DJ(8,2)*XI1*XI2
   DJAC(3,1)=DJ(2,3) + DJ(5,3)*XI2 + DJ(6,3)*XI3 + DJ(8,3)*XI2*XI3
   DJAC(3,2)=DJ(3,3) + DJ(5,3)*XI1 + DJ(7,3)*XI3 + DJ(8,3)*XI1*XI3
   DJAC(3,3)=DJ(4,3) + DJ(6,3)*XI1 + DJ(7,3)*XI2 + DJ(8,3)*XI1*XI2
   
   DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3)) &
        -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3)) &
        +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
     
   XX=DJ(1,1) + DJAC(1,1)*XI1 + DJ(3,1)*XI2 + DJ(4,1)*XI3 + DJ(7,1)*XI2*XI3
   YY=DJ(1,2) + DJ(2,2)*XI1 + DJAC(2,2)*XI2 + DJ(4,2)*XI3 + DJ(6,2)*XI1*XI3
   ZZ=DJ(1,3) + DJ(2,3)*XI1 + DJ(3,3)*XI2 + DJAC(3,3)*XI3 + DJ(5,3)*XI1*XI2
   
   RF%dcoor(:,ii) = [xx,yy,zz]
     
!    write(*,*) RF%dcoor(:,ii)
   
   end do
 end do
end do

ii = 0
do i=1,nCBP-1
 do j=1,nCBP-1
  do k=1,nCBP-1
   
   ii = ii + 1
   RF%kvert(:,ii) = [(k-1)*nCBP*nCBP+(j-1)*nCBP + i,&
                     (k-1)*nCBP*nCBP+(j-1)*nCBP + i + 1 ,&
                     (k-1)*nCBP*nCBP+(j)*nCBP + i + 1 ,&
                     (k-1)*nCBP*nCBP+(j)*nCBP + i ,&
                     (k)*nCBP*nCBP+(j-1)*nCBP + i ,&
                     (k)*nCBP*nCBP+(j-1)*nCBP + i + 1 ,&
                     (k)*nCBP*nCBP+(j)*nCBP + i + 1 ,&
                     (k)*nCBP*nCBP+(j)*nCBP + i ]
   
!    write(*,*) RF%kvert(:,ii)
   
   end do
 end do
end do

! pause

END SUBROUTINE FillUpRefined
! ----------------------------------------------
SUBROUTINE FillUpRefinedElement(RF,E,iP,bV,bRecursive)
logical bRecursive
integer :: iP
type(RefinerMesh) RF
logical bV(8)
real*8 E(3,8)
real*8 :: Q8=0.125d0
real*8 DJ(8,3),XX,YY,ZZ
real*8 :: CBP(3,64),xi1,xi2,xi3,DJAC(3,3),DETJ
integer i,j,k,ii
real*8 :: dPhi = 0.0d0
character*256 cF
integer, allocatable :: kvert(:,:)

 cF = cTemplates(myTemplate)

open(file=TRIM(ADJUSTL(cF)),unit=1)

read(1,*)
read(1,*) 
read(1,*) RF%nOfElem, RF%nOfVert
read(1,*)

allocate(RF%kvert(8,RF%nOfElem))
allocate(RF%dcoor(3,RF%nOfVert))
allocate(RF%knpr(RF%nOfVert))
RF%knpr = 0

DO i=1,RF%nOfVert
 read(1,*) CBP(:,i)
end do

read(1,*)
 
DO i=1,RF%nOfElem
 read(1,*) RF%kvert(:,i)
end do
 
close(1)

if (bRecursive) then

 if (myTemplate.eq.5) then
   RF%knpr(20) = 1
   RF%knpr(27) = 1
   RF%knpr(33) = 1
   RF%knpr(22) = 1
   RF%knpr(28) = 1
 end if

 if (myTemplate.eq.10) then
   RF%knpr(31) = 1
   RF%knpr(71) = 1
   RF%knpr(47) = 1
   RF%knpr(59) = 1
   RF%knpr(72) = 1
   RF%knpr(56) = 1
   RF%knpr(87) = 1
 end if
 
 if (myTemplate.eq.15) then
   RF%knpr(5)  = 1
   RF%knpr(12) = 1
   RF%knpr(17) = 1
   RF%knpr(19) = 1
   RF%knpr(29) = 1
   RF%knpr(49) = 1
   RF%knpr(54) = 1
   RF%knpr(68) = 1
   RF%knpr(70) = 1
   RF%knpr(65) = 1
   RF%knpr(72) = 1
 end if
 
 if (myTemplate.eq.18) then
  RF%knpr(8 ) = 1
  RF%knpr(10) = 1
  RF%knpr(13) = 1
  RF%knpr(21) = 1
  RF%knpr(25) = 1
  RF%knpr(32) = 1
  RF%knpr(38) = 1
  RF%knpr(42) = 1
  RF%knpr(49) = 1
  RF%knpr(4 ) = 1
  RF%knpr(9 ) = 1
  RF%knpr(19) = 1
  RF%knpr(29) = 1
  RF%knpr(36) = 1
  RF%knpr(46) = 1
 END IF
 
 if (myTemplate.eq.21) then
  RF%knpr(87 ) = 1
  RF%knpr(88 ) = 1
  RF%knpr(89 ) = 1
  RF%knpr(92 ) = 1
  RF%knpr(93 ) = 1
  RF%knpr(97 ) = 1
  RF%knpr(100) = 1
  RF%knpr(103) = 1
  RF%knpr(104) = 1
  RF%knpr(65 ) = 1
  RF%knpr(67 ) = 1
  RF%knpr(70 ) = 1
  RF%knpr(71 ) = 1
  RF%knpr(75 ) = 1
  RF%knpr(77 ) = 1
  RF%knpr(72 ) = 1
  RF%knpr(73 ) = 1
  RF%knpr(75 ) = 1
  RF%knpr(77 ) = 1
  RF%knpr(79 ) = 1
  RF%knpr(85 ) = 1
 END IF
  
end if

DJ(1,1)=( E(1,1)+E(1,2)+E(1,3)+E(1,4)+E(1,5)+E(1,6)+E(1,7)+E(1,8))*Q8
DJ(1,2)=( E(2,1)+E(2,2)+E(2,3)+E(2,4)+E(2,5)+E(2,6)+E(2,7)+E(2,8))*Q8
DJ(1,3)=( E(3,1)+E(3,2)+E(3,3)+E(3,4)+E(3,5)+E(3,6)+E(3,7)+E(3,8))*Q8
DJ(2,1)=(-E(1,1)+E(1,2)+E(1,3)-E(1,4)-E(1,5)+E(1,6)+E(1,7)-E(1,8))*Q8
DJ(2,2)=(-E(2,1)+E(2,2)+E(2,3)-E(2,4)-E(2,5)+E(2,6)+E(2,7)-E(2,8))*Q8
DJ(2,3)=(-E(3,1)+E(3,2)+E(3,3)-E(3,4)-E(3,5)+E(3,6)+E(3,7)-E(3,8))*Q8
DJ(3,1)=(-E(1,1)-E(1,2)+E(1,3)+E(1,4)-E(1,5)-E(1,6)+E(1,7)+E(1,8))*Q8
DJ(3,2)=(-E(2,1)-E(2,2)+E(2,3)+E(2,4)-E(2,5)-E(2,6)+E(2,7)+E(2,8))*Q8
DJ(3,3)=(-E(3,1)-E(3,2)+E(3,3)+E(3,4)-E(3,5)-E(3,6)+E(3,7)+E(3,8))*Q8
DJ(4,1)=(-E(1,1)-E(1,2)-E(1,3)-E(1,4)+E(1,5)+E(1,6)+E(1,7)+E(1,8))*Q8
DJ(4,2)=(-E(2,1)-E(2,2)-E(2,3)-E(2,4)+E(2,5)+E(2,6)+E(2,7)+E(2,8))*Q8
DJ(4,3)=(-E(3,1)-E(3,2)-E(3,3)-E(3,4)+E(3,5)+E(3,6)+E(3,7)+E(3,8))*Q8
DJ(5,1)=( E(1,1)-E(1,2)+E(1,3)-E(1,4)+E(1,5)-E(1,6)+E(1,7)-E(1,8))*Q8
DJ(5,2)=( E(2,1)-E(2,2)+E(2,3)-E(2,4)+E(2,5)-E(2,6)+E(2,7)-E(2,8))*Q8
DJ(5,3)=( E(3,1)-E(3,2)+E(3,3)-E(3,4)+E(3,5)-E(3,6)+E(3,7)-E(3,8))*Q8
DJ(6,1)=( E(1,1)-E(1,2)-E(1,3)+E(1,4)-E(1,5)+E(1,6)+E(1,7)-E(1,8))*Q8
DJ(6,2)=( E(2,1)-E(2,2)-E(2,3)+E(2,4)-E(2,5)+E(2,6)+E(2,7)-E(2,8))*Q8
DJ(6,3)=( E(3,1)-E(3,2)-E(3,3)+E(3,4)-E(3,5)+E(3,6)+E(3,7)-E(3,8))*Q8
DJ(7,1)=( E(1,1)+E(1,2)-E(1,3)-E(1,4)-E(1,5)-E(1,6)+E(1,7)+E(1,8))*Q8
DJ(7,2)=( E(2,1)+E(2,2)-E(2,3)-E(2,4)-E(2,5)-E(2,6)+E(2,7)+E(2,8))*Q8
DJ(7,3)=( E(3,1)+E(3,2)-E(3,3)-E(3,4)-E(3,5)-E(3,6)+E(3,7)+E(3,8))*Q8
DJ(8,1)=(-E(1,1)+E(1,2)-E(1,3)+E(1,4)+E(1,5)-E(1,6)+E(1,7)-E(1,8))*Q8
DJ(8,2)=(-E(2,1)+E(2,2)-E(2,3)+E(2,4)+E(2,5)-E(2,6)+E(2,7)-E(2,8))*Q8
DJ(8,3)=(-E(3,1)+E(3,2)-E(3,3)+E(3,4)+E(3,5)-E(3,6)+E(3,7)-E(3,8))*Q8

DO ii=1,RF%nOfVert

   xi1 = CBP(1,ii)
   xi2 = CBP(2,ii)
   xi3 = CBP(3,ii)

   DJAC(1,1)=DJ(2,1) + DJ(5,1)*XI2 + DJ(6,1)*XI3 + DJ(8,1)*XI2*XI3
   DJAC(1,2)=DJ(3,1) + DJ(5,1)*XI1 + DJ(7,1)*XI3 + DJ(8,1)*XI1*XI3
   DJAC(1,3)=DJ(4,1) + DJ(6,1)*XI1 + DJ(7,1)*XI2 + DJ(8,1)*XI1*XI2
   DJAC(2,1)=DJ(2,2) + DJ(5,2)*XI2 + DJ(6,2)*XI3 + DJ(8,2)*XI2*XI3
   DJAC(2,2)=DJ(3,2) + DJ(5,2)*XI1 + DJ(7,2)*XI3 + DJ(8,2)*XI1*XI3
   DJAC(2,3)=DJ(4,2) + DJ(6,2)*XI1 + DJ(7,2)*XI2 + DJ(8,2)*XI1*XI2
   DJAC(3,1)=DJ(2,3) + DJ(5,3)*XI2 + DJ(6,3)*XI3 + DJ(8,3)*XI2*XI3
   DJAC(3,2)=DJ(3,3) + DJ(5,3)*XI1 + DJ(7,3)*XI3 + DJ(8,3)*XI1*XI3
   DJAC(3,3)=DJ(4,3) + DJ(6,3)*XI1 + DJ(7,3)*XI2 + DJ(8,3)*XI1*XI2
   
   DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3)) &
        -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3)) &
        +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
     
   XX=DJ(1,1) + DJAC(1,1)*XI1 + DJ(3,1)*XI2 + DJ(4,1)*XI3 + DJ(7,1)*XI2*XI3
   YY=DJ(1,2) + DJ(2,2)*XI1 + DJAC(2,2)*XI2 + DJ(4,2)*XI3 + DJ(6,2)*XI1*XI3
   ZZ=DJ(1,3) + DJ(2,3)*XI1 + DJ(3,3)*XI2 + DJAC(3,3)*XI3 + DJ(5,3)*XI1*XI2
   
   RF%dcoor(:,ii) = [xx,yy,zz]
   
END DO


END SUBROUTINE FillUpRefinedElement

END Module MeshRefRefine