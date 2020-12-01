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
logical, allocatable :: bRefField(:)

integer iCombination, nOfCombinations,nAdditionalElems

nRefScheme = 3
nNonRefScheme = 1

ilev = mg_Mesh%nlmax

nel = mg_mesh%level(ilev)%nel
nvt = mg_mesh%level(ilev)%nvt

allocate(level(nvt))
allocate(bRefField(nvt))

11 continue

nAdditionalElems = 0

bRefField = .false.
level = 3
mytempCount = 0

do i=1,nel
 if (MarkerE(i).EQ.1.or.MarkerE(i).EQ.2) then
   DO j=1,8
    ivt = mg_mesh%level(ilev)%kvert(j,i)
    bRefField(ivt) = .true.
   END DO
 end if
end do


!Determine the coupling nodes
do i=1,nel

   bVert = .FALSE.
   DO j=1,8
    ivt = mg_mesh%level(ilev)%kvert(j,i)
    if (bRefField(ivt)) bVert(j) = .TRUE.
   END DO
  
   iVert = mg_mesh%level(ilev)%kvert(:,i)
   CALL RotatePatch(bVert,iVert)
   CALL DetermineTemplate(bVert,myTemplate)
   
   if (myTemplate.eq.5.or.&
       myTemplate.eq.10.or.&
       myTemplate.eq.15.or.&
       myTemplate.eq.18.or.&
       myTemplate.eq.21) then

       if (myTemplate.eq.5 ) level(iVert(4)) = min(level(iVert(4)),2)
       if (myTemplate.eq.10) level(iVert(7)) = min(level(iVert(7)),1)
       if (myTemplate.eq.15) level(iVert(7)) = min(level(iVert(7)),1)
       if (myTemplate.eq.18) level(iVert(7)) = min(level(iVert(7)),2)
       if (myTemplate.eq.18) level(iVert(8)) = min(level(iVert(8)),2)
       if (myTemplate.eq.21) level(iVert(8)) = min(level(iVert(8)),1)
    end if
   
end do

allocate(myRF(nel))

!Resolve elements with decoupling nodes

do i=1,nel

   bFound=.false.
   
   bVert = .FALSE.
   DO j=1,8
      ivt = mg_mesh%level(ilev)%kvert(j,i)
      if (level(ivt).ne.3) bVert(j) = .TRUE.
   END DO

   ii = 0
   do j=1,8
    if (bVert(j))  ii = ii + 1
   end  do
   
   iVert = mg_mesh%level(ilev)%kvert(:,i)
   CALL RotatePatch(bVert,iVert)
   CALL DetermineTemplate(bVert,myTemplate)
   
  if (myTemplate.eq.1.or.myTemplate.eq.2.or.myTemplate.eq.5.or.myTemplate.eq.10) then
   
    if (myTemplate.eq.2) then
     if (level(iVert(1)).eq.1.or.level(iVert(2)).eq.1) then
      myTemplate = 3
      if (level(iVert(1)).eq.1) THEN
       bVert(2) = .false.
       bVert(4) = .true.
       bVert(5) = .true.
      END IF
      if (level(iVert(2)).eq.1) THEN
       bVert(1) = .false.
       bVert(3) = .true.
       bVert(6) = .true.
      END IF
      CALL RotatePatch(bVert,iVert)
     end if
    END IF
    
    do j=1,8
     E(:,j) = mg_mesh%level(ilev)%dcorvg(:,iVert(j))
    end do
    
    ! assign to the patch the real refinemnet mask
    bVert = bRefField(iVert)
    

    CALL FillUpRefinedElement(myRF(i),E,ii,bvert,.true.,1)
    myRF(i)%patchID = myTemplate
    bFound=.true.
    
  else

    if (myTemplate.eq.0) bFound=.true.

  END IF
  
  if (.not.bFound) THEN
   MarkerE(i) = 2
   nAdditionalElems = nAdditionalElems + 1
  end if
    
end do

IF (nAdditionalElems.gt.0) THEN
 deallocate(myRF)
 write(*,'(A,I0,A)') 'Decoupling node combination was not possible to resolve for ',nAdditionalElems,' element!'
 GOTO 11
END IF

nOfCombinations = 0
allocate(mytemp(8,22),mytempCount(22))
mytempCount = 0

do i=1,nel

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
   
   iVert = mg_mesh%level(ilev)%kvert(:,i)
   CALL RotatePatch(bVert,iVert)
   CALL DetermineTemplate(bVert,myTemplate)
   
   do j=1,8
    E(:,j) = mg_mesh%level(ilev)%dcorvg(:,iVert(j))
   end do
    
   if (allocated(myRF(i)%kvert)) then
   
    allocate(myRF(i)%myRF(myRF(i)%nOfElem))
    
    DO iel=1,myRF(i)%nOfElem
     
     DO j=1,8
      ivt = myRF(i)%kvert(j,iel)
      iVert(j) = ivt
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
     CALL DetermineTemplate(bVert,myTemplate)
     
     do j=1,8
      E(:,j) = myRF(i)%dcoor(:,iVert(j))
     end do
     
!     if (myTemplate.eq.1) then
     if (myTemplate.eq.0.or.myTemplate.eq.1.or.myTemplate.eq.2.or.myTemplate.eq.8.or.myTemplate.eq.22) then
       CALL FillUpRefinedElement(myRF(i)%myRF(iel),E,jj,bvert,.false.,0)
    
     END IF
     
    END DO
   
   else
   
   if (myTemplate.eq.0.or.myTemplate.eq.1.or.myTemplate.eq.2.or.myTemplate.eq.8.or.myTemplate.eq.22) then
!     if (myTemplate.eq.8) then
     
     CALL FillUpRefinedElement(myRF(i),E,ii,bvert,.false.,0)
     myRF(i)%patchID = myTemplate
     bFound=.true.
    end if
   
   end if
   
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
   
end do
   

 DO iCombination = 1,nOfCombinations
   j = 0
   do i=1,8
    if (mytemp(i,iCombination)) j = j + 1
   end do
   write(*,*) j,' :: ', mytemp(:,iCombination),mytempCount(iCombination)
 END DO
 
END SUBROUTINE CreateRefinedMesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CreateRefinedMesh_Fine
integer i,j,k,l,ii,ivt,jj
real*8 E(3,8)
logical bVert(8),bVertCopy(8)
integer iVert(8)
logical bFound,bUnique
logical, allocatable :: mytemp(:,:)
integer, allocatable :: mytempCount(:)

integer iCombination, nOfCombinations,nAdditionalElems

nRefScheme = 6
nNonRefScheme = 2

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
       myTemplate.eq.9.or.&
       myTemplate.eq.10.or.&
       myTemplate.eq.11.or.&
       myTemplate.eq.15.or.&
       myTemplate.eq.18.or.&
       myTemplate.eq.21.or.&
       myTemplate.eq.22) then
    else
     
     write(*,*) bVert
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
   
   write(*,*) i,markerE(i),mytemplate
   
   if (myTemplate.eq.1.or.myTemplate.eq.2.or.myTemplate.eq.8) then
   
    iVert = mg_mesh%level(ilev)%kvert(:,i)

    CALL RotatePatch(bVert,iVert)
   
    do j=1,8
     E(:,j) = mg_mesh%level(ilev)%dcorvg(:,iVert(j))
    end do
    
    write(*,*) 'before being done',ii
    CALL FillUpRefinedElementF(myRF(i),E,ii,bvert,.false.)
    write(*,*) 'done'
    myRF(i)%patchID = myTemplate
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
    myRF(i)%patchID = myTemplate
    bFound=.true.
    
   end if

   
   if (myTemplate.eq. 5.or.&
       myTemplate.eq. 9.or.&
       myTemplate.eq.10.or.&
       myTemplate.eq.11.or.&
       myTemplate.eq.15.or.&
       myTemplate.eq.18.or.&
       myTemplate.eq.21) then
!    if ((ii.eq.3.or.ii.eq.5.or.ii.eq.6.or.ii.eq.7)) then
   
    iVert = mg_mesh%level(ilev)%kvert(:,i)

!      write(*,'(8L,8((", "),I0))') bVert,iVert
     CALL RotatePatch(bVert,iVert)
!      write(*,'(8L,8((", "),I0))') bVert,iVert
   
    do j=1,8
     E(:,j) = mg_mesh%level(ilev)%dcorvg(:,iVert(j))
    end do
    
    CALL FillUpRefinedElementF(myRF(i),E,ii,bvert,.true.)
    myRF(i)%patchID = myTemplate
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
      CALL FillUpRefinedElementF(myRF(i)%myRF(iel),E,jj+8,bvert,.false.)
      myRF(i)%myRF(iel)%selection = jj
      myRF(i)%myRF(iel)%patchID = myTemplate
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
      myRF(i)%myRF(iel)%patchID = 0 !myTemplate
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
    CALL FillUpRefinedElementF(myRF(i),E,0,bVert,.false.)
    
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

 END SUBROUTINE CreateRefinedMesh_Fine
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
SUBROUTINE FillUpRefinedElement(RF,E,iP,bV,bRecursive,iSet)
logical bRecursive
integer :: iP
type(RefinerMesh) RF
logical bV(8)
real*8 E(3,8)
real*8 :: Q8=0.125d0
real*8 DJ(8,3),XX,YY,ZZ
real*8 :: CBP(3,64),xi1,xi2,xi3,DJAC(3,3),DETJ,DKNPR(8)
integer i,j,k,ii
real*8 :: dPhi = 0.0d0,dVal
character*256 cF
integer, allocatable :: kvert(:,:)

if (iSet.eq.1) then
 cF = cTemplatesX(myTemplate)
else
 cF = cTemplates(myTemplate)
end if


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

DO i=1,8
 IF (bV(i)) THEN
  dknpr(i) = 1d0
 else
  dknpr(i) = 0d0
 END IF
END DO

! !  if (myTemplate.eq.18) then
! !   RF%knpr( 6) = 1
! !   RF%knpr(11) = 1
! !   RF%knpr (2) = 1
! !   RF%knpr( 4) = 1
! !   RF%knpr( 9) = 1
! !   RF%knpr(10) = 1
! !  END IF
! ! !  
! !  if (myTemplate.eq.21) then
! !   RF%knpr(1) = 1
! !   RF%knpr(7) = 1
! !   RF%knpr(3) = 1
! !   RF%knpr(5) = 1
! !   RF%knpr(12) = 1
! !   RF%knpr(10) = 1
! !   RF%knpr(14) = 1
! !  END IF
  
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
   
   if (bRecursive) then
    CALL INT_E011(dVal)
    IF (dVAl.gt.0.99d0) THEN
     RF%knpr(ii)  = 1d0
    ELSE
     RF%knpr(ii)  = 0d0
    END IF
!     write(*,'(8E12.4,A,E12.4,I8)') dknpr(:),' : ',dVal,RF%knpr(ii)
!     pause
   end if
   
END DO

 CONTAINS
 
  SUBROUTINE INT_E011(dU)
  implicit none
  REAL*8 dU
  REAL*8 :: Q8 = 0.125d0
  REAL*8 DH(8)
  integer idofe
 
  DH(1)=Q8*(1D0-Xi1)*(1D0-Xi2)*(1D0-Xi3)
  DH(2)=Q8*(1D0+Xi1)*(1D0-Xi2)*(1D0-Xi3)
  DH(3)=Q8*(1D0+Xi1)*(1D0+Xi2)*(1D0-Xi3)
  DH(4)=Q8*(1D0-Xi1)*(1D0+Xi2)*(1D0-Xi3)
  DH(5)=Q8*(1D0-Xi1)*(1D0-Xi2)*(1D0+Xi3)
  DH(6)=Q8*(1D0+Xi1)*(1D0-Xi2)*(1D0+Xi3)
  DH(7)=Q8*(1D0+Xi1)*(1D0+Xi2)*(1D0+Xi3)
  DH(8)=Q8*(1D0-Xi1)*(1D0+Xi2)*(1D0+Xi3)

  dU=0D0
  DO idofe=1,8
     dU = dU + DKNPR(idofe)*DH(idofe)
  END DO

  RETURN
  
  END SUBROUTINE INT_E011

END SUBROUTINE FillUpRefinedElement
! ----------------------------------------------
SUBROUTINE FillUpRefinedElementF(RF,E,iP,bV,bRecursive)
logical bRecursive
integer :: iP
type(RefinerMesh) RF
logical bV(8)
real*8 E(3,8)
real*8 :: Q8=0.125d0
real*8 DJ(8,3),XX,YY,ZZ
real*8 :: xi1,xi2,xi3,DJAC(3,3),DETJ
real*8 , allocatable ::cBP(:,:)
integer i,j,k,ii
real*8 :: dPhi = 0.0d0
character*256 cF
integer, allocatable :: kvert(:,:)

 cF = cTemplatesF(myTemplate)

open(file=TRIM(ADJUSTL(cF)),unit=1)

read(1,*)
read(1,*) 
read(1,*) RF%nOfElem, RF%nOfVert
read(1,*)

allocate(cBP(3,RF%nOfVert))
allocate(RF%kvert(8,RF%nOfElem))
allocate(RF%dcoor(3,RF%nOfVert))
allocate(RF%knpr(RF%nOfVert))
RF%knpr = 0

if (iP.eq.4) write(*,*) 'here -- CRIT',RF%nOfElem,RF%nOfVert,ADJUSTL(TRIM(cF))
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
 
 if (myTemplate.eq.9) then
   RF%knpr(49) = 1
   RF%knpr(29) = 1
   RF%knpr(70) = 1
   RF%knpr(65) = 1
   RF%knpr(72) = 1
   RF%knpr(59) = 1
   RF%knpr(61) = 1
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
 
 if (myTemplate.eq.11) then
   RF%knpr(30) = 1
   RF%knpr(37) = 1
   RF%knpr(72) = 1
   RF%knpr(70) = 1
   RF%knpr(65) = 1
   RF%knpr(54) = 1
   RF%knpr(68) = 1
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

deallocate(cBP)

END SUBROUTINE FillUpRefinedElementF

END Module MeshRefRefine