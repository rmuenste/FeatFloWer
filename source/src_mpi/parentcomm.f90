  ! ----------------------------------------------
  ! ----------------------------------------------
  ! ----------------------------------------------
  SUBROUTINE PARENTCOMM(NAT,NEL,NVT,DCORVG,DCORAG,KAREA,KVERT) !ok
  USE PP3D_MPI
  
    implicit none
    REAL*8  DCORAG(3,*),DCORVG(3,*)
    INTEGER KAREA(6,*),KVERT(8,*)
    INTEGER NAT,NEL,NVT
    INTEGER I,J,JI,K,iaux
    REAL*8,ALLOCATABLE ::  pCOORDINATES(:,:),dCOORDINATES(:,:)
    REAL*8 DIST
    real*4 time0,time1
    CHARACTER*256 :: cFileOut

    INTEGER pNEL,pNAT,pNVT,pID,pJD,n,iFace,jFace
    INTEGER ,DIMENSION(:,:), ALLOCATABLE :: ParFind
    INTEGER ,DIMENSION(:,:), ALLOCATABLE :: myFACELINK,myFACEPINK
    INTEGER ,DIMENSION(:), ALLOCATABLE :: nFACELISTS
    INTEGER :: NodeTab(subnodes,subnodes),iCount
    CHARACTER*9 ccgcc
    integer, allocatable :: sendcounts(:),displs(:),gathered_data(:)

    TYPE TVector
      INTEGER :: i,Num
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: Mids,Aux
    END TYPE TVector

    TYPE TStructure
      INTEGER :: NeighNum
      TYPE(TVector), DIMENSION(:), ALLOCATABLE :: Face
    END TYPE TStructure
    TYPE(TStructure), DIMENSION(:), ALLOCATABLE :: NeighSt

    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

    !--------------------------------------------------------------
    !CREATING MAPPING STRUCTURE FOR MASTER --> ASSISTANT //ELEMENTS
    !--------------------------------------------------------------
    
    call ztime(time0)
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    If (myid.eq.1) write(*,*) "ParentComm start!"
    
    IF (myid.eq.MASTER) THEN
      coarse%pElem = NEL
      coarse%pFace = NAT
      coarse%pVert = NVT
      ALLOCATE (coarse%pElemLink(subnodes,coarse%pElem))
      ALLOCATE (coarse%pNEL(subnodes))
      ALLOCATE (coarse%pNVT(subnodes))
      ALLOCATE (coarse%pFaceLink(subnodes,coarse%pFace))
      ALLOCATE (coarse%pVERTLink(subnodes,coarse%pVert))
      ALLOCATE (coarse%pDX(coarse%pFace))
      ALLOCATE(coarse%myELEMLINK(NEL))
      ALLOCATE(coarse%myVERTLINK(NVT))
    ELSE
      coarse%pFace = NAT
      ALLOCATE (coarse%DX(coarse%pFace))
      ALLOCATE(coarse%myELEMLINK(NEL))
      ALLOCATE(coarse%myVERTLINK(NVT))
    END IF

    !
    ! --------------------------------------- update ---------------------------------------------------------
    !

    IF (myid.eq.master) THEN
     pNVT=NVT
     ALLOCATE (dCOORDINATES(3,pNVT))
     DO I=1,pNVT
      dCOORDINATES(1,I)=DCORVG(1,I)
      dCOORDINATES(2,I)=DCORVG(2,I)
      dCOORDINATES(3,I)=DCORVG(3,I)
     END DO
    END IF
    
    call MPI_Bcast(pNVT, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    IF (myid.ne.master) THEN
     ALLOCATE (dCOORDINATES(3,pNVT))
    END IF
    
    call MPI_Bcast(dCOORDINATES, 3*pNVT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    IF (myid.ne.master) THEN
      CALL InitOctTree(dCOORDINATES,pNVT)
      
       DO I=1,NVT
       
        CALL FindInOctTree(dCOORDINATES,pNVT,DCORVG(:,I),J,DIST)
        
        IF (J.lt.0) then
         WRITE(*,*) I,"PROBLEM of vert assignement ..."
        end if
        IF (DIST.LT.DEpsPrec) THEN 
         coarse%myVERTLINK(I)=J
        END IF
       
       END DO
    
      CALL FreeOctTree()
    END IF
    
    allocate(sendcounts(0:numnodes),displs(0:numnodes+1))
    
    call MPI_allgather(NVT, 1, MPI_INTEGER, sendcounts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    sendcounts(0) = 0

    displs = 0
    do i = 2, numnodes+1
      displs(i) = displs(i-1) + sendcounts(i-1)
    end do
    
    if (myid.eq.master) then
      allocate(gathered_data(displs(numnodes+1)))
      n=0
    else 
     n = nvt
    endif
    
    call MPI_Gatherv(coarse%myVERTLINK, n, MPI_INTEGER, &
                   gathered_data, sendcounts, displs, &
                   MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
                   
    IF (myid.eq.0) THEN
      coarse%pVERTLINK = 0
      DO pID=1,subnodes
       j = 0
       DO i=displs(pID)+1, displs(pID+1)
        j = j + 1
        coarse%pVERTLINK(pID,J)=gathered_data(i)
       END DO
      end do
    END IF

!     IF (myid.eq.0) THEN
!      DO pID=1,subnodes
!       write(cFileOut,'(A,I0,A)') 'Vertlink_0/Vertlink_',PID,'.txt'
!       OPEN(FILE=ADJUSTL(tRIM(cFileOut)),unit=14142)
!       WRITE(14142,*) coarse%pVERTLINK(pID,:)
!       CLOSE(14142) 
!      END DO
!     END IF
    
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    call ztime(time1)
    If (myid.eq.1) write(*,*) "ParentComm stage 0 :: ",time1-time0 
    time0 = time1
    
    
    !
    ! --------------------------------------------------------------------------------------------------------
    !
    IF (myid.eq.master) THEN
     pNEL=NEL
     IF (ALLOCATED(dCOORDINATES)) deallocate(dCOORDINATES)
     ALLOCATE (dCOORDINATES(3,pNEL))
     DO I=1,pNEL
       dCOORDINATES(1,I)=0d0
       dCOORDINATES(2,I)=0d0
       dCOORDINATES(3,I)=0d0

       DO J=1,8
        dCOORDINATES(1,I) = dCOORDINATES(1,I) + DCORVG(1,KVERT(J,I))
        dCOORDINATES(2,I) = dCOORDINATES(2,I) + DCORVG(2,KVERT(J,I))
        dCOORDINATES(3,I) = dCOORDINATES(3,I) + DCORVG(3,KVERT(J,I))
       END DO

       dCOORDINATES(1,I)=0.125d0*dCOORDINATES(1,I)
       dCOORDINATES(2,I)=0.125d0*dCOORDINATES(2,I)
       dCOORDINATES(3,I)=0.125d0*dCOORDINATES(3,I)
     END DO
    ELSE
     ALLOCATE (pCOORDINATES(3,NEL))
     DO I=1,NEL
       pCOORDINATES(1,I)=0d0
       pCOORDINATES(2,I)=0d0
       pCOORDINATES(3,I)=0d0

       DO J=1,8
        pCOORDINATES(1,I) = pCOORDINATES(1,I) + DCORVG(1,KVERT(J,I))
        pCOORDINATES(2,I) = pCOORDINATES(2,I) + DCORVG(2,KVERT(J,I))
        pCOORDINATES(3,I) = pCOORDINATES(3,I) + DCORVG(3,KVERT(J,I))
       END DO

       pCOORDINATES(1,I)=0.125d0*pCOORDINATES(1,I)
       pCOORDINATES(2,I)=0.125d0*pCOORDINATES(2,I)
       pCOORDINATES(3,I)=0.125d0*pCOORDINATES(3,I)
     END DO
    END IF
    
    call MPI_Bcast(pNEL, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    IF (myid.ne.master) THEN
     IF (ALLOCATED(dCOORDINATES)) deallocate(dCOORDINATES)
     ALLOCATE (dCOORDINATES(3,pNEL))
    END IF
    
    call MPI_Bcast(dCOORDINATES, 3*pNEL, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    IF (myid.ne.master) THEN
      CALL InitOctTree(dCOORDINATES,pNEL)
      
       DO I=1,NEL
       
        CALL FindInOctTree(dCOORDINATES,pNEL,pCOORDINATES(:,I),J,DIST)
        
        IF (J.lt.0) then
         WRITE(*,*) I,"PROBLEM of elem assignement ..."
        end if
        IF (DIST.LT.DEpsPrec) THEN 
         coarse%myELEMLINK(I)=J
        END IF
       
       END DO
    
      CALL FreeOctTree()
    END IF

    call MPI_allgather(NEL, 1, MPI_INTEGER, sendcounts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    sendcounts(0) = 0

    displs = 0
    do i = 2, numnodes+1
      displs(i) = displs(i-1) + sendcounts(i-1)
    end do
    
    if (myid.eq.master) then
      if (allocated(gathered_data)) deallocate(gathered_data)
      allocate(gathered_data(displs(numnodes+1)))
      n=0
    else 
     n = nel
    endif
    
    call MPI_Gatherv(coarse%myELEMLINK, n, MPI_INTEGER, &
                   gathered_data, sendcounts, displs, &
                   MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
                   
    IF (myid.eq.0) THEN
      coarse%pELEMLINK = 0
      DO pID=1,subnodes
       j = 0
       DO i=displs(pID)+1, displs(pID+1)
        j = j + 1
        coarse%pELEMLINK(pID,J)=gathered_data(i)
       END DO
      end do
    END IF

!     IF (myid.eq.0) THEN
!      DO pID=1,subnodes
!       write(cFileOut,'(A,I0,A)') 'Vertlink_0/Vertlink_',PID,'.txt'
!       OPEN(FILE=ADJUSTL(tRIM(cFileOut)),unit=14142)
!       WRITE(14142,*) coarse%pELEMLINK(pID,:)
!       CLOSE(14142) 
!      END DO
!     END IF
    
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    call ztime(time1)
    If (myid.eq.1) write(*,*) "ParentComm stage 1 :: ",time1-time0
    time0 = time1
    
    IF (ALLOCATED(pCOORDINATES)) DEALLOCATE(pCOORDINATES)
    IF (ALLOCATED(dCOORDINATES)) DEALLOCATE(dCOORDINATES)

    !--------------------------------------------------------------
    !CREATING MAPPING STRUCTURE FOR MASTER --> ASSISTANT //FACES
    !--------------------------------------------------------------

    if (myid.eq.master) THEN
     pNAT = NAT
     allocate(dCOORDINATES(3,NAT))
     DO I=1,pNAT
      dCOORDINATES(1,I)=DCORAG(1,I)
      dCOORDINATES(2,I)=DCORAG(2,I)
      dCOORDINATES(3,I)=DCORAG(3,I)
     END DO
    end if
    
   call MPI_Bcast(pNAT, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    IF (myid.ne.master) THEN
     ALLOCATE (dCOORDINATES(3,pNAT))
    END IF
    
    call MPI_Bcast(dCOORDINATES, 3*pNAT, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     
    IF (myid.ne.master) THEN
    
      allocate(myFACELINK(2,pNAT))
      allocate(myFACEPINK(2,pNAT))
      myFACELINK = 0
      
      CALL InitOctTree(dCOORDINATES,pNAT)
      
      DO I=1,NAT
      
       CALL FindInOctTree(dCOORDINATES,pNAT,DCORAG(:,I),J,DIST)
       
       IF (J.lt.0) then
        WRITE(*,*) I,"PROBLEM of vert assignement ..."
       end if
       IF (DIST.LT.DEpsPrec) THEN 
        myFACELINK(1,J)=I
        myFACELINK(2,J)=myid
       END IF
      
       CALL FindInPeriodicOctTree(dCOORDINATES,nat,DCORAG(:,I),J,DIST,dPeriodicity)
       IF (DIST.LT.DEpsPrec) THEN 
        myFACELINK(1,J)=I
        myFACELINK(2,J)=myid
       END IF
      
      END DO
    
      CALL FreeOctTree()
      
      call MPI_Allreduce(myFACELINK, myFACEPINK, 2*pNAT, MPI_INTEGER, MPI_SUM, MPI_COMM_SUBS, ierr)      
      
      ALLOCATE(nFACELISTS(subnodes))
      nFACELISTS = 0
      
      DO I=1,pNAT
       IF (myFACELINK(1,I).ne.0) THEN
        IF (myFACELINK(1,I).ne.myFACEPINK(1,I).ne.0) THEN
         pID = myid
         PJD = myFACEPINK(2,I) - myid
         iFace = myFACELINK(1,I)
         jFace = myFACEPINK(1,I) - myFACELINK(1,I)
         nFACELISTS(PJD) = nFACELISTS(PJD) + 1
        END IF
       END IF
      END DO
      
      ALLOCATE(mg_mpi(1:9))
      ALLOCATE(mg_mpi(1)%UE(1:NAT))
      mg_mpi(1)%NeighNum = 0
      
      DO i=1,subnodes
       IF (nFACELISTS(i).ne.0) THEN
        mg_mpi(1)%NeighNum = mg_mpi(1)%NeighNum + 1
       END IF
      END dO
      ALLOCATE(mg_mpi(1)%parST(mg_mpi(1)%NeighNum))
      
      pID = 0
      DO i=1,subnodes
       IF (nFACELISTS(i).ne.0) THEN
        pID = pID + 1
        mg_mpi(1)%parST(pID)%Neigh=i
        mg_mpi(1)%parST(pID)%Num=nFACELISTS(i)
        ALLOCATE(mg_mpi(1)%parST(pID)%FaceLink(2,1:mg_mpi(1)%parST(pID)%Num))        
        ALLOCATE(mg_mpi(1)%parST(pID)%SideLink(1:mg_mpi(1)%parST(pID)%Num))
        ALLOCATE(mg_mpi(1)%parST(pID)%ElemLink(2,1:mg_mpi(1)%parST(pID)%Num))
        ALLOCATE(mg_mpi(1)%parST(pID)%SDVect(1:mg_mpi(1)%parST(pID)%Num))
        ALLOCATE(mg_mpi(1)%parST(pID)%RDVect(1:mg_mpi(1)%parST(pID)%Num))
        ALLOCATE(mg_mpi(1)%parST(pID)%SVVect(1:mg_mpi(1)%parST(pID)%Num))
        ALLOCATE(mg_mpi(1)%parST(pID)%RVVect(1:mg_mpi(1)%parST(pID)%Num))
        ALLOCATE(mg_mpi(1)%parST(pID)%PE(1:mg_mpi(1)%parST(pID)%Num))
       END IF
      END dO
       
      nFACELISTS = 0
      DO pID=1,mg_mpi(1)%NeighNum
       pJD = mg_mpi(1)%parST(pID)%Neigh
       DO I=1,pNAT
        IF (myFACELINK(1,I).ne.0) THEN
         IF (myFACELINK(1,I).ne.myFACEPINK(1,I).ne.0) THEN
          IF (pJD.eq.myFACEPINK(2,I) - myid) THEN
           iFace = myFACELINK(1,I)
           jFace = myFACEPINK(1,I) - myFACELINK(1,I)
           nFACELISTS(pJD) = nFACELISTS(pJD) + 1
           mg_mpi(1)%parST(pID)%FaceLink(1,nFACELISTS(pJD)) = iFace
           mg_mpi(1)%parST(pID)%FaceLink(2,nFACELISTS(pJD)) = iFace
          END IF
         END IF
        END IF
       
       END DO
      ENDDO
      
      
      DO pID=1,mg_mpi(1)%NeighNum
       DO K=1,mg_mpi(1)%parST(pID)%Num
        DO I=1,NEL
         DO J=1,6
          JI = KAREA(J,I)
          IF (JI.EQ.mg_mpi(1)%parST(pID)%FaceLink(1,K)) THEN
            mg_mpi(1)%parST(pID)%ElemLink(1,K)=I
            mg_mpi(1)%parST(pID)%SideLink(K)=J
          END IF
         END DO
        END DO
       END DO
       DO K=1,mg_mpi(1)%parST(pID)%Num
        DO I=1,NEL
         DO J=1,6
          JI = KAREA(J,I)
          IF (JI.EQ.mg_mpi(1)%parST(pID)%FaceLink(2,K)) THEN
            mg_mpi(1)%parST(pID)%ElemLink(2,K)=I
          END IF
         END DO
        END DO
       END DO
      END DO

!       write(cFileOut,'(A,I0,A)') 'Vertlink_0/Vertlink_',myid,'.txt'
!       OPEN(FILE=ADJUSTL(tRIM(cFileOut)),unit=14142)
!       DO pID=1,mg_mpi(1)%NeighNum
!        pJD = mg_mpi(1)%parST(pID)%Neigh
!        write(14142,*) myid,' to ', pJD, ' :: ',mg_mpi(1)%parST(pID)%FaceLink(1,:), ' :: ',mg_mpi(1)%parST(pID)%FaceLink(2,:)
!       END DO
!       CLOSE(14142) 
      
!       WRITE(*,*) myid, ' - ', nFACELISTS
      
    END IF

    
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    call ztime(time1)
    If (myid.eq.1) write(*,*) "ParentComm stage 2 :: ",time1-time0
    time0 = time1
    IF (ALLOCATED(pCOORDINATES)) DEALLOCATE(pCOORDINATES)
    IF (ALLOCATED(dCOORDINATES)) DEALLOCATE(dCOORDINATES)

END 
