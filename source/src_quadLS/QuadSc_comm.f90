! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E012RecCOMM_Init()
use, INTRINSIC :: ISO_C_BINDING
use mpi
use shared_memory_module, only : get_shared_memory_INT,get_shared_memory_DBL,myPC
USE PP3D_MPI, ONLY : ierr,myid,master,numnodes,subnodes,MPI_COMM_SUBS,MPI_COMM_SUBGROUP
USE PP3D_MPI, ONLY : SENDD_myMPI,RECVD_myMPI,SENDK_myMPI,RECVK_myMPI,MGE013,SENDI_myMPI,RECVI_myMPI
USE def_feat, ONLY: ILEV
USE var_QuadScalar,ONLY:knvt,knet,knat,knel,myStat,nlmax
use var_QuadScalar, ONLY : myRecComm

implicit none
INTEGER STATUS(MPI_STATUS_SIZE)

INTEGER I,J,pID,pJD,nSIZE,nEIGH,iLOC,iSHIFT,iAUX,jAUX,nXX
INTEGER II,JJ
integer, allocatable :: NumberOfMyRecords(:),NumberOfAllRecords(:,:)!,StartOfAllRecords(:,:)
! integer(kind=MPI_ADDRESS_KIND) :: win_size
INTEGER :: dblesize=8,intsize=4

character(len=256) :: cFMT
REAL*4 tt0,tt1

IF (myid.eq.MASTER) return

if (.not.allocated(myPC)) allocate(myPC(1:nlmax))
 
! CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
CALL ztime(tt0)

allocate (NumberOfMyRecords(myRecComm%NumHosts),NumberOfAllRecords(myRecComm%NumHosts,myRecComm%NumNodes))
allocate (myPC(ILEV)%s(myRecComm%NumHosts),myPC(ILEV)%r(myRecComm%NumHosts))
allocate (myPC(ILEV)%winS(myRecComm%NumHosts),myPC(ILEV)%ptrS(myRecComm%NumHosts))
allocate (myPC(ILEV)%winR(myRecComm%NumHosts),myPC(ILEV)%ptrR(myRecComm%NumHosts))
allocate (myPC(ILEV)%StartOfAllRecords(myRecComm%NumHosts,myRecComm%NumNodes+1))
myPC(ILEV)%StartOfAllRecords = 0

allocate(myPC(ILEV)%CODECs_win(myRecComm%NumHosts),myPC(ILEV)%CODECsptr(myRecComm%NumHosts),myPC(ILEV)%CODECs(myRecComm%NumHosts))
allocate(myPC(ILEV)%CODECr_win(myRecComm%NumHosts),myPC(ILEV)%CODECrptr(myRecComm%NumHosts),myPC(ILEV)%CODECr(myRecComm%NumHosts))

! Lets collect the number of vaules which will have to be exchanged to the different node-groups.
! NumberOfMyRecords will collect the number of values from this particular process to be sent to the processes on the other hosts

NumberOfMyRecords = 0
NumberOfAllRecords = 0
DO pID=1,subnodes
 pJD = myRecComm%groupIDs(pID)
 IF (MGE013(ILEV)%SP(pID)%Num.GT.0) THEN
  NumberOfMyRecords(pJD) = NumberOfMyRecords(pJD) + 1
 END IF
END DO

call MPI_allgather(NumberOfMyRecords, myRecComm%NumHosts, MPI_INTEGER, NumberOfAllRecords, myRecComm%NumHosts, MPI_INTEGER, MPI_COMM_SUBGROUP, ierr)

CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

DO pJD=1,myRecComm%NumHosts
  myPC(ILEV)%StartOfAllRecords(pJD,1) = 0
  do pID=1,myRecComm%NumNodes
   myPC(ILEV)%StartOfAllRecords(pJD,pID+1) = myPC(ILEV)%StartOfAllRecords(pJD,pID) + NumberOfAllRecords(pJD,pID)
  end do
END DO

! if (myRecComm%myid.eq.0) then 
!  DO pID = 1,myRecComm%NumHosts
!   write(*,'(I5,A,100I5)') myRecComm%MyNodeGroup, 'AllRecords: ',myPC(ILEV)%StartOfAllRecords(pID,:)
!  end do
! end if

! StartOfAllRecords(:,myRecComm%NumNodes+1) is the size of the Records
! we allocate the CODEC buffer 1 element larger so that the total size can also be included aside of the displacements
DO pJD=1,myRecComm%NumHosts
  nSIZE = myPC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1)
  if (nSIZE.gt.0) then
   CALL get_shared_memory_INT(MPI_COMM_SUBGROUP,myPC(ILEV)%CODECs(pJD)%x,(nSIZE+1)*3, myPC(ILEV)%CODECs_win(pJD),myPC(ILEV)%CODECsptr(pJD),ierr)
   myPC(ILEV)%CODECs(pJD)%x = 0
  end if
!   if (nSIZE.gt.0.and.pJD.ne.myRecComm%MyNodeGroup) then
!    CALL get_shared_memory_INT(MPI_COMM_SUBGROUP,myPC(ILEV)%CODECr(pJD)%x,(nSIZE+1)*3, myPC(ILEV)%CODECr_win(pJD),myPC(ILEV)%CODECrptr(pJD),ierr)
!    myPC(ILEV)%CODECr(pJD)%x = 0
!   end if
!   if (myRecComm%myid.eq.0) then 
!    write(*,*) myRecComm%MyNodeGroup,' to ', pJD,' size: ',nSIZE,size(myPC(ILEV)%CODECs(pJD)%x)
!   end if
END DO

CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

! We need to to create a displacement info which will help us later to define the start-position this particular process to fill up 
! the node-collected exchange buffers.
! Fill up the CODECs with the [sender,receiver,datasize] information
DO pID=1,subnodes
 pJD = myRecComm%groupIDs(pID)
 IF (MGE013(ILEV)%SP(pID)%Num.GT.0) THEN
  iAux = myPC(ILEV)%StartOfAllRecords(pJD,myRecComm%myid+1)
  myPC(ILEV)%CODECs(pJD)%x(3*iAux+1) = myid
  myPC(ILEV)%CODECs(pJD)%x(3*iAux+2) = pID
  myPC(ILEV)%CODECs(pJD)%x(3*iAux+3) = MGE013(ILEV)%SP(pID)%nElems(1)
  myPC(ILEV)%StartOfAllRecords(pJD,myRecComm%myid+1) = myPC(ILEV)%StartOfAllRecords(pJD,myRecComm%myid+1) + 1
 END IF
END DO

CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

! ! ! ! ! ! ! ! if (myRecComm%myid.eq.0.and.myRecComm%MyNodeGroup.eq.1) then 
! ! ! ! ! ! ! if (myRecComm%myid.eq.0) then
! ! ! ! ! ! !  DO pJD=1,myRecComm%NumHosts
! ! ! ! ! ! !    write(*,*) "S ",myid,pjd, " : ",myPC(ILEV)%CODECs(pJD)%x
! ! ! ! ! ! ! !   DO i=0,size(myPC(ILEV)%CODECs(pJD)%x)/3-1
! ! ! ! ! ! ! !    write(*,*) pJD,"S",myPC(ILEV)%CODECs(pJD)%x(3*i+1:3*i+3)
! ! ! ! ! ! ! !   END DO
! ! ! ! ! ! !  END DO
! ! ! ! ! ! ! end if
! ! ! ! ! ! ! !pause

! convert the datasizes to displacements
if (myRecComm%myid.eq.0) then
 DO pJD=1,myRecComm%NumHosts
  if (myPC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0) then 
   nSIZE = size(myPC(ILEV)%CODECs(pJD)%x)/3-1
!   IF (nSIZE.gt.0) THEN
   iaux = myPC(ILEV)%CODECs(pJD)%x(3)
   myPC(ILEV)%CODECs(pJD)%x(3) = 0
   DO i=1,nSIZE
    jaux = myPC(ILEV)%CODECs(pJD)%x(3*i+3)
    myPC(ILEV)%CODECs(pJD)%x(3*i+3) = iaux
    iaux = iaux + jaux
   END DO
  END IF
 END DO
END IF

CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

! if (myRecComm%myid.eq.0.and.myRecComm%MyNodeGroup.eq.1) then 
! ! if (myid.eq.24) then 
!  DO pJD=1,myRecComm%NumHosts
! !   if (NumberOfAllRecords(pJD,myRecComm%myid).gt.0) then 
!    DO i=0,size(myPC(ILEV)%CODECs(pJD)%x)/3-1
!     write(*,*) pJD,"S",myPC(ILEV)%CODECs(pJD)%x(3*i+1:3*i+3)
!    END DO
! !   end if
!  END DO
! end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
 

CALL get_shared_memory_INT(MPI_COMM_SUBGROUP,myPC(ILEV)%sendSIZE,myRecComm%NumHosts, myPC(ILEV)%SIZEs_win,myPC(ILEV)%SIZEsptr,ierr)
CALL get_shared_memory_INT(MPI_COMM_SUBGROUP,myPC(ILEV)%recvSIZE,myRecComm%NumHosts, myPC(ILEV)%SIZEr_win,myPC(ILEV)%SIZErptr,ierr)

myPC(ILEV)%sendSIZE(1:myRecComm%NumHosts) = 0
myPC(ILEV)%recvSIZE(1:myRecComm%NumHosts) = 0

CALL MPI_BARRIER(MPI_COMM_SUBGROUP,IERR)

IF (myRecComm%myID.eq.0) THEN
 DO pID=1,myRecComm%NumHosts
  IF (myRecComm%myNodeGroup.NE.pID) THEN
   IF (myPC(ILEV)%StartOfAllRecords(pID,myRecComm%NumNodes+1).gt.0) then 
    iLoc = size(myPC(ILEV)%CODECs(pID)%x)/3
    myPC(ILEV)%sendSIZE(pID) = iLoc
    CALL SENDI_myMPI(iLoc,myRecComm%hostleaders(pID))
   END IF
  ELSE
   iLoc = size(myPC(ILEV)%CODECs(pID)%x)/3
   myPC(ILEV)%sendSIZE(pID) = iLoc
   DO pJD=1,myRecComm%NumHosts
    if (myPC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0.and.pJD.ne.myRecComm%myNodeGroup) then 
     CALL RECVI_myMPI(iLoc,myRecComm%hostleaders(pJD))
     myPC(ILEV)%recvSIZE(pJD) = iLoc
    END IF
   END DO
  END IF
 END DO
END IF

CALL MPI_BARRIER(MPI_COMM_SUBGROUP,IERR)

DO pJD=1,myRecComm%NumHosts
  nSIZE = myPC(ILEV)%recvSIZE(pJD)
  if (nSIZE.gt.0.and.pJD.ne.myRecComm%MyNodeGroup) then
   CALL get_shared_memory_INT(MPI_COMM_SUBGROUP,myPC(ILEV)%CODECr(pJD)%x,nSIZE*3, myPC(ILEV)%CODECr_win(pJD),myPC(ILEV)%CODECrptr(pJD),ierr)
   myPC(ILEV)%CODECr(pJD)%x = 0
  end if
END DO

IF (myRecComm%myID.eq.0) THEN
 ! Communication of the group-leaders
 DO pID=1,myRecComm%NumHosts
  IF (myRecComm%myNodeGroup.NE.pID) THEN
   IF (myPC(ILEV)%StartOfAllRecords(pID,myRecComm%NumNodes+1).gt.0) then 
    iaux = myPC(ILEV)%sendSIZE(pID)
    CALL SENDK_myMPI(myPC(ILEV)%CODECs(pID)%x,3*iaux,myRecComm%hostleaders(pID))
!     write(*,*) myid,PID,iaux,size(myPC(ILEV)%CODECs(pID)%x)!myPC(ILEV)%CODECs(pID)%x
   END IF
  ELSE
   DO pJD=1,myRecComm%NumHosts
    if (myPC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0.and.pJD.ne.myRecComm%myNodeGroup) then 
     iaux = myPC(ILEV)%recvSIZE(pJD)
     CALL RECVK_myMPI(myPC(ILEV)%CODECr(pJD)%x,3*iaux,myRecComm%hostleaders(pJD))
!      write(*,*) myid,PJD,myPC(ILEV)%CODECr(pJD)%x
    END IF
   END DO
  END IF
 END DO
END IF

CALL MPI_BARRIER(MPI_COMM_SUBGROUP,IERR)
! if (myRecComm%myid.eq.0) then
!  DO pJD=1,myRecComm%NumHosts
! !   write(*,*) "R ",myid,pjd, " : ",myPC(ILEV)%CODECr(pJD)%x
!    DO i=0,size(myPC(ILEV)%CODECr(pJD)%x)/3-1
!     write(*,*) myid,pJD,"R",myPC(ILEV)%CODECr(pJD)%x(3*i+1:3*i+3)
!    END DO
!  END DO
! end if
! pause

! Creating the node-exchange buffers
DO pJD=1,myRecComm%NumHosts
 if (myPC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0) then 
  iaux = myPC(ILEV)%sendSIZE(pJD)
  IF (iaux.gt.0) then
   nSIZE = myPC(ILEV)%CODECs(pJD)%x(3*iaux)
   CALL get_shared_memory_DBL(MPI_COMM_SUBGROUP,myPC(ILEV)%s(pJD)%x,4*nSIZE, myPC(ILEV)%winS(pJD),myPC(ILEV)%ptrS(pJD),ierr)
  END IF
  
  iaux = myPC(ILEV)%recvSIZE(pJD)
  IF (iaux.gt.0.and.pJD.ne.myRecComm%myNodeGroup) then
   nSIZE = myPC(ILEV)%CODECr(pJD)%x(3*iaux)
   CALL get_shared_memory_DBL(MPI_COMM_SUBGROUP,myPC(ILEV)%r(pJD)%x,4*nSIZE, myPC(ILEV)%winR(pJD),myPC(ILEV)%ptrR(pJD),ierr)
  END IF
 END IF
END DO

CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

if (allocated(NumberOfMyRecords)) deallocate (NumberOfMyRecords)
if (allocated(NumberOfAllRecords)) deallocate (NumberOfAllRecords)
myPC(ILEV)%bPrepared = .true.

END SUBROUTINE E012RecCOMM_Init
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013RecCOMM_Init() !ok
use, INTRINSIC :: ISO_C_BINDING
use mpi
use shared_memory_module, only : get_shared_memory_INT,get_shared_memory_DBL,myRC
USE PP3D_MPI, ONLY : ierr,myid,master,numnodes,subnodes,MPI_COMM_SUBS,MPI_COMM_SUBGROUP
USE PP3D_MPI, ONLY : SENDD_myMPI,RECVD_myMPI,SENDK_myMPI,RECVK_myMPI,MGE013
USE def_feat, ONLY: ILEV
USE var_QuadScalar,ONLY:knvt,knet,knat,knel,myStat,nlmax
use var_QuadScalar, ONLY : myRecComm

implicit none
INTEGER STATUS(MPI_STATUS_SIZE)

INTEGER I,J,pID,pJD,nSIZE,nEIGH,iLOC,iSHIFT,iAUX,jAUX,nXX
integer, allocatable :: NumberOfMyRecords(:),NumberOfAllRecords(:,:)
INTEGER :: dblesize=8,intsize=4

character(len=256) :: cFMT
REAL*4 tt0,tt1

IF (myid.eq.MASTER) return

if (.not.allocated(myRC)) allocate(myRC(1:nlmax))
 
CALL ztime(tt0)

allocate (NumberOfMyRecords(myRecComm%NumHosts),NumberOfAllRecords(myRecComm%NumHosts,myRecComm%NumNodes))
allocate (myRC(ILEV)%s(myRecComm%NumHosts),myRC(ILEV)%r(myRecComm%NumHosts))
allocate (myRC(ILEV)%winS(myRecComm%NumHosts),myRC(ILEV)%ptrS(myRecComm%NumHosts))
allocate (myRC(ILEV)%winR(myRecComm%NumHosts),myRC(ILEV)%ptrR(myRecComm%NumHosts))
allocate (myRC(ILEV)%StartOfAllRecords(myRecComm%NumHosts,myRecComm%NumNodes+1))
myRC(ILEV)%StartOfAllRecords = 0

allocate(myRC(ILEV)%CODECs_win(myRecComm%NumHosts),myRC(ILEV)%CODECsptr(myRecComm%NumHosts),myRC(ILEV)%CODECs(myRecComm%NumHosts))
allocate(myRC(ILEV)%CODECr_win(myRecComm%NumHosts),myRC(ILEV)%CODECrptr(myRecComm%NumHosts),myRC(ILEV)%CODECr(myRecComm%NumHosts))

! Lets collect the number of vaules which will have to be exchanged to the different node-groups.
! NumberOfMyRecords will collect the number of values from this particular process to be sent to the processes on the other hosts

NumberOfMyRecords = 0
NumberOfAllRecords = 0
DO pID=1,subnodes
 pJD = myRecComm%groupIDs(pID)
 IF (MGE013(ILEV)%ST(pID)%Num.GT.0) THEN
  NumberOfMyRecords(pJD) = NumberOfMyRecords(pJD) + 1
 END IF
END DO

call MPI_allgather(NumberOfMyRecords, myRecComm%NumHosts, MPI_INTEGER, NumberOfAllRecords, myRecComm%NumHosts, MPI_INTEGER, MPI_COMM_SUBGROUP, ierr)

CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

DO pJD=1,myRecComm%NumHosts
  myRC(ILEV)%StartOfAllRecords(pJD,1) = 0
  do pID=1,myRecComm%NumNodes
   myRC(ILEV)%StartOfAllRecords(pJD,pID+1) = myRC(ILEV)%StartOfAllRecords(pJD,pID) + NumberOfAllRecords(pJD,pID)
  end do
END DO

! if (myRecComm%myid.eq.0) then 
!  DO pID = 1,myRecComm%NumHosts
!   write(*,'(I5,A,100I5)') myRecComm%MyNodeGroup, 'AllRecords: ',myRC(ILEV)%StartOfAllRecords(pID,:)
!  end do
! end if

! StartOfAllRecords(:,myRecComm%NumNodes+1) is the size of the Records
! we allocate the CODEC buffer 1 element larger so that the total size can also be included aside of the displacements
DO pJD=1,myRecComm%NumHosts
  nSIZE = myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1)
  if (nSIZE.gt.0) then
   CALL get_shared_memory_INT(MPI_COMM_SUBGROUP,myRC(ILEV)%CODECs(pJD)%x,(nSIZE+1)*3, myRC(ILEV)%CODECs_win(pJD),myRC(ILEV)%CODECsptr(pJD),ierr)
   myRC(ILEV)%CODECs(pJD)%x = 0
  end if
  if (nSIZE.gt.0.and.pJD.ne.myRecComm%MyNodeGroup) then
   CALL get_shared_memory_INT(MPI_COMM_SUBGROUP,myRC(ILEV)%CODECr(pJD)%x,(nSIZE+1)*3, myRC(ILEV)%CODECr_win(pJD),myRC(ILEV)%CODECrptr(pJD),ierr)
   myRC(ILEV)%CODECr(pJD)%x = 0
  end if
!   if (myRecComm%myid.eq.0) then 
!    write(*,*) myRecComm%MyNodeGroup,' to ', pJD,' size: ',nSIZE,size(myRC(ILEV)%CODECs(pJD)%x)
!   end if
END DO

CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

! We need to to create a displacement info which will help us later to define the start-position this particular process to fill up 
! the node-collected exchange buffers.
! Fill up the CODECs with the [sender,receiver,datasize] information
DO pID=1,subnodes
 pJD = myRecComm%groupIDs(pID)
 IF (MGE013(ILEV)%ST(pID)%Num.GT.0) THEN
  iAux = myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%myid+1)
  myRC(ILEV)%CODECs(pJD)%x(3*iAux+1) = myid
  myRC(ILEV)%CODECs(pJD)%x(3*iAux+2) = pID
  myRC(ILEV)%CODECs(pJD)%x(3*iAux+3) = MGE013(ILEV)%ST(pID)%Num
  myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%myid+1) = myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%myid+1) + 1
 END IF
END DO

CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

! if (myRecComm%myid.eq.0.and.myRecComm%MyNodeGroup.eq.1) then 
!  DO pJD=1,myRecComm%NumHosts
!   DO i=0,size(myRC(ILEV)%CODECs(pJD)%x)/3-1
!    write(*,*) pJD,"S",myRC(ILEV)%CODECs(pJD)%x(3*i+1:3*i+3)
!   END DO
!  END DO
! end if

! convert the datasizes to displacements
if (myRecComm%myid.eq.0) then
 DO pJD=1,myRecComm%NumHosts
  if (myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0) then 
   nSIZE = size(myRC(ILEV)%CODECs(pJD)%x)/3-1
!   IF (nSIZE.gt.0) THEN
   iaux = myRC(ILEV)%CODECs(pJD)%x(3)
   myRC(ILEV)%CODECs(pJD)%x(3) = 0
   DO i=1,nSIZE
    jaux = myRC(ILEV)%CODECs(pJD)%x(3*i+3)
    myRC(ILEV)%CODECs(pJD)%x(3*i+3) = iaux
    iaux = iaux + jaux
   END DO
  END IF
 END DO
END IF

CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

! if (myRecComm%myid.eq.0.and.myRecComm%MyNodeGroup.eq.5) then 
! if (myid.eq.24) then 
!  DO pJD=1,myRecComm%NumHosts
! !   if (NumberOfAllRecords(pJD,myRecComm%myid).gt.0) then 
!    DO i=0,size(myRC(ILEV)%CODECs(pJD)%x)/3-1
!     write(*,*) pJD,"S",myRC(ILEV)%CODECs(pJD)%x(3*i+1:3*i+3)
!    END DO
! !   end if
!  END DO
! end if

! pause

! Creating the node-exchange buffers
DO pJD=1,myRecComm%NumHosts
 if (myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0) then 
  iaux = size(myRC(ILEV)%CODECs(pJD)%x)/3
  IF (iaux.gt.0) then
   nSIZE = myRC(ILEV)%CODECs(pJD)%x(3*iaux)
   CALL get_shared_memory_DBL(MPI_COMM_SUBGROUP,myRC(ILEV)%s(pJD)%x,3*nSIZE, myRC(ILEV)%winS(pJD),myRC(ILEV)%ptrS(pJD),ierr)
  END IF
  
  IF (iaux.gt.0.and.pJD.ne.myRecComm%myNodeGroup) then
   nSIZE = myRC(ILEV)%CODECs(pJD)%x(3*iaux)
   CALL get_shared_memory_DBL(MPI_COMM_SUBGROUP,myRC(ILEV)%r(pJD)%x,3*nSIZE, myRC(ILEV)%winR(pJD),myRC(ILEV)%ptrR(pJD),ierr)
  END IF
 END IF
END DO

CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

! Communication of the group-leaders
IF (myRecComm%myID.eq.0) THEN
 DO pID=1,myRecComm%NumHosts
  IF (myRecComm%myNodeGroup.NE.pID) THEN
   IF (myRC(ILEV)%StartOfAllRecords(pID,myRecComm%NumNodes+1).gt.0) then 
    iaux = size(myRC(ILEV)%CODECs(pID)%x)/3
    CALL SENDK_myMPI(myRC(ILEV)%CODECs(pID)%x,3*iaux,myRecComm%hostleaders(pID))
   END IF
  ELSE
   DO pJD=1,myRecComm%NumHosts
    if (myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0.and.pJD.ne.myRecComm%myNodeGroup) then 
     iaux = size(myRC(ILEV)%CODECs(pJD)%x)/3
     CALL RECVK_myMPI(myRC(ILEV)%CODECr(pJD)%x,3*iaux,myRecComm%hostleaders(pJD))
    END IF
   END DO
  END IF
 END DO
END IF

if (allocated(NumberOfMyRecords)) deallocate (NumberOfMyRecords)
if (allocated(NumberOfAllRecords)) deallocate (NumberOfAllRecords)

myRC(ILEV)%bPrepared = .TRUE.

END SUBROUTINE E013RecCOMM_Init
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013UVWSumRec(FX) !ok
use, INTRINSIC :: ISO_C_BINDING
use mpi
use shared_memory_module, only : get_shared_memory_INT,get_shared_memory_DBL,myRC
USE PP3D_MPI, ONLY : ierr,myid,master,numnodes,subnodes,MPI_COMM_SUBS,MPI_COMM_SUBGROUP
USE PP3D_MPI, ONLY : SENDD_myMPI,RECVD_myMPI,SENDK_myMPI,RECVK_myMPI,MGE013
USE def_feat, ONLY: ILEV
USE var_QuadScalar,ONLY:knvt,knet,knat,knel,myStat,nlmax
use var_QuadScalar, ONLY : myRecComm

implicit none
REAL*8  FX(*)

INTEGER send_req(numnodes),recv_req(numnodes)
INTEGER STATUS(MPI_STATUS_SIZE)

INTEGER I,J,pID,pJD,nSIZE,nEIGH,iLOC,iSHIFT,iAUX,jAUX,nXX
INTEGER MEQ,MEQ1,MEQ2,MEQ3
INTEGER LEQ,LEQ1,LEQ2,LEQ3
INTEGER :: dblesize=8,intsize=4

character(len=256) :: cFMT
REAL*4 tt0,tt1

IF (myid.eq.MASTER) return

send_req = MPI_REQUEST_NULL
recv_req = MPI_REQUEST_NULL

if (.not.allocated(myRC)) allocate(myRC(1:nlmax))

LEQ = KNVT(ILEV) + KNAT(ILEV) + KNET(ILEV) + KNEL(ILEV)
LEQ1 =0
LEQ2 =LEQ
LEQ3 =2*LEQ
 
! CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
CALL ztime(tt0)

IF (.not.myRC(ILEV)%bPrepared) CALL E013RecCOMM_Init()

! extract the data to be sent
! All processes from the given host upload the data packages to be sent
DO pJD=1,myRecComm%NumHosts

 if (myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0) then 
  iaux = size(myRC(ILEV)%CODECs(pJD)%x)/3-1
!  IF (iaux.gt.0) THEN
  
   DO j=1,iaux
    jaux = 3*(j-1)
    IF (myRC(ILEV)%CODECs(pJD)%x(jaux+1).eq.myid) THEN
     pID  = myRC(ILEV)%CODECs(pJD)%x(jaux+2) ! receiver
     iLoc = myRC(ILEV)%CODECs(pJD)%x(jaux+3) ! Start location for reading the data
     nSIZE = myRC(ILEV)%CODECs(pJD)%x(jaux+6) - myRC(ILEV)%CODECs(pJD)%x(jaux+3) !MGE013(ILEV)%ST(pID)%Num
 !     write(*,*) myid,pID,nSIZE,MGE013(ILEV)%ST(pID)%Num
     DO I=1,nSIZE
       iShift = MGE013(ILEV)%ST(pID)%VertLink(1,I)
       myRC(ILEV)%s(pJD)%x(3*(iLoc + I - 1) + 1) = FX(LEQ1 + iShift)
       myRC(ILEV)%s(pJD)%x(3*(iLoc + I - 1) + 2) = FX(LEQ2 + iShift)
       myRC(ILEV)%s(pJD)%x(3*(iLoc + I - 1) + 3) = FX(LEQ3 + iShift)
     END DO
    END IF
   END DO
 END IF
 
END DO

CALL MPI_BARRIER(MPI_COMM_SUBGROUP,IERR)

! Communication of the group-leaders
IF (myRecComm%myID.eq.0) THEN

 DO pID=1,myRecComm%NumHosts
 
  IF (myRecComm%myNodeGroup.NE.pID) THEN
  
   if (myRC(ILEV)%StartOfAllRecords(pID,myRecComm%NumNodes+1).gt.0) then 
     iaux = size(myRC(ILEV)%CODECs(pID)%x)/3
     nSIZE = myRC(ILEV)%CODECs(pID)%x(3*iaux)
     
!      CALL SENDD_myMPI(myRC(ILEV)%s(pID)%x,3*nSIZE,myRecComm%hostleaders(pID))
     CALL MPI_ISEND(myRC(ILEV)%s(pID)%x,3*nSIZE,MPI_DOUBLE_PRECISION,myRecComm%hostleaders(pID),1001,MPI_COMM_WORLD,send_req(pID),IERR)
    END IF
   
  ELSE
   DO pJD=1,myRecComm%NumHosts
   
    if (myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0.and.pJD.ne.myRecComm%myNodeGroup) then 
     iaux = size(myRC(ILEV)%CODECs(pJD)%x)/3
     nSIZE = myRC(ILEV)%CODECs(pJD)%x(3*iaux)
!      CALL RECVD_myMPI(myRC(ILEV)%r(pJD)%x,3*nSIZE,myRecComm%hostleaders(pJD))
     CALL MPI_IRECV(myRC(ILEV)%r(pJD)%x,3*nSIZE,MPI_DOUBLE_PRECISION,myRecComm%hostleaders(pJD),1001,MPI_COMM_WORLD,recv_req(pJD),IERR)
    END IF
   END DO
  END IF
  
 END DO
END IF

IF (myRecComm%myID.eq.0) THEN
 DO pID=1,myRecComm%NumHosts
  IF (myRC(ILEV)%StartOfAllRecords(pID,myRecComm%NumNodes+1).gt.0.and.myRecComm%myNodeGroup.NE.pID) then 
     CALL MPI_Wait(send_req(pID),STATUS, IERR )
     CALL MPI_Wait(recv_req(pID),STATUS, IERR )
  END IF
 END DO
END IF

CALL MPI_BARRIER(MPI_COMM_SUBGROUP,IERR)

! if (myRecComm%myid.eq.0) then 
!  DO pJD=1,myRecComm%NumHosts
!   DO i=0,size(myRC(ILEV)%CODECr(pJD)%x)/3-1
!    write(*,*) pJD,"X",myRC(ILEV)%CODECr(pJD)%x(3*i+1:3*i+3)
!   END DO
!  END DO
! end if
! 
! write(*,*) myid,'done'
! pause

! All processes from the given host extract their own data packages
DO pJD=1,myRecComm%NumHosts

 IF (pJD.ne.myRecComm%myNodeGroup) then
  if (myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0) then 
    iaux = size(myRC(ILEV)%CODECr(pJD)%x)/3-1
!   IF (iaux.gt.0) THEN
   
    DO j=1,iaux
     jaux = 3*(j-1)
     IF (myRC(ILEV)%CODECr(pJD)%x(jaux+2).eq.myid) THEN ! receiver is me
      pID  = myRC(ILEV)%CODECr(pJD)%x(jaux+1) ! sender
      iLoc = myRC(ILEV)%CODECr(pJD)%x(jaux+3) ! Start location for reading the data
      nSIZE = myRC(ILEV)%CODECr(pJD)%x(jaux+6) - myRC(ILEV)%CODECr(pJD)%x(jaux+3) !MGE013(ILEV)%ST(pID)%Num
  !     write(*,*) myid,pID,nSIZE,MGE013(ILEV)%ST(pID)%Num
      DO I=1,nSIZE
        iShift = MGE013(ILEV)%ST(pID)%VertLink(2,I)
        FX(LEQ1 + iShift) = FX(LEQ1 + iShift) + myRC(ILEV)%r(pJD)%x(3*(iLoc + I - 1) + 1)
        FX(LEQ2 + iShift) = FX(LEQ2 + iShift) + myRC(ILEV)%r(pJD)%x(3*(iLoc + I - 1) + 2)
        FX(LEQ3 + iShift) = FX(LEQ3 + iShift) + myRC(ILEV)%r(pJD)%x(3*(iLoc + I - 1) + 3)
      END DO
     END IF
    END DO
  END IF
  
 ELSE !pJD.eq.myRecComm%myNodeGroup)
 
  if (myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0) then 
   iaux = size(myRC(ILEV)%CODECs(pJD)%x)/3-1
!   IF (iaux.gt.0) THEN
   DO j=1,iaux
    jaux = 3*(j-1)
    IF (myRC(ILEV)%CODECs(pJD)%x(jaux+2).eq.myid) THEN ! receiver is me
     pID  = myRC(ILEV)%CODECs(pJD)%x(jaux+1) ! sender
     iLoc = myRC(ILEV)%CODECs(pJD)%x(jaux+3) ! Start location for reading the data
     nSIZE = myRC(ILEV)%CODECs(pJD)%x(jaux+6) - myRC(ILEV)%CODECs(pJD)%x(jaux+3) !MGE013(ILEV)%ST(pID)%Num
!     write(*,*) myid,pID,nSIZE,MGE013(ILEV)%ST(pID)%Num
     DO I=1,nSIZE
        iShift = MGE013(ILEV)%ST(pID)%VertLink(2,I)
        FX(LEQ1 + iShift) = FX(LEQ1 + iShift) + myRC(ILEV)%s(pJD)%x(3*(iLoc + I - 1) + 1)
        FX(LEQ2 + iShift) = FX(LEQ2 + iShift) + myRC(ILEV)%s(pJD)%x(3*(iLoc + I - 1) + 2)
        FX(LEQ3 + iShift) = FX(LEQ3 + iShift) + myRC(ILEV)%s(pJD)%x(3*(iLoc + I - 1) + 3)
     END DO
    END IF
   END DO
  END IF
   
 END IF
END DO

! CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

CALL ztime(tt1)
myStat%tCommV = myStat%tCommV + (tt1-tt0)
 
! write(*,*) myid,'done', (tt1-tt0)
! pause
 
END
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013SumRec(FX) !ok
use, INTRINSIC :: ISO_C_BINDING
use mpi
use shared_memory_module, only : get_shared_memory_INT,get_shared_memory_DBL,myRC
USE PP3D_MPI, ONLY : ierr,myid,master,numnodes,subnodes,MPI_COMM_SUBS,MPI_COMM_SUBGROUP
USE PP3D_MPI, ONLY : SENDD_myMPI,RECVD_myMPI,SENDK_myMPI,RECVK_myMPI,MGE013
USE def_feat, ONLY: ILEV
USE var_QuadScalar,ONLY:knvt,knet,knat,knel,myStat,nlmax
use var_QuadScalar, ONLY : myRecComm

implicit none
REAL*8  FX(*)

INTEGER send_req(numnodes),recv_req(numnodes)
INTEGER STATUS(MPI_STATUS_SIZE)
INTEGER I,J,pID,pJD,nSIZE,nEIGH,iLOC,iSHIFT,iAUX,jAUX,nXX
INTEGER MEQ,MEQ1,MEQ2,MEQ3
INTEGER LEQ,LEQ1,LEQ2,LEQ3
INTEGER :: dblesize=8,intsize=4

character(len=256) :: cFMT
REAL*4 tt0,tt1

IF (myid.eq.MASTER) return

send_req = MPI_REQUEST_NULL
recv_req = MPI_REQUEST_NULL

if (.not.allocated(myRC)) allocate(myRC(1:nlmax))

LEQ = KNVT(ILEV) + KNAT(ILEV) + KNET(ILEV) + KNEL(ILEV)
LEQ1 =0
LEQ2 =LEQ
LEQ3 =2*LEQ
 
! CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
CALL ztime(tt0)

IF (.not.myRC(ILEV)%bPrepared) CALL E013RecCOMM_Init()

! extract the data to be sent
! All processes from the given host upload the data packages to be sent
DO pJD=1,myRecComm%NumHosts

 if (myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0) then 
  iaux = size(myRC(ILEV)%CODECs(pJD)%x)/3-1
!  IF (iaux.gt.0) THEN
  
   DO j=1,iaux
    jaux = 3*(j-1)
    IF (myRC(ILEV)%CODECs(pJD)%x(jaux+1).eq.myid) THEN
     pID  = myRC(ILEV)%CODECs(pJD)%x(jaux+2) ! receiver
     iLoc = myRC(ILEV)%CODECs(pJD)%x(jaux+3) ! Start location for reading the data
     nSIZE = myRC(ILEV)%CODECs(pJD)%x(jaux+6) - myRC(ILEV)%CODECs(pJD)%x(jaux+3) !MGE013(ILEV)%ST(pID)%Num
 !     write(*,*) myid,pID,nSIZE,MGE013(ILEV)%ST(pID)%Num
     DO I=1,nSIZE
       iShift = MGE013(ILEV)%ST(pID)%VertLink(1,I)
       myRC(ILEV)%s(pJD)%x(iLoc + I ) = FX(iShift)
!        myRC(ILEV)%s(pJD)%x(3*(iLoc + I - 1) + 2) = FX(LEQ2 + iShift)
!        myRC(ILEV)%s(pJD)%x(3*(iLoc + I - 1) + 3) = FX(LEQ3 + iShift)
     END DO
    END IF
   END DO
 END IF
 
END DO

CALL MPI_BARRIER(MPI_COMM_SUBGROUP,IERR)

! Communication of the group-leaders
IF (myRecComm%myID.eq.0) THEN

 DO pID=1,myRecComm%NumHosts
 
  IF (myRecComm%myNodeGroup.NE.pID) THEN
  
   if (myRC(ILEV)%StartOfAllRecords(pID,myRecComm%NumNodes+1).gt.0) then 
     iaux = size(myRC(ILEV)%CODECs(pID)%x)/3
     nSIZE = myRC(ILEV)%CODECs(pID)%x(3*iaux)
     
!     CALL SENDD_myMPI(myRC(ILEV)%s(pID)%x,nSIZE,myRecComm%hostleaders(pID))
     CALL MPI_ISEND(myRC(ILEV)%s(pID)%x,nSIZE,MPI_DOUBLE_PRECISION,myRecComm%hostleaders(pID),1001,MPI_COMM_WORLD,send_req(pID),IERR)
    END IF
   
  ELSE
   DO pJD=1,myRecComm%NumHosts
   
    if (myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0.and.pJD.ne.myRecComm%myNodeGroup) then 
     iaux = size(myRC(ILEV)%CODECs(pJD)%x)/3
     nSIZE = myRC(ILEV)%CODECs(pJD)%x(3*iaux)
!     CALL RECVD_myMPI(myRC(ILEV)%r(pJD)%x,nSIZE,myRecComm%hostleaders(pJD))
     CALL MPI_IRECV(myRC(ILEV)%r(pJD)%x,nSIZE,MPI_DOUBLE_PRECISION,myRecComm%hostleaders(pJD),1001,MPI_COMM_WORLD,recv_req(pJD),IERR)
    END IF
   END DO
  END IF
  
 END DO
END IF

IF (myRecComm%myID.eq.0) THEN
 DO pID=1,myRecComm%NumHosts
  IF (myRC(ILEV)%StartOfAllRecords(pID,myRecComm%NumNodes+1).gt.0.and.myRecComm%myNodeGroup.NE.pID) then 
     CALL MPI_Wait(send_req(pID),STATUS, IERR )
     CALL MPI_Wait(recv_req(pID),STATUS, IERR )
  END IF
 END DO
END IF

CALL MPI_BARRIER(MPI_COMM_SUBGROUP,IERR)

! if (myRecComm%myid.eq.0) then 
!  DO pJD=1,myRecComm%NumHosts
!   DO i=0,size(myRC(ILEV)%CODECr(pJD)%x)/3-1
!    write(*,*) pJD,"X",myRC(ILEV)%CODECr(pJD)%x(3*i+1:3*i+3)
!   END DO
!  END DO
! end if
! 
! write(*,*) myid,'done'
! pause

! All processes from the given host extract their own data packages
DO pJD=1,myRecComm%NumHosts

 IF (pJD.ne.myRecComm%myNodeGroup) then
  if (myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0) then 
    iaux = size(myRC(ILEV)%CODECr(pJD)%x)/3-1
!   IF (iaux.gt.0) THEN
   
    DO j=1,iaux
     jaux = 3*(j-1)
     IF (myRC(ILEV)%CODECr(pJD)%x(jaux+2).eq.myid) THEN ! receiver is me
      pID  = myRC(ILEV)%CODECr(pJD)%x(jaux+1) ! sender
      iLoc = myRC(ILEV)%CODECr(pJD)%x(jaux+3) ! Start location for reading the data
      nSIZE = myRC(ILEV)%CODECr(pJD)%x(jaux+6) - myRC(ILEV)%CODECr(pJD)%x(jaux+3) !MGE013(ILEV)%ST(pID)%Num
  !     write(*,*) myid,pID,nSIZE,MGE013(ILEV)%ST(pID)%Num
      DO I=1,nSIZE
        iShift = MGE013(ILEV)%ST(pID)%VertLink(2,I)
        FX(iShift) = FX(iShift) + myRC(ILEV)%r(pJD)%x(iLoc + I)
      END DO
     END IF
    END DO
  END IF
  
 ELSE !pJD.eq.myRecComm%myNodeGroup)
 
  if (myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0) then 
   iaux = size(myRC(ILEV)%CODECs(pJD)%x)/3-1
!   IF (iaux.gt.0) THEN
   DO j=1,iaux
    jaux = 3*(j-1)
    IF (myRC(ILEV)%CODECs(pJD)%x(jaux+2).eq.myid) THEN ! receiver is me
     pID  = myRC(ILEV)%CODECs(pJD)%x(jaux+1) ! sender
     iLoc = myRC(ILEV)%CODECs(pJD)%x(jaux+3) ! Start location for reading the data
     nSIZE = myRC(ILEV)%CODECs(pJD)%x(jaux+6) - myRC(ILEV)%CODECs(pJD)%x(jaux+3) !MGE013(ILEV)%ST(pID)%Num
!     write(*,*) myid,pID,nSIZE,MGE013(ILEV)%ST(pID)%Num
     DO I=1,nSIZE
        iShift = MGE013(ILEV)%ST(pID)%VertLink(2,I)
        FX(iShift) = FX(iShift) + myRC(ILEV)%s(pJD)%x(iLoc + I)
     END DO
    END IF
   END DO
  END IF
   
 END IF
END DO

! CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

CALL ztime(tt1)
myStat%tCommV = myStat%tCommV + (tt1-tt0)
 
! write(*,*) myid,'done', (tt1-tt0)
 
END
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013UVWMAT_Rec(A11,A22,A33,KLDA,NU) !ok
use, INTRINSIC :: ISO_C_BINDING
use mpi
use shared_memory_module, only : get_shared_memory_INT,get_shared_memory_DBL,myRC
USE PP3D_MPI, ONLY : ierr,myid,master,numnodes,subnodes,MPI_COMM_SUBS,MPI_COMM_SUBGROUP
USE PP3D_MPI, ONLY : SENDD_myMPI,RECVD_myMPI,SENDK_myMPI,RECVK_myMPI,MGE013
USE def_feat, ONLY: ILEV
USE var_QuadScalar,ONLY:knvt,knet,knat,knel,myStat,nlmax
use var_QuadScalar, ONLY : myRecComm

implicit none
REAL*8 A11(*),A22(*),A33(*)
INTEGER KLDA(*),NU

INTEGER send_req(numnodes),recv_req(numnodes)
INTEGER STATUS(MPI_STATUS_SIZE)
INTEGER I,J,pID,pJD,nSIZE,nEIGH,iLOC,iSHIFT,iAUX,jAUX,nXX
INTEGER MEQ,MEQ1,MEQ2,MEQ3
INTEGER LEQ,LEQ1,LEQ2,LEQ3
INTEGER :: dblesize=8,intsize=4

character(len=256) :: cFMT
REAL*4 tt0,tt1

IF (myid.eq.MASTER) return

send_req = MPI_REQUEST_NULL
recv_req = MPI_REQUEST_NULL

if (.not.allocated(myRC)) allocate(myRC(1:nlmax))

LEQ = KNVT(ILEV) + KNAT(ILEV) + KNET(ILEV) + KNEL(ILEV)
LEQ1 =0
LEQ2 =LEQ
LEQ3 =2*LEQ
 
! CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
CALL ztime(tt0)

IF (.not.myRC(ILEV)%bPrepared) CALL E013RecCOMM_Init()

DO I=1,NU
 MGE013(ILEV)%UE11(I)=A11(KLDA(I))
 MGE013(ILEV)%UE22(I)=A22(KLDA(I))
 MGE013(ILEV)%UE33(I)=A33(KLDA(I))
ENDDO

LEQ = KNVT(ILEV) + KNAT(ILEV) + KNET(ILEV) + KNEL(ILEV)
LEQ1 =0
LEQ2 =LEQ
LEQ3 =2*LEQ

! extract the data to be sent
! All processes from the given host upload the data packages to be sent
DO pJD=1,myRecComm%NumHosts

 if (myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0) then 
  iaux = size(myRC(ILEV)%CODECs(pJD)%x)/3-1
!  IF (iaux.gt.0) THEN
  
   DO j=1,iaux
    jaux = 3*(j-1)
    IF (myRC(ILEV)%CODECs(pJD)%x(jaux+1).eq.myid) THEN
     pID  = myRC(ILEV)%CODECs(pJD)%x(jaux+2) ! receiver
     iLoc = myRC(ILEV)%CODECs(pJD)%x(jaux+3) ! Start location for reading the data
     nSIZE = myRC(ILEV)%CODECs(pJD)%x(jaux+6) - myRC(ILEV)%CODECs(pJD)%x(jaux+3) !MGE013(ILEV)%ST(pID)%Num
 !     write(*,*) myid,pID,nSIZE,MGE013(ILEV)%ST(pID)%Num
     DO I=1,nSIZE
       iShift = MGE013(ILEV)%ST(pID)%VertLink(1,I)
       myRC(ILEV)%s(pJD)%x(3*(iLoc + I - 1) + 1) = MGE013(ILEV)%UE11(iShift)
       myRC(ILEV)%s(pJD)%x(3*(iLoc + I - 1) + 2) = MGE013(ILEV)%UE22(iShift)
       myRC(ILEV)%s(pJD)%x(3*(iLoc + I - 1) + 3) = MGE013(ILEV)%UE33(iShift)
     END DO
    END IF
   END DO
 END IF
 
END DO

CALL MPI_BARRIER(MPI_COMM_SUBGROUP,IERR)

! Communication of the group-leaders
IF (myRecComm%myID.eq.0) THEN

 DO pID=1,myRecComm%NumHosts
 
  IF (myRecComm%myNodeGroup.NE.pID) THEN
  
   if (myRC(ILEV)%StartOfAllRecords(pID,myRecComm%NumNodes+1).gt.0) then 
     iaux = size(myRC(ILEV)%CODECs(pID)%x)/3
     nSIZE = myRC(ILEV)%CODECs(pID)%x(3*iaux)
!     CALL SENDD_myMPI(myRC(ILEV)%s(pID)%x,3*nSIZE,myRecComm%hostleaders(pID))
     CALL MPI_ISEND(myRC(ILEV)%s(pID)%x,3*nSIZE,MPI_DOUBLE_PRECISION,myRecComm%hostleaders(pID),1001,MPI_COMM_WORLD,send_req(pID),IERR)
    END IF
   
  ELSE
   DO pJD=1,myRecComm%NumHosts
   
    if (myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0.and.pJD.ne.myRecComm%myNodeGroup) then 
     iaux = size(myRC(ILEV)%CODECs(pJD)%x)/3
     nSIZE = myRC(ILEV)%CODECs(pJD)%x(3*iaux)
!     CALL RECVD_myMPI(myRC(ILEV)%r(pJD)%x,3*nSIZE,myRecComm%hostleaders(pJD))
     CALL MPI_IRECV(myRC(ILEV)%r(pJD)%x,3*nSIZE,MPI_DOUBLE_PRECISION,myRecComm%hostleaders(pJD),1001,MPI_COMM_WORLD,recv_req(pJD),IERR)
    END IF
   END DO
  END IF
  
 END DO
END IF

IF (myRecComm%myID.eq.0) THEN
 DO pID=1,myRecComm%NumHosts
  IF (myRC(ILEV)%StartOfAllRecords(pID,myRecComm%NumNodes+1).gt.0.and.myRecComm%myNodeGroup.NE.pID) then 
     CALL MPI_Wait(send_req(pID),STATUS, IERR )
     CALL MPI_Wait(recv_req(pID),STATUS, IERR )
  END IF
 END DO
END IF

CALL MPI_BARRIER(MPI_COMM_SUBGROUP,IERR)

! if (myRecComm%myid.eq.0) then 
!  DO pJD=1,myRecComm%NumHosts
!   DO i=0,size(myRC(ILEV)%CODECr(pJD)%x)/3-1
!    write(*,*) pJD,"X",myRC(ILEV)%CODECr(pJD)%x(3*i+1:3*i+3)
!   END DO
!  END DO
! end if
! 
! write(*,*) myid,'done'
! pause

! All processes from the given host extract their own data packages
DO pJD=1,myRecComm%NumHosts

 IF (pJD.ne.myRecComm%myNodeGroup) then
  if (myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0) then 
    iaux = size(myRC(ILEV)%CODECr(pJD)%x)/3-1
!   IF (iaux.gt.0) THEN
   
    DO j=1,iaux
     jaux = 3*(j-1)
     IF (myRC(ILEV)%CODECr(pJD)%x(jaux+2).eq.myid) THEN ! receiver is me
      pID  = myRC(ILEV)%CODECr(pJD)%x(jaux+1) ! sender
      iLoc = myRC(ILEV)%CODECr(pJD)%x(jaux+3) ! Start location for reading the data
      nSIZE = myRC(ILEV)%CODECr(pJD)%x(jaux+6) - myRC(ILEV)%CODECr(pJD)%x(jaux+3) !MGE013(ILEV)%ST(pID)%Num
  !     write(*,*) myid,pID,nSIZE,MGE013(ILEV)%ST(pID)%Num
      DO I=1,nSIZE
        iShift = MGE013(ILEV)%ST(pID)%VertLink(2,I)
        MGE013(ILEV)%UE11(iShift) = MGE013(ILEV)%UE11(iShift) + myRC(ILEV)%r(pJD)%x(3*(iLoc + I - 1) + 1)
        MGE013(ILEV)%UE22(iShift) = MGE013(ILEV)%UE22(iShift) + myRC(ILEV)%r(pJD)%x(3*(iLoc + I - 1) + 2)
        MGE013(ILEV)%UE33(iShift) = MGE013(ILEV)%UE33(iShift) + myRC(ILEV)%r(pJD)%x(3*(iLoc + I - 1) + 3)
      END DO
     END IF
    END DO
  END IF
  
 ELSE !pJD.eq.myRecComm%myNodeGroup)
 
  if (myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0) then 
   iaux = size(myRC(ILEV)%CODECs(pJD)%x)/3-1
!   IF (iaux.gt.0) THEN
   DO j=1,iaux
    jaux = 3*(j-1)
    IF (myRC(ILEV)%CODECs(pJD)%x(jaux+2).eq.myid) THEN ! receiver is me
     pID  = myRC(ILEV)%CODECs(pJD)%x(jaux+1) ! sender
     iLoc = myRC(ILEV)%CODECs(pJD)%x(jaux+3) ! Start location for reading the data
     nSIZE = myRC(ILEV)%CODECs(pJD)%x(jaux+6) - myRC(ILEV)%CODECs(pJD)%x(jaux+3) !MGE013(ILEV)%ST(pID)%Num
!     write(*,*) myid,pID,nSIZE,MGE013(ILEV)%ST(pID)%Num
     DO I=1,nSIZE
        iShift = MGE013(ILEV)%ST(pID)%VertLink(2,I)
        MGE013(ILEV)%UE11(iShift) = MGE013(ILEV)%UE11(iShift) + myRC(ILEV)%s(pJD)%x(3*(iLoc + I - 1) + 1)
        MGE013(ILEV)%UE22(iShift) = MGE013(ILEV)%UE22(iShift) + myRC(ILEV)%s(pJD)%x(3*(iLoc + I - 1) + 2)
        MGE013(ILEV)%UE33(iShift) = MGE013(ILEV)%UE33(iShift) + myRC(ILEV)%s(pJD)%x(3*(iLoc + I - 1) + 3)
     END DO
    END IF
   END DO
  END IF
   
 END IF
END DO

! CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

CALL ztime(tt1)
myStat%tCommV = myStat%tCommV + (tt1-tt0)

END SUBROUTINE E013UVWMAT_Rec
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E013MAT_Rec(A,KLDA,NU) !ok
! USE PP3D_MPI
! USE def_feat, ONLY: ILEV
! USE var_QuadScalar,ONLY:myStat

use, INTRINSIC :: ISO_C_BINDING
use mpi
use shared_memory_module, only : get_shared_memory_INT,get_shared_memory_DBL,myRC
USE PP3D_MPI, ONLY : ierr,myid,master,numnodes,subnodes,MPI_COMM_SUBS,MPI_COMM_SUBGROUP
USE PP3D_MPI, ONLY : SENDD_myMPI,RECVD_myMPI,SENDK_myMPI,RECVK_myMPI,MGE013
USE def_feat, ONLY: ILEV
USE var_QuadScalar,ONLY:knvt,knet,knat,knel,myStat,nlmax
use var_QuadScalar, ONLY : myRecComm

implicit none
REAL*8 A(*)
INTEGER KLDA(*),NU

INTEGER send_req(numnodes),recv_req(numnodes)
INTEGER STATUS(MPI_STATUS_SIZE)

INTEGER I,J,pID,pJD,nSIZE,nEIGH,iLOC,iSHIFT,iAUX,jAUX,nXX
INTEGER MEQ,MEQ1,MEQ2,MEQ3
INTEGER LEQ,LEQ1,LEQ2,LEQ3
integer, allocatable :: NumberOfMyRecords(:),NumberOfAllRecords(:,:)!,StartOfAllRecords(:,:)
! integer(kind=MPI_ADDRESS_KIND) :: win_size
INTEGER :: dblesize=8,intsize=4

character(len=256) :: cFMT
REAL*4 tt0,tt1

IF (myid.eq.MASTER) return

send_req = MPI_REQUEST_NULL
recv_req = MPI_REQUEST_NULL

if (.not.allocated(myRC)) allocate(myRC(1:nlmax))

IF (.not.myRC(ILEV)%bPrepared) CALL E013RecCOMM_Init()

CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
CALL ztime(tt0)

DO I=1,NU
 MGE013(ILEV)%UE(I)=A(KLDA(I))
ENDDO

! extract the data to be sent
! All processes from the given host upload the data packages to be sent
DO pJD=1,myRecComm%NumHosts

 if (myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0) then 
  iaux = size(myRC(ILEV)%CODECs(pJD)%x)/3-1
!  IF (iaux.gt.0) THEN
  
   DO j=1,iaux
    jaux = 3*(j-1)
    IF (myRC(ILEV)%CODECs(pJD)%x(jaux+1).eq.myid) THEN
     pID  = myRC(ILEV)%CODECs(pJD)%x(jaux+2) ! receiver
     iLoc = myRC(ILEV)%CODECs(pJD)%x(jaux+3) ! Start location for reading the data
     nSIZE = myRC(ILEV)%CODECs(pJD)%x(jaux+6) - myRC(ILEV)%CODECs(pJD)%x(jaux+3) !MGE013(ILEV)%ST(pID)%Num
 !     write(*,*) myid,pID,nSIZE,MGE013(ILEV)%ST(pID)%Num
     DO I=1,nSIZE
       iShift = MGE013(ILEV)%ST(pID)%VertLink(1,I)
       myRC(ILEV)%s(pJD)%x(iLoc + I ) = MGE013(ILEV)%UE(iShift)
!        myRC(ILEV)%s(pJD)%x(3*(iLoc + I - 1) + 2) = FX(LEQ2 + iShift)
!        myRC(ILEV)%s(pJD)%x(3*(iLoc + I - 1) + 3) = FX(LEQ3 + iShift)
     END DO
    END IF
   END DO
 END IF
 
END DO

CALL MPI_BARRIER(MPI_COMM_SUBGROUP,IERR)

! Communication of the group-leaders
IF (myRecComm%myID.eq.0) THEN

 DO pID=1,myRecComm%NumHosts
 
  IF (myRecComm%myNodeGroup.NE.pID) THEN
  
   if (myRC(ILEV)%StartOfAllRecords(pID,myRecComm%NumNodes+1).gt.0) then 
     iaux = size(myRC(ILEV)%CODECs(pID)%x)/3
     nSIZE = myRC(ILEV)%CODECs(pID)%x(3*iaux)
    CALL SENDD_myMPI(myRC(ILEV)%s(pID)%x,nSIZE,myRecComm%hostleaders(pID))
!      CALL MPI_ISEND(myRC(ILEV)%s(pID)%x,nSIZE,MPI_DOUBLE_PRECISION,myRecComm%hostleaders(pID),1001,MPI_COMM_WORLD,send_req(pID),IERR)
    END IF
   
  ELSE
   DO pJD=1,myRecComm%NumHosts
   
    if (myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0.and.pJD.ne.myRecComm%myNodeGroup) then 
     iaux = size(myRC(ILEV)%CODECs(pJD)%x)/3
     nSIZE = myRC(ILEV)%CODECs(pJD)%x(3*iaux)
    CALL RECVD_myMPI(myRC(ILEV)%r(pJD)%x,nSIZE,myRecComm%hostleaders(pJD))
!      CALL MPI_IRECV(myRC(ILEV)%r(pJD)%x,nSIZE,MPI_DOUBLE_PRECISION,myRecComm%hostleaders(pJD),1001,MPI_COMM_WORLD,recv_req(pJD),IERR)
    END IF
   END DO
  END IF
  
 END DO
END IF

! IF (myRecComm%myID.eq.0) THEN
!  DO pID=1,myRecComm%NumHosts
!   IF (myRC(ILEV)%StartOfAllRecords(pID,myRecComm%NumNodes+1).gt.0.and.myRecComm%myNodeGroup.NE.pID) then 
!      CALL MPI_Wait(send_req(pID),STATUS, IERR )
!      CALL MPI_Wait(recv_req(pID),STATUS, IERR )
!   END IF
!  END DO
! END IF

CALL MPI_BARRIER(MPI_COMM_SUBGROUP,IERR)

! if (myRecComm%myid.eq.0) then 
!  DO pJD=1,myRecComm%NumHosts
!   DO i=0,size(myRC(ILEV)%CODECr(pJD)%x)/3-1
!    write(*,*) pJD,"X",myRC(ILEV)%CODECr(pJD)%x(3*i+1:3*i+3)
!   END DO
!  END DO
! end if
! 
! write(*,*) myid,'done'
! pause

! All processes from the given host extract their own data packages
DO pJD=1,myRecComm%NumHosts

 IF (pJD.ne.myRecComm%myNodeGroup) then
  if (myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0) then 
    iaux = size(myRC(ILEV)%CODECr(pJD)%x)/3-1
!   IF (iaux.gt.0) THEN
   
    DO j=1,iaux
     jaux = 3*(j-1)
     IF (myRC(ILEV)%CODECr(pJD)%x(jaux+2).eq.myid) THEN ! receiver is me
      pID  = myRC(ILEV)%CODECr(pJD)%x(jaux+1) ! sender
      iLoc = myRC(ILEV)%CODECr(pJD)%x(jaux+3) ! Start location for reading the data
      nSIZE = myRC(ILEV)%CODECr(pJD)%x(jaux+6) - myRC(ILEV)%CODECr(pJD)%x(jaux+3) !MGE013(ILEV)%ST(pID)%Num
  !     write(*,*) myid,pID,nSIZE,MGE013(ILEV)%ST(pID)%Num
      DO I=1,nSIZE
        iShift = MGE013(ILEV)%ST(pID)%VertLink(2,I)
        MGE013(ILEV)%UE(iShift) = MGE013(ILEV)%UE(iShift) + myRC(ILEV)%r(pJD)%x(iLoc + I)
      END DO
     END IF
    END DO
  END IF
  
 ELSE !pJD.eq.myRecComm%myNodeGroup)
 
  if (myRC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0) then 
   iaux = size(myRC(ILEV)%CODECs(pJD)%x)/3-1
!   IF (iaux.gt.0) THEN
   DO j=1,iaux
    jaux = 3*(j-1)
    IF (myRC(ILEV)%CODECs(pJD)%x(jaux+2).eq.myid) THEN ! receiver is me
     pID  = myRC(ILEV)%CODECs(pJD)%x(jaux+1) ! sender
     iLoc = myRC(ILEV)%CODECs(pJD)%x(jaux+3) ! Start location for reading the data
     nSIZE = myRC(ILEV)%CODECs(pJD)%x(jaux+6) - myRC(ILEV)%CODECs(pJD)%x(jaux+3) !MGE013(ILEV)%ST(pID)%Num
!     write(*,*) myid,pID,nSIZE,MGE013(ILEV)%ST(pID)%Num
     DO I=1,nSIZE
        iShift = MGE013(ILEV)%ST(pID)%VertLink(2,I)
        MGE013(ILEV)%UE(iShift) = MGE013(ILEV)%UE(iShift) + myRC(ILEV)%s(pJD)%x(iLoc + I)
     END DO
    END IF
   END DO
  END IF
   
 END IF
END DO

! CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE E013MAT_Rec
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE GetParPressureRec(P,PP)
use, INTRINSIC :: ISO_C_BINDING
use mpi
use shared_memory_module, only : get_shared_memory_INT,get_shared_memory_DBL,myPC
USE PP3D_MPI, ONLY : ierr,myid,master,numnodes,subnodes,MPI_COMM_SUBS,MPI_COMM_SUBGROUP
USE PP3D_MPI, ONLY : SENDD_myMPI,RECVD_myMPI,SENDK_myMPI,RECVK_myMPI,MGE013,SENDI_myMPI,RECVI_myMPI
USE def_feat, ONLY: ILEV
USE var_QuadScalar,ONLY:knvt,knet,knat,knel,myStat,nlmax
use var_QuadScalar, ONLY : myRecComm

implicit none
REAL*8 P(*),PP(*)

INTEGER send_req(numnodes),recv_req(numnodes)
INTEGER STATUS(MPI_STATUS_SIZE)

INTEGER I,J,pID,pJD,nSIZE,nEIGH,iLOC,iSHIFT,iAUX,jAUX,nXX
INTEGER II,JJ
INTEGER :: dblesize=8,intsize=4

character(len=256) :: cFMT
REAL*4 tt0,tt1

IF (myid.eq.MASTER) return

send_req = MPI_REQUEST_NULL
recv_req = MPI_REQUEST_NULL

! IF (myid.eq.1) write(*,*) ilev, "communication..."

if (.not.allocated(myPC)) allocate(myPC(1:nlmax))
 
! CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
CALL ztime(tt0)

IF (.not.myPC(ILEV)%bPrepared) CALL E012RecCOMM_Init()

! extract the data to be sent
! All processes from the given host upload the data packages to be sent
DO pJD=1,myRecComm%NumHosts

 if (myPC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0) then 
  iaux = myPC(ILEV)%sendSIZE(pJD)-1
  
   DO j=1,iaux
    jaux = 3*(j-1)
    IF (myPC(ILEV)%CODECs(pJD)%x(jaux+1).eq.myid) THEN
     pID  = myPC(ILEV)%CODECs(pJD)%x(jaux+2) ! receiver
     iLoc = myPC(ILEV)%CODECs(pJD)%x(jaux+3) ! Start location for reading the data
     nSIZE = myPC(ILEV)%CODECs(pJD)%x(jaux+6) - myPC(ILEV)%CODECs(pJD)%x(jaux+3) !MGE013(ILEV)%SP(pID)%Num
     DO I=1,nSIZE
       JJ = 4*(MGE013(ILEV)%SP(pID)%VertLink(1,I)-1)+1
       myPC(ILEV)%s(pJD)%x(4*(iLoc + I - 1) + 1) = P(JJ+0)
       myPC(ILEV)%s(pJD)%x(4*(iLoc + I - 1) + 2) = P(JJ+1)
       myPC(ILEV)%s(pJD)%x(4*(iLoc + I - 1) + 3) = P(JJ+2)
       myPC(ILEV)%s(pJD)%x(4*(iLoc + I - 1) + 4) = P(JJ+3)
     END DO
    END IF
   END DO
 END IF
 
END DO

CALL MPI_BARRIER(MPI_COMM_SUBGROUP,IERR)

! Communication of the group-leaders
IF (myRecComm%myID.eq.0) THEN

 DO pID=1,myRecComm%NumHosts
 
  IF (myRecComm%myNodeGroup.NE.pID) THEN
  
   if (myPC(ILEV)%StartOfAllRecords(pID,myRecComm%NumNodes+1).gt.0) then 
     iaux = myPC(ILEV)%sendSIZE(pID)
     nSIZE = myPC(ILEV)%CODECs(pID)%x(3*iaux)
     
!      CALL SENDD_myMPI(myPC(ILEV)%s(pID)%x,3*nSIZE,myRecComm%hostleaders(pID))
     CALL MPI_ISEND(myPC(ILEV)%s(pID)%x,4*nSIZE,MPI_DOUBLE_PRECISION,myRecComm%hostleaders(pID),1001,MPI_COMM_WORLD,send_req(pID),IERR)
    END IF
   
  ELSE
   DO pJD=1,myRecComm%NumHosts
   
    if (myPC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0.and.pJD.ne.myRecComm%myNodeGroup) then 
     iaux = myPC(ILEV)%recvSIZE(pJD)
     nSIZE = myPC(ILEV)%CODECr(pJD)%x(3*iaux)
!      CALL RECVD_myMPI(myPC(ILEV)%r(pJD)%x,3*nSIZE,myRecComm%hostleaders(pJD))
     CALL MPI_IRECV(myPC(ILEV)%r(pJD)%x,4*nSIZE,MPI_DOUBLE_PRECISION,myRecComm%hostleaders(pJD),1001,MPI_COMM_WORLD,recv_req(pJD),IERR)
    END IF
   END DO
  END IF
  
 END DO
END IF

IF (myRecComm%myID.eq.0) THEN
 DO pID=1,myRecComm%NumHosts
  IF (myPC(ILEV)%StartOfAllRecords(pID,myRecComm%NumNodes+1).gt.0.and.myRecComm%myNodeGroup.NE.pID) then 
     CALL MPI_Wait(send_req(pID),STATUS, IERR )
     CALL MPI_Wait(recv_req(pID),STATUS, IERR )
  END IF
 END DO
END IF

CALL MPI_BARRIER(MPI_COMM_SUBGROUP,IERR)
! 
! ! if (myRecComm%myid.eq.0) then 
! !  DO pJD=1,myRecComm%NumHosts
! !   DO i=0,size(myPC(ILEV)%CODECr(pJD)%x)/3-1
! !    write(*,*) pJD,"X",myPC(ILEV)%CODECr(pJD)%x(3*i+1:3*i+3)
! !   END DO
! !  END DO
! ! end if
! ! 
! ! write(*,*) myid,'done'
! ! pause
! 
! All processes from the given host extract their own data packages
DO pJD=1,myRecComm%NumHosts

 IF (pJD.ne.myRecComm%myNodeGroup) then
  if (myPC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0) then 
    iaux = myPC(ILEV)%recvSIZE(pJD)-1
   
    DO j=1,iaux
     jaux = 3*(j-1)
     IF (myPC(ILEV)%CODECr(pJD)%x(jaux+2).eq.myid) THEN ! receiver is me
      pID  = myPC(ILEV)%CODECr(pJD)%x(jaux+1) ! sender
      iLoc = myPC(ILEV)%CODECr(pJD)%x(jaux+3) ! Start location for reading the data
      nSIZE = myPC(ILEV)%CODECr(pJD)%x(jaux+6) - myPC(ILEV)%CODECr(pJD)%x(jaux+3) !MGE013(ILEV)%SP(pID)%Num
!       write(*,*) myid,pID,nSIZE,MGE013(ILEV)%SP(pID)%nElems(2)
      DO I=1,nSIZE
        JJ = 4*(MGE013(ILEV)%SP(pID)%VertLink(2,I)-1)+1
        PP(JJ+0) = myPC(ILEV)%r(pJD)%x(4*(iLoc + I - 1) + 1)
        PP(JJ+1) = myPC(ILEV)%r(pJD)%x(4*(iLoc + I - 1) + 2)
        PP(JJ+2) = myPC(ILEV)%r(pJD)%x(4*(iLoc + I - 1) + 3)
        PP(JJ+3) = myPC(ILEV)%r(pJD)%x(4*(iLoc + I - 1) + 4)
      END DO
     END IF
    END DO
  END IF
  
 ELSE !pJD.eq.myRecComm%myNodeGroup)
 
  if (myPC(ILEV)%StartOfAllRecords(pJD,myRecComm%NumNodes+1).gt.0) then 
   iaux = myPC(ILEV)%sendSIZE(pJD)-1
!    write(*,*) myid,pJD,iaux
    
   DO j=1,iaux
    jaux = 3*(j-1)
    IF (myPC(ILEV)%CODECs(pJD)%x(jaux+2).eq.myid) THEN ! receiver is me
     pID  = myPC(ILEV)%CODECs(pJD)%x(jaux+1) ! sender
     iLoc = myPC(ILEV)%CODECs(pJD)%x(jaux+3) ! Start location for reading the data
     nSIZE = myPC(ILEV)%CODECs(pJD)%x(jaux+6) - myPC(ILEV)%CODECs(pJD)%x(jaux+3) !MGE013(ILEV)%SP(pID)%Num
!     write(*,*) myid,pID,nSIZE,MGE013(ILEV)%SP(pID)%Num
!       write(*,*) myid,pID,nSIZE,MGE013(ILEV)%SP(pID)%nElems(2)
     DO I=1,nSIZE
        JJ = 4*(MGE013(ILEV)%SP(pID)%VertLink(2,I)-1)+1
        PP(JJ+0) = myPC(ILEV)%s(pJD)%x(4*(iLoc + I - 1) + 1)
        PP(JJ+1) = myPC(ILEV)%s(pJD)%x(4*(iLoc + I - 1) + 2)
        PP(JJ+2) = myPC(ILEV)%s(pJD)%x(4*(iLoc + I - 1) + 3)
        PP(JJ+3) = myPC(ILEV)%s(pJD)%x(4*(iLoc + I - 1) + 4)
     END DO
    END IF
   END DO
  END IF
   
 END IF
END DO

! CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

CALL ztime(tt1)
myStat%tCommP = myStat%tCommP + (tt1-tt0)
 
! write(*,*) myid,'done', (tt1-tt0)
! pause
 
END

