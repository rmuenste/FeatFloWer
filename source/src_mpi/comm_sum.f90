! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE COMM_INLN_OLD(INLComplete)
USE PP3D_MPI, ONLY : ierr,myid,master,numnodes,subnodes
USE PP3D_MPI, ONLY : SENDI_myMPI,RECVI_myMPI

implicit none
INTEGER INLComplete
INTEGER iNL,piNL,pID

! write(*,*) "complete?",myid,INLComplete
IF (myid.eq.MASTER) THEN

  iNL=1
  DO pID=1,subnodes
  CALL RECVI_myMPI(piNL,pID)
  iNL=iNL*pINL
  END DO

  pINL=iNL
  DO pID=1,subnodes
  CALL SENDI_myMPI(piNL,pID)
  END DO
  INLComplete=iNL

ELSE
  piNL=INLComplete
  CALL SENDI_myMPI(piNL,0)
  CALL RECVI_myMPI(piNL,0)
  INLComplete=piNL
END IF

END SUBROUTINE COMM_INLN_OLD
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE COMM_INLN_NEW(INLComplete)
use mpi
USE PP3D_MPI, ONLY : ierr,myid,master
USE var_QuadScalar, ONLY : BaSynch

implicit none
INTEGER INLComplete
INTEGER INL,NN
INTEGER req
INTEGER STATUS(MPI_STATUS_SIZE)

req = MPI_REQUEST_NULL

if (myid.eq.master) INLComplete = 1
NN = 1

! Perform the asynchronous all-reduce operation to find the maximum value
IF (BaSynch) then
 call MPI_Iallreduce(INLComplete, INL, NN, MPI_INTEGER, MPI_PROD, MPI_COMM_WORLD, req, ierr)
 CALL MPI_Wait(req,STATUS, ierr )
ELSE
 call MPI_allreduce(INLComplete, INL, NN, MPI_INTEGER, MPI_PROD, MPI_COMM_WORLD, ierr)
END IF

INLComplete = INL

END SUBROUTINE COMM_INLN_NEW
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE COMM_MaximumN_OLD(DVAL,NN)
USE PP3D_MPI, ONLY : ierr,myid,master,numnodes,subnodes
USE PP3D_MPI, ONLY : SENDD_myMPI,RECVD_myMPI

implicit none
INTEGER NN
REAL*8 DVAL(NN)
REAL*8, ALLOCATABLE :: pVal(:)
INTEGER pID,i

ALLOCATE (pVal(NN))

IF (myid.eq.MASTER) THEN
  DVAL = -1d30
  DO pID=1,subnodes
  CALL RECVD_myMPI(pVAL,NN,pID)
  DO i=1,NN
  DVal(i) = MAX(DVAL(i),pVal(i))
  END DO
  END DO

  pVAL = DVAL
  DO pID=1,subnodes
  CALL SENDD_myMPI(pVAL,NN,pID)
  END DO

ELSE
  pVAL=DVAL
  CALL SENDD_myMPI(pVAL,NN,0)
  CALL RECVD_myMPI(pVAL,NN,0)
  DVAL = pVal
END IF

DEALLOCATE (pVal)

END SUBROUTINE COMM_MaximumN_OLD
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE COMM_MaximumN_NEW(DVAL,NN)
use mpi
USE PP3D_MPI, ONLY : ierr,myid,master,numnodes,subnodes
USE PP3D_MPI, ONLY : SENDD_myMPI,RECVD_myMPI
USE var_QuadScalar, ONLY : BaSynch

implicit none
INTEGER NN
REAL*8 DVAL(NN),DDVAL(NN)
INTEGER pID,i
INTEGER req
INTEGER STATUS(MPI_STATUS_SIZE)

req = MPI_REQUEST_NULL

if (myid.eq.master) DVAL(1:NN) = -1d99

! Perform the asynchronous all-reduce operation to find the maximum value
IF (BaSynch) then
 call MPI_Iallreduce(DVAL, DDVAL, NN, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, req, ierr)
 CALL MPI_Wait(req,STATUS, ierr )
ELSE
 call MPI_allreduce(DVAL, DDVAL, NN, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
END IF

DVAL = DDVAL

END SUBROUTINE COMM_MaximumN_NEW
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE COMM_SUMMN_OLD(DVAL,NN)
USE PP3D_MPI, ONLY : ierr,myid,master,numnodes,subnodes
USE PP3D_MPI, ONLY : SENDD_myMPI,RECVD_myMPI

implicit none
INTEGER NN
REAL*8 DVAL(NN)
REAL*8, ALLOCATABLE :: pVal(:)
INTEGER pID,i

ALLOCATE (pVal(NN))

IF (myid.eq.MASTER) THEN
  DVAL = 0d0
  DO pID=1,subnodes
  CALL RECVD_myMPI(pVAL,NN,pID)
  DVal = DVal + pVal
  END DO

  pVAL = DVAL
  DO pID=1,subnodes
  CALL SENDD_myMPI(pVAL,NN,pID)
  END DO

ELSE
  pVAL=DVAL
  CALL SENDD_myMPI(pVAL,NN,0)
  CALL RECVD_myMPI(pVAL,NN,0)
  DVAL = pVal
END IF

DEALLOCATE (pVal)

END SUBROUTINE COMM_SUMMN_OLD
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE COMM_SUMMN_NEW(DVAL,NN)
use mpi
USE PP3D_MPI, ONLY : ierr,myid,master,numnodes,subnodes
USE PP3D_MPI, ONLY : SENDD_myMPI,RECVD_myMPI
USE var_QuadScalar, ONLY : BaSynch

implicit none
INTEGER NN
REAL*8 DVAL(NN),DDVAL(NN)
INTEGER pID,i
INTEGER req
INTEGER STATUS(MPI_STATUS_SIZE)

req = MPI_REQUEST_NULL

if (myid.eq.master) DVAL(1:NN) = 0d0

! Perform the asynchronous all-reduce operation to find the maximum value
IF (BaSynch) then
 call MPI_Iallreduce(DVAL, DDVAL, NN, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, req, ierr)
 CALL MPI_Wait(req,STATUS, ierr )
ELSE
 call MPI_allreduce(DVAL, DDVAL, NN, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
END IF

DVAL = DDVAL

END SUBROUTINE COMM_SUMMN_NEW
