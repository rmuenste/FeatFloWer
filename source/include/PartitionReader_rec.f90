CMESH1="_mesh/                 "                     ! PARALLEL
LenFile = LEN((TRIM(ADJUSTL(cGridFileName))))
WRITE(CMESH1(7:7+LenFile),'(A,A1)') TRIM(ADJUSTL(cGridFileName)),'/'

IF (myid.ne.0) THEN                                  ! PARALLEL
  if (myRecComm%subIndex.le.0 .or. myRecComm%gridIndex.le.0) then
    write(*,*) 'Warning: missing sub/grid index for rank',myid
    stop 'FindNodes() did not populate the partition mapping.'
  end if

  WRITE(CMESH1(7+LenFile+1:7+LenFile+7+1),'(A3,I4.4,A1)') 'sub',myRecComm%subIndex,'/'  ! PARALLEL

  cProjectFolder = CMESH1
  WRITE(CMESH1(15+LenFile+1:15+11+LenFile+1),'(A4,I4.4,A4)') 'GRID',myRecComm%gridIndex,'.tri'  ! PARALLEL
  WRITE(cProjectNumber(1:4),'(I4.4)') myRecComm%gridIndex

ELSE                                                 ! PARALLEL
  cProjectFolder = CMESH1
  WRITE(CMESH1(7+LenFile+1:14+LenFile+1),'(A8)') 'GRID.tri'  ! PARALLEL
END IF                                               ! PARALLEL

CALL CommBarrier()
