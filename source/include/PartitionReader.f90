 CMESH1="_mesh/                 "                     ! PARALLEL
 LenFile = LEN((TRIM(ADJUSTL(cGridFileName))))
 WRITE(CMESH1(7:7+LenFile),'(A,A1)') TRIM(ADJUSTL(cGridFileName)),"/"
 IF (myid.ne.0) THEN                                  ! PARALLEL
   kSubPart = FLOOR(DBLE(subnodes)/DBLE(nSubCoarseMesh)-1d-10)+1
   iSubPart = FLOOR(DBLE(subnodes-myid+1)/DBLE(kSubPart)-1d-10)+1
   iPart    = (subnodes-myid+1) - (iSubPart-1)*kSubPart
   WRITE(CMESH1(7+LenFile+1:7+LenFile+7+1),'(A3,I4.4,A1)') "sub",iSubpart,"/"  ! PARALLEL

   cProjectFolder = CMESH1
   WRITE(CMESH1(15+LenFile+1:15+11+LenFile+1),'(A4,I4.4,A4)') "GRID",iPart,".tri"  ! PARALLEL
   WRITE(cProjectNumber(1:4),'(I4.4)') iPart
    
 ELSE                                                 ! PARALLEL
   cProjectFolder = CMESH1
   WRITE(CMESH1(7+LenFile+1:14+LenFile+1),'(A8)') "GRID.tri"  ! PARALLEL
 END IF                                               ! PARALLEL

 write(*,*) "Rank:",myid,"is assigned to the mesh:",ADJUSTL(TRIM(CMESH1))
 CALL CommBarrier()
