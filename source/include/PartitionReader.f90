! CMESH1="_mesh/                 "                     ! PARALLEL
! write(*,*)'length: ', len(CMESH1)
 LenFile = LEN((TRIM(ADJUSTL(cGridFileName))))
 CMESH1 = "_mesh/" // TRIM(ADJUSTL(cGridFileName)) // "/"
 write(*,*) 'mesh name: ', CMESH1
 !WRITE(CMESH1(7:7+LenFile),'(A,A1)') TRIM(ADJUSTL(cGridFileName)),"/"
 IF (myid.ne.0) THEN                                  ! PARALLEL
!   kSubPart = FLOOR(DBLE(subnodes)/DBLE(nSubCoarseMesh)-1d-10)+1
!   iSubPart = FLOOR(DBLE(subnodes-myid+1)/DBLE(kSubPart)-1d-10)+1
!   kSubPart = FLOOR(DBLE(subnodes)/DBLE(nSubCoarseMesh)-1d-10)+1
!   iSubPart = FLOOR(DBLE(myid)/DBLE(kSubPart)-1d-10)+1
!   iPart    = (subnodes-myid+1) - (iSubPart-1)*kSubPart
!   kSubPart = FLOOR(DBLE(subnodes)/DBLE(nSubCoarseMesh)-1d-10)+1
!   iSubPart = FLOOR(DBLE(myid)/DBLE(kSubPart)-1d-10)+1
!   iPart    = myid - (iSubPart-1)*kSubPart
   iPart = 1
   iSubpart = myid
   ctemp = "" 
   write(ctemp, '(A, I3.3, A)') "sub", iSubpart, "/"
   CMESH1 = CMESH1 // ctemp 

   WRITE(CMESH1(7+LenFile+1:7+LenFile+7+1),'(A3,I4.4,A1)') "sub",iSubpart,"/"  ! PARALLEL

   cProjectFolder = CMESH1
   WRITE(CMESH1(15+LenFile+1:15+11+LenFile+1),'(A4,I4.4,A4)') "GRID",iPart,".tri"  ! PARALLEL
   WRITE(cProjectNumber(1:4),'(I4.4)') iPart

   write(*,*)'Myid: ', myid , 'gets part: ',iPart, " sub: ", iSubpart, " c:",cmesh1 
    
 ELSE                                                 ! PARALLEL
   cProjectFolder = CMESH1
   WRITE(CMESH1(7+LenFile+1:14+LenFile+1),'(A8)') "GRID.tri"  ! PARALLEL
 END IF                                               ! PARALLEL

 write(*,*) "Rank:",myid,"is assigned to the mesh:",ADJUSTL(TRIM(CMESH1))
 CALL CommBarrier()
