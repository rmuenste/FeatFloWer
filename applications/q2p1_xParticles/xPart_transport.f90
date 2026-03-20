subroutine Transport_xParticles_MPI(log_unit)

use var_QuadScalar
USE PP3D_MPI, ONLY : myid,master,showid,subnodes,numnodes,MPI_COMM_SUBS,Barrier_myMPI
USE xPart_def, ONLY : xFactor,xChunks
use Sigma_User, only : myProcess
use omp_lib

implicit none

include 'mpif.h'

integer, intent(in) :: log_unit
integer nParticles,iChunk,jChunk,nChunks,ivt_min,ivt_max,pID,pJD,i,mParticles,iaux
integer iStatus(MPI_STATUS_SIZE)
INTEGER iErr,indice,n,Stats(6)
real*8 chunkTimeStart,chunkTimeEnd,chunkDuration,chunkTimeReport,chunkTime
real*8, allocatable :: daux(:,:),saux(:,:)
integer, allocatable :: iauxbuf(:)
character :: cCSV_File*(128)
integer, allocatable :: sendcounts(:),displs(:)
integer, allocatable :: sendcounts_idx(:),displs_idx(:)
real*8, allocatable :: gathered_data(:,:)
integer, allocatable :: gathered_indices(:)
real*8 BoundingBox(3,2),P(3),Q(3),dL(3),dist
logical bxOutput
integer ixOutput,iInflow,iInInflow

integer :: num_threads, thread_num,xivt_min,xivt_max,xnParticles
integer :: nLocalComplete



! Determine the start Package
ILEV = NLMAX
nParticles  = int(DBLE(mg_mesh%level(ILEV)%nel)*xFactor)
nChunks = xChunks*subnodes

! get the bounding box of the mesh
BoundingBox(3,1) =  1d8
BoundingBox(3,2) = -1d8

DO i=1,mg_mesh%level(ILEV)%nel
 P = mg_mesh%level(ILEV)%dcorvg(:,i)
 if (BoundingBox(1,1).gt.P(1)) BoundingBox(1,1) = P(1)
 if (BoundingBox(2,1).gt.P(2)) BoundingBox(2,1) = P(2)
 if (BoundingBox(3,1).gt.P(3)) BoundingBox(3,1) = P(3)
 if (BoundingBox(1,2).lt.P(1)) BoundingBox(1,2) = P(1)
 if (BoundingBox(2,2).lt.P(2)) BoundingBox(2,2) = P(2)
 if (BoundingBox(3,2).lt.P(3)) BoundingBox(3,2) = P(3)
END DO

dL = (BoundingBox(:,2) - BoundingBox(:,1))
BoundingBox(:,1) = BoundingBox(:,1) - 0.05*dL
BoundingBox(:,2) = BoundingBox(:,2) + 0.05*dL
if (myid.eq.1) then
 write(*,*) "The bounding box is : "
 write(*,'(A,2ES12.4)') "x: ", BoundingBox(1,:)
 write(*,'(A,2ES12.4)') "y: ", BoundingBox(2,:)
 write(*,'(A,2ES12.4)') "z: ", BoundingBox(3,:)
end if

 call PrepareDistanceToInflowOrdering()

! initialization of LostSet and CompleteSet and process the first data Chunk
if (myid.ne.master) then
 chunkTimeReport = 0d0
 ALLOCATE(myLostSet(1:nParticles))
 nLostSet = nParticles
 DO indice=1,nParticles
  myLostSet(indice)%coor   = 0d0
  myLostSet(indice)%time   = 0d0
  myLostSet(indice)%indice = 0
  myLostSet(indice)%id     = 0
  myLostSet(indice)%velo   = 0d0
 END DO
 
 ALLOCATE(myCompleteSet(1:nParticles))
 nCompleteSet = 0
 DO indice=1,nParticles
  myCompleteSet(indice)%coor   = 0d0
  myCompleteSet(indice)%time   = 0d0
  myCompleteSet(indice)%indice = 0
  myCompleteSet(indice)%id     = 0
  myCompleteSet(indice)%velo   = 0d0
 END DO

 iChunk = myid
 ivt_min = (iChunk-1)*(nParticles/nChunks) + 1
 ivt_max = (iChunk+0)*(nParticles/nChunks) + 0
 if (iChunk.eq.nChunks) ivt_max = nParticles
 
 Stats(6) = iChunk
 chunkTimeStart = MPI_WTIME()
 CALL Transport_xParticles(log_unit,iChunk,ivt_min,ivt_max,Stats)
 chunkTimeEnd = MPI_WTIME()
 chunkDuration = chunkTimeEnd - chunkTimeStart
 chunkTimeReport = chunkDuration
 mParticles = myParticleParam%nParticles
 
end if

if (myid.eq.0) THEN 
 write(*,'(9A10)') " -- ","Chunk","pID","#ALL","#ALLinFD","#Active","#Lost","Time[s]","Progress"
end if

CALL Barrier_myMPI()

! master is distributing the data Chunks
if (myid.eq.0) THEN 
 jChunk = 0
 DO iChunk = subnodes + 1, nChunks
  jChunk = jChunk + 1
  CALL MPI_RECV(Stats,6,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,iStatus,iErr)
  pID = Stats(1)
  CALL MPI_RECV(chunkTime,1,MPI_DOUBLE_PRECISION,pID,1,MPI_COMM_WORLD,iStatus,iErr)
  write(*,'(A10,6I10,F10.3,F8.1,A)') "done:",Stats(6),Stats(1:5),chunkTime,1d2*dble(jChunk)/dble(nChunks),"%"
  
  CALL MPI_SEND(iChunk,1,MPI_INT,pID,0,MPI_COMM_WORLD,iStatus,iErr)
 END DO
 
 DO pJD=1,subnodes
  iChunk = 0
  jChunk = jChunk + 1
  CALL MPI_RECV(Stats,6,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,iStatus,iErr)
  pID = Stats(1)
  CALL MPI_RECV(chunkTime,1,MPI_DOUBLE_PRECISION,pID,1,MPI_COMM_WORLD,iStatus,iErr)
  write(*,'(A10,6I10,F10.3,F8.1,A)') "done:",Stats(6),Stats(1:5),chunkTime,1d2*dble(jChunk)/dble(nChunks),"%"
  
  CALL MPI_SEND(iChunk,1,MPI_INT,pID,0,MPI_COMM_WORLD,iStatus,iErr)
 END DO
 
END IF

! workers receive the new data Chunks and process them
if (myid.ne.0) then
 DO 
  CALL MPI_SEND(Stats,6,MPI_INT,0,0,MPI_COMM_WORLD,iStatus,iErr)
  CALL MPI_SEND(chunkTimeReport,1,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,iStatus,iErr)
  CALL MPI_RECV(iChunk,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,iStatus,iErr)
  Stats(6) = iChunk
  if (iChunk.eq.0) GOTO 1
  
  ivt_min = (iChunk-1)*(nParticles/nChunks) + 1
  ivt_max = (iChunk+0)*(nParticles/nChunks) + 0
  if (iChunk.eq.nChunks) ivt_max = nParticles
  chunkTimeStart = MPI_WTIME()
  CALL Transport_xParticles(log_unit,iChunk,ivt_min,ivt_max,Stats)
  chunkTimeEnd = MPI_WTIME()
  chunkDuration = chunkTimeEnd - chunkTimeStart
  chunkTimeReport = chunkDuration
  mParticles = mParticles + myParticleParam%nParticles
 end do
end if

1 continue
! all data Chunks are processed 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Gather the Lost data set to see the partcles at the outflow 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (myid.ne.0) then

 call MPI_AllReduce(mParticles, iaux, 1, MPI_INT, MPI_SUM, MPI_COMM_SUBS, ierr)
 mParticles = iaux
 
 ALLOCATE(daux(3,nParticles),saux(3,nParticles))
 
 DO indice=1,nParticles
  daux(:,indice) = myLostSet(indice)%coor
 END DO
 call MPI_AllReduce(daux, saux, 3*nParticles, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
 DO indice=1,nParticles
  myLostSet(indice)%coor = saux(:,indice)
 END DO

 DO indice=1,nParticles
  daux(:,indice) = myLostSet(indice)%velo
 END DO
 call MPI_AllReduce(daux, saux, 3*nParticles, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
 DO indice=1,nParticles
  myLostSet(indice)%velo = saux(:,indice)
 END DO

 DO indice=1,nParticles
  daux(:,indice) = [dble(myLostSet(indice)%ID),dble(myLostSet(indice)%indice),myLostSet(indice)%time]
 END DO
 call MPI_AllReduce(daux, saux, 3*nParticles, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
 DO indice=1,nParticles
  myLostSet(indice)%ID     = int(saux(1,indice))
  myLostSet(indice)%indice = int(saux(2,indice))
  myLostSet(indice)%time   = saux(3,indice)
 END DO
 DeALLOCATE(daux,saux)
end if
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Lets collect all the particles being still inside of the FD  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 nLocalComplete = nCompleteSet
 allocate(sendcounts(0:numnodes),displs(0:numnodes+1))
 allocate(sendcounts_idx(0:numnodes),displs_idx(0:numnodes+1))
 sendcounts = 0; displs = 0
 sendcounts_idx = 0; displs_idx = 0
 
 call MPI_allgather(3*nLocalComplete, 1, MPI_INTEGER, sendcounts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

 displs = 0
 do i = 2, numnodes+1
   displs(i) = displs(i-1) + sendcounts(i-1)
 end do
 do i=0,numnodes
   sendcounts_idx(i) = sendcounts(i)/3
 end do
 displs_idx = 0
 do i = 2, numnodes+1
   displs_idx(i) = displs_idx(i-1) + sendcounts_idx(i-1)
 end do
 
 ALLOCATE(daux(3,max(1,nLocalComplete)))
 ALLOCATE(iauxbuf(max(1,nLocalComplete)))
 if (nLocalComplete.gt.0) then
  DO i=1,nLocalComplete
   daux(:,i) = myCompleteSet(i)%coor
   iauxbuf(i) = myCompleteSet(i)%indice
  END DO
 else
  daux = 0d0
  iauxbuf = 0
 end if
 
 if (myid.eq.master) then
   allocate(gathered_data(3, max(1,displs(numnodes+1)/3)))
   allocate(gathered_indices(max(1,displs_idx(numnodes+1))))
 else
   allocate(gathered_data(3,1))
   allocate(gathered_indices(1))
 end if
 
 call MPI_Gatherv(daux, 3*nLocalComplete, MPI_DOUBLE_PRECISION, &
                  gathered_data, sendcounts, displs, &
                  MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
 call MPI_Gatherv(iauxbuf, nLocalComplete, MPI_INTEGER, &
                  gathered_indices, sendcounts_idx, displs_idx, &
                  MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
                  
 
 call MPI_AllReduce(nLocalComplete, iaux, 1, MPI_INT, MPI_SUM, MPI_COMM_World, ierr)
 nCompleteSet = iaux
 
if (myid.eq.master) then
 ALLOCATE(myCompleteSet(1:max(1,nCompleteSet)))
 DO i=1,nCompleteSet
  myCompleteSet(i)%coor = gathered_data(:,i)
  myCompleteSet(i)%indice = gathered_indices(i)
 END DO
endif
if (allocated(daux))      deallocate(daux)
if (allocated(iauxbuf))   deallocate(iauxbuf)
if (allocated(gathered_data)) deallocate(gathered_data)
if (allocated(gathered_indices)) deallocate(gathered_indices)
if (allocated(sendcounts)) deallocate(sendcounts)
if (allocated(displs)) deallocate(displs)
if (allocated(sendcounts_idx)) deallocate(sendcounts_idx)
if (allocated(displs_idx)) deallocate(displs_idx)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Release the remaining particles 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (myid.eq.master) then
  write(cCSV_File,'(A)') 'IN.csv'
  OPEN(FILE=adjustl(trim(cCSV_File)),UNIT=11)
  WRITE(11,*) 'x,y,z,indice'
  DO i=1,nCompleteSet
    WRITE(11,'(3(ES12.4,(",")),1(I0,(",")),I0)') myCompleteSet(i)%coor,i
  END DO
  write(*,*) "Overall number of remaining particles: ", nCompleteSet
  CLOSE(11)
  if (nCompleteSet.gt.0) then
   call WriteRemainingParticlesVTU('IN_particles.vtu',myCompleteSet,nCompleteSet,&
        mg_mesh%level(ILEV)%dcorvg,mg_mesh%level(ILEV)%nvt+mg_mesh%level(ILEV)%net+&
        mg_mesh%level(ILEV)%nat+mg_mesh%level(ILEV)%nel,&
        mg_mesh%level(ILEV)%kvert,mg_mesh%level(ILEV)%nel)
  end if
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Release the outflow particles 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (myid.ne.0) then
 if (myid.eq.1) then
  write(*,*) "Overall number of traced particles: ", mParticles
  iaux = 0
  ixOutput = 0
  write(cCSV_File,'(A,2(i4.4,A))') 'ValidOutflow.csv'
  OPEN(FILE=adjustl(trim(cCSV_File)),UNIT=11)
  WRITE(11,*) 'x,y,z,vx,vy,vz,t,id,indice'
  write(cCSV_File,'(A,2(i4.4,A))') 'InvalidOutflow.csv'
  OPEN(FILE=adjustl(trim(cCSV_File)),UNIT=12)
  WRITE(12,*) 'x,y,z,vx,vy,vz,t,id,indice'
  DO i=1,nParticles
   if (myLostSet(i)%indice.gt.0) then
    iaux = iaux + 1
    P = myLostSet(i)%coor
    bxOutput = .TRUE.
    if (BoundingBox(1,1).gt.P(1)) bxOutput = .FALSE.
    if (BoundingBox(2,1).gt.P(2)) bxOutput = .FALSE.
    if (BoundingBox(3,1).gt.P(3)) bxOutput = .FALSE.
    if (BoundingBox(1,2).lt.P(1)) bxOutput = .FALSE.
    if (BoundingBox(2,2).lt.P(2)) bxOutput = .FALSE.
    if (BoundingBox(3,2).lt.P(3)) bxOutput = .FALSE.
    IF (bxOutput) THEN   
     ixOutput = ixOutput + 1  
     WRITE(11,'(7(ES12.4,(",")),1(I0,(",")),I0)') myLostSet(i)%coor,myLostSet(i)%velo,myLostSet(i)%time,myLostSet(i)%id,myLostSet(i)%indice
    ELSE
     WRITE(12,'(7(ES12.4,(",")),1(I0,(",")),I0)') myLostSet(i)%coor,myLostSet(i)%velo,myLostSet(i)%time,myLostSet(i)%id,myLostSet(i)%indice
    end if
   end if
  END DO
  CLOSE(11)
  write(*,*) "Overall number of outflow particles: ", iaux
  write(*,*) "Overall number of VALID outflow particles: ", ixOutput
  write(*,*) "Overall number of INVALID outflow particles: ", iaux-ixOutput
 end if
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if (myid.ne.0) then
 
  iInInflow = 0
  allocate(MaterialDistribution(1:NLMAX))
  allocate(MaterialDistribution(NLMAX)%x(nParticles))
  
  DO i=1,nParticles
   if (myLostSet(i)%indice.gt.0) then
    P = myLostSet(i)%coor
    bxOutput = .TRUE.
    if (BoundingBox(1,1).gt.P(1)) bxOutput = .FALSE.
    if (BoundingBox(2,1).gt.P(2)) bxOutput = .FALSE.
    if (BoundingBox(3,1).gt.P(3)) bxOutput = .FALSE.
    if (BoundingBox(1,2).lt.P(1)) bxOutput = .FALSE.
    if (BoundingBox(2,2).lt.P(2)) bxOutput = .FALSE.
    if (BoundingBox(3,2).lt.P(3)) bxOutput = .FALSE.
    MaterialDistribution(NLMAX)%x(i) = -myLostSet(i)%id
    IF (bxOutput) THEN
     DO iInflow=1,myProcess%nOfInflows
      Q = myProcess%myInflow(iInflow)%Center
      dist = sqrt((P(1)-Q(1))**2d0 + (P(2)-Q(2))**2d0 + (P(3)-Q(3))**2d0)
      IF (dist.le.myProcess%myInflow(iInflow)%outerradius*1.05d0) THEN
       iInInflow = iInInflow + 1
       MaterialDistribution(NLMAX)%x(i) = myProcess%myInflow(iInflow)%Material
      END IF
     END DO
    end if
   else
    MaterialDistribution(NLMAX)%x(i) = myLostSet(i)%id
   end if
  END DO
  
!  write(*,*) "Overall number of identified outflow particles(",myid,"): ", iInInflow

 end if
 
 END subroutine Transport_xParticles_MPI
!
! ------------------------------------------------------------------------------
!
subroutine Transport_xParticles(log_unit,iChunk,ivt_min,ivt_max,Stats)
    
USE def_FEAT
USE types
use var_QuadScalar
USE PP3D_MPI, ONLY : myid,master,showid
USE var_QuadScalar, ONLY : mg_Mesh,QuadSc,myExport,Shell
USE OctTreeSearch
USE xPart_def, ONLY : d_CorrDist,dTimeStep,minDist,nTime
use omp_lib

implicit none
integer, intent(in) :: log_unit,iChunk,ivt_min,ivt_max
integer :: Stats(6)

real*8, allocatable, Dimension(:) :: xField,yField,zField
integer ndof,indice,i,iParticle,iTime,ld00,ld11
!!!!!!!!!!1
REAL*8 :: dTime,dStart
character :: cCSV_File*(128)
logical :: bCSV=.false.

integer iMinParticleIndex,iMaxParticleIndex,thread_num,num_threads,nStillActiveSet
integer, allocatable :: iStillActiveSet(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!! USER params !!!!!!!!!!!!!!!!!!!!!!!!
myParticleParam%OutflowZPos = 1d9
!!!!!!!!!!!!!!!!!!!!!!!!!!!! USER params !!!!!!!!!!!!!!!!!!!!!!!!

ILEV = NLMAX
CALL SETLEV(2)

iParticle = 0
DO i=ivt_min,ivt_max
 indice = mg_mesh%level(ILEV)%nvt + mg_mesh%level(ILEV)%net + mg_mesh%level(ILEV)%nat +  i 
 if (Shell(indice).gt.minDist) then
  iParticle = iParticle + 1
 end if
END DO
myParticleParam%nParticles = iParticle

ALLOCATE(myActiveSet  (1:myParticleParam%nParticles))

! write(*,*) myid,myParticleParam%nParticles,mg_mesh%level(ILEV)%nvt,mg_mesh%level(ILEV)%net, mg_mesh%level(ILEV)%nat,mg_mesh%level(ILEV)%nel
! pause

iParticle = 0
DO i=ivt_min,ivt_max
 indice = mg_mesh%level(ILEV)%nvt + mg_mesh%level(ILEV)%net + mg_mesh%level(ILEV)%nat +  i 
 if (Shell(indice).gt.minDist) then
  iParticle = iParticle + 1
  myActiveSet(iParticle)%coor(:) = mg_mesh%level(ILEV)%dcorvg(:,indice)
  myActiveSet(iParticle)%time    = 0d0
  myActiveSet(iParticle)%id      = 0
  myActiveSet(iParticle)%indice  = i
  myActiveSet(iParticle)%velo    = 0d0
 else 
  myLostSet(i)%id       = -myid
  myLostSet(i)%indice   = -myid
 end if
END DO

if (bCSV) THEN
 write(cCSV_File,'(A,2(i4.4,A))') "P_",iChunk,"__",0,'.csv'
 OPEN(FILE=adjustl(trim(cCSV_File)),UNIT=11)
 WRITE(11,*) 'x,y,z,vx,vy,vz,t,id,indice'
 DO i=1,myParticleParam%nParticles
  WRITE(11,'(7(ES12.4,(",")),1(I0,(",")),I0)') myActiveSet(i)%coor,myActiveSet(i)%velo,myActiveSet(i)%time,myActiveSet(i)%id,myActiveSet(i)%indice
 END DO
 CLOSE(11)
END IF
 
nActiveSet   = myParticleParam%nParticles
nLostSet = 0

CALL InitOctTree(mg_mesh%level(ILEV)%dcorvg,mg_mesh%level(ILEV)%nvt)

!$omp parallel private(thread_num, iMinParticleIndex ,iMaxParticleIndex,dStart,dTime,iTime)

dTime     = 0d0

num_threads = omp_get_num_threads()  ! Get the total number of threads
thread_num = omp_get_thread_num()    ! Get the thread number for the current thread
if (thread_num.eq.0) then
 allocate (iStillActiveSet(0:num_threads-1))
 iStillActiveSet = 0
end if
!$omp barrier
 
DO iTime = 1,nTime
 dStart = dTime
 dTime  = dTime + dTimeStep
 
 iMinParticleIndex = (thread_num  )*(nActiveSet/num_threads) + 1
 iMaxParticleIndex = (thread_num+1)*(nActiveSet/num_threads) + 0
 if (thread_num.eq.num_threads-1) iMaxParticleIndex = nActiveSet
!  print *, myid,"Thread number:", thread_num,iMinParticleIndex,iMaxParticleIndex,nActiveSet
 
 CALL Move_xParticle_Fast(mg_mesh%level(ILEV)%dcorvg,&
                          mg_mesh%level(ILEV)%kvert,&
                          mg_mesh%level(ILEV)%kedge,&
                          mg_mesh%level(ILEV)%karea,&
                          mg_mesh%level(ILEV)%elementsAtVertexIdx,&
                          mg_mesh%level(ILEV)%elementsAtVertex,&
                          QuadSc%ValU,QuadSc%ValV,QuadSc%ValW,&
                          QuadSc%ValU,QuadSc%ValV,QuadSc%ValW,&
                          Shell,&
                          mg_mesh%level(ILEV)%nvt,&
                          mg_mesh%level(ILEV)%net,&
                          mg_mesh%level(ILEV)%nat,&
                          mg_mesh%level(ILEV)%nel,&
                          dTime,dStart,&
                          iMinParticleIndex,iMaxParticleIndex,&
                          iStillActiveSet(thread_num))
 
!$omp barrier
if (thread_num.eq.0) then
 nStillActiveSet = 0
 do i=0,num_threads-1
  nStillActiveSet = nStillActiveSet + iStillActiveSet(i)
 end do
end if

!  if (iChunk.eq.1) write(*,*) thread_num,iStillActiveSet(thread_num),iMinParticleIndex,iMaxParticleIndex
!$omp barrier
 if (nStillActiveSet.eq.0) exit
 
!  if (nActiveSet.eq.0) EXIT
 
!  if (bCSV) THEN
!   write(cCSV_File,'(A,2(i4.4,A))') "P_",iChunk,"__",iTime,'.csv'
!   OPEN(FILE=adjustl(trim(cCSV_File)),UNIT=11)
!   WRITE(11,*) 'x,y,z,vx,vy,vz,t,id,indice'
!   DO i=1,nActiveSet
!    WRITE(11,'(7(ES12.4,(",")),1(I0,(",")),I0)') myActiveSet(i)%coor,myActiveSet(i)%velo,myActiveSet(i)%time,myActiveSet(i)%id,myActiveSet(i)%indice
!   END DO
!   CLOSE(11)
!  END IF
 
END DO                            

!$omp end parallel

! pause
!write(*,'(A,I8,A,1000I8)') "myid: ",myid, " :: ", iStillActiveSet

Stats(1) = myid                           !! signaure
Stats(2) = ivt_max-ivt_min + 1            !! all potential particles
Stats(3) = iParticle                      !! all particles in FD
Stats(4) = nStillActiveSet                     !! all particles remained in FD
Stats(5) = iParticle - nStillActiveSet         !! all particles which got lost

if (nStillActiveSet.ne.0) THEN
 DO i=1,nActiveSet
  if (myActiveSet(i)%id.eq.0) then
   nCompleteSet = nCompleteSet + 1
   myCompleteSet(nCompleteSet)%coor   = myActiveSet(i)%coor
   myCompleteSet(nCompleteSet)%time   = myActiveSet(i)%time
   myCompleteSet(nCompleteSet)%indice = myActiveSet(i)%indice
   myCompleteSet(nCompleteSet)%id     = myid
   myCompleteSet(nCompleteSet)%velo   = myActiveSet(i)%velo
  
   myLostSet(myActiveSet(i)%indice)%id       = -myid
   myLostSet(myActiveSet(i)%indice)%indice   = -myid
  end if
  
 END DO
END IF

CALL FreeOctTree()  

! write(*,*) "nCompleteSet = ", nCompleteSet, " : ", myid,nActiveSet,nStillActiveSet
deALLOCATE(myActiveSet)

end subroutine Transport_xParticles

!========================================================================================
!                           Sub: WriteRemainingParticlesVTU
!========================================================================================
subroutine WriteRemainingParticlesVTU(filename,particles,nParticles,dcorvg,nCoord,kvert,nElem)
 USE types, ONLY : tParticle
 implicit none
 CHARACTER(*), intent(in) :: filename
 TYPE(tParticle), intent(in) :: particles(*)
 integer, intent(in) :: nParticles
 REAL*8, intent(in) :: dcorvg(3,*)
 integer, intent(in) :: nCoord
 integer, intent(in) :: kvert(8,*)
 integer, intent(in) :: nElem
 character(len=512) :: baseName,pointsFile,hexFile
 integer :: lenName

 if (nParticles.le.0) return

 baseName = ADJUSTL(TRIM(filename))
 lenName = LEN_TRIM(baseName)
 if (lenName.gt.4 .and. baseName(lenName-3:lenName).eq.'.vtu') then
  pointsFile = baseName(1:lenName-4)//'_points.vtu'
  hexFile    = baseName(1:lenName-4)//'_hex.vtu'
 else
  pointsFile = baseName//'_points.vtu'
  hexFile    = baseName//'_hex.vtu'
 end if

 call WriteParticlePoints(pointsFile,particles,nParticles)
 call WriteParticleHex(hexFile,particles,nParticles,dcorvg,nCoord,kvert,nElem)

contains

 subroutine WriteParticlePoints(cFile,particles,nParticles)
  character(*), intent(in) :: cFile
  TYPE(tParticle), intent(in) :: particles(*)
  integer, intent(in) :: nParticles
  integer :: unitVTU,i,offset
  character(len=32) :: cTmpPoints,cTmpCells

  WRITE(cTmpPoints,'(I0)') nParticles
  WRITE(cTmpCells ,'(I0)') nParticles

  open(newunit=unitVTU,file=ADJUSTL(TRIM(cFile)),status='replace',action='write')
  write(unitVTU,'(A)') '<?xml version="1.0"?>'
  write(unitVTU,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
  write(unitVTU,'(A)') '  <UnstructuredGrid>'
  write(unitVTU,'(A)') '    <Piece NumberOfPoints="'//TRIM(cTmpPoints)//'" NumberOfCells="'//TRIM(cTmpCells)//'">'

  write(unitVTU,'(A)') '      <Points>'
  write(unitVTU,'(A)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
  do i=1,nParticles
   write(unitVTU,'(3(1X,ES23.16))') particles(i)%coor
  end do
  write(unitVTU,'(A)') '        </DataArray>'
  write(unitVTU,'(A)') '      </Points>'

  write(unitVTU,'(A)') '      <Cells>'
  write(unitVTU,'(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
  do i=1,nParticles
   write(unitVTU,'(1X,I0)') i-1
  end do
  write(unitVTU,'(A)') '        </DataArray>'

  write(unitVTU,'(A)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
  offset = 0
  do i=1,nParticles
   offset = offset + 1
   write(unitVTU,'(1X,I0)') offset
  end do
  write(unitVTU,'(A)') '        </DataArray>'

  write(unitVTU,'(A)') '        <DataArray type="UInt8" Name="types" format="ascii">'
  do i=1,nParticles
   write(unitVTU,'(1X,I0)') 1
  end do
  write(unitVTU,'(A)') '        </DataArray>'
  write(unitVTU,'(A)') '      </Cells>'

  write(unitVTU,'(A)') '      <CellData Scalars="pair_id">'
  write(unitVTU,'(A)') '        <DataArray type="Int32" Name="pair_id" format="ascii">'
  do i=1,nParticles
   write(unitVTU,'(1X,I0)') i
  end do
  write(unitVTU,'(A)') '        </DataArray>'
  write(unitVTU,'(A)') '        <DataArray type="Int32" Name="start_element" format="ascii">'
  do i=1,nParticles
   write(unitVTU,'(1X,I0)') particles(i)%indice
  end do
  write(unitVTU,'(A)') '        </DataArray>'
  write(unitVTU,'(A)') '      </CellData>'

  write(unitVTU,'(A)') '    </Piece>'
  write(unitVTU,'(A)') '  </UnstructuredGrid>'
  write(unitVTU,'(A)') '</VTKFile>'
  close(unitVTU)
 end subroutine WriteParticlePoints

 subroutine WriteParticleHex(cFile,particles,nParticles,dcorvg,nCoord,kvert,nElem)
  character(*), intent(in) :: cFile
  TYPE(tParticle), intent(in) :: particles(*)
  integer, intent(in) :: nParticles
  REAL*8, intent(in) :: dcorvg(3,*)
  integer, intent(in) :: nCoord
  integer, intent(in) :: kvert(8,*)
  integer, intent(in) :: nElem
  integer :: unitVTU,i,j,elem,offset,pointBase,ivt
  REAL*8 :: coords(3)
  logical :: hasHex
  character(len=32) :: cTmpPoints,cTmpCells

  WRITE(cTmpPoints,'(I0)') 8*nParticles
  WRITE(cTmpCells ,'(I0)') nParticles

  open(newunit=unitVTU,file=ADJUSTL(TRIM(cFile)),status='replace',action='write')
  write(unitVTU,'(A)') '<?xml version="1.0"?>'
  write(unitVTU,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
  write(unitVTU,'(A)') '  <UnstructuredGrid>'
  write(unitVTU,'(A)') '    <Piece NumberOfPoints="'//TRIM(cTmpPoints)//'" NumberOfCells="'//TRIM(cTmpCells)//'">'

  write(unitVTU,'(A)') '      <Points>'
  write(unitVTU,'(A)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
  do i=1,nParticles
   elem = particles(i)%indice
   hasHex = elem.ge.1 .and. elem.le.nElem
   do j=1,8
    if (hasHex) then
     ivt = kvert(j,elem)
     if (ivt.ge.1 .and. ivt.le.nCoord) then
      coords = dcorvg(:,ivt)
     else
      coords = particles(i)%coor
     end if
    else
     coords = particles(i)%coor
    end if
    write(unitVTU,'(3(1X,ES23.16))') coords
   end do
  end do
  write(unitVTU,'(A)') '        </DataArray>'
  write(unitVTU,'(A)') '      </Points>'

  write(unitVTU,'(A)') '      <Cells>'
  write(unitVTU,'(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
  do i=1,nParticles
   pointBase = (i-1)*8
   write(unitVTU,'(8(1X,I0))') (pointBase + j - 1,j=1,8)
  end do
  write(unitVTU,'(A)') '        </DataArray>'

  write(unitVTU,'(A)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
  offset = 0
  do i=1,nParticles
   offset = offset + 8
   write(unitVTU,'(1X,I0)') offset
  end do
  write(unitVTU,'(A)') '        </DataArray>'

  write(unitVTU,'(A)') '        <DataArray type="UInt8" Name="types" format="ascii">'
  do i=1,nParticles
   write(unitVTU,'(1X,I0)') 12
  end do
  write(unitVTU,'(A)') '        </DataArray>'
  write(unitVTU,'(A)') '      </Cells>'

  write(unitVTU,'(A)') '      <CellData Scalars="pair_id">'
  write(unitVTU,'(A)') '        <DataArray type="Int32" Name="pair_id" format="ascii">'
  do i=1,nParticles
   write(unitVTU,'(1X,I0)') i
  end do
  write(unitVTU,'(A)') '        </DataArray>'
  write(unitVTU,'(A)') '        <DataArray type="Int32" Name="start_element" format="ascii">'
  do i=1,nParticles
   write(unitVTU,'(1X,I0)') particles(i)%indice
  end do
  write(unitVTU,'(A)') '        </DataArray>'
  write(unitVTU,'(A)') '      </CellData>'

  write(unitVTU,'(A)') '    </Piece>'
  write(unitVTU,'(A)') '  </UnstructuredGrid>'
  write(unitVTU,'(A)') '</VTKFile>'
  close(unitVTU)
end subroutine WriteParticleHex

end subroutine WriteRemainingParticlesVTU

subroutine PrepareDistanceToInflowOrdering()
 use def_FEAT, ONLY : NLMAX
 use var_QuadScalar, ONLY : mg_mesh
 use Sigma_User, only : myProcess
 use PP3D_MPI, ONLY : myid,master
 use xPart_def, ONLY : DistanceToInflow,InflowOrdering,nInflowDistanceBins,InflowBins
 implicit none
 include 'mpif.h'
 integer :: ILEV,nElem,elemIdx,iInflow,binIdx,position,localIdx,iErr
 real*8 :: elemCenter(3),center(3),diffVec(3)
 real*8 :: radialDist,distCandidate,minDist
 real*8 :: minDistance,maxDistance,binWidth,innerRadius,outerRadius,centerlineRadius
 integer, allocatable :: binCounts(:),binOffsets(:),binPosition(:)
 real*8, allocatable :: localHist(:),globalHist(:)

 ILEV = NLMAX
 nElem = mg_mesh%level(ILEV)%nel

 if (myProcess%nOfInflows.le.0) then
  if (myid.eq.1) write(*,*) "ERROR: No inflows defined for q2p1_xParticles."
  call MPI_ABORT(MPI_COMM_WORLD, 1, iErr)
 end if

 if (allocated(DistanceToInflow)) then
  if (size(DistanceToInflow).ne.nElem) deallocate(DistanceToInflow)
 end if
 if (.not.allocated(DistanceToInflow)) allocate(DistanceToInflow(nElem))

 if (allocated(InflowOrdering%ids)) then
  if (size(InflowOrdering%ids).ne.nElem) deallocate(InflowOrdering%ids)
 end if
 if (.not.allocated(InflowOrdering%ids)) allocate(InflowOrdering%ids(nElem))

 if (allocated(InflowOrdering%distances)) then
  if (size(InflowOrdering%distances).ne.nElem) deallocate(InflowOrdering%distances)
 end if
 if (.not.allocated(InflowOrdering%distances)) allocate(InflowOrdering%distances(nElem))

 if (allocated(InflowBins)) then
  if (size(InflowBins).ne.nInflowDistanceBins) then
   do binIdx=1,size(InflowBins)
    if (allocated(InflowBins(binIdx)%ids)) deallocate(InflowBins(binIdx)%ids)
    if (allocated(InflowBins(binIdx)%distances)) deallocate(InflowBins(binIdx)%distances)
   end do
   deallocate(InflowBins)
  end if
 end if
 if (.not.allocated(InflowBins)) allocate(InflowBins(nInflowDistanceBins))

 do binIdx=1,nInflowDistanceBins
  InflowBins(binIdx)%count = 0
  if (allocated(InflowBins(binIdx)%ids)) deallocate(InflowBins(binIdx)%ids)
  if (allocated(InflowBins(binIdx)%distances)) deallocate(InflowBins(binIdx)%distances)
 end do

 do elemIdx=1,nElem
  elemCenter = mg_mesh%level(ILEV)%dcorvg(:,mg_mesh%level(ILEV)%nvt + &
               mg_mesh%level(ILEV)%net + mg_mesh%level(ILEV)%nat + elemIdx)
  minDist = huge(1d0)
  do iInflow=1,myProcess%nOfInflows
   center = myProcess%myInflow(iInflow)%Center
   diffVec = elemCenter - center
   radialDist = sqrt(diffVec(1)**2 + diffVec(2)**2 + diffVec(3)**2)
   select case (myProcess%myInflow(iInflow)%iBCtype)
   case (2)
    innerRadius = myProcess%myInflow(iInflow)%innerradius
    outerRadius = myProcess%myInflow(iInflow)%outerradius
    centerlineRadius = 0.5d0*(innerRadius + outerRadius)
    distCandidate = abs(radialDist - centerlineRadius)
   case default
    distCandidate = radialDist
   end select
   if (distCandidate.lt.minDist) minDist = distCandidate
  end do
  DistanceToInflow(elemIdx) = minDist
 end do

 if (myid.eq.master) then
  write(*,*) "Computed inflow distances for all elements."
 end if

 minDistance = minval(DistanceToInflow)
 maxDistance = maxval(DistanceToInflow)
 if (maxDistance.le.minDistance) then
  binWidth = 0d0
 else
  binWidth = (maxDistance - minDistance)/real(nInflowDistanceBins,8)
 end if

 allocate(localHist(nInflowDistanceBins))
 allocate(globalHist(nInflowDistanceBins))
 localHist = 0d0
 globalHist = 0d0

 do elemIdx=1,nElem
  if (binWidth.eq.0d0) then
   binIdx = 1
  else
   binIdx = int((DistanceToInflow(elemIdx) - minDistance)/binWidth) + 1
  end if
  binIdx = max(1,min(nInflowDistanceBins,binIdx))
  localHist(binIdx) = localHist(binIdx) + 1d0
 end do

 call MPI_Allreduce(localHist,globalHist,nInflowDistanceBins,MPI_DOUBLE_PRECISION,&
                    MPI_SUM,MPI_COMM_WORLD,iErr)

 if (myid.eq.master) then
  write(*,*) "Global inflow-distance histogram collected."
 end if

 allocate(binCounts(nInflowDistanceBins))
 allocate(binOffsets(nInflowDistanceBins+1))
 allocate(binPosition(nInflowDistanceBins))
 binCounts = 0

 do elemIdx=1,nElem
  if (binWidth.eq.0d0) then
   binIdx = 1
  else
   binIdx = int((DistanceToInflow(elemIdx) - minDistance)/binWidth) + 1
  end if
  binIdx = max(1,min(nInflowDistanceBins,binIdx))
  binCounts(binIdx) = binCounts(binIdx) + 1
 end do

 binOffsets(1) = 1
 do binIdx=1,nInflowDistanceBins
  binOffsets(binIdx+1) = binOffsets(binIdx) + binCounts(binIdx)
  binPosition(binIdx) = binOffsets(binIdx)
  InflowBins(binIdx)%count = binCounts(binIdx)
  if (binCounts(binIdx).gt.0) then
   allocate(InflowBins(binIdx)%ids(binCounts(binIdx)))
   allocate(InflowBins(binIdx)%distances(binCounts(binIdx)))
  end if
 end do

 do elemIdx=1,nElem
  if (binWidth.eq.0d0) then
   binIdx = 1
  else
   binIdx = int((DistanceToInflow(elemIdx) - minDistance)/binWidth) + 1
  end if
  binIdx = max(1,min(nInflowDistanceBins,binIdx))
  position = binPosition(binIdx)
  binPosition(binIdx) = binPosition(binIdx) + 1
  InflowOrdering%ids(position) = elemIdx
  InflowOrdering%distances(position) = DistanceToInflow(elemIdx)
  if (InflowBins(binIdx)%count.gt.0) then
   localIdx = position - binOffsets(binIdx) + 1
  InflowBins(binIdx)%ids(localIdx) = elemIdx
  InflowBins(binIdx)%distances(localIdx) = DistanceToInflow(elemIdx)
 end if
end do

 if (myid.eq.master) then
  write(*,'(A,2ES12.4)') "Inflow distance range:", minDistance,maxDistance
  write(*,*) "Prepared inflow-aware element ordering (no chunk coupling yet)."
  write(*,*) "Inflow-distance ordering arrays populated."
 end if

 deallocate(localHist,globalHist,binCounts,binOffsets,binPosition)

end subroutine PrepareDistanceToInflowOrdering
