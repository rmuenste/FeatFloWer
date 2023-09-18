subroutine Transport_xParticles_MPI(log_unit)

use var_QuadScalar
USE PP3D_MPI, ONLY : myid,master,showid,subnodes,numnodes,MPI_COMM_SUBS
USE xPart_def, ONLY : xFactor,xChunks

implicit none

include 'mpif.h'

integer, intent(in) :: log_unit
integer nParticles,iChunk,nChunks,ivt_min,ivt_max,pID,pJD,i,mParticles,iaux
integer iStatus(MPI_STATUS_SIZE)
INTEGER iErr,indice,n,Stats(6)
real*8, allocatable :: daux(:,:),saux(:,:)
character :: cCSV_File*(128)
integer, allocatable :: sendcounts(:),displs(:)
real*8, allocatable :: gathered_data(:,:)
real*8 BoundingBox(3,2),P(3),dL(3)
logical bxOutput
integer ixOutput

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

! initialization of LostSet and CompleteSet and process the first data Chunk
if (myid.ne.master) then
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
 CALL Transport_xParticles(log_unit,iChunk,ivt_min,ivt_max,Stats)
 mParticles = myParticleParam%nParticles
end if

! master is distributing the data Chunks
if (myid.eq.0) THEN
 write(*,'(7A10)') " -- ","Chunk","pID","#ALL","#ALLinFD","Active","Lost"
 DO iChunk = subnodes + 1, nChunks
  CALL MPI_RECV(Stats,6,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,iStatus,iErr)
  pID = Stats(1)
  write(*,'(A10,6I10)') "done:",Stats(6),Stats(1:5)
  
  CALL MPI_SEND(iChunk,1,MPI_INT,pID,0,MPI_COMM_WORLD,iStatus,iErr)
 END DO
 
 DO pJD=1,subnodes
  iChunk = 0
  CALL MPI_RECV(Stats,6,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,iStatus,iErr)
  pID = Stats(1)
  write(*,'(A10,6I10)') "done:",Stats(6),Stats(1:5)
  
  CALL MPI_SEND(iChunk,1,MPI_INT,pID,0,MPI_COMM_WORLD,iStatus,iErr)
 END DO
 
END IF

! workers receive the new data Chunks and process them
if (myid.ne.0) then
 DO 
  CALL MPI_SEND(Stats,6,MPI_INT,0,0,MPI_COMM_WORLD,iStatus,iErr)
  CALL MPI_RECV(iChunk,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,iStatus,iErr)
  Stats(6) = iChunk
  if (iChunk.eq.0) GOTO 1
  
  ivt_min = (iChunk-1)*(nParticles/nChunks) + 1
  ivt_max = (iChunk+0)*(nParticles/nChunks) + 0
  if (iChunk.eq.nChunks) ivt_max = nParticles
  CALL Transport_xParticles(log_unit,iChunk,ivt_min,ivt_max,Stats)
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
  daux(:,indice) = [dble(myLostSet(indice)%ID),dble(myLostSet(indice)%indice),0d0]
 END DO
 call MPI_AllReduce(daux, saux, 3*nParticles, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
 DO indice=1,nParticles
  myLostSet(indice)%ID     = saux(1,indice)
  myLostSet(indice)%indice = saux(2,indice)
 END DO
 DeALLOCATE(daux,saux)
end if
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Lets collect all the particles being still inside of the FD  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 allocate(sendcounts(0:numnodes),displs(0:numnodes+1))
 sendcounts = 0; displs = 0
 
 call MPI_allgather(3*nCompleteSet, 1, MPI_INTEGER, sendcounts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

 displs = 0
 do i = 2, numnodes+1
   displs(i) = displs(i-1) + sendcounts(i-1)
 end do
 
 if (myid.eq.master) then
   allocate(gathered_data(3,displs(numnodes+1)/3))
   nCompleteSet = 0 
   n=0
 else 
  n = 3*nCompleteSet
  ALLOCATE(daux(3,nCompleteSet))
  DO i=1,nCompleteSet
   daux(:,i) = myCompleteSet(i)%coor
  END DO
 endif
 
 call MPI_Gatherv(daux, n, MPI_DOUBLE_PRECISION, &
                  gathered_data, sendcounts, displs, &
                  MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
                  
 
 call MPI_AllReduce(nCompleteSet, iaux, 1, MPI_INT, MPI_SUM, MPI_COMM_World, ierr)
 nCompleteSet = iaux
 
 if (myid.eq.master) then
  ALLOCATE(myCompleteSet(1:nCompleteSet))
  DO i=1,nCompleteSet
   myCompleteSet(i)%coor = gathered_data(:,i)
  END DO
  DeALLOCATE(gathered_data)
 else
  DeALLOCATE(daux)
 endif
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
   if (myLostSet(i)%indice.ne.0) then
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

implicit none
integer, intent(in) :: log_unit,iChunk,ivt_min,ivt_max
integer :: Stats(6)

real*8, allocatable, Dimension(:) :: xField,yField,zField
integer ndof,indice,i,iParticle,iTime
!!!!!!!!!!1
REAL*8 :: dTime,dStart
character :: cCSV_File*(128)
logical :: bCSV=.false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!! USER params !!!!!!!!!!!!!!!!!!!!!!!!
myParticleParam%OutflowZPos = 1d9
!!!!!!!!!!!!!!!!!!!!!!!!!!!! USER params !!!!!!!!!!!!!!!!!!!!!!!!

dTime     = 0d0

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
 
iParticle = 0
DO i=ivt_min,ivt_max
 indice = mg_mesh%level(ILEV)%nvt + mg_mesh%level(ILEV)%net + mg_mesh%level(ILEV)%nat +  i 
 if (Shell(indice).gt.minDist) then
  iParticle = iParticle + 1
  myActiveSet(iParticle)%coor(:) = mg_mesh%level(ILEV)%dcorvg(:,indice)
  myActiveSet(iParticle)%time    = 0d0
  myActiveSet(iParticle)%id      = 1
  myActiveSet(iParticle)%indice  = i
  myActiveSet(iParticle)%velo    = 0d0
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

DO iTime = 1,nTime
 dStart = dTime
 dTime  = dTime + dTimeStep
 
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
                          dTime,dStart)
 
!  write(*,*) "nActiveSet@iChunk",iChunk," = ", nActiveSet
 if (nActiveSet.eq.0) EXIT
 
if (bCSV) THEN
 write(cCSV_File,'(A,2(i4.4,A))') "P_",iChunk,"__",iTime,'.csv'
 OPEN(FILE=adjustl(trim(cCSV_File)),UNIT=11)
 WRITE(11,*) 'x,y,z,vx,vy,vz,t,id,indice'
 DO i=1,nActiveSet
  WRITE(11,'(7(ES12.4,(",")),1(I0,(",")),I0)') myActiveSet(i)%coor,myActiveSet(i)%velo,myActiveSet(i)%time,myActiveSet(i)%id,myActiveSet(i)%indice
 END DO
 CLOSE(11)
END IF
 
END DO                            

Stats(1) = myid                           !! signaure
Stats(2) = ivt_max-ivt_min + 1            !! all potential particles
Stats(3) = iParticle                      !! all particles in FD
Stats(4) = nActiveSet                     !! all particles remained in FD
Stats(5) = iParticle - nActiveSet         !! all particles which got lost

if (nActiveSet.ne.0) THEN
 DO i=1,nActiveSet
  nCompleteSet = nCompleteSet + 1
  myCompleteSet(nCompleteSet)%coor   = myActiveSet(i)%coor
  myCompleteSet(nCompleteSet)%time   = myActiveSet(i)%time
  myCompleteSet(nCompleteSet)%indice = myActiveSet(i)%indice
  myCompleteSet(nCompleteSet)%id     = myid
  myCompleteSet(nCompleteSet)%velo   = myActiveSet(i)%velo
 END DO
END IF

CALL FreeOctTree()  

deALLOCATE(myActiveSet)

end subroutine Transport_xParticles
