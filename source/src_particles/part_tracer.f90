MODULE Particle

USE types
use var_QuadScalar
USE Sigma_User, ONLY: mySigma,myProcess,myMultiMat

USE PP3D_MPI, ONLY : myid,master,showid,Comm_Summ,Comm_SummN,subnodes,COMM_Maximum,COMM_Minimum
! USE Sigma_User, ONLY: myRTD
USE UMFPackSolver, ONLY : myUmfPack_Factorize,myUmfPack_Solve

USE app_initialization, ONLY : init_sol_same_level,init_sol_repart

USE particles_input

use particle_step

use solution_io, only : write_sol_to_file, read_sol_from_file


INTEGER nBuffer
PARAMETER (nBuffer=40)

REAL*8, ALLOCATABLE :: Lambda(:),Length(:,:)
REAL*8, ALLOCATABLE :: ZPosParticles(:,:,:)
LOGICAL, ALLOCATABLE :: bPosParticles(:,:)

LOGICAL, ALLOCATABLE :: pPresent(:)
LOGICAL bOutputLambda
INTEGER iOutputLambda

!------------------------------------------------------------

CONTAINS

SUBROUTINE Transport_Particle_xse(mfile)
USE Transport_Q2P1, only : Create_MMat,mg_mlmat
INTEGER mfile
INTEGER i,nExSum,nActSum,nActSum0,nActSumOld,iCycle,iFile,iLevel0,iLevel1
REAL*8  dTime,daux,dPeriod,dTimeStep,dMeltVolume,dStart,dBuffer(nBuffer)
REAL*8  dMeltFlowRate,dCharacteristicTime
CHARACTER*99 cFile
LOGICAL :: bOutput=.TRUE.
real*8 xmin, xmax, ymin, ymax, zmin, zmax
integer iPeriodicityShift,jFile,iAngle,iMat

dTime=0d0
dPeriod = 6d1/myParticleParam%f
dTimeStep = dPeriod/DBLE(myParticleParam%nTimeLevels)
iAngle = 360/myParticleParam%nTimeLevels
iPeriodicityShift = INT(myParticleParam%nTimeLevels/myParticleParam%nPeriodicity)
ilev=nlmax
IF (.not.allocated(mySegmentIndicator))   allocate(mySegmentIndicator(2,mg_mesh%level(ilev)%nvt+&
                                          mg_mesh%level(ilev)%net+&
                                          mg_mesh%level(ilev)%nat+&
                                          mg_mesh%level(ilev)%nel))
! GOTO 222

myMatrixRenewal%D = 0
myMatrixRenewal%K = 0
myMatrixRenewal%S = 0
myMatrixRenewal%M = 0
myMatrixRenewal%C = 0

!!!!!!!!!!!!!!!!!!!  ---- Velocity Fields are to be loaded -----   !!!!!!!!!!!!!!!!!!!
 ALLOCATE(myVelo(0:myParticleParam%nTimeLevels-1))
DO iFile=0,myParticleParam%nTimeLevels/myParticleParam%nPeriodicity-1 !myParticleParam%nTimeLevels
 myParticleParam%dump_in_file = iFile*iAngle
 if (myParticleParam%DumpFormat.eq.1) THEN
  WRITE(cFile,'(I0)') myParticleParam%dump_in_file
  CALL init_sol_same_level(cFile)
 END IF
 
 if (myParticleParam%DumpFormat.eq.2) CALL Load_ListFiles_PRT_Tracer(myParticleParam%dump_in_file)
 
 if (myParticleParam%DumpFormat.eq.4) CALL LoadMPIDumpFiles(iFile*iAngle,'v,x,s')
 
 ALLOCATE(myVelo(iFile)%x(QuadSc%ndof))
 ALLOCATE(myVelo(iFile)%y(QuadSc%ndof))
 ALLOCATE(myVelo(iFile)%z(QuadSc%ndof))
 IF (ADJUSTL(TRIM(mySigma%cType)).EQ."DIE") THEN
  if (myParticleParam%bBacktrace) THEN
   myVelo(iFile)%x = -QuadSc%ValU
   myVelo(iFile)%y = -QuadSc%ValV
   myVelo(iFile)%z = -QuadSc%ValW
  else
   myVelo(iFile)%x = QuadSc%ValU
   myVelo(iFile)%y = QuadSc%ValV
   myVelo(iFile)%z = QuadSc%ValW
  end if
 ELSE
  myVelo(iFile)%x = QuadSc%ValU
  myVelo(iFile)%y = QuadSc%ValV
  myVelo(iFile)%z = QuadSc%ValW
 END IF
 IF (myid.eq.1) WRITE(*,*) 'File ',iFile, ' is loaded... ==> angle :', iFile*iAngle

 if (myParticleParam%nPeriodicity.gt.1) then
  DO iPerio = 1,myParticleParam%nPeriodicity-1
   jFile = iFile + iPeriodicityShift*(iPerio) 
   ALLOCATE(myVelo(jFile)%x(QuadSc%ndof))
   ALLOCATE(myVelo(jFile)%y(QuadSc%ndof))
   ALLOCATE(myVelo(jFile)%z(QuadSc%ndof))
 IF (ADJUSTL(TRIM(mySigma%cType)).EQ."DIE") THEN
  if (myParticleParam%bBacktrace) THEN
   myVelo(jFile)%x = -QuadSc%ValU
   myVelo(jFile)%y = -QuadSc%ValV
   myVelo(jFile)%z = -QuadSc%ValW
  else
   myVelo(jFile)%x = QuadSc%ValU
   myVelo(jFile)%y = QuadSc%ValV
   myVelo(jFile)%z = QuadSc%ValW
  end if
 ELSE
  myVelo(jFile)%x = QuadSc%ValU
  myVelo(jFile)%y = QuadSc%ValV
  myVelo(jFile)%z = QuadSc%ValW
 END IF
   IF (myid.eq.1) WRITE(*,*) 'File ',jFile, ' is loaded... ==> angle :', jFile*iAngle
  END DO
 end if
END DO
! !!!!!!!!!!!!!!!!!!!  ---- Velocity Fields are to be loaded -----   !!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!! Volume estimation !!!!!!!!!!!!!!!!!!!!!!!!!!!! 
IF (ADJUSTL(TRIM(mySigma%cType)).EQ."DIE") THEN

 CALL Create_MMat()
 IF (myid.eq.showID) WRITE(MTERM,*) 
 IF (myid.eq.showID) WRITE(MFILE,*) 
 dMeltVolume = 0d0
 DO i=1,QuadSc%ndof
  IF (mySegmentIndicator(2,i).gt.0.5d0) then
   dMeltVolume = dMeltVolume + mg_mlmat(nlmax)%a(i)
  END IF
 END DO
 CALL Comm_Summ(dMeltVolume)

 dMeltFlowRate = 0d0
 DO iInflow=1,myProcess%nOfInflows
  iMat          = myProcess%myInflow(iInflow)%Material
  dMeltFlowRate = dMeltFlowRate + myProcess%myInflow(iInflow)%massflowrate * myMultiMat%Mat(iMat)%Thermodyn%density ! [kg/h]/[kg/l]=[l/h]
 END DO

 dCharacteristicTime = dMeltVolume*3.6d0/dMeltFlowRate ! [cm3]/[l/h] = 0.001*[l]/[l/3600s] = 3.6*[s]  

 IF (myid.eq.showID) WRITE(MTERM,'(A,3ES12.4)') "MeltVolume[l],MeltFlowrate[l/h],RT[s]:",1e-3*dMeltVolume,dMeltFlowRate,dCharacteristicTime
 IF (myid.eq.showID) WRITE(MFILE,'(A,3ES12.4)') "MeltVolume[l],MeltFlowrate[l/h],RT[s]:",1e-3*dMeltVolume,dMeltFlowRate,dCharacteristicTime

 dTimeStep = dCharacteristicTime/1d2
 
END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!! Volume estimation !!!!!!!!!!!!!!!!!!!!!!!!!!!! 

IF (myid.ne.0) THEN
 ILEV=NLMAX-1
 CALL SETLEV(2)
END IF
nStartActiveSet = 0
CALL Extract_Particle_xse(mg_mesh%level(ILEV)%dcorvg,&
                      mg_mesh%level(ILEV)%kvert,&
                      mg_mesh%level(ILEV)%kedge,&
                      mg_mesh%level(ILEV)%karea,&
                      mg_mesh%level(ILEV)%elementsAtVertexIdx,&
                      mg_mesh%level(ILEV)%elementsAtVertex,&
                      QuadSc%ValU,QuadSc%ValV,QuadSc%ValW,&
                      nvt,net,nat,nel,dTime,myParticleParam%d_CorrDist)

! Get the bounds of the mesh. Note: We need to set it to nlmax because else
! it won't work
IF (myid.ne.0) THEN
 ILEV=NLMAX
 CALL SETLEV(2)
END IF

call GetMeshBounds(mg_mesh%level(NLMAX)%dcorvg,nvt,xmin,xmax,ymin,ymax,zmin,zmax)
myMeshInfo%xmin = xmin
myMeshInfo%xmax = xmax
myMeshInfo%ymin = ymin
myMeshInfo%ymax = ymax
myMeshInfo%zmin = zmin
myMeshInfo%zmax = zmax
if (myid .eq. showid) then
  write(mfile,'(A,E12.4,A7,E12.4)') "MESH-Bounds in X-Dimension xmin: ",myMeshInfo%xmin, " xmax: ", myMeshInfo%xmax
  write(mfile,'(A,E12.4,A7,E12.4)') "MESH-Bounds in Y-Dimension ymin: ",myMeshInfo%ymin, " ymax: ", myMeshInfo%ymax
  write(mfile,'(A,E12.4,A7,E12.4)') "MESH-Bounds in Z-Dimension zmin: ",myMeshInfo%zmin, " zmax: ", myMeshInfo%zmax
  write(mterm,'(A,E12.4,A7,E12.4)') "MESH-Bounds in X-Dimension xmin: ",myMeshInfo%xmin, " xmax: ", myMeshInfo%xmax
  write(mterm,'(A,E12.4,A7,E12.4)') "MESH-Bounds in Y-Dimension ymin: ",myMeshInfo%ymin, " ymax: ", myMeshInfo%ymax
  write(mterm,'(A,E12.4,A7,E12.4)') "MESH-Bounds in Z-Dimension zmin: ",myMeshInfo%zmin, " zmax: ", myMeshInfo%zmax
end if ! Dump mesh-info

IF (myid.eq.1) OPEN(FILE='_RTD/RTD.csv',UNIT=947)
dBuffer = 0d0
nLostSet = 0

IF (myid.eq.1) WRITE(947,'(3(E12.4,A))') 0d0,', ',0d0,', ',0d0,' '

! Output initial positions of particles to csv.
IF (bOutput) THEN
 CALL OutputParticlesCSV(0)
END IF

DO iTimeSteps=1,myParticleParam%nRotation*myParticleParam%nTimeLevels

  nStartActiveSet = 0
  iCycle = 0
  iLevel0 = MOD(iTimeSteps-1,myParticleParam%nTimeLevels)
  iLevel1 = MOD(iTimeSteps  ,myParticleParam%nTimeLevels)

  dStart = dTime
  dTime = dTime + dTimeStep

  IF (myid.eq.1) WRITE(MFILE,'(A,I8,A,I8,A,I0,A,I0,ES14.6)') 'Timestep ',iTimeSteps, ' / ',myParticleParam%nRotation*myParticleParam%nTimeLevels,'  to be performed ...',iLevel0,'/',iLevel1,dTime
  IF (myid.eq.1) WRITE(MTERM,'(A,I8,A,I8,A,I0,A,I0,ES14.6)') 'Timestep ',iTimeSteps, ' / ',myParticleParam%nRotation*myParticleParam%nTimeLevels,'  to be performed ...',iLevel0,'/',iLevel1,dTime

  DO
   IF (myid.ne.0) THEN
    CALL Move_Particle_xse(mg_mesh%level(ILEV)%dcorvg,&
                       mg_mesh%level(ILEV)%kvert,&
                       mg_mesh%level(ILEV)%kedge,&
                       mg_mesh%level(ILEV)%karea,&
                       mg_mesh%level(ILEV)%elementsAtVertexIdx,&
                       mg_mesh%level(ILEV)%elementsAtVertex,&
                       myVelo(iLevel0)%x,myVelo(iLevel0)%y,myVelo(iLevel0)%z,&
                       myVelo(iLevel1)%x,myVelo(iLevel1)%y,myVelo(iLevel1)%z,&
                       nvt,net,nat,nel,dTime,dStart,myParticleParam%d_CorrDist)
   END IF

   daux = DBLE(nActiveSet)
   CALL Comm_Summ(daux)
   nActSum =INT(daux)
   IF (iTimeSteps.EQ.1) nActSumOld = nActSum
   IF (iTimeSteps.EQ.1) nActSum0   = nActSum

   daux = DBLE(nExchangeSet)
   CALL Comm_Summ(daux)
   nExSum =INT(daux)

   iCycle = iCycle + 1
   IF (myid.eq.1) WRITE(MFILE,'(A,I0,A,I0,A,I0)') 'Exchange_of_particles: Cycle: ',iCycle, ' Num_Of_Points: ',nExSum,'/',nActSum
   IF (myid.eq.1) WRITE(MTERM,'(A,I0,A,I0,A,I0)') 'Exchange_of_particles: Cycle: ',iCycle, ' Num_Of_Points: ',nExSum,'/',nActSum


   IF (nExSum.EQ.0) THEN
    IF (myid.eq.1) WRITE(MFILE,'(A)') '...........................................................................'
    IF (myid.eq.1) WRITE(MTERM,'(A)') '...........................................................................'
    EXIT
   ELSE

    CALL Exchange_Particle(nExSum)
    CALL Extract_Particle_xse(mg_mesh%level(ILEV)%dcorvg,&
                          mg_mesh%level(ILEV)%kvert,&
                          mg_mesh%level(ILEV)%kedge,&
                          mg_mesh%level(ILEV)%karea,&
                          mg_mesh%level(ILEV)%elementsAtVertexIdx,&
                          mg_mesh%level(ILEV)%elementsAtVertex,&
                          QuadSc%ValU,QuadSc%ValV,QuadSc%ValW,&
                          nvt,net,nat,nel,dTime,myParticleParam%d_CorrDist)

   END IF

  END DO

  IF (bOutput) THEN

   CALL OutputParticlesCSV(iTimeSteps)

   bOutputLambda = .FALSE.
   IF (MOD(iTimeSteps-1,myParticleParam%nTimeLevels).EQ.0) bOutputLambda=.TRUE.
   iOutputLambda = iTimeSteps/myParticleParam%nTimeLevels

   ! not the end of the simulation - so call getLambda with false
   CALL GetLambda(.false.)

  END IF

  dStat = 0d0
  DO i=1,nBuffer
   dStat = dStat +  dBuffer(i)
  END DO
  dStat = dStat/REAL(nBuffer)

  CALL GetCutplanes()

  IF (myid.eq.1) WRITE(947,'(3(E12.4,A))') (REAL(iTimeSteps)/REAL(myParticleParam%nTimeLevels))/(myParticleParam%f/6d1),', ',REAL(nActSum0-nActSum)/REAL(nActSum0),', ',dStat/REAL(nActSum0),' '
  !IF (myid.eq.1) WRITE(947,'(3(E12.4,A))') REAL(iTimeSteps),', ',REAL(nActSum0-nActSum)/REAL(nActSum0),', ',dStat/REAL(nActSum0),' '
  iBuffer = MOD(iTimeSteps,nBuffer)+1
  dBuffer(iBuffer) = (nActSumOld-nActSum)
  nActSumOld = nActSum

  IF (DBLE(nActSum0-nActSum)/DBLE(nActSum0).ge.myParticleParam%minFrac) EXIT
  

END DO

IF (bOutput) THEN

 CALL OutputLostParticlesCSV()

END IF

CALL OutputParticlesAtZtoCSV()

! Now its the end of the simulation - so output the final lambda
call GetLambda(.true.)

IF (myid.eq.1) CLOSE(947)

! 222 CONTINUE
!
! CALL PostProcessRTD()

END SUBROUTINE Transport_Particle_xse
!
! --------------------------------------------------------------------
!
SUBROUTINE Transport_Particle_Arch(mfile)
INTEGER mfile
INTEGER i,nExSum,nActSum,nActSum0,nActSumOld,iCycle,iFile,iLevel0,iLevel1
REAL*8  dTime,daux,dPeriod,dTimeStep,dStart,dBuffer(nBuffer)
CHARACTER*99 cFile
LOGICAL :: bOutput=.TRUE.
real*8 xmin, xmax, ymin, ymax, zmin, zmax
integer :: iOutputIndex=0
real*8, allocatable :: FlowField_U(:),FlowField_V(:),FlowField_W(:)
real*8 :: dPeriodicTime = 0.245d0,dPeriodicTimeUnit,dVeloTime_0,dVeloTime_1,dVeloTime,dVeloTime_F

dTime=0d0
dPeriod = 6d1/myParticleParam%f
dTimeStep = dPeriod/DBLE(myParticleParam%nTimeLevels)

! GOTO 222

!!!!!!!!!!!!!!!!!!  ---- Velocity Fields are to be loaded -----   !!!!!!!!!!!!!!!!!!!
ALLOCATE(myVelo(myParticleParam%nTimeLevels))
ALLOCATE(FlowField_U(QuadSc%ndof))
ALLOCATE(FlowField_V(QuadSc%ndof))
ALLOCATE(FlowField_W(QuadSc%ndof))
dPeriodicTimeUnit = dPeriodicTime / dble(myParticleParam%nTimeLevels)

DO iFile=1,myParticleParam%nTimeLevels
 WRITE(cFile,'(I0)')iFile
 CALL init_sol_same_level(CFile)
 ALLOCATE(myVelo(iFile)%x(QuadSc%ndof))
 ALLOCATE(myVelo(iFile)%y(QuadSc%ndof))
 ALLOCATE(myVelo(iFile)%z(QuadSc%ndof))
 myVelo(iFile)%x = QuadSc%ValU
 myVelo(iFile)%y = QuadSc%ValV
 myVelo(iFile)%z = QuadSc%ValW
END DO

!!!!!!!!!!!!!!!!!!!  ---- Velocity Fields are to be loaded -----   !!!!!!!!!!!!!!!!!!!

myMatrixRenewal%D = 0
myMatrixRenewal%K = 0
myMatrixRenewal%S = 0
myMatrixRenewal%M = 0
myMatrixRenewal%C = 0

! ALLOCATE(myVelo(1))
! 
! !!!!!!!!!!!!!!!!!!! choose one fof these !!!!!!!!!!!!!!!!!!!!!!
! if (myParticleParam%DumpFormat.eq.1) THEN
!  WRITE(cFile,'(I0)') myParticleParam%dump_in_file
!  CALL init_sol_same_level(cFile)
! END IF
! if (myParticleParam%DumpFormat.eq.2) CALL Load_ListFiles_PRT_Tracer(myParticleParam%dump_in_file)
! if (myParticleParam%DumpFormat.eq.3) THEN
!  WRITE(cFile,'(I0)') myParticleParam%dump_in_file
!  CALL init_sol_repart(cFile)
! END IF
! 
! ALLOCATE(myVelo(1)%x(QuadSc%ndof))
! ALLOCATE(myVelo(1)%y(QuadSc%ndof))
! ALLOCATE(myVelo(1)%z(QuadSc%ndof))
! myVelo(1)%x = QuadSc%ValU
! myVelo(1)%y = QuadSc%ValV
! myVelo(1)%z = QuadSc%ValW

IF (myid.ne.0) THEN
 ILEV=NLMAX
!  ILEV=NLMAX-1
 CALL SETLEV(2)
END IF
nStartActiveSet = 0
CALL Extract_Particle(mg_mesh%level(ILEV)%dcorvg,&
                      mg_mesh%level(ILEV)%kvert,&
                      mg_mesh%level(ILEV)%kedge,&
                      mg_mesh%level(ILEV)%karea,&
                      mg_mesh%level(ILEV)%elementsAtVertexIdx,&
                      mg_mesh%level(ILEV)%elementsAtVertex,&
                      QuadSc%ValU,QuadSc%ValV,QuadSc%ValW,&
                      nvt,net,nat,nel,dTime,myParticleParam%d_CorrDist)

! Get the bounds of the mesh. Note: We need to set it to nlmax because else
! it won't work
IF (myid.ne.0) THEN
 ILEV=NLMAX
 CALL SETLEV(2)
END IF

call GetMeshBounds(mg_mesh%level(NLMAX)%dcorvg,nvt,xmin,xmax,ymin,ymax,zmin,zmax)
myMeshInfo%xmin = xmin
myMeshInfo%xmax = xmax
myMeshInfo%ymin = ymin
myMeshInfo%ymax = ymax
myMeshInfo%zmin = zmin
myMeshInfo%zmax = zmax
if (myid .eq. showid) then
  write(mfile,'(A,E12.4,A7,E12.4)') "MESH-Bounds in X-Dimension xmin: ",myMeshInfo%xmin, " xmax: ", myMeshInfo%xmax
  write(mfile,'(A,E12.4,A7,E12.4)') "MESH-Bounds in Y-Dimension ymin: ",myMeshInfo%ymin, " ymax: ", myMeshInfo%ymax
  write(mfile,'(A,E12.4,A7,E12.4)') "MESH-Bounds in Z-Dimension zmin: ",myMeshInfo%zmin, " zmax: ", myMeshInfo%zmax
  write(mterm,'(A,E12.4,A7,E12.4)') "MESH-Bounds in X-Dimension xmin: ",myMeshInfo%xmin, " xmax: ", myMeshInfo%xmax
  write(mterm,'(A,E12.4,A7,E12.4)') "MESH-Bounds in Y-Dimension ymin: ",myMeshInfo%ymin, " ymax: ", myMeshInfo%ymax
  write(mterm,'(A,E12.4,A7,E12.4)') "MESH-Bounds in Z-Dimension zmin: ",myMeshInfo%zmin, " zmax: ", myMeshInfo%zmax
end if ! Dump mesh-info

IF (myid.eq.1) OPEN(FILE='_RTD/RTD.csv',UNIT=947)
dBuffer = 0d0
nLostSet = 0

IF (myid.eq.1) WRITE(947,'(3(E12.4,A))') 0d0,', ',0d0,', ',0d0,' '

! Output initial positions of particles to csv.
IF (bOutput) THEN
 CALL OutputParticlesCSV(0)
END IF

DO iTimeSteps=1,myParticleParam%nRotation*myParticleParam%nTimeLevels

  nStartActiveSet = 0
  iCycle = 0
  iLevel0 = MOD(iTimeSteps-1,myParticleParam%nTimeLevels)+1
  iLevel1 = MOD(iTimeSteps  ,myParticleParam%nTimeLevels)+1
  
  dVeloTime          = mod(dStart+0.5d0*dTimeStep,dPeriodicTime)
  dVeloTime_0        = (dVeloTime - mod(dVeloTime,dPeriodicTimeUnit)) / dPeriodicTimeUnit
  dVeloTime_1        = dVeloTime_0 + 1.0
  dVeloTime_F        = mod(dVeloTime,dPeriodicTimeUnit) / dPeriodicTimeUnit
  if ((dVeloTime_1+1d0).gt.dble(myParticleParam%nTimeLevels)) dVeloTime_1 = 0d0
  FlowField_U = myVelo(INT(dVeloTime_1)+1)%x*(dVeloTime_F) + myVelo(INT(dVeloTime_0)+1)%x*(1d0-dVeloTime_F)
  FlowField_V = myVelo(INT(dVeloTime_1)+1)%y*(dVeloTime_F) + myVelo(INT(dVeloTime_0)+1)%y*(1d0-dVeloTime_F)
  FlowField_W = myVelo(INT(dVeloTime_1)+1)%z*(dVeloTime_F) + myVelo(INT(dVeloTime_0)+1)%z*(1d0-dVeloTime_F)
  
  if (myid.eq.1) write(*,*) dVeloTime_0,dVeloTime_1,dVeloTime_F
    
  dStart = dTime
  dTime = dTime + dTimeStep

  IF (myid.eq.1) WRITE(MFILE,'(A,I8,A,I8,A,I0,A,I0,ES14.6)') 'Timestep ',iTimeSteps, ' / ',myParticleParam%nRotation*myParticleParam%nTimeLevels,'  to be performed ...',iLevel0,'/',iLevel1,dTime
  IF (myid.eq.1) WRITE(MTERM,'(A,I8,A,I8,A,I0,A,I0,ES14.6)') 'Timestep ',iTimeSteps, ' / ',myParticleParam%nRotation*myParticleParam%nTimeLevels,'  to be performed ...',iLevel0,'/',iLevel1,dTime

  DO
   IF (myid.ne.0) THEN
    CALL Move_Particle(mg_mesh%level(ILEV)%dcorvg,&
                       mg_mesh%level(ILEV)%kvert,&
                       mg_mesh%level(ILEV)%kedge,&
                       mg_mesh%level(ILEV)%karea,&
                       mg_mesh%level(ILEV)%elementsAtVertexIdx,&
                       mg_mesh%level(ILEV)%elementsAtVertex,&
                       FlowField_U,FlowField_V,FlowField_W,&
                       FlowField_U,FlowField_V,FlowField_W,&
                       nvt,net,nat,nel,dTime,dStart,myParticleParam%d_CorrDist)
   END IF

   daux = DBLE(nActiveSet)
   CALL Comm_Summ(daux)
   nActSum =INT(daux)
   IF (iTimeSteps.EQ.1) nActSumOld = nActSum
   IF (iTimeSteps.EQ.1) nActSum0   = nActSum

   daux = DBLE(nExchangeSet)
   CALL Comm_Summ(daux)
   nExSum =INT(daux)

   iCycle = iCycle + 1
   IF (myid.eq.1) WRITE(MFILE,'(A,I0,A,I0,A,I0)') 'Exchange_of_particles: Cycle: ',iCycle, ' Num_Of_Points: ',nExSum,'/',nActSum
   IF (myid.eq.1) WRITE(MTERM,'(A,I0,A,I0,A,I0)') 'Exchange_of_particles: Cycle: ',iCycle, ' Num_Of_Points: ',nExSum,'/',nActSum


   IF (nExSum.EQ.0) THEN
    IF (myid.eq.1) WRITE(MFILE,'(A)') '...........................................................................'
    IF (myid.eq.1) WRITE(MTERM,'(A)') '...........................................................................'
    EXIT
   ELSE

    CALL Exchange_Particle(nExSum)
    CALL Extract_Particle(mg_mesh%level(ILEV)%dcorvg,&
                          mg_mesh%level(ILEV)%kvert,&
                          mg_mesh%level(ILEV)%kedge,&
                          mg_mesh%level(ILEV)%karea,&
                          mg_mesh%level(ILEV)%elementsAtVertexIdx,&
                          mg_mesh%level(ILEV)%elementsAtVertex,&
                          QuadSc%ValU,QuadSc%ValV,QuadSc%ValW,&
                          nvt,net,nat,nel,dTime,myParticleParam%d_CorrDist)

   END IF

  END DO

  !!! here we need to collide the particles against each other
  CALL CheckForCollision(nvt,net,nat,nel,myParticleParam%d_CorrDist)
  
  IF (bOutput.and.mod(iTimeSteps,myParticleParam%OutputFreq).eq.0) THEN
  
   iOutputIndex = iOutputIndex + 1

   CALL OutputParticlesCSV(iOutputIndex)

   bOutputLambda = .FALSE.
   IF (MOD(iTimeSteps-1,myParticleParam%nTimeLevels).EQ.0) bOutputLambda=.TRUE.
   iOutputLambda = iTimeSteps/myParticleParam%nTimeLevels

   ! not the end of the simulation - so call getLambda with false
   CALL GetLambda(.false.)

  END IF

  dStat = 0d0
  DO i=1,nBuffer
   dStat = dStat +  dBuffer(i)
  END DO
  dStat = dStat/REAL(nBuffer)

  CALL GetCutplanes()

  IF (myid.eq.1) WRITE(947,'(3(E12.4,A))') (REAL(iTimeSteps)/REAL(myParticleParam%nTimeLevels))/(myParticleParam%f/6d1),', ',REAL(nActSum0-nActSum)/REAL(nActSum0),', ',dStat/REAL(nActSum0),' '
  !IF (myid.eq.1) WRITE(947,'(3(E12.4,A))') REAL(iTimeSteps),', ',REAL(nActSum0-nActSum)/REAL(nActSum0),', ',dStat/REAL(nActSum0),' '
  iBuffer = MOD(iTimeSteps,nBuffer)+1
  dBuffer(iBuffer) = (nActSumOld-nActSum)
  nActSumOld = nActSum

  IF (DBLE(nActSum0-nActSum)/DBLE(nActSum0).ge.myParticleParam%minFrac) EXIT
  

END DO

IF (bOutput) THEN

 CALL OutputLostParticlesCSV()

END IF

CALL OutputParticlesAtZtoCSV()

! Now its the end of the simulation - so output the final lambda
call GetLambda(.true.)

IF (myid.eq.1) CLOSE(947)

! 222 CONTINUE
!
! CALL PostProcessRTD()

END SUBROUTINE Transport_Particle_Arch
!
! --------------------------------------------------------------------
!
SUBROUTINE Transport_Particle(mfile)
INTEGER mfile
INTEGER i,nExSum,nActSum,nActSum0,nActSumOld,iCycle,iFile,iLevel0,iLevel1
REAL*8  dTime,daux,dPeriod,dTimeStep,dStart,dBuffer(nBuffer)
CHARACTER*99 cFile
LOGICAL :: bOutput=.TRUE.
real*8 xmin, xmax, ymin, ymax, zmin, zmax
integer :: iOutputIndex=0,iOutputFreq=1

dTime=0d0
dPeriod = 6d1/myParticleParam%f
dTimeStep = dPeriod/DBLE(myParticleParam%nTimeLevels)

! GOTO 222

!!!!!!!!!!!!!!!!!!!  ---- Velocity Fields are to be loaded -----   !!!!!!!!!!!!!!!!!!!
! ALLOCATE(myVelo(myParticleParam%nTimeLevels))
! WRITE(cFile,'(A,I2.2)') '_dump/',1
! CALL SolFromFile(CFile,1)
!
! DO iFile=1,myParticleParam%nTimeLevels
! ! CALL LoadSmartDumpFiles(cFile,1)
!  ALLOCATE(myVelo(iFile)%x(QuadSc%ndof))
!  ALLOCATE(myVelo(iFile)%y(QuadSc%ndof))
!  ALLOCATE(myVelo(iFile)%z(QuadSc%ndof))
!  myVelo(iFile)%x = QuadSc%ValU
!  myVelo(iFile)%y = QuadSc%ValV
!  myVelo(iFile)%z = QuadSc%ValW
! ! IF (myid.eq.1) WRITE(*,*) 'File ',iFile, ' is loaded...'
! END DO
! !!!!!!!!!!!!!!!!!!!  ---- Velocity Fields are to be loaded -----   !!!!!!!!!!!!!!!!!!!

myMatrixRenewal%D = 0
myMatrixRenewal%K = 0
myMatrixRenewal%S = 0
myMatrixRenewal%M = 0
myMatrixRenewal%C = 0

ALLOCATE(myVelo(1))

! CALL SolFromFile(cFile,1)

!!!!!!!!!!!!!!!!!!! choose one fof these !!!!!!!!!!!!!!!!!!!!!!
if (myParticleParam%DumpFormat.eq.1) THEN
 WRITE(cFile,'(I0)') myParticleParam%dump_in_file
 CALL init_sol_same_level(cFile)
END IF
if (myParticleParam%DumpFormat.eq.2) CALL Load_ListFiles_PRT_Tracer(myParticleParam%dump_in_file)
if (myParticleParam%DumpFormat.eq.3) THEN
 WRITE(cFile,'(I0)') myParticleParam%dump_in_file
 CALL init_sol_repart(cFile)
END IF
if (myParticleParam%DumpFormat.eq.4) CALL LoadMPIDumpFiles(iFile*iAngle,'v,x')

ALLOCATE(myVelo(1)%x(QuadSc%ndof))
ALLOCATE(myVelo(1)%y(QuadSc%ndof))
ALLOCATE(myVelo(1)%z(QuadSc%ndof))
myVelo(1)%x = QuadSc%ValU
myVelo(1)%y = QuadSc%ValV
myVelo(1)%z = QuadSc%ValW

IF (myid.ne.0) THEN
 ILEV=NLMAX
!  ILEV=NLMAX-1
 CALL SETLEV(2)
END IF
nStartActiveSet = 0
CALL Extract_Particle(mg_mesh%level(ILEV)%dcorvg,&
                      mg_mesh%level(ILEV)%kvert,&
                      mg_mesh%level(ILEV)%kedge,&
                      mg_mesh%level(ILEV)%karea,&
                      mg_mesh%level(ILEV)%elementsAtVertexIdx,&
                      mg_mesh%level(ILEV)%elementsAtVertex,&
                      QuadSc%ValU,QuadSc%ValV,QuadSc%ValW,&
                      nvt,net,nat,nel,dTime,myParticleParam%d_CorrDist)

! Get the bounds of the mesh. Note: We need to set it to nlmax because else
! it won't work
IF (myid.ne.0) THEN
 ILEV=NLMAX
 CALL SETLEV(2)
END IF

call GetMeshBounds(mg_mesh%level(NLMAX)%dcorvg,nvt,xmin,xmax,ymin,ymax,zmin,zmax)
myMeshInfo%xmin = xmin
myMeshInfo%xmax = xmax
myMeshInfo%ymin = ymin
myMeshInfo%ymax = ymax
myMeshInfo%zmin = zmin
myMeshInfo%zmax = zmax
if (myid .eq. showid) then
  write(mfile,'(A,E12.4,A7,E12.4)') "MESH-Bounds in X-Dimension xmin: ",myMeshInfo%xmin, " xmax: ", myMeshInfo%xmax
  write(mfile,'(A,E12.4,A7,E12.4)') "MESH-Bounds in Y-Dimension ymin: ",myMeshInfo%ymin, " ymax: ", myMeshInfo%ymax
  write(mfile,'(A,E12.4,A7,E12.4)') "MESH-Bounds in Z-Dimension zmin: ",myMeshInfo%zmin, " zmax: ", myMeshInfo%zmax
  write(mterm,'(A,E12.4,A7,E12.4)') "MESH-Bounds in X-Dimension xmin: ",myMeshInfo%xmin, " xmax: ", myMeshInfo%xmax
  write(mterm,'(A,E12.4,A7,E12.4)') "MESH-Bounds in Y-Dimension ymin: ",myMeshInfo%ymin, " ymax: ", myMeshInfo%ymax
  write(mterm,'(A,E12.4,A7,E12.4)') "MESH-Bounds in Z-Dimension zmin: ",myMeshInfo%zmin, " zmax: ", myMeshInfo%zmax
end if ! Dump mesh-info

IF (myid.eq.1) OPEN(FILE='_RTD/RTD.csv',UNIT=947)
dBuffer = 0d0
nLostSet = 0

IF (myid.eq.1) WRITE(947,'(3(E12.4,A))') 0d0,', ',0d0,', ',0d0,' '

! Output initial positions of particles to csv.
IF (bOutput) THEN
 CALL OutputParticlesCSV(0)
END IF

DO iTimeSteps=1,myParticleParam%nRotation*myParticleParam%nTimeLevels

  nStartActiveSet = 0
  iCycle = 0
  iLevel0 = MOD(iTimeSteps-1,myParticleParam%nTimeLevels)+1
  iLevel1 = MOD(iTimeSteps  ,myParticleParam%nTimeLevels)+1

  dStart = dTime
  dTime = dTime + dTimeStep

  IF (myid.eq.1) WRITE(MFILE,'(A,I8,A,I8,A,I0,A,I0,ES14.6)') 'Timestep ',iTimeSteps, ' / ',myParticleParam%nRotation*myParticleParam%nTimeLevels,'  to be performed ...',iLevel0,'/',iLevel1,dTime
  IF (myid.eq.1) WRITE(MTERM,'(A,I8,A,I8,A,I0,A,I0,ES14.6)') 'Timestep ',iTimeSteps, ' / ',myParticleParam%nRotation*myParticleParam%nTimeLevels,'  to be performed ...',iLevel0,'/',iLevel1,dTime

  DO
   IF (myid.ne.0) THEN
    CALL Move_Particle(mg_mesh%level(ILEV)%dcorvg,&
                       mg_mesh%level(ILEV)%kvert,&
                       mg_mesh%level(ILEV)%kedge,&
                       mg_mesh%level(ILEV)%karea,&
                       mg_mesh%level(ILEV)%elementsAtVertexIdx,&
                       mg_mesh%level(ILEV)%elementsAtVertex,&
                       myVelo(1)%x,myVelo(1)%y,myVelo(1)%z,&
                       myVelo(1)%x,myVelo(1)%y,myVelo(1)%z,&
                       nvt,net,nat,nel,dTime,dStart,myParticleParam%d_CorrDist)
   END IF

   daux = DBLE(nActiveSet)
   CALL Comm_Summ(daux)
   nActSum =INT(daux)
   IF (iTimeSteps.EQ.1) nActSumOld = nActSum
   IF (iTimeSteps.EQ.1) nActSum0   = nActSum

   daux = DBLE(nExchangeSet)
   CALL Comm_Summ(daux)
   nExSum =INT(daux)

   iCycle = iCycle + 1
   IF (myid.eq.1) WRITE(MFILE,'(A,I0,A,I0,A,I0)') 'Exchange_of_particles: Cycle: ',iCycle, ' Num_Of_Points: ',nExSum,'/',nActSum
   IF (myid.eq.1) WRITE(MTERM,'(A,I0,A,I0,A,I0)') 'Exchange_of_particles: Cycle: ',iCycle, ' Num_Of_Points: ',nExSum,'/',nActSum


   IF (nExSum.EQ.0) THEN
    IF (myid.eq.1) WRITE(MFILE,'(A)') '...........................................................................'
    IF (myid.eq.1) WRITE(MTERM,'(A)') '...........................................................................'
    EXIT
   ELSE

    CALL Exchange_Particle(nExSum)
    CALL Extract_Particle(mg_mesh%level(ILEV)%dcorvg,&
                          mg_mesh%level(ILEV)%kvert,&
                          mg_mesh%level(ILEV)%kedge,&
                          mg_mesh%level(ILEV)%karea,&
                          mg_mesh%level(ILEV)%elementsAtVertexIdx,&
                          mg_mesh%level(ILEV)%elementsAtVertex,&
                          QuadSc%ValU,QuadSc%ValV,QuadSc%ValW,&
                          nvt,net,nat,nel,dTime,myParticleParam%d_CorrDist)

   END IF

  END DO

  !!! here we need to collide the particles against each other
!   CALL CheckForCollision(nvt,net,nat,nel,myParticleParam%d_CorrDist)
  
  IF (bOutput.and.mod(iTimeSteps,iOutputFreq).eq.0) THEN
  
   iOutputIndex = iOutputIndex + 1

   CALL OutputParticlesCSV(iOutputIndex)

   bOutputLambda = .FALSE.
   IF (MOD(iTimeSteps-1,myParticleParam%nTimeLevels).EQ.0) bOutputLambda=.TRUE.
   iOutputLambda = iTimeSteps/myParticleParam%nTimeLevels

   ! not the end of the simulation - so call getLambda with false
   CALL GetLambda(.false.)

  END IF

  dStat = 0d0
  DO i=1,nBuffer
   dStat = dStat +  dBuffer(i)
  END DO
  dStat = dStat/REAL(nBuffer)

  CALL GetCutplanes()

  IF (myid.eq.1) WRITE(947,'(3(E12.4,A))') (REAL(iTimeSteps)/REAL(myParticleParam%nTimeLevels))/(myParticleParam%f/6d1),', ',REAL(nActSum0-nActSum)/REAL(nActSum0),', ',dStat/REAL(nActSum0),' '
  !IF (myid.eq.1) WRITE(947,'(3(E12.4,A))') REAL(iTimeSteps),', ',REAL(nActSum0-nActSum)/REAL(nActSum0),', ',dStat/REAL(nActSum0),' '
  iBuffer = MOD(iTimeSteps,nBuffer)+1
  dBuffer(iBuffer) = (nActSumOld-nActSum)
  nActSumOld = nActSum

  IF (DBLE(nActSum0-nActSum)/DBLE(nActSum0).ge.myParticleParam%minFrac) EXIT
  

END DO

IF (bOutput) THEN

 CALL OutputLostParticlesCSV()

END IF

CALL OutputParticlesAtZtoCSV()

! Now its the end of the simulation - so output the final lambda
call GetLambda(.true.)

IF (myid.eq.1) CLOSE(947)

! 222 CONTINUE
!
! CALL PostProcessRTD()

END SUBROUTINE Transport_Particle
!
! --------------------------------------------------------------------
!
SUBROUTINE BackTransport_Particle(mfile)
USE Sigma_User, ONLY: myMultiMat
INTEGER mfile
INTEGER i,nExSum,nActSum,nActSum0,nActSumOld,iCycle,iFile,iLevel0,iLevel1
REAL*8  dTime,daux,dPeriod,dTimeStep,dStart,dBuffer(nBuffer)
CHARACTER*99 cFile
CHARACTER*60 cInFile

LOGICAL :: bOutput=.TRUE.
real*8 xmin, xmax, ymin, ymax, zmin, zmax
INTEGER, ALLOCATABLE :: MatDist(:)

dTime=0d0
dPeriod = 6d1/myParticleParam%f
dTimeStep = dPeriod/DBLE(myParticleParam%nTimeLevels)

ALLOCATE(myVelo(1))

!  IF (.not.ALLOCATED(MaterialDistribution)) ALLOCATE(MaterialDistribution(1:NLMAX))
!  IF (.not.ALLOCATED(MaterialDistribution(NLMAX)%x)) ALLOCATE(MaterialDistribution(NLMAX)%x(mg_mesh%level(NLMAX)%nel))
!  MaterialDistribution(NLMAX)%x = 0
!  cInFile = "1"
!  CALL read_sol_from_file(cInFile,1,timens)
!  
!  CALL Output_Profiles(0)
  
myMatrixRenewal%D = 0
myMatrixRenewal%K = 0
myMatrixRenewal%S = 0
myMatrixRenewal%M = 0
myMatrixRenewal%C = 0

WRITE(cFile,'(I0)') myParticleParam%dump_in_file

!!!!!!!!!!!!!!!!!!! choose one fof these !!!!!!!!!!!!!!!!!!!!!!
if (myParticleParam%DumpFormat.eq.1) CALL init_sol_same_level(cFile)
if (myParticleParam%DumpFormat.eq.2) CALL Load_ListFiles_PRT_Tracer(myParticleParam%dump_in_file)
if (myParticleParam%DumpFormat.eq.3) CALL init_sol_repart(cFile)

ALLOCATE(myVelo(1)%x(QuadSc%ndof))
ALLOCATE(myVelo(1)%y(QuadSc%ndof))
ALLOCATE(myVelo(1)%z(QuadSc%ndof))
myVelo(1)%x = -QuadSc%ValU
myVelo(1)%y = -QuadSc%ValV
myVelo(1)%z = -QuadSc%ValW

IF (myid.ne.0) THEN
 ILEV=NLMAX
 CALL SETLEV(2)
END IF
nStartActiveSet = 0
CALL Extract_Particle(mg_mesh%level(ILEV)%dcorvg,&
                      mg_mesh%level(ILEV)%kvert,&
                      mg_mesh%level(ILEV)%kedge,&
                      mg_mesh%level(ILEV)%karea,&
                      mg_mesh%level(ILEV)%elementsAtVertexIdx,&
                      mg_mesh%level(ILEV)%elementsAtVertex,&
                      myVelo(1)%x,myVelo(1)%y,myVelo(1)%z,&
                      nvt,net,nat,nel,dTime,myParticleParam%d_CorrDist)

! Get the bounds of the mesh. Note: We need to set it to nlmax because else
! it won't work
IF (myid.ne.0) THEN
 ILEV=NLMAX
 CALL SETLEV(2)
END IF

call GetMeshBounds(mg_mesh%level(NLMAX)%dcorvg,nvt,xmin,xmax,ymin,ymax,zmin,zmax)
myMeshInfo%xmin = xmin
myMeshInfo%xmax = xmax
myMeshInfo%ymin = ymin
myMeshInfo%ymax = ymax
myMeshInfo%zmin = zmin
myMeshInfo%zmax = zmax
if (myid .eq. showid) then
  write(mfile,'(A,E12.4,A7,E12.4)') "MESH-Bounds in X-Dimension xmin: ",myMeshInfo%xmin, " xmax: ", myMeshInfo%xmax
  write(mfile,'(A,E12.4,A7,E12.4)') "MESH-Bounds in Y-Dimension ymin: ",myMeshInfo%ymin, " ymax: ", myMeshInfo%ymax
  write(mfile,'(A,E12.4,A7,E12.4)') "MESH-Bounds in Z-Dimension zmin: ",myMeshInfo%zmin, " zmax: ", myMeshInfo%zmax
  write(mterm,'(A,E12.4,A7,E12.4)') "MESH-Bounds in X-Dimension xmin: ",myMeshInfo%xmin, " xmax: ", myMeshInfo%xmax
  write(mterm,'(A,E12.4,A7,E12.4)') "MESH-Bounds in Y-Dimension ymin: ",myMeshInfo%ymin, " ymax: ", myMeshInfo%ymax
  write(mterm,'(A,E12.4,A7,E12.4)') "MESH-Bounds in Z-Dimension zmin: ",myMeshInfo%zmin, " zmax: ", myMeshInfo%zmax
end if ! Dump mesh-info

dBuffer = 0d0
nLostSet = 0

! Output initial positions of particles to csv.
IF (bOutput) THEN
 CALL OutputParticlesCSV(0)
END IF

DO iTimeSteps=1,myParticleParam%nRotation*myParticleParam%nTimeLevels

  nStartActiveSet = 0
  iCycle = 0
  iLevel0 = MOD(iTimeSteps-1,myParticleParam%nTimeLevels)+1
  iLevel1 = MOD(iTimeSteps  ,myParticleParam%nTimeLevels)+1

  dStart = dTime
  dTime = dTime + dTimeStep

  IF (myid.eq.1) WRITE(MFILE,'(A,I8,A,I8,A,I0,A,I0,ES14.6)') 'Timestep ',iTimeSteps, ' / ',myParticleParam%nRotation*myParticleParam%nTimeLevels,'  to be performed ...',iLevel0,'/',iLevel1,dTime
  IF (myid.eq.1) WRITE(MTERM,'(A,I8,A,I8,A,I0,A,I0,ES14.6)') 'Timestep ',iTimeSteps, ' / ',myParticleParam%nRotation*myParticleParam%nTimeLevels,'  to be performed ...',iLevel0,'/',iLevel1,dTime

  DO
   IF (myid.ne.0) THEN
    CALL Move_Particle(mg_mesh%level(ILEV)%dcorvg,&
                       mg_mesh%level(ILEV)%kvert,&
                       mg_mesh%level(ILEV)%kedge,&
                       mg_mesh%level(ILEV)%karea,&
                       mg_mesh%level(ILEV)%elementsAtVertexIdx,&
                       mg_mesh%level(ILEV)%elementsAtVertex,&
                       myVelo(1)%x,myVelo(1)%y,myVelo(1)%z,&
                       myVelo(1)%x,myVelo(1)%y,myVelo(1)%z,&
                       nvt,net,nat,nel,dTime,dStart,myParticleParam%d_CorrDist)
   END IF

   daux = DBLE(nActiveSet)
   CALL Comm_Summ(daux)
   nActSum =INT(daux)
   IF (iTimeSteps.EQ.1) nActSumOld = nActSum
   IF (iTimeSteps.EQ.1) nActSum0   = nActSum

   daux = DBLE(nExchangeSet)
   CALL Comm_Summ(daux)
   nExSum =INT(daux)

   iCycle = iCycle + 1
   IF (myid.eq.1) WRITE(MFILE,'(A,I0,A,I0,A,I0)') 'Exchange_of_particles: Cycle: ',iCycle, ' Num_Of_Points: ',nExSum,'/',nActSum
   IF (myid.eq.1) WRITE(MTERM,'(A,I0,A,I0,A,I0)') 'Exchange_of_particles: Cycle: ',iCycle, ' Num_Of_Points: ',nExSum,'/',nActSum


   IF (nExSum.EQ.0) THEN
    IF (myid.eq.1) WRITE(MFILE,'(A)') '...........................................................................'
    IF (myid.eq.1) WRITE(MTERM,'(A)') '...........................................................................'
    EXIT
   ELSE

    CALL Exchange_Particle(nExSum)
    CALL Extract_Particle(mg_mesh%level(ILEV)%dcorvg,&
                          mg_mesh%level(ILEV)%kvert,&
                          mg_mesh%level(ILEV)%kedge,&
                          mg_mesh%level(ILEV)%karea,&
                          mg_mesh%level(ILEV)%elementsAtVertexIdx,&
                          mg_mesh%level(ILEV)%elementsAtVertex,&
                          myVelo(1)%x,myVelo(1)%y,myVelo(1)%z,&
                          nvt,net,nat,nel,dTime,myParticleParam%d_CorrDist)

   END IF

  END DO

  IF (bOutput) THEN

   CALL OutputParticlesCSV(iTimeSteps)

  END IF

END DO

IF (bOutput) THEN

 myMultiMat%nOfMaterials = myParticleParam%NumberOfInflowRegions
 GenLinScalar%cName = "Temper"
 nMaterials = myMultiMat%nOfMaterials+1
 GenLinScalar%prm%nOfFields = nMaterials
 GenLinScalar%nOfFields = nMaterials
 ALLOCATE(GenLinScalar%prm%cField(GenLinScalar%prm%nOfFields))
 ALLOCATE(GenLinScalar%Fld(nMaterials))
 
 DO iFld=1,nMaterials
  if (iFld.eq.1) GenLinScalar%prm%cField(iFld) = 'temp'
  if (iFld.gt.1) then
   write(GenLinScalar%prm%cField(iFld),'(A,I0)') 'alpha',iFld-1
  end if
 end do

 DO iFld = 1,nMaterials
  GenLinScalar%Fld(iFld)%cName = GenLinScalar%prm%cField(iFld)
 END DO
  
 ALLOCATE(MatDist(myParticleParam%nParticles))
 MatDist = 0
 CALL AssignInflowPropertyToParticles(MatDist)
 CALL ExtractMatrialProperties(MatDist)
 
 if (myid.ne.0) THEN
 
  ILEV=NLMAX
  CALL SETLEV(2)
  
  ALLOCATE(GenLinScalar%Fld(1)%aux(QuadSc%ndof))
  ALLOCATE(GenLinScalar%Fld(1)%val(QuadSc%ndof))
  GenLinScalar%Fld(1)%val = myParticleParam%MeltTemperature
  
  DO iFld = 2,nMaterials
  
   ALLOCATE(GenLinScalar%Fld(iFld)%aux(QuadSc%ndof))
   ALLOCATE(GenLinScalar%Fld(iFld)%val(QuadSc%ndof))
   
   call INT_ParP0toQ2(MaterialDistribution(ILEV)%x,&
                      GenLinScalar%Fld(iFld)%val,&
                      GenLinScalar%Fld(iFld)%aux,&
                      mg_mesh%level(ILEV)%dvol,&
                      mg_mesh%level(ILEV)%kvert,&
                      mg_mesh%level(ILEV)%kedge,&
                      mg_mesh%level(ILEV)%karea,&
                      mg_mesh%level(ILEV)%nvt,&
                      mg_mesh%level(ILEV)%net,&
                      mg_mesh%level(ILEV)%nat,&
                      mg_mesh%level(ILEV)%nel,iFld)
  END DO
 END IF
 
 CALL OutputLostParticlesCSV()
 
 myExport%Level = NLMAX
 myExport%LevelMax = myExport%Level
 DEALLOCATE(myExport%Fields)
 ALLOCATE(myExport%Fields(2))
 myExport%Fields(2) = "Material_E"
 myExport%Fields(1) = "GenScalar"
 CALL Output_Profiles(0)
  
 CALL Release_ListFiles_General(0,'q')
 
!  CALL Release_ListFiles_General(0,'s')
 
!  CALL write_sol_to_file(10,0d0,0)

END IF

END SUBROUTINE BackTransport_Particle
!
!-------------------------------------------------------------------
!
SUBROUTINE ExtractMatrialProperties(m)

implicit none
real*8, allocatable :: separator(:),daux(:,:)
integer m(*)
integer i,j,iMat,iel,pID,ind,jnd
real*8 distMin,dist,point(3)

ILEV = NLMAX

allocate(separator(0:subnodes))

separator(:) = 0d0
separator(myid) = DBLE(mg_mesh%level(ILEV)%nel)

CALL Comm_SummN(separator,subnodes+1)

separator(0) = 0d0

do pID=1,subnodes
 separator(pID) = separator(pID) + separator(pID-1)
end do

ALLOCATE(daux         (3,myParticleParam%nParticles))
daux = 0d0

if (myid.ne.0) then
 DO i=1,mg_mesh%level(ILEV)%nel
   j = INT(separator(myid-1)) + i
   daux(:,j) = mg_mesh%level(ILEV)%dcorvg(:,nvt+net+nat+i)
 END DO
END IF

CALL Comm_SummN(daux,3*myParticleParam%nParticles)

IF (.not.ALLOCATED(MaterialDistribution)) ALLOCATE(MaterialDistribution(1:NLMAX))
IF (.not.ALLOCATED(MaterialDistribution(NLMAX)%x)) ALLOCATE(MaterialDistribution(NLMAX)%x(mg_mesh%level(NLMAX)%nel))
MaterialDistribution(NLMAX)%x = 0

! write(*,*) 'separator = ',separator

if (myid.ne.0) then
 do iel=1,mg_mesh%level(NLMAX)%nel

  ind = INT(separator(myid-1)) + iel
  
  if (m(ind).eq.0) then
   distMin = 1d30
   iMat = 0
   do jnd=1,myParticleParam%nParticles
    if (m(jnd).gt.0) then
      dist = sqrt((daux(1,jnd)-daux(1,ind))**2d0 + (daux(2,jnd)-daux(2,ind))**2d0 + (daux(3,jnd)-daux(3,ind))**2d0)
      if (distMin.gt.dist) then
       distMin = dist
       iMat = m(jnd)
      end if
    end if
   end do
  else
   iMat = m(ind)
  end if
  MaterialDistribution(NLMAX)%x(iel) = iMat
 end do
end if

END SUBROUTINE ExtractMatrialProperties
!
! --------------------------------------------------------------------
!
SUBROUTINE Init_Particle_xse(mfile)
INTEGER i,j,iParticles
REAL*8 myRandomNumber(3),R_min,R_max,Y,X_box,Y_box,Z_min
REAL*8 dBuff(3)
INTEGER iO
INTEGER inittype

call prt_read_config(myParticleParam, mfile, mterm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!! CorrDistEstimation !!!!!!!!!!!!!!!!!!!!!!!!!!!! 
IF (ADJUSTL(TRIM(mySigma%cType)).EQ."DIE") THEN
 myParticleParam%d_CorrDist = myProcess%ExtrusionGapSize*0.20d0*0.10d0 !(0.10d0 :: mm ==> cm), (0.10d0 :: safe min distance estimation for the outflow
 myParticleParam%Epsilon = myParticleParam%d_CorrDist*myParticleParam%Epsilon !(0.10d0 :: mm ==> cm), (0.10d0 :: safe min distance estimation for the outflow
 IF (myid.eq.1) THEN
 ! write(*,*) ADJUSTL(TRIM(mySigma%cType)), myParticleParam%d_CorrDist,myParticleParam%Epsilon
  write(*,*) "d_Corrdist and Epsilon is adaptively adjsuted due to DIE simulation", myParticleParam%d_CorrDist,myParticleParam%Epsilon
 end if
END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!! CorrDistEstimation !!!!!!!!!!!!!!!!!!!!!!!!!!!! 

IF (myid.eq.1) THEN
 WRITE(mfile,*) 'myParticleParam%dump_in_file = ',myParticleParam%dump_in_file
 WRITE(mfile,*) 'myParticleParam%nTimeLevels = ', myParticleParam%nTimeLevels
 WRITE(mfile,*) 'myParticleParam%nParticles  = ', myParticleParam%nParticles
 WRITE(mfile,*) 'myParticleParam%nRotation   = ', myParticleParam%nRotation
 WRITE(mfile,*) 'myParticleParam%d_CorrDist  = ', myParticleParam%d_CorrDist
 WRITE(mfile,*) 'myParticleParam%minFrac     = ', myParticleParam%minFrac
 WRITE(mfile,*) 'myParticleParam%Raster      = ', myParticleParam%Raster
 WRITE(mfile,*) 'myParticleParam%dEps1       = ', myParticleParam%dEps1
 WRITE(mfile,*) 'myParticleParam%dEps2       = ', myParticleParam%dEps2
 WRITE(mfile,*) 'myParticleParam%D_Out       = ', myParticleParam%D_Out
 WRITE(mfile,*) 'myParticleParam%D_In        = ', myParticleParam%D_In
 WRITE(mfile,*) 'myParticleParam%Z_seed      = ', myParticleParam%Z_seed
 WRITE(mfile,*) 'myParticleParam%OutflowZPos = ', myParticleParam%OutflowZPos
 WRITE(mfile,*) 'myParticleParam%RotMovement = ', myParticleParam%bRotationalMovement
 WRITE(mfile,*) 'myParticleParam%DumpFormat  = ', myParticleParam%DumpFormat
 WRITE(mfile,*) 'myParticleParam%Timelevels  = ', myParticleParam%nTimeLevels
 WRITE(mfile,*) 'myParticleParam%Periodicity = ', myParticleParam%nPeriodicity
 WRITE(mfile,*) 'myParticleParam%f           = ', myParticleParam%f
 WRITE(mfile,*) 'myParticleParam%Epsilon     = ', myParticleParam%Epsilon
 WRITE(mfile,*) 'myParticleParam%hSize       = ', myParticleParam%hSize
 WRITE(mterm,*) 'myParticleParam%dump_in_file = ', myParticleParam%dump_in_file
 WRITE(mterm,*) 'myParticleParam%nTimeLevels = ', myParticleParam%nTimeLevels
 WRITE(mterm,*) 'myParticleParam%nParticles  = ', myParticleParam%nParticles
 WRITE(mterm,*) 'myParticleParam%nRotation   = ', myParticleParam%nRotation
 WRITE(mterm,*) 'myParticleParam%d_CorrDist  = ', myParticleParam%d_CorrDist
 WRITE(mterm,*) 'myParticleParam%minFrac     = ', myParticleParam%minFrac
 WRITE(mterm,*) 'myParticleParam%Raster      = ', myParticleParam%Raster
 WRITE(mterm,*) 'myParticleParam%dEps1       = ', myParticleParam%dEps1
 WRITE(mterm,*) 'myParticleParam%dEps2       = ', myParticleParam%dEps2
 WRITE(mterm,*) 'myParticleParam%D_Out       = ', myParticleParam%D_Out
 WRITE(mterm,*) 'myParticleParam%D_In        = ', myParticleParam%D_In
 WRITE(mterm,*) 'myParticleParam%Z_seed      = ', myParticleParam%Z_seed
 WRITE(mterm,*) 'myParticleParam%OutflowZPos = ', myParticleParam%OutflowZPos
 WRITE(mterm,*) 'myParticleParam%RotMovement = ', myParticleParam%bRotationalMovement
 WRITE(mterm,*) 'myParticleParam%DumpFormat  = ', myParticleParam%DumpFormat
 WRITE(mterm,*) 'myParticleParam%Timelevels  = ', myParticleParam%nTimeLevels
 WRITE(mterm,*) 'myParticleParam%Periodicity = ', myParticleParam%nPeriodicity
 WRITE(mterm,*) 'myParticleParam%f           = ', myParticleParam%f
 WRITE(mterm,*) 'myParticleParam%Epsilon     = ', myParticleParam%Epsilon
 WRITE(mterm,*) 'myParticleParam%hSize       = ', myParticleParam%hSize
END IF

if (myParticleParam%inittype .eq. ParticleSeed_CSVFILE) then
 call Init_Particles_from_csv_xse(mfile)
else if(myParticleParam%inittype .eq. ParticleSeed_E3D) then
 CALL Init_Particles_viaE3D(mfile)
ELSE
 write(*,*) 'no valid initialization is chosen'
 STOP
END IF
! Now the part that is identical for all initialisations
if (myParticleParam%nZposCutplanes .gt. 0 ) then
  ALLOCATE(ZPosParticles(myParticleParam%nZposCutplanes,5,1:myParticleParam%nParticles))
  ! Initialisation with zero - always good to have known values in there
  ZPosParticles = 0.0d0
  ALLOCATE(bPosParticles(myParticleParam%nZposCutplanes,1:myParticleParam%nParticles))
  bPosParticles = .TRUE.
end if

END SUBROUTINE Init_Particle_xse
!
! --------------------------------------------------------------------
!
SUBROUTINE Init_Particle(mfile)
INTEGER i,j,iParticles
REAL*8 myRandomNumber(3),R_min,R_max,Y,X_box,Y_box,Z_min
REAL*8 dBuff(3)
INTEGER iO

call prt_read_config(myParticleParam, mfile, mterm)

IF (myid.eq.1) THEN
 WRITE(mfile,*) 'myParticleParam%dump_in_file = ',myParticleParam%dump_in_file
 WRITE(mfile,*) 'myParticleParam%nTimeLevels = ', myParticleParam%nTimeLevels
 WRITE(mfile,*) 'myParticleParam%nParticles  = ', myParticleParam%nParticles
 WRITE(mfile,*) 'myParticleParam%nRotation   = ', myParticleParam%nRotation
 WRITE(mfile,*) 'myParticleParam%minFrac     = ', myParticleParam%minFrac
 WRITE(mfile,*) 'myParticleParam%Raster      = ', myParticleParam%Raster
 WRITE(mfile,*) 'myParticleParam%dEps1       = ', myParticleParam%dEps1
 WRITE(mfile,*) 'myParticleParam%dEps2       = ', myParticleParam%dEps2
 WRITE(mfile,*) 'myParticleParam%D_Out       = ', myParticleParam%D_Out
 WRITE(mfile,*) 'myParticleParam%D_In        = ', myParticleParam%D_In
 WRITE(mfile,*) 'myParticleParam%Z_seed      = ', myParticleParam%Z_seed
 WRITE(mfile,*) 'myParticleParam%OutflowZPos = ', myParticleParam%OutflowZPos
 WRITE(mfile,*) 'myParticleParam%RotMovement = ', myParticleParam%bRotationalMovement
 WRITE(mfile,*) 'myParticleParam%DumpFormat  = ', myParticleParam%DumpFormat
 WRITE(mfile,*) 'myParticleParam%Timelevels  = ', myParticleParam%nTimeLevels
 WRITE(mfile,*) 'myParticleParam%Periodicity = ', myParticleParam%nPeriodicity
 WRITE(mfile,*) 'myParticleParam%f           = ', myParticleParam%f
 WRITE(mfile,*) 'myParticleParam%Epsilon     = ', myParticleParam%Epsilon
 WRITE(mfile,*) 'myParticleParam%hSize       = ', myParticleParam%hSize
 WRITE(mfile,*) 'myParticleParam%MeltTemperature   = ', myParticleParam%MeltTemperature
 WRITE(mfile,*) 'myParticleParam%nOfInflows  = ', myParticleParam%NumberOfInflowRegions
 IF (myParticleParam%NumberOfInflowRegions.gt.1) then
  do i=1,myParticleParam%NumberOfInflowRegions
   WRITE(mfile,'(A,I0,A,3ES12.4)') 'myParticleParam%InflowCenter',i,' = ', myParticleParam%InflowRegion(i)%Center
   WRITE(mfile,'(A,I0,A,1ES12.4)') 'myParticleParam%InflowRadius',i,' = ', myParticleParam%InflowRegion(i)%Radius
  end do
 END IF
 WRITE(mfile,*) 'myParticleParam%bPhysParticles  = ', myParticleParam%bPhysParticles
 IF (myParticleParam%bPhysParticles) then
  WRITE(mfile,*) 'myParticleParam%PhysParticles%rho_l    = ', myParticleParam%PhysParticles%rho_l
  WRITE(mfile,*) 'myParticleParam%PhysParticles%rho_p    = ', myParticleParam%PhysParticles%rho_p
  WRITE(mfile,*) 'myParticleParam%PhysParticles%mu_l     = ', myParticleParam%PhysParticles%mu_l
  WRITE(mfile,*) 'myParticleParam%PhysParticles%d_p      = ', myParticleParam%PhysParticles%d_p
  WRITE(mfile,*) 'myParticleParam%PhysParticles%gravity  = ', myParticleParam%PhysParticles%gravity
 END IF

 WRITE(mterm,*) 'myParticleParam%dump_in_file = ', myParticleParam%dump_in_file
 WRITE(mterm,*) 'myParticleParam%nTimeLevels = ', myParticleParam%nTimeLevels
 WRITE(mterm,*) 'myParticleParam%nParticles  = ', myParticleParam%nParticles
 WRITE(mterm,*) 'myParticleParam%nRotation   = ', myParticleParam%nRotation
 WRITE(mterm,*) 'myParticleParam%d_CorrDist  = ', myParticleParam%d_CorrDist
 WRITE(mterm,*) 'myParticleParam%minFrac     = ', myParticleParam%minFrac
 WRITE(mterm,*) 'myParticleParam%Raster      = ', myParticleParam%Raster
 WRITE(mterm,*) 'myParticleParam%dEps1       = ', myParticleParam%dEps1
 WRITE(mterm,*) 'myParticleParam%dEps2       = ', myParticleParam%dEps2
 WRITE(mterm,*) 'myParticleParam%D_Out       = ', myParticleParam%D_Out
 WRITE(mterm,*) 'myParticleParam%D_In        = ', myParticleParam%D_In
 WRITE(mterm,*) 'myParticleParam%Z_seed      = ', myParticleParam%Z_seed
 WRITE(mterm,*) 'myParticleParam%OutflowZPos = ', myParticleParam%OutflowZPos
 WRITE(mterm,*) 'myParticleParam%RotMovement = ', myParticleParam%bRotationalMovement
 WRITE(mterm,*) 'myParticleParam%DumpFormat  = ', myParticleParam%DumpFormat
 WRITE(mterm,*) 'myParticleParam%Timelevels  = ', myParticleParam%nTimeLevels
 WRITE(mterm,*) 'myParticleParam%Periodicity = ', myParticleParam%nPeriodicity
 WRITE(mterm,*) 'myParticleParam%f           = ', myParticleParam%f
 WRITE(mterm,*) 'myParticleParam%Epsilon     = ', myParticleParam%Epsilon
 WRITE(mterm,*) 'myParticleParam%hSize       = ', myParticleParam%hSize
 WRITE(mterm,*) 'myParticleParam%MeltTemperature   = ', myParticleParam%MeltTemperature
 WRITE(mterm,*) 'myParticleParam%nOfInflows  = ', myParticleParam%NumberOfInflowRegions
 IF (myParticleParam%NumberOfInflowRegions.gt.1) then
  do i=1,myParticleParam%NumberOfInflowRegions
   WRITE(mterm,'(A,I0,A,3ES12.4)') 'myParticleParam%InflowCenter',i,' = ', myParticleParam%InflowRegion(i)%Center
   WRITE(mterm,'(A,I0,A,1ES12.4)') 'myParticleParam%InflowRadius',i,' = ', myParticleParam%InflowRegion(i)%Radius
  end do
 END IF
 WRITE(mterm,*) 'myParticleParam%bPhysParticles  = ', myParticleParam%bPhysParticles
 IF (myParticleParam%bPhysParticles) then
  WRITE(mterm,*) 'myParticleParam%PhysParticles%rho_l    = ', myParticleParam%PhysParticles%rho_l
  WRITE(mterm,*) 'myParticleParam%PhysParticles%rho_p    = ', myParticleParam%PhysParticles%rho_p
  WRITE(mterm,*) 'myParticleParam%PhysParticles%mu_l     = ', myParticleParam%PhysParticles%mu_l
  WRITE(mterm,*) 'myParticleParam%PhysParticles%d_p      = ', myParticleParam%PhysParticles%d_p
  WRITE(mterm,*) 'myParticleParam%PhysParticles%gravity  = ', myParticleParam%PhysParticles%gravity
 END IF
END IF

if (myParticleParam%inittype .eq. ParticleSeed_Parameterfile ) then
  call Init_Particles_from_parameterfile(mfile)
else if(myParticleParam%inittype .eq. ParticleSeed_CSVFILE) then
  call Init_Particles_from_csv(mfile)
else if(myParticleParam%inittype .eq. ParticleSeed_OUTPUTFILE) then
  call Init_Particles_from_old_output(mfile)
else if(myParticleParam%inittype .eq. ParticleSeed_ELEMCENTER) then
  call Init_Particles_in_all_elemCenters(mfile)
else if(myParticleParam%inittype .eq. ParticleSeed_PLANE) then
  call Init_Particles_inPlane(mfile)
else if(myParticleParam%inittype .eq. ParticleSeed_Volume) then
  call Init_Particles_inVolume(mfile)
else
  write(mterm,*) "ERROR: Unknown Starting procedure ", myParticleParam%inittype
  write(mfile,*) "ERROR: Unknown starting procedure ", myParticleParam%inittype
  stop
end if


! Now the part that is identical for all initialisations
if (myParticleParam%nZposCutplanes .gt. 0 ) then
  ALLOCATE(ZPosParticles(myParticleParam%nZposCutplanes,5,1:myParticleParam%nParticles))
  ! Initialisation with zero - always good to have known values in there
  ZPosParticles = 0.0d0
  ALLOCATE(bPosParticles(myParticleParam%nZposCutplanes,1:myParticleParam%nParticles))
  bPosParticles = .TRUE.
end if

END SUBROUTINE Init_Particle
!
!-------------------------------------------------------------------
!
SUBROUTINE Init_Particles_inVolume(mfile)
INTEGER nfile
REAL*8 xmax,ymax,zmax,xmin,ymin,zmin
REAL*8 myRandomNumber(3),xxmin,yymin,xxmax,yymax,myRandomX,myRandomY,myRandomZ

xmax = -1d30
ymax = -1d30
zmax = -1d30

xmin = +1d30
ymin = +1d30
zmin = +1d30

ILEV = NLMAX

DO i=1,mg_mesh%level(ILEV)%nvt
 if (xmin > mg_mesh%level(ILEV)%dcorvg(1,i)) xmin = mg_mesh%level(ILEV)%dcorvg(1,i)
 if (xmax < mg_mesh%level(ILEV)%dcorvg(1,i)) xmax = mg_mesh%level(ILEV)%dcorvg(1,i)
 
 if (ymin > mg_mesh%level(ILEV)%dcorvg(2,i)) ymin = mg_mesh%level(ILEV)%dcorvg(2,i)
 if (ymax < mg_mesh%level(ILEV)%dcorvg(2,i)) ymax = mg_mesh%level(ILEV)%dcorvg(2,i)
 
 if (zmin > mg_mesh%level(ILEV)%dcorvg(3,i)) zmin = mg_mesh%level(ILEV)%dcorvg(3,i)
 if (zmax < mg_mesh%level(ILEV)%dcorvg(3,i)) zmax = mg_mesh%level(ILEV)%dcorvg(3,i)
end do

call COMM_Maximum(xmax)
call COMM_Maximum(ymax)
call COMM_Maximum(zmax)
call COMM_Minimum(xmin)
call COMM_Minimum(ymin)
call COMM_Minimum(zmin)

iElement = 0

numberOfFluidElements = myParticleParam%VolumeParticles
myParticleParam%nParticles = myParticleParam%VolumeParticles*4

ALLOCATE(myCompleteSet(1:myParticleParam%nParticles))
ALLOCATE(myActiveSet  (1:myParticleParam%nParticles))
ALLOCATE(myExchangeSet(1:myParticleParam%nParticles))
ALLOCATE(myLostSet    (1:myParticleParam%nParticles))

ALLOCATE(Lambda(numberOfFluidElements),Length(myParticleParam%nParticles,3))
ALLOCATE(pPresent(myParticleParam%nParticles))

 write(*,'(6A12)')   'xmax','ymax','zmax','xmin','ymin','zmin'
 write(*,'(6ES12.4)') xmax,ymax,zmax,xmin,ymin,zmin
 write(*,*) 'n= ',myParticleParam%VolumeParticles

if (myid.eq.1) then
 write(*,'(6A12)')   'xmax','ymax','zmax','xmin','ymin','zmin'
 write(*,'(6ES12.4)') xmax,ymax,zmax,xmin,ymin,zmin
END IF

do 

 myRandomX = xmin + (xmax-xmin)*rand(0)
 myRandomY = ymin + (ymax-ymin)*rand(0)
 myRandomZ = zmin + (zmax-zmin)*rand(0)

 myRandomNumber = [myRandomX,myRandomY,myRandomZ]
 
 cdx = 1d1*myRandomNumber(1)
 cdy = 1d1*myRandomNumber(2)
 cdz = 1d1*myRandomNumber(3)

 CALL GetDistAndProjPToAllSTLs(cdx,cdy,cdz,cpx,cpy,cpz,dist_CGAL)
 
 if (dist_CGAL.gt.myParticleParam%d_CorrDist) then
!  CALL GetDistToSTL(cdx,cdy,cdz,1,dist_CGAL,.true.)
!  if ((     myParticleParam%bRotationalMovement.and.dist_CGAL.gt.myParticleParam%d_CorrDist).or.&
!      (.not.myParticleParam%bRotationalMovement.and.dist_CGAL.lt.myParticleParam%d_CorrDist)) then
    
   iElement = iElement+ 1
   
   ! Calculate the indices of the particles that we
   ! are generating. We create them elementwise:
   ! First those from element 1, then those from
   ! element two, then those from element 3, ...
   idx1 = 4*iElement - 3
   idx2 = 4*iElement - 2
   idx3 = 4*iElement - 1
   idx4 = 4*iElement - 0

   if (myid.eq.1) write(*,*) myRandomNumber(1),myRandomNumber(2),myRandomNumber(3),iElement
   
   myExchangeSet(idx1)%coor(:) = [myRandomNumber(1),myRandomNumber(2),myRandomNumber(3)]
   myExchangeSet(idx1)%time    = 0d0
   myExchangeSet(idx1)%indice  = idx1

   myExchangeSet(idx2)%coor(:) = [myRandomNumber(1)+myParticleParam%Epsilon,myRandomNumber(2),myRandomNumber(3)]
   myExchangeSet(idx2)%time    = 0d0
   myExchangeSet(idx2)%indice  = idx2

   myExchangeSet(idx3)%coor(:) = [myRandomNumber(1),myRandomNumber(2)+myParticleParam%Epsilon,myRandomNumber(3)]
   myExchangeSet(idx3)%time    = 0d0
   myExchangeSet(idx3)%indice  = idx3

   myExchangeSet(idx4)%coor(:) = [myRandomNumber(1),myRandomNumber(2),myRandomNumber(3)+myParticleParam%Epsilon]
   myExchangeSet(idx4)%time    = 0d0
   myExchangeSet(idx4)%indice  = idx4
   
 END IF
 
 if (iElement.ge.myParticleParam%VolumeParticles) exit
 
END DO

write(*,*) 'iElement=',iElement
CLOSE(745)

nExchangeSet = myParticleParam%nParticles

END SUBROUTINE Init_Particles_inVolume
!
!-------------------------------------------------------------------
!
SUBROUTINE Init_Particles_inPlane(mfile)
INTEGER nfile
REAL*8 xmax,ymax,zmax,xmin,ymin,zmin
REAL*8 myRandomNumber(3),xxmin,yymin,xxmax,yymax,myRandomX,myRandomY
integer i,idx1,idx2,idx3,idx4,iElement,numberOfFluidElements 
REAL*8 cdX0,cdY0,cdZ0,dist_CGAL
logical :: bProjection=.false.
integer iProjection

xmax = -1d30
ymax = -1d30
zmax = -1d30

xmin = +1d30
ymin = +1d30
zmin = +1d30

ILEV = NLMAX

DO i=1,mg_mesh%level(ILEV)%nvt
 if (xmin > mg_mesh%level(ILEV)%dcorvg(1,i)) xmin = mg_mesh%level(ILEV)%dcorvg(1,i)
 if (xmax < mg_mesh%level(ILEV)%dcorvg(1,i)) xmax = mg_mesh%level(ILEV)%dcorvg(1,i)
 
 if (ymin > mg_mesh%level(ILEV)%dcorvg(2,i)) ymin = mg_mesh%level(ILEV)%dcorvg(2,i)
 if (ymax < mg_mesh%level(ILEV)%dcorvg(2,i)) ymax = mg_mesh%level(ILEV)%dcorvg(2,i)
 
 if (zmin > mg_mesh%level(ILEV)%dcorvg(3,i)) zmin = mg_mesh%level(ILEV)%dcorvg(3,i)
 if (zmax < mg_mesh%level(ILEV)%dcorvg(3,i)) zmax = mg_mesh%level(ILEV)%dcorvg(3,i)
end do


call COMM_Maximum(xmax)
call COMM_Maximum(ymax)
call COMM_Maximum(zmax)
call COMM_Minimum(xmin)
call COMM_Minimum(ymin)
call COMM_Minimum(zmin)

IF (myParticleParam%Plane.eq.1) then
 xxmin = ymin
 xxmax = ymax
 yymin = zmin
 yymax = zmax
END IF

IF (myParticleParam%Plane.eq.2) then
 xxmin = xmin
 xxmax = xmax
 yymin = zmin
 yymax = zmax
END IF

IF (myParticleParam%Plane.eq.3) then
 xxmin = xmin
 xxmax = xmax
 yymin = ymin
 yymax = ymax
END IF

iElement = 0

numberOfFluidElements = myParticleParam%PlaneParticles
myParticleParam%nParticles = myParticleParam%PlaneParticles*4

ALLOCATE(myCompleteSet(1:myParticleParam%nParticles))
ALLOCATE(myActiveSet  (1:myParticleParam%nParticles))
ALLOCATE(myExchangeSet(1:myParticleParam%nParticles))
ALLOCATE(myLostSet    (1:myParticleParam%nParticles))

ALLOCATE(Lambda(numberOfFluidElements),Length(myParticleParam%nParticles,3))
ALLOCATE(pPresent(myParticleParam%nParticles))

if (myid.eq.1) then
write(*,'(6A12)')   'xmax','ymax','zmax','xmin','ymin','zmin'
write(*,'(6ES12.4)') xmax,ymax,zmax,xmin,ymin,zmin
write(*,'(6ES12.4)') xxmax,yymax,xxmin,yymin
END IF

do 

 myRandomX = xxmin + (xxmax-xxmin)*rand(0)
 myRandomY = yymin + (yymax-yymin)*rand(0)

 IF (myParticleParam%Plane.eq.1)  myRandomNumber = [myParticleParam%PlaneOffSet,myRandomX,myRandomY]
 IF (myParticleParam%Plane.eq.2)  myRandomNumber = [myRandomX,myParticleParam%PlaneOffSet,myRandomY]
 IF (myParticleParam%Plane.eq.3)  myRandomNumber = [myRandomX,myRandomY,myParticleParam%PlaneOffSet]
 
 cdx = 1d1*myRandomNumber(1)
 cdy = 1d1*myRandomNumber(2)
 cdz = 1d1*myRandomNumber(3)
 
 CALL GetDistAndProjPToAllSTLs(cdx,cdy,cdz,cpx,cpy,cpz,dist_CGAL)
 
 if (dist_CGAL.gt.myParticleParam%d_CorrDist) then
 
!  CALL GetDistToSTL(cdx,cdy,cdz,1,dist_CGAL,.true.)
!  if ((     myParticleParam%bRotationalMovement.and.dist_CGAL.gt.myParticleParam%d_CorrDist).or.&
!      (.not.myParticleParam%bRotationalMovement.and.dist_CGAL.lt.myParticleParam%d_CorrDist)) then
!     
   iElement = iElement+ 1
   
   ! Calculate the indices of the particles that we
   ! are generating. We create them elementwise:
   ! First those from element 1, then those from
   ! element two, then those from element 3, ...
   idx1 = 4*iElement - 3
   idx2 = 4*iElement - 2
   idx3 = 4*iElement - 1
   idx4 = 4*iElement - 0

   if (myid.eq.1) write(*,*) myRandomNumber(1),myRandomNumber(2),myRandomNumber(3),iElement
   
   myExchangeSet(idx1)%coor(:) = [myRandomNumber(1),myRandomNumber(2),myRandomNumber(3)]
   myExchangeSet(idx1)%time    = 0d0
   myExchangeSet(idx1)%indice  = idx1

   myExchangeSet(idx2)%coor(:) = [myRandomNumber(1)+myParticleParam%Epsilon,myRandomNumber(2),myRandomNumber(3)]
   myExchangeSet(idx2)%time    = 0d0
   myExchangeSet(idx2)%indice  = idx2

   myExchangeSet(idx3)%coor(:) = [myRandomNumber(1),myRandomNumber(2)+myParticleParam%Epsilon,myRandomNumber(3)]
   myExchangeSet(idx3)%time    = 0d0
   myExchangeSet(idx3)%indice  = idx3

   myExchangeSet(idx4)%coor(:) = [myRandomNumber(1),myRandomNumber(2),myRandomNumber(3)+myParticleParam%Epsilon]
   myExchangeSet(idx4)%time    = 0d0
   myExchangeSet(idx4)%indice  = idx4
   
  END IF
 
 if (iElement.ge.myParticleParam%PlaneParticles) exit
 
END DO

write(*,*) iElement
CLOSE(745)

nExchangeSet = myParticleParam%nParticles

END SUBROUTINE Init_Particles_inPlane
!
!-------------------------------------------------------------------
!
SUBROUTINE Init_Particles_in_all_elemCenters(mfile)

real*8, allocatable :: separator(:),daux(:,:)
integer i

ILEV = NLMAX

allocate(separator(0:subnodes))

separator(:) = 0d0
separator(myid) = DBLE(mg_mesh%level(ILEV)%nel)

CALL Comm_SummN(separator,subnodes+1)

separator(0) = 0d0

do pID=1,subnodes
 separator(pID) = separator(pID) + separator(pID-1)
end do

myParticleParam%nParticles = INT(separator(subnodes))

ALLOCATE(myCompleteSet(1:myParticleParam%nParticles))
ALLOCATE(myActiveSet  (1:myParticleParam%nParticles))
ALLOCATE(myExchangeSet(1:myParticleParam%nParticles))
ALLOCATE(myLostSet    (1:myParticleParam%nParticles))
ALLOCATE(daux         (3,myParticleParam%nParticles))
daux                     = 0d0

DO i=1,myParticleParam%nParticles
 myExchangeSet(i)%coor(:) = 0d0
 myExchangeSet(i)%time    = 0d0
 myExchangeSet(i)%indice  = i
END DO

if (myid.ne.0) then
 DO i=1,mg_mesh%level(ILEV)%nel

   j = INT(separator(myid-1)) + i
   daux(:,j) = mg_mesh%level(ILEV)%dcorvg(:,nvt+net+nat+i)
   myExchangeSet(j)%time    = 0d0
   myExchangeSet(j)%indice  = j

 END DO
END IF

CALL Comm_SummN(daux,3*myParticleParam%nParticles)

 DO i=1,myParticleParam%nParticles

   myExchangeSet(i)%coor(:) = daux(:,i)

 END DO

nExchangeSet =  myParticleParam%nParticles

deallocate(separator,daux)

! if (myid.eq.1) then
!  open(file='start.csv',unit=1278)
!  WRITE(1278,'(5A)') '"coor_X",','"coor_Y",','"coor_Z",', '"indice",','"time"'
!  DO i=1,myParticleParam%nParticles
! 
!   ! Now output the particles to the file
!    WRITE(1278,'(3(E16.7,A),I0,A,E16.7)') REAL(myExchangeSet(i)%coor(1)*myParticleParam%dFacUnitOut),',',&
!                                 REAL(myExchangeSet(i)%coor(2)*myParticleParam%dFacUnitOut),',',&
!                                 REAL(myExchangeSet(i)%coor(3)*myParticleParam%dFacUnitOut),',',&
!                                 myExchangeSet(i)%indice,',',REAL(myExchangeSet(i)%time)
! !  write(1278,*) myExchangeSet(i)%coor(1)i
! 
!  END DO
!  close(1278)
! end if

if (myid.eq.1) WRITE(*,*) "myParticleParam%nParticles = ",myParticleParam%nParticles 
!  pause

END SUBROUTINE Init_Particles_in_all_elemCenters
!
!-------------------------------------------------------------------
!
SUBROUTINE Init_Particles_from_parameterfile(mfile)
INTEGER i,j,iElement, iSegment, iElementInSegment
REAL*8 myRandomNumber(3),R_min,R_max,Y,X_box,Y_box,Z_min
REAL*8 dBuff(3)
INTEGER iO
real*8 :: myRandomRadius, myRandomAngle, myNumberOfsegments, particlesPerSegment
real*8 :: stepAngles, langle, uangle
real*8 :: dParticle
real*8 :: myPi = 4.D0*DATAN(1.D0)
integer :: numberOfFluidElements, numberOfFluidElementsPerSegment
integer :: idx1, idx2, idx3, idx4

! Recalculate how many elements we have:
numberOfFluidElements = myParticleParam%nParticles/4

myNumberOfSegments = myParticleParam%numberSegments
! Check if number of segments and elements fit to each other
! This means that in each segment we have the same amout of elements
! If it does not work out up the number of elements and recalculate the
! number of particles
do while (mod(numberOfFluidElements,int(myNumberOfSegments)) .ne. 0)
  numberOfFluidElements = numberOfFluidElements + 1.0d0
end do
myParticleParam%nParticles=numberOfFluidElements*4.0d0

if (myid .eq. 1) write(*,'(A,I0,A)') 'initialised arrays for ', myParticleParam%nParticles, ' particles'



ALLOCATE(myCompleteSet(1:myParticleParam%nParticles))
ALLOCATE(myActiveSet  (1:myParticleParam%nParticles))
ALLOCATE(myExchangeSet(1:myParticleParam%nParticles))
ALLOCATE(myLostSet    (1:myParticleParam%nParticles))

R_max = 0.5d0*myParticleParam%D_Out*myParticleParam%dEps2
R_min = 0.5d0*myParticleParam%D_in
X_box = 0.5d0*myParticleParam%D_Out*myParticleParam%dEps2
Y_box = 0.5d0*myParticleParam%D_Out*myParticleParam%dEps2
Z_Min = myParticleParam%Z_seed

ALLOCATE(Lambda(numberOfFluidElements),Length(myParticleParam%nParticles,3))
ALLOCATE(pPresent(myParticleParam%nParticles))

Lambda = -1d0

stepAngles = 2.0d0*myPi / myNumberOfSegments

! This works because we adjusted the number of fluid elements above
numberOfFluidElementsPerSegment = numberOfFluidElements/int(myNumberOfSegments)
iElement = 0
DO iSegment=1,int(myNumberOfSegments)
  ! Calculate the lower and upper bound for the angle
  langle = 0.0d0 + real(iSegment -1)*stepAngles
  uangle =  0.0d0 + real(iSegment)*stepAngles

  do iElementInSegment=1,numberOfFluidElementsPerSegment
    ! We create an element
    iElement = iElement+1
    ! Now create a random angle inbetween these bounds
    ! This works as rand(0) returns a random number in
    ! the range (0,1)
    myRandomAngle = langle + (uangle - langle)*rand(0)
    myRandomRadius = R_min + (R_max - R_min) * rand(0)
    ! This leads to the following point
    myRandomNumber = [myRandomRadius*cos(myRandomAngle), &
                      myRandomRadius*sin(myRandomAngle), &
                      Z_Min]
    ! Calculate the indices of the particles that we
    ! are generating. We create them elementwise:
    ! First those from element 1, then those from
    ! element two, then those from element 3, ...
    idx1 = 4*iElement - 3
    idx2 = 4*iElement - 2
    idx3 = 4*iElement - 1
    idx4 = 4*iElement - 0

    ! Now set the coordinates
    ! This way, the set is divided into four blocks:
    ! First all "initial" points, then all with a delta_X,
    ! then all with a delta_Y, then all with a delta_Z
    myExchangeSet(idx1)%coor(:) = [myRandomNumber(1),myRandomNumber(2),myRandomNumber(3)]
    myExchangeSet(idx1)%time    = 0d0
    myExchangeSet(idx1)%indice  = idx1

    myExchangeSet(idx2)%coor(:) = [myRandomNumber(1)+myParticleParam%Epsilon,myRandomNumber(2),myRandomNumber(3)]
    myExchangeSet(idx2)%time    = 0d0
    myExchangeSet(idx2)%indice  = idx2

    myExchangeSet(idx3)%coor(:) = [myRandomNumber(1),myRandomNumber(2)+myParticleParam%Epsilon,myRandomNumber(3)]
    myExchangeSet(idx3)%time    = 0d0
    myExchangeSet(idx3)%indice  = idx3

    myExchangeSet(idx4)%coor(:) = [myRandomNumber(1),myRandomNumber(2),myRandomNumber(3)+myParticleParam%Epsilon]
    myExchangeSet(idx4)%time    = 0d0
    myExchangeSet(idx4)%indice  = idx4
  end do ! Loop over elements in segment

END DO ! Loop over all segments

nExchangeSet = myParticleParam%nParticles

END SUBROUTINE Init_Particles_from_parameterfile
!
!-------------------------------------------------------------------
!
SUBROUTINE Init_Particles_viaE3D(mfile)
INTEGER nfile
REAL*8 xmax,ymax,zmax,xmin,ymin,zmin
REAL*8 myRandomNumber(3),myRandomX,myRandomY
integer i,idx1,idx2,idx3,idx4,iElement,numberOfFluidElements,iPElement
REAL*8 cdX0,cdY0,cdZ0,dist_CGAL
logical :: bProjection=.false.
integer iProjection, iInflow, iTryDir, iFoundDir
REAL*8 dCenter(3), dNormal(3),dTry(3),dRefNormal1(3),dRefNormal2(3),dLine(3),dPlaneConstant,dOuterRadius,dDist,dTryDist,ProjPdist
REAL*8 tX,tY
real*8 :: myPi = 4.D0*DATAN(1.D0)


ILEV = NLMAX

numberOfFluidElements = myProcess%nOfInflows*myParticleParam%PlaneParticles
myParticleParam%nParticles = myProcess%nOfInflows*myParticleParam%PlaneParticles*4

ALLOCATE(myCompleteSet(1:myParticleParam%nParticles))
ALLOCATE(myActiveSet  (1:myParticleParam%nParticles))
ALLOCATE(myExchangeSet(1:myParticleParam%nParticles))
ALLOCATE(myLostSet    (1:myParticleParam%nParticles))

ALLOCATE(Lambda(numberOfFluidElements),Length(myParticleParam%nParticles,3))
ALLOCATE(pPresent(myParticleParam%nParticles))

iElement = 0
 
!  pause
DO iInflow = 1,myProcess%nOfInflows

  dCenter       = myProcess%myInflow(iInflow)%center
  dNormal       = myProcess%myInflow(iInflow)%normal
  dDist = ((dNormal(1))**2d0 + (dNormal(2))**2d0 + (dNormal(3))**2d0 )**0.5d0
  dNormal = dNormal/dDist
  dOuterRadius  = myProcess%myInflow(iInflow)%outerradius
  dCenter = dCenter + 0.0125d0*dOuterRadius*dNormal

  ! plane equation
  dPlaneConstant = dNormal(1)*dCenter(1) + dNormal(2)*dCenter(2) + dNormal(3)*dCenter(3)
  dTryDist = 1d-8
  iFoundDir = 0
  do iTryDir=1,3
   dTry = dCenter
   dTry(iTryDir) = dTry(iTryDir) + dOuterRadius
   
   ProjPdist = dNormal(1)*dTry(1) + dNormal(2)*dTry(2) + dNormal(3)*dTry(3) - dPlaneConstant
   dTry(1) = dTry(1) - ProjPdist*dNormal(1)
   dTry(2) = dTry(2) - ProjPdist*dNormal(2)
   dTry(3) = dTry(3) - ProjPdist*dNormal(3)
   ProjPdist = dNormal(1)*dTry(1) + dNormal(2)*dTry(2) + dNormal(3)*dTry(3) - dPlaneConstant
   
   dLine = (dTry - dCenter)
   dDist = ((dLine(1))**2d0 + (dLine(2))**2d0 + (dLine(3))**2d0)**0.5d0
   
   if (dTryDist.lt.dDist) then
    dTryDist = dDist
    dRefNormal1 = dLine/dDist
   end if
  end do
  
  dTry = dCenter
  dTry(iFoundDir) = dTry(iFoundDir) + dOuterRadius
  
  dRefNormal2 = [ dRefNormal1(2)*dNormal(3) - dRefNormal1(3)*dNormal(2),&
                 -dRefNormal1(1)*dNormal(3) + dRefNormal1(3)*dNormal(1),&
                  dRefNormal1(1)*dNormal(2) - dRefNormal1(2)*dNormal(1)]
  dDist = ((dRefNormal2(1))**2d0 + (dRefNormal2(2))**2d0 + (dRefNormal2(3))**2d0 )**0.5d0
  dRefNormal2 = dRefNormal2 / dDist
  
  iPElement = 0
  do 
   
   tX     =  (-1d0 + 2d0*rand(0))*dOuterRadius
   tAngle =  (-myPI + 2d0*myPI*rand(0))
   tY     =  dOuterRadius*sin(tAngle)
   dDist = sqrt(tX**2d0 + tY**2d0)
   IF (dDist.lt.dOuterRadius) then
   
      myRandomNumber = dCenter + tX*dRefNormal1 + tY*dRefNormal2
   
      cdx = myParticleParam%dFacUnitSourcefile*myRandomNumber(1)
      cdy = myParticleParam%dFacUnitSourcefile*myRandomNumber(2)
      cdz = myParticleParam%dFacUnitSourcefile*myRandomNumber(3)

      CALL SearchPointsWRT_STLs(cdx,cdy,cdz,0d0,dist_CGAL,bProjection,cdX0,cdY0,cdZ0,iProjection)
      if (dist_CGAL.gt.1d0*myParticleParam%d_CorrDist) then
         
         iElement = iElement+ 1
         iPElement = iPElement+ 1
         if (myid.eq.1) write(*,'(3ES12.4,(" "),I0)') cdx,cdy,cdz,iElement
         
         ! Calculate the indices of the particles that we
         ! are generating. We create them elementwise:
         ! First those from element 1, then those from
         ! element two, then those from element 3, ...
         idx1 = 4*iElement - 3
         idx2 = 4*iElement - 2
         idx3 = 4*iElement - 1
         idx4 = 4*iElement - 0

         myExchangeSet(idx1)%coor(:) = [myRandomNumber(1),myRandomNumber(2),myRandomNumber(3)]
         myExchangeSet(idx1)%velo    = 0d0
         myExchangeSet(idx1)%time    = 0d0
         myExchangeSet(idx1)%indice  = idx1
         myExchangeSet(idx1)%id = iInflow

         myExchangeSet(idx2)%coor(:) = [myRandomNumber(1)+myParticleParam%Epsilon,myRandomNumber(2),myRandomNumber(3)]
         myExchangeSet(idx2)%velo    = 0d0
         myExchangeSet(idx2)%time    = 0d0
         myExchangeSet(idx2)%indice  = idx2
         myExchangeSet(idx2)%id = iInflow

         myExchangeSet(idx3)%coor(:) = [myRandomNumber(1),myRandomNumber(2)+myParticleParam%Epsilon,myRandomNumber(3)]
         myExchangeSet(idx3)%velo    = 0d0
         myExchangeSet(idx3)%time    = 0d0
         myExchangeSet(idx3)%indice  = idx3
         myExchangeSet(idx3)%id = iInflow

         myExchangeSet(idx4)%coor(:) = [myRandomNumber(1),myRandomNumber(2),myRandomNumber(3)+myParticleParam%Epsilon]
         myExchangeSet(idx4)%velo    = 0d0
         myExchangeSet(idx4)%time    = 0d0
         myExchangeSet(idx4)%indice  = idx4
         myExchangeSet(idx4)%id = iInflow
     
      END IF
   
    END IF
    
   if (ipElement.ge.myParticleParam%PlaneParticles) exit
   
  END DO
  
 
END DO !iInflow

nExchangeSet = myParticleParam%nParticles

END SUBROUTINE Init_Particles_viaE3D
!
!-------------------------------------------------------------------
!
SUBROUTINE Init_Particles_from_csv_xse(mfile)
USE PP3D_MPI, ONLY : myid,master
use Sigma_User, only: mySigma, myProcess

implicit none
INTEGER i,j,iParticles, iElements, iElement, mfile
integer idx1,idx2,idx3,idx4
REAL*8 myRandomNumber(3),R_min,R_max,Y,X_box,Y_box,Z_min,dist
REAL*8 dBuff(3),cdX0,cdY0,cdZ0

Character*256 cHeader
INTEGER nColumns,Columns(4)
REAL*8, allocatable :: dValues(:)
INTEGER iO
logical :: bProjection=.false.
integer iProjection

! Count how many particles are in the sourcefile!
! We allocate one element per particle that we find
iParticles = 0
! if (myid.eq.14) open(785,file='points.csv')
! if (myid.eq.14) WRITE(785,'(9A)') '"coor_X",','"coor_Y",','"coor_Z",', '"d1",','"d2",', '"d3",','"d"'

myProcess%Angle = 0d0

OPEN(FILE=myParticleParam%sourcefile,UNIT=745)
READ(745,'(A)') cHeader

!!! Estimate the number of columns in csv
CALL getSubstring(cHeader,',',nColumns,Columns)
ALLOCATE(dValues(nColumns))

iElements = 0
DO
   READ(745,*,IOSTAT=io)  dValues
   
   IF (io > 0) THEN
      WRITE(*,*) 'Check input.  Something was wrong'
      EXIT
   ELSE IF (io < 0) THEN
      IF (myid.eq.1) WRITE(*,*) 'iParticles= ',iParticles,iElements
      EXIT
   ELSE
      dBuff = [dValues(Columns(2)),dValues(Columns(3)),dValues(Columns(4))]
      cdx = myParticleParam%dFacUnitSourcefile*dBuff(1)
      cdy = myParticleParam%dFacUnitSourcefile*dBuff(2)
      cdz = myParticleParam%dFacUnitSourcefile*dBuff(3)
      CALL SearchPointsWRT_STLs(cdx,cdy,cdz,0d0,dist_CGAL,bProjection,cdX0,cdY0,cdZ0,iProjection)
!      CALL SearchPointsWRT_STLs(cdx,cdy,cdz,0d0,dist_CGAL,.true.,cdX0,cdY0,cdZ0,iProjection)
!       write(*,'(8es12.4)') cdx,cdy,cdz,cdX0,cdY0,cdZ0, dist_CGAL,myParticleParam%d_CorrDist
      if (dist_CGAL.gt.1d0*myParticleParam%d_CorrDist) then
        iElements = iElements + 1
      end if
   END IF
END DO
CLOSE(745)

! CLOSE(785)

! write(*,*) myid,iElements

! We will create 4 times the number of particles
! from the sourcefile - because we will create
! the elements that are neccesary for the elongation
! (Lambda).
myParticleParam%nParticles = 4*iElements


ALLOCATE(myCompleteSet(1:myParticleParam%nParticles))
ALLOCATE(myActiveSet  (1:myParticleParam%nParticles))
ALLOCATE(myExchangeSet(1:myParticleParam%nParticles))
ALLOCATE(myLostSet    (1:myParticleParam%nParticles))

R_max = 0.5d0*myParticleParam%D_Out*myParticleParam%dEps2
R_min = 0.5d0*myParticleParam%D_in
X_box = 0.5d0*myParticleParam%D_Out*myParticleParam%dEps2
Y_box = 0.5d0*myParticleParam%D_Out*myParticleParam%dEps2
Z_Min = myParticleParam%Z_seed


ALLOCATE(Lambda(iElements),Length(myParticleParam%nParticles,3))
ALLOCATE(pPresent(myParticleParam%nParticles))
Lambda = -1d0

OPEN(FILE=myParticleParam%sourcefile,UNIT=745)
! Read the header and ignore it
READ(745,*)

! We read the file in one round
iElement = 0

DO
!DO iElement=1,iElements

  READ(745,*,IOSTAT=io)  dValues
  dBuff = [dValues(Columns(2)),dValues(Columns(3)),dValues(Columns(4))]
  myRandomNumber = 1d0*myParticleParam%dFacUnitSourcefile*dBuff

   IF (io > 0) THEN
      WRITE(*,*) 'Check input.  Something was wrong'
      EXIT
   ELSE IF (io < 0) THEN
      IF (myid.eq.1) WRITE(*,*) 'iParticles= ',iParticles,iElement
      EXIT
   ELSE

!       if (myid.eq.1) write(*,*) dBuff
      cdx = myParticleParam%dFacUnitSourcefile*dBuff(1)
      cdy = myParticleParam%dFacUnitSourcefile*dBuff(2)
      cdz = myParticleParam%dFacUnitSourcefile*dBuff(3)
      CALL SearchPointsWRT_STLs(cdx,cdy,cdz,0d0,dist_CGAL,bProjection,cdX0,cdY0,cdZ0,iProjection)
!      CALL GetDistToSTL(cdx,cdy,cdz,1,dist_CGAL,.true.)
      if (dist_CGAL.gt.1d0*myParticleParam%d_CorrDist) then
         
         iElement = iElement+ 1
         
         ! Calculate the indices of the particles that we
         ! are generating. We create them elementwise:
         ! First those from element 1, then those from
         ! element two, then those from element 3, ...
         idx1 = 4*iElement - 3
         idx2 = 4*iElement - 2
         idx3 = 4*iElement - 1
         idx4 = 4*iElement - 0

         myExchangeSet(idx1)%coor(:) = [myRandomNumber(1),myRandomNumber(2),myRandomNumber(3)]
         myExchangeSet(idx1)%velo    = 0d0
         myExchangeSet(idx1)%time    = 0d0
         myExchangeSet(idx1)%indice  = idx1
         if (dValues(Columns(1)).ne.0) myExchangeSet(idx1)%id = int(dValues(Columns(1))) 

         myExchangeSet(idx2)%coor(:) = [myRandomNumber(1)+myParticleParam%Epsilon,myRandomNumber(2),myRandomNumber(3)]
         myExchangeSet(idx2)%velo    = 0d0
         myExchangeSet(idx2)%time    = 0d0
         myExchangeSet(idx2)%indice  = idx2
         if (dValues(Columns(1)).ne.0) myExchangeSet(idx2)%id = int(dValues(Columns(1))) 

         myExchangeSet(idx3)%coor(:) = [myRandomNumber(1),myRandomNumber(2)+myParticleParam%Epsilon,myRandomNumber(3)]
         myExchangeSet(idx3)%velo    = 0d0
         myExchangeSet(idx3)%time    = 0d0
         myExchangeSet(idx3)%indice  = idx3
         if (dValues(Columns(1)).ne.0) myExchangeSet(idx3)%id = int(dValues(Columns(1))) 

         myExchangeSet(idx4)%coor(:) = [myRandomNumber(1),myRandomNumber(2),myRandomNumber(3)+myParticleParam%Epsilon]
         myExchangeSet(idx4)%velo    = 0d0
         myExchangeSet(idx4)%time    = 0d0
         myExchangeSet(idx4)%indice  = idx4
         if (dValues(Columns(1)).ne.0) myExchangeSet(idx4)%id = int(dValues(Columns(1))) 
     
      END IF
   END IF
END DO

write(*,*) iElement
CLOSE(745)

nExchangeSet = myParticleParam%nParticles

END SUBROUTINE Init_Particles_from_csv_xse
!
!-------------------------------------------------------------------
!
SUBROUTINE Init_Particles_from_csv(mfile)
INTEGER i,j,iParticles, iElements, iElement
REAL*8 myRandomNumber(3),R_min,R_max,Y,X_box,Y_box,Z_min
Character*256 cHeader
REAL*8 dBuff(3)
INTEGER nColumns,Columns(4)
REAL*8, allocatable :: dValues(:)
INTEGER iO

! Count how many particles are in the sourcefile!
! We allocate one element per particle that we find
iParticles = 0
OPEN(FILE=myParticleParam%sourcefile,UNIT=745)
READ(745,'(A)') cHeader

!!! Estimate the number of columns in csv
CALL getSubstring(cHeader,',',nColumns,Columns)
ALLOCATE(dValues(nColumns))

iElements = 0
DO
   READ(745,*,IOSTAT=io)  dValues
!   READ(745,*,IOSTAT=io)  dCrap,dBuff
   
   IF (io > 0) THEN
      WRITE(*,*) 'Check input.  Something was wrong'
      EXIT
   ELSE IF (io < 0) THEN
      IF (myid.eq.1) WRITE(*,*) 'iParticles= ',iParticles,iElements
      EXIT
   ELSE
      dBuff = [dValues(Columns(2)),dValues(Columns(3)),dValues(Columns(4))]
      cdx = 1d1*dBuff(1)
      cdy = 1d1*dBuff(2)
      cdz = 1d1*dBuff(3)
      CALL GetDistAndProjPToAllSTLs(cdx,cdy,cdz,cpx,cpy,cpz,dist_CGAL)
      
      if (dist_CGAL.gt.myParticleParam%d_CorrDist) then
!       CALL GetDistToSTL(cdx,cdy,cdz,1,dist_CGAL,.true.)
!       if ((     myParticleParam%bRotationalMovement.and.dist_CGAL.gt.myParticleParam%d_CorrDist).or.&
!           (.not.myParticleParam%bRotationalMovement.and.dist_CGAL.lt.myParticleParam%d_CorrDist)) then
        iElements = iElements + 1
      end if
   END IF
END DO
CLOSE(745)


! We will create 4 times the number of particles
! from the sourcefile - because we will create
! the elements that are neccesary for the elongation
! (Lambda).
myParticleParam%nParticles = 4*iElements


ALLOCATE(myCompleteSet(1:myParticleParam%nParticles))
ALLOCATE(myActiveSet  (1:myParticleParam%nParticles))
ALLOCATE(myExchangeSet(1:myParticleParam%nParticles))
ALLOCATE(myLostSet    (1:myParticleParam%nParticles))

R_max = 0.5d0*myParticleParam%D_Out*myParticleParam%dEps2
R_min = 0.5d0*myParticleParam%D_in
X_box = 0.5d0*myParticleParam%D_Out*myParticleParam%dEps2
Y_box = 0.5d0*myParticleParam%D_Out*myParticleParam%dEps2
Z_Min = myParticleParam%Z_seed


ALLOCATE(Lambda(iElements),Length(myParticleParam%nParticles,3))
ALLOCATE(pPresent(myParticleParam%nParticles))
Lambda = -1d0

OPEN(FILE=myParticleParam%sourcefile,UNIT=745)
! Read the header and ignore it
READ(745,*)

! We read the file in one round
iElement = 0
DO
!DO iElement=1,iElements

   READ(745,*,IOSTAT=io)  dValues
   dBuff = [dValues(Columns(2)),dValues(Columns(3)),dValues(Columns(4))]
   
!  READ(745,*,IOSTAT=io)  dCrap,dBuff
  myRandomNumber = 1d0*myParticleParam%dFacUnitSourcefile*dBuff

   IF (io > 0) THEN
      WRITE(*,*) 'Check input.  Something was wrong'
      EXIT
   ELSE IF (io < 0) THEN
      IF (myid.eq.1) WRITE(*,*) 'iParticles= ',iParticles,iElement
      EXIT
   ELSE

!       if (myid.eq.1) write(*,*) dBuff
      cdx = 1d1*dBuff(1)
      cdy = 1d1*dBuff(2)
      cdz = 1d1*dBuff(3)
      
      CALL GetDistAndProjPToAllSTLs(cdx,cdy,cdz,cpx,cpy,cpz,dist_CGAL)
      
      if (dist_CGAL.gt.myParticleParam%d_CorrDist) then
      
!       CALL GetDistToSTL(cdx,cdy,cdz,1,dist_CGAL,.true.)
!       if ((     myParticleParam%bRotationalMovement.and.dist_CGAL.gt.myParticleParam%d_CorrDist).or.&
!           (.not.myParticleParam%bRotationalMovement.and.dist_CGAL.lt.myParticleParam%d_CorrDist)) then
         
         iElement = iElement+ 1
         
         ! Calculate the indices of the particles that we
         ! are generating. We create them elementwise:
         ! First those from element 1, then those from
         ! element two, then those from element 3, ...
         idx1 = 4*iElement - 3
         idx2 = 4*iElement - 2
         idx3 = 4*iElement - 1
         idx4 = 4*iElement - 0

         myExchangeSet(idx1)%coor(:) = [myRandomNumber(1),myRandomNumber(2),myRandomNumber(3)]
         myExchangeSet(idx1)%time    = 0d0
         myExchangeSet(idx1)%indice  = idx1
         if (dValues(Columns(1)).ne.0) myExchangeSet(idx1)%id = int(dValues(Columns(1))) 

         myExchangeSet(idx2)%coor(:) = [myRandomNumber(1)+myParticleParam%Epsilon,myRandomNumber(2),myRandomNumber(3)]
         myExchangeSet(idx2)%time    = 0d0
         myExchangeSet(idx2)%indice  = idx2
         if (dValues(Columns(1)).ne.0) myExchangeSet(idx2)%id = int(dValues(Columns(1))) 

         myExchangeSet(idx3)%coor(:) = [myRandomNumber(1),myRandomNumber(2)+myParticleParam%Epsilon,myRandomNumber(3)]
         myExchangeSet(idx3)%time    = 0d0
         myExchangeSet(idx3)%indice  = idx3
         if (dValues(Columns(1)).ne.0) myExchangeSet(idx3)%id = int(dValues(Columns(1))) 

         myExchangeSet(idx4)%coor(:) = [myRandomNumber(1),myRandomNumber(2),myRandomNumber(3)+myParticleParam%Epsilon]
         myExchangeSet(idx4)%time    = 0d0
         myExchangeSet(idx4)%indice  = idx4
         if (dValues(Columns(1)).ne.0) myExchangeSet(idx4)%id = int(dValues(Columns(1))) 
     
      END IF
   END IF
END DO

write(*,*) iElement
!pause
CLOSE(745)

nExchangeSet = myParticleParam%nParticles

END SUBROUTINE Init_Particles_from_csv
!-------------------------------------------------------------------

! Assumtion:
! This is the output from an old run. This also means
! that the particles are in correct order for calculating the
! element elongation.
! This routine uses the index the particles have as ordering
! mechanism, so the file itself might look wild.
SUBROUTINE Init_Particles_from_old_output(mfile)
INTEGER i,j,iParticles
REAL*8 myRandomNumber(3),R_min,R_max,Y,X_box,Y_box,Z_min
REAL*8 dBuff(3)
INTEGER idx, idxmax,iID
INTEGER iO

! Count how many particles are in the sourcefile!
! Attention: his does not mean standard counting but
! intelligent counting: If you count 31500 particles in
! the file but found one with index 32000, then it is really
! likely that already some particles were lost - but
! to keep the old indicies, we still need to allocate
! arrays for 32000 particles.
! We also need to cope with the case that we found
! 31999 particles in the file (and the  um index is
! 31999 - the calculation of the elongation will fail
! if we allocate arrays for only 31999 particles, we still
! need to allocate arrays for 32000 particles.)
iParticles = 0
OPEN(FILE=myParticleParam%sourcefile,UNIT=745)
READ(745,*)
iParticles = 0
idxmax = 0
DO
   READ(745,*,IOSTAT=io)  dBuff, idx
   IF (io > 0) THEN
      WRITE(*,*) 'Check input.  Something was wrong'
      EXIT
   ELSE IF (io < 0) THEN
      !IF (myid.eq.1) WRITE(*,*) 'iParticles= ',iParticles
      EXIT
   ELSE
      iParticles = iParticles + 1
      idxmax = max(idxmax,idx)
   END IF
END DO
CLOSE(745)

myParticleParam%nParticles = max(iParticles, idxmax)
! Now the tricky part that is done quite easily:
! We know (see other init routines) that we have
! 4 particles per element (for the calculation
! of the element elongation). So we just add
! as many particles unitl it can be diveded by 4
do while (mod(myParticleParam%nParticles,4) .ne. 0)
  myParticleParam%nParticles = myParticleParam%nParticles + 1
end do
if (myid .eq. 1) write(*,'(A,I0,A)') 'initialised arrays for ', myParticleParam%nParticles, ' particles'

ALLOCATE(myCompleteSet(1:myParticleParam%nParticles))
ALLOCATE(myActiveSet  (1:myParticleParam%nParticles))
ALLOCATE(myExchangeSet(1:myParticleParam%nParticles))
ALLOCATE(myLostSet    (1:myParticleParam%nParticles))

R_max = 0.5d0*myParticleParam%D_Out*myParticleParam%dEps2
R_min = 0.5d0*myParticleParam%D_in
X_box = 0.5d0*myParticleParam%D_Out*myParticleParam%dEps2
Y_box = 0.5d0*myParticleParam%D_Out*myParticleParam%dEps2
Z_Min = myParticleParam%Z_seed

ALLOCATE(Lambda(int(myParticleParam%nParticles/4)),Length(myParticleParam%nParticles,3))
ALLOCATE(pPresent(myParticleParam%nParticles))

Lambda = -1d0

OPEN(FILE=myParticleParam%sourcefile,UNIT=745)
READ(745,*)

DO i=1,iParticles

  READ(745,*,IOSTAT=io)  dBuff, idx, iID
  myRandomNumber = 1d0*myParticleParam%dFacUnitSourcefile*dBuff

  myExchangeSet(idx)%coor(:) = [myRandomNumber(1),myRandomNumber(2),myRandomNumber(3)]
  myExchangeSet(idx)%time    = 0d0
  myExchangeSet(idx)%indice  = idx
  myExchangeSet(idx)%id      = iID

END DO

CLOSE(745)

nExchangeSet = myParticleParam%nParticles

END SUBROUTINE Init_Particles_from_old_output

!
! --------------------------------------------------------------------
!
SUBROUTINE OutputLostParticlesCSV()
CHARACTER*99 cFile
INTEGER i,j,k,nExSum
REAL*8 daux
INTEGER :: nxGrid,nyGrid
REAL*8 xmin,xmax,ymin,ymax,xp,yp, zp, zMinReach
REAL*8, ALLOCATABLE :: xGrid(:),yGrid(:)
REAL*8, ALLOCATABLE :: table(:,:,:)
REAL*8 :: dC1, dC2, dC3

IF (myid.eq.1) THEN

  ! Create the file and write the header
  cFile= '_RTD/ParticlesAtOutflow.csv'
  OPEN(FILE=TRIM(ADJUSTL(cFile)),UNIT = 412)

  WRITE(412,'(6A)') '"coor_X",','"coor_Y",','"coor_Z",', '"indice",', '"ID",','"time"'
  ! Now output the particles to the file
  DO i=1,nLostSet
   WRITE(412,'(3(E16.7,A),2(I0,A),E16.7,3(A,E16.7))') REAL(myLostSet(i)%coor(1)*myParticleParam%dFacUnitOut),',',&
                                REAL(myLostSet(i)%coor(2)*myParticleParam%dFacUnitOut),',',&
                                REAL(myLostSet(i)%coor(3)*myParticleParam%dFacUnitOut),',',&
                                myLostSet(i)%indice,',',myLostSet(i)%id,',',REAL(myLostSet(i)%time),',',&
                                REAL(myCompleteSet(i)%velo(1)),',',REAL(myCompleteSet(i)%velo(2)),',',REAL(myCompleteSet(i)%velo(3))
  END DO

  CLOSE(412)
  ! We will only accept particles in the raster that reached
  ! at least 99% of the z-length
  zMinReach = min(myParticleParam%OutflowZPos,myMeshInfo%zmin + (myMeshInfo%zmax - myMeshInfo%zmin)*0.99d0)

!  zMinReach = myParticleParam%OutflowZPos !myMeshInfo%zmin + (myMeshInfo%zmax - myMeshInfo%zmin)*0.99d0
!   zMinReach = myMeshInfo%zmin + (myMeshInfo%zmax - myMeshInfo%zmin)*0.99d0

  nxGrid = myParticleParam%Raster

  ! Prepare everything for the raster. Attention:
  ! This code breaks if the lost set is empty (which should not happen
  ! very often, but can happen if the maximum number of rotations is set
  ! too small).
  if (nLostSet .ge. 1 ) then
    xmin=myLostSet(1)%coor(1)
    xmax=myLostSet(1)%coor(1)
    ymin=myLostSet(1)%coor(2)
    ymax=myLostSet(1)%coor(2)

    DO i=2,nLostSet
     xmin = min(xmin, myLostSet(i)%coor(1))
     xmax = max(xmax, myLostSet(i)%coor(1))
     ymin = min(ymin, myLostSet(i)%coor(2))
     ymax = max(ymax, myLostSet(i)%coor(2))
    END DO

    nyGrid = INT(DBLE(nxGrid)*(ymax-ymin)/(xmax-xmin))
    ALLOCATE(xGrid(nxGrid+1))
    xGrid(1) = xmin
    DO i=1,nxGrid
     xGrid(i+1) = xGrid(i) + (xmax-xmin)/DBLE(nxGrid)
    END DO

    ALLOCATE(yGrid(nyGrid+1))
    yGrid(1) = ymin
    DO i=1,nyGrid
     yGrid(i+1) = yGrid(i) + (ymax-ymin)/DBLE(nyGrid)
    END DO

    ALLOCATE(table(nxGrid,nyGrid,2))
    table = 0d0

    DO i=1,nLostSet
     xp = myLostSet(i)%coor(1)
     yp = myLostSet(i)%coor(2)
     zp = myLostSet(i)%coor(3)
     if (zp .ge. zMinReach ) then
       DO j=1,nxGrid
        IF (xp.gt.xGrid(j).and.xp.lt.xGrid(j+1)) THEN
         DO  k=1,nyGrid
          IF (yp.gt.yGrid(k).and.yp.lt.yGrid(k+1)) THEN
           IF (myLostSet(i)%indice.gt.myParticleParam%nParticles/2) THEN
            table(j,k,2) = table(j,k,2) + 1d0
           ELSE
            table(j,k,1) = table(j,k,1) + 1d0
           END IF
          END IF
         END DO
        END IF
       END DO
     end if ! zp > zMinReach
    END DO

    OPEN(FILE='_RTD/Outflow.csv',UNIT = 412)
    WRITE(412,'(20A)') '"coor_X",','"coor_Y",','"coor_Z",', '"color1"',', ','"color2"', ',' ,'"color3"'

    DO j=1,nxGrid
     DO  k=1,nyGrid
      xp = 0.5d0*(xGrid(j) + xGrid(j+1))*myParticleParam%dFacUnitOut
      yp = 0.5d0*(yGrid(k) + yGrid(k+1))*myParticleParam%dFacUnitOut
      IF (table(j,k,1)+table(j,k,2).GT.0d0) THEN
       dC1 = table(j,k,1) / (table(j,k,1)+table(j,k,2))
       dC2 = table(j,k,2) / (table(j,k,1)+table(j,k,2))
       dC3 = 2.0d0*min(dC1,dC2)
       WRITE(412,'(6(E16.7,A))') xp,',',yp,',',0.0,',',dC1,',',dC2, ',', dC3
    !   ELSE
    !    dC1 = 0d0
    !    dC2 = 0d0
      END IF
     END DO
    END DO
    CLOSE(412)

  else
    ! This means we have less than 1 particle in the lost set.
    ! This means no particle is lost or has reached the outflow.
    ! we should not really do anything here, but since some
    ! postprocessing tools expect the file to be present
    ! in every case we will produce the file - but with
    ! complete garbage data
    OPEN(FILE='_RTD/Outflow.csv',UNIT = 412)
    WRITE(412,'(20A)') '"coor_X",','"coor_Y",','"coor_Z",', '"color1"',', ','"color2"', ',' ,'"color3"'
    write(412,'(6(E16.7,A))') 0.0,',', 0.0,',',0.0,',',-1.0,',',-1.0,',',-1.0
    close(412)
  end if ! (nLostSet .ge. 1)


END IF

END SUBROUTINE OutputLostParticlesCSV
!
! --------------------------------------------------------------------
!
SUBROUTINE OutputParticlesCSV(iT)
CHARACTER*99 cFile
INTEGER i,iT,nExSum
REAL*8 daux

daux = DBLE(nActiveSet)
CALL Comm_Summ(daux)
nExSum =INT(daux)
CALL Gather_Particle(nExSum)

IF (myid.eq.1) THEN

WRITE(cFile,'(A,I8.8,A4)') '_RTD/Particles_',iT,'.csv'
OPEN(FILE=TRIM(ADJUSTL(cFile)),UNIT = 412)

WRITE(412,'(9A)') '"coor_X",','"coor_Y",','"coor_Z",', '"indice",', '"ID",','"time",','"velo_X",','"velo_Y",','"velo_Z"'
DO i=1,nCompleteSet
 WRITE(412,'(3(E16.7,A),2(I0,A),E16.7,3(A,E16.7))') &
                               REAL(myCompleteSet(i)%coor(1)*myParticleParam%dFacUnitOut),',', &
                               REAL(myCompleteSet(i)%coor(2)*myParticleParam%dFacUnitOut),',', &
                               REAL(myCompleteSet(i)%coor(3)*myParticleParam%dFacUnitOut),',', &
                               myCompleteSet(i)%indice,',',myCompleteSet(i)%id,',',REAL(myCompleteSet(i)%time),',',&
                               REAL(myCompleteSet(i)%velo(1)),',',REAL(myCompleteSet(i)%velo(2)),',',REAL(myCompleteSet(i)%velo(3))
END DO

CLOSE(412)

END IF

END SUBROUTINE OutputParticlesCSV
!
! --------------------------------------------------------------------
!
SUBROUTINE GetLambda(lEndOfSimulation)
INTEGEr nH,i,j,nElements
REAL*8 X,Y,Z,ddd
CHARACTER*99 cFile
logical, intent(in) :: lEndOfSimulation

IF (myid.eq.1) THEN

  ! Only do something if we output a file!
  if (bOutputLambda .or. lEndOfSimulation) then
    pPresent = .FALSE.

    ! Place all coordinates sorted into the array
    ! Length - so that the coordinates of particle
    ! i end up in Length(i,:)
    DO i=1,nCompleteSet
     j = myCompleteSet(i)%indice
     X = myCompleteSet(i)%coor(1)
     Y = myCompleteSet(i)%coor(2)
     Z = myCompleteSet(i)%coor(3)
     Length(j,:) = [X,Y,Z]
     pPresent(j) = .TRUE.
    END DO

    ! Now calculate the elongations. Because of
    ! the ordering of the particles we first need
    ! to get the number of elements
    nElements = Size(Lambda)
    ! Calculate the elongation for one element:
    DO i=1,nElements
     ! Remember the ordering:
     ! One element after the other, and each element
     ! has 4 particles. So 4*i-3 is the first particle
     ! of this element
     i1 = 4*i - 3
     i2 = 4*i - 2
     i3 = 4*i - 1
     i4 = 4*i - 0
     IF (pPresent(i1).and.pPresent(i2).and.pPresent(i3).and.pPresent(i4)) THEN
      X = Length(i1,1)-Length(i2,1)
      Y = Length(i1,2)-Length(i3,2)
      Z = Length(i1,3)-Length(i4,3)
      ! Save the elongation for this element in the Lambda-Array.
      ! It is entry i because we did the calculation for the element i
      Lambda(i) = sqrt(X*X + Y*Y + Z*Z)
     END IF
    END DO

  end if ! (boutputLambda .or. lEndOfSimulation)


  IF (bOutputLambda) THEN
   ddd = (3d0*myParticleParam%Epsilon**2d0)**0.5d0
   WRITE(cFile,'(A,I3.3,A)') '_RTD/Lambda_Rot',iOutputLambda,'.csv'
   OPEN(412,FILE=TRIM(ADJUSTL(cFile)))
   DO i=1,Size(Lambda)
    WRITE(412,'(I10,4ES16.8)') i,Lambda(i),Lambda(i)/ddd,log10(Lambda(i)/ddd)
   END DO
   CLOSE(412)
  END IF

  IF (lEndOfSimulation) THEN
   ddd = (3d0*myParticleParam%Epsilon**2d0)**0.5d0
   WRITE(cFile,'(A)') '_RTD/Lambda_final.csv'
   OPEN(412,FILE=TRIM(ADJUSTL(cFile)))
   DO i=1,Size(Lambda)
    WRITE(412,'(I10,4ES16.8)') i,Lambda(i),Lambda(i)/ddd,log10(Lambda(i)/ddd)
   END DO
   CLOSE(412)
  END IF

END IF ! myId == 1

END SUBROUTINE GetLambda
!
! --------------------------------------------------------------------
!
SUBROUTINE GetCutplanes()
INTEGEr nH,i,j, iZpos
REAL*8 X,Y,Z,ddd,t,dID

IF (myid.eq.1) THEN
!  DO i=1,nCompleteSet
!
!   j = myCompleteSet(i)%indice
!   X = myCompleteSet(i)%coor(1)
!   Y = myCompleteSet(i)%coor(2)
!   Z = myCompleteSet(i)%coor(3)
!
!
!   IF (z.gt.myParticleParam%Z1.and.bPosParticles(1,j)) then
!    ZPosParticles(1,:,j) = [X,Y,myParticleParam%Z1]
!    bPosParticles(1,j)   = .FALSE.
!   END IF
!
!   IF (z.gt.myParticleParam%Z2.and.bPosParticles(2,j)) then
!    ZPosParticles(2,:,j) = [X,Y,myParticleParam%Z2]
!    bPosParticles(2,j)   = .FALSE.
!   END IF
!
!  END DO

  do i=1,nCompleteSet
    j = myCompleteSet(i)%indice
    X = myCompleteSet(i)%coor(1)
    Y = myCompleteSet(i)%coor(2)
    Z = myCompleteSet(i)%coor(3)
    t = myCompleteSet(i)%time
    dID = DBLE(myCompleteSet(i)%ID)
    do iZpos = 1,myParticleParam%nZposCutplanes
      IF (z.gt.myParticleParam%cutplanePositions(iZpos).and.bPosParticles(iZpos,j)) then
       ZPosParticles(iZpos,1:3,j) = [X,Y,myParticleParam%cutplanePositions(iZpos)]
       ZPosParticles(iZpos,4,j)   = t
       ZPosParticles(iZpos,5,j)   = dID
       bPosParticles(iZpos,j)     = .FALSE.
      END IF
    end do ! loop over all z positions
  end do ! over the complete set

!   do i=1,nCompleteSet
!     j = myCompleteSet(i)%indice
!     X = myCompleteSet(i)%coor(1)
!     Y = myCompleteSet(i)%coor(2)
!     Z = myCompleteSet(i)%coor(3)
!     do iZpos = 1,myParticleParam%nZposCutplanes
!       IF (z.gt.myParticleParam%cutplanePositions(iZpos).and.bPosParticles(iZpos,j)) then
!        ZPosParticles(iZpos,:,j) = [X,Y,myParticleParam%cutplanePositions(iZpos)]
!        bPosParticles(iZpos,j)   = .FALSE.
!       END IF
!     end do ! loop over all z positions
!   end do ! over the complete set
END IF

END SUBROUTINE GetCutplanes
!
! --------------------------------------------------------------------
!
SUBROUTINE OutputParticlesAtZtoCSV
CHARACTER*99 cFile
integer :: iZpos, j

IF (myid.eq.1) THEN

!  cFile= '_RTD/ParticlesAtZ1.csv'
!  OPEN(FILE=TRIM(ADJUSTL(cFile)),UNIT = 412)
!
!  WRITE(412,'(4A)') '"coor_X",','"coor_Y",','"coor_Z",', '"indice"'
!  DO j=1,myParticleParam%nParticles
! !   j = myCompleteSet(i)%indice
!   IF (.NOT.bPosParticles(1,j)) THEN
!    WRITE(412,'(3(E16.7,A),8I)') REAL(ZPosParticles(1,1,j)*myParticleParam%dFacUnitOut),',', &
!                                 REAL(ZPosParticles(1,2,j)*myParticleParam%dFacUnitOut),',', &
!                                 REAL(ZPosParticles(1,3,j)*myParticleParam%dFacUnitOut),',',j
!   END IF
!  END DO
!
!  CLOSE(412)
!
!  cFile= '_RTD/ParticlesAtZ2.csv'
!  OPEN(FILE=TRIM(ADJUSTL(cFile)),UNIT = 412)
!
!  WRITE(412,'(4A)') '"coor_X",','"coor_Y",','"coor_Z",', '"indice"'
!  DO j=1,myParticleParam%nParticles
! !   j = myCompleteSet(i)%indice
!   IF (.NOT.bPosParticles(2,j)) THEN
!    WRITE(412,'(3(E16.7,A),8I)') REAL(ZPosParticles(2,1,j)*myParticleParam%dFacUnitOut),',',&
!                                 REAL(ZPosParticles(2,2,j)*myParticleParam%dFacUnitOut),',',&
!                                 REAL(ZPosParticles(2,3,j)*myParticleParam%dFacUnitOut),',',j
!   END IF
!  END DO
!
!  CLOSE(412)

  do izPos=1,myParticleParam%nZposCutplanes

    write(cFile,'(A18,I2.2,A4)' ) '_RTD/ParticlesAtZ_' , izPos , '.csv'
    OPEN(FILE=TRIM(ADJUSTL(cFile)),UNIT = 412)

    WRITE(412,'(6A)') '"coor_X",','"coor_Y",','"coor_Z",', '"indice",', '"ID",','"time"'
    DO j=1,myParticleParam%nParticles

      IF (.NOT.bPosParticles(iZpos,j)) THEN
       WRITE(412,'(3(E16.7,A),2(I0,A),ES16.7)') &
                                     REAL(ZPosParticles(iZpos,1,j)*myParticleParam%dFacUnitOut),',', &
                                     REAL(ZPosParticles(iZpos,2,j)*myParticleParam%dFacUnitOut),',', &
                                     REAL(ZPosParticles(iZpos,3,j)*myParticleParam%dFacUnitOut),',',j,',',&
                                     INT(REAL(ZPosParticles(iZpos,5,j))),',', REAL(ZPosParticles(iZpos,4,j))
      END IF
    END DO
    CLOSE(412)
  end do ! loop over all Z positions

END IF

END SUBROUTINE OutputParticlesAtZtoCSV
!
! --------------------------------------------------------------------
!
SUBROUTINE OutPutParticles(iT)
CHARACTER*99 cFile
INTEGER i,iT,nExSum
REAL*8 daux

daux = DBLE(nActiveSet)
CALL Comm_Summ(daux)
nExSum =INT(daux)
CALL Gather_Particle(nExSum)

IF (myid.eq.1) THEN

WRITE(cFile,'(A,I4.4,A4)') 'Particles/out_',iT,'.vtu'

OPEN(FILE=TRIM(ADJUSTL(cFile)),UNIT = 412)

WRITE(412,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
WRITE(412,'(A)') '  <UnstructuredGrid>'
WRITE(412,'(A,I0,A)') '    <Piece NumberOfPoints="',nCompleteSet,'" NumberOfCells="0">'
WRITE(412,'(A)') '      <Points>'
WRITE(412,'(A)') '        <DataArray name="Position" type="Float32" NumberOfComponents="3" format="ascii">'
DO i=1,nCompleteSet
 WRITE(412,'(A10,3E16.7)') "          ",REAL(myCompleteSet(i)%coor(1)*myParticleParam%dFacUnitOut),&
                                        REAL(myCompleteSet(i)%coor(2)*myParticleParam%dFacUnitOut),&
                                        REAL(myCompleteSet(i)%coor(3)*myParticleParam%dFacUnitOut)
END DO
WRITE(412,'(A)') '        </DataArray>'
WRITE(412,'(A)') '      </Points>'
WRITE(412,'(A)') '      <PointData  Vectors="vector">'
! WRITE(412,'(A)') '        <DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="ascii">'
!         4 4 4 4 0 0 2 2 -2
! WRITE(412,'(A)') '        </DataArray>'
WRITE(412,'(A)') '    <DataArray type="Int32" Name="Indice" format="ascii">'
DO i=1,nCompleteSet
 WRITE(412,'(A10,I0)') "          ",myCompleteSet(i)%indice
END DO
WRITE(412,'(A)') '        </DataArray>'
WRITE(412,'(A)') '      </PointData>'
WRITE(412,'(A)') '      <Cells>'
WRITE(412,'(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
WRITE(412,'(A)') '        </DataArray>'
WRITE(412,'(A)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
WRITE(412,'(A)') '        </DataArray>'
WRITE(412,'(A)') '        <DataArray type="UInt8" Name="types" format="ascii">'
WRITE(412,'(A)') '        </DataArray>'
WRITE(412,'(A)') '      </Cells>'
WRITE(412,'(A)') '    </Piece>'
WRITE(412,'(A)') '  </UnstructuredGrid>'
WRITE(412,'(A)') '</VTKFile>'

CLOSE(412)
END IF

END SUBROUTINE OutPutParticles
!
! --------------------------------------------------------------------
!
SUBROUTINE PostProcessRTD()
INTEGER IOstatus,nData
REAL*8, ALLOCATABLE :: RTD_data(:,:),fineF(:)
LOGICAL BNEGATIVE


IF (myid.eq.1) THEN

 WRITE(*,*) 'filtering data ...'

 OPEN(FILE='_RTD/RTD.csv',UNIT=947)
!  READ(947,*)
 i = 0
 DO
  READ(947,*,IOSTAT=IOstatus) daux,daux,daux
  IF (IOSTATUS.NE.0) EXIT
  i = i + 1
 END DO
 REWIND(947)
 nData = i
 ALLOCATE(RTD_data(3,nData))
!  READ(947,*)
 jMark = 0
 DO i=1,nData
  READ(947,*,IOSTAT=IOstatus) RTD_data(:,i)
  IF (jMark.eq.0.and.RTD_data(3,i).NE.0) jMark = i-40
 END DO
 CLOSE (947)

 ALLOCATE(fineF(nData))
 fineF = 0d0

 jMark = 1
!  WRITE(*,*) 'going to subpostprocess ...'
 CALL subPostProcessRTD(RTD_data(1,jMark:),RTD_data(2,jMark:),fineF(jMark:),nData-jMark+1,20)

!  DO i=2,nData-1
! !  WRITE(*,*) (fineF(i)-fineF(i-1))*(fineF(i+1)-fineF(i))
!   IF ((fineF(i)-fineF(i-1))*(fineF(i+1)-fineF(i)).LT.0) THEN
!    iMarker=i
!    EXIT
!   END IF
!  END DO
!
!  WRITE(*,*) 'iMarker',iMarker
!
!  DO i=1,nData
!   IF (i.LT.iMarker.OR.fineF(i).LT.0d0) fineF(i) = 0d0
!  END DO

 WRITE(*,*) 'done....'




 OPEN(FILE='_RTD/RTD_filtered.csv',UNIT=947)

 WRITE(947,'(3(E12.4,A))') RTD_data(1,1),0e0,0e0
 DO i=2,nData
  WRITE(947,'(3(E12.4,A))') RTD_data(1,i),', ',fineF(i),', ',fineF(i)-fineF(i-1),' '
 END DO
 CLOSE (947)

END IF

END SUBROUTINE PostProcessRTD
!
! --------------------------------------------------------------------
!
SUBROUTINE subPostProcessRTD(oFineT,oFineF,nFineF,nF,nC)
REAL*8 oFineT(*),oFineF(*),nFineF(*)
INTEGER nF,nC
REAL*8, ALLOCATABLE :: crsTime(:)
REAL*8, ALLOCATABLE :: A(:,:),M(:,:),b(:),sol(:)
REAL*8 dFactor

ALLOCATE(crsTime(nC))
ALLOCATE(A(nF,nC+2),M(nC+2,nC+2),b(nC+2),sol(nC+2))

! WRITE(*,*) 'inside subpostprocess 1 ...', nF,nC

DO i=1,nF
 A(i,1) = 1d0
 A(i,2) = oFineT(i)
END DO

dFactor = oFineT(nF)/DBLE(nC)
DO j=1,nC
 crsTime(j) = dFactor*(j-1)
END DO

DO i=1,nF
 DO j=1,nC
  IF (oFineT(i)-crsTime(j).GT.0d0) THEN
   A(i,j+2) = (oFineT(i)-crsTime(j))**3d0
  ELSE
   A(i,j+2) = 0d0
  END IF
 END DO
END DO

! WRITE(*,*) 'inside subpostprocess 2 ...', nF,nC

DO i=1,nC+2
 DO j=1,nC+2

  M(i,j) = 0d0

  DO k=1,nF
   M(i,j) = M(i,j) + A(k,i)*A(k,j)
  END DO

 END DO
END DO

! WRITE(*,*) 'inside subpostprocess 3 ...', nF,nC

DO j=1,nC+2
 b(j) = 0d0
 DO i=1,nF
  b(j) = b(j) + oFineF(i)*A(i,j)
 END DO
END DO

! WRITE(*,*) 'writing b ...'
! WRITE(*,*) b

! WRITE(*,*) 'writing M ...'
! WRITE(*,*) M

nA = (nC+2)*(nC+2)

ALLOCATE (UMF_CMat(nA))
!UMF_CMat = M
ALLOCATE (UMF_lMat%ColA(nA))
ALLOCATE (UMF_lMat%LdA(nC+2+1))

UMF_lMat%LdA(1) = 1
jj = 0
DO j=1,nC+2

 UMF_lMat%LdA(j+1) =  UMF_lMat%LdA(j) + (nC+2)
 DO k=1,nC+2
  jj = jj + 1
  UMF_lMat%ColA(jj) = k
  UMF_CMat(jj) = M(j,k)
 END DO

END DO
UMF_lMat%nu   = nC+2
UMF_lMat%na   = nA

! WRITE(*,*) 'factorization ...', nF,nC
CALL myUmfPack_Factorize(UMF_CMat,UMF_lMat)
! WRITE(*,*) 'solution ...', nF,nC
CALL myUmfPack_Solve(sol,b,UMF_CMat,UMF_lMat,1)
! WRITE(*,*) 'writing sol ...'
! WRITE(*,*) sol

DO i=1,nF
 nFineF(i) = 0d0
 DO j=1,nC+2
  nFineF(i) = nFineF(i) + A(i,j)*sol(j)
 END DO
END DO

! WRITE(*,*) 'getting back ...'

END SUBROUTINE subPostProcessRTD

SUBROUTINE getSubstring(cH,cS,nC,iC)
use iniparser, only : inip_toupper_replace
Character*256 cH
Character*60 cF
Character*1 cS
integer iC(4),nC

integer i,jPos1,jPos2
integer, allocatable :: iPos(:)
Character*1 cA
Integer nS

nC = 0
iC = 0

if (myid.eq.1) write(*,*) ",",adjustl(trim(cH)),","

nS = 0
Do i=1,256
 READ(cH(i:i),'(A1)') cA
 if (cA.eq.cS) nS = nS + 1
END DO

allocate(iPos(nS+1))
iPos = 0

!write(*,*) 'nS=',nS

nS = 0
Do i=1,256
 READ(cH(i:i),'(A1)') cA
 if (cA.eq.cS) THEN
  nS = nS + 1
  iPos(nS) = i
 end if
END DO

iPos(nS+1) = LEN(adjustl(trim(cH)))+1
nC = nS + 1

if (myid.eq.1) write(*,*) 'nS=',nS, 'iPos= ',iPos

jPos1 = 0
Do i=1,nS+1
 jPos2 = iPos(i)
 READ(cH(jPos1+1:jPos2-1),'(A)') cF
 CALL inip_toupper_replace(cF)
 if (INDEX(ADJUSTL(TRIM(cF)),"RESULT").ne.0) then
  if (myid.eq.1) write(*,*) 'res-column is :', i,ADJUSTL(TRIM(cF))
   iC(1) = i
 end if
 if (INDEX(ADJUSTL(TRIM(cF)),"POINTS:0").ne.0) then
  if (myid.eq.1) write(*,*) 'x-column is :', i,ADJUSTL(TRIM(cF))
  iC(2) = i
 end if
 if (INDEX(ADJUSTL(TRIM(cF)),"POINTS:1").ne.0) then
  if (myid.eq.1) write(*,*) 'y-column is :', i,ADJUSTL(TRIM(cF))
  iC(3) = i
 end if
 if (INDEX(ADJUSTL(TRIM(cF)),"POINTS:2").ne.0) then
  iC(4) = i
  if (myid.eq.1) write(*,*) 'z-column is :', i,ADJUSTL(TRIM(cF))
 end if
 jPos1=jPos2
END DO

! pause

END SUBROUTINE getSubstring
END
!
!-------------------------------------------------------
!
SUBROUTINE SearchPointsWRT_STLs(X,Y,Z,t,dist,bProjection,X0,Y0,Z0,iProj)
use geometry_processing, only: STLR_elem,STLL_elem,STL_elem
use Sigma_User, only: mySigma
USE PP3D_MPI, ONLY : myid,master

IMPLICIT NONE
REAL*8 X,Y,Z,t,X0,Y0,Z0,dist,ProjP(3)
LOGICAL :: bProjection
INTEGER :: k,iProj
!!!!!
REAL*8  :: dSeg0,dSeg1,dSeg2,D0,D1,D2,tt
REAL*8 XC,YC,ZC
INTEGER :: inpr
REAL*8 :: ProjP0(3),ProjP1(3),ProjP2(3)

XC = X
YC = Y
ZC = Z

D0 = 1d8
D1 = 1d8
D2 = 1d8

!----------------------------------------------------------

!  if (bProjection) then
!   XC = 0.0d0
!   YC = 1.5d0
!   ZC = 5.0d0
!  end if

DO k=1, mySigma%NumberOfSeg
 dSeg0=1d8
 dSeg1=1d8
 dSeg2=1d8
 IF (mySigma%mySegment(k)%ART.EQ.'STL_R')  CALL STLR_elem(XC,YC,ZC,tt,k,dSeg1,dSeg2,inpr,bProjection,ProjP)
 IF (mySigma%mySegment(k)%ART.EQ.'STL_L')  CALL STLL_elem(XC,YC,ZC,tt,k,dSeg1,dSeg2,inpr,bProjection,ProjP)
 IF (mySigma%mySegment(k)%ART.EQ.'STL'  )  CALL STL_elem(XC,YC,ZC,tt,k,dSeg0,dSeg1,dSeg2,inpr,bProjection,ProjP)
 
!  if (myid.eq.1.and.bProjection) then
!   write(*,*) mySigma%mySegment(k)%ART,inpr,dSeg0
!   write(*,*) X,Y,Z
!   write(*,*) ProjP
!  end if

 if (dSeg1.lt.d1.and.bProjection) then
  ProjP1 = ProjP
  if (d1.lt.0d0) then
   XC = ProjP(1)
   YC = ProjP(2)
   ZC = ProjP(3)
  end if
 end if
 if (dSeg2.lt.d2.and.bProjection) then
  ProjP2 = ProjP
  if (d2.lt.0d0) then
   XC = ProjP(1)
   YC = ProjP(2)
   ZC = ProjP(3)
  end if
 end if
 if (dSeg0.lt.d0.and.bProjection) then
  ProjP0 = ProjP
  if (d0.lt.0d0) then
   XC = ProjP(1)
   YC = ProjP(2)
   ZC = ProjP(3)
  end if
 end if
 
 D0 = max(-1d8,min(dSeg0,D0))
 D1 = max(-1d8,min(dSeg1,D1))
 D2 = max(-1d8,min(dSeg2,D2))
  
END DO

dist = min(d0,d1,d2)

if (bProjection) then
 if ((d0.lt.d1).and.(d0.lt.d2)) then
  iProj = 1
  X0= ProjP0(1)
  Y0= ProjP0(2)
  Z0= ProjP0(3)
 else
  iProj = -1
  if (d1.lt.d2) then
   X0= ProjP1(1)
   Y0= ProjP1(2)
   Z0= ProjP1(3)
  else
   X0= ProjP2(1)
   Y0= ProjP2(2)
   Z0= ProjP2(3)
  end if
 end if
end if

END SUBROUTINE SearchPointsWRT_STLs
!
!-------------------------------------------------------
!
SUBROUTINE GetDistAndProjPToAllSTLs(X,Y,Z,X0,Y0,Z0,dist)
use fbm, only : myFBM

implicit none
REAL*8 X,Y,Z,X0,Y0,Z0,dist,X1,Y1,Z1
real*8 dist_sign,d_temp,daux
integer ip,ipc,isin,idynType

 Dist = 1d8
 
 DO IP = 1,myFBM%nParticles
  
  ipc=ip-1
  isin = 0
  call get_dynamics_type(ipc, idynType) 
  
  dist_sign = +1d0
  if (idynType.eq.2) dist_sign = -1d0
  
  call isinelementid(X,Y,Z,ipc,isin)
  if(isin .gt. 0)then
   dist_sign = -1d0*dist_sign
  else
   dist_sign = +1d0*dist_sign
  end if
  call getclosestpointid(X,Y,Z,x1,y1,z1,d_temp,ipc);

  if (dist_sign * d_temp.lt.dist) then
   dist = dist_sign * d_temp
   x0 = x1
   y0 = y1
   z0 = z1
  end if
  
 end do

END SUBROUTINE GetDistAndProjPToAllSTLs
