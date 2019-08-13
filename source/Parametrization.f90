MODULE Parametrization

USE PP3D_MPI, ONLY:myid,showid,myMPI_Barrier,master
USE var_QuadScalar

IMPLICIT NONE

TYPE tBndr
 INTEGER nVerts,nFaces
 LOGICAL, ALLOCATABLE :: Vert(:)
 INTEGER, ALLOCATABLE :: Face(:,:)
 REAL*8, ALLOCATABLE :: CoorList(:,:)
END TYPE tBndr

TYPE tParBndr
 TYPE(tBndr), ALLOCATABLE :: Bndr(:)
 CHARACTER :: Names*200,Types*200,Parameters*200
 INTEGER :: nBndrPar,CGAL_ID,Dimens
 REAL*8, ALLOCATABLE :: dBndrPar(:)
END TYPE tParBndr

INTEGER nBnds,iBnds
TYPE (tParBndr), ALLOCATABLE :: myParBndr(:)

CONTAINS
!----------------------------------------------------------------------------------------
SUBROUTINE InitBoundaryStructure(kvert,kedge)
implicit none
INTEGER kvert(8,*),kedge(12,*)
INTEGER ndof,i,j,k,ivt1,ivt2,iLoc
INTEGER iSym(3)
INTEGER Neigh(2,12)
DATA Neigh/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
CHARACTER cAux*50,cType*50
INTEGER   iType, iPhase, jPhase, iInterface,iTemp
logical :: bOK
 
 ndof = KNVT(ilev) + KNET(ilev) + KNAT(ilev) + KNEL(ilev)
 ALLOCATE(myBoundary%bWall(ndof))
 ALLOCATE(myBoundary%iInflow(ndof))
 ALLOCATE(myBoundary%iTemperature(ndof))
 ALLOCATE(myBoundary%iPhase(ndof))
 ALLOCATE(myBoundary%bOutflow(ndof))
 ALLOCATE(myBoundary%bSymmetry(3,ndof))
 ALLOCATE(BndrForce(ndof))
 ALLOCATE(myBoundary%LS_zero(ndof))
 ALLOCATE(myBoundary%bDisp_DBC(ndof)) 

 myBoundary%bWall     = .FALSE.
 myBoundary%iInflow   = 0
 myBoundary%iTemperature   = 0
 myBoundary%iPhase    = 0
 myBoundary%bOutflow  = .FALSE.
 myBoundary%bSymmetry = .FALSE.
 myBoundary%LS_zero   = 0
 BndrForce            = .FALSE.
 myBoundary%bDisp_DBC = .FALSE. 

 jPhase = 0
 DO iBnds=1,nBnds

  cAux = ADJUSTL(TRIM(myParBndr(iBnds)%Types))
  iType = 0
  IF (cAux(1:6).EQ.'Inflow') THEN
   READ(cAux(7:),'(I3)') iType
  END IF
  iTemp = 0
  IF (cAux(1:11).EQ.'Temperature') THEN
   READ(cAux(12:),'(I3)') iTemp
  END IF
  iPhase = 0
  IF (cAux(1:5).EQ.'Phase') THEN
   READ(cAux(6:),'(I2)') iPhase
  END IF
  iInterface = 0
  IF (cAux(1:6).EQ.'Bubble') THEN
   READ(cAux(7:),'(I2)') iInterface
  END IF  
  
  DO i=1,NVT
   IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Wall') myBoundary%bWall(i) = .TRUE.     
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Disp_DBC') myBoundary%bDisp_DBC(i) = .TRUE.
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'WallF') THEN
      myBoundary%bWall(i) = .TRUE.
      BndrForce(i) = .TRUE.
    END IF
    IF (iType.NE.0) myBoundary%iInflow(i) = iType
    IF (iTemp.GT.0) myBoundary%iTemperature(i) = iTemp
    IF (iPhase.GT.0) myBoundary%iPhase(i) = iPhase
    IF (iInterface.GT.0) myBoundary%LS_zero(i) = iInterface
!     IF (iPhase.GT.0) jPhase = jPhase + 1
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Outflow') myBoundary%bOutFlow(i) = .TRUE.
    iLoc = INDEX(myParBndr(iBnds)%Types,'Symmetry')
    IF (ILOC.EQ.1) THEN
     READ(myParBndr(iBnds)%Types(9:11),'(3I1)') iSym
     IF (iSym(1).EQ.1) myBoundary%bSymmetry(1,i) = .TRUE.
     IF (iSym(2).EQ.1) myBoundary%bSymmetry(2,i) = .TRUE.
     IF (iSym(3).EQ.1) myBoundary%bSymmetry(3,i) = .TRUE.
    END IF
   END IF
  END DO

  k=1
  DO i=1,NEL
   DO j=1,12
    IF (k.eq.kedge(j,i)) THEN
     ivt1 = kvert(Neigh(1,j),i)
     ivt2 = kvert(Neigh(2,j),i)
     IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(ivt1).AND.myParBndr(iBnds)%Bndr(ILEV)%Vert(ivt2)) THEN
      IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Wall') myBoundary%bWall(nvt+k) = .TRUE.
      IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Disp_DBC') myBoundary%bDisp_DBC(nvt+k) = .TRUE.      
      IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'WallF') THEN
        myBoundary%bWall(nvt+k) = .TRUE.
        BndrForce(nvt+k) = .TRUE.
      END IF
      IF (iType.NE.0) myBoundary%iInflow(nvt+k) = iType
      IF (iTemp.GT.0) myBoundary%iTemperature(nvt+k) = iTemp
      IF (iPhase.GT.0) myBoundary%iPhase(nvt+k) = iPhase
      IF (iInterface.GT.0) myBoundary%LS_zero(nvt+k) = iInterface
!       IF (iPhase.GT.0) jPhase = jPhase + 1
      IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Outflow') myBoundary%bOutFlow(nvt+k) = .TRUE.
      iLoc = INDEX(myParBndr(iBnds)%Types,'Symmetry')
      IF (ILOC.EQ.1) THEN
       READ(myParBndr(iBnds)%Types(9:11),'(3I1)') iSym
       IF (iSym(1).EQ.1) myBoundary%bSymmetry(1,nvt+k) = .TRUE.
       IF (iSym(2).EQ.1) myBoundary%bSymmetry(2,nvt+k) = .TRUE.
       IF (iSym(3).EQ.1) myBoundary%bSymmetry(3,nvt+k) = .TRUE.
      END IF
     END IF
     k = k + 1
    END IF
   END DO
  END DO

  DO i=1,NAT
   IF (myParBndr(iBnds)%Bndr(ILEV)%Face(1,i).NE.0) THEN
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Wall') myBoundary%bWall(nvt+net+i) = .TRUE.
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Disp_DBC') myBoundary%bDisp_DBC(nvt+net+i) = .TRUE.    
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'WallF') THEN
      myBoundary%bWall(nvt+net+i) = .TRUE.
      BndrForce(nvt+net+i) = .TRUE.
    END IF
    IF (iType.NE.0) myBoundary%iInflow(nvt+net+i) = iType
    IF (iTemp.GT.0) myBoundary%iTemperature(nvt+net+i) = iTemp
    IF (iPhase.GT.0) myBoundary%iPhase(nvt+net+i) = iPhase
    IF (iInterface.GT.0) myBoundary%LS_zero(nvt+net+i) = iInterface
!     IF (iPhase.GT.0) jPhase = jPhase + 1
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Outflow') myBoundary%bOutFlow(nvt+net+i) = .TRUE.
    iLoc = INDEX(myParBndr(iBnds)%Types,'Symmetry')
    IF (ILOC.EQ.1) THEN
     READ(myParBndr(iBnds)%Types(9:11),'(3I1)') iSym
     IF (iSym(1).EQ.1) myBoundary%bSymmetry(1,nvt+net+i) = .TRUE.
     IF (iSym(2).EQ.1) myBoundary%bSymmetry(2,nvt+net+i) = .TRUE.
     IF (iSym(3).EQ.1) myBoundary%bSymmetry(3,nvt+net+i) = .TRUE.
    END IF
   END IF   
  END DO

  DO i=1,NEL
   bOK = myParBndr(iBnds)%Bndr(ILEV)%Vert(kvert(1,i))
   DO j=2,8
     bOK = (bOK.AND.myParBndr(iBnds)%Bndr(ILEV)%Vert(kvert(j,i)))
   END DO
   IF (bOK) THEN
     IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Wall') myBoundary%bWall(nvt+net+nat+i) = .TRUE.
     IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Disp_DBC') myBoundary%bDisp_DBC(nvt+net+nat+i) = .TRUE.      
     IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'WallF') THEN
       myBoundary%bWall(nvt+net+nat+i) = .TRUE.
       BndrForce(nvt+net+nat+i) = .TRUE.
     END IF
     IF (iType.NE.0) myBoundary%iInflow(nvt+net+nat+i) = iType
     IF (iTemp.GT.0) myBoundary%iTemperature(nvt+net+nat+i) = iTemp
     IF (iPhase.GT.0) myBoundary%iPhase(nvt+net+nat+i) = iPhase
     IF (iPhase.GT.0) jPhase = jPhase + 1
     IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Outflow') myBoundary%bOutFlow(nvt+net+nat+i) = .TRUE.
     iLoc = INDEX(myParBndr(iBnds)%Types,'Symmetry')
     IF (ILOC.EQ.1) THEN
      READ(myParBndr(iBnds)%Types(9:11),'(3I1)') iSym
      IF (iSym(1).EQ.1) myBoundary%bSymmetry(1,nvt+net+nat+i) = .TRUE.
      IF (iSym(2).EQ.1) myBoundary%bSymmetry(2,nvt+net+nat+i) = .TRUE.
      IF (iSym(3).EQ.1) myBoundary%bSymmetry(3,nvt+net+nat+i) = .TRUE.
     END IF
   END IF
  END DO
 
 END DO

END SUBROUTINE InitBoundaryStructure
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE ParametrizeQ2Nodes(dCoor)
REAL*8 dCoor(3,*)
INTEGER I1,I2

 ILEV=NLMAX
 CALL SETLEV(2)
 i1 = NVT+1
 i2 = NVT+NET+NAT

 NLMAX=NLMAX+1
 ILEV=NLMAX
 CALL SETLEV(2)

 DO iBnds = 1, nBnds
  
   CALL Parametrize(dCoor,i1,i2,ilev)

 END DO

 NLMAX=NLMAX-1
 ILEV=NLMAX
 CALL SETLEV(2)

END SUBROUTINE ParametrizeQ2Nodes
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE ProlongateParametrization_STRCT(mgMesh,ilevel)
USE PP3D_MPI, ONLY: myid,coarse,myMPI_Barrier
IMPLICIT NONE
integer :: ilevel
type(tMultiMesh) :: mgMesh

INTEGER I1,I2

 DO iBnds = 1, nBnds

  IF (.NOT.ALLOCATED(myParBndr(iBnds)%Bndr(ilevel)%Vert)) THEN
   ALLOCATE (myParBndr(iBnds)%Bndr(ilevel)%Vert(mgMesh%level(ilevel)%NVT))
   ALLOCATE (myParBndr(iBnds)%Bndr(ilevel)%Face(2,mgMesh%level(ilevel)%NAT))

  myParBndr(iBnds)%Bndr(ilevel)%nVerts =0
  myParBndr(iBnds)%Bndr(ilevel)%nFaces =0
  myParBndr(iBnds)%Bndr(ilevel)%Vert =.FALSE.
  myParBndr(iBnds)%Bndr(ilevel)%Face = 0

  I1 = ilevel-1
  I2 = ilevel

  CALL GetHigherStructures(mgMesh%level(i1)%kvert,&
                           mgMesh%level(i1)%kedge,&
                           mgMesh%level(i2)%kvert,&
                           mgMesh%level(i2)%karea,&
                           mgMesh%level(i2)%kadj,&
                           mgMesh%level(i1)%nel,&
                           mgMesh%level(i1)%nvt,&
                           mgMesh%level(i1)%nat,&
                           mgMesh%level(i1)%net)
  END IF

 END DO

END SUBROUTINE ProlongateParametrization_STRCT
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE ParametrizeBndr(mgMesh,ilevel)
USE PP3D_MPI, ONLY: myid,coarse,myMPI_Barrier
IMPLICIT NONE
integer :: ilevel
type(tMultiMesh) :: mgMesh

INTEGER I1,I2

 DO iBnds = 1, nBnds

  IF (.NOT.ALLOCATED(myParBndr(iBnds)%Bndr(ilevel)%Vert)) THEN
   ALLOCATE (myParBndr(iBnds)%Bndr(ilevel)%Vert(mgMesh%level(ilevel)%NVT))
   ALLOCATE (myParBndr(iBnds)%Bndr(ilevel)%Face(2,mgMesh%level(ilevel)%NAT))

  myParBndr(iBnds)%Bndr(ilevel)%nVerts =0
  myParBndr(iBnds)%Bndr(ilevel)%nFaces =0
  myParBndr(iBnds)%Bndr(ilevel)%Vert =.FALSE.
  myParBndr(iBnds)%Bndr(ilevel)%Face = 0

  I1 = ilevel-1
  I2 = ilevel

  CALL GetHigherStructures(mgMesh%level(i1)%kvert,&
                           mgMesh%level(i1)%kedge,&
                           mgMesh%level(i2)%kvert,&
                           mgMesh%level(i2)%karea,&
                           mgMesh%level(i2)%kadj,&
                           mgMesh%level(i1)%nel,&
                           mgMesh%level(i1)%nvt,&
                           mgMesh%level(i1)%nat,&
                           mgMesh%level(i1)%net)
  END IF

  CALL Parametrize(mgmesh%level(ilevel)%dcorvg,1,mgmesh%level(ilevel)%nvt,ilevel)

 END DO

END SUBROUTINE ParametrizeBndr
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE ParametrizeBndryPoints(mgMesh,ilevel)
USE PP3D_MPI, ONLY: myid,coarse,myMPI_Barrier
IMPLICIT NONE
integer :: ilevel
type(tMultiMesh) :: mgMesh

 DO iBnds = 1, nBnds

  CALL Parametrize(mgmesh%level(ilevel)%dcorvg,1,mgmesh%level(ilevel)%nvt,ilevel)

 END DO

END SUBROUTINE ParametrizeBndryPoints
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE Parametrize(DCORVG,NVT1,NVT2,ilevel)
 implicit none
 REAL*8 DCORVG(3,*)
 INTEGER NVT1,NVT2
 integer :: ilevel

 INTEGER i,j
 REAL*8 dx,dy,dz,dist,dFact
 REAL*8 px,py,pz,dScale
 REAL*8 RX,RY,RZ,RAD,DFX,DFY,DFZ
 REAL*8 DA,DB,DC,DD,DSQUARE,DSUM,RAD1,RAD2,Z1,Z2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 IF (myParBndr(iBnds)%Types(1:6).EQ.'Bubble') THEN
  IF (itns.ge.1) THEN
   DO i=NVT1,NVT2
    IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i).and.myid.ne.0) THEN
     j = j + 1
     dx = DCORVG(1,i)
     dy = DCORVG(2,i)
     dz = DCORVG(3,i)
     CALL ParametrizePointOnTheBubble(dx,dy,dz,i)
     DCORVG(1,i) = dx
     DCORVG(2,i) = dy
     DCORVG(3,i) = dz
    END IF
   END DO
   RETURN
  END IF
 END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 IF (myParBndr(iBnds)%nBndrPar.EQ.0) RETURN

 IF (myParBndr(iBnds)%nBndrPar.EQ.1) THEN
  IF (.NOT.ALLOCATED(myParBndr(iBnds)%Bndr(ILEV)%CoorList)) THEN
   j = 0
   DO i=NVT1,NVT2
    IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
     j = j + 1
    END IF
   END DO
   ALLOCATE (myParBndr(iBnds)%Bndr(ILEV)%CoorList(3,j))
   j = 0
   DO i=NVT1,NVT2
    IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
     j = j + 1
     dx = DCORVG(1,i)
     dy = DCORVG(2,i)
     dz = DCORVG(3,i)
     myParBndr(iBnds)%Bndr(ILEV)%CoorList(:,j) = [dx,dy,dz]
    END IF
   END DO
  ELSE
   j = 0
   DO i=NVT1,NVT2
    IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
     j = j + 1
     px = myParBndr(iBnds)%Bndr(ILEV)%CoorList(1,j)
     py = myParBndr(iBnds)%Bndr(ILEV)%CoorList(2,j)
     pz = myParBndr(iBnds)%Bndr(ILEV)%CoorList(3,j)
     DCORVG(1,i) = px
     DCORVG(2,i) = py
     DCORVG(3,i) = pz
    END IF
   END DO
  END IF
 END IF

 IF (myParBndr(iBnds)%nBndrPar.EQ.4) THEN
  dA = myParBndr(iBnds)%dBndrPar(1)
  dB = myParBndr(iBnds)%dBndrPar(2)
  dC = myParBndr(iBnds)%dBndrPar(3)
  dD = myParBndr(iBnds)%dBndrPar(4)
  DSquare = dA*dA + dB*dB + dC*dC

  j = 0
  DO i=NVT1,NVT2
   IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
    j = j + 1
    dx = DCORVG(1,i)
    dy = DCORVG(2,i)
    dz = DCORVG(3,i)
!     dist = dA*dx + dB*dY + dC*dZ + dD
    dSum = (dx*dA + dy*dB + dz*dC + dD)
    px = dx - dA*dSum/dSquare
    py = dy - dB*dSum/dSquare
    pz = dz - dC*dSum/dSquare
!     WRITE(*,'(6e12.4)') dx,dy,dz,px,py,pz
!     WRITE(*,'(6e12.4)') dA,dB,dC,dD,dist
!     pause
    DCORVG(1,i) = px
    DCORVG(2,i) = py
    DCORVG(3,i) = pz
   END IF
  END DO
 END IF

 IF (myParBndr(iBnds)%nBndrPar.EQ.7) THEN

  RX  = myParBndr(iBnds)%dBndrPar(1)
  RY  = myParBndr(iBnds)%dBndrPar(2)
  RZ  = myParBndr(iBnds)%dBndrPar(3)
  RAD = myParBndr(iBnds)%dBndrPar(4)
  DFX = myParBndr(iBnds)%dBndrPar(5)
  DFY = myParBndr(iBnds)%dBndrPar(6)
  DFZ = myParBndr(iBnds)%dBndrPar(7)

  j = 0
  dScale = 1d0 !+ 1.0e-6
  IF (NVT1.ne.1) dScale = 1d0 !+ 0.2e-6
  DO i=NVT1,NVT2
   IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
    j = j + 1
    dx = DCORVG(1,i)-RX
    dy = DCORVG(2,i)-RY
    dz = DCORVG(3,i)-RZ
!     IF (nvt1.ne.1) write(*,*) myid,i,DCORVG(:,i),SIZE(myParBndr(iBnds)%Bndr(ILEV)%Vert)
    dist = DSQRT(DFX*dx**2d0 + DFY*dy**2d0 + DFZ*dz**2d0)
    dFact = dScale*RAD/dist
    px = (1d0-DFX)*DCORVG(1,i) + DFX*(RX + dFact*dx)
    py = (1d0-DFY)*DCORVG(2,i) + DFY*(RY + dFact*dy)
    pz = (1d0-DFZ)*DCORVG(3,i) + DFZ*(RZ + dFact*dz)
    DCORVG(1,i) = px
    DCORVG(2,i) = py
    DCORVG(3,i) = pz
   END IF
  END DO
 END IF
 
 ! Parametrization with respect to a varying radius cylinder
 IF (myParBndr(iBnds)%nBndrPar.EQ.10) THEN

  RX  = myParBndr(iBnds)%dBndrPar(1)
  RY  = myParBndr(iBnds)%dBndrPar(2)
  RZ  = myParBndr(iBnds)%dBndrPar(3)
  RAD1 = myParBndr(iBnds)%dBndrPar(4)
  Z1 = myParBndr(iBnds)%dBndrPar(5)
  RAD2 = myParBndr(iBnds)%dBndrPar(6)
  Z2 = myParBndr(iBnds)%dBndrPar(7)
  DFX = myParBndr(iBnds)%dBndrPar(8)
  DFY = myParBndr(iBnds)%dBndrPar(9)
  DFZ = myParBndr(iBnds)%dBndrPar(10)

  j = 0
  DO i=NVT1,NVT2
   IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
    j = j + 1
    dx = DCORVG(1,i)-RX
    dy = DCORVG(2,i)-RY
    dz = DCORVG(3,i)-RZ
    dist = DSQRT(DFX*dx**2d0 + DFY*dy**2d0 + DFZ*dz**2d0)
    RAD = RAD1 + (RAD2-RAD1)*(DCORVG(3,i)-Z1)/(Z2-Z1)
!     WRITE(*,*) dz
    dFact = RAD/dist
    px = (1d0-DFX)*DCORVG(1,i) + DFX*(RX + dFact*dx)
    py = (1d0-DFY)*DCORVG(2,i) + DFY*(RY + dFact*dy)
    pz = (1d0-DFZ)*DCORVG(3,i) + DFZ*(RZ + dFact*dz)
    DCORVG(1,i) = px
    DCORVG(2,i) = py
    DCORVG(3,i) = pz
   END IF
  END DO
 END IF


END SUBROUTINE Parametrize
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE Parametrize_old(DCORVG,NVT,ilevel)
 REAL*8 DCORVG(3,*)
 INTEGER NVT,NVT1
 integer :: ilevel
 INTEGER i,j
 REAL*8 dx,dy,dz,dist,dFact
 REAL*8 px,py,pz
 REAL*8 RX,RY,RZ,RAD,DFX,DFY,DFZ
 REAL*8 DA,DB,DC,DD,DSQUARE,DSUM,RAD1,RAD2,Z1,Z2

 IF (myParBndr(iBnds)%nBndrPar.EQ.0) RETURN

 IF (myParBndr(iBnds)%nBndrPar.EQ.1) THEN
  IF (.NOT.ALLOCATED(myParBndr(iBnds)%Bndr(ilevel)%CoorList)) THEN
   j = 0
   DO i=1,NVT
    IF (myParBndr(iBnds)%Bndr(ilevel)%Vert(i)) THEN
     j = j + 1
    END IF
   END DO
   ALLOCATE (myParBndr(iBnds)%Bndr(ilevel)%CoorList(3,j))
   j = 0
   DO i=1,NVT
    IF (myParBndr(iBnds)%Bndr(ilevel)%Vert(i)) THEN
     j = j + 1
     dx = DCORVG(1,i)
     dy = DCORVG(2,i)
     dz = DCORVG(3,i)
     myParBndr(iBnds)%Bndr(ilevel)%CoorList(:,j) = [dx,dy,dz]
    END IF
   END DO
  ELSE
   j = 0
   DO i=1,NVT
    IF (myParBndr(iBnds)%Bndr(ilevel)%Vert(i)) THEN
     j = j + 1
     px = myParBndr(iBnds)%Bndr(ilevel)%CoorList(1,j)
     py = myParBndr(iBnds)%Bndr(ilevel)%CoorList(2,j)
     pz = myParBndr(iBnds)%Bndr(ilevel)%CoorList(3,j)
     DCORVG(1,i) = px
     DCORVG(2,i) = py
     DCORVG(3,i) = pz
    END IF
   END DO
  END IF
 END IF

 IF (myParBndr(iBnds)%nBndrPar.EQ.4) THEN
  dA = myParBndr(iBnds)%dBndrPar(1)
  dB = myParBndr(iBnds)%dBndrPar(2)
  dC = myParBndr(iBnds)%dBndrPar(3)
  dD = myParBndr(iBnds)%dBndrPar(4)
  DSquare = dA*dA + dB*dB + dC*dC

  j = 0
  DO i=1,NVT
   IF (myParBndr(iBnds)%Bndr(ilevel)%Vert(i)) THEN
    j = j + 1
    dx = DCORVG(1,i)
    dy = DCORVG(2,i)
    dz = DCORVG(3,i)

    dSum = (dx*dA + dy*dB + dz*dC + dD)
    px = dx - dA*dSum/dSquare
    py = dy - dB*dSum/dSquare
    pz = dz - dC*dSum/dSquare
    
    DCORVG(1,i) = px
    DCORVG(2,i) = py
    DCORVG(3,i) = pz
   END IF
  END DO
 END IF

 IF (myParBndr(iBnds)%nBndrPar.EQ.7) THEN

 IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Bubble'.AND.TIMENS.LT.1d-5) THEN
  RETURN
 END IF
 
! IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Bubble'.AND.TIMENS.GT.1d-5) THEN
!  DO i=NVT1,NVT2
!   IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i).and.myid.ne.0) THEN
!    j = j + 1
!    dx = DCORVG(1,i)
!    dy = DCORVG(2,i)
!    dz = DCORVG(3,i)
!    CALL ParametrizePointOnTheBubble(dx,dy,dz,i)
!    DCORVG(1,i) = dx
!    DCORVG(2,i) = dy
!    DCORVG(3,i) = dz
!   END IF
!  END DO
!
!  RETURN
! END IF
 
  RX  = myParBndr(iBnds)%dBndrPar(1)
  RY  = myParBndr(iBnds)%dBndrPar(2)
  RZ  = myParBndr(iBnds)%dBndrPar(3)
  RAD = myParBndr(iBnds)%dBndrPar(4)
  DFX = myParBndr(iBnds)%dBndrPar(5)
  DFY = myParBndr(iBnds)%dBndrPar(6)
  DFZ = myParBndr(iBnds)%dBndrPar(7)

  j = 0
  DO i=1,NVT
   IF (myParBndr(iBnds)%Bndr(ilevel)%Vert(i)) THEN
    j = j + 1
    dx = DCORVG(1,i)-RX
    dy = DCORVG(2,i)-RY
    dz = DCORVG(3,i)-RZ
    dist = DSQRT(DFX*dx**2d0 + DFY*dy**2d0 + DFZ*dz**2d0)
    dFact = RAD/dist
    px = (1d0-DFX)*DCORVG(1,i) + DFX*(RX + dFact*dx)
    py = (1d0-DFY)*DCORVG(2,i) + DFY*(RY + dFact*dy)
    pz = (1d0-DFZ)*DCORVG(3,i) + DFZ*(RZ + dFact*dz)
    DCORVG(1,i) = px
    DCORVG(2,i) = py
    DCORVG(3,i) = pz
   END IF
  END DO
 END IF
 
 ! Parametrization with respect to a varying radius cylinder
 IF (myParBndr(iBnds)%nBndrPar.EQ.10) THEN

  RX  = myParBndr(iBnds)%dBndrPar(1)
  RY  = myParBndr(iBnds)%dBndrPar(2)
  RZ  = myParBndr(iBnds)%dBndrPar(3)
  RAD1 = myParBndr(iBnds)%dBndrPar(4)
  Z1 = myParBndr(iBnds)%dBndrPar(5)
  RAD2 = myParBndr(iBnds)%dBndrPar(6)
  Z2 = myParBndr(iBnds)%dBndrPar(7)
  DFX = myParBndr(iBnds)%dBndrPar(8)
  DFY = myParBndr(iBnds)%dBndrPar(9)
  DFZ = myParBndr(iBnds)%dBndrPar(10)

  j = 0
  DO i=1,NVT
   IF (myParBndr(iBnds)%Bndr(ilevel)%Vert(i)) THEN
    j = j + 1
    dx = DCORVG(1,i)-RX
    dy = DCORVG(2,i)-RY
    dz = DCORVG(3,i)-RZ
    dist = DSQRT(DFX*dx**2d0 + DFY*dy**2d0 + DFZ*dz**2d0)
    RAD = RAD1 + (RAD2-RAD1)*(DCORVG(3,i)-Z1)/(Z2-Z1)
!     WRITE(*,*) dz
    dFact = RAD/dist
    px = (1d0-DFX)*DCORVG(1,i) + DFX*(RX + dFact*dx)
    py = (1d0-DFY)*DCORVG(2,i) + DFY*(RY + dFact*dy)
    pz = (1d0-DFZ)*DCORVG(3,i) + DFZ*(RZ + dFact*dz)
    DCORVG(1,i) = px
    DCORVG(2,i) = py
    DCORVG(3,i) = pz
   END IF
  END DO
 END IF


END SUBROUTINE Parametrize_old
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE GetHigherStructures(KVERT1,KEDGE1,KVERT2,KAREA2,KADJ2,NEL1,NVT1,NAT1,NET1)
 INTEGER NAT1,NVT1,NEL1,NET1
 INTEGER KADJ2(6,*),KAREA2(6,*),KVERT2(8,*)
 INTEGER KEDGE1(12,*),KVERT1(8,*)
 INTEGER i,j,k,iel,iat,kat
 INTEGER kel(8)
 INTEGER jel(4),jat(4),jarea(4),jvt(4)
 INTEGER map(4,6),Neigh(2,12)
 DATA Neigh/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
 DATA map/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

 DO i=1,NAT1
  iel = myParBndr(iBnds)%Bndr(ILEV-1)%Face(1,i)
  iat = myParBndr(iBnds)%Bndr(ILEV-1)%Face(2,i)
  IF (iel.NE.0) THEN

   kel(1)=iel
   kel(2)=KADJ2(3,kel(1))
   kel(3)=KADJ2(3,kel(2))
   kel(4)=KADJ2(3,kel(3))
   kel(5)=KADJ2(6,kel(1))
   kel(6)=KADJ2(3,kel(5))
   kel(7)=KADJ2(3,kel(6))
   kel(8)=KADJ2(3,kel(7))

   IF (iat.EQ.1) THEN
    jel(1) = kel(1); jat(1) = 1
    jel(2) = kel(2); jat(2) = 1
    jel(3) = kel(3); jat(3) = 1
    jel(4) = kel(4); jat(4) = 1
   END IF

   IF (iat.EQ.2) THEN
    jel(1) = kel(1); jat(1) = 2
    jel(2) = kel(2); jat(2) = 5
    jel(3) = kel(6); jat(3) = 5
    jel(4) = kel(5); jat(4) = 2
   END IF

   IF (iat.EQ.3) THEN
    jel(1) = kel(2); jat(1) = 2
    jel(2) = kel(3); jat(2) = 5
    jel(3) = kel(7); jat(3) = 5
    jel(4) = kel(6); jat(4) = 2
   END IF

   IF (iat.EQ.4) THEN
    jel(1) = kel(3); jat(1) = 2
    jel(2) = kel(4); jat(2) = 5
    jel(3) = kel(8); jat(3) = 5
    jel(4) = kel(7); jat(4) = 2
   END IF

   IF (iat.EQ.5) THEN
    jel(1) = kel(4); jat(1) = 2
    jel(2) = kel(1); jat(2) = 5
    jel(3) = kel(5); jat(3) = 5
    jel(4) = kel(8); jat(4) = 2
   END IF

   IF (iat.EQ.6) THEN
    jel(1) = kel(5); jat(1) = 1
    jel(2) = kel(6); jat(2) = 1
    jel(3) = kel(7); jat(3) = 1
    jel(4) = kel(8); jat(4) = 1
   END IF

   jarea(1) = KAREA2(jat(1),jel(1))
   jarea(2) = KAREA2(jat(2),jel(2))
   jarea(3) = KAREA2(jat(3),jel(3))
   jarea(4) = KAREA2(jat(4),jel(4))

   myParBndr(iBnds)%Bndr(ILEV)%Face(:,jarea(1)) = [jel(1),jat(1)]
   myParBndr(iBnds)%Bndr(ILEV)%Face(:,jarea(2)) = [jel(2),jat(2)]
   myParBndr(iBnds)%Bndr(ILEV)%Face(:,jarea(3)) = [jel(3),jat(3)]
   myParBndr(iBnds)%Bndr(ILEV)%Face(:,jarea(4)) = [jel(4),jat(4)]

   myParBndr(iBnds)%Bndr(ILEV)%nFaces = myParBndr(iBnds)%Bndr(ILEV)%nFaces + 4

   DO kat = 1,4
    DO j=1,4
     jvt(j) = KVERT2(map(j,jat(kat)),jel(kat))
     IF (.NOT.myParBndr(iBnds)%Bndr(ILEV)%Vert(jvt(j))) THEN
      myParBndr(iBnds)%Bndr(ILEV)%nVerts = myParBndr(iBnds)%Bndr(ILEV)%nVerts + 1
      myParBndr(iBnds)%Bndr(ILEV)%Vert(jvt(j)) = .TRUE.
     END IF
    END DO
   END DO

  END IF
 END DO

 DO i=1,NVT1
  IF (myParBndr(iBnds)%Bndr(ILEV-1)%Vert(i)) myParBndr(iBnds)%Bndr(ILEV)%Vert(i) = .TRUE.
 END DO

 k=1
 DO i=1,NEL1
  DO j=1,12
   IF (k.eq.KEDGE1(j,i)) THEN
    jvt(1) = KVERT1(Neigh(1,j),i)
    jvt(2) = KVERT1(Neigh(2,j),i)
    IF (myParBndr(iBnds)%Bndr(ILEV-1)%Vert(jvt(1)).AND.myParBndr(iBnds)%Bndr(ILEV-1)%Vert(jvt(2))) THEN
     myParBndr(iBnds)%Bndr(ILEV)%Vert(NVT1+k) = .TRUE.
    END IF
    k = k + 1
   END IF
  END DO
 END DO

END SUBROUTINE GetHigherStructures
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE InitParametrization_STRCT(mesh,ilevel)
 type(tMesh) :: mesh
 integer :: ilevel
 INTEGER i,iVert,iLong,iAux,iError
 CHARACTER cFile*200,string*200
 integer :: iunit = 333
 integer :: istat

 CALL GetFileList()
 
 DO iBnds = 1, nBnds
  ALLOCATE (myParBndr(iBnds)%Bndr(1)%Vert(&
    mesh%NVT))
  ALLOCATE (myParBndr(iBnds)%Bndr(1)%Face(2,&
    mesh%NAT))
  cFile = ADJUSTL(TRIM(cProjectFolder))
  iLong = LEN(ADJUSTL(TRIM(cFile)))+1

  string = myParBndr(iBnds)%Names
  IF (myid.lt.1)     THEN 
   WRITE(cFile(iLong:),"(2A)") ADJUSTL(TRIM(string)),".par"
  ELSE
   WRITE(cFile(iLong:),"(A,A1,A3,A)") ADJUSTL(TRIM(string)),"_",cProjectNumber,".par"
  END IF

  OPEN(UNIT = iunit, FILE = TRIM(ADJUSTL(cFile)), action='read',iostat=istat)
  if(istat .ne. 0)then
    write(*,*)"Could not open file for reading. ",TRIM(ADJUSTL(cFile))
    stop          
  end if

  READ(iunit,*)  myParBndr(iBnds)%Bndr(ilevel)%nVerts, myParBndr(iBnds)%Types
  READ(iunit,*)  myParBndr(iBnds)%Parameters
  
  myParBndr(iBnds)%Bndr(1)%Vert = .FALSE.
  myParBndr(iBnds)%Bndr(1)%Face = 0

  DO i=1, myParBndr(iBnds)%Bndr(1)%nVerts
   READ(iunit,*) iVert
   myParBndr(iBnds)%Bndr(1)%Vert(iVert) = .TRUE.
  END DO

  CLOSE(iunit)

!  WRITE(*,'(I,3A)') myid,"|",myParBndr(iBnds)%Parameters,"|"
  READ(myParBndr(iBnds)%Parameters,*,IOSTAT=iError) myParBndr(iBnds)%nBndrPar
  IF (iError.NE.0.OR.myParBndr(iBnds)%nBndrPar.EQ.0) THEN
    myParBndr(iBnds)%nBndrPar = 0
  ELSE
    IF (myParBndr(iBnds)%nBndrPar.GT.0) THEN
      ALLOCATE(myParBndr(iBnds)%dBndrPar(myParBndr(iBnds)%nBndrPar))
      READ(myParBndr(iBnds)%Parameters,*) iAux,(myParBndr(iBnds)%dBndrPar(i),i=1,myParBndr(iBnds)%nBndrPar)
      myParBndr(iBnds)%Dimens = 3 ! ==> default:surface | will need to be redifened if we introduce specific line/point/volume parametrization ... 
    ELSE
      READ(myParBndr(iBnds)%Parameters,*) iAux,myParBndr(iBnds)%CGAL_ID
      myParBndr(iBnds)%Dimens = ABS(iAux)
      myParBndr(iBnds)%nBndrPar = 0
!       if (myid.eq.1) WRITE(*,*) myParBndr(iBnds)%Dimens,myParBndr(iBnds)%CGAL_ID
    END IF
  END IF
!  WRITE(*,'(2I8,<myParBndr(iBnds)%nBndrPar>D12.4)') myid,myParBndr(iBnds)%nBndrPar,myParBndr(iBnds)%dBndrPar

  CALL GetFaces(mesh%kvert,mesh%karea,mesh%NEL)
  CALL GetEdges()

 END DO

!   PAUSE

END SUBROUTINE InitParametrization_STRCT
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE DeterminePointParametrization_STRCT(mgMesh,ilevel)
 integer :: ilevel
 type(tMultiMesh) :: mgMesh
 integer inode,i,j,k,iel,iat,ivt1,ivt2,ivt3,ivt4,nn
 REAL*8, ALLOCATABLE :: daux(:)
 INTEGER NeighA(4,6),NeighU(4,6)
 DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
 DATA NeighU/1,2,3,4,1,6,9,5,2,7,10,6,3,8,11,7,4,5,12,8,9,10,11,12/

 ALLOCATE(mgMesh%BndryNodes(mgMesh%level(ilevel+1)%nvt))
 ALLOCATE(daux(mgMesh%level(ilevel+1)%nvt))
 daux = 0d0
 
 DO iel=1,mgMesh%level(ilevel)%nel
  DO j=1,6
   k = mgMesh%level(ilevel)%karea(j,iel)
   daux(mgMesh%level(ilevel)%nvt+mgMesh%level(ilevel)%net+k) = &
   daux(mgMesh%level(ilevel)%nvt+mgMesh%level(ilevel)%net+k) + 1d0
  END DO
 END DO
 
 ILEV = ilevel
 if (bParallel) CALL E013SUM(daux)

 k = 1
 DO iel=1,mgMesh%level(ilevel)%nel
   DO j=1,6
     IF (k.eq.mgMesh%level(ilevel)%karea(j,iel)) THEN
       IF (daux(mgMesh%level(ilevel)%nvt+mgMesh%level(ilevel)%net+k).eq.1d0) THEN
         ! CornerPoint
         ivt1 = mgMesh%level(ilevel)%kvert(NeighA(1,j),iel)
         ivt2 = mgMesh%level(ilevel)%kvert(NeighA(2,j),iel)
         ivt3 = mgMesh%level(ilevel)%kvert(NeighA(3,j),iel)
         ivt4 = mgMesh%level(ilevel)%kvert(NeighA(4,j),iel)
         mgMesh%BndryNodes(ivt1)%bOuterPoint = .TRUE.
         mgMesh%BndryNodes(ivt2)%bOuterPoint = .TRUE.
         mgMesh%BndryNodes(ivt3)%bOuterPoint = .TRUE.
         mgMesh%BndryNodes(ivt4)%bOuterPoint = .TRUE.
         ! EdgePoint
         ivt1 = mgMesh%level(ilevel)%kedge(NeighU(1,j),iel)
         ivt2 = mgMesh%level(ilevel)%kedge(NeighU(2,j),iel)
         ivt3 = mgMesh%level(ilevel)%kedge(NeighU(3,j),iel)
         ivt4 = mgMesh%level(ilevel)%kedge(NeighU(4,j),iel)
         mgMesh%BndryNodes(mgMesh%level(ilevel)%nvt+ivt1)%bOuterPoint = .TRUE.
         mgMesh%BndryNodes(mgMesh%level(ilevel)%nvt+ivt2)%bOuterPoint = .TRUE.
         mgMesh%BndryNodes(mgMesh%level(ilevel)%nvt+ivt3)%bOuterPoint = .TRUE.
         mgMesh%BndryNodes(mgMesh%level(ilevel)%nvt+ivt4)%bOuterPoint = .TRUE.
         ! FacePoint
         mgMesh%BndryNodes(mgMesh%level(ilevel)%nvt+mgMesh%level(ilevel)%net+k)%bOuterPoint = .TRUE.
       ELSE
!          WRITE(*,*)
       END IF
       k = k + 1
     END IF
   END DO
 END DO

 IF (myid.ne.MASTER) then
   DO iBnds = 1, nBnds
    DO iNode=1,mgMesh%level(ilevel+1)%nvt
      IF (myParBndr(iBnds)%Bndr(ilevel+1)%Vert(iNode)) THEN
        IF (myParBndr(iBnds)%Dimens.eq.1) mgMesh%BndryNodes(iNode)%nPoint   = mgMesh%BndryNodes(iNode)%nPoint + 1
        IF (myParBndr(iBnds)%Dimens.eq.2) mgMesh%BndryNodes(iNode)%nLine    = mgMesh%BndryNodes(iNode)%nLine + 1
        IF (myParBndr(iBnds)%Dimens.eq.3) mgMesh%BndryNodes(iNode)%nSurf    = mgMesh%BndryNodes(iNode)%nSurf + 1
        IF (myParBndr(iBnds)%Dimens.eq.4) mgMesh%BndryNodes(iNode)%nVolume  = mgMesh%BndryNodes(iNode)%nVolume + 1
      END IF
    END DO
   END DO

   DO iNode=1,mgMesh%level(ilevel+1)%nvt
     IF (mgMesh%BndryNodes(iNode)%nPoint .ne.0) ALLOCATE(mgMesh%BndryNodes(iNode)%P(mgMesh%BndryNodes(iNode)%nPoint))
     IF (mgMesh%BndryNodes(iNode)%nLine  .ne.0) ALLOCATE(mgMesh%BndryNodes(iNode)%L(mgMesh%BndryNodes(iNode)%nLine))
     IF (mgMesh%BndryNodes(iNode)%nSurf  .ne.0) ALLOCATE(mgMesh%BndryNodes(iNode)%S(mgMesh%BndryNodes(iNode)%nSurf))
     IF (mgMesh%BndryNodes(iNode)%nVolume.ne.0) ALLOCATE(mgMesh%BndryNodes(iNode)%V(mgMesh%BndryNodes(iNode)%nVolume))
   END DO

   DO iNode=1,mgMesh%level(ilevel+1)%nvt
     mgMesh%BndryNodes(iNode)%ParamTypes(1:4)=.FALSE.
     IF (mgMesh%BndryNodes(iNode)%nPoint .ne.0) mgMesh%BndryNodes(iNode)%ParamTypes(1)=.TRUE.
     IF (mgMesh%BndryNodes(iNode)%nLine  .ne.0) mgMesh%BndryNodes(iNode)%ParamTypes(2)=.TRUE.
     IF (mgMesh%BndryNodes(iNode)%nSurf  .ne.0) mgMesh%BndryNodes(iNode)%ParamTypes(3)=.TRUE.
     IF (mgMesh%BndryNodes(iNode)%nVolume.ne.0) mgMesh%BndryNodes(iNode)%ParamTypes(4)=.TRUE.
   END DO
   
   mgMesh%BndryNodes(:)%nPoint  = 0
   mgMesh%BndryNodes(:)%nLine   = 0
   mgMesh%BndryNodes(:)%nSurf   = 0
   mgMesh%BndryNodes(:)%nVolume = 0
   
   DO iBnds = 1, nBnds
     DO iNode=1,mgMesh%level(ilevel+1)%nvt
       IF (myParBndr(iBnds)%Bndr(ilevel+1)%Vert(iNode)) THEN
         IF (myParBndr(iBnds)%Dimens.eq.1) THEN
           mgMesh%BndryNodes(iNode)%nPoint   = mgMesh%BndryNodes(iNode)%nPoint + 1
           nn = mgMesh%BndryNodes(iNode)%nPoint
           mgMesh%BndryNodes(iNode)%P(nn) = iBnds
         END IF
         IF (myParBndr(iBnds)%Dimens.eq.2) THEN
           mgMesh%BndryNodes(iNode)%nLine    = mgMesh%BndryNodes(iNode)%nLine + 1
           nn = mgMesh%BndryNodes(iNode)%nLine
           mgMesh%BndryNodes(iNode)%L(nn) = iBnds
         END IF
         IF (myParBndr(iBnds)%Dimens.eq.3) THEN 
           mgMesh%BndryNodes(iNode)%nSurf    = mgMesh%BndryNodes(iNode)%nSurf + 1
           nn = mgMesh%BndryNodes(iNode)%nSurf
           mgMesh%BndryNodes(iNode)%S(nn) = iBnds
         END IF
         IF (myParBndr(iBnds)%Dimens.eq.4) THEN
           mgMesh%BndryNodes(iNode)%nVolume  = mgMesh%BndryNodes(iNode)%nVolume + 1
           nn = mgMesh%BndryNodes(iNode)%nVolume
           mgMesh%BndryNodes(iNode)%V(nn) = iBnds
         END IF
       END IF
     END DO
   END DO
  
   IF (bBoundaryCheck) THEN
    DO iNode=1,mgMesh%level(ilevel+1)%nvt
      IF (mgMesh%BndryNodes(iNode)%bOuterPoint) THEN
       IF (.NOT.(mgMesh%BndryNodes(iNode)%ParamTypes(1).OR.&
                 mgMesh%BndryNodes(iNode)%ParamTypes(2).OR.&
                 mgMesh%BndryNodes(iNode)%ParamTypes(3).OR.&
                 mgMesh%BndryNodes(iNode)%ParamTypes(4))) then
           WRITE(*,'(I4,A,I8,A)') myid,'BNDRY POINT:',iNode,' is not assigned to any parametrizations ...!'
       END IF
      END IF
    END DO
   END IF

 END IF
 
 if (bParallel) then
  daux = 0d0
  do i=1,mgMesh%level(ilevel+1)%nvt
   if (mgMesh%BndryNodes(i)%bOuterPoint) then
    daux(i) = 1d0
   end if
  end do
  CALL E013SUM(daux)
  do i=1,mgMesh%level(ilevel+1)%nvt
   if (daux(i).ge.1d0) then
    mgMesh%BndryNodes(i)%bOuterPoint = .true.
   end if
  end do
 end if
 
 deallocate(daux)
 
END SUBROUTINE DeterminePointParametrization_STRCT
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE ParametrizeBndryPoints_STRCT(mgMesh,ii)
USE PP3D_MPI, ONLY: myid,coarse,myMPI_Barrier
IMPLICIT NONE
integer :: ii
type(tMultiMesh) :: mgMesh
REAL*8 x,y,z,cpx,cpy,cpz,d_temp,xiS(3)
INTEGER iNode,ipc,iBnds,iS,W
integer :: iRepeat,nRepeat=8

 IF (myid.eq.MASTER) return
 
 DO iNode = 1, mgmesh%level(ii)%nvt
   IF ((mgMesh%BndryNodes(iNode)%ParamTypes(1).OR.&
        mgMesh%BndryNodes(iNode)%ParamTypes(2).OR.&
        mgMesh%BndryNodes(iNode)%ParamTypes(3).OR.&
        mgMesh%BndryNodes(iNode)%ParamTypes(4))) then
     IF (mgMesh%BndryNodes(iNode)%ParamTypes(1)) then
     
       GOTO 1
     END IF
     IF (mgMesh%BndryNodes(iNode)%ParamTypes(2)) then
        iBnds = mgMesh%BndryNodes(iNode)%L(1)
        ipc   = myParBndr(iBnds)%CGAL_ID
        x = dCGALtoRealFactor*mg_mesh%level(ii)%dcorvg(1,iNode)
        y = dCGALtoRealFactor*mg_mesh%level(ii)%dcorvg(2,iNode)
        z = dCGALtoRealFactor*mg_mesh%level(ii)%dcorvg(3,iNode)
        call projectonboundaryid(x,y,z,cpx,cpy,cpz,d_temp,ipc)
        mg_mesh%level(ii)%dcorvg(1,iNode) = cpx/dCGALtoRealFactor
        mg_mesh%level(ii)%dcorvg(2,iNode) = cpy/dCGALtoRealFactor
        mg_mesh%level(ii)%dcorvg(3,iNode) = cpz/dCGALtoRealFactor
       GOTO 1
     END IF
     IF (mgMesh%BndryNodes(iNode)%ParamTypes(3)) then
       xiS = 0d0
       W = DBLE(mgMesh%BndryNodes(iNode)%nSurf)
       DO iRepeat=1,nRepeat
       DO iS=1,mgMesh%BndryNodes(iNode)%nSurf
         iBnds = mgMesh%BndryNodes(iNode)%S(iS)
         IF (myParBndr(iBnds)%nBndrPar.NE.0) THEN !analytical
           x = mg_mesh%level(ii)%dcorvg(1,iNode)
           y = mg_mesh%level(ii)%dcorvg(2,iNode)
           z = mg_mesh%level(ii)%dcorvg(3,iNode)
           CALL projectonanalyticplane(x,y,z,cpx,cpy,cpz,iBnds)
           mg_mesh%level(ii)%dcorvg(1,iNode) = cpx
           mg_mesh%level(ii)%dcorvg(2,iNode) = cpy
           mg_mesh%level(ii)%dcorvg(3,iNode) = cpz
         ELSE  !CGAL
           ipc   = myParBndr(iBnds)%CGAL_ID
           x = dCGALtoRealFactor*mg_mesh%level(ii)%dcorvg(1,iNode)
           y = dCGALtoRealFactor*mg_mesh%level(ii)%dcorvg(2,iNode)
           z = dCGALtoRealFactor*mg_mesh%level(ii)%dcorvg(3,iNode)
           call projectonboundaryid(x,y,z,cpx,cpy,cpz,d_temp,ipc)
           mg_mesh%level(ii)%dcorvg(1,iNode) = cpx/dCGALtoRealFactor
           mg_mesh%level(ii)%dcorvg(2,iNode) = cpy/dCGALtoRealFactor
           mg_mesh%level(ii)%dcorvg(3,iNode) = cpz/dCGALtoRealFactor
         END IF
       ENDDO
       ENDDO
       GOTO 1
     END IF

   END IF
1  CONTINUE 
 END DO

END SUBROUTINE ParametrizeBndryPoints_STRCT
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE projectonanalyticplane(x1,y1,z1,x2,y2,z2,iBnds)
REAL*8 x1,y1,z1,x2,y2,z2
REAL*8 dx,dy,dz,dist,dFact
REAL*8 px,py,pz
REAL*8 RX,RY,RZ,RAD,DFX,DFY,DFZ
REAL*8 DA,DB,DC,DD,DSQUARE,DSUM,RAD1,RAD2,DZ1,DZ2
INTEGER iBnds,iSec,nSec
REAL*8, ALLOCATABLE :: Sec(:,:)

 IF (myParBndr(iBnds)%nBndrPar.EQ.0.OR.myParBndr(iBnds)%nBndrPar.EQ.1) THEN
  x2 = x1
  y2 = y1
  z2 = z1
  RETURN
 END IF

 IF (myParBndr(iBnds)%nBndrPar.EQ.4) THEN
  dA = myParBndr(iBnds)%dBndrPar(1)
  dB = myParBndr(iBnds)%dBndrPar(2)
  dC = myParBndr(iBnds)%dBndrPar(3)
  dD = myParBndr(iBnds)%dBndrPar(4)
!------------------------------------------------
  DSquare = dA*dA + dB*dB + dC*dC
  dx = x1
  dy = y1
  dz = z1
  dSum = (dx*dA + dy*dB + dz*dC + dD)
  px = dx - dA*dSum/dSquare
  py = dy - dB*dSum/dSquare
  pz = dz - dC*dSum/dSquare
  x2 = px
  y2 = py
  z2 = pz
  RETURN
 END IF

 IF (myParBndr(iBnds)%nBndrPar.EQ.7) THEN
 
  RX  = myParBndr(iBnds)%dBndrPar(1)
  RY  = myParBndr(iBnds)%dBndrPar(2)
  RZ  = myParBndr(iBnds)%dBndrPar(3)
  RAD = myParBndr(iBnds)%dBndrPar(4)
  DFX = myParBndr(iBnds)%dBndrPar(5)
  DFY = myParBndr(iBnds)%dBndrPar(6)
  DFZ = myParBndr(iBnds)%dBndrPar(7)
!------------------------------------------------
  dx = x1-RX
  dy = y1-RY
  dz = z1-RZ
  dist = DSQRT(DFX*dx**2d0 + DFY*dy**2d0 + DFZ*dz**2d0)
  dFact = RAD/dist
  px = (1d0-DFX)*x1 + DFX*(RX + dFact*dx)
  py = (1d0-DFY)*y1 + DFY*(RY + dFact*dy)
  pz = (1d0-DFZ)*z1 + DFZ*(RZ + dFact*dz)
  x2 = px
  y2 = py
  z2 = pz
  RETURN
 END IF
 
 ! Parametrization with respect to a varying radius cylinder
 IF (myParBndr(iBnds)%nBndrPar.GE.10) THEN
 
  nSec = (myParBndr(iBnds)%nBndrPar - 6)/2
  allocate(Sec(2,nSec))
  RX  = myParBndr(iBnds)%dBndrPar(1)
  RY  = myParBndr(iBnds)%dBndrPar(2)
  RZ  = myParBndr(iBnds)%dBndrPar(3)
  do iSec=1,nSec
   Sec(2,iSec) = myParBndr(iBnds)%dBndrPar(3+iSec)
  end do
  do iSec=1,nSec
   Sec(1,iSec) = myParBndr(iBnds)%dBndrPar(3+nSec+iSec)
  end do
  DFX = myParBndr(iBnds)%dBndrPar(3+2*(nSec-1)+3)
  DFY = myParBndr(iBnds)%dBndrPar(3+2*(nSec-1)+4)
  DFZ = myParBndr(iBnds)%dBndrPar(3+2*(nSec-1)+5)
  
!------------------------------------------------

  dx = x1-RX
  dy = y1-RY
  dz = z1-RZ
  do iSec=1,nSec-1
   RAD1= Sec(1,iSec  )
   RAD2= Sec(1,iSec+1)
   DZ1= Sec(2,iSec  )
   DZ2= Sec(2,iSec+1)
   if (dz.ge.DZ1.and.dz.le.DZ2) then
    dist = DSQRT(DFX*dx**2d0 + DFY*dy**2d0 + DFZ*dz**2d0)
    RAD = RAD1 + (RAD2-RAD1)*(z1-DZ1)/(DZ2-DZ1)
    dFact = RAD/dist
    px = (1d0-DFX)*x1 + DFX*(RX + dFact*dx)
    py = (1d0-DFY)*y1 + DFY*(RY + dFact*dy)
    pz = (1d0-DFZ)*z1 + DFZ*(RZ + dFact*dz)
   end if
  end do
  x2 = px
  y2 = py
  z2 = pz
  deallocate(Sec)
  RETURN
 END IF

END SUBROUTINE projectonanalyticplane
!
!----------------------------------------------------------------------------------------
!
 SUBROUTINE InitParametrization(mesh,ilevel)
 type(tMesh) :: mesh
 integer :: ilevel
 INTEGER i,iVert,iLong,iAux,iError
 CHARACTER cFile*200,string*200
 integer :: iunit = 333
 integer :: istat

 CALL GetFileList()

 DO iBnds = 1, nBnds
  ALLOCATE (myParBndr(iBnds)%Bndr(1)%Vert(&
    mesh%NVT))
  ALLOCATE (myParBndr(iBnds)%Bndr(1)%Face(2,&
    mesh%NAT))
  cFile = ADJUSTL(TRIM(cProjectFolder))
  iLong = LEN(ADJUSTL(TRIM(cFile)))+1

  string = myParBndr(iBnds)%Names
  IF (myid.lt.1)     THEN 
   WRITE(cFile(iLong:),"(2A)") ADJUSTL(TRIM(string)),".par"
  ELSE
   WRITE(cFile(iLong:),"(A,A1,A3,A)") ADJUSTL(TRIM(string)),"_",cProjectNumber,".par"
  END IF

  OPEN(UNIT = iunit, FILE = TRIM(ADJUSTL(cFile)), action='read',iostat=istat)
  if(istat .ne. 0)then
    write(*,*)"Could not open file for reading. ",TRIM(ADJUSTL(cFile))
    stop          
  end if

  READ(iunit,*)  myParBndr(iBnds)%Bndr(ilevel)%nVerts, myParBndr(iBnds)%Types
  READ(iunit,*)  myParBndr(iBnds)%Parameters

  myParBndr(iBnds)%Bndr(1)%Vert = .FALSE.
  myParBndr(iBnds)%Bndr(1)%Face = 0

  DO i=1, myParBndr(iBnds)%Bndr(1)%nVerts
   READ(iunit,*) iVert
   myParBndr(iBnds)%Bndr(1)%Vert(iVert) = .TRUE.
  END DO

  CLOSE(iunit)

!  WRITE(*,'(I,3A)') myid,"|",myParBndr(iBnds)%Parameters,"|"
  READ(myParBndr(iBnds)%Parameters,*,IOSTAT=iError) myParBndr(iBnds)%nBndrPar
  IF (iError.NE.0.OR.myParBndr(iBnds)%nBndrPar.EQ.0) THEN
    myParBndr(iBnds)%nBndrPar = 0
  ELSE
    IF (myParBndr(iBnds)%nBndrPar.GT.0) THEN
      ALLOCATE(myParBndr(iBnds)%dBndrPar(myParBndr(iBnds)%nBndrPar))
      READ(myParBndr(iBnds)%Parameters,*) iAux,(myParBndr(iBnds)%dBndrPar(i),i=1,myParBndr(iBnds)%nBndrPar)
      myParBndr(iBnds)%Dimens = 3 ! ==> default:surface | will need to be redifened if we introduce specific line/point/volume parametrization ... 
    ELSE
      READ(myParBndr(iBnds)%Parameters,*) iAux,myParBndr(iBnds)%CGAL_ID
      myParBndr(iBnds)%Dimens = ABS(iAux)
    END IF
  END IF
!  WRITE(*,'(2I8,<myParBndr(iBnds)%nBndrPar>D12.4)') myid,myParBndr(iBnds)%nBndrPar,myParBndr(iBnds)%dBndrPar
!  PAUSE

  CALL GetFaces(mesh%kvert,mesh%karea,mesh%NEL)
  CALL GetEdges()

  CALL Parametrize(mesh%dcorvg,1,mesh%NVT,ilevel)

 END DO


END SUBROUTINE InitParametrization
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE GetFileList()
 CHARACTER(LEN=200) :: string,cFile
 CHARACTER cWD*200,cSub*200
 INTEGER lenCommand,i,iPos,LenStr,iEnd

 if (myid.eq.1.or.(.not.bParallel))  Write(*,*)'Project File: "'//ADJUSTL(TRIM(cProjectFile))//'"'
 
 nBnds = 0
 OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cProjectFile)))
 DO
  READ(1,FMT='(200A)',IOSTAT=iEnd) string
  IF (iEnd.EQ.-1) EXIT
  LenStr = LEN(ADJUSTL(TRIM(string)))
  IF (LenStr.gt.4) THEN
   cFile = ADJUSTL(TRIM(string))
   IF (cFile(LenStr-3:LenStr).EQ.".par") nBnds = nBnds + 1
  END IF
 END DO

 if (myid.eq.1.or.(.not.bParallel)) Write(*,*)'Number of boundary parametrizations: ',nBnds
 ALLOCATE (myParBndr(nBnds))

 if (myid.eq.1.or.(.not.bParallel)) Write(*,*)'Warning: Allocated an additional boundary level'
 DO iBnds = 1, nBnds
  ALLOCATE (myParBndr(iBnds)%Bndr(NLMIN:NLMAX+1))
 END DO

 nBnds = 0
 REWIND(1)

 DO
  READ(1,FMT='(200A)',IOSTAT=iEnd) string
  IF (iEnd.EQ.-1) EXIT
  LenStr = LEN(ADJUSTL(TRIM(string)))
  IF (LenStr.gt.4) THEN
   cFile = ADJUSTL(TRIM(string))
   IF (cFile(LenStr-3:LenStr).EQ.".par") THEN
    if (myid.eq.1.or.(.not.bParallel)) Write(*,*)'Boundary parametrization: "'//ADJUSTL(TRIM(cFile))//'"'
    nBnds = nBnds + 1
    READ(cFile(1:LenStr-4),"(A)") myParBndr(nBnds)%Names
   END IF
  END IF
 END DO
 CLOSE(1)

END SUBROUTINE GetFileList
!----------------------------------------------------------------------------------------
SUBROUTINE GetEdges()

END SUBROUTINE GetEdges
!----------------------------------------------------------------------------------------
SUBROUTINE GetFaces(KVERT,KAREA,NEL)
 INTEGER KVERT(8,*),KAREA(6,*),NEL
 INTEGER i,iel,iat,iarea
 INTEGER IVT1,IVT2,IVT3,IVT4
 INTEGER map(4,6)
 DATA map/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

 myParBndr(iBnds)%Bndr(ILEV)%nFaces = 0

 DO iel=1,NEL
  DO iat=1,6
   
   ivt1 = KVERT(map(1,iat),iel)
   ivt2 = KVERT(map(2,iat),iel)
   ivt3 = KVERT(map(3,iat),iel)
   ivt4 = KVERT(map(4,iat),iel)

   IF (myParBndr(iBnds)%Bndr(1)%Vert(ivt1).AND.myParBndr(iBnds)%Bndr(1)%Vert(ivt2).AND.&
       myParBndr(iBnds)%Bndr(1)%Vert(ivt3).AND.myParBndr(iBnds)%Bndr(1)%Vert(ivt4)) THEN
       myParBndr(iBnds)%Bndr(1)%nFaces = myParBndr(iBnds)%Bndr(1)%nFaces+ 1
       iarea = KAREA(iat,iel)
       myParBndr(iBnds)%Bndr(1)%Face(1,iarea) = iel
       myParBndr(iBnds)%Bndr(1)%Face(2,iarea) = iat
   END IF

  END DO
 END DO

END SUBROUTINE GetFaces
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE ParametrizePointOnTheBubble(px,py,pz,iPoint)
REAL*8  :: px,py,pz

INTEGER i,j,k,kk,jClose,iTriang,iDef,iPoint
INTEGER nXX,iAux,iMinDist
REAL*8 cAux(3),dAux,P(3),ProjP(3),dTriang(3,3),dTestVector(4),AX2,AY2,AZ2,AX3,AY3,AZ3,ProjPdist,dist
REAL*8 tx,ty,tz
PARAMETER (nXX = 8)
REAL*8  :: dMonitor(nXX+1),cMonitor(3,nXX+1)
INTEGER :: iMonitor(nXX+1),jSet(nXX)
!!!!!!!!!!!!!!!!!!!!!!!!!!1
REAL*8 PP(3,9),dBAS(9),T(3),S(3),R(3),PPP(3,3), dMin_dL,dL,distP,distE
INTEGER lineFace(3,4)
LOGICAL bFound,bSFound
DATA LineFace/1,3,2, 3,5,4, 5,7,6, 7,1,8/
!!!1
INTEGER iIF

iIF = myBoundary%LS_zero(iPoint)

P = [px,py,pz]

dMonitor = 1d30
iMonitor = 0

DO i=1,myTSurf(iIF)%nT
 dist = (myTSurf(iIF)%T(i)%C(1,9)-P(1))**2d0 + (myTSurf(iIF)%T(i)%C(2,9)-P(2))**2d0 + (myTSurf(iIF)%T(i)%C(3,9)-P(3))**2d0
 IF (dist.lt.dMonitor(nXX)) THEN
  dMonitor(nXX) = dist
  iMonitor(nXX) = i
  DO kk=nXX-1,1,-1
   IF (dMonitor(kk+1).lt.dMonitor(kk)) THEN
    iAux = iMonitor(kk)
    dAux = dMonitor(kk)
    iMonitor(kk) = iMonitor(kk+1)
    dMonitor(kk) = dMonitor(kk+1)
    iMonitor(kk+1) = iAux
    dMonitor(kk+1) = dAux
   END IF
  END DO
 END IF
END DO

bFound = .FALSE.

DO i = 1,nXX

iMinDist = iMonitor(i)

PP(:,1) = myTSurf(iIF)%T(iMinDist)%C(:,1)
PP(:,2) = myTSurf(iIF)%T(iMinDist)%C(:,3)
PP(:,3) = myTSurf(iIF)%T(iMinDist)%C(:,5)
PP(:,4) = myTSurf(iIF)%T(iMinDist)%C(:,7)
PP(:,5) = myTSurf(iIF)%T(iMinDist)%C(:,2)
PP(:,6) = myTSurf(iIF)%T(iMinDist)%C(:,4)
PP(:,7) = myTSurf(iIF)%T(iMinDist)%C(:,6)
PP(:,8) = myTSurf(iIF)%T(iMinDist)%C(:,8)
PP(:,9) = myTSurf(iIF)%T(iMinDist)%C(:,9)

CALL ProjectPointOntoInterphaseSurf(P,T,PP)

IF (bFound) THEN
!      DIST = SQRT((px-t(1))**2d0+(py-t(2))**2d0+(pz-t(3))**2d0)
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) P
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) T
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) DIST
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) '----  A ----'
 EXIT
END IF

END DO

IF (.NOT.bFound) THEN

 dMin_dL = 1d30

 DO i = 1,nXX

  iMinDist = iMonitor(i)

  DO j=1,4
   PPP(:,1) = myTSurf(iIF)%T(iMinDist)%C(:,LineFace(1,j))
   PPP(:,2) = myTSurf(iIF)%T(iMinDist)%C(:,LineFace(2,j))
   PPP(:,3) = myTSurf(iIF)%T(iMinDist)%C(:,LineFace(3,j))
   bSFound = .FALSE.
   CALL ProjectPointOntoInterphaseLine(P,T,PPP,dL)	
   IF (bSFound) THEN 
    IF (dMin_dL.gt.dL) THEN
      bFound = .TRUE.
      dMin_dL = dL
      S = T
!      DIST = SQRT((px-t(1))**2d0+(py-t(2))**2d0+(pz-t(3))**2d0)
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) P
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) T
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) DIST
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) '----  E ----'
    END IF
   END IF
  END DO
 END DO

 IF (bFound) THEN
  T = S
  distE = dMin_dL
 END IF
END IF

! IF (.NOT.bFound) THEN
!  write(*,*) 'funny things are happenning .. '
!  dMin_dL = 1d30
 bSFound = .FALSE.
 DO i = 1,nXX
  iMinDist = iMonitor(i)
  DO j=1,9
   R = myTSurf(iIF)%T(iMinDist)%C(:,j)
   dL = (R(1)-P(1))**2d0+(R(2)-P(2))**2d0+(R(3)-P(3))**2d0
   IF (dMin_dL.gt.dL) THEN
     dMin_dL = dL
     S = R
     bSFound = .TRUE.
     distP = dMin_dL
!      DIST = SQRT((px-t(1))**2d0+(py-t(2))**2d0+(pz-t(3))**2d0)
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) P
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) T
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) DIST
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) '----  P ----'
   END IF
  END DO
 END DO

 IF (distP.lt.distE.and.bSFound) THEN
  T = S
 END IF
! END IF

tx = T(1)
ty = T(2)
tz = T(3)

DIST = SQRT((px-tx)**2d0+(py-ty)**2d0+(pz-tz)**2d0)

px = Tx
py = Ty
pz = Tz

! if (myid.eq.10.and.iPoint.eq.394) write(*,*) P
! if (myid.eq.10.and.iPoint.eq.394) write(*,*) T
! if (myid.eq.10.and.iPoint.eq.394) write(*,*) DIST
! if (myid.eq.10.and.iPoint.eq.394) write(*,*) '==== ===== ===== ====== ====== ====== ===== ===== '

 CONTAINS

 SUBROUTINE ProjectPointOntoInterphaseSurf(P,T,C)
 REAL*8 pP(2),P(3),T(3),pPx,pPy
 REAL*8 C(3,9),Der(2),Val(3),ddd
 REAL*8 XXX1,XXX2,XXX3,YYY1,YYY2,YYY3,df1,df2,dg1,dg2
 REAL*8 updateX,updateY
 INTEGER :: nW=10,iW
 REAL*8 dEPS,dBoxEPS
 
 dEPS = 1d-5
 dBoxEPS = 1d-5
 pP = 0d0
 
 CALL Get2DQ2basis(pP(1),pP(2),P,T,C,Val(1))
!  WRITE(*,'(3ES12.4)') pP,Val(1) 
 DO iW = 1,nW

  XXX1 = pP(1) - dEPS
  XXX2 = pP(1) 
  XXX3 = pP(1) + dEPS
  YYY1 = pP(2) - dEPS
  YYY2 = pP(2) 
  YYY3 = pP(2) + dEPS
  
!   CALL Get2DQ2basis(pPx,pP(2),P,T,C,Val(2))
  CALL Get2DQ2basis(XXX1,YYY2,P,T,C,Val(1))
  CALL Get2DQ2basis(XXX2,YYY2,P,T,C,Val(2))
  CALL Get2DQ2basis(XXX3,YYY2,P,T,C,Val(3))
  df1 = (Val(3)              - Val(1)) / (2d0*DEPS)
  df2 = (Val(3) - 2d0*Val(2) + Val(1)) / (DEPS*DEPS)
  CALL Get2DQ2basis(XXX2,YYY1,P,T,C,Val(1))
  CALL Get2DQ2basis(XXX2,YYY3,P,T,C,Val(3))
  dg1 = (Val(3)              - Val(1)) / (2d0*DEPS)
  dg2 = (Val(3) - 2d0*Val(2) + Val(1)) / (DEPS*DEPS)
  
  updateX = 1.0d0*df1/df2
  updateY = 1.0d0*dg1/dg2
  IF (df2.eq.0d0) updateX = 0d0
  IF (dg2.eq.0d0) updateY = 0d0

  ddd = SQRT(updateX**2d0 +updateY**2d0)
  
!   IF (myid.eq.1.and.iPoint.eq.32) WRITE(*,'(I,10ES12.4)') iW,Val(2),updateX,updateY,pP,dg1,dg2

  if (ddd.lt.1d-3) DEPS = 1d-6
  if (ddd.lt.1d-4) DEPS = 1d-7
  if (ddd.lt.1d-5) DEPS = 1d-8
  if (ddd.lt.1d-6) DEPS = 1d-9
    
  pP(1) =   MAX(-1d0-dBoxEPS,MIN(1d0+dBoxEPS,pP(1) - updateX))
  pP(2) =   MAX(-1d0-dBoxEPS,MIN(1d0+dBoxEPS,pP(2) - updateY))
  
  if (ddd.lt.1d-5) THEN
   bFound = .TRUE.
   EXIT
  END IF

 END DO

 IF (bFound) THEN
  XXX2 = pP(1) 
  YYY2 = pP(2) 
  CALL Get2DQ2basis(XXX2,YYY2,P,T,C,Val(2))
 END IF
 
END SUBROUTINE ProjectPointOntoInterphaseSurf

SUBROUTINE Get2DQ2basis(X1,X2,P,T,C,DIST)
 REAL*8 DHELP(9),X1,X2,P(3),C(3,9),T(3),DIST
 REAL*8 Q4,Q2
 PARAMETER (Q4=.25D0,Q2=.5D0)
 integer j

 DHELP(1)= Q4*(1D0-X1)*(1D0-X2)*X1*X2
 DHELP(2)=-Q4*(1D0+X1)*(1D0-X2)*X1*X2
 DHELP(3)= Q4*(1D0+X1)*(1D0+X2)*X1*X2
 DHELP(4)=-Q4*(1D0-X1)*(1D0+X2)*X1*X2
 DHELP(5)=-Q2*(1D0-X1*X1)*(1D0-X2)*X2
 DHELP(6)= Q2*(1D0+X1)*(1D0-X2*X2)*X1
 DHELP(7)= Q2*(1D0-X1*X1)*(1D0+X2)*X2
 DHELP(8)=-Q2*(1D0-X1)*(1D0-X2*X2)*X1
 DHELP(9)= (1D0-X1*X1)*(1D0-X2*X2)

 T = 0d0
 DO j=1,9
  T(1) = T(1) + DHELP(j)*C(1,j)
  T(2) = T(2) + DHELP(j)*C(2,j)
  T(3) = T(3) + DHELP(j)*C(3,j)
 END DO

 DIST = (T(1)-P(1))**2d0+(T(2)-P(2))**2d0+(T(3)-P(3))**2d0
 
END SUBROUTINE Get2DQ2basis

SUBROUTINE ProjectPointOntoInterphaseLine(P,T,C,dLineDIST)
 REAL*8 pP(1),P(3),T(3),pPx,pPy,dLineDIST
 REAL*8 C(3,3),Der(1),Val(3),ddd,daux
 REAL*8 XXX1,XXX2,XXX3,df1,df2
 REAL*8 updateX
 INTEGER :: nW=10,iW
 REAL*8 dEPS,dBoxEPS
 
 dEPS = 1d-5
 dBoxEPS = 1d-5
 
 pP = 0d0
 
 CALL Get1DQ2basis(pP(1),P,T,C,Val(1))

 DO iW = 1,nW

  XXX1 = pP(1) - dEPS
  XXX2 = pP(1) 
  XXX3 = pP(1) + dEPS
  
  CALL Get1DQ2basis(XXX1,P,T,C,Val(1))
  CALL Get1DQ2basis(XXX2,P,T,C,Val(2))
  CALL Get1DQ2basis(XXX3,P,T,C,Val(3))
  df1 = (Val(3)              - Val(1)) / (2d0*DEPS)
  df2 = (Val(3) - 2d0*Val(2) + Val(1)) / (DEPS*DEPS)
  
  updateX = 1.0d0*df1/df2
  IF (df2.eq.0d0) updateX = 0d0

  ddd = ABS(updateX)
  
!   IF (myid.eq.1.and.iPoint.eq.32) WRITE(*,'(I,10ES12.4)') iW,Val(2),updateX,updateY,pP,dg1,dg2

  if (ddd.lt.1d-3) DEPS = 1d-6
  if (ddd.lt.1d-4) DEPS = 1d-7
  if (ddd.lt.1d-5) DEPS = 1d-8
  if (ddd.lt.1d-6) DEPS = 1d-9
!   if (ddd.lt.1d-3) DEPS = 0.1*deps
    
  pP(1) =   MAX(-1d0-dBoxEPS,MIN(1d0+dBoxEPS,pP(1) - updateX))
  
  if (ddd.lt.1d-5) THEN
!    WRITE(*,*) pP,iW,i
   bSFound = .TRUE.
   EXIT
  END IF

 END DO
 
 IF (bSFound) THEN
  XXX2 = pP(1)
  CALL Get1DQ2basis(XXX2,P,T,C,dLineDIST)
 END IF

END SUBROUTINE ProjectPointOntoInterphaseLine

SUBROUTINE Get1DQ2basis(X1,P,T,C,DIST)
 REAL*8 DHELP(3),X1,P(3),C(3,3),T(3),DIST
 REAL*8 Q4,Q2
 PARAMETER (Q4=.25D0,Q2=.5D0)
 integer j

 DHELP(1)= -Q2*(1D0-X1)*X1
 DHELP(2)=  Q2*(1D0+X1)*X1
 DHELP(3)= (1D0-X1)*(1d0+X1)

 T = 0d0
 DO j=1,3
  T(1) = T(1) + DHELP(j)*C(1,j)
  T(2) = T(2) + DHELP(j)*C(2,j)
  T(3) = T(3) + DHELP(j)*C(3,j)
 END DO

 DIST = (T(1)-P(1))**2d0+(T(2)-P(2))**2d0+(T(3)-P(3))**2d0
 
END SUBROUTINE Get1DQ2basis

END SUBROUTINE ParametrizePointOnTheBubble

END MODULE Parametrization
