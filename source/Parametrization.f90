MODULE Parametrization

USE PP3D_MPI, ONLY:myid,showid
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
 CHARACTER :: Names*40,Types*40,Parameters*200
 INTEGER :: nBndrPar
 REAL*8, ALLOCATABLE :: dBndrPar(:)
END TYPE tParBndr

INTEGER nBnds,iBnds
TYPE (tParBndr), ALLOCATABLE :: myParBndr(:)

CONTAINS
!----------------------------------------------------------------------------------------
SUBROUTINE InitBoundaryStructure(kvert,kedge)
INTEGER kvert(8,*),kedge(12,*)
INTEGER ndof,i,j,k,ivt1,ivt2,iLoc
INTEGER iSym(3)
INTEGER Neigh(2,12)
DATA Neigh/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
CHARACTER cAux*50,cType*50
INTEGER   iType
 
 ndof = NVT + NET + NAT + NEL
 ALLOCATE(myBoundary%bWall(ndof))
 ALLOCATE(myBoundary%iInflow(ndof))
 ALLOCATE(myBoundary%bOutflow(ndof))
 ALLOCATE(myBoundary%bSymmetry(3,ndof))
 ALLOCATE(BndrForce(ndof))

 myBoundary%bWall     = .FALSE.
 myBoundary%iInflow   = 0
 myBoundary%bOutflow  = .FALSE.
 myBoundary%bSymmetry = .FALSE.
 BndrForce            = .FALSE.

 DO iBnds=1,nBnds

  cAux = ADJUSTL(TRIM(myParBndr(iBnds)%Types))
  iType = 0
  IF (cAux(1:6).EQ.'Inflow') THEN
   READ(cAux(7:),'(I)') iType
  END IF

  DO i=1,NVT
   IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Wall') myBoundary%bWall(i) = .TRUE.
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'WallF') THEN
      myBoundary%bWall(i) = .TRUE.
      BndrForce(i) = .TRUE.
    END IF
    IF (iType.GT.0) myBoundary%iInflow(i) = iType
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
      IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'WallF') THEN
        myBoundary%bWall(nvt+k) = .TRUE.
        BndrForce(nvt+k) = .TRUE.
      END IF
      IF (iType.GT.0) myBoundary%iInflow(nvt+k) = iType
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
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'WallF') THEN
      myBoundary%bWall(nvt+net+i) = .TRUE.
      BndrForce(nvt+net+i) = .TRUE.
    END IF
    IF (iType.GT.0) myBoundary%iInflow(nvt+net+i) = iType
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

 END DO

END SUBROUTINE InitBoundaryStructure
!----------------------------------------------------------------------------------------
SUBROUTINE ParametrizeBndr()
INTEGER I1,I2

 DO iBnds = 1, nBnds

  IF (.NOT.ALLOCATED(myParBndr(iBnds)%Bndr(ILEV)%Vert)) THEN
   ALLOCATE (myParBndr(iBnds)%Bndr(ILEV)%Vert(NVT))
   ALLOCATE (myParBndr(iBnds)%Bndr(ILEV)%Face(2,NAT))
  END IF

  myParBndr(iBnds)%Bndr(ILEV)%nVerts =0
  myParBndr(iBnds)%Bndr(ILEV)%nFaces =0
  myParBndr(iBnds)%Bndr(ILEV)%Vert =.FALSE.
  myParBndr(iBnds)%Bndr(ILEV)%Face = 0

! KVERT1,KEDGE1,KVERT2,KAREA2,KADJ2,NEL1,NVT1,NAT1
  I1 = ILEV-1
  I2 = ILEV
  CALL GetHigherStructures(KWORK(L(KLVERT(I1))),KWORK(L(KLEDGE(I1))),&
       KWORK(L(KLVERT(I2))),KWORK(L(KLAREA(I2))),KWORK(L(KLADJ(I2))),&
       KNEL(I1),KNVT(I1),KNAT(I1),KNET(I1))

  CALL Parametrize(DWORK(L(LCORVG)),NVT)

 END DO

END SUBROUTINE ParametrizeBndr
!----------------------------------------------------------------------------------------
SUBROUTINE Parametrize(DCORVG,NVT)
 REAL*8 DCORVG(3,*)
 INTEGER NVT,NVT1
 INTEGER i,j
 REAL*8 dx,dy,dz,dist,dFact
 REAL*8 px,py,pz
 REAL*8 RX,RY,RZ,RAD,DFX,DFY,DFZ
 REAL*8 DA,DB,DC,DD,DSQUARE,DSUM,RAD1,RAD2,Z1,Z2

 IF (myParBndr(iBnds)%nBndrPar.EQ.0) RETURN

 IF (myParBndr(iBnds)%nBndrPar.EQ.1) THEN
  IF (.NOT.ALLOCATED(myParBndr(iBnds)%Bndr(ILEV)%CoorList)) THEN
   j = 0
   DO i=1,NVT
    IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
     j = j + 1
    END IF
   END DO
   ALLOCATE (myParBndr(iBnds)%Bndr(ILEV)%CoorList(3,j))
   j = 0
   DO i=1,NVT
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
   DO i=1,NVT
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
  DO i=1,NVT
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
  DO i=1,NVT
   IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
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
!----------------------------------------------------------------------------------------
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
!----------------------------------------------------------------------------------------
SUBROUTINE InitParametrization()
 INTEGER i,iVert,iLong,iAux,iError
 CHARACTER cFile*100,string*10

 CALL GetFileList()

 DO iBnds = 1, nBnds
  ALLOCATE (myParBndr(iBnds)%Bndr(1)%Vert(NVT))
  ALLOCATE (myParBndr(iBnds)%Bndr(1)%Face(2,NAT))
  cFile = ADJUSTL(TRIM(cProjectFolder))
  iLong = LEN(ADJUSTL(TRIM(cFile)))+1

  string = myParBndr(iBnds)%Names
  IF (myid.lt.1)     THEN 
   WRITE(cFile(iLong:),"(2A)") ADJUSTL(TRIM(string)),".par"
  ELSE
   WRITE(cFile(iLong:),"(A,A1,A3,A)") ADJUSTL(TRIM(string)),"_",cProjectNumber,".par"
  END IF

!  WRITE(*,*) myid,nBnds,cFile
  OPEN(FILE = TRIM(ADJUSTL(cFile)),UNIT = 333)

  READ(333,*)  myParBndr(iBnds)%Bndr(ILEV)%nVerts, myParBndr(iBnds)%Types
  READ(333,*)  myParBndr(iBnds)%Parameters

  myParBndr(iBnds)%Bndr(1)%Vert = .FALSE.
  myParBndr(iBnds)%Bndr(1)%Face = 0

  DO i=1, myParBndr(iBnds)%Bndr(1)%nVerts
   READ(333,*) iVert
   myParBndr(iBnds)%Bndr(1)%Vert(iVert) = .TRUE.
  END DO

  CLOSE(333)

!  WRITE(*,'(I,3A)') myid,"|",myParBndr(iBnds)%Parameters,"|"
  READ(myParBndr(iBnds)%Parameters,*,IOSTAT=iError) myParBndr(iBnds)%nBndrPar
  IF (myParBndr(iBnds)%nBndrPar.GT.0.AND.iError.EQ.0) THEN
   ALLOCATE(myParBndr(iBnds)%dBndrPar(myParBndr(iBnds)%nBndrPar))
   READ(myParBndr(iBnds)%Parameters,*) iAux,(myParBndr(iBnds)%dBndrPar(i),i=1,myParBndr(iBnds)%nBndrPar)
  ELSE
   myParBndr(iBnds)%nBndrPar = 0
  END IF
!  WRITE(*,'(2I8,<myParBndr(iBnds)%nBndrPar>D12.4)') myid,myParBndr(iBnds)%nBndrPar,myParBndr(iBnds)%dBndrPar
!  PAUSE

  CALL GetFaces(KWORK(L(LVERT)),KWORK(L(LAREA)),NEL)
  CALL GetEdges()

  CALL Parametrize(DWORK(L(LCORVG)),NVT)

 END DO


END SUBROUTINE InitParametrization
!----------------------------------------------------------------------------------------
SUBROUTINE GetFileList()
 CHARACTER(LEN=40) :: string,cFile
 CHARACTER cWD*40,cSub*40
 INTEGER lenCommand,i,iPos,LenStr,iEnd

 nBnds = 0
 OPEN(UNIT=1,NAME=ADJUSTL(TRIM(cProjectFile)))
 DO
  READ(1,FMT='(40A)',IOSTAT=iEnd) string
  IF (iEnd.EQ.-1) EXIT
  LenStr = LEN(ADJUSTL(TRIM(string)))
  IF (LenStr.gt.4) THEN
   cFile = ADJUSTL(TRIM(string))
   IF (cFile(LenStr-3:LenStr).EQ.".par") nBnds = nBnds + 1
  END IF
 END DO

 ALLOCATE (myParBndr(nBnds))

 DO iBnds = 1, nBnds
  ALLOCATE (myParBndr(iBnds)%Bndr(NLMIN:NLMAX))
 END DO

 nBnds = 0
 REWIND(1)

 DO
  READ(1,FMT='(40A)',IOSTAT=iEnd) string
  IF (iEnd.EQ.-1) EXIT
  LenStr = LEN(ADJUSTL(TRIM(string)))
  IF (LenStr.gt.4) THEN
   cFile = ADJUSTL(TRIM(string))
   IF (cFile(LenStr-3:LenStr).EQ.".par") THEN
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

END MODULE Parametrization
