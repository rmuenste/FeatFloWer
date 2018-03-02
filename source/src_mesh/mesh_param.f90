module mesh_param

contains

SUBROUTINE mp_ParametrizeBndr()
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

END SUBROUTINE mp_ParametrizeBndr
!----------------------------------------------------------------------------------------
SUBROUTINE mp_Parametrize(DCORVG,NVT)
 REAL*8 DCORVG(3,*)
 INTEGER NVT,NVT1
 INTEGER i,j
 REAL*8 dx,dy,dz,dist,dFact
 REAL*8 px,py,pz
 REAL*8 RX,RY,RZ,RAD,DFX,DFY,DFZ
 REAL*8 DA,DB,DC,DD,DSQUARE,DSUM,RAD1,RAD2,Z1,Z2
 REAL*8 dFrac,dAlpha,dBeta,dxx,dyy
 REAL*8 :: myPI=dATAN(1d0)*4d0

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
!     DCORVG(3,i) = pz
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

 IF (myParBndr(iBnds)%nBndrPar.EQ.8) THEN
  RX  = myParBndr(iBnds)%dBndrPar(1)
  RY  = myParBndr(iBnds)%dBndrPar(2)
  RZ  = myParBndr(iBnds)%dBndrPar(3)
  RAD = myParBndr(iBnds)%dBndrPar(4)
  DFX = myParBndr(iBnds)%dBndrPar(6)
  DFY = myParBndr(iBnds)%dBndrPar(7)
  DFZ = myParBndr(iBnds)%dBndrPar(8)

  j = 0
  DO i=1,NVT
   IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
    j = j + 1
    dx = DCORVG(1,i)-RX
    dy = DCORVG(2,i)-RY

    dBeta=2d0*(DCORVG(3,i)-0d0)*myPI/54d0
    dyy = dx*cos(dBeta) - dy*sin(dBeta)
    dxx = dx*sin(dBeta) + dy*cos(dBeta)

    dist = DSQRT(DFX*dxx**2d0 + DFY*dyy**2d0)
    dFact = RAD/dist
    px = RX + dFact*dx
    py = RY + dFact*dy
    pz = 0d0
    DCORVG(1,i) = px
    DCORVG(2,i) = py
!    DCORVG(3,i) = pz
   END IF
  END DO
 END IF

 !  IF (myParBndr(iBnds)%nBndrPar.EQ.8) THEN
!   RX  = myParBndr(iBnds)%dBndrPar(1)
!   RY  = myParBndr(iBnds)%dBndrPar(2)
!   RZ  = myParBndr(iBnds)%dBndrPar(3)
!   RAD1 = myParBndr(iBnds)%dBndrPar(4)
!   RAD2 = myParBndr(iBnds)%dBndrPar(5)
!   DFX = myParBndr(iBnds)%dBndrPar(6)
!   DFY = myParBndr(iBnds)%dBndrPar(7)
!   DFZ = myParBndr(iBnds)%dBndrPar(8)
! 
!   j = 0
!   DO i=1,NVT
!    IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
!     j = j + 1
!     dx = DCORVG(1,i)-RX
!     dy = DCORVG(2,i)-RY
!     dz = DCORVG(3,i)-RZ
! !     dist = DSQRT(DFX*dx**2d0 + DFY*dy**2d0 + DFZ*dz**2d0)
! !     dAlpha = atan((DCORVG(2,i)-RY)/(DCORVG(1,i)-RX))
! !     dFrac  = (0.5d0*myPI + dAlpha)/myPI
! !     RAD = (1d0-dFrac)*RAD1 + (dFrac)*RAD2
! !     dFact = RAD/dist
! !     RAD = RAD1
! !     px = (1d0-DFX)*DCORVG(1,i) + DFX*(RX + dFact*dx)
! !     py = (1d0-DFY)*DCORVG(2,i) + DFY*(RY + dFact*dy)
! !     pz = (1d0-DFZ)*DCORVG(3,i) + DFZ*(RZ + dFact*dz)
!     dist = DSQRT(dx**2d0 + dy**2d0)
!     dAlpha = atan((DCORVG(2,i)-RY)/(DCORVG(1,i)-RX))
!     dFrac  = abs(2d0*dAlpha/myPI)
!     RAD = (1d0-dFrac)*RAD1 + (dFrac)*RAD2
!     dFact = RAD/dist
! !    RAD = RAD1
!     px = 1d0*(RX + dFact*dx)
!     py = 1d0*(RY + dFact*dy)
! !    pz = (1d0-DFZ)*DCORVG(3,i) + DFZ*(RZ + dFact*dz)
!     DCORVG(1,i) = px
!     DCORVG(2,i) = py
! !    DCORVG(3,i) = pz
!    END IF
!   END DO
!  END IF

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


END SUBROUTINE mp_Parametrize
!----------------------------------------------------------------------------------------
SUBROUTINE mp_GetHigherStructures(KVERT1,KEDGE1,KVERT2,KAREA2,KADJ2,NEL1,NVT1,NAT1,NET1)
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

END SUBROUTINE mp_GetHigherStructures

end module mesh_param

