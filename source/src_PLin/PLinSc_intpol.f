      SUBROUTINE IntQ0toQ1(DF0,DF1,DAUX,KVERT,DCORVG,VVOL,NEL,NVT)
      use PP3D_MPI, ONLY: E011Sum
      REAL*4 VVOL(*)
      REAL*8 DF0(*),DF1(*),DAUX(*)
      INTEGER KVERT(8,*),NEL,NVT,IVT,IVERT
      INTEGER IEL

      DO IEL=1,NEL
       DVAL = DF0(IEL)
       DVOL = DBLE(VVOL(IEL))*0.125d0
       DO IVT=1,8
        IVERT = KVERT(IVT,IEL)
        DF1 (IVERT) = DF1 (IVERT) + DVAL*DVOL
        DAUX(IVERT) = DAUX(IVERT) + DVOL
       END DO
      END DO

      CALL E011SUM(DF1)
      CALL E011SUM(DAUX)

      DO IVT=1,NVT
       DF1(IVT) = DF1(IVT)/DAUX(IVT)
      END DO

      END


      SUBROUTINE IntP1toQ1(DPHIQ,DPHIP,DWEIGHT,KVERT,
     *           DCORVG,AVOL,ELE)
      use PP3D_MPI, ONLY: E011Sum
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)

      REAL*4  AVOL(*)
      REAL*8  DPHIP(*),DPHIQ(*),DWEIGHT(*),DCORVG(3,*)
      INTEGER KVERT(8,*)

      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)

      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,
     *                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),
     *                IEL,NDIM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD

      DO I= 1,NNDER
       BDER(I)=.FALSE.
      END DO

      BDER(1)=.TRUE.

      DO IEL=1,NEL

      DVOL=DBLE(AVOL(IEL))

! *** Evaluation of coordinates of the vertices
      DX0 = 0d0
      DY0 = 0d0
      DZ0 = 0d0
      DO IVE=1,NVE
       IP=KVERT(IVE,IEL)
       KVE(IVE)=IP
       DX(IVE)=DCORVG(1,IP)
       DY(IVE)=DCORVG(2,IP)
       DZ(IVE)=DCORVG(3,IP)
       DX0 = DX0 + 0.125d0*DX(IVE)
       DY0 = DY0 + 0.125d0*DY(IVE)
       DZ0 = DZ0 + 0.125d0*DZ(IVE)
      END DO

      DO IVE=1,8
       IP=KVE(IVE)
       CALL ELE(DX(IVE)-DX0, DY(IVE)-DY0, DZ(IVE)-DZ0, 0)

       DPHI=0D0
       DO JDOFE=1,4
        JDFL=JDOFE
        JDFG=JDOFE+4*(IEL -1)
        HBASJ=DBAS(1,JDFL,1)
        DPHI=DPHI+DPHIP(JDFG)*HBASJ
       END DO

       DWEIGHT(IP) = DWEIGHT(IP) + 0.125D0*DVOL
       DPHIQ(IP)   = DPHIQ(IP)   + 0.125D0*DVOL*DPHI

      END DO

      END DO ! IEL

      CALL E011SUM(DPHIQ)
      CALL E011SUM(DWEIGHT)

      DO IVT=1,NVT
       DPHIQ(IVT) = DPHIQ(IVT)/DWEIGHT(IVT)
      END DO

      END
C
C
C
      SUBROUTINE IntPolNormals(dNorm,dVal,Layer,kVert,dCorvg,aVol,
     *           dMid,nel,nvt)
      use PP3D_MPI, ONLY: E011Sum
      IMPLICIT NONE

      REAL*4  aVol(*)
      REAL*8  dNorm(3,nvt),dVal(4,*),dCorvg(3,*),dMid(3,*)
      INTEGER kVert(8,*),Layer(*),nel,nvt
      REAL*8  daux,ds,dV,dX(3),dWeight
      INTEGER iel,ivt,iVert

       ! Normal directions are weighted from the old time step 
      dWeight = 0d0
      DO iVert=1,nvt
       dNorm(1,iVert) = dWeight*dNorm(1,iVert)
       dNorm(2,iVert) = dWeight*dNorm(2,iVert)
       dNorm(3,iVert) = dWeight*dNorm(3,iVert)
      END DO

      DO iel=1,nel
       IF (ABS(Layer(iel)).NE.100) THEN
        daux=DBLE(aVol(iel))
!         IF (Layer(iel).NE.0) THEN
         daux = 1d0!DSIGN(daux,DBLE(Layer(iel)))
         DO ivt=1,8
          iVert = kVert(ivt,iel)
          dNorm(1,iVert) = dNorm(1,iVert) + 0.125D0*daux*dVal(2,iel)
          dNorm(2,iVert) = dNorm(2,iVert) + 0.125D0*daux*dVal(3,iel)
          dNorm(3,iVert) = dNorm(3,iVert) + 0.125D0*daux*dVal(4,iel)
         END DO !ivt
        ELSE
         DO ivt=1,8
          iVert = kVert(ivt,iel)
          dX(1) = dCorvg(1,iVert) - dMid(1,iel)
          dX(2) = dCorvg(2,iVert) - dMid(2,iel)
          dX(3) = dCorvg(3,iVert) - dMid(3,iel)
          dV = dVal(1,iel)       + dVal(2,iel)*dX(1)
     *       + dVal(3,iel)*dX(2) + dVal(4,iel)*dX(3)
          ds = DSIGN(daux,dV)
          dNorm(1,iVert) = dNorm(1,iVert) + 0.125D0*ds*dVal(2,iel)
          dNorm(2,iVert) = dNorm(2,iVert) + 0.125D0*ds*dVal(3,iel)
          dNorm(3,iVert) = dNorm(3,iVert) + 0.125D0*ds*dVal(4,iel)
         END DO !ivt
!         END IF
       END IF
      END DO ! iel

      CALL E011SUM(dNorm(1,:))
      CALL E011SUM(dNorm(2,:))
      CALL E011SUM(dNorm(3,:))

      DO iVert=1,nvt
       daux = SQRT(dNorm(1,iVert)**2d0 + dNorm(2,iVert)**2d0 +
     *             dNorm(3,iVert)**2d0)
       IF (daux.NE.0d0) THEN
        dNorm(1,iVert) = dNorm(1,iVert)/daux
        dNorm(2,iVert) = dNorm(2,iVert)/daux
        dNorm(3,iVert) = dNorm(3,iVert)/daux
       ELSE
        dNorm(1,iVert) = 0d0
        dNorm(2,iVert) = 0d0
        dNorm(3,iVert) = 0d0
       END IF
      END DO

      END
C
C
C
      SUBROUTINE IntQ1toQ0toQ1(DPHIQ,DWEIGHT,KVERT,AVOL,NEL,NVT)
      use PP3D_MPI, ONLY: E011Sum
      IMPLICIT NONE

      REAL*4  AVOL(*)
      REAL*8  DPHIQ(*),DWEIGHT(*)
      INTEGER KVERT(8,*),NEL,NVT,IEL,IVT,IVE,IVERT
      REAL*8  DVAL(NEL),DVOL

      DO IEL=1,NEL
       DVAL(IEL) = 0d0
       DO IVE=1,8
        IVERT=KVERT(IVE,IEL)
        DVAL(IEL) = DVAL(IEL) + 0.125d0*DPHIQ(IVERT)
       END DO
      END DO

      DO IVT=1,NVT
       DPHIQ(IVT) = 0d0
       DWEIGHT(IVT) = 0d0
      END DO

      DO IEL=1,NEL
       DVOL = 0.125d0*DBLE(AVOL(IEL))
!        write(*,*) dvol,dval(iel)
       DO IVE=1,8
        IVERT=KVERT(IVE,IEL)
        DPHIQ  (IVERT) = DPHIQ  (IVERT) + DVOL*DVAL(IEL)
        DWEIGHT(IVERT) = DWEIGHT(IVERT) + DVOL
       END DO
      END DO

      CALL E011SUM(DPHIQ)
      CALL E011SUM(DWEIGHT)

      DO IVT=1,NVT
       DPHIQ(IVT) = DPHIQ(IVT)/DWEIGHT(IVT)
      END DO

      END
