!************************************************************************
      SUBROUTINE GetTorqueMixer(U1,U2,U3,P,ALPHA,KVERT,KAREA,KEDGE,&
                                DCORVG,DVISC,DTrqForce,ELE,iref)
!************************************************************************
      USE PP3D_MPI, ONLY:myid,showID,COMM_SUMMN
      USE var_QuadScalar, ONLY : myFBM
      USE var_QuadScalar, ONLY: my1DTorque
      USE Sigma_User, ONLY: myOutput,mySigma

      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
!
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,&
                 NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
!
      REAL*8  U1(*),U2(*),U3(*),P(*),DCORVG(NNDIM,*),DVISC(*)
      INTEGER ALPHA(*)
      INTEGER IREF
      INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      INTEGER KDFG(NNBAS),KDFL(NNBAS)
!
      REAL*8 DResForceX,DResForceY,DResForceZ
      REAL*8 DTrqForce(3)
      REAL*8 Center(3),dForce(6)
!
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,&
                      DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),&
                      IEL,NDIM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,&
                      NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
!
! *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)
      EQUIVALENCE (IAUSAV,VIPARM)
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,&
                    IMASS,IMASSL,IUPW,IPRECA,IPRECB,&
                    ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,&
                    INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,&
                    ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,&
                    IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA

      SAVE

      if (.not.allocated(my1DTorque)) THEN
       IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE") THEN 
        allocate(my1DTorque(2,myOutput%nOf1DLayers))
       ELSE
        allocate(my1DTorque(1,myOutput%nOf1DLayers))
       END IF
      END IF
            
      IF (iref.eq.101) my1DTorque(1,:) = 0d0
      IF (iref.eq.102) my1DTorque(2,:) = 0d0
      IF (iref.eq.103) my1DTorque(1,:) = 0d0
      
      IF (myid.eq.0) GOTO 999

      IF (iref.eq.101) dSignedAxisDist = +mySigma%a/2d0
      IF (iref.eq.102) dSignedAxisDist = -mySigma%a/2d0
      IF (iref.eq.103) dSignedAxisDist = 0d0

      DO 1 I= 1,NNDER

1     BDER(I)=.FALSE.

      DO 2 I=1,4
2     BDER(I)=.TRUE.

      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)!set number of element
      IDFL=NDFL(IELTYP)


      ICUB=9
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999

      IP = IREF

      DTrqForce(1) = 0d0
      DTrqForce(2) = 0d0
      DTrqForce(3) = 0d0

! *** Loop over all elements
      nnel = 0
      DO 100 IEL=1,NEL
!
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
!
!-----------------------------------------------------------------
       NJALFA=0
       NIALFA=0
       DO I=1,IDFL
         IG=KDFG(I)
         IF (ALPHA(IG).EQ.iref) THEN
          NJALFA=NJALFA+1
         ELSE
          NIALFA=NIALFA+1
         ENDIF
       ENDDO
!--------ONLY FOR THOSE ELEMENTS WHICH HAVE NONZERO GRAD ----------
      IF(NJALFA.EQ.27.OR.NIALFA.EQ.27) GOTO 100
!--------ONLY FOR THOSE ELEMENTS WHICH HAVE NONZERO GRAD ----------
!
      nnel = nnel + 1
!
! *** Evaluation of coordinates of the vertices
      DX0 = 0d0
      DY0 = 0d0
      DZ0 = 0d0
      DO 120 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
      DZ(IVE)=DCORVG(3,JP)
      DX0 = DX0 + 0.125d0*DX(IVE)
      DY0 = DY0 + 0.125d0*DY(IVE)
      DZ0 = DZ0 + 0.125d0*DZ(IVE)
120   CONTINUE
!
      DJ11=( DX(1)+DX(2)+DX(3)+DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
      DJ12=( DY(1)+DY(2)+DY(3)+DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
      DJ13=( DZ(1)+DZ(2)+DZ(3)+DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
      DJ21=(-DX(1)+DX(2)+DX(3)-DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
      DJ22=(-DY(1)+DY(2)+DY(3)-DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
      DJ23=(-DZ(1)+DZ(2)+DZ(3)-DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
      DJ31=(-DX(1)-DX(2)+DX(3)+DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
      DJ32=(-DY(1)-DY(2)+DY(3)+DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
      DJ33=(-DZ(1)-DZ(2)+DZ(3)+DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
      DJ41=(-DX(1)-DX(2)-DX(3)-DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
      DJ42=(-DY(1)-DY(2)-DY(3)-DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
      DJ43=(-DZ(1)-DZ(2)-DZ(3)-DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
      DJ51=( DX(1)-DX(2)+DX(3)-DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
      DJ52=( DY(1)-DY(2)+DY(3)-DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
      DJ53=( DZ(1)-DZ(2)+DZ(3)-DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
      DJ61=( DX(1)-DX(2)-DX(3)+DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
      DJ62=( DY(1)-DY(2)-DY(3)+DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
      DJ63=( DZ(1)-DZ(2)-DZ(3)+DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
      DJ71=( DX(1)+DX(2)-DX(3)-DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
      DJ72=( DY(1)+DY(2)-DY(3)-DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
      DJ73=( DZ(1)+DZ(2)-DZ(3)-DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
      DJ81=(-DX(1)+DX(2)-DX(3)+DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
      DJ82=(-DY(1)+DY(2)-DY(3)+DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
      DJ83=(-DZ(1)+DZ(2)-DZ(3)+DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
!
      CALL ELE(0D0,0D0,0D0,-2)
      IF (IER.LT.0) GOTO 99999
!
! *** Loop over all cubature points
      DO 200 ICUBP=1,NCUBP
!
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
!
! *** Jacobian of the bilinear mapping onto the reference element
      DJAC(1,1)=DJ21+DJ51*XI2+DJ61*XI3+DJ81*XI2*XI3
      DJAC(1,2)=DJ31+DJ51*XI1+DJ71*XI3+DJ81*XI1*XI3
      DJAC(1,3)=DJ41+DJ61*XI1+DJ71*XI2+DJ81*XI1*XI2
      DJAC(2,1)=DJ22+DJ52*XI2+DJ62*XI3+DJ82*XI2*XI3
      DJAC(2,2)=DJ32+DJ52*XI1+DJ72*XI3+DJ82*XI1*XI3
      DJAC(2,3)=DJ42+DJ62*XI1+DJ72*XI2+DJ82*XI1*XI2
      DJAC(3,1)=DJ23+DJ53*XI2+DJ63*XI3+DJ83*XI2*XI3
      DJAC(3,2)=DJ33+DJ53*XI1+DJ73*XI3+DJ83*XI1*XI3
      DJAC(3,3)=DJ43+DJ63*XI1+DJ73*XI2+DJ83*XI1*XI2
      DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))&
           -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))&
           +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      OM=DOMEGA(ICUBP)*ABS(DETJ)
!
      XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
      YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
      ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2
      dRR = sqrt(XX*XX + YY*YY)
!
!       dddu1 = - (250d0/60d0)*2d0*PI*YY
!       dddu2 = + (250d0/60d0)*2d0*PI*XX
!      dddpX = + 0.5d0*(dConst*dRR)**2d0

      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
!
!     Evaluate the solution values and derivatives in the cubature point
       DNY =0d0                  ! Viscosity
!
       DU1X=0D0                  ! U1 x deriv
       DU1Y=0d0                  ! U1 y deriv
       DU1Z=0D0                  ! U1 z deriv
!
       DU2X=0d0                  ! U2 x deriv
       DU2Y=0D0                  ! U2 y deriv
       DU2Z=0D0                  ! U2 z deriv
!
       DU3X=0D0                  ! U3 x deriv
       DU3Y=0D0                  ! U3 y deriv
       DU3Z=0D0                  ! U3 z deriv
!
       DALX=0D0                  ! ALFA x deriv
       DALY=0D0                  ! ALFA y deriv
       DALZ=0D0                  ! ALFA z deriv

!        Press =0d0

       DO I=1,IDFL
         IG=KDFG(I)
         DBI1=DBAS(1,KDFL(I),1)
         DBI2=DBAS(1,KDFL(I),2)
         DBI3=DBAS(1,KDFL(I),3)
         DBI4=DBAS(1,KDFL(I),4)
!---------------Viscosity----------------
         DNY=DNY + DVISC(IG)*DBI1
! !---------------Pressure----------------
!          Press=Press + P(IG)*DBI1
!---------------FOR U1----------------
         DU1X=DU1X + U1(IG)*DBI2
         DU1Y=DU1Y + U1(IG)*DBI3
         DU1Z=DU1Z + U1(IG)*DBI4
!---------------FOR U2----------------
         DU2X=DU2X + U2(IG)*DBI2
         DU2Y=DU2Y + U2(IG)*DBI3
         DU2Z=DU2Z + U2(IG)*DBI4
!---------------FOR U3----------------
         DU3X=DU3X + U3(IG)*DBI2
         DU3Y=DU3Y + U3(IG)*DBI3
         DU3Z=DU3Z + U3(IG)*DBI4
!---------------FOR ALFA----------------
         IF (ALPHA(IG).EQ.iref) THEN
          DALPHA = 1d0
         ELSE
          DALPHA = 0d0
         END IF
         DALX=DALX + DALPHA*DBI2
         DALY=DALY + DALPHA*DBI3
         DALZ=DALZ + DALPHA*DBI4
       ENDDO
!
!---------------------------------------------------------
       JJ = 4*(IEL-1) + 1
       Press =          P(JJ  ) + (XX-DX0)*P(JJ+1) + &
               (YY-DY0)*P(JJ+2) + (ZZ-DZ0)*P(JJ+3)
!--------------------------------------------------------
!-----------Form the integrand------------------
       DN1=-DALX
       DN2=-DALY
       DN3=-DALZ
!-------------------Acting force-------------------------
!--------------Deformation calculation-------------
!
!        AH1=-Press*DN1 + dny*(DU1X*DN1 + DU1Y*DN2 + DU1Z*DN3)
!        AH2=-Press*DN2 + dny*(DU2X*DN1 + DU2Y*DN2 + DU2Z*DN3)
!        AH3=-Press*DN3 + dny*(DU3X*DN1 + DU3Y*DN2 + DU3Z*DN3)

       AH1=-Press*DN1 + DNY*((DU1X+DU1X)*DN1+(DU1Y+DU2X)*DN2+(DU1Z+DU3X)*DN3)
       AH2=-Press*DN2 + DNY*((DU2X+DU1Y)*DN1+(DU2Y+DU2Y)*DN2+(DU2Z+DU3Y)*DN3)
       AH3=-Press*DN3 + DNY*((DU3X+DU1Z)*DN1+(DU3Y+DU2Z)*DN2+(DU3Z+DU3Z)*DN3)

!-------------------Torque force------------------------- 
       IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE") THEN

        IF (ADJUSTL(TRIM(mySigma%RotationAxis)).EQ."PARALLEL") THEN

         XTORQUE = XX
         YTORQUE = YY - dSignedAxisDist
         ZTORQUE = 0D0
        
        ELSE
         IF (iref.eq.101) THEN
          CALL TransformPointToNonparallelRotAxis(xx,yy,zz,XXX,YYY,ZZZ,-1d0)
         END IF
         IF (iref.eq.102) THEN
          CALL TransformPointToNonparallelRotAxis(xx,yy,zz,XXX,YYY,ZZZ,+1d0)
         END IF
       
         XTORQUE = XXX
         YTORQUE = YYY
         ZTORQUE = 0D0
        END IF
       END IF

       IF (ADJUSTL(TRIM(mySigma%cType)).EQ."SSE") THEN
        XTORQUE = XX
        YTORQUE = YY
        ZTORQUE = 0D0
       END IF
!
       ATQX = YTORQUE*AH3 - ZTORQUE*AH2                        ! full3D
       ATQY = ZTORQUE*AH1 - XTORQUE*AH3                        ! full3D
       ATQZ = XTORQUE*AH2 - YTORQUE*AH1
!
       DTrqForce(1) = DTrqForce(1) + ATQX*OM
       DTrqForce(2) = DTrqForce(2) + ATQY*OM
       DTrqForce(3) = DTrqForce(3) + ATQZ*OM
       
       CALL FindSegmentForTorque(ZZ,mySigma%L,myZPosition,myOutput%nOf1DLayers)
       IF (iref.eq.102) THEN
        my1DTorque(2,myZPosition) = my1DTorque(2,myZPosition) + ATQZ*OM
       ELSE
        my1DTorque(1,myZPosition) = my1DTorque(1,myZPosition) + ATQZ*OM
       END IF
      
200   CONTINUE
100   CONTINUE

!      END DO ! 1 and 2 screw

999   CALL COMM_SUMMN(DTrqForce,3)
      
      IF (iref.eq.102) THEN
       CALL COMM_SUMMN(my1DTorque(2,:),myOutput%nOf1DLayers)
      ELSE
       CALL COMM_SUMMN(my1DTorque(1,:),myOutput%nOf1DLayers)
      END IF

99999 CONTINUE
!
!       write(*,*)  myOutput%nOf1DLayers,myid,mySigma%L, ADJUSTL(TRIM(mySigma%cType))
!       if (myid.eq.1) then
!        write(*,*)  my1DTorque(1,:)
!        write(*,*)  
!        write(*,*)  my1DTorque(2,:)
!       end if
!       pause
      
      END
!      
      SUBROUTINE FindSegmentForTorque(zPos,zLen,iPos,nPos)
      implicit none
      INTEGER iPos,nPos
      REAL*8  zPos,zLen
      INTEGER i
      REAL*8  zMin,zMax
      
      iPos = -1
      DO i=1,nPos
       zMin = DBLE(i-1)*(1.0001d0*zLen)/DBLE(nPos)
       zMax = DBLE(i  )*(1.0001d0*zLen)/DBLE(nPos)
       IF (zPos.ge.zMin.and.zPos.lt.zMax) THEN
        iPos = i
        EXIT
       END IF
      END DO
       
      IF (iPos.lt.0) THEN
       WRITE(*,*) "Torque Value cannot be stored ... ",zPos,zLen,iPos,nPos
       STOP
      END IF
      
      END 
