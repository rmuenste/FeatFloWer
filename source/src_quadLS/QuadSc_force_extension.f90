!************************************************************************
      SUBROUTINE GetForces(U1,U2,U3,P,ALPHA,DVISC,KVERT,KAREA,KEDGE,&
                           DCORVG,ELE)
!************************************************************************
!*     Discrete convection operator: Q1~ elements (nonparametric)
!*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid,showID,COMM_SUMMN
      USE var_QuadScalar, ONLY : myFBM,myExport,Properties
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
!
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,&
                 NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
!
      REAL*8  U1(*),U2(*),U3(*),P(*),DVISC(*),DCORVG(NNDIM,*)
      INTEGER ALPHA(*)
      INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      INTEGER KDFG(NNBAS),KDFL(NNBAS)
!
      REAL*8 DResForceX,DResForceY,DResForceZ
      REAL*8 DTrqForceX,DTrqForceY,DTrqForceZ
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
!
      SAVE
!
      IF (myid.eq.0) GOTO 999

      if(myid.eq.1)then
        write(*,*)'Number of particles: ',myFBM%nParticles
      end if
!
      DO 1 I= 1,NNDER
1     BDER(I)=.FALSE.
!
      DO 2 I=1,4
2     BDER(I)=.TRUE.
!
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
!
      ICUB=9
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
!
      DO IP = 1,myFBM%nParticles

      Center = myFBM%particleNew(IP)%Position

      if(myid.eq.1)then
        if(IP.eq.50)then
        write(*,*)'Position: ',myFBM%particleNew(IP)%Position(1),&
                               myFBM%particleNew(IP)%Position(2),&
                               myFBM%particleNew(IP)%Position(3)
        end if
      end if

      DResForceX = 0D0
      DResForceY = 0D0
      DResForceZ = 0D0
      DTrqForceX = 0d0
      DTrqForceY = 0d0
      DTrqForceZ = 0d0
!
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
         ! Count the number of dofs outside
         IF((ALPHA(IG).EQ.0).or.(ALPHA(IG).NE.IP))THEN
          NJALFA=NJALFA+1
         ENDIF
         ! Count the number of dofs inside
         IF (ALPHA(IG).EQ.IP) THEN
          NIALFA=NIALFA+1
         ENDIF
       ENDDO
!--------ONLY FOR THOSE ELEMENTS WHICH HAVE NONZERO GRAD ----------
      IF(NJALFA.EQ.27.OR.NIALFA.EQ.27) GOTO 100
!--------ONLY FOR THOSE ELEMENTS WHICH HAVE NONZERO GRAD ----------
!

!      myExport%p_DataScalarCell(1)%pData(iel)=2
      !write(*,*)'adding IELEval: ',IEL
      nnel = nnel + 1
      DNY = DVISC(IEL)
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
!
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
!
!     Evaluate the solution values and derivatives in the cubature point

       DU1V=0D0     ! U1 value
       DU1X=0D0     ! U1 x deriv
       DU1Y=0D0     ! U1 y deriv
       DU1Z=0D0     ! U1 z deriv
!
       DU2V=0D0     ! U2 value
       DU2X=0D0     ! U2 x deriv
       DU2Y=0D0     ! U2 y deriv
       DU2Z=0D0     ! U2 z deriv
!
       DU3V=0D0     ! U3 value
       DU3X=0D0     ! U3 x deriv
       DU3Y=0D0     ! U3 y deriv
       DU3Z=0D0     ! U3 z deriv
!
       DALV=0D0     ! ALFA value
       DALX=0D0     ! ALFA x deriv
       DALY=0D0     ! ALFA y deriv
       DALZ=0D0     ! ALFA z deriv

       DO I=1,IDFL
         IG=KDFG(I)
         DBI1=DBAS(1,KDFL(I),1) ! value
         DBI2=DBAS(1,KDFL(I),2) ! d/dx
         DBI3=DBAS(1,KDFL(I),3) ! d/dy
         DBI4=DBAS(1,KDFL(I),4) ! d/dz
!---------------FOR U1----------------
         DU1V=DU1V+U1(IG)*DBI1
         DU1X=DU1X+U1(IG)*DBI2
         DU1Y=DU1Y+U1(IG)*DBI3
         DU1Z=DU1Z+U1(IG)*DBI4
!---------------FOR U2----------------
         DU2V=DU2V+U2(IG)*DBI1
         DU2X=DU2X+U2(IG)*DBI2
         DU2Y=DU2Y+U2(IG)*DBI3
         DU2Z=DU2Z+U2(IG)*DBI4
!---------------FOR U3----------------
         DU3V=DU3V+U3(IG)*DBI1
         DU3X=DU3X+U3(IG)*DBI2
         DU3Y=DU3Y+U3(IG)*DBI3
         DU3Z=DU3Z+U3(IG)*DBI4
!---------------FOR ALFA----------------
         IF (ALPHA(IG).EQ.IP) THEN
          DALPHA = 1d0
         ELSE
          DALPHA = 0d0
         END IF
         DALV=DALV+DALPHA*DBI1
         DALX=DALX+DALPHA*DBI2
         DALY=DALY+DALPHA*DBI3
         DALZ=DALZ+DALPHA*DBI4
       ENDDO
!C
!C---------------------------------------------------------
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
       AH1=-Press*DN1 + DNY*(DU1X*DN1 + DU1Y*DN2 + DU1Z*DN3)
       AH2=-Press*DN2 + DNY*(DU2X*DN1 + DU2Y*DN2 + DU2Z*DN3)
       AH3=-Press*DN3 + DNY*(DU3X*DN1 + DU3Y*DN2 + DU3Z*DN3)
!
!        AH1=-P(IEL)*DN1+DNY*((DU1X+DU1X)*DN1+(DU1Y+DU2X)*DN2 + ! full3D
!      *     (DU1Z+DU3X)*DN3)
!        AH2=-P(IEL)*DN2+DNY*((DU2X+DU1Y)*DN1+(DU2Y+DU2Y)*DN2 + ! full3D
!      *     (DU2Z+DU3Y)*DN3)
!        AH3=-P(IEL)*DN3+DNY*((DU3X+DU1Z)*DN1+(DU3Y+DU2Z)*DN2 + ! full3D
!      *     (DU3Z+DU3Z)*DN3)
!
!--------------------------------------------------------------
       DResForceX = DResForceX + AH1*OM
       DResForceY = DResForceY + AH2*OM
       DResForceZ = DResForceZ + AH3*OM
!-------------------Torque force------------------------- 
       XTORQUE = XX - Center(1)
       YTORQUE = YY - Center(2)
       ZTORQUE = ZZ - Center(3)
!
       ATQX = YTORQUE*AH3 - ZTORQUE*AH2                        ! full3D
       ATQY = ZTORQUE*AH1 - XTORQUE*AH3                        ! full3D
       ATQZ = XTORQUE*AH2 - YTORQUE*AH1
!
       DTrqForceX = DTrqForceX + ATQX*OM
       DTrqForceY = DTrqForceY + ATQY*OM
       DTrqForceZ = DTrqForceZ + ATQZ*OM
!
!===============================================================
!--------------UPDATE VELOCITY AND POSITION---------------------
!
200   CONTINUE
!
100   CONTINUE
!
      iPointer = 6*(IP-1)
      myFBM%Force(iPointer+1) = DResForceX
      myFBM%Force(iPointer+2) = DResForceY
      myFBM%Force(iPointer+3) = DResForceZ
      myFBM%Force(iPointer+4) = DTrqForceX
      myFBM%Force(iPointer+5) = DTrqForceY
      myFBM%Force(iPointer+6) = DTrqForceZ

      END DO ! nParticles

999   CALL COMM_SUMMN(myFBM%Force,6*myFBM%nParticles)

      DO IP = 1,myFBM%nParticles
       iPointer = 6*(IP-1)+1

       ! Gather translational force
       myFBM%ParticleNew(IP)%ResistanceForce(1)= &
       Properties%ForceScale(1) * &
       myFBM%Force(iPointer)

       myFBM%ParticleNew(IP)%ResistanceForce(2) = &
       Properties%ForceScale(2) * &
       myFBM%Force(iPointer+1)

       myFBM%ParticleNew(IP)%ResistanceForce(3) = &
       Properties%ForceScale(3) * &
       myFBM%Force(iPointer+2)

       ! Gather rotational force
       myFBM%ParticleNew(IP)%TorqueForce(1) = &
       myFBM%Force(iPointer+3)
       myFBM%ParticleNew(IP)%TorqueForce(2) = &
       myFBM%Force(iPointer+4)
       myFBM%ParticleNew(IP)%TorqueForce(3) = &
       myFBM%Force(iPointer+5)
      END DO

      if(myid.eq.1)then
        write(*,*)'Force(50): ',myFBM%particleNew(50)%ResistanceForce(1),&
                               myFBM%particleNew(50)%ResistanceForce(2),&
                               myFBM%particleNew(50)%ResistanceForce(3)
      end if
      
!     write(*,*)'z-force: ',myFBM%ParticleNew(1)%ResistanceForce(3)
      
99999 CONTINUE

      END
!
!-----------------------------------------------------------------------
!
      SUBROUTINE updateFBM(DensityL,dTime,simTime,Gravity,mfile,myid)
      USE var_QuadScalar, ONLY : myFBM
      USE cinterface, only : FBM_SetParticles 
      INTEGER IP,ipc,myip
      INTEGER myid
      INTEGER mfile
      REAL*8 DensityL,dTime,simTime,Gravity(3)
      REAL*8 pi,volume,mass,massR,radius,dimomir
      PARAMETER (PI=3.1415926535897931D0)
      REAL*8 RForce(3),dVelocity(3),dOmega(3),timecoll
      INTEGER :: iSubSteps
      REAL*8 :: dSubStep,x,y,z

      x = 0
      y = 0
      z = 0
      
      iSubSteps = 1
      call settimestep(dTime)

!     communicate new force + torque
      DO IP = 1,myFBM%nParticles
      ipc = ip-1
      ! Communicate the force
       call setforce(myFBM%particleOld(IP)%ResistanceForce(1),&
       myFBM%particleOld(IP)%ResistanceForce(2),&
       myFBM%particleOld(IP)%ResistanceForce(3),ipc)

      ! Communicate the torque
       call settorque(myFBM%particleOld(IP)%TorqueForce(1),&
       myFBM%particleOld(IP)%TorqueForce(2), &
       myFBM%particleOld(IP)%TorqueForce(3),ipc)

      END DO ! all particles

      ! update velocities by the force determined in the time step
      call settimestep(dTime)      
      call velocityupdate_soft()     
            
      call settimestep(dTime/real(iSubSteps))
      call starttiming()      
      ! start the collision handling
      do iStep=1,iSubSteps
        call stepsoftbody(x,y,z,dtime)
      end do
      call gettiming(timecoll)

      call FBM_SetParticles()

      !-------------------------------------------------------------------------------------
      ! Backup the particle parameters so the old values are available in the next timestep
      !-------------------------------------------------------------------------------------
      DO IP = 1,myFBM%nParticles
      ipc = ip-1

      ! Backup the forces
      myFBM%particleOld(IP)%ResistanceForce = &
      myFBM%particleNew(IP)%ResistanceForce
         
      ! Backup the torque
      myFBM%particleOld(IP)%TorqueForce = &
      myFBM%particleNew(IP)%TorqueForce

      ! Backup the positions
      myFBM%particleOld(IP)%Position = &
      myFBM%particleNew(IP)%Position

      ! Backup the velocities
      myFBM%particleOld(IP)%Velocity = &
      myFBM%particleNew(IP)%Velocity

      ! Backup the angles
      myFBM%particleOld(IP)%Angle = &
      myFBM%particleNew(IP)%Angle

      ! Backup the angular velocity
      myFBM%particleOld(IP)%AngularVelocity = &
      myFBM%particleNew(IP)%AngularVelocity

      END DO ! all particles

      END
!
!-----------------------------------------------------------
!
      SUBROUTINE GetFictKnpr(X,Y,Z,iBndr,inpr,Dist)
      USE var_QuadScalar, ONLY : myFBM
      IMPLICIT NONE
      REAL*8 X,Y,Z,Dist,daux
      REAL*8 PX,PY,PZ,RAD
      INTEGER iBndr,inpr,iP,iaux,ipc,isin

       inpr = 0
       !Dist = 0.05d0
       DO IP = 1,myFBM%nParticles
        PX = myFBM%particleNew(iP)%Position(1)
        PY = myFBM%particleNew(iP)%Position(2)
        PZ = myFBM%particleNew(iP)%Position(3)
        RAD =myFBM%particleNew(iP)%sizes(1)
        ipc=ip-1
        isin = 0

        call insidesoftbody(x,y,z,ipc,isin)

        if(isin .gt. 0)then
          inpr = isin
        end if
       end do

      END SUBROUTINE GetFictKnpr
!
!-----------------------------------------------------------
!
      SUBROUTINE GetVeloFictBCVal(X,Y,Z,ValU,ValV,ValW,IP,t)
      USE var_QuadScalar, ONLY : myFBM,bRefFrame
      IMPLICIT NONE
      INTEGER iP,ipc
      REAL*8 X,Y,Z,t,ValU,ValV,ValW
      REAL*8 PX,PY,PZ,RAD
      REAL*8 Velo(3),Pos(3),Omega(3)
      REAL*8 DVELZ_X,DVELZ_Y,DVELY_Z,DVELY_X,DVELX_Y,DVELX_Z

      ValU = 0d0
      ValV = 0d0
      ValW = 0d0

      if(IP .ne. 0)then
        call getsoftvel(x,y,z,valu,valv,valw,ip)
      end if

      END SUBROUTINE GetVeloFictBCVal
!
!-----------------------------------------------------------
!
      SUBROUTINE GetDistanceALE(X,Y,Z,ale_val,ipc)
      USE var_QuadScalar, ONLY : myFBM
      IMPLICIT NONE
      REAL*8 X,Y,Z
      REAL*8 ale_val
      INTEGER ipc

        call getdistanceid(x,y,z,ale_val,ipc)

      END SUBROUTINE GetDistanceALE
!
!-----------------------------------------------------------
!
      subroutine communicateforce(fx,fy,fz,tx,ty,tz)
      use var_QuadScalar, only : myFBM
      implicit none
      real*8 fx(*),fy(*),fz(*),tx(*),ty(*),tz(*)
      integer IP

      do IP = 1,myFBM%nParticles
        fx(IP)=myFBM%particleNew(IP)%ResistanceForce(1)
        fy(IP)=myFBM%particleNew(IP)%ResistanceForce(2)
        fz(IP)=myFBM%particleNew(IP)%ResistanceForce(3)
        tx(IP)=myFBM%particleNew(IP)%TorqueForce(1)
        ty(IP)=myFBM%particleNew(IP)%TorqueForce(2)
        tz(IP)=myFBM%particleNew(IP)%TorqueForce(3)
      end do
      end subroutine
!
!-----------------------------------------------------------------------
!
