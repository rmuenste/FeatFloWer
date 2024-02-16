!#define OUTPUT_LEVEL2 
!#define SED_BENCH
!************************************************************************

function calculate_l2_norm(fx, fy, fz)
  real(kind=8), intent(in) :: fx, fy, fz
  real(kind=8) :: calculate_l2_norm

  calculate_l2_norm = sqrt(fx**2 + fy**2 + fz**2)

end function calculate_l2_norm

!************************************************************************
!
!************************************************************************
SUBROUTINE ForcesLocalParticles(factors,U1,U2,U3,P,ALPHA,DVISC,KVERT,KAREA,KEDGE,&
                                DCORVG,ELE, maxLocal)
!************************************************************************
!*     Discrete convection operator: Q1~ elements (nonparametric)
!*-----------------------------------------------------------------------
USE PP3D_MPI, ONLY:myid,showID,COMM_SUMMN
USE var_QuadScalar, ONLY : myExport,Properties,FictKNPR_uint64, FictKNPR
use cinterface
use dem_query
IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
CHARACTER SUB*6,FMT*15,CPARAM*120

PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,&
           NNDIM=3,NNCOF=10)
PARAMETER (Q2=0.5D0,Q8=0.125D0)

! An array of scaling factors
Real*8 :: factors(*)

! U/V/W velocity components
REAL*8 :: U1(*),U2(*),U3(*)

! Pressure, viscosity
Real*8 :: P(*),DVISC(*)

! Coordinates
Real*8 :: DCORVG(NNDIM,*)

! Alpha function
integer :: ALPHA(*)

! The maximum force magnitude
REAL*8, intent(inout) :: maxLocal

INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
INTEGER KDFG(NNBAS),KDFL(NNBAS)

REAL*8 :: DResForceX,DResForceY,DResForceZ
REAL*8 :: DTrqForceX,DTrqForceY,DTrqForceZ
REAL*8 :: Center(3),dForce(6)

type(tParticleData), dimension(:), allocatable :: theParticles
integer :: numParticles, particleId

logical :: crit1, crit2

real*8 :: theNorm, totalMax
!integer, dimension(1) :: processRanks
!integer :: MPI_W0, MPI_EX0
!integer :: MPI_Comm_EX0, new_comm
!integer :: error_indicator

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

!user COMMON blocks
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

!  processRanks(1) = 0
!  CALL MPI_COMM_GROUP(MPI_COMM_WORLD, MPI_W0, error_indicator)
!  CALL MPI_GROUP_EXCL(MPI_W0, 1, processRanks, MPI_EX0, error_indicator)
!  CALL MPI_COMM_CREATE(MPI_COMM_WORLD, MPI_EX0, MPI_Comm_EX0, error_indicator)
 
  maxLocal = 0.0

  IF (myid.eq.0) return! GOTO 999

#ifdef SED_BENCH
  if(myid.eq.1)then
    write(*,*) 'Force with SED BENCH settings!'
  end if
#endif
  
  DO I= 1,NNDER
    BDER(I)=.FALSE.
  end do
  
  DO I=1,4
    BDER(I)=.TRUE.
  end do
  
  IELTYP=-1
  CALL ELE(0D0,0D0,0D0,IELTYP)
  IDFL=NDFL(IELTYP)
  
  ICUB=9
  CALL CB3H(ICUB)
  IF (IER.ne.0)then
    call ExitError('Error in GetForces',1702)
  end if

  numParticles = numLocalParticles()
  if (numParticles .eq. 0)return

  if(.not. allocated(theParticles)) then
    allocate(theParticles(numLocalParticles())) 
  else if ((allocated(theParticles)).and.(size(theParticles) .ne. numLocalParticles()))then
    deallocate(theParticles)
    allocate(theParticles(numLocalParticles())) 
  end if

  !=====================================================================
  ! We loop over all particles first
  !=====================================================================
#ifdef OUTPUT_LEVEL2 
  if (numParticles > 0) then
    write(*,'(I3,A,I5,A)')myid,'> FBM Force Calculation for:', numParticles, 'local particles'
  end if
#endif

  call getAllParticles(theParticles)

  DO IP = 1,numParticles
#ifdef OUTPUT_LEVEL2 
  write(*,'(I3,A,I5)')myid,'>Particle: ',theParticles(IP)%bytes(1) + 1

#endif

  particleId = theParticles(IP)%localIdx
  particleId = theParticles(IP)%bytes(1) + 1
  
  Center = theParticles(IP)%position 
  
  DResForceX = 0D0
  DResForceY = 0D0
  DResForceZ = 0D0
  DTrqForceX = 0d0
  DTrqForceY = 0d0
  DTrqForceZ = 0d0
  
  nnel = 0
  !=====================================================================
  ! loop over all elements 
  !=====================================================================
  DO IEL=1,NEL
  
  CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
  IF (IER.LT.0)then
    call ExitError('Error in GetForces',1728)
  end if 
  
   NJALFA=0
   NIALFA=0
   DO I=1,IDFL
     IG=KDFG(I)
     ! map particle Id to systemId, then
     ! map IG to see if IG is in systemId
     
     crit1 = .false.
     crit2 = .false.
     IF((ALPHA(IG).EQ.0).or.(ALPHA(IG).NE.particleId))THEN
       crit1 = .true.
     end if

     ! The preferred method of checking is to:
     ! compare the particle system ID (which is a 64bit unsigned integer) to
     ! the fortran long representation of the FictKNPR array
     IF((ALPHA(IG).EQ.0).or.(.not. longIdMatch(IG, theParticles(IP)%bytes)))THEN
     !IF((ALPHA(IG).EQ.0).or.(.not. map_local_to_system(particleId,IG)))THEN
     !IF((ALPHA(IG).EQ.0))THEN
      crit2 = .true.
      NJALFA=NJALFA+1
     ENDIF

!     if(.not. (crit1 .eqv. crit2))then
!       write(*,*) 'myid:',myid,'> a(ig)=',alpha(ig), ' partId=',particleId, ' partLongId=',theParticles(ip)%bytes,' uintFict',&
!         FictKNPR_uint64(ig)%bytes, ' F=',FictKNPR(ig)
!!       call ExitError('Error in GetForces',157)
!     end if

     crit1 = .false.
     crit2 = .false.
     IF (ALPHA(IG).eq.particleId) THEN
       crit1 = .true.
     end if
     IF (longIdMatch(IG, theParticles(IP)%bytes)) THEN
      crit2 = .true.
     !IF ((ALPHA(IG).EQ.1)) THEN
      NIALFA=NIALFA+1
     ENDIF

 !    if(.not. (crit1 .eqv. crit2))then
 !      write(*,*) 'myid:',myid,'> a(ig)=',alpha(ig), ' partId=',particleId, ' partLongId=',theParticles(ip)%bytes,' uintFict',&
 !        FictKNPR_uint64(ig)%bytes, ' F=',FictKNPR(ig)
 !!      call ExitError('Error in GetForces',174)
 !    end if

   ENDDO

  ! Skip elements where all dofs are inside or
  ! all dofs are outside
  IF(NJALFA.EQ.27.OR.NIALFA.EQ.27)then
!    write(*,*)myid,'>Particle: ',IP, '| njalfa: ', njalfa, ' nialfa :', nialfa
    cycle
  end if
!  write(*,*)myid,'>Particle: ',IP, ' iel=', iel, '| njalfa: ', njalfa, ' nialfa :', nialfa

  nnel = nnel + 1
  DNY = DVISC(IEL)

  ! *** Evaluation of coordinates of the vertices
  DX0 = 0d0
  DY0 = 0d0
  DZ0 = 0d0

  ! Calculate the element center in d0
  ! and get the element vertices in dx
  DO IVE=1,NVE
    JP=KVERT(IVE,IEL)
    KVE(IVE)=JP
    DX(IVE)=DCORVG(1,JP)
    DY(IVE)=DCORVG(2,JP)
    DZ(IVE)=DCORVG(3,JP)
    DX0 = DX0 + 0.125d0*DX(IVE)
    DY0 = DY0 + 0.125d0*DY(IVE)
    DZ0 = DZ0 + 0.125d0*DZ(IVE)
  end do

!==============================================================    
!                  Calculation of the jacobian
!==============================================================    
!     Precompute averaged values for later use 
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

  ! Initialize the ELE
  CALL ELE(0D0,0D0,0D0,-2)

  IF (IER.LT.0)then
    call ExitError('Error in GetForces',1800)
  end if 

    ! Loop over all cubature points
    DO ICUBP=1,NCUBP
  
    ! Get the coordinate of the cubature point
    XI1=DXI(ICUBP,1)
    XI2=DXI(ICUBP,2)
    XI3=DXI(ICUBP,3)

    ! Jacobian of the bilinear mapping onto the reference element
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
  
    ! Now calculate the real coordinates of the cubature point
    XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
    YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
    ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2
  
    ! Call ELE with the cubature point to SET the 
    ! derivatives in the DBAS-Array
    CALL ELE(XI1,XI2,XI3,-3)

    IF (IER.LT.0)then
      call ExitError('Error in GetForces',1837)
    end if 

    ! Evaluate the solution values and derivatives in the cubature point
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
      DBI1=DBAS(1,KDFL(I),1)
      DBI2=DBAS(1,KDFL(I),2)
      DBI3=DBAS(1,KDFL(I),3)
      DBI4=DBAS(1,KDFL(I),4)
  
      !------FOR U1--------
      DU1V=DU1V+U1(IG)*DBI1
      DU1X=DU1X+U1(IG)*DBI2
      DU1Y=DU1Y+U1(IG)*DBI3
      DU1Z=DU1Z+U1(IG)*DBI4
  
      !------FOR U2--------
      DU2V=DU2V+U2(IG)*DBI1
      DU2X=DU2X+U2(IG)*DBI2
      DU2Y=DU2Y+U2(IG)*DBI3
      DU2Z=DU2Z+U2(IG)*DBI4
  
      !------FOR U3--------
      DU3V=DU3V+U3(IG)*DBI1
      DU3X=DU3X+U3(IG)*DBI2
      DU3Y=DU3Y+U3(IG)*DBI3
      DU3Z=DU3Z+U3(IG)*DBI4
  
      !------FOR ALFA------
      !IF (ALPHA(IG).EQ.particleId) THEN
      IF (longIdMatch(IG, theParticles(IP)%bytes)) THEN
       DALPHA = 1d0
      ELSE
       DALPHA = 0d0
      END IF
      DALV=DALV+DALPHA*DBI1
      DALX=DALX+DALPHA*DBI2
      DALY=DALY+DALPHA*DBI3
      DALZ=DALZ+DALPHA*DBI4
    ENDDO
  
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
    end do ! end loop cubature points

  end do ! end loop elements

#ifdef SED_BENCH
  DResForceX = 0.0
  DResForceY = 0.0
  DResForceZ = 4.0 * DResForceZ
  theParticles(ip)%force(:) = (/DResForceX, DResForceY, DResForceZ/)
  theParticles(ip)%torque(:) = (/0.0, 0.0, 0.0/)
#else
  theParticles(ip)%force(:) = (/DResForceX, DResForceY, DResForceZ/)
  theParticles(ip)%torque(:) = (/DTrqForceX, DTrqForceY, DTrqForceZ/)
#endif

  theNorm = calculate_l2_norm(DResForceX, DResForceY, DResForceZ)
  if( theNorm >= maxLocal ) then
    maxLocal = theNorm
  end if

#ifdef OUTPUT_LEVEL2 
  write(*,'(I3,A,I5)')myid,'>Particle: ',theParticles(IP)%bytes(1) + 1
  write(*,'(I3,A,I5,A,3D12.4,A,I10)')myid,' pidx=', theParticles(ip)%bytes(1) + 1, ' theForce: ', (/DResForceX, DResForceY, DResForceZ/), ' elems: ', nnel
!  write(*,*)'TheForce: ', (/DResForceX, DResForceY, DResForceZ/)
!  write(*,*)myid,' nnel: ',nnel 
  !write(*,*)'TheTorque: ', theParticles(ip)%torque
#endif

  ! This function is in the dem_query module
  call setForcesMapped(theParticles(ip))

  END DO ! nParticles

!  call MPI_Reduce(totalMax, maxForce, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_EX0, ierr)
!  if (myid .eq. 1)then
!    write(*,*)'maxFluidForce1: ', maxLocal 
!  end if

END SUBROUTINE ForcesLocalParticles
!************************************************************************
!
!************************************************************************
SUBROUTINE ForcesRemoteParticles(factors,U1,U2,U3,P,ALPHA,DVISC,KVERT,KAREA,KEDGE,&
                        DCORVG,ELE, maxLocal)
!************************************************************************
!*     Discrete convection operator: Q1~ elements (nonparametric)
!*-----------------------------------------------------------------------
USE PP3D_MPI, ONLY:myid,showID,COMM_SUMMN
USE var_QuadScalar, ONLY : myExport,Properties,FictKNPR_uint64, FictKNPR
use cinterface
use dem_query
IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
CHARACTER SUB*6,FMT*15,CPARAM*120

PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,&
           NNDIM=3,NNCOF=10)
PARAMETER (Q2=0.5D0,Q8=0.125D0)

! An array of scaling factors
Real*8 :: factors(*)

! U/V/W velocity components
REAL*8 :: U1(*),U2(*),U3(*)

! Pressure, viscosity
Real*8 :: P(*),DVISC(*)

! Coordinates
Real*8 :: DCORVG(NNDIM,*)

! Alpha function
integer :: ALPHA(*)

! The maximum force magnitude
REAL*8, intent(inout) :: maxLocal

INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
INTEGER KDFG(NNBAS),KDFL(NNBAS)

REAL*8 :: DResForceX,DResForceY,DResForceZ
REAL*8 :: DTrqForceX,DTrqForceY,DTrqForceZ
REAL*8 :: Center(3),dForce(6)

type(tParticleData), dimension(:), allocatable :: theParticles
integer :: numParticles, particleId
real*8 :: pressSum, momSum
logical :: crit1, crit2

real*8 :: theNorm, totalMax
!integer, dimension(1) :: processRanks
!integer :: MPI_W0, MPI_EX0, ierr
!integer :: MPI_Comm_EX0, new_comm
!integer :: error_indicator

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

!user COMMON blocks
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
 
  localMax = 0.0
!  processRanks(1) = 0
!  CALL MPI_COMM_GROUP(MPI_COMM_WORLD, MPI_W0, error_indicator)
!  CALL MPI_GROUP_EXCL(MPI_W0, 1, processRanks, MPI_EX0, error_indicator)
!  CALL MPI_COMM_CREATE(MPI_COMM_WORLD, MPI_EX0, MPI_Comm_EX0, error_indicator)

  IF (myid.eq.0) return! GOTO 999
  
  pressSum = 0.0
  DO I= 1,NNDER
    BDER(I)=.FALSE.
  end do
  
  DO I=1,4
    BDER(I)=.TRUE.
  end do
  
  IELTYP=-1
  CALL ELE(0D0,0D0,0D0,IELTYP)
  IDFL=NDFL(IELTYP)
  
  ICUB=9
  CALL CB3H(ICUB)
  IF (IER.ne.0)then
    call ExitError('Error in GetForces',1702)
  end if

  numParticles = numRemParticles()
  if (numParticles .eq. 0)return

  if(allocated(theParticles))then
    deallocate(theParticles)
  end if
  allocate(theParticles(numParticles)) 

!  if(.not. allocated(theParticles)) then
!    allocate(theParticles(numParticles)) 
!  else if ((allocated(theParticles)).and.(size(theParticles) .ne. numParticles))then
!    deallocate(theParticles)
!    allocate(theParticles(numParticles)) 
!  end if

  !=====================================================================
  ! We loop over all particles first
  !=====================================================================
#ifdef OUTPUT_LEVEL2 
  write(*,'(I3,A,I5,A)')myid, '> FBM Force Calculation for:', numParticles, 'remote particle(s)'
#endif

  call getAllRemoteParticles(theParticles)

  DO IP = 1,numParticles
#ifdef OUTPUT_LEVEL2 
  write(*,'(I3,A,I5,A,I5,A,I5)')myid,'>Rem Particle: #',theParticles(IP)%localIdx, ', of: ', size(theParticles), ' shortId', theParticles(IP)%bytes(1) + 1
#endif

  particleId = theParticles(IP)%localIdx
  particleId = theParticles(IP)%bytes(1) + 1
  
  theParticles(IP)%uniqueIdx = 0

  Center = theParticles(IP)%position 
  
  DResForceX = 0D0
  DResForceY = 0D0
  DResForceZ = 0D0
  DTrqForceX = 0d0
  DTrqForceY = 0d0
  DTrqForceZ = 0d0
  
  nnel = 0
  !=====================================================================
  ! loop over all elements 
  !=====================================================================
  DO IEL=1,NEL
  
  CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
  IF (IER.LT.0)then
    call ExitError('Error in GetForces',1728)
  end if 
  
   NJALFA=0
   NIALFA=0
   DO I=1,IDFL
     IG=KDFG(I)

     crit1 = .false.
     crit2 = .false.
     IF((ALPHA(IG).EQ.0).or.(.not. Alpha(ig) .eq. particleId))THEN
       crit1 = .true.
     end if

     ! The preferred method of checking is to:
     ! compare the particle system ID (which is a 64bit unsigned integer) to
     ! the fortran long representation of the FictKNPR array
     IF((ALPHA(IG).EQ.0).or.(.not. longIdMatch(IG, theParticles(IP)%bytes)))THEN
!     IF((ALPHA(IG).EQ.0).or.(.not. map_local_to_system(particleId,IG)))THEN
      crit2 = .true.
      NJALFA=NJALFA+1
     ENDIF

 !    ! This condition represents a serious error
 !    if(.not. (crit1 .eqv. crit2))then
 !      write(*,*) 'Rem: myid:',myid,'> a(ig)=',alpha(ig), ' partId=',particleId, ' partLongId=',theParticles(ip)%bytes,' uintFict',&
 !        FictKNPR_uint64(ig)%bytes, ' F=',FictKNPR(ig)
 !!      call ExitError('Error in GetForces',157)
 !    end if


     crit1 = .false.
     crit2 = .false.

     IF (Alpha(ig) .eq. particleId) THEN
       crit1 = .true.
     end if

     IF (longIdMatch(IG, theParticles(IP)%bytes)) THEN
      crit2 = .true.
      NIALFA=NIALFA+1
     ENDIF

!     ! This condition represents a serious error
!     if(.not. (crit1 .eqv. crit2))then
!       write(*,*) 'Rem: myid:',myid,'> a(ig)=',alpha(ig), ' partId=',particleId, ' partLongId=',theParticles(ip)%bytes,' uintFict',&
!         FictKNPR_uint64(ig)%bytes, ' F=',FictKNPR(ig)
!!       call ExitError('Error in GetForces',174)
!     end if

   ENDDO


  ! Skip elements where all dofs are inside or
  ! all dofs are outside
  IF(NJALFA.EQ.27.OR.NIALFA.EQ.27) cycle

  nnel = nnel + 1
  DNY = DVISC(IEL)

  ! *** Evaluation of coordinates of the vertices
  DX0 = 0d0
  DY0 = 0d0
  DZ0 = 0d0

  ! Calculate the element center in d0
  ! and get the element vertices in dx
  DO IVE=1,NVE
    JP=KVERT(IVE,IEL)
    KVE(IVE)=JP
    DX(IVE)=DCORVG(1,JP)
    DY(IVE)=DCORVG(2,JP)
    DZ(IVE)=DCORVG(3,JP)
    DX0 = DX0 + 0.125d0*DX(IVE)
    DY0 = DY0 + 0.125d0*DY(IVE)
    DZ0 = DZ0 + 0.125d0*DZ(IVE)
  end do

!==============================================================    
!                  Calculation of the jacobian
!==============================================================    
!     Precompute averaged values for later use 
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

  ! Initialize the ELE
  CALL ELE(0D0,0D0,0D0,-2)

  IF (IER.LT.0)then
    call ExitError('Error in GetForces',1800)
  end if 

    pressSum = 0.0
    momSum = 0.0
    ! Loop over all cubature points
    DO ICUBP=1,NCUBP
  
    ! Get the coordinate of the cubature point
    XI1=DXI(ICUBP,1)
    XI2=DXI(ICUBP,2)
    XI3=DXI(ICUBP,3)

    ! Jacobian of the bilinear mapping onto the reference element
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
  
    ! Now calculate the real coordinates of the cubature point
    XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
    YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
    ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2
  
    ! Call ELE with the cubature point to SET the 
    ! derivatives in the DBAS-Array
    CALL ELE(XI1,XI2,XI3,-3)

    IF (IER.LT.0)then
      call ExitError('Error in GetForces',1837)
    end if 

    ! Evaluate the solution values and derivatives in the cubature point
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
      DBI1=DBAS(1,KDFL(I),1)
      DBI2=DBAS(1,KDFL(I),2)
      DBI3=DBAS(1,KDFL(I),3)
      DBI4=DBAS(1,KDFL(I),4)
  
      !------FOR U1--------
      DU1V=DU1V+U1(IG)*DBI1
      DU1X=DU1X+U1(IG)*DBI2
      DU1Y=DU1Y+U1(IG)*DBI3
      DU1Z=DU1Z+U1(IG)*DBI4
  
      !------FOR U2--------
      DU2V=DU2V+U2(IG)*DBI1
      DU2X=DU2X+U2(IG)*DBI2
      DU2Y=DU2Y+U2(IG)*DBI3
      DU2Z=DU2Z+U2(IG)*DBI4
  
      !------FOR U3--------
      DU3V=DU3V+U3(IG)*DBI1
      DU3X=DU3X+U3(IG)*DBI2
      DU3Y=DU3Y+U3(IG)*DBI3
      DU3Z=DU3Z+U3(IG)*DBI4
  
      !------FOR ALFA------
      !IF((map_local_to_system(particleId,IG)))THEN
      IF (longIdMatch(IG, theParticles(IP)%bytes)) THEN
      !IF (Alpha(ig) .eq. particleId) THEN
       DALPHA = 1d0
      ELSE
       DALPHA = 0d0
      END IF
      DALV=DALV+DALPHA*DBI1
      DALX=DALX+DALPHA*DBI2
      DALY=DALY+DALPHA*DBI3
      DALZ=DALZ+DALPHA*DBI4
    ENDDO
  
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
    momSum = momSum + (DU3Z)
    pressSum = pressSum + press 
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
    end do ! end loop cubature points

  end do ! end loop elements

#ifdef SED_BENCH
  DResForceX = 0.0
  DResForceY = 0.0
  DResForceZ = 4.0 * DResForceZ
  theParticles(ip)%force(:) = (/DResForceX, DResForceY, DResForceZ/)
  theParticles(ip)%torque(:) = (/0.0, 0.0, 0.0/)
#else
  theParticles(ip)%force(:) = (/DResForceX, DResForceY, DResForceZ/)
  theParticles(ip)%torque(:) = (/DTrqForceX, DTrqForceY, DTrqForceZ/)
#endif

  theNorm = calculate_l2_norm(DResForceX, DResForceY, DResForceZ)
  if( theNorm >= localMax ) then
    localMax = theNorm
  end if

#ifdef OUTPUT_LEVEL2 
  write(*,'(I3,A,I3,A,3D12.4,A,I10)')myid,' pidx=', theParticles(ip)%bytes(1) + 1, ' theRemoteForce(x,y,z): ', (/DResForceX, DResForceY, DResForceZ/), ' elems: ', nnel
#endif
!  write(*,*)myid,' nnel: ',nnel 
  !write(*,*)'TheTorque: ', theParticles(ip)%torque

  ! This function is in the dem_query module
  call setRemoteForcesMapped(theParticles(ip))

  END DO ! nParticles

!  call MPI_Reduce(totalMax, maxForce, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_EX0, ierr)
!  if (myid .eq. 1)then
!    write(*,*)'maxRemoteFluidForce: ', maxForce 
!    write(*,*)'maxFluidForce2: ', maxLocal 
!  end if

END SUBROUTINE ForcesRemoteParticles
!************************************************************************
!
!************************************************************************
SUBROUTINE GetForcesFC2(factors,U1,U2,U3,P,ALPHA,DVISC,KVERT,KAREA,KEDGE,&
                        DCORVG,ELE)
!************************************************************************
!*     Discrete convection operator: Q1~ elements (nonparametric)
!*-----------------------------------------------------------------------
USE PP3D_MPI, ONLY:myid,showID,COMM_SUMMN,Reduce_myMPI 
USE var_QuadScalar, ONLY : myExport,Properties
use cinterface
use dem_query


IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
CHARACTER SUB*6,FMT*15,CPARAM*120

PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,&
           NNDIM=3,NNCOF=10)
PARAMETER (Q2=0.5D0,Q8=0.125D0)

! An array of scaling factors
Real*8 :: factors(*)

! U/V/W velocity components
REAL*8 :: U1(*),U2(*),U3(*)

! Pressure, viscosity
Real*8 :: P(*),DVISC(*)

! Coordinates
Real*8 :: DCORVG(NNDIM,*)

! Alpha function
integer :: ALPHA(*)
INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
INTEGER KDFG(NNBAS),KDFL(NNBAS)

REAL*8 :: DResForceX,DResForceY,DResForceZ, localMax, localMaxRemote, totalMax
REAL*8 :: DTrqForceX,DTrqForceY,DTrqForceZ
REAL*8 :: Center(3),dForce(6)

type(tParticleData), dimension(:), allocatable :: theParticles
integer :: numParticles, particleId, ierr

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

!user COMMON blocks
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

  localMax       = 0.0
  localMaxRemote = 0.0
  totalMax       = 0.0

  IF (myid /= 0) then

  ! Local particles first
  call ForcesLocalParticles(factors,U1,U2,U3,P,ALPHA,&
                            DVISC,KVERT,KAREA,KEDGE,&
                            DCORVG,ELE, localMax)

  ! Remote particles second
  call ForcesRemoteParticles(factors,U1,U2,U3,P,ALPHA,&
                             DVISC,KVERT,KAREA,KEDGE,&
                             DCORVG,ELE, localMaxRemote)

  ! Now we synchronize the forces
  call sync_forces()
  call debug_output_force() 
  end if

  localMax = MAX(localMax, localMaxRemote)

!  call MPI_Reduce(totalMax, localMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_ALL, ierr)
!  call MPI_Barrier(MPI_COMM_ALL, ierr)
  call Reduce_myMPI(localMax, totalMax)

!  IF (myid == 0)then
!    write(*,*)'totalMaxFluidForce: ', totalMax 
!  end if
END
