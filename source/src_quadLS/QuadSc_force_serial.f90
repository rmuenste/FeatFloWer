!************************************************************************
!
! QuadSc_force_serial.f90
!
! Serial PE mode implementation of ForcesLocalParticles
! This file provides a clean, focused implementation for math students
! and researchers to understand the core fictitious boundary method.
!
! In serial PE mode:
! - Each CFD domain maintains full particle information
! - All particles are "local" (no remote/shadow copies)
! - Force synchronization happens at CFD layer via MPI
!
! Structure:
!   ForcesLocalParticlesSerial_Standard - Brute-force (all elements)
!   ForcesLocalParticlesSerial_KVEL     - KVEL-accelerated (candidate elements)
!   ForcesLocalParticlesSerial          - Wrapper (same external interface)
!
!************************************************************************

!========================================================================
! Standard (brute-force) force computation
! Loops over ALL elements for each particle.
!========================================================================
SUBROUTINE ForcesLocalParticlesSerial_Standard(factors,U1,U2,U3,P,ALPHA,DVISC,KVERT,KAREA,KEDGE,&
                                DCORVG,ELE, theParticles, numParticles, &
                                forceArray_out, boundaryElems_out, nBoundaryElems_out, maxBdryPerPart)
USE PP3D_MPI, ONLY:myid,showID,COMM_SUMMN
USE var_QuadScalar, ONLY : myExport,Properties,FictKNPR_uint64, FictKNPR
USE var_QuadScalar, ONLY : AlphaRelax
use cinterface
use dem_query
IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
CHARACTER SUB*6,FMT*15,CPARAM*120

PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,&
           NNDIM=3,NNCOF=10)
PARAMETER (Q2=0.5D0,Q8=0.125D0)

! Input arrays (same as original)
Real*8 :: factors(*)
REAL*8 :: U1(*),U2(*),U3(*)
Real*8 :: P(*),DVISC(*)
Real*8 :: DCORVG(NNDIM,*)
integer :: ALPHA(*)
INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
EXTERNAL ELE

! Particle data (passed from wrapper)
integer, intent(in) :: numParticles
type(tParticleData), intent(in) :: theParticles(numParticles)

! Output arrays
integer, intent(in) :: maxBdryPerPart
real*8, intent(out) :: forceArray_out(6*numParticles)
integer, intent(out) :: boundaryElems_out(maxBdryPerPart, numParticles)
integer, intent(out) :: nBoundaryElems_out(numParticles)

! Local variables
INTEGER KDFG(NNBAS),KDFL(NNBAS)
REAL*8 :: DResForceX,DResForceY,DResForceZ
REAL*8 :: DTrqForceX,DTrqForceY,DTrqForceZ
REAL*8 :: Center(3), omega(3)
integer :: resi
real*8 :: sliX
real*8, dimension(:), allocatable :: forceArray
integer :: iPointer

COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
COMMON /ERRCTL/ IER,ICHECK
COMMON /CHAR/   SUB,FMT(3),CPARAM
COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,&
                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),&
                IEL,NDIM
COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,&
                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,ICYCLE,KPRSM,KPOSM
COMMON /COAUX1/ KDFG,KDFL,IDFL

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

! Initialize outputs
nBoundaryElems_out = 0
boundaryElems_out = 0
omega = (/0.0, 0.0, 0.0/)

allocate(forceArray(6*numParticles))
forceArray = 0.0d0

IF (myid /= 0) THEN

  !========================================================================
  ! Element setup (workers only)
  !========================================================================
  DO I= 1,NNDER
    BDER(I)=.false.
  end do
  DO I=1,4
    BDER(I)=.true.
  end do

  IELTYP=-1
  CALL ELE(0D0,0D0,0D0,IELTYP)
  IDFL=NDFL(IELTYP)

  ICUB=9
  CALL CB3H(ICUB)
  IF (IER /= 0)then
    call ExitError('Error in ForcesLocalParticlesSerial_Standard',1702)
  end if

  !========================================================================
  ! MAIN PARTICLE LOOP
  !========================================================================
  DO IP = 1,numParticles

  Center = theParticles(IP)%position

  resi = 0
  sliX = 0.0
#ifdef ENABLE_LUBRICATION
  call sliding_wall_force(center,theParticles(IP)%velocity, omega, resi, sliX)
#endif

  DResForceX = 0D0
  DResForceY = 0D0
  DResForceZ = 0D0
  DTrqForceX = 0d0
  DTrqForceY = 0d0
  DTrqForceZ = 0d0

  !========================================================================
  ! ELEMENT LOOP - ALL elements (brute-force)
  !========================================================================
  DO IEL = 1, NEL

  CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
  IF (IER < 0)then
    call ExitError('Error in ForcesLocalParticlesSerial_Standard',1728)
  end if

   !----------------------------------------------------------------------
   ! Check if element intersects particle boundary
   !----------------------------------------------------------------------
   NJALFA=0
   NIALFA=0
   DO I=1,IDFL
     IG=KDFG(I)
     IF((ALPHA(IG) == 0).or.(.not. longIdMatch(IG, theParticles(IP)%bytes)))THEN
      NJALFA=NJALFA+1
     ENDIF
     IF (longIdMatch(IG, theParticles(IP)%bytes)) THEN
      NIALFA=NIALFA+1
     ENDIF
   ENDDO

  ! Skip elements where all dofs are inside or all dofs are outside
  IF(NJALFA == 27.OR.NIALFA == 27) cycle

  ! Track boundary element
  nBoundaryElems_out(IP) = nBoundaryElems_out(IP) + 1
  if (nBoundaryElems_out(IP) <= maxBdryPerPart) then
    boundaryElems_out(nBoundaryElems_out(IP), IP) = IEL
  end if

  DNY = DVISC(IEL)

  !----------------------------------------------------------------------
  ! Element vertex coordinates
  !----------------------------------------------------------------------
  DX0 = 0d0
  DY0 = 0d0
  DZ0 = 0d0

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

  !======================================================================
  ! Jacobian precomputation (for trilinear mapping)
  !======================================================================
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

  CALL ELE(0D0,0D0,0D0,-2)
  IF (IER < 0)then
    call ExitError('Error in ForcesLocalParticlesSerial_Standard',1800)
  end if

    !====================================================================
    ! CUBATURE POINT LOOP
    !====================================================================
    DO ICUBP=1,NCUBP

    XI1=DXI(ICUBP,1)
    XI2=DXI(ICUBP,2)
    XI3=DXI(ICUBP,3)

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

    XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
    YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
    ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2

    CALL ELE(XI1,XI2,XI3,-3)
    IF (IER < 0)then
      call ExitError('Error in ForcesLocalParticlesSerial_Standard',1837)
    end if

    DU1V=0D0; DU1X=0D0; DU1Y=0D0; DU1Z=0D0
    DU2V=0D0; DU2X=0D0; DU2Y=0D0; DU2Z=0D0
    DU3V=0D0; DU3X=0D0; DU3Y=0D0; DU3Z=0D0
    DALV=0D0; DALX=0D0; DALY=0D0; DALZ=0D0

    DO I=1,IDFL
      IG=KDFG(I)
      DBI1=DBAS(1,KDFL(I),1)
      DBI2=DBAS(1,KDFL(I),2)
      DBI3=DBAS(1,KDFL(I),3)
      DBI4=DBAS(1,KDFL(I),4)

      DU1V=DU1V+U1(IG)*DBI1
      DU1X=DU1X+U1(IG)*DBI2
      DU1Y=DU1Y+U1(IG)*DBI3
      DU1Z=DU1Z+U1(IG)*DBI4

      DU2V=DU2V+U2(IG)*DBI1
      DU2X=DU2X+U2(IG)*DBI2
      DU2Y=DU2Y+U2(IG)*DBI3
      DU2Z=DU2Z+U2(IG)*DBI4

      DU3V=DU3V+U3(IG)*DBI1
      DU3X=DU3X+U3(IG)*DBI2
      DU3Y=DU3Y+U3(IG)*DBI3
      DU3Z=DU3Z+U3(IG)*DBI4

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

    JJ = 4*(IEL-1) + 1
    Press =          P(JJ  ) + (XX-DX0)*P(JJ+1) + &
            (YY-DY0)*P(JJ+2) + (ZZ-DZ0)*P(JJ+3)

    DN1=-DALX
    DN2=-DALY
    DN3=-DALZ

    AH1=-Press*DN1+DNY*((DU1X+DU1X)*DN1+(DU1Y+DU2X)*DN2 +&
    (DU1Z+DU3X)*DN3)
    AH2=-Press*DN2+DNY*((DU2X+DU1Y)*DN1+(DU2Y+DU2Y)*DN2 +&
    (DU2Z+DU3Y)*DN3)
    AH3=-Press*DN3+DNY*((DU3X+DU1Z)*DN1+(DU3Y+DU2Z)*DN2 +&
    (DU3Z+DU3Z)*DN3)

    DResForceX = DResForceX + AH1*OM
    DResForceY = DResForceY + AH2*OM
    DResForceZ = DResForceZ + AH3*OM

    XTORQUE = XX - Center(1)
    YTORQUE = YY - Center(2)
    ZTORQUE = ZZ - Center(3)
    ATQX = YTORQUE*AH3 - ZTORQUE*AH2
    ATQY = ZTORQUE*AH1 - XTORQUE*AH3
    ATQZ = XTORQUE*AH2 - YTORQUE*AH1
    DTrqForceX = DTrqForceX + ATQX*OM
    DTrqForceY = DTrqForceY + ATQY*OM
    DTrqForceZ = DTrqForceZ + ATQZ*OM

    end do ! cubature points

  end do ! elements

  !========================================================================
  ! Post-processing: benchmark-specific modifications
  !========================================================================
#ifdef SED_BENCH
  DResForceX = 0.0
  DResForceY = 0.0
  DResForceZ = 4.0 * DResForceZ
  DTrqForceX = 0.0
  DTrqForceY = 0.0
  DTrqForceZ = 0.0
#else
#ifdef ENABLE_LUBRICATION
  if (resi > 0) then
    DResForceX = DResForceX + AlphaRelax * sliX
  end if
#endif
#endif

  ! Pack into forceArray
  iPointer = 6*(IP-1)
  forceArray(iPointer+1) = DResForceX
  forceArray(iPointer+2) = DResForceY
  forceArray(iPointer+3) = DResForceZ
  forceArray(iPointer+4) = DTrqForceX
  forceArray(iPointer+5) = DTrqForceY
  forceArray(iPointer+6) = DTrqForceZ

  END DO ! particles

END IF ! myid /= 0

! MPI force summation (ALL ranks)
CALL COMM_SUMMN(forceArray, 6*numParticles)
forceArray_out = forceArray
deallocate(forceArray)

END SUBROUTINE ForcesLocalParticlesSerial_Standard


!========================================================================
! KVEL-accelerated force computation
! Uses ParticleVertexCache + KVEL/KEEL/KAAL to find candidate elements.
!========================================================================
SUBROUTINE ForcesLocalParticlesSerial_KVEL(factors,U1,U2,U3,P,ALPHA,DVISC,KVERT,KAREA,KEDGE,&
                                DCORVG,ELE, theParticles, numParticles, &
                                forceArray_out, boundaryElems_out, nBoundaryElems_out, maxBdryPerPart)
USE PP3D_MPI, ONLY:myid,showID,COMM_SUMMN
USE var_QuadScalar, ONLY : myExport,Properties,FictKNPR_uint64, FictKNPR
USE var_QuadScalar, ONLY : AlphaRelax
USE var_QuadScalar, ONLY : ParticleVertexCache, bUseKVEL_Accel, myKVEL_Stats, mg_mesh
use cinterface
use dem_query
IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
CHARACTER SUB*6,FMT*15,CPARAM*120

PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,&
           NNDIM=3,NNCOF=10)
PARAMETER (Q2=0.5D0,Q8=0.125D0)

! Input arrays (same as original)
Real*8 :: factors(*)
REAL*8 :: U1(*),U2(*),U3(*)
Real*8 :: P(*),DVISC(*)
Real*8 :: DCORVG(NNDIM,*)
integer :: ALPHA(*)
INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
EXTERNAL ELE

! Particle data (passed from wrapper)
integer, intent(in) :: numParticles
type(tParticleData), intent(in) :: theParticles(numParticles)

! Output arrays
integer, intent(in) :: maxBdryPerPart
real*8, intent(out) :: forceArray_out(6*numParticles)
integer, intent(out) :: boundaryElems_out(maxBdryPerPart, numParticles)
integer, intent(out) :: nBoundaryElems_out(numParticles)

! Local variables
INTEGER KDFG(NNBAS),KDFL(NNBAS)
REAL*8 :: DResForceX,DResForceY,DResForceZ
REAL*8 :: DTrqForceX,DTrqForceY,DTrqForceZ
REAL*8 :: Center(3), omega(3)
integer :: resi
real*8 :: sliX
real*8, dimension(:), allocatable :: forceArray
integer :: iPointer
real*8 :: dbuf1(2)

! KVEL acceleration variables
LOGICAL, ALLOCATABLE :: bCandidateElement(:)
INTEGER, ALLOCATABLE :: CandidateList(:)
INTEGER :: nCandidates, iCand, iVtx, j, iedge, iface

COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
COMMON /ERRCTL/ IER,ICHECK
COMMON /CHAR/   SUB,FMT(3),CPARAM
COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,&
                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),&
                IEL,NDIM
COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,&
                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,ICYCLE,KPRSM,KPOSM
COMMON /COAUX1/ KDFG,KDFL,IDFL

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

! Initialize outputs
nBoundaryElems_out = 0
boundaryElems_out = 0
omega = (/0.0, 0.0, 0.0/)

! Initialize KVEL statistics
myKVEL_Stats%nCachedVertices = 0
myKVEL_Stats%nCandidateElements = 0
myKVEL_Stats%nBoundaryElements = 0

allocate(forceArray(6*numParticles))
forceArray = 0.0d0

IF (myid /= 0) THEN

  !========================================================================
  ! Element setup (workers only)
  !========================================================================
  DO I= 1,NNDER
    BDER(I)=.false.
  end do
  DO I=1,4
    BDER(I)=.true.
  end do

  IELTYP=-1
  CALL ELE(0D0,0D0,0D0,IELTYP)
  IDFL=NDFL(IELTYP)

  ICUB=9
  CALL CB3H(ICUB)
  IF (IER /= 0)then
    call ExitError('Error in ForcesLocalParticlesSerial_KVEL',1702)
  end if

  !========================================================================
  ! MAIN PARTICLE LOOP
  !========================================================================
  DO IP = 1,numParticles

  Center = theParticles(IP)%position

  resi = 0
  sliX = 0.0
#ifdef ENABLE_LUBRICATION
  call sliding_wall_force(center,theParticles(IP)%velocity, omega, resi, sliX)
#endif

  DResForceX = 0D0
  DResForceY = 0D0
  DResForceZ = 0D0
  DTrqForceX = 0d0
  DTrqForceY = 0d0
  DTrqForceZ = 0d0

  !========================================================================
  ! Build KVEL Candidate Element Set
  !========================================================================
  allocate(bCandidateElement(NEL))
  allocate(CandidateList(NEL))

  bCandidateElement = .FALSE.
  nCandidates = 0

#ifdef ENABLE_FBM_ACCELERATION
  if (bUseKVEL_Accel .and. allocated(ParticleVertexCache)) then
#else
  if (.FALSE.) then
#endif
    if (IP <= size(ParticleVertexCache)) then
      if (ParticleVertexCache(IP)%nVertices > 0) then

        DO iVtx = 1, ParticleVertexCache(IP)%nVertices
          ivt = ParticleVertexCache(IP)%dofIndices(iVtx)

          IF (ivt <= NVT) THEN
            DO j = 1, mg_mesh%level(NLMAX)%nvel
              IEL = mg_mesh%level(NLMAX)%kvel(j, ivt)
              IF (IEL == 0) EXIT
              IF (.not. bCandidateElement(IEL)) THEN
                bCandidateElement(IEL) = .TRUE.
                nCandidates = nCandidates + 1
                CandidateList(nCandidates) = IEL
              END IF
            END DO

          ELSE IF (ivt <= NVT + NET) THEN
            iedge = ivt - NVT
            DO j = 1, mg_mesh%level(NLMAX)%neel
              IEL = mg_mesh%level(NLMAX)%keel(j, iedge)
              IF (IEL == 0) EXIT
              IF (.not. bCandidateElement(IEL)) THEN
                bCandidateElement(IEL) = .TRUE.
                nCandidates = nCandidates + 1
                CandidateList(nCandidates) = IEL
              END IF
            END DO

          ELSE IF (ivt <= NVT + NET + NAT) THEN
            iface = ivt - NVT - NET
            DO j = 1, mg_mesh%level(NLMAX)%naal
              IEL = mg_mesh%level(NLMAX)%kaal(j, iface)
              IF (IEL == 0) EXIT
              IF (.not. bCandidateElement(IEL)) THEN
                bCandidateElement(IEL) = .TRUE.
                nCandidates = nCandidates + 1
                CandidateList(nCandidates) = IEL
              END IF
            END DO
          END IF
        END DO

      end if
    end if
  end if

  ! Handle case when no candidates found
  if (nCandidates == 0) then
    if (bUseKVEL_Accel .and. allocated(ParticleVertexCache)) then
      ! Cache is populated but particle has no DOFs on this rank
      deallocate(bCandidateElement)
      deallocate(CandidateList)
      ! Pack zero forces for this particle
      iPointer = 6*(IP-1)
      forceArray(iPointer+1) = 0.0d0
      forceArray(iPointer+2) = 0.0d0
      forceArray(iPointer+3) = 0.0d0
      forceArray(iPointer+4) = 0.0d0
      forceArray(iPointer+5) = 0.0d0
      forceArray(iPointer+6) = 0.0d0
      cycle  ! next particle
    else
      ! No cache available, fall back to all elements
      nCandidates = NEL
      DO IEL = 1, NEL
        CandidateList(IEL) = IEL
      END DO
    end if
  end if

  ! Accumulate statistics
  myKVEL_Stats%nCandidateElements = myKVEL_Stats%nCandidateElements + nCandidates

  !========================================================================
  ! ELEMENT LOOP - Process CANDIDATE elements only
  !========================================================================
  DO iCand = 1, nCandidates
  IEL = CandidateList(iCand)

  CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
  IF (IER < 0)then
    call ExitError('Error in ForcesLocalParticlesSerial_KVEL',1728)
  end if

   NJALFA=0
   NIALFA=0
   DO I=1,IDFL
     IG=KDFG(I)
     IF((ALPHA(IG) == 0).or.(.not. longIdMatch(IG, theParticles(IP)%bytes)))THEN
      NJALFA=NJALFA+1
     ENDIF
     IF (longIdMatch(IG, theParticles(IP)%bytes)) THEN
      NIALFA=NIALFA+1
     ENDIF
   ENDDO

  IF(NJALFA == 27.OR.NIALFA == 27) cycle

  ! Track boundary element
  nBoundaryElems_out(IP) = nBoundaryElems_out(IP) + 1
  if (nBoundaryElems_out(IP) <= maxBdryPerPart) then
    boundaryElems_out(nBoundaryElems_out(IP), IP) = IEL
  end if

  DNY = DVISC(IEL)

  DX0 = 0d0
  DY0 = 0d0
  DZ0 = 0d0

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

  CALL ELE(0D0,0D0,0D0,-2)
  IF (IER < 0)then
    call ExitError('Error in ForcesLocalParticlesSerial_KVEL',1800)
  end if

    DO ICUBP=1,NCUBP

    XI1=DXI(ICUBP,1)
    XI2=DXI(ICUBP,2)
    XI3=DXI(ICUBP,3)

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

    XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
    YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
    ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2

    CALL ELE(XI1,XI2,XI3,-3)
    IF (IER < 0)then
      call ExitError('Error in ForcesLocalParticlesSerial_KVEL',1837)
    end if

    DU1V=0D0; DU1X=0D0; DU1Y=0D0; DU1Z=0D0
    DU2V=0D0; DU2X=0D0; DU2Y=0D0; DU2Z=0D0
    DU3V=0D0; DU3X=0D0; DU3Y=0D0; DU3Z=0D0
    DALV=0D0; DALX=0D0; DALY=0D0; DALZ=0D0

    DO I=1,IDFL
      IG=KDFG(I)
      DBI1=DBAS(1,KDFL(I),1)
      DBI2=DBAS(1,KDFL(I),2)
      DBI3=DBAS(1,KDFL(I),3)
      DBI4=DBAS(1,KDFL(I),4)

      DU1V=DU1V+U1(IG)*DBI1
      DU1X=DU1X+U1(IG)*DBI2
      DU1Y=DU1Y+U1(IG)*DBI3
      DU1Z=DU1Z+U1(IG)*DBI4

      DU2V=DU2V+U2(IG)*DBI1
      DU2X=DU2X+U2(IG)*DBI2
      DU2Y=DU2Y+U2(IG)*DBI3
      DU2Z=DU2Z+U2(IG)*DBI4

      DU3V=DU3V+U3(IG)*DBI1
      DU3X=DU3X+U3(IG)*DBI2
      DU3Y=DU3Y+U3(IG)*DBI3
      DU3Z=DU3Z+U3(IG)*DBI4

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

    JJ = 4*(IEL-1) + 1
    Press =          P(JJ  ) + (XX-DX0)*P(JJ+1) + &
            (YY-DY0)*P(JJ+2) + (ZZ-DZ0)*P(JJ+3)

    DN1=-DALX
    DN2=-DALY
    DN3=-DALZ

    AH1=-Press*DN1+DNY*((DU1X+DU1X)*DN1+(DU1Y+DU2X)*DN2 +&
    (DU1Z+DU3X)*DN3)
    AH2=-Press*DN2+DNY*((DU2X+DU1Y)*DN1+(DU2Y+DU2Y)*DN2 +&
    (DU2Z+DU3Y)*DN3)
    AH3=-Press*DN3+DNY*((DU3X+DU1Z)*DN1+(DU3Y+DU2Z)*DN2 +&
    (DU3Z+DU3Z)*DN3)

    DResForceX = DResForceX + AH1*OM
    DResForceY = DResForceY + AH2*OM
    DResForceZ = DResForceZ + AH3*OM

    XTORQUE = XX - Center(1)
    YTORQUE = YY - Center(2)
    ZTORQUE = ZZ - Center(3)
    ATQX = YTORQUE*AH3 - ZTORQUE*AH2
    ATQY = ZTORQUE*AH1 - XTORQUE*AH3
    ATQZ = XTORQUE*AH2 - YTORQUE*AH1
    DTrqForceX = DTrqForceX + ATQX*OM
    DTrqForceY = DTrqForceY + ATQY*OM
    DTrqForceZ = DTrqForceZ + ATQZ*OM

    end do ! cubature points

  end do ! candidate elements

  ! Cleanup candidate lists
  deallocate(bCandidateElement)
  deallocate(CandidateList)

  !========================================================================
  ! Post-processing
  !========================================================================
#ifdef SED_BENCH
  DResForceX = 0.0
  DResForceY = 0.0
  DResForceZ = 4.0 * DResForceZ
  DTrqForceX = 0.0
  DTrqForceY = 0.0
  DTrqForceZ = 0.0
#else
#ifdef ENABLE_LUBRICATION
  if (resi > 0) then
    DResForceX = DResForceX + AlphaRelax * sliX
  end if
#endif
#endif

  ! Pack into forceArray
  iPointer = 6*(IP-1)
  forceArray(iPointer+1) = DResForceX
  forceArray(iPointer+2) = DResForceY
  forceArray(iPointer+3) = DResForceZ
  forceArray(iPointer+4) = DTrqForceX
  forceArray(iPointer+5) = DTrqForceY
  forceArray(iPointer+6) = DTrqForceZ

  END DO ! particles

END IF ! myid /= 0

! Output KVEL acceleration statistics (aggregated across all ranks)
#ifdef ENABLE_FBM_ACCELERATION
if (bUseKVEL_Accel) then
#else
if (.FALSE.) then
#endif
  dbuf1(1) = dble(myKVEL_Stats%nCandidateElements)
  dbuf1(2) = dble(NEL) * dble(numParticles)
  call COMM_SUMMN(dbuf1, 2)
  if (myid == 0) then
    WRITE(*,'(A,ES12.4,A,ES12.4,A,F8.1,A)') &
      'KVEL: ', dbuf1(1), ' candidates vs ', &
      dbuf1(2), ' brute-force (', &
      dbuf1(2)/max(dbuf1(1),1.0d0), 'x speedup)'
  end if
end if

! MPI force summation (ALL ranks)
CALL COMM_SUMMN(forceArray, 6*numParticles)
forceArray_out = forceArray
deallocate(forceArray)

END SUBROUTINE ForcesLocalParticlesSerial_KVEL


!************************************************************************
! Wrapper: ForcesLocalParticlesSerial
! Keeps the EXACT SAME external interface as the original subroutine.
! In normal mode: calls _Standard or _KVEL based on flags.
! In debug mode (DEBUG_FBM_OPTIMIZATION): calls BOTH, compares, writes files.
!************************************************************************
SUBROUTINE ForcesLocalParticlesSerial(factors,U1,U2,U3,P,ALPHA,DVISC,KVERT,KAREA,KEDGE,&
                                DCORVG,ELE, maxLocal)
USE PP3D_MPI, ONLY:myid,showID,COMM_SUMMN,COMM_Maximumn
USE var_QuadScalar, ONLY : myExport,Properties,FictKNPR_uint64, FictKNPR, total_lubrication
USE var_QuadScalar, ONLY : AlphaRelax
USE var_QuadScalar, ONLY : ParticleVertexCache, bUseKVEL_Accel, myKVEL_Stats, mg_mesh
use cinterface
use dem_query
IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
CHARACTER SUB*6,FMT*15,CPARAM*120

PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,&
           NNDIM=3,NNCOF=10)
PARAMETER (Q2=0.5D0,Q8=0.125D0)

! Same signature as original
Real*8 :: factors(*)
REAL*8 :: U1(*),U2(*),U3(*)
Real*8 :: P(*),DVISC(*)
Real*8 :: DCORVG(NNDIM,*)
integer :: ALPHA(*)
REAL*8, intent(inout) :: maxLocal
INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
EXTERNAL ELE

! Local variables
type(tParticleData), dimension(:), allocatable :: theParticles
integer :: numParticles
real*8 :: dbuf1(2)
real*8 :: theNorm
integer :: iPointer
character(len=*), parameter :: fmt_sed = '(A,1X,A,ES14.6,1X,A,I6,1X,3ES14.6)'
real*8 :: time_out

! Force arrays and boundary element tracking
integer, parameter :: maxBdryPerPart = 500
real*8, dimension(:), allocatable :: forceArray_result
integer, dimension(:,:), allocatable :: boundaryElems_result
integer, dimension(:), allocatable :: nBoundaryElems_result

#ifdef DEBUG_FBM_OPTIMIZATION
real*8, dimension(:), allocatable :: forceArray_std, forceArray_kvel
integer, dimension(:,:), allocatable :: boundaryElems_std, boundaryElems_kvel
integer, dimension(:), allocatable :: nBoundaryElems_std, nBoundaryElems_kvel
real*8 :: maxDiff, diff_val
integer :: nMismatchedParticles, iDebugUnit
character(len=40) :: debugFileName
logical :: bSetsMatch
integer :: iElem, jElem
#endif

COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
COMMON /ERRCTL/ IER,ICHECK
COMMON /CHAR/   SUB,FMT(3),CPARAM
COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,&
                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,ICYCLE,KPRSM,KPOSM

SAVE

!========================================================================
! Initialization
!========================================================================
maxLocal = 0.0

if(myid == 1) write(*,*) 'Force calculation: SERIAL PE mode'

!========================================================================
! Get particle data (ALL RANKS need particle count for MPI)
!========================================================================
numParticles = numTotalParticles()
dbuf1(1) = dble(numParticles)
call COMM_Maximumn(dbuf1, 1)
numParticles = int(dbuf1(1))
if (numParticles == 0) then
  return
end if

! Allocate and fetch particle data
if(allocated(theParticles)) deallocate(theParticles)
allocate(theParticles(numParticles))

IF (myid /= 0) THEN
  call getAllParticles(theParticles)
#ifdef SED_BENCH
  if(myid == 1) write(*,*) 'Force with SED BENCH settings!'
#endif
END IF

!========================================================================
! Call force computation subroutine(s)
!========================================================================

#ifdef DEBUG_FBM_OPTIMIZATION
  !----------------------------------------------------------------------
  ! DEBUG MODE: Run BOTH methods and compare
  !----------------------------------------------------------------------
  allocate(forceArray_std(6*numParticles))
  allocate(boundaryElems_std(maxBdryPerPart, numParticles))
  allocate(nBoundaryElems_std(numParticles))

  allocate(forceArray_kvel(6*numParticles))
  allocate(boundaryElems_kvel(maxBdryPerPart, numParticles))
  allocate(nBoundaryElems_kvel(numParticles))

  ! Run standard (brute-force)
  if (myid == 1) write(*,'(A)') 'DEBUG_FBM: Running STANDARD (brute-force) force computation...'
  call ForcesLocalParticlesSerial_Standard(factors,U1,U2,U3,P,ALPHA,DVISC,KVERT,KAREA,KEDGE,&
       DCORVG,ELE, theParticles, numParticles, &
       forceArray_std, boundaryElems_std, nBoundaryElems_std, maxBdryPerPart)

  ! Run KVEL-accelerated
  if (myid == 1) write(*,'(A)') 'DEBUG_FBM: Running KVEL-accelerated force computation...'
  call ForcesLocalParticlesSerial_KVEL(factors,U1,U2,U3,P,ALPHA,DVISC,KVERT,KAREA,KEDGE,&
       DCORVG,ELE, theParticles, numParticles, &
       forceArray_kvel, boundaryElems_kvel, nBoundaryElems_kvel, maxBdryPerPart)

  !----------------------------------------------------------------------
  ! Compare forces and report
  !----------------------------------------------------------------------
  if (myid == 1) then
    write(*,'(A)') '=================================================='
    write(*,'(A)') 'DEBUG_FBM: Force Comparison (Standard vs KVEL)'
    write(*,'(A)') '=================================================='

    nMismatchedParticles = 0
    DO IP = 1, numParticles
      iPointer = 6*(IP-1)
      maxDiff = 0.0d0
      DO I = 1, 6
        diff_val = abs(forceArray_std(iPointer+I) - forceArray_kvel(iPointer+I))
        if (diff_val > maxDiff) maxDiff = diff_val
      END DO

      if (maxDiff > 1.0d-10) then
        nMismatchedParticles = nMismatchedParticles + 1
        write(*,'(A,I5,A,ES14.6)') '  Particle ', theParticles(IP)%bytes(1) + 1, &
          ' max force diff: ', maxDiff
        write(*,'(A,6ES14.6)') '    STD:  ', forceArray_std(iPointer+1:iPointer+6)
        write(*,'(A,6ES14.6)') '    KVEL: ', forceArray_kvel(iPointer+1:iPointer+6)
      end if

      ! Compare boundary element sets
      bSetsMatch = (nBoundaryElems_std(IP) == nBoundaryElems_kvel(IP))
      if (bSetsMatch) then
        ! Check if same elements (may be in different order)
        DO I = 1, min(nBoundaryElems_std(IP), maxBdryPerPart)
          iElem = boundaryElems_std(I, IP)
          bSetsMatch = .FALSE.
          DO J = 1, min(nBoundaryElems_kvel(IP), maxBdryPerPart)
            if (boundaryElems_kvel(J, IP) == iElem) then
              bSetsMatch = .TRUE.
              exit
            end if
          END DO
          if (.not. bSetsMatch) exit
        END DO
      end if
      if (.not. bSetsMatch) then
        write(*,'(A,I5,A,I5,A,I5)') '  Particle ', theParticles(IP)%bytes(1) + 1, &
          ' boundary elem mismatch: std=', nBoundaryElems_std(IP), &
          ' kvel=', nBoundaryElems_kvel(IP)
      end if
    END DO

    if (nMismatchedParticles == 0) then
      write(*,'(A)') 'DEBUG_FBM: Force comparison PASSED - all particles match'
    else
      write(*,'(A,I5,A)') 'DEBUG_FBM: Force comparison FAILED - ', &
        nMismatchedParticles, ' particles differ'
    end if
    write(*,'(A)') '=================================================='
  end if

  !----------------------------------------------------------------------
  ! Write per-timestep debug files (rank 1 only)
  !----------------------------------------------------------------------
  if (myid == 1) then
    ! Standard forces file
    write(debugFileName, '(A,I5.5,A)') 'debug_forces_standard_', ITNS, '.dat'
    iDebugUnit = 899
    open(unit=iDebugUnit, file=debugFileName, status='replace', action='write')
    DO IP = 1, numParticles
      iPointer = 6*(IP-1)
      write(iDebugUnit, '(I6,6ES20.12,I6)', advance='no') &
        theParticles(IP)%bytes(1) + 1, &
        forceArray_std(iPointer+1), forceArray_std(iPointer+2), &
        forceArray_std(iPointer+3), forceArray_std(iPointer+4), &
        forceArray_std(iPointer+5), forceArray_std(iPointer+6), &
        nBoundaryElems_std(IP)
      DO I = 1, min(nBoundaryElems_std(IP), maxBdryPerPart)
        write(iDebugUnit, '(I8)', advance='no') boundaryElems_std(I, IP)
      END DO
      write(iDebugUnit, *)
    END DO
    close(iDebugUnit)

    ! KVEL forces file
    write(debugFileName, '(A,I5.5,A)') 'debug_forces_kvel_', ITNS, '.dat'
    open(unit=iDebugUnit, file=debugFileName, status='replace', action='write')
    DO IP = 1, numParticles
      iPointer = 6*(IP-1)
      write(iDebugUnit, '(I6,6ES20.12,I6)', advance='no') &
        theParticles(IP)%bytes(1) + 1, &
        forceArray_kvel(iPointer+1), forceArray_kvel(iPointer+2), &
        forceArray_kvel(iPointer+3), forceArray_kvel(iPointer+4), &
        forceArray_kvel(iPointer+5), forceArray_kvel(iPointer+6), &
        nBoundaryElems_kvel(IP)
      DO I = 1, min(nBoundaryElems_kvel(IP), maxBdryPerPart)
        write(iDebugUnit, '(I8)', advance='no') boundaryElems_kvel(I, IP)
      END DO
      write(iDebugUnit, *)
    END DO
    close(iDebugUnit)

    write(*,'(A,I5.5,A)') 'DEBUG_FBM: Wrote debug files for timestep ', ITNS, &
      ' (debug_forces_standard_*.dat, debug_forces_kvel_*.dat)'
  end if

  ! Use standard (trusted) forces as the result
  allocate(forceArray_result(6*numParticles))
  forceArray_result = forceArray_std

  deallocate(forceArray_std)
  deallocate(boundaryElems_std)
  deallocate(nBoundaryElems_std)
  deallocate(forceArray_kvel)
  deallocate(boundaryElems_kvel)
  deallocate(nBoundaryElems_kvel)

#else
  !----------------------------------------------------------------------
  ! NORMAL MODE: Call one method
  !----------------------------------------------------------------------
  allocate(forceArray_result(6*numParticles))
  allocate(boundaryElems_result(maxBdryPerPart, numParticles))
  allocate(nBoundaryElems_result(numParticles))

#ifdef ENABLE_FBM_ACCELERATION
  if (bUseKVEL_Accel) then
    call ForcesLocalParticlesSerial_KVEL(factors,U1,U2,U3,P,ALPHA,DVISC,KVERT,KAREA,KEDGE,&
         DCORVG,ELE, theParticles, numParticles, &
         forceArray_result, boundaryElems_result, nBoundaryElems_result, maxBdryPerPart)
  else
    call ForcesLocalParticlesSerial_Standard(factors,U1,U2,U3,P,ALPHA,DVISC,KVERT,KAREA,KEDGE,&
         DCORVG,ELE, theParticles, numParticles, &
         forceArray_result, boundaryElems_result, nBoundaryElems_result, maxBdryPerPart)
  end if
#else
  call ForcesLocalParticlesSerial_Standard(factors,U1,U2,U3,P,ALPHA,DVISC,KVERT,KAREA,KEDGE,&
       DCORVG,ELE, theParticles, numParticles, &
       forceArray_result, boundaryElems_result, nBoundaryElems_result, maxBdryPerPart)
#endif

  deallocate(boundaryElems_result)
  deallocate(nBoundaryElems_result)
#endif

!========================================================================
! Unpack forces and call setForcesMapped
!========================================================================
if (myid /= 0) then
  DO IP = 1, numParticles
    iPointer = 6*(IP-1)
    theParticles(IP)%force(1) =  forceArray_result(iPointer+1)
    theParticles(IP)%force(2) =  forceArray_result(iPointer+2)
    theParticles(IP)%force(3) =  forceArray_result(iPointer+3)
    theParticles(IP)%torque(1) = forceArray_result(iPointer+4)
    theParticles(IP)%torque(2) = forceArray_result(iPointer+5)
    theParticles(IP)%torque(3) = forceArray_result(iPointer+6)

    theNorm = calculate_l2_norm(theParticles(IP)%force(1), &
              theParticles(IP)%force(2), theParticles(IP)%force(3))
    if (theNorm >= maxLocal) maxLocal = theNorm

    call setForcesMapped(theParticles(ip))
  END DO

#ifdef SED_BENCH
  if (myid == 1) then
    time_out = dble(itns - 1) * tstep
    DO IP = 1, numParticles
      write(*,fmt_sed) 'SED_BENCH_FORCE', 'time=', time_out, 'ip=', IP, &
        theParticles(IP)%force(:)
      write(*,fmt_sed) 'SED_BENCH_POS  ', 'time=', time_out, 'ip=', IP, &
        theParticles(IP)%position(:)
      write(*,fmt_sed) 'SED_BENCH_VEL  ', 'time=', time_out, 'ip=', IP, &
        theParticles(IP)%velocity(:)
    END DO
  end if
#endif

end if

deallocate(forceArray_result)

! Deallocate vertex cache (single-timestep lifetime)
if (allocated(ParticleVertexCache)) then
  DO IP = 1, size(ParticleVertexCache)
    if (allocated(ParticleVertexCache(IP)%dofIndices)) then
      deallocate(ParticleVertexCache(IP)%dofIndices)
    end if
  END DO
  deallocate(ParticleVertexCache)
end if

total_lubrication = 0.0d0

END SUBROUTINE ForcesLocalParticlesSerial
