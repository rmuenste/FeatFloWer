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
!************************************************************************
SUBROUTINE ForcesLocalParticlesSerial(factors,U1,U2,U3,P,ALPHA,DVISC,KVERT,KAREA,KEDGE,&
                                DCORVG,ELE, maxLocal)
!************************************************************************
!*     Fictitious Boundary Method force calculation: Q2/P1 elements
!*     Serial PE mode: handles all particles on each domain
!*-----------------------------------------------------------------------
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

! An array of scaling factors
Real*8 :: factors(*)

! U/V/W velocity components
REAL*8 :: U1(*),U2(*),U3(*)

! Pressure, viscosity
Real*8 :: P(*),DVISC(*)

! Coordinates
Real*8 :: DCORVG(NNDIM,*)

! Alpha function (particle indicator field)
integer :: ALPHA(*)

! The maximum force magnitude
REAL*8, intent(inout) :: maxLocal

INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
INTEGER KDFG(NNBAS),KDFL(NNBAS)

REAL*8 :: DResForceX,DResForceY,DResForceZ
REAL*8 :: DTrqForceX,DTrqForceY,DTrqForceZ
REAL*8 :: Center(3),dForce(6), omega(3)

type(tParticleData), dimension(:), allocatable :: theParticles
integer :: numParticles, resi, totalSliding
real*8 :: dbuf1(1)

real*8 :: theNorm, totalMax, localSliding, accumulatedSliding,sliX, localHydro
real*8, dimension(:), allocatable :: forceArray
integer :: iPointer
character(len=*), parameter :: fmt_sed = '(A,1X,A,ES14.6,1X,A,I6,1X,3ES14.6)'
real*8 :: time_out

! KVEL acceleration variables
LOGICAL, ALLOCATABLE :: bCandidateElement(:)
INTEGER, ALLOCATABLE :: CandidateList(:)
INTEGER :: nCandidates, iCand, iVtx, j

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

!========================================================================
! Initialization
!========================================================================
maxLocal = 0.0
omega = (/0.0, 0.0, 0.0/)

localSliding = 0.0
localHydro = 0.0
accumulatedSliding = 0.0

! Initialize KVEL statistics
myKVEL_Stats%nCachedVertices = 0
myKVEL_Stats%nCandidateElements = 0
myKVEL_Stats%nBoundaryElements = 0

!========================================================================
! Get particle data (ALL RANKS need particle count for MPI)
!========================================================================
if(myid == 1) write(*,*) 'Force calculation: SERIAL PE mode'

! Get total number of particles (ALL ranks need the same value).
! Rank 0 does not run PE and would see 0 here. Use a collective
! max across ranks so everyone agrees on the nonzero count.
numParticles = numTotalParticles()
dbuf1(1) = dble(numParticles)
call COMM_Maximumn(dbuf1, 1)
numParticles = int(dbuf1(1))
if (numParticles == 0) then
  return
end if

! Allocate theParticles (ALL ranks need this for unpacking after MPI)
if(allocated(theParticles)) deallocate(theParticles)
allocate(theParticles(numParticles))

IF (myid /= 0) THEN
  ! Only workers fetch particle data and compute forces
  call getAllParticles(theParticles)

#ifdef SED_BENCH
  if(myid == 1)then
    write(*,*) 'Force with SED BENCH settings!'
  end if
#endif

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
    call ExitError('Error in ForcesLocalParticlesSerial',1702)
  end if

#ifdef OUTPUT_LEVEL2
if (numParticles > 0) then
  write(*,'(A,I5,A,I3)')'FBM Force Calculation for ', numParticles, ' total particle(s) in dom ', myid
end if
#endif

!========================================================================
! MAIN PARTICLE LOOP
!========================================================================
totalSliding = 0
DO IP = 1,numParticles
#ifdef OUTPUT_LEVEL3
write(*,'(I3,A,I5)')myid,'>Particle: ',theParticles(IP)%bytes(1) + 1
#endif

Center = theParticles(IP)%position

resi = 0
sliX = 0.0
#ifdef ENABLE_LUBRICATION
call sliding_wall_force(center,theParticles(IP)%velocity, omega, resi, sliX)
#endif

localSliding = localSliding + sliX
totalSliding = totalSliding + resi

DResForceX = 0D0
DResForceY = 0D0
DResForceZ = 0D0
DTrqForceX = 0d0
DTrqForceY = 0d0
DTrqForceZ = 0d0

nnel = 0

!========================================================================
! Build KVEL Candidate Element Set
!========================================================================
allocate(bCandidateElement(NEL))
allocate(CandidateList(NEL))

bCandidateElement = .FALSE.
nCandidates = 0

! Build candidate set using vertex cache and KVEL
if (bUseKVEL_Accel .and. allocated(ParticleVertexCache)) then
  if (IP <= size(ParticleVertexCache)) then
    if (ParticleVertexCache(IP)%nVertices > 0) then

      DO iVtx = 1, ParticleVertexCache(IP)%nVertices
        ivt = ParticleVertexCache(IP)%dofIndices(iVtx)

        ! Only process corner vertices (Q2 vertex DOFs) for KVEL lookup
        IF (ivt <= NVT) THEN
          ! Use KVEL to find elements touching this vertex
          DO j = 1, mg_mesh%level(NLMAX)%nvel
            IEL = mg_mesh%level(NLMAX)%kvel(j, ivt)
            IF (IEL == 0) EXIT  ! No more elements

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

! Fallback if no candidates found
if (nCandidates == 0) then
  if (myid == 1) WRITE(*,'(A,I0,A)') 'WARNING: No KVEL candidates for particle ', IP, ' - using all elements'
  nCandidates = NEL
  DO IEL = 1, NEL
    CandidateList(IEL) = IEL
  END DO
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
  call ExitError('Error in ForcesLocalParticlesSerial',1728)
end if

 !----------------------------------------------------------------------
 ! Check if element intersects particle boundary
 !----------------------------------------------------------------------
 NJALFA=0
 NIALFA=0
 DO I=1,IDFL
   IG=KDFG(I)
   ! Check if DOF belongs to current particle using 64-bit ID matching
   IF((ALPHA(IG) == 0).or.(.not. longIdMatch(IG, theParticles(IP)%bytes)))THEN
    NJALFA=NJALFA+1
   ENDIF

   IF (longIdMatch(IG, theParticles(IP)%bytes)) THEN
    NIALFA=NIALFA+1
   ENDIF
 ENDDO

! Skip elements where all dofs are inside or all dofs are outside
IF(NJALFA == 27.OR.NIALFA == 27)then
  cycle
end if

nnel = nnel + 1
DNY = DVISC(IEL)

!----------------------------------------------------------------------
! Element vertex coordinates
!----------------------------------------------------------------------
DX0 = 0d0
DY0 = 0d0
DZ0 = 0d0

! Calculate the element center in d0 and get the element vertices in dx
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

! Initialize the element
CALL ELE(0D0,0D0,0D0,-2)

IF (IER < 0)then
  call ExitError('Error in ForcesLocalParticlesSerial',1800)
end if

  !====================================================================
  ! CUBATURE POINT LOOP - Core integration
  !====================================================================
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

  ! Call ELE with the cubature point to SET the derivatives in the DBAS-Array
  CALL ELE(XI1,XI2,XI3,-3)

  IF (IER < 0)then
    call ExitError('Error in ForcesLocalParticlesSerial',1837)
  end if

  !------------------------------------------------------------------
  ! Evaluate solution values and derivatives at cubature point
  !------------------------------------------------------------------
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

  !------------------------------------------------------------------
  ! Pressure interpolation (Q1~ discontinuous)
  !------------------------------------------------------------------
  JJ = 4*(IEL-1) + 1
  Press =          P(JJ  ) + (XX-DX0)*P(JJ+1) + &
          (YY-DY0)*P(JJ+2) + (ZZ-DZ0)*P(JJ+3)

  !------------------------------------------------------------------
  ! Normal vector from alpha field gradient
  !------------------------------------------------------------------
  DN1=-DALX
  DN2=-DALY
  DN3=-DALZ

  !------------------------------------------------------------------
  ! Stress tensor contribution: -p*I + 2*nu*D(u)
  ! where D(u) is the symmetric gradient (rate of deformation tensor)
  !------------------------------------------------------------------
  AH1=-Press*DN1+DNY*((DU1X+DU1X)*DN1+(DU1Y+DU2X)*DN2 +& ! full3D
  (DU1Z+DU3X)*DN3)
  AH2=-Press*DN2+DNY*((DU2X+DU1Y)*DN1+(DU2Y+DU2Y)*DN2 +& ! full3D
  (DU2Z+DU3Y)*DN3)
  AH3=-Press*DN3+DNY*((DU3X+DU1Z)*DN1+(DU3Y+DU2Z)*DN2 +& ! full3D
  (DU3Z+DU3Z)*DN3)

  !------------------------------------------------------------------
  ! Integrate force (surface integral over particle boundary)
  !------------------------------------------------------------------
  DResForceX = DResForceX + AH1*OM
  DResForceY = DResForceY + AH2*OM
  DResForceZ = DResForceZ + AH3*OM

  !------------------------------------------------------------------
  ! Integrate torque: r × F
  !------------------------------------------------------------------
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
  end do ! end loop cubature points

end do ! end loop elements

! Cleanup candidate lists
deallocate(bCandidateElement)
deallocate(CandidateList)

  !========================================================================
  ! Post-processing: benchmark-specific modifications
  !========================================================================
#ifdef SED_BENCH
  DResForceX = 0.0
  DResForceY = 0.0
  DResForceZ = 4.0 * DResForceZ
  theParticles(ip)%force(:) = (/DResForceX, DResForceY, DResForceZ/)
  theParticles(ip)%torque(:) = (/0.0, 0.0, 0.0/)
#else
#ifdef ENABLE_LUBRICATION
  if (resi > 0) then
    DResForceX = DResForceX + AlphaRelax * sliX
  end if
#endif
  theParticles(ip)%force(:) = (/DResForceX, DResForceY, DResForceZ/)
  theParticles(ip)%torque(:) = (/DTrqForceX, DTrqForceY, DTrqForceZ/)
#endif
  
  theNorm = calculate_l2_norm(DResForceX, DResForceY, DResForceZ)
  if( theNorm >= maxLocal ) then
    maxLocal = theNorm
  end if
  
#ifdef OUTPUT_LEVEL2
  write(*,'(A,I5,A,3D12.4,A,3D12.4,I3)')'pidx=', theParticles(ip)%bytes(1) + 1, ' theForce    : ', (/DResForceX, DResForceY, DResForceZ/),&
  ' tau: ', (/DTrqForceX, DTrqForceY, DTrqForceZ/), myid
#endif
  
  if (resi > 0) then
    localHydro = localHydro + DResForceX
  end if

END DO ! nParticles

! Output KVEL acceleration statistics
if (bUseKVEL_Accel .and. myid == 1) then
  WRITE(*,'(A,I0,A,I0,A,F8.1,A)') &
    'KVEL: ', myKVEL_Stats%nCandidateElements, ' candidates vs ', &
    NEL*numParticles, ' brute-force (', &
    real(NEL*numParticles)/max(real(myKVEL_Stats%nCandidateElements),1.0), 'x speedup)'
end if

END IF  ! myid /= 0 (end of force calculation section)

!========================================================================
! MPI Force Summation (CFD layer) - ALL RANKS PARTICIPATE
!========================================================================
! In serial PE mode, each CFD domain computes partial hydrodynamic forces
! from its local mesh elements. These contributions must be summed across
! all domains to get the total force on each particle.
! We use the CFD's MPI communicator (not PE's MPI) for this reduction.
! IMPORTANT: COMM_SUMMN is a collective operation requiring ALL ranks!

! Allocate temporary array for MPI reduction (ALL ranks)
allocate(forceArray(6*numParticles))

! Pack forces from theParticles into flat array (workers only)
IF (myid /= 0) THEN
  DO IP = 1, numParticles
    iPointer = 6*(IP-1)
    forceArray(iPointer+1) = theParticles(IP)%force(1)
    forceArray(iPointer+2) = theParticles(IP)%force(2)
    forceArray(iPointer+3) = theParticles(IP)%force(3)
    forceArray(iPointer+4) = theParticles(IP)%torque(1)
    forceArray(iPointer+5) = theParticles(IP)%torque(2)
    forceArray(iPointer+6) = theParticles(IP)%torque(3)
  END DO
ELSE
  ! Rank 0 contributes zero (no force calculation)
  forceArray = 0.0d0
END IF

! MPI sum across all CFD domains (ALL RANKS)
CALL COMM_SUMMN(forceArray, 6*numParticles)

if (myid /= 0)then
  ! Unpack summed forces back into theParticles with scaling factors (ALL ranks)
  DO IP = 1, numParticles
    iPointer = 6*(IP-1)
    theParticles(IP)%force(1) =  forceArray(iPointer+1)
    theParticles(IP)%force(2) =  forceArray(iPointer+2)
    theParticles(IP)%force(3) =  forceArray(iPointer+3)
    theParticles(IP)%torque(1) = forceArray(iPointer+4)
    theParticles(IP)%torque(2) = forceArray(iPointer+5)
    theParticles(IP)%torque(3) = forceArray(iPointer+6)
  
    ! Write the summed+scaled forces to PE bodies (ALL ranks)
    call setForcesMapped(theParticles(ip))
  END DO
  
#ifdef SED_BENCH
  if (myid == 1) then
    time_out = dble(itns - 1) * tstep
    DO IP = 1, numParticles
      iPointer = 6*(IP-1)
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

deallocate(forceArray)

! Deallocate vertex cache (single-timestep lifetime)
if (allocated(ParticleVertexCache)) then
  DO IP = 1, size(ParticleVertexCache)
    if (allocated(ParticleVertexCache(IP)%dofIndices)) then
      deallocate(ParticleVertexCache(IP)%dofIndices)
    end if
  END DO
  deallocate(ParticleVertexCache)
end if

total_lubrication = localSliding

END SUBROUTINE ForcesLocalParticlesSerial
