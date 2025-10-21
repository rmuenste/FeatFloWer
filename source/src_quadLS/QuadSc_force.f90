!************************************************************************
      SUBROUTINE GetForces(factors,U1,U2,U3,P,ALPHA,DVISC,KVERT,KAREA,KEDGE,&
                           DCORVG,ELE)
!************************************************************************
!*     Discrete convection operator: Q1~ elements (nonparametric)
!*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid,showID,COMM_SUMMN
      USE var_QuadScalar, ONLY : myFBM,myExport,Properties
      use cinterface
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
!
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,&
                 NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)

      ! An array of scaling factors
      Real*8 :: factors(*)
!
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
!
      REAL*8 :: DResForceX,DResForceY,DResForceZ
      REAL*8 :: DTrqForceX,DTrqForceY,DTrqForceZ
      REAL*8 :: Center(3),dForce(6)
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

! Main computation (skipped for master process)
      IF (myid /= 0) THEN

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
        call ExitError('Error in GetForces',605)
      end if

      if(myid == 1) write(*,*)'> FBM Force Calculation'
      DO IP = 1,myFBM%nParticles
 
      Center = myFBM%particleNew(IP)%Position
 
      DResForceX = 0D0
      DResForceY = 0D0
      DResForceZ = 0D0
      DTrqForceX = 0d0
      DTrqForceY = 0d0
      DTrqForceZ = 0d0
 
! *** Loop over all elements
      nnel = 0
      DO IEL=1,NEL
!
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER < 0) RETURN
!
!-----------------------------------------------------------------
       NJALFA=0
       NIALFA=0
       DO I=1,IDFL
         IG=KDFG(I)
         IF((ALPHA(IG) == 0).or.(ALPHA(IG) /= IP))THEN
          NJALFA=NJALFA+1
         ENDIF
         IF (ALPHA(IG) == IP) THEN
          NIALFA=NIALFA+1
         ENDIF
       ENDDO

      ! Skip elements where all dofs are inside or
      ! all dofs are outside
      IF(NJALFA == 27.OR.NIALFA == 27) cycle

!      myExport%p_DataScalarCell(1)%pData(iel)=2
      !write(*,*)'adding IELEval: ',IEL
      nnel = nnel + 1
      DNY = DVISC(IEL)

!-----------------------------------------------------------------------
! Setup for trilinear geometry mapping from reference element [-1,1]^3
! to a physical hexahedral element with 8 corner vertices.
!
! Each coordinate (x,y,z) is interpolated as:
!   x(ξ) = ∑_{k=1}^8 N_k(ξ1,ξ2,ξ3) * x_k
!   y(ξ) = ∑_{k=1}^8 N_k(ξ1,ξ2,ξ3) * y_k
!   z(ξ) = ∑_{k=1}^8 N_k(ξ1,ξ2,ξ3) * z_k
!
! The Q1 shape functions N_k(ξ1, ξ2, ξ3) for k = 1..8 are:
!   N_k(ξ) = (1 + ξ1*ξ1k)(1 + ξ2*ξ2k)(1 + ξ3*ξ3k) / 8
! where (ξ1k, ξ2k, ξ3k) ∈ {±1} are the reference coordinates of vertex k.
!
! We precompute the coefficients of the full trilinear polynomial:
!   f(ξ1,ξ2,ξ3) = a₀ + a₁ ξ1 + a₂ ξ2 + a₃ ξ3
!                + a₄ ξ1 ξ2 + a₅ ξ1 ξ3 + a₆ ξ2 ξ3 + a₇ ξ1 ξ2 ξ3
!
! These 8 coefficients are stored for each coordinate:
!   - x: DJm1
!   - y: DJm2
!   - z: DJm3
!
!-----------------------------------------------------------------------      

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

!-----------------------------------------------------------------------
! Precompute monomial coefficients for mapping:
! These encode:
!   a₀      → DJ11/12/13  (constant term)
!   a₁      → DJ21/22/23  (∝ ξ1)
!   a₂      → DJ31/32/33  (∝ ξ2)
!   a₃      → DJ41/42/43  (∝ ξ3)
!   a₄      → DJ51/52/53  (∝ ξ1 ξ2)
!   a₅      → DJ61/62/63  (∝ ξ1 ξ3)
!   a₆      → DJ71/72/73  (∝ ξ2 ξ3)
!   a₇      → DJ81/82/83  (∝ ξ1 ξ2 ξ3)
!
! Used for fast evaluation of:
!   x(ξ1,ξ2,ξ3) = DJ11 + DJ21*ξ1 + DJ31*ξ2 + ... + DJ81*ξ1*ξ2*ξ3
!   (and same for y, z)
!-----------------------------------------------------------------------      
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
      IF (IER < 0) RETURN
!
! *** Loop over all cubature points
      DO ICUBP=1,NCUBP
!
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
!
!-----------------------------------------------------------------------
! Evaluate the Jacobian matrix J = ∂(x,y,z) / ∂(ξ1, ξ2, ξ3)
! at the current cubature (Gauss) point (XI1, XI2, XI3)
!
! This is needed to:
!   - transform basis gradients: ∇_x N = J^{-T} ∇_ξ N
!   - compute integration weight |J| for volume element
!
! The entries are obtained by differentiating the trilinear mapping:
!   x(ξ) = a₀ + a₁ ξ1 + a₂ ξ2 + a₃ ξ3 + ... + a₇ ξ1 ξ2 ξ3
! ⇒ ∂x/∂ξ1 = a₁ + a₄ ξ2 + a₅ ξ3 + a₇ ξ2 ξ3
!    (and similarly for y, z, and other ξ-variables)
!
! Each DJAC(i,j) = ∂(coord i) / ∂(ξ_j)
!-----------------------------------------------------------------------

! Derivatives w.r.t. ξ1 → 1st column of Jacobian      
      DJAC(1,1)=DJ21+DJ51*XI2+DJ61*XI3+DJ81*XI2*XI3
      DJAC(1,2)=DJ31+DJ51*XI1+DJ71*XI3+DJ81*XI1*XI3
      DJAC(1,3)=DJ41+DJ61*XI1+DJ71*XI2+DJ81*XI1*XI2
      DJAC(2,1)=DJ22+DJ52*XI2+DJ62*XI3+DJ82*XI2*XI3
      DJAC(2,2)=DJ32+DJ52*XI1+DJ72*XI3+DJ82*XI1*XI3
      DJAC(2,3)=DJ42+DJ62*XI1+DJ72*XI2+DJ82*XI1*XI2
      DJAC(3,1)=DJ23+DJ53*XI2+DJ63*XI3+DJ83*XI2*XI3
      DJAC(3,2)=DJ33+DJ53*XI1+DJ73*XI3+DJ83*XI1*XI3
      DJAC(3,3)=DJ43+DJ63*XI1+DJ73*XI2+DJ83*XI1*XI2

!-----------------------------------------------------------------------
! Compute determinant of Jacobian |J| at Gauss point
! (used to convert volume integrals from reference to physical space)
! Also needed for transforming gradients: ∇x = J^{-T} ∇ξ
!-----------------------------------------------------------------------      
      DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))&
           -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))&
           +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))

! Multiply with Gauss weight to get physical-space integration weight           
      OM=DOMEGA(ICUBP)*ABS(DETJ)
!
!-----------------------------------------------------------------------
! Evaluate the physical coordinates (x, y, z) of the current Gauss point
! by partially reconstructing the trilinear mapping from reference point
! (XI1, XI2, XI3) ∈ [-1,1]³ to physical coordinates (XX, YY, ZZ).
!
! The full trilinear mapping has 8 monomial terms:
!   f(ξ) = a₀ + a₁ ξ1 + a₂ ξ2 + a₃ ξ3 + a₄ ξ1ξ2 + a₅ ξ1ξ3 + a₆ ξ2ξ3 + a₇ ξ1ξ2ξ3
!
! Here, the code **reuses DJAC(i,i)*ξᵢ** to implicitly insert:
!   a₁ ξ1 + a₄ ξ1ξ2 + a₅ ξ1ξ3 + a₇ ξ1ξ2ξ3  for x (via DJAC(1,1)*XI1)
!   a₂ ξ2 + a₄ ξ1ξ2 + a₆ ξ2ξ3 + a₇ ξ1ξ2ξ3  for y (via DJAC(2,2)*XI2)
!   a₃ ξ3 + a₅ ξ1ξ3 + a₆ ξ2ξ3 + a₇ ξ1ξ2ξ3  for z (via DJAC(3,3)*XI3)
!
! It then adds the remaining terms explicitly:
!   a₂ ξ2, a₃ ξ3, a₆ ξ2ξ3 for x
!   a₁ ξ1, a₃ ξ3, a₅ ξ1ξ3 for y
!   a₁ ξ1, a₂ ξ2, a₄ ξ1ξ2 for z
!
!-----------------------------------------------------------------------      
      XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
      YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
      ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2

!-----------------------------------------------------------------------
! Evaluate Q2 (27-node) basis functions and their derivatives at the
! current Gauss point (XI1, XI2, XI3) ∈ [-1,1]^3 in the reference element.
!
! What this call does:
!   - Shape functions N_k(ξ) are defined for k = 1..27 (Q2 element)
!   - Computes:
!       * Value         :     N_k(ξ)
!       * Gradient in x : ∂N_k/∂x
!       * Gradient in y : ∂N_k/∂y
!       * Gradient in z : ∂N_k/∂z
!
!   - Applies chain rule using the inverse Jacobian:
!       ∇_x N_k = J^{-T} · ∇_ξ N_k
!
! Output is stored in:
!   DBAS(dim,k,deriv) where:
!     dim   = 1 (scalar field)
!     k     = 1..27 (local basis function index)
!     deriv = 1: value         → N_k(ξ)
!             2: ∂N_k/∂x
!             3: ∂N_k/∂y
!             4: ∂N_k/∂z
!
! Usage example:
!   DO I = 1, IDFL   (IDFL = 27 for Q2)
!     k  = KDFL(I)        ! local shape function index
!     DBI1 = DBAS(1,k,1)  ! N_k(ξ)
!     DBI2 = DBAS(1,k,2)  ! ∂N_k/∂x
!     DBI3 = DBAS(1,k,3)  ! ∂N_k/∂y
!     DBI4 = DBAS(1,k,4)  ! ∂N_k/∂z
!   END DO
!
! Note:
!   - These are the scalar-valued Q2 basis functions (e.g. for α, U1, U2, U3)
!   - The geometry mapping (earlier) is Q1 with only 8 nodes → don’t confuse!
!-----------------------------------------------------------------------
!
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER < 0) RETURN
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

!-----------------------------------------------------------------------
! Interpolate field values and spatial derivatives at the current
! Gauss (cubature) point using Q2 shape functions.
!
! This loop performs the **standard finite element interpolation** of
! multiple scalar fields (U1, U2, U3, ALPHA) and their gradients.
!
! Each field is represented in the element as:
!     f_h(x) = ∑_{k=1}^{27} f_k · N_k(x)
! and its gradient as:
!     ∇f_h(x) = ∑_{k=1}^{27} f_k · ∇N_k(x)
!
! What the loop does:
!   - IDFL = 27 is the number of Q2 basis functions per element
!   - KDFG(I) = global DOF index for I-th local basis function
!   - KDFL(I) = local index k ∈ 1..27 used to index DBAS
!
! For each local basis function N_k:
!   - DBI1 = N_k(ξ)       (basis function value at Gauss point)
!   - DBI2 = ∂N_k/∂x      (basis function x-derivative)
!   - DBI3 = ∂N_k/∂y      (basis function y-derivative)
!   - DBI4 = ∂N_k/∂z      (basis function z-derivative)
!
! These are used to interpolate:
!   - U1, U2, U3: components of velocity vector field
!   - ALPHA: indicator function αᵢ(x) for particle i
!
! Result:
!   After the loop, you obtain:
!     - DU1V = u₁(x_q),   DU1X = ∂u₁/∂x, etc.
!     - DALV = αᵢ(x_q),   DALX = ∂αᵢ/∂x, etc.
!   at the current Gauss point
!
! Notes:
!   - The interpolated gradients are used to compute ∇u and ∇α
!   - ALPHA is stored as an integer indicator (per node):
!       ALPHA(IG) == IP → node belongs to current particle
!   - DALPHA is set to 1 only if current particle matches node marking
!     ⇒ gives piecewise-constant αᵢ, smoothed to Gauss point by shape fcts
!-----------------------------------------------------------------------       
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
         IF (ALPHA(IG) == IP) THEN
          DALPHA = 1d0
         ELSE
          DALPHA = 0d0
         END IF
         DALV=DALV+DALPHA*DBI1
         DALX=DALX+DALPHA*DBI2
         DALY=DALY+DALPHA*DBI3
         DALZ=DALZ+DALPHA*DBI4
       ENDDO

!-----------------------------------------------------------------------
! Interpolate the pressure at the current Gauss point using the
! element-wise linear pressure representation (P1-discontinuous basis).
!
! Each pressure field is stored with 4 coefficients per element:
!   - P(JJ)     : constant term p₀
!   - P(JJ+1)   : linear coefficient in x-direction (∝ x - x₀)
!   - P(JJ+2)   : linear coefficient in y-direction (∝ y - y₀)
!   - P(JJ+3)   : linear coefficient in z-direction (∝ z - z₀)
!
! The pressure is evaluated at the Gauss point coordinates (XX, YY, ZZ)
! using the element center (DX0, DY0, DZ0) as the expansion point:
!
!     p(x,y,z) ≈ p₀ + (x - x₀) * px + (y - y₀) * py + (z - z₀) * pz
!
! where (x₀, y₀, z₀) = element centroid
!
! This linear field is element-local and does not enforce continuity
! between elements (which is appropriate for pressure in incompressible
! velocity-pressure formulations using e.g. Q2/P1_disc pairing).
!-----------------------------------------------------------------------
!C---------------------------------------------------------
       JJ = 4*(IEL-1) + 1
       Press =          P(JJ  ) + (XX-DX0)*P(JJ+1) + &
               (YY-DY0)*P(JJ+2) + (ZZ-DZ0)*P(JJ+3)
!--------------------------------------------------------
!-----------------------------------------------------------------------
! Evaluate the integrand for the hydrodynamic force and torque at the
! current Gauss point (XX, YY, ZZ).
!
! Based on the volume-integrated expression from the paper:
!     Fᵢ = - ∫ σ_h · ∇αᵢ dV
! where
!     σ_h = -p I + μ (∇u + ∇uᵀ)     -- full Newtonian Cauchy stress
!
! Quantities needed at the quadrature point:
!   - Velocity u = (u₁, u₂, u₃) and gradients ∇u
!   - Pressure p
!   - Indicator function αᵢ and gradient ∇αᵢ
!   - Viscosity μ = DNY or DMU0 (interpolated or element-based)
!
! At this point in the code:
!   - DU1X, DU1Y, DU1Z = ∂u₁/∂x, ∂u₁/∂y, ∂u₁/∂z
!   - DU2X, DU2Y, DU2Z = ∂u₂/∂x, ∂u₂/∂y, ∂u₂/∂z
!   - DU3X, DU3Y, DU3Z = ∂u₃/∂x, ∂u₃/∂y, ∂u₃/∂z
!   - DALX, DALY, DALZ = ∂αᵢ/∂x, ∂αᵢ/∂y, ∂αᵢ/∂z
!
! The integrand is:
!     f = σ_h · (-∇αᵢ) = -p (-∇αᵢ) + μ (∇u + ∇uᵀ) · (-∇αᵢ)
!
! This code uses the **symmetric** (physically correct) form:
!     σ_h = -p I + μ (∇u + ∇uᵀ)
!
! NOTE:
!     The non-symmetric approximation σ_h ≈ -p I + μ ∇u
!     is available in the commented-out block below, but is
!     less accurate for computing hydrodynamic forces.
!-----------------------------------------------------------------------

!-----------Form the integrand------------------
       DN1=-DALX
       DN2=-DALY
       DN3=-DALZ
!-------------------Acting force-------------------------
!--------------Symmetric stress tensor calculation-------------
!
! Non-symmetric form (less accurate, provided for reference):
!       AH1=-Press*DN1 + DNY*(DU1X*DN1 + DU1Y*DN2 + DU1Z*DN3)
!       AH2=-Press*DN2 + DNY*(DU2X*DN1 + DU2Y*DN2 + DU2Z*DN3)
!       AH3=-Press*DN3 + DNY*(DU3X*DN1 + DU3Y*DN2 + DU3Z*DN3)
!
! Symmetric form (active, physically correct):
        AH1=-Press*DN1+DNY*((DU1X+DU1X)*DN1+(DU1Y+DU2X)*DN2 + & ! full3D
           (DU1Z+DU3X)*DN3)
        AH2=-Press*DN2+DNY*((DU2X+DU1Y)*DN1+(DU2Y+DU2Y)*DN2 + & ! full3D
           (DU2Z+DU3Y)*DN3)
        AH3=-Press*DN3+DNY*((DU3X+DU1Z)*DN1+(DU3Y+DU2Z)*DN2 + & ! full3D
           (DU3Z+DU3Z)*DN3)

!--------------------------------------------------------------
       DResForceX = DResForceX + AH1*OM
       DResForceY = DResForceY + AH2*OM
       DResForceZ = DResForceZ + AH3*OM
!-----------------------------------------------------------------------
! Evaluate the torque contribution at this Gauss point via:
!
!     τ_i = ∫ (x - Xᵢ) × (σ_h · ∇αᵢ) dV
!
! At this point:
!   - AH = σ_h · (−∇αᵢ) has just been computed
!   - (XX, YY, ZZ) = physical coordinates of Gauss point
!   - Center(:) = center of mass of particle i
!
! The torque integrand is the cross product:
!     τ = r × f = (x - Xᵢ) × AH
!
! This gives:
!     τ₁ = (y - Yᵢ)*AH₃ - (z - Zᵢ)*AH₂
!     τ₂ = (z - Zᵢ)*AH₁ - (x - Xᵢ)*AH₃
!     τ₃ = (x - Xᵢ)*AH₂ - (y - Yᵢ)*AH₁
!
! The result is then accumulated using the quadrature weight OM.
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
!
      end do ! end loop elements
!
      iPointer = 6*(IP-1)

      myFBM%Force(iPointer+1) = DResForceX
      myFBM%Force(iPointer+2) = DResForceY
      myFBM%Force(iPointer+3) = DResForceZ

      myFBM%Force(iPointer+4) = DTrqForceX
      myFBM%Force(iPointer+5) = DTrqForceY
      myFBM%Force(iPointer+6) = DTrqForceZ

      END DO ! nParticles

      END IF ! myid /= 0

999   CALL COMM_SUMMN(myFBM%Force,6*myFBM%nParticles)

      DO IP = 1,myFBM%nParticles
       iPointer = 6*(IP-1)+1

       ! Gather translational force
       myFBM%ParticleNew(IP)%ResistanceForce(1)= &
       factors(1) * &
       myFBM%Force(iPointer)

       myFBM%ParticleNew(IP)%ResistanceForce(2) = &
       factors(2) * &
       myFBM%Force(iPointer+1)

       myFBM%ParticleNew(IP)%ResistanceForce(3) = &
       factors(3) * &
       myFBM%Force(iPointer+2)

       ! Gather rotational force
       myFBM%ParticleNew(IP)%TorqueForce(1) = &
       factors(4) * myFBM%Force(iPointer+3)

       myFBM%ParticleNew(IP)%TorqueForce(2) = &
       factors(5) * myFBM%Force(iPointer+4)

       myFBM%ParticleNew(IP)%TorqueForce(3) = &
       factors(6) * myFBM%Force(iPointer+5)
      END DO

      END
!
!-----------------------------------------------------------------------
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
!
!
      SUBROUTINE ExtractBoundaryNormals_Q2(NU,NV,NW,DD,KVERT,&
                 KAREA,KEDGE,DCORVG,ELE)
      USE PP3D_MPI, ONLY:myid
      USE var_QuadScalar, ONLY : myBoundary
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
!
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,&
                 NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
!

      REAL*8 NU(*),NV(*),NW(*)
      REAL*8 DD(*), DCORVG(NNDIM,*)
      INTEGER  KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      REAL*8 DMyOmgP(NNCUBP),DMyCubP(NNCUBP,NNAE,NNDIM)
      REAL*8    TSTEP,dN(3)
      REAL*8    DHELP(NNBAS,4,NNCUBP),DPP(NNDIM),dNorm(3)
      TYPE tF
       REAL*8   DHELP(NNBAS,4,NNCUBP)
      END TYPE
      TYPE(tF) F(6)
      REAL*8    dFluidNormal,PointC(3),PointA(3),dLine(3),dVELO(3)
      INTEGER NeighA(4,6)
      DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

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
      ICUB = 3 
      CALL SetUpMyCub(DMyOmgP,DMyCubP,NCUBP,ICUB)
!
      DO IAT=1,NNAE
      DO ICUBP=1,NCUBP
       XI1=DMyCubP(ICUBP,IAT,1)
       XI2=DMyCubP(ICUBP,IAT,2)
       XI3=DMyCubP(ICUBP,IAT,3)
       CALL E013A(XI1,XI2,XI3,DHELP,ICUBP)
      END DO
      F(IAT)%DHELP = DHELP
      END DO
!
!***** ******************************************************************
! *** Calculation of the matrix - storage technique 7 or 8
!***** ******************************************************************
!
! *** Loop over all elements
      DO 100 IEL=1,NEL
!      
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
!
      IF (IER.LT.0) GOTO 99999
!
!       IF (myBoundary%iPhase(NVT+NET+NAT+IEL).ne.0) THEN
       dFluidNormal = +1d0
!       ELSE
!        dFluidNormal = 0d0
!       END IF
! C     
      PointC(:) = DCORVG(:,NVT+NET+NAT+IEL)
!      
      DO 150 IAT=1,6
!
      ivt1 = kvert(NeighA(1,IAT),IEL)
      ivt2 = kvert(NeighA(2,IAT),IEL)
      ivt3 = kvert(NeighA(3,IAT),IEL)
      ivt4 = kvert(NeighA(4,IAT),IEL)
      
      IF (myBoundary%bSlip(ivt1).and.&
          myBoundary%bSlip(ivt2).and.&
          myBoundary%bSlip(ivt3).and.&
          myBoundary%bSlip(ivt4)) THEN
!      
      PointA(:) = DCORVG(:,NVT+NET+KAREA(IAT,IEL))
      dLine(:)  = PointC(:) - PointA(:)
!      
      DHELP = F(IAT)%DHELP
!      
      DO 200 ICUBP=1,NCUBP
!
      XI1=DMyCubP(ICUBP,IAT,1)
      XI2=DMyCubP(ICUBP,IAT,2)
      XI3=DMyCubP(ICUBP,IAT,3)
      OM=DMyOmgP (ICUBP)
!
      DJAC=0d0
      DO JDOFE=1,IDFL
       JDFL=KDFL(JDOFE)
       JDFG=KDFG(JDOFE)
       DPP(:) = DCORVG(:,JDFG)
       DJAC(1,1)= DJAC(1,1) +  DPP(1)*DHELP(JDFL,2,ICUBP)
       DJAC(2,1)= DJAC(2,1) +  DPP(2)*DHELP(JDFL,2,ICUBP)
       DJAC(3,1)= DJAC(3,1) +  DPP(3)*DHELP(JDFL,2,ICUBP)
       DJAC(1,2)= DJAC(1,2) +  DPP(1)*DHELP(JDFL,3,ICUBP)
       DJAC(2,2)= DJAC(2,2) +  DPP(2)*DHELP(JDFL,3,ICUBP)
       DJAC(3,2)= DJAC(3,2) +  DPP(3)*DHELP(JDFL,3,ICUBP)
       DJAC(1,3)= DJAC(1,3) +  DPP(1)*DHELP(JDFL,4,ICUBP)
       DJAC(2,3)= DJAC(2,3) +  DPP(2)*DHELP(JDFL,4,ICUBP)
       DJAC(3,3)= DJAC(3,3) +  DPP(3)*DHELP(JDFL,4,ICUBP)
      END DO
!
      DETJ = DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))&
            -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))&
            +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
!
      IF (IAT.eq.2.or.iat.eq.4) CALL SurfDet(DJAC,2,dA,dN)
      IF (IAT.eq.1.or.iat.eq.6) CALL SurfDet(DJAC,3,dA,dN)
      IF (IAT.eq.3.or.iat.eq.5) CALL SurfDet(DJAC,1,dA,dN)
!
      CALL ELE(XI1,XI2,XI3,0)
      IF (IER.LT.0) GOTO 99999
!
      dNorm = dFluidNormal*dN
      dDotProd = dNorm(1)*dLine(1)+dNorm(2)*dLine(2)+dNorm(3)*dLine(3)
      IF (dDotProd.lt.0d0) dNorm = -dNorm

      DO I=1,IDFL
        IG = KDFG(I)
        IL = KDFL(I)
        HBAS = DBAS(1,IL,1)
        NU(IG) = NU(IG) - dA*OM*dNorm(1)*HBAS
        NV(IG) = NV(IG) - dA*OM*dNorm(2)*HBAS
        NW(IG) = NW(IG) - dA*OM*dNorm(3)*HBAS
        DD(IG) = DD(IG) + dA*OM*HBAS
      ENDDO

200   CONTINUE
      END IF
150   CONTINUE
100   CONTINUE
99999 END
!
!
!
      SUBROUTINE ExtractBoundaryNormals_Q1(NU,NV,NW,DD,KVERT,&
                 KAREA,KEDGE,DCORVG,ELE)
      USE PP3D_MPI, ONLY:myid
      USE var_QuadScalar, ONLY : myBoundary
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
!
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,&
                 NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
!

      REAL*8 NU(*),NV(*),NW(*)
      REAL*8 DD(*), DCORVG(NNDIM,*)
      INTEGER  KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      REAL*8 DMyOmgP(NNCUBP),DMyCubP(NNCUBP,NNAE,NNDIM)
      REAL*8    TSTEP,dN(3)
      REAL*8    DHELP(8,4,NNCUBP),DPP(NNDIM),dNorm(3)
      TYPE tF
       REAL*8   DHELP(8,4,NNCUBP)
      END TYPE
      TYPE(tF) F(6)
      REAL*8    dFluidNormal,PointC(3),PointA(3),dLine(3),dVELO(3)
      INTEGER NeighA(4,6)
      DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

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
      ICUB = 3 
      CALL SetUpMyCub(DMyOmgP,DMyCubP,NCUBP,ICUB)
!
      DO IAT=1,NNAE
      DO ICUBP=1,NCUBP
       XI1=DMyCubP(ICUBP,IAT,1)
       XI2=DMyCubP(ICUBP,IAT,2)
       XI3=DMyCubP(ICUBP,IAT,3)
       CALL E011A(XI1,XI2,XI3,DHELP,ICUBP)
      END DO
      F(IAT)%DHELP = DHELP
      END DO
!
!***** ******************************************************************
! *** Calculation of the matrix - storage technique 7 or 8
!***** ******************************************************************
!
! *** Loop over all elements
      DO 100 IEL=1,NEL
!      
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
!
      IF (IER.LT.0) GOTO 99999
!
!       IF (myBoundary%iPhase(NVT+NET+NAT+IEL).ne.0) THEN
       dFluidNormal = +1d0
!       ELSE
!        dFluidNormal = 0d0
!       END IF
! C     
      PointC(:) = DCORVG(:,NVT+NET+NAT+IEL)
!      
      DO 150 IAT=1,6
!
      ivt1 = kvert(NeighA(1,IAT),IEL)
      ivt2 = kvert(NeighA(2,IAT),IEL)
      ivt3 = kvert(NeighA(3,IAT),IEL)
      ivt4 = kvert(NeighA(4,IAT),IEL)
      
      IF (myBoundary%bSlip(ivt1).and.&
          myBoundary%bSlip(ivt2).and.&
          myBoundary%bSlip(ivt3).and.&
          myBoundary%bSlip(ivt4)) THEN
!      
      PointA(:) = DCORVG(:,NVT+NET+KAREA(IAT,IEL))
      dLine(:)  = PointC(:) - PointA(:)
!      
      DHELP = F(IAT)%DHELP
!      
      DO 200 ICUBP=1,NCUBP
!
      XI1=DMyCubP(ICUBP,IAT,1)
      XI2=DMyCubP(ICUBP,IAT,2)
      XI3=DMyCubP(ICUBP,IAT,3)
      OM=DMyOmgP (ICUBP)
!
      DJAC=0d0
      DO JDOFE=1,8
       JDFL=KDFL(JDOFE)
       JDFG=KDFG(JDOFE)
       DPP(:) = DCORVG(:,JDFG)
       DJAC(1,1)= DJAC(1,1) +  DPP(1)*DHELP(JDFL,2,ICUBP)
       DJAC(2,1)= DJAC(2,1) +  DPP(2)*DHELP(JDFL,2,ICUBP)
       DJAC(3,1)= DJAC(3,1) +  DPP(3)*DHELP(JDFL,2,ICUBP)
       DJAC(1,2)= DJAC(1,2) +  DPP(1)*DHELP(JDFL,3,ICUBP)
       DJAC(2,2)= DJAC(2,2) +  DPP(2)*DHELP(JDFL,3,ICUBP)
       DJAC(3,2)= DJAC(3,2) +  DPP(3)*DHELP(JDFL,3,ICUBP)
       DJAC(1,3)= DJAC(1,3) +  DPP(1)*DHELP(JDFL,4,ICUBP)
       DJAC(2,3)= DJAC(2,3) +  DPP(2)*DHELP(JDFL,4,ICUBP)
       DJAC(3,3)= DJAC(3,3) +  DPP(3)*DHELP(JDFL,4,ICUBP)
      END DO
!
      DETJ = DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))&
            -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))&
            +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
!
      IF (IAT.eq.2.or.iat.eq.4) CALL SurfDet(DJAC,2,dA,dN)
      IF (IAT.eq.1.or.iat.eq.6) CALL SurfDet(DJAC,3,dA,dN)
      IF (IAT.eq.3.or.iat.eq.5) CALL SurfDet(DJAC,1,dA,dN)
!
      CALL ELE(XI1,XI2,XI3,0)
      IF (IER.LT.0) GOTO 99999
!
      dNorm = dFluidNormal*dN
      dDotProd = dNorm(1)*dLine(1)+dNorm(2)*dLine(2)+dNorm(3)*dLine(3)
      IF (dDotProd.lt.0d0) dNorm = -dNorm

      DO I=1,IDFL
        IG = KDFG(I)
        IL = KDFL(I)
        HBAS = DBAS(1,IL,1)
        NU(IG) = NU(IG) - dA*OM*dNorm(1)*HBAS
        NV(IG) = NV(IG) - dA*OM*dNorm(2)*HBAS
        NW(IG) = NW(IG) - dA*OM*dNorm(3)*HBAS
        DD(IG) = DD(IG) + dA*OM*HBAS
      ENDDO

200   CONTINUE
      END IF
150   CONTINUE
100   CONTINUE
99999 END



      
