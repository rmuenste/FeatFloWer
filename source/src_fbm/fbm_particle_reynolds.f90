MODULE fbm_particle_reynolds

  USE PP3D_MPI, ONLY: myid, COMM_SUMMN
  USE var_QuadScalar, ONLY: myFBM, mg_mesh

  IMPLICIT NONE

CONTAINS

#ifdef HAVE_PE
  SUBROUTINE fbm_compute_particle_reynolds(U1, U2, U3, DVISC, RHOFLUID, mfile, ELE)
    USE fbmaux, ONLY: fbmaux_PointInHex
    use dem_query, only: getAllParticles, tParticleData
    INTEGER, PARAMETER :: NNBAS = 27, NNDER = 10, NNVE = 8, NNDIM = 3

    REAL*8, DIMENSION(:), INTENT(IN) :: U1, U2, U3
    REAL*8, DIMENSION(:), INTENT(IN) :: DVISC
    REAL*8, INTENT(IN) :: RHOFLUID
    INTEGER, INTENT(IN) :: mfile
    EXTERNAL :: ELE  ! Element evaluation function (e.g., E013)

    INTEGER :: ip, ig, il, IELTYP, ilev, nel_local
    REAL*8 :: xi1, xi2, xi3
    REAL*8 :: vel_sample(3), slip(3), speed, diameter, mu_loc
    REAL*8 :: xverts(8), yverts(8), zverts(8)
    REAL*8 :: xmin, xmax, ymin, ymax, zmin, zmax, eps_box
    LOGICAL :: found

    TYPE(tParticleData), DIMENSION(:), ALLOCATABLE :: theParticles
    REAL*8, ALLOCATABLE :: re_local(:), re_weight(:)

    INTEGER KDFG(NNBAS), KDFL(NNBAS)
    REAL*8 DX(NNVE), DY(NNVE), DZ(NNVE), DJAC(3, 3), DETJ
    REAL*8 DBAS(NNDIM, NNBAS, NNDER)
    LOGICAL BDER(NNDER)
    INTEGER KVE(NNVE), NDIM, IEL
    INTEGER :: ive  ! Loop variable for element vertices

    COMMON/ELEM/DX, DY, DZ, DJAC, DETJ, DBAS, BDER, KVE, IEL, NDIM
    COMMON/COAUX1/KDFG, KDFL, IDFL
    INTEGER :: IDFL

    EXTERNAL :: NDFGL, SETLEV
    INTEGER :: NDFL

    IF (.NOT. ALLOCATED(myFBM%ParticleRe)) RETURN
    IF (myFBM%nParticles .EQ. 0) THEN
      myFBM%ParticleRe = 0D0
      RETURN
    END IF

    ilev = mg_mesh%nlmax
    CALL SETLEV(2)
    nel_local = mg_mesh%level(ilev)%nel
    IF (SIZE(DVISC) .LT. nel_local) THEN
      myFBM%ParticleRe = 0D0
      RETURN
    END IF

    ALLOCATE (theParticles(myFBM%nParticles))
    ALLOCATE (re_local(myFBM%nParticles))
    ALLOCATE (re_weight(myFBM%nParticles))
    re_local = 0D0
    re_weight = 0D0

    DO ive = 1, NNDER
      BDER(ive) = .FALSE.
    END DO
    BDER(1) = .TRUE.

    IELTYP = -1
    CALL ELE(0D0, 0D0, 0D0, IELTYP)
    IDFL = NDFL(IELTYP)

    eps_box = 1D-10

    IF (myid .NE. 0) THEN
      ! Fetch particle data from PE library
      CALL getAllParticles(theParticles)

      DO ip = 1, myFBM%nParticles
        vel_sample = 0D0
        found = .FALSE.
        xi1 = 0D0
        xi2 = 0D0
        xi3 = 0D0

        DO iel = 1, nel_local
          DO ive = 1, NNVE
            ig = mg_mesh%level(ilev)%kvert(ive, iel)
            xverts(ive) = mg_mesh%level(ilev)%dcorvg(1, ig)
            yverts(ive) = mg_mesh%level(ilev)%dcorvg(2, ig)
            zverts(ive) = mg_mesh%level(ilev)%dcorvg(3, ig)
          END DO

          xmin = MINVAL(xverts)
          xmax = MAXVAL(xverts)
          ymin = MINVAL(yverts)
          ymax = MAXVAL(yverts)
          zmin = MINVAL(zverts)
          zmax = MAXVAL(zverts)

          IF (myFBM%ParticleNew(ip)%Position(1) .LT. xmin - eps_box) CYCLE
          IF (myFBM%ParticleNew(ip)%Position(1) .GT. xmax + eps_box) CYCLE
          IF (myFBM%ParticleNew(ip)%Position(2) .LT. ymin - eps_box) CYCLE
          IF (myFBM%ParticleNew(ip)%Position(2) .GT. ymax + eps_box) CYCLE
          IF (myFBM%ParticleNew(ip)%Position(3) .LT. zmin - eps_box) CYCLE
          IF (myFBM%ParticleNew(ip)%Position(3) .GT. zmax + eps_box) CYCLE

          xi1 = 0D0
          xi2 = 0D0
          xi3 = 0D0
          found = fbmaux_PointInHex(myFBM%ParticleNew(ip)%Position(1), &
                                    myFBM%ParticleNew(ip)%Position(2), &
                                    myFBM%ParticleNew(ip)%Position(3), &
                                    xverts, yverts, zverts, xi1, xi2, xi3, iel)
          IF (.NOT. found) CYCLE

          DO ive = 1, NNDER
            BDER(ive) = .FALSE.
          END DO
          BDER(1) = .TRUE.

          CALL NDFGL(iel, 1, IELTYP, mg_mesh%level(ilev)%kvert, &
                     mg_mesh%level(ilev)%kedge, mg_mesh%level(ilev)%karea, &
                     KDFG, KDFL)

          DO ive = 1, NNVE
            ig = mg_mesh%level(ilev)%kvert(ive, iel)
            KVE(ive) = ig
            DX(ive) = xverts(ive)
            DY(ive) = yverts(ive)
            DZ(ive) = zverts(ive)
          END DO

          CALL ELE(xi1, xi2, xi3, 0)

          vel_sample = 0D0
          DO il = 1, IDFL
            ig = KDFG(il)
            ive = KDFL(il)
            vel_sample(1) = vel_sample(1) + U1(ig)*DBAS(1, ive, 1)
            vel_sample(2) = vel_sample(2) + U2(ig)*DBAS(1, ive, 1)
            vel_sample(3) = vel_sample(3) + U3(ig)*DBAS(1, ive, 1)
          END DO

          mu_loc = DVISC(iel)
          slip = vel_sample - myFBM%ParticleNew(ip)%Velocity
          speed = SQRT(slip(1)**2 + slip(2)**2 + slip(3)**2)
          ! Get diameter from particle AABB (largest extent)
          diameter = MAXVAL(theParticles(ip)%aabb)

          IF (mu_loc .GT. 0D0 .AND. diameter .GT. 0D0) THEN
            re_local(ip) = RHOFLUID*speed*diameter/mu_loc
          ELSE
            re_local(ip) = 0D0
          END IF
          re_weight(ip) = 1D0
          IF (ALLOCATED(myFBM%iel_ug)) THEN
            IF (ip .LE. SIZE(myFBM%iel_ug)) myFBM%iel_ug(ip) = iel
          END IF
          EXIT
        END DO
      END DO
    END IF

    CALL COMM_SUMMN(re_local, myFBM%nParticles)
    CALL COMM_SUMMN(re_weight, myFBM%nParticles)

    DO ip = 1, myFBM%nParticles
      IF (re_weight(ip) .GT. 0D0) THEN
        myFBM%ParticleRe(ip) = re_local(ip)/re_weight(ip)
      ELSE
        myFBM%ParticleRe(ip) = 0D0
      END IF
    END DO

    IF (myid .EQ. 1) THEN
      IF (myFBM%nParticles .GT. 0) THEN
        WRITE (mfile, '(A,2E16.8)') 'Particle Reynolds number range: ', &
          MINVAL(myFBM%ParticleRe), MAXVAL(myFBM%ParticleRe)
      END IF
    END IF

    DEALLOCATE (theParticles)
    DEALLOCATE (re_local)
    DEALLOCATE (re_weight)

  END SUBROUTINE fbm_compute_particle_reynolds
#else
  ! Stub implementation when PE library is not available
  SUBROUTINE fbm_compute_particle_reynolds(U1, U2, U3, DVISC, RHOFLUID, mfile, ELE)
    REAL*8, DIMENSION(:), INTENT(IN) :: U1, U2, U3
    REAL*8, DIMENSION(:), INTENT(IN) :: DVISC
    REAL*8, INTENT(IN) :: RHOFLUID
    INTEGER, INTENT(IN) :: mfile
    EXTERNAL :: ELE
    ! No-op: PE library not available, cannot compute particle Reynolds numbers
    IF (ALLOCATED(myFBM%ParticleRe)) myFBM%ParticleRe = 0D0
    RETURN
  END SUBROUTINE fbm_compute_particle_reynolds
#endif

!=========================================================================
! fbm_compute_particle_reynolds_interface - Interface-based Reynolds calculation
!
! Computes particle Reynolds numbers by sampling fluid velocity at interface
! elements (elements that cross the particle boundary) rather than at the
! particle center. This approach is suitable for fully resolved DNS where
! the particle center velocity is constrained to match particle motion.
!
! Method: For each particle, identify all interface elements (elements with
! some DOFs inside and some outside the particle), evaluate velocity at
! element centers, average these velocities, and compute Reynolds number
! from the slip velocity.
!
! Input:
!   U1, U2, U3  - Velocity field components (Q2)
!   ALPHA       - Particle indicator field (0 for fluid, particle ID for solid)
!   DVISC       - Element-wise dynamic viscosity
!   RHOFLUID    - Fluid density
!   mfile       - Output file handle
!
! Output:
!   myFBM%ParticleRe(:) - Reynolds number for each particle
!=========================================================================
#ifdef HAVE_PE
  SUBROUTINE fbm_compute_particle_reynolds_interface(U1, U2, U3, ALPHA, DVISC, RHOFLUID, &
                                                     mfile, ELE)
    !USE cinterface, ONLY: tParticleData
    use dem_query, only: numTotalParticles, getAllParticles, tParticleData 
    USE PP3D_MPI, ONLY: COMM_Maximumn
    INTEGER, PARAMETER :: NNBAS = 27, NNDER = 10, NNVE = 8, NNDIM = 3, NNEE = 12, NNAE = 6

    REAL*8, DIMENSION(:), INTENT(IN) :: U1, U2, U3
    INTEGER, DIMENSION(:), INTENT(IN) :: ALPHA
    REAL*8, DIMENSION(:), INTENT(IN) :: DVISC
    REAL*8, INTENT(IN) :: RHOFLUID
    INTEGER, INTENT(IN) :: mfile
    EXTERNAL :: ELE  ! Element evaluation function (e.g., E013)

    INTEGER :: ip, ig, il, IELTYP, ilev, nel_local, numParticles
    INTEGER :: NJALFA, NIALFA, n_interface, i
    REAL*8 :: vel_sum(3), vel_sample(3), vel_avg(3), slip(3)
    REAL*8 :: speed, diameter, mu_avg, weight_sum
    REAL*8 :: dbuf1(1)

    TYPE(tParticleData), DIMENSION(:), ALLOCATABLE :: theParticles
    REAL*8, ALLOCATABLE :: re_local(:), re_weight(:)

    INTEGER KDFG(NNBAS), KDFL(NNBAS)
    REAL*8 DX(NNVE), DY(NNVE), DZ(NNVE), DJAC(3, 3), DETJ
    REAL*8 DBAS(NNDIM, NNBAS, NNDER)
    LOGICAL BDER(NNDER)
    INTEGER KVE(NNVE), NDIM, IEL
    INTEGER :: ive  ! Loop variable for element vertices

    COMMON/ELEM/DX, DY, DZ, DJAC, DETJ, DBAS, BDER, KVE, IEL, NDIM
    COMMON/COAUX1/KDFG, KDFL, IDFL
    INTEGER :: IDFL

    EXTERNAL :: NDFGL, SETLEV
    INTEGER :: NDFL

    ! Get total particle count (ALL ranks need this for MPI)
    ! Rank 0 does not run PE and would see 0. Use collective max so all ranks agree.
    numParticles = numTotalParticles()
    dbuf1(1) = DBLE(numParticles)
    CALL COMM_Maximumn(dbuf1, 1)
    numParticles = INT(dbuf1(1))

    ! Early return if no particles
    IF (numParticles .EQ. 0) RETURN

    ! Allocate or reallocate ParticleRe array if needed
    IF (.NOT. ALLOCATED(myFBM%ParticleRe)) THEN
      ALLOCATE (myFBM%ParticleRe(numParticles))
      myFBM%ParticleRe = 0D0
    ELSE IF (SIZE(myFBM%ParticleRe) .NE. numParticles) THEN
      DEALLOCATE (myFBM%ParticleRe)
      ALLOCATE (myFBM%ParticleRe(numParticles))
      myFBM%ParticleRe = 0D0
    END IF

    ! Allocate particle data and result arrays (ALL ranks)
    ALLOCATE (theParticles(numParticles))
    ALLOCATE (re_local(numParticles))
    ALLOCATE (re_weight(numParticles))
    re_local = 0D0
    re_weight = 0D0

    ilev = mg_mesh%nlmax
    CALL SETLEV(2)
    nel_local = mg_mesh%level(ilev)%nel
    IF (SIZE(DVISC) .LT. nel_local) THEN
      myFBM%ParticleRe = 0D0
      DEALLOCATE (theParticles, re_local, re_weight)
      RETURN
    END IF

    ! Setup element type - only need function values (BDER(1))
    DO ive = 1, NNDER
      BDER(ive) = .FALSE.
    END DO
    BDER(1) = .TRUE.

    IELTYP = -1
    CALL ELE(0D0, 0D0, 0D0, IELTYP)
    IDFL = NDFL(IELTYP)

    ! Worker ranks only
    IF (myid .NE. 0) THEN

      ! Fetch particle data from PE library
      CALL getAllParticles(theParticles)

      ! ===== PARTICLE LOOP =====
      DO ip = 1, numParticles

        vel_sum = 0D0
        weight_sum = 0D0
        n_interface = 0

        ! ===== ELEMENT LOOP =====
        DO iel = 1, nel_local

          ! Get DOF mapping for this element
          CALL NDFGL(iel, 1, IELTYP, &
                     mg_mesh%level(ilev)%kvert, &
                     mg_mesh%level(ilev)%kedge, &
                     mg_mesh%level(ilev)%karea, &
                     KDFG, KDFL)

          ! ===== INTERFACE ELEMENT DETECTION =====
          ! Count DOFs inside (NIALFA) and outside (NJALFA) the particle
          NJALFA = 0
          NIALFA = 0
          DO i = 1, IDFL
            ig = KDFG(i)
            IF (ALPHA(ig) .EQ. 0) THEN
              NJALFA = NJALFA + 1  ! Outside particle (fluid)
            ELSE IF (ALPHA(ig) .EQ. ip) THEN
              NIALFA = NIALFA + 1  ! Inside particle ip
            END IF
          END DO

          ! Skip elements that are fully inside or fully outside
          ! Interface elements have: 0 < NIALFA < 27 AND 0 < NJALFA < 27
          IF (NJALFA .EQ. 27 .OR. NIALFA .EQ. 27) CYCLE

          ! This is an interface element!
          n_interface = n_interface + 1

          ! Setup element geometry for basis function evaluation
          DO ive = 1, NNVE
            ig = mg_mesh%level(ilev)%kvert(ive, iel)
            KVE(ive) = ig
            DX(ive) = mg_mesh%level(ilev)%dcorvg(1, ig)
            DY(ive) = mg_mesh%level(ilev)%dcorvg(2, ig)
            DZ(ive) = mg_mesh%level(ilev)%dcorvg(3, ig)
          END DO

          ! Evaluate basis functions at element center (0,0,0) in reference coords
          CALL ELE(0D0, 0D0, 0D0, 0)

          ! Interpolate velocity at element center
          vel_sample = 0D0
          DO il = 1, IDFL
            ig = KDFG(il)
            ive = KDFL(il)
            vel_sample(1) = vel_sample(1) + U1(ig)*DBAS(1, ive, 1)
            vel_sample(2) = vel_sample(2) + U2(ig)*DBAS(1, ive, 1)
            vel_sample(3) = vel_sample(3) + U3(ig)*DBAS(1, ive, 1)
          END DO

          ! Accumulate velocity contribution
          vel_sum = vel_sum + vel_sample
          weight_sum = weight_sum + 1D0

        END DO  ! elements

        ! ===== COMPUTE REYNOLDS NUMBER FOR THIS PARTICLE =====
        IF (weight_sum .GT. 0D0 .AND. n_interface .GT. 0) THEN

          ! Average velocity from all interface elements
          vel_avg = vel_sum/weight_sum

          ! Compute slip velocity
          slip = vel_avg - theParticles(ip)%velocity
          speed = SQRT(slip(1)**2 + slip(2)**2 + slip(3)**2)

          ! Get diameter from particle AABB (largest extent)
          diameter = MAXVAL(theParticles(ip)%aabb) 

          ! Use average viscosity (could improve by weighted average from interface elements)
          mu_avg = SUM(DVISC(1:nel_local))/DBLE(nel_local)

          ! Compute Reynolds number: Re = ρ * U * D / μ
          IF (mu_avg .GT. 0D0 .AND. diameter .GT. 0D0) THEN
            re_local(ip) = RHOFLUID*speed*diameter/mu_avg
          ELSE
            re_local(ip) = 0D0
          END IF
          re_weight(ip) = 1D0

        ELSE
          ! No interface elements found for this particle
          re_local(ip) = 0D0
          re_weight(ip) = 0D0
        END IF

      END DO  ! particles

    END IF  ! myid /= 0

    ! MPI reduction across all ranks
    CALL COMM_SUMMN(re_local, numParticles)
    CALL COMM_SUMMN(re_weight, numParticles)

    ! Finalize Reynolds numbers and store in myFBM array
    DO ip = 1, numParticles
      IF (re_weight(ip) .GT. 0D0) THEN
        myFBM%ParticleRe(ip) = re_local(ip)/re_weight(ip)
      ELSE
        myFBM%ParticleRe(ip) = 0D0
      END IF
    END DO

    ! Output logging
    IF (myid .EQ. 1) THEN
      IF (numParticles .GT. 0) THEN
        WRITE (mfile, '(A,2E16.8)') 'Particle Re (interface method): ', &
          MINVAL(myFBM%ParticleRe(1:numParticles)), MAXVAL(myFBM%ParticleRe(1:numParticles))
        WRITE (*, '(A,2E16.8)') 'Particle Re (interface method): ', &
          MINVAL(myFBM%ParticleRe(1:numParticles)), MAXVAL(myFBM%ParticleRe(1:numParticles))
      END IF
    END IF

    ! Cleanup
    DEALLOCATE (theParticles)
    DEALLOCATE (re_local)
    DEALLOCATE (re_weight)

  END SUBROUTINE fbm_compute_particle_reynolds_interface
#else
  ! Stub implementation when PE library is not available
  SUBROUTINE fbm_compute_particle_reynolds_interface(U1, U2, U3, ALPHA, DVISC, RHOFLUID, &
                                                     mfile, ELE)
    REAL*8, DIMENSION(:), INTENT(IN) :: U1, U2, U3
    INTEGER, DIMENSION(:), INTENT(IN) :: ALPHA
    REAL*8, DIMENSION(:), INTENT(IN) :: DVISC
    REAL*8, INTENT(IN) :: RHOFLUID
    INTEGER, INTENT(IN) :: mfile
    EXTERNAL :: ELE
    ! No-op: PE library not available
    IF (ALLOCATED(myFBM%ParticleRe)) myFBM%ParticleRe = 0D0
    RETURN
  END SUBROUTINE fbm_compute_particle_reynolds_interface
#endif

!=========================================================================
! fbm_compute_particle_reynolds_interface_extended - Extended interface method
!
! Enhanced version of interface-based Reynolds calculation with second-layer
! sampling for improved accuracy in DNS simulations.
!
! Improvements over standard interface method:
!   - Samples velocity from interface + second-layer elements
!   - Uses distance-based weighting (closer elements weighted more)
!   - Filters second-layer to fluid-dominated elements only
!   - O(1) element classification using integer array
!
! Method: Two-pass algorithm
!   Pass 1: Classify elements as inside/interface using ALPHA field
!   Pass 2: Expand from interface via face neighbors (kadj)
!
! Input:
!   U1, U2, U3  - Velocity field components (Q2)
!   ALPHA       - Particle indicator field (0 for fluid, particle ID for solid)
!   DVISC       - Element-wise dynamic viscosity
!   RHOFLUID    - Fluid density
!   mfile       - Output file handle
!   ELE         - Element evaluation function
!   n_layers    - Number of layers to sample (optional, default=2)
!
! Output:
!   myFBM%ParticleRe(:) - Reynolds number for each particle
!=========================================================================
#ifdef HAVE_PE
  SUBROUTINE fbm_compute_particle_reynolds_interface_extended(U1, U2, U3, ALPHA, DVISC, &
                                                              RHOFLUID, mfile, ELE, n_layers)
    use dem_query, only: numTotalParticles, getAllParticles, tParticleData
    USE PP3D_MPI, ONLY: COMM_Maximumn
    INTEGER, PARAMETER :: NNBAS = 27, NNDER = 10, NNVE = 8, NNDIM = 3

    REAL*8, DIMENSION(:), INTENT(IN) :: U1, U2, U3
    INTEGER, DIMENSION(:), INTENT(IN) :: ALPHA
    REAL*8, DIMENSION(:), INTENT(IN) :: DVISC
    REAL*8, INTENT(IN) :: RHOFLUID
    INTEGER, INTENT(IN) :: mfile
    EXTERNAL :: ELE
    INTEGER, INTENT(IN), OPTIONAL :: n_layers

    INTEGER :: ip, ig, il, IELTYP, ilev, nel_local, numParticles
    INTEGER :: NJALFA, NIALFA, n_interface, n_second_layer, i, iface, neighbor
    INTEGER :: num_layers, NJALFA_neighbor
    REAL*8 :: vel_sum(3), vel_sample(3), vel_avg(3), slip(3)
    REAL*8 :: speed, diameter, mu_avg, weight_sum, weight
    REAL*8 :: elem_center(3), dist_to_particle
    REAL*8 :: dbuf1(1)

    TYPE(tParticleData), DIMENSION(:), ALLOCATABLE :: theParticles
    REAL*8, ALLOCATABLE :: re_local(:), re_weight(:)
    INTEGER, ALLOCATABLE :: elem_class(:)  ! 0=unclassified, 1=inside, 2=interface, 3=second_layer

    INTEGER KDFG(NNBAS), KDFL(NNBAS)
    REAL*8 DX(NNVE), DY(NNVE), DZ(NNVE), DJAC(3, 3), DETJ
    REAL*8 DBAS(NNDIM, NNBAS, NNDER)
    LOGICAL BDER(NNDER)
    INTEGER KVE(NNVE), NDIM, IEL
    INTEGER :: ive

    COMMON/ELEM/DX, DY, DZ, DJAC, DETJ, DBAS, BDER, KVE, IEL, NDIM
    COMMON/COAUX1/KDFG, KDFL, IDFL
    INTEGER :: IDFL

    EXTERNAL :: NDFGL, SETLEV
    INTEGER :: NDFL

    ! Set number of layers (default = 2: interface + second layer)
    IF (PRESENT(n_layers)) THEN
      num_layers = n_layers
    ELSE
      num_layers = 2
    END IF

    ! Get total particle count (ALL ranks need this for MPI)
    numParticles = numTotalParticles()
    dbuf1(1) = DBLE(numParticles)
    CALL COMM_Maximumn(dbuf1, 1)
    numParticles = INT(dbuf1(1))

    ! Early return if no particles
    IF (numParticles .EQ. 0) RETURN

    ! Allocate or reallocate ParticleRe array if needed
    IF (.NOT. ALLOCATED(myFBM%ParticleRe)) THEN
      ALLOCATE (myFBM%ParticleRe(numParticles))
      myFBM%ParticleRe = 0D0
    ELSE IF (SIZE(myFBM%ParticleRe) .NE. numParticles) THEN
      DEALLOCATE (myFBM%ParticleRe)
      ALLOCATE (myFBM%ParticleRe(numParticles))
      myFBM%ParticleRe = 0D0
    END IF

    ! Allocate particle data and result arrays
    ALLOCATE (theParticles(numParticles))
    ALLOCATE (re_local(numParticles))
    ALLOCATE (re_weight(numParticles))
    re_local = 0D0
    re_weight = 0D0

    ilev = mg_mesh%nlmax
    CALL SETLEV(2)
    nel_local = mg_mesh%level(ilev)%nel
    IF (SIZE(DVISC) .LT. nel_local) THEN
      myFBM%ParticleRe = 0D0
      DEALLOCATE (theParticles, re_local, re_weight)
      RETURN
    END IF

    ! Allocate element classification array
    ALLOCATE (elem_class(nel_local))

    ! Setup element type
    DO ive = 1, NNDER
      BDER(ive) = .FALSE.
    END DO
    BDER(1) = .TRUE.

    IELTYP = -1
    CALL ELE(0D0, 0D0, 0D0, IELTYP)
    IDFL = NDFL(IELTYP)

    ! Worker ranks only
    IF (myid .NE. 0) THEN

      ! Fetch particle data
      CALL getAllParticles(theParticles)

      ! ===== PARTICLE LOOP =====
      DO ip = 1, numParticles

        ! Initialize element classification
        elem_class = 0

        n_interface = 0
        n_second_layer = 0

        ! ===== PASS 1: CLASSIFY ELEMENTS =====
        DO iel = 1, nel_local

          ! Get DOF mapping for this element
          CALL NDFGL(iel, 1, IELTYP, &
                     mg_mesh%level(ilev)%kvert, &
                     mg_mesh%level(ilev)%kedge, &
                     mg_mesh%level(ilev)%karea, &
                     KDFG, KDFL)

          ! Count DOFs inside and outside the particle
          NJALFA = 0  ! Outside particle (fluid)
          NIALFA = 0  ! Inside particle ip
          DO i = 1, IDFL
            ig = KDFG(i)
            IF (ALPHA(ig) .EQ. 0) THEN
              NJALFA = NJALFA + 1
            ELSE IF (ALPHA(ig) .EQ. ip) THEN
              NIALFA = NIALFA + 1
            END IF
          END DO

          ! Classify element
          IF (NJALFA .EQ. 0 .AND. NIALFA .GT. 0) THEN
            elem_class(iel) = 1  ! Fully inside particle
          ELSE IF (NJALFA .GT. 0 .AND. NIALFA .GT. 0) THEN
            elem_class(iel) = 2  ! Interface element
            n_interface = n_interface + 1
          END IF

        END DO  ! elements (Pass 1)

        ! ===== PASS 2: EXPAND TO SECOND LAYER =====
        IF (num_layers .GE. 2) THEN
          DO iel = 1, nel_local

            ! For each interface element, check face neighbors
            IF (elem_class(iel) .EQ. 2) THEN

              DO iface = 1, 6  ! 6 faces of hexahedron
                neighbor = mg_mesh%level(ilev)%kadj(iface, iel)

                ! Check if neighbor exists, is not on boundary, and not yet classified
                IF (neighbor .GT. 0 .AND. neighbor .LE. nel_local) THEN
                  IF (elem_class(neighbor) .EQ. 0) THEN

                    ! Check if neighbor is fluid-dominated
                    CALL NDFGL(neighbor, 1, IELTYP, &
                               mg_mesh%level(ilev)%kvert, &
                               mg_mesh%level(ilev)%kedge, &
                               mg_mesh%level(ilev)%karea, &
                               KDFG, KDFL)

                    NJALFA_neighbor = 0
                    DO i = 1, IDFL
                      ig = KDFG(i)
                      IF (ALPHA(ig) .EQ. 0) NJALFA_neighbor = NJALFA_neighbor + 1
                    END DO

                    ! Only add if mostly fluid (threshold: 20 out of 27 DOFs)
                    IF (NJALFA_neighbor .GE. 20) THEN
                      elem_class(neighbor) = 3  ! Second layer
                      n_second_layer = n_second_layer + 1
                    END IF

                  END IF
                END IF

              END DO  ! faces

            END IF

          END DO  ! elements (Pass 2)
        END IF  ! num_layers >= 2

        ! ===== SAMPLE VELOCITY FROM CLASSIFIED ELEMENTS =====
        vel_sum = 0D0
        weight_sum = 0D0

        DO iel = 1, nel_local

          ! Sample from interface (2) and second-layer (3) elements
          IF (elem_class(iel) .EQ. 2 .OR. elem_class(iel) .EQ. 3) THEN

            ! Setup element geometry
            DO ive = 1, NNVE
              ig = mg_mesh%level(ilev)%kvert(ive, iel)
              KVE(ive) = ig
              DX(ive) = mg_mesh%level(ilev)%dcorvg(1, ig)
              DY(ive) = mg_mesh%level(ilev)%dcorvg(2, ig)
              DZ(ive) = mg_mesh%level(ilev)%dcorvg(3, ig)
            END DO

            ! Compute element center
            elem_center(1) = SUM(DX(1:NNVE))/DBLE(NNVE)
            elem_center(2) = SUM(DY(1:NNVE))/DBLE(NNVE)
            elem_center(3) = SUM(DZ(1:NNVE))/DBLE(NNVE)

            ! Distance from element center to particle center
            dist_to_particle = SQRT((elem_center(1) - theParticles(ip)%position(1))**2 + &
                                    (elem_center(2) - theParticles(ip)%position(2))**2 + &
                                    (elem_center(3) - theParticles(ip)%position(3))**2)

            ! Distance-based weight (closer elements have more influence)
            weight = 1D0/(1D0 + dist_to_particle)

            ! Evaluate basis functions at element center
            CALL ELE(0D0, 0D0, 0D0, 0)

            ! Get DOF mapping
            CALL NDFGL(iel, 1, IELTYP, &
                       mg_mesh%level(ilev)%kvert, &
                       mg_mesh%level(ilev)%kedge, &
                       mg_mesh%level(ilev)%karea, &
                       KDFG, KDFL)

            ! Interpolate velocity at element center
            vel_sample = 0D0
            DO il = 1, IDFL
              ig = KDFG(il)
              ive = KDFL(il)
              vel_sample(1) = vel_sample(1) + U1(ig)*DBAS(1, ive, 1)
              vel_sample(2) = vel_sample(2) + U2(ig)*DBAS(1, ive, 1)
              vel_sample(3) = vel_sample(3) + U3(ig)*DBAS(1, ive, 1)
            END DO

            ! Accumulate weighted velocity
            vel_sum = vel_sum + weight*vel_sample
            weight_sum = weight_sum + weight

          END IF

        END DO  ! elements (sampling)

        ! ===== COMPUTE REYNOLDS NUMBER =====
        IF (weight_sum .GT. 0D0) THEN

          ! Weighted average velocity
          vel_avg = vel_sum/weight_sum

          ! Slip velocity
          slip = vel_avg - theParticles(ip)%velocity
          speed = SQRT(slip(1)**2 + slip(2)**2 + slip(3)**2)

          ! Get diameter from particle AABB (largest extent)
          diameter = MAXVAL(theParticles(ip)%aabb)

          ! Average viscosity
          mu_avg = SUM(DVISC(1:nel_local))/DBLE(nel_local)

          ! Compute Reynolds number
          IF (mu_avg .GT. 0D0 .AND. diameter .GT. 0D0) THEN
            re_local(ip) = RHOFLUID*speed*diameter/mu_avg
          ELSE
            re_local(ip) = 0D0
          END IF
          re_weight(ip) = 1D0

        ELSE
          re_local(ip) = 0D0
          re_weight(ip) = 0D0
        END IF

      END DO  ! particles

    END IF  ! myid /= 0

    ! MPI reduction across all ranks
    CALL COMM_SUMMN(re_local, numParticles)
    CALL COMM_SUMMN(re_weight, numParticles)

    ! Finalize Reynolds numbers
    DO ip = 1, numParticles
      IF (re_weight(ip) .GT. 0D0) THEN
        myFBM%ParticleRe(ip) = re_local(ip)/re_weight(ip)
      ELSE
        myFBM%ParticleRe(ip) = 0D0
      END IF
    END DO

    ! Output logging
    IF (myid .EQ. 1) THEN
      IF (numParticles .GT. 0) THEN
        WRITE (mfile, '(A,2E16.8)') 'Particle Re (extended interface, 2-layer): ', &
          MINVAL(myFBM%ParticleRe(1:numParticles)), MAXVAL(myFBM%ParticleRe(1:numParticles))
        WRITE (*, '(A,2E16.8)') 'Particle Re (extended interface, 2-layer): ', &
          MINVAL(myFBM%ParticleRe(1:numParticles)), MAXVAL(myFBM%ParticleRe(1:numParticles))
      END IF
    END IF

    ! Cleanup
    DEALLOCATE (theParticles)
    DEALLOCATE (re_local)
    DEALLOCATE (re_weight)
    DEALLOCATE (elem_class)

  END SUBROUTINE fbm_compute_particle_reynolds_interface_extended
#else
  ! Stub implementation when PE library is not available
  SUBROUTINE fbm_compute_particle_reynolds_interface_extended(U1, U2, U3, ALPHA, DVISC, &
                                                              RHOFLUID, mfile, ELE, n_layers)
    REAL*8, DIMENSION(:), INTENT(IN) :: U1, U2, U3
    INTEGER, DIMENSION(:), INTENT(IN) :: ALPHA
    REAL*8, DIMENSION(:), INTENT(IN) :: DVISC
    REAL*8, INTENT(IN) :: RHOFLUID
    INTEGER, INTENT(IN) :: mfile
    EXTERNAL :: ELE
    INTEGER, INTENT(IN), OPTIONAL :: n_layers
    ! No-op: PE library not available
    IF (ALLOCATED(myFBM%ParticleRe)) myFBM%ParticleRe = 0D0
    RETURN
  END SUBROUTINE fbm_compute_particle_reynolds_interface_extended
#endif

!=========================================================================
! fbm_compute_particle_reynolds_farfield - Far-field Reynolds calculation
!
! Simplified far-field Reynolds number estimation by sampling velocity at
! a point upstream of the particle. This provides a "free-stream" reference
! velocity for comparison with near-field methods.
!
! Method: For each particle, find the mesh vertex closest to the target point
!         (x_particle - diameter, y_particle, z_particle) and use its velocity
!         DOF value as the far-field velocity.
!
! Note: This is a simplified implementation for testing purposes, not optimized
!       for general production use. The sampling point is always upstream in
!       the x-direction.
!
! Input:
!   U1, U2, U3  - Velocity field components (Q2)
!   DVISC       - Element-wise dynamic viscosity
!   RHOFLUID    - Fluid density
!   mfile       - Output file handle
!
! Output:
!   myFBM%ParticleRe(:) - Far-field Reynolds number for each particle
!=========================================================================
#ifdef HAVE_PE
  SUBROUTINE fbm_compute_particle_reynolds_farfield(U1, U2, U3, DVISC, RHOFLUID, mfile)
    USE PP3D_MPI, ONLY: COMM_Maximumn
    use dem_query, only: numTotalParticles, getAllParticles, tParticleData

    REAL*8, DIMENSION(:), INTENT(IN) :: U1, U2, U3
    REAL*8, DIMENSION(:), INTENT(IN) :: DVISC
    REAL*8, INTENT(IN) :: RHOFLUID
    INTEGER, INTENT(IN) :: mfile

    INTEGER :: ip, ig, ilev, nel_local, numParticles, nvt_local
    INTEGER :: closest_vertex
    REAL*8 :: target_point(3), dist, min_dist, diameter
    REAL*8 :: vel_farfield(3), slip(3), speed, mu_avg
    REAL*8 :: dbuf1(1)

    TYPE(tParticleData), DIMENSION(:), ALLOCATABLE :: theParticles
    REAL*8, ALLOCATABLE :: re_local(:), re_weight(:)

    EXTERNAL :: SETLEV

    ! Get total particle count
    numParticles = numTotalParticles()
    dbuf1(1) = DBLE(numParticles)
    CALL COMM_Maximumn(dbuf1, 1)
    numParticles = INT(dbuf1(1))

    ! Early return if no particles
    IF (numParticles .EQ. 0) RETURN

    ! Allocate or reallocate ParticleRe array
    IF (.NOT. ALLOCATED(myFBM%ParticleRe)) THEN
      ALLOCATE (myFBM%ParticleRe(numParticles))
      myFBM%ParticleRe = 0D0
    ELSE IF (SIZE(myFBM%ParticleRe) .NE. numParticles) THEN
      DEALLOCATE (myFBM%ParticleRe)
      ALLOCATE (myFBM%ParticleRe(numParticles))
      myFBM%ParticleRe = 0D0
    END IF

    ! Allocate particle data and result arrays
    ALLOCATE (theParticles(numParticles))
    ALLOCATE (re_local(numParticles))
    ALLOCATE (re_weight(numParticles))
    re_local = 0D0
    re_weight = 0D0

    ilev = mg_mesh%nlmax
    CALL SETLEV(2)
    nel_local = mg_mesh%level(ilev)%nel
    nvt_local = mg_mesh%level(ilev)%nvt

    IF (SIZE(DVISC) .LT. nel_local) THEN
      myFBM%ParticleRe = 0D0
      DEALLOCATE (theParticles, re_local, re_weight)
      RETURN
    END IF

    ! NOTE ON VALIDITY:
    !   This implementation is calibrated for the single-particle sedimentation
    !   benchmark: the "far-field" sample is taken one particle diameter upstream
    !   in the negative x-direction, viscosity is assumed constant, and the
    !   nearest mesh vertex is accepted without checking the surrounding phase.
    !   For other flow configurations (multiple particles, arbitrary inflow
    !   direction, non-uniform viscosity) the caller must extend the sampling
    !   logic to (a) configure the offset direction/magnitude, and (b) verify
    !   that the sampled location lies in fluid. Without these checks the
    !   reported Reynolds numbers can be misleading.

    ! Worker ranks only
    IF (myid .NE. 0) THEN

      ! Fetch particle data
      CALL getAllParticles(theParticles)

      ! ===== PARTICLE LOOP =====
      DO ip = 1, numParticles

        ! Get diameter from particle AABB (largest extent)
        diameter = MAXVAL(theParticles(ip)%aabb)

        ! Target point: one diameter upstream in x-direction
        target_point(1) = theParticles(ip)%position(1) - diameter
        target_point(2) = theParticles(ip)%position(2)
        target_point(3) = theParticles(ip)%position(3)

        ! Find closest mesh vertex to target point
        min_dist = 1D10
        closest_vertex = -1

        DO ig = 1, nvt_local
          dist = SQRT((mg_mesh%level(ilev)%dcorvg(1, ig) - target_point(1))**2 + &
                      (mg_mesh%level(ilev)%dcorvg(2, ig) - target_point(2))**2 + &
                      (mg_mesh%level(ilev)%dcorvg(3, ig) - target_point(3))**2)

          IF (dist .LT. min_dist) THEN
            min_dist = dist
            closest_vertex = ig
          END IF
        END DO

        ! If we found a vertex, evaluate velocity there
        IF (closest_vertex .GT. 0) THEN

          ! Get velocity at closest vertex (Q2 DOF)
          vel_farfield(1) = U1(closest_vertex)
          vel_farfield(2) = U2(closest_vertex)
          vel_farfield(3) = U3(closest_vertex)

          ! Compute slip velocity
          slip = vel_farfield - theParticles(ip)%velocity
          speed = SQRT(slip(1)**2 + slip(2)**2 + slip(3)**2)

          ! Average viscosity
          mu_avg = SUM(DVISC(1:nel_local))/DBLE(nel_local)

          ! Compute Reynolds number
          IF (mu_avg .GT. 0D0 .AND. diameter .GT. 0D0) THEN
            re_local(ip) = RHOFLUID*speed*diameter/mu_avg
          ELSE
            re_local(ip) = 0D0
          END IF
          re_weight(ip) = 1D0

        ELSE
          re_local(ip) = 0D0
          re_weight(ip) = 0D0
        END IF

      END DO  ! particles

    END IF  ! myid /= 0

    ! MPI reduction across all ranks
    CALL COMM_SUMMN(re_local, numParticles)
    CALL COMM_SUMMN(re_weight, numParticles)

    ! Finalize Reynolds numbers
    DO ip = 1, numParticles
      IF (re_weight(ip) .GT. 0D0) THEN
        myFBM%ParticleRe(ip) = re_local(ip)/re_weight(ip)
      ELSE
        myFBM%ParticleRe(ip) = 0D0
      END IF
    END DO

    ! Output logging
    IF (myid .EQ. 1) THEN
      IF (numParticles .GT. 0) THEN
        WRITE (mfile, '(A,2E16.8)') 'Particle Re (far-field, 1-diameter upstream): ', &
          MINVAL(myFBM%ParticleRe(1:numParticles)), MAXVAL(myFBM%ParticleRe(1:numParticles))
        WRITE (*, '(A,2E16.8)') 'Particle Re (far-field, 1-diameter upstream): ', &
          MINVAL(myFBM%ParticleRe(1:numParticles)), MAXVAL(myFBM%ParticleRe(1:numParticles))
      END IF
    END IF

    ! Cleanup
    DEALLOCATE (theParticles)
    DEALLOCATE (re_local)
    DEALLOCATE (re_weight)

  END SUBROUTINE fbm_compute_particle_reynolds_farfield
#else
  ! Stub implementation when PE library is not available
  SUBROUTINE fbm_compute_particle_reynolds_farfield(U1, U2, U3, DVISC, RHOFLUID, mfile)
    REAL*8, DIMENSION(:), INTENT(IN) :: U1, U2, U3
    REAL*8, DIMENSION(:), INTENT(IN) :: DVISC
    REAL*8, INTENT(IN) :: RHOFLUID
    INTEGER, INTENT(IN) :: mfile
    ! No-op: PE library not available
    IF (ALLOCATED(myFBM%ParticleRe)) myFBM%ParticleRe = 0D0
    RETURN
  END SUBROUTINE fbm_compute_particle_reynolds_farfield
#endif

!=========================================================================
! fbm_compute_particle_reynolds_el - Euler-Lagrange particle Reynolds number
!
! Computes particle Reynolds numbers for Euler-Lagrange simulations where
! particles are treated as subgrid point particles embedded in a volume-
! averaged fluid field. Designed for semi-dense suspensions (volume
! fractions of order 5-10%) where drag correlations (Schiller-Naumann,
! Wen-Yu, Gidaspow, Di Felice, etc.) require Re_p and the local void
! fraction epsilon_f.
!
! Standard definition:
!   Re_p = rho_f * |u_f - v_p| * d_p / mu
!
! The void fraction does NOT enter Re_p itself; it is stored separately
! for the drag correlation to use (e.g. Wen-Yu multiplies C_D by
! epsilon_f^{-2.65}).
!
! Method:
!   For each particle, locate the containing mesh element via bounding-
!   box culling + trilinear inversion (fbmaux_PointInHex), interpolate
!   the Q2 velocity field at the particle position, and read the element-
!   wise void fraction. The slip velocity vector is stored for direct
!   use in the drag force computation.
!
! Input:
!   U1, U2, U3 - Velocity field components (Q2 DOFs)
!   DVISC       - Element-wise dynamic viscosity
!   RHOFLUID    - Fluid density
!   EPSI        - Element-wise fluid void fraction, epsilon_f = 1 - epsilon_s
!   mfile       - Output file handle
!   ELE         - Element evaluation function (e.g., E013)
!
! Output (stored in myFBM):
!   ParticleRe(:)        - Particle Reynolds number
!   ParticleSlipVel(:,:) - Slip velocity components (3, nParticles)
!   ParticleEpsilon(:)   - Local void fraction at each particle
!=========================================================================
#ifdef HAVE_PE
  SUBROUTINE fbm_compute_particle_reynolds_el(U1, U2, U3, DVISC, RHOFLUID, &
                                               EPSI, mfile, ELE)
    USE fbmaux, ONLY: fbmaux_PointInHex
    use dem_query, only: numTotalParticles, getAllParticles, tParticleData
    USE PP3D_MPI, ONLY: COMM_Maximumn
    INTEGER, PARAMETER :: NNBAS = 27, NNDER = 10, NNVE = 8, NNDIM = 3

    REAL*8, DIMENSION(:), INTENT(IN) :: U1, U2, U3
    REAL*8, DIMENSION(:), INTENT(IN) :: DVISC
    REAL*8, INTENT(IN) :: RHOFLUID
    REAL*8, DIMENSION(:), INTENT(IN) :: EPSI
    INTEGER, INTENT(IN) :: mfile
    EXTERNAL :: ELE

    INTEGER :: ip, ig, il, IELTYP, ilev, nel_local, numParticles
    REAL*8 :: xi1, xi2, xi3
    REAL*8 :: vel_sample(3), slip(3), speed, diameter, mu_loc, epsi_loc
    REAL*8 :: xverts(8), yverts(8), zverts(8)
    REAL*8 :: xmin, xmax, ymin, ymax, zmin, zmax, eps_box
    REAL*8 :: dbuf1(1)
    LOGICAL :: found

    TYPE(tParticleData), DIMENSION(:), ALLOCATABLE :: theParticles
    REAL*8, ALLOCATABLE :: re_local(:), re_weight(:)
    REAL*8, ALLOCATABLE :: slip_local(:,:), epsi_local(:)

    INTEGER KDFG(NNBAS), KDFL(NNBAS)
    REAL*8 DX(NNVE), DY(NNVE), DZ(NNVE), DJAC(3, 3), DETJ
    REAL*8 DBAS(NNDIM, NNBAS, NNDER)
    LOGICAL BDER(NNDER)
    INTEGER KVE(NNVE), NDIM, IEL
    INTEGER :: ive

    COMMON/ELEM/DX, DY, DZ, DJAC, DETJ, DBAS, BDER, KVE, IEL, NDIM
    COMMON/COAUX1/KDFG, KDFL, IDFL
    INTEGER :: IDFL

    EXTERNAL :: NDFGL, SETLEV
    INTEGER :: NDFL

    ! --- Determine global particle count (all ranks must agree) ---
    numParticles = numTotalParticles()
    dbuf1(1) = DBLE(numParticles)
    CALL COMM_Maximumn(dbuf1, 1)
    numParticles = INT(dbuf1(1))

    IF (numParticles .EQ. 0) RETURN

    ! --- (Re-)allocate result arrays in myFBM ---
    IF (.NOT. ALLOCATED(myFBM%ParticleRe)) THEN
      ALLOCATE (myFBM%ParticleRe(numParticles))
    ELSE IF (SIZE(myFBM%ParticleRe) .NE. numParticles) THEN
      DEALLOCATE (myFBM%ParticleRe)
      ALLOCATE (myFBM%ParticleRe(numParticles))
    END IF
    myFBM%ParticleRe = 0D0

    IF (.NOT. ALLOCATED(myFBM%ParticleSlipVel)) THEN
      ALLOCATE (myFBM%ParticleSlipVel(3, numParticles))
    ELSE IF (SIZE(myFBM%ParticleSlipVel, 2) .NE. numParticles) THEN
      DEALLOCATE (myFBM%ParticleSlipVel)
      ALLOCATE (myFBM%ParticleSlipVel(3, numParticles))
    END IF
    myFBM%ParticleSlipVel = 0D0

    IF (.NOT. ALLOCATED(myFBM%ParticleEpsilon)) THEN
      ALLOCATE (myFBM%ParticleEpsilon(numParticles))
    ELSE IF (SIZE(myFBM%ParticleEpsilon) .NE. numParticles) THEN
      DEALLOCATE (myFBM%ParticleEpsilon)
      ALLOCATE (myFBM%ParticleEpsilon(numParticles))
    END IF
    myFBM%ParticleEpsilon = 1D0  ! default: pure fluid

    ! --- Allocate work arrays ---
    ALLOCATE (theParticles(numParticles))
    ALLOCATE (re_local(numParticles))
    ALLOCATE (re_weight(numParticles))
    ALLOCATE (slip_local(3, numParticles))
    ALLOCATE (epsi_local(numParticles))
    re_local = 0D0
    re_weight = 0D0
    slip_local = 0D0
    epsi_local = 0D0

    ilev = mg_mesh%nlmax
    CALL SETLEV(2)
    nel_local = mg_mesh%level(ilev)%nel
    IF (SIZE(DVISC) .LT. nel_local .OR. SIZE(EPSI) .LT. nel_local) THEN
      myFBM%ParticleRe = 0D0
      myFBM%ParticleSlipVel = 0D0
      myFBM%ParticleEpsilon = 1D0
      DEALLOCATE (theParticles, re_local, re_weight, slip_local, epsi_local)
      RETURN
    END IF

    ! Setup basis function evaluation (values only, no derivatives)
    DO ive = 1, NNDER
      BDER(ive) = .FALSE.
    END DO
    BDER(1) = .TRUE.

    IELTYP = -1
    CALL ELE(0D0, 0D0, 0D0, IELTYP)
    IDFL = NDFL(IELTYP)

    eps_box = 1D-10

    IF (myid .NE. 0) THEN

      CALL getAllParticles(theParticles)

      ! ===== PARTICLE LOOP =====
      DO ip = 1, numParticles
        found = .FALSE.

        ! ===== ELEMENT SEARCH (bounding-box cull + trilinear inversion) =====
        DO iel = 1, nel_local
          DO ive = 1, NNVE
            ig = mg_mesh%level(ilev)%kvert(ive, iel)
            xverts(ive) = mg_mesh%level(ilev)%dcorvg(1, ig)
            yverts(ive) = mg_mesh%level(ilev)%dcorvg(2, ig)
            zverts(ive) = mg_mesh%level(ilev)%dcorvg(3, ig)
          END DO

          xmin = MINVAL(xverts); xmax = MAXVAL(xverts)
          ymin = MINVAL(yverts); ymax = MAXVAL(yverts)
          zmin = MINVAL(zverts); zmax = MAXVAL(zverts)

          IF (theParticles(ip)%position(1) .LT. xmin - eps_box) CYCLE
          IF (theParticles(ip)%position(1) .GT. xmax + eps_box) CYCLE
          IF (theParticles(ip)%position(2) .LT. ymin - eps_box) CYCLE
          IF (theParticles(ip)%position(2) .GT. ymax + eps_box) CYCLE
          IF (theParticles(ip)%position(3) .LT. zmin - eps_box) CYCLE
          IF (theParticles(ip)%position(3) .GT. zmax + eps_box) CYCLE

          xi1 = 0D0; xi2 = 0D0; xi3 = 0D0
          found = fbmaux_PointInHex(theParticles(ip)%position(1), &
                                    theParticles(ip)%position(2), &
                                    theParticles(ip)%position(3), &
                                    xverts, yverts, zverts, xi1, xi2, xi3, iel)
          IF (.NOT. found) CYCLE

          ! --- Evaluate Q2 basis at particle position ---
          DO ive = 1, NNDER
            BDER(ive) = .FALSE.
          END DO
          BDER(1) = .TRUE.

          CALL NDFGL(iel, 1, IELTYP, mg_mesh%level(ilev)%kvert, &
                     mg_mesh%level(ilev)%kedge, mg_mesh%level(ilev)%karea, &
                     KDFG, KDFL)

          DO ive = 1, NNVE
            ig = mg_mesh%level(ilev)%kvert(ive, iel)
            KVE(ive) = ig
            DX(ive) = xverts(ive)
            DY(ive) = yverts(ive)
            DZ(ive) = zverts(ive)
          END DO

          CALL ELE(xi1, xi2, xi3, 0)

          ! --- Interpolate fluid velocity at particle position ---
          vel_sample = 0D0
          DO il = 1, IDFL
            ig = KDFG(il)
            ive = KDFL(il)
            vel_sample(1) = vel_sample(1) + U1(ig)*DBAS(1, ive, 1)
            vel_sample(2) = vel_sample(2) + U2(ig)*DBAS(1, ive, 1)
            vel_sample(3) = vel_sample(3) + U3(ig)*DBAS(1, ive, 1)
          END DO

          ! --- Read element-wise quantities ---
          mu_loc = DVISC(iel)
          epsi_loc = EPSI(iel)

          ! --- Slip velocity ---
          slip = vel_sample - theParticles(ip)%velocity
          speed = SQRT(slip(1)**2 + slip(2)**2 + slip(3)**2)

          ! Particle diameter from AABB (largest extent)
          diameter = MAXVAL(theParticles(ip)%aabb)

          ! --- Standard particle Reynolds number ---
          ! Re_p = rho_f * |u_f - v_p| * d_p / mu
          IF (mu_loc .GT. 0D0 .AND. diameter .GT. 0D0) THEN
            re_local(ip) = RHOFLUID * speed * diameter / mu_loc
          ELSE
            re_local(ip) = 0D0
          END IF

          slip_local(:, ip) = slip
          epsi_local(ip) = epsi_loc
          re_weight(ip) = 1D0

          EXIT  ! particle found in this element, move to next particle
        END DO  ! elements
      END DO  ! particles

    END IF  ! myid /= 0

    ! --- MPI reduction ---
    CALL COMM_SUMMN(re_local, numParticles)
    CALL COMM_SUMMN(re_weight, numParticles)
    CALL COMM_SUMMN(slip_local, 3*numParticles)
    CALL COMM_SUMMN(epsi_local, numParticles)

    ! --- Store results ---
    DO ip = 1, numParticles
      IF (re_weight(ip) .GT. 0D0) THEN
        myFBM%ParticleRe(ip) = re_local(ip) / re_weight(ip)
        myFBM%ParticleSlipVel(:, ip) = slip_local(:, ip) / re_weight(ip)
        myFBM%ParticleEpsilon(ip) = epsi_local(ip) / re_weight(ip)
      ELSE
        myFBM%ParticleRe(ip) = 0D0
        myFBM%ParticleSlipVel(:, ip) = 0D0
        myFBM%ParticleEpsilon(ip) = 1D0  ! pure fluid if particle not found
      END IF
    END DO

    ! --- Logging ---
    IF (myid .EQ. 1) THEN
      IF (numParticles .GT. 0) THEN
        WRITE (mfile, '(A,2E16.8)') 'Particle Re (E-L): ', &
          MINVAL(myFBM%ParticleRe(1:numParticles)), MAXVAL(myFBM%ParticleRe(1:numParticles))
        WRITE (mfile, '(A,2E16.8)') 'Void fraction range at particles: ', &
          MINVAL(myFBM%ParticleEpsilon(1:numParticles)), MAXVAL(myFBM%ParticleEpsilon(1:numParticles))
      END IF
    END IF

    DEALLOCATE (theParticles, re_local, re_weight, slip_local, epsi_local)

  END SUBROUTINE fbm_compute_particle_reynolds_el
#else
  SUBROUTINE fbm_compute_particle_reynolds_el(U1, U2, U3, DVISC, RHOFLUID, &
                                               EPSI, mfile, ELE)
    REAL*8, DIMENSION(:), INTENT(IN) :: U1, U2, U3
    REAL*8, DIMENSION(:), INTENT(IN) :: DVISC
    REAL*8, INTENT(IN) :: RHOFLUID
    REAL*8, DIMENSION(:), INTENT(IN) :: EPSI
    INTEGER, INTENT(IN) :: mfile
    EXTERNAL :: ELE
    IF (ALLOCATED(myFBM%ParticleRe)) myFBM%ParticleRe = 0D0
    IF (ALLOCATED(myFBM%ParticleSlipVel)) myFBM%ParticleSlipVel = 0D0
    IF (ALLOCATED(myFBM%ParticleEpsilon)) myFBM%ParticleEpsilon = 1D0
    RETURN
  END SUBROUTINE fbm_compute_particle_reynolds_el
#endif

!=========================================================================
! fbm_compute_particle_stokes - Compute Stokes number from Reynolds number
!
! Computes the Stokes number for each particle based on previously
! calculated Reynolds numbers and particle/fluid densities.
!
! The Stokes number characterizes particle inertia relative to fluid flow:
!   St = (ρ_particle / ρ_fluid) * (Re_p / 18)
!
! Physical interpretation:
!   St << 1: Particle closely follows fluid streamlines
!   St ~ 1:  Particle responds to flow with some lag
!   St >> 1: Particle dominated by inertia, crosses streamlines
!
! Prerequisites:
!   - myFBM%ParticleRe must be populated (call Reynolds calculation first)
!   - Particle data must be available from PE (density via getAllParticles)
!
! Input:
!   RHOFLUID - Fluid density
!   mfile    - Output file handle for logging
!
! Output:
!   myFBM%ParticleSt(:) - Stokes number for each particle
!=========================================================================
#ifdef HAVE_PE
  SUBROUTINE fbm_compute_particle_stokes(RHOFLUID, mfile)
    USE PP3D_MPI, ONLY: COMM_Maximumn
    use dem_query, only: numTotalParticles, getAllParticles, tParticleData

    REAL*8, INTENT(IN) :: RHOFLUID
    INTEGER, INTENT(IN) :: mfile

    INTEGER :: ip, numParticles
    REAL*8 :: rho_particle
    REAL*8 :: dbuf1(1)

    TYPE(tParticleData), DIMENSION(:), ALLOCATABLE :: theParticles

    ! Get total particle count
    numParticles = numTotalParticles()
    dbuf1(1) = DBLE(numParticles)
    CALL COMM_Maximumn(dbuf1, 1)
    numParticles = INT(dbuf1(1))

    ! Early return if no particles
    IF (numParticles .EQ. 0) RETURN

    ! Check that ParticleRe exists and is populated
    IF (.NOT. ALLOCATED(myFBM%ParticleRe)) THEN
      IF (myid .EQ. 1) THEN
        WRITE (mfile, '(A)') 'WARNING: ParticleRe not allocated. Compute Reynolds number first.'
        WRITE (*, '(A)') 'WARNING: ParticleRe not allocated. Compute Reynolds number first.'
      END IF
      RETURN
    END IF

    IF (SIZE(myFBM%ParticleRe) .NE. numParticles) THEN
      IF (myid .EQ. 1) THEN
        WRITE (mfile, '(A)') 'WARNING: ParticleRe size mismatch. Compute Reynolds number first.'
        WRITE (*, '(A)') 'WARNING: ParticleRe size mismatch. Compute Reynolds number first.'
      END IF
      RETURN
    END IF

    ! Allocate or reallocate ParticleSt array
    IF (.NOT. ALLOCATED(myFBM%ParticleSt)) THEN
      ALLOCATE (myFBM%ParticleSt(numParticles))
      myFBM%ParticleSt = 0D0
    ELSE IF (SIZE(myFBM%ParticleSt) .NE. numParticles) THEN
      DEALLOCATE (myFBM%ParticleSt)
      ALLOCATE (myFBM%ParticleSt(numParticles))
      myFBM%ParticleSt = 0D0
    END IF

    ! Allocate particle data
    ALLOCATE (theParticles(numParticles))

    ! Worker ranks only
    IF (myid .NE. 0) THEN
      ! Fetch particle data from PE library
      CALL getAllParticles(theParticles)
    END IF

    ! Compute Stokes number for each particle
    ! St = (ρ_particle / ρ_fluid) * (Re_p / 18)
    DO ip = 1, numParticles
      IF (myid .NE. 0) THEN
        rho_particle = theParticles(ip)%density
      ELSE
        rho_particle = 0D0
      END IF

      IF (RHOFLUID .GT. 0D0) THEN
        myFBM%ParticleSt(ip) = (rho_particle/RHOFLUID) * (myFBM%ParticleRe(ip)/18D0)
      ELSE
        myFBM%ParticleSt(ip) = 0D0
      END IF
    END DO

    ! Output logging
    IF (myid .EQ. 1) THEN
      IF (numParticles .GT. 0) THEN
        WRITE (mfile, '(A,2E16.8)') 'Particle Stokes number range: ', &
          MINVAL(myFBM%ParticleSt(1:numParticles)), MAXVAL(myFBM%ParticleSt(1:numParticles))
        WRITE (*, '(A,2E16.8)') 'Particle Stokes number range: ', &
          MINVAL(myFBM%ParticleSt(1:numParticles)), MAXVAL(myFBM%ParticleSt(1:numParticles))
      END IF
    END IF

    ! Cleanup
    DEALLOCATE (theParticles)

  END SUBROUTINE fbm_compute_particle_stokes
#else
  ! Stub implementation when PE library is not available
  SUBROUTINE fbm_compute_particle_stokes(RHOFLUID, mfile)
    REAL*8, INTENT(IN) :: RHOFLUID
    INTEGER, INTENT(IN) :: mfile
    ! No-op: PE library not available
    IF (ALLOCATED(myFBM%ParticleSt)) myFBM%ParticleSt = 0D0
    RETURN
  END SUBROUTINE fbm_compute_particle_stokes
#endif

END MODULE fbm_particle_reynolds
