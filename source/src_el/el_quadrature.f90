MODULE EL_QUADRATURE

  USE MPI
  USE TYPES, ONLY: tMultiMesh
  USE EL_HALO, ONLY: tElParticleRecord
  USE EL_KERNEL_FUNCTIONS, ONLY: EL_KERNEL_VALUE
  USE EL_FIELDS, ONLY: tElFieldStorage
  USE EL_CONFIG, ONLY: el_kernel_width_factor
  USE PP3D_MPI, ONLY: myid, MPI_COMM_SUBS

  IMPLICIT NONE

  INTEGER, PARAMETER :: EL_SAMPLE_SIZE = 18
  INTEGER, PARAMETER :: EL_SAMPLE_NORMALIZATION = 1
  INTEGER, PARAMETER :: EL_SAMPLE_UF_BEGIN = 2
  INTEGER, PARAMETER :: EL_SAMPLE_UF_END = 4
  INTEGER, PARAMETER :: EL_SAMPLE_GRAD_U_BEGIN = 5
  INTEGER, PARAMETER :: EL_SAMPLE_GRAD_U_END = 13
  INTEGER, PARAMETER :: EL_SAMPLE_PRESSURE = 14
  INTEGER, PARAMETER :: EL_SAMPLE_GRAD_P_BEGIN = 15
  INTEGER, PARAMETER :: EL_SAMPLE_GRAD_P_END = 17
  INTEGER, PARAMETER :: EL_SAMPLE_EPSILON_F = 18
  INTEGER, PARAMETER :: NNBAS = 27
  INTEGER, PARAMETER :: NNDER = 10
  INTEGER, PARAMETER :: NNCUBP = 36
  INTEGER, PARAMETER :: NNVE = 8

CONTAINS

  SUBROUTINE EL_INTEGRATE_PARTICLE(mesh, ilev, velocity_u, velocity_v, velocity_w, &
      pressure, epsilon_f, particle, sample)

    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    REAL*8, INTENT(IN) :: velocity_u(:), velocity_v(:), velocity_w(:)
    REAL*8, INTENT(IN) :: pressure(:), epsilon_f(:)
    TYPE(tElParticleRecord), INTENT(IN) :: particle
    REAL*8, INTENT(OUT) :: sample(EL_SAMPLE_SIZE)

    INTEGER :: iel, icubp, i, ig, jj
    INTEGER :: kdfl(NNBAS), kdfg(NNBAS), ieltyp, idfl
    REAL*8 :: xi(3), point(3), jac(3,3), detj, omega, kernel_weight
    REAL*8 :: basis, grad_basis(3), velocity(3), grad_velocity(3,3)
    REAL*8 :: pressure_value, pressure_grad(3), center(3), width
    LOGICAL :: supported

    REAL*8 :: dx(NNVE), dy(NNVE), dz(NNVE), djac(3,3), detj_common
    REAL*8 :: dbas(3,NNBAS,NNDER), dxi(NNCUBP,3), domega(NNCUBP)
    LOGICAL :: bder(NNDER)
    INTEGER :: kve(NNVE), iel_common, ndim, idfl_common
    INTEGER :: nel_common, nvt_common, net_common, nat_common
    INTEGER :: nve_common, nee_common, nae_common, nvel_common, neel_common
    INTEGER :: nved_common, nvar_common, near_common, nbct_common
    INTEGER :: nvbd_common, nebd_common, nabd_common, ncubp
    INTEGER :: ier, icheck

    COMMON /ERRCTL/ ier,icheck
    COMMON /ELEM/ dx,dy,dz,djac,detj_common,dbas,bder,kve,iel_common,ndim
    COMMON /TRIAD/ nel_common,nvt_common,net_common,nat_common,nve_common, &
      nee_common,nae_common,nvel_common,neel_common,nved_common,nvar_common, &
      near_common,nbct_common,nvbd_common,nebd_common,nabd_common
    COMMON /CUB/ dxi,domega,ncubp,icubp
    COMMON /COAUX1/ kdfg,kdfl,idfl_common

    INTEGER, EXTERNAL :: NDFL
    EXTERNAL E013

    sample = 0.0d0
    width = 2.0d0*el_kernel_width_factor*particle%radius
    IF (width.LE.0.0d0) RETURN

    bder = .FALSE.
    bder(1:4) = .TRUE.
    ieltyp = -1
    CALL E013(0.0d0,0.0d0,0.0d0,ieltyp)
    idfl = NDFL(ieltyp)
    idfl_common = idfl
    CALL CB3H(9)
    IF (ier.NE.0) RETURN

    DO iel=1,mesh%level(ilev)%nel
      supported = EL_ELEMENT_INTERSECTS_SUPPORT(mesh, ilev, iel, &
        particle%position, width)
      IF (.NOT.supported) CYCLE

      iel_common = iel
      center = 0.0d0
      DO i=1,NNVE
        ig = mesh%level(ilev)%kvert(i,iel)
        kve(i) = ig
        dx(i) = mesh%level(ilev)%dcorvg(1,ig)
        dy(i) = mesh%level(ilev)%dcorvg(2,ig)
        dz(i) = mesh%level(ilev)%dcorvg(3,ig)
        center = center + 0.125d0*mesh%level(ilev)%dcorvg(:,ig)
      END DO
      CALL NDFGL(iel,1,ieltyp,mesh%level(ilev)%kvert, &
        mesh%level(ilev)%kedge,mesh%level(ilev)%karea,kdfg,kdfl)
      IF (ier.LT.0) RETURN
      CALL EL_ASSERT_KDFG_IN_BOUNDS('EL_INTEGRATE_PARTICLE', ilev, iel, &
        ieltyp, MIN(SIZE(velocity_u), SIZE(velocity_v), SIZE(velocity_w)), &
        kdfg, idfl)
      CALL E013(0.0d0,0.0d0,0.0d0,-2)

      DO icubp=1,ncubp
        xi = dxi(icubp,:)
        CALL EL_Q1_MAP(dx,dy,dz,xi,point,jac,detj)
        djac = jac
        detj_common = detj
        omega = domega(icubp)*ABS(detj)
        kernel_weight = EL_KERNEL_VALUE(SQRT(SUM((point-particle%position)**2)), &
                                        width)
        IF (kernel_weight.LE.0.0d0) CYCLE

        CALL E013(xi(1),xi(2),xi(3),-3)
        IF (ier.LT.0) RETURN
        velocity = 0.0d0
        grad_velocity = 0.0d0
        DO i=1,idfl
          ig = kdfg(i)
          basis = dbas(1,kdfl(i),1)
          grad_basis = dbas(1,kdfl(i),2:4)
          velocity(1) = velocity(1) + velocity_u(ig)*basis
          velocity(2) = velocity(2) + velocity_v(ig)*basis
          velocity(3) = velocity(3) + velocity_w(ig)*basis
          grad_velocity(1,:) = grad_velocity(1,:) + velocity_u(ig)*grad_basis
          grad_velocity(2,:) = grad_velocity(2,:) + velocity_v(ig)*grad_basis
          grad_velocity(3,:) = grad_velocity(3,:) + velocity_w(ig)*grad_basis
        END DO

        jj = 4*(iel-1)+1
        pressure_grad = pressure(jj+1:jj+3)
        pressure_value = pressure(jj) + DOT_PRODUCT(point-center,pressure_grad)

        sample(1) = sample(1) + kernel_weight*omega
        sample(2:4) = sample(2:4) + kernel_weight*omega*velocity
        sample(5:7) = sample(5:7) + kernel_weight*omega*grad_velocity(1,:)
        sample(8:10) = sample(8:10) + kernel_weight*omega*grad_velocity(2,:)
        sample(11:13) = sample(11:13) + kernel_weight*omega*grad_velocity(3,:)
        sample(14) = sample(14) + kernel_weight*omega*pressure_value
        sample(15:17) = sample(15:17) + kernel_weight*omega*pressure_grad
        sample(18) = sample(18) + kernel_weight*omega*epsilon_f(iel)
      END DO
    END DO

  END SUBROUTINE EL_INTEGRATE_PARTICLE

  !----------------------------------------------------------------------------
  ! EL_INTEGRATE_FLUID_MOMENTUM
  !
  ! Physical element-cubature integral of the fluid momentum over all LOCALLY
  ! OWNED elements (each cell counted exactly once), using the same Q2 (E013)
  ! basis/quadrature as the solver. This avoids the partition/periodic shared-DOF
  ! double counting of the lumped  Sum(MlRhoPmat*u)  diagnostic.
  !   momentum     = integral rho * u            dV
  !   momentum_eps = integral rho * epsilon_f * u dV   (void-fraction weighted)
  ! Local (per-rank) result; the caller MPI_Allreduce's over MPI_COMM_SUBS.
  !----------------------------------------------------------------------------
  SUBROUTINE EL_INTEGRATE_FLUID_MOMENTUM(mesh, ilev, velocity_u, velocity_v, &
      velocity_w, epsilon_f, density, momentum, momentum_eps)

    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    REAL*8, INTENT(IN) :: velocity_u(:), velocity_v(:), velocity_w(:)
    REAL*8, INTENT(IN) :: epsilon_f(:), density
    REAL*8, INTENT(OUT) :: momentum(3), momentum_eps(3)

    INTEGER :: iel, icubp, i, ig, ieltyp, idfl
    INTEGER :: kdfl(NNBAS), kdfg(NNBAS)
    REAL*8 :: xi(3), point(3), jac(3,3), detj, omega, basis, velocity(3)

    REAL*8 :: dx(NNVE), dy(NNVE), dz(NNVE), djac(3,3), detj_common
    REAL*8 :: dbas(3,NNBAS,NNDER), dxi(NNCUBP,3), domega(NNCUBP)
    LOGICAL :: bder(NNDER)
    INTEGER :: kve(NNVE), iel_common, ndim, idfl_common
    INTEGER :: nel_common, nvt_common, net_common, nat_common
    INTEGER :: nve_common, nee_common, nae_common, nvel_common, neel_common
    INTEGER :: nved_common, nvar_common, near_common, nbct_common
    INTEGER :: nvbd_common, nebd_common, nabd_common, ncubp
    INTEGER :: ier, icheck

    COMMON /ERRCTL/ ier,icheck
    COMMON /ELEM/ dx,dy,dz,djac,detj_common,dbas,bder,kve,iel_common,ndim
    COMMON /TRIAD/ nel_common,nvt_common,net_common,nat_common,nve_common, &
      nee_common,nae_common,nvel_common,neel_common,nved_common,nvar_common, &
      near_common,nbct_common,nvbd_common,nebd_common,nabd_common
    COMMON /CUB/ dxi,domega,ncubp,icubp
    COMMON /COAUX1/ kdfg,kdfl,idfl_common

    INTEGER, EXTERNAL :: NDFL
    EXTERNAL E013

    momentum = 0.0d0
    momentum_eps = 0.0d0

    bder = .FALSE.
    bder(1) = .TRUE.
    ieltyp = -1
    CALL E013(0.0d0,0.0d0,0.0d0,ieltyp)
    idfl = NDFL(ieltyp)
    idfl_common = idfl
    CALL CB3H(9)
    IF (ier.NE.0) RETURN

    DO iel=1,mesh%level(ilev)%nel
      iel_common = iel
      DO i=1,NNVE
        ig = mesh%level(ilev)%kvert(i,iel)
        kve(i) = ig
        dx(i) = mesh%level(ilev)%dcorvg(1,ig)
        dy(i) = mesh%level(ilev)%dcorvg(2,ig)
        dz(i) = mesh%level(ilev)%dcorvg(3,ig)
      END DO
      CALL NDFGL(iel,1,ieltyp,mesh%level(ilev)%kvert, &
        mesh%level(ilev)%kedge,mesh%level(ilev)%karea,kdfg,kdfl)
      IF (ier.LT.0) RETURN
      CALL EL_ASSERT_KDFG_IN_BOUNDS('EL_INTEGRATE_FLUID_MOMENTUM', ilev, &
        iel, ieltyp, MIN(SIZE(velocity_u), SIZE(velocity_v), &
        SIZE(velocity_w)), kdfg, idfl)
      CALL E013(0.0d0,0.0d0,0.0d0,-2)

      DO icubp=1,ncubp
        xi = dxi(icubp,:)
        CALL EL_Q1_MAP(dx,dy,dz,xi,point,jac,detj)
        djac = jac
        detj_common = detj
        omega = domega(icubp)*ABS(detj)
        CALL E013(xi(1),xi(2),xi(3),-3)
        IF (ier.LT.0) RETURN
        velocity = 0.0d0
        DO i=1,idfl
          ig = kdfg(i)
          basis = dbas(1,kdfl(i),1)
          velocity(1) = velocity(1) + velocity_u(ig)*basis
          velocity(2) = velocity(2) + velocity_v(ig)*basis
          velocity(3) = velocity(3) + velocity_w(ig)*basis
        END DO
        momentum = momentum + density*omega*velocity
        momentum_eps = momentum_eps + density*epsilon_f(iel)*omega*velocity
      END DO
    END DO

  END SUBROUTINE EL_INTEGRATE_FLUID_MOMENTUM

  !----------------------------------------------------------------------------
  ! EL_INTEGRATE_UZ_LANES
  !
  ! Element-cubature accumulation of  integral u_z dV  and  integral dV  into
  ! nbins uniform slabs in x and (independently) in y, over all LOCALLY OWNED
  ! elements. Feeds the mesoscale lane-mode filter: the caller Allreduce's the
  ! bins over MPI_COMM_SUBS and forms the zero-mean slab-average profiles
  ! ubar_z(x), ubar_z(y). Same E013 basis/quadrature as the momentum audit.
  !----------------------------------------------------------------------------
  SUBROUTINE EL_INTEGRATE_UZ_LANES(mesh, ilev, velocity_w, nbins, &
      xmin, xlen, ymin, ylen, num_x, den_x, num_y, den_y)

    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev, nbins
    REAL*8, INTENT(IN) :: velocity_w(:)
    REAL*8, INTENT(IN) :: xmin, xlen, ymin, ylen
    REAL*8, INTENT(OUT) :: num_x(nbins), den_x(nbins)
    REAL*8, INTENT(OUT) :: num_y(nbins), den_y(nbins)

    INTEGER :: iel, icubp, i, ig, ieltyp, idfl, ibx, iby
    INTEGER :: kdfl(NNBAS), kdfg(NNBAS)
    REAL*8 :: xi(3), point(3), jac(3,3), detj, omega, wvel

    REAL*8 :: dx(NNVE), dy(NNVE), dz(NNVE), djac(3,3), detj_common
    REAL*8 :: dbas(3,NNBAS,NNDER), dxi(NNCUBP,3), domega(NNCUBP)
    LOGICAL :: bder(NNDER)
    INTEGER :: kve(NNVE), iel_common, ndim, idfl_common
    INTEGER :: nel_common, nvt_common, net_common, nat_common
    INTEGER :: nve_common, nee_common, nae_common, nvel_common, neel_common
    INTEGER :: nved_common, nvar_common, near_common, nbct_common
    INTEGER :: nvbd_common, nebd_common, nabd_common, ncubp
    INTEGER :: ier, icheck

    COMMON /ERRCTL/ ier,icheck
    COMMON /ELEM/ dx,dy,dz,djac,detj_common,dbas,bder,kve,iel_common,ndim
    COMMON /TRIAD/ nel_common,nvt_common,net_common,nat_common,nve_common, &
      nee_common,nae_common,nvel_common,neel_common,nved_common,nvar_common, &
      near_common,nbct_common,nvbd_common,nebd_common,nabd_common
    COMMON /CUB/ dxi,domega,ncubp,icubp
    COMMON /COAUX1/ kdfg,kdfl,idfl_common

    INTEGER, EXTERNAL :: NDFL
    EXTERNAL E013

    num_x = 0.0d0
    den_x = 0.0d0
    num_y = 0.0d0
    den_y = 0.0d0
    IF (xlen.LE.0.0d0 .OR. ylen.LE.0.0d0) RETURN

    bder = .FALSE.
    bder(1) = .TRUE.
    ieltyp = -1
    CALL E013(0.0d0,0.0d0,0.0d0,ieltyp)
    idfl = NDFL(ieltyp)
    idfl_common = idfl
    CALL CB3H(9)
    IF (ier.NE.0) RETURN

    DO iel=1,mesh%level(ilev)%nel
      iel_common = iel
      DO i=1,NNVE
        ig = mesh%level(ilev)%kvert(i,iel)
        kve(i) = ig
        dx(i) = mesh%level(ilev)%dcorvg(1,ig)
        dy(i) = mesh%level(ilev)%dcorvg(2,ig)
        dz(i) = mesh%level(ilev)%dcorvg(3,ig)
      END DO
      CALL NDFGL(iel,1,ieltyp,mesh%level(ilev)%kvert, &
        mesh%level(ilev)%kedge,mesh%level(ilev)%karea,kdfg,kdfl)
      IF (ier.LT.0) RETURN
      CALL EL_ASSERT_KDFG_IN_BOUNDS('EL_INTEGRATE_UZ_LANES', ilev, &
        iel, ieltyp, SIZE(velocity_w), kdfg, idfl)
      CALL E013(0.0d0,0.0d0,0.0d0,-2)

      DO icubp=1,ncubp
        xi = dxi(icubp,:)
        CALL EL_Q1_MAP(dx,dy,dz,xi,point,jac,detj)
        djac = jac
        detj_common = detj
        omega = domega(icubp)*ABS(detj)
        CALL E013(xi(1),xi(2),xi(3),-3)
        IF (ier.LT.0) RETURN
        wvel = 0.0d0
        DO i=1,idfl
          wvel = wvel + velocity_w(kdfg(i))*dbas(1,kdfl(i),1)
        END DO
        ibx = MIN(nbins, MAX(1, 1 + INT((point(1)-xmin)/xlen*DBLE(nbins))))
        iby = MIN(nbins, MAX(1, 1 + INT((point(2)-ymin)/ylen*DBLE(nbins))))
        num_x(ibx) = num_x(ibx) + omega*wvel
        den_x(ibx) = den_x(ibx) + omega
        num_y(iby) = num_y(iby) + omega*wvel
        den_y(iby) = den_y(iby) + omega
      END DO
    END DO

  END SUBROUTINE EL_INTEGRATE_UZ_LANES

  SUBROUTINE EL_DEPOSIT_PARTICLE(mesh, ilev, particle, normalization, force, fields, &
      drag_b)

    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    TYPE(tElParticleRecord), INTENT(IN) :: particle
    REAL*8, INTENT(IN) :: normalization, force(3)
    TYPE(tElFieldStorage), INTENT(INOUT) :: fields
    REAL*8, INTENT(IN), OPTIONAL :: drag_b

    INTEGER :: iel, icubp, i, ig
    INTEGER :: kdfl(NNBAS), kdfg(NNBAS), ieltyp, idfl
    REAL*8 :: xi(3), point(3), jac(3,3), detj, omega, kernel_weight
    REAL*8 :: width, particle_volume, element_integral, basis_integral(NNBAS)
    LOGICAL :: supported

    REAL*8 :: dx(NNVE), dy(NNVE), dz(NNVE), djac(3,3), detj_common
    REAL*8 :: dbas(3,NNBAS,NNDER), dxi(NNCUBP,3), domega(NNCUBP)
    LOGICAL :: bder(NNDER)
    INTEGER :: kve(NNVE), iel_common, ndim, idfl_common
    INTEGER :: nel_common, nvt_common, net_common, nat_common
    INTEGER :: nve_common, nee_common, nae_common, nvel_common, neel_common
    INTEGER :: nved_common, nvar_common, near_common, nbct_common
    INTEGER :: nvbd_common, nebd_common, nabd_common, ncubp
    INTEGER :: ier, icheck

    COMMON /ERRCTL/ ier,icheck
    COMMON /ELEM/ dx,dy,dz,djac,detj_common,dbas,bder,kve,iel_common,ndim
    COMMON /TRIAD/ nel_common,nvt_common,net_common,nat_common,nve_common, &
      nee_common,nae_common,nvel_common,neel_common,nved_common,nvar_common, &
      near_common,nbct_common,nvbd_common,nebd_common,nabd_common
    COMMON /CUB/ dxi,domega,ncubp,icubp
    COMMON /COAUX1/ kdfg,kdfl,idfl_common

    INTEGER, EXTERNAL :: NDFL
    EXTERNAL E013

    IF (normalization.LE.0.0d0) RETURN
    width = 2.0d0*el_kernel_width_factor*particle%radius
    particle_volume = 4.0d0*ACOS(-1.0d0)*particle%radius**3/3.0d0

    bder = .FALSE.
    bder(1) = .TRUE.
    ieltyp = -1
    CALL E013(0.0d0,0.0d0,0.0d0,ieltyp)
    idfl = NDFL(ieltyp)
    idfl_common = idfl
    CALL CB3H(9)
    IF (ier.NE.0) RETURN

    DO iel=1,mesh%level(ilev)%nel
      supported = EL_ELEMENT_INTERSECTS_SUPPORT(mesh, ilev, iel, &
        particle%position, width)
      IF (.NOT.supported) CYCLE

      iel_common = iel
      DO i=1,NNVE
        ig = mesh%level(ilev)%kvert(i,iel)
        kve(i) = ig
        dx(i) = mesh%level(ilev)%dcorvg(1,ig)
        dy(i) = mesh%level(ilev)%dcorvg(2,ig)
        dz(i) = mesh%level(ilev)%dcorvg(3,ig)
      END DO
      CALL NDFGL(iel,1,ieltyp,mesh%level(ilev)%kvert, &
        mesh%level(ilev)%kedge,mesh%level(ilev)%karea,kdfg,kdfl)
      IF (ier.LT.0) RETURN
      CALL EL_ASSERT_KDFG_IN_BOUNDS('EL_DEPOSIT_PARTICLE', ilev, iel, &
        ieltyp, SIZE(fields%force_rhs,2), kdfg, idfl)

      CALL E013(0.0d0,0.0d0,0.0d0,-2)

      element_integral = 0.0d0
      basis_integral = 0.0d0
      DO icubp=1,ncubp
        xi = dxi(icubp,:)
        CALL EL_Q1_MAP(dx,dy,dz,xi,point,jac,detj)
        djac = jac
        detj_common = detj
        omega = domega(icubp)*ABS(detj)
        kernel_weight = EL_KERNEL_VALUE(SQRT(SUM((point-particle%position)**2)), &
                                        width)
        IF (kernel_weight.LE.0.0d0) CYCLE
        element_integral = element_integral + kernel_weight*omega
        CALL E013(xi(1),xi(2),xi(3),-3)
        IF (ier.LT.0) RETURN
        DO i=1,idfl
          basis_integral(i) = basis_integral(i) + &
            dbas(1,kdfl(i),1)*kernel_weight*omega
        END DO
      END DO

      fields%alpha_p(iel) = fields%alpha_p(iel) + &
        particle_volume*element_integral/(normalization*mesh%level(ilev)%dvol(iel))
      DO i=1,idfl
        ig = kdfg(i)
        fields%force_rhs(:,ig) = fields%force_rhs(:,ig) - &
          force*basis_integral(i)/normalization
        IF (PRESENT(drag_b)) THEN
          fields%drag_B_ml(ig) = fields%drag_B_ml(ig) + &
            drag_b*basis_integral(i)/normalization
        END IF
      END DO
    END DO

  END SUBROUTINE EL_DEPOSIT_PARTICLE

  SUBROUTINE EL_ASSERT_KDFG_IN_BOUNDS(context, ilev, iel, ieltyp, ndof, kdfg, idfl)

    CHARACTER(*), INTENT(IN) :: context
    INTEGER, INTENT(IN) :: ilev, iel, ieltyp, ndof, idfl
    INTEGER, INTENT(IN) :: kdfg(:)

    INTEGER :: ierr, kmin, kmax

    IF (idfl.LE.0) THEN
      WRITE(*,*) 'Fatal EL NDFGL bounds: rank=', myid, &
        ' context=', TRIM(context), ' ilev=', ilev, ' iel=', iel, &
        ' ieltyp=', ieltyp, ' ndof=', ndof, ' idfl=', idfl
      CALL MPI_Abort(MPI_COMM_SUBS, 905, ierr)
    END IF

    kmin = MINVAL(kdfg(1:idfl))
    kmax = MAXVAL(kdfg(1:idfl))
    IF (kmin.LT.1 .OR. kmax.GT.ndof) THEN
      WRITE(*,*) 'Fatal EL NDFGL bounds: rank=', myid, &
        ' context=', TRIM(context), ' ilev=', ilev, ' iel=', iel, &
        ' ieltyp=', ieltyp, ' ndof=', ndof, ' min_kdfg=', kmin, &
        ' max_kdfg=', kmax
      CALL MPI_Abort(MPI_COMM_SUBS, 905, ierr)
    END IF

  END SUBROUTINE EL_ASSERT_KDFG_IN_BOUNDS

  LOGICAL FUNCTION EL_ELEMENT_INTERSECTS_SUPPORT(mesh, ilev, iel, position, radius)

    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev, iel
    REAL*8, INTENT(IN) :: position(3), radius
    REAL*8 :: lower(3), upper(3), gap(3)
    INTEGER :: i, vertex

    lower = HUGE(1.0d0)
    upper = -HUGE(1.0d0)
    DO i=1,NNVE
      vertex = mesh%level(ilev)%kvert(i,iel)
      lower = MIN(lower,mesh%level(ilev)%dcorvg(:,vertex))
      upper = MAX(upper,mesh%level(ilev)%dcorvg(:,vertex))
    END DO
    gap = MAX(0.0d0,MAX(lower-position,position-upper))
    EL_ELEMENT_INTERSECTS_SUPPORT = SUM(gap*gap).LE.radius*radius

  END FUNCTION EL_ELEMENT_INTERSECTS_SUPPORT

  ! Geometric isoparametric map and Jacobian for the cubature points.
  ! NOTE: this uses the trilinear (Q1) hex map even though the fluid fields
  ! live in the triquadratic Q2 (E013) space. That is exact for the affine /
  ! straight-edged hexes produced by the structured mesh hierarchy used here
  ! (the geometry is genuinely Q1), so the cubature weights omega = |detJ|*w
  ! and the kernel-volume conservation are correct. If curved/Q2-geometry
  ! hexes are ever introduced, this map must be upgraded to the full Q2
  ! isoparametric Jacobian to preserve conservation.
  SUBROUTINE EL_Q1_MAP(x, y, z, xi, point, jac, detj)

    REAL*8, INTENT(IN) :: x(8), y(8), z(8), xi(3)
    REAL*8, INTENT(OUT) :: point(3), jac(3,3), detj
    REAL*8, PARAMETER :: signs(3,8) = RESHAPE((/ &
      -1.0d0,-1.0d0,-1.0d0,  1.0d0,-1.0d0,-1.0d0, &
       1.0d0, 1.0d0,-1.0d0, -1.0d0, 1.0d0,-1.0d0, &
      -1.0d0,-1.0d0, 1.0d0,  1.0d0,-1.0d0, 1.0d0, &
       1.0d0, 1.0d0, 1.0d0, -1.0d0, 1.0d0, 1.0d0 /), (/3,8/))
    REAL*8 :: shape, derivative(3), coordinates(3)
    INTEGER :: i

    point = 0.0d0
    jac = 0.0d0
    DO i=1,8
      coordinates = (/x(i),y(i),z(i)/)
      shape = 0.125d0*(1.0d0+signs(1,i)*xi(1))* &
        (1.0d0+signs(2,i)*xi(2))*(1.0d0+signs(3,i)*xi(3))
      derivative(1) = 0.125d0*signs(1,i)*(1.0d0+signs(2,i)*xi(2))* &
        (1.0d0+signs(3,i)*xi(3))
      derivative(2) = 0.125d0*signs(2,i)*(1.0d0+signs(1,i)*xi(1))* &
        (1.0d0+signs(3,i)*xi(3))
      derivative(3) = 0.125d0*signs(3,i)*(1.0d0+signs(1,i)*xi(1))* &
        (1.0d0+signs(2,i)*xi(2))
      point = point + shape*coordinates
      jac(:,1) = jac(:,1) + derivative(1)*coordinates
      jac(:,2) = jac(:,2) + derivative(2)*coordinates
      jac(:,3) = jac(:,3) + derivative(3)*coordinates
    END DO
    detj = jac(1,1)*(jac(2,2)*jac(3,3)-jac(3,2)*jac(2,3)) - &
      jac(2,1)*(jac(1,2)*jac(3,3)-jac(3,2)*jac(1,3)) + &
      jac(3,1)*(jac(1,2)*jac(2,3)-jac(2,2)*jac(1,3))

  END SUBROUTINE EL_Q1_MAP

END MODULE EL_QUADRATURE
