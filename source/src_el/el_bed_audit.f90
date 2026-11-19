!=========================================================================
! EL_BED_AUDIT - bed and column audit line for fluidised / packed beds (case 11)
!
! SimPar@ELBedAuditFreq = N (EL_CONFIG; default 0 = off): every N steps, after the fluid solve,
! one line on the protocol and the terminal (rank showid):
!
!   EL_BED_AUDIT t= step= p_in_out= p_in p_out area_in_out= area_in area_out
!                dp_hydro_dyn= dp dp_hydro dp_dyn n= n weight_z_floor_z= weight_z floor_z
!                z50_95_99_max= z50 z95 z99 zmax eps_bed= eps_bed window= z_lo z_hi
!
!   p_in, p_out   area-mean pressure on the bottom (z = zmin of the EL domain box) and the top
!                 (z = zmax) plane. The P1 (discontinuous linear) pressure of each element is
!                 p0 + grad p . (x - c) with c the vertex mean (the form the EL sampler uses,
!                 EL_QUADRATURE); it is integrated over each element face that lies in the
!                 plane (2 x 2 Gauss points on the bilinear face, exact for this integrand),
!                 summed over the ranks and divided by the summed face area. No smoothing, no
!                 projection to the vertices.
!   area_in/out   summed face areas (the cross-section A of the column).
!   dp            p_in - p_out.
!   dp_hydro      the hydrostatic part contained in dp: -rho_f g_z (zmax - zmin) when the fluid
!                 carries gravity (ELFluidGravity = Yes), 0 otherwise. Every momentum model
!                 implies ELFluidGravity = No (dynamic pressure, buoyancy on the particle side),
!                 so dp_hydro = 0 there and the pressure contains no hydrostatic part.
!   dp_dyn        dp - dp_hydro; in a fluidised bed dp_dyn A = buoyant bed weight.
!   n             number of particles (owned, summed over the ranks).
!   weight_z      sum over the particles of (rho_p - rho_f) V_p g_z (the buoyant weight; the
!                 EL_FORCE_BUDGET grav_buoy z with the same sign).
!   floor_z       force of the bottom z-wall of the particle engine on the particles,
!                 -(contact impulse on wall 0)_z / dt of the last macro step (pe
!                 getElWallContactImpulse; 0 without z walls).
!   z50/z95/z99   heights below which 50 / 95 / 99 % of the particle centres lie (histogram of
!                 NBINS bins over [zmin, zmax], linear inside a bin: resolution (zmax-zmin)/NBINS).
!   zmax          highest particle centre.
!   eps_bed       porosity of the bed window [z_lo, z_hi] = [zmin + 2 d, z99 + r - 2 d] from the
!                 particle positions (exact sphere-slab overlap, d and r the mean diameter and
!                 radius): 1 - V_solid(window) / (area_in (z_hi - z_lo)); -1 when the window is
!                 empty. Laterally periodic and walled boxes alike (the slab volume of a sphere
!                 does not depend on x, y).
!
! The hydrodynamic forces on the particles (drag, pressure force) are in the EL_FORCE_BUDGET
! line of the transfer pass of the same step (evaluated with the fluid state before the
! solve of the step); the post tool of case 11 joins the two by step.
!
! Collectives on MPI_COMM_SUBS: the routine is called by every worker (never by the master),
! with or without particles.
!=========================================================================
MODULE EL_BED_AUDIT

  USE MPI
  USE PP3D_MPI, ONLY: myid, showid, MPI_COMM_SUBS
  USE TYPES, ONLY: tMultiMesh
  USE EL_CONFIG, ONLY: el_bed_audit_freq, el_fluid_gravity
  USE EL_GEOMETRY, ONLY: el_domain_box, el_domain_box_set
  USE EL_QUADRATURE, ONLY: EL_Q1_MAP
#ifdef HAVE_PE
  USE DEM_QUERY, ONLY: tParticleData, numLocalParticles, getAllParticles, &
                       getElWallContactImpulse
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: EL_WRITE_BED_AUDIT
  ! helpers, public for the unit test
  PUBLIC :: EL_SPHERE_SLAB_VOLUME, EL_HISTOGRAM_QUANTILE, EL_PLANE_PRESSURE_LOCAL

  INTEGER, PARAMETER :: NBINS = 4000

CONTAINS

  SUBROUTINE EL_WRITE_BED_AUDIT(pressure, time, dt, mfile, istep, mesh, ilev, density, gravity)
    REAL*8, INTENT(IN) :: pressure(:)
    REAL*8, INTENT(IN) :: time, dt
    INTEGER, INTENT(IN) :: mfile, istep
    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    REAL*8, INTENT(IN) :: density, gravity(3)

    REAL*8 :: zmin, zmax, tol, loc(7), glob(7), p_in, p_out, dp, dp_hydro
    REAL*8 :: zq(3), zlo, zhi, eps_bed, dmean, rmean, zmax_loc, zmax_glob, floor_z, imp(3)
    REAL*8 :: vwin_loc, vwin
    INTEGER :: hist_loc(NBINS), hist(NBINS), nloc, nglob, ierr, i, ib
#ifdef HAVE_PE
    TYPE(tParticleData), ALLOCATABLE :: particles(:)
#endif

    IF (el_bed_audit_freq.LE.0) RETURN
    IF (MOD(ABS(istep), el_bed_audit_freq).NE.0) RETURN
    IF (.NOT.el_domain_box_set) RETURN
    zmin = el_domain_box(5)
    zmax = el_domain_box(6)
    tol = 1.0d-6*(zmax - zmin)

    ! loc: 1 area_in, 2 int p in, 3 area_out, 4 int p out, 5 weight_z, 6 sum r,
    !      7 floor contact impulse z
    loc = 0.0d0
    CALL EL_PLANE_PRESSURE_LOCAL(mesh, ilev, pressure, zmin, tol, loc(1), loc(2))
    CALL EL_PLANE_PRESSURE_LOCAL(mesh, ilev, pressure, zmax, tol, loc(3), loc(4))

    hist_loc = 0
    nloc = 0
    zmax_loc = -HUGE(1.0d0)
#ifdef HAVE_PE
    nloc = numLocalParticles()
    IF (nloc.GT.0) THEN
      ALLOCATE(particles(nloc))
      CALL getAllParticles(particles)
      DO i = 1, nloc
        loc(5) = loc(5) + (particles(i)%density - density)* &
                 4.0d0*ACOS(-1.0d0)*particles(i)%radius**3/3.0d0*gravity(3)
        loc(6) = loc(6) + particles(i)%radius
        ib = 1 + INT(DBLE(NBINS)*(particles(i)%position(3) - zmin)/(zmax - zmin))
        ib = MIN(NBINS, MAX(1, ib))
        hist_loc(ib) = hist_loc(ib) + 1
        zmax_loc = MAX(zmax_loc, particles(i)%position(3))
      END DO
    END IF
    imp = 0.0d0
    CALL getElWallContactImpulse(0, imp)
    loc(7) = imp(3)
#endif
    CALL MPI_Allreduce(loc, glob, 7, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(hist_loc, hist, NBINS, MPI_INTEGER, MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(nloc, nglob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(zmax_loc, zmax_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_SUBS, ierr)

    p_in = 0.0d0
    p_out = 0.0d0
    IF (glob(1).GT.0.0d0) p_in = glob(2)/glob(1)
    IF (glob(3).GT.0.0d0) p_out = glob(4)/glob(3)
    dp = p_in - p_out
    dp_hydro = 0.0d0
    IF (el_fluid_gravity) dp_hydro = -density*gravity(3)*(zmax - zmin)
    floor_z = 0.0d0
    IF (dt.GT.0.0d0) floor_z = -glob(7)/dt

    ! bed heights and the porosity of the bed window (identical histogram on every rank)
    zq = zmin
    eps_bed = -1.0d0
    zlo = zmin
    zhi = zmin
    IF (nglob.GT.0) THEN
      zq(1) = EL_HISTOGRAM_QUANTILE(hist, zmin, zmax, 0.50d0)
      zq(2) = EL_HISTOGRAM_QUANTILE(hist, zmin, zmax, 0.95d0)
      zq(3) = EL_HISTOGRAM_QUANTILE(hist, zmin, zmax, 0.99d0)
      rmean = glob(6)/DBLE(nglob)
      dmean = 2.0d0*rmean
      zlo = zmin + 2.0d0*dmean
      zhi = zq(3) + rmean - 2.0d0*dmean
      vwin_loc = 0.0d0
#ifdef HAVE_PE
      IF (zhi.GT.zlo) THEN
        DO i = 1, nloc
          vwin_loc = vwin_loc + EL_SPHERE_SLAB_VOLUME(particles(i)%position(3), &
                                                      particles(i)%radius, zlo, zhi)
        END DO
      END IF
#endif
      CALL MPI_Allreduce(vwin_loc, vwin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
      IF (zhi.GT.zlo .AND. glob(1).GT.0.0d0) eps_bed = 1.0d0 - vwin/(glob(1)*(zhi - zlo))
    END IF
#ifdef HAVE_PE
    IF (ALLOCATED(particles)) DEALLOCATE(particles)
#endif

    IF (myid.EQ.showid) THEN
      CALL WRITE_LINE(6)
      CALL WRITE_LINE(mfile)
    END IF

  CONTAINS

    SUBROUTINE WRITE_LINE(unit)
      INTEGER, INTENT(IN) :: unit
      WRITE(unit,'(A,ES14.6,A,I0,A,2ES16.8,A,2ES14.6,A,3ES16.8,A,I0,A,2ES14.6,A,4ES14.6,' // &
                 'A,ES14.6,A,2ES14.6)') &
        'EL_BED_AUDIT t= ', time, ' step= ', istep, ' p_in_out= ', p_in, p_out, &
        ' area_in_out= ', glob(1), glob(3), ' dp_hydro_dyn= ', dp, dp_hydro, dp - dp_hydro, &
        ' n= ', nglob, ' weight_z_floor_z= ', glob(5), floor_z, &
        ' z50_95_99_max= ', zq, zmax_glob, ' eps_bed= ', eps_bed, ' window= ', zlo, zhi
    END SUBROUTINE WRITE_LINE

  END SUBROUTINE EL_WRITE_BED_AUDIT

  ! Area and area integral of the P1 pressure over the element faces of this rank that lie in
  ! the plane z = zplane (all four face vertices within tol). Pressure layout of the solver:
  ! 4 values per element (value at the vertex mean c, then the gradient), p(x) = p0 +
  ! grad p . (x - c).
  SUBROUTINE EL_PLANE_PRESSURE_LOCAL(mesh, ilev, pressure, zplane, tol, area, p_integral)
    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    REAL*8, INTENT(IN) :: pressure(:), zplane, tol
    REAL*8, INTENT(INOUT) :: area, p_integral

    REAL*8, PARAMETER :: GP = 0.57735026918962576d0
    REAL*8 :: x(8), y(8), z(8), center(3), xi(3), point(3), jac(3,3), detj
    REAL*8 :: ta(3), tb(3), cn(3), da, pval
    INTEGER :: iel, i, ig, nface, k, sf, kf, sface, ia, ib, ga, gb, jj
    LOGICAL :: found

    DO iel = 1, mesh%level(ilev)%nel
      nface = 0
      center = 0.0d0
      DO i = 1, 8
        ig = mesh%level(ilev)%kvert(i,iel)
        x(i) = mesh%level(ilev)%dcorvg(1,ig)
        y(i) = mesh%level(ilev)%dcorvg(2,ig)
        z(i) = mesh%level(ilev)%dcorvg(3,ig)
        center = center + 0.125d0*mesh%level(ilev)%dcorvg(:,ig)
        IF (ABS(z(i) - zplane).LT.tol) nface = nface + 1
      END DO
      IF (nface.LT.4) CYCLE

      ! the reference face xi(kf) = sface whose centre maps onto the plane
      found = .FALSE.
      kf = 0
      sface = 0
      DO k = 1, 3
        DO sf = -1, 1, 2
          xi = 0.0d0
          xi(k) = DBLE(sf)
          CALL EL_Q1_MAP(x, y, z, xi, point, jac, detj)
          IF (ABS(point(3) - zplane).LT.tol .AND. .NOT.found) THEN
            kf = k
            sface = sf
            found = .TRUE.
          END IF
        END DO
      END DO
      IF (.NOT.found) CYCLE
      ia = MOD(kf, 3) + 1
      ib = MOD(kf + 1, 3) + 1

      jj = 4*(iel - 1) + 1
      DO ga = -1, 1, 2
        DO gb = -1, 1, 2
          xi(kf) = DBLE(sface)
          xi(ia) = DBLE(ga)*GP
          xi(ib) = DBLE(gb)*GP
          CALL EL_Q1_MAP(x, y, z, xi, point, jac, detj)
          ta = jac(:,ia)
          tb = jac(:,ib)
          cn(1) = ta(2)*tb(3) - ta(3)*tb(2)
          cn(2) = ta(3)*tb(1) - ta(1)*tb(3)
          cn(3) = ta(1)*tb(2) - ta(2)*tb(1)
          da = SQRT(SUM(cn*cn))
          pval = pressure(jj) + DOT_PRODUCT(point - center, pressure(jj+1:jj+3))
          area = area + da
          p_integral = p_integral + da*pval
        END DO
      END DO
    END DO
  END SUBROUTINE EL_PLANE_PRESSURE_LOCAL

  ! Volume of the part of a sphere (centre height zc, radius r) between the planes za < zb.
  PURE REAL*8 FUNCTION EL_SPHERE_SLAB_VOLUME(zc, r, za, zb)
    REAL*8, INTENT(IN) :: zc, r, za, zb
    REAL*8 :: lo, hi
    lo = MAX(za - zc, -r)
    hi = MIN(zb - zc, r)
    IF (hi.LE.lo) THEN
      EL_SPHERE_SLAB_VOLUME = 0.0d0
    ELSE
      ! integral of pi (r^2 - s^2) ds from lo to hi
      EL_SPHERE_SLAB_VOLUME = ACOS(-1.0d0)*(r*r*(hi - lo) - (hi**3 - lo**3)/3.0d0)
    END IF
  END FUNCTION EL_SPHERE_SLAB_VOLUME

  ! Height below which the fraction q of the histogram counts lies (bins of equal width over
  ! [z0, z1]; linear inside the bin that crosses q).
  PURE REAL*8 FUNCTION EL_HISTOGRAM_QUANTILE(hist, z0, z1, q)
    INTEGER, INTENT(IN) :: hist(:)
    REAL*8, INTENT(IN) :: z0, z1, q
    REAL*8 :: total, target, cum, width
    INTEGER :: ib

    total = DBLE(SUM(hist))
    width = (z1 - z0)/DBLE(SIZE(hist))
    EL_HISTOGRAM_QUANTILE = z0
    IF (total.LE.0.0d0) RETURN
    target = q*total
    cum = 0.0d0
    DO ib = 1, SIZE(hist)
      IF (hist(ib).GT.0 .AND. cum + DBLE(hist(ib)).GE.target) THEN
        EL_HISTOGRAM_QUANTILE = z0 + width*(DBLE(ib - 1) + (target - cum)/DBLE(hist(ib)))
        RETURN
      END IF
      cum = cum + DBLE(hist(ib))
    END DO
    EL_HISTOGRAM_QUANTILE = z1
  END FUNCTION EL_HISTOGRAM_QUANTILE

END MODULE EL_BED_AUDIT
