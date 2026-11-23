!=========================================================================
! EL_BED_AUDIT - bed and column audit line for fluidised / packed beds (case 11)
!
! SimPar@ELBedAuditFreq = N (EL_CONFIG; default 0 = off): every N steps, after the fluid solve,
! one line on the protocol and the terminal (rank showid):
!
!   EL_BED_AUDIT t= step= p_in_out= p_in p_out area_in_out= area_in area_out
!                dp_hydro_dyn= dp dp_hydro dp_dyn n= n weight_z_floor_z= weight_z floor_z
!                z50_95_99_max= z50 z95 z99 zmax eps_bed= eps_bed window= z_lo z_hi
!                gradz_win= G floor_bal_z= F eps_kern_slip_win= e s win_margin= m
!                gradz_win_dyn= Gd
!
! Unavailable values are printed as NaN (never 0 or -HUGE): a plane without an element face
! (p_in / p_out and with them dp, dp_dyn), no particle (heights, window, eps_bed, gradz_win),
! an empty window, a floor balance outside its assumptions (floor_bal_z), no particle centre
! in the window (eps_kern_slip_win).
!
! Geometry: the audit needs the box domain (SimPar@ELDomainType = box); EL_CONFIG refuses a
! positive ELBedAuditFreq with any other domain at start-up on every rank.
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
!   eps_bed       porosity of the bed window [z_lo, z_hi] = [zmin + m, z99 + r - m] from the
!                 particle positions (exact sphere-slab overlap, d and r the mean diameter and
!                 radius, m = win_margin): 1 - V_solid(window) / (area_in (z_hi - z_lo)); NaN when
!                 the window is empty. Laterally periodic and walled boxes alike (the slab volume
!                 of a sphere does not depend on x, y).
!   win_margin    m = MAX(2 d, kernel support radius of the mean particle) (EL_CONFIG:
!                 EL_KERNEL_SUPPORT_RADIUS, the radius the transfer maps with): the window
!                 excludes two diameters and the mapping support at both ends.
!   gradz_win     -<dp/dz> over the same window: the element-wise constant z gradient of the P1
!                 pressure weighted with the part of the element volume inside [z_lo, z_hi]
!                 (overlap of the element z range with the window). EXACT ONLY for elements of
!                 constant horizontal cross-section between two z planes (prisms along z, the
!                 case-11 mesh); for any other element, even a tapered hex with horizontal top and
!                 bottom faces, the height fraction is not the volume fraction and the value is an
!                 estimate. The total gradient, including the
!                 hydrostatic part when the fluid carries gravity (gate F2 uses gradz_win_dyn).
!   gradz_win_dyn the window gradient without the hydrostatic part, gradz_win - dp_hydro /
!                 (zmax - zmin) (= gradz_win when the fluid carries no gravity, as in every
!                 momentum model): the counterpart of dp_dyn.
!   floor_bal_z   force of the walls on the particles from the momentum balance of the particles
!                 over the last PE advance: (dp - expected)_z / dt with dp the measured change of
!                 the total particle momentum and expected the impulse of the forces the CFD side
!                 applied (EL_DIAGNOSTICS: EL_NEWTON_PAIR_END, kept every step). It is the floor
!                 force when the floor is the only wall the particles touch (case 11: the top wall
!                 is far above the bed). Independent of the contact bookkeeping of pe; floor_z
!                 (pe wall-contact impulse) reports the unrelaxed PGS impulse, 1/0.9 of the
!                 impulse pe applies (design note section 8). Assumptions: the CFD side arms pe
!                 with the hydrodynamic forces and the buoyant weight (ELApplyParticleForces =
!                 Yes; expected is then the impulse of exactly those forces) and pe adds no
!                 external force of its own (no pe gravity). The value then contains the force of
!                 ALL walls. NaN when the balance is not available: forces not armed
!                 (ELApplyParticleForces = No, nothing in expected acted on the particles) or no
!                 balance recorded for this step.
!   eps_kern_slip_win  mean kernel-sampled porosity (the value the drag closure uses, clipped as
!                 there) and mean sampled slip u_f,z - v_z of the particles whose centre lies in
!                 the window (bin centre in [z_lo, z_hi]), from the post-advance transfer pass of
!                 this step (EL_BED_AUDIT_RECORD_SAMPLE); NaN / NaN without particles there.
!
! The hydrodynamic forces on the particles (drag, pressure force) are in the EL_FORCE_BUDGET
! line of the transfer pass of the same step (evaluated with the fluid state before the
! solve of the step); the post tool of case 11 joins the two by step.
!
! Window audit (SimPar@ELBedWindowAudit = Yes, EL_CONFIG; experiment E2 of case 11), written on the
! same steps after the audit line, by EL_WRITE_BED_WINDOWS (called through the fluid adapter, which
! supplies the vertex-mode sums). One line per window of the table WINDOWS below:
!
!   EL_BED_WINDOW t= step= win= NAME z_lo_hi= zl zh gradz= G n= N eps_geo= e eps_kern= e
!                 slip= sx sy sz drag_mean= Fd drag_of_means= Fm pres_mean= Fp visc_mean= Fv
!
!   windows      'm2d' / 'm2.5d': [zmin + m, z99 + r - m] with m = 2 d / 2.5 d (d, r the mean
!                diameter and radius, z99 of the audit line), 'interior': [zmin + 5 mm, zmin + 20 mm]
!   gradz        -<dp/dz> of the broken P1 pressure over the window (EL_WINDOW_GRADIENT_LOCAL)
!   n            particles with the centre in the window (bin centre, as eps_kern_slip_win)
!   eps_geo      1 - V_solid(window) / (area_in (zh - zl)) from the particle positions
!   eps_kern, slip   means over those particles of the kernel-sampled porosity (clipped as in the
!                closure) and of the sampled slip u_f - v (post-advance transfer pass)
!   drag_mean    mean z closure drag PER PARTICLE [N], <B(eps_i, |s_i|) s_z,i> (EL_FORCE_BUDGET drag)
!   drag_of_means  z drag of one particle of mean radius and density at the mean porosity and the
!                mean slip vector [N], B(<eps>, |<s>|) <s_z>
!   pres_mean, visc_mean   mean z pressure / viscous force per particle [N]
!
!   EL_BED_VERTEX_MODE t= step= z_lo_hi= zl zh vertex= n mean std min max edge= ... face= ...
!                 centre= ...
!   w (m/s) over the Q2 dofs of each class with zmin + 5 mm <= z <= zmin + 20 mm, every dof counted
!   once (weight 1 / number of copies); n is the dof count of the class in the window.
!
! Collectives on MPI_COMM_SUBS: the routine is called by every worker (never by the master),
! with or without particles.
!=========================================================================
MODULE EL_BED_AUDIT

  USE MPI
  USE PP3D_MPI, ONLY: myid, showid, MPI_COMM_SUBS
  USE TYPES, ONLY: tMultiMesh
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_VALUE, IEEE_QUIET_NAN
  USE EL_CONFIG, ONLY: el_bed_audit_freq, el_fluid_gravity, el_apply_particle_forces, &
                       EL_KERNEL_SUPPORT_RADIUS, el_bed_window_audit
  USE EL_FORCES, ONLY: EL_DRAG_CLOSURE
  USE EL_MOMENTUM_MODEL, ONLY: EL_MM_DRAG_SCALE
  USE EL_GEOMETRY, ONLY: el_domain_box, el_domain_box_set
  USE EL_QUADRATURE, ONLY: EL_Q1_MAP
  USE EL_DIAGNOSTICS, ONLY: el_pair_last_step, el_pair_last_dp, el_pair_last_expected
#ifdef HAVE_PE
  USE DEM_QUERY, ONLY: tParticleData, numLocalParticles, getAllParticles, &
                       getElWallContactImpulse
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: EL_WRITE_BED_AUDIT, EL_BED_AUDIT_BEGIN_SAMPLES, EL_BED_AUDIT_RECORD_SAMPLE
  PUBLIC :: EL_BED_AUDIT_RECORD_FORCES, EL_WRITE_BED_WINDOWS, EL_BED_WINDOW_BOUNDS
  PUBLIC :: NWIN, BW_VM_LO, BW_VM_HI
  ! helpers, public for the unit test
  PUBLIC :: EL_SPHERE_SLAB_VOLUME, EL_HISTOGRAM_QUANTILE, EL_PLANE_PRESSURE_LOCAL, &
            EL_WINDOW_GRADIENT_LOCAL, EL_BED_WINDOW_MARGIN, EL_BED_FLOOR_BALANCE, &
            EL_BED_WINDOW_SAMPLE_SUMS, EL_BED_UNAVAILABLE

  INTEGER, PARAMETER :: NBINS = 4000
  ! per height bin of the particle centre: sum of the kernel-sampled porosity, of the sampled slip
  ! (z) and the count, of the owned particles of this rank (EL_BED_AUDIT_RECORD_SAMPLE)
  REAL*8, SAVE :: sample_bins(3, NBINS) = 0.0d0
  ! window audit: per height bin, sums of slip x, slip y, closure drag z, pressure force z, viscous
  ! force z, radius, particle density of the owned particles (EL_BED_AUDIT_RECORD_FORCES; the
  ! porosity, slip z and count are those of sample_bins)
  INTEGER, PARAMETER :: NFB = 7
  REAL*8, SAVE :: force_bins(NFB, NBINS) = 0.0d0
  ! table of the windows: kind 1 margin of win_a / win_b mean diameters from the floor / bed top,
  ! kind 2 fixed heights win_a / win_b above zmin [m]
  INTEGER, PARAMETER :: NWIN = 3
  CHARACTER(LEN=8), PARAMETER :: win_name(NWIN) = (/ 'm2d     ', 'm2.5d   ', 'interior' /)
  INTEGER, PARAMETER :: win_kind(NWIN) = (/ 1, 1, 2 /)
  REAL*8, PARAMETER :: win_a(NWIN) = (/ 2.0d0, 2.5d0, 0.005d0 /), win_b(NWIN) = (/ 2.0d0, 2.5d0, 0.020d0 /)
  ! the vertex-mode window (heights above zmin, m)
  REAL*8, PARAMETER :: BW_VM_LO = 0.005d0, BW_VM_HI = 0.020d0
  ! bed top and mean radius of the last audit line (window bounds of the window audit)
  INTEGER, SAVE :: last_step = -1
  REAL*8, SAVE :: last_z99 = 0.0d0, last_rmean = 0.0d0, last_area = 0.0d0

CONTAINS

  SUBROUTINE EL_WRITE_BED_AUDIT(pressure, time, dt, mfile, istep, mesh, ilev, density, gravity)
    REAL*8, INTENT(IN) :: pressure(:)
    REAL*8, INTENT(IN) :: time, dt
    INTEGER, INTENT(IN) :: mfile, istep
    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    REAL*8, INTENT(IN) :: density, gravity(3)

    REAL*8 :: zmin, zmax, tol, loc(7), glob(7), p_in, p_out, dp, dp_hydro
    REAL*8 :: zq(3), zlo, zhi, eps_bed, rmean, margin, zmax_loc, zmax_glob, floor_z, imp(3)
    REAL*8 :: win_loc(3), win(3), gradz_win, gradz_dyn, floor_bal, kern(3), kern_g(3)
    REAL*8 :: eps_kern, slip_kern
    INTEGER :: hist_loc(NBINS), hist(NBINS), nloc, nglob, ierr, i, ib
#ifdef HAVE_PE
    TYPE(tParticleData), ALLOCATABLE :: particles(:)
#endif

    IF (el_bed_audit_freq.LE.0) RETURN
    IF (MOD(ABS(istep), el_bed_audit_freq).NE.0) RETURN
    ! EL_CONFIG refuses the audit without a box domain at start-up on every rank, so this
    ! cannot be reached; a rank-local return would leave the other workers in the collectives
    IF (.NOT.el_domain_box_set) THEN
      WRITE(*,'(A,I0)') 'Fatal: EL_BED_AUDIT without a box domain on rank ', myid
      CALL MPI_Abort(MPI_COMM_SUBS, 914, ierr)
    END IF
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

    ! plane means (unavailable when the plane holds no element face of the mesh)
    p_in = EL_BED_UNAVAILABLE()
    p_out = EL_BED_UNAVAILABLE()
    IF (glob(1).GT.0.0d0) p_in = glob(2)/glob(1)
    IF (glob(3).GT.0.0d0) p_out = glob(4)/glob(3)
    dp = p_in - p_out
    dp_hydro = 0.0d0
    IF (el_fluid_gravity) dp_hydro = -density*gravity(3)*(zmax - zmin)
    floor_z = 0.0d0
    IF (dt.GT.0.0d0) floor_z = -glob(7)/dt

    ! bed heights and the porosity of the bed window (identical histogram on every rank; the
    ! branch depends on the global count only)
    zq = EL_BED_UNAVAILABLE()
    zmax_glob = MERGE(zmax_glob, EL_BED_UNAVAILABLE(), nglob.GT.0)
    eps_bed = EL_BED_UNAVAILABLE()
    zlo = EL_BED_UNAVAILABLE()
    zhi = EL_BED_UNAVAILABLE()
    margin = EL_BED_UNAVAILABLE()
    gradz_win = EL_BED_UNAVAILABLE()
    IF (nglob.GT.0) THEN
      zq(1) = EL_HISTOGRAM_QUANTILE(hist, zmin, zmax, 0.50d0)
      zq(2) = EL_HISTOGRAM_QUANTILE(hist, zmin, zmax, 0.95d0)
      zq(3) = EL_HISTOGRAM_QUANTILE(hist, zmin, zmax, 0.99d0)
      rmean = glob(6)/DBLE(nglob)
      margin = EL_BED_WINDOW_MARGIN(rmean)
      zlo = zmin + margin
      zhi = zq(3) + rmean - margin
      ! win: 1 solid volume in the window, 2 fluid-mesh volume in the window, 3 its integral
      ! of -dp/dz
      win_loc = 0.0d0
      IF (zhi.GT.zlo) THEN
#ifdef HAVE_PE
        DO i = 1, nloc
          win_loc(1) = win_loc(1) + EL_SPHERE_SLAB_VOLUME(particles(i)%position(3), &
                                                          particles(i)%radius, zlo, zhi)
        END DO
#endif
        CALL EL_WINDOW_GRADIENT_LOCAL(mesh, ilev, pressure, zlo, zhi, win_loc(2), win_loc(3))
      END IF
      CALL MPI_Allreduce(win_loc, win, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
      IF (zhi.GT.zlo .AND. glob(1).GT.0.0d0) eps_bed = 1.0d0 - win(1)/(glob(1)*(zhi - zlo))
      IF (win(2).GT.0.0d0) gradz_win = win(3)/win(2)
    END IF
#ifdef HAVE_PE
    IF (ALLOCATED(particles)) DEALLOCATE(particles)
#endif
    ! the window gradient without the hydrostatic part (the counterpart of dp_dyn)
    gradz_dyn = gradz_win - dp_hydro/(zmax - zmin)
    ! bed top of this step for the window audit (EL_WRITE_BED_WINDOWS)
    last_step = istep
    last_z99 = zq(3)
    last_rmean = EL_BED_UNAVAILABLE()
    IF (nglob.GT.0) last_rmean = glob(6)/DBLE(nglob)
    last_area = glob(1)

    ! wall force from the particle momentum balance of this step (forces actually armed only)
    floor_bal = EL_BED_FLOOR_BALANCE(el_pair_last_step, istep, el_pair_last_dp(3), &
                                     el_pair_last_expected(3), dt, el_apply_particle_forces)

    ! kernel-sampled porosity and slip of the particles with the centre in the window (no
    ! particle: zlo, zhi unavailable, nothing selected; every worker enters the reduction)
    kern = 0.0d0
    IF (nglob.GT.0) CALL EL_BED_WINDOW_SAMPLE_SUMS(zmin, zmax, zlo, zhi, kern)
    CALL MPI_Allreduce(kern, kern_g, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    eps_kern = EL_BED_UNAVAILABLE()
    slip_kern = EL_BED_UNAVAILABLE()
    IF (kern_g(3).GT.0.0d0) THEN
      eps_kern = kern_g(1)/kern_g(3)
      slip_kern = kern_g(2)/kern_g(3)
    END IF

    IF (myid.EQ.showid) THEN
      CALL WRITE_LINE(6)
      CALL WRITE_LINE(mfile)
    END IF

  CONTAINS

    SUBROUTINE WRITE_LINE(unit)
      INTEGER, INTENT(IN) :: unit
      WRITE(unit,'(A,ES14.6,A,I0,A,2ES16.8,A,2ES14.6,A,3ES16.8,A,I0,A,2ES14.6,A,4ES14.6,' // &
                 'A,ES14.6,A,2ES14.6,A,ES16.8,A,ES16.8,A,2ES16.8,A,ES14.6,A,ES16.8)') &
        'EL_BED_AUDIT t= ', time, ' step= ', istep, ' p_in_out= ', p_in, p_out, &
        ' area_in_out= ', glob(1), glob(3), ' dp_hydro_dyn= ', dp, dp_hydro, dp - dp_hydro, &
        ' n= ', nglob, ' weight_z_floor_z= ', glob(5), floor_z, &
        ' z50_95_99_max= ', zq, zmax_glob, ' eps_bed= ', eps_bed, ' window= ', zlo, zhi, &
        ' gradz_win= ', gradz_win, ' floor_bal_z= ', floor_bal, &
        ' eps_kern_slip_win= ', eps_kern, slip_kern, ' win_margin= ', margin, &
        ' gradz_win_dyn= ', gradz_dyn
    END SUBROUTINE WRITE_LINE

  END SUBROUTINE EL_WRITE_BED_AUDIT

  ! The value printed for an unavailable quantity of the audit line (quiet NaN).
  REAL*8 FUNCTION EL_BED_UNAVAILABLE()
    EL_BED_UNAVAILABLE = IEEE_VALUE(1.0d0, IEEE_QUIET_NAN)
  END FUNCTION EL_BED_UNAVAILABLE

  ! Margin of the bed window at both ends: two mean diameters (the spec) or the kernel support
  ! radius of the mean particle (the radius the transfer maps with), whichever is larger.
  REAL*8 FUNCTION EL_BED_WINDOW_MARGIN(rmean)
    REAL*8, INTENT(IN) :: rmean
    EL_BED_WINDOW_MARGIN = MAX(4.0d0*rmean, EL_KERNEL_SUPPORT_RADIUS(rmean))
  END FUNCTION EL_BED_WINDOW_MARGIN

  ! Wall force on the particles (z) from their momentum balance over the PE advance of step
  ! istep: (dp - expected)/dt, with dp the measured change of the particle momentum and expected
  ! the impulse of the forces the CFD side armed. Unavailable (NaN) unless the CFD side armed pe
  ! with the forces (forces_armed: ELApplyParticleForces = Yes; otherwise nothing of expected
  ! acted and the difference is no wall force), the balance belongs to this step and dt > 0.
  REAL*8 FUNCTION EL_BED_FLOOR_BALANCE(last_step, istep, dp_z, expected_z, dt, forces_armed)
    INTEGER, INTENT(IN) :: last_step, istep
    REAL*8, INTENT(IN) :: dp_z, expected_z, dt
    LOGICAL, INTENT(IN) :: forces_armed
    EL_BED_FLOOR_BALANCE = EL_BED_UNAVAILABLE()
    IF (forces_armed .AND. last_step.EQ.istep .AND. dt.GT.0.0d0) &
      EL_BED_FLOOR_BALANCE = (dp_z - expected_z)/dt
  END FUNCTION EL_BED_FLOOR_BALANCE

  ! Sums (porosity, slip z, count) of the recorded particle samples of this rank whose height bin
  ! has its centre in [zlo, zhi] (bins of EL_BED_AUDIT_RECORD_SAMPLE over [zmin, zmax]).
  SUBROUTINE EL_BED_WINDOW_SAMPLE_SUMS(zmin, zmax, zlo, zhi, sums)
    REAL*8, INTENT(IN) :: zmin, zmax, zlo, zhi
    REAL*8, INTENT(OUT) :: sums(3)
    REAL*8 :: zb
    INTEGER :: ib
    sums = 0.0d0
    DO ib = 1, NBINS
      zb = zmin + (DBLE(ib) - 0.5d0)*(zmax - zmin)/DBLE(NBINS)
      IF (zb.GE.zlo .AND. zb.LE.zhi) sums = sums + sample_bins(:, ib)
    END DO
  END SUBROUTINE EL_BED_WINDOW_SAMPLE_SUMS

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

  ! Kernel-sampled quantities of the particles for the window means of the audit: reset at the
  ! start of the post-advance transfer pass, one record per owned particle (EL_TRANSFER, only
  ! with SimPar@ELBedAuditFreq > 0). Rank-local, no collective.
  SUBROUTINE EL_BED_AUDIT_BEGIN_SAMPLES()
    sample_bins = 0.0d0
    IF (el_bed_window_audit) force_bins = 0.0d0
  END SUBROUTINE EL_BED_AUDIT_BEGIN_SAMPLES

  ! Window audit: the forces and slip of one owned particle, next to EL_BED_AUDIT_RECORD_SAMPLE
  ! (EL_TRANSFER, post-advance pass, only with SimPar@ELBedWindowAudit = Yes). Rank-local.
  SUBROUTINE EL_BED_AUDIT_RECORD_FORCES(z, slip, drag, pressure, viscous, radius, rho_p)
    REAL*8, INTENT(IN) :: z, slip(3), drag(3), pressure(3), viscous(3), radius, rho_p
    INTEGER :: ib
    IF (.NOT.el_domain_box_set) RETURN
    ib = 1 + INT(DBLE(NBINS)*(z - el_domain_box(5))/(el_domain_box(6) - el_domain_box(5)))
    ib = MIN(NBINS, MAX(1, ib))
    force_bins(1, ib) = force_bins(1, ib) + slip(1)
    force_bins(2, ib) = force_bins(2, ib) + slip(2)
    force_bins(3, ib) = force_bins(3, ib) + drag(3)
    force_bins(4, ib) = force_bins(4, ib) + pressure(3)
    force_bins(5, ib) = force_bins(5, ib) + viscous(3)
    force_bins(6, ib) = force_bins(6, ib) + radius
    force_bins(7, ib) = force_bins(7, ib) + rho_p
  END SUBROUTINE EL_BED_AUDIT_RECORD_FORCES

  ! Bounds [zl, zh] of window k (table above) for floor zmin, bed top z99 + r and mean radius r;
  ! unavailable (NaN) when the bed top is.
  PURE SUBROUTINE EL_BED_WINDOW_BOUNDS(k, zmin, z99, rmean, zl, zh)
    INTEGER, INTENT(IN) :: k
    REAL*8, INTENT(IN) :: zmin, z99, rmean
    REAL*8, INTENT(OUT) :: zl, zh
    IF (win_kind(k).EQ.1) THEN
      zl = zmin + win_a(k)*2.0d0*rmean
      zh = z99 + rmean - win_b(k)*2.0d0*rmean
    ELSE
      zl = zmin + win_a(k)
      zh = zmin + win_b(k)
    END IF
  END SUBROUTINE EL_BED_WINDOW_BOUNDS

  ! The window lines of this step (see the header). Collective over the workers; called by every
  ! worker on the audit steps when SimPar@ELBedWindowAudit = Yes, after EL_WRITE_BED_AUDIT of the
  ! step. vm_sum(3,4) / vm_min(4) / vm_max(4): the vertex-mode sums of this rank (sum 1/m, sum w/m,
  ! sum w^2/m), minima and maxima per dof class (vertex, edge, face, centre), formed by the caller.
  SUBROUTINE EL_WRITE_BED_WINDOWS(pressure, time, mfile, istep, mesh, ilev, density, viscosity, &
                                  vm_sum, vm_min, vm_max)
    REAL*8, INTENT(IN) :: pressure(:), time, density, viscosity
    INTEGER, INTENT(IN) :: mfile, istep, ilev
    TYPE(tMultiMesh), INTENT(IN) :: mesh
    REAL*8, INTENT(IN) :: vm_sum(3,4), vm_min(4), vm_max(4)

    INTEGER, PARAMETER :: NPW = 3 + 10
    REAL*8 :: loc(NPW*NWIN + 12), glob(NPW*NWIN + 12), vmin(4), vmax(4), zmin, zl(NWIN), zh(NWIN)
    REAL*8 :: zb, state(8), fm(3), bm, cnt, mean_eps, slip(3), rmean_w, rho_w, eps_geo, gz, nan
    REAL*8 :: vm_mean(4), vm_std(4)
    INTEGER :: k, ib, i, nloc, ierr, j0
    CHARACTER(LEN=8), PARAMETER :: cname(4) = (/ 'vertex  ', 'edge    ', 'face    ', 'centre  ' /)
#ifdef HAVE_PE
    TYPE(tParticleData), ALLOCATABLE :: particles(:)
#endif

    IF (.NOT.el_bed_window_audit) RETURN
    IF (el_bed_audit_freq.LE.0) RETURN
    IF (MOD(ABS(istep), el_bed_audit_freq).NE.0) RETURN
    nan = EL_BED_UNAVAILABLE()
    zmin = el_domain_box(5)
    DO k = 1, NWIN
      CALL EL_BED_WINDOW_BOUNDS(k, zmin, last_z99, last_rmean, zl(k), zh(k))
    END DO
    IF (last_step.NE.istep) THEN
      zl = nan
      zh = nan
    END IF

    ! per window: 1 mesh volume, 2 int -dp/dz, 3 solid volume, 4..13 particle sums (count, eps,
    ! slip x y z, drag z, pressure z, viscous z, radius, rho_p); then the vertex-mode sums
    loc = 0.0d0
#ifdef HAVE_PE
    nloc = numLocalParticles()
    ALLOCATE(particles(MAX(nloc, 0)))
    IF (nloc.GT.0) CALL getAllParticles(particles)
#endif
    DO k = 1, NWIN
      j0 = NPW*(k - 1)
      IF (.NOT.(zh(k).GT.zl(k))) CYCLE
      CALL EL_WINDOW_GRADIENT_LOCAL(mesh, ilev, pressure, zl(k), zh(k), loc(j0+1), loc(j0+2))
#ifdef HAVE_PE
      DO i = 1, nloc
        loc(j0+3) = loc(j0+3) + EL_SPHERE_SLAB_VOLUME(particles(i)%position(3), &
                                                      particles(i)%radius, zl(k), zh(k))
      END DO
#endif
      DO ib = 1, NBINS
        zb = zmin + (DBLE(ib) - 0.5d0)*(el_domain_box(6) - zmin)/DBLE(NBINS)
        IF (zb.LT.zl(k) .OR. zb.GT.zh(k)) CYCLE
        loc(j0+4) = loc(j0+4) + sample_bins(3, ib)
        loc(j0+5) = loc(j0+5) + sample_bins(1, ib)
        loc(j0+6) = loc(j0+6) + force_bins(1, ib)
        loc(j0+7) = loc(j0+7) + force_bins(2, ib)
        loc(j0+8) = loc(j0+8) + sample_bins(2, ib)
        loc(j0+9:j0+13) = loc(j0+9:j0+13) + force_bins(3:7, ib)
      END DO
    END DO
#ifdef HAVE_PE
    DEALLOCATE(particles)
#endif
    loc(NPW*NWIN+1:NPW*NWIN+12) = RESHAPE(vm_sum, (/12/))
    CALL MPI_Allreduce(loc, glob, NPW*NWIN + 12, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(vm_min, vmin, 4, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(vm_max, vmax, 4, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_SUBS, ierr)
    IF (myid.NE.showid) RETURN

    DO k = 1, 4
      vm_mean(k) = nan
      vm_std(k) = nan
      cnt = glob(NPW*NWIN + 3*(k-1) + 1)
      IF (cnt.GT.0.0d0) THEN
        vm_mean(k) = glob(NPW*NWIN + 3*(k-1) + 2)/cnt
        vm_std(k) = SQRT(MAX(0.0d0, glob(NPW*NWIN + 3*(k-1) + 3)/cnt - vm_mean(k)**2))
      ELSE
        vmin(k) = nan
        vmax(k) = nan
      END IF
    END DO

    DO k = 1, NWIN
      j0 = NPW*(k - 1)
      gz = nan
      IF (glob(j0+1).GT.0.0d0) gz = glob(j0+2)/glob(j0+1)
      eps_geo = nan
      IF (zh(k).GT.zl(k) .AND. last_area.GT.0.0d0) eps_geo = 1.0d0 - glob(j0+3)/(last_area*(zh(k) - zl(k)))
      cnt = glob(j0+4)
      mean_eps = nan
      slip = nan
      fm = nan
      bm = nan
      IF (cnt.GT.0.0d0) THEN
        mean_eps = glob(j0+5)/cnt
        slip = glob(j0+6:j0+8)/cnt
        rmean_w = glob(j0+12)/cnt
        rho_w = glob(j0+13)/cnt
        ! one particle of mean radius and density at rest in the mean slip: B(<eps>, |<s>|) <s>
        state = 0.0d0
        state(7) = rmean_w
        state(8) = rho_w
        CALL EL_DRAG_CLOSURE(state, slip, density, viscosity, mean_eps, fm, bm)
        fm = EL_MM_DRAG_SCALE(mean_eps)*fm
      END IF
      CALL WRITE_WINDOW(6)
      CALL WRITE_WINDOW(mfile)
    END DO
    CALL WRITE_VMODE(6)
    CALL WRITE_VMODE(mfile)

  CONTAINS

    SUBROUTINE WRITE_WINDOW(unit)
      INTEGER, INTENT(IN) :: unit
      REAL*8 :: per(4)
      per = nan
      IF (cnt.GT.0.0d0) per = glob(j0+9:j0+12)/cnt
      WRITE(unit,'(A,ES16.8,A,I0,A,A,A,2ES16.8,A,ES16.8,A,I0,2(A,ES16.8),A,3ES16.8,4(A,ES16.8))') &
        'EL_BED_WINDOW t= ', time, ' step= ', istep, ' win= ', TRIM(win_name(k)), &
        ' z_lo_hi= ', zl(k), zh(k), ' gradz= ', gz, ' n= ', NINT(cnt), ' eps_geo= ', eps_geo, &
        ' eps_kern= ', mean_eps, ' slip= ', slip, ' drag_mean= ', per(1), ' drag_of_means= ', fm(3), &
        ' pres_mean= ', per(2), ' visc_mean= ', per(3)
    END SUBROUTINE WRITE_WINDOW

    SUBROUTINE WRITE_VMODE(unit)
      INTEGER, INTENT(IN) :: unit
      INTEGER :: c
      WRITE(unit,'(A,ES16.8,A,I0,A,2ES16.8,4(A,A,5ES16.8))') 'EL_BED_VERTEX_MODE t= ', time, &
        ' step= ', istep, ' z_lo_hi= ', zmin + BW_VM_LO, zmin + BW_VM_HI, &
        (' '//TRIM(cname(c))//'= ', '', glob(NPW*NWIN + 3*(c-1) + 1), vm_mean(c), vm_std(c), &
         vmin(c), vmax(c), c = 1, 4)
    END SUBROUTINE WRITE_VMODE

  END SUBROUTINE EL_WRITE_BED_WINDOWS

  SUBROUTINE EL_BED_AUDIT_RECORD_SAMPLE(z, eps, slip_z)
    REAL*8, INTENT(IN) :: z, eps, slip_z
    INTEGER :: ib
    IF (.NOT.el_domain_box_set) RETURN
    ib = 1 + INT(DBLE(NBINS)*(z - el_domain_box(5))/(el_domain_box(6) - el_domain_box(5)))
    ib = MIN(NBINS, MAX(1, ib))
    sample_bins(1, ib) = sample_bins(1, ib) + eps
    sample_bins(2, ib) = sample_bins(2, ib) + slip_z
    sample_bins(3, ib) = sample_bins(3, ib) + 1.0d0
  END SUBROUTINE EL_BED_AUDIT_RECORD_SAMPLE

  ! Volume of the elements of this rank inside the slab za < z < zb and its integral of -dp/dz
  ! (P1 pressure: the gradient is constant per element, entry 4 of its block). The part of an
  ! element in the slab is taken as the overlap of its vertex z range with the slab times its
  ! volume: exact for an element of constant horizontal cross-section between two z planes (a
  ! prism along z, the case-11 mesh); an estimate otherwise, also for an element bounded by two
  ! z planes whose cross-section varies with z (a tapered hex).
  SUBROUTINE EL_WINDOW_GRADIENT_LOCAL(mesh, ilev, pressure, za, zb, volume, gradient_integral)
    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    REAL*8, INTENT(IN) :: pressure(:), za, zb
    REAL*8, INTENT(INOUT) :: volume, gradient_integral
    REAL*8 :: ze0, ze1, overlap, v
    INTEGER :: iel, i

    DO iel = 1, mesh%level(ilev)%nel
      ze0 = HUGE(1.0d0)
      ze1 = -HUGE(1.0d0)
      DO i = 1, 8
        ze0 = MIN(ze0, mesh%level(ilev)%dcorvg(3, mesh%level(ilev)%kvert(i,iel)))
        ze1 = MAX(ze1, mesh%level(ilev)%dcorvg(3, mesh%level(ilev)%kvert(i,iel)))
      END DO
      overlap = MIN(ze1, zb) - MAX(ze0, za)
      IF (overlap.LE.0.0d0 .OR. ze1.LE.ze0) CYCLE
      v = mesh%level(ilev)%dvol(iel)*overlap/(ze1 - ze0)
      volume = volume + v
      gradient_integral = gradient_integral - v*pressure(4*(iel-1)+4)
    END DO
  END SUBROUTINE EL_WINDOW_GRADIENT_LOCAL

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
