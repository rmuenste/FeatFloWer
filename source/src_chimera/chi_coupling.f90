!=========================================================================
! CHI_COUPLING - Chimera-S coupling state and hook implementations,
! Layer H (design: chimera-integration-design.md v3, sections 3, 4, 6, 7;
! milestone M1: static steady Chimera-S).
!
! This module owns every call from the Chimera component into legacy
! global state (var_QuadScalar, def_FEAT COMMON, PP3D_MPI, E013 sync).
! Existing solver code never USEs it directly - only through CHIMERA_API.
!
! Data flow per coupling update (CHI_COUPLING_BEGIN_STEP, workers only):
!   cached outer quadrature points -> CHI_EXCHANGE_BG_EVAL (collective on
!   MPI_COMM_SUBS) -> Robin data h = sigma(u,p).n - alpha (u.n) u (design
!   section 5) -> owner worker solves the submesh (Picard, direct) ->
!   solution broadcast -> fringe values evaluated from the replicated
!   submesh solution (bit-identical on all workers, so shared interface
!   dofs need no arbitration) -> forces (design section 5 conventions),
!   reported once as "ChimeraForce:" protocol lines.
!
! Markers (two-array scheme, v2-review blocker 2): marker_kind
! (0 free / 1 fringe / 2 hole) and marker_pid, classified geometrically
! on the finest background level, synchronised across partitions with
! the production E013Max_SUPER (numeric MAX = hole > fringe > free) and a
! second conditional pass for the body id.
!
! Hole-region background pressure dofs are left free (FBM precedent).
! The matrix row filter follows the FictKNPR practice of indexing the
! finest-level marker array on every multigrid level (coarse levels
! affect only the preconditioner, never the converged solution).
!=========================================================================
MODULE CHI_COUPLING

  USE PP3D_MPI, ONLY: myid, master, showid, subnodes, MPI_COMM_SUBS
  USE var_QuadScalar, ONLY: mg_mesh, myQ2Coor, Properties, postParams, mg_qMat
  USE def_FEAT, ONLY: NLMIN, NLMAX, ILEV, TIMENS, TSTEP
  USE CHIMERA_CONFIG, ONLY: chimera_outer_bc, chimera_particle_file, &
    chimera_submesh_file, chimera_submesh_nlmax, chimera_robin_alpha, &
    chimera_sub_nl, chimera_write_vtk, bChimeraW, chimera_gamma_max, &
    chimera_penalty_lumped, chimera_proj_cap, chimera_beta_full, chimera_beta_zero, &
    chimera_coupling_relax
  USE CHI_SUBMESH, ONLY: tChimeraSubmesh, CHI_SHAPE_CYLINDER_Z, &
    CHI_SHAPE_SPHERE, CHI_SURF_INNER, CHI_SURF_OUTER, CHI_SURF_ZLO, &
    CHI_SURF_ZHI, CHI_SUBMESH_RELEASE
  USE CHI_LEGACY_MESH_ADAPTER, ONLY: CHI_LOAD_SUBMESH, CHI_RELEASE_SUBMESH_MESH
  USE CHI_SOLVER, ONLY: tChiSubSolver, CHI_SOLVER_INIT, CHI_SOLVE_STEADY_TAB, &
    CHI_SOLVER_ADVANCE, CHI_SOLVER_RELEASE, CHI_DIR_W, CHI_DIR_ALL
  USE CHI_FORCES, ONLY: CHI_COMPUTE_FORCES
  USE CHI_LOCATOR, ONLY: tChimeraLocator, CHI_LOCATOR_BUILD, CHI_LOCATE, &
    CHI_LOCATE_NEAREST, CHI_LOCATOR_RELEASE
  USE CHI_FEM_EVAL, ONLY: CHI_EVAL_FIELD_AT
  USE CHI_KERNELS, ONLY: CHI_ROBIN_POINTS
  USE CHI_EXCHANGE, ONLY: CHI_EXCHANGE_BG_EVAL, CHI_EXCHANGE_BCAST, &
    CHI_EXCHANGE_RANK, CHI_EXCHANGE_MAX_INT, CHI_BG_NVAL, CHI_EXCHANGE_ALLSUM
  USE CHI_MARKERS, ONLY: tChiBody, CHI_BODY_CYLINDER_Z, CHI_BODY_SPHERE, &
    CHI_MARK_FREE, CHI_MARK_FRINGE, CHI_MARK_HOLE, CHI_CLASSIFY_MARKERS
  USE CHI_OUTPUT, ONLY: CHI_WRITE_SUBMESH_VTK
  USE CHI_PENALTY, ONLY: tChiPenaltyTab, CHI_PENALTY_NQ, CHI_PENALTY_TABULATE, &
    CHI_PENALTY_RELEASE, CHI_PENALTY_ASSEMBLE_D, CHI_PENALTY_ASSEMBLE_G, &
    CHI_PENALTY_DIAG, CHI_PENALTY_MATVEC, CHI_PCG3, chi_filter3_iface, &
    CHI_PENALTY_LUMP, CHI_PENALTY_NODAL

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: CHI_COUPLING_INIT
  PUBLIC :: CHI_COUPLING_BEGIN_STEP
  PUBLIC :: CHI_COUPLING_APPLY_DEF
  PUBLIC :: CHI_COUPLING_APPLY_VAL
  PUBLIC :: CHI_COUPLING_FILTER_MAT
  PUBLIC :: CHI_COUPLING_FILTER_MAT_9
  PUBLIC :: CHI_COUPLING_FINALIZE
  PUBLIC :: CHI_COUPLING_ACTIVE
  PUBLIC :: CHI_COUPLING_ADD_MAT
  PUBLIC :: CHI_COUPLING_ADD_DEFECT
  PUBLIC :: CHI_COUPLING_ADD_RHS
  PUBLIC :: CHI_COUPLING_CORRECT
  PUBLIC :: CHI_COUPLING_WRITE_RESTART
  PUBLIC :: CHI_COUPLING_READ_RESTART
  PUBLIC :: CHI_COUPLING_ADD_PRESSURE_MASS

  EXTERNAL E013Max_SUPER, E013Sum, E013Sum3

  ! ---- private SAVEd state (workers only; the master holds nothing) ----
  LOGICAL, SAVE :: active = .FALSE.
  INTEGER, SAVE :: mfile_unit = 0
  INTEGER, SAVE :: nstep_done = 0
  INTEGER, SAVE :: myrank = -1            ! rank in MPI_COMM_SUBS
  LOGICAL, SAVE :: dirichlet_mode = .FALSE.
  REAL*8,  SAVE :: rho = 1d0, mu = 1d0

  INTEGER, SAVE :: nsub = 0
  TYPE(tChimeraSubmesh), ALLOCATABLE, SAVE :: atm(:)
  TYPE(tChiSubSolver),   ALLOCATABLE, SAVE :: slv(:)
  TYPE(tChimeraLocator), ALLOCATABLE, SAVE :: subloc(:)
  TYPE(tChiBody),        ALLOCATABLE, SAVE :: bodies(:)
  INTEGER, ALLOCATABLE, SAVE :: owner(:)      ! owning rank in MPI_COMM_SUBS

  ! background finest level
  TYPE(tChimeraLocator), SAVE :: bgloc
  INTEGER, SAVE :: nbgdof = 0
  INTEGER, ALLOCATABLE, SAVE :: marker_kind(:), marker_pid(:)
  REAL*8,  ALLOCATABLE, SAVE :: fringeU(:), fringeV(:), fringeW(:)

  ! constrained background dofs with donor cache (iel = 0: hole node)
  INTEGER, SAVE :: ncon = 0
  INTEGER, ALLOCATABLE, SAVE :: con_dof(:), con_sub(:), con_iel(:)
  REAL*8,  ALLOCATABLE, SAVE :: con_xi(:,:)

  ! Robin quadrature points, concatenated over submeshes
  INTEGER, SAVE :: nqtot = 0
  INTEGER, ALLOCATABLE, SAVE :: qoff(:)        ! (nsub+1), face offsets
  REAL*8,  ALLOCATABLE, SAVE :: qpts(:,:), qnrm(:,:)   ! (3, 9*nfaces)

  ! outer Q2 nodes (Dirichlet diagnostic mode), concatenated
  INTEGER, SAVE :: ndtot = 0
  INTEGER, ALLOCATABLE, SAVE :: doff(:), ddof(:)
  REAL*8,  ALLOCATABLE, SAVE :: dpts(:,:)

  ! ---- weak variant (Chimera-W, Phase 4): penalty operator per level ----
  ! D shares the Q2 pattern mg_qMat(lev) (design section 3, name avoids
  ! the viscous DMat); g and the Jacobi/weight vectors live on NLMAX.
  TYPE tPenLevel
    REAL*8, ALLOCATABLE :: d(:)       ! consistent D on the level pattern (rank-partial)
    REAL*8, ALLOCATABLE :: dl(:)      ! row-sum lumped D_L, rank-partial
    REAL*8, ALLOCATABLE :: dlg(:)     ! D_L assembly-summed (global)
  END TYPE tPenLevel
  LOGICAL, SAVE :: lumped = .TRUE.
  REAL*8,  SAVE :: pen_dt_cmat = -1d0   ! dt frozen into the projection operator
  ! nodal (lumped) variant: finest-level penalised nodes and their donors
  INTEGER, SAVE :: nnod = 0
  INTEGER, ALLOCATABLE, SAVE :: nd_dof(:), nd_sub(:), nd_iel(:)
  REAL*8,  ALLOCATABLE, SAVE :: nd_xi(:,:)
  LOGICAL, SAVE :: weak_mode = .FALSE.
  INTEGER, SAVE :: pen_lmin = 0, pen_lmax = -1
  TYPE(tChiPenaltyTab), ALLOCATABLE, SAVE :: ptab(:)
  TYPE(tPenLevel),      ALLOCATABLE, SAVE :: pmat(:)
  INTEGER, ALLOCATABLE, SAVE :: pd_iel(:,:)          ! donors of the finest points
  REAL*8,  ALLOCATABLE, SAVE :: pd_xi(:,:,:)
  REAL*8,  ALLOCATABLE, SAVE :: gU(:), gV(:), gW(:)  ! rank-partial g
  REAL*8,  ALLOCATABLE, SAVE :: pdiag(:), pwts(:)    ! global diag(D), 1/share count
  LOGICAL, SAVE :: g_valid = .FALSE.
  REAL*8,  ALLOCATABLE, SAVE :: gU_prev(:), gV_prev(:), gW_prev(:)
  LOGICAL, SAVE :: fringe_valid = .FALSE.
  INTEGER, SAVE :: ncorr_it = 0
  REAL*8,  SAVE :: corr_resid = 0d0

CONTAINS

  LOGICAL FUNCTION CHI_COUPLING_ACTIVE()
    CHI_COUPLING_ACTIVE = active
  END FUNCTION CHI_COUPLING_ACTIVE

  !=======================================================================
  ! Initialization (hook H5, app-local).  Workers build the replicated
  ! submesh state and the background markers; the master returns.
  !=======================================================================
  SUBROUTINE CHI_COUPLING_INIT(mfile)
    INTEGER, INTENT(IN) :: mfile

    INTEGER :: k, i, lev, nvt, net, nat, nel, nholes, nfringe, nmax
    INTEGER :: nf, ifc, q, e, fl, nd, ilev_save
    LOGICAL :: ok, found
    REAL*8 :: excess, maxexcess
    INTEGER, ALLOCATABLE :: kind0(:), tmp(:)

    mfile_unit = mfile
    nstep_done = 0
    IF (myid .EQ. master) RETURN

    myrank = CHI_EXCHANGE_RANK(MPI_COMM_SUBS)
    dirichlet_mode = (TRIM(chimera_outer_bc) .EQ. 'dirichlet')
    weak_mode = bChimeraW
    rho = Properties%Density(1)
    mu  = Properties%Density(1)*Properties%Viscosity(1)   ! Prop@Viscosity is kinematic

    !---- bodies ---------------------------------------------------------
    CALL read_particle_file(chimera_particle_file)
    IF (nsub .LE. 0) THEN
      WRITE(*,'(A)') 'CHI_COUPLING error: no bodies in ' // TRIM(chimera_particle_file)
      STOP 1
    END IF

    !---- submeshes (replicated on every worker) --------------------------
    ALLOCATE(slv(nsub), subloc(nsub), owner(nsub))
    DO k = 1, nsub
      atm(k)%nlmax = chimera_submesh_nlmax
      CALL CHI_LOAD_SUBMESH(atm(k), TRIM(chimera_submesh_file), ok)
      IF (.NOT. ok) THEN
        WRITE(*,'(A,I0)') 'CHI_COUPLING error: cannot load submesh for body ', k
        STOP 1
      END IF
      lev = atm(k)%nlmax
      CALL CHI_LOCATOR_BUILD(subloc(k), atm(k)%mesh%level(lev)%dcorvg, &
        atm(k)%mesh%level(lev)%kvert, atm(k)%mesh%level(lev)%nel, &
        atm(k)%mesh%level(lev)%nvt)
      CALL CHI_SOLVER_INIT(slv(k), atm(k), ok)
      IF (.NOT. ok) THEN
        WRITE(*,'(A,I0)') 'CHI_COUPLING error: solver init failed for body ', k
        STOP 1
      END IF
      owner(k) = MOD(k-1, subnodes)
    END DO

    !---- background finest level -----------------------------------------
    lev = NLMAX
    nvt = mg_mesh%level(lev)%nvt
    net = mg_mesh%level(lev)%net
    nat = mg_mesh%level(lev)%nat
    nel = mg_mesh%level(lev)%nel
    nbgdof = nvt + net + nat + nel
    CALL CHI_LOCATOR_BUILD(bgloc, mg_mesh%level(lev)%dcorvg, &
      mg_mesh%level(lev)%kvert, nel, nvt)

    ALLOCATE(marker_kind(nbgdof), marker_pid(nbgdof))
    ALLOCATE(fringeU(nbgdof), fringeV(nbgdof), fringeW(nbgdof))
    fringeU = 0d0; fringeV = 0d0; fringeW = 0d0

    IF (weak_mode) THEN
      ! Chimera-W: no nodal constraints; the coupling is the penalty.
      marker_kind = CHI_MARK_FREE
      marker_pid = 0
    ELSE
    CALL CHI_CLASSIFY_MARKERS(nel, nvt, net, nat, mg_mesh%level(lev)%kvert, &
      mg_mesh%level(lev)%kedge, mg_mesh%level(lev)%karea, myQ2Coor, &
      nsub, bodies, marker_kind, marker_pid)

    ! Parallel synchronisation: MAX on the kind (hole > fringe > free),
    ! then the body id among partitions that carry the winning kind.
    ALLOCATE(kind0(nbgdof), tmp(nbgdof))
    kind0 = marker_kind
    ilev_save = ILEV
    ILEV = NLMAX
    CALL E013Max_SUPER(marker_kind)
    DO i = 1, nbgdof
      IF (marker_kind(i) .GT. CHI_MARK_FREE .AND. kind0(i) .EQ. marker_kind(i)) THEN
        tmp(i) = marker_pid(i)
      ELSE
        tmp(i) = 0
      END IF
    END DO
    CALL E013Max_SUPER(tmp)
    marker_pid = tmp
    ILEV = ilev_save
    DEALLOCATE(kind0, tmp)
    END IF

    !---- constraint list + donor cache -----------------------------------
    ncon = 0
    DO i = 1, nbgdof
      IF (marker_kind(i) .NE. CHI_MARK_FREE) ncon = ncon + 1
    END DO
    ALLOCATE(con_dof(MAX(ncon,1)), con_sub(MAX(ncon,1)), con_iel(MAX(ncon,1)), &
             con_xi(3,MAX(ncon,1)))
    ncon = 0
    nholes = 0
    nfringe = 0
    maxexcess = 0d0
    DO i = 1, nbgdof
      IF (marker_kind(i) .EQ. CHI_MARK_FREE) CYCLE
      ncon = ncon + 1
      con_dof(ncon) = i
      con_sub(ncon) = marker_pid(i)
      con_iel(ncon) = 0
      con_xi(:,ncon) = 0d0
      IF (marker_kind(i) .EQ. CHI_MARK_HOLE) THEN
        nholes = nholes + 1
        CYCLE
      END IF
      nfringe = nfringe + 1
      k = marker_pid(i)
      IF (k .LT. 1 .OR. k .GT. nsub) THEN
        WRITE(*,'(A,I0)') 'CHI_COUPLING error: fringe node without body id, dof ', i
        STOP 1
      END IF
      lev = atm(k)%nlmax
      CALL CHI_LOCATE(subloc(k), atm(k)%mesh%level(lev)%dcorvg, &
        atm(k)%mesh%level(lev)%kvert, myQ2Coor(:,i), con_iel(ncon), &
        con_xi(:,ncon), found)
      IF (.NOT. found) THEN
        CALL CHI_LOCATE_NEAREST(subloc(k), atm(k)%mesh%level(lev)%dcorvg, &
          atm(k)%mesh%level(lev)%kvert, myQ2Coor(:,i), con_iel(ncon), &
          con_xi(:,ncon), excess, found)
        IF (.NOT. found) THEN
          WRITE(*,'(A,I0,A,3ES14.6)') 'CHI_COUPLING error: fringe dof ', i, &
            ' not inside its atmosphere, x =', myQ2Coor(:,i)
          STOP 1
        END IF
        maxexcess = MAX(maxexcess, excess)
      END IF
    END DO

    !---- Robin quadrature point cache ------------------------------------
    ALLOCATE(qoff(nsub+1))
    qoff(1) = 0
    DO k = 1, nsub
      qoff(k+1) = qoff(k) + SIZE(atm(k)%outerFaces,2)
    END DO
    nqtot = 9*qoff(nsub+1)
    ALLOCATE(qpts(3,MAX(nqtot,1)), qnrm(3,MAX(nqtot,1)))
    DO k = 1, nsub
      lev = atm(k)%nlmax
      nf = SIZE(atm(k)%outerFaces,2)
      IF (nf .GT. 0) CALL CHI_ROBIN_POINTS(atm(k)%outerFaces, nf, &
        atm(k)%mesh%level(lev)%kvert, atm(k)%mesh%level(lev)%dcorvg, &
        qpts(:,9*qoff(k)+1:9*qoff(k+1)), qnrm(:,9*qoff(k)+1:9*qoff(k+1)))
    END DO

    !---- outer Q2 nodes (Dirichlet diagnostic mode) ----------------------
    ALLOCATE(doff(nsub+1))
    doff(1) = 0
    DO k = 1, nsub
      nd = 0
      DO i = 1, atm(k)%ndof
        IF (IAND(atm(k)%dofmask(i), CHI_SURF_OUTER) .NE. 0) nd = nd + 1
      END DO
      doff(k+1) = doff(k) + nd
    END DO
    ndtot = doff(nsub+1)
    ALLOCATE(ddof(MAX(ndtot,1)), dpts(3,MAX(ndtot,1)))
    nd = 0
    DO k = 1, nsub
      DO i = 1, atm(k)%ndof
        IF (IAND(atm(k)%dofmask(i), CHI_SURF_OUTER) .EQ. 0) CYCLE
        nd = nd + 1
        ddof(nd) = i
        dpts(:,nd) = atm(k)%q2coor(:,i)
      END DO
    END DO

    IF (weak_mode) CALL setup_penalty()

    active = .TRUE.

    nmax = CHI_EXCHANGE_MAX_INT(MPI_COMM_SUBS, ncon)
    IF (myid .EQ. showid) THEN
      IF (weak_mode) THEN
        WRITE(*,'(A,I0,A,I0,A,ES10.2)') 'Chimera: ', nsub, ' body/bodies, ', &
          subnodes, ' worker(s); weak (Chimera-W) coupling, gamma_max = ', &
          chimera_gamma_max
      ELSE
        WRITE(*,'(A,I0,A,I0,A)') 'Chimera: ', nsub, ' body/bodies, ', &
          subnodes, ' worker(s); strong (Chimera-S) coupling'
      END IF
      WRITE(*,'(A,I0,A,I0,A,I0)') 'Chimera: submesh Q2 dofs = ', atm(1)%ndof, &
        ', unknowns = ', slv(1)%n, ', levels = ', atm(1)%nlmax
      WRITE(*,'(A,I0,A,I0,A,I0,A,ES10.2)') 'Chimera: rank ', myid, &
        ' hole dofs = ', nholes, ', fringe dofs = ', nfringe, &
        ', max donor extrapolation = ', maxexcess
      WRITE(mfile,'(A,I0,A,I0,A,I0)') 'Chimera: bodies = ', nsub, &
        ', submesh levels = ', atm(1)%nlmax, ', submesh Q2 dofs = ', atm(1)%ndof
      WRITE(mfile,'(A,A)') 'Chimera: outer BC = ', TRIM(chimera_outer_bc)
      WRITE(mfile,'(A,I0)') 'Chimera: max constrained dofs per worker = ', nmax
    END IF
  END SUBROUTINE CHI_COUPLING_INIT

  !-----------------------------------------------------------------------
  ! Particle file: '#' comment lines; first data line = number of
  ! bodies; then one line per body:
  !   cylinder_z  cx cy cz  radius  H  zlo zhi
  !   sphere      cx cy cz  radius  H
  ! (H = atmosphere width, outer radius = radius + H; cz unused for
  ! cylinders, which are treated as infinite along z).
  !-----------------------------------------------------------------------
  SUBROUTINE read_particle_file(fname)
    CHARACTER(*), INTENT(IN) :: fname
    INTEGER :: iu, ios, k
    CHARACTER(LEN=512) :: line
    CHARACTER(LEN=32) :: shape
    REAL*8 :: c(3), r, h, zlo, zhi

    iu = 772
    OPEN(UNIT=iu, FILE=TRIM(fname), STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios .NE. 0) THEN
      WRITE(*,'(A)') 'CHI_COUPLING error: cannot open ChimeraParticleFile ' // TRIM(fname)
      STOP 1
    END IF
    CALL next_data_line(iu, line, ios)
    IF (ios .NE. 0) THEN
      WRITE(*,'(A)') 'CHI_COUPLING error: empty ChimeraParticleFile'
      STOP 1
    END IF
    READ(line, *, IOSTAT=ios) nsub
    IF (ios .NE. 0 .OR. nsub .LT. 1) THEN
      WRITE(*,'(A)') 'CHI_COUPLING error: bad body count in ChimeraParticleFile'
      STOP 1
    END IF
    ALLOCATE(bodies(nsub), atm(nsub))
    DO k = 1, nsub
      CALL next_data_line(iu, line, ios)
      IF (ios .NE. 0) THEN
        WRITE(*,'(A,I0)') 'CHI_COUPLING error: missing body line ', k
        STOP 1
      END IF
      READ(line, *, IOSTAT=ios) shape
      CALL lowercase(shape)
      SELECT CASE (TRIM(shape))
      CASE ('cylinder_z')
        READ(line, *, IOSTAT=ios) shape, c, r, h, zlo, zhi
        IF (ios .NE. 0) THEN
          WRITE(*,'(A,I0)') 'CHI_COUPLING error: cylinder_z line needs cx cy cz r H zlo zhi, body ', k
          STOP 1
        END IF
        bodies(k)%shape = CHI_BODY_CYLINDER_Z
        atm(k)%shape = CHI_SHAPE_CYLINDER_Z
        atm(k)%zlo = zlo
        atm(k)%zhi = zhi
      CASE ('sphere')
        READ(line, *, IOSTAT=ios) shape, c, r, h
        IF (ios .NE. 0) THEN
          WRITE(*,'(A,I0)') 'CHI_COUPLING error: sphere line needs cx cy cz r H, body ', k
          STOP 1
        END IF
        bodies(k)%shape = CHI_BODY_SPHERE
        atm(k)%shape = CHI_SHAPE_SPHERE
      CASE DEFAULT
        WRITE(*,'(A,A)') 'CHI_COUPLING error: unknown body shape ', TRIM(shape)
        STOP 1
      END SELECT
      IF (r .LE. 0d0 .OR. h .LE. 0d0) THEN
        WRITE(*,'(A,I0)') 'CHI_COUPLING error: radius and H must be > 0, body ', k
        STOP 1
      END IF
      bodies(k)%center = c
      bodies(k)%radius = r
      atm(k)%center = c
      atm(k)%radius_inner = r
      atm(k)%radius_outer = r + h
    END DO
    CLOSE(iu)
  END SUBROUTINE read_particle_file

  SUBROUTINE next_data_line(iu, line, ios)
    INTEGER, INTENT(IN) :: iu
    CHARACTER(*), INTENT(OUT) :: line
    INTEGER, INTENT(OUT) :: ios
    DO
      READ(iu, '(A)', IOSTAT=ios) line
      IF (ios .NE. 0) RETURN
      line = ADJUSTL(line)
      IF (LEN_TRIM(line) .EQ. 0) CYCLE
      IF (line(1:1) .EQ. '#') CYCLE
      RETURN
    END DO
  END SUBROUTINE next_data_line

  SUBROUTINE lowercase(s)
    CHARACTER(*), INTENT(INOUT) :: s
    INTEGER :: i
    DO i = 1, LEN_TRIM(s)
      IF (s(i:i) .GE. 'A' .AND. s(i:i) .LE. 'Z') s(i:i) = CHAR(IACHAR(s(i:i)) + 32)
    END DO
  END SUBROUTINE lowercase

  !=======================================================================
  ! Coupling update (hook H1), workers only.
  !=======================================================================
  SUBROUTINE CHI_COUPLING_BEGIN_STEP(valU, valV, valW, valP)
    REAL*8, INTENT(IN) :: valU(*), valV(*), valW(*), valP(*)

    INTEGER :: lev, k, i, j, ip, a, b, nmissing, nf, ndof, nel, ic, ierr
    INTEGER, ALLOCATABLE :: qowner(:), dirmask(:), downer(:)
    REAL*8,  ALLOCATABLE :: qvals(:,:), hq(:,:,:), ubc(:,:), dvals(:,:)
    REAL*8 :: uq(3), g(3,3), pq, sig(3,3), n(3), un, resid, F(3), T(3), fac
    REAL*8 :: uval(3), gradu(3,3), pval
    LOGICAL :: ok
    CHARACTER(LEN=256) :: vtkname

    IF (.NOT. active) RETURN
    nstep_done = nstep_done + 1
    lev = NLMAX

    !---- background data at the Robin quadrature points -----------------
    ALLOCATE(qvals(CHI_BG_NVAL, MAX(nqtot,1)), qowner(MAX(nqtot,1)))
    IF (.NOT. dirichlet_mode) THEN
      CALL CHI_EXCHANGE_BG_EVAL(MPI_COMM_SUBS, nqtot, qpts, bgloc, &
        mg_mesh%level(lev)%dcorvg, mg_mesh%level(lev)%kvert, &
        mg_mesh%level(lev)%kedge, mg_mesh%level(lev)%karea, &
        mg_mesh%level(lev)%nvt, mg_mesh%level(lev)%net, mg_mesh%level(lev)%nat, &
        valU, valV, valW, valP, qvals, qowner, nmissing)
      IF (nmissing .GT. 0) THEN
        WRITE(*,'(A,I0,A)') 'CHI_COUPLING error: ', nmissing, &
          ' outer quadrature point(s) lie outside the background mesh'
        STOP 1
      END IF
    END IF

    !---- background velocity at outer nodes (Dirichlet mode) -------------
    ALLOCATE(dvals(CHI_BG_NVAL, MAX(ndtot,1)), downer(MAX(ndtot,1)))
    IF (dirichlet_mode) THEN
      CALL CHI_EXCHANGE_BG_EVAL(MPI_COMM_SUBS, ndtot, dpts, bgloc, &
        mg_mesh%level(lev)%dcorvg, mg_mesh%level(lev)%kvert, &
        mg_mesh%level(lev)%kedge, mg_mesh%level(lev)%karea, &
        mg_mesh%level(lev)%nvt, mg_mesh%level(lev)%net, mg_mesh%level(lev)%nat, &
        valU, valV, valW, valP, dvals, downer, nmissing)
      IF (nmissing .GT. 0) THEN
        WRITE(*,'(A,I0,A)') 'CHI_COUPLING error: ', nmissing, &
          ' outer node(s) lie outside the background mesh'
        STOP 1
      END IF
    END IF

    !---- per submesh: boundary data, owner solve, broadcast --------------
    DO k = 1, nsub
      nf = qoff(k+1) - qoff(k)
      ndof = atm(k)%ndof
      ALLOCATE(hq(3,9,MAX(nf,1)), dirmask(ndof), ubc(3,ndof))
      hq = 0d0
      dirmask = 0
      ubc = 0d0

      ! Robin data h = sigma.n - alpha (u.n) u, sigma = -p I + mu (G + G^T)
      IF (.NOT. dirichlet_mode) THEN
        DO j = 1, nf
          DO ip = 1, 9
            ic = 9*qoff(k) + 9*(j-1) + ip
            uq = qvals(1:3, ic)
            DO b = 1, 3
              g(:,b) = qvals(3+3*(b-1)+1:3+3*b, ic)   ! g(a,b) = du_a/dx_b
            END DO
            pq = qvals(13, ic)
            DO a = 1, 3
              DO b = 1, 3
                sig(a,b) = mu*(g(a,b) + g(b,a))
              END DO
              sig(a,a) = sig(a,a) - pq
            END DO
            n = qnrm(:, ic)
            un = uq(1)*n(1) + uq(2)*n(2) + uq(3)*n(3)
            hq(:,ip,j) = MATMUL(sig, n) - chimera_robin_alpha*un*uq
          END DO
        END DO
      END IF

      ! Dirichlet rows: inner surface = rigid-body velocity (static: 0);
      ! cylinder z faces = symmetry (w = 0 only); outer surface = background
      ! velocity in the diagnostic mode.
      DO i = 1, ndof
        IF (IAND(atm(k)%dofmask(i), CHI_SURF_INNER) .NE. 0) THEN
          dirmask(i) = CHI_DIR_ALL
          ubc(:,i) = 0d0
        ELSE IF (IAND(atm(k)%dofmask(i), CHI_SURF_ZLO+CHI_SURF_ZHI) .NE. 0) THEN
          dirmask(i) = CHI_DIR_W
          ubc(3,i) = 0d0
        END IF
      END DO
      IF (dirichlet_mode) THEN
        DO j = doff(k)+1, doff(k+1)
          i = ddof(j)
          IF (dirmask(i) .EQ. CHI_DIR_ALL) CYCLE      ! rim shared with inner
          dirmask(i) = CHI_DIR_ALL
          ubc(:,i) = dvals(1:3, j)
        END DO
      END IF

      IF (owner(k) .EQ. myrank) THEN
        ! Same time discretization as the background step (backward
        ! Euler, step TSTEP): the submesh problem is the time-discrete
        ! one, advanced from its own previous level.
        CALL CHI_SOLVE_STEADY_TAB(slv(k), atm(k), rho, mu, 1d0/TSTEP, &
          chimera_robin_alpha, dirichlet_mode, dirmask, ubc, hq, &
          chimera_sub_nl, resid, ok)
        IF (.NOT. ok) THEN
          WRITE(*,'(A,I0)') 'CHI_COUPLING error: submesh solve failed, body ', k
          STOP 1
        END IF
        WRITE(*,'(A,I0,A,I0,A,ES10.3)') 'Chimera: body ', k, ' solved on worker ', &
          myid, ', last Picard update = ', resid
      END IF
      nel = atm(k)%mesh%level(atm(k)%nlmax)%nel
      CALL CHI_EXCHANGE_BCAST(MPI_COMM_SUBS, owner(k), slv(k)%u, ndof)
      CALL CHI_EXCHANGE_BCAST(MPI_COMM_SUBS, owner(k), slv(k)%v, ndof)
      CALL CHI_EXCHANGE_BCAST(MPI_COMM_SUBS, owner(k), slv(k)%w, ndof)
      CALL CHI_EXCHANGE_BCAST(MPI_COMM_SUBS, owner(k), slv(k)%p, 4*nel)
      CALL CHI_SOLVER_ADVANCE(slv(k))     ! identical on all workers
      DEALLOCATE(hq, dirmask, ubc)
    END DO
    DEALLOCATE(qvals, qowner, dvals, downer)

    !---- weak variant: Dirichlet data at the penalty points -> g ---------
    IF (weak_mode) CALL update_penalty_rhs()

    !---- fringe values from the replicated submesh solutions -------------
    DO ic = 1, ncon
      i = con_dof(ic)
      IF (con_iel(ic) .EQ. 0) THEN            ! hole: rigid-body velocity
        fringeU(i) = 0d0; fringeV(i) = 0d0; fringeW(i) = 0d0
        CYCLE
      END IF
      k = con_sub(ic)
      lev = atm(k)%nlmax
      CALL CHI_EVAL_FIELD_AT(con_iel(ic), con_xi(:,ic), &
        atm(k)%mesh%level(lev)%kvert, atm(k)%mesh%level(lev)%kedge, &
        atm(k)%mesh%level(lev)%karea, atm(k)%mesh%level(lev)%nvt, &
        atm(k)%mesh%level(lev)%net, atm(k)%mesh%level(lev)%nat, &
        atm(k)%mesh%level(lev)%dcorvg, slv(k)%u, slv(k)%v, slv(k)%w, slv(k)%p, &
        uval, gradu, pval, ok)
      IF (chimera_coupling_relax .LT. 1d0 .AND. fringe_valid) THEN
        fringeU(i) = chimera_coupling_relax*uval(1) + (1d0-chimera_coupling_relax)*fringeU(i)
        fringeV(i) = chimera_coupling_relax*uval(2) + (1d0-chimera_coupling_relax)*fringeV(i)
        fringeW(i) = chimera_coupling_relax*uval(3) + (1d0-chimera_coupling_relax)*fringeW(i)
      ELSE
        fringeU(i) = uval(1); fringeV(i) = uval(2); fringeW(i) = uval(3)
      END IF
    END DO
    fringe_valid = (ncon .GT. 0)

    !---- forces (design section 5), reported once ------------------------
    fac = 2d0/(postParams%U_mean*postParams%U_mean*postParams%D*postParams%H)
    DO k = 1, nsub
      lev = atm(k)%nlmax
      CALL CHI_COMPUTE_FORCES(atm(k)%innerFaces, SIZE(atm(k)%innerFaces,2), &
        atm(k)%mesh%level(lev)%nvt, atm(k)%mesh%level(lev)%net, &
        atm(k)%mesh%level(lev)%nat, atm(k)%mesh%level(lev)%nel, &
        atm(k)%mesh%level(lev)%kvert, atm(k)%mesh%level(lev)%kedge, &
        atm(k)%mesh%level(lev)%karea, atm(k)%mesh%level(lev)%dcorvg, &
        slv(k)%u, slv(k)%v, slv(k)%w, slv(k)%p, mu, atm(k)%center, F, T)
      IF (myid .EQ. showid) THEN
        WRITE(mfile_unit,'(A,I0,A,7ES15.7E2)') 'ChimeraForce', k, ': ', &
          timens, fac*F(1), fac*F(2), F(1), F(2), F(3), T(3)
        WRITE(*,'(A,I0,A,7ES15.7E2)') 'ChimeraForce', k, ': ', &
          timens, fac*F(1), fac*F(2), F(1), F(2), F(3), T(3)
        IF (chimera_write_vtk) THEN
          WRITE(vtkname,'(A,I0,A,I0,A)') '_vtk/chimera_body', k, '_', nstep_done, '.vtk'
          CALL CHI_WRITE_SUBMESH_VTK(TRIM(vtkname), atm(k)%mesh%level(lev)%nvt, &
            atm(k)%mesh%level(lev)%nel, atm(k)%mesh%level(lev)%kvert, &
            atm(k)%mesh%level(lev)%dcorvg, slv(k)%u, slv(k)%v, slv(k)%w, slv(k)%p)
        END IF
      END IF
    END DO
  END SUBROUTINE CHI_COUPLING_BEGIN_STEP

  !=======================================================================
  ! Hook bodies (H2-H4): siblings of the FictKNPR branches.
  !=======================================================================
  SUBROUTINE CHI_COUPLING_APPLY_DEF(defU, defV, defW, ndof)
    REAL*8, INTENT(INOUT) :: defU(*), defV(*), defW(*)
    INTEGER, INTENT(IN) :: ndof
    INTEGER :: i
    IF (.NOT. active) RETURN
    DO i = 1, MIN(ndof, nbgdof)
      IF (marker_kind(i) .NE. CHI_MARK_FREE) THEN
        defU(i) = 0d0
        defV(i) = 0d0
        defW(i) = 0d0
      END IF
    END DO
  END SUBROUTINE CHI_COUPLING_APPLY_DEF

  SUBROUTINE CHI_COUPLING_APPLY_VAL(valU, valV, valW, ndof)
    REAL*8, INTENT(INOUT) :: valU(*), valV(*), valW(*)
    INTEGER, INTENT(IN) :: ndof
    INTEGER :: i
    IF (.NOT. active) RETURN
    DO i = 1, MIN(ndof, nbgdof)
      IF (marker_kind(i) .NE. CHI_MARK_FREE) THEN
        valU(i) = fringeU(i)
        valV(i) = fringeV(i)
        valW(i) = fringeW(i)
      END IF
    END DO
  END SUBROUTINE CHI_COUPLING_APPLY_VAL

  SUBROUTINE CHI_COUPLING_FILTER_MAT(DA11, DA22, DA33, KLD, ndof)
    REAL*8, INTENT(INOUT) :: DA11(*), DA22(*), DA33(*)
    INTEGER, INTENT(IN) :: KLD(*), ndof
    INTEGER :: i, icol
    IF (.NOT. active) RETURN
    DO i = 1, MIN(ndof, nbgdof)
      IF (marker_kind(i) .EQ. CHI_MARK_FREE) CYCLE
      DO icol = KLD(i)+1, KLD(i+1)-1
        DA11(icol) = 0d0
        DA22(icol) = 0d0
        DA33(icol) = 0d0
      END DO
    END DO
  END SUBROUTINE CHI_COUPLING_FILTER_MAT

  SUBROUTINE CHI_COUPLING_FILTER_MAT_9(DA11, DA22, DA33, DA12, DA13, DA23, &
                                       DA21, DA31, DA32, KLD, ndof)
    REAL*8, INTENT(INOUT) :: DA11(*), DA22(*), DA33(*), DA12(*), DA13(*), &
                             DA23(*), DA21(*), DA31(*), DA32(*)
    INTEGER, INTENT(IN) :: KLD(*), ndof
    INTEGER :: i, icol
    IF (.NOT. active) RETURN
    DO i = 1, MIN(ndof, nbgdof)
      IF (marker_kind(i) .EQ. CHI_MARK_FREE) CYCLE
      icol = KLD(i)
      DA12(icol) = 0d0; DA13(icol) = 0d0; DA23(icol) = 0d0
      DA21(icol) = 0d0; DA31(icol) = 0d0; DA32(icol) = 0d0
      DO icol = KLD(i)+1, KLD(i+1)-1
        DA11(icol) = 0d0; DA22(icol) = 0d0; DA33(icol) = 0d0
        DA12(icol) = 0d0; DA13(icol) = 0d0; DA23(icol) = 0d0
        DA21(icol) = 0d0; DA31(icol) = 0d0; DA32(icol) = 0d0
      END DO
    END DO
  END SUBROUTINE CHI_COUPLING_FILTER_MAT_9

  !=======================================================================
  ! Release (hook H6); idempotent, partial-init-safe.
  !=======================================================================
  SUBROUTINE CHI_COUPLING_FINALIZE()
    INTEGER :: k
    IF (ALLOCATED(slv)) THEN
      DO k = 1, SIZE(slv)
        CALL CHI_SOLVER_RELEASE(slv(k))
      END DO
      DEALLOCATE(slv)
    END IF
    IF (ALLOCATED(subloc)) THEN
      DO k = 1, SIZE(subloc)
        CALL CHI_LOCATOR_RELEASE(subloc(k))
      END DO
      DEALLOCATE(subloc)
    END IF
    IF (ALLOCATED(atm)) THEN
      DO k = 1, SIZE(atm)
        CALL CHI_RELEASE_SUBMESH_MESH(atm(k))
        CALL CHI_SUBMESH_RELEASE(atm(k))
      END DO
      DEALLOCATE(atm)
    END IF
    IF (ALLOCATED(bodies)) DEALLOCATE(bodies)
    IF (ALLOCATED(owner)) DEALLOCATE(owner)
    CALL CHI_LOCATOR_RELEASE(bgloc)
    IF (ALLOCATED(marker_kind)) DEALLOCATE(marker_kind)
    IF (ALLOCATED(marker_pid)) DEALLOCATE(marker_pid)
    IF (ALLOCATED(fringeU)) DEALLOCATE(fringeU, fringeV, fringeW)
    IF (ALLOCATED(con_dof)) DEALLOCATE(con_dof, con_sub, con_iel, con_xi)
    IF (ALLOCATED(qoff)) DEALLOCATE(qoff)
    IF (ALLOCATED(qpts)) DEALLOCATE(qpts, qnrm)
    IF (ALLOCATED(doff)) DEALLOCATE(doff)
    IF (ALLOCATED(ddof)) DEALLOCATE(ddof, dpts)
    CALL release_penalty()
    nsub = 0
    ncon = 0
    nqtot = 0
    ndtot = 0
    nbgdof = 0
    active = .FALSE.
  END SUBROUTINE CHI_COUPLING_FINALIZE

  !=======================================================================
  ! Weak variant (Chimera-W): penalty operator setup and hook bodies
  ! (design sections 3 and 5; H8/H10/H11).
  !=======================================================================
  SUBROUTINE setup_penalty()
    INTEGER :: lev, k, a, q, na, nact, nmiss, npts, ilev_save, sl
    REAL*8, ALLOCATABLE :: hw(:), xn(:,:)
    INTEGER, ALLOCATABLE :: nb_of(:)
    LOGICAL, ALLOCATABLE :: inb(:)
    REAL*8 :: excess, maxexcess
    LOGICAL :: found

    nact = 0
    npts = 0
    maxexcess = 0d0
    pen_lmin = NLMIN
    pen_lmax = NLMAX
    ALLOCATE(ptab(pen_lmin:pen_lmax), pmat(pen_lmin:pen_lmax), hw(nsub))
    DO k = 1, nsub
      hw(k) = atm(k)%radius_outer - atm(k)%radius_inner
    END DO
    lumped = chimera_penalty_lumped
    DO lev = pen_lmin, pen_lmax
      na = mg_qMat(lev)%LdA(mg_qMat(lev)%nu+1) - 1
      ALLOCATE(pmat(lev)%dl(mg_qMat(lev)%nu), pmat(lev)%dlg(mg_qMat(lev)%nu))
      IF (lumped) THEN
        ! nodal-quadrature lumped penalty (positive, constant-exact)
        ALLOCATE(nb_of(mg_qMat(lev)%nu), inb(mg_qMat(lev)%nu), xn(3, mg_qMat(lev)%nu))
        CALL CHI_PENALTY_NODAL(mg_mesh%level(lev)%nel, mg_mesh%level(lev)%nvt, &
          mg_mesh%level(lev)%net, mg_mesh%level(lev)%nat, mg_mesh%level(lev)%kvert, &
          mg_mesh%level(lev)%kedge, mg_mesh%level(lev)%karea, mg_mesh%level(lev)%dcorvg, &
          nsub, bodies, hw, chimera_gamma_max, mg_qMat(lev)%nu, pmat(lev)%dl, &
          nb_of, inb, xn, chimera_beta_full, chimera_beta_zero)
        IF (lev .EQ. NLMAX) THEN
          ! finest level: penalised node list with donors in the atmosphere
          nnod = COUNT(pmat(lev)%dl .GT. 0d0)
          ALLOCATE(nd_dof(MAX(nnod,1)), nd_sub(MAX(nnod,1)), nd_iel(MAX(nnod,1)), &
                   nd_xi(3, MAX(nnod,1)))
          nnod = 0
          npts = 0
          nmiss = 0
          maxexcess = 0d0
          DO a = 1, mg_qMat(lev)%nu
            IF (pmat(lev)%dl(a) .LE. 0d0) CYCLE
            nnod = nnod + 1
            nd_dof(nnod) = a
            nd_sub(nnod) = nb_of(a)
            nd_iel(nnod) = 0
            nd_xi(:,nnod) = 0d0
            IF (inb(a)) CYCLE                       ! hole node: rigid velocity
            k = nb_of(a)
            sl = atm(k)%nlmax
            npts = npts + 1
            CALL CHI_LOCATE(subloc(k), atm(k)%mesh%level(sl)%dcorvg, &
              atm(k)%mesh%level(sl)%kvert, xn(:,a), nd_iel(nnod), nd_xi(:,nnod), found)
            IF (.NOT. found) THEN
              CALL CHI_LOCATE_NEAREST(subloc(k), atm(k)%mesh%level(sl)%dcorvg, &
                atm(k)%mesh%level(sl)%kvert, xn(:,a), nd_iel(nnod), nd_xi(:,nnod), &
                excess, found)
              IF (found) THEN
                maxexcess = MAX(maxexcess, excess)
              ELSE
                nmiss = nmiss + 1
              END IF
            END IF
          END DO
          IF (nmiss .GT. 0) THEN
            WRITE(*,'(A,I0,A)') 'CHI_COUPLING error: ', nmiss, &
              ' penalised node(s) have no donor element in their atmosphere'
            STOP 1
          END IF
        END IF
        DEALLOCATE(nb_of, inb, xn)
      ELSE
        ! consistent penalty matrix on the level pattern (paper scheme)
        CALL CHI_PENALTY_TABULATE(mg_mesh%level(lev)%nel, mg_mesh%level(lev)%kvert, &
          mg_mesh%level(lev)%dcorvg, nsub, bodies, hw, ptab(lev), &
          chimera_beta_full, chimera_beta_zero)
        ALLOCATE(pmat(lev)%d(na))
        CALL CHI_PENALTY_ASSEMBLE_D(ptab(lev), chimera_gamma_max, &
          mg_mesh%level(lev)%kvert, mg_mesh%level(lev)%kedge, mg_mesh%level(lev)%karea, &
          mg_mesh%level(lev)%nvt, mg_mesh%level(lev)%net, mg_mesh%level(lev)%nat, &
          mg_qMat(lev)%LdA, mg_qMat(lev)%ColA, mg_qMat(lev)%nu, pmat(lev)%d)
        CALL CHI_PENALTY_LUMP(mg_qMat(lev)%LdA, mg_qMat(lev)%nu, pmat(lev)%d, pmat(lev)%dl)
      END IF
      pmat(lev)%dlg = pmat(lev)%dl
      ilev_save = ILEV
      ILEV = lev
      CALL E013Sum(pmat(lev)%dlg)
      ILEV = ilev_save
    END DO
    DEALLOCATE(hw)
    lev = NLMAX
    IF (mg_qMat(lev)%nu .NE. nbgdof) THEN
      WRITE(*,'(A)') 'CHI_COUPLING error: Q2 pattern size differs from the dof count'
      STOP 1
    END IF
    ALLOCATE(gU(nbgdof), gV(nbgdof), gW(nbgdof))
    gU = 0d0; gV = 0d0; gW = 0d0
    g_valid = .FALSE.
    IF (.NOT. lumped) THEN
      ! donors of the finest-level points with beta > 0 outside the bodies
      nact = ptab(lev)%nact
      ALLOCATE(pd_iel(CHI_PENALTY_NQ, MAX(nact,1)), pd_xi(3, CHI_PENALTY_NQ, MAX(nact,1)))
      pd_iel = 0
      pd_xi = 0d0
      nmiss = 0
      npts = 0
      maxexcess = 0d0
      DO a = 1, nact
        DO q = 1, CHI_PENALTY_NQ
          IF (ptab(lev)%w(q,a) .EQ. 0d0 .OR. ptab(lev)%inbody(q,a)) CYCLE
          k = ptab(lev)%body(q,a)
          sl = atm(k)%nlmax
          npts = npts + 1
          CALL CHI_LOCATE(subloc(k), atm(k)%mesh%level(sl)%dcorvg, &
            atm(k)%mesh%level(sl)%kvert, ptab(lev)%x(:,q,a), pd_iel(q,a), pd_xi(:,q,a), found)
          IF (.NOT. found) THEN
            CALL CHI_LOCATE_NEAREST(subloc(k), atm(k)%mesh%level(sl)%dcorvg, &
              atm(k)%mesh%level(sl)%kvert, ptab(lev)%x(:,q,a), pd_iel(q,a), pd_xi(:,q,a), &
              excess, found)
            IF (found) THEN
              maxexcess = MAX(maxexcess, excess)
            ELSE
              nmiss = nmiss + 1
            END IF
          END IF
        END DO
      END DO
      IF (nmiss .GT. 0) THEN
        WRITE(*,'(A,I0,A)') 'CHI_COUPLING error: ', nmiss, &
          ' penalty quadrature point(s) have no donor element in their atmosphere'
        STOP 1
      END IF
      ALLOCATE(pdiag(nbgdof), pwts(nbgdof))
      ! global diagonal of D and partition-of-unity weights (design section 5)
      CALL CHI_PENALTY_DIAG(mg_qMat(lev)%LdA, nbgdof, pmat(lev)%d, pdiag)
      pwts = 1d0
      ilev_save = ILEV
      ILEV = NLMAX
      CALL E013Sum(pdiag)
      CALL E013Sum(pwts)
      ILEV = ilev_save
      pwts = 1d0/pwts
    END IF
    IF (myid .EQ. showid) THEN
      IF (lumped) THEN
        WRITE(*,'(A,I0,A,I0,A,ES10.2)') 'Chimera: penalised nodes (finest) = ', &
          nnod, ', donor nodes = ', npts, ', max donor extrapolation = ', maxexcess
      ELSE
        WRITE(*,'(A,I0,A,I0,A,ES10.2)') 'Chimera: penalty active elements (finest) = ', &
          nact, ', donor points = ', npts, ', max donor extrapolation = ', maxexcess
      END IF
      WRITE(mfile_unit,'(A,ES10.2,A,I0,A,I0,A,L1)') 'Chimera: weak coupling gamma_max = ', &
        chimera_gamma_max, ', penalty levels = ', pen_lmin, '..', pen_lmax, &
        ', lumped = ', lumped
    END IF
  END SUBROUTINE setup_penalty

  SUBROUTINE release_penalty()
    INTEGER :: lev
    IF (ALLOCATED(ptab)) THEN
      DO lev = LBOUND(ptab,1), UBOUND(ptab,1)
        CALL CHI_PENALTY_RELEASE(ptab(lev))
      END DO
      DEALLOCATE(ptab)
    END IF
    IF (ALLOCATED(pmat)) DEALLOCATE(pmat)
    IF (ALLOCATED(pd_iel)) DEALLOCATE(pd_iel, pd_xi)
    IF (ALLOCATED(nd_dof)) DEALLOCATE(nd_dof, nd_sub, nd_iel, nd_xi)
    nnod = 0
    IF (ALLOCATED(gU)) DEALLOCATE(gU, gV, gW)
    IF (ALLOCATED(gU_prev)) DEALLOCATE(gU_prev, gV_prev, gW_prev)
    IF (ALLOCATED(pdiag)) DEALLOCATE(pdiag, pwts)
    weak_mode = .FALSE.
    g_valid = .FALSE.
    pen_lmin = 0
    pen_lmax = -1
    pen_dt_cmat = -1d0
  END SUBROUTINE release_penalty

  ! Dirichlet data of the penalty at the finest-level points: rigid
  ! velocity inside a body (static: 0), the replicated submesh solution
  ! in the atmosphere; then g.
  SUBROUTINE update_penalty_rhs()
    INTEGER :: lev, a, q, k, sl, i
    REAL*8 :: uval(3), gradu(3,3), pval
    LOGICAL :: ok
    lev = NLMAX
    IF (lumped) THEN
      ! g_i = D_L(i) * uhat(x_i) (rank-partial D_L, like every assembled rhs)
      gU = 0d0; gV = 0d0; gW = 0d0
      DO a = 1, nnod
        i = nd_dof(a)
        IF (nd_iel(a) .EQ. 0) CYCLE               ! hole node: uhat = 0 (static)
        k = nd_sub(a)
        sl = atm(k)%nlmax
        CALL CHI_EVAL_FIELD_AT(nd_iel(a), nd_xi(:,a), &
          atm(k)%mesh%level(sl)%kvert, atm(k)%mesh%level(sl)%kedge, &
          atm(k)%mesh%level(sl)%karea, atm(k)%mesh%level(sl)%nvt, &
          atm(k)%mesh%level(sl)%net, atm(k)%mesh%level(sl)%nat, &
          atm(k)%mesh%level(sl)%dcorvg, slv(k)%u, slv(k)%v, slv(k)%w, slv(k)%p, &
          uval, gradu, pval, ok)
        gU(i) = pmat(lev)%dl(i)*uval(1)
        gV(i) = pmat(lev)%dl(i)*uval(2)
        gW(i) = pmat(lev)%dl(i)*uval(3)
      END DO
      CALL relax_g()
      g_valid = .TRUE.
      RETURN
    END IF
    DO a = 1, ptab(lev)%nact
      DO q = 1, CHI_PENALTY_NQ
        IF (ptab(lev)%w(q,a) .EQ. 0d0) CYCLE
        IF (ptab(lev)%inbody(q,a)) THEN
          ptab(lev)%uhat(:,q,a) = 0d0
          CYCLE
        END IF
        k = ptab(lev)%body(q,a)
        sl = atm(k)%nlmax
        CALL CHI_EVAL_FIELD_AT(pd_iel(q,a), pd_xi(:,q,a), &
          atm(k)%mesh%level(sl)%kvert, atm(k)%mesh%level(sl)%kedge, &
          atm(k)%mesh%level(sl)%karea, atm(k)%mesh%level(sl)%nvt, &
          atm(k)%mesh%level(sl)%net, atm(k)%mesh%level(sl)%nat, &
          atm(k)%mesh%level(sl)%dcorvg, slv(k)%u, slv(k)%v, slv(k)%w, slv(k)%p, &
          uval, gradu, pval, ok)
        ptab(lev)%uhat(:,q,a) = uval
      END DO
    END DO
    CALL CHI_PENALTY_ASSEMBLE_G(ptab(lev), chimera_gamma_max, &
      mg_mesh%level(lev)%kvert, mg_mesh%level(lev)%kedge, mg_mesh%level(lev)%karea, &
      mg_mesh%level(lev)%nvt, mg_mesh%level(lev)%net, mg_mesh%level(lev)%nat, &
      nbgdof, gU, gV, gW)
    CALL relax_g()
    g_valid = .TRUE.
  END SUBROUTINE update_penalty_rhs

  ! Coupling under-relaxation of g (ChimeraCouplingRelax); gU_prev holds
  ! the previously applied data (first call: none -> full update).
  SUBROUTINE relax_g()
    REAL*8 :: th
    th = chimera_coupling_relax
    IF (th .LT. 1d0) THEN
      IF (.NOT. ALLOCATED(gU_prev)) THEN
        ALLOCATE(gU_prev(nbgdof), gV_prev(nbgdof), gW_prev(nbgdof))
      ELSE
        gU = th*gU + (1d0-th)*gU_prev
        gV = th*gV + (1d0-th)*gV_prev
        gW = th*gW + (1d0-th)*gW_prev
      END IF
      gU_prev = gU; gV_prev = gV; gW_prev = gW
    END IF
  END SUBROUTINE relax_g

  ! H8: A11/A22/A33 += coef * D on level ilev (coef = tstep, design section 5).
  SUBROUTINE CHI_COUPLING_ADD_MAT(DA11, DA22, DA33, KLD, nu, ilev, coef)
    REAL*8, INTENT(INOUT) :: DA11(*), DA22(*), DA33(*)
    INTEGER, INTENT(IN) :: KLD(*), nu, ilev
    REAL*8, INTENT(IN) :: coef
    INTEGER :: na, p
    IF (.NOT. (active .AND. weak_mode)) RETURN
    IF (ilev .LT. pen_lmin .OR. ilev .GT. pen_lmax) RETURN
    na = KLD(nu+1) - 1
    IF (nu .NE. SIZE(pmat(ilev)%dl)) THEN
      WRITE(*,'(A,I0)') 'CHI_COUPLING error: penalty dof count mismatch on level ', ilev
      STOP 1
    END IF
    IF (.NOT. lumped) THEN
      IF (na .NE. SIZE(pmat(ilev)%d)) THEN
        WRITE(*,'(A,I0)') 'CHI_COUPLING error: penalty pattern mismatch on level ', ilev
        STOP 1
      END IF
    END IF
    IF (lumped) THEN
      DO p = 1, nu
        DA11(KLD(p)) = DA11(KLD(p)) + coef*pmat(ilev)%dl(p)
        DA22(KLD(p)) = DA22(KLD(p)) + coef*pmat(ilev)%dl(p)
        DA33(KLD(p)) = DA33(KLD(p)) + coef*pmat(ilev)%dl(p)
      END DO
    ELSE
      DO p = 1, na
        DA11(p) = DA11(p) + coef*pmat(ilev)%d(p)
        DA22(p) = DA22(p) + coef*pmat(ilev)%d(p)
        DA33(p) = DA33(p) + coef*pmat(ilev)%d(p)
      END DO
    END IF
  END SUBROUTINE CHI_COUPLING_ADD_MAT

  ! Defect contribution def += coef * D u on the finest level (the branches
  ! of Matdef that rebuild the defect from parts instead of the assembled A).
  SUBROUTINE CHI_COUPLING_ADD_DEFECT(valU, valV, valW, defU, defV, defW, &
                                     KLD, KCOL, nu, coef)
    REAL*8, INTENT(IN) :: valU(*), valV(*), valW(*)
    REAL*8, INTENT(INOUT) :: defU(*), defV(*), defW(*)
    INTEGER, INTENT(IN) :: KLD(*), KCOL(*), nu
    REAL*8, INTENT(IN) :: coef
    INTEGER :: i
    IF (.NOT. (active .AND. weak_mode)) RETURN
    IF (nu .NE. nbgdof) RETURN
    IF (lumped) THEN
      DO i = 1, nu
        defU(i) = defU(i) + coef*pmat(NLMAX)%dl(i)*valU(i)
        defV(i) = defV(i) + coef*pmat(NLMAX)%dl(i)*valV(i)
        defW(i) = defW(i) + coef*pmat(NLMAX)%dl(i)*valW(i)
      END DO
    ELSE
      CALL CHI_PENALTY_MATVEC(KLD, KCOL, nu, pmat(NLMAX)%d, coef, valU(1:nu), defU(1:nu))
      CALL CHI_PENALTY_MATVEC(KLD, KCOL, nu, pmat(NLMAX)%d, coef, valV(1:nu), defV(1:nu))
      CALL CHI_PENALTY_MATVEC(KLD, KCOL, nu, pmat(NLMAX)%d, coef, valW(1:nu), defW(1:nu))
    END IF
  END SUBROUTINE CHI_COUPLING_ADD_DEFECT

  ! H10: rhs += coef * g (rank-partial, like every assembled right-hand side).
  SUBROUTINE CHI_COUPLING_ADD_RHS(defU, defV, defW, ndof, coef)
    REAL*8, INTENT(INOUT) :: defU(*), defV(*), defW(*)
    INTEGER, INTENT(IN) :: ndof
    REAL*8, INTENT(IN) :: coef
    INTEGER :: i
    IF (.NOT. (active .AND. weak_mode .AND. g_valid)) RETURN
    DO i = 1, MIN(ndof, nbgdof)
      defU(i) = defU(i) + coef*gU(i)
      defV(i) = defV(i) + coef*gV(i)
      defW(i) = defW(i) + coef*gW(i)
    END DO
  END SUBROUTINE CHI_COUPLING_ADD_RHS

  ! H11, paper eq. (12): [ML + tstep D] delta = def, u = u~ - delta.
  ! def is the assembly-summed, Dirichlet-filtered dt*B(dp); ml the global
  ! lumped mass of the correction; filter3 the caller's defect filter.
  ! applied = .FALSE. leaves the caller's diagonal loop (fast path).
  SUBROUTINE CHI_COUPLING_CORRECT(valU, valV, valW, defU, defV, defW, ml, ndof, &
                                  filter3, applied)
    REAL*8, INTENT(INOUT) :: valU(*), valV(*), valW(*)
    REAL*8, INTENT(IN) :: defU(*), defV(*), defW(*), ml(*)
    INTEGER, INTENT(IN) :: ndof
    PROCEDURE(chi_filter3_iface) :: filter3
    LOGICAL, INTENT(OUT) :: applied
    REAL*8, ALLOCATABLE :: rhs(:,:), x(:,:), dglob(:)
    INTEGER :: n, i
    REAL*8 :: meff

    applied = .FALSE.
    IF (.NOT. (active .AND. weak_mode)) RETURN
    IF (ndof .NE. nbgdof) THEN
      WRITE(*,'(A)') 'CHI_COUPLING error: correction called with a foreign dof count'
      STOP 1
    END IF
    n = nbgdof
    IF (lumped) THEN
      ! diagonal system [M_L + dt D_L]; the projection operator (hook H15)
      ! was built with the same dt, which the weak variant requires to be
      ! constant.
      IF (pen_dt_cmat .GT. 0d0 .AND. ABS(TSTEP - pen_dt_cmat) .GT. 1d-12*pen_dt_cmat) THEN
        WRITE(*,'(A)') 'CHI_COUPLING error: the time step changed after the projection ' // &
          'operator was built; the weak variant requires a constant TimeStep.'
        STOP 1
      END IF
      ! same capped operator as the projection (hook H15): exactly
      ! discretely divergence-free corrected velocity
      DO i = 1, n
        meff = ml(i) + MIN(TSTEP*pmat(NLMAX)%dlg(i), chimera_proj_cap*ml(i))
        valU(i) = valU(i) - defU(i)/meff
        valV(i) = valV(i) - defV(i)/meff
        valW(i) = valW(i) - defW(i)/meff
      END DO
      applied = .TRUE.
      RETURN
    END IF
    ALLOCATE(rhs(n,3), x(n,3), dglob(n))
    rhs(:,1) = defU(1:n); rhs(:,2) = defV(1:n); rhs(:,3) = defW(1:n)
    dglob = ml(1:n) + TSTEP*pdiag
    CALL CHI_PCG3(n, mg_qMat(NLMAX)%LdA, mg_qMat(NLMAX)%ColA, pmat(NLMAX)%d, TSTEP, &
      ml(1:n), dglob, pwts, rhs(:,1), rhs(:,2), rhs(:,3), x(:,1), x(:,2), x(:,3), &
      1d-10, 500, pen_sum3, filter3, pen_allsum, ncorr_it, corr_resid)
    valU(1:n) = valU(1:n) - x(:,1)
    valV(1:n) = valV(1:n) - x(:,2)
    valW(1:n) = valW(1:n) - x(:,3)
    DEALLOCATE(rhs, x, dglob)
    applied = .TRUE.
    IF (myid .EQ. showid) WRITE(*,'(A,I0,A,ES10.3)') &
      'Chimera: correction CG iterations = ', ncorr_it, ', rel. residual = ', corr_resid
  END SUBROUTINE CHI_COUPLING_CORRECT

  ! H15: penalised lumped mass of the projection operator on level ilev
  ! (lumped variant only; the consistent variant keeps the paper's M_L).
  SUBROUTINE CHI_COUPLING_ADD_PRESSURE_MASS(meff, nu, ilev, coef)
    REAL*8, INTENT(INOUT) :: meff(*)
    INTEGER, INTENT(IN) :: nu, ilev
    REAL*8, INTENT(IN) :: coef
    INTEGER :: i
    IF (.NOT. (active .AND. weak_mode .AND. lumped)) RETURN
    IF (ilev .LT. pen_lmin .OR. ilev .GT. pen_lmax) RETURN
    IF (nu .NE. SIZE(pmat(ilev)%dlg)) THEN
      WRITE(*,'(A,I0)') 'CHI_COUPLING error: projection dof count mismatch on level ', ilev
      STOP 1
    END IF
    DO i = 1, nu
      meff(i) = meff(i) + MIN(coef*pmat(ilev)%dlg(i), chimera_proj_cap*meff(i))
    END DO
    IF (ilev .EQ. NLMAX) pen_dt_cmat = coef
  END SUBROUTINE CHI_COUPLING_ADD_PRESSURE_MASS

  SUBROUTINE pen_sum3(y1, y2, y3)
    REAL*8, INTENT(INOUT) :: y1(*), y2(*), y3(*)
    CALL E013Sum3(y1, y2, y3)
  END SUBROUTINE pen_sum3

  SUBROUTINE pen_allsum(vals, n)
    INTEGER, INTENT(IN) :: n
    REAL*8, INTENT(INOUT) :: vals(n)
    CALL CHI_EXCHANGE_ALLSUM(MPI_COMM_SUBS, vals, n)
  END SUBROUTINE pen_allsum

  !=======================================================================
  ! Restart (hook H12, app-driven): the replicated submesh states are
  ! written once (first worker) beside the flow dump and read by every
  ! worker; donor caches, markers and g are rebuilt, never restored.
  !=======================================================================
  SUBROUTINE CHI_COUPLING_WRITE_RESTART(idx)
    INTEGER, INTENT(IN) :: idx
    CHARACTER(LEN=64) :: dir
    CHARACTER(LEN=20) :: tag
    INTEGER :: iu, k, ios
    IF (.NOT. active) RETURN
    IF (myid .NE. 1) RETURN
    WRITE(dir,'(A,I0)') '_dump/', idx
    CALL EXECUTE_COMMAND_LINE('mkdir -p ' // TRIM(dir), WAIT=.TRUE.)
    iu = 4711
    OPEN(iu, FILE=TRIM(dir)//'/chimera.dmp', FORM='UNFORMATTED', STATUS='REPLACE', &
      ACTION='WRITE', IOSTAT=ios)
    IF (ios .NE. 0) THEN
      WRITE(*,'(A)') 'CHI_COUPLING error: cannot write ' // TRIM(dir) // '/chimera.dmp'
      STOP 1
    END IF
    tag = 'CHIMERA_RESTART_V1'
    WRITE(iu) tag
    WRITE(iu) nsub, nstep_done, TIMENS
    DO k = 1, nsub
      WRITE(iu) atm(k)%ndof, atm(k)%mesh%level(atm(k)%nlmax)%nel
      WRITE(iu) slv(k)%u, slv(k)%v, slv(k)%w
      WRITE(iu) slv(k)%uold, slv(k)%vold, slv(k)%wold
      WRITE(iu) slv(k)%p
    END DO
    CLOSE(iu)
    WRITE(*,'(A)') 'Chimera: restart state written to ' // TRIM(dir) // '/chimera.dmp'
  END SUBROUTINE CHI_COUPLING_WRITE_RESTART

  SUBROUTINE CHI_COUPLING_READ_RESTART(name)
    CHARACTER(LEN=*), INTENT(IN) :: name
    CHARACTER(LEN=256) :: fname
    CHARACTER(LEN=20) :: tag
    INTEGER :: iu, k, ios, ns, nd, ne
    REAL*8 :: t
    IF (.NOT. active) RETURN
    fname = '_dump/' // TRIM(ADJUSTL(name)) // '/chimera.dmp'
    iu = 4711
    OPEN(iu, FILE=TRIM(fname), FORM='UNFORMATTED', STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios .NE. 0) THEN
      WRITE(*,'(A)') 'CHI_COUPLING error: restart requested but ' // TRIM(fname) // &
        ' is missing'
      STOP 1
    END IF
    READ(iu) tag
    IF (tag .NE. 'CHIMERA_RESTART_V1') THEN
      WRITE(*,'(A)') 'CHI_COUPLING error: unknown restart format in ' // TRIM(fname)
      STOP 1
    END IF
    READ(iu) ns, nstep_done, t
    IF (ns .NE. nsub) THEN
      WRITE(*,'(A)') 'CHI_COUPLING error: restart body count differs from the body table'
      STOP 1
    END IF
    DO k = 1, nsub
      READ(iu) nd, ne
      IF (nd .NE. atm(k)%ndof .OR. ne .NE. atm(k)%mesh%level(atm(k)%nlmax)%nel) THEN
        WRITE(*,'(A,I0)') 'CHI_COUPLING error: restart submesh size differs, body ', k
        STOP 1
      END IF
      READ(iu) slv(k)%u, slv(k)%v, slv(k)%w
      READ(iu) slv(k)%uold, slv(k)%vold, slv(k)%wold
      READ(iu) slv(k)%p
    END DO
    CLOSE(iu)
    IF (myid .EQ. showid) WRITE(*,'(A,I0,A,ES12.5)') 'Chimera: restart state read from ' // &
      TRIM(fname) // ', steps done = ', nstep_done, ', dump time = ', t
  END SUBROUTINE CHI_COUPLING_READ_RESTART

END MODULE CHI_COUPLING

