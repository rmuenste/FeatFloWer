PROGRAM TEST_EL_PERIODIC_MPI

  ! Distributed test of the periodic image treatment of the particle-mesh coupling
  ! (SimPar@ELPeriodicKernel = wrap) on 4 or 8 ranks.
  !
  ! The fluid cell is a cube of side L = 0.02 with n = 8 cells per axis, periodic in x and y,
  ! walls on the z planes. Every rank holds one block of a px x py partition (2 x 2 on four
  ! ranks, 4 x 2 on eight), all workers, as the workers of the application; the master of the
  ! application has no part in the coupling. The test runs the routines of the application:
  !
  !   EL_REFRESH_COUPLING_HALO   neighbour graph (MPI_Dist_graph_create_adjacent) and record
  !                              exchange, with the lattice active
  !   EL_SAMPLE_RECORD           kernel sample of every record on every rank
  !   EL_SAMPLE_RECORD_VISCOUS   additive pieces of the kernel-averaged Laplacian
  !   EL_REDUCE_TO_OWNERS        sum of the pieces at the owner
  !   EL_BROADCAST_FROM_OWNERS   normalisation back to every rank that holds the record
  !   EL_DEPOSIT_RECORD          solid volume, its gradient, reaction force, drag coefficient
  !
  ! and compares with the same particles coupled on one rank that holds the whole cell (every
  ! rank computes this reference for itself on the undivided mesh): owner samples, viscous
  ! records, and the deposited fields gathered over the ranks element by element and dof by
  ! dof.
  !
  ! Cases in the configuration
  !   * a pair of ranks that are neighbours directly and across a periodic plane (the two
  !     blocks in y; on four ranks also in x)
  !   * a particle at a periodic corner, one exactly on a periodic plane, one whose position
  !     is outside the cell (not wrapped), one next to a wall
  !   * workers that own no particle
  !   * on eight ranks the kernel radius equals the gap between next-nearest blocks in x up
  !     to rounding (0.005 = L/4 written as 2 x 2.5 x 0.001), the tie on which two ranks must
  !     reach the same verdict: the neighbour lists of all ranks are gathered and must be
  !     symmetric
  !
  ! What it does not cover: the dof pairing of the fluid across a periodic plane (the fields
  ! are compared before any parallel sum of the solver), EL_PARTICLE_MESH_PASS itself (needs
  ! the particle engine and the solver arrays) and the force closures.

  USE MPI
  USE TYPES, ONLY: tMultiMesh
  USE PP3D_MPI, ONLY: myid, master, showid, MPI_COMM_SUBS
  USE DEM_QUERY, ONLY: tParticleData
  USE EL_CONFIG, ONLY: el_kernel, el_kernel_width_factor, EL_KERNEL_SUPPORT_RADIUS
  USE EL_FIELDS, ONLY: el_field_data, EL_ENSURE_FIELDS, EL_BEGIN_FIELD_UPDATE, EL_FINALIZE
  USE EL_HALO, ONLY: tElParticleRecord, halo, EL_REFRESH_COUPLING_HALO, EL_REDUCE_TO_OWNERS, &
                     EL_BROADCAST_FROM_OWNERS, EL_HALO_LOCAL_BOX, EL_HALO_DOMAIN_BOX, &
                     EL_FINALIZE_HALO
  USE EL_PERIODIC, ONLY: EL_SET_PERIODIC_LATTICE, EL_CLEAR_PERIODIC_LATTICE, &
                         EL_PERIODIC_MARKED_ON_PLANES
  USE EL_IMAGE_TRANSFER, ONLY: EL_SAMPLE_RECORD, EL_SAMPLE_RECORD_VISCOUS, EL_DEPOSIT_RECORD
  USE EL_QUADRATURE, ONLY: EL_SAMPLE_SIZE, EL_SAMPLE_NORMALIZATION
  USE EL_VISCOUS, ONLY: EL_VIS_SIZE
  USE EL_TEST_HELPERS, ONLY: EL_TEST_INIT, EL_TEST_CHECK, EL_TEST_REPORT, &
                             EL_TEST_RELEASE_MESH, BUILD_DOF_COORDINATES, &
                             nel, nvt, net, nat, nve, nee, nae, nvel, neel, nved, nvar, &
                             near, nbct, nvbd, nebd, nabd

  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 8                        ! cells per axis of the whole cell
  INTEGER, PARAMETER :: npart = 6
  REAL*8, PARAMETER :: period = 0.02d0
  REAL*8, PARAMETER :: particle_radius = 0.001d0     ! kernel radius 0.005 = L/4 = 2 cells
  REAL*8, PARAMETER :: tol = 1.0d-11

  TYPE(tMultiMesh) :: mesh
  REAL*8, ALLOCATABLE :: coor(:,:), val_u(:), val_v(:), val_w(:), val_p(:), eps(:)
  LOGICAL, ALLOCATABLE :: wall(:)
  INTEGER, ALLOCATABLE :: global_dof(:,:,:)          ! half-step lattice index -> dof of the cell
  INTEGER :: ndof, rank, nproc, ierr, failures, mfile, px, py
  INTEGER :: bx, by, nc(3), i0(3)                    ! block of this rank
  INTEGER :: mx, my, mz, ex, ey, fx, fy              ! cells and index offsets of the current mesh block
  REAL*8 :: dh, pi, positions(3,npart), forces(3,npart), drag_b(npart)
  LOGICAL :: lattice_ok

  ! the reference: all particles coupled on one rank holding the whole cell
  REAL*8 :: ref_sample(EL_SAMPLE_SIZE,npart), ref_vis(EL_VIS_SIZE,npart)
  REAL*8, ALLOCATABLE :: ref_alpha(:), ref_grad(:,:), ref_force(:,:), ref_drag(:)
  INTEGER :: nel_cell, ndof_cell

  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  CALL EL_TEST_INIT()
  failures = 0

  IF (nproc.NE.4 .AND. nproc.NE.8) THEN
    IF (rank.EQ.0) WRITE(*,'(A,I0,A)') 'EL periodic kernel (distributed): needs 4 or 8 ranks, got ', &
      nproc, '; nothing tested.'
    CALL MPI_Finalize(ierr)
    STOP
  END IF

  ! every rank is a worker of the coupling; the roles of the application's communicator
  MPI_COMM_SUBS = MPI_COMM_WORLD
  master = -1
  myid = rank + 1
  showid = 1
  OPEN(NEWUNIT=mfile, STATUS='SCRATCH', ACTION='READWRITE')

  pi = ACOS(-1.0d0)
  dh = period/DBLE(n)
  el_kernel = 'deen_poly'
  el_kernel_width_factor = 2.5d0
  px = nproc/2
  py = 2
  bx = MOD(rank, px)
  by = rank/px

  CALL set_particles()
  CALL EL_SET_PERIODIC_LATTICE((/ .TRUE., .TRUE., .FALSE. /), (/ period, period, period /), &
                               lattice_ok)

  CALL compute_reference()
  CALL test_distributed(failures)

  CALL EL_FINALIZE_HALO()
  CALL EL_CLEAR_PERIODIC_LATTICE()
  CALL EL_TEST_REPORT(failures, rank, nproc, 'EL periodic kernel (distributed)')
  CALL EL_TEST_RELEASE_MESH(mesh)
  CALL EL_FINALIZE()
  CLOSE(mfile)
  CALL MPI_Finalize(ierr)

  IF (failures.GT.0) ERROR STOP 1

CONTAINS

  !----------------------------------------------------------------------
  ! Particles (positions in periods), forces and drag coefficients
  !----------------------------------------------------------------------
  SUBROUTINE set_particles()
    positions(:,1) = (/ 0.03d0, 0.04d0, 0.50d0 /)     ! periodic corner
    positions(:,2) = (/ 0.97d0, 0.49d0, 0.30d0 /)     ! x plane and the block boundary in y
    positions(:,3) = (/ 0.30d0, 0.30d0, 0.60d0 /)     ! interior
    positions(:,4) = (/ 1.02d0, 0.98d0, 0.50d0 /)     ! outside the cell, not wrapped
    positions(:,5) = (/ 0.51d0, 0.01d0, 0.10d0 /)     ! wall z = 0 and the y plane
    positions(:,6) = (/ 0.00d0, 0.75d0, 0.45d0 /)     ! exactly on the x plane
    positions = positions*period
    forces(:,1) = (/ 0.3d0, -1.1d0, 0.7d0 /)
    forces(:,2) = (/ -0.4d0, 0.2d0, 0.9d0 /)
    forces(:,3) = (/ 1.0d0, 0.5d0, -0.3d0 /)
    forces(:,4) = (/ -0.8d0, -0.6d0, 0.1d0 /)
    forces(:,5) = (/ 0.2d0, 0.7d0, -1.2d0 /)
    forces(:,6) = (/ 0.6d0, -0.9d0, 0.4d0 /)
    drag_b = (/ 2.5d0, 1.5d0, 0.5d0, 3.0d0, 1.0d0, 2.0d0 /)
  END SUBROUTINE set_particles

  FUNCTION record_of(ip) RESULT(p)
    INTEGER, INTENT(IN) :: ip
    TYPE(tElParticleRecord) :: p
    p%position = positions(:,ip)
    p%velocity = (/ 0.1d0, -0.2d0, 0.3d0 /)
    p%radius = particle_radius
    p%density = 2.0d0
    p%system_id = ip
  END FUNCTION record_of

  ! Block of the partition that owns a particle: the block that holds the position wrapped
  ! into the cell (the particle engine owns by the wrapped position as well).
  LOGICAL FUNCTION owned_here(ip)
    INTEGER, INTENT(IN) :: ip
    REAL*8 :: w(2)
    INTEGER :: ib(2)
    w = positions(1:2,ip)/period
    w = w - FLOOR(w)
    ib(1) = MIN(px-1, INT(w(1)*DBLE(px)))
    ib(2) = MIN(py-1, INT(w(2)*DBLE(py)))
    owned_here = (ib(1).EQ.bx .AND. ib(2).EQ.by)
  END FUNCTION owned_here

  !----------------------------------------------------------------------
  ! The reference on the undivided cell
  !----------------------------------------------------------------------
  SUBROUTINE compute_reference()
    TYPE(tElParticleRecord) :: p
    REAL*8 :: cell_box(6)
    INTEGER :: ip, i, ix(3)

    CALL build_block(mesh, (/ n, n, n /), (/ 0, 0, 0 /))
    CALL set_block_fields()
    nel_cell = nel
    ndof_cell = ndof
    cell_box = (/ 0.0d0, period, 0.0d0, period, 0.0d0, period /)

    ALLOCATE(global_dof(0:2*n,0:2*n,0:2*n))
    global_dof = 0
    DO i = 1, ndof
      ix = NINT(2.0d0*coor(:,i)/dh)
      global_dof(ix(1),ix(2),ix(3)) = i
    END DO

    CALL EL_ENSURE_FIELDS(nel, ndof)
    CALL EL_BEGIN_FIELD_UPDATE(.TRUE.)
    DO ip = 1, npart
      p = record_of(ip)
      CALL EL_SAMPLE_RECORD(mesh, 1, val_u, val_v, val_w, val_p, eps, p, cell_box, &
                            ref_sample(:,ip))
      CALL EL_SAMPLE_RECORD_VISCOUS(mesh, 1, coor, val_u, val_v, val_w, wall, p, cell_box, &
                                    ref_vis(:,ip))
    END DO
    DO ip = 1, npart
      CALL EL_DEPOSIT_RECORD(mesh, 1, record_of(ip), cell_box, &
        ref_sample(EL_SAMPLE_NORMALIZATION,ip), forces(:,ip), el_field_data, &
        drag_b=drag_b(ip))
    END DO
    ref_alpha = el_field_data%alpha_p
    ref_grad = el_field_data%grad_alpha_p
    ref_force = el_field_data%force_rhs
    ref_drag = el_field_data%drag_B_ml
  END SUBROUTINE compute_reference

  !----------------------------------------------------------------------
  ! The distributed pass
  !----------------------------------------------------------------------
  SUBROUTINE test_distributed(nf)
    INTEGER, INTENT(INOUT) :: nf
    TYPE(tParticleData), ALLOCATABLE :: owned(:)
    TYPE(tElParticleRecord), ALLOCATABLE :: records(:)
    REAL*8, ALLOCATABLE :: local_samples(:,:), owned_samples(:,:), local_vis(:,:), owned_vis(:,:)
    REAL*8, ALLOCATABLE :: owned_norm(:,:), record_norm(:,:)
    REAL*8, ALLOCATABLE :: part(:), total(:)
    INTEGER, ALLOCATABLE :: adj_local(:), adj(:)
    REAL*8 :: local_box(6), domain_box(6), err, err_all, scale, delta
    INTEGER :: ip, i, j, k, n_owned, n_records, count_local, count_all, min_owned, max_owned
    INTEGER :: iel, ix(3), g, asym, deg_min, deg_max, max_records, n_marked, n_marked_all, c
    INTEGER :: n_both, n_both_all

    nc = (/ n/px, n/py, n /)
    i0 = (/ bx*nc(1), by*nc(2), 0 /)
    CALL build_block(mesh, nc, i0)
    CALL set_block_fields()
    CALL EL_ENSURE_FIELDS(nel, ndof)
    CALL EL_BEGIN_FIELD_UPDATE(.TRUE.)

    n_owned = 0
    DO ip = 1, npart
      IF (owned_here(ip)) n_owned = n_owned + 1
    END DO
    ALLOCATE(owned(n_owned))
    k = 0
    DO ip = 1, npart
      IF (.NOT.owned_here(ip)) CYCLE
      k = k + 1
      owned(k)%position = positions(:,ip)
      owned(k)%velocity = (/ 0.1d0, -0.2d0, 0.3d0 /)
      owned(k)%angvel = 0.0d0
      owned(k)%time = 0.0d0
      owned(k)%density = 2.0d0
      owned(k)%aabb = 0.0d0
      owned(k)%radius = particle_radius
      owned(k)%typeId = 0
      owned(k)%localIdx = k - 1
      owned(k)%uniqueIdx = ip - 1
      owned(k)%systemIdx = ip
      owned(k)%bytes = 0
      owned(k)%bytes(1) = INT(ip, KIND(owned(k)%bytes))
    END DO
    CALL MPI_Allreduce(n_owned, count_all, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(n_owned, min_owned, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(n_owned, max_owned, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL EL_TEST_CHECK(count_all.EQ.npart .AND. min_owned.EQ.0, &
      'setup: every particle has one owner, at least one worker owns none', nf, rank)

    ! --- neighbour graph and record exchange of the application
    CALL EL_REFRESH_COUPLING_HALO(mesh, 1, owned, records, mfile)
    n_records = SIZE(records)
    local_box = EL_HALO_LOCAL_BOX()
    domain_box = EL_HALO_DOMAIN_BOX()
    delta = EL_KERNEL_SUPPORT_RADIUS(particle_radius)

    ALLOCATE(adj_local(nproc*nproc), adj(nproc*nproc))
    adj_local = 0
    DO j = 1, halo%neighbor_count
      adj_local(1 + rank + nproc*halo%neighbors(j)) = 1
    END DO
    CALL MPI_Allreduce(adj_local, adj, nproc*nproc, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    asym = 0
    deg_min = nproc
    deg_max = 0
    DO i = 0, nproc-1
      k = 0
      DO j = 0, nproc-1
        IF (adj(1+i+nproc*j).NE.adj(1+j+nproc*i)) asym = asym + 1
        k = k + adj(1+i+nproc*j)
      END DO
      deg_min = MIN(deg_min, k)
      deg_max = MAX(deg_max, k)
    END DO
    IF (rank.EQ.0) WRITE(*,'(A,I0,A,I0,A,I0,A,ES22.15,A,ES22.15)') '    neighbour ranks per worker: ', &
      deg_min, ' to ', deg_max, ' of ', nproc-1, '; kernel radius ', delta, ', block width in x ', &
      local_box(2)-local_box(1)
    CALL EL_TEST_CHECK(asym.EQ.0, &
      'graph: the neighbour lists of all workers are symmetric (also at the tie distance)', nf, rank)
    ! every block touches the blocks next to it in x and the other row in y: 3 on four ranks,
    ! at least 5 on eight (the next-nearest blocks in x are the tie)
    CALL EL_TEST_CHECK(deg_min.GE.MERGE(3, 5, nproc.EQ.4), &
      'graph: direct and wrap neighbours are all linked', nf, rank)
    CALL EL_TEST_CHECK(MAXVAL(ABS(domain_box - (/ 0.0d0, period, 0.0d0, period, 0.0d0, period /))) &
      .LE.1.0d-15, 'graph: the gathered boxes span the lattice cell', nf, rank)
    ! the pair of blocks in y is linked directly and across the plane: both gaps are zero
    n_both = 0
    DO i = 1, nproc
      IF (i-1.EQ.rank) CYCLE
      IF (direct_gap(local_box, halo%boxes(:,i)).LE.0.0d0 .AND. &
          wrap_only_gap(local_box, halo%boxes(:,i)).LE.0.0d0) n_both = n_both + 1
    END DO
    CALL MPI_Allreduce(n_both, n_both_all, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    CALL EL_TEST_CHECK(n_both_all.GE.1, &
      'graph: every worker has a neighbour that is both a direct and a wrap neighbour', nf, rank)

    count_local = n_records
    CALL MPI_Allreduce(count_local, count_all, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(count_local, max_records, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    count_local = MERGE(1, 0, n_owned.EQ.0 .AND. n_records.GT.0)
    CALL MPI_Allreduce(count_local, c, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    IF (rank.EQ.0) WRITE(*,'(A,I0,A,I0,A,I0)') '    records over all workers: ', count_all, &
      ', most on one worker ', max_records, ', workers with records but no particle of their own: ', c
    CALL EL_TEST_CHECK(count_all.GT.npart .AND. max_records.LE.npart .AND. c.GE.1, &
      'exchange: records reach other workers, one per particle and worker, also particle-free ones', &
      nf, rank)

    ! --- sample, viscous pieces, sum at the owner
    ALLOCATE(local_samples(EL_SAMPLE_SIZE,n_records), owned_samples(EL_SAMPLE_SIZE,n_owned))
    ALLOCATE(local_vis(EL_VIS_SIZE,n_records), owned_vis(EL_VIS_SIZE,n_owned))
    n_marked = EL_PERIODIC_MARKED_ON_PLANES(coor, wall, domain_box)
    CALL MPI_Allreduce(n_marked, n_marked_all, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(n_marked, c, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    CALL EL_TEST_CHECK(n_marked_all.EQ.0 .AND. c.EQ.0, &
      'viscous: the wall mask of every block holds no face on a periodic plane', nf, rank)
    DO i = 1, n_records
      CALL EL_SAMPLE_RECORD(mesh, 1, val_u, val_v, val_w, val_p, eps, records(i), local_box, &
                            local_samples(:,i))
      CALL EL_SAMPLE_RECORD_VISCOUS(mesh, 1, coor, val_u, val_v, val_w, wall, records(i), &
                                    local_box, local_vis(:,i))
    END DO
    owned_samples = 0.0d0
    owned_vis = 0.0d0
    CALL EL_REDUCE_TO_OWNERS(local_samples, owned_samples)
    CALL EL_REDUCE_TO_OWNERS(local_vis, owned_vis)

    err = 0.0d0
    DO i = 1, n_owned
      ip = records(i)%system_id
      scale = MAXVAL(ABS(ref_sample(:,ip)))
      err = MAX(err, MAXVAL(ABS(owned_samples(:,i) - ref_sample(:,ip)))/scale)
    END DO
    CALL MPI_Allreduce(err, err_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    IF (rank.EQ.0) WRITE(*,'(A,ES10.2)') '    owner samples against one rank, largest rel. difference ', err_all
    CALL EL_TEST_CHECK(err_all.LE.tol, &
      'sample: the sum at the owner equals the sample on one rank (all slots, all particles)', &
      nf, rank)
    err = 0.0d0
    DO i = 1, n_owned
      ip = records(i)%system_id
      scale = MAXVAL(ABS(ref_vis(:,ip)))
      err = MAX(err, MAXVAL(ABS(owned_vis(:,i) - ref_vis(:,ip)))/scale)
    END DO
    CALL MPI_Allreduce(err, err_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    IF (rank.EQ.0) WRITE(*,'(A,ES10.2)') '    owner viscous records against one rank, largest rel. difference ', err_all
    CALL EL_TEST_CHECK(err_all.LE.tol, &
      'viscous: the sum at the owner equals the record on one rank', nf, rank)

    ! --- normalisation back to the record holders, deposit
    ALLOCATE(owned_norm(1,n_owned), record_norm(1,n_records))
    owned_norm(1,:) = owned_samples(EL_SAMPLE_NORMALIZATION,:)
    record_norm = 0.0d0
    CALL EL_BROADCAST_FROM_OWNERS(owned_norm, record_norm)
    err = 0.0d0
    DO i = 1, n_records
      ip = records(i)%system_id
      err = MAX(err, ABS(record_norm(1,i)/ref_sample(EL_SAMPLE_NORMALIZATION,ip) - 1.0d0))
    END DO
    CALL MPI_Allreduce(err, err_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL EL_TEST_CHECK(err_all.LE.tol, &
      'normalisation: every holder of a record receives the kernel integral over all images', &
      nf, rank)
    DO i = 1, n_records
      ip = records(i)%system_id
      CALL EL_DEPOSIT_RECORD(mesh, 1, records(i), local_box, record_norm(1,i), forces(:,ip), &
                             el_field_data, drag_b=drag_b(ip))
    END DO

    ! --- deposited fields, gathered on the numbering of the whole cell
    ALLOCATE(part(4*nel_cell), total(4*nel_cell))
    part = 0.0d0
    DO iel = 1, nel
      ix = NINT((coor(:, nvt+net+nat+iel)/dh - 0.5d0))
      g = 1 + ix(1) + n*(ix(2) + n*ix(3))
      part(g) = el_field_data%alpha_p(iel)
      part(nel_cell+3*(g-1)+1:nel_cell+3*g) = el_field_data%grad_alpha_p(:,iel)
    END DO
    CALL MPI_Allreduce(part, total, 4*nel_cell, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    err = MAXVAL(ABS(total(1:nel_cell) - ref_alpha))/MAXVAL(ABS(ref_alpha))
    CALL EL_TEST_CHECK(err.LE.tol, 'deposit: solid volume fraction of every element as on one rank', &
      nf, rank)
    err = MAXVAL(ABS(total(nel_cell+1:) - RESHAPE(ref_grad, (/ 3*nel_cell /))))/MAXVAL(ABS(ref_grad))
    CALL EL_TEST_CHECK(err.LE.tol, 'deposit: gradient of the solid volume fraction as on one rank', &
      nf, rank)
    DEALLOCATE(part, total)

    ALLOCATE(part(4*ndof_cell), total(4*ndof_cell))
    part = 0.0d0
    DO i = 1, ndof
      ix = NINT(2.0d0*coor(:,i)/dh)
      g = global_dof(ix(1),ix(2),ix(3))
      part(3*(g-1)+1:3*g) = part(3*(g-1)+1:3*g) + el_field_data%force_rhs(:,i)
      part(3*ndof_cell+g) = part(3*ndof_cell+g) + el_field_data%drag_B_ml(i)
    END DO
    CALL MPI_Allreduce(part, total, 4*ndof_cell, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    err = MAXVAL(ABS(total(1:3*ndof_cell) - RESHAPE(ref_force, (/ 3*ndof_cell /))))/MAXVAL(ABS(ref_force))
    IF (rank.EQ.0) WRITE(*,'(A,ES10.2)') '    reaction force per dof against one rank, largest rel. difference ', err
    CALL EL_TEST_CHECK(err.LE.tol, 'deposit: reaction force on every dof as on one rank', nf, rank)
    err = MAXVAL(ABS(total(3*ndof_cell+1:) - ref_drag))/MAXVAL(ABS(ref_drag))
    CALL EL_TEST_CHECK(err.LE.tol, 'deposit: drag coefficient on every dof as on one rank', nf, rank)
    err = 0.0d0
    DO c = 1, 3
      err = MAX(err, ABS(SUM(total(c:3*ndof_cell:3)) + SUM(forces(c,:))))
    END DO
    CALL EL_TEST_CHECK(err.LE.tol*MAXVAL(ABS(forces)), &
      'deposit: the reaction forces of all particles sum to minus the forces', nf, rank)
    err = ABS(SUM(total(3*ndof_cell+1:)) - SUM(drag_b))/SUM(drag_b)
    CALL EL_TEST_CHECK(err.LE.tol, 'deposit: the drag coefficients are conserved', nf, rank)
    DEALLOCATE(part, total)
  END SUBROUTINE test_distributed

  ! gap between two boxes without the lattice (zero when they touch)
  REAL*8 FUNCTION direct_gap(a, b)
    REAL*8, INTENT(IN) :: a(6), b(6)
    REAL*8 :: g(3)
    INTEGER :: k
    DO k = 1, 3
      g(k) = MAX(0.0d0, a(2*k-1)-b(2*k), b(2*k-1)-a(2*k))
    END DO
    direct_gap = SQRT(SUM(g*g))
  END FUNCTION direct_gap

  ! gap between box a and the images of box b shifted by one period in x or y (not b itself)
  REAL*8 FUNCTION wrap_only_gap(a, b)
    REAL*8, INTENT(IN) :: a(6), b(6)
    REAL*8 :: shifted(6)
    INTEGER :: sx, sy
    wrap_only_gap = HUGE(1.0d0)
    DO sy = -1, 1
      DO sx = -1, 1
        IF (sx.EQ.0 .AND. sy.EQ.0) CYCLE
        shifted = b
        shifted(1:2) = shifted(1:2) + DBLE(sx)*period
        shifted(3:4) = shifted(3:4) + DBLE(sy)*period
        wrap_only_gap = MIN(wrap_only_gap, direct_gap(a, shifted))
      END DO
    END DO
  END FUNCTION wrap_only_gap

  !----------------------------------------------------------------------
  ! Mesh block and fields
  !----------------------------------------------------------------------

  ! Structured Q2/P1 block of nc(1) x nc(2) x nc(3) cells of size dh whose lowest vertex has
  ! the lattice index i0 (the layout of EL_TEST_HELPERS: BUILD_STRUCTURED_MESH, for a block
  ! that need not be a cube and need not start at the origin). Sets /TRIAD/, the dof
  ! coordinates, and the mask of the wall faces (z planes of the cell).
  SUBROUTINE build_block(m, cells, first)
    TYPE(tMultiMesh), INTENT(INOUT) :: m
    INTEGER, INTENT(IN) :: cells(3), first(3)
    INTEGER :: ie, je, ke, iel, i

    CALL EL_TEST_RELEASE_MESH(m)
    mx = cells(1)
    my = cells(2)
    mz = cells(3)
    nel = mx*my*mz
    nvt = (mx+1)*(my+1)*(mz+1)
    ex = mx*(my+1)*(mz+1)
    ey = (mx+1)*my*(mz+1)
    net = ex + ey + (mx+1)*(my+1)*mz
    fx = (mx+1)*my*mz
    fy = mx*(my+1)*mz
    nat = fx + fy + mx*my*(mz+1)
    nve = 8
    nee = 12
    nae = 6
    nvel = 0
    neel = 0
    nved = 0
    nvar = 0
    near = 0
    nbct = 1
    nvbd = 0
    nebd = 0
    nabd = 0

    ALLOCATE(m%level(1))
    m%nlmin = 1
    m%nlmax = 1
    m%maxlevel = 1
    m%level(1)%nel = nel
    m%level(1)%nvt = nvt
    m%level(1)%net = net
    m%level(1)%nat = nat
    ALLOCATE(m%level(1)%dcorvg(3,nvt), m%level(1)%kvert(8,nel), m%level(1)%kedge(12,nel))
    ALLOCATE(m%level(1)%karea(6,nel), m%level(1)%dvol(nel))

    DO ke = 0, mz
      DO je = 0, my
        DO ie = 0, mx
          m%level(1)%dcorvg(:, v(ie,je,ke)) = (/ DBLE(first(1)+ie)*dh, DBLE(first(2)+je)*dh, &
                                                 DBLE(first(3)+ke)*dh /)
        END DO
      END DO
    END DO
    DO ke = 0, mz-1
      DO je = 0, my-1
        DO ie = 0, mx-1
          iel = 1 + ie + mx*(je + my*ke)
          m%level(1)%dvol(iel) = dh*dh*dh
          m%level(1)%kvert(:,iel) = (/ v(ie,je,ke), v(ie+1,je,ke), v(ie+1,je+1,ke), v(ie,je+1,ke), &
            v(ie,je,ke+1), v(ie+1,je,ke+1), v(ie+1,je+1,ke+1), v(ie,je+1,ke+1) /)
          m%level(1)%kedge(:,iel) = (/ e1(ie,je,ke), e2(ie+1,je,ke), e1(ie,je+1,ke), e2(ie,je,ke), &
            e3(ie,je,ke), e3(ie+1,je,ke), e3(ie+1,je+1,ke), e3(ie,je+1,ke), &
            e1(ie,je,ke+1), e2(ie+1,je,ke+1), e1(ie,je+1,ke+1), e2(ie,je,ke+1) /)
          m%level(1)%karea(:,iel) = (/ f3(ie,je,ke), f2(ie,je,ke), f1(ie+1,je,ke), f2(ie,je+1,ke), &
            f1(ie,je,ke), f3(ie,je,ke+1) /)
        END DO
      END DO
    END DO

    ndof = nvt + net + nat + nel
    IF (ALLOCATED(coor)) DEALLOCATE(coor, val_u, val_v, val_w, val_p, eps, wall)
    ALLOCATE(coor(3,ndof), val_u(ndof), val_v(ndof), val_w(ndof), val_p(4*nel), eps(nel), wall(ndof))
    CALL BUILD_DOF_COORDINATES(m, coor)
    ! boundary-face mask of the particle viscous force: face-centre dofs on the walls
    wall = .FALSE.
    DO i = nvt+net+1, nvt+net+nat
      wall(i) = (coor(3,i).LT.1.0d-12 .OR. coor(3,i).GT.period-1.0d-12)
    END DO

  END SUBROUTINE build_block

  ! Vertex, edge and face numbers of the current block (mx, my, mz cells; build_block)
    INTEGER FUNCTION v(i, j, k)
      INTEGER, INTENT(IN) :: i, j, k
      v = 1 + i + (mx+1)*(j + (my+1)*k)
    END FUNCTION v
    INTEGER FUNCTION e1(i, j, k)          ! edges along x
      INTEGER, INTENT(IN) :: i, j, k
      e1 = 1 + i + mx*(j + (my+1)*k)
    END FUNCTION e1
    INTEGER FUNCTION e2(i, j, k)          ! edges along y
      INTEGER, INTENT(IN) :: i, j, k
      e2 = ex + 1 + i + (mx+1)*(j + my*k)
    END FUNCTION e2
    INTEGER FUNCTION e3(i, j, k)          ! edges along z
      INTEGER, INTENT(IN) :: i, j, k
      e3 = ex + ey + 1 + i + (mx+1)*(j + (my+1)*k)
    END FUNCTION e3
    INTEGER FUNCTION f1(i, j, k)          ! faces normal to x
      INTEGER, INTENT(IN) :: i, j, k
      f1 = 1 + i + (mx+1)*(j + my*k)
    END FUNCTION f1
    INTEGER FUNCTION f2(i, j, k)          ! faces normal to y
      INTEGER, INTENT(IN) :: i, j, k
      f2 = fx + 1 + i + mx*(j + (my+1)*k)
    END FUNCTION f2
    INTEGER FUNCTION f3(i, j, k)          ! faces normal to z
      INTEGER, INTENT(IN) :: i, j, k
      f3 = fx + fy + 1 + i + mx*(j + my*k)
    END FUNCTION f3

  ! Smooth fields of period L in x and y on the dofs and cells of the current block: the same
  ! function of the coordinates on every block and on the whole cell.
  SUBROUTINE set_block_fields()
    REAL*8 :: x(3), tp
    INTEGER :: i, iel, jj

    tp = 2.0d0*pi/period
    DO i = 1, ndof
      x = coor(:,i)
      val_u(i) = SIN(tp*x(1))*COS(tp*x(2)) + 0.5d0*SIN(tp*x(3))
      val_v(i) = COS(tp*x(1)) + SIN(tp*x(2))*SIN(tp*x(3))
      val_w(i) = 1.0d0 + 0.3d0*COS(tp*(x(1)+x(2))) + 0.2d0*SIN(tp*x(3))
    END DO
    DO iel = 1, nel
      x = coor(:, nvt+net+nat+iel)
      jj = 4*(iel-1) + 1
      val_p(jj) = COS(tp*x(1))*SIN(tp*x(2)) + 0.1d0*COS(tp*x(3))
      val_p(jj+1) = -tp*SIN(tp*x(1))*SIN(tp*x(2))
      val_p(jj+2) = tp*COS(tp*x(1))*COS(tp*x(2))
      val_p(jj+3) = -0.1d0*tp*SIN(tp*x(3))
      eps(iel) = 0.8d0 + 0.1d0*SIN(tp*x(1))*COS(tp*x(3)) + 0.05d0*COS(tp*x(2))
    END DO
  END SUBROUTINE set_block_fields

END PROGRAM TEST_EL_PERIODIC_MPI
