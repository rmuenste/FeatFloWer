MODULE EL_HALO

  USE ISO_C_BINDING, ONLY: c_int
  USE MPI
  USE TYPES, ONLY: tMultiMesh
  USE PP3D_MPI, ONLY: myid, master, showid, MPI_COMM_SUBS
  USE DEM_QUERY, ONLY: tParticleData
  USE EL_CONFIG, ONLY: el_kernel_width_factor

  IMPLICIT NONE

  INTEGER, PARAMETER :: EL_RECORD_REAL_COUNT = 8
  INTEGER, PARAMETER :: EL_RECORD_INT_COUNT = 3

  TYPE tElParticleRecord
    REAL*8 :: position(3) = 0.0d0
    REAL*8 :: velocity(3) = 0.0d0
    REAL*8 :: radius = 0.0d0
    REAL*8 :: density = 0.0d0
    INTEGER(c_int) :: system_id = 0
    INTEGER :: owner_rank = -1
    INTEGER :: owner_slot = 0
    LOGICAL :: owned = .FALSE.
  END TYPE tElParticleRecord

  TYPE tElHaloState
    LOGICAL :: initialized = .FALSE.
    LOGICAL :: warned_wide_kernel = .FALSE.
    INTEGER :: communicator = MPI_COMM_NULL
    INTEGER :: rank = -1
    INTEGER :: size = 0
    INTEGER :: neighbor_count = 0
    INTEGER, ALLOCATABLE :: neighbors(:)
    REAL*8, ALLOCATABLE :: boxes(:,:)
    REAL*8 :: delta_max = -1.0d0
    INTEGER, ALLOCATABLE :: send_counts(:), recv_counts(:)
    INTEGER, ALLOCATABLE :: send_displs(:), recv_displs(:)
    INTEGER, ALLOCATABLE :: send_owner_slots(:)
  END TYPE tElHaloState

  TYPE(tElHaloState), SAVE :: halo

CONTAINS

  SUBROUTINE EL_REFRESH_COUPLING_HALO(mesh, ilev, owned_particles, records, mfile)

    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev, mfile
    TYPE(tParticleData), INTENT(IN) :: owned_particles(:)
    TYPE(tElParticleRecord), ALLOCATABLE, INTENT(OUT) :: records(:)

    REAL*8 :: local_radius, global_radius, requested_delta
    INTEGER :: ierr

    IF (myid.EQ.master) THEN
      ALLOCATE(records(0))
      RETURN
    END IF

    local_radius = 0.0d0
    IF (SIZE(owned_particles).GT.0) local_radius = MAXVAL(owned_particles%radius)
    CALL MPI_Allreduce(local_radius, global_radius, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, MPI_COMM_SUBS, ierr)
    requested_delta = 2.0d0*el_kernel_width_factor*global_radius

    IF (.NOT.halo%initialized .OR. &
        ABS(requested_delta-halo%delta_max).GT.1.0d-12*MAX(1.0d0,requested_delta)) THEN
      CALL EL_BUILD_HALO_GRAPH(mesh, ilev, requested_delta, mfile)
    END IF

    CALL EL_EXCHANGE_PARTICLE_RECORDS(owned_particles, records)

  END SUBROUTINE EL_REFRESH_COUPLING_HALO

  SUBROUTINE EL_BUILD_HALO_GRAPH(mesh, ilev, delta_max, mfile)

    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev, mfile
    REAL*8, INTENT(IN) :: delta_max

    REAL*8 :: local_box(6), min_width, global_min_width
    INTEGER :: i, ierr, degree, affected_local, affected_global
    INTEGER :: max_layers_local, max_layers_global
    INTEGER, ALLOCATABLE :: candidate_neighbors(:)
    LOGICAL :: overlap

    CALL EL_RELEASE_GRAPH()

    CALL MPI_Comm_rank(MPI_COMM_SUBS, halo%rank, ierr)
    CALL MPI_Comm_size(MPI_COMM_SUBS, halo%size, ierr)

    local_box = (/ MINVAL(mesh%level(ilev)%dcorvg(1,:)), &
                   MAXVAL(mesh%level(ilev)%dcorvg(1,:)), &
                   MINVAL(mesh%level(ilev)%dcorvg(2,:)), &
                   MAXVAL(mesh%level(ilev)%dcorvg(2,:)), &
                   MINVAL(mesh%level(ilev)%dcorvg(3,:)), &
                   MAXVAL(mesh%level(ilev)%dcorvg(3,:)) /)
    ALLOCATE(halo%boxes(6,halo%size))
    CALL MPI_Allgather(local_box, 6, MPI_DOUBLE_PRECISION, halo%boxes, 6, &
                       MPI_DOUBLE_PRECISION, MPI_COMM_SUBS, ierr)

    ALLOCATE(candidate_neighbors(MAX(1,halo%size-1)))
    degree = 0
    DO i=1,halo%size
      IF (i-1.EQ.halo%rank) CYCLE
      overlap = EL_BOXES_WITHIN_DELTA(local_box, halo%boxes(:,i), delta_max)
      IF (overlap) THEN
        degree = degree + 1
        candidate_neighbors(degree) = i-1
      END IF
    END DO

    halo%neighbor_count = degree
    ALLOCATE(halo%neighbors(MAX(1,degree)))
    IF (degree.GT.0) halo%neighbors = candidate_neighbors(1:degree)
    DEALLOCATE(candidate_neighbors)

    CALL MPI_Dist_graph_create_adjacent(MPI_COMM_SUBS, degree, halo%neighbors, &
      MPI_UNWEIGHTED, degree, halo%neighbors, MPI_UNWEIGHTED, MPI_INFO_NULL, &
      .FALSE., halo%communicator, ierr)

    ALLOCATE(halo%send_counts(MAX(1,degree)), halo%recv_counts(MAX(1,degree)))
    ALLOCATE(halo%send_displs(MAX(1,degree)), halo%recv_displs(MAX(1,degree)))
    halo%send_counts = 0
    halo%recv_counts = 0
    halo%send_displs = 0
    halo%recv_displs = 0
    halo%delta_max = delta_max
    halo%initialized = .TRUE.

    min_width = MIN(local_box(2)-local_box(1), local_box(4)-local_box(3), &
                    local_box(6)-local_box(5))
    CALL MPI_Allreduce(min_width, global_min_width, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MIN, MPI_COMM_SUBS, ierr)
    affected_local = MERGE(1, 0, delta_max.GE.min_width)
    max_layers_local = CEILING(delta_max/MAX(min_width,TINY(1.0d0)))
    CALL MPI_Allreduce(affected_local, affected_global, 1, MPI_INTEGER, &
                       MPI_SUM, MPI_COMM_SUBS, ierr)
    CALL MPI_Allreduce(max_layers_local, max_layers_global, 1, MPI_INTEGER, &
                       MPI_MAX, MPI_COMM_SUBS, ierr)

    IF (affected_global.GT.0 .AND. .NOT.halo%warned_wide_kernel) THEN
      IF (myid.EQ.showid) THEN
        WRITE(*,'(A)') 'WARNING: E-L kernel support spans at least one CFD partition width.'
        WRITE(*,'(A,ES14.6)') '  maximum kernel radius = ', delta_max
        WRITE(*,'(A,ES14.6)') '  minimum partition width = ', global_min_width
        WRITE(*,'(A,I0)') '  affected CFD ranks = ', affected_global
        WRITE(*,'(A,I0)') '  maximum partition layers crossed = ', max_layers_global
        WRITE(*,'(A)') '  Transfer remains conservative, but this decomposition is inefficient.'
        WRITE(mfile,'(A)') 'WARNING: E-L kernel support spans at least one CFD partition width.'
        WRITE(mfile,'(A,2ES14.6,2I8)') '  delta/min-width/ranks/layers: ', &
          delta_max, global_min_width, affected_global, max_layers_global
      END IF
      halo%warned_wide_kernel = .TRUE.
    END IF

  END SUBROUTINE EL_BUILD_HALO_GRAPH

  SUBROUTINE EL_EXCHANGE_PARTICLE_RECORDS(owned_particles, records)

    TYPE(tParticleData), INTENT(IN) :: owned_particles(:)
    TYPE(tElParticleRecord), ALLOCATABLE, INTENT(OUT) :: records(:)

    REAL*8, ALLOCATABLE :: send_real(:), recv_real(:)
    INTEGER, ALLOCATABLE :: send_int(:), recv_int(:), cursor(:)
    INTEGER :: i, j, total_send, total_recv, ierr, slot

    halo%send_counts = 0
    DO j=1,halo%neighbor_count
      DO i=1,SIZE(owned_particles)
        IF (EL_SPHERE_INTERSECTS_BOX(owned_particles(i)%position, &
            2.0d0*el_kernel_width_factor*owned_particles(i)%radius, &
            halo%boxes(:,halo%neighbors(j)+1))) THEN
          halo%send_counts(j) = halo%send_counts(j) + 1
        END IF
      END DO
    END DO

    CALL MPI_Neighbor_alltoall(halo%send_counts, 1, MPI_INTEGER, &
      halo%recv_counts, 1, MPI_INTEGER, halo%communicator, ierr)
    CALL EL_SET_DISPLACEMENTS(halo%send_counts, halo%send_displs, total_send)
    CALL EL_SET_DISPLACEMENTS(halo%recv_counts, halo%recv_displs, total_recv)

    ALLOCATE(send_real(MAX(1,EL_RECORD_REAL_COUNT*total_send)))
    ALLOCATE(recv_real(MAX(1,EL_RECORD_REAL_COUNT*total_recv)))
    ALLOCATE(send_int(MAX(1,EL_RECORD_INT_COUNT*total_send)))
    ALLOCATE(recv_int(MAX(1,EL_RECORD_INT_COUNT*total_recv)))
    IF (ALLOCATED(halo%send_owner_slots)) DEALLOCATE(halo%send_owner_slots)
    ALLOCATE(halo%send_owner_slots(MAX(1,total_send)))
    ALLOCATE(cursor(MAX(1,halo%neighbor_count)))
    IF (halo%neighbor_count.GT.0) cursor(1:halo%neighbor_count) = halo%send_displs + 1

    DO j=1,halo%neighbor_count
      DO i=1,SIZE(owned_particles)
        IF (.NOT.EL_SPHERE_INTERSECTS_BOX(owned_particles(i)%position, &
            2.0d0*el_kernel_width_factor*owned_particles(i)%radius, &
            halo%boxes(:,halo%neighbors(j)+1))) CYCLE
        slot = cursor(j)
        CALL EL_PACK_RECORD(owned_particles(i), i, slot, send_real, send_int)
        halo%send_owner_slots(slot) = i
        cursor(j) = cursor(j) + 1
      END DO
    END DO

    CALL EL_NEIGHBOR_ALLTOALLV_REAL(send_real, recv_real, EL_RECORD_REAL_COUNT)
    CALL EL_NEIGHBOR_ALLTOALLV_INT(send_int, recv_int, EL_RECORD_INT_COUNT)

    ALLOCATE(records(SIZE(owned_particles)+total_recv))
    DO i=1,SIZE(owned_particles)
      CALL EL_COPY_OWNED_RECORD(owned_particles(i), i, records(i))
    END DO
    DO i=1,total_recv
      CALL EL_UNPACK_RECORD(i, recv_real, recv_int, records(SIZE(owned_particles)+i))
    END DO
    CALL EL_ASSERT_UNIQUE_RECORDS(records)

    DEALLOCATE(send_real, recv_real, send_int, recv_int, cursor)

  END SUBROUTINE EL_EXCHANGE_PARTICLE_RECORDS

  SUBROUTINE EL_ASSERT_UNIQUE_RECORDS(records)

    TYPE(tElParticleRecord), INTENT(IN) :: records(:)
    INTEGER :: i, j, ierr

    DO i=1,SIZE(records)
      IF (records(i)%owner_rank.LT.0) THEN
        WRITE(*,'(A,I0,A,I0)') 'Fatal E-L halo error on rank ', halo%rank, &
          ': particle has invalid owner rank, system ID ', records(i)%system_id
        CALL MPI_Abort(MPI_COMM_SUBS, 902, ierr)
      END IF
      DO j=i+1,SIZE(records)
        IF (records(i)%system_id.EQ.records(j)%system_id) THEN
          WRITE(*,'(A,I0,A,I0)') 'Fatal E-L halo error on rank ', halo%rank, &
            ': duplicate particle record, system ID ', records(i)%system_id
          CALL MPI_Abort(MPI_COMM_SUBS, 903, ierr)
        END IF
      END DO
    END DO

  END SUBROUTINE EL_ASSERT_UNIQUE_RECORDS

  SUBROUTINE EL_REDUCE_TO_OWNERS(local_values, owned_values)

    REAL*8, INTENT(IN) :: local_values(:,:)
    REAL*8, INTENT(OUT) :: owned_values(:,:)

    REAL*8, ALLOCATABLE :: send_buffer(:), recv_buffer(:)
    INTEGER :: ncomp, n_owned, total_send, total_recv, i, slot

    ncomp = SIZE(local_values,1)
    n_owned = SIZE(owned_values,2)
    total_send = SUM(halo%send_counts)
    total_recv = SUM(halo%recv_counts)
    owned_values = local_values(:,1:n_owned)

    ALLOCATE(send_buffer(MAX(1,ncomp*total_recv)))
    ALLOCATE(recv_buffer(MAX(1,ncomp*total_send)))
    DO i=1,total_recv
      send_buffer(ncomp*(i-1)+1:ncomp*i) = local_values(:,n_owned+i)
    END DO
    CALL EL_NEIGHBOR_ALLTOALLV_VALUES(send_buffer, recv_buffer, ncomp, &
      halo%recv_counts, halo%recv_displs, halo%send_counts, halo%send_displs)

    DO i=1,total_send
      slot = halo%send_owner_slots(i)
      owned_values(:,slot) = owned_values(:,slot) + &
        recv_buffer(ncomp*(i-1)+1:ncomp*i)
    END DO
    DEALLOCATE(send_buffer, recv_buffer)

  END SUBROUTINE EL_REDUCE_TO_OWNERS

  SUBROUTINE EL_BROADCAST_FROM_OWNERS(owned_values, record_values)

    REAL*8, INTENT(IN) :: owned_values(:,:)
    REAL*8, INTENT(OUT) :: record_values(:,:)

    REAL*8, ALLOCATABLE :: send_buffer(:), recv_buffer(:)
    INTEGER :: ncomp, n_owned, total_send, total_recv, i, slot

    ncomp = SIZE(owned_values,1)
    n_owned = SIZE(owned_values,2)
    total_send = SUM(halo%send_counts)
    total_recv = SUM(halo%recv_counts)
    record_values(:,1:n_owned) = owned_values

    ALLOCATE(send_buffer(MAX(1,ncomp*total_send)))
    ALLOCATE(recv_buffer(MAX(1,ncomp*total_recv)))
    DO i=1,total_send
      slot = halo%send_owner_slots(i)
      send_buffer(ncomp*(i-1)+1:ncomp*i) = owned_values(:,slot)
    END DO
    CALL EL_NEIGHBOR_ALLTOALLV_VALUES(send_buffer, recv_buffer, ncomp, &
      halo%send_counts, halo%send_displs, halo%recv_counts, halo%recv_displs)
    DO i=1,total_recv
      record_values(:,n_owned+i) = recv_buffer(ncomp*(i-1)+1:ncomp*i)
    END DO
    DEALLOCATE(send_buffer, recv_buffer)

  END SUBROUTINE EL_BROADCAST_FROM_OWNERS

  SUBROUTINE EL_FINALIZE_HALO()

    CALL EL_RELEASE_GRAPH()
    halo%warned_wide_kernel = .FALSE.

  END SUBROUTINE EL_FINALIZE_HALO

  SUBROUTINE EL_RELEASE_GRAPH()

    INTEGER :: ierr

    IF (halo%communicator.NE.MPI_COMM_NULL) CALL MPI_Comm_free(halo%communicator, ierr)
    IF (ALLOCATED(halo%neighbors)) DEALLOCATE(halo%neighbors)
    IF (ALLOCATED(halo%boxes)) DEALLOCATE(halo%boxes)
    IF (ALLOCATED(halo%send_counts)) DEALLOCATE(halo%send_counts)
    IF (ALLOCATED(halo%recv_counts)) DEALLOCATE(halo%recv_counts)
    IF (ALLOCATED(halo%send_displs)) DEALLOCATE(halo%send_displs)
    IF (ALLOCATED(halo%recv_displs)) DEALLOCATE(halo%recv_displs)
    IF (ALLOCATED(halo%send_owner_slots)) DEALLOCATE(halo%send_owner_slots)
    halo%communicator = MPI_COMM_NULL
    halo%initialized = .FALSE.

  END SUBROUTINE EL_RELEASE_GRAPH

  LOGICAL FUNCTION EL_BOXES_WITHIN_DELTA(a, b, delta)

    REAL*8, INTENT(IN) :: a(6), b(6), delta
    REAL*8 :: gap(3)

    gap(1) = MAX(0.0d0, MAX(a(1)-b(2), b(1)-a(2)))
    gap(2) = MAX(0.0d0, MAX(a(3)-b(4), b(3)-a(4)))
    gap(3) = MAX(0.0d0, MAX(a(5)-b(6), b(5)-a(6)))
    EL_BOXES_WITHIN_DELTA = SQRT(SUM(gap*gap)).LE.delta

  END FUNCTION EL_BOXES_WITHIN_DELTA

  LOGICAL FUNCTION EL_SPHERE_INTERSECTS_BOX(position, radius, box)

    REAL*8, INTENT(IN) :: position(3), radius, box(6)
    REAL*8 :: gap(3)

    gap(1) = MAX(0.0d0, MAX(box(1)-position(1), position(1)-box(2)))
    gap(2) = MAX(0.0d0, MAX(box(3)-position(2), position(2)-box(4)))
    gap(3) = MAX(0.0d0, MAX(box(5)-position(3), position(3)-box(6)))
    EL_SPHERE_INTERSECTS_BOX = SUM(gap*gap).LE.radius*radius

  END FUNCTION EL_SPHERE_INTERSECTS_BOX

  SUBROUTINE EL_SET_DISPLACEMENTS(counts, displacements, total)

    INTEGER, INTENT(IN) :: counts(:)
    INTEGER, INTENT(OUT) :: displacements(:), total
    INTEGER :: i

    total = 0
    DO i=1,SIZE(counts)
      displacements(i) = total
      total = total + counts(i)
    END DO

  END SUBROUTINE EL_SET_DISPLACEMENTS

  SUBROUTINE EL_PACK_RECORD(particle, owner_slot, slot, real_buffer, int_buffer)

    TYPE(tParticleData), INTENT(IN) :: particle
    INTEGER, INTENT(IN) :: owner_slot, slot
    REAL*8, INTENT(INOUT) :: real_buffer(:)
    INTEGER, INTENT(INOUT) :: int_buffer(:)
    INTEGER :: r0, i0

    r0 = EL_RECORD_REAL_COUNT*(slot-1)
    i0 = EL_RECORD_INT_COUNT*(slot-1)
    real_buffer(r0+1:r0+3) = particle%position
    real_buffer(r0+4:r0+6) = particle%velocity
    real_buffer(r0+7) = particle%radius
    real_buffer(r0+8) = particle%density
    int_buffer(i0+1) = particle%systemIdx
    int_buffer(i0+2) = halo%rank
    int_buffer(i0+3) = owner_slot

  END SUBROUTINE EL_PACK_RECORD

  SUBROUTINE EL_COPY_OWNED_RECORD(particle, owner_slot, record)

    TYPE(tParticleData), INTENT(IN) :: particle
    INTEGER, INTENT(IN) :: owner_slot
    TYPE(tElParticleRecord), INTENT(OUT) :: record

    record%position = particle%position
    record%velocity = particle%velocity
    record%radius = particle%radius
    record%density = particle%density
    record%system_id = particle%systemIdx
    record%owner_rank = halo%rank
    record%owner_slot = owner_slot
    record%owned = .TRUE.

  END SUBROUTINE EL_COPY_OWNED_RECORD

  SUBROUTINE EL_UNPACK_RECORD(slot, real_buffer, int_buffer, record)

    INTEGER, INTENT(IN) :: slot
    REAL*8, INTENT(IN) :: real_buffer(:)
    INTEGER, INTENT(IN) :: int_buffer(:)
    TYPE(tElParticleRecord), INTENT(OUT) :: record
    INTEGER :: r0, i0

    r0 = EL_RECORD_REAL_COUNT*(slot-1)
    i0 = EL_RECORD_INT_COUNT*(slot-1)
    record%position = real_buffer(r0+1:r0+3)
    record%velocity = real_buffer(r0+4:r0+6)
    record%radius = real_buffer(r0+7)
    record%density = real_buffer(r0+8)
    record%system_id = int_buffer(i0+1)
    record%owner_rank = int_buffer(i0+2)
    record%owner_slot = int_buffer(i0+3)
    record%owned = .FALSE.

  END SUBROUTINE EL_UNPACK_RECORD

  SUBROUTINE EL_NEIGHBOR_ALLTOALLV_REAL(send_buffer, recv_buffer, width)

    REAL*8, INTENT(IN) :: send_buffer(:)
    REAL*8, INTENT(OUT) :: recv_buffer(:)
    INTEGER, INTENT(IN) :: width
    INTEGER, ALLOCATABLE :: sc(:), rc(:), sd(:), rd(:)
    INTEGER :: ierr

    ALLOCATE(sc(halo%neighbor_count), rc(halo%neighbor_count))
    ALLOCATE(sd(halo%neighbor_count), rd(halo%neighbor_count))
    sc = width*halo%send_counts
    rc = width*halo%recv_counts
    sd = width*halo%send_displs
    rd = width*halo%recv_displs
    CALL MPI_Neighbor_alltoallv(send_buffer, sc, sd, MPI_DOUBLE_PRECISION, &
      recv_buffer, rc, rd, MPI_DOUBLE_PRECISION, halo%communicator, ierr)
    DEALLOCATE(sc, rc, sd, rd)

  END SUBROUTINE EL_NEIGHBOR_ALLTOALLV_REAL

  SUBROUTINE EL_NEIGHBOR_ALLTOALLV_INT(send_buffer, recv_buffer, width)

    INTEGER, INTENT(IN) :: send_buffer(:)
    INTEGER, INTENT(OUT) :: recv_buffer(:)
    INTEGER, INTENT(IN) :: width
    INTEGER, ALLOCATABLE :: sc(:), rc(:), sd(:), rd(:)
    INTEGER :: ierr

    ALLOCATE(sc(halo%neighbor_count), rc(halo%neighbor_count))
    ALLOCATE(sd(halo%neighbor_count), rd(halo%neighbor_count))
    sc = width*halo%send_counts
    rc = width*halo%recv_counts
    sd = width*halo%send_displs
    rd = width*halo%recv_displs
    CALL MPI_Neighbor_alltoallv(send_buffer, sc, sd, MPI_INTEGER, &
      recv_buffer, rc, rd, MPI_INTEGER, halo%communicator, ierr)
    DEALLOCATE(sc, rc, sd, rd)

  END SUBROUTINE EL_NEIGHBOR_ALLTOALLV_INT

  SUBROUTINE EL_NEIGHBOR_ALLTOALLV_VALUES(send_buffer, recv_buffer, width, &
      send_counts, send_displs, recv_counts, recv_displs)

    REAL*8, INTENT(IN) :: send_buffer(:)
    REAL*8, INTENT(OUT) :: recv_buffer(:)
    INTEGER, INTENT(IN) :: width
    INTEGER, INTENT(IN) :: send_counts(:), send_displs(:)
    INTEGER, INTENT(IN) :: recv_counts(:), recv_displs(:)
    INTEGER, ALLOCATABLE :: sc(:), rc(:), sd(:), rd(:)
    INTEGER :: ierr

    ALLOCATE(sc(halo%neighbor_count), rc(halo%neighbor_count))
    ALLOCATE(sd(halo%neighbor_count), rd(halo%neighbor_count))
    sc = width*send_counts
    rc = width*recv_counts
    sd = width*send_displs
    rd = width*recv_displs
    CALL MPI_Neighbor_alltoallv(send_buffer, sc, sd, MPI_DOUBLE_PRECISION, &
      recv_buffer, rc, rd, MPI_DOUBLE_PRECISION, halo%communicator, ierr)
    DEALLOCATE(sc, rc, sd, rd)

  END SUBROUTINE EL_NEIGHBOR_ALLTOALLV_VALUES

END MODULE EL_HALO
