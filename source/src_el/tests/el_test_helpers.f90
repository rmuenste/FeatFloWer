MODULE EL_TEST_HELPERS

  USE MPI
  USE TYPES, ONLY: tMultiMesh

  IMPLICIT NONE

  INTEGER :: ier, icheck
  INTEGER :: nel, nvt, net, nat, nve, nee, nae, nvel, neel, nved, &
             nvar, near, nbct, nvbd, nebd, nabd
  COMMON /ERRCTL/ ier, icheck
  COMMON /TRIAD/ nel, nvt, net, nat, nve, nee, nae, nvel, neel, nved, &
                 nvar, near, nbct, nvbd, nebd, nabd
  SAVE /ERRCTL/, /TRIAD/

  INTEGER :: helper_ncell = 0

CONTAINS

  SUBROUTINE EL_TEST_INIT()

    icheck = 0
    ier = 0

  END SUBROUTINE EL_TEST_INIT

  SUBROUTINE EL_TEST_RELEASE_MESH(m)

    TYPE(tMultiMesh), INTENT(INOUT) :: m

    IF (ALLOCATED(m%level)) THEN
      IF (SIZE(m%level).GE.1) THEN
        IF (ASSOCIATED(m%level(1)%dcorvg)) DEALLOCATE(m%level(1)%dcorvg)
        IF (ALLOCATED(m%level(1)%kvert)) DEALLOCATE(m%level(1)%kvert)
        IF (ALLOCATED(m%level(1)%kedge)) DEALLOCATE(m%level(1)%kedge)
        IF (ALLOCATED(m%level(1)%karea)) DEALLOCATE(m%level(1)%karea)
        IF (ALLOCATED(m%level(1)%dvol)) DEALLOCATE(m%level(1)%dvol)
      END IF
      DEALLOCATE(m%level)
    END IF

  END SUBROUTINE EL_TEST_RELEASE_MESH

  ! Synthetic structured Q2/P1 hex mesh on the unit cube.
  ! Sets /TRIAD/ and fills the connectivity the EL transfer routines and
  ! NDFGL consume: kvert, kedge, karea, dcorvg, dvol.
  SUBROUTINE BUILD_STRUCTURED_MESH(m, n)

    TYPE(tMultiMesh), INTENT(INOUT) :: m
    INTEGER, INTENT(IN) :: n
    INTEGER :: ie, je, ke, iel, ex, ey, ez, fx, fy, fz
    REAL*8  :: dh

    CALL EL_TEST_RELEASE_MESH(m)
    helper_ncell = n
    dh = 1.0d0 / DBLE(n)

    nel = n*n*n
    nvt = (n+1)*(n+1)*(n+1)
    ex = n*(n+1)*(n+1)
    ey = (n+1)*n*(n+1)
    ez = (n+1)*(n+1)*n
    net = ex + ey + ez
    fx = (n+1)*n*n
    fy = n*(n+1)*n
    fz = n*n*(n+1)
    nat = fx + fy + fz
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
    ALLOCATE(m%level(1)%dcorvg(3, nvt))
    ALLOCATE(m%level(1)%kvert(8, nel))
    ALLOCATE(m%level(1)%kedge(12, nel))
    ALLOCATE(m%level(1)%karea(6, nel))
    ALLOCATE(m%level(1)%dvol(nel))

    DO ke = 0, n
      DO je = 0, n
        DO ie = 0, n
          m%level(1)%dcorvg(:, vid(ie, je, ke)) = &
            (/DBLE(ie)*dh, DBLE(je)*dh, DBLE(ke)*dh/)
        END DO
      END DO
    END DO

    DO ke = 0, n-1
      DO je = 0, n-1
        DO ie = 0, n-1
          iel = 1 + ie + n*(je + n*ke)
          m%level(1)%dvol(iel) = dh*dh*dh

          m%level(1)%kvert(1, iel) = vid(ie,   je,   ke  )
          m%level(1)%kvert(2, iel) = vid(ie+1, je,   ke  )
          m%level(1)%kvert(3, iel) = vid(ie+1, je+1, ke  )
          m%level(1)%kvert(4, iel) = vid(ie,   je+1, ke  )
          m%level(1)%kvert(5, iel) = vid(ie,   je,   ke+1)
          m%level(1)%kvert(6, iel) = vid(ie+1, je,   ke+1)
          m%level(1)%kvert(7, iel) = vid(ie+1, je+1, ke+1)
          m%level(1)%kvert(8, iel) = vid(ie,   je+1, ke+1)

          m%level(1)%kedge(1,  iel) = edx(ie, je,   ke  )
          m%level(1)%kedge(2,  iel) = edy(ie+1, je, ke  )
          m%level(1)%kedge(3,  iel) = edx(ie, je+1, ke  )
          m%level(1)%kedge(4,  iel) = edy(ie, je,   ke  )
          m%level(1)%kedge(5,  iel) = edz(ie,   je,   ke)
          m%level(1)%kedge(6,  iel) = edz(ie+1, je,   ke)
          m%level(1)%kedge(7,  iel) = edz(ie+1, je+1, ke)
          m%level(1)%kedge(8,  iel) = edz(ie,   je+1, ke)
          m%level(1)%kedge(9,  iel) = edx(ie, je,   ke+1)
          m%level(1)%kedge(10, iel) = edy(ie+1, je, ke+1)
          m%level(1)%kedge(11, iel) = edx(ie, je+1, ke+1)
          m%level(1)%kedge(12, iel) = edy(ie, je,   ke+1)

          m%level(1)%karea(1, iel) = faz(ie, je, ke  )
          m%level(1)%karea(2, iel) = fay(ie, je, ke  )
          m%level(1)%karea(3, iel) = fax(ie+1, je, ke)
          m%level(1)%karea(4, iel) = fay(ie, je+1, ke)
          m%level(1)%karea(5, iel) = fax(ie, je, ke  )
          m%level(1)%karea(6, iel) = faz(ie, je, ke+1)
        END DO
      END DO
    END DO

  END SUBROUTINE BUILD_STRUCTURED_MESH

  SUBROUTINE BUILD_DOF_COORDINATES(m, c)

    TYPE(tMultiMesh), INTENT(IN) :: m
    REAL*8, INTENT(OUT) :: c(:,:)
    INTEGER :: iel, i, g, a, b
    INTEGER, PARAMETER :: ed(2,12) = RESHAPE((/ &
      1,2, 2,3, 4,3, 1,4, 1,5, 2,6, 3,7, 4,8, 5,6, 6,7, 8,7, 5,8 /), (/2,12/))
    INTEGER, PARAMETER :: fc(4,6) = RESHAPE((/ &
      1,2,3,4,  1,2,6,5,  2,3,7,6,  4,3,7,8,  1,4,8,5,  5,6,7,8 /), (/4,6/))
    REAL*8 :: v(3,8)

    c = 0.0d0
    DO iel = 1, m%level(1)%nel
      DO i = 1, 8
        v(:, i) = m%level(1)%dcorvg(:, m%level(1)%kvert(i, iel))
        c(:, m%level(1)%kvert(i, iel)) = v(:, i)
      END DO
      DO i = 1, 12
        a = ed(1, i)
        b = ed(2, i)
        g = m%level(1)%nvt + m%level(1)%kedge(i, iel)
        c(:, g) = 0.5d0*(v(:, a) + v(:, b))
      END DO
      DO i = 1, 6
        g = m%level(1)%nvt + m%level(1)%net + m%level(1)%karea(i, iel)
        c(:, g) = 0.25d0*(v(:, fc(1,i)) + v(:, fc(2,i)) + &
                          v(:, fc(3,i)) + v(:, fc(4,i)))
      END DO
      g = m%level(1)%nvt + m%level(1)%net + m%level(1)%nat + iel
      c(:, g) = 0.125d0*SUM(v, DIM=2)
    END DO

  END SUBROUTINE BUILD_DOF_COORDINATES

  SUBROUTINE EL_TEST_CHECK(ok, label, nf, rank)

    LOGICAL, INTENT(IN) :: ok
    CHARACTER(LEN=*), INTENT(IN) :: label
    INTEGER, INTENT(INOUT) :: nf
    INTEGER, INTENT(IN) :: rank

    IF (.NOT.ok) THEN
      nf = nf + 1
      IF (rank.EQ.0) WRITE(*,'(A,A)') '  [FAIL] ', label
    ELSE
      IF (rank.EQ.0) WRITE(*,'(A,A)') '  [ ok ] ', label
    END IF

  END SUBROUTINE EL_TEST_CHECK

  SUBROUTINE EL_TEST_REPORT(nf, rank, nproc, suite_name)

    INTEGER, INTENT(IN) :: nf, rank, nproc
    CHARACTER(LEN=*), INTENT(IN) :: suite_name
    INTEGER :: total, ierr

    CALL MPI_Allreduce(nf, total, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    IF (rank.EQ.0) THEN
      IF (total.EQ.0) THEN
        WRITE(*,'(A,A,I0,A)') TRIM(suite_name), ' tests passed on ', nproc, &
          ' rank(s).'
      ELSE
        WRITE(*,'(A,A,I0,A,I0,A)') TRIM(suite_name), ' tests FAILED: ', total, &
          ' failure(s) across ', nproc, ' rank(s).'
      END IF
    END IF

  END SUBROUTINE EL_TEST_REPORT

  INTEGER FUNCTION vid(i, j, k)
    INTEGER, INTENT(IN) :: i, j, k
    vid = 1 + i + (helper_ncell+1)*(j + (helper_ncell+1)*k)
  END FUNCTION vid

  INTEGER FUNCTION edx(i, j, k)
    INTEGER, INTENT(IN) :: i, j, k
    edx = 1 + i + helper_ncell*(j + (helper_ncell+1)*k)
  END FUNCTION edx

  INTEGER FUNCTION edy(i, j, k)
    INTEGER, INTENT(IN) :: i, j, k
    edy = helper_ncell*(helper_ncell+1)*(helper_ncell+1) + &
      1 + i + (helper_ncell+1)*(j + helper_ncell*k)
  END FUNCTION edy

  INTEGER FUNCTION edz(i, j, k)
    INTEGER, INTENT(IN) :: i, j, k
    edz = helper_ncell*(helper_ncell+1)*(helper_ncell+1) + &
      (helper_ncell+1)*helper_ncell*(helper_ncell+1) + &
      1 + i + (helper_ncell+1)*(j + (helper_ncell+1)*k)
  END FUNCTION edz

  INTEGER FUNCTION fax(i, j, k)
    INTEGER, INTENT(IN) :: i, j, k
    fax = 1 + i + (helper_ncell+1)*(j + helper_ncell*k)
  END FUNCTION fax

  INTEGER FUNCTION fay(i, j, k)
    INTEGER, INTENT(IN) :: i, j, k
    fay = (helper_ncell+1)*helper_ncell*helper_ncell + &
      1 + i + helper_ncell*(j + (helper_ncell+1)*k)
  END FUNCTION fay

  INTEGER FUNCTION faz(i, j, k)
    INTEGER, INTENT(IN) :: i, j, k
    faz = (helper_ncell+1)*helper_ncell*helper_ncell + &
      helper_ncell*(helper_ncell+1)*helper_ncell + &
      1 + i + helper_ncell*(j + helper_ncell*k)
  END FUNCTION faz

END MODULE EL_TEST_HELPERS
