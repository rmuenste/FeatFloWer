!=========================================================================
! CHI_EXCHANGE - MPI service for the Chimera coupling, Layer M
! (design: chimera-integration-design.md v3, sections 3 and 4).
!
! Contract: every rank of the communicator calls with the IDENTICAL
! replicated query list (equal counts by construction).  Each rank
! locates the points in its own partition; the owner of a point is the
! lowest rank that found it (MPI_Allreduce MIN over proposed ranks);
! values are reduced by MPI_Allreduce SUM over owner-masked
! contributions, so every rank ends up with the same bits.  Points that
! no rank found are reported (nmissing) - the caller decides (fatal for
! the coupling: a submesh outer surface must lie inside the background
! domain).
!
! The communicator is passed in (the coupling layer uses the worker
! communicator MPI_COMM_SUBS; the master never calls here).  This
! module is COMMON-free and var_QuadScalar-free.
!=========================================================================
MODULE CHI_EXCHANGE

  USE CHI_LOCATOR, ONLY: tChimeraLocator, CHI_LOCATE
  USE CHI_FEM_EVAL, ONLY: CHI_EVAL_FIELD_AT

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  PRIVATE

  PUBLIC :: CHI_EXCHANGE_BG_EVAL
  PUBLIC :: CHI_EXCHANGE_BCAST
  PUBLIC :: CHI_EXCHANGE_RANK
  PUBLIC :: CHI_EXCHANGE_SIZE
  PUBLIC :: CHI_EXCHANGE_MAX_INT
  PUBLIC :: CHI_BG_NVAL

  ! value layout per point: u(3), grad u (3x3, column-major: (a,b) at
  ! 3 + a + 3*(b-1)), p
  INTEGER, PARAMETER :: CHI_BG_NVAL = 13

CONTAINS

  INTEGER FUNCTION CHI_EXCHANGE_RANK(comm)
    INTEGER, INTENT(IN) :: comm
    INTEGER :: ierr
    CALL MPI_COMM_RANK(comm, CHI_EXCHANGE_RANK, ierr)
  END FUNCTION CHI_EXCHANGE_RANK

  INTEGER FUNCTION CHI_EXCHANGE_SIZE(comm)
    INTEGER, INTENT(IN) :: comm
    INTEGER :: ierr
    CALL MPI_COMM_SIZE(comm, CHI_EXCHANGE_SIZE, ierr)
  END FUNCTION CHI_EXCHANGE_SIZE

  !-----------------------------------------------------------------------
  ! Evaluate a Q2/P1 field, distributed over the communicator, at the
  ! replicated point list pts(3,npts).  On return vals(CHI_BG_NVAL,npts)
  ! is identical on all ranks; owner(npts) is the owning rank (or the
  ! communicator size for points nobody found); nmissing counts those.
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_EXCHANGE_BG_EVAL(comm, npts, pts, loc, dcorvg, &
                                  kvert, kedge, karea, nvt, net, nat, &
                                  u, v, w, p, vals, owner, nmissing)
    INTEGER, INTENT(IN) :: comm, npts
    REAL*8,  INTENT(IN) :: pts(3,*)
    TYPE(tChimeraLocator), INTENT(IN) :: loc
    REAL*8,  INTENT(IN) :: dcorvg(3,*)
    INTEGER, INTENT(IN) :: kvert(8,*), kedge(12,*), karea(6,*)
    INTEGER, INTENT(IN) :: nvt, net, nat
    REAL*8,  INTENT(IN) :: u(*), v(*), w(*), p(*)
    REAL*8,  INTENT(OUT) :: vals(CHI_BG_NVAL,*)
    INTEGER, INTENT(OUT) :: owner(*)
    INTEGER, INTENT(OUT) :: nmissing

    INTEGER :: myrank, nranks, ierr, ip, i
    INTEGER, ALLOCATABLE :: prop(:), iel(:)
    REAL*8,  ALLOCATABLE :: xi(:,:)
    REAL*8 :: uval(3), gradu(3,3), pval
    LOGICAL :: found, ok

    CALL MPI_COMM_RANK(comm, myrank, ierr)
    CALL MPI_COMM_SIZE(comm, nranks, ierr)

    nmissing = 0
    IF (npts .LE. 0) RETURN

    ALLOCATE(prop(npts), iel(npts), xi(3,npts))
    DO ip = 1, npts
      CALL CHI_LOCATE(loc, dcorvg, kvert, pts(:,ip), iel(ip), xi(:,ip), found)
      IF (found) THEN
        prop(ip) = myrank
      ELSE
        prop(ip) = nranks
      END IF
    END DO

    CALL MPI_ALLREDUCE(prop, owner, npts, MPI_INTEGER, MPI_MIN, comm, ierr)

    DO ip = 1, npts
      vals(:,ip) = 0d0
      IF (owner(ip) .EQ. myrank) THEN
        CALL CHI_EVAL_FIELD_AT(iel(ip), xi(:,ip), kvert, kedge, karea, &
          nvt, net, nat, dcorvg, u, v, w, p, uval, gradu, pval, ok)
        IF (.NOT. ok) THEN
          WRITE(*,*) 'CHI_EXCHANGE_BG_EVAL: singular Jacobian at point', ip
          STOP 1
        END IF
        vals(1:3,ip) = uval
        DO i = 1, 3
          vals(3+3*(i-1)+1:3+3*i, ip) = gradu(:,i)
        END DO
        vals(13,ip) = pval
      END IF
    END DO

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, vals, CHI_BG_NVAL*npts, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)

    DO ip = 1, npts
      IF (owner(ip) .GE. nranks) nmissing = nmissing + 1
    END DO
    DEALLOCATE(prop, iel, xi)
  END SUBROUTINE CHI_EXCHANGE_BG_EVAL

  !-----------------------------------------------------------------------
  ! Broadcast a real array from root (rank in comm).
  !-----------------------------------------------------------------------
  SUBROUTINE CHI_EXCHANGE_BCAST(comm, root, arr, n)
    INTEGER, INTENT(IN) :: comm, root, n
    REAL*8, INTENT(INOUT) :: arr(*)
    INTEGER :: ierr
    IF (n .LE. 0) RETURN
    CALL MPI_BCAST(arr, n, MPI_DOUBLE_PRECISION, root, comm, ierr)
  END SUBROUTINE CHI_EXCHANGE_BCAST

  !-----------------------------------------------------------------------
  ! Global integer maximum (diagnostics / consistency checks).
  !-----------------------------------------------------------------------
  INTEGER FUNCTION CHI_EXCHANGE_MAX_INT(comm, val)
    INTEGER, INTENT(IN) :: comm, val
    INTEGER :: ierr, res
    CALL MPI_ALLREDUCE(val, res, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
    CHI_EXCHANGE_MAX_INT = res
  END FUNCTION CHI_EXCHANGE_MAX_INT

END MODULE CHI_EXCHANGE
