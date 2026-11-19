!=========================================================================
! EL_VELOCITY_AUDIT - where is the fastest fluid, and what is the porosity there (case 11)
!
! SimPar@ELVelocityAuditFreq = N (EL_CONFIG; default 0 = off): every N steps, after the fluid
! solve, one line on the protocol and the terminal (rank showid):
!
!   EL_UMAX_AUDIT t= step= umax= |u| u v w at= x y z eps_h= e dragB_fsrc= B f
!                 eps_min= e at= x y z inlet_eps_min= e inlet_wmax= w
!
!   umax      largest Euclidean norm of the Q2 velocity over all dofs of all ranks, its
!             components and the coordinates of that dof
!   eps_h     transport porosity (Q1 reconstruction embedded in Q2) at that dof, -1 when the
!             field is not built (it is built under the conservative continuity form)
!   dragB     lumped semi-implicit drag coefficient of the dof (drag_B_ml), fsrc the z
!             component of the deposited explicit particle force at the dof (force_rhs); -1 / 0
!             when the arrays do not have the dof layout
!   eps_min   smallest element porosity epsilon_f and the centre of that element
!   inlet_*   smallest transport porosity and largest w on the bottom plane z = zmin of the EL
!             domain box (the inlet of case 11: w = U / eps_h there)
!
! Diagnostic only. Collectives on MPI_COMM_SUBS: called by every worker (never by the master).
!=========================================================================
MODULE EL_VELOCITY_AUDIT

  USE MPI
  USE PP3D_MPI, ONLY: myid, showid, MPI_COMM_SUBS
  USE TYPES, ONLY: tMultiMesh
  USE EL_CONFIG, ONLY: el_velocity_audit_freq
  USE EL_GEOMETRY, ONLY: el_domain_box, el_domain_box_set
  USE EL_FIELDS, ONLY: el_field_data

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: EL_WRITE_VELOCITY_AUDIT

CONTAINS

  SUBROUTINE EL_WRITE_VELOCITY_AUDIT(u, v, w, coor, ndof, time, mfile, istep, mesh, ilev)
    INTEGER, INTENT(IN) :: ndof, mfile, istep, ilev
    REAL*8, INTENT(IN) :: u(:), v(:), w(:), coor(:,:), time
    TYPE(tMultiMesh), INTENT(IN) :: mesh

    REAL*8 :: pair(2), gpair(2), rec(10), grec(10), zmin, tol, s, emin_loc(2), emin(2)
    REAL*8 :: erec(3), gerec(3), inl(2), ginl(2)
    INTEGER :: i, imax, ierr, owner, nel, iemin

    IF (el_velocity_audit_freq.LE.0) RETURN
    IF (MOD(ABS(istep), el_velocity_audit_freq).NE.0) RETURN

    ! fastest dof of this rank
    pair = (/ -1.0d0, DBLE(myid) /)
    imax = 0
    DO i = 1, ndof
      s = u(i)*u(i) + v(i)*v(i) + w(i)*w(i)
      IF (s.GT.pair(1)) THEN
        pair(1) = s
        imax = i
      END IF
    END DO
    CALL MPI_Allreduce(pair, gpair, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, MPI_COMM_SUBS, ierr)
    owner = NINT(gpair(2))
    rec = 0.0d0
    IF (myid.EQ.owner .AND. imax.GT.0) THEN
      rec(1) = SQRT(pair(1))
      rec(2:4) = (/ u(imax), v(imax), w(imax) /)
      rec(5:7) = coor(1:3, imax)
      rec(8) = -1.0d0
      IF (ALLOCATED(el_field_data%epsilon_transport_q2)) THEN
        IF (SIZE(el_field_data%epsilon_transport_q2).GE.imax) rec(8) = el_field_data%epsilon_transport_q2(imax)
      END IF
      rec(9) = -1.0d0
      IF (ALLOCATED(el_field_data%drag_B_ml)) THEN
        IF (SIZE(el_field_data%drag_B_ml).GE.imax) rec(9) = el_field_data%drag_B_ml(imax)
      END IF
      IF (ALLOCATED(el_field_data%force_rhs)) THEN
        IF (SIZE(el_field_data%force_rhs,2).GE.imax) rec(10) = el_field_data%force_rhs(3,imax)
      END IF
    END IF
    CALL MPI_Allreduce(rec, grec, 10, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)

    ! smallest element porosity
    nel = mesh%level(ilev)%nel
    emin_loc = (/ HUGE(1.0d0), DBLE(myid) /)
    iemin = 0
    IF (ALLOCATED(el_field_data%epsilon_f)) THEN
      DO i = 1, MIN(nel, SIZE(el_field_data%epsilon_f))
        IF (el_field_data%epsilon_f(i).LT.emin_loc(1)) THEN
          emin_loc(1) = el_field_data%epsilon_f(i)
          iemin = i
        END IF
      END DO
    END IF
    CALL MPI_Allreduce(emin_loc, emin, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_SUBS, ierr)
    erec = 0.0d0
    IF (myid.EQ.NINT(emin(2)) .AND. iemin.GT.0) THEN
      DO i = 1, 8
        erec = erec + 0.125d0*mesh%level(ilev)%dcorvg(:, mesh%level(ilev)%kvert(i, iemin))
      END DO
    END IF
    CALL MPI_Allreduce(erec, gerec, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_SUBS, ierr)

    ! inlet plane
    inl = (/ -HUGE(1.0d0), -HUGE(1.0d0) /)   ! (-eps_min, w_max)
    IF (el_domain_box_set) THEN
      zmin = el_domain_box(5)
      tol = 1.0d-6*(el_domain_box(6) - el_domain_box(5))
      DO i = 1, ndof
        IF (ABS(coor(3,i) - zmin).GE.tol) CYCLE
        inl(2) = MAX(inl(2), w(i))
        IF (ALLOCATED(el_field_data%epsilon_transport_q2)) THEN
          IF (SIZE(el_field_data%epsilon_transport_q2).GE.i) &
            inl(1) = MAX(inl(1), -el_field_data%epsilon_transport_q2(i))
        END IF
      END DO
    END IF
    CALL MPI_Allreduce(inl, ginl, 2, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_SUBS, ierr)

    IF (myid.EQ.showid) THEN
      CALL WRITE_LINE(6)
      CALL WRITE_LINE(mfile)
    END IF

  CONTAINS

    SUBROUTINE WRITE_LINE(unit)
      INTEGER, INTENT(IN) :: unit
      WRITE(unit,'(A,ES14.6,A,I0,A,4ES14.6,A,3ES14.6,A,ES14.6,A,2ES14.6,A,ES14.6,A,3ES14.6,A,ES14.6,A,ES14.6)') &
        'EL_UMAX_AUDIT t= ', time, ' step= ', istep, ' umax= ', grec(1:4), ' at= ', grec(5:7), &
        ' eps_h= ', grec(8), ' dragB_fsrc= ', grec(9:10), ' eps_min= ', emin(1), ' at= ', gerec, &
        ' inlet_eps_min= ', -ginl(1), ' inlet_wmax= ', ginl(2)
    END SUBROUTINE WRITE_LINE

  END SUBROUTINE EL_WRITE_VELOCITY_AUDIT

END MODULE EL_VELOCITY_AUDIT
