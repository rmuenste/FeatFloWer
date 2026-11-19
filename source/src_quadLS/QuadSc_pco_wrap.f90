!=========================================================================
! QuadSc_pco_wrap.f90 - master-side assembly of the periodic coarse pressure operator
!
! Part of the periodic coarse pressure operator (overview: QuadSc_pco_operator.f90).
!
! The own matrices of the master B (velocity rows, P1 pressure columns), B^T (pressure rows,
! velocity columns) and lumped masses live on the unpartitioned mesh, where the two images
! of a periodic velocity dof are separate dofs. With Q the map that copies the value of a
! periodic class to all its images (classes: pco_periodic_topology), the periodic operator is
!     C = (B^T Q) M_cls^-1 (Q^T B):
! the B rows of the images of a class are summed into one class row, the B^T columns into one
! class column, and every class has one mass and one constraint mask. Per class k:
!   ordinary mass        m_k = sum_i M_i                     (M = mg_MlRhoPmat)
!   weighted EL mass     m_k = (sum_i U_i)^2 / (sum_i W_i)   (U = rho Ml, W = lump(eps rho);
!                        the workers sum U and W over coupled dofs first and transform once,
!                        QuadSc_el_mm.f90 - transforming per image and summing is wrong at a
!                        seam with unequal porosities)
!   mask                 the Dirichlet flags of the class (equal over the class, checked by the
!                        coordinator before assembly)
!   C(p,q) += sum_c BT_c(p,k) B_c(k,q) / m_k  over the unconstrained components c.
! All pressure columns of a class row are merged before the product (no first-match loss of
! terms when one pressure element touches two images of the same class), every image
! counts. The slip projection is in B and B^T as the master stores them (SetSlipOnBandBT).
!
! Structure (PCO_WRAP_BUILD, once per pattern): per class the distinct pressure elements of
! its B rows and B^T columns and their block positions in the master layout; the column
! index of B^T per velocity dof. Values (PCO_WRAP_ASSEMBLE, every renewal): O(nnz(B)).
! Master only. The result is the oracle-checked alternative to the gathered operator.
!=========================================================================
MODULE pco_wrap

  USE var_QuadScalar, ONLY: mg_qlMat, mg_lqMat, mg_BXMat, mg_BYMat, mg_BZMat, &
                            mg_BTXMat, mg_BTYMat, mg_BTZMat
  USE pco_periodic_topology, ONLY: tPcoClasses
  USE pco_coarse_solver, ONLY: PCO_NEIGHBOUR_POS, PCO_BLOCK_SLOT

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: PCO_MASS_ORDINARY = 1, PCO_MASS_EL_WEIGHTED = 2
  ! check-only variant of the weighted mass (sum_i U_i^2/W_i, the seam mass of the diagnostic
  ! prototype); used by the check of the coordinator to show the check discriminates
  INTEGER, PARAMETER, PUBLIC :: PCO_MASS_EL_PER_IMAGE = 3

  TYPE, PUBLIC :: tPcoWrap
    LOGICAL :: built = .FALSE.
    INTEGER :: lev = 0, ndof = 0
    INTEGER, ALLOCATABLE :: bt_ptr(:), bt_pos(:)     ! velocity dof -> slots of its B^T column
    INTEGER, ALLOCATABLE :: bt_prow(:)               ! ... and their pressure rows
    INTEGER :: maxne = 0                             ! largest element count of a class
    INTEGER, ALLOCATABLE :: el_ptr(:), el(:)         ! class -> pressure elements touched
    INTEGER, ALLOCATABLE :: bp_ptr(:), bp(:)         ! class -> block position (a,b) row-major
    INTEGER :: nbad = 0                              ! element pairs outside the layout
  END TYPE tPcoWrap

  PUBLIC :: PCO_WRAP_BUILD, PCO_WRAP_ASSEMBLE

CONTAINS

  !-----------------------------------------------------------------------
  ! Structure of the wrapped assembly on level lev (master). w%nbad > 0: the classes couple
  ! elements the layout does not hold (refused by the coordinator).
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_WRAP_BUILD(w, lev, cls)
    TYPE(tPcoWrap), INTENT(INOUT) :: w
    INTEGER, INTENT(IN) :: lev
    TYPE(tPcoClasses), INTENT(IN) :: cls
    INTEGER :: ndof, nprow, i, k, p, j, m, ne, a, b, pos, ip
    INTEGER, ALLOCATABLE :: fill(:), stamp(:), buf(:)

    w%lev = lev
    ndof = cls%ndof
    w%ndof = ndof
    nprow = mg_lqMat(lev)%nu

    ! column access of B^T: slots of the entries with velocity column i
    IF (ALLOCATED(w%bt_ptr)) DEALLOCATE(w%bt_ptr, w%bt_pos, w%bt_prow)
    ALLOCATE(w%bt_ptr(ndof+1), fill(ndof+1))
    fill = 0
    DO k = 1, mg_lqMat(lev)%LdA(nprow+1)-1
      i = mg_lqMat(lev)%ColA(k)
      fill(i+1) = fill(i+1) + 1
    END DO
    w%bt_ptr(1) = 1
    DO i = 1, ndof
      w%bt_ptr(i+1) = w%bt_ptr(i) + fill(i+1)
    END DO
    ALLOCATE(w%bt_pos(MAX(w%bt_ptr(ndof+1)-1,1)), w%bt_prow(MAX(w%bt_ptr(ndof+1)-1,1)))
    fill(1:ndof) = w%bt_ptr(1:ndof)
    DO p = 1, nprow
      DO k = mg_lqMat(lev)%LdA(p), mg_lqMat(lev)%LdA(p+1)-1
        i = mg_lqMat(lev)%ColA(k)
        w%bt_pos(fill(i)) = k
        w%bt_prow(fill(i)) = p
        fill(i) = fill(i) + 1
      END DO
    END DO
    DEALLOCATE(fill)

    ! pressure elements per class: union over the members of the B row and B^T column
    IF (ALLOCATED(w%el_ptr)) DEALLOCATE(w%el_ptr, w%el, w%bp_ptr, w%bp)
    ALLOCATE(stamp(nprow/4), buf(nprow/4), w%el_ptr(cls%ncls+1), w%bp_ptr(cls%ncls+1))
    stamp = 0
    DO ip = 1, 2
      w%el_ptr(1) = 1
      w%bp_ptr(1) = 1
      DO k = 1, cls%ncls
        ne = 0
        DO j = cls%ptr(k), cls%ptr(k+1)-1
          i = cls%mem(j)
          DO m = mg_qlMat(lev)%LdA(i), mg_qlMat(lev)%LdA(i+1)-1
            CALL add_elem((mg_qlMat(lev)%ColA(m)-1)/4+1)
          END DO
          DO m = w%bt_ptr(i), w%bt_ptr(i+1)-1
            CALL add_elem((w%bt_prow(m)-1)/4+1)
          END DO
        END DO
        IF (ip.EQ.2) THEN
          w%el(w%el_ptr(k):w%el_ptr(k)+ne-1) = buf(1:ne)
          DO a = 1, ne
            DO b = 1, ne
              IF (a.EQ.b) THEN
                pos = 0
              ELSE
                pos = PCO_NEIGHBOUR_POS(buf(a), buf(b))
                IF (pos.EQ.0) THEN
                  w%nbad = w%nbad + 1
                  pos = -1
                END IF
              END IF
              w%bp(w%bp_ptr(k) + (a-1)*ne + b - 1) = pos
            END DO
          END DO
        END IF
        w%el_ptr(k+1) = w%el_ptr(k) + ne
        w%bp_ptr(k+1) = w%bp_ptr(k) + ne*ne
      END DO
      IF (ip.EQ.1) THEN
        ALLOCATE(w%el(MAX(w%el_ptr(cls%ncls+1)-1,1)), w%bp(MAX(w%bp_ptr(cls%ncls+1)-1,1)))
        w%nbad = 0
      END IF
    END DO
    DEALLOCATE(stamp, buf)
    w%maxne = 0
    DO k = 1, cls%ncls
      w%maxne = MAX(w%maxne, w%el_ptr(k+1) - w%el_ptr(k))
    END DO
    w%built = .TRUE.

  CONTAINS

    SUBROUTINE add_elem(e)
      INTEGER, INTENT(IN) :: e
      IF (stamp(e).EQ.(ip-1)*cls%ncls + k) RETURN
      stamp(e) = (ip-1)*cls%ncls + k
      ne = ne + 1
      buf(ne) = e
    END SUBROUTINE add_elem

  END SUBROUTINE PCO_WRAP_BUILD

  !-----------------------------------------------------------------------
  ! raw = (B^T Q) M_cls^-1 (Q^T B) on level w%lev (master; raw in the pco layout).
  ! mass_kind: PCO_MASS_ORDINARY (m = mass), PCO_MASS_EL_WEIGHTED (u, wgt raw masses),
  ! PCO_MASS_EL_PER_IMAGE (check variant). ku/kv/kw: velocity constraint flags of the level.
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_WRAP_ASSEMBLE(w, cls, mass_kind, mass, u, wgt, ku, kv, kw, raw)
    TYPE(tPcoWrap), INTENT(IN) :: w
    TYPE(tPcoClasses), INTENT(IN) :: cls
    INTEGER, INTENT(IN) :: mass_kind
    REAL*8, INTENT(IN) :: mass(:), u(:), wgt(:)
    INTEGER, INTENT(IN) :: ku(:), kv(:), kw(:)
    REAL*8, INTENT(INOUT) :: raw(:)
    INTEGER :: k, j, i, m, ne, a, b, mp, mq, e, mode, lev, slot, pos, r
    REAL*8, ALLOCATABLE :: bl(:,:,:), btl(:,:,:)
    REAL*8 :: s1, s2, mk, dm(3), val

    lev = w%lev
    raw = 0d0
    ALLOCATE(bl(3,4,MAX(w%maxne,1)), btl(3,4,MAX(w%maxne,1)))
    DO k = 1, cls%ncls
      ne = w%el_ptr(k+1) - w%el_ptr(k)
      IF (ne.EQ.0) CYCLE
      bl(:,:,1:ne) = 0d0
      btl(:,:,1:ne) = 0d0

      ! class mass
      s1 = 0d0
      s2 = 0d0
      DO j = cls%ptr(k), cls%ptr(k+1)-1
        i = cls%mem(j)
        SELECT CASE (mass_kind)
        CASE (PCO_MASS_ORDINARY)
          s1 = s1 + mass(i)
        CASE (PCO_MASS_EL_WEIGHTED)
          s1 = s1 + u(i)
          s2 = s2 + wgt(i)
        CASE (PCO_MASS_EL_PER_IMAGE)
          s1 = s1 + u(i)*u(i)/wgt(i)
        END SELECT
      END DO
      IF (mass_kind.EQ.PCO_MASS_EL_WEIGHTED) THEN
        mk = s1*s1/s2
      ELSE
        mk = s1
      END IF
      r = cls%mem(cls%ptr(k))
      dm(1) = MERGE(0d0, 1d0/mk, ku(r).EQ.1)
      dm(2) = MERGE(0d0, 1d0/mk, kv(r).EQ.1)
      dm(3) = MERGE(0d0, 1d0/mk, kw(r).EQ.1)

      ! class row of B and class column of B^T, merged per (element, mode)
      DO j = cls%ptr(k), cls%ptr(k+1)-1
        i = cls%mem(j)
        DO m = mg_qlMat(lev)%LdA(i), mg_qlMat(lev)%LdA(i+1)-1
          e = (mg_qlMat(lev)%ColA(m)-1)/4 + 1
          mode = mg_qlMat(lev)%ColA(m) - 4*(e-1)
          a = local(e)
          bl(1,mode,a) = bl(1,mode,a) + mg_BXMat(lev)%a(m)
          bl(2,mode,a) = bl(2,mode,a) + mg_BYMat(lev)%a(m)
          bl(3,mode,a) = bl(3,mode,a) + mg_BZMat(lev)%a(m)
        END DO
        DO m = w%bt_ptr(i), w%bt_ptr(i+1)-1
          slot = w%bt_pos(m)
          e = (w%bt_prow(m)-1)/4 + 1
          mode = w%bt_prow(m) - 4*(e-1)
          a = local(e)
          btl(1,mode,a) = btl(1,mode,a) + mg_BTXMat(lev)%a(slot)
          btl(2,mode,a) = btl(2,mode,a) + mg_BTYMat(lev)%a(slot)
          btl(3,mode,a) = btl(3,mode,a) + mg_BTZMat(lev)%a(slot)
        END DO
      END DO

      ! C(p,q) += sum_c BT_c(p) dm_c B_c(q)
      DO a = 1, ne
        DO b = 1, ne
          pos = w%bp(w%bp_ptr(k) + (a-1)*ne + b - 1)
          IF (pos.LT.0) CYCLE
          DO mp = 1, 4
            DO mq = 1, 4
              val = btl(1,mp,a)*dm(1)*bl(1,mq,b) + btl(2,mp,a)*dm(2)*bl(2,mq,b) + &
                    btl(3,mp,a)*dm(3)*bl(3,mq,b)
              slot = PCO_BLOCK_SLOT(w%el(w%el_ptr(k)+a-1), mp, pos, mq)
              raw(slot) = raw(slot) + val
            END DO
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(bl, btl)

  CONTAINS

    INTEGER FUNCTION local(e0)
      INTEGER, INTENT(IN) :: e0
      INTEGER :: t
      DO t = 1, ne
        IF (w%el(w%el_ptr(k)+t-1).EQ.e0) THEN
          local = t
          RETURN
        END IF
      END DO
      local = 0
    END FUNCTION local

  END SUBROUTINE PCO_WRAP_ASSEMBLE

END MODULE pco_wrap
