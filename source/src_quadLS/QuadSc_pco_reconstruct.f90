!=========================================================================
! QuadSc_pco_reconstruct.f90 - P0 -> vertex averaging of the type-3/4 coarse correction at
! periodic planes
!
! Part of the periodic coarse pressure operator (overview: QuadSc_pco_operator.f90).
!
! Coarse solver types 3 and 4 solve for the element means only, average them to the vertices
! (INTPVB in QuadSc_mg.f90: value of a vertex = sum_e (vol_e/8) p_e / sum_e (vol_e/8) over the
! elements of the vertex) and project the vertex field back to P1 (IntQ1toP1). On the
! unpartitioned mesh the two images of a periodic vertex see only the elements of
! their own side; the periodic average sums numerators AND denominators over the vertex
! class, divides once and gives the value to every image (the geometry of the projection
! stays the physical one of each element). For a vertex without images the formula is the
! one of INTPVB (equal to rounding; the compiler may contract the products differently).
!
! The classes are installed by the coordinator (PCO_SET_VERTEX_CLASSES); the shared type-3/4
! tail of the master coarse solve (crs_P0_Correction in QuadSc_mg.f90) calls PCO_P0_TO_VERTEX
! when a periodic operator is active and INTPVB otherwise. The vertex values and the
! denominators live in arrays of this module, sized nvt (the multigrid scratch vectors of the
! legacy tail are sized for the pressure vector, 4 nel, which a strip of elements with more
! vertices than 4 nel would overrun); pco_vertex_values holds the result.
!=========================================================================
MODULE pco_p0_reconstruct

  IMPLICIT NONE
  PRIVATE

  ! vertex classes with more than one member (CSR, ascending members)
  INTEGER, SAVE :: nvcls = 0
  INTEGER, ALLOCATABLE, SAVE :: vptr(:), vmem(:)
  ! result (vertex values) and denominators, nvt each
  REAL*8, ALLOCATABLE, SAVE, TARGET, PUBLIC :: pco_vertex_values(:)
  REAL*8, ALLOCATABLE, SAVE :: vden(:)

  PUBLIC :: PCO_SET_VERTEX_CLASSES, PCO_P0_TO_VERTEX

CONTAINS

  ! classes over the dofs (cls_of/ptr/mem of pco_periodic_topology); only the vertex classes
  ! (members <= nvt) with more than one member are kept
  SUBROUTINE PCO_SET_VERTEX_CLASSES(nvt, ncls, ptr, mem)
    INTEGER, INTENT(IN) :: nvt, ncls, ptr(:), mem(:)
    INTEGER :: k, n, m
    IF (ALLOCATED(vptr)) DEALLOCATE(vptr, vmem)
    n = 0
    m = 0
    DO k = 1, ncls
      IF (ptr(k+1)-ptr(k).LT.2 .OR. mem(ptr(k)).GT.nvt) CYCLE
      n = n + 1
      m = m + ptr(k+1) - ptr(k)
    END DO
    nvcls = n
    ALLOCATE(vptr(n+1), vmem(MAX(m,1)))
    vptr(1) = 1
    n = 0
    DO k = 1, ncls
      IF (ptr(k+1)-ptr(k).LT.2 .OR. mem(ptr(k)).GT.nvt) CYCLE
      n = n + 1
      vptr(n+1) = vptr(n) + ptr(k+1) - ptr(k)
      vmem(vptr(n):vptr(n+1)-1) = mem(ptr(k):ptr(k+1)-1)
    END DO
  END SUBROUTINE PCO_SET_VERTEX_CLASSES

  !-----------------------------------------------------------------------
  ! INTPVB with periodic vertex classes: pco_vertex_values(v), v = 1..nvt.
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_P0_TO_VERTEX(dp, dvol, kvert, nel, nvt)
    INTEGER, INTENT(IN) :: nel, nvt, kvert(8,*)
    REAL*8, INTENT(IN) :: dp(*), dvol(*)
    INTEGER :: iel, j, iv, k, m
    REAL*8 :: w, num, den

    IF (ALLOCATED(pco_vertex_values)) THEN
      IF (SIZE(pco_vertex_values).LT.nvt) DEALLOCATE(pco_vertex_values, vden)
    END IF
    IF (.NOT.ALLOCATED(pco_vertex_values)) ALLOCATE(pco_vertex_values(nvt), vden(nvt))
    ASSOCIATE(dpl => pco_vertex_values, daux => vden)

    dpl(1:nvt) = 0d0
    daux(1:nvt) = 0d0
    DO iel = 1, nel
      w = 0.125d0*dvol(iel)
      DO j = 1, 8
        iv = kvert(j,iel)
        dpl(iv) = dpl(iv) + w*dp(iel)
        daux(iv) = daux(iv) + w
      END DO
    END DO
    ! numerators and denominators summed over each class, the sums given to every image
    DO k = 1, nvcls
      num = 0d0
      den = 0d0
      DO m = vptr(k), vptr(k+1)-1
        num = num + dpl(vmem(m))
        den = den + daux(vmem(m))
      END DO
      DO m = vptr(k), vptr(k+1)-1
        dpl(vmem(m)) = num
        daux(vmem(m)) = den
      END DO
    END DO
    DO iv = 1, nvt
      dpl(iv) = dpl(iv)/daux(iv)
    END DO
    END ASSOCIATE
  END SUBROUTINE PCO_P0_TO_VERTEX

END MODULE pco_p0_reconstruct
