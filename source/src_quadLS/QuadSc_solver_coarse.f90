!===============================================================================
! Module: QuadSc_solver_coarse
!
! Purpose: UMFPACK coarse grid solver interface for QuadScalar
!          Extracted from QuadSc_def.f90 as part of Phase 1 refactoring
!          (See docs/md_docs/refactoring_roadmap.md)
!
! History:
!   Phase 1 (2025-12-04): Initial extraction from QuadSc_def.f90
!
! Notes:
!   - This module encapsulates coarse grid direct solver setup logic
!   - Supports UMFPACK types 2, 3, and 4
!   - Type 4 and 3 use /16 geometric coarsening
!   - Only executed on myid==0 (coarse grid process)
!===============================================================================

MODULE QuadSc_solver_coarse

USE PP3D_MPI, ONLY: myid, showID
USE var_QuadScalar
USE UMFPackSolver, ONLY: myUmfPack_Factorize

IMPLICIT NONE
PRIVATE

! Public interface
PUBLIC :: Setup_UMFPACK_CoarseSolver

CONTAINS

!===============================================================================
! Subroutine: Setup_UMFPACK_CoarseSolver
!
! Purpose: Set up UMFPACK direct solver for coarse grid pressure system
!
! Description:
!   - Handles three solver types:
!     * Type 2: Full coarse grid matrix (no geometric coarsening)
!     * Type 3: Coarse grid with /16 geometric coarsening (no factorization)
!     * Type 4: Coarse grid with /16 geometric coarsening + LU factorization
!   - Type 4/3 apply stride-4 extraction to reduce matrix size
!   - Only executes on myid==0 (coarse grid process)
!
! Arguments:
!   coarse_solver - INTEGER: solver type (1-4)
!   coarse_lev    - INTEGER: coarse level index (typically NLMIN)
!   bNoOutFlow    - LOGICAL: singular configuration flag
!
! Side Effects:
!   - Allocates and fills UMF_CMat, UMF_lMat (type 2)
!   - Allocates and fills crsSTR arrays (type 3/4)
!   - Calls myUmfPack_Factorize for LU decomposition (type 2/4)
!   - Prints diagnostic if at NLMIN
!===============================================================================
SUBROUTINE Setup_UMFPACK_CoarseSolver(knprP, coarse_solver, coarse_lev, bNoOutFlow)
  TYPE(mg_kVector), INTENT(IN) :: knprP(*)
  INTEGER, INTENT(IN) :: coarse_solver
  INTEGER, INTENT(IN) :: coarse_lev
  LOGICAL, INTENT(IN) :: bNoOutFlow

  INTEGER :: I, J, jCol, iEntry
  TYPE(TMatrix), POINTER :: lMat_work
  REAL*8, POINTER :: CMat_work(:)

  ! Only coarse grid process (myid==0) performs this setup
  IF (myid /= 0) RETURN

  ! Validate solver type
  IF (coarse_solver < 1 .OR. coarse_solver > 4) RETURN

  ! Set active level to coarse grid
  ILEV = coarse_lev
  CALL SETLEV(2)
  lMat_work => mg_lMat(ILEV)
  CMat_work => mg_CMat(ILEV)%a

  ! Handle singular (no outflow) boundary condition
  IF (bNoOutFlow) THEN
    DO I = lMat_work%LdA(1)+1, lMat_work%LdA(2)-1
      CMat_work(I) = 0d0
    END DO
  END IF

  ! ===========================================================================
  ! Solver Type 2: Full UMFPACK (no geometric coarsening)
  ! ===========================================================================
  IF (coarse_solver == 2) THEN
    ! Allocate UMFPACK arrays
    IF (.NOT. ALLOCATED(UMF_CMat)) ALLOCATE(UMF_CMat(lMat_work%na))
    UMF_CMat = CMat_work

    IF (.NOT. ALLOCATED(UMF_lMat%ColA)) ALLOCATE(UMF_lMat%ColA(lMat_work%na))
    IF (.NOT. ALLOCATED(UMF_lMat%LdA))  ALLOCATE(UMF_lMat%LdA(lMat_work%nu+1))

    ! Copy matrix structure
    UMF_lMat%ColA = lMat_work%ColA
    UMF_lMat%LdA  = lMat_work%LdA
    UMF_lMat%nu   = lMat_work%nu
    UMF_lMat%na   = lMat_work%na

    ! Perform LU factorization
    CALL myUmfPack_Factorize(UMF_CMat, UMF_lMat)
  END IF

  ! ===========================================================================
  ! Solver Type 4 & 3: UMFPACK with /16 geometric coarsening
  ! ===========================================================================
  IF (coarse_solver == 4 .OR. coarse_solver == 3) THEN
    ! Apply /4 reduction for unknowns, /16 for matrix entries
    crsSTR%A%nu = lMat_work%nu / 4
    crsSTR%A%na = lMat_work%na / 16

    ! Allocate coarse solver arrays
    IF (.NOT. ALLOCATED(crsSTR%A%LdA))  ALLOCATE(crsSTR%A%LdA(crsSTR%A%nu+1))
    IF (.NOT. ALLOCATED(crsSTR%A%ColA)) ALLOCATE(crsSTR%A%ColA(crsSTR%A%na))
    IF (.NOT. ALLOCATED(crsSTR%A_MAT))  ALLOCATE(crsSTR%A_MAT(crsSTR%A%na))
    IF (.NOT. ALLOCATED(crsSTR%A_RHS))  ALLOCATE(crsSTR%A_RHS(crsSTR%A%nu))
    IF (.NOT. ALLOCATED(crsSTR%A_SOL))  ALLOCATE(crsSTR%A_SOL(crsSTR%A%nu))

    ! Build coarse matrix structure (row pointers)
    crsSTR%A%LdA(1) = 1
    DO I = 1, crsSTR%A%nu
      J = 4*I - 3  ! Stride-4 access to fine grid
      crsSTR%A%LdA(I+1) = crsSTR%A%LdA(I) + (lMat_work%LdA(J+1) - lMat_work%LdA(J)) / 4
    END DO

    ! Extract coarse matrix entries (stride-4 access)
    iEntry = 0
    DO I = 1, crsSTR%A%nu
      J = 4*(I-1) + 1  ! Stride-4 access to fine grid
      DO jCol = lMat_work%LdA(J), lMat_work%LdA(J+1)-1, 4
        iEntry = iEntry + 1
        crsSTR%A%ColA(iEntry) = (lMat_work%ColA(jCol) - 1) / 4 + 1
        crsSTR%A_MAT(iEntry) = CMat_work(jCol)
      END DO
    END DO

    ! ========== Diagnostic Logging (Phase 0.4) ==========
    ! Verify /16 geometric coarsening efficiency
    IF (ILEV == NLMIN) THEN
      WRITE(MTERM,'(A)') ''
      WRITE(MTERM,'(A)') '========== UMFPACK Coarse Grid Diagnostic =========='
      WRITE(MTERM,'(A,I0)') 'Level:                            ', ILEV
      WRITE(MTERM,'(A,I0)') 'Fine grid unknowns (lMat%nu):     ', lMat_work%nu
      WRITE(MTERM,'(A,I0)') 'Fine grid non-zeros (lMat%na):    ', lMat_work%na
      WRITE(MTERM,'(A,I0)') 'Coarse unknowns (lMat%nu/4):      ', lMat_work%nu/4
      WRITE(MTERM,'(A,I0)') 'Coarse allocated (lMat%na/16):    ', lMat_work%na/16
      WRITE(MTERM,'(A,I0)') 'Coarse actual (iEntry):           ', iEntry
      IF (lMat_work%nu > 0) THEN
        WRITE(MTERM,'(A,F10.4)') 'nu ratio (fine/coarse):           ', &
                                  REAL(lMat_work%nu,8)/REAL(lMat_work%nu/4,8)
      END IF
      IF (iEntry > 0) THEN
        WRITE(MTERM,'(A,F10.4)') 'na ratio (fine/actual):           ', &
                                  REAL(lMat_work%na,8)/REAL(iEntry,8)
        WRITE(MTERM,'(A,F10.2,A)') 'Allocation efficiency:            ', &
                                     100.0*REAL(iEntry,8)/REAL(lMat_work%na/16,8), '%'
      END IF
      WRITE(MTERM,'(A)') '===================================================='
      WRITE(MTERM,'(A)') ''
    END IF
    ! =====================================================

    ! Apply boundary conditions (set diagonal to small value for constrained DOFs)
    DO I = 1, crsSTR%A%nu
      J = 4*(I-1) + 1
      IF (knprP(ILEV)%x(I) > 0) THEN
        crsSTR%A_MAT(crsSTR%A%LDA(I)) = 1d-16
      END IF
    END DO

    ! Perform LU factorization (type 4 only)
    IF (coarse_solver == 4) THEN
      CALL myUmfPack_Factorize(crsSTR%A_MAT, crsSTR%A)
    END IF
  END IF

END SUBROUTINE Setup_UMFPACK_CoarseSolver

END MODULE QuadSc_solver_coarse
