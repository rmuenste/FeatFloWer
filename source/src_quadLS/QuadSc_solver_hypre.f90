!===============================================================================
! Module: QuadSc_solver_hypre
!
! Purpose: HYPRE solver interface for QuadScalar
!          Extracted from QuadSc_def.f90 as part of Phase 1 refactoring
!          (See docs/md_docs/refactoring_roadmap.md)
!
! History:
!   Phase 1 (2025-12-04): Initial extraction from QuadSc_def.f90
!
! Notes:
!   - This module encapsulates HYPRE-specific matrix setup logic
!   - Two setup routines: Full coarse (type 7) and Geometric coarse (type 8)
!   - Both operate on MinLev (coarse level), then restore to NLMAX
!   - Uses global myHYPRE structure from types.f90
!===============================================================================

MODULE QuadSc_solver_hypre

USE PP3D_MPI, ONLY: myid, showID
USE var_QuadScalar
use, intrinsic :: ieee_arithmetic

IMPLICIT NONE
PRIVATE

! Public interface
PUBLIC :: Setup_HYPRE_CoarseLevel_Full
PUBLIC :: Setup_HYPRE_CoarseLevel_Geometric

CONTAINS

!===============================================================================
! Subroutine: Setup_HYPRE_CoarseLevel_Full
!
! Purpose: Set up HYPRE structures for coarse level (MinLev) WITHOUT geometric coarsening
!          (CrsSolverType == 7)
!
! Description:
!   - Creates global numbering for HYPRE on coarse level
!   - Fills HYPRE matrix structures with full CSR data (no /16 reduction)
!   - Handles both local (lMat, CMat) and parallel (lPMat, CPMat) contributions
!   - Converts to zero-based indexing if required
!   - Restores level to NLMAX at end
!
! Arguments: None (uses module-level pointers and global myHYPRE)
!
! Side Effects:
!   - Modifies global myHYPRE structure
!   - Allocates myHYPRE arrays (Numbering, ncols, sol, rhs, rows, cols, values)
!===============================================================================
SUBROUTINE Setup_HYPRE_CoarseLevel_Full(lScalar_in)
  TYPE(TLinScalar), INTENT(IN) :: lScalar_in
  INTEGER :: IEQ, IA, ICOL, II, III, NDOF_p, MaxDofs, NU, NEL
  INTEGER, ALLOCATABLE :: iDofs(:)
  REAL*8, ALLOCATABLE :: dDofs(:)

  IF (myid.ne.0) THEN
    ! Set active level to MinLev (coarsest)
    ILEV = lScalar_in%prm%MGprmIn%MinLev
    CALL SETLEV(2)
    lMat      => mg_lMat(ILEV)
    CMat      => mg_CMat(ILEV)%a
    lPMat     => mg_lPMat(ILEV)
    CPMat     => mg_CPMat(ILEV)%a

    ! Allocate numbering array
    myHYPRE%nrows = lPMat%nu
    allocate(myHYPRE%Numbering(myHYPRE%nrows))
  END IF

  ! Get global numbering limits for this process
  CALL GetMyHYPRENumberingLimits(myHYPRE%ilower, myHYPRE%iupper, NEL)

  IF (myid.ne.0) THEN
    ! Create consecutive global numbering starting from ilower
    DO IEQ = 1, myHYPRE%nrows
      myHYPRE%Numbering(IEQ) = myHYPRE%ilower + IEQ - 1
    END DO
  END IF

  ! Fill up the HYPRE structures
  IF (myid.ne.0) THEN
    ! Determine maximum column index from parallel matrix
    NDOF_p = 0
    DO IEQ=1, lPMat%nu
      DO IA=lPMat%LdA(IEQ), lPMat%LdA(IEQ+1)-1
        ICOL = lPMat%ColA(IA)
        NDOF_p = max(ICOL, NDOF_p)
      END DO
    END DO

    ! Allocate and get off-partition numbering
    IF (ALLOCATED(myHYPRE%OffPartitionNumbering)) DEALLOCATE(myHYPRE%OffPartitionNumbering)
    ALLOCATE(myHYPRE%OffPartitionNumbering(NDOF_p))
    CALL GetHYPREParPressureIndices(myHYPRE%OffPartitionNumbering)

    ! Allocate main HYPRE arrays
    MaxDofs = 0
    myHYPRE%nonzeros = lPMat%na + lMat%na
    allocate(myHYPRE%ncols(myHYPRE%nrows))
    allocate(myHYPRE%sol(myHYPRE%nrows))
    allocate(myHYPRE%rhs(myHYPRE%nrows))
    allocate(myHYPRE%rows(myHYPRE%nrows))
    allocate(myHYPRE%cols(myHYPRE%nonzeros))
    allocate(myHYPRE%values(myHYPRE%nonzeros))

    ! Count columns per row
    DO IEQ=1, myHYPRE%nrows
      NU = (lMat%LdA(IEQ+1) - lMat%LdA(IEQ)) + (lPMat%LdA(IEQ+1) - lPMat%LdA(IEQ))
      myHYPRE%ncols(IEQ) = NU
      MaxDofs = max(MaxDofs, NU)
    END DO

    ! Set up row indices
    DO IEQ=1, myHYPRE%nrows
      myHYPRE%rows(IEQ) = myHYPRE%Numbering(IEQ)
    END DO

    ! Allocate work arrays
    allocate(iDofs(MaxDofs), dDofs(MaxDofs))

    ! Assemble column indices and values
    III = 0
    DO IEQ=1, myHYPRE%nrows
      II = 0
      ! Local matrix contributions
      DO IA=lMat%LdA(IEQ), lMat%LdA(IEQ+1)-1
        ICOL = lMat%ColA(IA)
        II = II + 1
        iDofs(II) = myHYPRE%Numbering(ICOL)
        dDofs(II) = CMat(IA)
      END DO
      ! Parallel matrix contributions
      DO IA=lPMat%LdA(IEQ), lPMat%LdA(IEQ+1)-1
        ICOL = lPMat%ColA(IA)
        II = II + 1
        iDofs(II) = myHYPRE%OffPartitionNumbering(ICOL)
        dDofs(II) = CPMat(IA)
      END DO

      ! Sort DOFs by column index
      CALL SORT_DOFs(iDofs, dDofs, II)

      ! Copy sorted data into HYPRE arrays
      DO IA=1, II
        III = III + 1
        myHYPRE%cols(III)   = iDofs(IA)
        myHYPRE%values(III) = dDofs(IA)
      END DO
    END DO

    ! Convert to zero-based indexing if needed
    IF (myHYPRE%ZeroBased) THEN
      myHYPRE%ilower = myHYPRE%ilower - 1
      myHYPRE%iupper = myHYPRE%iupper - 1
      myHYPRE%rows = myHYPRE%rows - 1
      myHYPRE%cols = myHYPRE%cols - 1
    END IF

    ! Clean up work arrays
    deallocate(iDofs, dDofs)

    ! Restore level to NLMAX
    ILEV = NLMAX
    CALL SETLEV(2)
    lMat      => mg_lMat(ILEV)
    CMat      => mg_CMat(ILEV)%a
    lPMat     => mg_lPMat(ILEV)
    CPMat     => mg_CPMat(ILEV)%a
  END IF

END SUBROUTINE Setup_HYPRE_CoarseLevel_Full

!===============================================================================
! Subroutine: Setup_HYPRE_CoarseLevel_Geometric
!
! Purpose: Set up HYPRE structures for coarse level (MinLev) WITH /16 geometric coarsening
!          (CrsSolverType == 8)
!
! Description:
!   - Creates global numbering for HYPRE on coarse level
!   - Applies /4 reduction for unknowns and /16 for matrix entries
!   - Stride-4 loops extract coarse grid stencil
!   - Fills HYPRE matrix structures with CSR data
!   - Optional diagnostic output for verification
!   - Restores level to NLMAX at end
!
! Arguments: None (uses module-level pointers and global myHYPRE)
!
! Side Effects:
!   - Modifies global myHYPRE structure
!   - Allocates myHYPRE arrays (Numbering, ncols, sol, rhs, rows, cols, values)
!   - Prints diagnostic if at NLMIN and myid == showID
!===============================================================================
SUBROUTINE Setup_HYPRE_CoarseLevel_Geometric(lScalar_in)
  TYPE(TLinScalar), INTENT(IN) :: lScalar_in
  INTEGER :: IEQ, IA, ICOL, II, III, NDOF_p, MaxDofs, NU, NEL, JEQ
  INTEGER, ALLOCATABLE :: iDofs(:)
  REAL*8, ALLOCATABLE :: dDofs(:)
  LOGICAL :: bNoOutFlow

  IF (myid.ne.0) THEN
    ! Set active level to MinLev (coarsest)
    ILEV = lScalar_in%prm%MGprmIn%MinLev
    CALL SETLEV(2)

    lMat      => mg_lMat(ILEV)
    CMat      => mg_CMat(ILEV)%a
    lPMat     => mg_lPMat(ILEV)
    CPMat     => mg_CPMat(ILEV)%a

    ! Allocate numbering array
    myHYPRE%nrows = lPMat%nu
    allocate(myHYPRE%Numbering(myHYPRE%nrows))
  END IF

  ! Handle singular (no outflow) configuration
  bNoOutFlow = .false.  ! TODO: Get from appropriate module/parameter
  IF (myid.ne.0 .AND. bNoOutFlow) THEN
    IF (myid.eq.1) THEN
      WRITE(*,*) 'Imposing Dirichlet pressure for the singular (no outflow) configuration'
      DO IEQ=lMat%LdA(1)+1, lMat%LdA(2)-1
        CMat(IEQ) = 0d0
      END DO
    END IF
  END IF

  ! Get global numbering limits
  CALL GetMyHYPRENumberingLimits(myHYPRE%ilower, myHYPRE%iupper, NEL)

  IF (myid.ne.0) THEN
    ! Create consecutive global numbering
    DO IEQ = 1, myHYPRE%nrows
      myHYPRE%Numbering(IEQ) = myHYPRE%ilower + IEQ - 1
    END DO
  END IF

  ! Fill up the HYPRE structures with /16 coarsening
  IF (myid.ne.0) THEN
    ! Determine maximum column index from parallel matrix
    NDOF_p = 0
    DO IEQ=1, lPMat%nu
      DO IA=lPMat%LdA(IEQ), lPMat%LdA(IEQ+1)-1
        ICOL = lPMat%ColA(IA)
        NDOF_p = max(ICOL, NDOF_p)
      END DO
    END DO

    ! Allocate and get off-partition numbering
    IF (ALLOCATED(myHYPRE%OffPartitionNumbering)) DEALLOCATE(myHYPRE%OffPartitionNumbering)
    ALLOCATE(myHYPRE%OffPartitionNumbering(NDOF_p))
    CALL GetHYPREParPressureIndices(myHYPRE%OffPartitionNumbering)

    ! Apply /4 reduction for unknowns
    MaxDofs = 0
    myHYPRE%nrows = lPMat%nu / 4
    myHYPRE%ilower = (myHYPRE%ilower + 3) / 4
    myHYPRE%iupper = myHYPRE%iupper / 4

    ! Apply /16 reduction for matrix entries
    myHYPRE%nonzeros = lPMat%na/16 + lMat%na/16
    allocate(myHYPRE%ncols(myHYPRE%nrows))
    allocate(myHYPRE%sol(myHYPRE%nrows))
    allocate(myHYPRE%rhs(myHYPRE%nrows))
    allocate(myHYPRE%rows(myHYPRE%nrows))
    allocate(myHYPRE%cols(myHYPRE%nonzeros))
    allocate(myHYPRE%values(myHYPRE%nonzeros))

    ! Count columns per row (stride-4 access)
    DO IEQ=1, myHYPRE%nrows
      JEQ = 4*(IEQ-1) + 1
      NU = ((lMat%LdA(JEQ+1) - lMat%LdA(JEQ)) + (lPMat%LdA(JEQ+1) - lPMat%LdA(JEQ))) / 4
      myHYPRE%ncols(IEQ) = NU
      MaxDofs = max(MaxDofs, NU)
    END DO

    ! Set up row indices (stride-4 access)
    DO IEQ=1, myHYPRE%nrows
      JEQ = 4*(IEQ-1) + 1
      myHYPRE%rows(IEQ) = (myHYPRE%Numbering(JEQ) + 3) / 4
    END DO

    ! Allocate work arrays
    allocate(iDofs(MaxDofs), dDofs(MaxDofs))

    ! Assemble column indices and values (stride-4 access)
    III = 0
    DO IEQ=1, myHYPRE%nrows
      II = 0
      JEQ = 4*(IEQ-1) + 1
      ! Local matrix contributions (stride-4)
      DO IA=lMat%LdA(JEQ), lMat%LdA(JEQ+1)-1, 4
        ICOL = lMat%ColA(IA)
        II = II + 1
        iDofs(II) = (myHYPRE%Numbering(ICOL) + 3) / 4
        dDofs(II) = CMat(IA)
      END DO
      ! Parallel matrix contributions (stride-4)
      DO IA=lPMat%LdA(JEQ), lPMat%LdA(JEQ+1)-1, 4
        ICOL = lPMat%ColA(IA)
        II = II + 1
        iDofs(II) = (myHYPRE%OffPartitionNumbering(ICOL) + 3) / 4
        dDofs(II) = CPMat(IA)
      END DO

      ! Sort DOFs by column index
      CALL SORT_DOFs(iDofs, dDofs, II)

      ! Copy sorted data into HYPRE arrays
      DO IA=1, II
        III = III + 1
        myHYPRE%cols(III)   = iDofs(IA)
        myHYPRE%values(III) = dDofs(IA)
      END DO
    END DO

    ! ========== Diagnostic Logging (Phase 0.4) ==========
    ! Verify /16 geometric coarsening efficiency
    IF (myid == showID .and. ILEV == NLMIN) THEN
      WRITE(MTERM,'(A)') ''
      WRITE(MTERM,'(A)') '========== HYPRE Coarse Grid Diagnostic =========='
      WRITE(MTERM,'(A,I0)') 'Level:                            ', ILEV
      WRITE(MTERM,'(A,I0)') 'lMat%nu:                          ', lMat%nu
      WRITE(MTERM,'(A,I0)') 'lMat%na:                          ', lMat%na
      WRITE(MTERM,'(A,I0)') 'lPMat%nu:                         ', lPMat%nu
      WRITE(MTERM,'(A,I0)') 'lPMat%na:                         ', lPMat%na
      WRITE(MTERM,'(A,I0)') 'HYPRE nrows (lPMat%nu/4):         ', myHYPRE%nrows
      WRITE(MTERM,'(A,I0)') 'HYPRE allocated (lPMat%na/16 + lMat%na/16): ', myHYPRE%nonzeros
      WRITE(MTERM,'(A,I0)') '  = lPMat%na/16:                  ', lPMat%na/16
      WRITE(MTERM,'(A,I0)') '  + lMat%na/16:                   ', lMat%na/16
      WRITE(MTERM,'(A,I0)') 'HYPRE actual (III):               ', III
      IF (III > 0) THEN
        WRITE(MTERM,'(A,F10.2,A)') 'Allocation efficiency:            ', &
                                     100.0*REAL(III,8)/REAL(myHYPRE%nonzeros,8), '%'
      END IF
      WRITE(MTERM,'(A)') '===================================================='
      WRITE(MTERM,'(A)') ''
    END IF
    ! =====================================================

    ! Convert to zero-based indexing if needed
    IF (myHYPRE%ZeroBased) THEN
      myHYPRE%ilower = myHYPRE%ilower - 1
      myHYPRE%iupper = myHYPRE%iupper - 1
      myHYPRE%rows = myHYPRE%rows - 1
      myHYPRE%cols = myHYPRE%cols - 1
    END IF

    ! Clean up work arrays
    deallocate(iDofs, dDofs)

    ! Restore level to NLMAX
    ILEV = NLMAX
    CALL SETLEV(2)
    lMat      => mg_lMat(ILEV)
    CMat      => mg_CMat(ILEV)%a
    lPMat     => mg_lPMat(ILEV)
    CPMat     => mg_CPMat(ILEV)%a
  END IF

END SUBROUTINE Setup_HYPRE_CoarseLevel_Geometric

!===============================================================================
! Subroutine: SORT_DOFs
!
! Purpose: Sort DOF column indices and corresponding matrix values
!
! Description:
!   Simple insertion sort for small arrays
!   Sorts iDofs (column indices) in ascending order
!   Maintains correspondence with dDofs (matrix values)
!
! Arguments:
!   LW(N) - INTEGER array of column indices (in/out)
!   RW(N) - REAL*8 array of matrix values (in/out)
!   N     - Number of entries to sort
!===============================================================================
SUBROUTINE SORT_DOFs(LW, RW, N)
  INTEGER, INTENT(IN) :: N
  INTEGER :: LW(N), LWA
  REAL*8  :: RW(N), RWA
  INTEGER :: I, J

  ! Insertion sort
  DO I = 2, N
    LWA = LW(I)
    RWA = RW(I)
    J = I - 1
    DO WHILE (J >= 1)
      IF (LW(J) <= LWA) EXIT
      LW(J+1) = LW(J)
      RW(J+1) = RW(J)
      J = J - 1
    END DO
    LW(J+1) = LWA
    RW(J+1) = RWA
  END DO

END SUBROUTINE SORT_DOFs

END MODULE QuadSc_solver_hypre
