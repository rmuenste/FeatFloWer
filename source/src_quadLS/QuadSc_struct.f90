!===============================================================================
! Module: QuadSc_struct
!
! Purpose: Matrix structure allocation routines for QuadScalar solver
!          Extracted from QuadSc_def.f90 as part of Phase 3.1 refactoring
!          (See docs/md_docs/refactoring_roadmap.md)
!
! History:
!   Phase 3.1 (2025-12-06): Initial extraction from QuadSc_def.f90
!
! Notes:
!   - This module contains pure structure allocation routines (no coefficients)
!   - Creates sparsity patterns for multigrid levels (NLMIN:NLMAX)
!   - Used during initialization before coefficient assembly
!
! Routines:
!   - Create_QuadMatStruct: Q2 velocity matrix structure (M, K, D, A matrices)
!   - Create_QuadLinMatStruct: Q2-P1 coupling matrices (B, B^T)
!   - Create_LinMatStruct: P1 pressure matrix structure (C matrix)
!   - Create_ParLinMatStruct: Parallel P1 pressure matrix structure
!===============================================================================

MODULE QuadSc_struct

USE PP3D_MPI, ONLY: myid, showID
USE var_QuadScalar

IMPLICIT NONE
PRIVATE

! Public interface
PUBLIC :: Create_QuadMatStruct
PUBLIC :: Create_QuadLinMatStruct
PUBLIC :: Create_LinMatStruct
PUBLIC :: Create_ParLinMatStruct

CONTAINS

!===============================================================================
! Subroutine: Create_QuadMatStruct
!
! Purpose: Create sparsity pattern for Q2 velocity-space matrices
!
! Description:
!   Allocates and initializes the matrix structure descriptor mg_qMat for all
!   multigrid levels. This structure is shared by mass matrix (M), stiffness
!   matrix (K), diffusion matrix (D), and advection matrix (A).
!
!   Q2 elements: Quadratic basis functions (27 DOFs per hexahedron)
!
! Side Effects:
!   - Allocates mg_qMat(NLMIN:NLMAX)
!   - For each level:
!     * Allocates mg_qMat(ILEV)%LdA (row pointers)
!     * Allocates mg_qMat(ILEV)%ColA (column indices)
!     * Sets mg_qMat(ILEV)%nu (number of unknowns)
!     * Sets mg_qMat(ILEV)%na (number of non-zeros)
!   - Sets global pointer qMat => mg_qMat(NLMAX)
!
! Performance: O(NDOF) per level, NDOF ~ vertices + edges + faces + elements
!===============================================================================
SUBROUTINE Create_QuadMatStruct()
INTEGER iSymm,nERow
INTEGER , DIMENSION(:)  , ALLOCATABLE :: TempColA
INTEGER I,J,MatSize,NDOF
REAL(KIND=8) :: t_start, t_end  ! Phase 0.5: Performance timing
EXTERNAL E013,coefst

! ============== Performance Timing (Phase 0.5) ==============
CALL CPU_TIME(t_start)
! =============================================================

 ALLOCATE(mg_qMat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX
   CALL SETLEV(2)

   ndof = mg_mesh%level(ilev)%nvt+&
          mg_mesh%level(ilev)%net+&
          mg_mesh%level(ilev)%nat+&
          mg_mesh%level(ilev)%nel

   MatSize = 300*NDOF

   ALLOCATE(TempColA(MatSize))
   ALLOCATE(mg_qMat(ILEV)%LdA(NDOF+1))
   mg_qMat(ILEV)%nu = NDOF
   mg_qMat(ILEV)%na = MatSize
   iSymm =   0
   nERow =   300

   CALL AP7(TempColA,mg_qMat(ILEV)%LdA,mg_qMat(ILEV)%na,&
            mg_qMat(ILEV)%nu,E013,iSymm,nERow,&
                      mg_mesh%level(ILEV)%kvert,&
                      mg_mesh%level(ILEV)%kedge,&
                      mg_mesh%level(ILEV)%karea)

   IF (myid.eq.showID) WRITE(MTERM,'(A40,2I10)') &
      "M,K,D,A matrix structure created",mg_qMat(ILEV)%nu,mg_qMat(ILEV)%na

   ALLOCATE(mg_qMat(ILEV)%ColA(mg_qMat(ILEV)%na))
   mg_qMat(ILEV)%ColA(:) = TempColA(1:mg_qMat(ILEV)%na)
   DEALLOCATE(TempColA)

 END DO

 ILEV=NLMAX
 CALL SETLEV(2)
 qMat => mg_qMat(NLMAX)

! ============== Performance Timing (Phase 0.5) ==============
CALL CPU_TIME(t_end)
IF (myid == showID) THEN
  WRITE(MTERM,'(A,F10.4,A)') 'Create_QuadMatStruct time:        ', t_end - t_start, ' sec'
END IF
! =============================================================

END SUBROUTINE Create_QuadMatStruct
!
! ----------------------------------------------
!
!===============================================================================
! Subroutine: Create_QuadLinMatStruct
!
! Purpose: Create sparsity patterns for Q2-P1 coupling matrices (gradient operators)
!
! Description:
!   Allocates and initializes two matrix structure descriptors:
!   - mg_qlMat: Q2 → P1 coupling (gradient B: velocity → pressure)
!   - mg_lqMat: P1 → Q2 coupling (divergence B^T: pressure → velocity)
!
!   These matrices couple the quadratic velocity space to the linear pressure
!   space in the Q2/P1 Stokes discretization.
!
! Side Effects:
!   - Allocates mg_qlMat(NLMIN:NLMAX), mg_lqMat(NLMIN:NLMAX)
!   - For each level:
!     * mg_qlMat: nu ~ Q2 DOFs, structure for B matrices (Bx, By, Bz)
!     * mg_lqMat: nu ~ 4*nel (P1 DOFs), structure for B^T matrices
!   - Sets global pointers qlMat => mg_qlMat(NLMAX), lqMat => mg_lqMat(NLMAX)
!
! Performance: O(nel) per level, nel = number of elements
!===============================================================================
SUBROUTINE Create_QuadLinMatStruct()
INTEGER iSymm,nERow
INTEGER , DIMENSION(:)  , ALLOCATABLE :: TempColA,TempLdA
INTEGER I,J,MatSize,NDOF
REAL(KIND=8) :: t_start, t_end  ! Phase 0.5: Performance timing
CHARACTER*10 myFile
EXTERNAL E011,E013,E010,coefst

! ============== Performance Timing (Phase 0.5) ==============
CALL CPU_TIME(t_start)
! =============================================================

 ALLOCATE (mg_qlMat(NLMIN:NLMAX))
 ALLOCATE (mg_lqMat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX
  CALL SETLEV(2)

  NDOF = mg_mesh%level(ilev)%nvt+&
         mg_mesh%level(ilev)%net+&
         mg_mesh%level(ilev)%nat+&
         mg_mesh%level(ilev)%nel

  MatSize = 16*27*mg_mesh%level(ilev)%nel

  ALLOCATE(TempColA(MatSize))
  ALLOCATE(mg_qlMat(ILEV)%LdA(NDOF+1))
  mg_qlMat(ILEV)%nu = NDOF
  mg_qlMat(ILEV)%na = MatSize
  iSymm =   0
  nERow =   16

  CALL AP9(TempColA,mg_qlMat(ILEV)%LdA,mg_qlMat(ILEV)%na,&
           mg_qlMat(ILEV)%nu,E013,E010,nERow,&
           mg_mesh%level(ILEV)%kvert,&
           mg_mesh%level(ILEV)%kedge,&
           mg_mesh%level(ILEV)%karea)

  mg_qlMat(ILEV)%na = 4*mg_qlMat(ILEV)%na

  ALLOCATE(mg_qlMat(ILEV)%ColA(mg_qlMat(ILEV)%na))

  CALL MatStructQ2P1(TempColA,mg_qlMat(ILEV)%ColA,mg_qlMat(ILEV)%LdA,&
       mg_qlMat(ILEV)%na,mg_qlMat(ILEV)%nu)

  IF (myid.eq.showID) WRITE(MTERM,'(A40,2I10)') &
   "B matrix structure created",mg_qlMat(ILEV)%nu,mg_qlMat(ILEV)%na

!   CALL OutputMatrixStuct("MatB",mg_qlMat(ILEV))

  MatSize = 4*27*mg_mesh%level(ilev)%nel
  ALLOCATE (TempLdA(4*mg_mesh%level(ilev)%nel+1))
  mg_lqMat(ILEV)%nu = mg_mesh%level(ilev)%nel
  mg_lqMat(ILEV)%na = MatSize
  iSymm =   0
  nERow =   27

  CALL AP9(TempColA,TempLdA,mg_lqMat(ILEV)%na,mg_lqMat(ILEV)%nu,E010,E013,nERow,&
           mg_mesh%level(ILEV)%kvert,&
           mg_mesh%level(ILEV)%kedge,&
           mg_mesh%level(ILEV)%karea)

  mg_lqMat(ILEV)%nu = 4*mg_lqMat(ILEV)%nu
  mg_lqMat(ILEV)%na = 4*mg_lqMat(ILEV)%na

  ALLOCATE(mg_lqMat(ILEV)%LdA(mg_lqMat(ILEV)%nu+1))
  ALLOCATE(mg_lqMat(ILEV)%ColA(mg_lqMat(ILEV)%na))

  CALL MatStructP1Q2(TempLdA,mg_lqMat(ILEV)%LdA,&
                     TempColA,&
                     mg_lqMat(ILEV)%ColA,&
                     MatSize,mg_mesh%level(ilev)%nel)

  IF (myid.eq.showID) WRITE(MTERM,'(A40,2I10)') &
  "BT matrix structure created",mg_lqMat(ILEV)%nu,mg_lqMat(ILEV)%na

  DEALLOCATE(TempColA,TempLdA)

 END DO

 ILEV=NLMAX
 CALL SETLEV(2)

 qlMat => mg_qlMat(NLMAX)
 lqMat => mg_lqMat(NLMAX)

! ============== Performance Timing (Phase 0.5) ==============
CALL CPU_TIME(t_end)
IF (myid == showID) THEN
  WRITE(MTERM,'(A,F10.4,A)') 'Create_QuadLinMatStruct time:     ', t_end - t_start, ' sec'
END IF
! =============================================================

END SUBROUTINE Create_QuadLinMatStruct
!
! ----------------------------------------------
!
!===============================================================================
! Subroutine: Create_LinMatStruct
!
! Purpose: Create sparsity pattern for P1 pressure matrix (local part)
!
! Description:
!   Allocates and initializes the matrix structure descriptor mg_lMat for the
!   P1 pressure Laplacian (C matrix). This is the local (on-process) part of
!   the pressure matrix used in the Schur complement.
!
!   P1 elements: Linear basis functions (4 DOFs per hexahedron)
!
! Side Effects:
!   - Allocates mg_lMat(NLMIN:NLMAX)
!   - For each level:
!     * Calls Get_CMatLen to determine matrix size
!     * Calls Get_CMatStruct to build sparsity pattern
!     * Allocates mg_lMat(ILEV)%LdA, mg_lMat(ILEV)%ColA
!   - Sets global pointer lMat => mg_lMat(NLMAX)
!
! Dependencies:
!   - Requires mg_qlMat already allocated (created by Create_QuadLinMatStruct)
!   - Uses Get_CMatLen, Get_CMatStruct assembly kernels
!===============================================================================
SUBROUTINE Create_LinMatStruct()

 ALLOCATE (mg_lMat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)
  qlMat => mg_qlMat(ILEV)
  qlPMat => mg_qlPMat(ILEV)

  CALL Get_CMatLen(qlMat%LdA,qlMat%ColA,&
                   mg_mesh%level(ILEV)%kvert,&
                   mg_mesh%level(ilev)%nel,&
                   mg_lMat(ILEV)%na,1)



  mg_lMat(ILEV)%nu = 4*mg_mesh%level(ilev)%nel
  ALLOCATE (mg_lMat(ILEV)%LdA(mg_lMat(ILEV)%nu+1),mg_lMat(ILEV)%ColA(mg_lMat(ILEV)%na))

  mg_lMat(ILEV)%LdA=0
  mg_lMat(ILEV)%ColA=0

  CALL Get_CMatStruct(mg_lMat(ILEV)%LdA,mg_lMat(ILEV)%ColA,qlMat%LdA,qlMat%ColA,&
       mg_mesh%level(ilev)%kvert,mg_mesh%level(ilev)%nel,mg_lMat(ILEV)%na,1)

  ! CALL OutputMatrixStuct("MatC",mg_lMat(ILEV))
  IF (myid.eq.showID) WRITE(MTERM,'(A40,2I10)') &
  "C matrix structure created",mg_lMat(ILEV)%nu,mg_lMat(ILEV)%na

 END DO

 ILEV=NLMAX
 CALL SETLEV(2)

 qlMat => mg_qlMat(NLMAX)
 lMat  => mg_lMat(NLMAX)

END SUBROUTINE Create_LinMatStruct
!
! ----------------------------------------------
!
!===============================================================================
! Subroutine: Create_ParLinMatStruct
!
! Purpose: Create sparsity pattern for P1 pressure matrix (parallel part)
!
! Description:
!   Allocates and initializes the matrix structure descriptor mg_lPMat for the
!   parallel (off-process) part of the P1 pressure Laplacian. This handles
!   contributions from DOFs owned by neighboring MPI processes.
!
!   Used in conjunction with mg_lMat for complete pressure system assembly.
!
! Side Effects:
!   - Allocates mg_lPMat(NLMIN:NLMAX)
!   - For each level:
!     * Calls Get_CMatLen with mode=2 (parallel mode)
!     * Calls Get_CMatStruct to build parallel sparsity pattern
!     * Allocates mg_lPMat(ILEV)%LdA, mg_lPMat(ILEV)%ColA
!   - Sets global pointers qlPMat, lPMat => NLMAX level
!
! Dependencies:
!   - Requires mg_qlPMat already allocated (parallel B matrix structure)
!   - Uses Get_CMatLen, Get_CMatStruct with parallel mode
!
! Notes:
!   - Only relevant for multi-process runs (myid != 0)
!   - Line 256: Bug? Sets qlPMat => mg_lqMat (should be mg_qlPMat?)
!===============================================================================
SUBROUTINE Create_ParLinMatStruct()

 ALLOCATE (mg_lPMat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)
  qlPMat => mg_qlPMat(ILEV)




  CALL Get_CMatLen(qlPMat%LdA,qlPMat%ColA,&
                   mg_mesh%level(ilev)%kvert,&
                   mg_mesh%level(ilev)%nel,&
                   mg_lPMat(ILEV)%na,2)

  mg_lPMat(ILEV)%nu = 4*NEL
  ALLOCATE (mg_lPMat(ILEV)%LdA(mg_lPMat(ILEV)%nu+1),mg_lPMat(ILEV)%ColA(mg_lPMat(ILEV)%na))
  mg_lPMat(ILEV)%LdA=0
  mg_lPMat(ILEV)%ColA=0

  CALL Get_CMatStruct(mg_lPMat(ILEV)%LdA,mg_lPMat(ILEV)%ColA,qlPMat%LdA,qlPMat%ColA,&
       mg_mesh%level(ilev)%kvert,mg_mesh%level(ilev)%nel,&
       mg_lPMat(ILEV)%na,2)

  ! CALL OutputMatrixStuct("MaPC",mg_lPMat(ILEV))
  IF (myid.eq.showID) WRITE(MTERM,'(A40,2I10)') &
  "Parallel C matrix structure created",mg_lPMat(ILEV)%nu,mg_lPMat(ILEV)%na
  ! pause
 END DO

 ILEV=NLMAX
 CALL SETLEV(2)

 qlPMat => mg_lqMat(NLMAX)
 lPMat => mg_lPMat(NLMAX)

END SUBROUTINE Create_ParLinMatStruct
!
! ----------------------------------------------
!
!===============================================================================
! Subroutine: MatStructQ2P1
!
! Purpose: Build column indices for Q2-P1 coupling matrix structure
!
! Description:
!   Helper routine for Create_QuadLinMatStruct. Transforms the sparsity
!   pattern from single-component to 4-component (for x,y,z velocity + scalar).
!
!   Takes the preliminary column structure (ColAo) and expands it to account
!   for the 4 components of the coupled system.
!
! Arguments:
!   ColAo - INTEGER(*): Input column indices (preliminary pattern)
!   ColAn - INTEGER(*): Output column indices (expanded pattern)
!   LdA   - INTEGER(*): Row pointers (modified in-place)
!   na    - INTEGER: Number of non-zeros
!   nu    - INTEGER: Number of unknowns
!
! Side Effects:
!   - Modifies LdA array in-place (expands row pointers)
!   - Fills ColAn with expanded column indices
!===============================================================================
SUBROUTINE MatStructQ2P1(ColAo,ColAn,LdA,na,nu)
IMPLICIT NONE
INTEGER ColAo(*),ColAn(*),LdA(*),na,nu
INTEGER i,j,k,nn,ijEntry

nn = 0
DO i=1,nu
 DO j=LdA(i),LdA(i+1)-1
  ijEntry = ColAo(j)
  DO k=1,4
   nn = nn + 1
   ColAn(nn) = 4*(ijEntry-1) + k
  END DO
 END DO
END DO

DO i=2,nu+1
 LdA(i) = 4*(LdA(i)-1)+1
END DO

END SUBROUTINE MatStructQ2P1
!
! ----------------------------------------------
!
!===============================================================================
! Subroutine: MatStructP1Q2
!
! Purpose: Build structure for P1-Q2 coupling matrix (transpose pattern)
!
! Description:
!   Helper routine for Create_QuadLinMatStruct. Builds the transpose pattern
!   (B^T) from the preliminary Q2-P1 pattern.
!
!   Expands the pattern to 4 components and adjusts row pointers for the
!   transposed structure.
!
! Arguments:
!   LdAo  - INTEGER(*): Input row pointers (preliminary)
!   LdAn  - INTEGER(*): Output row pointers (expanded)
!   ColAo - INTEGER(*): Input column indices
!   ColAn - INTEGER(*): Output column indices (expanded)
!   na    - INTEGER: Number of non-zeros
!   nu    - INTEGER: Number of unknowns
!
! Side Effects:
!   - Fills LdAn with expanded row pointers
!   - Fills ColAn with expanded column indices
!===============================================================================
SUBROUTINE MatStructP1Q2(LdAo,LdAn,ColAo,ColAn,na,nu)
IMPLICIT NONE
INTEGER ColAo(*),ColAn(*),LdAo(*),LdAn(*),na,nu
INTEGER i,j,k,nn,iPos

nn = 0
DO i=1,nu
 DO k=1,4
  DO j=LdAo(i),LdAo(i+1)-1
   nn = nn + 1
   ColAn(nn) = ColAo(j)
  END DO
 END DO
END DO

nn = 0
iPos = -26
DO i=1,nu
 DO k=1,4
  nn = nn + 1
  iPos = iPos + 27
  LdAn(nn) = iPos
 END DO
END DO
nn = nn + 1
LdAn(nn) = LdAn(nn-1) + 27

END SUBROUTINE MatStructP1Q2
!
! ----------------------------------------------
!
END MODULE QuadSc_struct
