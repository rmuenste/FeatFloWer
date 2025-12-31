!===============================================================================
! MODULE: QuadSc_assembly
!
! DESCRIPTION:
!   Matrix assembly routines for quadratic/linear scalar systems.
!   Extracted from def_QuadScalar module as part of Phase 3.2 structural split.
!
! PURPOSE:
!   - Centralizes all matrix coefficient assembly logic
!   - Generic assembly routines (from Phase 2 deduplication)
!   - Specific matrix builders (mass, diffusion, convection, coupling)
!   - Reduces complexity of monolithic def_QuadScalar module
!
! CONTAINS (12 PUBLIC ROUTINES):
!   Generic Assembly (Phase 2):
!     * Assemble_Mass_Generic              - Generic mass matrix assembly
!     * Assemble_Diffusion_Alpha_Generic   - Generic diffusion with alpha
!     * Assemble_ParallelMatrix_Generic    - Generic parallel matrix setup
!
!   Mass Matrices:
!     * Create_MRhoMat                     - Mass matrix with density (rho*M)
!     * Create_MMat                        - Standard mass matrix (M)
!
!   Diffusion Matrices:
!     * Create_hDiffMat                    - h-dependent diffusion
!     * Create_ConstDiffMat                - Constant diffusion
!     * Create_DiffMat                     - User-supplied scalar diffusion
!
!   Velocity-Pressure Coupling:
!     * Create_BMat                        - Gradient matrix (Q2→P1)
!     * Create_CMat                        - Divergence + pressure Laplacian
!
!   Convection Matrices:
!     * Create_SMat                        - Streamline diffusion matrix
!     * Create_KMat                        - Convection matrix
!
! DEPENDENCIES:
!   - var_QuadScalar: Data structures (mg_qMat, mg_lMat, etc.)
!   - PP3D_MPI: Parallel utilities (myid, showID)
!   - QuadSc_struct: Structure allocation routines (Create_LinMatStruct, etc.)
!
! AUTHOR: Extracted during def_QuadScalar refactoring Phase 3.2
! DATE: 2025
!===============================================================================

MODULE QuadSc_assembly

  USE PP3D_MPI, ONLY: myid, showID
  USE var_QuadScalar
  USE QuadSc_struct, ONLY: Create_LinMatStruct, Create_QuadLinMatStruct
  USE QuadSc_solver_coarse, ONLY: Setup_UMFPACK_CoarseSolver

  IMPLICIT NONE
  PRIVATE

  ! Public interface - Generic assembly routines
  PUBLIC :: Assemble_Mass_Generic
  PUBLIC :: Assemble_Diffusion_Alpha_Generic
  PUBLIC :: Assemble_ParallelMatrix_Generic

  ! Public interface - Specific matrix builders
  PUBLIC :: Create_MRhoMat
  PUBLIC :: Create_MMat
  PUBLIC :: Create_CMat
  PUBLIC :: Create_BMat
  PUBLIC :: Create_hDiffMat
  PUBLIC :: Create_ConstDiffMat
  PUBLIC :: Create_DiffMat
  PUBLIC :: Create_SMat
  PUBLIC :: Create_KMat

CONTAINS

!===============================================================================
! GENERIC ASSEMBLY ROUTINES (Phase 2)
!===============================================================================

SUBROUTINE Assemble_Mass_Generic(use_density, mg_MlMatrix, mg_MlPMatrix, &
                                  density_opt, label)
  LOGICAL, INTENT(IN) :: use_density
  TYPE(mg_Matrix), INTENT(INOUT) :: mg_MlMatrix(NLMIN:NLMAX)
  TYPE(mg_Matrix), INTENT(INOUT) :: mg_MlPMatrix(NLMIN:NLMAX)
  TYPE(mg_dVector), INTENT(IN), OPTIONAL :: density_opt(NLMIN:NLMAX)
  CHARACTER(LEN=*), INTENT(IN) :: label

  EXTERNAL E013
  REAL*8  DML
  INTEGER I, J

  ! Ensure mg_Mmat is allocated (shared by both MMat and MRhoMat)
  ! Note: mg_MlMatrix and mg_MlPMatrix are passed by caller, so they handle allocation
  IF (.not.ALLOCATED(mg_Mmat)) ALLOCATE(mg_Mmat(NLMIN:NLMAX))

  ! Loop over all multigrid levels
  DO ILEV=NLMIN,NLMAX
    CALL SETLEV(2)
    qMat => mg_qMat(ILEV)

    ! Allocate mass matrix if needed
    IF (.not.ALLOCATED(mg_Mmat(ILEV)%a)) THEN
      ALLOCATE(mg_Mmat(ILEV)%a(qMat%na))
    END IF

    mg_Mmat(ILEV)%a = 0d0

    ! Progress indicator
    IF (myid.eq.showID) THEN
      IF (ILEV.EQ.NLMIN) THEN
        WRITE(MTERM,'(A,I1,A)', advance='no') " "//TRIM(label)//": [", ILEV,"]"
      ELSE
        WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
      END IF
    END IF

    ! Call appropriate matrix builder
    IF (use_density) THEN
      IF (.NOT. PRESENT(density_opt)) THEN
        WRITE(*,*) 'ERROR: Assemble_Mass_Generic called with use_density=.TRUE. but no density array!'
        STOP
      END IF
      CALL BuildMRhoMat(density_opt(ILEV)%x, mg_Mmat(ILEV)%a, qMat%na, &
                        qMat%ColA, qMat%LdA, &
                        mg_mesh%level(ILEV)%kvert, &
                        mg_mesh%level(ILEV)%karea, &
                        mg_mesh%level(ILEV)%kedge, &
                        mg_mesh%level(ILEV)%dcorvg, &
                        E013)
    ELSE
      CALL BuildMMat(mg_Mmat(ILEV)%a, qMat%na, &
                     qMat%ColA, qMat%LdA, &
                     mg_mesh%level(ILEV)%kvert, &
                     mg_mesh%level(ILEV)%karea, &
                     mg_mesh%level(ILEV)%kedge, &
                     mg_mesh%level(ILEV)%dcorvg, &
                     E013)
    END IF

    ! Compute lumped mass matrix (row sum) - identical for both variants
    IF (.not.ALLOCATED(mg_MlMatrix(ILEV)%a)) ALLOCATE(mg_MlMatrix(ILEV)%a(qMat%nu))

    DO I=1,qMat%nu
      DML = 0d0
      DO J=qMat%LdA(I),qMat%LdA(I+1)-1
        DML = DML + mg_Mmat(ILEV)%a(J)
      END DO
      mg_MlMatrix(ILEV)%a(I) = DML
    END DO

    ! Parallel synchronization - identical for both variants
    IF (.not.ALLOCATED(mg_MlPMatrix(ILEV)%a)) ALLOCATE(mg_MlPMatrix(ILEV)%a(qMat%nu))
    mg_MlPMatrix(ILEV)%a = mg_MlMatrix(ILEV)%a
    CALL E013SUM(mg_MlPMatrix(ILEV)%a)

  END DO

  ! Restore to NLMAX and set pointers
  ILEV=NLMAX
  CALL SETLEV(2)
  qMat => mg_qMat(NLMAX)

  IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='no') " |"

END SUBROUTINE Assemble_Mass_Generic

!-------------------------------------------------------------------------------

SUBROUTINE Assemble_Diffusion_Alpha_Generic(mg_OutMatrix, alpha, label, &
                                             require_worker)
  TYPE(mg_Matrix), INTENT(INOUT) :: mg_OutMatrix(NLMIN:NLMAX)
  REAL*8, INTENT(IN) :: alpha
  CHARACTER(LEN=*), INTENT(IN) :: label
  LOGICAL, INTENT(IN) :: require_worker

  EXTERNAL E013

  CALL ZTIME(myStat%t0)

  ! Loop over multigrid levels
  DO ILEV = NLMIN, NLMAX

    CALL SETLEV(2)
    qMat => mg_qMat(ILEV)

    ! Allocate per-level matrix
    IF (.NOT. ALLOCATED(mg_OutMatrix(ILEV)%a)) THEN
      ALLOCATE(mg_OutMatrix(ILEV)%a(qMat%na))
    END IF

    mg_OutMatrix(ILEV)%a = 0d0

    ! Progress indicator
    IF (myid.eq.showID) THEN
      IF (ILEV.EQ.NLMIN) THEN
        WRITE(MTERM,'(A,I1,A)', advance='no') " "//TRIM(label)//": [", ILEV,"]"
      ELSE
        WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
      END IF
    END IF

    ! Assembly kernel (with optional worker-only guard)
    IF (require_worker) THEN
      IF (myid.ne.0) THEN
        CALL DIFFQ2_alpha(mg_OutMatrix(ILEV)%a, qMat%na, qMat%ColA, &
                          qMat%LdA, &
                          mg_mesh%level(ILEV)%kvert, &
                          mg_mesh%level(ILEV)%karea, &
                          mg_mesh%level(ILEV)%kedge, &
                          mg_mesh%level(ILEV)%dcorvg, &
                          E013, alpha)
      END IF
    ELSE
      CALL DIFFQ2_alpha(mg_OutMatrix(ILEV)%a, qMat%na, qMat%ColA, &
                        qMat%LdA, &
                        mg_mesh%level(ILEV)%kvert, &
                        mg_mesh%level(ILEV)%karea, &
                        mg_mesh%level(ILEV)%kedge, &
                        mg_mesh%level(ILEV)%dcorvg, &
                        E013, alpha)
    END IF

  END DO

  IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='no') " |"

  ! Restore to NLMAX level
  ILEV = NLMAX
  CALL SETLEV(2)

  CALL ZTIME(myStat%t1)
  myStat%tDMat = myStat%tDMat + (myStat%t1-myStat%t0)

END SUBROUTINE Assemble_Diffusion_Alpha_Generic

!-------------------------------------------------------------------------------

SUBROUTINE Assemble_ParallelMatrix_Generic(allocate_structure, myPLinSc_opt)
  LOGICAL, INTENT(IN) :: allocate_structure
  TYPE(TParLinScalar), INTENT(INOUT), OPTIONAL :: myPLinSc_opt

  ! Local variables
  INTEGER :: pNEL, MatSize
  INTEGER, ALLOCATABLE :: TempLdB(:)

  ! Top-level allocation (only if creating structure)
  IF (allocate_structure) THEN
    ALLOCATE(mg_BXPMat(NLMIN:NLMAX))
    ALLOCATE(mg_BYPMat(NLMIN:NLMAX))
    ALLOCATE(mg_BZPMat(NLMIN:NLMAX))
    ALLOCATE(mg_qlPMat(NLMIN:NLMAX))
  END IF

  ! Process all multigrid levels
  DO ILEV = NLMIN, NLMAX

    CALL SETLEV(2)

    ! Set up matrix pointers for this level
    qlMat => mg_qlMat(ILEV)
    BXMat => mg_BXMat(ILEV)%a
    BYMat => mg_BYMat(ILEV)%a
    BZMat => mg_BZMat(ILEV)%a

    ! Structure allocation (only if creating)
    IF (allocate_structure) THEN
      ! Initialize parallel pressure communication
      CALL ParPresComm_Init(qlMat%ColA, qlMat%LdA, qlMat%nu, &
                            mg_mesh%level(ILEV)%nel, ILEV)

      ! Allocate parallel matrix structure
      mg_qlPMat(ILEV)%nu = qlMat%nu
      ALLOCATE(mg_qlPMat(ILEV)%LdA(mg_qlPMat(ILEV)%nu+1))
      mg_qlPMat(ILEV)%LdA = 0

      ! Get parallel matrix row pointer structure
      CALL Create_ParB_LD(mg_qlPMat(ILEV)%LdA, qlMat%LdA, qlMat%nu, ILEV)

      MatSize = mg_qlPMat(ILEV)%LdA(mg_qlPMat(ILEV)%nu+1)

      IF (myid.eq.showID) THEN
        WRITE(MTERM,*) "Parallel B Matrix and Structure created", &
                       mg_qlPMat(ILEV)%nu, MatSize
      END IF

      ! Allocate column index and coefficient arrays
      ALLOCATE(mg_qlPMat(ILEV)%ColA(MatSize))
      mg_qlPMat(ILEV)%ColA = 0

      ALLOCATE(mg_BXPMat(ILEV)%a(MatSize))
      ALLOCATE(mg_BYPMat(ILEV)%a(MatSize))
      ALLOCATE(mg_BZPMat(ILEV)%a(MatSize))
    END IF

    ! Fill coefficient arrays (always executed for both create and fill)
    ALLOCATE(TempLdB(mg_qlPMat(ILEV)%nu+1))
    TempLdB = mg_qlPMat(ILEV)%LdA

    ! Get parallel matrix column structure and coefficients
    CALL Create_ParB_COLMAT(BXMat, BYMat, BZMat, &
                            mg_BXPMat(ILEV)%a, mg_BYPMat(ILEV)%a, mg_BZPMat(ILEV)%a, &
                            mg_qlPMat(ILEV)%LdA, mg_qlPMat(ILEV)%ColA, TempLdB, &
                            qlMat%LdA, qlMat%ColA, qlMat%nu, NEL, pNEL, ILEV)

    DEALLOCATE(TempLdB)

  END DO

  ! Set up output parameter structure (only if creating and parameter present)
  IF (allocate_structure .AND. PRESENT(myPLinSc_opt)) THEN
    myPLinSc_opt%ndof = pNEL
    IF (ALLOCATED(myPLinSc_opt%Val)) DEALLOCATE(myPLinSc_opt%Val)
    ALLOCATE(myPLinSc_opt%Val(4*pNEL))
  END IF

  ! Restore pointers to finest level (NLMAX)
  ILEV = NLMAX
  CALL SETLEV(2)

  qlMat  => mg_qlMat(NLMAX)
  qlPMat => mg_qlPMat(NLMAX)
  BXPMat => mg_BXPMat(NLMAX)%a
  BYPMat => mg_BYPMat(NLMAX)%a
  BZPMat => mg_BZPMat(NLMAX)%a
  BXMat  => mg_BXMat(NLMAX)%a
  BYMat  => mg_BYMat(NLMAX)%a
  BZMat  => mg_BZMat(NLMAX)%a

END SUBROUTINE Assemble_ParallelMatrix_Generic

!===============================================================================
! MASS MATRIX ROUTINES
!===============================================================================

SUBROUTINE Create_MRhoMat()
! Phase 2.1: Refactored to use Assemble_Mass_Generic
INTEGER :: ILEV_save

 if (.not.bMasterTurnedOn) return

 CALL ZTIME(myStat%t0)

 ! ========== New Implementation (Phase 2.1) ==========
 ! Allocate top-level arrays (must be done by caller before passing to generic)
 IF (.not.ALLOCATED(mg_MlRhomat))  ALLOCATE(mg_MlRhomat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_MlRhoPmat)) ALLOCATE(mg_MlRhoPmat(NLMIN:NLMAX))

 ! Call generic assembly routine with density
 CALL Assemble_Mass_Generic(.TRUE., mg_MlRhomat, mg_MlRhoPmat, &
                             mgDensity, "[MRho] & [MlRho]")

 ! Handle bSteadyState special case (only in MRhoMat, not MMat)
 if(bSteadyState)then
   DO ILEV_save=NLMIN,NLMAX
     mg_MMat(ILEV_save)%a = 0d0
   END DO
 end if
 ! ====================================================

 ! Set final pointers to NLMAX level
 qMat      => mg_qMat(NLMAX)
 Mmat      => mg_Mmat(NLMAX)%a
 MlRhomat  => mg_MlRhomat(NLMAX)%a
 MlRhoPmat => mg_MlRhoPmat(NLMAX)%a

 CALL ZTIME(myStat%t1)
 myStat%tMMat = myStat%tMMat + (myStat%t1-myStat%t0)

END SUBROUTINE Create_MRhoMat

!-------------------------------------------------------------------------------

SUBROUTINE Create_MMat()
! Phase 2.1: Refactored to use Assemble_Mass_Generic

 if (.not.bMasterTurnedOn) return

 CALL ZTIME(myStat%t0)

 ! ========== New Implementation (Phase 2.1) ==========
 ! Allocate top-level arrays (must be done by caller before passing to generic)
 IF (.not.ALLOCATED(mg_MlMat))  ALLOCATE(mg_MlMat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_MlPMat))  ALLOCATE(mg_MlPMat(NLMIN:NLMAX))

 ! Call generic assembly routine without density
 CALL Assemble_Mass_Generic(.FALSE., mg_MlMat, mg_MlPMat, &
                             label="[MRho] & [MlRho]")
 ! ====================================================

 ! Set final pointers to NLMAX level
 qMat      => mg_qMat(NLMAX)
 Mmat      => mg_Mmat(NLMAX)%a
 MlMat     => mg_MlMat(NLMAX)%a
 MlPMat    => mg_MlPMat(NLMAX)%a

 CALL ZTIME(myStat%t1)
 myStat%tMMat = myStat%tMMat + (myStat%t1-myStat%t0)

END SUBROUTINE Create_MMat

!===============================================================================
! DIFFUSION MATRIX ROUTINES
!===============================================================================

SUBROUTINE Create_hDiffMat()
! Phase 2.2: Refactored to use Assemble_Diffusion_Alpha_Generic

 ! ========== New Implementation (Phase 2.2) ==========
 ! Check if matrix already exists (early exit)
 IF (.NOT. ALLOCATED(mg_hDmat)) THEN
   ALLOCATE(mg_hDmat(NLMIN:NLMAX))
 ELSE
   IF (myid.eq.showID) THEN
     WRITE(MTERM,'(A)', advance='no') " [hD]: Exists |"
   END IF
   RETURN
 END IF

 ! Call generic assembly routine with alpha=1.0, worker-only guard
 CALL Assemble_Diffusion_Alpha_Generic(mg_hDmat, 1.0d0, "[hD]", .TRUE.)
 ! ====================================================

 ! Set final pointer to NLMAX level
 qMat  => mg_qMat(NLMAX)
 hDMat => mg_hDMat(NLMAX)%a

END SUBROUTINE Create_hDiffMat

!-------------------------------------------------------------------------------

SUBROUTINE Create_ConstDiffMat()
! Phase 2.2: Refactored to use Assemble_Diffusion_Alpha_Generic

 ! ========== New Implementation (Phase 2.2) ==========
 ! Check if matrix already exists (early exit)
 IF (.NOT. ALLOCATED(mg_ConstDMat)) THEN
   ALLOCATE(mg_ConstDMat(NLMIN:NLMAX))
 ELSE
   IF (myid.eq.showID) THEN
     WRITE(MTERM,'(A)', advance='no') " [VD]: Exists |"
   END IF
   RETURN
 END IF

 ! Call generic assembly routine with alpha=0.0, no worker guard
 CALL Assemble_Diffusion_Alpha_Generic(mg_ConstDMat, 0.0d0, "[VD]", .FALSE.)
 ! ====================================================

 ! Set final pointer to NLMAX level
 qMat      => mg_qMat(NLMAX)
 ConstDMat => mg_ConstDMat(NLMAX)%a

END SUBROUTINE Create_ConstDiffMat

!-------------------------------------------------------------------------------

SUBROUTINE Create_DiffMat(myScalar)
TYPE(TQuadScalar) myScalar
EXTERNAL E013

 if (.not.bMasterTurnedOn) return

 CALL ZTIME(myStat%t0)

 IF (.not.ALLOCATED(mg_Dmat)) ALLOCATE(mg_Dmat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)
  qMat => mg_qMat(ILEV)

  IF (.not.ALLOCATED(mg_Dmat(ILEV)%a)) THEN
   ALLOCATE(mg_Dmat(ILEV)%a(qMat%na))
  END IF

  mg_Dmat(ILEV)%a=0d0

  IF (myid.eq.showID) THEN
   IF (ILEV.EQ.NLMIN) THEN
    WRITE(MTERM,'(A,I1,A)', advance='no') " [D]: [", ILEV,"]"
   ELSE
    WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
   END IF
  END IF

  if(bNonNewtonian) THEN
   if (bMultiMat) then
    CALL DIFFQ2_AlphaNNEWT(myScalar%valU, myScalar%valV,myScalar%valW, &
         MaterialDistribution(ILEV)%x,&
         mg_Dmat(ILEV)%a,qMat%na,qMat%ColA,qMat%LdA,&
         mg_mesh%level(ILEV)%kvert,&
         mg_mesh%level(ILEV)%karea,&
         mg_mesh%level(ILEV)%kedge,&
         mg_mesh%level(ILEV)%dcorvg,&
         E013)
   else
    CALL DIFFQ2_NNEWT(myScalar%valU, myScalar%valV,myScalar%valW, &
         Temperature,& !MaterialDistribution(ILEV)%x,&
         mg_Dmat(ILEV)%a,qMat%na,qMat%ColA,qMat%LdA,&
         mg_mesh%level(ILEV)%kvert,&
         mg_mesh%level(ILEV)%karea,&
         mg_mesh%level(ILEV)%kedge,&
         mg_mesh%level(ILEV)%dcorvg,&
         E013)
   end if
  else
    CALL DIFFQ2_NEWT(mg_Dmat(ILEV)%a,qMat%na,qMat%ColA,&
         qMat%LdA,&
         mg_mesh%level(ILEV)%kvert,&
         mg_mesh%level(ILEV)%karea,&
         mg_mesh%level(ILEV)%kedge,&
         mg_mesh%level(ILEV)%dcorvg,&
         E013)
  end if

 END DO

 IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='no') " |"

 ILEV=NLMAX
 CALL SETLEV(2)

 qMat  => mg_qMat(NLMAX)
 DMat => mg_DMat(NLMAX)%a

 CALL ZTIME(myStat%t1)
 myStat%tDMat = myStat%tDMat + (myStat%t1-myStat%t0)

END SUBROUTINE Create_DiffMat

!===============================================================================
! VELOCITY-PRESSURE COUPLING MATRICES
!===============================================================================

SUBROUTINE Create_BMat() !(B)
INTEGER nERow,pNEL
INTEGER I,J
real*8 ddx,ddy,ddz
CHARACTER*10 myFile
EXTERNAL E011,E013

 if (.not.bMasterTurnedOn) return

 IF (.NOT.ALLOCATED(mg_BXMat)) ALLOCATE(mg_BXMat(NLMIN:NLMAX))
 IF (.NOT.ALLOCATED(mg_BYMat)) ALLOCATE(mg_BYMat(NLMIN:NLMAX))
 IF (.NOT.ALLOCATED(mg_BZMat)) ALLOCATE(mg_BZMat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)

  qlMat => mg_qlMat(ILEV)
  IF (.NOT.ALLOCATED(mg_BXMat(ILEV)%a)) ALLOCATE(mg_BXMat(ILEV)%a(qlMat%na))
  IF (.NOT.ALLOCATED(mg_BYMat(ILEV)%a)) ALLOCATE(mg_BYMat(ILEV)%a(qlMat%na))
  IF (.NOT.ALLOCATED(mg_BZMat(ILEV)%a)) ALLOCATE(mg_BZMat(ILEV)%a(qlMat%na))
  mg_BXMat(ILEV)%a=0d0
  mg_BYMat(ILEV)%a=0d0
  mg_BZMat(ILEV)%a=0d0

  IF (myid.eq.showID) THEN
   IF (ILEV.EQ.NLMIN) THEN
    WRITE(MTERM,'(A,I1,A)', advance='no') " [B]: [", ILEV,"]"
   ELSE
    WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
   END IF
  END IF
  CALL Build_BMatP1(mg_BXMat(ILEV)%a,mg_BYMat(ILEV)%a,&
       mg_BZMat(ILEV)%a,qlMat%LdA,qlMat%ColA,&
       mg_mesh%level(ILEV)%kvert,&
       mg_mesh%level(ILEV)%karea,&
       mg_mesh%level(ILEV)%kedge,&
       mg_mesh%level(ILEV)%dcorvg,&
       qlMat%na,E013)

 END DO

 IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='no') " |"

 IF (.NOT.ALLOCATED(mg_BTXMat)) ALLOCATE(mg_BTXMat(NLMIN:NLMAX))
 IF (.NOT.ALLOCATED(mg_BTYMat)) ALLOCATE(mg_BTYMat(NLMIN:NLMAX))
 IF (.NOT.ALLOCATED(mg_BTZMat)) ALLOCATE(mg_BTZMat(NLMIN:NLMAX))

DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)

  lqMat => mg_lqMat(ILEV)
  IF (.NOT.ALLOCATED(mg_BTXMat(ILEV)%a)) ALLOCATE(mg_BTXMat(ILEV)%a(lqMat%na))
  IF (.NOT.ALLOCATED(mg_BTYMat(ILEV)%a)) ALLOCATE(mg_BTYMat(ILEV)%a(lqMat%na))
  IF (.NOT.ALLOCATED(mg_BTZMat(ILEV)%a)) ALLOCATE(mg_BTZMat(ILEV)%a(lqMat%na))
  mg_BTXMat(ILEV)%a=0d0
  mg_BTYMat(ILEV)%a=0d0
  mg_BTZMat(ILEV)%a=0d0


  IF (myid.eq.showID) THEN
   IF (ILEV.EQ.NLMIN) THEN
    WRITE(MTERM,'(A,I1,A)', advance='no') " [B{T}]: [", ILEV,"]"
   ELSE
    WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
   END IF
  END IF
  CALL Build_BTMatP1(mg_BTXMat(ILEV)%a,mg_BTYMat(ILEV)%a,&
       mg_BTZMat(ILEV)%a,lqMat%LdA,lqMat%ColA,&
       mg_mesh%level(ILEV)%kvert,&
       mg_mesh%level(ILEV)%karea,&
       mg_mesh%level(ILEV)%kedge,&
       mg_mesh%level(ILEV)%dcorvg,&
       lqMat%na,E013)

 END DO

 IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='no') " |"

55 CONTINUE

 ILEV=NLMAX
 CALL SETLEV(2)

 qlMat  => mg_qlMat(NLMAX)
 BXMat => mg_BXMat(NLMAX)%a
 BYMat => mg_BYMat(NLMAX)%a
 BZMat => mg_BZMat(NLMAX)%a

 lqMat  => mg_lqMat(NLMAX)
 BTXMat => mg_BTXMat(NLMAX)%a
 BTYMat => mg_BTYMat(NLMAX)%a
 BTZMat => mg_BTZMat(NLMAX)%a

END SUBROUTINE Create_BMat

!-------------------------------------------------------------------------------

SUBROUTINE Create_CMat(knprU,knprV,knprW,knprP,coarse_lev,coarse_solver) !(C)
INTEGER coarse_lev,coarse_solver
TYPE(mg_kVector) :: knprU(*),knprV(*),knprW(*),knprP(*)
INTEGER i,j,iEntry,jCol

 if (.not.bMasterTurnedOn) return

 CALL ZTIME(myStat%t0)

 IF (.NOT.ALLOCATED(mg_CMat)) ALLOCATE(mg_CMat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)
  qMat      => mg_qMat(ILEV)
  lMat      => mg_lMat(ILEV)
  qlMat     => mg_qlMat(ILEV)
  lqMat     => mg_lqMat(ILEV)
  BXMat     => mg_BXMat(ILEV)%a
  BYMat     => mg_BYMat(ILEV)%a
  BZMat     => mg_BZMat(ILEV)%a
  BTXMat    => mg_BTXMat(ILEV)%a
  BTYMat    => mg_BTYMat(ILEV)%a
  BTZMat    => mg_BTZMat(ILEV)%a
  MlRhoPmat => mg_MlRhoPmat(ILEV)%a

  IF (.NOT.ALLOCATED(mg_CMat(ILEV)%a)) ALLOCATE(mg_CMat(ILEV)%a(lMat%na))
  mg_CMat(ILEV)%a=0d0
  CALL Get_CMat(MlRhoPmat,mg_CMat(ILEV)%a,lMat%LdA,lMat%ColA,&
       BXMat,BYMat,BZMat,qlMat%LdA,qlMat%ColA,&
       BTXMat,BTYMat,BTZMat,lqMat%LdA,lqMat%ColA, &
       knprU(ILEV)%x,knprV(ILEV)%x,knprW(ILEV)%x,lMat%nu,qMat%nu)

  IF (myid.eq.showID) THEN
   IF (ILEV.EQ.NLMIN) THEN
    WRITE(MTERM,'(A,I1,A)', advance='no') " [B{T} MRho{-1} B]: [", ILEV,"]"
   ELSE
    WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
   END IF
  END IF
 END DO

 ILEV=NLMAX
 CALL SETLEV(2)

 qMat      => mg_qMat(NLMAX)
 lMat      => mg_lMat(NLMAX)
 qlMat     => mg_qlMat(NLMAX)
 lqMat     => mg_lqMat(NLMAX)
 BXMat     => mg_BXMat(NLMAX)%a
 BYMat     => mg_BYMat(NLMAX)%a
 BZMat     => mg_BZMat(NLMAX)%a
 BTXMat    => mg_BTXMat(NLMAX)%a
 BTYMat    => mg_BTYMat(NLMAX)%a
 BTZMat    => mg_BTZMat(NLMAX)%a
 CMat      => mg_CMat(NLMAX)%a
 MlRhoPmat => mg_MlRhoPmat(NLMAX)%a

 ! Phase 1: UMFPACK coarse grid solver setup (types 2, 3, 4)
 ! Call extracted UMFPACK setup routine
 CALL Setup_UMFPACK_CoarseSolver(knprP, coarse_solver, coarse_lev, bNoOutflow)

 CALL ZTIME(myStat%t1)
 myStat%tCMat = myStat%tCMat + (myStat%t1-myStat%t0)
 IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='no') " |"


END SUBROUTINE Create_CMat

!===============================================================================
! CONVECTION MATRIX ROUTINES
!===============================================================================

SUBROUTINE Create_SMat(myScalar)
TYPE(TQuadScalar) myScalar
EXTERNAL E013
INTEGER i

 CALL ZTIME(myStat%t0)

 IF (.not.ALLOCATED(mg_S11mat)) ALLOCATE(mg_S11mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_S22mat)) ALLOCATE(mg_S22mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_S33mat)) ALLOCATE(mg_S33mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_S12mat)) ALLOCATE(mg_S12mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_S13mat)) ALLOCATE(mg_S13mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_S23mat)) ALLOCATE(mg_S23mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_S21mat)) ALLOCATE(mg_S21mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_S31mat)) ALLOCATE(mg_S31mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_S32mat)) ALLOCATE(mg_S32mat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)
  qMat => mg_qMat(ILEV)

  IF (.not.ALLOCATED(mg_S11mat(ILEV)%a)) ALLOCATE(mg_S11mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_S22mat(ILEV)%a)) ALLOCATE(mg_S22mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_S33mat(ILEV)%a)) ALLOCATE(mg_S33mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_S12mat(ILEV)%a)) ALLOCATE(mg_S12mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_S13mat(ILEV)%a)) ALLOCATE(mg_S13mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_S23mat(ILEV)%a)) ALLOCATE(mg_S23mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_S21mat(ILEV)%a)) ALLOCATE(mg_S21mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_S31mat(ILEV)%a)) ALLOCATE(mg_S31mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_S32mat(ILEV)%a)) ALLOCATE(mg_S32mat(ILEV)%a(qMat%na))

  mg_S11mat(ILEV)%a=0d0
  mg_S22mat(ILEV)%a=0d0
  mg_S33mat(ILEV)%a=0d0
  mg_S12mat(ILEV)%a=0d0
  mg_S13mat(ILEV)%a=0d0
  mg_S23mat(ILEV)%a=0d0
  mg_S21mat(ILEV)%a=0d0
  mg_S31mat(ILEV)%a=0d0
  mg_S32mat(ILEV)%a=0d0

  IF (myid.eq.showID) THEN
   IF (ILEV.EQ.NLMIN) THEN
    WRITE(MTERM,'(A,I1,A)', advance='no') " [S]: [", ILEV,"]"
   ELSE
    WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
   END IF
  END IF

  CALL CUBATURESTRESS(myScalar%valU, myScalar%valV,myScalar%valW, &
       Temperature,&
       mg_S11mat(ILEV)%a,mg_S22mat(ILEV)%a,mg_S33mat(ILEV)%a,&
       mg_S12mat(ILEV)%a,mg_S13mat(ILEV)%a,mg_S23mat(ILEV)%a,&
       mg_S21mat(ILEV)%a,mg_S31mat(ILEV)%a,mg_S32mat(ILEV)%a,&
       qMat%na,qMat%ColA,qMat%LdA,&
       mg_mesh%level(ILEV)%kvert,&
       mg_mesh%level(ILEV)%karea,&
       mg_mesh%level(ILEV)%kedge,&
       mg_mesh%level(ILEV)%dcorvg,&
       E013)

 END DO

 ILEV=NLMAX
 CALL SETLEV(2)

 qMat  => mg_qMat(NLMAX)
 S11Mat => mg_S11Mat(NLMAX)%a
 S22Mat => mg_S22Mat(NLMAX)%a
 S33Mat => mg_S33Mat(NLMAX)%a
 S12Mat => mg_S12Mat(NLMAX)%a
 S13Mat => mg_S13Mat(NLMAX)%a
 S23Mat => mg_S23Mat(NLMAX)%a
 S21Mat => mg_S21Mat(NLMAX)%a
 S31Mat => mg_S31Mat(NLMAX)%a
 S32Mat => mg_S32Mat(NLMAX)%a

 CALL ZTIME(myStat%t1)
 myStat%tSMat = myStat%tSMat + (myStat%t1-myStat%t0)

END SUBROUTINE Create_SMat

!-------------------------------------------------------------------------------

SUBROUTINE Create_KMat(myScalar)
TYPE(TQuadScalar) myScalar
INTEGER LINT
EXTERNAL E013
! Assembly for Convection Kmat

 if (.not.bMasterTurnedOn) return

 CALL ZTIME(myStat%t0)

 IF (.not.ALLOCATED(mg_Kmat)) ALLOCATE(mg_Kmat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)
  qMat => mg_qMat(ILEV)

  IF (.not.ALLOCATED(mg_Kmat(ILEV)%a)) THEN
   ALLOCATE(mg_Kmat(ILEV)%a(qMat%na))
  END IF

  mg_Kmat(ILEV)%a=0d0

  IF (myid.eq.showID) THEN
   IF (ILEV.EQ.NLMIN) THEN
    WRITE(MTERM,'(A,I1,A)', advance='no') " [KRho]: [", ILEV,"]"
   ELSE
    WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
   END IF
  END IF

  LINT = KLINT(ILEV)

  CALL CONVQ2(mgDensity(ILEV)%x,myScalar%valU,myScalar%valV,myScalar%valW,&
  myALE%MeshVelo,&
  mg_Kmat(ILEV)%a,qMat%nu,qMat%ColA,qMat%LdA,&
  mg_mesh%level(ILEV)%kvert,&
  mg_mesh%level(ILEV)%karea,&
  mg_mesh%level(ILEV)%kedge,&
  mg_mesh%level(ILEV)%dcorvg,&
  E013)

 END DO

 IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='no') " |"

 ILEV=NLMAX
 CALL SETLEV(2)

 qMat  => mg_qMat(NLMAX)
 KMat => mg_KMat(NLMAX)%a

 CALL ZTIME(myStat%t1)
 myStat%tKMat = myStat%tKMat + (myStat%t1-myStat%t0)

END SUBROUTINE Create_KMat

!===============================================================================
END MODULE QuadSc_assembly
