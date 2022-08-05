MODULE HypreSolver
USE PP3D_MPI
IMPLICIT NONE

include 'HYPREf.h'

CONTAINS

SUBROUTINE myHypre_Solve(x,b,Mat)
USE var_QuadScalar
IMPLICIT NONE

REAL*8, intent(in) :: x(*), b(*)
TYPE(TMatrix), intent(in) :: Mat

integer*8 hypre_A, parcsr_A, hypre_b, par_b, hypre_x, par_x
integer*8 hypre_solver

integer ilower, iupper, local_size
integer hypreCommunicator, color
integer num_iterations
real*8  final_res_norm
integer ierr

! DELETE LATER; PLACEHOLDER
integer :: nrows
integer, dimension(:), allocatable :: ncols, rows, cols
real*8,  dimension(:), allocatable :: values, rhs_values, x_values


! initialize hypre structure
call HYPRE_Init(ierr)

call HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE, ierr)
call HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE, ierr)

call HYPRE_SetSpGemmUseVendor(0, ierr)

! set up parallel structure for hypre solvers
! ---------------------------------------------

! split communicator without id 0
if (myid.eq.0) then
  color = MPI_undefined
else
  color = 1
end if

call MPI_COMM_split(MPI_COMM_WORLD, color, myid-1, hypreCommunicator)

! get numbering of partitions
call hypre_factorize(ilower, iupper)

call HYPRE_IJMatrixCreate(hypreCommunicator, ilower, iupper,&
                         ilower, iupper, hypre_A, ierr)
call HYPRE_IJMatrixSetObjectType(hypre_A, HYPRE_PARCSR, ierr)
call HYPRE_IJMatrixInitialize(hypre_A, ierr)


! set up hypre matrix structure
! ---------------------------------------------

! get matrix structure for hypre partitions
call hypre_getMatrixStruct(nrows, ncols, rows, cols, values)
call HYPRE_IJMatrixSetValues(hypre_A, nrows, ncols, rows, cols, values, ierr)

call HYPRE_IJMatrixAssemble(hypre_A, ierr)
call HYPRE_IJMatrixGetObject(hypre_A, parcsr_A, ierr)

! set up parallel structure for RHS and solution vector
! ---------------------------------------------
call HYPRE_IJVectorCreate(hypreCommunicator,ilower, iupper, hypre_b, ierr)
call HYPRE_IJVectorSetObjectType(hypre_b, HYPRE_PARCSR, ierr)
call HYPRE_IJVectorInitialize(hypre_b, ierr)

call HYPRE_IJVectorCreate(hypreCommunicator,ilower, iupper, hypre_x, ierr)
call HYPRE_IJVectorSetObjectType(hypre_x, HYPRE_PARCSR, ierr)
call HYPRE_IJVectorInitialize(hypre_x, ierr)

! set vector values
call hypre_getVectorValues(rhs_values, x_values)
local_size = iupper - ilower + 1
call HYPRE_IJVectorSetValues(hypre_b, local_size, rows, rhs_values, ierr)
call HYPRE_IJVectorSetValues(hypre_x, local_size, rows, x_values, ierr)

call HYPRE_IJVectorAssemble(hypre_b, ierr)
call HYPRE_IJVectorAssemble(hypre_x, ierr)

call HYPRE_IJVectorGetObject(hypre_b, par_b, ierr)
call HYPRE_IJVectorGetObject(hypre_x, par_x, ierr)

!=============================================== TO DO: setup parameter values
! solve the system
! ---------------------------------------------
!        Create solver
call HYPRE_BoomerAMGCreate(hypre_solver, ierr)


!        Set some parameters (See Reference Manual for more parameters)

!        print solve info + parameters
call HYPRE_BoomerAMGSetPrintLevel(hypre_solver, 3, ierr)
!        old defaults, Falgout coarsening, mod. class. interpolation
call HYPRE_BoomerAMGSetOldDefault(hypre_solver, ierr)
!        G-S/Jacobi hybrid relaxation
call HYPRE_BoomerAMGSetRelaxType(hypre_solver, 3, ierr)
!        C/F relaxation
call HYPRE_BoomerAMGSetRelaxOrder(hypre_solver, 1, ierr)
!        Sweeeps on each level
call HYPRE_BoomerAMGSetNumSweeps(hypre_solver, 1, ierr)
!         maximum number of levels
call HYPRE_BoomerAMGSetMaxLevels(hypre_solver, 20, ierr)
!        conv. tolerance
call HYPRE_BoomerAMGSetTol(hypre_solver, 1.0d-7, ierr)
!        Keep local transposes
call HYPRE_BoomerAMGSetKeepTransp(hypre_solver, 1, ierr)

!        Now setup and solve!
call HYPRE_BoomerAMGSetup(hypre_solver, parcsr_A, par_b, par_x, ierr)
call HYPRE_BoomerAMGSolve(hypre_solver, parcsr_A, par_b, par_x, ierr)


!        Run info - needed logging turned on
call HYPRE_BoomerAMGGetNumIterations(hypre_solver, num_iterations,ierr)
call HYPRE_BoomerAMGGetFinalReltvRes(hypre_solver, final_res_norm,ierr)


if ( myid .eq. 0 ) then
  print *
  print '(A,I2)', " Iterations = ", num_iterations
  print '(A,ES16.8)'," Final Relative Residual Norm = ", final_res_norm
  print *
endif

!        Destroy solver
call HYPRE_BoomerAMGDestroy(hypre_solver, ierr)

end subroutine myHypre_Solve
!
! -----------------------------------------------------
!
subroutine hypre_factorize(ilower, iupper)
implicit none

integer, intent(inout) :: ilower,iupper
end subroutine
!
! -----------------------------------------------------
!
subroutine hypre_getMatrixStruct(nrows, ncols, rows, cols, values)
implicit none

integer, intent(in) :: nrows
integer, dimension(:), intent(inout) :: ncols, rows, cols
real*8,  dimension(:), intent(inout) :: values


end subroutine
!
! -----------------------------------------------------
!
subroutine hypre_getVectorValues(rhs_values, x_values)
implicit none

real*8,  dimension(:), intent(inout) :: rhs_values, x_values


end subroutine
!
! -----------------------------------------------------
!




END MODULE HypreSolver