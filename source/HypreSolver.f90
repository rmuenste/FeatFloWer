MODULE HypreSolver
USE PP3D_MPI
IMPLICIT NONE

include 'HYPREf.h'

#ifdef HYPRE_TIMING
integer :: hypre_amg_call_count = 0
integer :: hypre_gmres_call_count = 0
integer :: hypre_pcg_call_count = 0
#endif

CONTAINS

SUBROUTINE myHypre_Solve(num_iterations)
USE var_QuadScalar
IMPLICIT NONE

integer local_size
integer color
real*8  final_res_norm
integer ierr
integer status(MPI_status_size)

#ifdef HYPRE_TIMING
real*8 MPI_WTIME
real*8 total_start, total_end, setup_start, setup_end, solve_start, solve_end
#endif

integer, intent(inout) :: num_iterations


if (.not.myHypre%solverIsSet) then
  ! initialize hypre structure
  call HYPRE_Init(ierr)

  call HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE, ierr)
  call HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE, ierr)

  call HYPRE_SetSpGemmUseVendor(0, ierr)

  ! set up parallel structure for hypre solvers
  ! ---------------------------------------------
  if (myid.eq.0) then
    color = MPI_undefined
  else
    color = 1
  end if

  call MPI_COMM_split(MPI_COMM_WORLD, color, myid-1, myHypre%communicator, ierr)

  if (myid.ne.0) then

  call HYPRE_IJMatrixCreate(myHypre%communicator, myHypre%ilower, myHypre%iupper,&
                         myHypre%ilower, myHypre%iupper, myHypre%A, ierr)
  call HYPRE_IJMatrixSetObjectType(myHypre%A, HYPRE_PARCSR, ierr)
  call HYPRE_IJMatrixInitialize(myHypre%A, ierr)

  ! set up hypre matrix structure
  ! ---------------------------------------------

  ! get matrix structure for hypre partitions

  call HYPRE_IJMatrixSetValues(myHypre%A, myHypre%nrows, myHypre%ncols, myHypre%rows, myHypre%cols, myHypre%values, ierr)

  call HYPRE_IJMatrixAssemble(myHypre%A, ierr)
  call HYPRE_IJMatrixGetObject(myHypre%A, myHypre%parcsr_A, ierr)

  ! set up parallel structure for RHS and solution vector
  ! ---------------------------------------------
  call HYPRE_IJVectorCreate(myHypre%communicator,myHypre%ilower, myHypre%iupper, myHypre%b, ierr)
  call HYPRE_IJVectorSetObjectType(myHypre%b, HYPRE_PARCSR, ierr)
  call HYPRE_IJVectorInitialize(myHypre%b, ierr)

  call HYPRE_IJVectorCreate(myHypre%communicator,myHypre%ilower, myHypre%iupper, myHypre%x, ierr)
  call HYPRE_IJVectorSetObjectType(myHypre%x, HYPRE_PARCSR, ierr)
  call HYPRE_IJVectorInitialize(myHypre%x, ierr)


  call HYPRE_IJVectorGetObject(myHypre%b, myHypre%par_b, ierr)
  call HYPRE_IJVectorGetObject(myHypre%x, myHypre%par_x, ierr)

  !=============================================== TO DO: setup parameter values
  ! solve the system
  ! ---------------------------------------------
  !        Create solver
  call HYPRE_BoomerAMGCreate(myHypre%solver, ierr)

  !print solve info + parameters
  call HYPRE_BoomerAMGSetPrintLevel(myHypre%solver, 0, ierr)

  call HYPRE_BoomerAMGSetStrongThrshld(myHypre%solver, 0.25d0, ierr)
#ifdef HYPRE_GPU_AMG
  ! PMIS and symmetric hybrid Gauss-Seidel/Jacobi are GPU-capable.
  call HYPRE_BoomerAMGSetCoarsenType(myHypre%solver, 8, ierr)
  call HYPRE_BoomerAMGSetRelaxType(myHypre%solver, 6, ierr)
#else
  !use HMIS coarsening (robust for 3D Laplace systems)
  call HYPRE_BoomerAMGSetCoarsenType(myHypre%solver, 10, ierr)
  !G-S/Jacobi hybrid relaxation
  call HYPRE_BoomerAMGSetRelaxType(myHypre%solver, 8, ierr)
#endif
  !C/F relaxation
  call HYPRE_BoomerAMGSetRelaxOrder(myHypre%solver, 1, ierr)
  !Sweeps on each level
  call HYPRE_BoomerAMGSetNumSweeps(myHypre%solver, 2, ierr)
  !maximum number of levels
  call HYPRE_BoomerAMGSetMaxLevels(myHypre%solver, 20, ierr)
#ifdef HYPRE_GPU_AMG
  ! Extended+i interpolation with truncation is GPU-capable.
  call Hypre_BoomerAMGSetInterpType(myHypre%solver, 6, ierr)
  call HYPRE_BoomerAMGSetTruncFactor(myHypre%solver, 0.25d0, ierr)
#else
  !set interpolation type
  call Hypre_BoomerAMGSetInterpType(myHypre%solver, 13, ierr)
#endif
  !Max numbers per rows
  call HYPRE_BoomerAMGSetPMaxElmts(myHypre%solver, 5, ierr)

!   call HYPRE_BoomerAMGSetNumFunctions(myHypre%solver, 4, ierr)
!   call HYPRE_BoomerAMGSetNodal(myHypre%solver, 3, ierr)

  call HYPRE_BoomerAMGSetCycleType(myHypre%solver, 1, ierr)

  call Hypre_BoomerAMGSetMaxIter(myHypre%solver, 50, ierr)


  !conv. tolerance
  call HYPRE_BoomerAMGSetTol(myHypre%solver, 1.0d-7, ierr)

  end if
end if ! solver is set
myHypre%solverIsSet = .true.

#ifdef HYPRE_TIMING
hypre_amg_call_count = hypre_amg_call_count + 1
#endif

  !setup and solve
if (myid.ne.0) then
#ifdef HYPRE_TIMING
  total_start = MPI_WTIME()
#endif
  ! set vector values
  local_size = myHypre%iupper - myHypre%ilower + 1
  call HYPRE_IJVectorSetValues(myHypre%b, local_size, myHypre%rows, myHypre%rhs, ierr)
  call HYPRE_IJVectorSetValues(myHypre%x, local_size, myHypre%rows, myHypre%sol, ierr)

  call HYPRE_IJVectorAssemble(myHypre%b, ierr)
  call HYPRE_IJVectorAssemble(myHypre%x, ierr)

#ifdef HYPRE_TIMING
  setup_start = MPI_WTIME()
#endif
  call HYPRE_BoomerAMGSetup(myHypre%solver, myHypre%parcsr_A, myHypre%par_b, myHypre%par_x, ierr)
#ifdef HYPRE_TIMING
  setup_end = MPI_WTIME()
  solve_start = MPI_WTIME()
#endif
  call HYPRE_BoomerAMGSolve(myHypre%solver, myHypre%parcsr_A, myHypre%par_b, myHypre%par_x, ierr)
#ifdef HYPRE_TIMING
  solve_end = MPI_WTIME()
#endif

  !        Run info - needed logging turned on
  call HYPRE_BoomerAMGGetNumIterations(myHypre%solver, num_iterations,ierr)
  call HYPRE_BoomerAMGGetFinalReltvRes(myHypre%solver, final_res_norm,ierr)
end if


if (myid.ne.0) then
 ! Recover the values from HYPRE back to "x"
 call HYPRE_IJVectorGetValues(myHypre%x, local_size, myHypre%rows, myHypre%sol, ierr)
#ifdef HYPRE_TIMING
 total_end = MPI_WTIME()
 call ReportHypreTiming('amg', hypre_amg_call_count, setup_end - setup_start, &
                        solve_end - solve_start, total_end - total_start, &
                        num_iterations, final_res_norm)
#endif
end if


  if (myid.eq.1) then
   call MPI_send(num_iterations, 1, MPI_int, 0, 1, MPI_COMM_WORLD, ierr)
  else if (myid.eq.0) then
   call MPI_recv(num_iterations, 1, MPI_int, 1, 1, MPI_COMM_WORLD, status, ierr)
  end if

! pause
end subroutine myHypre_Solve



SUBROUTINE myHypre_GMRESDestroy()
USE var_QuadScalar
IMPLICIT NONE
integer ierr

if (myid.ne.0) then
  call HYPRE_IJVectorDestroy(myHypre%b, ierr)
  call HYPRE_IJVectorDestroy(myHypre%x, ierr)
  call HYPRE_IJMatrixDestroy(myHypre%A, ierr)
  call HYPRE_ParCSRGMRESDestroy(myHypre%solver, ierr)
  call HYPRE_BoomerAMGDestroy(myHypre%precond, ierr)
end if

END SUBROUTINE myHypre_GMRESDestroy



SUBROUTINE myHypreGMRES_Solve(num_iterations)
USE var_QuadScalar
IMPLICIT NONE


integer local_size
integer color
real*8  final_res_norm
integer ierr
integer status(MPI_status_size)

#ifdef HYPRE_TIMING
real*8 MPI_WTIME
real*8 total_start, total_end, setup_start, setup_end, solve_start, solve_end
#endif

integer, intent(inout) :: num_iterations


if (.not.myHypre%solverIsSet) then
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

  call MPI_COMM_split(MPI_COMM_WORLD, color, myid-1, myHypre%communicator, ierr)
  
  if (myid.ne.0) then

  call HYPRE_IJMatrixCreate(myHypre%communicator, myHypre%ilower, myHypre%iupper,&
                           myHypre%ilower, myHypre%iupper, myHypre%A, ierr)
  call HYPRE_IJMatrixSetObjectType(myHypre%A, HYPRE_PARCSR, ierr)
  call HYPRE_IJMatrixInitialize(myHypre%A, ierr)

  ! set up hypre matrix structure
  ! ---------------------------------------------

  call HYPRE_IJMatrixSetValues(myHypre%A, myHypre%nrows, myHypre%ncols, myHypre%rows, myHypre%cols, myHypre%values, ierr)

  call HYPRE_IJMatrixAssemble(myHypre%A, ierr)
  call HYPRE_IJMatrixGetObject(myHypre%A, myHypre%parcsr_A, ierr)

  ! set up parallel structure for RHS and solution vector
  ! ---------------------------------------------
  call HYPRE_IJVectorCreate(myHypre%communicator,myHypre%ilower, myHypre%iupper, myHypre%b, ierr)
  call HYPRE_IJVectorSetObjectType(myHypre%b, HYPRE_PARCSR, ierr)
  call HYPRE_IJVectorInitialize(myHypre%b, ierr)

  call HYPRE_IJVectorCreate(myHypre%communicator,myHypre%ilower, myHypre%iupper, myHypre%x, ierr)
  call HYPRE_IJVectorSetObjectType(myHypre%x, HYPRE_PARCSR, ierr)
  call HYPRE_IJVectorInitialize(myHypre%x, ierr)

  ! Setup the GMRES Solver
  ! ---------------------------------------------
  call HYPRE_ParCSRGMRESCreate(myHypre%communicator, myHypre%solver, ierr)
  call HYPRE_ParCSRGMRESSetMaxIter(myHypre%solver, 80, ierr)
  call HYPRE_ParCSRGMRESSetTol(myHypre%solver, 1.0d-5, ierr)
  call HYPRE_ParCSRGMRESSetPrintLevel(myHypre%solver, 0, ierr)
  call HYPRE_ParCSRGMRESSetLogging(myHypre%solver, 0, ierr)

   
  ! setup the preconditioner
  ! ---------------------------------------------
  !Create solver
  call HYPRE_BoomerAMGCreate(myHypre%precond, ierr)

  call HYPRE_BoomerAMGSetStrongThrshld(myHypre%precond, 0.25d0, ierr)
#ifdef HYPRE_GPU_AMG
  call HYPRE_BoomerAMGSetCoarsenType(myHypre%precond, 8, ierr)
  call HYPRE_BoomerAMGSetRelaxType(myHypre%precond, 6, ierr)
#else
  !set coarsening type (HMIS for scalar Laplace)
  call HYPRE_BoomerAMGSetCoarsenType(myHypre%precond, 10, ierr)
  !G-S/Jacobi hybrid relaxation
  call HYPRE_BoomerAMGSetRelaxType(myHypre%precond, 8, ierr)
#endif
  !C/F relaxation
  call HYPRE_BoomerAMGSetRelaxOrder(myHypre%precond, 1, ierr)
  ! Sweeeps on each level
  call HYPRE_BoomerAMGSetNumSweeps(myHypre%precond, 2, ierr)
  !maximum number of levels
  call HYPRE_BoomerAMGSetMaxLevels(myHypre%precond, 20, ierr)
#ifdef HYPRE_GPU_AMG
  call Hypre_BoomerAMGSetInterpType(myHypre%precond, 6, ierr)
  call HYPRE_BoomerAMGSetTruncFactor(myHypre%precond, 0.25d0, ierr)
#else
  !set interpolation type
  call Hypre_BoomerAMGSetInterpType(myHypre%precond, 13, ierr)
#endif
  !Max numbers per rows
  call HYPRE_BoomerAMGSetPMaxElmts(myHypre%precond, 5, ierr)
  
  CALL HYPRE_BoomerAMGSetMinCoarseSize(myHypre%precond, 60, ierr)

  CALL HYPRE_BoomerAMGSetMaxCoarseSize(myHypre%precond, 240, ierr)
  
  CALL HYPRE_BoomerAMGSetCoarsenType(myHypre%precond,8,ierr)

#ifndef HYPRE_GPU_AMG
  ! The Hypre 2.33 scalar C/F expansion for nodal systems is host-only.
  call HYPRE_BoomerAMGSetNumFunctions(myHypre%precond, 4, ierr)
  call HYPRE_BoomerAMGSetNodal(myHypre%precond, 3, ierr)
#endif

  call HYPRE_BoomerAMGSetCycleType(myHypre%precond, 1, ierr)

  call Hypre_BoomerAMGSetMaxIter(myHypre%precond, 1, ierr)

  call HYPRE_BoomerAMGSetTol(myHypre%precond, 0.0d0, ierr)



  !set amg as the pcg preconditioner, precond_id = 2 -> AMG
  call HYPRE_ParCSRGMRESSetPrecond(myHypre%solver, 2, myHypre%precond, ierr)

  end if
end if ! solver is Setup
myHypre%solverIsSet = .true.

#ifdef HYPRE_TIMING
hypre_gmres_call_count = hypre_gmres_call_count + 1
#endif

if (myid.ne.0) then
#ifdef HYPRE_TIMING
  total_start = MPI_WTIME()
#endif
  ! set vector values
  local_size = myHypre%iupper - myHypre%ilower + 1
  call HYPRE_IJVectorSetValues(myHypre%b, local_size, myHypre%rows, myHypre%rhs, ierr)
  call HYPRE_IJVectorSetValues(myHypre%x, local_size, myHypre%rows, myHypre%sol, ierr)

  call HYPRE_IJVectorAssemble(myHypre%b, ierr)
  call HYPRE_IJVectorAssemble(myHypre%x, ierr)

  call HYPRE_IJVectorGetObject(myHypre%b, myHypre%par_b, ierr)
  call HYPRE_IJVectorGetObject(myHypre%x, myHypre%par_x, ierr)
  
  ! setup the solver and solve the system
#ifdef HYPRE_TIMING
  setup_start = MPI_WTIME()
#endif
  call HYPRE_ParCSRGMRESSetup(myHypre%solver, myHypre%parcsr_A, myHypre%par_b,myHypre%par_x, ierr)
#ifdef HYPRE_TIMING
  setup_end = MPI_WTIME()
  solve_start = MPI_WTIME()
#endif
  call HYPRE_ParCSRGMRESSolve(myHypre%solver, myHypre%parcsr_A, myHypre%par_b,myHypre%par_x, ierr)
#ifdef HYPRE_TIMING
  solve_end = MPI_WTIME()
#endif

  call HYPRE_ParCSRGMRESGetNumIteratio(myHypre%solver, num_iterations, ierr)
  call HYPRE_ParCSRGMRESGetFinalRelati(myHypre%solver, final_res_norm, ierr)
  
end if


if (myid.ne.0) then
 ! Recover the values from HYPRE back to "x"
 call HYPRE_IJVectorGetValues(myHypre%x, local_size, myHypre%rows, myHypre%sol, ierr)
#ifdef HYPRE_TIMING
 total_end = MPI_WTIME()
 call ReportHypreTiming('gmres-amg', hypre_gmres_call_count, setup_end - setup_start, &
                        solve_end - solve_start, total_end - total_start, &
                        num_iterations, final_res_norm)
#endif
end if


  if (myid.eq.1) then
   call MPI_send(num_iterations, 1, MPI_int, 0, 1, MPI_COMM_WORLD, ierr)
  else if (myid.eq.0) then
   call MPI_recv(num_iterations, 1, MPI_int, 1, 1, MPI_COMM_WORLD, status, ierr)
  end if
end subroutine myHypreGMRES_Solve






SUBROUTINE myHyprePCG_Solve
USE var_QuadScalar
IMPLICIT NONE


integer local_size
integer color
integer num_iterations
real*8  final_res_norm
integer ierr
integer status(MPI_status_size)

#ifdef HYPRE_TIMING
real*8 MPI_WTIME
real*8 total_start, total_end, setup_start, setup_end, solve_start, solve_end
#endif


if (.not.myHypre%solverIsSet) then
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

  call MPI_COMM_split(MPI_COMM_WORLD, color, myid-1, myHypre%communicator, ierr)

  if (myid.ne.0) then

  call HYPRE_IJMatrixCreate(myHypre%communicator, myHypre%ilower, myHypre%iupper,&
                           myHypre%ilower, myHypre%iupper, myHypre%A, ierr)
  call HYPRE_IJMatrixSetObjectType(myHypre%A, HYPRE_PARCSR, ierr)
  call HYPRE_IJMatrixInitialize(myHypre%A, ierr)

  ! set up hypre matrix structure
  ! ---------------------------------------------

  ! get matrix structure for hypre partitions

  call HYPRE_IJMatrixSetValues(myHypre%A, myHypre%nrows, myHypre%ncols, myHypre%rows, myHypre%cols, myHypre%values, ierr)

  call HYPRE_IJMatrixAssemble(myHypre%A, ierr)
  call HYPRE_IJMatrixGetObject(myHypre%A, myHypre%parcsr_A, ierr)

  ! set up parallel structure for RHS and solution vector
  ! ---------------------------------------------
  call HYPRE_IJVectorCreate(myHypre%communicator,myHypre%ilower, myHypre%iupper, myHypre%b, ierr)
  call HYPRE_IJVectorSetObjectType(myHypre%b, HYPRE_PARCSR, ierr)
  call HYPRE_IJVectorInitialize(myHypre%b, ierr)

  call HYPRE_IJVectorCreate(myHypre%communicator,myHypre%ilower, myHypre%iupper, myHypre%x, ierr)
  call HYPRE_IJVectorSetObjectType(myHypre%x, HYPRE_PARCSR, ierr)
  call HYPRE_IJVectorInitialize(myHypre%x, ierr)

  ! Setup the PCG Solver
  ! ---------------------------------------------

  call HYPRE_ParCSRPCGCreate(myHypre%communicator, myHypre%solver, ierr)
  call HYPRE_ParCSRPCGSetMaxIter(myHypre%solver, 10, ierr)
  call HYPRE_ParCSRPCGSetTol(myHypre%solver, 1.0d-5, ierr)
  call HYPRE_ParCSRPCGSetPrintLevel(myHypre%solver, 0, ierr)
  call HYPRE_ParCSRPCGSetLogging(myHypre%solver, 0, ierr)

   
  ! setup the preconditioner
  ! ---------------------------------------------
  !        Create solver
  call HYPRE_BoomerAMGCreate(myHypre%precond, ierr)

  !        print solve info + parameters
  ! call HYPRE_BoomerAMGSetPrintLevel(myHypre%precond, 3, ierr)
  !        set 3d problem
  call HYPRE_BoomerAMGSetStrongThrshld(myHypre%precond, 0.25d0, ierr)
#ifdef HYPRE_GPU_AMG
  call HYPRE_BoomerAMGSetCoarsenType(myHypre%precond, 8, ierr)
  call HYPRE_BoomerAMGSetRelaxType(myHypre%precond, 6, ierr)
#else
  !        use HMIS coarsening
  call HYPRE_BoomerAMGSetCoarsenType(myHypre%precond, 10, ierr)
  !        G-S/Jacobi hybrid relaxation
  call HYPRE_BoomerAMGSetRelaxType(myHypre%precond, 8, ierr)
#endif
  !        C/F relaxation
  call HYPRE_BoomerAMGSetRelaxOrder(myHypre%precond, 1, ierr)
  !        Sweeeps on each level
  call HYPRE_BoomerAMGSetNumSweeps(myHypre%precond, 2, ierr)
  !        maximum number of levels
  call HYPRE_BoomerAMGSetMaxLevels(myHypre%precond, 20, ierr)
#ifdef HYPRE_GPU_AMG
  call Hypre_BoomerAMGSetInterpType(myHypre%precond, 6, ierr)
  call HYPRE_BoomerAMGSetTruncFactor(myHypre%precond, 0.25d0, ierr)
#else
  !        set interpolation type
  call Hypre_BoomerAMGSetInterpType(myHypre%precond, 13, ierr)
#endif
  !        Max numbers per rows
  call HYPRE_BoomerAMGSetPMaxElmts(myHypre%precond, 5, ierr)

#ifndef HYPRE_GPU_AMG
  call HYPRE_BoomerAMGSetNumFunctions(myHypre%precond, 4, ierr)
  call HYPRE_BoomerAMGSetNodal(myHypre%precond, 3, ierr)
#endif

  call HYPRE_BoomerAMGSetCycleType(myHypre%precond, 1, ierr)

  call Hypre_BoomerAMGSetMaxIter(myHypre%precond, 1, ierr)

  call HYPRE_BoomerAMGSetTol(myHypre%precond, 0.0d0, ierr)



  !        set amg as the pcg preconditioner, precond_id = 2 -> AMG
  call HYPRE_ParCSRPCGSetPrecond(myHypre%solver, 2, myHypre%precond, ierr)

  end if
end if ! solver is Setup
myHypre%solverIsSet = .true.

#ifdef HYPRE_TIMING
hypre_pcg_call_count = hypre_pcg_call_count + 1
#endif

if (myid.ne.0) then
  !        Now setup and solve!
#ifdef HYPRE_TIMING
  total_start = MPI_WTIME()
#endif
  ! set vector values
  local_size = myHypre%iupper - myHypre%ilower + 1
  call HYPRE_IJVectorSetValues(myHypre%b, local_size, myHypre%rows, myHypre%rhs, ierr)
  call HYPRE_IJVectorSetValues(myHypre%x, local_size, myHypre%rows, myHypre%sol, ierr)

  call HYPRE_IJVectorAssemble(myHypre%b, ierr)
  call HYPRE_IJVectorAssemble(myHypre%x, ierr)

  call HYPRE_IJVectorGetObject(myHypre%b, myHypre%par_b, ierr)
  call HYPRE_IJVectorGetObject(myHypre%x, myHypre%par_x, ierr)

#ifdef HYPRE_TIMING
  setup_start = MPI_WTIME()
#endif
  call HYPRE_ParCSRPCGSetup(myHypre%solver, myHypre%parcsr_A, myHypre%par_b,myHypre%par_x, ierr)
#ifdef HYPRE_TIMING
  setup_end = MPI_WTIME()
  solve_start = MPI_WTIME()
#endif
  call HYPRE_ParCSRPCGSolve(myHypre%solver, myHypre%parcsr_A, myHypre%par_b,myHypre%par_x, ierr)
  call HYPRE_ParCSRPCGGetNumIterations(myHypre%solver, num_iterations, ierr)
  call HYPRE_ParCSRPCGGetFinalRelative(myHypre%solver, final_res_norm, ierr)
#ifdef HYPRE_TIMING
  solve_end = MPI_WTIME()
#endif
end if


if (myid.ne.0) then
 ! Recover the values from HYPRE back to "x"
 call HYPRE_IJVectorGetValues(myHypre%x, local_size, myHypre%rows, myHypre%sol, ierr)
#ifdef HYPRE_TIMING
 total_end = MPI_WTIME()
 call ReportHypreTiming('pcg-amg', hypre_pcg_call_count, setup_end - setup_start, &
                        solve_end - solve_start, total_end - total_start, &
                        num_iterations, final_res_norm)
#endif
end if



  if (myid.eq.1) then
   call MPI_send(num_iterations, 1, MPI_int, 0, 1, MPI_COMM_WORLD, ierr)
  else if (myid.eq.0) then
   call MPI_recv(num_iterations, 1, MPI_int, 1, 1, MPI_COMM_WORLD, status, ierr)
  end if

end subroutine myHyprePCG_Solve


#ifdef HYPRE_TIMING
SUBROUTINE ReportHypreTiming(label, call_index, setup_time, solve_time, total_time, &
                             num_iterations, final_res_norm)
USE var_QuadScalar
IMPLICIT NONE

character(len=*), intent(in) :: label
integer, intent(in) :: call_index, num_iterations
real*8, intent(in) :: setup_time, solve_time, total_time, final_res_norm
integer :: ierr, subrank
real*8 :: local_times(3), max_times(3)

local_times = (/ setup_time, solve_time, total_time /)
call MPI_Comm_rank(myHypre%communicator, subrank, ierr)
call MPI_Reduce(local_times, max_times, 3, MPI_DOUBLE_PRECISION, MPI_MAX, 0, &
                myHypre%communicator, ierr)

if (subrank.eq.0) then
  write(*,'(A,A,A,I0,A,ES14.6,A,ES14.6,A,ES14.6,A,I0,A,ES14.6)') &
    'HYPRE_TIMING solver=', trim(label), ' call=', call_index, &
    ' setup_s=', max_times(1), ' solve_s=', max_times(2), &
    ' total_s=', max_times(3), ' iterations=', num_iterations, &
    ' rel_res=', final_res_norm
end if

END SUBROUTINE ReportHypreTiming
#endif







END MODULE HypreSolver
