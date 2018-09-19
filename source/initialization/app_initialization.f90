module app_initialization

use var_QuadScalar,only:knvt,knet,knat,knel
!-------------------------------------------------------------------------------------------------
! A module for saving the solution values to 
! a file. This output(dump) is mainly done
! to a binary file.
!-------------------------------------------------------------------------------------------------

contains

!========================================================================================
!                              Sub: init_q2p1_app
!========================================================================================
! @param log_unit An integer unit for the log/protocol file
subroutine init_q2p1_app(log_unit)
! Default routine to initialize a q2p1 application 
  USE def_FEAT
  USE PLinScalar, ONLY : Init_PLinScalar,InitCond_PLinLS, &
    UpdateAuxVariables,Transport_PLinLS,Reinitialize_PLinLS, &
    Reinit_Interphase,dMaxSTF
  USE Transport_Q2P1, ONLY : Init_QuadScalar_Stuctures, &
    InitCond_QuadScalar,ProlongateSolution, &
    ResetTimer,bTracer,bViscoElastic,StaticMeshAdaptation,&
    LinScalar_InitCond
  USE ViscoScalar, ONLY : Init_ViscoScalar_Stuctures, &
    Transport_ViscoScalar,IniProf_ViscoScalar,ProlongateViscoSolution
  USE Transport_Q1, ONLY : Init_LinScalar,InitCond_LinScalar, &
    Transport_LinScalar
  USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
  USE var_QuadScalar, ONLY : myStat,cFBM_File,mg_Mesh

  integer, intent(in) :: log_unit

  !-------INIT PHASE-------

  ! Initialization for FEATFLOW
  call General_init_ext(79,log_unit)

  call Init_QuadScalar_Stuctures(log_unit)

  IF(bViscoElastic)call Init_ViscoScalar_Stuctures(log_unit)

  CALL Init_LinScalar(log_unit)

  call InitCond_LinScalar()

  ! Normal start from inital configuration
  if (istart.eq.0) then
    if (myid.ne.0) call CreateDumpStructures(1)
    call InitCond_QuadScalar()
    IF(bViscoElastic)call IniProf_ViscoScalar()

  ! Start from a solution on the same lvl
  ! with the same number of partitions
  elseif (istart.eq.1) then
    if (myid.ne.0) call CreateDumpStructures(1)
    call SolFromFile(CSTART,1)

  ! Start from a solution on a lower lvl
  ! with the same number of partitions
  elseif (istart.eq.2)then
    ! In order to read in from a lower level
    ! the lower level structures are needed
    if (myid.ne.0) call CreateDumpStructures(0)
    call SolFromFile(CSTART,0)
    call ProlongateSolution()

    ! Now generate the structures for the actual level 
    if (myid.ne.0) call CreateDumpStructures(1)

  ! Start from a solution on the same lvl
  ! with a different number of partitions
  elseif (istart.eq.3) then
    IF (myid.ne.0) CALL CreateDumpStructures(1)
    call SolFromFileRepart(CSTART,1)
  end if

end subroutine init_q2p1_app

!========================================================================================
!                             Sub: init_q2p1_particle_tracer
!========================================================================================
! @param log_unit An integer unit for the log/protocol file
subroutine init_q2p1_particle_tracer(log_unit)
! Routine to initialize a particle tracer application 
  USE def_FEAT
  USE PLinScalar, ONLY : Init_PLinScalar,InitCond_PLinLS, &
    UpdateAuxVariables,Transport_PLinLS,Reinitialize_PLinLS, &
    Reinit_Interphase,dMaxSTF
  USE Transport_Q2P1, ONLY : Init_QuadScalar_Stuctures, &
    InitCond_QuadScalar,ProlongateSolution, &
    ResetTimer,bTracer,bViscoElastic,StaticMeshAdaptation,&
    LinScalar_InitCond, Init_QuadScalar, QuadSc 
  USE ViscoScalar, ONLY : Init_ViscoScalar_Stuctures, &
    Transport_ViscoScalar,IniProf_ViscoScalar,ProlongateViscoSolution
  USE Transport_Q1, ONLY : Init_LinScalar,InitCond_LinScalar, &
    Transport_LinScalar
  USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
  USE var_QuadScalar, ONLY : myStat,cFBM_File,mg_Mesh

  integer, intent(in) :: log_unit

  integer :: i

  !-------INIT PHASE-------

  ! Initialization for FEATFLOW
  call General_init_ext(79,log_unit)

  call Init_QuadScalar_Stuctures(log_unit)

  ! The particle tracer by default starts with a single 
  ! combined solution file 
  if (istart.eq.3) then
    IF (myid.ne.0) CALL CreateDumpStructures(1)
    !call SolFromFileRepartPartTracer(CSTART,1)
  else
    if(myid .eq. showid)then
      write(*,*)'Error: A particle tracer application needs setting istart = 3'
    end if
    call myMPI_Barrier()
    stop
  end if

  do i = 1, mg_mesh%level(maxlevel)%NVT 
    mg_mesh%level(maxlevel)%dcorvg(1,i) = QuadSc%auxU(i)
    mg_mesh%level(maxlevel)%dcorvg(2,i) = QuadSc%auxV(i)
    mg_mesh%level(maxlevel)%dcorvg(3,i) = QuadSc%auxW(i)
  end do

end subroutine init_q2p1_particle_tracer

end module app_initialization
