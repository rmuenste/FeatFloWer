module app_initialization

use var_QuadScalar,only:knvt,knet,knat,knel
USE var_QuadScalar, ONLY : myStat,cFBM_File,mg_Mesh
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
USE Transport_Q2P1, ONLY: QuadSc,ProlongateSolution,OperatorRegenaration 
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
  USE Transport_Q2P1, ONLY : Init_QuadScalar_Stuctures, &
    InitCond_QuadScalar, &
    bTracer,bViscoElastic,StaticMeshAdaptation,&
    LinScalar_InitCond
  USE ViscoScalar, ONLY : Init_ViscoScalar_Stuctures, &
    Transport_ViscoScalar,IniProf_ViscoScalar,ProlongateViscoSolution
  USE Transport_Q1, ONLY : Init_LinScalar,InitCond_LinScalar, &
    Transport_LinScalar
  USE checkpoint_config, ONLY: checkpoint_dispatch_enabled, &
    selected_checkpoint_format, CHECKPOINT_FORMAT_MPI_PRF, &
    checkpoint_slot_from_start_file

  integer, intent(in) :: log_unit
  integer :: checkpoint_slot, checkpoint_error

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
    if (checkpoint_dispatch_enabled .and. &
        selected_checkpoint_format.eq.CHECKPOINT_FORMAT_MPI_PRF) then
      checkpoint_slot = checkpoint_slot_from_start_file(CSTART,checkpoint_error)
      call validate_checkpoint_slot(checkpoint_slot,checkpoint_error,CSTART)
      call init_sol_same_level_mpi(checkpoint_slot)
    else
      call init_sol_same_level(CSTART)
    end if

  ! Start from a solution on a lower lvl
  ! with the same number of partitions
  elseif (istart.eq.2)then
    if (checkpoint_dispatch_enabled .and. &
        selected_checkpoint_format.eq.CHECKPOINT_FORMAT_MPI_PRF) then
      checkpoint_slot = checkpoint_slot_from_start_file(CSTART,checkpoint_error)
      call validate_checkpoint_slot(checkpoint_slot,checkpoint_error,CSTART)
      call init_sol_lower_level_mpi(checkpoint_slot,log_unit)
    else
      call init_sol_lower_level(CSTART)
    end if

  ! Start from a solution on the same lvl
  ! with a different number of partitions
  elseif (istart.eq.3) then
    if (checkpoint_dispatch_enabled .and. &
        selected_checkpoint_format.eq.CHECKPOINT_FORMAT_MPI_PRF) then
      checkpoint_slot = checkpoint_slot_from_start_file(CSTART,checkpoint_error)
      call validate_checkpoint_slot(checkpoint_slot,checkpoint_error,CSTART)
      call init_sol_same_level_mpi(checkpoint_slot)
    else
      call init_sol_repart(CSTART)
    end if

  end if


end subroutine init_q2p1_app

subroutine validate_checkpoint_slot(slot, parse_error, start_file)
  USE PP3D_MPI, ONLY: myid,showid,MPI_COMM_WORLD
  USE var_QuadScalar, ONLY: myErrorCode
  implicit none
  integer, intent(in) :: slot, parse_error
  character(len=*), intent(in) :: start_file
  integer :: mpi_error

  if (parse_error.eq.0 .and. slot.ge.0) return
  if (myid.eq.showid) then
    write(*,'(A,A)') 'Invalid MPI-PRF SimPar@StartFile: ',trim(start_file)
    write(*,'(A)') 'Use a numeric slot or a path such as _dump/1.'
  end if
  call MPI_Abort(MPI_COMM_WORLD,myErrorCode%PARAM_FILE_READ_ERROR,mpi_error)
end subroutine validate_checkpoint_slot

subroutine init_sol_same_level_mpi(slot)
  USE Transport_Q2P1, ONLY: OperatorRegenaration
  USE checkpoint_service, ONLY: checkpoint_mpi_read
  USE checkpoint_types, ONLY: t_checkpoint_context, checkpoint_fields_fc
  implicit none
  integer, intent(in) :: slot
  type(t_checkpoint_context) :: checkpoint_context

  checkpoint_context%slot = slot
  checkpoint_context%clock_valid = .true.
  call checkpoint_mpi_read(checkpoint_context,checkpoint_fields_fc())
  call OperatorRegenaration(1)
  call OperatorRegenaration(2)
  call OperatorRegenaration(3)
end subroutine init_sol_same_level_mpi

subroutine init_sol_lower_level_mpi(slot,log_unit)
  USE Transport_Q2P1, ONLY: ProlongateSolution,OperatorRegenaration
  USE def_QuadScalar, ONLY: ProlongateSingleFieldQ2
  USE var_QuadScalar, ONLY: GenLinScalar, mg_mesh
  USE checkpoint_service, ONLY: checkpoint_mpi_read
  USE checkpoint_types, ONLY: t_checkpoint_context, checkpoint_fields_fc, &
    CHECKPOINT_PROLONGATE
  implicit none
  integer, intent(in) :: slot,log_unit
  integer :: field_index, coarse_dofs, fine_dofs, fine_level
  type(t_checkpoint_context) :: checkpoint_context

  checkpoint_context%slot = slot
  checkpoint_context%restart_mode = CHECKPOINT_PROLONGATE
  checkpoint_context%clock_valid = .true.
  call checkpoint_mpi_read(checkpoint_context,checkpoint_fields_fc(),log_unit)
  call ProlongateSolution()
  if (allocated(GenLinScalar%Fld)) then
    fine_level = mg_mesh%maxlevel
    coarse_dofs = KNVT(fine_level-1) + KNAT(fine_level-1) + &
      KNET(fine_level-1) + KNEL(fine_level-1)
    fine_dofs = KNVT(fine_level) + KNAT(fine_level) + &
      KNET(fine_level) + KNEL(fine_level)
    do field_index=1,GenLinScalar%nOfFields
      call ProlongateSingleFieldQ2(GenLinScalar%Fld(field_index)%Val, &
        coarse_dofs,fine_dofs,adjustl(trim(GenLinScalar%Fld(field_index)%cName)))
    end do
  end if
  call OperatorRegenaration(1)
  call OperatorRegenaration(2)
  call OperatorRegenaration(3)
end subroutine init_sol_lower_level_mpi

!========================================================================================
!                             Sub: init_q2p1_particle_tracer
!========================================================================================
! @param log_unit An integer unit for the log/protocol file
subroutine init_q2p1_particle_tracer(log_unit)
! Routine to initialize a particle tracer application 
  USE def_FEAT
  USE Transport_Q2P1, ONLY : Init_QuadScalar_ReducedStuctures, &
    InitCond_QuadScalar,ProlongateSolution, &
    bTracer,bViscoElastic,StaticMeshAdaptation,&
    LinScalar_InitCond, Init_QuadScalar 
  USE ViscoScalar, ONLY : Init_ViscoScalar_Stuctures, &
    Transport_ViscoScalar,IniProf_ViscoScalar,ProlongateViscoSolution
  USE Transport_Q1, ONLY : Init_LinScalar,InitCond_LinScalar, &
    Transport_LinScalar

  integer, intent(in) :: log_unit

  integer :: i

  !-------INIT PHASE-------

  ! Initialization for FEATFLOW
  call General_init_ext(79,log_unit)

!   call Init_QuadScalar_Stuctures(log_unit)

  call Init_QuadScalar_ReducedStuctures(log_unit)
  
  if (myid.ne.0) call CreateDumpStructures(1)
  
  ! The particle tracer by default starts with a single 
  ! combined solution file 
  ! Normal start from inital configuration
!   if (istart.eq.0) then
!   
!    if (myid.ne.0) call CreateDumpStructures(1)
! 
!   ! Start from a solution on the same lvl
!   ! with the same number of partitions
!   elseif (istart.eq.1) then
! 
!   call init_sol_same_level(CSTART)
! 
!   ! Start from a solution on a lower lvl
!   ! with the same number of partitions
!   elseif (istart.eq.2)then
! 
!     call init_sol_lower_level(CSTART)
! 
!   ! Start from a solution on the same lvl
!   ! with a different number of partitions
!   elseif (istart.eq.3) then
! 
!     call init_sol_repart(CSTART)
! 
!   end if

end subroutine init_q2p1_particle_tracer
!========================================================================================
!                             Sub: init_sol_same_level_heat
!========================================================================================
subroutine init_sol_same_level_heat(start_file)
implicit none
character(len=*), intent(in) :: start_file


! Locals
integer :: i, ilev

if (myid.ne.0) call CreateDumpStructures(1)

call SolFromFile_heat(start_file,1)

if (myid .ne. 0) then
  do i = 1, mg_mesh%level(mg_Mesh%maxlevel)%NVT 

    mg_mesh%level(mg_Mesh%maxlevel)%dcorvg(1,i) = QuadSc%auxU(i)
    mg_mesh%level(mg_Mesh%maxlevel)%dcorvg(2,i) = QuadSc%auxV(i)
    mg_mesh%level(mg_Mesh%maxlevel)%dcorvg(3,i) = QuadSc%auxW(i)

  end do
end if

ilev = mg_Mesh%nlmin

call ExchangeNodeValuesOnCoarseLevel(&
  mg_mesh%level(ilev)%dcorvg,&
  mg_mesh%level(ilev)%kvert,&
  mg_mesh%level(ilev)%nvt,&
  mg_mesh%level(ilev)%nel)

! call OperatorRegenaration(1)
! call OperatorRegenaration(2)
! call OperatorRegenaration(3)

end subroutine init_sol_same_level_heat
!========================================================================================
!                             Sub: init_sol_same_level
!========================================================================================
subroutine init_sol_same_level(start_file)
use prov_dump_config, only: use_prov_dump_io
implicit none
character(len=*), intent(in) :: start_file


! Locals
integer :: i, ilev

if (myid.ne.0) call CreateDumpStructures(1)

if (use_prov_dump_io) then
  call SolFromFileProv(start_file,1)
else
  call SolFromFile(start_file,1)
end if

if (myid .ne. 0) then
  do i = 1, mg_mesh%level(mg_Mesh%maxlevel)%NVT 

    mg_mesh%level(mg_Mesh%maxlevel)%dcorvg(1,i) = QuadSc%auxU(i)
    mg_mesh%level(mg_Mesh%maxlevel)%dcorvg(2,i) = QuadSc%auxV(i)
    mg_mesh%level(mg_Mesh%maxlevel)%dcorvg(3,i) = QuadSc%auxW(i)

  end do
end if

ilev = mg_Mesh%nlmin

call ExchangeNodeValuesOnCoarseLevel(&
  mg_mesh%level(ilev)%dcorvg,&
  mg_mesh%level(ilev)%kvert,&
  mg_mesh%level(ilev)%nvt,&
  mg_mesh%level(ilev)%nel)

call OperatorRegenaration(1)
call OperatorRegenaration(2)
call OperatorRegenaration(3)

end subroutine init_sol_same_level
!========================================================================================
!                             Sub: init_sol_lower_level
!========================================================================================
subroutine init_sol_lower_level(start_file)
USE Parametrization, ONLY: ParametrizeBndr
use var_QuadScalar,only:ilev,nlmax
use prov_dump_config, only: use_prov_dump_io
implicit none

character(len=*), intent(in) :: start_file

! Locals
integer :: i

! In order to read in from a lower level
! the lower level structures are needed
if (myid.ne.0) call CreateDumpStructures(0)

if (use_prov_dump_io) then
  call SolFromFileProv(start_file,0)
else
  call SolFromFile(start_file,0)
end if


if (myid .ne. 0) then
  do i = 1, mg_mesh%level(mg_Mesh%maxlevel-1)%NVT 

    mg_mesh%level(mg_Mesh%maxlevel-1)%dcorvg(1,i) = QuadSc%auxU(i)
    mg_mesh%level(mg_Mesh%maxlevel-1)%dcorvg(2,i) = QuadSc%auxV(i)
    mg_mesh%level(mg_Mesh%maxlevel-1)%dcorvg(3,i) = QuadSc%auxW(i)

    if (abs(mg_mesh%level(mg_Mesh%maxlevel)%dcorvg(1,i) - QuadSc%auxU(i)) > 1.0E-5) then
      write(*,*)"myid: ", myid
      write(*,*)"idx: ", i
      write(*,*)"computed: " , mg_mesh%level(mg_Mesh%maxlevel-1)%dcorvg(1,i)
      write(*,*)"read: ",QuadSc%auxU(i)
    end if

  end do
end if

call ProlongateCoordinates(&
  mg_mesh%level(mg_Mesh%maxlevel-1)%dcorvg,&
  mg_mesh%level(mg_Mesh%maxlevel)%dcorvg,&
  mg_mesh%level(mg_Mesh%maxlevel-1)%karea,&
  mg_mesh%level(mg_Mesh%maxlevel-1)%kvert,&
  mg_mesh%level(mg_Mesh%maxlevel-1)%kedge,&
  mg_mesh%level(mg_Mesh%maxlevel-1)%nel,&
  mg_mesh%level(mg_Mesh%maxlevel-1)%nvt,&
  mg_mesh%level(mg_Mesh%maxlevel-1)%net,&
  mg_mesh%level(mg_Mesh%maxlevel-1)%nat)

call ExchangeNodeValuesOnCoarseLevel(&
  mg_mesh%level(1)%dcorvg,&
  mg_mesh%level(1)%kvert,&
  mg_mesh%level(1)%nvt,&
  mg_mesh%level(1)%nel)

call ProlongateSolution()

if(myid.ne.0)then
  NLMAX = NLMAX + 1
  ILEV=NLMAX
  CALL SETLEV(2)
  CALL ParametrizeBndr(mg_mesh,ILEV)
  NLMAX = NLMAX - 1
endif
    
call OperatorRegenaration(1)
call OperatorRegenaration(2)
call OperatorRegenaration(3)

! Now generate the structures for the actual level 
if (myid.ne.0) call CreateDumpStructures(1)


end subroutine init_sol_lower_level
!========================================================================================
!                             Sub: init_sol_repart
!========================================================================================
subroutine init_sol_repart(start_file)
use prov_dump_config, only: use_prov_dump_io
implicit none

character(len=*), intent(in) :: start_file

! Locals
integer :: i,ilev

IF (myid.ne.0) CALL CreateDumpStructures(1)

if (use_prov_dump_io) then
  call SolFromFileRepartProv(start_file,1)
else
  call SolFromFileRepart(start_file,1)
end if

if (myid .ne. 0) then
 do i = 1, mg_mesh%level(mg_Mesh%maxlevel)%NVT 

   mg_mesh%level(mg_Mesh%maxlevel)%dcorvg(1,i) = QuadSc%auxU(i)
   mg_mesh%level(mg_Mesh%maxlevel)%dcorvg(2,i) = QuadSc%auxV(i)
   mg_mesh%level(mg_Mesh%maxlevel)%dcorvg(3,i) = QuadSc%auxW(i)

 end do
end if

ilev = mg_Mesh%nlmin

call ExchangeNodeValuesOnCoarseLevel(&
  mg_mesh%level(ilev)%dcorvg,&
  mg_mesh%level(ilev)%kvert,&
  mg_mesh%level(ilev)%nvt,&
  mg_mesh%level(ilev)%nel)

call OperatorRegenaration(1)
call OperatorRegenaration(2)
call OperatorRegenaration(3)

end subroutine init_sol_repart

end module app_initialization
