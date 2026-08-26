module checkpoint_service

use checkpoint_types, only: t_checkpoint_context, &
  t_checkpoint_field_selection, checkpoint_field_codes, CHECKPOINT_REPLACE, &
  CHECKPOINT_MERGE_FIELDS, CHECKPOINT_SAME_LEVEL, CHECKPOINT_PROLONGATE, &
  CHECKPOINT_REPARTITION

implicit none

private
public :: checkpoint_mpi_write, checkpoint_mpi_read

contains

subroutine checkpoint_mpi_write(context,fields)
  type(t_checkpoint_context), intent(in) :: context
  type(t_checkpoint_field_selection), intent(in) :: fields
  character(len=32) :: codes

  codes = checkpoint_field_codes(fields)
  if (len_trim(codes) == 0) return

  select case (context%write_mode)
  case (CHECKPOINT_REPLACE)
    call ReleaseMPIDumpFilesAtStep(context%slot,trim(codes), &
      int(context%completed_step))
  case (CHECKPOINT_MERGE_FIELDS)
    call ReleaseMPIDumpFilesMerge(context%slot,trim(codes))
  case default
    error stop 'Invalid checkpoint write mode.'
  end select
end subroutine checkpoint_mpi_write

subroutine checkpoint_mpi_read(context,fields,log_unit)
  type(t_checkpoint_context), intent(in) :: context
  type(t_checkpoint_field_selection), intent(in) :: fields
  integer, intent(in), optional :: log_unit
  character(len=32) :: codes

  codes = checkpoint_field_codes(fields)
  if (len_trim(codes) == 0) return

  select case (context%restart_mode)
  case (CHECKPOINT_SAME_LEVEL,CHECKPOINT_REPARTITION)
    if (context%clock_valid) then
      call LoadMPIDumpFilesRestoreClock(context%slot,trim(codes))
    else
      call LoadMPIDumpFiles(context%slot,trim(codes))
    end if
  case (CHECKPOINT_PROLONGATE)
    if (.not.present(log_unit)) &
      error stop 'Checkpoint prolongation requires a log unit.'
    if (context%clock_valid) then
      call LoadMPIDumpFilesProlongateRestoreClock(context%slot,trim(codes), &
        log_unit)
    else
      call LoadMPIDumpFilesProlongateSSE(context%slot,trim(codes),log_unit)
    end if
  case default
    error stop 'Invalid checkpoint restart mode.'
  end select
end subroutine checkpoint_mpi_read

end module checkpoint_service
