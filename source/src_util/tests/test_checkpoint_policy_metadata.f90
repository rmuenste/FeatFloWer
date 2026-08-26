program test_checkpoint_policy_metadata

use, intrinsic :: iso_fortran_env, only: int64, real64
use checkpoint_config
use checkpoint_mpi_prf
use checkpoint_types, only: CHECKPOINT_LAYOUT_Q2

implicit none

integer :: ierr, status, unit_number
character(len=256) :: message
logical :: exists
real(real64) :: payload(4)
type(t_prf_metadata) :: written, loaded

call checkpoint_register_application('test', &
  CHECKPOINT_MASK_MPI_PRF + CHECKPOINT_MASK_DMP)

call checkpoint_reset_selection()
call checkpoint_finalize_selection(ierr,message)
call assert_true(ierr == 0, 'default policy is valid')
call assert_true(selected_checkpoint_format == CHECKPOINT_FORMAT_MPI_PRF, &
  'default policy selects MPI_PRF')

call checkpoint_reset_selection()
call checkpoint_set_legacy_format(.false.)
call checkpoint_finalize_selection(ierr,message)
call assert_true(ierr == 0, 'legacy DMP policy is valid')
call assert_true(selected_checkpoint_format == CHECKPOINT_FORMAT_DMP, &
  'UseProvDump=No maps to DMP')

call checkpoint_reset_selection()
call checkpoint_set_format_token('MPI',ierr,message)
call assert_true(ierr == 0, 'MPI alias parses')
call checkpoint_finalize_selection(ierr,message)
call assert_true(selected_checkpoint_format == CHECKPOINT_FORMAT_MPI_PRF, &
  'MPI alias maps to MPI_PRF')

call checkpoint_reset_selection()
call checkpoint_set_format_token('MPI_PRF',ierr,message)
call checkpoint_set_legacy_format(.false.)
call checkpoint_finalize_selection(ierr,message)
call assert_true(ierr /= 0, 'conflicting old and new parameters fail')

payload = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
call prf_metadata_initialize(written,1.25_real64,42_int64,4,1_int64,3,.true.)
written%generation = 1001_int64
call prf_metadata_add_field(written,'test_field',CHECKPOINT_LAYOUT_Q2, &
  1,1,4_int64,32_int64)
call prf_begin_write(1,.true.,written%generation,ierr,message)
call assert_true(ierr == 0, 'transaction begins')
open(newunit=unit_number,file='_dump/1/test_field_comp1_chunk_1.prf', &
  status='replace',access='stream',form='unformatted',action='write')
write(unit_number) payload
close(unit_number)
call prf_commit_write(1,written,ierr,message)
call assert_true(ierr == 0, 'metadata commits')
call prf_read_metadata(1,loaded,status,message)
call assert_true(status == PRF_STATUS_VALID, 'committed metadata validates')
call assert_true(loaded%completed_step == 42_int64, &
  'completed step round-trips')
call assert_true(loaded%simulation_time == 1.25_real64, &
  'simulation time round-trips bitwise')

open(newunit=unit_number,file='_dump/2/time.prf',status='replace', &
  action='write')
write(unit_number,'(A)') 'timens'
write(unit_number,'(D16.8)') 1.25_real64
close(unit_number)
call prf_read_metadata(2,loaded,status,message)
call assert_true(status == PRF_STATUS_SMARTDUMP, &
  'legacy SmartDump time.prf is recognized')

open(newunit=unit_number,file='_dump/3/legacy_key.prf',status='replace', &
  access='stream',form='unformatted',action='write')
close(unit_number)
call prf_read_metadata(3,loaded,status,message)
call assert_true(status == PRF_STATUS_LEGACY, &
  'marker-less MPI_PRF is accepted as legacy')

open(newunit=unit_number,file='_dump/4/checkpoint_incomplete.prf', &
  status='replace',action='write')
write(unit_number,'(A)') 'generation 1'
close(unit_number)
call prf_read_metadata(4,loaded,status,message)
call assert_true(status == PRF_STATUS_INCOMPLETE, &
  'incomplete transaction is rejected')

open(newunit=unit_number,file='_dump/5/velocity_comp1_chunk_9.prf', &
  status='replace',access='stream',form='unformatted',action='write')
write(unit_number) payload
close(unit_number)
call prf_begin_write(5,.true.,5_int64,ierr,message)
call assert_true(ierr == 0, 'replacement transaction begins')
inquire(file='_dump/5/velocity_comp1_chunk_9.prf',exist=exists)
call assert_true(.not.exists, 'replacement removes orphan payload files')

call prf_metadata_initialize(written,2.5_real64,84_int64,4,1_int64,3,.true.)
call prf_metadata_add_field(written,'short_field',CHECKPOINT_LAYOUT_Q2, &
  1,1,4_int64,32_int64)
call prf_begin_write(6,.true.,6_int64,ierr,message)
open(newunit=unit_number,file='_dump/6/short_field_comp1_chunk_1.prf', &
  status='replace',access='stream',form='unformatted',action='write')
write(unit_number) payload(1:2)
close(unit_number)
call prf_commit_write(6,written,ierr,message)
call assert_true(ierr == 0, 'truncated-payload metadata commits')
call prf_read_metadata(6,loaded,status,message)
call assert_true(status == PRF_STATUS_MALFORMED, &
  'truncated payload is rejected before MPI-IO')

open(newunit=unit_number,file='_dump/7/alpha1_comp1_chunk_1.prf', &
  status='replace',access='stream',form='unformatted',action='write')
close(unit_number)
open(newunit=unit_number,file='_dump/7/alpha10_comp1_chunk_1.prf', &
  status='replace',access='stream',form='unformatted',action='write')
close(unit_number)
call prf_cleanup_field(7,'alpha1',ierr,message)
call assert_true(ierr == 0, 'field cleanup succeeds')
inquire(file='_dump/7/alpha1_comp1_chunk_1.prf',exist=exists)
call assert_true(.not.exists, 'field cleanup removes the selected field')
inquire(file='_dump/7/alpha10_comp1_chunk_1.prf',exist=exists)
call assert_true(exists, 'field cleanup preserves prefix-related fields')

contains

subroutine assert_true(condition,description)
  logical, intent(in) :: condition
  character(len=*), intent(in) :: description

  if (.not.condition) then
    write(*,'(A,A)') 'FAILED: ',trim(description)
    error stop 1
  end if
end subroutine assert_true

end program test_checkpoint_policy_metadata
