module checkpoint_mpi_prf

use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
use, intrinsic :: iso_fortran_env, only: int32, int64, real64
use checkpoint_types, only: CHECKPOINT_LAYOUT_P1, CHECKPOINT_LAYOUT_Q2, &
  CHECKPOINT_LAYOUT_Q2_COORDINATES

implicit none

integer, parameter :: PRF_METADATA_VERSION = 1
integer, parameter :: PRF_STATUS_VALID = 0
integer, parameter :: PRF_STATUS_LEGACY = 1
integer, parameter :: PRF_STATUS_INCOMPLETE = 2
integer, parameter :: PRF_STATUS_MALFORMED = 3
integer, parameter :: PRF_STATUS_MISSING = 4
integer, parameter :: PRF_STATUS_SMARTDUMP = 5

type t_prf_field_metadata
  character(len=64) :: name = ''
  integer :: layout = 0
  integer :: components = 0
  integer :: chunks_per_component = 0
  integer(int64) :: global_entries_per_component = 0_int64
  integer(int64) :: payload_bytes_per_component = 0_int64
end type t_prf_field_metadata

type t_prf_metadata
  integer :: version = PRF_METADATA_VERSION
  integer(int64) :: generation = 0_int64
  logical :: clock_valid = .false.
  real(real64) :: simulation_time = 0.0_real64
  integer(int64) :: completed_step = 0_int64
  integer :: writer_ranks = 0
  integer :: integer_bytes = 0
  integer :: real_bytes = 0
  character(len=8) :: byte_order = ''
  integer(int64) :: coarse_elements = 0_int64
  integer :: source_level = 0
  type(t_prf_field_metadata), allocatable :: field(:)
end type t_prf_metadata

interface
  integer(c_int) function ff_checkpoint_remove_file(path) &
      bind(C, name='ff_checkpoint_remove_file')
    import :: c_char, c_int
    character(kind=c_char), intent(in) :: path(*)
  end function ff_checkpoint_remove_file

  integer(c_int) function ff_checkpoint_rename_file(old_path, new_path) &
      bind(C, name='ff_checkpoint_rename_file')
    import :: c_char, c_int
    character(kind=c_char), intent(in) :: old_path(*)
    character(kind=c_char), intent(in) :: new_path(*)
  end function ff_checkpoint_rename_file

  integer(c_int) function ff_checkpoint_remove_payloads(directory, prefix) &
      bind(C, name='ff_checkpoint_remove_payloads')
    import :: c_char, c_int
    character(kind=c_char), intent(in) :: directory(*)
    character(kind=c_char), intent(in) :: prefix(*)
  end function ff_checkpoint_remove_payloads

  integer(c_int) function ff_checkpoint_has_mpi_prf_payload(directory) &
      bind(C, name='ff_checkpoint_has_mpi_prf_payload')
    import :: c_char, c_int
    character(kind=c_char), intent(in) :: directory(*)
  end function ff_checkpoint_has_mpi_prf_payload
end interface

private

public :: PRF_STATUS_VALID, PRF_STATUS_LEGACY, PRF_STATUS_INCOMPLETE
public :: PRF_STATUS_MALFORMED, PRF_STATUS_MISSING, PRF_STATUS_SMARTDUMP
public :: t_prf_field_metadata, t_prf_metadata
public :: prf_metadata_initialize, prf_metadata_add_field
public :: prf_metadata_remove_field, prf_metadata_field_chunks
public :: prf_begin_write, prf_cleanup_field, prf_commit_write
public :: prf_read_metadata, prf_status_name, prf_checkpoint_directory

contains

subroutine prf_metadata_initialize(metadata, simulation_time, completed_step, &
    writer_ranks, coarse_elements, source_level, clock_valid)
  type(t_prf_metadata), intent(out) :: metadata
  real(real64), intent(in) :: simulation_time
  integer(int64), intent(in) :: completed_step, coarse_elements
  integer, intent(in) :: writer_ranks, source_level
  logical, intent(in) :: clock_valid
  integer :: clock_count

  call system_clock(clock_count)
  metadata%generation = int(clock_count, int64)
  metadata%clock_valid = clock_valid
  metadata%simulation_time = simulation_time
  metadata%completed_step = completed_step
  metadata%writer_ranks = writer_ranks
  metadata%integer_bytes = storage_size(0)/8
  metadata%real_bytes = storage_size(0.0_real64)/8
  metadata%byte_order = native_byte_order()
  metadata%coarse_elements = coarse_elements
  metadata%source_level = source_level
  allocate(metadata%field(0))
end subroutine prf_metadata_initialize

subroutine prf_metadata_add_field(metadata, name, layout, components, chunks, &
    global_entries, payload_bytes)
  type(t_prf_metadata), intent(inout) :: metadata
  character(len=*), intent(in) :: name
  integer, intent(in) :: layout, components, chunks
  integer(int64), intent(in) :: global_entries, payload_bytes
  type(t_prf_field_metadata), allocatable :: expanded(:)
  integer :: count

  count = size(metadata%field)
  allocate(expanded(count+1))
  if (count > 0) expanded(1:count) = metadata%field
  expanded(count+1)%name = adjustl(trim(name))
  expanded(count+1)%layout = layout
  expanded(count+1)%components = components
  expanded(count+1)%chunks_per_component = chunks
  expanded(count+1)%global_entries_per_component = global_entries
  expanded(count+1)%payload_bytes_per_component = payload_bytes
  call move_alloc(expanded, metadata%field)
end subroutine prf_metadata_add_field

subroutine prf_metadata_remove_field(metadata, name)
  type(t_prf_metadata), intent(inout) :: metadata
  character(len=*), intent(in) :: name
  type(t_prf_field_metadata), allocatable :: retained(:)
  integer :: i, retained_count

  retained_count = count([(trim(metadata%field(i)%name) /= trim(name), &
    i=1,size(metadata%field))])
  allocate(retained(retained_count))
  retained_count = 0
  do i = 1, size(metadata%field)
    if (trim(metadata%field(i)%name) == trim(name)) cycle
    retained_count = retained_count + 1
    retained(retained_count) = metadata%field(i)
  end do
  call move_alloc(retained, metadata%field)
end subroutine prf_metadata_remove_field

integer function prf_metadata_field_chunks(metadata, name, legacy_value)
  type(t_prf_metadata), intent(in) :: metadata
  character(len=*), intent(in) :: name
  integer, intent(in) :: legacy_value
  integer :: i

  prf_metadata_field_chunks = legacy_value
  do i = 1, size(metadata%field)
    if (trim(metadata%field(i)%name) == trim(name)) then
      prf_metadata_field_chunks = metadata%field(i)%chunks_per_component
      return
    end if
  end do
end function prf_metadata_field_chunks

subroutine prf_begin_write(slot, replace_payloads, generation, ierr, message)
  integer, intent(in) :: slot
  logical, intent(in) :: replace_payloads
  integer(int64), intent(in) :: generation
  integer, intent(out) :: ierr
  character(len=*), intent(out) :: message
  character(len=256) :: directory, marker, temporary, sentinel
  integer :: unit_number, io_status

  ierr = 0
  message = ''
  call checkpoint_paths(slot, directory, marker, temporary, sentinel)

  open(newunit=unit_number, file=trim(sentinel), status='replace', &
    action='write', iostat=io_status)
  if (io_status /= 0) then
    ierr = io_status
    write(message,'(A,A)') 'Cannot create checkpoint sentinel: ', trim(sentinel)
    return
  end if
  write(unit_number,'(A,1X,I0)') 'generation', generation
  close(unit_number)

  io_status = ff_checkpoint_remove_file(trim(marker)//c_null_char)
  if (io_status == 0) then
    io_status = ff_checkpoint_remove_file(trim(temporary)//c_null_char)
  end if
  if (io_status == 0 .and. replace_payloads) then
    io_status = ff_checkpoint_remove_payloads( &
      trim(directory)//c_null_char, c_null_char)
  end if
  if (io_status /= 0) then
    ierr = io_status
    write(message,'(A,A)') 'Cannot prepare checkpoint directory: ', &
      trim(directory)
  end if
end subroutine prf_begin_write

subroutine prf_cleanup_field(slot, field_name, ierr, message)
  integer, intent(in) :: slot
  character(len=*), intent(in) :: field_name
  integer, intent(out) :: ierr
  character(len=*), intent(out) :: message
  character(len=256) :: directory

  directory = prf_checkpoint_directory(slot)
  ierr = ff_checkpoint_remove_payloads(trim(directory)//c_null_char, &
    trim(field_name)//c_null_char)
  message = ''
  if (ierr /= 0) then
    write(message,'(A,A)') 'Cannot remove old payloads for field: ', &
      trim(field_name)
  end if
end subroutine prf_cleanup_field

subroutine prf_commit_write(slot, metadata, ierr, message)
  integer, intent(in) :: slot
  type(t_prf_metadata), intent(in) :: metadata
  integer, intent(out) :: ierr
  character(len=*), intent(out) :: message
  character(len=256) :: directory, marker, temporary, sentinel
  integer :: unit_number, io_status, i

  ierr = 0
  message = ''
  call checkpoint_paths(slot, directory, marker, temporary, sentinel)

  open(newunit=unit_number, file=trim(temporary), status='replace', &
    action='write', iostat=io_status)
  if (io_status /= 0) then
    ierr = io_status
    write(message,'(A,A)') 'Cannot write checkpoint metadata: ', trim(temporary)
    return
  end if

  write(unit_number,'(A)') 'FF_MPI_PRF_CHECKPOINT'
  write(unit_number,'(A,1X,I0)') 'version', metadata%version
  write(unit_number,'(A,1X,I0)') 'generation', metadata%generation
  write(unit_number,'(A,1X,I0)') 'clock_valid', merge(1, 0, metadata%clock_valid)
  write(unit_number,'(A,1X,ES25.17E3)') &
    'simulation_time', metadata%simulation_time
  write(unit_number,'(A,1X,I0)') 'completed_step', metadata%completed_step
  write(unit_number,'(A,1X,I0)') 'writer_ranks', metadata%writer_ranks
  write(unit_number,'(A,1X,I0)') 'integer_bytes', metadata%integer_bytes
  write(unit_number,'(A,1X,I0)') 'real_bytes', metadata%real_bytes
  write(unit_number,'(A,1X,A)') 'byte_order', trim(metadata%byte_order)
  write(unit_number,'(A,1X,I0)') 'coarse_elements', metadata%coarse_elements
  write(unit_number,'(A,1X,I0)') 'source_level', metadata%source_level
  write(unit_number,'(A,1X,I0)') 'field_count', size(metadata%field)
  do i = 1, size(metadata%field)
    write(unit_number,'(A,1X,A,1X,3(I0,1X),2(I0,1X))') 'field', &
      trim(metadata%field(i)%name), metadata%field(i)%layout, &
      metadata%field(i)%components, &
      metadata%field(i)%chunks_per_component, &
      metadata%field(i)%global_entries_per_component, &
      metadata%field(i)%payload_bytes_per_component
  end do
  close(unit_number, iostat=io_status)
  if (io_status /= 0) then
    ierr = io_status
    message = 'Cannot close checkpoint metadata temporary file.'
    return
  end if

  io_status = ff_checkpoint_rename_file(trim(temporary)//c_null_char, &
    trim(marker)//c_null_char)
  if (io_status == 0) then
    io_status = ff_checkpoint_remove_file(trim(sentinel)//c_null_char)
  end if
  if (io_status /= 0) then
    ierr = io_status
    message = 'Cannot commit checkpoint metadata marker.'
  end if
end subroutine prf_commit_write

subroutine prf_read_metadata(slot, metadata, status, message)
  integer, intent(in) :: slot
  type(t_prf_metadata), intent(out) :: metadata
  integer, intent(out) :: status
  character(len=*), intent(out) :: message
  character(len=256) :: directory, marker, temporary, sentinel
  character(len=512) :: line
  character(len=64) :: key, field_name
  character(len=8) :: byte_order
  logical :: marker_exists, sentinel_exists, time_exists
  integer :: unit_number, io_status, clock_value, expected_fields
  integer :: layout, components, chunks
  integer(int64) :: global_entries, payload_bytes

  call initialize_empty_metadata(metadata)
  message = ''
  call checkpoint_paths(slot, directory, marker, temporary, sentinel)
  inquire(file=trim(marker), exist=marker_exists)
  inquire(file=trim(sentinel), exist=sentinel_exists)

  if (.not. marker_exists) then
    if (sentinel_exists) then
      status = PRF_STATUS_INCOMPLETE
      message = 'Checkpoint has an incomplete-generation sentinel.'
      return
    end if
    if (ff_checkpoint_has_mpi_prf_payload( &
        trim(directory)//c_null_char) /= 0) then
      status = PRF_STATUS_LEGACY
      message = 'Legacy MPI-PRF checkpoint has no metadata marker.'
      return
    end if
    inquire(file=trim(directory)//'/time.prf', exist=time_exists)
    if (time_exists .and. has_smartdump_header(trim(directory)//'/time.prf')) then
      status = PRF_STATUS_SMARTDUMP
      message = 'SmartDump time.prf is not an MPI-PRF checkpoint.'
      return
    end if
    status = PRF_STATUS_MISSING
    message = 'MPI-PRF checkpoint does not exist.'
    return
  end if

  open(newunit=unit_number, file=trim(marker), status='old', &
    action='read', iostat=io_status)
  if (io_status /= 0) then
    status = PRF_STATUS_MALFORMED
    message = 'Cannot open checkpoint metadata marker.'
    return
  end if

  read(unit_number,'(A)',iostat=io_status) line
  if (io_status /= 0 .or. trim(line) /= 'FF_MPI_PRF_CHECKPOINT') then
    close(unit_number)
    status = PRF_STATUS_MALFORMED
    message = 'Checkpoint metadata has an invalid magic string.'
    return
  end if

  expected_fields = -1
  do
    read(unit_number,'(A)',iostat=io_status) line
    if (io_status < 0) exit
    if (io_status > 0) then
      status = PRF_STATUS_MALFORMED
      message = 'Cannot read checkpoint metadata.'
      close(unit_number)
      return
    end if
    key = ''
    read(line,*,iostat=io_status) key
    if (io_status /= 0) cycle
    select case (trim(key))
    case ('version')
      read(line,*,iostat=io_status) key, metadata%version
    case ('generation')
      read(line,*,iostat=io_status) key, metadata%generation
    case ('clock_valid')
      read(line,*,iostat=io_status) key, clock_value
      metadata%clock_valid = clock_value /= 0
    case ('simulation_time')
      read(line,*,iostat=io_status) key, metadata%simulation_time
    case ('completed_step')
      read(line,*,iostat=io_status) key, metadata%completed_step
    case ('writer_ranks')
      read(line,*,iostat=io_status) key, metadata%writer_ranks
    case ('integer_bytes')
      read(line,*,iostat=io_status) key, metadata%integer_bytes
    case ('real_bytes')
      read(line,*,iostat=io_status) key, metadata%real_bytes
    case ('byte_order')
      read(line,*,iostat=io_status) key, byte_order
      metadata%byte_order = byte_order
    case ('coarse_elements')
      read(line,*,iostat=io_status) key, metadata%coarse_elements
    case ('source_level')
      read(line,*,iostat=io_status) key, metadata%source_level
    case ('field_count')
      read(line,*,iostat=io_status) key, expected_fields
    case ('field')
      read(line,*,iostat=io_status) key, field_name, layout, components, &
        chunks, global_entries, payload_bytes
      if (io_status == 0) call prf_metadata_add_field(metadata, field_name, &
        layout, components, chunks, global_entries, payload_bytes)
    case default
      io_status = 0
    end select
    if (io_status /= 0) then
      status = PRF_STATUS_MALFORMED
      write(message,'(A,A)') 'Malformed checkpoint metadata entry: ', &
        trim(line)
      close(unit_number)
      return
    end if
  end do
  close(unit_number)

  if (metadata%version /= PRF_METADATA_VERSION) then
    status = PRF_STATUS_MALFORMED
    message = 'Unsupported checkpoint metadata version.'
  else if (expected_fields < 0 .or. &
      expected_fields /= size(metadata%field)) then
    status = PRF_STATUS_MALFORMED
    message = 'Checkpoint metadata field count does not match its manifest.'
  else if (.not.manifest_is_valid(metadata)) then
    status = PRF_STATUS_MALFORMED
    message = 'Checkpoint metadata contains an invalid field manifest.'
  else if (metadata%integer_bytes /= storage_size(0)/8 .or. &
      metadata%real_bytes /= storage_size(0.0_real64)/8 .or. &
      trim(metadata%byte_order) /= trim(native_byte_order())) then
    status = PRF_STATUS_MALFORMED
    message = 'Checkpoint native integer/real representation is incompatible.'
  else
    call validate_payload_files(slot, metadata, status, message)
  end if
end subroutine prf_read_metadata

logical function manifest_is_valid(metadata)
  type(t_prf_metadata), intent(in) :: metadata
  integer :: i, j

  manifest_is_valid = .false.
  if (metadata%writer_ranks <= 0 .or. metadata%coarse_elements <= 0_int64 .or. &
      metadata%source_level <= 0) return
  if (size(metadata%field) == 0) return

  do i = 1, size(metadata%field)
    if (len_trim(metadata%field(i)%name) == 0) return
    if (index(metadata%field(i)%name,'/') /= 0 .or. &
        index(metadata%field(i)%name,achar(92)) /= 0) return
    if (metadata%field(i)%layout < CHECKPOINT_LAYOUT_P1 .or. &
        metadata%field(i)%layout > CHECKPOINT_LAYOUT_Q2_COORDINATES) return
    if (metadata%field(i)%components <= 0 .or. &
        metadata%field(i)%chunks_per_component <= 0) return
    if (metadata%field(i)%global_entries_per_component <= 0_int64 .or. &
        metadata%field(i)%payload_bytes_per_component /= &
        metadata%field(i)%global_entries_per_component * &
        int(metadata%real_bytes,int64)) return
    do j = 1, i-1
      if (trim(metadata%field(j)%name) == trim(metadata%field(i)%name)) return
    end do
  end do
  manifest_is_valid = .true.
end function manifest_is_valid

subroutine validate_payload_files(slot, metadata, status, message)
  integer, intent(in) :: slot
  type(t_prf_metadata), intent(in) :: metadata
  integer, intent(out) :: status
  character(len=*), intent(out) :: message
  character(len=256) :: directory, file_name
  logical :: exists
  integer :: i, component, chunk
  integer(int64) :: file_size, component_size

  directory = prf_checkpoint_directory(slot)
  do i = 1, size(metadata%field)
    do component = 1, metadata%field(i)%components
      component_size = 0_int64
      do chunk = 1, metadata%field(i)%chunks_per_component
        write(file_name,'(A,A,A,I0,A,I0,A)') trim(directory), '/', &
          trim(metadata%field(i)%name)//'_comp', component, '_chunk_', &
          chunk, '.prf'
        inquire(file=trim(file_name), exist=exists, size=file_size)
        if (.not. exists) then
          status = PRF_STATUS_MALFORMED
          write(message,'(A,A)') 'Checkpoint payload is missing: ', &
            trim(file_name)
          return
        end if
        component_size = component_size + file_size
      end do
      if (component_size /= &
          metadata%field(i)%payload_bytes_per_component) then
        status = PRF_STATUS_MALFORMED
        write(message,'(A,A)') 'Checkpoint payload size mismatch for field: ', &
          trim(metadata%field(i)%name)
        return
      end if
    end do
  end do
  status = PRF_STATUS_VALID
  message = ''
end subroutine validate_payload_files

subroutine initialize_empty_metadata(metadata)
  type(t_prf_metadata), intent(out) :: metadata

  metadata%version = 0
  metadata%integer_bytes = 0
  metadata%real_bytes = 0
  allocate(metadata%field(0))
end subroutine initialize_empty_metadata

subroutine checkpoint_paths(slot, directory, marker, temporary, sentinel)
  integer, intent(in) :: slot
  character(len=*), intent(out) :: directory, marker, temporary, sentinel

  directory = prf_checkpoint_directory(slot)
  marker = trim(directory)//'/checkpoint_meta.prf'
  temporary = trim(directory)//'/checkpoint_meta.prf.tmp'
  sentinel = trim(directory)//'/checkpoint_incomplete.prf'
end subroutine checkpoint_paths

function prf_checkpoint_directory(slot) result(directory)
  integer, intent(in) :: slot
  character(len=256) :: directory

  write(directory,'(A,I0)') '_dump/', slot
end function prf_checkpoint_directory

logical function has_smartdump_header(path)
  character(len=*), intent(in) :: path
  character(len=32) :: header
  integer :: unit_number, io_status

  has_smartdump_header = .false.
  open(newunit=unit_number, file=trim(path), status='old', &
    action='read', iostat=io_status)
  if (io_status /= 0) return
  read(unit_number,*,iostat=io_status) header
  close(unit_number)
  if (io_status == 0) has_smartdump_header = trim(header) == 'timens'
end function has_smartdump_header

function native_byte_order() result(byte_order)
  character(len=8) :: byte_order
  integer(int32) :: one
  character(len=4) :: bytes

  one = 1_int32
  bytes = transfer(one, bytes)
  if (iachar(bytes(1:1)) == 1) then
    byte_order = 'little'
  else
    byte_order = 'big'
  end if
end function native_byte_order

function prf_status_name(status) result(name)
  integer, intent(in) :: status
  character(len=16) :: name

  select case (status)
  case (PRF_STATUS_VALID)
    name = 'valid'
  case (PRF_STATUS_LEGACY)
    name = 'legacy'
  case (PRF_STATUS_INCOMPLETE)
    name = 'incomplete'
  case (PRF_STATUS_MALFORMED)
    name = 'malformed'
  case (PRF_STATUS_MISSING)
    name = 'missing'
  case (PRF_STATUS_SMARTDUMP)
    name = 'smartdump'
  case default
    name = 'unknown'
  end select
end function prf_status_name

end module checkpoint_mpi_prf
