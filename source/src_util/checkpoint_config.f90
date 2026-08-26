module checkpoint_config

implicit none

integer, parameter :: CHECKPOINT_FORMAT_MPI_PRF   = 1
integer, parameter :: CHECKPOINT_FORMAT_DMP       = 2
integer, parameter :: CHECKPOINT_FORMAT_LST       = 3
integer, parameter :: CHECKPOINT_FORMAT_PROVENANCE = 4

integer, parameter :: CHECKPOINT_MASK_MPI_PRF = 1
integer, parameter :: CHECKPOINT_MASK_DMP = 2
integer, parameter :: CHECKPOINT_MASK_LST = 4
integer, parameter :: CHECKPOINT_MASK_PROVENANCE = 8

integer, parameter :: CHECKPOINT_SOURCE_DEFAULT = 1
integer, parameter :: CHECKPOINT_SOURCE_EXPLICIT = 2
integer, parameter :: CHECKPOINT_SOURCE_LEGACY = 3

integer, save :: selected_checkpoint_format = CHECKPOINT_FORMAT_MPI_PRF
integer, save :: checkpoint_selection_source = CHECKPOINT_SOURCE_DEFAULT
integer, save :: registered_format_mask = 0
logical, save :: checkpoint_dispatch_enabled = .false.

logical, save :: explicit_format_seen = .false.
logical, save :: legacy_format_seen = .false.
integer, save :: explicit_format = CHECKPOINT_FORMAT_MPI_PRF
integer, save :: legacy_format = CHECKPOINT_FORMAT_DMP
character(len=32), save :: registered_application = 'unregistered'

private

public :: CHECKPOINT_FORMAT_MPI_PRF, CHECKPOINT_FORMAT_DMP
public :: CHECKPOINT_FORMAT_LST, CHECKPOINT_FORMAT_PROVENANCE
public :: CHECKPOINT_MASK_MPI_PRF, CHECKPOINT_MASK_DMP
public :: CHECKPOINT_MASK_LST, CHECKPOINT_MASK_PROVENANCE
public :: CHECKPOINT_SOURCE_DEFAULT, CHECKPOINT_SOURCE_EXPLICIT
public :: CHECKPOINT_SOURCE_LEGACY
public :: selected_checkpoint_format, checkpoint_selection_source
public :: checkpoint_dispatch_enabled
public :: checkpoint_register_application, checkpoint_reset_selection
public :: checkpoint_set_format_token, checkpoint_set_legacy_format
public :: checkpoint_finalize_selection, checkpoint_format_name
public :: checkpoint_source_name, checkpoint_format_is_supported
public :: checkpoint_slot_from_start_file

contains

subroutine checkpoint_register_application(name, supported_mask)
  character(len=*), intent(in) :: name
  integer, intent(in) :: supported_mask

  registered_application = adjustl(trim(name))
  registered_format_mask = supported_mask
  checkpoint_dispatch_enabled = .true.
end subroutine checkpoint_register_application

subroutine checkpoint_reset_selection()
  selected_checkpoint_format = CHECKPOINT_FORMAT_MPI_PRF
  checkpoint_selection_source = CHECKPOINT_SOURCE_DEFAULT
  explicit_format_seen = .false.
  legacy_format_seen = .false.
  explicit_format = CHECKPOINT_FORMAT_MPI_PRF
  legacy_format = CHECKPOINT_FORMAT_DMP
end subroutine checkpoint_reset_selection

subroutine checkpoint_set_format_token(token, ierr, message)
  character(len=*), intent(in) :: token
  integer, intent(out) :: ierr
  character(len=*), intent(out) :: message
  character(len=32) :: normalized

  normalized = uppercase(adjustl(trim(token)))
  ierr = 0
  message = ''

  select case (trim(normalized))
  case ('MPI_PRF', 'MPI')
    explicit_format = CHECKPOINT_FORMAT_MPI_PRF
  case ('DMP')
    explicit_format = CHECKPOINT_FORMAT_DMP
  case ('LST')
    explicit_format = CHECKPOINT_FORMAT_LST
  case ('PROVENANCE')
    explicit_format = CHECKPOINT_FORMAT_PROVENANCE
  case default
    ierr = 1
    write(message,'(A,A,A)') 'Unknown DumpFormat "', trim(token), &
      '"; expected MPI_PRF, DMP, LST, or PROVENANCE.'
    return
  end select

  explicit_format_seen = .true.
end subroutine checkpoint_set_format_token

subroutine checkpoint_set_legacy_format(use_provenance)
  logical, intent(in) :: use_provenance

  legacy_format_seen = .true.
  if (use_provenance) then
    legacy_format = CHECKPOINT_FORMAT_PROVENANCE
  else
    legacy_format = CHECKPOINT_FORMAT_DMP
  end if
end subroutine checkpoint_set_legacy_format

subroutine checkpoint_finalize_selection(ierr, message)
  integer, intent(out) :: ierr
  character(len=*), intent(out) :: message

  ierr = 0
  message = ''

  if (explicit_format_seen) then
    selected_checkpoint_format = explicit_format
    checkpoint_selection_source = CHECKPOINT_SOURCE_EXPLICIT
    if (legacy_format_seen .and. legacy_format /= explicit_format) then
      ierr = 1
      write(message,'(A,A,A,A,A)') 'DumpFormat=', &
        trim(checkpoint_format_name(explicit_format)), &
        ' conflicts with legacy UseProvDump=', &
        trim(checkpoint_format_name(legacy_format)), '.'
      return
    end if
  else if (legacy_format_seen) then
    selected_checkpoint_format = legacy_format
    checkpoint_selection_source = CHECKPOINT_SOURCE_LEGACY
  else
    selected_checkpoint_format = CHECKPOINT_FORMAT_MPI_PRF
    checkpoint_selection_source = CHECKPOINT_SOURCE_DEFAULT
  end if

  if (checkpoint_dispatch_enabled .and. &
      .not. checkpoint_format_is_supported(selected_checkpoint_format)) then
    ierr = 1
    write(message,'(A,A,A,A,A)') 'DumpFormat ', &
      trim(checkpoint_format_name(selected_checkpoint_format)), &
      ' is not supported by application ', trim(registered_application), '.'
  end if
end subroutine checkpoint_finalize_selection

logical function checkpoint_format_is_supported(format_id)
  integer, intent(in) :: format_id
  integer :: format_mask

  if (format_id < CHECKPOINT_FORMAT_MPI_PRF .or. &
      format_id > CHECKPOINT_FORMAT_PROVENANCE) then
    checkpoint_format_is_supported = .false.
    return
  end if

  format_mask = 2**(format_id-1)
  checkpoint_format_is_supported = iand(registered_format_mask, format_mask) /= 0
end function checkpoint_format_is_supported

function checkpoint_format_name(format_id) result(name)
  integer, intent(in) :: format_id
  character(len=16) :: name

  select case (format_id)
  case (CHECKPOINT_FORMAT_MPI_PRF)
    name = 'MPI_PRF'
  case (CHECKPOINT_FORMAT_DMP)
    name = 'DMP'
  case (CHECKPOINT_FORMAT_LST)
    name = 'LST'
  case (CHECKPOINT_FORMAT_PROVENANCE)
    name = 'PROVENANCE'
  case default
    name = 'UNKNOWN'
  end select
end function checkpoint_format_name

function checkpoint_source_name(source_id) result(name)
  integer, intent(in) :: source_id
  character(len=16) :: name

  select case (source_id)
  case (CHECKPOINT_SOURCE_DEFAULT)
    name = 'default'
  case (CHECKPOINT_SOURCE_EXPLICIT)
    name = 'explicit'
  case (CHECKPOINT_SOURCE_LEGACY)
    name = 'legacy'
  case default
    name = 'unknown'
  end select
end function checkpoint_source_name

integer function checkpoint_slot_from_start_file(start_file, ierr)
  character(len=*), intent(in) :: start_file
  integer, intent(out) :: ierr
  character(len=256) :: value
  integer :: slash_pos

  value = adjustl(trim(start_file))
  slash_pos = index(trim(value), '/', back=.true.)
  if (slash_pos > 0) value = value(slash_pos+1:)

  read(value,*,iostat=ierr) checkpoint_slot_from_start_file
  if (ierr /= 0) checkpoint_slot_from_start_file = -1
end function checkpoint_slot_from_start_file

function uppercase(value) result(result_value)
  character(len=*), intent(in) :: value
  character(len=len(value)) :: result_value
  integer :: i, code

  result_value = value
  do i = 1, len(result_value)
    code = iachar(result_value(i:i))
    if (code >= iachar('a') .and. code <= iachar('z')) then
      result_value(i:i) = achar(code - iachar('a') + iachar('A'))
    end if
  end do
end function uppercase

end module checkpoint_config
