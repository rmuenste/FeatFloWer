module checkpoint_types

use, intrinsic :: iso_fortran_env, only: int64, real64

implicit none

integer, parameter :: CHECKPOINT_LAYOUT_P1 = 1
integer, parameter :: CHECKPOINT_LAYOUT_Q2 = 2
integer, parameter :: CHECKPOINT_LAYOUT_Q2_COORDINATES = 3

integer, parameter :: CHECKPOINT_REPLACE = 1
integer, parameter :: CHECKPOINT_MERGE_FIELDS = 2

integer, parameter :: CHECKPOINT_SAME_LEVEL = 1
integer, parameter :: CHECKPOINT_PROLONGATE = 2
integer, parameter :: CHECKPOINT_REPARTITION = 3

type t_checkpoint_component_view
  real(real64), pointer :: values(:) => null()
end type t_checkpoint_component_view

type t_checkpoint_field
  character(len=64) :: name = ''
  integer :: layout = 0
  integer :: component_count = 0
  logical :: required = .true.
  type(t_checkpoint_component_view), allocatable :: component(:)
end type t_checkpoint_field

type t_checkpoint_context
  integer :: slot = 0
  integer :: communicator = 0
  integer :: source_level = 0
  integer :: target_level = 0
  integer(int64) :: coarse_element_count = 0_int64
  real(real64) :: simulation_time = 0.0_real64
  integer(int64) :: completed_step = 0_int64
  logical :: clock_valid = .false.
  integer :: restart_mode = CHECKPOINT_SAME_LEVEL
  integer :: write_mode = CHECKPOINT_REPLACE
end type t_checkpoint_context

type t_checkpoint_field_selection
  logical :: pressure = .false.
  logical :: velocity = .false.
  logical :: coordinates = .false.
  logical :: distance = .false.
  logical :: segment = .false.
  logical :: shell = .false.
  logical :: mixer_boundary = .false.
  logical :: temperature = .false.
  logical :: generic_scalars = .false.
end type t_checkpoint_field_selection

contains

subroutine checkpoint_validate_fields(fields, ierr, message)
  type(t_checkpoint_field), intent(in) :: fields(:)
  integer, intent(out) :: ierr
  character(len=*), intent(out) :: message
  integer :: i, component

  ierr = 0
  message = ''

  do i = 1, size(fields)
    if (len_trim(fields(i)%name) == 0) then
      ierr = 1
      message = 'Checkpoint field has no name.'
      return
    end if
    if (fields(i)%layout < CHECKPOINT_LAYOUT_P1 .or. &
        fields(i)%layout > CHECKPOINT_LAYOUT_Q2_COORDINATES) then
      ierr = 1
      write(message,'(A,A)') 'Checkpoint field has invalid layout: ', &
        trim(fields(i)%name)
      return
    end if
    if (.not. allocated(fields(i)%component)) then
      if (fields(i)%required) then
        ierr = 1
        write(message,'(A,A)') 'Required checkpoint field is unallocated: ', &
          trim(fields(i)%name)
        return
      end if
      cycle
    end if
    if (size(fields(i)%component) /= fields(i)%component_count) then
      ierr = 1
      write(message,'(A,A)') 'Checkpoint component count mismatch: ', &
        trim(fields(i)%name)
      return
    end if
    do component = 1, fields(i)%component_count
      if (.not. associated(fields(i)%component(component)%values) .and. &
          fields(i)%required) then
        ierr = 1
        write(message,'(A,A)') 'Required checkpoint field is unassociated: ', &
          trim(fields(i)%name)
        return
      end if
    end do
  end do
end subroutine checkpoint_validate_fields

function checkpoint_fields_fc() result(fields)
  type(t_checkpoint_field_selection) :: fields

  fields%pressure = .true.
  fields%velocity = .true.
  fields%coordinates = .true.
  fields%temperature = .true.
  fields%generic_scalars = .true.
end function checkpoint_fields_fc

function checkpoint_fields_sse() result(fields)
  type(t_checkpoint_field_selection) :: fields

  fields%pressure = .true.
  fields%velocity = .true.
  fields%coordinates = .true.
  fields%distance = .true.
  fields%segment = .true.
  fields%shell = .true.
  fields%mixer_boundary = .true.
  fields%temperature = .true.
  fields%generic_scalars = .true.
end function checkpoint_fields_sse

function checkpoint_fields_temperature() result(fields)
  type(t_checkpoint_field_selection) :: fields

  fields%temperature = .true.
end function checkpoint_fields_temperature

function checkpoint_field_codes(fields) result(codes)
  type(t_checkpoint_field_selection), intent(in) :: fields
  character(len=32) :: codes

  codes = ''
  if (fields%pressure) call append_code(codes,'p')
  if (fields%velocity) call append_code(codes,'v')
  if (fields%coordinates) call append_code(codes,'x')
  if (fields%distance) call append_code(codes,'d')
  if (fields%segment) call append_code(codes,'s')
  if (fields%shell) call append_code(codes,'y')
  if (fields%mixer_boundary) call append_code(codes,'z')
  if (fields%temperature) call append_code(codes,'t')
  if (fields%generic_scalars) call append_code(codes,'q')
end function checkpoint_field_codes

subroutine append_code(codes,code)
  character(len=*), intent(inout) :: codes
  character(len=1), intent(in) :: code

  if (len_trim(codes) > 0) codes = trim(codes)//','
  codes = trim(codes)//code
end subroutine append_code

end module checkpoint_types
