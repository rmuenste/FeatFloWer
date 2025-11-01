subroutine werr(code, routine)
  implicit none
  integer, intent(in) :: code
  character(len=*), intent(in) :: routine

  write(*, '(A,I0,2A)') 'WERR stub invoked: code=', code, ' routine=', trim(routine)
  stop 'FEAT error encountered in sandbox'
end subroutine werr

subroutine otrc(name, stamp)
  implicit none
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: stamp

  write(*, '(A,2A)') 'OTRC stub: ', trim(name), ' ', trim(stamp)
end subroutine otrc

