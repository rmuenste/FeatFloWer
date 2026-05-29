module prov_dump_config

implicit none

logical :: use_prov_dump_io = .false.

contains

subroutine set_use_prov_dump(flag)
  logical, intent(in) :: flag
  use_prov_dump_io = flag
end subroutine set_use_prov_dump

end module prov_dump_config
