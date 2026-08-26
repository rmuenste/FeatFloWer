module prov_dump_config

use checkpoint_config, only: checkpoint_set_legacy_format, &
  selected_checkpoint_format, CHECKPOINT_FORMAT_PROVENANCE

implicit none

logical :: use_prov_dump_io = .false.

contains

subroutine set_use_prov_dump(flag)
  logical, intent(in) :: flag
  use_prov_dump_io = flag
  call checkpoint_set_legacy_format(flag)
end subroutine set_use_prov_dump

subroutine sync_use_prov_dump()
  use_prov_dump_io = selected_checkpoint_format == CHECKPOINT_FORMAT_PROVENANCE
end subroutine sync_use_prov_dump

end module prov_dump_config
