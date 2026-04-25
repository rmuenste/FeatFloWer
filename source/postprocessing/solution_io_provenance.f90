module solution_io_provenance

use var_QuadScalar, only: knvt, knet, knat, knel

implicit none

integer :: ifile_prov = 0

type tProvQ2Meta
  integer :: nelem_global = 0
  integer :: nslots = 0
  integer :: owner_count = 0
  integer, allocatable :: global_ids(:,:)
  integer, allocatable :: owner_elem(:)
  integer, allocatable :: owner_slot(:)
  integer, allocatable :: owner_index_by_gid(:)
  integer, allocatable :: duplicate_count(:)
end type tProvQ2Meta

contains

subroutine write_sol_to_file_prov(imax_out, time_ns, output_idx)
  use def_FEAT, only: NLMAX
  use var_QuadScalar, only: QuadSc, LinSc, Temperature, MaterialDistribution
  use var_QuadScalar, only: bTracer, Tracer, myDump, istep_ns, fieldPtr, mg_mesh
  use var_QuadScalar, only: GenLinScalar, GlobalNumberingQ2
  use pp3d_mpi, only: myid, coarse

  implicit none

  integer, intent(in) :: imax_out
  real*8, intent(in) :: time_ns
  integer, optional :: output_idx

  integer :: iout, ndof, nelem, iFld
  character(len=256) :: out_dir
  character(len=60) :: field_name
  type(fieldPtr), dimension(3) :: packed
  type(tProvQ2Meta) :: meta

  if (.not.present(output_idx)) then
    ifile_prov = ifile_prov + 1
    iout = mod(ifile_prov + imax_out - 1, imax_out) + 1
  else
    iout = output_idx
  end if

  nelem = knel(NLMAX)
  ndof = knvt(NLMAX) + knat(NLMAX) + knet(NLMAX) + knel(NLMAX)

  call build_prov_output_dir(iout, out_dir)
  call gather_q2_metadata(meta)
  call write_manifest(out_dir, meta)
  call write_q2_ownership_files(out_dir, meta)

  packed(1)%p => QuadSc%ValU
  packed(2)%p => QuadSc%ValV
  packed(3)%p => QuadSc%ValW
  call write_q2_field(out_dir, 'velocity', 3, packed, meta)

  QuadSc%auxU = mg_mesh%level(nlmax+1)%dcorvg(1,:)
  QuadSc%auxV = mg_mesh%level(nlmax+1)%dcorvg(2,:)
  QuadSc%auxW = mg_mesh%level(nlmax+1)%dcorvg(3,:)
  packed(1)%p => QuadSc%auxU
  packed(2)%p => QuadSc%auxV
  packed(3)%p => QuadSc%auxW
  call write_q2_field(out_dir, 'coordinates', 3, packed, meta)

  if (bTracer .and. allocated(Temperature)) then
    if (myid.ne.0) QuadSc%auxU(1:ndof) = Tracer%Val(NLMAX+1)%x(1:ndof)
    packed(1)%p => QuadSc%auxU
    call write_q2_field(out_dir, 'temperature', 1, packed(1:1), meta)
  end if

  if (allocated(MaterialDistribution)) then
    QuadSc%auxU = 0d0
    QuadSc%auxU((knvt(NLMAX)+knat(NLMAX)+knet(NLMAX))+1:) = MaterialDistribution(NLMAX)%x(1:knel(NLMAX))
    packed(1)%p => QuadSc%auxU
    call write_q2_field(out_dir, 'MaterialDistribution', 1, packed(1:1), meta)
  end if

  if (allocated(GenLinScalar%Fld)) then
    if (myid.ne.0) then
      do iFld = 1, GenLinScalar%nOfFields
        field_name = adjustl(trim(GenLinScalar%prm%cField(iFld)))
        QuadSc%auxU = GenLinScalar%Fld(iFld)%Val
        packed(1)%p => QuadSc%auxU
        call write_q2_field(out_dir, field_name, 1, packed(1:1), meta)
      end do
    end if
  end if

  call write_pressure_field(out_dir, nelem, LinSc%ValP(NLMAX)%x)
  call write_time_field(out_dir, istep_ns, time_ns)

  call free_q2_meta(meta)
end subroutine write_sol_to_file_prov

subroutine read_sol_from_file_prov(startFrom, iLevel, time_ns)
  character(len=*), intent(in) :: startFrom
  integer, intent(in) :: iLevel
  real*8, intent(inout) :: time_ns

  call read_sol_from_file_prov_common(startFrom, iLevel, time_ns)
end subroutine read_sol_from_file_prov

subroutine read_sol_from_file_repart_prov(startFrom, iLevel, time_ns)
  character(len=*), intent(in) :: startFrom
  integer, intent(in) :: iLevel
  real*8, intent(inout) :: time_ns

  call read_sol_from_file_prov_common(startFrom, iLevel, time_ns)
end subroutine read_sol_from_file_repart_prov

subroutine read_sol_from_file_prov_common(startFrom, iLevel, time_ns)
  use def_FEAT, only: NLMIN, NLMAX
  use var_QuadScalar, only: QuadSc, LinSc, Temperature, MaterialDistribution
  use var_QuadScalar, only: bTracer, Tracer, myDump, istep_ns, fieldPtr
  use var_QuadScalar, only: GenLinScalar
  use pp3d_mpi, only: myid, coarse

  implicit none

  character(len=*), intent(in) :: startFrom
  integer, intent(in) :: iLevel
  real*8, intent(inout) :: time_ns

  integer :: ndof, nelem, nslots_q2, nslots_p1, nelem_global, owner_count
  integer, allocatable :: owner_index_map(:,:)
  real*8, allocatable :: pressure_map(:,:,:)
  integer :: iFld
  character(len=256) :: base_dir
  character(len=60) :: field_name
  type(fieldPtr), dimension(3) :: packed

  nelem = knel(NLMAX)
  ndof = knvt(NLMAX) + knat(NLMAX) + knet(NLMAX) + knel(NLMAX)

  base_dir = '_dump_prov/'//trim(adjustl(startFrom))
  call read_manifest(base_dir, nelem_global, nslots_q2, nslots_p1, owner_count)
  call read_q2_ownership_map(base_dir, nelem_global, nslots_q2, owner_index_map)
  call read_pressure_map(base_dir, nelem_global, nslots_p1, pressure_map)

  packed(1)%p => QuadSc%ValU
  packed(2)%p => QuadSc%ValV
  packed(3)%p => QuadSc%ValW
  call read_q2_field(base_dir, 'velocity', 3, packed, owner_index_map, owner_count)

  packed(1)%p => QuadSc%auxU
  packed(2)%p => QuadSc%auxV
  packed(3)%p => QuadSc%auxW
  call read_q2_field(base_dir, 'coordinates', 3, packed, owner_index_map, owner_count)

  call read_pressure_field(base_dir, pressure_map, LinSc%ValP(NLMAX)%x)
  call read_time_field(base_dir, istep_ns, time_ns)

  if (bTracer .and. allocated(Temperature)) then
    packed(1)%p => QuadSc%auxU
    if (q2_field_exists(base_dir, 'temperature')) then
      call read_q2_field(base_dir, 'temperature', 1, packed(1:1), owner_index_map, owner_count)
      if (myid.ne.0) Tracer%Val(NLMAX+1)%x(1:ndof) = QuadSc%auxU(1:ndof)
      Temperature = QuadSc%auxU
    end if
  end if

  if (allocated(MaterialDistribution)) then
    packed(1)%p => QuadSc%auxU
    if (q2_field_exists(base_dir, 'MaterialDistribution')) then
      call read_q2_field(base_dir, 'MaterialDistribution', 1, packed(1:1), owner_index_map, owner_count)
      MaterialDistribution(NLMAX)%x(1:knel(NLMAX)) = nint(QuadSc%auxU((knvt(NLMAX)+knat(NLMAX)+knet(NLMAX))+1:))
    end if
  end if

  if (allocated(GenLinScalar%Fld)) then
    do iFld = 1, GenLinScalar%nOfFields
      field_name = adjustl(trim(GenLinScalar%prm%cField(iFld)))
      if (q2_field_exists(base_dir, field_name)) then
        packed(1)%p => QuadSc%auxU
        call read_q2_field(base_dir, field_name, 1, packed(1:1), owner_index_map, owner_count)
        GenLinScalar%fld(iFld)%Val = QuadSc%auxU
      end if
    end do
  end if

  deallocate(owner_index_map, pressure_map)
end subroutine read_sol_from_file_prov_common

subroutine build_prov_output_dir(iout, out_dir)
  implicit none
  integer, intent(in) :: iout
  character(len=*), intent(out) :: out_dir
  character(len=32) :: idx_str
  integer :: ierr

  write(idx_str,'(I0)') iout
  out_dir = '_dump_prov/'//trim(adjustl(idx_str))
  call execute_command_line('mkdir -p _dump_prov', wait=.true., exitstat=ierr)
  call execute_command_line('mkdir -p '//trim(adjustl(out_dir)), wait=.true., exitstat=ierr)
end subroutine build_prov_output_dir

subroutine gather_q2_metadata(meta)
  use def_FEAT, only: NLMIN
  use var_QuadScalar, only: myDump, GlobalNumberingQ2
  use pp3d_mpi, only: myid, subnodes, coarse
  use pp3d_mpi, only: SENDI_myMPI, RECVI_myMPI, SENDK_myMPI, RECVK_myMPI

  implicit none

  type(tProvQ2Meta), intent(out) :: meta

  integer :: nslots, nlocal, nelem_global
  integer, allocatable :: local_ids(:,:), recv_ids(:,:)
  integer :: iel, ivt, pID, global_iel, gid, max_gid

  nslots = size(myDump%Vertices,2)
  if (myid.ne.0 .and. .not.allocated(GlobalNumberingQ2)) then
    write(*,*) 'solution_io_provenance requires GlobalNumberingQ2 to be allocated.'
    stop 1
  end if

  if (myid.ne.0) then
    nlocal = knel(NLMIN)
    allocate(local_ids(nslots,nlocal))
    do iel = 1, nlocal
      do ivt = 1, nslots
        local_ids(ivt,iel) = GlobalNumberingQ2(myDump%Vertices(iel,ivt))
      end do
    end do
    call SENDI_myMPI(nslots, 0)
    call SENDI_myMPI(nlocal, 0)
    call SENDK_myMPI(local_ids, nslots*nlocal, 0)
    deallocate(local_ids)
  else
    nelem_global = knel(NLMIN)
    meta%nelem_global = nelem_global
    meta%nslots = nslots
    allocate(meta%global_ids(nslots, nelem_global))
    meta%global_ids = 0
    do pID = 1, subnodes
      call RECVI_myMPI(nslots, pID)
      call RECVI_myMPI(nlocal, pID)
      allocate(recv_ids(nslots,nlocal))
      call RECVK_myMPI(recv_ids, nslots*nlocal, pID)
      do iel = 1, nlocal
        global_iel = coarse%pELEMLINK(pID, iel)
        meta%global_ids(:,global_iel) = recv_ids(:,iel)
      end do
      deallocate(recv_ids)
    end do

    max_gid = maxval(meta%global_ids)
    allocate(meta%owner_elem(max_gid), meta%owner_slot(max_gid), meta%owner_index_by_gid(max_gid), meta%duplicate_count(max_gid))
    meta%owner_elem = huge(1)
    meta%owner_slot = huge(1)
    meta%owner_index_by_gid = 0
    meta%duplicate_count = 0

    do global_iel = 1, nelem_global
      do ivt = 1, nslots
        gid = meta%global_ids(ivt,global_iel)
        if (gid.le.0) cycle
        meta%duplicate_count(gid) = meta%duplicate_count(gid) + 1
        if (global_iel.lt.meta%owner_elem(gid)) then
          meta%owner_elem(gid) = global_iel
          meta%owner_slot(gid) = ivt
        else if (global_iel.eq.meta%owner_elem(gid) .and. ivt.lt.meta%owner_slot(gid)) then
          meta%owner_slot(gid) = ivt
        end if
      end do
    end do

    meta%owner_count = 0
    do global_iel = 1, nelem_global
      do ivt = 1, nslots
        gid = meta%global_ids(ivt,global_iel)
        if (gid.le.0) cycle
        if (meta%owner_elem(gid).eq.global_iel .and. meta%owner_slot(gid).eq.ivt) then
          meta%owner_count = meta%owner_count + 1
          meta%owner_index_by_gid(gid) = meta%owner_count
        end if
      end do
    end do
  end if
end subroutine gather_q2_metadata

subroutine free_q2_meta(meta)
  type(tProvQ2Meta), intent(inout) :: meta
  if (allocated(meta%global_ids)) deallocate(meta%global_ids)
  if (allocated(meta%owner_elem)) deallocate(meta%owner_elem)
  if (allocated(meta%owner_slot)) deallocate(meta%owner_slot)
  if (allocated(meta%owner_index_by_gid)) deallocate(meta%owner_index_by_gid)
  if (allocated(meta%duplicate_count)) deallocate(meta%duplicate_count)
  meta%nelem_global = 0
  meta%nslots = 0
  meta%owner_count = 0
end subroutine free_q2_meta

subroutine write_manifest(out_dir, meta)
  use def_FEAT, only: NLMAX
  use pp3d_mpi, only: myid
  implicit none
  character(len=*), intent(in) :: out_dir
  type(tProvQ2Meta), intent(in) :: meta
  integer :: unit_id

  if (myid.ne.0) return
  unit_id = 701
  open(unit=unit_id, file=trim(adjustl(out_dir))//'/manifest.txt', status='replace')
  write(unit_id,'(A)') 'format=ff_prov_dump_v1'
  write(unit_id,'(A,I0)') 'coarse_elements=', meta%nelem_global
  write(unit_id,'(A,I0)') 'q2_slots_per_coarse=', meta%nslots
  write(unit_id,'(A,I0)') 'p1_slots_per_coarse=', 8**(NLMAX-1)
  write(unit_id,'(A,I0)') 'q2_owner_count=', meta%owner_count
  close(unit_id)
end subroutine write_manifest

subroutine write_q2_ownership_files(out_dir, meta)
  use pp3d_mpi, only: myid
  implicit none
  character(len=*), intent(in) :: out_dir
  type(tProvQ2Meta), intent(in) :: meta
  integer :: unit_meta, unit_audit
  integer :: gid, iel, ivt, owner_idx

  if (myid.ne.0) return
  unit_meta = 702
  unit_audit = 703

  open(unit=unit_meta, file=trim(adjustl(out_dir))//'/q2_ownership.csv', status='replace')
  write(unit_meta,'(A)') 'coarse_elem_global,slot,global_q2_id,owner_coarse_elem_global,owner_slot,owner_index,is_owner,duplicate_count'
  do iel = 1, meta%nelem_global
    do ivt = 1, meta%nslots
      gid = meta%global_ids(ivt,iel)
      owner_idx = meta%owner_index_by_gid(gid)
      write(unit_meta,'(8(I0,:,","))') iel, ivt, gid, meta%owner_elem(gid), meta%owner_slot(gid), owner_idx, &
        merge(1,0,meta%owner_elem(gid).eq.iel .and. meta%owner_slot(gid).eq.ivt), meta%duplicate_count(gid)
    end do
  end do
  close(unit_meta)

  open(unit=unit_audit, file=trim(adjustl(out_dir))//'/q2_ownership_audit.csv', status='replace')
  write(unit_audit,'(A)') 'owner_index,global_q2_id,owner_coarse_elem_global,owner_slot,duplicate_count'
  do gid = 1, size(meta%owner_index_by_gid)
    if (meta%owner_index_by_gid(gid).gt.0) then
      write(unit_audit,'(5(I0,:,","))') meta%owner_index_by_gid(gid), gid, meta%owner_elem(gid), meta%owner_slot(gid), meta%duplicate_count(gid)
    end if
  end do
  close(unit_audit)
end subroutine write_q2_ownership_files

subroutine write_q2_field(out_dir, field_name, ncomp, field_pack, meta)
  use def_FEAT, only: NLMIN
  use pp3d_mpi, only: myid, subnodes, coarse
  use pp3d_mpi, only: SENDI_myMPI, RECVI_myMPI, SENDD_myMPI, RECVD_myMPI
  use var_QuadScalar, only: myDump, fieldPtr

  implicit none

  character(len=*), intent(in) :: out_dir
  character(len=*), intent(in) :: field_name
  integer, intent(in) :: ncomp
  type(fieldPtr), dimension(:), intent(in) :: field_pack
  type(tProvQ2Meta), intent(in) :: meta

  integer :: nslots, nlocal, pID, iel, ivt, comp, global_iel, owner_idx, gid
  real*8, allocatable :: local_vals(:,:,:), recv_vals(:,:,:), owner_vals(:,:)
  character(len=256) :: file_name
  integer :: unit_id

  nslots = size(myDump%Vertices,2)
  if (myid.ne.0) then
    nlocal = knel(NLMIN)
    allocate(local_vals(nslots,nlocal,ncomp))
    do comp = 1, ncomp
      do iel = 1, nlocal
        do ivt = 1, nslots
          local_vals(ivt,iel,comp) = field_pack(comp)%p(myDump%Vertices(iel,ivt))
        end do
      end do
    end do
    call SENDI_myMPI(nslots, 0)
    call SENDI_myMPI(nlocal, 0)
    call SENDI_myMPI(ncomp, 0)
    call SENDD_myMPI(local_vals, nslots*nlocal*ncomp, 0)
    deallocate(local_vals)
  else
    allocate(owner_vals(meta%owner_count,ncomp))
    owner_vals = 0d0
    do pID = 1, subnodes
      call RECVI_myMPI(nslots, pID)
      call RECVI_myMPI(nlocal, pID)
      call RECVI_myMPI(comp, pID)
      allocate(recv_vals(nslots,nlocal,comp))
      call RECVD_myMPI(recv_vals, nslots*nlocal*comp, pID)
      do iel = 1, nlocal
        global_iel = coarse%pELEMLINK(pID, iel)
        do ivt = 1, nslots
          gid = meta%global_ids(ivt,global_iel)
          if (meta%owner_elem(gid).eq.global_iel .and. meta%owner_slot(gid).eq.ivt) then
            owner_idx = meta%owner_index_by_gid(gid)
            owner_vals(owner_idx,1:comp) = recv_vals(ivt,iel,1:comp)
          end if
        end do
      end do
      deallocate(recv_vals)
    end do

    file_name = trim(adjustl(out_dir))//'/'//trim(adjustl(field_name))//'.csv'
    unit_id = 710
    open(unit=unit_id, file=file_name, status='replace')
    write(unit_id,'(A)') 'owner_index'//repeat(',value', ncomp)
    do owner_idx = 1, meta%owner_count
      select case(ncomp)
      case(1)
        write(unit_id,'(I0,",",ES24.16E3)') owner_idx, owner_vals(owner_idx,1)
      case(2)
        write(unit_id,'(I0,",",ES24.16E3,",",ES24.16E3)') owner_idx, owner_vals(owner_idx,1), owner_vals(owner_idx,2)
      case default
        write(unit_id,'(I0,",",ES24.16E3,",",ES24.16E3,",",ES24.16E3)') owner_idx, owner_vals(owner_idx,1), owner_vals(owner_idx,2), owner_vals(owner_idx,3)
      end select
    end do
    close(unit_id)
    deallocate(owner_vals)
  end if
end subroutine write_q2_field

subroutine write_pressure_field(out_dir, nelem, pres)
  use def_FEAT, only: NLMIN
  use pp3d_mpi, only: myid, subnodes, coarse
  use pp3d_mpi, only: SENDI_myMPI, RECVI_myMPI, SENDD_myMPI, RECVD_myMPI
  use var_QuadScalar, only: myDump

  implicit none

  character(len=*), intent(in) :: out_dir
  integer, intent(in) :: nelem
  real*8, dimension(:), intent(in) :: pres

  integer :: nslots, nlocal, pID, iel, islot, global_iel, idx
  real*8, allocatable :: local_vals(:,:,:), recv_vals(:,:,:), all_vals(:,:,:)
  integer :: unit_id

  nslots = size(myDump%Elements,2)
  if (myid.ne.0) then
    nlocal = knel(NLMIN)
    allocate(local_vals(nslots,nlocal,4))
    do iel = 1, nlocal
      do islot = 1, nslots
        idx = myDump%Elements(iel,islot)
        local_vals(islot,iel,1) = pres(4*(idx-1)+1)
        local_vals(islot,iel,2) = pres(4*(idx-1)+2)
        local_vals(islot,iel,3) = pres(4*(idx-1)+3)
        local_vals(islot,iel,4) = pres(4*(idx-1)+4)
      end do
    end do
    call SENDI_myMPI(nslots, 0)
    call SENDI_myMPI(nlocal, 0)
    call SENDD_myMPI(local_vals, nslots*nlocal*4, 0)
    deallocate(local_vals)
  else
    allocate(all_vals(nslots,knel(NLMIN),4))
    all_vals = 0d0
    do pID = 1, subnodes
      call RECVI_myMPI(nslots, pID)
      call RECVI_myMPI(nlocal, pID)
      allocate(recv_vals(nslots,nlocal,4))
      call RECVD_myMPI(recv_vals, nslots*nlocal*4, pID)
      do iel = 1, nlocal
        global_iel = coarse%pELEMLINK(pID, iel)
        all_vals(:,global_iel,:) = recv_vals(:,iel,:)
      end do
      deallocate(recv_vals)
    end do

    unit_id = 711
    open(unit=unit_id, file=trim(adjustl(out_dir))//'/pressure.csv', status='replace')
    write(unit_id,'(A)') 'coarse_elem_global,slot,val,val_dx,val_dy,val_dz'
    do global_iel = 1, knel(NLMIN)
      do islot = 1, nslots
        write(unit_id,'(I0,",",I0,4(",",ES24.16E3))') global_iel, islot, all_vals(islot,global_iel,1), &
          all_vals(islot,global_iel,2), all_vals(islot,global_iel,3), all_vals(islot,global_iel,4)
      end do
    end do
    close(unit_id)
    deallocate(all_vals)
  end if
end subroutine write_pressure_field

subroutine write_time_field(out_dir, istep, sim_time)
  use pp3d_mpi, only: myid
  implicit none
  character(len=*), intent(in) :: out_dir
  integer, intent(in) :: istep
  real*8, intent(in) :: sim_time
  integer :: unit_id

  if (myid.ne.0 .and. myid.ne.1) return
  if (myid.eq.1) then
    unit_id = 712
    open(unit=unit_id, file=trim(adjustl(out_dir))//'/time.txt', status='replace')
    write(unit_id,'(ES24.16E3)') sim_time
    write(unit_id,'(I0)') istep
    close(unit_id)
  end if
end subroutine write_time_field

subroutine read_manifest(base_dir, nelem_global, nslots_q2, nslots_p1, owner_count)
  use pp3d_mpi, only: myid
  implicit none
  character(len=*), intent(in) :: base_dir
  integer, intent(out) :: nelem_global, nslots_q2, nslots_p1, owner_count
  integer :: unit_id, ios, pos
  character(len=256) :: line

  nelem_global = 0
  nslots_q2 = 0
  nslots_p1 = 0
  owner_count = 0

  unit_id = 720 + myid
  open(unit=unit_id, file=trim(adjustl(base_dir))//'/manifest.txt', status='old', action='read', iostat=ios)
  if (ios.ne.0) then
    write(*,*) 'Could not open provenance manifest: ', trim(adjustl(base_dir))
    stop 1
  end if

  do
    read(unit_id,'(A)',iostat=ios) line
    if (ios.ne.0) exit
    pos = index(line,'=')
    if (pos.le.0) cycle
    select case (trim(adjustl(line(1:pos-1))))
    case ('coarse_elements')
      read(line(pos+1:),*) nelem_global
    case ('q2_slots_per_coarse')
      read(line(pos+1:),*) nslots_q2
    case ('p1_slots_per_coarse')
      read(line(pos+1:),*) nslots_p1
    case ('q2_owner_count')
      read(line(pos+1:),*) owner_count
    end select
  end do
  close(unit_id)
end subroutine read_manifest

subroutine read_q2_ownership_map(base_dir, nelem_global, nslots_q2, owner_index_map)
  use pp3d_mpi, only: myid
  implicit none
  character(len=*), intent(in) :: base_dir
  integer, intent(in) :: nelem_global, nslots_q2
  integer, allocatable, intent(out) :: owner_index_map(:,:)
  integer :: unit_id, ios
  character(len=512) :: line
  integer :: coarse_elem_global, slot, global_q2_id, owner_coarse, owner_slot, owner_index, is_owner, duplicate_count

  allocate(owner_index_map(nslots_q2, nelem_global))
  owner_index_map = 0
  unit_id = 730 + myid
  open(unit=unit_id, file=trim(adjustl(base_dir))//'/q2_ownership.csv', status='old', action='read', iostat=ios)
  if (ios.ne.0) then
    write(*,*) 'Could not open q2 ownership file in ', trim(adjustl(base_dir))
    stop 1
  end if
  read(unit_id,'(A)') line
  do
    read(unit_id,'(A)',iostat=ios) line
    if (ios.ne.0) exit
    read(line,*) coarse_elem_global, slot, global_q2_id, owner_coarse, owner_slot, owner_index, is_owner, duplicate_count
    owner_index_map(slot, coarse_elem_global) = owner_index
  end do
  close(unit_id)
end subroutine read_q2_ownership_map

subroutine read_pressure_map(base_dir, nelem_global, nslots_p1, pressure_map)
  use pp3d_mpi, only: myid
  implicit none
  character(len=*), intent(in) :: base_dir
  integer, intent(in) :: nelem_global, nslots_p1
  real*8, allocatable, intent(out) :: pressure_map(:,:,:)
  integer :: unit_id, ios
  character(len=512) :: line
  integer :: coarse_elem_global, slot
  real*8 :: v0, v1, v2, v3

  allocate(pressure_map(4,nslots_p1,nelem_global))
  pressure_map = 0
  unit_id = 740 + myid
  open(unit=unit_id, file=trim(adjustl(base_dir))//'/pressure.csv', status='old', action='read', iostat=ios)
  if (ios.ne.0) then
    write(*,*) 'Could not open pressure provenance file in ', trim(adjustl(base_dir))
    stop 1
  end if
  read(unit_id,'(A)') line
  do
    read(unit_id,'(A)',iostat=ios) line
    if (ios.ne.0) exit
    read(line,*) coarse_elem_global, slot, v0, v1, v2, v3
    pressure_map(1,slot,coarse_elem_global) = v0
    pressure_map(2,slot,coarse_elem_global) = v1
    pressure_map(3,slot,coarse_elem_global) = v2
    pressure_map(4,slot,coarse_elem_global) = v3
  end do
  close(unit_id)
end subroutine read_pressure_map

subroutine read_q2_field(base_dir, field_name, ncomp, field_pack, owner_index_map, owner_count)
  use def_FEAT, only: NLMIN
  use pp3d_mpi, only: myid, coarse
  use var_QuadScalar, only: myDump, fieldPtr

  implicit none

  character(len=*), intent(in) :: base_dir
  character(len=*), intent(in) :: field_name
  integer, intent(in) :: ncomp
  type(fieldPtr), dimension(:), intent(inout) :: field_pack
  integer, intent(in) :: owner_index_map(:,:)
  integer, intent(in) :: owner_count

  integer :: unit_id, ios, idx, iel, ivt, owner_idx
  character(len=512) :: line
  real*8, allocatable :: owner_vals(:,:)

  if (myid.eq.0) return

  allocate(owner_vals(owner_count,ncomp))
  owner_vals = 0d0

  unit_id = 750 + myid
  open(unit=unit_id, file=trim(adjustl(base_dir))//'/'//trim(adjustl(field_name))//'.csv', status='old', action='read', iostat=ios)
  if (ios.ne.0) then
    write(*,*) 'Could not open provenance field: ', trim(adjustl(field_name))
    stop 1
  end if
  read(unit_id,'(A)') line
  do
    read(unit_id,'(A)',iostat=ios) line
    if (ios.ne.0) exit
    select case(ncomp)
    case(1)
      read(line,*) idx, owner_vals(idx,1)
    case(2)
      read(line,*) idx, owner_vals(idx,1), owner_vals(idx,2)
    case default
      read(line,*) idx, owner_vals(idx,1), owner_vals(idx,2), owner_vals(idx,3)
    end select
  end do
  close(unit_id)

  do iel = 1, knel(NLMIN)
    do ivt = 1, size(myDump%Vertices,2)
      owner_idx = owner_index_map(ivt, coarse%myELEMLINK(iel))
      do idx = 1, ncomp
        field_pack(idx)%p(myDump%Vertices(iel,ivt)) = owner_vals(owner_idx,idx)
      end do
    end do
  end do

  deallocate(owner_vals)
end subroutine read_q2_field

subroutine read_pressure_field(base_dir, pressure_map, pres)
  use def_FEAT, only: NLMIN
  use pp3d_mpi, only: myid, coarse
  use var_QuadScalar, only: myDump

  implicit none

  character(len=*), intent(in) :: base_dir
  real*8, intent(in) :: pressure_map(:,:,:)
  real*8, dimension(:), intent(inout) :: pres

  integer :: iel, islot, idx

  if (myid.eq.0) return
  do iel = 1, knel(NLMIN)
    do islot = 1, size(myDump%Elements,2)
      idx = myDump%Elements(iel,islot)
      pres(4*(idx-1)+1) = pressure_map(1,islot,coarse%myELEMLINK(iel))
      pres(4*(idx-1)+2) = pressure_map(2,islot,coarse%myELEMLINK(iel))
      pres(4*(idx-1)+3) = pressure_map(3,islot,coarse%myELEMLINK(iel))
      pres(4*(idx-1)+4) = pressure_map(4,islot,coarse%myELEMLINK(iel))
    end do
  end do
end subroutine read_pressure_field

subroutine read_time_field(base_dir, istep, sim_time)
  use pp3d_mpi, only: myid
  implicit none
  character(len=*), intent(in) :: base_dir
  integer, intent(out) :: istep
  real*8, intent(out) :: sim_time
  integer :: unit_id, ios

  unit_id = 760 + myid
  open(unit=unit_id, file=trim(adjustl(base_dir))//'/time.txt', status='old', action='read', iostat=ios)
  if (ios.ne.0) then
    write(*,*) 'Could not open provenance time file in ', trim(adjustl(base_dir))
    stop 1
  end if
  read(unit_id,*) sim_time
  read(unit_id,*) istep
  close(unit_id)
  istep = istep + 1
end subroutine read_time_field

logical function q2_field_exists(base_dir, field_name)
  implicit none
  character(len=*), intent(in) :: base_dir
  character(len=*), intent(in) :: field_name
  integer :: unit_id, ios

  q2_field_exists = .false.
  unit_id = 799
  open(unit=unit_id, file=trim(adjustl(base_dir))//'/'//trim(adjustl(field_name))//'.csv', status='old', action='read', iostat=ios)
  if (ios.eq.0) then
    q2_field_exists = .true.
    close(unit_id)
  end if
end function q2_field_exists

end module solution_io_provenance
