module el_frozen_driver

  use iso_c_binding, only: c_int
  use PP3D_MPI, only: myid, showid, coarse
  use var_QuadScalar, only: mg_mesh, myDump
  use Transport_Q2P1, only: QuadSc
  use dem_query, only: tParticleData, getAllParticles, getAllRemoteParticles, &
                       numLocalParticles, numRemParticles, numTotalParticles, &
                       setForcesMapped, setRemoteForcesMapped
  use fbmaux, only: fbmaux_PointInHex

  implicit none

  integer, parameter :: EL_KERNEL_NONE = 0
  integer, parameter :: EL_KERNEL_TRACER = 1
  integer, parameter :: EL_KERNEL_STOKES = 2
  integer, parameter :: EL_KERNEL_SCHILLER_NAUMANN = 3
  integer, parameter :: EL_MAX_DEBUG_IDS = 64

  type tElClosureCapabilities
    logical :: supports_force = .false.
    logical :: supports_torque = .false.
    logical :: supports_pressure = .false.
    logical :: supports_slip_velocity = .false.
    logical :: supports_reynolds = .false.
  end type tElClosureCapabilities

  character(len=32) :: el_force_kernel = 'none'
  logical :: el_write_diagnostics = .false.
  logical :: el_sampling_debug_enabled = .false.
  logical :: el_write_sampling_summary_debug = .false.
  logical :: el_write_sampling_dof_debug = .false.
  logical :: el_write_sampling_connectivity_debug = .false.
  logical :: el_write_sampling_topology_debug = .false.
  logical :: el_write_sampling_field_debug = .false.
  logical :: el_write_sampling_dumpmap_debug = .false.
  logical :: el_write_sampling_dof_sources_debug = .false.
  logical :: el_apply_forces = .true.
  logical :: el_enable_buoyancy = .false.
  real*8 :: el_fluid_density = 1.0d0
  real*8 :: el_kinematic_viscosity = 1.0d-3
  real*8 :: el_particle_density = 1.0d0
  real*8 :: el_gravity(3) = 0.0d0
  integer :: el_debug_ids(EL_MAX_DEBUG_IDS) = 0
  integer :: el_num_debug_ids = 0

  interface
    subroutine step_el_frozen_trace_c() bind(C, name="step_el_frozen_trace_")
    end subroutine step_el_frozen_trace_c
    subroutine sync_el_frozen_forces_c() bind(C, name="sync_el_frozen_forces_")
    end subroutine sync_el_frozen_forces_c
  end interface

contains

  subroutine el_run_frozen_particle_pass(log_unit, istep)

    integer, intent(in) :: log_unit, istep

    integer, parameter :: NNBAS = 27, NNDER = 10, NNVE = 8, NNDIM = 3

    integer :: num_local_particles, num_particles, ierr, ip
    integer :: ilev, nel_local, found_total, duplicate_total, missing_total
    integer :: local_elem
    integer :: local_num(1), global_num(1)
    integer, allocatable :: found_local(:), found_global(:)
    integer, allocatable :: elem_local(:), elem_global(:)
    integer, allocatable :: owner_local(:), owner_global(:)
    real*8, allocatable :: vel_local(:, :), vel_global(:, :)
    real*8, allocatable :: force_global(:, :), torque_global(:, :)
    real*8, allocatable :: slip_global(:, :), re_global(:)
    type(tParticleData), allocatable :: particles(:)
    type(tElClosureCapabilities) :: caps
    integer :: m, mt, mkeyb, mterm, merr, mprot, msys, mtrc, irecl8
    real*8 :: xi_debug(3), bbox_debug(6)
    integer :: dof_count_debug, kdfg_debug(27), kdfl_debug(27), direct_q2_debug(27)
    real*8 :: basis_debug(27), nodal_u_debug(27), nodal_v_debug(27), nodal_w_debug(27)
    integer :: kvert_debug(8), kedge_debug(12), karea_debug(6)
    real*8 :: xverts_debug(8), yverts_debug(8), zverts_debug(8)

    common /OUTPUT/ m, mt, mkeyb, mterm, merr, mprot, msys, mtrc, irecl8

    include 'mpif.h'

#ifndef PE_SERIAL_MODE
    call el_run_frozen_particle_pass_parallel(log_unit, istep)
    return
#endif

    local_num(1) = 0
    if (myid .ne. 0) local_num(1) = numTotalParticles()
    call MPI_Allreduce(local_num, global_num, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    num_particles = global_num(1)

    if (num_particles .le. 0) then
      if (myid .eq. showid) then
        write(mterm,'(A)') 'Frozen-field particle pass: no PE particles available.'
        write(log_unit,'(A)') 'Frozen-field particle pass: no PE particles available.'
      end if
      return
    end if

    allocate(found_local(num_particles), found_global(num_particles))
    allocate(elem_local(num_particles), elem_global(num_particles))
    allocate(owner_local(num_particles), owner_global(num_particles))
    allocate(vel_local(3, num_particles), vel_global(3, num_particles))
    allocate(force_global(3, num_particles), torque_global(3, num_particles))
    allocate(slip_global(3, num_particles), re_global(num_particles))

    found_local = 0
    elem_local = 0
    owner_local = 0
    vel_local = 0.0d0
    vel_global = 0.0d0
    force_global = 0.0d0
    torque_global = 0.0d0
    slip_global = 0.0d0
    re_global = 0.0d0

    if (myid .ne. 0) then
      num_local_particles = numTotalParticles()
      allocate(particles(num_particles))
      if (num_local_particles .gt. 0) then
        call getAllParticles(particles(1:num_local_particles))
      end if
    end if

    ilev = mg_mesh%nlmax
    nel_local = mg_mesh%level(ilev)%nel

    if (myid .ne. 0 .and. num_local_particles .gt. 0) then
      call el_setup_q2_sampling(ilev)

      do ip = 1, num_particles
        local_elem = 0
        call el_sample_particle_velocity_local(particles(ip), ilev, nel_local, vel_local(:, ip), local_elem, &
                                               xi_debug, bbox_debug, dof_count_debug, kdfg_debug, kdfl_debug, &
                                               basis_debug, nodal_u_debug, nodal_v_debug, nodal_w_debug, direct_q2_debug, &
                                               kvert_debug, kedge_debug, karea_debug, xverts_debug, yverts_debug, zverts_debug)
        if (local_elem .gt. 0) then
          found_local(ip) = 1
          elem_local(ip) = local_elem
          owner_local(ip) = myid
        end if
        call el_write_sampling_debug_record(istep, particles(ip), myid, 'serial', local_elem, found_local(ip), &
                                            xi_debug, bbox_debug, vel_local(:, ip))
        call el_write_sampling_dof_debug_record(istep, particles(ip), myid, 'serial', local_elem, dof_count_debug, &
                                                kdfg_debug, kdfl_debug, basis_debug, nodal_u_debug, nodal_v_debug, nodal_w_debug)
        call el_write_sampling_connectivity_debug_record(istep, particles(ip), myid, 'serial', local_elem, dof_count_debug, &
                                                         kdfg_debug, direct_q2_debug)
        call el_write_sampling_topology_debug_record(istep, particles(ip), myid, 'serial', local_elem, &
                                                     kvert_debug, kedge_debug, karea_debug, xverts_debug, yverts_debug, zverts_debug)
        call el_write_sampling_field_debug_record(istep, particles(ip), myid, 'serial', local_elem, direct_q2_debug, &
                                                  QuadSc%ValU, QuadSc%ValV, QuadSc%ValW)
        call el_write_sampling_dumpmap_debug_record(istep, particles(ip), myid, 'serial', local_elem)
        call el_write_sampling_dof_sources_debug_record(istep, particles(ip), myid, 'serial', local_elem, direct_q2_debug)
      end do
    end if

    call MPI_Allreduce(found_local, found_global, num_particles, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(elem_local, elem_global, num_particles, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(owner_local, owner_global, num_particles, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(vel_local, vel_global, 3*num_particles, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    found_total = 0
    duplicate_total = 0
    missing_total = 0

    do ip = 1, num_particles
      if (found_global(ip) .gt. 0) then
        vel_global(:, ip) = vel_global(:, ip) / dble(found_global(ip))
        found_total = found_total + 1
      else
        missing_total = missing_total + 1
      end if
      if (found_global(ip) .gt. 1) duplicate_total = duplicate_total + 1
    end do

    call el_get_closure_capabilities(caps)

    if (myid .ne. 0 .and. num_local_particles .gt. 0) then
      do ip = 1, num_particles
        call el_evaluate_closure(particles(ip), vel_global(:, ip), found_global(ip), &
                                 force_global(:, ip), torque_global(:, ip), &
                                 slip_global(:, ip), re_global(ip))
        if (el_apply_forces) then
          particles(ip)%force = force_global(:, ip)
          particles(ip)%torque = torque_global(:, ip)
          call setForcesMapped(particles(ip))
        end if
      end do
    end if

    if (myid .eq. showid) then
      write(mterm,'(A,I0,A,I0,A,I0)') 'Frozen-field particle pass: located ', found_total, &
        ' / ', num_particles, ' particles; duplicates=', duplicate_total
      write(log_unit,'(A,I0,A,I0,A,I0)') 'Frozen-field particle pass: located ', found_total, &
        ' / ', num_particles, ' particles; duplicates=', duplicate_total
      if (missing_total .gt. 0) then
        write(mterm,'(A,I0)') '  particles missing from all local meshes: ', missing_total
        write(log_unit,'(A,I0)') '  particles missing from all local meshes: ', missing_total
      end if
      write(mterm,'(A,A)') '  closure kernel: ', trim(el_force_kernel)
      write(log_unit,'(A,A)') '  closure kernel: ', trim(el_force_kernel)
      call el_print_pass_summary(log_unit, num_particles, found_global, force_global, re_global)
      call el_write_particle_diagnostics(istep, particles, num_particles, found_global, elem_global, &
                                         owner_global, vel_global, slip_global, force_global, re_global, caps)
    end if

    if (allocated(particles)) deallocate(particles)
    deallocate(found_local, found_global, elem_local, elem_global, owner_local, owner_global)
    deallocate(vel_local, vel_global, force_global, torque_global, slip_global, re_global)

  end subroutine el_run_frozen_particle_pass

  subroutine el_run_frozen_particle_pass_parallel(log_unit, istep)

    integer, intent(in) :: log_unit, istep

    integer :: num_local_particles, num_remote_particles, num_particles, ierr, ip
    integer :: ilev, nel_local, local_elem
    integer :: found_owned_local, found_owned_global
    integer :: missing_owned_local, missing_owned_global
    integer :: owned_local, owned_global
    integer :: local_found
    type(tParticleData), allocatable :: local_particles(:), remote_particles(:)
    type(tParticleData), allocatable :: diag_particles(:)
    integer, allocatable :: diag_found(:), diag_elem(:), diag_owner(:)
    real*8, allocatable :: diag_vel(:, :), diag_force(:, :), diag_torque(:, :)
    real*8, allocatable :: diag_slip(:, :), diag_re(:)
    real*8 :: vel_sample(3), force(3), torque(3), slip(3), re_p
    real*8 :: local_sum_re, global_sum_re
    real*8 :: local_min_re, global_min_re, local_max_re, global_max_re
    real*8 :: local_sum_force, global_sum_force
    real*8 :: local_min_force, global_min_force, local_max_force, global_max_force
    real*8 :: force_mag
    type(tElClosureCapabilities) :: caps
    integer :: m, mt, mkeyb, mterm, merr, mprot, msys, mtrc, irecl8
    real*8 :: xi_debug(3), bbox_debug(6)
    integer :: dof_count_debug, kdfg_debug(27), kdfl_debug(27), direct_q2_debug(27)
    real*8 :: basis_debug(27), nodal_u_debug(27), nodal_v_debug(27), nodal_w_debug(27)
    integer :: kvert_debug(8), kedge_debug(12), karea_debug(6)
    real*8 :: xverts_debug(8), yverts_debug(8), zverts_debug(8)

    common /OUTPUT/ m, mt, mkeyb, mterm, merr, mprot, msys, mtrc, irecl8

    include 'mpif.h'

    num_local_particles = 0
    num_remote_particles = 0
    if (myid .ne. 0) then
      num_local_particles = numLocalParticles()
      num_remote_particles = numRemParticles()
    end if
    num_particles = num_local_particles + num_remote_particles

    if (num_local_particles .gt. 0) then
      allocate(local_particles(num_local_particles))
      call getAllParticles(local_particles)
    end if

    if (num_remote_particles .gt. 0) then
      allocate(remote_particles(num_remote_particles))
      call getAllRemoteParticles(remote_particles)
    end if

    ilev = 0
    nel_local = 0
    if (myid .ne. 0) then
      ilev = mg_mesh%nlmax
      nel_local = mg_mesh%level(ilev)%nel
      call el_setup_q2_sampling(ilev)
    end if
    call el_get_closure_capabilities(caps)

    found_owned_local = 0
    missing_owned_local = 0
    owned_local = num_local_particles
    local_sum_re = 0.0d0
    local_min_re = huge(local_min_re)
    local_max_re = -huge(local_max_re)
    local_sum_force = 0.0d0
    local_min_force = huge(local_min_force)
    local_max_force = -huge(local_max_force)

    if (num_local_particles .gt. 0) then
      allocate(diag_particles(num_local_particles))
      allocate(diag_found(num_local_particles), diag_elem(num_local_particles), diag_owner(num_local_particles))
      allocate(diag_vel(3, num_local_particles), diag_force(3, num_local_particles))
      allocate(diag_torque(3, num_local_particles), diag_slip(3, num_local_particles), diag_re(num_local_particles))

      do ip = 1, num_local_particles
        call el_sample_particle_velocity_local(local_particles(ip), ilev, nel_local, vel_sample, local_elem, &
                                               xi_debug, bbox_debug, dof_count_debug, kdfg_debug, kdfl_debug, &
                                                   basis_debug, nodal_u_debug, nodal_v_debug, nodal_w_debug, direct_q2_debug, &
                                                   kvert_debug, kedge_debug, karea_debug, xverts_debug, yverts_debug, zverts_debug)
        local_found = merge(1, 0, local_elem .gt. 0)
        call el_evaluate_closure(local_particles(ip), vel_sample, local_found, force, torque, slip, re_p)
        call el_write_sampling_debug_record(istep, local_particles(ip), myid, 'parallel_owned', local_elem, local_found, &
                                            xi_debug, bbox_debug, vel_sample)
        call el_write_sampling_dof_debug_record(istep, local_particles(ip), myid, 'parallel_owned', local_elem, dof_count_debug, &
                                                kdfg_debug, kdfl_debug, basis_debug, nodal_u_debug, nodal_v_debug, nodal_w_debug)
        call el_write_sampling_connectivity_debug_record(istep, local_particles(ip), myid, 'parallel_owned', local_elem, &
                                                         dof_count_debug, kdfg_debug, direct_q2_debug)
        call el_write_sampling_topology_debug_record(istep, local_particles(ip), myid, 'parallel_owned', local_elem, &
                                                     kvert_debug, kedge_debug, karea_debug, xverts_debug, yverts_debug, zverts_debug)
        call el_write_sampling_field_debug_record(istep, local_particles(ip), myid, 'parallel_owned', local_elem, direct_q2_debug, &
                                                  QuadSc%ValU, QuadSc%ValV, QuadSc%ValW)
        call el_write_sampling_dumpmap_debug_record(istep, local_particles(ip), myid, 'parallel_owned', local_elem)
        call el_write_sampling_dof_sources_debug_record(istep, local_particles(ip), myid, 'parallel_owned', local_elem, direct_q2_debug)

        if (el_apply_forces) then
          local_particles(ip)%force = force
          local_particles(ip)%torque = torque
          call setForcesMapped(local_particles(ip))
        end if

        if (local_found .gt. 0) then
          found_owned_local = found_owned_local + 1
        else
          missing_owned_local = missing_owned_local + 1
        end if

        force_mag = sqrt(sum(force*force))
        local_sum_re = local_sum_re + re_p
        local_min_re = min(local_min_re, re_p)
        local_max_re = max(local_max_re, re_p)
        local_sum_force = local_sum_force + force_mag
        local_min_force = min(local_min_force, force_mag)
        local_max_force = max(local_max_force, force_mag)

        diag_particles(ip) = local_particles(ip)
        diag_found(ip) = local_found
        diag_elem(ip) = local_elem
        diag_owner(ip) = myid
        diag_vel(:, ip) = vel_sample
        diag_force(:, ip) = force
        diag_torque(:, ip) = torque
        diag_slip(:, ip) = slip
        diag_re(ip) = re_p
      end do
    end if

    do ip = 1, num_remote_particles
      call el_sample_particle_velocity_local(remote_particles(ip), ilev, nel_local, vel_sample, local_elem, &
                                             xi_debug, bbox_debug, dof_count_debug, kdfg_debug, kdfl_debug, &
                                         basis_debug, nodal_u_debug, nodal_v_debug, nodal_w_debug, direct_q2_debug, &
                                         kvert_debug, kedge_debug, karea_debug, xverts_debug, yverts_debug, zverts_debug)
      local_found = merge(1, 0, local_elem .gt. 0)
      call el_evaluate_closure(remote_particles(ip), vel_sample, local_found, force, torque, slip, re_p)

      call el_write_sampling_debug_record(istep, remote_particles(ip), myid, 'parallel_remote', local_elem, local_found, &
                                          xi_debug, bbox_debug, vel_sample)
      call el_write_sampling_dof_debug_record(istep, remote_particles(ip), myid, 'parallel_remote', local_elem, dof_count_debug, &
                                              kdfg_debug, kdfl_debug, basis_debug, nodal_u_debug, nodal_v_debug, nodal_w_debug)
      call el_write_sampling_connectivity_debug_record(istep, remote_particles(ip), myid, 'parallel_remote', local_elem, &
                                                       dof_count_debug, kdfg_debug, direct_q2_debug)
      call el_write_sampling_topology_debug_record(istep, remote_particles(ip), myid, 'parallel_remote', local_elem, &
                                                   kvert_debug, kedge_debug, karea_debug, xverts_debug, yverts_debug, zverts_debug)
      call el_write_sampling_field_debug_record(istep, remote_particles(ip), myid, 'parallel_remote', local_elem, direct_q2_debug, &
                                                QuadSc%ValU, QuadSc%ValV, QuadSc%ValW)
      call el_write_sampling_dumpmap_debug_record(istep, remote_particles(ip), myid, 'parallel_remote', local_elem)
      call el_write_sampling_dof_sources_debug_record(istep, remote_particles(ip), myid, 'parallel_remote', local_elem, direct_q2_debug)

      if (el_apply_forces) then
        remote_particles(ip)%force = force
        remote_particles(ip)%torque = torque
        call setRemoteForcesMapped(remote_particles(ip))
      end if
    end do

    if (num_local_particles .eq. 0) then
      local_min_re = huge(local_min_re)
      local_max_re = -huge(local_max_re)
      local_min_force = huge(local_min_force)
      local_max_force = -huge(local_max_force)
    end if

    if (el_apply_forces .and. myid .ne. 0) call sync_el_frozen_forces_c()

    call MPI_Allreduce(owned_local, owned_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(found_owned_local, found_owned_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(missing_owned_local, missing_owned_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(local_sum_re, global_sum_re, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(local_min_re, global_min_re, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(local_max_re, global_max_re, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(local_sum_force, global_sum_force, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(local_min_force, global_min_force, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(local_max_force, global_max_force, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    if (owned_global .le. 0) then
      if (myid .eq. showid) then
        write(mterm,'(A)') 'Frozen-field particle pass: no locally owned PE particles available.'
        write(log_unit,'(A)') 'Frozen-field particle pass: no locally owned PE particles available.'
      end if
    else if (myid .eq. showid) then
      write(mterm,'(A,I0,A,I0,A)') 'Frozen-field particle pass: located ', found_owned_global, &
        ' / ', owned_global, ' locally owned particles'
      write(log_unit,'(A,I0,A,I0,A)') 'Frozen-field particle pass: located ', found_owned_global, &
        ' / ', owned_global, ' locally owned particles'
      if (missing_owned_global .gt. 0) then
        write(mterm,'(A,I0)') '  locally owned particles missing from all local meshes: ', missing_owned_global
        write(log_unit,'(A,I0)') '  locally owned particles missing from all local meshes: ', missing_owned_global
      end if
      write(mterm,'(A,A)') '  closure kernel: ', trim(el_force_kernel)
      write(log_unit,'(A,A)') '  closure kernel: ', trim(el_force_kernel)
      write(mterm,'(A,I0)') '  particle count = ', owned_global
      write(log_unit,'(A,I0)') '  particle count = ', owned_global
      write(mterm,'(A,3(ES14.6,1X))') '  Re_p    avg/min/max = ', &
        global_sum_re / dble(owned_global), global_min_re, global_max_re
      write(log_unit,'(A,3(ES14.6,1X))') '  Re_p    avg/min/max = ', &
        global_sum_re / dble(owned_global), global_min_re, global_max_re
      write(mterm,'(A,3(ES14.6,1X))') '  |F|     avg/min/max = ', &
        global_sum_force / dble(owned_global), global_min_force, global_max_force
      write(log_unit,'(A,3(ES14.6,1X))') '  |F|     avg/min/max = ', &
        global_sum_force / dble(owned_global), global_min_force, global_max_force
      write(mterm,'(A,ES14.6)') '  found_count avg = ', dble(found_owned_global) / dble(owned_global)
      write(log_unit,'(A,ES14.6)') '  found_count avg = ', dble(found_owned_global) / dble(owned_global)
    end if

    if (num_local_particles .gt. 0) then
      call el_write_particle_diagnostics_parallel(istep, diag_particles, num_local_particles, diag_found, &
                                                  diag_elem, diag_owner, diag_vel, diag_slip, &
                                                  diag_force, diag_re, caps)
    end if

    if (allocated(local_particles)) deallocate(local_particles)
    if (allocated(remote_particles)) deallocate(remote_particles)
    if (allocated(diag_particles)) deallocate(diag_particles)
    if (allocated(diag_found)) deallocate(diag_found, diag_elem, diag_owner)
    if (allocated(diag_vel)) deallocate(diag_vel, diag_force, diag_torque, diag_slip, diag_re)

  end subroutine el_run_frozen_particle_pass_parallel

  subroutine el_step_frozen_particles(log_unit, istep)

    integer, intent(in) :: log_unit, istep

    integer :: ierr
    integer :: m, mt, mkeyb, mterm, merr, mprot, msys, mtrc, irecl8
    integer :: local_before_max(1), local_before_min(1)
    integer :: global_before_max(1), global_before_min(1)
    integer :: local_after_max(1), local_after_min(1)
    integer :: global_after_max(1), global_after_min(1)
    integer :: removed_particles

    common /OUTPUT/ m, mt, mkeyb, mterm, merr, mprot, msys, mtrc, irecl8

    include 'mpif.h'

#ifndef PE_SERIAL_MODE
    call el_step_frozen_particles_parallel(log_unit, istep)
    return
#endif

    local_before_max(1) = -huge(local_before_max(1))
    local_before_min(1) =  huge(local_before_min(1))
    if (myid .ne. 0) then
      local_before_max(1) = numTotalParticles()
      local_before_min(1) = local_before_max(1)
    end if
    call MPI_Allreduce(local_before_max, global_before_max, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(local_before_min, global_before_min, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)

    if (global_before_max(1) .ne. global_before_min(1)) then
      if (myid .eq. showid) then
        write(mterm,'(A,I0,A,I0)') 'Frozen-field PE step inconsistency before step: min/max particle count = ', &
          global_before_min(1), ' / ', global_before_max(1)
        write(log_unit,'(A,I0,A,I0)') 'Frozen-field PE step inconsistency before step: min/max particle count = ', &
          global_before_min(1), ' / ', global_before_max(1)
      end if
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    if (global_before_max(1) .le. 0) then
      if (myid .eq. showid) then
        write(mterm,'(A,I0,A)') 'Frozen-field PE step ', istep, ': no particles to advance.'
        write(log_unit,'(A,I0,A)') 'Frozen-field PE step ', istep, ': no particles to advance.'
      end if
      return
    end if

    if (myid .ne. 0) then
      call step_el_frozen_trace_c()
    end if

    local_after_max(1) = -huge(local_after_max(1))
    local_after_min(1) =  huge(local_after_min(1))
    if (myid .ne. 0) then
      local_after_max(1) = numTotalParticles()
      local_after_min(1) = local_after_max(1)
    end if
    call MPI_Allreduce(local_after_max, global_after_max, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(local_after_min, global_after_min, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)

    if (global_after_max(1) .ne. global_after_min(1)) then
      if (myid .eq. showid) then
        write(mterm,'(A,I0,A,I0)') 'Frozen-field PE step inconsistency after step: min/max particle count = ', &
          global_after_min(1), ' / ', global_after_max(1)
        write(log_unit,'(A,I0,A,I0)') 'Frozen-field PE step inconsistency after step: min/max particle count = ', &
          global_after_min(1), ' / ', global_after_max(1)
      end if
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    removed_particles = max(0, global_before_max(1) - global_after_max(1))

    if (myid .eq. showid) then
      write(mterm,'(A,I0,A,I0,A,I0,A,I0)') 'Frozen-field PE step ', istep, &
        ': particles before/after = ', global_before_max(1), ' / ', global_after_max(1), &
        ', outflow removed = ', removed_particles
      write(log_unit,'(A,I0,A,I0,A,I0,A,I0)') 'Frozen-field PE step ', istep, &
        ': particles before/after = ', global_before_max(1), ' / ', global_after_max(1), &
        ', outflow removed = ', removed_particles
    end if

  end subroutine el_step_frozen_particles

  subroutine el_step_frozen_particles_parallel(log_unit, istep)

    integer, intent(in) :: log_unit, istep

    integer :: ierr
    integer :: local_before, global_before
    integer :: local_after, global_after
    integer :: removed_particles
    integer :: m, mt, mkeyb, mterm, merr, mprot, msys, mtrc, irecl8

    common /OUTPUT/ m, mt, mkeyb, mterm, merr, mprot, msys, mtrc, irecl8

    include 'mpif.h'

    local_before = 0
    if (myid .ne. 0) local_before = numLocalParticles()
    call MPI_Allreduce(local_before, global_before, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    if (global_before .le. 0) then
      if (myid .eq. showid) then
        write(mterm,'(A,I0,A)') 'Frozen-field PE step ', istep, ': no particles to advance.'
        write(log_unit,'(A,I0,A)') 'Frozen-field PE step ', istep, ': no particles to advance.'
      end if
      return
    end if

    if (myid .ne. 0) then
      call step_el_frozen_trace_c()
    end if

    local_after = 0
    if (myid .ne. 0) local_after = numLocalParticles()
    call MPI_Allreduce(local_after, global_after, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    removed_particles = max(0, global_before - global_after)

    if (myid .eq. showid) then
      write(mterm,'(A,I0,A,I0,A,I0,A,I0)') 'Frozen-field PE step ', istep, &
        ': owned particles before/after = ', global_before, ' / ', global_after, &
        ', outflow removed = ', removed_particles
      write(log_unit,'(A,I0,A,I0,A,I0,A,I0)') 'Frozen-field PE step ', istep, &
        ': owned particles before/after = ', global_before, ' / ', global_after, &
        ', outflow removed = ', removed_particles
    end if

  end subroutine el_step_frozen_particles_parallel

  subroutine el_setup_q2_sampling(ilev_target)

    integer, intent(in) :: ilev_target

    integer, parameter :: NNBAS = 27, NNDER = 10, NNVE = 8, NNDIM = 3
    integer :: ieltyp, i
    real*8 :: dx(nnve), dy(nnve), dz(nnve), djac(3, 3), detj
    real*8 :: dbas(nndim, nnbas, nnder)
    logical :: bder(nnder)
    integer :: kve(nnve), ndim, iel
    integer :: kdfg(nnbas), kdfl(nnbas), idfl
    integer :: ilev, nlev, nlmin, nlmax
    integer :: ielstp, ielmax, ielmul, ielold, ielnew

    common /ELEM/ dx, dy, dz, djac, detj, dbas, bder, kve, iel, ndim
    common /COAUX1/ kdfg, kdfl, idfl
    common /MGPAR/  ilev, nlev, nlmin, nlmax, ielstp, ielmax, ielmul, ielold, ielnew

    integer, external :: ndfl
    external :: e013, setlev

    ilev = ilev_target
    call setlev(2)
    do i = 1, nnder
      bder(i) = .false.
    end do
    bder(1) = .true.

    ieltyp = -1
    call e013(0d0, 0d0, 0d0, ieltyp)
    idfl = ndfl(ieltyp)

  end subroutine el_setup_q2_sampling

  subroutine el_sample_particle_velocity_local(particle, ilev, nel_local, vel_sample, elem_id, xi_out, bbox_out, &
                                               dof_count_out, kdfg_out, kdfl_out, basis_out, nodal_u_out, nodal_v_out, nodal_w_out, direct_q2_out, &
                                               kvert_out, kedge_out, karea_out, xverts_out, yverts_out, zverts_out)

    type(tParticleData), intent(in) :: particle
    integer, intent(in) :: ilev, nel_local
    real*8, intent(out) :: vel_sample(3)
    integer, intent(out) :: elem_id
    real*8, intent(out), optional :: xi_out(3), bbox_out(6)
    integer, intent(out), optional :: dof_count_out, kdfg_out(27), kdfl_out(27)
    real*8, intent(out), optional :: basis_out(27), nodal_u_out(27), nodal_v_out(27), nodal_w_out(27)
    integer, intent(out), optional :: direct_q2_out(27)
    integer, intent(out), optional :: kvert_out(8), kedge_out(12), karea_out(6)
    real*8, intent(out), optional :: xverts_out(8), yverts_out(8), zverts_out(8)

    integer, parameter :: NNBAS = 27, NNDER = 10, NNVE = 8, NNDIM = 3

    integer :: iel, ive, ig, il
    integer :: kdfg(nnbas), kdfl(nnbas), idfl
    integer :: kve(nnve), ndim, iel_common
    real*8 :: dx(nnve), dy(nnve), dz(nnve), djac(3, 3), detj
    real*8 :: dbas(nndim, nnbas, nnder)
    logical :: bder(nnder)
    real*8 :: xverts(8), yverts(8), zverts(8)
    real*8 :: xi1, xi2, xi3
    real*8 :: xmin, xmax, ymin, ymax, zmin, zmax
    real*8 :: eps_box
    logical :: found

    common /ELEM/ dx, dy, dz, djac, detj, dbas, bder, kve, iel_common, ndim
    common /COAUX1/ kdfg, kdfl, idfl

    external :: e013, ndfgl

    vel_sample = 0.0d0
    elem_id = 0
    eps_box = 1.0d-10
    if (present(xi_out)) xi_out = 0.0d0
    if (present(bbox_out)) bbox_out = 0.0d0
    if (present(dof_count_out)) dof_count_out = 0
    if (present(kdfg_out)) kdfg_out = 0
    if (present(kdfl_out)) kdfl_out = 0
    if (present(direct_q2_out)) direct_q2_out = 0
    if (present(basis_out)) basis_out = 0.0d0
    if (present(nodal_u_out)) nodal_u_out = 0.0d0
    if (present(nodal_v_out)) nodal_v_out = 0.0d0
    if (present(nodal_w_out)) nodal_w_out = 0.0d0
    if (present(kvert_out)) kvert_out = 0
    if (present(kedge_out)) kedge_out = 0
    if (present(karea_out)) karea_out = 0
    if (present(xverts_out)) xverts_out = 0.0d0
    if (present(yverts_out)) yverts_out = 0.0d0
    if (present(zverts_out)) zverts_out = 0.0d0

    do iel = 1, nel_local
      do ive = 1, nnve
        ig = mg_mesh%level(ilev)%kvert(ive, iel)
        xverts(ive) = mg_mesh%level(ilev)%dcorvg(1, ig)
        yverts(ive) = mg_mesh%level(ilev)%dcorvg(2, ig)
        zverts(ive) = mg_mesh%level(ilev)%dcorvg(3, ig)
      end do

      xmin = minval(xverts)
      xmax = maxval(xverts)
      ymin = minval(yverts)
      ymax = maxval(yverts)
      zmin = minval(zverts)
      zmax = maxval(zverts)

      if (particle%position(1) .lt. xmin - eps_box) cycle
      if (particle%position(1) .gt. xmax + eps_box) cycle
      if (particle%position(2) .lt. ymin - eps_box) cycle
      if (particle%position(2) .gt. ymax + eps_box) cycle
      if (particle%position(3) .lt. zmin - eps_box) cycle
      if (particle%position(3) .gt. zmax + eps_box) cycle

      xi1 = 0.0d0
      xi2 = 0.0d0
      xi3 = 0.0d0
      found = fbmaux_PointInHex(particle%position(1), particle%position(2), particle%position(3), &
                                xverts, yverts, zverts, xi1, xi2, xi3, iel)
      if (.not. found) cycle

      call ndfgl(iel, 1, 13, mg_mesh%level(ilev)%kvert, mg_mesh%level(ilev)%kedge, &
                 mg_mesh%level(ilev)%karea, kdfg, kdfl)

      if (present(direct_q2_out)) then
        direct_q2_out(1:8) = mg_mesh%level(ilev)%kvert(:, iel)
        direct_q2_out(9:20) = mg_mesh%level(ilev)%nvt + mg_mesh%level(ilev)%kedge(:, iel)
        direct_q2_out(21:26) = mg_mesh%level(ilev)%nvt + mg_mesh%level(ilev)%net + mg_mesh%level(ilev)%karea(:, iel)
        direct_q2_out(27) = mg_mesh%level(ilev)%nvt + mg_mesh%level(ilev)%net + mg_mesh%level(ilev)%nat + iel
      end if

      do ive = 1, nnve
        ig = mg_mesh%level(ilev)%kvert(ive, iel)
        kve(ive) = ig
        dx(ive) = xverts(ive)
        dy(ive) = yverts(ive)
        dz(ive) = zverts(ive)
      end do

      call e013(xi1, xi2, xi3, 0)

      do il = 1, idfl
        ig = kdfg(il)
        ive = kdfl(il)
        if (present(kdfg_out)) kdfg_out(il) = ig
        if (present(kdfl_out)) kdfl_out(il) = ive
        if (present(basis_out)) basis_out(il) = dbas(1, ive, 1)
        if (present(nodal_u_out)) nodal_u_out(il) = QuadSc%ValU(ig)
        if (present(nodal_v_out)) nodal_v_out(il) = QuadSc%ValV(ig)
        if (present(nodal_w_out)) nodal_w_out(il) = QuadSc%ValW(ig)
        vel_sample(1) = vel_sample(1) + QuadSc%ValU(ig) * dbas(1, ive, 1)
        vel_sample(2) = vel_sample(2) + QuadSc%ValV(ig) * dbas(1, ive, 1)
        vel_sample(3) = vel_sample(3) + QuadSc%ValW(ig) * dbas(1, ive, 1)
      end do

      elem_id = iel
      if (present(xi_out)) xi_out = (/xi1, xi2, xi3/)
      if (present(bbox_out)) bbox_out = (/xmin, xmax, ymin, ymax, zmin, zmax/)
      if (present(dof_count_out)) dof_count_out = idfl
      if (present(kvert_out)) kvert_out = mg_mesh%level(ilev)%kvert(:, iel)
      if (present(kedge_out)) kedge_out = mg_mesh%level(ilev)%kedge(:, iel)
      if (present(karea_out)) karea_out = mg_mesh%level(ilev)%karea(:, iel)
      if (present(xverts_out)) xverts_out = xverts
      if (present(yverts_out)) yverts_out = yverts
      if (present(zverts_out)) zverts_out = zverts
      return
    end do

  end subroutine el_sample_particle_velocity_local

  logical function el_any_sampling_debug_enabled()

    el_any_sampling_debug_enabled = el_sampling_debug_enabled .or. &
      el_write_sampling_summary_debug .or. el_write_sampling_dof_debug .or. &
      el_write_sampling_connectivity_debug .or. el_write_sampling_topology_debug .or. &
      el_write_sampling_field_debug .or. el_write_sampling_dumpmap_debug .or. &
      el_write_sampling_dof_sources_debug

  end function el_any_sampling_debug_enabled

  logical function el_is_debug_particle(particle_id)

    integer, intent(in) :: particle_id
    integer :: i

    el_is_debug_particle = .false.
    if (.not. el_any_sampling_debug_enabled()) return

    do i = 1, el_num_debug_ids
      if (el_debug_ids(i) .eq. particle_id) then
        el_is_debug_particle = .true.
        return
      end if
    end do

  end function el_is_debug_particle

  subroutine el_write_sampling_debug_record(istep, particle, rank_id, sample_kind, elem_id, found_count, xi, bbox, carrier_vel)

    integer, intent(in) :: istep, rank_id, elem_id, found_count
    type(tParticleData), intent(in) :: particle
    character(len=*), intent(in) :: sample_kind
    real*8, intent(in) :: xi(3), bbox(6), carrier_vel(3)

    integer :: unit_id, ios
    logical :: file_exists
    character(len=128) :: filename

    if (.not. el_write_sampling_summary_debug) return
    if (.not. el_is_debug_particle(particle%uniqueIdx)) return

    write(filename,'(A,I4.4,A,I6.6,A)') 'el_frozen_sampling_debug_rank', rank_id, '_step', istep, '.csv'
    inquire(file=trim(filename), exist=file_exists)
    open(newunit=unit_id, file=trim(filename), status='unknown', position='append', action='write', iostat=ios)
    if (ios .ne. 0) return

    if (.not. file_exists) then
      write(unit_id,'(A)') 'id,system_id,rank,sample_kind,found_count,elem_id,x,y,z,ux,uy,uz,xi1,xi2,xi3,xmin,xmax,ymin,ymax,zmin,zmax'
    end if

    write(unit_id,'(3(I0,","),A,",",2(I0,","),14(ES18.10,","),ES18.10)') particle%uniqueIdx, &
      particle%systemIdx, rank_id, trim(sample_kind), found_count, elem_id, &
      particle%position(1), particle%position(2), particle%position(3), &
      carrier_vel(1), carrier_vel(2), carrier_vel(3), &
      xi(1), xi(2), xi(3), bbox(1), bbox(2), bbox(3), bbox(4), bbox(5), bbox(6)

    close(unit_id)

  end subroutine el_write_sampling_debug_record

  subroutine el_write_sampling_dof_debug_record(istep, particle, rank_id, sample_kind, elem_id, dof_count, &
                                                kdfg, kdfl, basis, nodal_u, nodal_v, nodal_w)

    integer, intent(in) :: istep, rank_id, elem_id, dof_count
    type(tParticleData), intent(in) :: particle
    character(len=*), intent(in) :: sample_kind
    integer, intent(in) :: kdfg(:), kdfl(:)
    real*8, intent(in) :: basis(:), nodal_u(:), nodal_v(:), nodal_w(:)

    integer :: unit_id, ios, il
    logical :: file_exists
    character(len=128) :: filename

    if (.not. el_write_sampling_dof_debug) return
    if (.not. el_is_debug_particle(particle%uniqueIdx)) return

    write(filename,'(A,I4.4,A,I6.6,A)') 'el_frozen_sampling_dofs_rank', rank_id, '_step', istep, '.csv'
    inquire(file=trim(filename), exist=file_exists)
    open(newunit=unit_id, file=trim(filename), status='unknown', position='append', action='write', iostat=ios)
    if (ios .ne. 0) return

    if (.not. file_exists) then
      write(unit_id,'(A)') 'id,system_id,rank,sample_kind,elem_id,dof_idx,kdfg,kdfl,basis,u,v,w'
    end if

    do il = 1, dof_count
      write(unit_id,'(3(I0,","),A,",",4(I0,","),3(ES18.10,","),ES18.10)') particle%uniqueIdx, particle%systemIdx, &
        rank_id, trim(sample_kind), elem_id, il, kdfg(il), kdfl(il), basis(il), nodal_u(il), nodal_v(il), nodal_w(il)
    end do

    close(unit_id)

  end subroutine el_write_sampling_dof_debug_record

  subroutine el_write_sampling_connectivity_debug_record(istep, particle, rank_id, sample_kind, elem_id, dof_count, &
                                                         ndfgl_kdfg, direct_q2)

    integer, intent(in) :: istep, rank_id, elem_id, dof_count
    type(tParticleData), intent(in) :: particle
    character(len=*), intent(in) :: sample_kind
    integer, intent(in) :: ndfgl_kdfg(:), direct_q2(:)

    integer :: unit_id, ios, il
    logical :: file_exists
    character(len=128) :: filename

    if (.not. el_write_sampling_connectivity_debug) return
    if (.not. el_is_debug_particle(particle%uniqueIdx)) return

    write(filename,'(A,I4.4,A,I6.6,A)') 'el_frozen_sampling_connectivity_rank', rank_id, '_step', istep, '.csv'
    inquire(file=trim(filename), exist=file_exists)
    open(newunit=unit_id, file=trim(filename), status='unknown', position='append', action='write', iostat=ios)
    if (ios .ne. 0) return

    if (.not. file_exists) then
      write(unit_id,'(A)') 'id,system_id,rank,sample_kind,elem_id,slot,ndfgl_kdfg,direct_q2'
    end if

    do il = 1, max(dof_count, size(direct_q2))
      write(unit_id,'(3(I0,","),A,",",4(I0,","),I0)') particle%uniqueIdx, particle%systemIdx, rank_id, &
        trim(sample_kind), elem_id, il, ndfgl_kdfg(il), direct_q2(il)
    end do

    close(unit_id)

  end subroutine el_write_sampling_connectivity_debug_record

  subroutine el_write_sampling_topology_debug_record(istep, particle, rank_id, sample_kind, elem_id, &
                                                     kvert, kedge, karea, xverts, yverts, zverts)

    integer, intent(in) :: istep, rank_id, elem_id
    type(tParticleData), intent(in) :: particle
    character(len=*), intent(in) :: sample_kind
    integer, intent(in) :: kvert(8), kedge(12), karea(6)
    real*8, intent(in) :: xverts(8), yverts(8), zverts(8)

    integer :: unit_id, ios
    logical :: file_exists
    character(len=128) :: filename

    if (.not. el_write_sampling_topology_debug) return
    if (.not. el_is_debug_particle(particle%uniqueIdx)) return

    write(filename,'(A,I4.4,A,I6.6,A)') 'el_frozen_sampling_topology_rank', rank_id, '_step', istep, '.csv'
    inquire(file=trim(filename), exist=file_exists)
    open(newunit=unit_id, file=trim(filename), status='unknown', position='append', action='write', iostat=ios)
    if (ios .ne. 0) return

    if (.not. file_exists) then
      write(unit_id,'(A)') 'id,system_id,rank,sample_kind,elem_id,kvert1,kvert2,kvert3,kvert4,kvert5,kvert6,kvert7,kvert8,' // &
                           'kedge1,kedge2,kedge3,kedge4,kedge5,kedge6,kedge7,kedge8,kedge9,kedge10,kedge11,kedge12,' // &
                           'karea1,karea2,karea3,karea4,karea5,karea6,' // &
                           'x1,x2,x3,x4,x5,x6,x7,x8,y1,y2,y3,y4,y5,y6,y7,y8,z1,z2,z3,z4,z5,z6,z7,z8'
    end if

    write(unit_id,'(3(I0,","),A,27(",",I0),24(",",ES18.10))') particle%uniqueIdx, particle%systemIdx, rank_id, trim(sample_kind), &
      elem_id, kvert(1), kvert(2), kvert(3), kvert(4), kvert(5), kvert(6), kvert(7), kvert(8), &
      kedge(1), kedge(2), kedge(3), kedge(4), kedge(5), kedge(6), kedge(7), kedge(8), kedge(9), kedge(10), kedge(11), kedge(12), &
      karea(1), karea(2), karea(3), karea(4), karea(5), karea(6), &
      xverts(1), xverts(2), xverts(3), xverts(4), xverts(5), xverts(6), xverts(7), xverts(8), &
      yverts(1), yverts(2), yverts(3), yverts(4), yverts(5), yverts(6), yverts(7), yverts(8), &
      zverts(1), zverts(2), zverts(3), zverts(4), zverts(5), zverts(6), zverts(7), zverts(8)

    close(unit_id)

  end subroutine el_write_sampling_topology_debug_record

  subroutine el_write_sampling_field_debug_record(istep, particle, rank_id, sample_kind, elem_id, direct_q2, &
                                                  val_u, val_v, val_w)

    integer, intent(in) :: istep, rank_id, elem_id
    type(tParticleData), intent(in) :: particle
    character(len=*), intent(in) :: sample_kind
    integer, intent(in) :: direct_q2(27)
    real*8, intent(in) :: val_u(:), val_v(:), val_w(:)

    integer :: unit_id, ios, il, jl
    integer :: sorted_idx(27), tmp
    logical :: file_exists
    character(len=128) :: filename

    if (.not. el_write_sampling_field_debug) return
    if (.not. el_is_debug_particle(particle%uniqueIdx)) return
    if (elem_id .le. 0) return

    sorted_idx = direct_q2
    do il = 1, 26
      do jl = il + 1, 27
        if (sorted_idx(jl) .lt. sorted_idx(il)) then
          tmp = sorted_idx(il)
          sorted_idx(il) = sorted_idx(jl)
          sorted_idx(jl) = tmp
        end if
      end do
    end do

    write(filename,'(A,I4.4,A,I6.6,A)') 'el_frozen_sampling_field_rank', rank_id, '_step', istep, '.csv'
    inquire(file=trim(filename), exist=file_exists)
    open(newunit=unit_id, file=trim(filename), status='unknown', position='append', action='write', iostat=ios)
    if (ios .ne. 0) return

    if (.not. file_exists) then
      write(unit_id,'(A)') 'id,system_id,rank,sample_kind,elem_id,slot,global_dof,u,v,w'
    end if

    do il = 1, 27
      write(unit_id,'(3(I0,","),A,",",3(I0,","),3(ES18.10,","),ES18.10)') particle%uniqueIdx, particle%systemIdx, rank_id, &
        trim(sample_kind), elem_id, il, sorted_idx(il), val_u(sorted_idx(il)), val_v(sorted_idx(il)), val_w(sorted_idx(il))
    end do

    close(unit_id)

  end subroutine el_write_sampling_field_debug_record

  subroutine el_find_dump_coarse_mapping(elem_id, coarse_row, coarse_elem_global, subelem_slot)

    integer, intent(in) :: elem_id
    integer, intent(out) :: coarse_row, coarse_elem_global, subelem_slot

    integer :: iel, jel, nrow, nsub, nv

    coarse_row = 0
    coarse_elem_global = 0
    subelem_slot = 0

    if (.not. allocated(myDump%Elements)) return

    nrow = size(myDump%Elements, 1)
    nsub = size(myDump%Elements, 2)

    do iel = 1, nrow
      do jel = 1, nsub
        if (myDump%Elements(iel, jel) .eq. elem_id) then
          coarse_row = iel
          subelem_slot = jel
          if (allocated(coarse%myELEMLINK)) coarse_elem_global = coarse%myELEMLINK(iel)
          return
        end if
      end do
    end do

  end subroutine el_find_dump_coarse_mapping

  subroutine el_write_sampling_dumpmap_debug_record(istep, particle, rank_id, sample_kind, elem_id)

    integer, intent(in) :: istep, rank_id, elem_id
    type(tParticleData), intent(in) :: particle
    character(len=*), intent(in) :: sample_kind

    integer :: unit_id, ios
    integer :: coarse_row, coarse_elem_global, subelem_slot
    integer :: unit_map, unit_elem, ios_map, ios_elem, il
    logical :: file_exists
    character(len=128) :: filename
    logical :: file_exists_elem

    if (.not. el_write_sampling_dumpmap_debug) return
    if (.not. el_is_debug_particle(particle%uniqueIdx)) return
    if (elem_id .le. 0) return

    call el_find_dump_coarse_mapping(elem_id, coarse_row, coarse_elem_global, subelem_slot)

    write(filename,'(A,I4.4,A,I6.6,A)') 'el_frozen_sampling_dumpmap_entries_rank', rank_id, '_step', istep, '.csv'
    inquire(file=trim(filename), exist=file_exists)
    open(newunit=unit_map, file=trim(filename), status='unknown', position='append', action='write', iostat=ios_map)
    if (ios_map .ne. 0) return

    if (.not. file_exists) then
      write(unit_map,'(A)') 'id,system_id,rank,sample_kind,elem_id,coarse_row,coarse_elem_global,subelem_slot,map_slot,map_dof'
    end if

    if (coarse_row .gt. 0 .and. allocated(myDump%Vertices)) then
      do il = 1, size(myDump%Vertices, 2)
        write(unit_map,'(3(I0,","),A,",",6(I0,","),I0)') particle%uniqueIdx, particle%systemIdx, rank_id, trim(sample_kind), &
          elem_id, coarse_row, coarse_elem_global, subelem_slot, il, myDump%Vertices(coarse_row, il)
      end do
    end if

    close(unit_map)

    write(filename,'(A,I4.4,A,I6.6,A)') 'el_frozen_sampling_dumpelem_rank', rank_id, '_step', istep, '.csv'
    inquire(file=trim(filename), exist=file_exists_elem)
    open(newunit=unit_elem, file=trim(filename), status='unknown', position='append', action='write', iostat=ios_elem)
    if (ios_elem .ne. 0) return

    if (.not. file_exists_elem) then
      write(unit_elem,'(A)') 'id,system_id,rank,sample_kind,elem_id,coarse_row,coarse_elem_global,subelem_slot,elem_slot,fine_elem'
    end if

    if (coarse_row .gt. 0 .and. allocated(myDump%Elements)) then
      do il = 1, size(myDump%Elements, 2)
        write(unit_elem,'(3(I0,","),A,",",6(I0,","),I0)') particle%uniqueIdx, particle%systemIdx, rank_id, trim(sample_kind), &
          elem_id, coarse_row, coarse_elem_global, subelem_slot, il, myDump%Elements(coarse_row, il)
      end do
    end if

    close(unit_elem)

  end subroutine el_write_sampling_dumpmap_debug_record

  subroutine el_write_sampling_dof_sources_debug_record(istep, particle, rank_id, sample_kind, elem_id, direct_q2)

    integer, intent(in) :: istep, rank_id, elem_id
    type(tParticleData), intent(in) :: particle
    character(len=*), intent(in) :: sample_kind
    integer, intent(in) :: direct_q2(27)

    integer :: unit_id, ios
    integer :: idof, irow, islot, isource
    logical :: file_exists
    character(len=128) :: filename

    if (.not. el_write_sampling_dof_sources_debug) return
    if (.not. el_is_debug_particle(particle%uniqueIdx)) return
    if (elem_id .le. 0) return
    if (.not. allocated(myDump%Vertices)) return

    write(filename,'(A,I4.4,A,I6.6,A)') 'el_frozen_sampling_dof_sources_rank', rank_id, '_step', istep, '.csv'
    inquire(file=trim(filename), exist=file_exists)
    open(newunit=unit_id, file=trim(filename), status='unknown', position='append', action='write', iostat=ios)
    if (ios .ne. 0) return

    if (.not. file_exists) then
      write(unit_id,'(A)') 'id,system_id,rank,sample_kind,elem_id,global_dof,source_idx,coarse_row,coarse_elem_global,map_slot'
    end if

    do idof = 1, size(direct_q2)
      if (direct_q2(idof) .le. 0) cycle
      isource = 0
      do irow = 1, size(myDump%Vertices, 1)
        do islot = 1, size(myDump%Vertices, 2)
          if (myDump%Vertices(irow, islot) .eq. direct_q2(idof)) then
            isource = isource + 1
            write(unit_id,'(3(I0,","),A,",",6(I0,","),I0)') particle%uniqueIdx, particle%systemIdx, rank_id, trim(sample_kind), &
              elem_id, direct_q2(idof), isource, irow, coarse%myELEMLINK(irow), islot
          end if
        end do
      end do
    end do

    close(unit_id)

  end subroutine el_write_sampling_dof_sources_debug_record

  subroutine el_get_closure_capabilities(caps)

    type(tElClosureCapabilities), intent(out) :: caps
    integer :: kernel_id

    caps = tElClosureCapabilities()
    kernel_id = el_kernel_id(trim(el_force_kernel))

    select case (kernel_id)
    case (EL_KERNEL_NONE)
      continue
    case (EL_KERNEL_TRACER)
      caps%supports_slip_velocity = .true.
    case (EL_KERNEL_STOKES, EL_KERNEL_SCHILLER_NAUMANN)
      caps%supports_force = .true.
      caps%supports_slip_velocity = .true.
      caps%supports_reynolds = .true.
    case default
      continue
    end select

  end subroutine el_get_closure_capabilities

  subroutine el_evaluate_closure(particle, carrier_vel, found_count, force, torque, slip, re_p)

    type(tParticleData), intent(in) :: particle
    real*8, intent(in) :: carrier_vel(3)
    integer, intent(in) :: found_count
    real*8, intent(out) :: force(3), torque(3), slip(3), re_p

    integer :: kernel_id
    real*8 :: mu, diameter, slip_norm, correction, particle_volume

    force = 0.0d0
    torque = 0.0d0
    slip = 0.0d0
    re_p = 0.0d0

    if (found_count .le. 0) return

    kernel_id = el_kernel_id(trim(el_force_kernel))
    slip = carrier_vel - particle%velocity

    select case (kernel_id)
    case (EL_KERNEL_NONE)
      return
    case (EL_KERNEL_TRACER)
      return
    case (EL_KERNEL_STOKES, EL_KERNEL_SCHILLER_NAUMANN)
      mu = el_fluid_density * el_kinematic_viscosity
      diameter = el_particle_diameter(particle)
      if (mu .le. 0.0d0 .or. diameter .le. 0.0d0) return

      slip_norm = sqrt(sum(slip*slip))
      if (slip_norm .gt. 0.0d0) re_p = el_fluid_density * slip_norm * diameter / mu

      correction = 1.0d0
      if (kernel_id .eq. EL_KERNEL_SCHILLER_NAUMANN) then
        if (re_p .gt. 0.0d0 .and. re_p .lt. 1000.0d0) then
          correction = 1.0d0 + 0.15d0 * re_p**0.687d0
        else if (re_p .ge. 1000.0d0) then
          correction = 0.44d0 * re_p / 24.0d0
        end if
      end if

      force = 3.0d0 * acos(-1.0d0) * mu * diameter * correction * slip

      if (el_enable_buoyancy) then
        particle_volume = 4.0d0/3.0d0 * acos(-1.0d0) * particle%radius**3
        force = force + particle_volume * (el_particle_density - el_fluid_density) * el_gravity
      end if
    case default
      return
    end select

  end subroutine el_evaluate_closure

  subroutine el_print_pass_summary(log_unit, num_particles, found_count, force, re_p)

    integer, intent(in) :: log_unit, num_particles
    integer, intent(in) :: found_count(:)
    real*8, intent(in) :: force(:, :), re_p(:)

    integer :: ip
    integer :: m, mt, mkeyb, mterm, merr, mprot, msys, mtrc, irecl8
    real*8, allocatable :: force_mag(:)
    real*8 :: avg_re, min_re, max_re
    real*8 :: avg_force, min_force, max_force
    real*8 :: avg_found

    common /OUTPUT/ m, mt, mkeyb, mterm, merr, mprot, msys, mtrc, irecl8

    if (num_particles .le. 0) return

    allocate(force_mag(num_particles))
    do ip = 1, num_particles
      force_mag(ip) = sqrt(sum(force(:, ip)*force(:, ip)))
    end do

    avg_re = sum(re_p) / dble(num_particles)
    min_re = minval(re_p)
    max_re = maxval(re_p)

    avg_force = sum(force_mag) / dble(num_particles)
    min_force = minval(force_mag)
    max_force = maxval(force_mag)

    avg_found = sum(dble(found_count)) / dble(num_particles)

    write(mterm,'(A,I0)') '  particle count = ', num_particles
    write(log_unit,'(A,I0)') '  particle count = ', num_particles
    write(mterm,'(A,3(ES14.6,1X))') '  Re_p    avg/min/max = ', avg_re, min_re, max_re
    write(log_unit,'(A,3(ES14.6,1X))') '  Re_p    avg/min/max = ', avg_re, min_re, max_re
    write(mterm,'(A,3(ES14.6,1X))') '  |F|     avg/min/max = ', avg_force, min_force, max_force
    write(log_unit,'(A,3(ES14.6,1X))') '  |F|     avg/min/max = ', avg_force, min_force, max_force
    write(mterm,'(A,ES14.6)') '  found_count avg = ', avg_found
    write(log_unit,'(A,ES14.6)') '  found_count avg = ', avg_found

    deallocate(force_mag)

  end subroutine el_print_pass_summary

  subroutine el_write_particle_diagnostics(istep, particles, num_particles, found_count, elem_sum, owner_sum, &
                                           carrier_vel, slip, force, re_p, caps)

    integer, intent(in) :: istep, num_particles
    type(tParticleData), intent(in) :: particles(:)
    integer, intent(in) :: found_count(:), elem_sum(:), owner_sum(:)
    real*8, intent(in) :: carrier_vel(:, :), slip(:, :), force(:, :), re_p(:)
    type(tElClosureCapabilities), intent(in) :: caps

    integer :: ip, unit_id, ios
    character(len=128) :: filename
    integer :: m, mt, mkeyb, mterm, merr, mprot, msys, mtrc, irecl8

    common /OUTPUT/ m, mt, mkeyb, mterm, merr, mprot, msys, mtrc, irecl8

    if (.not. el_write_diagnostics) return

    write(filename,'(A,I6.6,A)') 'el_frozen_particles_step', istep, '.csv'
    open(newunit=unit_id, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios .ne. 0) then
      write(mterm,'(A,A)') 'Failed to open frozen particle diagnostics: ', trim(filename)
      return
    end if

    write(unit_id,'(A)') 'id,system_id,x,y,z,radius,pvx,pvy,pvz,ux,uy,uz,slipx,slipy,slipz,fx,fy,fz,re_p,found_count,owner_rank_sum,elem_sum'
    do ip = 1, num_particles
      write(unit_id,'(2(I0,","),17(ES18.10,","),I0,",",I0,",",I0)') particles(ip)%uniqueIdx, &
        particles(ip)%systemIdx, particles(ip)%position(1), particles(ip)%position(2), particles(ip)%position(3), &
        particles(ip)%radius, particles(ip)%velocity(1), particles(ip)%velocity(2), particles(ip)%velocity(3), &
        carrier_vel(1, ip), carrier_vel(2, ip), carrier_vel(3, ip), &
        slip(1, ip), slip(2, ip), slip(3, ip), &
        force(1, ip), force(2, ip), force(3, ip), &
        re_p(ip), found_count(ip), owner_sum(ip), elem_sum(ip)
    end do
    close(unit_id)

    if (caps%supports_force) then
      write(mterm,'(A,A)') '  wrote particle diagnostics to ', trim(filename)
    else
      write(mterm,'(A,A)') '  wrote sampling diagnostics to ', trim(filename)
    end if

  end subroutine el_write_particle_diagnostics

  subroutine el_write_particle_diagnostics_parallel(istep, particles, num_particles, found_count, elem_sum, owner_sum, &
                                                    carrier_vel, slip, force, re_p, caps)

    integer, intent(in) :: istep, num_particles
    type(tParticleData), intent(in) :: particles(:)
    integer, intent(in) :: found_count(:), elem_sum(:), owner_sum(:)
    real*8, intent(in) :: carrier_vel(:, :), slip(:, :), force(:, :), re_p(:)
    type(tElClosureCapabilities), intent(in) :: caps

    integer :: ip, unit_id, ios
    character(len=128) :: filename
    integer :: m, mt, mkeyb, mterm, merr, mprot, msys, mtrc, irecl8

    common /OUTPUT/ m, mt, mkeyb, mterm, merr, mprot, msys, mtrc, irecl8

    if (.not. el_write_diagnostics) return
    if (num_particles .le. 0) return

    write(filename,'(A,I4.4,A,I6.6,A)') 'el_frozen_particles_rank', myid, '_step', istep, '.csv'
    open(newunit=unit_id, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios .ne. 0) then
      write(mterm,'(A,A)') 'Failed to open frozen particle diagnostics: ', trim(filename)
      return
    end if

    write(unit_id,'(A)') 'id,system_id,x,y,z,radius,pvx,pvy,pvz,ux,uy,uz,slipx,slipy,slipz,fx,fy,fz,re_p,found_count,owner_rank_sum,elem_sum'
    do ip = 1, num_particles
      write(unit_id,'(2(I0,","),17(ES18.10,","),I0,",",I0,",",I0)') particles(ip)%uniqueIdx, &
        particles(ip)%systemIdx, particles(ip)%position(1), particles(ip)%position(2), particles(ip)%position(3), &
        particles(ip)%radius, particles(ip)%velocity(1), particles(ip)%velocity(2), particles(ip)%velocity(3), &
        carrier_vel(1, ip), carrier_vel(2, ip), carrier_vel(3, ip), &
        slip(1, ip), slip(2, ip), slip(3, ip), &
        force(1, ip), force(2, ip), force(3, ip), &
        re_p(ip), found_count(ip), owner_sum(ip), elem_sum(ip)
    end do
    close(unit_id)

    if (myid .eq. showid) then
      if (caps%supports_force) then
        write(mterm,'(A)') '  wrote rank-local particle diagnostics el_frozen_particles_rankXXXX_stepXXXXXX.csv'
      else
        write(mterm,'(A)') '  wrote rank-local sampling diagnostics el_frozen_particles_rankXXXX_stepXXXXXX.csv'
      end if
    end if

  end subroutine el_write_particle_diagnostics_parallel

  subroutine el_print_configuration(log_unit)

    integer, intent(in) :: log_unit
    type(tElClosureCapabilities) :: caps
    integer :: m, mt, mkeyb, mterm, merr, mprot, msys, mtrc, irecl8

    common /OUTPUT/ m, mt, mkeyb, mterm, merr, mprot, msys, mtrc, irecl8

    call el_get_closure_capabilities(caps)

    if (myid .eq. showid) then
      write(mterm,'(A)') 'Frozen-field E-L configuration:'
      write(log_unit,'(A)') 'Frozen-field E-L configuration:'
      write(mterm,'(A,A)') '  ELForceKernel = ', trim(el_force_kernel)
      write(log_unit,'(A,A)') '  ELForceKernel = ', trim(el_force_kernel)
      write(mterm,'(A,L1)') '  ELWriteDiagnostics = ', el_write_diagnostics
      write(log_unit,'(A,L1)') '  ELWriteDiagnostics = ', el_write_diagnostics
      write(mterm,'(A,L1)') '  ELWriteSamplingDebug = ', el_sampling_debug_enabled
      write(log_unit,'(A,L1)') '  ELWriteSamplingDebug = ', el_sampling_debug_enabled
      write(mterm,'(A,L1)') '  ELWriteSamplingSummaryDebug = ', el_write_sampling_summary_debug
      write(log_unit,'(A,L1)') '  ELWriteSamplingSummaryDebug = ', el_write_sampling_summary_debug
      write(mterm,'(A,L1)') '  ELWriteSamplingDofDebug = ', el_write_sampling_dof_debug
      write(log_unit,'(A,L1)') '  ELWriteSamplingDofDebug = ', el_write_sampling_dof_debug
      write(mterm,'(A,L1)') '  ELWriteSamplingConnectivityDebug = ', el_write_sampling_connectivity_debug
      write(log_unit,'(A,L1)') '  ELWriteSamplingConnectivityDebug = ', el_write_sampling_connectivity_debug
      write(mterm,'(A,L1)') '  ELWriteSamplingTopologyDebug = ', el_write_sampling_topology_debug
      write(log_unit,'(A,L1)') '  ELWriteSamplingTopologyDebug = ', el_write_sampling_topology_debug
      write(mterm,'(A,L1)') '  ELWriteSamplingFieldDebug = ', el_write_sampling_field_debug
      write(log_unit,'(A,L1)') '  ELWriteSamplingFieldDebug = ', el_write_sampling_field_debug
      write(mterm,'(A,L1)') '  ELWriteSamplingDumpmapDebug = ', el_write_sampling_dumpmap_debug
      write(log_unit,'(A,L1)') '  ELWriteSamplingDumpmapDebug = ', el_write_sampling_dumpmap_debug
      write(mterm,'(A,L1)') '  ELWriteSamplingDofSourcesDebug = ', el_write_sampling_dof_sources_debug
      write(log_unit,'(A,L1)') '  ELWriteSamplingDofSourcesDebug = ', el_write_sampling_dof_sources_debug
      write(mterm,'(A,L1)') '  ELApplyForces = ', el_apply_forces
      write(log_unit,'(A,L1)') '  ELApplyForces = ', el_apply_forces
      write(mterm,'(A,L1)') '  ELEnableBuoyancy = ', el_enable_buoyancy
      write(log_unit,'(A,L1)') '  ELEnableBuoyancy = ', el_enable_buoyancy
      write(mterm,'(A,ES12.4)') '  ELFluidDensity = ', el_fluid_density
      write(log_unit,'(A,ES12.4)') '  ELFluidDensity = ', el_fluid_density
      write(mterm,'(A,ES12.4)') '  ELKinematicViscosity = ', el_kinematic_viscosity
      write(log_unit,'(A,ES12.4)') '  ELKinematicViscosity = ', el_kinematic_viscosity
      write(mterm,'(A,ES12.4)') '  ELParticleDensity = ', el_particle_density
      write(log_unit,'(A,ES12.4)') '  ELParticleDensity = ', el_particle_density
      write(mterm,'(A,3(ES12.4,1X))') '  ELGravity = ', el_gravity
      write(log_unit,'(A,3(ES12.4,1X))') '  ELGravity = ', el_gravity
      write(mterm,'(A,I0)') '  ELDebugIDCount = ', el_num_debug_ids
      write(log_unit,'(A,I0)') '  ELDebugIDCount = ', el_num_debug_ids
      write(mterm,'(A,L1,A,L1,A,L1)') '  capabilities: force=', caps%supports_force, &
        ', torque=', caps%supports_torque, ', reynolds=', caps%supports_reynolds
      write(log_unit,'(A,L1,A,L1,A,L1)') '  capabilities: force=', caps%supports_force, &
        ', torque=', caps%supports_torque, ', reynolds=', caps%supports_reynolds
    end if

  end subroutine el_print_configuration

  integer function el_kernel_id(name)

    character(len=*), intent(in) :: name
    character(len=32) :: kernel_name

    kernel_name = el_normalize_string(name)
    select case (trim(kernel_name))
    case ('none')
      el_kernel_id = EL_KERNEL_NONE
    case ('tracer')
      el_kernel_id = EL_KERNEL_TRACER
    case ('stokes_drag')
      el_kernel_id = EL_KERNEL_STOKES
    case ('schiller_naumann')
      el_kernel_id = EL_KERNEL_SCHILLER_NAUMANN
    case default
      el_kernel_id = EL_KERNEL_NONE
    end select

  end function el_kernel_id

  real*8 function el_particle_diameter(particle)

    type(tParticleData), intent(in) :: particle
    integer(c_int), parameter :: PE_SPHERE_TYPE = 1_c_int

    if (particle%typeId .eq. PE_SPHERE_TYPE .and. particle%radius .gt. 0d0) then
      el_particle_diameter = 2d0 * particle%radius
    else
      el_particle_diameter = maxval(particle%aabb)
    end if

  end function el_particle_diameter

  logical function el_parse_yes_no(value)

    character(len=*), intent(in) :: value
    character(len=32) :: normalized

    normalized = el_normalize_string(value)
    el_parse_yes_no = (trim(normalized) .eq. 'yes' .or. trim(normalized) .eq. 'true' .or. &
                       trim(normalized) .eq. 'on' .or. trim(normalized) .eq. '1')

  end function el_parse_yes_no

  character(len=32) function el_normalize_string(value)

    character(len=*), intent(in) :: value
    integer :: i, ascii_code

    el_normalize_string = adjustl(trim(value))
    do i = 1, len_trim(el_normalize_string)
      ascii_code = iachar(el_normalize_string(i:i))
      if (ascii_code .ge. iachar('A') .and. ascii_code .le. iachar('Z')) then
        el_normalize_string(i:i) = achar(ascii_code + 32)
      end if
    end do

  end function el_normalize_string

  subroutine el_set_debug_ids(value)

    character(len=*), intent(in) :: value
    character(len=512) :: work
    integer :: i, start_pos, end_pos, comma_pos

    el_num_debug_ids = 0
    el_debug_ids = 0

    work = trim(adjustl(value))
    if (len_trim(work) .le. 0) return

    start_pos = 1
    do
      comma_pos = index(work(start_pos:), ',')
      if (comma_pos .gt. 0) then
        end_pos = start_pos + comma_pos - 2
      else
        end_pos = len_trim(work)
      end if

      if (end_pos .ge. start_pos) then
        if (el_num_debug_ids .lt. EL_MAX_DEBUG_IDS) then
          el_num_debug_ids = el_num_debug_ids + 1
          read(work(start_pos:end_pos), *) el_debug_ids(el_num_debug_ids)
        end if
      end if

      if (comma_pos .le. 0) exit
      start_pos = end_pos + 2
      if (start_pos .gt. len_trim(work)) exit
    end do

  end subroutine el_set_debug_ids

end module el_frozen_driver
