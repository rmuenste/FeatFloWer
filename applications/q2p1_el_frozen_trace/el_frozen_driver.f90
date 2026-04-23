module el_frozen_driver

  use iso_c_binding, only: c_int
  use PP3D_MPI, only: myid, showid
  use var_QuadScalar, only: mg_mesh
  use Transport_Q2P1, only: QuadSc
  use dem_query, only: tParticleData, getAllParticles, numTotalParticles, setForcesMapped
  use fbmaux, only: fbmaux_PointInHex

  implicit none

  integer, parameter :: EL_KERNEL_NONE = 0
  integer, parameter :: EL_KERNEL_TRACER = 1
  integer, parameter :: EL_KERNEL_STOKES = 2
  integer, parameter :: EL_KERNEL_SCHILLER_NAUMANN = 3

  type tElClosureCapabilities
    logical :: supports_force = .false.
    logical :: supports_torque = .false.
    logical :: supports_pressure = .false.
    logical :: supports_slip_velocity = .false.
    logical :: supports_reynolds = .false.
  end type tElClosureCapabilities

  character(len=32) :: el_force_kernel = 'none'
  logical :: el_write_diagnostics = .true.
  logical :: el_apply_forces = .true.
  logical :: el_enable_buoyancy = .false.
  real*8 :: el_fluid_density = 1.0d0
  real*8 :: el_kinematic_viscosity = 1.0d-3
  real*8 :: el_particle_density = 1.0d0
  real*8 :: el_gravity(3) = 0.0d0

  interface
    subroutine step_el_frozen_trace_c() bind(C, name="step_el_frozen_trace_")
    end subroutine step_el_frozen_trace_c
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

    common /OUTPUT/ m, mt, mkeyb, mterm, merr, mprot, msys, mtrc, irecl8

    include 'mpif.h'

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
      call el_setup_q2_sampling()

      do ip = 1, num_particles
        local_elem = 0
        call el_sample_particle_velocity_local(particles(ip), ilev, nel_local, vel_local(:, ip), local_elem)
        if (local_elem .gt. 0) then
          found_local(ip) = 1
          elem_local(ip) = local_elem
          owner_local(ip) = myid
        end if
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

  subroutine el_setup_q2_sampling()

    integer, parameter :: NNBAS = 27, NNDER = 10, NNVE = 8, NNDIM = 3
    integer :: ieltyp, i
    real*8 :: dx(nnve), dy(nnve), dz(nnve), djac(3, 3), detj
    real*8 :: dbas(nndim, nnbas, nnder)
    logical :: bder(nnder)
    integer :: kve(nnve), ndim, iel
    integer :: kdfg(nnbas), kdfl(nnbas), idfl

    common /ELEM/ dx, dy, dz, djac, detj, dbas, bder, kve, iel, ndim
    common /COAUX1/ kdfg, kdfl, idfl

    integer, external :: ndfl
    external :: e013, setlev

    call setlev(2)
    do i = 1, nnder
      bder(i) = .false.
    end do
    bder(1) = .true.

    ieltyp = -1
    call e013(0d0, 0d0, 0d0, ieltyp)
    idfl = ndfl(ieltyp)

  end subroutine el_setup_q2_sampling

  subroutine el_sample_particle_velocity_local(particle, ilev, nel_local, vel_sample, elem_id)

    type(tParticleData), intent(in) :: particle
    integer, intent(in) :: ilev, nel_local
    real*8, intent(out) :: vel_sample(3)
    integer, intent(out) :: elem_id

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
        vel_sample(1) = vel_sample(1) + QuadSc%ValU(ig) * dbas(1, ive, 1)
        vel_sample(2) = vel_sample(2) + QuadSc%ValV(ig) * dbas(1, ive, 1)
        vel_sample(3) = vel_sample(3) + QuadSc%ValW(ig) * dbas(1, ive, 1)
      end do

      elem_id = iel
      return
    end do

  end subroutine el_sample_particle_velocity_local

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

    write(unit_id,'(A)') 'id,x,y,z,radius,pvx,pvy,pvz,ux,uy,uz,slipx,slipy,slipz,fx,fy,fz,re_p,found_count,owner_rank_sum,elem_sum'
    do ip = 1, num_particles
      write(unit_id,'(I0,17(",",ES18.10),3(",",I0))') particles(ip)%uniqueIdx, &
        particles(ip)%position(1), particles(ip)%position(2), particles(ip)%position(3), &
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

end module el_frozen_driver
