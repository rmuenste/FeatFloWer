!##############################################################################
!# hashgrid_gate_driver
!#
!# Stage sequencing for the Layer-2 FBM classification gate. This module only
!# ORDERS calls: everything that decides a verdict lives in module fbm
!# (fbm_query_gate_reset / fbm_gate_preflight / fbm_query_gate_verdict), and
!# everything that knows the geometry - the HashGrid span, pe's contact
!# threshold, each body's real AABB - lives in pe
!# (libs/pe/src/interface/fbm_gate_probe.cpp, design in its header).
!#
!# Every routine here is called on EVERY rank, in the same order, from the
!# program's loop. That is what makes the collectives inside the verdict safe;
!# pe itself is only ever touched on the workers (the master has no pe world).
!##############################################################################
module hashgrid_gate_driver

  use iso_c_binding, only: c_int, c_double
  use PP3D_MPI, only: myid, MPI_COMM_WORLD
  use var_QuadScalar, only: cGateStageList, bGateRefreshIndex, bFbmOwnerTable, myErrorCode
  use Transport_Q2P1, only: updateFBMGeometry
  use fbm, only: fbm_query_gate_reset, fbm_query_gate_verdict, fbm_gate_preflight, &
                 FBM_GATE_MUST_PASS, FBM_GATE_MUST_DIFFER
  use dem_query, only: fbm_gate_place_stage, fbm_gate_nudge
  use fbm_ownership, only: fbm_own_report

  implicit none
  private
  public :: gate_preflight, gate_run_stages

  ! Stage codes and shape codes of fbm_gate_place_stage_ (pe/interface/fbm_gate_probe.h)
  integer(c_int), parameter :: STAGE_PARK = 0, STAGE_LATTICE = 1, STAGE_OVERLAP = 2, STAGE_PARTITION = 3
  integer(c_int), parameter :: ALL_SHAPES = -1
  integer, parameter :: NUM_SHAPES = 5
  character(len=9), parameter :: shapeName(0:NUM_SHAPES-1) = &
    [character(len=9) :: 'sphere', 'box', 'capsule', 'ellipsoid', 'cylinder']

  ! The nudge along x: ONE Q2 DOF PITCH of the QBOX9Z3 level-2 mesh (18 cells, pitch 1/36).
  !
  ! Layer 1 found the defect with a nudge of 0.035 % of a diameter because it chose its query
  ! points just inside the far surface. On the real mesh the query points are the Q2 DOFs, on
  ! planes x = m/36. A probe parked with its AABB minimum delta below plane k*span reaches to
  ! (k+1)*span - 3 delta; a small nudge leaves its far surface inside cell k, every contained
  ! DOF hashes to cell k, and the 27-cell stencil of cell k still covers the stale cell k-1 -
  ! no miss is possible. The first run of this harness proved that: A2 compared 53391 DOFs and
  ! found nothing. For a miss, a DOF plane BEYOND (k+1)*span must lie inside the body, and the
  ! smallest nudge that guarantees one is the DOF pitch. The minimum then sits in cell k
  ! (stale index says k-1), the far DOFs in cell k+1 with stencil {k, k+1, k+2}: a miss.
  ! Four of the five planes gain such a DOF plane inside the domain; the fifth lies beyond
  ! x = 1. If the mesh changes, A2 fails as a control and says so - which is the point.
  ! (named GATE_NUDGE, not NUDGE: Fortran is case-insensitive and there is a subroutine nudge)
  real(c_double), parameter :: GATE_NUDGE(3) = [1.0d0/36.0d0, 0.0d0, 0.0d0]

  ! Fewer compared classifications than this and a stage is judged vacuous. The lattice
  ! sweep compares ~17 800 DOFs per worker, so this only has to exclude "nothing ran".
  integer, parameter :: MIN_COMPARED = 1000

  ! pe's refresh_spatial_index_(); FeatFloWer calls it as an implicit external everywhere
  ! (gfortran's trailing-underscore mangling is the C name), there is no interface block.
  external :: refresh_spatial_index

contains

  subroutine gate_preflight(mfile)
    integer, intent(in) :: mfile
    call fbm_gate_preflight(mfile)
  end subroutine gate_preflight

  ! Re-file the spatial index - unless the deck turned that off, which makes every
  ! verdict stage a control (a debugging knob, see SimPar@GateRefreshIndex).
  subroutine refresh()
    if (myid /= 0 .and. bGateRefreshIndex) call refresh_spatial_index()
  end subroutine refresh

  subroutine place(stage, shapeFilter)
    integer(c_int), intent(in) :: stage, shapeFilter
    integer(c_int) :: nPlaced
    nPlaced = 0
    if (myid /= 0) call fbm_gate_place_stage(stage, shapeFilter, nPlaced)
  end subroutine place

  subroutine nudge()
    if (myid /= 0) call fbm_gate_nudge(GATE_NUDGE)
  end subroutine nudge

  ! The real classification sweep (QuadScalar_FictKnpr over every Q2 DOF, then
  ! E013Max_SUPER), followed by the collective verdict. Both are rank-symmetric.
  ! expectMulti: the stage was built to overlap bodies on real DOFs, so the ownership
  ! table must have rows - otherwise the overlap construction reached nothing.
  subroutine sweep_and_judge(simTime, name, mode, expectMulti)
    real*8,           intent(in) :: simTime
    character(len=*), intent(in) :: name
    integer,          intent(in) :: mode
    logical, optional, intent(in) :: expectMulti
    logical :: em
    em = .false.
    if (present(expectMulti)) em = expectMulti
    call updateFBMGeometry()
    call fbm_query_gate_verdict(simTime, name, mode, MIN_COMPARED)
    if (bFbmOwnerTable) call fbm_own_report(simTime, name, em)
  end subroutine sweep_and_judge

  subroutine run_one(simTime, token)
    real*8,           intent(in) :: simTime
    character(len=*), intent(in) :: token
    integer :: s, ierr

    select case (trim(token))

    case ('A1')
      ! Baseline: every probe on the lattice, AABB minima delta below the planes, index
      ! current. The move from the parking positions is large, so this also proves the
      ! refresh handles a body changing cells by a long way.
      call fbm_query_gate_reset(); call place(STAGE_LATTICE, ALL_SHAPES); call refresh()
      call sweep_and_judge(simTime, 'A1_lattice', FBM_GATE_MUST_PASS)

    case ('A2')
      ! Teeth check: the same construction nudged one DOF pitch past the planes, NO refresh. If this
      ! does not differ, the placement missed the lattice and nothing below means anything.
      call fbm_query_gate_reset(); call place(STAGE_LATTICE, ALL_SHAPES); call refresh()
      call nudge()
      call sweep_and_judge(simTime, 'A2_teeth', FBM_GATE_MUST_DIFFER)

    case ('G')
      ! The adversarial construction, one shape at a time, stale index: the control pass.
      do s = 0, NUM_SHAPES-1
        call fbm_query_gate_reset(); call place(STAGE_LATTICE, int(s, c_int)); call refresh()
        call nudge()
        call sweep_and_judge(simTime, 'G_'//trim(shapeName(s)), FBM_GATE_MUST_DIFFER)
      end do

    case ('R')
      ! The same construction with the index refreshed after the nudge: the verdict pass.
      do s = 0, NUM_SHAPES-1
        call fbm_query_gate_reset(); call place(STAGE_LATTICE, int(s, c_int)); call refresh()
        call nudge(); call refresh()
        call sweep_and_judge(simTime, 'R_'//trim(shapeName(s)), FBM_GATE_MUST_PASS)
      end do

    case ('E')
      ! Overlapping pairs across a plane. Both query paths own an overlap by the minimum
      ! system id since pe 8ea90af, so the full 8-byte id must agree on every DOF - and with
      ! the ownership table on, the full owner SETS must agree and the table must be non-empty.
      call fbm_query_gate_reset(); call place(STAGE_OVERLAP, ALL_SHAPES); call refresh()
      call sweep_and_judge(simTime, 'E_overlap', FBM_GATE_MUST_PASS, expectMulti=.true.)

    case ('P')
      ! Probes centred on the z = 1/3 and z = 2/3 subdomain cuts: one body's interior DOFs
      ! are classified by two ranks. Only a real partitioned mesh can ask this.
      call fbm_query_gate_reset(); call place(STAGE_PARTITION, ALL_SHAPES); call refresh()
      call sweep_and_judge(simTime, 'P_partition', FBM_GATE_MUST_PASS)

    case default
      ! A typo must not become a silently skipped stage - that is the vacuity trap again.
      if (myid == 1) then
        write(*,'(A,A,A)') 'DNS_QUERY_GATE FAILED: unknown stage "', trim(token), &
          '" in SimPar@GateStageList (known: A1 A2 G R E P)'
        flush(6)
      end if
      call MPI_Abort(MPI_COMM_WORLD, myErrorCode%PARAMETER_CONFLICT, ierr)
    end select
  end subroutine run_one

  ! Run the comma-separated stage list from the deck, in order.
  subroutine gate_run_stages(simTime, mfile)
    real*8,  intent(in) :: simTime
    integer, intent(in) :: mfile
    character(len=64) :: rest
    character(len=16) :: token
    integer :: ic, nStages

    rest = adjustl(cGateStageList)
    nStages = 0
    do while (len_trim(rest) > 0)
      ic = index(rest, ',')
      if (ic == 0) then
        token = adjustl(rest); rest = ''
      else
        token = adjustl(rest(1:ic-1)); rest = adjustl(rest(ic+1:))
      end if
      if (len_trim(token) == 0) cycle
      nStages = nStages + 1
      if (myid == 1) write(mfile,'(A,A)') 'gate stage: ', trim(token)
      call run_one(simTime, token)
    end do

    if (myid == 1) then
      write(*,'(A,I0,A)') 'DNS_QUERY_GATE_SUMMARY stages_judged= ', nStages, ' verdict= green'
      flush(6)
    end if
  end subroutine gate_run_stages

end module hashgrid_gate_driver
