!=========================================================================
! QuadSc_el_bedwin.f90
!
! Window audit of fixed and fluidised beds (SimPar@ELBedWindowAudit; EL_BED_AUDIT writes the
! lines and documents them): fluid-side adapter. Included in Transport_Q2P1 (host association to
! QuadSc, LinSc, mg_mesh, Properties, E013Sum). One call after the bed audit of the step, on every
! worker; a no-op unless the deck key is set and the step is a bed-audit step.
!
! The vertex mode: the vertical gas velocity w over the Q2 dofs whose height lies in the interior
! window (EL_BED_AUDIT: BW_VM_LO .. BW_VM_HI above zmin), per dof class (vertex, edge, face, centre,
! the layout of EL_VA_ELEMENT_DOFS), each dof counted once over the ranks with the weight 1/m, m the
! number of its copies (a parallel sum of ones on the finest level, formed once per run on every
! worker alike). Read-only.
!=========================================================================
SUBROUTINE EL_BedWin_Adapter(mfile, istep)

  USE EL_BED_AUDIT, ONLY: EL_WRITE_BED_WINDOWS, BW_VM_LO, BW_VM_HI
  USE EL_GEOMETRY, ONLY: el_domain_box
  USE EL_CONFIG, ONLY: el_bed_window_audit
  USE PP3D_MPI, ONLY: myid, master

  INTEGER, INTENT(IN) :: mfile, istep
  REAL*8, ALLOCATABLE, SAVE :: inv_mult(:)
  LOGICAL, SAVE :: mult_ready = .FALSE.
  REAL*8, ALLOCATABLE :: ones(:)
  REAL*8 :: vm_sum(3,4), vm_min(4), vm_max(4), zl, zh, tol, w, z
  INTEGER :: n, i, c, nv, ne, na, ilev_save

  IF (.NOT.el_bed_window_audit) RETURN
  IF (myid.EQ.master) RETURN
  IF (el_bed_audit_freq.LE.0) RETURN
  IF (MOD(ABS(istep), el_bed_audit_freq).NE.0) RETURN

  n = QuadSc%ndof
  ! copies of every dof: once per run, on every worker at the same (first audit) call
  IF (.NOT.mult_ready) THEN
    ALLOCATE(ones(n), inv_mult(n))
    ones = 1.0d0
    ilev_save = ILEV
    ILEV = NLMAX
    CALL E013Sum(ones)
    ILEV = ilev_save
    inv_mult = 1.0d0/ones
    DEALLOCATE(ones)
    mult_ready = .TRUE.
  END IF

  nv = mg_mesh%level(NLMAX)%nvt
  ne = mg_mesh%level(NLMAX)%net
  na = mg_mesh%level(NLMAX)%nat
  zl = el_domain_box(5) + BW_VM_LO
  zh = el_domain_box(5) + BW_VM_HI
  tol = 1.0d-9*(el_domain_box(6) - el_domain_box(5))
  vm_sum = 0.0d0
  vm_min = HUGE(1.0d0)
  vm_max = -HUGE(1.0d0)
  DO i = 1, n
    z = mg_mesh%level(NLMAX)%dcorvg(3, i)
    IF (z.LT.zl - tol .OR. z.GT.zh + tol) CYCLE
    IF (i.LE.nv) THEN
      c = 1
    ELSE IF (i.LE.nv + ne) THEN
      c = 2
    ELSE IF (i.LE.nv + ne + na) THEN
      c = 3
    ELSE
      c = 4
    END IF
    w = QuadSc%valW(i)
    vm_sum(1,c) = vm_sum(1,c) + inv_mult(i)
    vm_sum(2,c) = vm_sum(2,c) + inv_mult(i)*w
    vm_sum(3,c) = vm_sum(3,c) + inv_mult(i)*w*w
    vm_min(c) = MIN(vm_min(c), w)
    vm_max(c) = MAX(vm_max(c), w)
  END DO

  CALL EL_WRITE_BED_WINDOWS(LinSc%valP(NLMAX)%x, timens, mfile, istep, mg_mesh, NLMAX, &
    Properties%Density(1), EL_CLOSURE_VISCOSITY(Properties%Viscosity(1)), vm_sum, vm_min, vm_max)

END SUBROUTINE EL_BedWin_Adapter
