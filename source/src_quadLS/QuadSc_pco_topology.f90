!=========================================================================
! QuadSc_pco_topology.f90 - periodic topology of the coarse pressure level of the master
!
! Part of the periodic coarse pressure operator (overview: QuadSc_pco_operator.f90). This
! module answers two questions and owns no solver state:
!
! 1. Is periodicity active? (PCO_LINK_CENSUS, collective)
!    A finite period alone does not say so (q2p1_dns_drag sets unit periods on every box).
!    The evidence is the provenance of the couplings the workers apply: the communication
!    set-up (E013_CreateComm*, QuadSc_mpi.f90) accepts a velocity dof link between two ranks
!    when the coordinates coincide or lie one period apart, axis by axis, and counts the
!    accepted links that matched through the period of axis a (MGE013(lev)%nPerLink). Every
!    periodic coupling of the worker operator goes through such links (the pressure ghost
!    elements of CP are the elements at the linked dofs, ParPresComm_Init). Axis a is
!    active when the global count is positive. This is not re-derived from geometry: an
!    element pair that shares an ordinary vertex AND a periodic image (two coarse cells along
!    a periodic axis, one per rank) still contributes its periodic links, and a pair joined
!    through several images contributes every axis. The master then has to find the images
!    of every active axis in its classes (checked by the coordinator).
!    In addition every coupled element pair of CP is checked for consistency: the two
!    elements must share a vertex unshifted or after some shift s P (PCO_PAIR_SHIFT); a pair
!    that shares none in any image (nbad) means periods that do not fit the mesh.
!
! 2. Which dofs of the unpartitioned coarse mesh are periodic images of each other?
!    (PCO_BUILD_CLASSES, master only)
!    Velocity (Q2) dofs of level lev, numbered vertices, edges, faces, elements, with the
!    coordinates of level lev+1 (where they are vertices). Partners are searched with every
!    combined shift s P, s in {-1,0,1}^3 on the active axes (a vertex on a periodic edge has
!    3 images, in a triply periodic corner 7), by a hash of the coordinates with the tolerance
!    the workers use for their links (DEpsPrec). The partner relation is closed by union-find
!    into classes. Refused (counted, reported by the caller):
!      ambiguous  more than one dof within the tolerance of one shifted probe
!      incomplete a class whose members span m periodic axes must have 2^m members and every
!                 member 2^m - 1 direct partners
!      mixed      a class mixing vertex, edge, face or element dofs
!    The planes are not assumed to be the faces of a bounding box: any surface pair related
!    by the translation P is found, a missing partner on a periodic plane shows up later as a
!    worker coupling outside the master structure (pco_gather) or as an incomplete class.
!
! Also here: the conflict count of the velocity constraint masks over a class
! (PCO_CLASS_MASK_CONFLICTS) and the element graph through vertex classes
! (PCO_ELEMENT_GRAPH: elements sharing a vertex or one of its periodic images), which is
! the block pattern of the periodic coarse pressure matrix.
!=========================================================================
MODULE pco_periodic_topology

  USE PP3D_MPI, ONLY: myid, master, MPI_COMM_WORLD, MPI_INTEGER, MPI_SUM, MPI_LOGICAL, MPI_LAND, &
                      dPeriodicity, DEpsPrec, MGE013
  USE var_QuadScalar, ONLY: mg_mesh, mg_lPMat
  USE def_FEAT, ONLY: ILEV

  IMPLICIT NONE
  PRIVATE

  ! dof classes of one mesh level; singletons are classes too
  TYPE, PUBLIC :: tPcoClasses
    INTEGER :: ndof = 0, nvt = 0, ncls = 0, nmulti = 0, maxsize = 0
    INTEGER, ALLOCATABLE :: cls_of(:)        ! dof -> class
    INTEGER, ALLOCATABLE :: ptr(:), mem(:)   ! members of class k: mem(ptr(k):ptr(k+1)-1), ascending
    ! refusal counters of the construction
    INTEGER :: nambiguous = 0, nincomplete = 0, nmixed = 0
    INTEGER :: first_bad = 0                 ! a dof of the first refused class (0: none)
    LOGICAL :: axis_found(3) = .FALSE.       ! some class spans the periodic axis a
  END TYPE tPcoClasses

  TYPE, PUBLIC :: tPcoLinkCensus
    INTEGER :: nshift(3) = 0                 ! accepted dof links through the period of axis a
    INTEGER :: nbad = 0                      ! couplings without a shared vertex in any image
    LOGICAL :: axis_active(3) = .FALSE.
    LOGICAL :: possible = .TRUE.             ! every worker has the pressure structure of lev
  END TYPE tPcoLinkCensus

  PUBLIC :: PCO_LINK_CENSUS, PCO_BUILD_CLASSES, PCO_CLASS_MASK_CONFLICTS, PCO_ELEMENT_GRAPH
  PUBLIC :: PCO_GHOST_BUFFER_SIZE, PCO_MASKS_DIFFER, PCO_MASTER_HAS_DOF_COORDINATES, PCO_PAIR_SHIFT

  ! coordinate hash for the partner search
  TYPE tPointHash
    REAL*8 :: h = 0d0, x0(3) = 0d0
    INTEGER :: n = 0
    INTEGER(8), ALLOCATABLE :: key(:)
    INTEGER, ALLOCATABLE :: idx(:)
  END TYPE tPointHash
  INTEGER(8), PARAMETER :: HB = 2_8**20      ! cells per axis of the hash key

CONTAINS

  !-----------------------------------------------------------------------
  ! Size of the ghost buffer GetParPressure fills on level lev (workers): four values per
  ! ghost element slot of the exchange lists, and at least every column of CP.
  !-----------------------------------------------------------------------
  INTEGER FUNCTION PCO_GHOST_BUFFER_SIZE(lev) RESULT(nbuf)
    INTEGER, INTENT(IN) :: lev
    INTEGER :: pjd, i
    nbuf = 4
    IF (myid.EQ.master) RETURN
    IF (ALLOCATED(MGE013(lev)%SP)) THEN
      DO pjd = 1, SIZE(MGE013(lev)%SP)
        DO i = 1, MGE013(lev)%SP(pjd)%Num
          nbuf = MAX(nbuf, 4*MGE013(lev)%SP(pjd)%VertLink(2,i))
        END DO
      END DO
    END IF
    IF (ALLOCATED(mg_lPMat(lev)%ColA)) THEN
      IF (SIZE(mg_lPMat(lev)%ColA).GT.0) nbuf = MAX(nbuf, MAXVAL(mg_lPMat(lev)%ColA))
    END IF
  END FUNCTION PCO_GHOST_BUFFER_SIZE

  !-----------------------------------------------------------------------
  ! Census of level lev. Collective (every rank; the master contributes nothing).
  ! Activity: nshift(a) = global number of accepted velocity dof links through the period of
  ! axis a on level lev (MGE013(lev)%nPerLink, counted by the communication set-up; see the
  ! module header). A walled mesh has none: adjacent cells of any width match coincident
  ! dofs, and a period that is not a distance between interface dofs matches nothing.
  ! Consistency: the workers exchange the 24 vertex coordinates of their elements (six
  ! pressure exchanges of four values each); an element pair coupled by a ghost column of CP
  ! that shares no vertex unshifted nor after any shift s P (s in {-1,0,1}^3 on the axes with
  ! a finite period) is counted in nbad.
  ! census%possible is .FALSE. on every rank when some worker has no pressure structure on
  ! level lev (then nothing is counted).
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_LINK_CENSUS(lev, census)
    INTEGER, INTENT(IN) :: lev
    TYPE(tPcoLinkCensus), INTENT(OUT) :: census
    REAL*8, ALLOCATABLE :: p(:), pp(:), vx(:,:), gx(:,:)
    INTEGER :: nel, iel, k, ia, ic, ierr, ilev_save, nbuf, chunk, g, shift(3), kind
    INTEGER :: nloc(4), nglob(4)
    LOGICAL :: can, can_all
    INTEGER, PARAMETER :: NV = 24

    can = .TRUE.
    IF (myid.NE.master) THEN
      can = ALLOCATED(mg_lPMat) .AND. ALLOCATED(mg_mesh%level)
      IF (can) can = (lev.GE.LBOUND(mg_lPMat,1) .AND. lev.LE.UBOUND(mg_lPMat,1) .AND. &
                      lev.LE.SIZE(mg_mesh%level))
      IF (can) can = ALLOCATED(mg_lPMat(lev)%LdA)
    END IF
    CALL MPI_Allreduce(can, can_all, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    census%possible = can_all
    IF (.NOT.can_all) RETURN

    nloc = 0
    ilev_save = ILEV
    IF (myid.NE.master) THEN
      nel = mg_mesh%level(lev)%nel
      nbuf = PCO_GHOST_BUFFER_SIZE(lev)
      ALLOCATE(p(4*nel), pp(nbuf), vx(NV,nel), gx(NV,(nbuf+3)/4))
      DO iel = 1, nel
        DO k = 1, 8
          vx(3*(k-1)+1:3*k, iel) = mg_mesh%level(lev)%dcorvg(:, mg_mesh%level(lev)%kvert(k,iel))
        END DO
      END DO
      gx = 0d0
    ELSE
      ALLOCATE(p(4), pp(4))
    END IF
    ILEV = lev
    DO chunk = 1, NV/4
      IF (myid.NE.master) THEN
        DO iel = 1, nel
          p(4*(iel-1)+1:4*iel) = vx(4*(chunk-1)+1:4*chunk, iel)
        END DO
        pp = 0d0
      END IF
      CALL GetParPressure(p, pp)
      IF (myid.NE.master) THEN
        DO g = 1, SIZE(gx,2)
          gx(4*(chunk-1)+1:4*chunk, g) = pp(4*(g-1)+1:4*g)
        END DO
      END IF
    END DO
    ILEV = ilev_save

    IF (myid.NE.master) THEN
      DO iel = 1, nel
        ic = 4*(iel-1)+1
        DO ia = mg_lPMat(lev)%LdA(ic), mg_lPMat(lev)%LdA(ic+1)-1
          k = mg_lPMat(lev)%ColA(ia)
          IF (MOD(k-1,4).NE.0) CYCLE            ! one entry per coupled ghost element
          g = (k-1)/4 + 1
          CALL PCO_PAIR_SHIFT(vx(:,iel), gx(:,g), dPeriodicity, DEpsPrec, shift, kind)
          IF (kind.LT.0) nloc(4) = nloc(4) + 1
        END DO
      END DO
      DEALLOCATE(vx, gx)
      ! the provenance of the periodic couplings (communication set-up)
      IF (ALLOCATED(MGE013)) THEN
        IF (lev.GE.LBOUND(MGE013,1) .AND. lev.LE.UBOUND(MGE013,1)) nloc(1:3) = MGE013(lev)%nPerLink
      END IF
    END IF
    DEALLOCATE(p, pp)

    CALL MPI_Allreduce(nloc, nglob, 4, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    census%nshift = nglob(1:3)
    census%nbad = nglob(4)
    census%axis_active = (census%nshift.GT.0)


  END SUBROUTINE PCO_LINK_CENSUS

  !-----------------------------------------------------------------------
  ! Consistency of a coupled element pair from the 8 vertex coordinates of each (ve, vf:
  ! x1,y1,z1, x2, ...): kind 0 when they share a vertex unshifted, 1 when they share one only
  ! after the shift shift*period (shift in {-1,0,1}^3, nonzero only on axes with a finite
  ! period, the first such shift found), -1 when they share none in any image. tol: the
  ! coincidence tolerance (DEpsPrec in the census). Used for nbad only: kind 0 does not mean
  ! that the pair has no periodic coupling (it may share a periodic image as well).
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_PAIR_SHIFT(ve, vf, period, tol, shift, kind)
    REAL*8, INTENT(IN) :: ve(24), vf(24), period(3), tol
    INTEGER, INTENT(OUT) :: shift(3), kind
    INTEGER :: sx, sy, sz
    kind = 0
    shift = 0
    IF (share((/0, 0, 0/))) RETURN
    DO sx = -1, 1
    DO sy = -1, 1
    DO sz = -1, 1
      shift = (/sx, sy, sz/)
      IF (ALL(shift.EQ.0)) CYCLE
      IF (ANY(shift.NE.0 .AND. period.GE.1d8)) CYCLE
      IF (share(shift)) THEN
        kind = 1
        RETURN
      END IF
    END DO
    END DO
    END DO
    shift = 0
    kind = -1
  CONTAINS
    LOGICAL FUNCTION share(s)
      INTEGER, INTENT(IN) :: s(3)
      INTEGER :: i, j
      REAL*8 :: q(3)
      share = .FALSE.
      DO i = 1, 8
        q = ve(3*(i-1)+1:3*i) + DBLE(s)*MERGE(period, 0d0, s.NE.0)
        DO j = 1, 8
          IF (NORM2(q - vf(3*(j-1)+1:3*j)).LT.tol) THEN
            share = .TRUE.
            RETURN
          END IF
        END DO
      END DO
    END FUNCTION share
  END SUBROUTINE PCO_PAIR_SHIFT

  !-----------------------------------------------------------------------
  ! The master holds the coordinates of the Q2 dofs of level lev (the vertices of level
  ! lev+1, shared storage of the finest master level) - not every application refines the
  ! master mesh beyond its MedLev (General_init of the LES driver does not).
  !-----------------------------------------------------------------------
  LOGICAL FUNCTION PCO_MASTER_HAS_DOF_COORDINATES(lev) RESULT(ok)
    INTEGER, INTENT(IN) :: lev
    INTEGER :: ndof
    ok = ALLOCATED(mg_mesh%level)
    IF (ok) ok = (lev.GE.1 .AND. SIZE(mg_mesh%level).GE.lev+1)
    IF (ok) ok = ASSOCIATED(mg_mesh%level(lev+1)%dcorvg)
    IF (ok) THEN
      ndof = mg_mesh%level(lev)%nvt + mg_mesh%level(lev)%net + mg_mesh%level(lev)%nat + &
             mg_mesh%level(lev)%nel
      ok = (SIZE(mg_mesh%level(lev+1)%dcorvg,2).GE.ndof .AND. ndof.GT.0)
    END IF
  END FUNCTION PCO_MASTER_HAS_DOF_COORDINATES

  !-----------------------------------------------------------------------
  ! Periodic classes of the Q2 dofs of level lev on the master (see module header).
  ! axis_active / period: the axes found active by the census and their periods.
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_BUILD_CLASSES(lev, axis_active, period, cls)
    INTEGER, INTENT(IN) :: lev
    LOGICAL, INTENT(IN) :: axis_active(3)
    REAL*8, INTENT(IN) :: period(3)
    TYPE(tPcoClasses), INTENT(INOUT) :: cls

    TYPE(tPointHash) :: hash
    INTEGER :: ndof, nvt, net, nat, i, j, k, sx, sy, sz, nhit, hits(2), a, m, sz_k
    INTEGER, ALLOCATABLE :: parent(:), npartner(:), cnt(:), root2cls(:)
    REAL*8, ALLOCATABLE :: x(:,:)
    REAL*8 :: q(3), span(3), tol
    INTEGER :: shift(3)

    nvt = mg_mesh%level(lev)%nvt
    net = mg_mesh%level(lev)%net
    nat = mg_mesh%level(lev)%nat
    ndof = nvt + net + nat + mg_mesh%level(lev)%nel
    tol = DEpsPrec
    ALLOCATE(x(3,ndof))
    x = mg_mesh%level(lev+1)%dcorvg(:,1:ndof)

    CALL hash_build(hash, x, ndof, tol)

    ALLOCATE(parent(ndof), npartner(ndof))
    DO i = 1, ndof
      parent(i) = i
    END DO
    npartner = 0
    cls%nambiguous = 0
    cls%nincomplete = 0
    cls%nmixed = 0
    cls%first_bad = 0

    DO i = 1, ndof
      DO sx = -1, 1
      DO sy = -1, 1
      DO sz = -1, 1
        shift = (/sx, sy, sz/)
        IF (ALL(shift.EQ.0)) CYCLE
        IF (ANY(shift.NE.0 .AND. .NOT.axis_active)) CYCLE
        q = x(:,i) + DBLE(shift)*MERGE(period, 0d0, axis_active)
        CALL hash_query(hash, x, q, tol, hits, nhit)
        IF (nhit.GT.1) THEN
          cls%nambiguous = cls%nambiguous + 1
          IF (cls%first_bad.EQ.0) cls%first_bad = i
        END IF
        IF (nhit.GE.1) THEN
          npartner(i) = npartner(i) + 1
          CALL uf_union(parent, i, hits(1))
        END IF
      END DO
      END DO
      END DO
    END DO

    ! classes in the order of their smallest dof
    ALLOCATE(root2cls(ndof), cnt(ndof))
    root2cls = 0
    IF (ALLOCATED(cls%cls_of)) DEALLOCATE(cls%cls_of)
    ALLOCATE(cls%cls_of(ndof))
    cls%ncls = 0
    DO i = 1, ndof
      j = uf_find(parent, i)
      IF (root2cls(j).EQ.0) THEN
        cls%ncls = cls%ncls + 1
        root2cls(j) = cls%ncls
      END IF
      cls%cls_of(i) = root2cls(j)
    END DO
    IF (ALLOCATED(cls%ptr)) DEALLOCATE(cls%ptr)
    IF (ALLOCATED(cls%mem)) DEALLOCATE(cls%mem)
    ALLOCATE(cls%ptr(cls%ncls+1), cls%mem(ndof))
    cnt = 0
    DO i = 1, ndof
      cnt(cls%cls_of(i)) = cnt(cls%cls_of(i)) + 1
    END DO
    cls%ptr(1) = 1
    DO k = 1, cls%ncls
      cls%ptr(k+1) = cls%ptr(k) + cnt(k)
    END DO
    cnt(1:cls%ncls) = cls%ptr(1:cls%ncls)
    DO i = 1, ndof
      k = cls%cls_of(i)
      cls%mem(cnt(k)) = i
      cnt(k) = cnt(k) + 1
    END DO

    ! completeness and type consistency
    cls%nmulti = 0
    cls%maxsize = 1
    cls%axis_found = .FALSE.
    DO k = 1, cls%ncls
      sz_k = cls%ptr(k+1) - cls%ptr(k)
      IF (sz_k.EQ.1) CYCLE
      cls%nmulti = cls%nmulti + 1
      cls%maxsize = MAX(cls%maxsize, sz_k)
      m = 0
      DO a = 1, 3
        IF (.NOT.axis_active(a)) CYCLE
        span(a) = MAXVAL(x(a, cls%mem(cls%ptr(k):cls%ptr(k+1)-1))) - &
                  MINVAL(x(a, cls%mem(cls%ptr(k):cls%ptr(k+1)-1)))
        IF (span(a).GT.0.5d0*period(a)) THEN
          m = m + 1
          cls%axis_found(a) = .TRUE.
        END IF
      END DO
      IF (sz_k.NE.2**m .OR. ANY(npartner(cls%mem(cls%ptr(k):cls%ptr(k+1)-1)).NE.2**m-1)) THEN
        cls%nincomplete = cls%nincomplete + 1
        IF (cls%first_bad.EQ.0) cls%first_bad = cls%mem(cls%ptr(k))
      END IF
      IF (dof_kind(cls%mem(cls%ptr(k)), nvt, net, nat).NE. &
          dof_kind(cls%mem(cls%ptr(k+1)-1), nvt, net, nat)) THEN
        cls%nmixed = cls%nmixed + 1
        IF (cls%first_bad.EQ.0) cls%first_bad = cls%mem(cls%ptr(k))
      END IF
    END DO
    ! members of a class are ascending, so kinds of the first and the last member bound
    ! the kinds of all (dof_kind is monotone in the dof number)

    cls%ndof = ndof
    cls%nvt = nvt
    DEALLOCATE(x, parent, npartner, cnt, root2cls)
    CALL hash_free(hash)

  END SUBROUTINE PCO_BUILD_CLASSES

  INTEGER FUNCTION dof_kind(i, nvt, net, nat)
    INTEGER, INTENT(IN) :: i, nvt, net, nat
    IF (i.LE.nvt) THEN
      dof_kind = 1
    ELSE IF (i.LE.nvt+net) THEN
      dof_kind = 2
    ELSE IF (i.LE.nvt+net+nat) THEN
      dof_kind = 3
    ELSE
      dof_kind = 4
    END IF
  END FUNCTION dof_kind

  !-----------------------------------------------------------------------
  ! Number of classes whose members carry different velocity constraint flags (Dirichlet
  ! per component, the test of Get_CMat: flag == 1). On the workers the flags of coupled
  ! dofs are the same; a class with different flags on the master is a periodic pair with
  ! different physical conditions on its two sides.
  !-----------------------------------------------------------------------
  INTEGER FUNCTION PCO_CLASS_MASK_CONFLICTS(cls, ku, kv, kw, first) RESULT(nconf)
    TYPE(tPcoClasses), INTENT(IN) :: cls
    INTEGER, INTENT(IN) :: ku(:), kv(:), kw(:)
    INTEGER, INTENT(OUT), OPTIONAL :: first      ! a dof of the first conflicting class
    INTEGER :: k, j, r
    nconf = 0
    IF (PRESENT(first)) first = 0
    DO k = 1, cls%ncls
      IF (cls%ptr(k+1)-cls%ptr(k).EQ.1) CYCLE
      r = cls%mem(cls%ptr(k))
      DO j = cls%ptr(k)+1, cls%ptr(k+1)-1
        IF (PCO_MASKS_DIFFER(ku(cls%mem(j)), kv(cls%mem(j)), kw(cls%mem(j)), ku(r), kv(r), kw(r))) THEN
          nconf = nconf + 1
          IF (PRESENT(first)) THEN
            IF (first.EQ.0) first = r
          END IF
          EXIT
        END IF
      END DO
    END DO
  END FUNCTION PCO_CLASS_MASK_CONFLICTS

  !-----------------------------------------------------------------------
  ! .TRUE. when two dofs carry different Dirichlet masks (flag == 1 per component, the test
  ! of Get_CMat). Each component compared on its own: .NEQV. binds weaker than .OR.
  !-----------------------------------------------------------------------
  LOGICAL FUNCTION PCO_MASKS_DIFFER(u1, v1, w1, u2, v2, w2)
    INTEGER, INTENT(IN) :: u1, v1, w1, u2, v2, w2
    LOGICAL :: m1(3), m2(3)
    m1 = (/u1.EQ.1, v1.EQ.1, w1.EQ.1/)
    m2 = (/u2.EQ.1, v2.EQ.1, w2.EQ.1/)
    PCO_MASKS_DIFFER = ANY(m1 .NEQV. m2)
  END FUNCTION PCO_MASKS_DIFFER

  !-----------------------------------------------------------------------
  ! Element graph of level lev through vertex classes: f is a neighbour of e (f /= e) when a
  ! vertex of e and a vertex of f are in one class. eptr/enbr: CSR, neighbours ascending.
  ! The vertex classes are the classes of the first nvt dofs (a class never mixes kinds).
  !-----------------------------------------------------------------------
  SUBROUTINE PCO_ELEMENT_GRAPH(lev, cls, eptr, enbr)
    INTEGER, INTENT(IN) :: lev
    TYPE(tPcoClasses), INTENT(IN) :: cls
    INTEGER, ALLOCATABLE, INTENT(OUT) :: eptr(:), enbr(:)
    INTEGER :: nel, nvt, iel, j, v, w, k, f, n, ipos, mw
    INTEGER, ALLOCATABLE :: vptr(:), vel(:), fill(:), stamp(:), buf(:)

    nel = mg_mesh%level(lev)%nel
    nvt = mg_mesh%level(lev)%nvt
    ALLOCATE(vptr(nvt+1), fill(nvt+1), vel(8*nel), stamp(nel), buf(nel))
    fill = 0
    DO iel = 1, nel
      DO j = 1, 8
        v = mg_mesh%level(lev)%kvert(j,iel)
        fill(v+1) = fill(v+1) + 1
      END DO
    END DO
    vptr(1) = 1
    DO v = 1, nvt
      vptr(v+1) = vptr(v) + fill(v+1)
    END DO
    fill(1:nvt) = vptr(1:nvt)
    DO iel = 1, nel
      DO j = 1, 8
        v = mg_mesh%level(lev)%kvert(j,iel)
        vel(fill(v)) = iel
        fill(v) = fill(v) + 1
      END DO
    END DO

    stamp = 0
    ALLOCATE(eptr(nel+1))
    ! two passes: count, then fill
    DO ipos = 1, 2
      eptr(1) = 1
      DO iel = 1, nel
        n = 0
        DO j = 1, 8
          v = mg_mesh%level(lev)%kvert(j,iel)
          k = cls%cls_of(v)
          DO mw = cls%ptr(k), cls%ptr(k+1)-1
            w = cls%mem(mw)
            DO f = vptr(w), vptr(w+1)-1
              IF (vel(f).EQ.iel) CYCLE
              IF (stamp(vel(f)).EQ.iel + (ipos-1)*nel) CYCLE
              stamp(vel(f)) = iel + (ipos-1)*nel
              n = n + 1
              buf(n) = vel(f)
            END DO
          END DO
        END DO
        IF (ipos.EQ.2) THEN
          CALL isort(buf, n)
          enbr(eptr(iel):eptr(iel)+n-1) = buf(1:n)
        END IF
        eptr(iel+1) = eptr(iel) + n
      END DO
      IF (ipos.EQ.1) ALLOCATE(enbr(MAX(eptr(nel+1)-1,1)))
    END DO
    DEALLOCATE(vptr, fill, vel, stamp, buf)

  END SUBROUTINE PCO_ELEMENT_GRAPH

  SUBROUTINE isort(a, n)
    INTEGER, INTENT(INOUT) :: a(:)
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, j, t
    DO i = 2, n
      t = a(i)
      j = i - 1
      DO WHILE (j.GE.1)
        IF (a(j).LE.t) EXIT
        a(j+1) = a(j)
        j = j - 1
      END DO
      a(j+1) = t
    END DO
  END SUBROUTINE isort

  !-----------------------------------------------------------------------
  ! union-find
  !-----------------------------------------------------------------------
  INTEGER FUNCTION uf_find(parent, i) RESULT(r)
    INTEGER, INTENT(INOUT) :: parent(:)
    INTEGER, INTENT(IN) :: i
    INTEGER :: j, t
    r = i
    DO WHILE (parent(r).NE.r)
      r = parent(r)
    END DO
    j = i
    DO WHILE (parent(j).NE.r)
      t = parent(j)
      parent(j) = r
      j = t
    END DO
  END FUNCTION uf_find

  SUBROUTINE uf_union(parent, i, j)
    INTEGER, INTENT(INOUT) :: parent(:)
    INTEGER, INTENT(IN) :: i, j
    INTEGER :: ri, rj
    ri = uf_find(parent, i)
    rj = uf_find(parent, j)
    IF (ri.EQ.rj) RETURN
    IF (ri.LT.rj) THEN
      parent(rj) = ri
    ELSE
      parent(ri) = rj
    END IF
  END SUBROUTINE uf_union

  !-----------------------------------------------------------------------
  ! coordinate hash: cells of size h >= 4 tol, keys sorted; a query visits the 27 cells
  ! around the probe and returns the points within tol (at most two reported, nhit counts all)
  !-----------------------------------------------------------------------
  SUBROUTINE hash_build(hs, x, n, tol)
    TYPE(tPointHash), INTENT(INOUT) :: hs
    INTEGER, INTENT(IN) :: n
    REAL*8, INTENT(IN) :: x(3,n), tol
    INTEGER :: i
    REAL*8 :: ext
    hs%n = n
    hs%x0 = MINVAL(x, DIM=2)
    ext = MAXVAL(MAXVAL(x, DIM=2) - hs%x0)
    hs%h = MAX(4d0*tol, 2d0*ext/DBLE(HB/2))
    hs%x0 = hs%x0 - 2d0*hs%h
    ALLOCATE(hs%key(n), hs%idx(n))
    DO i = 1, n
      hs%key(i) = cell_key(hs, x(:,i))
      hs%idx(i) = i
    END DO
    CALL heapsort_key(hs%key, hs%idx, n)
  END SUBROUTINE hash_build

  SUBROUTINE hash_free(hs)
    TYPE(tPointHash), INTENT(INOUT) :: hs
    IF (ALLOCATED(hs%key)) DEALLOCATE(hs%key)
    IF (ALLOCATED(hs%idx)) DEALLOCATE(hs%idx)
    hs%n = 0
  END SUBROUTINE hash_free

  INTEGER(8) FUNCTION cell_key(hs, p) RESULT(key)
    TYPE(tPointHash), INTENT(IN) :: hs
    REAL*8, INTENT(IN) :: p(3)
    key = index_key(FLOOR((p - hs%x0)/hs%h, KIND=8))
  END FUNCTION cell_key

  INTEGER(8) FUNCTION index_key(ic) RESULT(key)
    INTEGER(8), INTENT(IN) :: ic(3)
    IF (ANY(ic.LT.0_8) .OR. ANY(ic.GE.HB)) THEN
      key = -1_8
    ELSE
      key = (ic(1)*HB + ic(2))*HB + ic(3)
    END IF
  END FUNCTION index_key

  SUBROUTINE hash_query(hs, x, q, tol, hits, nhit)
    TYPE(tPointHash), INTENT(IN) :: hs
    REAL*8, INTENT(IN) :: x(:,:), q(3), tol
    INTEGER, INTENT(OUT) :: hits(2), nhit
    INTEGER :: dx, dy, dz, lo, hi, m
    INTEGER(8) :: key, ic0(3)
    nhit = 0
    hits = 0
    ! the 27 cells around the cell of q, by integer offsets (distinct cells, no rounding)
    ic0 = FLOOR((q - hs%x0)/hs%h, KIND=8)
    DO dx = -1, 1
    DO dy = -1, 1
    DO dz = -1, 1
      key = index_key(ic0 + INT((/dx, dy, dz/), KIND=8))
      IF (key.LT.0_8) CYCLE
      CALL key_range(hs%key, hs%n, key, lo, hi)
      DO m = lo, hi
        IF (NORM2(x(:,hs%idx(m)) - q).LT.tol) THEN
          nhit = nhit + 1
          IF (nhit.LE.2) hits(nhit) = hs%idx(m)
        END IF
      END DO
    END DO
    END DO
    END DO
  END SUBROUTINE hash_query

  ! positions lo..hi of key in the sorted array (hi < lo: absent)
  SUBROUTINE key_range(keys, n, key, lo, hi)
    INTEGER(8), INTENT(IN) :: keys(:), key
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(OUT) :: lo, hi
    INTEGER :: l, r, mid
    l = 1
    r = n + 1
    DO WHILE (l.LT.r)
      mid = (l + r)/2
      IF (keys(mid).LT.key) THEN
        l = mid + 1
      ELSE
        r = mid
      END IF
    END DO
    lo = l
    hi = l - 1
    DO WHILE (hi+1.LE.n)
      IF (keys(hi+1).NE.key) EXIT
      hi = hi + 1
    END DO
  END SUBROUTINE key_range

  SUBROUTINE heapsort_key(key, idx, n)
    INTEGER(8), INTENT(INOUT) :: key(:)
    INTEGER, INTENT(INOUT) :: idx(:)
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, iend
    DO i = n/2, 1, -1
      CALL sift(i, n)
    END DO
    DO iend = n, 2, -1
      CALL swap(1, iend)
      CALL sift(1, iend-1)
    END DO
  CONTAINS
    SUBROUTINE sift(start, last)
      INTEGER, INTENT(IN) :: start, last
      INTEGER :: root, child
      root = start
      DO
        child = 2*root
        IF (child.GT.last) EXIT
        IF (child+1.LE.last) THEN
          IF (key(child+1).GT.key(child)) child = child + 1
        END IF
        IF (key(root).GE.key(child)) EXIT
        CALL swap(root, child)
        root = child
      END DO
    END SUBROUTINE sift
    SUBROUTINE swap(a, b)
      INTEGER, INTENT(IN) :: a, b
      INTEGER(8) :: tk
      INTEGER :: ti
      tk = key(a); key(a) = key(b); key(b) = tk
      ti = idx(a); idx(a) = idx(b); idx(b) = ti
    END SUBROUTINE swap
  END SUBROUTINE heapsort_key

END MODULE pco_periodic_topology
