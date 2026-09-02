!=========================================================================
! test_chi_markers - Phase 3 "two-partition cut-cell marker test"
! (chimera-integration-design.md v3, sections 3 and 9).
!
! A 2x1x1 box mesh (elements A: x in [0,1], B: x in [1,2]) is classified
! (a) as one mesh and (b) as two single-element partitions whose shared
! interface nodes (x = 1: 4 vertices, 4 edge midpoints, 1 face centre)
! are merged with the design's MAX rule (CHI_MERGE_MARKER, the semantics
! of the production E013Max_SUPER two-pass synchronisation).  The merged
! partition result must equal the single-mesh classification, for
!   case 1: a sphere inside A only  -> A is a cut cell, the interface
!           nodes are FRINGE seen from A and FREE seen from B (the
!           v2-review counterexample: a signed encoding would have lost
!           the fringe under MAX);
!   case 2: a sphere crossing the interface -> the face centre is a HOLE
!           node, both elements are cut cells, interface corners fringe.
! Also checked: body ids accompany every non-free marker; free nodes
! carry id 0.  Meshes are written as .tri files and loaded through the
! legacy mesh adapter path (readTriCoarse/genMeshStructures), so the Q2
! numbering is the production one.
!=========================================================================
PROGRAM test_chi_markers

  USE types, ONLY: tMultiMesh
  USE mesh_structures, ONLY: readTriCoarse, genMeshStructures, &
    getNumberOfEdgesOnVerts
  USE CHI_MARKERS, ONLY: tChiBody, CHI_BODY_SPHERE, CHI_MARK_FREE, &
    CHI_MARK_FRINGE, CHI_MARK_HOLE, CHI_CLASSIFY_MARKERS, CHI_MERGE_MARKER
  USE CHI_FEM_EVAL, ONLY: CHI_Q2_DOFMAP, CHI_Q2_REFNODES
  USE CHI_GEOMETRY, ONLY: CHI_Q1_MAP

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: ierr, nfail, icase
  TYPE(tChiBody) :: body(1)

  nfail = 0
  CALL MPI_INIT(ierr)

  CALL write_box('tcm_global.tri', 0d0, 2d0, 2)
  CALL write_box('tcm_partA.tri', 0d0, 1d0, 1)
  CALL write_box('tcm_partB.tri', 1d0, 2d0, 1)

  DO icase = 1, 2
    body(1)%shape = CHI_BODY_SPHERE
    IF (icase .EQ. 1) THEN
      body(1)%center = (/ 0.5d0, 0.5d0, 0.5d0 /)
      body(1)%radius = 0.35d0
    ELSE
      body(1)%center = (/ 0.75d0, 0.5d0, 0.5d0 /)
      body(1)%radius = 0.3d0
    END IF
    CALL run_case(icase, body, nfail)
  END DO

  CALL MPI_FINALIZE(ierr)
  IF (nfail .GT. 0) THEN
    WRITE(*,*) 'test_chi_markers: ', nfail, ' failure(s)'
    STOP 1
  END IF
  WRITE(*,*) 'test_chi_markers: PASS'

CONTAINS

  SUBROUTINE write_box(fname, x0, x1, nx)
    CHARACTER(*), INTENT(IN) :: fname
    REAL*8, INTENT(IN) :: x0, x1
    INTEGER, INTENT(IN) :: nx
    INTEGER :: iu, ix, iy, iz, i
    iu = 91
    OPEN(iu, FILE=fname, STATUS='REPLACE')
    WRITE(iu,'(A)') 'test_chi_markers box'
    WRITE(iu,'(A)') 'Parametrisierung PARXC, PARYC, TMAXC'
    WRITE(iu,'(I0,1X,I0,A)') nx, (nx+1)*4, ' 1 8 12 6     NEL NVT NBCT NVE NEE NAE'
    WRITE(iu,'(A)') 'DCORVG'
    DO iz = 0, 1
      DO iy = 0, 1
        DO ix = 0, nx
          WRITE(iu,'(3ES24.16)') x0 + (x1-x0)*ix/nx, DBLE(iy), DBLE(iz)
        END DO
      END DO
    END DO
    WRITE(iu,'(A)') 'KVERT'
    DO ix = 0, nx-1
      WRITE(iu,'(8(I0,1X))') vid(nx,ix,0,0), vid(nx,ix+1,0,0), vid(nx,ix+1,1,0), vid(nx,ix,1,0), &
                             vid(nx,ix,0,1), vid(nx,ix+1,0,1), vid(nx,ix+1,1,1), vid(nx,ix,1,1)
    END DO
    WRITE(iu,'(A)') 'KNPR'
    DO i = 1, (nx+1)*4
      WRITE(iu,'(I0)') 1
    END DO
    CLOSE(iu)
  END SUBROUTINE write_box

  INTEGER FUNCTION vid(nx, ix, iy, iz)
    INTEGER, INTENT(IN) :: nx, ix, iy, iz
    vid = 1 + ix + (nx+1)*(iy + 2*iz)
  END FUNCTION vid

  SUBROUTINE load(fname, m)
    CHARACTER(*), INTENT(IN) :: fname
    TYPE(tMultiMesh), INTENT(INOUT) :: m
    INTEGER :: noe
    ALLOCATE(m%level(1))
    m%nlmin = 1; m%nlmax = 1; m%maxlevel = 1
    m%level(1)%nel = 0
    CALL readTriCoarse(fname, m)
    CALL getNumberOfEdgesOnVerts(m%level(1), noe)
    CALL genMeshStructures(m, .FALSE., 1, noe)
  END SUBROUTINE load

  SUBROUTINE q2coords(m, q2c, ndof)
    TYPE(tMultiMesh), INTENT(IN) :: m
    REAL*8, ALLOCATABLE, INTENT(OUT) :: q2c(:,:)
    INTEGER, INTENT(OUT) :: ndof
    REAL*8 :: refxi(3,27), nodes(3,8), x(3), jac(3,3), detj
    INTEGER :: e, i, idx(27)
    ndof = m%level(1)%nvt + m%level(1)%net + m%level(1)%nat + m%level(1)%nel
    ALLOCATE(q2c(3,ndof))
    CALL CHI_Q2_REFNODES(refxi)
    DO e = 1, m%level(1)%nel
      DO i = 1, 8
        nodes(:,i) = m%level(1)%dcorvg(:,m%level(1)%kvert(i,e))
      END DO
      CALL CHI_Q2_DOFMAP(e, m%level(1)%kvert, m%level(1)%kedge, m%level(1)%karea, &
        m%level(1)%nvt, m%level(1)%net, m%level(1)%nat, idx)
      DO i = 1, 27
        CALL CHI_Q1_MAP(nodes, refxi(:,i), x, jac, detj)
        q2c(:,idx(i)) = x
      END DO
    END DO
  END SUBROUTINE q2coords

  SUBROUTINE classify(m, body, kind, pid, q2c, ndof)
    TYPE(tMultiMesh), INTENT(IN) :: m
    TYPE(tChiBody), INTENT(IN) :: body(1)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: kind(:), pid(:)
    REAL*8, ALLOCATABLE, INTENT(OUT) :: q2c(:,:)
    INTEGER, INTENT(OUT) :: ndof
    CALL q2coords(m, q2c, ndof)
    ALLOCATE(kind(ndof), pid(ndof))
    CALL CHI_CLASSIFY_MARKERS(m%level(1)%nel, m%level(1)%nvt, m%level(1)%net, &
      m%level(1)%nat, m%level(1)%kvert, m%level(1)%kedge, m%level(1)%karea, &
      q2c, 1, body, kind, pid)
  END SUBROUTINE classify

  ! index of the dof of mesh (q2c,ndof) at coordinate x (exact match)
  INTEGER FUNCTION find_dof(q2c, ndof, x)
    REAL*8, INTENT(IN) :: q2c(:,:), x(3)
    INTEGER, INTENT(IN) :: ndof
    INTEGER :: i
    find_dof = 0
    DO i = 1, ndof
      IF (MAXVAL(ABS(q2c(:,i) - x)) .LT. 1d-12) THEN
        find_dof = i
        RETURN
      END IF
    END DO
  END FUNCTION find_dof

  SUBROUTINE run_case(icase, body, cnt)
    INTEGER, INTENT(IN) :: icase
    TYPE(tChiBody), INTENT(IN) :: body(1)
    INTEGER, INTENT(INOUT) :: cnt
    TYPE(tMultiMesh) :: mg, ma, mb
    INTEGER, ALLOCATABLE :: kg(:), pg(:), ka(:), pa(:), kb(:), pb(:)
    REAL*8, ALLOCATABLE :: cg(:,:), ca(:,:), cb(:,:)
    INTEGER :: ng, na, nb, i, ig, ib, nshared, nh, nfr, nmism, kmerged, pmerged
    LOGICAL :: is_iface

    CALL load('tcm_global.tri', mg)
    CALL load('tcm_partA.tri', ma)
    CALL load('tcm_partB.tri', mb)
    CALL classify(mg, body, kg, pg, cg, ng)
    CALL classify(ma, body, ka, pa, ca, na)
    CALL classify(mb, body, kb, pb, cb, nb)

    ! sanity of the single-mesh classification
    nh = COUNT(kg .EQ. CHI_MARK_HOLE)
    nfr = COUNT(kg .EQ. CHI_MARK_FRINGE)
    WRITE(*,'(A,I0,A,I0,A,I0)') ' case ', icase, ': global holes = ', nh, ', fringe = ', nfr
    IF (icase .EQ. 1) THEN
      ! only A's centre node is inside; all other 26 nodes of A are fringe
      IF (nh .NE. 1 .OR. nfr .NE. 26) THEN
        WRITE(*,*) 'FAIL: case 1 global counts'
        cnt = cnt + 1
      END IF
      ! the interface face centre must be fringe in the global view
      ig = find_dof(cg, ng, (/1d0, 0.5d0, 0.5d0/))
      IF (kg(ig) .NE. CHI_MARK_FRINGE) THEN
        WRITE(*,*) 'FAIL: case 1 interface face centre not fringe'
        cnt = cnt + 1
      END IF
    ELSE
      ig = find_dof(cg, ng, (/1d0, 0.5d0, 0.5d0/))
      IF (kg(ig) .NE. CHI_MARK_HOLE) THEN
        WRITE(*,*) 'FAIL: case 2 interface face centre not hole'
        cnt = cnt + 1
      END IF
      IF (nh .LT. 2) THEN
        WRITE(*,*) 'FAIL: case 2 expects at least two hole nodes'
        cnt = cnt + 1
      END IF
    END IF
    DO i = 1, ng
      IF ((kg(i) .NE. CHI_MARK_FREE) .NEQV. (pg(i) .EQ. 1)) THEN
        WRITE(*,*) 'FAIL: body id inconsistent with kind at dof', i
        cnt = cnt + 1
        EXIT
      END IF
    END DO

    ! partition view: every global dof is found in A and/or B; merge
    nshared = 0
    nmism = 0
    DO i = 1, ng
      ig = find_dof(ca, na, cg(:,i))
      ib = find_dof(cb, nb, cg(:,i))
      IF (ig .EQ. 0 .AND. ib .EQ. 0) THEN
        WRITE(*,*) 'FAIL: global dof missing from both partitions', i
        cnt = cnt + 1
        CYCLE
      END IF
      is_iface = (ig .NE. 0 .AND. ib .NE. 0)
      IF (is_iface) THEN
        nshared = nshared + 1
        kmerged = ka(ig); pmerged = pa(ig)
        CALL CHI_MERGE_MARKER(kmerged, pmerged, kb(ib), pb(ib))
      ELSE IF (ig .NE. 0) THEN
        kmerged = ka(ig); pmerged = pa(ig)
      ELSE
        kmerged = kb(ib); pmerged = pb(ib)
      END IF
      IF (kmerged .NE. kg(i) .OR. pmerged .NE. pg(i)) THEN
        nmism = nmism + 1
        IF (nmism .LE. 3) WRITE(*,'(A,3F6.2,A,2I3,A,2I3)') '   mismatch at x=', &
          cg(:,i), ' merged (kind,pid)=', kmerged, pmerged, ' global=', kg(i), pg(i)
      END IF
    END DO
    IF (nshared .NE. 9) THEN
      WRITE(*,*) 'FAIL: expected 9 shared interface dofs, got', nshared
      cnt = cnt + 1
    END IF
    IF (nmism .GT. 0) THEN
      WRITE(*,*) 'FAIL: merged partition markers differ from the single-mesh result:', nmism
      cnt = cnt + 1
    END IF
    ! in case 1 the interface nodes are FREE in B's own view (the trap)
    IF (icase .EQ. 1) THEN
      ib = find_dof(cb, nb, (/1d0, 0.5d0, 0.5d0/))
      ig = find_dof(ca, na, (/1d0, 0.5d0, 0.5d0/))
      IF (kb(ib) .NE. CHI_MARK_FREE .OR. ka(ig) .NE. CHI_MARK_FRINGE) THEN
        WRITE(*,*) 'FAIL: case 1 partition views of the interface face centre', kb(ib), ka(ig)
        cnt = cnt + 1
      END IF
    END IF
    WRITE(*,'(A,I0,A)') ' case ', icase, ': merged partition markers == single-mesh markers'
  END SUBROUTINE run_case

END PROGRAM test_chi_markers
