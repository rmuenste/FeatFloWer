!=========================================================================
! EL_PERIODIC - periodic images of the particle-mesh coupling
!
! The one place that knows about periodic images (SimPar@ELPeriodicKernel = wrap, EL_CONFIG).
! On a periodic axis the fluid domain is one cell of a lattice; a particle at x stands for the
! whole row of images x + n L (n integer, L the period of the axis). The kernel of the particle
! is the sum of the kernels of its images, so that
!   * a particle within one kernel radius of a periodic plane deposits and samples on both
!     sides of the plane,
!   * its weights still sum to one over the whole support (the normalisation is the sum of the
!     kernel integrals of all images on all ranks, as in the interior),
!   * the fields do not change when the particle is carried across the plane and its position
!     jumps by L (translation invariance).
!
! Conventions
!   * Lattice: per axis a flag and a length (EL_SET_PERIODIC_LATTICE). Only the LENGTH enters,
!     never an origin: an image is position + shift with shift(a) = n(a) L(a). The routines are
!     therefore independent of where the particle engine keeps a wrapped particle ([0,L), (0,L],
!     or a position slightly outside the box during a step); a particle exactly on a periodic
!     plane has images on both sides like any other particle near it.
!   * A record (EL_HALO) carries the particle position as its owner knows it, once per rank.
!     Images are enumerated locally by the rank that integrates: EL_IMAGE_SHIFTS returns every
!     lattice shift whose kernel sphere touches the box of that rank, the unshifted particle
!     included when it does. A rank that touches both periodic planes of an axis (few ranks on
!     that axis, or one) simply finds two or three images of the same record; a periodic edge
!     or corner gives up to 4 or 8.
!   * The halo exchange needs two predicates, both periodic distances: rank boxes are neighbours
!     when the gap between one box and the nearest image of the other is within the kernel
!     radius (EL_PERIODIC_BOXES_WITHIN_DELTA), and a particle is sent to a neighbour when one
!     of its images touches that box (EL_PERIODIC_SPHERE_REACHES_BOX).
!   * Rank boxes are axis-aligned bounding boxes inside one lattice cell (extent <= L on a
!     periodic axis); the kernel radius must be smaller than every period
!     (EL_PERIODIC_KERNEL_FITS, checked by the halo before the first exchange). Together they
!     bound the images of one particle on one rank by 3 per periodic axis, EL_MAX_IMAGES = 27.
!     A kernel wider than a subdomain (EL_HALO: warned_wide_kernel) or wider than half a period
!     (a particle overlapping its own image) is covered.
!
! Stateless apart from the lattice, which is set once at start-up; no MPI. With the lattice
! inactive (the default, SimPar@ELPeriodicKernel = off) no caller reaches the image routines:
! EL_HALO and EL_TRANSFER branch on EL_PERIODIC_ACTIVE and run the code of before.
!
! Extension points: a sheared lattice (Lees-Edwards) would add a position- and time-dependent
! offset to the shift of the sheared axis in EL_IMAGE_SHIFTS and a velocity offset to the
! record of an image (EL_TRANSFER: EL_SAMPLE_RECORD / EL_DEPOSIT_RECORD are the only callers).
!=========================================================================
MODULE EL_PERIODIC

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: EL_MAX_IMAGES, EL_SET_PERIODIC_LATTICE, EL_CLEAR_PERIODIC_LATTICE, &
            EL_PERIODIC_ACTIVE, EL_PERIODIC_AXES, EL_PERIODIC_LENGTHS, &
            EL_PERIODIC_KERNEL_FITS, EL_IMAGE_SHIFTS, EL_PERIODIC_MARKED_ON_PLANES, &
            EL_PERIODIC_SPHERE_REACHES_BOX, EL_PERIODIC_BOXES_WITHIN_DELTA

  ! most images of one particle on one rank: 3 per periodic axis (module header)
  INTEGER, PARAMETER :: EL_MAX_IMAGES = 27

  LOGICAL, SAVE :: lattice_axis(3) = .FALSE.
  REAL*8, SAVE :: lattice_length(3) = 0.0d0

CONTAINS

  ! Set the lattice: axis(a) says whether axis a is periodic, length(a) is its period (read
  ! only where axis(a) is set). Called once at start-up, on every rank, with deck values.
  SUBROUTINE EL_SET_PERIODIC_LATTICE(axis, length)

    LOGICAL, INTENT(IN) :: axis(3)
    REAL*8, INTENT(IN) :: length(3)
    INTEGER :: a

    DO a = 1, 3
      IF (axis(a) .AND. .NOT.(length(a).GT.0.0d0)) THEN
        WRITE(*,'(A,I0,A,ES14.6)') 'EL_SET_PERIODIC_LATTICE: period of axis ', a, &
          ' must be positive, got ', length(a)
        ERROR STOP 1
      END IF
    END DO
    lattice_axis = axis
    lattice_length = 0.0d0
    WHERE (axis) lattice_length = length

  END SUBROUTINE EL_SET_PERIODIC_LATTICE

  SUBROUTINE EL_CLEAR_PERIODIC_LATTICE()
    lattice_axis = .FALSE.
    lattice_length = 0.0d0
  END SUBROUTINE EL_CLEAR_PERIODIC_LATTICE

  ! True when the coupling wraps on at least one axis. The same on every rank.
  PURE LOGICAL FUNCTION EL_PERIODIC_ACTIVE()
    EL_PERIODIC_ACTIVE = ANY(lattice_axis)
  END FUNCTION EL_PERIODIC_ACTIVE

  PURE FUNCTION EL_PERIODIC_AXES() RESULT(axis)
    LOGICAL :: axis(3)
    axis = lattice_axis
  END FUNCTION EL_PERIODIC_AXES

  ! Periods of the three axes; 0 on an axis that is not periodic.
  PURE FUNCTION EL_PERIODIC_LENGTHS() RESULT(length)
    REAL*8 :: length(3)
    length = lattice_length
  END FUNCTION EL_PERIODIC_LENGTHS

  ! The image bound of the module header holds for this kernel radius.
  PURE LOGICAL FUNCTION EL_PERIODIC_KERNEL_FITS(radius)
    REAL*8, INTENT(IN) :: radius
    INTEGER :: a
    EL_PERIODIC_KERNEL_FITS = .TRUE.
    DO a = 1, 3
      IF (lattice_axis(a) .AND. radius.GE.lattice_length(a)) EL_PERIODIC_KERNEL_FITS = .FALSE.
    END DO
  END FUNCTION EL_PERIODIC_KERNEL_FITS

  ! Lattice shifts of the images of a particle that reach a box: every shift s = n L with
  ! |gap(position + s, box)| <= radius, where gap is the distance vector from the box to the
  ! point (zero inside). box = (xmin, xmax, ymin, ymax, zmin, zmax). The unshifted particle is
  ! one of them when it reaches the box itself; n = 0 when no image does. Without a lattice the
  ! only candidate is the particle itself.
  SUBROUTINE EL_IMAGE_SHIFTS(position, radius, box, shifts, n)

    REAL*8, INTENT(IN) :: position(3), radius, box(6)
    REAL*8, INTENT(OUT) :: shifts(3,EL_MAX_IMAGES)
    INTEGER, INTENT(OUT) :: n

    INTEGER :: lo(3), hi(3), a, i1, i2, i3
    REAL*8 :: s(3), gap2_3, gap2_23, gap2

    ! per axis, the integers m with box_min - radius <= position + m L <= box_max + radius
    lo = 0
    hi = 0
    DO a = 1, 3
      IF (.NOT.lattice_axis(a)) CYCLE
      lo(a) = CEILING((box(2*a-1) - radius - position(a))/lattice_length(a))
      hi(a) = FLOOR((box(2*a) + radius - position(a))/lattice_length(a))
    END DO

    n = 0
    shifts = 0.0d0
    DO i3 = lo(3), hi(3)
      s(3) = DBLE(i3)*lattice_length(3)
      gap2_3 = axis_gap(3, s(3))**2
      DO i2 = lo(2), hi(2)
        s(2) = DBLE(i2)*lattice_length(2)
        gap2_23 = gap2_3 + axis_gap(2, s(2))**2
        DO i1 = lo(1), hi(1)
          s(1) = DBLE(i1)*lattice_length(1)
          gap2 = gap2_23 + axis_gap(1, s(1))**2
          IF (gap2.GT.radius*radius) CYCLE
          IF (n.GE.EL_MAX_IMAGES) THEN
            ! excluded by EL_PERIODIC_KERNEL_FITS and a rank box inside one lattice cell
            WRITE(*,'(A,I0,A,ES14.6,A,3ES14.6)') 'EL_IMAGE_SHIFTS: more than ', EL_MAX_IMAGES, &
              ' images of one particle reach a rank box; kernel radius ', radius, &
              ' periods ', lattice_length
            ERROR STOP 1
          END IF
          n = n + 1
          shifts(:,n) = s
        END DO
      END DO
    END DO

  CONTAINS

    ! distance along axis a from the box to the image shifted by sa (zero inside the box)
    PURE REAL*8 FUNCTION axis_gap(a, sa)
      INTEGER, INTENT(IN) :: a
      REAL*8, INTENT(IN) :: sa
      axis_gap = MAX(0.0d0, box(2*a-1) - (position(a) + sa), (position(a) + sa) - box(2*a))
    END FUNCTION axis_gap

  END SUBROUTINE EL_IMAGE_SHIFTS

  ! One of the images of the particle reaches the box (periodic form of the sphere-box test
  ! of the halo send list).
  LOGICAL FUNCTION EL_PERIODIC_SPHERE_REACHES_BOX(position, radius, box)

    REAL*8, INTENT(IN) :: position(3), radius, box(6)
    REAL*8 :: shifts(3,EL_MAX_IMAGES)
    INTEGER :: n

    CALL EL_IMAGE_SHIFTS(position, radius, box, shifts, n)
    EL_PERIODIC_SPHERE_REACHES_BOX = (n.GT.0)

  END FUNCTION EL_PERIODIC_SPHERE_REACHES_BOX

  ! The gap between box a and the nearest image of box b is within delta (periodic form of the
  ! neighbour test of the halo graph; symmetric in a and b). Both boxes lie inside one lattice
  ! cell, so the nearest image along an axis is b itself or b shifted by one period either way.
  PURE LOGICAL FUNCTION EL_PERIODIC_BOXES_WITHIN_DELTA(a, b, delta)

    REAL*8, INTENT(IN) :: a(6), b(6), delta
    REAL*8 :: gap(3), s
    INTEGER :: k, m

    DO k = 1, 3
      gap(k) = MAX(0.0d0, MAX(a(2*k-1)-b(2*k), b(2*k-1)-a(2*k)))
      IF (.NOT.lattice_axis(k)) CYCLE
      DO m = -1, 1, 2
        s = DBLE(m)*lattice_length(k)
        gap(k) = MIN(gap(k), MAX(0.0d0, MAX(a(2*k-1)-(b(2*k)+s), (b(2*k-1)+s)-a(2*k))))
      END DO
    END DO
    EL_PERIODIC_BOXES_WITHIN_DELTA = SQRT(SUM(gap*gap)).LE.delta

  END FUNCTION EL_PERIODIC_BOXES_WITHIN_DELTA

  ! Number of marked points on a periodic plane of the domain box (the two faces of the box
  ! normal to each periodic axis; within 1e-8 of the period). Used to verify that a mask of
  ! kernel-cutting boundary faces does not contain a periodic plane (EL_TRANSFER).
  PURE INTEGER FUNCTION EL_PERIODIC_MARKED_ON_PLANES(coor, marked, domain_box) RESULT(count)

    REAL*8, INTENT(IN) :: coor(:,:), domain_box(6)
    LOGICAL, INTENT(IN) :: marked(:)
    INTEGER :: i, a
    REAL*8 :: tol

    count = 0
    DO i = 1, MIN(SIZE(marked), SIZE(coor,2))
      IF (.NOT.marked(i)) CYCLE
      DO a = 1, 3
        IF (.NOT.lattice_axis(a)) CYCLE
        tol = 1.0d-8*lattice_length(a)
        IF (ABS(coor(a,i)-domain_box(2*a-1)).LE.tol .OR. ABS(coor(a,i)-domain_box(2*a)).LE.tol) THEN
          count = count + 1
          EXIT
        END IF
      END DO
    END DO

  END FUNCTION EL_PERIODIC_MARKED_ON_PLANES

END MODULE EL_PERIODIC
