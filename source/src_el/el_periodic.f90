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
!   * Rank boxes are axis-aligned bounding boxes that do not span a period (extent < L on a
!     periodic axis; EL_PERIODIC_BOX_SPANS_AXIS, refused by the halo), the lattice is the
!     fluid domain (EL_PERIODIC_LATTICE_FITS_DOMAIN, refused at start-up otherwise) and the
!     kernel radius is smaller than every period (EL_PERIODIC_KERNEL_FITS, refused by the
!     halo). Then the images of one particle that reach one rank are among the three lattice
!     shifts per periodic axis next to the box, EL_MAX_IMAGES = 27; EL_IMAGE_SHIFTS looks at
!     exactly those, so the bound holds whatever the input. A kernel wider than a subdomain
!     (EL_HALO: warned_wide_kernel) or wider than half a period (a particle overlapping its
!     own image) is covered.
!   * The neighbour test is symmetric bit for bit (EL_PERIODIC_BOXES_WITHIN_DELTA): two ranks
!     evaluate it for the same pair and must agree, or the neighbour graph is inconsistent.
!
! Stateless apart from the lattice, which is set once at start-up; no MPI, no stop: a
! refusal is returned to the caller, which ends the run in the way of its layer. With the
! lattice inactive (the default, SimPar@ELPeriodicKernel = off) no caller reaches the image
! routines: EL_HALO, EL_IMAGE_TRANSFER, EL_TRANSFER and EL_GEOMETRY branch on
! EL_PERIODIC_ACTIVE and run the code of before.
!
! Extension points: a sheared lattice (Lees-Edwards) would add a position- and time-dependent
! offset to the shift of the sheared axis in EL_IMAGE_SHIFTS and a velocity offset to the
! record of an image (EL_IMAGE_TRANSFER: EL_RECORD_IMAGES is the only caller).
!=========================================================================
MODULE EL_PERIODIC

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: EL_MAX_IMAGES, EL_SET_PERIODIC_LATTICE, EL_CLEAR_PERIODIC_LATTICE, &
            EL_PERIODIC_ACTIVE, EL_PERIODIC_AXES, EL_PERIODIC_LENGTHS, &
            EL_PERIODIC_KERNEL_FITS, EL_IMAGE_SHIFTS, EL_PERIODIC_MARKED_ON_PLANES, &
            EL_PERIODIC_SPHERE_REACHES_BOX, EL_PERIODIC_BOXES_WITHIN_DELTA, &
            EL_PERIODIC_LATTICE_FITS_DOMAIN, EL_PERIODIC_BOX_SPANS_AXIS, EL_PERIODIC_TOLERANCE

  ! most images of one particle on one rank: 3 per periodic axis (module header)
  INTEGER, PARAMETER :: EL_MAX_IMAGES = 27
  ! relative tolerance (in periods) of the geometric start-up checks: lattice against domain,
  ! rank box against period, points on a periodic plane
  REAL*8, PARAMETER :: EL_PERIODIC_TOLERANCE = 1.0d-8

  LOGICAL, SAVE :: lattice_axis(3) = .FALSE.
  REAL*8, SAVE :: lattice_length(3) = 0.0d0

CONTAINS

  ! Set the lattice: axis(a) says whether axis a is periodic, length(a) is its period (read
  ! only where axis(a) is set). Called once at start-up, on every rank, with deck values.
  ! ok = .FALSE. and no lattice when the period of a periodic axis is not a positive finite
  ! number; the caller stops the run.
  SUBROUTINE EL_SET_PERIODIC_LATTICE(axis, length, ok)

    LOGICAL, INTENT(IN) :: axis(3)
    REAL*8, INTENT(IN) :: length(3)
    LOGICAL, INTENT(OUT) :: ok
    INTEGER :: a

    ok = .TRUE.
    DO a = 1, 3
      IF (axis(a) .AND. .NOT.(length(a).GT.0.0d0 .AND. length(a).LE.HUGE(1.0d0))) ok = .FALSE.
    END DO
    lattice_axis = .FALSE.
    lattice_length = 0.0d0
    IF (.NOT.ok) RETURN
    lattice_axis = axis
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

  ! The lattice is the fluid domain: on every periodic axis the period equals the extent of
  ! domain_box = (xmin, xmax, ymin, ymax, zmin, zmax) within EL_PERIODIC_TOLERANCE periods.
  ! bad_axis is the first axis that fails (0 when the function is true).
  LOGICAL FUNCTION EL_PERIODIC_LATTICE_FITS_DOMAIN(domain_box, bad_axis) RESULT(fits)
    REAL*8, INTENT(IN) :: domain_box(6)
    INTEGER, INTENT(OUT) :: bad_axis
    INTEGER :: a
    fits = .TRUE.
    bad_axis = 0
    DO a = 1, 3
      IF (.NOT.lattice_axis(a)) CYCLE
      ! written so that a NaN extent fails
      IF (.NOT.(ABS((domain_box(2*a)-domain_box(2*a-1)) - lattice_length(a)).LE. &
                EL_PERIODIC_TOLERANCE*lattice_length(a))) THEN
        fits = .FALSE.
        IF (bad_axis.EQ.0) bad_axis = a
      END IF
    END DO
  END FUNCTION EL_PERIODIC_LATTICE_FITS_DOMAIN

  ! A rank box reaches from one periodic plane of an axis to the other (extent >= period,
  ! within EL_PERIODIC_TOLERANCE periods): the rank would hold both planes of that axis.
  PURE LOGICAL FUNCTION EL_PERIODIC_BOX_SPANS_AXIS(box) RESULT(spans)
    REAL*8, INTENT(IN) :: box(6)
    INTEGER :: a
    spans = .FALSE.
    DO a = 1, 3
      IF (.NOT.lattice_axis(a)) CYCLE
      IF (.NOT.(box(2*a)-box(2*a-1).LT.lattice_length(a)*(1.0d0-EL_PERIODIC_TOLERANCE))) &
        spans = .TRUE.
    END DO
  END FUNCTION EL_PERIODIC_BOX_SPANS_AXIS

  ! Lattice shifts of the images of a particle that reach a box: every shift s = n L with
  ! |gap(position + s, box)| <= radius, where gap is the distance vector from the box to the
  ! point (zero inside). box = (xmin, xmax, ymin, ymax, zmin, zmax). The unshifted particle is
  ! one of them when it reaches the box itself; n = 0 when no image does. Without a lattice the
  ! only candidate is the particle itself.
  !
  ! Per periodic axis the candidates are the image nearest to the centre of the box and its
  ! two lattice neighbours. With a box that does not span the period and a radius below the
  ! period (module header) no other image can reach the box: the nearest image is within L/2
  ! of the centre, the next ones at least 3L/2 away, the reach is below L/2 + L. The nearest
  ! image is found in real arithmetic (ANINT), so a position far outside the cell needs no
  ! integer conversion, and n <= EL_MAX_IMAGES by construction. A position that is not finite,
  ! or so far outside the cell that neighbouring lattice indices are no longer distinct
  ! numbers (2**52 periods), has no image (n = 0; the owner then finds no kernel support and
  ! stops the run, EL_TRANSFER: error 901).
  SUBROUTINE EL_IMAGE_SHIFTS(position, radius, box, shifts, n)

    REAL*8, INTENT(IN) :: position(3), radius, box(6)
    REAL*8, INTENT(OUT) :: shifts(3,EL_MAX_IMAGES)
    INTEGER, INTENT(OUT) :: n

    REAL*8, PARAMETER :: EXACT_INDEX_LIMIT = 4.0d15     ! below 2**52: base + 1 is exact
    INTEGER :: m(3), a, i1, i2, i3
    REAL*8 :: base(3), s(3), gap2_3, gap2_23, gap2

    n = 0
    shifts = 0.0d0
    IF (.NOT.ALL(ABS(position).LE.HUGE(1.0d0))) RETURN

    ! base(a): lattice index (a whole number, kept as a real) of the image nearest to the
    ! centre of the box; m(a): how far to look on either side of it
    base = 0.0d0
    m = 0
    DO a = 1, 3
      IF (.NOT.lattice_axis(a)) CYCLE
      base(a) = -ANINT((position(a) - 0.5d0*(box(2*a-1)+box(2*a)))/lattice_length(a))
      IF (.NOT.(ABS(base(a)).LT.EXACT_INDEX_LIMIT)) RETURN
      m(a) = 1
    END DO

    DO i3 = -m(3), m(3)
      s(3) = (base(3)+DBLE(i3))*lattice_length(3)
      gap2_3 = axis_gap(3, s(3))**2
      DO i2 = -m(2), m(2)
        s(2) = (base(2)+DBLE(i2))*lattice_length(2)
        gap2_23 = gap2_3 + axis_gap(2, s(2))**2
        DO i1 = -m(1), m(1)
          s(1) = (base(1)+DBLE(i1))*lattice_length(1)
          gap2 = gap2_23 + axis_gap(1, s(1))**2
          ! written so that a gap that is not a number is not an image
          IF (.NOT.(gap2.LE.radius*radius)) CYCLE
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
  ! neighbour test of the halo graph). Both boxes lie inside one lattice cell, so the nearest
  ! image along an axis is b itself or b shifted by one period either way.
  !
  ! Symmetric bit for bit: f(a,b) and f(b,a) perform the same floating-point operations. Per
  ! axis d1 = a_lo - b_hi and d2 = b_lo - a_hi are formed first; the gap to the image of b
  ! shifted by s is MAX(0, d1 - s, d2 + s). Exchanging a and b exchanges d1 and d2, which maps
  ! the shift s onto -s, and both signs are taken. (Shifting a box before subtracting, as in
  ! a_lo - (b_hi + s), rounds differently on the two sides: two ranks could disagree on a
  ! pair whose gap equals delta.)
  PURE LOGICAL FUNCTION EL_PERIODIC_BOXES_WITHIN_DELTA(a, b, delta)

    REAL*8, INTENT(IN) :: a(6), b(6), delta
    REAL*8 :: gap(3), d1, d2, len
    INTEGER :: k

    DO k = 1, 3
      d1 = a(2*k-1) - b(2*k)
      d2 = b(2*k-1) - a(2*k)
      gap(k) = MAX(0.0d0, d1, d2)
      IF (.NOT.lattice_axis(k)) CYCLE
      len = lattice_length(k)
      gap(k) = MIN(gap(k), MAX(0.0d0, d1 - len, d2 + len), MAX(0.0d0, d1 + len, d2 - len))
    END DO
    EL_PERIODIC_BOXES_WITHIN_DELTA = SQRT(SUM(gap*gap)).LE.delta

  END FUNCTION EL_PERIODIC_BOXES_WITHIN_DELTA

  ! Number of marked points on a periodic plane of the domain: the two faces of domain_box
  ! normal to each periodic axis, within EL_PERIODIC_TOLERANCE periods. Used to verify that a
  ! mask of kernel-cutting boundary faces holds no periodic plane (EL_TRANSFER). domain_box is
  ! the bounding box of the whole fluid mesh. -1 when coor does not cover the mask.
  PURE INTEGER FUNCTION EL_PERIODIC_MARKED_ON_PLANES(coor, marked, domain_box) RESULT(count)

    REAL*8, INTENT(IN) :: coor(:,:), domain_box(6)
    LOGICAL, INTENT(IN) :: marked(:)
    INTEGER :: i, a
    REAL*8 :: tol

    count = -1
    IF (SIZE(coor,1).LT.3 .OR. SIZE(coor,2).LT.SIZE(marked)) RETURN
    count = 0
    DO i = 1, SIZE(marked)
      IF (.NOT.marked(i)) CYCLE
      DO a = 1, 3
        IF (.NOT.lattice_axis(a)) CYCLE
        tol = EL_PERIODIC_TOLERANCE*lattice_length(a)
        IF (ABS(coor(a,i)-domain_box(2*a-1)).LE.tol .OR. ABS(coor(a,i)-domain_box(2*a)).LE.tol) THEN
          count = count + 1
          EXIT
        END IF
      END DO
    END DO

  END FUNCTION EL_PERIODIC_MARKED_ON_PLANES

END MODULE EL_PERIODIC
