MODULE EL_GEOMETRY

  IMPLICIT NONE

  REAL*8 :: el_domain_box(6) = 0.0d0
  LOGICAL :: el_domain_box_set = .FALSE.
  REAL*8 :: el_domain_cylinder_center(3) = 0.0d0
  REAL*8 :: el_domain_cylinder_radius = -1.0d0
  CHARACTER(LEN=8) :: el_domain_cylinder_axis = 'x'
  LOGICAL :: el_domain_cylinder_set = .FALSE.
  LOGICAL :: el_domain_set = .FALSE.

CONTAINS

  SUBROUTINE EL_SET_DOMAIN_BOX(xmin, xmax, ymin, ymax, zmin, zmax)

    REAL*8, INTENT(IN) :: xmin, xmax, ymin, ymax, zmin, zmax

    el_domain_box = (/xmin, xmax, ymin, ymax, zmin, zmax/)
    el_domain_box_set = .TRUE.
    el_domain_cylinder_set = .FALSE.
    el_domain_set = .TRUE.

  END SUBROUTINE EL_SET_DOMAIN_BOX

  SUBROUTINE EL_SET_DOMAIN_CYLINDER(center, radius, axis)

    REAL*8, INTENT(IN) :: center(3), radius
    CHARACTER(LEN=*), INTENT(IN) :: axis

    el_domain_cylinder_center = center
    el_domain_cylinder_radius = radius
    el_domain_cylinder_axis = axis
    CALL EL_LOWERCASE(el_domain_cylinder_axis)
    el_domain_cylinder_set = .TRUE.
    el_domain_box_set = .FALSE.
    el_domain_set = .TRUE.

  END SUBROUTINE EL_SET_DOMAIN_CYLINDER

  REAL*8 FUNCTION EL_WALL_DISTANCE(position)

    REAL*8, INTENT(IN) :: position(3)
    REAL*8 :: normal(3)

    CALL EL_WALL_DISTANCE_AND_NORMAL(position, EL_WALL_DISTANCE, normal)

  END FUNCTION EL_WALL_DISTANCE

  SUBROUTINE EL_WALL_DISTANCE_AND_NORMAL(position, distance, normal)

    REAL*8, INTENT(IN) :: position(3)
    REAL*8, INTENT(OUT) :: distance, normal(3)
    REAL*8 :: distances(6), radial(3), rnorm
    INTEGER :: iwall

    normal = 0.0d0
    IF (.NOT.el_domain_set) THEN
      distance = HUGE(1.0d0)
      RETURN
    END IF

    IF (el_domain_cylinder_set) THEN
      radial = position - el_domain_cylinder_center
      SELECT CASE (TRIM(el_domain_cylinder_axis))
      CASE ('x')
        radial(1) = 0.0d0
      CASE ('y')
        radial(2) = 0.0d0
      CASE DEFAULT
        radial(3) = 0.0d0
      END SELECT
      rnorm = SQRT(SUM(radial**2))
      distance = el_domain_cylinder_radius - rnorm
      IF (rnorm.GT.0.0d0) THEN
        normal = -radial/rnorm
      ELSE
        SELECT CASE (TRIM(el_domain_cylinder_axis))
        CASE ('x')
          normal(2) = -1.0d0
        CASE ('y')
          normal(1) = -1.0d0
        CASE DEFAULT
          normal(1) = -1.0d0
        END SELECT
      END IF
    ELSE
      distances = (/position(1)-el_domain_box(1), el_domain_box(2)-position(1), &
                    position(2)-el_domain_box(3), el_domain_box(4)-position(2), &
                    position(3)-el_domain_box(5), el_domain_box(6)-position(3)/)
      iwall = MINLOC(distances, DIM=1)
      distance = distances(iwall)

      SELECT CASE (iwall)
      CASE (1)
        normal(1) = 1.0d0
      CASE (2)
        normal(1) = -1.0d0
      CASE (3)
        normal(2) = 1.0d0
      CASE (4)
        normal(2) = -1.0d0
      CASE (5)
        normal(3) = 1.0d0
      CASE (6)
        normal(3) = -1.0d0
      END SELECT
    END IF

  END SUBROUTINE EL_WALL_DISTANCE_AND_NORMAL

  REAL*8 FUNCTION EL_DOMAIN_Z_CENTER()

    IF (el_domain_box_set) THEN
      EL_DOMAIN_Z_CENTER = 0.5d0*(el_domain_box(5) + el_domain_box(6))
    ELSE
      EL_DOMAIN_Z_CENTER = 0.5d0
    END IF

  END FUNCTION EL_DOMAIN_Z_CENTER

  SUBROUTINE EL_LOWERCASE(text)

    CHARACTER(LEN=*), INTENT(INOUT) :: text
    INTEGER :: i, code

    DO i=1,LEN_TRIM(text)
      code = IACHAR(text(i:i))
      IF (code.GE.IACHAR('A') .AND. code.LE.IACHAR('Z')) THEN
        text(i:i) = ACHAR(code + IACHAR('a') - IACHAR('A'))
      END IF
    END DO

  END SUBROUTINE EL_LOWERCASE

  ! More complex geometries should route through the PE/CGAL distance backend,
  ! e.g. getdistanceid in the FullC0ntact Fortran/C++ interface.

END MODULE EL_GEOMETRY
