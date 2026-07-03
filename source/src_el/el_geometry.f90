MODULE EL_GEOMETRY

  IMPLICIT NONE

  REAL*8 :: el_domain_box(6) = 0.0d0
  LOGICAL :: el_domain_box_set = .FALSE.

CONTAINS

  SUBROUTINE EL_SET_DOMAIN_BOX(xmin, xmax, ymin, ymax, zmin, zmax)

    REAL*8, INTENT(IN) :: xmin, xmax, ymin, ymax, zmin, zmax

    el_domain_box = (/xmin, xmax, ymin, ymax, zmin, zmax/)
    el_domain_box_set = .TRUE.

  END SUBROUTINE EL_SET_DOMAIN_BOX

  REAL*8 FUNCTION EL_WALL_DISTANCE(position)

    REAL*8, INTENT(IN) :: position(3)
    REAL*8 :: normal(3)

    CALL EL_WALL_DISTANCE_AND_NORMAL(position, EL_WALL_DISTANCE, normal)

  END FUNCTION EL_WALL_DISTANCE

  SUBROUTINE EL_WALL_DISTANCE_AND_NORMAL(position, distance, normal)

    REAL*8, INTENT(IN) :: position(3)
    REAL*8, INTENT(OUT) :: distance, normal(3)
    REAL*8 :: distances(6)
    INTEGER :: iwall

    normal = 0.0d0
    IF (.NOT.el_domain_box_set) THEN
      distance = HUGE(1.0d0)
      RETURN
    END IF

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

  END SUBROUTINE EL_WALL_DISTANCE_AND_NORMAL

  REAL*8 FUNCTION EL_DOMAIN_Z_CENTER()

    IF (el_domain_box_set) THEN
      EL_DOMAIN_Z_CENTER = 0.5d0*(el_domain_box(5) + el_domain_box(6))
    ELSE
      EL_DOMAIN_Z_CENTER = 0.5d0
    END IF

  END FUNCTION EL_DOMAIN_Z_CENTER

  ! Future non-box geometries should route this module through the PE/CGAL
  ! distance backend, e.g. getdistanceid in the FullC0ntact Fortran/C++ interface.

END MODULE EL_GEOMETRY
