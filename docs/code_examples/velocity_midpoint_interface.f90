! Simple interface for velocity midpoint evaluation
! Easy-to-use wrappers for common use cases

MODULE VelocityMidpoint_Interface
!
! Purpose: Provide easy-to-use interfaces for velocity evaluation at element centers
!          and other common post-processing tasks
!

USE var_QuadScalar, ONLY: TQuadScalar
USE mg_mesh_TypeDef, ONLY: mg_mesh
USE var_LScalars, ONLY: ILEV

IMPLICIT NONE

CONTAINS

!===========================================================================

SUBROUTINE Call_VelocityMidpointEvaluation(QuadSc, output_vtk, mfile)
!
! Purpose: Complete velocity midpoint evaluation with optional VTK output
!          This is the main interface that most users will want to call
!
IMPLICIT NONE

TYPE(TQuadScalar), INTENT(IN) :: QuadSc
LOGICAL, INTENT(IN) :: output_vtk
INTEGER, INTENT(IN) :: mfile

! Local arrays
REAL*8, ALLOCATABLE :: velocity_midpoints(:,:)  ! (3, NEL)
REAL*8, ALLOCATABLE :: element_centers(:,:)     ! (3, NEL)
INTEGER :: NEL

!---------------------------------------------------------------------------

NEL = mg_mesh%level(ILEV)%NEL

! Allocate arrays
ALLOCATE(velocity_midpoints(3, NEL))
ALLOCATE(element_centers(3, NEL))

! Evaluate velocities at element centers
CALL Evaluate_Velocity_Midpoint(QuadSc, velocity_midpoints, element_centers, mfile)

! Optional VTK output
IF (output_vtk) THEN
  CALL Write_Velocity_Midpoints_VTK(velocity_midpoints, element_centers, &
                                   'velocity_midpoints', mfile)
END IF

! Report statistics
CALL Report_VelocityMidpoint_Statistics(velocity_midpoints, mfile)

! Clean up
DEALLOCATE(velocity_midpoints)
DEALLOCATE(element_centers)

END SUBROUTINE Call_VelocityMidpointEvaluation

!===========================================================================

SUBROUTINE Report_VelocityMidpoint_Statistics(velocity_midpoints, mfile)
!
! Purpose: Report basic statistics about element-centered velocities
!
IMPLICIT NONE

REAL*8, INTENT(IN) :: velocity_midpoints(:,:)
INTEGER, INTENT(IN) :: mfile

REAL*8 :: u_min, u_max, u_avg, u_rms
REAL*8 :: v_min, v_max, v_avg, v_rms
REAL*8 :: w_min, w_max, w_avg, w_rms
REAL*8 :: mag_min, mag_max, mag_avg, mag_rms
REAL*8 :: magnitude
INTEGER :: NEL, IEL

!---------------------------------------------------------------------------

NEL = SIZE(velocity_midpoints, 2)

! Initialize
u_min = HUGE(1D0); u_max = -HUGE(1D0); u_avg = 0D0; u_rms = 0D0
v_min = HUGE(1D0); v_max = -HUGE(1D0); v_avg = 0D0; v_rms = 0D0
w_min = HUGE(1D0); w_max = -HUGE(1D0); w_avg = 0D0; w_rms = 0D0
mag_min = HUGE(1D0); mag_max = -HUGE(1D0); mag_avg = 0D0; mag_rms = 0D0

! Compute statistics
DO IEL = 1, NEL
  ! U component
  u_min = MIN(u_min, velocity_midpoints(1, IEL))
  u_max = MAX(u_max, velocity_midpoints(1, IEL))
  u_avg = u_avg + velocity_midpoints(1, IEL)
  u_rms = u_rms + velocity_midpoints(1, IEL)**2

  ! V component
  v_min = MIN(v_min, velocity_midpoints(2, IEL))
  v_max = MAX(v_max, velocity_midpoints(2, IEL))
  v_avg = v_avg + velocity_midpoints(2, IEL)
  v_rms = v_rms + velocity_midpoints(2, IEL)**2

  ! W component
  w_min = MIN(w_min, velocity_midpoints(3, IEL))
  w_max = MAX(w_max, velocity_midpoints(3, IEL))
  w_avg = w_avg + velocity_midpoints(3, IEL)
  w_rms = w_rms + velocity_midpoints(3, IEL)**2

  ! Magnitude
  magnitude = SQRT(SUM(velocity_midpoints(1:3, IEL)**2))
  mag_min = MIN(mag_min, magnitude)
  mag_max = MAX(mag_max, magnitude)
  mag_avg = mag_avg + magnitude
  mag_rms = mag_rms + magnitude**2
END DO

! Normalize
u_avg = u_avg / NEL
v_avg = v_avg / NEL
w_avg = w_avg / NEL
mag_avg = mag_avg / NEL

u_rms = SQRT(u_rms / NEL)
v_rms = SQRT(v_rms / NEL)
w_rms = SQRT(w_rms / NEL)
mag_rms = SQRT(mag_rms / NEL)

! Report
WRITE(mfile, '(/A)') '=== Element-Centered Velocity Statistics ==='
WRITE(mfile, '(A,I0)') 'Number of elements: ', NEL
WRITE(mfile, '(/A)') 'U-component:'
WRITE(mfile, '(A,3(E12.4,2X))') '  Min/Max/Avg: ', u_min, u_max, u_avg
WRITE(mfile, '(A,E12.4)') '  RMS:         ', u_rms
WRITE(mfile, '(/A)') 'V-component:'
WRITE(mfile, '(A,3(E12.4,2X))') '  Min/Max/Avg: ', v_min, v_max, v_avg
WRITE(mfile, '(A,E12.4)') '  RMS:         ', v_rms
WRITE(mfile, '(/A)') 'W-component:'
WRITE(mfile, '(A,3(E12.4,2X))') '  Min/Max/Avg: ', w_min, w_max, w_avg
WRITE(mfile, '(A,E12.4)') '  RMS:         ', w_rms
WRITE(mfile, '(/A)') 'Magnitude:'
WRITE(mfile, '(A,3(E12.4,2X))') '  Min/Max/Avg: ', mag_min, mag_max, mag_avg
WRITE(mfile, '(A,E12.4)') '  RMS:         ', mag_rms
WRITE(mfile, '(A/)') '============================================='

END SUBROUTINE Report_VelocityMidpoint_Statistics

!===========================================================================

SUBROUTINE Find_Elements_With_HighVelocity(QuadSc, velocity_threshold, &
                                           high_velocity_elements, num_found, mfile)
!
! Purpose: Find elements where velocity magnitude exceeds a threshold
!          Useful for identifying regions of high shear, recirculation, etc.
!
IMPLICIT NONE

TYPE(TQuadScalar), INTENT(IN) :: QuadSc
REAL*8, INTENT(IN) :: velocity_threshold
INTEGER, INTENT(OUT) :: high_velocity_elements(:)
INTEGER, INTENT(OUT) :: num_found
INTEGER, INTENT(IN) :: mfile

! Local arrays
REAL*8, ALLOCATABLE :: velocity_midpoints(:,:)
REAL*8, ALLOCATABLE :: element_centers(:,:)
REAL*8 :: magnitude
INTEGER :: NEL, IEL

!---------------------------------------------------------------------------

NEL = mg_mesh%level(ILEV)%NEL

! Allocate arrays
ALLOCATE(velocity_midpoints(3, NEL))
ALLOCATE(element_centers(3, NEL))

! Evaluate velocities
CALL Evaluate_Velocity_Midpoint(QuadSc, velocity_midpoints, element_centers, mfile)

! Find high velocity elements
num_found = 0
DO IEL = 1, NEL
  magnitude = SQRT(SUM(velocity_midpoints(1:3, IEL)**2))
  
  IF (magnitude .GT. velocity_threshold) THEN
    num_found = num_found + 1
    IF (num_found .LE. SIZE(high_velocity_elements)) THEN
      high_velocity_elements(num_found) = IEL
    END IF
  END IF
END DO

! Report results
WRITE(mfile, '(A,E12.4)') 'Velocity threshold: ', velocity_threshold
WRITE(mfile, '(A,I0,A,I0,A)') 'Found ', num_found, ' elements (out of ', NEL, &
                               ') with velocity above threshold'

IF (num_found .GT. SIZE(high_velocity_elements)) THEN
  WRITE(mfile, '(A,I0,A)') 'WARNING: Output array too small, only first ', &
                           SIZE(high_velocity_elements), ' elements returned'
END IF

! Clean up
DEALLOCATE(velocity_midpoints)
DEALLOCATE(element_centers)

END SUBROUTINE Find_Elements_With_HighVelocity

!===========================================================================

SUBROUTINE Evaluate_Velocity_Along_Line(QuadSc, start_point, end_point, &
                                        num_points, velocities, positions, mfile)
!
! Purpose: Evaluate velocity along a line through the domain
!          Useful for creating velocity profiles, line plots, etc.
!
! Note: This is a simplified version that assumes the line passes through
!       elements in a straightforward way. A more robust implementation
!       would need element search and intersection algorithms.
!
IMPLICIT NONE

TYPE(TQuadScalar), INTENT(IN) :: QuadSc
REAL*8, INTENT(IN) :: start_point(3), end_point(3)
INTEGER, INTENT(IN) :: num_points
REAL*8, INTENT(OUT) :: velocities(:,:)   ! (3, num_points)
REAL*8, INTENT(OUT) :: positions(:,:)    ! (3, num_points)
INTEGER, INTENT(IN) :: mfile

! Local variables
REAL*8 :: direction(3), step_size, t
INTEGER :: i_point
REAL*8 :: current_position(3), current_velocity(3)

!---------------------------------------------------------------------------

! Calculate direction and step size
direction = end_point - start_point
step_size = 1D0 / DBLE(num_points - 1)

WRITE(mfile, '(A)') 'Evaluating velocity along line...'
WRITE(mfile, '(A,3F10.4)') 'Start point: ', start_point
WRITE(mfile, '(A,3F10.4)') 'End point:   ', end_point
WRITE(mfile, '(A,I0)') 'Number of evaluation points: ', num_points

! Evaluate at each point along the line
DO i_point = 1, num_points
  t = DBLE(i_point - 1) * step_size
  current_position = start_point + t * direction
  
  ! Store position
  positions(1:3, i_point) = current_position
  
  ! Find element containing this point and evaluate velocity
  ! This is simplified - in practice would need robust element search
  CALL Find_Element_Containing_Point(current_position, current_velocity)
  
  ! Store velocity
  velocities(1:3, i_point) = current_velocity
END DO

WRITE(mfile, '(A)') 'Line evaluation complete.'

END SUBROUTINE Evaluate_Velocity_Along_Line

!===========================================================================

SUBROUTINE Find_Element_Containing_Point(point, velocity)
!
! Purpose: Find element containing a given point and evaluate velocity there
!          Simplified placeholder implementation
!
IMPLICIT NONE

REAL*8, INTENT(IN) :: point(3)
REAL*8, INTENT(OUT) :: velocity(3)

! Placeholder implementation
! In practice, this would need:
! 1. Efficient element search (e.g., using bounding boxes)
! 2. Inverse isoparametric mapping to find reference coordinates
! 3. Proper velocity evaluation using those coordinates

velocity = 0D0  ! Default to zero if not implemented

END SUBROUTINE Find_Element_Containing_Point

END MODULE VelocityMidpoint_Interface