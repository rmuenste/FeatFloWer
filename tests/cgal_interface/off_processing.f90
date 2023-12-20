PROGRAM OFF_PROCESSING
  IMPLICIT NONE

  ! Main program
  INTEGER :: processId, id, in

  DOUBLE PRECISION :: dist

  ! Set the coordinates of the point (x, y, z)
  DOUBLE PRECISION :: testPoint(3)

  ! Initialize processId
  processId = 1

  ! Initialize the C++ library
  CALL init_fc_easy(processId)

  testPoint = [0.0, 0.0, 1.5]

  ! Initialize variables
  dist = -1.0
  id = 0

  write(*,'(A)') "|==Test Setup=========================================================|"
  ! Print the test point coordinates
  write(*,'(A, 3F8.4)')"Test Point Coordinates: ", testPoint
  write(*,'(A)') "Testing distance and containment..."

  ! Print section header for distance computation
  write(*,'(A)') "|==Distance===========================================================|"

  ! Call library function to compute distance to the mesh and get the corresponding id
  CALL getdistanceid(testPoint(1), testPoint(2), testPoint(3), dist, id)

  ! Print the computed distance
  write(*,'(A, F8.4)')"Distance should be 1.0 and is = ", dist

  ! Print section header for containment test
  write(*,'(A)') "|==Containment Test===================================================|"

  ! Call library function to check if the point is contained within the mesh
  CALL isinelementid(testPoint(1), testPoint(2), testPoint(3), id, in)

  ! Print the result of the containment test
  write(*,'(A, I0)') "Containment test result should be 0 and is = ", in

END PROGRAM
