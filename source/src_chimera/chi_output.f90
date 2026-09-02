!=========================================================================
! CHI_OUTPUT - legacy-VTK dump of a submesh solution, Layer M.
! Vertex velocities (the first nvt Q2 dofs) and the element-centroid
! pressure (first P1 dof of each element).  Enabled per deck via
! SimPar@ChimeraWriteVTK.
!=========================================================================
MODULE CHI_OUTPUT

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: CHI_WRITE_SUBMESH_VTK

CONTAINS

  SUBROUTINE CHI_WRITE_SUBMESH_VTK(fname, nvt, nel, kvert, dcorvg, u, v, w, p)
    CHARACTER(*), INTENT(IN) :: fname
    INTEGER, INTENT(IN) :: nvt, nel, kvert(8,*)
    REAL*8,  INTENT(IN) :: dcorvg(3,*), u(*), v(*), w(*), p(*)

    INTEGER :: iu, i, e, ios

    iu = 771
    OPEN(UNIT=iu, FILE=TRIM(fname), STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    IF (ios .NE. 0) THEN
      WRITE(*,*) 'CHI_WRITE_SUBMESH_VTK: cannot open ', TRIM(fname)
      RETURN
    END IF
    WRITE(iu,'(A)') '# vtk DataFile Version 3.0'
    WRITE(iu,'(A)') 'Chimera submesh solution'
    WRITE(iu,'(A)') 'ASCII'
    WRITE(iu,'(A)') 'DATASET UNSTRUCTURED_GRID'
    WRITE(iu,'(A,I0,A)') 'POINTS ', nvt, ' double'
    DO i = 1, nvt
      WRITE(iu,'(3ES22.14)') dcorvg(1:3,i)
    END DO
    WRITE(iu,'(A,I0,1X,I0)') 'CELLS ', nel, 9*nel
    DO e = 1, nel
      WRITE(iu,'(I0,8(1X,I0))') 8, kvert(1:8,e) - 1
    END DO
    WRITE(iu,'(A,I0)') 'CELL_TYPES ', nel
    DO e = 1, nel
      WRITE(iu,'(I0)') 12
    END DO
    WRITE(iu,'(A,I0)') 'POINT_DATA ', nvt
    WRITE(iu,'(A)') 'VECTORS velocity double'
    DO i = 1, nvt
      WRITE(iu,'(3ES22.14)') u(i), v(i), w(i)
    END DO
    WRITE(iu,'(A,I0)') 'CELL_DATA ', nel
    WRITE(iu,'(A)') 'SCALARS pressure double 1'
    WRITE(iu,'(A)') 'LOOKUP_TABLE default'
    DO e = 1, nel
      WRITE(iu,'(ES22.14)') p(4*(e-1)+1)
    END DO
    CLOSE(iu)
  END SUBROUTINE CHI_WRITE_SUBMESH_VTK

END MODULE CHI_OUTPUT
