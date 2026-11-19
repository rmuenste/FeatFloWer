!=========================================================================
! EL_INLET_TABLE - piecewise-linear inlet velocity table U(t) of boundary type Inflow975
!
! SimPar@ELInletVelocityTable = "<file>" (EL_CONFIG) names a text file of rows "t U"
! (seconds, m/s; the SUPERFICIAL velocity like the ramp keys). Blank lines and lines starting
! with '#' or '!' are skipped; the two numbers may be separated by blanks or a comma; anything
! after them on a row is ignored. Refused (ierr /= 0 with a message): a missing or unreadable
! file, a row that does not hold two numbers, a non-finite value, times that do not increase
! strictly, an empty table. U(t) is the linear interpolation between the rows and constant
! beyond the first and the last row.
!
! Pure table code: no MPI, no deck access. EL_CONFIG loads the table on every rank from the
! same file during the deck validation (every rank parses the deck itself), so every rank
! holds the same rows or every rank stops with the same message.
!=========================================================================
MODULE EL_INLET_TABLE

  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: tElInletTable, EL_INLET_TABLE_READ, EL_INLET_TABLE_VALUE, EL_INLET_TABLE_PATH

  TYPE tElInletTable
    INTEGER :: n = 0
    REAL*8, ALLOCATABLE :: t(:), u(:)
  END TYPE tElInletTable

  INTEGER, PARAMETER :: MAX_LINE = 512

CONTAINS

  ! The file name of the deck value: the text between quotes ("..." or '...'), else the first
  ! blank-delimited word (a list-directed READ would stop at the first '/' of a path).
  FUNCTION EL_INLET_TABLE_PATH(text) RESULT(path)
    CHARACTER(LEN=*), INTENT(IN) :: text
    CHARACTER(LEN=256) :: path
    CHARACTER(LEN=LEN(text)) :: work
    INTEGER :: iend

    path = ' '
    work = ADJUSTL(text)
    IF (LEN_TRIM(work).EQ.0) RETURN
    IF (work(1:1).EQ.'"' .OR. work(1:1).EQ."'") THEN
      iend = INDEX(work(2:), work(1:1))
      IF (iend.EQ.0) THEN
        path = work(2:)
      ELSE
        path = work(2:iend)
      END IF
    ELSE
      iend = INDEX(work, ' ')
      IF (iend.EQ.0) iend = LEN(work) + 1
      path = work(1:iend-1)
    END IF
  END FUNCTION EL_INLET_TABLE_PATH

  ! Read and check the table. On success ierr = 0 and table holds the rows; otherwise ierr > 0,
  ! table%n = 0 and message says why (file name and row number included).
  SUBROUTINE EL_INLET_TABLE_READ(filename, table, ierr, message)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    TYPE(tElInletTable), INTENT(INOUT) :: table
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message

    CHARACTER(LEN=MAX_LINE) :: line
    INTEGER :: iunit, ios, nrow, pass, lineno, i
    REAL*8 :: tv, uv
    LOGICAL :: exists

    ierr = 0
    message = ' '
    table%n = 0
    IF (ALLOCATED(table%t)) DEALLOCATE(table%t)
    IF (ALLOCATED(table%u)) DEALLOCATE(table%u)

    INQUIRE(FILE=TRIM(filename), EXIST=exists)
    IF (.NOT.exists .OR. LEN_TRIM(filename).EQ.0) THEN
      ierr = 1
      message = 'ELInletVelocityTable: file not found: '//TRIM(filename)
      RETURN
    END IF

    ! pass 1 counts and checks the rows, pass 2 stores them
    DO pass = 1, 2
      OPEN(NEWUNIT=iunit, FILE=TRIM(filename), STATUS='OLD', ACTION='READ', IOSTAT=ios)
      IF (ios.NE.0) THEN
        ierr = 2
        message = 'ELInletVelocityTable: cannot open '//TRIM(filename)
        RETURN
      END IF
      nrow = 0
      lineno = 0
      DO
        READ(iunit, '(A)', IOSTAT=ios) line
        IF (ios.NE.0) EXIT
        lineno = lineno + 1
        line = ADJUSTL(line)
        IF (LEN_TRIM(line).EQ.0) CYCLE
        IF (line(1:1).EQ.'#' .OR. line(1:1).EQ.'!') CYCLE
        DO i = 1, LEN_TRIM(line)
          IF (line(i:i).EQ.',' .OR. line(i:i).EQ.CHAR(9)) line(i:i) = ' '
        END DO
        READ(line, *, IOSTAT=ios) tv, uv
        IF (ios.NE.0) THEN
          ierr = 3
          WRITE(message,'(A,A,A,I0,A)') 'ELInletVelocityTable: ', TRIM(filename), ' row ', &
            lineno, ' does not hold two numbers "t U"'
          CLOSE(iunit)
          RETURN
        END IF
        IF (.NOT.(IEEE_IS_FINITE(tv) .AND. IEEE_IS_FINITE(uv))) THEN
          ierr = 4
          WRITE(message,'(A,A,A,I0,A)') 'ELInletVelocityTable: ', TRIM(filename), ' row ', &
            lineno, ' holds a non-finite value'
          CLOSE(iunit)
          RETURN
        END IF
        nrow = nrow + 1
        IF (pass.EQ.2) THEN
          table%t(nrow) = tv
          table%u(nrow) = uv
          IF (nrow.GT.1) THEN
            IF (.NOT.(table%t(nrow).GT.table%t(nrow-1))) THEN
              ierr = 5
              WRITE(message,'(A,A,A,I0,A)') 'ELInletVelocityTable: ', TRIM(filename), &
                ' row ', lineno, ': times must increase strictly'
              CLOSE(iunit)
              DEALLOCATE(table%t, table%u)
              RETURN
            END IF
          END IF
        END IF
      END DO
      CLOSE(iunit)
      IF (pass.EQ.1) THEN
        IF (nrow.EQ.0) THEN
          ierr = 6
          message = 'ELInletVelocityTable: no rows in '//TRIM(filename)
          RETURN
        END IF
        ALLOCATE(table%t(nrow), table%u(nrow))
      END IF
    END DO
    table%n = nrow

  END SUBROUTINE EL_INLET_TABLE_READ

  ! U(t): linear between the rows, constant beyond the ends (the caller ensures table%n >= 1).
  PURE REAL*8 FUNCTION EL_INLET_TABLE_VALUE(table, t)
    TYPE(tElInletTable), INTENT(IN) :: table
    REAL*8, INTENT(IN) :: t
    INTEGER :: lo, hi, mid

    IF (t.LE.table%t(1)) THEN
      EL_INLET_TABLE_VALUE = table%u(1)
      RETURN
    END IF
    IF (t.GE.table%t(table%n)) THEN
      EL_INLET_TABLE_VALUE = table%u(table%n)
      RETURN
    END IF
    ! table%t(lo) < t < table%t(hi)
    lo = 1
    hi = table%n
    DO WHILE (hi - lo.GT.1)
      mid = (lo + hi)/2
      IF (table%t(mid).LE.t) THEN
        lo = mid
      ELSE
        hi = mid
      END IF
    END DO
    EL_INLET_TABLE_VALUE = table%u(lo) + (table%u(hi) - table%u(lo))* &
                           (t - table%t(lo))/(table%t(hi) - table%t(lo))
  END FUNCTION EL_INLET_TABLE_VALUE

END MODULE EL_INLET_TABLE
