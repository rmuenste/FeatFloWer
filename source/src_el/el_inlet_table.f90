!=========================================================================
! EL_INLET_TABLE - piecewise-linear inlet velocity table U(t) of boundary type Inflow975
!
! SimPar@ELInletVelocityTable = "<file>" (EL_CONFIG) names a text file of rows "t U"
! (seconds, m/s; the SUPERFICIAL velocity like the ramp keys). Blank lines (spaces and tabs
! only) and lines whose first non-blank character is '#' or '!' are skipped (a tab is a blank
! everywhere); the two numbers may be separated by blanks, tabs or a comma; anything
! after them on a row is ignored. Refused (ierr /= 0 with a message): a missing or unreadable
! file, a read error (distinguished from the end of the file), a row longer than MAX_LINE, a
! row whose first two tokens are not both numbers (a token is a number when it is an optional
! sign, digits with at most one decimal point and an optional exponent e/E/d/D with digits;
! "0.5 /" or "0.5" is refused, not read list-directed with a value left over from an earlier
! row), a non-finite value, times that do not increase strictly, an empty table. U(t) is the
! linear interpolation between the rows and constant beyond the first and the last row.
!
! Three separate parts:
!   EL_INLET_TABLE_READ         rank-local: one snapshot of the file (opened once, every line
!                               read into memory, closed), then parsed and checked; no MPI
!   EL_INLET_TABLE_LOAD_SHARED  distribution: rank 0 of the communicator reads, the status,
!                               the message and the rows are broadcast, so every rank holds the
!                               same rows or the same refusal (EL_CONFIG stops every rank then)
!   EL_INLET_TABLE_VALUE        evaluation U(t)
! The table is a start-up setting: EL_CONFIG loads it once during the deck validation of the
! application; a reload of the parameter file (Reload_SimPar, generic parser) does not touch
! it (stated in the echo of the table).
!=========================================================================
MODULE EL_INLET_TABLE

  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  USE MPI

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: tElInletTable, EL_INLET_TABLE_READ, EL_INLET_TABLE_LOAD_SHARED, &
            EL_INLET_TABLE_VALUE, EL_INLET_TABLE_PATH, EL_INLET_TABLE_IS_NUMBER

  TYPE tElInletTable
    INTEGER :: n = 0
    REAL*8, ALLOCATABLE :: t(:), u(:)
  END TYPE tElInletTable

  INTEGER, PARAMETER :: MAX_LINE = 1024

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

  ! True when tok (no blanks) is a decimal number: [sign] digits [. digits] [exponent], at least
  ! one digit in the mantissa, exponent letter e/E/d/D followed by [sign] and at least one digit.
  PURE LOGICAL FUNCTION EL_INLET_TABLE_IS_NUMBER(tok)
    CHARACTER(LEN=*), INTENT(IN) :: tok
    INTEGER :: i, n, ndig
    EL_INLET_TABLE_IS_NUMBER = .FALSE.
    n = LEN_TRIM(tok)
    IF (n.EQ.0) RETURN
    i = 1
    IF (tok(i:i).EQ.'+' .OR. tok(i:i).EQ.'-') i = i + 1
    ndig = 0
    DO WHILE (i.LE.n)
      IF (.NOT.IS_DIGIT(tok(i:i))) EXIT
      ndig = ndig + 1
      i = i + 1
    END DO
    IF (i.LE.n) THEN
      IF (tok(i:i).EQ.'.') THEN
        i = i + 1
        DO WHILE (i.LE.n)
          IF (.NOT.IS_DIGIT(tok(i:i))) EXIT
          ndig = ndig + 1
          i = i + 1
        END DO
      END IF
    END IF
    IF (ndig.EQ.0) RETURN
    IF (i.LE.n) THEN
      IF (INDEX('eEdD', tok(i:i)).EQ.0) RETURN
      i = i + 1
      IF (i.LE.n) THEN
        IF (tok(i:i).EQ.'+' .OR. tok(i:i).EQ.'-') i = i + 1
      END IF
      ndig = 0
      DO WHILE (i.LE.n)
        IF (.NOT.IS_DIGIT(tok(i:i))) EXIT
        ndig = ndig + 1
        i = i + 1
      END DO
      IF (ndig.EQ.0) RETURN
    END IF
    EL_INLET_TABLE_IS_NUMBER = (i.GT.n)
  CONTAINS
    PURE LOGICAL FUNCTION IS_DIGIT(c)
      CHARACTER(LEN=1), INTENT(IN) :: c
      IS_DIGIT = (c.GE.'0' .AND. c.LE.'9')
    END FUNCTION IS_DIGIT
  END FUNCTION EL_INLET_TABLE_IS_NUMBER

  ! Read and check the table (rank-local, no MPI). The file is opened once and every line is
  ! read into memory before anything is parsed (one snapshot: a file changed while it is read
  ! cannot give two different row counts). On success ierr = 0 and table holds the rows;
  ! otherwise ierr > 0, table%n = 0 and message says why (file name and row number included).
  SUBROUTINE EL_INLET_TABLE_READ(filename, table, ierr, message)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    TYPE(tElInletTable), INTENT(INOUT) :: table
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message

    CHARACTER(LEN=MAX_LINE), ALLOCATABLE :: lines(:)
    INTEGER :: nlines, nrow, k
    REAL*8 :: tv, uv
    LOGICAL :: data_row
    REAL*8, ALLOCATABLE :: t(:), u(:)

    ierr = 0
    message = ' '
    CALL CLEAR(table)

    CALL READ_SNAPSHOT(filename, lines, nlines, ierr, message)
    IF (ierr.NE.0) RETURN

    ALLOCATE(t(nlines), u(nlines))
    nrow = 0
    DO k = 1, nlines
      CALL PARSE_ROW(lines(k), data_row, tv, uv, ierr)
      IF (ierr.EQ.3) THEN
        WRITE(message,'(A,A,A,I0,A)') 'ELInletVelocityTable: ', TRIM(filename), ' row ', k, &
          ' does not hold two numbers "t U"'
        RETURN
      ELSE IF (ierr.EQ.4) THEN
        WRITE(message,'(A,A,A,I0,A)') 'ELInletVelocityTable: ', TRIM(filename), ' row ', k, &
          ' holds a non-finite value'
        RETURN
      END IF
      IF (.NOT.data_row) CYCLE
      nrow = nrow + 1
      t(nrow) = tv
      u(nrow) = uv
      IF (nrow.GT.1) THEN
        IF (.NOT.(t(nrow).GT.t(nrow-1))) THEN
          ierr = 5
          WRITE(message,'(A,A,A,I0,A)') 'ELInletVelocityTable: ', TRIM(filename), &
            ' row ', k, ': times must increase strictly'
          RETURN
        END IF
      END IF
    END DO
    IF (nrow.EQ.0) THEN
      ierr = 6
      message = 'ELInletVelocityTable: no rows in '//TRIM(filename)
      RETURN
    END IF
    ALLOCATE(table%t(nrow), table%u(nrow))
    table%t = t(1:nrow)
    table%u = u(1:nrow)
    table%n = nrow
  END SUBROUTINE EL_INLET_TABLE_READ

  ! Rank 0 of comm reads the table (EL_INLET_TABLE_READ), every rank of comm receives the
  ! status, the message and on success the rows. Collective over comm: every rank calls it.
  SUBROUTINE EL_INLET_TABLE_LOAD_SHARED(filename, table, ierr, message, comm)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    TYPE(tElInletTable), INTENT(INOUT) :: table
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    INTEGER, INTENT(IN) :: comm
    INTEGER :: rank, mpierr, head(2)

    CALL MPI_Comm_rank(comm, rank, mpierr)
    ierr = 0
    message = ' '
    IF (rank.EQ.0) THEN
      CALL EL_INLET_TABLE_READ(filename, table, ierr, message)
    ELSE
      CALL CLEAR(table)
    END IF
    head = (/ ierr, table%n /)
    CALL MPI_Bcast(head, 2, MPI_INTEGER, 0, comm, mpierr)
    CALL MPI_Bcast(message, LEN(message), MPI_CHARACTER, 0, comm, mpierr)
    ierr = head(1)
    IF (ierr.NE.0) THEN
      CALL CLEAR(table)
      RETURN
    END IF
    IF (rank.NE.0) THEN
      table%n = head(2)
      ALLOCATE(table%t(table%n), table%u(table%n))
    END IF
    CALL MPI_Bcast(table%t, table%n, MPI_DOUBLE_PRECISION, 0, comm, mpierr)
    CALL MPI_Bcast(table%u, table%n, MPI_DOUBLE_PRECISION, 0, comm, mpierr)
  END SUBROUTINE EL_INLET_TABLE_LOAD_SHARED

  SUBROUTINE CLEAR(table)
    TYPE(tElInletTable), INTENT(INOUT) :: table
    table%n = 0
    IF (ALLOCATED(table%t)) DEALLOCATE(table%t)
    IF (ALLOCATED(table%u)) DEALLOCATE(table%u)
  END SUBROUTINE CLEAR

  ! Every line of the file in memory (one open). ierr 1 missing file, 2 cannot open, 7 read
  ! error (not the end of the file), 8 a line of MAX_LINE characters or more. Non-advancing
  ! reads: the end of a record ends a line, a full buffer without it is a line too long.
  SUBROUTINE READ_SNAPSHOT(filename, lines, nlines, ierr, message)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=MAX_LINE), ALLOCATABLE, INTENT(OUT) :: lines(:)
    INTEGER, INTENT(OUT) :: nlines, ierr
    CHARACTER(LEN=*), INTENT(INOUT) :: message
    CHARACTER(LEN=MAX_LINE), ALLOCATABLE :: grown(:)
    CHARACTER(LEN=MAX_LINE) :: line
    INTEGER :: iunit, ios, nchar
    LOGICAL :: exists

    nlines = 0
    ALLOCATE(lines(64))
    IF (LEN_TRIM(filename).EQ.0) THEN
      ierr = 1
      message = 'ELInletVelocityTable: empty file name'
      RETURN
    END IF
    INQUIRE(FILE=TRIM(filename), EXIST=exists)
    IF (.NOT.exists) THEN
      ierr = 1
      message = 'ELInletVelocityTable: file not found: '//TRIM(filename)
      RETURN
    END IF
    OPEN(NEWUNIT=iunit, FILE=TRIM(filename), STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios.NE.0) THEN
      ierr = 2
      message = 'ELInletVelocityTable: cannot open '//TRIM(filename)
      RETURN
    END IF
    DO
      line = ' '
      nchar = 0
      READ(iunit, '(A)', ADVANCE='NO', SIZE=nchar, IOSTAT=ios) line
      IF (IS_IOSTAT_END(ios)) THEN
        IF (nchar.EQ.0) EXIT
        ios = 0                               ! last line without a line end
      ELSE IF (IS_IOSTAT_EOR(ios)) THEN
        ios = 0
      ELSE IF (ios.EQ.0) THEN
        ierr = 8
        WRITE(message,'(A,A,A,I0,A,I0,A)') 'ELInletVelocityTable: ', TRIM(filename), ' row ', &
          nlines + 1, ' is longer than ', MAX_LINE - 1, ' characters'
        CLOSE(iunit)
        RETURN
      END IF
      IF (ios.NE.0) THEN
        ierr = 7
        WRITE(message,'(A,A,A,I0,A,I0)') 'ELInletVelocityTable: read error in ', &
          TRIM(filename), ' after row ', nlines, ', iostat ', ios
        CLOSE(iunit)
        RETURN
      END IF
      nlines = nlines + 1
      IF (nlines.GT.SIZE(lines)) THEN
        ALLOCATE(grown(2*SIZE(lines)))
        grown(1:SIZE(lines)) = lines
        CALL MOVE_ALLOC(grown, lines)
      END IF
      lines(nlines) = line
    END DO
    CLOSE(iunit)
    ierr = 0
  END SUBROUTINE READ_SNAPSHOT

  ! One line of the snapshot: data_row is .FALSE. for a blank or comment line. ierr 3: the
  ! first two tokens are not both numbers; 4: a non-finite value.
  SUBROUTINE PARSE_ROW(text, data_row, tv, uv, ierr)
    CHARACTER(LEN=*), INTENT(IN) :: text
    LOGICAL, INTENT(OUT) :: data_row
    REAL*8, INTENT(OUT) :: tv, uv
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=LEN(text)) :: line
    CHARACTER(LEN=LEN(text)) :: tok(2)
    REAL*8 :: val(2)
    INTEGER :: i, k, ios, pos, iend

    ierr = 0
    data_row = .FALSE.
    tv = 0.0d0
    uv = 0.0d0
    ! horizontal whitespace first: a tab is a blank for the blank-line and comment tests too
    ! (ADJUSTL and LEN_TRIM treat only the space as blank), the comma separator afterwards
    line = text
    DO i = 1, LEN_TRIM(line)
      IF (line(i:i).EQ.CHAR(9)) line(i:i) = ' '
    END DO
    line = ADJUSTL(line)
    IF (LEN_TRIM(line).EQ.0) RETURN
    IF (line(1:1).EQ.'#' .OR. line(1:1).EQ.'!') RETURN
    data_row = .TRUE.
    DO i = 1, LEN_TRIM(line)
      IF (line(i:i).EQ.',') line(i:i) = ' '
    END DO
    pos = 1
    DO k = 1, 2
      tok(k) = ' '
      DO WHILE (pos.LE.LEN(line))
        IF (line(pos:pos).NE.' ') EXIT
        pos = pos + 1
      END DO
      IF (pos.GT.LEN(line)) THEN
        ierr = 3
        RETURN
      END IF
      iend = INDEX(line(pos:), ' ')
      IF (iend.EQ.0) THEN
        iend = LEN(line)
      ELSE
        iend = pos + iend - 2
      END IF
      tok(k) = line(pos:iend)
      pos = iend + 1
      IF (.NOT.EL_INLET_TABLE_IS_NUMBER(tok(k))) THEN
        ierr = 3
        RETURN
      END IF
      READ(tok(k), *, IOSTAT=ios) val(k)
      IF (ios.NE.0) THEN
        ierr = 4
        RETURN
      END IF
    END DO
    IF (.NOT.(IEEE_IS_FINITE(val(1)) .AND. IEEE_IS_FINITE(val(2)))) THEN
      ierr = 4
      RETURN
    END IF
    tv = val(1)
    uv = val(2)
  END SUBROUTINE PARSE_ROW

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
