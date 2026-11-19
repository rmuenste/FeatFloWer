!=========================================================================
! test_el_inlet_table - inlet velocity table of boundary type Inflow975 (EL_INLET_TABLE,
! EL_CONFIG:EL_INLET_VELOCITY; SimPar@ELInletVelocityTable). Case 11; docs/md_docs/
! el_min_fluidization_prep_2026-11-19.md.
!
! Checks (exit code = the failing check):
!   10-16  table reading: comments, blank lines, commas, quotes in the deck value
!   20-27  interpolation: rows exact, midpoints, constant beyond both ends, one-row table
!   30-35  refusals (ierr, no stop): missing file, malformed row, non-finite value, equal and
!          decreasing times, empty file
!   40-42  EL_INLET_VELOCITY: ramp unchanged without a table (same expression), table with one
!=========================================================================
PROGRAM TEST_EL_INLET_TABLE

  USE EL_INLET_TABLE, ONLY: tElInletTable, EL_INLET_TABLE_READ, EL_INLET_TABLE_VALUE, &
                            EL_INLET_TABLE_PATH
  USE EL_CONFIG, ONLY: el_inlet_u0, el_inlet_u1, el_inlet_tramp, el_inlet_table, &
                       el_inlet_table_file, EL_INLET_VELOCITY, EL_INLET_CONFIGURED

  IMPLICIT NONE

  TYPE(tElInletTable) :: tab
  INTEGER :: ierr, i
  CHARACTER(LEN=320) :: msg
  REAL*8 :: v

  ! ---- 10-16: reading --------------------------------------------------------------------
  CALL WRITE_FILE('tab_ok.txt', [CHARACTER(LEN=40) :: '# t U (protocol)', '', &
    '0.0  0.30', '0.5, 0.30   ! end of preparation', '  0.52 0.28', '', '1.0 0.28'])
  CALL EL_INLET_TABLE_READ('tab_ok.txt', tab, ierr, msg)
  IF (ierr.NE.0) THEN
    WRITE(*,*) TRIM(msg)
    STOP 10
  END IF
  IF (tab%n.NE.4) STOP 11
  IF (ANY(tab%t.NE.(/0.0d0, 0.5d0, 0.52d0, 1.0d0/))) STOP 12
  IF (ANY(tab%u.NE.(/0.30d0, 0.30d0, 0.28d0, 0.28d0/))) STOP 13
  IF (TRIM(EL_INLET_TABLE_PATH(' "runs/x/inlet table.txt" ! c')).NE.'runs/x/inlet table.txt') STOP 14
  IF (TRIM(EL_INLET_TABLE_PATH(" 'a/b.txt'")).NE.'a/b.txt') STOP 15
  IF (TRIM(EL_INLET_TABLE_PATH('  _data/inlet.txt  ')).NE.'_data/inlet.txt') STOP 16

  ! ---- 20-27: interpolation ----------------------------------------------------------------
  IF (EL_INLET_TABLE_VALUE(tab, 0.5d0).NE.0.30d0) STOP 20
  IF (EL_INLET_TABLE_VALUE(tab, 0.52d0).NE.0.28d0) STOP 21
  IF (ABS(EL_INLET_TABLE_VALUE(tab, 0.51d0) - 0.29d0).GT.1.0d-15) STOP 22
  IF (ABS(EL_INLET_TABLE_VALUE(tab, 0.505d0) - 0.295d0).GT.1.0d-15) STOP 23
  IF (EL_INLET_TABLE_VALUE(tab, -1.0d0).NE.0.30d0) STOP 24
  IF (EL_INLET_TABLE_VALUE(tab, 7.5d0).NE.0.28d0) STOP 25
  IF (EL_INLET_TABLE_VALUE(tab, 0.25d0).NE.0.30d0) STOP 26
  CALL WRITE_FILE('tab_one.txt', [CHARACTER(LEN=40) :: '0.1 0.2'])
  CALL EL_INLET_TABLE_READ('tab_one.txt', tab, ierr, msg)
  IF (ierr.NE.0 .OR. tab%n.NE.1) STOP 27
  IF (EL_INLET_TABLE_VALUE(tab, 0.0d0).NE.0.2d0 .OR. EL_INLET_TABLE_VALUE(tab, 9.0d0).NE.0.2d0) STOP 27

  ! ---- 30-35: refusals -----------------------------------------------------------------------
  CALL EL_INLET_TABLE_READ('no_such_table.txt', tab, ierr, msg)
  CALL EXPECT_REFUSAL(ierr, tab, msg, 30)
  CALL WRITE_FILE('tab_bad.txt', [CHARACTER(LEN=40) :: '0.0 0.3', '0.5'])
  CALL EL_INLET_TABLE_READ('tab_bad.txt', tab, ierr, msg)
  CALL EXPECT_REFUSAL(ierr, tab, msg, 31)
  CALL WRITE_FILE('tab_nan.txt', [CHARACTER(LEN=40) :: '0.0 0.3', '0.5 NaN'])
  CALL EL_INLET_TABLE_READ('tab_nan.txt', tab, ierr, msg)
  CALL EXPECT_REFUSAL(ierr, tab, msg, 32)
  CALL WRITE_FILE('tab_eq.txt', [CHARACTER(LEN=40) :: '0.0 0.3', '0.5 0.3', '0.5 0.28'])
  CALL EL_INLET_TABLE_READ('tab_eq.txt', tab, ierr, msg)
  CALL EXPECT_REFUSAL(ierr, tab, msg, 33)
  CALL WRITE_FILE('tab_dec.txt', [CHARACTER(LEN=40) :: '0.0 0.3', '0.5 0.3', '0.4 0.28'])
  CALL EL_INLET_TABLE_READ('tab_dec.txt', tab, ierr, msg)
  CALL EXPECT_REFUSAL(ierr, tab, msg, 34)
  CALL WRITE_FILE('tab_empty.txt', [CHARACTER(LEN=40) :: '# only a comment', ''])
  CALL EL_INLET_TABLE_READ('tab_empty.txt', tab, ierr, msg)
  CALL EXPECT_REFUSAL(ierr, tab, msg, 35)

  ! ---- 40-42: EL_INLET_VELOCITY -------------------------------------------------------------
  el_inlet_u0 = 0.002d0
  el_inlet_u1 = 0.02d0
  el_inlet_tramp = 0.4d0
  IF (.NOT.EL_INLET_CONFIGURED()) STOP 40
  DO i = 0, 10
    v = 0.05d0*DBLE(i)
    IF (EL_INLET_VELOCITY(v).NE.el_inlet_u0 + (el_inlet_u1 - el_inlet_u0)* &
        MIN(1.0d0, MAX(0.0d0, v)/el_inlet_tramp)) STOP 41
  END DO
  el_inlet_u0 = 0.0d0
  el_inlet_u1 = 0.0d0
  el_inlet_table_file = 'tab_ok.txt'
  IF (.NOT.EL_INLET_CONFIGURED()) STOP 42
  CALL EL_INLET_TABLE_READ(el_inlet_table_file, el_inlet_table, ierr, msg)
  IF (ierr.NE.0) STOP 42
  IF (ABS(EL_INLET_VELOCITY(0.51d0) - 0.29d0).GT.1.0d-15) STOP 42

  WRITE(*,'(A)') 'test_el_inlet_table: all checks passed'

CONTAINS

  SUBROUTINE WRITE_FILE(name, lines)
    CHARACTER(LEN=*), INTENT(IN) :: name
    CHARACTER(LEN=*), INTENT(IN) :: lines(:)
    INTEGER :: u, k
    OPEN(NEWUNIT=u, FILE=name, STATUS='REPLACE', ACTION='WRITE')
    DO k = 1, SIZE(lines)
      WRITE(u, '(A)') TRIM(lines(k))
    END DO
    CLOSE(u)
  END SUBROUTINE WRITE_FILE

  SUBROUTINE EXPECT_REFUSAL(ierr_in, table, message, code)
    INTEGER, INTENT(IN) :: ierr_in, code
    TYPE(tElInletTable), INTENT(IN) :: table
    CHARACTER(LEN=*), INTENT(IN) :: message
    WRITE(*,'(A,I0,A,A)') 'refusal ', code, ': ', TRIM(message)
    IF (ierr_in.EQ.0 .OR. table%n.NE.0 .OR. LEN_TRIM(message).EQ.0) STOP code
  END SUBROUTINE EXPECT_REFUSAL

END PROGRAM TEST_EL_INLET_TABLE
