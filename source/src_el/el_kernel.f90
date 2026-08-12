MODULE EL_KERNEL_FUNCTIONS

  USE EL_CONFIG, ONLY: el_kernel

  IMPLICIT NONE

CONTAINS

  PURE REAL*8 FUNCTION EL_KERNEL_VALUE(distance, width)

    REAL*8, INTENT(IN) :: distance, width
    REAL*8 :: q

    EL_KERNEL_VALUE = 0.0d0
    IF (width.LE.0.0d0 .OR. distance.LT.0.0d0) RETURN

    q = distance / width
    IF (q.GE.1.0d0) RETURN

    SELECT CASE (TRIM(el_kernel))
    CASE ('gaussian')
      EL_KERNEL_VALUE = EXP(-4.5d0*q*q) / width**3
    CASE DEFAULT
      ! Positive compact Lucy polynomial. Discrete transfer normalizes the
      ! weights, so the continuous prefactor is intentionally omitted.
      EL_KERNEL_VALUE = (1.0d0 + 3.0d0*q) * (1.0d0 - q)**3 / width**3
    END SELECT

  END FUNCTION EL_KERNEL_VALUE

END MODULE EL_KERNEL_FUNCTIONS
