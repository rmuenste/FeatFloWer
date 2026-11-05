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

  PURE REAL*8 FUNCTION EL_KERNEL_DERIVATIVE(distance, width)

    ! dW/dr of EL_KERNEL_VALUE (same prefactor conventions).
    REAL*8, INTENT(IN) :: distance, width
    REAL*8 :: q

    EL_KERNEL_DERIVATIVE = 0.0d0
    IF (width.LE.0.0d0 .OR. distance.LT.0.0d0) RETURN

    q = distance / width
    IF (q.GE.1.0d0) RETURN

    SELECT CASE (TRIM(el_kernel))
    CASE ('gaussian')
      EL_KERNEL_DERIVATIVE = -9.0d0*q*EXP(-4.5d0*q*q) / width**4
    CASE DEFAULT
      ! d/dq [(1+3q)(1-q)^3] = -12 q (1-q)^2
      EL_KERNEL_DERIVATIVE = -12.0d0*q*(1.0d0 - q)**2 / width**4
    END SELECT

  END FUNCTION EL_KERNEL_DERIVATIVE

END MODULE EL_KERNEL_FUNCTIONS
