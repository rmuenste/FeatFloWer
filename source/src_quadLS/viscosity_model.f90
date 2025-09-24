module viscosity_model
  implicit none
  private
  public :: mu_eff, set_mu_f

  ! Runtime-settable viscosity
  real*8 :: mu_f = 1.0      ! Default; to be set from input
  real*8, parameter :: phi_max   = 0.64
  real*8, parameter :: kd_exp    = 0.605
  real*8, parameter :: phi_low   = 0.05
  real*8, parameter :: gdot_min  = 0.001
  real*8, parameter :: gdot_max  = 1000.0

contains

  subroutine set_mu_f(mu_input)
    real*8, intent(in) :: mu_input
    mu_f = mu_input
  end subroutine set_mu_f

  pure function mu_low(phi) result(mu)
    real*8, intent(in) :: phi
    real*8 :: mu
    mu = mu_f * (1.0 + 2.5 * phi)
  end function mu_low

  pure function mu_fit(phi) result(mu)
    real*8, intent(in) :: phi
    real*8 :: mu, denom
    denom = 1.0 - phi / phi_max
    if (denom <= 0.0) then
      mu = 1.0e10
    else
      mu = mu_f * denom ** (-2.5 * kd_exp)
    end if
  end function mu_fit

  pure function mu_shear(phi, gdot) result(mu)
    real*8, intent(in) :: phi, gdot
    real*8 :: mu, beta, g_ref
    beta = 0.10
    g_ref = 1.0
    mu = mu_fit(phi) * (1.0 + beta * log10(gdot / g_ref))
  end function mu_shear

  pure function mu_eff(phi, gdot) result(mu)
    real*8, intent(in) :: phi, gdot
    real*8 :: mu
    if (phi < phi_low) then
      mu = mu_low(phi)
    else if (gdot >= gdot_min .and. gdot <= gdot_max) then
      mu = mu_shear(phi, gdot)
    else
      mu = mu_fit(phi)
    end if
  end function mu_eff

end module viscosity_model
