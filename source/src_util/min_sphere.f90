subroutine min_sphere(dx, n, center, radius)
  implicit none
  integer, intent(in) :: n
  real(kind=8), dimension(3,n), intent(in) :: dx
  real(kind=8), dimension(3), intent(out) :: center
  real(kind=8), intent(out) :: radius
  real(kind=8), dimension(3) :: q
  real(kind=8), dimension(3,3) :: b, z
  real(kind=8), dimension(n) :: d, u
  real(kind=8) :: alpha, beta, gamma, t, delta
  integer :: i, j, k, l, m
  logical :: modified

  ! Initialization
  center = 0.0
  do i = 1, n
    center = center + dx(:,i)
  end do
  center = center / n
  
  d = 0.0
  do i = 1, n
    d(i) = sqrt(sum((dx(:,i) - center)**2))
  end do
  
  radius = maxval(d)
  q = dx(:,1)

  do i = 1, 1000  ! maximum number of iterations
    modified = .false.

    do j = 1, n
      if (d(j) > radius * (1.0 + 1.0e-12)) then  ! if point j is outside the sphere
        t = sqrt(sum((dx(:,j) - q)**2))
        if (t > 1.0e-12) then  ! if the point is not coincident with q
          u = (dx(:,j) - q) / t
          alpha = dot_product(u, (dx(:,j) - center))**2
          beta = radius**2 - t**2
          gamma = dot_product((dx(:,j) - center), q)**2
          delta = beta * gamma - alpha * (gamma - beta)
          if (delta >= 0.0) then  ! if the discriminant is non-negative
            t = beta * dot_product(u, (center - dx(:,j))) / delta
            q = q + t * ((center - dx(:,j)) + (q - dx(:,j)) * gamma / beta)
            radius = sqrt(sum((q - center)**2))
            modified = .true.
          end if
        end if
      end if
    end do

    if (.not. modified) exit  ! exit if no point has been moved

  end do
  
  DO j = 1,n
   q = dx(:,j)
   if (radius < sqrt(sum((q - center)**2))) then
    radius = sqrt(sum((q - center)**2))
   end if
  END DO

end subroutine min_sphere



