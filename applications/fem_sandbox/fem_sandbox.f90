program fem_sandbox
  use def_FEAT
  implicit none

  integer, parameter :: nvt_loc = 8
  integer, parameter :: net_loc = 12
  integer, parameter :: nat_loc = 6
  integer, parameter :: nel_loc = 1
  integer, parameter :: ndof_loc = 27

  real*8 :: dcorvg_local(3, nvt_loc)
  integer :: kvert_local(8, nel_loc)
  real*8 :: q2coords(3, ndof_loc)
  real*8 :: xi(3)
  real*8 :: shapeder(ndof_loc, 4, 1)
  real*8 :: jac(3, 3), jac_det
  real*8 :: invJ(3, 3)
  real*8 :: grads(ndof_loc, 3)
  integer :: i
  external :: E013A

  call initialise_globals()
  call build_hexahedron(dcorvg_local, kvert_local)
  call compute_q2_coordinates(dcorvg_local, kvert_local, q2coords)

  xi = (/0.2d0, -0.3d0, 0.4d0/)
  call evaluate_q2_at_point(xi, q2coords, shapeder, jac, jac_det, invJ, grads)

  call print_results(xi, q2coords, shapeder(:, 1, 1), jac, jac_det, grads)

contains

  subroutine initialise_globals()
    implicit none
    integer :: i

    NEL = nel_loc
    NVT = nvt_loc
    NET = net_loc
    NAT = nat_loc
    NVE = 8
    NEE = 12
    NAE = 6
    NDIM = 3
    IEL = 1
    ICHECK = 0

    do i = 1, NNDER
      BDER(i) = .false.
    end do
    BDER(1) = .true.
    BDER(2) = .true.
    BDER(3) = .true.
    BDER(4) = .true.
  end subroutine initialise_globals

  subroutine build_hexahedron(dcor, kv)
    implicit none
    real*8, intent(out) :: dcor(3, nvt_loc)
    integer, intent(out) :: kv(8, nel_loc)

    dcor(:, 1) = (/0.0d0, 0.0d0, 0.0d0/)
    dcor(:, 2) = (/2.0d0, 0.1d0, 0.0d0/)
    dcor(:, 3) = (/2.2d0, 1.0d0, 0.1d0/)
    dcor(:, 4) = (/0.1d0, 1.2d0, -0.1d0/)
    dcor(:, 5) = (/0.0d0, 0.0d0, 1.1d0/)
    dcor(:, 6) = (/2.0d0, 0.1d0, 1.2d0/)
    dcor(:, 7) = (/2.3d0, 1.1d0, 1.0d0/)
    dcor(:, 8) = (/0.2d0, 1.3d0, 1.3d0/)

    kv(:, 1) = (/1, 2, 3, 4, 5, 6, 7, 8/)
  end subroutine build_hexahedron

  subroutine compute_q2_coordinates(dcor, kv, coords)
    implicit none
    real*8, intent(in)  :: dcor(3, nvt_loc)
    integer, intent(in) :: kv(8, nel_loc)
    real*8, intent(out) :: coords(3, ndof_loc)
    integer, parameter :: neighE(2, 12) = reshape( &
         (/1,2, 2,3, 3,4, 4,1, 1,5, 2,6, 3,7, 4,8, 5,6, 6,7, 7,8, 8,5/), (/2,12/))
    integer, parameter :: neighA(4, 6) = reshape( &
         (/1,2,3,4, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 5,6,7,8/), (/4,6/))
    integer :: j
    integer :: ivt1, ivt2, ivt3, ivt4
    real*8 :: px, py, pz

    do j = 1, nvt_loc
      coords(:, j) = dcor(:, j)
    end do

    do j = 1, net_loc
      ivt1 = kv(neighE(1, j), 1)
      ivt2 = kv(neighE(2, j), 1)
      coords(:, nvt_loc + j) = 0.5d0 * (dcor(:, ivt1) + dcor(:, ivt2))
    end do

    do j = 1, nat_loc
      ivt1 = kv(neighA(1, j), 1)
      ivt2 = kv(neighA(2, j), 1)
      ivt3 = kv(neighA(3, j), 1)
      ivt4 = kv(neighA(4, j), 1)
      coords(:, nvt_loc + net_loc + j) = 0.25d0 * (dcor(:, ivt1) + dcor(:, ivt2) + dcor(:, ivt3) + dcor(:, ivt4))
    end do

    px = 0.0d0
    py = 0.0d0
    pz = 0.0d0
    do j = 1, nvt_loc
      px = px + 0.125d0 * dcor(1, kv(j, 1))
      py = py + 0.125d0 * dcor(2, kv(j, 1))
      pz = pz + 0.125d0 * dcor(3, kv(j, 1))
    end do
    coords(:, nvt_loc + net_loc + nat_loc + 1) = (/px, py, pz/)
  end subroutine compute_q2_coordinates

  subroutine evaluate_q2_at_point(refXi, coords, dhelp, jacobian, detjac, invjac, gradients)
    implicit none
    real*8, intent(in)  :: refXi(3)
    real*8, intent(in)  :: coords(3, ndof_loc)
    real*8, intent(out) :: dhelp(ndof_loc, 4, 1)
    real*8, intent(out) :: jacobian(3, 3)
    real*8, intent(out) :: detjac
    real*8, intent(out) :: invjac(3, 3)
    real*8, intent(out) :: gradients(ndof_loc, 3)
    integer :: a
    real*8 :: refgrad(3)

    jacobian = 0.0d0
    gradients = 0.0d0

    call set_derivative_flags()
    call E013A(refXi(1), refXi(2), refXi(3), dhelp, 1)

    do a = 1, ndof_loc
      jacobian(1, 1) = jacobian(1, 1) + coords(1, a) * dhelp(a, 2, 1)
      jacobian(2, 1) = jacobian(2, 1) + coords(2, a) * dhelp(a, 2, 1)
      jacobian(3, 1) = jacobian(3, 1) + coords(3, a) * dhelp(a, 2, 1)

      jacobian(1, 2) = jacobian(1, 2) + coords(1, a) * dhelp(a, 3, 1)
      jacobian(2, 2) = jacobian(2, 2) + coords(2, a) * dhelp(a, 3, 1)
      jacobian(3, 2) = jacobian(3, 2) + coords(3, a) * dhelp(a, 3, 1)

      jacobian(1, 3) = jacobian(1, 3) + coords(1, a) * dhelp(a, 4, 1)
      jacobian(2, 3) = jacobian(2, 3) + coords(2, a) * dhelp(a, 4, 1)
      jacobian(3, 3) = jacobian(3, 3) + coords(3, a) * dhelp(a, 4, 1)
    end do

    detjac = determinant3x3(jacobian)
    invjac = inverse3x3(jacobian, detjac)

    do a = 1, ndof_loc
      refgrad = (/dhelp(a, 2, 1), dhelp(a, 3, 1), dhelp(a, 4, 1)/)
      gradients(a, 1) = invjac(1, 1) * refgrad(1) + invjac(1, 2) * refgrad(2) + invjac(1, 3) * refgrad(3)
      gradients(a, 2) = invjac(2, 1) * refgrad(1) + invjac(2, 2) * refgrad(2) + invjac(2, 3) * refgrad(3)
      gradients(a, 3) = invjac(3, 1) * refgrad(1) + invjac(3, 2) * refgrad(2) + invjac(3, 3) * refgrad(3)
    end do
  end subroutine evaluate_q2_at_point

  subroutine set_derivative_flags()
    implicit none
    integer :: i
    do i = 1, NNDER
      BDER(i) = .false.
    end do
    BDER(1) = .true.
    BDER(2) = .true.
    BDER(3) = .true.
    BDER(4) = .true.
  end subroutine set_derivative_flags

  real*8 function determinant3x3(A)
    implicit none
    real*8, intent(in) :: A(3, 3)
    determinant3x3 = A(1, 1) * (A(2, 2) * A(3, 3) - A(2, 3) * A(3, 2)) &
                   - A(1, 2) * (A(2, 1) * A(3, 3) - A(2, 3) * A(3, 1)) &
                   + A(1, 3) * (A(2, 1) * A(3, 2) - A(2, 2) * A(3, 1))
  end function determinant3x3

  function inverse3x3(A, detA) result(invA)
    implicit none
    real*8, intent(in) :: A(3, 3)
    real*8, intent(in) :: detA
    real*8 :: invA(3, 3)
    real*8 :: invDet

    if (abs(detA) < 1.0d-12) then
      write(*, *) 'Jacobian determinant is near zero, mapping is degenerate.'
      stop
    end if

    invDet = 1.0d0 / detA

    invA(1, 1) =  (A(2, 2) * A(3, 3) - A(2, 3) * A(3, 2)) * invDet
    invA(1, 2) = -(A(1, 2) * A(3, 3) - A(1, 3) * A(3, 2)) * invDet
    invA(1, 3) =  (A(1, 2) * A(2, 3) - A(1, 3) * A(2, 2)) * invDet

    invA(2, 1) = -(A(2, 1) * A(3, 3) - A(2, 3) * A(3, 1)) * invDet
    invA(2, 2) =  (A(1, 1) * A(3, 3) - A(1, 3) * A(3, 1)) * invDet
    invA(2, 3) = -(A(1, 1) * A(2, 3) - A(1, 3) * A(2, 1)) * invDet

    invA(3, 1) =  (A(2, 1) * A(3, 2) - A(2, 2) * A(3, 1)) * invDet
    invA(3, 2) = -(A(1, 1) * A(3, 2) - A(1, 2) * A(3, 1)) * invDet
    invA(3, 3) =  (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) * invDet
  end function inverse3x3

  subroutine print_results(refXi, coords, shapeVals, jacobian, detjac, gradients)
    implicit none
    real*8, intent(in) :: refXi(3)
    real*8, intent(in) :: coords(3, ndof_loc)
    real*8, intent(in) :: shapeVals(ndof_loc)
    real*8, intent(in) :: jacobian(3, 3)
    real*8, intent(in) :: detjac
    real*8, intent(in) :: gradients(ndof_loc, 3)
    integer :: i

    write(*, '(A)') '--- Q2 sandbox evaluation ---'
    write(*, '(A,3F12.6)') 'Reference coordinates (xi): ', refXi
    write(*, '(A)') 'First eight Q2 node coordinates:'
    do i = 1, 8
      write(*, '(I3,3F12.6)') i, coords(1, i), coords(2, i), coords(3, i)
    end do
    write(*, '(A)') 'First eight Q2 shape values at xi:'
    do i = 1, 8
      write(*, '(I3,F18.10)') i, shapeVals(i)
    end do
    write(*, '(A)') 'Element Jacobian:'
    do i = 1, 3
      write(*, '(3F12.6)') jacobian(i, 1), jacobian(i, 2), jacobian(i, 3)
    end do
    write(*, '(A,F18.10)') 'det(J) = ', detjac
    write(*, '(A)') 'Gradients of first four basis functions in physical space:'
    do i = 1, 4
      write(*, '(I3,3F18.10)') i, gradients(i, 1), gradients(i, 2), gradients(i, 3)
    end do
  end subroutine print_results

end program fem_sandbox
