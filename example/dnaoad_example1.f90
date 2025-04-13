!> Example #1
!!
!! Program demonstrating capabilities of the dnaoad library.
!!
program dnaoad_example1
  use iso_fortran_env, only : real64
  use dnaoad
  implicit none

  integer, parameter :: prec = real64
  integer, parameter :: order = 5
  integer :: kk
  real(prec) :: r0
  complex(prec) :: z0
  type(dual_cmplx64) :: z, fz
  type(dual_real64) :: r, fr

  z0 = (1.1_prec,2.2_prec)
  r0 = 1.1_prec

  ! We initialize the dual number to 0 and up to the 5th derivative
  call initialize_dual(z, order)

  !we set the 0-th and 1-th components. If dual numbers are used to
  !calculate D^n f(z0) then r must be of the form r = r0 + 1*eps_1
  z%f(0) = z0
  z%f(1) = 1

  fz = sin(z)**log(z*z)

  !writing the derivatives, from the 0th derivative up to the
  !order-th derivative.
  print*,"Complex derivatives"
  do kk = 0, order
    print*, fz%f(kk)
  end do

  ! We initialize the dual number to 0 and up to the 5th derivative
  call initialize_dual(r, order)
  r%f(0) = r0
  r%f(1) = 1

  fr = sin(r)**log(r*r)

  print*,"real derivatives"
  do kk = 0, order
    print*, fr%f(kk)
  end do

end program dnaoad_example1
