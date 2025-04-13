!> Example #2
!!
!! Program demonstrating capabilities of the dnaoad library.
!!
program dnaoad_example2
  use iso_fortran_env, only : real64
  use dnaoad
  use src_multivariablefunctions
  implicit none

  integer, parameter :: prec = real64
  integer, parameter :: order = 5
  complex(prec) :: v(3), q(3), fr(3), v5D(5)
  complex(prec), allocatable :: H(:,:)
  type(dual_cmplx64) :: fval
  real :: t1, t2

  ! since a dualzn numbers is an allocatable entity, do not forget to
  ! initialize it

  v(:) = [(1,0),(2,0),(3,0)]
  q(:) = v / 10.0 + (0,1)

  call cpu_time(t1)
  fval = dnfscalar_dir(fstest, v, q, 1)
  call cpu_time(t2)

  write(*,*)fval%f

  print*,"elapsed time (s):",t2-t1

  call cpu_time(t1)
  fr(:) = gradient(fstest, q)
  call cpu_time(t2)

  write(*,*)fr

  print*,"elapsed time (s):",t2-t1

  v5D(:) = (1,.5)
  H = hessian(fstest, v5D)
  block
    integer :: ii
    do ii = 1, size(H,dim=2)
      write(*,*)H(:,ii)
    end do
  end block

contains

  !> Example of scalar function f = sin(x*targety*z) + cos(x*y*z)
  function fstest(r) result(fr)

    type(dual_cmplx64), intent(in) :: r(:)
    type(dual_cmplx64) :: fr

    fr = sin(product(r)) + cos(product(r))

  end function fstest

  !> Example of vector function f = [f1,f2,f3]
  !! f = fvectest(r) is a function f:D^m ---> Dn
  function fvectest(r) result(fr)

    type(dual_cmplx64), intent(in) :: r(:)
    type(dual_cmplx64), allocatable :: fr(:)

    type(dual_cmplx64) :: f1,f2,f3

    f1 = sin(product(r))
    f2 = cos(product(r))*sqrt(r(4)/r(2) - r(1)/r(3))
    f3 = sin(log(product(r)))

    allocate(fr(3))
    fr = [f1,f2,f3]

  end function fvectest

end program dnaoad_example2
