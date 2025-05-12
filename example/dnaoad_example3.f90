!> Example #2
!!
!! Program demonstrating capabilities of the dnaoad library.
!!
program dnaoad_example2
  use iso_fortran_env, only : real64
  use dnaoad
  use src_multivariablefunctions, only : jacobian, dnfvector, hessian, gradient
  implicit none

  integer, parameter :: prec = real64
  integer, parameter :: order = 5
  integer, parameter :: nf = 3, mq= 4
  integer :: ii
  complex(prec) :: vec(mq), q(mq), Jmat(nf, mq), JV(nf), GV(3), Hmat(3,3)
  complex(prec), parameter :: im = cmplx(0,1,kind=prec)

  vec(:) = [(1,0),(2,0),(3,0),(4,0)]
  q(:) = vec / 10.0 + im

  print*,"Jv using matmul"
  Jmat = Jacobian(fvectest, q, nf)
  JV = matmul(Jmat, vec)
  do ii = 1, nf
    write(*,*) JV(ii)
  end do
  write(*,*)

  print *, "Jv using vector directional derivative"
  JV = f_part(dnfvector(fvectest,vec,q,nf,1), 1)
  do ii = 1, nf
    write(*,*) JV(ii)
  end do
  write(*,*)

  print *, "---Hessian matrix---"
  Hmat(:,:) = Hessian(fstest,q(1:3))
  do ii = 1, 3
    write(*,"(A,i0)") "row:",ii
    write(*,*) Hmat(:, ii)
  end do
  write(*,*)

  print *, "---Gradient---"
  GV(:) = gradient(fstest,q(1:3))
  do ii = 1, 3
    write(*,*) GV(ii)
  end do

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
    integer :: order

    order = r(1)%order()
    call initialize_dual(f1, order)
    call initialize_dual(f2, order)
    call initialize_dual(f3, order)
    f1 = sin(product(r))
    f2 = cos(product(r))*sqrt(r(4)/r(2) - r(1)/r(3))
    f3 = sin(log(product(r)))

    allocate(fr(3))
    call initialize_dual(fr, order)
    fr(:) = [f1,f2,f3]

  end function fvectest

end program dnaoad_example2
