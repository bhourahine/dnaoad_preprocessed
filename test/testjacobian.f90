!> Test app for DNAOAD to drive Fortuno unit tests
module test_jacobian
  use iso_fortran_env, only: real64
  use dnaoad
  use src_multivariablefunctions, only : jacobian, dnfvector
  use fortuno_serial, only : is_equal, test => serial_case_item, check => serial_check, test_list
  implicit none

contains

  function tests()
    type(test_list) :: tests

    tests = test_list([&
      test("test_jacobian", test_jacobian_43)&
    ])

  end function tests

  ! Test: tests jacobian times direction matches directional derivative along direction
  subroutine test_jacobian_43()
    integer, parameter :: prec = real64
    integer, parameter :: nf = 3, mq= 4
    integer :: ii
    complex(prec) :: vec(mq), q(mq), Jmat(nf, mq), JV_1(nf), JV_2(nf)
    complex(prec), parameter :: im = cmplx(0,1,kind=prec)

    vec(:) = [(1,0),(2,0),(3,0),(4,0)]
    q(:) = vec / 10.0 + im

    print*,"Jv using matmul of jacobian with direction"
    Jmat(:,:) = Jacobian(fvectest, q, nf)
    JV_1(:) = matmul(Jmat, vec)
    do ii = 1, nf
      write(*,*) JV_1(ii)
    end do
    write(*,*)

    print *, "Jv using vector directional derivative"
    JV_2(:) = f_part(dnfvector(fvectest,vec,q,nf,1), 1)
    do ii = 1, nf
      write(*,*) JV_2(ii)
    end do
    write(*,*)

    call check(all(abs(JV_1 - JV_2) < 10.0_real64 * epsilon(0.0_real64)))

  end subroutine test_jacobian_43

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

end module test_jacobian


program testmatrix
  use fortuno_serial, only : execute_serial_cmd_app
  use test_jacobian, only : tests
  implicit none

  call execute_serial_cmd_app(tests())

end program testmatrix
