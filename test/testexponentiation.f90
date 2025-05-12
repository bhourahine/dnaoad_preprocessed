!> Test app for DNAOAD to drive Fortuno unit tests
module test_exponentiation
  use iso_fortran_env, only: real64
  use dnaoad
  use fortuno_serial, only : is_equal, test => serial_case_item, check => serial_check, test_list
  implicit none

contains

  function tests()
    type(test_list) :: tests

    tests = test_list([&
      test("power_real_duals", power_real_duals),&
      test("power_real_dualint", power_real_dualint),&
      test("power_real_dualint", power_real_intdual),&
      test("power_real_zero", power_real_0),&
      test("power_cmplx_duals", power_cmplx_duals),&
      test("power_cmplx_dualint", power_cmplx_dualint),&
      test("power_cmplx_dualint", power_cmplx_intdual),&
      test("power_cmplx_zero", power_cmplx_0)&
    ])

  end function tests

  ! Tests dual raised to a dual : x^x
  subroutine power_real_duals()
    type(dual_real64) :: r, s, fr
    real(real64), parameter :: x = 2.0_real64
    real(real64), parameter :: one = 1.0_real64
    real(real64), parameter :: zero = 0.0_real64
    call initialize_dual(r, 2)
    call initialize_dual(s, 2)
    call initialize_dual(fr, 2)
    call r%set_derivatives([x,one,zero])
    s = r
    fr = r ** s
    write(*,*)r%get_derivatives()
    write(*,*)s%get_derivatives()
    write(*,*)fr%get_derivatives()
    write(*,*)"r  :", r%get_derivatives()
    write(*,*)"r^2:", fr%get_derivatives()
    write(*,*)x*x, x*x*(log(x)+one), x*x*(1/x + (one + log(x))**2)
    call check(abs(fr%get_derivative(0) - x*x) < epsilon(0.0))
    call check(abs(fr%get_derivative(1) - x*x*(log(x)+one)) < epsilon(0.0))
    call check(abs(fr%get_derivative(2) - x*x*(1/x + (one + log(x))**2)) < epsilon(0.0))
  end subroutine power_real_duals

  ! Test: tests dual raised to an integer x^2
  subroutine power_real_dualint()
    type(dual_real64) :: r, fr
    real(real64), parameter :: x = 2.0_real64
    real(real64), parameter :: one = 1.0_real64
    real(real64), parameter :: zero = 0.0_real64
    call initialize_dual(r, 2)
    call initialize_dual(fr, 2)
    call r%set_derivatives([x,one,zero])
    fr = r ** 2
    write(*,*)"r  :", r%get_derivatives()
    write(*,*)"r^2:", fr%get_derivatives()
    write(*,*)x*x, 2*x, 2
    call check(abs(fr%get_derivative(0) - x*x) < epsilon(0.0))
    call check(abs(fr%get_derivative(1) - 2.0_real64 * x) < epsilon(0.0))
    call check(abs(fr%get_derivative(2) - 2.0_real64) < epsilon(0.0))
  end subroutine power_real_dualint

  ! Test: tests integer raised to a dual 2^x
  subroutine power_real_intdual()
    type(dual_real64) :: r, fr
    real(real64), parameter :: x = 2.0_real64
    real(real64), parameter :: one = 1.0_real64
    real(real64), parameter :: zero = 0.0_real64
    call initialize_dual(r, 2)
    call initialize_dual(fr, 2)
    call r%set_derivatives([x,one,zero])
    fr = 2 ** r
    write(*,*)"r  :", r%get_derivatives()
    write(*,*)"2^r:", fr%get_derivatives()
    write(*,*)x*x, log(2.0_real64) * 2**x, (log(2.0_real64)**2) * 2**x
    call check(abs(fr%get_derivative(0) - x*x) < epsilon(0.0))
    call check(abs(fr%get_derivative(1) - log(2.0_real64) * 2**x) < epsilon(0.0))
    call check(abs(fr%get_derivative(2) - (log(2.0_real64)**2) * 2**x) < epsilon(0.0))
  end subroutine power_real_intdual

  ! Tests variants of 0 raised to or raising as a power
  subroutine power_real_0()
    type(dual_real64) :: r, fr
    integer, parameter :: order = 4
    integer :: ref(0:order)
    call initialize_dual(r, order)
    call initialize_dual(fr, order)
    ref(:) = 0
    ref(0) = 1
    fr = 2 ** r
    call check(all(abs(fr%get_derivatives() - ref) < epsilon(0.0)))
    fr = r ** 2
    call check(all(abs(fr%get_derivatives()) < epsilon(0.0)))
    fr = 0 ** r
    call check(all(abs(fr%get_derivatives() - ref) < epsilon(0.0)))
    fr = r ** 0
    call check(all(abs(fr%get_derivatives() - ref) < epsilon(0.0)))
    fr = 0.0_real64 ** r
    call check(all(abs(fr%get_derivatives() - ref) < epsilon(0.0)))
    fr = r ** 0.0_real64
    call check(all(abs(fr%get_derivatives() - ref) < epsilon(0.0)))
    fr = r ** r
    call check(all(abs(fr%get_derivatives() - ref) < epsilon(0.0)))
  end subroutine power_real_0

  ! Tests dual raised to a dual : x^x
  subroutine power_cmplx_duals()
    type(dual_cmplx64) :: r, s, fr
    complex(real64), parameter :: x = 2.0_real64
    complex(real64), parameter :: one = 1.0_real64
    complex(real64), parameter :: zero = 0.0_real64
    call initialize_dual(r, 2)
    call initialize_dual(s, 2)
    call initialize_dual(fr, 2)
    call r%set_derivatives([x,one,zero])
    s = r
    fr = r ** s
    write(*,*)r%get_derivatives()
    write(*,*)s%get_derivatives()
    write(*,*)fr%get_derivatives()
    write(*,*)"r  :", r%get_derivatives()
    write(*,*)"r^2:", fr%get_derivatives()
    write(*,*)x*x, x*x*(log(x)+one), x*x*(1/x + (one + log(x))**2)
    call check(abs(fr%get_derivative(0) - x*x) < epsilon(0.0))
    call check(abs(fr%get_derivative(1) - x*x*(log(x)+one)) < epsilon(0.0))
    call check(abs(fr%get_derivative(2) - x*x*(1/x + (one + log(x))**2)) < epsilon(0.0))
  end subroutine power_cmplx_duals

  ! Test: tests dual raised to an integer x^2
  subroutine power_cmplx_dualint()
    type(dual_cmplx64) :: r, fr
    complex(real64), parameter :: x = 2.0_real64
    complex(real64), parameter :: one = 1.0_real64
    complex(real64), parameter :: zero = 0.0_real64
    call initialize_dual(r, 2)
    call initialize_dual(fr, 2)
    call r%set_derivatives([x,one,zero])
    fr = r ** 2
    write(*,*)"r  :", r%get_derivatives()
    write(*,*)"r^2:", fr%get_derivatives()
    write(*,*)x*x, 2*x, 2
    call check(abs(fr%get_derivative(0) - x*x) < epsilon(0.0))
    call check(abs(fr%get_derivative(1) - 2.0_real64 * x) < epsilon(0.0))
    call check(abs(fr%get_derivative(2) - 2.0_real64) < epsilon(0.0))
  end subroutine power_cmplx_dualint

  ! Test: tests integer raised to a dual 2^x
  subroutine power_cmplx_intdual()
    type(dual_cmplx64) :: r, fr
    complex(real64), parameter :: x = 2.0_real64
    complex(real64), parameter :: one = 1.0_real64
    complex(real64), parameter :: zero = 0.0_real64
    call initialize_dual(r, 2)
    call initialize_dual(fr, 2)
    call r%set_derivatives([x,one,zero])
    fr = 2 ** r
    write(*,*)"r  :", r%get_derivatives()
    write(*,*)"2^r:", fr%get_derivatives()
    write(*,*)x*x, log(2.0_real64) * 2**x, (log(2.0_real64)**2) * 2**x
    call check(abs(fr%get_derivative(0) - x*x) < epsilon(0.0))
    call check(abs(fr%get_derivative(1) - log(2.0_real64) * 2**x) < epsilon(0.0))
    call check(abs(fr%get_derivative(2) - (log(2.0_real64)**2) * 2**x) < epsilon(0.0))
  end subroutine power_cmplx_intdual

  ! Tests variants of 0 raised to or raising as a power
  subroutine power_cmplx_0()
    type(dual_cmplx64) :: r, fr
    integer, parameter :: order = 4
    integer :: ref(0:order)
    call initialize_dual(r, order)
    call initialize_dual(fr, order)
    ref(:) = 0
    ref(0) = 1
    fr = 2 ** r
    call check(all(abs(fr%get_derivatives() - ref) < epsilon(0.0)))
    fr = r ** 2
    call check(all(abs(fr%get_derivatives()) < epsilon(0.0)))
    fr = 0 ** r
    call check(all(abs(fr%get_derivatives() - ref) < epsilon(0.0)))
    fr = r ** 0
    call check(all(abs(fr%get_derivatives() - ref) < epsilon(0.0)))
    fr = (0.0_real64, 0.0_real64) ** r
    call check(all(abs(fr%get_derivatives() - ref) < epsilon(0.0)))
    fr = r ** (0.0_real64, 0.0_real64)
    call check(all(abs(fr%get_derivatives() - ref) < epsilon(0.0)))
    fr = r ** r
    call check(all(abs(fr%get_derivatives() - ref) < epsilon(0.0)))
  end subroutine power_cmplx_0

end module test_exponentiation


program testexponentiation
  use fortuno_serial, only : execute_serial_cmd_app
  use test_exponentiation, only : tests
  implicit none

  call execute_serial_cmd_app(tests())

end program testexponentiation
