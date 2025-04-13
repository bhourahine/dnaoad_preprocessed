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

  ! Test: tests dual raised to a dual
  subroutine power_real_duals()
    type(dual_real64) :: r, s, fr
    call initialize_dual(r, 4)
    call initialize_dual(s, 4)
    call initialize_dual(fr, 4)
    r%f(0) = 2.0
    r%f(1) = 1.0
    s%f(0) = 2.0
    s%f(1) = 1.0
    fr = r ** s
    write(*,*)r%f
    write(*,*)s%f
    write(*,*)fr%f
    call check(abs(fr%f(0) - 4.0) < epsilon(0.0))
  end subroutine power_real_duals

  ! Test: tests dual raised to an integer
  subroutine power_real_dualint()
    type(dual_real64) :: r, fr
    call initialize_dual(r, 4)
    call initialize_dual(fr, 4)
    r%f(0) = 2.0
    r%f(1) = 1.0
    fr = r ** 2
    write(*,*)r%f
    write(*,*)fr%f
    call check(abs(fr%f(0) - 4.0) < epsilon(0.0))
  end subroutine power_real_dualint

  ! Test: tests integer raised to a dual
  subroutine power_real_intdual()
    type(dual_real64) :: r, fr
    call initialize_dual(r, 4)
    call initialize_dual(fr, 4)
    r%f(0) = 2.0
    r%f(1) = 1.0
    fr = 2 ** r
    write(*,*)r%f
    write(*,*)fr%f
    call check(abs(fr%f(0) - 4.0) < epsilon(0.0))
  end subroutine power_real_intdual

  ! Test: tests integer raised to a dual
  subroutine power_real_0()
    type(dual_real64) :: r, fr
    call initialize_dual(r, 4)
    call initialize_dual(fr, 4)
    fr = 2 ** r
    call check(abs(fr%f(0) - 1) < epsilon(0.0))
    fr = r ** 2
    call check(abs(fr%f(0)) < epsilon(0.0))
    fr = 0 ** r
    call check(abs(fr%f(0) - 1) < epsilon(0.0))
    fr = r ** 0
    call check(abs(fr%f(0) - 1) < epsilon(0.0))
    fr = 0.0_real64 ** r
    call check(abs(fr%f(0) - 1) < epsilon(0.0))
    fr = r ** 0.0_real64
    call check(abs(fr%f(0) - 1) < epsilon(0.0))
    fr = r ** r
    call check(abs(fr%f(0) - 1) < epsilon(0.0))
  end subroutine power_real_0

  ! Test: tests dual raised to a dual
  subroutine power_cmplx_duals()
    type(dual_cmplx64) :: r, s, fr
    call initialize_dual(r, 4)
    call initialize_dual(s, 4)
    call initialize_dual(fr, 4)
    r%f(0) = 2.0
    r%f(1) = 1.0
    s%f(0) = 2.0
    s%f(1) = 1.0
    fr = r ** s
    write(*,*)r%f
    write(*,*)s%f
    write(*,*)fr%f
    call check(abs(fr%f(0) - 4.0) < epsilon(0.0))
  end subroutine power_cmplx_duals

  ! Test: tests dual raised to an integer
  subroutine power_cmplx_dualint()
    type(dual_cmplx64) :: r, fr
    call initialize_dual(r, 4)
    call initialize_dual(fr, 4)
    r%f(0) = 2.0
    r%f(1) = 1.0
    fr = r ** 2
    write(*,*)r%f
    write(*,*)fr%f
    call check(abs(fr%f(0) - 4.0) < epsilon(0.0))
  end subroutine power_cmplx_dualint

  ! Test: tests integer raised to a dual
  subroutine power_cmplx_intdual()
    type(dual_cmplx64) :: r, fr
    call initialize_dual(r, 4)
    call initialize_dual(fr, 4)
    r%f(0) = 2.0
    r%f(1) = 1.0
    fr = 2 ** r
    write(*,*)r%f
    write(*,*)fr%f
    call check(abs(fr%f(0) - 4.0) < epsilon(0.0))
  end subroutine power_cmplx_intdual

  ! Test: tests integer raised to a dual
  subroutine power_cmplx_0()
    type(dual_cmplx64) :: r, fr
    call initialize_dual(r, 4)
    call initialize_dual(fr, 4)
    fr = 2 ** r
    call check(abs(fr%f(0) - 1) < epsilon(0.0))
    fr = r ** 2
    call check(abs(fr%f(0)) < epsilon(0.0))
    fr = 0 ** r
    call check(abs(fr%f(0) - 1) < epsilon(0.0))
    fr = r ** 0
    call check(abs(fr%f(0) - 1) < epsilon(0.0))
    fr = cmplx(0,0,kind=real64) ** r
    call check(abs(fr%f(0) - 1) < epsilon(0.0))
    fr = r ** cmplx(0,0,kind=real64)
    call check(abs(fr%f(0) - 1) < epsilon(0.0))
    fr = 0_real64 ** r
    call check(abs(fr%f(0) - 1) < epsilon(0.0))
    fr = r ** 0_real64
    call check(abs(fr%f(0) - 1) < epsilon(0.0))
    fr = r ** r
    call check(abs(fr%f(0) - 1) < epsilon(0.0))
  end subroutine power_cmplx_0

end module test_exponentiation


program testexponentiation
  use fortuno_serial, only : execute_serial_cmd_app
  use test_exponentiation, only : tests
  implicit none

  call execute_serial_cmd_app(tests())

end program testexponentiation
