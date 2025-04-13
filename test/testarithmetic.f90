!> Test app for DNAOAD to drive Fortuno unit tests
module test_arithmetic
  use iso_fortran_env, only: real64
  use dnaoad
  use fortuno_serial, only : is_equal, test => serial_case_item, check => serial_check, test_list
  implicit none

contains

  function tests()
    type(test_list) :: tests

    tests = test_list([&
      test("add_0", add_real_0),&
      test("add_1", add_cmplx_1),&
      test("add_1", sub_real_0),&
      test("add_1", sub_cmplx_1),&
      test("div_0", div_real_0),&
      test("div_1", div_cmplx_1),&
      test("int_1", test_isint),&
      test("inv_r", test_invr),&
      test("inv_z", test_invz)&
    ])

  end function tests

  ! Test: tests 2 * x = x + x
  subroutine add_real_0()
    type(dual_real64) :: r
    call initialize_dual(r, 4)
    r%f(0) = 2.0
    r%f(1) = 1.0
    call check(2 * r == r + r)
  end subroutine add_real_0

  ! Test: tests 2 * x = x + x
  subroutine add_cmplx_1()
    type(dual_cmplx64) :: z
    call initialize_dual(z, 4)
    z%f(0) = (2.0,-3.0)
    z%f(1) = 1.0
    call check(2 * z == z + z)
  end subroutine add_cmplx_1

  ! Test: tests 2 * x - x =  x
  subroutine sub_real_0()
    type(dual_real64) :: r
    call initialize_dual(r, 4)
    r%f(0) = 2.0
    r%f(1) = 1.0
    call check(2 * r - r == r)
  end subroutine sub_real_0

  ! Test: tests 2 * x = x + x
  subroutine sub_cmplx_1()
    type(dual_cmplx64) :: z
    call initialize_dual(z, 4)
    z%f(0) = (2.0,-3.0)
    z%f(1) = 1.0
    call check(2 * z - z == z)
  end subroutine sub_cmplx_1

  ! Test: tests x**2 / x = x
  subroutine div_real_0()
    type(dual_real64) :: x
    call initialize_dual(x, 4)
    x%f(0) = 2.0
    x%f(1) = 1.0
    call check(x**2 / x == x)
  end subroutine div_real_0

  ! Test: tests x**2 / x = x
  subroutine div_cmplx_1()
    type(dual_cmplx64) :: z, fz
    call initialize_dual(z, 4)
    z%f(0) = cmplx(2.0,-3.0)
    z%f(1) = 1.0
    fz = z**2 / z - z
    call check(all( abs(fz%f) < 10.0_real64 * epsilon(0.0_real64) ))
  end subroutine div_cmplx_1

  ! Test: tests r * r^-1 = 1
  subroutine test_invr()
    type(dual_real64) :: r, fr
    call initialize_dual(r, 4)
    r%f(0) = 2.0
    r%f(1) = 1.0
    fr = r * inv(r)
    call check(all(abs(fr%f - [1,0,0,0,0]) < 10.0_real64 * epsilon(0.0_real64) ))
  end subroutine test_invr

  ! Test: tests z * z^-1 = 1
  subroutine test_invz()
    type(dual_cmplx64) :: z, fz
    call initialize_dual(z, 4)
    z%f(0) = (2.0,-1.0)
    z%f(1) = 1.0
    fz = z * inv(z)
    call check(all(abs(fz%f - [1,0,0,0,0]) < 10.0_real64 * epsilon(0.0_real64) ))
  end subroutine test_invz

  ! Test: tests if a dual number is equivalent to an integer
  subroutine test_isint()
    type(dual_real64) :: r
    type(dual_cmplx64) :: z
    call initialize_dual(r, 4)
    call initialize_dual(z, 4)
    r%f(0) = 1
    z%f(0) = 2
    call check(isInteger(r))
    call check(isInteger(z))
    r = 2
    z = 3
    call check(isInteger(r))
    call check(isInteger(z))
    r = 2
    r%f(1) = 1
    z = 3
    z%f(1) = 1
    call check(.not.isInteger(r))
    call check(.not.isInteger(z))
    z = 1
    call check(isInteger(z))
    z = cmplx(1,2,real64)
    call check(.not.isInteger(z))
  end subroutine test_isint

end module test_arithmetic


program testarithmetic
  use fortuno_serial, only : execute_serial_cmd_app
  use test_arithmetic, only : tests
  implicit none

  call execute_serial_cmd_app(tests())

end program testarithmetic
