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
      test("sub_0", sub_real_0),&
      test("sub_1", sub_cmplx_1),&
      test("div_0", div_real_0),&
      test("div_1", div_cmplx_1),&
      test("inv_r", test_invx),&
      test("inv_z", test_invz),&
      test("int_1", test_isint)&
    ])

  end function tests

  ! Test: tests 2 * x = x + x
  subroutine add_real_0()
    type(dual_real64) :: x
    call initialize_dual(x, 2)
    call x%set_derivatives([2.0_real64, 1.0_real64, 0.0_real64])
    call check(2 * x == x + x)
  end subroutine add_real_0

  ! Test: tests 2 * z = z + z
  subroutine add_cmplx_1()
    type(dual_cmplx64) :: z
    call initialize_dual(z, 2)
    call z%set_derivatives([cmplx(2.0,-3.0, real64), cmplx(1.0,0.0, real64), cmplx(0,0, real64)])
    call check(2 * z == z + z)
  end subroutine add_cmplx_1

  ! Test: tests 2 * x - x =  x
  subroutine sub_real_0()
    type(dual_real64) :: x
    call initialize_dual(x, 2)
    x = [2,1]
    call check(2 * x - x == x)
  end subroutine sub_real_0

  ! Test: tests 2 * z - z =  z
  subroutine sub_cmplx_1()
    type(dual_cmplx64) :: z
    call initialize_dual(z, 2)
    call z%set_derivatives([cmplx(2.0,-3.0, real64), cmplx(1.0,0.0, real64), cmplx(0,0, real64)])
    call check(2 * z - z == z)
  end subroutine sub_cmplx_1

  ! Test: tests x**2 / x = x
  subroutine div_real_0()
    type(dual_real64) :: x
    call initialize_dual(x, 4)
    x = [2,1]
    call check(x**2 / x == x)
  end subroutine div_real_0

  ! Test: tests z**2 / z = z
  subroutine div_cmplx_1()
    type(dual_cmplx64) :: z, fz
    call initialize_dual(z, 4)
    z = [cmplx(2.0,-3.0, real64), cmplx(1.0,0.0, real64)]
    fz = z**2 / z - z
    call check(all( abs(fz%get_derivatives()) < 10.0_real64 * epsilon(0.0_real64) ))
  end subroutine div_cmplx_1

  ! Test: tests x * x^-1 = 1
  subroutine test_invx()
    type(dual_real64) :: x, fx
    complex(real64), parameter :: one = 1.0_real64
    complex(real64), parameter :: zero = 0.0_real64
    call initialize_dual(x, 4)
    x = [2,1]
    fx = x * inv(x)
    call check(all(abs(fx%get_derivatives() - [one,zero,zero,zero,zero])&
        & < 10.0_real64 * epsilon(0.0_real64) ))
  end subroutine test_invx

  ! Test: tests z * z^-1 = 1
  subroutine test_invz()
    type(dual_cmplx64) :: z, fz
    call initialize_dual(z, 4)
    z = [cmplx(2.0,-1.0, real64), cmplx(1.0,0.0, real64)]
    fz = z * inv(z)
    call check(all(abs(fz%get_derivatives() - [1,0,0,0,0]) < 10.0_real64 * epsilon(0.0_real64) ))
  end subroutine test_invz

  ! Test: tests if a dual number is equivalent to an integer
  subroutine test_isint()
    type(dual_real64) :: x
    type(dual_cmplx64) :: z
    call initialize_dual(x, 4)
    call initialize_dual(z, 4)
    x = 1
    call check(isInteger(x))
    z = 2
    call check(isInteger(z))
    x = 1.2_real64
    call check(.not.isInteger(x))
    z = (2.2_real64,0.0_real64)
    call check(.not.isInteger(z))
    x = [2,1]
    call check(.not.isInteger(x))
    z = [3,1]
    call check(.not.isInteger(z))
    x = 1
    call check(isInteger(x))
    z = 2
    call check(isInteger(z))
    z = cmplx(1,0,real64)
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
