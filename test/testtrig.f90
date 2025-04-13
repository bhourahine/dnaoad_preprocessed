!> Test app for DNAOAD to drive Fortuno unit tests
module test_trig
  use iso_fortran_env, only: real64
  use dnaoad
  use fortuno_serial, only : is_equal, test => serial_case_item, check => serial_check, test_list
  implicit none

contains

  function tests()
    type(test_list) :: tests

    tests = test_list([&
        test("sin(x)", test_sin_real),&
        test("cos(x)", test_cos_real),&
        test("tan(x)", test_tan_real),&
        test("sin(z)", test_sin_cmplx),&
        test("cos(z)", test_cos_cmplx),&
        test("tan(z)", test_tan_cmplx),&
        test("sinh(x)", test_sinh_real),&
        test("cosh(x)", test_cosh_real),&
        test("tanh(x)", test_tanh_real),&
        test("sinh(z)", test_sinh_cmplx),&
        test("cosh(z)", test_cosh_cmplx),&
        test("tanh(z)", test_tanh_cmplx)&
    ])

  end function tests

  ! Test: sin and asin as inverses
  subroutine test_sin_real()
    type(dual_real64) :: fr, r
    call initialize_dual(r, 4)
    call initialize_dual(fr, 4)
    r%f(0) = 0.85
    r%f(1) = 1.0
    fr = sin(r)
    write(*,*)"sin()", fr%f
    fr = asin(fr)
    write(*,*)"asin(sin())", fr%f
    fr = fr - r
    write(*,*)"difference", fr%f
    call check(all( abs(fr%f) < 100.0_real64 * epsilon(0.0_real64) ))
  end subroutine test_sin_real

  ! Test: cos and acos as inverses
  subroutine test_cos_real()
    type(dual_real64) :: fr, r
    call initialize_dual(r, 4)
    call initialize_dual(fr, 4)
    r%f(0) = 0.85
    r%f(1) = 1.0
    fr = cos(r)
    write(*,*)"cos()", fr%f
    fr = acos(fr)
    write(*,*)"acos(cos())", fr%f
    fr = fr - r
    write(*,*)"difference", fr%f
    call check(all( abs(fr%f) < 100.0_real64 * epsilon(0.0_real64) ))
  end subroutine test_cos_real

  ! Test: tan and atan as inverses
  subroutine test_tan_real()
    type(dual_real64) :: fr, r
    call initialize_dual(r, 4)
    call initialize_dual(fr, 4)
    r%f(0) = 0.85
    r%f(1) = 1.0
    fr = tan(r)
    write(*,*)"tan()", fr%f
    fr = atan(fr)
    write(*,*)"atan(tan())", fr%f
    fr = fr - r
    write(*,*)"difference", fr%f
    call check(all( abs(fr%f) < 100.0_real64 * epsilon(0.0_real64) ))
  end subroutine test_tan_real

  ! Test: sin and asin as inverses
  subroutine test_sin_cmplx()
    type(dual_cmplx64) :: fz, z
    call initialize_dual(z, 4)
    call initialize_dual(fz, 4)
    z%f(0) = (0.85,0.15)
    z%f(1) = 1.0
    fz = sin(z)
    write(*,*)"sin()", fz%f
    fz = asin(fz)
    write(*,*)"asin(sin())", fz%f
    fz = fz - z
    write(*,*)"difference", fz%f
    call check(all( abs(fz%f) < 100.0_real64 * epsilon(0.0_real64) ))
  end subroutine test_sin_cmplx

  ! Test: cos and acos as inverses
  subroutine test_cos_cmplx()
    type(dual_cmplx64) :: fz, z
    call initialize_dual(z, 4)
    call initialize_dual(fz, 4)
    z%f(0) = (0.85,0.15)
    z%f(1) = 1.0
    fz = cos(z)
    write(*,*)"cos()", fz%f
    fz = acos(fz)
    write(*,*)"acos(cos())", fz%f
    fz = fz - z
    write(*,*)"difference", fz%f
    call check(all( abs(fz%f) < 100.0_real64 * epsilon(0.0_real64) ))
  end subroutine test_cos_cmplx

  ! Test: tan and atan as inverses
  subroutine test_tan_cmplx()
    type(dual_cmplx64) :: fz, z
    call initialize_dual(z, 4)
    call initialize_dual(fz, 4)
    z%f(0) = (0.85,0.15)
    z%f(1) = 1.0
    fz = tan(z)
    write(*,*)"tan()", fz%f
    fz = atan(fz)
    write(*,*)"atan(tan())", fz%f
    fz = fz - z
    write(*,*)"difference", fz%f
    call check(all( abs(fz%f) < 1000.0_real64 * epsilon(0.0_real64) ))
  end subroutine test_tan_cmplx

  ! Hyperbolic functions

  ! Test: sinh and asinh as inverses
  subroutine test_sinh_real()
    type(dual_real64) :: fr, r
    call initialize_dual(r, 4)
    call initialize_dual(fr, 4)
    r%f(0) = 0.85
    r%f(1) = 1.0
    fr = sinh(r)
    write(*,*)"sinh()", fr%f
    fr = asinh(fr)
    write(*,*)"asinh(sinh())", fr%f
    fr = fr - r
    write(*,*)"difference", fr%f
    call check(all( abs(fr%f) < 100.0_real64 * epsilon(0.0_real64) ))
  end subroutine test_sinh_real

  ! Test: cosh and acosh as inverses
  subroutine test_cosh_real()
    type(dual_real64) :: fr, r
    call initialize_dual(r, 4)
    call initialize_dual(fr, 4)
    r%f(0) = 0.85
    r%f(1) = 1.0
    fr = cosh(r)
    write(*,*)"cosh()", fr%f
    fr = acosh(fr)
    write(*,*)"acosh(cosh())", fr%f
    fr = fr - r
    write(*,*)"difference", fr%f
    call check(all( abs(fr%f) < 100.0_real64 * epsilon(0.0_real64) ))
  end subroutine test_cosh_real

  ! Test: tanh and atanh as inverses
  subroutine test_tanh_real()
    type(dual_real64) :: fr, r
    call initialize_dual(r, 4)
    call initialize_dual(fr, 4)
    r%f(0) = 0.85
    r%f(1) = 1.0
    fr = tanh(r)
    write(*,*)"tanh()", fr%f
    fr = atanh(fr)
    write(*,*)"atanh(tanh())", fr%f
    fr = fr - r
    write(*,*)"difference", fr%f
    call check(all( abs(fr%f) < 100.0_real64 * epsilon(0.0_real64) ))
  end subroutine test_tanh_real

  ! Test: sinh and asinh as inverses
  subroutine test_sinh_cmplx()
    type(dual_cmplx64) :: fz, z
    call initialize_dual(z, 4)
    call initialize_dual(fz, 4)
    z%f(0) = (0.85,0.15)
    z%f(1) = 1.0
    fz = sinh(z)
    write(*,*)"sinh()", fz%f
    fz = asinh(fz)
    write(*,*)"asinh(sinh())", fz%f
    fz = fz - z
    write(*,*)"difference", fz%f
    call check(all( abs(fz%f) < 100.0_real64 * epsilon(0.0_real64) ))
  end subroutine test_sinh_cmplx

  ! Test: cosh and acosh as inverses
  subroutine test_cosh_cmplx()
    type(dual_cmplx64) :: fz, z
    call initialize_dual(z, 4)
    call initialize_dual(fz, 4)
    z%f(0) = (0.85,0.15)
    z%f(1) = 1.0
    fz = cosh(z)
    write(*,*)"cosh()", fz%f
    fz = acosh(fz)
    write(*,*)"acosh(cosh())", fz%f
    fz = fz - z
    write(*,*)"difference", fz%f
    call check(all( abs(fz%f) < 100.0_real64 * epsilon(0.0_real64) ))
  end subroutine test_cosh_cmplx

  ! Test: tanh and atanh as inverses
  subroutine test_tanh_cmplx()
    type(dual_cmplx64) :: fz, z
    call initialize_dual(z, 4)
    call initialize_dual(fz, 4)
    z%f(0) = (0.85,0.15)
    z%f(1) = 1.0
    fz = tanh(z)
    write(*,*)"tanh()", fz%f
    fz = atanh(fz)
    write(*,*)"atanh(tanh())", fz%f
    fz = fz - z
    write(*,*)"difference", fz%f
    call check(all( abs(fz%f) < 1000.0_real64 * epsilon(0.0_real64) ))
  end subroutine test_tanh_cmplx

end module test_trig

program testtrig
  use fortuno_serial, only : execute_serial_cmd_app
  use test_trig, only : tests
  implicit none

  call execute_serial_cmd_app(tests())

end program testtrig
