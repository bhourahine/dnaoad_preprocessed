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
    logical, allocatable :: status(:)
    call initialize_dual(r, 4)
    call initialize_dual(fr, 4)
    r = [0.85_real64, 1.0_real64]
    fr = sin(r)
    write(*,*)"sin(x)", fr%get_derivatives()
    fr = asin(fr)
    write(*,*)"asin(sin(x))", fr%get_derivatives()
    fr = fr - r
    write(*,*)"difference", fr%get_derivatives()
    status = abs(fr%get_derivatives()) < 100.0_real64 * epsilon(0.0_real64)
    write(*,*)status
    call check(all(status))
  end subroutine test_sin_real

  ! Test: cos and acos as inverses
  subroutine test_cos_real()
    type(dual_real64) :: fr, r
    logical, allocatable :: status(:)
    call initialize_dual(r, 4)
    call initialize_dual(fr, 4)
    r = [0.85_real64, 1.0_real64]
    fr = cos(r)
    write(*,*)"cos(x)", fr%get_derivatives()
    fr = acos(fr)
    write(*,*)"acos(cos(x))", fr%get_derivatives()
    fr = fr - r
    write(*,*)"difference", fr%get_derivatives()
    status = abs(fr%get_derivatives()) < 100.0_real64 * epsilon(0.0_real64)
    write(*,*)status
    call check(all(status))
  end subroutine test_cos_real

  ! Test: tan and atan as inverses
  subroutine test_tan_real()
    type(dual_real64) :: fr, r
    logical, allocatable :: status(:)
    call initialize_dual(r, 4)
    call initialize_dual(fr, 4)
    r = [0.85_real64, 1.0_real64]
    fr = tan(r)
    write(*,*)"tan(x)", fr%get_derivatives()
    fr = atan(fr)
    write(*,*)"atan(tan(x))", fr%get_derivatives()
    fr = fr - r
    write(*,*)"difference", fr%get_derivatives()
    status = abs(fr%get_derivatives()) < 1000.0_real64 * epsilon(0.0_real64)
    write(*,*)status
    call check(all(status))
  end subroutine test_tan_real

  ! Test: sin and asin as inverses
  subroutine test_sin_cmplx()
    type(dual_cmplx64) :: fz, z
    logical, allocatable :: status(:)
    call initialize_dual(z, 4)
    call initialize_dual(fz, 4)
    z = [(0.85_real64,0.15_real64), (1.0_real64,0.0_real64)]
    fz = sin(z)
    write(*,*)"sin(z)", fz%get_derivatives()
    fz = asin(fz)
    write(*,*)"asin(sin(z))", fz%get_derivatives()
    fz = fz - z
    write(*,*)"difference", fz%get_derivatives()
    status = abs(fz%get_derivatives()) < 1000.0_real64 * epsilon(0.0_real64)
    write(*,*)status
    call check(all(status))
  end subroutine test_sin_cmplx

  ! Test: cos and acos as inverses
  subroutine test_cos_cmplx()
    type(dual_cmplx64) :: fz, z
    logical, allocatable :: status(:)
    call initialize_dual(z, 4)
    call initialize_dual(fz, 4)
    z = [(0.85_real64,0.15_real64), (1.0_real64,0.0_real64)]
    fz = cos(z)
    write(*,*)"cos(z)", fz%get_derivatives()
    fz = acos(fz)
    write(*,*)"acos(cos(z))", fz%get_derivatives()
    fz = fz - z
    write(*,*)"difference", fz%get_derivatives()
    status = abs(fz%get_derivatives()) < 100.0_real64 * epsilon(0.0_real64)
    write(*,*)status
    call check(all(status))
  end subroutine test_cos_cmplx

  ! Test: tan and atan as inverses
  subroutine test_tan_cmplx()
    type(dual_cmplx64) :: fz, z
    logical, allocatable :: status(:)
    call initialize_dual(z, 4)
    call initialize_dual(fz, 4)
    z = [(0.85_real64,0.15_real64), (1.0_real64,0.0_real64)]
    fz = tan(z)
    write(*,*)"tan(z)", fz%get_derivatives()
    fz = atan(fz)
    write(*,*)"atan(tan(z))", fz%get_derivatives()
    fz = fz - z
    write(*,*)"difference", fz%get_derivatives()
    status = abs(fz%get_derivatives()) < 1000.0_real64 * epsilon(0.0_real64)
    write(*,*)status
    call check(all(status))
  end subroutine test_tan_cmplx

  ! Hyperbolic functions

  ! Test: sinh and asinh as inverses
  subroutine test_sinh_real()
    type(dual_real64) :: fr, r
    logical, allocatable :: status(:)
    call initialize_dual(r, 4)
    call initialize_dual(fr, 4)
    r = [0.85_real64, 1.0_real64]
    fr = sinh(r)
    write(*,*)"sinh(x)", fr%get_derivatives()
    fr = asinh(fr)
    write(*,*)"asinh(sinh(x))", fr%get_derivatives()
    fr = fr - r
    write(*,*)"difference", fr%get_derivatives()
    status = abs(fr%get_derivatives()) < 100.0_real64 * epsilon(0.0_real64)
    write(*,*)status
    call check(all(status))
  end subroutine test_sinh_real

  ! Test: cosh and acosh as inverses
  subroutine test_cosh_real()
    type(dual_real64) :: fr, r
    logical, allocatable :: status(:)
    call initialize_dual(r, 4)
    call initialize_dual(fr, 4)
    r = [0.85_real64, 1.0_real64]
    fr = cosh(r)
    write(*,*)"cosh(x)", fr%get_derivatives()
    fr = acosh(fr)
    write(*,*)"acosh(cosh(x))", fr%get_derivatives()
    fr = fr - r
    write(*,*)"difference", fr%get_derivatives()
    status = abs(fr%get_derivatives()) < 100.0_real64 * epsilon(0.0_real64)
    write(*,*)status
    call check(all(status))
  end subroutine test_cosh_real

  ! Test: tanh and atanh as inverses
  subroutine test_tanh_real()
    type(dual_real64) :: fr, r
    logical, allocatable :: status(:)
    call initialize_dual(r, 4)
    call initialize_dual(fr, 4)
    r = [0.85_real64, 1.0_real64]
    fr = tanh(r)
    write(*,*)"tanh(x)", fr%get_derivatives()
    fr = atanh(fr)
    write(*,*)"atanh(tanh(x))", fr%get_derivatives()
    fr = fr - r
    write(*,*)"difference", fr%get_derivatives()
    status = abs(fr%get_derivatives()) < 1000.0_real64 * epsilon(0.0_real64)
    write(*,*)status
    call check(all(status))
  end subroutine test_tanh_real

  ! Test: sinh and asinh as inverses
  subroutine test_sinh_cmplx()
    type(dual_cmplx64) :: fz, z
    logical, allocatable :: status(:)
    call initialize_dual(z, 4)
    call initialize_dual(fz, 4)
    z = [(0.85_real64,0.15_real64), (1.0_real64,0.0_real64)]
    fz = sinh(z)
    write(*,*)"sinh(z)", fz%get_derivatives()
    fz = asinh(fz)
    write(*,*)"asinh(sinh(z))", fz%get_derivatives()
    fz = fz - z
    write(*,*)"difference", fz%get_derivatives()
    status = abs(fz%get_derivatives()) < 100.0_real64 * epsilon(0.0_real64)
    write(*,*)status
    call check(all(status))
  end subroutine test_sinh_cmplx

  ! Test: cosh and acosh as inverses
  subroutine test_cosh_cmplx()
    type(dual_cmplx64) :: fz, z
    logical, allocatable :: status(:)
    call initialize_dual(z, 4)
    call initialize_dual(fz, 4)
    z = [(0.85_real64,0.15_real64), (1.0_real64,0.0_real64)]
    fz = cosh(z)
    write(*,*)"cosh(z)", fz%get_derivatives()
    fz = acosh(fz)
    write(*,*)"acosh(cosh(z))", fz%get_derivatives()
    fz = fz - z
    write(*,*)"difference", fz%get_derivatives()
    status = abs(fz%get_derivatives()) < 100.0_real64 * epsilon(0.0_real64)
    write(*,*)status
    call check(all(status))
  end subroutine test_cosh_cmplx

  ! Test: tanh and atanh as inverses
  subroutine test_tanh_cmplx()
    type(dual_cmplx64) :: fz, z
    logical, allocatable :: status(:)
    call initialize_dual(z, 4)
    call initialize_dual(fz, 4)
    z = [(0.85_real64,0.15_real64), (1.0_real64,0.0_real64)]
    fz = tanh(z)
    write(*,*)"tanh(z)", fz%get_derivatives()
    fz = atanh(fz)
    write(*,*)"atanh(tanh(z))", fz%get_derivatives()
    fz = fz - z
    write(*,*)"difference", fz%get_derivatives()
    status = abs(fz%get_derivatives()) < 1000.0_real64 * epsilon(0.0_real64)
    write(*,*)status
    call check(all(status))
  end subroutine test_tanh_cmplx

end module test_trig

program testtrig
  use fortuno_serial, only : execute_serial_cmd_app
  use test_trig, only : tests
  implicit none

  call execute_serial_cmd_app(tests())

end program testtrig
