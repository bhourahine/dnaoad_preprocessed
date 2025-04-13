!> Test app for DNAOAD to drive Fortuno unit tests
module test_matrix
  use iso_fortran_env, only: real64
  use dnaoad
  use fortuno_serial, only : is_equal, test => serial_case_item, check => serial_check, test_list
  implicit none

contains

  function tests()
    type(test_list) :: tests

    tests = test_list([&
      test("sum_1real", sum_real_1D),&
      test("sum_2real", sum_real_2D),&
      test("sum_3real", sum_real_3D),&
      test("product_1real", prod_real_1D)&
    ])

  end function tests

  ! Test: tests sum of matrix elements matches product
  subroutine sum_real_1D()
    integer, parameter :: n = 10
    integer, parameter :: order = 2
    type(dual_real64) :: r, rArray(n)
    call initialize_dual(r, order)
    r%f(0) = 2.0
    r%f(1) = 1.0
    rArray(:) = r
    call check(sum(rArray) == n * r)
  end subroutine sum_real_1D

  ! Test: tests sum of matrix elements matches product
  subroutine sum_real_2D()
    integer, parameter :: n = 4
    integer, parameter :: order = 2
    type(dual_real64) :: r, rArray(n,n)
    call initialize_dual(r, order)
    r%f(0) = 2.0
    r%f(1) = 1.0
    rArray(:,:) = r
    call check(sum(rArray) == n * n * r)
  end subroutine sum_real_2D

  ! Test: tests sum of matrix elements matches product
  subroutine sum_real_3D()
    integer, parameter :: n = 4
    integer, parameter :: order = 2
    type(dual_real64) :: r, rArray(n,n,n)
    call initialize_dual(r, order)
    r%f(0) = 2.0
    r%f(1) = 1.0
    rArray(:,:,:) = r
    call check(sum(rArray) == n * n * n * r)
  end subroutine sum_real_3D

  ! Test: tests product of matrix elements matches power
  subroutine prod_real_1D()
    integer, parameter :: n = 10
    integer, parameter :: order = 4
    type(dual_real64) :: r, rArray(n)
    call initialize_dual(r, order)
    r%f(0) = 2.0
    r%f(1) = 1.0
    rArray(:) = r
    call check(product(rArray) == r**n)
  end subroutine prod_real_1D

  ! Test: tests product of matrix elements matches power
  subroutine prod_real_2D()
    integer, parameter :: n = 4
    integer, parameter :: order = 4
    type(dual_real64) :: r, rArray(n,n)
    call initialize_dual(r, order)
    r%f(0) = 2.0
    r%f(1) = 1.0
    rArray(:,:) = r
    call check(product(rArray) == r**(n*n))
  end subroutine prod_real_2D

end module test_matrix


program testmatrix
  use fortuno_serial, only : execute_serial_cmd_app
  use test_matrix, only : tests
  implicit none

  call execute_serial_cmd_app(tests())

end program testmatrix
