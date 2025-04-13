!> Test app for DNAOAD to drive Fortuno unit tests
module test_typednaoad
  use dnaoad
  use fortuno_serial, only : is_equal, test => serial_case_item, check => serial_check, test_list
  implicit none

contains

  function tests()
    type(test_list) :: tests

    tests = test_list([&
      test("initialize_0", test_initialize_0),&
      test("initialize_1", test_initialize_1),&
      test("eqv", test_eqv),&
      test("neqv", test_neqv)&
    ])

  end function tests

  ! Test: initializing scalar value
  subroutine test_initialize_0()
    type(dual_real64) :: r
    call initialize_dual(r, 4)
    call check(size(r%f) == 5)
  end subroutine test_initialize_0

  ! Test: initializing vector of values
  subroutine test_initialize_1()
    type(dual_real64) :: r(4)
    call initialize_dual(r, 4)
    call check(size(r) == 4)
    call check(size(r(1)%f) == 5)
  end subroutine test_initialize_1

  ! Test: equivalence between dual numbers
  subroutine test_eqv()
    type(dual_real64) :: r, s, t
    call initialize_dual(r, 4)
    call initialize_dual(s, 4)
    call initialize_dual(t, 3)
    r%f(:) = [0,1,2,3,4]
    s%f(:) = [0,1,2,3,4]
    call check(r == s)
    call check(.not.(r == t))
  end subroutine test_eqv

  ! Test: inequivalence between dual numbers
  subroutine test_neqv()
    type(dual_real64) :: r, s, t
    call initialize_dual(r, 4)
    call initialize_dual(s, 4)
    call initialize_dual(t, 3)
    r%f(:) = [0,1,2,3,4]
    s%f(:) = [1,2,3,4,5]
    call check(r /= s)
    call check(r /= t)
  end subroutine test_neqv

end module test_typednaoad


program testtype
  use fortuno_serial, only : execute_serial_cmd_app
  use test_typednaoad, only : tests
  implicit none

  call execute_serial_cmd_app(tests())

end program testtype
