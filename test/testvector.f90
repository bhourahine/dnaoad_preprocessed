!> Test app for DNAOAD to drive Fortuno unit tests
module test_vector
  use iso_fortran_env, only: real64
  use dnaoad
  use src_multivariablefunctions, only : divergence, curl, dnfvector
  use fortuno_serial, only : is_equal, test => serial_case_item, check => serial_check, test_list
  implicit none

  integer, parameter :: prec = real64

contains

  function tests()
    type(test_list) :: tests

    tests = test_list([&
      test("test_vector", test_vector_div_real)&
    ])

  end function tests

  ! Test: Divergence and curl of sourceless field
  subroutine test_vector_div_real()
    real(prec) :: q(3), div, rotor(3)
    integer :: ii

    q(:) = [1,2,3]

    div = divergence(fdivfree, q)
    write(*,*)'Divergence', div
    rotor(:) = curl(fcurlfree, q)
    write(*,*)'Curl', rotor

    call check(div < 10.0_real64 * epsilon(0.0_real64))
    call check(all(rotor < 10.0_real64 * epsilon(0.0_real64)))

  contains

    !> Example vector function with no divergence
    function fdivfree(r) result(fr)

      type(dual_real64), intent(in) :: r(:)
      type(dual_real64), allocatable :: fr(:)

      type(dual_real64) :: f1,f2,f3
      integer :: order

      order = r(1)%order()

      allocate(fr(3))
      call initialize_dual(fr, order)
      fr(:) = [r(2),r(3),r(1)]

    end function fdivfree

    !> Example vector function with no curl
    function fcurlfree(r) result(fr)

      type(dual_real64), intent(in) :: r(:)
      type(dual_real64), allocatable :: fr(:)

      type(dual_real64) :: f1,f2,f3
      integer :: order

      order = r(1)%order()

      allocate(fr(3))
      call initialize_dual(fr, order)
      fr(:) = [r(1),r(2),r(3)]

    end function fcurlfree

  end subroutine test_vector_div_real

!   ! wavenumbers (normalised such that wavenumber of external medium is 1)
!    ks=1.0_O ! external medium k
!
!    epirel = epirel / epsrel
!    muirel = muirel / musrel
!
!  scaleFactor = real(sqrt(epsrel*musrel),O)
!  ki=sqrt(epirel*muirel) Internal
!  Ci=sqrt(epirel/muirel)
!  Cs=sqrt(epsrel/musrel)


  !> Hansen vectors based on plane waves, cartesian components of wave with k-vector along
  !! (theta,phi) and polarization angle of kappa
  subroutine planewave(M, N, k, r, theta, phi, kappa)

    !> Hansen vectors
    complex(prec), intent(out) :: M(3), N(3)

    !> Dielectric constant k
    complex(prec), intent(in)  :: k

    !> Point at which to evaluate field
    real(prec), intent(in) :: r(3)

    real(prec), intent(in) :: theta, phi, kappa

    real(prec) :: M0(3), N0(3), ksvect(3)
    complex(prec) :: Mcrt(3), Ncrt(3)

    real(prec) :: st, ct, sp, cp, x(3)

    st = sin(theta)
    ct = cos(theta)
    sp = sin(phi)
    cp = cos(phi)

    x(3) = r(1) * cos(r(2))
    x(2) = r(1) * sin(r(2)) * sin(r(3))
    x(1) = r(1) * sin(r(2)) * cos(r(3))

    ! Cartesian representation
    ksvect = [st*cp, st*sp, ct]
    M0 = [ct*cp, ct*sp, -st]
    N0 = [-sp, cp, 0.0_prec]

    Mcrt = M0 * sin(k * dot_product(ksvect,x))
    Ncrt = N0 * cos(k * dot_product(ksvect,x))

  end subroutine planewave


  !> C = A X B
  pure function cross_real(A,B) result(C)

    type(dual_real64), intent(in) :: A(3)

    type(dual_real64), intent(in) :: B(3)

    type(dual_real64) :: C(3)

    C(:) = [A(2)*B(3)-A(3)*B(2), A(3)*B(1)-A(1)*B(3), A(1)*B(2)-A(2)*B(1)]

  end function cross_real


end module test_vector


program testmatrix
  use fortuno_serial, only : execute_serial_cmd_app
  use test_vector, only : tests
  implicit none

  call execute_serial_cmd_app(tests())

end program testmatrix
