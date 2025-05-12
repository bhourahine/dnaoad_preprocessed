!> Dual numbers for automatic differentiation up to order n
!! This can be used to compute first, second, ... nth order derivatives
!! F. Pe~nu~nuri
!! UADY, Merida Yucatan Mexico, 2024

#:include 'common.fypp'

!> Data types for dual values
module src_multivariablefunctions
  use iso_fortran_env, only : ${ISO_CONSTANTS}$
  use src_types
  implicit none
  private
  public :: dnfscalar, dnfvector, gradient, hessian, jacobian

  !> Abstract interfaces for functions passed as arguments to other
  !! functions.

#:for name, type, label in FLAVOURS
#:for bits in REALMODEL
  !> Interface for a scalar dual function f: D^m --> D (similar to f: R^m --> R)
  !! Interface to generically define application of function to a ${type}$ array, returning a dual
  !! variable result of type ${type}$ and ${bits}$ precision.
  abstract interface
    function fsdual_${label}$${bits}$(xd) result(frsd)
      use iso_fortran_env, only: real${bits}$
      import :: dual_${name}$${bits}$
      !> Array of dual values to apply the function to
      type(dual_${name}$${bits}$), intent(in) :: xd(:)
      !> Resulting dual variable
      type(dual_${name}$${bits}$) :: frsd
    end function fsdual_${label}$${bits}$
  end interface

  !> Interface for a vector valued dual function f: D^m --> D^n
  !! Interface to generically define application of function to a ${type}$ array, returning a vector
  !! of dual variables of type ${type}$ and ${bits}$ precision.
  interface
    function fvecdual_${label}$${bits}$(xd) result(frd)
      use iso_fortran_env, only: real${bits}$
      import :: dual_${name}$${bits}$
      !> Array of dual values to apply the function to
      type(dual_${name}$${bits}$), intent(in) :: xd(:)
      !> Array valued results of the function
      type(dual_${name}$${bits}$), allocatable :: frd(:)
    end function fvecdual_${label}$${bits}$
  end interface

#:endfor
#:endfor

  interface dnfscalar
#:for name, type, label in FLAVOURS
#:for bits in REALMODEL
    module procedure dnfscalar_${label}$${bits}$
#:endfor
#:endfor
  end interface dnfscalar

    interface dnfvector
#:for name, type, label in FLAVOURS
#:for bits in REALMODEL
    module procedure dnfvector_${label}$${bits}$
#:endfor
#:endfor
  end interface dnfvector

  interface d2fdudv
#:for _, _, label in FLAVOURS
#:for bits in REALMODEL
    module procedure d2fdudv_${label}$${bits}$
#:endfor
#:endfor
  end interface d2fdudv

  interface gradient
#:for name, type, label in FLAVOURS
#:for bits in REALMODEL
    module procedure gradient_${label}$${bits}$
#:endfor
#:endfor
  end interface gradient

  interface hessian
#:for name, type, label in FLAVOURS
#:for bits in REALMODEL
    module procedure hessian_${label}$${bits}$
#:endfor
#:endfor
  end interface hessian

    interface jacobian
#:for name, type, label in FLAVOURS
#:for bits in REALMODEL
    module procedure jacobian_${label}$${bits}$
#:endfor
#:endfor
  end interface jacobian

contains

#:for name, type, label in FLAVOURS
#:for bits in REALMODEL

  !> Directional derivative to order n of a scalar function along vector v, evaluated at point q
  function dnfscalar_${label}$${bits}$(fsd, v, q, n) result(fr)

    !> Function to differentiate
    procedure(fsdual_${label}$${bits}$) :: fsd

    !> Direction vector
    ${type}$(real${bits}$), intent(in) :: v(:)

    !> Point at which to evaluate
    ${type}$(real${bits}$), intent(in) :: q(:)

    !> Order of derivative
    integer, intent(in) :: n

    !> Resulting scalar directional derivative
    type(dual_${name}$${bits}$) :: fr

    type(dual_${name}$${bits}$) :: eps

    ${type}$(real${bits}$), parameter :: one = 1.0_real${bits}$

    call initialize_dual(eps, n)
    call eps%set_derivative(1, one)

    fr = fsd(q + eps * v)

  end function dnfscalar_${label}$${bits}$

!  !> Mixed derivative at values q (unchecked)
!  function dnfmixed_${label}$${bits}$(fsd, q) result(fr)
!
!    !> Function to differentiate
!    procedure(fsdual_${label}$${bits}$) :: fsd
!
!    !> Points at which to evaluate derivatives
!    ${type}$(real${bits}$), intent(in) :: q(:)
!
!    !> Resulting derivative
!    type(dual_${name}$${bits}$) :: fr
!
!    type(dual_${name}$${bits}$) :: eps(size(q))
!
!    call initialize_dual(eps, n)
!    eps(:)%f(1) = 1
!
!    fr = fsd(q + sum(eps)) -
!
!  end function dnfmixed_${label}$${bits}$

  !> Directional derivative, up to order, of a n-dimensional vector-valued function along vector v,
  !! evaluated at point q
  function dnfvector_${label}$${bits}$(fvecd, v, q, n, order) result(fr)

    !> Function to differentiate
    procedure(fvecdual_${label}$${bits}$) :: fvecd

    !> Direction vector
    ${type}$(real${bits}$), intent(in) :: v(:)

    !> Point at which to evaluate
    ${type}$(real${bits}$), intent(in) :: q(:)

    !> Size of the vector returned by fvecd
    integer, intent(in) :: n

    !> Order of derivative
    integer, intent(in) :: order

    !> Resulting scalar directional derivative
    type(dual_${name}$${bits}$), allocatable :: fr(:)

    type(dual_${name}$${bits}$) :: eps

    ${type}$(real${bits}$), parameter :: one = 1.0_real${bits}$

    allocate(fr(n))
    call initialize_dual(fr, order)
    call initialize_dual(eps, order)

    call eps%set_derivative(1, one)
    fr(:) = fvecd(q + eps * v)

  end function dnfvector_${label}$${bits}$

  !> Mixed directional 2nd derivative of a scalar valued function along vectors u and v, evaluated
  !! at point q
  function d2fdudv_${label}$${bits}$(fsd, u, v, q) result(fr)

    !> Function to differentiate
    procedure(fsdual_${label}$${bits}$) :: fsd

    !> First directional vector
    ${type}$(real${bits}$), intent(in) :: u(:)

    !> Second directional vector
    ${type}$(real${bits}$), intent(in) :: v(:)

    !> Point at which to evaluate function derivatives
    ${type}$(real${bits}$), intent(in) :: q(:)

    !> Resulting scalar directional derivative
    ${type}$(real${bits}$) :: fr

    type(dual_${name}$${bits}$) :: ftmp

    ftmp = dnfscalar(fsd, u + v, q, 2) - dnfscalar(fsd, u, q, 2) - dnfscalar(fsd, v, q, 2)
    fr = 0.5_real${bits}$ * f_part(ftmp, 2)

  end function d2fdudv_${label}$${bits}$

  !> Evaluate gradients (i.e. 1st derivative) of scalar function fsd at point q
  function gradient_${label}$${bits}$(fsd, q) result(fr)

    !> Scalar valued function to evaluate the gradient
    procedure(fsdual_${label}$${bits}$) :: fsd

    !> Point at which to evaluate the derivative
    ${type}$(real${bits}$), intent(in) :: q(:)

    ${type}$(real${bits}$) :: fr(size(q)), ei(size(q))
    integer :: ii

    do ii = 1, size(q)
      ei = 0
      ei(ii) = 1
      fr(ii) = f_part(dnfscalar(fsd, ei, q, 1), 1)
    end do

  end function gradient_${label}$${bits}$

  !> Evaluate the hessian matrix of a scalar function fsd at point q
  function hessian_${label}$${bits}$(fsd, q) result(Hmat)

    !> Scalar valued function to evaluate the gradient
    procedure(fsdual_${label}$${bits}$) :: fsd

    !> Point at which to evaluate the derivative
    ${type}$(real${bits}$), intent(in) :: q(:)

    ${type}$(real${bits}$) :: Hmat(size(q), size(q))

    ${type}$(real${bits}$) :: ei(size(q)), ej(size(q))
    integer :: ii, jj

    Hmat(:,:) = 0
    do ii = 1, size(q)
      ei(:) = 0
      ei(ii) = 1
      Hmat(ii, ii) = f_part(dnfscalar(fsd, ei, q, 2), 2)
      do jj = 1, size(q)
        if (ii == jj) cycle
        ej(:) = 0
        ej(jj) = 1
        Hmat(jj, ii) = d2fdudv(fsd, ei, ej, q)
      end do
    end do

  end function hessian_${label}$${bits}$

  !> Evaluate the jacobian matrix of a vector function fvecd at point q
  function jacobian_${label}$${bits}$(fvecd, q, n) result(Jmat)

    !> Scalar valued function to evaluate the gradient
    procedure(fvecdual_${label}$${bits}$) :: fvecd

    !> Point at which to evaluate the derivative
    ${type}$(real${bits}$), intent(in) :: q(:)

    !> Number of dependent variables (i.e. return size of fvecd array)
    integer, intent(in) :: n

    ${type}$(real${bits}$) :: Jmat(n, size(q)), ei(size(q))
    integer :: ii

    Jmat(:,:) = 0
    do ii = 1, size(q)
      ei(:) = 0
      ei(ii) = 1
      Jmat(:, ii) = f_part(dnfvector(fvecd, ei, q, n, 1), 1)
    end do

  end function jacobian_${label}$${bits}$

#:endfor
#:endfor

end module src_multivariablefunctions
