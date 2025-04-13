!> Dual numbers for automatic differentiation up to order n
!! This can be used to compute first, second, ... nth order derivatives
!! F. Pe~nu~nuri
!! UADY, Merida Yucatan Mexico, 2024

#:include 'common.fypp'

!> Data types for dual values
module src_types
  use iso_fortran_env, only : ${ISO_CONSTANTS}$
  implicit none
  private
  public :: initialize_dual

#:for name, type, _ in FLAVOURS
#:for bits in REALMODEL
  !> Dual ${type}$ variable type for ${bits}$ precision
  type, public :: dual_${name}$${bits}$
    ${type}$(kind=real${bits}$), allocatable :: f(:)
  end type dual_${name}$${bits}$
#:endfor

#:for bits in DUMMYREALMODEL
  !> Dummy ${type}$ variable type for ${bits}$ precision
  type, public :: dual_${name}$${bits}$
    ${type}$, allocatable :: f(:)
  end type dual_${name}$${bits}$
#:endfor
#:endfor

  interface initialize_dual
#:for _, _, label in FLAVOURS
#:for bits in REALMODEL
    module procedure initialize_dual_${label}$${bits}$
#:endfor
#:endfor
  end interface initialize_dual

contains

#:for name, _, label in FLAVOURS
#:for bits in REALMODEL

  !> Initialize a dual number of type ${name}$${bits}$ to be zero.
  !! Optional order of derivatives can be set (defaults to 1 otherwise).
  elemental subroutine initialize_dual_${label}$${bits}$(num, order)

    !> Dual number
    type(dual_${name}$${bits}$), intent(out) :: num

    !> Order to which the derivatives are evaluated with this value
    integer, intent(in), optional :: order

    integer :: order_

    order_ = 1
    if (present(order)) order_ = order

    allocate(num%f(0:order_)) ! notice, order+1 components
    num%f(:) = 0.0_real${bits}$

  end subroutine initialize_dual_${label}$${bits}$

#:endfor
#:endfor

end module src_types
