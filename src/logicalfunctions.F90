!> Dual numbers for automatic differentiation up to order n
!! This can be used to compute first, second, ... nth order derivatives
!! F. Pe~nu~nuri
!! UADY, Merida Yucatan Mexico, 2024

#:include 'common.fypp'

!> Logical functions between dual numbers
module src_logicalfunctions
  use iso_fortran_env, only : ${ISO_CONSTANTS}$
  use src_types
  implicit none
  private
  public :: operator(==), operator(/=)

  !> Equality comparison
  interface operator (==)
#:for _, _, label in FLAVOURS
#:for bits in REALMODEL
    module procedure equal_dual${label}$${bits}$
#:endfor
#:endfor
  end interface operator (==)

  !> Not equality comparison
  interface operator (/=)
#:for _, _, label in FLAVOURS
#:for bits in REALMODEL
    module procedure notequal_dual${label}$${bits}$
#:endfor
#:endfor
  end interface operator (/=)

contains


#:for name, type, label in FLAVOURS
#:for bits in REALMODEL

  !> Check equality between two dual numbers of type ${label}$${bits}$
  elemental function equal_dual${label}$${bits}$(lhs, rhs) result (eqv)

    !> First dual number
    type(dual_${name}$${bits}$), intent(in) :: lhs

    !> Second dual number
    type(dual_${name}$${bits}$), intent(in) :: rhs

    !> Are these numbers equivalent
    logical :: eqv

    eqv = size(lhs%f) == size(rhs%f)
    if (eqv) eqv = all(lhs%f == rhs%f)

  end function equal_dual${label}$${bits}$

  !> Check inequality between two dual numbers of type ${label}$${bits}$
  elemental function notequal_dual${label}$${bits}$(lhs, rhs) result (neqv)

    !> First dual number
    type(dual_${name}$${bits}$), intent(in) :: lhs

    !> Second dual number
    type(dual_${name}$${bits}$), intent(in) :: rhs

    !> Are these numbers equivalent
    logical :: neqv

    neqv = size(lhs%f) /= size(rhs%f)
    if (.not. neqv) neqv = any(lhs%f /= rhs%f)

  end function notequal_dual${label}$${bits}$

#:endfor
#:endfor

end module src_logicalfunctions
