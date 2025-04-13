!> Dual numbers for automatic differentiation up to order n
!! This can be used to compute first, second, ... nth order derivatives
!! F. Pe~nu~nuri
!! UADY, Merida Yucatan Mexico, 2024

#:include 'common.fypp'

module src_constants
  use iso_fortran_env, only : ${ISO_CONSTANTS}$
  implicit none

  public

#:for bits in REALMODEL
  !> Value of pi for real${bits}$ type
  real(real${bits}$), public, parameter :: Pi_${bits}$ = 4.0_real${bits}$ * atan(1.0_real${bits}$)

#:endfor
#:for bits in DUMMYREALMODEL
  !> Dummy value of pi for an unsupported real${bits}$ type
  real, public, parameter :: Pi_${bits}$ = -1.0

#:endfor

end module src_constants
