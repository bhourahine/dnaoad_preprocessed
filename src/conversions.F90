!> Dual numbers for automatic differentiation up to order n
!! This can be used to compute first, second, ... nth order derivatives
!! F. Pe~nu~nuri
!! UADY, Merida Yucatan Mexico, 2024

#:include 'common.fypp'

!> Conversions between types of dual numbers
module src_conversions
  use iso_fortran_env, only : ${ISO_CONSTANTS}$
  use src_types

  implicit none
  private
  public :: aimag
#:for outbits in REALMODEL
  public :: real_${outbits}$, cmplx_${outbits}$
#:endfor
#:for outbits in REALMODEL

  !> Interface for conversion to a dual number of type real_${outbits}$
  interface real_${outbits}$
#:for _, _, label in FLAVOURS
#:for inbits in REALMODEL
    module procedure real${outbits}$_${label}$${inbits}$
#:endfor
#:endfor
  end interface real_${outbits}$

#:endfor
#:for outbits in REALMODEL

  !> Interface for conversion to a dual number of type cmplx_${outbits}$
  interface cmplx_${outbits}$
#:for inbits in REALMODEL
    module procedure cmplx${outbits}$_c${inbits}$
    module procedure cmplx${outbits}$_r${inbits}$
#:for bitsopt in REALMODEL
    module procedure cmplx${outbits}$_r${inbits}$_opt${bitsopt}$
#:endfor
#:endfor
  end interface cmplx_${outbits}$

#:endfor

  !> Interface for extracting imaginary part of a complex dual number as a real dual
  interface aimag
#:for bits in REALMODEL
    module procedure aimag_c${bits}$
#:endfor
  end interface aimag

contains

#:for outbits in REALMODEL
#:for name, type, label in FLAVOURS
#:for inbits in REALMODEL

  !> Converts a ${type}$ ${inbits}$-bit dual value to a real ${inbits}$-bit dual
  elemental function real${outbits}$_${label}$${inbits}$(A) result(B)

    !> Dual number
    type(dual_${name}$${inbits}$), intent(in) :: A

    !> Dual number
    type(dual_real${outbits}$) :: B

    integer :: order

    order = size(A%f) - 1
    call initialize_dual(B, order)

    B%f(:) = real(A%f, kind=real${outbits}$)

  end function real${outbits}$_${label}$${inbits}$

#:endfor
#:endfor
#:endfor

#:for outbits in REALMODEL

  !> Converts a complecx ${outbits}$-bit dual value to a real ${outbits}$-bit dual
  elemental function aimag_c${outbits}$(A) result(B)

    !> Dual number
    type(dual_cmplx${outbits}$), intent(in) :: A

    !> Dual number
    type(dual_real${outbits}$) :: B

    integer :: order

    order = size(A%f) - 1
    call initialize_dual(B, order)

    B%f(:) = aimag(A%f)

  end function aimag_c${outbits}$

#:for inbits in REALMODEL

  !> Converts a complex ${inbits}$-bit dual value to a complex ${outbits}$-bit dual
  elemental function cmplx${outbits}$_c${inbits}$(A) result(B)

    !> Dual number
    type(dual_cmplx${inbits}$), intent(in) :: A

    !> Dual number
    type(dual_cmplx${outbits}$) :: B

    integer :: order

    order = size(A%f) - 1
    call initialize_dual(B, order)

    B%f(:) = cmplx(A%f, kind=real${outbits}$)

  end function cmplx${outbits}$_c${inbits}$

  !> Converts a real ${inbits}$-bit dual value to a complex ${outbits}$-bit dual
  elemental function cmplx${outbits}$_r${inbits}$(A) result(C)

    !> Dual number for real part
    type(dual_real${inbits}$), intent(in) :: A

    !> Dual number
    type(dual_cmplx${outbits}$) :: C

    integer :: order

    order = size(A%f) - 1
    call initialize_dual(C, order)

    C%f(:) = cmplx(A%f, kind=real${outbits}$)

  end function cmplx${outbits}$_r${inbits}$

#:for bitsopt in REALMODEL

  !> Converts a real ${inbits}$-bit dual value to a complex ${outbits}$-bit dual
  elemental function cmplx${outbits}$_r${inbits}$_opt${bitsopt}$(A, B) result(C)

    !> Dual number for real part
    type(dual_real${inbits}$), intent(in) :: A

    !> Dual number for imaginary part
    type(dual_real${bitsopt}$), intent(in) :: B

    !> Dual number
    type(dual_cmplx${outbits}$) :: C

    integer :: order

    order = min(size(A%f), size(B%f)) - 1
    call initialize_dual(C, order)

    C%f(:) = cmplx(A%f(0:order), B%f(0:order), kind=real${outbits}$)

  end function cmplx${outbits}$_r${inbits}$_opt${bitsopt}$

#:endfor

#:endfor

#:endfor

end module src_conversions
