!> Dual numbers for automatic differentiation up to order n
!! This can be used to compute first, second, ... nth order derivatives
!! F. Pe~nu~nuri
!! UADY, Merida Yucatan Mexico, 2024

#:include 'common.fypp'

!> Chain rule differentiation
module src_chainrule
  use iso_fortran_env, only : ${ISO_CONSTANTS}$
  use src_helperfunctions, only : BellY
  use src_types
  implicit none
  private
  public :: Dnd

  !> Driver for differentiation of a function
  interface Dnd
#:for _, _, label in FLAVOURS
#:for bits in REALMODEL
    module procedure Dnd_${label}$${bits}$
#:endfor
#:endfor
  end interface Dnd

#:for name, type, label in FLAVOURS
#:for bits in REALMODEL
  !> Interface to generically define application of function to a ${type}$ variable, returning a
  !! dual variable result of type ${type}$ and ${bits}$ precision.
  abstract interface
    pure function fun_dual_${label}$${bits}$(val, order) result(f_result)
      use iso_fortran_env, only: real${bits}$
      import :: dual_${name}$${bits}$
      !> Value to apply the function to
      ${type}$(real${bits}$), intent(in) :: val
      !> Order of dual variable being processed
      integer, intent(in) :: order
      !> Resulting dual variable
      type(dual_${name}$${bits}$) :: f_result
    end function fun_dual_${label}$${bits}$
  end interface
#:endfor
#:endfor

contains


#:for name, type, label in FLAVOURS
#:for bits in REALMODEL

  !> Chain rule, D^n (f(g)), to differentiate functions of dual numbers
  pure function Dnd_${label}$${bits}$(fc, gdual, n) result(dnfc)

    !> Procedure being applied to the number
    procedure(fun_dual_${label}$${bits}$) :: fc

    !> Dual value
    type(dual_${name}$${bits}$), intent(in) :: gdual

    !> Order of derivative
    integer, intent(in) :: n

    !> Resulting term for order n derivative
    ${type}$(real${bits}$) :: dnfc

    ${type}$(real${bits}$) :: g0, suma
    type(dual_${name}$${bits}$) :: fvd
    ${type}$(real${bits}$), allocatable :: xvg(:)
    integer :: kk, jj, order

    order = size(gdual%f) - 1
    g0 = gdual%f(0)
    fvd = fc(g0, order)
    if(n == 0) then
      dnfc = fvd%f(0)
    else
      suma = 0
      do kk = 1, n
        allocate(xvg(1:n-kk+1))
        do jj=1,n-kk+1
          xvg(jj) = gdual%f(jj)
        end do
        suma = suma + fvd%f(kk) * BellY(n,kk,xvg)
        deallocate(xvg)
      end do
      dnfc = suma
    end if

  end function Dnd_${label}$${bits}$

#:endfor
#:endfor


end module src_chainrule
