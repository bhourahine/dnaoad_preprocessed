!> Dual numbers for automatic differentiation up to order n
!! This can be used to compute first, second, ... nth order derivatives
!! F. Pe~nu~nuri
!! UADY, Merida Yucatan Mexico, 2024

#:include 'common.fypp'

!> Helper functions for differentiation
module src_helperfunctions
  use iso_fortran_env, only : ${ISO_CONSTANTS}$
  use src_types
  implicit none
  private
  public :: f_part, BellY
#:for bits in REALMODEL
  public :: binomial${bits}$
#:endfor

  !> Extract specified order part from a dual number
  interface f_part
#:for _, _, label in FLAVOURS
#:for bits in REALMODEL
    module procedure f_part_${label}$${bits}$
#:endfor
#:endfor
  end interface f_part

  !> Partial Bell polynomial
  interface BellY
#:for _, _, label in FLAVOURS
#:for bits in REALMODEL
    module procedure BellY_${label}$${bits}$
#:endfor
#:endfor
  end interface BellY

contains

#:for bits in REALMODEL

  !> Binomial combinations for real${bits}$
  !! m! / ((m-n)! * n!)
  pure function binomial${bits}$(m, n) result(binom)

    !> Choose n from m
    integer, intent(in) :: m, n

    !> Resulting binomial
    real(real${bits}$) :: binom

    integer :: jj

    if (n == 0 .or. n == m) then
      binom = 1.0_real${bits}$
    else
      binom = 1.0_real${bits}$
      do jj = 1, n
        binom = binom * real(m - jj + 1, real${bits}$) / real(jj, real${bits}$)
      end do
    endif

  end function binomial${bits}$

#:endfor


#:for name, type, label in FLAVOURS
#:for bits in REALMODEL

  !> Extract the k-th derivative part of a dual number or an array of duals
  elemental function f_part_${label}$${bits}$(A, k) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: A

    integer, intent(in) :: k

    !> Resulting ${type}$ value
    ${type}$(real${bits}$) :: fr

    fr = A%f(k)

  end function f_part_${label}$${bits}$


  !> Partial Bell polynomial function, returning ${name}$ result with ${bits}$ bits
  pure function BellY_${label}$${bits}$(n, k, x) result(result_value)

    !> Orders of the polynomial
    integer, intent(in) :: n, k

    !> Value at which to calculate the polynomial
    ${type}$(real${bits}$), dimension(:), intent(in) :: x

    !> Resulting polynomial value
    ${type}$(real${bits}$) :: result_value

    ${type}$(real${bits}$) :: dp(n+1,k+1), sum_val
    ${type}$(real${bits}$), allocatable :: newx(:)
    integer :: nn, kk, ii, LX, nx

    ! Special cases
    if (n == 0 .and. k == 0) then
      result_value = 1.0_real${bits}$
      return
    elseif (size(x) == 0) then
      result_value = 0.0_real${bits}$
      return
    end if

    dp(:, :) = 0.0_real${bits}$
    dp(1, 1) = 1.0_real${bits}$

    ! Main loop to compute BellY[n, k, x]
    do nn = 1, n
      do kk = 1, k
        LX = size(x)
        nx = max(nn - kk + 1, LX)
        if (nx > 0) then
          allocate(newx(nx))
          newx = 0.0_real${bits}$
          newx(1:LX) = x
          sum_val = 0.0_real${bits}$
          do ii = 0, nn - kk
            sum_val = sum_val + binomial${bits}$(nn - 1, ii) * newx(ii + 1) * dp(nn - ii, kk)
          end do
          dp(nn + 1, kk + 1) = sum_val
          deallocate(newx)
        end if
      end do
    end do

    result_value = dp(n + 1, k + 1)

  end function BellY_${label}$${bits}$

#:endfor
#:endfor

end module src_helperfunctions
