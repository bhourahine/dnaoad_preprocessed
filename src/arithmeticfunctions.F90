!> Dual numbers for automatic differentiation up to order n
!! This can be used to compute first, second, ... nth order derivatives
!! F. Pe~nu~nuri
!! UADY, Merida Yucatan Mexico, 2024

#:include 'common.fypp'

!> Arithmetic functions between numbers and dual numbers
module src_arithmeticfunctions
  use iso_fortran_env, only : ${ISO_CONSTANTS}$
  use src_chainrule, only : Dnd
  use src_constants
#:for bits in REALMODEL
  use src_conversions, only : cmplx_${bits}$
#:endfor
  use src_helperfunctions
  use src_logicalfunctions
  use src_types

  implicit none
  private
  public :: assignment (=)
  public :: operator(+), operator(-), operator(*), operator(/), operator(**)
  public :: inv, sin, cos, tan, exp, log, sqrt, asin, acos, atan, asinh, acosh
  public :: atanh, sinh, cosh, tanh, absx, atan2
  public :: sind, cosd, tand, asind, acosd, atand, log10
  public :: conjg
  public :: isInteger

  ! To do: mixed precision for operators, bessel_*, erf, erfc, erfc_scaled, gamma, log_gamma

  !> Assignment operator
  interface assignment (=)
#:for _, _, label in FLAVOURS
#:for bits in REALMODEL
#:for intbits in INTMODEL
    module procedure assign_int${intbits}$_to_dual${label}$${bits}$
#:endfor
    module procedure assign_${label}$${bits}$_to_dual${label}$${bits}$
    module procedure assign_dual${label}$${bits}$_to_dual${label}$${bits}$
#:endfor
#:endfor
  end interface assignment (=)

#! Loop over add, minus, times and div interfaces
#:for op, name, Name in [("+","add", "Plus"), ("-","minus", "Minus"),&
  & ("*","times", "Times"), ("/","div", "Divide"), ("**", "power", "Exponentiate")]

  !> ${Name}$ operation involving dual numbers and integers, floating point or other dual numbers
  interface operator (${op}$)
#:for _, _, label in FLAVOURS
#:for bits in REALMODEL
#:if op == "+" or op == "-"
    module procedure unary${name}$_${label}$${bits}$
#:endif
#:for intbits in INTMODEL
    module procedure ${name}$_int${intbits}$_to_dual${label}$${bits}$
    module procedure ${name}$_dual${label}$${bits}$_to_int${intbits}$
#:endfor
    module procedure ${name}$_${label}$${bits}$_to_dual${label}$${bits}$
    module procedure ${name}$_dual${label}$${bits}$_to_${label}$${bits}$
#:for _, _, label2 in FLAVOURS
    module procedure ${name}$Duals_${label}$${bits}$_${label2}$${bits}$
#:endfor
#:endfor
#:endfor
  end interface operator (${op}$)

#:endfor


#:for op in [("exp"), ("log"), ("sin"), ("cos"), ("tan"), ("sqrt"), ("absx"), ("sinh"), ("cosh"),&
  & ("tanh"), ("asin"), ("acos"), ("asinh"), ("acosh"), ("atanh"), ("atan2"), ("sind"), ("cosd"),&
      & ("tand"), ("asind"), ("acosd"), ("atand"), ("log10")]
  !> Interface for ${op}$
  interface ${op}$
#:for _, _, label in FLAVOURS
#:for bits in REALMODEL
    module procedure ${op}$_${label}$${bits}$
#:endfor
#:endfor
  end interface ${op}$
#:endfor

  !> Interface for atan
  interface atan
#:for _, _, label in FLAVOURS
#:for bits in REALMODEL
    module procedure atan_${label}$${bits}$
    ! Fortran 2008 standard for 2 argument atan:
    module procedure atan2_${label}$${bits}$
#:endfor
#:endfor
  end interface atan


  !> Conjugation operation.
  !! Notice that the conjugation operation is not forward differentiable. We define this operation
  !! as (df)* not d(f*).
  interface conjg
#:for bits in REALMODEL
    module procedure conjg_${bits}$
#:endfor
  end interface conjg


#:for describe, op in [("Sum", "add"), ("Difference", "minus"), ("Product", "times"),&
  & ("Ratio", "div"), ("Exponentiation", "power")]
  !> ${describe}$ between two dual numbers
  interface ${op}$Duals
#:for _, _, label in FLAVOURS
#:for _, _, label2 in FLAVOURS
#:for bits in REALMODEL
    module procedure ${op}$Duals_${label}$${bits}$_${label2}$${bits}$
#:endfor
#:endfor
#:endfor

  end interface ${op}$Duals

#:endfor


  !> Sub-term helper for multiplying two numbers together
  interface times_subterm
#:for _, _, label in FLAVOURS
#:for _, _, label2 in FLAVOURS
#:for bits in REALMODEL
    module procedure times_subterm_${label}$${bits}$_${label2}$${bits}$
#:endfor
#:endfor
#:endfor
  end interface times_subterm

  !> Inverse of a dual value (i.e., 1/value)
  interface inv
#:for _, _, label in FLAVOURS
#:for bits in REALMODEL
    module procedure inv_${label}$${bits}$
#:endfor
#:endfor
  end interface inv

  !> Inverse evaluation sub-term helper
  interface inv_subterm
#:for _, _, label in FLAVOURS
#:for bits in REALMODEL
    module procedure inv_subterm_${label}$${bits}$
#:endfor
#:endfor
  end interface inv_subterm

  !> Tests if a dual number is equivalent to an integer
  interface isInteger
#:for _, _, label in FLAVOURS
#:for bits in REALMODEL
    module procedure isInteger_dual${label}$${bits}$
    module procedure isInteger_${label}$${bits}$
#:endfor
#:endfor
  end interface isInteger

  !> Internal routine for complex atan2
  interface atan2_z
#:for bits in REALMODEL
    module procedure atan2_z_${bits}$
#:endfor
  end interface atan2_z

contains

#:for name, type, label in FLAVOURS
#:for bits in REALMODEL

#! Integers with dual numbers
#:for intbits in INTMODEL
  !> Assign an int${intbits}$ value to an already initialized dual number of type ${name}$${bits}$
  !! Note, the order of the dual number is retained in the operation
  elemental subroutine assign_int${intbits}$_to_dual${label}$${bits}$(A, ii)

    !> Resulting dual number
    type(dual_${name}$${bits}$), intent(inout) :: A

    !> Integer value to copy into the dual number
    integer(int${intbits}$), intent(in) :: ii

    A%f(0) = ${name}$(ii)
    A%f(1:) = ${name}$(0)

  end subroutine assign_int${intbits}$_to_dual${label}$${bits}$
#:endfor

  !> Assign a ${name}$${bits}$ value to an already initialized dual number of type ${name}$${bits}$
  !! Note, the order of the output dual number is retained in the operation
  elemental subroutine assign_${label}$${bits}$_to_dual${label}$${bits}$(A, x)

    !> Resulting dual number
    type(dual_${name}$${bits}$), intent(inout) :: A

    !> Value to copy into the dual number
    ${type}$(real${bits}$), intent(in) :: x

    A%f(0) = x
    A%f(1:) = 0

  end subroutine assign_${label}$${bits}$_to_dual${label}$${bits}$

  !> Assign a dual number of type ${name}$${bits}$ value to a dual number of type
  !! ${name}$${bits}$
  !! Note, the order of the input dual number is retained in the output value
  elemental subroutine assign_dual${label}$${bits}$_to_dual${label}$${bits}$(A, B)

    !> Resulting dual number
    type(dual_${name}$${bits}$), intent(out) :: A

    !> Dual number to copy
    type(dual_${name}$${bits}$), intent(in) :: B

    A%f = B%f

  end subroutine assign_dual${label}$${bits}$_to_dual${label}$${bits}$

#:endfor
#:endfor


#! Loop over add, minus, times and divide operations between dual number and another variable
#:for comment, op in [("+","add"), ("-","minus"), ("*","times"), ("/","div")]

#:for name, type, label in FLAVOURS
#:for bits in REALMODEL

#:if comment == "+" or comment == "-"

  !> Unary ${op}$ operation on a number
  elemental function unary${op}$_${label}$${bits}$(A) result (fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: A

    !> Resulting dual number
    type(dual_${name}$${bits}$) :: fr

    integer :: order

    order = size(A%f) - 1
    call initialize_dual(fr, order)

    fr%f(:) = ${comment}$ A%f

  end function unary${op}$_${label}$${bits}$

#:endif

#! Integers with dual numbers
#:for intbits in INTMODEL

  !> Integer ${comment}$ dual value
  elemental function ${op}$_int${intbits}$_to_dual${label}$${bits}$(ii, A) result(fr)

    !> Integer value
    integer(int${intbits}$), intent(in) :: ii

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: A

    !> Resulting dual number
    type(dual_${name}$${bits}$) :: fr

    type(dual_${name}$${bits}$) :: tmpDual
    integer :: order

    order = size(A%f) - 1
    call initialize_dual(tmpDual, order)
    call initialize_dual(fr, order)

    tmpDual%f(0) = ${name}$(ii)
    fr = ${op}$Duals(tmpDual, A)

  end function ${op}$_int${intbits}$_to_dual${label}$${bits}$


  !> Dual value ${comment}$ Integer
  elemental function ${op}$_dual${label}$${bits}$_to_int${intbits}$(A, ii) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: A

    !> Integer value
    integer(int${intbits}$), intent(in) :: ii

    !> Resulting dual number
    type(dual_${name}$${bits}$) :: fr

    type(dual_${name}$${bits}$) :: tmpDual
    integer :: order

    order = size(A%f) - 1
    call initialize_dual(tmpDual, order)
    call initialize_dual(fr, order)

    tmpDual%f(0) = ${name}$(ii)
    fr = ${op}$Duals(A, tmpDual)

  end function ${op}$_dual${label}$${bits}$_to_int${intbits}$

#:endfor

  !> ${op}$ ${label}$${bits}$ and ${comment}$ dual value
  elemental function ${op}$_${label}$${bits}$_to_dual${label}$${bits}$(x, A) result(fr)

    !> ${type}$ value
    ${type}$(real${bits}$), intent(in) :: x

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: A

    !> Resulting dual number
    type(dual_${name}$${bits}$) :: fr

    type(dual_${name}$${bits}$) :: tmpDual
    integer :: order

    order = size(A%f) - 1
    call initialize_dual(tmpDual, order)
    call initialize_dual(fr, order)

    tmpDual%f(0) = x
    fr = ${op}$Duals(tmpDual, A)

  end function ${op}$_${label}$${bits}$_to_dual${label}$${bits}$


  !> ${op}$ dual value ${comment}$ and ${label}$${bits}$
  elemental function ${op}$_dual${label}$${bits}$_to_${label}$${bits}$(A, x) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: A

    !> ${type}$ value
    ${type}$(real${bits}$), intent(in) :: x

    !> Resulting dual number
    type(dual_${name}$${bits}$) :: fr

    type(dual_${name}$${bits}$) :: tmpDual
    integer :: order

    order = size(A%f) - 1
    call initialize_dual(tmpDual, order)
    call initialize_dual(fr, order)

    tmpDual%f(0) = x
    fr = ${op}$Duals(A, tmpDual)

  end function ${op}$_dual${label}$${bits}$_to_${label}$${bits}$

#:endfor
#:endfor

#:endfor


#:for name, type, label in FLAVOURS
#:for bits in REALMODEL

#! Integers with dual numbers
#:for intbits in INTMODEL

  !> Raise integer ** dual value
  elemental function power_int${intbits}$_to_dual${label}$${bits}$(ii, A) result(fr)

    !> Integer value
    integer(int${intbits}$), intent(in) :: ii

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: A

    !> Resulting dual number
    type(dual_${name}$${bits}$) :: fr

    type(dual_${name}$${bits}$) :: tmpDual
    integer :: order

    order = size(A%f) - 1
    call initialize_dual(tmpDual, order)
    call initialize_dual(fr, order)

    tmpDual%f(0) = ${name}$(ii)
    fr = powerDuals(tmpDual, A)

  end function power_int${intbits}$_to_dual${label}$${bits}$


  !> Raise dual value ** Integer
  elemental function power_dual${label}$${bits}$_to_int${intbits}$(A, n) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: A

    !> Integer value
    integer(int${intbits}$), intent(in) :: n

    !> Resulting dual number
    type(dual_${name}$${bits}$) :: fr

    type(dual_${name}$${bits}$) :: tmpDual
    integer(int${intbits}$) :: kk
    integer :: order

    order = size(A%f) - 1
    call initialize_dual(tmpDual, order)
    call initialize_dual(fr, order)

    tmpDual = 1
    if (A == tmpDual) then
      fr = A
      return
    elseif (n > 0) then
      tmpDual = 0
      if (A == tmpDual) then
        fr = A
        return
      end if
    end if

    ! 0**0 ---> 1
    fr = 1

    if (n /= 0) then
      do kk = 1, abs(n)
        fr = fr * A
      end do
    end if
    if (n < 0) then
      fr = inv(fr)
    end if

  end function power_dual${label}$${bits}$_to_int${intbits}$

#:endfor

  !> Raise ${label}$${bits}$ ** dual value
  elemental function power_${label}$${bits}$_to_dual${label}$${bits}$(x, B) result(fr)

    !> ${type}$ value
    ${type}$(real${bits}$), intent(in) :: x

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: B

    !> Resulting dual number
    type(dual_${name}$${bits}$) :: fr

    type(dual_${name}$${bits}$) :: tmpDual
    integer :: order

    order = size(B%f) - 1
    call initialize_dual(tmpDual, order)
    call initialize_dual(fr, order)

    tmpDual%f(0) = x
    fr = powerDuals(tmpDual, B)

  end function power_${label}$${bits}$_to_dual${label}$${bits}$


  !> Raise dual value ** ${label}$${bits}$
  elemental function power_dual${label}$${bits}$_to_${label}$${bits}$(A, x) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: A

    !> ${type}$ value
    ${type}$(real${bits}$), intent(in) :: x

    !> Resulting dual number
    type(dual_${name}$${bits}$) :: fr

    type(dual_${name}$${bits}$) :: tmpDual
    integer :: order

    order = size(A%f) - 1
    call initialize_dual(fr, order)

    if (isInteger(x)) then

      fr = A ** int(x)

    else

      call initialize_dual(tmpDual, order)
      tmpDual%f(0) = x
      fr = powerDuals(A, tmpDual)

    end if

  end function power_dual${label}$${bits}$_to_${label}$${bits}$

#:endfor
#:endfor



#:for name, type, label in FLAVOURS
#:for bits in REALMODEL

#:for name2, type2, label2 in FLAVOURS

#:for op, txt, routine in [("+", "Add together", "add"), ("-", "Difference between", "minus")]

  !> ${txt}$ two dual variables
  elemental function ${routine}$Duals_${label}$${bits}$_${label2}$${bits}$(A, B) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: A

    !> Dual number
    type(dual_${name2}$${bits}$), intent(in) :: B

    !> Dual number (= A ${op}$ B)
#:if (label == 'c' or label2 == 'c')
    type(dual_cmplx${bits}$) :: fr
#:else
    type(dual_real${bits}$) :: fr
#:endif

    integer :: ii, order

    order = min(size(A%f), size(B%f)) - 1
    call initialize_dual(fr, order)
    do ii = 0, order
      fr%f(ii) = A%f(ii) ${op}$ B%f(ii)
    end do

  end function ${routine}$Duals_${label}$${bits}$_${label2}$${bits}$

#:endfor

  !> Product between two dual numbers
  elemental function timesDuals_${label}$${bits}$_${label2}$${bits}$(A, B) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: A

    !> Dual number
    type(dual_${name2}$${bits}$), intent(in) :: B

    !> Dual number (= A * B)
#:if (label == 'c' or label2 == 'c')
    type(dual_cmplx${bits}$) :: fr
#:else
    type(dual_real${bits}$) :: fr
#:endif

    integer :: order, ii

    order = min(size(A%f), size(B%f)) - 1
    call initialize_dual(fr, order)

    do ii = 0, order
      fr%f(ii) = times_subterm(A, B, ii)
    end do

  end function timesDuals_${label}$${bits}$_${label2}$${bits}$


  !> Dual number divided by another dual number
  elemental function divDuals_${label}$${bits}$_${label2}$${bits}$(A, B) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: A

    !> Dual number
    type(dual_${name2}$${bits}$), intent(in) :: B

    !> Dual number (= A / B)
#:if (label == 'c' or label2 == 'c')
    type(dual_cmplx${bits}$) :: fr
#:else
    type(dual_real${bits}$) :: fr
#:endif

    integer :: order

    order = min(size(A%f), size(B%f)) - 1
    call initialize_dual(fr, order)

    fr = A * inv(B)

  end function divDuals_${label}$${bits}$_${label2}$${bits}$


  !> Exponentiate a dual variable with another
#:if type == "real"
  !! For real values, convert to complex for conversion, as analytical for negative quadrant
#:endif
  elemental function powerDuals_${label}$${bits}$_${label2}$${bits}$(A, B) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: A

    !> Dual number
    type(dual_${name2}$${bits}$), intent(in) :: B

    !> Dual number (= A ** B)
#:if (label == 'c' or label2 == 'c')
    type(dual_cmplx${bits}$) :: fr
#:else
    type(dual_real${bits}$) :: fr
#:endif

    integer :: order

#:if type == "real" and type2 == "real"

    type(dual_cmplx${bits}$) :: cA, cB, cfr

#:endif

    order = min(size(A%f), size(B%f)) - 1
    call initialize_dual(fr, order)

#:if type == "real" and type2 == "real"

    if (isInteger(B)) then

      fr = A**int(B%f(0))

    else

      call initialize_dual(cA, order)
      call initialize_dual(cB, order)
      call initialize_dual(cfr, order)
      cA%f = A%f
      cB%f = B%f

      cfr = exp(log(cA)*cB)
      fr%f(:) = cfr%f%re

    end if

#:else

    if (isInteger(B)) then
#:if type2 == "real"
      fr = A**int(B%f(0))
#:else
      fr = cmplx_${bits}$(A)**int(B%f(0)%re)
#:endif
    else
      fr = exp(log(A)*B)
    end if

#:endif

  end function powerDuals_${label}$${bits}$_${label2}$${bits}$

#:endfor

  !> Determines if a ${name}$ value is equivalent to a real integer
  elemental function isInteger_${label}$${bits}$(x) result(isIt)

    !> ${type}$ value
    ${type}$(real${bits}$), intent(in) :: x

    !> Is the dual variable an integer?
    logical :: isIt

    isIt = (x == int(x))

  end function isInteger_${label}$${bits}$


  !> Determines if a dual value is equivalent to a real integer
  elemental function isInteger_dual${label}$${bits}$(A) result(isIt)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: A

    !> Is the dual variable an integer?
    logical :: isIt

#:if type == "complex"

    isIt = (A%f(0)%re == int(A%f(0)%re) .and. A%f(0)%im == 0 .and. all(A%f(1:) == 0))

#:else

    isIt = (A%f(0) == int(A%f(0)) .and. all(A%f(1:) == 0))

#:endif

  end function isInteger_dual${label}$${bits}$

#:for name2, type2, label2 in FLAVOURS

  !> Sub-terms in product between dual numbers
  pure function times_subterm_${label}$${bits}$_${label2}$${bits}$(A, B, kk) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: A

    !> Dual number
    type(dual_${name2}$${bits}$), intent(in) :: B

    !> Term being evaluated
    integer, intent(in) :: kk

#:if (label == 'c' or label2 == 'c')
    complex(real${bits}$) :: fr
#:else
    real(real${bits}$) :: fr
#:endif

    integer :: ii

    fr = 0
    do ii = 0, kk
      fr = fr + binomial${bits}$(kk, ii) * A%f(ii) * B%f(kk-ii)
    end do

  end function times_subterm_${label}$${bits}$_${label2}$${bits}$

#:endfor

#:for comment, op in [("Inverse", "inv"), ("Exponential", "exp"), ("Logarithm", "log"),&
  & ("sin", "sin"), ("cos", "cos")]

  !> ${comment}$ function of a dual number
  elemental function ${op}$_${label}$${bits}$(g) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: g

    !> Resulting value
    type(dual_${name}$${bits}$) :: fr

    integer :: kk, order

    order = size(g%f) - 1
    call initialize_dual(fr, order)

    do kk = 0, order
      fr%f(kk) = Dnd(${op}$_subterm_${label}$${bits}$, g, kk)
    end do

  end function ${op}$_${label}$${bits}$

#:endfor


  !> Sub-term in the 1/z evaluation
  pure function inv_subterm_${label}$${bits}$(z, order) result(fr)

    !> Number from part of dual value
    ${type}$(real${bits}$), intent(in) :: z

    !> Overall order of the dual value
    integer, intent(in) :: order

    !> Resulting term
    type(dual_${name}$${bits}$) :: fr

    integer :: kk, signo

    call initialize_dual(fr, order)

    signo=1
    do kk = 0, order
      fr%f(kk) = signo * gamma(real(kk+1, kind=real${bits}$)) / (z**(kk+1))
      signo = - signo
    end do

  end function inv_subterm_${label}$${bits}$


  !> Sub-term in the exp evaluation
  pure function exp_subterm_${label}$${bits}$(z, order) result(fr)

    !> Number from part of dual value
    ${type}$(real${bits}$), intent(in) :: z

    !> Overall order of the dual value
    integer, intent(in) :: order

    !> Resulting term
    type(dual_${name}$${bits}$) :: fr

    integer :: kk

    call initialize_dual(fr, order)

    do kk = 0, order
      fr%f(kk) = exp(z)
    end do

  end function exp_subterm_${label}$${bits}$


  !> Sub-term in the log evaluation
  pure function log_subterm_${label}$${bits}$(z, order) result(fr)

    !> Number from part of dual value
    ${type}$(real${bits}$), intent(in) :: z

    !> Overall order of the dual value
    integer, intent(in) :: order

    !> Resulting term
    type(dual_${name}$${bits}$) :: fr

    integer :: kk, signo

    call initialize_dual(fr, order)
    fr%f(0) = log(z)
    signo=1
    do kk = 1, order
      fr%f(kk) = signo * gamma(real(kk, kind=real${bits}$)) / (z**kk)
      signo = - signo
    end do

  end function log_subterm_${label}$${bits}$


  !> Sub-term in the sin evaluation
  pure function sin_subterm_${label}$${bits}$(z, order) result(fr)

    !> Number from part of dual value
    ${type}$(real${bits}$), intent(in) :: z

    !> Overall order of the dual value
    integer, intent(in) :: order

    !> Resulting term
    type(dual_${name}$${bits}$) :: fr

    integer :: kk

    call initialize_dual(fr, order)
    do kk = 0, order
      fr%f(kk) = sin(z + kk * Pi_${bits}$ / 2)
    end do

  end function sin_subterm_${label}$${bits}$


  !> Sub-term in the cos evaluation
  pure function cos_subterm_${label}$${bits}$(z, order) result(fr)

    !> Number from part of dual value
    ${type}$(real${bits}$), intent(in) :: z

    !> Overall order of the dual value
    integer, intent(in) :: order

    !> Resulting term
    type(dual_${name}$${bits}$) :: fr

    integer :: kk

    call initialize_dual(fr, order)
    do kk = 0, order
      fr%f(kk) = cos(z + kk * Pi_${bits}$ / 2)
    end do

  end function cos_subterm_${label}$${bits}$

  !> Square root
  pure function sqrt_${label}$${bits}$(g) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: g

    !> Resulting term
    type(dual_${name}$${bits}$) :: fr

    integer :: order

    order = size(g%f) - 1
    call initialize_dual(fr, order)
#:if type == "complex"
    fr = g ** cmplx(0.5, 0.0, kind=real${bits}$)
#:else
    fr = g ** 0.5_real${bits}$
#:endif

  end function sqrt_${label}$${bits}$

#:if type == "complex"
  !> absx = sqrt(z*z), NOT sqrt(z*conjg(z))
  !! Note: convention matches behavior of the JAX framework, see Table 1 of
  !! https://arxiv.org/abs/2409.06752v3
#:else
  !> absx = sqrt(z*z)
#:endif
  elemental function absx_${label}$${bits}$(g) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: g

    !> Resulting term
    type(dual_${name}$${bits}$) :: fr

    ${type}$(real${bits}$) :: g0
    integer :: order

    order = size(g%f) - 1
    call initialize_dual(fr, order)
    g0 = g%f(0)
    fr%f(:) = g%f * g0 / sqrt(g0*g0)

  end function absx_${label}$${bits}$

  !> Hyperbolic sin function
  elemental function sinh_${label}$${bits}$(g) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: g

    !> Resulting term
    type(dual_${name}$${bits}$) :: fr

    integer :: order

    order = size(g%f) - 1
    call initialize_dual(fr, order)
    fr = (exp(g) - exp(-g))
    fr%f(:) = 0.5_real${bits}$ * fr%f

  end function sinh_${label}$${bits}$

  !> Hyperbolic cosine function
  elemental function cosh_${label}$${bits}$(g) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: g

    !> Resulting term
    type(dual_${name}$${bits}$) :: fr

    integer :: order

    order = size(g%f) - 1
    call initialize_dual(fr, order)
    fr = (exp(g) + exp(-g))
    fr%f(:) = 0.5_real${bits}$ * fr%f

  end function cosh_${label}$${bits}$

  !> Hyperbolic tangent function
  elemental function tanh_${label}$${bits}$(g) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: g

    !> Resulting term
    type(dual_${name}$${bits}$) :: fr

    integer :: order

    order = size(g%f) - 1
    call initialize_dual(fr, order)
    fr = (exp(g) - exp(-g)) / (exp(g) + exp(-g))

  end function tanh_${label}$${bits}$


    !> tangent function
  elemental function tan_${label}$${bits}$(g) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: g

    !> Resulting term
    type(dual_${name}$${bits}$) :: fr

    integer :: order

    order = size(g%f) - 1
    call initialize_dual(fr, order)
    fr = sin(g) * inv(cos(g))

  end function tan_${label}$${bits}$

  !> asin function of a dual number
  elemental function asin_${label}$${bits}$(g) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: g

    !> Resulting value
    type(dual_${name}$${bits}$) :: fr

    integer :: order, kk

    order = size(g%f) - 1
    call initialize_dual(fr, order)

    do kk = 0, order
      fr%f(kk) = Dnd(asin_subterm_${label}$${bits}$, g, kk)
    end do

  end function asin_${label}$${bits}$

  !> Sub-term in the asin evaluation
  pure function asin_subterm_${label}$${bits}$(z, order) result(fr)

    !> Number from part of dual value
    ${type}$(real${bits}$), intent(in) :: z

    !> Overall order of the dual value
    integer, intent(in) :: order

    !> Resulting term
    type(dual_${name}$${bits}$) :: fr

    integer :: ii, kk
    ${type}$(real${bits}$) :: sumterm, sqz
    type(dual_${name}$${bits}$) :: sqzdn, zdn

    ${type}$(real${bits}$), parameter :: zero = 0.0_real${bits}$
    ${type}$(real${bits}$), parameter :: one = 1.0_real${bits}$

    call initialize_dual(fr, order)
    fr%f(0) = asin(z)
    if (order == 0) return
    fr%f(1) = one / sqrt(one - z**2)
    if (order == 1) return

    call initialize_dual(zdn, max(order, 1))
    zdn%f(0) = z
    zdn%f(1) = one

    sqz = sqrt(one - z * z)
    call initialize_dual(sqzdn, max(order, 1))
    sqzdn%f(0) = one
    sqzdn = sqrt(sqzdn - zdn * zdn)

    do kk = 2, order
       sumterm = zero
       do ii = 1, kk - 1
          sumterm = sumterm + binomial${bits}$(kk-1, ii) * sqzdn%f(ii) * fr%f(kk-ii) / sqz
       end do
       fr%f(kk) = -sumterm
     end do

   end function asin_subterm_${label}$${bits}$

   !> acos function of a dual number
  elemental function acos_${label}$${bits}$(g) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: g

    !> Resulting value
    type(dual_${name}$${bits}$) :: fr

    integer :: order, kk

    order = size(g%f) - 1
    call initialize_dual(fr, order)

    do kk = 0, order
      fr%f(kk) = Dnd(acos_subterm_${label}$${bits}$, g, kk)
    end do

  end function acos_${label}$${bits}$

  !> Sub-term in the acos evaluation
  pure function acos_subterm_${label}$${bits}$(z, order) result(fr)

    !> Number from part of dual value
    ${type}$(real${bits}$), intent(in) :: z

    !> Overall order of the dual value
    integer, intent(in) :: order

    !> Resulting term
    type(dual_${name}$${bits}$) :: fr

    integer :: ii, kk
    ${type}$(real${bits}$) :: sumterm, sqz
    type(dual_${name}$${bits}$) :: sqzdn, zdn

    ${type}$(real${bits}$), parameter :: zero = 0.0_real${bits}$
    ${type}$(real${bits}$), parameter :: one = 1.0_real${bits}$

    call initialize_dual(fr, order)
    fr%f(0) = acos(z)
    if (order == 0) return
    fr%f(1) = -one / sqrt(one - z**2)
    if (order == 1) return

    call initialize_dual(zdn, max(order, 1))
    zdn%f(0) = z
    zdn%f(1) = one

    sqz = sqrt(one - z * z)
    call initialize_dual(sqzdn, max(order, 1))
    sqzdn%f(0) = one
    sqzdn = sqrt(sqzdn - zdn * zdn)

    do kk = 2, order
       sumterm = zero
       do ii = 1, kk - 1
          sumterm = sumterm + binomial${bits}$(kk-1, ii) * sqzdn%f(ii) * fr%f(kk-ii) / sqz
       end do
       fr%f(kk) = -sumterm
     end do

   end function acos_subterm_${label}$${bits}$

   !> atan function of a dual number
   elemental function atan_${label}$${bits}$(g) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: g

    !> Resulting value
    type(dual_${name}$${bits}$) :: fr

    integer :: order, kk

    order = size(g%f) - 1
    call initialize_dual(fr, order)

    do kk = 0, order
      fr%f(kk) = Dnd(atan_subterm_${label}$${bits}$, g, kk)
    end do

  end function atan_${label}$${bits}$

  !> Sub-term in the atan evaluation
  pure function atan_subterm_${label}$${bits}$(z, order) result(fr)

    !> Number from part of dual value
    ${type}$(real${bits}$), intent(in) :: z

    !> Overall order of the dual value
    integer, intent(in) :: order

    !> Resulting term
    type(dual_${name}$${bits}$) :: fr

    integer :: ii, kk
    ${type}$(real${bits}$) :: sumterm, den
    type(dual_${name}$${bits}$) :: auxdn

    ${type}$(real${bits}$), parameter :: zero = 0.0_real${bits}$
    ${type}$(real${bits}$), parameter :: one = 1.0_real${bits}$
    ${type}$(real${bits}$), parameter :: two = 2.0_real${bits}$

    call initialize_dual(fr, order)

    fr%f(0) = atan(z)
    if (order == 0) return

    fr%f(1) = one / (one + z**2)
    if (order == 1) return
    den = one + z * z

    call initialize_dual(auxdn, order)
    auxdn%f(0) = den
    auxdn%f(1) = two * z
    auxdn%f(2) = two

    do kk = 2, order
      sumterm = zero
      do ii = 1, kk - 1
        sumterm = sumterm + binomial${bits}$(kk-1, ii) * auxdn%f(ii) * fr%f(kk-ii) / den
      end do
      fr%f(kk) = -sumterm
    end do

  end function atan_subterm_${label}$${bits}$

  !> asinh function of a dual number
  elemental function asinh_${label}$${bits}$(g) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: g

    !> Resulting value
    type(dual_${name}$${bits}$) :: fr

    integer :: order, kk

    order = size(g%f) - 1
    call initialize_dual(fr, order)

    do kk = 0, order
      fr%f(kk) = Dnd(asinh_subterm_${label}$${bits}$, g, kk)
    end do

  end function asinh_${label}$${bits}$

  !> Sub-term in the asinh evaluation
  pure function asinh_subterm_${label}$${bits}$(z, order) result(fr)

    !> Number from part of dual value
    ${type}$(real${bits}$), intent(in) :: z

    !> Overall order of the dual value
    integer, intent(in) :: order

    !> Resulting term
    type(dual_${name}$${bits}$) :: fr

    integer :: ii, kk
    ${type}$(real${bits}$) :: sumterm, den
    type(dual_${name}$${bits}$) :: auxdn, zdn

    ${type}$(real${bits}$), parameter :: zero = 0.0_real${bits}$
    ${type}$(real${bits}$), parameter :: one = 1.0_real${bits}$
    ${type}$(real${bits}$), parameter :: two = 2.0_real${bits}$

    call initialize_dual(fr, order)

    fr%f(0) = asinh(z)
    if (order == 0) return

    den = sqrt(one + z * z)
    fr%f(1) = one / den
    if (order == 1) return

    call initialize_dual(auxdn, order)
    call initialize_dual(zdn, order)

    zdn%f(0) = z
    zdn%f(1) = one

    auxdn = sqrt(one + zdn * zdn)

    do kk = 2, order
      sumterm = zero
      do ii = 1, kk - 1
        sumterm = sumterm + binomial${bits}$(kk-1, ii) * auxdn%f(ii) * fr%f(kk-ii) / den
      end do
      fr%f(kk) = -sumterm
    end do

  end function asinh_subterm_${label}$${bits}$

  !> acosh function of a dual number
  elemental function acosh_${label}$${bits}$(g) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: g

    !> Resulting value
    type(dual_${name}$${bits}$) :: fr

    integer :: order, kk

    order = size(g%f) - 1
    call initialize_dual(fr, order)

    do kk = 0, order
      fr%f(kk) = Dnd(acosh_subterm_${label}$${bits}$, g, kk)
    end do

  end function acosh_${label}$${bits}$

  !> Sub-term in the acosh evaluation
  pure function acosh_subterm_${label}$${bits}$(z, order) result(fr)

    !> Number from part of dual value
    ${type}$(real${bits}$), intent(in) :: z

    !> Overall order of the dual value
    integer, intent(in) :: order

    !> Resulting term
    type(dual_${name}$${bits}$) :: fr

    integer :: ii, kk
    ${type}$(real${bits}$) :: sumterm, den
    type(dual_${name}$${bits}$) :: auxdn, zdn

    ${type}$(real${bits}$), parameter :: zero = 0.0_real${bits}$
    ${type}$(real${bits}$), parameter :: one = 1.0_real${bits}$
    ${type}$(real${bits}$), parameter :: two = 2.0_real${bits}$

    call initialize_dual(fr, order)

    fr%f(0) = acosh(z)
    if (order == 0) return

    den = sqrt(z - one) * sqrt(one + z)
    fr%f(1) = one / den
    if (order == 1) return

    call initialize_dual(auxdn, order)
    call initialize_dual(zdn, order)

    zdn%f(0) = z
    zdn%f(1) = one

    auxdn = sqrt(zdn - one) * sqrt(zdn + one)
    do kk = 2, order
      sumterm = zero
      do ii = 1, kk - 1
        sumterm = sumterm + binomial${bits}$(kk-1, ii) * auxdn%f(ii) * fr%f(kk-ii) / den
      end do
      fr%f(kk) = -sumterm
    end do

  end function acosh_subterm_${label}$${bits}$

  !> atanh function of a dual number
  elemental function atanh_${label}$${bits}$(g) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: g

    !> Resulting value
    type(dual_${name}$${bits}$) :: fr

    integer :: order, kk

    order = size(g%f) - 1
    call initialize_dual(fr, order)

    do kk = 0, order
      fr%f(kk) = Dnd(atanh_subterm_${label}$${bits}$, g, kk)
    end do

  end function atanh_${label}$${bits}$

  !> Sub-term in the atanh evaluation
  pure function atanh_subterm_${label}$${bits}$(z, order) result(fr)

    !> Number from part of dual value
    ${type}$(real${bits}$), intent(in) :: z

    !> Overall order of the dual value
    integer, intent(in) :: order

    !> Resulting term
    type(dual_${name}$${bits}$) :: fr

    integer :: ii, kk
    ${type}$(real${bits}$) :: sumterm, den
    type(dual_${name}$${bits}$) :: auxdn

    ${type}$(real${bits}$), parameter :: zero = 0.0_real${bits}$
    ${type}$(real${bits}$), parameter :: one = 1.0_real${bits}$
    ${type}$(real${bits}$), parameter :: two = 2.0_real${bits}$

    call initialize_dual(fr, order)

    fr%f(0) = atanh(z)
    if (order == 0) return

    ! do not 'simplify' they are complex
    den = one - z*z
    fr%f(1) = one / den
    if (order == 1) return

    call initialize_dual(auxdn, order)
    auxdn%f(0) = den
    auxdn%f(1) = -two * z
    auxdn%f(2) = -two

    do kk = 2, order
      sumterm = zero
      do ii = 1, kk - 1
        sumterm = sumterm + binomial${bits}$(kk-1, ii) * auxdn%f(ii) * fr%f(kk-ii) / den
      end do
      fr%f(kk) = -sumterm
    end do

  end function atanh_subterm_${label}$${bits}$

  !> atan2 function of a dual number
  elemental function atan2_${label}$${bits}$(y, x) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: y

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: x

    !> Resulting value
    type(dual_${name}$${bits}$) :: fr

    integer :: order
    ${type}$(real${bits}$) :: x0, y0

    order = min(size(x%f), size(y%f)) - 1
    call initialize_dual(fr, order)

    fr = atan(y / x)
    x0 = x%f(0)
    y0 = x%f(0)
#:if type == "complex"
    fr%f(0) = atan2_z(y0,x0)
#:else
    fr%f(0) = atan2(y0, x0)
#:endif

  end function atan2_${label}$${bits}$

#:if type == "complex"
  !> Complex atan2 evaluation
  pure function atan2_z_${bits}$(zy, zx) result(f_res)

    !> y value for working out quadrant aware atan
    complex(real${bits}$), intent(in) :: zy

    !> x value for working out quadrant aware atan
    complex(real${bits}$), intent(in) :: zx

    !> Resulting value
    complex(real${bits}$) :: f_res

    complex(real${bits}$), parameter :: im = cmplx(0, 1, kind=real${bits}$)
    complex(real${bits}$) :: num, den, divnd, t1, t2
    real(real${bits}$) :: x1, x2, y1, y2
    complex(real${bits}$) :: x1c, x2c, y1c, y2c

    x1 = zx%re
    x2 = zx%im

    y1 = zy%re
    y2 = zy%im

    x1c = x1
    x2c = x2
    y1c = y1
    y2c = y2

    num = x1c + im*x2c + im*(y1c + im*y2c)
    den = sqrt((x1c + im*x2c)**2 + (y1c + im*y2c)**2)
    divnd = num/den;
    t1 = atan2(divnd%im,divnd%re)
    t2 = im * log( sqrt((x2c + y1c)**2 + (x1c - y2c)**2)&
        & / ( (2.0_real${bits}$*x1c* x2c +  2.0_real${bits}$*y1c*y2c)**2&
        & + (x1c**2 - x2c**2 + y1c**2 - y2c**2)**2 )**0.25_real${bits}$ )

    f_res = t1 - t2

  end function atan2_z_${bits}$
#:endif


#:for op in [("sin"), ("cos"), ("tan")]

  !> ${op}$ function of a dual number in degrees
  elemental function ${op}$d_${label}$${bits}$(g) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: g

    !> Resulting value
    type(dual_${name}$${bits}$) :: fr

    integer :: order

    order = size(g%f) - 1
    call initialize_dual(fr, order)

#:if type == "complex"
    fr = ${op}$(g * cmplx(Pi_${bits}$ / 180.0_real${bits}$, 0.0, kind=real${bits}$))
#:else
    fr = ${op}$(Pi_${bits}$ * g / 180.0_real${bits}$)
#:endif

  end function ${op}$d_${label}$${bits}$

#:endfor

#:for op in [("asin"), ("acos"), ("atan")]

  !> ${op}$ function of a dual number, returning in degrees
  elemental function ${op}$d_${label}$${bits}$(g) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: g

    !> Resulting value
    type(dual_${name}$${bits}$) :: fr

    integer :: order

    order = size(g%f) - 1
    call initialize_dual(fr, order)

    fr = ${op}$(g)

#:if type == "complex"
    fr = fr * cmplx(Pi_${bits}$ / 180.0_real${bits}$, 0.0, kind=real${bits}$)
#:else
    fr = fr * Pi_${bits}$ / 180.0_real${bits}$
#:endif

  end function ${op}$d_${label}$${bits}$

#:endfor

  !> log10 function of a dual number
  elemental function log10_${label}$${bits}$(g) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: g

    !> Resulting value
    type(dual_${name}$${bits}$) :: fr

    integer :: order

    order = size(g%f) - 1
    call initialize_dual(fr, order)

#:if type == "complex"
    fr = log(g) / log(cmplx(10.0_real${bits}$, 0.0, kind=real${bits}$))
#:else
    fr = log(g) / log(10.0_real${bits}$)
#:endif

  end function log10_${label}$${bits}$

#:endfor
#:endfor


#:for bits in REALMODEL
  !> Conjg for (df)* NOT d(f*)
  elemental function conjg_${bits}$(g) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: g

    !> Resulting dual number
    type(dual_cmplx${bits}$) :: fr

    integer :: order

    order = size(g%f) - 1
    call initialize_dual(fr, order)
    fr%f(:) = conjg(g%f)

  end function conjg_${bits}$

#:endfor

end module src_arithmeticfunctions
