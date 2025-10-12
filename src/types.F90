!> Dual numbers for automatic differentiation up to order n
!! This can be used to compute first, second, ... nth order derivatives
!! F. Pe~nu~nuri
!! UADY, Merida Yucatan Mexico, 2024

#:include 'common.fypp'

!> Data types for dual values
module src_types
  use iso_fortran_env, only : ${ISO_CONSTANTS}$
  use src_constants
  implicit none
  private
#:for name, type, label in FLAVOURS
#:for bits in REALMODEL
  public :: dual_${name}$${bits}$
#:endfor
#:endfor
  public :: initialize_dual
  public :: assignment (=)
  public :: operator(==), operator(/=)
  public :: operator(+), operator(-), operator(*), operator(/), operator(**)
  public :: inv, sin, cos, tan, exp, log, sqrt, asin, acos, atan, asinh, acosh
  public :: atanh, sinh, cosh, tanh, absx, atan2
  public :: sind, cosd, tand, asind, acosd, atand, log10
  public :: conjg
  public :: isInteger
  public :: aimag, f_part
#:for bits in REALMODEL
  public :: real_${bits}$, cmplx_${bits}$
#:endfor
  public :: Dnd

  ! To do: mixed precision for operators, bessel_*, erf, erfc, erfc_scaled, gamma, log_gamma

#:for name, type, label in FLAVOURS

#:for bits in REALMODEL
  !> Dual ${type}$ variable type for ${bits}$ precision
  type :: dual_${name}$${bits}$
    ${type}$(kind=real${bits}$), private, allocatable :: f(:)
  contains
    procedure :: order => dual_${name}$${bits}$_order
    procedure :: get_derivative => f_part_${label}$${bits}$
    procedure :: get_derivatives => f_${label}$${bits}$
    procedure :: set_derivative => setDeriv_${label}$${bits}$
    procedure :: set_derivatives => setDerivs_${label}$${bits}$
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

  !> Assignment operator
  interface assignment (=)
#:for _, _, label in FLAVOURS
#:for bits in REALMODEL
#:for intbits in INTMODEL
    module procedure assign_int${intbits}$_to_dual${label}$${bits}$
    module procedure assign_ints${intbits}$_to_dual${label}$${bits}$
#:endfor
    module procedure assign_${label}$${bits}$_to_dual${label}$${bits}$
    module procedure assign_${label}$${bits}$s_to_dual${label}$${bits}$
    module procedure assign_dual${label}$${bits}$_to_dual${label}$${bits}$
#:endfor
#:endfor
  end interface assignment (=)

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

  !> Driver for differentiation of a function
  interface Dnd
#:for _, _, label in FLAVOURS
#:for bits in REALMODEL
    module procedure Dnd_${label}$${bits}$
#:endfor
#:endfor
  end interface Dnd

#:for bits in REALMODEL

#:for name, type, label in FLAVOURS
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

  !> Interface for conversion to a dual number of type real_${bits}$
  interface real_${bits}$
#:for _, _, label in FLAVOURS
#:for inbits in REALMODEL
    module procedure real${bits}$_${label}$${inbits}$
#:endfor
#:endfor
  end interface real_${bits}$

  !> Interface for conversion to a dual number of type cmplx_${bits}$
  interface cmplx_${bits}$
#:for inbits in REALMODEL
    module procedure cmplx${bits}$_c${inbits}$
    module procedure cmplx${bits}$_r${inbits}$
#:for bitsopt in REALMODEL
    module procedure cmplx${bits}$_r${inbits}$_opt${bitsopt}$
#:endfor
#:endfor
  end interface cmplx_${bits}$

#:endfor


  !> Interface for extracting imaginary part of a complex dual number as a real dual
  interface aimag
#:for bits in REALMODEL
    module procedure aimag_c${bits}$
#:endfor
  end interface aimag

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

  !> Converts a complecx ${bits}$-bit dual value to a real ${bits}$-bit dual
  elemental function aimag_c${bits}$(A) result(B)

    !> Dual number
    type(dual_cmplx${bits}$), intent(in) :: A

    !> Dual number
    type(dual_real${bits}$) :: B

    integer :: order

    order = A%order()
    call initialize_dual(B, order)

    B%f(:) = aimag(A%f)

  end function aimag_c${bits}$

#:for inbits in REALMODEL

  !> Converts a complex ${inbits}$-bit dual value to a complex ${bits}$-bit dual
  elemental function cmplx${bits}$_c${inbits}$(A) result(B)

    !> Dual number
    type(dual_cmplx${inbits}$), intent(in) :: A

    !> Dual number
    type(dual_cmplx${bits}$) :: B

    integer :: order

    order = A%order()
    call initialize_dual(B, order)

    B%f(:) = cmplx(A%f, kind=real${bits}$)

  end function cmplx${bits}$_c${inbits}$

  !> Converts a real ${inbits}$-bit dual value to a complex ${bits}$-bit dual
  elemental function cmplx${bits}$_r${inbits}$(A) result(C)

    !> Dual number for real part
    type(dual_real${inbits}$), intent(in) :: A

    !> Dual number
    type(dual_cmplx${bits}$) :: C

    integer :: order

    order = A%order()
    call initialize_dual(C, order)

    C%f(:) = cmplx(A%f, kind=real${bits}$)

  end function cmplx${bits}$_r${inbits}$

#:for bitsopt in REALMODEL

  !> Converts a real ${inbits}$-bit dual value to a complex ${bits}$-bit dual
  elemental function cmplx${bits}$_r${inbits}$_opt${bitsopt}$(A, B) result(C)

    !> Dual number for real part
    type(dual_real${inbits}$), intent(in) :: A

    !> Dual number for imaginary part
    type(dual_real${bitsopt}$), intent(in) :: B

    !> Dual number
    type(dual_cmplx${bits}$) :: C

    integer :: order

    order = min(A%order(), B%order())
    call initialize_dual(C, order)

    C%f(:) = cmplx(A%f(0:order), B%f(0:order), kind=real${bits}$)

  end function cmplx${bits}$_r${inbits}$_opt${bitsopt}$

#:endfor

#:endfor

#:for name, type, label in FLAVOURS

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

  pure function dual_${name}$${bits}$_order(this) result(order_value)

    class(dual_${name}$${bits}$), intent(in) :: this

    integer :: order_value

    order_value = size(this%f) - 1

  end function dual_${name}$${bits}$_order

  !> Extract the k-th derivative part of a dual number or an array of duals
  elemental function f_part_${label}$${bits}$(A, k) result(fr)

    !> Dual number
    class(dual_${name}$${bits}$), intent(in) :: A

    !> Derivative order to extract
    integer, intent(in) :: k

    !> Resulting ${type}$ value
    ${type}$(real${bits}$) :: fr

    fr = A%f(k)

  end function f_part_${label}$${bits}$

  !> Extract the derivatives part of a dual number or an array of duals
  pure function f_${label}$${bits}$(A) result(fr)

    !> Dual number
    class(dual_${name}$${bits}$), intent(in) :: A

    !> Resulting ${type}$ value
    ${type}$(real${bits}$) :: fr(size(A%f))

    fr(:) = A%f

  end function f_${label}$${bits}$

  !> Sets the k-th derivative part of a dual number
  pure subroutine setDeriv_${label}$${bits}$(A, k, fr)

    !> Dual number
    class(dual_${name}$${bits}$), intent(inout) :: A

    !> Derivative order
    integer, intent(in) :: k

    !> ${type}$ value
    ${type}$(real${bits}$), intent(in) :: fr

    A%f(k) = fr

  end subroutine setDeriv_${label}$${bits}$

  !> Sets all orders of the derivative parts of a dual number
  pure subroutine setDerivs_${label}$${bits}$(A, fr)

    !> Dual number
    class(dual_${name}$${bits}$), intent(inout) :: A

    !> ${type}$ value
    ${type}$(real${bits}$), intent(in) :: fr(:)

    A%f(:) = fr

  end subroutine setDerivs_${label}$${bits}$

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

  !> Assign an array of int${intbits}$ values to an already initialized dual number of type
  !! ${name}$${bits}$
  !! Note, the order of the dual number is retained in the operation
  subroutine assign_ints${intbits}$_to_dual${label}$${bits}$(A, ii)

    !> Resulting dual number
    type(dual_${name}$${bits}$), intent(inout) :: A

    !> Integer value to copy into the dual number
    integer(int${intbits}$), intent(in) :: ii(0:)

    integer :: order

    order = min(size(ii) - 1, A%order())
    A%f(:) = ${name}$(0)
    A%f(0:order) = ${name}$(ii(0:order))

  end subroutine assign_ints${intbits}$_to_dual${label}$${bits}$

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

  !> Assign an array of ${name}$${bits}$ values to an already initialized dual number of type
  !! ${name}$${bits}$ as its derivatives.
  !! Note, the order of the output dual number is retained in the operation
  subroutine assign_${label}$${bits}$s_to_dual${label}$${bits}$(A, x)

    !> Resulting dual number
    type(dual_${name}$${bits}$), intent(inout) :: A

    !> Values to copy into the dual number
    ${type}$(real${bits}$), intent(in) :: x(0:)

    integer :: order

    order = min(size(x) - 1, A%order())

    A%f(:) = 0
    A%f(0:order) = x(0:order)

  end subroutine assign_${label}$${bits}$s_to_dual${label}$${bits}$

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

  !> Check equality between two dual numbers of type ${label}$${bits}$
  elemental function equal_dual${label}$${bits}$(lhs, rhs) result (eqv)

    !> First dual number
    type(dual_${name}$${bits}$), intent(in) :: lhs

    !> Second dual number
    type(dual_${name}$${bits}$), intent(in) :: rhs

    !> Are these numbers equivalent
    logical :: eqv

    eqv = lhs%order() == rhs%order()
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

    neqv = lhs%order() /= rhs%order()
    if (.not. neqv) neqv = any(lhs%f /= rhs%f)

  end function notequal_dual${label}$${bits}$

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

    order = gdual%order()
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

#:for inbits in REALMODEL

  !> Converts a ${type}$ ${inbits}$-bit dual value to a real ${inbits}$-bit dual
  elemental function real${bits}$_${label}$${inbits}$(A) result(B)

    !> Dual number
    type(dual_${name}$${inbits}$), intent(in) :: A

    !> Dual number
    type(dual_real${bits}$) :: B

    integer :: order

    order = A%order()
    call initialize_dual(B, order)

    B%f(:) = real(A%f, kind=real${bits}$)

  end function real${bits}$_${label}$${inbits}$

#:endfor

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

    call initialize_dual(fr, A%order())

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

    order = A%order()
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

    order = A%order()
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

    order = A%order()
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

    order = A%order()
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

    order = A%order()
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

    order = A%order()
    call initialize_dual(tmpDual, order)
    call initialize_dual(fr, order)

    tmpDual = 1
    if (A == tmpDual) then
      fr = A
      return
    end if

    tmpDual = 0
    if (n > 0 .and. A == tmpDual) then
      fr = A
      return
    end if

    fr = 1 ! also covers 0**0 ---> 1
    if (n /= 0) then
      do kk = 1, abs(n)
        fr = fr * A
      end do
      if (n < 0) then
        fr = inv(fr)
      end if
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

    order = B%order()
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

    order = A%order()
    call initialize_dual(fr, order)

    if (isInteger(x)) then

    #:if type == "complex"
      fr = A ** nint(x%re)
    #:else
      fr = A ** nint(x)
    #:endif

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

    order = min(A%order(), B%order())
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

    order = min(A%order(), B%order())
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

    order = min(A%order(), B%order())
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

    order = min(A%order(), B%order())
    call initialize_dual(fr, order)

#:if type == "real" and type2 == "real"

    if (isInteger(B)) then

      fr = A ** nint(B%f(0))

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
      fr = A ** nint(B%f(0))
#:else
      fr = cmplx_${bits}$(A) ** nint(B%f(0)%re)
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

    !> Is the dual variable a real valued integer?
    logical :: isIt

  #:if type == "complex"
    isIt = (x == nint(x%re))
  #:else
    isIt = (x == nint(x))
  #:endif

  end function isInteger_${label}$${bits}$


  !> Determines if a dual value is equivalent to a real integer
  elemental function isInteger_dual${label}$${bits}$(A) result(isIt)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: A

    !> Is the dual variable an integer?
    logical :: isIt

#:if type == "complex"

    isIt = (A%f(0)%re == nint(A%f(0)%re) .and. A%f(0)%im == 0 .and. all(A%f(1:) == 0))

#:else

    isIt = (A%f(0) == nint(A%f(0)) .and. all(A%f(1:) == 0))

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

    order = g%order()
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

    order = g%order()
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

    order = g%order()
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

    order = g%order()
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

    order = g%order()
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

    order = g%order()
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

    order = g%order()
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

    order = g%order()
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

    order = g%order()
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

    order = g%order()
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

    order = g%order()
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

    order = g%order()
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

    order = g%order()
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

    order = min(x%order(), y%order())
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


#:for op in [("sind"), ("cosd"), ("tand")]
#:set func = op[:-1]

  !> ${func}$ function of a dual number in degrees
  elemental function ${op}$_${label}$${bits}$(g) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: g

    !> Resulting value
    type(dual_${name}$${bits}$) :: fr

    integer :: order

    order = g%order()
    call initialize_dual(fr, order)

#:if type == "complex"
    fr = ${func}$(g * cmplx(Pi_${bits}$ / 180.0_real${bits}$, 0.0, kind=real${bits}$))
#:else
    fr = ${func}$(Pi_${bits}$ * g / 180.0_real${bits}$)
#:endif

  end function ${op}$_${label}$${bits}$

#:endfor

#:for op in [("asind"), ("acosd"), ("atand")]
#:set func = op[:-1]

  !> ${func}$ function of a dual number, returning in degrees
  elemental function ${op}$_${label}$${bits}$(g) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: g

    !> Resulting value
    type(dual_${name}$${bits}$) :: fr

    integer :: order

    order = g%order()
    call initialize_dual(fr, order)

    fr = ${func}$(g)

#:if type == "complex"
    fr = fr * cmplx(Pi_${bits}$ / 180.0_real${bits}$, 0.0, kind=real${bits}$)
#:else
    fr = fr * Pi_${bits}$ / 180.0_real${bits}$
#:endif

  end function ${op}$_${label}$${bits}$

#:endfor

  !> log10 function of a dual number
  elemental function log10_${label}$${bits}$(g) result(fr)

    !> Dual number
    type(dual_${name}$${bits}$), intent(in) :: g

    !> Resulting value
    type(dual_${name}$${bits}$) :: fr

    integer :: order

    order = g%order()
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

    order = g%order()
    call initialize_dual(fr, order)
    fr%f(:) = conjg(g%f)

  end function conjg_${bits}$

#:endfor

end module src_types
