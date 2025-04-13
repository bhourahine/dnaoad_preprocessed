!> Dual numbers for automatic differentiation up to order n
!! This can be used to compute first, second, ... nth order derivatives
!! F. Pe~nu~nuri
!! UADY, Merida Yucatan Mexico, 2024

#:include 'common.fypp'

module dnaoad

  use src_chainrule, only : Dnd

  use src_arithmeticfunctions, only : inv, sin, cos, tan, exp, log, sqrt, asin, acos, atan, asinh
  use src_arithmeticfunctions, only : acosh, atanh, sinh, cosh, tanh, absx, atan2
  use src_arithmeticfunctions, only : conjg
  use src_arithmeticfunctions, only : sind, cosd, tand, asind, acosd, atand
  use src_arithmeticfunctions, only : log10

  use src_matrixfunctions, only : sum, product, matmul

  use src_arithmeticfunctions, only : assignment(=)
  use src_arithmeticfunctions, only : operator(*), operator(/), operator(+), operator(-)
  use src_arithmeticfunctions, only : operator(**)
  use src_logicalfunctions, only : operator(==), operator(/=)

  use src_arithmeticfunctions, only : isInteger
  use src_helperfunctions, only : f_part
  use src_types
  use src_conversions

  implicit none

  public

end module dnaoad
