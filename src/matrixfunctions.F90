!> Dual numbers for automatic differentiation up to order n
!! This can be used to compute first, second, ... nth order derivatives
!! F. Pe~nu~nuri
!! UADY, Merida Yucatan Mexico, 2024

! To do: dot_product, BLAS and ScaLAPACK bindings for selected operations

#:include 'common.fypp'

#:set DIMS = [("0", ""), ("1", "(:)"), ("2", "(:,:)"), ("3", "(:,:,:)"), ("4", "(:,:,:,:)")]
#:set INDICES = [('ii'),('jj'),('kk'),('ll')]

!> Functions for matrices of dual numbers
module src_matrixfunctions
  use iso_fortran_env, only : ${ISO_CONSTANTS}$
  use src_arithmeticfunctions, only : assignment(=), operator(+), operator(*)
  use src_types

  implicit none
  private
  public :: sum, product, matmul

  ! To do: transpose, dot_product, unpack, norm2, optional GEMM and GEMV (note is_contiguous test
  ! required)

#:for subname, descript in [("sum", "sum"), ("prod", "product")]

  !> The ${descript}$ of elements of a matrix of dual numbers
  interface ${descript}$
#:for _, _, label in FLAVOURS
#:for bits in REALMODEL
#:for dim, _ in DIMS[1:]
    module procedure ${subname}$${dim}$D_${label}$${bits}$
#:endfor
#:for dim, _ in DIMS[2:]
    module procedure ${subname}$${dim}$D_dim_${label}$${bits}$
#:endfor
#:endfor
#:endfor
  end interface ${descript}$

#:endfor


  !> Matrix multiplication
  interface matmul
#:for _, _, label in FLAVOURS
#:for bits in REALMODEL
#:for intbits in INTMODEL
    module procedure matmul_int${intbits}$_with_dual${label}$${bits}$
    module procedure matmul_dual${label}$${bits}$_with_int${intbits}$
#:endfor
    module procedure matmul_${label}$${bits}$_with_dual${label}$${bits}$
    module procedure matmul_dual${label}$${bits}$_with_${label}$${bits}$
    module procedure matmulDuals_${label}$${bits}$
#:endfor
#:endfor
  end interface matmul

contains

#:for subname, descript, op in [("sum", "sum", "+"), ("prod", "product", "*")]

#:for dim, indices in DIMS[1:]
#:set firstelement = indices.replace(':', '1')
#:set iCount = int(dim)
#:set contrCount = int(dim) - 1
#:set loopindices = ', '.join(INDICES[:iCount])
#:set loopindices2 = ', '.join(["indx("+str(i)+")" for i in range(1,iCount)])
#:set contractedout = DIMS[iCount-1][1]

#:for name, type, label in FLAVOURS
#:for bits in REALMODEL

  ! pairwise sum ...

  !> The ${descript}$() operation for ${type}$ matrices of dual numbers with ${bits}$ bits
  function ${subname}$${dim}$D_${label}$${bits}$(A, mask) result(fr)

    !> Matrix of dual numbers
    type(dual_${name}$${bits}$), intent(in) :: A${indices}$

    !> Optional logical mask of elements
    logical, intent(in), optional :: mask${indices}$

    !> Resulting ${descript}$
    type(dual_${name}$${bits}$) :: fr

    integer :: order, ${loopindices}$
    ${type}$(real${bits}$), parameter :: one = 1.0_real${bits}$

    order = size(A${firstelement}$%f) - 1
    call initialize_dual(fr, order)
#:if subname == "prod"
    fr = one
#:endif

    if (present(mask)) then

#:for ind in range(iCount, 0, -1)
#:set jCounter = INDICES[ind-1]
#:set ws = "  " * (iCount-ind)
      ${ws}$do ${jCounter}$ = 1, size(A, dim=${ind}$)
#:endfor
      ${ws}$  if (mask(${loopindices}$)) fr = fr ${op}$ A(${loopindices}$)
#:for ind in range(iCount, 0, -1)
#:set ws = "  "*(ind-1)
      ${ws}$enddo
#:endfor

    else

#:for ind in range(iCount, 0, -1)
#:set jCounter = INDICES[ind-1]
#:set ws = "  " * (iCount-ind)
      ${ws}$do ${jCounter}$ = 1, size(A, dim=${ind}$)
#:endfor
      ${ws}$  fr = fr ${op}$ A(${loopindices}$)
#:for ind in range(iCount, 0, -1)
#:set ws = "  "*(ind-1)
      ${ws}$enddo
#:endfor

    end if

  end function ${subname}$${dim}$D_${label}$${bits}$

#:if contrCount > 0

  !> The ${descript}$() operation for ${type}$ matrices of dual numbers with ${bits}$ bits
  function ${subname}$${dim}$D_dim_${label}$${bits}$(A, dim, mask) result(fr)

    !> Matrix of dual numbers
    type(dual_${name}$${bits}$), intent(in) :: A${indices}$

    !> Dimension along which to sum
    integer, intent(in) :: dim

    !> Optional logical mask of elements
    logical, intent(in), optional :: mask${indices}$

    !> Resulting ${descript}$
    type(dual_${name}$${bits}$), allocatable :: fr${contractedout}$

    integer :: order, indx(${contrCount}$)
    integer :: ${loopindices}$
    type(dual_${name}$${bits}$), allocatable :: frTmp(:)
    ${type}$(real${bits}$), parameter :: one = 1.0_real${bits}$

    order = size(A${firstelement}$%f) - 1
    indx(:) = skipDim(dim,shape(A))
    allocate(frTmp(product(indx)))
    fr = reshape(frTmp, indx)
    deallocate(frTmp)
    call initialize_dual(fr, order)
#:if subname == "prod"
    fr${contractedout}$ = one
#:endif

    if (present(mask)) then

#:for ind in range(iCount, 0, -1)
#:set jCounter = INDICES[ind-1]
#:set ws = "  " * (iCount-ind)
      ${ws}$do ${jCounter}$ = 1, size(A, dim=${ind}$)
#:endfor
      ${ws}$  if (mask(${loopindices}$)) then
      ${ws}$    indx(:) = skipDim(dim,[${loopindices}$])
      ${ws}$    fr(${loopindices2}$) = fr(${loopindices2}$) ${op}$ A(${loopindices}$)
      ${ws}$  endif
#:for ind in range(iCount, 0, -1)
#:set ws = "  "*(ind-1)
      ${ws}$enddo
#:endfor

    else

#:for ind in range(iCount, 0, -1)
#:set jCounter = INDICES[ind-1]
#:set ws = "  " * (iCount-ind)
      ${ws}$do ${jCounter}$ = 1, size(A, dim=${ind}$)
#:endfor
      ${ws}$  indx(:) = skipDim(dim,[${loopindices}$])
      ${ws}$  fr(${loopindices2}$) = fr(${loopindices2}$) ${op}$ A(${loopindices}$)
#:for ind in range(iCount, 0, -1)
#:set ws = "  "*(ind-1)
      ${ws}$enddo
#:endfor

    end if

  end function ${subname}$${dim}$D_dim_${label}$${bits}$

#:endif

#:endfor
#:endfor
#:endfor

#:endfor


#:for name, type, label in FLAVOURS
#:for bits in REALMODEL

#! Integers with dual numbers
#:for intbits in INTMODEL

  !> Multiply a matrix of integers with a matrix of dual numbers
  function matmul_int${intbits}$_with_dual${label}$${bits}$(A, B) result(fr)

    !> Matrix of integers
    integer(int${intbits}$), intent(in) :: A(:,:)

    !> Matrix of dual numbers
    type(dual_${name}$${bits}$), intent(in) :: B(:,:)

    !> Dual number matrix product (= A (x) B)
    type(dual_${name}$${bits}$) :: fr(size(A, dim=1), size(B, dim=2))

    type(dual_${name}$${bits}$) :: Xd(size(A, dim=1), size(A, dim=2))

    Xd = A
    fr = matmul(Xd, A)

  end function matmul_int${intbits}$_with_dual${label}$${bits}$

  !> Multiply a matrix of dual numbers with a matrix of integers
  function matmul_dual${label}$${bits}$_with_int${intbits}$(A, B) result(fr)

    !> Matrix of dual numbers
    type(dual_${name}$${bits}$), intent(in) :: A(:,:)

    !> Matrix of integers
    integer(int${intbits}$), intent(in) :: B(:,:)

    !> Dual number matrix product (= A (x) B)
    type(dual_${name}$${bits}$) :: fr(size(A, dim=1), size(B, dim=2))

    type(dual_${name}$${bits}$) :: Xd(size(A, dim=1), size(A, dim=2))

    Xd = B
    fr = matmul(A, Xd)

  end function matmul_dual${label}$${bits}$_with_int${intbits}$

#:endfor

#! Floating point with dual numbers

  !> Multiply a matrix of ${label}$${bits}$ with a matrix of dual numbers
  function matmul_${label}$${bits}$_with_dual${label}$${bits}$(A, B) result(fr)

    !> Matrix of ${label}$ values
    ${type}$(real${bits}$), intent(in) :: A(:,:)

    !> Matrix of dual numbers
    type(dual_${name}$${bits}$), intent(in) :: B(:,:)

    !> Dual number matrix product (= A (x) B)
    type(dual_${name}$${bits}$) :: fr(size(A, dim=1), size(B, dim=2))

    type(dual_${name}$${bits}$) :: Xd(size(A, dim=1), size(A, dim=2))

    Xd = A
    fr = matmul(Xd, A)

  end function matmul_${label}$${bits}$_with_dual${label}$${bits}$

  !> Multiply a matrix of dual numbers with a matrix of ${label}$${bits}$
  function matmul_dual${label}$${bits}$_with_${label}$${bits}$(A, B) result(fr)

    !> Matrix of dual numbers
    type(dual_${name}$${bits}$), intent(in) :: A(:,:)

    !> Matrix of ${label}$ values
    ${type}$(real${bits}$), intent(in) :: B(:,:)

    !> Dual number matrix product (= A (x) B)
    type(dual_${name}$${bits}$) :: fr(size(A, dim=1), size(B, dim=2))

    type(dual_${name}$${bits}$) :: Xd(size(A, dim=1), size(A, dim=2))

    Xd = B
    fr = matmul(A, Xd)

  end function matmul_dual${label}$${bits}$_with_${label}$${bits}$

#! Duals with duals

  !> Multiply together two matrices of dual variables
  !! Note: simple loop based multiplication
  function matmulDuals_${label}$${bits}$(A, B) result(fr)

    !> Dual number matrix
    type(dual_${name}$${bits}$), intent(in) :: A(:,:)

    !> Dual number matrix
    type(dual_${name}$${bits}$), intent(in) :: B(:,:)

    !> Dual number matrix product (= A (x) B)
    type(dual_${name}$${bits}$) :: fr(size(A, dim=1), size(B, dim=2))

    integer :: ii, jj, kk, m, n, p, q, order

    m = size(A,1)
    n = size(A,2)
    p = size(B,1)
    q = size(B,2)

    if (n /= p) then
      print*, "Error: Matrix dimensions do not align for multiplication"
      error stop
    end if

    order = min(size(A(1,1)%f), size(B(1,1)%f)) - 1
    call initialize_dual(fr, order)

    do ii = 1, m
      do jj = 1, q
        do kk = 1, n
          fr(ii,jj) = fr(ii,jj) + A(ii,kk) * B(kk,jj)
        end do
      end do
    end do

  end function matmulDuals_${label}$${bits}$

#:endfor
#:endfor

  !> Dimension along which to sum
  !integer, intent(in), optional :: dim

  function skipDim(ii, indx)

    integer, intent(in) :: ii

    integer, intent(in) :: indx(:)

    integer, allocatable :: skipDim(:)

    integer :: iDim, jj, n

    n = size(indx)
    allocate(skipDim(n-1))
    iDim = 1
    do jj = 1, n
      if (ii == jj) cycle
      skipDim(iDim) = indx(jj)
      iDim = iDim + 1
    end do

  end function skipDim

end module src_matrixfunctions
