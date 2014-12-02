/*
 * optimpack-linalg.h --
 *
 * Definitions for linear algebra routines in OptimPack library.
 *
 *-----------------------------------------------------------------------------
 *
 * The OptimPack library is licensed under the MIT "Expat" License:
 *
 * Copyright (c) 2003-2014: Éric Thiébaut
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 *
 *-----------------------------------------------------------------------------
 */

#ifndef _OPTIMPACK_LINALG_H
#define _OPTIMPACK_LINALG_H 1

#ifndef _OPTIMPACK_H
# include "optimpack.h"
#endif

OPK_BEGIN_C_DECLS

/*---------------------------------------------------------------------------*/
/* Group: BLAS Constants */
/* SBLAS - Simple Basic Linear Algebra Subroutines */

/* Following BLAS/LAPACK conventions, the prefixes 's' and 'd' are
 * used to distinguish between the single/double precision floating
 * point versions of some numerical functions.
 */

#ifdef HAVE_CBLAS

# include <cblas.h>
# define OPK_BLAS_ROW_MAJOR   CblasRowMajor
# define OPK_BLAS_COL_MAJOR   CblasColMajor
# define OPK_BLAS_NO_TRANS    CblasNoTrans
# define OPK_BLAS_TRANS       CblasTrans
# define OPK_BLAS_CONJ_TRANS  CblasConjTrans
# define OPK_BLAS_UPPER       CblasUpper
# define OPK_BLAS_LOWER       CblasLower
# define OPK_BLAS_NON_UNIT    CblasNonUnit
# define OPK_BLAS_UNIT        CblasUnit
# define OPK_BLAS_LEFT        CblasLeft
# define OPK_BLAS_RIGHT       CblasRight
# define opk_blas_order_t     CBLAS_ORDER
# define opk_blas_trans_t     CBLAS_TRANSPOSE
# define opk_blas_uplo_t      CBLAS_UPLO
# define opk_blas_diag_t      CBLAS_DIAG
# define opk_blas_side_t      CBLAS_SIDE

#else /* not HAVE_CBLAS */

/*
 * Enum: opk_blas_order_t
 *
 * OPK_BLAS_ROW_MAJOR - Matrix elements stored in row major order.
 * OPK_BLAS_COL_MAJOR - Matrix elements stored in column major order.
 */
typedef enum _opk_blas_order {
  OPK_BLAS_ROW_MAJOR  = 101,
  OPK_BLAS_COL_MAJOR  = 102
} opk_blas_order_t;

/*
 * Enum: opk_blas_trans_t
 *
 * OPK_BLAS_NO_TRANS   - Do not transpose matrix.
 * OPK_BLAS_TRANS      - Transpose matrix.
 * OPK_BLAS_CONJ_TRANS - Conjugate transpose.
 */
typedef enum _opk_blas_trans {
  OPK_BLAS_NO_TRANS   = 111,
  OPK_BLAS_TRANS      = 112,
  OPK_BLAS_CONJ_TRANS = 113
} opk_blas_trans_t;

/*
 * Enum: opk_blas_uplo_t
 *
 * OPK_BLAS_UPPER - Use upper triangular part.
 * OPK_BLAS_LOWER - Use lower triangular part.
 */
typedef enum _opk_blas_uplo {
  OPK_BLAS_UPPER = 121,
  OPK_BLAS_LOWER = 122
} opk_blas_uplo_t;

/*
 * Enum: opk_blas_diag_t
 *
 *   OPK_BLAS_NON_UNIT - Non unit diagonal.
 *   OPK_BLAS_UNIT     - Unit diagonal.
 */
typedef enum _opk_blas_diag {
  OPK_BLAS_NON_UNIT = 131,
  OPK_BLAS_UNIT     = 132
} opk_blas_diag_t;

/*
 * Enum: opk_blas_side_t
 *
 *   OPK_BLAS_LEFT  - (CblasLeft).
 *   OPK_BLAS_RIGHT - (CblasRight).
 */
typedef enum _opk_blas_side {
  OPK_BLAS_LEFT  = 141,
  OPK_BLAS_RIGHT = 142
} opk_blas_side_t;

#endif /* HAVE_CBLAS */

/*---------------------------------------------------------------------------*/
/* Group: Vector operations */
/* Level 1 BLAS-like routines */

/*
 * Function: opk_damax
 *   Maximum absolute value of a double precision vector.
 *
 * Description:
 *   Get maximum absolute value of a vector.  This is also the infinite norm of
 *   the vector.
 *
 * Parameters:
 *   n    - The number of elements to consider in _x_.
 *   x    - A vector of _n_ * _incx_ values.
 *   incx - The index increment (greater or equal 1).
 *
 * Returns:
 *   The maximum absolute value among the _n_ elements of _x_.
 *
 * See Also:
 *   <opk_samax>, <opk_idamax>.
 */
extern double
opk_damax(opk_index_t n, const double x[],
          opk_index_t incx);

/*
 * Function: opk_samax
 *   Maximum absolute value of a single precision vector.
 *
 * See Also:
 *   <opk_damax>, <opk_isamax>.
 */
extern float
opk_samax(opk_index_t n, const float x[],
          opk_index_t incx);

/*
 * Function: opk_dasum
 *   Sum of the absolute values of a vector.
 *
 * Description:
 *   This function computes the sum of the absolute values of a vector.  This
 *   is also the L-1 norm of the vector.
 *
 * Parameters:
 *   n    - The number of elements to consider in _x_.
 *   x    - A vector of _n_ * _incx_ values.
 *   incx - The index increment (greater or equal 1).
 *
 * Returns:
 *   Returns the sum of the absolute values of _x_.
 *
 * See Also:
 *   <opk_sasum>.
 */
extern double
opk_dasum(opk_index_t n,
          const double x[], opk_index_t incx);

/*
 * Function: opk_sasum
 *   Sum of the absolute values of a vector.
 *
 * See Also:
 *   <opk_dasum>.
 */
extern float
opk_sasum(opk_index_t n,
          const float x[], opk_index_t incx);

/*
 * Function: opk_daxpy
 *   Linear combination of two vectors.
 *
 * Description:
 *   This function stores _a_ * _x_[ix] + _y_[iy] into _y_[iy]. The index
 *   increments may be negative to consider elements in reverse order.  The
 *   code is optimized for _a_ = +/-1 and _a_ = 0.
 *
 * Parameters:
 *   n    - The number of elements to consider in _x_ and _y_.
 *   a    - The scalar factor.
 *   x    - A vector of _n_ * | _incx_ | values.
 *   incx - The index increment for _x_.
 *   y    - A vector of _n_ * | _incy_ | values.
 *   incy - The index increment for _y_.
 *
 * See Also:
 *   <opk_saxpy>.
 */
extern void
opk_daxpy(opk_index_t n, double a,
          const double x[], opk_index_t incx,
          double y[], opk_index_t incy);

/*
 * Function: opk_saxpy
 *   Linear combination of two vectors.
 *
 * See Also:
 *   <opk_daxpy>.
 */
extern void
opk_saxpy(opk_index_t n, float a,
          const float x[], opk_index_t incx,
          float y[], opk_index_t incy);

/*
 * Function: opk_dcopy
 *   Copy a vector into another one.
 *
 * Description:
 *   This function copies _n_ elements of vectors _x_ to vector _y_.  The index
 *   increments may be negative to consider elements in reverse order.
 *
 * Parameters:
 *   n    - The number of elements to consider in _x_ and _y_.
 *   x    - A vector of _n_ * | _incx_ | values.
 *   incx - The index increment for _x_.
 *   y    - A vector of _n_ * | _incy_ | values.
 *   incy - The index increment for _y_.
 *
 * See Also:
 *   <opk_scopy>, <opk_dswap>.
 */
extern void
opk_dcopy(opk_index_t n,
          const double x[], opk_index_t incx,
          double y[], opk_index_t incy);

/*
 * Function: opk_scopy
 *   Copy a vector into another one.
 *
 * See Also:
 *   <opk_dcopy>, <opk_sswap>.
 */
extern void
opk_scopy(opk_index_t n,
          const float x[], opk_index_t incx,
          float y[], opk_index_t incy);

/*
 * Function: opk_ddot
 *   Dot product of two vectors.
 *
 * Description:
 *   This function computes the dot product of two vectors. The index
 *   increments may be negative to consider elements in reverse order.
 *
 * Parameters:
 *   n    - The number of elements to consider in _x_ and _y_.
 *   x    - A vector of _n_ * | _incx_ | values.
 *   incx - The index increment for _x_.
 *   y    - A vector of _n_ * | _incy_ | values.
 *   incy - The index increment for _y_.
 *
 * Returns:
 *   The dot product of _x_ by _y_.
 *
 * See Also:
 *   <opk_sdot>.
 */
extern double
opk_ddot(opk_index_t n,
         const double x[], opk_index_t incx,
         const double y[], opk_index_t incy);

/*
 * Function: opk_sdot
 *   Dot product of two vectors.
 *
 * See Also:
 *   <opk_ddot>.
 */
extern float
opk_sdot(opk_index_t n,
         const float x[], opk_index_t incx,
         const float y[], opk_index_t incy);

/*
 * Function: opk_dnrm2
 *   Euclidean norm of a vector.
 *
 * Description:
 *   This function computes the Euclidean norm of a vector, avoiding
 *   overflows.  This is also the L-2 norm of the vector.
 *
 * Parameters:
 *   n    - The number of elements to consider in _x_.
 *   x    - A vector of _n_ * _incx_ values.
 *   incx - The index increment (greater or equal 1).
 *
 * Returns:
 *   The Euclidean norm of vector _x_: sqrt(x[0]*x[0] + x[incx]*x[incx] + ...).
 *
 * See Also:
 *   <opk_snrm2>.
 */
extern double
opk_dnrm2(opk_index_t n,
          const double x[], opk_index_t incx);

/*
 * Function: opk_snrm2
 *   Euclidean norm of a vector.
 *
 * See Also:
 *   <opk_dnrm2>.
 */
extern float
opk_snrm2(opk_index_t n,
          const float x[], opk_index_t incx);

/*
 * Function: opk_dscal
 *   Scaling of a vector.
 *
 * Description:
 *   This function scales a vector by a scalar, the operation is done in
 *   place.  The _n_ elements of _x_ get multiplied by _a_.  Does nothing if
 *   _n_ or _incx_ are less than 1.  The code is optimized for scalar _a_ =
 *   +/-1 and _a_ = 0.
 *
 * Parameters:
 *   n    - The number of elements to consider in _x_.
 *   a    - The scalar factor.
 *   x    - A vector of _n_ * _incx_ values.
 *   incx - The index increment (greater or equal 1).
 *
 * See Also:
 *   <opk_sscal>.
 */
extern void
opk_dscal(opk_index_t n, double a,
          double x[], opk_index_t incx);

/*
 * Function: opk_sscal
 *   Scaling of a vector.
 *
 * See Also:
 *   <opk_dscal>.
 */
extern void
opk_sscal(opk_index_t n, float a,
          float  x[], opk_index_t incx);

/*
 * Function:  opk_dsum
 *   Sum of the values of a vector.
 *
 * Description:
 *   This function computes the sum of the values of a vector.
 *
 * Parameters:
 *   n    - The number of elements to consider in _x_.
 *   x    - A vector of _n_ * _incx_ values.
 *   incx - The index increment (greater or equal 1).
 *
 * Returns:
 *   Returns the sum of the values of _x_.
 *
 * See Also:
 *   <opk_ssum>, <opk_dasum>.
 */
extern double
 opk_dsum(opk_index_t n,
          const double x[], opk_index_t incx);

/*
 * Function:  opk_ssum
 *   Sum of the absolute values of a vector.
 *
 * See Also:
 *   <opk_dsum>, <opk_sasum>.
 */
extern float
opk_ssum(opk_index_t n,
         const float x[], opk_index_t incx);

/*
 * Function: opk_dswap
 *   Exchanging contents of two vectors.
 *
 * Description:
 *   This function exchanges _n_ elements of vectors _x_ and _y_.  The index
 *   increments may be negative to consider elements in reverse order.
 *
 * Parameters:
 *   n    - The number of elements to consider in _x_ and _y_.
 *   x    - A vector of _n_ * | _incx_ | values.
 *   incx - The index increment for _x_.
 *   y    - A vector of _n_ * | _incy_ | values.
 *   incy - The index increment for _y_.
 *
 * See Also:
 *   <opk_sswap>, <opk_dcopy>.
 */
extern void
opk_dswap(opk_index_t n,
          double x[], opk_index_t incx,
          double y[], opk_index_t incy);

/*
 * Function: opk_sswap
 *   Exchanging contents of two vectors.
 *
 * See Also:
 *   <opk_dswap>, <opk_scopy>.
 */
extern void
opk_sswap(opk_index_t n,
          float x[], opk_index_t incx,
          float y[], opk_index_t incy);

/*
 * Function: opk_dzero
 *   Fill an array with zeros.
 *
 * Description:
 *   This function fills the _n_ elements of vector _x_ by step of _incx_ with
 *   zeros.
 *
 * Parameters:
 *   n    - The number of elements to consider in _x_.
 *   x    - A vector of _n_ * _incx_ values.
 *   incx - The index increment (greater or equal 1).
 *
 * See Also:
 *   <opk_szero>.
 */

extern void
opk_dzero(opk_index_t n,
          double x[], opk_index_t incx);

/*
 * Function: opk_szero
 *   Fill an array with zeros.
 *
 * Description:
 *   This function fills the _n_ elements of vector _x_ by step of _incx_ with
 *   zeros.
 *
 * Parameters:
 *   n    - The number of elements to consider in _x_.
 *   x    - A vector of _n_ * _incx_ values.
 *   incx - The index increment (greater or equal 1).
 *
 * See Also:
 *   <opk_dzero>.
 */
extern void
opk_szero(opk_index_t n,
          float  x[], opk_index_t incx);

/*
 * Function: opk_idamax
 *   Get index of maximum absolute value of a vector.
 *
 * Description:
 *   This function returns the index of the maximum absolute value of a vector.
 *   Following FORTRAN conventions, indices are 1-based.  If _n_ <= 0 or
 *   _incx_ <= 0, the returned value is 0.
 *
 * Parameters:
 *   n    - The number of elements to consider in _x_.
 *   x    - A vector of _n_ * _incx_ values.
 *   incx - The index increment (greater or equal 1).
 *
 * Returns:
 *   Returns 1-based index of maximum absolute value of _x_.
 *   Returns 0 if _n_ <= 0 or _incx_ <= 0.
 *
 * See Also:
 *   <opk_isamax>, <opk_damax>.
 */
extern opk_index_t
opk_idamax(opk_index_t n,
           const double x[], opk_index_t incx);

/*
 * Function: opk_isamax
 *   Get index of maximum absolute value of a vector.
 *
 * See Also:
 *   <opk_idamax>, <opk_samax>.
 */
extern opk_index_t
opk_isamax(opk_index_t n,
           const float x[], opk_index_t incx);

/*---------------------------------------------------------------------------*/
/* Group: Matrix-vector operations */
/* Level 2 BLAS-like routines */

/*
 * Function: opk_dgemv
 *   Matrix-vector operation.
 *
 * Description:
 *   This function performs one of the matrix-vector operations
 *
 *   > y := alpha A.x + beta y,
 *
 *   or
 *
 *   >  y := alpha A'.x + beta y,
 *
 *  where _alpha_ and _beta_ are scalars, _x_ and _y_ are vectors and _A_ is
 *  an _m_ by _n_ matrix.
 *
 * Parameters:
 *   trans - Transpose flag, specifies the matrix multiplication to be
 *           performed as follows: OPK_BLAS_NO_TRANS for A.x, or
 *           OPK_BLAS_TRANS for A'.x.
 *   m     - The number of rows of matrix _a_.
 *   n     - The number of columns of matrix _a_.
 *   alpha - The first scalar parameter.
 *   a     - Matrix of dimensions (_dla_, _n_).
 *   lda   - The leading dimension of _a_; _lda_ >= max(1,_m_).
 *   x     - The source vector of length at least 1+(_n_-1)*abs(_incx_)
 *           when _trans_ =  OPK_BLAS_NO_TRANS or 1+(_m_-1)*abs(_incx_)
 *           otherwise.
 *   incx  - The index increment for vector _x_ (must be non-zero).
 *   beta  - The second scalar parameter.
 *   y     - The destination vector of length at least 1+(_m_-1)*abs(_incy_)
 *           when _trans_ =  OPK_BLAS_NO_TRANS or 1+(_n_-1)*abs(_incy_)
 *           otherwise.
 *   incy  - The index increment for vector _y_ (must be non-zero).
 *
 * Results:
 *   Non-zero result _k_ means invalid _k_-th argument.
 *
 * See Also:
 *   <opk_sgemv>.
 */
extern int
opk_dgemv(opk_blas_trans_t trans,
          opk_index_t m, opk_index_t n,
          double alpha,
          const double a[], opk_index_t lda,
          const double x[], opk_index_t incx,
          double beta,
          double y[], opk_index_t incy);

/*
 * Function: opk_sgemv
 *   Matrix-vector operation.
 *
 * See Also:
 *   <opk_dgemv>.
 */
extern int
opk_sgemv(opk_blas_trans_t trans,
          opk_index_t m, opk_index_t n,
          float alpha,
          const float a[], opk_index_t lda,
          const float x[], opk_index_t incx,
          float beta,
          float y[], opk_index_t incy);

/*
 * Function: opk_dtrmv
 *   Multiplication of a vector by a triangular matrix.
 *
 * Description:
 *   This function performs one of the matrix-vector operations:
 *   > x := A.x,
 *   or
 *   > x := A'.x,
 *   where _x_ is an _n_ element vector and _A_ is an _n_ by _n_ unit, or
 *   non-unit, upper or lower triangular matrix.
 *
 * Parameters:
 *   uplo  - Specifies whether the matrix is an upper or lower triangular
 *           matrix as follows: _uplo_ = <OPK_BLAS_UPPER> if _A_ is an upper
 *           triangular matrix, and _uplo_ = <OPK_BLAS_LOWER> if _A_ is a lower
 *           triangular matrix.
 *   trans - Specifies the operation to be performed as follows: if _trans_ =
 *           <OPK_BLAS_NO_TRANS>, then x := A.x; otherwise if _trans_ =
 *           <OPK_BLAS_TRANS> or <OPK_BLAS_CONJ_TRANS>, then x := A'.x.
 *   diag  - Specifies whether or not _A_ is unit triangular as follows: if
 *           _diag_ = <OPK_BLAS_UNIT>, _A_ is assumed to be unit triangular;
 *           else if _diag_ = <OPK_BLAS_NON_UNIT>, _A_ is not assumed to be unit
 *           triangular.
 *   n     - The order of the matrix _A_, must be at least zero.
 *   a     - An array of dimension _lda_ by _n_.
 *           With _uplo_ = <OPK_BLAS_UPPER>, the leading _n_ by _n_
 *           upper triangular part of the array _A_ must contain the upper
 *           triangular matrix and the strictly lower triangular part of
 *           _A_ is not used.
 *           With _uplo_ = <OPK_BLAS_LOWER>, the leading _n_ by _n_
 *           lower triangular part of the array _A_ must contain the lower
 *           triangular matrix and the strictly upper triangular part of
 *           _A_ is not used.
 *           Note that when _diag_ = <OPK_BLAS_UNIT>, the diagonal elements of
 *           _A_ are not referenced either, but are assumed to be unity.
 *   lda   - The leading dimension of _A_, must be at least max(_n_,1).
 *   x     - An array of at least (1+(_n_-1)*abs(_incx_) elements.  Before
 *           entry, the incremented array must contain the _n_ element vector
 *           _x_. On exit, it is overwritten with the tranformed vector _x_.
 *  incx   - The increment for the elements of _x_, must not be zero.
 *
 * Results:
 *   Non-zero result _k_ means invalid _k_-th argument.
 *
 * See Also:
 *   <opk_strmv>, <opk_dtrsv>.
 */
extern int
opk_dtrmv(opk_blas_uplo_t uplo,
          opk_blas_trans_t trans,
          opk_blas_diag_t diag,
          opk_index_t n,
          const double a[], opk_index_t lda,
          double x[], opk_index_t incx);

/*
 * Function: opk_strmv
 *   Multiplication of a vector by a triangular matrix.
 *
 * See Also:
 *   <opk_dtrmv>, <opk_strsv>.
 */
extern int
opk_strmv(opk_blas_uplo_t uplo,
          opk_blas_trans_t trans,
          opk_blas_diag_t diag,
          opk_index_t n,
          const float a[], opk_index_t lda,
          float x[], opk_index_t incx);

/*
 * Function: opk_dtrsv
 *   Solves a tringular linear system of equations.
 *
 * Description:
 *   This function solves one of the matrix-vector operations:
 *   > A.x = b,
 *   or
 *   > A'.x = b,
 *   where _x_ and _b_ are _n_ element vectors and _A_ is an _n_ by _n_ unit,
 *   or non-unit, upper or lower triangular matrix.
 *
 *   No test for singularity or near-singularity is included in this
 *   routine. Such tests must be performed before calling this routine.
 *
 * Parameters:
 *   uplo  - Specifies whether the matrix is an upper or lower triangular
 *           matrix as follows: _uplo_ = OPK_BLAS_UPPER if _A_ is an upper
 *           triangular matrix, and _uplo_ = OPK_BLAS_LOWER if _A_ is a lower
 *           triangular matrix.
 *   trans - Specifies the operation to be performed as follows: if _trans_ =
 *           OPK_BLAS_NO_TRANS, then solves A.x = b; otherwise if _trans_ =
 *           OPK_BLAS_TRANS or OPK_BLAS_CONJ_TRANS, then solves A'.x = b.
 *   diag  - Specifies whether or not _A_ is unit triangular as follows: if
 *           _diag_ = OPK_BLAS_UNIT, _A_ is assumed to be unit triangular;
 *           else if _diag_ = OPK_BLAS_NON_UNIT, _A_ is not assumed to be unit
 *           triangular.
 *   n     - The order of the matrix _A_, must be at least zero.
 *   a     - An array of dimension _lda_ by _n_.
 *           With _uplo_ = OPK_BLAS_UPPER, the leading _n_ by _n_
 *           upper triangular part of the array _A_ must contain the upper
 *           triangular matrix and the strictly lower triangular part of
 *           _A_ is not used.
 *           With _uplo_ = OPK_BLAS_LOWER, the leading _n_ by _n_
 *           lower triangular part of the array _A_ must contain the lower
 *           triangular matrix and the strictly upper triangular part of
 *           _A_ is not used.
 *           Note that when _diag_ = OPK_BLAS_UNIT, the diagonal elements of
 *           _A_ are not referenced either, but are assumed to be unity.
 *   lda   - The leading dimension of _A_, must be at least max(_n_,1).
 *   x -     An array of at least (1+(_n_-1)*abs(_incx_) elements.  Before
 *           entry, the incremented array _x_ must contain the _n_ element
 *           vector _b_. On exit, _x_ is overwritten with the solution vector
 *           _x_.
 *  incx   - The increment for the elements of _x_, must not be zero.
 *
 * Results:
 *   Non-zero result _k_ means invalid _k_-th argument.
 *
 * See Also:
 *   <opk_strsv>, <opk_dtrmv>.
 */
extern int
opk_dtrsv(opk_blas_uplo_t uplo,
          opk_blas_trans_t trans,
          opk_blas_diag_t diag,
          opk_index_t n,
          const double a[], opk_index_t lda,
          double x[], opk_index_t incx);

/*
 * Function: opk_strsv
 *   Solves a tringular linear system of equations.
 *
 * See Also:
 *   <opk_dtrsv>, <opk_strmv>.
 */
extern int
opk_strsv(opk_blas_uplo_t uplo,
          opk_blas_trans_t trans,
          opk_blas_diag_t diag,
          opk_index_t n,
          const float a[], opk_index_t lda,
          float x[], opk_index_t incx);

/*---------------------------------------------------------------------------*/
/* Group: Matrix-matrix operations */
/* Level 3 BLAS-like routines */

/*
 * Function: opk_dgemm
 *   Performs a matrix-matrix operation.
 *
 * Description:
 *   This function performs one of the matrix-matrix operations
 *   > C := alpha op(A).op(B) + beta C,
 *   where op(X) = X or op(X) = X', _alpha_ and _beta_ are scalars, and _A_,
 *   _B_ and _C_ are matrices, with op(_A_) an _m_ by _k_ matrix, op(_B_) a
 *   _k_ by _n_ matrix and C an _m_ by _n_ matrix.
 *
 * Parameters:
 *   transa - Specifies the form of op(_A_) to be used in the matrix
 *            multiplication as follows: if _transa_ = OPK_BLAS_NO_TRANS, then
 *            op(_A_) = _A_; otherwise if _trans_ = OPK_BLAS_TRANS or
 *            OPK_BLAS_CONJ_TRANS, then op(_A_) = _A_'.
 *   transb - Specifies the form of op(_B_) to be used in the matrix
 *            multiplication; see _transa_.
 *   m      - The number of rows of the matrix op(_A_) and of the matrix _C_,
 *            must be at least zero.
 *   n      - The number of columns of the matrix op(_B_) and of the matrix
 *            _C_, must be at least zero.
 *   k      - The number of columns of the matrix op(_A_) and the number of
 *            rows of the matrix op(_B_), must be at least zero.
 *   alpha  - The scalar _alpha_.
 *   a      - An array of dimension _lda_ by _ka_, where _ka_ is _k_ when
 *            _transa_ = OPK_BLAS_NO_TRANS, and is _m_ otherwise.  Before entry
 *            with _transa_ = OPK_BLAS_NO_TRANS, the leading _m_ by _k_ part of
 *            the array _a_ must contain the matrix _A_, otherwise the leading
 *            _k_ by _m_ part of the array _a_ must contain the matrix _A_.
 *   lda    - The first dimension of _a_. When _transa_ = OPK_BLAS_NO_TRANS
 *            then _lda_ must be at least max(_m_,1), otherwise _lda_ must be
 *            at least max(_k_,1).
 *   b      - An array of dimension _ldb_ by _kb_, where _kb_ is _n_ when
 *            _transb_ = OPK_BLAS_NO_TRANS, and is _k_ otherwise.  Before
 *            entry with _transb_ = OPK_BLAS_NO_TRANS, the leading _k_ by _n_
 *            part of the array _b_ must contain the matrix _B_, otherwise the
 *            leading _n_ by _k_ part of the array _b_ must contain the matrix
 *            _B_.
 *   ldb    - The first dimension of _b_. When _transb_ = OPK_BLAS_NO_TRANS
 *            then _ldb_ must be at least max(_k_,1), otherwise _ldb_ must be
 *            at least max(_n_,1).
 *   beta   - The scalar _beta_.  When _beta_ is supplied as zero then _C_ need
 *            not be set on input.
 *   c      - An array of dimension _ldc_ by _n_.
 *            Before entry, the leading _m_ by _n_ part of the array _c_ must
 *            contain the matrix _C_, except when _beta_ is zero, in which
 *            case _c_ need not be set on entry.
 *            On exit, the array _c_ is overwritten by the _m_ by _n_ matrix
 *            _alpha_ op(_A_)*op(_B_) + _beta_ _C_.
 *   ldc    - The first dimension of _c_. _ldc_  must  be  at  least max(_m_,1).
 *
 * Results:
 *   Non-zero result _k_ means invalid _k_-th argument.
 *
 * See Also:
 *   <opk_sgemm>.
 */
extern int
opk_dgemm(opk_blas_trans_t transa,
          opk_blas_trans_t transb,
          opk_index_t m,
          opk_index_t n,
          opk_index_t k,
          double alpha,
          const double a[], opk_index_t lda,
          const double b[], opk_index_t ldb,
          double beta,
          double c[], opk_index_t ldc);

/*
 * Function: opk_sgemm
 *   Performs a matrix-matrix operation.
 *
 * See Also:
 *   <opk_dgemm>.
 */
extern int
opk_sgemm(opk_blas_trans_t transa,
          opk_blas_trans_t transb,
          opk_index_t m,
          opk_index_t n,
          opk_index_t k,
          float alpha,
          const float a[], opk_index_t lda,
          const float b[], opk_index_t ldb,
          float beta,
          float c[], opk_index_t ldc);

/*
 * Function: opk_dsyrk
 *   Performs a symmetric rank k operation.
 *
 * Description:
 *   This function performs one of the symmetric rank k operations
 *   > C := alpha*A*A' + beta*C,
 *   or
 *   > C := alpha*A'*A + beta*C,
 *   where _alpha_ and _beta_ are scalars, _C_ is an _n_ by _n_ symmetric
 *   matrix and _A_ is an _n_ by _k_ matrix in the first case and a _k_ by _n_
 *   matrix in the second case.
 *
 * Parameters:
 *   uplo  - Specifies whether the matrix is the upper or lower triangular part
 *           of the array _C_ is to be referenced as follows.  If _uplo_ =
 *           <OPK_BLAS_UPPER>, only the upper triangular part of _C_ is to be
 *           referenced.  If _uplo_ = <OPK_BLAS_LOWER>, only the lower
 *           triangular part of _C_ is to be referenced.
 *   trans - Specifies the operation to be performed as follows: C :=
 *           alpha*A*A' + beta*C, if _trans_ = <OPK_BLAS_NO_TRANS>; or C :=
 *           alpha*A'*A + beta*C, _trans_ = <OPK_BLAS_TRANS> or
 *           <OPK_BLAS_CONJ_TRANS>.
 *   n     - The order of the matrix _C_; _N_ must be at least zero.
 *   k     - On entry with _trans_ = <OPK_BLAS_NO_TRANS>, _k_ specifies the
 *           number of columns of the matrix _A_, and on entry with _trans_ =
 *           <OPK_BLAS_TRANS> or <OPK_BLAS_CONJ_TRANS>, _k_ specifies the
 *           number of rows of the matrix _A_.  _K_ must be at least zero.
 *   alpha - The scalar _alpha_.
 *   a     - Real array of dimension _lda_ by _ka_, where _ka_ is _k_ when
 *           _trans_ = <OPK_BLAS_NO_TRANS>, and is _n_ otherwise.  Before
 *           entry with _trans_ = <OPK_BLAS_NO_TRANS>, the leading _n_ by _k_
 *           part of the array _a_ must contain the matrix _A_, otherwise the
 *           leading _k_ by _n_ part of the array _a_ must contain the matrix
 *           _A_.
 *   lda   - The leading dimension of _a_.  When _trans_ = <OPK_BLAS_NO_TRANS>,
 *           then _lda_ must be at least max(_n_,1), otherwise _LDA_ must be
 *           at least max(_k_,1).
 *   beta  - The scalar _beta_.
 *   c     - Real array of dimension _ldc_ by _n_.
 *           Before entry with _uplo_ = <OPK_BLAS_UPPER>, the leading _n_ by
 *           _n_ upper triangular part of the array _c_ must contain the upper
 *           triangular part of the symmetric matrix and the strictly lower
 *           triangular part of _c_ is not referenced.  On exit, the upper
 *           triangular part of the array _c_ is overwritten by the upper
 *           triangular part of the updated matrix.
 *           Before entry with _uplo_ = <OPK_BLAS_LOWER>, the leading _n_ by
 *           _n_ lower triangular part of the array _c_ must contain the lower
 *           triangular part of the symmetric matrix and the strictly upper
 *           triangular part of _c_ is not referenced.  On exit, the lower
 *           triangular part of the array _c_ is overwritten by the lower
 *           triangular part of the updated matrix.
 *   ldc   - The leading dimension of _c_.  _ldc_  must  be  at  least
 *           max(_n_,1).
 *
 * Results:
 *   Non-zero result _k_ means invalid _k_-th argument.
 *
 * See Also:
 *   <opk_ssyrk>.
 */
extern int
opk_dsyrk(opk_blas_uplo_t uplo,
          opk_blas_trans_t trans,
          opk_index_t n,
          opk_index_t k,
          double alpha,
          const double a[], opk_index_t lda,
          double beta,
          double c[], opk_index_t ldc);

/*
 * Function: opk_ssyrk
 *   Performs a symmetric rank k operation.
 *
 * See Also:
 *   <opk_dsyrk>.
 */
extern int
opk_ssyrk(opk_blas_uplo_t uplo,
          opk_blas_trans_t trans,
          opk_index_t n,
          opk_index_t k,
          float alpha,
          const float a[], opk_index_t lda,
          float beta,
          float c[], opk_index_t ldc);

/*---------------------------------------------------------------------------*/
/* Group: LAPACK-like routines */

/*
 * Function: opk_dpotf2
 *   Cholesky factorization of a real symmetric positive definite matrix.
 *
 * Description:
 *   This function computes the Cholesky factorization of a real symmetric
 *   positive definite matrix _A_.
 *
 *   The factorization has the form _A_ = _U_'._U_, if _uplo_ =
 *   <OPK_BLAS_UPPER>, or _A_ = _L_._L_', if _uplo_ = <OPK_BLAS_LOWER> where
 *   _U_ is an upper triangular matrix and _L_ is lower triangular.
 *
 *  This is the unblocked version of the algorithm, calling Level 2 BLAS.
 *
 * Parameters:
 *
 *   uplo - Specifies whether the upper or lower triangular part of the
 *          symmetric matrix _A_ is stored. _uplo_ is one of: <OPK_BLAS_UPPER>
 *          or <OPK_BLAS_LOWER>.
 *   n    - The order of the matrix _A_, must be greater of equal zero.
 *   a    - Input/output array of dimension _lda_ by _n_.
 *          On entry, _a_ contains the symmetric matrix _A_.  If _uplo_ =
 *          <OPK_BLAS_UPPER>, the leading _n_ by _n_ upper triangular part of
 *          _a_ contains the upper triangular part of the matrix _A_, and the
 *          strictly lower triangular part of _a_ is not referenced.  If
 *          _uplo_ = <OPK_BLAS_LOWER>, the leading _n_ by _n_ lower triangular
 *          part of _a_ contains the lower triangular part of the matrix _A_,
 *          and the strictly upper triangular part of _a_ is not referenced.
 *          On successful exit (with returned value 0), _a_ contains the
 *          factor _U_ or _L_ from the Cholesky factorization _A_ = _U'_ . _U_
 *          or _A_ = _L_ . _L'_ .
 *
 *   lda   - The leading dimension of the array _a_. _lda_ >= max(_n_,1).
 *
 * Returns:
 *   The returned value, says _info_ is:
 *   * _info_ = 0 - successful exit.
 *   * _info_ < 0 - the _k_-th argument had an illegal value, with
 *                  _k_ = -_info_.
 *   * _info_ > 0 - the leading minor of order _k_ is not positive definite,
 *                  and the factorization could not be completed, with _k_ =
 *                  _info_.
 *
 * See Also:
 *   <opk_spotf2>.
 */
extern opk_index_t
opk_dpotf2(opk_blas_uplo_t uplo,
           opk_index_t n,
           double a[], opk_index_t lda);

/*
 * Function: opk_spotf2
 *   Cholesky factorization of a real symmetric positive definite matrix.
 *
 * See Also:
 *   <opk_dpotf2>.
 */
extern opk_index_t
opk_spotf2(opk_blas_uplo_t uplo,
           opk_index_t n,
           float a[], opk_index_t lda);

/*---------------------------------------------------------------------------*/

/* Group: Linear Conjugate Gradient */

/*
 * This package provides methods to solve a linear system of equations with a
 * positive definite left hand-side matrix by means of iterative conjugate
 * gradient method.
 */

/*
 * Enum: opk_cg_state_t
 *   The possible values of Conjugate Gradient state.
 *
 *   OPK_CG_ERROR      - error;
 *   OPK_CG_START      - start with no initial variable (all zero);
 *   OPK_CG_RESTART    - start or restart with initial variables;
 *   OPK_CG_NEWX       - new step taken;
 *   OPK_CG_PRECOND    - caller has to to compute z = Q.r;
 *   OPK_CG_AP         - caller has to to compute q = A.p;
 *   OPK_CG_FINISH     - convergence or rounding errors prevent further
 *                       progress;
 *   OPK_CG_NON_CONVEX - non-positive definitiveness has been detected;
 *   OPK_CG_TRUNCATED  - step truncated at trust region boundary.
 */
typedef enum {
  OPK_CG_ERROR      = -1,
  OPK_CG_START      =  0,
  OPK_CG_RESTART    =  1,
  OPK_CG_NEWX       =  2,
  OPK_CG_PRECOND    =  3,
  OPK_CG_AP         =  4,
  OPK_CG_FINISH     =  5,
  OPK_CG_NON_CONVEX =  6,
  OPK_CG_TRUNCATED  =  7
} opk_cg_state_t;

/*
 * Function: opk_plcg_d
 *   Preconditioned linear conjugate gradient (double precision).
 *
 * Description:
 *   Iteratively solve a linear system _A.x = b_ where _A_ is a symmetric
 *   positive definite matrix by preconditioned conjugate gradient method with
 *   reverse communication.
 *
 * Parameters:
 *   n - The number of variables.
 *   p - A work array of size _n_; see parameter _state_ for explanations.
 *   q - A work array of size _n_; on return with _state_ set to <OPK_CG_AP>,
 *       the caller must compute the matrix vector multiplication _q=A.p_
 *       (i.e. store into vector _q_ the result of multiplying vector _p_ by
 *       the matrix _A_); on return with _state_ set to <OPK_CG_NEWX>, _q_ is
 *       filled with the unscaled residuals step.
 *   r - An array of size _n_; on initialization or restart (with state set to
 *       <OPK_CG_START> or <OPK_CG_RESTART>), _r_ stores the right hand side
 *       vector _b_; otherwise, _r_ stores the current residuals.
 *   x - An array of size _n_ to store the unknowns.  On initialization with
 *       given parameter values or restart, i.e. with _state_ set to
 *       <OPK_CG_RESTART>, _x_ stores the initial unknowns _x0_; on initial
 *       entry with _state_ set to <OPK_CG_START>, the contents of _x_ is
 *       ignored (_x_ will be filled with zeros); otherwise, _x_ stores the
 *       current solution.
 *   z - Optional work array of size _n_; on return with _state_ set to
 *       <OPK_CG_PRECOND>, the caller must compute and store the
 *       preconditioned residuals into _z_ given the current residuals which
 *       are in _r_.  If no preconditioning is to be used, _z_ must be set to
 *       _NULL_ (in which case, a _state_ with value <OPK_CG_PRECOND> will
 *       never get returned).  If not _NULL_, _z_ can be the same array as _q_
 *       to save memory.
 *   rho - A 4-element array used to store the dot product _r'.z_ of the
 *       current and previous iterations and to store some parameters of the
 *       algorithm.  _rho[0]_ is the former value of _rho[1]_; _rho[1]_ is the
 *       dot product of the residuals _r_ by the (preconditioned) residuals
 *       _z_; _rho[2]_ is the optimal step size (_alpha_); _rho[3]_ is the
 *       weight of the former search direction (_beta_).
 *   state - Address of integer variable to store the current stage of the
 *       algorithm (must not be altered by the caller between calls except to
 *       restart the algorithm).  On initialization or restart of the
 *       conjugate gradient iterations, _state_ must be set to <OPK_CG_START>
 *       or <OPK_CG_RESTART> and _r_ must be set to the right hand side
 *       vector, that is: _r=b_; if _state_ is <OPK_CG_START>, _x_ is filled
 *       with zeros by the PLCG routine; if _state_ is <OPK_CG_RESTART>, _x_
 *       must be set with some initial values of the parameters, that is:
 *       _x=x0_. Upon return of the function, the value stored in _state_ is
 *       one of:
 *
 *      - <OPK_CG_NEWX>: _x_ stores the current solution, _r_ stores the
 *        current residuals, the change in parameters is equal to _rho[2]*p_
 *        and the corresponding residuals change is equal to _-rho[2]*q_.  The
 *        caller can decide to terminate the iterations, to pursue the
 *        iterations or to restart the algorithm.
 *
 *      - <OPK_CG_AP>: The caller must compute the matrix vector
 *        multiplication _q=A.p_ (i.e., store into vector _q_ the result of
 *        multiplying vector _p_ by the matrix _A_) and call PLCG again.
 *
 *      - <OPK_CG_PRECOND>: The caller must compute and store the
 *        preconditioned residuals into _z_ (given the current residuals which
 *        are in _r_) and call PLCG again.
 *
 *      - <OPK_CG_FINISH>: No further progress are possible either because the
 *        exact solution has been found or because rounding errors prevent
 *        further progress.
 *
 *      - <OPK_CG_ERROR>: An error has been detected by the algorithm: invalid
 *        parameter or non-positive definitiveness of the matrix _A_ or of the
 *        preconditioner.
 *
 *     For consistency reasons, <OPK_CG_NEWX> is returned after internal
 *     initialization, i.e. _x_ is just _x0_ and _r_ is just _r = b - A.x0_.
 *
 *
 * Example:
 *   In the following piece of code, we assume that the arrays _a_ and _c_
 *   store respectively the coefficients of the matrix _A_ and of the
 *   preconditioner.
 *
 *   (code)
 *   for (i = 0; i < n; ++i) {
 *     x[i] = ...; // initial set of parameters
 *   }
 *   for (i = 0; i < n; ++i) {
 *     r[i] = b[i];        // copy b into r
 *   }
 *   state = OPK_CG_RESTART;
 *   for (;;) {
 *     opk_plcg_d(n, p, q, r, x, z, rho, &state);
 *     if (state == OPK_CG_AP) {
 *       for (i = 0; i < n; ++i) {
 *         double s = 0.0;
 *         for (j = 0; j < n; ++j) {
 *           s += a[i*n + j]*p[j];
 *         }
 *         q[i] = s;
 *       }
 *     } else if (state == OPK_CG_PRECOND) {
 *       for (i = 0; i < n; ++i) {
 *         double s = 0.0;
 *         for (j = 0; j < n; ++j) {
 *           s += c[i*n + j]*r[j];
 *         }
 *         z[i] = s;
 *       }
 *     } else if (state == OPK_CG_NEWX) {
 *       double s = 0.0;
 *       for (i = 0; i < n; ++i) {
 *         s += r[i]*r[i];
 *       }
 *       if (s < 1e-3) {
 *         break; // convergence
 *       }
 *     } else {
 *       // unexpected (must be an error)
 *       fprintf(stderr, "error in PLCG\n");
 *       break;
 *     }
 *   }
 *   (end code)
 *
 * See Also:
 *   <opk_lcg_d()>, <opk_plcg_f()>.
 */
extern void
opk_plcg_d(opk_index_t n, double p[], double q[],
           double r[], double x[], double z[],
           double rho[4], opk_cg_state_t *state);

/*
 * Function: opk_plcg_f
 *   Preconditioned linear conjugate gradient (single precision).
 *
 * See Also:
 *   <opk_lcg_f()>, <opk_plcg_d()> for the meaning of the arguments.
 */
extern void
opk_plcg_f(opk_index_t n, float p[], float q[],
           float r[], float x[], float z[],
           float rho[4], opk_cg_state_t *state);

/*
 * Function: opk_lcg_d
 *   Linear conjugate gradient (double precision).
 *
 * Description:
 *   This function implements reverse communication linear conjugate gradient
 *   without predconditioner for double precision floating point.  See
 *   <opk_plcg_d()> for the meaning of the arguments and <opk_lcg_f()> for a
 *   single precision version.
 *
 * See Also:
 *   <opk_lcg_f()>, <opk_plcg_d()>.
 */
extern void
opk_lcg_d(opk_index_t n, double p[], double q[], double r[],
          double x[], double rho[4], opk_cg_state_t *state);

/*
 * Function: opk_lcg_f
 *   Linear conjugate gradient (single precision).
 *
 * Description:
 *   This function implements reverse communication linear conjugate gradient
 *   without predconditioner for single precision floating point.  See
 *   <opk_plcg_d()> for the meaning of the arguments and <opk_lcg_d()> for a
 *   double precision version.
 *
 * See Also:
 *   <opk_lcg_d()>, <opk_plcg_d()>.
 */
extern void
opk_lcg_f(opk_index_t n, float p[], float q[], float r[],
          float x[], float rho[4], opk_cg_state_t *state);

/*
 * Function: opk_trcg_d
 *   Trust region conjugate gradient (double precision).
 *
 * Description:
 *   This function implements a reverse communication version of a trust
 *   region linear conjugate gradient algorithm with optional preconditioning.
 *   The trust-region (or truncated) conjugate gradient method is due to
 *   Steihaug (see References below).  This version is for double precision
 *   variables, see <opk_trcg_f()> for a single precision version.  See
 *   <opk_plcg_d()> for an example of using reverse communication.
 *
 * Parameters:
 *   n - The number of variables.
 *   p - A work array of size _n_; see parameter _state_ for explanations.
 *   q - A work array of size _n_; on return with _state_ set to <OPK_CG_AP>,
 *       the caller must compute the matrix vector multiplication _q=A.p_
 *       (i.e. store into vector _q_ the result of multiplying vector _p_ by
 *       the matrix _A_); on return with _state_ set to <OPK_CG_NEWX>, _q_ is
 *       filled with the unscaled residuals step.
 *   r - An array of size _n_; on initialization or restart (i.e. with _state_
 *       set to <OPK_CG_START> or <OPK_CG_RESTART>), _r_ stores the right hand
 *       side vector _b_; otherwise, _r_ stores the current residuals.
 *   x - An array of size _n_ to store the unknowns.  On initialization with
 *       given parameter values or restart, i.e. with _state_ set to
 *       <OPK_CG_RESTART>, _x_ stores the initial unknowns _x0_; on initial
 *       entry with _state_ set to <OPK_CG_START>, the contents of _x_ is
 *       ignored (_x_ will be filled with zeroes); otherwise, _x_ stores the
 *       current solution.
 *   z - Optional work array of size _n_; on return with _state_ set to
 *       <OPK_CG_PRECOND>, the caller must compute and store the
 *       preconditioned residuals into _z_ given the current residuals which
 *       are in _r_.  If no preconditioning is to be used, _z_ must be set to
 *       _NULL_ (in which case, a _state_ with value <OPK_CG_PRECOND> will
 *       never get returned).
 *   delta - The maximum length (Euclidean norm) of vector _x_.  This value
 *       must not be changed once the algorithm is started.
 *   rho - A 5-element array used to store the dot product _r'.z_ of the
 *       current and previous iterations and to store some parameters of the
 *       algorithm.  _rho[0]_ is the former value of _rho[1]_; _rho[1]_ is the
 *       dot product of the residuals _r_ by the (preconditioned) residuals
 *       _z_; _rho[2]_ is the optimal step size (_alpha_); _rho[3]_ is the
 *       weight of the former search direction (_beta_);  _rho[4]_ is the
 *       Euclidean norm of _x_.
 *   state - Address of integer variable to store the current stage of the
 *       algorithm (must not be altered by the caller between calls except to
 *       restart the algorithm).  On initialization or restart of the
 *       conjugate gradient iterations, _state_ must be set to <OPK_CG_START>
 *       or <OPK_CG_RESTART> and _r_ must be set to the right hand side
 *       vector, that is: _r=b_; if _state_ is <OPK_CG_START>, _x_ is filled
 *       with zeroes by the PLCG routine; if _state_ is <OPK_CG_RESTART>, _x_
 *       must be set with some initial values of the parameters, that is:
 *       _x=x0_. Upon return of the function, the value stored in _state_ is
 *       one of:
 *
 *      - <OPK_CG_NEWX>: _x_ stores the current solution, _r_ stores the
 *        current residuals, the change in parameters is equal to _rho[2]*p_
 *        and the corresponding residuals change is equal to _-rho[2]*q_.  The
 *        caller can decide to terminate the iterations, to pursue the
 *        iterations or to restart the algorithm.
 *
 *      - <OPK_CG_PRECOND>: The caller must compute and store the
 *        preconditioned residuals into _z_ (given the current residuals which
 *        are in _r_) and call PLCG again.
 *
 *      - <OPK_CG_AP>: The caller must compute the matrix vector
 *        multiplication _q=A.p_ (i.e., store into vector _q_ the result of
 *        multiplying vector _p_ by the matrix _A_) and call PLCG again.
 *
 *      - <OPK_CG_FINISH>: No further progress are possible either because the
 *        exact solution has been found or because rounding errors prevent
 *        further progress.  The Euclidean norm of _x_ is less than _delta_.
 *
 *      - <OPK_CG_NON_CONVEX>: A further conjugate gradient step cannot be
 *        taken because non positive definitiveness of the preconditioner has
 *        been detected.  The Euclidean norm of _x_ is less than _delta_.
 *
 *      - <OPK_CG_TRUNCATED>: A truncated gradient step has been taken.  The
 *        Euclidean norm of _x_ is equal to _delta_. No further iterations
 *        should be taken (unless with a larger _delta_ and _state_ set to
 *        <OPK_CG_RESTART>).
 *
 *      - <OPK_CG_ERROR>: An error has been detected by the algorithm: invalid
 *        parameter, or corrupted workspace _rho_, or non-positive
 *        definitiveness of the matrix _A_ or of the preconditioner.
 *
 *   For consistency reasons, <OPK_CG_NEWX> is returned after internal
 *   initialization, i.e. _x_ is just _x0_ and _r_ is just _r = b - A.x0_.
 *
 * Restrictions:
 *   It is a bad idea to restart the algorithm with an initial _x_ with norm
 *   greater than _delta_, in this case, _state_ is set to <OPK_CG_TRUNCATED>
 *   and the initial _x_ is rescaled to have a Euclidean norm equals to
 *   _delta_.
 *
 * References:
 *   - J. Nocedal and S. J. Wright, "Numerical Optimization", Springer Verlag,
 *     2006.
 *   - T. Steihaug, "The conjugate gradient method and trust regions in large
 *     scale optimization", SIAM Journal on Numerical Analysis, vol. 20,
 *     pp. 626-637, 1983.
 *
 * See Also:
 *   <opk_trcg_f()>, <opk_plcg_d()>.
 */
extern void
opk_trcg_d(opk_index_t n, double p[], const double q[],
           double r[], double x[], const double z[],
           double delta, double rho[5],
           opk_cg_state_t *state);

/*
 * Function: opk_trcg_f
 *   Trust region conjugate gradient (single precision).
 *
 * Description:
 *   This function implements reverse communication trust region linear
 *   conjugate gradient with optional predconditioning.  This version is for
 *   single precision variables, see <opk_trcg_f()> for a double precision
 *   version and for explanations.
 *
 * See Also:
 *   <opk_trcg_d()>, <opk_plcg_d()>.
 */
extern void
opk_trcg_f(opk_index_t n, float p[], const float q[],
           float r[], float x[], const float z[],
           float delta, float rho[5],
           opk_cg_state_t *state);

/*---------------------------------------------------------------------------*/

OPK_END_C_DECLS

#endif /* _OPTIMPACK_LINALG_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * fill-column: 79
 * coding: utf-8
 * ispell-local-dictionary: "american"
 * End:
 */
