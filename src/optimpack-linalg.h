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

/**
 * @addtogroup LinearAlgebra
 *
 * This part of the library provides simple versions of some BLAS/LINPACK/LAPACK
 * linear algebra routines.  This routines are intended for *small* problems.
 *
 * Following BLAS, LINPACK and LAPACK conventions, the prefixes `s` and `d` are
 * used to distinguish between the single/double precision floating point
 * versions of some numerical functions.
 *
 * For large scale problems, variables are stored in so-called *vectors* which
 * have nothing to do with the moderate size arrays targeted by the functions
 * of this package.  Prefix 'v' is used to disantangle between the routines of
 * the different packages.  For instance, opk_ddot() and opk_sdot() are the
 * functions to compute the dot product beween two vectors (respectively with
 * double and single precision floating point elements) stored in conventional
 * memory, while opk_vdot() computes the dot product of two vectors (of the
 * same vector space) that may be stored in any suitable way.
 */

/*---------------------------------------------------------------------------*/
/* BLAS CONSTANTS  */

/**
 * @defgroup BLASConstants      BLAS Constants
 * @ingroup LinearAlgebra
 * @{
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

/** Matrix storage order. */
typedef enum _opk_blas_order {
  OPK_BLAS_ROW_MAJOR  = 101, /**< Matrix elements stored in row major order. */
  OPK_BLAS_COL_MAJOR  = 102  /**< Matrix elements stored in column major order. */
} opk_blas_order_t;

/** Matrix transpose. */
typedef enum _opk_blas_trans {
  OPK_BLAS_NO_TRANS   = 111, /**< Do not transpose matrix. */
  OPK_BLAS_TRANS      = 112, /**< Transpose matrix. */
  OPK_BLAS_CONJ_TRANS = 113  /**< Conjugate transpose. */
} opk_blas_trans_t;

/** Part of a triangular matrix to use. */
typedef enum _opk_blas_uplo {
  OPK_BLAS_UPPER = 121, /**< Use upper triangular part. */
  OPK_BLAS_LOWER = 122  /**< Use lower triangular part. */
} opk_blas_uplo_t;

/** Is diagonal unit? */
typedef enum _opk_blas_diag {
  OPK_BLAS_NON_UNIT = 131, /**< Non unit diagonal. */
  OPK_BLAS_UNIT     = 132  /**< Unit diagonal. */
} opk_blas_diag_t;

/** Side of matrix multiply. */
typedef enum _opk_blas_side {
  OPK_BLAS_LEFT  = 141, /**< Left multiply. */
  OPK_BLAS_RIGHT = 142  /**< Right multiply. */
} opk_blas_side_t;

#endif /* HAVE_CBLAS */

/** @} */

/*---------------------------------------------------------------------------*/
/* VECTOR OPERATIONS */

/**
 * @defgroup BLAS1      Level 1 BLAS-like routines
 * @ingroup LinearAlgebra
 * @{
 */

/**
 * Maximum absolute value of a double precision vector.
 *
 * Get maximum absolute value of a vector.  This is also the infinite norm of
 * the vector.
 *
 * @param n    - The number of elements to consider in `x`.
 * @param x    - A vector of `n*incx` values.
 * @param incx - The index increment (greater or equal 1).
 *
 * @return The maximum absolute value among the `n` elements of `x`.
 *
 * @see opk_samax, opk_idamax.
 */
extern double
opk_damax(opk_index_t n, const double x[],
          opk_index_t incx);

/**
 * Maximum absolute value of a single precision vector.
 *
 * @see opk_damax, opk_isamax.
 */
extern float
opk_samax(opk_index_t n, const float x[],
          opk_index_t incx);

/**
 * Sum of the absolute values of a vector.
 *
 * This function computes the sum of the absolute values of a vector.  This is
 * also the L-1 norm of the vector.
 *
 * @param n    - The number of elements to consider in `x`.
 * @param x    - A vector of `n*incx` values.
 * @param incx - The index increment (greater or equal 1).
 *
 * @return Returns the sum of the absolute values of `x`.
 *
 * @see opk_sasum.
 */
extern double
opk_dasum(opk_index_t n,
          const double x[], opk_index_t incx);

/**
 * Sum of the absolute values of a vector.
 *
 * @see opk_dasum.
 */
extern float
opk_sasum(opk_index_t n,
          const float x[], opk_index_t incx);

/**
 * Linear combination of two vectors.
 *
 * This function stores `a*x[ix] + y[iy]` into `y[iy]`. The index increments
 * may be negative to consider elements in reverse order.  The code is
 * optimized for `a = +/-1` and `a = 0`.
 *
 * @param n    - The number of elements to consider in `x` and `y`.
 * @param a    - The scalar factor.
 * @param x    - A vector of `n*|incx|` values.
 * @param incx - The index increment for `x`.
 * @param y    - A vector of `n*|incy|` values.
 * @param incy - The index increment for `y`.
 *
 * @see opk_saxpy.
 */
extern void
opk_daxpy(opk_index_t n, double a,
          const double x[], opk_index_t incx,
          double y[], opk_index_t incy);

/**
 * Linear combination of two vectors.
 *
 * @see opk_daxpy.
 */
extern void
opk_saxpy(opk_index_t n, float a,
          const float x[], opk_index_t incx,
          float y[], opk_index_t incy);

/**
 * Copy a vector into another one.
 *
 * This function copies `n` elements of vectors `x` to vector `y`.  The index
 * increments may be negative to consider elements in reverse order.
 *
 * @param n    - The number of elements to consider in `x` and `y`.
 * @param x    - A vector of `n*|incx|` values.
 * @param incx - The index increment for `x`.
 * @param y    - A vector of `n*|incy|` values.
 * @param incy - The index increment for `y`.
 *
 * @see opk_scopy, opk_dswap.
 */
extern void
opk_dcopy(opk_index_t n,
          const double x[], opk_index_t incx,
          double y[], opk_index_t incy);

/**
 * Copy a vector into another one.
 *
 * @see opk_dcopy, opk_sswap.
 */
extern void
opk_scopy(opk_index_t n,
          const float x[], opk_index_t incx,
          float y[], opk_index_t incy);

/**
 * Dot product of two vectors.
 *
 * This function computes the dot product of two vectors. The index increments
 * may be negative to consider elements in reverse order.
 *
 * @param n    - The number of elements to consider in `x` and `y`.
 * @param x    - A vector of `n*|incx|` values.
 * @param incx - The index increment for `x`.
 * @param y    - A vector of `n*|incy|` values.
 * @param incy - The index increment for `y`.
 *
 * @return The dot product of `x` by `y`.
 *
 * @see opk_sdot.
 */
extern double
opk_ddot(opk_index_t n,
         const double x[], opk_index_t incx,
         const double y[], opk_index_t incy);

/**
 * Dot product of two vectors.
 *
 * @see opk_ddot.
 */
extern float
opk_sdot(opk_index_t n,
         const float x[], opk_index_t incx,
         const float y[], opk_index_t incy);

/**
 * Euclidean norm of a vector.
 *
 * This function computes the Euclidean norm of a vector, avoiding overflows.
 * This is also the L-2 norm of the vector.
 *
 * @param n    - The number of elements to consider in `x`.
 * @param x    - A vector of `n*incx` values.
 * @param incx - The index increment (greater or equal 1).
 *
 * @return The Euclidean norm of vector `x`:
 *
 *     sqrt(x[0]*x[0] + x[incx]*x[incx] + x[2*incx]*x[2*incx] + ...).
 *
 * @see opk_snrm2.
 */
extern double
opk_dnrm2(opk_index_t n,
          const double x[], opk_index_t incx);

/**
 * Euclidean norm of a vector.
 *
 * @see opk_dnrm2.
 */
extern float
opk_snrm2(opk_index_t n,
          const float x[], opk_index_t incx);

/**
 * Scaling of a vector.
 *
 * This function scales a vector by a scalar, the operation is done in place.
 * The `n` elements of `x` get multiplied by `a`.  Does nothing if `n` or
 * `incx` are less than 1.  The code is optimized for scalar `a = +/-1` and
 * `a = 0`.
 *
 * Parameters:
 * @param n    - The number of elements to consider in `x`.
 * @param a    - The scalar factor.
 * @param x    - A vector of `n*incx` values.
 * @param incx - The index increment (greater or equal 1).
 *
 * @see opk_sscal.
 */
extern void
opk_dscal(opk_index_t n, double a,
          double x[], opk_index_t incx);

/**
 * Scaling of a vector.
 *
 * @see opk_dscal.
 */
extern void
opk_sscal(opk_index_t n, float a,
          float  x[], opk_index_t incx);

/**
 * Sum of the values of a vector.
 *
 * This function computes the sum of the values of a vector.
 *
 * @param n    - The number of elements to consider in `x`.
 * @param x    - A vector of `n*incx` values.
 * @param incx - The index increment (greater or equal 1).
 *
 * @return Returns the sum of the values of `x`.
 *
 * @see opk_ssum, opk_dasum.
 */
extern double
 opk_dsum(opk_index_t n,
          const double x[], opk_index_t incx);

/**
 * Sum of the absolute values of a vector.
 *
 * @see opk_dsum, opk_sasum.
 */
extern float
opk_ssum(opk_index_t n,
         const float x[], opk_index_t incx);

/**
 * Exchanging contents of two vectors.
 *
 * This function exchanges `n` elements of vectors `x` and `y`.  The index
 * increments may be negative to consider elements in reverse order.
 *
 * @param n    - The number of elements to consider in `x` and `y`.
 * @param x    - A vector of `n*|incx|` values.
 * @param incx - The index increment for `x`.
 * @param y    - A vector of `n*|incy|` values.
 * @param incy - The index increment for `y`.
 *
 * @see opk_sswap, opk_dcopy.
 */
extern void
opk_dswap(opk_index_t n,
          double x[], opk_index_t incx,
          double y[], opk_index_t incy);

/**
 * Exchanging contents of two vectors.
 *
 * @see opk_dswap, opk_scopy.
 */
extern void
opk_sswap(opk_index_t n,
          float x[], opk_index_t incx,
          float y[], opk_index_t incy);

/**
 * Fill an array with zeros.
 *
 * This function fills the `n` elements of vector `x` by step of `incx` with
 * zeros.
 *
 * @param n    - The number of elements to consider in `x`.
 * @param x    - A vector of `n*incx` values.
 * @param incx - The index increment (greater or equal 1).
 *
 * @see opk_szero.
 */

extern void
opk_dzero(opk_index_t n,
          double x[], opk_index_t incx);

/**
 * Fill an array with zeros.
 *
 * This function fills the `n` elements of vector `x` by step of `incx` with
 * zeros.
 *
 * @param n    - The number of elements to consider in `x`.
 * @param x    - A vector of `n*incx` values.
 * @param incx - The index increment (greater or equal 1).
 *
 * @see opk_dzero.
 */
extern void
opk_szero(opk_index_t n,
          float  x[], opk_index_t incx);

/**
 * Get index of maximum absolute value of a vector.
 *
 * This function returns the index of the maximum absolute value of a vector.
 * Following FORTRAN conventions, indices are 1-based.  If `n = 0` or
 * `incx = 0`, the returned value is 0.
 *
 * @param n    - The number of elements to consider in `x`.
 * @param x    - A vector of `n*incx` values.
 * @param incx - The index increment (greater or equal 1).
 *
 * @return Returns the 1-based index of maximum absolute value of `x`.
 *         Returns 0 if `n` = 0 or `incx` = 0.
 *
 * @see opk_isamax, opk_damax.
 */
extern opk_index_t
opk_idamax(opk_index_t n,
           const double x[], opk_index_t incx);

/**
 * Get index of maximum absolute value of a vector.
 *
 * @see opk_idamax, opk_samax.
 */
extern opk_index_t
opk_isamax(opk_index_t n,
           const float x[], opk_index_t incx);


/** @} */

/*---------------------------------------------------------------------------*/
/* MATRIX-VECTOR OPERATIONS */

/**
 * @defgroup BLAS2      Level 2 BLAS-like routines
 * @ingroup LinearAlgebra
 * @{
 */

/**
 * Matrix-vector operation.
 *
 * This function performs one of the matrix-vector operations
 *
 *     y := alpha A.x + beta y,
 *
 * or
 *
 *     y := alpha A'.x + beta y,
 *
 *  where `alpha` and `beta` are scalars, `x` and `y` are vectors and `A` is
 *  an `m` by `n` matrix.
 *
 * Parameters:
 * @param trans - Transpose flag, specifies the matrix multiplication to be
 *                performed as follows: `OPK_BLAS_NO_TRANS` for <tt>A.x</tt>,
 *                or `OPK_BLAS_TRANS` for <tt>A'.x</tt>.
 * @param m     - The number of rows of matrix `A`.
 * @param n     - The number of columns of matrix `A`.
 * @param alpha - The first scalar parameter.
 * @param a     - Matrix of dimensions `(dla,n)`.
 * @param lda   - The leading dimension of `a`; `lda >= max(1,m)`.
 * @param x     - The source vector of length at least `1+(n-1)*abs(incx)`
 *                when `trans = OPK_BLAS_NO_TRANS` or `1+(m-1)*abs(incx)`
 *                otherwise.
 * @param incx  - The index increment for vector `x` (must be non-zero).
 * @param beta  - The second scalar parameter.
 * @param y     - The destination vector of length at least `1+(m-1)*abs(incy)`
 *                when `trans = OPK_BLAS_NO_TRANS` or `1+(n-1)*abs(incy)`
 *                otherwise.
 * @param incy  - The index increment for vector `y` (must be non-zero).
 *
 * @return Non-zero result `k` means invalid `k`-th argument.
 *
 * @see opk_sgemv.
 */
extern int
opk_dgemv(opk_blas_trans_t trans,
          opk_index_t m, opk_index_t n,
          double alpha,
          const double a[], opk_index_t lda,
          const double x[], opk_index_t incx,
          double beta,
          double y[], opk_index_t incy);

/**
 * Matrix-vector operation.
 *
 * @see opk_dgemv.
 */
extern int
opk_sgemv(opk_blas_trans_t trans,
          opk_index_t m, opk_index_t n,
          float alpha,
          const float a[], opk_index_t lda,
          const float x[], opk_index_t incx,
          float beta,
          float y[], opk_index_t incy);

/**
 * Multiplication of a vector by a triangular matrix.
 *
 * This function performs one of the matrix-vector operations:
 *
 *     x := A.x,
 *
 * or
 *
 *     x := A'.x,
 *
 * where `x` is an `n` element vector and `A` is an `n` by `n` unit, or
 * non-unit, upper or lower triangular matrix.
 *
 * @param uplo  - Specifies whether the matrix is an upper or lower triangular
 *                matrix as follows: `uplo = OPK_BLAS_UPPER` if `A` is an upper
 *                triangular matrix, and `uplo = OPK_BLAS_LOWER` if `A` is a
 *                lower triangular matrix.
 * @param trans - Specifies the operation to be performed as follows: if
 *                `trans = OPK_BLAS_NO_TRANS`, then <tt>x := A.x</tt>; otherwise
 *                if `trans = OPK_BLAS_TRANS` or `trans = OPK_BLAS_CONJ_TRANS`,
 *                then <tt>x := A'.x</tt>.
 * @param diag  - Specifies whether or not `A` is unit triangular as follows: if
 *                `diag = OPK_BLAS_UNIT`, `A` is assumed to be unit triangular;
 *                else if `diag = OPK_BLAS_NON_UNIT`, `A` is not assumed to be
 *                unit triangular.
 * @param n     - The order of the matrix `A`, must be at least zero.
 * @param a     - An array of dimension `lda` by `n`.  With `uplo =
 *                OPK_BLAS_UPPER`, the leading `n` by `n` upper triangular part
 *                of the array `A` must contain the upper triangular matrix and
 *                the strictly lower triangular part of `A` is not used.  With
 *                `uplo = OPK_BLAS_LOWER`, the leading `n` by `n` lower
 *                triangular part of the array `A` must contain the lower
 *                triangular matrix and the strictly upper triangular part of
 *                `A` is not used.  Note that when `diag = OPK_BLAS_UNIT`, the
 *                diagonal elements of `A` are not referenced either, but are
 *                assumed to be unity.
 * @param lda   - The leading dimension of `A`, must be at least `max(n,1)`.
 * @param x     - An array of at least `(1+(n-1)*abs(incx)` elements.  Before
 *                entry, the incremented array must contain the `n` element vector
 *                `x`. On exit, it is overwritten with the tranformed vector `x`.
 * @param incx  - The increment for the elements of `x`, must not be zero.
 *
 * @result Non-zero result `k` means invalid `k`-th argument.
 *
 * @see opk_strmv, opk_dtrsv.
 */
extern int
opk_dtrmv(opk_blas_uplo_t uplo,
          opk_blas_trans_t trans,
          opk_blas_diag_t diag,
          opk_index_t n,
          const double a[], opk_index_t lda,
          double x[], opk_index_t incx);

/**
 * Multiplication of a vector by a triangular matrix.
 *
 * @see opk_dtrmv, opk_strsv.
 */
extern int
opk_strmv(opk_blas_uplo_t uplo,
          opk_blas_trans_t trans,
          opk_blas_diag_t diag,
          opk_index_t n,
          const float a[], opk_index_t lda,
          float x[], opk_index_t incx);

/**
 * Solves a tringular linear system of equations.
 *
 * This function solves one of the matrix-vector operations:
 *
 *     A.x = b,
 *
 * or
 *
 *     A'.x = b,
 *
 * where `x` and `b` are `n` element vectors and `A` is an `n` by `n` unit, or
 * non-unit, upper or lower triangular matrix.
 *
 * No test for singularity or near-singularity is included in this routine.
 * Such tests must be performed before calling this routine.
 *
 * @param uplo  - Specifies whether the matrix is an upper or lower triangular
 *                matrix as follows: `uplo = OPK_BLAS_UPPER` if `A` is an upper
 *                triangular matrix, and `uplo = OPK_BLAS_LOWER` if `A` is a
 *                lower triangular matrix.
 * @param trans - Specifies the operation to be performed as follows: if
 *                `trans = OPK_BLAS_NO_TRANS`, then solves <tt>A.x = b</tt>;
 *                otherwise if `trans = OPK_BLAS_TRANS` or
 *                `trans = OPK_BLAS_CONJ_TRANS`, then solves <tt>A'.x = b</tt>.
 * @param diag  - Specifies whether or not `A` is unit triangular as follows: if
 *                `diag = OPK_BLAS_UNIT`, `A` is assumed to be unit triangular;
 *                else if `diag = OPK_BLAS_NON_UNIT`, `A` is not assumed to be
 *                unit triangular.
 * @param n     - The order of the matrix `A`, must be at least zero.
 * @param a     - An array of dimension `lda` by `n`.  With `uplo =
 *                OPK_BLAS_UPPER`, the leading `n` by `n` upper triangular part
 *                of the array `A` must contain the upper triangular matrix and
 *                the strictly lower triangular part of `A` is not used.  With
 *                `uplo = OPK_BLAS_LOWER`, the leading `n` by `n` lower
 *                triangular part of the array `A` must contain the lower
 *                triangular matrix and the strictly upper triangular part of
 *                `A` is not used.  Note that when `diag = OPK_BLAS_UNIT`, the
 *                diagonal elements of `A` are not referenced either, but are
 *                assumed to be unity.
 * @param lda   - The leading dimension of `A`, must be at least `max(n,1)`.
 * @param x     - An array of at least `(1+(n-1)*abs(incx)` elements.  Before
 *                entry, the incremented array `x` must contain the `n`
 *                element vector `b`. On exit, `x` is overwritten with the
 *                solution vector `x`.
 * @param incx  - The increment for the elements of `x`, must not be zero.
 *
 * @result Non-zero result `k` means invalid `k`-th argument.
 *
 * @see opk_strsv, opk_dtrmv.
 */
extern int
opk_dtrsv(opk_blas_uplo_t uplo,
          opk_blas_trans_t trans,
          opk_blas_diag_t diag,
          opk_index_t n,
          const double a[], opk_index_t lda,
          double x[], opk_index_t incx);

/**
 * Solves a tringular linear system of equations.
 *
 * @see opk_dtrsv, opk_strmv.
 */
extern int
opk_strsv(opk_blas_uplo_t uplo,
          opk_blas_trans_t trans,
          opk_blas_diag_t diag,
          opk_index_t n,
          const float a[], opk_index_t lda,
          float x[], opk_index_t incx);

/** @} */

/*---------------------------------------------------------------------------*/
/* MATRIX-MATRIX OPERATIONS */

/**
 * @defgroup BLAS3      Level 3 BLAS-like routines
 * @ingroup LinearAlgebra
 *
 * These routines are a sub-set of the Level 3 BLAS routines for matrix-matrix
 * operataions.
 *
 * @{
 */

/**
 * Performs a matrix-matrix operation.
 *
 * This function performs one of the matrix-matrix operations
 *
 *      C := alpha op(A).op(B) + beta C,
 *
 * where <tt>op(X) = X</tt> or <tt>op(X) = X'</tt>, `alpha` and `beta` are
 * scalars, and `A`, `B` and `C` are matrices, with `op(A)` an `m` by `k`
 * matrix, `op(B)` a `k` by `n` matrix and `C` an `m` by `n` matrix.
 *
 * @param transa - Specifies the form of `op(A)` to be used in the matrix
 *                 multiplication as follows: if `transa` = OPK_BLAS_NO_TRANS,
 *                 then `op(A)` = `A`; otherwise if `trans = OPK_BLAS_TRANS` or
 *                 `OPK_BLAS_CONJ_TRANS`, then <tt>op(A) = A'</tt>.
 * @param transb - Specifies the form of `op(B)` to be used in the matrix
 *                 multiplication; see `transa`.
 * @param m      - The number of rows of the matrix `op(A)` and of the matrix
 *                 `C`, must be at least zero.
 * @param n      - The number of columns of the matrix `op(B)` and of the matrix
 *                 `C`, must be at least zero.
 * @param k      - The number of columns of the matrix `op(A)` and the number of
 *                 rows of the matrix `op(B)`, must be at least zero.
 * @param alpha  - The scalar `alpha`.
 * @param a      - An array of dimension `lda` by `ka`, where `ka` is `k` when
 *                 `transa = OPK_BLAS_NO_TRANS`, and is `m` otherwise.  Before
 *                 entry with `transa = OPK_BLAS_NO_TRANS`, the leading `m` by
 *                 `k` part of the array `a` must contain the matrix `A`,
 *                 otherwise the leading `k` by `m` part of the array `a` must
 *                 contain the matrix `A`.
 * @param lda    - The first dimension of `a`. When `transa = OPK_BLAS_NO_TRANS`
 *                 then `lda` must be at least `max(m,1)`, otherwise `lda` must
 *                 be at least `max(k,1)`.
 * @param b      - An array of dimension `ldb` by `kb`, where `kb` is `n` when
 *                 `transb = OPK_BLAS_NO_TRANS`, and is `k` otherwise.  Before
 *                 entry with `transb = OPK_BLAS_NO_TRANS`, the leading `k` by
 *                 `n` part of the array `b` must contain the matrix `B`,
 *                 otherwise the leading `n` by `k` part of the array `b` must
 *                 contain the matrix `B`.
 * @param ldb    - The first dimension of `b`. When `transb = OPK_BLAS_NO_TRANS`
 *                 then `ldb` must be at least `max(k,1)`, otherwise `ldb` must
 *                 be at least `max(n,1)`.
 * @param beta   - The scalar `beta`.  When `beta` is supplied as zero then `C`
 *                 need not be set on input.
 * @param c      - An array of dimension `ldc` by `n`.  Before entry, the
 *                 leading `m` by `n` part of the array `c` must contain the
 *                 matrix `C`, except when `beta` is zero, in which case `c`
 *                 need not be set on entry.  On exit, the array `c` is
 *                 overwritten by the `m` by `n` matrix
 *                 <tt>alpha op(A).op(B) + beta C</tt>.
 * @param ldc    - The first dimension of `c`. `ldc` must be at least
 *                 `max(m,1)`.
 *
 * @return Non-zero result `k` means invalid `k`-th argument.
 *
 * @see opk_sgemm.
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

/**
 * Performs a matrix-matrix operation.
 *
 * @see opk_dgemm.
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

/**
 * Performs a symmetric rank k operation.
 *
 * This function performs one of the symmetric rank `k` operations
 *
 *     C := alpha A.A' + beta C,
 *
 * or
 *
 *     C := alpha A'.A + beta C,
 *
 * where `alpha` and `beta` are scalars, `C` is an `n` by `n` symmetric matrix
 * and `A` is an `n` by `k` matrix in the first case and a `k` by `n` matrix in
 * the second case.
 *
 * @param uplo  - Specifies whether the matrix is the upper or lower triangular
 *                part of the array `C` is to be referenced as follows.  If
 *                `uplo = OPK_BLAS_UPPER`, only the upper triangular part of
 *                `C` is to be referenced.  If `uplo = OPK_BLAS_LOWER`, only
 *                the lower triangular part of `C` is to be referenced.
 * @param trans - Specifies the operation to be performed as follows:
 *                <tt>C := alpha A.A' + beta C</tt>, if
 *                `trans = OPK_BLAS_NO_TRANS`; or
 *                <tt>C := alpha A'.A + beta C</tt>,
 *                if `trans = OPK_BLAS_TRANS` or `OPK_BLAS_CONJ_TRANS`.
 * @param n     - The order of the matrix `C`; `N` must be at least zero.
 * @param k     - On entry with `trans = OPK_BLAS_NO_TRANS`, `k` specifies the
 *                number of columns of the matrix `A`, and on entry with
 *                `trans` = OPK_BLAS_TRANS` or `OPK_BLAS_CONJ_TRANS`, `k`
 *                specifies the number of rows of the matrix `A`.  `k` must be
 *                at least zero.
 * @param alpha - The scalar `alpha`.
 * @param a     - Real array of dimension `lda` by `ka`, where `ka` is `k`
 *                when `trans = OPK_BLAS_NO_TRANS`, and is `n` otherwise.
 *                Before entry with `trans = OPK_BLAS_NO_TRANS`, the leading
 *                `n` by `k` part of the array `a` must contain the matrix `A`,
 *                otherwise the leading `k` by `n` part of the array `a` must
 *                contain the matrix `A`.
 * @param lda   - The leading dimension of `a`.  When
 *                `trans = OPK_BLAS_NO_TRANS`, then `lda` must be at least
 *                `max(n,1)`, otherwise `LDA` must be at least `max(k,1)`.
 * @param beta  - The scalar `beta`.
 * @param c     - Real array of dimension `ldc` by `n`.
 *                Before entry with `uplo = OPK_BLAS_UPPER`, the leading `n` by
 *                `n` upper triangular part of the array `c` must contain the
 *                upper triangular part of the symmetric matrix and the
 *                strictly lower triangular part of `c` is not referenced.  On
 *                exit, the upper triangular part of the array `c` is
 *                overwritten by the upper triangular part of the updated
 *                matrix.  Before entry with `uplo = OPK_BLAS_LOWER`, the
 *                leading `n` by `n` lower triangular part of the array `c`
 *                must contain the lower triangular part of the symmetric
 *                matrix and the strictly upper triangular part of `c` is not
 *                referenced.  On exit, the lower triangular part of the array
 *                `c` is overwritten by the lower triangular part of the
 *                updated matrix.
 * @param ldc   - The leading dimension of `c`.  `ldc` must be at least
 *                `max(n,1)`.
 *
 * @return Non-zero result `k` means invalid `k`-th argument.
 *
 * @see opk_ssyrk.
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

/**
 * Performs a symmetric rank k operation.
 *
 * @see opk_dsyrk.
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

/** @} */

/*---------------------------------------------------------------------------*/
/* LAPACK-LIKE ROUTINES */

/**
 * @defgroup LAPACK      LAPACK-like routines
 * @ingroup LinearAlgebra
 *
 * This part of the library implements the non-blocked version of a few LINPACK
 * or LAPACK routines needed by other part of the softaware.  These routines
 * are not intended for large scale problems.
 *
 * @{
 */

/**
 * Cholesky factorization of a real symmetric positive definite matrix.
 *
 * This function computes the Cholesky factorization of a real symmetric
 * positive definite matrix `A`.
 *
 * The factorization has the form <tt>A = U'.U</tt>, if <tt>uplo =
 * OPK_BLAS_UPPER</tt>, or <tt>A = L.L'</tt>, if <tt>uplo = OPK_BLAS_LOWER</tt>
 * where `U` is an upper triangular matrix and `L` is lower triangular.
 *
 * This is the unblocked version of the algorithm, calling Level 2 BLAS.
 *
 * @param uplo - Specifies whether the upper or lower triangular part of the
 *               symmetric matrix `A` is stored. `uplo` is one of:
 *               `OPK_BLAS_UPPER` or `OPK_BLAS_LOWER`.
 * @param n    - The order of the matrix `A`, must be greater of equal zero.
 * @param a    - Input/output array of dimension `lda` by `n`.  On entry, `a`
 *               contains the symmetric matrix `A`.  If <tt>uplo =
 *               OPK_BLAS_UPPER</tt>, the leading `n` by `n` upper triangular
 *               part of `a` contains the upper triangular part of the matrix
 *               `A`, and the strictly lower triangular part of `a` is not
 *               referenced.  If <tt>uplo = OPK_BLAS_LOWER</tt>, the leading
 *               `n` by `n` lower triangular part of `a` contains the lower
 *               triangular part of the matrix `A`, and the strictly upper
 *               triangular part of `a` is not referenced.  On successful exit
 *               (with returned value 0), `a` contains the factor `U` or `L`
 *               from the Cholesky factorization <tt>A = U'.U</tt> or <tt>A =
 *               L.L'</tt>.
 * @param lda  - The leading dimension of the array `a`. <tt>lda >=
 *               max(n,1)</tt>.
 *
 * @return
 *   The returned value, says `info` is:
 *   * `info = 0` - successful exit.
 *   * `info < 0` - the `k`-th argument had an illegal value, with
 *                  `k = -info`.
 *   * `info > 0` - the leading minor of order `k` is not positive definite,
 *                  and the factorization could not be completed, with
 *                  `k = info`.
 *
 * @see
 *   opk_spotf2.
 */
extern opk_index_t
opk_dpotf2(opk_blas_uplo_t uplo,
           opk_index_t n,
           double a[], opk_index_t lda);

/**
 * Cholesky factorization of a real symmetric positive definite matrix.
 *
 * @see opk_dpotf2.
 */
extern opk_index_t
opk_spotf2(opk_blas_uplo_t uplo,
           opk_index_t n,
           float a[], opk_index_t lda);

/** @} */

/*---------------------------------------------------------------------------*/
/* LINEAR CONJUGATE GRADIENT */

/**
 * @defgroup LCG      Linear conjugate gradient
 * @ingroup LinearAlgebra
 *
 * This package provides methods to solve a linear system of equations with a
 * positive definite left hand-side matrix by means of iterative conjugate
 * gradient method (Hestenes & Stiefel, 1952).  A reverse communication version
 * of the linear conjugate gradient with and without a preconditioner is
 * implemented.
 *
 * For computing trust region-like step or if left hand-side matrix is not
 * positive definite, the truncated conjugate gradient method of Steihaug
 * (1983) is provided.
 *
 *
 * ### References
 *
 * 1. M.R. Hestenes & E. Stiefel, "Methods of Conjugate Gradients for Solving
 *    Linear Systems," Journal of Research of the National Bureau of Standards
 *    **49**, 409-436 (1952).
 *
 * 2. T. Steihaug, "The conjugate gradient method and trust regions in large
 *    scale optimization," SIAM Journal on Numerical Analysis **20**, 626-637
 *    (1983).
 *
 * @{
 */

/**
 * Possible values of conjugate gradient state.
 */
typedef enum {
  OPK_CG_ERROR      = -1, /**< error */
  OPK_CG_START      =  0, /**< start with no initial variable (all zero) */
  OPK_CG_RESTART    =  1, /**< start or restart with initial variables */
  OPK_CG_NEWX       =  2, /**< new step taken */
  OPK_CG_PRECOND    =  3, /**< caller has to to compute z = Q.r */
  OPK_CG_AP         =  4, /**< caller has to to compute q = A.p */
  OPK_CG_FINISH     =  5, /**< convergence or rounding errors prevent further
                               progress */
  OPK_CG_NON_CONVEX =  6, /**< non-positive definitiveness has been detected */
  OPK_CG_TRUNCATED  =  7  /**< step truncated at trust region boundary */
} opk_cg_state_t;

/**
 * Preconditioned linear conjugate gradient (double precision).
 *
 * Iteratively solve a linear system <tt>A.x = b</tt> where `A` is a symmetric
 * positive definite matrix by preconditioned conjugate gradient method with
 * reverse communication.
 *
 * @param n     - The number of variables.
 * @param p     - A work array of size `n`; see parameter `state` for
 *                explanations.
 * @param q     - A work array of size `n`; on return with `state` set to
 *                `OPK_CG_AP`, the caller must compute the matrix vector
 *                multiplication <tt>q=A.p</tt> (i.e. store into vector `q` the
 *                result of multiplying vector `p` by the matrix `A`); on
 *                return with `state` set to `OPK_CG_NEWX`, `q` is filled with
 *                the unscaled residuals step.
 * @param r     - An array of size `n`; on initialization or restart (with
 *                `state` set to `OPK_CG_START` or `OPK_CG_RESTART`), `r`
 *                stores the right hand side vector `b`; otherwise, `r` stores
 *                the current residuals.
 * @param x     - An array of size `n` to store the unknowns.  On
 *                initialization with given parameter values or restart,
 *                i.e. with `state` set to `OPK_CG_RESTART`, `x` stores the
 *                initial unknowns `x0`; on initial entry with `state` set to
 *                `OPK_CG_START`, the contents of `x` is ignored (`x` will be
 *                filled with zeros); otherwise, `x` stores the current
 *                solution.
 * @param z     - Optional work array of size `n`; on return with `state` set
 *                to `OPK_CG_PRECOND`, the caller must compute and store the
 *                preconditioned residuals into `z` given the current residuals
 *                which are in `r`.  If no preconditioning is to be used, `z`
 *                must be set to `NULL` (in which case, a `state` with value
 *                `OPK_CG_PRECOND` will never get returned).  If not `NULL`,
 *                `z` can be the same array as `q` to save memory.
 * @param rho   - A 4-element array used to store the dot product <tt>r'.z</tt>
 *                of the current and previous iterations and to store some
 *                parameters of the algorithm.  <tt>rho[0]</tt> is the former
 *                value of <tt>rho[1]</tt>; <tt>rho[1]</tt> is the dot product of
 *                the residuals `r` by the (preconditioned) residuals `z`;
 *                <tt>rho[2]</tt> is the optimal step size (`alpha`);
 *                <tt>rho[3]</tt> is the weight of the former search direction
 *                (`beta`).
 * @param state - Address of integer variable to store the current stage of the
 *                algorithm (must not be altered by the caller between calls
 *                except to restart the algorithm).  On initialization or
 *                restart of the conjugate gradient iterations, `state` must be
 *                set to `OPK_CG_START` or `OPK_CG_RESTART` and `r` must be set
 *                to the right hand side vector, that is: <tt>r = b</tt>; if
 *                `state` is `OPK_CG_START`, `x` is filled with zeros by the
 *                PLCG routine; if `state` is `OPK_CG_RESTART`, `x` must be set
 *                with some initial values of the parameters, that is: <tt>x
 *                = x0</tt>. Upon return of the function, the value stored in
 *                `state` is one of:
 *
 * - `OPK_CG_NEWX`: `x` stores the current solution, `r` stores the current
 *   residuals, the change in parameters is equal to <tt>rho[2]*p</tt> and the
 *   corresponding residuals change is equal to <tt>-rho[2]*q</tt>.  The caller
 *   can decide to terminate the iterations, to pursue the iterations or to
 *   restart the algorithm.
 *
 * - `OPK_CG_AP`: The caller must compute the matrix vector multiplication
 *   <tt>q = A.p</tt> (*i.e.*, store into vector `q` the result of multiplying
 *   vector `p` by the matrix `A`) and call the conjugate gradient routine
 *   again.
 *
 * - `OPK_CG_PRECOND`: The caller must compute and store the preconditioned
 *   residuals into `z` (given the current residuals which are in `r`) and call
 *   the conjugate gradient routine again.
 *
 * - `OPK_CG_FINISH`: No further progress are possible either because the exact
 *   solution has been found or because rounding errors prevent further
 *   progress.
 *
 * - `OPK_CG_ERROR`: An error has been detected by the algorithm: invalid
 *   parameter or non-positive definitiveness of the matrix `A` or of the
 *   preconditioner.
 *
 * For consistency reasons, `OPK_CG_NEWX` is returned after internal
 * initialization, *i.e.*, `x` is just `x0` and `r` is just
 * <tt>r = b - A.x0</tt>.
 *
 *
 * ### Example
 *
 * In the following piece of code, we assume that the arrays `a` and `c` store
 * respectively the coefficients of the matrix `A` and of the preconditioner.
 *
 *     for (i = 0; i < n; ++i) {
 *       x[i] = ...; // initial set of parameters
 *     }
 *     for (i = 0; i < n; ++i) {
 *       r[i] = b[i];        // copy b into r
 *     }
 *     state = OPK_CG_RESTART;
 *     for (;;) {
 *       opk_dplcg(n, p, q, r, x, z, rho, &state);
 *       if (state == OPK_CG_AP) {
 *         for (i = 0; i < n; ++i) {
 *           double s = 0.0;
 *           for (j = 0; j < n; ++j) {
 *             s += a[i*n + j]*p[j];
 *           }
 *           q[i] = s;
 *         }
 *       } else if (state == OPK_CG_PRECOND) {
 *         for (i = 0; i < n; ++i) {
 *           double s = 0.0;
 *           for (j = 0; j < n; ++j) {
 *             s += c[i*n + j]*r[j];
 *           }
 *           z[i] = s;
 *         }
 *       } else if (state == OPK_CG_NEWX) {
 *         double s = 0.0;
 *         for (i = 0; i < n; ++i) {
 *           s += r[i]*r[i];
 *         }
 *         if (s < 1e-3) {
 *           break; // convergence
 *         }
 *       } else {
 *         // unexpected (must be an error)
 *         fprintf(stderr, "error in PLCG\n");
 *         break;
 *       }
 *     }
 *
 * @see opk_dlcg, opk_splcg.
 */
extern void
opk_dplcg(opk_index_t n, double p[], double q[],
          double r[], double x[], double z[],
          double rho[4], opk_cg_state_t *state);

/**
 * Preconditioned linear conjugate gradient (single precision).
 *
 * @see opk_slcg(), opk_dplcg() for the meaning of the arguments.
 */
extern void
opk_splcg(opk_index_t n, float p[], float q[],
          float r[], float x[], float z[],
          float rho[4], opk_cg_state_t *state);

/**
 * Linear conjugate gradient (double precision).
 *
 * This function implements reverse communication linear conjugate gradient
 * without predconditioner for double precision floating point.  See
 * opk_dplcg() for the meaning of the arguments and opk_slcg() for a single
 * precision version.
 *
 * @see opk_slcg, opk_dplcg.
 */
extern void
opk_dlcg(opk_index_t n, double p[], double q[], double r[],
         double x[], double rho[4], opk_cg_state_t *state);

/**
 * Linear conjugate gradient (single precision).
 *
 * This function implements reverse communication linear conjugate gradient
 * without predconditioner for single precision floating point.  See
 * opk_dplcg() for the meaning of the arguments and opk_dlcg() for a double
 * precision version.
 *
 * @see opk_dlcg(), opk_dplcg().
 */
extern void
opk_slcg(opk_index_t n, float p[], float q[], float r[],
         float x[], float rho[4], opk_cg_state_t *state);

/**
 * Trust region conjugate gradient (double precision).
 *
 * This function implements a reverse communication version of a trust region
 * linear conjugate gradient algorithm with optional preconditioning.  The
 * trust-region (or truncated) conjugate gradient method is due to Steihaug
 * (see References below).  This version is for double precision variables, see
 * opk_strcg() for a single precision version.  See opk_dplcg() for an
 * example of using reverse communication.
 *
 * @param n - The number of variables.
 * @param p - A work array of size `n`; see parameter `state` for explanations.
 * @param q - A work array of size `n`; on return with `state` set to
 *            `OPK_CG_AP`, the caller must compute the matrix vector
 *            multiplication <tt>q=A.p</tt> (i.e., store into vector `q` the
 *            result of multiplying vector `p` by the matrix `A`); on return
 *            with `state` set to `OPK_CG_NEWX`, `q` is filled with the
 *            unscaled residuals step.
 * @param r - An array of size `n`; on initialization or restart (i.e., with
 *            `state` set to `OPK_CG_START` or `OPK_CG_RESTART`), `r` stores
 *            the right hand side vector `b`; otherwise, `r` stores the current
 *            residuals.
 * @param x - An array of size `n` to store the unknowns.  On initialization
 *            with given parameter values or restart, i.e. with `state` set to
 *            `OPK_CG_RESTART`, `x` stores the initial unknowns `x0`; on
 *            initial entry with `state` set to `OPK_CG_START`, the contents of
 *            `x` is ignored (`x` will be filled with zeroes); otherwise, `x`
 *            stores the current solution.
 * @param z - Optional work array of size `n`; on return with `state` set to
 *            `OPK_CG_PRECOND`, the caller must compute and store the
 *            preconditioned residuals into `z` given the current residuals
 *            which are in `r`.  If no preconditioning is to be used, `z` must
 *            be set to `NULL` (in which case, a `state` with value
 *            `OPK_CG_PRECOND` will never get returned).
 * @param delta - The maximum length (Euclidean norm) of vector `x`.  This
 *            value must not be changed once the algorithm is started.
 * @param rho - A 5-element array used to store the dot product <tt>r'.z</tt>
 *            of the current and previous iterations and to store some
 *            parameters of the algorithm.  <tt>rho[0]</tt> is the former value
 *            of <tt>rho[1]</tt>; <tt>rho[1]</tt> is the dot product of the
 *            residuals `r` by the (preconditioned) residuals `z`;
 *            <tt>rho[2]</tt> is the optimal step size (`alpha`);
 *            <tt>rho[3]</tt> is the weight of the former search direction
 *            (`beta`); <tt>rho[4]</tt> is the Euclidean norm of `x`.
 * @param state - Address of integer variable to store the current stage of the
 *            algorithm (must not be altered by the caller between calls except
 *            to restart the algorithm).  On initialization or restart of the
 *            conjugate gradient iterations, `state` must be set to
 *            `OPK_CG_START` or `OPK_CG_RESTART` and `r` must be set to the
 *            right hand side vector, that is: <tt>r=b</tt>; if `state` is
 *            `OPK_CG_START`, `x` is filled with zeroes by the conjugate
 *            gradient routine; if `state` is `OPK_CG_RESTART`, `x` must be set
 *            with some initial values of the parameters, that is:
 *            <tt>x=x0</tt>. Upon return of the function, the value stored in
 *            `state` is one of:
 *
 * - `OPK_CG_NEWX`: `x` stores the current solution, `r` stores the current
 *   residuals, the change in parameters is equal to <tt>rho[2]*p</tt> and the
 *   corresponding residuals change is equal to <tt>-rho[2]*q</tt>.  The caller
 *   can decide to terminate the iterations, to pursue the iterations or to
 *   restart the algorithm.
 *
 * - `OPK_CG_PRECOND`: The caller must compute and store the preconditioned
 *   residuals into `z` (given the current residuals which are in `r`) and call
 *   the conjugate gradient routine again.
 *
 * - `OPK_CG_AP`: The caller must compute the matrix vector multiplication
 *   <tt>q=A.p</tt> (i.e., store into vector `q` the result of multiplying
 *   vector `p` by the matrix `A`) and call the conjugate gradient routine
 *   again.
 *
 * - `OPK_CG_FINISH`: No further progress are possible either because the exact
 *   solution has been found or because rounding errors prevent further
 *   progress.  The Euclidean norm of `x` is less than `delta`.
 *
 * - `OPK_CG_NON_CONVEX`: A further conjugate gradient step cannot be taken
 *   because non positive definitiveness of the preconditioner has been
 *   detected.  The Euclidean norm of `x` is less than `delta`.
 *
 * - `OPK_CG_TRUNCATED`: A truncated gradient step has been taken.  The
 *   Euclidean norm of `x` is equal to `delta`. No further iterations should be
 *   taken (unless with a larger `delta` and `state` set to `OPK_CG_RESTART`).
 *
 * - `OPK_CG_ERROR`: An error has been detected by the algorithm: invalid
 *   parameter, or corrupted workspace `rho`.
 *
 * For consistency reasons, `OPK_CG_NEWX` is returned after internal
 * initialization, i.e., `x` is just `x0` and `r` is just <tt>r = b - A.x0<tt>.
 *
 * It is a bad idea to restart the algorithm with an initial `x` with norm
 * greater than `delta`, in this case, `state` is set to `OPK_CG_TRUNCATED` and
 * the initial `x` is rescaled to have a Euclidean norm equals to `delta`.
 *
 *
 * ### References
 *
 *   - J. Nocedal and S. J. Wright, "Numerical Optimization", Springer Verlag,
 *     2006.
 *
 *   - T. Steihaug, "The conjugate gradient method and trust regions in large
 *     scale optimization", SIAM Journal on Numerical Analysis, vol. **20**,
 *     pp. 626-637, 1983.
 *
 * @see opk_strcg, opk_dplcg.
 */
extern void
opk_dtrcg(opk_index_t n, double p[], const double q[],
          double r[], double x[], const double z[],
          double delta, double rho[5],
          opk_cg_state_t *state);

/**
 * Trust region conjugate gradient (single precision).
 *
 * This function implements reverse communication trust region linear conjugate
 * gradient with optional predconditioning.  This version is for single
 * precision variables, see opk_strcg() for a double precision version and for
 * explanations.
 *
 * @see opk_dtrcg, opk_dplcg.
 */
extern void
opk_strcg(opk_index_t n, float p[], const float q[],
          float r[], float x[], const float z[],
          float delta, float rho[5],
          opk_cg_state_t *state);

/** @} */

/*---------------------------------------------------------------------------*/

OPK_END_C_DECLS

#endif /* _OPTIMPACK_LINALG_H */
