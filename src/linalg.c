/*
 * linalg.c --
 *
 * Linear algebra routines. Based on code from LINPACK and LAPACK libraries. These
 * routines are designed for "small" size problems.
 *
 *----------------------------------------------------------------------------------------
 *
 * The OptimPack library is licensed under the MIT "Expat" License:
 *
 * Copyright (c) 2008-2014: Éric Thiébaut
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this
 * software and associated documentation files (the "Software"), to deal in the Software
 * without restriction, including without limitation the rights to use, copy, modify,
 * merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies
 * or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 *----------------------------------------------------------------------------------------
 */

#ifndef OPK_LINALG_C_
#define OPK_LINALG_C_ 1

#include <string.h>
#include <math.h>

#include "optimpack.h"


/* LOOP over variable (0-based index) */
#define LOOP(VARIABLE, NUMBER) for (VARIABLE=0; VARIABLE<NUMBER; ++VARIABLE)

/* START_0 (START_1) get 0-based (1-based) starting iteration index NUMBER is the number
   of steps to take and STEP is the, possibly negative, increment. */
#define START_0(NUMBER, STEP)    ((STEP) >= 0 ? 0 : (1 - (NUMBER))*(STEP))
#define START_1(NUMBER, STEP)    ((STEP) >= 0 ? 1 : 1 + (1 - (NUMBER))*(STEP))

/* SWAP interchanges two variables VAR1 and VAR2 using temporary variable TMP (must be
   enclosed into a block). */
#define SWAP(VAR1, VAR2, TMP) TMP = VAR1; VAR1 = VAR2; VAR2 = TMP

/* Macros for loop unrolling with 4 or 8 levels of unrolling.
 *
 * The 'a' versions do not use any temporary variables. The 'b' versions use temporary
 * variables: t0, t1, ... which must have been correctly initialized. All these macros
 * require the macro UNROLL_OP to be defined as:
 *
 *     UNROLL_OP(i)    for type 'a' macros
 *     UNROLL_OP(t, i) for type 'b' macros
 *
 * The arguments of the macros are:
 *
 *     I - integer variable for element indexing
 *     M - integer variable
 *     N - number of elements (*must* be > 0)
 */
#define UNROLL_4a(i,m,n)                        \
    do {                                        \
        m = (n & 3);                            \
        for (i = 0; i < m; ++i) {               \
            UNROLL_OP(i);                       \
        }                                       \
        for (i = m; i < n; i += 4) {            \
            UNROLL_OP(i);                       \
            UNROLL_OP(i + 1);                   \
            UNROLL_OP(i + 2);                   \
            UNROLL_OP(i + 3);                   \
        }                                       \
    } while (0)

#define UNROLL_4b(i,m,n)                        \
    do {                                        \
        m = (n & 3);                            \
        for (i = 0; i < m; ++i) {               \
            UNROLL_OP(t0, i);                   \
        }                                       \
        for (i = m; i < n; i += 4) {            \
            UNROLL_OP(t0, i);                   \
            UNROLL_OP(t1, i + 1);               \
            UNROLL_OP(t2, i + 2);               \
            UNROLL_OP(t3, i + 3);               \
        }                                       \
    } while (0)

#define UNROLL_8a(i,m,n)                        \
    do {                                        \
        m = (n & 7);                            \
        for (i = 0; i < m; ++i) {               \
            UNROLL_OP(i);                       \
        }                                       \
        for (i = m; i < n; i += 8) {            \
            UNROLL_OP(i);                       \
            UNROLL_OP(i + 1);                   \
            UNROLL_OP(i + 2);                   \
            UNROLL_OP(i + 3);                   \
            UNROLL_OP(i + 4);                   \
            UNROLL_OP(i + 5);                   \
            UNROLL_OP(i + 6);                   \
            UNROLL_OP(i + 7);                   \
        }                                       \
    } while (0)

#define UNROLL_8b(i,m,n)                        \
    do {                                        \
        m = (n & 7);                            \
        for (i = 0; i < m; ++i) {               \
            UNROLL_OP(t0, i);                   \
        }                                       \
        for (i = m; i < n; i += 8) {            \
            UNROLL_OP(t0, i);                   \
            UNROLL_OP(t1, i + 1);               \
            UNROLL_OP(t2, i + 2);               \
            UNROLL_OP(t3, i + 3);               \
            UNROLL_OP(t4, i + 4);               \
            UNROLL_OP(t5, i + 5);               \
            UNROLL_OP(t6, i + 6);               \
            UNROLL_OP(t7, i + 7);               \
        }                                       \
    } while (0)

/* First pass for single precision version of the routines defined in this file. */
#include "linalg-single.h"
#include __FILE__

/* Second pass for double precision version of the routines defined in this file. */
#include "linalg-double.h"
#include __FILE__


#else /* OPK_LINALG_C_ is defined ------------------------------------------------------*/

/*---[ Level 1 BLAS-like routines ]-----------------------------------------------------*/

/* AMAX */

#ifdef OPK_AMAX
real_t
OPK_AMAX(opk_index n, const real_t x[],
         opk_index incx)
{
    real_t a, b;
    opk_index i, m;

    a = OPK_REALCONST(0.0);
    if (n < 1 || incx < 1) {
        return a;
    }
    b = a;
    if (incx == 1) {
        for (i = 0; i < n; ++i) {
            if (x[i] < a) a = x[i];
            if (x[i] > b) b = x[i];
        }
    } else {
        m = n*incx;
        for (i = 0; i < m; i += incx) {
            if (x[i] < a) a = x[i];
            if (x[i] > b) b = x[i];
        }
    }
    return (((a = -a) >= b) ? a : b);
}
#endif /* OPK_AMAX */

/*--------------------------------------------------------------------------------------*/
/* SWAP */

#ifdef OPK_SWAP
void
OPK_SWAP(opk_index n,
         real_t x[], opk_index incx,
         real_t y[], opk_index incy)
{
    real_t t0, t1, t2, t3;
    opk_index i, m, ix, iy;

    if (n <= 0 || x == y) {
        return;
    }
    if (incx == incy) {
        if (incx == 1) {
#define UNROLL_OP(t,j)  SWAP(x[j], y[j], t)
            UNROLL_4b(i,m,n);
#undef UNROLL_OP
        } else {
            ix = START_0(n, incx);
            for (i = 0; i < n; ++i) {
                SWAP(x[ix], y[ix], t0);
                ix += incx;
            }
        }
    } else if (incx == 1) {
        iy = START_0(n, incy);
        for (i = 0; i < n; ++i) {
            SWAP(x[i], y[iy], t0);
            iy += incy;
        }
    } else if (incy == 1) {
        ix = START_0(n, incx);
        for (i = 0; i < n; ++i) {
            SWAP(x[ix], y[i], t0);
            ix += incx;
        }
    } else {
        ix = START_0(n, incx);
        iy = START_0(n, incy);
        for (i = 0; i < n; ++i) {
            SWAP(x[ix], y[iy], t0);
            ix += incx;
            iy += incy;
        }
    }
}
#endif /* OPK_SWAP */

/*--------------------------------------------------------------------------------------*/
/* SCALE */

#ifdef OPK_SCAL
void
OPK_SCAL(opk_index n, real_t a,
         real_t x[], opk_index incx)
{
    const real_t ZERO = OPK_REALCONST(0.0);
    opk_index i, m;

    if (n <= 0 || a == OPK_REALCONST(1.0)) {
        return;
    }
    if (incx == 1) {
        if (a == ZERO) {
#define UNROLL_OP(i) x[i] = ZERO
            UNROLL_8a(i,m,n);
#undef UNROLL_OP
        } else if (a == OPK_REALCONST(-1.0)) {
#define UNROLL_OP(i) x[i] = -x[i]
            UNROLL_8a(i,m,n);
#undef UNROLL_OP
        } else {
#define UNROLL_OP(i) x[i] *= a
            UNROLL_8a(i,m,n);
#undef UNROLL_OP
        }
    } else if (incx > 0) {
        m = n*incx;
        if (a == ZERO) {
            for (i = 0; i < m; i += incx) {
                x[i] = ZERO;
            }
        } else if (a == OPK_REALCONST(-1.0)) {
            for (i = 0; i < m; i += incx) {
                x[i] = -x[i];
            }
        } else {
            for (i = 0; i < m; i += incx) {
                x[i] *= a;
            }
        }
    }
}
#endif /* OPK_SCAL */

/*--------------------------------------------------------------------------------------*/
/* COPY */

#ifdef OPK_COPY
void
OPK_COPY(opk_index n,
         const real_t x[], opk_index incx,
         real_t y[], opk_index incy)
{
    opk_index i, m, ix, iy;

    if (n <= 0 || (incx == incy && x == y)) {
        return;
    }
    if (incx == 1) {
        if (incy == 1) {
            /* Even if we have `memcpy` we cannot use it in case of overlapping arrays and
               because the copy direction in `memcpy` is unspecified. */
# define UNROLL_OP(i) y[i] = x[i]
            UNROLL_8a(i,m,n);
# undef UNROLL_OP
        } else {
            iy = START_0(n, incy);
            for (ix = 0; ix < n; ++ix) {
                y[iy] = x[ix];
                iy += incy;
            }
        }
    } else {
        if (incy == 1) {
            ix = START_0(n, incx);
            for (iy = 0; iy < n; ++iy) {
                y[iy] = x[ix];
                ix += incx;
            }
        } else {
            ix = START_0(n, incx);
            iy = START_0(n, incy);
            for (i = 0; i < n; ++i) {
                y[iy] = x[ix];
                ix += incx;
                iy += incy;
            }
        }
    }
}
#endif /* OPK_COPY */

/*--------------------------------------------------------------------------------------*/
/* AXPY */

#ifdef OPK_AXPY
void
OPK_AXPY(opk_index n, real_t a,
         const real_t x[], opk_index incx,
         real_t y[], opk_index incy)
{
    opk_index m, i, ix, iy;

    if (n <= 0 || a == OPK_REALCONST(0.0)) {
        return;
    }
    if (a == OPK_REALCONST(1.0)) {
        if (incx == 1) {
            if (incy == 1) {
#define UNROLL_OP(i) y[i] += x[i]
                UNROLL_4a(i,m,n);
#undef UNROLL_OP
            } else {
                iy = START_0(n, incy);
                for (i = 0; i < n; ++i) {
                    y[iy] += x[i];
                    iy += incy;
                }
            }
        } else {
            if (incy == 1) {
                ix = START_0(n, incx);
                for (i = 0; i < n; ++i) {
                    y[i] += x[ix];
                    ix += incx;
                }
            } else {
                ix = START_0(n, incx);
                iy = START_0(n, incy);
                for (i = 0; i < n; ++i) {
                    y[iy] += x[ix];
                    ix += incx;
                    iy += incy;
                }
            }
        }
    } else if (a == OPK_REALCONST(-1.0)) {
        if (incx == 1) {
            if (incy == 1) {
#define UNROLL_OP(i) y[i] -= x[i]
                UNROLL_4a(i,m,n);
#undef UNROLL_OP
            } else {
                iy = START_0(n, incy);
                for (i = 0; i < n; ++i) {
                    y[iy] -= x[i];
                    iy += incy;
                }
            }
        } else {
            if (incy == 1) {
                ix = START_0(n, incx);
                for (i = 0; i < n; ++i) {
                    y[i] -= x[ix];
                    ix += incx;
                }
            } else {
                ix = START_0(n, incx);
                iy = START_0(n, incy);
                for (i = 0; i < n; ++i) {
                    y[iy] -= x[ix];
                    ix += incx;
                    iy += incy;
                }
            }
        }
    } else {
        if (incx == 1) {
            if (incy == 1) {
#define UNROLL_OP(i) y[i] += a*x[i]
                UNROLL_4a(i,m,n);
#undef UNROLL_OP
            } else {
                iy = START_0(n, incy);
                for (i = 0; i < n; ++i) {
                    y[iy] += a*x[i];
                    iy += incy;
                }
            }
        } else {
            if (incy == 1) {
                ix = START_0(n, incx);
                for (i = 0; i < n; ++i) {
                    y[i] += a*x[ix];
                    ix += incx;
                }
            } else {
                ix = START_0(n, incx);
                iy = START_0(n, incy);
                for (i = 0; i < n; ++i) {
                    y[iy] += a*x[ix];
                    ix += incx;
                    iy += incy;
                }
            }
        }
    }
}
#endif /* OPK_AXPY */

/*--------------------------------------------------------------------------------------*/
/* DOT */

#ifdef OPK_DOT
real_t
OPK_DOT(opk_index n,
        const real_t x[], opk_index incx,
        const real_t y[], opk_index incy)
{
    real_t t0, t1, t2, t3;
    opk_index m, i, ix, iy;

    t0 = OPK_REALCONST(0.0);
    if (n <= 0) {
        return t0;
    }
    if (incx == 1) {
        if (incy == 1) {
            t1 = t2 = t3 = t0;
#define UNROLL_OP(t,i) t += x[i]*y[i]
            UNROLL_4b(i,m,n);
#undef UNROLL_OP
            return (t0 + t1) + (t2 + t3);
        } else {
            iy = START_0(n, incy);
            for (ix = 0; ix < n; ++ix, iy += incy) {
                t0 += x[ix]*y[iy];
            }
        }
    } else {
        if (incy == 1) {
            ix = START_0(n, incx);
            for (iy = 0; iy < n; ++iy, ix += incx) {
                t0 += x[ix]*y[iy];
            }
        } else {
            ix = START_0(n, incx);
            iy = START_0(n, incy);
            for (i = 0; i < n; ++i, ix += incx, iy += incy) {
                t0 += x[ix]*y[iy];
            }
        }
    }

    return t0;
}
#endif /* OPK_DOT */

/*--------------------------------------------------------------------------------------*/
/* ASUM */

#ifdef OPK_ASUM
real_t
OPK_ASUM(opk_index n, const real_t x[],
         opk_index incx)
{
    real_t t0, t1, t2, t3;
    opk_index i, m;

    t0 = OPK_REALCONST(0.0);
    if (n < 1 || incx < 1) {
        return t0;
    }
    if (incx == 1) {
        t1 = t2 = t3 = t0;
# define UNROLL_OP(t,i) t += FABS(x[i])
        UNROLL_4b(i,m,n);
# undef UNROLL_OP
        return (t0 + t1) + (t2 + t3);
    }
    m = n*incx;
    for (i = 0; i < m; i += incx) {
        t0 += FABS(x[i]);
    }
    return t0;
}
#endif /* OPK_ASUM */

/*--------------------------------------------------------------------------------------*/
/* SUM */

#ifdef OPK_SUM
real_t
OPK_SUM(opk_index n, const real_t x[],
        opk_index incx)
{
    real_t t0, t1, t2, t3;
    opk_index i, m;

    t0 = OPK_REALCONST(0.0);
    if (n < 1 || incx < 1) {
        return t0;
    }
    if (incx == 1) {
        t1 = t2 = t3 = t0;
# define UNROLL_OP(t,i) t += x[i]
        UNROLL_4b(i,m,n);
# undef UNROLL_OP
        return (t0 + t1) + (t2 + t3);
    }
    m = n*incx;
    for (i = 0; i < m; i += incx) {
        t0 += x[i];
    }
    return t0;
}
#endif /* OPK_SUM */

/*--------------------------------------------------------------------------------------*/
/* NRM2 (Euclidean norm) */

#ifdef OPK_NRM2
real_t
OPK_NRM2(opk_index n, const real_t x[],
         opk_index incx)
{
    const real_t ZERO = OPK_REALCONST(0.0);
    const real_t ONE = OPK_REALCONST(1.0);
    real_t a, b, s, t;
    opk_index i, m;

    if (n < 1 || incx < 1) {
        return ZERO;
    }

    /* Get maximum absolute value of X. */
    a = b = x[0];
    if (incx == 1) {
        for (i = 1; i < n; ++i) {
            if (x[i] < a) a = x[i];
            if (x[i] > b) b = x[i];
        }
    } else {
        m = n*incx;
        for (i = incx; i < m; i += incx) {
            if (x[i] < a) a = x[i];
            if (x[i] > b) b = x[i];
        }
    }
    if ((a = -a) < b) {
        a = b;
    }
    if (a <= ZERO) {
        return ZERO;
    }

    /* Compute sum of squares avoiding overflows. */
    s = ZERO;
    b = ONE/a;
    if (incx == 1) {
        for (i = 0; i < n; ++i) {
            t = b*x[i];
            s += t*t;
        }
    } else {
        m = n*incx;
        for (i = 0; i < m; i += incx) {
            t = b*x[i];
            s += t*t;
        }
    }
    return a*SQRT(s);
}
#endif /* OPK_NRM2 */

/*--------------------------------------------------------------------------------------*/
/* IAMAX */

#ifdef OPK_IAMAX
opk_index
OPK_IAMAX(opk_index n, const real_t x[],
          opk_index incx)
{
    real_t t, amax;
    opk_index i, ix, imax;

    if (n < 1L || incx < 1L) return 0L;
    if (n == 1L) return 1L;
    amax = FABS(x[0]);
    imax = 0L;
    if (incx == 1L) {
        for (i = 1L; i < n; ++i) {
            t = FABS(x[i]);
            if (t > amax) {
                amax = t;
                imax = i;
            }
        }
    } else /* incx > 1 */ {
        ix = incx;
        for (i = 1L; i < n; ++i) {
            t = FABS(x[ix]);
            if (t > amax) {
                amax = t;
                imax = i;
            }
            ix += incx;
        }
    }
    return imax + 1L; /* FORTRAN indexing */
}
#endif /* OPK_IAMAX */

#if 0 /* code below is intended to be smarter (and faster) but still need to be
         tested... */
#ifdef OPK_IAMAX
opk_index
OPK_IAMAX(opk_index n,
          const real_t x[],
          opk_index incx)
{
    real_t a, b;
    opk_index i, m, ia, ib;

    if (n < 1L || incx < 1L) return 0L;
    if (n == 1L) return 1L;
    if (incx == 1L) {
        ia = ib = -1L;
        a = b = OPK_REALCONST(0.0);
        for (i = 0L; i < n; ++i) {
            if (x[i] < a) {
                a = x[i];
                ia = i;
            }
            if (x[i] > b) {
                b = x[i];
                ib = i;
            }
        }
    } else /* incx > 1 */ {
        ia = ib = -incx;
        a = b = OPK_REALCONST(0.0);
        m = n*incx;
        for (i = 0L; i < m; i += incx) {
            if (x[i] < a) {
                a = x[i];
                ia = i;
            }
            if (x[i] > b) {
                b = x[i];
                ib = i;
            }
        }
    }
    /* FORTRAN indexing */
    return ((a = -a) > b ? ia : (a < b ? ib : (ia <= ib ? ia : ib))) + 1L;
}
#endif /* OPK_IAMAX */
#endif

/*--------------------------------------------------------------------------------------*/

/* For IEEE floating-points, it is sufficient to fill the memory with null bytes. This may
   however depends on data representation. Hence we cannot use `memset` and we really copy
   floating point zeroes with some loop unrolling to achieve some speedup. */

#ifdef OPK_ZERO
void
OPK_ZERO(opk_index n, real_t x[], opk_index incx)
{
    const real_t ZERO = OPK_REALCONST(0.0);
    opk_index i, m;

    if (n >= 1) {
        if (incx == 1) {
# define UNROLL_OP(i) x[i] = ZERO
            UNROLL_8a(i,m,n);
# undef UNROLL_OP
        } else if (incx > 1) {
            m = n*incx;
            for (i = 0; i < m; i += incx) {
                x[i] = ZERO;
            }
        }
    }
}
#endif /* OPK_ZERO */

/*---[ Level 2 BLAS-like routines ]-----------------------------------------------------*/

#ifdef OPK_GEMV
/*
 *  DGEMV: Matrix-vector multiplication.
 *
 *  Level 2 Blas routine.
 *
 *  -- Written on 22-October-1986.
 *     Jack Dongarra, Argonne National Lab.
 *     Jeremy Du Croz, Nag Central Office.
 *     Sven Hammarling, Nag Central Office.
 *     Richard Hanson, Sandia National Labs.
 *
 *  -- C-version on 25 January 2006 by Eric Thiébaut (CRAL);
 *     non-zero integer result, says K, means invalid K-th argument.
 */
int
OPK_GEMV(opk_blas_trans trans, opk_index m, opk_index n,
         real_t alpha, const real_t a_[], opk_index lda,
         const real_t x_[], opk_index incx, real_t beta,
         real_t y_[], opk_index incy)
{
    /* Define some macros to mimic FORTRAN indexing. */
#undef y
#define y(a1) y_[a1]
#undef x
#define x(a1) x_[a1]
#undef a
#define a(a1,a2) a_[(a2)*lda + (a1)]

    /* Local constants and variables. */
    const real_t ZERO = OPK_REALCONST(0.0);
    const real_t ONE  = OPK_REALCONST(1.0);
    real_t temp;
    opk_index i, ix, iy, j, jx, jy, kx, ky, lenx, leny;
    int notrans;

    /* Check arguments and set LENX and LENY, the lengths
       of the vectors x and y. */
    notrans = (trans == OPK_BLAS_NO_TRANS);
    if (notrans) {
        lenx = n;
        leny = m;
    } else if (trans == OPK_BLAS_TRANS || trans == OPK_BLAS_CONJ_TRANS) {
        lenx = m;
        leny = n;
    } else {
        return 1;
    }
    if (m < 0) {
        return 2;
    }
    if (n < 0) {
        return 3;
    }
    if (lda < m || lda < 1) {
        return 6;
    }
    if (incx == 0) {
        return 8;
    }
    if (incy == 0) {
        return 11;
    }

    /* Quick return if possible. */
    if ((m == 0) || (n == 0) || ((alpha == ZERO) && (beta == ONE))) {
        return 0;
    }

    /* Set up the start points in X and Y. */
    kx = START_0(lenx, incx);
    ky = START_0(leny, incy);

    /*
     * Start the operations. In this version the elements of A are accessed sequentially
     * with ONE pass through A along the storage order of the elements of A.
     */

    /*
     * First form  y := beta*y.
     */
    if (beta != ONE) {
        if (incy == 1) {
            if (beta == ZERO) {
                LOOP(i, leny) {
                    y(i) = ZERO;
                }
            } else {
                LOOP(i, leny) {
                    y(i) *= beta;
                }
            }
        } else {
            iy = ky;
            if (beta == ZERO) {
                LOOP(i, leny) {
                    y(iy) = ZERO;
                    iy += incy;
                }
            } else {
                LOOP(i, leny) {
                    y(iy) *= beta;
                    iy += incy;
                }
            }
        }
    }

    if (alpha != ZERO) {
        if (notrans) {
            /*
             * Form  y := alpha*A*x + y.
             */
            jx = kx;
            if (incy == 1) {
                LOOP(j, n) {
                    if (x(jx) != ZERO) {
                        temp = alpha*x(jx);
                        LOOP(i, m) {
                            y(i) += temp*a(i, j);
                        }
                    }
                    jx += incx;
                }
            } else {
                LOOP(j, n) {
                    if (x(jx) != ZERO) {
                        temp = alpha*x(jx);
                        iy = ky;
                        LOOP(i, m) {
                            y(iy) += temp*a(i, j);
                            iy += incy;
                        }
                    }
                    jx += incx;
                }
            }
        } else {
            /*
             * Form  y := alpha*A'*x + y.
             */
            jy = ky;
            if (incx == 1) {
                LOOP(j, n) {
                    temp = ZERO;
                    LOOP(i, m) {
                        temp += a(i, j)*x(i);
                    }
                    y(jy) += alpha*temp;
                    jy += incy;
                }
            } else {
                LOOP(j, n) {
                    temp = ZERO;
                    ix = kx;
                    LOOP(i, m) {
                        temp += a(i, j)*x(ix);
                        ix += incx;
                    }
                    y(jy) += alpha*temp;
                    jy += incy;
                }
            }
        }
    }
    return 0;
}
#endif /* OPK_GEMV */

#ifdef OPK_TRMV
/*
 *  DTRMV: Multiplication of a vector by a triangular matrix.
 *
 *  Level 2 Blas routine.
 *
 *  -- Written on 22-October-1986.
 *     Jack Dongarra, Argonne National Lab.
 *     Jeremy Du Croz, Nag Central Office.
 *     Sven Hammarling, Nag Central Office.
 *     Richard Hanson, Sandia National Labs.
 *
 *  -- C-version on 25 January 2006 by Eric Thiébaut (CRAL);
 *     non-zero integer result, says K, means invalid K-th argument.
 */
int
OPK_TRMV(opk_blas_uplo uplo,
         opk_blas_trans trans,
         opk_blas_diag diag,
         opk_index n, const real_t a_[], opk_index lda,
         real_t x_[], opk_index incx)
{
    /* Define some macros to mimic FORTRAN indexing. */
#undef x
#define x(a1) x_[a1]
#undef a
#define a(a1,a2) a_[(a2)*lda + a1]

    /* Local constants and variables. */
    const real_t ZERO  = OPK_REALCONST(0.0);
    real_t temp;
    opk_index i, ix, j, jx, kx;
    int nounit, notrans, upper;

    /* Check arguments. */
    if (uplo == OPK_BLAS_UPPER) {
        upper = 1;
    } else if (uplo == OPK_BLAS_LOWER) {
        upper = 0;
    } else {
        return 1;
    }
    if (trans == OPK_BLAS_NO_TRANS) {
        notrans = 1;
    } else if (trans == OPK_BLAS_TRANS || trans == OPK_BLAS_CONJ_TRANS) {
        notrans = 0;
    } else {
        return 2;
    }
    if (diag == OPK_BLAS_NON_UNIT) {
        nounit = 1;
    } else if (diag == OPK_BLAS_UNIT) {
        nounit = 0;
    } else {
        return 3;
    }
    if (n < 0) {
        return 4;
    }
    if (lda < n || lda < 1) {
        return 6;
    }
    if (incx == 0) {
        return 8;
    }

    /* Quick return if possible. */
    if (n == 0) {
        return 0;
    }

    /* Set up the start point in X if the increment is not unity. This will be
       (N - 1)*INCX too small for descending loops. */
    kx = START_0(n, incx);

    /* Start the operations. In this version the elements of A are accessed sequentially
       with one pass through A. */
    if (notrans) {
        /* Form  x := A*x. */
        if (upper) {
            if (incx == 1) {
                for (j = 0; j < n; ++j) {
                    if (x(j) != ZERO) {
                        temp = x(j);
                        for (i = 0; i < j; ++i) {
                            x(i) += temp*a(i,j);
                        }
                        if (nounit) x(j) *= a(j,j);
                    }
                }
            } else {
                jx = kx;
                for (j = 0; j < n; ++j) {
                    if (x(jx) != ZERO) {
                        temp = x(jx);
                        ix = kx;
                        for (i = 0; i < j; ++i) {
                            x(ix) += temp*a(i,j);
                            ix += incx;
                        }
                        if (nounit) x(jx) *= a(j,j);
                    }
                    jx += incx;
                }
            }
        } else {
            if (incx == 1) {
                for (j = n - 1; j >= 0; --j) {
                    if (x(j) != ZERO) {
                        temp = x(j);
                        for (i = n - 1; i > j; --i) {
                            x(i) += temp*a(i,j);
                        }
                        if (nounit) x(j) *= a(j,j);
                    }
                }
            } else {
                kx += (n - 1)*incx;
                jx = kx;
                for (j = n - 1; j >= 0; --j) {
                    if (x(jx) != ZERO) {
                        temp = x(jx);
                        ix = kx;
                        for (i = n - 1; i > j; --i) {
                            x(ix) += temp*a(i,j);
                            ix -= incx;
                        }
                        if (nounit) x(jx) *= a(j,j);
                    }
                    jx -= incx;
                }
            }
        }
    } else {
        /*
         * Form  x := A'*x.
         */
        if (upper) {
            if (incx == 1) {
                for (j = n - 1; j >= 0; --j) {
                    temp = x(j);
                    if (nounit) temp *= a(j,j);
                    for (i = j - 1; i >= 0; --i) {
                        temp += a(i,j)*x(i);
                    }
                    x(j) = temp;
                }
            } else {
                jx = kx + (n - 1)*incx;
                for (j = n - 1; j >= 0; --j) {
                    temp = x(jx);
                    if (nounit) temp *= a(j,j);
                    ix = jx;
                    for (i = j - 1; i >= 0; --i) {
                        ix -= incx;
                        temp += a(i,j)*x(ix);
                    }
                    x(jx) = temp;
                    jx -= incx;
                }
            }
        } else {
            if (incx == 1) {
                for (j = 0; j < n; ++j) {
                    temp = x(j);
                    if (nounit) temp *= a(j,j);
                    for (i = j + 1; i < n; ++i) {
                        temp += a(i,j)*x(i);
                    }
                    x(j) = temp;
                }
            } else {
                jx = kx;
                for (j = 0; j < n; ++j) {
                    temp = x(jx);
                    if (nounit) temp *= a(j,j);
                    ix = jx;
                    for (i = j + 1; i < n; ++i) {
                        ix += incx;
                        temp += a(i,j)*x(ix);
                    }
                    x(jx) = temp;
                    jx += incx;
                }
            }
        }
    }

    return 0;
}
#endif /* OPK_TRMV */

#ifdef OPK_TRSV
/*
 * DTRSV: Solve a tringular linear system of equations.
 *
 * Level 2 Blas routine.
 *
 * -- Written on 22-October-1986.
 *    Jack Dongarra, Argonne National Lab.
 *    Jeremy Du Croz, Nag Central Office.
 *    Sven Hammarling, Nag Central Office.
 *    Richard Hanson, Sandia National Labs.
 *
 * -- C-version on 25 January 2006 by Eric Thiébaut (CRAL);
 *    non-zero integer result K, means invalid K-th argument.
 */
int
OPK_TRSV(opk_blas_uplo uplo,
         opk_blas_trans trans,
         opk_blas_diag diag,
         opk_index n, const real_t a_[], opk_index lda,
         real_t x_[], opk_index incx)
{
    /* Define some macros to mimic FORTRAN indexing. */
#undef x
#define x(a1) x_[a1-1]
#undef a
#define a(a1,a2) a_[a1-1+lda*(a2-1)]

    /* Local constants and variables. */
    const real_t ZERO = OPK_REALCONST(0.0);
    real_t temp;
    opk_index i, ix, j, jx, kx;
    int nounit, notrans, upper;

    /* Check arguments. */
    if (uplo == OPK_BLAS_UPPER) {
        upper = 1;
    } else if (uplo == OPK_BLAS_LOWER) {
        upper = 0;
    } else {
        return 1;
    }
    if (trans == OPK_BLAS_NO_TRANS) {
        notrans = 1;
    } else if (trans == OPK_BLAS_TRANS || trans == OPK_BLAS_CONJ_TRANS) {
        notrans = 0;
    } else {
        return 2;
    }
    if (diag == OPK_BLAS_NON_UNIT) {
        nounit = 1;
    } else if (diag == OPK_BLAS_UNIT) {
        nounit = 0;
    } else {
        return 3;
    }
    if (n < 0) {
        return 4;
    }
    if (lda < MAX(1, n)) {
        return 6;
    }
    if (incx == 0) {
        return 8;
    }

    /* Quick return if possible. */
    if (n == 0) {
        return 0;
    }

    /* Set up the start point in X if the increment is not unity. This will be
       (N - 1)*INCX too small for descending loops. */
    kx = (incx >= 0 ? 1 : 1 + (1 - n) * incx);

    /* Start the operations. In this version the elements of A are accessed sequentially
       with one pass through A. */
    if (notrans) {
        /* Form  x := inv(A)*x. */
        if (upper) {
            if (incx == 1) {
                for (j = n; j >= 1; --j) {
                    if (x(j) != ZERO) {
                        if (nounit) x(j) /= a(j,j);
                        temp = x(j);
                        for (i = j - 1; i >= 1; --i) {
                            x(i) -= temp*a(i,j);
                        }
                    }
                }
            } else {
                jx = kx + (n - 1)*incx;
                for (j = n; j >= 1; --j) {
                    if (x(jx) != ZERO) {
                        if (nounit) x(jx) /= a(j,j);
                        temp = x(jx);
                        ix = jx;
                        for (i = j - 1; i >= 1; --i) {
                            ix -= incx;
                            x(ix) -= temp*a(i,j);
                        }
                    }
                    jx -= incx;
                }
            }
        } else {
            if (incx == 1) {
                for (j = 1; j <= n; ++j) {
                    if (x(j) != ZERO) {
                        if (nounit) x(j) /= a(j,j);
                        temp = x(j);
                        for (i = j + 1; i <= n; ++i) {
                            x(i) -= temp*a(i,j);
                        }
                    }
                }
            } else {
                jx = kx;
                for (j = 1; j <= n; ++j) {
                    if (x(jx) != ZERO) {
                        if (nounit) x(jx) /= a(j,j);
                        temp = x(jx);
                        ix = jx;
                        for (i = j + 1; i <= n; ++i) {
                            ix += incx;
                            x(ix) -= temp*a(i,j);
                        }
                    }
                    jx += incx;
                }
            }
        }
    } else {
        /* Form  x := inv(A')*x. */
        if (upper) {
            if (incx == 1) {
                for (j = 1; j <= n; ++j) {
                    temp = x(j);
                    for (i = 1; i < j; ++i) {
                        temp -= a(i,j)*x(i);
                    }
                    x(j) = (nounit ? temp/a(j,j) : temp);
                }
            } else {
                jx = kx;
                for (j = 1; j <= n; ++j) {
                    temp = x(jx);
                    ix = kx;
                    for (i = 1; i < j; ++i) {
                        temp -= a(i,j)*x(ix);
                        ix += incx;
                    }
                    x(jx) = (nounit ? temp/a(j,j) : temp);
                    jx += incx;
                }
            }
        } else {
            if (incx == 1) {
                for (j = n; j >= 1; --j) {
                    temp = x(j);
                    for (i = n; i > j; --i) {
                        temp -= a(i,j)*x(i);
                    }
                    x(j) = (nounit ? temp/a(j,j) : temp);
                }
            } else {
                kx += (n - 1)*incx;
                jx = kx;
                for (j = n; j >= 1; --j) {
                    temp = x(jx);
                    ix = kx;
                    for (i = n; i > j; --i) {
                        temp -= a(i,j)*x(ix);
                        ix -= incx;
                    }
                    x(jx) = (nounit ? temp/a(j,j) : temp);
                    jx -= incx;
                }
            }
        }
    }

    return 0;
}
#endif /* OPK_TRSV */

/*---[ Level 3 BLAS-like routines ]------------------------------------------*/

#ifdef OPK_GEMM
/*
 * DGEMM: Performs a matrix-matrix operation.
 *
 * Level 3 Blas routine.
 *
 * -- Written on 8-February-1989.
 *    Jack Dongarra, Argonne National Laboratory.
 *    Iain Duff, AERE Harwell.
 *    Jeremy Du Croz, Numerical Algorithms Group Ltd.
 *    Sven Hammarling, Numerical Algorithms Group Ltd.
 *
 * -- C-version on 25 January 2006 by Eric Thiébaut (CRAL).
 *    non-zero integer result K, means invalid K-th argument.
 */
int
OPK_GEMM(opk_blas_trans transa,
         opk_blas_trans transb,
         opk_index m, opk_index n, opk_index k,
         real_t alpha, const real_t a_[], opk_index lda,
         const real_t b_[], opk_index ldb, real_t beta,
         real_t c_[], opk_index ldc)
{
    /* Define some macros to mimic FORTRAN indexing. */
#undef c
#define c(a1,a2) c_[a1 + ldc*(a2)]
#undef b
#define b(a1,a2) b_[a1 + ldb*(a2)]
#undef a
#define a(a1,a2) a_[a1 + lda*(a2)]

    /* Local constants and variables. */
    const real_t ZERO = OPK_REALCONST(0.0);
    const real_t ONE = OPK_REALCONST(1.0);
    real_t temp;
    opk_index i, j, l, /* ncola,*/ nrowa, nrowb;
    int nota, notb;

    /* Test the input parameters. Set NOTA and NOTB as true if A and B respectively are
       not transposed and set NROWA, NCOLA and NROWB as the number of rows and columns of
       A and the number of rows of B respectively. */
    if (transa == OPK_BLAS_NO_TRANS) {
        nota = 1;
        nrowa = m;
        /*ncola = k;*/
    } else if (transa == OPK_BLAS_TRANS || transa == OPK_BLAS_CONJ_TRANS) {
        nota = 0;
        nrowa = k;
        /*ncola = m;*/
    } else {
        return 1;
    }
    if (transb == OPK_BLAS_NO_TRANS) {
        notb = 1;
        nrowb = k;
    } else if (transb == OPK_BLAS_TRANS || transb == OPK_BLAS_CONJ_TRANS) {
        notb = 0;
        nrowb = n;
    } else {
        return 1;
    }
    if (m < 0) return 3;
    if (n < 0) return 4;
    if (k < 0) return 5;
    if (lda < 1 || lda < nrowa) return 8;
    if (ldb < 1 || ldb < nrowb) return 10;
    if (ldc < 1 || ldc < m) return 13;

    /* Quick return if possible. */
    if (m == 0 || n == 0 || (beta == ONE && (alpha == ZERO || k == 0))) {
        return 0;
    }

    /* And if alpha = ZERO. */
    if (alpha == ZERO) {
        if (beta == ZERO) {
            for (j = 0; j < n; ++j) {
                for (i = 0; i < m; ++i) {
                    c(i,j) = ZERO;
                }
            }
        } else {
            for (j = 0; j < n; ++j) {
                for (i = 0; i < m; ++i) {
                    c(i,j) *= beta;
                }
            }
        }
        return 0;
    }

    /* Start the operations. */
    if (notb) {
        if (nota) {
            /* Form C := alpha*A*B + beta*C. */
            for (j = 0; j < n; ++j) {
                if (beta == ZERO) {
                    for (i = 0; i < m; ++i) {
                        c(i,j) = ZERO;
                    }
                } else if (beta != ONE) {
                    for (i = 0; i < m; ++i) {
                        c(i,j) *= beta;
                    }
                }
                for (l = 0; l < k; ++l) {
                    if (b(l,j) != ZERO) {
                        temp = alpha*b(l,j);
                        for (i = 0; i < m; ++i) {
                            c(i,j) += temp*a(i,l);
                        }
                    }
                }
            }
        } else {
            /* Form C := alpha*A'*B + beta*C */
            for (j = 0; j < n; ++j) {
                for (i = 0; i < m; ++i) {
                    temp = ZERO;
                    for (l = 0; l < k; ++l) {
                        temp += a(l,i)*b(l,j);
                    }
                    if (beta == ZERO) {
                        c(i,j) = alpha*temp;
                    } else {
                        c(i,j) = alpha*temp + beta*c(i,j);
                    }
                }
            }
        }
    } else {
        if (nota) {
            /* Form C := alpha*A*B' + beta*C */
            for (j = 0; j < n; ++j) {
                if (beta == ZERO) {
                    for (i = 0; i < m; ++i) {
                        c(i,j) = ZERO;
                    }
                } else if (beta != ONE) {
                    for (i = 0; i < m; ++i) {
                        c(i,j) *= beta;
                    }
                }
                for (l = 0; l < k; ++l) {
                    if (b(j,l) != ZERO) {
                        temp = alpha*b(j, l);
                        for (i = 0; i < m; ++i) {
                            c(i,j) += temp*a(i,l);
                        }
                    }
                }
            }
        } else {
            /* Form C := alpha*A'*B' + beta*C */
            for (j = 0; j < n; ++j) {
                for (i = 0; i < m; ++i) {
                    temp = ZERO;
                    for (l = 0; l < k; ++l) {
                        temp += a(l,i)*b(j,l);
                    }
                    if (beta == ZERO) {
                        c(i,j) = alpha*temp;
                    } else {
                        c(i,j) = alpha*temp + beta*c(i,j);
                    }
                }
            }
        }
    }

    return 0;
}
#endif /* OPK_GEMM */

#ifdef OPK_SYRK
/*
 * DSYRK performs one of the symmetric rank k operations
 *
 *   C := alpha*A*A' + beta*C,
 *
 * or
 *
 *   C := alpha*A'*A + beta*C,
 *
 * where alpha and beta are scalars, C is an n by n symmetric matrix and A is an n by k
 * matrix in the first case and a k by n matrix in the second case.
 *
 * Level 3 Blas routine.
 *
 * -- Written on 8-February-1989.
 *    Jack Dongarra, Argonne National Laboratory.
 *    Iain Duff, AERE Harwell.
 *    Jeremy Du Croz, Numerical Algorithms Group Ltd.
 *    Sven Hammarling, Numerical Algorithms Group Ltd.
 *
 * -- C-version on 25 January 2006 by Eric Thiébaut (CRAL).
 *    non-zero integer result K, means invalid K-th argument.
 */
int
OPK_SYRK(opk_blas_uplo uplo,
         opk_blas_trans trans,
         opk_index n, opk_index k, real_t alpha,
         const real_t a_[], opk_index lda,
         real_t beta, real_t c_[], opk_index ldc)
{
    /* Define some macros to mimic FORTRAN indexing. */
#undef c
#define c(a1,a2) c_[(a2)*ldc + a1]
#undef a
#define a(a1,a2) a_[(a2)*ldc + a1]

    /* Local constants and variables. */
    const real_t ZERO = OPK_REALCONST(0.0);
    const real_t ONE = OPK_REALCONST(1.0);
    real_t temp;
    opk_index i, j, l, nrowa;
    int upper, notrans;

    /* Test the input parameters. */
    if (uplo == OPK_BLAS_UPPER) {
        upper = 1;
    } else if (uplo == OPK_BLAS_LOWER) {
        upper = 0;
    } else {
        return 1;
    }
    if (trans == OPK_BLAS_NO_TRANS) {
        notrans = 1;
        nrowa = n;
    } else if (trans == OPK_BLAS_TRANS || trans == OPK_BLAS_CONJ_TRANS) {
        notrans = 0;
        nrowa = k;
    } else {
        return 2;
    }
    if (n < 0) return 3;
    if (k < 0) return 4;
    if (lda < 1 || lda < nrowa) return 7;
    if (ldc < 1 || ldc < n) return 10;

    /* Quick return if possible. */
    if (n == 0 || (beta == ONE && (alpha == ZERO || k == 0))) {
        return 0;
    }

    /* And when alpha = ZERO. */
    if (alpha == ZERO) {
        if (upper) {
            if (beta == ZERO) {
                for (j = 0; j < n; ++j) {
                    for (i = 0; i <= j; ++i) {
                        c(i,j) = ZERO;
                    }
                }
            } else {
                for (j = 0; j < n; ++j) {
                    for (i = 0; i <= j; ++i) {
                        c(i,j) *= beta;
                    }
                }
            }
        } else {
            if (beta == ZERO) {
                for (j = 0; j < n; ++j) {
                    for (i = j; i < n; ++i) {
                        c(i,j) = ZERO;
                    }
                }
            } else {
                for (j = 0; j < n; ++j) {
                    for (i = j; i < n; ++i) {
                        c(i,j) *= beta;
                    }
                }
            }
        }
        return 0;
    }

    /* Start the operations. */
    if (notrans) {
        /* Form  C := alpha*A*A' + beta*C. */
        if (upper) {
            for (j = 0; j < n; ++j) {
                if (beta == ZERO) {
                    for (i = 0; i <= j; ++i) {
                        c(i,j) = ZERO;
                    }
                } else if (beta != ONE) {
                    for (i = 0; i <= j; ++i) {
                        c(i,j) *= beta;
                    }
                }
                for (l = 0; l < k; ++l) {
                    if (a(j,l) != ZERO) {
                        temp = alpha*a(j,l);
                        for (i = 0; i <= j; ++i) {
                            c(i,j) += temp*a(i,l);
                        }
                    }
                }
            }
        } else {
            for (j = 0; j < n; ++j) {
                if (beta == ZERO) {
                    for (i = j; i < n; ++i) {
                        c(i,j) = ZERO;
                    }
                } else if (beta != ONE) {
                    for (i = j; i < n; ++i) {
                        c(i,j) *= beta;
                    }
                }
                for (l = 0; l < k; ++l) {
                    if (a(j,l) != ZERO) {
                        temp = alpha*a(j,l);
                        for (i = j; i < n; ++i) {
                            c(i,j) += temp*a(i,l);
                        }
                    }
                }
            }
        }
    } else {
        /* Form  C := alpha*A'*A + beta*C. */
        if (upper) {
            for (j = 0; j < n; ++j) {
                for (i = 0; i <= j; ++i) {
                    temp = ZERO;
                    for (l = 0; l < k; ++l) {
                        temp += a(l,i)*a(l,j);
                    }
                    if (beta == ZERO) {
                        c(i,j) = alpha*temp;
                    } else {
                        c(i,j) = alpha*temp + beta*c(i,j);
                    }
                }
            }
        } else {
            for (j = 0; j < n; ++j) {
                for (i = j; i < n; ++i) {
                    temp = ZERO;
                    for (l = 0; l < k; ++l) {
                        temp += a(l,i)*a(l,j);
                    }
                    if (beta == ZERO) {
                        c(i,j) = alpha*temp;
                    } else {
                        c(i,j) = alpha*temp + beta*c(i,j);
                    }
                }
            }
        }
    }

    return 0;
}
#endif /* OPK_SYRK */

/*---[ LAPACK-like routines ]------------------------------------------------*/

#ifdef OPK_POTF2
/*
 * DPOTF2: Computes the Cholesky factorization of a real symmetric positive definite
 *         matrix A.
 *
 * -- LAPACK routine (version 3.0) --
 *    Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
 *    Courant Institute, Argonne National Lab, and Rice University
 *    February 29, 1992
 *
 * -- C-version on 25 January 2006 by Eric Thiébaut (CRAL).
 */
opk_index
OPK_POTF2(opk_blas_uplo uplo,
          opk_index n,
          real_t a_[],
          opk_index lda)
{
    /* Define some macros to mimic FORTRAN indexing. */
    /* FIXME:*/
#undef a
#define a(a1,a2) a_[a1-1+lda*(a2-1)]
#undef one

    /* Local constants and variables. */
    const real_t ZERO = OPK_REALCONST(0.0);
    const real_t ONE = OPK_REALCONST(1.0);
    real_t ajj;
    opk_index j;

    /* Test the input parameters. */
    if (n < 0) {
        return -2;
    }
    if (lda < 1 || lda < n) {
        return -4;
    }

    /* Quick return if possible. */
    if (n == 0) {
        return 0;
    }

    if (uplo == OPK_BLAS_UPPER) {

        /* Compute the Cholesky factorization A = U'*U. */

        for (j = 1; j <= n; ++j) {

            /* Compute U(J,J) and test for non-positive-definiteness. */
            ajj = a(j, j) - OPK_DOT(j - 1, &a(1, j), 1, &a(1, j), 1);
            if (ajj <= ZERO) {
                a(j, j) = ajj;
                return j;
            }
            ajj = SQRT(ajj);
            a(j, j) = ajj;

            /* Compute elements J+1:N of row J. */
            if (j < n) {
                OPK_GEMV(OPK_BLAS_TRANS, j - 1, n - j, -ONE, &a(1, j + 1),
                         lda, &a(1, j), 1, ONE, &a(j, j + 1), lda);
                OPK_SCAL(n - j, ONE/ajj, &a(j, j + 1), lda);
            }

        }

    } else if (uplo == OPK_BLAS_LOWER) {

        /* Compute the Cholesky factorization A = L*L'. */

        for (j = 1; j <= n; ++j) {

            /* Compute L(J,J) and test for non-positive-definiteness. */
            ajj = a(j, j) - OPK_DOT(j - 1, &a(j, 1), lda, &a(j, 1), lda);
            if (ajj <= ZERO) {
                a(j, j) = ajj;
                return j;
            }
            ajj = SQRT(ajj);
            a(j, j) = ajj;

            /* Compute elements J+1:N of column J. */
            if (j < n) {
                OPK_GEMV(OPK_BLAS_NO_TRANS, n - j, j - 1, -ONE, &a(j + 1, 1),
                         lda, &a(j, 1), lda, ONE, &a(j + 1, j), 1);
                OPK_SCAL(n - j, ONE/ajj, &a(j + 1, j), 1);
            }

        }

    } else {

        /* Bad first parameter. */
        return -1;

    }
    return 0;
}
#endif /* OPK_POTF2 */

#endif /* OPK_LINALG_C_ */
