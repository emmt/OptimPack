/*
 * linalg-common.h --
 *
 * Private macro definitions for single/double precision version of numerical
 * routines.
 *
 *-----------------------------------------------------------------------------
 *
 * The OptimPack library is licensed under the MIT "Expat" License:
 *
 * Copyright (c) 2008-2014: Éric Thiébaut
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

#ifndef _OPK_LINALG_COMMON_H
#define _OPK_LINALG_COMMON_H 1

#ifndef _OPTIMPACK_H
# include "optimpack.h"
#endif
#ifndef _OPTIMPACK_LINALG_H
# include "optimpack-linalg.h"
#endif

/* Macros to replace some Fortran intrinsics. */
#define MIN(x,y)  ((x)<=(y)?(x):(y))
#define MAX(x,y)  ((x)>=(y)?(x):(y))
#define SIGN(x,y) ((((x)<0)==((y)<0))?(x):-(x))

#endif /* _OPK_LINALG_COMMON_H */

/* Undefined macros prior to proper definitions. */
#undef real_t
#undef OPK_REALCONST
#undef OPK_MATHFUNC
#undef OPK_BLAS
#undef OPK_LAPACK
#undef OPK_MINPACK2
#undef OPK_FUNC

/* Setup definitions to generate single precision version
   of the code. */
#if defined(OPK_SINGLE_PRECISION) && OPK_SINGLE_PRECISION
# define real_t float
# define OPK_BLAS(name)       OPK_JOIN(opk_s,name)
# define OPK_LAPACK(name)     OPK_JOIN(opk_s,name)
# define OPK_MINPACK2(name)   OPK_JOIN(opk_s,name)
# define OPK_FUNC(name)       OPK_JOIN(opk_,OPK_JOIN(name,_f))
# define OPK_REALCONST(expr)  OPK_JOIN(expr,F)
# define OPK_MATHFUNC(func)   OPK_JOIN(func,f)
#endif /* OPK_SINGLE_PRECISION */

/* Setup definitions to generate double precision version
   of the code. */
#if defined(OPK_DOUBLE_PRECISION) && OPK_DOUBLE_PRECISION
# define real_t double
# define OPK_BLAS(name)       OPK_JOIN(opk_d,name)
# define OPK_LAPACK(name)     OPK_JOIN(opk_d,name)
# define OPK_MINPACK2(name)   OPK_JOIN(opk_d,name)
# define OPK_FUNC(name)       OPK_JOIN(OPK_JOIN(opk_,name),_d)
# define OPK_REALCONST(expr)  expr
# define OPK_MATHFUNC(func)   func
#endif /* OPK_DOUBLE_PRECISION */

/* Apply rules defined above. */
#ifdef real_t

/* OPK_LCG: linear conjugate gradient */
# undef OPK_LCG
# define OPK_LCG     OPK_FUNC(lcg)

/* OPK_PLCG: preconditioned linear conjugate gradient */
# undef OPK_PLCG
# define OPK_PLCG    OPK_FUNC(plcg)

/* OPK_TRCG: trust region conjugate gradient */
# undef OPK_TRCG
# define OPK_TRCG    OPK_FUNC(trcg)

/* OPK_NLCG: non-linear conjugate gradient */
# undef OPK_NLCG
# define OPK_NLCG    OPK_FUNC(nlcg)

/* Level 1 BLAS-like functions. */
# undef OPK_AMAX
# define OPK_AMAX    OPK_BLAS(amax)
# undef OPK_ASUM
# define OPK_ASUM    OPK_BLAS(asum)
# undef OPK_AXPY
# define OPK_AXPY    OPK_BLAS(axpy)
# undef OPK_COPY
# define OPK_COPY    OPK_BLAS(copy)
# undef OPK_DOT
# define OPK_DOT     OPK_BLAS(dot)
# undef OPK_NRM2
# define OPK_NRM2    OPK_BLAS(nrm2)
# undef OPK_SCAL
# define OPK_SCAL    OPK_BLAS(scal)
# undef OPK_SUM
# define OPK_SUM     OPK_BLAS(sum)
# undef OPK_SWAP
# define OPK_SWAP    OPK_BLAS(swap)
# undef OPK_ZERO
# define OPK_ZERO    OPK_BLAS(zero)
# undef OPK_IAMAX
# if defined(OPK_SINGLE_PRECISION) && OPK_SINGLE_PRECISION
#  define OPK_IAMAX   opk_isamax
# else
#  define OPK_IAMAX   opk_idamax
# endif

/* Level 2 BLAS-like functions. */
# undef OPK_GEMV
# define OPK_GEMV    OPK_BLAS(gemv)
# undef OPK_TRMV
# define OPK_TRMV    OPK_BLAS(trmv)
# undef OPK_TRSV
# define OPK_TRSV    OPK_BLAS(trsv)

/* Level 3 BLAS-like functions. */
# undef OPK_GEMM
# define OPK_GEMM    OPK_BLAS(gemm)
# undef OPK_SYRK
# define OPK_SYRK    OPK_BLAS(syrk)

/* LAPACK-like functions. */
# undef OPK_POTF2
# define OPK_POTF2   OPK_LAPACK(potf2)

/* MINPACK-2 functions. */
# undef OPK_ESTSV
# define OPK_ESTSV   OPK_MINPACK2(estsv)
# undef OPK_GQT
# define OPK_GQT     OPK_MINPACK2(gqt)

/* Mathematical functions. */
# undef LOG
# define LOG       OPK_MATHFUNC(log)
# undef LOG10
# define LOG10     OPK_MATHFUNC(log10)
# undef EXP
# define EXP       OPK_MATHFUNC(exp)
# undef SQRT
# define SQRT      OPK_MATHFUNC(sqrt)
# undef CBRT
# define CBRT      OPK_MATHFUNC(cbrt)
# undef POW
# define POW       OPK_MATHFUNC(pow)
# undef FABS
# define FABS      OPK_MATHFUNC(fabs)
# undef HYPOT
# define HYPOT     OPK_MATHFUNC(hypot)
# undef COS
# define COS       OPK_MATHFUNC(cos)
# undef SIN
# define SIN       OPK_MATHFUNC(sin)
# undef TAN
# define TAN       OPK_MATHFUNC(tan)
# undef ACOS
# define ACOS      OPK_MATHFUNC(acos)
# undef ASIN
# define ASIN      OPK_MATHFUNC(asin)
# undef ATAN
# define ATAN      OPK_MATHFUNC(atan)
# undef ATAN2
# define ATAN2     OPK_MATHFUNC(atan2)
# undef COSH
# define COSH      OPK_MATHFUNC(cosh)
# undef SINH
# define SINH      OPK_MATHFUNC(sinh)
# undef TANH
# define TANH      OPK_MATHFUNC(tanh)
# undef ACOSH
# define ACOSH     OPK_MATHFUNC(acosh)
# undef ASINH
# define ASINH     OPK_MATHFUNC(asinh)
# undef ATANH
# define ATANH     OPK_MATHFUNC(atanh)
# undef FLOOR
# define FLOOR     OPK_MATHFUNC(floor)
# undef CEIL
# define CEIL      OPK_MATHFUNC(ceil)
# undef FMOD
# define FMOD      OPK_MATHFUNC(fmod)

#endif /* real_t */

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
