/*
 * mgh-wrappers.c --
 *
 * Implementation of wrappers to call FORTRAN subroutines of the unconstrained
 * minimization test suite from the MINPACK-1 Project.
 *
 * References:
 * [1] J. J. Moré, B. S. Garbow and K. E. Hillstrom, "Testing unconstrained
 *     optimization software," ACM Trans. Math. Software 7, 17-41 (1981).
 * [2] J. J. Moré, B. S. Garbow and K. E. Hillstrom, "Fortran subroutines for
 *     testing unconstrained optimization software," ACM Trans. Math. Software
 *     7, 136-140 (1981).
 *
 * History:
 *  - Argonne National Laboratory. MINPACK-1 Project.  March 1980.
 *    Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More.
 *  - Conversion to Yorick.  November 2001. Éric Thiébaut.
 *  - Conversion to C.  february 2014. Éric Thiébaut.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2014 Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>.
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

#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "mgh-wrappers.h"

/* The following macros should be set to match the data types for the Fortran
   compiler used to build BLAS and/or LAPACK libraries. */
#ifndef CONST
# define CONST     const
#endif
#ifndef CHARACTER
# define CHARACTER char
#endif
#ifndef LOGICAL
# define LOGICAL int
#endif
#ifndef INTEGER
# define INTEGER int
#endif
#ifndef REAL
# define REAL double
#endif
#ifndef VOID
# define VOID      int /* type returned by FORTRAN subroutines */
#endif


#define JOIN2_(a,b)           a##b
#define JOIN3_(a,b,c)         a##b##c
#define JOIN4_(a,b,c,d)       a##b##c##d

#define JOIN(a,b)             JOIN2_(a,b)
#define JOIN2(a,b)            JOIN2_(a,b)
#define JOIN3(a,b,c)          JOIN3_(a,b,c)
#define JOIN4(a,b,c,d)        JOIN4_(a,b,c,d)

#define STRINGIFY_(x)         #x
#define STRINGIFY(x)          STRINGIFY_(x)

/*
 * Fortran compilers mangle names differently the name of a FORTRAN symbol may
 * be in upper or lower case letters and may be prefixed or post fixed by an
 * underscore letter.  To cope with this variety, we define the macro
 * FORTRAN_NAME to build the name of a FORTRAN symbol for the C compiler.
 * This macros take two arguments: the FORTRAN symbol name in lowercase and in
 * uppercase letters.  The value of the macro FORTRAN_STYLE may be set to
 * select the appropriate macro definition.  The value of FORTRAN_STYLE is an
 * integer value in the range 0:7 as follows:
 *
 *   - 1st bit set to use uppercase letters, otherwise lowercase;
 *   - 2nd bit set to use an underscore prefix;
 *   - 3rd bit set to use an underscore suffix.
 */
#ifndef FORTRAN_STYLE
# ifdef _WIN32
#  define FORTRAN_STYLE 1
# else
#  define FORTRAN_STYLE 2
# endif
#endif
#if (FORTRAN_STYLE == 0)
# define FORTRAN_NAME(name,NAME) name
#elif (FORTRAN_STYLE == 1)
# define FORTRAN_NAME(name,NAME) NAME
#elif (FORTRAN_STYLE == 2)
# define FORTRAN_NAME(name,NAME) name##_
#elif (FORTRAN_STYLE == 3)
# define FORTRAN_NAME(name,NAME) NAME##_
#elif (FORTRAN_STYLE == 4)
# define FORTRAN_NAME(name,NAME) _##name
#elif (FORTRAN_STYLE == 5)
# define FORTRAN_NAME(name,NAME) _##NAME
#elif (FORTRAN_STYLE == 6)
# define FORTRAN_NAME(name,NAME) _##name##_
#elif (FORTRAN_STYLE == 7)
# define FORTRAN_NAME(name,NAME) _##NAME##_
#else
# error illegal FORTRAN_STYLE
#endif

#define INITPT FORTRAN_NAME(initpt, INITPT)
#define OBJFCN FORTRAN_NAME(objfcn, OBJFCN)
#define GRDFCN FORTRAN_NAME(grdfcn, GRDFCN)
#define HESFCN FORTRAN_NAME(hesfcn, HESFCN)
#define EXPAND FORTRAN_NAME(expand, EXPAND)
#define SSQFCN FORTRAN_NAME(ssqfcn, SSQFCN)
#define SSQJAC FORTRAN_NAME(ssqjac, SSQJAC)
#define VECFCN FORTRAN_NAME(vecfcn, VECFCN)
#define VECJAC FORTRAN_NAME(vecjac, VECJAC)

extern VOID INITPT(INTEGER* n, REAL* x, INTEGER* nprob, REAL* factor);
extern VOID OBJFCN(INTEGER* n, REAL* x, REAL* f, INTEGER* nprob);
extern VOID GRDFCN(INTEGER* n, REAL* x, REAL* g, INTEGER* nprob);
extern VOID HESFCN(INTEGER* n, REAL* x, REAL* h, INTEGER* ldh, INTEGER* nprob);
extern VOID EXPAND(INTEGER* n, REAL* alin, REAL* amat, INTEGER* lda);
extern VOID SSQFCN(INTEGER* m, INTEGER* n, REAL* x, REAL* fvec, INTEGER* nprob);
extern VOID SSQJAC(INTEGER* m, INTEGER* n, REAL* x, REAL* fjac, INTEGER* ldfjac, INTEGER* nprob);
extern VOID VECFCN(INTEGER* n, REAL* x, REAL* fvec, INTEGER* nprob);
extern VOID VECJAC(INTEGER* n, REAL* x, REAL* fjac, INTEGER* ldfjac, INTEGER* nprob);

void mgh_initpt(int n, double* x, int nprob, double factor)
{
  REAL _factor = factor;
  INTEGER _n = n, _nprob = nprob;
  assert(sizeof(REAL) == sizeof(*x));
  INITPT(&_n, (REAL*)x, &_nprob, &_factor);
}

double mgh_objfcn(int n, const double* x, int nprob)
{
  REAL _f;
  INTEGER _n = n, _nprob = nprob;
  assert(sizeof(REAL) == sizeof(*x));
  OBJFCN(&_n, (REAL*)x, &_f, &_nprob);
  return (double)_f;
}

void mgh_grdfcn(int n, const double* x, double* g, int nprob)
{
  INTEGER _n = n, _nprob = nprob;
  assert(sizeof(REAL) == sizeof(*x));
  GRDFCN(&_n, (REAL*)x, (REAL*)g, &_nprob);
}

void mgh_hesfcn(int n, const double* x, double* h, int ldh, int nprob)
{
  INTEGER _n = n, _ldh = ldh, _nprob = nprob;
  assert(sizeof(REAL) == sizeof(*x));
  HESFCN(&_n, (REAL*)x, (REAL*)h, &_ldh, &_nprob);
}

void mgh_expand(int n, const double* alin, double* amat, int lda)
{
  INTEGER _n = n, _lda = lda;
  assert(sizeof(REAL) == sizeof(*alin));
  EXPAND(&_n, (REAL*)alin, (REAL*)amat, &_lda);
}

void mgh_ssqfcn(int m, int n, const double* x, double* fvec, int nprob)
{
  INTEGER _m = m, _n = n, _nprob = nprob;
  assert(sizeof(REAL) == sizeof(*x));
  SSQFCN(&_m, &_n, (REAL*)x, (REAL*)fvec, &_nprob);
}

void mgh_ssqjac(int m, int n, const double* x, double* fjac, int ldfjac, int nprob)
{
  INTEGER _m = m, _n = n, _ldfjac = ldfjac, _nprob = nprob;
  assert(sizeof(REAL) == sizeof(*x));
  SSQJAC(&_m, &_n, (REAL*)x, (REAL*)fjac, &_ldfjac, &_nprob);
}

void mgh_vecfcn(int n, const double* x, double* fvec, int nprob)
{
  INTEGER _n = n, _nprob = nprob;
  assert(sizeof(REAL) == sizeof(*x));
  VECFCN(&_n, (REAL*)x, (REAL*)fvec, &_nprob);
}

void mgh_vecjac(int n, const double* x, double* fjac, int ldfjac, int nprob)
{
  INTEGER _n = n, _ldfjac = ldfjac, _nprob = nprob;
  assert(sizeof(REAL) == sizeof(*x));
  VECJAC(&_n, (REAL*)x, (REAL*)fjac, &_ldfjac, &_nprob);
}

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: nil
 * c-basic-offset: 2
 * fill-column: 79
 * coding: utf-8
 * End:
 */
