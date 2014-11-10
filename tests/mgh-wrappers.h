/*
 * mgh-wrappers.h --
 *
 * Definitions of wrappers to call FORTRAN subroutines of the unconstrained
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

#ifndef _MGH_WRAPPERS_H
#define _MGH_WRAPPERS_H 1

#ifdef __cplusplus
extern "C" {
#endif
#if 0 /* fool Emacs auto-indent */
}
#endif

extern void mgh_initpt(int n, double* x, int nprob, double factor);
extern double mgh_objfcn(int n, const double* x, int nprob);
extern void mgh_grdfcn(int n, const double* x, double* g, int nprob);
extern void mgh_hesfcn(int n, const double* x, double* h, int ldh, int nprob);
extern void mgh_expand(int n, const double* alin, double* amat, int lda);
extern void mgh_ssqfcn(int m, int n, const double* x, double* fvec, int nprob);
extern void mgh_ssqjac(int m, int n, const double* x, double* fjac, int ldfjac, int nprob);
extern void mgh_vecfcn(int n, const double* x, double* fvec, int nprob);
extern void mgh_vecjac(int n, const double* x, double* fjac, int ldfjac, int nprob);

#if 0 /* fool Emacs auto-indent */
{
#endif
#ifdef __cplusplus
}
#endif

#endif /* _MGH_WRAPPERS_H */

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
