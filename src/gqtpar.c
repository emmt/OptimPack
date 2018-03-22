/*
 * gqtpar.c --
 *
 * GQTPAR algorithm for "Computing A Trust Region Step" by Moré,
 * J.J. & Sorensen, D.C. (SIAM J. Sci. Stat. Comp., 1983, vol. 4, pp. 553-572).
 *
 *-----------------------------------------------------------------------------
 *
 * Original FORTRAN code dgqt.f and destsv.f from MINPACK-2 project.
 *
 * Copyright (C) 1993, Brett M. Averick and Jorge J. Moré (MINPACK-2 Project,
 * Argonne National Laboratory and University of Minnesota).
 *
 * Copyright (C) 1994, Brett M. Averick, Richard Carter and Jorge J. Moré
 * (MINPACK-2 Project, Argonne National Laboratory and University of
 * Minnesota).
 *
 * C-version of GQTPAR is part of OptimPack library is licensed under the MIT
 * "Expat" License.
 *
 * Copyright (c) 2008-2015, Éric Thiébaut
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

#ifndef _OPK_GQTPAR_C
#define _OPK_GQTPAR_C 1

#include <string.h>
#include <math.h>
#include "optimpack.h"

#ifndef NULL
# define NULL ((void *)0)
#endif
#define TRUE  1
#define FALSE 0
#define UPDATE_MAX(VAR, TEMP, EXPR) if ((TEMP = (EXPR)) > VAR) VAR = TEMP

/* First pass for single precision version of the
   routines defined in this file. */
#include "linalg-single.h"
#include __FILE__

/* Second pass for double precision version of the
   routines defined in this file. */
#include "linalg-double.h"
#include __FILE__

#else /* _OPK_GQTPAR_C is defined --------------------------------------------*/

#ifdef OPK_GQT

int
OPK_GQT(const opk_index_t n, real_t a_[], const opk_index_t lda,
        const real_t b_[], const real_t delta,
        const real_t rtol, const real_t atol,
        const opk_index_t itmax, real_t *par_ptr, real_t *fx_ptr,
        real_t x_[], opk_index_t *iter_ptr, real_t z_[],
        real_t wa1_[], real_t wa2_[])
{
  /* Define some macros to mimic FORTRAN indexing. */
#undef wa2
#define wa2(a1) wa2_[a1-1]
#undef wa1
#define wa1(a1) wa1_[a1-1]
#undef z
#define z(a1) z_[a1-1]
#undef x
#define x(a1) x_[a1-1]
#undef b
#define b(a1) b_[a1-1]
#undef a
#define a(a1,a2) a_[a1-1+lda*(a2-1)]

  /* Local constants and variables. */
  const real_t ZERO  = OPK_REALCONST(0.0);
  const real_t ONE   = OPK_REALCONST(1.0);
  const real_t HALF  = OPK_REALCONST(0.5);
  const real_t SMALL = OPK_REALCONST(0.001);
  real_t alpha, anorm, bnorm, prod, rxnorm, rznorm, xnorm;
  real_t par, parc, parf, parl, pars, paru;
  real_t temp, temp1, temp2, temp3;
  opk_index_t indef, j, iter;
  int rednc, info;

  /* Initialization. */
  par = (par_ptr != NULL ? *par_ptr : ZERO);
  parf = ZERO;
  xnorm = ZERO;
  rxnorm = ZERO;
  rednc = FALSE;
  alpha = ZERO; /* avoids compiler warnings */
  rznorm = ZERO; /* avoids compiler warnings */
  for (j = 1; j <= n; ++j) {
    x(j) = ZERO;
    z(j) = ZERO;
  }

  /* Copy the diagonal and save A in its lower triangle. */
  OPK_COPY(n, a_, lda+1, wa1_, 1);
  for (j = 1; j < n; ++j) {
    OPK_COPY(n-j, &a(j,j+1), lda, &a(j+1,j), 1);
  }

  /* Calculate the l1-norm of A, the Gershgorin row sums,
     and the l2-norm of b. */
  anorm = ZERO;
  for (j = 1 ; j <= n ; ++j) {
    temp = OPK_ASUM(n, &a(1,j), 1);
    if (anorm < temp) anorm = temp;
    wa2(j) = temp - FABS(wa1(j));
  }
  bnorm = OPK_NRM2(n, b_, 1);

  /* Calculate a lower bound, PARS, for the domain of the problem.
     Also calculate an upper bound, PARU, and a lower bound, PARL,
     for the Lagrange multiplier. */
  pars = parl = paru = -anorm;
  for (j = 1 ; j <= n ; ++j) {
    UPDATE_MAX(pars, temp1, -wa1(j));
    UPDATE_MAX(parl, temp2, wa2(j) + wa1(j));
    UPDATE_MAX(paru, temp3, wa2(j) - wa1(j));
  }
  /* Compute: parl = max(zero, bnorm/delta - parl, pars) */
  temp = bnorm/delta;
  parl = temp - parl;
  if (parl < pars) parl = pars;
  if (parl < ZERO) parl = ZERO;
  /* Compute: paru = max(zero, bnorm/delta + paru) */
  paru = temp + paru;
  if (paru < ZERO) paru = ZERO;

  /* If the input PAR lies outside of the interval (PARL,PARU),
     set PAR to the closer endpoint. */
  if (par < parl) par = parl;
  if (par > paru) par = paru;

  /* Special case: parl = paru. */
  UPDATE_MAX(paru, temp, (ONE + rtol)*parl);

  /* Beginning of an iteration. */
  info = 0;
  iter = 0;
  while (TRUE) {

    /* Safeguard PAR. */
    if (par <= pars && paru > ZERO) {
      temp = SQRT(parl/paru);
      par = MAX(temp, SMALL)*paru;
    }

    /* Copy the lower triangle of A into its upper triangle and
       compute A + par*I. */
    for (j = 1 ; j < n ; ++j) {
      OPK_COPY(n-j, &a(j+1,j), 1, &a(j,j+1), lda);
    }
    for (j = 1 ; j <= n ; ++j) {
      a(j,j) = wa1(j) + par;
    }

    /* Attempt the Cholesky factorization of A without referencing
       the lower triangular part. */
    /* call dpotrf('u', n, a, lda, indef) -- note that dpotf2 is the
       unblocked version of dpotrf */
    indef = OPK_POTF2(OPK_BLAS_UPPER, n, a_, lda);

    if (indef == 0) {

      /***********************************************
       **  Case 1: A + par*I is positive definite.  **
       ***********************************************/

      /* Compute an approximate solution x and save the
         last value of par with A + par*I positive definite. */
      parf = par;
      OPK_COPY(n, b_, 1, wa2_, 1);
      OPK_TRSV(OPK_BLAS_UPPER, OPK_BLAS_TRANS, OPK_BLAS_NON_UNIT,
               n, a_, lda, wa2_, 1);
      rxnorm = OPK_NRM2(n, wa2_, 1);
      OPK_TRSV(OPK_BLAS_UPPER, OPK_BLAS_NO_TRANS, OPK_BLAS_NON_UNIT,
               n, a_, lda, wa2_, 1);
      OPK_COPY(n, wa2_, 1, x_, 1);
      OPK_SCAL(n, -ONE, x_, 1);
      xnorm = OPK_NRM2(n, x_, 1);

      /* Test for convergence.*/
      if (FABS(xnorm - delta) <= rtol*delta
          || (par == ZERO && xnorm <= (ONE + rtol)*delta)) {
        info = 1;
      }

      /* Compute a direction of negative curvature and use this
         information to improve pars. */
      rznorm = OPK_ESTSV(n, a_, lda, z_);
      UPDATE_MAX(pars, temp, par - rznorm*rznorm);

      /* Compute a negative curvature solution of the form
         x + alpha*z where norm(x + alpha*z) = delta. */
      rednc = FALSE;
      if (xnorm < delta) {

        /* Compute ALPHA. */
        prod = OPK_DOT(n, z_, 1, x_, 1)/delta;
        temp = (delta - xnorm)*((delta + xnorm)/delta);
        alpha = temp/(FABS(prod) + SQRT(prod*prod + temp/delta));
        alpha = SIGN(alpha, prod);

        /* Test to decide if the negative curvature step
           produces a larger reduction than with z = 0. */
        rznorm *= FABS(alpha);
        temp1 = rznorm/delta;
        temp1 *= temp1; /* temp1 = (rznorm/delta)**2 */
        temp2 = xnorm/delta;
        if (temp1 + par*temp2*temp2 <= par) rednc = TRUE;

        /* Test for convergence. */
        temp2 = rxnorm/delta;
        temp2 = par + temp2*temp2; /* temp2 = par + (rxnorm/delta)**2 */
        if (HALF*temp1 <= rtol*(ONE - HALF*rtol)*temp2) {
          info = 1;
        } else if (info == 0 && HALF*temp2 <= (atol/delta)/delta) {
          info = 2;
        } else if (xnorm == ZERO) {
          info = 1;
        }
      }

      /* Compute the Newton correction PARC to PAR. */
      if (xnorm == ZERO) {
        parc = -par;
      } else {
        OPK_COPY(n, x_, 1, wa2_, 1);
        OPK_SCAL(n, ONE/xnorm, wa2_, 1);
        OPK_TRSV(OPK_BLAS_UPPER, OPK_BLAS_TRANS, OPK_BLAS_NON_UNIT,
                 n, a_, lda, wa2_, 1);
        temp = OPK_NRM2(n, wa2_, 1);
        parc = (((xnorm - delta)/delta)/temp)/temp;
      }

      /* Update PARL or PARU.*/
      if (xnorm > delta && parl < par) parl = par;
      if (xnorm < delta && paru > par) paru = par;

    } else {

      /***************************************************
       **  Case 2: A + par*I is not positive definite.  **
       ***************************************************/

      /* Use the rank information from the Cholesky
         decomposition to update PAR. */
      if (indef > 1) {
        /* Restore column indef to A + par*I. */
        OPK_COPY(indef-1, &a(indef,1), lda, &a(1,indef), 1);
        a(indef,indef) = wa1(indef) + par;

        /* Compute PARC. */
        OPK_COPY(indef-1, &a(1,indef), 1, wa2_, 1);
        OPK_TRSV(OPK_BLAS_UPPER, OPK_BLAS_TRANS, OPK_BLAS_NON_UNIT,
                 indef-1, a_, lda, wa2_, 1);
        temp = OPK_NRM2(indef-1, wa2_, 1);
        a(indef,indef) -= temp*temp;
        OPK_TRSV(OPK_BLAS_UPPER, OPK_BLAS_NO_TRANS, OPK_BLAS_NON_UNIT,
                 indef-1, a_, lda, wa2_, 1);
      }
      wa2(indef) = -ONE;
      temp = OPK_NRM2(indef, wa2_, 1);
      parc = -(a(indef,indef)/temp)/temp;
      temp = (parc > ZERO ? par + parc : par);
      if (pars < temp) pars = temp; /* pars = max(pars, par, par+parc) */

      /* If necessary, increase PARU slightly.  This is needed because in
         some exceptional situations PARU is the optimal value of PAR. */
      temp = (ONE + rtol)*pars;
      if (paru < temp) paru = temp; /* paru = max(paru,(one+rtol)*pars) */
    }

    /* Use PARS to update PARL. */
    if (parl < pars) parl = pars;

    /* Test for termination. */
    ++iter;
    if (info == 0) {
      if (paru == ZERO) {
        info = 2;
      } else if (paru <= (ONE + HALF*rtol)*pars) {
        info = 3;
      } else if (iter >= itmax) {
        info = 4;
      }
    }

    /* If exiting, store the best approximation and restore
       the upper triangle of A. */
    if (info != 0) {

      /* Compute the best current estimates for x and f. */
      par = parf;
      if (rednc) {
        if (fx_ptr != NULL) {
          *fx_ptr = -HALF*((rxnorm*rxnorm + par*delta*delta) - rznorm*rznorm);
        }
        OPK_AXPY(n, alpha, z_, 1, x_, 1);
      } else {
        if (fx_ptr != NULL) {
          *fx_ptr = -HALF*(rxnorm*rxnorm + par*xnorm*xnorm);
        }
      }

      /* Restore the upper triangle of A. */
      for (j = 1 ; j < n ; ++j) {
        OPK_COPY(n-j, &a(j+1,j), 1, &a(j,j+1), lda);
      }
      OPK_COPY(n, wa1_, 1, a_, lda+1);

      if (par_ptr != NULL) {
        *par_ptr = par;
      }
      if (iter_ptr != NULL) {
        *iter_ptr = iter;
      }
      return info;
    }

    /* Compute an improved estimate for PAR. */
    if ((par += parc) < parl) par = parl;

  } /* End of an iteration. */

#undef wa2
#undef wa1
#undef z
#undef x
#undef b
#undef a

}

#endif /* OPK_GQT */


#ifdef OPK_ESTSV

real_t
OPK_ESTSV(opk_index_t n, const real_t r_[], opk_index_t ldr, real_t z_[])
{
  /* Define some macros to mimic FORTRAN indexing. */
#undef z
#define z(a1) z_[a1-1]
#undef r
#define r(a1,a2) r_[a1-1+ldr*(a2-1)]

  /* Local constants and variables. */
  const real_t ZERO  = OPK_REALCONST(0.0);
  const real_t ONE   = OPK_REALCONST(1.0);
  const real_t SMALL = OPK_REALCONST(0.01);
  real_t e, s, sm, temp, w, wm, ynorm, znorm;
  opk_index_t i, j;

  /* Clear vector z. */
  if (n > 0) {
    memset(z_, 0, n*sizeof(*z_));
  }

  /* This choice of e makes the algorithm scale invariant. */
  e = FABS(r(1,1));
  if (e == ZERO) {
    z(1) = ONE;
    return ZERO;
  }

  /* Solve R'*y = e.*/
  for (i = 1; i <= n; ++i) {

    /* Scale y. The factor of 0.01 reduces the number of scalings. */
    e = SIGN(e, -z(i));
    if (FABS(e - z(i)) > FABS(r(i,i))) {
      temp = FABS(r(i,i)/(e - z(i)));
      if (temp > SMALL) temp = SMALL;
      OPK_SCAL(n, temp, z_, 1);
      e *= temp;
    }

    /* Determine the two possible choices of y(i). */
    if (r(i,i) == ZERO) {
      w = ONE;
      wm = ONE;
    } else {
      w  =  (e - z(i))/r(i,i);
      wm = -(e + z(i))/r(i,i);
    }

    /* Choose y(i) based on the predicted value of y(j) for j > i. */
    s  = FABS(e - z(i));
    sm = FABS(e + z(i));
    for (j = i + 1; j <= n; ++j) {
      sm += FABS(z(j) + wm*r(i,j));
    }
    if (i < n) {
      OPK_AXPY(n-i, w, &r(i,i+1), ldr, &z(i+1), 1);
      s += OPK_ASUM(n-i, &z(i+1), 1);
    }
    if (s < sm) {
      temp = wm - w;
      w = wm;
      if (i < n && temp != ZERO) {
        OPK_AXPY(n-i, temp, &r(i,i+1), ldr, &z(i+1), 1);
      }
    }
    z(i) = w;
  }
  ynorm = OPK_NRM2(n, &z(1), 1);

  /* Solve R**z = y. */
  for (j = n; j >= 1; --j) {
    if (FABS(z(j)) > FABS(r(j,j))) {
      /* Scale z. */
      temp = FABS(r(j,j)/z(j));
      if (temp > SMALL) temp = SMALL;
      OPK_SCAL(n, temp, &z(1), 1);
      ynorm *= temp;
    }
    if (r(j,j) == ZERO) {
      z(j) = ONE;
    } else {
      z(j) /= r(j,j);
    }
    OPK_AXPY(j-1, -z(j), &r(1,j), 1, &z(1), 1);
  }

  /* Normalize z and return svmin.*/
  znorm = ONE/OPK_NRM2(n, &z(1), 1);
  OPK_SCAL(n, znorm, &z(1), 1);
  return ynorm*znorm;

#undef z
#undef r

}

#endif /* OPK_ESTSV */

#endif /* _OPK_GQTPAR_C */
