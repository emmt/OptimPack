/*
 * lcg.c --
 *
 * Iteratively solve a linear system A.x = b where A is a symmetric positive definite
 * matrix by preconditioned conjugate gradient method with reverse communication.
 *
 *----------------------------------------------------------------------------------------
 *
 * The OptimPack library is licensed under the MIT "Expat" License:
 *
 * Copyright (c) 2003-2014: Éric Thiébaut
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

#ifndef OPK_LCG_C_
#define OPK_LCG_C_ 1

#include <math.h>

#include "optimpack-linalg.h"

#define COPY(n, x, y)    OPK_COPY(n, x, 1, y, 1)
#define AXPY(n, a, x, y) OPK_AXPY(n, a, x, 1, y, 1)
#define DOT(n, x, y)     OPK_DOT(n, x, 1, y, 1)
#define NRM2(n, x)       OPK_NRM2(n, x, 1)
#define SCAL(n, a, x)    OPK_SCAL(n, a, x, 1);
#define ZERO(n, x)       OPK_ZERO(n, x, 1);

#ifndef NULL
# define NULL ((void *)0)
#endif

/* First pass for single precision version of the routines defined in this file. */
#include "linalg-single.h"
#include __FILE__

/* Second pass for double precision version of the routines defined in this file. */
#include "linalg-double.h"
#include __FILE__

#else /* OPK_LCG_C_ defined */

#ifdef OPK_PLCG
void
OPK_PLCG(const opk_index n, real_t p[], real_t q[], real_t r[],
         real_t x[], real_t z[], real_t rho[4], opk_cg_state *state)
{
    const real_t zero = OPK_REALCONST(0.0);
    real_t pq, alpha, beta;
    opk_index i;

    switch (*state) {

    case OPK_CG_START:

        /* Start with no initial guess: fill x with zeros, there is no needs to compute q
           = A.x0 since it is 0. */
        ZERO(n, x);
        rho[0] = rho[1] = rho[2] = rho[3] = zero;
        *state = OPK_CG_NEWX;
        return;

    case OPK_CG_RESTART:

        /* Start or restart with initial x given: copy initial x into p and request caller
           to compute: q = A.p */
        COPY(n, x, p);
        rho[0] = rho[1] = rho[2] = rho[3] = zero;
        *state = OPK_CG_AP;
        return;

    case OPK_CG_NEWX:

        if (z == NULL) {
            /* No preconditioning. Take z = r and jump to next case to compute conjugate
               gradient direction. */
            z = r;
        } else {
            /* Preconditioned version of the algorithm. Request caller to compute
               preconditioned residuals. */
            *state = OPK_CG_PRECOND;
            return;
        }

    case OPK_CG_PRECOND:

        /* Caller has been requested to compute preconditioned residuals. Use conjugate
           gradients recurrence to compute next search direction p. */
        rho[1] = DOT(n, r, z);
        if (rho[1] <= zero) {
            /* If r'.z too small, then algorithm has converged or A is not positive
               definite. */
            *state = (rho[1] < zero ? OPK_CG_NON_CONVEX : OPK_CG_FINISH);
            return;
        }
        if (rho[0] > zero) {
            beta = rho[1]/rho[0];
            for (i = 0; i < n; ++i) {
                p[i] = z[i] + beta*p[i];
            }
        } else {
            beta = zero;
            COPY(n, z, p);
        }
        rho[3] = beta;

        /* Request caller to compute: q = A.p */
        *state = OPK_CG_AP;
        return;

    case OPK_CG_AP:

        if (rho[1] > zero) {
            /* Caller has been requested to compute q = A.p. Compute optimal step size and
               update the variables and the residuals. */
            pq = DOT(n, p, q);
            if (pq <= zero) {
                /* Operator A is not positive definite. */
                *state = OPK_CG_NON_CONVEX;
                return;
            }
            alpha = rho[1]/pq; /* optimal step size */
            rho[2] = alpha; /* memorize optimal step size */
            if (alpha == zero) {
                /* If alpha too small, then algorithm has converged. */
                *state = OPK_CG_FINISH;
                return;
            }
            AXPY(n,  alpha, p, x);
            AXPY(n, -alpha, q, r);
            rho[0] = rho[1];
        } else {
            /* Caller has been requested to compute q = A.x0, update the residuals. */
            AXPY(n, -1.0, q, r);
        }

        /* Variables X and residuals R available for inspection. */
        *state = OPK_CG_NEWX;
        return;

    case OPK_CG_FINISH:
    case OPK_CG_NON_CONVEX:

        /* If state is OPK_CG_FINISH, the caller can restart the algorithm with final x,
           state set to OPK_CG_RESTART and r set to b. */
        return;

    default:

        /* There must be something wrong... */
        *state = OPK_CG_ERROR;
        return;
    }

}
#endif /* OPK_PLCG */

#if defined(OPK_PLCG) && defined(OPK_LCG)
void
OPK_LCG(const opk_index n, real_t p[], real_t q[], real_t r[],
        real_t x[], real_t rho[4], opk_cg_state *state)
{
    OPK_PLCG(n, p, q, r, x, (real_t *)0, rho, state);
}
#endif /* OPK_PLCG and OPK_LCG */


#if defined(OPK_TRCG)

/*
 *  ... RHO[0] = RHO_PREV
 *  ... RHO[1] = RHO
 *  ... RHO[2] = ALPHA
 *  ... RHO[3] = BETA
 *  ... RHO[4] = |X|
 */
void
OPK_TRCG(const opk_index n, real_t p[], const real_t q[], real_t r[],
         real_t x[], const real_t z[], const real_t delta,
         real_t rho[5], opk_cg_state *state)
{
    const real_t zero = OPK_REALCONST(0.0);
    const real_t one  = OPK_REALCONST(1.0);
    real_t a, b, c, d, e, xn, sum, tmp;
    real_t pq, alpha, beta;
    long i;

    if (delta <= zero) {
        *state = OPK_CG_ERROR;
        return;
    }

    switch (*state) {

    case OPK_CG_START:

        /* Start with no initial guess: fill x with zeros, there is no needs to compute
           q = A.x0 since it is 0. */
        ZERO(n, x);
        rho[0] = rho[1] = rho[2] = rho[3] = rho[4] = zero;
        *state = OPK_CG_NEWX;
        return;

    case OPK_CG_RESTART:

        /* Start or restart with initial x given: copy initial x into p and request caller
           to compute: q = A.p */
        rho[0] = rho[1] = rho[2] = rho[3] = zero;
        xn = NRM2(n, x);
        if (xn >= delta) {
            if (xn > delta) {
                SCAL(n, delta/xn, x);
            }
            rho[4] = delta; /* |X| */
            *state = OPK_CG_TRUNCATED;
        } else {
            rho[4] = xn; /* |X| */
            COPY(n, x, p);
            *state = OPK_CG_AP;
        }
        return;

    case OPK_CG_NEWX:

        if (z == NULL) {
            /* No preconditioning. Take z = r and jump to next case to compute conjugate
               gradient direction. */
            z = r;
        } else {
            /* Preconditioned version of the algorithm. Request caller to compute
               preconditioned residuals. */
            *state = OPK_CG_PRECOND;
            return;
        }

    case OPK_CG_PRECOND:

        /* Caller has been requested to compute preconditioned residuals. Use conjugate
           gradients recurrence to compute next search direction p. */
        rho[1] = DOT(n, r, z);
        if (rho[1] <= zero) {
            /* If r'.z too small, then algorithm has converged or preconditioner is not
               positive definite. */
            *state = (rho[1] < zero ? OPK_CG_NON_CONVEX : OPK_CG_FINISH);
            return;
        }
        if (rho[0] > zero) {
            beta = rho[1]/rho[0];
            for (i = 0; i < n; ++i) {
                p[i] = z[i] + beta*p[i];
            }
        } else {
            beta = zero;
            COPY(n, z, p);
        }
        rho[3] = beta;

        /* Request caller to compute: q = A.p */
        *state = OPK_CG_AP;
        return;

    case OPK_CG_AP:

        if (rho[1] > zero) {
            /* Caller has been requested to compute q = A.p. Compute optimal step size and
               update the variables and the residuals. */
            pq = DOT(n, p, q);
            if (pq > zero) {
                alpha = rho[1]/pq; /* optimal step size */
                rho[2] = alpha; /* memorize optimal step size */
                if (alpha == zero) {
                    /* If alpha too small, then algorithm has converged. */
                    *state = OPK_CG_FINISH;
                    return;
                }
                sum = zero;
                for (i = 0; i < n; ++i) {
                    tmp = x[i] + alpha*p[i];
                    sum += tmp*tmp;
                }
                xn = SQRT(sum);
                if (xn <= delta) {
                    /* Optimal step not too long, take it. */
                    AXPY(n,  alpha, p, x);
                    AXPY(n, -alpha, q, r);
                    rho[0] = rho[1];
                    rho[4] = xn;
                    *state = (xn < delta ? OPK_CG_NEWX : OPK_CG_TRUNCATED);
                    return;
                }
            }

            /* Operator A is not positive definite (P'.A.P < 0) or optimal step leads us
               outside the trust region. In these cases, we take a truncated step X +
               ALPHA*P along P so that |X + ALPHA*P| = DELTA with ALPHA >= 0. This amounts
               to find the positive root of: A*ALPHA^2 + 2*B*ALPHA + C with: A = |P|^2, B
               = X'.P, and C = |X|^2 - DELTA^2. Note that the reduced discriminant is D =
               B*B - A*C which expands to D = (DELTA^2 - |X|^2*sin(THETA)^2)*|P|^2 with
               THETA the angle between X and P, as DELTA > |X|, then D > 0 must hold. */
            a = DOT(n, p, p);
            if (a <= zero) {
                /* This can only occurs if P = 0 (and thus A = 0). It is probably due to
                   rounding errors and the algorithm should be restarted. */
                *state = OPK_CG_FINISH;
                return;
            }
            b = DOT(n, x, p);
            c = (rho[4] + delta)*(rho[4] - delta); /* RHO[4] = |X| */
            if (c >= zero) {
                /* This can only occurs if the caller has modified DELTA or RHO[4], which
                   stores the norm of X, and thus is considered as an error. */
                *state = OPK_CG_ERROR;
                return;
            }
            /* Normalize the coefficients to avoid overflows. */
            d = FABS(a);
            if ((e = FABS(b)) > d) d = e;
            if ((e = FABS(c)) > d) d = e;
            d = one/d;
            a *= d;
            b *= d;
            c *= d;
            /* Compute the reduced discriminant. */
            d = b*b - a*c;
            if (d > zero) {
                /* The polynomial has two real roots of opposite signs (because C/A < 0 by
                   construction), ALPHA is the positive one. Compute this root avoiding
                   numerical errors (by adding numbers of same sign). */
                e = SQRT(d);
                if (b >= zero) {
                    alpha = -c/(e + b);
                } else {
                    alpha = (e - b)/a;
                }
            } else {
                /* D > 0 must hold.  There must be something wrong... */
                *state = OPK_CG_ERROR;
                return;
            }
            if (alpha > zero) {
                AXPY(n,  alpha, p, x);
                AXPY(n, -alpha, q, r);
            }
            rho[0] = rho[1];
            rho[2] = alpha; /* memorize optimal step size */
            rho[4] = delta;
            *state = OPK_CG_TRUNCATED;
            return;

        } else {

            /* Caller has been requested to compute q = A.x0, update the residuals. */
            for (i = 0; i < n; ++i) {
                r[i] -= q[i];
            }
            rho[2] = zero; /* ALPHA = 0 (no step has been taken yet) */
            rho[3] = zero; /* BETA = 0 (ditto) */

        }

        /* Variables X and residuals R available for inspection. */
        *state = OPK_CG_NEWX;
        return;


    case OPK_CG_FINISH:
    case OPK_CG_NON_CONVEX:
    case OPK_CG_TRUNCATED:

        /* The caller can restart the algorithm with final 'x', 'state' set to
           OPK_CG_RESTART and 'r' set to 'b' and, maybe, a larger 'delta'. */
        return;

    default:

        /* There must be something wrong... */
        *state = OPK_CG_ERROR;
        return;
    }

}
#endif /* OPK_TRCG */

#endif /* OPK_LCG_C_ defined */
