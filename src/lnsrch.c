/*
 * lnsrch.c --
 *
 * Linesearch routines for OptimPack library.  Implements Armijo linesearch,
 * Moré & Thuente inexact linesearch and nonmonotone linesearch methods.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2003-2014 Éric Thiébaut
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
#include <string.h>
#include <math.h>
#include <errno.h>

#include "optimpack-private.h"

#define TRUE   OPK_TRUE
#define FALSE  OPK_FALSE

#define MAX(a,b)  OPK_MAX(a,b)
#define MIN(a,b)  OPK_MIN(a,b)

#define ROUND_UP(a,b)   OPK_ROUND_UP(a,b)

#if (OPK_LNSRCH_SEARCH != 0)
# error OPK_LNSRCH_SEARCH != 0
#endif

/*---------------------------------------------------------------------------*/
/* PRIVATE ROUTINES */

static int
non_finite(double x)
{
  return (isnan(x) || isinf(x));
}

/*---------------------------------------------------------------------------*/
/* UNIFIED INTERFACE FOR LINE SEARCH */

void
finalize_line_search(opk_object_t* obj)
{
  opk_lnsrch_t* ws = (opk_lnsrch_t*)obj;
  if (ws->ops->finalize != NULL) {
    ws->ops->finalize(ws);
  }
}

opk_lnsrch_t*
opk_allocate_line_search(opk_lnsrch_operations_t *ops,
                         size_t size)
{
  opk_lnsrch_t* ls;

  if (ops == NULL || ops->start == NULL || ops->iterate == NULL) {
    errno = EFAULT;
    return NULL;
  }
  if (ops->start == NULL || ops->iterate == NULL) {
    errno = EINVAL;
    return NULL;
  }
  size = OPK_MAX(size, sizeof(opk_lnsrch_t));
  ls = (opk_lnsrch_t*)opk_allocate_object(finalize_line_search, size);
  if (ls != NULL) {
    ls->ops = ops;
    ls->status = OPK_LNSRCH_ERROR_NOT_STARTED;
  }
  return ls;
}

/* after an error or convergence, you must call opk_lnsrch_start */

int
opk_lnsrch_start(opk_lnsrch_t* ls, double f0, double g0,
                 double stp, double stpmin, double stpmax)
{
  if (ls == NULL) {
    return OPK_LNSRCH_ERROR_ILLEGAL_ADDRESS;
  }
  if (stpmin < 0.0) {
    ls->status = OPK_LNSRCH_ERROR_STPMIN_LT_ZERO;
  } else if (stpmin > stpmax) {
    ls->status = OPK_LNSRCH_ERROR_STPMIN_GT_STPMAX;
  } else if (stp < stpmin) {
    ls->status = OPK_LNSRCH_ERROR_STP_LT_STPMIN;
  } else if (stp > stpmax) {
    ls->status = OPK_LNSRCH_ERROR_STP_GT_STPMAX;
  } else if (g0 >= 0.0) {
    ls->status = OPK_LNSRCH_ERROR_INITIAL_DERIVATIVE_GE_ZERO;
  } else {
    ls->stp = stp;
    ls->stpmin = stpmin;
    ls->stpmax = stpmax;
    ls->finit = f0;
    ls->ginit = g0;
    ls->status = ls->ops->start(ls);
  }
  return ls->status;
}

int
opk_lnsrch_iterate(opk_lnsrch_t* ls, double* stp_ptr,
                   double f1, double g1)
{
  if (ls == NULL || stp_ptr == NULL) {
    return OPK_LNSRCH_ERROR_ILLEGAL_ADDRESS;
  }
  if (ls->status == OPK_LNSRCH_SEARCH) {
    if (*stp_ptr != ls->stp) {
      ls->status = OPK_LNSRCH_ERROR_STP_CHANGED;
    } else {
      ls->status = ls->ops->iterate(ls, stp_ptr, f1, g1);
      if (*stp_ptr > ls->stpmax) {
        if (ls->stp >= ls->stpmax) {
          ls->status = OPK_LNSRCH_WARNING_STP_EQ_STPMAX;
        }
        *stp_ptr = ls->stpmax;
      } else if (*stp_ptr < ls->stpmin) {
        if (ls->stp <= ls->stpmin) {
          ls->status = OPK_LNSRCH_WARNING_STP_EQ_STPMIN;
        }
        *stp_ptr = ls->stpmin;
      }
      ls->stp = *stp_ptr;
    }
  } else {
    ls->status = OPK_LNSRCH_ERROR_NOT_STARTED;
  }
  return ls->status;
}

double
opk_lnsrch_get_step(const opk_lnsrch_t* ls)
{
  return (ls != NULL && ls->status == OPK_LNSRCH_SEARCH ? ls->stp : -1.0);
}

int
opk_lnsrch_get_status(const opk_lnsrch_t* ls)
{
  return (ls != NULL ? ls->status : OPK_LNSRCH_ERROR_ILLEGAL_ADDRESS);
}

const char*
opk_lnsrch_message(int status)
{
  switch (status) {
  case OPK_LNSRCH_ERROR_ILLEGAL_ADDRESS:
    return "Illegal address in line search";
  case OPK_LNSRCH_ERROR_CORRUPTED_WORKSPACE:
    return "Corrupted line search workspace";
  case OPK_LNSRCH_ERROR_BAD_WORKSPACE:
    return "Bad line search workspace";
  case OPK_LNSRCH_ERROR_STP_CHANGED:
    return "Line search step modified by caller";
  case OPK_LNSRCH_ERROR_STP_OUTSIDE_BRACKET:
    return "Line search step outside bracket";
  case OPK_LNSRCH_ERROR_NOT_A_DESCENT:
    return "Line search direction is not a descent";
  case OPK_LNSRCH_ERROR_STPMIN_GT_STPMAX:
    return "Minimum line search step greater than maximum step";
  case OPK_LNSRCH_ERROR_STPMIN_LT_ZERO:
    return "Minimum line search step less than zero";
  case OPK_LNSRCH_ERROR_STP_LT_STPMIN:
    return "Line search step less than minimum step";
  case OPK_LNSRCH_ERROR_STP_GT_STPMAX:
    return "Line search step greater than maximum step";
  case OPK_LNSRCH_ERROR_INITIAL_DERIVATIVE_GE_ZERO:
    return "Initial derivative along line search greater or equal zero";
  case OPK_LNSRCH_ERROR_NOT_STARTED:
    return "Line search not started";
  case OPK_LNSRCH_SEARCH:
    return "Line search in progress";
  case OPK_LNSRCH_CONVERGENCE:
    return "Line search has converged";
  case OPK_LNSRCH_WARNING_ROUNDING_ERRORS_PREVENT_PROGRESS:
    return "Rounding errors prevent progress in line search";
  case OPK_LNSRCH_WARNING_XTOL_TEST_SATISFIED:
    return "OPK_LNSRCH_WARNING_XTOL_TEST_SATISFIED";
  case OPK_LNSRCH_WARNING_STP_EQ_STPMAX:
    return "Line search step at upper bound";
  case OPK_LNSRCH_WARNING_STP_EQ_STPMIN:
    return "Line search step at lower bound";
  default:
    return NULL;
  }
}

opk_bool_t
opk_lnsrch_has_errors(const opk_lnsrch_t* ls)
{
  return (opk_lnsrch_get_status(ls) < 0);
}

opk_bool_t
opk_lnsrch_has_warnings(const opk_lnsrch_t* ls)
{
  return (opk_lnsrch_get_status(ls) > OPK_LNSRCH_CONVERGENCE);
}

opk_bool_t
opk_lnsrch_converged(const opk_lnsrch_t* ls)
{
  return (opk_lnsrch_get_status(ls) == OPK_LNSRCH_CONVERGENCE);
}

opk_bool_t
opk_lnsrch_finished(const opk_lnsrch_t* ls)
{
  return (opk_lnsrch_get_status(ls) != OPK_LNSRCH_SEARCH);
}

opk_bool_t
opk_lnsrch_use_deriv(const opk_lnsrch_t* ls)
{
  return ls->ops->use_deriv;
}

/*---------------------------------------------------------------------------*/
/* ARMIJO (BACKTRACKING) LINE SEARCH */

/* Sub-type for Armijo line search. */
typedef struct _backtrack_lnsrch backtrack_lnsrch_t;
struct _backtrack_lnsrch {
  opk_lnsrch_t base; /* Base type (must be the first member). */
  double ftol;
};

static int
backtrack_start(opk_lnsrch_t* _ws)
{
  return OPK_LNSRCH_SEARCH;
}

static int
backtrack_iterate(opk_lnsrch_t* ls,
                  double* stp_ptr, double f, double g)
{
  backtrack_lnsrch_t* bls = (backtrack_lnsrch_t*)ls;
  int status = OPK_LNSRCH_SEARCH;

  if (f <= ls->finit + bls->ftol*(*stp_ptr)*ls->ginit) {
    /* First Wolfe conditions satisfied. */
    status = OPK_LNSRCH_CONVERGENCE;
  } else {
    /* Take a bisection step unless already at the lower bound. */
    if (*stp_ptr <= ls->stpmin) {
      status = OPK_LNSRCH_WARNING_STP_EQ_STPMIN;
    } else {
      *stp_ptr = (*stp_ptr + ls->stpmin)*0.5;
    }
    /* Safeguard the step. */
    if (*stp_ptr < ls->stpmin) {
      *stp_ptr = ls->stpmin;
    }
  }
  return status;
}

static opk_lnsrch_operations_t backtrack_operations = {
  NULL,
  backtrack_start,
  backtrack_iterate,
  FALSE
};

opk_lnsrch_t*
opk_lnsrch_new_backtrack(double ftol)
{
  opk_lnsrch_t* ls;

  if (ftol <= 0.0) {
    errno = EINVAL;
    return NULL;
  }
  ls = opk_allocate_line_search(&backtrack_operations,
                                sizeof(backtrack_lnsrch_t));
  if (ls != NULL) {
    backtrack_lnsrch_t* bls = (backtrack_lnsrch_t*)ls;
    bls->ftol = ftol;
  }
  return ls;
}

/*---------------------------------------------------------------------------*/
/* NONMONOTONE LINE SEARCH */

/* Sub-type for nonmonotone line search. */
typedef struct _nonmonotone_lnsrch nonmonotone_lnsrch_t;
struct _nonmonotone_lnsrch {
  opk_lnsrch_t base; /**< Base type (must be the first member). */
  double sigma1;     /**< Lower steplength bound to trigger bisection. */
  double sigma2;     /**< Upper steplength relative bound to trigger bisection. */
  double ftol;       /**< Parameter for the function reduction criterion. */
  double fmax;       /**< Maximum function value for the past M steps. */
  double* fsav;      /**< Function values for M last accepted steps. */
  opk_index_t m;     /**< Number of previous steps to remember. */
  opk_index_t mp;    /**< Number of steps since starting. */
};

static int
nonmonotone_start(opk_lnsrch_t* ls)
{
  nonmonotone_lnsrch_t* nmls = (nonmonotone_lnsrch_t*)ls;
  opk_index_t j, n;

  /* Save function value. */
  nmls->fsav[nmls->mp%nmls->m] = ls->finit;
  ++nmls->mp;

  /* Get the worst function value among the N last steps. */
  n = MIN(nmls->mp, nmls->m);
  nmls->fmax = nmls->fsav[0];
  for (j = 1; j < n; ++j) {
    if (nmls->fsav[j] > nmls->fmax) {
      nmls->fmax = nmls->fsav[j];
    }
  }

  return OPK_LNSRCH_SEARCH;
}

#if 0
static void
nonmonotone_reset(opk_lnsrch_t* ls)
{
  nonmonotone_lnsrch_t* nmls = (nonmonotone_lnsrch_t*)ls;
  nmls->mp = 0;
}
#endif

static int
nonmonotone_iterate(opk_lnsrch_t* ls,
                    double* stp_ptr, double f, double g)
{
  nonmonotone_lnsrch_t* nmls = (nonmonotone_lnsrch_t*)ls;
  double alpha, delta, q, r;

  /* Check whether Armijo-like condition satisfied. */
  alpha = *stp_ptr;   /* current steplength */
  delta = ls->ginit;  /* directional derivative at alpha=0 */
  if (f <= nmls->fmax + nmls->ftol*alpha*delta) {
    /* Convergence criterion satisfied. */
    return OPK_LNSRCH_CONVERGENCE;
  }

  /* Check whether step is already at the lower bound. */
  if (alpha <= ls->stpmin) {
    *stp_ptr = ls->stpmin;
    return OPK_LNSRCH_WARNING_STP_EQ_STPMIN;
  }

  /* Attempt to use safeguarded quadratic interpolation to find a better step.
     The optimal steplength estimated by quadratic interpolation is q/r and
     r > 0 must hold for the quadratic approximation to be strictly convex. */
  q = -delta*alpha*alpha;
  r = (f - ls->finit - alpha*delta)*2.0;
  if (r > 0.0 && nmls->sigma1*r <= q && q <= nmls->sigma2*alpha*r) {
    /* Quadratic approximation is strictly convex and its minimum is within
       the bounds.  Take the quadratic interpolation step. */
    alpha = q/r;
  } else {
    /* Take the bisection step. */
    alpha = (alpha + ls->stpmin)/2.0;
  }

  /* Safeguard the step. */
  alpha = MAX(alpha, ls->stpmin);
  *stp_ptr = alpha;
  return (alpha > 0.0 ? OPK_LNSRCH_SEARCH : OPK_LNSRCH_WARNING_STP_EQ_STPMIN);
}

static opk_lnsrch_operations_t nonmonotone_operations = {
  NULL,
  nonmonotone_start,
  nonmonotone_iterate,
  FALSE
};

opk_lnsrch_t*
opk_lnsrch_new_nonmonotone(opk_index_t m, double ftol,
                           double sigma1, double sigma2)
{
  opk_lnsrch_t* ls;
  size_t size, offset;

  if (non_finite(ftol) || non_finite(sigma1) || non_finite(sigma1) ||
      ftol <= 0.0 || sigma1 <= 0.0 || sigma1 >= sigma2 || sigma2 >= 1.0 ||
      m < 1) {
    errno = EINVAL;
    return NULL;
  }
  offset = ROUND_UP(sizeof(nonmonotone_lnsrch_t), sizeof(double));
  size = offset + m*sizeof(double);
  ls = opk_allocate_line_search(&nonmonotone_operations, size);
  if (ls != NULL) {
    nonmonotone_lnsrch_t* nmls = (nonmonotone_lnsrch_t*)ls;
    nmls->ftol = ftol;
    nmls->sigma1 = sigma1;
    nmls->sigma2 = sigma2;
    nmls->fsav = (double*)(((char*)ls) + offset);
    nmls->m = m;
    nmls->mp = 0;
  }
  return ls;
}

/*---------------------------------------------------------------------------*/
/* MORÉ AND THUENTE CUBIC LINE SEARCH */

/* Sub-type for Moré and Thuente cubic line search. */
typedef struct _csrch_lnsrch csrch_lnsrch_t;
struct _csrch_lnsrch {
  opk_lnsrch_t base; /* Base type (must be the first member). */

  /* Convergence parameters. */
  double ftol;
  double gtol;
  double xtol;

  /* GTEST is used to check for Wolfe conditions. */
  double gtest;

  /* The variables STX, FX, GX contain the values of the step, function, and
     derivative at the best step. */
  double stx, fx, gx;

  /* The variables STY, FY, GY contain the value of the step, function, and
     derivative at STY. */
  double sty, fy, gy;

  /* Parameters to track the interval where to seek for the step. */
  double stmin;
  double stmax;
  double width;
  double width1;
  opk_bool_t brackt;

  /* The algorithm has two different stages. */
  int stage;
};

static csrch_lnsrch_t*
csrch_get_workspace(opk_lnsrch_t* ls);

static int
csrch_start(opk_lnsrch_t* ls)
{
  csrch_lnsrch_t* cls;

  cls = csrch_get_workspace(ls);
  if (cls == NULL) {
    return OPK_LNSRCH_ERROR_CORRUPTED_WORKSPACE;
  }

  /* Convergence threshold for this step. */
  cls->gtest = cls->ftol*cls->base.ginit;

  /* Initialize parameters for the interval of search. */
  cls->stmin = cls->base.stpmin; /* FIXME: ??? */
  cls->stmax = cls->base.stpmax;
  cls->width = cls->base.stpmax - cls->base.stpmin;
  cls->width1 = cls->width/0.5;
  cls->brackt = FALSE;

  /* The variables STX, FX, GX contain the values of the step,
     function, and derivative at the best step. */
  cls->stx = 0.0;
  cls->fx = cls->base.finit;
  cls->gx = cls->base.ginit;

  /* The variables STY, FY, GY contain the value of the step,
     function, and derivative at STY. */
  cls->sty = 0.0;
  cls->fy = cls->base.finit;
  cls->gy = cls->base.ginit;

  cls->stage = 1;
  return OPK_LNSRCH_SEARCH;
}

static int
csrch_iterate(opk_lnsrch_t* _ws,
              double* stp_ptr, double f1, double g1)
{
  double ftest;
  csrch_lnsrch_t* cls;
  int result;

  cls = csrch_get_workspace(_ws);
  if (cls == NULL) {
    return OPK_LNSRCH_ERROR_CORRUPTED_WORKSPACE;
  }

  /* Test for convergence. */
  ftest = cls->base.finit + (*stp_ptr)*cls->gtest;
  if (f1 <= ftest && fabs(g1) <= -cls->gtol*cls->base.ginit) {
    /* Strong Wolfe conditions satisfied. */
    return OPK_LNSRCH_CONVERGENCE;
  }

  /* Test for warnings. */
  if (*stp_ptr == cls->base.stpmin && (f1 > ftest || g1 >= cls->gtest)) {
    return OPK_LNSRCH_WARNING_STP_EQ_STPMIN;
  }
  if (*stp_ptr == cls->base.stpmax && f1 <= ftest && g1 <= cls->gtest) {
    return OPK_LNSRCH_WARNING_STP_EQ_STPMAX;
  }
  if (cls->brackt && cls->stmax - cls->stmin <= cls->xtol * cls->stmax) {
    return OPK_LNSRCH_WARNING_XTOL_TEST_SATISFIED;
  }
  if (cls->brackt && (*stp_ptr <= cls->stmin || *stp_ptr >= cls->stmax)) {
    return OPK_LNSRCH_WARNING_ROUNDING_ERRORS_PREVENT_PROGRESS;
  }

  /* If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
     algorithm enters the second stage. */
  if (cls->stage == 1 && f1 <= ftest && g1 >= 0.0) {
    cls->stage = 2;
  }

  /* A modified function is used to predict the step during the first stage if
     a lower function value has been obtained but the decrease is not
     sufficient. */
  if (cls->stage == 1 && f1 <= cls->fx && f1 > ftest) {
    /* Define the modified function and derivative values and call CSTEP to
       update STX, STY, and to compute the new step.  Then restore the
       function and derivative values for F.*/
    double fm = f1 - *stp_ptr * cls->gtest;
    double fxm = cls->fx - cls->stx * cls->gtest;
    double fym = cls->fy - cls->sty * cls->gtest;
    double gm = g1 - cls->gtest;
    double gxm = cls->gx - cls->gtest;
    double gym = cls->gy - cls->gtest;
    result = opk_cstep(&cls->stx, &fxm, &gxm,
                       &cls->sty, &fym, &gym,
                       stp_ptr, fm, gm,
                       &cls->brackt, cls->stmin, cls->stmax);
    if (result < 0) {
      return result;
    }
    cls->fx = fxm + cls->stx * cls->gtest;
    cls->fy = fym + cls->sty * cls->gtest;
    cls->gx = gxm + cls->gtest;
    cls->gy = gym + cls->gtest;
  } else {
    /* Call CSTEP to update STX, STY, and to compute the new step. */
    result = opk_cstep(&cls->stx, &cls->fx, &cls->gx,
                       &cls->sty, &cls->fy, &cls->gy,
                       stp_ptr, f1, g1,
                       &cls->brackt, cls->stmin, cls->stmax);
    if (result < 0) {
      return result;
    }
  }

  /* Decide if a bisection step is needed. */
  if (cls->brackt) {
    double new_width = fabs(cls->sty - cls->stx);
    if (new_width >= 0.66*cls->width1) {
      *stp_ptr = cls->stx + 0.5*(cls->sty - cls->stx);
    }
    cls->width1 = cls->width;
    cls->width = new_width;
  }

  /* Set the minimum and maximum steps allowed for stp. */
  if (cls->brackt) {
    cls->stmin = MIN(cls->stx, cls->sty);
    cls->stmax = MAX(cls->stx, cls->sty);
  } else {
    cls->stmin = *stp_ptr + (*stp_ptr - cls->stx)*1.1;
    cls->stmax = *stp_ptr + (*stp_ptr - cls->stx)*4.0;
  }

  /* Force the step to be within the bounds stpmax and stpmin. */
  *stp_ptr = MAX(*stp_ptr, cls->base.stpmin);
  *stp_ptr = MIN(*stp_ptr, cls->base.stpmax);

  /* If further progress is not possible, let stp be the best
     point obtained during the search. */
  if (cls->brackt && (*stp_ptr <= cls->stmin || *stp_ptr >= cls->stmax
                     || cls->stmax - cls->stmin <= cls->xtol * cls->stmax)) {
    *stp_ptr = cls->stx;
  }

  /* Obtain another function and derivative. */
  return OPK_LNSRCH_SEARCH;
}

static opk_lnsrch_operations_t csrch_operations = {
  NULL,
  csrch_start,
  csrch_iterate,
  TRUE
};


static csrch_lnsrch_t*
csrch_get_workspace(opk_lnsrch_t* ls)
{
  if (ls->ops == &csrch_operations) {
    return (csrch_lnsrch_t*)ls;
  } else {
    return NULL;
  }
}

opk_lnsrch_t*
opk_lnsrch_new_csrch(double ftol, double gtol, double xtol)
{
  opk_lnsrch_t* _ws;

  if (ftol < 0.0) {
    /* ERROR: FTOL .LT. ZERO */
    errno = EINVAL;
    return NULL;
  }
  if (gtol < 0.0) {
    /* ERROR: GTOL .LT. ZERO */
    errno = EINVAL;
    return NULL;
  }
  if (xtol < 0.0) {
    /* ERROR: XTOL .LT. ZERO */
    errno = EINVAL;
    return NULL;
  }
  _ws = opk_allocate_line_search(&csrch_operations,
                                 sizeof(csrch_lnsrch_t));
  if (_ws != NULL) {
    csrch_lnsrch_t* ws = (csrch_lnsrch_t*)_ws;
    ws->ftol = ftol;
    ws->gtol = gtol;
    ws->xtol = xtol;
    ws->stage = 0;
  }
  return _ws;
}

/**
 * Find a step that satisfies a sufficient decrease condition and a curvature
 * condition.
 *
 * This subroutine finds a step that satisfies a sufficient decrease condition
 * and a curvature condition.
 *
 * Each call of the subroutine updates an interval with endpoints stx and
 * sty. The interval is initially chosen so that it contains a minimizer of the
 * modified function
 *
 *     psi(stp) = f(stp) - f(0) - ftol*stp*f'(0).
 *
 * If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the interval is chosen
 * so that it contains a minimizer of f.
 *
 * The algorithm is designed to find a step that satisfies the sufficient
 * decrease condition
 *
 *     f(stp) <= f(0) + ftol*stp*f'(0),
 *
 * and the curvature condition
 *
 *     abs(f'(stp)) <= gtol*abs(f'(0)).
 *
 * If ftol is less than gtol and if, for example, the function is bounded
 * below, then there is always a step which satisfies both conditions.
 *
 * If no step can be found that satisfies both conditions, then the algorithm
 * stops with a warning. In this case stp only satisfies the sufficient
 * decrease condition.
 *
 * A typical invocation of dcsrch has the following outline:
 *
 *     Evaluate the function at stp = 0.0d0; store in f.
 *     Evaluate the gradient at stp = 0.0d0; store in g.
 *     Choose a starting step stp.
 *
 *     task = 'START'
 *  10 continue
 *        call dcsrch(stp,f,g,ftol,gtol,xtol,task,stpmin,stpmax,
 *    +               isave,dsave)
 *        if (task .eq. 'FG') then
 *           Evaluate the function and the gradient at stp
 *           go to 10
 *           end if
 *
 *     NOTE: The user must not alter work arrays between calls.
 *
 *     The subroutine statement is
 *
 *       subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,
 *                         task,isave,dsave)
 *     where
 *
 *       stp is a double precision variable.
 *         On entry stp is the current estimate of a satisfactory
 *            step. On initial entry, a positive initial estimate
 *            must be provided.
 *         On exit stp is the current estimate of a satisfactory step
 *            if task = 'FG'. If task = 'CONV' then stp satisfies
 *            the sufficient decrease and curvature condition.
 *
 *       f is a double precision variable.
 *         On initial entry f is the value of the function at 0.
 *            On subsequent entries f is the value of the
 *            function at stp.
 *         On exit f is the value of the function at stp.
 *
 *       g is a double precision variable.
 *         On initial entry g is the derivative of the function at 0.
 *            On subsequent entries g is the derivative of the
 *            function at stp.
 *         On exit g is the derivative of the function at stp.
 *
 *       ftol is a double precision variable.
 *         On entry ftol specifies a nonnegative tolerance for the
 *            sufficient decrease condition.
 *         On exit ftol is unchanged.
 *
 *       gtol is a double precision variable.
 *         On entry gtol specifies a nonnegative tolerance for the
 *            curvature condition.
 *         On exit gtol is unchanged.
 *
 *       xtol is a double precision variable.
 *         On entry xtol specifies a nonnegative relative tolerance
 *            for an acceptable step. The subroutine exits with a
 *            warning if the relative difference between sty and stx
 *            is less than xtol.
 *         On exit xtol is unchanged.
 *
 *       task is a character variable of length at least 60.
 *         On initial entry task must be set to 'START'.
 *         On exit task indicates the required action:
 *
 *            If task(1:2) = 'FG' then evaluate the function and
 *            derivative at stp and call dcsrch again.
 *
 *            If task(1:4) = 'CONV' then the search is successful.
 *
 *            If task(1:4) = 'WARN' then the subroutine is not able
 *            to satisfy the convergence conditions. The exit value of
 *            stp contains the best point found during the search.
 *
 *            If task(1:5) = 'ERROR' then there is an error in the
 *            input arguments.
 *
 *         On exit with convergence, a warning or an error, the
 *            variable task contains additional information.
 *
 *       stpmin is a double precision variable.
 *         On entry stpmin is a nonnegative lower bound for the step.
 *         On exit stpmin is unchanged.
 *
 *       stpmax is a double precision variable.
 *         On entry stpmax is a nonnegative upper bound for the step.
 *         On exit stpmax is unchanged.
 *
 *       isave is an integer work array of dimension 2.
 *
 *       dsave is a double precision work array of dimension 13.
 *
 *     Subprograms called
 *
 *       MINPACK-2 ... dcstep
 *
 *     MINPACK-1 Project. June 1983.
 *     Argonne National Laboratory.
 *     Jorge J. More' and David J. Thuente.
 *
 *     MINPACK-2 Project. November 1993.
 *     Argonne National Laboratory and University of Minnesota.
 *     Brett M. Averick, Richard G. Carter, and Jorge J. More'.
 */

/*---------------------------------------------------------------------------*/

static double max3(double val1, double val2, double val3)
{
  double result = val1;
  if (val2 > result) result = val2;
  if (val3 > result) result = val3;
  return result;
}

/* FIXME: add a reference to the documentation and check history. */

/**
 * Compute a safeguarded step for a line search.
 *
 * This routine computes a safeguarded step for a search procedure and updates
 * an interval that contains a step that satisfies a sufficient decrease and a
 * curvature condition.
 *
 * The parameter stx contains the step with the least function value. If
 * brackt is true, then a minimizer has been bracketed in an interval
 * with endpoints stx and sty.  The parameter stp contains the current step.
 * The subroutine assumes that if brackt is true then:
 *
 *       min(stx,sty) < stp < max(stx,sty),
 *
 * and that the derivative at stx is negative in the direction of the step;
 * that is:
 *
 *      dx*(stp - stx) < 0
 *
 *
 * @param stx_ptr
 *        A pointer to stx.  On entry, stx is the best step obtained so far
 *        and is an endpoint of the interval that contains the minimizer.  On
 *        exit, stx is the updated best step.
 *
 * @param fx_ptr
 *        A pointer to fx.  On entry, fx is the function at stx.  On exit, fx
 *        is the function at stx.
 *
 * @param dx_ptr
 *        A pointer to dx.  On entry, dx is the derivative of the function at
 *        stx. The derivative must be negative in the direction of the step,
 *        that is, dx and stp - stx must have opposite signs.  On exit, dx is
 *        the derivative of the function at stx.
 *
 * @param sty_ptr
 *        A pointer to sty.  On entry, sty is the second endpoint of the
 *        interval that contains the minimizer.  On exit, sty is the updated
 *        endpoint of the interval that contains the minimizer.
 *
 * @param fy_ptr
 *        A pointer to fy.  On entry, fy is the function at sty.  On exit, fy
 *        is the function at sty.
 *
 * @param dy_ptr
 *        A pointer to dy.  On entry, dy is the derivative of the function at
 *        sty.  On exit, dy is the derivative of the function at sty.
 *
 * @param stp_ptr
 *        A pointer to stp. On entry, stp is the current step. If the value at
 *        brackt_ptr is true, then on input, stp must be between stx and sty.
 *        On exit, stp is a new trial step.
 *
 * @param fp
 *        The value of the function at stp on entry.
 *
 * @param dp
 *        The the derivative of the function at stp on entry.
 *
 * @param brackt_ptr
 *        A pointer to the boolean variable brackt.  On entry, brackt
 *        specifies if a minimizer has been bracketed.  Initially brackt must
 *        be set to false On exit, brackt specifies if a minimizer has been
 *        bracketed.  When a minimizer is bracketed brackt is set to true.
 *
 * @return A strictly negative value on error, a strictly positive value on
 * success.  The value returned on success is between 1 and 4 and corresponds
 * to one of the four possible cases.  The value returned on error is one of:
 *
 *    OPK_LNSRCH_ERROR_STP_OUTSIDE_BRACKET
 *    if brackt is true but the step is outside the bracket endpoints.
 *
 *    OPK_LNSRCH_ERROR_NOT_A_DESCENT
 *    if the descent condition is violated.
 *
 *    OPK_LNSRCH_ERROR_STPMIN_GT_STPMAX
 *    if stpmin > stpmax.
 *
 * @history
 * MINPACK-1 Project. June 1983
 * Argonne National Laboratory.
 * Jorge J. Moré and David J. Thuente.
 *
 * MINPACK-2 Project. November 1993.
 * Argonne National Laboratory and University of Minnesota.
 * Brett M. Averick and Jorge J. Moré.
 */
int opk_cstep(double* stx_ptr, double* fx_ptr, double* dx_ptr,
              double* sty_ptr, double* fy_ptr, double* dy_ptr,
              double* stp_ptr, double  fp,     double  dp,
              int* brackt_ptr, double stpmin, double stpmax)
{
  /* Constants. */
  const double ZERO = 0.0;
  const double TWO = 2.0;
  const double THREE = 3.0;

  /* Get values of input/output variables. */
  double stx = *stx_ptr, fx = *fx_ptr, dx = *dx_ptr;
  double sty = *sty_ptr, fy = *fy_ptr, dy = *dy_ptr;
  double stp = *stp_ptr;

  /* Local variables. */
  double gamma, theta, p, q, r, s, temp;
  double stpc; /* cubic step */
  double stpq; /* quadratic step */
  double stpf;
  int opposite, result;

  /* Check the input parameters for errors. */
  if ((*brackt_ptr) && (stx < sty ? (stp <= stx || stp >= sty)
                                  : (stp <= sty || stp >= stx))) {
    return OPK_LNSRCH_ERROR_STP_OUTSIDE_BRACKET;
  } else if (dx*(stp - stx) >= ZERO) {
    return OPK_LNSRCH_ERROR_NOT_A_DESCENT;
  } else if (stpmin > stpmax) {
    return OPK_LNSRCH_ERROR_STPMIN_GT_STPMAX;
  }

  /* Determine if the derivatives have opposite signs. */
  opposite = ((dp < ZERO && dx > ZERO) || (dp > ZERO && dx < ZERO));

  if (fp > fx) {
    /* First case.  A higher function value.  The minimum is bracketed.  If
       the cubic step is closer to STX than the quadratic step, the cubic step
       is taken, otherwise the average of the cubic and quadratic steps is
       taken. */
    result = 1;
    *brackt_ptr = TRUE;
    theta = THREE*(fx - fp)/(stp - stx) + dx + dp;
    s = max3(fabs(theta), fabs(dx), fabs(dp));
    temp = theta/s;
    gamma = s*sqrt(temp*temp - (dx/s)*(dp/s));
    if (stp < stx) gamma = -gamma;
    p =  (gamma - dx) + theta;
    q = ((gamma - dx) + gamma) + dp;
    r = p/q;
    stpc = stx + r*(stp - stx);
    stpq = stx + ((dx/((fx - fp)/(stp - stx) + dx))/TWO)*(stp - stx);
    if (fabs(stpc - stx) < fabs(stpq - stx)) {
      stpf = stpc;
    } else {
      stpf = stpc + (stpq - stpc)/TWO;
    }
  } else if (opposite) {
    /* Second case.  A lower function value and derivatives of opposite sign.
       The minimum is bracketed.  If the cubic step is farther from STP than
       the secant (quadratic) step, the cubic step is taken, otherwise the
       secant step is taken. */
    result = 2;
    *brackt_ptr = TRUE;
    theta = THREE*(fx - fp)/(stp - stx) + dx + dp;
    s = max3(fabs(theta), fabs(dx), fabs(dp));
    temp = theta/s;
    gamma = s*sqrt(temp*temp - (dx/s)*(dp/s));
    if (stp > stx) gamma = -gamma;
    p =  (gamma - dp) + theta;
    q = ((gamma - dp) + gamma) + dx;
    r = p/q;
    stpc = stp + r*(stx - stp);
    stpq = stp + (dp/(dp - dx))*(stx - stp);
    if (fabs(stpc - stp) > fabs(stpq - stp)) {
      stpf = stpc;
    } else {
      stpf = stpq;
    }
  } else if (fabs(dp) < fabs(dx)) {
    /* Third case.  A lower function value, derivatives of the same sign, and
       the magnitude of the derivative decreases.  The cubic step is computed
       only if the cubic tends to infinity in the direction of the step or if
       the minimum of the cubic is beyond STP.  Otherwise the cubic step is
       defined to be the secant step.  The case GAMMA = 0 only arises if the
       cubic does not tend to infinity in the direction of the step. */
    result = 3;
    theta = THREE*(fx - fp)/(stp - stx) + dx + dp;
    s = max3(fabs(theta), fabs(dx), fabs(dp));
    temp = theta/s;
    temp = temp*temp - (dx/s)*(dp/s);
    if (temp > ZERO) {
      gamma = s*sqrt(temp);
      if (stp > stx) gamma = -gamma;
    } else {
      gamma = ZERO;
    }
    p = (gamma - dp) + theta;
    q = (gamma + (dx - dp)) + gamma;
    r = p/q;
    if (r < ZERO && gamma != ZERO) {
      stpc = stp + r*(stx - stp);
    } else if (stp > stx) {
      stpc = stpmax;
    } else {
      stpc = stpmin;
    }
    stpq = stp + (dp/(dp - dx))*(stx - stp);

    if (*brackt_ptr) {
      /* A minimizer has been bracketed.  If the cubic step is closer to STP
         than the secant step, the cubic step is taken, otherwise the secant
         step is taken. */
      if (fabs(stpc - stp) < fabs(stpq - stp)) {
        stpf = stpc;
      } else {
        stpf = stpq;
      }
      temp = stp + 0.66*(sty - stp);
      if (stp > stx ? stpf > temp : stpf < temp) {
        stpf = temp;
      }
    } else {
      /* A minimizer has not been bracketed. If the cubic step is farther from
         stp than the secant step, the cubic step is taken, otherwise the
         secant step is taken. */
      if (fabs(stpc - stp) > fabs(stpq - stp)) {
        stpf = stpc;
      } else {
        stpf = stpq;
      }
      if (stpf > stpmax) stpf = stpmax;
      if (stpf < stpmin) stpf = stpmin;
    }
  } else {
    /* Fourth case.  A lower function value, derivatives of the same sign, and
       the magnitude of the derivative does not decrease.  If the minimum is
       not bracketed, the step is either STPMIN or STPMAX, otherwise the cubic
       step is taken. */
    result = 4;
    if (*brackt_ptr) {
      theta = THREE*(fp - fy)/(sty - stp) + dy + dp;
      s = max3(fabs(theta), fabs(dy), fabs(dp));
      temp = theta/s;
      gamma = s*sqrt(temp*temp - (dy/s)*(dp/s));
      if (stp > sty) gamma = -gamma;
      p =  (gamma - dp) + theta;
      q = ((gamma - dp) + gamma) + dy;
      r = p/q;
      stpc = stp + r*(sty - stp);
      stpf = stpc;
    } else if (stp > stx) {
      stpf = stpmax;
    } else {
      stpf = stpmin;
    }
  }

  /* Update the interval which contains a minimizer. */
  if (fp > fx) {
    *sty_ptr = stp;
    *fy_ptr = fp;
    *dy_ptr = dp;
  } else {
    if (opposite) {
      *sty_ptr = stx;
      *fy_ptr = fx;
      *dy_ptr = dx;
    }
    *stx_ptr = stp;
    *fx_ptr = fp;
    *dx_ptr = dp;
  }

  /* Store the new safeguarded step. */
  *stp_ptr = stpf;
  return result;
}
