/*
 * lbfgs.c --
 *
 * Limited memory variable metric method for OptimPack library.  The method has
 * been first described by Nocedal (1980) under the name VMLM (for Variable
 * Metric Limited Memory) and popularized under the name LBFGS by Byrd et
 * al. (1994).
 *
 *  - Nocedal, J. "Updating Quasi-Newton Matrices with Limited Storage,"
 *    Mathematics of Computation, Vol. 35, pp. 773-782 (1980)
 *
 *  - Byrd, R. H., Nocedal, J. & Schnabel, R. B. "Representations of
 *    quasi-Newton matrices and their use in limited memory methods,"
 *    Mathematical Programming, Vol. 63, pp. 129-156 (1994).
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of OptimPack (https://github.com/emmt/OptimPack).
 *
 * Copyright (C) 2015 Éric Thiébaut
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

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>

#include "optimpack-private.h"

#define SAVE_MEMORY 1 /* save memory to use 2 + 2m vectors instead of 4 + 2m */

#define TRUE                          OPK_TRUE
#define FALSE                         OPK_FALSE
#define MAX(a, b)                     OPK_MAX(a, b)
#define MIN(a, b)                     OPK_MIN(a, b)
#define ROUND_UP(a, b)                OPK_ROUND_UP(a,b )
#define ADDRESS(type, base, offset)   OPK_ADDRESS(type, base, offset)

/*---------------------------------------------------------------------------*/
/* VARIABLE METRIC OPTIMIZER */

/* Default parameters for the line search.  These values are from original
   LBFGS algorithm by Jorge Nocedal but other values may be more suitable. */
static const double STPMIN = 1E-20;
static const double STPMAX = 1E+20;
static const double SFTOL = 1E-4;
static const double SGTOL = 1E-1;
static const double SXTOL = DBL_EPSILON;

/* Default parameters for the global convergence. */
static const double GRTOL = 1E-6;
static const double GATOL = 0.0;

/* Other default parameters. */
static const double DELTA   = 5e-2;
static const double EPSILON = 1e-2;

struct _opk_lbfgs {
  opk_object_t base;       /**< Base type (must be the first member). */
  double delta;            /**< Relative size for a small step. */
  double epsilon;          /**< Threshold to accept descent direction. */
  double grtol;            /**< Relative threshold for the norm or the gradient
                                (relative to GINIT the norm of the initial
                                gradient) for convergence. */
  double gatol;            /**< Absolute threshold for the norm or the gradient
                                for convergence. */
  double ginit;            /**< Euclidean norm or the initial gradient. */
  double f0;               /**< Function value at x0. */
  double gnorm;            /**< Euclidean norm of the gradient at the last
                                tested step. */
  double stp;              /**< Current step length. */
  double stpmin;           /**< Relative mimimum step length. */
  double stpmax;           /**< Relative maximum step length. */
  opk_vspace_t* vspace;    /**< Variable space. */
  opk_lnsrch_t* lnsrch;    /**< Line search method. */
  opk_vector_t* x0;        /**< Variables at the start of the line search. */
  opk_vector_t* g0;        /**< Gradient at x0. */
  opk_vector_t* d;         /**< Anti-search direction; a iterate is computed
                                as: x1 = x0 - stp*d with stp > 0. */
  opk_index_t evaluations; /**< Number of functions (and gradients)
                                evaluations. */
  opk_index_t iterations;  /**< Number of iterations (successful steps
                                taken). */
  opk_index_t restarts;    /**< Number of LBFGS recurrence restarts. */
  opk_status_t status;     /**< Last error. */
  opk_task_t task;         /**< Pending task. */

  /* Limited memory BFGS approximation of the Hessian of the objective
     function. */
  double gamma;            /**< Scale factor to approximate inverse Hessian. */
  opk_vector_t** s;        /**< Storage for variable differences. */
  opk_vector_t** y;        /**< Storage for gradient differences. */
  double* beta;            /**< Workspace to save <d,s>/<s,y> */
  double* rho;             /**< Workspace to save 1/<s,y> */
  opk_index_t m;           /**< Maximum number of memorized steps. */
  opk_index_t mp;          /**< Actual number of memorized steps
                                (0 <= mp <= m). */
  opk_index_t updates;     /**< Number of BFGS updates since start. */
};

static double
max3(double a1, double a2, double a3)
{
  if (a3 >= a2) {
    return (a3 >= a1 ? a3 : a1);
  } else {
    return (a2 >= a1 ? a2 : a1);
  }
}

static int
non_finite(double x)
{
  return (isnan(x) || isinf(x));
}

static void
finalize_lbfgs(opk_object_t* obj)
{
  opk_lbfgs_t* opt = (opk_lbfgs_t*)obj;
  opk_index_t k;

  /* Drop references to all objects (neither opt->x0 nor opt->g0 which are weak
     references to specific vectors in opt->s and opt->y). */
  OPK_DROP(opt->vspace);
  OPK_DROP(opt->lnsrch);
#if ! SAVE_MEMORY
  OPK_DROP(opt->x0);
  OPK_DROP(opt->g0);
#endif
  OPK_DROP(opt->d);
  if (opt->s != NULL) {
    for (k = 0; k < opt->m; ++k) {
      OPK_DROP(opt->s[k]);
    }
  }
  if (opt->y != NULL) {
    for (k = 0; k < opt->m; ++k) {
      OPK_DROP(opt->y[k]);
    }
  }
}

/* The following function returns the index where is stored the (updates-j)-th
   correction pair.  Argument j must be in the inclusive range [0:mp] with mp
   the actual number of saved corrections.  At any moment,
   0 ≤ mp ≤ min(m,updates); thus updates - j ≥ 0. */
static opk_index_t
slot(const opk_lbfgs_t* opt, opk_index_t j)
{
  return (opt->updates - j)%opt->m;
}

#if 0
opk_vector_t*
opk_get_lbfgs_s(opk_lbfgs_t* opt, opk_index_t k)
{
  return (0 <= k && k <= opt->mp ? opt->s[slot(opt, k)] : NULL);
}

opk_vector_t*
opk_get_lbfgs_y(opk_lbfgs_t* opt, opk_index_t k)
{
  return (0 <= k && k <= opt->mp ? opt->y[slot(opt, k)] : NULL);
}
#endif

static opk_task_t
success(opk_lbfgs_t* opt, opk_task_t task)
{
  opt->status = OPK_SUCCESS;
  opt->task = task;
  return task;
}

static opk_task_t
failure(opk_lbfgs_t* opt, opk_status_t status)
{
  opt->status = status;
  opt->task = OPK_TASK_ERROR;
  return opt->task;
}

opk_lbfgs_t*
opk_new_lbfgs_optimizer_with_line_search(opk_vspace_t* space,
                                         opk_index_t m,
                                         opk_lnsrch_t* lnsrch)
{
  opk_lbfgs_t* opt;
  size_t s_offset, y_offset, beta_offset, rho_offset, size;
  opk_index_t k;

  /* Check the input arguments for errors. */
  if (space == NULL || lnsrch == NULL) {
    errno = EFAULT;
    return NULL;
  }
  if (space->size < 1 || m < 1) {
    errno = EINVAL;
    return NULL;
  }
  if (m > space->size) {
    m = space->size;
  }

  /* Allocate enough memory for the workspace and its arrays. */
  s_offset = ROUND_UP(sizeof(opk_lbfgs_t), sizeof(opk_vector_t*));
  y_offset = s_offset + m*sizeof(opk_vector_t*);
  beta_offset = ROUND_UP(y_offset + m*sizeof(opk_vector_t*), sizeof(double));
  rho_offset = beta_offset + m*sizeof(double);
  size = rho_offset + m*sizeof(double);
  opt = (opk_lbfgs_t*)opk_allocate_object(finalize_lbfgs, size);
  if (opt == NULL) {
    return NULL;
  }
  opt->s    = ADDRESS(opk_vector_t*, opt,    s_offset);
  opt->y    = ADDRESS(opk_vector_t*, opt,    y_offset);
  opt->beta = ADDRESS(double,        opt, beta_offset);
  opt->rho  = ADDRESS(double,        opt,  rho_offset);
  opt->m = m;
  for (k = 0; k < m; ++k) {
    opt->s[k] = opk_vcreate(space);
    if (opt->s[k] == NULL) {
      goto error;
    }
    opt->y[k] = opk_vcreate(space);
    if (opt->y[k] == NULL) {
      goto error;
    }
  }
  opt->vspace = OPK_HOLD_VSPACE(space);
  opt->lnsrch = OPK_HOLD_LNSRCH(lnsrch);
#if ! SAVE_MEMORY
  opt->x0 = opk_vcreate(space);
  if (opt->x0 == NULL) {
    goto error;
  }
  opt->g0 = opk_vcreate(space);
  if (opt->g0 == NULL) {
    goto error;
  }
#endif
  opt->d = opk_vcreate(space);
  if (opt->d == NULL) {
    goto error;
  }
  opk_set_lbfgs_options(opt, NULL);
  opt->gamma = 1.0;
  return opt;

 error:
  OPK_DROP(opt);
  return NULL;
}

opk_lbfgs_t*
opk_new_lbfgs_optimizer(opk_vspace_t* space, opk_index_t m)
{
  opk_lnsrch_t* lnsrch;
  opk_lbfgs_t* opt;

  lnsrch = opk_lnsrch_new_csrch(SFTOL, SGTOL, SXTOL);
  if (lnsrch == NULL) {
    return NULL;
  }
  opt = opk_new_lbfgs_optimizer_with_line_search(space, m, lnsrch);
  OPK_DROP(lnsrch); /* the line search is now owned by the optimizer */
  return opt;
}

opk_task_t
opk_start_lbfgs(opk_lbfgs_t* opt, opk_vector_t* x)
{
  opt->iterations = 0;
  opt->evaluations = 0;
  opt->restarts = 0;
  opt->updates = 0;
  opt->mp = 0;
  return success(opt, OPK_TASK_COMPUTE_FG);
}


#define S(k)     opt->s[k]
#define Y(k)     opt->y[k]
#define BETA(k)  opt->beta[k]
#define RHO(k)   opt->rho[k]
#define UPDATE(dst, alpha, x)            opk_vaxpby(dst, 1, dst, alpha, x)
#define COMBINE(dst, alpha, x, beta, y)  opk_vaxpby(dst, alpha, x, beta, y)
#define SCALE(x, alpha)                  opk_vscale(x, alpha, x)
#define DOT(x, y)                        opk_vdot(x, y)

opk_task_t
opk_iterate_lbfgs(opk_lbfgs_t* opt, opk_vector_t* x,
                  double f, opk_vector_t* g)
{
  double dtg, yty, sty;
  opk_index_t j, k;
  opk_status_t status;
  opk_lnsrch_task_t lnsrch_task;

  switch (opt->task) {

  case OPK_TASK_COMPUTE_FG:

    /* Caller has computed the function value and the gradient at the current
       point. */
    ++opt->evaluations;

    /* Check for global convergence. */
    opt->gnorm = opk_vnorm2(g);
    if (opt->evaluations == 1) {
      opt->ginit = opt->gnorm;
    }
    if (opt->gnorm <= max3(0.0, opt->gatol, opt->grtol*opt->ginit)) {
      return success(opt, OPK_TASK_FINAL_X);
    }

    if (opt->evaluations > 1) {
      /* A line search is in progress, check whether it has converged. */
      if (opk_lnsrch_use_deriv(opt->lnsrch)) {
        /* Compute directional derivative. */
        dtg = -opk_vdot(opt->d, g);
      } else {
        /* Line search does not need directional derivative. */
        dtg = 0;
      }
      lnsrch_task = opk_lnsrch_iterate(opt->lnsrch, &opt->stp, f, dtg);
      if (lnsrch_task == OPK_LNSRCH_SEARCH) {
        /* Line search has not yet converged, break to compute a new trial
           point along the search direction. */
        break;
      }
      if (lnsrch_task != OPK_LNSRCH_CONVERGENCE) {
        status = opk_lnsrch_get_status(opt->lnsrch);
        if (lnsrch_task != OPK_LNSRCH_WARNING ||
            status != OPK_ROUNDING_ERRORS_PREVENT_PROGRESS) {
          return failure(opt, status);
        }
      }
      ++opt->iterations;
    }
    return success(opt, OPK_TASK_NEW_X);

  case OPK_TASK_NEW_X:
  case OPK_TASK_FINAL_X:

    if (opt->iterations >= 1) {
      /* Update L-BFGS approximation of the Hessian. */
      k = slot(opt, 0);
      opk_vaxpby(S(k), 1, x, -1, opt->x0);
      opk_vaxpby(Y(k), 1, g, -1, opt->g0);
      sty = DOT(Y(k), S(k));
      if (sty <= 0) {
        RHO(k) = 0;
      } else {
        RHO(k) = 1/sty;
        yty = DOT(Y(k), Y(k));
        if (yty > 0) {
          opt->gamma = sty/yty;
        }
      }
      ++opt->updates;
      if (opt->mp < opt->m) {
        ++opt->mp;
      }
    }

    /* Compute a search direction.  We take care of checking whether -D is a
       sufficient descent direction (here D is a sufficient ascent direction).
       As shown by Zoutendijk, this is true if cos(theta) = (D/|D|)'.(G/|G|) is
       larger or equal EPSILON > 0, where G is the gradient at X and D the
       descent direction. */
    dtg = 0;
    if (opt->mp >= 1) {
      /* Apply the L-BFGS Strang's two-loop recursion to compute a search
         direction. */
      opt->gamma = 0;
      opk_vcopy(opt->d, g);
      for (j = 1; j <= opt->mp; ++j) {
        k = slot(opt, j);
        BETA(k) = RHO(k)*DOT(opt->d, S(k));
        UPDATE(opt->d, -BETA(k), Y(k));
      }
      if (opt->gamma > 0 && opt->gamma != 1) {
        /* Apply initial inverse Hessian approximation. */
        SCALE(opt->d, opt->gamma);
      }
      for (j = opt->mp; j >= 1; --j) {
        k = slot(opt, j);
        if (RHO(k) > 0) {
          UPDATE(opt->d, BETA(k) - RHO(k)*DOT(opt->d, Y(k)), S(k));
        }
      }

      /* Check whether the algorithm has produced a sufficient descent
         direction. */
      dtg = -DOT(opt->d, g);
      if (opt->epsilon > 0 &&
          dtg > -opt->epsilon*opk_vnorm2(opt->d)*opt->gnorm) {
        /* Set DTG to zero to indicate that we do not have a sufficient
           descent direction. */
        dtg = 0;
      }
    }
    if (dtg < 0) {
      /* Use search direction produced by L-BFGS recursion and an initial unit
         step. */
      opt->stp = 1.0;
    } else {
      /* Use steepest ascent. */
      if (opt->mp > 0) {
        /* L-BFGS recursion did not produce a sufficient descent direction. */
        ++opt->restarts;
        opt->mp = 0;
      }
      dtg = -opt->gnorm*opt->gnorm;
      if (f != 0) {
        opt->stp = 2*fabs(f/dtg);
      } else {
        double dnorm = opt->gnorm;
        double xnorm = opk_vnorm2(x);
        if (xnorm > 0) {
          opt->stp = opt->delta*xnorm/dnorm;
        } else {
          opt->stp = opt->delta/dnorm;
        }
      }
    }

    /* Save current point. */
#if SAVE_MEMORY
    k = slot(opt, 0);
    opt->x0 = S(k); /* weak reference */
    opt->g0 = Y(k); /* weak reference */
    if (opt->mp == opt->m) {
      --opt->mp;
    }
#endif
    opk_vcopy(opt->x0, x);
    opk_vcopy(opt->g0, g);
    opt->f0 = f;

    /* Start the line search and break to compute the first trial point along
       the line search. */
    if (opk_lnsrch_start(opt->lnsrch, f, dtg, opt->stp,
                         opt->stp*opt->stpmin,
                         opt->stp*opt->stpmax) != OPK_LNSRCH_SEARCH) {
      return failure(opt, opk_lnsrch_get_status(opt->lnsrch));
    }
    break;

  default:

    /* There is probably something wrong. */
    return opt->task;

  }

  /* Compute a trial point along the line search. */
  opk_vaxpby(x, 1, opt->x0, -opt->stp, opt->d);
  return success(opt, OPK_TASK_COMPUTE_FG);
}

opk_task_t
opk_get_lbfgs_task(opk_lbfgs_t* opt)
{
  return opt->task;
}

opk_status_t
opk_get_lbfgs_status(opk_lbfgs_t* opt)
{
  return opt->status;
}

opk_index_t
opk_get_lbfgs_iterations(opk_lbfgs_t* opt)
{
  return opt->iterations;
}

opk_index_t
opk_get_lbfgs_evaluations(opk_lbfgs_t* opt)
{
  return opt->evaluations;
}

opk_index_t
opk_get_lbfgs_restarts(opk_lbfgs_t* opt)
{
  return opt->restarts;
}

double
opk_get_lbfgs_step(opk_lbfgs_t* opt)
{
  return (opt == NULL ? -1.0 : opt->stp);
}

opk_status_t
opk_get_lbfgs_options(opk_lbfgs_options_t* dst, const opk_lbfgs_t* src)
{
  if (dst == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  if (src == NULL) {
    /* Get default options. */
    dst->delta = DELTA;
    dst->epsilon = EPSILON;
    dst->grtol = GRTOL;
    dst->gatol = GATOL;
    dst->stpmin = STPMIN;
    dst->stpmax = STPMAX;
  } else {
    /* Get current options. */
    dst->delta = src->delta;
    dst->epsilon = src->epsilon;
    dst->grtol = src->grtol;
    dst->gatol = src->gatol;
    dst->stpmin = src->stpmin;
    dst->stpmax = src->stpmax;
  }
  return OPK_SUCCESS;
}

opk_status_t
opk_set_lbfgs_options(opk_lbfgs_t* dst, const opk_lbfgs_options_t* src)
{
  if (dst == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  if (src == NULL) {
    /* Set default options. */
    dst->delta = DELTA;
    dst->epsilon = EPSILON;
    dst->grtol = GRTOL;
    dst->gatol = GATOL;
    dst->stpmin = STPMIN;
    dst->stpmax = STPMAX;
  } else {
    /* Check and set given options. */
    if (non_finite(src->gatol) || src->gatol < 0 ||
        non_finite(src->grtol) || src->grtol < 0 ||
        non_finite(src->delta) || src->delta <= 0 ||
        non_finite(src->epsilon) || src->epsilon < 0 ||
        non_finite(src->stpmin) || src->stpmin < 0 ||
        non_finite(src->stpmax) || src->stpmax <=  src->stpmin) {
      return OPK_INVALID_ARGUMENT;
    }
    dst->delta = src->delta;
    dst->epsilon = src->epsilon;
    dst->grtol = src->grtol;
    dst->gatol = src->gatol;
    dst->stpmin = src->stpmin;
    dst->stpmax = src->stpmax;
  }
  return OPK_SUCCESS;
}
