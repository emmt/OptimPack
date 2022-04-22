/*
 * vmlmb.c --
 *
 * Limited memory variable metric method possibly with optional bound
 * constraints for OptimPack library.
 *
 * Three distinct methods are implemented: VMLM (also known as L-BFGS), VMLM-B
 * and BLMVM.  All methods use a limited memory BFGS approximation of the
 * inverse Hessian which has been first described by Nocedal (1980) under the
 * name VMLM (for Variable Metric Limited Memory) and popularized under the
 * name L-BFGS (for Limited memory BFGS) by Lyu and Nocedal (1989).  VMLM-B and
 * BLMVM implement bound constraints.  BLMVM, invented by Benson & Moré (2001),
 * uses the projected gradient to update the BFGS approximation.  VMLM-B,
 * described by Thiébaut (2002), apply the L-BFGS recursion in the subspace of
 * the free variables.  These two latter methods are the competitors of
 * L-BFGS-B.
 *
 * References:
 *  - J. Nocedal, "Updating Quasi-Newton Matrices with Limited Storage",
 *    Mathematics of Computation, Vol. 35, pp. 773-782 (1980)
 *
 *  - D. Liu and J. Nocedal, "On the limited memory BFGS method for large scale
 *    optimization", Mathematical Programming B 45, 503-528 (1989).
 *
 *  - S.J. Benson, & J. Moré, "A Limited-Memory Variable-Metric Algorithm for
 *    Bound-Constrained Minimization", Tech. Report ANL/MCS-P909-0901,
 *    Mathematics and Computer Science Division, Argonne National Laboratory
 *    (2001).
 *
 *  - É. Thiébaut, "Optimization issues in blind deconvolution algorithms",
 *    SPIE Conf. Astronomical Data Analysis II, Vol. 4847, 174-183 (2002).
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of OptimPack (https://github.com/emmt/OptimPack).
 *
 * Copyright (C) 2002, 2015-2019 Éric Thiébaut
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

#define OPK_UPDATE_AS_SOON_AS_POSSIBLE (1 << 0) /**< Update LBFGS approximation
                                                     as soon as possible
                                                     (before returnning
                                                     OPK_FINAL_X or OPK_NEW_X).
                                                     This is useful for
                                                     hierarchical optimization
                                                     in an ADMM scheme. */
#define OPK_EMULATE_BLMVM              (1 << 1) /**< Emulate Benson & Moré
                                                     BLMVM method. */



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
static const double SFTOL = 0.0001;
static const double SGTOL = 0.9;
static const double SXTOL = DBL_EPSILON;
static const double SAMIN = 0.1;

/* Default parameters for the global convergence. */
static const double GRTOL = 1E-6;
static const double GATOL = 0.0;

/* Other default parameters. */
static const double DELTA   = 5e-2;
static const double EPSILON = 0.0;
#if 0 /* unused stuff */
static const opk_bfgs_scaling SCALING = OPK_SCALING_OREN_SPEDICATO;
#endif

struct opk_vmlmb_ {
  opk_object base;       /**< Base type (must be the first member). */
  double delta;            /**< Relative size for a small step. */
  double epsilon;          /**< Threshold to accept descent direction. */
  double grtol;            /**< Relative threshold for the norm or the
                                (projected) gradient (relative to GINIT the
                                norm of the initial (projected) gradient) for
                                convergence. */
  double gatol;            /**< Absolute threshold for the norm or the
                                (projected) gradient for convergence. */
  double ginit;            /**< Euclidean norm or the initial (projected)
                            * gradient. */
  double f0;               /**< Function value at x0. */
  double gnorm;            /**< Euclidean norm of the (projected) gradient at
                                the last tested step. */
  double stp;              /**< Current step length. */
  double stpmin;           /**< Relative minimum step length. */
  double stpmax;           /**< Relative maximum step length. */
  double bsmin;            /**< Step size to the first encountered bound. */
  opk_vspace* vspace;    /**< Variable space. */
  opk_lnsrch* lnsrch;    /**< Line search method. */
  opk_vector* x0;        /**< Variables at the start of the line search. */
  opk_vector* g0;        /**< Gradient at x0. */
  opk_vector* d;         /**< Anti-search direction; a iterate is computed
                                as: x1 = x0 - stp*d with stp > 0. */
  opk_vector* w;         /**< Vector whose elements are set to 1 for free
                                variables and 0 otherwise. */
  opk_vector* gp;        /**< Work vector used to store the effective step
                                size and the projected gradient (for the BLMVM
                                method). */
  opk_convexset* box;    /**< Convex set implementing the box constraints (or
                            *   `NULL`). */
  int method;              /**< The method to use/emulate. */
  opk_index evaluations; /**< Number of functions (and gradients)
                                evaluations. */
  opk_index iterations;  /**< Number of iterations (successful steps
                                taken). */
  opk_index restarts;    /**< Number of LBFGS recurrence restarts. */
  opk_status status;     /**< Last error. */
  opk_task task;         /**< Pending task. */
  opk_bool save_memory;  /**< Save some memory? */

  /* Limited memory BFGS approximation of the Hessian of the objective
     function. */
  double gamma;            /**< Scale factor to approximate inverse Hessian. */
  opk_vector** s;        /**< Storage for variable differences. */
  opk_vector** y;        /**< Storage for gradient differences. */
  double* alpha;           /**< Workspace to save <d,s>/<s,y> */
  double* rho;             /**< Workspace to save 1/<s,y> */
  opk_index m;           /**< Maximum number of memorized steps. */
  opk_index mp;          /**< Actual number of memorized steps
                                (0 <= mp <= m). */
  opk_index updates;     /**< Number of BFGS updates since start. */
};

static double
max3(double a1, double a2, double a3)
{
  return (a3 >= a2 ? MAX(a1, a3) : MAX(a1, a2));
}

static int
non_finite(double x)
{
  return (isnan(x) || isinf(x));
}

static void
finalize_vmlmb(opk_object* obj)
{
  opk_vmlmb* opt = (opk_vmlmb*)obj;
  opk_index k;

  /* Drop references to all objects (neither opt->x0 nor opt->g0 which are weak
     references to specific vectors in opt->s and opt->y). */
  OPK_DROP(opt->vspace);
  OPK_DROP(opt->lnsrch);
  if (! opt->save_memory) {
    OPK_DROP(opt->x0);
    OPK_DROP(opt->g0);
  }
  OPK_DROP(opt->d);
  OPK_DROP(opt->w);
  OPK_DROP(opt->gp);
  OPK_DROP(opt->box);
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

/* The following function returns the index where is stored the
   (updates-j)-th correction pair.  Argument j must be in the
   inclusive range [0:mp] with mp the actual number of saved
   corrections.  At any moment, `0 ≤ mp ≤ min(m,updates)`; thus
   `updates - j ≥ 0`. */
static opk_index
slot(const opk_vmlmb* opt, opk_index j)
{
  return (opt->updates - j)%opt->m;
}

static opk_task
success(opk_vmlmb* opt, opk_task task)
{
  opt->status = OPK_SUCCESS;
  opt->task = task;
  return task;
}

static opk_task
failure(opk_vmlmb* opt, opk_status status)
{
  opt->status = status;
  opt->task = OPK_TASK_ERROR;
  return opt->task;
}

opk_vmlmb*
opk_new_vmlmb_optimizer(const opk_vmlmb_options* opts,
                        opk_vspace* space,
                        opk_lnsrch* lnsrch,
                        opk_convexset* box)
{
  opk_vmlmb_options options;
  opk_vmlmb* opt;
  size_t s_offset, y_offset, alpha_offset, rho_offset, size;
  opk_index k, m;

  /* Check options. */
  if (opts == NULL) {
    /* Use default options. */
    opk_get_vmlmb_default_options(&options);
    opts = &options;
  }
  if (opk_check_vmlmb_options(opts) != OPK_SUCCESS) {
    errno = EINVAL;
    return NULL;
  }
  m = opts->mem;

  /* Check other input arguments for errors. */
  if (space == NULL) {
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
  if (box != NULL && box->space != space) {
    errno = EINVAL;
    return NULL;
  }

  /* Allocate enough memory for the workspace and its arrays. */
  s_offset = ROUND_UP(sizeof(opk_vmlmb), sizeof(opk_vector*));
  y_offset = s_offset + m*sizeof(opk_vector*);
  alpha_offset = ROUND_UP(y_offset + m*sizeof(opk_vector*), sizeof(double));
  rho_offset = alpha_offset + m*sizeof(double);
  size = rho_offset + m*sizeof(double);
  opt = (opk_vmlmb*)opk_allocate_object(finalize_vmlmb, size);
  if (opt == NULL) {
    return NULL;
  }
  opt->s           = ADDRESS(opk_vector*, opt,     s_offset);
  opt->y           = ADDRESS(opk_vector*, opt,     y_offset);
  opt->alpha       = ADDRESS(double,        opt, alpha_offset);
  opt->rho         = ADDRESS(double,        opt,   rho_offset);
  opt->m           = m;
  opt->gamma       = 1.0;
  opt->delta       = opts->delta;
  opt->epsilon     = opts->epsilon;
  opt->grtol       = opts->grtol;
  opt->gatol       = opts->gatol;
  opt->stpmin      = opts->stpmin;
  opt->stpmax      = opts->stpmax;
  opt->save_memory = opts->save_memory;
  if (box == NULL) {
    opt->method = OPK_LBFGS;
  } else if (opts->blmvm) {
    /* For the BLMVM method, the scratch vector is used to store the
       projected gradient. */
    opt->method = OPK_BLMVM;
    opt->gp = opk_vcreate(space);
    if (opt->gp == NULL) {
      goto error;
    }
  } else {
    opt->method = OPK_VMLMB;
  }

  /* Allocate work vectors.  If saving memory, x0 and g0 will be weak
     references to one of the saved vectors in the LBFGS operator. */
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
  if (lnsrch != NULL) {
    opt->lnsrch = OPK_HOLD_LNSRCH(lnsrch);
  } else {
    if (box != NULL) {
      opt->lnsrch = opk_lnsrch_new_backtrack(SFTOL, SAMIN);
    } else {
      opt->lnsrch = opk_lnsrch_new_csrch(SFTOL, SGTOL, SXTOL);
    }
    if (opt->lnsrch == NULL) {
      goto error;
    }
  }
  if (! opt->save_memory) {
    opt->x0 = opk_vcreate(space);
    if (opt->x0 == NULL) {
      goto error;
    }
    opt->g0 = opk_vcreate(space);
    if (opt->g0 == NULL) {
      goto error;
    }
  }
  opt->d = opk_vcreate(space);
  if (opt->d == NULL) {
    goto error;
  }
  if (box != NULL) {
    opt->box = OPK_HOLD_CONVEXSET(box);
    opt->w = opk_vcreate(space);
    if (opt->w == NULL) {
      goto error;
    }
  }

  /* Enforce calling opk_vmlmb_start and return the optimizer. */
  failure(opt, OPK_NOT_STARTED);
  return opt;

 error:
  OPK_DROP(opt);
  return NULL;
}

opk_task
opk_start_vmlmb(opk_vmlmb* opt, opk_vector* x)
{
  opt->iterations = 0;
  opt->evaluations = 0;
  opt->restarts = 0;
  opt->updates = 0;
  opt->mp = 0;
  if (opt->box != NULL) {
    opk_status status = opk_project_variables(x, x, opt->box);
    if (status != OPK_SUCCESS) {
      return failure(opt, status);
    }
  }
  return success(opt, OPK_TASK_COMPUTE_FG);
}

/* Define a few macros to make the code more readable. */
#define S(k)                           opt->s[k]
#define Y(k)                           opt->y[k]
#define ALPHA(k)                       opt->alpha[k]
#define RHO(k)                         opt->rho[k]
#define SLOT(j)                        slot(opt, j)
#define UPDATE(dst, alpha, x)          AXPBY(dst, 1, dst, alpha, x)
#define AXPBY(dst, alpha, x, beta, y)  opk_vaxpby(dst, alpha, x, beta, y)
#define COPY(dst, src)                 opk_vcopy(dst, src)
#define SCALE(x, alpha)                opk_vscale(x, alpha, x)
#define DOT(x, y)                      opk_vdot(x, y)
#define WDOT(x, y)                     opk_vdot3(opt->w, x, y)
#define NORM2(x)                       opk_vnorm2(x)
#define WNORM2(x)                      sqrt(WDOT(x, x))

/* Update L-BFGS approximation of the Hessian. */
static void
update(opk_vmlmb* opt,
       const opk_vector* x,
       const opk_vector* g)
{
  double sty, yty;
  opk_index k;

  k = SLOT(0);
  AXPBY(S(k), 1, x, -1, opt->x0);
  AXPBY(Y(k), 1, g, -1, opt->g0);
  if (opt->method != OPK_VMLMB) {
    /* Compute initial inverse Hessian approximation. */
    sty = DOT(S(k), Y(k));
    if (sty <= 0) {
      /* This pair will be skipped.  This may however indicate a problem, see
         Nocedal & Wright "Numerical Optimization", section 8.1, p. 201 (1999).
         FIXME: restart? */
      RHO(k) = 0;
    } else {
      /* Compute RHO(k) and GAMMA. */
      RHO(k) = 1/sty;
      yty = DOT(Y(k), Y(k));
      if (yty > 0) {
        opt->gamma = sty/yty;
      }
    }
  }
  ++opt->updates;
  if (opt->mp < opt->m) {
    ++opt->mp;
  }
}

/* Apply the L-BFGS Strang's two-loop recursion to compute a search direction.
   Returned value indicates whether the operation was successful otherwise the
   direction is just the gradient (steepest ascent). */
static opk_status
apply(opk_vmlmb* opt, const opk_vector* g)
{
  double sty, yty;
  opk_index j, k;
  opk_status result;

  COPY(opt->d, g);
  result = OPK_NOT_POSITIVE_DEFINITE;
  if (opt->method != OPK_VMLMB) {
    /* Apply the original L-BFGS Strang's two-loop recursion. */
    for (j = 1; j <= opt->mp; ++j) {
      k = SLOT(j);
      if (RHO(k) > 0) {
        ALPHA(k) = RHO(k)*DOT(opt->d, S(k));
        UPDATE(opt->d, -ALPHA(k), Y(k));
        result = OPK_SUCCESS;
      }
    }
    if (result == OPK_SUCCESS) {
      if (opt->gamma != 1) {
        /* Apply initial inverse Hessian approximation. */
        SCALE(opt->d, opt->gamma);
      }
      for (j = opt->mp; j >= 1; --j) {
        k = SLOT(j);
        if (RHO(k) > 0) {
          double beta = RHO(k)*DOT(opt->d, Y(k));
          UPDATE(opt->d, ALPHA(k) - beta, S(k));
        }
      }
    }
  } else {
    /* Apply L-BFGS Strang's two-loop recursion restricted to the subspace of
       free variables. */
    opt->gamma = 0;
    for (j = 1; j <= opt->mp; ++j) {
      k = SLOT(j);
      sty = WDOT(Y(k), S(k));
      if (sty <= 0) {
        RHO(k) = 0;
      } else {
        RHO(k) = 1/sty;
        ALPHA(k) = RHO(k)*WDOT(opt->d, S(k));
        UPDATE(opt->d, -ALPHA(k), Y(k));
        result = OPK_SUCCESS;
        if (opt->gamma == 0) {
          yty = WDOT(Y(k), Y(k));
          opt->gamma = sty/yty;
        }
      }
    }
    if (result == OPK_SUCCESS) {
      if (opt->gamma != 1) {
        SCALE(opt->d, opt->gamma);
      }
      for (j = opt->mp; j >= 1; --j) {
        k = SLOT(j);
        if (RHO(k) > 0) {
          double beta = RHO(k)*WDOT(opt->d, Y(k));
          UPDATE(opt->d, ALPHA(k) - beta, S(k));
        }
      }
      /* Enforce search direction to belong to the subset of the free
         variables. */
      opk_vproduct(opt->d, opt->w, opt->d);
    }
  }
  return result;
}

opk_task
opk_iterate_vmlmb(opk_vmlmb* opt, opk_vector* x,
                  double f, opk_vector* g)
{
  double dtg, gtest, stpmin, stpmax;
  opk_index k;
  opk_status status;
  opk_lnsrch_task lnsrch_task;
  opk_bool bounded;

  if (opt == NULL) {
    return OPK_TASK_ERROR;
  }
  bounded = (opt->box != NULL);

  switch (opt->task) {

  case OPK_TASK_COMPUTE_FG:

    /* Caller has computed the function value and the gradient at the current
       point. */
    ++opt->evaluations;

    if (opt->evaluations > 1) {
      /* A line search is in progress, check whether it has converged. */
      if (opk_lnsrch_use_deriv(opt->lnsrch)) {
        if (bounded) {
          /* Compute the directional derivative as the inner product between
             the effective step and the gradient divided by the step length. */
#if 0
          if (opt->tmp == NULL &&
              (opt->tmp = opk_vcreate(opt->vspace)) == NULL) {
            return failure(opt, OPK_INSUFFICIENT_MEMORY);
          }
          AXPBY(opt->tmp, 1, x, -1, opt->x0);
          dtg = DOT(opt->tmp, g)/opt->stp;
#else
          dtg = (DOT(x, g) - DOT(opt->x0, g))/opt->stp;
#endif
        } else {
          /* Compute the directional derivative. */
          dtg = -opk_vdot(opt->d, g);
        }
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
        /* An error may have occurred during the line search.  Figure out
           whether this error can be safely ignored. */
        status = opk_lnsrch_get_status(opt->lnsrch);
        if (lnsrch_task != OPK_LNSRCH_WARNING ||
            status != OPK_ROUNDING_ERRORS_PREVENT_PROGRESS) {
          return failure(opt, status);
        }
      }
      ++opt->iterations;
    }

    /* The current step is acceptable.  Check for global convergence. */
    if (bounded) {
      /* Determine the set of free variables. */
      status = opk_get_free_variables(opt->w, x, opt->box, g, OPK_ASCENT);
      if (status != OPK_SUCCESS) {
        return failure(opt, status);
      }
    }
    if (opt->method == OPK_VMLMB) {
      /* Compute the Euclidean norm of the projected gradient. */
      opt->gnorm = WNORM2(g);
    } else if (opt->method == OPK_BLMVM) {
      /* Compute the projected gradient and its norm. */
      opk_vproduct(opt->gp, opt->w, g);
      opt->gnorm = NORM2(opt->gp);
    } else {
      /* Compute the Euclidean norm of the gradient. */
      opt->gnorm = NORM2(g);
    }
    if (opt->evaluations == 1) {
      opt->ginit = opt->gnorm;
    }
    gtest = max3(0.0, opt->gatol, opt->grtol*opt->ginit);
    return success(opt, (opt->gnorm <= gtest
                         ? OPK_TASK_FINAL_X
                         : OPK_TASK_NEW_X));

  case OPK_TASK_NEW_X:
  case OPK_TASK_FINAL_X:

    /* Compute a new search direction. */
    if (opt->iterations >= 1) {
      /* Update L-BFGS approximation of the Hessian. */
      update(opt, x, (opt->method == OPK_BLMVM ? opt->gp : g));
    }
    status = apply(opt, g);
    if (status == OPK_SUCCESS) {
      /* The L-BFGS approximation produces a search direction D.  To warrant
         convergence, we have to check whether -D is a sufficient descent
         direction (that is to say that D is a sufficient ascent direction).
         As shown by Zoutendijk, this is true if cos(theta) = (D/|D|)'.(G/|G|)
         is larger or equal EPSILON > 0, where G is the gradient at X and D
         the, ascent for us, search direction. */
      if (bounded) {
        /* Project the search direction produced by the L-BFGS recursion. */
        status = opk_project_direction(opt->d, x, opt->box, opt->d, OPK_ASCENT);
        if (status != OPK_SUCCESS) {
          return failure(opt, status);
        }
      }
      dtg = -DOT(opt->d, g);
      if (opt->epsilon > 0 && -dtg < opt->epsilon*NORM2(opt->d)*opt->gnorm) {
        /* -D is not a sufficient descent direction.  Set the directional
           derivative to zero to force using the steepest descent direction. */
        dtg = 0.0;
      }
    } else {
      /* The L-BFGS approximation is not available (first iteration or just
         after a reset) or failed to produce a direction.  Set the directional
         derivative to zero to use the steepest descent direction. */
      dtg = 0.0;
    }

    /* Determine the initial step length. */
    if (dtg < 0) {
      /* A sufficient descent direction has been produced by L-BFGS recursion.
         An initial unit step will be used. */
      opt->stp = 1.0;
    } else {
      /* First iteration or L-BFGS recursion failed to produce a sufficient
         descent direction, use the (projected) gradient as a search
         direction. */
      if (opt->mp > 0) {
        /* L-BFGS recursion did not produce a sufficient descent direction. */
        ++opt->restarts;
        opt->mp = 0;
      }
      if (opt->method == OPK_VMLMB) {
        /* Use the projected gradient. */
        opk_vproduct(opt->d, opt->w, g);
      } else if (opt->method == OPK_BLMVM) {
        /* Use the projected gradient (which has already been computed and
         * stored in the scratch vector). */
        opk_vcopy(opt->d, opt->gp);
      } else {
        /* Use the gradient. */
        opk_vcopy(opt->d, g);
      }
      dtg = -opt->gnorm*opt->gnorm;
      if (f != 0) {
        opt->stp = 2*fabs(f/dtg);
      } else {
        /* Use a small step compared to X. */
        double dnorm = opt->gnorm;
        double xnorm = (bounded ? WNORM2(x) : NORM2(x));
        if (xnorm > 0) {
          opt->stp = opt->delta*xnorm/dnorm;
        } else {
          opt->stp = opt->delta/dnorm;
        }
      }
    }

    stpmin = opt->stp*opt->stpmin;
    stpmax = opt->stp*opt->stpmax;
    if (bounded) {
      /* Shortcut the step length. */
      double bsmin1, bsmin2, bsmax;
      status = opk_get_step_limits(&bsmin1, &bsmin2, &bsmax,
                                   x, opt->box, opt->d, OPK_ASCENT);
      if (bsmin1 < 0) {
        fprintf(stderr, "FIXME: SMIN1 =%g, SMIN2 =%g, SMAX =%g\n",
                bsmin1, bsmin2, bsmax);
      }
      if (status != OPK_SUCCESS) {
        return failure(opt, status);
      }
      if (bsmax <= 0) {
        return failure(opt, OPK_WOULD_BLOCK);
      }
      if (opt->stp > bsmax) {
        opt->stp = bsmax;
      }
      if (stpmax > bsmax) {
        stpmax = bsmax;
      }
      opt->bsmin = bsmin2;
    }

    /* Save current point. */
    if (opt->save_memory) {
      k = SLOT(0);
      opt->x0 = S(k); /* weak reference */
      opt->g0 = Y(k); /* weak reference */
      if (opt->mp == opt->m) {
        --opt->mp;
      }
    }
    COPY(opt->x0, x);
    COPY(opt->g0, (opt->method == OPK_BLMVM ? opt->gp : g));
    opt->f0 = f;

    /* Start the line search and break to take the first step along the line
       search. */
    if (opk_lnsrch_start(opt->lnsrch, f, dtg, opt->stp,
                         stpmin, stpmax) != OPK_LNSRCH_SEARCH) {
      return failure(opt, opk_lnsrch_get_status(opt->lnsrch));
    }
    break;

  default:

    /* There must be something wrong. */
    return opt->task;

  }

  /* Compute a new trial point along the line search. */
  opk_vaxpby(x, 1, opt->x0, -opt->stp, opt->d);
  if (bounded && opt->stp > opt->bsmin) {
    opk_status status = opk_project_variables(x, x, opt->box);
    if (status != OPK_SUCCESS) {
      return failure(opt, status);
    }
  }
  return success(opt, OPK_TASK_COMPUTE_FG);
}

opk_vmlmb_method
opk_get_vmlmb_method(const opk_vmlmb* opt)
{
  return opt->method;
}

const char*
opk_get_vmlmb_method_name(const opk_vmlmb* opt)
{
  switch (opt->method) {
  case OPK_LBFGS: return "VMLM/L-BFGS";
  case OPK_VMLMB: return "VMLMB";
  case OPK_BLMVM: return "BLMVM";
  default: return "*** unknown method ***";
  }
}

size_t
opk_get_vmlmb_description(char* buf, size_t size, const opk_vmlmb* opt)
{
  char str[80];

  switch (opt->method) {

  case OPK_LBFGS:
    sprintf(str, "variable metric method with %ld memorized step(s)",
            (long)opt->m);
    break;

  case OPK_VMLMB:
  case OPK_BLMVM:
    sprintf(str, "variable metric method with %ld memorized step(s) and bounds",
            (long)opt->m);
    break;

  default:
    strcat(str, "*** unknown method ***");
  }

  return opk_copy_string(buf, size, str);
}

opk_task
opk_get_vmlmb_task(const opk_vmlmb* opt)
{
  return opt->task;
}

opk_status
opk_get_vmlmb_status(const opk_vmlmb* opt)
{
  return opt->status;
}

void
opk__set_vmlmb_status(opk_vmlmb* opt, opk_status status)
{
  opt->status = status;
  if (status != OPK_SUCCESS) {
    opt->task = OPK_TASK_ERROR;
  }
}

opk_index
opk_get_vmlmb_iterations(const opk_vmlmb* opt)
{
  return opt->iterations;
}

opk_index
opk_get_vmlmb_evaluations(const opk_vmlmb* opt)
{
  return opt->evaluations;
}

opk_index
opk_get_vmlmb_restarts(const opk_vmlmb* opt)
{
  return opt->restarts;
}

double
opk_get_vmlmb_step(const opk_vmlmb* opt)
{
  return (opt == NULL ? -1.0 : opt->stp);
}

double
opk_get_vmlmb_gnorm(const opk_vmlmb* opt)
{
  return (opt == NULL ? -1.0 : opt->gnorm);
}

opk_index
opk_get_vmlmb_mp(const opk_vmlmb* opt)
{
  return (opt == NULL ? 0 : opt->mp);
}

opk_vector*
opk_get_vmlmb_s(const opk_vmlmb* opt, opk_index j)
{
  return (0 <= j && opt != NULL && j <= opt->mp ? opt->s[slot(opt, j)] : NULL);
}

opk_vector*
opk_get_vmlmb_y(const opk_vmlmb* opt, opk_index j)
{
  return (0 <= j && opt != NULL && j <= opt->mp ? opt->y[slot(opt, j)] : NULL);
}

void
opk_get_vmlmb_default_options(opk_vmlmb_options* opts)
{
  if (opts != NULL) {
    opts->delta = DELTA;
    opts->epsilon = EPSILON;
    opts->grtol = GRTOL;
    opts->gatol = GATOL;
    opts->stpmin = STPMIN;
    opts->stpmax = STPMAX;
    opts->mem = 5;
    opts->blmvm = FALSE;
    opts->save_memory = TRUE;
  }
}

opk_status
opk_check_vmlmb_options(const opk_vmlmb_options* opts)
{
  if (opts == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  if (non_finite(opts->gatol) || opts->gatol < 0 ||
      non_finite(opts->grtol) || opts->grtol < 0 ||
      non_finite(opts->delta) || opts->delta <= 0 ||
      non_finite(opts->epsilon) || opts->epsilon < 0 || opts->epsilon >= 1 ||
      non_finite(opts->stpmin) || opts->stpmin < 0 ||
      non_finite(opts->stpmax) || opts->stpmax <= opts->stpmin ||
      opts->mem < 1) {
    return OPK_INVALID_ARGUMENT;
  }
  return OPK_SUCCESS;
}
