/*
 * vmlmn.c --
 *
 * Limited memory variable metric method possibly with bound constraints for
 * OptimPack library.
 *
 * Three distinct methods are implemented: VMLM (also known as L-BFGS), VMLMN
 * and BLMVM.  All methods use a limited memory BFGS approximation of the
 * inverse Hessian which has been first described by Nocedal (1980) under the
 * name VMLM (for Variable Metric Limited Memory) and popularized under the
 * name L-BFGS (for Limited memory BFGS) by Lyu and Nocedal (1989).  VMLMN and
 * BLMVM implement bound constraints.  BLMVM, invented by Benson & Moré (2001),
 * uses the projected gradient to update the BFGS approximation.  VMLMN,
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
 * Copyright (C) 2002, 2015 Éric Thiébaut
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
static const opk_bfgs_scaling_t SCALING = OPK_SCALING_OREN_SPEDICATO;

struct _opk_vmlmn {
  opk_object_t base;       /**< Base type (must be the first member). */
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
  double stpmin;           /**< Relative mimimum step length. */
  double stpmax;           /**< Relative maximum step length. */
  double bsmin;            /**< Step size to the first encountered bound. */
  opk_vspace_t* vspace;    /**< Variable space. */
  opk_lnsrch_t* lnsrch;    /**< Line search method. */
  opk_vector_t* x0;        /**< Variables at the start of the line search. */
  opk_vector_t* g0;        /**< Gradient at x0. */
  opk_vector_t* d;         /**< Anti-search direction; a iterate is computed
                                as: x1 = x0 - stp*d with stp > 0. */
  opk_vector_t* w;         /**< Vector whose elements are set to 1 for free
                                variables and 0 otherwise. */
  opk_vector_t* tmp;       /**< Scratch work vector sued to store the effective
                                step size and the projected gradient (for the
                                BLMVM method). */
  opk_bound_t* xl;         /**< Lower bound (or `NULL`). */
  opk_bound_t* xu;         /**< Upper bound (or `NULL`). */
  int bounds;              /**< Type of bounds. */
  int method;              /**< The method to use/emulate. */
  unsigned int flags;      /**< Bitwise options. */
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
  return (a3 >= a2 ? MAX(a1, a3) : MAX(a1, a2));
}

static int
non_finite(double x)
{
  return (isnan(x) || isinf(x));
}

static void
finalize_vmlmn(opk_object_t* obj)
{
  opk_vmlmn_t* opt = (opk_vmlmn_t*)obj;
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
  OPK_DROP(opt->w);
  OPK_DROP(opt->tmp);
  OPK_DROP(opt->xl);
  OPK_DROP(opt->xu);
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
static opk_index_t
slot(const opk_vmlmn_t* opt, opk_index_t j)
{
  return (opt->updates - j)%opt->m;
}

static opk_task_t
success(opk_vmlmn_t* opt, opk_task_t task)
{
  opt->status = OPK_SUCCESS;
  opt->task = task;
  return task;
}

static opk_task_t
failure(opk_vmlmn_t* opt, opk_status_t status)
{
  opt->status = status;
  opt->task = OPK_TASK_ERROR;
  return opt->task;
}

opk_vmlmn_t*
opk_new_vmlmn_optimizer(opk_vspace_t* space,
                        opk_index_t m,
                        unsigned int flags,
                        opk_bound_t* xl,
                        opk_bound_t* xu,
                        opk_lnsrch_t* lnsrch)
{
  opk_vmlmn_t* opt;
  size_t s_offset, y_offset, beta_offset, rho_offset, size;
  opk_index_t k;
  int bounds;

  /* Check the input arguments for errors. */
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
  bounds = 0;
  if (xl != NULL) {
    if (xl->owner != space) {
      errno = EINVAL;
      return NULL;
    }
    if (xl->type == OPK_BOUND_SCALAR) {
      bounds += 1;
    } else if (xl->type == OPK_BOUND_VECTOR) {
      bounds += 2;
    } else {
      /* Discard unused bound. */
      xl = NULL;
    }
  }
  if (xu != NULL) {
    if (xu->owner != space) {
      errno = EINVAL;
      return NULL;
    }
    if (xu->type == OPK_BOUND_SCALAR) {
      bounds += 3;
    } else if (xu->type == OPK_BOUND_VECTOR) {
      bounds += 6;
    } else {
      /* Discard unused bound. */
      xu = NULL;
    }
  }

  /* Allocate enough memory for the workspace and its arrays. */
  s_offset = ROUND_UP(sizeof(opk_vmlmn_t), sizeof(opk_vector_t*));
  y_offset = s_offset + m*sizeof(opk_vector_t*);
  beta_offset = ROUND_UP(y_offset + m*sizeof(opk_vector_t*), sizeof(double));
  rho_offset = beta_offset + m*sizeof(double);
  size = rho_offset + m*sizeof(double);
  opt = (opk_vmlmn_t*)opk_allocate_object(finalize_vmlmn, size);
  if (opt == NULL) {
    return NULL;
  }
  opt->s      = ADDRESS(opk_vector_t*, opt,    s_offset);
  opt->y      = ADDRESS(opk_vector_t*, opt,    y_offset);
  opt->beta   = ADDRESS(double,        opt, beta_offset);
  opt->rho    = ADDRESS(double,        opt,  rho_offset);
  opt->m      = m;
  opt->gamma  = 1.0;
  opt->bounds = bounds;
  opt->flags  = flags;
  if (opt->bounds == 0) {
    opt->method = OPK_LBFGS;
  } else if ((flags & OPK_EMULATE_BLMVM) != 0) {
    /* For the BLMVM method, the scratch vector is used to store the
       projected gradient. */
    opt->method = OPK_BLMVM;
    opt->tmp = opk_vcreate(space);
    if (opt->tmp == NULL) {
      goto error;
    }
  } else {
    opt->method = OPK_VMLMN;
  }
  opk_set_vmlmn_options(opt, NULL);

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
    if (bounds != 0) {
      opt->lnsrch = opk_lnsrch_new_backtrack(SFTOL, SAMIN);
    } else {
      opt->lnsrch = opk_lnsrch_new_csrch(SFTOL, SGTOL, SXTOL);
    }
    if (opt->lnsrch == NULL) {
      goto error;
    }
  }
  if (xl != NULL) {
    opt->xl = OPK_HOLD_BOUND(xl);
  }
  if (xu != NULL) {
    opt->xu = OPK_HOLD_BOUND(xu);
  }
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
  if (opt->bounds != 0) {
    opt->w = opk_vcreate(space);
    if (opt->w == NULL) {
      goto error;
    }
  }

  /* Enforce calling opk_vmlmn_start and return the optimizer. */
  failure(opt, OPK_NOT_STARTED);
  return opt;

 error:
  OPK_DROP(opt);
  return NULL;
}

opk_task_t
opk_start_vmlmn(opk_vmlmn_t* opt, opk_vector_t* x)
{
  opt->iterations = 0;
  opt->evaluations = 0;
  opt->restarts = 0;
  opt->updates = 0;
  opt->mp = 0;
  if (opt->bounds != 0) {
    opk_status_t status = opk_box_project_variables(x, x, opt->xl, opt->xu);
    if (status != OPK_SUCCESS) {
      return failure(opt, status);
    }
  }
  return success(opt, OPK_TASK_COMPUTE_FG);
}

/* Define a few macros to make the code more readable. */
#define S(k)                           opt->s[k]
#define Y(k)                           opt->y[k]
#define BETA(k)                        opt->beta[k]
#define RHO(k)                         opt->rho[k]
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
update(opk_vmlmn_t* opt,
       const opk_vector_t* x,
       const opk_vector_t* g)
{
  double sty, yty;
  opk_index_t k;

  k = slot(opt, 0);
  AXPBY(S(k), 1, x, -1, opt->x0);
  AXPBY(Y(k), 1, g, -1, opt->g0);
  if (opt->method != OPK_VMLMN) {
    /* Compute initial inverse Hessian approximation. */
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
  }
  ++opt->updates;
  if (opt->mp < opt->m) {
    ++opt->mp;
  }
}

/* Apply the L-BFGS Strang's two-loop recursion to compute a search
   direction. */
static opk_status_t
apply(opk_vmlmn_t* opt, const opk_vector_t* g)
{
  double sty, yty;
  opk_index_t j, k;

  if (opt->mp < 1) {
    /* Will use the steepest descent direction. */
    return OPK_NOT_POSITIVE_DEFINITE;
  }
  COPY(opt->d, g);
  if (opt->method != OPK_VMLMN) {
    /* Apply the original L-BFGS Strang's two-loop recursion. */
    for (j = 1; j <= opt->mp; ++j) {
      k = slot(opt, j);
      if (RHO(k) > 0) {
        BETA(k) = RHO(k)*DOT(opt->d, S(k));
        UPDATE(opt->d, -BETA(k), Y(k));
      }
    }
    if (opt->gamma != 1) {
      /* Apply initial inverse Hessian approximation. */
      SCALE(opt->d, opt->gamma);
    }
    for (j = opt->mp; j >= 1; --j) {
      k = slot(opt, j);
      if (RHO(k) > 0) {
        UPDATE(opt->d, BETA(k) - RHO(k)*DOT(opt->d, Y(k)), S(k));
      }
    }
  } else {
    /* Apply L-BFGS Strang's two-loop recursion restricted to the subspace of
       free variables. */
    opt->gamma = 0;
    for (j = 1; j <= opt->mp; ++j) {
      k = slot(opt, j);
      sty = WDOT(Y(k), S(k));
      if (sty <= 0) {
        RHO(k) = 0;
        continue;
      }
      RHO(k) = 1/sty;
      BETA(k) = RHO(k)*WDOT(opt->d, S(k));
      UPDATE(opt->d, -BETA(k), Y(k));
      if (opt->gamma == 0) {
        yty = WDOT(Y(k), Y(k));
        if (yty > 0) {
          opt->gamma = sty/yty;
        }
      }
    }
    if (opt->gamma != 1) {
      if (opt->gamma <= 0) {
        /* Force using the steepest descent direction. */
        return OPK_NOT_POSITIVE_DEFINITE;
      }
      SCALE(opt->d, opt->gamma);
    }
    for (j = opt->mp; j >= 1; --j) {
      k = slot(opt, j);
      if (RHO(k) > 0) {
        UPDATE(opt->d, BETA(k) - RHO(k)*WDOT(opt->d, Y(k)), S(k));
      }
    }
  }

  if (opt->bounds != 0) {
    /* Enforce search direction to belong to the subset of the free
       variables. */
    opk_vproduct(opt->d, opt->w, opt->d);
  }

  return OPK_SUCCESS;
}

opk_task_t
opk_iterate_vmlmn(opk_vmlmn_t* opt, opk_vector_t* x,
                  double f, opk_vector_t* g)
{
  double dtg;
  opk_index_t k;
  opk_status_t status;
  opk_lnsrch_task_t lnsrch_task;
  opk_bool_t bounded, final;

  bounded = (opt->bounds != 0);

  switch (opt->task) {

  case OPK_TASK_COMPUTE_FG:

    /* Caller has computed the function value and the gradient at the current
       point. */
    ++opt->evaluations;

    if (opt->evaluations > 1) {
      /* A line search is in progress, check whether it has converged. */
      if (opk_lnsrch_use_deriv(opt->lnsrch)) {
        if (bounded) {
          /* Compute effective step and directional derivative. */
          if (opt->tmp == NULL &&
              (opt->tmp = opk_vcreate(opt->vspace)) == NULL) {
            return failure(opt, OPK_INSUFFICIENT_MEMORY);
          }
          AXPBY(opt->tmp, 1, x, -1, opt->x0);
          dtg = DOT(opt->tmp, g)/opt->stp;
        } else {
          /* Compute directional derivative. */
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
        status = opk_lnsrch_get_status(opt->lnsrch);
        if (lnsrch_task != OPK_LNSRCH_WARNING ||
            status != OPK_ROUNDING_ERRORS_PREVENT_PROGRESS) {
          return failure(opt, status);
        }
      }
      ++opt->iterations;
    }

    if (bounded) {
      /* Determine the set of free variables. */
      status = opk_box_get_free_variables(opt->w, x, opt->xl, opt->xu,
                                          g, OPK_ASCENT);
      if (status != OPK_SUCCESS) {
        return failure(opt, status);
      }
    }

    /* Check for global convergence. */
    if (opt->method == OPK_VMLMN) {
      /* Compute the Euclidean norm of the projected gradient. */
      opt->gnorm = WNORM2(g);
    } else if (opt->method == OPK_BLMVM) {
      /* Compute the projected gradient and its norm. */
      opk_vproduct(opt->tmp, opt->w, g);
      opt->gnorm = NORM2(opt->tmp);
    } else {
      /* Compute the Euclidean norm of the gradient. */
      opt->gnorm = NORM2(g);
    }
    if (opt->evaluations == 1) {
      opt->ginit = opt->gnorm;
    }
    final = (opt->gnorm <= max3(0.0, opt->gatol, opt->grtol*opt->ginit));
    return success(opt, (final ? OPK_TASK_FINAL_X : OPK_TASK_NEW_X));

  case OPK_TASK_NEW_X:
  case OPK_TASK_FINAL_X:

    /* Compute a new search direction. */
    if (opt->iterations >= 1) {
      /* Update L-BFGS approximation of the Hessian. */
      update(opt, x, (opt->method == OPK_BLMVM ? opt->tmp : g));
    }
    if (apply(opt, g) == OPK_SUCCESS) {
      /* We take care of checking whether -D is a sufficient descent direction
         (that is to say that D is a sufficient ascent direction).  As shown by
         Zoutendijk, this is true if cos(theta) = (D/|D|)'.(G/|G|) is larger or
         equal EPSILON > 0, where G is the gradient at X and D the (ascent for
         us) search direction. */
      dtg = -DOT(opt->d, g);
      if (opt->epsilon > 0 &&
          dtg > -opt->epsilon*NORM2(opt->d)*opt->gnorm) {
        /* We do not have a sufficient descent direction.  Set the directional
           derivative to zero to force using the steepest descent direction. */
        dtg = 0.0;
      }
    } else {
      /* The L-BFGS approximation is unset (first iteration) or failed to
         produce a direction.  Set the directional derivative to zero to use
         the steepest descent direction. */
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
      if (opt->method == OPK_VMLMN) {
        /* Use the projected gradient. */
        opk_vproduct(opt->d, opt->w, g);
      } else if (opt->method == OPK_BLMVM) {
        /* Use the projected gradient (which has aready been computed and
         * stored in the scratch vector). */
        opk_vcopy(opt->d, opt->tmp);
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

    if (bounded) {
      /* Shortcut the step length. */
      double bsmin, bsmax, wolfe;
      status = opk_box_get_step_limits(&bsmin, &wolfe, &bsmax,
                                       x, opt->xl, opt->xu,
                                       opt->d, OPK_ASCENT);
      if (status != OPK_SUCCESS) {
        return failure(opt, status);
      }
      if (bsmax <= 0) {
        return failure(opt, OPK_WOULD_BLOCK);
      }
      if (opt->stp > bsmax) {
        opt->stp = bsmax;
      }
      opt->bsmin = bsmin;
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
    COPY(opt->x0, x);
    COPY(opt->g0, (opt->method == OPK_BLMVM ? opt->tmp : g));
    opt->f0 = f;

    /* Start the line search and break to take the first step along the line
       search. */
    if (opk_lnsrch_start(opt->lnsrch, f, dtg, opt->stp,
                         opt->stp*opt->stpmin,
                         opt->stp*opt->stpmax) != OPK_LNSRCH_SEARCH) {
      return failure(opt, opk_lnsrch_get_status(opt->lnsrch));
    }
    break;

  default:

    /* There must be something wrong. */
    return opt->task;

  }

  /* Compute a new trial point along the line search. */
  opk_vaxpby(x, 1, opt->x0, -opt->stp, opt->d);
  if (opt->bounds != 0 && opt->stp > opt->bsmin) {
    opk_status_t status = opk_box_project_variables(x, x, opt->xl, opt->xu);
    if (status != OPK_SUCCESS) {
      return failure(opt, status);
    }
  }
  return success(opt, OPK_TASK_COMPUTE_FG);
}

opk_task_t
opk_get_vmlmn_task(opk_vmlmn_t* opt)
{
  return opt->task;
}

opk_status_t
opk_get_vmlmn_status(opk_vmlmn_t* opt)
{
  return opt->status;
}

opk_vmlmn_method_t
opk_get_vmlmn_method(opk_vmlmn_t* opt)
{
  return opt->method;
}

const char*
opk_get_vmlmn_method_name(opk_vmlmn_t* opt)
{
  switch (opt->method) {
  case OPK_LBFGS: return "VMLM/L-BFGS";
  case OPK_VMLMN: return "VMLMB";
  case OPK_BLMVM: return "BLMVM";
  default: return "*** unknown method ***";
  }
}

opk_index_t
opk_get_vmlmn_iterations(opk_vmlmn_t* opt)
{
  return opt->iterations;
}

opk_index_t
opk_get_vmlmn_evaluations(opk_vmlmn_t* opt)
{
  return opt->evaluations;
}

opk_index_t
opk_get_vmlmn_restarts(opk_vmlmn_t* opt)
{
  return opt->restarts;
}

double
opk_get_vmlmn_step(opk_vmlmn_t* opt)
{
  return (opt == NULL ? -1.0 : opt->stp);
}

opk_vector_t*
opk_get_vmlmn_s(opk_vmlmn_t* opt, opk_index_t k)
{
  return (0 <= k && k <= opt->mp ? opt->s[slot(opt, k)] : NULL);
}

opk_vector_t*
opk_get_vmlmn_y(opk_vmlmn_t* opt, opk_index_t k)
{
  return (0 <= k && k <= opt->mp ? opt->y[slot(opt, k)] : NULL);
}

opk_status_t
opk_get_vmlmn_options(opk_vmlmn_options_t* dst, const opk_vmlmn_t* src)
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
opk_set_vmlmn_options(opk_vmlmn_t* dst, const opk_vmlmn_options_t* src)
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
