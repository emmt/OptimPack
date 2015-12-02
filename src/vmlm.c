/*
 * vmlm.c --
 *
 * Limited memory variable metric methods for OptimPack library.
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of OptimPack (https://github.com/emmt/OptimPack).
 *
 * Copyright (C) 2014, 2015 Éric Thiébaut
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

#define TRUE   OPK_TRUE
#define FALSE  OPK_FALSE

#define MAX(a,b)  OPK_MAX(a,b)
#define MIN(a,b)  OPK_MIN(a,b)

#define ROUND_UP(a,b)   OPK_ROUND_UP(a,b)

/*---------------------------------------------------------------------------*/
/* LIMITED MEMORY BFGS OPERATOR */

struct _opk_lbfgs_operator {
  opk_operator_t base; /**< Base type (must be the first member). */
  double gamma;        /**< Scale factor to approximate inverse Hessian. */
  opk_vector_t** s;    /**< Storage for variable differences. */
  opk_vector_t** y;    /**< Storage for gradient differences. */
  opk_vector_t* tmp;   /**< Scratch vector, needed if there is a
                            preconditioner. */
  opk_operator_t* H0;  /**< Preconditioner or NULL. */
  double* alpha;       /**< Workspace to save <s[j],d>/<s[j],y[j]> */
  double* rho;         /**< Workspace to save 1/<s[j],y[j]> */
  opk_index_t m;       /**< Maximum number of memorized steps. */
  opk_index_t mp;      /**< Actual number of memorized steps (0 <= mp <= m). */
  opk_index_t updates; /**< Number of updates since start. */
  int scaling;
};

static opk_index_t
slot(opk_lbfgs_operator_t* op, opk_index_t j)
{
  return (op->updates - j)%op->m;
}

static void
finalize_lbfgs_operator(opk_operator_t* self)
{
  opk_lbfgs_operator_t* op = (opk_lbfgs_operator_t*)self;
  opk_index_t k;

  if (op->s != NULL) {
    for (k = 0; k < op->m; ++k) {
      OPK_DROP(op->s[k]);
    }
  }
  if (op->y != NULL) {
    for (k = 0; k < op->m; ++k) {
      OPK_DROP(op->y[k]);
    }
  }
  OPK_DROP(op->tmp);
  OPK_DROP(op->H0);
}

static opk_status_t
apply_lbfgs_operator(opk_operator_t* self, opk_vector_t* dst,
                     const opk_vector_t* src)
{
  opk_lbfgs_operator_t* op = (opk_lbfgs_operator_t*)self;
  opk_vector_t* tmp;
  opk_index_t j, k;
  opk_status_t status;

  /* Initialize work vectors. */
  if (op->H0 == NULL) {
    /* No preconditioner, no needs for a scratch vector. */
    tmp = dst;
  } else {
    /* With a preconditioner, a scratch vector is needed. */
    if (op->tmp == NULL) {
      op->tmp = opk_vcreate(op->base.outspace);
      if (op->tmp == NULL) {
        return OPK_INSUFFICIENT_MEMORY;
      }
    }
    tmp = op->tmp;
  }

  /* First loop of the recursion (from the newest saved pair to the oldest
     one). */
  opk_vcopy(tmp, src);
  for (j = 1; j <= op->mp; ++j) {
    k = slot(op, j);
    if (op->rho[k] > 0.0) {
      op->alpha[k] = op->rho[k]*opk_vdot(tmp, op->s[k]);
      opk_vaxpby(tmp, 1.0, tmp, -op->alpha[k], op->y[k]);
    } else {
      op->alpha[k] = 0.0;
    }
  }

  /* Apply approximation of inverse Hessian. */
  if (op->H0 != NULL) {
    status = opk_apply_direct(op->H0, dst, tmp);
    if (status != OPK_SUCCESS) {
      return status;
    }
  }
  if (op->gamma != 1.0) {
    opk_vscale(dst, op->gamma, dst);
  }

  /* Second loop of the recursion (from the oldest saved pair to the newest
     one). */
  for (j = op->mp; j >= 1; --j) {
    k = slot(op, j);
    if (op->rho[k] > 0.0) {
      double beta = op->rho[k]*opk_vdot(dst, op->y[k]);
      opk_vaxpby(dst, 1.0, dst, op->alpha[k] - beta, op->s[k]);
    }
  }
  return OPK_SUCCESS;
}

static opk_operator_operations_t lbfgs_operations = {
  finalize_lbfgs_operator,
  apply_lbfgs_operator,
  apply_lbfgs_operator,
  NULL,
};

opk_lbfgs_operator_t*
opk_new_lbfgs_operator(opk_vspace_t* vspace, opk_index_t m,
                       int scaling)
{
  opk_lbfgs_operator_t* op;
  size_t s_offset, y_offset, alpha_offset, rho_offset, size;
  opk_index_t k;

  /* Check arguments. */
  if (vspace == NULL) {
    errno = EFAULT;
    return NULL;
  }
  if (vspace->size < 1 || m < 1) {
    errno = EINVAL;
    return NULL;
  }
  if (m > vspace->size) {
    m = vspace->size;
  }

  /* Allocate enough memory for the workspace and its arrays. */
  s_offset = ROUND_UP(sizeof(opk_lbfgs_operator_t), sizeof(opk_vector_t*));
  y_offset = s_offset + m*sizeof(opk_vector_t*);
  alpha_offset = ROUND_UP(y_offset + m*sizeof(opk_vector_t*), sizeof(double));
  rho_offset = alpha_offset + m*sizeof(double);
  size = rho_offset + m*sizeof(double);
  op = (opk_lbfgs_operator_t*)opk_allocate_operator(&lbfgs_operations,
                                                    vspace, vspace, size);
  if (op == NULL) {
    return NULL;
  }
  op->s = (opk_vector_t**)(((char*)op) + s_offset);
  op->y = (opk_vector_t**)(((char*)op) + y_offset);
  op->alpha = (double*)(((char*)op) + alpha_offset);
  op->rho = (double*)(((char*)op) + rho_offset);
  op->gamma = 1.0;
  op->m = m;
  op->mp = 0;
  op->updates = 0;
  op->scaling = scaling;
  for (k = 0; k < m; ++k) {
    op->s[k] = opk_vcreate(vspace);
    if (op->s[k] == NULL) {
      goto error;
    }
    op->y[k] = opk_vcreate(vspace);
    if (op->y[k] == NULL) {
      goto error;
    }
  }
  return op;
 error:
  OPK_DROP(op);
  return NULL;
}

void
opk_reset_lbfgs_operator(opk_lbfgs_operator_t* op)
{
  op->mp = 0;
  op->gamma = 1.0;
}

opk_vector_t*
opk_get_lbfgs_s(opk_lbfgs_operator_t* op, opk_index_t k)
{
  return (0 <= k && k <= op->mp ? op->s[slot(op, k)] : NULL);
}

opk_vector_t*
opk_get_lbfgs_y(opk_lbfgs_operator_t* op, opk_index_t k)
{
  return (0 <= k && k <= op->mp ? op->y[slot(op, k)] : NULL);
}

void
opk_set_lbfgs_operator_preconditioner(opk_lbfgs_operator_t* op,
                                      opk_operator_t* H0)
{
  if (op->H0 != H0) {
    opk_operator_t* old = op->H0;
    op->H0 = OPK_HOLD_OPERATOR(H0);
    OPK_DROP(old);
  }
}

void
opk_update_lbfgs_operator(opk_lbfgs_operator_t* op,
                          const opk_vector_t* x,
                          const opk_vector_t* x0,
                          const opk_vector_t* g,
                          const opk_vector_t* g0)
{
  double sty;
  opk_index_t k;

  /* Store the variables and gradient differences in the slot just after the
     last update. */
  k = slot(op, 0);
  opk_vaxpby(op->s[k], 1, x, -1, x0);
  opk_vaxpby(op->y[k], 1, g, -1, g0);

  /* Compute RHO[k] and GAMMA.  If the update formula for GAMMA does not yield
     a strictly positive value, the strategy is to keep the previous value. */
  sty = opk_vdot(op->s[k], op->y[k]);
  if (sty <= 0.0) {
    /* This pair will be skipped. */
    op->rho[k] = 0.0;
    if (op->mp == op->m) {
      --op->mp;
    }
  } else {
    /* Compute RHO[k] and GAMMA. */
    op->rho[k] = 1.0/sty;
    if (op->scaling == OPK_SCALING_OREN_SPEDICATO) {
      double ynorm = opk_vnorm2(op->y[k]);
      op->gamma = (sty/ynorm)/ynorm;
    } else if (op->scaling == OPK_SCALING_BARZILAI_BORWEIN) {
      double snorm = opk_vnorm2(op->s[k]);
      op->gamma = (snorm/sty)*snorm;
    }
    /* Update the mark and the number of saved pairs. */
    ++op->updates;
    if (op->mp < op->m) {
      ++op->mp;
    }
  }
}

/*---------------------------------------------------------------------------*/
/* VARIABLE METRIC OPTIMIZER */

/* Default parameters for the line search.  These values are from original
   LBFGS algorithm by Jorge Nocedal but other values may be more suitable. */
static const double STPMIN = 1E-20;
static const double STPMAX = 1E+20;
static const double SFTOL = 1E-4;
static const double SGTOL = 0.9;
static const double SXTOL = DBL_EPSILON;

/* Default parameters for the global convergence. */
static const double GRTOL = 1E-6;
static const double GATOL = 0.0;

/* Other default parameters. */
static const double DELTA   = 5e-2;
static const double EPSILON = 1e-2;
static const int    SCALING = OPK_SCALING_OREN_SPEDICATO;

struct _opk_vmlm {
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
  double g0norm;           /**< Euclidean norm of g0 the gradient at x0. */
  double gnorm;            /**< Euclidean norm of the gradient at the last
                                tried x. */
  double stp;              /**< Current step length. */
  double stpmin;           /**< Relative mimimum step length. */
  double stpmax;           /**< Relative maximum step length. */
  opk_vspace_t* vspace;    /**< Variable space. */
  opk_lnsrch_t* lnsrch;    /**< Line search method. */
  opk_lbfgs_operator_t* H; /**< L-BFGS approximation of inverse Hessian. */
  opk_vector_t* x0;        /**< Variables at the start of the line search. */
  opk_vector_t* g0;        /**< Gradient at x0. */
  opk_vector_t* d;         /**< Anti-search direction; a iterate is computed
                                as: x1 = x0 - stp*d with stp > 0. */
  opk_index_t evaluations; /**< Number of functions (and gradients)
                                evaluations. */
  opk_index_t iterations;  /**< Number of iterations (successful steps
                                taken). */
  opk_index_t restarts;    /**< Number of restarts. */
  opk_status_t reason;     /**< Last error. */
  opk_task_t task;         /**< Pending task. */
  opk_bool_t save_memory;  /**< To save space, the variable and gradient at the
                                start of a line search are weak references to
                                the (s,y) pair of vectors of the LBFGS operator
                                just after the mark. */
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
finalize_vmlm(opk_object_t* obj)
{
  opk_vmlm_t* opt = (opk_vmlm_t*)obj;
  OPK_DROP(opt->vspace);
  OPK_DROP(opt->lnsrch);
  OPK_DROP(opt->H);
  if (! opt->save_memory) {
    OPK_DROP(opt->x0);
    OPK_DROP(opt->g0);
  }
  OPK_DROP(opt->d);
}

static opk_task_t
success(opk_vmlm_t* opt, opk_task_t task)
{
  opt->reason = OPK_SUCCESS;
  opt->task = task;
  return task;
}

static opk_task_t
failure(opk_vmlm_t* opt, opk_status_t status)
{
  opt->reason = status;
  opt->task = OPK_TASK_ERROR;
  return opt->task;
}

opk_vmlm_t*
opk_new_vmlm_optimizer_with_line_search(opk_vspace_t* vspace,
                                        opk_index_t m,
                                        opk_lnsrch_t* lnsrch)
{
  opk_vmlm_t* opt;

  /* Check the input arguments for errors. */
  if (vspace == NULL || lnsrch == NULL) {
    errno = EFAULT;
    return NULL;
  }
  if (vspace->size < 1 || m < 1) {
    errno = EINVAL;
    return NULL;
  }

  /* Allocate and instanciate the workspace (not the part which is done by
     opk_start_vmlm). */
  opt = (opk_vmlm_t*)opk_allocate_object(finalize_vmlm, sizeof(opk_vmlm_t));
  if (opt == NULL) {
    return NULL;
  }
  opt->vspace = OPK_HOLD_VSPACE(vspace);
  opt->lnsrch = OPK_HOLD_LNSRCH(lnsrch);
  opt->H = opk_new_lbfgs_operator(vspace, m, SCALING);
  if (opt->H == NULL) {
    goto error;
  }
  opt->save_memory = TRUE;
  opt->grtol = GRTOL;
  opt->gatol = GATOL;
  opt->stpmin = STPMIN;
  opt->stpmax = STPMAX;
  opt->delta = DELTA;
  opt->epsilon = EPSILON;

  /* Allocate work vectors.  If saving memory, x0 and g0 will be weak references
     to one of the saved vectors in the LBFGS operator. */
  if (! opt->save_memory) {
    opt->x0 = opk_vcreate(opt->vspace);
    if (opt->x0 == NULL) {
      goto error;
    }
    opt->g0 = opk_vcreate(vspace);
    if (opt->g0 == NULL) {
      goto error;
    }
  }
  opt->d = opk_vcreate(vspace);
  if (opt->d == NULL) {
    goto error;
  }
  opk_start_vmlm(opt);
  return opt;

 error:
  OPK_DROP(opt);
  return NULL;
}

opk_vmlm_t*
opk_new_vmlm_optimizer(opk_vspace_t* vspace, opk_index_t m)
{
  opk_lnsrch_t* lnsrch;
  opk_vmlm_t* opt;

  lnsrch = opk_lnsrch_new_csrch(SFTOL, SGTOL, SXTOL);
  if (lnsrch == NULL) {
    return NULL;
  }
  opt = opk_new_vmlm_optimizer_with_line_search(vspace, m, lnsrch);
  OPK_DROP(lnsrch); /* the line search is now owned by the optimizer */
  return opt;
}

opk_task_t
opk_start_vmlm(opk_vmlm_t* opt)
{
  opt->iterations = 0;
  opt->evaluations = 0;
  opt->restarts = 0;
  opt->H->updates = 0;
  opk_reset_lbfgs_operator(opt->H);
  return success(opt, OPK_TASK_COMPUTE_FG);
}

opk_task_t
opk_iterate_vmlm(opk_vmlm_t* opt, opk_vector_t* x,
                 double f, opk_vector_t* g)
{
  double gtest, dtg;
  opk_status_t status;
  opk_lnsrch_task_t lnsrch_task;

  switch (opt->task) {

  case OPK_TASK_COMPUTE_FG:

    /* Caller has computed the function value and the gradient at the current
       point. */
    ++opt->evaluations;
    if (opt->evaluations > 1) {
      /* A line search is in progress, check whether it has converged. */
      if (opk_lnsrch_use_deriv(opt->lnsrch)) {
        /* Compute effective step and directional derivative. */
        dtg = -opk_vdot(opt->d, g);
      } else {
        /* Line search does not need directional derivative. */
        dtg = 0;
      }
      lnsrch_task = opk_lnsrch_iterate(opt->lnsrch, &opt->stp, f, dtg);
      if (lnsrch_task == OPK_LNSRCH_SEARCH) {
        /* Line search has not yet converged. */
        goto new_try;
      }
      if (lnsrch_task != OPK_LNSRCH_CONVERGENCE) {
        status = opk_lnsrch_get_status(opt->lnsrch);
        if (status != OPK_ROUNDING_ERRORS_PREVENT_PROGRESS) {
          return failure(opt, status);
        }
      }
      ++opt->iterations;
    }

    /* The current step is acceptable. Check for global convergence. */
    opt->gnorm = opk_vnorm2(g);
    if (opt->evaluations == 1) {
      opt->ginit = opt->gnorm;
    }
    gtest = max3(0.0, opt->gatol, opt->grtol*opt->ginit);
    return success(opt, (opt->gnorm <= gtest
                         ? OPK_TASK_FINAL_X
                         : OPK_TASK_NEW_X));

  case OPK_TASK_NEW_X:
  case OPK_TASK_FINAL_X:

    if (opt->evaluations > 1) {
      /* Update the LBFGS matrix. */
      opk_update_lbfgs_operator(opt->H, x, opt->x0, g, opt->g0);
    }

    /* Compute a search direction.  We take care of checking whether -D is a
       sufficient descent direction (here D is a sufficient ascent direction).
       As shown by Zoutendijk, this is true if cos(theta) = (D/|D|)'.(G/|G|) is
       larger or equal EPSILON > 0, where G is the gradient at X and D the
       descent direction. */

    while (TRUE) {
      opk_apply_direct((opk_operator_t*)opt->H, opt->d, g);
      dtg = -opk_vdot(opt->d, g);
      if (opt->epsilon > 0 &&
          dtg > -opt->epsilon*opk_vnorm2(opt->d)*opt->gnorm) {
        /* Set DTG to zero to indicate that we do not have a sufficient
           descent direction. */
        dtg = 0;
      }
      if (dtg < 0) {
        /* Accept D as a sufficient ascent direction. */
        break;
      }
      if (opt->H->mp < 1) {
        /* Initial iteration or recursion has just been restarted.  This means
           that the initial inverse Hessian approximation is not positive
           definite. */
        return failure(opt, OPK_BAD_PRECONDITIONER);
      }
      /* Reset the LBFGS recursion and loop to use H0 to compute an initial
         search direction. */
      opk_reset_lbfgs_operator(opt->H);
      ++opt->restarts;
    }

    /* Save current variables X0, gradient G0 and function value F0. */
    if (opt->save_memory) {
      /* Use the slot just after the mark to store X0 and G0.  Note that this
         is a weak reference: we do not "hold" the vectors. */
      opt->x0 = opk_get_lbfgs_s(opt->H, 0);
      opt->g0 = opk_get_lbfgs_y(opt->H, 0);
      if (opt->H->mp == opt->H->m) {
        --opt->H->mp;
      }
    }
    opk_vcopy(opt->x0, x);
    opk_vcopy(opt->g0, g);
    opt->f0 = f;
    opt->g0norm = opt->gnorm;

    /* Estimate the length of the first step. */
    if (opt->H->mp >= 1 || opt->H->scaling == OPK_SCALING_NONE) {
      opt->stp = 1.0;
    } else {
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

    /* Start the line search and take the first step along the search
       direction. */
    if (opk_lnsrch_start(opt->lnsrch, f, dtg, opt->stp,
                         opt->stpmin*opt->stp,
                         opt->stpmax*opt->stp) != OPK_LNSRCH_SEARCH) {
      return failure(opt, opk_lnsrch_get_status(opt->lnsrch));
    }
    goto new_try;

  default:

    /* There must be something wrong. */
    return opt->task;

  }

 new_try:
  opk_vaxpby(x, 1, opt->x0, -opt->stp, opt->d);
  return success(opt, OPK_TASK_COMPUTE_FG);
}

opk_task_t
opk_get_vmlm_task(opk_vmlm_t* opt)
{
  return opt->task;
}

opk_status_t
opk_get_vmlm_reason(opk_vmlm_t* opt)
{
  return opt->reason;
}

opk_index_t
opk_get_vmlm_iterations(opk_vmlm_t* opt)
{
  return opt->iterations;
}

opk_index_t
opk_get_vmlm_evaluations(opk_vmlm_t* opt)
{
  return opt->evaluations;
}

opk_index_t
opk_get_vmlm_restarts(opk_vmlm_t* opt)
{
  return opt->restarts;
}

int
opk_get_vmlm_scaling(opk_vmlm_t* opt)
{
  return (opt == NULL ? SCALING : opt->H->scaling);
}

opk_status_t
opk_set_vmlm_scaling(opk_vmlm_t* opt, int scaling)
{
  if (opt == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  switch (scaling) {
  case OPK_SCALING_NONE:
  case OPK_SCALING_OREN_SPEDICATO:
  case OPK_SCALING_BARZILAI_BORWEIN:
    opt->H->scaling = scaling;
    return OPK_SUCCESS;
  default:
    return OPK_INVALID_ARGUMENT;
  }
}

double
opk_get_vmlm_gatol(opk_vmlm_t* opt)
{
  return (opt == NULL ? GATOL : opt->gatol);
}

opk_status_t
opk_set_vmlm_gatol(opk_vmlm_t* opt, double gatol)
{
  if (opt == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  if (non_finite(gatol) || gatol < 0) {
    return OPK_INVALID_ARGUMENT;
  }
  opt->gatol = gatol;
  return OPK_SUCCESS;
}

double
opk_get_vmlm_grtol(opk_vmlm_t* opt)
{
  return (opt == NULL ? GRTOL : opt->grtol);
}

opk_status_t
opk_set_vmlm_grtol(opk_vmlm_t* opt, double grtol)
{
  if (opt == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  if (non_finite(grtol) || grtol < 0) {
    return OPK_INVALID_ARGUMENT;
  }
  opt->grtol = grtol;
  return OPK_SUCCESS;
}

double
opk_get_vmlm_step(opk_vmlm_t* opt)
{
  return (opt == NULL ? -1.0 : opt->stp);
}

double
opk_get_vmlm_stpmin(opk_vmlm_t* opt)
{
  return (opt == NULL ? STPMIN : opt->stpmin);
}

double
opk_get_vmlm_stpmax(opk_vmlm_t* opt)
{
  return (opt == NULL ? STPMAX : opt->stpmax);
}

opk_status_t
opk_set_vmlm_stpmin_and_stpmax(opk_vmlm_t* opt,
                               double stpmin, double stpmax)
{
  if (opt == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  if (non_finite(stpmin) || non_finite(stpmax) ||
      stpmin < 0 || stpmax <= stpmin) {
    return OPK_INVALID_ARGUMENT;
  }
  opt->stpmin = stpmin;
  opt->stpmax = stpmax;
  return OPK_SUCCESS;
}
