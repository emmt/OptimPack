/*
 * vmlm.c --
 *
 * Limited memory variable metric methods for OptimPack library.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (c) 2014 Éric Thiébaut
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
  double epsilon;      /**< Threshold to skip update. */
  double gamma;        /**< Scale factor to approximate inverse Hessian. */
  opk_vector_t** s;    /**< Storage for variable differences. */
  opk_vector_t** y;    /**< Storage for gradient differences. */
  opk_vector_t* tmp;   /**< Scratch vector, needed if there is a
                            preconditioner. */
  opk_operator_t* H0;  /**< Preconditioner or NULL. */
  double* beta;        /**< Workspace to save 1/<s[j],y[j]> */
  double* rho;         /**< Workspace to save 1/<s[j],y[j]> */
  opk_index_t m;       /**< Maximum number of memorized steps. */
  opk_index_t mp;      /**< Actual number of memorized steps (0 <= mp <= m). */
  opk_index_t mark;    /**< Index of oldest saved step. */
  int scaling;
};

static opk_index_t
slot(opk_lbfgs_operator_t* op, opk_index_t k)
{
  return (op->m + op->mark + k)%op->m;
}

static void
finalize_lbfgs_operator(opk_operator_t* self)
{
  opk_lbfgs_operator_t* op = (opk_lbfgs_operator_t*)self;
  opk_index_t k, m = op->m;

  for (k = 0; k < m; ++k) {
    OPK_DROP(op->s[k]);
    OPK_DROP(op->y[k]);
  }
  OPK_DROP(op->tmp);
  OPK_DROP(op->H0);
}

static int
apply_lbfgs_operator(opk_operator_t* self, opk_vector_t* dst,
                     const opk_vector_t* src)
{
  opk_lbfgs_operator_t* op = (opk_lbfgs_operator_t*)self;
  opk_vector_t* tmp;
  const opk_index_t newest = 0;
  opk_index_t oldest = 1 - op->mp;
  opk_index_t j, k;

  /* Initialize work vectors. */
  if (op->H0 == NULL) {
    /* No preconditioner, no needs for a scratch vector. */
    tmp = dst;
  } else {
    /* With a preconditioner, a scratch vector is needed. */
    if (op->tmp == NULL) {
      op->tmp = opk_vcreate(op->base.outspace);
      if (op->tmp == NULL) {
        return OPK_FAILURE;
      }
    }
    tmp = op->tmp;
  }

  /* First loop of the recursion (from the newest saved pair to the oldest
     one). */
  opk_vcopy(tmp, src);
  for (k = newest; k >= oldest; --k) {
    j = slot(op, k);
    if (op->rho[j] > 0.0) {
      op->beta[j] = op->rho[j]*opk_vdot(tmp, op->s[j]);
      opk_vaxpby(tmp, 1.0, tmp, -op->beta[j], op->y[j]);
    } else {
      op->beta[j] = 0.0;
    }
  }

  /* Apply approximation of inverse Hessian. */
  if (op->H0 != NULL) {
    if (opk_apply_direct(op->H0, dst, tmp) != OPK_SUCCESS) {
      return OPK_FAILURE;
    }
  }
  if (op->gamma != 1.0) {
    opk_vscale(dst, op->gamma, dst);
  }

  /* Second loop of the recursion (from the oldest saved pair to the
   * newest one). */
  for (k = oldest; k <= newest; ++k) {
    j = slot(op, k);
    if (op->rho[j] > 0.0) {
      double phi = op->rho[j]*opk_vdot(dst, op->y[j]);
      opk_vaxpby(dst, 1.0, dst, op->beta[j] - phi, op->s[j]);
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
  size_t s_offset, y_offset, beta_offset, rho_offset, size;
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
  s_offset = ROUND_UP(sizeof(opk_lbfgs_operator_t),
                      sizeof(opk_vector_t*));
  y_offset = ROUND_UP(s_offset + m*sizeof(opk_vector_t*),
                      sizeof(opk_vector_t*));
  beta_offset = ROUND_UP(y_offset + m*sizeof(opk_vector_t*),
                         sizeof(double));
  rho_offset = ROUND_UP(beta_offset + m*sizeof(double),
                        sizeof(double));
  size = rho_offset + m*sizeof(double);
  op = (opk_lbfgs_operator_t*)opk_allocate_operator(&lbfgs_operations,
                                                    vspace, vspace, size);
  if (op == NULL) {
    return NULL;
  }
  op->s = (opk_vector_t**)(((char*)op) + s_offset);
  op->y = (opk_vector_t**)(((char*)op) + y_offset);
  op->beta = (double*)(((char*)op) + beta_offset);
  op->rho = (double*)(((char*)op) + rho_offset);
  op->epsilon = 1e-6;
  op->gamma = 1.0;
  op->m = m;
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
  opk_reset_lbfgs_operator(op);
  return op;
 error:
  OPK_DROP(op);
  return NULL;
}

void
opk_reset_lbfgs_operator(opk_lbfgs_operator_t* op)
{
  op->mp = 0;
  op->mark = -1;
  op->gamma = 1.0;
}

opk_vector_t*
opk_get_lbfgs_s(opk_lbfgs_operator_t* op, opk_index_t k)
{
  return op->s[slot(op, k)];
}

opk_vector_t*
opk_get_lbfgs_y(opk_lbfgs_operator_t* op, opk_index_t k)
{
  return op->y[slot(op, k)];
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
                          const opk_vector_t* x1,
                          const opk_vector_t* x0,
                          const opk_vector_t* g1,
                          const opk_vector_t* g0)
{
  double sty, ynorm, snorm;

  /* Store the variables and gradient differences in the slot just after the
     mark (which is the index of the last saved pair or -1 if none). */
  opk_index_t j = slot(op, 1);
  opk_vaxpby(op->s[j], 1.0, x1, -1.0, x0);
  snorm = opk_vnorm2(op->s[j]);
  opk_vaxpby(op->y[j], 1.0, g1, -1.0, g0);
  ynorm = opk_vnorm2(op->y[j]);

  /* Compute RHO[j] and GAMMA.  If the update formula for GAMMA does not yield
     a strictly positive value, the strategy is to keep the previous value. */
  sty = opk_vdot(op->s[j], op->y[j]);
  if (sty <= op->epsilon*snorm*ynorm) {
    /* This pair will be skipped. */
    op->rho[j] = 0.0;
  } else {
    /* Compute RHO[j] and GAMMA. */
    op->rho[j] = 1.0/sty;
    if (op->scaling == OPK_SCALING_OREN_SPEDICATO) {
      op->gamma = (sty/ynorm)/ynorm;
    } else if (op->scaling == OPK_SCALING_BARZILAI_BORWEIN) {
      op->gamma = (snorm/sty)*snorm;
    }
    /* Update the mark and the number of saved pairs. */
    op->mark = j;
    if (op->mp < op->m) {
      ++op->mp;
    }
  }
}

double
opk_get_lbfgs_operator_update_threshold(opk_lbfgs_operator_t* op)
{
  return op->epsilon;
}

opk_status_t
opk_set_lbfgs_operator_update_threshold(opk_lbfgs_operator_t* op,
                                        double epsilon)
{
  if (epsilon < 0.0 || epsilon >= 1.0) {
    errno = EINVAL;
    return OPK_FAILURE;
  }
  op->epsilon = epsilon;
  return OPK_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/* VARIABLE METRIC OPTIMIZER */

/* FIXME: make a common super-type for NLCG and VMLM to share API. */

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
static const double DELTA = 0.01;
static const double EPSILON = 1e-3;
static const int    SCALING = OPK_SCALING_OREN_SPEDICATO;

struct _opk_vmlm {
  opk_object_t base; /**< Base type (must be the first member). */
  double alpha;      /**< Current step size. */
  double delta;      /**< Threshold to accept descent direction. */
  double epsilon;    /**< Relative size for a small step. */
  double grtol;      /**< Relative threshold for the norm or the gradient
                          (relative to GINIT the norm of the initial gradient)
                          for convergence. */
  double gatol;      /**< Absolute threshold for the norm or the gradient for
                          convergence. */
  double ginit;      /**< Euclidean norm or the initial gradient. */
  double f0;         /**< Function value at x0. */
  double dg0;        /**< Directional derivative at x0. */
  double g0norm;     /**< Euclidean norm of g0 the gradient at x0. */
  double g1norm;     /**< Euclidean norm of the gradient at the last accepted
                          step. */
  double pnorm;      /**< Euclidean norm of p the anti-search direction */
  double stpmin;
  double stpmax;

  opk_vspace_t* vspace;
  opk_lnsrch_t* lnsrch;
  opk_lbfgs_operator_t* H;
  opk_vector_t*  x0; /**< Variables at the start of the line search. */
  opk_vector_t*  g0; /**< Gradient at x0. */
  opk_vector_t*  p;  /**< Anti-search direction; a iterate is computed
                          as: x1 = x0 - alpha*p with alpha > 0. */

  opk_index_t evaluations; /**< Number of functions (and gradients)
                                evaluations. */
  opk_index_t iterations;  /**< Number of iterations (successful steps
                                taken). */
  opk_index_t restarts;
  opk_task_t task;
  int reason; /* some details about the error */
  opk_bool_t starting;
  opk_bool_t save_memory; /**< To save space, the variable and gradient at the
                               start of a line search are weak references to
                               the (s,y) pair of vectors of the LBFGS operator
                               just after the mark. */
};

typedef enum {
  OPK_NO_PROBLEMS = 0,
  OPK_BAD_PRECONDITIONER,  /**< Preconditioner is not positive definite. */
  OPK_LINE_SEARCH_WARNING, /**< Warning in line search. */
  OPK_LINE_SEARCH_ERROR,   /**< Error in line search. */
} opk_reason_t;

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
  OPK_DROP(opt->p);
}

static opk_task_t
optimizer_success(opk_vmlm_t* opt, opk_task_t task)
{
  opt->reason = OPK_NO_PROBLEMS;
  opt->task = task;
  return task;
}

static opk_task_t
optimizer_failure(opk_vmlm_t* opt, opk_reason_t reason)
{
  opt->reason = reason;
  opt->task = OPK_TASK_ERROR;
  return opt->task;
}

static opk_task_t
line_search_failure(opk_vmlm_t* opt)
{
  opk_reason_t reason;
  if (opk_lnsrch_has_errors(opt->lnsrch)) {
    reason = OPK_LINE_SEARCH_ERROR;
  } else {
    reason = OPK_LINE_SEARCH_WARNING;
  }
  return optimizer_failure(opt, reason);
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
  opt->p = opk_vcreate(vspace);
  if (opt->p == NULL) {
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
  opt->starting = TRUE;
  opk_reset_lbfgs_operator(opt->H);
  return optimizer_success(opt, OPK_TASK_COMPUTE_FG);
}

/* FIXME: There should be some thread safe global data used to store last error
   information.   For now, I use the thread safe global variable errno. */

/* FIXME: check order of arguments. */


/* Notes:
 *
 *  - during line search, ws->s[ws->mark] is the (anti-)search direction p, and
 *    ws->y[ws->mark] is g0, the gradient at the start of the line search;
 *
 *  - gd and gd0 are the dot product of the gradient and search direction;
 *
 * To do:
 *
 *  - To save a copy operation it may be better to store the search direction p
 *    in a separate slot and to store x0 in ws->s[ws->mark] and g0 in
 *    ws->y[ws->mark].  This is also advantageous if the effective step (rather
 *    than the scaled search direction) is used to compute the variables change
 *    stored in ws->s.  Finally this is more logical for the organization of the
 *    code (e.g. the LBFGS operator is applied to p not to x0).
 */

static opk_task_t
next_step(opk_vmlm_t* opt, opk_vector_t* x1)
{
  opk_vaxpby(x1, 1.0, opt->x0, -opt->alpha, opt->p);
  return optimizer_success(opt, OPK_TASK_COMPUTE_FG);
}

opk_task_t
opk_iterate_vmlm(opk_vmlm_t* opt, opk_vector_t* x1,
                 double f1, opk_vector_t* g1)
{
  double gtest, pg1;
  int status;

  switch (opt->task) {

  case OPK_TASK_COMPUTE_FG:

    /* Caller has computed the function value and the gradient at the current
       point. */
    ++opt->evaluations;
    if (opt->evaluations > 1) {
      /* A line search is in progress.  Compute directional derivative and
         check whether line search has converged. */
      pg1 = opk_vdot(opt->p, g1);
      status = opk_lnsrch_iterate(opt->lnsrch, &opt->alpha, f1, -pg1);
      if (status == OPK_LNSRCH_SEARCH) {
        /* Line search has not yet converged. */
        return next_step(opt, x1);
      }
      if (status != OPK_LNSRCH_CONVERGENCE &&
          status != OPK_LNSRCH_WARNING_ROUNDING_ERRORS_PREVENT_PROGRESS) {
        return line_search_failure(opt);
      }
      ++opt->iterations;
    }

    /* The current step is acceptable. Check for global convergence. */
    opt->g1norm = opk_vnorm2(g1);
    if (opt->evaluations == 1) {
      opt->ginit = opt->g1norm;
    }
    gtest = max3(0.0, opt->gatol, opt->grtol*opt->ginit);
    return optimizer_success(opt, ((opt->g1norm <= gtest)
                                   ? OPK_TASK_FINAL_X
                                   : OPK_TASK_NEW_X));

  case OPK_TASK_NEW_X:

    if (opt->evaluations > 1) {
      /* Update the LBFGS matrix. */
      opk_update_lbfgs_operator(opt->H,
                                x1, opt->x0,
                                g1, opt->g0);
    }

  case OPK_TASK_FINAL_X:

    /* Compute a search direction.  We take care of checking whether D = -P is
       a sufficient descent direction.  As shown by Zoutendijk, this is true
       if: cos(theta) = -(D/|D|)'.(G/|G|) >= EPSILON > 0 where G is the
       gradient and D the descent direction. */
    while (TRUE) {
      opk_apply_direct((opk_operator_t*)opt->H, opt->p, g1);
      opt->pnorm = opk_vnorm2(opt->p);
      pg1 = opk_vdot(opt->p, g1);
      if (pg1 >= opt->delta*opt->pnorm*opt->g1norm) {
        /* Accept P (respectively D = -P) as a sufficient ascent (respectively
           descent) direction and set the directional derivative. */
        opt->dg0 = -pg1;
        break;
      }
      if (opt->H->mp < 1) {
        /* Initial iteration or recursion has just been restarted.  This means
           that the initial inverse Hessian approximation is not positive
           definite. */
        return optimizer_failure(opt, OPK_BAD_PRECONDITIONER);
      }
      /* Reset the LBFGS recursion and loop to use H0 to compute an initial
         search direction. */
      opk_reset_lbfgs_operator(opt->H);
      ++opt->restarts;
    }

    /* Save current variables X0, gradient G0 and function value F0.  The
       directional derivative DG0 = <D,G0> has already been set in the above
       loop.  */
    if (opt->save_memory) {
      /* Use the slot just after the mark to store X0 and G0.  Note that this
         is a weak reference: we do not "hold" the vectors. */
      opt->x0 = opk_get_lbfgs_s(opt->H, 1);
      opt->g0 = opk_get_lbfgs_y(opt->H, 1);
      if (opt->H->mp > opt->H->m - 1) {
        opt->H->mp = opt->H->m - 1;
      }
    }
    opk_vcopy(opt->x0, x1);
    opk_vcopy(opt->g0, g1);
    opt->g0norm = opt->g1norm;
    opt->f0 = f1;

    /* Estimate the length of the first step, start the line search and take
       the first step along the search direction. */
    if (opt->H->mp >= 1 || opt->H->scaling == OPK_SCALING_NONE) {
      opt->alpha = 1.0;
    } else if (0.0 < opt->epsilon && opt->epsilon < 1.0) {
      double x1norm = opk_vnorm2(x1);
      if (x1norm > 0.0) {
        opt->alpha = (x1norm/opt->g1norm)*opt->epsilon;
      } else {
        opt->alpha = 1.0/opt->g1norm;
      }
    } else {
      opt->alpha = 1.0/opt->g1norm;
    }
    status = opk_lnsrch_start(opt->lnsrch, opt->f0, opt->dg0, opt->alpha,
                              opt->stpmin*opt->alpha,
                              opt->stpmax*opt->alpha);
    if (status != OPK_LNSRCH_SEARCH) {
      return line_search_failure(opt);
    }
    return next_step(opt, x1);

  default:

    /* There must be something wrong. */
    return opt->task;

  }

}

opk_task_t
opk_get_vmlm_task(opk_vmlm_t* opt)
{
  return opt->task;
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

int
opk_set_vmlm_scaling(opk_vmlm_t* opt, int scaling)
{
  if (opt == NULL) {
    errno = EFAULT;
    return OPK_FAILURE;
  }
  switch (scaling) {
  case OPK_SCALING_NONE:
  case OPK_SCALING_OREN_SPEDICATO:
  case OPK_SCALING_BARZILAI_BORWEIN:
    opt->H->scaling = scaling;
    return OPK_SUCCESS;
  default:
    errno = EINVAL;
    return OPK_FAILURE;
  }
}

double
opk_get_vmlm_gatol(opk_vmlm_t* opt)
{
  return (opt == NULL ? GATOL : opt->gatol);
}

int
opk_set_vmlm_gatol(opk_vmlm_t* opt, double gatol)
{
  if (opt == NULL) {
    errno = EFAULT;
    return OPK_FAILURE;
  }
  if (non_finite(gatol) || gatol < 0.0) {
    errno = EINVAL;
    return OPK_FAILURE;
  }
  opt->gatol = gatol;
  return OPK_SUCCESS;
}

double
opk_get_vmlm_grtol(opk_vmlm_t* opt)
{
  return (opt == NULL ? GRTOL : opt->grtol);
}

int
opk_set_vmlm_grtol(opk_vmlm_t* opt, double grtol)
{
  if (opt == NULL) {
    errno = EFAULT;
    return OPK_FAILURE;
  }
  if (non_finite(grtol) || grtol < 0.0) {
    errno = EINVAL;
    return OPK_FAILURE;
  }
  opt->grtol = grtol;
  return OPK_SUCCESS;
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

int
opk_set_vmlm_stpmin_and_stpmax(opk_vmlm_t* opt, double stpmin, double stpmax)
{
  if (opt == NULL) {
    errno = EFAULT;
    return OPK_FAILURE;
  }
  if (non_finite(stpmin) || non_finite(stpmax) ||
      stpmin < 0.0 || stpmax <= stpmin) {
    errno = EINVAL;
    return OPK_FAILURE;
  }
  opt->stpmin = stpmin;
  opt->stpmax = stpmax;
  return OPK_SUCCESS;
}

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
