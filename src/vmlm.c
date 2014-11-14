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

/**
 * @defgroup VariableMetric     Variable metric methods.
 * @addtogroup VariableMetric
 * @{
 */

typedef struct _opk_vmlm opk_vmlm_t;
typedef struct _opk_lbfgs_operator opk_lbfgs_operator_t;

typedef enum {
  OPK_CUSTOM_APPROX = 0,
  OPK_BARZILAI_BORWEIN_1, /**< gamma = <s,y>/<y,y> */
  OPK_BARZILAI_BORWEIN_2, /**< gamma = <s,s>/<s,y> */
} opk_inverse_hessian_rule_t;

/**
 * Create a new limited memory BFGS operator.
 *
 * @param vspace - The input and output vector space of the operator.
 * @param m - The maximum number of previous steps to memorize.
 * @param rule - The rule for updating the scale of the approximation of the
 *               inverse Hessian.
 */
extern opk_lbfgs_operator_t*
opk_new_lbfgs_operator(opk_vspace_t* vspace, opk_index_t m,
                       opk_inverse_hessian_rule_t rule);

/** Forget all memorized steps in limited memory BFGS operator. */
extern void
opk_reset_lbfgs_operator(opk_lbfgs_operator_t* op);

/**
 * Query a memorized variable difference from a limited memory BFGS operator.
 *
 * It is the caller responsibility to use proper arguments.
 *
 * @param op - A limited memory BFGS operator.
 * @param k - The index of the memorized step to consider.
 *
 * @return s_k
 */
extern opk_vector_t*
opk_get_lbfgs_s(opk_lbfgs_operator_t* op, opk_index_t k);

/**
 * Query a memorized gradient difference from a limited memory BFGS operator.
 *
 * It is the caller responsibility to use proper arguments.
 *
 * @param op - A limited memory BFGS operator.
 * @param k - The index of the memorized step to consider.
 *
 * @return y_k
 */
extern opk_vector_t*
opk_get_lbfgs_y(opk_lbfgs_operator_t* op, opk_index_t k);

extern void
opk_set_lbfgs_operator_preconditioner(opk_lbfgs_operator_t* op,
                                      opk_operator_t* B0);

/**
 * Update LBFGS operator with a new pair of variables and gradient
 * differences.
 * @param x1 - The new variables.
 * @param x0 - The previous variables.
 * @param g1 - The gradient at {@code x1}.
 * @param g0 - The gradient at {@code x0}.
 * @throws IncorrectSpaceException
 */
extern void
opk_lbfgs_operator_update(opk_lbfgs_operator_t* op,
                          const opk_vector_t* x1,
                          const opk_vector_t* x0,
                          const opk_vector_t* g1,
                          const opk_vector_t* g0);

extern opk_vmlm_t*
opk_new_vmlm_optimizer(opk_vspace_t* vspace, opk_index_t m,
                       double frtol, double fatol, double fmin);

extern opk_task_t
opk_vmlm_start(opk_vmlm_t* opt);

/** @} */

/*---------------------------------------------------------------------------*/
/* LIMITED MEMORY BFGS OPERATOR */

struct _opk_lbfgs_operator {
  opk_operator_t base; /**< Base type (must be the first member). */
  double gamma;        /**< Scale factor to approximate inverse Hessian. */
  opk_vector_t** s;
  opk_vector_t** y;
  opk_vector_t* tmp;   /**< Scratch vector needed if there is a
                          preconditioner. */
  opk_operator_t* B0;  /**< Preconditioner or NULL. */
  double* alpha;
  double* rho;
  opk_index_t m;      /* maximum number of memorized (s,y) pairs */
  opk_index_t mp;     /* actual number of saved (s,y) pairs (0 <= mp <= m) */
  opk_index_t mark;   /* index of oldest saved pair */
  opk_inverse_hessian_rule_t rule;
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
  OPK_DROP(op->B0);
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
  if (op->B0 == NULL) {
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
      op->alpha[j] = op->rho[j]*opk_vdot(tmp, op->s[j]);
      opk_vaxpby(tmp, 1.0, tmp, -op->alpha[j], op->y[j]);
    } else {
      op->alpha[j] = 0.0;
    }
  }

  /* Apply approximation of inverse Hessian. */
  if (op->B0 != NULL) {
    if (opk_apply_direct(op->B0, tmp, dst) != OPK_SUCCESS) {
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
      double beta = op->rho[j]*opk_vdot(dst, op->y[j]);
      opk_vaxpby(dst, 1.0, dst, op->alpha[j] - beta, op->s[j]);
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
                       opk_inverse_hessian_rule_t rule)
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
  s_offset = ROUND_UP(sizeof(opk_lbfgs_operator_t),
                      sizeof(opk_vector_t*));
  y_offset = ROUND_UP(s_offset + m*sizeof(opk_vector_t*),
                      sizeof(opk_vector_t*));
  alpha_offset = ROUND_UP(y_offset + m*sizeof(opk_vector_t*),
                          sizeof(double));
  rho_offset = ROUND_UP(alpha_offset + m*sizeof(double),
                        sizeof(double));
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
  op->m = m;
  opk_reset_lbfgs_operator(op);
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
  op->rule = rule;
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
  return op->s[k];
}

opk_vector_t*
opk_get_lbfgs_y(opk_lbfgs_operator_t* op, opk_index_t k)
{
  return op->y[k];
}

void
opk_set_lbfgs_operator_preconditioner(opk_lbfgs_operator_t* op,
                                      opk_operator_t* B0)
{
  if (op->B0 != B0) {
    opk_operator_t* old = op->B0;
    op->B0 = OPK_HOLD_OPERATOR(B0);
    OPK_DROP(old);
  }
}

void
opk_lbfgs_operator_update(opk_lbfgs_operator_t* op,
                          const opk_vector_t* x1,
                          const opk_vector_t* x0,
                          const opk_vector_t* g1,
                          const opk_vector_t* g0)
{
  double sty;

  /* Store the variables and gradient differences in the slot just after the
     mark (which is the index of the last saved pair or -1 if none). */
  opk_index_t j = slot(op, 1);
  opk_vaxpby(op->s[j], 1.0, x1, -1.0, x0);
  opk_vaxpby(op->y[j], 1.0, g1, -1.0, g0);

  /* Compute RHO[j] and GAMMA.  If the update formula for GAMMA does not yield
     a strictly positive value, the strategy is to keep the previous value. */
  sty = opk_vdot(op->s[j], op->y[j]);
  if (sty > 0.0) {
    op->rho[j] = 1.0/sty;
    if (op->rule == OPK_BARZILAI_BORWEIN_1) {
      double yty = opk_vdot(op->y[j], op->y[j]);
      if (yty > 0.0) {
        op->gamma = sty/yty;
      }
    } else if (op->rule == OPK_BARZILAI_BORWEIN_2) {
      double sts = opk_vdot(op->s[j], op->s[j]);
      if (sts > 0.0) {
        op->gamma = sts/sty;
      }
    }
  } else {
    /* This pair will be skipped. */
    op->rho[j] = 0.0;
  }

  /* Update the mark and the number of saved pairs. */
  op->mark = j;
  if (op->mp < op->m) {
    ++op->mp;
  }
}

/*---------------------------------------------------------------------------*/
/* VARIABLE METRIC OPTIMIZER */

/* FIXME: make a common super-type for NLCG and VMLM to share API. */

struct _opk_vmlm {
  opk_object_t base; /**< Base type (must be the first member). */
  double alpha;      /**< Current step size. */
  double epsilon;    /**< Threshold to accept descent direction. */
  double grtol;      /**< Relative threshold for the norm or the gradient
                          (relative to GTEST the norm of the initial gradient)
                          for convergence. */
  double gatol;      /**< Absolute threshold for the norm or the gradient for
                          convergence. */
  double gtest;      /**< Euclidean norm or the initial gradient FIXME:
                          GINIT */
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
  opk_lbfgs_operator_t* B;
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
  opk_bool_t starting;
  opk_bool_t save_memory; /**< To save space, the variable and gradient at the
                               start of a line search are references to the
                               (s,y) pair of vectors of the LBFGS operator just
                               after the mark. */
};

static void
finalize_vmlm(opk_object_t* obj)
{
  opk_vmlm_t* opt = (opk_vmlm_t*)obj;
  OPK_DROP(opt->vspace);
  OPK_DROP(opt->lnsrch);
  OPK_DROP(opt->B);
  OPK_DROP(opt->x0);
  OPK_DROP(opt->g0);
  OPK_DROP(opt->x0);
  OPK_DROP(opt->p);
}

opk_vmlm_t*
opk_new_vmlm_optimizer_with_line_search(opk_vspace_t* vspace,
                                        opk_index_t m,
                                        opk_lnsrch_t* lnsrch,
                                        double frtol,
                                        double fatol,
                                        double fmin)
{
  size_t s_offset, y_offset, alpha_offset, size;
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
  if (frtol <= 0.0) {
    errno = EINVAL;
    return NULL;
  }
  if (fatol <= 0.0) {
    errno = EINVAL;
    return NULL;
  }

  /* Allocate enough memory for the workspave and its arrays. */
  s_offset = ROUND_UP(sizeof(opk_vmlm_t),
                      sizeof(opk_vector_t*));
  y_offset = ROUND_UP(s_offset + m*sizeof(opk_vector_t*),
                      sizeof(opk_vector_t*));
  alpha_offset = ROUND_UP(y_offset + m*sizeof(opk_vector_t*),
                          sizeof(double));
  size = alpha_offset + m*sizeof(double);
  opt = (opk_vmlm_t*)opk_allocate_object(finalize_vmlm, size);
  if (opt == NULL) {
    return NULL;
  }

  /* Instanciate the workspace (not the part which is done by
     opk_vmlm_start). */
  opt->vspace = OPK_HOLD_VSPACE(vspace);
  opt->lnsrch = OPK_HOLD_LNSRCH(lnsrch);
  opt->B = opk_new_lbfgs_operator(vspace, m, OPK_BARZILAI_BORWEIN_2);
  if (opt->B == NULL) {
    goto error;
  }
#if 0
  opt->frtol = frtol;
  opt->fatol = fatol;
  opt->fmin = fmin;
#endif
  opt->x0 = opk_vcreate(vspace);
  if (opt->x0 == NULL) {
    goto error;
  }
  opt->g0 = opk_vcreate(vspace);
  if (opt->g0 == NULL) {
    goto error;
  }
  opt->p = opk_vcreate(vspace);
  if (opt->p == NULL) {
    goto error;
  }

  opk_vmlm_start(opt);
  return opt;

error:
  OPK_DROP(opt);
  return NULL;
}

opk_vmlm_t*
opk_new_vmlm_optimizer(opk_vspace_t* vspace,
                       opk_index_t m,
                       double frtol,
                       double fatol,
                       double fmin)
{
  opk_lnsrch_t* lnsrch;
  opk_vmlm_t* opt;
  lnsrch = opk_lnsrch_new_csrch(/* sftol  */  1e-4,
                                /* sgtol  */  1e-1,
                                /* sxtol  */  1E-17);
  if (lnsrch == NULL) {
    return NULL;
  }
  opt = opk_new_vmlm_optimizer_with_line_search(vspace, m, lnsrch,
                                                frtol, fatol, fmin);
  OPK_DROP(lnsrch); /* the line search is now owned by the optimizer */
  return opt;
}

opk_task_t
opk_vmlm_start(opk_vmlm_t* opt)
{
  opt->iterations = 0;
  opt->starting = TRUE;
  opk_reset_lbfgs_operator(opt->B);
  opt->task = OPK_TASK_COMPUTE_FG;
  return opt->task;
}

#if 0

/* The three possible stages. */
#define START   0 /* just (re)started */
#define SEARCH  1 /* search in progress */
#define FINISH  2 /* convergence, errors or warnings */


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
next_step(opk_vector_t* x1)
{
  double alpha = opk_lnsrch_get_step(ws->lnsrch);
  opk_vaxpby(x1, 1.0, ws->x0, -alpha, ws->p);
  ws->task = COMPUTE_FG;
  return ws->task;
}

opk_task_t
opk_vmlm_iterate(opk_vmlm_t* ws, opk_vector_t* x,
                 double f1, opk_vector_t* g)
{
  double df, df1, df2;
  opk_vspace_t* vspace = ws->vspace;
  int corrupted;

  /* New iterate is computed as:
   *
   *   x(alpha) = x0 + alpha*sigma*p
   *
   * where x0 is the point at the start of the line search, p is the search
   * direction and alpha > 0 the step length.  Hence the derivative along
   * the search direction is:
   *
   *   d f(x0 + alpha*sigma*p) / d alpha = sigma*dot(p,g(x(alpha)))
   */
  const double SIGMA = -1.0; /* +/- 1 */

  /* Check consistency of internal state and return if not searching.
     Making this early simplifies the various alternatives in the
     remaining part of the code. */
  if (ws->task == COMPUTE_FG || ws->task == NEW_X) {
    corrupted = (ws->stage != START && ws->stage != SEARCH);
  } else if (ws->stage == FINISH) {
    return SUCCESS;
  } else {
    corrupted = TRUE;
  }
  if (corrupted) {
    ws->reason = CORRUPTED_WORKSPACE;
    return (ws->task = OPK_TASK_ERROR);
  }


  if (ws->task == OPK_TASK_COMPUTE_FG) {
    /* Caller has computed the function value and the gradient at the current
       point. */
    ++ws->evaluations;
    if (evaluations > 1) {
      /* A line search is in progress.  Compute directional derivative and
       * check whether line search has converged. */
      int status = lnsrch.iterate(alpha, f1, -p.dot(g1));
      if (status == LineSearch.SEARCH) {
        return nextStep(x1);
      }
      if (status != LineSearch.CONVERGENCE &&
          status != LineSearch.WARNING_ROUNDING_ERRORS_PREVENT_PROGRESS) {
        return lineSearchFailure(status);
      }
      ++ws->iterations;
    }

    /* The current step is acceptable. Check for global convergence. */
    g1norm = opk_vnorm2(g1);
    if (ws->evaluations == 1) {
      ws->gtest = g1norm;
    }
    reason = NO_PROBLEM;
    task = (g1norm <= getGradientThreshold() ? OptimTask.FINAL_X : OptimTask.NEW_X);

  } else if (task == OptimTask.NEW_X || task == OptimTask.FINAL_X) {



    opk_bool_t accept;

    if (ws->starting) {
      /* Make the initial variables and function value available for
         inspection. */
      ws->f0 = f1;
      ws->stp = 1.0;
      ws->gd0 = opk_dot(g1, g1); /* FIXME: sign? */
      accept = TRUE;
    } else {
      /* Check whether line search has converged. */
      double df1 = SIGMA*opk_dot(ws->p, g);
      status = opk_lnsrch_iterate(ws->lnsrch, &ws->stp, f1, df1);
      if (status == OPK_LNSRCH_CONVERGENCE ||
          status == OPK_LNSRCH_WARNING_ROUNDING_ERRORS_PREVENT_PROGRESS) {
        /* Line search has converged, a new success step is available
           for inspection. */
        ++ws->iter;
        accept = TRUE;
      } else if (status == OPK_LNSRCH_SEARCH) {
        /* Line search has not converged, proceed with next step. */
        accept = FALSE;
      } else {
        /* Some error has occured. */
        if (status < 0) {
          ws->reason = OPK_VMLM_ERROR_LINE_SEARCH;
          ws->task = OPK_TASK_ERROR;
        } else {
          ws->reason = OPK_VMLM_WARNING_LINE_SEARCH;
          ws->task = OPK_TASK_WARNING;
          /* FIXME: shall we restore x0? */
        }
        return ws->task;
      }
    }
    if (accept) {
      /* Check for global convergence.  Otherwise a new iterate is
         available for inspection. */
      double t1 = fabs(f1 - ws->f0);
      double t2 = ws->stp*fabs(ws->gd0);
      double df = MAX(t1, t2);
      if (df <= ws->frtol*fabs(ws->f0)) {
        ws->reason = OPK_VMLM_CONVERGENCE_FRTOL_TEST_SATISFIED;
        ws->task = OPK_TASK_FINAL_X;
      } else if (df <= ws->fatol) {
        ws->reason = OPK_VMLM_CONVERGENCE_FATOL_TEST_SATISFIED;
        ws->task = OPK_TASK_FINAL_X;
      } else {
        /* A new iterate is available for inspection. */
        ws->reason = OPK_VMLM_NONE;
        ws->task = OPK_TASK_NEW_X;
      }
      return ws->task;
    }

  } else if (ws->task == OPK_TASK_NEW_X) {

    /* Caller returned after an initial or a successful step.  We
       first check for convergence, then compute a new search
       direction, start the line search and proceed with the next
       step. */
    double stp, stpmin, stpmax;
    if (ws->starting) {
      /* Initialize for very first search direction. */
      if (f1 <= ws->fmin) {
        ws->task = OPK_VMLM_ERROR_INITIAL_F_LE_FMIN;
        ws->stage = FINISH;
        return OPK_FAILURE;
      }
      ws->mp = 0; /* FIXME: not needed? */
      ws->scale = V_NORM_2(g);
      V_SCALE(1.0/ws->scale, g); /* FIXME: why? */
      stp = ...;
      stpmin = ...;
      stmpax = ..;
    } else {
      /* Store new (s,y) pair and apply recurence to compute a new
         search direction. */
    }

    opk_lnsrch_start(ws->lnsrch, f1, df, stp, stpmin, stpmax);
  }

  /* Compute new variables to try. */
  V_AXPBY(1.0, ws->x0, SIGMA*(*stp_ptr), ws->p, x);
  ws->task = COMPUTE_FG;
  ws->stage = SEARCH;
  return OPL_SUCCESS;


  if (ws->stage == SEARCH) {

    /* Search is in progress, possible values for task are NEW_X or
       COMPUTE_FG. */
    int status;

    if (ws->task == NEW_X) {

      ws->on_new_x(ws, x, fx_ptr, gx);
      if (ws->stage == FINISH) {
        return (ws->task < 0 ? FAILURE : SUCCESS);
      }

      status = ws->lnsrch->start(ws->lnsrch, &ws->stp, ws->f0, ws->d0);

    } else {

      if (ws->task != COMPUTE_FG) {
        return corrupted(ws);
      }

      /* Compute derivative along search direction:
         d f(x0 + alpha*sigma*p) / d alpha. */
      double dx = SIGMA*opk_dot(ws->p, gx);
      status = ws->lnsrch->iterate(ws->lnsrch, &ws->stp, fx, dx);

    }

    /* Check line search state. */
    if (status < 0) {
      /* An error occured. */
      if (ws->task != NEW_X /* FIXME: *fx_ptr > ws->f0 */) {
        /* Restore best solution found so far. */
        V_COPY(ws->x0, x);
        *fx_ptr = ws->f0;
        V_COPY(ws->g0, gx);
      }
      ws->stage = FINISH;
      ws->task = ERROR_LINE_SEARCH; /* FIXME: status */
      return FAILURE;
    } else if (status > 0) {
      /* Line search converged (FIXME: check warnings). */
      ++iter;
      ws->task = NEW_X;
    } else {
      /* A new iterate is required by the line search procedure. */
      V_AXPBY(1.0, ws->x0, SIGMA*ws->stp, p, x);
    }
    ws->stage = SEARCH; /* FIXME: not needed */
    return SUCCESS;

  } /* ws->task == SEARCH */






  if (ws->stage == START) {

   if (f <= ws->fmin) {
     /* ERROR: INITIAL F .LE. FMIN */
     errno = EINVAL;
     return OPK_FAILURE;
   }
   ws->iter = 1; /* FIXME: not needed? */

    /* Initialize step information. */
    ws->scale = V_NORM_2(g);
    V_SCALE(1.0/ws->scale, g);

    /* Set work to start the search. */
    opk_lnsrch_start(ws->lnsrch, f1, df, stp, stpmin, stpmax);
line_search_task = LINE_SEARCH_TASK_START_SEARCH; /* "START SEARCH" */

  } else {

    /* Restore local variables. */
    if (ws->state == 1) {
      ws->line_search_task = LINE_SEARCH_TASK_SEARCH; /* "SEARCH" */
    }
    if (ws->state == 2) {
      ws->line_search_task = LINE_SEARCH_TASK_SEARCH_DIRECTION; /* "SEARCH DIRECTION" */
    }
  }
L20:

  if (ws->line_search_task = LINE_SEARCH_TASK_START_SEARCH) { /* "START SEARCH" */
    /* Initialize the line search subroutine. */
    double stpmin, stpmax;
    ws->f0 = f;
    ws->gd0 = -opk_dot(g, ws->s[ws->mark]); /* FIXME: use ws->p (see my notes) */
    stpmin = 0.0;
    stpmax = (ws->fmin - ws->f0)/(ws->sgtol*ws->gd0);
    ws->stp = MIN(1.0, stpmax);
    V_COPY(x, ws->x0); /* FIXME: use ws->s[ws->mark] (see my notes) */
    V_COPY(g, ws->y[ws->mark]);
    opk_set_csrch_bounds(ws->lnsrch, stpmin, stpmax);
    ws->task = START_SEARCH;
    lnsrch_task = OPK_LNSRCH_TASK_SEARCH;
    //opk_set_csrch_task(ws->lnsrch, OPK_LNSRCH_TASK_SEARCH);
  }
  if (lnsrch_task == OPK_LNSRCH_TASK_SEARCH) {
    /* Determine the line search parameter. */
    if (f < ws->fmin) {
      ws->task = OPK_VMLM_WARNING_F_LT_FMIN;
      goto done;
    }
    ws->gd = -opk_dot(g, ws->s[ws->mark]);
    res = opk_csrch(ws->lnsrch, &ws->stp, f, &ws->gd);

    /* Compute the new iterate. */
    V_AXPBY(1.0, ws->x0, -ws->stp, ws->s[ws->mark], x);

    /* Continue if the line search has converged. */
    if (s_cmp(task, "CONV") != 0 && s_cmp(task, "WARNING: XTOL TEST SATISFIED") != 0) {
      goto done;
    }

    /* Compute the step and gradient change. */
    ++ws->iter;
    V_AXPY(-1.0, g, ws->y[ws->mark]);
    V_SCALE(ws->stp, ws->s[ws->mark]);
    rho[ws->mark] = opk_dot(ws->y[ws->mark], ws->s[ws->mark]);

    /* Compute the scale. */
    if (rho[ws->mark] > 0.0) {
      ws->scale = rho[ws->mark]/opk_dot(ws->y[ws->mark], ws->y[ws->mark]);
    } else {
      ws->scale = 1.0;
    }

    /* Set task to signal a new iterate.
       Set work to compute a new search direction. */

    /* Test for convergence. */
    df1 = fabs(f - ws->f0);
    df2 = ws->stp*fabs(ws->gd0);
    df = MAX(df1, df2);
    if (df <= ws->frtol*fabs(ws->f0)) {
      ws->task = OPK_VMLM_CONVERGENCE_FRTOL_TEST_SATISFIED;
    } else if (df <= ws->fatol) {
      ws->task = OPK_VMLM_CONVERGENCE_FATOL_TEST_SATISFIED;
    } else {
      ws->task = OPK_VMLM_NEW_X;
    }
    s_copy(work, "SEARCH DIRECTION");

    goto done;
  }
  if (s_cmp(work, "SEARCH DIRECTION") == 0) {
    /* Compute -H*g. */
    V_COPY(g, ws->x0);
    if ((mm = ws->iter - 1) > m) mm = m;
    opk_apply_lbfgs(vspace, mm, ws->s, ws->y, ws->rho,
                    ws->scale, ws->mark, ws->x0, ws->alpha);
    if (++ws->mark >= m) {
      ws->mark = 0;
    }
    V_COPY(ws->x0, ws->s[ws->mark]);
    /* Set task and work to initialize the line search. */
    s_copy(task, "START SEARCH");
    s_copy(work, "START SEARCH");
    goto L20;
  }
done:
  /* Save local variables. */
  if (s_cmp(work, "SEARCH") == 0) {
    ws->state = 1;
  }
  if (s_cmp(work, "SEARCH DIRECTION") == 0) {
    ws->state = 2;
  }
  return 0;
}

#endif

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
