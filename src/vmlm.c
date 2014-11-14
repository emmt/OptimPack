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
  return op->s[slot(op, k)];
}

opk_vector_t*
opk_get_lbfgs_y(opk_lbfgs_operator_t* op, opk_index_t k)
{
  return op->y[slot(op, k)];
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
opk_update_lbfgs_operator(opk_lbfgs_operator_t* op,
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
  double tiny;       /**< Relative size for a small step. */
  double grtol;      /**< Relative threshold for the norm or the gradient
                          (relative to GINIT the norm of the initial gradient)
                          for convergence. */
  double gatol;      /**< Absolute threshold for the norm or the gradient for
                          convergence. */
  double ginit;      /**< Euclidean norm or the initial gradient FIXME:
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

static void
finalize_vmlm(opk_object_t* obj)
{
  opk_vmlm_t* opt = (opk_vmlm_t*)obj;
  OPK_DROP(opt->vspace);
  OPK_DROP(opt->lnsrch);
  OPK_DROP(opt->B);
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
  opt->save_memory = TRUE;
#if 0
  opt->frtol = frtol;
  opt->fatol = fatol;
  opt->fmin = fmin;
#endif

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
opk_vmlm_iterate(opk_vmlm_t* opt, opk_vector_t* x1,
                 double f1, opk_vector_t* g1)
{
  double gtest, pg1;
  int status;

  if (opt->task == OPK_TASK_COMPUTE_FG) {
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
    gtest = opt->gatol + opt->grtol*opt->ginit;
    return optimizer_success(opt, ((opt->g1norm <= gtest)
                                   ? OPK_TASK_FINAL_X
                                   : OPK_TASK_NEW_X));
  }

  if (opt->task != OPK_TASK_NEW_X && opt->task != OPK_TASK_FINAL_X) {
    /* There must be something wrong. */
    return opt->task;
  }

  if (opt->task == OPK_TASK_NEW_X && opt->evaluations > 1) {
    /* Update the LBFGS matrix. */
    opk_update_lbfgs_operator(opt->B,
                              x1, opt->x0,
                              g1, opt->g0);
  }

  /* Compute a search direction.  We take care of checking whether D = -P is
     a sufficient descent direction.  As shown by Zoutendijk, this is true
     if: cos(theta) = -(D/|D|)'.(G/|G|) >= EPSILON > 0 where G is the
     gradient and D the descent direction. */
  while (TRUE) {
    opk_apply_direct((opk_operator_t*)opt->B, opt->p, g1);
    opt->pnorm = opk_vnorm2(opt->p);
    pg1 = opk_vdot(opt->p, g1);
    if (pg1 >= opt->epsilon*opt->pnorm*opt->g1norm) {
      /* Accept P (respectively D = -P) as a sufficient ascent (respectively
         descent) direction and set the directional derivative. */
      opt->dg0 = -pg1;
      break;
    }
    if (opt->B->mp < 1) {
      /* Initial iteration or recursion has just been restarted.  This means
         that the initial inverse Hessian approximation is not positive
         definite. */
      return optimizer_failure(opt, OPK_BAD_PRECONDITIONER);
    }
    /* Reset the LBFGS recursion and loop to use B0 to compute an initial
       search direction. */
    opk_reset_lbfgs_operator(opt->B);
    ++opt->restarts;
  }

  /* Save current variables X0, gradient G0 and function value F0.  The
     directional derivative DG0 = <D,G0> has already been set in the above
     loop.  */
  if (opt->save_memory) {
    /* Use the slot just after the mark to store X0 and G0.  Note that this
       is a weak reference: we do not "hold" the vectors. */
    opt->x0 = opk_get_lbfgs_s(opt->B, 1);
    opt->g0 = opk_get_lbfgs_y(opt->B, 1);
    if (opt->B->mp > opt->B->m - 1) {
      opt->B->mp = opt->B->m - 1;
    }
  }
  opk_vcopy(opt->x0, x1);
  opk_vcopy(opt->g0, g1);
  opt->g0norm = opt->g1norm;
  opt->f0 = f1;

  /* Estimate the length of the first step, start the line search and take
     the first step along the search direction. */
  if (opt->B->mp >= 1 || opt->B->rule == OPK_CUSTOM_APPROX) {
    opt->alpha = 1.0;
  } else if (0.0 < opt->tiny && opt->tiny < 1.0) {
    double x1norm = opk_vnorm2(x1);
    if (x1norm > 0.0) {
      opt->alpha = (x1norm/opt->g1norm)*opt->tiny;
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
