/*
 * vmlmc.c --
 *
 * Implementation of VMLMC, a variable metric limited memory optimization
 * algorithm with convex set contraints, for OptimPack library.
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

extern opk_vmlmc_t*
opk_new_vmlmc_optimizer(opk_vspace_t* vspace,
                        opk_index_t m,
                        double xsmall);

extern opk_vmlmc_t*
opk_new_vmlmc_optimizer_with_line_search(opk_vspace_t* vspace,
                                         opk_index_t m,
                                         double xsmall,
                                         opk_lnsrch_t* lnsrch);

extern opk_task_t
opk_start_vmlmc(opk_vmlmc_t* opt);


#define TRUE   OPK_TRUE
#define FALSE  OPK_FALSE

#define MAX(a,b)  OPK_MAX(a,b)
#define MIN(a,b)  OPK_MIN(a,b)

#define ROUND_UP(a,b)   OPK_ROUND_UP(a,b)

typedef enum {
  OPK_NO_PROBLEMS = 0,
  OPK_UNEXPECTED_TASK,     /**< Unexpected value for the pending optimizer
                                task.  This is most likely due to a misuse of
                                the optimizer. */
  OPK_BAD_DIRECTION,       /**< The opposite of the projected gradient is not a
                                descent direction or is not feasible.  This is
                                most likely due to a bug in the implementation
                                of the projector. */
  OPK_LINE_SEARCH_WARNING, /**< Warning in line search. */
  OPK_LINE_SEARCH_ERROR,   /**< Error in line search. */
} opk_reason_t;


/*---------------------------------------------------------------------------*/
/* OPTIMIZER WORKSPACE AND SETTINGS */

/* Default parameters for the line search.  These values are from the original
   SPG algorithm of Birgin et al. (2000) but other values may be more
   suitable. */
static const double STPMIN = 1E-30;
static const double STPMAX = 1E+30;
static const double SFTOL = 1E-4;
static const double SIGMA1 = 0.1;
static const double SIGMA2 = 0.2;

/* Default parameters for the global convergence. */
static const double GRTOL = 1E-6;
static const double GATOL = 0.0;

/* Other default parameters. */
static const int    SCALING = OPK_SCALING_OREN_SPEDICATO;
static const double STPSIZ = 1.0;

struct _opk_vmlmc {
  opk_object_t base;       /**< Base type (must be the first member). */
  double gamma;            /**< Scale factor to approximate inverse Hessian. */
  double grtol;            /**< Relative threshold for the norm or the gradient
                                (relative to PGINIT the norm of the initial
                                projected gradient) for convergence. */
  double gatol;            /**< Absolute threshold for the norm or the gradient
                                for convergence. */
  double pginit;           /**< Euclidean norm or the initial projected
                                gradient. */
  double pgnorm;           /**< Euclidean norm of the projected gradient at
                                the last accepted step. */
  double f0;               /**< Function value at x0. */
  double stp;              /**< Current step length. */
  double stpmin;           /**< Maximum relative step length. */
  double stpmax;           /**< Minimum relative step length. */
  double stpsiz;           /**< Euclidean norm of a small change of variables,
                                must be strictly positive.  The first step at
                                start or after a restart is the steepest
                                descent scaled to have this length. */
  opk_vspace_t* vspace;    /**< Vector space for variables of the problem. */
  opk_lnsrch_t* lnsrch;    /**< Line search method. */
  opk_vector_t* x0;        /**< Variables at the start of the line search. */
  opk_vector_t* g0;        /**< Gradient at x0. */
  opk_vector_t* pg;        /**< Projected gradient. */
  opk_vector_t* h0;        /**< Diagonal preconditioner or NULL. */
  opk_vector_t** s;        /**< Storage for variable differences. */
  opk_vector_t** y;        /**< Storage for gradient differences. */
  double* alpha;           /**< Workspace to save 1/<s[j],y[j]> */
  double* rho;             /**< Workspace to save 1/<s[j],y[j]> */
  opk_index_t m;           /**< Maximum number of memorized steps. */
  opk_index_t mp;          /**< Actual number of memorized steps
                                (0 <= mp <= m). */
  opk_index_t mark;        /**< Index of oldest saved step. */
  opk_index_t evaluations; /**< Number of functions (and gradients)
                                evaluations. */
  opk_index_t iterations;  /**< Number of iterations (successful steps
                                taken). */
  opk_index_t restarts;    /**< Number of restarts. */
  opk_index_t projections; /** Number of projections. */
  opk_task_t task;         /**< Pending task. */
  int reason;              /**< Some details about the error, if any. */
  int scaling;             /**< Method to compute GAMMA. */
  int stage;               /**< STAGE = 0 for the starting solution, STAGE = 1
                                for the first line search step, STAGE = 2 for
                                the other steps. */
  opk_bool_t save_memory;  /**< To save space, the variable and gradient at the
                                start of a line search are weak references to
                                the (s,y) pair of vectors of the LBFGS operator
                                just after the mark. */
};

/*---------------------------------------------------------------------------*/
/* L-BFGS OPERATOR */

static void
lbfgs_reset(opk_vmlmc_t* opt)
{
  opt->mp = 0;
}

/* Define a a bunch of macros to make the code easier to read. */
#define DOT(x,y)              opk_vdot(x, y)
#define DOT3(w,x,y)           opk_vdot3(w, x, y)
#define S(j)                  opt->s[j]
#define Y(j)                  opt->y[j]
#define RHO(j)                opt->rho[j]
#define ALPHA(j)              opt->alpha[j]
#define AXPBY(dst,a,x,b,y)    opk_vaxpby(dst, a, x, b, y)
#define SCALE(dst,alpha,src)  opk_vscale(dst, alpha, src)

/* Get index of k-th (s,y) pair memorized by L-BFGS operator.  Argument k is
   relative to the last saved pair (which corresponds to k = 0); thus k = -1
   for the pair just saved before the last one, and k = 1 for the slot just
   after the last one.  Argument k must be in the range 1 - MP <= k <= 1, with
   MP <= M the actual number of saved pairs. */
static opk_index_t
lbfgs_slot(opk_vmlmc_t* opt, opk_index_t k)
{
  return (opt->m + opt->mark + k)%opt->m;
}

/* Apply L-BFGS operator in place.  The operation is split in two stages
   corresponding to the first and second loop of the two-loop recursion
   described by Nocedal (1980) and due to Strang. */

static void
lbfgs_loop1(opk_vmlmc_t* opt, opk_vector_t* d)
{
  opk_index_t j, k;

  /* First loop of the recursion from the newest saved (s,y) pair to the
     oldest one. */
  for (k = 0; k < opt->mp; ++k) {
    j = lbfgs_slot(opt, -k);
    if (RHO(j) > 0.0) {
      ALPHA(j) = RHO(j)*DOT(S(j), d);
      AXPBY(d, 1.0, d, -ALPHA(j), Y(j));
    } else {
      ALPHA(j) = 0.0;
    }
  }

  /* Apply scaling. */
  if (opt->gamma > 0.0 && opt->gamma != 1.0) {
    SCALE(d, opt->gamma, d);
  }
}

static void
lbfgs_loop2(opk_vmlmc_t* opt, opk_vector_t* d)
{
  double beta;
  opk_index_t j, k;

  /* Second loop of the recursion from the oldest saved (s,y) pair to the
     newest one. */
  for (k = opt->mp - 1; k >= 0; --k) {
    j = lbfgs_slot(opt, -k);
    if (RHO(j) > 0.0) {
      beta = RHO(j)*DOT(d, Y(j));
      opk_vaxpby(d, 1.0, d, ALPHA(j) - beta, S(j));
    }
  }
}

static void
lbfgs_update(opk_vmlmc_t* opt,
             const opk_vector_t* x1,
             const opk_vector_t* x0,
             const opk_vector_t* g1,
             const opk_vector_t* g0)
{
  double sty;
  opk_index_t j;

  /* Store the variables and gradient differences in the slot just after the
     mark (which is the index of the last saved pair or -1 if none). */
  j = lbfgs_slot(opt, 1);
  AXPBY(S(j), 1.0, x1, -1.0, x0);
  AXPBY(Y(j), 1.0, g1, -1.0, g0);

  /* Compute RHO[j] and GAMMA.  If the update formula for GAMMA does not yield
     a strictly positive value, the strategy is to keep the previous value and
     discard the (s,y) pair. */
  sty = DOT(S(j), Y(j));
  if (sty <= 0.0) {
    /* This pair will be skipped. */
    RHO(j) = 0.0;
    if (opt->mp == opt->m) {
      --opt->mp;
    }
    return;
  }
  RHO(j) = 1.0/sty;
  if (opt->scaling == OPK_SCALING_OREN_SPEDICATO) {
    opt->gamma = sty/DOT(Y(j), Y(j));
  } else if (opt->scaling == OPK_SCALING_BARZILAI_BORWEIN) {
    opt->gamma = DOT(S(j), S(j))/sty;
  }

  /* Update the mark and the number of saved pairs. */
  opt->mark = j;
  if (opt->mp < opt->m) {
    ++opt->mp;
  }
}

#undef DOT
#undef DOT3
#undef S
#undef Y
#undef RHO
#undef ALPHA
#undef AXPBY
#undef SCALE

/*---------------------------------------------------------------------------*/

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
finalize_vmlmc(opk_object_t* obj)
{
  opk_vmlmc_t* opt = (opk_vmlmc_t*)obj;
  opk_index_t k;

  OPK_DROP(opt->vspace);
  OPK_DROP(opt->lnsrch);
  if (! opt->save_memory) {
    OPK_DROP(opt->x0);
    OPK_DROP(opt->g0);
  }
  OPK_DROP(opt->pg);
  OPK_DROP(opt->h0);

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

opk_vmlmc_t*
opk_new_vmlmc_optimizer_with_line_search(opk_vspace_t* vspace,
                                         opk_index_t m,
                                         double stpsiz,
                                         opk_lnsrch_t* lnsrch)
{
  opk_vmlmc_t* opt;
  size_t s_offset, y_offset, alpha_offset, rho_offset, size;
  opk_index_t k;

  /* Check the input arguments for errors. */
  if (vspace == NULL || lnsrch == NULL) {
    errno = EFAULT;
    return NULL;
  }
  if (vspace->size < 1 || m < 1 || non_finite(stpsiz) || stpsiz <= 0.0) {
    errno = EINVAL;
    return NULL;
  }

  /* Allocate and instanciate the workspace (not the part which is done by
     opk_start_vmlmc). */
  s_offset = OPK_ROUND_UP(sizeof(opk_vmlmc_t), sizeof(opk_vector_t*));
  y_offset = s_offset + m*sizeof(opk_vector_t*);
  alpha_offset = OPK_ROUND_UP(y_offset + m*sizeof(opk_vector_t*),
                              sizeof(double));
  rho_offset = alpha_offset + m*sizeof(double);
  size = rho_offset + m*sizeof(double);
  opt = (opk_vmlmc_t*)opk_allocate_object(finalize_vmlmc, size);
  if (opt == NULL) {
    return NULL;
  }
  opt->s = (opk_vector_t**)(((unsigned char*)opt) + s_offset);
  opt->y = (opk_vector_t**)(((unsigned char*)opt) + y_offset);
  opt->alpha =    (double*)(((unsigned char*)opt) + alpha_offset);
  opt->rho =      (double*)(((unsigned char*)opt) + rho_offset);
  opt->gamma = 1.0;
  opt->grtol = GRTOL;
  opt->gatol = GATOL;
  opt->stpmin = STPMIN;
  opt->stpmax = STPMAX;
  opt->stpsiz = stpsiz;
  opt->vspace = OPK_HOLD_VSPACE(vspace);
  opt->lnsrch = OPK_HOLD_LNSRCH(lnsrch);
  opt->m = m;
  opt->scaling = SCALING;
  opt->save_memory = TRUE;

  /* Allocate work vectors.  If saving memory, x0 and g0 will be weak
     references to one of the saved vectors in the LBFGS operator. */
  if (! opt->save_memory) {
    opt->x0 = opk_vcreate(vspace);
    if (opt->x0 == NULL) {
      goto error;
    }
    opt->g0 = opk_vcreate(vspace);
    if (opt->g0 == NULL) {
      goto error;
    }
  }
  opt->pg = opk_vcreate(vspace);
  if (opt->pg == NULL) {
    goto error;
  }
  for (k = 0; k < opt->m; ++k) {
    opt->s[k] = opk_vcreate(vspace);
    if (opt->s[k] == NULL) {
      goto error;
    }
    opt->y[k] = opk_vcreate(vspace);
    if (opt->y[k] == NULL) {
      goto error;
    }
  }
  (void)opk_start_vmlmc(opt);
  return opt;

 error:
  OPK_DROP(opt);
  return NULL;
}

opk_vmlmc_t*
opk_new_vmlmc_optimizer(opk_vspace_t* vspace,
                        opk_index_t m,
                        double stpsiz)
{
  opk_lnsrch_t* lnsrch;
  opk_vmlmc_t* opt;

#if 0
  /* Create nonmonotone line search with same parameters as in Birgin et
     al. (2000) but with a memory of 1 to mimic monotone line search with
     backtracking and quadratic interpolation. */
  lnsrch = opk_lnsrch_new_nonmonotone(1, SFTOL, SIGMA1, SIGMA2);
#else
  lnsrch = opk_lnsrch_new_csrch(1E-4, 0.9, 2E-17);
#endif
  if (lnsrch == NULL) {
    return NULL;  }
  opt = opk_new_vmlmc_optimizer_with_line_search(vspace, m, stpsiz, lnsrch);
  OPK_DROP(lnsrch); /* the line search is now owned by the optimizer */
  return opt;
}

static opk_task_t
optimizer_success(opk_vmlmc_t* opt, opk_task_t task)
{
  opt->reason = OPK_NO_PROBLEMS;
  opt->task = task;
  return task;
}

static opk_task_t
optimizer_failure(opk_vmlmc_t* opt, opk_reason_t reason)
{
  opt->reason = reason;
  opt->task = OPK_TASK_ERROR;
  return opt->task;
}

/* Make first line search step. */
static opk_task_t
first_step(opk_vmlmc_t* opt, opk_vector_t* x,
           double stp, const opk_vector_t* d)
{
  opt->stp = fabs(stp);
  opk_vaxpby(x, 1.0, opt->x0, stp, d);
  opt->stage = 1;
  return optimizer_success(opt, OPK_TASK_PROJECT_X);
}

static opk_task_t
line_search_failure(opk_vmlmc_t* opt)
{
  opk_reason_t reason;
  if (opk_lnsrch_has_errors(opt->lnsrch)) {
    reason = OPK_LINE_SEARCH_ERROR;
  } else {
    reason = OPK_LINE_SEARCH_WARNING;
  }
  return optimizer_failure(opt, reason);
}

opk_task_t
opk_start_vmlmc(opk_vmlmc_t* opt)
{
  lbfgs_reset(opt);
  opt->evaluations = 0;
  opt->iterations = 0;
  opt->restarts = 0;
  opt->projections = 0;
  opt->stp = 0.0;
  opt->stage = 0;
  return optimizer_success(opt, OPK_TASK_PROJECT_X);
}

static void
save_iterate(opk_vmlmc_t* opt, const opk_vector_t* x,
             double f, const opk_vector_t* g)
{
  /* Save current variables X0, gradient G0 and function value F0. */
  if (opt->save_memory) {
    /* Use the slot just after the mark in the LBFGS operator to store X0 and
       G0.  Note that this is a "weak" reference: we do not "hold" the
       vectors. */
    opk_index_t j = lbfgs_slot(opt, 1);
    opt->x0 = opt->s[j];
    opt->g0 = opt->y[j];
    if (opt->mp > opt->m - 1) {
      opt->mp = opt->m - 1;
    }
  }
  opk_vcopy(opt->x0, x);
  opk_vcopy(opt->g0, g);
  opt->f0 = f;
}

/*
 * During line search, the k-th iterate writes:
 *
 *     x_{k} = x_{0} + stp_{k} d
 *
 * with x_{0} the point at the start of the line search, d the search direction,
 * and stp_{k} the step length (equal to zero for k=0 and to one for k=1).
 * The search direction is thus:
 *
 *     d = x_{1} - x_{0}
 *
 * with x_{1} the first tried point (after projection).  It is then possible to
 * rewrite the k-th iterate as:
 *
 *     x_{k} = x_{0} + (stp_{k}/stp_{k-1}) (x_{k-1} - x_{0})
 *           = (1 - phi_{k}) x_{0} + phi_{k} x_{k-1}
 *
 * with:
 *
 *     phi_{k} = stp_{k}/stp_{k-1}
 *
 * this is an economical way to update the step because it does not require to
 * save the search direction which is only explicitely needed for the second
 * Wolfe condition.  Thus this formula may be considered when only Armijo-like
 * condition is used to accept a step.  Another advantage is that if the
 * variables x have to be restricted to a convex feasible set, then x_{k} is
 * automatically feasible if x_{0} and x_{k-1} are feasible and if 0 <= phi_{k}
 * <= 1 which is the case during backtracking (moreover, in this case, phi_k is
 * a constant, e.g. 1/2).  Using this formula may however be more subject to
 * the accumulation of rounding errors.
 */
opk_task_t
opk_iterate_vmlmc(opk_vmlmc_t* opt, opk_vector_t* x, double f,
                  opk_vector_t* g, opk_vector_t* d)
{
  int status;

  switch (opt->task) {

  case OPK_TASK_PROJECT_X:

    /* Caller has projected the variables x to the feasible set. */
    ++opt->projections;

    if (opt->stage == 1) {
      /* This is the first step along a new search direction.  Form actual
         search direction and check whether is is a descent direction.  See
         Nocedal & Wright ("Numerical Optimization", 2006) for a justification
         that just the sign of the directional derivative has to be checked
         (i.e. not a threshold). */
      double dg0; /* FIXME: initial directional derivative */
      opk_vaxpby(d, 1.0, x, -1.0, opt->x0);
      dg0 = opk_vdot(d, opt->g0);
      if (dg0 >= 0.0) {
        /* Initial step is not along a descent direction.  If the search
           direction has been produced by the L-BFGS recursion, it means that
           this approximation is not positive definite in the local sub-space
           of feasible directions.  In that case, we restart L-BFGS recursion
           and manage to use the projected gradient as the next search
           direction. */
        if (opt->mp == 0) {
          /* We were already using the steepest descent projected direction.
             Thus there must be an error. */
          fprintf(stderr, "bad steepest descent projected direction\n");
          return optimizer_failure(opt, OPK_BAD_DIRECTION);
        }
        lbfgs_reset(opt);
        ++opt->restarts;
        return first_step(opt, x, -(opt->stpsiz/opt->pgnorm), opt->pg);
      } else {
        /* Initial step is along a descent direction.  Initialize line search
           with safeguard bounds to only allow for bracktracking. */
        opt->stp = 1.0;
        status = opk_lnsrch_start(opt->lnsrch, opt->f0, dg0, opt->stp,
                                  opt->stpmin*opt->stp, opt->stp);
        if (status != OPK_LNSRCH_SEARCH) {
          return line_search_failure(opt);
        }
        opt->stage = 2;
      }
    }

    /* Request caller to compute objective function and gradient at x. */
    return optimizer_success(opt, OPK_TASK_COMPUTE_FG);


  case OPK_TASK_COMPUTE_FG:

    /* Caller has computed the function value and the gradient at the current
       iterate. */
    ++opt->evaluations;
    if (opt->evaluations == 1) {
      /* Save the initial iterate. */
      save_iterate(opt, x, f, g);
    }

    if (opt->stage == 2) {
      /* A line search is in progress.  Compute directional derivative and check
         whether line search has converged. */
      double dg; /* FIXME: directional derivative */
      if (opk_lnsrch_use_deriv(opt->lnsrch)) {
        dg = opk_vdot(d, g);
      } else {
        dg = 0.0;
      }
      switch (opk_lnsrch_iterate(opt->lnsrch, &opt->stp, f, dg)) {
      case OPK_LNSRCH_SEARCH:
        /* Line search has not yet converged.  Compute the new iterate to try
           and request caller to compute function value and gradient.  (There
           are no needs to project the variables as the new iterate is
           guaranteed to be feasible.) */
        opk_vaxpby(x, 1.0, opt->x0, opt->stp, d);
        return optimizer_success(opt, OPK_TASK_COMPUTE_FG);
      case OPK_LNSRCH_CONVERGENCE:
      case OPK_LNSRCH_WARNING_ROUNDING_ERRORS_PREVENT_PROGRESS:
      case OPK_LNSRCH_WARNING_STP_EQ_STPMAX:
        /* Line search has converged. */
        ++opt->iterations;
        break;
      default:
        /* There are errors. */
        return line_search_failure(opt);
      }
    }

    /* Starting solution or line search has converged.  We require the caller
       to compute the projected gradient which is needed to check for global
       convergence and to compute the next search direction. */
    opk_vcopy(d, g);
    return optimizer_success(opt, OPK_TASK_PROJECT_D);

  case OPK_TASK_PROJECT_D:

    ++opt->projections;
    if (opt->stage < 3) {
      /* The vector d contains the projected gradient.  Check for global
         convergence. */
      opk_task_t next_task;
      opk_vcopy(opt->pg, d); /* save the projected gradient */
      opt->pgnorm = opk_vnorm2(opt->pg);
      if (opt->evaluations == 1) {
        opt->pginit = opt->pgnorm;
      }
      if (opt->pgnorm <= max3(0.0, opt->gatol, opt->grtol*opt->pginit)) {
        next_task = OPK_TASK_FINAL_X;
      } else {
        next_task = OPK_TASK_NEW_X;
      }
      return optimizer_success(opt, next_task);
    } else {
      /* The caller has projected the vector d resulting from the first loop of
         the L-BFGS recursion.  Apply the 2nd loop of L-BFGS recursion and
         compute the first iterate of the next line search. */
      lbfgs_loop2(opt, d);
      save_iterate(opt, x, f, g);
      return first_step(opt, x, -1.0, d);
    }

  case OPK_TASK_FINAL_X:
  case OPK_TASK_NEW_X:

    if (opt->stage > 0) {
      /* Update L-BFGS operator and apply the first loop of the L-BFGS
         recursion to vector d which already contains the projected
         gradient. */
      lbfgs_update(opt, x, opt->x0, g, opt->g0);
      opk_vcopy(d, g);
      lbfgs_loop1(opt, d);
      opt->stage = 3;
      return optimizer_success(opt, OPK_TASK_PROJECT_D);
    } else {
      /* Initial iterate or algorithm has been restarted.  The search direction
         is the scaled projected gradient. */
      return first_step(opt, x, -(opt->stpsiz/opt->pgnorm), opt->pg);
    }

  default:

    /* There must be something wrong. */
    return optimizer_failure(opt, OPK_UNEXPECTED_TASK);

  }

}

opk_task_t
opk_get_vmlmc_task(opk_vmlmc_t* opt)
{
  return opt->task;
}

opk_index_t
opk_get_vmlmc_iterations(opk_vmlmc_t* opt)
{
  return opt->iterations;
}

opk_index_t
opk_get_vmlmc_evaluations(opk_vmlmc_t* opt)
{
  return opt->evaluations;
}

opk_index_t
opk_get_vmlmc_restarts(opk_vmlmc_t* opt)
{
  return opt->restarts;
}

opk_index_t
opk_get_vmlmc_projections(opk_vmlmc_t* opt)
{
  return opt->projections;
}

int
opk_get_vmlmc_scaling(opk_vmlmc_t* opt)
{
  return (opt == NULL ? SCALING : opt->scaling);
}

int
opk_set_vmlmc_scaling(opk_vmlmc_t* opt, int scaling)
{
  if (opt == NULL) {
    errno = EFAULT;
    return OPK_FAILURE;
  }
  switch (scaling) {
  case OPK_SCALING_NONE:
  case OPK_SCALING_OREN_SPEDICATO:
  case OPK_SCALING_BARZILAI_BORWEIN:
    opt->scaling = scaling;
    return OPK_SUCCESS;
  default:
    errno = EINVAL;
    return OPK_FAILURE;
  }
}

double
opk_get_vmlmc_gatol(opk_vmlmc_t* opt)
{
  return (opt == NULL ? GATOL : opt->gatol);
}

int
opk_set_vmlmc_gatol(opk_vmlmc_t* opt, double gatol)
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
opk_get_vmlmc_grtol(opk_vmlmc_t* opt)
{
  return (opt == NULL ? GRTOL : opt->grtol);
}

int
opk_set_vmlmc_grtol(opk_vmlmc_t* opt, double grtol)
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
opk_get_vmlmc_stpmin(opk_vmlmc_t* opt)
{
  return (opt == NULL ? STPMIN : opt->stpmin);
}

double
opk_get_vmlmc_stpmax(opk_vmlmc_t* opt)
{
  return (opt == NULL ? STPMAX : opt->stpmax);
}

int
opk_set_vmlmc_stpmin_and_stpmax(opk_vmlmc_t* opt, double stpmin, double stpmax)
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
