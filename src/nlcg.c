/*
 * nlcg.c --
 *
 * Non-linear conjugate gradient methods for OptimPack library.
 *
 * References:
 *
 * [1] Hestenes, M.R. & Stiefel, E., "Methods of Conjugate Gradients for
 *     Solving Linear Systems," Journal of Research of the National Bureau of
 *     Standards 49, 409-436 (1952).
 *
 * [2] Hager, W.W. & Zhang, H., "A survey of nonlinear conjugate gradient
 *     methods," Pacific Journal of Optimization, Vol. 2, pp. 35-58 (2006).
 *
 * [3] Hager, W. W. & Zhang, H. "A New Conjugate Gradient Method with
 *     Guaranteed Descent and an Efficient Line Search," SIAM J. Optim.,
 *     Vol. 16, pp. 170-192 (2005).
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (c) 2003-2014 Éric Thiébaut
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

/*---------------------------------------------------------------------------*/
/* PRIVATE DEFINITIONS */

/* Default parameters for the line search.  These values are from CG+ but
   other values may be more suitable. */
static const double STPMIN = 1E-20;
static const double STPMAX = 1E+20;
static const double SFTOL = 1E-4;
static const double SGTOL = 1E-1;
static const double SXTOL = DBL_EPSILON;

/* Default parameters for the global convergence. */
static const double GRTOL = 1E-6;
static const double GATOL = 0.0;

struct _opk_nlcg {
  opk_object_t base;     /* Base type (must be the first member). */
  double f0;             /* Function value at the start of the line search. */
  double g0norm;         /* Euclidean norm of G0, the gradient at the start of
                            the line search. */
  double g1norm;         /* Euclidean norm of G1, the gradient of the end of
                            the line search / last accepted point. */
  double dg0;            /* Directional derivative at the start of the line
                            search; given by the inner product:
                            <d,g0> = - <p,g0> */
  double dg1;            /* Directional derivative at the end or during the
                            line search; given by the inner product:
                            <d,g1> = - <p,g1> */
#if 0 /* not used */
  double alpha0;         /* Scale factor for the initial step step size. */
#endif
  double grtol;          /* Relative threshold for the norm or the gradient
                            (relative to the initial gradient) for
                            convergence. */
  double gatol;          /* Absolute threshold for the norm or the gradient for
                            convergence. */
  double ginit;          /* Euclidean norm of the initial gradient. */
  double fmin;           /* Minimal function value if provided. */
  double alpha;          /* Current step length. */
  double beta;           /* Current parameter in conjugate gradient update rule
                            (for information). */
  double stpmin;         /* Relative lower bound for the step length. */
  double stpmax;         /* Relative upper bound for the step length. */
  int (*update)(opk_nlcg_t* opt,
                const opk_vector_t* x1,
                const opk_vector_t* g1);
                         /* The update "method" is called to update the search
                            direction.  The returned value indicates whether
                            the updating rule has been successful, otherwise a
                            restart is needed. */
  opk_vspace_t* vspace;  /* Vector space of the variables of the problem. */
  opk_lnsrch_t* lnsrch;  /* Line search method. */
  opk_vector_t* x0;      /* Variables at start of line search. */
  opk_vector_t* g0;	 /* Gradient at start of line search. */
  opk_vector_t* p;	 /* (Anti-)search direction, new iterate is searched
                            as: x1 = x0 - alpha*p, for alpha >= 0. */
  opk_vector_t* y;       /* Work vector (e.g., to store the gradient
                            difference: Y = G1 - G0). */
  opk_index_t iter;      /* Iteration number. */
  opk_index_t nrestarts; /* Number of algorithm restarts. */
  opk_index_t nevals;    /* Number of function and gradient evaluations. */
  unsigned int method;   /* Conjugate gradient method. */
  opk_task_t task;       /* Current caller task. */
  opk_bool_t start;      /* Indicate whether algorithm is starting */
  opk_bool_t fmin_given; /* Indicate whether FMIN is specified. */
  opk_bool_t update_Hager_Zhang_orig;
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

/*
 * Most non-linear conjugate gradient methods, update the new search direction
 * by the following rule:
 *
 *     d' = -g1 + beta*d
 *
 * with d' the new search direction, g1 the current gradient, d the previous
 * search direction and beta a parameter which depends on the method.  For us,
 * the anti-search direction is:
 *
 *     p' = -d' = g1 + beta*p
 *
 * with p = -d the previous anti-search direction.  Some methods (e.g., Perry
 * & Shanno) implement the following rule:
 *
 *     d' = (-g1 + beta*d + gamma*y)*delta
 *
 * with y = g1 - g0.
 */

/* Helper function to compute search direction as: p' = g1 + beta*p. */
static int
update0(opk_nlcg_t* opt,
        const opk_vector_t* g1,
        double beta)
{
  if ((opt->beta = beta) != 0.0) {
    opk_vaxpby(opt->p, 1.0, g1, beta, opt->p);
    return OPK_SUCCESS;
  } else {
    return OPK_FAILURE;
  }
}

/* Helper function to compute search direction as: p' = g1 + beta*p
   possibly with the constraint that beta > 0. */
static int
update1(opk_nlcg_t* opt,
        const opk_vector_t* g1,
        double beta)
{
  if ((opt->method & OPK_NLCG_POWELL) != 0 && beta < 0.0) {
    ++opt->nrestarts;
    beta = 0.0;
  }
  if ((opt->beta = beta) != 0.0) {
    opk_vaxpby(opt->p, 1.0, g1, beta, opt->p);
    return OPK_SUCCESS;
  } else {
    return OPK_FAILURE;
  }
}

/* Form: Y = G1 - G0 */
static void
form_y(opk_nlcg_t* opt,
       const opk_vector_t* g1)
{
  opk_vaxpby(opt->y, 1.0, g1, -1.0, opt->g0);
}

/*
 * For Hestenes & Stiefel method:
 *
 *     beta = <g1,y>/<d,y> = - <g1,y>/<p,y>
 *
 * with y = g1 - g0.
 */
static int
update_Hestenes_Stiefel(opk_nlcg_t* opt,
                        const opk_vector_t* x1,
                        const opk_vector_t* g1)
{
  double g1y, dy, beta;
  form_y(opt, g1);
  g1y = opk_vdot(g1, opt->y);    /* Compute: g1y = <g1,y> */
  dy = -opk_vdot(opt->p, opt->y); /* Compute: dy = <d,y> = - <p,y> */
  beta = (dy != 0.0 ? g1y/dy : 0.0);
  return update1(opt, g1, beta);
}

/*
 * For Fletcher & Reeves method:
 *
 *     beta = <g1,g1>/<g0,g0>
 *
 * (this value is always >= 0 and can only be zero at a stationary point).
 */
static int
update_Fletcher_Reeves(opk_nlcg_t* opt,
                       const opk_vector_t* x1,
                       const opk_vector_t* g1)
{
  double r = opt->g1norm/opt->g0norm;
  return update0(opt, g1, r*r);
}

/*
 * For Polak-Ribière-Polyak method:
 *
 *     beta = <g1,y>/<g0,g0>
 */
static int
update_Polak_Ribiere_Polyak(opk_nlcg_t* opt,
                            const opk_vector_t* x1,
                            const opk_vector_t* g1)
{
  double beta;
  form_y(opt, g1);
  beta = opk_vdot(g1, opt->y)/(opt->g0norm*opt->g0norm);
  return update1(opt, g1, beta);
}

/*
 * For Fletcher "Conjugate Descent" method:
 *
 *     beta = <g1,g1>/(-<d,g0>)
 *
 * (this value is always >= 0 and can only be zero at a stationnary point).
 */
static int
update_Fletcher(opk_nlcg_t* opt,
                const opk_vector_t* x1,
                const opk_vector_t* g1)
{
  double beta = opt->g1norm*(opt->g1norm/(-opt->dg0));
  return update0(opt, g1, beta);
}

/*
 * For Liu & Storey method:
 *
 *     beta = <g1,y>/(-<d,g0>)
 */
static int
update_Liu_Storey(opk_nlcg_t* opt,
                  const opk_vector_t* x1,
                  const opk_vector_t* g1)
{
  double g1y, beta;
  form_y(opt, g1);
  g1y =  opk_vdot(g1, opt->y);    /* Compute: g1y = <g1,y> */
  beta = g1y/(-opt->dg0);
  return update1(opt, g1, beta);
}

/*
 * For Dai & Yuan method:
 *
 *     beta = <g1,g1>/<d,y>
 */
static int
update_Dai_Yuan(opk_nlcg_t* opt,
                const opk_vector_t* x1,
                const opk_vector_t* g1)
{
  double dy, beta;
  form_y(opt, g1);
  dy = -opk_vdot(opt->p, opt->y); /* Compute: dy = <d,y> = - <p,y> */
  beta = (dy != 0.0 ? opt->g1norm*(opt->g1norm/dy) : 0.0);
  return update1(opt, g1, beta);
}

/*
 * For Hager & Zhang method:
 *
 *     beta = <y - (2*<y,y>/<d,y>)*d,g1>
 *          = (<g1,y> - 2*<y,y>*<d,g1>/<d,y>)/<d,y>
 */
static int
update_Hager_Zhang(opk_nlcg_t* opt,
                   const opk_vector_t* x1,
                   const opk_vector_t* g1)
{
  double dy, beta;
  form_y(opt, g1);
  dy = -opk_vdot(opt->p, opt->y);
  if (dy != 0.0) {
    if (opt->update_Hager_Zhang_orig) {
      /* Original formulation. */
      double q = 1.0/dy;
      double r = q*opk_vnorm2(opt->y);
      opk_vaxpby(opt->y, q, opt->y, 2.0*r*r, opt->p);
      beta = opk_vdot(opt->y, g1);
    } else {
      /* Improved formulation which spares one linear combination and thus has
         less overhead (only 3 scalar products plus 2 linear combinations
         instead of 3 scalar products and 3 linear combinations).  The
         rounding errors are however different, so one or the other
         formulation can be by chance more efficient.  Though there is no
         systematic trend. */
      double yg = opk_vdot(opt->y, g1);
      double dg = opt->dg1;
      double r = opk_vnorm2(opt->y)/dy;
      beta = yg/dy - 2.0*r*r*dg;
    }
  } else {
    beta = 0.0;
  }
  return update1(opt, g1, beta);
}

/* Perry & Shanno, update rule (used in CONMIN and see Eq. (1.4) in [3])
 * writes:
 *
 *     d' = alpha*(-c1*g1 + c2*d - c3*y)  ==>   p' = c1*g1 + c2*p + c3*y
 *
 *     c1 = (1/alpha)*<s,y>/<y,y>
 *        =  <d,y>/<y,y>
 *        = -<p,y>/<y,y>
 *
 *     c2 = <g1,y>/<y,y> - 2*<s,g1>/<s,y>
 *	  = <g1,y>/<y,y> - 2*<d,g1>/<d,y>
 *	  = <g1,y>/<y,y> - 2*<p,g1>/<p,y>
 *
 *     c3 = -(1/alpha)*<s,g1>/<y,y>
 *        = -<d,g1>/<y,y>
 *        =  <p,g1>/<y,y>
 *
 * with alpha the step length, s = x1 - x0 = alpha*d = -alpha*p.  For this
 * method, beta = c2/c1.
 */
static int
update_Perry_Shanno(opk_nlcg_t* opt,
                    const opk_vector_t* x1,
                    const opk_vector_t* g1)
{
  double yy, dy, g1y, dg1, c1, c2, c3;
  form_y(opt, g1);
  yy = opk_vdot(opt->y, opt->y);
  if (yy <= 0.0) return OPK_FAILURE;
  dy = -opk_vdot(opt->p, opt->y);
  if (dy == 0.0) return OPK_FAILURE;
  g1y = opk_vdot(g1, opt->y);
  dg1 = opt->dg1;
  c1 = dy/yy;
  c2 = g1y/yy - 2.0*dg1/dy;
  c3 = -dg1/yy;
  opt->beta = c2/c1;
  opk_vaxpbypcz(opt->p, c1, g1, c2, opt->p, c3, opt->y);
  return OPK_SUCCESS;
}

static void
finalize_nlcg(opk_object_t* obj)
{
  opk_nlcg_t* opt = (opk_nlcg_t*)obj;
  OPK_DROP(opt->vspace);
  OPK_DROP(opt->lnsrch);
  OPK_DROP(opt->x0);
  OPK_DROP(opt->g0);
  OPK_DROP(opt->p);
  OPK_DROP(opt->y);
}

/*---------------------------------------------------------------------------*/
/* PUBLIC INTERFACE */

opk_nlcg_t*
opk_new_nlcg_optimizer_with_line_search(opk_vspace_t* vspace,
                                        unsigned int method,
                                        opk_lnsrch_t* lnsrch)
{
  opk_nlcg_t* opt;
  opk_bool_t g0_needed, y_needed;
  int (*update)(opk_nlcg_t* opt,
                const opk_vector_t* x1,
                const opk_vector_t* g1);

  /* Check the input arguments for errors. */
  if (vspace == NULL || lnsrch == NULL) {
    errno = EFAULT;
    return NULL;
  }
  if (method == 0) {
    method = OPK_NLCG_DEFAULT;
  }
  switch ((method & 0xff)) {
  case OPK_NLCG_FLETCHER_REEVES:
    update = update_Fletcher_Reeves;
    g0_needed = FALSE;
    y_needed = FALSE;
    break;
  case OPK_NLCG_HESTENES_STIEFEL:
    update = update_Hestenes_Stiefel;
    g0_needed = TRUE;
    y_needed = TRUE;
    break;
  case OPK_NLCG_POLAK_RIBIERE_POLYAK:
    update = update_Polak_Ribiere_Polyak;
    g0_needed = TRUE;
    y_needed = TRUE;
    break;
  case OPK_NLCG_FLETCHER:
    update = update_Fletcher;
    g0_needed = FALSE;
    y_needed = FALSE;
    break;
  case OPK_NLCG_LIU_STOREY:
    update = update_Liu_Storey;
    g0_needed = TRUE;
    y_needed = TRUE;
    break;
  case OPK_NLCG_DAI_YUAN:
    update = update_Dai_Yuan;
    g0_needed = TRUE;
    y_needed = TRUE;
    break;
  case OPK_NLCG_PERRY_SHANNO:
    update = update_Perry_Shanno;
    g0_needed = TRUE;
    y_needed = TRUE;
    break;
  case OPK_NLCG_HAGER_ZHANG:
    update = update_Hager_Zhang;
    g0_needed = TRUE;
    y_needed = TRUE;
    break;
  default:
    errno = EINVAL;
    return NULL;
  }

  /* We allocate enough memory for the workspace and instanciate it. */
  opt = (opk_nlcg_t*)opk_allocate_object(finalize_nlcg, sizeof(opk_nlcg_t));
  if (opt == NULL) {
    return NULL;
  }
  opt->update = update;
  opt->vspace = OPK_HOLD_VSPACE(vspace);
  opt->lnsrch = OPK_HOLD_LNSRCH(lnsrch);
  opt->fmin_given = FALSE;
  opt->method = method;
  opt->grtol = GRTOL;
  opt->gatol = GATOL;
  opt->stpmin = STPMIN;
  opt->stpmax = STPMAX;
  opt->x0 = opk_vcreate(vspace);
  if (opt->x0 == NULL) {
    goto error;
  }
  if (g0_needed) {
    opt->g0 = opk_vcreate(vspace);
    if (opt->g0 == NULL) {
      goto error;
    }
  }
  opt->p = opk_vcreate(vspace);
  if (opt->p == NULL) {
    goto error;
  }
  if (y_needed) {
    opt->y = opk_vcreate(vspace);
    if (opt->y == NULL) {
      goto error;
    }
  }
  opt->update_Hager_Zhang_orig = FALSE;
  opk_start_nlcg(opt);
  return opt;

 error:
  OPK_DROP(opt);
  return NULL;
}

opk_nlcg_t*
opk_new_nlcg_optimizer(opk_vspace_t* vspace, unsigned int method)
{
  opk_lnsrch_t* lnsrch;
  opk_nlcg_t* opt;

  lnsrch = opk_lnsrch_new_csrch(SFTOL, SGTOL, SXTOL);
  if (lnsrch == NULL) {
    return NULL;
  }
  opt = opk_new_nlcg_optimizer_with_line_search(vspace, method, lnsrch);
  OPK_DROP(lnsrch); /* the line search is now owned by the optimizer */
  return opt;
}

double
opk_get_nlcg_gatol(opk_nlcg_t* opt)
{
  return (opt == NULL ? GATOL : opt->gatol);
}

int
opk_set_nlcg_gatol(opk_nlcg_t* opt, double gatol)
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
opk_get_nlcg_grtol(opk_nlcg_t* opt)
{
  return (opt == NULL ? GRTOL : opt->grtol);
}

int
opk_set_nlcg_grtol(opk_nlcg_t* opt, double grtol)
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
opk_get_nlcg_stpmin(opk_nlcg_t* opt)
{
  return (opt == NULL ? STPMIN : opt->stpmin);
}

double
opk_get_nlcg_stpmax(opk_nlcg_t* opt)
{
  return (opt == NULL ? STPMAX : opt->stpmax);
}

int
opk_set_nlcg_stpmin_and_stpmax(opk_nlcg_t* opt, double stpmin, double stpmax)
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

int
opk_get_nlcg_fmin(opk_nlcg_t* opt, double* fmin)
{
  if (opt == NULL) {
    errno = EFAULT;
    return OPK_FAILURE;
  }
  if (opt->fmin_given) {
    if (fmin != NULL) {
      *fmin = opt->fmin;
    }
    return OPK_SUCCESS;
  } else {
    return OPK_FAILURE;
  }
}

int
opk_set_nlcg_fmin(opk_nlcg_t* opt, double fmin)
{
  if (opt == NULL) {
    errno = EFAULT;
    return OPK_FAILURE;
  }
  if (isnan(fmin) || isinf(fmin)) {
    errno = EINVAL;
    return OPK_FAILURE;
  }
  opt->fmin = fmin;
  opt->fmin_given = TRUE;
  return OPK_SUCCESS;
}

int
opk_unset_nlcg_fmin(opk_nlcg_t* opt)
{
  if (opt == NULL) {
    errno = EFAULT;
    return OPK_FAILURE;
  }
  opt->fmin_given = FALSE;
  return OPK_SUCCESS;
}

opk_index_t
opk_get_nlcg_iterations(opk_nlcg_t* opt)
{
  return opt->iter;
}

opk_index_t
opk_get_nlcg_restarts(opk_nlcg_t* opt)
{
  return opt->nrestarts;
}

opk_index_t
opk_get_nlcg_evaluations(opk_nlcg_t* opt)
{
  return opt->nevals;
}

unsigned int
opk_get_nlcg_method(opk_nlcg_t* opt)
{
  return opt->method;
}

opk_bool_t
opk_get_nlcg_starting(opk_nlcg_t* opt)
{
  return opt->start;
}

opk_task_t
opk_get_nlcg_task(opk_nlcg_t* opt)
{
  return opt->task;
}

double
opk_get_nlcg_alpha(opk_nlcg_t* opt)
{
  return opt->alpha;
}

double
opk_get_nlcg_beta(opk_nlcg_t* opt)
{
  return opt->beta;
}

opk_task_t
opk_start_nlcg(opk_nlcg_t* opt)
{
  opt->iter = 0;
  opt->nevals = 0;
  opt->nrestarts = 0;
  opt->task = OPK_TASK_COMPUTE_FG;
  opt->start = TRUE;
  return opt->task;
}

opk_task_t
opk_iterate_nlcg(opk_nlcg_t* opt, opk_vector_t* x1,
                 double f1, opk_vector_t* g1)
{
  /*
   * The new iterate is:
   *    x_{k+1} = x_{k} + \alpha_{k} d_{k}
   *            = x_{k} - \alpha_{k} p_{k}
   * as we consider the anti-search direction p = -d here.
   */
  int status = 0;

  if (opt->task == OPK_TASK_COMPUTE_FG) {
    opk_bool_t accept;
    ++opt->nevals;
    if (opt->start) {
      opt->g1norm = opk_vnorm2(g1);
      opt->ginit = opt->g1norm;
      accept = TRUE;
    } else {
      /* Compute directional derivative and check whether line search has
         converged. */
      opt->dg1 = -opk_vdot(opt->p, g1);
      status = opk_lnsrch_iterate(opt->lnsrch, &opt->alpha, f1, opt->dg1);
      if (status != 0) {
        /* Line search is finished, hopefully because of convergence
           (FIXME: check for warnings). */
        if (status < 0) {
          goto line_search_error;
        }
        ++opt->iter;
        opt->g1norm = opk_vnorm2(g1);
        accept = TRUE;
      } else {
        accept = FALSE;
      }
    }
    if (accept) {
      /* Check for global convergence. */
      if (opt->g1norm <= max3(0.0, opt->gatol, opt->grtol*opt->ginit)) {
        opt->task = OPK_TASK_FINAL_X;
      } else {
        opt->task = OPK_TASK_NEW_X;
      }
      return opt->task;
    }
  } else if (opt->task == OPK_TASK_NEW_X) {
    /* Compute a search direction and start line search. */
    int restart;
    if (opt->start) {
      restart = TRUE;
      opt->beta = 0.0;
    } else {
      restart = (opt->update(opt, x1, g1) != OPK_SUCCESS);
      if (! restart) {
        double dg = -opk_vdot(opt->p, g1);
        if (dg >= 0.0) {
          /* Restart if not a descent direction, not all updates warrant that.
             (FIXME: Generate an error instead?) */
          restart = TRUE;
        } else {
          /* Compute an initial step size ALPHA. */
          if ((opt->method & OPK_NLCG_SHANNO_PHUA) != 0) {
            /* Initial step size is such that:
               <alpha_{k+1}*d_{k+1},g_{k+1}> = <alpha_{k}*d_{k},g_{k}> */
            opt->alpha *= opt->dg0/dg;
#if 0 /* keep same ALPHA */
          } else {
            /* FIXME: check this */
            opt->alpha = (-opt->dg0)/p0norm/p0norm/opt->alpha;
#endif
          }
        }
        /* Save directional derivative. */
        opt->dg0 = dg;
      }
    }
    if (restart) {
      /* Initial search direction or recurrence has been restarted.  FIXME:
         other possibility is to use Fletcher's formula, see BGLS p. 39) */
      if (! opt->start) {
        ++opt->nrestarts;
      }
      opt->beta = 0.0;
#if 0
      double x1norm = opk_vnorm2(x1);
      if (x1norm > 0.0) {
        opt->alpha = (x1norm/opt->g1norm)*opt->tiny;
      } else {
        opt->alpha = (1e-3*MAX(fabs(f1), 1.0)/opt->g1norm)/opt->g1norm;
      }
#else
      opt->alpha = 1.0/opt->g1norm;
#endif
      opk_vcopy(opt->p, g1);
      opt->dg0 = -opt->g1norm*opt->g1norm;
    }

    /* Store current position as X0, f0, etc. */
    opk_vcopy(opt->x0, x1);
    opt->f0 = f1;
    if (opt->g0 != NULL) {
      opk_vcopy(opt->g0, g1);
    }
    opt->g0norm = opt->g1norm;

    /* Start the line search. */
    opt->dg1 = opt->dg0;
    status = opk_lnsrch_start(opt->lnsrch, opt->f0, opt->dg0, opt->alpha,
                              opt->stpmin*opt->alpha,
                              opt->stpmax*opt->alpha);
    if (status < 0) {
      goto line_search_error;
    }
  } else {
    return opt->task;
  }

  /* Build a new step to try. */
  opk_vaxpby(x1, 1.0, opt->x0, -opt->alpha, opt->p);
  opt->start = FALSE;
  opt->task = OPK_TASK_COMPUTE_FG;
  return opt->task;

  /* Some errors occurred in line search. */
 line_search_error:
  fprintf(stderr, "opk_lnsrch_start error: %d\n", status);
  /* FIXME: add line search error state */
  opt->task = OPK_TASK_ERROR;
  return opt->task;
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
