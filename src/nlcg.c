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
 * This file is part of OptimPack (https://github.com/emmt/OptimPack).
 *
 * Copyright (C) 2003-2015 Éric Thiébaut
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

/* Code returned by the `update` procedures. */
#define UPDATE_SUCCESS   0
#define UPDATE_FAILURE  -1

/*---------------------------------------------------------------------------*/
/* PRIVATE DEFINITIONS */

/* Default parameters for the line search.  These values are from CG+ but
   other values may be more suitable. */
static const double STPMIN = 1E-20;
static const double STPMAX = 1E+20;
static const double SFTOL = 0.0001;
static const double SGTOL = 0.1;
static const double SXTOL = DBL_EPSILON;

/* Default parameters for the global convergence. */
static const double GRTOL = 1E-6;
static const double GATOL = 0.0;

/* Other default parameters. */
static const double DELTA   = 5e-2;
static const double EPSILON = 0.0;

struct _opk_nlcg {
  opk_object_t base;      /* Base type (must be the first member). */
  double f0;              /* Function value at the start of the line search. */
  double g0norm;          /* Euclidean norm of G0, the gradient at the start of
                             the line search. */
  double gnorm;           /* Euclidean norm of G, the gradient of the last
                             accepted point. */
  double dtg0;            /* Directional derivative at the start of the line
                             search; given by the inner product: -<d,g0> */
  double dtg;             /* Directional derivative at the last trial point;
                             given by the inner product: -<d,g> */
  double grtol;           /* Relative threshold for the norm or the gradient
                             (relative to the initial gradient) for
                             convergence. */
  double gatol;           /* Absolute threshold for the norm or the gradient
                             for convergence. */
  double ginit;           /* Euclidean norm of the initial gradient. */
  double fmin;            /* Minimal function value if provided. */
  double delta;           /* Relative size for a small step. */
  double epsilon;         /* Threshold to accept descent direction. */
  double alpha;           /* Current step length. */
  double beta;            /* Current parameter in conjugate gradient update
                             rule (for information). */
  double stpmin;          /* Relative lower bound for the step length. */
  double stpmax;          /* Relative upper bound for the step length. */
  int (*update)(opk_nlcg_t* opt,
                const opk_vector_t* x,
                const opk_vector_t* g);
                          /* The update "method" is called to update the search
                             direction.  The returned value indicates whether
                             the updating rule has been successful, otherwise a
                             restart is needed. */
  opk_vspace_t* vspace;   /* Vector space of the variables of the problem. */
  opk_lnsrch_t* lnsrch;   /* Line search method. */
  opk_vector_t* x0;       /* Variables at start of line search. */
  opk_vector_t* g0;       /* Gradient at start of line search. */
  opk_vector_t* d;        /* (Anti-)search direction, new iterate is searched
                             as: x = x0 - alpha*d, for alpha >= 0. */
  opk_vector_t* y;        /* Work vector (e.g., to store the gradient
                             difference: Y = G - G0). */
  opk_index_t iterations; /* Number of iterations. */
  opk_index_t restarts;   /* Number of algorithm restarts. */
  opk_index_t evaluations;/* Number of function and gradient evaluations. */
  unsigned int flags;     /* Conjugate gradient method and options. */
  opk_status_t status;    /* Current error status. */
  opk_task_t task;        /* Pending caller task. */
  opk_bool_t fmin_given;  /* Indicate whether FMIN is specified. */
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

static opk_task_t
success(opk_nlcg_t* opt, opk_task_t task)
{
  opt->status = OPK_SUCCESS;
  opt->task = task;
  return opt->task;
}

static opk_task_t
failure(opk_nlcg_t* opt, opk_status_t status)
{
  opt->status = status;
  opt->task = OPK_TASK_ERROR;
  return opt->task;
}

/*
 * Most non-linear conjugate gradient methods, update the new search direction
 * by the following rule:
 *
 *     d' = -g + beta*d
 *
 * with d' the new search direction, g the current gradient, d the previous
 * search direction and beta a parameter which depends on the method.
 *
 * Some methods (e.g., Perry & Shanno) implement the following rule:
 *
 *     d' = (-g + beta*d + gamma*y)*delta
 *
 * with y = g - g0.
 *
 * For us, the anti-search direction is used instead, thus:
 *
 *     d' = g + beta*d
 */

/* Helper function to compute search direction as: d' = g + beta*d. */
static int
update0(opk_nlcg_t* opt,
        const opk_vector_t* g,
        double beta)
{
  if ((opt->beta = beta) == 0) {
    return UPDATE_FAILURE;
  }
  opk_vaxpby(opt->d, 1, g, beta, opt->d);
  return UPDATE_SUCCESS;
}

/* Helper function to compute search direction as: d' = g + beta*d
   possibly with the constraint that beta > 0. */
static int
update1(opk_nlcg_t* opt,
        const opk_vector_t* g,
        double beta)
{
  if ((opt->flags & OPK_NLCG_POWELL) == OPK_NLCG_POWELL && beta < 0) {
    beta = 0;
  }
  if ((opt->beta = beta) == 0) {
    return UPDATE_FAILURE;
  }
  opk_vaxpby(opt->d, 1, g, beta, opt->d);
  return UPDATE_SUCCESS;
}

/* Form: Y = G - G0 */
static void
form_y(opk_nlcg_t* opt,
       const opk_vector_t* g)
{
  opk_vaxpby(opt->y, 1, g, -1, opt->g0);
}

/*
 * For Hestenes & Stiefel method:
 *
 *     beta = <g,y>/<d,y>
 *
 * with y = g - g0.
 */
static int
update_Hestenes_Stiefel(opk_nlcg_t* opt,
                        const opk_vector_t* x,
                        const opk_vector_t* g)
{
  double gty, dty, beta;
  form_y(opt, g);
  gty =  opk_vdot(g, opt->y);
  dty = -opk_vdot(opt->d, opt->y);
  beta = (dty != 0.0 ? gty/dty : 0.0);
  return update1(opt, g, beta);
}

/*
 * For Fletcher & Reeves method:
 *
 *     beta = <g,g>/<g0,g0>
 *
 * (this value is always >= 0 and can only be zero at a stationary point).
 */
static int
update_Fletcher_Reeves(opk_nlcg_t* opt,
                       const opk_vector_t* x,
                       const opk_vector_t* g)
{
  double r = opt->gnorm/opt->g0norm;
  return update0(opt, g, r*r);
}

/*
 * For Polak-Ribière-Polyak method:
 *
 *     beta = <g,y>/<g0,g0>
 */
static int
update_Polak_Ribiere_Polyak(opk_nlcg_t* opt,
                            const opk_vector_t* x,
                            const opk_vector_t* g)
{
  double beta;
  form_y(opt, g);
  beta = (opk_vdot(g, opt->y)/opt->g0norm)/opt->g0norm;
  return update1(opt, g, beta);
}

/*
 * For Fletcher "Conjugate Descent" method:
 *
 *     beta = -<g,g>/<d,g0>
 *
 * (this value is always >= 0 and can only be zero at a stationnary point).
 */
static int
update_Fletcher(opk_nlcg_t* opt,
                const opk_vector_t* x,
                const opk_vector_t* g)
{
  double beta = -opt->gnorm*(opt->gnorm/opt->dtg0);
  return update0(opt, g, beta);
}

/*
 * For Liu & Storey method:
 *
 *     beta = -<g,y>/<d,g0>
 */
static int
update_Liu_Storey(opk_nlcg_t* opt,
                  const opk_vector_t* x,
                  const opk_vector_t* g)
{
  double gty, beta;
  form_y(opt, g);
  gty =  opk_vdot(g, opt->y);
  beta = -gty/opt->dtg0;
  return update1(opt, g, beta);
}

/*
 * For Dai & Yuan method:
 *
 *     beta = <g,g>/<d,y>
 */
static int
update_Dai_Yuan(opk_nlcg_t* opt,
                const opk_vector_t* x,
                const opk_vector_t* g)
{
  double dty, beta;
  form_y(opt, g);
  dty = -opk_vdot(opt->d, opt->y);
  beta = (dty != 0.0 ? opt->gnorm*(opt->gnorm/dty) : 0.0);
  return update1(opt, g, beta);
}

/*
 * For Hager & Zhang method:
 *
 *     beta = <y - (2*<y,y>/<d,y>)*d,g>/<d,y>
 *          = (<g,y> - 2*<y,y>*<d,g>/<d,y>)/<d,y>
 */
static int
update_Hager_Zhang(opk_nlcg_t* opt,
                   const opk_vector_t* x,
                   const opk_vector_t* g)
{
  double dty, beta;
  form_y(opt, g);
  dty = -opk_vdot(opt->d, opt->y);
  if (dty != 0) {
    if (opt->update_Hager_Zhang_orig) {
      /* Original formulation, using Y as a scratch vector. */
      double q = 1.0/dty;
      double r = q*opk_vnorm2(opt->y);
      opk_vaxpby(opt->y, q, opt->y, 2.0*r*r, opt->d);
      beta = opk_vdot(opt->y, g);
    } else {
      /* Improved formulation which spares one linear combination and thus has
         less overhead (only 3 scalar products plus 2 linear combinations
         instead of 3 scalar products and 3 linear combinations).  The
         rounding errors are however different, so one or the other
         formulation can be by chance more efficient.  Though there is no
         systematic trend. */
      double ytg = opk_vdot(opt->y, g);
      double dtg = opt->dtg;
      double ynorm = opk_vnorm2(opt->y);
      beta = (ytg - 2.0*(ynorm/dty)*ynorm*dtg)/dty;
    }
  } else {
    beta = 0.0;
  }
  return update1(opt, g, beta);
}

/* Perry & Shanno, update rule (used in CONMIN and see Eq. (1.4) in [3])
 * writes:
 *
 *     d' = alpha*(-c1*g + c2*d - c3*y)  ==>   d' = c1*g + c2*d + c3*y
 *
 *     c1 = (1/alpha)*<s,y>/<y,y>
 *        =  <d,y>/<y,y>
 *        = -<d,y>/<y,y>
 *
 *     c2 = <g,y>/<y,y> - 2*<s,g>/<s,y>
 *        = <g,y>/<y,y> - 2*<d,g>/<d,y>
 *        = <g,y>/<y,y> - 2*<d,g>/<d,y>
 *
 *     c3 = -(1/alpha)*<s,g>/<y,y>
 *        = -<d,g>/<y,y>
 *        =  <d,g>/<y,y>
 *
 * with alpha the step length, s = x - x0 = alpha*d = -alpha*d.  For this
 * method, beta = c2/c1.
 */
static int
update_Perry_Shanno(opk_nlcg_t* opt,
                    const opk_vector_t* x,
                    const opk_vector_t* g)
{
  double yty, dty, gty, dtg, c1, c2, c3;
  form_y(opt, g);
  yty = opk_vdot(opt->y, opt->y);
  if (yty <= 0) {
    return UPDATE_FAILURE;
  }
  dty = -opk_vdot(opt->d, opt->y);
  if (dty == 0) {
    return UPDATE_FAILURE;
  }
  gty = opk_vdot(g, opt->y);
  dtg = opt->dtg;
  c1 = dty/yty;
  c2 = gty/yty - 2.0*dtg/dty;
  c3 = -dtg/yty;
  opt->beta = c2/c1;
  opk_vaxpbypcz(opt->d, c1, g, c2, opt->d, c3, opt->y);
  return UPDATE_SUCCESS;
}

static void
finalize_nlcg(opk_object_t* obj)
{
  opk_nlcg_t* opt = (opk_nlcg_t*)obj;
  OPK_DROP(opt->vspace);
  OPK_DROP(opt->lnsrch);
  OPK_DROP(opt->x0);
  OPK_DROP(opt->g0);
  OPK_DROP(opt->d);
  OPK_DROP(opt->y);
}

/*---------------------------------------------------------------------------*/
/* PUBLIC INTERFACE */

opk_nlcg_t*
opk_new_nlcg_optimizer(opk_vspace_t* vspace,
                       unsigned int flags,
                       opk_lnsrch_t* lnsrch)
{
  opk_nlcg_t* opt;
  opk_bool_t g0_needed, y_needed;
  int (*update)(opk_nlcg_t* opt,
                const opk_vector_t* x,
                const opk_vector_t* g);

  /* Check the input arguments for errors. */
  if (vspace == NULL) {
    errno = EFAULT;
    return NULL;
  }
  if (flags == 0) {
    flags = OPK_NLCG_DEFAULT;
  }
  switch ((flags & 0xff)) {
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
  if (lnsrch != NULL) {
    opt->lnsrch = OPK_HOLD_LNSRCH(lnsrch);
  } else {
    opt->lnsrch = opk_lnsrch_new_csrch(SFTOL, SGTOL, SXTOL);
    if (opt->lnsrch == NULL) {
      goto error;
    }
  }
  opt->fmin_given = FALSE;
  opt->flags = flags;
  opt->update_Hager_Zhang_orig = FALSE;
  opk_set_nlcg_options(opt, NULL);

  /* Allocate work vectors. */
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
  opt->d = opk_vcreate(vspace);
  if (opt->d == NULL) {
    goto error;
  }
  if (y_needed) {
    opt->y = opk_vcreate(vspace);
    if (opt->y == NULL) {
      goto error;
    }
  }

  /* Enforce calling opk_nlcg_start and return the optimizer. */
  failure(opt, OPK_NOT_STARTED);
  return opt;

 error:
  OPK_DROP(opt);
  return NULL;
}

double
opk_get_nlcg_step(opk_nlcg_t* opt)
{
  return (opt == NULL ? -1.0 : opt->alpha);
}

double
opk_get_nlcg_gnorm(opk_nlcg_t* opt)
{
  return (opt == NULL ? -1.0 : opt->gnorm);
}

opk_status_t
opk_get_nlcg_fmin(opk_nlcg_t* opt, double* fmin)
{
  if (opt == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  if (opt->fmin_given) {
    if (fmin != NULL) {
      *fmin = opt->fmin;
    }
    return OPK_SUCCESS;
  } else {
    return OPK_UNDEFINED_VALUE;
  }
}

opk_status_t
opk_set_nlcg_fmin(opk_nlcg_t* opt, double fmin)
{
  if (opt == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  if (isnan(fmin) || isinf(fmin)) {
    return OPK_INVALID_ARGUMENT;
  }
  opt->fmin = fmin;
  opt->fmin_given = TRUE;
  return OPK_SUCCESS;
}

opk_status_t
opk_unset_nlcg_fmin(opk_nlcg_t* opt)
{
  if (opt == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  opt->fmin_given = FALSE;
  return OPK_SUCCESS;
}

opk_index_t
opk_get_nlcg_iterations(opk_nlcg_t* opt)
{
  return opt->iterations;
}

opk_index_t
opk_get_nlcg_restarts(opk_nlcg_t* opt)
{
  return opt->restarts;
}

opk_index_t
opk_get_nlcg_evaluations(opk_nlcg_t* opt)
{
  return opt->evaluations;
}

unsigned int
opk_get_nlcg_flags(opk_nlcg_t* opt)
{
  return (opt != NULL ? opt->flags : -1);
}

size_t
opk_get_nlcg_description(char* buf, size_t size, const opk_nlcg_t* opt)
{
  char str[80];
  int and;

  if (opt == NULL || (buf == NULL && size > 0)) {
    return 0;
  }
  switch (opt->flags & 0xff) {
  case OPK_NLCG_FLETCHER_REEVES:
    strcpy(str, "Fletcher & Reeves");
    break;
  case OPK_NLCG_HESTENES_STIEFEL:
    strcpy(str, "Hestenes & Stiefel");
    break;
  case OPK_NLCG_POLAK_RIBIERE_POLYAK:
    strcpy(str, "Polak, Ribière & Polyak");
    break;
  case OPK_NLCG_FLETCHER:
    strcpy(str, "Fletcher conjugate descent");
    break;
  case OPK_NLCG_LIU_STOREY:
    strcpy(str, "Liu & Storey");
    break;
  case OPK_NLCG_DAI_YUAN:
    strcpy(str, "Dai & Yuan");
    break;
  case OPK_NLCG_PERRY_SHANNO:
    strcpy(str, "Perry & Shanno");
    break;
  case OPK_NLCG_HAGER_ZHANG:
    strcpy(str, "Hager & Zhang");
    break;
  default:
    str[0] = '\0';
    return OPK_CORRUPTED_WORKSPACE;
  }
  strcat(str, " updates");
  if ((opt->flags & OPK_NLCG_POWELL) == OPK_NLCG_POWELL) {
    and = TRUE;
    strcat(str, " with Powell restarts");
  } else {
    and = FALSE;
  }
  if ((opt->flags & OPK_NLCG_SHANNO_PHUA) == OPK_NLCG_SHANNO_PHUA) {
    strcat(str, (and ? " and" : " with"));
    strcat(str, " Shanno & Phua step size");
  }
  return opk_copy_string(buf, size, str);
}

opk_task_t
opk_get_nlcg_task(opk_nlcg_t* opt)
{
  return opt->task;
}

opk_status_t
opk_get_nlcg_status(opk_nlcg_t* opt)
{
  return opt->status;
}

double
opk_get_nlcg_beta(opk_nlcg_t* opt)
{
  return opt->beta;
}

opk_status_t
opk_get_nlcg_options(opk_nlcg_options_t* dst, const opk_nlcg_t* src)
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
opk_set_nlcg_options(opk_nlcg_t* dst, const opk_nlcg_options_t* src)
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

opk_task_t
opk_start_nlcg(opk_nlcg_t* opt, opk_vector_t* x)
{
  opt->iterations = 0;
  opt->evaluations = 0;
  opt->restarts = 0;
  return success(opt, OPK_TASK_COMPUTE_FG);
}

opk_task_t
opk_iterate_nlcg(opk_nlcg_t* opt, opk_vector_t* x,
                 double f, opk_vector_t* g)
{
  /*
   * The new iterate is:
   *    x_{k+1} = x_{k} - \alpha_{k} d_{k}
   * as we consider the anti-search direction here.
   */
  opk_status_t status = 0;
  opk_task_t next_task;
  opk_lnsrch_task_t lnsrch_task = 0;

  switch (opt->task) {

  case OPK_TASK_COMPUTE_FG:

    ++opt->evaluations;
    if (opt->evaluations > 1) {
      /* Line search in progress. Compute directional derivative and check
         whether line search has converged. */
      opt->dtg = -opk_vdot(opt->d, g);
      lnsrch_task = opk_lnsrch_iterate(opt->lnsrch, &opt->alpha, f, opt->dtg);
      if (lnsrch_task != OPK_LNSRCH_CONVERGENCE) {
        if (lnsrch_task == OPK_LNSRCH_SEARCH) {
          /* Line search has not converged, break to compute a new trial point
             along the search direction. */
          break;
        }
        status = opk_lnsrch_get_status(opt->lnsrch);
        if (lnsrch_task != OPK_LNSRCH_WARNING ||
            status != OPK_ROUNDING_ERRORS_PREVENT_PROGRESS) {
          return failure(opt, status);
        }
      }
      /* Line search has converged. */
      ++opt->iterations;
    }

    /* The current step is acceptable.  Check for global convergence. */
    opt->gnorm = opk_vnorm2(g);
    if (opt->evaluations <= 1) {
      opt->ginit = opt->gnorm;
    }
    if (opt->gnorm <= max3(0.0, opt->gatol, opt->grtol*opt->ginit)) {
      next_task = OPK_TASK_FINAL_X;
    } else {
      next_task = OPK_TASK_NEW_X;
    }
    return success(opt, next_task);

  case OPK_TASK_NEW_X:
  case OPK_TASK_FINAL_X:

    /* Compute search direction and initial step size. */
    if (opt->evaluations <= 1 || opt->update(opt, x, g) != UPDATE_SUCCESS) {
      /* First evaluation or update failed, set DTG to zero to use the steepest
         descent direction. */
      opt->dtg = 0;
    } else {
      opt->dtg = -opk_vdot(opt->d, g);
      if (opt->epsilon > 0 &&
          opt->dtg > -opt->epsilon*opk_vnorm2(opt->d)*opt->gnorm) {
        /* Set DTG to zero to indicate that we do not have a sufficient
           descent direction. */
        opt->dtg = 0;
      }
    }
    if (opt->dtg < 0) {
      /* The recursion yields a sufficient descent direction (not all methods
         warrant that).  Compute an initial step size ALPHA along the new
         direction. */
      if ((opt->flags & OPK_NLCG_SHANNO_PHUA) == OPK_NLCG_SHANNO_PHUA) {
        /* Initial step size is such that:
           <alpha_{k+1}*d_{k+1},g_{k+1}> = <alpha_{k}*d_{k},g_{k}> */
        opt->alpha *= (opt->dtg0/opt->dtg);
      }
    } else {
      /* Initial search direction or recurrence has been restarted.  FIXME:
         other possibility is to use Fletcher's formula, see BGLS p. 39) */
      if (opt->evaluations > 1) {
        ++opt->restarts;
      }
      opk_vcopy(opt->d, g);
      opt->dtg = -opt->gnorm*opt->gnorm;
      if (opt->fmin_given && opt->fmin < f) {
        opt->alpha = 2*(opt->fmin - f)/opt->dtg;
      } else if (f != 0) {
        opt->alpha = 2*fabs(f/opt->dtg);
      } else {
        double dnorm = opt->gnorm;
        double xnorm = opk_vnorm2(x);
        if (xnorm > 0) {
          opt->alpha = opt->delta*xnorm/dnorm;
        } else {
          opt->alpha = opt->delta/dnorm;
        }
      }
      opt->beta = 0;
    }

    /* Store current position as X0, f0, etc. */
    opk_vcopy(opt->x0, x);
    opt->f0 = f;
    if (opt->g0 != NULL) {
      opk_vcopy(opt->g0, g);
    }
    opt->g0norm = opt->gnorm;
    opt->dtg0 = opt->dtg;

    /* Start the line search and break to compute the first trial point along
       the line search. */
    if (opk_lnsrch_start(opt->lnsrch, opt->f0, opt->dtg0, opt->alpha,
                         opt->stpmin*opt->alpha,
                         opt->stpmax*opt->alpha) != OPK_LNSRCH_SEARCH) {
      return failure(opt, opk_lnsrch_get_status(opt->lnsrch));
    }
    break;

  default:

    /* There is probably something wrong. */
    return opt->task;
  }

  /* Compute a trial point along the line search. */
  opk_vaxpby(x, 1.0, opt->x0, -opt->alpha, opt->d);
  return success(opt, OPK_TASK_COMPUTE_FG);
}
