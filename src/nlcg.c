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
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>

#include "optimpack-private.h"

/*---------------------------------------------------------------------------*/
/* PRIVATE DEFINITIONS */

static const double STPMIN = 1E-20;
static const double STPMAX = 1E+20;

struct _opk_nlcg_workspace {
  double f0;     /* Function value at the start of the line search. */
  double g0norm; /* Euclidean norm of G0, the gradient at the start of the
                    line search. */
  double g1norm; /* Euclidean norm of G1, the gradient of the end of the line
                    search / last accepted point. */
  double dg0;    /* Directional derivative at the start of the line search;
                    given by the inner product: <d,g0> = - <p,g0> */
  double dg1;    /* Directional derivative at the end or during the line
                    search; given by the inner product: <d,g1> = - <p,g1> */
#if 0 /* not used */
  double alpha0; /* Scale factor for the initial step step size. */
#endif
  double grtol;  /* Relative threshold for the norm or the gradient (relative
                    to the initial gradient) for convergence. */
  double gatol;  /* Absolute threshold for the norm or the gradient for
                    convergence. */
  double gtest;  /* Threshold for the norm or the gradient for
                    convergence: GTEST = MAX(0, GATOL, GRTOL*NORM(GINIT)) */
  double fmin;   /* Minimal function value if provided. */
  double alpha;  /* Current step length. */
  double beta;   /* Current parameter in conjugate gradient update rule (for
		    information). */
  double stpmin; /* Relative lower bound for the step length. */
  double stpmax; /* Relative upper bound for the step length. */
  opk_vspace_t* vspace;
  int (*update)(opk_nlcg_workspace_t* ws,
                const opk_vector_t* x1,
                const opk_vector_t* g1); /* This "method" is called to update
                                            the search direction.  The
                                            returned value indicates whether
                                            the updating rule has been
                                            successful, otherwise a restart is
                                            needed. */
  opk_lnsrch_workspace_t* lnsrch;
  opk_vector_t* x0;      /* Variables at start of line search. */
  opk_vector_t* g0;	 /* Gradient at start of line search. */
  opk_vector_t* p;	 /* (Anti-)search direction, new iterate is searched
                            as: x1 = x0 - alpha*p, for alpha >= 0. */
  opk_vector_t* y;       /* Work vector (e.g., to store the gradient
                            difference: Y = G1 - G0). */
  int iter;              /* Iteration number. */
  int nrestarts;         /* Number of algorithm restarts. */
  int nevals;            /* Number of function and gradient evaluations. */
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
update0(opk_nlcg_workspace_t* ws,
        const opk_vector_t* g1,
        double beta)
{
  if ((ws->beta = beta) != 0.0) {
    opk_vaxpby(ws->p, 1.0, g1, beta, ws->p);
    return OPK_SUCCESS;
  } else {
    return OPK_FAILURE;
  }
}

/* Helper function to compute search direction as: p' = g1 + beta*p
   possibly with the constraint that beta > 0. */
static int
update1(opk_nlcg_workspace_t* ws,
        const opk_vector_t* g1,
        double beta)
{
  if ((ws->method & OPK_NLCG_POWELL) != 0 && beta < 0.0) {
    ++ws->nrestarts;
    beta = 0.0;
  }
  if ((ws->beta = beta) != 0.0) {
    opk_vaxpby(ws->p, 1.0, g1, beta, ws->p);
    return OPK_SUCCESS;
  } else {
    return OPK_FAILURE;
  }
}

/* Form: Y = G1 - G0 */
static void
form_y(opk_nlcg_workspace_t* ws,
       const opk_vector_t* g1)
{
  opk_vaxpby(ws->y, 1.0, g1, -1.0, ws->g0);
}

/*
 * For Hestenes & Stiefel method:
 *
 *     beta = <g1,y>/<d,y> = - <g1,y>/<p,y>
 *
 * with y = g1 - g0.
 */
static int
update_Hestenes_Stiefel(opk_nlcg_workspace_t* ws,
			const opk_vector_t* x1,
			const opk_vector_t* g1)
{
  double g1y, dy, beta;
  form_y(ws, g1);
  g1y = opk_vdot(g1, ws->y);    /* Compute: g1y = <g1,y> */
  dy = -opk_vdot(ws->p, ws->y); /* Compute: dy = <d,y> = - <p,y> */
  beta = (dy != 0.0 ? g1y/dy : 0.0);
  return update1(ws, g1, beta);
}

/*
 * For Fletcher & Reeves method:
 *
 *     beta = <g1,g1>/<g0,g0>
 *
 * (this value is always >= 0 and can only be zero at a stationary point).
 */
static int
update_Fletcher_Reeves(opk_nlcg_workspace_t* ws,
		       const opk_vector_t* x1,
		       const opk_vector_t* g1)
{
  double r = ws->g1norm/ws->g0norm;
  return update0(ws, g1, r*r);
}

/*
 * For Polak-Ribière-Polyak method:
 *
 *     beta = <g1,y>/<g0,g0>
 */
static int
update_Polak_Ribiere_Polyak(opk_nlcg_workspace_t* ws,
			    const opk_vector_t* x1,
			    const opk_vector_t* g1)
{
  double beta;
  form_y(ws, g1);
  beta = opk_vdot(g1, ws->y)/(ws->g0norm*ws->g0norm);
  return update1(ws, g1, beta);
}

/*
 * For Fletcher "Conjugate Descent" method:
 *
 *     beta = <g1,g1>/(-<d,g0>)
 *
 * (this value is always >= 0 and can only be zero at a stationnary point).
 */
static int
update_Fletcher(opk_nlcg_workspace_t* ws,
		const opk_vector_t* x1,
		const opk_vector_t* g1)
{
  double beta = ws->g1norm*(ws->g1norm/(-ws->dg0));
  return update0(ws, g1, beta);
}

/*
 * For Liu & Storey method:
 *
 *     beta = <g1,y>/(-<d,g0>)
 */
static int
update_Liu_Storey(opk_nlcg_workspace_t* ws,
		  const opk_vector_t* x1,
		  const opk_vector_t* g1)
{
  double g1y, beta;
  form_y(ws, g1);
  g1y =  opk_vdot(g1, ws->y);    /* Compute: g1y = <g1,y> */
  beta = g1y/(-ws->dg0);
  return update1(ws, g1, beta);
}

/*
 * For Dai & Yuan method:
 *
 *     beta = <g1,g1>/<d,y>
 */
static int
update_Dai_Yuan(opk_nlcg_workspace_t* ws,
		const opk_vector_t* x1,
		const opk_vector_t* g1)
{
  double dy, beta;
  form_y(ws, g1);
  dy = -opk_vdot(ws->p, ws->y); /* Compute: dy = <d,y> = - <p,y> */
  beta = (dy != 0.0 ? ws->g1norm*(ws->g1norm/dy) : 0.0);
  return update1(ws, g1, beta);
}

/*
 * For Hager & Zhang method:
 *
 *     beta = <y - (2*<y,y>/<d,y>)*d,g1>
 *          = (<g1,y> - 2*<y,y>*<d,g1>/<d,y>)/<d,y>
 */
static int
update_Hager_Zhang(opk_nlcg_workspace_t* ws,
		   const opk_vector_t* x1,
		   const opk_vector_t* g1)
{
  double dy, beta;
  form_y(ws, g1);
  dy = -opk_vdot(ws->p, ws->y);
  if (dy != 0.0) {
    if (ws->update_Hager_Zhang_orig) {
      /* Original formulation. */
      double q = 1.0/dy;
      double r = q*opk_vnorm2(ws->y);
      opk_vaxpby(ws->y, q, ws->y, 2.0*r*r, ws->p);
      beta = opk_vdot(ws->y, g1);
    } else {
      /* Improved formulation which spares one linear combination and thus has
         less overhead (only 3 scalar products plus 2 linear combinations
         instead of 3 scalar products and 3 linear combinations).  The
         rounding errors are however different, so one or the other
         formulation can be by chance more efficient.  Though there is no
         systematic trend. */
      double yg = opk_vdot(ws->y, g1);
      double dg = ws->dg1;
      double r = opk_vnorm2(ws->y)/dy;
      beta = yg/dy - 2.0*r*r*dg;
    }
  } else {
    beta = 0.0;
  }
  return update1(ws, g1, beta);
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
update_Perry_Shanno(opk_nlcg_workspace_t* ws,
                    const opk_vector_t* x1,
                    const opk_vector_t* g1)
{
  double yy, dy, g1y, dg1, c1, c2, c3;
  form_y(ws, g1);
  yy = opk_vdot(ws->y, ws->y);
  if (yy <= 0.0) return OPK_FAILURE;
  dy = -opk_vdot(ws->p, ws->y);
  if (dy == 0.0) return OPK_FAILURE;
  g1y = opk_vdot(g1, ws->y);
  dg1 = ws->dg1;
  c1 = dy/yy;
  c2 = g1y/yy - 2.0*dg1/dy;
  c3 = -dg1/yy;
  ws->beta = c2/c1;
  opk_vaxpbypcz(ws->p, c1, g1, c2, ws->p, c3, ws->y);
  return OPK_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/* PUBLIC ROUTINES */

opk_nlcg_workspace_t*
opk_nlcg_new_with_line_search(opk_vspace_t* vspace, unsigned int method,
			      opk_lnsrch_workspace_t* lnsrch)
{
  opk_nlcg_workspace_t* ws;
  int (*update)(opk_nlcg_workspace_t* ws,
                const opk_vector_t* x1,
                const opk_vector_t* g1);
  opk_bool_t g0_needed, y_needed;

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
  ws = NEW(opk_nlcg_workspace_t);
  if (ws == NULL) {
    return NULL;
  }
  ws->vspace = vspace; /* must be set early in case of errors */
  ws->update = update;
  ws->fmin_given = FALSE;
  ws->method = method;
  ws->lnsrch = lnsrch;
  ws->grtol = 1e-6;
  ws->gatol = 0.0;
  ws->stpmin = STPMIN;
  ws->stpmax = STPMAX;
  ws->x0 = opk_vcreate(vspace);
  if (ws->x0 == NULL) {
    goto error;
  }
  if (g0_needed) {
    ws->g0 = opk_vcreate(vspace);
    if (ws->g0 == NULL) {
      goto error;
    }
  }
  ws->p = opk_vcreate(vspace);
  if (ws->p == NULL) {
    goto error;
  }
  if (y_needed) {
    ws->y = opk_vcreate(vspace);
    if (ws->y == NULL) {
      goto error;
    }
  }
  ws->update_Hager_Zhang_orig = OPK_FALSE;

  opk_nlcg_start(ws);
  return ws;

error:
  opk_nlcg_delete(ws);
  return NULL;
}


opk_nlcg_workspace_t*
opk_nlcg_new(opk_vspace_t* vspace, unsigned int method)
{
  opk_lnsrch_workspace_t* lnsrch;
  /* FIXME: choose more suitable values (e.g., in CG+: FTOL=1E-4, GTOL=1E-1,
     not less than 1E-4, XTOL=1E-17, STPMIN=1E-20 STPMAX=1E+20 and MAXFEV=40) */
  lnsrch = opk_lnsrch_new_csrch(/* sftol  */  1e-4,
				/* sgtol  */  1e-1,
				/* sxtol  */  1E-17);
  if (lnsrch == NULL) {
    return NULL;
  }
  return opk_nlcg_new_with_line_search(vspace, method, lnsrch);
}

void
opk_nlcg_delete(opk_nlcg_workspace_t* ws)
{
  if (ws != NULL) {
    /* Delete the vectors of the workspace. */
    opk_vspace_t* vspace = ws->vspace;
    if (vspace != NULL && vspace->delete != NULL) {
      opk_vdelete(ws->x0);
      opk_vdelete(ws->g0);
      opk_vdelete(ws->p);
      opk_vdelete(ws->y);
    }

    /* Delete the line search workspace. */
    if (ws->lnsrch != NULL) {
      opk_lnsrch_delete(ws->lnsrch);
    }

    /* Delete the workspace. */
    free(ws);
  }
}

int
opk_nlcg_get_fmin(opk_nlcg_workspace_t* ws, double* fmin)
{
  if (ws == NULL) {
    errno = EFAULT;
    return OPK_FAILURE;
  }
  if (ws->fmin_given) {
    if (fmin != NULL) {
      *fmin = ws->fmin;
    }
    return OPK_SUCCESS;
  } else {
    return OPK_FAILURE;
  }
}

int
opk_nlcg_set_fmin(opk_nlcg_workspace_t* ws, double fmin)
{
  if (ws == NULL) {
    errno = EFAULT;
    return OPK_FAILURE;
  }
  if (isnan(fmin) || isinf(fmin)) {
    errno = EINVAL;
    return OPK_FAILURE;
  }
  ws->fmin = fmin;
  ws->fmin_given = TRUE;
  return OPK_SUCCESS;
}

int
opk_nlcg_unset_fmin(opk_nlcg_workspace_t* ws)
{
  if (ws == NULL) {
    errno = EFAULT;
    return OPK_FAILURE;
  }
  ws->fmin_given = FALSE;
  return OPK_SUCCESS;
}

int
opk_nlcg_get_iterations(opk_nlcg_workspace_t* ws)
{
  return ws->iter;
}

int
opk_nlcg_get_restarts(opk_nlcg_workspace_t* ws)
{
  return ws->nrestarts;
}

int
opk_nlcg_get_evaluations(opk_nlcg_workspace_t* ws)
{
  return ws->nevals;
}

unsigned int
opk_nlcg_get_method(opk_nlcg_workspace_t* ws)
{
  return ws->method;
}

opk_bool_t
opk_nlcg_get_starting(opk_nlcg_workspace_t* ws)
{
  return ws->start;
}

opk_task_t
opk_nlcg_get_task(opk_nlcg_workspace_t* ws)
{
  return ws->task;
}

double
opk_nlcg_get_alpha(opk_nlcg_workspace_t* ws)
{
  return ws->alpha;
}

double
opk_nlcg_get_beta(opk_nlcg_workspace_t* ws)
{
  return ws->beta;
}

opk_task_t
opk_nlcg_start(opk_nlcg_workspace_t* ws)
{
  ws->iter = 0;
  ws->nevals = 0;
  ws->nrestarts = 0;
  ws->task = OPK_TASK_COMPUTE_FG;
  ws->start = TRUE;
  return ws->task;
}

opk_task_t
opk_nlcg_iterate(opk_nlcg_workspace_t* ws, opk_vector_t* x1,
		 double f1, opk_vector_t* g1)
{
  /*
   * The new iterate is:
   *    x_{k+1} = x_{k} + \alpha_{k} d_{k}
   *            = x_{k} - \alpha_{k} p_{k}
   * as we consider the anti-search direction p = -d here.
   */
  int status = 0;

  if (ws->task == OPK_TASK_COMPUTE_FG) {
    opk_bool_t accept;
    ++ws->nevals;
    if (ws->start) {
      ws->g1norm = opk_vnorm2(g1);
      ws->gtest = max3(0.0, ws->gatol, ws->grtol*ws->g1norm);
      accept = TRUE;
    } else {
      /* Compute directional derivative and check whether line search has
         converged. */
      ws->dg1 = -opk_vdot(ws->p, g1);
      status = opk_lnsrch_iterate(ws->lnsrch, &ws->alpha, f1, ws->dg1);
      if (status != 0) {
	/* Line search is finished, hopefully because of convergence
	   (FIXME: check for warnings). */
	if (status < 0) {
	  goto line_search_error;
	}
        ++ws->iter;
	ws->g1norm = opk_vnorm2(g1);
        accept = TRUE;
      } else {
        accept = FALSE;
      }
    }
    if (accept) {
      /* Check for global convergence. */
      if (ws->g1norm <= ws->gtest) {
	ws->task = OPK_TASK_FINAL_X;
      } else {
	ws->task = OPK_TASK_NEW_X;
      }
      return ws->task;
    }
  } else if (ws->task == OPK_TASK_NEW_X) {
    /* Compute a search direction and start line search. */
    int restart;
    if (ws->start) {
      restart = TRUE;
      ws->beta = 0.0;
    } else {
      restart = (ws->update(ws, x1, g1) != OPK_SUCCESS);
      if (! restart) {
        double dg = -opk_vdot(ws->p, g1);
        if (dg >= 0.0) {
          /* Restart if not a descent direction, not all updates warrant that.
             (FIXME: Generate an error instead?) */
          restart = TRUE;
        } else {
          /* Compute an initial step size ALPHA. */
          if ((ws->method & OPK_NLCG_SHANNO_PHUA) != 0) {
            /* Initial step size is such that:
               <alpha_{k+1}*d_{k+1},g_{k+1}> = <alpha_{k}*d_{k},g_{k}> */
            ws->alpha *= ws->dg0/dg;
#if 0 /* keep same ALPHA */
          } else {
            /* FIXME: check this */
            ws->alpha = (-ws->dg0)/p0norm/p0norm/ws->alpha;
#endif
          }
        }
        /* Save directional derivative. */
        ws->dg0 = dg;
      }
    }
    if (restart) {
      /* Initial search direction or recurrence has been restarted.  FIXME:
	 other possibility is to use Fletcher's formula, see BGLS p. 39) */
      if (! ws->start) {
        ++ws->nrestarts;
      }
      ws->beta = 0.0;
#if 0
      double x1norm = opk_vnorm2(x1);
      if (x1norm > 0.0) {
	ws->alpha = (x1norm/ws->g1norm)*ws->tiny;
      } else {
	ws->alpha = (1e-3*MAX(fabs(f1), 1.0)/ws->g1norm)/ws->g1norm;
      }
#else
      ws->alpha = 1.0/ws->g1norm;
#endif
      opk_vcopy(ws->p, g1);
      ws->dg0 = -ws->g1norm*ws->g1norm;
    }

    /* Store current position as X0, f0, etc. */
    opk_vcopy(ws->x0, x1);
    ws->f0 = f1;
    if (ws->g0 != NULL) {
      opk_vcopy(ws->g0, g1);
    }
    ws->g0norm = ws->g1norm;

    /* Start the line search. */
    ws->dg1 = ws->dg0;
    status = opk_lnsrch_start(ws->lnsrch, ws->f0, ws->dg0, ws->alpha,
                              ws->stpmin*ws->alpha,
                              ws->stpmax*ws->alpha);
    if (status < 0) {
      goto line_search_error;
    }
  } else {
    return ws->task;
  }

  /* Build a new step to try. */
  opk_vaxpby(x1, 1.0, ws->x0, -ws->alpha, ws->p);
  ws->start = FALSE;
  ws->task = OPK_TASK_COMPUTE_FG;
  return ws->task;

  /* Some errors occurred in line search. */
 line_search_error:
  fprintf(stderr, "opk_lnsrch_start error: %d\n", status);
  /* FIXME: add line search error state */
  ws->task = OPK_TASK_ERROR;
  return ws->task;
}

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * c-basic-offset: 2
 * fill-column: 78
 * coding: utf-8
 * ispell-local-dictionary: "american"
 * End:
 */
