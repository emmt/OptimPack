/*
 * lnsrch.c --
 *
 * Line search routines for OptimPack library.  Implements Armijo and Moré &
 * Thuente inexact line search methods.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2003-2014 Éric Thiébaut
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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "optimpack-private.h"

#if (OPK_LNSRCH_SEARCH != 0)
# error OPK_LNSRCH_SEARCH != 0
#endif

/*---------------------------------------------------------------------------*/
/* UNIFIED INTERFACE FOR LINE SEARCH */

opk_lnsrch_workspace_t*
_opk_lnsrch_new(size_t size,
		int (*start)(opk_lnsrch_workspace_t* ws),
		int (*iterate)(opk_lnsrch_workspace_t* ws,
			       double* stp_ptr, double f1, double g1),
	        void (*delete)(opk_lnsrch_workspace_t* ws))
{
  opk_lnsrch_workspace_t* ws;

  if (start == NULL || iterate == NULL) {
    errno = EFAULT;
    return NULL;
  }
  size = MAX(size, sizeof(opk_lnsrch_workspace_t));
  ws = (opk_lnsrch_workspace_t*)malloc(size);
  if (ws != NULL) {
    memset(ws, 0, size);
    ws->start = start;
    ws->iterate = iterate;
    ws->delete = delete;
    ws->status = OPK_LNSRCH_ERROR_NOT_STARTED;
  }
  return ws;
}

/* after an error or convergence, you must call opk_lnsrch_start */

int
opk_lnsrch_start(opk_lnsrch_workspace_t* ws, double f0, double g0,
		 double stp, double stpmin, double stpmax)
{
  if (ws == NULL) {
    return OPK_LNSRCH_ERROR_ILLEGAL_ADDRESS;
  }
  if (stpmin < 0.0) {
    ws->status = OPK_LNSRCH_ERROR_STPMIN_LT_ZERO;
  } else if (stpmin > stpmax) {
    ws->status = OPK_LNSRCH_ERROR_STPMIN_GT_STPMAX;
  } else if (stp < stpmin) {
    ws->status = OPK_LNSRCH_ERROR_STP_LT_STPMIN;
  } else if (stp > stpmax) {
    ws->status = OPK_LNSRCH_ERROR_STP_GT_STPMAX;
  } else if (g0 >= 0.0) {
    ws->status = OPK_LNSRCH_ERROR_INITIAL_DERIVATIVE_GE_ZERO;
  } else {
    ws->stp = stp;
    ws->stpmin = stpmin;
    ws->stpmax = stpmax;
    ws->finit = f0;
    ws->ginit = g0;
    ws->status = ws->start(ws);
  }
  return ws->status;
}

int
opk_lnsrch_iterate(opk_lnsrch_workspace_t* ws, double* stp_ptr,
		   double f1, double g1)
{
  if (ws == NULL || stp_ptr == NULL) {
    return OPK_LNSRCH_ERROR_ILLEGAL_ADDRESS;
  }
  if (ws->status == OPK_LNSRCH_SEARCH) {
    if (*stp_ptr != ws->stp) {
      ws->status = OPK_LNSRCH_ERROR_STP_CHANGED;
    } else {
      ws->status = ws->iterate(ws, stp_ptr, f1, g1);
      if (*stp_ptr > ws->stpmax) {
	if (ws->stp >= ws->stpmax) {
	  ws->status = OPK_LNSRCH_WARNING_STP_EQ_STPMAX;
	}
	*stp_ptr = ws->stpmax;
      } else if (*stp_ptr < ws->stpmin) {
	if (ws->stp <= ws->stpmin) {
	  ws->status = OPK_LNSRCH_WARNING_STP_EQ_STPMIN;
	}
	*stp_ptr = ws->stpmin;
      }
      ws->stp = *stp_ptr;
    }
  } else {
    ws->status = OPK_LNSRCH_ERROR_NOT_STARTED;
  }
  return ws->status;
}

void
opk_lnsrch_delete(opk_lnsrch_workspace_t* ws)
{
  if (ws != NULL) {
    if (ws->delete != NULL) {
      ws->delete(ws);
    }
    free(ws);
  }
}

double
opk_lnsrch_get_step(const opk_lnsrch_workspace_t* ws)
{
  return (ws != NULL && ws->status == OPK_LNSRCH_SEARCH ? ws->stp : -1.0);
}

int
opk_lnsrch_get_status(const opk_lnsrch_workspace_t* ws)
{
  return (ws != NULL ? ws->status : OPK_LNSRCH_ERROR_ILLEGAL_ADDRESS);
}

int
opk_lnsrch_has_errors(const opk_lnsrch_workspace_t* ws)
{
  return (opk_lnsrch_get_status(ws) < 0);
}

int
opk_lnsrch_has_warnings(const opk_lnsrch_workspace_t* ws)
{
  return (opk_lnsrch_get_status(ws) > OPK_LNSRCH_CONVERGENCE);
}

int
opk_lnsrch_converged(const opk_lnsrch_workspace_t* ws)
{
  return (opk_lnsrch_get_status(ws) == OPK_LNSRCH_CONVERGENCE);
}

int
opk_lnsrch_finished(const opk_lnsrch_workspace_t* ws)
{
  return (opk_lnsrch_get_status(ws) != OPK_LNSRCH_SEARCH);
}

/*---------------------------------------------------------------------------*/
/* ARMIJO (BACKTRACKING) LINE SEARCH */

/* Workspace for Armijo line search. */
typedef struct _backtrack_workspace backtrack_workspace_t;
struct _backtrack_workspace {
  opk_lnsrch_workspace_t super;
  double ftol;
};

/* Initialize workspace. */
static int
backtrack_start(opk_lnsrch_workspace_t* _ws)
{
  return OPK_LNSRCH_SEARCH;
}

static int
backtrack_iterate(opk_lnsrch_workspace_t* _ws,
		  double* stp_ptr, double f1, double g1)
{
  backtrack_workspace_t* ws = (backtrack_workspace_t*)_ws;
  int status;

  if (f1 <= ws->super.finit + ws->ftol*(*stp_ptr)*ws->super.ginit) {
    /* First Wolfe conditions satisfied. */
    status = OPK_LNSRCH_CONVERGENCE;
  } else {
    /* Take a bisection step unless already at the lower bound. */
    if (*stp_ptr <= ws->super.stpmin) {
      status = OPK_LNSRCH_WARNING_STP_EQ_STPMIN;
    } else {
      *stp_ptr = (*stp_ptr + ws->super.stpmin)*0.5;
    }
    if (*stp_ptr < ws->super.stpmin) {
      *stp_ptr = ws->super.stpmin;
    }
  }
  return status;
}

opk_lnsrch_workspace_t*
opk_lnsrch_new_backtrack(double ftol)
{
  opk_lnsrch_workspace_t* _ws;

  if (ftol < 0.0) {
    errno = EINVAL;
    return NULL;
  }
  _ws = _opk_lnsrch_new(sizeof(backtrack_workspace_t),
			backtrack_start, backtrack_iterate, NULL);
  if (_ws != NULL) {
    backtrack_workspace_t* ws = (backtrack_workspace_t*)_ws;
    ws->ftol = ftol;
  }
  return _ws;
}

/*---------------------------------------------------------------------------*/
/* MORÉ AND THUENTE CUBIC LINE SEARCH */

/* Workspace for Moré and Thuente cubic line search. */
typedef struct _csrch_workspace csrch_workspace_t;
struct _csrch_workspace {
  opk_lnsrch_workspace_t super;

  /* Convergence parameters. */
  double ftol;
  double gtol;
  double xtol;

  /* GTEST is used to check for Wolfe conditions. */
  double gtest;

  /* The variables STX, FX, GX contain the values of the step, function, and
     derivative at the best step. */
  double stx, fx, gx;

  /* The variables STY, FY, GY contain the value of the step, function, and
     derivative at STY. */
  double sty, fy, gy;

  /* Parameters to track the interval where to seek for the step. */
  double stmin;
  double stmax;
  double width;
  double width1;
  opk_bool_t brackt;

  /* The algorithm has two different stages. */
  int stage;
};

static csrch_workspace_t*
csrch_get_workspace(opk_lnsrch_workspace_t* ws);

static int
csrch_start(opk_lnsrch_workspace_t* _ws)
{
  csrch_workspace_t* ws;

  ws = csrch_get_workspace(_ws);
  if (ws == NULL) {
    return OPK_LNSRCH_ERROR_CORRUPTED_WORKSPACE;
  }

  /* Convergence threshold for this step. */
  ws->gtest = ws->ftol*ws->super.ginit;

  /* Initialize parameters for the interval of search. */
  ws->stmin = ws->super.stpmin;
  ws->stmax = ws->super.stpmax;
  ws->width = ws->super.stpmax - ws->super.stpmin;
  ws->width1 = ws->width/0.5;
  ws->brackt = FALSE;

  /* The variables STX, FX, GX contain the values of the step,
     function, and derivative at the best step. */
  ws->stx = 0.0;
  ws->fx = ws->super.finit;
  ws->gx = ws->super.ginit;

  /* The variables STY, FY, GY contain the value of the step,
     function, and derivative at STY. */
  ws->sty = 0.0;
  ws->fy = ws->super.finit;
  ws->gy = ws->super.ginit;

  ws->stage = 1;
  return OPK_LNSRCH_SEARCH;
}

static int
csrch_iterate(opk_lnsrch_workspace_t* _ws,
	      double* stp_ptr, double f1, double g1)
{
  double ftest;
  csrch_workspace_t* ws;
  int result;

  ws = csrch_get_workspace(_ws);
  if (ws == NULL) {
    return OPK_LNSRCH_ERROR_CORRUPTED_WORKSPACE;
  }

  /* Test for convergence. */
  ftest = ws->super.finit + (*stp_ptr)*ws->gtest;
  if (f1 <= ftest && fabs(g1) <= -ws->gtol*ws->super.ginit) {
    /* Strong Wolfe conditions satisfied. */
    return OPK_LNSRCH_CONVERGENCE;
  }

  /* Test for warnings. */
  if (*stp_ptr == ws->super.stpmin && (f1 > ftest || g1 >= ws->gtest)) {
    return OPK_LNSRCH_WARNING_STP_EQ_STPMIN;
  }
  if (*stp_ptr == ws->super.stpmax && f1 <= ftest && g1 <= ws->gtest) {
    return OPK_LNSRCH_WARNING_STP_EQ_STPMAX;
  }
  if (ws->brackt && ws->stmax - ws->stmin <= ws->xtol * ws->stmax) {
    return OPK_LNSRCH_WARNING_XTOL_TEST_SATISFIED;
  }
  if (ws->brackt && (*stp_ptr <= ws->stmin || *stp_ptr >= ws->stmax)) {
    return OPK_LNSRCH_WARNING_ROUNDING_ERRORS_PREVENT_PROGRESS;
  }

  /* If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
     algorithm enters the second stage. */
  if (ws->stage == 1 && f1 <= ftest && g1 >= 0.0) {
    ws->stage = 2;
  }

  /* A modified function is used to predict the step during the first stage if
     a lower function value has been obtained but the decrease is not
     sufficient. */
  if (ws->stage == 1 && f1 <= ws->fx && f1 > ftest) {
    /* Define the modified function and derivative values and call CSTEP to
       update STX, STY, and to compute the new step.  Then restore the
       function and derivative values for F.*/
    double fm = f1 - *stp_ptr * ws->gtest;
    double fxm = ws->fx - ws->stx * ws->gtest;
    double fym = ws->fy - ws->sty * ws->gtest;
    double gm = g1 - ws->gtest;
    double gxm = ws->gx - ws->gtest;
    double gym = ws->gy - ws->gtest;
    result = opk_cstep(&ws->stx, &fxm, &gxm,
		       &ws->sty, &fym, &gym,
		       stp_ptr, fm, gm,
		       &ws->brackt, ws->stmin, ws->stmax);
    if (result < 0) {
      return result;
    }
    ws->fx = fxm + ws->stx * ws->gtest;
    ws->fy = fym + ws->sty * ws->gtest;
    ws->gx = gxm + ws->gtest;
    ws->gy = gym + ws->gtest;
  } else {
    /* Call CSTEP to update STX, STY, and to compute the new step. */
    result = opk_cstep(&ws->stx, &ws->fx, &ws->gx,
		       &ws->sty, &ws->fy, &ws->gy,
		       stp_ptr, f1, g1,
		       &ws->brackt, ws->stmin, ws->stmax);
    if (result < 0) {
      return result;
    }
  }

  /* Decide if a bisection step is needed. */
  if (ws->brackt) {
    double new_width = fabs(ws->sty - ws->stx);
    if (new_width >= 0.66*ws->width1) {
      *stp_ptr = ws->stx + 0.5*(ws->sty - ws->stx);
    }
    ws->width1 = ws->width;
    ws->width = new_width;
  }

  /* Set the minimum and maximum steps allowed for stp. */
  if (ws->brackt) {
    ws->stmin = MIN(ws->stx, ws->sty);
    ws->stmax = MAX(ws->stx, ws->sty);
  } else {
    ws->stmin = *stp_ptr + (*stp_ptr - ws->stx)*1.1;
    ws->stmax = *stp_ptr + (*stp_ptr - ws->stx)*4.0;
  }

  /* Force the step to be within the bounds stpmax and stpmin. */
  *stp_ptr = MAX(*stp_ptr, ws->super.stpmin);
  *stp_ptr = MIN(*stp_ptr, ws->super.stpmax);

  /* If further progress is not possible, let stp be the best
     point obtained during the search. */
  if (ws->brackt && (*stp_ptr <= ws->stmin || *stp_ptr >= ws->stmax
                     || ws->stmax - ws->stmin <= ws->xtol * ws->stmax)) {
    *stp_ptr = ws->stx;
  }

  /* Obtain another function and derivative. */
  return OPK_LNSRCH_SEARCH;
}

static csrch_workspace_t*
csrch_get_workspace(opk_lnsrch_workspace_t* ws)
{
  if (ws->start == csrch_start && ws->iterate == csrch_iterate) {
    return (csrch_workspace_t*)ws;
  } else {
    return NULL;
  }
}

opk_lnsrch_workspace_t*
opk_lnsrch_new_csrch(double ftol, double gtol, double xtol)
{
  opk_lnsrch_workspace_t* _ws;

  if (ftol < 0.0) {
    /* ERROR: FTOL .LT. ZERO */
    errno = EINVAL;
    return NULL;
  }
  if (gtol < 0.0) {
    /* ERROR: GTOL .LT. ZERO */
    errno = EINVAL;
    return NULL;
  }
  if (xtol < 0.0) {
    /* ERROR: XTOL .LT. ZERO */
    errno = EINVAL;
    return NULL;
  }
  _ws = _opk_lnsrch_new(sizeof(csrch_workspace_t),
			csrch_start, csrch_iterate, NULL);
  if (_ws != NULL) {
    csrch_workspace_t* ws = (csrch_workspace_t*)_ws;
    ws->ftol = ftol;
    ws->gtol = gtol;
    ws->xtol = xtol;
    ws->stage = 0;
  }
  return _ws;
}

/**
 * Find a step that satisfies a sufficient decrease condition and a curvature
 * condition.
 *
 * This subroutine finds a step that satisfies a sufficient decrease condition
 * and a curvature condition.
 *
 * Each call of the subroutine updates an interval with endpoints stx and
 * sty. The interval is initially chosen so that it contains a minimizer of the
 * modified function
 *
 *	     psi(stp) = f(stp) - f(0) - ftol*stp*f'(0).
 *
 * If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the interval is chosen
 * so that it contains a minimizer of f.
 *
 * The algorithm is designed to find a step that satisfies the sufficient
 * decrease condition
 *
 *	     f(stp) <= f(0) + ftol*stp*f'(0),
 *
 * and the curvature condition
 *
 *	     abs(f'(stp)) <= gtol*abs(f'(0)).
 *
 * If ftol is less than gtol and if, for example, the function is bounded
 * below, then there is always a step which satisfies both conditions.
 *
 * If no step can be found that satisfies both conditions, then the algorithm
 * stops with a warning. In this case stp only satisfies the sufficient
 * decrease condition.
 *
 * A typical invocation of dcsrch has the following outline:
 *
 *     Evaluate the function at stp = 0.0d0; store in f.
 *     Evaluate the gradient at stp = 0.0d0; store in g.
 *     Choose a starting step stp.
 *
 *     task = 'START'
 *  10 continue
 *	  call dcsrch(stp,f,g,ftol,gtol,xtol,task,stpmin,stpmax,
 *    +		      isave,dsave)
 *	  if (task .eq. 'FG') then
 *	     Evaluate the function and the gradient at stp
 *	     go to 10
 *	     end if
 *
 *     NOTE: The user must not alter work arrays between calls.
 *
 *     The subroutine statement is
 *
 *	 subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,
 *			   task,isave,dsave)
 *     where
 *
 *	 stp is a double precision variable.
 *	   On entry stp is the current estimate of a satisfactory
 *	      step. On initial entry, a positive initial estimate
 *	      must be provided.
 *	   On exit stp is the current estimate of a satisfactory step
 *	      if task = 'FG'. If task = 'CONV' then stp satisfies
 *	      the sufficient decrease and curvature condition.
 *
 *	 f is a double precision variable.
 *	   On initial entry f is the value of the function at 0.
 *	      On subsequent entries f is the value of the
 *	      function at stp.
 *	   On exit f is the value of the function at stp.
 *
 *	 g is a double precision variable.
 *	   On initial entry g is the derivative of the function at 0.
 *	      On subsequent entries g is the derivative of the
 *	      function at stp.
 *	   On exit g is the derivative of the function at stp.
 *
 *	 ftol is a double precision variable.
 *	   On entry ftol specifies a nonnegative tolerance for the
 *	      sufficient decrease condition.
 *	   On exit ftol is unchanged.
 *
 *	 gtol is a double precision variable.
 *	   On entry gtol specifies a nonnegative tolerance for the
 *	      curvature condition.
 *	   On exit gtol is unchanged.
 *
 *	 xtol is a double precision variable.
 *	   On entry xtol specifies a nonnegative relative tolerance
 *	      for an acceptable step. The subroutine exits with a
 *	      warning if the relative difference between sty and stx
 *	      is less than xtol.
 *	   On exit xtol is unchanged.
 *
 *	 task is a character variable of length at least 60.
 *	   On initial entry task must be set to 'START'.
 *	   On exit task indicates the required action:
 *
 *	      If task(1:2) = 'FG' then evaluate the function and
 *	      derivative at stp and call dcsrch again.
 *
 *	      If task(1:4) = 'CONV' then the search is successful.
 *
 *	      If task(1:4) = 'WARN' then the subroutine is not able
 *	      to satisfy the convergence conditions. The exit value of
 *	      stp contains the best point found during the search.
 *
 *	      If task(1:5) = 'ERROR' then there is an error in the
 *	      input arguments.
 *
 *	   On exit with convergence, a warning or an error, the
 *	      variable task contains additional information.
 *
 *	 stpmin is a double precision variable.
 *	   On entry stpmin is a nonnegative lower bound for the step.
 *	   On exit stpmin is unchanged.
 *
 *	 stpmax is a double precision variable.
 *	   On entry stpmax is a nonnegative upper bound for the step.
 *	   On exit stpmax is unchanged.
 *
 *	 isave is an integer work array of dimension 2.
 *
 *	 dsave is a double precision work array of dimension 13.
 *
 *     Subprograms called
 *
 *	 MINPACK-2 ... dcstep
 *
 *     MINPACK-1 Project. June 1983.
 *     Argonne National Laboratory.
 *     Jorge J. More' and David J. Thuente.
 *
 *     MINPACK-2 Project. November 1993.
 *     Argonne National Laboratory and University of Minnesota.
 *     Brett M. Averick, Richard G. Carter, and Jorge J. More'.
 */

/*---------------------------------------------------------------------------*/

static double max3(double val1, double val2, double val3)
{
  double result = val1;
  if (val2 > result) result = val2;
  if (val3 > result) result = val3;
  return result;
}

/* FIXME: add a reference to the documentation and check history. */

/**
 * Compute a safeguarded step for a line search.
 *
 * This routine computes a safeguarded step for a search procedure and updates
 * an interval that contains a step that satisfies a sufficient decrease and a
 * curvature condition.
 *
 * The parameter stx contains the step with the least function value. If
 * brackt is true, then a minimizer has been bracketed in an interval
 * with endpoints stx and sty.  The parameter stp contains the current step.
 * The subroutine assumes that if brackt is true then:
 *
 *       min(stx,sty) < stp < max(stx,sty),
 *
 * and that the derivative at stx is negative in the direction of the step;
 * that is:
 *
 *      dx*(stp - stx) < 0
 *
 *
 * @param stx_ptr
 *        A pointer to stx.  On entry, stx is the best step obtained so far
 *        and is an endpoint of the interval that contains the minimizer.  On
 *        exit, stx is the updated best step.
 *
 * @param fx_ptr
 *        A pointer to fx.  On entry, fx is the function at stx.  On exit, fx
 *        is the function at stx.
 *
 * @param dx_ptr
 *        A pointer to dx.  On entry, dx is the derivative of the function at
 *        stx. The derivative must be negative in the direction of the step,
 *        that is, dx and stp - stx must have opposite signs.  On exit, dx is
 *        the derivative of the function at stx.
 *
 * @param sty_ptr
 *        A pointer to sty.  On entry, sty is the second endpoint of the
 *        interval that contains the minimizer.  On exit, sty is the updated
 *        endpoint of the interval that contains the minimizer.
 *
 * @param fy_ptr
 *        A pointer to fy.  On entry, fy is the function at sty.  On exit, fy
 *        is the function at sty.
 *
 * @param dy_ptr
 *        A pointer to dy.  On entry, dy is the derivative of the function at
 *        sty.  On exit, dy is the derivative of the function at sty.
 *
 * @param stp_ptr
 *        A pointer to stp. On entry, stp is the current step. If the value at
 *	  brackt_ptr is true, then on input, stp must be between stx and sty.
 *	  On exit, stp is a new trial step.
 *
 * @param fp
 *	  The value of the function at stp on entry.
 *
 * @param dp
 *	  The the derivative of the function at stp on entry.
 *
 * @param brackt_ptr
 *	  A pointer to the boolean variable brackt.  On entry, brackt
 *	  specifies if a minimizer has been bracketed.  Initially brackt must
 *	  be set to false On exit, brackt specifies if a minimizer has been
 *	  bracketed.  When a minimizer is bracketed brackt is set to true.
 *
 * @return A strictly negative value on error, a strictly positive value on
 * success.  The value returned on success is between 1 and 4 and corresponds
 * to one of the four possible cases.  The value returned on error is one of:
 *
 *    OPK_LNSRCH_ERROR_STP_OUTSIDE_BRACKET
 *    if brackt is true but the step is outside the bracket endpoints.
 *
 *    OPK_LNSRCH_ERROR_NOT_A_DESCENT
 *    if the descent condition is violated.
 *
 *    OPK_LNSRCH_ERROR_STPMIN_GT_STPMAX
 *    if stpmin > stpmax.
 *
 * @history
 * MINPACK-1 Project. June 1983
 * Argonne National Laboratory.
 * Jorge J. Moré and David J. Thuente.
 *
 * MINPACK-2 Project. November 1993.
 * Argonne National Laboratory and University of Minnesota.
 * Brett M. Averick and Jorge J. Moré.
 */
int opk_cstep(double* stx_ptr, double* fx_ptr, double* dx_ptr,
	      double* sty_ptr, double* fy_ptr, double* dy_ptr,
	      double* stp_ptr, double  fp,     double  dp,
	      int* brackt_ptr, double stpmin, double stpmax)
{
  /* Constants. */
  const double ZERO = 0.0;
  const double TWO = 2.0;
  const double THREE = 3.0;

  /* Get values of input/output variables. */
  double stx = *stx_ptr, fx = *fx_ptr, dx = *dx_ptr;
  double sty = *sty_ptr, fy = *fy_ptr, dy = *dy_ptr;
  double stp = *stp_ptr;

  /* Local variables. */
  double gamma, theta, p, q, r, s, temp;
  double stpc; /* cubic step */
  double stpq; /* quadratic step */
  double stpf;
  int opposite, result;

  /* Check the input parameters for errors. */
  if ((*brackt_ptr) && (stx < sty ? (stp <= stx || stp >= sty)
		                  : (stp <= sty || stp >= stx))) {
    return OPK_LNSRCH_ERROR_STP_OUTSIDE_BRACKET;
  } else if (dx*(stp - stx) >= ZERO) {
    return OPK_LNSRCH_ERROR_NOT_A_DESCENT;
  } else if (stpmin > stpmax) {
    return OPK_LNSRCH_ERROR_STPMIN_GT_STPMAX;
  }

  /* Determine if the derivatives have opposite signs. */
  opposite = ((dp < ZERO && dx > ZERO) || (dp > ZERO && dx < ZERO));

  if (fp > fx) {
    /* First case.  A higher function value.  The minimum is bracketed.  If
       the cubic step is closer to STX than the quadratic step, the cubic step
       is taken, otherwise the average of the cubic and quadratic steps is
       taken. */
    result = 1;
    *brackt_ptr = TRUE;
    theta = THREE*(fx - fp)/(stp - stx) + dx + dp;
    s = max3(fabs(theta), fabs(dx), fabs(dp));
    temp = theta/s;
    gamma = s*sqrt(temp*temp - (dx/s)*(dp/s));
    if (stp < stx) gamma = -gamma;
    p =  (gamma - dx) + theta;
    q = ((gamma - dx) + gamma) + dp;
    r = p/q;
    stpc = stx + r*(stp - stx);
    stpq = stx + ((dx/((fx - fp)/(stp - stx) + dx))/TWO)*(stp - stx);
    if (fabs(stpc - stx) < fabs(stpq - stx)) {
      stpf = stpc;
    } else {
      stpf = stpc + (stpq - stpc)/TWO;
    }
  } else if (opposite) {
    /* Second case.  A lower function value and derivatives of opposite sign.
       The minimum is bracketed.  If the cubic step is farther from STP than
       the secant (quadratic) step, the cubic step is taken, otherwise the
       secant step is taken. */
    result = 2;
    *brackt_ptr = TRUE;
    theta = THREE*(fx - fp)/(stp - stx) + dx + dp;
    s = max3(fabs(theta), fabs(dx), fabs(dp));
    temp = theta/s;
    gamma = s*sqrt(temp*temp - (dx/s)*(dp/s));
    if (stp > stx) gamma = -gamma;
    p =  (gamma - dp) + theta;
    q = ((gamma - dp) + gamma) + dx;
    r = p/q;
    stpc = stp + r*(stx - stp);
    stpq = stp + (dp/(dp - dx))*(stx - stp);
    if (fabs(stpc - stp) > fabs(stpq - stp)) {
      stpf = stpc;
    } else {
      stpf = stpq;
    }
  } else if (fabs(dp) < fabs(dx)) {
    /* Third case.  A lower function value, derivatives of the same sign, and
       the magnitude of the derivative decreases.  The cubic step is computed
       only if the cubic tends to infinity in the direction of the step or if
       the minimum of the cubic is beyond STP.  Otherwise the cubic step is
       defined to be the secant step.  The case GAMMA = 0 only arises if the
       cubic does not tend to infinity in the direction of the step. */
    result = 3;
    theta = THREE*(fx - fp)/(stp - stx) + dx + dp;
    s = max3(fabs(theta), fabs(dx), fabs(dp));
    temp = theta/s;
    temp = temp*temp - (dx/s)*(dp/s);
    if (temp > ZERO) {
      gamma = s*sqrt(temp);
      if (stp > stx) gamma = -gamma;
    } else {
      gamma = ZERO;
    }
    p = (gamma - dp) + theta;
    q = (gamma + (dx - dp)) + gamma;
    r = p/q;
    if (r < ZERO && gamma != ZERO) {
      stpc = stp + r*(stx - stp);
    } else if (stp > stx) {
      stpc = stpmax;
    } else {
      stpc = stpmin;
    }
    stpq = stp + (dp/(dp - dx))*(stx - stp);

    if (*brackt_ptr) {
      /* A minimizer has been bracketed.  If the cubic step is closer to STP
	 than the secant step, the cubic step is taken, otherwise the secant
	 step is taken. */
      if (fabs(stpc - stp) < fabs(stpq - stp)) {
	stpf = stpc;
      } else {
	stpf = stpq;
      }
      temp = stp + 0.66*(sty - stp);
      if (stp > stx ? stpf > temp : stpf < temp) {
	stpf = temp;
      }
    } else {
      /* A minimizer has not been bracketed. If the cubic step is farther from
	 stp than the secant step, the cubic step is taken, otherwise the
	 secant step is taken. */
      if (fabs(stpc - stp) > fabs(stpq - stp)) {
	stpf = stpc;
      } else {
	stpf = stpq;
      }
      if (stpf > stpmax) stpf = stpmax;
      if (stpf < stpmin) stpf = stpmin;
    }
  } else {
    /* Fourth case.  A lower function value, derivatives of the same sign, and
       the magnitude of the derivative does not decrease.  If the minimum is
       not bracketed, the step is either STPMIN or STPMAX, otherwise the cubic
       step is taken. */
    result = 4;
    if (*brackt_ptr) {
      theta = THREE*(fp - fy)/(sty - stp) + dy + dp;
      s = max3(fabs(theta), fabs(dy), fabs(dp));
      temp = theta/s;
      gamma = s*sqrt(temp*temp - (dy/s)*(dp/s));
      if (stp > sty) gamma = -gamma;
      p =  (gamma - dp) + theta;
      q = ((gamma - dp) + gamma) + dy;
      r = p/q;
      stpc = stp + r*(sty - stp);
      stpf = stpc;
    } else if (stp > stx) {
      stpf = stpmax;
    } else {
      stpf = stpmin;
    }
  }

  /* Update the interval which contains a minimizer. */
  if (fp > fx) {
    *sty_ptr = stp;
    *fy_ptr = fp;
    *dy_ptr = dp;
  } else {
    if (opposite) {
      *sty_ptr = stx;
      *fy_ptr = fx;
      *dy_ptr = dx;
    }
    *stx_ptr = stp;
    *fx_ptr = fp;
    *dx_ptr = dp;
  }

  /* Store the new safeguarded step. */
  *stp_ptr = stpf;
  return result;
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
