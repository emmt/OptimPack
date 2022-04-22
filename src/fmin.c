/*
 * fmin.c --
 *
 * Minimization of an univariate function.  The implemented algorithms are the
 * golden section method and Brent's method [1] with some improvements to
 * account for various boundary conditions and initial search interval.  A
 * reverse communication version is also provided.
 *
 * [1] Brent, R.P., "Algorithms for Minimization without Derivatives,"
 *     Prentice-Hall, Inc. (1973).
 *
 * See code at: http://netlib.org/go/fmin.f
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of OptimPack (https://github.com/emmt/OptimPack).
 *
 * Copyright (C) 2009-2019 Éric Thiébaut
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

#ifndef OPK_FMIN_C_
#define OPK_FMIN_C_ 1

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <float.h>

#include "optimpack.h"

#define OPK_FMIN_BOUNDED (OPK_FMIN_BOUNDED_BY_A | OPK_FMIN_BOUNDED_BY_B)

/* The HUGE_VAL macro (defined in math.h) is an expression representing a
   particular very large number. On machines that use IEEE floating point
   format, the value is "infinity". On other machines, it's typically the
   largest positive number that can be represented. */

static const double FMIN_EPSILON = DBL_EPSILON;
static const double FMIN_INFINITY = HUGE_VAL;
static double FMIN_SQRT_EPSILON_ = -1.0;

#define FMIN_SQRT_EPSILON ((FMIN_SQRT_EPSILON_ > 0.0) ? FMIN_SQRT_EPSILON_ \
                         : (FMIN_SQRT_EPSILON_ = sqrt(FMIN_EPSILON)))

#define SWAP(a, b, temp)  temp = a; a = b; b = temp

/* Golden section constants.
 *   ALPHA = (3 - sqrt(5))/2
 *   BETA  = 1 - ALPHA = (sqrt(5) - 1)/2
 *   GAMMA = BETA/ALPHA = (sqrt(5) + 1)/2
 */
static const double FMIN_ALPHA = 3.81966011250105151795413165634362E-1;
#if 0 /* unused stuff */
static const double FMIN_BETA  = 6.18033988749894848204586834365638E-1;
#endif
static const double FMIN_GAMMA = 1.61803398874989484820458683436564E+0;

typedef struct {
  double a, b, u, fu, v, fv, w, fw, x;
  double prec;
  long nevals;
  const char *msg;
  unsigned int flags;
  opk_fmin_task status;
  int stage;
} opk_fmin_context;

/* DOCUMENT FMIN_EPSILON  = the smallest positive floating point x such
                            that 1 + x is numerically not equal to 1.
            FMIN_HUGE     = the largest representable finite floating-point
                            number.
            FMIN_INFINITY = the floating-point representation of infinite.
            FMIN_SQRT_EPSILON = the square root of FMIN_EPSILON.

   SEE ALSO: machine_constant, fmin_golden, fmin_brent.
 */

opk_fmin_context*
opk_fmin_new(int method)
{
  opk_fmin_context* ctx;
  ctx = (opk_fmin_context*)malloc(sizeof(opk_fmin_context));
  if (ctx != NULL) {
    memset(ctx, 0, sizeof(opk_fmin_context));
    ctx->status = OPK_FMIN_START - 1;
    ctx->prec = FMIN_SQRT_EPSILON;
  }
  return ctx;
}

void
opk_fmin_destroy(opk_fmin_context* ctx)
{
  if (ctx != NULL) {
    free((void*)ctx);
  }
}

opk_status
opk_fmin_set_precision(opk_fmin_context* ctx, double prec)
{
  if (ctx == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  if (prec < 0) {
    return OPK_INVALID_ARGUMENT;
  }
  ctx->prec = prec;
  return OPK_SUCCESS;
}

double
opk_fmin_get_precision(opk_fmin_context* ctx)
{
  if (ctx != NULL) return ctx->prec;
  errno = EFAULT;
  return FMIN_SQRT_EPSILON;
}

unsigned int
opk_fmin_get_flags(opk_fmin_context* ctx)
{
  if (ctx != NULL) return ctx->flags;
  errno = EFAULT;
  return 0;
}

opk_fmin_task
opk_fmin_get_status(opk_fmin_context* ctx)
{
  if (ctx != NULL) return ctx->status;
  errno = EFAULT;
  return OPK_FMIN_ERROR;
}

opk_status
opk_fmin_start(opk_fmin_context* ctx,
               double a, double b,
               unsigned int flags)
{
  if (ctx == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  if (a == b) {
    return OPK_INVALID_ARGUMENT;
  }
  ctx->a = a;
  ctx->b = b;
  ctx->nevals = 0;
  ctx->flags = (flags & (OPK_FMIN_BOUNDED | OPK_FMIN_SMOOTH));
  ctx->stage = 0;
  ctx->status = OPK_FMIN_START;
  return OPK_SUCCESS;
}

#define MAYBE_GET(ptr, expr)  if (ptr != NULL) *ptr = expr

opk_status
opk_fmin_get_initial(opk_fmin_context* ctx,
                     double* a, double* b,
                     unsigned int* flags)
{
  if (ctx == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  if (ctx->status < OPK_FMIN_START) {
    return OPK_INVALID_ARGUMENT;
  }
  MAYBE_GET(a, ctx->a);
  MAYBE_GET(b, ctx->b);
  MAYBE_GET(flags, ctx->flags);
  return OPK_SUCCESS;
}

opk_status
opk_fmin_get_final(opk_fmin_context* ctx,
                   double* xmin, double* xlo, double* xup,
                   double* fmin, double* flo, double* fup,
                   long* nevals)
{
  if (ctx == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  /* FIXME: bracket? */
  if (ctx->status < OPK_FMIN_START) {
    return OPK_INVALID_ARGUMENT;
  }
  MAYBE_GET(xmin, ctx->v);
  MAYBE_GET(fmin, ctx->fv);
  if (ctx->u > ctx->w) {
    MAYBE_GET(xlo, ctx->w);
    MAYBE_GET(flo, ctx->fw);
    MAYBE_GET(xup, ctx->v);
    MAYBE_GET(fup, ctx->fv);
  } else {
    MAYBE_GET(xlo, ctx->v);
    MAYBE_GET(flo, ctx->fv);
    MAYBE_GET(xup, ctx->w);
    MAYBE_GET(fup, ctx->fw);
  }
  MAYBE_GET(nevals, ctx->nevals);
  return OPK_SUCCESS;
}

/*
 * For golden section algorithm:
 *    ctx->stage = -1    error or not yet initialized
 *    ctx->stage =  0    initial interval set
 *    ctx->stage =  1    evaluation of f(a)
 *    ctx->stage =  2    evaluation of f(b)
 *    ctx->stage =  3    evaluation of f(x) for x inside (a,b)
 *    ctx->stage =  4    evaluation of f(x) while trying to bracket
 *    ctx->stage =  5    evaluation of f(x) while converging
 */

opk_fmin_task
opk_fmin_next(opk_fmin_context* ctx, double* xptr, double fx)
{
#define a     (ctx->a)
#define b     (ctx->b)
#define u     (ctx->u)
#define v     (ctx->v)
#define w     (ctx->w)
#define x     (ctx->x)
#define fu    (ctx->fu)
#define fv    (ctx->fv)
#define fw    (ctx->fw)
#define prec  (ctx->prec)
#define flags (ctx->flags)

  double r, t;
  /* int brent; */

  if (ctx == NULL) {
    errno = EFAULT;
    return OPK_FMIN_ERROR;
  }
  if (xptr == NULL) {
    errno = EFAULT;
    return (ctx->status = OPK_FMIN_ERROR);
  }

  /* The following macros are to facilitate reverse communication and to hide
     the necessary ugliness of the resulting code.  S is the stage number, N is
     the member name for the variable value and EXPR is the expression
     (evaluated only once) to compute the variable value. */
#define _EVAL_u(s, expr)    _EVAL_n(s, u, expr)
#define _EVAL_v(s, expr)    _EVAL_n(s, v, expr)
#define _EVAL_w(s, expr)    _EVAL_n(s, w, expr)
#define _EVAL_x(s, expr)    _EVAL_o(s, *xptr = x = (expr))
#define _EVAL_n(s, n, expr) _EVAL_o(s, *xptr = x = n = (expr)); f##n = fx
#define _EVAL_o(s, init)                        \
     init;                                      \
     ctx->stage = s;                            \
     return (ctx->status = OPK_FMIN_FX);        \
   stage_##s:                                   \
     ++ctx->nevals
#define EVAL(s, n, expr)   _EVAL_##n(s, expr)
#define JUMP(s)            case s: goto stage_##s

  if ((ctx->status == OPK_FMIN_FX) ||
      (ctx->status == OPK_FMIN_CONVERGENCE)) {
    /* brent = ((flags & OPK_FMIN_SMOOTH) != 0); */
    if (*xptr != x) {
      /* The workspace has been corrupted or the caller has changed the value
         of the variable. */
      return (ctx->status = OPK_FMIN_ERROR);
    }
    switch (ctx->stage) {
      JUMP(1);
      JUMP(2);
      JUMP(3);
      JUMP(4);
      JUMP(5);
    default:
      return (ctx->status = OPK_FMIN_ERROR);
    }
  }

#undef JUMP

  if (ctx->status != OPK_FMIN_START) {
    return (ctx->status = OPK_FMIN_ERROR);
  }

  /* Evaluation of function at U = A. */
  if ((flags & OPK_FMIN_BOUNDED_BY_A) != 0) {
    u = a;
    fu = FMIN_INFINITY;
  } else {
    EVAL(1, u, a);
  }

  /* Evaluation of function at W = B. */
  if ((flags & OPK_FMIN_BOUNDED_BY_B) != 0) {
    w = b;
    fw = FMIN_INFINITY;
  } else {
    EVAL(2, w, b);
  }

  /* Take a golden section interpolation step (choose a point which is closer
     to the end-point with the highest function value; that is, U in our
     settings). */
  if (fw > fu) {
    /* Make sure that f(u) >= f(w). */
    SWAP(u,  w,  t);
    SWAP(fu, fw, t);
  }
  EVAL(3, v, FMIN_ALPHA*(w - u) + u);

  /* Take a golden section extrapolation step until a bracket is found. */
  while (fv > fw) {
    u = v;
    fu = fv;
    v = w;
    fv = fw;
    EVAL(4, w, FMIN_GAMMA*(v - u) + v);
  }

  /* Reduce the bracketing interval to improve the precision. */
  for (;;) {
    /* Check for convergence. */
    r = fabs(u);
    if ((t = fabs(v)) > r) r = t;
    if ((t = fabs(w)) > r) r = t;
    if (fabs(w - v) <= prec*r) {
      return (ctx->status = OPK_FMIN_CONVERGENCE);
    }

    /* Take a golden section interpolation step. */
    EVAL(5, x, FMIN_ALPHA*(w - v) + v);

    /* Update the interval.  In case of a tie (fx = fv and fu = fw), the
       situation is symmetrical and the choice arbitrary. */
    if ((fx > fv) || ((fx == fv) && (fu >= fw))) { /* FIXME: check this secondary rule. */
      /* The new interval is (x,v,u). */
      w = u;
      fw = fu;
      u = x;
      fu = fx;
    } else {
      /* The new interval is (v,x,w). */
      u = v;
      fu = fv;
      v = x;
      fv = fx;
    }
  }

#undef a
#undef b
#undef u
#undef v
#undef w
#undef x
#undef fu
#undef fv
#undef fw
#undef prec
#undef flags
}

/* Usage:


   opk_fmin(f, a, b, flags, prec, maxeval)

   if OPK_FMIN_DATA is set then f(data, x) yields f(x)
   otherwise f(x)


   Reverse callback interface:

   double x, a, b;
   opk_fmin_context* c = opk_fmin_new();
   int stage;

   if (opk_fmin_start(c, a, b, flags) != OPK_SUCCESS) {
     error...
   }
   for (;;) {
     int stage = opk_fmin_next_stage(c, &x, fx);
     if (stage == OPK_FMIN_FX) {
       fx = ....;
     } else if (stage == OPK_FMIN_CONVERGENCE) {
       break;
     } else if (stage == OPK_FMIN_ERROR) {
       error...;
       break;
     }
   }


   opk_fmin_destroy(c);


   Any number of minimizations possible (use opk_fmin_start to initiate the
   search).

   Continuation possible : keep calling  opk_fmin_next_stage after
   OPK_FMIN_CONVERGENCE

   Some parameters may be changed on the fly:

      - relative precision
      - maximum number of evaluations

   Some parameters can only be changed right after opk_fmin_new or right after
   opk_fmin_start:

      - optimization method
      - starting interval

   Parameters that can be queried:

      - number of function evaluations
      - initial interval (flags and values)
      - optimization method
      - (relative) width of curent interval
      - parameters of current interval (u,fu,v,fv,w,fw and bracket status)
      - requested precision
      - maximum number of function evaluations
      - value of the parameter to compute the function
*/


#define OPK_FMIN_WITH_CONTEXT 0
#include __FILE__


#define OPK_FMIN_WITH_CONTEXT 1
#include __FILE__

#if 0
double fmin_test_f(double x)
{
  /* This function has 2 minima: a smooth one at
   * (123 - sqrt(24141))/30 ~ -1.07912 and
   * a non-smooth one at 12.7.
   */
  fx = fabs(x - 12.7)*(x + 5.2)*(x - 4.8);
  return fx;
}
#endif

#else /* OPK_FMIN_C_ *********************************************************/

#ifdef OPK_FMIN_WITH_CONTEXT

#if OPK_FMIN_WITH_CONTEXT
# define func  opk_fmin_with_context
# define args  double (*f_)(void* data, double x), double a, double b, \
               unsigned int flags, long maxeval, double prec, double out[7], \
               void* data
# define f(x)  f_(data, (x))
#else
# define func  opk_fmin
# define args  double (*f_)(double x), double a, double b, \
               unsigned int flags, long maxeval, double prec, double out[7]
# define f(x)  f_(x)
#endif

int func(args)
{
  double r, t, u, fu, v, fv, w, fw, x, fx;
  long nevals;
  int status;

  status = -1;
  if (prec < 0.0) {
    prec = FMIN_SQRT_EPSILON;
  }
  if (a == b) {
    /* A and B must be different. */
    errno = EINVAL;
    return -1;
  }

  /* Search a bracketing of the minimum. */
  nevals = 0;
  u = a;
  if ((flags & OPK_FMIN_BOUNDED_BY_A) == 0) {
    fu = f(u);
    ++nevals;
  } else {
    fu = FMIN_INFINITY;
  }
  w = b;
  if ((flags & OPK_FMIN_BOUNDED_BY_B) == 0) {
    fw = f(w);
    ++nevals;
  } else {
    fw = FMIN_INFINITY;
  }
  if (fw > fu) {
    /* Make sure that f(u) >= f(w). */
    t =  u;  u =  w;  w = t;
    t = fu; fu = fw; fw = t;
  }
  v = FMIN_ALPHA*(w - u) + u; /* Golden section interpolation step: choose a
                                 point which is closer to the end-point with
                                 the highest function value (that is, U in our
                                 settings). */
  fv = f(v);
  ++nevals;
  while (fv > fw) {
    if (maxeval > 0 && nevals >= maxeval) {
      status = 2; /* no bracket yet */
      goto done;
    }
    u = v;
    fu = fv;
    v = w;
    fv = fw;
    w = FMIN_GAMMA*(v - u) + v; /* Golden section extrapolation step. */
    fw = f(w);
    ++nevals;
  }

  /* Reduce the bracketing interval to improve the precision. */
  for (;;) {
    /* Check for convergence. */
    r = fabs(u);
    if ((t = fabs(v)) > r) r = t;
    if ((t = fabs(w)) > r) r = t;
    if (fabs(w - v) <= prec*r) {
      status = 0; /* convergence */
      goto done;
    } else if (maxeval > 0 && nevals >= maxeval) {
      status = 1; /* bracket but no convergence */
      goto done;
    }

    /* Take a golden section interpolation step. */
    x = FMIN_ALPHA*(w - v) + v;
#if 0
    /* FIXME: This test may be used to save one function evaluation when
       extreme precision is requested and rounding errors prevent further
       progresses. */
    if (x == u || x == v || x == w) {
      break;
    }
#endif
    fx = f(x);
    ++nevals;
    /* Update the interval.  In case of a tie (fx = fv and fu = fw), the
       situation is symmetrical and the choice arbitrary. */
    if ((fx > fv) || ((fx == fv) && (fu >= fw))) { /* FIXME: check this secondary rule. */
      /* The new interval is (x,v,u). */
      w = u;
      fw = fu;
      u = x;
      fu = fx;
    } else {
      /* The new interval is (v,x,w). */
      u = v;
      fu = fv;
      v = x;
      fv = fx;
    }
  }

  /* Return the result. */
 done:
  if (u > w) {
    out[0] = v;
    out[1] = w;
    out[2] = u;
    out[3] = fv;
    out[4] = fw;
    out[5] = fu;
  } else {
    out[0] = v;
    out[1] = u;
    out[2] = w;
    out[3] = fv;
    out[4] = fu;
    out[5] = fw;
  }
  out[6] = nevals;
  return status;
}

#undef func
#undef args
#undef f
#undef OPK_FMIN_WITH_CONTEXT

#endif /* OPK_FMIN_WITH_CONTEXT */

#endif /* OPK_FMIN_C_ */
