/*
 * simple_vector.c --
 *
 * Simple vector space implementation for OptimPack library.
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

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include "optimpack-private.h"

#define ROUND_UP(a,b)   OPK_ROUND_UP(a,b)


#ifndef SINGLE_PRECISION
#  error expecting macro SINGLE_PRECISION to be defined
#endif

#if SINGLE_PRECISION
#  define REAL     float
#  define ABS(x)   fabsf(x)
#  define SQRT(x)  sqrtf(x)
#  define ALPHA    _alpha
#  define BETA     _beta
#  define GAMMA    _gamma
#  define ZERO     0.0f
#  define ONE      1.0f
#  define FLOAT_CHOICE(a,b) a
#else
#  define REAL     double
#  define ABS(x)   fabs(x)
#  define SQRT(x)  sqrt(x)
#  define ALPHA    alpha
#  define BETA     beta
#  define GAMMA    gamma
#  define ZERO     0.0
#  define ONE      1.0
#  define FLOAT_CHOICE(a,b) b
#endif

typedef struct _simple_vector simple_vector_t;
struct _simple_vector {
  opk_vector_t base;  /* base type (must be the first member) */
  REAL* data;
  void* client_data;
  void (*free_client_data)(void* client_data);
};

#define DATA(v) ((simple_vector_t*)(v))->data

static opk_vector_t*
create(opk_vspace_t* vspace)
{
  const size_t offset = ROUND_UP(sizeof(simple_vector_t), sizeof(REAL));
  size_t size = offset + vspace->size*sizeof(REAL);
  opk_vector_t* v = opk_allocate_vector(vspace, size);
  if (v != NULL) {
    simple_vector_t* sv = (simple_vector_t*)v;
    sv->data = (REAL*)(((char*)v) + offset);
    sv->client_data = NULL;
    sv->free_client_data = NULL;
  }
  return v;
}

static void
finalize(opk_vspace_t* vspace,
         opk_vector_t* v)
{
  simple_vector_t* sv = (simple_vector_t*)v;
  if (sv->free_client_data != NULL) {
    sv->free_client_data(sv->client_data);
  }
}

static void
fill(opk_vspace_t* vspace, opk_vector_t* vect, double alpha)
{
  REAL* x = DATA(vect);
  opk_index_t j, n = vspace->size;
  if (alpha == 0.0) {
    memset(x, 0, n*sizeof(REAL));
  } else {
#if SINGLE_PRECISION
    REAL ALPHA = (REAL)alpha;
#endif
    for (j = 0; j < n; ++j) {
      x[j] = ALPHA;
    }
  }
}

static double
norm1(opk_vspace_t* vspace,
      const opk_vector_t* vx)
{
  REAL result = ZERO;
  const REAL* x = DATA(vx);
  opk_index_t j, n = vspace->size;
  for (j = 0; j < n; ++j) {
    result += ABS(x[j]);
  }
  return (double)result;
}

static double
norm2(opk_vspace_t* vspace,
      const opk_vector_t* vx)
{
  REAL result = ZERO;
  const REAL* x = DATA(vx);
  opk_index_t j, n = vspace->size;
  for (j = 0; j < n; ++j) {
    REAL xj = x[j];
    result += xj*xj;
  }
  return (double)SQRT(result);
}

static double
norminf(opk_vspace_t* vspace,
        const opk_vector_t* vx)
{
  REAL result = ZERO;
  const REAL* x = DATA(vx);
  opk_index_t j, n = vspace->size;
  for (j = 0; j < n; ++j) {
    REAL axj = ABS(x[j]);
    if (axj > result) {
      result = axj;
    }
  }
  return (double)result;
}

static double
dot(opk_vspace_t* vspace,
    const opk_vector_t* vx,
    const opk_vector_t* vy)
{
  REAL result = ZERO;
  const REAL* x = DATA(vx);
  const REAL* y = DATA(vy);
  opk_index_t j, n = vspace->size;
  for (j = 0; j < n; ++j) {
    result += x[j]*y[j];
  }
  return (double)result;
}

static double
dot3(opk_vspace_t* vspace,
     const opk_vector_t* vw,
     const opk_vector_t* vx,
     const opk_vector_t* vy)
{
  REAL result = ZERO;
  const REAL* w = DATA(vw);
  const REAL* x = DATA(vx);
  const REAL* y = DATA(vy);
  opk_index_t j, n = vspace->size;
  for (j = 0; j < n; ++j) {
    result += w[j]*x[j]*y[j];
  }
  return (double)result;
}

static void
copy(opk_vspace_t* vspace,
     opk_vector_t* vdst, const opk_vector_t* vsrc)
{
  REAL* dst = DATA(vdst);
  const REAL* src = DATA(vsrc);
  if (dst != src) {
    memcpy(dst, src, vspace->size*sizeof(REAL));
  }
}

static void
swap(opk_vspace_t* vspace,
     opk_vector_t* vx, opk_vector_t* vy)
{
  REAL* x = DATA(vx);
  REAL* y = DATA(vy);
  if (x != y) {
    opk_index_t j, n = vspace->size;
    for (j = 0; j < n; ++j) {
      REAL t = x[j];
      x[j] = y[j];
      y[j] = t;
    }
  }
}

static void
scale(opk_vspace_t* vspace, opk_vector_t* vdst,
      double alpha, const opk_vector_t* vsrc)
{
  /* Note: we already know that ALPHA is neither 0 nor 1. */
  REAL* dst = DATA(vdst);
  const REAL* src = DATA(vsrc);
  opk_index_t j, n = vspace->size;
#if SINGLE_PRECISION
  REAL ALPHA = (REAL)alpha;
#endif
  for (j = 0; j < n; ++j) {
    dst[j] = ALPHA*src[j];
  }
}

static void
product(opk_vspace_t* vspace, opk_vector_t* vdst,
        const opk_vector_t* vw, const opk_vector_t* vx)
{
  REAL* dst = DATA(vdst);
  const REAL* w = DATA(vw);
  const REAL* x = DATA(vx);
  opk_index_t j, n = vspace->size;
  for (j = 0; j < n; ++j) {
    dst[j] = w[j]*x[j];
  }
}

static void
axpby(opk_vspace_t* vspace, opk_vector_t* vdst,
      double alpha, const opk_vector_t* vx,
      double beta,  const opk_vector_t* vy)
{
  /* Note: we already know that neither ALPHA nor BETA is 0. */
  const REAL* x = DATA(vx);
  const REAL* y = DATA(vy);
  REAL* dst = DATA(vdst);
  opk_index_t j, n = vspace->size;
  if (alpha == 1.0) {
    if (beta == 1.0) {
      for (j = 0; j < n; ++j) {
        dst[j] = x[j] + y[j];
      }
    } else if (beta == -1.0) {
      for (j = 0; j < n; ++j) {
        dst[j] = x[j] - y[j];
      }
    } else {
#if SINGLE_PRECISION
      REAL BETA = (REAL)beta;
#endif
      for (j = 0; j < n; ++j) {
        dst[j] = x[j] + BETA*y[j];
      }
    }
  } else if (alpha == -1.0) {
    if (beta == 1.0) {
      for (j = 0; j < n; ++j) {
        dst[j] = y[j] - x[j];
      }
    } else if (beta == -1.0) {
      for (j = 0; j < n; ++j) {
        dst[j] = -y[j] - x[j];
      }
    } else {
#if SINGLE_PRECISION
      REAL BETA = (REAL)beta;
#endif
      for (j = 0; j < n; ++j) {
        dst[j] = BETA*y[j] - x[j];
      }
    }
  } else {
#if SINGLE_PRECISION
    REAL ALPHA = (REAL)alpha;
#endif
    if (beta == 1.0) {
      for (j = 0; j < n; ++j) {
        dst[j] = ALPHA*x[j] + y[j];
      }
    } else if (beta == -1.0) {
      for (j = 0; j < n; ++j) {
        dst[j] = ALPHA*x[j] - y[j];
      }
    } else {
#if SINGLE_PRECISION
      REAL BETA = (REAL)beta;
#endif
      for (j = 0; j < n; ++j) {
        dst[j] = ALPHA*x[j] + BETA*y[j];
      }
    }
  }
}

static void
axpbypcz(opk_vspace_t* vspace, opk_vector_t *vdst,
         double alpha, const opk_vector_t* vx,
         double beta,  const opk_vector_t* vy,
         double gamma, const opk_vector_t* vz)
{
  /* Note: we already know that neither ALPHA nor BETA nor GAMMA is 0. */
  const REAL* x = DATA(vx);
  const REAL* y = DATA(vy);
  const REAL* z = DATA(vz);
  REAL* dst = DATA(vdst);
  opk_index_t j, n = vspace->size;
#if SINGLE_PRECISION
  REAL ALPHA = (REAL)alpha;
  REAL BETA  = (REAL)beta;
  REAL GAMMA = (REAL)gamma;
#endif
  for (j = 0; j < n; ++j) {
    dst[j] = ALPHA*x[j] + BETA*y[j] + GAMMA*z[j];
  }
}
#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))

static int
boxprojvar(opk_vspace_t* space,
           opk_vector_t* dstvec,
           const opk_vector_t* srcvec,
           const void* lower,
           const void* upper,
           int bound)
{
  REAL* dst = DATA(dstvec);
  const REAL* x = DATA(srcvec);
  const REAL* xl;
  const REAL* xu;
  REAL a, b, t;
  opk_index_t i, n = space->size;

#define VALUE(addr) (*((double*)(addr)))
  switch (bound) {
  case 0:
    if (dst != x) {
      memcpy(dst, x, n*sizeof(REAL));
    }
    break;
  case 1:
    a = VALUE(lower);
    for (i = 0; i < n; ++i) {
      t = x[i];
      if (t < a) t = a;
      dst[i] = t;
    }
    break;
  case 2:
    xl = DATA(lower);
    for (i = 0; i < n; ++i) {
      a = xl[i];
      t = x[i];
      if (t < a) t = a;
      dst[i] = t;
    }
    break;
  case 3:
    b = VALUE(upper);
    for (i = 0; i < n; ++i) {
      t = x[i];
      if (t > b) t = b;
      dst[i] = t;
    }
    break;
  case 4:
    a = VALUE(lower);
    b = VALUE(upper);
    for (i = 0; i < n; ++i) {
      t = x[i];
      if (t < a) t = a;
      if (t > b) t = b;
      dst[i] = t;
    }
    break;
  case 5:
    xl = DATA(lower);
    b = VALUE(upper);
    for (i = 0; i < n; ++i) {
      a = xl[i];
      t = x[i];
      if (t < a) t = a;
      if (t > b) t = b;
      dst[i] = t;
    }
    break;
  case 6:
    xu = DATA(upper);
    for (i = 0; i < n; ++i) {
      t = x[i];
      b = xu[i];
      if (t > b) t = b;
      dst[i] = t;
    }
    break;
  case 7:
    a = VALUE(lower);
    xu = DATA(upper);
    for (i = 0; i < n; ++i) {
      t = x[i];
      b = xu[i];
      if (t < a) t = a;
      if (t > b) t = b;
      dst[i] = t;
    }
    break;
  case 8:
    xl = DATA(lower);
    xu = DATA(upper);
    for (i = 0; i < n; ++i) {
      a = xl[i];
      t = x[i];
      b = xu[i];
      if (t < a) t = a;
      if (t > b) t = b;
      dst[i] = t;
    }
    break;
  }
#undef VALUE
  return OPK_SUCCESS;
}

static int
boxprojdir(opk_vspace_t* space, opk_vector_t* dstvec,
           const opk_vector_t* srcvec,
           const void* lower, const void* upper, int bound,
           const opk_vector_t* dirvec, int orient)
{
  REAL* dst = DATA(dstvec);
  const REAL* x = DATA(srcvec);
  const REAL* d = DATA(dirvec);
  const REAL* xl;
  const REAL* xu;
  REAL a, b;
  opk_index_t i, n = space->size;

#define VALUE(addr)  (*((double*)(addr)))

#define BOXED(lo, hi)                                           \
  if (orient > 0) {                                             \
    for (i = 0; i < n; ++i) {                                   \
      dst[i] = (d[i] < 0 ? (x[i] > lo ? d[i] : 0) :             \
                (d[i] > 0 ? (x[i] < hi ? d[i] : 0) : 0));       \
    }                                                           \
  } else {                                                      \
    for (i = 0; i < n; ++i) {                                   \
      dst[i] = (d[i] > 0 ? (x[i] > lo ? d[i] : 0) :             \
                (d[i] < 0 ? (x[i] < hi ? d[i] : 0) : 0));       \
    }                                                           \
  }

#define LOWER(lo)                                               \
  if (orient > 0) {                                             \
    for (i = 0; i < n; ++i) {                                   \
      dst[i] = (d[i] >= 0 || x[i] > lo ? d[i] : 0);             \
    }                                                           \
  } else {                                                      \
    for (i = 0; i < n; ++i) {                                   \
      dst[i] = (d[i] <= 0 || x[i] > lo ? d[i] : 0);             \
    }                                                           \
  }

#define UPPER(hi)                                               \
  if (orient > 0) {                                             \
    for (i = 0; i < n; ++i) {                                   \
      dst[i] = (d[i] <= 0 || x[i] < hi ? d[i] : 0);             \
    }                                                           \
  } else {                                                      \
    for (i = 0; i < n; ++i) {                                   \
      dst[i] = (d[i] >= 0 || x[i] < hi ? d[i] : 0);             \
    }                                                           \
  }

  switch (bound) {
  case 0:
    if (dst != d) {
      memcpy(dst, d, n*sizeof(REAL));
    }
    break;
  case 1:
    a = VALUE(lower);
    LOWER(a);
    break;
  case 2:
    xl = DATA(lower);
    LOWER(xl[i]);
    break;
  case 3:
    b = VALUE(upper);
    UPPER(b);
    break;
  case 4:
    a = VALUE(lower);
    b = VALUE(upper);
    BOXED(a, b);
    break;
  case 5:
    xl = DATA(lower);
    b = VALUE(upper);
    BOXED(xl[i], b);
    break;
  case 6:
    xu = DATA(upper);
    UPPER(xu[i]);
    break;
  case 7:
    a = VALUE(lower);
    xu = DATA(upper);
    BOXED(a, xu[i]);
    break;
  case 8:
    xl = DATA(lower);
    xu = DATA(upper);
    BOXED(xl[i], xu[i]);
    break;
  }

#undef LOWER
#undef UPPER
#undef BOXED
#undef VALUE

  return OPK_SUCCESS;
}

static int
boxfreevar(opk_vspace_t* space, opk_vector_t* dstvec,
           const opk_vector_t* srcvec,
           const void* lower, const void* upper, int bound,
           const opk_vector_t* dirvec, int orient)
{
  REAL* dst = DATA(dstvec);
  const REAL* x = DATA(srcvec);
  const REAL* d = DATA(dirvec);
  const REAL* xl;
  const REAL* xu;
  REAL a, b;
  opk_index_t i, n = space->size;

#define VALUE(addr)    (*((double*)(addr)))

#define BOXED(lo, hi)                                   \
  if (orient > 0) {                                     \
    for (i = 0; i < n; ++i) {                           \
      dst[i] = d[i] < 0 ? (x[i] > lo ? 1 : 0) :         \
        (      d[i] > 0 ? (x[i] < hi ? 1 : 0) : 0);     \
    }                                                   \
  } else {                                              \
    for (i = 0; i < n; ++i) {                           \
      dst[i] = d[i] > 0 ? (x[i] > lo ? 1 : 0) :         \
        (      d[i] < 0 ? (x[i] < hi ? 1 : 0) : 0);     \
    }                                                   \
  }

#define LOWER(lo)                                       \
  if (orient > 0) {                                     \
    for (i = 0; i < n; ++i) {                           \
      dst[i] = d[i] < 0 ? (x[i] > lo ? 1 : 0) :         \
        (      d[i] > 0 ? 1 : 0);                       \
    }                                                   \
  } else {                                              \
    for (i = 0; i < n; ++i) {                           \
      dst[i] = d[i] > 0 ? (x[i] > lo ? 1 : 0) :         \
        (      d[i] < 0 ? 1 : 0);                       \
    }                                                   \
  }

#define UPPER(hi)                                       \
  if (orient > 0) {                                     \
    for (i = 0; i < n; ++i) {                           \
      dst[i] = d[i] < 0 ? 1 :                           \
        (      d[i] > 0 ? (x[i] < hi ? 1 : 0) : 0);     \
    }                                                   \
  } else {                                              \
    for (i = 0; i < n; ++i) {                           \
      dst[i] = d[i] > 0 ? 1 :                           \
        (      d[i] < 0 ? (x[i] < hi ? 1 : 0) : 0);     \
    }                                                   \
  }

  switch (bound) {
  case 0:
    for (i = 0; i < n; ++i) {
      dst[i] = 1;
    }
    break;
  case 1:
    a = VALUE(lower);
    LOWER(a);
    break;
  case 2:
    xl = DATA(lower);
    LOWER(xl[i]);
    break;
  case 3:
    b = VALUE(upper);
    UPPER(b);
    break;
  case 4:
    a = VALUE(lower);
    b = VALUE(upper);
    BOXED(a, b);
    break;
  case 5:
    xl = DATA(lower);
    b = VALUE(upper);
    BOXED(xl[i], b);
    break;
  case 6:
    xu = DATA(upper);
    UPPER(xu[i]);
    break;
  case 7:
    a = VALUE(lower);
    xu = DATA(upper);
    BOXED(a, xu[i]);
    break;
  case 8:
    xl = DATA(lower);
    xu = DATA(upper);
    BOXED(xl[i], xu[i]);
    break;
  }

#undef LOWER
#undef UPPER
#undef BOXED
#undef VALUE

  return OPK_SUCCESS;
}

static int
boxsteplimits(opk_vspace_t* space,
              double* smin, double* wolfe, double* smax,
              const opk_vector_t* xvec,
              const void* lower, const void* upper, int bound,
              const opk_vector_t* dvec, int orient)
{
  const REAL* x = DATA(xvec);
  const REAL* d = DATA(dvec);
  const REAL* xl;
  const REAL* xu;
  const REAL inf = FLOAT_CHOICE(FLT_MAX, DBL_MAX);
  REAL a, b;
  REAL s, s1 = inf, s2 = inf, s3 = 0;
  opk_index_t i, n = space->size;

#define VALUE(addr)    (*((double*)(addr)))

#define BOXED(lo, hi)                           \
  if (orient > 0) {                             \
    for (i = 0; i < n; ++i) {                   \
      REAL p = d[i];                            \
      if (p > 0) {                              \
        s = (hi - x[i])/p;                      \
      } else if (p < 0) {                       \
        s = (lo - x[i])/p;                      \
      } else {                                  \
        continue;                               \
      }                                         \
      if (s < s1) s1 = s;                       \
      if (s < s2 && s > 0) s2 = s;              \
      if (s > s3) s3 = s;                       \
    }                                           \
  } else {                                      \
    for (i = 0; i < n; ++i) {                   \
      REAL p = d[i];                            \
      if (p < 0) {                              \
        s = (x[i] - hi)/p;                      \
      } else if (p > 0) {                       \
        s = (x[i] - lo)/p;                      \
      } else {                                  \
        continue;                               \
      }                                         \
      if (s < s1) s1 = s;                       \
      if (s < s2 && s > 0) s2 = s;              \
      if (s > s3) s3 = s;                       \
    }                                           \
  }

#define LOWER(lo)                               \
  if (orient > 0) {                             \
    for (i = 0; i < n; ++i) {                   \
      REAL p = d[i];                            \
      if (p < 0) {                              \
        s = (lo - x[i])/p;                      \
        if (s < s1) s1 = s;                     \
        if (s < s2 && s > 0) s2 = s;            \
        if (s > s3) s3 = s;                     \
      } else if (p > 0) {                       \
        s3 = inf;                               \
      }                                         \
    }                                           \
  } else {                                      \
    for (i = 0; i < n; ++i) {                   \
      REAL p = d[i];                            \
      if (p > 0) {                              \
        s = (x[i] - lo)/p;                      \
        if (s < s1) s1 = s;                     \
        if (s < s2 && s > 0) s2 = s;            \
        if (s > s3) s3 = s;                     \
      } else if (p > 0) {                       \
        s3 = inf;                               \
      }                                         \
    }                                           \
  }

#define UPPER(hi)                               \
  if (orient > 0) {                             \
    for (i = 0; i < n; ++i) {                   \
      REAL p = d[i];                            \
      if (p > 0) {                              \
        s = (hi - x[i])/p;                      \
        if (s < s1) s1 = s;                     \
        if (s < s2 && s > 0) s2 = s;            \
        if (s > s3) s3 = s;                     \
      } else if (p < 0) {                       \
        s3 = inf;                               \
      }                                         \
    }                                           \
  } else {                                      \
    for (i = 0; i < n; ++i) {                   \
      REAL p = d[i];                            \
      if (p < 0) {                              \
        s = (x[i] - hi)/p;                      \
        if (s < s1) s1 = s;                     \
        if (s < s2 && s > 0) s2 = s;            \
        if (s > s3) s3 = s;                     \
      } else if (p > 0) {                       \
        s3 = inf;                               \
      }                                         \
    }                                           \
  }

  switch (bound) {
  case 0:
    s3 = inf;
    break;
  case 1:
    a = VALUE(lower);
    LOWER(a);
    break;
  case 2:
    xl = DATA(lower);
    LOWER(xl[i]);
    break;
  case 3:
    b = VALUE(upper);
    UPPER(b);
    break;
  case 4:
    a = VALUE(lower);
    b = VALUE(upper);
    BOXED(a, b);
    break;
  case 5:
    xl = DATA(lower);
    b = VALUE(upper);
    BOXED(xl[i], b);
    break;
  case 6:
    xu = DATA(upper);
    UPPER(xu[i]);
    break;
  case 7:
    a = VALUE(lower);
    xu = DATA(upper);
    BOXED(a, xu[i]);
    break;
  case 8:
    xl = DATA(lower);
    xu = DATA(upper);
    BOXED(xl[i], xu[i]);
    break;
  }

#undef LOWER
#undef UPPER
#undef BOXED
#undef VALUE

  if (smin  != NULL) *smin  = s1;
  if (wolfe != NULL) *wolfe = s2;
  if (smax  != NULL) *smax  = s3;
  return OPK_SUCCESS;
}

#define NEW_VECTOR_SPACE OPK_JOIN3(opk_new_simple_, REAL, _vector_space)
#define WRAP_VECTOR OPK_JOIN3(opk_wrap_simple_, REAL, _vector)
#define REWRAP_VECTOR OPK_JOIN3(opk_rewrap_simple_, REAL, _vector)
#define GET_DATA OPK_JOIN3(opk_get_simple_, REAL, _vector_data)
#define GET_CLIENT_DATA OPK_JOIN3(opk_get_simple_, REAL, _vector_client_data)
#define GET_FREE_CLIENT_DATA OPK_JOIN3(opk_get_simple_, REAL, _vector_free_client_data)

#if SINGLE_PRECISION
#  define NOUN "single"
#else
#  define NOUN "double"
#endif

static opk_vspace_operations_t operations = {
  "simple vector space for " NOUN " precision floating point values",
  NULL,
  create,
  finalize,
  fill,
  norm1,
  norm2,
  norminf,
  dot,
  dot3,
  copy,
  swap,
  scale,
  product,
  axpby,
  axpbypcz,
  boxprojvar,
  boxprojdir,
  boxfreevar,
  boxsteplimits
};

opk_vspace_t*
NEW_VECTOR_SPACE(opk_index_t size)
{
  return opk_allocate_vector_space(&operations, size, 0);
}

opk_vector_t*
WRAP_VECTOR(opk_vspace_t* vspace, REAL data[],
            void (*free_client_data)(void*), void* client_data)
{
  opk_vector_t* v;
  if (vspace->ops != &operations) {
    errno = EINVAL;
    return NULL;
  }
  if (data == NULL) {
    errno = EFAULT;
    return NULL;
  }
  v = opk_allocate_vector(vspace, sizeof(simple_vector_t));
  if (v != NULL) {
    simple_vector_t* sv = (simple_vector_t*)v;
    sv->data = data;
    sv->client_data = client_data;
    sv->free_client_data = free_client_data;
  }
  return v;
}

REAL*
GET_DATA(opk_vector_t* v)
{
  if (v == NULL) {
    errno = EFAULT;
    return NULL;
  }
  if (v->owner->ops != &operations) {
    errno = EINVAL;
    return NULL;
  }
  return ((simple_vector_t*)v)->data;
}

void*
GET_CLIENT_DATA(opk_vector_t* v)
{
  if (v == NULL) {
    errno = EFAULT;
    return NULL;
  }
  if (v->owner->ops != &operations) {
    errno = EINVAL;
    return NULL;
  }
  return ((simple_vector_t*)v)->client_data;
}

opk_free_proc*
GET_FREE_CLIENT_DATA(opk_vector_t* v)
{
  if (v == NULL) {
    errno = EFAULT;
    return NULL;
  }
  if (v->owner->ops != &operations) {
    errno = EINVAL;
    return NULL;
  }
  return ((simple_vector_t*)v)->free_client_data;
}

int
REWRAP_VECTOR(opk_vector_t* v, REAL new_data[],
              void (*new_free_client_data)(void*),
              void* new_client_data)
{
  simple_vector_t* sv;
  void *old_client_data;
  void (*old_free_client_data)(void*);

  /* Check arguments. */
  if (v == NULL) {
    errno = EFAULT;
    return OPK_FAILURE;
  }
  if (v->owner->ops != &operations) {
    errno = EINVAL;
    return OPK_FAILURE;
  }
  if (new_data == NULL) {
    errno = EFAULT;
    return OPK_FAILURE;
  }

  /* Get old members and make sure to not apply free_client_data again. */
  sv = (simple_vector_t*)v;
  old_client_data = sv->client_data;
  old_free_client_data = sv->free_client_data;
  sv->client_data = NULL;
  sv->free_client_data = NULL;
  if (old_free_client_data != NULL
      && (old_free_client_data != new_free_client_data
          || old_client_data != new_client_data)) {
    /* Apply old callback. */
    old_free_client_data(old_client_data);
  }

  /* Update contents. */
  sv->data = new_data;
  sv->client_data = new_client_data;
  sv->free_client_data = new_free_client_data;
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
