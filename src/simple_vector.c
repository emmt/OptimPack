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
#else
#  define REAL     double
#  define ABS(x)   fabs(x)
#  define SQRT(x)  sqrt(x)
#  define ALPHA    alpha
#  define BETA     beta
#  define GAMMA    gamma
#  define ZERO     0.0
#  define ONE      1.0
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
  axpbypcz
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
