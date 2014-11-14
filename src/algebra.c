/*
 * algebra.c --
 *
 * High level vector space implementation for OptimPack library.
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

/*---------------------------------------------------------------------------*/
/* BASIC OBJECTS */

opk_object_t*
opk_allocate_object(void (*finalize)(opk_object_t* self),
                    size_t nbytes)
{
  opk_object_t* object;
  if (nbytes < sizeof(opk_object_t)) {
    nbytes = sizeof(opk_object_t);
  }
  object = (opk_object_t*)malloc(nbytes);
  if (object != NULL) {
    memset(object, 0, nbytes);
    object->finalize = finalize;
    object->references = 1;
  }
  return object;
}

opk_object_t*
opk_hold_object(opk_object_t* obj)
{
  ++obj->references;
  return obj;
}

void
opk_drop_object(opk_object_t* obj)
{
  if (obj != NULL && --obj->references < 1) {
    if (obj->finalize != NULL) {
      obj->finalize(obj);
    }
    free((void*)obj);
  }
}

opk_index_t
opk_get_object_references(opk_object_t* obj)
{
  return (obj != NULL ? obj->references : 0);
}

/*---------------------------------------------------------------------------*/
/* VECTOR SPACES */

static void
finalize_vector_space(opk_object_t* obj)
{
  opk_vspace_t* vspace = (opk_vspace_t*)obj;
  if (vspace->ops->finalize_space != NULL) {
    vspace->ops->finalize_space(vspace);
  }
}

opk_vspace_t*
opk_allocate_vector_space(const opk_vspace_operations_t* ops,
                          opk_index_t nvariables,
                          size_t nbytes)
{
  opk_vspace_t* vspace;
  if (nvariables < 1) {
    errno = EINVAL;
    return NULL;
  }
  if (nbytes < sizeof(opk_vspace_t)) {
    nbytes = sizeof(opk_vspace_t);
  }
  vspace = (opk_vspace_t*)opk_allocate_object(finalize_vector_space, nbytes);
  if (vspace != NULL) {
    vspace->ops = ops;
    vspace->size = nvariables;
  }
  return vspace;
}

/*---------------------------------------------------------------------------*/
/* VECTORS */

/* Destroy a vector when no longer referenced. */
static void
finalize_vector(opk_object_t* obj)
{
  opk_vector_t* vector = (opk_vector_t*)obj;
  opk_vspace_t* vspace = vector->owner;
  vspace->ops->finalize_vector(vspace, vector);
  OPK_DROP(vspace);
}

/* Utility for allocating a new vector object. */
opk_vector_t*
opk_allocate_vector(opk_vspace_t* vspace, size_t nbytes)
{
  opk_vector_t* vect;

  if (nbytes < sizeof(opk_vector_t)) {
    nbytes = sizeof(opk_vector_t);
  }
  vect = (opk_vector_t*)opk_allocate_object(finalize_vector, nbytes);
  if (vect != NULL) {
    vect->owner = OPK_HOLD_VSPACE(vspace);
  }
  return vect;
}

#define BAD_VECTOR(name)                                                \
  opk_error("vector does not belong to the correct space in `opk_"      \
            name "`")
#define BAD_VECTORS(name)                                               \
  opk_error("vectors do not belong to the same space in `opk_"          \
            name "`")

/* Create a new vector. */
opk_vector_t*
opk_vcreate(opk_vspace_t* vspace)
{
  return vspace->ops->create(vspace);
}

/* Zero-fill a vector. */
void
opk_vzero(opk_vector_t* vect)
{
  opk_vspace_t* vspace = vect->owner;
  vspace->ops->fill(vspace, vect, 0.0);
}

/* Fill a vector with a scalar. */
void
opk_vfill(opk_vector_t* vect, double alpha)
{
  opk_vspace_t* vspace = vect->owner;
  vspace->ops->fill(vspace, vect, alpha);
}

/* Copy vector contents. */
void
opk_vcopy(opk_vector_t* dst, const opk_vector_t* src)
{
  if (src != dst) {
    opk_vspace_t* vspace = dst->owner;
    if (src->owner != vspace) {
      BAD_VECTORS("vcopy");
    } else {
      vspace->ops->copy(vspace, dst, src);
    }
  }
}

/* Scale vector contents. */
void
opk_vscale(opk_vector_t* dst, double alpha, const opk_vector_t* src)
{
  opk_vspace_t* vspace = dst->owner;
  if (src->owner != vspace) {
    BAD_VECTORS("vscale");
  } else {
    if (alpha == 1.0) {
      if (src != dst) {
        vspace->ops->copy(vspace, dst, src);
      }
    } else if (alpha == 0.0) {
      vspace->ops->fill(vspace, dst, 0.0);
    } else {
      vspace->ops->scale(vspace, dst, alpha, src);
    }
  }
}

/* Exchange vector contents. */
void
opk_vswap(opk_vector_t* x, opk_vector_t* y)
{
  if (x != y) {
    opk_vspace_t* vspace = x->owner;
    if (y->owner != vspace) {
      BAD_VECTORS("vswap");
    } else {
      vspace->ops->swap(vspace, x, y);
    }
  }
}

/* Compute inner product. */
double
opk_vdot(const opk_vector_t* x, const opk_vector_t* y)
{
  opk_vspace_t* vspace = x->owner;
  if (y->owner != vspace) {
    BAD_VECTORS("vdot");
    return 0.0;
  } else {
    return vspace->ops->dot(vspace, x, y);
  }
}

/* Compute L2 norm. */
double
opk_vnorm2(const opk_vector_t* x)
{
  opk_vspace_t* vspace = x->owner;
  return vspace->ops->norm2(vspace, x);
}

/* Compute L1 norm. */
double
opk_vnorm1(const opk_vector_t* x)
{
  opk_vspace_t* vspace = x->owner;
  return vspace->ops->norm1(vspace, x);
}

/* Compute infinite norm. */
double
opk_vnorminf(const opk_vector_t* x)
{
  opk_vspace_t* vspace = x->owner;
  return vspace->ops->norminf(vspace, x);
}

/* Compute linear combination of two vectors. */
void
opk_vaxpby(opk_vector_t* dst,
           double alpha, const opk_vector_t* x,
           double beta,  const opk_vector_t* y)
{
  opk_vspace_t* vspace = dst->owner;
  if (x->owner != vspace || y->owner != vspace) {
    BAD_VECTORS("vaxpby");
  } else {
    if (alpha == 0.0) {
      if (beta == 0.0) {
        vspace->ops->fill(vspace, dst, 0.0);
      } else {
        vspace->ops->scale(vspace, dst, beta, y);
      }
    } else {
      if (beta == 0.0) {
        vspace->ops->scale(vspace, dst, alpha, x);
      } else {
        vspace->ops->axpby(vspace, dst, alpha, x, beta, y);
      }
    }
  }
}

/* Compute linear combination of three vectors. */
void
opk_vaxpbypcz(opk_vector_t* dst,
              double alpha, const opk_vector_t* x,
              double beta,  const opk_vector_t* y,
              double gamma, const opk_vector_t* z)
{
  opk_vspace_t* vspace = dst->owner;
  if (x->owner != vspace || y->owner != vspace || z->owner != vspace) {
    BAD_VECTORS("vaxpbypcz");
  } else {
    if (alpha == 0.0) {
      if (beta == 0.0) {
        if (gamma == 0.0) {
          vspace->ops->fill(vspace, dst, 0.0);
        } else {
          vspace->ops->scale(vspace, dst, gamma, z);
        }
      } else {
        if (gamma == 0.0) {
          vspace->ops->scale(vspace, dst, beta, y);
        } else {
          vspace->ops->axpby(vspace, dst, beta, y, gamma, z);
        }
      }
    } else {
      if (beta == 0.0) {
        if (gamma == 0.0) {
          vspace->ops->scale(vspace, dst, alpha, x);
        } else {
          vspace->ops->axpby(vspace, dst, alpha, x, gamma, z);
        }
      } else {
        if (gamma == 0.0) {
          vspace->ops->axpby(vspace, dst, alpha, x, beta, y);
        } else {
          vspace->ops->axpbypcz(vspace, dst, alpha, x, beta, y, gamma, z);
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/* OPERATORS */

static void
finalize_operator(opk_object_t* obj)
{
  opk_operator_t* op = (opk_operator_t*)obj;
  if (op->ops->finalize != NULL) {
    op->ops->finalize(op);
  }
  OPK_DROP(op->inpspace);
  OPK_DROP(op->outspace);
}

opk_operator_t*
opk_allocate_operator(const opk_operator_operations_t* ops,
                      opk_vspace_t* inpspace,
                      opk_vspace_t* outspace,
                      size_t size)
{
  opk_operator_t* op;

  /* Check arguments. */
  if (inpspace == NULL || outspace == NULL) {
    errno = EFAULT;
    return NULL;
  }
  if (size < sizeof(opk_operator_t)) {
    size = sizeof(opk_operator_t);
  }
  op = (opk_operator_t*)opk_allocate_object(finalize_operator, size);
  if (op != NULL) {
    op->ops = ops;
    op->inpspace = OPK_HOLD_VSPACE(inpspace);
    op->outspace = OPK_HOLD_VSPACE(outspace);
  }
  return op;
}

int
opk_apply_direct(opk_operator_t* op, opk_vector_t* dst,
                 const opk_vector_t* src)
{
  if (op == NULL || dst == NULL || src == NULL) {
    errno = EFAULT;
    return OPK_FAILURE;
  }
  if (dst->owner != op->outspace || src->owner != op->inpspace) {
    errno = EINVAL;
    return OPK_FAILURE;
  }
  if (op->ops->apply_direct == NULL) {
    errno = EPERM; /* Operation not permitted */
    return OPK_FAILURE;
  }
  return op->ops->apply_direct(op, dst, src);
}

int
opk_apply_adjoint(opk_operator_t* op, opk_vector_t* dst,
                 const opk_vector_t* src)
{
  if (op == NULL || dst == NULL || src == NULL) {
    errno = EFAULT;
    return OPK_FAILURE;
  }
  if (dst->owner != op->inpspace || src->owner != op->outspace) {
    errno = EINVAL;
    return OPK_FAILURE;
  }
  if (op->ops->apply_adjoint == NULL) {
    errno = EPERM; /* Operation not permitted */
    return OPK_FAILURE;
  }
  return op->ops->apply_adjoint(op, dst, src);
}

int
opk_apply_inverse(opk_operator_t* op, opk_vector_t* dst,
                 const opk_vector_t* src)
{
  if (op == NULL || dst == NULL || src == NULL) {
    errno = EFAULT;
    return OPK_FAILURE;
  }
  if (dst->owner != op->inpspace || src->owner != op->outspace) {
    errno = EINVAL;
    return OPK_FAILURE;
  }
  if (op->ops->apply_inverse == NULL) {
    errno = EPERM; /* Operation not permitted */
    return OPK_FAILURE;
  }
  return op->ops->apply_inverse(op, dst, src);
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
