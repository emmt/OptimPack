/*
 * algebra.c --
 *
 * High level vector space implementation for OptimPack library.
 *
 *-----------------------------------------------------------------------------
 *
 * The OptimPack library is licensed under the MIT "Expat" License:
 *
 * Copyright (c) 2014: Éric Thiébaut
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

double
opk_vdot3(const opk_vector_t* w, const opk_vector_t* x, const opk_vector_t* y)
{
  opk_vspace_t* vspace = w->owner;
  if (w->owner != vspace || x->owner != vspace || y->owner != vspace) {
    BAD_VECTORS("vdot3");
    return 0.0;
  } else {
    return vspace->ops->dot3(vspace, w, x, y);
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

/* Compute the elementwise product of two vectors. */
void
opk_vproduct(opk_vector_t* dst,
             const opk_vector_t* x, const opk_vector_t* y)
{
  opk_vspace_t* vspace = dst->owner;
  if (x->owner != vspace || y->owner != vspace) {
    BAD_VECTORS("vproduct");
  } else {
    vspace->ops->product(vspace, dst, x, y);
  }
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
/* BOUNDS AND BOXED SETS */

static void
finalize_bound(opk_object_t* obj)
{
  opk_bound_t* bnd = (opk_bound_t*)obj;
  opk_unset_bound(bnd);
  OPK_DROP(bnd->owner);
}

opk_bound_t*
opk_new_bound(opk_vspace_t* space, opk_bound_type_t type,
              void* value)
{
  opk_bound_t* bnd;

  /* Check arguments. */
  if (space == NULL || (type != OPK_BOUND_NONE && value == NULL)) {
    errno = EFAULT;
    return NULL;
  }
  if (space->size < 1 || type < 0 || type > 2 ||
      (type == OPK_BOUND_VECTOR && ((opk_vector_t*)value)->owner != space)) {
    errno = EINVAL;
    return NULL;
  }

  /* Allocate enough memory for the object and instanciate it. */
  bnd = (opk_bound_t*)opk_allocate_object(finalize_bound, sizeof(opk_bound_t));
  if (bnd == NULL) {
    return NULL;
  }
  bnd->owner = (opk_vspace_t*)OPK_HOLD(space);
  bnd->type = type;
  if (type == OPK_BOUND_VECTOR) {
    bnd->value.vector = (opk_vector_t*)OPK_HOLD(value);
  } else if (type == OPK_BOUND_SCALAR) {
    bnd->value.scalar = *(double*)value;
  }
  return bnd;
}

opk_bound_type_t
opk_get_bound_type(const opk_bound_t* bnd)
{
  return (bnd == NULL ? OPK_BOUND_NONE : bnd->type);
}

void
opk_unset_bound(opk_bound_t* bnd)
{
  opk_vector_t* vec;
  if (bnd != NULL) {
    vec = (bnd->type == OPK_BOUND_VECTOR ? bnd->value.vector : NULL);
    bnd->type = OPK_BOUND_NONE;
    bnd->value.vector = NULL;
    if (vec != NULL) {
      OPK_DROP(vec);
    }
  }
}

void
opk_set_scalar_bound(opk_bound_t* bnd, double val)
{
  if (bnd != NULL) {
    opk_unset_bound(bnd);
    bnd->value.scalar = val;
    bnd->type = OPK_BOUND_SCALAR;
  }
}

int
opk_set_vector_bound(opk_bound_t* bnd, opk_vector_t* vec)
{
  if (bnd != NULL) {
    if (vec != NULL && bnd->owner != vec->owner) {
      errno = EINVAL;
      return OPK_FAILURE;
    }
    opk_unset_bound(bnd);
    if (vec != NULL) {
      bnd->value.vector = (opk_vector_t*)OPK_HOLD(vec);
      bnd->type = OPK_BOUND_VECTOR;
    }
  }
  return OPK_SUCCESS;
}

double opk_box_shortcut_step(double alpha,
                             const opk_vector_t* x,
                             const opk_bound_t* xl,
                             const opk_bound_t* xu,
                             const opk_vector_t* d,
                             opk_orientation_t orient)
{
  double stpmax;
  if (opk_box_get_step_limits(NULL, NULL, &stpmax, x, xl, xu,
                              d, orient) != OPK_SUCCESS) {
    return 0;
  }
  return OPK_MIN(alpha, stpmax);
}

static int
get_bound(const void** lower, const opk_bound_t* xl,
          const void** upper, const opk_bound_t* xu)
{
  int type;

  if (xl != NULL && xl->type == OPK_BOUND_SCALAR) {
    type = 1;
    *lower = &(xl->value.scalar);
  } else if (xl != NULL && xl->type == OPK_BOUND_VECTOR) {
    type = 2;
    *lower = xl->value.vector;
  } else {
    type = 0;
    *lower = NULL;
  }
  if (xu != NULL && xu->type == OPK_BOUND_SCALAR) {
    type += 3;
    *upper = &(xu->value.scalar);
  } else if (xu != NULL && xu->type == OPK_BOUND_VECTOR) {
    type += 6;
    *upper = xu->value.vector;
  } else {
    *upper = NULL;
  }
  return type;
}

int
opk_box_project_variables(opk_vector_t* dst,
                          const opk_vector_t* x,
                          const opk_bound_t* xl,
                          const opk_bound_t* xu)
{
  opk_vspace_t* space;
  const void* lower;
  const void* upper;
  int type;

  if (dst == NULL || x == NULL) {
    errno = EFAULT;
    return OPK_FAILURE;
  }
  space = dst->owner;
  if (x->owner != space
      || (xl != NULL && xl->owner != space)
      || (xu != NULL && xu->owner != space)) {
    errno = EINVAL;
    return OPK_FAILURE;
  }
  type = get_bound(&lower, xl, &upper, xu);
  if (space->ops->boxprojvar == NULL) {
    errno = ENOSYS;
    return OPK_FAILURE;
  }
  return space->ops->boxprojvar(space, dst, x, lower, upper, type);
}

int
opk_box_project_direction(opk_vector_t* dst,
                          const opk_vector_t* x,
                          const opk_bound_t* xl,
                          const opk_bound_t* xu,
                          const opk_vector_t* d,
                          opk_orientation_t orient)
{
  opk_vspace_t* space;
  const void* lower;
  const void* upper;
  int type;

  if (dst == NULL || x == NULL || d == NULL) {
    errno = EFAULT;
    return OPK_FAILURE;
  }
  space = dst->owner;
  if (x->owner != space || d->owner != space
      || (xl != NULL && xl->owner != space)
      || (xu != NULL && xu->owner != space)) {
    errno = EINVAL;
    return OPK_FAILURE;
  }
  type = get_bound(&lower, xl, &upper, xu);
  if (space->ops->boxprojdir == NULL) {
    errno = ENOSYS;
    return OPK_FAILURE;
  }
  return space->ops->boxprojdir(space, dst, x, lower, upper, type,
                                d, orient);
}

int
opk_box_get_free_variables(opk_vector_t* dst,
                           const opk_vector_t* x,
                           const opk_bound_t* xl,
                           const opk_bound_t* xu,
                           const opk_vector_t* d,
                           opk_orientation_t orient)
{
  opk_vspace_t* space;
  const void* lower;
  const void* upper;
  int type;

  if (dst == NULL || x == NULL || d == NULL) {
    errno = EFAULT;
    return OPK_FAILURE;
  }
  space = dst->owner;
  if (x->owner != space || d->owner != space
      || (xl != NULL && xl->owner != space)
      || (xu != NULL && xu->owner != space)) {
    errno = EINVAL;
    return OPK_FAILURE;
  }
  type = get_bound(&lower, xl, &upper, xu);
  if (space->ops->boxfreevar == NULL) {
    errno = ENOSYS;
    return OPK_FAILURE;
  }
  return space->ops->boxfreevar(space, dst, x, lower, upper, type,
                                d, orient);
}

extern int
opk_box_get_step_limits(double* smin, double* wolfe, double *smax,
                        const opk_vector_t* x,
                        const opk_bound_t* xl,
                        const opk_bound_t* xu,
                        const opk_vector_t* d,
                        opk_orientation_t orient)
{
  opk_vspace_t* space;
  const void* lower;
  const void* upper;
  int type;

  if (x == NULL || d == NULL) {
    errno = EFAULT;
    return OPK_FAILURE;
  }
  space = x->owner;
  if (d->owner != space
      || (xl != NULL && xl->owner != space)
      || (xu != NULL && xu->owner != space)) {
    errno = EINVAL;
    return OPK_FAILURE;
  }
  type = get_bound(&lower, xl, &upper, xu);
  if (space->ops->boxsteplimits == NULL) {
    errno = ENOSYS;
    return OPK_FAILURE;
  }
  return space->ops->boxsteplimits(space, smin, wolfe, smax,
                                   x, lower, upper, type, d, orient);
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
