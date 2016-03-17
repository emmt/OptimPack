/*
 * algebra.c --
 *
 * High level vector space implementation for OptimPack library.
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of OptimPack (https://github.com/emmt/OptimPack).
 *
 * Copyright (C) 2014-2016 Éric Thiébaut
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

opk_status_t
opk_vpeek(const opk_vector_t* vect, opk_index_t k, double* ptr)
{
  if (vect == NULL || ptr == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  if (k < 0 || k >= vect->owner->size) {
    return OPK_OUT_OF_BOUNDS_INDEX;
  }
  *ptr = vect->owner->ops->peek(vect->owner, vect, k);
  return OPK_SUCCESS;
}

opk_status_t
opk_vpoke(opk_vector_t* vect, opk_index_t k, double value)
{
  if (vect == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  if (k < 0 || k >= vect->owner->size) {
    return OPK_OUT_OF_BOUNDS_INDEX;
  }
  vect->owner->ops->poke(vect->owner, vect, k, value);
  return OPK_SUCCESS;
}

void
opk_vprint(FILE* file, const char* name, const opk_vector_t* vect,
           opk_index_t nmax)
{
  opk_vspace_t* vspace;
  opk_index_t i, n;

  if (vect == NULL) {
    if (name != NULL) {
      fputs(name, file);
      fputs(" = NULL;\n", file);
    } else {
      fputs("NULL;\n", file);
    }
    return;
  }
  if (file == NULL) {
    file = stdout;
  }
  if (name != NULL) {
    fputs(name, file);
    fputs(" = {", file);
  } else {
    fputs("{", file);
  }
  vspace = vect->owner;
  n = vspace->size;
  if (nmax > 0 && nmax < n) {
    n = nmax;
  }
  for (i = 0; i < n; ++i) {
    fprintf(file, "%g", vspace->ops->peek(vspace, vect, i));
    fputs((i < n - 1 ? "," : ""), file);
  }
  fputs((n < vspace->size ? ",...};\n" : "};\n"), file);
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
/* CONVEX SETS */

static void
finalize_convexset(opk_object_t* obj)
{
  opk_convexset_t* set = (opk_convexset_t*)obj;
  if (set->finalize != NULL) {
    set->finalize(set);
  }
  OPK_DROP(set->space);
}

opk_convexset_t*
opk_allocate_convexset(opk_vspace_t* space,
                       void (*finalize)(opk_convexset_t* self),
                       opk_status_t (*projvar)(opk_vector_t* dst,
                                               const opk_vector_t* x,
                                               const opk_convexset_t* set),
                       opk_status_t (*projdir)(opk_vector_t* dst,
                                               const opk_vector_t* x,
                                               const opk_convexset_t* set,
                                               const opk_vector_t* d,
                                               opk_orientation_t orient),
                       opk_status_t (*freevar)(opk_vector_t* dst,
                                               const opk_vector_t* x,
                                               const opk_convexset_t* set,
                                               const opk_vector_t* d,
                                               opk_orientation_t orient),
                       opk_status_t (*steplim)(double* smin1, double* smin2,
                                               double *smax,
                                               const opk_vector_t* x,
                                               const opk_convexset_t* set,
                                               const opk_vector_t* d,
                                               opk_orientation_t orient),
                       size_t nbytes)
{
  opk_convexset_t* set;

  if (space == NULL) {
    errno = EFAULT;
    return NULL;
  }
  if (nbytes < sizeof(opk_convexset_t)) {
    nbytes = sizeof(opk_convexset_t);
  }
  set = (opk_convexset_t*)opk_allocate_object(finalize_convexset, nbytes);
  if (set != NULL) {
    set->space = OPK_HOLD_VSPACE(space);
    set->finalize = finalize;
    set->projvar = projvar;
    set->projdir = projdir;
    set->freevar = freevar;
    set->steplim = steplim;
  }
  return set;
}

opk_bool_t
opk_can_project_variables(const opk_convexset_t* set)
{
  return (set != NULL && set->projvar != NULL);
}

opk_bool_t
opk_can_project_directions(const opk_convexset_t* set)
{
  return (set != NULL && set->projdir != NULL);
}

opk_bool_t
opk_can_get_free_variables(const opk_convexset_t* set)
{
  return (set != NULL && set->freevar != NULL);
}

opk_bool_t
opk_can_get_step_limits(const opk_convexset_t* set)
{
  return (set != NULL && set->steplim != NULL);
}

opk_status_t
opk_project_variables(opk_vector_t* dst,
                      const opk_vector_t* x,
                      const opk_convexset_t* set)
{
  opk_vspace_t* space;

  if (dst == NULL || x == NULL || set == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  space = set->space;
  if (dst->owner != space || x->owner != space) {
    return OPK_BAD_SPACE;
  }
  if (set->projvar == NULL) {
    return OPK_NOT_IMPLEMENTED;
  }
  return set->projvar(dst, x, set);
}

opk_status_t
opk_project_direction(opk_vector_t* dst,
                      const opk_vector_t* x,
                      const opk_convexset_t* set,
                      const opk_vector_t* d,
                      opk_orientation_t orient)
{
  opk_vspace_t* space;

  if (dst == NULL || x == NULL || set == NULL || d == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  space = set->space;
  if (dst->owner != space || x->owner != space || d->owner != space) {
    return OPK_BAD_SPACE;
  }
  if (set->projdir == NULL) {
    return OPK_NOT_IMPLEMENTED;
  }
  return set->projdir(dst, x, set, d, orient);
}

opk_status_t
opk_get_free_variables(opk_vector_t* dst,
                       const opk_vector_t* x,
                       const opk_convexset_t* set,
                       const opk_vector_t* d,
                       opk_orientation_t orient)
{
  opk_vspace_t* space;

  if (dst == NULL || x == NULL || set == NULL || d == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  space = set->space;
  if (dst->owner != space || x->owner != space || d->owner != space) {
    return OPK_BAD_SPACE;
  }
  if (set->freevar == NULL) {
    return OPK_NOT_IMPLEMENTED;
  }
  return set->freevar(dst, x, set, d, orient);
}

opk_status_t
opk_get_step_limits(double* smin1, double* smin2, double *smax,
                    const opk_vector_t* x,
                    const opk_convexset_t* set,
                    const opk_vector_t* d,
                    opk_orientation_t orient)
{
  opk_vspace_t* space;

  if (x == NULL || set == NULL || d == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  space = set->space;
  if (x->owner != space || d->owner != space) {
    return OPK_BAD_SPACE;
  }
  if (set->steplim == NULL) {
    return OPK_NOT_IMPLEMENTED;
  }
  return set->steplim(smin1, smin2, smax, x, set, d, orient);
}

/* A box set is a convex set with simple separable bounds.  This type is not
   publically available (the caller receives an `opk_convexset_t` type
   instead. */

typedef struct _opk_bound {
  union {
    opk_vector_t* vector; /**< The bound as a vector. */
    double        scalar; /**< The bound as a scalar. */
  } value;                /**< Value of the bound. */
  opk_bound_type_t type;  /**< Type of the bound. */
} opk_bound_t;

typedef struct _opk_boxset {
  opk_convexset_t base;  /**< Base type (must be the first member). */
  opk_bound_t     lower; /**< Lower bound. */
  opk_bound_t     upper; /**< Upper bound. */
} opk_boxset_t;

static void
finalize_boxset(opk_convexset_t* set)
{
  opk_boxset_t* box = (opk_boxset_t*)set;
  if (box->lower.type == OPK_BOUND_VECTOR) {
    OPK_DROP(box->lower.value.vector);
  }
  if (box->upper.type == OPK_BOUND_VECTOR) {
    OPK_DROP(box->upper.value.vector);
  }
}

static int
get_bounds(const void** lower, const void** upper, const opk_boxset_t* box)
{
  int type;

  if (box->lower.type == OPK_BOUND_SCALAR) {
    type = 1;
    *lower = &box->lower.value.scalar;
  } else if (box->lower.type == OPK_BOUND_VECTOR) {
    type = 2;
    *lower = box->lower.value.vector;
  } else {
    type = 0;
    *lower = NULL;
  }
  if (box->upper.type == OPK_BOUND_SCALAR) {
    type += 3;
    *upper = &box->upper.value.scalar;
  } else if (box->upper.type == OPK_BOUND_VECTOR) {
    type += 6;
    *upper = box->upper.value.vector;
  } else {
    *upper = NULL;
  }
  return type;
}

static opk_status_t
box_projvar(opk_vector_t* dst,
            const opk_vector_t* x,
            const opk_convexset_t* set)
{
  const void* lower;
  const void* upper;
  opk_vspace_t* space = set->space;
  int type = get_bounds(&lower, &upper, (opk_boxset_t*)set);
  return space->ops->boxprojvar(space, dst, x, lower, upper, type);
}

static opk_status_t
box_projdir(opk_vector_t* dst,
            const opk_vector_t* x,
            const opk_convexset_t* set,
            const opk_vector_t* d,
            opk_orientation_t orient)
{
  const void* lower;
  const void* upper;
  opk_vspace_t* space = set->space;
  int type = get_bounds(&lower, &upper, (opk_boxset_t*)set);
  return space->ops->boxprojdir(space, dst,
                                x, lower, upper, type, d, orient);
}

static opk_status_t
box_freevar(opk_vector_t* dst,
            const opk_vector_t* x,
            const opk_convexset_t* set,
            const opk_vector_t* d,
            opk_orientation_t orient)
{
  const void* lower;
  const void* upper;
  opk_vspace_t* space = set->space;
  int type = get_bounds(&lower, &upper, (opk_boxset_t*)set);
  return space->ops->boxfreevar(set->space, dst,
                                x, lower, upper, type, d, orient);
}

static opk_status_t
box_steplim(double* smin1, double* smin2, double *smax,
            const opk_vector_t* x,
            const opk_convexset_t* set,
            const opk_vector_t* d,
            opk_orientation_t orient)
{
  const void* lower;
  const void* upper;
  opk_vspace_t* space = set->space;
  int type = get_bounds(&lower, &upper, (opk_boxset_t*)set);
  return space->ops->boxsteplim(space, smin1, smin2, smax,
                                x, lower, upper, type, d, orient);
}

static opk_status_t
set_bound(const opk_vspace_t* space, opk_bound_t* bound,
          opk_bound_type_t type, void* value)
{
  memset(bound, 0, sizeof(opk_bound_t));
  if (type == OPK_BOUND_NONE) {
    if (value != NULL) {
      return OPK_INVALID_ARGUMENT;
    }
  } else if (type == OPK_BOUND_SCALAR) {
    if (value == NULL) {
      return OPK_ILLEGAL_ADDRESS;
    }
    bound->value.scalar = *(double*)value;
  } else if (type == OPK_BOUND_VECTOR) {
    if (value == NULL) {
      return OPK_ILLEGAL_ADDRESS;
    }
    if (((opk_vector_t*)value)->owner != space) {
      return OPK_BAD_SPACE;
    }
    bound->value.vector = (opk_vector_t*)OPK_HOLD(value);
  } else {
    return OPK_INVALID_ARGUMENT;
  }
  bound->type = type;
  return OPK_SUCCESS;
}

opk_convexset_t*
opk_new_boxset(opk_vspace_t* space,
               opk_bound_type_t lower_type, void* lower,
               opk_bound_type_t upper_type, void* upper)
{
  opk_boxset_t* box;
  const opk_vspace_operations_t* ops;
  opk_status_t (*projvar)(opk_vector_t* dst,
                          const opk_vector_t* x,
                          const opk_convexset_t* set);
  opk_status_t (*projdir)(opk_vector_t* dst,
                          const opk_vector_t* x,
                          const opk_convexset_t* set,
                          const opk_vector_t* d,
                          opk_orientation_t orient);
  opk_status_t (*freevar)(opk_vector_t* dst,
                          const opk_vector_t* x,
                          const opk_convexset_t* set,
                          const opk_vector_t* d,
                          opk_orientation_t orient);
  opk_status_t (*steplim)(double* smin1, double* smin2,
                          double *smax,
                          const opk_vector_t* x,
                          const opk_convexset_t* set,
                          const opk_vector_t* d,
                          opk_orientation_t orient);

  if (space == NULL) {
    errno = EFAULT;
    return NULL;
  }
  ops = space->ops;
  projvar = (ops->boxprojvar != NULL ? box_projvar : NULL);
  projdir = (ops->boxprojdir != NULL ? box_projdir : NULL);
  freevar = (ops->boxfreevar != NULL ? box_freevar : NULL);
  steplim = (ops->boxsteplim != NULL ? box_steplim : NULL);
  box = (opk_boxset_t*)opk_allocate_convexset(space, finalize_boxset,
                                              projvar, projdir,
                                              freevar, steplim,
                                              sizeof(opk_boxset_t));
  if (box != NULL) {
    if (set_bound(space, &box->lower, lower_type, lower) != OPK_SUCCESS ||
        set_bound(space, &box->upper, upper_type, upper) != OPK_SUCCESS) {
      OPK_DROP(box);
      errno = EINVAL;
      return NULL;
    }
  }
  return (opk_convexset_t*)box;
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

opk_status_t
opk_apply_direct(opk_operator_t* op, opk_vector_t* dst,
                 const opk_vector_t* src)
{
  if (op == NULL || dst == NULL || src == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  if (dst->owner != op->outspace || src->owner != op->inpspace) {
    return OPK_BAD_SPACE;
  }
  if (op->ops->apply_direct == NULL) {
    return OPK_NOT_IMPLEMENTED;
  }
  return op->ops->apply_direct(op, dst, src);
}

opk_status_t
opk_apply_adjoint(opk_operator_t* op, opk_vector_t* dst,
                 const opk_vector_t* src)
{
  if (op == NULL || dst == NULL || src == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  if (dst->owner != op->inpspace || src->owner != op->outspace) {
    return OPK_BAD_SPACE;
  }
  if (op->ops->apply_adjoint == NULL) {
    return OPK_NOT_IMPLEMENTED;
  }
  return op->ops->apply_adjoint(op, dst, src);
}

opk_status_t
opk_apply_inverse(opk_operator_t* op, opk_vector_t* dst,
                 const opk_vector_t* src)
{
  if (op == NULL || dst == NULL || src == NULL) {
    return OPK_ILLEGAL_ADDRESS;
  }
  if (dst->owner != op->inpspace || src->owner != op->outspace) {
    return OPK_BAD_SPACE;
  }
  if (op->ops->apply_inverse == NULL) {
    return OPK_NOT_IMPLEMENTED;
  }
  return op->ops->apply_inverse(op, dst, src);
}
