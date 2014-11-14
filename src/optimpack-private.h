/*
 * optimpack-private.h --
 *
 * Private definitions for building OptimPack library.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2014 Éric Thiébaut <thiebaut@obs.univ-lyon1.fr>
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

#ifndef _OPTIMPACK_PRIVATE_H
#define _OPTIMPACK_PRIVATE_H 1

#include "optimpack.h"


/*---------------------------------------------------------------------------*/
/* LOW LEVEL API FOR BASIC OBJECTS AND DERIVED TYPES */

/**
 * Basic object structure.
 *
 * The definition of this structure must be exposed so that others types
 * can "inherit" from it.  However, the contents must be left unchanged
 * otherwise unexpected results may occur.
 */
struct _opk_object {
  void (*finalize)(opk_object_t* self);
  long references;
};

/**
 * Allocate a new (empty) object.
 *
 * This function allocates enough bytes of memory to store a type derived
 * from the basic object type, initializes the basic type and zero-fill
 * the extra bytes of allocated memory.
 *
 * The caller of this function holds a reference on the returned object.  When
 * the object is no longer needed by the caller, he/she has to call
 * opk_drop_object() to release this reference.
 *
 * When the object is no longer in use (its last reference has been dropped),
 * it is effectively destroyed by first applying the finalize() method on it
 * (unless it is NULL) and then freeing the dynamic memory where the object
 * was stored.
 *
 * The typical usage consists in building a derived type as follows:
 * <pre>
 * typedef struct _sub_type {
 *     opk_object_t base; // base type (must be the first member)
 *     ...;               // other members
 * } sub_type_t;
 * </pre>
 * and provide a constructor like:
 * <pre>
 * sub_type_t* create_sub_type(args...) {
 *    sub_type_t* obj = (sub_type_t*)opk_allocate_object(finalize_sub_type,
 *                                                       sizeof(sub_type_t));
 *    ...; // any other initialization
 *    return obj;
 * }
 * </pre>
 * where finalize_sub_type() is in charge of releasing any specific ressources
 * of the object.  Note that this function must be able to deal with a
 * partially initialized object, although this is simplified by the fact that
 * the non-basic parts of the object are intially zero-filled.  The memory
 * allocated by opk_allocate_object() is automatically freed and must not be
 * freed by the finalize() method.
 *
 * The object model of OptimPack is very simple but has some drawbacks (C is
 * not C++): casts are needed to force the proper type of an object pointer,
 * memory management by reference counting forbid circular references.
 *
 * @param finalize - If not NULL, method to call when the object is no longer
 *                   referenced.  This method is in charge of freeing object
 *                   ressources specific to the sub-type.  This method is
 *                   *never* called with a NULL argument.
 *
 * @param nbytes   - Number of bytes to allocate for the object.  This value
 *                   may be adjusted so that at least a basic object can be
 *                   stored.
 *
 * @return A new basic object with one reference set or NULL in case of error.
 */
extern opk_object_t*
opk_allocate_object(void (*finalize)(opk_object_t* self),
                    size_t nbytes);

/*---------------------------------------------------------------------------*/
/* LOW LEVEL API FOR VECTOR SPACES AND DERIVED TYPES */

typedef struct _opk_vspace_operations opk_vspace_operations_t;

/* The _opk_vector structure is intentionally exposed to let different
 * implementations coexist (although in separate codes).   If one want
 * to implement a sub-type of the vector type, it is sufficient to define
 * a new structure whose first member is an `opk_vector_t`.  The function
 * `opk_valloc()` can be used to allocate the whole structure.
 * <pre>
 * struct _my_vector {
 *   opk_vector_t base;
 *   double* data;
 * };
 * </pre>
 *
 * OptimPack routines
 * only require the address of such vectors and treat them as opaque
 * structures. */
struct _opk_vector {
  opk_object_t  base;   /* base type (must be the first member) */
  opk_vspace_t* owner;  /* vector space to which the vector belongs */
};

struct _opk_vspace {
  opk_object_t base;    /* base type (must be the first member) */

  /* Table of operations on vectors of this space. */
  const opk_vspace_operations_t* ops;

  /* Number of elements of vectors of this vector space (real problems are in
     finite dimension!). */
  opk_index_t size;
};

struct _opk_vspace_operations {
  /* A short description of the vector space. */
  const char* description;

  /* Method called to effectively free anything else than the memory allocated
     to store the vector space itself.  It is guaranteed that this method is
     called with a non-NULL argument and that no references exist on the vector
     space (in particular all vectors of this vector space have been
     finalized). */
  void (*finalize_space)(opk_vspace_t* vspace);

  /* Create an element of this vector space (the contents of the vector is
     undefined). */
  opk_vector_t* (*create)(opk_vspace_t* vspace);

  /* Delete an element of this vector space.  Same semantics as the
     finalize_space() method. FIXME: improve explanations */
  void (*finalize_vector)(opk_vspace_t* vspace, opk_vector_t* v);

  /* Fill a vector with a scalar ALPHA.  When this method is called, it is
     guaranteed that X belong to VSPACE.  The implementation of this method
     may be optimized when ALPHA is equal to zero. */
  void (*fill)(opk_vspace_t* vspace, opk_vector_t* v, double alpha);

  /* Compute the L1 norm (sum of absolute values) of an element of this vector
     space.  When this method is called, it is guaranteed that X belongs to
     VSPACE. */
  double (*norm1)(opk_vspace_t* vspace, const opk_vector_t* x);

  /* Compute the L2 norm (square root of the sum of squared values, also known
     as the Euclidean norm) of an element of this vector space.  When this
     method is called, it is guaranteed that X belongs to VSPACE. */
  double (*norm2)(opk_vspace_t* vspace, const opk_vector_t* x);

  /* Compute the infinite norm (maximum absolute value) of an element of this
     vector space.  When this method is called, it is guaranteed that X
     belongs to VSPACE. */
  double (*norminf)(opk_vspace_t* vspace, const opk_vector_t* x);

  /* Compute the inner product between two elements of this vector space.
     When this method is called, it is guaranteed that X and Y belong to
     VSPACE. */
  double (*dot)(opk_vspace_t* vspace,
                const opk_vector_t* x,
                const opk_vector_t* y);

  /* Copy vector SRC to vector DST.  When this method is called, it is
     guaranteed that DST and SRC are different and both belong to VSPACE. */
  void (*copy)(opk_vspace_t* vspace,
               opk_vector_t* dst, const opk_vector_t* src);

  /* Exchange contents of vector SRC and vector DST.  When this method is
     called, it is guaranteed that DST and SRC are different and both belong
     to VSPACE. */
  void (*swap)(opk_vspace_t* vspace,
               opk_vector_t* x, opk_vector_t* y);

  /* Scale vector by a scalar.  This methods performs the operation:
     DST = ALPHA*SRC.   When this method is called, it is guaranteed that
     DST and SRC belong to VSPACE and that ALPHA is neither 0 nor 1. */
  void (*scale)(opk_vspace_t* vspace,opk_vector_t* dst,
                double alpha, const opk_vector_t* src);

  /* Compute the linear combination of two elements of this vector space
     (DST = ALPHA*X + BETA*Y).   When this method is called, it is
     guaranteed that all vectors (DST, X and Y) belong to VSPACE and that
     neither ALPHA, nor BETA are equal to zero. */
  void (*axpby)(opk_vspace_t* vspace, opk_vector_t* dst,
                double alpha, const opk_vector_t* x,
                double beta,  const opk_vector_t* y);

  /* Compute the linear combination of three elements of this vector space
     (DST = ALPHA*X + BETA*Y + GAMMA*Z).  When this method is called, it is
     guaranteed that all vectors (DST, X, Y and Z) belong to VSPACE and that
     none of ALPHA, BETA nor GAMMA are equal to zero. */
  void (*axpbypcz)(opk_vspace_t* vspace, opk_vector_t* dst,
                   double alpha, const opk_vector_t* x,
                   double beta,  const opk_vector_t* y,
                   double gamma, const opk_vector_t* z);

};

/* `opk_allocate_vector_space()` allocates `nbytes` bytes of memory to store a
   vector space derived type and initializes it as a basic `opk_vspace_t` type.
   The value `nbytes` is adjusted to be at least `sizeof(opk_vspace_t)`.
   `NULL` is returned in case of error.  The returned object must be
   unreferenced `opk_drop_object()` or macro `OPK_DROP()` when no longer
   needed. */
extern opk_vspace_t*
opk_allocate_vector_space(const opk_vspace_operations_t* vops,
                          opk_index_t nvariables,
                          size_t nbytes);

/**
 * Allocate a new vector object.
 *
 * This utility function is designed for constructors of vector sub-types,
 * that is the `create()` method of a vector space table of operations.  To
 * create vectors, the end-users use the `opk_vcreate()` function.
 *
 * This function allocates dynamic memory and instanciates it with as a basic
 * vector type.  The storage size is adjusted to be at least sufficient for a
 * basic vector of type `opk_vector_t`.  As part of its initialization, the
 * returned vector holds a reference on its vector space.  The returned vector
 * is managed as any other OptimPack object; in particular, the function
 * `opk_drop_object()` or the macro `OPK_DROP()` has to be called to drop the
 * reference on the vector.
 *
 * @param vspace - The owner of the vector.
 * @param nbytes - The minimum number of bytes to allocate.
 *
 * @return A new vector object instantiated as a basic vector type; `NULL` in
 *         case of error.
 */
extern opk_vector_t*
opk_allocate_vector(opk_vspace_t* vspace, size_t nbytes);


/*---------------------------------------------------------------------------*/
/* LOW LEVEL API FOR LINE SEARCH AND DERIVED TYPES */

typedef struct _opk_lnsrch_operations opk_lnsrch_operations_t;

/* The base structure for line search must be exposed for line search
   methods. */
struct _opk_lnsrch {
  opk_object_t base;            /* Base type (must be the first member). */
  opk_lnsrch_operations_t *ops; /* Table of line search methods. */
  double stp;                   /* Current step length. */
  double stpmin;                /* Lower bound for the step. */
  double stpmax;                /* Upper bound for the step. */
  double finit;                 /* Function value at the start of the
                                   search. */
  double ginit;                 /* Directional derivative value at the start of
                                   the search. */
  int status;                   /* Last value returned by line search
                                   methods. */
  int searching;                /* True if search is in progress. */
};

struct _opk_lnsrch_operations {
  /* Method used to release any ressources specifically allocated by the
     line search method (not the workspace itself). */
  void (*finalize)(opk_lnsrch_t* self);

  /* Method called by opk_lnsrch_start() to initiate a line search.  This
     method is called after having set member STP with the length of the first
     step to perform, members STPMIN and STPMAX with the lower and upper bounds
     of the step length, members FINIT and GINIT with the function value and
     the derivative the function along the search direction at the start of the
     search. */
  int (*start)(opk_lnsrch_t* self);

  /* Method to iterate during a line search.  STP_PTR, F_PTR and D_PTR are the
     addresses of variables which store STP, F and D.  On entry, STP is the
     current step length, F is the function value for this step and D is the
     corresponding derivative of the function along the search direction.
     FIXME: On exit, if convergence is achieved (or in case of error/warning)
     STP, F and D are left unchanged; otherwise STP is the new step to try. */
  int (*iterate)(opk_lnsrch_t* self,
                 double* stp_ptr, double f1, double d1);
};

extern opk_lnsrch_t*
opk_allocate_line_search(opk_lnsrch_operations_t *ops,
                         size_t size);

/*---------------------------------------------------------------------------*/
/* LOW LEVEL API FOR OPERATORS AND DERIVED TYPES */

typedef struct _opk_operator_operations opk_operator_operations_t;

extern opk_operator_t*
opk_allocate_operator(const opk_operator_operations_t* ops,
                      opk_vspace_t* inpspace,
                      opk_vspace_t* outspace,
                      size_t size);

struct _opk_operator {
  opk_object_t base;      /* base type (must be the first member) */
  const opk_operator_operations_t* ops; /* table of operator methods */
  opk_vspace_t* inpspace; /* input space of the operator */
  opk_vspace_t* outspace; /* output space of the operator */
};

struct _opk_operator_operations {

  /* If not NULL, the finalize() method is called to finalize any specific
     parts of the operator sub-type before finalizing the base operator
     type. */
  void (*finalize)(opk_operator_t* self);

  /* If not NULL, the apply_direct() method is called by opk_apply_direct() to
     apply the operator to the source vector `src` and store the result in the
     destination vector `dst`.  It is guaranteed that the arguments are valid
     (i.e., that they are non NULL and that they belong to the correct vector
     space).  If NULL, such operation is considered as not allowed by the
     operator. */
  int (*apply_direct)(opk_operator_t* self, opk_vector_t* dst,
                      const opk_vector_t* src);

  /* If not NULL, the apply_adjoint() method is called by opk_apply_adjoint()
     to apply the adjoint of the operator to the source vector `src` and store
     the result in the destination vector `dst`.  The same assumptions as for
     the apply_direct() method otherwise hold. */
  int (*apply_adjoint)(opk_operator_t* self, opk_vector_t* dst,
                       const opk_vector_t* src);

  /* If not NULL, the apply_inverse() method is called by opk_apply_adjoint()
     to apply the inverse of the operator to the source vector `src` and store
     the result in the destination vector `dst`.  The same assumptions as for
     the apply_direct() method otherwise hold. */
  int (*apply_inverse)(opk_operator_t* self, opk_vector_t* dst,
                       const opk_vector_t* src);
};

#endif /* _OPTIMPACK_PRIVATE_H */

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
