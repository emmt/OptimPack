/*
 * optimpack-private.h --
 *
 * Private definitions for building OptimPack library.
 *
 *----------------------------------------------------------------------------------------
 *
 * This file is part of OptimPack (https://github.com/emmt/OptimPack).
 *
 * Copyright (C) 2014-2016 Éric Thiébaut
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this
 * software and associated documentation files (the "Software"), to deal in the Software
 * without restriction, including without limitation the rights to use, copy, modify,
 * merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies
 * or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 *----------------------------------------------------------------------------------------
 */

#ifndef OPTIMPACK_PRIVATE_H_
#define OPTIMPACK_PRIVATE_H_ 1

#include "optimpack.h"


/*--------------------------------------------------------------------------------------*/
/* LOW LEVEL API FOR BASIC OBJECTS AND DERIVED TYPES */

/**
 * @defgroup LowObjects  Allocation of basic objects.
 * @ingroup LowLevel
 * @{
 */

/**
 * Basic object structure.
 *
 * The definition of this structure must be exposed so that others types can "inherit"
 * from it. However, the contents must be left unchanged otherwise unexpected results may
 * occur.
 */
struct opk_object_ {
  void (*finalize)(opk_object* self);
  opk_index references;
};

/**
 * Allocate a new (empty) object.
 *
 * This function allocates enough bytes of memory to store a type derived from the basic
 * object type, initializes the basic type and zero-fill the extra bytes of allocated
 * memory.
 *
 * The caller of this function holds a reference on the returned object. When the object
 * is no longer needed by the caller, he/she has to call `opk_drop_object()` to release
 * this reference.
 *
 * When the object is no longer in use (its last reference has been dropped), it is
 * effectively destroyed by first applying the `finalize()` method on it (unless it is
 * NULL) and then freeing the dynamic memory where the object was stored.
 *
 * The typical usage consists in building a derived type as follows:
 *
 * ~~~~~~~~~~{.c}
 * typedef struct {
 *     opk_object base; // base type (must be the first member)
 *     ...;               // other members
 * } sub_type;
 * ~~~~~~~~~~
 *
 * and provide a destructor and a constructor like:
 *
 * ~~~~~~~~~~{.c}
 * static void finalize_sub_type(opk_object* self)
 * {
 *     ...; // do whatever cleanup is needed
 * }
 *
 * sub_type* create_sub_type(args...) {
 *    sub_type* obj = (sub_type*)opk_allocate_object(finalize_sub_type,
 *                                                       sizeof(sub_type));
 *    ...; // any other initialization
 *    return obj;
 * }
 * ~~~~~~~~~~
 *
 * where `finalize_sub_type()` is in charge of releasing any specific ressources of the
 * object. Note that this function must be able to deal with a partially initialized
 * object, although this is simplified by the fact that the non-basic parts of the object
 * are intially zero-filled. The memory allocated by `opk_allocate_object()` is
 * automatically freed and must not be freed by the `finalize()` method.
 *
 * The object model of OptimPack is very simple but has some drawbacks (C is not C++):
 * casts are needed to force the proper type of an object pointer, memory management by
 * reference counting forbid circular references.
 *
 * @param finalize - If not NULL, method to call when the object is no longer referenced.
 *                   This method is in charge of freeing object ressources specific to the
 *                   sub-type. This method is *never* called with a NULL argument.
 *
 * @param nbytes   - Number of bytes to allocate for the object. This value may be
 *                   adjusted so that at least a basic object can be stored.
 *
 * @return A new basic object with one reference set or NULL in case of error.
 */
extern opk_object*
opk_allocate_object(void (*finalize)(opk_object* self),
                    size_t nbytes);

/** @} */

/*--------------------------------------------------------------------------------------*/
/* LOW LEVEL API FOR VECTOR SPACES AND DERIVED TYPES */

/**
 * @defgroup LowVectors     Implementing vectors based on specific memory model.
 * @ingroup LowLevel
 * @{
 */

/** Table of methods for vector spaces. */
typedef struct opk_vspace_operations_ opk_vspace_operations;

/**
 * Structure implementing a basic vector type.
 *
 * The `opk_vector` structure is intentionally exposed to let different implementations
 * coexist (although in separate codes). If one want to implement a sub-type of the vector
 * type, it is sufficient to define a new structure whose first member is an `opk_vector`.
 * The function `opk_allocate_vector()` **must** be used to allocate the whole structure.
 * For instance:
 *
 * ~~~~~~~~~~{.c}
 * typedef struct {
 *   opk_vector base;
 *   double* data;
 * } my_vector;
 *
 * void finalize_vector(my_vector* v)
 * {
 *   if (v->data != NULL) {
 *     free(v->data);
 *   }
 * }
 *
 * my_vector* new_vector(int n)
 * {
 *    my_vector* v = (my_vector*)opk_allocate_vector(space, sizeof(my_vector));
 *    if (v != NULL) {
 *       v->data = (double*)malloc(n*sizeof(double));
 *       if (v->data == NULL) {
 *         OPK_DROP(v);
 *         return NULL;
 *       }
 *    }
 *    return v;
 * }
 * ~~~~~~~~~~
 *
 * OptimPack routines only require the address of such vectors and treat them as opaque
 * structures.
 */
struct opk_vector_ {
  opk_object  base;   /**< Base type (must be the first member). */
  opk_vspace* owner;  /**< Vector space to which the vector belongs. */
};

struct opk_vspace_ {
  opk_object base;    /* base type (must be the first member) */

  /* Table of operations on vectors of this space. */
  const opk_vspace_operations* ops;

  /* Number of elements of vectors of this vector space (real problems are in finite
     dimension!). */
  opk_index size;
};

struct opk_vspace_operations_ {
  /* A short description of the vector space. */
  const char* description;

  /* Method called to effectively free anything else than the memory allocated to store
     the vector space itself. It is guaranteed that this method is called with a non-NULL
     argument and that no references exist on the vector space (in particular all vectors
     of this vector space have been finalized). */
  void (*finalize_space)(opk_vspace* vspace);

  /* Create an element of this vector space (the contents of the vector is undefined). */
  opk_vector* (*create)(opk_vspace* vspace);

  /* Delete an element of this vector space. Same semantics as the finalize_space()
     method. FIXME: improve explanations */
  void (*finalize_vector)(opk_vspace* vspace, opk_vector* v);

  /* Get a specific value in a vector. This function is only needed for debugging. When
     called, it is guaranteed, that `v` belongs to `vspace` and that index `k` is valid (0
     <= k < vspace->size). */
  double (*peek)(const opk_vspace* vspace, const opk_vector* v,
                 opk_index k);

  /* Set a specific value in a vector. This function is only needed for debugging. When
     called, it is guaranteed, that `v` belongs to `vspace` and that index `k` is valid (0
     <= k < vspace->size). */
  void (*poke)(const opk_vspace* vspace, opk_vector* v,
               opk_index k, double value);

  /* Copy the values of a conventional array into a vector. */
  void (*import)(const opk_vspace* space, opk_vector* dst,
                 const void* src, opk_eltype type);

  /* Copy the values of a vector in a conventional array. */
  void (*export)(const opk_vspace* space, void* dst, opk_eltype type,
                 const opk_vector* src);

  /* Fill a vector with a scalar ALPHA. When this method is called, it is guaranteed that
     X belong to VSPACE. The implementation of this method may be optimized when ALPHA is
     equal to zero. */
  void (*fill)(opk_vspace* vspace, opk_vector* v, double alpha);

  /* Compute the L1 norm (sum of absolute values) of an element of this vector space. When
     this method is called, it is guaranteed that X belongs to VSPACE. */
  double (*norm1)(opk_vspace* vspace, const opk_vector* x);

  /* Compute the L2 norm (square root of the sum of squared values, also known as the
     Euclidean norm) of an element of this vector space. When this method is called, it is
     guaranteed that X belongs to VSPACE. */
  double (*norm2)(opk_vspace* vspace, const opk_vector* x);

  /* Compute the infinite norm (maximum absolute value) of an element of this vector
     space. When this method is called, it is guaranteed that X belongs to VSPACE. */
  double (*norminf)(opk_vspace* vspace, const opk_vector* x);

  /* Compute the inner product between two elements of this vector space. When this method
     is called, it is guaranteed that X and Y belong to VSPACE. */
  double (*dot)(opk_vspace* vspace,
                const opk_vector* x,
                const opk_vector* y);

  /* Compute the inner product between three elements of this vector space defined as
     sum_i W[i]*X[i]*Y[i]. When this method is called, it is guaranteed that W, X and Y
     belong to VSPACE. */
  double (*dot3)(opk_vspace* vspace,
                 const opk_vector* w,
                 const opk_vector* x,
                 const opk_vector* y);

  /* Copy vector SRC to vector DST. When this method is called, it is guaranteed that DST
     and SRC are different and both belong to VSPACE. */
  void (*copy)(opk_vspace* vspace,
               opk_vector* dst, const opk_vector* src);

  /* Exchange contents of vector SRC and vector DST. When this method is called, it is
     guaranteed that DST and SRC are different and both belong to VSPACE. */
  void (*swap)(opk_vspace* vspace,
               opk_vector* x, opk_vector* y);

  /* Scale vector by a scalar. This methods performs the operation: DST = ALPHA*SRC. When
     this method is called, it is guaranteed that DST and SRC belong to VSPACE and that
     ALPHA is neither 0 nor 1. */
  void (*scale)(opk_vspace* vspace,opk_vector* dst,
                double alpha, const opk_vector* src);

  /* Compute the elementwise product of two vectors: DST[i] = X[i]*Y[i] for all i. When
     this method is called, it is guaranteed that DST, X and Y belong to VSPACE. */
  void (*product)(opk_vspace* vspace, opk_vector* dst,
                  const opk_vector* x, const opk_vector* y);

  /* Compute the linear combination of two elements of this vector space (DST = ALPHA*X +
     BETA*Y). When this method is called, it is guaranteed that all vectors (DST, X and Y)
     belong to VSPACE and that neither ALPHA, nor BETA are equal to zero. */
  void (*axpby)(opk_vspace* vspace, opk_vector* dst,
                double alpha, const opk_vector* x,
                double beta,  const opk_vector* y);

  /* Compute the linear combination of three elements of this vector space (DST = ALPHA*X
     + BETA*Y + GAMMA*Z). When this method is called, it is guaranteed that all vectors
     (DST, X, Y and Z) belong to VSPACE and that none of ALPHA, BETA nor GAMMA are equal
     to zero. */
  void (*axpbypcz)(opk_vspace* vspace, opk_vector* dst,
                   double alpha, const opk_vector* x,
                   double beta,  const opk_vector* y,
                   double gamma, const opk_vector* z);

  /* Project variables X such that DST = MAX(XL, MIN(X, XU)). There are 9 possibilities
     for the bounds. If (BOUND%3) = 0, there is no lower bound; if (BOUND%3) = 1, the
     lower bound is a scalar; if (BOUND%3) = 2, the lower bound is a vector. If (BOUND/3)
     = 0, there is no uppwer bound; if (BOUND/3) = 1, the upper bound is a scalar; if
     (BOUND/3) = 2, the upper bound is a vector. */
  opk_status (*boxprojvar)(opk_vspace* space, opk_vector* dst,
                             const opk_vector* x,
                             const void* xl, const void* xu, int bound);

  opk_status (*boxprojdir)(opk_vspace* space, opk_vector* dst,
                             const opk_vector* x,
                             const void* xl, const void* xu, int bound,
                             const opk_vector* d, opk_orientation orient);

  opk_status (*boxfreevar)(opk_vspace* space, opk_vector* dst,
                             const opk_vector* x,
                             const void* xl, const void* xu, int bound,
                             const opk_vector* d, opk_orientation orient);

  opk_status (*boxsteplim)(opk_vspace* space,
                             double* smin1, double* smin2, double* smax,
                             const opk_vector* x,
                             const void* xl, const void* xu, int bound,
                             const opk_vector* d,
                             opk_orientation orient);

};

/* `opk_allocate_vector_space()` allocates `nbytes` bytes of memory to store a vector
   space derived type and initializes it as a basic `opk_vspace` type. The value `nbytes`
   is adjusted to be at least `sizeof(opk_vspace)`. `NULL` is returned in case of error.
   The returned object must be unreferenced `opk_drop_object()` or macro `OPK_DROP()` when
   no longer needed. */
extern opk_vspace*
opk_allocate_vector_space(const opk_vspace_operations* vops,
                          opk_index nvariables,
                          size_t nbytes);

/**
 * Allocate a new vector object.
 *
 * This utility function is designed for constructors of vector sub-types, that is the
 * `create()` method of a vector space table of operations. To create vectors, the
 * end-users use the `opk_vcreate()` function.
 *
 * This function allocates dynamic memory and instanciates it with as a basic vector type.
 * The storage size is adjusted to be at least sufficient for a basic vector of type
 * `opk_vector`. As part of its initialization, the returned vector holds a reference on
 * its vector space. The returned vector is managed as any other OptimPack object; in
 * particular, the function `opk_drop_object()` or the macro `OPK_DROP()` has to be called
 * to drop the reference on the vector.
 *
 * @param vspace - The owner of the vector.
 * @param nbytes - The minimum number of bytes to allocate.
 *
 * @return A new vector object instantiated as a basic vector type; `NULL` in
 *         case of error.
 */
extern opk_vector*
opk_allocate_vector(opk_vspace* vspace, size_t nbytes);

/** @} */

/*--------------------------------------------------------------------------------------*/
/* LOW LEVEL API FOR SEPARABLE BOUND CONSTRAINTS */

/**
 * @defgroup LowConvexSet Implementing convex sets.
 * @ingroup LowLevel
 * @{
 */

/**
 * Private structure to store the base of an instance derived from a convex set.
 */
struct opk_convexset_ {
  opk_object base;       /**< Base type (must be the first member). */
  opk_vspace* space;     /**< Variable space. */
  void (*finalize)(opk_convexset* self);
  opk_status (*projvar)(opk_vector* dst,
                          const opk_vector* x,
                          const opk_convexset* set);
  opk_status (*projdir)(opk_vector* dst,
                          const opk_vector* x,
                          const opk_convexset* set,
                          const opk_vector* d,
                          opk_orientation orient);
  opk_status (*freevar)(opk_vector* dst,
                          const opk_vector* x,
                          const opk_convexset* set,
                          const opk_vector* d,
                          opk_orientation orient);
  opk_status (*steplim)(double* smin1, double* smin2, double *smax,
                          const opk_vector* x,
                          const opk_convexset* set,
                          const opk_vector* d,
                          opk_orientation orient);
};

/**
 * Allocate a new instance of a convex set.
 */
extern opk_convexset*
opk_allocate_convexset(opk_vspace* space,
                       void (*finalize)(opk_convexset* self),
                       opk_status (*projvar)(opk_vector* dst,
                                               const opk_vector* x,
                                               const opk_convexset* set),
                       opk_status (*projdir)(opk_vector* dst,
                                               const opk_vector* x,
                                               const opk_convexset* set,
                                               const opk_vector* d,
                                               opk_orientation orient),
                       opk_status (*freevar)(opk_vector* dst,
                                               const opk_vector* x,
                                               const opk_convexset* set,
                                               const opk_vector* d,
                                               opk_orientation orient),
                       opk_status (*steplim)(double* smin1, double* smin2,
                                               double *smax,
                                               const opk_vector* x,
                                               const opk_convexset* set,
                                               const opk_vector* d,
                                               opk_orientation orient),
                       size_t nbytes);

/** @} */

/*--------------------------------------------------------------------------------------*/
/* LOW LEVEL API FOR LINE SEARCH AND DERIVED TYPES */

/**
 * @defgroup LowLineSearch  Implementing line search methods.
 * @ingroup LowLevel
 * @{
 */

typedef struct opk_lnsrch_operations_ opk_lnsrch_operations;

/* The base structure for line search must be exposed for line search
   methods. */
struct opk_lnsrch_ {
  opk_object base;            /**< Base type (must be the first member). */
  opk_lnsrch_operations *ops; /**< Table of line search methods. */
  double stp;                 /**< Current step length. */
  double stpmin;              /**< Lower bound for the step. */
  double stpmax;              /**< Upper bound for the step. */
  double finit;               /**< Function value at the start of the search. */
  double ginit;               /**< Directional derivative value at the start of the search. */
  opk_status status;          /**< Last status. */
  opk_lnsrch_task task;       /**< Current pending task. */
  int searching;              /**< True if search is in progress. */
};

struct opk_lnsrch_operations_ {
  /* Method used to release any ressources specifically allocated by the line search
     method (not the workspace itself). */
  void (*finalize)(opk_lnsrch* self);

  /* Method called by opk_lnsrch_start() to initiate a line search. This method is called
     after having set member STP with the length of the first step to perform, members
     STPMIN and STPMAX with the lower and upper bounds of the step length, members FINIT
     and GINIT with the function value and the derivative the function along the search
     direction at the start of the search. */
  opk_lnsrch_task (*start)(opk_lnsrch* self);

  /* Method to iterate during a line search. STP_PTR is the addresse of the variable which
     store the current step size STP. F is the function value for this step size. DF is
     the corresponding directional derivative along the search direction. On exit, if
     convergence is achieved (or in case of error/warning) the step size is left
     unchanged; otherwise, the contents of STP_PTR is set with the new step to try. */
  opk_lnsrch_task (*iterate)(opk_lnsrch* self,
                               double* stp_ptr, double f, double df);

  /* Flag to indicate whether directional derivative are needed to check the line search
     convergence; else only function values are used (e.g. Armijo or 1st Wolfe condition
     only). The initial directional derivative at the start of each search is always
     needed. */
  opk_bool use_deriv;
};

extern opk_lnsrch*
opk_allocate_line_search(opk_lnsrch_operations *ops,
                         size_t size);

extern void
opk__set_nlcg_status(opk_nlcg* opt, opk_status status);

extern void
opk__set_vmlmb_status(opk_vmlmb* opt, opk_status status);

/** @} */

/*--------------------------------------------------------------------------------------*/
/* LOW LEVEL API FOR OPERATORS AND DERIVED TYPES */

/**
 * @defgroup LowOperators  Implementing operators acting on vectors.
 * @ingroup LowLevel
 * @{
 */

typedef struct opk_operator_operations_ opk_operator_operations;

extern opk_operator*
opk_allocate_operator(const opk_operator_operations* ops,
                      opk_vspace* inpspace,
                      opk_vspace* outspace,
                      size_t size);

struct opk_operator_ {
  opk_object base;      /* base type (must be the first member) */
  const opk_operator_operations* ops; /* table of operator methods */
  opk_vspace* inpspace; /* input space of the operator */
  opk_vspace* outspace; /* output space of the operator */
};

struct opk_operator_operations_ {

  /* If not NULL, the finalize() method is called to finalize any specific parts of the
     operator sub-type before finalizing the base operator type. */
  void (*finalize)(opk_operator* self);

  /* If not NULL, the apply_direct() method is called by opk_apply_direct() to apply the
     operator to the source vector `src` and store the result in the destination vector
     `dst`. It is guaranteed that the arguments are valid (i.e., that they are non NULL
     and that they belong to the correct vector space). If NULL, such operation is
     considered as not allowed by the operator. */
  opk_status (*apply_direct)(opk_operator* self, opk_vector* dst,
                               const opk_vector* src);

  /* If not NULL, the apply_adjoint() method is called by opk_apply_adjoint() to apply the
     adjoint of the operator to the source vector `src` and store the result in the
     destination vector `dst`. The same assumptions as for the apply_direct() method
     otherwise hold. */
  opk_status (*apply_adjoint)(opk_operator* self, opk_vector* dst,
                                const opk_vector* src);

  /* If not NULL, the apply_inverse() method is called by opk_apply_adjoint() to apply the
     inverse of the operator to the source vector `src` and store the result in the
     destination vector `dst`. The same assumptions as for the apply_direct() method
     otherwise hold. */
  opk_status  (*apply_inverse)(opk_operator* self, opk_vector* dst,
                                 const opk_vector* src);
};

/** @} */

#endif /* OPTIMPACK_PRIVATE_H_ */
