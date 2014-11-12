/*
 * optimpack.h --
 *
 * Common definitions for OptimPack library.
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

#ifndef _OPTIMPACK_H
#define _OPTIMPACK_H 1

#include <stddef.h>

/*---------------------------------------------------------------------------*/

/*
 * Defines: Miscellaneous macros.
 *
 *   OPK_BEGIN_C_DECLS - Start declarations of C elements (functions, variables,
 *                       types, structures, etc.).
 *   OPK_END_C_DECLS   - Terminate C declarations.
 *
 *   <OPK_BEGIN_C_DECLS> should be used at the beginning of your declarations,
 *   so that C++ compilers don't mangle their names.  Use <OPK_END_C_DECLS> at
 *   the end of C declarations.
 */
#ifdef __cplusplus
# define OPK_BEGIN_C_DECLS extern "C" {
# define OPK_END_C_DECLS   }
#else /* not __cplusplus */
# define OPK_BEGIN_C_DECLS
# define OPK_END_C_DECLS
#endif /* __cplusplus */


/*
 * Macro: OPK_JOIN
 *
 * OPK_JOIN(_a_,_b_) joins its two arguments, possibly after performing macro
 *  expansion of _a_ and _b_.
 */
#define OPK_JOIN(a,b)                _OPK_JOIN2(a,b)

#define OPK_JOIN2(a,b)               _OPK_JOIN2(a,b)
#define OPK_JOIN3(a,b,c)             _OPK_JOIN3(a,b,c)
#define OPK_JOIN4(a,b,c,d)           _OPK_JOIN4(a,b,c,d)
#define OPK_JOIN5(a,b,c,d,e)         _OPK_JOIN5(a,b,c,d,e)
#define OPK_JOIN6(a,b,c,d,e,f)       _OPK_JOIN6(a,b,c,d,e,f)
#define OPK_JOIN7(a,b,c,d,e,f,g)     _OPK_JOIN7(a,b,c,d,e,f,g)
#define OPK_JOIN8(a,b,c,d,e,f,g,h)   _OPK_JOIN8(a,b,c,d,e,f,g,h)
#define OPK_JOIN9(a,b,c,d,e,f,g,h,i) _OPK_JOIN8(a,b,c,d,e,f,g,h,i)

/*
 * Macro: OPK_STRINGIFY
 *
 * OPK_STRINGIFY(x) takes its argument _x_ and wraps it in "" (double
 * quotation marks), possibly after performing macro expansion of the
 * argument.
 */
#define OPK_STRINGIFY(x)  _OPK_STRINGIFY(x)

/* Stringify/concatenate arguments without expansion. */
#if defined(__STDC__) || defined(__cplusplus) || defined(c_plusplus)
# define _OPK_STRINGIFY(x)             # x
# define _OPK_JOIN2(a,b)               a##b
# define _OPK_JOIN3(a,b,c)             a##b##c
# define _OPK_JOIN4(a,b,c,d)           a##b##c##d
# define _OPK_JOIN5(a,b,c,d,e)         a##b##c##d##e
# define _OPK_JOIN6(a,b,c,d,e,f)       a##b##c##d##e##f
# define _OPK_JOIN7(a,b,c,d,e,f,g)     a##b##c##d##e##f##g
# define _OPK_JOIN8(a,b,c,d,e,f,g,h)   a##b##c##d##e##f##g##h
# define _OPK_JOIN9(a,b,c,d,e,f,g,h,i) a##b##c##d##e##f##g##h##i
#else
# define _OPK_STRINGIFY(x)             "x"
# define _OPK_JOIN2(a,b)               a/**/b
# define _OPK_JOIN3(a,b,c)             a/**/b/**/c
# define _OPK_JOIN4(a,b,c,d)           a/**/b/**/c/**/d
# define _OPK_JOIN5(a,b,c,d,e)         a/**/b/**/c/**/d/**/e
# define _OPK_JOIN6(a,b,c,d,e,f)       a/**/b/**/c/**/d/**/e/**/f
# define _OPK_JOIN7(a,b,c,d,e,f,g)     a/**/b/**/c/**/d/**/e/**/f/**/g
# define _OPK_JOIN8(a,b,c,d,e,f,g,h)   a/**/b/**/c/**/d/**/e/**/f/**/g/**/h
# define _OPK_JOIN9(a,b,c,d,e,f,g,h,i) a/**/b/**/c/**/d/**/e/**/f/**/g/**/h/**/i
#endif


/*---------------------------------------------------------------------------*/

/*
 * Macro: OPK_ABS
 *
 * OPK_ABS(_a_) gives the absolute value of _a_.
 */
#define OPK_ABS(a)   ((a) >= 0 ? (a) : -(a))

/*
 * Macro: OPK_MIN
 *
 * OPK_MIN(_a_, _b_) yields the minimum value between _a_ and _b_.
 */
#define OPK_MIN(a,b) ((a) <= (b) ? (a) : (b))

/*
 * Macro: OPK_MAX
 *
 * OPK_MAX(_a_, _b_) yields the maximum value between _a_ and _b_.
 */
#define OPK_MAX(a,b) ((a) >= (b) ? (a) : (b))

/*
 * Macro: OPK_HOW_MANY
 *
 * OPK_HOW_MANY(_a_, _b_) yields the minimal number of chunks with _b_
 * elements needed to store _a_ elements.
 */
#define OPK_HOW_MANY(a, b) ((((b) - 1) + (a))/(b))

/*
 * Macro: OPK_ROUND_UP
 *
 * OPK_ROUND_UP(_a_, _b_) yields the integer _a_ rounded up to a multiple of
 * integer _b_.
 */
#define OPK_ROUND_UP(a, b) (OPK_HOW_MANY(a, b)*(b))

/*
 * Macro: OPK_LOOP
 *
 * OPK_LOOP(_var_, _cnt_) yields a simple loop over variable _var_ from 0 to
 * _cnt_-1.
 */
#define OPK_LOOP(var, cnt) for (var = 0; var < cnt; ++var)


/* Yields |a|*sign(b). */
#define opk_sign(a, b) ({ typeof(a) _1 = (a); typeof(b) _2 = (b); \
                          ((_1 < 0) == ((b) < 0)) ? _1 : -_1;})

/* Yield absolute value. */
#define opk_abs(a)    ({ typeof(a) _1 = (a); _1 >= 0 ? _1 : -_1; })

/* Yield minimum of two values. */
#define opk_min(a, b) ({ typeof(a) _1 = (a); typeof(b) _2 = (b); \
                         _1 <= _2 ? _1 : _2;})

/* Yield maximum of two values. */
#define opk_max(a, b) ({ typeof(a) _1 = (a); typeof(b) _2 = (b); \
                         _1 >= _2 ? _1 : _2;})
/*---------------------------------------------------------------------------*/

/*
 * Macro: OPK_NEW
 *
 * OPK_NEW(_type_) allocates an object of structure/class _type_.
 *
 * See Also:
 *   <OPK_ALLOC>, <OPK_FREE>, malloc(3), free(3).
 */
#define OPK_NEW(type)             ((type *)malloc(sizeof(type)))

/*
 * Macro: OPK_ALLOC
 *
 * OPK_ALLOC(_type_, _number_) allocates an array of _number_ items of
 * structure/class _type_.
 *
 * See Also:
 *   <OPK_NEW>, <OPK_FREE>, malloc(3), free(3).
 */
#define OPK_ALLOC(type, number)   ((type *)malloc((number)*sizeof(type)))

/*
 * Macro: OPK_FREE
 *
 * OPK_FREE(_ptr_) frees pointer _ptr_ if non-null.
 *
 * See Also:
 *   <OPK_ALLOC>, <OPK_NEW>, malloc(3), free(3).
 */
#define OPK_FREE(ptr)   if (! (ptr)) ; else free(ptr)

/*
 * Macro: OPK_OFFSET
 *
 * OPK_OFFSET(_type_, _member_) yields offset (in bytes) of _member_ in
 * structure/class _type_.
 */
#define OPK_OFFSET(type, member)  ((char*)&((type*)0)->member - (char*)0)

/*---------------------------------------------------------------------------*/

OPK_BEGIN_C_DECLS

/*
 * Enum:  opk_status_t
 *
 *  OPK_FAILURE   - Status value returned upon success.
 *  OPK_SUCCESS   - Status value returned upon failure.
 *
 *  Values returned by SPL routines.
 */
typedef enum {
  OPK_FAILURE = -1,
  OPK_SUCCESS = 0
} opk_status_t;

/*---------------------------------------------------------------------------*/
/* DATA TYPES */

/*
 * Type: opk_index_t
 *
 *   Integer data type used for array indices in SPL.
 */
typedef ptrdiff_t opk_index_t;

/*
 * Type: opk_bool_t
 *
 *   Data type for boolean (logical) values.
 */
typedef int opk_bool_t;

#define OPK_TRUE  1
#define OPK_FALSE 0

/*---------------------------------------------------------------------------*/
/* VECTOR SPACES */

typedef struct _opk_vector opk_vector_t; /* a vector */
typedef struct _opk_vspace opk_vspace_t; /* a vector space */

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
  opk_vspace_t* owner;  /* the vector space to which the vector belongs */
};

struct _opk_vspace {
  /* Any unique (for instance the address of a static variable) address
     which can be used to check the type of the vector space. */
  const void* ident;

  /* Number of elements of vectors of this vector space (real problems are in
     finite dimension!). */
  opk_index_t size;

  /* Method called to effectively destroy the vector space.  For now, the user
     is responsible of not destroying the vector space until all its elements
     have been deleted (this could be handled by reference counting, but this
     is not yet implemented). */
  void (*finalize)(opk_vspace_t* vspace);

  /* Create an element of this vector space (the contents of the vector is
     undefined). */
  opk_vector_t* (*create)(opk_vspace_t* vspace);

  /* Delete an element of this vector space. */
  void (*delete)(opk_vspace_t* vspace, opk_vector_t* v);

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

/*---------------------------------------------------------------------------*/
/* LOW-LEVEL ROUTINES FOR VECTOR SPACES */

/* `opk_valloc()` allocates `size` bytes of memory to store a vector type and
   initializes it as a basic `opk_vector_t` type.  The value `size` is adjusted
   to be at least `sizeof(opk_vector_t)`.  `NULL` is returned in case of error.
   The allocated memory must be freed with `opk_vfree()`. */
extern opk_vector_t*
opk_valloc(opk_vspace_t* vspace, size_t size);

/* `opk_vfree()` must be used to free the memory allocated by `opk_valloc()`.
   It does not free other ressources. */
extern void
opk_vfree(opk_vector_t* vect);

/*---------------------------------------------------------------------------*/
/* HIGH-LEVEL ROUTINES FOR VECTOR SPACES */

typedef void opk_free_proc(void*);

/**
 * Create a vector space for array of double's in conventional memory.
 *
 * This particular type of vector space deals with arrays of values stored
 * contiguously and accessible as conventional arrays.  This include arrays
 * allocated from the heap, dynamically allocated with `malloc()`, arrays in
 * shared memory, memory mapped files, etc.
 *
 * To create vectors belonging to this kind of vector space, as for any type
 * of vector spaces, it is possible to call `opk_vcreate()` but it is also
 * possible to call `opk_wrap_simple_double_vector()` to wrap an existing
 * array (of the correct size and type of course) into a vector.
 *
 * @param size - The number of elements of the vectors of the space.
 *
 * @return A new vector space or `NULL` in case of errors.
 */
extern opk_vspace_t*
opk_new_simple_double_vector_space(opk_index_t size);

/**
 * Wrap an existing array into a simple vector.
 *
 * This function creates a new vector whose elements are stored into an array
 * provided by the caller.  The caller is responsible of ensuring that the
 * memory is sufficiently large (the array has at least `vspace->size`
 * elements) and correctly aligned.
 *
 * When the vector is destroyed, the function `free_client_data()`, if not
 * `NULL`, is called with argument `client_data` to release ressources.  Then
 * the container is freed.  If function `free_client_data()` is `NULL`, it is
 * assumed that the caller is responsible of releasing the data when no longer
 * needed.
 *
 * A typical usage is:
 * <pre>
 *     #define N 1000
 *     opk_vspace_t* vspace = opk_new_simple_double_vector_space(N);
 *     double heap_array[N];
 *     opk_vector_t* v1 = opk_wrap_simple_double_vector(vspace, heap_array,
 *                                                      NULL, NULL);
 *     double* dynamic_array = (double*)malloc(N*sizeof(double));
 *     opk_vector_t* v2 = opk_wrap_simple_double_vector(vspace, dynamic_array,
 *                                                      dynamic_array, free);
 * </pre>
 *
 * which creates two vectors `v1` and `v2` which are respectively wrapped
 * around an array allocated on the heap and around a dynamically allocated
 * array.
 *
 * In the above example, the `client_data` and the `data` are the same but the
 * possible distinction is needed to allow for using of various kind of
 * objects which contains an array of values that can be wrapped into a
 * vector.  For objects of type `object_t`, we can do somthing like:
 * <pre>
 *     object_t* obj = ...;
 *     opk_vspace_t* vspace = opk_new_simple_double_vector_space(get_number(obj));
 *     opk_vector_t* v = opk_wrap_simple_double_vector(vspace, get_data(obj),
 *                                                     (void*)obj,
 *                                                     delete_object);
 * </pre>
 * where `get_number()` returns the number of elements stored in the data part
 * of the object, `get_data()` returns the address of these elements, and
 * `delete_object()` delete the object.  Of course, if one prefers to keep the
 * control on the object management, passing `NULL` for the
 * `free_client_data()` function is always possible.
 *
 * @param vspace - The vector space which will own the vector.
 * @param data   - The array of values, must have at least `vspace->size`
 *                 elements.
 * @param client_data - Anything required by the `free_client_data()` method.
 * @param free_client_data - Function called to release ressources.  If not
 *                 `NULL`, it is called with argument `client_data` when the
 *                 vector is destroyed.
 *
 * @return A new vector of `vspace`, `NULL` in case of error.
 */
extern opk_vector_t*
opk_wrap_simple_double_vector(opk_vspace_t* vspace, double data[],
                              void* client_data,
                              void (*free_client_data)(void*));

extern double*
opk_get_simple_double_vector_data(opk_vector_t* v);

extern void*
opk_get_simple_double_vector_client_data(opk_vector_t* v);

extern opk_free_proc*
opk_get_simple_double_vector_free_client_data(opk_vector_t* v);

/**
 * Re-wrap an array into an existing simple vector.
 *
 * This functions replaces the contents of a simple wrapped vector.  It is
 * assumed that the vector `v` is a wrapped vector, that the new data
 * `new_data` is correctly aligned and large enough.  If the former
 * `free_client_data()` method of the wrapped vector `v` is not `NULL` and if
 * either the new `free_client_data()` method or the new `client_data` differ
 * from the former ones, then the former `free_client_data()` method is
 * applied to the former `client_data`.
 *
 * Re-wrapping is considered as a hack which merely saves the time needed to
 * allocate a container for a wrapped vector.  It is the caller responsibility
 * to ensure that all the assumptions hold.  In many cases deleting the old
 * vector and wrapping the arguments into a new vector is safer.
 *
 * @param v - The vector to re-wrap.
 * @param new_data - The new array of values.
 * @param new_client_data - The new client data.
 * @param new_free_client_data - The new method to free client data.
 *
 * @return `OPK_SUCCESS` or `OPK_FAILURE`.  In case of failure, global
 *         variable `errno` is set to: `EFAULT` if `v` or `new_data` are
 *         `NULL`, `EINVAL` if `v` is not a vector of the correct kind.
 */
extern int
opk_rewrap_simple_double_vector(opk_vector_t* v, double new_data[],
                                void* new_client_data,
                                void (*new_free_client_data)(void*));

/**
 * Create a vector space for array of float's in conventional memory.
 *
 * See `opk_new_simple_double_vector()` for a description.
 */
extern opk_vspace_t*
opk_new_simple_float_vector_space(opk_index_t size);

/**
 * Wrap an existing array into a simple vector.
 *
 * See `opk_wrap_simple_double_vector()` for a description.
 */
extern opk_vector_t*
opk_wrap_simple_float_vector(opk_vspace_t* vspace, float data[],
                             void* client_data,
                             void (*free_client_data)(void*));


extern float*
opk_get_simple_float_vector_data(opk_vector_t* v);

extern void*
opk_get_simple_float_vector_client_data(opk_vector_t* v);

extern opk_free_proc*
opk_get_simple_float_vector_free_client_data(opk_vector_t* v);

extern int
opk_rewrap_simple_float_vector(opk_vector_t* v, float new_data[],
                               void* new_client_data,
                               void (*new_free_client_data)(void*));


extern opk_vspace_t*
opk_allocate_vector_space(const void* ident,
                          opk_index_t nvariables,
                          opk_index_t nbytes);

extern void
opk_delete_vector_space(opk_vspace_t* vspace);

/**
 * Delete a vector.
 *
 * @param vect - The vector to delete (can be NULL).
 */
extern void
opk_vdelete(opk_vector_t* vect);

extern opk_vector_t*
opk_vcreate(opk_vspace_t* vspace);

extern void
opk_vzero(opk_vector_t* vect);

extern void
opk_vfill(opk_vector_t* vect, double alpha);

extern void
opk_vcopy(opk_vector_t* dst, const opk_vector_t* src);

extern void
opk_vscale(opk_vector_t* dst,
           double alpha, const opk_vector_t* src);

extern void
opk_vswap(opk_vector_t* x, opk_vector_t* y);

extern double
opk_vdot(const opk_vector_t* x, const opk_vector_t* y);

extern double
opk_vnorm2(const opk_vector_t* v);

extern double
opk_vnorm1(const opk_vector_t* v);

extern double
opk_vnorminf(const opk_vector_t* v);

extern void
opk_vaxpby(opk_vector_t* dst,
           double alpha, const opk_vector_t* x,
           double beta,  const opk_vector_t* y);

extern void
opk_vaxpbypcz(opk_vector_t* dst,
              double alpha, const opk_vector_t* x,
              double beta,  const opk_vector_t* y,
              double gamma, const opk_vector_t* z);

/*---------------------------------------------------------------------------*/
/* OPERATORS */

typedef struct _opk_operator opk_operator_t;

extern opk_operator_t*
opk_allocate_operator(opk_vspace_t* inpspace,
                      opk_vspace_t* outspace,
                      size_t size);

extern void
opk_delete_operator(opk_operator_t* op);

extern int
opk_apply_direct(opk_operator_t* op, opk_vector_t* dst,
                 const opk_vector_t* src);

extern int
opk_apply_adjoint(opk_operator_t* op, opk_vector_t* dst,
                  const opk_vector_t* src);

extern int
opk_apply_inverse(opk_operator_t* op, opk_vector_t* dst,
                  const opk_vector_t* src);

struct _opk_operator {
  opk_vspace_t* inpspace; /* The input space of the operator. */
  opk_vspace_t* outspace; /* The output space of the operator. */

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

  /* If not NULL, the delete() method is called by opk_delete_operator() to
     release any specific ressources.  The memory occupied by instance itself,
     which has been allocated by opk_allocate_operator(), is freed by
     opk_delete_operator(). */
  void (*delete)(opk_operator_t* self);
};

/*---------------------------------------------------------------------------*/
/* ERROR MANAGEMENT */

typedef void opk_error_handler(const char* message);

extern opk_error_handler*
opk_get_error_handler(void);

/**
 * Set the error handler.
 *
 * @param handler - The new error handler or NULL to restore the default
 *                  handler.
 */
extern opk_error_handler*
opk_set_error_handler(opk_error_handler* handler);

/**
 * Throw an error.
 *
 * This function calls the current error handler.  It is used in OptimPack
 * library to throw errors corresponding to a misuse of the library.  For
 * instance, when one attempts to compute the dot product of two vectors which
 * do not belong to the same vector space.
 *
 * @param resaon - The error message indicating the reason of the failure.
 */
extern void
opk_error(const char* reason);

/*---------------------------------------------------------------------------*/
/* LIMITED MEMORY VARIABLE METRIC METHODS */

extern void
opk_apply_lbfgs(opk_vspace_t* vspace,
                int m,
                const opk_vector_t* s[],
                const opk_vector_t* y[],
                const double rho[],
                double gamma,
                int mark,
                opk_vector_t* v,
                double alpha[]);

/*---------------------------------------------------------------------------*/
/* REVERSE COMMUNICATION */

/**
 * @brief Code returned by the reverse communication verison of optimzation
 * algorithms.
 */
typedef enum {
  OPK_TASK_ERROR      = -1, /**< An error has ocurred. */
  OPK_TASK_COMPUTE_FG =  0, /**< Caller shall compute f(x) and g(x). */
  OPK_TASK_NEW_X      =  1, /**< A new iterate is available. */
  OPK_TASK_FINAL_X    =  2, /**< Algorithm has converged, solution is available. */
  OPK_TASK_WARNING    =  3  /**< Algorithm terminated with a warning. */
} opk_task_t;


/*---------------------------------------------------------------------------*/
/* LINE SEARCH METHODS */

/* Workspace for line search methods. */
typedef struct _opk_lnsrch_workspace opk_lnsrch_workspace_t;

/* Create a workspace for Moré and Thuente cubic line search. */
extern opk_lnsrch_workspace_t*
opk_lnsrch_new_csrch(double ftol, double gtol, double xtol);

/* Create a workspace for Armijo (backtracking) line search. */
extern opk_lnsrch_workspace_t*
opk_lnsrch_new_backtrack(double ftol);

/* Possible values returned by opk_lnsrch_start and opk_lnsrch_iterate. */
#define OPK_LNSRCH_ERROR_ILLEGAL_ADDRESS                    -12
#define OPK_LNSRCH_ERROR_CORRUPTED_WORKSPACE                -11
#define OPK_LNSRCH_ERROR_BAD_WORKSPACE                      -10
#define OPK_LNSRCH_ERROR_STP_CHANGED                         -9
#define OPK_LNSRCH_ERROR_STP_OUTSIDE_BRACKET                 -8
#define OPK_LNSRCH_ERROR_NOT_A_DESCENT                       -7
#define OPK_LNSRCH_ERROR_STPMIN_GT_STPMAX                    -6
#define OPK_LNSRCH_ERROR_STPMIN_LT_ZERO                      -5
#define OPK_LNSRCH_ERROR_STP_LT_STPMIN                       -4
#define OPK_LNSRCH_ERROR_STP_GT_STPMAX                       -3
#define OPK_LNSRCH_ERROR_INITIAL_DERIVATIVE_GE_ZERO          -2
#define OPK_LNSRCH_ERROR_NOT_STARTED                         -1
#define OPK_LNSRCH_SEARCH                                     0
#define OPK_LNSRCH_CONVERGENCE                                1
#define OPK_LNSRCH_WARNING_ROUNDING_ERRORS_PREVENT_PROGRESS   2
#define OPK_LNSRCH_WARNING_XTOL_TEST_SATISFIED                3
#define OPK_LNSRCH_WARNING_STP_EQ_STPMAX                      4
#define OPK_LNSRCH_WARNING_STP_EQ_STPMIN                      5

/* Start a new line search.  The returned value is strictly negative to
   indicate an error; it is equal to zero otherwise. */
extern int
opk_lnsrch_start(opk_lnsrch_workspace_t* ws, double f0, double g0,
                 double stp, double stpmin, double stpmax);

/* Check whether line search has converged or update the step size.  The
   returned value is strictly negative to indicate an error; it is equal to
   zero when searching is in progress; it is strictly positive when line search
   has converged or cannot make any more progresses.*/
extern int
opk_lnsrch_iterate(opk_lnsrch_workspace_t* ws, double* stp_ptr,
                   double f1, double g1);

/* Destroy a line search workspace. */
extern void
opk_lnsrch_delete(opk_lnsrch_workspace_t* ws);

/* Get current step lenght. Returned value should be >= 0; -1 is returned in
   case of error. */
extern double
opk_lnsrch_get_step(const opk_lnsrch_workspace_t* ws);

extern int
opk_lnsrch_get_status(const opk_lnsrch_workspace_t* ws);

extern int
opk_lnsrch_has_errors(const opk_lnsrch_workspace_t* ws);

extern int
opk_lnsrch_has_warnings(const opk_lnsrch_workspace_t* ws);

extern int
opk_lnsrch_converged(const opk_lnsrch_workspace_t* ws);

extern int
opk_lnsrch_finished(const opk_lnsrch_workspace_t* ws);


/* Moré & Thuente method to perform a cubic safeguarded step. */
extern int
opk_cstep(double *stx_ptr, double *fx_ptr, double *dx_ptr,
          double *sty_ptr, double *fy_ptr, double *dy_ptr,
          double *stp_ptr, double  fp,     double  dp,
          int *brackt, double stpmin, double stpmax);

/*---------------------------------------------------------------------------*/
/* NON-LINEAR CONJUGATE GRADIENT METHODS */

typedef struct _opk_nlcg_workspace opk_nlcg_workspace_t;

/**
 * Create a new optimizer instance for non-linear conjugate gradient method.
 *
 * This function creates an optimizer instance for minimization by a
 * non-linear conjugate gradient method over a given vector space.  The
 * returned instance must be destroyed by calling opk_nlcg_delete when no
 * longer needed.
 *
 * @param vspace   The vector space of the unknowns.
 * @param method   A bitwise combination of the non-linear conjugate gradient
 *                 update method and options.
 *
 * @return The address of a new optimizer instance, or NULL in case of error.
 *         Global variable errno may be ENOMEM if there is not enough memory
 *         or EINVAL if one of the arguments is invalid or EFAULT if vspace is
 *         NULL.
 */
extern opk_nlcg_workspace_t*
opk_nlcg_new(opk_vspace_t* vspace, unsigned int method);

extern opk_nlcg_workspace_t*
opk_nlcg_new_with_line_search(opk_vspace_t* vspace, unsigned int method,
                              opk_lnsrch_workspace_t* lnsrch);
extern void
opk_nlcg_delete(opk_nlcg_workspace_t* ws);

extern opk_task_t
opk_nlcg_start(opk_nlcg_workspace_t* ws);

extern opk_task_t
opk_nlcg_iterate(opk_nlcg_workspace_t* ws, opk_vector_t* x1,
                 double f1, opk_vector_t* g1);

extern int
opk_nlcg_get_ftol(opk_nlcg_workspace_t* ws, double* frtol,
                  double* fatol);

extern int
opk_nlcg_get_gtol(opk_nlcg_workspace_t* ws, double* grtol,
                  double* gatol);

extern int
opk_nlcg_get_fmin(opk_nlcg_workspace_t* ws, double* fmin);

extern int
opk_nlcg_set_fmin(opk_nlcg_workspace_t* ws, double fmin);

extern int
opk_nlcg_unset_fmin(opk_nlcg_workspace_t* ws);

extern int
opk_nlcg_get_iterations(opk_nlcg_workspace_t* ws);

extern int
opk_nlcg_get_restarts(opk_nlcg_workspace_t* ws);

extern int
opk_nlcg_get_evaluations(opk_nlcg_workspace_t* ws);

extern unsigned int
opk_nlcg_get_method(opk_nlcg_workspace_t* ws);

extern opk_bool_t
opk_nlcg_get_starting(opk_nlcg_workspace_t* ws);

extern opk_task_t
opk_nlcg_get_task(opk_nlcg_workspace_t* ws);

extern double
opk_nlcg_get_alpha(opk_nlcg_workspace_t* ws);

extern double
opk_nlcg_get_beta(opk_nlcg_workspace_t* ws);


#define OPK_NLCG_FLETCHER_REEVES        1
#define OPK_NLCG_HESTENES_STIEFEL       2
#define OPK_NLCG_POLAK_RIBIERE_POLYAK   3
#define OPK_NLCG_FLETCHER               4
#define OPK_NLCG_LIU_STOREY             5
#define OPK_NLCG_DAI_YUAN               6
#define OPK_NLCG_PERRY_SHANNO           7
#define OPK_NLCG_HAGER_ZHANG            8
#define OPK_NLCG_POWELL              (1<<8) /* force beta >= 0 */
#define OPK_NLCG_SHANNO_PHUA         (1<<9) /* compute scale from previous
                                               iteration */

/* For instance: (OPK_NLCG_POLAK_RIBIERE_POLYAK | OPK_NLCG_POWELL) merely
   corresponds to PRP+ (Polak, Ribiere, Polyak) while (OPK_NLCG_PERRY_SHANNO |
   OPK_NLCG_SHANNO_PHUA) merely corresponds to the conjugate gradient method
   implemented in CONMIN. */

/* Default settings for non linear conjugate gradient (should correspond to
   the method which is, in general, the most successful). */
#define OPK_NLCG_DEFAULT (OPK_NLCG_HAGER_ZHANG | OPK_NLCG_SHANNO_PHUA)

/*---------------------------------------------------------------------------*/
/* MINIMIZATION OF AN UNIVARIATE FUNCTION */

/* Defines: Options for opk_fmin routines.
 *
 * OPK_FMIN_BOUNDED_BY_A  - point A of the initial interval is an
 *                          exlusive bound;
 * OPK_FMIN_BOUNDED_BY_B  - point B of the initial interval is an
 *                          exlusive bound;
 * OPK_FMIN_SMOOTH        - the function is smooth (will use Brent's
 *                          algorithm, otherwise will use golden section
 *                          algorithm).
 */
#define OPK_FMIN_BOUNDED_BY_A  1
#define OPK_FMIN_BOUNDED_BY_B  2
#define OPK_FMIN_SMOOTH        4

/* Values returned by the opk_fmin_next routine. */
#define OPK_FMIN_ERROR         (-1)
#define OPK_FMIN_START           0
#define OPK_FMIN_FX              1
#define OPK_FMIN_NEWX            2
#define OPK_FMIN_CONVERGENCE     3

/*
 * Function: opk_fmin
 *
 *   Search for the minimum of a function.
 *
 *
 * Description:
 *
 *   This function searches for the minimum _xmin_ of the univariate
 *   function f(_x_).
 *
 *   The algorithm requires an initial interval (_a_, _b_).  If the bit
 *   <OPK_FMIN_BOUNDED_BY_A> is set in _flags_, then the value _a_ is a strict
 *   bound for the search.  Similarly, if the bit <OPK_FMIN_BOUNDED_BY_B> is
 *   set in _flags_, then the value _b_ is a strict bound for the search.  If
 *   _a_ and/or _b_ are not exclusive bounds, their values are tried first by
 *   the algorithm.
 *
 *   If the bit <OPK_FMIN_SMOOTH> is set in _flags_, then the function is
 *   assumed to be smooth and Brent's algorithm [1] is used to find the
 *   minimum; otherwise, the golden section method is used.
 *
 *   The result is stored into array _out_ as follows
 *
 *    - out[0] = the approximative solution _xmin_
 *    - out[1] = the lower bound _xlo_of the final interval
 *    - out[2] = the upper bound _xhi_ of the final interval
 *    - out[3] = f(_xmin_)
 *    - out[4] = f(_xlo_)
 *    - out[5] = f(_xhi_)
 *    - out[6] = the number of function evaluations
 *
 *   Depending whether the input interval has exclusive bounds, the minimum
 *   number of function evaluations is between 1 and 3 whatever is the value of
 *   _maxeval_.  If the minimum has not been bracketted, the result is:
 *
 *    - out[0:2] = {_u_, _v_, _w_} and
 *    - out[3:5] = {f(_u_), f(_v_), f(_w_)}
 *
 *   with {_u_, _v_, _w_} the 3 last positions tried by the algorithm (in no
 *   particular order).
 *
 * Parameters:
 *   f       - The function to minimize.
 *   a       - The first point of the initial interval.
 *   b       - The other point of the initial interval.
 *   flags   - A bitwise combination of flags (see description).
 *   maxeval - If non-negative, the maximum number of function evaluations;
 *             otherwise, no limits.
 *   prec    - The relative precision: the algorithm stop when the uncertainty
 *             is less or equal _prec_ times the magnitude of the solution.
 *             If _prec_ is strictly negative, then a default precision of
 *             sqrt(_epsilon_) is used with _epsilon_ the machine relative
 *             precision.
 *   out     - A 7-element array to store the result.
 *
 * Returns:
 *    0 on convergence, 1 if too many iterations but minimum was bracketted, 2
 *    if too many iterations but minimum was *not* bracketted, and -1 on error
 *    (invalid input arguments).
 *
 * References:
 *    [1] Brent, R.P. 1973, "Algorithms for Minimization without Derivatives"
 *        (Englewood Cliffs, NJ: Prentice-Hall), Chapter 5.
 *
 * See Also:
 *    <opk_fmin_with_context>.
 */
extern int
opk_fmin(double (*f)(double x), double a, double b,
         unsigned int flags, long maxeval, double prec,
         double out[7]);


/*
 * Function: opk_fmin_with_context
 *
 *   Search for the minimum of a function.
 *
 * Description:
 *
 *   This function is identical to <opk_fmin> except that the user-defined
 *   function is called as f(_data_, _x_) to evaluate the funtion at _x_.
 *   The _data_ argument is simply passed to the user-defined function and may
 *   be used to store any parameters (but the variable value) needed by the
 *   function.
 *
 * Parameters:
 *   f         - The function to minimize.
 *   a         - The first point of the initial interval.
 *   b         - The other point of the initial interval.
 *   flags     - A bitwise combination of flags.
 *   maxeval   - If non-negative, the maximum number of function
 *               evaluations; otherwise no limits.
 *   prec      - The relative precision.
 *   out       - A 7-element array to store the result.
 *   data      - Anything needed by the user-defined function.
 *
 * Returns:
 *   Same value as <spf_fmin>.
 *
 * See Also:
 *   <opk_fmin>.
 */
extern int
opk_fmin_with_context(double (*f)(void *data, double x),
                      double a, double b, unsigned int flags,
                      long maxeval, double prec,
                      double out[7], void *data);

OPK_END_C_DECLS

#endif /* _OPTIMPACK_H */

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
