/*
 * optimpack.h --
 *
 * Common definitions for OptimPack library.
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

#ifndef _OPTIMPACK_H
#define _OPTIMPACK_H 1

#include <stddef.h>

/**
 * @defgroup Objects         Management of objects and derived types.
 * @defgroup Vectors         Vectors to store variables and vector spaces.
 * @defgroup Operators       Operators acting on vectors.
 * @defgroup LineSearch      Line search methods.
 * @defgroup NLCG            Non-linear conjugate gradient methods.
 * @defgroup VariableMetric  Variable metric methods.
 * @defgroup ConstrainedVariableMetric  Variable metric methods with convex constraints.
 * @defgroup Error           Reporting of errors.
 * @defgroup TrustRegion     Trust region methods.
 * @defgroup Utilities       Miscellaneous utility functions and macros.
 * @defgroup LowLevel        Low-level API to implement objects.
 * @defgroup LinearAlgebra   BLAS/LAPACK/LINPACK-like linear algebra routines for small arrays.
 *
 * @mainpage Optimization of large scale problems
 *
 * OptimPack is a library for large scale optimization problems.  Its primary
 * focus is solving large inverse problems as image reconstruction.
 *
 * @section ManagementOfObjects Management of objects and derived types.
 *
 * @subsection BasicObjects Basic objects.
 * @ref Objects
 *
 * @subsection VectorSpaces Vectors and vector spaces.
 * @ref Vectors
 *
 * @section Linear Algebra
 * @ref LinearAlgebra
 */

/*---------------------------------------------------------------------------*/

/*
 * Miscellaneous macros.
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


/**
 * Macro `OPK_JOIN(a,b)` joins its two arguments, possibly after performing
 * macro expansion of `a` and `b`.
 */
#define OPK_JOIN(a,b)                _OPK_JOIN2(a,b)

#define OPK_JOIN2(a,b)               _OPK_JOIN2(a,b)
#define OPK_JOIN3(a,b,c)             _OPK_JOIN3(a,b,c)
#define OPK_JOIN4(a,b,c,d)           _OPK_JOIN4(a,b,c,d)
#define OPK_JOIN5(a,b,c,d,e)         _OPK_JOIN5(a,b,c,d,e)
#define OPK_JOIN6(a,b,c,d,e,f)       _OPK_JOIN6(a,b,c,d,e,f)
#define OPK_JOIN7(a,b,c,d,e,f,g)     _OPK_JOIN7(a,b,c,d,e,f,g)
#define OPK_JOIN8(a,b,c,d,e,f,g,h)   _OPK_JOIN8(a,b,c,d,e,f,g,h)
#define OPK_JOIN9(a,b,c,d,e,f,g,h,i) _OPK_JOIN9(a,b,c,d,e,f,g,h,i)

/**
 * Macro `OPK_STRINGIFY(x)` wraps its argument in "" (double quotation marks),
 * possibly after performing macro expansion of the argument.
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

/**
 * Macro `OPK_ABS(a)` gives the absolute value of `a`.
 *
 * Beware that `a` is evaluated twice.
 */
#define OPK_ABS(a)   ((a) >= 0 ? (a) : -(a))

/**
 * Macro `OPK_MIN(a,b)` yields the minimum value between `a` and `b`.
 */
#define OPK_MIN(a,b) ((a) <= (b) ? (a) : (b))

/**
 * Macro `OPK_MAX(a,b)` yields the maximum value between `a` and `b`.
 */
#define OPK_MAX(a,b)  ((a) >= (b) ? (a) : (b))

/**
 * Macro `OPK_HOW_MANY(a,b)` yields the minimal number of chunks with `b`
 * elements needed to store `a` elements.  Both `a` and `b` must be integers.
 */
#define OPK_HOW_MANY(a,b)  ((((b) - 1) + (a))/(b))

/**
 * Macro `OPK_ROUND_UP(a,b)` yields the integer `a` rounded up to a multiple of
 * integer `b`.
 */
#define OPK_ROUND_UP(a,b)  (OPK_HOW_MANY(a,b)*(b))

/**
 * Macro `OPK_ADDRESS(type,base,offset)` yields a `type*` address at `offset`
 * (in bytes) from `base` address.
 */
#define OPK_ADDRESS(type, base, offset) ((type*)((char*)(base) + (offset)))

/**
 * Macro `OPK_OFFSET(type, field)` yields the offset (in bytes) of member
 * `field` in structure/class `type`.
 */
#ifdef offsetof
#  define OPK_OFFSET(type, field)  offsetof(type, field)
#elif 1
#  define OPK_OFFSET(type, field)  ((size_t)((char*)&((type*)0)->field))
#else
#  define OPK_OFFSET(type, field)  ((char*)&((type*)0)->field - (char*)0)
#endif

/**
 * Macro `OPK_LOOP(var,cnt)` yields a simple loop over variable `var` from 0 to
 * `cnt` - 1.
 */
#define OPK_LOOP(var, cnt)  for (var = 0; var < cnt; ++var)


/** Yields |a|*sign(b). */
#define opk_sign(a, b) ({ typeof(a) _1 = (a); typeof(b) _2 = (b); \
                          ((_1 < 0) == ((b) < 0)) ? _1 : -_1;})

/** Yield absolute value. */
#define opk_abs(a)    ({ typeof(a) _1 = (a); _1 >= 0 ? _1 : -_1; })

/** Yield minimum of two values. */
#define opk_min(a, b) ({ typeof(a) _1 = (a); typeof(b) _2 = (b); \
                         _1 <= _2 ? _1 : _2;})

/** Yield maximum of two values. */
#define opk_max(a, b) ({ typeof(a) _1 = (a); typeof(b) _2 = (b); \
                         _1 >= _2 ? _1 : _2;})
/*---------------------------------------------------------------------------*/

/**
 * Macro `OPK_NEW(type)` allocates an object of structure/class `type`.
 */
#define OPK_NEW(type)             ((type *)malloc(sizeof(type)))

/**
 * Macro `OPK_ALLOC(type,number)` allocates an array of `number` items of
 * structure/class `type`.
 */
#define OPK_ALLOC(type, number)   ((type *)malloc((number)*sizeof(type)))

/**
 * Macro `OPK_FREE(ptr)` frees pointer `ptr` if non-null.
 */
#define OPK_FREE(ptr)   if (! (ptr)) ; else free(ptr)

/*---------------------------------------------------------------------------*/

OPK_BEGIN_C_DECLS

/**
 *  Values returned by OptimPack routines.
 *
 *  OPK_FAILURE   - Status value returned upon success.
 *  OPK_SUCCESS   - Status value returned upon failure.
 */
typedef enum {
  OPK_FAILURE = -1,
  OPK_SUCCESS = 0
} opk_status_t;

/*---------------------------------------------------------------------------*/
/* DATA TYPES */

/**
 * Integer data type used for array indices in OptimPack.
 */
typedef ptrdiff_t opk_index_t;

/**
 * Data type for boolean (logical) values.
 */
typedef int opk_bool_t;

#define OPK_TRUE  1
#define OPK_FALSE 0

/*---------------------------------------------------------------------------*/
/* BASIC OBJECTS */

/**
 * @addtogroup Objects
 * @{
 */

/**
 * Opaque basic object type.
 *
 * Users of OptimPack library normally do not need to known exactly the
 * contents of an object.  They simply use the functions of the library to
 * manipulate object pointers.
 *
 * Developpers who need to implement derived types or specific variants of
 * existing types must include the header `optimpack-private.h` to unveil
 * the definitions of the structures representing the OptimPack objects.
 */
typedef struct _opk_object opk_object_t;

/**
 * Hold a reference on an object.
 *
 * This function increments the reference count of its argument and returns it.
 * If the argument is `NULL`, nothing is done.  Every call to this function must
 * be balanced by a call to `opk_drop_object()`.
 *
 * It is a good practice to hold a reference whenever a persistent structure
 * (e.g., another object) remembers the object.  To limit the number of casts,
 * the macro `OPK_HOLD()` can be used instead.
 *
 * @param obj - The object to lock (can be NULL).
 *
 * @return Its argument (with one more reference if not NULL).
 */
extern opk_object_t*
opk_hold_object(opk_object_t* obj);

/**
 * Drop a reference on an object.
 *
 * This function decrements the reference count of its argument and delete it
 * if there are no other references.  If the argument is `NULL`, nothing is
 * done.  To effectively delete an object, there must be a call to this
 * function for every call to `opk_hold_object()` on this object plus one (to
 * release the reference by the creator of the object).
 *
 * @param obj - The object to release (can be `NULL`).
 */
extern void
opk_drop_object(opk_object_t* obj);

/**
 * Get the number of references on an object.
 *
 * @param obj - The object (can be `NULL`).
 *
 * @return The number of references set on the object.  If the object address
 *         is `NULL`, the result is zero; otherwise, the result is greater of
 *         equal one.
 */
extern opk_index_t
opk_get_object_references(opk_object_t* obj);

/**
 * Macro to set a reference on an object whatever its type.  Compared to the
 * function `opk_hold_object()`, this macro casts its argument to be an
 * `opk_object_t` pointer.  The result is the argument cast as a basic object
 * pointer.
 */
#define OPK_HOLD(obj)       opk_hold_object((opk_object_t*)(obj))

/**
 * Macro to drop an object whatever its type.  Compared to the function
 * `opk_drop_object()`, this macro casts its argument to be an `opk_object_t`
 * pointer.
 */
#define OPK_DROP(obj)       opk_drop_object((opk_object_t*)(obj))

/**
 * Macro to query the number of references on an object whatever its type.
 * Compared to the function `opk_get_object_references()`, this macro casts its
 * argument to be a `opk_object_t` pointer.
 */
#define OPK_REFS(obj)       opk_get_object_references((opk_object_t*)(obj))

#define OPK_HOLD_VSPACE(vsp)  ((opk_vspace_t*)OPK_HOLD(vsp))
#define OPK_HOLD_VECTOR(vec)  ((opk_vector_t*)OPK_HOLD(vec))
#define OPK_HOLD_LNSRCH(vec)  ((opk_lnsrch_t*)OPK_HOLD(vec))
#define OPK_HOLD_OPERATOR(op) ((opk_operator_t*)OPK_HOLD(op))

/** @} */

/*---------------------------------------------------------------------------*/
/* VECTORS AND VECTOR SPACES */

/**
 * @addtogroup Vectors
 * @{
 */

/** Opaque vector type.  This sub-type inherits from `opk_object_t`. */
typedef struct _opk_vector opk_vector_t;

/** Opaque vector space type.  This sub-type inherits from `opk_object_t`. */
typedef struct _opk_vspace opk_vspace_t;

/** Prototype of function to release client data. */
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
 * ~~~~~~~~~~{.c}
 * #define N 1000
 * opk_vspace_t* vspace = opk_new_simple_double_vector_space(N);
 *
 * double heap_array[N];
 * opk_vector_t* v1 = opk_wrap_simple_double_vector(vspace, heap_array,
 *                                                  NULL, NULL);
 *
 * double* dynamic_array = (double*)malloc(N*sizeof(double));
 * opk_vector_t* v2 = opk_wrap_simple_double_vector(vspace, dynamic_array,
 *                                                  dynamic_array, free);
 * ~~~~~~~~~~
 *
 * which creates two vectors, `v1` and `v2`, which are respectively wrapped
 * around an array allocated on the heap and around a dynamically allocated
 * array.
 *
 * In the above example, the `client_data` and the `data` are the same but the
 * possible distinction is needed to allow for using of various kind of
 * objects which contain an array of values that can be wrapped into a
 * vector.  For objects of type `object_t`, we can do something like:
 * ~~~~~~~~~~{.c}
 * object_t* obj = ...;
 * opk_vspace_t* vspace = opk_new_simple_double_vector_space(get_number(obj));
 * opk_vector_t* v = opk_wrap_simple_double_vector(vspace, get_data(obj),
 *                                                 (void*)obj,
 *                                                 delete_object);
 * ~~~~~~~~~~
 * where `get_number()` returns the number of elements stored in the data part
 * of the object, `get_data()` returns the address of these elements, and
 * `delete_object()` delete the object.  Of course, if one prefers to keep the
 * control on the object management, passing `NULL` for the
 * `free_client_data()` function is always possible.
 *
 * @param vspace - The vector space which will own the vector.
 * @param data   - The array of values, must have at least `vspace->size`
 *                 elements.
 * @param free_client_data - Function called to release ressources.  If not
 *                 `NULL`, it is called with argument `client_data` when the
 *                 vector is destroyed.
 * @param client_data - Anything required by the `free_client_data()` method.
 *
 * @return A new vector of `vspace`, `NULL` in case of error.
 */
extern opk_vector_t*
opk_wrap_simple_double_vector(opk_vspace_t* vspace, double data[],
                              void (*free_client_data)(void*),
                              void* client_data);

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
 * assumed that the vector `vect` is a wrapped vector, that the new data
 * `new_data` is correctly aligned and large enough.  If the former
 * `free_client_data()` method of the wrapped vector `vect` is not `NULL` and
 * if either the new `free_client_data()` method or the new `client_data`
 * differ from the former ones, then the former `free_client_data()` method is
 * applied to the former `client_data`.
 *
 * Re-wrapping is considered as a hack which merely saves the time needed to
 * allocate a container for a wrapped vector.  It is the caller responsibility
 * to ensure that all the assumptions hold.  In many cases dropping the old
 * vector and wrapping the arguments into a new vector is safer.
 *
 * @param vect - The vector to re-wrap.
 * @param new_data - The new array of values.
 * @param new_free_client_data - The new method to free client data.
 * @param new_client_data - The new client data.
 *
 * @return `OPK_SUCCESS` or `OPK_FAILURE`.  In case of failure, global variable
 *         `errno` is set to: `EFAULT` if `vect` or `new_data` are `NULL`,
 *         `EINVAL` if `vect` is not a vector of the correct kind.
 */
extern int
opk_rewrap_simple_double_vector(opk_vector_t* vect, double new_data[],
                                void (*new_free_client_data)(void*),
                                void* new_client_data);

/**
 * Create a vector space for array of float's in conventional memory.
 *
 * See `opk_new_simple_double_vector_space()` for a description.
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
                             void (*free_client_data)(void*),
                             void* client_data);


extern float*
opk_get_simple_float_vector_data(opk_vector_t* v);

extern void*
opk_get_simple_float_vector_client_data(opk_vector_t* v);

extern opk_free_proc*
opk_get_simple_float_vector_free_client_data(opk_vector_t* v);

extern int
opk_rewrap_simple_float_vector(opk_vector_t* v, float new_data[],
                               void (*new_free_client_data)(void*),
                               void* new_client_data);

/**
 * Create a vector instance.
 *
 * This functions creates a new vector of a given vector space.  The contents
 * of the vector is undefined.  The caller holds a reference on the returned
 * vector which has to be released with the function `opk_drop_object()` or
 * with the macro `OPK_DROP()`.
 *
 * @param vspace - The vector space which owns the vector to create.
 *
 * @return A new vector of the vector space; `NULL` in case of error.
 */
extern opk_vector_t*
opk_vcreate(opk_vspace_t* vspace);

/**
 * Fill a vector with zeros.
 *
 * This functions set all elements of a vector to zero.
 *
 * @param vect - The vector to fill.
 */
extern void
opk_vzero(opk_vector_t* vect);

/**
 * Fill a vector with a given value.
 *
 * This functions set all elements of a vector to the given value.
 *
 * @param vect  - The vector to fill.
 * @param alpha - The value.
 */
extern void
opk_vfill(opk_vector_t* vect, double alpha);

/**
 * Copy vector contents.
 *
 * This functions copies the contents of the source vector into the destination
 * vector.  Both vectors must belong to the same vector space.
 *
 * @param dst - The destination vector.
 * @param src - The source vector.
 */
extern void
opk_vcopy(opk_vector_t* dst, const opk_vector_t* src);

/**
 * Scale a vector by a scalar.
 *
 * This functions multiplies the elements of the source vector by the given
 * value and stores the result into the destination vector.  Both vectors must
 * belong to the same vector space.  The operation is optimized for specfific
 * values of `alpha`: with `alpha = 1`, the operation is the same as
 * `opk_vcopy()`; with `alpha = 0`, the operation is the same as `opk_vzero()`.
 *
 * @param dst   - The destination vector.
 * @param alpha - The scale factor.
 * @param src   - The source vector.
 */
extern void
opk_vscale(opk_vector_t* dst,
           double alpha, const opk_vector_t* src);

/**
 * Exchange vector contents.
 *
 * This functions exchanges the contents of two vectors.  Both vectors must
 * belong to the same vector space.
 *
 * @param x - A vector.
 * @param y - Another vector.
 */
extern void
opk_vswap(opk_vector_t* x, opk_vector_t* y);

/**
 * Compute the inner product of two vectors.
 *
 * This functions computes the inner product, also known as scalar product, of
 * two vectors.  Both vectors must belong to the same vector space.
 *
 * @param x - A vector.
 * @param y - Another vector.
 *
 * @return The inner product of the two vectors, that is the sum of the product
 *         of their elements.
 */
extern double
opk_vdot(const opk_vector_t* x, const opk_vector_t* y);

/**
 * Compute the inner product of three vectors.
 *
 * This functions computes the sum of the componentwise product of the elements
 * of three vectors.  All three vectors must belong to the same vector space.
 *
 * @param w - A vector.
 * @param x - Another vector.
 * @param y - Yet another vector.
 *
 * @return The sum of the componentwise product of the elements of three
 *         vectors.
 */
extern double
opk_vdot3(const opk_vector_t* w, const opk_vector_t* x, const opk_vector_t* y);

/**
 * Compute the L2 norm of a vector.
 *
 * This functions computes the L2 norm, also known as the Euclidean norm, of a
 * vector.
 *
 * @param v - A vector.
 *
 * @return The Euclidean norm of the vector, that is the square root of the sum
 *         of its squared elements.
 */
extern double
opk_vnorm2(const opk_vector_t* v);

/**
 * Compute the L1 norm of a vector.
 *
 * This functions computes the L1 norm of a vector.
 *
 * @param v - A vector.
 *
 * @return The L1 norm of the vector, that is the sum of the absolute values of
 *         its elements.
 */
extern double
opk_vnorm1(const opk_vector_t* v);

/**
 * Compute the infinite norm of a vector.
 *
 * This functions computes the infinite norm of a vector.
 *
 * @param v - A vector.
 *
 * @return The infinite norm of the vector, that is the maximum absolute value
 *         of its elements.
 */
extern double
opk_vnorminf(const opk_vector_t* v);

/**
 * Compute the elementwise product of two vectors.
 *
 * All vectors must belong to the same vector space.
 *
 * @param dst   - The destination vector.
 * @param x     - A vector.
 * @param y     - Another vector.
 */
extern void
opk_vproduct(opk_vector_t* dst,
             const opk_vector_t* x,
             const opk_vector_t* y);

/**
 * Compute the linear combination of two vectors.
 *
 * This functions stores in the destination vector `dst` the linear combination
 * `alpha*x + beta*y` where `alpha` and `beta` are two scalars while `x` and
 * `y` are two vectors.  All vectors must belong to the same vector space.
 *
 * @param dst   - The destination vector.
 * @param alpha - The factor for the vector `x`.
 * @param x     - A vector.
 * @param beta  - The factor for the vector `y`.
 * @param y     - Another vector.
 */
extern void
opk_vaxpby(opk_vector_t* dst,
           double alpha, const opk_vector_t* x,
           double beta,  const opk_vector_t* y);

/**
 * Compute the linear combination of three vectors.
 *
 * This functions stores in the destination vector `dst` the linear combination
 * `alpha*x + beta*y + gamma*z` where `alpha`, `beta` and `gamma` are three
 * scalars while `x`, `y` and `z` are three vectors.  All vectors must belong
 * to the same vector space.
 *
 * @param dst   - The destination vector.
 * @param alpha - The factor for the vector `x`.
 * @param x     - A vector.
 * @param beta  - The factor for the vector `y`.
 * @param y     - Another vector.
 * @param gamma - The factor for the vector `z`.
 * @param z     - Yet another vector.
 */
extern void
opk_vaxpbypcz(opk_vector_t* dst,
              double alpha, const opk_vector_t* x,
              double beta,  const opk_vector_t* y,
              double gamma, const opk_vector_t* z);

/** @} */

/*---------------------------------------------------------------------------*/
/* OPERATORS */

/**
 * @addtogroup Operators
 * @{
 */

/** Opaque operator type.  This sub-type inherits from `opk_object_t`. */
typedef struct _opk_operator opk_operator_t;

extern int
opk_apply_direct(opk_operator_t* op, opk_vector_t* dst,
                 const opk_vector_t* src);

extern int
opk_apply_adjoint(opk_operator_t* op, opk_vector_t* dst,
                  const opk_vector_t* src);

extern int
opk_apply_inverse(opk_operator_t* op, opk_vector_t* dst,
                  const opk_vector_t* src);

/** @} */

/*---------------------------------------------------------------------------*/
/* ERROR MANAGEMENT */

/**
 * @addtogroup Error
 * @{
 */

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
 * @param reason - The error message indicating the reason of the failure.
 */
extern void
opk_error(const char* reason);

/** @} */

/*---------------------------------------------------------------------------*/
/* REVERSE COMMUNICATION */

/**
 * @brief Code returned by the reverse communication version of optimzation
 * algorithms.
 */
typedef enum {
  OPK_TASK_ERROR       = -1, /**< An error has ocurred. */
  OPK_TASK_START       =  0, /**< Caller must call `start` method. */
  OPK_TASK_COMPUTE_FG  =  1, /**< Caller must compute f(x) and g(x). */
  OPK_TASK_NEW_X       =  2, /**< A new iterate is available. */
  OPK_TASK_FINAL_X     =  3, /**< Algorithm has converged, solution is
                                  available. */
  OPK_TASK_WARNING     =  4  /**< Algorithm terminated with a warning. */
} opk_task_t;


/*---------------------------------------------------------------------------*/
/* LINE SEARCH METHODS */

/**
 * @addtogroup LineSearch
 * @{
 */

/** Opaque linesearch type.  This sub-type inherits from `opk_object_t`. */
typedef struct _opk_lnsrch opk_lnsrch_t;

/** Create a Moré and Thuente cubic line search. */
extern opk_lnsrch_t*
opk_lnsrch_new_csrch(double ftol, double gtol, double xtol);

/** Create a backtracking (Armijo) line search. */
extern opk_lnsrch_t*
opk_lnsrch_new_backtrack(double ftol);

/**
 * Create a nonmonotone line search.
 *
 * Nonmonotone line search is described in SPG2 paper:
 *
 * E.G. Birgin, J.M. Martínez, & M. Raydan, "<i>Nonmonotone
 * Spectral Projected Gradient Methods on Convex Sets</i>," SIAM
 * J. Optim. <b>10</b>, 1196-1211 (2000).
 *
 * The parameters used in the SPG2 paper:
 * <pre>
 * m = 10
 * ftol = 1E-4
 * sigma1 = 0.1
 * sigma2 = 0.9
 * </pre>
 *
 * With {@code m = 1}, this line search method is equivalent to Armijo's line
 * search except that it attempts to use quadratic interpolation rather than
 * systematically use bisection to reduce the step size.
 *
 * @param m      - Number of previous steps to remember.
 * @param ftol   - Parameter for the function reduction criterion.
 * @param sigma1 - Lower steplength bound to trigger bissection.
 * @param sigma2 - Upper steplength relative bound to trigger bissection.
 *
 * @return A new line search object.
 */
extern opk_lnsrch_t*
opk_lnsrch_new_nonmonotone(opk_index_t m, double ftol,
                           double sigma1, double sigma2);

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

/**
 * Get the description of a line search status.
 *
 * @param status - The status code (e.g., as returned by
 *                 `opk_lnsrch_get_status()`).
 *
 * @return A pointer to a string describing the status or `NULL` if
 *         the status does not correspond to any known status.
 */
extern const char*
opk_lnsrch_message(int status);

/**
 * Start a new linesearch.
 *
 * @param ls     - The linesearch object.
 * @param f0     - The function value at the start of the line search (that
 *                 is, for a step of length 0).
 * @param df0    - The dirctional derivative at the start of the line search.
 * @param stp1   - The length of the first step to try (must be between
 *                 `stpmin` and `stpmax`).
 * @param stpmin - The minimum allowed step length (must be nonnegative).
 * @param stpmax - The maximum allowed step length (must be greater than
 *                 `stpmin`).
 *
 * @return The linesearch task, which is normally `OPK_LNSRCH_SEARCH` (zero).
 *         A different value (strictly negative) indicate an error.
 */
extern int
opk_lnsrch_start(opk_lnsrch_t* ls, double f0, double df0,
                 double stp1, double stpmin, double stpmax);

/**
 * Check whether line search has converged or update the step size.
 *
 * @param ls      - The linesearch object.
 * @param stp_ptr - The address at which the step length is stored.  On entry,
 *                  it must be the current step length; on exit, it is updated
 *                  with the next step to try (unless the line search is
 *                  finished).
 * @param f       - The function value for the current step length.
 * @param df      - The directional derivative for the current step length.

 * @return The returned value is strictly negative to indicate an error; it is
 * equal to zero when searching is in progress; it is strictly positive when
 * line search has converged or cannot make any more progresses.
 */
extern int
opk_lnsrch_iterate(opk_lnsrch_t* ls, double* stp_ptr,
                   double f, double df);

/** Get current step lenght. Returned value should be >= 0; -1 is returned in
   case of error. */
extern double
opk_lnsrch_get_step(const opk_lnsrch_t* ls);

extern int
opk_lnsrch_get_status(const opk_lnsrch_t* ls);

extern opk_bool_t
opk_lnsrch_has_errors(const opk_lnsrch_t* ls);

extern opk_bool_t
opk_lnsrch_has_warnings(const opk_lnsrch_t* ls);

extern opk_bool_t
opk_lnsrch_converged(const opk_lnsrch_t* ls);

extern opk_bool_t
opk_lnsrch_finished(const opk_lnsrch_t* ls);

extern opk_bool_t
opk_lnsrch_use_deriv(const opk_lnsrch_t* ls);

/** Moré & Thuente method to perform a cubic safeguarded step. */
extern int
opk_cstep(double *stx_ptr, double *fx_ptr, double *dx_ptr,
          double *sty_ptr, double *fy_ptr, double *dy_ptr,
          double *stp_ptr, double  fp,     double  dp,
          int *brackt, double stpmin, double stpmax);

/** @} */

/*---------------------------------------------------------------------------*/
/* NON-LINEAR CONJUGATE GRADIENT METHODS */

/**
 * @addtogroup NLCG
 * @{
 */

/** Opaque type for non-linear conjugate gradient optimizers. */
typedef struct _opk_nlcg opk_nlcg_t;

/**
 * Create a new optimizer instance for non-linear conjugate gradient method.
 *
 * This function creates an optimizer instance for minimization by a non-linear
 * conjugate gradient method over a given vector space.  The returned instance
 * must be unreferenced by calling the function `opk_drop_object()`, or th
 * macro `OPK_DROP()` when no longer needed.
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
extern opk_nlcg_t*
opk_new_nlcg_optimizer(opk_vspace_t* vspace, unsigned int method);

/*
 * Note that the optimizer will hold a reference to the line search object.
 */
extern opk_nlcg_t*
opk_new_nlcg_optimizer_with_line_search(opk_vspace_t* vspace,
                                        unsigned int method,
                                        opk_lnsrch_t* lnsrch);
extern opk_task_t
opk_start_nlcg(opk_nlcg_t* opt);

extern opk_task_t
opk_iterate_nlcg(opk_nlcg_t* opt, opk_vector_t* x,
                 double f, opk_vector_t* g);

/**
 * Get the absolute threshold for the norm or the gradient for convergence.
 *
 * This function retrieves the absolute threshold for the norm or the gradient
 * for convergence.  The convergence of the non-linear convergence gradient
 * (NLCG) method is defined by:
 * <pre>
 * ||g|| <= max(0, gatol, grtol*||ginit||)
 * </pre>
 * where `||g||` is the Euclidean norm of the current gradient `g`, `||ginit||`
 * is the Euclidean norm of the initial gradient `ginit` while `gatol` and
 * `grtol` are absolute and relative thresholds.

 * @param opt - The NLCG optimizer or `NULL` to get the default value.
 *
 * @return The value of `gatol` for the optimzer `opt`, if not `NULL`; the
 *         default value of `gatol`, if `opt` is `NULL`.
 */
extern double
opk_get_nlcg_gatol(opk_nlcg_t* opt);

/**
 * Set the absolute threshold for the norm or the gradient for convergence.
 *
 * @param opt - The NLCG optimizer.
 * @param gatol - The new value of the absolute threshold.
 *
 * @return `OPK_SUCCESS`, or `OPK_FAILURE` on error with global variable
 * `errno` set to `EFAULT` if `opt` is `NULL` and to `EINVAL` if the value of
 * `gatol` is not valid.
 *
 * @see opk_get_nlcg_gatol()
 */
extern int
opk_set_nlcg_gatol(opk_nlcg_t* opt, double gatol);

/**
 * Get the relative threshold for the norm or the gradient for convergence.
 *
 * @param opt - The NLCG optimizer or `NULL` to get the default value.
 *
 * @return The value of `grtol` for the optimzer `opt`, if not `NULL`; the
 *         default value of `grtol`, if `opt` is `NULL`.
 *
 * @see opk_get_nlcg_gatol()
 */
extern double
opk_get_nlcg_grtol(opk_nlcg_t* opt);

/**
 * Set the relative threshold for the norm or the gradient for convergence.
 *
 * @param opt - The NLCG optimizer.
 * @param grtol - The new value of the relative threshold.
 *
 * @return `OPK_SUCCESS`, or `OPK_FAILURE` on error with global variable
 * `errno` set to `EFAULT` if `opt` is `NULL` and to `EINVAL` if the value of
 * `grtol` is not valid.
 *
 * @see opk_get_nlcg_grtol()
 */
extern int
opk_set_nlcg_grtol(opk_nlcg_t* opt, double grtol);

/**
 * Get the current step length.
 *
 * This function retrieves the value of the current step size.
 *
 * @param opt - The NLCG optimizer.
 *
 * @return The value of the current step size, should be strictly positive.
 *
 * @see opk_get_nlcg_stpmax(), opk_set_nlcg_stpmin_and_stpmax().
 */
extern double
opk_get_nlcg_step(opk_nlcg_t* opt);

/**
 * Get the minimum relative step size.
 *
 * This function retrieves the value of the minimum relative step size `stpmin`
 * for the line search.
 *
 * @param opt - The NLCG optimizer or `NULL` to get the default value.
 *
 * @return The value of `stpmin` for the optimzer `opt`, if not `NULL`; the
 *         default value of `stpmin`, if `opt` is `NULL`.
 *
 * @see opk_get_nlcg_stpmax(), opk_set_nlcg_stpmin_and_stpmax().
 */
extern double
opk_get_nlcg_stpmin(opk_nlcg_t* opt);

/**
 * Get the maximum relative step size.
 *
 * During a line search, the step is constrained to be within `stpmin` and
 * `stpmax` times the lenght of the first step.  The relative bounds must be
 * such that:
 * <pre>
 * 0 <= stpmin < stpmax
 * </pre>
 * This function retrieves the value of the maximum relative step size `stpmax`
 * for the line search.
 *
 * @param opt - The NLCG optimizer or `NULL` to get the default value.
 *
 * @return The value of `stpmax` for the optimzer `opt`, if not `NULL`; the
 *         default value of `stpmax`, if `opt` is `NULL`.
 *
 * @see opk_get_nlcg_stpmin(), opk_set_nlcg_stpmin_and_stpmax().
 */
extern double
opk_get_nlcg_stpmax(opk_nlcg_t* opt);

/**
 * Set the relative step limits.
 *
 * @param opt - The NLCG optimizer or `NULL` to get the default value.
 * @param stpmin - The minimum relative step size.
 * @param stpmax - The maximum relative step size.
 *
 * @return `OPK_SUCCESS`, or `OPK_FAILURE` on error with global variable
 * `errno` set to `EFAULT` if `opt` is `NULL` and to `EINVAL` if the values of
 * `stpmin` or `stpmax` are not valid.
 *
 * @see opk_get_nlcg_stpmin(), opk_get_nlcg_stpmax().
 */
extern int
opk_set_nlcg_stpmin_and_stpmax(opk_nlcg_t* opt, double stpmin, double stpmax);

extern int
opk_get_nlcg_fmin(opk_nlcg_t* opt, double* fmin);

extern int
opk_set_nlcg_fmin(opk_nlcg_t* opt, double fmin);

extern int
opk_unset_nlcg_fmin(opk_nlcg_t* opt);

extern opk_index_t
opk_get_nlcg_iterations(opk_nlcg_t* opt);

extern opk_index_t
opk_get_nlcg_restarts(opk_nlcg_t* opt);

extern opk_index_t
opk_get_nlcg_evaluations(opk_nlcg_t* opt);

extern unsigned int
opk_get_nlcg_method(opk_nlcg_t* opt);

extern opk_bool_t
opk_get_nlcg_starting(opk_nlcg_t* opt);

extern opk_task_t
opk_get_nlcg_task(opk_nlcg_t* opt);

extern double
opk_get_nlcg_beta(opk_nlcg_t* opt);

/* Rules to compute the search direction in NLCG: */
#define OPK_NLCG_FLETCHER_REEVES        1
#define OPK_NLCG_HESTENES_STIEFEL       2
#define OPK_NLCG_POLAK_RIBIERE_POLYAK   3
#define OPK_NLCG_FLETCHER               4
#define OPK_NLCG_LIU_STOREY             5
#define OPK_NLCG_DAI_YUAN               6
#define OPK_NLCG_PERRY_SHANNO           7
#define OPK_NLCG_HAGER_ZHANG            8

/* The rule can be combined (bitwise or'ed) with the following bits to force
   beta >= 0 (according to Powell's prescription): */
#define OPK_NLCG_POWELL              (1<<8)

/* The rule can be combined (bitwise or'ed) with one of the following bits to
   compute the initial step size from the previous iteration: */
#define OPK_NLCG_SHANNO_PHUA         (1<<9)
#define OPK_NLCG_OREN_SPEDICATO      (2<<9)
#define OPK_NLCG_BARZILAI_BORWEIN    (3<<9)

/* For instance: (OPK_NLCG_POLAK_RIBIERE_POLYAK | OPK_NLCG_POWELL) merely
   corresponds to PRP+ (Polak, Ribiere, Polyak) while (OPK_NLCG_PERRY_SHANNO |
   OPK_NLCG_SHANNO_PHUA) merely corresponds to the conjugate gradient method
   implemented in CONMIN. */

/* Default settings for non linear conjugate gradient (should correspond to
   the method which is, in general, the most successful). */
#define OPK_NLCG_DEFAULT (OPK_NLCG_HAGER_ZHANG | OPK_NLCG_SHANNO_PHUA)


/** @} */

/*---------------------------------------------------------------------------*/
/* LIMITED MEMORY BFGS OPERATOR */

/**
 * @defgroup LBGFS      Limited memory BFGS operator
 * @ingroup VariableMetric
 * @{
 *
 * Implement limited memory quasi-Newton approximation of the inverse Hessian.
 *
 * The approximation of the inverse of Hessian is based on BFGS (Broyden,
 * Fletcher, Goldfarb & Shanno) updates using the 2-loop recursive algorithm of
 * Strang (described in Nocedal, 1980) combined with a preconditioner (initial
 * approximation of the inverse Hessian) or automatic scalings (along the ideas
 * of Gilbert & Lemaréchal (1989); and Shanno).
 */

/** Opaque type to store a limited memory BFGS operator.  This type inherits
    from the `opk_operator_t` object.  */
typedef struct _opk_lbfgs_operator opk_lbfgs_operator_t;

/** Rules for scaling the inverse Hessian approximation. */
#define OPK_SCALING_NONE               0
#define OPK_SCALING_OREN_SPEDICATO     1  /**< gamma = <s,y>/<y,y> */
#define OPK_SCALING_BARZILAI_BORWEIN   2  /**< gamma = <s,s>/<s,y> */

/**
 * Create a new limited memory BFGS operator.
 *
 * @param vspace - The input and output vector space of the operator.
 * @param m       - The maximum number of previous steps to memorize.
 * @param scaling - The rule for scaling the approximation of the
 *                   inverse Hessian.
 */
extern opk_lbfgs_operator_t*
opk_new_lbfgs_operator(opk_vspace_t* vspace, opk_index_t m,
                       int rule);

/** Forget all memorized steps in limited memory BFGS operator. */
extern void
opk_reset_lbfgs_operator(opk_lbfgs_operator_t* op);

/**
 * Query a memorized variable difference from a limited memory BFGS operator.
 *
 * It is the caller responsibility to use proper arguments.
 *
 * @param op - A limited memory BFGS operator.
 * @param k - The index of the memorized step to consider.
 *
 * @return s_k
 */
extern opk_vector_t*
opk_get_lbfgs_s(opk_lbfgs_operator_t* op, opk_index_t k);

/**
 * Query a memorized gradient difference from a limited memory BFGS operator.
 *
 * It is the caller responsibility to use proper arguments.
 *
 * @param op - A limited memory BFGS operator.
 * @param k - The index of the memorized step to consider.
 *
 * @return y_k
 */
extern opk_vector_t*
opk_get_lbfgs_y(opk_lbfgs_operator_t* op, opk_index_t k);

extern void
opk_set_lbfgs_operator_preconditioner(opk_lbfgs_operator_t* op,
                                      opk_operator_t* B0);

/**
 * Update LBFGS operator with a new pair of variables and gradient
 * differences.
 * @param x1 - The new variables.
 * @param x0 - The previous variables.
 * @param g1 - The gradient at {@code x1}.
 * @param g0 - The gradient at {@code x0}.
 * @throws IncorrectSpaceException
 */
extern void
opk_update_lbfgs_operator(opk_lbfgs_operator_t* op,
                          const opk_vector_t* x1,
                          const opk_vector_t* x0,
                          const opk_vector_t* g1,
                          const opk_vector_t* g0);

extern void
opk_set_lbfgs_operator_preconditioner(opk_lbfgs_operator_t* op,
                                      opk_operator_t* B0);

/** @} */

/*---------------------------------------------------------------------------*/
/* VARIABLE METRIC OPTIMIZATION METHOD */

/**
 * @addtogroup VariableMetric
 * @{
 */

/** Opaque type for a variable metric optimizer. */
typedef struct _opk_vmlm opk_vmlm_t;

extern opk_vmlm_t*
opk_new_vmlm_optimizer(opk_vspace_t* vspace, opk_index_t m);

extern opk_vmlm_t*
opk_new_vmlm_optimizer_with_line_search(opk_vspace_t* vspace,
                                        opk_index_t m,
                                        opk_lnsrch_t* lnsrch);

extern opk_task_t
opk_start_vmlm(opk_vmlm_t* opt);

extern opk_task_t
opk_iterate_vmlm(opk_vmlm_t* opt, opk_vector_t* x,
                 double f, opk_vector_t* g);

extern opk_task_t
opk_get_vmlm_task(opk_vmlm_t* opt);

extern opk_index_t
opk_get_vmlm_iterations(opk_vmlm_t* opt);

extern opk_index_t
opk_get_vmlm_evaluations(opk_vmlm_t* opt);

extern opk_index_t
opk_get_vmlm_restarts(opk_vmlm_t* opt);

extern int
opk_get_vmlm_scaling(opk_vmlm_t* opt);

extern int
opk_set_vmlm_scaling(opk_vmlm_t* opt, int scaling);

extern double
opk_get_vmlm_gatol(opk_vmlm_t* opt);

extern int
opk_set_vmlm_gatol(opk_vmlm_t* opt, double gatol);

extern double
opk_get_vmlm_grtol(opk_vmlm_t* opt);

extern int
opk_set_vmlm_grtol(opk_vmlm_t* opt, double grtol);

extern double
opk_get_vmlm_step(opk_vmlm_t* opt);

extern double
opk_get_vmlm_stpmin(opk_vmlm_t* opt);

extern double
opk_get_vmlm_stpmax(opk_vmlm_t* opt);

extern int
opk_set_vmlm_stpmin_and_stpmax(opk_vmlm_t* opt, double stpmin, double stpmax);

/** @} */

/*---------------------------------------------------------------------------*/
/* SEPARABLE BOUND CONSTRAINTS */

/**
 * @addtogroup BoxConstraints
 * @{
 */

typedef struct _opk_bound opk_bound_t;

typedef enum {
  OPK_BOUND_NONE = 0,
  OPK_BOUND_SCALAR = 1,
  OPK_BOUND_VECTOR = 2
} opk_bound_type_t;

extern opk_bound_t*
opk_new_bound(opk_vspace_t* space, opk_bound_type_t type,
              void* value);
extern void
opk_unset_bound(opk_bound_t* bnd);

extern void
opk_set_scalar_bound(opk_bound_t* bnd, double val);

extern int
opk_set_vector_bound(opk_bound_t* bnd, opk_vector_t* vec);

extern opk_bound_type_t
opk_get_bound_type(const opk_bound_t* bnd);

extern int
opk_box_project_variables(opk_vector_t* dst,
                          const opk_vector_t* x,
                          const opk_bound_t* xl,
                          const opk_bound_t* xu);

typedef enum {
  OPK_ASCENT  = -1,
  OPK_DESCENT =  1
} opk_orientation_t;

extern int
opk_box_project_direction(opk_vector_t* dst,
                          const opk_vector_t* x,
                          const opk_bound_t* xl,
                          const opk_bound_t* xu,
                          const opk_vector_t* d,
                          opk_orientation_t orientation);

extern int
opk_box_get_free_variables(opk_vector_t* dst,
                           const opk_vector_t* x,
                           const opk_bound_t* xl,
                           const opk_bound_t* xu,
                           const opk_vector_t* d,
                           opk_orientation_t orientation);

extern int
opk_box_get_step_limits(double* smin, double* wolfe, double *smax,
                        const opk_vector_t* x,
                        const opk_bound_t* xl,
                        const opk_bound_t* xu,
                        const opk_vector_t* d,
                        opk_orientation_t orient);

/** @} */

/*---------------------------------------------------------------------------*/
/* VARIABLE METRIC OPTIMIZATION METHOD WITH CONVEX CONSTRAINTS */

/**
 * @addtogroup ConstrainedVariableMetric
 * @{
 */

/** Opaque type for a variable metric optimizer. */
typedef struct _opk_vmlmc opk_vmlmc_t;


/**
 * Create a reverse communication optimizer implementing a limited memory
 * quasi-Newton method with convex constraints and nonmonotone line search.
 *
 * @param vspace - The space to which belong the variables.
 * @param m      - The number of previous steps to memorize.
 * @param stpsiz - The length of the first variable step at start or after a
 *                 restart is the steepest descent scaled to have this
 *                 length.
 */
extern opk_vmlmc_t*
opk_new_vmlmc_optimizer(opk_vspace_t* vspace,
                        opk_index_t m,
                        double stpsiz);

/**
 * Create a reverse communication optimizer implementing a limited memory
 * quasi-Newton method with convex constraints and specific line search.
 *
 * @param vspace - The space to which belong the variables.
 * @param m      - The number of previous steps to memorize.
 * @param stpsiz - The length of the first variable step at start or after a
 *                 restart is the steepest descent scaled to have this
 *                 length.
 * @param lnsrch - The line search method to use.
 */
extern opk_vmlmc_t*
opk_new_vmlmc_optimizer_with_line_search(opk_vspace_t* vspace,
                                         opk_index_t m,
                                         double stpsiz,
                                         opk_lnsrch_t* lnsrch);

extern opk_task_t
opk_start_vmlmc(opk_vmlmc_t* opt);

extern opk_task_t
opk_iterate_vmlmc(opk_vmlmc_t* opt, opk_vector_t* x,
                  double f, opk_vector_t* g, opk_vector_t* d);

extern opk_task_t
opk_get_vmlmc_task(opk_vmlmc_t* opt);

extern opk_index_t
opk_get_vmlmc_iterations(opk_vmlmc_t* opt);

extern opk_index_t
opk_get_vmlmc_evaluations(opk_vmlmc_t* opt);

extern opk_index_t
opk_get_vmlmc_restarts(opk_vmlmc_t* opt);

extern opk_index_t
opk_get_vmlmc_projections(opk_vmlmc_t* opt);

extern int
opk_get_vmlmc_scaling(opk_vmlmc_t* opt);

extern int
opk_set_vmlmc_scaling(opk_vmlmc_t* opt, int scaling);

extern double
opk_get_vmlmc_gatol(opk_vmlmc_t* opt);

extern int
opk_set_vmlmc_gatol(opk_vmlmc_t* opt, double gatol);

extern double
opk_get_vmlmc_grtol(opk_vmlmc_t* opt);

extern int
opk_set_vmlmc_grtol(opk_vmlmc_t* opt, double grtol);

extern double
opk_get_vmlmc_step(opk_vmlmc_t* opt);

extern double
opk_get_vmlmc_stpmin(opk_vmlmc_t* opt);

extern double
opk_get_vmlmc_stpmax(opk_vmlmc_t* opt);

extern int
opk_set_vmlmc_stpmin_and_stpmax(opk_vmlmc_t* opt, double stpmin, double stpmax);

/** @} */

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

/*---------------------------------------------------------------------------*/
/* TRUST REGION STEP */

/**
 * @addtogroup TrustRegion
 * @{
 */

/**
 * Computes a trust region step.
 *
 * Given an `n` by `n` symmetric matrix `A`, an `n`-vector `b`, and a positive
 * number `delta`, this subroutine determines a vector `x` which approximately
 * minimizes the quadratic function:
 *
 *     f(x) = (1/2) x'.A.x + b'.x
 *
 * subject to the Euclidean norm constraint
 *
 *     norm(x) <= delta.
 *
 * This subroutine computes an approximation `x` and a Lagrange multiplier
 * `par` such that either `par` is zero and
 *
 *     norm(x) <= (1 + rtol)*delta,
 *
 * or `par` is positive and
 *
 *     abs(norm(x) - delta) <= rtol*delta.
 *
 * If `xsol` is the exact solution to the problem, the approximation `x`
 * satisfies
 *
 *     f(x) <= ((1 - rtol)^2)*f(xsol)
 *
 * (where `^` means exponentiation).
 *
 *
 * @param n - The order of `A`.
 *
 * @param a - A real array of dimensions `lda` by `n`.  On entry the full upper
 *            triangle of `a` must contain the full upper triangle of the
 *            symmetric matrix `A`.  On exit the array contains the matrix `A`.
 *
 * @param lda - The leading dimension of the array `a`.
 *
 * @param b - A real array of dimension `n`.  On entry `b` specifies the linear
 *            term in the quadratic.  On exit `b` is unchanged.
 *
 * @param   delta - The bound on the Euclidean norm of `x`.
 *
 * @param rtol - The relative accuracy desired in the solution. Convergence
 *            occurs if `f(x) <= ((1-rtol)^2)*f(xsol)` where `xsol` is the
 *            exact solution.
 *
 * @param atol - The absolute accuracy desired in the solution. Convergence
 *            occurs when `norm(x) <= (1+rtol)*delta` and
 *            `max(-f(x),-f(xsol)) <= atol`.
 *
 * @param itmax - The maximum number of iterations.
 *
 * @param par_ptr - If non `NULL`, the address of an integer variable used to
 *            store the Lagrange multiplier.  On entry `*par_ptr` is an initial
 *            estimate of the Lagrange multiplier for the constraint
 *            `norm(x) <= delta`.  On exit `*par_ptr` contains the final estimate
 *            of the multiplier. If `NULL`, the initial Lagrange parameter is 0.
 *
 * @param f_ptr - If non `NULL`, the address of a floating point variable used
 *            to store the function value. On entry `*f_ptr` needs not be
 *            specified.  On exit `*f_ptr` is set to `f(x)` at the output `x`.
 *
 * @param x - A real array of dimension `n`.  On entry `x` need not be
 *            specified.  On exit `x` is set to the final estimate of the
 *            solution.
 *
 * @param iter_ptr - If non `NULL`, the address of an integer variable used to
 *            store the number of iterations.
 *
 * @param z - A real work array of dimension `n`.
 *
 * @param wa1 - A real work array of dimension `n`.
 *
 * @param wa2 - A real work array of dimension `n`.
 *
 *
 * @return
 * The function returns one of the following integer values:
 *
 * * 1 - The function value `f(x)` has the relative accuracy specified by
 *       `rtol`.
 *
 * * 2 - The function value `f(x)` has the absolute accuracy specified by
 *       `atol`.
 *
 * * 3 - Rounding errors prevent further progress.  On exit `x` is the best
 *       available approximation.
 *
 * * 4 - Failure to converge after `itmax` iterations.  On exit `x` is the
 *       best available approximation.
 *
 *
 * @see
 *
 *  * MINPACK-2: `opk_destsv`, `opk_sgqt`.
 *
 *  * LAPACK: `opk_dpotrf`.
 *
 *  * Level 1 BLAS: `opk_dasum`, `opk_daxpy`, `opk_dcopy`, `opk_ddot`,
 *   `opk_dnrm2`, `opk_dscal`.
 *
 *  * Level 2 BLAS: `opk_dtrmv`, `opk_dtrsv`.
 *
 *
 * ### History
 *
 *  * MINPACK-2 Project. July 1994.  Argonne National Laboratory and University
 *    of Minnesota.  Brett M. Averick, Richard Carter, and Jorge J. Moré
 *
 *  * C-version on 30 January 2006 by Éric Thiébaut (CRAL); `info` is the value
 *    returned by the function.
 */
extern int
opk_dgqt(opk_index_t n, double a[], opk_index_t lda, const double b[],
         double delta, double rtol, double atol,
         opk_index_t itmax, double *par_ptr, double *f_ptr,
         double x[], opk_index_t *iter_ptr, double z[],
         double wa1[], double wa2[]);

/**
 * Computes a trust region step.
 *
 * This function is the single precision version of `opk_dgqt`.
 *
 * @see
 *
 *  * MINPACK-2: `opk_sestsv`, `opk_dgqt`.
 *
 *  * LAPACK: `opk_spotrf`
 *
 *  * Level 1 BLAS: `opk_sasum`, `opk_saxpy`, `opk_scopy`, `opk_sdot`,
 *    `opk_snrm2`, `opk_sscal`.
 *
 *  * Level 2 BLAS: `opk_strmv`, `opk_strsv`.
 *
 *
 * ### History
 *
 *  * MINPACK-2 Project. July 1994.  Argonne National Laboratory and University
 *    of Minnesota.  Brett M. Averick, Richard Carter, and Jorge J. Moré
 *
 *  * C-version on 30 January 2006 by Éric Thiébaut (CRAL); `info` is the value
 *    returned by the function.
 */
extern
int opk_sgqt(opk_index_t n, float a[], opk_index_t lda, const float b[],
             float delta, float rtol, float atol,
             opk_index_t itmax, float *par_ptr, float *f_ptr,
             float x[], opk_index_t *iter_ptr, float z[],
             float wa1[], float wa2[]);

/**
 * Computes smallest singular value and corresponding vector from an upper
 * triangular matrix.
 *
 * Given an `n` by `n` upper triangular matrix `R`, this subroutine estimates
 * the smallest singular value and the associated singular vector of `R`.
 *
 * In the algorithm a vector `e` is selected so that the solution `y` to the
 * system {@code R'.y = e} is large. The choice of sign for the components
 * of `e` cause maximal local growth in the components of `y` as the forward
 * substitution proceeds. The vector `z` is the solution of the system
 * {@code R.z = y}, and the estimate `svmin` is {@code norm(y)/norm(z)} in the
 * Euclidean norm.
 *
 * @param n - The order of `R`.
 *
 * @param r - A real array of dimension `ldr` by `n`.  On entry the full upper
 *            triangle must contain the full upper triangle of the matrix `R`.
 *            On exit `r` is unchanged.

 * @param ldr - The leading dimension of `r`.
 *
 * @param z - A real array of dimension `n`.  On entry `z` need not be
 *            specified.  On exit `z` contains a singular vector associated
 *            with the estimate `svmin` such that {@code norm(R.z) = svmin}
 *            and {@code norm(z) = 1} in the Euclidean norm.
 *
 * @return
 * This function returns `svmin`, an estimate for the smallest singular value
 * of `R`.
 *
 * @see
 *  * MINPACK-2: `opk_dgqt`, `opk_sestsv`.
 *
 *  * Level 1 BLAS: `opk_dasum`, `opk_daxpy`, `opk_dnrm2`, `opk_dscal`.
 *
 * ### History
 *
 *  * MINPACK-2 Project. October 1993.
 *     Argonne National Laboratory
 *     Brett M. Averick and Jorge J. Moré.
 *
 *  * C-version on 30 January 2006 by Éric Thiébaut (CRAL);
 *    `svmin` is the value returned by the function.
 */
extern double
opk_destsv(opk_index_t n, const double r[], opk_index_t ldr, double z[]);

/**
 * Computes smallest singular value and corresponding vector from an upper
 * triangular matrix.
 *
 * This function is the single precision version of `opk_destsv`.
 *
 * @see
 *  * MINPACK-2: `opk_sgqt`, `opk_destsv`.
 *
 *  * Level 1 BLAS: `opk_sasum`, `opk_saxpy`, `opk_snrm2`, `opk_sscal`.
 */
extern float
opk_sestsv(opk_index_t n, const float r[], opk_index_t ldr, float z[]);

/** @} */

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
