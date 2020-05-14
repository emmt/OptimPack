/*
 * optimpack.h --
 *
 * Common definitions for OptimPack library.
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

#ifndef _OPTIMPACK_H
#define _OPTIMPACK_H 1

#include <stdio.h>
#include <stddef.h>

/**
 * @defgroup LimitedMemory   Limited memory methods.
 * @defgroup NLCG            Non-linear conjugate gradient methods.
 * @defgroup VariableMetric  Variable metric methods.
 * @defgroup LineSearch      Line search methods.
 * @defgroup ReverseCommunication Reverse communication model.
 * @defgroup Objects         Management of objects and derived types.
 * @defgroup Vectors         Vectors to store variables and vector spaces.
 * @defgroup Operators       Operators acting on vectors.
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
 *
 * @section LargeScale Large scale optimization
 *
 * OptimPack implements limited memory optimization methods dedicated to the
 * minimization of function of a large number of variables.  We routinely use
 * these method to solve image restoration problems with millions, of even
 * billions, of variables.  Two families of methods are provided: non-linear
 * conjugate gradients and variable metric methods.  These methods require the
 * objective function and its gradient thus the objective function has to be \e
 * smooth that is to say \e differentiable.
 *
 * @subsection ReverseCommunicationDriver General reverse communication driver
 *
 * @ref LimitedMemory
 *
 * @ref NLCG
 *
 * @ref VariableMetric
 *
 * @ref LineSearch
 *
 *
 * @section ManagementOfObjects Management of objects and derived types
 *
 * OptimPack implements a simple object-oriented API with various object
 * classes.  The most basic reason of this choice is to correctly manage data
 * and structures which may be shared between different piece of code.  All
 * OptimPack objects inherit from @ref Objects::opk_object_t "opk_object_t"
 * whose instances maintain a reference count and take care of releasing
 * ressources when objects are no longer in use (understand \e referenced).
 * This memory management model is lightweight and very effective and but is
 * not as sophisticated as, say a garbage collector which allows for circular
 * references (circular references must be avoided in OptimPack).
 *
 * @subsection BasicObjects Basic objects
 * @ref Objects
 *
 * @subsection VectorSpaces Vectors and vector spaces
 *
 * A high level driver is provided by OptimPack to the implemented limited
 * memory optimization methods.  At a lower level, these methods operate on
 * variables which can be stored in any form suitable to the application needs.
 * OptimPack will never directly access the components of these variables and
 * only applies a limited set of operations (dot product, scaling, linear
 * combination, *etc.*) to the variables.  These operations are very similar to
 * the one available to \e vectors in mathematics and the notions of vectors
 * and associated vector spaces, operators and linear algebra are central to
 * OptimPack (and to optimization in general).
 *
 * @ref Vectors
 *
 * @ref Operators
 *
 * OptimPack provides a simple memory model to store variables (as conventional
 * C arrays of float/double values) but no memory model is imposed.  For
 * instance, you may have your variables stored in non-conventional memory (GPU)
 * or distributed between several machines.  The only requirement is to provide
 * OptimPack with a set of functions to manipulate your own implementation of
 * vectors.
 *
 *
 * @section Linear Algebra
 *
 * OptimPack has a number of linear algebra routines for small and simple
 * arrays mostly ported from the LINPACK library (in FORTRAN).
 *
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

#define _OPK_STATUS_LIST                                                \
  _OPK_STATUS( 0, OPK_SUCCESS, "Success")                               \
  _OPK_STATUS( 1, OPK_INVALID_ARGUMENT, "Invalid argument")             \
  _OPK_STATUS( 2, OPK_INSUFFICIENT_MEMORY, "Insufficient memory")       \
  _OPK_STATUS( 3, OPK_ILLEGAL_ADDRESS, "Illegal address")               \
  _OPK_STATUS( 4, OPK_NOT_IMPLEMENTED, "Not implemented")               \
  _OPK_STATUS( 5, OPK_CORRUPTED_WORKSPACE, "Corrupted workspace")       \
  _OPK_STATUS( 6, OPK_BAD_SPACE, "Bad variable space")                  \
  _OPK_STATUS( 7, OPK_OUT_OF_BOUNDS_INDEX, "Out of bounds index")       \
  _OPK_STATUS( 8, OPK_NOT_STARTED, "Line search not started")           \
  _OPK_STATUS( 9, OPK_NOT_A_DESCENT, "Not a descent direction")         \
  _OPK_STATUS(10, OPK_STEP_CHANGED, "Step changed")                     \
  _OPK_STATUS(11, OPK_STEP_OUTSIDE_BRACKET, "Step outside bracket")     \
  _OPK_STATUS(12, OPK_STPMIN_GT_STPMAX,                                 \
               "Lower step bound larger than upper bound")              \
  _OPK_STATUS(13, OPK_STPMIN_LT_ZERO,                                   \
              "Minimal step length less than zero")                     \
  _OPK_STATUS(14, OPK_STEP_LT_STPMIN, "Step lesser than lower bound")   \
  _OPK_STATUS(15, OPK_STEP_GT_STPMAX, "Step greater than upper bound")  \
  _OPK_STATUS(16, OPK_FTOL_TEST_SATISFIED,                              \
              "Convergence within variable tolerance")                  \
  _OPK_STATUS(17, OPK_GTOL_TEST_SATISFIED,                              \
                "Convergence within function tolerance")                \
  _OPK_STATUS(18, OPK_XTOL_TEST_SATISFIED,                              \
                "Convergence within gradient tolerance")                \
  _OPK_STATUS(19, OPK_STEP_EQ_STPMAX, "Step blocked at upper bound")    \
  _OPK_STATUS(20, OPK_STEP_EQ_STPMIN, "Step blocked at lower bound")    \
  _OPK_STATUS(21, OPK_ROUNDING_ERRORS_PREVENT_PROGRESS,                 \
              "Rounding errors prevent progress")                       \
  _OPK_STATUS(22, OPK_NOT_POSITIVE_DEFINITE,                            \
              "Operator is not positive definite")                      \
  _OPK_STATUS(23, OPK_BAD_PRECONDITIONER,                               \
              "Preconditioner is not positive definite")                \
  _OPK_STATUS(24, OPK_INFEASIBLE_BOUNDS, "Box set is infeasible")       \
  _OPK_STATUS(25, OPK_WOULD_BLOCK,                                      \
              "Variables cannot be improved (would block)")             \
  _OPK_STATUS(26, OPK_UNDEFINED_VALUE, "Undefined value")               \
  _OPK_STATUS(27, OPK_TOO_MANY_EVALUATIONS, "Too many evaluations")     \
  _OPK_STATUS(28, OPK_TOO_MANY_ITERATIONS, "Too many iterations")


/**
 *  Values returned by OptimPack routines.
 *
 *  `OPK_SUCCESS` indicates that the routine was successfull, any other value
 *  indicate a failure or a warning.
 */
typedef enum {
#define _OPK_STATUS(a,b,c) b = a,
  _OPK_STATUS_LIST
#undef _OPK_STATUS
  OPK_MAX_STATUS
} opk_status_t;

/**
 * Retrieve a textual description for a given status.
 *
 * @param status - The status code.
 *
 * @return A pointer to a string describing the status or an empty string, "",
 *         if the status does not correspond to any known status.
 */
extern const char*
opk_get_reason(opk_status_t status);

/**
 * Retrieve OptimPack status from `errno`.
 *
 * This function is needed to figure out the kind of errors for the few
 * routines which do not return a status (mostly the ones which create
 * objects).
 */
extern opk_status_t
opk_guess_status();

/**
 * Copy a string.
 *
 * @param dst  - The destination buffer to copy the soruce (can be `NULL`).
 * @param size - The number of available bytes in `buf`.
 * @param src  - The source string; `NULL` is considered as being the
 *               same as an empty string "".
 * @return The minimum number of bytes required to store the source
 *         string (including the terminating '\0' character).
 */
extern size_t
opk_copy_string(char* dst, size_t size, const char* src);

/**
 * Get a constant given its name.
 *
 * This function is mostly need for OptimPack wrappers in other languages than
 * C to avoid hardcoding the constants of the library.  For now this function
 * is not particularly efficient so it should be sparsely used.
 *
 * @param name  - The name of the constant; for instance: `"OPK_SUCCESS"`.
 * @param ptr   - The address where to store the value of the constant.
 *
 * @return `OPK_SUCCESS` or `OPK_INVALID_ARGUMENT` if the name is unknown.
 */
extern opk_status_t
opk_get_integer_constant(const char* name, long* ptr);

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

/** Type of the variables in a conventional array. */
typedef enum {
  OPK_FLOAT, /**< Variables are single precision floating point numbers. */
  OPK_DOUBLE /**< Variables are double precision floating point numbers. */
} opk_type_t;

/*---------------------------------------------------------------------------*/
/* BASIC OBJECTS */

/**
 * @addtogroup Objects
 * @{
 *
 * Users of OptimPack library normally do not need to known exactly the
 * contents of an object.  They simply use the functions of the library to
 * manipulate object pointers as opaque structure pointers.
 *
 * Developpers who need to implement derived types or specific variants of
 * existing types must include the header `optimpack-private.h` to unveil
 * the definitions of the structures representing the OptimPack objects.
 */

/**
 * Opaque basic object type.
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

#define OPK_HOLD_VSPACE(vsp)     ((opk_vspace_t*)OPK_HOLD(vsp))
#define OPK_HOLD_VECTOR(vec)     ((opk_vector_t*)OPK_HOLD(vec))
#define OPK_HOLD_LNSRCH(lns)     ((opk_lnsrch_t*)OPK_HOLD(lns))
#define OPK_HOLD_OPERATOR(op)    ((opk_operator_t*)OPK_HOLD(op))
#define OPK_HOLD_CONVEXSET(set)  ((opk_convexset_t*)OPK_HOLD(set))

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
 *                                                  free, dynamic_array);
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
 *                                                 delete_object,
 *                                                 (void*)obj);
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
 * Print vector contents.
 *
 * @param file  - The output file stream, `stdout` is used if `NULL`.
 * @param name  - The name of the vector, can be `NULL`.
 * @param nmax  - The maximum number of elements to print.  The vector size is
 *                used if this parameter is not strictly positive.
 */
extern void
opk_vprint(FILE* file, const char* name, const opk_vector_t* vect,
           opk_index_t nmax);

/**
 * Fetch a specific vector component.
 *
 * This function is by no means intended to be efficient and should be avoided
 * except for debugging purposes.
 *
 * @param vect - A vector.
 * @param k    - The index of the compoent to peek.
 * @param ptr  - The address to store the component value (as a double
 *               precision floating point).
 *
 * @return A standard status.
 */
extern opk_status_t
opk_vpeek(const opk_vector_t* vect, opk_index_t k, double* ptr);

/**
 * Set the value of a specific vector component.
 *
 * This function is by no means intended to be efficient and should be avoided
 * except for debugging purposes.
 *
 * @param vect  - A vector.
 * @param k     - The index of the component to set.
 * @param value - The value to store in the component (as a double precision
 *                floating point).
 *
 * @return A standard status.
 */
extern opk_status_t
opk_vpoke(opk_vector_t* vect, opk_index_t k, double value);

/**
 * Copy the values of a conventional array into a vector.
 *
 * @param dst  - The destination vector.
 * @param src  - The source array.
 * @param type - The type of the elements of the source array.
 * @param n    - The number of elements in the source array.
 *
 * @return A standard status.  The number of elements of the source must match
 *         those of the destination.
 */
extern opk_status_t
opk_vimport(opk_vector_t* dst, const void* src, opk_type_t type, opk_index_t n);

/**
 * Copy the values of a vector into a conventional array.
 *
 * @param dst  - The destination array.
 * @param type - The type of the elements of the destination array.
 * @param n    - The number of elements in the destination array.
 * @param src  - The source vector.
 *
 * @return A standard status.  The number of elements of the source must match
 *         those of the destination.
 */
opk_status_t
opk_vexport(void* dst, opk_type_t type, opk_index_t n, const opk_vector_t* src);

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

extern opk_status_t
opk_apply_direct(opk_operator_t* op, opk_vector_t* dst,
                 const opk_vector_t* src);

extern opk_status_t
opk_apply_adjoint(opk_operator_t* op, opk_vector_t* dst,
                  const opk_vector_t* src);

extern opk_status_t
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
 *
 * @return The former error handler.
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
 * @addtogroup ReverseCommunication
 * @{
 *
 * For maximum portability and integration in any other language, OptimPack
 * limited memory optimizers exploit a reverse communication model to minimize
 * the objective function.
 */

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

/** @} */

/*---------------------------------------------------------------------------*/
/* LINE SEARCH METHODS */

/**
 * @addtogroup LineSearch
 * @{
 *
 * Line search consists in optimizing the variables along a given search
 * direction.
 *
 * The variables `x` are updated along the search direction `d` as follows:
 * ~~~~~{.cpp}
 * x = x0 +/- alpha*d
 * ~~~~~
 * with `x0` the variables at the start of the line search, `alpha` the step
 * size and +/- a plus if `d` is a descent direction and a minus if it is an
 * ascent direction.  The line search approximately solves:
 * ~~~~~{.cpp}
 * min f(x0 +/- alpha*d)
 * ~~~~~
 * with respect to `alpha > 0` and where `f(x)` is the objective function.
 */

/** Opaque line search type.  This sub-type inherits from `opk_object_t`. */
typedef struct _opk_lnsrch opk_lnsrch_t;

/**
 * Create a Moré and Thuente cubic line search.
 *
 * Moré & Thuente cubic line search method is designed to find a step `stp`
 * that satisfies the sufficient decrease condition (a.k.a. first Wolfe
 * condition):
 *
 *     f(stp) ≤ f(0) + ftol⋅stp⋅f'(0)
 *
 * and the curvature condition (a.k.a. second strong Wolfe condition):
 *
 *     abs(f'(stp)) ≤ gtol⋅abs(f'(0))
 *
 * where `f(stp)` is the value of the objective function for a step `stp` along
 * the search direction while `f'(stp)` is the derivative of this function.
 *
 * The algorithm is described in:
 *
 * - J.J. Moré and D.J. Thuente, "Line search algorithms with guaranteed
 *   sufficient decrease" in ACM Transactions on Mathematical Software,
 *   vol. 20, pp. 286–307 (1994).
 *
 * @param ftol specifies the nonnegative tolerance for the sufficient decrease
 *             condition.
 *
 * @param gtol specifies the nonnegative tolerance for the curvature condition.
 *
 * @param xtol specifies a nonnegative relative tolerance for an acceptable
 *             step.  The method exits with a warning if the relative size of
 *             the bracketting interval is less than `xtol`.
 *
 * @return  A line search object.
 */
extern opk_lnsrch_t*
opk_lnsrch_new_csrch(double ftol, double gtol, double xtol);

/**
 * Create a backtracking line search.
 *
 * @param ftol - Parameter of the first Wolfe condition.  Must be in the
 *               range (0,1/2); however a small value is recommended.
 * @param amin - Smallest parameter of the first Wolfe condition.  Must be in
 *               the range (0,1); if larger of equal 1/2, a bisection step is
 *               always take (as in Armijo's rule).
 *
 * @return  A line search object.
 */
extern opk_lnsrch_t*
opk_lnsrch_new_backtrack(double ftol, double amin);

/**
 * Create a nonmonotone line search.
 *
 * Nonmonotone line search is described in SPG2 paper:
 *
 * > E.G. Birgin, J.M. Martínez, & M. Raydan, "Nonmonotone Spectral Projected
 * > Gradient Methods on Convex Sets," SIAM J. Optim. **10**, 1196-1211 (2000).
 *
 * The parameters used in the SPG2 paper:
 * ~~~~~{.cpp}
 * m = 10;
 * ftol = 1E-4;
 * sigma1 = 0.1;
 * sigma2 = 0.9;
 * ~~~~~
 *
 * With `m = 1`, this line search method is equivalent to Armijo's line search
 * except that it attempts to use quadratic interpolation rather than
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

/**
 * Possible values returned by {@link opk_lnsrch_start} and {@link
 * opk_lnsrch_iterate}.
 *
 * These values are for the caller to decide what to do next.  In case of error
 * or warning, {@link opk_lnsrch_get_status} can be used to query more
 * information.
 */
typedef enum {
  OPK_LNSRCH_ERROR       = -1, /**< An error occurred. */
  OPK_LNSRCH_SEARCH      =  0, /**< Line search in progress. */
  OPK_LNSRCH_CONVERGENCE =  1, /**< Line search has converged. */
  OPK_LNSRCH_WARNING     =  2  /**< Line search terminated with warnings. */
} opk_lnsrch_task_t;

/**
 * Start a new line search.
 *
 * @param ls     - The line search object.
 * @param f0     - The function value at the start of the line search (that
 *                 is, for a step of length 0).
 * @param df0    - The directional derivative at the start of the line search.
 *                 (see {@link opk_lnsrch_use_deriv} for the definition).
 * @param stp1   - The length of the first step to try (must be between
 *                 `stpmin` and `stpmax`).
 * @param stpmin - The minimum allowed step length (must be nonnegative).
 * @param stpmax - The maximum allowed step length (must be greater than
 *                 `stpmin`).
 *
 * @return The line search task, which is normally @link OPK_LNSRCH_SEARCH}.
 *         A different value indicates an error.
 */
extern int
opk_lnsrch_start(opk_lnsrch_t* ls, double f0, double df0,
                 double stp1, double stpmin, double stpmax);

/**
 * Check whether line search has converged or update the step size.
 *
 * @param ls      - The line search object.
 * @param stp_ptr - The address at which the step length is stored.  On entry,
 *                  it must be the current step length; on exit, it is updated
 *                  with the next step to try (unless the line search is
 *                  finished).
 * @param f       - The function value for the current step length.
 * @param df      - The directional derivative for the current step length (see
 *                  {@link opk_lnsrch_use_deriv} for the definition).
 *
 * @return The returned value is strictly negative to indicate an error; it is
 * equal to zero when searching is in progress; it is strictly positive when
 * line search has converged or cannot make any more progresses.
 */
extern int
opk_lnsrch_iterate(opk_lnsrch_t* ls, double* stp_ptr,
                   double f, double df);

/**
 * Get current step lenght.
 *
 * @param ls - The line search object.
 *
 * @return Returned value should be >= 0; -1 is returned in case of error.
 */
extern double
opk_lnsrch_get_step(const opk_lnsrch_t* ls);

/**
 * Get current line search pending task.
 *
 * In case of error or warning, {@link opk_lnsrch_get_status} can be used to
 * query more information.
 *
 * @param ls - The line search object.
 *
 * @return The current line search pending task.
 */
extern opk_lnsrch_task_t
opk_lnsrch_get_task(const opk_lnsrch_t* ls);

/**
 * Get current line search status.
 *
 * @param ls - The line search object.
 *
 * @return The current line search status.
 */
extern opk_status_t
opk_lnsrch_get_status(const opk_lnsrch_t* ls);

/**
 * Check whether line search has errors.
 *
 * @param ls - The line search object.
 *
 * @return A boolean value, true if there are any errors.
 */
extern opk_bool_t
opk_lnsrch_has_errors(const opk_lnsrch_t* ls);

/**
 * Check whether line search has warnings.
 *
 * @param ls - The line search object.
 *
 * @return A boolean value, true if there are any warnings.
 */
extern opk_bool_t
opk_lnsrch_has_warnings(const opk_lnsrch_t* ls);

/**
 * Check whether line search has converged.
 *
 * @param ls - The line search object.
 *
 * @return A boolean value, true if line search has converged.
 */
extern opk_bool_t
opk_lnsrch_converged(const opk_lnsrch_t* ls);

/**
 * Check whether line search has finished.
 *
 * Note that line search termination may be due to convergence, possibly with
 * warnings, or to errors.
 *
 * @param ls - The line search object.
 *
 * @return A boolean value, true if line search has finished.
 */
extern opk_bool_t
opk_lnsrch_finished(const opk_lnsrch_t* ls);

/**
 * Check whether line search requires the directional derivative.
 *
 * The directional derivative is the derivative with respect to the step size
 * of the objective function during the line search.  The directional
 * derivative is given by:
 * ~~~~~{.cpp}
 *  +/-d'.g(x0 +\- alpha*d)
 * ~~~~~
 * with `d` the search direction, `g(x)` the gradient of the objective function
 * (the `+/-` signs are `+` if `d` is a descent direction and `-` if it is an
 * ascent direction) and assuming that, during the line search, the variables
 * change as:
 * ~~~~~{.cpp}
 * x = x0 +/- alpha*d
 * ~~~~~
 * with `x0` the variables at the start of the line search and `alpha` the step
 * size.
 *
 * @param ls - The line search object.
 *
 * @return A boolean value, true if line search makes use of the directional
 *         derivative.
 */
extern opk_bool_t
opk_lnsrch_use_deriv(const opk_lnsrch_t* ls);

/** Moré & Thuente method to perform a cubic safeguarded step. */
extern opk_status_t
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
 * Structure used to specify the settings of a NLCG optimizer.
 *
 * The absolute threshold for the norm or the gradient for convergence are
 * specified by the members `gatol` and `grtol` of this structure.  The
 * convergence of the non-linear convergence gradient (NLCG) method is defined
 * by:
 * ~~~~~{.cpp}
 * ||g|| <= max(0, gatol, grtol*||ginit||)
 * ~~~~~
 * where `||g||` is the Euclidean norm of the current gradient `g`, `||ginit||`
 * is the Euclidean norm of the initial gradient `ginit` while `gatol` and
 * `grtol` are absolute and relative thresholds.
 *
 * During a line search, the step is constrained to be within `stpmin` and
 * `stpmax` times the lenght of the first step.  The relative bounds must be
 * such that:
 * ~~~~~{.cpp}
 * 0 <= stpmin < stpmax
 * ~~~~~
 */
typedef struct _opk_nlcg_options {
  double          delta; /**< Relative size for a small step. */
  double        epsilon; /**< Threshold to accept descent direction. */
  double          grtol; /**< Relative threshold for the norm or the gradient
                              (relative to the norm of the initial gradient)
                              for convergence. */
  double          gatol; /**< Absolute threshold for the norm or the gradient
                              for convergence. */
  double         stpmin; /**< Relative minimum step length. */
  double         stpmax; /**< Relative maximum step length. */
  double           fmin; /**< Minimal function value if provided. */
  unsigned int    flags; /**< A bitwise combination of the non-linear conjugate
                              gradient update method and options. */
  opk_bool_t fmin_given; /**< Minimal function value is provided? */
} opk_nlcg_options_t;

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

/* For instance: (OPK_NLCG_POLAK_RIBIERE_POLYAK | OPK_NLCG_POWELL) merely
   corresponds to PRP+ (Polak, Ribiere, Polyak) while (OPK_NLCG_PERRY_SHANNO |
   OPK_NLCG_SHANNO_PHUA) merely corresponds to the conjugate gradient method
   implemented in CONMIN. */

/* Default settings for non linear conjugate gradient (should correspond to
   the method which is, in general, the most successful). */
#define OPK_NLCG_DEFAULT (OPK_NLCG_POLAK_RIBIERE_POLYAK | \
                          OPK_NLCG_POWELL | OPK_NLCG_SHANNO_PHUA)

/**
 * Query default nonlinear conjugate gradient optimizer parameters.
 *
 * @param opts - The address of the structure where to store the parameters.
 */
extern void
opk_get_nlcg_default_options(opk_nlcg_options_t* opts);

/**
 * Check nonlinear conjugate gradient optimizer parameters.
 *
 * @param opts - The address of the structure with the parameters to check.
 *
 * @return A standard status.
 */
extern opk_status_t
opk_check_nlcg_options(const opk_nlcg_options_t* opts);

/**
 * Create a new optimizer instance for non-linear conjugate gradient method.
 *
 * This function creates an optimizer instance for minimization by a non-linear
 * conjugate gradient method over a given vector space.  The returned instance
 * must be unreferenced by calling the function `opk_drop_object()`, or th
 * macro `OPK_DROP()` when no longer needed.
 *
 * @param vspace   The vector space of the unknowns.
 * @param lnsrch - Optional line search method to use; can be `NULL` to use a
 *                 default one.  Note that the optimizer will hold a reference
 *                 to the line search object.
 *
 * @return The address of a new optimizer instance, or NULL in case of error.
 *         Global variable errno may be ENOMEM if there is not enough memory
 *         or EINVAL if one of the arguments is invalid or EFAULT if vspace is
 *         NULL.
 */
extern opk_nlcg_t*
opk_new_nlcg_optimizer(const opk_nlcg_options_t* opts,
                       opk_vspace_t* vspace,
                       opk_lnsrch_t* lnsrch);

extern opk_task_t
opk_start_nlcg(opk_nlcg_t* opt, opk_vector_t* x);

extern opk_task_t
opk_iterate_nlcg(opk_nlcg_t* opt, opk_vector_t* x,
                 double f, opk_vector_t* g);

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
opk_get_nlcg_step(const opk_nlcg_t* opt);

extern double
opk_get_nlcg_gnorm(const opk_nlcg_t* opt);

extern opk_index_t
opk_get_nlcg_evaluations(const opk_nlcg_t* opt);

extern opk_index_t
opk_get_nlcg_iterations(const opk_nlcg_t* opt);

extern opk_index_t
opk_get_nlcg_restarts(const opk_nlcg_t* opt);

extern opk_task_t
opk_get_nlcg_task(const opk_nlcg_t* opt);

extern opk_status_t
opk_get_nlcg_status(const opk_nlcg_t* opt);

extern size_t
opk_get_nlcg_name(char* buf, size_t size, const opk_nlcg_t* opt);

/**
 * Get description of nonlinear conjugate gradient method.
 *
 * @param buf - A string buffer to copy the description (can be `NULL`).
 * @param size - The number of available bytes in `buf`.
 * @param opt - The optimizer.
 *
 * @return The minimum number of bytes required to store the description
 *         (including the terminating '\0' character).
 */
extern size_t
opk_get_nlcg_description(char* buf, size_t size, const opk_nlcg_t* opt);

extern double
opk_get_nlcg_beta(const opk_nlcg_t* opt);

/** @} */

/*---------------------------------------------------------------------------*/
/* CONVEX SETS AND BOX SETS */

/**
 * @addtogroup BoxConstraints
 * @{
 */

/** Opaque structure to represent an instance of a convex set. */
typedef struct _opk_convexset opk_convexset_t;

/**
 * Type of bounds.
 *
 * Variables can be bounded or unbounded from below and from above.  The bounds
 * are specified by two parameters: a type and an address.  If the variables
 * are unbounded, then the corresponding bound type is `OPK_BOUND_NONE` and the
 * associated address must be `NULL`.  If all variables have the same bound, it
 * is more efficient to specify a scalar bound with type
 * `OPK_BOUND_SCALAR_FLOAT` or `OPK_BOUND_SCALAR_DOUBLE` and the address of a
 * single or double precision variable which stores the bound.  It is also
 * possible to specify a componentwise bound in three different ways depending
 * how the bounds are stored.  If the bounds are in an OptimPack vector (of the
 * same vector space of the variables) use type `OPK_BOUND_VECTOR` and provide
 * the address of the vector.  If the bounds are in a conventional array (with
 * the same number of elements as the variables), then the address of the array
 * must be provided with type `OPK_BOUND_STATIC_FLOAT` or
 * `OPK_BOUND_STATIC_DOUBLE` if the array will not be released while the bound
 * is in use or `OPK_BOUND_VOLATILE_FLOAT` or `OPK_BOUND_VOLATILE_DOUBLE`
 * otherwise.
 */
typedef enum {
  OPK_BOUND_NONE            = 0, /**< No-bound (associated value must be
                                  *   `NULL`). */
  OPK_BOUND_SCALAR_FLOAT    = 1, /**< Scalar bound (associated value is the
                                  *   address of a float). */
  OPK_BOUND_SCALAR_DOUBLE   = 2, /**< Scalar bound (associated value is the
                                  *   address of a double). */
  OPK_BOUND_STATIC_FLOAT    = 3, /**< Bounds are stored in a static array of
                                  *   float's (associated value is the address
                                  *   of the array). */
  OPK_BOUND_STATIC_DOUBLE   = 4, /**< Bounds are stored in a static array of
                                  *   double's (associated value is the address
                                  *   of the array). */
  OPK_BOUND_VOLATILE_FLOAT  = 5, /**< Bounds are stored in a volatile array of
                                  *   floats's (associated value is the address
                                  *   of the array). */
  OPK_BOUND_VOLATILE_DOUBLE = 6, /**< Bounds are stored in a volatile array of
                                  *   double's (associated value is the address
                                  *   of the array). */
  OPK_BOUND_VECTOR          = 7  /**< Vector bound (associated value is the
                                  *   address of an `opk_vector_t`). */
} opk_bound_type_t;

/**
 * Project the variables to the feasible set.
 *
 * Given input variables `x`, the projection produces feasible output variables
 * `dst` which belongs to the convex set.  The input and output variables can
 * be stored in the same vector (*i.e.* the method can be applied *in-place*).
 *
 * The function `opk_can_project_variables()` can be used to determine whether
 * this functionality is implemented by `set`.
 *
 * @param dst - The output projected variables.
 * @param x   - The input variables.
 * @param set - The convex set which implements the constraints.
 *
 * @return `OPK_SUCCESS` or `OPK_ILLEGAL_ADDRESS` if one argumeant has an
 *         invalid (`NULL`) address, `OPK_BAD_SPACE` if not all vectors
 *         and bounds belongs to the same space, `OPK_NOT_IMPLEMENTED` if
 *         this functionality is not implemented.
 */
extern opk_status_t
opk_project_variables(opk_vector_t* dst,
                      const opk_vector_t* x,
                      const opk_convexset_t* set);

/**
 * Check whether the projection of the variables is implemented.
 *
 * @param set - The convex set which implements the constraints.
 *
 * @return A boolean value, true if the projection of the variables is
 *         implemented by `set`.
 */
extern opk_bool_t
opk_can_project_variables(const opk_convexset_t* set);

/**
 * Orientation of a search direction.
 *
 * If orientation is `OPK_DESCENT` (or strictly positive), the search
 * direction `d` is a descent direction and the variables are updated as:
 * ~~~~~~~~~~{.c}
 *     x[i] + alpha*d[i]
 * ~~~~~~~~~~
 * otherwise, `d` is considered as an ascent disrection and the variables
 * are updated as:
 * ~~~~~~~~~~{.c}
 *     x[i] - alpha*d[i]
 * ~~~~~~~~~~
 */
typedef enum {
  OPK_ASCENT  = -1, /**< Ascent search direction. */
  OPK_DESCENT =  1  /**< Descent search direction. */
} opk_orientation_t;

/**
 * Project a direction.
 *
 * This function projects the direction `d` so that:
 * ~~~~~~~~~~{.c}
 * x + orient*alpha*d
 * ~~~~~~~~~~
 * yields a feasible position for `alpha > 0` sufficiently small.
 *
 * The function `opk_can_project_directions()` can be used to determine whether
 * this functionality is implemented by `set`.
 *
 * @param dst - The resulting projected direction.
 * @param x   - The current variables.
 * @param set - The convex set which implements the constraints.
 * @param d   - The direction to project.
 * @param orient - The orientation of the direction `d`.  Strictly positive if
 *              `d` is a descent direction, strictly negative if `d` is an
 *              ascent direction.  For convenience, the constants {@link
 *              #OPK_DESCENT} and {@link #OPK_ASCENT} can be used to specify
 *              the orientation.
 *
 * @return `OPK_SUCCESS` or `OPK_ILLEGAL_ADDRESS` if one argumeant has an
 *         invalid (`NULL`) address, `OPK_BAD_SPACE` if not all vectors
 *         and bounds belongs to the same space, `OPK_NOT_IMPLEMENTED` if
 *         this functionality is not implemented.
 */
extern opk_status_t
opk_project_direction(opk_vector_t* dst,
                      const opk_vector_t* x,
                      const opk_convexset_t* set,
                      const opk_vector_t* d,
                      opk_orientation_t orient);

/**
 * Check whether projection of the search direction is implemented.
 *
 * @param set - The convex set which implements the constraints.
 *
 * @return A boolean value, true if projection of the search direction is
 *         implemented by `set`.
 */
extern opk_bool_t
opk_can_project_directions(const opk_convexset_t* set);

/**
 * Find the non-binding constraints.
 *
 * The function `opk_can_get_free_variables()` can be used to determine whether
 * this functionality is implemented by `set`.
 *
 * @param dst - The resulting mask whose elements are set to 1 (resp. 0) if
 *              the corresponding variables are free (resp. binded).
 * @param x   - The current variables.
 * @param set - The convex set which implements the constraints.
 * @param d   - The search direction.
 * @param orient - The orientation of the search direction (see {@link
 *              opk_project_direction}).
 *
 * @return `OPK_SUCCESS` or `OPK_ILLEGAL_ADDRESS` if one argumeant has an
 *         invalid (`NULL`) address, `OPK_BAD_SPACE` if not all vectors
 *         and bounds belongs to the same space, `OPK_NOT_IMPLEMENTED` if
 *         this functionality is not implemented.
 */
extern opk_status_t
opk_get_free_variables(opk_vector_t* dst,
                       const opk_vector_t* x,
                       const opk_convexset_t* set,
                       const opk_vector_t* d,
                       opk_orientation_t orient);

/**
 * Check whether determining the free variables is implemented.
 *
 * @param set - The convex set which implements the constraints.
 *
 * @return A boolean value, true if determining the free variables is
 *         implemented by `set`.
 */
extern opk_bool_t
opk_can_get_free_variables(const opk_convexset_t* set);

/**
 * Find the limits of the step size.
 *
 * Along the search direction the new variables are computed as:
 * ~~~~~~~~~~{.c}
 * proj(x +/- alpha*d)
 * ~~~~~~~~~~
 * where `proj` is the projection onto the feasible set, `alpha` is the step
 * length and +/- is a plus for a descent direction and a minus otherwise.
 * This method computes 3 specific step lengths: `smin1` which is the step
 * length for the first bound reached by the variables along the search
 * direction; `smin2` which is the step length for the first bound reached by
 * the variables along the search direction and with a non-zero step length;
 * `smax` which is the step length for the last bound reached by the variables
 * along the search direction.
 *
 * In other words, for any `0 <= alpha <= smin1`, no variables may overreach a
 * bound, while, for any `alpha >= smax`, the projected variables are the same
 * as those obtained with `alpha = smax`.  If `d` has been properly projected
 * (e.g. by {@link opk_project_direction}), then `smin1 = smin2` otherwise
 * `0 <= smin1 <= smin2` and `smin2 > 0`.
 *
 * The function `opk_can_get_step_limits()` can be used to determine whether
 * this functionality is implemented by `set`.
 *
 * @param smin1  - The address to store the value of `smin1` or `NULL`.
 * @param smin2  - The address to store the value of `smin2` or `NULL`.
 * @param smax   - The address to store the value of `smax` or `NULL`.
 * @param x      - The current variables (assumed to be feasible).
 * @param set - The convex set which implements the constraints.
 * @param d      - The search direction.
 * @param orient - The orientation of the search direction (see
 *                 {@link opk_project_direction}).
 *
 * @return `OPK_SUCCESS` or `OPK_ILLEGAL_ADDRESS` if one argumeant has an
 *         invalid (`NULL`) address, `OPK_BAD_SPACE` if not all vectors
 *         and bounds belongs to the same space, `OPK_NOT_IMPLEMENTED` if
 *         this functionality is not implemented.
 */
extern opk_status_t
opk_get_step_limits(double* smin1, double* smin2, double *smax,
                    const opk_vector_t* x,
                    const opk_convexset_t* set,
                    const opk_vector_t* d,
                    opk_orientation_t orient);

/**
 * Check whether determining the step limits is implemented.
 *
 * @param set - The convex set which implements the constraints.
 *
 * @return A boolean value, true if determining the step limits is implemented
 *         by `set`.
 */
extern opk_bool_t
opk_can_get_step_limits(const opk_convexset_t* set);

/**
 * Create a new box set.
 *
 * A box set is a convex set with separable bounds on the variables.
 *
 * @param space - The vector space of the variables.  The convex set is a
 *                subset of this vector space.
 * @param lower_type - The type of the lower bound(s).
 * @param lower - The value of the lower bound(s).
 * @param upper_type - The type of the upper bound(s).
 * @param upper - The value of the upper bound(s).
 *
 * @return The address of a new box set, `NULL` in case of error (global
 *         variable `errno` can be used to figure out the reason of the
 *         failure).
 */
extern opk_convexset_t*
opk_new_boxset(opk_vspace_t* space,
               opk_bound_type_t lower_type, void* lower,
               opk_bound_type_t upper_type, void* upper);

/** @} */

/*---------------------------------------------------------------------------*/
/* VARIABLE METRIC OPTIMIZATION METHOD WITH OPTIONAL BOUND CONSTRAINTS */

/**
 * @addtogroup VariableMetric
 * @{
 *
 * Variable metric methods, also known as quasi-Newton methods, make use of
 * the previous steps and gradient changes to estimate an approximation of the
 * inverse Hessian of the objective function.
 */

/** Opaque type for a variable metric optimizer. */
typedef struct _opk_vmlmb opk_vmlmb_t;

/** Rules for scaling the inverse Hessian approximation. */
typedef enum {
  OPK_SCALING_NONE             = 0, /**< No-scaling. */
  OPK_SCALING_OREN_SPEDICATO   = 1, /**< Scaling by:
                                         {@code gamma = (s'.y)/(y'.y)} */
  OPK_SCALING_BARZILAI_BORWEIN = 2  /**< Scaling by:
                                         {@code gamma = (s'.s)/(s'.y)} */
} opk_bfgs_scaling_t;

/** Structure used to store the settings of a VMLMB optimizer. */
typedef struct _opk_vmlmb_options {
  double           delta; /**< Relative size for a small step. */
  double         epsilon; /**< Threshold to accept descent direction. */
  double           grtol; /**< Relative threshold for the norm or the projected
                               gradient (relative to the norm of the initial
                               projected gradient) for convergence. */
  double           gatol; /**< Absolute threshold for the norm or the projected
                               gradient for convergence. */
  double          stpmin; /**< Relative minimum step length. */
  double          stpmax; /**< Relative maximum step length. */
  opk_index_t        mem; /**< Maximum number of memorized steps. */
  opk_bool_t       blmvm; /**< Emulate Benson & Moré BLMVM method? */
  opk_bool_t save_memory; /**< Save some memory? */
} opk_vmlmb_options_t;

/**
 * Query default VMLMB optimizer parameters.
 *
 * @param opts - The address of the structure where to store the parameters.
 */
extern void
opk_get_vmlmb_default_options(opk_vmlmb_options_t* opts);

/**
 * Check VMLMB optimizer parameters.
 *
 * @param opts - The address of the structure with the parameters to check.
 *
 * @return A standard status.
 */
extern opk_status_t
opk_check_vmlmb_options(const opk_vmlmb_options_t* opts);

/**
 * Create a reverse communication optimizer implementing a limited memory
 * quasi-Newton method.
 *
 * The optimizer may account for bound constraints if argument `box` is
 * non-`NULL`.
 *
 * @param opts   - The options to sue (can be `NULL` to use default options).
 * @param space  - The space to which belong the variables.
 * @param lnsrch - Optional line search method to use; can be `NULL` to use a
 *                 default one.  Note that the optimizer will hold a reference
 *                 to the line search object.
 * @param box    - An optional box set implementing the bound constraints (can
 *                 be `NULL`).
 *
 * @return The address of a new optimizer instance, or NULL in case of error.
 *         Global variable errno may be ENOMEM if there is not enough memory
 *         or EINVAL if one of the arguments is invalid or EFAULT if vspace is
 *         NULL.
 */
extern opk_vmlmb_t*
opk_new_vmlmb_optimizer(const opk_vmlmb_options_t* opts,
                        opk_vspace_t* space,
                        opk_lnsrch_t* lnsrch,
                        opk_convexset_t* box);

/** The variants implemented by VMLMB. */
typedef enum { OPK_LBFGS, OPK_VMLMB, OPK_BLMVM } opk_vmlmb_method_t;
extern opk_vmlmb_method_t
opk_get_vmlmb_method(const opk_vmlmb_t* opt);

extern const char*
opk_get_vmlmb_method_name(const opk_vmlmb_t* opt);

extern opk_task_t
opk_start_vmlmb(opk_vmlmb_t* opt, opk_vector_t* x);

extern opk_task_t
opk_iterate_vmlmb(opk_vmlmb_t* opt, opk_vector_t* x,
                  double f, opk_vector_t* g);

extern opk_task_t
opk_get_vmlmb_task(const opk_vmlmb_t* opt);

extern opk_status_t
opk_get_vmlmb_status(const opk_vmlmb_t* opt);

/**
 * Get description of algorithm implemented by VMLMB.
 *
 * @param buf - A string buffer to copy the description (can be `NULL`).
 * @param size - The number of available bytes in `buf`.
 * @param opt - The optimizer.
 *
 * @return The minimum number of bytes required to store the description
 *         (including the terminating '\0' character).
 */
extern size_t
opk_get_vmlmb_description(char* buf, size_t size, const opk_vmlmb_t* opt);

extern opk_index_t
opk_get_vmlmb_evaluations(const opk_vmlmb_t* opt);

extern opk_index_t
opk_get_vmlmb_iterations(const opk_vmlmb_t* opt);

extern opk_index_t
opk_get_vmlmb_restarts(const opk_vmlmb_t* opt);

extern double
opk_get_vmlmb_step(const opk_vmlmb_t* opt);

extern double
opk_get_vmlmb_gnorm(const opk_vmlmb_t* opt);

/**
 * Get actual number of memorized steps.
 *
 * @param opt - The VMLMB optimizer.
 *
 * @return The actual number of memorized steps which is in the range `[0,m]`,
 *         with `m` the maximum number of memorized steps.
 */
extern opk_index_t
opk_get_vmlmb_mp(const opk_vmlmb_t* opt);

/**
 * Get a given memorized variable change.
 *
 * Variable metric methods store variable and gradient changes for the few last
 * steps to measure the effect of the Hessian.  Using pseudo-code notation the
 * following `(s,y)` pairs are memorized:
 * ~~~~~{.cpp}
 * s[k-j] = x[k-j+1] - x[k-j]     // variable change
 * y[k-j] = g[k-j+1] - g[k-j]     // gradient change
 * ~~~~~
 * with `x[k]` and `g[k]` the variables and corresponding gradient at `k`-th
 * iteration and `j=1,...,mp` the relative index of the saved pair.
 *
 * @param opt - The VMLMB optimizer.
 *
 * @param j   - The index of the memorized pair to consider relative to the
 *              last one.  Index `j` must be in the inclusive range `[0,mp]`
 *              with `mp` the actual number of saved corrections.  The special
 *              case `j = 0` corresponds to the next saved pair (which will be
 *              overwritten at the end of the current iteration).  The other
 *              cases, `j = 1,...,mp`, correspond to actually saved pairs.
 *              Because the algorithm may be automatically restarted or may try
 *              to save memory, the actual number of saved pairs `mp` may
 *              change between iterations and has to be retrieved each time a
 *              given save pair is queried.
 *
 * @return The variable difference `s[k-j]` where `k` is the current iteration
 *         number, with `m` the maximum number of memorized steps.  `NULL` is
 *         returned if `j` is out of bounds.
 *
 * @see opk_get_vmlmb_y, opk_get_vmlmb_mb.
 */
extern opk_vector_t*
opk_get_vmlmb_s(const opk_vmlmb_t* opt, opk_index_t j);

/**
 * Get a given memorized gradient change.
 *
 * @param opt - The VMLMB optimizer.
 * @param j   - The index of the memorized pair to consider relative to the
 *              last one.  Index `j` must be in the inclusive range `[0,mp]`
 *              with `mp` the actual number of saved corrections.
 *
 * @return The gradient difference `y[k-j]` where `k` is the current iteration
 *         number, with `m` the maximum number of memorized steps.  `NULL` is
 *         returned if `j` is out of bounds.
 *
 * @see opk_get_vmlmb_s, opk_get_vmlmb_mb.
 */
extern opk_vector_t*
opk_get_vmlmb_y(const opk_vmlmb_t* opt, opk_index_t j);

/** @} */

/*---------------------------------------------------------------------------*/
/* SIMPLE DRIVER FOR LIMITED MEMORY OPTIMIZATION */

/**
 * @addtogroup LimitedMemory
 * @{
 */

/** Limited memory optimization algorithm. */
typedef enum {
  OPK_ALGORITHM_NLCG, /**< Nonlinear conjugate gradient. */
  OPK_ALGORITHM_VMLMB /**< Limited memory variable metric (possibly with
                       *   bounds). */
} opk_algorithm_t;

/** Opaque structure for limited memory optimizer. */
typedef struct _opk_optimizer opk_optimizer_t;

/**
 * Create a reverse communication optimizer implementing a limited memory
 * optimization method.
 *
 * This simple driver assumes that the variables and the gradient are stored as
 * flat arrays of floating point values (`float` or `double`).  The implemented
 * algorithms are suitable for large problems and smooth (that is to say
 * differentiable) objective functions.
 *
 * Depending on the settings, optimization can be performed by an instance of
 * the nonlinear conjugate gradient methods or an instance of the limited
 * memory quasi-Newton (a.k.a. variable metric) methods.  Simple bound
 * constraints can be taken into account (providing the variable metric method
 * is selected).
 *
 * When no longer needed, the optimizer must be released with
 * {@link opk_destroy_optimizer}.
 *
 * Typical usage is:
 * ~~~~~{.cpp}
 * const int n = 100000;
 * const int type = OPK_FLOAT;
 * double fx;
 * float x[n];
 * float gx[n];
 * opk_optimizer_t* opt = opk_new_optimizer(OPK_ALGORITHM_VMLMB,
 *                                          NULL, type, n,
 *                                          OPK_BOUND_NONE, NULL,
 *                                          OPK_BOUND_NONE, NULL,
 *                                          NULL);
 * task = opk_start(opt, x);
 * for (;;) {
 *     if (task == OPK_TASK_COMPUTE_FG) {
 *          fx = f(x);
 *          for (i = 0; i < n; ++i) {
 *              gx[i] = ...;
 *          }
 *      } else if (task == OPK_TASK_NEW_X) {
 *          // a new iterate is available
 *          fprintf(stdout, "iter=%ld, f(x)=%g, |g(x)|=%g\n",
 *                  opk_get_iterations(opt), fx,
 *                  opk_get_gnorm(opt));
 *      } else {
 *          break;
 *      }
 *      task = opk_iterate(opt, x, fx, gx);
 * }
 * if (task != OPK_TASK_FINAL_X) {
 *     fprintf(stderr, "some error occured (%s)",
 *             opk_get_reason(opk_get_status(opt)));
 * }
 * opk_destroy_optimizer(opt);
 * ~~~~~
 *
 * @param algorithm - The limited memory algorithm to use.
 * @param opts   - Address of structure with algorithm options (can be `NULL`).
 * @param type   - The type of the variable limited memory algorithm to use
 *                 ({@link OPK_FLOAT} or {@link OPK_DOUBLE}).
 * @param n      - The number of variables (`n` > 0).
 * @param lower_type - The type of the lower bound.
 * @param lower  - Optional lower bound for the variables.  Can be `NULL` if
 *                 there are no lower bounds and `lower_type` is {@link
 *                 OPK_BOUND_NONE}; otherwise must have the same type as the
 *                 variables.
 * @param upper_type - The type of the upper bound.
 * @param upper  - Optional upper bound for the variables.  Can be `NULL` if
 *                 there are no upper bounds and `upper_type` is {@link
 *                 OPK_BOUND_NONE}; otherwise must have the same type as the
 *                 variables.
 * @param lnsrch - The line search method to use, can be `NULL` to use a default
 *                 line search which depends on the optimization algorithm.
 *
 * @return The address of a new optimizer instance, or `NULL` in case of error.
 *         In case of error, global variable `errno` may be `ENOMEM` if there
 *         is not enough memory or `EINVAL` if one of the arguments is invalid
 *         or `EFAULT` if the bounds are unexpectedly `NULL`.
 *
 */
extern opk_optimizer_t *
opk_new_optimizer(opk_algorithm_t algorithm, const void* opts,
                  opk_type_t type, opk_index_t n,
                  opk_bound_type_t lower_type, void* lower,
                  opk_bound_type_t upper_type, void* upper,
                  opk_lnsrch_t* lnschr);

/**
 * Destroy a reverse communication optimizer implementing a limited memory
 * optimization method.
 *
 * This function must be called when the optimizer is no longer in use.  It
 * reduces the reference count of the optimizer eventually freeing any
 * associated ressources.
 *
 * @param opt - An optimizer created by {@link #opk_new_optimizer}.
 */
extern void
opk_destroy_optimizer(opk_optimizer_t *opt);

/**
 * Start the optimization with given initial variables.
 *
 * @param opt  - An optimizer created by {@link #opk_new_optimizer}.
 * @param x    - The variables which must be of the same type (`float` or
 *               `double`) as assumed when creating the optimizer and have the
 *               correct number of elements.
 *
 * @return An integer indicating the next thing to do for the caller.
 *         Unless an error occured, it should be {@link OPK_TASK_COMPUTE_FG}.
 */
extern opk_task_t
opk_start(opk_optimizer_t *opt, void* x);

/**
 * Proceed with next optimization step.
 *
 * Note that the variables must not be changed by the caller after calling
 * {@link opk_start} and between calls to {@link opk_iterate}.  Arrays `x` and
 * `g` must be of the same type (`float` or `double`) as assumed when creating
 * the optimizer and have the correct number of elements.
 *
 * @param opt  - An optimizer created by {@link #opk_new_optimizer}.
 * @param x    - The current variables.
 * @param f    - The value of the objective function at the current variables.
 * @param g    - The gradient of the objective function at the current
 *               variables.
 *
 * @return An integer indicating the next thing to do for the caller.
 */
extern opk_task_t
opk_iterate(opk_optimizer_t *opt, void* x, double f, void* g);

/**
 * Get the current pending task.
 *
 * @param opt  - An optimizer created by {@link #opk_new_optimizer}.
 *
 * @return The current pending task.
 */
extern opk_task_t
opk_get_task(const opk_optimizer_t* opt);

/**
 * Get the current optimizer status.
 *
 * This function is useful to figure out which kind of problem occured when
 * the pending task is {@link OPK_TASK_WARNING} or {@link OPK_TASK_ERROR}.
 *
 * @param opt  - An optimizer created by {@link #opk_new_optimizer}.
 *
 * @return The current optimizer status.
 */
extern opk_status_t
opk_get_status(const opk_optimizer_t* opt);

extern opk_index_t
opk_get_evaluations(const opk_optimizer_t* opt);

extern opk_index_t
opk_get_iterations(const opk_optimizer_t* opt);

extern opk_index_t
opk_get_restarts(const opk_optimizer_t* opt);

extern size_t
opk_get_name(char* buf, size_t size, const opk_optimizer_t* opt);

extern size_t
opk_get_description(char* buf, size_t size, const opk_optimizer_t* opt);

extern double
opk_get_step(const opk_optimizer_t* opt);

extern double
opk_get_gnorm(const opk_optimizer_t* opt);

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
 *            estimate of the Lagrange multiplier for the constraint `norm(x)
 *            <= delta`.  On exit `*par_ptr` contains the final estimate of the
 *            multiplier.  If `NULL`, the initial Lagrange parameter is 0.
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
