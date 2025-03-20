/*
 * bobyqa.h -
 *
 * Definitions for Mike Powell's BOBYQA algorithm for minimizing a function of many
 * variables. The method is "derivatives free" (only the function values are needed) and
 * accounts for bound constraints on the variables. The algorithm is described in:
 *
 *     M.J.D. Powell, "The BOBYQA Algorithm for Bound Constrained Optimization Without
 *     Derivatives." Technical report, Department of Applied Mathematics and Theoretical
 *     Physics, University of Cambridge (2009).
 *
 * The present code is based on the original FORTRAN version written by Mike Powell who
 * kindly provides his code on demand (at mjdp@cam.ac.uk) and has been converted to C by
 * É. Thiébaut.
 *
 * Copyright (c) 2009, Mike Powell (FORTRAN version).
 * Copyright (c) 2015-2025, Éric Thiébaut (C version).
 *
 * Read the accompanying `LICENSE` file for details.
 */

#ifndef BOBYQA_H_
#define  BOBYQA_H_ 1

#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Prototype of the objective function assumed by the BOBYQA routine. The returned value
   is the function value at X the current variables, N is the number of variables and DATA
   is anything needed by the function (unused by BOBYQA itself). */
typedef double bobyqa_objfun(const ptrdiff_t n, const double* x, void* data);

/**
 * @brief Status for BOBYQA routines.
 *
 * This type enumerate the possible values returned by bobyqa(),
 * bobyqa_get_status() and bobyqa_iterate().
 */
typedef enum {
    BOBYQA_SUCCESS              =  0, /**< Algorithm converged */
    BOBYQA_BAD_NVARS            = -1, /**< Bad number of variables */
    BOBYQA_BAD_NPT              = -2, /**< NPT is not in the required interval */
    BOBYQA_BAD_RHO_RANGE        = -3, /**< Bad trust region radius parameters */
    BOBYQA_BAD_SCALING          = -4, /**< Bad scaling factor(s) */
    BOBYQA_TOO_CLOSE            = -5, /**< Insufficient space between the bounds */
    BOBYQA_ROUNDING_ERRORS      = -6, /**< Too much cancellation in a denominator */
    BOBYQA_TOO_MANY_EVALUATIONS = -7, /**< Maximum number of function evaluations exceeded */
    BOBYQA_STEP_FAILED          = -8, /**< A trust region step has failed to reduce Q */
} bobyqa_status;

/* BOBYQA seeks the least value of a function of many variables, by applying a trust
   region method that forms quadratic models by interpolation. There is usually some
   freedom in the interpolation conditions, which is taken up by minimizing the Frobenius
   norm of the change to the second derivative of the model, beginning with the zero
   matrix. The values of the variables are constrained by upper and lower bounds. The
   arguments of the subroutine are as follows.

   N must be set to the number of variables and must be at least two. NPT is the number of
   interpolation conditions. Its value must be in the interval [N+2,(N+1)(N+2)/2]. Choices
   that exceed 2*N+1 are not recommended.

   OBJFUN is provided by the user to compute the objective function value at the values of
   the variables X(1),X(2),...,X(N), which are generated automatically by BOBYQA in a way
   that satisfies the bounds given in XL and XU. DATA is anything needed by the function
   and which is passed as is to OBJFUN by BOBYQA.

   Initial values of the variables must be set in X(1),X(2),...,X(N). They will be changed
   to the values that give the least calculated F. For I=1,2,...,N, XL(I) and XU(I) must
   provide the lower and upper bounds, respectively, on X(I). The construction of
   quadratic models requires XL(I) to be strictly less than XU(I) for each I. Further, the
   contribution to a model from changes to the I-th variable is damaged severely by
   rounding errors if XU(I)-XL(I) is too small.

   RHOBEG and RHOEND must be set to the initial and final values of a trust region radius,
   so both must be positive with RHOEND no greater than RHOBEG. Typically, RHOBEG should
   be about one tenth of the greatest expected change to a variable, while RHOEND should
   indicate the accuracy that is required in the final values of the variables. An error
   return occurs if any of the differences XU(I)-XL(I), I=1,...,N, is less than 2*RHOBEG.

   The value of IPRINT should be set to 0, 1, 2 or 3, which controls the amount of
   printing. Specifically, there is no output if IPRINT=0 and there is output only at the
   return if IPRINT=1. Otherwise, each new value of RHO is printed, with the best vector
   of variables so far and the corresponding value of the objective function. Further,
   each new value of F with its variables are output if IPRINT=3.

   MAXFUN must be set to an upper bound on the number of calls of OBJFUN.

   The array W will be used for working space. Its length must be at least
   (NPT+5)*(NPT+N)+3*N*(N+5)/2. Upon successful return, the first element of W will be set
   to the function value at the solution. */
extern bobyqa_status bobyqa(
    const ptrdiff_t n, const ptrdiff_t npt, bobyqa_objfun* objfun, void* data, double* x,
    const double* xl, const double* xu, const double rhobeg, const double rhoend,
    const ptrdiff_t iprint, const ptrdiff_t maxfun, double* w);

/*
  SCL is an array of N scaling factors or `NULL`.

  The array W will be used for working space. Its length must be at least
  (NPT+5)*(NPT+N)+3*N*(N+5)/2 if SCL is `NULL` and at least (NPT+5)*(NPT+N)+3*N*(N+5)/2 +
  3*N otherwise. Upon successful return, the first element of W will be set to the
  function value at the solution.
*/
extern bobyqa_status bobyqa_optimize(
    const ptrdiff_t n, const ptrdiff_t npt, bool maximize, bobyqa_objfun* objfun,
    void* data, double* x, const double* xl, const double* xu, const double* scl,
    const double rhobeg, const double rhoend, const ptrdiff_t iprint,
    const ptrdiff_t maxfun, double* w);

extern const char* bobyqa_reason(bobyqa_status status);

/* Test problem for BOBYQA, the objective function being the sum of the reciprocals of all
   pairwise distances between the points P_I, I=1,2,...,M in two dimensions, where M=N/2
   and where the components of P_I are X(2*I-1) and X(2*I). Thus each vector X of N
   variables defines the M points P_I. The initial X gives equally spaced points on a
   circle. Four different choices of the pairs (N,NPT) are tried, namely (10,16), (10,21),
   (20,26) and (20,41). Convergence to a local minimum that is not global occurs in both
   the N=10 cases. The details of the results are highly sensitive to computer rounding
   errors. The choice IPRINT=2 provides the current X and optimal F so far whenever RHO is
   reduced. The bound constraints of the problem require every component of X to be in the
   interval [-1,1]. */
extern void bobyqa_test(void);

/*--------------------------------------------------------------------------------------*/
/* FORTRAN SUPPORT */

/* Depending on your FORTRAN compiler, the names of the compiled functions may have to be
   modified. The various possibilities can be chosen via the macro FORTRAN_LINKAGE:

   -UFORTRAN_LINKAGE (or FORTRAN_LINKAGE undefined)
                         No support for FORTRAN will be compiled.

   -DFORTRAN_LINKAGE=0   FORTRAN link name is the same as with the C compiler.

   -DFORTRAN_LINKAGE=1   FORTRAN link name is is the function name in upper case letters
                         (for instance, `foo` yields `FOO`).

   -DFORTRAN_LINKAGE=2   FORTRAN link name is the function name suffixed with an
                         underscore (for instance, `foo` yields `foo_`).

   -DFORTRAN_LINKAGE=3   FORTRAN link name is the function name in upper case letters and
                         suffixed with an underscore (for instance, `foo` yields `FOO_`).
*/

#ifdef FORTRAN_LINKAGE

# if FORTRAN_LINKAGE == 0
#   define FORTRAN_NAME(a,A) a
#   error names will clash
# elif FORTRAN_LINKAGE == 1
#   define FORTRAN_NAME(a,A) A
# elif FORTRAN_LINKAGE == 2
#   define FORTRAN_NAME(a,A) a##_
# elif FORTRAN_LINKAGE == 3
#   define FORTRAN_NAME(a,A) A##_
# else
#   error unsupported FORTRAN linkage
# endif

/* This subroutine is a version of BOBYQA that is callable from FORTRAN code. The main
   difference with the C version is that the objective function must be provided by the
   following external subroutine.

   SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set F to the value of
   the objective function for the current values of the variables X(1),X(2),...,X(N),
   which are generated automatically in a way that satisfies the bounds given in XL and
   XU.
*/
extern int
FORTRAN_NAME(bobyqa,BOBYQA)(const ptrdiff_t* n, const ptrdiff_t* npt, double* x,
                            const double* xl, const double* xu, const double* rhobeg,
                            const double* rhoend, const ptrdiff_t* iprint,
                            const ptrdiff_t* maxfun, double* w);

/* Wrapper function to emulate `bobyqa_objfun` objective function calling the user-defined
   `calfun_` subroutine. */
extern double
bobyqa_calfun_wrapper(const ptrdiff_t n, const double* x, void* data);

/* Subroutine that must be defined by the application to use the FORTRAN wrapper to
   BOBYQA. */
extern int
FORTRAN_NAME(calfun,CALFUN)(const ptrdiff_t* n, double* x, double* f);

#endif /* FORTRAN_LINKAGE */

#ifdef __cplusplus
}
#endif

#endif /* BOBYQA_H_ */
