/*
 * newuoa.h -
 *
 * Definitions for Mike Powell's NEWUOA algorithm for minimizing a function of many
 * variables. The method is "derivatives free" (only the function values are needed). The
 * algorithm is described in:
 *
 *     M.J.D. Powell, "The NEWUOA software for unconstrained minimization without
 *     derivatives", in Large-Scale Nonlinear Optimization, editors G. Di Pillo and M.
 *     Roma, Springer (2006), pages 255-297.
 *
 * The present code is based on the original FORTRAN version written by Mike Powell who
 * kindly provides his code on demand (at mjdp@cam.ac.uk) and has been converted to C by
 * É. Thiébaut.
 *
 * Copyright (c) 2004, Mike Powell (FORTRAN version).
 * Copyright (c) 2015, Éric Thiébaut (C version).
 *
 * Read the accompanying `LICENSE` file for details.
 */

#ifndef NEWUOA_H_
#define  NEWUOA_H_ 1

#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Prototype of the objective function assumed by the NEWUOA routine. The returned value
   is the function value at X, the current variables, N is the number of variables and
   DATA is anything else needed by the objective function (unused by NEWUOA itself). */
typedef double newuoa_objfun(const ptrdiff_t n, const double *x, void* data);

/**
 * @brief Status for NEWUOA routines.
 *
 * This type enumerate the possible values returned by newuoa(),
 * newuoa_get_status() and newuoa_iterate().
 */
typedef enum {
    NEWUOA_INITIAL_ITERATE      =  2, /**< Only used internaly */
    NEWUOA_ITERATE              =  1, /**< Caller is requested to evaluate the objective
                                       *   function and call newuoa_iterate() */
    NEWUOA_SUCCESS              =  0, /**< Algorithm converged */
    NEWUOA_BAD_NVARS            = -1, /**< Bad number of variables */
    NEWUOA_BAD_NPT              = -2, /**< NPT is not in the required interval */
    NEWUOA_BAD_RHO_RANGE        = -3, /**< Invalid RHOBEG/RHOEND */
    NEWUOA_BAD_SCALING          = -4, /**< Bad scaling factor(s) */
    NEWUOA_ROUNDING_ERRORS      = -5, /**< Too much cancellation in a denominator */
    NEWUOA_TOO_MANY_EVALUATIONS = -6, /**< Maximum number of function evaluations
                                       *   exceeded */
    NEWUOA_STEP_FAILED          = -7, /**< Trust region step has failed to reduce
                                       *   quadratic approximation */
    NEWUOA_BAD_ADDRESS          = -8, /**< Illegal NULL address */
    NEWUOA_CORRUPTED            = -9, /**< Corrupted or misused workspace */
} newuoa_status;

/* This subroutine seeks the least value of a function of many variables, by a trust
   region method that forms quadratic models by interpolation. There can be some freedom
   in the interpolation conditions, which is taken up by minimizing the Frobenius norm of
   the change to the second derivative of the quadratic model, beginning with a zero
   matrix. The arguments of the subroutine are as follows.

   N must be set to the number of variables and must be at least two. NPT is the number of
   interpolation conditions. Its value must be in the interval [N+2,(N+1)(N+2)/2].

   Function OBJFUN(N,X,DATA) must be provided by the user to compute the value of the
   objective function for the variables X(1),X(2),...,X(N). DATA is anything else needed
   by the objective function.

   Initial values of the variables must be set in X(1),X(2),...,X(N). They will be changed
   to the values that give the least calculated F.

   RHOBEG and RHOEND must be set to the initial and final values of a trust region radius,
   so both must be positive with RHOEND<=RHOBEG. Typically RHOBEG should be about one
   tenth of the greatest expected change to a variable, and RHOEND should indicate the
   accuracy that is required in the final values of the variables.

   The value of IPRINT should be set to 0, 1, 2 or 3, which controls the amount of
   printing. Specifically, there is no output if IPRINT=0 and there is output only at the
   return if IPRINT=1. Otherwise, each new value of RHO is printed, with the best vector
   of variables so far and the corresponding value of the objective function. Further,
   each new value of F with its variables are output if IPRINT=3.

   MAXFUN must be set to an upper bound on the number of calls of OBJFUN.

   The array W will be used for working space. Its length must be at least
   (NPT+13)*(NPT+N)+3*N*(N+3)/2.

   The returned value should be NEWUOA_SUCCESS, but a different value can be returned upon
   error (see `newuoa_reason` for an explanatory message). */
extern newuoa_status newuoa(
    const ptrdiff_t n, const ptrdiff_t npt,
    newuoa_objfun* objfun, void* data,
    double* x, const double rhobeg, const double rhoend,
    const ptrdiff_t iprint, const ptrdiff_t maxfun,
    double* work);

/**
 * Optimize a function of many variables without derivatives.
 *
 * This function seeks the least (or the most) value of a function `f(x)` of many
 * variables `x[0]`, `x[1]`, ..., `x[n-1], by a trust region method that forms quadratic
 * models by interpolation. There can be some freedom in the interpolation conditions,
 * which is taken up by minimizing the Frobenius norm of the change to the second
 * derivative of the quadratic model, beginning with a zero matrix.
 *
 * Arguments `rhobeg` and `rhoend` must be set to the initial and final values of a trust
 * region radius, so both must be positive with `0 < rhoend <= rhobeg`. Typically `rhobeg`
 * should be about one tenth of the greatest expected change to a variable, and `rhoend`
 * should indicate the accuracy that is required in the final values of the variables. The
 * proper scaling of the variables is important for the success of the algorithm and the
 * optional `scale` argument should be specified if the typical precision is not the same
 * for all variables. If specified, `scale` is an array of same dimensions as `x` with
 * strictly nonnegative values, such that `scale[i]*rho` (with `rho` the trust region
 * radius) is the size of the trust region for the i-th variable. If `scale` is not
 * specified, a unit scaling for all the variables is assumed.
 *
 * @param n      - The number of variables which must be at least 2.
 *
 * @param npt    - The number of interpolation conditions. Its value must be in the
 *                 interval `[n + 2, (n + 1)*(n + 2)/2]`. The recommended value is `2*n +
 *                 1`.
 *
 * @param maximize - If true, maximize the function; otherwise, minimize it.
 *
 * @param objfun - The objective function. Called as `objfun(n,x,data)`, it returns the
 *                 value of the objective function for the variables `x[0]`, `x[1]`, ...,
 *                 `x[n-1]`. Argument `data` is anything else needed by the objective
 *                 function.
 *
 * @param data   - Anything needed by the objective function. This address is passed to
 *                 the objective fucntion each time it is called.
 *
 * @param x      - On entry, the initial variables; on return, the solution.
 *
 * @param scale  - Scaling factors for the variables. May be `NULL` to use the same unit
 *                 scaling factors for all variables. Otherwise, must all be strictly
 *                 positive.
 *
 * @param rhobeg - The initial radius of the trust region.
 *
 * @param rhoend - The final radius of the trust region.
 *
 * @param iprint - The amount of printing, its value should be set to 0, 1, 2 or 3.
 *                 Specifically, there is no output if `iprint = 0` and there is output
 *                 only at the return if `iprint = 1`. Otherwise, each new value of `rho`
 *                 is printed, with the best vector of variables so far and the
 *                 corresponding value of the objective function. Further, each new value
 *                 of `f(x)` with its variables are output if `iprint = 3`.
 *
 * @param maxfun - The maximum number of calls to the objective function.
 *
 * @param work   - A workspace array. If `scl` is `NULL`, its length must be at least
 *                 `(npt + 13)*(npt + n) + 3*n*(n + 3)/2` and at least `(npt + 13)*(npt +
 *                 n) + 3*n*(n + 3)/2 + n` otherwise. On exit, `work[0]` is set with
 *                 `f(x)` the value of the objective function at the solution.
 *
 * @return On success, the returned value is `NEWUOA_SUCCESS`; a different value is
 *         returned on error (see `newuoa_reason` for an explanatory message).
 */
extern newuoa_status newuoa_optimize(
    ptrdiff_t n, ptrdiff_t npt, bool maximize, newuoa_objfun* objfun, void* data,
    double x[], const double scl[], double rhobeg, double rhoend, ptrdiff_t iprint,
    ptrdiff_t maxfun, double* work);

/* Get a textual explanation of the status returned by NEWUOA. */
extern const char* newuoa_reason(newuoa_status status);

/*--------------------------------------------------------------------------------------*/
/* REVERSE COMMUNICATION VERSION */

/* Opaque structure used by the reverse communication version of NEWUOA. */
typedef struct newuoa_context_ newuoa_context;

/* Allocate a new reverse communication workspace for NEWUOA algorithm. The returned
   address is `NULL` to indicate an error: either invalid parameters (external variable
   `errno` set to `EINVAL`), or insufficient memory (external variable `errno` set to
   `ENOMEM`). The arguments correspond to those of `newuoa` (except that the variables X
   are omitted because they need not be specified until the first iteration with
   `newuoa_iterate` and that the workspace W is automatically allocated). The initial
   status of the created workspace is guaranteed to be `NEWUOA_ITERATE`.

   When no longer needed, the workspace must be deleted with `newuoa_delete`.

   A typical usage is:
   ```
   double x[N], f;
   newuoa_context* ctx;
   x[...] = ...; // initial solution
   ctx = newuoa_create(N, NPT, RHOBEG, RHOEND, IPRINT, MAXFUN);
   status = newuoa_get_status(ctx);
   while (status == NEWUOA_ITERATE) {
     f = ...; // compute function value at X
     status = newuoa_iterate(ctx, f, x);
   }
   newuoa_delete(ctx);
   if (status != NEWUOA_SUCCESS) {
     fprintf(stderr, "Something wrong occured in NEWUOA: %s\n",
             newuoa_reason(status));
   }
   ```
 */
extern newuoa_context*
newuoa_create(const ptrdiff_t n, const ptrdiff_t npt, const double rhobeg,
              const double rhoend, const ptrdiff_t iprint, const ptrdiff_t maxfun);

/* Release resource allocated for NEWUOA reverse communication workspace. Argument can be
   `NULL`. */
extern void newuoa_delete(newuoa_context* ctx);

/* Perform the next iteration of the reverse communication version of the NEWUOA
   algorithm. On entry, the workspace status must be `NEWUOA_ITERATE`, `f` is the function
   value at `x`. On exit, the returned value (the new workspace status) is:
   `NEWUOA_ITERATE` if a new trial point has been stored in `x` and if user is requested
   to compute the function value for the new point; `NEWUOA_SUCCESS` if algorithm has
   converged; anything else indicates an error (see `newuoa_reason` for an explanatory
   message). */
extern newuoa_status newuoa_iterate(newuoa_context* ctx, double f, double* x);

/* Restart NEWUOA algorithm using the same parameters. The return value is the new status
   of the algorithm, see `newuoa_get_status` for details. */
extern newuoa_status newuoa_restart(newuoa_context* ctx);

/* Get the current status of the algorithm. Result is: `NEWUOA_ITERATE` if user is
   requested to compute F(X); `NEWUOA_SUCCESS` if algorithm has converged; anything else
   indicates an error (see `newuoa_reason` for an explanatory message). */
extern newuoa_status newuoa_get_status(const newuoa_context* ctx);

/* Get the current number of function evaluations. Result is -1 if something is wrong
   (e.g. CTX is NULL), nonnegative otherwise. */
extern ptrdiff_t newuoa_get_nevals(const newuoa_context* ctx);

/* Get the current size of the trust region. Result is 0 if algorithm not yet started
   (before first iteration), -1 if something is wrong (e.g. CTX is NULL), strictly
   positive otherwise. */
extern double newuoa_get_rho(const newuoa_context* ctx);

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

     -DFORTRAN_LINKAGE=3   FORTRAN link name is the function name in upper case letters
                           and suffixed with an underscore (for instance, `foo` yields
                           `FOO_`).
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

/* This subroutine is a version of NEWUOA that is callable from FORTRAN code. The main
   difference with the C version is that the objective function must be provided by the
   following external subroutine.

   SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set F to the value of
   the objective function for the current values of the variables X(1),X(2),...,X(N),
   which are generated automatically by NEWUOA. */
extern int
FORTRAN_NAME(newuoa,NEWUOA)(const ptrdiff_t* n, const ptrdiff_t* npt, double* x,
                            const double* rhobeg, const double* rhoend,
                            const ptrdiff_t* iprint, const ptrdiff_t* maxfun, double* w);

/* Wrapper function to emulate `newuoa_calfun` function calling
   the user-defined `calfun_` subroutine. */
extern double newuoa_calfun_wrapper(const ptrdiff_t n, const double* x, void* data);

/* Subroutine that must be defined by the application to use the FORTRAN
   wrapper to NEWUOA. */
extern int calfun_(const ptrdiff_t* n, double *x, double *f);

#endif /* FORTRAN_LINKAGE */

/*--------------------------------------------------------------------------------------*/
/* TESTS */

/* Test problem for NEWUOA, the objective function being the sum of the reciprocals of all
   pairwise distances between the points P_I, I=1,2,...,M in two dimensions, where M=N/2
   and where the components of P_I are X(2*I-1) and X(2*I). Thus each vector X of N
   variables defines the M points P_I. The initial X gives equally spaced points on a
   circle. Four different choices of the pairs (N,NPT) are tried, namely (10,16), (10,21),
   (20,26) and (20,41). Convergence to a local minimum that is not global occurs in both
   the N=10 cases. The details of the results are highly sensitive to computer rounding
   errors. The choice IPRINT=2 provides the current X and optimal F so far whenever RHO is
   reduced. The bound constraints of the problem require every component of X to be in the
   interval [-1,1]. */
extern void newuoa_test(void);

#ifdef __cplusplus
}
#endif

#endif /* NEWUOA_H_ */
