/*
 * cobyla.h -
 *
 * Definitions for Mike Powell's COBYLA algorithm for minimizing a function of
 * a few variables.  The method is "derivatives free" (only the function values
 * are needed) and accounts for constraints on the variables.  The algorithm is
 * described in:
 *
 *   M.J.D. Powell, "A direct search optimization method that models the
 *   objective and constraint functions by linear interpolation," in Advances
 *   in Optimization and Numerical Analysis Mathematics and Its Applications,
 *   vol. 275 (eds. Susana Gomez and Jean-Pierre Hennart), Kluwer Academic
 *   Publishers, pp. 51-67 (1994).
 *
 * The present code is based on the original FORTRAN version written by Mike
 * Powell who kindly provides his code on demand and has been converted to C by
 * É. Thiébaut.
 *
 * Copyright (c) 1992, Mike Powell (FORTRAN version).
 * Copyright (c) 2015, Éric Thiébaut (C version).
 */

#ifndef COBYLA_H_
#define  COBYLA_H_ 1

#include "optimpack.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Prototype of the function assumed by the COBYLA algorithm.
 *
 * Prototype of the objective function assumed by the COBYLA routine.  The
 * returned value is the function value at `x`, `n` is the number of variables,
 * `m` is the number of constraints, `x` are the current values of the
 * variables and `con` is to store the `m` constraints.  `data` is anything
 * needed by the function (unused by COBYLA itself).
 */
typedef double
cobyla_calcfc(opk_index n, opk_index m, const double x[],
              double con[], void* data);

/**
 * @brief Status for COBYLA routines.
 *
 * This type enumerate the possible values returned by cobyla(),
 * cobyla_get_status() and cobyla_iterate().
 */
typedef enum {
    COBYLA_INITIAL_ITERATE      =  2, /**< Only used internally */
    COBYLA_ITERATE              =  1, /**< User requested to compute
                                       *   F(X) and C(X) */
    COBYLA_SUCCESS              =  0, /**< Algorithm converged */
    COBYLA_BAD_NVARS            = -1, /**< Bad number of variables */
    COBYLA_BAD_NCONS            = -2, /**< Bad number of constraints */
    COBYLA_BAD_RHO_RANGE        = -3, /**< Invalid trust region parameters */
    COBYLA_BAD_SCALING          = -4, /**< Bad scaling factor(s) */
    COBYLA_ROUNDING_ERRORS      = -5, /**< Rounding errors prevent progress */
    COBYLA_TOO_MANY_EVALUATIONS = -6, /**< Too many evaluations */
    COBYLA_BAD_ADDRESS          = -7, /**< Illegal address */
    COBYLA_CORRUPTED            = -8, /**< Corrupted workspace */
} cobyla_status;

/**
 * Minimize a function of many variables subject to inequality constraints.
 *
 * The `cobyla` algorithm minimizes an objective function `f(x)` subject to `m`
 * inequality constraints on `x`, where `x` is a vector of variables that has
 * `n` components.  The algorithm employs linear approximations to the
 * objective and constraint functions, the approximations being formed by
 * linear interpolation at `n+1` points in the space of the variables.  We
 * regard these interpolation points as vertices of a simplex.  The parameter
 * `rho` controls the size of the simplex and it is reduced automatically from
 * `rhobeg` to `rhoend`.  For each `rho` the subroutine tries to achieve a good
 * vector of variables for the current size, and then `rho` is reduced until
 * the value `rhoend` is reached.  Therefore `rhobeg` and `rhoend` should be
 * set to reasonable initial changes to and the required accuracy in the
 * variables respectively, but this accuracy should be viewed as a subject for
 * experimentation because it is not guaranteed.  The subroutine has an
 * advantage over many of its competitors, however, which is that it treats
 * each constraint individually when calculating a change to the variables,
 * instead of lumping the constraints together into a single penalty function.
 * The name of the subroutine is derived from the phrase "Constrained
 * Optimization BY Linear Approximations".
 *
 * The user must set the values of `n`, `m`, `rhobeg` and `rhoend`, and must
 * provide an initial vector of variables in `x`.  Further, the value of
 * `iprint` should be set to 0, 1, 2 or 3, which controls the amount of
 * printing during the calculation.  Specifically, there is no output if
 * `iprint=0` and there is output only at the end of the calculation if
 * `iprint=1`.  Otherwise each new value of `rho` and `sigma` is printed.
 * Further, the vector of variables and some function information are given
 * either when `rho` is reduced or when each new value of `f(x)` is computed in
 * the cases `iprint=2` or `iprint=3` respectively.  Here `sigma` is a penalty
 * parameter, it being assumed that a change to `x` is an improvement if it
 * reduces the merit function:
 *
 *     f(x) + sigma*max(0.0, -C1(x), -C2(x), ..., -CM(x)),
 *
 * where `C1`, `C2`, ..., `CM` denote the constraint functions that should
 * become nonnegative eventually, at least to the precision of `rhoend`.  In
 * the printed output the displayed term that is multiplied by `sigma` is
 * called `maxcv`, which stands for "MAXimum Constraint Violation".  The
 * argument `maxfun` is the address of an integer variable that must be set by
 * the user to a limit on the number of calls of `fc`, the purpose of this
 * routine being given below.  The value of `maxfun` will be altered to the
 * number of calls of `fc` that are made.  The arguments `work` and `iact`
 * provide real and integer arrays that are used as working space.  Their
 * lengths must be at least `n*(3*n+2*m+11)+4*m+6` and `m+1` respectively.
 *
 * In order to define the objective and constraint functions, we require a
 * function `fc` that has the following prototype:
 *
 *     REAL fc(INTEGER n, INTEGER m, const REAL x[], REAL con[], void* data);
 *
 * The values of `n` and `m` are fixed and have been defined already, while `x`
 * is now the current vector of variables.  The function should return the
 * value of the objective function at `x` and store constraint functions at `x`
 * in `con[0]`, `con[1]`, ..., `con[m-1]`.  Argument `data` will be set with
 * whatever value has been provided when `cobyla` was called.  Note that we are
 * trying to adjust `x` so that `f(x)` is as small as possible subject to the
 * constraint functions being nonnegative.
 *
 * @param n      - The number of variables.
 *
 * @param m      - The number of constraints.
 *
 * @param fc     - The objective function.
 *
 * @param data   - Anything needed by the objective function.
 *
 * @param x      - On entry, the initial variables; on exit, the final
 *                 variables.
 *
 * @param rhobeg - The initial trust region radius.
 *
 * @param rhoend - The final trust region radius.
 *
 * @param maxfun - On entry, the maximum number of calls to `fc`; on exit, the
 *                 actual number of calls to `fc`.
 *
 * @param iprint - The level of verbosity.
 *
 * @param maxfun - The maximum number of calls to `fc`.
 *
 * @param work   - Workspace array with at least `n*(3*n+2*m+11)+4*m+6`
 *                 elements.  On successful exit, the value of the objective
 *                 function and of the worst constraint at the final `x` are
 *                 stored in `work[0]` and `work[1]` respectively.
 *
 * @param iact   - Workspace array with at least `m+1` elements.  On successful
 *                 exit, the actual number of calls to `fc` is stored in
 *                 `iact[0]`.
 *
 * @return `COBYLA_SUCCESS` is returned when the algorithm is successful; any
 *         other value indicates an error (use `cobyla_reason` to have an
 *         explanation).
 */
extern cobyla_status cobyla(
    opk_index n, opk_index m,
    cobyla_calcfc* fc, void* data,
    double x[], double rhobeg, double rhoend,
    opk_index iprint, opk_index maxfun,
    double work[], opk_index iact[]);

/**
 * Minimize or maximize a function of many variables subject to inequality
 * constraints with optional scaling.
 *
 * This function is a variant of `cobyla` which attempts to minimize or
 * maximize an objective function of many variables subject to inequality
 * constraints.  The scaling of the variables is important for the success of
 * the algorithm and the `scale` argument (if not `NULL`) should specify the
 * relative size of each variable.  If specified, `scale` is an array of `n`
 * strictly nonnegative values, such that `scale[i]*rho` (with `rho` the trust
 * region radius) is the size of the trust region for the `i`-th variable.
 * Thus `scale[i]*rhobeg` is the typical step size for the `i`-th variable at
 * the beginning of the algorithm and `scale[i]*rhoend` is the typical
 * precision on the `i`-th variable at the end.  If `scale` is not specified, a
 * unit scaling for all the variables is assumed.
 *
 * @param n      - The number of variables.
 *
 * @param m      - The number of constraints.
 *
 * @param maximize - Non-zero to attempt to maximize the objective function;
 *                 zero to attempt to minimize it.
 *
 * @param fc     - The objective function.
 *
 * @param data   - Anything needed by the objective function.
 *
 * @param x      - On entry, the initial variables; on exit, the final
 *                 variables.
 *
 * @param scale  - An array of `n` scaling factors, all strictly positive.  Can
 *                 be `NULL` to perform no scaling of the variables.
 *
 * @param rhobeg - The initial trust region radius.
 *
 * @param rhoend - The final trust region radius.
 *
 * @param maxfun - The maximum number of calls to `fc`.
 *
 * @param iprint - The level of verbosity.
 *
 * @param maxfun - On entry, the maximum number of calls to `fc`; on exit, the
 *                 actual number of calls to `fc`.
 *
 * @param work   - Workspace array with at least `n*(3*n+2*m+11)+4*m+6`
 *                 elements.  On successful exit, the value of the objective
 *                 function and of the worst constraint at the final `x` are
 *                 stored in `work[0]` and `work[1]` respectively.
 *
 * @param iact   - Workspace array with at least `m+1` elements.  On successful
 *                 exit, the actual number of calls to `fc` is stored in
 *                 `iact[0]`.
 *
 * @return `COBYLA_SUCCESS` is returned when the algorithm is successful; any
 *         other value indicates an error (use `cobyla_reason` to have an
 *         explanation).
 */
extern cobyla_status cobyla_optimize(
    opk_index n, opk_index m,
    opk_bool maximize, cobyla_calcfc* fc, void* data,
    double x[], const double scl[], double rhobeg, double rhoend,
    opk_index iprint, opk_index maxfun,
    double work[], opk_index iact[]);

/* Opaque structure used by the reverse communication variant of COBYLA. */
typedef struct cobyla_context_ cobyla_context;

/* Allocate a new reverse communication workspace for COBYLA algorithm.  The
   returned address is `NULL` to indicate an error due to either invalid
   parameters (external variable `errno` set to `EINVAL`), or insufficient
   memory (external variable `errno` set to `ENOMEM`).  When no longer needed,
   the workspace must be deleted with `cobyla_delete`.

   A typical usage is:
   ```
   double x[N], c[M], f;
   cobyla_context* ctx;
   x[...] = ...; // initial solution
   ctx = cobyla_create(N, M, RHOBEG, RHOEND, IPRINT, MAXFUN);
   status = cobyla_get_status(ctx);
   while (status == COBYLA_ITERATE) {
     f = ...; // compute function value at X
     c[...] = ...; // compute constraints at X
     status = cobyla_iterate(ctx, f, x, c);
   }
   cobyla_delete(ctx);
   if (status != COBYLA_SUCCESS) {
     fprintf(stderr, "Something wrong occured in COBYLA: %s\n",
             cobyla_reason(status));
   }
   ```
 */
extern cobyla_context*
cobyla_create(opk_index n, opk_index m, double rhobeg, double rhoend,
              opk_index iprint, opk_index maxfun);

/* Release ressources allocated for COBYLA reverse communication workspace.
   Argument can be `NULL`. */
extern void
cobyla_delete(cobyla_context* ctx);

/* Perform the next iteration of the reverse communication variant of the
   COBYLA algorithm.  On entry, the workspace status must be `COBYLA_ITERATE`,
   `f` and `c` are the function value and the constraints at `x`.  On exit, the
   returned value (the new workspace status) is: `COBYLA_ITERATE` if a new
   trial point has been stored in `x` and if user is requested to compute the
   function value and the constraints at the new point; `COBYLA_SUCCESS` if
   algorithm has converged and `x` has been set with the variables at the
   solution (the corresponding function value can be retrieved with
   `cobyla_get_last_f`); anything else indicates an error (see `cobyla_reason`
   for an explanatory message). */
extern cobyla_status cobyla_iterate(
    cobyla_context* ctx, double f, double x[], double c[]);

/* Restart COBYLA algorithm using the same parameters.  The return value is
   the new status of the algorithm, see `cobyla_get_status` for details. */
extern cobyla_status cobyla_restart(cobyla_context* ctx);

/* Get the current status of the algorithm.  Result is: `COBYLA_ITERATE` if
   user is requested to compute F(X) and C(X); `COBYLA_SUCCESS` if algorithm
   has converged; anything else indicates an error (see `cobyla_reason` for an
   explanatory message). */
extern cobyla_status cobyla_get_status(const cobyla_context* ctx);

/* Get the current number of function evaluations.  Result is -1 if something
   is wrong (e.g. CTX is NULL), nonnegative otherwise. */
extern opk_index
cobyla_get_nevals(const cobyla_context* ctx);

/* Get the current size of the trust region.  Result is 0 if algorithm has not
   yet been started (before first iteration), -1 if something is wrong
   (e.g. CTX is NULL), strictly positive otherwise. */
extern double
cobyla_get_rho(const cobyla_context* ctx);

/* Get the last function value.  Upon convergence of `cobyla_iterate`
   (i.e. return with status `COBYLA_SUCCESS`), this value corresponds to the
   function at the solution; otherwise, this value corresponds to the previous
   set of variables. */
extern double
cobyla_get_last_f(const cobyla_context* ctx);

/* Get a textual explanation of the status returned by `cobyla`,
   `cobyla_get_status` and `cobyla_iterate`. */
extern const char* cobyla_reason(cobyla_status status);

/*---------------------------------------------------------------------------*/
/* FORTRAN SUPPORT */

/* Depending on your FORTRAN compiler, the names of the compiled functions
   may have to be modified.  The various possibilities can be chosen via the
   macro FORTRAN_LINKAGE:

     -UFORTRAN_LINKAGE  (or FORTRAN_LINKAGE undefined)
           No support for FORTRAN will be compiled.

     -DFORTRAN_LINKAGE=0   FORTRAN link name is the same as with the C
                           compiler.

     -DFORTRAN_LINKAGE=1   FORTRAN link name is is the function name in upper
                           case letters (for instance, `foo` yields `FOO`).

     -DFORTRAN_LINKAGE=2   FORTRAN link name is the function name suffixed
                           with an underscore (for instance, `foo` yields
                           `foo_`).

     -DFORTRAN_LINKAGE=3   FORTRAN link name is the function name in upper
                           case letters and suffixed with an underscore
                           (for instance, `foo` yields `FOO_`).
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

/*
 * Subroutine variant designed to be callable from FORTRAN code.  Usage is
 * similar to that of `cobyla` except that everything is passed by address and
 * that, in order to define the objective and constraint functions, we require
 * a subroutine that has the name and arguments:
 *
 *            SUBROUTINE CALCFC (N,M,X,F,CON)
 *            DIMENSION X(*),CON(*)
 *
 * The values of N and M are fixed and have been defined already, while X is
 * now the current vector of variables.  The subroutine should return the
 * objective and constraint functions at X in F and CON(1),CON(2), ...,CON(M).
 * Note that we are trying to adjust X so that F(X) is as small as possible
 * subject to the constraint functions being nonnegative.
 */
extern int
FORTRAN_NAME(cobyla,COBYLA)(opk_index* n, opk_index* m, double x[],
                            double* rhobeg, double* rhoend,
                            opk_index* iprint, opk_index* maxfun,
                            double w[], opk_index iact[]);

#endif /* FORTRAN_LINKAGE */

#ifdef __cplusplus
}
#endif

#endif /* COBYLA_H_ */
