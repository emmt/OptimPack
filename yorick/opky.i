/*
 * opky.i --
 *
 * OPKY is OptimPack for Yorick.
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of OptimPack (https://github.com/emmt/OptimPack).
 *
 * Copyright (C) 2014, 2015 Éric Thiébaut
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

if (is_func(plug_in)) plug_in, "opky";

local OPK_NLCG_FLETCHER_REEVES, OPK_NLCG_HESTENES_STIEFEL;
local OPK_NLCG_POLAK_RIBIERE_POLYAK, OPK_NLCG_FLETCHER;
local OPK_NLCG_LIU_STOREY, OPK_NLCG_DAI_YUAN;
local OPK_NLCG_PERRY_SHANNO, OPK_NLCG_HAGER_ZHANG;
local OPK_NLCG_POWELL, OPK_NLCG_SHANNO_PHUA;
local OPK_NLCG_DEFAULT;
extern opk_nlcg;
/* DOCUMENT opt = opk_nlcg(dims, flags=..., single=...);

     The function opk_nlcg() creates a new instance of OPKY reverse
     communication optimizer implementing a non-linear conjugate gradient
     method.  DIMS is the dimension list of the variables of the problem.

     Keyword FLAGS may be set to specify a specific variant of the non-linear
     conjugate gradient method (see below for more details).

     Keyword SINGLE may be set true to use single precision floating point
     variables.  The default is to use double precision floating point
     variables.

     The dot notation can be used to query some members of the optimizer
     instance:
         opt.flags       - the name of the optimization method;
         opt.size        - the size of the problem;
         opt.dims        - the dimensions of the problem;
         opt.single      - single precision?
         opt.task        - current pending task;
         opt.status      - last optimizer status;
         opt.reason      - textual explanation about last failure;
         opt.iterations  - number of iterations;
         opt.evaluations - number of function (and gradient) evaluations;
         opt.restarts    - number of restarts;

     The FLAGS keyword can be set with one of the following rules to compute
     the search direction:

         OPK_NLCG_FLETCHER_REEVES
         OPK_NLCG_HESTENES_STIEFEL
         OPK_NLCG_POLAK_RIBIERE_POLYAK
         OPK_NLCG_FLETCHER
         OPK_NLCG_LIU_STOREY
         OPK_NLCG_DAI_YUAN
         OPK_NLCG_PERRY_SHANNO
         OPK_NLCG_HAGER_ZHANG

     The rule can be combined (bitwise or'ed) with:

         OPK_NLCG_POWELL

     to force beta >= 0 (according to Powell's prescription) and with one of:

         OPK_NLCG_SHANNO_PHUA
         OPK_NLCG_OREN_SPEDICATO
         OPK_NLCG_BARZILAI_BORWEIN

     to specify how to compute the initial step size from the previous
     iteration.


   SEE ALSO: opk_iterate, opk_task, opk_start, opk_vmlmb.
 */

local OPK_EMULATE_BLMVM;
extern opk_vmlmb;
/* DOCUMENT opt = opk_vmlmb(dims, mem=..., lower=..., upper=..., flags=...,
                                  single=...);

     The function opk_vmlmb() creates a new instance of OPKY reverse
     communication optimizer implementing VMLM-B algorithm.  This algorithm is
     a limited memory variant of the BFGS variable metric (quasi-Newton) method
     for solving an optimization problem with optional bound constraints.  DIMS
     is the dimension list of the variables of the problem.

     Keyword MEM can be set to specify the number of previous steps to
     memorize.  By defalust MEM=5 is used

     Keywords LOWER and UPPER may be used to set lower and/or upper bounds for
     the variables.  Their value can be nothing (no bound), a scalar (same
     bound value for all variables), or an array (to specify a different bound
     value for each variables).

     Keyword SINGLE may be set true to use single precision floating point
     variables.  The default is to use double precision floating point
     variables.

     The dot notation can be used to query some members of the optimizer
     instance:
         opt.size        - the size of the problem;
         opt.dims        - the dimensions of the problem;
         opt.single      - single precision?
         opt.task        - current pending task;
         opt.status      - last optimizer status;
         opt.reason      - textual explanation about last failure;
         opt.iterations  - number of iterations;
         opt.evaluations - number of function (and gradient) evaluations;
         opt.restarts    - number of restarts;
         opt.projections - number of projections;

   SEE ALSO: opk_iterate, opk_task, opk_start, opk_nlcg.
 */


extern opk_iterate;
/* DOCUMENT task = opk_iterate(opt, x, fx, gx);
     Proceed with next iteration for the optimizer OPT.  X stores the current
     variables, FX is the function value at X and GX is the gradient of the
     function at X.  For constrained optimization, D is an additional work
     array.  See opk_task() for the interpretation of the returned value.

   SEE ALSO: opk_nlcg, opk_vmlmb, opk_task, opk_start.
 */

extern opk_start;
/* DOCUMENT task = opk_start(opt, x);
     Start or re-start the reverse communication optimizer OPT with initial
     variables X.  See opk_task() for the interpretaion of the returned value.

   SEE ALSO: opk_nlcg, opk_vmlmb, opk_task, opk_iterate.
*/

local OPK_TASK_ERROR, OPK_TASK_WARNING;
local OPK_TASK_COMPUTE_FG;
local OPK_TASK_NEW_X, OPK_TASK_FINAL_X;
extern opk_get_task;
/* DOCUMENT task = opk_get_task(opt);
     Query the current pending task for the reverse communication optimizer
     OPT.  The possible values for the returned value are:

        OPK_TASK_ERROR      - an error has occured;
        OPK_TASK_COMPUTE_FG - caller must compute the function value and its
                              gradient for the current variable and call
                              opk_iterate again;
        OPK_TASK_NEW_X      - new improved variables are available for
                              examination before calling opk_iterate again;
        OPK_TASK_FINAL_X    - the method has converged;
        OPK_TASK_WARNING    - the method has terminated with a warning;

   SEE ALSO: opk_nlcg, opk_vmlmb, opk_iterate, opk_start.
*/

extern opk_get_status;
extern opk_get_reason;
/* DOCUMENT status = opk_get_status(opt);
         or reason = opk_get_reason(status);

     Query the last status of optimizer OPT and the textual reason for a given
     status value.


   SEE ALSO: opk_nlcg, opk_vmlmb, opk_iterate, opk_start.
*/

func opk_minimize(fg, x0, mem=, nlcg=, flags=, single=, verb=, lower=, upper=)
/* DOCUMENT x = opk_minimize(fg, x0);

     This driver minimizes a smooth multi-variate function.  FG is a function
     which takes two arguments X and GX and which, for the given variables X,
     stores the gradient of the function in GX and returns the value of the
     function:

        fx = fg(x, gx);

     X0 are the initial variables, they are left unchanged.

     Keyword SINGLE may be set true to use single precision floating point
     variables.  The default is to use double precision floating point
     variables.  Beware that this have an incidence on the allocated array to
     store the gradient and may result in costly conversions if the function FG
     is not designed to work at the assumed precision.

     Keyword NLCG can be set with a non-zero value indicating which non-linear
     conjugate gradient method to use.  By default, a variable metric method is
     used.

   SEE ALSO: opk_nlcg, opk_vmlmb.
 */
{
  TRUE = 1n;
  FALSE = 0n;
  if (identof(x0) > Y_DOUBLE) {
    error, "bad data type for X0";
  }
  if (nlcg && !(is_void(mem) && is_void(lower) && is_void(upper))) {
    error, "keywords MEM, LOWER and UPPER cannot be used for a non-linear conjugate gradient method";
  }
  if (single) {
    type = float;
    single = TRUE;
  } else {
    type = double;
    single = FALSE;
  }
  x = type(unref(x0));
  gx = array(type, dimsof(x));
  n = numberof(x);
  if (nlcg) {
    opt = opk_nlcg(n, flags=flags, single=single);
  } else {
    opt = opk_vmlmb(n, mem=mem, single=single, lower=lower, upper=upper,
                    flags=flags);
  }
  task = opk_start(opt, x);
  while (TRUE) {
    if (task == OPK_TASK_COMPUTE_FG) {
      fx = fg(x, gx);
    } else if (task == OPK_TASK_NEW_X) {
      if (verb) {
        write, format="%4d  %4d  %.16E\n",
          opt.iterations, opt.evaluations, fx;
      }
    } else {
      break;
    }
    task = opk_iterate(opt, x, fx, gx);
  }
  if (task != OPK_TASK_FINAL_X) {
    if (task == OPK_TASK_WARNING) {
      write, format="WARNING: %s\n", "some warning occured";
    } else if (task == OPK_TASK_ERROR) {
      error, "some error occured";
    } else {
      error, "unexpected task";
    }
  }
  return x;
}

extern opk_init;
/* DOCUMENT opk_init;
     This procedure initializes internals of OPKY.  It can safely be called to
     re-initialize and restore values of global variables.  It is automatically
     called when the plug-in is loaded.
 */
opk_init;

/*---------------------------------------------------------------------------*/
/* Yorick interface to Mike Powell's COBYLA algorithm. */

local cobyla_minimize, cobyla_maximize;
/* DOCUMENT xmin = cobyla_minimize(f, c, x0, rhobeg, rhoend);
         or  obj = cobyla_minimize(f, c, x0, rhobeg, rhoend, all=1);
         or xmax = cobyla_maximize(f, c, x0, rhobeg, rhoend);
         or  obj = cobyla_maximize(f, c, x0, rhobeg, rhoend, all=1);

     Minimize or maximize the multi-variate function F under the inequality
     constraints implemented by C and starting at the initial point X0.  RHOBEG
     and RHOEND are the initial and final values of the trust region radius (0
     < RHOEND <= RHOBEG).

     The objective is to find X which solves one of the problems:

          min f(x)    s.t.   c(x) <= 0   (cobyla_minimize)
          max f(x)    s.t.   c(x) <= 0   (cobyla_maximize)

     Arguments F and C are user defined functions which take a single argument,
     the variables X, and return the function value and the constraints
     respectively.  If there are no contraints, C can be empty.

     Note that the proper scaling of the variables is important for the success
     of the algorithm.  RHOBEG should be set to the typical size of the region
     to explorate and RHOEND should be set to the typical precision.

     Keyword NPT sets the number of of interpolation conditions.  Its default
     value is equal to 2*N+1 (the recommended value).

     Keyword MAXFUN sets the maximum number of function evaluations.  Its
     default value is equal to 30*N.

     If keyword ALL is true, the result is a structured object.  For
     `cobyla_minimize`, the members of the returned object are:

        obj.fmin   = minimal function value found
        obj.cmin   = constraints at the minimum
        obj.xmin   = corresponding parameters
        obj.nevals = number of function evaluations
        obj.status = status of the algorithm upon return
        obj.rho    = final radius of the trust region

     For `cobyla_maximize`, the two first members are:

        obj.fmax   = minimal function value found
        obj.cmax   = constraints at the maximum
        obj.xmax   = corresponding parameters

     Keyword ERR sets the behavior in case of abnormal termination.  If ERR=0,
     anything but a success throws an error (this is also the default
     behavior); if ERR > 0, non-fatal errors are reported by a warning message;
     if ERR < 0, non-fatal errors are silently ignored.

     Keyword VERB set the verbosity level.


   SEE ALSO: cobyla_create, cobyla_error.
 */

func cobyla_minimize(f, c, x0, rhobeg, rhoend, npt=, maxfun=, all=, verb=, err=)
{
  local cx, fx, fmin, cmin, xmin;
  x = double(unref(x0));
  n = numberof(x);
  constrained = (! is_void(c));
  if (constrained) cx = c(x);
  ctx = cobyla_create(n, numberof(cx), rhobeg, rhoend,
                      (is_void(verb) ? 0 : verb),
                      (is_void(maxfun) ? 50*n : maxfun));
  start = 1n;
  do {
    if (constrained && ! start) cx = c(x);
    fx = f(x);
    if (start || fmin > fx) {
      fmin = fx;
      cmin = cx;
      xmin = x;
      start = 0n;
    }
    status = cobyla_iterate(ctx, fx, x, cx);
  } while (status == COBYLA_ITERATE);
  if (status != COBYLA_SUCCESS) cobyla_error, status, err;
  return (all ? save(fmin, cmin, xmin, nevals = ctx.nevals,
                     status, rho = ctx.rho) : x);
}

func cobyla_maximize(f, c, x0, rhobeg, rhoend, npt=, maxfun=, all=, verb=, err=)
{
  local cx, fx, fmax, cmax, xmax;
  x = double(unref(x0));
  n = numberof(x);
  constrained = (! is_void(c));
  if (constrained) cx = c(x);
  ctx = cobyla_create(n, numberof(cx), rhobeg, rhoend,
                      (is_void(verb) ? 0 : verb),
                      (is_void(maxfun) ? 50*n : maxfun));
  start = 1n;
  do {
    if (constrained && ! start) cx = c(x);
    fx = f(x);
    if (start || fmax < fx) {
      fmax = fx;
      cmax = cx;
      xmax = x;
      start = 0n;
    }
    status = cobyla_iterate(ctx, -fx, x, cx);
  } while (status == COBYLA_ITERATE);
  if (status != COBYLA_SUCCESS) cobyla_error, status, err;
  return (all ? save(fmax, cmax, xmax, nevals = ctx.nevals,
                     status, rho = ctx.rho) : x);
}

func cobyla_error(status, errmode)
/* DOCUMENT cobyla_error, status;
         or cobyla_error, status, errmode;

     Report an error in COBYLA according to the value of STATUS.  Nothing is
     done if STATUS is COBYLA_SUCCESS; otherwise, the optional argument ERRMODE
     determines the behavior.  If ERRMODE = 0, the routine throws an error
     (this is also the default behavior); if ERRMODE > 0, non-fatal errors are
     reported by a warning message; if ERRMODE < 0, non-fatal errors are
     silently ignored.

   SEE ALSO: cobyla_reason, error.
 */
{
  if (status != COBYLA_SUCCESS) {
    if (errmode && (status == COBYLA_ROUNDING_ERRORS ||
                    status == COBYLA_TOO_MANY_EVALUATIONS ||
                    status == COBYLA_STEP_FAILED)) {
      if (errmode > 0) {
        write, format="WARNING: Something wrong occured in COBYLA: %s",
          cobyla_reason(status);
      }
    } else {
      error, swrite(format="Something wrong occured in COBYLA: %s",
                    cobyla_reason(status));
    }
  }
}
errs2caller, cobyla_error;

local COBYLA_SUCCESS, COBYLA_ITERATE, COBYLA_ROUNDING_ERRORS;
local COBYLA_TOO_MANY_EVALUATIONS, COBYLA_BAD_ADDRESS;
local COBYLA_CORRUPTED;
extern cobyla_create;
extern cobyla_iterate;
/* DOCUMENT ctx = cobyla_create(n, m, rhobeg, rhoend, iprint, maxfun);
         or status = cobyla_iterate(ctx, f, x, c);

     The function `cobyla_create` makes a new instance for Mike Powell's COBYLA
     algorithm for minimizing a function of a few variables.  The method is
     "derivatives free" (only the function values are needed) and accounts for
     inequality constraints on the variables.  N is the number of variables, M
     is the number of constraints, RHOBEG and RHOEND are the initial and final
     size of the trust region, IPRINT control the verbosity of the method (see
     below) and MAXFUN is the maximum number of function evaluations.

     The objective is to find X which solves the problem:

          min f(x)    s.t.   c(x) <= 0

     Given a penalty parameter SIGMA > 0, COBYLA assumes that a change to X is
     an improvement if it reduces the merit function:

          f(x) + sigma*max(0.0,-min(c(x))),

     where c(x) are the constraint functions that should become nonnegative
     eventually, at least to the precision of RHOEND.  In the printed output
     the displayed term that is multiplied by SIGMA is called MAXCV, which
     stands for 'MAXimum Constraint Violation'.

     IPRINT should be set to 0, 1, 2 or 3, which controls the amount of
     printing during the calculation.  Specifically, there is no output if
     IPRINT=0 and there is output only at the end of the calculation if
     IPRINT=1.  Otherwise each new value of RHO and SIGMA is printed.  Further,
     the vector of variables and some function information are given either
     when RHO is reduced or when each new value of F(X) is computed in the
     cases IPRINT=2 or IPRINT=3 respectively.

     The function `cobyla_iterate` performs an iteration of the algorithm given
     F the function value at X the current variables and C the constraints.

     Typical usage is:

     >   ctx = cobyla_create(n, m, rhobeg, rhoend, iprint, maxfun);
     >   x = ...; // initial solution
     >   while (ctx.status == COBYLA_ITERATE) {
     >     f = ...; // compute function value at X
     >     c = ...; // compute constraints at X, if any
     >     cobyla_iterate, ctx, f, x, c;
     >   }
     >   if (ctx.status != COBYLA_SUCCESS) {
     >     error, swrite(format="Something wrong occured in COBYLA: %s",
     >                   ctx.reason);

     If there are no constraints (M = 0), argument C must be [] (void) or can
     be omitted.

   REFERENCES
     The algorithm is described in:

         M.J.D. Powell, "A direct search optimization method that models the
         objective and constraint functions by linear interpolation," in
         Advances in Optimization and Numerical Analysis Mathematics and Its
         Applications, vol.  275 (eds.  Susana Gomez and Jean-Pierre Hennart),
         Kluwer Academic Publishers, pp. 51-67 (1994).

   SEE ALSO: cobyla_restart, cobyla_reason.
 */

extern cobyla_restart;
/* DOCUMENT cobyla_restart(ctx);
     Restart COBYLA algorithm using the same parameters.  The return value is
     the new status of the algorithm, see `cobyla_get_status` for details.

   SEE ALSO: cobyla_create.
 */

extern cobyla_reason;
/* DOCUMENT cobyla_reason(status);
     Get a textual explanation of the status returned by COBYLA.
   SEE ALSO: cobyla_create.
 */

extern cobyla_init;
/* DOCUMENT cobyla_init;
     Initialize COBYLA interface.  It is automatically called when COBYLA
     plugin is loaded but it can safely be called again to reinitialze
     constants.
   SEE ALSO: cobyla_create.
 */
cobyla_init;


/*---------------------------------------------------------------------------*/
/* Yorick interface to Mike Powell's NEWUOA algorithm. */

local newuoa_minimize, newuoa_maximize;
/* DOCUMENT xmin = newuoa_minimize(f, x0, rhobeg, rhoend);
         or  obj = newuoa_minimize(f, x0, rhobeg, rhoend, all=1);
         or xmax = newuoa_maximize(f, x0, rhobeg, rhoend);
         or  obj = newuoa_maximize(f, x0, rhobeg, rhoend, all=1);

     Minimize or maximize the multi-variate function F starting at the initial
     point X0.  RHOBEG and RHOEND are the initial and final values of the trust
     region radius (0 < RHOEND <= RHOBEG).

     Note that the proper scaling of the variables is important for the success
     of the algorithm.  RHOBEG should be set to the typical size of the region
     to explorate and RHOEND should be set to the typical precision.

     Keyword NPT sets the number of of interpolation conditions.  Its default
     value is equal to 2*N+1 (the recommended value).

     Keyword MAXFUN sets the maximum number of function evaluations.  Its
     default value is equal to 30*N.

     If keyword ALL is true, the result is a structured object.  For
     `newuoa_minimize`, the members of the returned object are:

        obj.fmin   = minimal function value found
        obj.xmin   = corresponding parameters
        obj.nevals = number of function evaluations
        obj.status = status of the algorithm upon return
        obj.rho    = final radius of the trust region

     For `newuoa_maximize`, the two first members are:

        obj.fmax   = minimal function value found
        obj.xmax   = corresponding parameters

     Keyword ERR sets the behavior in case of abnormal termination.  If ERR=0,
     anything but a success throws an error (this is also the default
     behavior); if ERR > 0, non-fatal errors are reported by a warning message;
     if ERR < 0, non-fatal errors are silently ignored.

     Keyword VERB set the verbosity level.


   SEE ALSO: newuoa_create, newuoa_error.
 */

func newuoa_minimize(f, x0, rhobeg, rhoend, npt=, maxfun=, all=, verb=, err=)
{
  x = double(unref(x0));
  n = numberof(x);
  ctx = newuoa_create(n, (is_void(npt) ? 2*n + 1 : npt), rhobeg, rhoend,
                      (is_void(verb) ? 0 : verb),
                      (is_void(maxfun) ? 30*n : maxfun));
  fmin = [];
  do {
    fx = f(x);
    if (is_void(fmin) || fmin > fx) {
      fmin = fx;
      xmin = x;
    }
    status = newuoa_iterate(ctx, fx, x);
  } while (status == NEWUOA_ITERATE);
  if (status != NEWUOA_SUCCESS) newuoa_error, status, err;
  return (all ? save(fmin, xmin, nevals = ctx.nevals,
                     status, rho = ctx.rho) : x);
}

func newuoa_maximize(f, x0, rhobeg, rhoend, npt=, maxfun=, all=, verb=, err=)
{
  x = double(unref(x0));
  n = numberof(x);
  ctx = newuoa_create(n, (is_void(npt) ? 2*n + 1 : npt), rhobeg, rhoend,
                      (is_void(verb) ? 0 : verb),
                      (is_void(maxfun) ? 30*n : maxfun));
  fmax = [];
  do {
    fx = f(x);
    if (is_void(fmax) || fmax < fx) {
      fmax = fx;
      xmax = x;
    }
    status = newuoa_iterate(ctx, -fx, x);
  } while (status == NEWUOA_ITERATE);
  if (status != NEWUOA_SUCCESS) newuoa_error, status, err;
  return (all ? save(fmax, xmax, nevals = ctx.nevals,
                     status, rho = ctx.rho) : x);
}

func newuoa_error(status, errmode)
/* DOCUMENT newuoa_error, status;
         or newuoa_error, status, errmode;

     Report an error in NEWUOA according to the value of STATUS.  Nothing is
     done if STATUS is NEWUOA_SUCCESS; otherwise, the optional argument ERRMODE
     determines the behavior.  If ERRMODE = 0, the routine throws an error
     (this is also the default behavior); if ERRMODE > 0, non-fatal errors are
     reported by a warning message; if ERRMODE < 0, non-fatal errors are
     silently ignored.

   SEE ALSO: newuoa_reason, error.
 */
{
  if (status != NEWUOA_SUCCESS) {
    if (errmode && (status == NEWUOA_ROUNDING_ERRORS ||
                    status == NEWUOA_TOO_MANY_EVALUATIONS ||
                    status == NEWUOA_STEP_FAILED)) {
      if (errmode > 0) {
        write, format="WARNING: Something wrong occured in NEWUOA: %s",
          newuoa_reason(status);
      }
    } else {
      error, swrite(format="Something wrong occured in NEWUOA: %s",
                    newuoa_reason(status));
    }
  }
}
errs2caller, newuoa_error;

local NEWUOA_ITERATE, NEWUOA_SUCCESS, NEWUOA_BAD_NPT, NEWUOA_ROUNDING_ERRORS;
local NEWUOA_TOO_MANY_EVALUATIONS, NEWUOA_STEP_FAILED, NEWUOA_BAD_ADDRESS;
local NEWUOA_CORRUPTED;
extern newuoa_create;
extern newuoa_iterate;
/* DOCUMENT ctx = newuoa_create(n, npt, rhobeg, rhoend, iprint, maxfun);
         or status = newuoa_iterate(ctx, f, x, c);

     The function `newuoa_create` makes a new instance for Mike Powell's NEWUOA
     algorithm for minimizing a function of many variables.  The method is
     "derivatives free" (only the function values are needed).

     N is the number of variables, NPT is the number of interpolation
     conditions. Its value must be in the interval [N+2,(N+1)(N+2)/2].  The
     recommended number of points for building the quadratic model is NPT=2*N+1.

     RHOBEG and RHOEND are the initial and final values of a trust region
     radius, so both must be positive with RHOEND <= RHOBEG.  Typically RHOBEG
     should be about one tenth of the greatest expected change to a variable,
     and RHOEND should indicate the accuracy that is required in the final
     values of the variables.

     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
     amount of printing. Specifically, there is no output if IPRINT=0 and there
     is output only at the return if IPRINT=1. Otherwise, each new value of RHO
     is printed, with the best vector of variables so far and the corresponding
     value of the objective function.  Further, each new value of F with its
     variables are output if IPRINT=3.

     MAXFUN must be set to an upper bound on the number of objective function
     calls.

     The function `newuoa_iterate` performs an iteration of the algorithm given
     F the function value at the current variables X.

     Typical usage is:

     >   ctx = newuoa_create(n, npt, rhobeg, rhoend, iprint, maxfun);
     >   x = ...; // initial solution
     >   while (ctx.status == NEWUOA_ITERATE) {
     >     f = ...; // compute function value at X
     >     newuoa_iterate, ctx, f, x;
     >   }
     >   if (ctx.status != NEWUOA_SUCCESS) {
     >     error, swrite(format="Something wrong occured in NEWUOA: %s",
     >                   ctx.reason);
     >   }

     The context object CTX returned by the function `newuoa_create` has the
     following members:

         ctx.n       number of variables
         ctx.npt     number of intepolation points
         ctx.rho     radius of the trust region
         ctx.status  current status
         ctx.nevals  number of function evaluations so far
         ctx.reason  textual description of current status


   REFERENCES
     The NEWUOA algorithm is described in:

         M.J.D. Powell, "The NEWUOA software for unconstrained minimization
         without derivatives", in Large-Scale Nonlinear Optimization, editors
         G. Di Pillo and M. Roma, Springer (2006), pages 255-297.

   SEE ALSO: newuoa_restart, newuoa_reason.
 */

extern newuoa_restart;
/* DOCUMENT newuoa_restart(ctx);
     Restart NEWUOA algorithm using the same parameters.  The return value is
     the new status of the algorithm, see `newuoa_get_status` for details.

   SEE ALSO: newuoa_create.
 */

extern newuoa_reason;
/* DOCUMENT newuoa_reason(status);
     Get a textual explanation of the status returned by NEWUOA.

   SEE ALSO: newuoa_create.
 */

extern newuoa_init;
/* DOCUMENT newuoa_init;
     Initialize NEWUOA interface.  It is automatically called when NEWUOA
     plugin is loaded but it can safely be called again to reinitialze
     constants.

   SEE ALSO: newuoa_create.
 */
newuoa_init;

/*---------------------------------------------------------------------------*/
