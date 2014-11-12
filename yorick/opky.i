/*
 * opky.i --
 *
 * OPKY is OptimPack for Yorick.
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

if (is_func(plug_in)) plug_in, "opky";

local OPK_NLCG_FLETCHER_REEVES, OPK_NLCG_HESTENES_STIEFEL;
local OPK_NLCG_POLAK_RIBIERE_POLYAK, OPK_NLCG_FLETCHER;
local OPK_NLCG_LIU_STOREY, OPK_NLCG_DAI_YUAN;
local OPK_NLCG_PERRY_SHANNO, OPK_NLCG_HAGER_ZHANG;
local OPK_NLCG_POWELL, OPK_NLCG_SHANNO_PHUA;
local OPK_NLCG_DEFAULT;
extern opk_nlcg;
/* DOCUMENT opt = opk_nlcg(n);
         or opt = opk_nlcg(n, m);

     The function opk_nlcg() creates a new instance of OPKY reverse
     communication optimizer implementing a non-linear conjugate gradient
     method.  N is the size of the problem.

     Keyword METHOD may be set with the non-linear conjugate gradient method to
     use.

     Keyword SINGLE may be set true to use single precision floating point
     variables.  The default is to use double precision floating point
     variables.

     The dot notation can be used to query some members of the optimizer
     instance:
         opt.method      - the name of the optimization method;
         opt.size        - the size of the problem;
         opt.single      - single precision?
         opt.task        - current pending task;
         opt.iterations  - number of iterations;
         opt.evaluations - number of function (and gradient) evaluations;
         opt.restarts    - number of restarts;

   SEE ALSO: opk_iterate, opk_lbfgs, opk_task, opk_start.
 */

extern opk_lbfgs;
/* DOCUMENT opt = opk_lbfgs(n, m);

     The function opk_lbfgs() creates a new instance of OPKY reverse
     communication optimizer implementing a limited memory variant of the BFGS
     variable metric method.  N is the size of the problem and M is the number
     of previous steps to memorize.

     Keyword SINGLE may be set true to use single precision floating point
     variables.  The default is to use double precision floating point
     variables.

 */

extern opk_iterate;
/* DOCUMENT task = opk_iterate(opt, x, fx, gx);
     Proceed with next iteration for the optimizer OPT.  X stores the current
     variables, FX is the function value at X and GX is the gradient of the
     function at X.  See opk_task() for the interpretaion of the returned
     value.

   SEE ALSO: opk_nlcg, opk_lbfgs, opk_task, opk_start.
 */

extern opk_start;
/* DOCUMENT task = opk_start(opt);
     Start or re-start the reverse communication optimizer OPT.  See
     opk_task() for the interpretaion of the returned value.

   SEE ALSO: opk_nlcg, opk_lbfgs, opk_task, opk_iterate.
*/

local OPK_TASK_ERROR;
local OPK_TASK_COMPUTE_FG;
local OPK_TASK_NEW_X;
local OPK_TASK_FINAL_X;
local OPK_TASK_WARNING;
extern opk_task;
/* DOCUMENT task = opk_task(opt);
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

   SEE ALSO: opk_nlcg, opk_lbfgs, opk_iterate, opk_start.
*/

func opk_minimize(fg, x0, m, lbfgs=, nlcg=, single=, verb=)
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

   SEE ALSO: opk_nlcg, opk_lbfgs.
 */
{
  TRUE = 1n;
  FALSE = 0n;
  if (identof(x0) > Y_DOUBLE) {
    error, "bad data type for X0";
  }
  if (lbfgs && nlcg) {
    error, "only one of the keywords LBFG or NLCG can be set";
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
  if (lbfgs) {
    if (is_void(m)) m = 5;
    opt = opk_lbfgs(n, m, single=single);
  } else {
    if (is_void(m)) m = OPK_NLCG_DEFAULT;
    opt = opk_nlcg(n, m, single=single);
  }
  task = opk_start(opt);
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

/*
 * Local Variables:
 * mode: Yorick
 * tab-width: 8
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * fill-column: 79
 * coding: utf-8
 * ispell-local-dictionary: "american"
 * End:
 */
