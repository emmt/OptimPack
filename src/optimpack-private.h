/*
 * optimpack-private.h --
 *
 * Private macro definitions for building OptimPack library.
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

#define TRUE        1
#define FALSE       0

#define MAX(a,b)   ((a) >= (b) ? (a) : (b))
#define MIN(a,b)   ((a) <= (b) ? (a) : (b))

#define ROUND_UP(a, b)        ((((a) + ((b) - 1))/(b))*(b))
#define NEW(TYPE)             ((TYPE*)calloc(1, sizeof(TYPE)))

/* Helper macros for simple loops. */
#define LOOP_0(var,num)    for (var = 0; var < num; ++var)  /* C-like */
#define LOOP_1(var,num)    for (var = 1; var <= num; ++var) /* FORTRAN-like */

/*---------------------------------------------------------------------------*/
/* PRIVATE API FOR LINE SEARCH */

/* The base structure for line search must be exposed for line search
   methods. */
struct _opk_lnsrch_workspace {
  double stp;    /* current step length */
  double stpmin; /* lower bound for the step */
  double stpmax; /* upper bound for the step */
  double finit; /* function value at the start of the search */
  double ginit; /* directional derivative value at the start of the search */

  /* Method called by opk_lnsrch_start() to initiate a line search.  This
     method is called after having set member STP with the length of the first
     step to perform, members STPMIN and STPMAX with the lower and upper bounds
     of the step length, members FINIT and GINIT with the function value and
     the derivative the function along the search direction at the start of the
     search. */
  int (*start)(opk_lnsrch_workspace_t* ws);

  /* Method to iterate during a line search.  STP_PTR, F_PTR and D_PTR are the
     addresses of variables which store STP, F and D.  On entry, STP is the
     current step length, F is the function value for this step and D is the
     corresponding derivative of the function along the search direction.
     FIXME: On exit, if convergence is achieved (or in case of error/warning)
     STP, F and D are left unchanged; otherwise STP is the new step to try. */
  int (*iterate)(opk_lnsrch_workspace_t* ws,
                 double* stp_ptr, double f1, double d1);

  /* Method used to release any ressources speciffically allocated by the
     line search method (not the workspace itself). */
  void (*delete)(opk_lnsrch_workspace_t* ws);

  int status; /* last value returned by line search methods */
  int searching; /* true if search is in progress */
};

extern opk_lnsrch_workspace_t*
_opk_lnsrch_new(size_t size,
                int (*start)(opk_lnsrch_workspace_t* ws),
                int (*iterate)(opk_lnsrch_workspace_t* ws,
                               double* stp_ptr, double f1, double g1),
                void (*delete)(opk_lnsrch_workspace_t* ws));

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
