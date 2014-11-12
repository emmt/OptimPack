/*
 * mgh-tests.c --
 *
 * Tests of routines for unconstrained minimization test suite from the MINPACK-1
 * project.  Note: "MGH" stands for "Moré, Garbow and Hillstrom".
 *
 * References:
 * [1] J. J. Moré, B. S. Garbow and K. E. Hillstrom, "Testing unconstrained
 *     optimization software," ACM Trans. Math. Software 7, 17-41 (1981).
 * [2] J. J. Moré, B. S. Garbow and K. E. Hillstrom, "Fortran subroutines for
 *     testing unconstrained optimization software," ACM Trans. Math. Software
 *     7, 136-140 (1981).
 *
 * History:
 *  - Argonne National Laboratory. MINPACK-1 Project.  March 1980.
 *    Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More.
 *  - Conversion to Yorick.  November 2001. Éric Thiébaut.
 *  - Conversion to C.  february 2014. Éric Thiébaut.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2014 Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>.
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mgh-wrappers.h"
#include "mgh-port.h"

#define MAX(a,b)   ((a) >= (b) ? (a) : (b))
#define MIN(a,b)   ((a) <= (b) ? (a) : (b))

#define TRUE  1
#define FALSE 0

#define MALLOC(type, number)  ((type*)malloc((number)*sizeof(type)))

/* Helper macros for simple loops. */
#define LOOP_0(var,num)    for (var = 0; var < num; ++var)  /* C-like */
#define LOOP_1(var,num)    for (var = 1; var <= num; ++var) /* FORTRAN-like */


/* Table of parameters for the MINPACK1 tests. */
static struct {
  int prob; /* problem number */
  int n;    /* problem size */
} prob_list[] = {
  { 1,  3},
  { 2,  6},
  { 3,  3},
  { 4,  2},
  { 5,  3},
  { 6,  5}, /* N >= 1 */
  { 7,  6}, /* N >= 2, usually 6 or 9 */
  { 7,  9}, /* N >= 2, usually 6 or 9 */
  { 8,  6}, /* N >= 1 */
  { 9,  8}, /* N >= 1 */
  {10,  2},
  {11,  4},
  {12,  3},
  {13,  8}, /* N >= 1 */
  {14,  8}, /* N >= 1 and a multiple of 2 */
  {15, 12}, /* N >= 1 and a multiple of 4 */
  {16,  2},
  {17,  4},
  {18, 30} /* 1 <= N <= MGH_TEST18_NMAX */
};
static const int NPROBLEMS = (sizeof(prob_list)/sizeof(prob_list[0]));

/* Compare two vectors (V0 is the "reference"). */
static void
v_compare(int n, const double* v0, const double* v1, double* amax_ptr, double* rmax_ptr)
{
  double amax = 0.0, rmax = 0.0;
  int j;
  for (j = 0; j < n; ++j) {
    double a0 = v0[j];
    double a1 = v1[j];
    double q = fabs(a0);
    double adif = fabs(a1 - a0);
    double rdif = (q > 0.0 ? adif/q : adif);
    amax = MAX(amax, adif);
    rmax = MAX(rmax, rdif);
  }
  if (rmax_ptr != NULL) *rmax_ptr = rmax;
  if (amax_ptr != NULL) *amax_ptr = amax;
}

#if 0
static void
v_print(FILE* output, const char* title, int n, const double v[])
{
  int i;

  if (output != NULL) {
    if (title != NULL) {
      fprintf(output, "%s\n", title);
    }
    for (i = 0; i < n; ++i) {
      fprintf(output, "  %10.8E", v[i]);
      if (i == n - 1 || (i + 1)%4 == 0) {
        fprintf(output, "\n");
      }
    }
    fflush(output);
  }
}
#endif

int
main(int argc, const char* argv[])
{
  const char* descr;
  FILE* output = stdout;
  int prob, n, nmax, k, j;
  double amax, rmax, factor;
  double  f0,  f1;
  double *x0, *x1;
  double *g0, *g1;

  /* Allocate memory. */
  nmax = 1;
  for (k = 0; k < NPROBLEMS; ++k) {
    n = prob_list[k].n;
    nmax = MAX(n, nmax);
  }
  x0 = MALLOC(double, 4*nmax);
  if (x0 == NULL) {
    fprintf(stderr, "error: insufficient memory\n");
    exit(1);
  }
  g0 = x0 +   nmax;
  x1 = x0 + 2*nmax;
  g1 = x0 + 3*nmax;

  /* Compare results for different problemes and factors. */
  for (k = 0; k < NPROBLEMS; ++k) {
    n = prob_list[k].n;
    prob = prob_list[k].prob;
    if (mgh_umck(&descr, n, prob) != MGH_SUCCESS) {
      fprintf(stderr, "error: %s\n", descr);
      exit(1);
    }
    for (factor = 1.0, j = 0; j < 4; ++j, factor *= 10.0) {
      /* Initialize the variables. */
      mgh_initpt(n, x0, prob, factor);
      mgh_umipt(n, x1, prob, factor);
      v_compare(n, x0, x1, &amax, &rmax);
      fprintf(output, "  problem: %s (n = %d) init: factor = %.0f, rmax = %g, amax = %g\n",
              descr, n, factor, rmax, amax);

      /* Compute the objective functions. */
      f0 = mgh_objfcn(n, x0, prob);
      f1 = mgh_umobj(n, x1, prob);
      v_compare(1, &f0, &f1, &amax, &rmax);
      fprintf(output, "  problem: %s (n = %d) objf: factor = %.0f, rmax = %g, amax = %g\n",
              descr, n, factor, rmax, amax);

      /* Compute the gradients. */
      mgh_grdfcn(n, x0, g0, prob);
      mgh_umgrd(n, x1, g1, prob);
      v_compare(n, g0, g1, &amax, &rmax);
      fprintf(output, "  problem: %s (n = %d) grad: factor = %.0f, rmax = %g, amax = %g\n",
              descr, n, factor, rmax, amax);
    }
  }

  /* Release ressources. */
  free((void*)x0);
  return 0;
}

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: nil
 * c-basic-offset: 2
 * fill-column: 79
 * coding: utf-8
 * End:
 */
