/*
 * nlcg-test.c --
 *
 * Tests of non-linear conjugate gradient methods for OptimPack library.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (c) 2003-2014 Éric Thiébaut
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
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include <stdio.h>

#include "optimpack.h"
#include "mgh-port.h"

#define MALLOC(type, number)  ((type*)malloc((number)*sizeof(type)))

/* Table of parameters for the MINPACK1 tests. */
static struct {
  int prob; /* problem number */
  int size;    /* problem size */
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
const int NPROBLEMS = (sizeof(prob_list)/sizeof(prob_list[0]));

static void
fatal_error(const char* reason, ...)
{
  opk_index_t length;
  va_list ap;
  length = (reason != NULL ? strlen(reason) : 0);
  fprintf(stderr, "error: ");
  va_start(ap, reason);
  vfprintf(stderr, reason, ap);
  va_end(ap);
  if (length < 1 || reason[length - 1] != '\n') {
    fprintf(stderr, "\n");
  }
  fflush(stderr);
  exit(1);
}

int
main(int argc, const char* argv[])
{
  FILE* output = stdout;
  const char* reason;
  opk_vspace_t* vspace;
  opk_vector_t* x1;
  opk_vector_t* g1;
  opk_task_t task;
  opk_nlcg_workspace_t* ws;
  double* data;
  double* x1buf;
  double* g1buf;
  int k, n, prob, mxfun;
  unsigned int method;
  const char* descr;
  double f0, f1, g0norm, g1norm, factor;
  opk_bool_t verbose = OPK_FALSE;


  /* Get problem number and size. */
  if (argc < 2 || argc > 4) {
    fprintf(stderr, "usage: %s PROB [[FACTOR] METHOD]\n", argv[0]);
    exit(1);
  }
  k = atoi(argv[1]);
  if (k < 0 || k >= NPROBLEMS) {
    fatal_error("bad problem number (%d)", k);
  }
  if (argc == 2) {
    factor = 1.0;
    method = (OPK_NLCG_HAGER_ZHANG|OPK_NLCG_SHANNO_PHUA);
  } else if (argc == 3) {
    factor = 1.0;
    method = atoi(argv[2]);
  } else {
    char *end;
    factor = strtod(argv[2], &end);
    if (end == argv[2] || *end != '\0') {
      fatal_error("bad value for FACTOR parameter (%s)", argv[2]);
    }
    method = atoi(argv[3]);
  }
  n = prob_list[k].size;
  prob = prob_list[k].prob;
  if (mgh_umck(&descr, n, prob) != MGH_SUCCESS) {
    fatal_error("%s", descr);
  }

  /* Allocate workspaces. */
  vspace = opk_new_simple_double_vector_space(n);
  if (vspace == NULL) {
    fatal_error("failed to create vector space");
  }
  data = MALLOC(double, 2*n);
  if (data == NULL) {
    fatal_error("failed to allocate memory for %ld reals", (long)(2*n));
  }
  x1buf = data;
  g1buf = data + n;
  x1 = opk_wrap_simple_double_vector(vspace, x1buf, NULL, NULL);
  if (x1 == NULL) {
    fatal_error("failed to create vector wrapper X1");
  }
  g1 = opk_wrap_simple_double_vector(vspace, g1buf, NULL, NULL);
  if (g1 == NULL) {
    fatal_error("failed to create vector wrapper G1");
  }

  ws = opk_nlcg_new(vspace, method);
  if (ws == NULL) {
    fatal_error("failed to create optimizer workspace");
  }

  fprintf(output, "Problem: %s (n = %d)  method: non-linear conjugate gradient (%u)\n",
          descr, n, method);
  mxfun = n*500;

  /* Fill X1 with initial solution. */
  mgh_umipt(n, x1buf, prob, factor);
  f0 = mgh_umobj(n, x1buf, prob);
  mgh_umgrd(n, x1buf, g1buf, prob);
  g0norm = opk_vnorm2(g1);

  task = opk_nlcg_start(ws);
  while (OPK_TRUE) {
    if (verbose) {
      fprintf(output,
              "ITER = %3d / NEVALS = %3d / TASK = %d / START = %d%s",
              opk_nlcg_get_iterations(ws),
              opk_nlcg_get_evaluations(ws),
              opk_nlcg_get_task(ws),
              opk_nlcg_get_starting(ws),
              (opk_nlcg_get_task(ws) == OPK_TASK_COMPUTE_FG ? " / " : "\n"));
    }
    if (task == OPK_TASK_COMPUTE_FG) {
      f1 = mgh_umobj(n, x1buf, prob);
      mgh_umgrd(n, x1buf, g1buf, prob);
      g1norm = opk_vnorm2(g1);
      if (verbose) {
        fprintf(output,
                "ALPHA =%9.2e / BETA = %+9.2e / F = %+15.8E / |G| =%9.2E\n",
                opk_nlcg_get_alpha(ws), opk_nlcg_get_beta(ws), f1, g1norm);
      }
      if (f1 != f1 || g1norm != g1norm) {
        /* Exit if f(x) or g(x) is NaN (Not a Number). */
        task = OPK_TASK_ERROR;
        break;
      }
    } else {
#if 0
      if (task == OPK_TASK_NEW_X || task >= OPK_TASK_FINAL_X) {
        /* print new iterate. */
        printf("%g\n", f1);
      }
#endif
      if (task != OPK_TASK_NEW_X) {
        /* Algorithm finished. */
        break;
      }
    }
    if (opk_nlcg_get_evaluations(ws) >= mxfun) {
      /* Exit if too many iterartions. */
      task = OPK_TASK_ERROR;
      break;
    }
    task = opk_nlcg_iterate(ws, x1, f1, g1);
  }

  /* Deal with the reason of stopping the algorithm. */
  if (task == OPK_TASK_FINAL_X) {
    reason = "convergence";
  } else if (task > OPK_TASK_FINAL_X) {
    reason = "*** WARNING ***";
  } else {
    reason = "*** ERROR ***";
  }
  fprintf(output, (" - %s (%d) / factor = %g\n"
                   "   start: f0 =%16.8E, ||g0|| =%10.2E\n"
                   "   final:  f =%16.8E,  ||g|| =%10.2E, %d iterations, %d function calls\n"),
          reason, task,
          factor, f0, g0norm, f1, opk_vnorm2(g1), opk_nlcg_get_iterations(ws),
          opk_nlcg_get_evaluations(ws));

#if 0
  /* Write final X and G vectors. */
  v_print(output, " FINAL X.", n, x);
  v_print(output, " FINAL G.", n, g);
#endif

  /* Release ressources. */
  opk_vdelete(x1);
  opk_vdelete(g1);
  free(data);
  opk_nlcg_delete(ws);
  opk_delete_vector_space(vspace);

  return 0;
}

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
