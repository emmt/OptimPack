/*
 * opky.c --
 *
 * Yorick interface for OptimPack library.
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
#include <errno.h>
#include <stdio.h>

#include <pstdlib.h>
#include <yapi.h>

#include "optimpack.h"

#define TRUE  1
#define FALSE 0

/* Define some macros to get rid of some GNU extensions when not compiling
   with GCC. */
#if ! (defined(__GNUC__) && __GNUC__ > 1)
#   define __attribute__(x)
#   define __inline__
#   define __FUNCTION__        ""
#   define __PRETTY_FUNCTION__ ""
#endif

PLUG_API void y_error(const char *) __attribute__ ((noreturn));

/*---------------------------------------------------------------------------*/
/* OPTIMIZER OBJECT */

enum {
  TYPE_NONE = 0,
  TYPE_NLCG,   /* non-linear conjugate gradient */
  TYPE_LBFGS   /* limited memory BFGS */
};

typedef struct _yopt_instance yopt_instance_t;
struct _yopt_instance {
  opk_vspace_t* vspace;
  void (*rewrap)(opk_vector_t* vect, int iarg);
  void* ws;
  opk_vector_t* x;
  opk_vector_t* gx;
  int type;
  int single;
};

static void yopt_free(void *);
static void yopt_print(void *);
/*static void yopt_eval(void *, int);*/
/*static void yopt_extract(void *, char *);*/

static y_userobj_t yopt_type = {
  "line search",
  yopt_free,
  yopt_print,
  /*yopt_eval*/NULL,
  /*yopt_extract*/NULL,
  NULL
};

static void yopt_free(void *ptr)
{
  yopt_instance_t* opt = (yopt_instance_t*)ptr;
  if (opt->ws != NULL) {
    switch (opt->type) {
    case TYPE_NLCG:
      opk_nlcg_delete((opk_nlcg_workspace_t*)opt->ws);
      break;
    }
  }
  if (opt->x != NULL) {
    opk_vdelete(opt->x);
  }
  if (opt->gx != NULL) {
    opk_vdelete(opt->gx);
  }
  if (opt->vspace != NULL) {
    opk_delete_vector_space(opt->vspace);
  }
}

static void yopt_print(void *ptr)
{
}
/*static void yopt_eval(void *, int);*/
/*static void yopt_extract(void *, char *);*/

/*---------------------------------------------------------------------------*/

/* Indices of keywords. */
static long method_index = -1L;
static long single_index = -1L;
static long task_index = -1L;

static void
error_handler(const char* message)
{
  y_error(message);
}

static void
rewrap_float(opk_vector_t* vect, int iarg)
{
  long ntot;
  float* data = ygeta_f(iarg, &ntot, NULL);
  if (ntot != vect->owner->size) {
    y_error("bad number of elements");
  }
  if (opk_rewrap_simple_float_vector(vect, data, NULL, NULL) != OPK_SUCCESS) {
    y_error("failed to wrap vector");
  }
}

static void
rewrap_double(opk_vector_t* vect, int iarg)
{
  long ntot;
  double* data = ygeta_d(iarg, &ntot, NULL);
  if (ntot != vect->owner->size) {
    y_error("bad number of elements");
  }
  if (opk_rewrap_simple_double_vector(vect, data, NULL, NULL) != OPK_SUCCESS) {
    y_error("failed to wrap vector");
  }
}

static void set_global_int(const char* name, int value)
{
  ypush_int(value);
  yput_global(yget_global(name, 0), 0);
  yarg_drop(1);
}

void Y_opk_init(int argc)
{
  opk_set_error_handler(error_handler);

#define GET_GLOBAL(a) a##_index = yget_global(#a, 0)
  GET_GLOBAL(method);
  GET_GLOBAL(single);
  GET_GLOBAL(task);
#undef GET_GLOBAL

#define SET_GLOBAL_INT(a) set_global_int(#a, a)
  SET_GLOBAL_INT(OPK_TASK_ERROR);
  SET_GLOBAL_INT(OPK_TASK_COMPUTE_FG);
  SET_GLOBAL_INT(OPK_TASK_NEW_X);
  SET_GLOBAL_INT(OPK_TASK_FINAL_X);
  SET_GLOBAL_INT(OPK_TASK_WARNING);
  SET_GLOBAL_INT(OPK_NLCG_FLETCHER_REEVES);
  SET_GLOBAL_INT(OPK_NLCG_HESTENES_STIEFEL);
  SET_GLOBAL_INT(OPK_NLCG_POLAK_RIBIERE_POLYAK);
  SET_GLOBAL_INT(OPK_NLCG_FLETCHER);
  SET_GLOBAL_INT(OPK_NLCG_LIU_STOREY);
  SET_GLOBAL_INT(OPK_NLCG_DAI_YUAN);
  SET_GLOBAL_INT(OPK_NLCG_PERRY_SHANNO);
  SET_GLOBAL_INT(OPK_NLCG_HAGER_ZHANG);
  SET_GLOBAL_INT(OPK_NLCG_POWELL);
  SET_GLOBAL_INT(OPK_NLCG_SHANNO_PHUA);
  SET_GLOBAL_INT(OPK_NLCG_DEFAULT);
#undef SET_GLOBAL_INT
  ypush_nil();
}

void Y_opk_nlcg(int argc)
{
  yopt_instance_t* opt;
  long n = 0;
  int iarg, position = 0, single = FALSE;
  unsigned int method = OPK_NLCG_DEFAULT;

  for (iarg = argc - 1; iarg >= 0; --iarg) {
    long index = yarg_key(iarg);
    if (index < 0L) {
      /* Non-keyword argument. */
      switch (++position) {
      case 1:
        n = ygets_l(iarg);
        if (n <= 0L) {
          y_error("illegal number of variables");
        }
        break;
      case 2:
        method = (unsigned int)ygets_l(iarg);
        break;
      default:
        y_error("too many arguments");
      }
    } else {
      /* Keyword argument. */
      --iarg;
      if (index == single_index) {
        single = yarg_true(iarg);
      } else {
        y_error("unknown keyword");
      }
    }
  }
  if (position < 1) {
    y_error("not enough arguments");
  }
  opt = (yopt_instance_t*)ypush_obj(&yopt_type, sizeof(yopt_instance_t));
  opt->type = TYPE_NLCG;
  opt->single = single;
  if (single) {
    opt->vspace = opk_new_simple_float_vector_space(n);
    opt->rewrap = rewrap_float;
  } else {
    opt->vspace = opk_new_simple_double_vector_space(n);
    opt->rewrap = rewrap_double;
  }
  if (opt->vspace == NULL) {
    y_error("failed to create vector space");
  }
  if (single) {
    float dummy = 0.0f;
    opt->x = opk_wrap_simple_float_vector(opt->vspace, &dummy, NULL, NULL);
    opt->gx = opk_wrap_simple_float_vector(opt->vspace, &dummy, NULL, NULL);
  } else {
    double dummy = 0.0;
    opt->x = opk_wrap_simple_double_vector(opt->vspace, &dummy, NULL, NULL);
    opt->gx = opk_wrap_simple_double_vector(opt->vspace, &dummy, NULL, NULL);
  }
  if (opt->x == NULL || opt->gx == NULL) {
    y_error("failed to create working vectors");
  }
  opt->ws = opk_nlcg_new(opt->vspace, method);
  if (opt->ws == NULL) {
    y_error("failed to create NLCG optimizer");
  }
}

void Y_opk_lbfgs(int argc)
{
  y_error("not yet implemented");
}

void Y_opk_task(int argc)
{
  yopt_instance_t* opt;

  if (argc != 1) {
    y_error("expecting exactly 1 argument");
  }
  opt = yget_obj(0, &yopt_type);
  switch (opt->type) {
  case TYPE_NLCG:
    ypush_int(opk_nlcg_get_task((opk_nlcg_workspace_t*)opt->ws));
    break;
  default:
    y_error("unexpected optimizer");
  }
}

void Y_opk_start(int argc)
{
  yopt_instance_t* opt;

  if (argc != 1) {
    y_error("expecting exactly 1 argument");
  }
  opt = yget_obj(0, &yopt_type);
  switch (opt->type) {
  case TYPE_NLCG:
    ypush_int(opk_nlcg_start((opk_nlcg_workspace_t*)opt->ws));
    break;
  default:
    y_error("unexpected optimizer");
  }
}

void Y_opk_iterate(int argc)
{
  yopt_instance_t* opt;
  double fx;

  if (argc != 4) {
    y_error("expecting exactly 4 arguments");
  }
  opt = yget_obj(0, &yopt_type);
  opt->rewrap(opt->x, 1);
  fx = ygets_d(2);
  opt->rewrap(opt->gx, 3);
  switch (opt->type) {
  case TYPE_NLCG:
    ypush_int(opk_nlcg_iterate((opk_nlcg_workspace_t*)opt->ws, opt->x, fx, opt->gx));
    break;
  default:
    y_error("unexpected optimizer");
  }
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
