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

#include "optimpack-private.h"

#define TRUE    OPK_TRUE
#define FALSE   OPK_FALSE

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
/* PRIVATE DATA AND DEFINTIONS */

/* Indices of keywords. */
static long evaluations_index = -1L;
static long iterations_index = -1L;
static long method_index = -1L;
static long restarts_index = -1L;
static long single_index = -1L;
static long size_index = -1L;
static long task_index = -1L;

static void push_string(const char *value);

/*---------------------------------------------------------------------------*/
/* OPTIMIZER OBJECT */

typedef struct _yopt_instance yopt_instance_t;
typedef struct _yopt_operations yopt_operations_t;

struct _yopt_operations {
  const char* method;
  void (*start)(opk_object_t* optimizer);
  void (*iterate)(opk_object_t* optimizer, opk_vector_t* x,
                  double fx, opk_vector_t* gx);
  void (*get_task)(opk_object_t* optimizer);
  void (*get_iterations)(opk_object_t* optimizer);
  void (*get_evaluations)(opk_object_t* optimizer);
  void (*get_restarts)(opk_object_t* optimizer);
};

struct _yopt_instance {
  yopt_operations_t* ops;
  opk_object_t* optimizer;
  opk_vspace_t* vspace;
  opk_vector_t* x;
  opk_vector_t* gx;
  void (*rewrap)(opk_vector_t* vect, int iarg);
  int single;
};

static void yopt_free(void *);
static void yopt_print(void *);
/*static void yopt_eval(void *, int);*/
static void yopt_extract(void *, char *);

static y_userobj_t yopt_type = {
  "OPKY optimizer",
  yopt_free,
  yopt_print,
  /*yopt_eval*/NULL,
  yopt_extract,
  NULL
};

static void
yopt_free(void* ptr)
{
  yopt_instance_t* opt = (yopt_instance_t*)ptr;
  OPK_DROP(opt->optimizer);
  OPK_DROP(opt->vspace);
  OPK_DROP(opt->x);
  OPK_DROP(opt->gx);
}

static void
yopt_print(void* ptr)
{
  yopt_instance_t* opt = (yopt_instance_t*)ptr;
  char buffer[100];
  y_print(yopt_type.type_name, FALSE);
  y_print(" implementing ", FALSE);
  y_print(opt->ops->method, FALSE);
  sprintf(buffer, " method (size=%ld, type=%s)",
          (long)opt->vspace->size,
          (opt->single ? "float" : "double"));
  y_print(buffer, TRUE);
}

/*static void yopt_eval(void *, int);*/

static void
yopt_extract(void* ptr, char* member)
{
  yopt_instance_t* opt = (yopt_instance_t*)ptr;
  long index = yget_global(member, 0);
  if (index == method_index) {
    push_string(opt->ops->method);
  } else if (index == task_index) {
    opt->ops->get_task(opt->optimizer);
  } else if (index == size_index) {
    ypush_long(opt->vspace->size);
  } else if (index == iterations_index) {
    opt->ops->get_iterations(opt->optimizer);
  } else if (index == evaluations_index) {
    opt->ops->get_evaluations(opt->optimizer);
  } else if (index == restarts_index) {
    opt->ops->get_restarts(opt->optimizer);
  } else if (index == single_index) {
    ypush_int(opt->single);
  } else {
    ypush_nil();
  }
}


/*---------------------------------------------------------------------------*/
/* NON-LINEAR CONJUGATE GRADIENT (NLCG) METHOD */

#define NLCG(obj) ((opk_nlcg_t*)(obj))

static void
nlcg_start(opk_object_t* optimizer)
{
  ypush_int(opk_start_nlcg(NLCG(optimizer)));
}

static void
nlcg_iterate(opk_object_t* optimizer,
             opk_vector_t* x, double fx, opk_vector_t* gx)
{
  ypush_int(opk_iterate_nlcg(NLCG(optimizer), x, fx, gx));
}

static void
nlcg_get_task(opk_object_t* optimizer)
{
  ypush_int(opk_get_nlcg_task(NLCG(optimizer)));
}

static void
nlcg_get_iterations(opk_object_t* optimizer)
{
  ypush_long(opk_get_nlcg_iterations(NLCG(optimizer)));
}

static void
nlcg_get_evaluations(opk_object_t* optimizer)
{
  ypush_long(opk_get_nlcg_evaluations(NLCG(optimizer)));
}

static void
nlcg_get_restarts(opk_object_t* optimizer)
{
  ypush_long(opk_get_nlcg_restarts(NLCG(optimizer)));
}

static yopt_operations_t nlcg_ops = {
  "non-linear conjugate gradient (NLCG)",
  nlcg_start,
  nlcg_iterate,
  nlcg_get_task,
  nlcg_get_iterations,
  nlcg_get_evaluations,
  nlcg_get_restarts
};

/*---------------------------------------------------------------------------*/
/* LIMITED MEMORY VARIABLE METRIC (VMLM/LBFGS) METHOD */

#define VMLM(obj) ((opk_vmlm_t*)(obj))

static void
vmlm_start(opk_object_t* optimizer)
{
  ypush_int(opk_start_vmlm(VMLM(optimizer)));
}

static void
vmlm_iterate(opk_object_t* optimizer,
             opk_vector_t* x, double fx, opk_vector_t* gx)
{
  ypush_int(opk_iterate_vmlm(VMLM(optimizer), x, fx, gx));
}

static void
vmlm_get_task(opk_object_t* optimizer)
{
  ypush_int(opk_get_vmlm_task(VMLM(optimizer)));
}

static void
vmlm_get_iterations(opk_object_t* optimizer)
{
  ypush_long(opk_get_vmlm_iterations(VMLM(optimizer)));
}

static void
vmlm_get_evaluations(opk_object_t* optimizer)
{
  ypush_long(opk_get_vmlm_evaluations(VMLM(optimizer)));
}

static void
vmlm_get_restarts(opk_object_t* optimizer)
{
  ypush_long(opk_get_vmlm_restarts(VMLM(optimizer)));
}

static yopt_operations_t vmlm_ops = {
  "non-linear conjugate gradient (VMLM)",
  vmlm_start,
  vmlm_iterate,
  vmlm_get_task,
  vmlm_get_iterations,
  vmlm_get_evaluations,
  vmlm_get_restarts
};

/*---------------------------------------------------------------------------*/
/* PRIVATE FUNCTIONS */

static void push_string(const char *value)
{
  ypush_q((long *)NULL)[0] = (value ? p_strcpy((char *)value) : NULL);
}

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

/*---------------------------------------------------------------------------*/
/* BUILTIN FUNCTIONS */

void Y_opk_init(int argc)
{
  opk_set_error_handler(error_handler);

#define GET_GLOBAL(a) a##_index = yget_global(#a, 0)
  GET_GLOBAL(evaluations);
  GET_GLOBAL(iterations);
  GET_GLOBAL(method);
  GET_GLOBAL(restarts);
  GET_GLOBAL(single);
  GET_GLOBAL(size);
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
  opt->ops = &nlcg_ops;
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
  opt->optimizer = (opk_object_t*)opk_new_nlcg_optimizer(opt->vspace, method);
  if (opt->optimizer == NULL) {
    y_error("failed to create NLCG optimizer");
  }
}

void Y_opk_vmlm(int argc)
{
  yopt_instance_t* opt;
  long n = 0, m = 5;
  int iarg, position = 0, single = FALSE;
  /*unsigned int rule = OPK_BARZILAI_BORWEIN_2;*/

  for (iarg = argc - 1; iarg >= 0; --iarg) {
    long index = yarg_key(iarg);
    if (index < 0L) {
      /* Non-keyword argument. */
      switch (++position) {
      case 1:
        n = ygets_l(iarg);
        if (n < 1L) {
          y_error("illegal number of variables");
        }
        break;
      case 2:
        m = ygets_l(iarg);
        if (m < 1L) {
          y_error("illegal number of memorized steps");
        }
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
  opt->ops = &vmlm_ops;
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
  opt->optimizer = (opk_object_t*)opk_new_vmlm_optimizer(opt->vspace, m);
  if (opt->optimizer == NULL) {
    y_error("failed to create VMLM optimizer");
  }
}

void Y_opk_task(int argc)
{
  yopt_instance_t* opt;

  if (argc != 1) {
    y_error("expecting exactly 1 argument");
  }
  opt = yget_obj(0, &yopt_type);
  opt->ops->get_task(opt->optimizer);
}

void Y_opk_start(int argc)
{
  yopt_instance_t* opt;

  if (argc != 1) {
    y_error("expecting exactly 1 argument");
  }
  opt = yget_obj(0, &yopt_type);
  opt->ops->start(opt->optimizer);
}

void Y_opk_iterate(int argc)
{
  yopt_instance_t* opt;
  double fx;

  if (argc != 4) {
    y_error("expecting exactly 4 arguments");
  }
  opt = yget_obj(3, &yopt_type);
  opt->rewrap(opt->x, 2);
  fx = ygets_d(1);
  opt->rewrap(opt->gx, 0);
  opt->ops->iterate(opt->optimizer, opt->x, fx, opt->gx);
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
