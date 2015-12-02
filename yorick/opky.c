/*
 * opky.c --
 *
 * Yorick interface for OptimPack library.
 *-----------------------------------------------------------------------------
 *
 * This file is part of OptimPack (https://github.com/emmt/OptimPack).
 *
 * Copyright (c) 2014, 2015 Éric Thiébaut
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
static long projections_index = -1L;
static long restarts_index = -1L;
static long single_index = -1L;
static long size_index = -1L;
static long task_index = -1L;
static long xmax_index = -1L;
static long xmin_index = -1L;

static void push_string(const char *value);

/*---------------------------------------------------------------------------*/
/* OPTIMIZER OBJECT */

typedef struct _yopt_instance yopt_instance_t;
typedef struct _yopt_operations yopt_operations_t;

struct _yopt_operations {
  const char* method;
  void (*start)(yopt_instance_t* opt);
  void (*iterate)(yopt_instance_t* opt, double fx);
  void (*get_task)(yopt_instance_t* opt);
  void (*get_iterations)(yopt_instance_t* opt);
  void (*get_evaluations)(yopt_instance_t* opt);
  void (*get_restarts)(yopt_instance_t* opt);
  void (*get_projections)(yopt_instance_t* opt);
};

struct _yopt_instance {
  yopt_operations_t* ops;
  opk_object_t* optimizer;
  opk_vspace_t* vspace;
  opk_vector_t* x;
  opk_vector_t* gx;
  opk_bound_t* xl;
  opk_bound_t* xu;
  void (*rewrap)(int iarg, opk_vector_t* vect);
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
  OPK_DROP(opt->xl);
  OPK_DROP(opt->xu);
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
    opt->ops->get_task(opt);
  } else if (index == size_index) {
    ypush_long(opt->vspace->size);
  } else if (index == iterations_index) {
    opt->ops->get_iterations(opt);
  } else if (index == evaluations_index) {
    opt->ops->get_evaluations(opt);
  } else if (index == restarts_index) {
    opt->ops->get_restarts(opt);
  } else if (index == projections_index) {
    opt->ops->get_projections(opt);
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
nlcg_start(yopt_instance_t* opt)
{
  ypush_int(opk_start_nlcg(NLCG(opt->optimizer), opt->x));
}

static void
nlcg_iterate(yopt_instance_t* opt, double fx)
{
  ypush_int(opk_iterate_nlcg(NLCG(opt->optimizer), opt->x, fx, opt->gx));
}

static void
nlcg_get_task(yopt_instance_t* opt)
{
  ypush_int(opk_get_nlcg_task(NLCG(opt->optimizer)));
}

static void
nlcg_get_iterations(yopt_instance_t* opt)
{
  ypush_long(opk_get_nlcg_iterations(NLCG(opt->optimizer)));
}

static void
nlcg_get_evaluations(yopt_instance_t* opt)
{
  ypush_long(opk_get_nlcg_evaluations(NLCG(opt->optimizer)));
}

static void
nlcg_get_restarts(yopt_instance_t* opt)
{
  ypush_long(opk_get_nlcg_restarts(NLCG(opt->optimizer)));
}

static void
nlcg_get_projections(yopt_instance_t* opt)
{
  ypush_long(0);
}

static yopt_operations_t nlcg_ops = {
  "non-linear conjugate gradient (NLCG)",
  nlcg_start,
  nlcg_iterate,
  nlcg_get_task,
  nlcg_get_iterations,
  nlcg_get_evaluations,
  nlcg_get_restarts,
  nlcg_get_projections
};

/*---------------------------------------------------------------------------*/
/* LIMITED MEMORY VARIABLE METRIC (VMLM) METHOD */

#define VMLM(obj) ((opk_vmlm_t*)(obj))

static void
vmlm_start(yopt_instance_t* opt)
{
  ypush_int(opk_start_vmlm(VMLM(opt->optimizer), opt->x));
}

static void
vmlm_iterate(yopt_instance_t* opt, double fx)
{
  ypush_int(opk_iterate_vmlm(VMLM(opt->optimizer), opt->x, fx, opt->gx));
}

static void
vmlm_get_task(yopt_instance_t* opt)
{
  ypush_int(opk_get_vmlm_task(VMLM(opt->optimizer)));
}

static void
vmlm_get_iterations(yopt_instance_t* opt)
{
  ypush_long(opk_get_vmlm_iterations(VMLM(opt->optimizer)));
}

static void
vmlm_get_evaluations(yopt_instance_t* opt)
{
  ypush_long(opk_get_vmlm_evaluations(VMLM(opt->optimizer)));
}

static void
vmlm_get_restarts(yopt_instance_t* opt)
{
  ypush_long(opk_get_vmlm_restarts(VMLM(opt->optimizer)));
}

static void
vmlm_get_projections(yopt_instance_t* opt)
{
  ypush_long(0);
}

static yopt_operations_t vmlm_ops = {
  "limited memory quasi-Newton method (VMLM)",
  vmlm_start,
  vmlm_iterate,
  vmlm_get_task,
  vmlm_get_iterations,
  vmlm_get_evaluations,
  vmlm_get_restarts,
  vmlm_get_projections
};

/*---------------------------------------------------------------------------*/
/* LIMITED MEMORY VARIABLE METRIC (LBFGS) METHOD */

#define LBFGS(obj) ((opk_lbfgs_t*)(obj))

static void
lbfgs_start(yopt_instance_t* opt)
{
  ypush_int(opk_start_lbfgs(LBFGS(opt->optimizer), opt->x));
}

static void
lbfgs_iterate(yopt_instance_t* opt, double fx)
{
  ypush_int(opk_iterate_lbfgs(LBFGS(opt->optimizer), opt->x, fx, opt->gx));
}

static void
lbfgs_get_task(yopt_instance_t* opt)
{
  ypush_int(opk_get_lbfgs_task(LBFGS(opt->optimizer)));
}

static void
lbfgs_get_iterations(yopt_instance_t* opt)
{
  ypush_long(opk_get_lbfgs_iterations(LBFGS(opt->optimizer)));
}

static void
lbfgs_get_evaluations(yopt_instance_t* opt)
{
  ypush_long(opk_get_lbfgs_evaluations(LBFGS(opt->optimizer)));
}

static void
lbfgs_get_restarts(yopt_instance_t* opt)
{
  ypush_long(opk_get_lbfgs_restarts(LBFGS(opt->optimizer)));
}

static void
lbfgs_get_projections(yopt_instance_t* opt)
{
  ypush_long(0);
}

static yopt_operations_t lbfgs_ops = {
  "limited memory quasi-Newton method (LBFGS)",
  lbfgs_start,
  lbfgs_iterate,
  lbfgs_get_task,
  lbfgs_get_iterations,
  lbfgs_get_evaluations,
  lbfgs_get_restarts,
  lbfgs_get_projections
};

/*---------------------------------------------------------------------------*/
/* LIMITED MEMORY VARIABLE METRIC METHOD WITH BOUNDS */

#define VMLMB(obj) ((opk_vmlmb_t*)(obj))

static void
vmlmb_start(yopt_instance_t* opt)
{
  ypush_int(opk_start_vmlmb(VMLMB(opt->optimizer), opt->x,
                            opt->xl, opt->xu));
}

static void
vmlmb_iterate(yopt_instance_t* opt, double fx)
{
  ypush_int(opk_iterate_vmlmb(VMLMB(opt->optimizer), opt->x, fx, opt->gx,
                              opt->xl, opt->xu));
}

static void
vmlmb_get_task(yopt_instance_t* opt)
{
  ypush_int(opk_get_vmlmb_task(VMLMB(opt->optimizer)));
}

static void
vmlmb_get_iterations(yopt_instance_t* opt)
{
  ypush_long(opk_get_vmlmb_iterations(VMLMB(opt->optimizer)));
}

static void
vmlmb_get_evaluations(yopt_instance_t* opt)
{
  ypush_long(opk_get_vmlmb_evaluations(VMLMB(opt->optimizer)));
}

static void
vmlmb_get_restarts(yopt_instance_t* opt)
{
  ypush_long(opk_get_vmlmb_restarts(VMLMB(opt->optimizer)));
}

static void
vmlmb_get_projections(yopt_instance_t* opt)
{
  ypush_long(opk_get_vmlmb_evaluations(VMLMB(opt->optimizer)));
}

static yopt_operations_t vmlmb_ops = {
  "limited memory quasi-Newton method with convex constraints (VMLMB)",
  vmlmb_start,
  vmlmb_iterate,
  vmlmb_get_task,
  vmlmb_get_iterations,
  vmlmb_get_evaluations,
  vmlmb_get_restarts,
  vmlmb_get_projections
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

static opk_vector_t*
create_wrapper(opk_vspace_t* vspace, int single)
{
  if (single) {
    float dummy = 0.0f;
    return opk_wrap_simple_float_vector(vspace, &dummy, NULL, NULL);
  } else {
    double dummy = 0.0;
    return opk_wrap_simple_double_vector(vspace, &dummy, NULL, NULL);
  }
}

static void
rewrap_float(int iarg, opk_vector_t* vect)
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
rewrap_double(int iarg, opk_vector_t* vect)
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

static void
set_bound(int iarg, yopt_instance_t* opt, opk_bound_t** bnd)
{
  int rank;

  if (*bnd != NULL) {
    opk_bound_t* tmp = *bnd;
    *bnd = NULL;
    OPK_DROP(tmp);
  }
  rank = yarg_rank(iarg);
  if (rank == 0) {
    /* Got a scalar. */
    double value = ygets_d(iarg);
    *bnd = opk_new_bound(opt->vspace, OPK_BOUND_SCALAR, &value);
  } else {
    /* Must be an array.  FIXME: for now we do a copy*/
    long ntot;
    void* src;
    opk_vector_t* vector;
    if (opt->single) {
      src = ygeta_f(iarg, &ntot, NULL);
    } else {
      src = ygeta_d(iarg, &ntot, NULL);
    }
    if (ntot != opt->vspace->size) {
      y_error("bad number of elements for the bound");
    }
    vector = opk_vcreate(opt->vspace);
    if (vector == NULL) {
      y_error("failed to create a \"vector\" for the bound");
    }
    *bnd = opk_new_bound(opt->vspace, OPK_BOUND_VECTOR, vector);
    OPK_DROP(vector);
    if (*bnd != NULL) {
      if (opt->single) {
        float* dst = opk_get_simple_float_vector_data(vector);
        memcpy(dst, src, ntot*sizeof(float));
      } else {
        double* dst = opk_get_simple_double_vector_data(vector);
        memcpy(dst, src, ntot*sizeof(double));
      }
    }
  }
  if (*bnd == NULL) {
    y_error("failed to create the bound");
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
  GET_GLOBAL(projections);
  GET_GLOBAL(single);
  GET_GLOBAL(size);
  GET_GLOBAL(task);
  GET_GLOBAL(xmax);
  GET_GLOBAL(xmin);
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
  opt->x = create_wrapper(opt->vspace, single);
  opt->gx = create_wrapper(opt->vspace, single);
  if (opt->x == NULL || opt->gx == NULL) {
    y_error("failed to create working vectors");
  }
  opt->optimizer = (opk_object_t*)opk_new_nlcg_optimizer(opt->vspace, method);
  if (opt->optimizer == NULL) {
    y_error("failed to create NLCG optimizer");
  }
}

void Y_opk_lbfgs(int argc)
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
  opt->ops = &lbfgs_ops;
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
  opt->x = create_wrapper(opt->vspace, single);
  opt->gx = create_wrapper(opt->vspace, single);
  if (opt->x == NULL || opt->gx == NULL) {
    y_error("failed to create working vectors");
  }
  opt->optimizer = (opk_object_t*)opk_new_lbfgs_optimizer(opt->vspace, m);
  if (opt->optimizer == NULL) {
    y_error("failed to create LBFGS optimizer");
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
  opt->x = create_wrapper(opt->vspace, single);
  opt->gx = create_wrapper(opt->vspace, single);
  if (opt->x == NULL || opt->gx == NULL) {
    y_error("failed to create working vectors");
  }
  opt->optimizer = (opk_object_t*)opk_new_vmlm_optimizer(opt->vspace, m);
  if (opt->optimizer == NULL) {
    y_error("failed to create VMLM optimizer");
  }
}

void Y_opk_vmlmb(int argc)
{
  yopt_instance_t* opt;
  long n = 0, m = 5;
  int iarg, position = 0, single = FALSE;
  int xmin = -1, xmax = -1;

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
      } else if (index == xmin_index) {
        xmin = iarg;
      } else if (index == xmax_index) {
        xmax = iarg;
      } else {
        y_error("unknown keyword");
      }
    }
  }
  opt = (yopt_instance_t*)ypush_obj(&yopt_type, sizeof(yopt_instance_t));
  opt->ops = &vmlmb_ops;
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
  opt->x = create_wrapper(opt->vspace, single);
  opt->gx = create_wrapper(opt->vspace, single);
  if (opt->x == NULL || opt->gx == NULL) {
    y_error("failed to create working vectors");
  }
  if (xmin >= 0) {
    set_bound(xmin + 1, opt, &opt->xl);
  }
  if (xmax >= 0) {
    set_bound(xmax + 1, opt, &opt->xu);
  }
  opt->optimizer = (opk_object_t*)opk_new_vmlmb_optimizer(opt->vspace, m);
  if (opt->optimizer == NULL) {
    y_error("failed to create VMLMB optimizer");
  }
}

void Y_opk_task(int argc)
{
  yopt_instance_t* opt;

  if (argc != 1) {
    y_error("expecting exactly 1 argument");
  }
  opt = yget_obj(0, &yopt_type);
  opt->ops->get_task(opt);
}

void Y_opk_start(int argc)
{
  yopt_instance_t* opt;

  if (argc != 2) {
    y_error("expecting exactly 2 arguments");
  }
  opt = yget_obj(1, &yopt_type);
  opt->rewrap(0, opt->x);
  opt->ops->start(opt);
}

void Y_opk_iterate(int argc)
{
  yopt_instance_t* opt;
  double fx;

  if (argc != 4) {
    y_error("expecting 4 arguments");
  }
  opt = yget_obj(3, &yopt_type);
  opt->rewrap(2, opt->x);
  fx = ygets_d(1);
  opt->rewrap(0, opt->gx);
  opt->ops->iterate(opt, fx);
}
