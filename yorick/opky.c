/*
 * opky.c --
 *
 * Yorick interface for OptimPack library.
 *
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

#define IS_INTEGER(id)   (Y_CHAR <= (id) && (id) <= Y_LONG)
#define IS_REAL(id)      (Y_CHAR <= (id) && (id) <= Y_DOUBLE)

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
static long atol_index = -1L;
static long description_index = -1L;
static long dims_index = -1L;
static long evaluations_index = -1L;
static long flags_index = -1L;
static long gnorm_index = -1L;
static long iterations_index = -1L;
static long lower_index = -1L;
static long maxiter_index = -1L;
static long mem_index = -1L;
static long name_index = -1L;
static long projections_index = -1L;
static long reason_index = -1L;
static long restarts_index = -1L;
static long rtol_index = -1L;
static long single_index = -1L;
static long size_index = -1L;
static long status_index = -1L;
static long step_index = -1L;
static long task_index = -1L;
static long upper_index = -1L;

static void push_string(const char *value);
static void copy_dims(long dst[], const long src[]);
static int same_dims(const long adims[], const long bdims[]);
static long get_dims(int iarg, long dims[]);

static void set_global_int(const char* name, int value);

static unsigned int get_optional_uint(int iarg, unsigned int def);
static long         get_optional_long(int iarg, long def);
static double       get_optional_double(int iarg, double def);

static void assign_int(long index, int value);
static void assign_long(long index, long value);
static void assign_double(long index, double value);

/*---------------------------------------------------------------------------*/
/* OPTIMIZER OBJECT */

typedef struct _yopt_instance yopt_instance_t;
typedef struct _yopt_operations yopt_operations_t;

struct _yopt_operations {
  opk_task_t   (*start)(yopt_instance_t* opt);
  opk_task_t   (*iterate)(yopt_instance_t* opt, double fx);
  opk_task_t   (*get_task)(yopt_instance_t* opt);
  opk_status_t (*get_status)(yopt_instance_t* opt);
  unsigned int (*get_flags)(yopt_instance_t* opt);
  long         (*get_iterations)(yopt_instance_t* opt);
  long         (*get_evaluations)(yopt_instance_t* opt);
  long         (*get_restarts)(yopt_instance_t* opt);
  long         (*get_projections)(yopt_instance_t* opt);
  size_t       (*get_name)(char* buf, size_t size, yopt_instance_t* opt);
  size_t       (*get_description)(char* buf, size_t size, yopt_instance_t* opt);
  double       (*get_step)(yopt_instance_t* opt);
  double       (*get_gnorm)(yopt_instance_t* opt);
};

#define START(opt)                       ((opt)->ops->start(opt))
#define ITERATE(opt, fx)                 ((opt)->ops->iterate(opt, fx))
#define GET_TASK(opt)                    ((opt)->ops->get_task(opt))
#define GET_FLAGS(opt)                   ((opt)->ops->get_flags(opt))
#define GET_STATUS(opt)                  ((opt)->ops->get_status(opt))
#define GET_ITERATIONS(opt)              ((opt)->ops->get_iterations(opt))
#define GET_EVALUATIONS(opt)             ((opt)->ops->get_evaluations(opt))
#define GET_RESTARTS(opt)                ((opt)->ops->get_restarts(opt))
#define GET_PROJECTIONS(opt)             ((opt)->ops->get_projections(opt))
#define GET_NAME(opt, buf, siz)          ((opt)->ops->get_name(buf, siz, opt))
#define GET_DESCRIPTION(opt, buf, siz)   ((opt)->ops->get_description(buf, siz, opt))
#define GET_GNORM(opt)                   ((opt)->ops->get_gnorm(opt))
#define GET_STEP(opt)                    ((opt)->ops->get_step(opt))

struct _yopt_instance {
  long dims[Y_DIMSIZE];
  yopt_operations_t* ops;
  opk_object_t* optimizer;
  opk_vspace_t* vspace;
  opk_vector_t* x;
  opk_vector_t* gx;
  opk_bound_t* xl;
  opk_bound_t* xu;
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
  GET_NAME(opt, buffer, sizeof(buffer));
  y_print(buffer, FALSE);
  y_print(" implementing ", FALSE);
  GET_DESCRIPTION(opt, buffer, sizeof(buffer));
  y_print(buffer, FALSE);
  sprintf(buffer, " (size=%ld, type=%s)",
          (long)opt->vspace->size,
          (opt->single ? "float" : "double"));
  y_print(buffer, TRUE);
}

/*static void yopt_eval(void *, int);*/

static void
yopt_extract(void* ptr, char* member)
{
  char buffer[100];
  yopt_instance_t* opt = (yopt_instance_t*)ptr;
  long index = yget_global(member, 0);
  if (index == task_index) {
    ypush_int(GET_TASK(opt));
  } else if (index == status_index) {
    ypush_int(GET_STATUS(opt));
  } else if (index == reason_index) {
    push_string(opk_get_reason(GET_STATUS(opt)));
  } else if (index == iterations_index) {
    ypush_long(GET_ITERATIONS(opt));
  } else if (index == evaluations_index) {
    ypush_long(GET_EVALUATIONS(opt));
  } else if (index == restarts_index) {
    ypush_long(GET_RESTARTS(opt));
  } else if (index == projections_index) {
    ypush_long(GET_PROJECTIONS(opt));
  } else if (index == gnorm_index) {
    ypush_double(GET_GNORM(opt));
  } else if (index == step_index) {
    ypush_double(GET_STEP(opt));
  } else if (index == flags_index) {
    ypush_long(GET_FLAGS(opt));
  } else if (index == description_index) {
    GET_DESCRIPTION(opt, buffer, sizeof(buffer));
    push_string(buffer);
  } else if (index == name_index) {
    GET_NAME(opt, buffer, sizeof(buffer));
    push_string(buffer);
  } else if (index == size_index) {
    ypush_long(opt->vspace->size);
  } else if (index == dims_index) {
    long ndims, dims[2];
    ndims = opt->dims[0];
    dims[0] = 1;
    dims[1] = ndims + 1;
    copy_dims(ypush_l(dims), opt->dims);
  } else if (index == single_index) {
    ypush_int(opt->single);
  } else {
    ypush_nil();
  }
}

/*---------------------------------------------------------------------------*/
/* NON-LINEAR CONJUGATE GRADIENT (NLCG) METHOD */

#define NLCG(obj) ((opk_nlcg_t*)(obj))

static opk_task_t
nlcg_start(yopt_instance_t* opt)
{
  return opk_start_nlcg(NLCG(opt->optimizer), opt->x);
}

static opk_task_t
nlcg_iterate(yopt_instance_t* opt, double fx)
{
  return opk_iterate_nlcg(NLCG(opt->optimizer), opt->x, fx, opt->gx);
}

static opk_task_t
nlcg_get_task(yopt_instance_t* opt)
{
  return opk_get_nlcg_task(NLCG(opt->optimizer));
}

static opk_status_t
nlcg_get_status(yopt_instance_t* opt)
{
  return opk_get_nlcg_status(NLCG(opt->optimizer));
}

static unsigned int
nlcg_get_flags(yopt_instance_t* opt)
{
  return opk_get_nlcg_flags(NLCG(opt->optimizer));
}

static long
nlcg_get_iterations(yopt_instance_t* opt)
{
  return opk_get_nlcg_iterations(NLCG(opt->optimizer));
}

static long
nlcg_get_evaluations(yopt_instance_t* opt)
{
  return opk_get_nlcg_evaluations(NLCG(opt->optimizer));
}

static long
nlcg_get_restarts(yopt_instance_t* opt)
{
  return opk_get_nlcg_restarts(NLCG(opt->optimizer));
}

static long
nlcg_get_projections(yopt_instance_t* opt)
{
  return 0L;
}

static size_t
nlcg_get_name(char* buf, size_t size, yopt_instance_t* opt)
{
  return opk_copy_string(buf, size, "NLCG");
}

static size_t
nlcg_get_description(char* buf, size_t size, yopt_instance_t* opt)
{
  return opk_get_nlcg_description(buf, size, NLCG(opt->optimizer));
}

static double
nlcg_get_step(yopt_instance_t* opt)
{
  return opk_get_nlcg_step(NLCG(opt->optimizer));
}

static double
nlcg_get_gnorm(yopt_instance_t* opt)
{
  return opk_get_nlcg_gnorm(NLCG(opt->optimizer));
}

static yopt_operations_t nlcg_ops = {
  nlcg_start,
  nlcg_iterate,
  nlcg_get_task,
  nlcg_get_status,
  nlcg_get_flags,
  nlcg_get_iterations,
  nlcg_get_evaluations,
  nlcg_get_restarts,
  nlcg_get_projections,
  nlcg_get_name,
  nlcg_get_description,
  nlcg_get_step,
  nlcg_get_gnorm
};

/*---------------------------------------------------------------------------*/
/* LIMITED MEMORY VARIABLE METRIC METHOD WITH BOUNDS */

#define VMLMB(obj) ((opk_vmlmb_t*)(obj))

static opk_task_t
vmlmb_start(yopt_instance_t* opt)
{
  return opk_start_vmlmb(VMLMB(opt->optimizer), opt->x);
}

static opk_task_t
vmlmb_iterate(yopt_instance_t* opt, double fx)
{
  return opk_iterate_vmlmb(VMLMB(opt->optimizer), opt->x, fx, opt->gx);
}

static opk_task_t
vmlmb_get_task(yopt_instance_t* opt)
{
  return opk_get_vmlmb_task(VMLMB(opt->optimizer));
}

static opk_status_t
vmlmb_get_status(yopt_instance_t* opt)
{
  return opk_get_vmlmb_status(VMLMB(opt->optimizer));
}

static unsigned int
vmlmb_get_flags(yopt_instance_t* opt)
{
  return opk_get_vmlmb_flags(VMLMB(opt->optimizer));
}

static long
vmlmb_get_iterations(yopt_instance_t* opt)
{
  return opk_get_vmlmb_iterations(VMLMB(opt->optimizer));
}

static long
vmlmb_get_evaluations(yopt_instance_t* opt)
{
  return opk_get_vmlmb_evaluations(VMLMB(opt->optimizer));
}

static long
vmlmb_get_restarts(yopt_instance_t* opt)
{
  return opk_get_vmlmb_restarts(VMLMB(opt->optimizer));
}

static long
vmlmb_get_projections(yopt_instance_t* opt)
{
  return opk_get_vmlmb_evaluations(VMLMB(opt->optimizer));
}

static size_t
vmlmb_get_name(char* buf, size_t size, yopt_instance_t* opt)
{
  return opk_copy_string(buf, size,
                         opk_get_vmlmb_method_name(VMLMB(opt->optimizer)));
}

static size_t
vmlmb_get_description(char* buf, size_t size, yopt_instance_t* opt)
{
  return opk_get_vmlmb_description(buf, size, VMLMB(opt->optimizer));
}

static double
vmlmb_get_step(yopt_instance_t* opt)
{
  return opk_get_vmlmb_step(VMLMB(opt->optimizer));
}

static double
vmlmb_get_gnorm(yopt_instance_t* opt)
{
  return opk_get_vmlmb_gnorm(VMLMB(opt->optimizer));
}

static yopt_operations_t vmlmb_ops = {
  vmlmb_start,
  vmlmb_iterate,
  vmlmb_get_task,
  vmlmb_get_status,
  vmlmb_get_flags,
  vmlmb_get_iterations,
  vmlmb_get_evaluations,
  vmlmb_get_restarts,
  vmlmb_get_projections,
  vmlmb_get_name,
  vmlmb_get_description,
  vmlmb_get_step,
  vmlmb_get_gnorm
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

static int same_dims(const long adims[], const long bdims[])
{
  long i, ndims = adims[0];
  for (i = 0; i <= ndims; ++i) {
    if (bdims[i] != adims[i]) {
      return OPK_FALSE;
    }
  }
  return OPK_TRUE;
}

static void copy_dims(long dst[], const long src[])
{
  long i, ndims = src[0];
  for (i = 0; i <= ndims; ++i) {
    dst[i] = src[i];
  }
}

static long get_dims(int iarg, long dims[])
{
  long dim;
  int rank, type;

  type = yarg_typeid(iarg);
  if (type == Y_VOID) {
    dims[0] = 0;
    return 1L;
  }
  if (type > Y_LONG) {
    y_error("invalid type for dimension list");
  }
  rank = yarg_rank(iarg);
  if (rank == 0) {
    /* Got a scalar. */
    dim = ygets_l(iarg);
    if (dim < 1) y_error("invalid dimension length");
    dims[0] = 1;
    dims[1] = dim;
    return dim;
  }
  if (rank == 1) {
    /* Got a vector. */
    long ntot, ndims, i, number = 1L;
    long* vals = ygeta_l(iarg, &ntot, NULL);
    if (ntot > Y_DIMSIZE) {
      y_error("too many dimensions");
    }
    ndims = ntot - 1;
    if (vals[0] == ndims) {
      dims[0] = ndims;
      for (i = 1; i <= ndims; ++i) {
        dim = vals[i];
        if (dim < 1) y_error("invalid dimension length");
        dims[i] = dim;
        number *= dim;
      }
      return number;
    }
  }

  y_error("invalid dimension list");
  return -1L;
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
rewrap(int iarg, yopt_instance_t* opt, opk_vector_t* vect)
{
  long dims[Y_DIMSIZE];
  if (opt->single) {
    float* data = ygeta_f(iarg, NULL, dims);
    if (! same_dims(opt->dims, dims)) {
      y_error("bad dimensions");
    }
    if (opk_rewrap_simple_float_vector(vect, data, NULL, NULL) != OPK_SUCCESS) {
      y_error("failed to wrap vector");
    }
  } else {
    double* data = ygeta_d(iarg, NULL, dims);
    if (! same_dims(opt->dims, dims)) {
      y_error("bad dimensions");
    }
    if (opk_rewrap_simple_double_vector(vect, data, NULL, NULL) != OPK_SUCCESS) {
      y_error("failed to wrap vector");
    }
  }
}

static opk_bound_t*
get_bound(int iarg, yopt_instance_t* opt)
{
  opk_bound_t* bnd = NULL;
  int rank;

  rank = yarg_rank(iarg);
  if (rank == 0) {
    /* Got a scalar. */
    double value = ygets_d(iarg);
    bnd = opk_new_bound(opt->vspace, OPK_BOUND_SCALAR, &value);
  } else if (rank > 0) {
    /* Got an array.  FIXME: for now we do a copy. */
    long ntot, dims[Y_DIMSIZE];
    void* src;
    opk_vector_t* vector;
    if (opt->single) {
      src = ygeta_f(iarg, &ntot, dims);
    } else {
      src = ygeta_d(iarg, &ntot, dims);
    }
    if (! same_dims(opt->dims, dims)) {
      y_error("bad bound dimensions");
    }
    if (ntot != opt->vspace->size) {
      y_error("bad number of elements for the bound");
    }
    vector = opk_vcreate(opt->vspace);
    if (vector == NULL) {
      y_error("failed to create a \"vector\" for the bound");
    }
    bnd = opk_new_bound(opt->vspace, OPK_BOUND_VECTOR, vector);
    OPK_DROP(vector); /* must be done before error checking */
    if (bnd == NULL) {
      y_error("failed to create the bound");
    }
    if (opt->single) {
      float* dst = opk_get_simple_float_vector_data(vector);
      memcpy(dst, src, ntot*sizeof(float));
    } else {
      double* dst = opk_get_simple_double_vector_data(vector);
      memcpy(dst, src, ntot*sizeof(double));
    }
  } else if (! yarg_nil(iarg)) {
    y_error("invalid bound");
  }
  return bnd;
}

static void assign_int(long index, int value)
{
  ypush_int(value);
  yput_global(index, 0);
  yarg_drop(1);
}

static void assign_long(long index, long value)
{
  ypush_long(value);
  yput_global(index, 0);
  yarg_drop(1);
}

static void assign_double(long index, double value)
{
  ypush_double(value);
  yput_global(index, 0);
  yarg_drop(1);
}

static void set_global_int(const char* name, int value)
{
  assign_int(yget_global(name, 0), value);
}

static unsigned int
get_optional_uint(int iarg, unsigned int def)
{
  int type = yarg_typeid(iarg);
  if (type == Y_VOID) {
    return def;
  }
  if (type > Y_LONG || yarg_rank(iarg) != 0) {
    y_error("expecting nothing or an integer scalar");
  }
  return (unsigned int)ygets_l(iarg);
}

static long
get_optional_long(int iarg, long def)
{
  int type = yarg_typeid(iarg);
  if (type == Y_VOID) {
    return def;
  }
  if (type > Y_LONG || yarg_rank(iarg) != 0) {
    y_error("expecting nothing or an integer scalar");
  }
  return ygets_l(iarg);
}

static double
get_optional_double(int iarg, double def)
{
  int type = yarg_typeid(iarg);
  if (type == Y_VOID) {
    return def;
  }
  if (type > Y_DOUBLE || yarg_rank(iarg) != 0) {
    y_error("expecting nothing or an integer scalar");
  }
  return ygets_d(iarg);
}

/*---------------------------------------------------------------------------*/
/* BUILTIN FUNCTIONS */

void Y_opk_init(int argc)
{
  opk_set_error_handler(error_handler);

#define GET_GLOBAL(a) a##_index = yget_global(#a, 0)
  GET_GLOBAL(atol);
  GET_GLOBAL(description);
  GET_GLOBAL(dims);
  GET_GLOBAL(evaluations);
  GET_GLOBAL(flags);
  GET_GLOBAL(gnorm);
  GET_GLOBAL(iterations);
  GET_GLOBAL(lower);
  GET_GLOBAL(maxiter);
  GET_GLOBAL(mem);
  GET_GLOBAL(name);
  GET_GLOBAL(projections);
  GET_GLOBAL(reason);
  GET_GLOBAL(restarts);
  GET_GLOBAL(rtol);
  GET_GLOBAL(single);
  GET_GLOBAL(size);
  GET_GLOBAL(status);
  GET_GLOBAL(step);
  GET_GLOBAL(task);
  GET_GLOBAL(upper);
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

  SET_GLOBAL_INT(OPK_EMULATE_BLMVM);
#undef SET_GLOBAL_INT

  ypush_nil();
}

void Y_opk_nlcg(int argc)
{
  yopt_instance_t* opt;
  long n = -1;
  long dims[Y_DIMSIZE];
  int iarg, single = FALSE;
  unsigned int flags = OPK_NLCG_DEFAULT;

  for (iarg = argc - 1; iarg >= 0; --iarg) {
    long index = yarg_key(iarg);
    if (index < 0L) {
      /* Non-keyword argument. */
      if (n != -1) {
        y_error("too many arguments");
      }
      n = get_dims(iarg, dims);
    } else {
      /* Keyword argument. */
      --iarg;
      if (index == single_index) {
        single = yarg_true(iarg);
      } else if (index == flags_index) {
        flags = get_optional_uint(iarg, flags);
      } else {
        y_error("unknown keyword");
      }
    }
  }
  if (n == -1) {
    y_error("not enough arguments");
  }
  opt = (yopt_instance_t*)ypush_obj(&yopt_type, sizeof(yopt_instance_t));
  copy_dims(opt->dims, dims);
  opt->ops = &nlcg_ops;
  opt->single = single;
  if (single) {
    opt->vspace = opk_new_simple_float_vector_space(n);
  } else {
    opt->vspace = opk_new_simple_double_vector_space(n);
  }
  if (opt->vspace == NULL) {
    y_error("failed to create vector space");
  }
  opt->x = create_wrapper(opt->vspace, single);
  opt->gx = create_wrapper(opt->vspace, single);
  if (opt->x == NULL || opt->gx == NULL) {
    y_error("failed to create working vectors");
  }
  opt->optimizer = (opk_object_t*)opk_new_nlcg_optimizer(opt->vspace,
                                                         flags, NULL);
  if (opt->optimizer == NULL) {
    y_error(opk_get_reason(opk_guess_status(errno)));
  }
}

void Y_opk_vmlmb(int argc)
{
  yopt_instance_t* opt;
  long n = -1, mem = 5;
  long dims[Y_DIMSIZE];
  int iarg, single = FALSE;
  int lower = -1, upper = -1;
  unsigned int flags = 0;

  for (iarg = argc - 1; iarg >= 0; --iarg) {
    long index = yarg_key(iarg);
    if (index < 0L) {
      /* Non-keyword argument. */
      if (n != -1) {
        y_error("too many arguments");
      }
      n = get_dims(iarg, dims);
    } else {
      /* Keyword argument. */
      --iarg;
      if (index == mem_index) {
        mem = get_optional_long(iarg, mem);
        if (mem <= 0) y_error("invalid value for MEM keyword");
      } else if (index == flags_index) {
        flags = get_optional_uint(iarg, flags);
      } else if (index == single_index) {
        single = yarg_true(iarg);
      } else if (index == lower_index) {
        lower = iarg;
      } else if (index == upper_index) {
        upper = iarg;
      } else {
        y_error("unknown keyword");
      }
    }
  }
  if (n == -1) {
    y_error("not enough arguments");
  }
  opt = (yopt_instance_t*)ypush_obj(&yopt_type, sizeof(yopt_instance_t));
  copy_dims(opt->dims, dims);
  opt->ops = &vmlmb_ops;
  opt->single = single;
  if (single) {
    opt->vspace = opk_new_simple_float_vector_space(n);
  } else {
    opt->vspace = opk_new_simple_double_vector_space(n);
  }
  if (opt->vspace == NULL) {
    y_error("failed to create vector space");
  }
  opt->x = create_wrapper(opt->vspace, single);
  opt->gx = create_wrapper(opt->vspace, single);
  if (opt->x == NULL || opt->gx == NULL) {
    y_error("failed to create working vectors");
  }
  opt->xl = (lower >= 0 ? get_bound(lower + 1, opt) : NULL);
  opt->xu = (upper >= 0 ? get_bound(upper + 1, opt) : NULL);
  opt->optimizer = (opk_object_t*)opk_new_vmlmb_optimizer(opt->vspace,
                                                          mem, flags,
                                                          opt->xl, opt->xu,
                                                          NULL);
  if (opt->optimizer == NULL) {
    y_error(opk_get_reason(opk_guess_status(errno)));
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
  opk_task_t task;

  if (argc != 2) {
    y_error("expecting exactly 2 arguments");
  }
  opt = yget_obj(1, &yopt_type);
  rewrap(0, opt, opt->x);
  task = START(opt);
  ypush_int(task);
}

void Y_opk_iterate(int argc)
{
  yopt_instance_t* opt;
  double fx;
  opk_task_t task;

  if (argc != 4) {
    y_error("expecting 4 arguments");
  }
  opt = yget_obj(3, &yopt_type);
  rewrap(2, opt, opt->x);
  fx = ygets_d(1);
  rewrap(0, opt, opt->gx);
  task = ITERATE(opt, fx);
  ypush_int(task);
}

void Y_opk_get_status(int argc)
{
  yopt_instance_t* opt;
  opk_status_t status;

  if (argc != 1) y_error("expecting exactly one argument");
  opt = yget_obj(0, &yopt_type);
  status = GET_STATUS(opt);
  ypush_int(status);
}

void Y_opk_get_task(int argc)
{
  yopt_instance_t* opt;
  opk_task_t task;

  if (argc != 1) y_error("expecting exactly one argument");
  opt = yget_obj(0, &yopt_type);
  task = GET_TASK(opt);
  ypush_int(task);
}

void Y_opk_get_reason(int argc)
{
  if (argc != 1) y_error("expecting exactly one argument");
  push_string(opk_get_reason(ygets_l(0)));
}

/*---------------------------------------------------------------------------*/
/* TRUST REGION */

void Y_opk_gqtpar(int argc)
{
  double atol = 1E-8, rtol = 5E-2;
  double delta = 0.0, lambda = 0.0, f = 0.0;
  long lambda_index = -1, f_index = -1;
  long iter_index = -1, status_index = -1;
  long n, elsize, iter, maxiter = 7;
  void *A, *b, *z, *wa1, *wa2;
  long dims[Y_DIMSIZE];
  int np = 0, type, iarg, status;
  int A_iarg = -1, A_type, A_temp, A_full;
  int b_iarg = -1, b_type;

  /* Parse arguments. */
  for (iarg = argc - 1; iarg >= 0; --iarg) {
    long index = yarg_key(iarg);
    if (index < 0L) {
      /* Non-keyword argument. */
      if (++np == 1) {
        A_iarg = iarg;
      } else if (np == 2) {
        b_iarg = iarg;
      } else if (np == 3) {
        if ((delta = ygets_d(iarg)) <= 0.0) {
          y_error("DELTA must be strictly greater than zero");
        }
      } else if (np == 4) {
        lambda_index = yget_ref(iarg);
        lambda = get_optional_double(iarg, lambda);
        if (lambda  < 0.0) {
          y_error("LAMBDA must be nonnegative");
        }
      } else if (np == 5) {
        if ((f_index = yget_ref(iarg)) < 0L) {
          y_error("argument F must be set with a variable reference");
        }
      } else if (np == 6) {
        if ((iter_index = yget_ref(iarg)) < 0L) {
          y_error("argument ITER must be set with a variable reference");
        }
      } else if (np == 7) {
        if ((status_index = yget_ref(iarg)) < 0L) {
          y_error("argument STATUS must be set with a variable reference");
        }
      } else {
        y_error("too many arguments");
      }
    } else {
      /* Keyword argument. */
      --iarg;
      if (index == atol_index) {
        atol = get_optional_double(iarg, atol);
        if (atol < 0) {
          y_error("invalid value for ATOL keyword");
        }
      } else if (index == rtol_index) {
        rtol = get_optional_double(iarg, rtol);
        if (rtol < 0) {
          y_error("invalid value for RTOL keyword");
        }
      } else if (index == maxiter_index) {
        maxiter = get_optional_long(iarg, maxiter);
        if (maxiter < 1) {
          y_error("invalid value for MAXITER keyword");
        }
      } else {
        y_error("unknown keyword");
      }
    }
  }
  if (np < 3) {
    y_error("too few arguments");
  }

  /* Get the A and B arrays. */
  A_type = yarg_typeid(A_iarg);
  b_type = yarg_typeid(b_iarg);
  if (A_type == Y_DOUBLE || b_type == Y_DOUBLE) {
    type = Y_DOUBLE;
    elsize = sizeof(double);
  } else {
    type = Y_FLOAT;
    elsize = sizeof(float);
  }
  if (A_type != type) {
    if (! IS_REAL(A_type)) {
      y_error("bad data type for matrix A");
    }
    A_temp = TRUE;
  } else {
    A_temp = (yget_ref(A_iarg) < 0L);
  }
  if (b_type != type && ! IS_REAL(b_type)) {
    y_error("bad data type for vector B");
  }
  if (type == Y_DOUBLE) {
    b = ygeta_d(b_iarg, NULL, dims);
  } else {
    b = ygeta_f(b_iarg, NULL, dims);
  }
  if (dims[0] == 1L) {
    n = dims[1];
  } else if (dims[0] == 0L) {
    n = 1L;
  } else {
    y_error("bad dimensions for vector B");
    n = 0L; /* avoids compiler warnings */
  }
  if (type == Y_DOUBLE) {
    A = ygeta_d(A_iarg, NULL, dims);
  } else {
    A = ygeta_f(A_iarg, NULL, dims);
  }
  if ((dims[0] == 0L && n == 1L) ||
      (dims[0] == 2L && dims[1] == n && dims[2] == n)) {
    A_full = TRUE;
  } else {
    A_full = FALSE;
    if (dims[0] != 1L || dims[1] != n*(n + 1L)/2L) {
      y_error("bad dimensions for matrix A (or incompatible with vector B)");
    }
  }
  if (! A_temp || ! A_full) {
    /* Allocate a new array A and copy its upper triangle part.
     * Note:
     *   - lower triangle: A(i,j)  i >= j
     *   - upper triangle: A(i,j)  i <= j
     */
    long j1, j2;
    dims[0] = 2;
    dims[1] = n;
    dims[2] = n;
    if (type == Y_DOUBLE) {
      double* dst = ypush_d(dims);
      const double* src = (const double*)A;
      if (A_full) {
        for (j2 = 0; j2 < n; ++j2, dst += n, src += n) {
          for (j1 = 0; j1 <= j2; ++j1) {
            dst[j1] = src[j1];
          }
        }
      } else {
        for (j2 = 0; j2 < n; ++j2, dst += n, src += j2) {
          for (j1 = 0; j1 <= j2; ++j1) {
            dst[j1] = src[j1];
          }
        }
      }
      A = (void*)dst;
    } else {
      float* dst = ypush_f(dims);
      const float* src = (const float*)A;
      if (A_full) {
        for (j2 = 0; j2 < n; ++j2, dst += n, src += n) {
          for (j1 = 0; j1 <= j2; ++j1) {
            dst[j1] = src[j1];
          }
        }
      } else {
        for (j2 = 0; j2 < n; ++j2, dst += n, src += j2) {
          for (j1 = 0; j1 <= j2; ++j1) {
            dst[j1] = src[j1];
          }
        }
      }
      A = (void*)dst;
    }
    yarg_swap(A_iarg + 1, 0);
    yarg_drop(1);
  }

  /* Allocate workspaces on the stack. */
  z = ypush_scratch(3*elsize*n, NULL);
  wa1 = (void*)(((char*)z) +   elsize*n);
  wa2 = (void*)(((char*)z) + 2*elsize*n);

  /* Call More & Sorenson routine. */
  dims[0] = 1;
  dims[1] = n;
  if (type == Y_FLOAT) {
    float temp[2];
    temp[0] = (float)lambda;
    temp[1] = (float)f;
    status = opk_sgqt(n, A, n, b, delta, rtol, atol, maxiter, &temp[0],
                      (f_index >= 0L ? &temp[1] : NULL),
                      ypush_f(dims), &iter, z, wa1, wa2);
    lambda = temp[0];
    f = temp[1];
  } else {
    status = opk_dgqt(n, A, n, b, delta, rtol, atol, maxiter, &lambda,
                      (f_index >= 0L ? &f : NULL),
                      ypush_d(dims), &iter, z, wa1, wa2);
  }

  /* Assign output variables. */
  if (f_index >= 0L) {
    assign_double(f_index, f);
  }
  if (lambda_index >= 0L) {
    assign_double(lambda_index, lambda);
  }
  if (iter_index >= 0L) {
    assign_long(iter_index, iter);
  }
  if (status_index >= 0L) {
    assign_int(status_index, status);
  }
}

/*---------------------------------------------------------------------------*/
