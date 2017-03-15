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
/* PRIVATE DATA AND DEFINITIONS */

/* Indices of keywords. */
static long atol_index = -1L;
static long delta_index = -1L;
static long description_index = -1L;
static long dims_index = -1L;
static long epsilon_index = -1L;
static long evaluations_index = -1L;
static long flags_index = -1L;
static long fmin_index = -1L;
static long gatol_index = -1L;
static long gnorm_index = -1L;
static long grtol_index = -1L;
static long iterations_index = -1L;
static long lower_index = -1L;
static long maxiter_index = -1L;
static long mem_index = -1L;
static long name_index = -1L;
static long reason_index = -1L;
static long restarts_index = -1L;
static long rtol_index = -1L;
static long single_index = -1L;
static long size_index = -1L;
static long status_index = -1L;
static long step_index = -1L;
static long stpmin_index = -1L;
static long stpmax_index = -1L;
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
/* YORICK VALUE */

/* Structure to store a Yorick value. */
typedef struct _value  value_t;

#if 0
/* Get a reference for a Yorick value from the stack. */
static void hold_value(int iarg, value_t* v);
#endif

/* Push a Yorick value to the stack. */
static void push_value(value_t* v);

/* Finalize a Yorick value. */
static void drop_value(value_t* v);

struct _value {
  union {
    int i;
    long l;
    double d;
    void* p;
  } data;
# define VALUE_NONE    0
# define VALUE_INT     1
# define VALUE_LONG    2
# define VALUE_DOUBLE  3
# define VALUE_HANDLE   4
  int type;
};

#if 0
static void
hold_value(int iarg, value_t* v)
{
  int type = yarg_typeid(iarg);
  if (type == Y_VOID) {
    memset(v, 0, sizeof(*v));
    return;
  } else if (type == Y_INT) {
    if (yarg_rank(iarg) == 0) {
      v->data.i = ygets_i(iarg);
      v->type = VALUE_INT;
      return;
    }
  } else if (type == Y_LONG) {
    if (yarg_rank(iarg) == 0) {
      v->data.l = ygets_l(iarg);
      v->type = VALUE_LONG;
      return;
    }
  } else if (type == Y_DOUBLE) {
    if (yarg_rank(iarg) == 0) {
      v->data.d = ygets_d(iarg);
      v->type = VALUE_DOUBLE;
      return;
    }
  }
  v->data.p = yget_use(iarg);
  v->type = VALUE_HANDLE;
}
#endif

static void
push_value(value_t* v)
{
  int type = (v == NULL ? VALUE_NONE : v->type);
  if (type == VALUE_INT) {
    ypush_int(v->data.i);
  } else if (type == VALUE_LONG) {
    ypush_long(v->data.l);
   } else if (type == VALUE_DOUBLE) {
    ypush_double(v->data.d);
  } else if (type == VALUE_HANDLE && v->data.p != NULL) {
    ykeep_use(v->data.p);
  } else {
    ypush_nil();
  }
}

static void
drop_value(value_t* v)
{
  if (v != NULL && v->type == VALUE_HANDLE && v->data.p != NULL) {
    v->type = VALUE_NONE;
    ydrop_use(v->data.p);
  }
}

/*---------------------------------------------------------------------------*/
/* OPTIMIZER OBJECT */

typedef struct _yopt_instance yopt_instance_t;

struct _yopt_instance {
  opk_optimizer_t* optimizer;
  value_t lower;
  value_t upper;
  long dims[Y_DIMSIZE];
  long size;
  opk_algorithm_t algorithm;
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
  NULL,
  yopt_extract,
  NULL
};

static void
yopt_free(void* ptr)
{
  yopt_instance_t* opt = (yopt_instance_t*)ptr;
  opk_destroy_optimizer(opt->optimizer);
  drop_value(&opt->lower);
  drop_value(&opt->upper);
}

static void
yopt_print(void* ptr)
{
  yopt_instance_t* opt = (yopt_instance_t*)ptr;
  char buffer[100];
  opk_get_name(buffer, sizeof(buffer), opt->optimizer);
  y_print(buffer, FALSE);
  y_print(" implementing ", FALSE);
  opk_get_description(buffer, sizeof(buffer), opt->optimizer);
  y_print(buffer, FALSE);
  sprintf(buffer, " (size=%ld, type=%s)",
          (long)opt->size,
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
    ypush_int(opk_get_task(opt->optimizer));
  } else if (index == status_index) {
    ypush_int(opk_get_status(opt->optimizer));
  } else if (index == reason_index) {
    push_string(opk_get_reason(opk_get_status(opt->optimizer)));
  } else if (index == iterations_index) {
    ypush_long(opk_get_iterations(opt->optimizer));
  } else if (index == evaluations_index) {
    ypush_long(opk_get_evaluations(opt->optimizer));
  } else if (index == restarts_index) {
    ypush_long(opk_get_restarts(opt->optimizer));
  } else if (index == gnorm_index) {
    ypush_double(opk_get_gnorm(opt->optimizer));
  } else if (index == step_index) {
    ypush_double(opk_get_step(opt->optimizer));
  } else if (index == description_index) {
    opk_get_description(buffer, sizeof(buffer), opt->optimizer);
    push_string(buffer);
  } else if (index == name_index) {
    opk_get_name(buffer, sizeof(buffer), opt->optimizer);
    push_string(buffer);
  } else if (index == lower_index) {
    push_value(&opt->lower);
  } else if (index == upper_index) {
    push_value(&opt->upper);
  } else if (index == size_index) {
    ypush_long(opt->size);
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

static int
same_dims(const long adims[], const long bdims[])
{
  long i, ndims = adims[0];
  for (i = 0; i <= ndims; ++i) {
    if (bdims[i] != adims[i]) {
      return OPK_FALSE;
    }
  }
  return OPK_TRUE;
}

static void
copy_dims(long dst[], const long src[])
{
  long i, ndims = src[0];
  for (i = 0; i <= ndims; ++i) {
    dst[i] = src[i];
  }
}

static long
get_dims(int iarg, long dims[])
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

static void*
get_array(int iarg, const yopt_instance_t* opt)
{
  int type;
  long ntot, dims[Y_DIMSIZE];
  void* arr = ygeta_any(iarg, &ntot, dims, &type);
  if (type != (opt->single ? Y_FLOAT : Y_DOUBLE)) {
    y_error((opt->single ? "expecting an array of floats"
             : "expecting an array of doubles"));
  }
  if (! same_dims(opt->dims, dims)) {
    y_error("bad dimensions");
  }
  return arr;
}

static void
get_bound(int iarg, yopt_instance_t* opt, value_t* val,
          opk_bound_type_t* type_ptr, void** addr_ptr)
{
  int type = yarg_typeid(iarg);
  if (type <= Y_DOUBLE) {
    if (yarg_rank(iarg) == 0) {
      /* Got a scalar. */
      val->data.d = ygets_d(iarg);
      val->type = VALUE_DOUBLE;
      *type_ptr = OPK_BOUND_SCALAR_DOUBLE;
      *addr_ptr = &val->data.d;
    } else {
      /* Got an array. */
      long ntot, dims[Y_DIMSIZE];
      void* arr = ygeta_any(iarg, &ntot, dims, &type);
      if (! same_dims(opt->dims, dims)) {
        y_error("bad bound dimensions");
      }
      if (ntot != opt->size) {
        y_error("bad number of elements for the bound");
      }
      if (type != (opt->single ? Y_FLOAT : Y_DOUBLE)) {
        arr = ygeta_coerce(iarg, arr, ntot, dims, type,
                           (opt->single ? Y_FLOAT : Y_DOUBLE));
      }
      val->data.p = yget_use(iarg);
      val->type = VALUE_HANDLE;
      *type_ptr = (opt->single ? OPK_BOUND_STATIC_FLOAT
                   : OPK_BOUND_STATIC_DOUBLE);
      *addr_ptr = arr;
    }
  } else if (type == Y_VOID) {
    val->type = VALUE_NONE;
    *type_ptr = OPK_BOUND_NONE;
    *addr_ptr = NULL;
  } else {
    y_error("bad bound data type");
  }
}

static void
assign_int(long index, int value)
{
  ypush_int(value);
  yput_global(index, 0);
  yarg_drop(1);
}

static void
assign_long(long index, long value)
{
  ypush_long(value);
  yput_global(index, 0);
  yarg_drop(1);
}

static void
assign_double(long index, double value)
{
  ypush_double(value);
  yput_global(index, 0);
  yarg_drop(1);
}

static void
set_global_int(const char* name, int value)
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
  GET_GLOBAL(delta);
  GET_GLOBAL(description);
  GET_GLOBAL(dims);
  GET_GLOBAL(epsilon);
  GET_GLOBAL(evaluations);
  GET_GLOBAL(flags);
  GET_GLOBAL(fmin);
  GET_GLOBAL(gatol);
  GET_GLOBAL(gnorm);
  GET_GLOBAL(grtol);
  GET_GLOBAL(iterations);
  GET_GLOBAL(lower);
  GET_GLOBAL(maxiter);
  GET_GLOBAL(mem);
  GET_GLOBAL(name);
  GET_GLOBAL(reason);
  GET_GLOBAL(restarts);
  GET_GLOBAL(rtol);
  GET_GLOBAL(single);
  GET_GLOBAL(size);
  GET_GLOBAL(status);
  GET_GLOBAL(step);
  GET_GLOBAL(stpmin);
  GET_GLOBAL(stpmax);
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
#undef SET_GLOBAL_INT

  ypush_nil();
}

void Y_opk_get_constant(int argc)
{
  const char* name;
  long value;
  if (argc != 1) {
    y_error("expecting exactly one argument");
  }
  name = ygets_q(0);
  if (opk_get_integer_constant(name, &value) != OPK_SUCCESS) {
    ypush_nil();
  } else {
    ypush_int((int)value);
  }
}

void Y_opk_nlcg(int argc)
{
  opk_nlcg_options_t options;
  yopt_instance_t* opt;
  long size = -1;
  long dims[Y_DIMSIZE];
  int iarg, single = FALSE;

  /* Parse arguments. */
  opk_get_nlcg_default_options(&options);
  for (iarg = argc - 1; iarg >= 0; --iarg) {
    long index = yarg_key(iarg);
    if (index < 0L) {
      /* Non-keyword argument. */
      if (size != -1) {
        y_error("too many arguments");
      }
      size = get_dims(iarg, dims);
    } else {
      /* Keyword argument. */
      --iarg;
      if (index == delta_index) {
        options.delta = get_optional_double(iarg, options.delta);
      } else if (index == epsilon_index) {
        options.epsilon = get_optional_double(iarg, options.epsilon);
      } else if (index == gatol_index) {
        options.gatol = get_optional_double(iarg, options.gatol);
      } else if (index == grtol_index) {
        options.grtol = get_optional_double(iarg, options.grtol);
      } else if (index == stpmin_index) {
        options.stpmin = get_optional_double(iarg, options.stpmin);
      } else if (index == stpmax_index) {
        options.stpmax = get_optional_double(iarg, options.stpmax);
      } else if (index == fmin_index) {
        options.fmin = get_optional_double(iarg, options.fmin);
      } else if (index == single_index) {
        single = yarg_true(iarg);
      } else if (index == flags_index) {
        options.flags = get_optional_uint(iarg, options.flags);
      } else {
        y_error("unknown keyword");
      }
    }
  }
  if (size == -1) {
    y_error("not enough arguments");
  }
  if (opk_check_nlcg_options(&options) != OPK_SUCCESS) {
    y_error("invalid option(s)");
  }

  /* Create instance. */
  opt = (yopt_instance_t*)ypush_obj(&yopt_type, sizeof(yopt_instance_t));
  copy_dims(opt->dims, dims);
  opt->algorithm = OPK_ALGORITHM_NLCG;
  opt->size = size;
  opt->single = single;
  opt->optimizer = opk_new_optimizer(opt->algorithm, &options,
                                     (single ? OPK_FLOAT : OPK_DOUBLE), size,
                                     OPK_BOUND_NONE, NULL,
                                     OPK_BOUND_NONE, NULL, NULL);
  if (opt->optimizer == NULL) {
    y_error(opk_get_reason(opk_guess_status(errno)));
  }
}

static void
new_vmlmb(int argc, int blmvm)
{
  opk_vmlmb_options_t options;
  yopt_instance_t* opt;
  long size = -1;
  long dims[Y_DIMSIZE];
  int iarg, single = FALSE;
  int lower_iarg = -1;
  int upper_iarg = -1;
  opk_bound_type_t lower_type = OPK_BOUND_NONE;
  opk_bound_type_t upper_type = OPK_BOUND_NONE;
  void* lower_addr = NULL;
  void* upper_addr = NULL;

  /* Parse arguments. */
  opk_get_vmlmb_default_options(&options);
  options.blmvm = (blmvm != 0);
  for (iarg = argc - 1; iarg >= 0; --iarg) {
    long index = yarg_key(iarg);
    if (index < 0L) {
      /* Non-keyword argument. */
      if (size != -1) {
        y_error("too many arguments");
      }
      size = get_dims(iarg, dims);
    } else {
      /* Keyword argument. */
      --iarg;
      if (index == mem_index) {
        options.mem = get_optional_long(iarg, options.mem);
        if (options.mem <= 0) {
          y_error("invalid value for MEM keyword");
        }
      } else if (index == delta_index) {
        options.delta = get_optional_double(iarg, options.delta);
      } else if (index == epsilon_index) {
        options.epsilon = get_optional_double(iarg, options.epsilon);
      } else if (index == gatol_index) {
        options.gatol = get_optional_double(iarg, options.gatol);
      } else if (index == grtol_index) {
        options.grtol = get_optional_double(iarg, options.grtol);
      } else if (index == stpmin_index) {
        options.stpmin = get_optional_double(iarg, options.stpmin);
      } else if (index == stpmax_index) {
        options.stpmax = get_optional_double(iarg, options.stpmax);
      } else if (index == single_index) {
        single = yarg_true(iarg);
      } else if (index == lower_index) {
        lower_iarg = iarg;
      } else if (index == upper_index) {
        upper_iarg = iarg;
      } else {
        y_error("unknown keyword");
      }
    }
  }
  if (size == -1) {
    y_error("not enough arguments");
  }
  if (opk_check_vmlmb_options(&options) != OPK_SUCCESS) {
    y_error("invalid option(s)");
  }

  /* Create instance. */
  opt = (yopt_instance_t*)ypush_obj(&yopt_type, sizeof(yopt_instance_t));
  copy_dims(opt->dims, dims);
  opt->algorithm = OPK_ALGORITHM_VMLMB;
  opt->size = size;
  opt->single = single;

  /* Increment argulents indexes because of new instance pushed on top of the
   * stack. */
  if (lower_iarg >= 0) {
    ++lower_iarg;
  }
  if (upper_iarg >= 0) {
    ++upper_iarg;
  }

  /* Get the bounds. */
  if (lower_iarg >= 0) {
    get_bound(lower_iarg, opt, &opt->lower, &lower_type, &lower_addr);
  }
  if (upper_iarg >= 0) {
    get_bound(upper_iarg, opt, &opt->upper, &upper_type, &upper_addr);
  }

  /* Create optimizer. */
  opt->optimizer = opk_new_optimizer(opt->algorithm, &options,
                                     (single ? OPK_FLOAT : OPK_DOUBLE), size,
                                     lower_type, lower_addr,
                                     upper_type, upper_addr, NULL);
  if (opt->optimizer == NULL) {
    y_error(opk_get_reason(opk_guess_status(errno)));
  }
}

void Y_opk_vmlmb(int argc)
{
  new_vmlmb(argc, 0);
}

void Y_opk_blmvm(int argc)
{
  new_vmlmb(argc, 1);
}

#define get_opt_instance(iarg)  ((yopt_instance_t*)yget_obj(iarg, &yopt_type))

void Y_opk_start(int argc)
{
  yopt_instance_t* opt;
  void* x;
  opk_task_t task;

  if (argc != 2) y_error("expecting exactly 2 arguments");
  opt = get_opt_instance(1);
  x = get_array(0, opt);
  task = opk_start(opt->optimizer, x);
  ypush_int(task);
}

void Y_opk_iterate(int argc)
{
  yopt_instance_t* opt;
  double f;
  void* x;
  void* g;
  opk_task_t task;

  if (argc != 4) y_error("expecting exactly 4 arguments");
  opt = get_opt_instance(3);
  x = get_array(2, opt);
  f = ygets_d(1);
  g = get_array(0, opt);
  task = opk_iterate(opt->optimizer, x, f, g);
  ypush_int(task);
}

void Y_opk_get_status(int argc)
{
  if (argc != 1) y_error("expecting exactly one argument");
  ypush_int(opk_get_status(get_opt_instance(0)->optimizer));
}

void Y_opk_get_task(int argc)
{
  if (argc != 1) y_error("expecting exactly one argument");
  ypush_int(opk_get_task(get_opt_instance(0)->optimizer));
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
