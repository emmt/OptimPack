/*
 * driver.c --
 *
 * Simple driver interface to limited memory optimization methods for OptimPack
 * library.
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of OptimPack (https://github.com/emmt/OptimPack).
 *
 * Copyright (C) 2002, 2015-2016 Éric Thiébaut
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
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>

#include "optimpack-private.h"

/*---------------------------------------------------------------------------*/

typedef struct _operations operations_t;

struct _opk_optimizer {
  opk_object_t base;  /* base type (must be the first member) */
  operations_t* ops;
  opk_object_t* optimizer;
  opk_vspace_t* vspace;
  opk_vector_t* x;
  opk_vector_t* g;
  opk_convexset_t* box;
  opk_index_t n;
  opk_algorithm_t algorithm;
  int single;
};

struct _operations {
  opk_task_t   (*start)(opk_optimizer_t* opt);
  opk_task_t   (*iterate)(opk_optimizer_t* opt, double f);
  opk_task_t   (*get_task)(const opk_optimizer_t* opt);
  opk_status_t (*get_status)(const opk_optimizer_t* opt);
  opk_status_t (*get_options)(void* dst, const opk_optimizer_t* src);
  opk_status_t (*set_options)(opk_optimizer_t* dst, const void* src);
  unsigned int (*get_flags)(const opk_optimizer_t* opt);
  opk_index_t  (*get_iterations)(const opk_optimizer_t* opt);
  opk_index_t  (*get_evaluations)(const opk_optimizer_t* opt);
  opk_index_t  (*get_restarts)(const opk_optimizer_t* opt);
  opk_index_t  (*get_projections)(const opk_optimizer_t* opt);
  size_t       (*get_name)(char* buf, size_t size,
                           const opk_optimizer_t* opt);
  size_t       (*get_description)(char* buf, size_t size,
                                  const opk_optimizer_t* opt);
  double       (*get_step)(const opk_optimizer_t* opt);
  double       (*get_gnorm)(const opk_optimizer_t* opt);
};

/*---------------------------------------------------------------------------*/
/* NON-LINEAR CONJUGATE GRADIENT (NLCG) METHOD */

#define NLCG(obj) ((opk_nlcg_t*)(obj))

static opk_task_t
nlcg_start(opk_optimizer_t* opt)
{
  return opk_start_nlcg(NLCG(opt->optimizer), opt->x);
}

static opk_task_t
nlcg_iterate(opk_optimizer_t* opt, double f)
{
  return opk_iterate_nlcg(NLCG(opt->optimizer), opt->x, f, opt->g);
}

static opk_task_t
nlcg_get_task(const opk_optimizer_t* opt)
{
  return opk_get_nlcg_task(NLCG(opt->optimizer));
}

static opk_status_t
nlcg_get_status(const opk_optimizer_t* opt)
{
  return opk_get_nlcg_status(NLCG(opt->optimizer));
}

static opk_status_t
nlcg_get_options(void* dst, const opk_optimizer_t* src)
{
  return opk_get_nlcg_options((opk_nlcg_options_t*)dst,
                              NLCG(src->optimizer));
}

static opk_status_t
nlcg_set_options(opk_optimizer_t* dst, const void* src)
{
  return opk_set_nlcg_options(NLCG(dst->optimizer),
                              (const opk_nlcg_options_t*)src);
}

static unsigned int
nlcg_get_flags(const opk_optimizer_t* opt)
{
  return opk_get_nlcg_flags(NLCG(opt->optimizer));
}

static opk_index_t
nlcg_get_iterations(const opk_optimizer_t* opt)
{
  return opk_get_nlcg_iterations(NLCG(opt->optimizer));
}

static opk_index_t
nlcg_get_evaluations(const opk_optimizer_t* opt)
{
  return opk_get_nlcg_evaluations(NLCG(opt->optimizer));
}

static opk_index_t
nlcg_get_restarts(const opk_optimizer_t* opt)
{
  return opk_get_nlcg_restarts(NLCG(opt->optimizer));
}

static opk_index_t
nlcg_get_projections(const opk_optimizer_t* opt)
{
  return 0L;
}

static size_t
nlcg_get_name(char* buf, size_t size, const opk_optimizer_t* opt)
{
  return opk_copy_string(buf, size, "NLCG");
}

static size_t
nlcg_get_description(char* buf, size_t size, const opk_optimizer_t* opt)
{
  return opk_get_nlcg_description(buf, size, NLCG(opt->optimizer));
}

static double
nlcg_get_step(const opk_optimizer_t* opt)
{
  return opk_get_nlcg_step(NLCG(opt->optimizer));
}

static double
nlcg_get_gnorm(const opk_optimizer_t* opt)
{
  return opk_get_nlcg_gnorm(NLCG(opt->optimizer));
}

static operations_t nlcg_ops = {
  nlcg_start,
  nlcg_iterate,
  nlcg_get_task,
  nlcg_get_status,
  nlcg_get_options,
  nlcg_set_options,
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
vmlmb_start(opk_optimizer_t* opt)
{
  return opk_start_vmlmb(VMLMB(opt->optimizer), opt->x);
}

static opk_task_t
vmlmb_iterate(opk_optimizer_t* opt, double f)
{
  return opk_iterate_vmlmb(VMLMB(opt->optimizer), opt->x, f, opt->g);
}

static opk_task_t
vmlmb_get_task(const opk_optimizer_t* opt)
{
  return opk_get_vmlmb_task(VMLMB(opt->optimizer));
}

static opk_status_t
vmlmb_get_status(const opk_optimizer_t* opt)
{
  return opk_get_vmlmb_status(VMLMB(opt->optimizer));
}

static opk_status_t
vmlmb_get_options(void* dst, const opk_optimizer_t* src)
{
  return opk_get_vmlmb_options((opk_vmlmb_options_t*)dst,
                               VMLMB(src->optimizer));
}

static opk_status_t
vmlmb_set_options(opk_optimizer_t* dst, const void* src)
{
  return opk_set_vmlmb_options(VMLMB(dst->optimizer),
                               (const opk_vmlmb_options_t*)src);
}

static unsigned int
vmlmb_get_flags(const opk_optimizer_t* opt)
{
  return opk_get_vmlmb_flags(VMLMB(opt->optimizer));
}

static opk_index_t
vmlmb_get_iterations(const opk_optimizer_t* opt)
{
  return opk_get_vmlmb_iterations(VMLMB(opt->optimizer));
}

static opk_index_t
vmlmb_get_evaluations(const opk_optimizer_t* opt)
{
  return opk_get_vmlmb_evaluations(VMLMB(opt->optimizer));
}

static opk_index_t
vmlmb_get_restarts(const opk_optimizer_t* opt)
{
  return opk_get_vmlmb_restarts(VMLMB(opt->optimizer));
}

static opk_index_t
vmlmb_get_projections(const opk_optimizer_t* opt)
{
  return opk_get_vmlmb_evaluations(VMLMB(opt->optimizer));
}

static size_t
vmlmb_get_name(char* buf, size_t size, const opk_optimizer_t* opt)
{
  return opk_copy_string(buf, size,
                         opk_get_vmlmb_method_name(VMLMB(opt->optimizer)));
}

static size_t
vmlmb_get_description(char* buf, size_t size, const opk_optimizer_t* opt)
{
  return opk_get_vmlmb_description(buf, size, VMLMB(opt->optimizer));
}

static double
vmlmb_get_step(const opk_optimizer_t* opt)
{
  return opk_get_vmlmb_step(VMLMB(opt->optimizer));
}

static double
vmlmb_get_gnorm(const opk_optimizer_t* opt)
{
  return opk_get_vmlmb_gnorm(VMLMB(opt->optimizer));
}

static operations_t vmlmb_ops = {
  vmlmb_start,
  vmlmb_iterate,
  vmlmb_get_task,
  vmlmb_get_status,
  vmlmb_get_options,
  vmlmb_set_options,
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
/* PUBLIC INTERFACE */

#define WRAP_FLOAT(vsp, ptr) \
  opk_wrap_simple_float_vector(vsp, ptr, NULL, NULL)

#define WRAP_DOUBLE(vsp, ptr) \
  opk_wrap_simple_double_vector(vsp, ptr, NULL, NULL)

#define REWRAP_FLOAT(vec, ptr) \
  opk_rewrap_simple_float_vector(vec, ptr, NULL, NULL)

#define REWRAP_DOUBLE(vec, ptr) \
  opk_rewrap_simple_double_vector(vec, ptr, NULL, NULL)



static void
finalize_optimizer(opk_object_t* obj)
{
  opk_optimizer_t *opt = (opk_optimizer_t*)obj;
  OPK_DROP(opt->optimizer);
  OPK_DROP(opt->vspace);
  OPK_DROP(opt->x);
  OPK_DROP(opt->g);
  OPK_DROP(opt->box);
}

#define CHECK_BOUND_TYPE(t) (OPK_BOUND_NONE <= (t) &&  \
                             (t) <= OPK_BOUND_VECTOR)

extern opk_optimizer_t *
opk_new_optimizer(opk_algorithm_t algorithm, /* optimization algorithm */
                  opk_type_t type, /* type of variables: OPK_FLOAT or OPK_DOUBLE */
                  opk_index_t n, /* number of variables */
                  opk_index_t m, /* number of memorized directions (m > 0, for
                                    quasi-Newton, m = 0 for non-linear
                                    conjugate gradient) */
                  unsigned int flags, /* algorithm flags */
                  opk_bound_type_t lower_type, void* lower,
                  opk_bound_type_t upper_type, void* upper,
                  opk_lnsrch_t* lnschr)
{
  opk_optimizer_t* opt = NULL;
  opk_bool_t single, bounded;

  /* Check some generic parameters. */
  if (n < 1) {
    errno = EINVAL;
    return NULL;
  }
  single = (type == OPK_FLOAT);
  if (! single && type != OPK_DOUBLE) {
    errno = EINVAL;
    return NULL;
  }

  /* Check consistency of bound settings. */
  if (! CHECK_BOUND_TYPE(lower_type) || ! CHECK_BOUND_TYPE(upper_type)) {
    errno = EINVAL;
    return NULL;
  }
  if ((lower_type == OPK_BOUND_NONE) != (lower == NULL) ||
      (upper_type == OPK_BOUND_NONE) != (upper == NULL)) {
    errno = EFAULT;
    return NULL;
  }

  /* Check algorithm and bounds. */
  bounded = (lower_type != OPK_BOUND_NONE ||
             upper_type != OPK_BOUND_NONE);
  if (algorithm == OPK_ALGORITHM_NLCG) {
    if (bounded) {
      errno = EINVAL;
      return NULL;
    }
  } else if (algorithm == OPK_ALGORITHM_VMLMB) {
    if (m < 1) {
      /* Default memory parameter. */
      m = 3;
    }
    if (m > n) {
      m = n;
    }
  } else {
    /* Illegal algorithm. */
    errno = EINVAL;
    return NULL;
  }

  /* Allocate optimizer and instanciate it. */
  opt = (opk_optimizer_t*)opk_allocate_object(finalize_optimizer,
                                              sizeof(opk_optimizer_t));
  opt->single = single;
  opt->n = n;
  opt->algorithm = algorithm;
  if (single) {
    float dummy = 0.0f;
    opt->vspace = opk_new_simple_float_vector_space(n);
    if (opt->vspace != NULL) {
      opt->x = WRAP_FLOAT(opt->vspace, &dummy);
      opt->g = WRAP_FLOAT(opt->vspace, &dummy);
    }
  } else {
    double dummy = 0.0;
    opt->vspace = opk_new_simple_double_vector_space(n);
    if (opt->vspace != NULL) {
      opt->x = WRAP_DOUBLE(opt->vspace, &dummy);
      opt->g = WRAP_DOUBLE(opt->vspace, &dummy);
    }
  }
  if (opt->vspace == NULL || opt->x == NULL || opt->g == NULL) {
    goto failure;
  }
  if (bounded) {
    /* Take care of simply wrapping static arrays if the bounds are provided as
       such. */
    opk_vector_t* lower_vector = NULL;
    opk_vector_t* upper_vector = NULL;
    if (lower_type == (single ? OPK_BOUND_STATIC_FLOAT
                       : OPK_BOUND_STATIC_DOUBLE)) {
      if (single) {
        lower_vector = WRAP_FLOAT(opt->vspace, lower);
      } else {
        lower_vector = WRAP_DOUBLE(opt->vspace, lower);
      }
      if (lower_vector == NULL) {
        goto failure;
      }
      lower = lower_vector;
      lower_type = OPK_BOUND_VECTOR;
    }
    if (upper_type == (single ? OPK_BOUND_STATIC_FLOAT
                       : OPK_BOUND_STATIC_DOUBLE)) {
      if (single) {
        upper_vector = WRAP_FLOAT(opt->vspace, upper);
      } else {
        upper_vector = WRAP_DOUBLE(opt->vspace, upper);
      }
      if (upper_vector == NULL) {
        OPK_DROP(lower_vector);
        goto failure;
      }
      upper = upper_vector;
      upper_type = OPK_BOUND_VECTOR;
    }
    opt->box = opk_new_boxset(opt->vspace,
                              lower_type, lower,
                              upper_type, upper);
    OPK_DROP(lower_vector);
    OPK_DROP(upper_vector);
    if (opt->box == NULL) {
      goto failure;
    }
  }
  if (algorithm == OPK_ALGORITHM_NLCG) {
    opt->optimizer = (opk_object_t*)opk_new_nlcg_optimizer(opt->vspace,
                                                           flags, lnschr);
    opt->ops = &nlcg_ops;
  } else {
    opt->optimizer = (opk_object_t*)opk_new_vmlmb_optimizer(opt->vspace,
                                                            m, flags,
                                                            opt->box,
                                                            lnschr);
    opt->ops = &vmlmb_ops;
  }
  if (opt->optimizer != NULL) {
    return opt;
  }

 failure:
  OPK_DROP(opt);
  return NULL;
}

void
opk_destroy_optimizer(opk_optimizer_t *opt)
{
  OPK_DROP(opt);
}

opk_task_t
opk_start(opk_optimizer_t *opt, opk_type_t type, opk_index_t n,  void* x)
{
  if (opt == NULL || x == NULL ||
      type != (opt->single ? OPK_FLOAT : OPK_DOUBLE) ||
      opt->n != n) {
    return OPK_TASK_ERROR;
  }
  if (opt->single) {
    if (REWRAP_FLOAT(opt->x, x) != OPK_SUCCESS) {
      return OPK_TASK_ERROR;
    }
  } else {
    if (REWRAP_DOUBLE(opt->x, x) != OPK_SUCCESS) {
      return OPK_TASK_ERROR;
    }
  }
  return opt->ops->start(opt);
}

opk_task_t
opk_iterate(opk_optimizer_t *opt, opk_type_t type, opk_index_t n,
            void* x, double f, void* g)
{
  if (opt == NULL || x == NULL || g == NULL ||
      type != (opt->single ? OPK_FLOAT : OPK_DOUBLE) ||
      opt->n != n) {
    return OPK_TASK_ERROR;
  }
  if (opt->single) {
    if (REWRAP_FLOAT(opt->x, x) != OPK_SUCCESS ||
        REWRAP_FLOAT(opt->g, g) != OPK_SUCCESS) {
      return OPK_TASK_ERROR;
    }
  } else {
    if (REWRAP_DOUBLE(opt->x, x) != OPK_SUCCESS ||
        REWRAP_DOUBLE(opt->g, g) != OPK_SUCCESS) {
      return OPK_TASK_ERROR;
    }
  }
  return opt->ops->iterate(opt, f);
}

opk_task_t
opk_get_task(const opk_optimizer_t* opt)
{
  return (opt == NULL ? OPK_TASK_ERROR : opt->ops->get_task(opt));
}

opk_status_t
opk_get_status(const opk_optimizer_t* opt)
{
  return (opt == NULL ? OPK_ILLEGAL_ADDRESS : opt->ops->get_status(opt));
}

opk_status_t
opk_get_options(void* dst, const opk_optimizer_t* src)
{
  return (src == NULL || dst == NULL ? OPK_ILLEGAL_ADDRESS :
          src->ops->get_options(dst, src));
}

opk_status_t
opk_set_options(opk_optimizer_t* dst, const void* src)
{
  return (src == NULL || dst == NULL ? OPK_ILLEGAL_ADDRESS :
          dst->ops->set_options(dst, src));
}

unsigned int
opk_get_flags(const opk_optimizer_t* opt)
{
  return (opt == NULL ? 0 : opt->ops->get_flags(opt));
}

opk_index_t
opk_get_iterations(const opk_optimizer_t* opt)
{
  return (opt == NULL ? 0 : opt->ops->get_iterations(opt));
}

opk_index_t
opk_get_restarts(const opk_optimizer_t* opt)
{
  return (opt == NULL ? 0 : opt->ops->get_restarts(opt));
}

opk_index_t
opk_get_evaluations(const opk_optimizer_t* opt)
{
  return (opt == NULL ? 0 : opt->ops->get_evaluations(opt));
}

opk_index_t
opk_get_projections(const opk_optimizer_t* opt)
{
  return (opt == NULL ? 0 : opt->ops->get_projections(opt));
}

size_t
opk_get_name(char* buf, size_t size, const opk_optimizer_t* opt)
{
  if (opt == NULL) {
    return opk_copy_string(buf, size, "");
  } else {
    return opt->ops->get_name(buf, size, opt);
  }
}

size_t
opk_get_description(char* buf, size_t size, const opk_optimizer_t* opt)
{
  return opt->ops->get_description(buf, size, opt);
}

double
opk_get_step(const opk_optimizer_t* opt)
{
  return (opt == NULL ? -1 : opt->ops->get_step(opt));
}

double
opk_get_gnorm(const opk_optimizer_t* opt)
{
  return (opt == NULL ? 0 : opt->ops->get_gnorm(opt));
}
