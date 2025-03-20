/*
 * driver.c --
 *
 * Simple driver interface to limited memory optimization methods for OptimPack library.
 *
 *----------------------------------------------------------------------------------------
 *
 * This file is part of OptimPack (https://github.com/emmt/OptimPack).
 *
 * Copyright (c) 2002-2025, Éric Thiébaut.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this
 * software and associated documentation files (the "Software"), to deal in the Software
 * without restriction, including without limitation the rights to use, copy, modify,
 * merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies
 * or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 *----------------------------------------------------------------------------------------
 */

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>

#include "optimpack-private.h"

/*--------------------------------------------------------------------------------------*/

typedef struct operations_ operations;

struct opk_optimizer_ {
    opk_object base;  /* base type (must be the first member) */
    operations* ops;
    opk_object* optimizer;
    opk_vspace* vspace;
    opk_vector* x;
    opk_vector* g;
    opk_convexset* box;
    opk_index n;
    opk_algorithm algorithm;
    int single;
};

struct operations_ {
    opk_task   (*start)(opk_optimizer* opt);
    opk_task   (*iterate)(opk_optimizer* opt, double f);
    opk_task   (*get_task)(const opk_optimizer* opt);
    opk_status (*get_status)(const opk_optimizer* opt);
    void         (*set_status)(opk_optimizer* opt, opk_status status);
    opk_index  (*get_evaluations)(const opk_optimizer* opt);
    opk_index  (*get_iterations)(const opk_optimizer* opt);
    opk_index  (*get_restarts)(const opk_optimizer* opt);
    size_t       (*get_name)(char* buf, size_t size,
                             const opk_optimizer* opt);
    size_t       (*get_description)(char* buf, size_t size,
                                    const opk_optimizer* opt);
    double       (*get_step)(const opk_optimizer* opt);
    double       (*get_gnorm)(const opk_optimizer* opt);
};

/*--------------------------------------------------------------------------------------*/
/* NON-LINEAR CONJUGATE GRADIENT (NLCG) METHOD */

#define NLCG(obj) ((opk_nlcg*)(obj))

static opk_task
nlcg_start(opk_optimizer* opt)
{
    return opk_start_nlcg(NLCG(opt->optimizer), opt->x);
}

static opk_task
nlcg_iterate(opk_optimizer* opt, double f)
{
    return opk_iterate_nlcg(NLCG(opt->optimizer), opt->x, f, opt->g);
}

static opk_task
nlcg_get_task(const opk_optimizer* opt)
{
    return opk_get_nlcg_task(NLCG(opt->optimizer));
}

static opk_status
nlcg_get_status(const opk_optimizer* opt)
{
    return opk_get_nlcg_status(NLCG(opt->optimizer));
}

static void
nlcg_set_status(opk_optimizer* opt, opk_status status)
{
    opk__set_nlcg_status(NLCG(opt->optimizer), status);
}

static opk_index
nlcg_get_evaluations(const opk_optimizer* opt)
{
    return opk_get_nlcg_evaluations(NLCG(opt->optimizer));
}

static opk_index
nlcg_get_iterations(const opk_optimizer* opt)
{
    return opk_get_nlcg_iterations(NLCG(opt->optimizer));
}

static opk_index
nlcg_get_restarts(const opk_optimizer* opt)
{
    return opk_get_nlcg_restarts(NLCG(opt->optimizer));
}

static size_t
nlcg_get_name(char* buf, size_t size, const opk_optimizer* opt)
{
    return opk_get_nlcg_name(buf, size, NLCG(opt->optimizer));
}

static size_t
nlcg_get_description(char* buf, size_t size, const opk_optimizer* opt)
{
    return opk_get_nlcg_description(buf, size, NLCG(opt->optimizer));
}

static double
nlcg_get_step(const opk_optimizer* opt)
{
    return opk_get_nlcg_step(NLCG(opt->optimizer));
}

static double
nlcg_get_gnorm(const opk_optimizer* opt)
{
    return opk_get_nlcg_gnorm(NLCG(opt->optimizer));
}

static operations nlcg_ops = {
    nlcg_start,
    nlcg_iterate,
    nlcg_get_task,
    nlcg_get_status,
    nlcg_set_status,
    nlcg_get_evaluations,
    nlcg_get_iterations,
    nlcg_get_restarts,
    nlcg_get_name,
    nlcg_get_description,
    nlcg_get_step,
    nlcg_get_gnorm
};

/*--------------------------------------------------------------------------------------*/
/* LIMITED MEMORY VARIABLE METRIC METHOD WITH BOUNDS */

#define VMLMB(obj) ((opk_vmlmb*)(obj))

static opk_task
vmlmb_start(opk_optimizer* opt)
{
    return opk_start_vmlmb(VMLMB(opt->optimizer), opt->x);
}

static opk_task
vmlmb_iterate(opk_optimizer* opt, double f)
{
    return opk_iterate_vmlmb(VMLMB(opt->optimizer), opt->x, f, opt->g);
}

static opk_task
vmlmb_get_task(const opk_optimizer* opt)
{
    return opk_get_vmlmb_task(VMLMB(opt->optimizer));
}

static opk_status
vmlmb_get_status(const opk_optimizer* opt)
{
    return opk_get_vmlmb_status(VMLMB(opt->optimizer));
}

static void
vmlmb_set_status(opk_optimizer* opt, opk_status status)
{
    opk__set_vmlmb_status(VMLMB(opt->optimizer), status);
}

static opk_index
vmlmb_get_evaluations(const opk_optimizer* opt)
{
    return opk_get_vmlmb_evaluations(VMLMB(opt->optimizer));
}

static opk_index
vmlmb_get_iterations(const opk_optimizer* opt)
{
    return opk_get_vmlmb_iterations(VMLMB(opt->optimizer));
}

static opk_index
vmlmb_get_restarts(const opk_optimizer* opt)
{
    return opk_get_vmlmb_restarts(VMLMB(opt->optimizer));
}

static size_t
vmlmb_get_name(char* buf, size_t size, const opk_optimizer* opt)
{
    return opk_copy_string(buf, size,
                           opk_get_vmlmb_method_name(VMLMB(opt->optimizer)));
}

static size_t
vmlmb_get_description(char* buf, size_t size, const opk_optimizer* opt)
{
    return opk_get_vmlmb_description(buf, size, VMLMB(opt->optimizer));
}

static double
vmlmb_get_step(const opk_optimizer* opt)
{
    return opk_get_vmlmb_step(VMLMB(opt->optimizer));
}

static double
vmlmb_get_gnorm(const opk_optimizer* opt)
{
    return opk_get_vmlmb_gnorm(VMLMB(opt->optimizer));
}

static operations vmlmb_ops = {
    vmlmb_start,
    vmlmb_iterate,
    vmlmb_get_task,
    vmlmb_get_status,
    vmlmb_set_status,
    vmlmb_get_evaluations,
    vmlmb_get_iterations,
    vmlmb_get_restarts,
    vmlmb_get_name,
    vmlmb_get_description,
    vmlmb_get_step,
    vmlmb_get_gnorm
};

/*--------------------------------------------------------------------------------------*/
/* PUBLIC INTERFACE */

#define WRAP_FLOAT(vsp, ptr)                            \
    opk_wrap_simple_float_vector(vsp, ptr, NULL, NULL)

#define WRAP_DOUBLE(vsp, ptr)                           \
    opk_wrap_simple_double_vector(vsp, ptr, NULL, NULL)

#define REWRAP_FLOAT(vec, ptr)                                  \
    opk_rewrap_simple_float_vector(vec, ptr, NULL, NULL)

#define REWRAP_DOUBLE(vec, ptr)                                 \
    opk_rewrap_simple_double_vector(vec, ptr, NULL, NULL)



static void
finalize_optimizer(opk_object* obj)
{
    opk_optimizer *opt = (opk_optimizer*)obj;
    OPK_DROP(opt->optimizer);
    OPK_DROP(opt->vspace);
    OPK_DROP(opt->x);
    OPK_DROP(opt->g);
    OPK_DROP(opt->box);
}

#define CHECK_BOUND_TYPE(t) (OPK_BOUND_NONE <= (t) &&   \
                             (t) <= OPK_BOUND_VECTOR)

extern opk_optimizer *
opk_new_optimizer(opk_algorithm algorithm, /* optimization algorithm */
                  const void* opts, /* options */
                  opk_eltype type, /* type of variables: OPK_FLOAT or OPK_DOUBLE */
                  opk_index n, /* number of variables */
                  opk_bound_type lower_type, void* lower_addr,
                  opk_bound_type upper_type, void* upper_addr,
                  opk_lnsrch* lnsrch)
{
    opk_optimizer* opt = NULL;
    opk_bool single, bounded;

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
    if ((lower_type == OPK_BOUND_NONE) != (lower_addr == NULL) ||
        (upper_type == OPK_BOUND_NONE) != (upper_addr == NULL)) {
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
    } else if (algorithm != OPK_ALGORITHM_VMLMB) {
        /* Illegal algorithm. */
        errno = EINVAL;
        return NULL;
    }

    /* Allocate optimizer and instantiate it. */
    opt = (opk_optimizer*)opk_allocate_object(finalize_optimizer,
                                              sizeof(opk_optimizer));
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
        opk_vector* lower_vector = NULL;
        opk_vector* upper_vector = NULL;
        if (lower_type == (single ? OPK_BOUND_STATIC_FLOAT
                           : OPK_BOUND_STATIC_DOUBLE)) {
            if (single) {
                lower_vector = WRAP_FLOAT(opt->vspace, lower_addr);
            } else {
                lower_vector = WRAP_DOUBLE(opt->vspace, lower_addr);
            }
            if (lower_vector == NULL) {
                goto failure;
            }
            lower_addr = lower_vector;
            lower_type = OPK_BOUND_VECTOR;
        }
        if (upper_type == (single ? OPK_BOUND_STATIC_FLOAT
                           : OPK_BOUND_STATIC_DOUBLE)) {
            if (single) {
                upper_vector = WRAP_FLOAT(opt->vspace, upper_addr);
            } else {
                upper_vector = WRAP_DOUBLE(opt->vspace, upper_addr);
            }
            if (upper_vector == NULL) {
                OPK_DROP(lower_vector);
                goto failure;
            }
            upper_addr = upper_vector;
            upper_type = OPK_BOUND_VECTOR;
        }
        opt->box = opk_new_boxset(opt->vspace,
                                  lower_type, lower_addr,
                                  upper_type, upper_addr);
        OPK_DROP(lower_vector);
        OPK_DROP(upper_vector);
        if (opt->box == NULL) {
            goto failure;
        }
    }
    if (algorithm == OPK_ALGORITHM_NLCG) {
        opt->optimizer = (opk_object*)opk_new_nlcg_optimizer(opts, opt->vspace, lnsrch);
        opt->ops = &nlcg_ops;
    } else {
        opt->optimizer = (opk_object*)opk_new_vmlmb_optimizer(opts, opt->vspace, lnsrch, opt->box);
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
opk_destroy_optimizer(opk_optimizer *opt)
{
    OPK_DROP(opt);
}

opk_task
opk_start(opk_optimizer *opt,  void* x)
{
    opk_status status;
    if (opt == NULL) {
        return OPK_TASK_ERROR;
    }
    if (x == NULL) {
        opt->ops->set_status(opt, OPK_ILLEGAL_ADDRESS);
        return OPK_TASK_ERROR;
    }
    if (opt->single) {
        status = REWRAP_FLOAT(opt->x, x);
    } else {
        status = REWRAP_DOUBLE(opt->x, x);
    }
    if (status != OPK_SUCCESS) {
        opt->ops->set_status(opt, status);
        return OPK_TASK_ERROR;
    }
    return opt->ops->start(opt);
}

opk_task
opk_iterate(opk_optimizer *opt, void* x, double f, void* g)
{
    opk_status status;
    if (opt == NULL) {
        return OPK_TASK_ERROR;
    }
    if (x == NULL || g == NULL) {
        opt->ops->set_status(opt, OPK_ILLEGAL_ADDRESS);
        return OPK_TASK_ERROR;
    }
    if (opt->single) {
        status = REWRAP_FLOAT(opt->x, x);
        if (status == OPK_SUCCESS) {
            status = REWRAP_FLOAT(opt->g, g);
        }
    } else {
        status = REWRAP_DOUBLE(opt->x, x);
        if (status == OPK_SUCCESS) {
            status = REWRAP_DOUBLE(opt->g, g);
        }
    }
    if (status != OPK_SUCCESS) {
        opt->ops->set_status(opt, status);
        return OPK_TASK_ERROR;
    }
    return opt->ops->iterate(opt, f);
}

opk_task
opk_get_task(const opk_optimizer* opt)
{
    return (opt == NULL ? OPK_TASK_ERROR : opt->ops->get_task(opt));
}

opk_status
opk_get_status(const opk_optimizer* opt)
{
    return (opt == NULL ? OPK_ILLEGAL_ADDRESS : opt->ops->get_status(opt));
}

opk_index
opk_get_evaluations(const opk_optimizer* opt)
{
    return (opt == NULL ? 0 : opt->ops->get_evaluations(opt));
}

opk_index
opk_get_iterations(const opk_optimizer* opt)
{
    return (opt == NULL ? 0 : opt->ops->get_iterations(opt));
}

opk_index
opk_get_restarts(const opk_optimizer* opt)
{
    return (opt == NULL ? 0 : opt->ops->get_restarts(opt));
}

size_t
opk_get_name(char* buf, size_t size, const opk_optimizer* opt)
{
    if (opt == NULL) {
        return opk_copy_string(buf, size, "");
    } else {
        return opt->ops->get_name(buf, size, opt);
    }
}

size_t
opk_get_description(char* buf, size_t size, const opk_optimizer* opt)
{
    return opt->ops->get_description(buf, size, opt);
}

double
opk_get_step(const opk_optimizer* opt)
{
    return (opt == NULL ? -1 : opt->ops->get_step(opt));
}

double
opk_get_gnorm(const opk_optimizer* opt)
{
    return (opt == NULL ? 0 : opt->ops->get_gnorm(opt));
}
