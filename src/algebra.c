/*
 * algebra.c --
 *
 * High level vector space implementation for OptimPack library.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (c) 2014 Éric Thiébaut
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

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <assert.h>

#include "optimpack-private.h"

opk_vector_t*
opk_valloc(opk_vspace_t* vspace, size_t size)
{
  opk_vector_t* vect;

  if (size < sizeof(opk_vector_t)) {
    size = sizeof(opk_vector_t);
  }
  vect = (opk_vector_t*)malloc(size);
  if (vect != NULL) {
    vect->owner = vspace;
  }
  return vect;
}

void
opk_vfree(opk_vector_t* vect)
{
  if (vect != NULL) {
    free((void*)vect);
  }
}

opk_vspace_t*
opk_allocate_vector_space(const void* ident,
                          opk_index_t nvariables,
                          opk_index_t nbytes)
{
  opk_vspace_t* vspace;
  if (nvariables < 1) {
    errno = EINVAL;
    return NULL;
  }
  if (nbytes < sizeof(opk_vspace_t)) {
    nbytes = sizeof(opk_vspace_t);
  }
  vspace = (opk_vspace_t*)malloc(nbytes);
  if (vspace != NULL) {
    memset(vspace, 0, nbytes);
    vspace->ident = ident;
    vspace->size = nvariables;
  }
  return vspace;
}

void
opk_delete_vector_space(opk_vspace_t* vspace)
{
  if (vspace != NULL && vspace->finalize != NULL) {
    vspace->finalize(vspace);
  }
}

#define BAD_VECTOR(name)                                                \
  opk_error("vector does not belong to the correct space in `opk_"      \
            name "`")
#define BAD_VECTORS(name)                                               \
  opk_error("vectors do not belong to the same space in `opk_"          \
            name "`")

void
opk_vdelete(opk_vector_t* vect)
{
  if (vect != NULL) {
    opk_vspace_t* vspace = vect->owner;
    vspace->delete(vspace, vect);
  }
}

opk_vector_t*
opk_vcreate(opk_vspace_t* vspace)
{
  return vspace->create(vspace);
}

void
opk_vzero(opk_vector_t* vect)
{
  opk_vspace_t* vspace = vect->owner;
  vspace->fill(vspace, vect, 0.0);
}

void
opk_vfill(opk_vector_t* vect, double alpha)
{
  opk_vspace_t* vspace = vect->owner;
  vspace->fill(vspace, vect, alpha);
}

void
opk_vcopy(opk_vector_t* dst, const opk_vector_t* src)
{
  if (src != dst) {
    opk_vspace_t* vspace = dst->owner;
    if (src->owner != vspace) {
      BAD_VECTORS("vcopy");
    } else {
      vspace->copy(vspace, dst, src);
    }
  }
}

void
opk_vscale(opk_vector_t* dst, double alpha, const opk_vector_t* src)
{
  opk_vspace_t* vspace = dst->owner;
  if (src->owner != vspace) {
    BAD_VECTORS("vscale");
  } else if (alpha == 1.0) {
    if (src != dst) {
      vspace->copy(vspace, dst, src);
    }
  } else if (alpha == 0.0) {
    vspace->fill(vspace, dst, 0.0);
  } else {
    vspace->scale(vspace, dst, alpha, src);
  }
}

void
opk_vswap(opk_vector_t* x, opk_vector_t* y)
{
  if (x != y) {
    opk_vspace_t* vspace = x->owner;
    if (y->owner != vspace) {
      BAD_VECTORS("vswap");
    } else {
      vspace->swap(vspace, x, y);
    }
  }
}

double
opk_vdot(const opk_vector_t* x, const opk_vector_t* y)
{
  opk_vspace_t* vspace = x->owner;
  if (y->owner != vspace) {
    BAD_VECTORS("vdot");
    return 0.0;
  } else {
    return vspace->dot(vspace, x, y);
  }
}

double
opk_vnorm2(const opk_vector_t* x)
{
  opk_vspace_t* vspace = x->owner;
  return vspace->norm2(vspace, x);
}

double
opk_vnorm1(const opk_vector_t* x)
{
  opk_vspace_t* vspace = x->owner;
  return vspace->norm1(vspace, x);
}

double
opk_vnorminf(const opk_vector_t* x)
{
  opk_vspace_t* vspace = x->owner;
  return vspace->norminf(vspace, x);
}

void
opk_vaxpby(opk_vector_t* dst,
           double alpha, const opk_vector_t* x,
           double beta,  const opk_vector_t* y)
{
  opk_vspace_t* vspace = dst->owner;
  if (x->owner != vspace || y->owner != vspace) {
    BAD_VECTORS("vaxpby");
  } else if (alpha == 0.0) {
    opk_vscale(dst, beta, y);
  } else if (beta == 0.0) {
    opk_vscale(dst, alpha, x);
  } else {
    vspace->axpby(vspace, dst, alpha, x, beta, y);
  }
}

void
opk_vaxpbypcz(opk_vector_t* dst,
              double alpha, const opk_vector_t* x,
              double beta,  const opk_vector_t* y,
              double gamma, const opk_vector_t* z)
{
  opk_vspace_t* vspace = dst->owner;
  if (x->owner != vspace || y->owner != vspace || z->owner != vspace) {
    BAD_VECTORS("vaxpbypcz");
  } else if (alpha == 0.0) {
    opk_vaxpby(dst, beta, y, gamma, z);
  } else if (beta == 0.0) {
    opk_vaxpby(dst, alpha, x, gamma, z);
  } else if (gamma == 0.0) {
    opk_vaxpby(dst, alpha, x, beta, y);
  } else {
    vspace->axpbypcz(vspace, dst, alpha, x, beta, y, gamma, z);
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
