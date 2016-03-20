/*
 * utils.c --
 *
 * Error management and miscellaneaous functions in OptimPack library.
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
#include <stdio.h>
#include <errno.h>

#include "optimpack.h"

extern size_t
opk_copy_string(char* dst, size_t size, const char* src)
{
  size_t length = (src == NULL || src[0] == '\0' ? 0 : strlen(src));
  if (dst != NULL) {
    if (size > length) {
      if (length > 0) {
        memcpy(dst, src, length + 1);
      } else {
        dst[0] = '\0';
      }
    } else {
      if (size > 1) {
        memcpy(dst, src, size - 1);
      }
      if (size > 0) {
        dst[size - 1] = '\0';
      }
    }
  }
  return length + 1;
}

static void
default_error_handler(const char* message)
{
  fprintf(stderr, "ERROR: %s\n", message);
  abort();
}

static void (*error_handler)(const char* message) = default_error_handler;

opk_error_handler*
opk_get_error_handler(void)
{
  return error_handler;
}

opk_error_handler*
opk_set_error_handler(opk_error_handler* new_handler)
{
  opk_error_handler* old_handler = error_handler;
  error_handler = (new_handler == NULL ? default_error_handler : new_handler);
  return old_handler;
}

void
opk_error(const char* reason)
{
  error_handler(reason);
}

const char*
opk_get_reason(opk_status_t status)
{
#define _OPK_STATUS(a,b,c) case b: return c;
  switch (status) {
    _OPK_STATUS_LIST
  default: return "";
  }
#undef _OPK_STATUS
}

opk_status_t
opk_guess_status()
{
  switch (errno) {
#ifdef ENOMEM
  case ENOMEM: return OPK_INSUFFICIENT_MEMORY;
#endif
#ifdef EFAULT
  case EFAULT: return OPK_ILLEGAL_ADDRESS;
#endif
#ifdef EINVAL
  case EINVAL: return OPK_INVALID_ARGUMENT;
#endif
  default: return -1;
  }
}

extern opk_status_t
opk_get_integer_constant(const char* name, long* ptr)
{
  long dummy;
  if (ptr == NULL) {
    ptr = &dummy;
  }

  if (name != NULL && name[0] != 0) {
#define TEST(sym) if (strcmp(name, #sym) == 0) { *ptr = (long)(sym); \
                                                  return OPK_SUCCESS; }
#define _OPK_STATUS(val,sym,msg) TEST(sym);
    _OPK_STATUS_LIST;
#undef _OPK_STATUS

    TEST(OPK_TRUE);
    TEST(OPK_FALSE);

    TEST(OPK_FLOAT);
    TEST(OPK_DOUBLE);

    TEST(OPK_TASK_ERROR);
    TEST(OPK_TASK_START);
    TEST(OPK_TASK_COMPUTE_FG);
    TEST(OPK_TASK_NEW_X);
    TEST(OPK_TASK_FINAL_X);
    TEST(OPK_TASK_WARNING);

    TEST(OPK_LNSRCH_ERROR);
    TEST(OPK_LNSRCH_SEARCH);
    TEST(OPK_LNSRCH_CONVERGENCE);
    TEST(OPK_LNSRCH_WARNING);

    TEST(OPK_NLCG_FLETCHER_REEVES);
    TEST(OPK_NLCG_HESTENES_STIEFEL);
    TEST(OPK_NLCG_POLAK_RIBIERE_POLYAK);
    TEST(OPK_NLCG_FLETCHER);
    TEST(OPK_NLCG_LIU_STOREY);
    TEST(OPK_NLCG_DAI_YUAN);
    TEST(OPK_NLCG_PERRY_SHANNO);
    TEST(OPK_NLCG_HAGER_ZHANG);
    TEST(OPK_NLCG_POWELL);
    TEST(OPK_NLCG_SHANNO_PHUA);
    TEST(OPK_NLCG_POLAK_RIBIERE_POLYAK);

    TEST(OPK_BOUND_NONE);
    TEST(OPK_BOUND_SCALAR_FLOAT);
    TEST(OPK_BOUND_SCALAR_DOUBLE);
    TEST(OPK_BOUND_STATIC_FLOAT);
    TEST(OPK_BOUND_STATIC_DOUBLE);
    TEST(OPK_BOUND_VOLATILE_FLOAT);
    TEST(OPK_BOUND_VOLATILE_DOUBLE);
    TEST(OPK_BOUND_VECTOR);

    TEST(OPK_ASCENT);
    TEST(OPK_DESCENT);

    TEST(OPK_SCALING_NONE);
    TEST(OPK_SCALING_OREN_SPEDICATO);
    TEST(OPK_SCALING_BARZILAI_BORWEIN);

    TEST(OPK_LBFGS);
    TEST(OPK_VMLMB);
    TEST(OPK_BLMVM);

    TEST(OPK_ALGORITHM_NLCG);
    TEST(OPK_ALGORITHM_VMLMB);

    TEST(OPK_FMIN_BOUNDED_BY_A);
    TEST(OPK_FMIN_BOUNDED_BY_B);
    TEST(OPK_FMIN_SMOOTH);
    TEST(OPK_FMIN_ERROR);
    TEST(OPK_FMIN_START);
    TEST(OPK_FMIN_FX);
    TEST(OPK_FMIN_NEWX);
    TEST(OPK_FMIN_CONVERGENCE);

#undef TEST
  }
  return OPK_INVALID_ARGUMENT;
}
