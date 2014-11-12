/*
 * mgh-port.c --
 *
 * Implementation of unconstrained minimization test suite from the MINPACK-1
 * project.  Note: "MGH" stands for "Moré, Garbow and Hillstrom".
 *
 * References:
 * [1] J. J. Moré, B. S. Garbow and K. E. Hillstrom, "Testing unconstrained
 *     optimization software," ACM Trans. Math. Software 7, 17-41 (1981).
 * [2] J. J. Moré, B. S. Garbow and K. E. Hillstrom, "Fortran subroutines for
 *     testing unconstrained optimization software," ACM Trans. Math. Software
 *     7, 136-140 (1981).
 *
 * History:
 *  - Argonne National Laboratory. MINPACK-1 Project.  March 1980.
 *    Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. Moré.
 *  - Conversion to Yorick.  November 2001.  Éric Thiébaut.
 *  - Conversion to C.  February 2014. Éric Thiébaut.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2014 Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>.
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
#include <math.h>

#include "mgh-port.h"

#define MALLOC(type, number)  ((type*)malloc((number)*sizeof(type)))

#define JOIN2_(a,b)           a##b
#define JOIN3_(a,b,c)         a##b##c
#define JOIN4_(a,b,c,d)       a##b##c##d

#define JOIN(a,b)             JOIN2_(a,b)
#define JOIN2(a,b)            JOIN2_(a,b)
#define JOIN3(a,b,c)          JOIN3_(a,b,c)
#define JOIN4(a,b,c,d)        JOIN4_(a,b,c,d)

#define STRINGIFY_(x)         #x
#define STRINGIFY(x)          STRINGIFY_(x)

static int
check(int success, const char** str, const char* descr,
      const char* reason)
{
  if (success) {
    if (str != NULL) *str = descr;
    return MGH_SUCCESS;
  } else {
    if (str != NULL) *str = reason;
    return MGH_FAILURE;
  }
}

int
mgh_umck(const char** str, int n, int prob)
{
  int status = MGH_SUCCESS;

  /* Check compatibility of arguments N and PROB. */
  switch (prob) {
  case 1:
    status = check(n == 3, str, "Helical valley function.",
                   "N must be 3 for problem #1");
    break;
  case 2:
    status = check(n == 6, str, "Biggs exp6 function.",
                   "N must be 6 for problem #2");
    break;
  case 3:
    status = check(n == 3, str, "Gaussian function.",
                   "N must be 3 for problem #3");
    break;
  case 4:
    status = check(n == 2, str, "Powell badly scaled function.",
                   "N must be 2 for problem #4");
    break;
  case 5:
    status = check(n == 3, str, "Box 3-dimensional function.",
                   "N must be 3 for problem #5");
    break;
  case 6:
    status = check(n >= 1, str, "Variably dimensioned function.",
                   "N must be >= 1 in problem #6");
    break;
  case 7:
    status = check(n >= 2, str, "Watson function.",
                   "N may be 2 or greater but is usually 6 or 9 for problem #7");
    break;
  case 8:
    status = check(n >= 1, str, "Penalty function I.",
                   "N must be >= 1 in problem #8");
    break;
  case 9:
    status = check(n >= 1, str, "Penalty function II.",
                   "N must be >= 1 in problem #9");
    break;
  case 10:
    status = check(n == 2, str, "Brown badly scaled function.",
                   "N must be 2 for problem #10");
    break;
  case 11:
    status = check(n == 4, str, "Brown and Dennis function.",
                   "N must be 4 for problem #11");
    break;
  case 12:
    status = check(n == 3, str, "Gulf research and development function.",
                   "N must be 3 for problem #12");
    break;
  case 13:
    status = check(n >= 1, str, "Trigonometric function.",
                   "N must be >= 1 in problem #13");
    break;
  case 14:
    status = check(n >= 1 && n%2 == 0, str, "Extended Rosenbrock function.",
                   "N must be a multiple of 2 in problem #14");
    break;
  case 15:
    status = check(n >= 1 && n%4 == 0, str, "Extended Powell function.",
                   "N must be a multiple of 4 in problem #15");
    break;
  case 16:
    status = check(n == 2, str, "Beale function.",
                   "N must be 2 for problem #16");
    break;
  case 17:
    status = check(n == 4, str, "Wood function.",
                   "N must be 4 for problem #17");
    break;
  case 18:
    status = check(n >= 1 && n <= MGH_TEST18_NMAX, str,
                   "Chebyquad function.",
                   "N must be <= " STRINGIFY(MGH_TEST18_NMAX) \
                   " for problem #18");
    break;
  default:
    status = MGH_FAILURE;
    if (str != NULL) {
      *str = "PROB must be an integer between 1 and 18";
    }
  }
  return status;
}

/*---------------------------------------------------------------------------*/

static const double um_y[] = {9.0e-4,   4.4e-3,   1.75e-2,  5.4e-2,   1.295e-1,
                              2.42e-1,  3.521e-1, 3.989e-1, 3.521e-1, 2.42e-1,
                              1.295e-1, 5.4e-2,   1.75e-2,  4.4e-3,   9.0e-4};
void
mgh_umipt(int n, double x[], int prob, double factor)
{
  double h;
  int j;

  /* Selection of initial point. */
  switch (prob) {
  case 1:
    /* Helical valley function. */
    x[0] = -1.0;
    x[1] = 0.0;
    x[2] = 0.0;
    break;
  case 2:
    /* Biggs exp6 function. */
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 1.0;
    x[3] = 1.0;
    x[4] = 1.0;
    x[5] = 1.0;
    break;
  case 3:
    /* Gaussian function. */
    x[0] = 0.4;
    x[1] = 1.0;
    x[2] = 0.0;
    break;
  case 4:
    /* Powell badly scaled function. */
    x[0] = 0.0;
    x[1] = 1.0;
    break;
  case 5:
    /* Box 3-dimensional function. */
    x[0] = 0.0;
    x[1] = 10.0;
    x[2] = 20.0;
    break;
  case 6:
    /* Variably dimensioned function. */
    h = 1.0/((double)n);
    for (j = 0; j < n; ++j) {
      x[j] = 1.0 - ((double)(j+1))*h;
    }
    break;
  case 7:
    /* Watson function. */
    for (j = 0; j < n; ++j) {
      x[j] = 0.0;
    }
    break;
  case 8:
    /* Penalty function I. */
    for (j = 0; j < n; ++j) {
      x[j] = ((double)(j+1));
    }
    break;
  case 9:
    /* Penalty function II. */
    for (j = 0; j < n; ++j) {
      x[j] = 0.5;
    }
    break;
  case 10:
    /* Brown badly scaled function. */
    x[0] = 1.0;
    x[1] = 1.0;
    break;
  case 11:
    /* Brown and Dennis function. */
    x[0] = 25.0;
    x[1] =  5.0;
    x[2] = -5.0;
    x[3] = -1.0;
    break;
  case 12:
    /* Gulf research and development function. */
    x[0] = 5.0;
    x[1] = 2.5;
    x[2] = 0.15;
    break;
  case 13:
    /* Trigonometric function. */
    h = 1.0/((double)n);
    for (j = 0; j < n; ++j) {
      x[j] = h;
    }
    break;
  case 14:
    /* Extended Rosenbrock function. */
    for (j = 0; j < n; j += 2) {
      x[j]   = -1.2;
      x[j+1] =  1.0;
    }
    break;
  case 15:
    /* Extended Powell singular function. */
    for (j = 0; j < n; j += 4) {
      x[j]   =  3.0;
      x[j+1] = -1.0;
      x[j+2] =  0.0;
      x[j+3] =  1.0;
    }
    break;
  case 16:
    /* Beale function. */
    x[0] = 1.0;
    x[1] = 1.0;
    break;
  case 17:
    /* Wood function. */
    x[0] = -3.0;
    x[1] = -1.0;
    x[2] = -3.0;
    x[3] = -1.0;
    break;
  case 18:
    /* Chebyquad function. */
    h = 1.0/((double)(n + 1));
    for (j = 0; j < n; ++j) {
      x[j] = ((double)(j+1))*h;
    }
  }

  /* Compute multiple of initial point. */
  if (factor != 1.0) {
    if (prob == 7) {
      for (j = 0; j < n; ++j) {
        x[j] = factor;
      }
    } else {
      for (j = 0; j < n; ++j) {
        x[j] *= factor;
      }
    }
  }
}

double
mgh_umobj(int n, const double x[], int prob)
{
  const double TPI = 2.0*M_PI;
  const double AP = 1e-5, BP = 1.0;
  double f, arg, d1, d2, r, s1, s2, s3, t, t1, t2, t3, th;
  double fvec[MGH_TEST18_NMAX];
  int i, j;
  int iev;

  /* Function routine selector. */
  switch (prob) {
  case 1:
    /* Helical valley function. */
    if      (x[0] > 0.0) th = atan(x[1]/x[0])/TPI;
    else if (x[0] < 0.0) th = atan(x[1]/x[0])/TPI + 0.5;
    else                 th = (x[1] >= 0.0 ? 0.25 : -0.25);
    arg = x[0]*x[0] + x[1]*x[1];
    r = sqrt(arg) - 1.0;
    t = x[2] - 10.0*th;
    f = 100.0*(t*t + r*r) + x[2]*x[2];
    break;
  case 2:
    /* Biggs exp6 function. */
    f = 0.0;
    for (i = 1; i <= 13; ++i) {
      d1 = ((double)i)/10.0;
      d2 = exp(-d1) - 5.0*exp(-10.0*d1) + 3.0*exp(-4.0*d1);
      s1 = exp(-d1*x[0]);
      s2 = exp(-d1*x[1]);
      s3 = exp(-d1*x[4]);
      t = x[2]*s1 - x[3]*s2 + x[5]*s3 - d2;
      f += t*t;
    }
    break;
  case 3:
    /* Gaussian function. */
    f = 0.0;
    for (i = 0; i < 15; ++i) {
      d1 = 0.5*((double)i);
      d2 = 3.5 - d1 - x[2];
      arg = -0.5*x[1]*(d2*d2);
      r = exp(arg);
      t = x[0]*r - um_y[i];
      f += t*t;
    }
    break;
  case 4:
    /* Powell badly scaled function. */
    t1 = 1e4*x[0]*x[1] - 1.0;
    s1 = exp(-x[0]);
    s2 = exp(-x[1]);
    t2 = s1 + s2 - 1.0001;
    f = t1*t1 + t2*t2;
    break;
  case 5:
    /* Box 3-dimensional function. */
    f = 0.0;
    for (i = 1; i <= 10; ++i) {
      d1 = (double)i;
      d2 = d1/10.0;
      s1 = exp(-d2*x[0]);
      s2 = exp(-d2*x[1]);
      s3 = exp(-d2) - exp(-d1);
      t = s1 - s2 - s3*x[2];
      f += t*t;
    }
    break;
  case 6:
    /* Variably dimensioned function. */
    t1 = 0.0;
    t2 = 0.0;
    for (j = 0; j < n; ++j) {
      t1 += ((double)(j + 1))*(x[j] - 1.0);
      t = x[j] - 1.0;
      t2 += t*t;
    }
    t = t1*t1;
    f = t2 + t*(1.0 + t);
    break;
  case 7:
    /* Watson function. */
    f = 0.0;
    for (i = 1; i <= 29; ++i) {
      d1 = ((double)i)/29.0;
      s1 = 0.0;
      d2 = 1.0;
      for (j = 1; j < n; ++j) {
        s1 += ((double)j)*d2*x[j];
        d2 *= d1;
      }
      s2 = 0.0;
      d2 = 1.0;
      for (j = 0; j < n; ++j) {
        s2 += d2*x[j];
        d2 *= d1;
      }
      t = s1 - s2*s2 - 1.0;
      f += t*t;
    }
    t = x[0]*x[0];
    t1 = x[1] - t - 1.0;
    f += t + t1*t1;
    break;
  case 8:
    /* Penalty function I. */
    t1 = -0.25;
    t2 = 0.0;
    for (j = 0; j < n; ++j) {
      t1 += x[j]*x[j];
      t = x[j] - 1.0;
      t2 += t*t;
    }
    f = AP*t2 + BP*(t1*t1);
    break;
  case 9:
    /* Penalty function II. */
    t1 = -1.0;
    t2 =  0.0;
    t3 =  0.0;
    d1 =  exp(0.1);
    d2 =  1.0;
    s2 =  0.0; /* avoid compiler warning about s2 used uninitialized */
    for (j = 0; j < n; ++j) {
      t1 += ((double)(n - j))*(x[j]*x[j]);
      s1 = exp(x[j]/10.0);
      if (j > 0) {
        s3 = s1 + s2 - d2*(d1 + 1.0);
        t2 += s3*s3;
        t = s1 - 1.0/d1;
        t3 += t*t;
      }
      s2 = s1;
      d2 *= d1;
    }
    t = x[0] - 0.2;
    f = AP*(t2 + t3) + BP*(t1*t1 + t*t);
    break;
  case 10:
    /* Brown badly scaled function. */
    t1 = x[0] - 1e6;
    t2 = x[1] - 2e-6;
    t3 = x[0]*x[1] - 2.0;
    f = t1*t1 + t2*t2 + t3*t3;
    break;
  case 11:
    /* Brown and Dennis function. */
    f = 0.0;
    for (i = 1; i <= 20; ++i) {
      d1 = ((double)i)/5.0;
      d2 = sin(d1);
      t1 = x[0] + d1*x[1] - exp(d1);
      t2 = x[2] + d2*x[3] - cos(d1);
      t = t1*t1 + t2*t2;
      f += t*t;
    }
    break;
  case 12:
    /* Gulf research and development function. */
    f = 0.0;
    d1 = 2.0/3.0;
    for (i = 1; i <= 99; ++i) {
      arg = ((double)i)/100.0;
      r = fabs(pow(-50.0*log(arg), d1) + 25.0 - x[1]);
      t1 = pow(r, x[2])/x[0];
      t2 = exp(-t1);
      t = t2 - arg;
      f += t*t;
    }
    break;
  case 13:
    /* Trigonometric function. */
    s1 = 0.0;
    for (j = 0; j < n; ++j) {
      s1 += cos(x[j]);
    }
    f = 0.0;
    for (j = 0; j < n; ++j) {
      t = ((double)(n + 1 + j)) - sin(x[j]) - s1
        - ((double)(1 + j))*cos(x[j]);
      f += t*t;
    }
    break;
  case 14:
    /* Extended Rosenbrock function. */
    f = 0.0;
    for (j = 0; j < n; j += 2) {
      t1 = 1.0 - x[j];
      t2 = 10.0*(x[j+1] - x[j]*x[j]);
      f += t1*t1 + t2*t2;
    }
    break;
  case 15:
    /* Extended Powell function. */
    f = 0.0;
    for (j = 0; j < n; j += 4) {
      t = x[j] + 10.0*x[j+1];
      t1 = x[j+2] - x[j+3];
      s1 = 5.0*t1;
      t2 = x[j+1] - 2.0*x[j+2];
      s2 = t2*t2*t2;
      t3 = x[j] - x[j+3];
      s3 = 10.0*(t3*t3*t3);
      f += t*t + s1*t1 + s2*t2 + s3*t3;
    }
    break;
  case 16:
    /* Beale function. */
    s1 = 1.0 - x[1];
    t1 = 1.5 - x[0]*s1;
    s2 = 1.0 - x[1]*x[1];
    t2 = 2.25 - x[0]*s2;
    s3 = 1.0 - x[1]*x[1]*x[1];
    t3 = 2.625 - x[0]*s3;
    f = t1*t1 + t2*t2 + t3*t3;
    break;
  case 17:
    /* Wood function. */
    s1 = x[1] - x[0]*x[0];
    s2 = 1.0 - x[0];
    s3 = x[1] - 1.0;
    t1 = x[3] - x[2]*x[2];
    t2 = 1.0 - x[2];
    t3 = x[3] - 1.0;
    d1 = s3 + t3;
    d2 = s3 - t3;
    f = (100.0*(s1*s1) + (s2*s2) + 90.0*(t1*t1) + (t2*t2)
         + 10.0*(d1*d1) + (d2*d2)/10.0);
    break;
  case 18:
    /* Chebyquad function. */
    for (i = 0; i < n; ++i) {
      fvec[i] = 0.0;
    }
    for (j = 0; j < n; ++j) {
      t1 = 1.0;
      t2 = 2.0*x[j] - 1.0;
      t = 2.0*t2;
      for (i = 0; i < n; ++i) {
        fvec[i] += t2;
        th = t*t2 - t1;
        t1 = t2;
        t2 = th;
      }
    }
    f = 0.0;
    d1 = 1.0/((double) n);
    iev = -1;
    for (i = 0; i < n; ++i) {
      t = d1*fvec[i];
      if (iev > 0) {
        r = (double)(i + 1);
        t += 1.0/(r*r - 1.0);
      }
      f += t*t;
      iev = -iev;
    }
    break;
  default:
    f = 0.0;
  }
  return f;
}

void
mgh_umgrd(int n, const double x[], double g[], int prob)
{
  const double TPI = 2.0*M_PI;
  const double AP = 1e-5, BP = 1.0;
  double arg, d1, d2, r, s1, s2, s3, t, t1, t2, t3, th;
  double fvec[MGH_TEST18_NMAX];
  int i, j;
  int iev;

  /* Gradient routine selector. */
  switch (prob) {
  case 1:
    /* Helical valley function. */
    if      (x[0] > 0.0) th = atan(x[1]/x[0])/TPI;
    else if (x[0] < 0.0) th = atan(x[1]/x[0])/TPI + 0.5;
    else                 th = (x[1] >= 0.0 ? 0.25 : -0.25);
    arg = x[0]*x[0] + x[1]*x[1];
    r = sqrt(arg);
    t = x[2] - 10.0*th;
    s1 = 10.0*t/(TPI*arg);
    g[0] = 200.0*(x[0] - x[0]/r + x[1]*s1);
    g[1] = 200.0*(x[1] - x[1]/r - x[0]*s1);
    g[2] = 2.0*(100.0*t + x[2]);
    break;
  case 2:
    /* Biggs exp6 function. */
    for (j = 0; j < 6; ++j) g[j] = 0.0;
    for (i = 1; i <= 13; ++i) {
      d1 = ((double)i)/10.0;
      d2 = exp(-d1) - 5.0*exp(-10.0*d1) + 3.0*exp(-4.0*d1);
      s1 = exp(-d1*x[0]);
      s2 = exp(-d1*x[1]);
      s3 = exp(-d1*x[4]);
      t = x[2]*s1 - x[3]*s2 + x[5]*s3 - d2;
      th = d1*t;
      g[0] -= s1*th;
      g[1] += s2*th;
      g[2] += s1*t;
      g[3] -= s2*t;
      g[4] -= s3*th;
      g[5] += s3*t;
    }
    g[0] *= 2.0*x[2];
    g[1] *= 2.0*x[3];
    g[2] *= 2.0;
    g[3] *= 2.0;
    g[4] *= 2.0*x[5];
    g[5] *= 2.0;
    break;
  case 3:
    /* Gaussian function. */
    g[0] = 0.0;
    g[1] = 0.0;
    g[2] = 0.0;
    for (i = 0; i < 15; ++i) {
      d1 = 0.5*((double)i);
      d2 = 3.5 - d1 - x[2];
      arg = -0.5*x[1]*(d2*d2);
      r = exp(arg);
      t = x[0]*r - um_y[i];
      s1 = r*t;
      s2 = d2*s1;
      g[0] += s1;
      g[1] -= d2*s2;
      g[2] += s2;
    }
    g[0] *= 2.0;
    g[1] *= x[0];
    g[2] *= 2.0*x[0]*x[1];
    break;
  case 4:
    /* Powell badly scaled function. */
    t1 = 1e4*x[0]*x[1] - 1.0;
    s1 = exp(-x[0]);
    s2 = exp(-x[1]);
    t2 = s1 + s2 - 1.0001;
    g[0] = 2.0*(1e4*x[1]*t1 - s1*t2);
    g[1] = 2.0*(1e4*x[0]*t1 - s2*t2);
    break;
  case 5:
    /* Box 3-dimensional function. */
    g[0] = 0.0;
    g[1] = 0.0;
    g[2] = 0.0;
    for (i = 1; i <= 10; ++i) {
      d1 = (double)i;
      d2 = d1/10.0;
      s1 = exp(-d2*x[0]);
      s2 = exp(-d2*x[1]);
      s3 = exp(-d2) - exp(-d1);
      t = s1 - s2 - s3*x[2];
      th = d2*t;
      g[0] -= s1*th;
      g[1] += s2*th;
      g[2] -= s3*t;
    }
    g[0] *= 2.0;
    g[1] *= 2.0;
    g[2] *= 2.0;
    break;
  case 6:
    /* Variably dimensioned function. */
    t1 = 0.0;
    for (j = 0; j < n; ++j) {
      t1 += ((double)(j+1))*(x[j] - 1.0);
    }
    t = t1*(1.0 + 2.0*(t1*t1));
    for (j = 0; j < n; ++j) {
      g[j] = 2.0*(x[j] - 1.0 + ((double)(j+1))*t);
    }
    break;
  case 7:
    /* Watson function. */
    for (j = 0; j < n; ++j) {
      g[j] = 0.0;
    }
    for (i = 1; i <= 29; ++i) {
      d1 = ((double)i)/29.0;
      s1 = 0.0;
      d2 = 1.0;
      for (j = 1; j < n; ++j) {
        s1 += ((double)j)*d2*x[j];
        d2 *= d1;
      }
      s2 = 0.0;
      d2 = 1.0;
      for (j = 0; j < n; ++j) {
        s2 += d2*x[j];
        d2 *= d1;
      }
      t = s1 - s2*s2 - 1.0;
      s3 = 2.0*d1*s2;
      d2 = 2.0/d1;
      for (j = 0; j < n; ++j) {
        g[j] += d2*(((double)j) - s3)*t;
        d2 *= d1;
      }
    }
    t1 = x[1] - x[0]*x[0] - 1.0;
    g[0] += x[0]*(2.0 - 4.0*t1);
    g[1] += 2.0*t1;
    break;
  case 8:
    /* Penalty function I. */
    t1 = -0.25;
    for (j = 0; j < n; ++j) {
      t1 += x[j]*x[j];
    }
    d1 = 2.0*AP;
    th = 4.0*BP*t1;
    for (j = 0; j < n; ++j) {
      g[j] = d1*(x[j] - 1.0) + x[j]*th;
    }
    break;
  case 9:
    /* Penalty function II. */
    s2 = 0.0; /* avoid compiler warning about s2 used uninitialized */
    t1 = -1.0;
    for (j = 0; j < n; ++j) {
      t1 += ((double)(n - j))*(x[j]*x[j]);
    }
    d1 = exp(0.1);
    d2 = 1.0;
    th = 4.0*BP*t1;
    for (j = 0; j < n; ++j) {
      g[j] = ((double)(n - j))*x[j]*th;
      s1 = exp(x[j]/10.0);
      if (j > 0) {
        s3 = s1 + s2 - d2*(d1 + 1.0);
        g[j] += AP*s1*(s3 + s1 - 1.0/d1)/5.0;
        g[j-1] += AP*s2*s3/5.0;
      }
      s2 = s1;
      d2 *= d1;
    }
    g[0] += 2.0*BP*(x[0] - 0.2);
    break;
  case 10:
    /* Brown badly scaled function. */
    t1 = x[0] - 1e6;
    t2 = x[1] - 2e-6;
    t3 = x[0]*x[1] - 2.0;
    g[0] = 2.0*(t1 + x[1]*t3);
    g[1] = 2.0*(t2 + x[0]*t3);
    break;
  case 11:
    /* Brown and Dennis function. */
    g[0] = 0.0;
    g[1] = 0.0;
    g[2] = 0.0;
    g[3] = 0.0;
    for (i = 1; i <= 20; ++i) {
      d1 = ((double)i)/5.0;
      d2 = sin(d1);
      t1 = x[0] + d1*x[1] - exp(d1);
      t2 = x[2] + d2*x[3] - cos(d1);
      t = t1*t1 + t2*t2;
      s1 = t1*t;
      s2 = t2*t;
      g[0] += s1;
      g[1] += d1*s1;
      g[2] += s2;
      g[3] += d2*s2;
    }
    g[0] *= 4.0;
    g[1] *= 4.0;
    g[2] *= 4.0;
    g[3] *= 4.0;
    break;
  case 12:
    /* Gulf research and development function. */
    g[0] = 0.0;
    g[1] = 0.0;
    g[2] = 0.0;
    d1 = 2.0/3.0;
    for (i = 1; i <= 99; ++i) {
      arg = ((double)i)/100.0;
      r = fabs(pow(-50.0*log(arg), d1) + 25.0 - x[1]);
      t1 = pow(r, x[2])/x[0];
      t2 = exp(-t1);
      t = t2 - arg;
      s1 = t1*t2*t;
      g[0] += s1;
      g[1] += s1/r;
      g[2] -= s1*log(r);
    }
    g[0] *= 2.0/x[0];
    g[1] *= 2.0*x[2];
    g[2] *= 2.0;
    break;
  case 13:
    /* Trigonometric function. */
    s1 = 0.0;
    for (j = 0; j < n; ++j) {
      g[j] = cos(x[j]);
      s1 += g[j];
    }
    s2 = 0.0;
    for (j = 0; j < n; ++j) {
      th = sin(x[j]);
      t = ((double)(n+1+j)) - th - s1 - ((double)(1+j))*g[j];
      s2 += t;
      g[j] = (((double)(1+j))*th - g[j])*t;
    }
    for (j = 0; j < n; ++j) {
      g[j] = 2.0*(g[j] + sin(x[j])*s2);
    }
    break;
  case 14:
    /* Extended Rosenbrock function. */
    for (j = 0; j < n; j += 2) {
      t1 = 1.0 - x[j];
      g[j+1] = 200.0*(x[j+1] - x[j]*x[j]);
      g[j] = -2.0*(x[j]*g[j+1] + t1);
    }
    break;
  case 15:
    /* Extended Powell function. */
    for (j = 0; j < n; j += 4) {
      t = x[j] + 10.0*x[j+1];
      t1 = x[j+2] - x[j+3];
      s1 = 5.0*t1;
      t2 = x[j+1] - 2.0*x[j+2];
      s2 = 4.0*(t2*t2*t2);
      t3 = x[j] - x[j+3];
      s3 = 20.0*(t3*t3*t3);
      g[j] = 2.0*(t + s3);
      g[j+1] = 20.0*t + s2;
      g[j+2] = 2.0*(s1 - s2);
      g[j+3] = -2.0*(s1 + s3);
    }
    break;
  case 16:
    /* Beale function. */
    s1 = 1.0 - x[1];
    t1 = 1.5 - x[0]*s1;
    s2 = 1.0 - x[1]*x[1];
    t2 = 2.25 - x[0]*s2;
    s3 = 1.0 - x[1]*x[1]*x[1];
    t3 = 2.625 - x[0]*s3;
    g[0] = -2.0*(s1*t1 + s2*t2 + s3*t3);
    g[1] = 2.0*x[0]*(t1 + x[1]*(2.0*t2 + 3.0*x[1]*t3));
    break;
  case 17:
    /* Wood function. */
    s1 = x[1] - x[0]*x[0];
    s2 = 1.0 - x[0];
    s3 = x[1] - 1.0;
    t1 = x[3] - x[2]*x[2];
    t2 = 1.0 - x[2];
    t3 = x[3] - 1.0;
    g[0] = -2.0*(200.0*x[0]*s1 + s2);
    g[1] = 200.0*s1 + 20.2*s3 + 19.8*t3;
    g[2] = -2.0*(180.0*x[2]*t1 + t2);
    g[3] = 180.0*t1 + 20.2*t3 + 19.8*s3;
    break;
  case 18:
    /* Chebyquad function. */
    for (i = 0; i < n; ++i) {
      fvec[i] = 0.0;
    }
    for (j = 0; j < n; ++j) {
      t1 = 1.0;
      t2 = 2.0*x[j] - 1.0;
      t = 2.0*t2;
      for (i = 0; i < n; ++i) {
        fvec[i] += t2;
        th = t*t2 - t1;
        t1 = t2;
        t2 = th;
      }
    }
    d1 = 1.0/((double)n);
    iev = -1;
    for (i = 0; i < n; ++i) {
      fvec[i] *= d1;
      if (iev > 0) {
        r = (double)(i + 1);
        fvec[i] += 1.0/(r*r - 1.0);
      }
      iev = -iev;
    }
    for (j = 0; j < n; ++j) {
      g[j] = 0.0;
      t1 = 1.0;
      t2 = 2.0*x[j] - 1.0;
      t = 2.0*t2;
      s1 = 0.0;
      s2 = 2.0;
      for (i = 0; i < n; ++i) {
        g[j] += fvec[i]*s2;
        th = 4.0*t2 + t*s2 - s1;
        s1 = s2;
        s2 = th;
        th = t*t2 - t1;
        t1 = t2;
        t2 = th;
      }
    }
    d2 = 2.0*d1;
    for (j = 0; j < n; ++j) {
      g[j] *= d2;
    }
    break;
  default:
    for (j = 0; j < n; ++j) {
      g[j] = 0.0;
    }
  }
}

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: nil
 * c-basic-offset: 2
 * fill-column: 79
 * coding: utf-8
 * End:
 */
