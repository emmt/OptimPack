/*
 * mgh-port.h --
 *
 * Definitions for unconstrained minimization test suite from the MINPACK-1
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

#define MGH_TEST18_NMAX  50

#define MGH_SUCCESS      (0)
#define MGH_FAILURE     (-1)

extern int mgh_umck(const char** str, int n, int prob);
/* Check validity of parameters for a given MINPACK1 unconstrained minimization
   problem.  N and PROB are the size and number of the problem.  If not NULL,
   STR is used to store the description of the problem on success or an error
   message on failure.  A standard result is returned: MGH_SUCCESS or
   MGH_FAILURE. */

extern void mgh_umipt(int n, double x[], int prob, double factor);
/* Stores the standard starting points for the functions defined by subroutine
   mgh_umobj.  Argument X is a vector of N elements, X is a multiple times
   FACTOR of the standard starting point.  For the seventh function the
   standard starting point is 0.0, so in this case, if FACTOR is not unity,
   then the function returns X filled with FACTOR.  PROB has the same meaning
   as in mgh_umobj. */

extern double mgh_umobj(int n, const double x[], int prob);
/* Returns the objective functions of eighteen nonlinear unconstrained
   minimization problems.  X is the parameter array: a vector of length N, PROB
   is the problem number (a positive integer between 1 and 18).

   The values of N for functions 1,2,3,4,5,10,11,12,16 and 17 are
   3,6,3,2,3,2,4,3,2 and 4, respectively.  For function 7, N may be 2 or
   greater but is usually 6 or 9.  For functions 6,8,9,13,14,15 and 18, N may
   be variable, however it must be even for function 14, a multiple of 4 for
   function 15, and not greater than 50 for function 18. */

extern void mgh_umgrd(int n, const double x[], double g[], int prob);
/* Computes the gradient vectors of eighteen nonlinear unconstrained
   minimization problems. The problem dimensions are as described in the
   prologue comments of mgh_umobj. */

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
