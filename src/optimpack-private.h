/*
 * optimpack-private.h --
 *
 * Private macro definitions for building OptimPack library.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2014 Éric Thiébaut <thiebaut@obs.univ-lyon1.fr>
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can use, modify
 * and/or redistribute the software under the terms of the CeCILL-C license as
 * circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty and the software's author, the holder of the
 * economic rights, and the successive licensors have only limited liability.
 *
 * In this respect, the user's attention is drawn to the risks associated with
 * loading, using, modifying and/or developing or reproducing the software by
 * the user in light of its specific status of free software, that may mean
 * that it is complicated to manipulate, and that also therefore means that it
 * is reserved for developers and experienced professionals having in-depth
 * computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling
 * the security of their systems and/or data to be ensured and, more
 * generally, to use and operate it in the same conditions as regards
 * security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 *
 *-----------------------------------------------------------------------------
 */

#ifndef _OPTIMPACK_PRIVATE_H
#define _OPTIMPACK_PRIVATE_H 1

#include "optimpack.h"

#define TRUE        1
#define FALSE       0

#define MAX(a,b)   ((a) >= (b) ? (a) : (b))
#define MIN(a,b)   ((a) <= (b) ? (a) : (b))

#define ROUND_UP(a, b)        ((((a) + ((b) - 1))/(b))*(b))
#define NEW(TYPE)             ((TYPE*)calloc(1, sizeof(TYPE)))

/* Helper macros for simple loops. */
#define LOOP_0(var,num)    for (var = 0; var < num; ++var)  /* C-like */
#define LOOP_1(var,num)    for (var = 1; var <= num; ++var) /* FORTRAN-like */

/*---------------------------------------------------------------------------*/
/* PRIVATE API FOR LINE SEARCH */

/* The base structure for line search must be exposed for line search
   methods. */
struct _opk_lnsrch_workspace {
  double stp;    /* current step length */
  double stpmin; /* lower bound for the step */
  double stpmax; /* upper bound for the step */
  double finit; /* function value at the start of the search */
  double ginit; /* directional derivative value at the start of the search */

  /* Method called by opk_lnsrch_start() to initiate a line search.  This
     method is called after having set member STP with the length of the first
     step to perform, members STPMIN and STPMAX with the lower and upper bounds
     of the step length, members FINIT and GINIT with the function value and
     the derivative the function along the search direction at the start of the
     search. */
  int (*start)(opk_lnsrch_workspace_t* ws);

  /* Method to iterate during a line search.  STP_PTR, F_PTR and D_PTR are the
     addresses of variables which store STP, F and D.  On entry, STP is the
     current step length, F is the function value for this step and D is the
     corresponding derivative of the function along the search direction.
     FIXME: On exit, if convergence is achieved (or in case of error/warning)
     STP, F and D are left unchanged; otherwise STP is the new step to try. */
  int (*iterate)(opk_lnsrch_workspace_t* ws,
		 double* stp_ptr, double f1, double d1);

  /* Method used to release any ressources speciffically allocated by the
     line search method (not the workspace itself). */
  void (*delete)(opk_lnsrch_workspace_t* ws);

  int status; /* last value returned by line search methods */
  int searching; /* true if search is in progress */
};

extern opk_lnsrch_workspace_t*
_opk_lnsrch_new(size_t size,
		int (*start)(opk_lnsrch_workspace_t* ws),
		int (*iterate)(opk_lnsrch_workspace_t* ws,
			       double* stp_ptr, double f1, double g1),
	        void (*delete)(opk_lnsrch_workspace_t* ws));

#endif /* _OPTIMPACK_PRIVATE_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * c-basic-offset: 2
 * fill-column: 79
 * coding: utf-8
 * ispell-local-dictionary: "american"
 * End:
 */
