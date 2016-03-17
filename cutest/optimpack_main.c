/* ====================================================
 * CUTEst interface for OptimPack
 *
 * Éric Thiébaut, November 12th, 2015
 *
 * (Based on CUTEr genc_main.c of D. Orban, Feb 3, 2003)
 * (CUTEst evolution, Nick Gould, Apr 2, 2014)
 * ====================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <time.h>

#define DESCRIPTION_MAX_SIZE 128
#define MAXLINE 512
#define VERBOSE 0

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

#include "cutest.h"
#include "optimpack.h"

typedef enum {
  NLCG = 0, VMLM, VMLMB
} algorithm_t;

static void getinfo(integer n, integer m, doublereal *bl, doublereal *bu,
                    doublereal *cl, doublereal *cu, logical *equatn,
                    logical *linear, VarTypes *vartypes)
{
  int i;

  vartypes->nlin = 0;
  vartypes->neq = 0;
  vartypes->nbnds = 0;
  vartypes->nrange = 0;
  vartypes->nlower = 0;
  vartypes->nupper = 0;
  vartypes->nineq = 0;
  vartypes->nineq_lin = 0;
  vartypes->nineq_nlin = 0;
  vartypes->neq_lin = 0;
  vartypes->neq_nlin = 0;

  for (i = 0; i < n; ++i) {
    if (bl[i] > -CUTE_INF || bu[i] < CUTE_INF) {
      ++vartypes->nbnds;
    }
  }
  for (i = 0; i < m; ++i) {
    if (linear[i]) {
      ++vartypes->nlin;
    }
    if (equatn[i]) {
      ++vartypes->neq;
      if (linear[i]) {
        ++vartypes->neq_lin;
      } else {
        ++vartypes->neq_nlin;
      }
    } else {
      if (cl[i] > -CUTE_INF) {
        if (cu[i] < CUTE_INF) {
          vartypes->nrange++;
        } else {
          vartypes->nlower++;
          vartypes->nineq++;
        }
      } else {
        if (cu[i] < CUTE_INF) {
          ++vartypes->nupper;
          ++vartypes->nineq;
        }
      }
      if (!equatn[i] && linear[i]) {
        ++vartypes->nineq_lin;
      } else {
        ++vartypes->nineq_nlin;
      }
    }
  }
}

/* Return the (projected) gradient norm as in CG_DESCENT or ASA_CG */
static double
gradnorm(const opk_vector_t* x, const opk_vector_t* g, opk_vector_t* gp,
         const opk_convexset_t* box)
{
  if (box != NULL) {
    /* compute ||x_k - Projection(x_k - g_k)||_infty */
    opk_vaxpby(gp, 1, x, -1, g);
    if (opk_project_variables(gp, gp, box) != OPK_SUCCESS) {
      printf("# Failed to project variables\n");
      exit(-1);
    }
    opk_vaxpby(gp, 1, x, -1, gp);
    return opk_vnorminf(gp);
  } else {
    return opk_vnorminf(g);
  }
}

/* global variables */
integer CUTEst_nvar;	    /* number of variables */
integer CUTEst_ncon;	    /* number of constraints */

/* main program */
int MAINENTRY(void)
{
  char *fname = "OUTSDIF.d"; /* CUTEst data file */
  integer funit = 42;        /* FORTRAN unit number for OUTSDIF.d */
  integer io_buffer = 11;    /* FORTRAN unit for internal i/o */
  integer iout = 6;          /* FORTRAN unit number for error output */
  integer ierr;              /* Exit flag from OPEN and CLOSE */
  integer status;            /* Exit flag from CUTEst tools */
  VarTypes vartypes;
  integer    ncon_dummy;
  doublereal *x, *g, *bl, *bu;
  doublereal *v = NULL, *cl = NULL, *cu = NULL;
  logical *equatn = NULL, *linear = NULL;
  char *pname, *vnames, *gnames, *cptr;
  char **Vnames, **Gnames; /* vnames and gnames as arrays of strings */
  integer e_order = 1, l_order = 0, v_order = 0;
  logical     constrained = FALSE_;
  doublereal  calls[7], cpu[2];
  integer     i, j;
  FILE       *spec;

  /* Variables for OptimPack methods. */
  opk_vspace_t* vspace;
  opk_vector_t* vx;
  opk_vector_t* vg;
  opk_vector_t* vgp;
  opk_vector_t* vbl;
  opk_vector_t* vbu;
  opk_nlcg_t* nlcg  = NULL; /* non-linear conjugate gradient optimizer */
  opk_vmlmb_t* vmlmb = NULL; /* idem with bound constraints */
  opk_lnsrch_t* lnsrch = NULL; /* Line search method */
  opk_task_t task;
  opk_status_t final_status;
#define NLCG  1
#define VMLMB 2
  int algorithm = VMLMB;
  unsigned int vmlmb_flags = 0;
  unsigned int nlcg_flags = OPK_NLCG_DEFAULT;
  const char* algorithm_name;
  opk_convexset_t* box;
  void* lower_value;
  void* upper_value;
  opk_bound_type_t lower_type;
  opk_bound_type_t upper_type;
  double f, gnorm, finalf, finalgnorm, gtest;
  double gatol = 1.0e-6;    /* required gradient absolute tolerance */
  double grtol = 0.0;       /* required gradient relative tolerance */
  double delta = 5e-2;
  double epsilon = 1e-2;
  double sftol = 1e-4;
  double sgtol = 0.9;
  double sxtol = 1e-15;
  double samin = 0.1;
  int delta_given = FALSE_;
  int epsilon_given = FALSE_;
  int maxiter = -1;
  int maxeval = -1;
  int verbose = 0;
  int evaluations, iterations, bounds;
  int mem = 5; /* number of saved (Y,S) pairs memorized by the quasi-Newton
                  method to approximate the Hessian, if zero non-linear
                  conjugate gradient is used */


  /* Open problem description file OUTSDIF.d */
  ierr = 0;
  FORTRAN_open( &funit, fname, &ierr );
  if (ierr) {
    printf("# Error opening file OUTSDIF.d.\nAborting.\n");
    exit(1);
  }

  /* Determine problem size */
  CUTEST_cdimen(&status, &funit, &CUTEst_nvar, &CUTEst_ncon);
  if (status != 0) {
    printf("# CUTEst error, status = %d, aborting\n", status);
    exit(status);
  }

  /* Determine whether to call constrained or unconstrained tools */
  constrained = (CUTEst_ncon > 0);

  /* Seems to be needed for some Solaris C compilers */
  ncon_dummy = CUTEst_ncon + 1;

  /* Reserve memory for variables, gradient, bounds, and multipliers and
     call appropriate initialization routine for CUTEst */
  x  = (double*)malloc(CUTEst_nvar*sizeof(double));
  g  = (double*)malloc(CUTEst_nvar*sizeof(double));
  bl = (double*)malloc(CUTEst_nvar*sizeof(double));
  bu = (double*)malloc(CUTEst_nvar*sizeof(double));
  if (constrained) {
    MALLOC( equatn, CUTEst_ncon, logical    );
    MALLOC( linear, CUTEst_ncon, logical    );
    MALLOC( v,      CUTEst_ncon, doublereal );
    MALLOC( cl,     CUTEst_ncon, doublereal );
    MALLOC( cu,     CUTEst_ncon, doublereal );
    CUTEST_csetup(&status, &funit, &iout, &io_buffer,
                  &CUTEst_nvar, &CUTEst_ncon, x, bl, bu,
                  v, cl, cu, equatn, linear,
                  &e_order, &l_order, &v_order);
  } else {
    CUTEST_usetup(&status, &funit, &iout, &io_buffer, &CUTEst_nvar, x, bl, bu);
  }
  if (status != 0) {
    printf("# CUTEst error, status = %d, aborting\n", status);
    exit(status);
  }

  /* Get problem, variables and constraints names */
  MALLOC(pname, FSTRING_LEN + 1, char);
  MALLOC(vnames, CUTEst_nvar * FSTRING_LEN, char);        /* For Fortran */
  MALLOC(Vnames, CUTEst_nvar, char *);               /* Array of strings */
  for (i = 0; i < CUTEst_nvar; i++) {
    MALLOC(Vnames[i], FSTRING_LEN + 1, char);
  }
  if (constrained) {
    MALLOC(gnames, CUTEst_ncon * FSTRING_LEN, char);   /* For Fortran */
    MALLOC(Gnames, CUTEst_ncon, char *);          /* Array of strings */
    for (i = 0; i < CUTEst_ncon; i++) {
      MALLOC(Gnames[i], FSTRING_LEN + 1, char);
    }
    CUTEST_cnames(&status, &CUTEst_nvar, &CUTEst_ncon,
                  pname, vnames, gnames);
  } else {
    CUTEST_unames(&status, &CUTEst_nvar, pname, vnames);
  }
  if (status != 0) {
    printf("# CUTEst error, status = %d, aborting\n", status);
    exit(status);
  }

  /* Make sure to null-terminate problem name */
  i = FSTRING_LEN;
  do {
     pname[i] = '\0';
  } while (--i >= 0 && pname[i] == ' ');

  /* Transfer variables and constraint names into arrays of
   * null-terminated strings.
   * If you know of a simpler way to do this portably, let me know!
   */
  for (i = 0; i < CUTEst_nvar; i++) {
    cptr = vnames + i * FSTRING_LEN;
    for (j = 0; j < FSTRING_LEN; j++) {
      Vnames[i][j] = *cptr;
      cptr++;
    }
    Vnames[i][FSTRING_LEN] = '\0';
  }

  for (i = 0; i < CUTEst_ncon; i++) {
    cptr = gnames + i * FSTRING_LEN;
    for (j = 0; j < FSTRING_LEN; j++) {
      Gnames[i][j] = *cptr;
      cptr++;
    }
    Gnames[i][FSTRING_LEN] = '\0';
  }

  /* Fortran strings no longer needed */
  FREE(vnames);
  if (constrained) FREE(gnames);

#if VERBOSE
  printf("# Variable names:\n");
  for (i = 0; i < CUTEst_nvar; i++) {
    printf("#  %s\n", Vnames[i]);
  }
#endif

  /* Free memory for variable names */
  for (i = 0; i < CUTEst_nvar; i++) {
    FREE(Vnames[i]);
  }
  FREE(Vnames);

#if VERBOSE
  if (constrained) {
    printf("# Constraint names:\n");
    for (i = 0; i < CUTEst_ncon; i++) {
      printf("#  %s\n", Gnames[i]);
    }
  }
#endif

  if (constrained) {
    /* Free memory for constraint names */
    for (i = 0; i < CUTEst_ncon; i++) {
      FREE(Gnames[i]);
    }
    FREE(Gnames);
  }

  if (constrained) {
    printf("# the problem %s has %i constraints\n",
            pname,  CUTEst_ncon);
    printf("# current OptimPack is for unconstrained optimization\n");
    exit(-1);
  }

  /* Obtain basic info on problem */
  getinfo(CUTEst_nvar, CUTEst_ncon, bl, bu, NULL, NULL, NULL, NULL, &vartypes);
  bounds = 0;
  for (i = 0; i < CUTEst_nvar; ++i) {
    /*if (i < 10) printf("# bl[%d] = %g / bu[%d] = %g\n", i, bl[i], i, bu[i]);*/
    if (bl[i] > -CUTE_INF) {
      bounds |= 1;
    }
    if (bu[i] < CUTE_INF) {
      bounds |= 2;
    }
  }


  /*printf("# Problem: %s (n = %i)\n", pname, CUTEst_nvar );*/

  /* MALLOC(vnames, CUTEst_nvar*FSTRING_LEN, char);
     CUTEST_unames( &status, &CUTEst_nvar, pname, vnames);
     if( status ) {
     printf("# CUTEst error, status = %d, aborting\n", status);
     exit(status);
     }
     FREE(vnames);
  */

  /* Set any parameter values here
     First read in the default parameter values */

  /* Parameter values are overwritten using any values stored in the
     OPTIMPACK.SPC file. The format of the file is parameter name at
     the start of the line followed by one or more spaces and then
     the parameter value.  */

  spec = fopen("OPTIMPACK.SPC", "r");
  if (spec != NULL) {
    char line[MAXLINE+1];
    char str[MAXLINE+1];
    int ival;
    double dval;
    int powell = FALSE_;
    unsigned int autostep = 0;
    nlcg_flags = OPK_NLCG_POLAK_RIBIERE_POLYAK;
    while (fgets(line, MAXLINE, spec) != (char*)NULL) {
      /* determine the parameter and its value */
      line[MAXLINE] = 0;
      if (sscanf(line, " algorithm %s", str) == 1) {
        if (strcasecmp(str, "nlcg") == 0) {
          algorithm = NLCG;
        } else if (strcasecmp(str, "vmlmb") == 0) {
          algorithm = VMLMB;
        } else if (strcasecmp(str, "blmvm") == 0) {
          algorithm = VMLMB;
          vmlmb_flags = OPK_EMULATE_BLMVM;
        } else {
          printf("# Unknown algorithm\n");
          exit(-1);
        }
        continue;
      }
      if (sscanf(line, " linesearch %s", str) == 1) {
        if (strcasecmp(str, "quadratic") == 0) {
          lnsrch = opk_lnsrch_new_backtrack(sftol, samin);
        } else if (strcasecmp(str, "Armijo") == 0) {
          lnsrch = opk_lnsrch_new_backtrack(sftol, 0.5);
        } else if (strcasecmp(str, "cubic") == 0) {
          lnsrch = opk_lnsrch_new_csrch(sftol, sgtol, sxtol);
        } else {
          printf("# Unknown line search method\n");
          exit(-1);
        }
        continue;
      }
      if (sscanf(line, " Powell %d", &ival) == 1) {
        if (ival != 0) {
          powell = TRUE_;
        }
        continue;
      }
      if (sscanf(line, " ShannoPhua %d", &ival) == 1) {
        if (ival != 0) {
          autostep = OPK_NLCG_SHANNO_PHUA;
        }
        continue;
      }
#if 0
      if (sscanf(line, " OrenSpedicato %d", &ival) == 1) {
        if (ival != 0) {
          autostep = OPK_NLCG_OREN_SPEDICATO;
        }
        continue;
      }
      if (sscanf(line, " BarzilaiBorwein %d", &ival) == 1) {
        if (ival != 0) {
          autostep = OPK_NLCG_BARZILAI_BORWEIN;
        }
        continue;
      }
#endif
      if (sscanf(line, " FletcherReeves %d", &ival) == 1) {
        if (ival != 0) {
          nlcg_flags = OPK_NLCG_FLETCHER_REEVES;
        }
        continue;
      }
      if (sscanf(line, " HestenesStiefel %d", &ival) == 1) {
        if (ival != 0) {
          nlcg_flags = OPK_NLCG_HESTENES_STIEFEL;
        }
        continue;
      }
      if (sscanf(line, " PolakRibierePolyak %d", &ival) == 1) {
        if (ival != 0) {
          nlcg_flags = OPK_NLCG_POLAK_RIBIERE_POLYAK;
        }
        continue;
      }
      if (sscanf(line, " Fletcher %d", &ival) == 1) {
        if (ival != 0) {
          nlcg_flags = OPK_NLCG_FLETCHER;
        }
        continue;
      }
      if (sscanf(line, " LiuStorey %d", &ival) == 1) {
        if (ival != 0) {
          nlcg_flags = OPK_NLCG_LIU_STOREY;
        }
        continue;
      }
      if (sscanf(line, " DaiYuan %d", &ival) == 1) {
        if (ival != 0) {
          nlcg_flags = OPK_NLCG_DAI_YUAN;
        }
        continue;
      }
      if (sscanf(line, " PerryShanno %d", &ival) == 1) {
        if (ival != 0) {
          nlcg_flags = OPK_NLCG_PERRY_SHANNO;
        }
        continue;
      }
      if (sscanf(line, " HagerZhang %d", &ival) == 1) {
        if (ival != 0) {
          nlcg_flags = OPK_NLCG_HAGER_ZHANG;
        }
        continue;
      }
      if (sscanf(line, " mem %d", &ival) == 1) {
        if (mem <= 0) {
          printf("# Illegal value for option MEM\n");
          exit(-1);
        }
        mem = ival;
        continue;
      }
      if (sscanf(line, " maxiter %d", &ival) == 1) {
        maxiter = ival;
        continue;
      }
      if (sscanf(line, " maxeval %d", &ival) == 1) {
        maxeval = ival;
        continue;
      }
      if (sscanf(line, " verb %d", &ival) == 1) {
        verbose = ival;
        continue;
      }
      if (sscanf(line, " delta %lf", &dval) == 1) {
        delta = dval;
        delta_given = TRUE_;
        continue;
      }
      if (sscanf(line, " epsilon %lf", &dval) == 1) {
        epsilon = dval;
        epsilon_given = TRUE_;
        continue;
      }
      if (sscanf(line, " gatol %lf", &dval) == 1) {
        gatol = dval;
        continue;
      }
      if (sscanf(line, " grtol %lf", &dval) == 1) {
        grtol = dval;
        continue;
      }
    }
    if (powell) {
      nlcg_flags |= OPK_NLCG_POWELL;
    }
    if (autostep != 0) {
      nlcg_flags |= autostep;
    }
    fclose(spec);
  }
  if (gatol < 0) {
    printf("# Bad value for GATOL\n");
    exit(-1);
  }
  if (grtol < 0) {
    printf("# Bad value for GRTOL\n");
    exit(-1);
  }

  /* Check whether the problem has constraints */
  if (bounds != 0 && mem == 0) {
    printf("# Use VMLMB or BLMVM for bound constrained optimization\n");
    exit(-1);
  }

  /* Build vector space and populate it with work vectors. */
  vspace = opk_new_simple_double_vector_space(CUTEst_nvar);
  if (vspace == NULL) {
    printf("# Failed to allocate vector space\n");
    exit(-1);
  }
  vx  = opk_wrap_simple_double_vector(vspace, x,  free, x);
  vg  = opk_wrap_simple_double_vector(vspace, g,  free, g);
  vbl = opk_wrap_simple_double_vector(vspace, bl, free, bl);
  vbu = opk_wrap_simple_double_vector(vspace, bu, free, bu);
  if (vx == NULL || vg == NULL || vbl == NULL || vbu == NULL) {
    printf("# Failed to wrap vectors\n");
    exit(-1);
  }
  if ((bounds & 1) == 0) {
    lower_value = NULL;
    lower_type = OPK_BOUND_NONE;
  } else {
    lower_value = vbl;
    lower_type = OPK_BOUND_VECTOR;
  }
  if ((bounds & 2) == 0) {
    upper_value = NULL;
    upper_type = OPK_BOUND_NONE;
  } else {
    upper_value = vbu;
    upper_type = OPK_BOUND_VECTOR;
  }
  if ((bounds & 3) == 0) {
    box = NULL;
  } else {
    box = opk_new_boxset(vspace,
                         lower_type, lower_value,
                         upper_type, upper_value);
    if (box == NULL) {
      printf("# Failed to create box set\n");
      exit(-1);
    }
  }

  /* Create the optimizer */
  vgp = NULL;
  if (algorithm == VMLMB) {
    opk_vmlmb_options_t options;
    vmlmb = opk_new_vmlmb_optimizer(vspace, mem, vmlmb_flags,
                                    box, lnsrch);
    if (vmlmb == NULL) {
      printf("# Failed to create VMLMB optimizer\n");
      exit(-1);
    }
    algorithm_name = opk_get_vmlmb_method_name(vmlmb);
    if (box != NULL) {
      vgp = opk_vcreate(vspace);
      if (vgp == NULL) {
        printf("# Failed to create projected gradient vector\n");
        exit(-1);
      }
    }
    opk_get_vmlmb_options(&options, vmlmb);
    options.gatol = 0.0;
    options.grtol = 0.0;
    if (delta_given) {
      options.delta = delta;
    } else {
      delta = options.delta;
    }
    if (epsilon_given) {
      options.epsilon = epsilon;
    } else {
      epsilon = options.epsilon;
    }
    if (opk_set_vmlmb_options(vmlmb, &options) != OPK_SUCCESS) {
      printf("# Bad VMLMB options\n");
      exit(-1);
    }
    task = opk_start_vmlmb(vmlmb, vx);
  } else if (algorithm == NLCG) {
    opk_nlcg_options_t options;
    algorithm_name = "NLCG";
    if (bounds != 0) {
      printf("# Algorithm %s cannot be used with bounds\n",
             algorithm_name);
      exit(-1);
    }
    nlcg = opk_new_nlcg_optimizer(vspace, nlcg_flags, lnsrch);
    if (nlcg == NULL) {
      printf("# Failed to create NLCG optimizer\n");
      exit(-1);
    }
    opk_get_nlcg_options(&options, nlcg);
    options.gatol = 0.0;
    options.grtol = 0.0;
    if (delta_given) {
      options.delta = delta;
    } else {
      delta = options.delta;
    }
    if (epsilon_given) {
      options.epsilon = epsilon;
    } else {
      epsilon = options.epsilon;
    }
    if (opk_set_nlcg_options(nlcg, &options) != OPK_SUCCESS) {
      printf("# Bad NLCG options\n");
      exit(-1);
    }
    task = opk_start_nlcg(nlcg, vx);
  } else {
    printf("# Bad algorithm\n");
    exit(-1);
  }

  /* Iteratively call the optimizer */
  evaluations = 0;
  iterations = 0;
  f = -1;
  gnorm = -1;
  finalf = -1;
  finalgnorm = -1;
  gtest = -1;
  if (verbose > 0) {
    char str[256];
    int len1, len2, off;
    strcpy(str, "--------------------------------------------------------------------");
    len1 = strlen(pname);
    len2 = strlen(str);
    off = (len2/2) - (len1/2);
    if (off < 3) {
      off = 3;
    }
    str[off] = '{';
    str[off+1] = ' ';
    strcpy(&str[off+2], pname);
    str[off+len1+2] = ' ';
    str[off+len1+3] = '}';
    printf("# %s\n", str);
    printf("# Iter.  Eval.  Reset            F               ||G||         Step\n");
    printf("# --------------------------------------------------------------------\n");
  }
  final_status = OPK_SUCCESS;
  while (TRUE_) {
    if (task == OPK_TASK_COMPUTE_FG) {
      logical grad = TRUE_;
      if (maxeval > 0 && evaluations > maxeval) {
        final_status = OPK_TOO_MANY_EVALUATIONS;
        task = OPK_TASK_FINAL_X; /* to force terminating */
      } else {
        CUTEST_uofg(&status, &CUTEst_nvar, x, &f, g, &grad);
        if ((status == 1) || (status == 2)) {
          printf("# CUTEst error, status = %d, aborting\n", status);
          exit(status);
        }
        gnorm = gradnorm(vx, vg, vgp, box); /* FIXME: this adds some overhead */
        ++evaluations;
        if (evaluations == 1) {
          gtest = gatol + grtol*gnorm;
        }
        if (gnorm <= gtest) {
          task = OPK_TASK_FINAL_X;
        }
        if (evaluations == 1 || f < finalf) {
          finalf = f;
          finalgnorm = gnorm;
        }
      }
    } else if (task != OPK_TASK_NEW_X) {
      /* Either algorithm converged, or an exception (error or warning)
         occured. */
      if (vmlmb != NULL) {
        final_status = opk_get_vmlmb_status(vmlmb);
      } else {
        final_status = opk_get_nlcg_status(nlcg);
      }
      task = OPK_TASK_FINAL_X; /* to force terminating */
    }
    if (task == OPK_TASK_NEW_X || task == OPK_TASK_FINAL_X) {
      if (evaluations > 1) {
        ++iterations;
      }
      if (verbose > 0 && (task == OPK_TASK_FINAL_X ||
                          (iterations%verbose) == 0)) {
        double alpha;
        opk_index_t restarts;
        if (vmlmb != NULL) {
          alpha = opk_get_vmlmb_step(vmlmb);
          restarts = opk_get_vmlmb_restarts(vmlmb);
        } else {
          alpha = opk_get_nlcg_step(nlcg);
          restarts = 0;
        }
        printf("  %6d %6d %6ld %23.15e %11.3e %11.3e\n",
               iterations, evaluations, (long)restarts, f, gnorm, alpha);
      }
      if (task == OPK_TASK_FINAL_X) {
        break;
      }
      if (maxiter > 0 && iterations >= maxiter) {
        final_status = OPK_TOO_MANY_ITERATIONS;
        break;
      }
    }
    if (vmlmb != NULL) {
      task = opk_iterate_vmlmb(vmlmb, vx, f, vg);
    } else {
      task = opk_iterate_nlcg(nlcg, vx, f, vg);
    }
  }
#if 0
  opk_vprint(stderr, " x", vx, 10);
  opk_vprint(stderr, " g", vg, 10);
  opk_vprint(stderr, "xl", vbl, 10);
  opk_vprint(stderr, "xu", vbu, 10);
#endif

  /* Get CUTEst statistics */
  CUTEST_creport(&status, calls, cpu);
  if (status != 0) {
    printf("# CUTEst error, status = %d, aborting\n", status);
    exit(status);
  }

  /* Print statistics */
  printf("# *********************** CUTEst statistics ************************\n");
  printf("#                Algorithm = OptimPack/%s ", algorithm_name);
  if (algorithm == VMLMB) {
    printf("(mem=%d)\n", OPK_MIN(CUTEst_nvar, mem));
  } else {
    char temp[DESCRIPTION_MAX_SIZE];
    opk_get_nlcg_description(temp, sizeof(temp), nlcg);
    printf("(%s)\n", temp);
  }
  printf("#                  Problem = %-s\n", pname);
  printf("#                Variables = %-10d\n", CUTEst_nvar);
  printf("#        Bound constraints = %-10d\n", vartypes.nbnds);
  printf("#               Iterations = %d\n", iterations);
  printf("# Objective function calls = %.0f\n", calls[0]);
  printf("# Objective gradient calls = %.0f\n", calls[1]);
  printf("#       Objective Hessians = %.0f\n", calls[2]);
  printf("#  Hessian-vector products = %.0f\n", calls[3]);
  printf("#                Exit code = %d (%s)\n", final_status,
         opk_get_reason(final_status));
  printf("#                  Final f = %-23.15e\n", finalf);
  printf("#              Final ||g|| = %-15.7e\n", finalgnorm);
  printf("#              Set up time = %#-14.6f seconds\n", cpu[0]);
  printf("#               Solve time = %#-14.6f seconds\n", cpu[1]);
  printf("#      Relative small step = %g\n", delta);
  printf("#        Descent threshold = %g\n", epsilon);
  printf("# ******************************************************************\n");

  ierr = 0;
  FORTRAN_close(&funit, &ierr);
  if (ierr != 0) {
    printf("# Error closing %s on unit %d.\n", fname, funit);
    printf("# Trying not to abort.\n");
  }

  /* Free workspace */
  FREE(pname);
  OPK_DROP(vx);
  OPK_DROP(vg);
  OPK_DROP(vgp);
  OPK_DROP(vbl);
  OPK_DROP(vbu);
  OPK_DROP(vspace);
  OPK_DROP(nlcg);
  OPK_DROP(vmlmb);
  OPK_DROP(box);
  CUTEST_uterminate(&status);

  return 0;
}

#ifdef __cplusplus
}    /* Closing brace for  extern "C"  block */
#endif
