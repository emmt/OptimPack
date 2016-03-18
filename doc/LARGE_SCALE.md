# Solving large scale smooth problems

Suppose you want to solve the following problem:
```
min f(x)    s.t.    l <= x <= u             (P1)
```
where `f(x)` is a smooth function of the variables `x`, `l` are lower bounds
and `u` are upper bounds (the inequalities above are componentwise).  If all
values of the lower bound `l` are -∞ and all values of the upper bound `u` are
+∞, then the problem is unconstrained; otherwise, the problem is constrained.
Here a *smooth function* means that it is differentiable.

With a large number of variables `x`, the methods of choice are limited memory
algorithms such as the nonlinear conjugate gradients or quasi-Newton methods.
These methods (with several variants) are implemented in OptimpPack to solve
(P1).  Nonlinear conjugate gradient (NLCG) methods have the smallest memory
footprint but are not as fast as quasi-Newton (though it may depend on the
specific problem and on the settings) and, as implemented in OptimpPack, can
only deal with unconstrained problems.  Limited memory variable metric methods
may use different amount of memory for more efficiency and the version
implemented in OptimpPack can cope with unconstrained and bound constrained
problems.


## Simple driver

OptimPack manages the variables in a way which is agnostic to the actual layout
of the data.  This is a strength of OptimPack but result in some additional
complexity.  A simple driver is implemented in OptimPack to provide limited
memory optimization methods when the variables are *flat* arrays of floating
point values (`float` or `double`) stored contiguously in conventional memory.

This driver uses reverse communication to facilitate its use from any
programming languages: only scalar values and arrays are exchanged (no function
pointer or things like that).

Depending on the settings, optimization can be performed by an instance of
the nonlinear conjugate gradient methods or an instance of the limited
memory quasi-Newton (a.k.a. variable metric) methods.  Simple bound
constraints can be taken into account (providing the variable metric method
is selected).


### Unconstrained minimization example

The following example shows how to solve an unconstrained problem by
a limited memory variable metric method (`OPK_ALGORITHM_VMLMB`):
~~~~~{.cpp}
const int n = 100000;       // size of the problem
const int type = OPK_FLOAT; // type of variables
int mem = 0;                // amount of memory to use, a small default
                            // number is taken if mem <= 0
unigned flags = 0;          // options
double fx;                  // the computed objective function
float x[n];                 // the variables
for (int i = 0; i < n; ++i) x[i] = 0; // initialize variables
float gx[n];                // the gradient of the objective function
opk_optimizer_t* opt = opk_new_optimizer(OPK_ALGORITHM_VMLMB,
                                         type, n, mem, flags,
                                         OPK_BOUND_NONE, NULL,
                                         OPK_BOUND_NONE, NULL,
                                         NULL);
opk_task_t task = opk_start(opt, x);
for (;;) {
    if (task == OPK_TASK_COMPUTE_FG) {
         // compute function value and gradient
         fx = f(x);
         for (i = 0; i < n; ++i) {
             gx[i] = ...;
         }
     } else if (task == OPK_TASK_NEW_X) {
         // a new iterate is available
         fprintf(stdout, "iter=%ld, f(x)=%g, |g(x)|=%g\n",
                 opk_get_iterations(opt), fx,
                 opk_get_gnorm(opt));
     } else {
         break;
     }
     task = opk_iterate(opt, x, fx, gx);
}
if (task != OPK_TASK_FINAL_X) {
    fprintf(stderr, "some error occured (%s)",
            opk_get_reason(opk_get_status(opt)));
}
opk_destroy_optimizer(opt); // finalize the optimizer
~~~~~
It is sufficient to use `OPK_ALGORITHM_NLCG` as the first argument of
`opk_new_optimizer()` to use nonlinear conjugate gradient instead (the value of
`mem` is ignored in this case).


### Bound constrained minimization example

In order to solve a bound constrained problem, you have to provide the upper
and lower bounds to `opk_new_optimizer()`.  Each bound is specified by its type
and by an address.  If the variables are unbounded from above and/or from
below, then the corresponding bound type is `OPK_BOUND_NONE` and the associated
address must be `NULL` (as in the previous example).  If all variables have the
same bound, it is more efficient to specify a scalar bound with type
`OPK_BOUND_SCALAR_FLOAT` or `OPK_BOUND_SCALAR_DOUBLE` along with the address of
a single or double precision variable which stores the bound.  It is also
possible to specify a componentwise bound by providing the address of a
conventional array (with the same number of elements as the variables) and type
`OPK_BOUND_STATIC_FLOAT` or `OPK_BOUND_STATIC_DOUBLE` if the array will not be
released while the bound is in use, `OPK_BOUND_VOLATILE_FLOAT` or
`OPK_BOUND_VOLATILE_DOUBLE` otherwise.

The following example shows how to solve constrained problem bounded below by
zero and above by different upper bounds:
~~~~~{.cpp}
const int n = 100000;       // size of the problem
const int type = OPK_FLOAT; // type of variables
int mem = 0;                // amount of memory to use, a small default
                            // number is taken if mem <= 0
unigned flags = 0;          // options
double lower = 0.0;         // lower bound
float upper[n];             // upper bound
for (int i = 0; i < n; ++i) upper[i] = i; // initialize bounds
double fx;                  // the computed objective function
float x[n];                 // the variables
for (int i = 0; i < n; ++i) x[i] = 0; // initialize variables
float gx[n];                // the gradient of the objective function
opk_optimizer_t* opt = opk_new_optimizer(OPK_ALGORITHM_VMLMB,
                                         type, n, mem, flags,
                                         OPK_BOUND_SCALAR_DOUBLE, &lower,
                                         OPK_BOUND_STATIC_FLOAT, upper,
                                         NULL);
opk_task_t task = opk_start(opt, x);
for (;;) {
    if (task == OPK_TASK_COMPUTE_FG) {
         // compute function value and gradient
         fx = f(x);
         for (i = 0; i < n; ++i) {
             gx[i] = ...;
         }
     } else if (task == OPK_TASK_NEW_X) {
         // a new iterate is available
         fprintf(stdout, "iter=%ld, f(x)=%g, |g(x)|=%g\n",
                 opk_get_iterations(opt), fx,
                 opk_get_gnorm(opt));
     } else {
         break;
     }
     task = opk_iterate(opt, x, fx, gx);
}
if (task != OPK_TASK_FINAL_X) {
    fprintf(stderr, "some error occured (%s)",
            opk_get_reason(opk_get_status(opt)));
}
opk_destroy_optimizer(opt); // finalize the optimizer
~~~~~

*Remarks*:

* The bounds are only specified when calling `opk_new_optimizer()`, the rest of
  the code is unchanged.

* If there are bounds, you must use `OPK_ALGORITHM_VMLMB` (not
  `OPK_ALGORITHM_NLCG`).

* The initial variables

* If the data type of the bounds does not match that of the variable, proper
  conversions are performed but you loose the advantages of having the bounds
  in *static* memory (*i.e.* some memory have to be allocated to store the
  bounds).

