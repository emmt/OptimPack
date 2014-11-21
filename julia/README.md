# OptimPack.jl

OptimPack.jl is Julia interface for OptimPack library.


## Unconstrained Minimization of a Nonlinear Smooth Function

There are two methods in OptimPack to minimize a nonlinear smooth
multivariate function without constraints: non-linear conjugate gradient
(NLCG) and limited memory variable metric method (VMLM).

The easiest way to use these minimizers is to provide a Julia function, say
`fg!`, which is in charge of computing the value of the function and its
gradient for given variables.  This function must have the form:
```julia
function fg!(x, g)
   g[...] = ... # store the gradient of the function
   f = ...      # compute the function value
   return f     # return the function value
end
```
where the arguments `x` and `g` are Julia arrays (same types and
dimensions) with, on entry, `x` storing the variables and, on exit, `g`
storing the gradient.  The user defined function shall return the function
value.


## Nonlinear Conjugate Gradient (NLCG)

The solution `x` can be computed by one of the implemented nonlinear
conjugate gradient methods with:
```julia
x = OptimPack.nlcg(fg!, x0, method)
```
where `x0` gives the initial value of the variables (as well as the data
type and dimensions of the solution).  `x0` is a Julia dense array with any
dimensions and with elements of type `Float64` or `Float32`.  Argument
`method` is optional and can be used to choose among the different implemented
methods (see below).

The keyword `verb` can be set true to print information at each iteration.
Other keywords are described in the following sub-sections.


### Method Settings

The different nonlinear conjugate gradient methods mainly differ by the way
they compute the search direction.  The conjugate gradient iteration
writes:
```julia
x_{k+1} = x_{k} + alpha_{k} * d_{k}
```
with `alpha_{k}` the step length and where the search direction `d_{k}` is
derived from the gradient `g(x_{k})` of the objective function at the
current point `x_{k}` and from the previous search direction `d_{k-1}` by
an *update rule* which depends on the specific method.  Typically:
```julia
d_{k} = -g(x_{k}) + beta_{k} * d_{k-1}
```
where `beta_{k}` is computed following different recipes.  To choose which
recipe to use, the value of the `method` argument can be set to one of the
following values:

- `OPK_NLCG_FLETCHER_REEVES` for Fletcher & Reeve method;
- `OPK_NLCG_HESTENES_STIEFEL` for Hestenes & Stiefel method;
- `OPK_NLCG_POLAK_RIBIERE_POLYAK` for Polak, Ribière & Polyak method;
- `OPK_NLCG_FLETCHER` for Fletcher "*Conjugate Descent*" method;
- `OPK_NLCG_LIU_STOREY` for Liu & Storey method;
- `OPK_NLCG_DAI_YUAN` for Dai & Yuan method;
- `OPK_NLCG_PERRY_SHANNO` for Perry & Shanno update rule;
- `OPK_NLCG_HAGER_ZHANG` for Hager & Zhang method.

The above values can be bitwise or'ed with the following bits:

- `OPK_NLCG_POWELL` to force parameter `beta` to be nonnegative;
- `OPK_NLCG_SHANNO_PHUA` to guess the stpe length following the
  prescription of Shanno & Phua.

For instance:
```julia
method = OPK_NLCG_POLAK_RIBIERE_POLYAK | OPK_NLCG_POWELL
```
merely corresponds to PRP+ algorithm by Polak, Ribière & Polyak; while:
```julia
method = OPK_NLCG_PERRY_SHANNO | OPK_NLCG_SHANNO_PHUA
```
merely corresponds to the nonlinear conjugate gradient method implemented
in CONMIN.

The default settings for nonlinear conjugate gradient is:
```julia
const OPK_NLCG_DEFAULT  = (OPK_NLCG_HAGER_ZHANG | OPK_NLCG_SHANNO_PHUA)
```


### Convergence Settings

The nonlinear conjugate gradient methods are iterative algorithms, the
convergence is assumed to be achieved when the Euclidean norm of the
gradient is smaller than a threshold.  In pseudo-code, the criterion is:
```julia
||g(x)|| <= max(0, gatol, grtol*||g(x0)||)
```
where `||g(x)||` is the Euclidean norm of the gradient at the current
solution `x`, `||g(x0)||` is the Euclidean norm of the initial gradient at
`x0`, `gatol` is an absolute threshold parameter and `grtol` is a relative
threshold parameter.  The keywords `gatol` and `grtol` can be used to
specify other values for these parameters than the default ones which are
`gatol = 0.0` and `grtol = 1E-6`.


### Linesearch Settings

The keyword `lnsrch` can be used to specify another linesearch method than
the default one:
```julia
x = OptimPack.nlcg(fg!, x0, method, lnsrch=ls)
```
where `ls` is one of the implemented line search methods:
```julia
ls = OptimPack.OptimPackArmijoLineSearch(ftol)
ls = OptimPack.OptimPackMoreThuenteLineSearch(ftol, gtol, xtol)
ls = OptimPack.OptimPackNonmonotoneLineSearch(ftol, m)
```
with `ftol` the tolerance on the function reduction for the Armijo or first
Wolfe condition, `gtol` the tolerance on the gradient for the second
(strong) Wolfe condition, `xtol` the relative precision for the step length
(set to the machine relative precision by default) and `m` the number of
previous steps to remember for the nonmonotone linesearch.

The linesearch is safeguarded by imposing lower and upper bounds on the
step.  Keywords `stpmin` and `stpmax` can be used to specify the step
bounds relatively to the size of the first step (for each linesearch).
Their default values are: `stpmin = 1E-20` and `stpmax = 1E+20`; if
specified, they must be such that: `0 <= stpmin < stpmax`.


## Variable Metric with Limited Memory (VMLM)

Alternatively, the solution `x` can be computed by a limited memory version
of the variable metric method (implementing BFGS updates) with:
```julia
x = OptimPack.vmlm(fg!, x0, m)
```
where the optional argument `m` is the number of previous steps to memorize
(by default `m=3`) while other arguments have the same meaning as for
`OptimPack.nlcg`.

Keywords `verb`, `gatol`, `grtol`, `lnsrch`, `stpmin` and `stpmax` can also
be specified for `OptimPack.vmlm` and have the same meaning as for
`OptimPack.nlcg`.


## Low-level Interface

### Operations on Vectors

To create a vector space for vectors of dimensions `dims` and element type
`T`:
```julia
space = OptimPack.OptimPackShapedVectorSpace(T, dims)
```
where `T` is `Float32` or `Float64` (or any type alias of these,
e.g. `Cfloat` or `Cdouble`) and `dims` is a tuple of the dimensions.

It is also possible to *wrap* a vector around a specific Julia array:
```julia
vect = OptimPack.wrap(space, arr)
```
where `space` is an OptimPack *shaped* vector space and `arr` is a Julia
array.  The element type and dimension list of the array must match those
of the vector space.  A method is available to change the contents of such
a vector:
```julia
OptimPack.wrap!(vect, arr2)
```
where `arr2` is another Julia array (the same constraints on the element
type and dimensions apply).

Note that `arr` must be a **dense array** (of type `DenseArray`) as the
elements of *shaped* vectors are supposed to be stored contiguously in
memory.  OptimPack offers the possibility to create custom vector spaces
and this will be exploited in a near futur to allow for other flavors of
Julia arrays.


### Error Management

Run-time errors throw Julia exception.

