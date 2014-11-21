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


## Non-Linear Conjugate Gradient (NLCG)

The solution `x` can be computed by one of the implemented non-linear
conjugate gradient methods with:
```julia
x = OptimPack.nlcg(fg!, x0, method)
```
where `x0` gives the initial value of the variables (as well as the data
type and dimensions of the solution).  `x0` is a Julia dense array with any
dimensions and with elements of type `Float64` or `Float32`.  Argument
`method` is optional and can be used to choose among the different implemented
methods.

The non-linear conjugate gradient is an iterative algorithm, the
convergence is assumed when the Euclidean norm of the gradient is
smaller than a threshold:
```julia
||g|| <= max(0, gatol, grtol*||g0||)
```
where `||g||` is the Euclidean norm of the current gradient, `||g0||`
is Euclidean norm of the current gradient, `gatol` is an absolute
threshold parameter and `grtol` is a relative threshold parameter.
The keywords `gatol` and `grtol` can be used to specify other values
for these parameters than the default ones which are `gatol = 0.0` and
`grtol = 1E-6`.

The keyword `verb` can be set true to print information at each
iteration.

The keyword `lnsrch` can be used to specify another linesearch method
than the default one:
```julia
x = OptimPack.nlcg(fg!, x0, method, lnsrch=ls)
```
where `ls` is one of the implemented line search methods:
```julia
ls = OptimPack.OptimPackArmijoLineSearch(ftol)
ls = OptimPack.OptimPackMoreThuenteLineSearch(ftol, gtol, xtol)
ls = OptimPack.OptimPackNonmonotoneLineSearch(ftol, m)
```
with `ftol` the tolerance on the function reduction for the Armijo or
first Wolfe condition, `gtol` the tolerance on the gradient the second
(strong) Wolfe condition, `xtol` the relative precision for the step
length (set to the machine relative precision by default) and `m` the
number of previous steps to remember for the nonmonotone linesearch.

The linesearch is safeguarded by imposing a lower and an upper bounds
on the step.  Keywords `stpmin` and `stpmax` can be used to specify
the step bounds relatively to the size of the first step (for each
linesearch).  Their default values are: `stpmin = 1E-20` and `stpmax =
1E+20`; they must be given such that: `0 <= stpmin < stpmax`.


## Limited Memory Variable Metric (VMLM)

Alternatively, the solution `x` can be computed by a limited memory version
of the variable metric method (implementing BFGS updates) with:
```julia
x = OptimPack.vmlm(fg!, x0, m)
```
where the optional argument `m` is the number of previous steps to memorize
(by default `m=3`) while other arguments have the same meaning as for
`OptimPack.nlcg`.

Keywords `verb`, `gatol`, `grtol`, `lnsrch`, `stpmin` and `stpmax` can
also be specified for `OptimPack.vmlm` and have the same meaning as
for `OptimPack.nlcg`.


## Low-level Interface

### Operations on Vectors

To create a vector space for vectors of dimensions `dims` and element type `T`:
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

