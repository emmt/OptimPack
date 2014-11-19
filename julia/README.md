OptimPack.jl is Julia interface for OptimPack library.


Unconstrained Minimization of a Nonlinear Smooth Function
=========================================================

There are two methods in OptimPack which can be used to minimize a
nonlinear smooth multivariate function without constraints: non-linear
conjugate gradient (NLCG) and limited memory variable metric method (VMLM).

The easiest way to use these minimizers is to provide a Julia function, say
`fg!`, which is in charge of computing the value of the function and its
gradient for given variables.  This function must have the form:
```
function fg!(x, g)
   g[...] = ... # store the gradient of the function
   f = ...      # compute the function value
   return f     # return the function value
end
```
where the arguments `x` and `g` are Julia arrays (same types and
dimensions) respectively with (on entry) the variables and (on exit) the
gradient.  The user defined function shall return the function value.

The solution `x` can be computed by one of the implemented non-linear
conjugate gradient methods with:
```
x = OptimPack.nlcg(fg!, x0, method)
```
where `x0` gives the initial value of the variables (as well as the data
type and dimensions of the solution).  `x0` is a Julia dense array with any
dimensions and with elements of type `Float64` or `Float32`.  Argument
`method` is optional and can be used to achoose among the different methods
implemented.

Alternatively, the solution `x` can be computed by a limited memory version
of the variable metric method (with BFGS updates) with:
```
x = OptimPack.vmlm(fg!, x0, m)
```
where the optional argument `m` is the number of previous steps to memorize
(by deafult `m=3`) while other arguments have the same meaning as for
`OptimPack.nlcg`.


Operations on Vectors
=====================

To create a vector space for vectors of length `len` and element type `T`:
```
s = OptimPackShapedVectorSpace(T, dims)
```
where `T` is `Float32` or `Float64` (or any type alias of these,
e.g. `Cfloat` or `Cdouble`) and `dims` gives the list of dimensions.

It is also possible to create a vector space to wrap specific Julia arrays.


Error Management
================

Run-time errors throw Julia exception.

