# OptimPack.jl

OptimPack.jl is Julia interface for OptimPack library.


## Unconstrained Minimization of a Nonlinear Smooth Function

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
```julia
x = OptimPack.nlcg(fg!, x0, method)
```
where `x0` gives the initial value of the variables (as well as the data
type and dimensions of the solution).  `x0` is a Julia dense array with any
dimensions and with elements of type `Float64` or `Float32`.  Argument
`method` is optional and can be used to achoose among the different methods
implemented.

Alternatively, the solution `x` can be computed by a limited memory version
of the variable metric method (with BFGS updates) with:
```julia
x = OptimPack.vmlm(fg!, x0, m)
```
where the optional argument `m` is the number of previous steps to memorize
(by default `m=3`) while other arguments have the same meaning as for
`OptimPack.nlcg`.


## Low-level Interface

### Operations on Vectors

To create a vector space for vectors of dimensions `dims` and element type `T`:
```julia
space = OptimPack.OptimPackShapedVectorSpace(T, dims)
```
where `T` is `Float32` or `Float64` (or any type alias of these,
e.g. `Cfloat` or `Cdouble`) and `dims` is a tuple of the dimensions.

It is also possible to *wrap* a vector around a specific Julia arrays:
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

