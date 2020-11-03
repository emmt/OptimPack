[![Build Status](https://travis-ci.com/emmt/OptimPack.svg?branch=master)](https://travis-ci.com/emmt/OptimPack)

# OptimPack (version 3.1.0)

This is **OptimPack**, a C library for solving optimization problems.  The
library is mostly targeted at very large problems (*e.g.* as the ones
encountered in image restoration) but also provide routines for problems of
smaller size.

This document provides a general overview of **OptimPack**, for more
specific information, see:

- [Software installation.](./doc/INSTALL.md)
- [Developer notes.](./doc/NOTES.md)
- [Software changes and history.](./doc/CHANGES.md)
- [Solving large scale smooth problems.](./doc/LARGE_SCALE.md)

Most of the documentation is in the header files, *e.g.*
[src/optimpack.h](src/optimpack.h), in Doxygen format.


## Large scale problems

For large scale problems involving millions of variables (or even more),
**OptimPack** provides:

- several non-linear conjugate gradient (NLCG) methods (see refs. [1-3]);

- limited memory variable metric (LBFGS, see Ref. [4]) possibly with bound
  constraints and/or preconditioning (VMLMB, see Ref. [5], or BLMVM, see
  Ref. [6]);

- inexact monotone and nonmonotone line searches (see Ref. [7,8]);

- linear conjugate gradients [1].

The large scale optimizers of the **OptimPack** library can work with the
unknowns stored in almost any form (provided a minimal set of functions to
manipulate them are implemented).  This feature may be used to exploit hardware
acceleration, multi-threading or to distribute the storage and computation
across multiple machines.

See [Solving large scale smooth problems](./doc/LARGE_SCALE.md) for examples
and more information about using Optimpack to solve large scale problems.


## Problems of small to moderate size

For problems of small to moderate size, **OptimPack** provides:

- Moré & Sorensen method to compute a trust region step (see Ref. [13]);

- Mike Powell's **COBYLA** (see Ref. [10]), **NEWUOA** (see Ref. [11]), and
  **BOBYQA** (see Ref. [12]) algorithms for minimizing a function of many
  variables.  These methods are *derivatives free* (only the function values
  are needed).  **NEWUOA** is for unconstrained optimization.  **COBYLA**
  accounts for general inequality constraints.  **BOBYQA** accounts for bound
  constraints on the variables.

- Brent's method for the minimization of an univariate function (see
  Ref. [14]).


## OptimPack bindings

OptimPack library is written in C but, in order to make embedding OptimPack
into another language as easy as possible, most routines use **reverse
communication**: all local variables needed by the optimizers get saved into
workspaces created by the library and the optimizers never directly call the
function to optimize.

The following language bindings allow OptimPack to be used in other programming
languages:

* Directory [yorick](https://github.com/emmt/OptimPack/tree/master/yorick)
  contains an implementation of OptimPack support in Yorick
  (https://github.com/emmt/OptimPack/tree/master/yorick).

* The [OptimPack.jl](https://github.com/emmt/OptimPack.jl) project implements
  of OptimPack support for [Julia](http://julialang.org/).

* [OptimPackNextGen.jl](https://github.com/emmt/OptimPackNextGen.jl) is a new
  project to provide the same features as OptimPack for
  [Julia](http://julialang.org/) but with most code written in pure Julia.  For
  now, pure Julia version of the methods devoted to large scale problems are
  available.  To use Powell algorithms (devoted to moderate size problems) the
  dynamic libraries of OptimPack are still needed.

* The [TiPi](https://github.com/emmt/TiPi) project provides a framework for
  solving inverse image reconstruction problems in **Java**.  The optimization
  package of TiPi is a Java version of the OptimPack library.


## References

1. M.R. Hestenes & E. Stiefel, "*Methods of Conjugate Gradients for Solving
   Linear Systems*," Journal of Research of the National Bureau of Standards
   49, pp. 409-436 (1952).

2. W.W. Hager & H. Zhang, "*A survey of nonlinear conjugate gradient methods*,"
   Pacific Journal of Optimization, Vol. 2, pp. 35-58 (2006).

3. W.W. Hager & H. Zhang, "*A New Conjugate Gradient Method with Guaranteed
   Descent and an Efficient Line Search*," SIAM J. Optim., Vol. 16, pp. 170-192
   (2005).

4. D. Liu and J. Nocedal, "*On the limited memory BFGS method for large scale
   optimization*", Mathematical Programming B **45**, 503-528 (1989).

5. É. Thiébaut, "*Optimization issues in blind deconvolution algorithms*," in
   Astronomical Data Analysis II, SPIE Proc. **4847**, 174-183 (2002).

6. S.J. Benson & J.J. Moré, "*A limited memory variable metric method in
   subspaces and bound constrained optimization problems*", in Subspaces and
   Bound Constrained Optimization Problems, (2001).

7. E.G. Birgin, J.M. Martínez & M. Raydan, "*Nonmonotone Spectral Projected
   Gradient Methods on Convex Sets*," SIAM J. Optim. **10**, 1196-1211 (2000).

8. Jorge J. Moré and David J. Thuente, "*Line search algorithms with guaranteed
   sufficient decrease*" in ACM Transactions on Mathematical Software (TOMS)
   Volume 20, Issue 3, Pages 286-307 (1994).

9. T. Steihaug, "*The conjugate gradient method and trust regions in large
   scale optimization*", SIAM Journal on Numerical Analysis, vol. **20**,
   pp. 626-637 (1983).

10. M.J.D. Powell, "*A direct search optimization method that models the
    objective and constraint functions by linear interpolation*," in
    Advances in Optimization and Numerical Analysis Mathematics and Its
    Applications, vol. **275** (eds. Susana Gomez and Jean-Pierre Hennart),
    Kluwer Academic Publishers, pp. 51-67 (1994).

11. M.J.D. Powell, "*The NEWUOA software for unconstrained minimization without
    derivatives*", in Large-Scale Nonlinear Optimization, editors G. Di Pillo
    and M. Roma, Springer, pp. 255-297 (2006).

12. M.J.D. Powell, "*The BOBYQA Algorithm for Bound Constrained Optimization
    Without Derivatives*."  Technical report, Department of Applied Mathematics
    and Theoretical Physics, University of Cambridge (2009).

13. J.J. Moré & D.C. Sorensen, "*Computing A Trust Region Step*," SIAM
    J. Sci. Stat. Comp. **4**, 553-572 (1983).

14. R.P. Brent, "*Algorithms for Minimization without Derivatives*,"
    Prentice-Hall, Inc. (1973).


## Authors

* Éric Thiébaut (https://github.com/emmt)


## Credits

The development of OptimPack was supported by the
[MiTiV](http://mitiv-univ-lyon1.fr) project funded by the French
[*Agence Nationale pour la Recherche*](http://www.agence-nationale-recherche.fr)
(Ref. ANR-09-EMER-008).

A simpler version of OptimPack is provided by
[OptimPackLegacy](https://github.com/emmt/OptimPackLegacy) project.


## License

The OptimPack library is released under the
[MIT "Expat" License](LICENSE.md).
