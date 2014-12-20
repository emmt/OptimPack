# OptimPack
# (version 2.0)

***Although it can already be used, this version is a work in progress: not
   all routines have been tested or implemented yet.  The older version of
   OptimPack is available
   [here](http://cral.univ-lyon1.fr/labo/perso/eric.thiebaut/?Software/OptimPack)***

This is OptimPack, a library for solving large scale optimization
problems.  This version implements:

- several non-linear conjugate gradient (NLCG) methods (see refs. [1-3]);

- limited memory variable metric (LBFGS, see ref. [4]) possibly with bound
  constraints and/or preconditioning (VMLMB, see ref. [5], or BLMVM, see ref. [6]);

- spectral project gradient (SPG, seee ref. [7]) method;

- inexact monotone and nonmonotone line searches (see ref. [7,8]);

- computaion of a trust region step;

- linear conjugate gradients [1] and trust region conjugate gradient [9].

Most of the documention is in the header files, *e.g.*
[src/optimpack.h](src/optimpack.h), in Doxygen format.


## OptimPack Bindings

OptimPack library is written in C but, in order to make embedding
OptimPack into another language as easy as possible, the routines use
reverse communication: all local variables needed by the optimizers
get saved into workspaces created by the library and the optimizers
never explicitely call the penalty function to optimize.

The following language bindings allow OptimPack to be used in other
programming languages:

* Directory
  [yorick](https://github.com/emmt/OptimPack/tree/master/yorick)
  contains an implementation of OptimPack support in Yorick
  (https://github.com/emmt/OptimPack/tree/master/yorick).

* The [OptimPack.jl](https://github.com/emmt/OptimPack.jl) project
  implements of OptimPack support for [Julia](http://julialang.org/).

* Directory `idl` contains an implementation of OptimPack support in
  [IDL](http://en.wikipedia.org/wiki/IDL_%28programming_language%29)/[PV-Wave](http://www.roguewave.com/products-services/pv-wave)/[GDL](http://gnudatalanguage.sourceforge.net/)
  (using `CALL_EXTERNAL`).

* The [TiPi](https://github.com/emmt/TiPi) project provides a
  framework for solving inverse image reconstruction problems in
  **Java**.  The optimization package of TiPi is a Java version of the
  OptimPack library.


## References

1. M.R. Hestenes & E. Stiefel, "*Methods of Conjugate Gradients for
   Solving Linear Systems*," Journal of Research of the National Bureau
   of Standards 49, 409-436 (1952).

2. W.W. Hager & H. Zhang, "*A survey of nonlinear conjugate gradient
   methods*," Pacific Journal of Optimization, Vol. 2, pp. 35-58
   (2006).

3. W.W. Hager & H. Zhang, "*A New Conjugate Gradient Method with
   Guaranteed Descent and an Efficient Line Search*," SIAM J. Optim.,
   Vol. 16, pp. 170-192 (2005).

4. D. Liu and J. Nocedal, "*On the limited memory BFGS method for
   large scale optimization*", Mathematical Programming B **45**,
   503-528 (1989).

5. É. Thiébaut, "*Optimization issues in blind deconvolution
   algorithms*," in Astronomical Data Analysis II, SPIE
   Proc. **4847**, 174-183 (2002).

6. S.J. Benson & J.J. Moré, "*A limited memory variable metric method
   in subspaces and bound constrained optimization problems*", in
   Subspaces and Bound Constrained Optimization Problems, (2001).

7. E.G. Birgin, J.M. Martínez & M. Raydan, "*Nonmonotone Spectral
   Projected Gradient Methods on Convex Sets*," SIAM J. Optim. **10**,
   1196-1211 (2000).

8. Jorge J. Moré and David J. Thuente, "*Line search algorithms with
   guaranteed sufficient decrease*" in ACM Transactions on Mathematical
   Software (TOMS) Volume 20, Issue 3, Pages 286-307 (September 1994).

9. T. Steihaug, "*The conjugate gradient method and trust regions in large
   scale optimization*", SIAM Journal on Numerical Analysis, vol. **20**,
   pp. 626-637, 1983.


## Installation

1. Go to directory `src` and edit the `Makefile`.

2. Compile the library:

        make


## Authors

* Éric Thiébaut (https://github.com/emmt)


## License

The OptimPack library is released under under the
[MIT "Expat" License](LICENSE.md).
