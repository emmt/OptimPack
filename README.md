# OptimPack
# (version 2.0)

![Travis build status](http://travis-ci.org/emmt/OptimPack.png)

***Although it can already be used, this version is a work in progress: not
   all routines have been tested or implemented yet.  The older version of
   OptimPack is available
   [here](http://cral.univ-lyon1.fr/labo/perso/eric.thiebaut/?Software/OptimPack)***

This is **OptimPack**, a library for solving optimization problems.  For
problems of small to moderate size, **OptimPack** provides:

- computation of a trust region step;

- Mike Powell's **COBYLA** (see ref. [10]), **NEWUOA** (see ref. [11]), and
  **BOBYQA** (see ref. [11]) algorithms for minimizing a function of many
  variables.  These methods are *derivatives free* (only the function
  values are needed).  **COBYLA** accounts for general inequality
  constraints.  **BOBYQA** accounts for bound constraints on the variables.

For large scale problems involving millions of variables (or more),
**OptimPack** provides:

- several non-linear conjugate gradient (NLCG) methods (see refs. [1-3]);

- limited memory variable metric (LBFGS, see ref. [4]) possibly with bound
  constraints and/or preconditioning (VMLMB, see ref. [5], or BLMVM, see
  ref. [6]);

- spectral project gradient (SPG, see ref. [7]) method;

- inexact monotone and nonmonotone line searches (see ref. [7,8]);

- linear conjugate gradients [1] and trust region conjugate gradient [9].

Most of the documentation is in the header files, *e.g.*
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

10. M.J.D. Powell, "*A direct search optimization method that models the
    objective and constraint functions by linear interpolation*," in
    Advances in Optimization and Numerical Analysis Mathematics and Its
    Applications, vol. **275** (eds. Susana Gomez and Jean-Pierre Hennart),
    Kluwer Academic Publishers, pp. 51-67 (1994).

11. M.J.D. Powell, "*The NEWUOA software for unconstrained minimization
    without derivatives*", in Large-Scale Nonlinear Optimization, editors
    G. Di Pillo and M. Roma, Springer (2006), pages 255-297.

12. M.J.D. Powell, "*The BOBYQA Algorithm for Bound Constrained
    Optimization Without Derivatives*."  Technical report, Department of
    Applied Mathematics and Theoretical Physics, University of Cambridge
    (2009).


## Installation

### OptimPack library

OptimPack uses standard GNU *autotools* to ease the building and
installation of the the library (either shared and static versions are
supported).  The whole process can be as short as the following three
steps:

1. Building can be done from the OptimPack distribution directory or from
   any *build* directory.  From the build directory run the `configure`
   script of OptimPack:

        path_to_distribution_directory/configure

   There are many options, so you may run `configure --help` first to
   figure out.

   If you directly clone OptimPack from GitHub, you'll have to generate the
   `configure` script before:

        ./autogen.sh

   Note that this requires that you have installed a recent version of the
   *autotools* (`autoconf`, `automake` and `libtool`).

2. Compile the library:

        make

3. Install the library:

        make install


### Yorick plugin

To install Yorick plugin you must go to the `yorick` subdirectory of the
distribution and follow instrutions there.


## Authors

* Éric Thiébaut (https://github.com/emmt)


## Credits

The development of OptimPack was supported by the
[MiTiV](http://mitiv-univ-lyon1.fr) project funded by the French *Agence
Nationale pour la Recherche* (ref. ANR-09-EMER-008).


## License

The OptimPack library is released under under the
[MIT "Expat" License](LICENSE.md).
