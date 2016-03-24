# Developer notes for OptimPack

## Publishing a new release

To prepare a new public release, follow these steps:

1. If the API has changed, update the interface version number in macro
   `libopk_la_LDFLAGS` in [src/Makefile.am](../src/Makefile.am) see
   [libtool documentation](http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html)
   for explanations.  The interface version number is `c:r:a` (for *current*,
   *revision* and *age*), revisions between `c` and `c-a` are supposed to be
   supported by the library.  If the library source code has changed at all
   since the last update, then increment revision (`c:r:a` becomes `c:r+1:a`).
   If any interfaces have been added, removed, or changed since the last
   update, increment current, and set revision to 0.  If any interfaces have
   been added since the last public release, then increment age.  If any
   interfaces have been removed or changed since the last public release, then
   set age to 0.  **Never** try to set the interface numbers so that they
   correspond to the release number of the package.

2. Bump version number in [README.md](../README.md) and
   [configure.ac](../configure.ac) and commit the changes.

3. Update the [CHANGES.md](./CHANGES.md) file and commit the changes.

4. Rebuild [configure](./configure) script and make the archive:
   ```{.sh}
   cd build
   make
   make dist-gzip
   mv optimpack-${VERSION}.tar.gz ../releases/.
   ```

5. Push the changes on GitHub:
   ```{.sh}
   git push
   ```

6. Open [OptimPack](https://github.com/emmt/OptimPack) page at GitHub and click
   on the [releases](https://github.com/emmt/OptimPack/releases) tab and then
   on the [Draft a new release](https://github.com/emmt/OptimPack/releases/new)
   button.

   - Tag the version.
   - Add a brief description and a detailed description of the most important
     changes.
   - Check/uncheck **This is a pre-release**.
   - Attach the release archive.
   - Click on the green **Publish release** button.

7. In the local git repository, pull the tag with:
   ```{.sh}
   git pull --all
   ```


## Initial step size/scaling

The `k`-th iteration of a conjugate gradient method consists in updating the
variables as:

    x[k+1] = x[k] + alpha[k]*d[k]

where `alpha[k] > 0` is the step length and `d[k]` is the search direction
given by:

    d[k] = -g[k] + beta[k]*d[k-1]

where `g[k]` is the gradient of the objective function at `x[k]`.  Different
formulae have been proposed to compute `beta[k]`.

A line search method is used to find a suitable `alpha` along the search
path.  The initial value of `alpha` can be:

    alpha[k] = alpha[k-1]*<d[k-1],g[k-1]>/<d[k],g[k]>

which is the rule used by Shanno & Phua (1978) in their CONMIN algorithm.


In quasi-Newton methods, the initial approximation of the inverse Hessian is
often given by a simple scaling by a parameter `gamma > 0`.  This yields to
the following search direction:

    d[k] = -gamma[k]*g[k]

Possible rules to guess the value `gamma` from the previous iterations are:

    gamma[k] = <s[k-1],y[k-1]>/<y[k-1],y[k-1]>    (Oren & Spedicato)

    alpha[k] = <s[k-1],s[k-1]>/<s[k-1],y[k-1]>    (Barzilai & Borwein)

where:

    s[k] = x[k+1] - x[k]
    y[k] = g[k+1] - g[k]

are the difference between successive variables and gradients respectively.
