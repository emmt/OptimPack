# Developer notes for OptimPack

## Publishing a new release

To prepare a new public release, follow these steps:

1. If the API has changed, update the interface version number in macros `libopk_version`,
   `libcobyla_version`, `libbobyqa_version`, and `libnewuoa_version` `libopk_la_LDFLAGS`
   in [configure.ac](../configure.ac). See [libtool
   documentation](http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html)
   for explanations. The interface version number is `c:r:a` (for *current*, *revision*
   and *age*), revisions between `c` and `c-a` are supposed to be supported by the
   library. If the library source code has changed at all since the last update, then
   increment revision (`c:r:a` becomes `c:r+1:a`). If any interfaces have been added,
   removed, or changed since the last update, increment current and set revision to 0
   (`c:r:a` becomes `c+1:0:0`). If any interfaces have been added since the last public
   release, then increment age. If any interfaces have been removed or changed since the
   last public release, then set age to 0. **Never** try to set the interface numbers so
   that they correspond to the release number of the package.

2. Bump release numbers in [README.md](../README.md) and in the macro `project_version` in
   file [configure.ac](../configure.ac).

3. Commit the changes made in steps 1 and 2.

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

6. Open [OptimPack](https://github.com/emmt/OptimPack) page at GitHub and click on the
   [releases](https://github.com/emmt/OptimPack/releases) tab and then on the [Draft a new
   release](https://github.com/emmt/OptimPack/releases/new) button.

   - Tag the version.
   - Add a brief description and a detailed description of the most important changes.
   - Check/uncheck **This is a pre-release**.
   - Attach the release archive.
   - Click on the green **Publish release** button.

7. In the local git repository, pull the tag with:
   ```{.sh}
   git pull --all
   ```

To have binaries automatically build by [Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil):

1. Fork the Yggdrasil project on GitHub:

2. Clone the fork, make a branch, and edit the file `Yggdrasil/O/OptimPack/build_tarballs.jl`:
   ```{.sh}
   git clone git@github.com:emmt/Yggdrasil.git
   cd Yggdrasil
   git co -b optimpack
   vi O/OptimPack/build_tarballs.jl
   ```
   In this file, at least update version number and SHA256 code for the source tarball:
   ```{.sh}
   sha256sum optimpack-${VERSION}.tar.gz
   ```

3. Check that products can be built (from the `Yggdrasil/O/OptimPack` directory):
   ```{.sh}
   julia-1.7 build_tarballs.jl
   ```
   You may upload the products (in the directory `Yggdrasil/O/OptimPack/products`) to
   the release page on GitHub.

4. Commit changes, push the branch to GitHub:
   ```{.sh}
   git checkin -m "Update for OptimPack version ${VERSION}" build_tarballs.jl
   git push --all
   ```

5. As suggested by this last command, create a pull request at
   https://github.com/emmt/Yggdrasil/pull/new/optimpack.


## Prerequisites

Developers under Debian/Ubuntu need the following packages: `build-essential`, `gcc` or
`clang`, `autoconf`, `automake`, `autotools-dev` and `libtool`.

## Initial step size/scaling

The `k`-th iteration of a conjugate gradient method consists in updating the variables as:

    x[k+1] = x[k] + alpha[k]*d[k]

where `alpha[k] > 0` is the step length and `d[k]` is the search direction given by:

    d[k] = -g[k] + beta[k]*d[k-1]

where `g[k]` is the gradient of the objective function at `x[k]`. Different formulae have
been proposed to compute `beta[k]`.

A line search method is used to find a suitable `alpha` along the search path. The initial
value of `alpha` can be:

    alpha[k] = alpha[k-1]*<d[k-1],g[k-1]>/<d[k],g[k]>

which is the rule used by Shanno & Phua (1978) in their CONMIN algorithm.

In quasi-Newton methods, the initial approximation of the inverse Hessian is often given
by a simple scaling by a parameter `gamma > 0`. This yields to the following search
direction:

    d[k] = -gamma[k]*g[k]

Possible rules to guess the value `gamma` from the previous iterations are:

    gamma[k] = <s[k-1],y[k-1]>/<y[k-1],y[k-1]>    (Oren & Spedicato)

    alpha[k] = <s[k-1],s[k-1]>/<s[k-1],y[k-1]>    (Barzilai & Borwein)

where:

    s[k] = x[k+1] - x[k]
    y[k] = g[k+1] - g[k]

are the difference between successive variables and gradients respectively.
