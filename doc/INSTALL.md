## Installation of OptimPack

OptimPack consists in four independent libraries which can be linked together
in a common executable without any conflicts, their functions use specific
prefixes.  The following table summarizes their usage.

| Headers          | Link flags     | Prefix    | Description                          |
|:-----------------|:---------------|:----------|:-------------------------------------|
| `<optimpack*.h>` | `-lopk -lm`    | `opk_`    | OptimPack methods (VMLMB, NLCG, ...) |
| `<bobyqa.h>`     | `-lbobyqa -lm` | `bobyqa_` | Powell's BOBYQA method               |
| `<cobyla.h>`     | `-lcobyla -lm` | `cobyla_` | Powell's COBYLA method               |
| `<newuoa.h>`     | `-lnewuoa -lm` | `newuoa_` | Powell's NEWUOA method               |


### Use pre-compiled libraries

Pre-compiled libraries are available for a variety of platforms and may be
downloaded from [this repository's GitHub releases
page](https://github.com/emmt/OptimPackBuilder/releases).


### Getting the source files

You may download the sources of the last version
([optimpack-3.1.0.tar.gz](https://github.com/emmt/OptimPack/releases/download/v3.1.0/optimpack-3.1.0.tar.gz))
or clone the [OptimPack repository](https://github.com/emmt/OptimPack). In that
latter case, you'll have to generate the `configure` script by executing the
following command in the repository directory:

    ./autogen.sh

Note that this requires that you have installed a recent version of the
*autotools* software (`autoconf`, `automake` and `libtool`).


### Build OptimPack libraries

OptimPack uses standard GNU *autotools* to ease the building and installation
of the the libraries (either shared and static versions are supported).  The
whole process can be as short as the following three steps:

1. Building can be done from the OptimPack distribution directory or from any
   *build* directory.  From the build directory run the `configure` script of
   OptimPack:

        path_to_distribution_directory/configure

   There are many options, so you may run `configure --help` first to figure
   out.

2. Compile the libraries:

        make

3. Install the libraries and header files:

        make install


### Yorick plugin

To install Yorick plugin you must go to the `yorick` sub-directory of the
distribution and follow the instructions [here](../yorick/README.md).
