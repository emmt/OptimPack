## Installation of OptimPack

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
distribution and follow the instrutions [there](./yorick/README.md).
