#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# PROJECT AND LIBRARIES VERSIONS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# The project version follows semantic versioning rules (https://semver.org).
#
m4_define([project_version], [3.3.0])
#
# The `-version-info` flag of LIBTOOL takes an argument of the form `current:revision:age`
# where:
#
# - `current` uniquely identifies the latest interface implemented, it can only grow.
#
# - `revision` is the implementation number of `current`, that is the number of bugs
#    fixed. This number shall reset to zero when `current` is incremented.
#
# - `age` is used to specify the backward compatibility. `age + 1` is the number of
#   interfaces implemented by `current`. The library supports all interfaces between
#   `current − age` and `current`.
#
# Here are a set of rules to help you update your library version information:
#
# 1. Start with version information of `0:0:0` for each LIBTOOL library.
#
# 2. Update the version information only immediately before a public release of
#    your software. More frequent updates are unnecessary, and only guarantee
#    that the current interface number gets larger faster.
#
# 3. If the library source code has changed at all since the last update, then
#    increment revision.
#
# 4. If any interfaces have been added, removed, or changed since the last update,
#    increment current, and set revision to 0.
#
# 5. If any interfaces have been added since the last public release, then increment age.
#
# 6. If any interfaces have been removed or changed since the last public release, then
#    set age to 0.
#
m4_define([libopk_version], [3:0:0])
m4_define([libbobyqa_version], [2:1:0])
m4_define([libcobyla_version], [2:1:0])
m4_define([libnewuoa_version], [2:1:0])

AC_PREREQ([2.71])
AC_INIT([OptimPack],[project_version],[eric.thiebaut@univ-lyon1.fr])
AC_CONFIG_SRCDIR([src/optimpack.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign -Wall -Werror dist-bzip2])
AC_SUBST([project_version], [project_version])
AC_SUBST([libopk_version], [libopk_version])
AC_SUBST([libbobyqa_version], [libbobyqa_version])
AC_SUBST([libcobyla_version], [libcobyla_version])
AC_SUBST([libnewuoa_version], [libnewuoa_version])

# Checks for programs.
AM_PROG_AR
AC_PROG_CC

# Turn on all warnings for gcc.  (Options -Wall -Werror in AM_INIT_AUTOMAKE
# have no effects.)
if test "$GCC" = "yes"; then
    #CFLAGS="$CFLAGS -Wall -Werror -Wsign-compare -Wconversion"
    CFLAGS="$CFLAGS -Wall -Werror"
fi

# Use `libtool` for compiling libraries.
LT_INIT([shared])

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([assert.h float.h stdarg.h stddef.h stdlib.h stdio.h string.h errno.h math.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([memset sqrt])

# Output files.
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT
