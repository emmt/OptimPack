#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.69])
AC_INIT([OptimPack], [3.2.0], [eric.thiebaut@univ-lyon1.fr])
AC_CONFIG_SRCDIR([src/optimpack.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign -Wall -Werror dist-bzip2])

# Checks for programs.
AM_PROG_AR
AC_PROG_CC_C99

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
