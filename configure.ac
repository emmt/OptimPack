#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([OptimPack], [2.0.0a], [eric.thiebaut@univ-lyon1.fr])
AC_CONFIG_SRCDIR([src/optimpack.h])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB
AM_PROG_AR

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
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
