# Specific compilation flags.
AM_CFLAGS = -Wall
#AM_LDFLAGS =

# Additional files to clean.
CLEANFILES = *~

# Headers to install:
include_HEADERS = optimpack.h optimpack-linalg.h bobyqa.h cobyla.h newuoa.h

# Additional files to include in the distribution.
EXTRA_DIST = simple_vector.c

# Target library.
lib_LTLIBRARIES = libopk.la
libopk_la_SOURCES = algebra.c error.c lbfgs.c lcg.c gqtpar.c lnsrch.c linalg.c nlcg.c fmin.c vmlm.c vmlmb.c fltvector.c dblvector.c $(include_HEADERS) optimpack-private.h linalg-common.h linalg-single.h linalg-double.h bobyqa.c bobyqa.h cobyla.c cobyla.h newuoa.c newuoa.h
libopk_la_LDFLAGS = -version-info 0:0:0
# (See https://www.lrde.epita.fr/~adl/autotools.html p. 125-127 for the
# rationale of libtool version numbers.)
