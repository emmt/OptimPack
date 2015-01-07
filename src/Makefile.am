# Specific compilation flags.
#AM_CFLAGS = --pedantic -Wall -std=c99 -O2
#AM_LDFLAGS =

# Headers to install:
include_HEADERS = optimpack.h optimpack-linalg.h

# Additional files to include in the distribution.
EXTRA_DIST = simple_vector.c

# Target library.
lib_LTLIBRARIES = libopk.la
libopk_la_SOURCES = algebra.c error.c lcg.c gqtpar.c lnsrch.c linalg.c nlcg.c fmin.c vmlm.c fltvector.c dblvector.c $(include_HEADERS) optimpack-private.h linalg-common.h linalg-single.h linalg-double.h
