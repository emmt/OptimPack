# Specific compilation flags.
AM_CFLAGS = -Wall
#AM_LDFLAGS =

# Additional files to clean.
CLEANFILES = *~

# Headers to install:
OPTIMPACK_INCS = optimpack.h optimpack-linalg.h optimpack-private.h
include_HEADERS = $(OPTIMPACK_INCS) bobyqa.h cobyla.h newuoa.h

# Additional files to include in the distribution.
EXTRA_DIST = simple_vector.c

# Target library.
lib_LTLIBRARIES = libopk.la libcobyla.la libbobyqa.la libnewuoa.la
libopk_la_SOURCES = \
    algebra.c \
    utils.c \
    lcg.c \
    gqtpar.c \
    lnsrch.c \
    linalg.c \
    nlcg.c \
    fmin.c \
    vmlmb.c \
    fltvector.c \
    dblvector.c \
    driver.c \
    $(OPTIMPACK_INCS) \
    linalg-common.h \
    linalg-single.h \
    linalg-double.h
libopk_la_LDFLAGS = -version-info 3:0:0
libcobyla_la_SOURCES = optimpack.h cobyla.c cobyla.h
libcobyla_la_LDFLAGS = -version-info 2:0:0
libbobyqa_la_SOURCES = optimpack.h bobyqa.c bobyqa.h
libbobyqa_la_LDFLAGS = -version-info 2:0:0
libnewuoa_la_SOURCES = optimpack.h newuoa.c newuoa.h
libnewuoa_la_LDFLAGS = -version-info 2:0:0
# (See doc/NOTES.md for the rationale of libtool version numbers.)
