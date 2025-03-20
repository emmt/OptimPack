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

libopk_la_LDFLAGS = -version-info @libopk_version@
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


libcobyla_la_LDFLAGS = -version-info @libcobyla_version@
libcobyla_la_SOURCES = cobyla.c cobyla.h

libbobyqa_la_LDFLAGS = -version-info @libbobyqa_version@
libbobyqa_la_SOURCES = bobyqa.c bobyqa.h

libnewuoa_la_LDFLAGS = -version-info @libnewuoa_version@
libnewuoa_la_SOURCES = newuoa.c newuoa.h

noinst_PROGRAMS = $(TESTS)

TESTS = \
    bobyqa-test \
    cobyla-test \
    newuoa-test

bobyqa_test_SOURCES = bobyqa-test.c
bobyqa_test__CPPFLAGS = -I$(srcdir)
bobyqa_test_LDADD = -lm

cobyla_test_SOURCES = cobyla-test.c
cobyla_test__CPPFLAGS = -I$(srcdir)
cobyla_test_LDADD = -lm

newuoa_test_SOURCES = newuoa-test.c
newuoa_test__CPPFLAGS = -I$(srcdir)
newuoa_test_LDADD = -lm
