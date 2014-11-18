LIBNAME = libOptimPack2
VERSION_MAJOR = 1
VERSION_MINOR = 0
PATCH_LEVEL = 0

VERSION = $(VERSION_MAJOR).$(VERSION_MINOR).$(PATCH_LEVEL)
CC = gcc

CFLAGS_SHARED = -fpic
CFLAGS_OPTIM = -O2 -DNDEBUG
CFLAGS_DEBUG = -O2 -g

CFLAGS_TEST = -I. -I$(MGH_DIR) -Wall

CFLAGS = -Wall $(CFLAGS_OPTIM) $(CFLAGS_SHARED)
CPPFLAGS = -I.

MGH_DIR = ../tests

HEADERS = optimpack.h optimpack-private.h
OBJS = algebra.o fltvector.o dblvector.o error.o lnsrch.o nlcg.o fmin.o vmlm.o

STATIC_LIBRARY = $(LIBNAME).a
SHARED_LIBRARY = $(LIBNAME).so.$(VERSION)

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin) 
WL_SONAME := install_name
else
WL_SONAME := soname
endif

RM = rm -f

default: $(STATIC_LIBRARY) $(SHARED_LIBRARY)

clean:
	rm -f *~ *.o core $(STATIC_LIBRARY) $(SHARED_LIBRARY)

%.o: %.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $<

$(SHARED_LIBRARY): $(OBJS)
	$(CC) -shared -Wl,-$(WL_SONAME),$(LIBNAME).so.$(VERSION_MAJOR) -o $@ $^ -lm

$(STATIC_LIBRARY): $(OBJS)
	$(AR) rv $@ $^

error.o: error.c optimpack.h
algebra.o: algebra.c $(HEADERS)

fltvector.o: simple_vector.c $(HEADERS)
	$(CC) -c $(CPPFLAGS) $(CFLAGS) -DSINGLE_PRECISION=1 $< -o $@

dblvector.o: simple_vector.c $(HEADERS)
	$(CC) -c $(CPPFLAGS) $(CFLAGS) -DSINGLE_PRECISION=0 $< -o $@

fmin.o: fmin.c $(HEADERS)
lnsrch.o: lnsrch.c $(HEADERS)
nlcg.o: nlcg.c $(HEADERS)
vmlm.o: vmlm.c $(HEADERS)
mgh-port.o: $(MGH_DIR)/mgh-port.c $(MGH_DIR)/mgh-port.h
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@

nlcg-test: nlcg-test.c optimpack.h $(MGH_DIR)/mgh-port.h mgh-port.o $(STATIC_LIBRARY)
	$(CC) $(CFLAGS_TEST) nlcg-test.c mgh-port.o $(STATIC_LIBRARY) -o $@ -lm

run-nlcg-test: nlcg-test
	@for rule in 1 2 3 4 5 6 7 8; do \
	  for opt in 0 256 512; do \
	    meth=`echo $$rule $$opt | awk '{ print $$1 + $$2; }'`; \
	    dest=nlcg-test.out.$$meth; \
	    $(RM) $$dest; \
	    for prob in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18; do \
	      ./nlcg-test $$prob $$meth >> $$dest; \
	    done; \
	  done; \
	done
