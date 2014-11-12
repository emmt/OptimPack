CC = gcc

CFLAGS_OPTIM = -O2 -DNDEBUG
CFLAGS_DEBUG = -O2 -g

CFLAGS_TEST = -I. -I$(MGH_DIR) -Wall

CFLAGS = -Wall $(CFLAGS_OPTIM)
CPPFLAGS = -I.

MGH_DIR = ../tests

HEADERS = optimpack.h optimpack-private.h
OBJS = algebra.o fltvector.o dblvector.o error.o lnsrch.o nlcg.o fmin.o
LIBNAME = libOptimPack2.a

RM = rm -f

default: $(LIBNAME)

clean:
	rm -f *~ *.o core

%.o: %.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $<

$(LIBNAME): $(OBJS)
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
vmlm.o: lnsrch.c $(HEADERS)
mgh-port.o: $(MGH_DIR)/mgh-port.c $(MGH_DIR)/mgh-port.h
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@

nlcg-test: nlcg-test.c optimpack.h $(MGH_DIR)/mgh-port.h mgh-port.o $(LIBNAME)
	$(CC) $(CFLAGS_TEST) nlcg-test.c mgh-port.o $(LIBNAME) -o $@ -lm

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
