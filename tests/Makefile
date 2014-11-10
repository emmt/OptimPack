

MGH_SRCDIR = ../orig/mgh
MGH_OBJS = grdfcn.o lhesfcn.o neqipt.o ssqfcn.o umipt.o vecjac.o \
           hesfcn.o lsqipt.o objfcn.o ssqjac.o vecfcn.o

RM = rm -f

CC = gcc
CFLAGS = -Wall -O2
CPPFLAGS = -I.

FC = gfortran
FFLAGS = $(CFLAGS)

LD = $(CC)
LDFLAGS =
LDLIBS = -lm

FILTER_NUMBERS = sed -e 's/\([.0-9]\)[eE]\([-+]\)/\1D\2/g'

OBJS = lnsrch.o nlcg.o

default: $(OBJS)

# Default rules to compile C code.
%.o: %.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $<

# Default rules to compile FORTRAN code.
%.o: %.f
	$(FC) -c $(FFLAGS) $< -o $@

default: check

clean:
	$(RM) *~ *.o core
	$(RM) mgh-tests

check: conmin-test.out conmin-test.ref
	colordiff conmin-test.ref conmin-test.out

#------------------------------------------------------------------------------
# Routines from MINPACK-1 project.

mgh-wrappers.o: mgh-wrappers.c mgh-wrappers.h
	$(CC) -c -I. $(CPPFLAGS) $(CFLAGS) $< -o $@

%.o: $(MGH_SRCDIR)/%.f
	$(FC) -c $(FFLAGS) $< -o $@

mgh-objs: $(MGH_OBJS)
umipt.o: $(MGH_SRCDIR)/umipt.f
grdfcn.o: $(MGH_SRCDIR)/grdfcn.f
objfcn.o: $(MGH_SRCDIR)/objfcn.f
lhesfcn.o: $(MGH_SRCDIR)/lhesfcn.f
hesfcn.o: $(MGH_SRCDIR)/hesfcn.f
neqipt.o: $(MGH_SRCDIR)/neqipt.f
lsqipt.o: $(MGH_SRCDIR)/lsqipt.f
ssqfcn.o: $(MGH_SRCDIR)/ssqfcn.f
ssqjac.o: $(MGH_SRCDIR)/ssqjac.f
vecfcn.o: $(MGH_SRCDIR)/vecfcn.f
vecjac.o: $(MGH_SRCDIR)/vecjac.f

MGH_TESTS_OBJS = umipt.o grdfcn.o objfcn.o \
                 hesfcn.o ssqfcn.o ssqjac.o vecfcn.o vecjac.o \
                 mgh-wrappers.o mgh-port.o mgh-tests.o

mgh-tests.o: mgh-tests.c mgh-wrappers.h mgh-port.h

mgh-tests: $(MGH_TESTS_OBJS)
	$(LD) $(LDFLAGS) $(MGH_TESTS_OBJS) -o $@ $(LDLIBS)

run-mgh-tests: mgh-tests
	@echo "*** Running MGH tests (comparison between FORTRAN and C codes)"
	@./mgh-tests | grep -v 'rmax = 0, amax = 0'

#------------------------------------------------------------------------------

conmin-test-orig: woods.f conmin.f
	$(FC) $(FFLAGS) woods.f conmin.f -o $@ -lm

minpack1-tests.o: $(MGH)/minpack1-tests.c $(MGH)/minpack1-tests.c
	$(CC) -c -I$(MGH) $(CPPFLAGS) $(CFLAGS) $< -o $@

conmin-test: conmin.c $(MGH)
	$(CC) -I$(MGH) $(CPPFLAGS) -DTESTING -DMGH $(CFLAGS) \
	  minpack1-tests.o conmin.c -o $@ -lm

conmin-test-orig.out: conmin-test-orig
	./conmin-test-orig > $@

conmin-test.out: conmin-test
	./conmin-test > $@
#	./conmin-test | $(FILTER_NUMBERS) > $@

#------------------------------------------------------------------------------
