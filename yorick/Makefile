# where are the sources? (automatically filled in by configure script)
srcdir=.

# these values filled in by "yorick -batch make.i" or configure script
Y_MAKEDIR=
Y_EXE=
Y_EXE_PKGS=
Y_EXE_HOME=
Y_EXE_SITE=
Y_HOME_PKG=

# ----------------------------------------------------- optimization flags

# options for make command line, e.g.-   make COPT=-g TGT=exe
COPT=$(COPT_DEFAULT)
TGT=$(DEFAULT_TGT)

# ------------------------------------------------ macros for this package

PKG_NAME=opky
PKG_I=${srcdir}/opky.i

OBJS=opky.o algebra.o fltvector.o dblvector.o utils.o gqtpar.o linalg.o \
     nlcg.o lnsrch.o vmlmb.o cobyla.o ycobyla.o newuoa.o ynewuoa.o

# change to give the executable a name other than yorick
PKG_EXENAME=yorick

PREFIX=/usr

# PKG_DEPLIBS=-Lsomedir -lsomelib   for dependencies of this package
PKG_DEPLIBS=
# set compiler (or rarely loader) flags specific to this package
PKG_CFLAGS=
PKG_LDFLAGS=

# list of additional package names you want in PKG_EXENAME
# (typically $(Y_EXE_PKGS) should be first here)
EXTRA_PKGS=$(Y_EXE_PKGS)

# list of additional files for clean
PKG_CLEAN=

# autoload file for this package, if any
PKG_I_START=${srcdir}/opky-start.i
# non-pkg.i include files for this package, if any
PKG_I_EXTRA=

# -------------------------------- standard targets and rules (in Makepkg)

# set macros Makepkg uses in target and dependency names
# DLL_TARGETS, LIB_TARGETS, EXE_TARGETS
# are any additional targets (defined below) prerequisite to
# the plugin library, archive library, and executable, respectively
PKG_I_DEPS=$(PKG_I)
Y_DISTMAKE=distmake

ifeq (,$(strip $(Y_MAKEDIR)))
$(info *** WARNING: Y_MAKEDIR not defined, you may run 'yorick -batch make.i' first)
else
include $(Y_MAKEDIR)/Make.cfg
include $(Y_MAKEDIR)/Makepkg
include $(Y_MAKEDIR)/Make$(TGT)
endif

# override macros Makepkg sets for rules and other macros
# see comments in Y_HOME/Makepkg for a list of possibilities

# if this package built with mpy: 1. be sure mpy appears in EXTRA_PKGS,
# 2. set TGT=exe, and 3. uncomment following two lines
# Y_MAIN_O=$(Y_LIBEXE)/mpymain.o
# include $(Y_MAKEDIR)/Makempy

# configure script for this package may produce make macros:
# include output-makefile-from-package-configure

# reduce chance of yorick-1.5 corrupting this Makefile
MAKE_TEMPLATE = protect-against-1.5

# ------------------------------------- targets and rules for this package

# Dummy default target in case Y_MAKEDIR was not defined:
dummy-default:
	@echo >&2 "*** ERROR: Y_MAKEDIR not defined, aborting..."; false

%.o: ${srcdir}/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

# Settings for OptimPack:
OPK_SRC = $(srcdir)/../src
OPK_CFLAGS = $(CPPFLAGS) $(CFLAGS) -DNDEBUG -I$(OPK_SRC)
OPK_HEADERS =  $(OPK_SRC)/optimpack.h $(OPK_SRC)/optimpack-private.h
OPK_LINALG_HEADERS = $(OPK_HEADERS) $(OPK_SRC)/optimpack-linalg.h \
    $(OPK_SRC)/linalg-single.h $(OPK_SRC)/linalg-double.h \
    $(OPK_SRC)/linalg-common.h

opky.o: $(srcdir)/opky.c $(OPK_SRC)/optimpack.h
	$(CC) $(OPK_CFLAGS) -o $@ -c $<

algebra.o: $(OPK_SRC)/algebra.c $(OPK_HEADERS)
	$(CC) $(OPK_CFLAGS) -o $@ -c $<

fltvector.o: $(OPK_SRC)/simple_vector.c $(OPK_HEADERS)
	$(CC) $(OPK_CFLAGS) -DSINGLE_PRECISION=1 -o $@ -c $<

dblvector.o: $(OPK_SRC)/simple_vector.c $(OPK_HEADERS)
	$(CC) $(OPK_CFLAGS) -DSINGLE_PRECISION=0 -o $@ -c $<

utils.o: $(OPK_SRC)/utils.c $(OPK_HEADERS)
	$(CC) $(OPK_CFLAGS) -o $@ -c $<

gqtpar.o:$(OPK_SRC)/gqtpar.c $(OPK_LINALG_HEADERS)
	$(CC) $(OPK_CFLAGS) -o $@ -c $<

linalg.o:$(OPK_SRC)/linalg.c $(OPK_LINALG_HEADERS)
	$(CC) $(OPK_CFLAGS) -o $@ -c $<

nlcg.o: $(OPK_SRC)/nlcg.c $(OPK_HEADERS)
	$(CC) $(OPK_CFLAGS) -o $@ -c $<

lnsrch.o: $(OPK_SRC)/lnsrch.c $(OPK_HEADERS)
	$(CC) $(OPK_CFLAGS) -o $@ -c $<

vmlmb.o: $(OPK_SRC)/vmlmb.c $(OPK_HEADERS)
	$(CC) $(OPK_CFLAGS) -o $@ -c $<

cobyla.o: $(OPK_SRC)/cobyla.c $(OPK_SRC)/optimpack.h $(OPK_SRC)/cobyla.h
	$(CC) $(OPK_CFLAGS) -o $@ -c $<

ycobyla.o: $(srcdir)/ycobyla.c $(OPK_SRC)/optimpack.h $(OPK_SRC)/cobyla.h
	$(CC) $(OPK_CFLAGS) -o $@ -c $<

newuoa.o: $(OPK_SRC)/newuoa.c $(OPK_SRC)/optimpack.h $(OPK_SRC)/newuoa.h
	$(CC) $(OPK_CFLAGS) -o $@ -c $<

ynewuoa.o: $(srcdir)/ynewuoa.c $(OPK_SRC)/optimpack.h $(OPK_SRC)/newuoa.h
	$(CC) $(OPK_CFLAGS) -o $@ -c $<


.PHONY: clean release

# -------------------------------------------------------- end of Makefile
