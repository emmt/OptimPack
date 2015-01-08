# Additional files to include in the distribution.
EXTRA_DIST = README mgh-tests.c mgh-wrappers.c mgh-wrappers.h

# Additional files to clean.
CLEANFILES = *~ nlcg-test.out*

# MinPack tests.
check_PROGRAMS = nlcg-test
nlcg_test_SOURCES = nlcg-test.c mgh-port.c mgh-port.h
nlcg_test_CPPFLAGS = -I$(srcdir)/../src
nlcg_test_LDADD = -L$(srcdir)/../src -lopk -lm
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
