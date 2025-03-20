AUTOMAKE_OPTIONS = foreign
SUBDIRS = src tests
ACLOCAL_AMFLAGS = -I m4

# Additional files to include in the distribution.
EXTRA_DIST = doc yorick .dir-locals.el \
             LICENSE.md README.md doc/NOTES.md doc/INSTALL.md doc/CHANGES.md \
             cutest/install cutest/optimpack cutest/README.md \
             cutest/makemaster cutest/optimpack_main.c \
             cutest/OPTIMPACK.SPC cutest/README.optimpack \
             src/bobyqa-test.out src/cobyla-test.out src/newuoa-test.out \
             src/newuoa-test.out~bigfloat src/newuoa-test.out~old \
             src/newuoa-test.out~orig

dist-hook:
	find "$(distdir)" -type f -name "*~" -exec rm -f {} \;
	rm -rf "$(distdir)/doc/latex"
	rm -rf "$(distdir)/doc/html"
	rm -rf "$(distdir)/yorick/build"
	(cd "$(distdir)/yorick"; rm -f *~ ywrap.* *.o *.a *.so *.dll *.dylib)
