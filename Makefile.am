AUTOMAKE_OPTIONS = foreign
SUBDIRS = src

# Additional files to include in the distribution.
EXTRA_DIST = doc yorick LICENSE.md README.md

dist-hook:
	find $(distdir)/doc -type f -name "*~" -exec rm -f {} \;
	(cd $(distdir)/yorick; rm -f *~ ywrap.* *.o *.a *.so *.dll *.dylib)

