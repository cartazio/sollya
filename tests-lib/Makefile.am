TESTFILES = check00.c check01.c

EXTRA_DIST = $(TESTFILES)

check_PROGRAMS = force_compilation $(TESTFILES:.c=)
TESTS= $(TESTFILES:.c=)

LDADD = $(top_builddir)/.libs/libsollya.la

.PHONY: force_compilation

force_compilation:
	@echo "Recompiling sollya";cd $(top_builddir); $(MAKE) all

