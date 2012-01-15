TESTFILES= check00.c check01.c

EXTRA_DIST= $(TESTFILES)

check_PROGRAMS = $(TESTFILES:.c=)
TESTS= $(check_PROGRAMS)

LDADD = $(top_builddir)/.libs/libsollya.la

.PHONY: $(top_builddir)/sollya

$(top_builddir)/sollya:
	@echo "Recompiling sollya";cd $(top_builddir); $(MAKE) all

$(TESTS): $(top_builddir)/sollya