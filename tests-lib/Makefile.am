TESTFILES= check00.c check01.c

EXTRA_DIST= $(TESTFILES)

check_PROGRAMS = $(TESTFILES:.c=)
TESTS= $(check_PROGRAMS)

LDADD = $(top_builddir)/.libs/libsollya.la

