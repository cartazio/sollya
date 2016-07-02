SUBDIRS = doc tests-tool tests-lib
EXTRA_DIST = libsollyaexample.c externalfun.c externalproc.c sollya.h2m

AM_YFLAGS = -d
AM_LFLAGS = -o$(LEX_OUTPUT_ROOT).c

# The following line is added to follow a suggestion from libtoolize.
ACLOCAL_AMFLAGS = -I m4

BUILT_SOURCES = \
miniparser.h \
parser.h \
internparser.h \
lexer.c \
minilexer.c \
internlexer.c

bin_PROGRAMS = sollya

sollya_SOURCES = \
$(libsollya_la_SOURCES) \
main.c \
main.h

sollya_CFLAGS = $(AM_CFLAGS)
sollya_CXXFLAGS = $(AM_CXXFLAGS)

include_HEADERS = sollya.h sollya-messages.h

lib_LTLIBRARIES = libsollya.la
libsollya_la_SOURCES = \
sollya.h \
mpfi-compat.h \
mpfi-compat.c \
hash.c \
hash.h \
hooks.c \
hooks.h \
polynomials.c \
polynomials.h \
sollya-library-wrappers.c \
sollya-library-wrappers.h \
sollya-messaging.h \
sollya-messages.h \
sollya-messaging.c \
bitfields.h \
bitfields.c \
printf.c \
printf.h \
internparser.y \
miniparser.y \
parser.y \
internlexer.l \
lexer.l \
minilexer.l \
help.h \
version.h \
assignment.c \
assignment.h \
autodiff.c \
autodiff.h \
chain.c \
chain.h \
double.c \
double.h \
execute.c \
execute.h \
expression.c \
expression.h \
external.c \
external.h \
general.c \
general.h \
signalhandling.c \
signalhandling.h \
implement.c \
implement.h \
implementconst.c \
implementconst.h \
infnorm.c \
infnorm.h \
integral.c \
integral.h \
library.c \
library.h \
plot.c \
plot.h \
proof.c \
proof.h \
fpminimax.cpp \
fpminimax.h \
remez.c \
remez.h \
match.c \
match.h \
taylorform.c \
taylorform.h \
chebyshevform.c \
chebyshevform.h \
chebyshevformaux.c \
chebyshevformaux.h \
supnorm.c \
supnorm.h \
sturm.c \
sturm.h \
taylor.c \
taylor.h \
worstcase.c \
worstcase.h \
xml.c \
xml.h

libsollya_la_LDFLAGS = \
	-version-info 5:0:0

##LIBTOOL_DEPS = @LIBTOOL_DEPS@
##libtool: $(LIBTOOL_DEPS)
##	$(SHELL) ./config.status --recheck

check-lib: $(BUILT_SOURCES)
	@fail= failcom='exit 1'; \
	for f in x $$MAKEFLAGS; do \
	  case $$f in \
	    *=* | --[!k]*);; \
	    *k*) failcom='fail=yes';; \
	  esac; \
	done; \
	echo "Making check in tests-lib"; \
	($(am__cd) tests-lib && $(MAKE) $(AM_MAKEFLAGS) check) \
	|| eval $$failcom; \
	$(MAKE) $(AM_MAKEFLAGS) check-am || exit 1; \
	test -z "$$fail"

check-tool: $(BUILT_SOURCES)
	@fail= failcom='exit 1'; \
	for f in x $$MAKEFLAGS; do \
	  case $$f in \
	    *=* | --[!k]*);; \
	    *k*) failcom='fail=yes';; \
	  esac; \
	done; \
	echo "Making check in tests-tool"; \
	($(am__cd) tests-tool && $(MAKE) $(AM_MAKEFLAGS) check) \
	|| eval $$failcom; \
	$(MAKE) $(AM_MAKEFLAGS) check-am || exit 1; \
	test -z "$$fail"

man_MANS = sollya.1

AM_H2MFLAGS = \
	--local=C.utf8 \
	--libtool \
	--no-info

# Short descriptions for manpages
sollya_DESCRIPTION =  "tool environment for safe floting-point code development -- command-line interface"

%.1: %
	 SOLLYA_BANNER_OFF=YES $(HELP2MAN) \
			-s 1 \
			$(AM_H2MFLAGS) \
			-I $(top_srcdir)/$*.h2m \
			-n $(if $($(subst -,_,$*)_DESCRIPTION), $($(subst -,_,$*)_DESCRIPTION), "manual page for $*") \
			-o $@ \
		$(top_builddir)/$<

#EXTRA_DIST += $(man_MANS)
