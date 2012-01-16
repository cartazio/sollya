TESTFILES = tbuild_function_abs.c tbuild_function_acos.c tbuild_function_acosh.c tbuild_function_add.c tbuild_function_asin.c tbuild_function_asinh.c tbuild_function_atan.c tbuild_function_atanh.c tbuild_function_ceil.c tbuild_function_cos.c tbuild_function_cosh.c tbuild_function_div.c tbuild_function_double.c tbuild_function_double_double.c tbuild_function_doubleextended.c tbuild_function_erf.c tbuild_function_erfc.c tbuild_function_exp.c tbuild_function_expm1.c tbuild_function_floor.c tbuild_function_free_variable.c tbuild_function_halfprecision.c tbuild_function_log10.c tbuild_function_log1p.c tbuild_function_log2.c tbuild_function_log.c tbuild_function_mul.c tbuild_function_nearestint.c tbuild_function_neg.c tbuild_function_pi.c tbuild_function_pow.c tbuild_function_quad.c tbuild_function_sin.c tbuild_function_single.c tbuild_function_sinh.c tbuild_function_sqrt.c tbuild_function_sub.c tbuild_function_tan.c tbuild_function_tanh.c tbuild_function_triple_double.c

EXTRA_DIST = check.proto $(TESTFILES) $(TESTFILES:.c=.reference)

check_PROGRAMS = force_compilation $(TESTFILES:.c=)
TESTS= $(TESTFILES:.c=.sh)

LDADD = $(top_builddir)/.libs/libsollya.la

.PHONY: force_compilation

force_compilation:
	@echo "Recompiling sollya";cd $(top_builddir); $(MAKE) all

$(TESTS): check.proto
	@echo "Building "$@
	@echo "#!/bin/sh" > $@
	@echo "TESTFILES=\""$(@:.sh=)"\"" >> $@
	@cat check.proto >> $@
	@chmod +x $@
