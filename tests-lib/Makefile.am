TESTFILES = tabs.c tabsolute.c tacos.c tacosh.c tadd.c tasin.c tasinh.c tatan.c tatanh.c tbinary.c tbuild_function_abs.c tbuild_function_acos.c tbuild_function_acosh.c tbuild_function_add.c tbuild_function_asin.c tbuild_function_asinh.c tbuild_function_atan.c tbuild_function_atanh.c tbuild_function_ceil.c tbuild_function_cos.c tbuild_function_cosh.c tbuild_function_div.c tbuild_function_double.c tbuild_function_double_double.c tbuild_function_doubleextended.c tbuild_function_erf.c tbuild_function_erfc.c tbuild_function_exp.c tbuild_function_expm1.c tbuild_function_floor.c tbuild_function_free_variable.c tbuild_function_halfprecision.c tbuild_function_log10.c tbuild_function_log1p.c tbuild_function_log2.c tbuild_function_log.c tbuild_function_mul.c tbuild_function_nearestint.c tbuild_function_neg.c tbuild_function_pi.c tbuild_function_pow.c tbuild_function_quad.c tbuild_function_sin.c tbuild_function_single.c tbuild_function_sinh.c tbuild_function_sqrt.c tbuild_function_sub.c tbuild_function_tan.c tbuild_function_tanh.c tbuild_function_triple_double.c tceil.c tcos.c tcosh.c tdecimal.c tdefault.c tdiv.c tdouble.c tdouble_double.c tdouble_double_obj.c tdoubleextended.c tdoubleextended_obj.c tdouble_obj.c tdyadic.c terf.c terfc.c terror.c texp.c texpm1.c tfalse.c tfile.c tfixed.c tfloating.c tfloor.c tfree_variable.c thalfprecision.c thalfprecision_obj.c thexadecimal.c thonorcoeffprec.c tis_absolute.c tis_binary.c tis_decimal.c tis_default.c tis_double_double_obj.c tis_doubleextended_obj.c tis_double_obj.c tis_dyadic.c tis_false.c tis_file.c tis_fixed.c tis_floating.c tis_halfprecision_obj.c tis_hexadecimal.c tis_honorcoeffprec.c tis_off.c tis_on.c tis_perturb.c tis_pi.c tis_postscript.c tis_postscriptfile.c tis_powers.c tis_quad_obj.c tis_relative.c tis_round_down.c tis_round_to_nearest.c tis_round_towards_zero.c tis_round_up.c tis_single_obj.c tis_triple_double_obj.c tis_true.c tis_void.c tlog10.c tlog1p.c tlog2.c tlog.c tmin.c tminus.c tmax.c tmul.c tnearestint.c tobj_is_end_elliptic_list.c tobj_is_error.c tobj_is_function.c tobj_is_list.c tobj_is_range.c tobj_is_string.c toff.c ton.c tperturb.c tpi.c tpostscript.c tpostscriptfile.c tpow.c tpowers.c tquad.c tquad_obj.c trelative.c tround_down.c tround_to_nearest.c tround_towards_zero.c tround_up.c tsin.c tsingle.c tsingle_obj.c tsinh.c tsqrt.c tsub.c ttan.c ttanh.c ttriple_double.c ttriple_double_obj.c ttrue.c tvoid.c

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
