/*

Copyright 2011 by

Laboratoire d'Informatique de Paris 6, equipe PEQUAN,
UPMC Universite Paris 06 - CNRS - UMR 7606 - LIP6, Paris, France,

Contributors Ch. Lauter

christoph.lauter@ens-lyon.org

This software is a computer program whose purpose is to provide an
environment for safe floating-point code development. It is
particularily targeted to the automatized implementation of
mathematical floating-point libraries (libm). Amongst other features,
it offers a certified infinity norm, an automatic polynomial
implementer and a fast Remez algorithm.

This software is governed by the CeCILL-C license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-C
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
herefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-C license and that you accept its terms.

This program is distributed WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

*/

#include <stdarg.h>
#include <gmp.h>
#include <mpfr.h>
#include <stdint.h>
#include <stdio.h>
#include "expression.h"
#include "execute.h"
#include "chain.h"
#include "mpfi-compat.h"
#include "sollya-library-wrappers.h"


int sollya_lib_init() {
  return 0; 
}

int sollya_lib_close() {
  return 0; 
}

int sollya_lib_install_msg_callback(int (*callback_func) (int)) {
  return 0; 
}

int sollya_lib_uninstall_msg_callback() {
  return 0; 
}

int sollya_lib_printf(const char *format, ...) {
  va_list varlist;
  int res;

  va_start(varlist,format);

  res = sollyaVfprintf(stdout,format,varlist);

  va_end(varlist);

  return res;
}

int sollya_lib_fprintf(FILE *fd, const char *format, ...) {
  va_list varlist;
  int res;

  va_start(varlist,format);

  res = sollyaVfprintf(fd,format,varlist);

  va_end(varlist);

  return res;
}

void sollya_lib_obj_clear(sollya_obj_t obj1) {
  freeThing(obj1);
}

int sollya_lib_cmp_objs_structurally(sollya_obj_t obj1, sollya_obj_t obj2) {
  return isEqualThing(obj1, obj2); 
}

void sollya_lib_plot(sollya_obj_t obj1, sollya_obj_t obj2, ...) {
 
}

void sollya_lib_printdouble(sollya_obj_t obj1) {
 
}

void sollya_lib_printsingle(sollya_obj_t obj1) {
 
}

void sollya_lib_printexpansion(sollya_obj_t obj1) {
 
}

void sollya_lib_implementconst(sollya_obj_t obj1, ...) {
 
}

void sollya_lib_bashexecute(sollya_obj_t obj1) {
 
}

void sollya_lib_externalplot(sollya_obj_t obj1, sollya_obj_t obj2, sollya_obj_t obj3, sollya_obj_t obj4, sollya_obj_t obj5, ...) {
 
}

void sollya_lib_asciiplot(sollya_obj_t obj1, sollya_obj_t obj2) {
 
}

void sollya_lib_printxml(sollya_obj_t obj1) {
 
}

void sollya_lib_printxml_newfile(sollya_obj_t obj1, sollya_obj_t obj2) {
 
}

void sollya_lib_printxml_appendfile(sollya_obj_t obj1, sollya_obj_t obj3) {
 
}

void sollya_lib_worstcase(sollya_obj_t obj1, sollya_obj_t obj2, sollya_obj_t obj3, sollya_obj_t obj4, sollya_obj_t obj5, ...) {
 
}

void sollya_lib_autoprint(sollya_obj_t obj1, ...) {
 
}

void sollya_lib_set_prec(sollya_obj_t obj1) {
 
}

void sollya_lib_set_points(sollya_obj_t obj1) {
 
}

void sollya_lib_set_diam(sollya_obj_t obj1) {
 
}

void sollya_lib_set_display(sollya_obj_t obj1) {
 
}

void sollya_lib_set_verbosity(sollya_obj_t obj1) {
 
}

void sollya_lib_set_canonical(sollya_obj_t obj1) {
 
}

void sollya_lib_set_autosimplify(sollya_obj_t obj1) {
 
}

void sollya_lib_set_taylorrecursions(sollya_obj_t obj1) {
 
}

void sollya_lib_set_timing(sollya_obj_t obj1) {
 
}

void sollya_lib_set_midpointmode(sollya_obj_t obj1) {
 
}

void sollya_lib_set_dieonerrormode(sollya_obj_t obj1) {
 
}

void sollya_lib_set_rationalmode(sollya_obj_t obj1) {
 
}

void sollya_lib_set_roundingwarnings(sollya_obj_t obj1) {
 
}

void sollya_lib_set_hopitalrecursions(sollya_obj_t obj1) {
 
}

void sollya_lib_set_prec_silent(sollya_obj_t obj1) {
 
}

void sollya_lib_set_points_silent(sollya_obj_t obj1) {
 
}

void sollya_lib_set_diam_silent(sollya_obj_t obj1) {
 
}

void sollya_lib_set_display_silent(sollya_obj_t obj1) {
 
}

void sollya_lib_set_verbosity_silent(sollya_obj_t obj1) {
 
}

void sollya_lib_set_canonical_silent(sollya_obj_t obj1) {
 
}

void sollya_lib_set_autosimplify_silent(sollya_obj_t obj1) {
 
}

void sollya_lib_set_taylorrecursions_silent(sollya_obj_t obj1) {
 
}

void sollya_lib_set_timing_silent(sollya_obj_t obj1) {
 
}

void sollya_lib_set_midpointmode_silent(sollya_obj_t obj1) {
 
}

void sollya_lib_set_dieonerrormode_silent(sollya_obj_t obj1) {
 
}

void sollya_lib_set_rationalmode_silent(sollya_obj_t obj1) {
 
}

void sollya_lib_set_roundingwarnings_silent(sollya_obj_t obj1) {
 
}

void sollya_lib_set_hopitalrecursions_silent(sollya_obj_t obj1) {
 
}

sollya_obj_t sollya_lib_free_variable() {
  return makeVariable(); 
}

sollya_obj_t sollya_lib_and(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_or(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_negate(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_cmp_equal(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_cmp_in(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_cmp_less(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_cmp_greater(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_cmp_less_equal(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_cmp_greater_equal(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_cmp_not_equal(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_add(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_sub(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_concat(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_append(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_prepend(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_apply(sollya_obj_t obj1, sollya_obj_t obj2, ...) {
  return NULL; 
}

sollya_obj_t sollya_lib_approx(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_mul(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_div(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_pow(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_minus(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_sup(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_mid(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_inf(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_diff(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_simplify(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_bashevaluate(sollya_obj_t obj1, ...) {
  return NULL; 
}

sollya_obj_t sollya_lib_remez(sollya_obj_t obj1, sollya_obj_t obj2, sollya_obj_t obj3, ...) {
  return NULL; 
}

sollya_obj_t sollya_lib_min(sollya_obj_t obj1, ...) {
  return NULL; 
}

sollya_obj_t sollya_lib_max(sollya_obj_t obj1, ...) {
  return NULL; 
}

sollya_obj_t sollya_lib_fpminimax(sollya_obj_t obj1, sollya_obj_t obj2, sollya_obj_t obj3, sollya_obj_t obj4, ...) {
  return NULL; 
}

sollya_obj_t sollya_lib_horner(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_canonical(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_expand(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_simplifysafe(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_taylor(sollya_obj_t obj1, sollya_obj_t obj2, sollya_obj_t obj3) {
  return NULL; 
}

sollya_obj_t sollya_lib_taylorform(sollya_obj_t obj1, sollya_obj_t obj2, sollya_obj_t obj3, ...) {
  return NULL; 
}

sollya_obj_t sollya_lib_autodiff(sollya_obj_t obj1, sollya_obj_t obj2, sollya_obj_t obj3) {
  return NULL; 
}

sollya_obj_t sollya_lib_degree(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_numerator(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_denominator(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_substitute(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_composepolynomials(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_coeff(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_subpoly(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_roundcoefficients(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_rationalapprox(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_round(sollya_obj_t obj1, sollya_obj_t obj2, sollya_obj_t obj3) {
  return NULL; 
}

sollya_obj_t sollya_lib_evaluate(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_parse(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_readxml(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_infnorm(sollya_obj_t obj1, sollya_obj_t obj2, ...) {
  return NULL; 
}

sollya_obj_t sollya_lib_supnorm(sollya_obj_t obj1, sollya_obj_t obj2, sollya_obj_t obj3, sollya_obj_t obj4, sollya_obj_t obj5) {
  return NULL; 
}

sollya_obj_t sollya_lib_findzeros(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_dirtyinfnorm(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_numberroots(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_integral(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_dirtyintegral(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_implementpoly(sollya_obj_t obj1, sollya_obj_t obj2, sollya_obj_t obj3, sollya_obj_t obj4, sollya_obj_t obj5, sollya_obj_t obj6, ...) {
  return NULL; 
}

sollya_obj_t sollya_lib_checkinfnorm(sollya_obj_t obj1, sollya_obj_t obj2, sollya_obj_t obj3) {
  return NULL; 
}

sollya_obj_t sollya_lib_zerodenominators(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_searchgal(sollya_obj_t obj1, ...) {
  return NULL; 
}

sollya_obj_t sollya_lib_guessdegree(sollya_obj_t obj1, sollya_obj_t obj2, sollya_obj_t obj3, ...) {
  return NULL; 
}

sollya_obj_t sollya_lib_dirtyfindzeros(sollya_obj_t obj1, sollya_obj_t obj2, sollya_obj_t obj3) {
  return NULL; 
}

sollya_obj_t sollya_lib_head(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_roundcorrectly(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_revert(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_sort(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_mantissa(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_exponent(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_tail(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_range(sollya_obj_t obj1, sollya_obj_t obj2) {
  return NULL; 
}

sollya_obj_t sollya_lib_sqrt(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_exp(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_log(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_log2(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_log10(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_sin(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_cos(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_tan(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_asin(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_acos(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_atan(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_sinh(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_cosh(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_tanh(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_asinh(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_acosh(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_atanh(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_abs(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_erf(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_erfc(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_log1p(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_expm1(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_double(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_single(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_quad(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_halfprecision(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_double_double(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_triple_double(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_doubleextended(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_ceil(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_floor(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_nearestint(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_length(sollya_obj_t obj1) {
  return NULL; 
}

sollya_obj_t sollya_lib_get_prec() {
  return NULL; 
}

sollya_obj_t sollya_lib_get_points() {
  return NULL; 
}

sollya_obj_t sollya_lib_get_diam() {
  return NULL; 
}

sollya_obj_t sollya_lib_get_display() {
  return NULL; 
}

sollya_obj_t sollya_lib_get_verbosity() {
  return NULL; 
}

sollya_obj_t sollya_lib_get_canonical() {
  return NULL; 
}

sollya_obj_t sollya_lib_get_autosimplify() {
  return NULL; 
}

sollya_obj_t sollya_lib_get_taylorrecursions() {
  return NULL; 
}

sollya_obj_t sollya_lib_get_timing() {
  return NULL; 
}

sollya_obj_t sollya_lib_get_midpointmode() {
  return NULL; 
}

sollya_obj_t sollya_lib_get_dieonerrormode() {
  return NULL; 
}

sollya_obj_t sollya_lib_get_rationalmode() {
  return NULL; 
}

sollya_obj_t sollya_lib_get_roundingwarnings() {
  return NULL; 
}

sollya_obj_t sollya_lib_get_hopitalrecursions() {
  return NULL; 
}

sollya_obj_t sollya_lib_on() {
  return makeOn();
}

sollya_obj_t sollya_lib_off() {
  return makeOff(); 
}

sollya_obj_t sollya_lib_dyadic() {
  return makeDyadic(); 
}

sollya_obj_t sollya_lib_powers() {
  return makePowers(); 
}

sollya_obj_t sollya_lib_binary() {
  return makeBinaryThing(); 
}

sollya_obj_t sollya_lib_hexadecimal() {
  return makeHexadecimalThing(); 
}

sollya_obj_t sollya_lib_file() {
  return makeFile(); 
}

sollya_obj_t sollya_lib_postscript() {
  return makePostscript(); 
}

sollya_obj_t sollya_lib_postscriptfile() {
  return makePostscriptFile(); 
}

sollya_obj_t sollya_lib_perturb() {
  return makePerturb(); 
}

sollya_obj_t sollya_lib_round_down() {
  return makeRoundDown(); 
}

sollya_obj_t sollya_lib_round_up() {
  return makeRoundUp(); 
}

sollya_obj_t sollya_lib_round_towards_zero() {
  return makeRoundToZero(); 
}

sollya_obj_t sollya_lib_round_to_nearest() {
  return makeRoundToNearest(); 
}

sollya_obj_t sollya_lib_honorcoeffprec() {
  return makeHonorCoeff(); 
}

sollya_obj_t sollya_lib_true() {
  return makeTrue(); 
}

sollya_obj_t sollya_lib_false() {
  return makeFalse(); 
}

sollya_obj_t sollya_lib_void() {
  return makeUnit(); 
}

sollya_obj_t sollya_lib_default() {
  return makeDefault(); 
}

sollya_obj_t sollya_lib_decimal() {
  return makeDecimal(); 
}

sollya_obj_t sollya_lib_absolute() {
  return makeAbsolute(); 
}

sollya_obj_t sollya_lib_relative() {
  return makeRelative(); 
}

sollya_obj_t sollya_lib_fixed() {
  return makeFixed(); 
}

sollya_obj_t sollya_lib_floating() {
  return makeFloating(); 
}

sollya_obj_t sollya_lib_error() {
  return makeError(); 
}

sollya_obj_t sollya_lib_double_obj() {
  return makeDoubleSymbol(); 
}

sollya_obj_t sollya_lib_single_obj() {
  return makeSingleSymbol(); 
}

sollya_obj_t sollya_lib_quad_obj() {
  return makeQuadSymbol(); 
}

sollya_obj_t sollya_lib_halfprecision_obj() {
  return makeHalfPrecisionSymbol(); 
}

sollya_obj_t sollya_lib_doubleextended_obj() {
  return makeDoubleextendedSymbol(); 
}

sollya_obj_t sollya_lib_double_double_obj() {
  return makeDoubleDoubleSymbol(); 
}

sollya_obj_t sollya_lib_triple_double_obj() {
  return makeTripleDoubleSymbol(); 
}

sollya_obj_t sollya_lib_pi() {
  return makePi(); 
}

sollya_obj_t sollya_lib_parse_string(char *str) {
  return parseString(str);
}

sollya_obj_t sollya_lib_string(char *str) {
  return makeString(str);
}

sollya_obj_t sollya_lib_range_from_interval(sollya_mpfi_t interval) {
  sollya_obj_t temp;
  mp_prec_t prec;
  mpfr_t left, right;

  prec = sollya_mpfi_get_prec(interval);
  mpfr_init2(left,prec);
  mpfr_init2(right,prec);
  sollya_mpfi_get_left(left,interval);
  sollya_mpfi_get_right(right,interval);
  temp = makeRange(makeConstant(left),makeConstant(right));
  mpfr_clear(left);
  mpfr_clear(right);

  return temp;
}

sollya_obj_t sollya_lib_range_from_bounds(mpfr_t left, mpfr_t right) {
  return makeRange(makeConstant(left),makeConstant(right)); 
}

sollya_obj_t sollya_lib_constant(mpfr_t value) {
  return makeConstant(value); 
}

sollya_obj_t sollya_lib_constant_from_double(double value) {
  sollya_obj_t temp;
  mpfr_t valueAsMpfr;

  mpfr_init2(valueAsMpfr,64); /* On some systems, doubles actually are double-extended */
  mpfr_set_d(valueAsMpfr, value, GMP_RNDN);
  temp = makeConstant(valueAsMpfr);
  mpfr_clear(valueAsMpfr);

  return temp; 
}

sollya_obj_t sollya_lib_constant_from_int(int value) {
  sollya_obj_t temp;
  mpfr_t valueAsMpfr;

  mpfr_init2(valueAsMpfr,8 * sizeof(int) + 5); 
  mpfr_set_si(valueAsMpfr, value, GMP_RNDN);
  temp = makeConstant(valueAsMpfr);
  mpfr_clear(valueAsMpfr);

  return temp; 
}

sollya_obj_t sollya_lib_constant_from_int64(int64_t value) {
  sollya_obj_t temp;
  double valueDoubleHi, valueDoubleLo;
  int64_t valueHi, valueLo, tempInt64;
  mpfr_t valueMpfr, tempMpfr;

  /* Cut the 64 bit signed value into two pieces
     such that 

     value = valueHi * 2^32 + valueLo.

     Then convert both values to double, which will always 
     hold. Adjust the high value. Finally convert to MPFR.
  */
  valueHi = value;
  valueHi >>= 32;       /* Performs floor, even signed */
  tempInt64 = valueHi;
  tempInt64 <<= 32;     /* Exact */
  valueLo = value - tempInt64;

  valueDoubleHi = (double) valueHi;
  valueDoubleLo = (double) valueLo;
  valueDoubleHi *= 4294967296.0;   /* Multiply by 2^32 */

  mpfr_init2(valueMpfr, 64);
  mpfr_init2(tempMpfr, 64);
  mpfr_set_d(valueMpfr, valueDoubleHi, GMP_RNDN);
  mpfr_set_d(tempMpfr, valueDoubleLo, GMP_RNDN);
  mpfr_add(valueMpfr, valueMpfr, tempMpfr, GMP_RNDN); /* exact */
  temp = makeConstant(valueMpfr);
  mpfr_clear(valueMpfr);
  mpfr_clear(valueMpfr);

  return temp; 
}

sollya_obj_t sollya_lib_constant_from_uint64(uint64_t value) {
  sollya_obj_t temp;
  double valueDoubleHi, valueDoubleLo;
  uint64_t valueHi, valueLo, tempInt64;
  mpfr_t valueMpfr, tempMpfr;

  /* Cut the 64 bit unsigned value into two pieces
     such that 

     value = valueHi * 2^32 + valueLo.

     Then convert both values to double, which will always 
     hold. Adjust the high value. Finally convert to MPFR.
  */
  valueHi = value;
  valueHi >>= 32;       /* Performs floor */
  tempInt64 = valueHi;
  tempInt64 <<= 32;     /* Exact */
  valueLo = value - tempInt64;

  valueDoubleHi = (double) valueHi;
  valueDoubleLo = (double) valueLo;
  valueDoubleHi *= 4294967296.0;   /* Multiply by 2^32 */

  mpfr_init2(valueMpfr, 64);
  mpfr_init2(tempMpfr, 64);
  mpfr_set_d(valueMpfr, valueDoubleHi, GMP_RNDN);
  mpfr_set_d(tempMpfr, valueDoubleLo, GMP_RNDN);
  mpfr_add(valueMpfr, valueMpfr, tempMpfr, GMP_RNDN); /* exact */
  temp = makeConstant(valueMpfr);
  mpfr_clear(valueMpfr);
  mpfr_clear(valueMpfr);

  return temp; 
}

int sollya_lib_get_interval_from_range(sollya_mpfi_t interval, sollya_obj_t obj1) {
  mpfr_t a, b;

  mpfr_init2(a,tools_precision); /* evaluateThingToRange will change the precision afterwards */
  mpfr_init2(b,tools_precision); /* evaluateThingToRange will change the precision afterwards */
  if (evaluateThingToRange(a, b, obj1)) {
    sollya_mpfi_interv_fr(interval, a, b); 
    mpfr_clear(a);
    mpfr_clear(b);
    return 1;
  } 

  mpfr_clear(a);
  mpfr_clear(b);
  return 0;
}

int sollya_lib_get_bounds_from_range(mpfr_t left, mpfr_t right, sollya_obj_t obj1) {
  mpfr_t a, b;
  mp_prec_t p, pp;
  sollya_mpfi_t temp;

  mpfr_init2(a,tools_precision); /* evaluateThingToRange will change the precision afterwards */
  mpfr_init2(b,tools_precision); /* evaluateThingToRange will change the precision afterwards */
  if (evaluateThingToRange(a, b, obj1)) {
    p = mpfr_get_prec(a);
    pp = mpfr_get_prec(b);
    if (pp > p) p = pp;
    sollya_mpfi_init2(temp,p);
    sollya_mpfi_interv_fr(temp, a, b); 
    sollya_mpfi_get_left(left, temp);
    sollya_mpfi_get_right(right, temp);
    sollya_mpfi_clear(temp);
    mpfr_clear(a);
    mpfr_clear(b);
    return 1;
  } 

  mpfr_clear(a);
  mpfr_clear(b);
  return 0;
}

int sollya_lib_get_string(char **str, sollya_obj_t obj1) {
  return evaluateThingToString(str, obj1);
}

int sollya_lib_get_constant(mpfr_t value, sollya_obj_t obj1) {
  sollya_obj_t evaluatedObj, simplifiedObj;

  evaluatedObj = evaluateThing(obj1);
  if (!isPureTree(evaluatedObj)) {
    freeThing(evaluatedObj);
    return 0;
  }

  simplifiedObj = simplifyTreeErrorfree(evaluatedObj);
  if (!isConstant(simplifiedObj)) {
    freeThing(evaluatedObj);
    freeThing(simplifiedObj);
    return 0;
  }

  if (evaluateThingToConstant(value, simplifiedObj, NULL, 1)) {
    freeThing(evaluatedObj);
    freeThing(simplifiedObj);
    return 1;
  }

  freeThing(evaluatedObj);
  freeThing(simplifiedObj);
  return 0;
}

int sollya_lib_get_constant_as_double(double *value, sollya_obj_t obj1) {
  mpfr_t temp;

  mpfr_init2(temp,64); /* sollya_lib_get_constant may change the precision afterwards */
  if (sollya_lib_get_constant(temp, obj1)) {
    *value = mpfr_get_d(temp, GMP_RNDN);
    mpfr_clear(temp);
    return 1;
  } 

  mpfr_clear(temp);
  return 0;
}

int sollya_lib_get_constant_as_int(int *value, sollya_obj_t obj1) {
  mpfr_t temp;

  mpfr_init2(temp,64); /* sollya_lib_get_constant may change the precision afterwards */
  if (sollya_lib_get_constant(temp, obj1)) {
    *value = mpfr_get_si(temp, GMP_RNDN);
    mpfr_clear(temp);
    return 1;
  } 

  mpfr_clear(temp);
  return 0;
}

int sollya_lib_get_constant_as_int64(int64_t *value, sollya_obj_t obj1) {
  return 0; /* TODO */
}

int sollya_lib_get_constant_as_uint64(uint64_t *value, sollya_obj_t obj1) {
  return 0; /* TODO */
}

sollya_obj_t sollya_lib_list(sollya_obj_t objects[], int num) {
  int i;
  chain *tempChain;

  if (num < 1) return makeEmptyList();
  tempChain = NULL;
  for (i=num-1;i>=0;i--) {
    tempChain = addElement(tempChain, copyThing(objects[i]));
  }
  return makeList(tempChain);
}

sollya_obj_t sollya_lib_end_elliptic_list(sollya_obj_t objects[], int num) {
  int i;
  chain *tempChain;

  if (num < 1) return makeEmptyList();
  tempChain = NULL;
  for (i=num-1;i>=0;i--) {
    tempChain = addElement(tempChain, copyThing(objects[i]));
  }
  return makeFinalEllipticList(tempChain);
}

int sollya_lib_get_list_elements(sollya_obj_t *objects[], int *num, int *end_elliptic, sollya_obj_t obj1) {
  sollya_obj_t evaluatedObj;
  int tempVal, i;
  chain *curr;

  evaluatedObj = evaluateThing(obj1);
  if (isEmptyList(evaluatedObj)) {
    *num = 0;
    *end_elliptic = 0;
    freeThing(evaluatedObj);
    return 1;
  }
  if (isPureList(evaluatedObj) || (tempVal = isPureFinalEllipticList(evaluatedObj))) {
    *num = lengthChain(evaluatedObj->arguments);
    *objects = (sollya_obj_t *) safeCalloc(*num,sizeof(sollya_obj_t));
    for (curr=evaluatedObj->arguments,i=0;curr!=NULL;curr=curr->next,i++) {
      (*objects)[i] = copyThing((node *) (curr->value));
    }
    *end_elliptic = tempVal;
    freeThing(evaluatedObj);
    return 1;
  } 
  
  freeThing(evaluatedObj);
  return 0;
}

int sollya_lib_obj_is_function(sollya_obj_t obj1) {
  return isPureTree(obj1);
}

fp_eval_result_t sollya_lib_evaluate_function_at_point(mpfr_t y, sollya_obj_t obj1, mpfr_t x, mpfr_t *cutoff) {
  return FP_EVAL_FAILURE; // TODO
}

ia_eval_result_t sollya_lib_evaluate_funtion_over_interval(sollya_mpfi_t y, sollya_obj_t obj1, sollya_mpfi_t x) {
  return INT_EVAL_FAILURE; // TODO
}

sollya_obj_t sollya_lib_build_function_free_variable() {
  return makeVariable();
}

sollya_obj_t sollya_lib_build_function_add(sollya_obj_t obj1, sollya_obj_t obj2) {
  return makeAdd(obj1,obj2);
}

sollya_obj_t sollya_lib_build_function_sub(sollya_obj_t obj1, sollya_obj_t obj2) {
  return makeSub(obj1,obj2);
}

sollya_obj_t sollya_lib_build_function_mul(sollya_obj_t obj1, sollya_obj_t obj2) {
  return makeMul(obj1,obj2);
}

sollya_obj_t sollya_lib_build_function_div(sollya_obj_t obj1, sollya_obj_t obj2) {
  return makeDiv(obj1,obj2);
}

sollya_obj_t sollya_lib_build_function_sqrt(sollya_obj_t obj1) {
  return makeSqrt(obj1);
}

sollya_obj_t sollya_lib_build_function_exp(sollya_obj_t obj1) {
  return makeExp(obj1); 
}

sollya_obj_t sollya_lib_build_function_log(sollya_obj_t obj1) {
  return makeLog(obj1); 
}

sollya_obj_t sollya_lib_build_function_log2(sollya_obj_t obj1) {
  return makeLog2(obj1); 
}

sollya_obj_t sollya_lib_build_function_log10(sollya_obj_t obj1) {
  return makeLog10(obj1); 
}

sollya_obj_t sollya_lib_build_function_sin(sollya_obj_t obj1) {
  return makeSin(obj1); 
}

sollya_obj_t sollya_lib_build_function_cos(sollya_obj_t obj1) {
  return makeCos(obj1); 
}

sollya_obj_t sollya_lib_build_function_tan(sollya_obj_t obj1) {
  return makeTan(obj1); 
}

sollya_obj_t sollya_lib_build_function_asin(sollya_obj_t obj1) {
  return makeAsin(obj1); 
}

sollya_obj_t sollya_lib_build_function_acos(sollya_obj_t obj1) {
  return makeAcos(obj1); 
}

sollya_obj_t sollya_lib_build_function_atan(sollya_obj_t obj1) {
  return makeAtan(obj1); 
}

sollya_obj_t sollya_lib_build_function_pow(sollya_obj_t obj1, sollya_obj_t obj2) {
  return makePow(obj1, obj2); 
}

sollya_obj_t sollya_lib_build_function_neg(sollya_obj_t obj1) {
  return makeNeg(obj1); 
}

sollya_obj_t sollya_lib_build_function_abs(sollya_obj_t obj1) {
  return makeAbs(obj1); 
}

sollya_obj_t sollya_lib_build_function_double(sollya_obj_t obj1) {
  return makeDouble(obj1); 
}

sollya_obj_t sollya_lib_build_function_single(sollya_obj_t obj1) {
  return makeSingle(obj1); 
}

sollya_obj_t sollya_lib_build_function_quad(sollya_obj_t obj1) {
  return makeQuad(obj1); 
}

sollya_obj_t sollya_lib_build_function_halfprecision(sollya_obj_t obj1) {
  return makeHalfPrecision(obj1); 
}

sollya_obj_t sollya_lib_build_function_double_double(sollya_obj_t obj1) {
  return makeDoubledouble(obj1); 
}

sollya_obj_t sollya_lib_build_function_triple_double(sollya_obj_t obj1) {
  return makeTripledouble(obj1); 
}

sollya_obj_t sollya_lib_build_function_erf(sollya_obj_t obj1) {
  return makeErf(obj1); 
}

sollya_obj_t sollya_lib_build_function_erfc(sollya_obj_t obj1) {
  return makeErfc(obj1); 
}

sollya_obj_t sollya_lib_build_function_log1p(sollya_obj_t obj1) {
  return makeLog1p(obj1); 
}

sollya_obj_t sollya_lib_build_function_expm1(sollya_obj_t obj1) {
  return makeExpm1(obj1); 
}

sollya_obj_t sollya_lib_build_function_doubleextended(sollya_obj_t obj1) {
  return makeDoubleextended(obj1); 
}

sollya_obj_t sollya_lib_build_function_ceil(sollya_obj_t obj1) {
  return makeCeil(obj1); 
}

sollya_obj_t sollya_lib_build_function_floor(sollya_obj_t obj1) {
  return makeFloor(obj1); 
}

sollya_obj_t sollya_lib_build_function_nearestint(sollya_obj_t obj1) {
  return makeNearestInt(obj1); 
}

sollya_obj_t sollya_lib_build_function_sinh(sollya_obj_t obj1) {
  return makeSinh(obj1); 
}

sollya_obj_t sollya_lib_build_function_cosh(sollya_obj_t obj1) {
  return makeCosh(obj1); 
}

sollya_obj_t sollya_lib_build_function_tanh(sollya_obj_t obj1) {
  return makeTanh(obj1); 
}

sollya_obj_t sollya_lib_build_function_asinh(sollya_obj_t obj1) {
  return makeAsinh(obj1); 
}

sollya_obj_t sollya_lib_build_function_acosh(sollya_obj_t obj1) {
  return makeAcosh(obj1); 
}

sollya_obj_t sollya_lib_build_function_atanh(sollya_obj_t obj1) {
  return makeAtanh(obj1); 
}

sollya_obj_t sollya_lib_build_function_pi() {
  return makePi(); 
}

