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
therefore means  that it is reserved for developers  and  experienced
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

#ifndef SOLLYA_LIBRARY_WRAPPERS_H
#define SOLLYA_LIBRARY_WRAPPERS_H

#include <gmp.h>
#include <mpfr.h>
#include <stdint.h>
#include <stdio.h>
#include "expression.h"
#include "execute.h"
#include "mpfi-compat.h"

/* Define a type for all Sollya objects */
typedef node * sollya_obj_t;

/* Define an enumeration type for the status
   of floating-point evaluation
*/
typedef enum fp_eval_result_enum_t fp_eval_result_t;
enum fp_eval_result_enum_t {
  FP_EVAL_OBJ_NO_FUNCTION,
  FP_EVAL_FAITHFUL,
  FP_EVAL_BELOW_CUTOFF,
  FP_EVAL_NOT_FAITHFUL_ZERO,
  FP_EVAL_NOT_FAITHFUL_NOT_ZERO,
  FP_EVAL_FAILURE
};

/* Define an enumeration type for the status
   of floating-point evaluation
*/
typedef enum ia_eval_result_enum_t ia_eval_result_t;
enum ia_eval_result_enum_t {
  INT_EVAL_OBJ_NO_FUNCTION,
  INT_EVAL_BOUNDED,
  INT_EVAL_UNBOUNDED,
  INT_EVAL_FAILURE
};

/* Initialization and finitialization functions */
int sollya_lib_init();
int sollya_lib_close();

/* Function to install and uninstall a call-back for the messages
   emitted by the Sollya core.  
*/
int sollya_lib_install_msg_callback(int (*) (int));
int sollya_lib_uninstall_msg_callback();

/* Functions to print anything, including Sollya objects */
int sollya_lib_printf(const char *, ...);
int sollya_lib_fprintf(FILE *, const char *, ...);

/* A function to clear Sollya_objects */
void sollya_lib_clear_obj(sollya_obj_t);

/* A function to structurally compare two Sollya objects */
int sollya_lib_cmp_objs_structurally(sollya_obj_t, sollya_obj_t);

/* A function to copy Sollya objects */
sollya_obj_t sollya_lib_copy_obj(sollya_obj_t);

/* Functions corresponding to Sollya commands */
void sollya_lib_plot(sollya_obj_t, sollya_obj_t, ...);
void sollya_lib_printdouble(sollya_obj_t);
void sollya_lib_printsingle(sollya_obj_t);
void sollya_lib_printexpansion(sollya_obj_t);
void sollya_lib_implementconst(sollya_obj_t, ...);
void sollya_lib_bashexecute(sollya_obj_t);
void sollya_lib_externalplot(sollya_obj_t, sollya_obj_t, sollya_obj_t, sollya_obj_t, sollya_obj_t, ...);
void sollya_lib_asciiplot(sollya_obj_t, sollya_obj_t);
void sollya_lib_printxml(sollya_obj_t);
void sollya_lib_printxml_newfile(sollya_obj_t, sollya_obj_t);
void sollya_lib_printxml_appendfile(sollya_obj_t, sollya_obj_t);
void sollya_lib_worstcase(sollya_obj_t, sollya_obj_t, sollya_obj_t, sollya_obj_t, sollya_obj_t, ...);
void sollya_lib_autoprint(sollya_obj_t, ...);
void sollya_lib_set_prec(sollya_obj_t);
void sollya_lib_set_points(sollya_obj_t);
void sollya_lib_set_diam(sollya_obj_t);
void sollya_lib_set_display(sollya_obj_t);
void sollya_lib_set_verbosity(sollya_obj_t);
void sollya_lib_set_canonical(sollya_obj_t);
void sollya_lib_set_autosimplify(sollya_obj_t);
void sollya_lib_set_taylorrecursions(sollya_obj_t);
void sollya_lib_set_timing(sollya_obj_t);
void sollya_lib_set_midpointmode(sollya_obj_t);
void sollya_lib_set_dieonerrormode(sollya_obj_t);
void sollya_lib_set_rationalmode(sollya_obj_t);
void sollya_lib_set_roundingwarnings(sollya_obj_t);
void sollya_lib_set_hopitalrecursions(sollya_obj_t);
void sollya_lib_set_prec_silent(sollya_obj_t);
void sollya_lib_set_points_silent(sollya_obj_t);
void sollya_lib_set_diam_silent(sollya_obj_t);
void sollya_lib_set_display_silent(sollya_obj_t);
void sollya_lib_set_verbosity_silent(sollya_obj_t);
void sollya_lib_set_canonical_silent(sollya_obj_t);
void sollya_lib_set_autosimplify_silent(sollya_obj_t);
void sollya_lib_set_taylorrecursions_silent(sollya_obj_t);
void sollya_lib_set_timing_silent(sollya_obj_t);
void sollya_lib_set_midpointmode_silent(sollya_obj_t);
void sollya_lib_set_dieonerrormode_silent(sollya_obj_t);
void sollya_lib_set_rationalmode_silent(sollya_obj_t);
void sollya_lib_set_roundingwarnings_silent(sollya_obj_t);
void sollya_lib_set_hopitalrecursions_silent(sollya_obj_t);

/* Functions corresponding to Sollya built-in procedures */
sollya_obj_t sollya_lib_free_variable();
sollya_obj_t sollya_lib_and(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_or(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_negate(sollya_obj_t);
sollya_obj_t sollya_lib_cmp_equal(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_cmp_in(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_cmp_less(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_cmp_greater(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_cmp_less_equal(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_cmp_greater_equal(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_cmp_not_equal(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_add(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_sub(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_concat(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_append(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_prepend(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_apply(sollya_obj_t, sollya_obj_t, ...);
sollya_obj_t sollya_lib_approx(sollya_obj_t);
sollya_obj_t sollya_lib_mul(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_div(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_pow(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_minus(sollya_obj_t);
sollya_obj_t sollya_lib_sup(sollya_obj_t);
sollya_obj_t sollya_lib_mid(sollya_obj_t);
sollya_obj_t sollya_lib_inf(sollya_obj_t);
sollya_obj_t sollya_lib_diff(sollya_obj_t);
sollya_obj_t sollya_lib_simplify(sollya_obj_t);
sollya_obj_t sollya_lib_bashevaluate(sollya_obj_t, ...);
sollya_obj_t sollya_lib_remez(sollya_obj_t, sollya_obj_t, sollya_obj_t, ...);
sollya_obj_t sollya_lib_min(sollya_obj_t, ...);
sollya_obj_t sollya_lib_max(sollya_obj_t, ...);
sollya_obj_t sollya_lib_fpminimax(sollya_obj_t, sollya_obj_t, sollya_obj_t, sollya_obj_t, ...);
sollya_obj_t sollya_lib_horner(sollya_obj_t);
sollya_obj_t sollya_lib_canonical(sollya_obj_t);
sollya_obj_t sollya_lib_expand(sollya_obj_t);
sollya_obj_t sollya_lib_simplifysafe(sollya_obj_t);
sollya_obj_t sollya_lib_taylor(sollya_obj_t, sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_taylorform(sollya_obj_t, sollya_obj_t, sollya_obj_t, ...);
sollya_obj_t sollya_lib_autodiff(sollya_obj_t, sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_degree(sollya_obj_t);
sollya_obj_t sollya_lib_numerator(sollya_obj_t);
sollya_obj_t sollya_lib_denominator(sollya_obj_t);
sollya_obj_t sollya_lib_substitute(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_composepolynomials(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_coeff(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_subpoly(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_roundcoefficients(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_rationalapprox(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_round(sollya_obj_t, sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_evaluate(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_parse(sollya_obj_t);
sollya_obj_t sollya_lib_readxml(sollya_obj_t);
sollya_obj_t sollya_lib_infnorm(sollya_obj_t, sollya_obj_t, ...);
sollya_obj_t sollya_lib_supnorm(sollya_obj_t, sollya_obj_t, sollya_obj_t, sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_findzeros(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_dirtyinfnorm(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_numberroots(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_integral(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_dirtyintegral(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_implementpoly(sollya_obj_t, sollya_obj_t, sollya_obj_t, sollya_obj_t, sollya_obj_t, sollya_obj_t, ...);
sollya_obj_t sollya_lib_checkinfnorm(sollya_obj_t, sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_zerodenominators(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_searchgal(sollya_obj_t, ...);
sollya_obj_t sollya_lib_guessdegree(sollya_obj_t, sollya_obj_t, sollya_obj_t, ...);
sollya_obj_t sollya_lib_dirtyfindzeros(sollya_obj_t, sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_head(sollya_obj_t);
sollya_obj_t sollya_lib_roundcorrectly(sollya_obj_t);
sollya_obj_t sollya_lib_revert(sollya_obj_t);
sollya_obj_t sollya_lib_sort(sollya_obj_t);
sollya_obj_t sollya_lib_mantissa(sollya_obj_t);
sollya_obj_t sollya_lib_exponent(sollya_obj_t);
sollya_obj_t sollya_lib_tail(sollya_obj_t);
sollya_obj_t sollya_lib_range(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_sqrt(sollya_obj_t);
sollya_obj_t sollya_lib_exp(sollya_obj_t);
sollya_obj_t sollya_lib_log(sollya_obj_t);
sollya_obj_t sollya_lib_log2(sollya_obj_t);
sollya_obj_t sollya_lib_log10(sollya_obj_t);
sollya_obj_t sollya_lib_sin(sollya_obj_t);
sollya_obj_t sollya_lib_cos(sollya_obj_t);
sollya_obj_t sollya_lib_tan(sollya_obj_t);
sollya_obj_t sollya_lib_asin(sollya_obj_t);
sollya_obj_t sollya_lib_acos(sollya_obj_t);
sollya_obj_t sollya_lib_atan(sollya_obj_t);
sollya_obj_t sollya_lib_sinh(sollya_obj_t);
sollya_obj_t sollya_lib_cosh(sollya_obj_t);
sollya_obj_t sollya_lib_tanh(sollya_obj_t);
sollya_obj_t sollya_lib_asinh(sollya_obj_t);
sollya_obj_t sollya_lib_acosh(sollya_obj_t);
sollya_obj_t sollya_lib_atanh(sollya_obj_t);
sollya_obj_t sollya_lib_abs(sollya_obj_t);
sollya_obj_t sollya_lib_erf(sollya_obj_t);
sollya_obj_t sollya_lib_erfc(sollya_obj_t);
sollya_obj_t sollya_lib_log1p(sollya_obj_t);
sollya_obj_t sollya_lib_expm1(sollya_obj_t);
sollya_obj_t sollya_lib_double(sollya_obj_t);
sollya_obj_t sollya_lib_single(sollya_obj_t);
sollya_obj_t sollya_lib_quad(sollya_obj_t);
sollya_obj_t sollya_lib_halfprecision(sollya_obj_t);
sollya_obj_t sollya_lib_double_double(sollya_obj_t);
sollya_obj_t sollya_lib_triple_double(sollya_obj_t);
sollya_obj_t sollya_lib_doubleextended(sollya_obj_t);
sollya_obj_t sollya_lib_ceil(sollya_obj_t);
sollya_obj_t sollya_lib_floor(sollya_obj_t);
sollya_obj_t sollya_lib_nearestint(sollya_obj_t);
sollya_obj_t sollya_lib_length(sollya_obj_t);
sollya_obj_t sollya_lib_get_prec();
sollya_obj_t sollya_lib_get_points();
sollya_obj_t sollya_lib_get_diam();
sollya_obj_t sollya_lib_get_display();
sollya_obj_t sollya_lib_get_verbosity();
sollya_obj_t sollya_lib_get_canonical();
sollya_obj_t sollya_lib_get_autosimplify();
sollya_obj_t sollya_lib_get_taylorrecursions();
sollya_obj_t sollya_lib_get_timing();
sollya_obj_t sollya_lib_get_midpointmode();
sollya_obj_t sollya_lib_get_dieonerrormode();
sollya_obj_t sollya_lib_get_rationalmode();
sollya_obj_t sollya_lib_get_roundingwarnings();
sollya_obj_t sollya_lib_get_hopitalrecursions();

/* Functions creating Sollya objects */
sollya_obj_t sollya_lib_on();
sollya_obj_t sollya_lib_off();
sollya_obj_t sollya_lib_dyadic();
sollya_obj_t sollya_lib_powers();
sollya_obj_t sollya_lib_binary();
sollya_obj_t sollya_lib_hexadecimal();
sollya_obj_t sollya_lib_file();
sollya_obj_t sollya_lib_postscript();
sollya_obj_t sollya_lib_postscriptfile();
sollya_obj_t sollya_lib_perturb();
sollya_obj_t sollya_lib_round_down();
sollya_obj_t sollya_lib_round_up();
sollya_obj_t sollya_lib_round_towards_zero();
sollya_obj_t sollya_lib_round_to_nearest();
sollya_obj_t sollya_lib_honorcoeffprec();
sollya_obj_t sollya_lib_true();
sollya_obj_t sollya_lib_false();
sollya_obj_t sollya_lib_void();
sollya_obj_t sollya_lib_default();
sollya_obj_t sollya_lib_decimal();
sollya_obj_t sollya_lib_absolute();
sollya_obj_t sollya_lib_relative();
sollya_obj_t sollya_lib_fixed();
sollya_obj_t sollya_lib_floating();
sollya_obj_t sollya_lib_error();
sollya_obj_t sollya_lib_double_obj();
sollya_obj_t sollya_lib_single_obj();
sollya_obj_t sollya_lib_quad_obj();
sollya_obj_t sollya_lib_halfprecision_obj();
sollya_obj_t sollya_lib_doubleextended_obj();
sollya_obj_t sollya_lib_double_double_obj();
sollya_obj_t sollya_lib_triple_double_obj();
sollya_obj_t sollya_lib_pi();

/* A function to parse expressions and evaluate them */
sollya_obj_t sollya_lib_parse_string(char *);

/* Functions to convert from constants to Sollya objects */
sollya_obj_t sollya_lib_string(char *);
sollya_obj_t sollya_lib_range_from_interval(sollya_mpfi_t);
sollya_obj_t sollya_lib_range_from_bounds(mpfr_t, mpfr_t);
sollya_obj_t sollya_lib_constant(mpfr_t);
sollya_obj_t sollya_lib_constant_from_double(double);
sollya_obj_t sollya_lib_constant_from_int(int);
sollya_obj_t sollya_lib_constant_from_int64(int64_t);
sollya_obj_t sollya_lib_constant_from_uint64(uint64_t);

/* Functions to get values contained in Sollya objects */
int sollya_lib_get_interval_from_range(sollya_mpfi_t, sollya_obj_t);
int sollya_lib_get_bounds_from_range(mpfr_t, mpfr_t, sollya_obj_t);
int sollya_lib_get_string(char **, sollya_obj_t);
int sollya_lib_get_constant_as_double(double *, sollya_obj_t);
int sollya_lib_get_constant_as_int(int *, sollya_obj_t);
int sollya_lib_get_constant_as_int64(int64_t *, sollya_obj_t);
int sollya_lib_get_constant_as_uint64(uint64_t *, sollya_obj_t);

/* The following function, in contrast to all others, 
   not only assigns a new value to the mpfr_t argument
   in case of success but also adjusts its precision
   in order to store the constant in the Sollya 
   object exactly (without any rounding).

   Rounding may nevertheless happen if the Sollya object 
   is not a constant by itself but a constant expression
   that needs to be evaluated. 
*/
int sollya_lib_get_constant(mpfr_t, sollya_obj_t);

/* Functions to build up Sollya lists from arrays of objects and 
   to get arrays of Sollya objects out of Sollya lists 
*/
sollya_obj_t sollya_lib_list(sollya_obj_t[], int);
sollya_obj_t sollya_lib_end_elliptic_list(sollya_obj_t[], int);
int sollya_lib_get_list_elements(sollya_obj_t *[], int *, int *, sollya_obj_t);

/* A function to check if a Sollya object represents a mathematical
   function 
*/
int sollya_lib_obj_is_function(sollya_obj_t);

/* Functions to evaluate Sollya objects that are mathematical
   functions at points or over intervals 
*/
fp_eval_result_t sollya_lib_evaluate_function_at_point(mpfr_t, sollya_obj_t, mpfr_t, mpfr_t *);
ia_eval_result_t sollya_lib_evaluate_funtion_over_interval(sollya_mpfi_t, sollya_obj_t, sollya_mpfi_t);

/* Functions for building Sollya objects representing 
   mathematical functions.

   Attention: in contrast to all other functions in
   the Sollya library, these functions "use up" the
   objects they take as an argument.

 */
sollya_obj_t sollya_lib_build_function_free_variable();
sollya_obj_t sollya_lib_build_function_add(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_build_function_sub(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_build_function_mul(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_build_function_div(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_build_function_sqrt(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_exp(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_log(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_log2(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_log10(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_sin(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_cos(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_tan(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_asin(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_acos(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_atan(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_pow(sollya_obj_t, sollya_obj_t);
sollya_obj_t sollya_lib_build_function_neg(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_abs(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_double(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_single(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_quad(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_halfprecision(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_double_double(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_triple_double(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_erf(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_erfc(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_log1p(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_expm1(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_doubleextended(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_ceil(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_floor(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_nearestint(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_sinh(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_cosh(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_tanh(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_asinh(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_acosh(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_atanh(sollya_obj_t);
sollya_obj_t sollya_lib_build_function_pi();

#endif /* ifdef SOLLYA_LIBRARY_WRAPPERS_H*/
