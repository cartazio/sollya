/*

Copyright 2008-2010 by 

Laboratoire de l'Informatique du Parallélisme, 
UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668

and by

LORIA (CNRS, INPL, INRIA, UHP, U-Nancy 2)

Contributors S. Chevillard, M. Joldes, Ch. Lauter

sylvain.chevillard@ens-lyon.org
mioara.joldes@ens-lyon.fr
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

*/

#ifndef AUTODIFF_H
#define AUTODIFF_H

#include <mpfi.h>
#include "expression.h"

void symbolic_poly_diff(mpfi_t *res, mpfi_t *coeff_array, int degree);
void symbolic_poly_evaluation_horner(mpfi_t res, mpfi_t *coeffs_array, mpfi_t x, int degree);
void symbolic_poly_evaluation_powers(mpfi_t res, mpfi_t *coeffs_array, mpfi_t *powers_array, mpfi_t x, int degree);

void exp_diff(mpfi_t *res, mpfi_t x0, int n);
void expm1_diff(mpfi_t *res, mpfi_t x0, int n);
void log1p_diff(mpfi_t *res, mpfi_t x0, int n);
void log_diff(mpfi_t *res, mpfi_t x0, int n);
void log2_diff(mpfi_t *res, mpfi_t x0, int n);
void log10_diff(mpfi_t *res, mpfi_t x0, int n);
void sin_diff(mpfi_t *res, mpfi_t x0, int n);
void cos_diff(mpfi_t *res, mpfi_t x0, int n);
void sinh_diff(mpfi_t *res, mpfi_t x0, int n);
void cosh_diff(mpfi_t *res, mpfi_t x0, int n);
void tan_diff(mpfi_t *res, mpfi_t x0, int n);
void tanh_diff(mpfi_t *res, mpfi_t x0, int n);
void atan_diff(mpfi_t *res, mpfi_t x0, int n);
void atanh_diff(mpfi_t *res, mpfi_t x0, int n);
void asin_diff(mpfi_t *res, mpfi_t x0, int n);
void acos_diff(mpfi_t *res, mpfi_t x0, int n);
void asinh_diff(mpfi_t *res, mpfi_t x0, int n);
void acosh_diff(mpfi_t *res, mpfi_t x0, int n);
void erf_diff(mpfi_t *res, mpfi_t x0, int n);
void erfc_diff(mpfi_t *res, mpfi_t x0, int n);
void abs_diff(mpfi_t *res, mpfi_t x0, int n);
void ceil_diff(mpfi_t *res, mpfi_t x, int n);
void double_diff(mpfi_t *res, mpfi_t x, int n);
void double_double_diff(mpfi_t *res, mpfi_t x, int n);
void double_extended_diff(mpfi_t *res, mpfi_t x, int n);
void floor_diff(mpfi_t *res, mpfi_t x, int n);
void nearestint_diff(mpfi_t *res, mpfi_t x, int n);
void single_diff(mpfi_t *res, mpfi_t x, int n);
void triple_double_diff(mpfi_t *res, mpfi_t x, int n);

void powerFunction_diff(mpfi_t *res, mpfr_t p, mpfi_t x0, int n);
void constantPower_diff(mpfi_t *res, mpfi_t x0, mpfr_t p, int n);
void baseFunction_diff(mpfi_t *res, int nodeType, mpfi_t x0, int n);

void multiplication_AD(mpfi_t *res, mpfi_t *f, mpfi_t *g, int n);
void composition_AD(mpfi_t *res, mpfi_t *g, mpfi_t *f, int n);
void auto_diff_scaled(mpfi_t* res, node *f, mpfi_t x0, int n);
void auto_diff(mpfi_t* res, node *f, mpfi_t x0, int n);

#endif /* AUTODIFF_H */
