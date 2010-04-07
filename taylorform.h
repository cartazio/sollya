/*

Copyright 2009 by 

Laboratoire de l'Informatique du Parallélisme, 
UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668

Contributors Ch. Lauter, S. Chevillard, N. Jourdan

christoph.lauter@ens-lyon.fr
sylvain.chevillard@ens-lyon.fr
nicolas.jourdan@ens-lyon.fr

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

#ifndef TAYLORFORM_H
#define TAYLORFORM_H

#include "expression.h"
#include "chain.h"

/* Autodiff related functions */

void exp_diff(mpfi_t *res, mpfi_t x, int n);
void expm1_diff(mpfi_t *res, mpfi_t x, int n);
void log_diff(mpfi_t *res, mpfi_t x, int n);
void log2_diff(mpfi_t *res, mpfi_t x, int n);
void log10_diff(mpfi_t *res, mpfi_t x, int n);
void sin_diff(mpfi_t *res, mpfi_t x, int n);
void cos_diff(mpfi_t *res, mpfi_t x, int n);
void sinh_diff(mpfi_t *res, mpfi_t x, int n);
void cosh_diff(mpfi_t *res, mpfi_t x, int n);
void tan_diff(mpfi_t *res, mpfi_t x, int n);
void tanh_diff(mpfi_t *res, mpfi_t x, int n);
void atan_diff(mpfi_t *res, mpfi_t x, int n);
void atanh_diff(mpfi_t *res, mpfi_t x, int n);
void asin_diff(mpfi_t *res, mpfi_t x, int n);
void acos_diff(mpfi_t *res, mpfi_t x, int n);
void asinh_diff(mpfi_t *res, mpfi_t x, int n);
void acosh_diff(mpfi_t *res, mpfi_t x, int n);

void powerFunction_diff(mpfi_t *res, mpfr_t p, mpfi_t x, int n);
void constantPower_diff(mpfi_t *res, mpfi_t x, mpfr_t p, int n);
void baseFunction_diff(mpfi_t *res, int nodeType, mpfi_t x, int n);


void taylorform(node **T, chain **errors, mpfi_t **delta,
		node *f, int n, mpfi_t *x0, mpfi_t *d, int mode);

#endif /* TAYLORFORM_H */
