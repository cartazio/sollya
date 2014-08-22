/*

  Copyright 2014 by

  Centre de recherche INRIA Sophia-Antipolis Mediterranee, equipe APICS,
  Sophia Antipolis, France.

  Contributors S. Chevillard

  sylvain.chevillard@ens-lyon.org

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

#include <stdlib.h>
#include "general.h"
#include "expression.h"
#include "base-functions.h"
#include "autodiff.h"
#include "mpfi-compat.h"
#include "infnorm.h"
#include "double.h"


/******************************************************************************/
/*                                                                            */
/*                    AUTOMATIC DIFFERENTIATION FUNCTIONS                     */
/*                                                                            */
/******************************************************************************/

void sqrt_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  mpfr_t oneHalf;
  mpfr_init2(oneHalf, 12);
  mpfr_set_d(oneHalf, 0.5, GMP_RNDN);
  constantPower_diff(res, x, oneHalf, n, silent);
  mpfr_clear(oneHalf);
  return;
}

void exp_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i;
  sollya_mpfi_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  sollya_mpfi_init2(temp, prec);

  sollya_mpfi_exp(temp, x);
  for(i=0;i<=n;i++) {
    sollya_mpfi_set(res[i], temp);
    sollya_mpfi_div_ui(temp, temp, i+1);
  }

  sollya_mpfi_clear(temp);
  return;
}

void expm1_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  exp_diff(res, x, n, silent);
  sollya_mpfi_sub_ui(res[0], res[0], 1);
  return;
}


void log_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  mpfr_t minusOne;
  mp_prec_t prec;
  int i;

  prec = getToolPrecision();

  sollya_mpfi_log(res[0], x);

  if(n>=1) {
    mpfr_init2(minusOne, prec);
    mpfr_set_si(minusOne, -1, GMP_RNDN);
    constantPower_diff(res+1, x, minusOne, n-1, silent);
    mpfr_clear(minusOne);
  }
  for(i=1;i<=n;i++) sollya_mpfi_div_ui(res[i], res[i], i);
  return;
}

void log1p_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  mpfr_t minusOne;
  sollya_mpfi_t u;
  int i;
  mp_prec_t prec;

  prec = getToolPrecision();

  sollya_mpfi_log1p(res[0], x);

  if(n>=1) {
    sollya_mpfi_init2(u, prec);
    sollya_mpfi_add_ui(u, x, 1);
    mpfr_init2(minusOne, prec);
    mpfr_set_si(minusOne, -1, GMP_RNDN);
    constantPower_diff(res+1, u, minusOne, n-1, silent);
    mpfr_clear(minusOne);
    sollya_mpfi_clear(u);
  }

  for(i=1;i<=n;i++) sollya_mpfi_div_ui(res[i], res[i], i);

  return;
}

/* log2(x) = log(x) * (1/log(2)) */
void log2_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i;
  sollya_mpfi_t log2;
  mp_prec_t prec;

  prec = getToolPrecision();
  sollya_mpfi_init2(log2, prec);

  sollya_mpfi_set_ui(log2, 2); sollya_mpfi_log(log2, log2);
  log_diff(res,x,n,silent);
  for(i=0;i<=n;i++) sollya_mpfi_div(res[i], res[i], log2);

  sollya_mpfi_clear(log2);
  return;
}

/* idem */
void log10_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i;
  sollya_mpfi_t log10;
  mp_prec_t prec;

  prec = getToolPrecision();
  sollya_mpfi_init2(log10, prec);

  sollya_mpfi_set_ui(log10, 10); sollya_mpfi_log(log10, log10);
  log_diff(res,x,n,silent);
  for(i=0;i<=n;i++) sollya_mpfi_div(res[i], res[i], log10);

  sollya_mpfi_clear(log10);
  return;
}

void sin_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i;

  UNUSED_PARAM(silent);

  sollya_mpfi_sin(res[0], x);
  for(i=2; i<=n; i+=2) sollya_mpfi_div_ui(res[i], res[i-2], (i-1)*i);
  for(i=2; i<=n; i +=4) sollya_mpfi_neg(res[i], res[i]);

  if(n>=1) {
    sollya_mpfi_cos(res[1], x);
    for(i=3; i<=n; i+=2) sollya_mpfi_div_ui(res[i], res[i-2], (i-1)*i);
    for(i=3; i<=n; i +=4) sollya_mpfi_neg(res[i], res[i]);
  }

  return;
}

void cos_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i;

  UNUSED_PARAM(silent);

  sollya_mpfi_cos(res[0], x);
  for(i=2; i<=n; i+=2) sollya_mpfi_div_ui(res[i], res[i-2], (i-1)*i);
  for(i=2; i<=n; i +=4) sollya_mpfi_neg(res[i], res[i]);

  if(n>=1) {
    sollya_mpfi_sin(res[1], x);
    for(i=3; i<=n; i+=2) sollya_mpfi_div_ui(res[i], res[i-2], (i-1)*i);
    for(i=1; i<=n; i +=4) sollya_mpfi_neg(res[i], res[i]);
  }

  return;
}

void sinh_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i;

  UNUSED_PARAM(silent);

  sollya_mpfi_sinh(res[0], x);
  for(i=2; i<=n; i+=2) sollya_mpfi_div_ui(res[i], res[i-2], (i-1)*i);

  if(n>=1) {
    sollya_mpfi_cosh(res[1], x);
    for(i=3; i<=n; i+=2) sollya_mpfi_div_ui(res[i], res[i-2], (i-1)*i);
  }

  return;
}

void cosh_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i;

  UNUSED_PARAM(silent);

  sollya_mpfi_cosh(res[0], x);
  for(i=2; i<=n; i+=2) sollya_mpfi_div_ui(res[i], res[i-2], (i-1)*i);

  if(n>=1) {
    sollya_mpfi_sinh(res[1], x);
    for(i=3; i<=n; i+=2) sollya_mpfi_div_ui(res[i], res[i-2], (i-1)*i);
  }

  return;
}

/*  u=tan(x), tan^(n) / n! = p_(n)(u) with
    p_0 = u;

    recurrence formula: p_(n+1)(u) = (p_n(u))' / (n+1) = p'_n(u) * (1+u^2) / (n+1)
    -> p_n of degree n+1
*/

void tan_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i,index;
  sollya_mpfi_t *coeffs_array;
  sollya_mpfi_t u;
  mp_prec_t prec;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  coeffs_array = (sollya_mpfi_t *)safeCalloc( (n+2),sizeof(sollya_mpfi_t));

  for (index=0; index<=n+1; index++) {
    sollya_mpfi_init2(coeffs_array[index], prec);
    sollya_mpfi_set_ui(coeffs_array[index], 0);
  }
  sollya_mpfi_init2(u, prec);

  sollya_mpfi_tan(u, x);
  sollya_mpfi_set_ui(coeffs_array[0], 0);
  sollya_mpfi_set_ui(coeffs_array[1], 1);

  sollya_mpfi_set(res[0], u);

  for(index=1; index<=n; index++) {
    /* coeffs_array represents p_(index-1) */

    symbolic_poly_diff(coeffs_array, coeffs_array, index);
    sollya_mpfi_set_ui(coeffs_array[index], 0);
    /* now it represents p_(index-1)' */

    for(i=index+1; i>=2; i--) {
      sollya_mpfi_add(coeffs_array[i], coeffs_array[i], coeffs_array[i-2]);
      sollya_mpfi_div_ui(coeffs_array[i], coeffs_array[i], index);
    }
    sollya_mpfi_div_ui(coeffs_array[1], coeffs_array[1], index);
    sollya_mpfi_div_ui(coeffs_array[0], coeffs_array[0], index);
    /* now it represents p_(index) */

    symbolic_poly_evaluation_horner(res[index], coeffs_array, u, index+1);
  }

  for (index=0; index<=n+1; index++){
    sollya_mpfi_clear(coeffs_array[index]);
  }
  sollya_mpfi_clear(u);
  safeFree(coeffs_array);

  return;
}

/*  u=tanh(x), tanh^(n) / n! = p_(n)(u) with
    p_0 = u;

    recurrence formula: p_(n+1)(u) = (p_n(u))' / (n+1) = p'_n(u) * (1-u^2) / (n+1)
    -> p_n of degree n+1
*/

void tanh_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i,index;
  sollya_mpfi_t *coeffs_array;
  sollya_mpfi_t u;
  mp_prec_t prec;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  coeffs_array = (sollya_mpfi_t *)safeCalloc( (n+2),sizeof(sollya_mpfi_t));

  for (index=0; index<=n+1; index++) {
    sollya_mpfi_init2(coeffs_array[index], prec);
    sollya_mpfi_set_ui(coeffs_array[index], 0);
  }
  sollya_mpfi_init2(u, prec);

  sollya_mpfi_tanh(u, x);
  sollya_mpfi_set_ui(coeffs_array[0], 0);
  sollya_mpfi_set_ui(coeffs_array[1], 1);

  sollya_mpfi_set(res[0], u);

  for(index=1; index<=n; index++) {
    /* coeffs_array represents p_(index-1) */

    symbolic_poly_diff(coeffs_array, coeffs_array, index);
    sollya_mpfi_set_ui(coeffs_array[index], 0);
    /* now it represents p_(index-1)' */

    for(i=index+1; i>=2; i--) {
      sollya_mpfi_sub(coeffs_array[i], coeffs_array[i], coeffs_array[i-2]);
      sollya_mpfi_div_ui(coeffs_array[i], coeffs_array[i], index);
    }
    sollya_mpfi_div_ui(coeffs_array[1], coeffs_array[1], index);
    sollya_mpfi_div_ui(coeffs_array[0], coeffs_array[0], index);
    /* now it represents p_(index) */

    symbolic_poly_evaluation_horner(res[index], coeffs_array, u, index+1);
  }

  for (index=0; index<=n+1; index++){
    sollya_mpfi_clear(coeffs_array[index]);
  }
  sollya_mpfi_clear(u);
  safeFree(coeffs_array);

  return;
}

/*
   Some thoughts:
   For all holonomic functions, but especially for those below, it would certainly be more interesting to use the recurrence relation that exists between the successive values of the derivative at point x, instead of using a recurrence on the polynomials giving the closed form of the n-th derivative.

   The fact is that the recurrence on the polynomials is more general, as it gives an expression for the derivative *at any point*, but since we are only interested at evaluating the successive derivatives at one single point, it is probably not a good idea as each polynomial requires O(n) operations to be constructed and evaluated, hence a total O(degree^2) operations. This should be compared with O(degree) operations with the recurrence on the coefficients, but it might be that it is not suitable for interval computations, as it may involve recurrence with two or more terms.

   Moreover, the code could be automatically generated (this of few interest here, but could turn very interesting if it were to be recoded inside COQ for instance) for any holonomic function. Here is something I began to write on this topic, on the example of asin:
   deq := holexprtodiffeq(arcsin(x+x0),y(x));
   rec := diffeqtorec(deq,y(x),a(n));
   recc := op(1, rec); cond := [seq(op(i, rec), i=2..nops(rec))];
   recc := collect(recc,a(n)); # a(n) is asin^(n)(x0)/n! and recc is the recurrence it satisfies
   recc := add(factor(op(i,recc)),i=1..nops(recc));
*/



/* atan_diff : reccurence formula: p_(n+1) = (p'_n * (1+x^2) - 2nx * p_n) / (n+1)
   atan^(0) = atan(x)
   atan^(n) / n! = p_(n)(x)/((1+x^2)^n)
   p_1=1;

   --> degree of p_n is (n-1)
*/
void atan_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i,index;
  sollya_mpfi_t *coeffs_array, *coeffs_array_diff;
  sollya_mpfi_t u, temp;

  mp_prec_t prec;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  coeffs_array = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));
  coeffs_array_diff = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));

  for (index=0; index<=n-1; index++) {
    sollya_mpfi_init2(coeffs_array[index], prec);
    sollya_mpfi_init2(coeffs_array_diff[index], prec);

    sollya_mpfi_set_ui(coeffs_array[index], 0);
    sollya_mpfi_set_ui(coeffs_array_diff[index], 0);
  }

  sollya_mpfi_init2(u, prec);
  sollya_mpfi_init2(temp, prec);

  sollya_mpfi_atan(res[0], x);

  if(n>=1) {
    sollya_mpfi_sqr(u, x);
    sollya_mpfi_add_ui(u, u, 1);

    sollya_mpfi_inv(res[1], u);

    sollya_mpfi_set_ui(coeffs_array[0], 1);
  }

  for(index=2; index<=n; index++) {
    /* coeffs_array represents p_(index-1) */

    symbolic_poly_diff(coeffs_array_diff, coeffs_array, index-2);
    sollya_mpfi_set_ui(coeffs_array_diff[index-2], 0);
    /* now it represents p_(index-1)' */

    for(i=index-1; i>=2; i--) {
      sollya_mpfi_add(coeffs_array[i], coeffs_array_diff[i], coeffs_array_diff[i-2]);
      sollya_mpfi_mul_ui(temp, coeffs_array[i-1], 2*(index-1));
      sollya_mpfi_sub(coeffs_array[i], coeffs_array[i], temp);
      sollya_mpfi_div_ui(coeffs_array[i], coeffs_array[i], index);
    }

    sollya_mpfi_mul_ui(temp, coeffs_array[0], 2*(index-1));
    sollya_mpfi_sub(coeffs_array[1], coeffs_array_diff[1], temp);
    sollya_mpfi_div_ui(coeffs_array[1], coeffs_array[1], index);

    sollya_mpfi_set(coeffs_array[0], coeffs_array_diff[0]);
    sollya_mpfi_div_ui(coeffs_array[0], coeffs_array[0], index);
    /* now it represents p_(index) */

    symbolic_poly_evaluation_horner(res[index], coeffs_array, x, index-1);
    sollya_mpfi_set_ui(temp, index);
    sollya_mpfi_pow(temp, u, temp);
    sollya_mpfi_div(res[index], res[index], temp);
  }

  for (index=0; index<=n-1; index++){
    sollya_mpfi_clear(coeffs_array[index]);
    sollya_mpfi_clear(coeffs_array_diff[index]);
  }
  sollya_mpfi_clear(u);
  sollya_mpfi_clear(temp);
  safeFree(coeffs_array);
  safeFree(coeffs_array_diff);

  return;
}


/* atanh_diff : reccurence formula: p_(n+1) = (p'_n * (1-x^2) + 2nx * p_n)/ (n+1)
   atanh^(0) = atanh(x)
   atanh^(n)/n! = p_(n)(x)/((1-x^2)^n)
   p_1=1;

   --> degree of p_n is (n-1)
*/
void atanh_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i,index;
  sollya_mpfi_t *coeffs_array, *coeffs_array_diff;
  sollya_mpfi_t u, temp;

  mp_prec_t prec;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  coeffs_array = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));
  coeffs_array_diff = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));

  for (index=0; index<=n-1; index++) {
    sollya_mpfi_init2(coeffs_array[index], prec);
    sollya_mpfi_init2(coeffs_array_diff[index], prec);

    sollya_mpfi_set_ui(coeffs_array[index], 0);
    sollya_mpfi_set_ui(coeffs_array_diff[index], 0);
  }

  sollya_mpfi_init2(u, prec);
  sollya_mpfi_init2(temp, prec);

  sollya_mpfi_atanh(res[0], x);

  if(n>=1) {
    sollya_mpfi_sqr(u, x);
    sollya_mpfi_sub_ui(u, u, 1); /* TODO: FIX IT when MPFI is patched */
    sollya_mpfi_neg(u, u);

    sollya_mpfi_inv(res[1], u);

    sollya_mpfi_set_ui(coeffs_array[0], 1);
  }

  for(index=2; index<=n; index++) {
    /* coeffs_array represents p_(index-1) */

    symbolic_poly_diff(coeffs_array_diff, coeffs_array, index-2);
    sollya_mpfi_set_ui(coeffs_array_diff[index-2], 0);
    /* now it represents p_(index-1)' */

    for(i=index-1; i>=2; i--) {
      sollya_mpfi_sub(coeffs_array[i], coeffs_array_diff[i], coeffs_array_diff[i-2]);
      sollya_mpfi_mul_ui(temp, coeffs_array[i-1], 2*(index-1));
      sollya_mpfi_add(coeffs_array[i], coeffs_array[i], temp);
      sollya_mpfi_div_ui(coeffs_array[i], coeffs_array[i], index);
    }

    sollya_mpfi_mul_ui(temp, coeffs_array[0], 2*(index-1));
    sollya_mpfi_add(coeffs_array[1], coeffs_array_diff[1], temp);
    sollya_mpfi_div_ui(coeffs_array[1], coeffs_array[1], index);

    sollya_mpfi_set(coeffs_array[0], coeffs_array_diff[0]);
    sollya_mpfi_div_ui(coeffs_array[0], coeffs_array[0], index);
    /* now it represents p_(index) */

    symbolic_poly_evaluation_horner(res[index], coeffs_array, x, index-1);
    sollya_mpfi_set_ui(temp, index);
    sollya_mpfi_pow(temp, u, temp);
    sollya_mpfi_div(res[index], res[index], temp);
  }

  for (index=0; index<=n-1; index++){
    sollya_mpfi_clear(coeffs_array[index]);
    sollya_mpfi_clear(coeffs_array_diff[index]);
  }
  sollya_mpfi_clear(u);
  sollya_mpfi_clear(temp);
  safeFree(coeffs_array);
  safeFree(coeffs_array_diff);

  return;
}



/* asin_diff : recurrence formula: p_(n+1)= (p'_n * (1-x^2) + (2n-1)x * p_n)/(n+1)
   asin^(0) = asin(x)
   asin^(n) / n! = p_(n)(x) / (1-x^2)^((2n-1)/2)
   p_1=1;

   --> degree of p_n is (n-1)
*/
void asin_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i,index;
  sollya_mpfi_t *coeffs_array, *coeffs_array_diff;
  sollya_mpfi_t u, temp;

  mp_prec_t prec;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  coeffs_array = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));
  coeffs_array_diff = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));

  for (index=0; index<=n-1; index++) {
    sollya_mpfi_init2(coeffs_array[index], prec);
    sollya_mpfi_init2(coeffs_array_diff[index], prec);

    sollya_mpfi_set_ui(coeffs_array[index], 0);
    sollya_mpfi_set_ui(coeffs_array_diff[index], 0);
  }

  sollya_mpfi_init2(u, prec);
  sollya_mpfi_init2(temp, prec);

  sollya_mpfi_asin(res[0], x);

  if(n>=1) {
    sollya_mpfi_sqr(u, x);
    sollya_mpfi_sub_ui(u, u, 1); /* TODO: FIX IT when MPFI is patched */
    sollya_mpfi_neg(u, u);
    sollya_mpfi_sqrt(u, u);

    sollya_mpfi_inv(res[1], u);

    sollya_mpfi_set_ui(coeffs_array[0], 1);
  }

  for(index=2; index<=n; index++) {
    /* coeffs_array represents p_(index-1) */

    symbolic_poly_diff(coeffs_array_diff, coeffs_array, index-2);
    sollya_mpfi_set_ui(coeffs_array_diff[index-2], 0);
    /* now it represents p_(index-1)' */

    for(i=index-1; i>=2; i--) {
      sollya_mpfi_sub(coeffs_array[i], coeffs_array_diff[i], coeffs_array_diff[i-2]);
      sollya_mpfi_mul_ui(temp, coeffs_array[i-1], 2*(index-1)-1);
      sollya_mpfi_add(coeffs_array[i], coeffs_array[i], temp);
      sollya_mpfi_div_ui(coeffs_array[i], coeffs_array[i], index);
    }

    sollya_mpfi_mul_ui(temp, coeffs_array[0], 2*(index-1)-1);
    sollya_mpfi_add(coeffs_array[1], coeffs_array_diff[1], temp);
    sollya_mpfi_div_ui(coeffs_array[1], coeffs_array[1], index);

    sollya_mpfi_set(coeffs_array[0], coeffs_array_diff[0]);
    sollya_mpfi_div_ui(coeffs_array[0], coeffs_array[0], index);
    /* now it represents p_(index) */

    symbolic_poly_evaluation_horner(res[index], coeffs_array, x, index-1);
    sollya_mpfi_set_ui(temp, 2*index-1);
    sollya_mpfi_pow(temp, u, temp);
    sollya_mpfi_div(res[index], res[index], temp);
  }

  for (index=0; index<=n-1; index++){
    sollya_mpfi_clear(coeffs_array[index]);
    sollya_mpfi_clear(coeffs_array_diff[index]);
  }
  sollya_mpfi_clear(u);
  sollya_mpfi_clear(temp);
  safeFree(coeffs_array);
  safeFree(coeffs_array_diff);

  return;
}


/* acos_diff : except for the res[0], all the terms are equal to -asin^(n)(x)/n! */
void acos_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i;

  asin_diff(res,x,n,silent);

  sollya_mpfi_acos(res[0], x);

  for (i=1; i<=n;i++)  sollya_mpfi_neg(res[i], res[i]);

  return;
}


/* asinh_diff : recurrence formula: p_(n+1) = (p'_n * (1+x^2) - (2n-1)x * p_n) / (n+1)
   asinh^(0) = asinh(x)
   asinh^(n)/n! = p_(n)(x) / (1+x^2)^((2n-1)/2)
   p_1=1;

   --> degree of p_n is (n-1)
*/
void asinh_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i,index;
  sollya_mpfi_t *coeffs_array, *coeffs_array_diff;
  sollya_mpfi_t u, temp;

  mp_prec_t prec;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  coeffs_array = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));
  coeffs_array_diff = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));

  for (index=0; index<=n-1; index++) {
    sollya_mpfi_init2(coeffs_array[index], prec);
    sollya_mpfi_init2(coeffs_array_diff[index], prec);

    sollya_mpfi_set_ui(coeffs_array[index], 0);
    sollya_mpfi_set_ui(coeffs_array_diff[index], 0);
  }

  sollya_mpfi_init2(u, prec);
  sollya_mpfi_init2(temp, prec);

  sollya_mpfi_asinh(res[0], x);

  if(n>=1) {
    sollya_mpfi_sqr(u, x);
    sollya_mpfi_add_ui(u, u, 1);
    sollya_mpfi_sqrt(u, u);

    sollya_mpfi_inv(res[1], u);

    sollya_mpfi_set_ui(coeffs_array[0], 1);
  }

  for(index=2; index<=n; index++) {
    /* coeffs_array represents p_(index-1) */

    symbolic_poly_diff(coeffs_array_diff, coeffs_array, index-2);
    sollya_mpfi_set_ui(coeffs_array_diff[index-2], 0);
    /* now it represents p_(index-1)' */

    for(i=index-1; i>=2; i--) {
      sollya_mpfi_add(coeffs_array[i], coeffs_array_diff[i], coeffs_array_diff[i-2]);
      sollya_mpfi_mul_ui(temp, coeffs_array[i-1], 2*(index-1)-1);
      sollya_mpfi_sub(coeffs_array[i], coeffs_array[i], temp);
      sollya_mpfi_div_ui(coeffs_array[i], coeffs_array[i], index);
    }

    sollya_mpfi_mul_ui(temp, coeffs_array[0], 2*(index-1)-1);
    sollya_mpfi_sub(coeffs_array[1], coeffs_array_diff[1], temp);
    sollya_mpfi_div_ui(coeffs_array[1], coeffs_array[1], index);

    sollya_mpfi_set(coeffs_array[0], coeffs_array_diff[0]);
    sollya_mpfi_div_ui(coeffs_array[0], coeffs_array[0], index);
    /* now it represents p_(index) */

    symbolic_poly_evaluation_horner(res[index], coeffs_array, x, index-1);
    sollya_mpfi_set_ui(temp, 2*index-1);
    sollya_mpfi_pow(temp, u, temp);
    sollya_mpfi_div(res[index], res[index], temp);
  }

  for (index=0; index<=n-1; index++){
    sollya_mpfi_clear(coeffs_array[index]);
    sollya_mpfi_clear(coeffs_array_diff[index]);
  }
  sollya_mpfi_clear(u);
  sollya_mpfi_clear(temp);
  safeFree(coeffs_array);
  safeFree(coeffs_array_diff);

  return;
}


/* acosh_diff : recurrence formula: p_(n+1) = (p'_n * (x^2-1) - (2n-1)x * p_n) / (n+1)
   acosh^(0) = acosh(x)
   acosh^(n)/n! = p_(n)(x) / (x^2-1)^((2n-1)/2)
   p_1=1;

   --> degree of p_n is (n-1)
*/
void acosh_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i,index;
  sollya_mpfi_t *coeffs_array, *coeffs_array_diff;
  sollya_mpfi_t u, temp;

  mp_prec_t prec;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  coeffs_array = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));
  coeffs_array_diff = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));

  for (index=0; index<=n-1; index++) {
    sollya_mpfi_init2(coeffs_array[index], prec);
    sollya_mpfi_init2(coeffs_array_diff[index], prec);

    sollya_mpfi_set_ui(coeffs_array[index], 0);
    sollya_mpfi_set_ui(coeffs_array_diff[index], 0);
  }

  sollya_mpfi_init2(u, prec);
  sollya_mpfi_init2(temp, prec);

  sollya_mpfi_acosh(res[0], x);

  if(n>=1) {
    sollya_mpfi_sqr(u, x);
    sollya_mpfi_sub_ui(u, u, 1);
    sollya_mpfi_sqrt(u, u);

    sollya_mpfi_inv(res[1], u);

    sollya_mpfi_set_ui(coeffs_array[0], 1);
  }

  for(index=2; index<=n; index++) {
    /* coeffs_array represents p_(index-1) */

    symbolic_poly_diff(coeffs_array_diff, coeffs_array, index-2);
    sollya_mpfi_set_ui(coeffs_array_diff[index-2], 0);
    /* now it represents p_(index-1)' */

    for(i=index-1; i>=2; i--) {
      sollya_mpfi_sub(coeffs_array[i], coeffs_array_diff[i-2], coeffs_array_diff[i]);
      sollya_mpfi_mul_ui(temp, coeffs_array[i-1], 2*(index-1)-1);
      sollya_mpfi_sub(coeffs_array[i], coeffs_array[i], temp);
      sollya_mpfi_div_ui(coeffs_array[i], coeffs_array[i], index);
    }

    sollya_mpfi_mul_ui(temp, coeffs_array[0], 2*(index-1)-1);
    sollya_mpfi_add(coeffs_array[1], temp, coeffs_array_diff[1]);
    sollya_mpfi_neg(coeffs_array[1], coeffs_array[1]);
    sollya_mpfi_div_ui(coeffs_array[1], coeffs_array[1], index);

    sollya_mpfi_neg(coeffs_array[0], coeffs_array_diff[0]);
    sollya_mpfi_div_ui(coeffs_array[0], coeffs_array[0], index);
    /* now it represents p_(index) */

    symbolic_poly_evaluation_horner(res[index], coeffs_array, x, index-1);
    sollya_mpfi_set_ui(temp, 2*index-1);
    sollya_mpfi_pow(temp, u, temp);
    sollya_mpfi_div(res[index], res[index], temp);
  }

  for (index=0; index<=n-1; index++){
    sollya_mpfi_clear(coeffs_array[index]);
    sollya_mpfi_clear(coeffs_array_diff[index]);
  }
  sollya_mpfi_clear(u);
  sollya_mpfi_clear(temp);
  safeFree(coeffs_array);
  safeFree(coeffs_array_diff);

  return;
}

/* erf^(n)(x)/n! = p_n(x)*e^(-x^2)             */
/* with p_1(x) = 2/sqrt(pi)                    */
/* and p_(n+1)(x) = (p_n'(x) - 2xp_n(x))/(n+1) */
/*  -> degree of p_n is n-1                    */
void erf_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i,index;
  sollya_mpfi_t *coeffs_array, *coeffs_array_diff;
  sollya_mpfi_t u, temp;

  mp_prec_t prec;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  coeffs_array = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));
  coeffs_array_diff = (sollya_mpfi_t *)safeCalloc( n,sizeof(sollya_mpfi_t));

  for (index=0; index<=n-1; index++) {
    sollya_mpfi_init2(coeffs_array[index], prec);
    sollya_mpfi_init2(coeffs_array_diff[index], prec);

    sollya_mpfi_set_ui(coeffs_array[index], 0);
    sollya_mpfi_set_ui(coeffs_array_diff[index], 0);
  }

  sollya_mpfi_init2(u, prec);
  sollya_mpfi_init2(temp, prec);

  sollya_mpfi_erf(res[0], x);

  if(n>=1) {
    sollya_mpfi_const_pi(temp);
    sollya_mpfi_sqrt(temp, temp);
    sollya_mpfi_ui_div(temp, 2, temp);

    sollya_mpfi_sqr(u, x);
    sollya_mpfi_neg(u, u);
    sollya_mpfi_exp(u, u);

    sollya_mpfi_mul(res[1], temp, u);

    sollya_mpfi_set(coeffs_array[0], temp);
  }

  for(index=2; index<=n; index++) {
    /* coeffs_array represents p_(index-1) */

    symbolic_poly_diff(coeffs_array_diff, coeffs_array, index-2);
    sollya_mpfi_set_ui(coeffs_array_diff[index-2], 0);
    /* now it represents p_(index-1)' */

    for(i=index-1; i>=1; i--) {
      sollya_mpfi_mul_ui(temp, coeffs_array[i-1], 2);
      sollya_mpfi_sub(coeffs_array[i], coeffs_array_diff[i], temp);
      sollya_mpfi_div_ui(coeffs_array[i], coeffs_array[i], index);
    }

    sollya_mpfi_set(coeffs_array[0], coeffs_array_diff[0]);
    sollya_mpfi_div_ui(coeffs_array[0], coeffs_array[0], index);
    /* now it represents p_(index) */

    symbolic_poly_evaluation_horner(res[index], coeffs_array, x, index-1);
    sollya_mpfi_mul(res[index], res[index], u);
  }

  for (index=0; index<=n-1; index++){
    sollya_mpfi_clear(coeffs_array[index]);
    sollya_mpfi_clear(coeffs_array_diff[index]);
  }
  sollya_mpfi_clear(u);
  sollya_mpfi_clear(temp);
  safeFree(coeffs_array);
  safeFree(coeffs_array_diff);

  return;
}

void erfc_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i;

  erf_diff(res, x, n, silent);

  sollya_mpfi_erfc(res[0], x);

  for (i=1; i<=n;i++)  sollya_mpfi_neg(res[i],res[i]);

  return;
}

void abs_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  int i;
  mpfr_t temp2;
  mp_prec_t prec;

  prec = getToolPrecision();

  sollya_mpfi_abs(res[0], x);
  if(n >= 1) {
    if (sollya_mpfi_has_zero(x))  sollya_mpfi_interv_si(res[1], -1, 1);
    else sollya_mpfi_set_si(res[1], sollya_mpfi_is_nonneg(x) ? 1 : (-1));
  }

  if(n >= 2) {
    mpfr_init2(temp2, prec);
    mpfr_set_nan(temp2);

    if (!(*silent)) {
      *silent = 1;
      printMessage(1, SOLLYA_MSG_ABS_NOT_TWICE_DIFFERENTIABLE, "Warning: the absolute value is not twice differentiable.\n");
      printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
    }
    for(i=2;i<=n;i++) sollya_mpfi_set_fr(res[i], temp2);
    mpfr_clear(temp2);
  }

  return;
}

void single_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
  int i;
  mpfr_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  mpfr_init2(temp, prec);
  mpfr_set_nan(temp);

  if (!(*silent)) {
    *silent = 1;
    printMessage(1, SOLLYA_MSG_SINGLE_NOT_DIFFERENTIABLE, "Warning: the single rounding operator is not differentiable.\n");
    printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
  }
  for(i=0;i<=n;i++) sollya_mpfi_set_fr(res[i], temp);
  mpfr_clear(temp);
}

void quad_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
  int i;
  mpfr_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  mpfr_init2(temp, prec);
  mpfr_set_nan(temp);

  if (!(*silent)) {
    *silent = 1;
    printMessage(1, SOLLYA_MSG_QUAD_NOT_DIFFERENTIABLE, "Warning: the quad rounding operator is not differentiable.\n");
    printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
  }
  for(i=0;i<=n;i++) sollya_mpfi_set_fr(res[i], temp);
  mpfr_clear(temp);
}

void halfprecision_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
  int i;
  mpfr_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  mpfr_init2(temp, prec);
  mpfr_set_nan(temp);

  if (!(*silent)) {
    *silent = 1;
    printMessage(1, SOLLYA_MSG_HALF_NOT_DIFFERENTIABLE, "Warning: the half-precision rounding operator is not differentiable.\n");
    printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
  }
  for(i=0;i<=n;i++) sollya_mpfi_set_fr(res[i], temp);
  mpfr_clear(temp);
}

void double_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
  int i;
  mpfr_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  mpfr_init2(temp, prec);
  mpfr_set_nan(temp);

  if (!(*silent)) {
    *silent = 1;
    printMessage(1, SOLLYA_MSG_DOUBLE_NOT_DIFFERENTIABLE, "Warning: the double rounding operator is not differentiable.\n");
    printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
  }
  for(i=0;i<=n;i++) sollya_mpfi_set_fr(res[i], temp);
  mpfr_clear(temp);
}

void doubledouble_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
  int i;
  mpfr_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  mpfr_init2(temp, prec);
  mpfr_set_nan(temp);

  if (!(*silent)) {
    *silent = 1;
    printMessage(1, SOLLYA_MSG_DOUBLE_DOUBLE_NOT_DIFFERENTIABLE, "Warning: the doubledouble rounding operator is not differentiable.\n");
    printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
  }
  for(i=0;i<=n;i++) sollya_mpfi_set_fr(res[i], temp);
  mpfr_clear(temp);
}

void tripledouble_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
  int i;
  mpfr_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  mpfr_init2(temp, prec);
  mpfr_set_nan(temp);

  if (!(*silent)) {
    *silent = 1;
    printMessage(1, SOLLYA_MSG_TRIPLE_DOUBLE_NOT_DIFFERENTIABLE, "Warning: the tripledouble rounding operator is not differentiable.\n");
    printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
  }
  for(i=0;i<=n;i++) sollya_mpfi_set_fr(res[i], temp);
  mpfr_clear(temp);
}

void doubleextended_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
  int i;
  mpfr_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  mpfr_init2(temp, prec);
  mpfr_set_nan(temp);

  if (!(*silent)) {
    *silent = 1;
    printMessage(1, SOLLYA_MSG_DOUBLEEXTENDED_NOT_DIFFERENTIABLE, "Warning: the doubleextended rounding operator is not differentiable.\n");
    printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
  }
  for(i=0;i<=n;i++) sollya_mpfi_set_fr(res[i], temp);
  mpfr_clear(temp);
}

void ceil_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
  int i;
  mpfr_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  mpfr_init2(temp, prec);
  mpfr_set_nan(temp);

  if (!(*silent)) {
    *silent = 1;
    printMessage(1, SOLLYA_MSG_CEIL_NOT_DIFFERENTIABLE, "Warning: the ceil operator is not differentiable.\n");
    printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
  }
  for(i=0;i<=n;i++) sollya_mpfi_set_fr(res[i], temp);
  mpfr_clear(temp);
}

void floor_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
  int i;
  mpfr_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  mpfr_init2(temp, prec);
  mpfr_set_nan(temp);

  if (!(*silent)) {
    *silent = 1;
    printMessage(1, SOLLYA_MSG_FLOOR_NOT_DIFFERENTIABLE, "Warning: the floor operator is not differentiable.\n");
    printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
  }
  for(i=0;i<=n;i++) sollya_mpfi_set_fr(res[i], temp);
  mpfr_clear(temp);
}

void nearestint_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
  int i;
  mpfr_t temp;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  mpfr_init2(temp, prec);
  mpfr_set_nan(temp);

  if (!(*silent)) {
    *silent = 1;
    printMessage(1, SOLLYA_MSG_NEARESTINT_NOT_DIFFERENTIABLE, "Warning: the nearestint operator is not differentiable.\n");
    printMessage(1, SOLLYA_MSG_CONTINUATION, "Will return [@NaN@, @NaN@].\n");
  }
  for(i=0;i<=n;i++) sollya_mpfi_set_fr(res[i], temp);
  mpfr_clear(temp);
}






/******************************************************************************/
/*                                                                            */
/*                    EXPRESSION DIFFERENTIATION FUNCTIONS                    */
/*                                                                            */
/******************************************************************************/

node *sqrt_diff_expr(node *g) {
  return makeDiv( differentiateUnsimplified(g),
                  makeMul(makeConstantDouble(2.0),
                          makeSqrt(copyTree(g))
                          )
                  );
}

node *exp_diff_expr(node *g) {
  return makeMul( makeExp(copyTree(g)),
                  differentiateUnsimplified(g)
                  );
}

node *log_diff_expr(node *g) {
  return makeMul( makeDiv(makeConstantDouble(1.0), copyTree(g)),
                  differentiateUnsimplified(g)
                  );
}

node *log2_diff_expr(node *g) {
  return makeMul( makeDiv( makeConstantDouble(1.0),
                           makeMul(copyTree(g), makeLog(makeConstantDouble(2.0)))
                           ),
                  differentiateUnsimplified(g)
                  );
}

node *log10_diff_expr(node *g) {
  return makeMul( makeDiv( makeConstantDouble(1.0),
                           makeMul(copyTree(g), makeLog(makeConstantDouble(10.0)))
                           ),
                  differentiateUnsimplified(g)
                  );
}

node *sin_diff_expr(node *g) {
  return makeMul( makeCos(copyTree(g)),
                  differentiateUnsimplified(g)
                  );
}

node *cos_diff_expr(node *g) {
  return makeMul( makeNeg(makeSin(copyTree(g))),
                  differentiateUnsimplified(g)
                  );
}

node *tan_diff_expr(node *g) {
  return makeMul( makeAdd(makeConstantDouble(1.0),
                          makePow(makeTan(copyTree(g)), makeConstantDouble(2.0))
                          ),
                  differentiateUnsimplified(g)
                  );
}

node *asin_diff_expr(node *g) {
  return makeMul( makeDiv(makeConstantDouble(1.0),
                          makeSqrt(makeSub(makeConstantDouble(1.0),
                                           makePow(copyTree(g), makeConstantDouble(2.0))
                                           )
                                   )
                          ),
                  differentiateUnsimplified(g)
                  );
}

node *acos_diff_expr(node *g) {
  return makeMul( makeDiv(makeConstantDouble(-1.0),
                          makeSqrt(makeSub(makeConstantDouble(1.0),
                                           makePow(copyTree(g), makeConstantDouble(2.0))
                                           )
                                   )
                          ),
                  differentiateUnsimplified(g)
                  );
}

node *atan_diff_expr(node *g) {
  return makeMul( makeDiv(makeConstantDouble(1.0),
                          makeAdd(makeConstantDouble(1.0),
                                  makeMul(copyTree(g), copyTree(g))
                                  )
                          ),
                  differentiateUnsimplified(g)
                  );
}

node *sinh_diff_expr(node *g) {
  return makeMul( makeCosh(copyTree(g)),
                  differentiateUnsimplified(g)
                  );
}

node *cosh_diff_expr(node *g) {
  return makeMul( makeSinh(copyTree(g)),
                  differentiateUnsimplified(g)
                  );
}

node *tanh_diff_expr(node *g) {
  return makeMul( makeSub(makeConstantDouble(1.0),
                          makePow(makeTanh(copyTree(g)), makeConstantDouble(2.0))
                          ),
                  differentiateUnsimplified(g)
                  );
}

node *asinh_diff_expr(node *g) {
  return makeMul( makeDiv(makeConstantDouble(1.0),
                          makeSqrt(makeAdd(makeConstantDouble(1.0),
                                           makeMul(copyTree(g), copyTree(g))
                                           )
                                   )
                          ),
                  differentiateUnsimplified(g)
                  );
}

node *acosh_diff_expr(node *g) {
  return makeMul( makeDiv(makeConstantDouble(1.0),
                          makeSqrt(makeAdd(makeConstantDouble(-1.0),
                                           makeMul(copyTree(g), copyTree(g))
                                           )
                                   )
                          ),
                  differentiateUnsimplified(g)
                  );
}

node *atanh_diff_expr(node *g) {
  return makeMul( makeDiv(makeConstantDouble(1.0),
                          makeSub(makeConstantDouble(1.0),
                                  makeMul(copyTree(g), copyTree(g))
                                  )
                          ),
                  differentiateUnsimplified(g)
                  );
}

node *erf_diff_expr(node *g) {
  return makeMul( makeDiv(makeExp(makeNeg(makePow(copyTree(g), makeConstantDouble(2.0)))),
                          makeSqrt(makeAtan(makeConstantDouble(1.0)))
                          ),
                  differentiateUnsimplified(g)
                  );
}

node *erfc_diff_expr(node *g) {
  return makeMul( makeDiv(makeNeg(makeExp(makeNeg(makePow(copyTree(g), makeConstantDouble(2.0))))),
                          makeSqrt(makeAtan(makeConstantDouble(1.0)))
                          ),
                  differentiateUnsimplified(g)
                  );
}

node *log1p_diff_expr(node *g) {
  return makeMul( makeDiv(makeConstantDouble(1.0),
                          makeAdd(makeConstantDouble(1.0),
                                  copyTree(g)
                                  )
                          ),
                  differentiateUnsimplified(g)
                  );
}

node *expm1_diff_expr(node *g) {
  return makeMul( makeExp(copyTree(g)),
                  differentiateUnsimplified(g)
                  );
}

node *abs_diff_expr(node *g) {
  return makeMul( makeDiv(copyTree(g),
                          makeAbs(copyTree(g))
                          ),
                  differentiateUnsimplified(g)
                  );
}

node *double_diff_expr(node *g) {
  UNUSED_PARAM(g);
  printMessage(1,SOLLYA_MSG_DOUBLE_NOT_DIFFERENTIABLE,
               "Warning: the double rounding operator is not differentiable.\nReplacing it by a constant function when differentiating.\n");
  return makeConstantDouble(0.0);
}

node *single_diff_expr(node *g) {
  UNUSED_PARAM(g);
  printMessage(1,SOLLYA_MSG_SINGLE_NOT_DIFFERENTIABLE,
               "Warning: the single rounding operator is not differentiable.\nReplacing it by a constant function when differentiating.\n");
  return makeConstantDouble(0.0);
}

node *quad_diff_expr(node *g) {
  UNUSED_PARAM(g);
  printMessage(1,SOLLYA_MSG_QUAD_NOT_DIFFERENTIABLE,
               "Warning: the quad rounding operator is not differentiable.\nReplacing it by a constant function when differentiating.\n");
  return makeConstantDouble(0.0);
}

node *halfprecision_diff_expr(node *g) {
  UNUSED_PARAM(g);
  printMessage(1,SOLLYA_MSG_HALF_NOT_DIFFERENTIABLE,
               "Warning: the half-precision rounding operator is not differentiable.\nReplacing it by a constant function when differentiating.\n");
  return makeConstantDouble(0.0);
}

node *doubledouble_diff_expr(node *g) {
  UNUSED_PARAM(g);
  printMessage(1,SOLLYA_MSG_DOUBLE_DOUBLE_NOT_DIFFERENTIABLE,
               "Warning: the double-double rounding operator is not differentiable.\nReplacing it by a constant function when differentiating.\n");
  return makeConstantDouble(0.0);
}

node *tripledouble_diff_expr(node *g) {
  UNUSED_PARAM(g);
  printMessage(1,SOLLYA_MSG_TRIPLE_DOUBLE_NOT_DIFFERENTIABLE,
               "Warning: the triple-double rounding operator is not differentiable.\nReplacing it by a constant function when differentiating.\n");
  return makeConstantDouble(0.0);
}

node *doubleextended_diff_expr(node *g) {
  UNUSED_PARAM(g);
  printMessage(1,SOLLYA_MSG_DOUBLEEXTENDED_NOT_DIFFERENTIABLE,
               "Warning: the double-extended rounding operator is not differentiable.\nReplacing it by a constant function when differentiating.\n");
  return makeConstantDouble(0.0);
}

node *ceil_diff_expr(node *g) {
  UNUSED_PARAM(g);	printMessage(1,SOLLYA_MSG_CEIL_NOT_DIFFERENTIABLE,
                                     "Warning: the ceil operator is not differentiable.\nReplacing it by a constant function when differentiating.\n");
  return makeConstantDouble(0.0);
}

node *floor_diff_expr(node *g) {
  UNUSED_PARAM(g);
  printMessage(1,SOLLYA_MSG_FLOOR_NOT_DIFFERENTIABLE,
               "Warning: the floor operator is not differentiable.\nReplacing it by a constant function when differentiating.\n");
  return makeConstantDouble(0.0);
}

node *nearestint_diff_expr(node *g) {
  UNUSED_PARAM(g);
  printMessage(1,SOLLYA_MSG_NEARESTINT_NOT_DIFFERENTIABLE,
               "Warning: the nearestint operator is not differentiable.\nReplacing it by a constant function when differentiating.\n");
  return makeConstantDouble(0.0);
}






/******************************************************************************/
/*                                                                            */
/*                   FUNCTIONS FOR SIMPLIFYING EXPRESSIONS                    */
/*                                                                            */
/******************************************************************************/

/* Generic algorithm for simplifying an expression of the form f(g) where g is already simplified */
node *simplify_generic(baseFunction *f, node *g) {
  node *simplified;
  mpfr_t *value;

  simplified = (node*) safeMalloc(sizeof(node));

  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value, tools_precision);
    if ( (f->point_eval(*value, *(accessThruMemRef(g)->value), GMP_RNDN) == 0) &&
         mpfr_number_p(*value) )  {
      simplified->nodeType = CONSTANT;
      simplified->value = value;
      free_memory(g);
      return simplified;
    }
    else {
      mpfr_clear(*value);
      safeFree(value);
    }
  }

  /* Default case: we do nothing but constructing f(g) */
  simplified->nodeType = UNARY_BASE_FUNC;
  simplified->baseFun = f;
  simplified->child1 = g;
  return simplified;
}

/* Algorithm for simplifying an expression of the form f(g) where f is nearestint|floor|ceil and g is already simplified */
node *simplify_integer_rounding_functions(baseFunction *f, node *g) {
  node *simplified;
  mpfr_t *value;
  rangetype xrange, yrange;
  mp_prec_t prec;
  int ok = 0;

  simplified = (node*) safeMalloc(sizeof(node));

  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    prec = mpfr_get_prec(*(accessThruMemRef(g)->value));
    prec = (tools_precision<prec)?prec:tools_precision;
    mpfr_init2(*value, prec);
    if ( (f->point_eval(*value, *(accessThruMemRef(g)->value), GMP_RNDN) == 0) &&
         mpfr_number_p(*value) )  {
      simplified->nodeType = CONSTANT;
      simplified->value = value;
      free_memory(g);
      return simplified;
    }
    else {
      mpfr_clear(*value);
      safeFree(value);
    }
  }
  else if (isConstant(g)) {
    xrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    xrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    yrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    yrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*(yrange.a),tools_precision);
    mpfr_init2(*(yrange.b),tools_precision);
    mpfr_init2(*(xrange.a),4 * tools_precision);
    mpfr_init2(*(xrange.b),4 * tools_precision);
    mpfr_set_ui(*(xrange.a),1,GMP_RNDD);
    mpfr_set_ui(*(xrange.b),1,GMP_RNDU);
    evaluateRangeFunction(yrange, g, xrange, 8 * tools_precision);
    f->point_eval(*(xrange.a),*(yrange.a),GMP_RNDD);
    f->point_eval(*(xrange.b),*(yrange.b),GMP_RNDU);
    if ( mpfr_number_p(*(xrange.a)) &&
         mpfr_number_p(*(xrange.b)) &&
         (mpfr_cmp(*(xrange.a),*(xrange.b)) == 0) ) {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      mpfr_set(*value,*(xrange.a),GMP_RNDN); /* Exact. TODO: why should it be exact? */
      free_memory(g);
      ok = 1;
    }
    mpfr_clear(*(xrange.a));
    mpfr_clear(*(xrange.b));
    mpfr_clear(*(yrange.a));
    mpfr_clear(*(yrange.b));
    safeFree(xrange.a);
    safeFree(xrange.b);
    safeFree(yrange.a);
    safeFree(yrange.b);
    if (ok) return simplified;
  }

  /* Default case: we do nothing but constructing f(g) */
  simplified->nodeType = UNARY_BASE_FUNC;
  simplified->baseFun = f;
  simplified->child1 = g;
  return simplified;
}

node *simplify_exp(node *g) { return simplify_generic(basefun_exp, g); }
node *simplify_log2(node *g) { return simplify_generic(basefun_log2, g); }
node *simplify_log10(node *g) { return simplify_generic(basefun_log10, g); }
node *simplify_sin(node *g) { return simplify_generic(basefun_sin, g); }
node *simplify_cos(node *g) { return simplify_generic(basefun_cos, g); }
node *simplify_tan(node *g) { return simplify_generic(basefun_tan, g); }
node *simplify_sinh(node *g) { return simplify_generic(basefun_sinh, g); }
node *simplify_cosh(node *g) { return simplify_generic(basefun_cosh, g); }
node *simplify_tanh(node *g) { return simplify_generic(basefun_tanh, g); }
node *simplify_asinh(node *g) { return simplify_generic(basefun_asinh, g); }
node *simplify_acosh(node *g) { return simplify_generic(basefun_acosh, g); }
node *simplify_atanh(node *g) { return simplify_generic(basefun_atanh, g); }
node *simplify_erf(node *g) { return simplify_generic(basefun_erf, g); }
node *simplify_erfc(node *g) { return simplify_generic(basefun_erfc, g); }
node *simplify_log1p(node *g) { return simplify_generic(basefun_log1p, g); }
node *simplify_expm1(node *g) { return simplify_generic(basefun_expm1, g); }

node *simplify_sqrt(node *g) {
  node *simplified;
  mpfr_t *value;
  mp_prec_t prec, p;

  simplified = (node*) safeMalloc(sizeof(node));
  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    simplified->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    prec = tools_precision;
    p = mpfr_get_prec(*(accessThruMemRef(g)->value));
    if (p > prec) prec = p;
    prec += 10;
    if (prec > 256 * tools_precision) prec = 256 * tools_precision;
    mpfr_init2(*value,prec);
    simplified->value = value;
    if ((mpfr_sqrt(*value, *(accessThruMemRef(g)->value), GMP_RNDN) != 0) ||
        (!mpfr_number_p(*value))) {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_sqrt;
      simplified->child1 = g;
      mpfr_clear(*value);
      safeFree(value);
    } else {
      free_memory(g);
    }
  } else {
    if ((accessThruMemRef(g)->nodeType == POW) &&
        (accessThruMemRef(accessThruMemRef(g)->child2)->nodeType == CONSTANT) &&
        (mpfr_cmp_d(*(accessThruMemRef(accessThruMemRef(g)->child2)->value),2.0) == 0.0) &&
        (!mpfr_nan_p(*(accessThruMemRef(accessThruMemRef(g)->child2)->value)))) {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_abs;
      simplified->child1 = copyTree(accessThruMemRef(g)->child1);
      free_memory(g);
    } else {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_sqrt;
      simplified->child1 = g;
    }
  }
  return simplified;
}


node *simplify_log(node *g) {
  node *simplified;
  mpfr_t *value;

  simplified = (node*) safeMalloc(sizeof(node));
  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    simplified->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value,tools_precision);
    simplified->value = value;
    if ((mpfr_log(*value, *(accessThruMemRef(g)->value), GMP_RNDN) != 0) ||
        (!mpfr_number_p(*value))) {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_log;
      simplified->child1 = g;
      mpfr_clear(*value);
      safeFree(value);
    } else {
      free_memory(g);
    }
  } else {
    if ((accessThruMemRef(g)->nodeType == ADD) &&
        (accessThruMemRef(accessThruMemRef(g)->child1)->nodeType == CONSTANT) &&
        (mpfr_cmp_d(*(accessThruMemRef(accessThruMemRef(g)->child1)->value),1.0) == 0.0) &&
        (!mpfr_nan_p(*(accessThruMemRef(accessThruMemRef(g)->child1)->value)))) {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_log1p;
      simplified->child1 = copyTree(accessThruMemRef(g)->child2);
      free_memory(g);
    } else {
      if ((accessThruMemRef(g)->nodeType == ADD) &&
          (accessThruMemRef(accessThruMemRef(g)->child2)->nodeType == CONSTANT) &&
          (mpfr_cmp_d(*(accessThruMemRef(accessThruMemRef(g)->child2)->value),1.0) == 0.0) &&
          (!mpfr_nan_p(*(accessThruMemRef(accessThruMemRef(g)->child2)->value)))) {
        simplified->nodeType = UNARY_BASE_FUNC;
        simplified->baseFun = basefun_log1p;
        simplified->child1 = copyTree(accessThruMemRef(g)->child1);
        free_memory(g);
      } else {
        if ( (accessThruMemRef(g)->nodeType == UNARY_BASE_FUNC) &&
             (accessThruMemRef(g)->baseFun->baseFunctionCode == EXP) ) {
          safeFree(simplified);
          simplified = copyTree(accessThruMemRef(g)->child1);
          free_memory(g);
        } else {
          simplified->nodeType = UNARY_BASE_FUNC;
          simplified->baseFun = basefun_log;
          simplified->child1 = g;
        }
      }
    }
  }
  return simplified;
}




node *simplify_asin(node *g) {
  node *simplified;
  mpfr_t *value;
  mpfr_t temp;

  simplified = (node*) safeMalloc(sizeof(node));
  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    mpfr_init2(temp,53);
    mpfr_set_ui(temp,1,GMP_RNDN);
    if ((mpfr_cmp(temp, *(accessThruMemRef(g)->value)) == 0) && (!mpfr_nan_p(*(accessThruMemRef(g)->value)))) {
      free_memory(g);
      simplified->child2 = (node *) safeMalloc(sizeof(node));
      simplified->child2->nodeType = CONSTANT;
      simplified->child2->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*(simplified->child2->value),53);
      mpfr_set_si(*(simplified->child2->value),2,GMP_RNDN);
      simplified->nodeType = DIV;
      simplified->child1 = (node *) safeMalloc(sizeof(node));
      simplified->child1->nodeType = PI_CONST;
    } else {
      mpfr_set_si(temp,-1,GMP_RNDN);
      if ((mpfr_cmp(temp, *(accessThruMemRef(g)->value)) == 0) && (!mpfr_nan_p(*(accessThruMemRef(g)->value)))) {
        free_memory(g);
        simplified->child2 = (node *) safeMalloc(sizeof(node));
        simplified->child2->nodeType = CONSTANT;
        simplified->child2->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*(simplified->child2->value),53);
        mpfr_set_si(*(simplified->child2->value),-2,GMP_RNDN);
        simplified->nodeType = DIV;
        simplified->child1 = (node *) safeMalloc(sizeof(node));
        simplified->child1->nodeType = PI_CONST;
      } else {
        simplified->nodeType = CONSTANT;
        value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*value,tools_precision);
        simplified->value = value;
        if ((mpfr_asin(*value, *(accessThruMemRef(g)->value), GMP_RNDN) != 0) ||
            (!mpfr_number_p(*value))) {
          simplified->nodeType = UNARY_BASE_FUNC;
          simplified->baseFun = basefun_asin;
          simplified->child1 = g;
          mpfr_clear(*value);
          safeFree(value);
        } else {
          free_memory(g);
        }
      }
    }
    mpfr_clear(temp);
  } else {
    simplified->nodeType = UNARY_BASE_FUNC;
    simplified->baseFun = basefun_asin;
    simplified->child1 = g;
  }
  return simplified;
}


node *simplify_acos(node *g) {
  node *simplified;
  mpfr_t *value;
  mpfr_t temp;

  simplified = (node*) safeMalloc(sizeof(node));
  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    mpfr_init2(temp,53);
    mpfr_set_si(temp,-1,GMP_RNDN);
    if ((mpfr_cmp(temp, *(accessThruMemRef(g)->value)) == 0) && (!mpfr_nan_p(*(accessThruMemRef(g)->value)))) {
      free_memory(g);
      simplified->nodeType = PI_CONST;
    } else {
      simplified->nodeType = CONSTANT;
      value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*value,tools_precision);
      simplified->value = value;
      if ((mpfr_acos(*value, *(accessThruMemRef(g)->value), GMP_RNDN) != 0) ||
          (!mpfr_number_p(*value))) {
        simplified->nodeType = UNARY_BASE_FUNC;
        simplified->baseFun = basefun_acos;
        simplified->child1 = g;
        mpfr_clear(*value);
        safeFree(value);
      } else {
        free_memory(g);
      }
    }
    mpfr_clear(temp);
  } else {
    simplified->nodeType = UNARY_BASE_FUNC;
    simplified->baseFun = basefun_acos;
    simplified->child1 = g;
  }
  return simplified;
}


node *simplify_atan(node *g) {
  node *simplified;
  mpfr_t *value;
  mpfr_t temp;

  simplified = (node*) safeMalloc(sizeof(node));
  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    mpfr_init2(temp,53);
    mpfr_set_ui(temp,1,GMP_RNDN);
    if ((mpfr_cmp(temp, *(accessThruMemRef(g)->value)) == 0) && (!mpfr_nan_p(*(accessThruMemRef(g)->value)))) {
      free_memory(g);
      simplified->child2 = (node *) safeMalloc(sizeof(node));
      simplified->child2->nodeType = CONSTANT;
      simplified->child2->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*(simplified->child2->value),53);
      mpfr_set_si(*(simplified->child2->value),4,GMP_RNDN);
      simplified->nodeType = DIV;
      simplified->child1 = (node *) safeMalloc(sizeof(node));
      simplified->child1->nodeType = PI_CONST;
    } else {
      mpfr_set_si(temp,-1,GMP_RNDN);
      if ((mpfr_cmp(temp, *(accessThruMemRef(g)->value)) == 0) && (!mpfr_nan_p(*(accessThruMemRef(g)->value)))) {
        free_memory(g);
        simplified->child2 = (node *) safeMalloc(sizeof(node));
        simplified->child2->nodeType = CONSTANT;
        simplified->child2->value = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*(simplified->child2->value),53);
        mpfr_set_si(*(simplified->child2->value),-4,GMP_RNDN);
        simplified->nodeType = DIV;
        simplified->child1 = (node *) safeMalloc(sizeof(node));
        simplified->child1->nodeType = PI_CONST;
      } else {
        simplified->nodeType = CONSTANT;
        value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*value,tools_precision);
        simplified->value = value;
        if ((mpfr_atan(*value, *(accessThruMemRef(g)->value), GMP_RNDN) != 0) ||
            (!mpfr_number_p(*value))) {
          simplified->nodeType = UNARY_BASE_FUNC;
          simplified->baseFun = basefun_atan;
          simplified->child1 = g;
          mpfr_clear(*value);
          safeFree(value);
        } else {
          free_memory(g);
        }
      }
    }
    mpfr_clear(temp);
  } else {
    simplified->nodeType = UNARY_BASE_FUNC;
    simplified->baseFun = basefun_atan;
    simplified->child1 = g;
  }
  return simplified;
}



node *simplify_abs(node *g) {
  node *simplified;
  mpfr_t *value;

  simplified = (node*) safeMalloc(sizeof(node));

  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value, tools_precision);
    if ( (mpfr_abs(*value, *(accessThruMemRef(g)->value), GMP_RNDN) == 0) &&
         mpfr_number_p(*value) )  {
      simplified->nodeType = CONSTANT;
      simplified->value = value;
      free_memory(g);
      return simplified;
    }
    else {
      mpfr_clear(*value);
      safeFree(value);
    }
  }

  /* abs(abs(something)) --> abs(something) */
  if ( (accessThruMemRef(g)->nodeType == UNARY_BASE_FUNC) &&
       (accessThruMemRef(g)->baseFun->baseFunctionCode == ABS) ) {
    simplified->nodeType = UNARY_BASE_FUNC;
    simplified->baseFun = basefun_abs;
    simplified->child1 = copyTree(accessThruMemRef(g)->child1);
    free_memory(g);
    return simplified;
  }

  /* Default case: we do nothing but constructing f(g) */
  simplified->nodeType = UNARY_BASE_FUNC;
  simplified->baseFun = basefun_abs;
  simplified->child1 = g;
  return simplified;
}


node *simplify_double(node *g) {
  node *simplified;
  mpfr_t *value;
  rangetype xrange, yrange;

  simplified = (node*) safeMalloc(sizeof(node));
  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    simplified->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value,((tools_precision > 64)? tools_precision : 64));
    simplified->value = value;
    mpfr_round_to_double(*value, *(accessThruMemRef(g)->value));
    if (!mpfr_number_p(*value)) {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_double;
      simplified->child1 = g;
      mpfr_clear(*value);
      safeFree(value);
    } else {
      free_memory(g);
    }
  } else {
    if (isConstant(g)) {
      xrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      xrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      yrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      yrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*(yrange.a),4* tools_precision);
      mpfr_init2(*(yrange.b),4 * tools_precision);
      mpfr_init2(*(xrange.a),((tools_precision > 64)? tools_precision : 64));
      mpfr_init2(*(xrange.b),((tools_precision > 64)? tools_precision : 64));
      mpfr_set_ui(*(xrange.a),1,GMP_RNDD);
      mpfr_set_ui(*(xrange.b),1,GMP_RNDU);
      evaluateRangeFunction(yrange, g, xrange, 8 * tools_precision);
      mpfr_round_to_double(*(xrange.a),*(yrange.a));
      mpfr_round_to_double(*(xrange.b),*(yrange.b));
      if (mpfr_number_p(*(xrange.a)) &&
          mpfr_number_p(*(xrange.b)) &&
          (mpfr_cmp(*(xrange.a),*(xrange.b)) == 0)) {
        simplified->nodeType = CONSTANT;
        value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*value,((tools_precision > 64)? tools_precision : 64));
        simplified->value = value;
        mpfr_set(*value,*(xrange.a),GMP_RNDN); /* Exact */
        free_memory(g);
      } else {
        simplified->nodeType = UNARY_BASE_FUNC;
        simplified->baseFun = basefun_double;
        simplified->child1 = g;
      }
      mpfr_clear(*(xrange.a));
      mpfr_clear(*(xrange.b));
      mpfr_clear(*(yrange.a));
      mpfr_clear(*(yrange.b));
      safeFree(xrange.a);
      safeFree(xrange.b);
      safeFree(yrange.a);
      safeFree(yrange.b);
    } else {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_double;
      simplified->child1 = g;
    }
  }
  return simplified;
}


node *simplify_single(node *g) {
  node *simplified;
  mpfr_t *value;
  rangetype xrange, yrange;

  simplified = (node*) safeMalloc(sizeof(node));
  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    simplified->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value,((tools_precision > 64)? tools_precision : 64));
    simplified->value = value;
    mpfr_round_to_single(*value, *(accessThruMemRef(g)->value));
    if (!mpfr_number_p(*value)) {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_single;
      simplified->child1 = g;
      mpfr_clear(*value);
      safeFree(value);
    } else {
      free_memory(g);
    }
  } else {
    if (isConstant(g)) {
      xrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      xrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      yrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      yrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*(yrange.a),4* tools_precision);
      mpfr_init2(*(yrange.b),4 * tools_precision);
      mpfr_init2(*(xrange.a),((tools_precision > 64)? tools_precision : 64));
      mpfr_init2(*(xrange.b),((tools_precision > 64)? tools_precision : 64));
      mpfr_set_ui(*(xrange.a),1,GMP_RNDD);
      mpfr_set_ui(*(xrange.b),1,GMP_RNDU);
      evaluateRangeFunction(yrange, g, xrange, 8 * tools_precision);
      mpfr_round_to_single(*(xrange.a),*(yrange.a));
      mpfr_round_to_single(*(xrange.b),*(yrange.b));
      if (mpfr_number_p(*(xrange.a)) &&
          mpfr_number_p(*(xrange.b)) &&
          (mpfr_cmp(*(xrange.a),*(xrange.b)) == 0)) {
        simplified->nodeType = CONSTANT;
        value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*value,((tools_precision > 64)? tools_precision : 64));
        simplified->value = value;
        mpfr_set(*value,*(xrange.a),GMP_RNDN); /* Exact */
        free_memory(g);
      } else {
        simplified->nodeType = UNARY_BASE_FUNC;
        simplified->baseFun = basefun_single;
        simplified->child1 = g;
      }
      mpfr_clear(*(xrange.a));
      mpfr_clear(*(xrange.b));
      mpfr_clear(*(yrange.a));
      mpfr_clear(*(yrange.b));
      safeFree(xrange.a);
      safeFree(xrange.b);
      safeFree(yrange.a);
      safeFree(yrange.b);
    } else {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_single;
      simplified->child1 = g;
    }
  }
  return simplified;
}


node *simplify_quad(node *g) {
  node *simplified;
  mpfr_t *value;
  rangetype xrange, yrange;

  simplified = (node*) safeMalloc(sizeof(node));
  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    simplified->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value,((tools_precision > 128)? tools_precision : 128));
    simplified->value = value;
    mpfr_round_to_quad(*value, *(accessThruMemRef(g)->value));
    if (!mpfr_number_p(*value)) {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_quad;
      simplified->child1 = g;
      mpfr_clear(*value);
      safeFree(value);
    } else {
      free_memory(g);
    }
  } else {
    if (isConstant(g)) {
      xrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      xrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      yrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      yrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*(yrange.a),4* tools_precision);
      mpfr_init2(*(yrange.b),4 * tools_precision);
      mpfr_init2(*(xrange.a),((tools_precision > 128)? tools_precision : 128));
      mpfr_init2(*(xrange.b),((tools_precision > 128)? tools_precision : 128));
      mpfr_set_ui(*(xrange.a),1,GMP_RNDD);
      mpfr_set_ui(*(xrange.b),1,GMP_RNDU);
      evaluateRangeFunction(yrange, g, xrange, 8 * tools_precision);
      mpfr_round_to_quad(*(xrange.a),*(yrange.a));
      mpfr_round_to_quad(*(xrange.b),*(yrange.b));
      if (mpfr_number_p(*(xrange.a)) &&
          mpfr_number_p(*(xrange.b)) &&
          (mpfr_cmp(*(xrange.a),*(xrange.b)) == 0)) {
        simplified->nodeType = CONSTANT;
        value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*value,((tools_precision > 128)? tools_precision : 128));
        simplified->value = value;
        mpfr_set(*value,*(xrange.a),GMP_RNDN); /* Exact */
        free_memory(g);
      } else {
        simplified->nodeType = UNARY_BASE_FUNC;
        simplified->baseFun = basefun_quad;
        simplified->child1 = g;
      }
      mpfr_clear(*(xrange.a));
      mpfr_clear(*(xrange.b));
      mpfr_clear(*(yrange.a));
      mpfr_clear(*(yrange.b));
      safeFree(xrange.a);
      safeFree(xrange.b);
      safeFree(yrange.a);
      safeFree(yrange.b);
    } else {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_quad;
      simplified->child1 = g;
    }
  }
  return simplified;
}


node *simplify_halfprecision(node *g) {
  node *simplified;
  mpfr_t *value;
  rangetype xrange, yrange;

  simplified = (node*) safeMalloc(sizeof(node));
  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    simplified->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value,((tools_precision > 64)? tools_precision : 64));
    simplified->value = value;
    mpfr_round_to_halfprecision(*value, *(accessThruMemRef(g)->value));
    if (!mpfr_number_p(*value)) {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_halfprecision;
      simplified->child1 = g;
      mpfr_clear(*value);
      safeFree(value);
    } else {
      free_memory(g);
    }
  } else {
    if (isConstant(g)) {
      xrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      xrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      yrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      yrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*(yrange.a),4* tools_precision);
      mpfr_init2(*(yrange.b),4 * tools_precision);
      mpfr_init2(*(xrange.a),((tools_precision > 64)? tools_precision : 64));
      mpfr_init2(*(xrange.b),((tools_precision > 64)? tools_precision : 64));
      mpfr_set_ui(*(xrange.a),1,GMP_RNDD);
      mpfr_set_ui(*(xrange.b),1,GMP_RNDU);
      evaluateRangeFunction(yrange, g, xrange, 8 * tools_precision);
      mpfr_round_to_halfprecision(*(xrange.a),*(yrange.a));
      mpfr_round_to_halfprecision(*(xrange.b),*(yrange.b));
      if (mpfr_number_p(*(xrange.a)) &&
          mpfr_number_p(*(xrange.b)) &&
          (mpfr_cmp(*(xrange.a),*(xrange.b)) == 0)) {
        simplified->nodeType = CONSTANT;
        value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*value,((tools_precision > 64)? tools_precision : 64));
        simplified->value = value;
        mpfr_set(*value,*(xrange.a),GMP_RNDN); /* Exact */
        free_memory(g);
      } else {
        simplified->nodeType = UNARY_BASE_FUNC;
        simplified->baseFun = basefun_halfprecision;
        simplified->child1 = g;
      }
      mpfr_clear(*(xrange.a));
      mpfr_clear(*(xrange.b));
      mpfr_clear(*(yrange.a));
      mpfr_clear(*(yrange.b));
      safeFree(xrange.a);
      safeFree(xrange.b);
      safeFree(yrange.a);
      safeFree(yrange.b);
    } else {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_halfprecision;
      simplified->child1 = g;
    }
  }
  return simplified;
}


node *simplify_doubledouble(node *g) {
  node *simplified;
  mpfr_t *value;
  rangetype xrange, yrange;

  simplified = (node*) safeMalloc(sizeof(node));
  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    simplified->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value,((tools_precision > 129)? tools_precision : 129));
    simplified->value = value;
    mpfr_round_to_doubledouble(*value, *(accessThruMemRef(g)->value));
    if (!mpfr_number_p(*value)) {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_doubledouble;
      simplified->child1 = g;
      mpfr_clear(*value);
      safeFree(value);
    } else {
      free_memory(g);
    }
  } else {
    if (isConstant(g)) {
      xrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      xrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      yrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      yrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*(yrange.a),4* tools_precision);
      mpfr_init2(*(yrange.b),4 * tools_precision);
      mpfr_init2(*(xrange.a),((tools_precision > 129)? tools_precision : 129));
      mpfr_init2(*(xrange.b),((tools_precision > 129)? tools_precision : 129));
      mpfr_set_ui(*(xrange.a),1,GMP_RNDD);
      mpfr_set_ui(*(xrange.b),1,GMP_RNDU);
      evaluateRangeFunction(yrange, g, xrange, 8 * tools_precision);
      mpfr_round_to_doubledouble(*(xrange.a),*(yrange.a));
      mpfr_round_to_doubledouble(*(xrange.b),*(yrange.b));
      if (mpfr_number_p(*(xrange.a)) &&
          mpfr_number_p(*(xrange.b)) &&
          (mpfr_cmp(*(xrange.a),*(xrange.b)) == 0)) {
        simplified->nodeType = CONSTANT;
        value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*value,((tools_precision > 129)? tools_precision : 129));
        simplified->value = value;
        mpfr_set(*value,*(xrange.a),GMP_RNDN); /* Exact */
        free_memory(g);
      } else {
        simplified->nodeType = UNARY_BASE_FUNC;
        simplified->baseFun = basefun_doubledouble;
        simplified->child1 = g;
      }
      mpfr_clear(*(xrange.a));
      mpfr_clear(*(xrange.b));
      mpfr_clear(*(yrange.a));
      mpfr_clear(*(yrange.b));
      safeFree(xrange.a);
      safeFree(xrange.b);
      safeFree(yrange.a);
      safeFree(yrange.b);
    } else {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_doubledouble;
      simplified->child1 = g;
    }
  }
  return simplified;
}


node *simplify_tripledouble(node *g) {
  node *simplified;
  mpfr_t *value;
  rangetype xrange, yrange;

  simplified = (node*) safeMalloc(sizeof(node));
  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    simplified->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value,((tools_precision > 200)? tools_precision : 200));
    simplified->value = value;
    mpfr_round_to_tripledouble(*value, *(accessThruMemRef(g)->value));
    if (!mpfr_number_p(*value)) {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_tripledouble;
      simplified->child1 = g;
      mpfr_clear(*value);
      safeFree(value);
    } else {
      free_memory(g);
    }
  } else {
    if (isConstant(g)) {
      xrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      xrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      yrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      yrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*(yrange.a),4* tools_precision);
      mpfr_init2(*(yrange.b),4 * tools_precision);
      mpfr_init2(*(xrange.a),((tools_precision > 200)? tools_precision : 200));
      mpfr_init2(*(xrange.b),((tools_precision > 200)? tools_precision : 200));
      mpfr_set_ui(*(xrange.a),1,GMP_RNDD);
      mpfr_set_ui(*(xrange.b),1,GMP_RNDU);
      evaluateRangeFunction(yrange, g, xrange, 8 * tools_precision);
      mpfr_round_to_tripledouble(*(xrange.a),*(yrange.a));
      mpfr_round_to_tripledouble(*(xrange.b),*(yrange.b));
      if (mpfr_number_p(*(xrange.a)) &&
          mpfr_number_p(*(xrange.b)) &&
          (mpfr_cmp(*(xrange.a),*(xrange.b)) == 0)) {
        simplified->nodeType = CONSTANT;
        value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*value,((tools_precision > 200)? tools_precision : 200));
        simplified->value = value;
        mpfr_set(*value,*(xrange.a),GMP_RNDN); /* Exact */
        free_memory(g);
      } else {
        simplified->nodeType = UNARY_BASE_FUNC;
        simplified->baseFun = basefun_tripledouble;
        simplified->child1 = g;
      }
      mpfr_clear(*(xrange.a));
      mpfr_clear(*(xrange.b));
      mpfr_clear(*(yrange.a));
      mpfr_clear(*(yrange.b));
      safeFree(xrange.a);
      safeFree(xrange.b);
      safeFree(yrange.a);
      safeFree(yrange.b);
    } else {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_tripledouble;
      simplified->child1 = g;
    }
  }
  return simplified;
}


node *simplify_doubleextended(node *g) {
  node *simplified;
  mpfr_t *value;
  rangetype xrange, yrange;

  simplified = (node*) safeMalloc(sizeof(node));
  if (accessThruMemRef(g)->nodeType == CONSTANT) {
    simplified->nodeType = CONSTANT;
    value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*value,((tools_precision > 128)? tools_precision : 128));
    simplified->value = value;
    mpfr_round_to_doubleextended(*value, *(accessThruMemRef(g)->value));
    if (!mpfr_number_p(*value)) {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_doubleextended;
      simplified->child1 = g;
      mpfr_clear(*value);
      safeFree(value);
    } else {
      free_memory(g);
    }
  } else {
    if (isConstant(g)) {
      xrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      xrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      yrange.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      yrange.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
      mpfr_init2(*(yrange.a),4* tools_precision);
      mpfr_init2(*(yrange.b),4 * tools_precision);
      mpfr_init2(*(xrange.a),((tools_precision > 128)? tools_precision : 128));
      mpfr_init2(*(xrange.b),((tools_precision > 128)? tools_precision : 128));
      mpfr_set_ui(*(xrange.a),1,GMP_RNDD);
      mpfr_set_ui(*(xrange.b),1,GMP_RNDU);
      evaluateRangeFunction(yrange, g, xrange, 8 * tools_precision);
      mpfr_round_to_doubleextended(*(xrange.a),*(yrange.a));
      mpfr_round_to_doubleextended(*(xrange.b),*(yrange.b));
      if (mpfr_number_p(*(xrange.a)) &&
          mpfr_number_p(*(xrange.b)) &&
          (mpfr_cmp(*(xrange.a),*(xrange.b)) == 0)) {
        simplified->nodeType = CONSTANT;
        value = (mpfr_t*) safeMalloc(sizeof(mpfr_t));
        mpfr_init2(*value,((tools_precision > 128)? tools_precision : 128));
        simplified->value = value;
        mpfr_set(*value,*(xrange.a),GMP_RNDN); /* Exact */
        free_memory(g);
      } else {
        simplified->nodeType = UNARY_BASE_FUNC;
        simplified->baseFun = basefun_doubleextended;
        simplified->child1 = g;
      }
      mpfr_clear(*(xrange.a));
      mpfr_clear(*(xrange.b));
      mpfr_clear(*(yrange.a));
      mpfr_clear(*(yrange.b));
      safeFree(xrange.a);
      safeFree(xrange.b);
      safeFree(yrange.a);
      safeFree(yrange.b);
    } else {
      simplified->nodeType = UNARY_BASE_FUNC;
      simplified->baseFun = basefun_doubleextended;
      simplified->child1 = g;
    }
  }
  return simplified;
}

node *simplify_floor(node *g) { return simplify_integer_rounding_functions(basefun_floor, g); }
node *simplify_ceil(node *g) { return simplify_integer_rounding_functions(basefun_ceil, g); }
node *simplify_nearestint(node *g) { return simplify_integer_rounding_functions(basefun_nearestint, g); }




/******************************************************************************/
/*                                                                            */
/*                     FUNCTIONS FOR EVALUATING THE SIGN                      */
/*                                                                            */
/******************************************************************************/

int evaluateSignTrigoUnsafe(int *s, node *child, int baseFunctionCode) {
  mpfr_t value, value2;
  mpfr_t piHalf;
  mpfr_t dummyX;
  int okay, res;
  node *tempNode;
  int signA;

  okay = 0;

  mpfr_init2(value,defaultprecision);
  mpfr_init2(piHalf,defaultprecision);
  mpfr_init2(dummyX,12);
  mpfr_set_ui(dummyX,1,GMP_RNDN);
  if (evaluateFaithful(value, child, dummyX, defaultprecision) &&
      mpfr_number_p(value)) {
    mpfr_const_pi(piHalf,GMP_RNDN);
    mpfr_div_2ui(piHalf,piHalf,1,GMP_RNDN);
    mpfr_div(value,value,piHalf,GMP_RNDN);
    mpfr_rint(value,value,GMP_RNDN);
    mpfr_div_2ui(value,value,1,GMP_RNDN);
    /* Here, diff is approximately value * pi
       and value * 2 is an integer
    */
    tempNode = makeMul(makeConstant(value),makePi());
    if (compareConstant(&signA, child, tempNode, NULL, 0)) {
      if (signA == 0) {
	/* Here, we have proven that child is equal to value * pi
	 */
	mpfr_init2(value2,defaultprecision);
	mpfr_rint(value2,value,GMP_RNDN);      /* exact, same precision */
	mpfr_sub(value,value,value2,GMP_RNDN); /* exact, Sterbenz */
	/* Here, we know that child is equal to (n + value) * pi for
	   some integer n. We know that value can only be 0 or +/- 0.5
	*/
	switch (baseFunctionCode) {
	case SIN:
	  /* sin is zero for all n * pi, n in Z */
	  if (mpfr_zero_p(value)) {
	    okay = 1;
	    res = 0;
	  }
	  break;
	case COS:
	  /* cos is zero for all (n + 1/2) * pi, n in Z */
	  if (!mpfr_zero_p(value)) {
	    okay = 1;
	    res = 0;
	  }
	  break;
	case TAN:
	  /* tan is zero for all n * pi, n in Z */
	  if (mpfr_zero_p(value)) {
	    okay = 1;
	    res = 0;
	  }
	  break;
	default:
	  sollyaFprintf(stderr,"Error: evaluateSignTrigoUnsafe: unknown identifier (%d) in the tree\n",baseFunctionCode);
	  exit(1);
	}
	mpfr_clear(value2);
      }
    }
    free_memory(tempNode);
  }
  mpfr_clear(dummyX);
  mpfr_clear(piHalf);
  mpfr_clear(value);

  if (okay) *s = res;
  return okay;
}

int sqrt_evalsign(int *sign, node *c) {
  int okayA, signA;
  okayA = evaluateSign(&signA, c);
  if (okayA && (signA >= 0)) {
    *sign = signA;
    return 1;
  }
  return 0;
}

int log_evalsign(int *sign, node *c) {
  int okayA, signA, okayB, signB;
  node * tempNode;
  tempNode = makeConstantDouble(1.0);
  okayA = compareConstant(&signA, c, tempNode, NULL, 0);
  free_memory(tempNode);
  okayB = evaluateSign(&signB, c);
  if (okayA && okayB && (signB > 0)) {
    *sign = signA;
    return 1;
  }
  return 0;
}

int log2_evalsign(int *sign, node *c) { return log_evalsign(sign, c); }
int log10_evalsign(int *sign, node *c) { return log_evalsign(sign, c); }

int sin_evalsign(int *sign, node *c) {
  int okayA, signA;
  okayA = evaluateSign(&signA, c);
  if (okayA && (signA == 0)) {
    *sign = 0;
    return 1;
  }
  /* else */
  return evaluateSignTrigoUnsafe(sign, c, SIN);
}

int cos_evalsign(int *sign, node *c) {
  int okayA, signA;
  okayA = evaluateSign(&signA, c);
  if (okayA && (signA == 0)) {
    *sign = 0;
    return 1;
  }
  /* else */
  return evaluateSignTrigoUnsafe(sign, c, COS);
}

int tan_evalsign(int *sign, node *c) {
  int okayA, signA;
  okayA = evaluateSign(&signA, c);
  if (okayA && (signA == 0)) {
    *sign = 0;
    return 1;
  }
  /* else */
  return evaluateSignTrigoUnsafe(sign, c, TAN);
}


int asin_evalsign(int *sign, node *c) {
  int okayA, signA, okayB, signB;
  node *tempNode, *tempNode2;
  okayA = evaluateSign(&signA, c);
  tempNode = makeAbs(copyTree(c));
  tempNode2 = makeConstantDouble(1.0);
  okayB = compareConstant(&signB, tempNode, tempNode2, NULL, 0);
  free_memory(tempNode);
  free_memory(tempNode2);

  if (okayA && okayB && (signB <= 0)) {
    *sign = signA;
    return 1;
  }
  return 0;
}

int acos_evalsign(int *sign, node *c) {
  int okayA, signA, okayB, signB, okayC, signC;
  node *tempNode, *tempNode2;
  okayA = evaluateSign(&signA, c);
  tempNode = makeAbs(copyTree(c));
  tempNode2 = makeConstantDouble(1.0);
  okayB = compareConstant(&signB, tempNode, tempNode2, NULL, 0);
  okayC = compareConstant(&signC, c, tempNode2, NULL, 0);
  free_memory(tempNode);
  free_memory(tempNode2);

  if (okayA && okayB && okayC && (signB <= 0)) {
    if (signC == 0) *sign = 0; else *sign = 1;
    return 1;
  }
  return 0;
}

int odd_increasing_function_evalsign(int *sign, node *c) {
  return evaluateSign(sign, c);
}

int atan_evalsign(int *sign, node *c) { return odd_increasing_function_evalsign(sign, c); }
int sinh_evalsign(int *sign, node *c) { return odd_increasing_function_evalsign(sign, c); }
int tanh_evalsign(int *sign, node *c) { return odd_increasing_function_evalsign(sign, c); }
int asinh_evalsign(int *sign, node *c) { return odd_increasing_function_evalsign(sign, c); }
int erf_evalsign(int *sign, node *c) { return odd_increasing_function_evalsign(sign, c); }
int expm1_evalsign(int *sign, node *c) { return odd_increasing_function_evalsign(sign, c); }

int positive_function_evalsign(int *sign, node *c) {
 int signA;
  if (evaluateSign(&signA, c)) {
    *sign = 1;
    return 1;
  }
  return 0;
}

int cosh_evalsign(int *sign, node *c) { return positive_function_evalsign(sign, c); }
int erfc_evalsign(int *sign, node *c) { return positive_function_evalsign(sign, c); }
int exp_evalsign(int *sign, node *c) {  return positive_function_evalsign(sign, c); }

int acosh_evalsign(int *sign, node *c) {
  int okayA, signA;
  node *tempNode;
  tempNode = makeConstantDouble(1.0);
  okayA = compareConstant(&signA, c, tempNode, NULL, 0);
  free_memory(tempNode);
  if (okayA && (signA >= 0)) {
    *sign = 1;
    return 1;
  }
  return 0;
}

int atanh_evalsign(int *sign, node *c) {
  int okayA, signA, okayB, signB;
  node *tempNode, *tempNode2;
  okayA = evaluateSign(&signA, c);
  tempNode = makeAbs(copyTree(c));
  tempNode2 = makeConstantDouble(1.0);
  okayB = compareConstant(&signB, tempNode, tempNode2, NULL, 0);
  free_memory(tempNode);
  free_memory(tempNode2);
  if (okayA && okayB && (signB < 0)) {
    *sign = signA;
    return 1;
  }
  return 0;
}

int abs_evalsign(int *sign, node *c) {
  int okayA, signA;
  okayA = evaluateSign(&signA, c);
  if (okayA) {
    if (signA == 0) *sign = 0; else *sign = 1;
    return 1;
  }
  return 0;
}

int void_evalsign(int *sign, node *c) {
  UNUSED_PARAM(sign);
  UNUSED_PARAM(c);
  return 0;
}

int double_evalsign(int *sign, node *c) { return void_evalsign(sign, c); }
int single_evalsign(int *sign, node *c) { return void_evalsign(sign, c); }
int quad_evalsign(int *sign, node *c) { return void_evalsign(sign, c); }
int halfprecision_evalsign(int *sign, node *c) { return void_evalsign(sign, c); }
int doubledouble_evalsign(int *sign, node *c) { return void_evalsign(sign, c); }
int tripledouble_evalsign(int *sign, node *c) { return void_evalsign(sign, c); }
int doubleextended_evalsign(int *sign, node *c) { return void_evalsign(sign, c); }

int log1p_evalsign(int *sign, node *c) {
  int okayA, signA, okayB, signB;
  node *tempNode;
  tempNode = makeConstantDouble(-1.0);
  okayA = compareConstant(&signA, c, tempNode, NULL, 0);
  okayB = evaluateSign(&signB, c);
  free_memory(tempNode);
  if (okayA && okayB && (signA > 0)) {
    *sign = signB;
    return 1;
  }
  return 0;
}

int ceil_evalsign(int *sign, node *c) {
  int okay, s, okayA, signA, okayB, signB;
  node *tempNode;
  okay = 0;
  okayA = evaluateSign(&signA, c);
  tempNode = makeConstantDouble(-1.0);
  if (okayA)
    okayB = compareConstant(&signB, c, tempNode, NULL, 0);
  else
    okayB = 0;
  if (okayA && okayB) {
    okay = 1;
    if (signB <= 0) {
      s = -1;
    } else {
      if (signA <= 0) {
        s = 0;
      } else {
        s = 1;
      }
    }
  }
  free_memory(tempNode);
  if (okay) *sign = s;
  return okay;
}

int floor_evalsign(int *sign, node *c) {
  int okay, s, okayA, signA, okayB, signB;
  node *tempNode;
  okay = 0;
  okayA = evaluateSign(&signA, c);
  tempNode = makeConstantDouble(1.0);
  if (okayA)
    okayB = compareConstant(&signB, c, tempNode, NULL, 0);
  else
    okayB = 0;
  if (okayA && okayB) {
    okay = 1;
    if (signA < 0) {
      s = -1;
    } else {
      if (signB < 0) {
        s = 0;
      } else {
        s = 1;
      }
    }
  }
  free_memory(tempNode);
  if (okay) *sign = s;
  return okay;
}

int nearestint_evalsign(int *sign, node *c) {
  int okay, s, okayA, signA, okayB, signB;
  node *tempNode;
  okay = 0;
  okayA = evaluateSign(&signA, c);
  if (okayA) {
    if (signA == 0) {
      okayB = 1;
      signB = 0;
    }
    else if (signA > 0) {
      tempNode = makeConstantDouble(0.5);
      okayB = compareConstant(&signB, c, tempNode, NULL, 0);
    }
    else {
      tempNode = makeConstantDouble(-0.5);
      okayB = compareConstant(&signB, c, tempNode, NULL, 0);
    }
  }
  else okayB = 0;

  if (okayA && okayB) {
    okay = 1;
    if ( (signA < 0) && (signB < 0) ) s = -1;
    else if ( (signA > 0) && (signB > 0) ) s = 1;
    else s = 0;
  }
  free_memory(tempNode);
  if (okay) *sign = s;
  return okay;
}




/******************************************************************************/
/*                                                                            */
/*                       POINTWISE EVALUATION FUNCTIONS                       */
/*                                                                            */
/******************************************************************************/

int double_point_eval(mpfr_t res, mpfr_t x, mp_rnd_t rnd) {
  mpfr_t myResult;
  int r;
  mpfr_init2(myResult, ((mpfr_get_prec(res) > 64)? mpfr_get_prec(res) : 64));
  mpfr_round_to_double(myResult, x);
  r = mpfr_set(res, myResult, rnd);
  mpfr_clear(myResult);
  return r;
}

int single_point_eval(mpfr_t res, mpfr_t x, mp_rnd_t rnd) {
  mpfr_t myResult;
  int r;
  mpfr_init2(myResult, ((mpfr_get_prec(res) > 64)? mpfr_get_prec(res) : 64));
  mpfr_round_to_single(myResult, x);
  r = mpfr_set(res, myResult, rnd);
  mpfr_clear(myResult);
  return r;
}

int quad_point_eval(mpfr_t res, mpfr_t x, mp_rnd_t rnd) {
  mpfr_t myResult;
  int r;
  mpfr_init2(myResult, ((mpfr_get_prec(res) > 128)? mpfr_get_prec(res) : 128));
  mpfr_round_to_quad(myResult, x);
  r = mpfr_set(res, myResult, rnd);
  mpfr_clear(myResult);
  return r;
}

int halfprecision_point_eval(mpfr_t res, mpfr_t x, mp_rnd_t rnd) {
  mpfr_t myResult;
  int r;
  mpfr_init2(myResult, ((mpfr_get_prec(res) > 64)? mpfr_get_prec(res) : 64));
  mpfr_round_to_halfprecision(myResult, x);
  r = mpfr_set(res, myResult, rnd);
  mpfr_clear(myResult);
  return r;
}

int doubledouble_point_eval(mpfr_t res, mpfr_t x, mp_rnd_t rnd) {
  mpfr_t myResult;
  int r;
  mpfr_init2(myResult, ((mpfr_get_prec(res) > 129)? mpfr_get_prec(res) : 129));
  mpfr_round_to_doubledouble(myResult, x);
  r = mpfr_set(res, myResult, rnd);
  mpfr_clear(myResult);
  return r;
}

int tripledouble_point_eval(mpfr_t res, mpfr_t x, mp_rnd_t rnd) {
  mpfr_t myResult;
  int r;
  mpfr_init2(myResult, ((mpfr_get_prec(res) > 200)? mpfr_get_prec(res) : 200));
  mpfr_round_to_tripledouble(myResult, x);
  r = mpfr_set(res, myResult, rnd);
  mpfr_clear(myResult);
  return r;
}

int doubleextended_point_eval(mpfr_t res, mpfr_t x, mp_rnd_t rnd) {
  mpfr_t myResult;
  int r;
  mpfr_init2(myResult, ((mpfr_get_prec(res) > 128)? mpfr_get_prec(res) : 128));
  mpfr_round_to_doubleextended(myResult, x);
  r = mpfr_set(res, myResult, rnd);
  mpfr_clear(myResult);
  return r;
}






/******************************************************************************/
/*                                                                            */
/*                     DEFINITIONS OF THE BASEFUN OBJECTS                     */
/*                                                                            */
/******************************************************************************/

baseFunction basefun_sqrt_obj = {
  .baseFunctionCode = SQRT,
  .functionName = "sqrt",
  .xmlString = "<root/>\n",
  .mpfrName = "mpfr_sqrt",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 0,
  .onlyZeroIsZero = 1,
  .doesNotVanish = 0,
  .monotonicity = INCREASING,
  .baseAutodiff = sqrt_diff,
  .interval_eval = sollya_mpfi_sqrt,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_sqrt),
  .diff_expr = sqrt_diff_expr,
  .simplify = simplify_sqrt,
  .evalsign = sqrt_evalsign
};
baseFunction *basefun_sqrt = &basefun_sqrt_obj;

baseFunction basefun_exp_obj = {
  .baseFunctionCode = EXP,
  .functionName = "exp",
  .xmlString = "<exp/>\n",
  .mpfrName = "mpfr_exp",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 1,
  .onlyZeroIsZero = 0,
  .doesNotVanish = 1,
  .monotonicity = INCREASING,
  .baseAutodiff = exp_diff,
  .interval_eval = sollya_mpfi_exp,
  .point_eval =  (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_exp),
  .diff_expr = exp_diff_expr,
  .simplify = simplify_exp,
  .evalsign = exp_evalsign
};
baseFunction *basefun_exp = &basefun_exp_obj;

baseFunction basefun_log_obj = {
  .baseFunctionCode = LOG,
  .functionName = "log",
  .xmlString = "<ln/>\n",
  .mpfrName = "mpfr_log",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 0,
  .onlyZeroIsZero = 0,
  .doesNotVanish = 0,
  .monotonicity = INCREASING,
  .baseAutodiff = log_diff,
  .interval_eval = sollya_mpfi_log,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_log),
  .diff_expr = log_diff_expr,
  .simplify = simplify_log,
  .evalsign = log_evalsign
};
baseFunction *basefun_log = &basefun_log_obj;

baseFunction basefun_log2_obj = {
  .baseFunctionCode = LOG_2,
  .functionName = "log2",
  .xmlString = "<log/><logbase><cn>2</cn></logbase>\n",
  .mpfrName = "mpfr_log2",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 0,
  .onlyZeroIsZero = 0,
  .doesNotVanish = 0,
  .monotonicity = INCREASING,
  .baseAutodiff = log2_diff,
  .interval_eval = sollya_mpfi_log2,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_log2),
  .diff_expr = log2_diff_expr,
  .simplify = simplify_log2,
  .evalsign = log2_evalsign
};
baseFunction *basefun_log2 = &basefun_log2_obj;

baseFunction basefun_log10_obj = {
  .baseFunctionCode = LOG_10,
  .functionName = "log10",
  .xmlString = "<log/>\n",
  .mpfrName = "mpfr_log10",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 0,
  .onlyZeroIsZero = 0,
  .doesNotVanish = 0,
  .monotonicity = INCREASING,
  .baseAutodiff = log10_diff,
  .interval_eval = sollya_mpfi_log10,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_log10),
  .diff_expr = log10_diff_expr,
  .simplify = simplify_log10,
  .evalsign = log10_evalsign
};
baseFunction *basefun_log10 = &basefun_log10_obj;

baseFunction basefun_sin_obj = {
  .baseFunctionCode = SIN,
  .functionName = "sin",
  .xmlString = "<sin/>\n",
  .mpfrName = "mpfr_sin",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 1,
  .onlyZeroIsZero = 0,
  .doesNotVanish = 0,
  .monotonicity = MONOTONICITY_NONE,
  .baseAutodiff = sin_diff,
  .interval_eval = sollya_mpfi_sin,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_sin),
  .diff_expr = sin_diff_expr,
  .simplify = simplify_sin,
  .evalsign = sin_evalsign
};
baseFunction *basefun_sin = &basefun_sin_obj;

baseFunction basefun_cos_obj = {
  .baseFunctionCode = COS,
  .functionName = "cos",
  .xmlString = "<cos/>\n",
  .mpfrName = "mpfr_cos",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 1,
  .onlyZeroIsZero = 0,
  .doesNotVanish = 0,
  .monotonicity = MONOTONICITY_NONE,
  .baseAutodiff = cos_diff,
  .interval_eval = sollya_mpfi_cos,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_cos),
  .diff_expr = cos_diff_expr,
  .simplify = simplify_cos,
  .evalsign = cos_evalsign
};
baseFunction *basefun_cos = &basefun_cos_obj;

baseFunction basefun_tan_obj = {
  .baseFunctionCode = TAN,
  .functionName = "tan",
  .xmlString = "<tan/>\n",
  .mpfrName = "mpfr_tan",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 0,
  .onlyZeroIsZero = 0,
  .doesNotVanish = 0,
  .monotonicity = MONOTONICITY_NONE,
  .baseAutodiff = tan_diff,
  .interval_eval = sollya_mpfi_tan,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_tan),
  .diff_expr = tan_diff_expr,
  .simplify = simplify_tan,
  .evalsign = tan_evalsign
};
baseFunction *basefun_tan = &basefun_tan_obj;

baseFunction basefun_asin_obj = {
  .baseFunctionCode = ASIN,
  .functionName = "asin",
  .xmlString = "<arcsin/>\n",
  .mpfrName = "mpfr_asin",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 0,
  .onlyZeroIsZero = 1,
  .doesNotVanish = 0,
  .monotonicity = INCREASING,
  .baseAutodiff = asin_diff,
  .interval_eval = sollya_mpfi_asin,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_asin),
  .diff_expr = asin_diff_expr,
  .simplify = simplify_asin,
  .evalsign = asin_evalsign
};
baseFunction *basefun_asin = &basefun_asin_obj;

baseFunction basefun_acos_obj = {
  .baseFunctionCode = ACOS,
  .functionName = "acos",
  .xmlString = "<arccos/>\n",
  .mpfrName = "mpfr_acos",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 0,
  .onlyZeroIsZero = 0,
  .doesNotVanish = 0,
  .monotonicity = DECREASING,
  .baseAutodiff = acos_diff,
  .interval_eval = sollya_mpfi_acos,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_acos),
  .diff_expr = acos_diff_expr,
  .simplify = simplify_acos,
  .evalsign = acos_evalsign
};
baseFunction *basefun_acos = &basefun_acos_obj;

baseFunction basefun_atan_obj = {
  .baseFunctionCode = ATAN,
  .functionName = "atan",
  .xmlString = "<arctan/>\n",
  .mpfrName = "mpfr_atan",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 1,
  .onlyZeroIsZero = 1,
  .doesNotVanish = 0,
  .monotonicity = INCREASING,
  .baseAutodiff = atan_diff,
  .interval_eval = sollya_mpfi_atan,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_atan),
  .diff_expr = atan_diff_expr,
  .simplify = simplify_atan,
  .evalsign = atan_evalsign
};
baseFunction *basefun_atan = &basefun_atan_obj;

baseFunction basefun_sinh_obj = {
  .baseFunctionCode = SINH,
  .functionName = "sinh",
  .xmlString = "<sinh/>\n",
  .mpfrName = "mpfr_sinh",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 1,
  .onlyZeroIsZero = 1,
  .doesNotVanish = 0,
  .monotonicity = INCREASING,
  .baseAutodiff = sinh_diff,
  .interval_eval = sollya_mpfi_sinh,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_sinh),
  .diff_expr = sinh_diff_expr,
  .simplify = simplify_sinh,
  .evalsign = sinh_evalsign
};
baseFunction *basefun_sinh = &basefun_sinh_obj;

baseFunction basefun_cosh_obj = {
  .baseFunctionCode = COSH,
  .functionName = "cosh",
  .xmlString = "<cosh/>\n",
  .mpfrName = "mpfr_cosh",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 1,
  .onlyZeroIsZero = 0,
  .doesNotVanish = 1,
  .monotonicity = MONOTONICITY_NONE,
  .baseAutodiff = cosh_diff,
  .interval_eval = sollya_mpfi_cosh,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_cosh),
  .diff_expr = cosh_diff_expr,
  .simplify = simplify_cosh,
  .evalsign = cosh_evalsign
};
baseFunction *basefun_cosh = &basefun_cosh_obj;

baseFunction basefun_tanh_obj = {
  .baseFunctionCode = TANH,
  .functionName = "tanh",
  .xmlString = "<tanh/>\n",
  .mpfrName = "mpfr_tanh",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 1,
  .onlyZeroIsZero = 1,
  .doesNotVanish = 0,
  .monotonicity = INCREASING,
  .baseAutodiff = tanh_diff,
  .interval_eval = sollya_mpfi_tanh,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_tanh),
  .diff_expr = tanh_diff_expr,
  .simplify = simplify_tanh,
  .evalsign = tanh_evalsign
};
baseFunction *basefun_tanh = &basefun_tanh_obj;

baseFunction basefun_asinh_obj = {
  .baseFunctionCode = ASINH,
  .functionName = "asinh",
  .xmlString = "<arcsinh/>\n",
  .mpfrName = "mpfr_asinh",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 1,
  .onlyZeroIsZero = 1,
  .doesNotVanish = 0,
  .monotonicity = INCREASING,
  .baseAutodiff = asinh_diff,
  .interval_eval = sollya_mpfi_asinh,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_asinh),
  .diff_expr = asinh_diff_expr,
  .simplify = simplify_asinh,
  .evalsign = asinh_evalsign
};
baseFunction *basefun_asinh = &basefun_asinh_obj;

baseFunction basefun_acosh_obj = {
  .baseFunctionCode = ACOSH,
  .functionName = "acosh",
  .xmlString = "<arccosh/>\n",
  .mpfrName = "mpfr_acosh",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 0,
  .onlyZeroIsZero = 0,
  .doesNotVanish = 0,
  .monotonicity = INCREASING,
  .baseAutodiff = acosh_diff,
  .interval_eval = sollya_mpfi_acosh,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_acosh),
  .diff_expr = acosh_diff_expr,
  .simplify = simplify_acosh,
  .evalsign = acosh_evalsign
};
baseFunction *basefun_acosh = &basefun_acosh_obj;

baseFunction basefun_atanh_obj = {
  .baseFunctionCode = ATANH,
  .functionName = "atanh",
  .xmlString = "<arctanh/>\n",
  .mpfrName = "mpfr_atanh",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 0,
  .onlyZeroIsZero = 1,
  .doesNotVanish = 0,
  .monotonicity = INCREASING,
  .baseAutodiff = atanh_diff,
  .interval_eval = sollya_mpfi_atanh,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_atanh),
  .diff_expr = atanh_diff_expr,
  .simplify = simplify_atanh,
  .evalsign = atanh_evalsign
};
baseFunction *basefun_atanh = &basefun_atanh_obj;

baseFunction basefun_abs_obj = {
  .baseFunctionCode = ABS,
  .functionName = "abs",
  .xmlString = "<abs/>\n",
  .mpfrName = "mpfr_abs",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 1,
  .onlyZeroIsZero = 1,
  .doesNotVanish = 0,
  .monotonicity = MONOTONICITY_NONE,
  .baseAutodiff = abs_diff,
  .interval_eval = sollya_mpfi_abs,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_abs),
  .diff_expr = abs_diff_expr,
  .simplify = simplify_abs,
  .evalsign = abs_evalsign
};
baseFunction *basefun_abs = &basefun_abs_obj;

baseFunction basefun_erf_obj = {
  .baseFunctionCode = ERF,
  .functionName = "erf",
  .xmlString = "<csymbol definitionURL=\"http://www.openmath.org/CDs/errorFresnelInts.ocd\" encoding=\"OpenMath\">erf</csymbol>\n",
  .mpfrName = "mpfr_erf",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 1,
  .onlyZeroIsZero = 1,
  .doesNotVanish = 0,
  .monotonicity = INCREASING,
  .baseAutodiff = erf_diff,
  .interval_eval = sollya_mpfi_erf,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_erf),
  .diff_expr = erf_diff_expr,
  .simplify = simplify_erf,
  .evalsign = erf_evalsign
};
baseFunction *basefun_erf = &basefun_erf_obj;

baseFunction basefun_erfc_obj = {
  .baseFunctionCode = ERFC,
  .functionName = "erfc",
  .xmlString = "<csymbol definitionURL=\"http://www.openmath.org/CDs/errorFresnelInts.ocd\" encoding=\"OpenMath\">erfc</csymbol>\n",
  .mpfrName = "mpfr_erfc",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 1,
  .onlyZeroIsZero = 0,
  .doesNotVanish = 1,
  .monotonicity = DECREASING,
  .baseAutodiff = erfc_diff,
  .interval_eval = sollya_mpfi_erfc,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_erfc),
  .diff_expr = erfc_diff_expr,
  .simplify = simplify_erfc,
  .evalsign = erfc_evalsign
};
baseFunction *basefun_erfc = &basefun_erfc_obj;

baseFunction basefun_log1p_obj = {
  .baseFunctionCode = LOG_1P,
  .functionName = "log1p",
  .xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">log1p</csymbol>\n",
  .mpfrName = "mpfr_log1p",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 0,
  .onlyZeroIsZero = 0,
  .doesNotVanish = 0,
  .monotonicity = INCREASING,
  .baseAutodiff = log1p_diff,
  .interval_eval = sollya_mpfi_log1p,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_log1p),
  .diff_expr = log1p_diff_expr,
  .simplify = simplify_log1p,
  .evalsign = log1p_evalsign
};
baseFunction *basefun_log1p = &basefun_log1p_obj;

baseFunction basefun_expm1_obj = {
  .baseFunctionCode = EXP_M1,
  .functionName = "expm1",
  .xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">expm1</csymbol>\n",
  .mpfrName = "mpfr_expm1",
  .handledByImplementconst = 1,
  .isDefinedEverywhere = 1,
  .onlyZeroIsZero = 1,
  .doesNotVanish = 0,
  .monotonicity = INCREASING,
  .baseAutodiff = expm1_diff,
  .interval_eval = sollya_mpfi_expm1,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_expm1),
  .diff_expr = expm1_diff_expr,
  .simplify = simplify_expm1,
  .evalsign = expm1_evalsign
};
baseFunction *basefun_expm1 = &basefun_expm1_obj;

baseFunction basefun_double_obj = {
  .baseFunctionCode = DOUBLE,
  .functionName = "double",
  .xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">double</csymbol>\n",
  .mpfrName = "",
  .handledByImplementconst = 0,
  .isDefinedEverywhere = 0, /* because of overflows: for large numbers, it becomes +Inf */
  .onlyZeroIsZero = 0,
  .doesNotVanish = 0,
  .monotonicity = NONDECREASING,
  .baseAutodiff = double_diff,
  .interval_eval = sollya_mpfi_round_to_double,
  .point_eval = double_point_eval,
  .diff_expr = double_diff_expr,
  .simplify = simplify_double,
  .evalsign = double_evalsign
};
baseFunction *basefun_double = &basefun_double_obj;

baseFunction basefun_single_obj = {
  .baseFunctionCode = SINGLE,
  .functionName = "single",
  .xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">single</csymbol>\n",
  .mpfrName = "",
  .handledByImplementconst = 0,
  .isDefinedEverywhere = 0, /* because of overflows: for large numbers, it becomes +Inf */
  .onlyZeroIsZero = 0,
  .doesNotVanish = 0,
  .monotonicity = NONDECREASING,
  .baseAutodiff = single_diff,
  .interval_eval = sollya_mpfi_round_to_single,
  .point_eval = single_point_eval,
  .diff_expr = single_diff_expr,
  .simplify = simplify_single,
  .evalsign = single_evalsign
};
baseFunction *basefun_single = &basefun_single_obj;

baseFunction basefun_halfprecision_obj = {
  .baseFunctionCode = HALFPRECISION,
  .functionName = "halfprecision",
  .xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">halfprecision</csymbol>\n",
  .mpfrName = "",
  .handledByImplementconst = 0,
  .isDefinedEverywhere = 0, /* because of overflows: for large numbers, it becomes +Inf */
  .onlyZeroIsZero = 0,
  .doesNotVanish = 0,
  .monotonicity = NONDECREASING,
  .baseAutodiff = halfprecision_diff,
  .interval_eval = sollya_mpfi_round_to_halfprecision,
  .point_eval = halfprecision_point_eval,
  .diff_expr = halfprecision_diff_expr,
  .simplify = simplify_halfprecision,
  .evalsign = halfprecision_evalsign
};
baseFunction *basefun_halfprecision = &basefun_halfprecision_obj;

baseFunction basefun_quad_obj = {
  .baseFunctionCode = QUAD,
  .functionName = "quad",
  .xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">quad</csymbol>\n",
  .mpfrName = "",
  .handledByImplementconst = 0,
  .isDefinedEverywhere = 0, /* because of overflows: for large numbers, it becomes +Inf */
  .onlyZeroIsZero = 0,
  .doesNotVanish = 0,
  .monotonicity = NONDECREASING,
  .baseAutodiff = quad_diff,
  .interval_eval = sollya_mpfi_round_to_quad,
  .point_eval = quad_point_eval,
  .diff_expr = quad_diff_expr,
  .simplify = simplify_quad,
  .evalsign = quad_evalsign
};
baseFunction *basefun_quad = &basefun_quad_obj;

baseFunction basefun_doubledouble_obj = {
  .baseFunctionCode = DOUBLEDOUBLE,
  .functionName = "doubledouble",
  .xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">doubledouble</csymbol>\n",
  .mpfrName = "",
  .handledByImplementconst = 0,
  .isDefinedEverywhere = 0, /* because of overflows: for large numbers, it becomes +Inf */
  .onlyZeroIsZero = 0,
  .doesNotVanish = 0,
  .monotonicity = NONDECREASING,
  .baseAutodiff = doubledouble_diff,
  .interval_eval = sollya_mpfi_round_to_doubledouble,
  .point_eval = doubledouble_point_eval,
  .diff_expr = doubledouble_diff_expr,
  .simplify = simplify_doubledouble,
  .evalsign = doubledouble_evalsign
};
baseFunction *basefun_doubledouble = &basefun_doubledouble_obj;

baseFunction basefun_tripledouble_obj = {
  .baseFunctionCode = TRIPLEDOUBLE,
  .functionName = "tripledouble",
  .xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">tripledouble</csymbol>\n",
  .mpfrName = "",
  .handledByImplementconst = 0,
  .isDefinedEverywhere = 0, /* because of overflows: for large numbers, it becomes +Inf */
  .onlyZeroIsZero = 0,
  .doesNotVanish = 0,
  .monotonicity = NONDECREASING,
  .baseAutodiff = tripledouble_diff,
  .interval_eval = sollya_mpfi_round_to_tripledouble,
  .point_eval = tripledouble_point_eval,
  .diff_expr = tripledouble_diff_expr,
  .simplify = simplify_tripledouble,
  .evalsign = tripledouble_evalsign
};
baseFunction *basefun_tripledouble = &basefun_tripledouble_obj;

baseFunction basefun_doubleextended_obj = {
  .baseFunctionCode = DOUBLEEXTENDED,
  .functionName = "doubleextended",
  .xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">doubleextended</csymbol>\n",
  .mpfrName = "",
  .handledByImplementconst = 0,
  .isDefinedEverywhere = 0, /* because of overflows: for large numbers, it becomes +Inf */
  .onlyZeroIsZero = 0,
  .doesNotVanish = 0,
  .monotonicity = NONDECREASING,
  .baseAutodiff = doubleextended_diff,
  .interval_eval = sollya_mpfi_round_to_doubleextended,
  .point_eval = doubleextended_point_eval,
  .diff_expr = doubleextended_diff_expr,
  .simplify = simplify_doubleextended,
  .evalsign = doubleextended_evalsign
};
baseFunction *basefun_doubleextended = &basefun_doubleextended_obj;

baseFunction basefun_ceil_obj = {
  .baseFunctionCode = CEIL,
  .functionName = "ceil",
  .xmlString = "<ceiling/>\n",
  .mpfrName = "",
  .handledByImplementconst = 0,
  .isDefinedEverywhere = 1,
  .onlyZeroIsZero = 0,
  .doesNotVanish = 0,
  .monotonicity = NONDECREASING,
  .baseAutodiff = ceil_diff,
  .interval_eval = sollya_mpfi_ceil,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_ceil),
  .diff_expr = ceil_diff_expr,
  .simplify = simplify_ceil,
  .evalsign = ceil_evalsign
};
baseFunction *basefun_ceil = &basefun_ceil_obj;

baseFunction basefun_floor_obj = {
  .baseFunctionCode = FLOOR,
  .functionName = "floor",
  .xmlString = "<floor/>\n",
  .mpfrName = "",
  .handledByImplementconst = 0,
  .isDefinedEverywhere = 1,
  .doesNotVanish = 0,
  .monotonicity = NONDECREASING,
  .baseAutodiff = floor_diff,
  .onlyZeroIsZero = 0,
  .interval_eval = sollya_mpfi_floor,
  .point_eval = (int (*)(mpfr_t, mpfr_t, mp_rnd_t))(mpfr_floor),
  .diff_expr = floor_diff_expr,
  .simplify = simplify_floor,
  .evalsign = floor_evalsign
};
baseFunction *basefun_floor = &basefun_floor_obj;

baseFunction basefun_nearestint_obj = {
  .baseFunctionCode = NEARESTINT,
  .functionName = "nearestint",
  .xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">nearestint</csymbol>\n",
  .mpfrName = "",
  .handledByImplementconst = 0,
  .isDefinedEverywhere = 1,
  .onlyZeroIsZero = 0,
  .doesNotVanish = 0,
  .monotonicity = NONDECREASING,
  .baseAutodiff = nearestint_diff,
  .interval_eval = sollya_mpfi_nearestint,
  .point_eval = sollya_mpfr_rint_nearestint,
  .diff_expr = nearestint_diff_expr,
  .simplify = simplify_nearestint,
  .evalsign = nearestint_evalsign
};
baseFunction *basefun_nearestint = &basefun_nearestint_obj;






/******************************************************************************/
/*                                                                            */
/*                     SHORTCUTS FOR CONSTRUCTING FUNCTIONS                   */
/*                                                                            */
/******************************************************************************/

node *makeSqrt(node *op1) {
  return makeUnary(op1,basefun_sqrt);
}

node *makeExp(node *op1) {
  return makeUnary(op1,basefun_exp);
}

node *makeLog(node *op1) {
  return makeUnary(op1,basefun_log);
}

node *makeLog2(node *op1) {
  return makeUnary(op1,basefun_log2);
}

node *makeLog10(node *op1) {
  return makeUnary(op1,basefun_log10);
}

node *makeSin(node *op1) {
  return makeUnary(op1,basefun_sin);
}

node *makeCos(node *op1) {
  return makeUnary(op1,basefun_cos);
}

node *makeTan(node *op1) {
  return makeUnary(op1,basefun_tan);
}

node *makeAsin(node *op1) {
  return makeUnary(op1,basefun_asin);
}

node *makeAcos(node *op1) {
  return makeUnary(op1,basefun_acos);
}

node *makeAtan(node *op1) {
  return makeUnary(op1,basefun_atan);
}

node *makeAbs(node *op1) {
  return makeUnary(op1,basefun_abs);
}

node *makeDouble(node *op1) {
  return makeUnary(op1,basefun_double);
}

node *makeSingle(node *op1) {
  return makeUnary(op1,basefun_single);
}

node *makeQuad(node *op1) {
  return makeUnary(op1,basefun_quad);
}

node *makeHalfPrecision(node *op1) {
  return makeUnary(op1,basefun_halfprecision);
}

node *makeDoubledouble(node *op1) {
  return makeUnary(op1,basefun_doubledouble);
}

node *makeTripledouble(node *op1) {
  return makeUnary(op1,basefun_tripledouble);
}

node *makeErf(node *op1 ) {
  return makeUnary(op1,basefun_erf);
}

node *makeErfc(node *op1) {
  return makeUnary(op1,basefun_erfc);
}

node *makeLog1p(node *op1) {
  return makeUnary(op1,basefun_log1p);
}

node *makeExpm1(node *op1) {
  return makeUnary(op1,basefun_expm1);
}

node *makeDoubleextended(node *op1) {
  return makeUnary(op1,basefun_doubleextended);
}

node *makeCeil(node *op1) {
  return makeUnary(op1,basefun_ceil);
}

node *makeFloor(node *op1) {
  return makeUnary(op1,basefun_floor);
}

node *makeNearestInt(node *op1) {
  return makeUnary(op1,basefun_nearestint);
}

node *makeSinh(node *op1) {
  return makeUnary(op1,basefun_sinh);
}

node *makeCosh(node *op1) {
  return makeUnary(op1,basefun_cosh);
}

node *makeTanh(node *op1) {
  return makeUnary(op1,basefun_tanh);
}

node *makeAsinh(node *op1) {
  return makeUnary(op1,basefun_asinh);
}

node *makeAcosh(node *op1) {
  return makeUnary(op1,basefun_acosh);
}

node *makeAtanh(node *op1) {
  return makeUnary(op1,basefun_atanh);
}
