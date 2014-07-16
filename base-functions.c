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

#include "general.h"
#include "expression.h"
#include "base-functions.h"
#include "autodiff.h"
#include "mpfi-compat.h"


void sqrt_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent) {
  mpfr_t oneHalf;
  mpfr_init2(oneHalf, prec);
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

void double_double_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
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

void triple_double_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
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

void double_extended_diff(sollya_mpfi_t *res, sollya_mpfi_t x, int n, int *silent){
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

/*

functionName = "sqrt";
functionName = "exp";
functionName = "log";
functionName = "log2";
functionName = "log10";
functionName = "sin";
functionName = "cos";
functionName = "tan";
functionName = "asin";
functionName = "acos";
functionName = "atan";
functionName = "sinh";
functionName = "cosh";
functionName = "tanh";
functionName = "asinh";
functionName = "acosh";
functionName = "atanh";
functionName = "-";
functionName = "abs";
functionName = "double";
functionName = "single";
functionName = "halfprecision";
functionName = "quad";
functionName = "doubledouble";
functionName = "tripledouble";
functionName = "erf";
functionName = "erfc";
functionName = "log1p";
functionName = "expm1";
functionName = "doubleextended";
functionName = "ceil";
functionName = "floor";
functionName = "nearestint";



xmlString = "<root/>\n";
xmlString = "<exp/>\n";
xmlString = "<ln/>\n";
xmlString = "<log/><logbase><cn>2</cn></logbase>\n";
xmlString = "<log/>\n";
xmlString = "<sin/>\n";
xmlString = "<cos/>\n";
xmlString = "<tan/>\n";
xmlString = "<arcsin/>\n";
xmlString = "<arccos/>\n";
xmlString = "<arctan/>\n";
xmlString = "<sinh/>\n";
xmlString = "<cosh/>\n";
xmlString = "<tanh/>\n";
xmlString = "<arcsinh/>\n";
xmlString = "<arccosh/>\n";
xmlString = "<arctanh/>\n";
xmlString = "<abs/>\n";
xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">double</csymbol>\n";
xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">single</csymbol>\n";
xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">quad</csymbol>\n";
xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">halfprecision</csymbol>\n";
xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">doubledouble</csymbol>\n";
xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">tripledouble</csymbol>\n";
xmlString = "<csymbol definitionURL=\"http://www.openmath.org/CDs/errorFresnelInts.ocd\" encoding=\"OpenMath\">erf</csymbol>\n";
xmlString = "<csymbol definitionURL=\"http://www.openmath.org/CDs/errorFresnelInts.ocd\" encoding=\"OpenMath\">erfc</csymbol>\n";
xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">log1p</csymbol>\n";
xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">expm1</csymbol>\n";
xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">doubleextended</csymbol>\n";
xmlString = "<ceiling/>\n";
xmlString = "<floor/>\n";
xmlString = "<csymbol definitionURL=\"http://www.google.com/\" encoding=\"OpenMath\">nearestint</csymbol>\n";

interval_eval = sollya_mpfi_exp;
interval_eval = sollya_mpfi_log;
interval_eval = sollya_mpfi_log2;
interval_eval = sollya_mpfi_log10;
interval_eval = sollya_mpfi_sin;
interval_eval = sollya_mpfi_cos;
interval_eval = sollya_mpfi_tan;
interval_eval = sollya_mpfi_asin;
interval_eval = sollya_mpfi_acos;
interval_eval = sollya_mpfi_atan;
interval_eval = sollya_mpfi_sinh;
interval_eval = sollya_mpfi_cosh;
interval_eval = sollya_mpfi_tanh;
interval_eval = sollya_mpfi_asinh;
interval_eval = sollya_mpfi_acosh;
interval_eval = sollya_mpfi_atanh;
interval_eval = sollya_mpfi_abs;
interval_eval = sollya_mpfi_round_to_double;
interval_eval = sollya_mpfi_round_to_single;
interval_eval = sollya_mpfi_round_to_halfprecision;
interval_eval = sollya_mpfi_round_to_quad;
interval_eval = sollya_mpfi_round_to_doubledouble;
interval_eval = sollya_mpfi_round_to_tripledouble;
interval_eval = sollya_mpfi_erf;
interval_eval = sollya_mpfi_erfc;
interval_eval = sollya_mpfi_log1p;
interval_eval = sollya_mpfi_expm1;
interval_eval = sollya_mpfi_round_to_doubleextended;
interval_eval = sollya_mpfi_ceil;
interval_eval = sollya_mpfi_floor;
interval_eval = sollya_mpfi_nearestint;
*/

baseFunction *basefun_sqrt;
baseFunction *basefun_exp;
baseFunction *basefun_log;
baseFunction *basefun_log2;
baseFunction *basefun_log10;
baseFunction *basefun_sin;
baseFunction *basefun_cos;
baseFunction *basefun_tan;
baseFunction *basefun_asin;
baseFunction *basefun_acos;
baseFunction *basefun_atan;
baseFunction *basefun_sinh;
baseFunction *basefun_cosh;
baseFunction *basefun_tanh;
baseFunction *basefun_asinh;
baseFunction *basefun_acosh;
baseFunction *basefun_atanh;
baseFunction *basefun_abs;
baseFunction *basefun_double;
baseFunction *basefun_single;
baseFunction *basefun_halfprecision;
baseFunction *basefun_quad;
baseFunction *basefun_doubledouble;
baseFunction *basefun_tripledouble;
baseFunction *basefun_erf;
baseFunction *basefun_erfc;
baseFunction *basefun_log1p;
baseFunction *basefun_expm1;
baseFunction *basefun_doubleextended;
baseFunction *basefun_ceil;
baseFunction *basefun_floor;
baseFunction *basefun_nearestint;

/* Shortcuts */
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
