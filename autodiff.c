/*

  Copyright 2008-2012 by

  Laboratoire de l'Informatique du Parallelisme,
  UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668,

  LORIA (CNRS, INPL, INRIA, UHP, U-Nancy 2),

  Laboratoire d'Informatique de Paris 6, equipe PEQUAN,
  UPMC Universite Paris 06 - CNRS - UMR 7606 - LIP6, Paris, France

  and by

  Centre de recherche INRIA Sophia-Antipolis Mediterranee, equipe APICS,
  Sophia Antipolis, France.

  Contributors S. Chevillard, M. Joldes, Ch. Lauter

  sylvain.chevillard@ens-lyon.org
  joldes@laas.fr
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

#include "autodiff.h"
#include "general.h"
#include "infnorm.h"
#include "execute.h"
#include <stdlib.h>

/* The functions in AD manipulate arrays of size (n+1): [u0...un] */
/* Each array is supposed to represent a function u at a given    */
/* point x0. The value ui is a small sollya_mpfi_t such that             */
/*            u^(i)(x0)/i!  belongs to ui                         */


/* Apply Leibniz' formula: (uv)_p = sum(i=0..p, u_i * v_(p-i))    */
void multiplication_AD(sollya_mpfi_t *res, sollya_mpfi_t *f, sollya_mpfi_t *g, int n) {
  int i,j,p;
  sollya_mpfi_t temp;
  mp_prec_t prec;
  sollya_mpfi_t *temp_array;

  prec = getToolPrecision();
  sollya_mpfi_init2(temp, prec);

  temp_array = (sollya_mpfi_t *)safeCalloc((n+1),sizeof(sollya_mpfi_t));
  for(p=0;p<=n;p++) sollya_mpfi_init2(temp_array[p], prec);

  for(p=0; p<=n; p++) {
    sollya_mpfi_set_ui(temp_array[p], 0);
    i=0; j=p;
    while(i<=p) {
      sollya_mpfi_mul(temp, f[i], g[j]);
      sollya_mpfi_add(temp_array[p], temp_array[p], temp);
      i++; j--;
    }
  }

  for(p=0; p<=n; p++) {
    sollya_mpfi_set(res[p], temp_array[p]);
    sollya_mpfi_clear(temp_array[p]);
  }
  safeFree(temp_array);
  sollya_mpfi_clear(temp);
  return;
}

/* Generic recursive algorithm for the successive derivation of (g o f)  */
/* The array [f0...fn] represents a function f at point x0               */
/* The array [g0...gn] represents a function g at point f(x0)            */
/* Algo:                                                                 */
/*    If n==0, return [g0]                                               */
/*    Else, (g o f)^(i+1) = ((g o f)')^(i) = (f' * (g' o f))^(i) = h^(i) */
/*      So (g o f)^(i+1) / (i+1)! = (1/(i+1)) * h^(i)/i!                 */
/*      So we compute the array [w0 ... w(n-1)] corresponding to         */
/*      (g' o f), up to order (n-1) (by a recursive call)                */
/*      We apply multiplication_AD to w and [(1*f1) ... (n*fn)]          */
/*      (we remark that [(1*f1) ... (n*fn)] corresponds to f')           */
/*      This leads to an array [h0...h(n-1)] corresponding to h          */
/*      Finally, we return [g0 (h0/1) ... (h(n-1)/n)]                    */
void composition_AD(sollya_mpfi_t *res, sollya_mpfi_t *g, sollya_mpfi_t *f, int n) {
  sollya_mpfi_t *fprime, *gprime;
  sollya_mpfi_t *temp_array;
  int i;
  mp_prec_t prec;

  prec = getToolPrecision();
  if(n==0) sollya_mpfi_set(res[0], g[0]);
  else {
    temp_array = (sollya_mpfi_t *)safeCalloc(n,sizeof(sollya_mpfi_t));
    fprime = (sollya_mpfi_t *)safeCalloc(n,sizeof(sollya_mpfi_t));
    gprime = (sollya_mpfi_t *)safeCalloc(n,sizeof(sollya_mpfi_t));
    for(i=0;i<=n-1;i++) {
      sollya_mpfi_init2(temp_array[i], prec);
      sollya_mpfi_init2(fprime[i], prec);
      sollya_mpfi_init2(gprime[i], prec);

      sollya_mpfi_mul_ui(fprime[i], f[i+1], i+1);
      sollya_mpfi_mul_ui(gprime[i], g[i+1], i+1);
    }

    composition_AD(temp_array, gprime, f, n-1);
    multiplication_AD(res+1, temp_array, fprime, n-1);

    sollya_mpfi_set(res[0], g[0]);
    for(i=1; i<=n; i++) sollya_mpfi_div_ui(res[i], res[i], i);
    for(i=0; i<=n-1; i++) {
      sollya_mpfi_clear(temp_array[i]);
      sollya_mpfi_clear(fprime[i]);
      sollya_mpfi_clear(gprime[i]);
    }

    safeFree(temp_array);
    safeFree(fprime);
    safeFree(gprime);
  }

  return ;
}


void binary_function_diff(sollya_mpfi_t *res, int nodeType, sollya_mpfi_t x0, node *f, node *g, int n, int *silent) {
  int i;
  sollya_mpfi_t *res1, *res2, *temp_array;
  mpfr_t minusOne;
  mp_prec_t prec;

  prec = getToolPrecision();
  res1 = (sollya_mpfi_t *)safeCalloc((n+1),sizeof(sollya_mpfi_t));
  res2 = (sollya_mpfi_t *)safeCalloc((n+1),sizeof(sollya_mpfi_t));
  for(i=0;i<=n;i++) {
    sollya_mpfi_init2(res1[i], prec);
    sollya_mpfi_init2(res2[i], prec);
  }
  auto_diff_scaled(res1, f, x0, n);
  auto_diff_scaled(res2, g, x0, n);

  switch(nodeType) {
  case ADD:
    for(i=0; i<=n; i++) sollya_mpfi_add(res[i], res1[i], res2[i]);
    break;
  case SUB:
    for(i=0; i<=n; i++) sollya_mpfi_sub(res[i], res1[i], res2[i]);
    break;
  case MUL:
    multiplication_AD(res, res1, res2, n);
    break;
  case DIV: /* We compute it by g/h = g * h^{-1} */
    temp_array = (sollya_mpfi_t *)safeCalloc((n+1),sizeof(sollya_mpfi_t));
    for(i=0;i<=n;i++) sollya_mpfi_init2(temp_array[i], prec);

    /* temp_array corresponds to x->1/x at point h(x0) */
    mpfr_init2(minusOne, prec);  mpfr_set_si(minusOne, -1, GMP_RNDN);
    constantPower_diff(temp_array, res2[0], minusOne, n, silent);
    mpfr_clear(minusOne);

    /* temp_array corresponds to (x->1/x)(h) = 1/h */
    composition_AD(temp_array, temp_array, res2, n);

    /* res corresponds to g * 1/h */
    multiplication_AD(res, res1, temp_array, n);

    for(i=0;i<=n;i++) sollya_mpfi_clear(temp_array[i]);
    safeFree(temp_array);
    break;

  default:
    sollyaFprintf(stderr, "Error in autodiff: unknown binary operator (%d)\n", nodeType);
    return;
  }

  for(i=0;i<=n;i++) {
    sollya_mpfi_clear(res1[i]);
    sollya_mpfi_clear(res2[i]);
  }
  safeFree(res1);
  safeFree(res2);
}


/* Computes the successive derivatives of y -> y^p at point x          */
/* [x^p/0!    p*x^(p-1)/1!   ...   p*(p-1)*...*(p-n+1)*x^(p-n)/n! ]    */
void constantPower_diff(sollya_mpfi_t *res, sollya_mpfi_t x, mpfr_t p, int n, int *silent) {
  sollya_mpfi_t expo, acc;
  mp_prec_t prec_expo, prec;
  int i;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  prec_expo = (prec > mpfr_get_prec(p))?prec:mpfr_get_prec(p);

  sollya_mpfi_init2(expo, prec_expo);
  sollya_mpfi_init2(acc, prec);

  sollya_mpfi_set_fr(expo, p);
  sollya_mpfi_set_ui(acc, 1);

  for(i=0; i<=n; i++) {
    if (sollya_mpfi_is_zero(acc)) sollya_mpfi_set_ui(res[i],0);
    else {
      sollya_mpfi_pow(res[i], x, expo);
      sollya_mpfi_mul(res[i], res[i], acc);

      sollya_mpfi_mul(acc, acc, expo);
      sollya_mpfi_div_ui(acc, acc, i+1);
      sollya_mpfi_sub_ui(expo, expo, 1);
    }
  }

  sollya_mpfi_clear(expo);
  sollya_mpfi_clear(acc);

  return;
}


/* the power function is: p^x, where p is a positive constant */
/* [p^x/0!, log(p)p^x/1!, ... , log(p)^n p^x / n! ] */
void powerFunction_diff(sollya_mpfi_t *res, mpfr_t p, sollya_mpfi_t x, int n, int *silent) {
  int i;
  sollya_mpfi_t temp1,temp2;
  mp_prec_t prec;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  sollya_mpfi_init2(temp1, prec);
  sollya_mpfi_init2(temp2, prec);

  sollya_mpfi_set_fr(temp1,p);

  sollya_mpfi_pow(temp2, temp1, x); /* temp2 = p^x */
  sollya_mpfi_log(temp1,temp1); /* temp1 = log(p) */

  for(i=0;i<=n;i++) {
    sollya_mpfi_set(res[i], temp2);
    sollya_mpfi_mul(temp2,temp2,temp1);
    sollya_mpfi_div_ui(temp2, temp2, i+1);
  }

  sollya_mpfi_clear(temp1);
  sollya_mpfi_clear(temp2);
  return;
}


/* Takes a polynomial given by the array of its coefficients [p0...pn]
   and differentiates it (returns the array of the derivative) */

/* It IS safe to use the same pointer for res and coeff_array */
/* (in-place computation) */
void symbolic_poly_diff(sollya_mpfi_t *res, sollya_mpfi_t *coeff_array, int degree) {
  int i;

  for(i=0;i<=degree-1;i++) sollya_mpfi_mul_ui(res[i], coeff_array[i+1], i+1);
}

/* Evaluates a symbolic polynomial at point x by Horner scheme */
void symbolic_poly_evaluation_horner(sollya_mpfi_t res, sollya_mpfi_t *coeffs_array, sollya_mpfi_t x, int degree) {
  int i;
  sollya_mpfi_t temp;
  mp_prec_t prec;

  prec = getToolPrecision();
  sollya_mpfi_init2(temp, prec);

  sollya_mpfi_set(temp, coeffs_array[degree]);
  for(i=degree-1;i>=0;i--) {
    sollya_mpfi_mul(temp, temp, x);
    sollya_mpfi_add(temp, temp, coeffs_array[i]);
  }
  sollya_mpfi_set(res, temp);
  sollya_mpfi_clear(temp);
}

/* Evaluates a symbolic polynomial at point x by computing successive powers */
void symbolic_poly_evaluation_powers(sollya_mpfi_t res, sollya_mpfi_t *coeffs_array, sollya_mpfi_t *powers_array, sollya_mpfi_t x, int degree) {
  int i;
  sollya_mpfi_t temp, acc;
  mp_prec_t prec;

  UNUSED_PARAM(x);

  prec = getToolPrecision();
  sollya_mpfi_init2(temp, prec);
  sollya_mpfi_init2(acc, prec);

  sollya_mpfi_set_ui(acc, 0);
  for(i=0;i<=degree;i++) {
    sollya_mpfi_mul(temp, coeffs_array[i], powers_array[i]);
    sollya_mpfi_add(acc, acc, temp);
  }
  sollya_mpfi_set(res, acc);

  sollya_mpfi_clear(temp);
  sollya_mpfi_clear(acc);
}


void libraryFunction_diff(sollya_mpfi_t *res, node *f, sollya_mpfi_t x, int n, int *silent) {
  sollya_mpfi_t fact, temp, temp2;
  mp_prec_t prec;
  int i;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  sollya_mpfi_init2(fact, prec);
  sollya_mpfi_set_ui(fact, 1);
  sollya_mpfi_init2(temp, prec);

  for(i=0;i<=n;i++) {
    accessThruMemRef(f)->libFun->code(temp, x, accessThruMemRef(f)->libFunDeriv + i);
    sollya_init_and_convert_interval(temp2, temp);
    sollya_mpfi_div(res[i], temp2, fact);
    sollya_mpfi_clear(temp2);
    sollya_mpfi_mul_ui(fact, fact, i+1);
  }
  sollya_mpfi_clear(fact);
  sollya_mpfi_clear(temp);
}

void procedureFunction_diff(sollya_mpfi_t *res, node *f, sollya_mpfi_t x, int n, int *silent) {
  sollya_mpfi_t fact;
  mp_prec_t prec;
  int i;

  UNUSED_PARAM(silent);

  prec = getToolPrecision();
  sollya_mpfi_init2(fact, prec);
  sollya_mpfi_set_ui(fact, 1);

  for(i=0;i<=n;i++) {
    computeFunctionWithProcedure(res[i], accessThruMemRef(f)->child2, x, (unsigned int) (accessThruMemRef(f)->libFunDeriv + i));
    sollya_mpfi_div(res[i], res[i], fact);
    sollya_mpfi_mul_ui(fact, fact, i+1);
  }
  sollya_mpfi_clear(fact);
}


/* res is a reserved space for n+1 sollya_mpfi_t such that: */
/*               res_i = f^(i)(x0)/i!                */
/* We proceed recursively on the structure.          */
void auto_diff_scaled(sollya_mpfi_t* res, node *f, sollya_mpfi_t x0, int n) {
  int i;
  sollya_mpfi_t *res1, *res2;
  node *simplifiedChild1, *simplifiedChild2, *tempTree;
  sollya_mpfi_t temp1, temp2;
  mp_prec_t prec;
  int silent = 0;

  prec = getToolPrecision();
  switch (accessThruMemRef(f)->nodeType) {
  case VARIABLE:
    sollya_mpfi_set(res[0], x0);
    if(n>=1) {
      sollya_mpfi_set_ui(res[1], 1);
      for(i=2; i<=n; i++) sollya_mpfi_set_ui(res[i], 0);
    }
    break;

  case PI_CONST:
    sollya_mpfi_const_pi(res[0]);
    for(i=1; i<=n; i++) sollya_mpfi_set_ui(res[i], 0);
    break;

  case LIBRARYCONSTANT:
    libraryConstantToInterval(res[0], accessThruMemRef(f));
    for(i=1; i<=n; i++) sollya_mpfi_set_ui(res[i], 0);
    break;

  case CONSTANT:
    sollya_mpfi_set_fr(res[0], *(accessThruMemRef(f)->value));
    for(i=1; i<=n; i++) sollya_mpfi_set_ui(res[i], 0);
    break;

  case NEG:
    auto_diff_scaled(res, accessThruMemRef(f)->child1, x0, n);
    for(i=0;i<=n;i++) sollya_mpfi_neg(res[i], res[i]);
    break;

  case ADD:
  case SUB:
  case MUL:
  case DIV:
    binary_function_diff(res, accessThruMemRef(f)->nodeType, x0, accessThruMemRef(f)->child1, accessThruMemRef(f)->child2, n, &silent);
    break;

  case UNARY_BASE_FUNC:
  case LIBRARYFUNCTION:
  case PROCEDUREFUNCTION:
    res1 = (sollya_mpfi_t *)safeCalloc((n+1),sizeof(sollya_mpfi_t));
    res2 = (sollya_mpfi_t *)safeCalloc((n+1),sizeof(sollya_mpfi_t));
    for(i=0;i<=n;i++) {
      sollya_mpfi_init2(res1[i], prec);
      sollya_mpfi_init2(res2[i], prec);
    }

    auto_diff_scaled(res1, accessThruMemRef(f)->child1, x0, n);
    if(accessThruMemRef(f)->nodeType==LIBRARYFUNCTION) libraryFunction_diff(res2, accessThruMemRef(f), res1[0], n, &silent);
    else if(accessThruMemRef(f)->nodeType==PROCEDUREFUNCTION) procedureFunction_diff(res2, accessThruMemRef(f), res1[0], n, &silent);
    else accessThruMemRef(f)->baseFun->baseAutodiff(res2, res1[0], n, &silent);
    composition_AD(res, res2, res1, n);

    for(i=0;i<=n;i++) {
      sollya_mpfi_clear(res1[i]);
      sollya_mpfi_clear(res2[i]);
    }
    safeFree(res1);
    safeFree(res2);
    break;

  case POW:
    simplifiedChild2 = simplifyTreeErrorfree(accessThruMemRef(f)->child2);
    simplifiedChild1 = simplifyTreeErrorfree(accessThruMemRef(f)->child1);

    /* x^p case */
    if ( (accessThruMemRef(simplifiedChild1)->nodeType == VARIABLE) &&
	 (accessThruMemRef(simplifiedChild2)->nodeType == CONSTANT) ) {
      constantPower_diff(res, x0, *(accessThruMemRef(simplifiedChild2)->value), n, &silent);
    }

    /* p^x case */
    else if ( (accessThruMemRef(simplifiedChild1)->nodeType == CONSTANT) &&
	      (accessThruMemRef(simplifiedChild2)->nodeType == VARIABLE) ) {
      powerFunction_diff(res, *(accessThruMemRef(simplifiedChild1)->value), x0, n, &silent);
    }

    /* p^q case */
    else if ( (accessThruMemRef(simplifiedChild1)->nodeType == CONSTANT) &&
	      (accessThruMemRef(simplifiedChild2)->nodeType == CONSTANT) ) {
      sollya_mpfi_init2(temp1, prec);
      sollya_mpfi_set_fr(temp1, *(accessThruMemRef(simplifiedChild1)->value));
      sollya_mpfi_init2(temp2, prec);
      sollya_mpfi_set_fr(temp2, *(accessThruMemRef(simplifiedChild2)->value));
      sollya_mpfi_pow(res[0], temp1, temp2);
      for(i=1; i<=n; i++) sollya_mpfi_set_ui(res[i], 0);

      sollya_mpfi_clear(temp1);
      sollya_mpfi_clear(temp2);
    }

    /* p^f or f^p case */
    else if ( (accessThruMemRef(simplifiedChild1)->nodeType==CONSTANT) ||
	      (accessThruMemRef(simplifiedChild2)->nodeType==CONSTANT) ) {

      res1 = (sollya_mpfi_t *)safeCalloc((n+1),sizeof(sollya_mpfi_t));
      res2 = (sollya_mpfi_t *)safeCalloc((n+1),sizeof(sollya_mpfi_t));
      for(i=0;i<=n;i++) {
	sollya_mpfi_init2(res1[i], prec);
	sollya_mpfi_init2(res2[i], prec);
      }

      if (accessThruMemRef(simplifiedChild1)->nodeType == CONSTANT) { /* p^f */
	auto_diff_scaled(res1, simplifiedChild2, x0, n);
	powerFunction_diff(res2, *(accessThruMemRef(simplifiedChild1)->value), res1[0], n, &silent);
      }
      else { /* f^p */
	auto_diff_scaled(res1, simplifiedChild1, x0, n);
	constantPower_diff(res2, res1[0], *(accessThruMemRef(simplifiedChild2)->value), n, &silent);
      }

      composition_AD(res, res2, res1, n);

      for(i=0; i<=n; i++) {
	sollya_mpfi_clear(res1[i]);
	sollya_mpfi_clear(res2[i]);
      }
      safeFree(res1);
      safeFree(res2);
    }

    /*  f^g case */
    /* f^g = exp(g*log(f)) */
    else {
      tempTree = makeExp(makeMul(copyTree(simplifiedChild2), makeLog(copyTree(simplifiedChild1))));
      auto_diff_scaled(res, tempTree, x0, n);
      free_memory(tempTree);
    }

    free_memory(simplifiedChild1);
    free_memory(simplifiedChild2);
    break;

  default:
    sollyaFprintf(stderr,"Error in autodiff: unknown identifier (%d) in the tree\n",accessThruMemRef(f)->nodeType);
    exit(1);
  }

  return;
}

/* res is a reserved space for n+1 sollya_mpfi_t such that: */
/*               res_i = f^(i)(x0)                   */
void auto_diff(sollya_mpfi_t* res, node *f, sollya_mpfi_t x0, int n) {
  int i;
  sollya_mpfi_t fact;
  mp_prec_t prec;

  prec = getToolPrecision();

  sollya_mpfi_init2(fact, prec);
  sollya_mpfi_set_ui(fact, 1);

  auto_diff_scaled(res, f, x0, n);
  for(i=1;i<=n;i++) {
    sollya_mpfi_mul_ui(fact, fact, i);
    sollya_mpfi_mul(res[i], res[i], fact);
  }

  sollya_mpfi_clear(fact);
}
