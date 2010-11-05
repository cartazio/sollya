/*

Copyright 2010 by 

Laboratoire d'Informatique de Paris 6, equipe PEQUAN,
UPMC Universite Paris 06 - CNRS - UMR 7606 - LIP6, Paris, France.

and

Laboratoire de l'Informatique du Parallélisme, 
UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668,

Contributors Ch. Lauter, M. Joldes

christoph.lauter@ens-lyon.org
mioara.joldes@ens-lyon.fr

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

#include "supnorm.h"
#include <mpfr.h>
#include "mpfi-compat.h"
#include "execute.h"
#include <stdio.h> 
#include <stdlib.h>
#include "expression.h"
#include "infnorm.h"
#include "autodiff.h"
#include "taylorform.h"
#include "chain.h"
#include "sturm.h"
#include "general.h"
#include "infnorm.h"

/* Add error codes here as needed. 

   When adding error codes, add warning messages below (in
   supremumNormBisect()).

*/
#define SUPNORM_NO_ERROR                       0
#define SUPNORM_SOME_ERROR                    -1
#define SUPNORM_NO_TAYLOR                      1
#define SUPNORM_NOT_ENOUGH_WORKING_PRECISION   2
#define SUPNORM_SINGULARITY_NOT_REMOVED        3
#define SUPNORM_COULD_NOT_SHOW_POSITIVITY      4

/* General remarks:

   (i) We will not handle removable singularities inside the
   expression tree for func right now. However, as we might in the
   future, it might be good idea not to use the Taylor form functions
   from taylorform.h directly but to wrap them into a function, which,
   in the future, may take these singularities into account.

   (ii) Do not use getToolPrecision(). Use the precision given in
   argument to the functions and propagate it as needed. Or: perform
   exact operations that determine their precision on their own.

   (iii) List of functions that might come handy:

   - node* addPolynomialsExactly(node *p1, node *p2)
   - node* subPolynomialsExactly(node *p1, node *p2)
   - node* scalePolynomialExactly(node *poly, mpfr_t scale)
   - int showPositivity(node * poly, sollya_mpfi_t dom, mp_prec_t prec)
   - void evaluateInterval(sollya_mpfi_t y, node *func, node *deriv, sollya_mpfi_t x)
     (remark that deriv may be set to NULL if it is not known)
   - void uncertifiedInfnorm(mpfr_t result, node *tree, mpfr_t a, mpfr_t b, unsigned long int points, mp_prec_t prec)
     (that's dirtyinfnorm)
   - void taylorform(node **T, chain **errors, sollya_mpfi_t **delta,
		node *f, int n, sollya_mpfi_t *x0, sollya_mpfi_t *d, int mode)
   - ...

 */


/* Exactly add two polynomials

   For each monomial degree do:

   In the case when both monomials are floating-point numbers,
   return a monomial that is a floating-point number.

   In the case when one of the monomials is a ratio of two
   floating-point numbers and the other is a floating-point number or
   a ratio of floating-point numbers, return a ratio of floating-point
   numbers.

   Otherwise, return a monomial representing the sum of the two
   expressions.

   This version is a first, crude version. We can do better if need
   be.
*/
node *addPolynomialsExactly(node *p1, node *p2) {
  node *temp, *res, *temp2;

  if (!(isPolynomial(p1) && isPolynomial(p2))) {
    temp = makeAdd(copyTree(p1),copyTree(p2));
    res = simplifyTreeErrorfree(temp);
    free_memory(temp);
    return res;
  }

  /* Here, we know that we have polynomials */
  temp = makeAdd(copyTree(p1),copyTree(p2));
  temp2 = horner(temp);
  res = simplifyRationalErrorfree(temp2);
  free_memory(temp);
  free_memory(temp2);

  return res;
}

/* Exactly subtract two polynomials

   For each monomial degree do:

   In the case when both monomials are floating-point numbers,
   return a monomial that is a floating-point number.

   In the case when one of the monomials is a ratio of two
   floating-point numbers and the other is a floating-point number or
   a ratio of floating-point numbers, return a ratio of floating-point
   numbers.

   Otherwise, return a monomial representing the sum of the two
   expressions.

   This version is a first, crude version. We can do better if need
   be.
*/
node *subPolynomialsExactly(node *p1, node *p2) {
  node *temp, *res, *temp2;

  if (!(isPolynomial(p1) && isPolynomial(p2))) {
    temp = makeSub(copyTree(p1),copyTree(p2));
    res = simplifyTreeErrorfree(temp);
    free_memory(temp);
    return res;
  }

  /* Here, we know that we have polynomials */
  temp = makeSub(copyTree(p1),copyTree(p2));
  temp2 = horner(temp);
  res = simplifyRationalErrorfree(temp2);
  free_memory(temp);
  free_memory(temp2);

  return res;
}

/* Exactly scale a polynomials

   Return sum (s * c_i) * x^i for s and p(x) = sum c_i * x^i 

   Wherever possible (even when more precision is needed), simplify
   the expressions s * c_i in the monomials.

   This version is a first, crude version. We can do better if need
   be.
*/
node *scalePolynomialExactly(node *poly, mpfr_t scale) {
  node *temp, *res, *temp2;

  if (!(isPolynomial(poly))) {
    temp = makeMul(copyTree(poly),makeConstant(scale));
    res = simplifyTreeErrorfree(temp);
    free_memory(temp);
    return res;
  }

  /* Here, we know that we have a polynomial */
  temp = makeMul(copyTree(poly),makeConstant(scale));
  temp2 = horner(temp);
  res = simplifyRationalErrorfree(temp2);
  free_memory(temp);
  free_memory(temp2);

  return res;
}


/* Show positivity of a polynomial using the Sturm algorithm 

   For a polynomial poly and a domain dom,

   return zero if there is a point in dom at which poly is non-positive,
   return a non-zero value otherwise.

   Return zero for expressions that are no polynomial.
   Return zero for domains that are not bounded by finite numbers.

   Return zero if something goes wrong in the Sturm routine.

   Generally, return zero in case of a doubt (not enough precision
   etc.)

*/
int showPositivity(node * poly, sollya_mpfi_t dom, mp_prec_t prec) {
  int res, positive, nbRoots;
  mpfr_t a, b, c, y, nbRootsMpfr;
  mp_prec_t pp;

  if (!isPolynomial(poly)) return 0;
  if (!sollya_mpfi_bounded_p(dom)) return 0;

  mpfr_init2(nbRootsMpfr,8 * sizeof(int));
  res = getNrRoots(nbRootsMpfr, poly, dom, prec);
  if (!mpfr_number_p(nbRootsMpfr)) {
    nbRoots = 1;
  } else {
    nbRoots = mpfr_get_si(nbRootsMpfr,GMP_RNDU);
  }
  mpfr_clear(nbRootsMpfr);
  if (!res) return 0;

  if (nbRoots != 0) return 0;

  /* Here, we know that we do not cross the abscissa

     We still have to establish that at some point (in the middle of
     the interval), we have a positive value.

  */
  pp = sollya_mpfi_get_prec(dom);
  mpfr_init2(a,pp);
  mpfr_init2(b,pp);
  mpfr_init2(c,pp+1);

  sollya_mpfi_get_left(a,dom);
  sollya_mpfi_get_right(b,dom);

  mpfr_add(c,a,b,GMP_RNDN);
  mpfr_div_2ui(c,c,1,GMP_RNDN);

  mpfr_init2(y,16);
  res = evaluateFaithful(y, poly, c, prec);

  positive = 1;
  if (!res) positive = 0;
  if (!mpfr_number_p(y)) positive = 0;
  if (mpfr_sgn(y) <= 0) positive = 0;

  mpfr_clear(a);
  mpfr_clear(b);
  mpfr_clear(c);
  mpfr_clear(y);

  return positive;
}


/* Compute the supremum norm on eps = p - f over dom

   The supremum norm is computed with an enclosure error less than accuracy,
   i.e. the obtained interval [l,u] satisfies in the end:

   abs(u-l/l) <= accuracy

   If everything works fine, result is affected with an interval
   safely enclosing the supremum norm. The return value is then ZERO
   (i.e. SUPNORM_NO_ERROR).

   Otherwise, if an error occurs, the return value is non-zero. The
   number of the return value then corresponds to a special error
   message meaning (see #defines above).

   No warning message is ever displayed by this function.

   The computing precision is prec.

   We are ensured that this function is called only if 

   * poly is a polynomial
   * dom is a closed non-empty interval containing only numbers that is not reduced to a point,
   * accuracy is a positive number.

   In the case when the computation fails but there is hope in
   obtaining a result by bisection, the algorithm may assign a point
   in the interior of dom to bisectPoint. The global bisection code
   will then try to bisect at that point. If no value is assigned, the
   global bisection will be performed at the midpoint of dom. This
   means: if you just want default behavior for the bisection (in the
   midpoint), then do not touch bisectPoint.

   To begin with, we do not care about removable singularities in the
   expression of func. In such a case, the absolute supnorm may simply
   fail for now.

*/
int supnormAbsolute(sollya_mpfi_t result, node *poly, node *func, sollya_mpfi_t dom, mpfr_t accuracy, mp_prec_t prec, mpfr_t bisectPoint) {

  /* TODO */

  return SUPNORM_SOME_ERROR;
}

/* Compute the supremum norm on eps = p/f - 1 over dom

   The supremum norm is computed with an enclosure error less than accuracy,
   i.e. the obtained interval [l,u] satisfies in the end:

   abs(u-l/l) <= accuracy

   If everything works fine, result is affected with an interval
   safely enclosing the supremum norm. The return value is then ZERO
   (i.e. SUPNORM_NO_ERROR).

   Otherwise, if an error occurs, the return value is non-zero. The
   number of the return value then corresponds to a special error
   message meaning (see #defines above).

   No warning message is ever displayed by this function.

   The computing precision is prec.

   We are ensured that this function is called only if 

   * poly is a polynomial
   * dom is a closed non-empty interval containing only numbers that is not reduced to a point,
   * accuracy is a positive number.

   In the case when the computation fails but there is hope in
   obtaining a result by bisection, the algorithm may assign a point
   in the interior of dom to bisectPoint. The global bisection code
   will then try to bisect at that point. If no value is assigned, the
   global bisection will be performed at the midpoint of dom. This
   means: if you just want default behavior for the bisection (in the
   midpoint), then do not touch bisectPoint.

   This function is supposed to detect and overcome false
   singularities. However, it is also supposed to perform a fast check
   first, i.e. it is not supposed to do a length detection of zeros of
   func if the (IA) image of func over dom does not contain zero. The
   use of void evaluateInterval(sollya_mpfi_t y, node *func, node
   *deriv, sollya_mpfi_t x); comes handy here (remark that deriv may
   be set to NULL).

   To begin with, we do not care about removable singularities in the
   expression of func. In such a case, the relative supnorm may simply
   fail for now.

   However, removable singularities in poly/func must be detected
   (after a fast check if there aren't any) and overcome. This must
   work for cases when there is one singularity in dom, though. If
   there are several singularities, bisection will eventually split
   the interval. If appropriate set bisectPoint to something
   reasonable. 

*/
int supnormRelative(sollya_mpfi_t result, node *poly, node *func, sollya_mpfi_t dom, mpfr_t accuracy, mp_prec_t prec, mpfr_t bisectPoint) {

  /* TODO */

  return SUPNORM_SOME_ERROR;
}

/* Compute the supremum norm on eps = p - f resp. eps = p/f - 1 over dom

   eps is defined according to the mode parameter:
   if mode = 0 then eps = p - f else eps = p/f -1

   The supremum norm is computed with an enclosure error less than accuracy,
   i.e. the obtained interval [l,u] satisfies in the end:

   abs(u-l/l) <= accuracy

   If everything works fine, result is affected with an interval
   safely enclosing the supremum norm. The return value is then ZERO
   (i.e. SUPNORM_NO_ERROR).

   Otherwise, if an error occurs, the return value is non-zero. The
   number of the return value then corresponds to a special error
   message meaning (see #defines above).

   No warning message is ever displayed by this function.

   The computing precision is prec.

   We are ensured that this function is called only if 

   * poly is a polynomial
   * dom is a closed non-empty interval containing only numbers that is not reduced to a point,
   * accuracy is a positive number.

   In the case when the computation fails but there is hope in
   obtaining a result by bisection, the algorithm may assign a point
   in the interior of dom to bisectPoint. The global bisection code
   will then try to bisect at that point. If no value is assigned, the
   global bisection will be performed at the midpoint of dom. This
   means: if you just want default behavior for the bisection (in the
   midpoint), then do not touch bisectPoint.

*/
int supremumNormInner(sollya_mpfi_t result, node *poly, node *func, sollya_mpfi_t dom, int mode, mpfr_t accuracy, mp_prec_t prec, mpfr_t bisectPoint) {
  int res;

  if (mode == 0) {
    res = supnormAbsolute(result,poly,func,dom,accuracy,prec,bisectPoint);
  } else {
    res = supnormRelative(result,poly,func,dom,accuracy,prec,bisectPoint);
  }

  return res;
}

/* Compute the supremum norm on eps = p - f resp. eps = p/f - 1 over [a,b]

   eps is defined according to the mode parameter:
   if mode = 0 then eps = p - f else eps = p/f -1

   The supremum norm is computed with an enclosure error less than accuracy,
   i.e. the obtained interval [l,u] satisfies in the end:

   abs(u-l/l) <= accuracy

   If everything works fine, result is affected with an interval
   safely enclosing the supremum norm. The return value is then ZERO
   (i.e. SUPNORM_NO_ERROR).

   Otherwise, if an error occurs, the return value reflects the last
   error message of the recursive calls.

   We are ensured that this function is called only if 

   * poly is a polynomial
   * a and b, a < b, form an interval [a,b] that is not reduced to a point
   * accuracy is a positive number,
   * diameter is a non-negative number.

   The algorithm is allowed to stop if a result has been found or if
   abs(b-a) is less than diameter.

*/
int supremumNormBisectInner(sollya_mpfi_t result, node *poly, node *func, mpfr_t a, mpfr_t b, int mode, mpfr_t accuracy, mpfr_t diameter, mp_prec_t prec) {
  int resFirst, res, resLeft, resRight;
  mp_prec_t p1, p2;
  sollya_mpfi_t dom, resultLeft, resultRight;
  mpfr_t c, width;
  mpfr_t ll, lr, ul, ur;

  p1 = mpfr_get_prec(a);
  p2 = mpfr_get_prec(b);
  if (p2 > p1) p1 = p2;
  sollya_mpfi_init2(dom,p1);
  sollya_mpfi_interv_fr(dom,a,b);

  mpfr_init2(c, p1 + 1);
  mpfr_add(c,a,b,GMP_RNDN);
  mpfr_div_2ui(c,c,1,GMP_RNDN);

  /* Call inner supnorm algorithm on whole interval */
  resFirst = supremumNormInner(result, poly, func, dom, mode, accuracy, prec, c);

  sollya_mpfi_clear(dom);

  /* If everything worked fine on the whole interval, we do not need to bisect */
  if (resFirst == SUPNORM_NO_ERROR) {
    mpfr_clear(c);
    return SUPNORM_NO_ERROR;
  }

  /* Some error occured, check if we still need to bisect or if we have reached diameter */
  mpfr_init2(width,p1);
  mpfr_sub(width,b,a,GMP_RNDU);
  if (mpfr_cmp(width,diameter) < 0) {
    /* We have reached diameter, we return the latest error code */
    mpfr_clear(c);
    mpfr_clear(width);
    return resFirst;
  }

  mpfr_clear(width);

  /* Here, we have to bisect.

     We check that the bisection point is a number in the interior of [a,b].
     If it is not, we set it to the middle of [a,b].

  */
  if ((!mpfr_number_p(c)) ||
      (mpfr_cmp(c,a) <= 0) ||
      (mpfr_cmp(b,c) <= 0)) {
    mpfr_add(c,a,b,GMP_RNDN);
    mpfr_div_2ui(c,c,1,GMP_RNDN);
  }

  /* Bisect */
  p2 = sollya_mpfi_get_prec(result);
  sollya_mpfi_init2(resultLeft,p2);

  resLeft = supremumNormBisectInner(resultLeft, poly, func, a, c, mode, accuracy, diameter, prec);

  if (resLeft != SUPNORM_NO_ERROR) {
    /* The bisection recursively failed on the left sub-interval [a,c] 
       Return the error code returned by that recursive call.
     */
    mpfr_clear(c);
    sollya_mpfi_clear(resultLeft);
    return resLeft;
  }
  /* Here, resLeft == SUPNORM_NO_ERROR */

  sollya_mpfi_init2(resultRight,p2);

  resRight = supremumNormBisectInner(resultRight, poly, func, c, b, mode, accuracy, diameter, prec);

  if (resRight != SUPNORM_NO_ERROR) {
    /* The bisection recursively failed on the right sub-interval [c,b] 
       Return the error code returned by that recursive call.
     */
    mpfr_clear(c);
    sollya_mpfi_clear(resultLeft);
    sollya_mpfi_clear(resultRight);
    return resRight;
  }
  
  /* Here, resLeft == SUPNORM_NO_ERROR and resRight == SUPNORM_NO_ERROR

     This means both recursive calls worked without error. 
     
     We combine the results by taking the max of the lower resp. the
     upper bounds.
     
   */
  mpfr_init2(ll,p2);
  mpfr_init2(lr,p2);
  mpfr_init2(ul,p2);
  mpfr_init2(ur,p2);
  
  sollya_mpfi_get_left(ll,resultLeft);
  sollya_mpfi_get_right(ul,resultLeft);
  sollya_mpfi_get_left(lr,resultRight);
  sollya_mpfi_get_right(ur,resultRight);

  if (mpfr_cmp(ll,lr) > 0) mpfr_set(lr,ll,GMP_RNDN); /* exact */
  if (mpfr_cmp(ul,ur) > 0) mpfr_set(ur,ul,GMP_RNDN); /* exact */

  sollya_mpfi_interv_fr(result,lr,ur);

  mpfr_clear(ll);
  mpfr_clear(lr);
  mpfr_clear(ul);
  mpfr_clear(ur);
  sollya_mpfi_clear(resultLeft);
  sollya_mpfi_clear(resultRight);
  mpfr_clear(c);
  return SUPNORM_NO_ERROR;
}

/* Compute the supremum norm on eps = p - f resp. eps = p/f - 1 over [a,b]

   eps is defined according to the mode parameter:
   if mode = 0 then eps = p - f else eps = p/f -1

   The supremum norm is computed with an enclosure error less than accuracy,
   i.e. the obtained interval [l,u] satisfies in the end:

   abs(u-l/l) <= accuracy

   If everything works fine, result is affected with an interval 
   safely enclosing the supremum norm. The return value is then non-zero.

   Otherwise, if an error occurs, the return value is 0.

   We are ensured that this function is called only if 

   * poly is a polynomial
   * a and b, a < b, form an interval [a,b] that is not reduced to a point
   * accuracy is a positive number,
   * diameter is a non-negative number.

   The algorithm is allowed to stop if a result has been found or if
   abs(b-a) is less than diameter.

*/
int supremumNormBisect(sollya_mpfi_t result, node *poly, node *func, mpfr_t a, mpfr_t b, int mode, mpfr_t accuracy, mpfr_t diameter) {
  int res;
  mp_prec_t prec, p;

  prec = getToolPrecision() + 25;
  p = mpfi_get_prec(result);
  if (p > prec) prec = p;

  res = supremumNormBisectInner(result, poly, func, a, b, mode, accuracy, diameter, prec);

  if (res == 0) return 1; /* everything's fine */
 
  /* In the following, perform error handling (messaging and return 1) */
  switch (res) {
  case SUPNORM_NO_TAYLOR:
    printMessage(1,"Warning: during supnorm computation, no suitable Taylor form could be found.\n");
    break;
  case SUPNORM_NOT_ENOUGH_WORKING_PRECISION:
    printMessage(1,"Warning: during supnorm computation, no result could be found as the working precision seems to be too low.\n");
    break;
  case SUPNORM_SINGULARITY_NOT_REMOVED:
    printMessage(1,"Warning: during supnorm computation, a singularity in the expression tree could not be removed.\n");
    break;
  case SUPNORM_COULD_NOT_SHOW_POSITIVITY:
    printMessage(1,"Warning: during supnorm computation, the positivity of a polynomial could not be established.\n");
    break;
  default:
    printMessage(1,"Warning: during supnorm computation, some generic error occured. No further description is available.\n");
  }

  return 0;
}

/* Assigns [NaN;NaN] to an interval 
*/
void assignNaNInterval(sollya_mpfi_t rop) {
  mpfr_t temp;
  
  mpfr_init2(temp,12);
  mpfr_set_nan(temp);
  sollya_mpfi_interv_fr(rop,temp,temp);
  mpfr_clear(temp);
}

/* Compute the supremum norm on eps = p - f resp. eps = p/f - 1 at a,
   i.e. evaluate abs(eps) at a.

   eps is defined according to the mode parameter:
   if mode = 0 then eps = p - f else eps = p/f -1

   The supremum norm is computed with an enclosure error less than accuracy,
   i.e. the obtained interval [l,u] satisfies in the end:

   abs(u-l/l) <= accuracy

   If everything works fine, result is affected with an interval 
   safely enclosing the supremum norm. The return value is then non-zero.

   Otherwise, if an error occurs, the return value is 0.

   We are ensured that this function is called only if 

   * poly is a polynomial
   * a is a number,
   * accuracy is a positive number,
   * diameter is a non-negative number.

*/
int supremumNormDegenerate(sollya_mpfi_t result, node *poly, node *func, mpfr_t a, int mode, mpfr_t accuracy) {
  node *absEps;
  int res;
  mpfr_t temp, y, ya, yb;
  unsigned int pr;
  mp_prec_t prec, pp;
  int tempRes;

  if (mode == 0) {
    /* Construct absEps = abs(poly - func) */
    absEps = makeAbs(makeSub(copyTree(poly),copyTree(func)));
  } else {
    /* Construct absEps = abs(poly/func - 1) */
    absEps = makeAbs(makeSub(makeDiv(copyTree(poly),copyTree(func)),makeConstantDouble(1.0)));
  }

  /* Compute pr = -floor(log2(accuracy)) to get the number of bits we need
     in the end
  */
  
  mpfr_init2(temp, 8 * sizeof(mp_prec_t) + 10);
  mpfr_log2(temp,accuracy,GMP_RNDD);
  mpfr_floor(temp,temp);
  mpfr_neg(temp,temp,GMP_RNDU);
  pr = mpfr_get_ui(temp,GMP_RNDD);
  mpfr_clear(temp);

  prec = getToolPrecision();
  if (pr > 2048 * prec) {
    printMessage(1,"Warning: the given accuracy depasses the current maximum precision of %d bits.\n",2048 * prec);
    printMessage(1,"Try to increase the precision of the tool.\n");
    assignNaNInterval(result);
    free_memory(absEps);
    return 0;
  }

  if (pr < prec) pp = prec; else pp = prec;
  pp += 10;
 
  mpfr_init2(y,pp);

  tempRes = evaluateFaithful(y, absEps, a, pp);

  res = 0;
  if (tempRes == 1) {
    res = 1;
    pp = mpfr_get_prec(y) - 5;
    mpfr_init2(ya,pp);
    mpfr_init2(yb,pp); 
    mpfr_set(ya,y,GMP_RNDD);
    mpfr_set(yb,y,GMP_RNDU);
    mpfr_nextbelow(ya);
    mpfr_nextbelow(ya);
    mpfr_nextabove(yb);
    mpfr_nextabove(yb);
    if (mpfr_sgn(ya) < 0) {
      mpfr_set_si(ya,0,GMP_RNDN);
    }

    sollya_mpfi_interv_fr(result,ya,yb);

    mpfr_clear(ya);
    mpfr_clear(yb);
  } else {
    printMessage(1,"Warning: could not perform a faithful evaluation of the error function between the given polynomial and the given function at the given point.\n");
    assignNaNInterval(result);
  }

  free_memory(absEps);
  mpfr_clear(y);

  return res;
}


/* Checks if all coefficients of poly can be written as ratios of
   floating-point numbers 

   Returns 0 if poly is not a polynomial or is a polynomial that does
   contain irrational coefficients.

   Returns a non-zero value otherwise.

*/
int hasOnlyMpqCoefficients(node *poly) {
  node **coefficients;
  int degree, i, res, okay;
  node *simplified;

  if (!isPolynomial(poly)) return 0;

  getCoefficients(&degree,&coefficients,poly);
  if (degree < 0) return 0;

  res = 1;
  for (i=0;i<=degree;i++) {
    if (coefficients[i] != NULL) {
      simplified = simplifyRationalErrorfree(coefficients[i]);
      okay = 0;
      if ((simplified->nodeType == CONSTANT) &&
	  (mpfr_number_p(*(simplified->value)))) okay = 1;
      else {
	if ((simplified->nodeType == DIV) &&
	    (((simplified->child1->nodeType == CONSTANT) && 
	      (mpfr_number_p(*(simplified->child1->value)))) && 
	     ((simplified->child2->nodeType == CONSTANT) && 
	      (mpfr_number_p(*(simplified->child2->value)))))) okay = 1;
      }
      free_memory(simplified);
      if (!okay) {
	res = 0;
	break;
      }
    }
  }

  for (i=0;i<=degree;i++) {
    if (coefficients[i] != NULL) free_memory(coefficients[i]);
  }
  free(coefficients);

  return res;
}



/* Compute the supremum norm on eps = p - f resp. eps = p/f - 1 over dom

   eps is defined according to the mode parameter:
   if mode = 0 then eps = p - f else eps = p/f -1

   The supremum norm is computed with an enclosure error less than accuracy,
   i.e. the interval [l,u] obtained satisfies in the end:

   abs(u-l/l) <= abs(accuracy)

   If everything works fine, result is affected with an interval 
   safely enclosing the supremum norm. The return value is then non-zero.

   Otherwise, if an error occurs, the return value is 0 and result is
   affected with [NaN;NaN]. A warning is printed in this case.
*/
int supremumnorm(sollya_mpfi_t result, node *poly, node *func, sollya_mpfi_t dom, int mode, mpfr_t accuracy) {
  mpfr_t a, b, diameter, temp, absAccuracy;
  mp_prec_t tempPrec;
  int res;

  if (!isPolynomial(poly)) {
    printMessage(1,"Warning: the given expression is not a polynomial.\n");
    assignNaNInterval(result);
    return 0;
  }

  tempPrec = sollya_mpfi_get_prec(dom);
  mpfr_init2(a,tempPrec);
  mpfr_init2(b,tempPrec);
  sollya_mpfi_get_left(a,dom);
  sollya_mpfi_get_right(b,dom);

  if (!(mpfr_number_p(a) && mpfr_number_p(b))) {
    printMessage(1,"Warning: the given domain is not a closed interval on the reals.\n");
    assignNaNInterval(result);
    mpfr_clear(a);
    mpfr_clear(b);
    return 0;
  }

  if (mpfr_cmp(a,b) > 0) {
    printMessage(1,"Warning: the given domain is empty.\n");
    assignNaNInterval(result);
    mpfr_clear(a);
    mpfr_clear(b);
    return 0;
  }

  if (mpfr_cmp(a,b) == 0) {
    printMessage(1,"Warning: the given domain is reduced to a point. Replacing the supremum norm with an evaluation.\n");
    res = supremumNormDegenerate(result,poly,func,a,mode,accuracy);
    if (!res) {
      printMessage(1,"Warning: could not evaluate the error function between the given polynomial and the given function at this point.\n");
      assignNaNInterval(result);
    } 
    mpfr_clear(a);
    mpfr_clear(b);
    return 1;
  }

  if (!mpfr_number_p(accuracy)) {
    printMessage(1,"Warning: the given accuracy is not a real number.\n");
    assignNaNInterval(result);
    mpfr_clear(a);
    mpfr_clear(b);
    return 0;
  }

  if (mpfr_zero_p(accuracy)) {
    printMessage(1,"Warning: the given accuracy is zero. In order to ensure the termination of the supremum norm algorithm, the accuracy parameter must be non-zero.\n");
    assignNaNInterval(result);
    mpfr_clear(a);
    mpfr_clear(b);
    return 0;
  }

  if (!hasOnlyMpqCoefficients(poly)) {
    printMessage(1,"Warning: the coefficients of the given polynomial cannot all be written as ratios of floating-point numbers.\nSupremum norm computation is only possible on such polynomials. Try to use roundcoefficients().\n");
    assignNaNInterval(result);
    mpfr_clear(a);
    mpfr_clear(b);
    return 0;
  }

  /* Here, we know that the interval is proper (no NaNs, no Infs) and
     that it is not reduced to a point. We know that accuracy is a
     non-zero number and that poly is a polynomial whose coefficients 
     can be written in floating-point numbers or ratios of floating-point
     numbers.

     We will call supremumNormBisect with diam * width(dom) and abs(accuracy).

  */

  mpfr_init2(temp,tempPrec * 4);
  mpfr_init2(diameter, tempPrec * 4 + 53);
  mpfr_sub(temp,b,a,GMP_RNDU);
  getToolDiameter(diameter);
  mpfr_mul(diameter,diameter,temp,GMP_RNDU);
  mpfr_abs(diameter,diameter,GMP_RNDN);

  mpfr_init2(absAccuracy,mpfr_get_prec(accuracy));
  mpfr_abs(absAccuracy,accuracy,GMP_RNDN); /* exact */

  res = supremumNormBisect(result,poly,func,a,b,mode,absAccuracy,diameter);

  if (!res) {
    printMessage(1,"Warning: an error occured during supremum norm computation. A safe enclosure of the supremum norm could not be computed.\n");
    assignNaNInterval(result);
  }  

  mpfr_clear(a);
  mpfr_clear(b);
  mpfr_clear(temp);
  mpfr_clear(diameter);
  mpfr_clear(absAccuracy);

  return res;
}




