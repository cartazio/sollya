/*

  Copyright 2014 by

  Laboratoire d'Informatique de Paris 6, equipe PEQUAN,
  UPMC Universite Paris 06 - CNRS - UMR 7606 - LIP6, Paris, France

  Contributor Ch. Lauter

  christoph.lauter@lip6.fr

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
#include "execute.h"
#include "infnorm.h"
#include "sturm.h"
#include "hooks.h"

int addEvaluationHook(eval_hook_t **hookPtr, 
		      void *data, 
		      int (*evaluateFunc)(sollya_mpfi_t, sollya_mpfi_t, mp_prec_t, void *), 
		      void (*freeFunc)(void *),
		      int (*compareFunc)(void *, void *),
		      void *(*copyFunc)(void *)) {
  eval_hook_t *newHook, *curr;

  /* Check if this hook has already been installed. If yes, deallocate
     the hook we have been given.
  */
  for (curr=*hookPtr;curr!=NULL;curr=curr->nextHook) {
    if (((curr->evaluateHook == evaluateFunc) &&
	 (curr->freeHook == freeFunc) &&
	 (curr->compareHook == compareFunc)) &&
	curr->compareHook(curr->data, data)) {
      freeFunc(data);
      return 0;
    }
  }

  /* Here, we know that the new hook does not yet exist 
     
     Create a new hook.
  */
  newHook = (eval_hook_t *) safeMalloc(sizeof(eval_hook_t));
  newHook->data = data;
  newHook->evaluateHook = evaluateFunc;
  newHook->freeHook = freeFunc;
  newHook->compareHook = compareFunc;
  newHook->copyHook = copyFunc;
  newHook->nextHook = *hookPtr;

  /* Assign the new hook */
  *hookPtr = newHook;
  
  /* Signal success */
  return 1;
}

int addEvaluationHookFromCopy(eval_hook_t **newHookPtr, eval_hook_t *hook) {
  int res, r;
  eval_hook_t *curr;

  res = 0;
  for (curr=hook; curr != NULL; curr=curr->nextHook) {
    r = addEvaluationHook(newHookPtr, 
			  curr->copyHook(curr->data),
			  curr->evaluateHook,
			  curr->freeHook,
			  curr->compareHook,
			  curr->copyHook);
    res = res || r;
  }
  return res;
}

void freeEvaluationHook(eval_hook_t **hookPtr) {
  eval_hook_t *curr, *next;
  
  curr = *hookPtr;
  while (curr != NULL) {
    next = curr->nextHook;
    curr->freeHook(curr->data);
    safeFree(curr);
    curr = next;
  }
  *hookPtr = NULL;
}

int evaluateWithEvaluationHook(sollya_mpfi_t y, sollya_mpfi_t x, mp_prec_t prec, eval_hook_t *hook) {
  eval_hook_t *curr;

  if (hook == NULL) return 0;

  for (curr=hook;curr!=NULL;curr=curr->nextHook) {
    if (curr->evaluateHook(y,x,prec,curr->data)) {
      return 1;
    }
  }
  
  return 0;
}

node_eval_hook_t *createNodeEvalHook(node *func, sollya_mpfi_t dom, sollya_mpfi_t delta, sollya_mpfi_t t) {
  node_eval_hook_t *newNodeEvalHook;

  newNodeEvalHook = (node_eval_hook_t *) safeMalloc(sizeof(node_eval_hook_t));
  sollya_mpfi_init2(newNodeEvalHook->domain, sollya_mpfi_get_prec(dom));
  sollya_mpfi_set(newNodeEvalHook->domain, dom);
  sollya_mpfi_init2(newNodeEvalHook->delta, sollya_mpfi_get_prec(delta));
  sollya_mpfi_set(newNodeEvalHook->delta, delta);
  sollya_mpfi_init2(newNodeEvalHook->t, sollya_mpfi_get_prec(t));
  sollya_mpfi_set(newNodeEvalHook->t, t);
  newNodeEvalHook->func = copyThing(func);

  return newNodeEvalHook;
}


int evaluateNodeEvalHook(sollya_mpfi_t y, sollya_mpfi_t x, mp_prec_t prec, void *data) {
  node_eval_hook_t *hook;
  mp_prec_t p, pY, pX;
  sollya_mpfi_t myY, myYRnd, myYRndWithDelta, redX;
  int okay;

  hook = (node_eval_hook_t *) data;

  if (sollya_mpfi_has_nan(x)) return 0;
  if (sollya_mpfi_has_infinity(x)) return 0;
  if (!sollya_mpfi_is_inside(x, hook->domain)) return 0;
  
  pY = sollya_mpfi_get_prec(y); 
  pX = sollya_mpfi_get_prec(x); 
  p = pY + 10;
  if (prec > p) p = prec;

  sollya_mpfi_init2(myY, p);
  sollya_mpfi_init2(redX, (p > pX ? p : pX));
  sollya_mpfi_sub(redX, x, hook->t);
  evaluateInterval(myY, hook->func, NULL, redX);

  okay = 0;
  sollya_mpfi_init2(myYRnd, pY + 5);
  sollya_mpfi_init2(myYRndWithDelta, pY + 5);

  sollya_mpfi_set(myYRnd, myY);
  sollya_mpfi_blow_1ulp(myYRnd);
  sollya_mpfi_add(myYRndWithDelta, myY, hook->delta);

  if (sollya_mpfi_is_inside(myYRndWithDelta, myYRnd) && 
      (!(sollya_mpfi_has_nan(myYRndWithDelta) || 
	 (sollya_mpfi_has_infinity(myYRndWithDelta) && 
	  (!sollya_mpfi_is_infinity(myYRndWithDelta)))))) okay = 1;
  
  if (okay) sollya_mpfi_set(y, myYRndWithDelta);

  sollya_mpfi_clear(myYRnd);
  sollya_mpfi_clear(myYRndWithDelta);
  sollya_mpfi_clear(myY);
  sollya_mpfi_clear(redX);

  return okay;
}

void freeNodeEvalHook(void *data) {
  node_eval_hook_t *hook;

  hook = (node_eval_hook_t *) data;
  freeThing(hook->func);
  sollya_mpfi_clear(hook->domain);
  sollya_mpfi_clear(hook->delta);
  sollya_mpfi_clear(hook->t);
  safeFree(hook);
}

void *copyNodeEvalHook(void *data) {
  node_eval_hook_t *hook;

  hook = (node_eval_hook_t *) data;
  return (void *) createNodeEvalHook(hook->func, hook->domain, hook->delta, hook->t);
}

int compareNodeEvalHook(void *data1, void *data2) {
  node_eval_hook_t *hook1, *hook2;

  hook1 = (node_eval_hook_t *) data1;
  hook2 = (node_eval_hook_t *) data2;

  if (!sollya_mpfi_equal_p(hook1->domain, hook2->domain)) return 0;
  if (!sollya_mpfi_equal_p(hook1->delta, hook2->delta)) return 0;
  if (!sollya_mpfi_equal_p(hook1->t, hook2->t)) return 0;
  if (!isEqualThing(hook1->func, hook2->func)) return 0;

  return 1;
}

int evaluatePolyEvalHook(sollya_mpfi_t, sollya_mpfi_t, mp_prec_t, void *);

poly_eval_hook_t *createPolyEvalHook(int degree, mpfr_t *coeffs, sollya_mpfi_t dom, sollya_mpfi_t delta, sollya_mpfi_t t) {
  poly_eval_hook_t *newPolyEvalHook;
  mpfr_t *derivCoeffs;
  int i, derivPolyCoeffsOkay;
  node *derivPoly;
  sollya_mpfi_t X, Y;
  mpfr_t ymax, deltamax;
  mpfr_t nbRoots;
  mp_prec_t maxPrec;
  unsigned long pp, ppp;
  mp_prec_t p;

  /* Allocate memory for the hook */
  newPolyEvalHook = (poly_eval_hook_t *) safeMalloc(sizeof(poly_eval_hook_t));

  /* Set some technical stuff */
  newPolyEvalHook->reusedVarMyYInit = 0;
  newPolyEvalHook->reusedVarXInit = 0;
  newPolyEvalHook->reusedVarTempInit = 0;
  newPolyEvalHook->reusedVarMyYBInit = 0;
  newPolyEvalHook->reusedVarXAInit = 0;
  newPolyEvalHook->reusedVarXBInit = 0;
  newPolyEvalHook->reusedVarMyYRndInit = 0;
  newPolyEvalHook->reusedVarMyYRndWithDeltaInit = 0;
  newPolyEvalHook->reusedVarAInit = 0;
  newPolyEvalHook->reusedVarBInit = 0;

  /* Set the easy stuff */
  sollya_mpfi_init2(newPolyEvalHook->domain, sollya_mpfi_get_prec(dom));
  sollya_mpfi_set(newPolyEvalHook->domain, dom);
  sollya_mpfi_init2(newPolyEvalHook->delta, sollya_mpfi_get_prec(delta));
  sollya_mpfi_set(newPolyEvalHook->delta, delta);
  sollya_mpfi_init2(newPolyEvalHook->t, sollya_mpfi_get_prec(t));
  sollya_mpfi_set(newPolyEvalHook->t, t);
  newPolyEvalHook->degree = degree;
  newPolyEvalHook->coefficients = (mpfr_t *) safeCalloc(degree + 1, sizeof(mpfr_t));
  for (i=0;i<=degree;i++) {
    mpfr_init2(newPolyEvalHook->coefficients[i], mpfr_get_prec(coeffs[i]));
    mpfr_set(newPolyEvalHook->coefficients[i], coeffs[i], GMP_RNDN); /* exact as precision the same */
  }

  /* Try to precompute information on whether the polynomial is monotone over the whole interval */
  newPolyEvalHook->polynomialIsMonotone = 0; 
  if (degree >= 1) {
    derivCoeffs = (mpfr_t *) safeCalloc(degree, sizeof(mpfr_t));
    derivPolyCoeffsOkay = 1;
    maxPrec = sollya_mpfi_get_prec(dom);
    if (sollya_mpfi_get_prec(t) > maxPrec) maxPrec = sollya_mpfi_get_prec(t);
    for (i=1;i<=degree;i++) {
      mpfr_init2(derivCoeffs[i-1],mpfr_get_prec(coeffs[i]) + 8 * sizeof(int) + 5);
      mpfr_mul_si(derivCoeffs[i-1],coeffs[i],i,GMP_RNDN); /* exact because precision enough */
      derivPolyCoeffsOkay = derivPolyCoeffsOkay && mpfr_number_p(derivCoeffs[i-1]);
      if (mpfr_get_prec(derivCoeffs[i-1]) > maxPrec) maxPrec = mpfr_get_prec(derivCoeffs[i-1]);
    }
    derivPoly = makePolynomial(derivCoeffs, degree-1);
    for (i=0;i<degree;i++) {
      mpfr_clear(derivCoeffs[i]);
    }
    safeFree(derivCoeffs);
    if (derivPolyCoeffsOkay) {
      sollya_mpfi_init2(X, 
			(sollya_mpfi_get_prec(dom) > sollya_mpfi_get_prec(t) ? 
			 sollya_mpfi_get_prec(dom) : sollya_mpfi_get_prec(t)));
      sollya_mpfi_init2(Y, sollya_mpfi_get_prec(X));
      sollya_mpfi_sub(X, dom, t);
      evaluateInterval(Y, derivPoly, NULL, X);
      if (!(sollya_mpfi_has_zero_inside(Y) || sollya_mpfi_has_nan(Y))) {
	/* The derivative of the polynomial does not change sign, so the 
	   polynomial is monotonous over the whole (shifted) domain.
	*/
	newPolyEvalHook->polynomialIsMonotone = 1; 
      } else {
	/* We still do not know */
	mpfr_init2(nbRoots, 8 * sizeof(int) + 5);
	if (getNrRoots(nbRoots, derivPoly, X, maxPrec, 1)) {
	  if (mpfr_zero_p(nbRoots)) {
	    /* We surely know that there is no root of the derivative polynomial in the domain */
	    newPolyEvalHook->polynomialIsMonotone = 1; 
	  }
	}
	mpfr_clear(nbRoots);
      }
      sollya_mpfi_clear(X);
      sollya_mpfi_clear(Y);
    }
    freeThing(derivPoly);
  } else {
    if (degree == 0) newPolyEvalHook->polynomialIsMonotone = 1;
  }

  /* Precompute information on whether we have an exact representation
     and/or what the maximum accuracy we can expect is 
  */
  newPolyEvalHook->exactRepresentation = 0;
  newPolyEvalHook->maxPrecKnown = 0;
  newPolyEvalHook->maxPrec = 0;
  if (sollya_mpfi_is_zero(newPolyEvalHook->delta)) {
    /* If the error interval delta is not zero, the representation is
       exact (i.e. the function we are hooking (up with ;-) is a
       polynomial.)
    */
    newPolyEvalHook->exactRepresentation = 1;
  } else {
    /* Here, the error interval delta is not zero 

       Use the hook we are just about to create to get bounds on the
       represented function. We just need some low precision.
       
    */
    if (!(sollya_mpfi_has_nan(newPolyEvalHook->delta) ||
	  sollya_mpfi_has_infinity(newPolyEvalHook->delta))) {
      sollya_mpfi_init2(Y, 12);
      if (evaluatePolyEvalHook(Y, dom, 15, (void *) newPolyEvalHook)) {
	if (!(sollya_mpfi_has_nan(Y) || 
	      sollya_mpfi_has_infinity(Y))) {
	  sollya_mpfi_abs(Y, Y);
	  mpfr_init2(ymax, 12);
	  mpfr_init2(deltamax, sollya_mpfi_get_prec(newPolyEvalHook->delta));
	  sollya_mpfi_get_right(ymax, Y);
	  sollya_mpfi_set_prec(Y, sollya_mpfi_get_prec(newPolyEvalHook->delta));
	  sollya_mpfi_abs(Y, newPolyEvalHook->delta);
	  sollya_mpfi_get_right(deltamax, Y);
	  if (mpfr_number_p(ymax) && 
	      mpfr_number_p(deltamax) && 
	      (!mpfr_zero_p(deltamax))) {
	    mpfr_div(ymax, deltamax, ymax, GMP_RNDU);
	    if (mpfr_sgn(ymax) > 0) {
	      mpfr_log2(ymax, ymax, GMP_RNDU);
	      mpfr_neg(ymax, ymax, GMP_RNDN);
	      if (mpfr_sgn(ymax) > 0) {
		mpfr_ceil(ymax, ymax);
		if (mpfr_fits_ulong_p(ymax, GMP_RNDD)) {
		  pp = mpfr_get_ui(ymax, GMP_RNDD) + 5; /* some margin accounting for the errors made here */
		  p = (mp_prec_t) pp;
		  ppp = (unsigned long) p;
		  if ((ppp == pp) && (p >= 7)) {
		    newPolyEvalHook->maxPrecKnown = 1;
		    newPolyEvalHook->maxPrec = p;
		  }
		}
	      }
	    }
	  }
	  mpfr_clear(ymax);
	  mpfr_clear(deltamax);
	  
	}
      }
      sollya_mpfi_clear(Y);
    }
  }


  return newPolyEvalHook;
}

void *copyPolyEvalHook(void *data) {
  poly_eval_hook_t *hook, *newPolyEvalHook;
  int i;

  hook = (poly_eval_hook_t *) data;

  /* Allocate memory for the hook */
  newPolyEvalHook = (poly_eval_hook_t *) safeMalloc(sizeof(poly_eval_hook_t));

  /* Set some technical stuff */
  newPolyEvalHook->reusedVarMyYInit = 0;
  newPolyEvalHook->reusedVarXInit = 0;
  newPolyEvalHook->reusedVarTempInit = 0;
  newPolyEvalHook->reusedVarMyYBInit = 0;
  newPolyEvalHook->reusedVarXAInit = 0;
  newPolyEvalHook->reusedVarXBInit = 0;
  newPolyEvalHook->reusedVarMyYRndInit = 0;
  newPolyEvalHook->reusedVarMyYRndWithDeltaInit = 0;
  newPolyEvalHook->reusedVarAInit = 0;
  newPolyEvalHook->reusedVarBInit = 0;

  /* Copy the basic stuff */
  sollya_mpfi_init2(newPolyEvalHook->domain, sollya_mpfi_get_prec(hook->domain));
  sollya_mpfi_set(newPolyEvalHook->domain, hook->domain);
  sollya_mpfi_init2(newPolyEvalHook->delta, sollya_mpfi_get_prec(hook->delta));
  sollya_mpfi_set(newPolyEvalHook->delta, hook->delta);
  sollya_mpfi_init2(newPolyEvalHook->t, sollya_mpfi_get_prec(hook->t));
  sollya_mpfi_set(newPolyEvalHook->t, hook->t);
  newPolyEvalHook->degree = hook->degree;
  newPolyEvalHook->coefficients = (mpfr_t *) safeCalloc(hook->degree + 1, sizeof(mpfr_t));
  for (i=0;i<=hook->degree;i++) {
    mpfr_init2(newPolyEvalHook->coefficients[i], mpfr_get_prec(hook->coefficients[i]));
    mpfr_set(newPolyEvalHook->coefficients[i], hook->coefficients[i], GMP_RNDN); /* exact as precision the same */
  }
  
  /* Copy the precomputed stuff */
  newPolyEvalHook->polynomialIsMonotone = hook->polynomialIsMonotone; 
  newPolyEvalHook->exactRepresentation = hook->exactRepresentation;
  newPolyEvalHook->maxPrecKnown = hook->maxPrecKnown;
  newPolyEvalHook->maxPrec = hook->maxPrec;

  /* Return the copy */
  return (void *) newPolyEvalHook;
}


int evaluatePolyEvalHook(sollya_mpfi_t y, sollya_mpfi_t x, mp_prec_t prec, void *data) {
  poly_eval_hook_t *hook;
  mp_prec_t p, pY, pX;
  int okay, i, polynomialIsMonotone;

  hook = (poly_eval_hook_t *) data;

  if (sollya_mpfi_has_nan(x)) return 0;
  if (sollya_mpfi_has_infinity(x)) return 0;
  if (!sollya_mpfi_is_inside(x, hook->domain)) return 0;
  
  pY = sollya_mpfi_get_prec(y); 
  pX = sollya_mpfi_get_prec(x); 
  p = pY + 10;
  if (prec > p) p = prec;

  if ((!(hook->exactRepresentation)) &&
      hook->maxPrecKnown &&
      (hook->maxPrec < pY)) {
    return 0;
  }

  if (hook->reusedVarMyYInit) { 
    sollya_mpfi_set_prec(hook->reusedVarMyY, p); 
  } else {
    sollya_mpfi_init2(hook->reusedVarMyY, p); 
    hook->reusedVarMyYInit = 1;
  }
  if (hook->reusedVarXInit) { 
    sollya_mpfi_set_prec(hook->reusedVarX, (p > pX ? p : pX)); 
  } else {
    sollya_mpfi_init2(hook->reusedVarX, (p > pX ? p : pX)); 
    hook->reusedVarXInit = 1;
  }
  sollya_mpfi_sub(hook->reusedVarX, x, hook->t);

  if (!sollya_mpfi_is_point_and_real(x)) {

    polynomialIsMonotone = hook->polynomialIsMonotone;
    if (!polynomialIsMonotone) {
      /* Try to determine if the polynomial is monotone over x - t 

	 Compute a bounding of its derivative and check if that
	 interval contains zero.

      */
      if (hook->reusedVarTempInit) { 
	sollya_mpfi_set_prec(hook->reusedVarTemp, p); 
      } else {
	sollya_mpfi_init2(hook->reusedVarTemp, p); 
	hook->reusedVarTempInit = 1;
      }
      sollya_mpfi_set_si(hook->reusedVarMyY, 0);
      for (i=hook->degree;i>=1;i--) {
	sollya_mpfi_mul(hook->reusedVarMyY, hook->reusedVarMyY, hook->reusedVarX);
	sollya_mpfi_set_fr(hook->reusedVarTemp, hook->coefficients[i]);
	sollya_mpfi_mul_ui(hook->reusedVarTemp, hook->reusedVarTemp, (unsigned int) i);
	sollya_mpfi_add(hook->reusedVarMyY, hook->reusedVarMyY, hook->reusedVarTemp);
      }
      if (!(sollya_mpfi_has_zero_inside(hook->reusedVarMyY) || 
	    sollya_mpfi_has_nan(hook->reusedVarMyY))) polynomialIsMonotone = 1;
    }

    if (polynomialIsMonotone) {
      /* Do a Horner evaluation in interval arithmetic over both
	 endpoint point-intervals and take the hull 
      */
      if (hook->reusedVarMyYBInit) { 
	sollya_mpfi_set_prec(hook->reusedVarMyYB, p); 
      } else {
	sollya_mpfi_init2(hook->reusedVarMyYB, p); 
	hook->reusedVarMyYBInit = 1;
      }
      if (hook->reusedVarXAInit) { 
	sollya_mpfi_set_prec(hook->reusedVarXA, sollya_mpfi_get_prec(hook->reusedVarX)); 
      } else {
	sollya_mpfi_init2(hook->reusedVarXA, sollya_mpfi_get_prec(hook->reusedVarX)); 
	hook->reusedVarXAInit = 1;
      }
      if (hook->reusedVarXBInit) { 
	sollya_mpfi_set_prec(hook->reusedVarXB, sollya_mpfi_get_prec(hook->reusedVarX)); 
      } else {
	sollya_mpfi_init2(hook->reusedVarXB, sollya_mpfi_get_prec(hook->reusedVarX)); 
	hook->reusedVarXBInit = 1;
      }
      if (hook->reusedVarAInit) { 
	mpfr_set_prec(hook->reusedVarA, sollya_mpfi_get_prec(hook->reusedVarX)); 
      } else {
	mpfr_init2(hook->reusedVarA, sollya_mpfi_get_prec(hook->reusedVarX)); 
	hook->reusedVarAInit = 1;
      }
      if (hook->reusedVarBInit) { 
	mpfr_set_prec(hook->reusedVarB, sollya_mpfi_get_prec(hook->reusedVarX)); 
      } else {
	mpfr_init2(hook->reusedVarB, sollya_mpfi_get_prec(hook->reusedVarX)); 
	hook->reusedVarBInit = 1;
      }
      sollya_mpfi_get_left(hook->reusedVarA, hook->reusedVarX);
      sollya_mpfi_get_right(hook->reusedVarB, hook->reusedVarX);
      sollya_mpfi_set_fr(hook->reusedVarXA, hook->reusedVarA);
      sollya_mpfi_set_fr(hook->reusedVarXB, hook->reusedVarB);
      
      sollya_mpfi_set_si(hook->reusedVarMyY, 0);
      sollya_mpfi_set_si(hook->reusedVarMyYB, 0);
      for (i=hook->degree;i>=0;i--) {
	sollya_mpfi_mul(hook->reusedVarMyY, hook->reusedVarMyY, hook->reusedVarXA);
	sollya_mpfi_mul(hook->reusedVarMyYB, hook->reusedVarMyYB, hook->reusedVarXB);
	sollya_mpfi_add_fr(hook->reusedVarMyY, hook->reusedVarMyY, hook->coefficients[i]);
	sollya_mpfi_add_fr(hook->reusedVarMyYB, hook->reusedVarMyYB, hook->coefficients[i]);
      }      
      sollya_mpfi_union(hook->reusedVarMyY, hook->reusedVarMyY, hook->reusedVarMyYB);

      okay = 0;
      if (hook->reusedVarMyYRndInit) { 
	sollya_mpfi_set_prec(hook->reusedVarMyYRnd, pY + 5); 
      } else {
	sollya_mpfi_init2(hook->reusedVarMyYRnd, pY + 5); 
	hook->reusedVarMyYRndInit = 1;
      }
      if (hook->reusedVarMyYRndWithDeltaInit) { 
	sollya_mpfi_set_prec(hook->reusedVarMyYRndWithDelta, pY + 5); 
      } else {
	sollya_mpfi_init2(hook->reusedVarMyYRndWithDelta, pY + 5); 
	hook->reusedVarMyYRndWithDeltaInit = 1;
      }
      sollya_mpfi_set(hook->reusedVarMyYRnd, hook->reusedVarMyY);
      sollya_mpfi_blow_1ulp(hook->reusedVarMyYRnd);
      sollya_mpfi_add(hook->reusedVarMyYRndWithDelta, hook->reusedVarMyY, hook->delta);

      if (sollya_mpfi_is_inside(hook->reusedVarMyYRndWithDelta, hook->reusedVarMyYRnd) && 
	  (!(sollya_mpfi_has_nan(hook->reusedVarMyYRndWithDelta) || 
	     (sollya_mpfi_has_infinity(hook->reusedVarMyYRndWithDelta) && 
	      (!sollya_mpfi_is_infinity(hook->reusedVarMyYRndWithDelta)))))) okay = 1;
  
      if (okay) sollya_mpfi_set(y, hook->reusedVarMyYRndWithDelta);

      return okay;
    }
    /* Fall through if we could not determine monotononicity */
  }

  /* Regular Horner evaluation with interval arithmetic */
  sollya_mpfi_set_fr(hook->reusedVarMyY, hook->coefficients[hook->degree]);
  for (i=hook->degree-1;i>=0;i--) {
    sollya_mpfi_mul(hook->reusedVarMyY, hook->reusedVarMyY, hook->reusedVarX);
    sollya_mpfi_add_fr(hook->reusedVarMyY, hook->reusedVarMyY, hook->coefficients[i]);
  }
  
  okay = 0;
  if (hook->reusedVarMyYRndInit) { 
    sollya_mpfi_set_prec(hook->reusedVarMyYRnd, pY + 5); 
  } else {
    sollya_mpfi_init2(hook->reusedVarMyYRnd, pY + 5); 
    hook->reusedVarMyYRndInit = 1;
  }
  if (hook->reusedVarMyYRndWithDeltaInit) { 
    sollya_mpfi_set_prec(hook->reusedVarMyYRndWithDelta, pY + 5); 
  } else {
    sollya_mpfi_init2(hook->reusedVarMyYRndWithDelta, pY + 5); 
    hook->reusedVarMyYRndWithDeltaInit = 1;
  }

  sollya_mpfi_set(hook->reusedVarMyYRnd, hook->reusedVarMyY);

  sollya_mpfi_blow_1ulp(hook->reusedVarMyYRnd);
  sollya_mpfi_add(hook->reusedVarMyYRndWithDelta, hook->reusedVarMyY, hook->delta);

  if (sollya_mpfi_is_inside(hook->reusedVarMyYRndWithDelta, hook->reusedVarMyYRnd) && 
      (!(sollya_mpfi_has_nan(hook->reusedVarMyYRndWithDelta) || 
	 (sollya_mpfi_has_infinity(hook->reusedVarMyYRndWithDelta) && 
	  (!sollya_mpfi_is_infinity(hook->reusedVarMyYRndWithDelta)))))) okay = 1;
  
  if (okay) sollya_mpfi_set(y, hook->reusedVarMyYRndWithDelta);

  return okay;
}

void freePolyEvalHook(void *data) {
  poly_eval_hook_t *hook;
  int i;

  hook = (poly_eval_hook_t *) data;
  sollya_mpfi_clear(hook->domain);
  sollya_mpfi_clear(hook->delta);
  sollya_mpfi_clear(hook->t);
  for (i=0;i<=hook->degree;i++) {
    mpfr_clear(hook->coefficients[i]);
  }
  safeFree(hook->coefficients);
  if (hook->reusedVarMyYInit) sollya_mpfi_clear(hook->reusedVarMyY);
  if (hook->reusedVarXInit) sollya_mpfi_clear(hook->reusedVarX);
  if (hook->reusedVarTempInit) sollya_mpfi_clear(hook->reusedVarTemp);
  if (hook->reusedVarMyYBInit) sollya_mpfi_clear(hook->reusedVarMyYB);
  if (hook->reusedVarXAInit) sollya_mpfi_clear(hook->reusedVarXA);
  if (hook->reusedVarXBInit) sollya_mpfi_clear(hook->reusedVarXB);
  if (hook->reusedVarMyYRndInit) sollya_mpfi_clear(hook->reusedVarMyYRnd);
  if (hook->reusedVarMyYRndWithDeltaInit) sollya_mpfi_clear(hook->reusedVarMyYRndWithDelta);
  if (hook->reusedVarAInit) mpfr_clear(hook->reusedVarA);
  if (hook->reusedVarBInit) mpfr_clear(hook->reusedVarB);
  safeFree(hook);
}

int comparePolyEvalHook(void *data1, void *data2) {
  poly_eval_hook_t *hook1, *hook2;
  int i;

  hook1 = (poly_eval_hook_t *) data1;
  hook2 = (poly_eval_hook_t *) data2;

  if (!sollya_mpfi_equal_p(hook1->domain, hook2->domain)) return 0;
  if (!sollya_mpfi_equal_p(hook1->delta, hook2->delta)) return 0;
  if (!sollya_mpfi_equal_p(hook1->t, hook2->t)) return 0;
  if (hook1->degree != hook2->degree) return 0;
  for (i=0;i<=hook1->degree;i++) {
    if (!mpfr_equal_p(hook1->coefficients[i], hook2->coefficients[i])) return 0;
  }

  return 1;
}


int addPolyEvaluationHook(eval_hook_t **hookPtr, node *func, sollya_mpfi_t dom, sollya_mpfi_t delta, sollya_mpfi_t t, mp_prec_t prec) {
  int degree, okay, i, coeffsOkay;
  node **coeffs;
  mpfr_t *evaluatedCoeffs;
  sollya_mpfi_t globalDelta, c, shiftedDom;
  mp_prec_t p, evalPrec, pp;
  mp_exp_t maxExp;
  
  if (!isPolynomial(func)) return 0;

  getCoefficients(&degree, &coeffs, func);
  evaluatedCoeffs = (mpfr_t *) safeCalloc(degree+1,sizeof(mpfr_t));
  evalPrec = prec;
  if (!(sollya_mpfi_is_zero(delta) ||
	sollya_mpfi_has_nan(delta) ||
	sollya_mpfi_has_infinity(delta))) {
    maxExp = sollya_mpfi_max_exp(delta);
    if (maxExp < 0) {
      pp = -maxExp;
      if ((pp > evalPrec) && (pp < 4 * evalPrec)) {
	evalPrec = pp;
	if (evalPrec < 12) evalPrec = 12;
      }
    }
  }
  sollya_mpfi_init2(globalDelta, (sollya_mpfi_get_prec(delta) > evalPrec ? sollya_mpfi_get_prec(delta) : evalPrec));
  p = prec;
  if (sollya_mpfi_get_prec(dom) > p) p = sollya_mpfi_get_prec(dom);
  if (sollya_mpfi_get_prec(t) > p) p = sollya_mpfi_get_prec(t);
  sollya_mpfi_init2(shiftedDom, p);
  sollya_mpfi_sub(shiftedDom, dom, t);
  sollya_mpfi_set_si(globalDelta, 0);
  sollya_mpfi_init2(c, evalPrec + 5);
  coeffsOkay = 1;
  for (i=degree;i>=0;i--) {
    if ((coeffs[i] != NULL) && (accessThruMemRef(coeffs[i])->nodeType == CONSTANT)) {
      mpfr_init2(evaluatedCoeffs[i],mpfr_get_prec(*(accessThruMemRef(coeffs[i])->value)));
      mpfr_set(evaluatedCoeffs[i],*(accessThruMemRef(coeffs[i])->value),GMP_RNDN); /* exact */
      sollya_mpfi_set_si(c, 0);
      freeThing(coeffs[i]);
    } else {
      mpfr_init2(evaluatedCoeffs[i],evalPrec);
      if (coeffs[i] == NULL) {
	sollya_mpfi_set_si(c, 0);
      } else {
	evaluateConstantExpressionToSharpInterval(c, coeffs[i]);
	freeThing(coeffs[i]);
      }
      sollya_mpfi_mid(evaluatedCoeffs[i], c);
      coeffsOkay = coeffsOkay && mpfr_number_p(evaluatedCoeffs[i]);
      sollya_mpfi_sub_fr(c, c, evaluatedCoeffs[i]);
    }
    sollya_mpfi_mul(globalDelta, globalDelta, shiftedDom);
    sollya_mpfi_add(globalDelta, globalDelta, c);
  }
  safeFree(coeffs);
  sollya_mpfi_add(globalDelta, globalDelta, delta);
  sollya_mpfi_clear(c);
  sollya_mpfi_clear(shiftedDom);
  
  if (coeffsOkay && (!(sollya_mpfi_has_nan(globalDelta) || sollya_mpfi_has_infinity(globalDelta)))) {
    okay = addEvaluationHook(hookPtr, (void *) createPolyEvalHook(degree, evaluatedCoeffs, dom, globalDelta, t),
			     evaluatePolyEvalHook, freePolyEvalHook, comparePolyEvalHook, copyPolyEvalHook); 
  } else {
    okay = 0;
  }

  sollya_mpfi_clear(globalDelta);
  for (i=0;i<=degree;i++) {
    mpfr_clear(evaluatedCoeffs[i]);
  }
  safeFree(evaluatedCoeffs);

  return okay;
}

int addNodeEvaluationHook(eval_hook_t **hookPtr, node *func, sollya_mpfi_t dom, sollya_mpfi_t delta, sollya_mpfi_t t, mp_prec_t prec) {
  UNUSED_PARAM(prec);
  return addEvaluationHook(hookPtr, (void *) createNodeEvalHook(func, dom, delta, t), 
			   evaluateNodeEvalHook, freeNodeEvalHook, compareNodeEvalHook, copyNodeEvalHook); 
}

int chooseAndAddEvaluationHook(eval_hook_t **hookPtr, node *func, sollya_mpfi_t dom, sollya_mpfi_t delta, sollya_mpfi_t t, mp_prec_t prec) {
  if (isPolynomial(func)) return addPolyEvaluationHook(hookPtr, func, dom, delta, t, prec);
  return addNodeEvaluationHook(hookPtr, func, dom, delta, t, prec);
}

int copyFunctionAndChooseAndAddEvaluationHook(node **copyPtr, node *orig, node *func, sollya_mpfi_t dom, sollya_mpfi_t delta, sollya_mpfi_t t, mp_prec_t prec) {
  node *copy;
  int okay;

  copy = addMemRef(copyThing(orig));
 
  if (copy->nodeType != MEMREF) {
    freeThing(copy);
    return 0;
  }

  okay = chooseAndAddEvaluationHook(&(copy->evaluationHook), func, dom, delta, t, prec);
    
  if (okay) {
    *copyPtr = copy;
  } else {
    freeThing(copy);
  }

  return okay;
}



