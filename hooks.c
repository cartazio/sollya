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
#include "hooks.h"

int addEvaluationHook(eval_hook_t **hookPtr, 
		      void *data, 
		      int (*evaluateFunc)(sollya_mpfi_t, sollya_mpfi_t, mp_prec_t, void *), 
		      void (*freeFunc)(void *),
		      int (*compareFunc)(void *, void *)) {
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
  newHook->nextHook = *hookPtr;

  /* Assign the new hook */
  *hookPtr = newHook;
  
  /* Signal success */
  return 1;
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

  for (curr=hook;curr!=NULL;curr=curr->nextHook) {
    if (curr->evaluateHook(y,x,prec,curr->data)) return 1;
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
  sollya_mpfi_init2(myYRnd, pY);
  sollya_mpfi_init2(myYRndWithDelta, pY);

  sollya_mpfi_set(myYRnd, myY);
  sollya_mpfi_add(myYRndWithDelta, myY, hook->delta);
  
  if (sollya_mpfi_is_inside(myYRnd, myYRndWithDelta)) okay = 1;
  
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


poly_eval_hook_t *createPolyEvalHook(int degree, mpfr_t *coeffs, sollya_mpfi_t dom, sollya_mpfi_t delta, sollya_mpfi_t t) {
  poly_eval_hook_t *newPolyEvalHook;
  int i;

  newPolyEvalHook = (poly_eval_hook_t *) safeMalloc(sizeof(poly_eval_hook_t));
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

  return newPolyEvalHook;
}

int evaluatePolyEvalHook(sollya_mpfi_t y, sollya_mpfi_t x, mp_prec_t prec, void *data) {

  return 0;
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





