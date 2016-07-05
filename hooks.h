/*

  Copyright 2014 by

  Laboratoire d'Informatique de Paris 6, equipe PEQUAN,
  UPMC Universite Paris 06 - CNRS - UMR 7606 - LIP6, Paris, France

  Contributor Ch. Lauter

  christoph.lauter@lip6.fr

  This software is a computer program whose purpose is to provide an
  environment for safe floating-point code development. It is
  particularly targeted to the automated implementation of
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

#ifndef HOOKS_H
#define HOOKS_H

#include <gmp.h>
#include <mpfr.h>
#include "mpfi-compat.h"

/* We need to know the nodeStruct structure */

struct nodeStruct;

/* General framework for evaluation hooks */

typedef struct __eval_hook_t_struct eval_hook_t;
struct __eval_hook_t_struct {
  void *data;
  int gettingUsed;
  int reuseMPFIInit;
  sollya_mpfi_t reuseMPFI;
  int (*evaluateHook)(sollya_mpfi_t, sollya_mpfi_t, mp_prec_t, int, void *);
  void (*freeHook)(void *);
  int (*compareHook)(void *, void *);
  void *(*copyHook)(void *);
  int (*composeHook)(eval_hook_t **, void *, struct nodeStruct *);
  eval_hook_t *nextHook;
};

int addEvaluationHook(eval_hook_t **, 
		      void *, 
		      int (*)(sollya_mpfi_t, sollya_mpfi_t, mp_prec_t, int, void *), 
		      void (*)(void *),
		      int (*)(void *, void *),
		      void *(*)(void *),
		      int (*)(eval_hook_t **, void *, struct nodeStruct *));

int addEvaluationHookFromCopy(eval_hook_t **, eval_hook_t *);

int addEvaluationHookFromComposition(eval_hook_t **, eval_hook_t *, struct nodeStruct *);

void freeEvaluationHook(eval_hook_t **);

int evaluateWithEvaluationHook(sollya_mpfi_t, sollya_mpfi_t, mp_prec_t, int, eval_hook_t *);



/* General one-replaces-one node-represented function evaluation hooks */

typedef struct __node_eval_hook_t_struct node_eval_hook_t;
struct __node_eval_hook_t_struct {
  sollya_mpfi_t domain;
  sollya_mpfi_t delta;
  sollya_mpfi_t t;
  struct nodeStruct *func;
};

node_eval_hook_t *createNodeEvalHook(struct nodeStruct *, sollya_mpfi_t, sollya_mpfi_t, sollya_mpfi_t);
int evaluateNodeEvalHook(sollya_mpfi_t, sollya_mpfi_t, mp_prec_t, int, void *);
void freeNodeEvalHook(void *);
int compareNodeEvalHook(void *, void *);
void *copyNodeEvalHook(void *);
int composeNodeEvalHook(eval_hook_t **, void *, struct nodeStruct *);

/* Polynomial replacement function evaluation hooks */

typedef struct __poly_eval_hook_t_struct poly_eval_hook_t;
struct __poly_eval_hook_t_struct {
  sollya_mpfi_t domain;
  sollya_mpfi_t delta;
  sollya_mpfi_t t;
  int degree;
  int polynomialIsMonotone;
  int polynomialHasZero;
  int maxPrecKnown;
  mp_prec_t maxPrec;
  int exactRepresentation;
  mpfr_t *coefficients;
  int reusedVarMyYInit;
  sollya_mpfi_t reusedVarMyY;
  int reusedVarXInit;
  sollya_mpfi_t reusedVarX;
  int reusedVarTempInit;
  sollya_mpfi_t reusedVarTemp;
  int reusedVarMyYBInit;
  sollya_mpfi_t reusedVarMyYB;
  int reusedVarXAInit;
  sollya_mpfi_t reusedVarXA;
  int reusedVarXBInit;
  sollya_mpfi_t reusedVarXB;
  int reusedVarMyYRndInit;
  sollya_mpfi_t reusedVarMyYRnd;
  int reusedVarMyYRndWithDeltaInit;
  sollya_mpfi_t reusedVarMyYRndWithDelta;
  int reusedVarAInit;
  mpfr_t reusedVarA;
  int reusedVarBInit;
  mpfr_t reusedVarB;
};

poly_eval_hook_t *createPolyEvalHook(int, mpfr_t *, sollya_mpfi_t, sollya_mpfi_t, sollya_mpfi_t);
int evaluatePolyEvalHook(sollya_mpfi_t, sollya_mpfi_t, mp_prec_t, int, void *);
void freePolyEvalHook(void *);
int comparePolyEvalHook(void *, void *);
void *copyPolyEvalHook(void *);
int composePolyEvalHook(eval_hook_t **, void *, struct nodeStruct *);

/* Composition evaluation hooks */

typedef struct __composition_eval_hook_t_struct composition_eval_hook_t;
struct __composition_eval_hook_t_struct {
  eval_hook_t *f;
  struct nodeStruct *g;
  sollya_mpfi_t reusedVarT;
  sollya_mpfi_t reusedVarTA;
  sollya_mpfi_t reusedVarTB;
  mpfr_t reusedVarTemp;
  int reusedVarTInit;
  int reusedVarTAInit;
  int reusedVarTBInit;
  int reusedVarTempInit;
};

composition_eval_hook_t *createCompositionEvalHook(eval_hook_t *, struct nodeStruct *);
int evaluateCompositionEvalHook(sollya_mpfi_t, sollya_mpfi_t, mp_prec_t, int, void *);
void freeCompositionEvalHook(void *);
int compareCompositionEvalHook(void *, void *);
void *copyCompositionEvalHook(void *);
int composeCompositionEvalHook(eval_hook_t **, void *, struct nodeStruct *);

/* Helper functions that install either a general or a polynomial replacement hook */

int chooseAndAddEvaluationHook(eval_hook_t **, struct nodeStruct *, sollya_mpfi_t, sollya_mpfi_t, sollya_mpfi_t, mp_prec_t);
int copyFunctionAndChooseAndAddEvaluationHook(struct nodeStruct **, struct nodeStruct *, struct nodeStruct *, sollya_mpfi_t, sollya_mpfi_t, sollya_mpfi_t, mp_prec_t);




#endif /* ifdef HOOKS_H*/
