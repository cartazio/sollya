/*

Copyright 2011 by 

Laboratoire d'Informatique de Paris 6, equipe PEQUAN,
UPMC Universite Paris 06 - CNRS - UMR 7606 - LIP6, Paris, France.

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
herefore means  that it is reserved for developers  and  experienced
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

#include "expression.h"
#include "execute.h"
#include "chain.h"
#include "assignment.h"
#include <string.h>
#include <stdlib.h>

/* Add association identifier -> thing to list 
   
   If identifier exists already:
      - And is associated to a thing that is equal to thing, return success (1)
      - And is associated to a different thing, return failure (0)

   Otherwise: create the association identifier -> thing and add it to the 
              list, return success (1)

*/
int associateThing(chain **associations, char *identifier, node *thing) {
  chain *curr;
  entry *newEntry;

  for (curr = *associations; curr != NULL; curr = curr->next) {
    if (!strcmp(identifier, ((char *) (((entry *) (curr->value))->name)))) {
      /* Identifier exists already in list */
      if (isEqualThing(thing, ((node *) (((entry *) (curr->value))->value)))) {
	/* But it is associated to the same thing */
	return 1;
      } else {
	/* But it is associated to another thing */
	return 0;
      }
    }
  }

  newEntry = (entry *) safeMalloc(sizeof(entry));
  newEntry->name = (char *) safeCalloc(strlen(identifier) + 1,sizeof(char));
  strcpy(newEntry->name, identifier);
  newEntry->value = copyThing(thing);
  *associations = addElement(*associations, newEntry);

  return 1;
}

int tryCombineAssociations(chain **associations, chain *assoc1, chain *assoc2) {
  int okay;
  chain *myAssociations = NULL;
  chain *curr;

  if ((assoc1 == NULL) && (assoc2 == NULL)) {
    *associations = NULL;
    return 1;
  }

  if (assoc1 == NULL) {
    *associations = copyChain(assoc2,copyEntryOnVoid);
    return 1;
  }

  if (assoc2 == NULL) {
    *associations = copyChain(assoc1,copyEntryOnVoid);
    return 1;
  }

  okay = 1;
  myAssociations = copyChain(assoc1,copyEntryOnVoid);

  for (curr=assoc2; curr != NULL; curr = curr->next) {
    if (!associateThing(&myAssociations, 
			((char *) ((entry *) (curr->value))->name), 
			((node *) ((entry *) (curr->value))->value))) {
      okay = 0;
      break;
    }
  }

  if (okay) {
    *associations = myAssociations;
  } else {
    if (myAssociations != NULL) {
      freeChain(myAssociations, freeEntryOnVoid);
    }
  }

  return okay;
}

int tryMatchExtendedPureTree(chain **associations, node *thingToMatch, node *possibleMatcher) {
  chain *leftAssoc, *rightAssoc;
  int okay;
  
  /* Special case: possibleMatcher is a free variable to bind
     Check if it is possible equal to the mathematical free 
     variable. If not, create an association 
  */
  if (possibleMatcher->nodeType == TABLEACCESS) {
    if ((variablename != NULL) &&
	(!strcmp(variablename, possibleMatcher->string))) {
      /* Here, the free variable to bind is actually equal
	 to the free mathematical variable. 
	 
	 We have a match if the thing to match is 
	 also equal to the free mathematical variable.

	 Otherwise, we have no match.

	 In none of the cases, we have to establish an 
	 association nor call ourselves for recursion.
      */
      return (thingToMatch->nodeType == VARIABLE);
    }
    /* Here, the free variable is not the mathematical one.
       We have to establish an association and return success
       (unless the association fails).
    */
    return associateThing(associations, possibleMatcher->string, thingToMatch);
  }

  /* Here, possibleMatcher is itself a pure tree, i.e. a function without
     free variables to bind.

     To start with, we cannot have a match if the head symbols are not
     the same.
  */
  if (possibleMatcher->nodeType != thingToMatch->nodeType) return 0;

  /* Here, the head symbols are always the same */
  switch (possibleMatcher->nodeType) {
  case VARIABLE:
    /* The free mathematical variable matches the free mathematical variable */
    return 1;
    break;
  case CONSTANT:
    /* Constants match if they are the same. We have to be careful with NaNs */
    if (mpfr_nan_p(*(possibleMatcher->value)) && mpfr_nan_p(*(thingToMatch->value))) return 1;
    if (mpfr_equal_p(*(possibleMatcher->value),*(thingToMatch->value))) return 1; 
    return 0;
    break;
  case ADD:
  case SUB:
  case MUL:
  case DIV:
  case POW:
    /* Binary functions match if both sub-expressions match and if we can combine ("unify") 
       both association lists.
    */
    leftAssoc = NULL;
    rightAssoc = NULL;
    if (!tryMatchExtendedPureTree(&leftAssoc, thingToMatch->child1, possibleMatcher->child1)) return 0;
    if (!tryMatchExtendedPureTree(&rightAssoc, thingToMatch->child2, possibleMatcher->child2)) return 0;
    okay = tryCombineAssociations(associations, leftAssoc, rightAssoc);
    if (leftAssoc != NULL) freeChain(leftAssoc, freeEntryOnVoid);
    if (rightAssoc != NULL) freeChain(rightAssoc, freeEntryOnVoid);
    return okay;
    break;
  case SQRT:
  case EXP:
  case LOG:
  case LOG_2:
  case LOG_10:
  case SIN:
  case COS:
  case TAN:
  case ASIN:
  case ACOS:
  case ATAN:
  case SINH:
  case COSH:
  case TANH:
  case ASINH:
  case ACOSH:
  case ATANH:
  case NEG:
  case ABS:
  case DOUBLE:
  case SINGLE:
  case DOUBLEDOUBLE:
  case TRIPLEDOUBLE:
  case ERF: 
  case ERFC:
  case LOG_1P:
  case EXP_M1:
  case DOUBLEEXTENDED:
  case CEIL:
  case FLOOR:
  case NEARESTINT:
    /* Unary functions: we have a match if we have a match on recursion
       Possible associations are established on recursion. 
    */
    return tryMatchExtendedPureTree(associations, thingToMatch->child1, possibleMatcher->child1);
    break;
  case PI_CONST:
    /* Pi matches pi */
    return 1;
    break;
  case LIBRARYFUNCTION:
    /* Library functions to not match if they are not bound to the same code
       or if they are at different differentiation levels
    */
    if ((possibleMatcher->libFun != thingToMatch->libFun) ||
	(possibleMatcher->libFunDeriv != thingToMatch->libFunDeriv)) return 0;
    /* Otherwise they match when recursion matches */
    return tryMatchExtendedPureTree(associations, thingToMatch->child1, possibleMatcher->child1);
    break;
  case LIBRARYCONSTANT:
    /* Library constants match iff they are bound to the same code */
    return (possibleMatcher->libFun == thingToMatch->libFun);
    break;
  case PROCEDUREFUNCTION:
    /* Procedure functions to not match if they are not based on the same
       procedure or if they are at different differentiation levels
    */
    if ((!isEqualThing(possibleMatcher->child2, thingToMatch->child2)) ||
	(possibleMatcher->libFunDeriv != thingToMatch->libFunDeriv)) return 0;
    /* Otherwise they match when recursion matches */
    return tryMatchExtendedPureTree(associations, thingToMatch->child1, possibleMatcher->child1);    
    break;
  default:
    sollyaFprintf(stderr,"Error: tryMatchExtendedPureTree: unknown identifier (%d) in the tree\n",possibleMatcher->nodeType);
    exit(1);
  }

  /* We should never get here */
  return 0;
}


int tryMatchInner(chain **associations, node *thingToMatch, node *possibleMatcher) {

  /* Default case: a match for everything, no association to perform */
  if (possibleMatcher->nodeType == DEFAULT) return 1;

  /* Base symbols: a match is obtained if the thing to match is equal
     to the matcher */
  if (isCorrectlyTypedBaseSymbol(possibleMatcher)) {
    if (isEqualThing(thingToMatch, possibleMatcher)) {
      return 1;
    }
    return 0;
  }

  /* Extended pure trees (functions with possible variables to bind to): 
     run specialized code if the thing to match is a pure tree (a function 
     without variables)
  */
  if (isExtendedPureTree(possibleMatcher)) {
    if (!isPureTree(thingToMatch)) return 0;
    return tryMatchExtendedPureTree(associations, thingToMatch, possibleMatcher);
  }

  /* TODO */


  return 0;
}



int tryMatch(chain **associations, node *thingToMatch, node *possibleMatcher) {
  int okay;
  chain *myAssociations = NULL;

  sollyaPrintf("Trying to match ");
  rawPrintThing(thingToMatch);
  sollyaPrintf(" with ");
  rawPrintThing(possibleMatcher);
  sollyaPrintf("\n");

  okay = tryMatchInner(&myAssociations, thingToMatch, possibleMatcher);

  sollyaPrintf("okay = %d\n",okay);

  if (okay) {
    *associations = myAssociations;
  } else {
    if (myAssociations != NULL) {
      freeChain(myAssociations, freeEntryOnVoid);
    }
  }

  return okay;
}
