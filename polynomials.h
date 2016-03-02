/*

  Copyright 2014-2015 by

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

#ifndef POLYNOMIALS_H
#define POLYNOMIALS_H

#include <stdint.h>
#include <gmp.h>
#include <mpfr.h>
#include <stdio.h>
#include "mpfi-compat.h"


/* An abstract type for polynomials */
typedef struct __polynomial_struct_t * polynomial_t;


/* Operations on polynomials */

/* Constructors */
polynomial_t polynomialFromMpfrConstant(mpfr_t);
polynomial_t polynomialFromMpzConstant(mpz_t);
polynomial_t polynomialFromMpqConstant(mpq_t);
polynomial_t polynomialFromIntConstant(int);
polynomial_t polynomialFromIdentity();
polynomial_t polynomialFromMpfrCoefficients(mpfr_t *, unsigned int);
int polynomialFromConstantExpressionCoefficients(polynomial_t *, struct nodeStruct **, unsigned int);
int polynomialFromExpression(polynomial_t *, struct nodeStruct *);
int polynomialFromExpressionOnlyRealCoeffs(polynomial_t *, struct nodeStruct *);

/* Copy-Constructor */
polynomial_t polynomialFromCopy(polynomial_t);

/* Destructor */
void polynomialFree(polynomial_t);

/* Comparisons */
int polynomialEqual(polynomial_t, polynomial_t, int);
int polynomialIsIdentity(polynomial_t, int);
int polynomialIsConstant(polynomial_t, int);
int polynomialStructurallyEqual(polynomial_t, polynomial_t, int);

/* Arithmetical operations */
polynomial_t polynomialAdd(polynomial_t, polynomial_t);
polynomial_t polynomialSub(polynomial_t, polynomial_t);
polynomial_t polynomialMul(polynomial_t, polynomial_t);
polynomial_t polynomialNeg(polynomial_t);
polynomial_t polynomialCompose(polynomial_t, polynomial_t);
void polynomialDiv(polynomial_t *, polynomial_t *, polynomial_t, polynomial_t);
int polynomialPow(polynomial_t *, polynomial_t, polynomial_t);
polynomial_t polynomialPowUnsignedInt(polynomial_t, unsigned int);
polynomial_t polynomialDeriv(polynomial_t);

/* Rewrite operations */
polynomial_t polynomialHornerize(polynomial_t);
polynomial_t polynomialCanonicalize(polynomial_t);

/* Tests for rewrite operations */
int polynomialIsHornerized(polynomial_t);
int polynomialIsCanonicalized(polynomial_t);

/* Accessors */
void polynomialGetDegree(mpz_t, polynomial_t);
int polynomialGetDegreeAsInt(polynomial_t);
struct nodeStruct *polynomialGetIthCoefficient(polynomial_t, mpz_t);
struct nodeStruct *polynomialGetIthCoefficientIntIndex(polynomial_t, int);
int polynomialGetCoefficients(struct nodeStruct ***, unsigned int *, polynomial_t); 
struct nodeStruct *polynomialGetExpression(polynomial_t);
struct nodeStruct *polynomialGetExpressionExplicit(polynomial_t);

/* Displaying and conversion to strings */
void polynomialFPrintf(FILE *, polynomial_t);
char *polynomialToString(polynomial_t);

/* Evaluation */
void polynomialEvalMpfr(mpfr_t, polynomial_t, mpfr_t);
void polynomialEvalMpfi(sollya_mpfi_t, polynomial_t, sollya_mpfi_t);

/* Rounding of coefficients */
int polynomialCoefficientsAreDyadic(polynomial_t, int);
int polynomialCoefficientsAreRational(polynomial_t, int);
polynomial_t polynomialRoundDyadic(polynomial_t, mp_prec_t);
polynomial_t polynomialRoundRational(polynomial_t, mp_prec_t);
polynomial_t polynomialRound(polynomial_t, mp_prec_t);

/* A function to prevent memory reference loops */
int polynomialReferencesExpression(polynomial_t, struct nodeStruct *);

/* A hashing function */
uint64_t polynomialHash(polynomial_t);


#endif /* ifdef POLYNOMIALS_H*/
