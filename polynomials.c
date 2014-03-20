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
#include "expression.h"
#include "polynomials.h"

typedef enum __constant_type_enum_t constant_type_t;
enum __constant_type_enum_t {
  EXPRESSION,
  MPFR,
  SCALEDMPQ
};

typedef struct __scaled_mpq_struct_t scaled_mpq_t;
struct __scaled_mpq_struct_t {
  mp_exp_t expo;
  mpq_t significand;
};

typedef struct __constant_struct_t * constant_t;
struct __constant_struct_t {
  unsigned int refCount;
  constant_type_t type;
  union { 
    node *expr;
    mpfr_t mpfr;
    scaled_mpq_t scaledMpq;
  } value;
};

struct __polynomial_struct_t {
  unsigned int refCount;
  mpz_t deg;
  unsigned int monomialCount;
  constant_t *coeffs;
  mpz_t *monomialDegrees;
};

constant_t constantFromExpression(node *c) {
  return NULL;
}

constant_t constantFromMpfr(mpfr_t c) {
  return NULL;
}

constant_t constantFromMpz(mpz_t c) {
  return NULL;
}

constant_t constantFromMpq(mpq_t c) {
  return NULL;
}

constant_t constantFromCopy(constant_t c) {
  return NULL;
}

void freeConstant(constant_t c) {
  
}

constant_t constantAdd(constant_t a, constant_t b) {
  return NULL;
}

constant_t constantSub(constant_t a, constant_t b) {
  return NULL;
}

constant_t constantMul(constant_t a, constant_t b) {
  return NULL;
}

constant_t constantDiv(constant_t a, constant_t b) {
  return NULL;
}

constant_t constantPow(constant_t a, constant_t b) {
  return NULL;
}

constant_t constantNeg(constant_t a) {
  return NULL;
}

int constantIsZero(constant_t a, int defVal) {
  return defVal;
}

int constantIsOne(constant_t a, int defVal) {
  return defVal;
}

node *constantToExpression(constant_t a) {
  return NULL;
}



polynomial_t polynomialFromExpression(node *p) {
  return NULL;
}

polynomial_t polynomialFromMpfrConstant(mpfr_t c) {
  return NULL;
}

polynomial_t polynomialFromMpzConstant(mpz_t c) {
  return NULL;
}

polynomial_t polynomialFromMpqConstant(mpq_t c) {
  return NULL;
}

polynomial_t polynomialFromIdentity() {
  return NULL;
}

polynomial_t polynomialFromMpfrCoefficients(mpfr_t *coeffs, int deg) {
  return NULL;
}

polynomial_t polynomialFromConstantExpressionCoefficients(node **coeffs, int deg) {
  return NULL;
}


polynomial_t polynomialCopy(polynomial_t p) {
  return NULL;
}

void freePolynomial(polynomial_t p) {
}

polynomial_t polynomialAdd(polynomial_t p, polynomial_t q) {
  return NULL;
}

polynomial_t polynomialSub(polynomial_t p, polynomial_t q) {
  return NULL;
}

polynomial_t polynomialMul(polynomial_t p, polynomial_t q) {
  return NULL;
}

polynomial_t polynomialNeg(polynomial_t p) {
  return NULL;
}

polynomial_t polynomialCompose(polynomial_t p, polynomial_t q) {
  return NULL;
}

void polynomialDiv(polynomial_t *quot, polynomial_t *rest, polynomial_t p, polynomial_t q) {
}

polynomial_t polynomialMpzPow(polynomial_t p, mpz_t c) {
  return NULL;
}

polynomial_t polynomialPow(polynomial_t p, polynomial_t q) {
  return NULL;
}

void polynomialGetDegree(mpz_t deg, polynomial_t p) {
}

int polynomialGetDegreeAsInt(polynomial_t p) {
  return -1;
}

node *polynomialGetIthCoefficient(polynomial_t p, mpz_t i) {
  return NULL;
}

node *polynomialGetIthCoefficientIntIndex(polynomial_t p, int i) {
  return NULL;
}

int polynomialGetCoefficients(node ***coeffs, int *deg, polynomial_t p) {
  return 0;
}

node *polynomialGetExpression(polynomial_t p, int canonical) {
  return NULL;
}

void polynomialFPrintf(FILE *fd, polynomial_t p, int canonical) {
}

char *polynomialToString(polynomial_t p, int canonical) {
  return NULL;
}

void polynomialEvalMpfr(mpfr_t y, polynomial_t p, mpfr_t x, mp_prec_t prec) {
}

void polynomialEvalMpfi(sollya_mpfi_t y, polynomial_t p, sollya_mpfi_t x, mp_prec_t prec) {
}

