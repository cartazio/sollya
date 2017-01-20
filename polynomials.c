/*

  Copyright 2014-2016 by

  Laboratoire d'Informatique de Paris 6 - Équipe PEQUAN
  Sorbonne Universités
  UPMC Univ Paris 06
  UMR 7606, LIP6
  Boîte Courrier 169
  4, place Jussieu
  F-75252 Paris Cedex 05
  France

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

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "general.h"
#include "execute.h"
#include "infnorm.h"
#include "expression.h"
#include "polynomials.h"
#include "printf.h"
#include "hash.h"

#if (!(defined(HAVE_MP_BITCNT_T) && (HAVE_MP_BITCNT_T)))
typedef unsigned long int mp_bitcnt_t;
#endif


/* Helper types */


struct __hash_result_struct_t {
  uint64_t hash;
  int hasHash;
};
typedef struct __hash_result_struct_t hash_result_t;

struct __boolean_result_cache_struct_t {
  int res;
  int cached;
};
typedef struct __boolean_result_cache_struct_t boolean_result_cache_t;



/* Types for constants, sparse polynomials and polynomials */

typedef enum __constant_type_enum_t constant_type_t;
enum __constant_type_enum_t {
  INTEGER = 0,
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
  boolean_result_cache_t isZero;
  boolean_result_cache_t isOne;
  boolean_result_cache_t isNonNegativeInteger;
  boolean_result_cache_t isPositive;
  boolean_result_cache_t isDyadic;
  boolean_result_cache_t isRational;
  hash_result_t hash;
  union { 
    int integer;
    node *expr;
    mpfr_t mpfr;
    scaled_mpq_t scaledMpq;
  } value;
};

typedef struct __sparse_polynomial_struct_t * sparse_polynomial_t;
struct __sparse_polynomial_struct_t {
  unsigned int refCount;
  constant_t deg;
  unsigned int monomialCount;
  hash_result_t hash;
  constant_t *coeffs;
  constant_t *monomialDegrees;
};

#define SPARSE_POLYNOMIAL_GCD_HEURISTIC_MAX_MONOMIAL_COUNT  ((unsigned int) ((((uint64_t) 1u) << 16) - ((uint64_t) 1u)))
#define SPARSE_POLYNOMIAL_GCD_HEURISTIC_MAX_DEGREE          ((int) ((((int64_t) 1) << 18) - ((int64_t) 1)))
#define SPARSE_POLYNOMIAL_GCD_HEURISTIC_MAX_INTEGER_BOUND   ((mp_bitcnt_t) ((((uint64_t) 1u) << 15) - ((uint64_t) 1u)))
#define SPARSE_POLYNOMIAL_GCD_HEURISTIC_TRIALS              (4)

typedef enum __polynomial_type_enum_t polynomial_type_t;
enum __polynomial_type_enum_t {
  SPARSE = 0,
  ADDITION,
  SUBTRACTION,
  MULTIPLICATION,
  COMPOSITION,
  NEGATE,
  POWER
};

typedef enum __polynomial_output_type_enum_t polynomial_output_type_t;
enum __polynomial_output_type_enum_t {
  ANY_FORM = 0,
  HORNER_FORM,
  CANONICAL_FORM
};

struct __polynomial_struct_t {
  unsigned int refCount;
  polynomial_type_t type;
  polynomial_output_type_t outputType;
  hash_result_t hash;
  boolean_result_cache_t usesExpressionConstant;
  union { 
    sparse_polynomial_t sparse;
    polynomial_t g;
    struct {
      polynomial_t g;
      polynomial_t h;
    } pair;
    struct {
      polynomial_t g;
      constant_t c;
    } powering;
  } value;
};

/* Globally defined constants and caches */

#define CONSTANT_INTEGER_CACHE_MIN (-16384)
#define CONSTANT_INTEGER_CACHE_MAX (16383)
#define CONSTANT_INTEGER_CACHE_SIZE ((CONSTANT_INTEGER_CACHE_MAX) - (CONSTANT_INTEGER_CACHE_MIN) + 1)

#define CONSTANT_MALLOC_CACHE_SIZE (4096)

typedef struct __cached_constant_struct_t __cached_constant_t;
struct __cached_constant_struct_t {
  constant_t constant;
  int initialized;
};

static int __constant_cache_initialized = 0;
static __cached_constant_t __constant_integer_cache[CONSTANT_INTEGER_CACHE_SIZE];

static constant_t __constant_malloc_cache[CONSTANT_MALLOC_CACHE_SIZE];
static int __constant_malloc_cache_index = 0;

/* Helper functions */

static inline int mpfr_is_machine_integer(int *intval, mpfr_t op) {
  long int t, ttt;
  int tt;

  if (!mpfr_number_p(op)) return 0;
  if (!mpfr_integer_p(op)) return 0;
  if (!mpfr_fits_slong_p(op, GMP_RNDN)) return 0;
  t = mpfr_get_si(op, GMP_RNDN); /* exact */
  tt = t;
  ttt = tt;
  if (t != ttt) return 0;
  *intval = tt;
  return 1;
}

static inline int mpfr_is_machine_unsigned_integer(unsigned int *intval, mpfr_t op) {
  unsigned long int t, ttt;
  unsigned int tt;

  if (!mpfr_number_p(op)) return 0;
  if (mpfr_sgn(op) < 0) return 0;
  if (!mpfr_integer_p(op)) return 0;
  if (!mpfr_fits_ulong_p(op, GMP_RNDN)) return 0;
  t = mpfr_get_ui(op, GMP_RNDN); /* exact */
  tt = t;
  ttt = tt;
  if (t != ttt) return 0;
  *intval = tt;
  return 1;
}

/* This function supposes that op is canonicalized, i.e. zero is
   represented as 0/1, that the numerator and the denominator have no
   common factors and that the denominator is positive. The function
   guarantees that op is canonicalized in output as well.
*/
static inline mp_exp_t mpq_remove_powers_of_two(mpq_t op) {
  mp_bitcnt_t dyadNum, dyadDen;
  mp_exp_t expo;

  if (mpq_sgn(op) == 0) return ((mp_exp_t) 0);

  dyadNum = mpz_scan1(mpq_numref(op), 0);
  dyadDen = mpz_scan1(mpq_denref(op), 0);
  mpz_tdiv_q_2exp(mpq_numref(op), mpq_numref(op), dyadNum);
  mpz_tdiv_q_2exp(mpq_denref(op), mpq_denref(op), dyadDen);
  expo = dyadNum - dyadDen;

  return expo;
}

static inline void scaledMpqAdd(mp_exp_t *EC, mpq_t c, 
				mp_exp_t EA, mpq_t a, 
				mp_exp_t EB, mpq_t b) {
  if (EB <= EA) {
    *EC = EB;
    mpq_mul_2exp(c, a, EA - EB); /* Exponent overflow possible */
    mpq_add(c, c, b);
  } else {
    *EC = EA;
    mpq_mul_2exp(c, b, EB - EA); /* Exponent overflow possible */
    mpq_add(c, c, a);
  }
  *EC += mpq_remove_powers_of_two(c); /* Exponent overflow possible */
}

static inline void scaledMpqAddInt(mp_exp_t *EC, mpq_t c, 
				   mp_exp_t EA, mpq_t a, 
				   int b) {
  mpq_t t;

  /* Can be optimized */
  mpq_init(t);
  mpq_set_si(t, b, 1u);
  mpq_canonicalize(t);
  scaledMpqAdd(EC, c, EA, a, 0, t);
  mpq_clear(t);
}

static inline void scaledMpqSub(mp_exp_t *EC, mpq_t c, 
				mp_exp_t EA, mpq_t a, 
				mp_exp_t EB, mpq_t b) {
  if (EB <= EA) {
    *EC = EB;
    mpq_mul_2exp(c, a, EA - EB); /* Exponent overflow possible */
    mpq_sub(c, c, b);
  } else {
    *EC = EA;
    mpq_mul_2exp(c, b, EB - EA); /* Exponent overflow possible */
    mpq_sub(c, a, c);
  }
  *EC += mpq_remove_powers_of_two(c); /* Exponent overflow possible */
}

static inline void scaledMpqSubInt(mp_exp_t *EC, mpq_t c, 
				   mp_exp_t EA, mpq_t a, 
				   int b) {
  mpq_t t;

  /* Can be optimized */
  mpq_init(t);
  mpq_set_si(t, b, 1u);
  mpq_canonicalize(t);
  scaledMpqSub(EC, c, EA, a, 0, t);
  mpq_clear(t);
}

static inline void scaledMpqIntSub(mp_exp_t *EC, mpq_t c, 
				   int a,
				   mp_exp_t EB, mpq_t b) {
  mpq_t t;

  /* Can be optimized */
  mpq_init(t);
  mpq_set_si(t, a, 1u);
  mpq_canonicalize(t);
  scaledMpqSub(EC, c, 0, t, EB, b);
  mpq_clear(t);
}

static inline void scaledMpqMul(mp_exp_t *EC, mpq_t c, 
				mp_exp_t EA, mpq_t a, 
				mp_exp_t EB, mpq_t b) {
  *EC = EA + EB; /* Exponent overflow possible */
  mpq_mul(c, a, b);
  *EC += mpq_remove_powers_of_two(c); /* Exponent overflow possible */
}

static inline void scaledMpqMulInt(mp_exp_t *EC, mpq_t c, 
				   mp_exp_t EA, mpq_t a, 
				   int b) {
  mpq_t t;

  /* Can be optimized */
  mpq_init(t);
  mpq_set_si(t, b, 1u);
  mpq_canonicalize(t);
  scaledMpqMul(EC, c, EA, a, 0, t);
  mpq_clear(t);
}

static inline int sollya_mpz_pow(mpz_t z, mpz_t x, mpz_t y) {
  unsigned long int k;
  int sgn, res;
  mpz_t yhi, ylo, xyhi, xylo, xt14;

  /* Get sign of y */
  sgn = mpz_sgn(y);

  /* Refuse work for negative y */
  if (sgn < 0) return 0;

  /* If y is zero, set z to 1 */
  if (sgn == 0) {
    mpz_set_ui(z, 1u);
    return 1;
  }
  
  /* Try to see if y fits into a ulong */
  if (mpz_fits_ulong_p(y)) {
    k = mpz_get_ui(y);
    mpz_pow_ui(z, x, k);
    return 1;
  }

  /* Decompose y into y = yhi * 2^14 + ylo */
  mpz_init(yhi);
  mpz_init(ylo);
  mpz_fdiv_qr_ui(yhi, ylo, y, (1u << 14));

  /* Compute x^y as 

     x^y = x^(2^14 * yhi + ylo) = (x^(2^14))^yhi * x^ylo.

  */
  mpz_init(xylo);
  mpz_init(xyhi);
  mpz_init(xt14);
  if (sollya_mpz_pow(xylo, x, ylo)) {
    mpz_pow_ui(xt14, x, (1u << 14));
    res = sollya_mpz_pow(xyhi, xt14, yhi); 
    mpz_mul(z, xyhi, xylo);
  } else {
    res = 0;
  }

  /* Clear temporaries */
  mpz_clear(xt14);
  mpz_clear(xyhi);
  mpz_clear(xylo);
  mpz_clear(ylo);
  mpz_clear(yhi);

  /* Return success indicator */
  return res;
}

static inline int tryScaledMpqDiv(mp_exp_t *EC, mpq_t c, 
				  mp_exp_t EA, mpq_t a, 
				  mp_exp_t EB, mpq_t b) {
  if (mpq_sgn(b) == 0) return 0;
  *EC = EA - EB; /* Exponent overflow possible */
  mpq_div(c, a, b);
  *EC += mpq_remove_powers_of_two(c); /* Exponent overflow possible */
  return 1;
}

static inline void sollya_mpq_gcd(mpq_t c, mpq_t a, mpq_t b) {
  mpz_t u, v;

  mpz_init(u);
  mpz_init(v);
  mpz_mul(mpq_denref(c), mpq_denref(a), mpq_denref(b));
  mpz_mul(u, mpq_numref(a), mpq_denref(b));
  mpz_mul(v, mpq_denref(a), mpq_numref(b));
  mpz_gcd(mpq_numref(c), u, v);
  mpz_clear(v);
  mpz_clear(u);
  mpq_canonicalize(c);
}

static inline int scaledMpqIsInteger(mp_exp_t E, mpq_t a) {
  mpq_t aa;
  mp_exp_t EE;
  int res;

  /* 0 is an integer */
  if (mpq_sgn(a) == 0) return 1;

  /* If the denominator q is one in magnitude and the exponent E
     non-negative, 2^E * p/q is integer
  */
  if ((E >= 0) && 
      (mpz_cmpabs_ui(mpq_denref(a),1u) == 0)) return 1;

  /* If the denominator q is one in magnitude and the exponent E
     non-negative, 2^E * p/q is integer
  */
  if ((E >= 0) && 
      (mpz_cmpabs_ui(mpq_denref(a),1u) == 0)) return 1;

  /* We can take short-cuts if the numerator p and the denominator q
     are both odd. Check that condition.
  */
  if (mpz_odd_p(mpq_numref(a)) &&
      mpz_odd_p(mpq_denref(a))) {
    /* Here p and q are co-prime and odd. If E is negative, 2^E * p/q
       is no integer. Otherwise 2^E * p/q is an integer iff abs(p) ==
       abs(q).
    */
    if (E < 0) return 0;
    return (mpz_cmpabs(mpq_numref(a),mpq_denref(a)) == 0);
  }

  /* Here we know that at least one of the numerator and denominator
     has 2 as a prime factor. 
  */
  mpq_init(aa);
  mpq_set(aa, a);
  EE = E;

  mpq_canonicalize(aa);
  EE += mpq_remove_powers_of_two(aa); /* Exponent overflow possible */

  /* Now we have 2^EE * aa = 2^E * a, aa in least factors and all
     prime factors 2 taken out. 
  */
  if (EE < 0) {
    mpq_clear(aa);
    return 0;
  }

  res = (mpz_cmpabs(mpq_numref(aa),mpq_denref(aa)) == 0);
  mpq_clear(aa);
  
  return res;
}

static inline int scaledMpqIsDyadic(mp_exp_t E, mpq_t a) {
  mpq_t aa;
  int res;

  UNUSED_PARAM(E);

  /* 0 is a dyadic */
  if (mpq_sgn(a) == 0) return 1;

  /* If the denominator q is one in magnitude, 2^E * p/q is dyadic. */
  if (mpz_cmpabs_ui(mpq_denref(a),1u) == 0) return 1;
  
  /* We can take short-cuts if the numerator p and the denominator q
     are both odd. Check that condition.
  */
  if (mpz_odd_p(mpq_numref(a)) &&
      mpz_odd_p(mpq_denref(a))) {
    /* Here p and q are co-prime and odd and q is different from 1. 
       The number is not dyadic.
    */
    return 0;
  }
  
  /* Here we know that at least one of the numerator and denominator
     has 2 as a prime factor. 
  */
  mpq_init(aa);
  mpq_set(aa, a);
  mpq_canonicalize(aa);
  mpq_remove_powers_of_two(aa);

  /* Now we have 2^EE * aa = 2^E * a, aa in least factors and all
     prime factors 2 taken out. 

     2^EE * pp/qq with pp and qq odd, co-prime. The number is 
     dyadic iff qq is one in magnitude.
  */
  res = (mpz_cmpabs_ui(mpq_denref(a),1u) == 0);
  mpq_clear(aa);
  return res;
}

/* Determine if 2^EA * a > 2^EB * b */
static inline int scaledMpqIsGreaterThan(mp_exp_t EA, mpq_t a,
					 mp_exp_t EB, mpq_t b) {
  mp_exp_t E, Emin, Emax, F, D;
  mpz_t p, q;
  mpq_t r;
  int res;
  
  /* Handle zero in input */
  if (mpq_sgn(a) == 0) {
    return (mpq_sgn(b) < 0);
  }
  if (mpq_sgn(b) == 0) {
    return (mpq_sgn(a) > 0);
  }

  /* Handle inputs of different signs */
  if (mpq_sgn(a) * mpq_sgn(b) < 0) {
    return (mpq_sgn(a) > 0);
  }

  /* Here 2^EA * a and 2^EB * b have the same sign 

     The comparison 2^EA * a > 2^EB * b is equivalent to

     *    2^(EA - EB) > b/a      if a is positive
     *    2^(EA - EB) < b/a      otherwise

     We start by computing Emin and Emax such that

     2^Emin < b/a < 2^Emax

     This way we can then perform an easy, first check
     based on the order of magnitudes.

  */
  E = (mp_exp_t) ((mpz_sizeinbase(mpq_numref(b),2) + 
		   mpz_sizeinbase(mpq_denref(a),2)) - 
		  (mpz_sizeinbase(mpq_denref(b),2) + 
		   mpz_sizeinbase(mpq_numref(a),2))); /* Exponent overflow possible */
  Emin = E - 2; /* Exponent overflow possible */
  Emax = E + 2; /* Exponent overflow possible */

  /* Perform an easy check based on the magnitudes, taking the sign of
     a into account. 
  */
  if (mpq_sgn(a) > 0) {
    /* a is positive, so we have to perform the comparison

       2^(EA - EB) > b/a 

       knowing that 2^Emin < b/a < 2^Emax.

       If EA - EB >= Emax, we have 
       
       2^(EA - EB) >= 2^Emax > b/a

       hence the answer surely is "yes, 2^(EA - EB) > b/a is
       satisfied"

       If EA - EB <= Emin, we have

       2^(EA - EB) <= 2^Emin < b/a

       hence the answer surely is "no, 2^(EA - EB) > b/a is not
       satisfied"
       
    */
    if (EA - EB >= Emax) return 1; /* Exponent overflow possible */
    if (EA - EB <= Emin) return 0; /* Exponent overflow possible */
  } else {
    /* a is positive, so we have to perform the comparison

       2^(EA - EB) < b/a 

       knowing that 2^Emin < b/a < 2^Emax.

       If EA - EB >= Emax, we have 
       
       2^(EA - EB) >= 2^Emax > b/a

       hence the answer surely is "no, 2^(EA - EB) < b/a is not
       satisfied"

       If EA - EB <= Emin, we have

       2^(EA - EB) <= 2^Emin < b/a

       hence the answer surely is "yes, 2^(EA - EB) < b/a is
       satisfied"
       
    */
    if (EA - EB >= Emax) return 0; /* Exponent overflow possible */
    if (EA - EB <= Emin) return 1; /* Exponent overflow possible */
  }

  /* Here, b/a is too close to 2^E to allow the comparison to be
     decided merely on the base of the order of magnitude.

     However, now it is reasonable to compute the ratio

     b/a = 2^F * p/q 

     where p and q are co-prime odd integers and q is positive.

     Then the comparison between 2^(EA - EB) and b/a 
     becomes a comparison between

     2^(EA - EB - F) * q and p.

     Let D = EA - EB - F

  */
  mpq_init(r);
  mpq_div(r, b, a);
  F = mpq_remove_powers_of_two(r);
  mpz_init(p);
  mpz_init(q);
  mpq_get_num(p, r);
  mpq_get_den(q, r);
  mpq_clear(r);

  /* Compute D = EA - EB - F */
  D = EA - EB - F; /* Exponent overflow possible */

  /* Now integrate 2^D into q or p, such that the comparison boils down
     to comparing p and q. 
  */
  if (D >= 0) {
    mpz_mul_2exp(q, q, (mp_bitcnt_t) D);
  } else {
    mpz_mul_2exp(p, p, (mp_bitcnt_t) (-D));
  }

  if (mpq_sgn(a) > 0) {
    /* a is positive. We have to perform the comparison

       2^(EA - EB) > b/a

       which is equivalent to 

       q > p.

    */
    res = (mpz_cmp(q, p) > 0);
    mpz_clear(p);
    mpz_clear(q);
    return res;
  } 

  /* a is negative. We have to perform the comparison

     2^(EA - EB) < b/a
     
     which is equivalent to 

     q < p.

  */
  res = (mpz_cmp(q, p) < 0);
  mpz_clear(p);
  mpz_clear(q);
  return res;
}

static inline void scaledMpqFloor(mp_exp_t *EB, mpq_t b,
				  mp_exp_t EA, mpq_t a) {
  mp_exp_t sizeNum, sizeDen, sizeOut, ER, ERP, Erest, Edelta;
  mp_prec_t precOut;
  mpfr_t t;
  mpq_t r, rPrime, one, rest, delta;
  mpz_t z, num, den;
  
  /* If 2^EA * a is an integer, we have nothing to do */
  if (scaledMpqIsInteger(EA, a)) {
    *EB = EA;
    mpq_set(b, a);
    return;
  }
  
  /* Here 2^EA * a is not an integer and hence non-zero 

     Get approximate sizes of the numerator and the denominator in
     bits.

  */
  sizeNum = (mp_exp_t) (mpz_sizeinbase(mpq_numref(a),2) + 1);
  sizeDen = (mp_exp_t) (mpz_sizeinbase(mpq_denref(a),2));

  /* Compute an upper bound on the size of the output */
  sizeOut = EA + sizeNum + 5 - sizeDen; /* Integer overflow possible */
  if (sizeOut < 12) sizeOut = 12;
 
  /* This upper bound on the size of the output gives us 
     an idea of the precision we need to compute a first guess 
     on the output.
  */
  precOut = sizeOut;
  if (precOut < 12) precOut = 12;

  /* Compute a first guess on the output */
  mpfr_init2(t, precOut);
  mpfr_set_z_2exp(t, mpq_numref(a), EA, GMP_RNDD);   
  mpfr_div_z(t, t, mpq_denref(a), GMP_RNDD);
  mpfr_floor(t, t);
  
  /* Represent the first guess as a scaled MPQ */
  mpq_init(r);
  mpz_init(z);
  if (mpfr_zero_p(t)) {
    ER = EA;
    mpq_set_si(r, 0, 1u);
  } else {
    ER = mpfr_get_z_2exp(z, t);
    mpq_set_z(r,z);
  }
  mpfr_clear(t);
  mpz_clear(z);
  mpq_canonicalize(r);
  ER += mpq_remove_powers_of_two(r); /* Exponent overflow possible */

  /* As we computed the first guess in round-down mode, we are sure
     that 2^ER * r < 2^EA * a.

     We now check if

     2^ER * r > 2^EA * a - 1,

     i.e. if 

     2^ER * r + 1 > 2^EA * a
     
     This condition should be satisfied in almost all cases. If it is
     not, the first guess was not accurate enough.

     Representing 2^ER * r + 1 is in general as memory intensive as 
     representing 2^ER * r.

  */
  mpq_init(rPrime);
  mpq_init(one);
  mpq_set_si(one, 1, 1u);
  scaledMpqAdd(&ERP, rPrime, ER, r, 0, one);
  mpq_clear(one);  

  if (scaledMpqIsGreaterThan(ERP, rPrime, EA, a)) {
    mpq_clear(rPrime);
    *EB = ER;
    mpq_set(b,r);
    mpq_clear(r);
    return;
  }
  mpq_clear(rPrime);
    
  /* The accuracy of the first guess was not enough. 
     We compute the rest 

     2^Erest * rest = 2^EA * a - 2^ER * r

     then we compute floor(2^Erest * rest) in MPQ/MPZ arithmetic.

  */
  mpq_init(rest);
  scaledMpqSub(&Erest, rest, EA, a, ER, r);
  mpz_init(num);
  mpz_init(den);
  mpq_get_num(num, rest);
  mpq_get_den(den, rest);
  mpq_clear(rest);
  if (Erest > 0) {
    mpz_mul_2exp(num, num, (mp_bitcnt_t) Erest);
  } else {
    mpz_mul_2exp(den, den, (mp_bitcnt_t) (-Erest)); /* Exponent overflow possible */
  }
  
  /* Now compute floor(num/den). That's a Euclidean division. */
  mpz_fdiv_q(num,num,den);
  mpz_clear(den);
  
  /* Represent the correction as a scaled MPQ */
  mpq_init(delta);
  mpq_set_z(delta, num);
  mpq_canonicalize(delta);
  Edelta = mpq_remove_powers_of_two(delta);
  mpz_clear(num);

  /* The result is 2^ER * r + 2^Edelta * delta */
  scaledMpqAdd(EB, b, ER, r, Edelta, delta);
  mpq_clear(r);
  mpq_clear(delta);
}

/* Checks if 2^E * a is an integer that can be represented as an unsigned int and 
   computes that integer.
*/
static inline int scaledMpqIsUnsignedInt(unsigned int *n, mp_exp_t E, mpq_t a) {
  mpfr_t t;
  unsigned long tI, ttI, uIAbs;
  unsigned int tt;
  signed long uI;
  mpz_t r, tmp, u, v;

  /* First check for zero */
  if (mpq_sgn(a) == 0) {
    *n = 0u;
    return 1;
  }

  /* Then check if 2^E * a is a non-negative integer */
  if (!scaledMpqIsInteger(E,a)) return 0;
  if (mpq_sgn(a) < 0) return 0;
  
  /* Here, 2^E * a is a positive integer 

     Compute a first guess on its integer value.

     We compute t = floor(2^E * p/q * (1 + eps))

     with abs(eps) <= 2 * e + e^2

     where abs(e) <= 2^(2 * -l - 7 + 1)

     where l is the number of bits in an unsigned int.

  */
  mpfr_init2(t, 2 * 8 * sizeof(unsigned int) + 7);
  mpfr_set_z_2exp(t, mpq_numref(a), E, GMP_RNDD);   
  mpfr_div_z(t, t, mpq_denref(a), GMP_RNDD);
  mpfr_floor(t, t);
  if (!mpfr_fits_ulong_p(t, GMP_RNDN)) {
    mpfr_clear(t);
    return 0;
  }
  tI = mpfr_get_ui(t, GMP_RNDN); /* exact */
  mpfr_clear(t);
  
  /* Now compute the "rest" 


                            / (2^E * p - q * t) / q  if E >= 0
     r/q = 2^E * p/q - t = | 
                            \ (p - 2^-E * q * t) / q if E < 0

     We can be sure that this computation is not too 
     memory intensive because we know that t is bounded.

  */
  mpz_init(r);
  if (E >= 0) {
    mpz_set(r, mpq_numref(a));
    mpz_mul_2exp(r, r, (mp_bitcnt_t) E);
    mpz_init(tmp);
    mpz_mul_ui(tmp, mpq_denref(a), tI);
    mpz_sub(r, r, tmp);
    mpz_clear(tmp);
  } else {
    mpz_set(r, mpq_denref(a));
    mpz_mul_2exp(r, r, (mp_bitcnt_t) (-E)); /* Exponent overflow possible */
    mpz_mul_ui(r, r, tI);
    mpz_sub(r, mpq_numref(a), r);
  }
  
  /* If the rest is zero, we know the integer value 

     tI = 2^E * a

  */
  if (mpz_sgn(r) == 0) {
    mpz_clear(r);
    tt = tI;
    ttI = tt;
    if (ttI != tI) return 0;
    *n = tt;
    return 1;
  }

  /* Here, the rest r/q is not zero. We perform 
     an Euclidean division.

     u = floor(r/q), q * u + v = r

  */
  mpz_init(u);
  mpz_init(v);
  mpz_tdiv_qr(u, v, r, mpq_denref(a));
  
  /* If v is non-zero, 2^E * a is no integer */
  if (mpz_sgn(v) != 0) {
    mpz_clear(v);
    mpz_clear(u);
    mpz_clear(r);
    return 0;
  }
  mpz_clear(v);
  mpz_clear(r);

  /* Here we know that u + tI = 2^E * a 

     We start by checking if u fits a signed long.

  */
  if (!mpz_fits_slong_p(u)) {
    mpz_clear(u);
    return 0;
  }
  uI = mpz_get_si(u);
  mpz_clear(u);

  /* Here we know that uI + tI = 2^E * a 

     where uI is a signed long and tI is an unsigned long.

  */
  if (uI < 0l) {
    uIAbs = -uI;
    if (uIAbs < tI) return 0;
    tI -= uIAbs;
    tt = tI;
    ttI = tt;
    if (ttI != tI) return 0;
    *n = tt;
    return 1;
  } 

  /* Here we know that uI + tI = 2^E * a and uI >= 0 */
  uIAbs = uI;
  tI += uIAbs;
  if (tI < uIAbs) return 0;
  tt = tI;
  ttI = tt;
  if (ttI != tI) return 0;
  *n = tt;
  return 1;
}

/* Checks if 2^E * a is an integer that can be represented as an int and 
   computes that integer.
*/
static inline int scaledMpqIsInt(int *n, mp_exp_t E, mpq_t a) {
  unsigned int nU, nUU, rMA;
  int r, rM;
  mpq_t aAbs;

  if (mpq_sgn(a) >= 0) {
    if (!scaledMpqIsUnsignedInt(&nU, E, a)) return 0;
    r = nU;
    if (r < 0) return 0;
    nUU = r;
    if (nUU != nU) return 0;
    *n = r;
    return 1;
  }

  mpq_init(aAbs);
  mpq_neg(aAbs, a);
  if (!scaledMpqIsUnsignedInt(&nU, E, aAbs)) {
    mpq_clear(aAbs);
    return 0;
  }
  mpq_clear(aAbs);
  r = nU;
  if (r < 0) return 0;
  nUU = r;
  if (nUU != nU) return 0;
  rM = -r;
  if (rM > 0) return 0;
  rMA = -rM;
  if (rMA != nU) return 0;
  *n = rM;
  return 1;
}

static inline int tryScaledMpqPow(mp_exp_t *EC, mpq_t c, 
				  mp_exp_t EA, mpq_t a, 
				  mp_exp_t EB, mpq_t b) {
  unsigned int n;
  mpz_t eaZ, nZ, num, den;
  signed long int eaI, ecI, ecII;
  mp_exp_t eaE, ecR;

  if (!scaledMpqIsUnsignedInt(&n, EB, b)) return 0;
  if (n == 0u) {
    *EC = 0;
    mpq_set_ui(c,1,1u);
    mpq_canonicalize(c);
    *EC += mpq_remove_powers_of_two(c); /* Exponent overflow possible */
    return 1;
  }
  if (n == 1u) {
    *EC = EA;
    mpq_set(c, a);
    mpq_canonicalize(c);
    *EC += mpq_remove_powers_of_two(c); /* Exponent overflow possible */
    return 1;
  }

  eaI = EA;
  eaE = eaI;
  if (eaE != EA) return 0;
  mpz_init_set_si(eaZ, eaI);
  mpz_init_set_ui(nZ, n);
  mpz_mul(eaZ, eaZ, nZ);
  if (!mpz_fits_slong_p(eaZ)) {
    mpz_clear(nZ);
    mpz_clear(eaZ);
    return 0;
  }
  ecI = mpz_get_si(eaZ);
  mpz_clear(nZ);
  mpz_clear(eaZ);
  ecR = ecI;
  ecII = ecR;
  if (ecII != ecI) return 0;
  
  mpz_init(num);
  mpz_init(den);

  mpq_get_num(num, a);
  mpq_get_den(den, a);

  mpz_pow_ui(num, num, n);
  mpz_pow_ui(den, den, n);

  *EC = ecR;
  mpq_set_num(c, num);
  mpq_set_den(c, den);
  mpq_canonicalize(c);
  *EC += mpq_remove_powers_of_two(c); /* Exponent overflow possible */

  mpz_clear(den);
  mpz_clear(num);

  return 1;
}

static inline int tryScaledMpqPowInt(mp_exp_t *EC, mpq_t c, 
				     mp_exp_t EA, mpq_t a, 
				     int b) {
  unsigned int n;
  mpz_t eaZ, nZ, num, den;
  signed long int eaI, ecI, ecII;
  mp_exp_t eaE, ecR;
  int bb;

  if (b < 0) return 0;
  n = b;
  bb = n;
  if (bb < 0) return 0;
  if (bb != b) return 0;

  if (n == 0u) {
    *EC = 0;
    mpq_set_ui(c,1,1u);
    mpq_canonicalize(c);
    *EC += mpq_remove_powers_of_two(c); /* Exponent overflow possible */
    return 1;
  }
  if (n == 1u) {
    *EC = EA;
    mpq_set(c, a);
    mpq_canonicalize(c);
    *EC += mpq_remove_powers_of_two(c); /* Exponent overflow possible */
    return 1;
  }

  eaI = EA;
  eaE = eaI;
  if (eaE != EA) return 0;
  mpz_init_set_si(eaZ, eaI);
  mpz_init_set_ui(nZ, n);
  mpz_mul(eaZ, eaZ, nZ);
  if (!mpz_fits_slong_p(eaZ)) {
    mpz_clear(nZ);
    mpz_clear(eaZ);
    return 0;
  }
  ecI = mpz_get_si(eaZ);
  mpz_clear(nZ);
  mpz_clear(eaZ);
  ecR = ecI;
  ecII = ecR;
  if (ecII != ecI) return 0;
  
  mpz_init(num);
  mpz_init(den);

  mpq_get_num(num, a);
  mpq_get_den(den, a);

  mpz_pow_ui(num, num, n);
  mpz_pow_ui(den, den, n);

  *EC = ecR;
  mpq_set_num(c, num);
  mpq_set_den(c, den);
  mpq_canonicalize(c);
  *EC += mpq_remove_powers_of_two(c); /* Exponent overflow possible */

  mpz_clear(den);
  mpz_clear(num);

  return 1;
}

static inline int tryExactIntAddition(int *c, int a, int b) {
  int aNeg, bNeg, cNeg, cT;
  unsigned int aAbs, bAbs, cAbs, cTAbs;

  if (a >= 0) {
    aAbs = a;
    aNeg = 0;
  } else {
    aAbs = -a;
    aNeg = 1;
  }
  if (b >= 0) {
    bAbs = b;
    bNeg = 0;
  } else {
    bAbs = -b;
    bNeg = 1;
  }
  if (aNeg == bNeg) {
    cNeg = aNeg;
    cAbs = aAbs + bAbs;
    if (cAbs < aAbs) return 0; /* Should never happen */
  } else {
    if (aAbs >= bAbs) {
      cNeg = aNeg;
      cAbs = aAbs - bAbs;
    } else {
      cNeg = bNeg;
      cAbs = bAbs - aAbs;
    }
  }
  if (cNeg) {
    cT = -cAbs;
    if (cT > 0) return 0;
    cTAbs = -cT;
    if (cTAbs != cAbs) return 0;
  } else {
    cT = cAbs;
    if (cT < 0) return 0;
    cTAbs = cT;
    if (cTAbs != cAbs) return 0;
  }
  *c = cT;
  return 1;
}

static inline int tryExactIntSubtraction(int *c, int a, int b) {
  int aNeg, bNeg, cNeg, cT;
  unsigned int aAbs, bAbs, cAbs, cTAbs;

  if (a >= 0) {
    aAbs = a;
    aNeg = 0;
  } else {
    aAbs = -a;
    aNeg = 1;
  }
  if (b >= 0) {
    bAbs = b;
    bNeg = 1;
  } else {
    bAbs = -b;
    bNeg = 0;
  }
  if (aNeg == bNeg) {
    cNeg = aNeg;
    cAbs = aAbs + bAbs;
    if (cAbs < aAbs) return 0; /* Should never happen */
  } else {
    if (aAbs >= bAbs) {
      cNeg = aNeg;
      cAbs = aAbs - bAbs;
    } else {
      cNeg = bNeg;
      cAbs = bAbs - aAbs;
    }
  }
  if (cNeg) {
    cT = -cAbs;
    if (cT > 0) return 0;
    cTAbs = -cT;
    if (cTAbs != cAbs) return 0;
  } else {
    cT = cAbs;
    if (cT < 0) return 0;
    cTAbs = cT;
    if (cTAbs != cAbs) return 0;
  }
  *c = cT;
  return 1;
}

static inline void exactUint64Mul(uint64_t *ch, 
				  uint64_t *cl, 
				  uint64_t a, 
				  uint64_t b) {
  uint64_t ah, al, bh, bl;
  uint64_t ahbh, ahbl, albh, albl;
  uint64_t ahblh, ahbll, albhh, albhl;
  uint64_t rh, rl;

  ah = a >> 32;        bh = a >> 32;
  al = a - (ah << 32); bl = b - (bh << 32);

  ahbh = ah * bh; ahbl = ah * bl;
  albh = al * bh; albl = al * bl;

  ahblh = ahbl >> 32; ahbll = ahbl - (ahblh << 32);
  albhh = albh >> 32; albhl = albh - (albhh << 32);

  ahbll <<= 32;
  albhl <<= 32;

  rl = albl + ahbll;
  if (rl < ahbll) ahblh++;
  
  rl += albhl;
  if (rl < albhl) ahblh++;

  rh = ahbh + ahblh + albhh;

  *ch = rh;
  *cl = rl;
}

static inline int tryExactUint64Multiplication(uint64_t *c, 
					       uint64_t a, 
					       uint64_t b) {
  uint64_t cth, ctl;

  exactUint64Mul(&cth, &ctl, a, b);
  
  if (cth == UINT64_C(0)) {
    *c = ctl;
    return 1;
  }
  return 0;
}

static inline int tryExactUnsignedIntMultiplication(unsigned int *c, 
						    unsigned int a, 
						    unsigned int b) {
  uint64_t a64, b64, c64, cT64;
  unsigned int aa, bb, cT;

  a64 = a;
  aa = a64;
  if (aa != a) return 0; /* should never happen */

  b64 = b;
  bb = b64;
  if (bb != b) return 0; /* should never happen */
  
  if (!tryExactUint64Multiplication(&c64, a64, b64)) return 0;

  cT = c64;
  cT64 = cT;
  if (c64 != cT64) return 0;

  *c = cT;
  return 1;
}

static inline int tryExactIntMultiplication(int *c, int a, int b) {
  int aNeg, bNeg, cNeg;
  unsigned int aAbs, bAbs, cAbs, cTAbs;
  int cT;
  
  if (a >= 0) {
    aAbs = a;
    aNeg = -1;
  } else {
    aAbs = -a;
    aNeg = 1;
  }
  if (b >= 0) {
    bAbs = b;
    bNeg = -1;
  } else {
    bAbs = -b;
    bNeg = 1;
  }
  cNeg = aNeg * bNeg;
  if (!tryExactUnsignedIntMultiplication(&cAbs, aAbs, bAbs)) return 0;
  if (cNeg < 0) {
    cT = -cAbs;
    if (cT > 0) return 0;
    cTAbs = -cT;
    if (cTAbs != cAbs) return 0;
  } else {
    cT = cAbs;
    if (cT < 0) return 0;
    cTAbs = cT;
    if (cTAbs != cAbs) return 0;
  }
  *c = cT;
  return 1;
}

static inline int tryExactUnsignedIntDivision(unsigned int *c, 
					      unsigned int a, 
					      unsigned int b) {
  if (b == 0u) return 0;
  if (a % b != 0u) return 0;
  *c = a / b;
  return 1;
}

static inline int tryExactIntDivision(int *c, int a, int b) {
  int aNeg, bNeg, cNeg;
  unsigned int aAbs, bAbs, cAbs, cTAbs;
  int cT;
  
  if (a >= 0) {
    aAbs = a;
    aNeg = -1;
  } else {
    aAbs = -a;
    aNeg = 1;
  }
  if (b >= 0) {
    bAbs = b;
    bNeg = -1;
  } else {
    bAbs = -b;
    bNeg = 1;
  }
  cNeg = aNeg * bNeg;
  if (!tryExactUnsignedIntDivision(&cAbs, aAbs, bAbs)) return 0;
  if (cNeg < 0) {
    cT = -cAbs;
    if (cT > 0) return 0;
    cTAbs = -cT;
    if (cTAbs != cAbs) return 0;
  } else {
    cT = cAbs;
    if (cT < 0) return 0;
    cTAbs = cT;
    if (cTAbs != cAbs) return 0;
  }
  *c = cT;
  return 1;
}

static inline int tryExactIntNegation(int *c, int a) {
  int r, aa;

  /* Handle zero */
  if (a == 0) {
    *c = a;
    return 1;
  }

  /* Try negation */
  r = -a;
  
  /* If r is zero or r and a have the same sign, something went
     wrong 
  */
  if (r == 0) return 0;
  if (!((r > 0) ^ (a > 0))) return 0;

  /* Now we know that r and a have opposite signs. Handle now clamping. */
  aa = -r;

  /* If aa is zero or aa and r have the same sign, something went
     wrong 
  */
  if (aa == 0) return 0;
  if (!((aa > 0) ^ (r > 0))) return 0;
  
  /* Now we know that aa and a have the same sign. 

     If aa is not equal to a, there has been some clamping.
     Otherwise, r is the opposite of a.
  */
  if (aa != a) return 0;

  *c = r;
  return 1;
}

static inline int tryExactIntPower(int *c, int a, int b) {
  unsigned int n;
  int r, t;

  if (b < 0) return 0;
  if (b == 0) {
    *c = 1;
    return 1;
  }
  n = b;

  t = a;
  r = 1;
  while (n > 0u) {
    if (n & 1u) {
      if (!tryExactIntMultiplication(&r,r,t)) return 0;
    }
    n >>= 1;
    if (n > 0u) {
      if (!tryExactIntMultiplication(&t,t,t)) return 0;
    }
  }

  *c = r;
  return 1;
}

/* Computes c = a + b exactly, adapting the precision of c */
static inline void mpfr_add_exact(mpfr_t c, mpfr_t a, mpfr_t b) {
  mp_exp_t Ea, Eb;
  mp_prec_t pa, pb, pc;
  

  if ((!mpfr_number_p(a)) ||
      (!mpfr_number_p(b))) {
    mpfr_add(c, a, b, GMP_RNDN); /* Producing NaN or Inf */
    return;
  }

  if (mpfr_zero_p(a)) {
    mpfr_set_prec(c, mpfr_get_prec(b));
    mpfr_set(c, b, GMP_RNDN); /* exact */
    return;
  }

  if (mpfr_zero_p(b)) {
    mpfr_set_prec(c, mpfr_get_prec(a));
    mpfr_set(c, a, GMP_RNDN); /* exact */
    return;
  }

  /* Here, both inputs are non-zero numbers */
  
  Ea = mpfr_get_exp(a);
  Eb = mpfr_get_exp(b);
  pa = mpfr_get_prec(a);
  pb = mpfr_get_prec(b);

  /* The ulps are the following

     ulp_pa(a) = 2^(Ea - pa)
     ulp_pb(b) = 2^(Eb - pb)

     So the ulp of the sum must be

     ulp_pc(c) = 2^min(Ea - pa, Eb - pb)
     
     The exponent of the result c is bounded
     by

     Ec <= max(Ea, Eb) + 1 

     The + 1 comes from the possible outgoing carry.

     So we need 

     pc = max(Ea, Eb) - min(Ea - pa, Eb - pb) + 1 bits
     
     for c.

  */
  pc = (Ea > Eb? Ea : Eb) - (Ea - pa < Eb - pb ? Ea - pa : Eb - pb) + 1; /* Exponent overflow possible */
  mpfr_set_prec(c, pc);

  /* Perform the addition */
  mpfr_add(c, a, b, GMP_RNDN); /* exact */

  /* It could be that we actually don't need all the precision of c to
     keep the result. So we reduce that precision to the minimum.
  */
  pc = mpfr_min_prec(c);
  if (pc < 12) pc = 12;
  
  /* Adjust the precision of c again. */
  mpfr_prec_round(c, pc, GMP_RNDN); /* exact */
}

/* Computes c = a + b exactly, adapting the precision of c */
static inline void mpfr_add_exact_int(mpfr_t c, mpfr_t a, int b) {
  mp_exp_t Ea, Eb;
  mp_prec_t pa, pb, pc;
  

  if (!mpfr_number_p(a)) {
    mpfr_add_si(c, a, b, GMP_RNDN); /* Producing NaN or Inf */
    return;
  }

  if (mpfr_zero_p(a)) {
    mpfr_set_prec(c, 8 * sizeof(int) + 5);
    mpfr_set_si(c, b, GMP_RNDN); /* exact */
    return;
  }

  if (b == 0) {
    mpfr_set_prec(c, mpfr_get_prec(a));
    mpfr_set(c, a, GMP_RNDN); /* exact */
    return;
  }

  /* Here, both inputs are non-zero numbers */
  
  Ea = mpfr_get_exp(a);
  Eb = 8 * sizeof(int);
  pa = mpfr_get_prec(a);
  pb = 8 * sizeof(int);

  /* The ulps are the following

     ulp_pa(a) = 2^(Ea - pa)
     ulp_pb(b) = 2^(Eb - pb)

     So the ulp of the sum must be

     ulp_pc(c) = 2^min(Ea - pa, Eb - pb)
     
     The exponent of the result c is bounded
     by

     Ec <= max(Ea, Eb) + 1 

     The + 1 comes from the possible outgoing carry.

     So we need 

     pc = max(Ea, Eb) - min(Ea - pa, Eb - pb) + 1 bits
     
     for c.

  */
  pc = (Ea > Eb? Ea : Eb) - (Ea - pa < Eb - pb ? Ea - pa : Eb - pb) + 1; /* Exponent overflow possible */
  mpfr_set_prec(c, pc);

  /* Perform the addition */
  mpfr_add_si(c, a, b, GMP_RNDN); /* exact */

  /* It could be that we actually don't need all the precision of c to
     keep the result. So we reduce that precision to the minimum.
  */
  pc = mpfr_min_prec(c);
  if (pc < 12) pc = 12;
  
  /* Adjust the precision of c again. */
  mpfr_prec_round(c, pc, GMP_RNDN); /* exact */
}


/* Computes c = a - b exactly, adapting the precision of c */
static inline void mpfr_sub_exact(mpfr_t c, mpfr_t a, mpfr_t b) {
  mp_exp_t Ea, Eb;
  mp_prec_t pa, pb, pc;
  

  if ((!mpfr_number_p(a)) ||
      (!mpfr_number_p(b))) {
    mpfr_sub(c, a, b, GMP_RNDN); /* Producing NaN or Inf */
    return;
  }

  if (mpfr_zero_p(a)) {
    mpfr_set_prec(c, mpfr_get_prec(b));
    mpfr_neg(c, b, GMP_RNDN); /* exact */
    return;
  }

  if (mpfr_zero_p(b)) {
    mpfr_set_prec(c, mpfr_get_prec(a));
    mpfr_set(c, a, GMP_RNDN); /* exact */
    return;
  }

  /* Here, both inputs are non-zero numbers */
  
  Ea = mpfr_get_exp(a);
  Eb = mpfr_get_exp(b);
  pa = mpfr_get_prec(a);
  pb = mpfr_get_prec(b);

  /* The ulps are the following

     ulp_pa(a) = 2^(Ea - pa)
     ulp_pb(b) = 2^(Eb - pb)

     So the ulp of the sum must be

     ulp_pc(c) = 2^min(Ea - pa, Eb - pb)
     
     The exponent of the result c is bounded
     by

     Ec <= max(Ea, Eb) + 1 

     The + 1 comes from the possible outgoing carry.

     So we need 

     pc = max(Ea, Eb) - min(Ea - pa, Eb - pb) + 1 bits
     
     for c.

  */
  pc = (Ea > Eb? Ea : Eb) - (Ea - pa < Eb - pb ? Ea - pa : Eb - pb) + 1; /* Exponent overflow possible */
  mpfr_set_prec(c, pc);

  /* Perform the subtraction */
  mpfr_sub(c, a, b, GMP_RNDN); /* exact */

  /* It could be that we actually don't need all the precision of c to
     keep the result. So we reduce that precision to the minimum.
  */
  pc = mpfr_min_prec(c);
  if (pc < 12) pc = 12;
  
  /* Adjust the precision of c again. */
  mpfr_prec_round(c, pc, GMP_RNDN); /* exact */
}

/* Computes c = a - b exactly, adapting the precision of c */
static inline void mpfr_sub_exact_int(mpfr_t c, mpfr_t a, int b) {
  mp_exp_t Ea, Eb;
  mp_prec_t pa, pb, pc;
  

  if (!mpfr_number_p(a)) {
    mpfr_sub_si(c, a, b, GMP_RNDN); /* Producing NaN or Inf */
    return;
  }

  if (mpfr_zero_p(a)) {
    mpfr_set_prec(c, 8 * sizeof(b) + 5);
    mpfr_set_si(c, b, GMP_RNDN); /* exact */
    mpfr_neg(c, c, GMP_RNDN); /* exact */
    return;
  }

  if (b == 0) {
    mpfr_set_prec(c, mpfr_get_prec(a));
    mpfr_set(c, a, GMP_RNDN); /* exact */
    return;
  }

  /* Here, both inputs are non-zero numbers */
  
  Ea = mpfr_get_exp(a);
  Eb = 8 * sizeof(int);
  pa = mpfr_get_prec(a);
  pb = 8 * sizeof(int);

  /* The ulps are the following

     ulp_pa(a) = 2^(Ea - pa)
     ulp_pb(b) = 2^(Eb - pb)

     So the ulp of the sum must be

     ulp_pc(c) = 2^min(Ea - pa, Eb - pb)
     
     The exponent of the result c is bounded
     by

     Ec <= max(Ea, Eb) + 1 

     The + 1 comes from the possible outgoing carry.

     So we need 

     pc = max(Ea, Eb) - min(Ea - pa, Eb - pb) + 1 bits
     
     for c.

  */
  pc = (Ea > Eb? Ea : Eb) - (Ea - pa < Eb - pb ? Ea - pa : Eb - pb) + 1; /* Exponent overflow possible */
  mpfr_set_prec(c, pc);

  /* Perform the subtraction */
  mpfr_sub_si(c, a, b, GMP_RNDN); /* exact */

  /* It could be that we actually don't need all the precision of c to
     keep the result. So we reduce that precision to the minimum.
  */
  pc = mpfr_min_prec(c);
  if (pc < 12) pc = 12;
  
  /* Adjust the precision of c again. */
  mpfr_prec_round(c, pc, GMP_RNDN); /* exact */
}

static inline void mpfr_int_sub_exact(mpfr_t c, int a, mpfr_t b) {
  mpfr_sub_exact_int(c, b, a);
  mpfr_neg(c, c, GMP_RNDN); /* exact */
}

/* Computes c = a * b exactly, adapting the precision of c */
static inline void mpfr_mul_exact(mpfr_t c, mpfr_t a, mpfr_t b) {
  mp_prec_t pa, pb, pc;
  
  if ((!mpfr_number_p(a)) ||
      (!mpfr_number_p(b))) {
    mpfr_mul(c, a, b, GMP_RNDN); /* Producing NaN or Inf */
    return;
  }

  if (mpfr_zero_p(a)) {
    mpfr_set_prec(c, mpfr_get_prec(a));
    mpfr_set(c, a, GMP_RNDN); /* exact */
    return;
  }

  if (mpfr_zero_p(b)) {
    mpfr_set_prec(c, mpfr_get_prec(b));
    mpfr_set(c, b, GMP_RNDN); /* exact */
    return;
  }

  if (mpfr_cmp_si(a,1) == 0) {
    mpfr_set_prec(c, mpfr_get_prec(b));
    mpfr_set(c, b, GMP_RNDN); /* exact */
    return;
  }
  
  if (mpfr_cmp_si(b,1) == 0) {
    mpfr_set_prec(c, mpfr_get_prec(a));
    mpfr_set(c, a, GMP_RNDN); /* exact */
    return;
  }

  /* Here, both inputs are non-zero numbers */
  pa = mpfr_get_prec(a);
  pb = mpfr_get_prec(b);

  /* Need pc = pa + pb bits of precision for an exact product. */
  pc = pa + pb; /* Precision overflow possible */
  mpfr_set_prec(c, pc);

  /* Perform the multiplication */
  mpfr_mul(c, a, b, GMP_RNDN); /* exact */

  /* It could be that we actually don't need all the precision of c to
     keep the result. So we reduce that precision to the minimum.
  */
  pc = mpfr_min_prec(c);
  if (pc < 12) pc = 12;
  
  /* Adjust the precision of c again. */
  mpfr_prec_round(c, pc, GMP_RNDN); /* exact */
}

/* Computes c = a * b exactly, adapting the precision of c */
static inline void mpfr_mul_exact_int(mpfr_t c, mpfr_t a, int b) {
  mp_prec_t pa, pb, pc;
  
  if (!mpfr_number_p(a)) {
    mpfr_mul_si(c, a, b, GMP_RNDN); /* Producing NaN or Inf */
    return;
  }

  if (mpfr_zero_p(a)) {
    mpfr_set_prec(c, mpfr_get_prec(a));
    mpfr_set(c, a, GMP_RNDN); /* exact */
    return;
  }

  if (b == 0) {
    mpfr_set_prec(c, 8 * sizeof(int) + 5);
    mpfr_set_si(c, b, GMP_RNDN); /* exact */
    return;
  }

  if (mpfr_cmp_si(a,1) == 0) {
    mpfr_set_prec(c, 8 * sizeof(int) + 5);
    mpfr_set_si(c, b, GMP_RNDN); /* exact */
    return;
  }
  
  if (b == 1) {
    mpfr_set_prec(c, mpfr_get_prec(a));
    mpfr_set(c, a, GMP_RNDN); /* exact */
    return;
  }

  /* Here, both inputs are non-zero numbers */
  pa = mpfr_get_prec(a);
  pb = 8 * sizeof(int) + 5;

  /* Need pc = pa + pb bits of precision for an exact product. */
  pc = pa + pb; /* Precision overflow possible */
  mpfr_set_prec(c, pc);

  /* Perform the multiplication */
  mpfr_mul_si(c, a, b, GMP_RNDN); /* exact */

  /* It could be that we actually don't need all the precision of c to
     keep the result. So we reduce that precision to the minimum.
  */
  pc = mpfr_min_prec(c);
  if (pc < 12) pc = 12;
  
  /* Adjust the precision of c again. */
  mpfr_prec_round(c, pc, GMP_RNDN); /* exact */
}

static inline int try_mpfr_pow_exact(mpfr_t c, mpfr_t a, unsigned long int b) {
  mp_prec_t pa, pc, paR;
  unsigned int paI, pcI, pcII; 
  unsigned int bI;
  unsigned long int bb;

  /* Special input NaN or Inf */
  if (!mpfr_number_p(a)) {
    /* Producing NaN, Inf or 1 */
    mpfr_set_prec(c,12);
    mpfr_pow_ui(c,a,b,GMP_RNDN);
    return 1;
  }

  /* Here, a is a real number. 

     Handle the easy cases b = 0 and b = 1.

  */
  if (b == 0) {
    mpfr_set_prec(c, 12);
    mpfr_set_si(c, 1, GMP_RNDN); /* exact */
    return 1;
  }

  if (b == 1) {
    mpfr_set_prec(c, mpfr_get_prec(a));
    mpfr_set(c, a, GMP_RNDN); /* exact */
    return 1;
  }

  /* It's so easy to handle b = 2, so let's do that */
  if (b == 2) {
    mpfr_set_prec(c, 2 * mpfr_get_prec(a));
    mpfr_mul(c, a, a, GMP_RNDN); /* exact */
    pc = mpfr_min_prec(c);
    if (pc < 12) pc = 12;
    mpfr_prec_round(c, pc, GMP_RNDN); /* exact */    
    return 1;    
  }

  /* Let's see if we can represent b on an unsigned int bI.
     If we can't, we refuse to do the job.
  */
  bI = b;
  bb = bI;
  if (bb != b) return 0;

  /* Handle the general case */
  pa = mpfr_get_prec(a);
  if (pa <= 1) return 0;

  /* pa >= 2 > 0 */
  paI = pa;
  paR = paI;
  if (paR != pa) return 0;

  /* Now we know that a is on paI bits, where paI is an unsigned int 

     At most, we need paI * bI bits of precision to represent a^bI.

     We try if we can represent paI * b on an unsigned int. If we
     can't, we refuse to do the job.

  */
  if (!tryExactUnsignedIntMultiplication(&pcI, paI, bI)) return 0;

  /* Now we know that we need pcI bits. We check if we can write pcI
     on an mp_prec_t. We first make sure that pcI is at least 12.
  */
  if (pcI < 12u) pcI = 12u;
  pc = pcI;
  if (pc < 12) pc = 12;
  pcII = pc;
  if (pcII != pcI) return 0;

  /* Now we know that we need pc bits for a^b */
  mpfr_set_prec(c, pc);
  mpfr_pow_ui(c,a,b,GMP_RNDN); /* exact */
  
  /* It could be that we actually don't need all the precision of c to
     keep the result. So we reduce that precision to the minimum.
  */
  pc = mpfr_min_prec(c);
  if (pc < 12) pc = 12;
  
  /* Adjust the precision of c again. */
  mpfr_prec_round(c, pc, GMP_RNDN); /* exact */
    
  return 1;
}

/* Part for constants */

static inline constant_t __polynomialGetIthCoefficientAsConstantIntIndex(polynomial_t, int);
static inline constant_t __polynomialGetIthCoefficientAsConstant(polynomial_t, mpz_t);

static inline void constantFree(constant_t);
static inline constant_t constantFromCopy(constant_t);
static inline void constantEvalMpfi(sollya_mpfi_t, constant_t); 
static inline constant_t constantAdd(constant_t, constant_t);
static inline constant_t constantSub(constant_t, constant_t);
static inline constant_t constantMul(constant_t, constant_t);
static inline constant_t constantDiv(constant_t, constant_t);
static inline constant_t constantPow(constant_t, constant_t);
static inline constant_t constantNeg(constant_t);

static inline int constantIsNonNegativeInteger(constant_t, int);
static inline int tryConstantToUnsignedInt(unsigned int *, constant_t);
static inline int tryConstantToMpz(mpz_t, constant_t);
static inline int constantIsGreater(constant_t, constant_t, int);

static inline void constantInitializeCaches() {
  int i;
  
  if (__constant_cache_initialized) return;
  for (i=0;i<CONSTANT_INTEGER_CACHE_SIZE;i++) {
    __constant_integer_cache[i].constant = NULL;
    __constant_integer_cache[i].initialized = 0;
  }
  __constant_malloc_cache_index = 0;
  __constant_cache_initialized = 1;
}

static inline void constantFreeCaches() {
  int i;

  if (!__constant_cache_initialized) return;
  for (i=0;i<CONSTANT_INTEGER_CACHE_SIZE;i++) {
    if (__constant_integer_cache[i].initialized) {
      constantFree(__constant_integer_cache[i].constant);
      __constant_integer_cache[i].initialized = 0;
    }
  }
  for (i=0;i<__constant_malloc_cache_index;i++) {
    safeFree(__constant_malloc_cache[i]);
  }
  __constant_malloc_cache_index = 0;
  __constant_cache_initialized = 0;
}

static inline constant_t __constantAllocate() {
  constantInitializeCaches();
  if (__constant_cache_initialized) {
    if ((0 < __constant_malloc_cache_index) && (__constant_malloc_cache_index <= CONSTANT_MALLOC_CACHE_SIZE)) {
      __constant_malloc_cache_index--;
      return __constant_malloc_cache[__constant_malloc_cache_index];
    }
  }
  return (constant_t) safeMalloc(sizeof(struct __constant_struct_t));
}

static inline void __constantFreeMem(constant_t c) {
  if (__constant_cache_initialized) {
    if ((0 <= __constant_malloc_cache_index) && (__constant_malloc_cache_index < CONSTANT_MALLOC_CACHE_SIZE)) {
      __constant_malloc_cache[__constant_malloc_cache_index] = c;
      __constant_malloc_cache_index++;
      return;
    }
  }
  safeFree(c);
}

static inline constant_t __constantAllocatePostTreatment(constant_t c) {
  return c;
}

static inline constant_t constantFromInt(int c) {
  constant_t res;

  if ((CONSTANT_INTEGER_CACHE_MIN <= c) && (c <= CONSTANT_INTEGER_CACHE_MAX)) {
    constantInitializeCaches();
    if (__constant_integer_cache[c - (CONSTANT_INTEGER_CACHE_MIN)].initialized) {
      return constantFromCopy(__constant_integer_cache[c - (CONSTANT_INTEGER_CACHE_MIN)].constant);
    }
  }
  
  res = __constantAllocate();
  res->refCount = 1;
  res->isZero.cached = 0;
  res->isOne.cached = 0;
  res->isNonNegativeInteger.cached = 0;
  res->isPositive.cached = 0;
  res->isDyadic.cached = 0;
  res->isRational.cached = 0;
  res->type = INTEGER;
  res->value.integer = c;
  res->hash.hasHash = 0;
  res = __constantAllocatePostTreatment(res);

  if ((CONSTANT_INTEGER_CACHE_MIN <= c) && (c <= CONSTANT_INTEGER_CACHE_MAX)) {
    constantInitializeCaches();
    if (!(__constant_integer_cache[c - (CONSTANT_INTEGER_CACHE_MIN)].initialized)) {
      __constant_integer_cache[c - (CONSTANT_INTEGER_CACHE_MIN)].constant = constantFromCopy(res);
      __constant_integer_cache[c - (CONSTANT_INTEGER_CACHE_MIN)].initialized = 1;
    }
  }
  
  return res;
}

static inline mp_prec_t mpfr_get_needed_prec(mpfr_t c) {
  mp_prec_t p;

  p = mpfr_min_prec(c);
  if (p < ((mp_prec_t) 12)) p = ((mp_prec_t) 12);

  return p;
}

static inline constant_t constantFromMpfr(mpfr_t c) {
  constant_t res;
  int intval;

  if (mpfr_is_machine_integer(&intval, c)) {
    return constantFromInt(intval);
  }
  
  res = __constantAllocate();
  res->refCount = 1;
  res->isZero.cached = 0;
  res->isOne.cached = 0;
  res->isNonNegativeInteger.cached = 0;
  res->isPositive.cached = 0;
  res->isDyadic.cached = 0;
  res->isRational.cached = 0;
  res->hash.hasHash = 0;
  res->type = MPFR;
  mpfr_init2(res->value.mpfr, mpfr_get_needed_prec(c));
  mpfr_set(res->value.mpfr, c, GMP_RNDN); /* exact, enough precision */
  
  return __constantAllocatePostTreatment(res);
}

static inline constant_t constantFromUnsignedInt(unsigned int c) {
  constant_t res;
  int cI;
  unsigned int cc;
  mpfr_t t;

  cI = c;
  cc = cI;
  if ((cI < 0) || (cc != c)) {
    mpfr_init2(t, 8 * sizeof(unsigned int) + 5);
    mpfr_set_ui(t, c, GMP_RNDN); /* exact */
    res = constantFromMpfr(t);
    mpfr_clear(t);
    return res;
  } 

  return constantFromInt(cI);
}

/* The argument c must be a canonicalized mpq_t number, i.e. zero 
   needs to be represented as 0/1, the numerator and the denominator 
   of c must not have any common factor and the denominator must be 
   positive.
*/
static inline constant_t constantFromScaledMpq(mp_exp_t e, mpq_t c) {
  mpz_t num, den;
  constant_t res;
  mp_bitcnt_t dyadNum, dyadDen;
  mp_exp_t expo;
  mp_prec_t p;
  mpfr_t mpfrval;

  if (mpq_sgn(c) == 0) {
    return constantFromInt(0);
  }

  mpz_init(num);
  mpz_init(den);
  mpq_get_num(num, c);
  mpq_get_den(den, c);
  dyadNum = mpz_scan1(num, 0);
  dyadDen = mpz_scan1(den, 0);
  mpz_tdiv_q_2exp(num, num, dyadNum);
  mpz_tdiv_q_2exp(den, den, dyadDen);
  expo = e + dyadNum - dyadDen; /* Exponent overflow possible */

  if (mpz_cmp_si(den, 1) == 0) {
    /* The denominator is one, so we can actually 
       use a MPFR (or integer) representation.
    */
    p = mpz_sizeinbase(num, 2);
    if (p < 12) p = 12;
    mpfr_init2(mpfrval, p);
    mpfr_set_z_2exp(mpfrval, num, expo, GMP_RNDN); /* exact as enough precision */
    if (mpfr_number_p(mpfrval)) {
      res = constantFromMpfr(mpfrval);
      mpfr_clear(mpfrval);
      mpz_clear(den);
      mpz_clear(num);
      return res;
    } else {
      mpfr_clear(mpfrval);
    }
  }

  /* Here, we are sure we must use a scaled MPQ representation */
  res = __constantAllocate();
  res->refCount = 1;
  res->isZero.cached = 0;
  res->isOne.cached = 0;
  res->isNonNegativeInteger.cached = 0;
  res->isPositive.cached = 0;
  res->isDyadic.cached = 0;
  res->isRational.cached = 0;
  res->type = SCALEDMPQ;
  res->value.scaledMpq.expo = expo;
  res->hash.hasHash = 0;
  mpq_init(res->value.scaledMpq.significand);
  mpq_set_num(res->value.scaledMpq.significand, num);
  mpq_set_den(res->value.scaledMpq.significand, den);
  mpz_clear(den);
  mpz_clear(num);
  
  return __constantAllocatePostTreatment(res);
}

static inline constant_t constantFromMpq(mpq_t c) {
  return constantFromScaledMpq((mp_exp_t) 0, c);
}

static inline constant_t constantFromMpz(mpz_t c) {
  mpq_t q;
  constant_t res;
  
  mpq_init(q);
  mpq_set_z(q, c);
  mpq_canonicalize(q);
  res = constantFromMpq(q);
  mpq_clear(q);
  
  return res;
}

static inline constant_t constantFromBinomialUnsignedInt(unsigned int n, unsigned int i) {
  mpz_t bin;
  constant_t res;

  mpz_init(bin);
  mpz_bin_uiui(bin, n, i);
  res = constantFromMpz(bin);
  mpz_clear(bin);

  return res;
}

static inline int __constantFromBinomialConstantAndUnsignedIntUnsafe(constant_t *res, constant_t n, unsigned int k) {
  mpz_t nz, bin;

  mpz_init(nz);
  if (!tryConstantToMpz(nz, n)) {
    mpz_clear(nz);
    return 0;
  }

  mpz_init(bin);
  mpz_bin_ui(bin, nz, k);
  *res = constantFromMpz(bin);
  mpz_clear(bin);
  mpz_clear(nz);
  return 1;
}

static inline int __constantFromBinomialUnsafe(constant_t *bin, constant_t n, constant_t k) {
  unsigned int nu, ku, ju;
  constant_t j;
  
  /* Check whether both n and k are representable as unsigned machine
     integers. 
  */
  
  if (tryConstantToUnsignedInt(&nu, n) &&
      tryConstantToUnsignedInt(&ku, k)) {
    /* Here we can use the special function on machine integers. */
    *bin = constantFromBinomialUnsignedInt(nu, ku);
    return 1;
  }

  /* Here, at least one of n or k is not representable as an unsigned
     machine integer. 

     Continue by checking whether at least k is representable as an
     unsigned machine integer.

  */
  if (tryConstantToUnsignedInt(&ku, k)) {
    /* Use the special function for bin(n, k) with k unsigned int */
    return __constantFromBinomialConstantAndUnsignedIntUnsafe(bin, n, ku);
  }

  /* Compute j = n - k and check if j (and n) are representable as
     unsigned machine integers. 
  */
  j = constantSub(n, k);
  if (!constantIsNonNegativeInteger(j, 0)) {
    constantFree(j);
    return 0;
  }

  if (!tryConstantToUnsignedInt(&ju, j)) {
    constantFree(j);
    return 0;    
  }
  constantFree(j);
  
  if (tryConstantToUnsignedInt(&nu, n)) {
    *bin = constantFromBinomialUnsignedInt(nu, ju);
    return 1;
  }

  return __constantFromBinomialConstantAndUnsignedIntUnsafe(bin, n, ju);
}

static inline int constantFromBinomial(constant_t *bin, constant_t n, constant_t k) {

  /* Handle stupid input */
  if (n == NULL) return 0;
  if (k == NULL) return 0;

  /* Check if n and k are both non-negative integers and if k is no
     greater than n. 
  */
  if (constantIsNonNegativeInteger(n, 0) &&
      constantIsNonNegativeInteger(k, 0) &&
      (!constantIsGreater(k, n, 1))) {
    /* Both n and k are non-negative integers.

       Use an auxilliary function and return its return value.
    */
    return __constantFromBinomialUnsafe(bin, n, k);
  }

  /* Here one of n or k is negative or no integer. Signal failure. */  
  return 0;
}

static inline constant_t constantFromCopy(constant_t c) {
  if (c == NULL) return NULL;
  c->refCount++;
  return c;
}

static inline constant_t constantFromExpression(node *c) {
  node *simplified, *copy;
  mpq_t rational;
  constant_t res;
  
  if (c == NULL) return NULL;
  if (!isConstant(c)) return NULL;
  if ((c->nodeType == MEMREF) &&
      (c->polynomialRepresentation != NULL)) {
    /* The expression is a polynomial but is constant. The constant of
       degree 0 does the job 
    */
    return __polynomialGetIthCoefficientAsConstantIntIndex(c->polynomialRepresentation, 0);
  }
  if (accessThruMemRef(c)->nodeType == CONSTANT) {
    return constantFromMpfr(*(accessThruMemRef(c)->value));
  }
  mpq_init(rational);
  if (tryEvaluateConstantTermToMpq(rational, c)) {
    mpq_canonicalize(rational);
    res = constantFromMpq(rational);
    mpq_clear(rational);
    return res;
  }
  mpq_clear(rational);
  simplified = simplifyRationalErrorfree(c);
  if (accessThruMemRef(simplified)->nodeType == CONSTANT) {
    res = constantFromMpfr(*(accessThruMemRef(simplified)->value));
    freeThing(simplified);
    return res;
  }
  if (simplified == c) {
    copy = addMemRef(copyThing(accessThruMemRef(simplified)));
    tryCopyTreeAnnotations(copy, simplified);
    freeThing(simplified);
    simplified = copy;
  }
  
  res = __constantAllocate();
  res->refCount = 1;
  res->isZero.cached = 0;
  res->isOne.cached = 0;
  res->isNonNegativeInteger.cached = 0;
  res->isPositive.cached = 0;
  res->isDyadic.cached = 0;
  res->isRational.cached = 0;
  res->type = EXPRESSION;
  res->value.expr = simplified;
  res->hash.hasHash = 0;
  
  return __constantAllocatePostTreatment(res);
}


static inline void constantFree(constant_t c) {
  if (c == NULL) return;
  c->refCount--;
  if (c->refCount > 0u) return;
  switch (c->type) {
  case INTEGER:
    break;
  case EXPRESSION:
    freeThing(c->value.expr);
    break;
  case MPFR:
    mpfr_clear(c->value.mpfr);
    break;
  case SCALEDMPQ:
    mpq_clear(c->value.scaledMpq.significand);
    break;
  }
  __constantFreeMem(c);
}

static inline int constantIsZero(constant_t a, int defVal) {
  int s;

  if (a == NULL) return defVal;
  if (a->isZero.cached) return a->isZero.res;
  if (a->isOne.cached && a->isOne.res) return 0;
  if (a->isPositive.cached && a->isPositive.res) return 0;
  switch (a->type) {
  case INTEGER:
    a->isZero.cached = 1;
    a->isZero.res = (a->value.integer == 0);
    return a->isZero.res;
    break;
  case EXPRESSION:
    if (accessThruMemRef(a->value.expr)->nodeType == CONSTANT) {
      if (!mpfr_number_p(*(accessThruMemRef(a->value.expr)->value))) return defVal;
      a->isZero.cached = 1;
      a->isZero.res = mpfr_zero_p(*(accessThruMemRef(a->value.expr)->value));
      return a->isZero.res;
    }
    if (evaluateSignFast(&s, a->value.expr)) {
      if (s == 0) {
	a->isZero.cached = 1;
	a->isZero.res = 1;	
	return a->isZero.res;
      } else {
	a->isZero.cached = 1;
	a->isZero.res = 0;	
	return a->isZero.res;
      }
    }
    return defVal;
    break;
  case MPFR:
    if (!mpfr_number_p(a->value.mpfr)) return defVal;
    a->isZero.cached = 1;
    a->isZero.res = mpfr_zero_p(a->value.mpfr);
    return a->isZero.res;
    break;
  case SCALEDMPQ:
    a->isZero.cached = 1;
    a->isZero.res = (mpq_sgn(a->value.scaledMpq.significand) == 0);
    return a->isZero.res;
    break;
  }
  return defVal;
}

static inline int constantIsOne(constant_t a, int defVal) {
  int s;

  if (a == NULL) return defVal;
  if (a->isOne.cached) return a->isOne.res;
  if (a->isZero.cached && a->isZero.res) return 0;
  switch (a->type) {
  case INTEGER:
    a->isOne.cached = 1;
    a->isOne.res = (a->value.integer == 1);
    return a->isOne.res;
    break;
  case EXPRESSION:
    if (accessThruMemRef(a->value.expr)->nodeType == CONSTANT) {
      if (!mpfr_number_p(*(accessThruMemRef(a->value.expr)->value))) return defVal;
      a->isOne.cached = 1;
      a->isOne.res = (mpfr_cmp_si(*(accessThruMemRef(a->value.expr)->value),1) == 0);
      return a->isOne.res;
    }
    if (evaluateSignFast(&s, a->value.expr)) {
      if (s <= 0) {
	a->isOne.cached = 1;
	a->isOne.res = 0;
	return a->isOne.res;
      }
    }
    return defVal;
    break;
  case MPFR:
    if (!mpfr_number_p(a->value.mpfr)) return defVal;
    a->isOne.cached = 1;
    a->isOne.res = (mpfr_cmp_si(a->value.mpfr,1) == 0);
    return a->isOne.res;
    break;
  case SCALEDMPQ:
    /* We have to compare 2^E * n/d == 1. This means checking 2^E * n == d */
    if (mpq_sgn(a->value.scaledMpq.significand) <= 0) {
      a->isOne.cached = 1;
      a->isOne.res = 0;
      return a->isOne.res;
    }
    if (a->value.scaledMpq.expo == ((mp_exp_t) 0)) {
      /* Here, E = 0, so we have to check whether n == d */
      a->isOne.cached = 1;
      if (mpz_cmp(mpq_numref(a->value.scaledMpq.significand), 
		  mpq_denref(a->value.scaledMpq.significand)) == 0) {
	a->isOne.res = 1;
      } else {
	a->isOne.res = 0;
      }
      return a->isOne.res;
    }
    /* Here E != 0 */
    a->value.scaledMpq.expo += mpq_remove_powers_of_two(a->value.scaledMpq.significand); /* Exponent overflow possible */
    if (a->value.scaledMpq.expo != ((mp_exp_t) 0)) {
      a->isOne.cached = 1;
      a->isOne.res = 0;
      return a->isOne.res;
    }
    if (mpz_cmp(mpq_numref(a->value.scaledMpq.significand), 
		mpq_denref(a->value.scaledMpq.significand)) == 0) {
      a->isOne.cached = 1;
    } else {
      a->isOne.cached = 0;
    }
    return a->isOne.res;
    break;
  }
  return defVal;
}

static inline int constantIsNonNegativeInteger(constant_t a, int defVal) {
  int s;
  node *t, *d;

  if (a == NULL) return defVal;
  if (a->isNonNegativeInteger.cached) return a->isNonNegativeInteger.res;
  if (a->isZero.cached && a->isZero.res) return 1;
  if (a->isOne.cached && a->isOne.res) return 1;
  switch (a->type) {
  case INTEGER:
    a->isNonNegativeInteger.cached = 1;
    a->isNonNegativeInteger.res = (a->value.integer >= 0);
    return a->isNonNegativeInteger.res;
    break;
  case EXPRESSION:
    if (accessThruMemRef(a->value.expr)->nodeType == CONSTANT) {
      if (!mpfr_number_p(*(accessThruMemRef(a->value.expr)->value))) return defVal;
      a->isNonNegativeInteger.cached = 1;
      a->isNonNegativeInteger.res = ((mpfr_sgn(*(accessThruMemRef(a->value.expr)->value)) >= 0) && 
				     mpfr_integer_p(*(accessThruMemRef(a->value.expr)->value)));
      return a->isNonNegativeInteger.res;
    }
    if (!evaluateSignFast(&s, a->value.expr)) return defVal;
    if (s < 0) {
      a->isNonNegativeInteger.cached = 1;
      a->isNonNegativeInteger.res = 0;
      return a->isNonNegativeInteger.res;
    }
    t = addMemRef(makeSub(copyThing(a->value.expr),
			  makeNearestInt(copyThing(a->value.expr))));
    d = simplifyRationalErrorfree(t);
    freeThing(t);
    if (!evaluateSignFast(&s, d)) {
      freeThing(d);
      return defVal;
    }
    freeThing(d);
    a->isNonNegativeInteger.cached = 1;
    a->isNonNegativeInteger.res = (s == 0);
    return a->isNonNegativeInteger.res;
    break;
  case MPFR:
    if (!mpfr_number_p(a->value.mpfr)) return defVal;
    a->isNonNegativeInteger.cached = 1;
    a->isNonNegativeInteger.res = (mpfr_integer_p(a->value.mpfr) && 
				   (mpfr_sgn(a->value.mpfr) >= 0));
    return a->isNonNegativeInteger.res;
    break;
  case SCALEDMPQ:
    if (mpq_sgn(a->value.scaledMpq.significand) < 0) return 0;
    a->isNonNegativeInteger.cached = 1;
    a->isNonNegativeInteger.res = scaledMpqIsInteger(a->value.scaledMpq.expo, 
						     a->value.scaledMpq.significand);
    return a->isNonNegativeInteger.res;
    break;
  }

  return defVal;
}

static inline int constantIsDyadic(constant_t a, int defVal) {
  if (a == NULL) return defVal;
  if (a->isDyadic.cached) return a->isDyadic.res;
  switch (a->type) {
  case INTEGER:
    a->isDyadic.cached = 1;
    a->isDyadic.res = 1;
    return a->isDyadic.res;
    break;
  case EXPRESSION:
    if (accessThruMemRef(a->value.expr)->nodeType == CONSTANT) {
      if (!mpfr_number_p(*(accessThruMemRef(a->value.expr)->value))) return defVal;
      a->isDyadic.cached = 1;
      a->isDyadic.res = 1;
      return a->isDyadic.res;
    }
    return defVal;
    break;
  case MPFR:
    if (!mpfr_number_p(a->value.mpfr)) return defVal;
    a->isDyadic.cached = 1;
    a->isDyadic.res = 1;
    return a->isDyadic.res;
    break;
  case SCALEDMPQ:
    a->isDyadic.cached = 1;
    a->isDyadic.res = scaledMpqIsDyadic(a->value.scaledMpq.expo, 
					a->value.scaledMpq.significand);
    return a->isDyadic.res;
    break;
  }
  return defVal;
}

static inline int constantIsRational(constant_t a, int defVal) {
  if (a == NULL) return defVal;
  if (a->isRational.cached) return a->isRational.res;
  switch (a->type) {
  case INTEGER:
    a->isRational.cached = 1;
    a->isRational.res = 1;
    return a->isRational.res;
    break;
  case EXPRESSION:
    if (accessThruMemRef(a->value.expr)->nodeType == CONSTANT) {
      if (!mpfr_number_p(*(accessThruMemRef(a->value.expr)->value))) return defVal;
      a->isRational.cached = 1;
      a->isRational.res = 1;
      return a->isRational.res;
    }
    return defVal;
    break;
  case MPFR:
    if (!mpfr_number_p(a->value.mpfr)) return defVal;
    a->isRational.cached = 1;
    a->isRational.res = 1;
    return a->isRational.res;
    break;
  case SCALEDMPQ:
    a->isRational.cached = 1;
    a->isRational.res = 1;
    return a->isRational.res;
    break;
  }
  return defVal;
}

static inline int constantHoldsOnPrecBits(constant_t a, mp_prec_t prec, int defVal) {
  sollya_mpfi_t t;
  int res;

  if (a == NULL) return defVal;
  if (prec < 2) return defVal;

  sollya_mpfi_init2(t, prec + 15);
  constantEvalMpfi(t, a);
  if (sollya_mpfi_has_nan(t) ||
      sollya_mpfi_has_infinity(t)) {
    res = defVal;
  } else {
    res = sollya_mpfi_is_point_and_real(t);
  }
  sollya_mpfi_clear(t);

  return res;
}

static inline int constantIsEqual(constant_t a, constant_t b, int defVal) {
  constant_t d;
  int res, sa, sb;
  mpq_t t;
  mp_exp_t G;

  /* Trivial answers */
  if (a == NULL) return defVal;
  if (b == NULL) return defVal;
  if (a == b) return 1;

  /* Using a hash here is actually a bad idea as we are talking about 
     mathematical identity not the structural one.
  */
  
  if (a->type != b->type) {
    /* If the constants are not of the same type, compute their
       difference and compare to zero. 
    */
    d = constantSub(a, b);
    res = constantIsZero(d, defVal);
    constantFree(d);
    return res;
  }

  /* Here a and b are of the same type */
  switch (a->type) {
  case INTEGER:
    return (a->value.integer == b->value.integer);
    break;
  case EXPRESSION:
    if ((accessThruMemRef(a->value.expr)->nodeType == CONSTANT) && 
	(accessThruMemRef(b->value.expr)->nodeType == CONSTANT)) {
      if (!mpfr_number_p(*(accessThruMemRef(a->value.expr)->value))) return defVal;
      if (!mpfr_number_p(*(accessThruMemRef(b->value.expr)->value))) return defVal;
      return mpfr_equal_p(*(accessThruMemRef(a->value.expr)->value),
			  *(accessThruMemRef(b->value.expr)->value));
    }
    if (evaluateSignFast(&sa, a->value.expr) &&
	evaluateSignFast(&sb, b->value.expr)) {
      if ((sa == 0) && (sb == 0)) {
	/* If both expressions are zero, they are equal */
	return 1;
      }
      if (sa * sb <= 0) {
	/* If the expressions are of different sign, they are not equal */
	return 0;
      }
    }
    return defVal;
    break;
  case MPFR:
    if (!mpfr_number_p(a->value.mpfr)) return defVal;
    if (!mpfr_number_p(b->value.mpfr)) return defVal;
    return mpfr_equal_p(a->value.mpfr,b->value.mpfr);
    break;
  case SCALEDMPQ:
    sa = mpq_sgn(a->value.scaledMpq.significand);
    sb = mpq_sgn(b->value.scaledMpq.significand);
    /* Reason upon the sign first */
    if ((sa == 0) && (sb == 0)) return 1;
    if (sa * sb <= 0) return 0;
    /* Here a and b are of the same sign and none is zero 

       If their exponents are equal and their significands are equal,
       they are equal.

    */
    if ((a->value.scaledMpq.expo == b->value.scaledMpq.expo) && 
	(mpq_cmp(a->value.scaledMpq.significand, 
		 b->value.scaledMpq.significand) == 0)) 
      return 1;
    /* Here the exponents are unequal or the significands are unequal 

       In order to compare 2^E * n/d to 2^F * p/q, compare

       2^(E - F) to (p/q) / (n/d).

    */
    mpq_init(t);
    mpq_div(t, b->value.scaledMpq.significand, a->value.scaledMpq.significand);
    G = mpq_remove_powers_of_two(t);

    /* Here, we have to compare 2^(E - F) to 2^G * t, 

       where t = r/s = 2^-G * (p/q) / (n/d)

       is a fraction such that both r and s are odd. 

       This means if 

       * E - F != G, a and b are unequal
       * E - F == G:
            - if t is equal to one, a and b are equal
	    - otherwise, a and b are unequal.
    */
    if (a->value.scaledMpq.expo - b->value.scaledMpq.expo != G) { /* Exponent overflow possible */
      res = 0;
    } else {
      res = (mpq_cmp_si(t, 1, 1u) == 0);
    }
    mpq_clear(t);
    return res;
    break;
  }
 
  return defVal;
}

static inline int constantIsPositive(constant_t a, int defVal) {
  int s;

  if (a == NULL) return defVal;
  if (a->isPositive.cached) return a->isPositive.res;
  if (a->isOne.cached && a->isOne.res) return 1;
  if (a->isZero.cached && a->isZero.res) return 0;
  switch (a->type) {
  case INTEGER:
    a->isPositive.cached = 1;
    a->isPositive.res = (a->value.integer > 0);
    return a->isPositive.res;
    break;
  case EXPRESSION:
    if (accessThruMemRef(a->value.expr)->nodeType == CONSTANT) {
      if (!mpfr_number_p(*(accessThruMemRef(a->value.expr)->value))) return defVal;
      a->isPositive.cached = 1;
      a->isPositive.res = (mpfr_sgn(*(accessThruMemRef(a->value.expr)->value)) > 0);
      return a->isPositive.res;
    }
    if (evaluateSignFast(&s, a->value.expr)) {
      if (s > 0) {
	a->isPositive.cached = 1;
	a->isPositive.res = 1;	
	return a->isPositive.res;
      }
    }
    return defVal;
    break;
  case MPFR:
    if (!mpfr_number_p(a->value.mpfr)) return defVal;
    a->isPositive.cached = 1;
    a->isPositive.res = (mpfr_sgn(a->value.mpfr) > 0);
    return a->isPositive.res;
    break;
  case SCALEDMPQ:
    a->isPositive.cached = 1;
    a->isPositive.res = (mpq_sgn(a->value.scaledMpq.significand) > 0);
    return a->isPositive.res;
    break;
  }
  return defVal;
}

static inline int constantIsGreater(constant_t a, constant_t b, int defVal) {
  constant_t d;
  int res;
  node *dn;
  int sa, sb;

  /* Handle stupid inputs */
  if (a == NULL) return defVal;
  if (b == NULL) return defVal;
  if (a == b) return 0;

  /* If the inputs don't have the same type, compute the difference
     and compare with 0. 
  */
  if (a->type != b->type) {
    d = constantSub(a, b);
    res = constantIsPositive(d, defVal);
    constantFree(d);
    return res;
  }

  /* Here, the inputs have the same representation type */
  switch (a->type) {
  case INTEGER:
    return (a->value.integer > b->value.integer);
    break;
  case EXPRESSION:
    if ((accessThruMemRef(a->value.expr)->nodeType == CONSTANT) && 
	(accessThruMemRef(b->value.expr)->nodeType == CONSTANT)) {
      if (!mpfr_number_p(*(accessThruMemRef(a->value.expr)->value))) return defVal;
      if (!mpfr_number_p(*(accessThruMemRef(b->value.expr)->value))) return defVal;
      return (mpfr_cmp(*(accessThruMemRef(a->value.expr)->value),
		       *(accessThruMemRef(b->value.expr)->value)) > 0);
    }
    if (evaluateSignFast(&sa, a->value.expr) &&
	evaluateSignFast(&sb, b->value.expr)) {
      if (sa * sb <= 0) {
	/* The inputs have different sign or are both zero 
	   
	   If they are both zero, a is not greater than b.
	   If a is positive, b is negative and a is greater than b.
	   If a is negative, b is positive and a is not greater than b.

	*/
	if ((sa == 0) && (sb == 0)) return 0;
	if (sa > 0) return 1;
	return 0;
      }
    }
    /* Construct the difference and try to evaluate its sign */
    dn = addMemRef(makeSub(copyThing(a->value.expr),
			   copyThing(b->value.expr)));
    if (evaluateSignFast(&sa, dn)) {
      freeThing(dn);
      return (sa > 0);
    }
    freeThing(dn);
    return defVal;
    break;
  case MPFR:
    if (!mpfr_number_p(a->value.mpfr)) return defVal;
    if (!mpfr_number_p(b->value.mpfr)) return defVal;
    return (mpfr_cmp(a->value.mpfr,b->value.mpfr) > 0);
    break;
  case SCALEDMPQ:
    return scaledMpqIsGreaterThan(a->value.scaledMpq.expo,
				  a->value.scaledMpq.significand,
				  b->value.scaledMpq.expo,
				  b->value.scaledMpq.significand);
    break;
  }

  return defVal;
}

static inline int constantIsGreaterOrEqual(constant_t a, constant_t b, int defVal) {
  int resGreater, resEqual;

  resGreater = constantIsGreater(a, b, 42);
  if (resGreater == 42) {
    resEqual = constantIsEqual(a, b, 42);
    if (resEqual == 42) return defVal;
    if (resEqual) return 1;
    return defVal;
  }
  if (resGreater) return 1;
  return constantIsEqual(a, b, defVal);
}

static inline node *constantToExpression(constant_t a) {
  mpfr_t num, den;
  mp_prec_t p;
  node *res;
  mp_exp_t ED, EN, EE;

  if (a == NULL) return NULL;
  switch (a->type) {
  case INTEGER:
    return addMemRef(makeConstantInt(a->value.integer));
    break;
  case EXPRESSION:
    return copyThing(a->value.expr);
    break;
  case MPFR:
    return addMemRef(makeConstant(a->value.mpfr));
    break;
  case SCALEDMPQ:
    if (mpq_sgn(a->value.scaledMpq.significand) == 0) 
      return addMemRef(makeConstantInt(0));
    if (mpz_cmp_si(mpq_denref(a->value.scaledMpq.significand),1) == 0) {
      /* The denominator is 1, so we do not have to take it into account */
      p = mpz_sizeinbase(mpq_numref(a->value.scaledMpq.significand), 2);
      if (p < 12) p = 12;
      mpfr_init2(num, p);
      mpfr_set_z_2exp(num, 
		      mpq_numref(a->value.scaledMpq.significand), 
		      a->value.scaledMpq.expo, 
		      GMP_RNDN); /* exact as enough precision */
      if (mpfr_number_p(num)) {
	res = addMemRef(makeConstant(num));
      } else {
	mpfr_set_z(num, mpq_numref(a->value.scaledMpq.significand), GMP_RNDN); /* exact as enough precision */
	res = addMemRef(makeMul(makePow(makeConstantInt(2),
					makeConstantInt(a->value.scaledMpq.expo)), 
				makeConstant(num)));
      }
      mpfr_clear(num);
      return res;
    }
    /* Here, we must take both the numerator and the denominator into
       account.
    */
    p = mpz_sizeinbase(mpq_numref(a->value.scaledMpq.significand), 2);
    if (p < 12) p = 12;
    if (a->value.scaledMpq.expo >= 0) {
      EN = a->value.scaledMpq.expo;
      ED = 0;
    } else {
      ED = -a->value.scaledMpq.expo; /* Exponent overflow possible */
      EN = 0;
    }
    mpfr_init2(num, p);
    mpfr_set_z_2exp(num, 
		    mpq_numref(a->value.scaledMpq.significand), 
		    EN,
		    GMP_RNDN); /* exact as enough precision */
    if (mpfr_number_p(num)) {
      EN = 0;
    } else {
      mpfr_set_z(num, 
		 mpq_numref(a->value.scaledMpq.significand), 
		 GMP_RNDN); /* exact as enough precision */      
    }
    p = mpz_sizeinbase(mpq_denref(a->value.scaledMpq.significand), 2);
    if (p < 12) p = 12;
    mpfr_init2(den, p);
    mpfr_set_z_2exp(den, 
		    mpq_denref(a->value.scaledMpq.significand),
		    ED, 
		    GMP_RNDN); /* exact as enough precision */
    if (mpfr_number_p(den)) {
      ED = 0;
    } else {
      mpfr_set_z(den, 
		 mpq_denref(a->value.scaledMpq.significand), 
		 GMP_RNDN); /* exact as enough precision */      
    }
    EE = EN - ED; /* Exponent overflow possible */
    if (EE == 0) {
      res = addMemRef(makeDiv(makeConstant(num),makeConstant(den)));
    } else {
      res = addMemRef(makeMul(makePow(makeConstantInt(2),makeConstantInt(EE)),makeDiv(makeConstant(num),makeConstant(den))));
    }
    mpfr_clear(num);
    mpfr_clear(den);
    return res;
    break;
  }
  return NULL;
}

static inline int tryConstantToScaledMpq(mp_exp_t *E, mpq_t rop, constant_t a) {
  mpq_t t;
  mpz_t mant;
  mp_exp_t expo;

  if (a == NULL) return 0;
  switch (a->type) {
  case INTEGER:
    mpq_set_si(rop, a->value.integer, 1u);
    mpq_canonicalize(rop);
    *E = mpq_remove_powers_of_two(rop);
    return 1;
    break;
  case EXPRESSION:
    mpq_init(t);
    if (tryEvaluateConstantTermToMpq(t, a->value.expr)) {
      mpq_set(rop, t);
      mpq_canonicalize(rop);
      *E = mpq_remove_powers_of_two(rop);      
      mpq_clear(t);
      return 1;
    }
    mpq_clear(t);
    return 0;
    break;
  case MPFR:
    if (!mpfr_number_p(a->value.mpfr)) return 0;
    if (mpfr_zero_p(a->value.mpfr)) {
      *E = 0;
      mpq_set_si(rop, 0, 1u);
      mpq_canonicalize(rop);
      return 1;
    }
    mpz_init(mant);
    expo = mpfr_get_z_2exp(mant, a->value.mpfr);
    mpq_set_z(rop, mant);
    mpq_canonicalize(rop);
    *E = expo + mpq_remove_powers_of_two(rop);  /* Exponent overflow possible */
    mpz_clear(mant);
    return 1;
    break;
  case SCALEDMPQ:
    *E = a->value.scaledMpq.expo;
    mpq_set(rop, a->value.scaledMpq.significand);
    mpq_canonicalize(rop);
    return 1;
    break;
  }
  return 0;
}

static inline int tryConstantToMpq(mpq_t rop, constant_t a) {
  mpq_t q;
  mp_exp_t E;

  mpq_init(q);
  if (!tryConstantToScaledMpq(&E, q, a)) {
    mpq_clear(q);
    return 0;
  }

  if (E >= ((mp_exp_t) 0)) {
    mpq_mul_2exp(rop, q, ((mp_bitcnt_t) E));
  } else {
    mpq_div_2exp(rop, q, ((mp_bitcnt_t) (-E)));
  }
  
  return 1;
}

static inline int tryConstantToMpz(mpz_t r, constant_t a) {
  mpz_t num, den;
  mpq_t q;
  mp_exp_t EQ;
  
  if (a == NULL) return 0;
  switch (a->type) {
  case INTEGER:
    mpz_set_si(r, a->value.integer);
    return 1;
    break;
  case EXPRESSION:
    if ((accessThruMemRef(a->value.expr)->nodeType == CONSTANT) &&
	mpfr_number_p(*(accessThruMemRef(a->value.expr)->value)) && 
	mpfr_integer_p(*(accessThruMemRef(a->value.expr)->value))) {
      mpfr_get_z(r, *(accessThruMemRef(a->value.expr)->value), GMP_RNDN); /* exact */
      return 1;      
    }    
    break;
  case MPFR:
    if (!mpfr_number_p(a->value.mpfr)) return 0;
    if (!mpfr_integer_p(a->value.mpfr)) return 0;
    mpfr_get_z(r, a->value.mpfr, GMP_RNDN); /* exact */
    return 1;
    break;    
  case SCALEDMPQ:
    if (!scaledMpqIsInteger(a->value.scaledMpq.expo, 
			    a->value.scaledMpq.significand)) return 0;
    mpz_init(num);
    mpz_init(den);
    mpq_get_num(num, a->value.scaledMpq.significand);
    mpq_get_den(den, a->value.scaledMpq.significand);
    if (a->value.scaledMpq.expo >= 0) {
      mpz_mul_2exp(num, num, (mp_bitcnt_t) a->value.scaledMpq.expo);
    } else {
      mpz_mul_2exp(den, den, (mp_bitcnt_t) (-a->value.scaledMpq.expo)); /* Exponent overflow possible */
    }
    mpz_fdiv_q(r, num, den);
    mpz_clear(den);
    mpz_clear(num);
    return 1;
    break;
  }

  mpq_init(q);
  if (!tryConstantToScaledMpq(&EQ, q, a)) {
    mpq_clear(q);
    return 0;
  }
  if (!scaledMpqIsInteger(EQ, q)) {
    mpq_clear(q);
    return 0;
  }
  mpz_init(num);
  mpz_init(den);
  mpq_get_num(num, a->value.scaledMpq.significand);
  mpq_get_den(den, a->value.scaledMpq.significand);
  if (a->value.scaledMpq.expo >= 0) {
    mpz_mul_2exp(num, num, (mp_bitcnt_t) a->value.scaledMpq.expo);
  } else {
    mpz_mul_2exp(den, den, (mp_bitcnt_t) (-a->value.scaledMpq.expo)); /* Exponent overflow possible */
  }
  mpz_fdiv_q(r, num, den);
  mpz_clear(den);
  mpz_clear(num);
  return 1;
}

static inline int tryConstantIntegerBound(mp_bitcnt_t *b, constant_t a) {
  mpz_t z;

  mpz_init(z);
  if (!tryConstantToMpz(z, a)) {
    mpz_clear(z);
    return 0;
  }

  *b = (mp_bitcnt_t) mpz_sizeinbase(z, 2) + 1;
  mpz_clear(z);
  return 1;
}

static inline int tryConstantDenominatorLcm(mpz_t lcm, constant_t a) {
  mpq_t q;
  mp_exp_t E;
  mpz_t n;
  
  mpq_init(q);
  if (!tryConstantToScaledMpq(&E, q, a)) {
    mpq_clear(q);
    return 0;
  }

  mpz_init(n);
  mpq_get_den(n, q);
  if (E < ((mp_exp_t) 0)) {
    mpz_mul_2exp(n, n, ((mp_bitcnt_t) (-E)));
  }
  mpz_lcm(lcm, lcm, n);
  mpz_clear(n);
  mpq_clear(q);
  
  return 1;
}


static inline int tryConstantToInt(int *r, constant_t a) {
  mpq_t q;
  mp_exp_t EQ;
  int res;

  if (a == NULL) return 0;
  switch (a->type) {
  case INTEGER:
    *r = a->value.integer;
    return 1;
    break;
  case EXPRESSION:
    if ((accessThruMemRef(a->value.expr)->nodeType == CONSTANT) &&
	mpfr_number_p(*(accessThruMemRef(a->value.expr)->value))) {
      return mpfr_is_machine_integer(r, *(accessThruMemRef(a->value.expr)->value));
    }    
    break;
  case MPFR:
    return mpfr_is_machine_integer(r, a->value.mpfr);
    break;
  default:
    break;
  }
  
  mpq_init(q);
  if (!tryConstantToScaledMpq(&EQ, q, a)) {
    mpq_clear(q);
    return 0;
  }
  if (!scaledMpqIsInt(&res, EQ, q)) {
    mpq_clear(q);
    return 0;
  }
  *r = res;
  return 1;
}

static inline int tryConstantToUnsignedInt(unsigned int *r, constant_t a) {
  mpq_t q;
  mp_exp_t EQ;
  unsigned int res;

  if (a == NULL) return 0;
  switch (a->type) {
  case INTEGER:
    if (a->value.integer < 0) return 0;
    *r = a->value.integer;
    return 1;
    break;
  case EXPRESSION:
    if ((accessThruMemRef(a->value.expr)->nodeType == CONSTANT) &&
	mpfr_number_p(*(accessThruMemRef(a->value.expr)->value))) {
      return mpfr_is_machine_unsigned_integer(r, *(accessThruMemRef(a->value.expr)->value));
    }    
    break;
  case MPFR:
    return mpfr_is_machine_unsigned_integer(r, a->value.mpfr);
    break;
  default:
    break;
  }
  
  mpq_init(q);
  if (!tryConstantToScaledMpq(&EQ, q, a)) {
    mpq_clear(q);
    return 0;
  }
  if (!scaledMpqIsUnsignedInt(&res, EQ, q)) {
    mpq_clear(q);
    return 0;
  }
  *r = res;
  return 1;
}

void constantFPrintf(FILE *fd, constant_t c) {
  if (c == NULL) {
    sollyaFprintf(fd, "(null)");
    return;
  }
  switch (c->type) {
  case INTEGER:
    sollyaFprintf(fd,"%d",c->value.integer);
    break;
  case EXPRESSION:
    sollyaFprintf(fd,"%b",c->value.expr);
    break;
  case MPFR:
    sollyaFprintf(fd,"%v",c->value.mpfr);
    break;
  case SCALEDMPQ:
    if (c->value.scaledMpq.expo == 0) {
      sollyaFprintf(fd,"%r",
		    c->value.scaledMpq.significand);
    } else {
      sollyaFprintf(fd,"2^(%lld) * %r",
		    (long long int) c->value.scaledMpq.expo, 
		    c->value.scaledMpq.significand);
    }
    break;
  }
}

char *constantToString(constant_t c) {
  char staticStr[8];
  char *str;
  int size, r;

  if (c == NULL) return NULL;
  switch (c->type) {
  case INTEGER:
    size = sollya_snprintf(staticStr,8,"%d",c->value.integer);
    break;
  case EXPRESSION:
    size = sollya_snprintf(staticStr,8,"%b",c->value.expr);
    break;
  case MPFR:
    size = sollya_snprintf(staticStr,8,"%v",c->value.mpfr);
    break;
  case SCALEDMPQ:
    if (c->value.scaledMpq.expo == 0) {
      size = sollya_snprintf(staticStr,8,"%r",
			     c->value.scaledMpq.significand);
    } else {
      size = sollya_snprintf(staticStr,8,"2^(%lld) * %r",
			     (long long int) c->value.scaledMpq.expo, 
			     c->value.scaledMpq.significand);
    }
    break;
  default:
    return NULL;
    break;
  }
  if (size < 0) return NULL;
  str = (char *) safeCalloc(size + 2, sizeof(char));
  switch (c->type) {
  case INTEGER:
    r = sollya_snprintf(str,size,"%d",c->value.integer);
    break;
  case EXPRESSION:
    r = sollya_snprintf(str,size,"%b",c->value.expr);
    break;
  case MPFR:
    r = sollya_snprintf(str,size,"%v",c->value.mpfr);
    break;
  case SCALEDMPQ:
    if (c->value.scaledMpq.expo == 0) {
      r = sollya_snprintf(str,size,"%r",
			  c->value.scaledMpq.significand);
    } else {
      r = sollya_snprintf(str,size,"2^(%lld) * %r",
			  (long long int) c->value.scaledMpq.expo, 
			  c->value.scaledMpq.significand);
    }
    break;
  default:
    safeFree(str);
    return NULL;
    break;
  }
  if (r < 0) {
    safeFree(str);
    return NULL;
  }
  return str;
}

static inline constant_t constantAdd(constant_t a, constant_t b) {
  node *aExpr, *bExpr, *cExpr;
  constant_t res;
  mpq_t aS, bS, cS;
  mp_exp_t EA, EB, EC;
  int cI;
  mpz_t cZ;
  unsigned long int bAbs;
  mpfr_t cM;

  /* Stupid inputs */
  if (a == NULL) return NULL;
  if (b == NULL) return NULL;

  /* Operation on neutral element */
  if (constantIsZero(a, 0)) return constantFromCopy(b);
  if (constantIsZero(b, 0)) return constantFromCopy(a);

  /* Handling expression representations */
  if ((a->type == EXPRESSION) ||
      (b->type == EXPRESSION)) {
    /* One of the two constants is an expression. 
       Get expression representations for both, 
       construct the expression for the operation
       and build the result constant.
    */
    aExpr = constantToExpression(a);
    bExpr = constantToExpression(b);
    cExpr = addMemRef(makeAdd(aExpr, bExpr));
    res = constantFromExpression(cExpr);
    freeThing(cExpr);
    return res;
  }

  /* Here, none of the constants a and b is represented by an
     expression

     Now handle the case when both constants are represented by the
     same type.

  */
  if (a->type == b->type) {
    switch (a->type) {
    case INTEGER:
      if (tryExactIntAddition(&cI, a->value.integer, b->value.integer)) {
	res = constantFromInt(cI);
	return res;
      } 
      mpz_init_set_si(cZ, a->value.integer);
      if (b->value.integer >= 0) {
	bAbs = b->value.integer;
	mpz_add_ui(cZ, cZ, bAbs);
      } else {
	bAbs = - b->value.integer;
	mpz_sub_ui(cZ, cZ, bAbs);
      }
      res = constantFromMpz(cZ);
      mpz_clear(cZ);
      return res;
    break;
    case MPFR:
      mpfr_init2(cM, 12);
      /* Next call may change the precision of cM */
      mpfr_add_exact(cM, a->value.mpfr, b->value.mpfr);
      if (mpfr_number_p(a->value.mpfr) &&
	  mpfr_number_p(b->value.mpfr) && 
	  (!mpfr_number_p(cM))) {
	cExpr = addMemRef(makeAdd(makeConstant(a->value.mpfr),
				  makeConstant(b->value.mpfr)));
	res = constantFromExpression(cExpr);
	freeThing(cExpr);
      } else {
	res = constantFromMpfr(cM);
      }
      mpfr_clear(cM);
      return res;
      break;
    case SCALEDMPQ:
      mpq_init(cS);
      scaledMpqAdd(&EC, cS, 
		   a->value.scaledMpq.expo, a->value.scaledMpq.significand, 
		   b->value.scaledMpq.expo, b->value.scaledMpq.significand);
      res = constantFromScaledMpq(EC, cS);
      mpq_clear(cS);
      return res;
      break;
    default:
      break;
    }
  }

  /* Handle the case when only one of the inputs is a machine integer */
  if (a->type == INTEGER) {
    switch (b->type) {
    case MPFR:
      mpfr_init2(cM, 12);
      /* Next call may change the precision of cM */
      mpfr_add_exact_int(cM, b->value.mpfr, a->value.integer);
      if (mpfr_number_p(b->value.mpfr) &&
	  (!mpfr_number_p(cM))) {
	cExpr = addMemRef(makeAdd(makeConstantInt(a->value.integer),
				  makeConstant(b->value.mpfr)));
	res = constantFromExpression(cExpr);
	freeThing(cExpr);
      } else {
	res = constantFromMpfr(cM);
      }
      mpfr_clear(cM);
      return res;      
      break;
    case SCALEDMPQ:
      mpq_init(cS);
      scaledMpqAddInt(&EC, cS, 
		      b->value.scaledMpq.expo, b->value.scaledMpq.significand, 
		      a->value.integer);
      res = constantFromScaledMpq(EC, cS);
      mpq_clear(cS);
      return res;      
      break;
    default:
      break;
    }
  }
  if (b->type == INTEGER) {
    switch (a->type) {
    case MPFR:
      mpfr_init2(cM, 12);
      /* Next call may change the precision of cM */
      mpfr_add_exact_int(cM, a->value.mpfr, b->value.integer);
      if (mpfr_number_p(a->value.mpfr) &&
	  (!mpfr_number_p(cM))) {
	cExpr = addMemRef(makeAdd(makeConstant(a->value.mpfr),
				  makeConstantInt(b->value.integer)));
	res = constantFromExpression(cExpr);
	freeThing(cExpr);
      } else {
	res = constantFromMpfr(cM);
      }
      mpfr_clear(cM);
      return res;      
      break;
    case SCALEDMPQ:
      mpq_init(cS);
      scaledMpqAddInt(&EC, cS, 
		      a->value.scaledMpq.expo, a->value.scaledMpq.significand, 
		      b->value.integer);
      res = constantFromScaledMpq(EC, cS);
      mpq_clear(cS);
      return res;      
      break;
    default:
      break;
    }
  }

  /* Here, none of the constants a and b is represented by an
     expression and they are of different types. 

     Convert both constants to scaled rational form, perform the
     operation and convert back.

     If one of the conversions to scaled rational form fails, convert
     to expression form and work with expression form.
  */
  mpq_init(aS);
  mpq_init(bS);
  if (tryConstantToScaledMpq(&EA, aS, a) &&
      tryConstantToScaledMpq(&EB, bS, b)) {
    mpq_init(cS);
    scaledMpqAdd(&EC, cS, EA, aS, EB, bS);
    res = constantFromScaledMpq(EC, cS);
    mpq_clear(aS);
    mpq_clear(bS);
    mpq_clear(cS);
    return res;
  } 
  mpq_clear(aS);
  mpq_clear(bS);

  /* Could not convert to scaled rational form. Convert to
     expressions. 
  */
  aExpr = constantToExpression(a);
  bExpr = constantToExpression(b);
  cExpr = addMemRef(makeAdd(aExpr, bExpr));
  res = constantFromExpression(cExpr);
  freeThing(cExpr);
  return res;  
}

static inline constant_t constantSub(constant_t a, constant_t b) {
  node *aExpr, *bExpr, *cExpr;
  constant_t res;
  mpq_t aS, bS, cS;
  mp_exp_t EA, EB, EC;
  int cI;
  mpz_t cZ;
  unsigned long int bAbs;
  mpfr_t cM;

  /* Stupid inputs */
  if (a == NULL) return NULL;
  if (b == NULL) return NULL;

  /* Operation on neutral element */
  if (constantIsZero(a, 0)) return constantNeg(b);
  if (constantIsZero(b, 0)) return constantFromCopy(a);

  /* Handling expression representations */
  if ((a->type == EXPRESSION) ||
      (b->type == EXPRESSION)) {
    /* One of the two constants is an expression. 
       Get expression representations for both, 
       construct the expression for the operation
       and build the result constant.
    */
    aExpr = constantToExpression(a);
    bExpr = constantToExpression(b);
    cExpr = addMemRef(makeSub(aExpr, bExpr));
    res = constantFromExpression(cExpr);
    freeThing(cExpr);
    return res;
  }

  /* Here, none of the constants a and b is represented by an
     expression

     Now handle the case when both constants are represented by the
     same type.

  */
  if (a->type == b->type) {
    switch (a->type) {
    case INTEGER:
      if (tryExactIntSubtraction(&cI, a->value.integer, b->value.integer)) {
	res = constantFromInt(cI);
	return res;
      } 
      mpz_init_set_si(cZ, a->value.integer);
      if (b->value.integer >= 0) {
	bAbs = b->value.integer;
	mpz_sub_ui(cZ, cZ, bAbs);
      } else {
	bAbs = - b->value.integer;
	mpz_add_ui(cZ, cZ, bAbs);
      }
      res = constantFromMpz(cZ);
      mpz_clear(cZ);
      return res;
    break;
    case MPFR:
      mpfr_init2(cM, 12);
      /* Next call may change the precision of cM */
      mpfr_sub_exact(cM, a->value.mpfr, b->value.mpfr);
      if (mpfr_number_p(a->value.mpfr) &&
	  mpfr_number_p(b->value.mpfr) && 
	  (!mpfr_number_p(cM))) {
	cExpr = addMemRef(makeSub(makeConstant(a->value.mpfr),
				  makeConstant(b->value.mpfr)));
	res = constantFromExpression(cExpr);
	freeThing(cExpr);
      } else {
	res = constantFromMpfr(cM);
      }
      mpfr_clear(cM);
      return res;
      break;
    case SCALEDMPQ:
      mpq_init(cS);
      scaledMpqSub(&EC, cS, 
		   a->value.scaledMpq.expo, a->value.scaledMpq.significand, 
		   b->value.scaledMpq.expo, b->value.scaledMpq.significand);
      res = constantFromScaledMpq(EC, cS);
      mpq_clear(cS);
      return res;
      break;
    default:
      break;
    }
  }

    /* Handle the case when only one of the inputs is a machine integer */
  if (a->type == INTEGER) {
    switch (b->type) {
    case MPFR:
      mpfr_init2(cM, 12);
      /* Next call may change the precision of cM */
      mpfr_int_sub_exact(cM, a->value.integer, b->value.mpfr);
      if (mpfr_number_p(b->value.mpfr) &&
	  (!mpfr_number_p(cM))) {
	cExpr = addMemRef(makeSub(makeConstantInt(a->value.integer),
				  makeConstant(b->value.mpfr)));
	res = constantFromExpression(cExpr);
	freeThing(cExpr);
      } else {
	res = constantFromMpfr(cM);
      }
      mpfr_clear(cM);
      return res;      
      break;
    case SCALEDMPQ:
      mpq_init(cS);
      scaledMpqIntSub(&EC, cS, 
		      a->value.integer,
		      b->value.scaledMpq.expo, b->value.scaledMpq.significand);
      res = constantFromScaledMpq(EC, cS);
      mpq_clear(cS);
      return res;      
      break;
    default:
      break;
    }
  }
  if (b->type == INTEGER) {
    switch (a->type) {
    case MPFR:
      mpfr_init2(cM, 12);
      /* Next call may change the precision of cM */
      mpfr_sub_exact_int(cM, a->value.mpfr, b->value.integer);
      if (mpfr_number_p(a->value.mpfr) &&
	  (!mpfr_number_p(cM))) {
	cExpr = addMemRef(makeSub(makeConstant(a->value.mpfr),
				  makeConstantInt(b->value.integer)));
	res = constantFromExpression(cExpr);
	freeThing(cExpr);
      } else {
	res = constantFromMpfr(cM);
      }
      mpfr_clear(cM);
      return res;      
      break;
    case SCALEDMPQ:
      mpq_init(cS);
      scaledMpqSubInt(&EC, cS, 
		      a->value.scaledMpq.expo, a->value.scaledMpq.significand, 
		      b->value.integer);
      res = constantFromScaledMpq(EC, cS);
      mpq_clear(cS);
      return res;      
      break;
    default:
      break;
    }
  }

  /* Here, none of the constants a and b is represented by an
     expression and they are of different types. 

     Convert both constants to scaled rational form, perform the
     operation and convert back.

     If one of the conversions to scaled rational form fails, convert
     to expression form and work with expression form.
  */
  mpq_init(aS);
  mpq_init(bS);
  if (tryConstantToScaledMpq(&EA, aS, a) &&
      tryConstantToScaledMpq(&EB, bS, b)) {
    mpq_init(cS);
    scaledMpqSub(&EC, cS, EA, aS, EB, bS);
    res = constantFromScaledMpq(EC, cS);
    mpq_clear(aS);
    mpq_clear(bS);
    mpq_clear(cS);
    return res;
  } 
  mpq_clear(aS);
  mpq_clear(bS);

  /* Could not convert to scaled rational form. Convert to
     expressions. 
  */
  aExpr = constantToExpression(a);
  bExpr = constantToExpression(b);
  cExpr = addMemRef(makeSub(aExpr, bExpr));
  res = constantFromExpression(cExpr);
  freeThing(cExpr);
  return res;  
}

static inline constant_t constantMul(constant_t a, constant_t b) {
  node *aExpr, *bExpr, *cExpr;
  constant_t res;
  mpq_t aS, bS, cS;
  mp_exp_t EA, EB, EC;
  int cI;
  mpz_t cZ;
  mpfr_t cM;

  /* Stupid inputs */
  if (a == NULL) return NULL;
  if (b == NULL) return NULL;

  /* Operation on neutral element */
  if (constantIsOne(a, 0)) return constantFromCopy(b);
  if (constantIsOne(b, 0)) return constantFromCopy(a);

  /* Handling expression representations */
  if ((a->type == EXPRESSION) ||
      (b->type == EXPRESSION)) {
    /* One of the two constants is an expression. 
       Get expression representations for both, 
       construct the expression for the operation
       and build the result constant.
    */
    aExpr = constantToExpression(a);
    bExpr = constantToExpression(b);
    cExpr = addMemRef(makeMul(aExpr, bExpr));
    res = constantFromExpression(cExpr);
    freeThing(cExpr);
    return res;
  }

  /* Here, none of the constants a and b is represented by an
     expression

     Now handle the case when both constants are represented by the
     same type.

  */
  if (a->type == b->type) {
    switch (a->type) {
    case INTEGER:
      if (tryExactIntMultiplication(&cI, a->value.integer, b->value.integer)) {
	res = constantFromInt(cI);
	return res;
      } 
      mpz_init_set_si(cZ, a->value.integer);
      mpz_mul_si(cZ, cZ, b->value.integer);
      res = constantFromMpz(cZ);
      mpz_clear(cZ);
      return res;
    break;
    case MPFR:
      mpfr_init2(cM, 12);
      /* Next call may change the precision of cM */
      mpfr_mul_exact(cM, a->value.mpfr, b->value.mpfr);
      if (mpfr_number_p(a->value.mpfr) &&
	  mpfr_number_p(b->value.mpfr) && 
	  (!mpfr_number_p(cM))) {
	cExpr = addMemRef(makeMul(makeConstant(a->value.mpfr),
				  makeConstant(b->value.mpfr)));
	res = constantFromExpression(cExpr);
	freeThing(cExpr);
      } else {
	res = constantFromMpfr(cM);
      }
      mpfr_clear(cM);
      return res;
      break;
    case SCALEDMPQ:
      mpq_init(cS);
      scaledMpqMul(&EC, cS, 
		   a->value.scaledMpq.expo, a->value.scaledMpq.significand, 
		   b->value.scaledMpq.expo, b->value.scaledMpq.significand);
      res = constantFromScaledMpq(EC, cS);
      mpq_clear(cS);
      return res;
      break;
    default:
      break;
    }
  }

  /* Handle the case when only one of the inputs is a machine integer */
  if (a->type == INTEGER) {
    switch (b->type) {
    case MPFR:
      mpfr_init2(cM, 12);
      /* Next call may change the precision of cM */
      mpfr_mul_exact_int(cM, b->value.mpfr, a->value.integer);
      if (mpfr_number_p(b->value.mpfr) && 
	  (!mpfr_number_p(cM))) {
	cExpr = addMemRef(makeMul(makeConstantInt(a->value.integer),
				  makeConstant(b->value.mpfr)));
	res = constantFromExpression(cExpr);
	freeThing(cExpr);
      } else {
	res = constantFromMpfr(cM);
      }
      mpfr_clear(cM);
      return res;      
      break;
    case SCALEDMPQ:
      mpq_init(cS);
      scaledMpqMulInt(&EC, cS, 
		      b->value.scaledMpq.expo, b->value.scaledMpq.significand, 
		      a->value.integer);
      res = constantFromScaledMpq(EC, cS);
      mpq_clear(cS);
      return res;      
      break;
    default:
      break;
    }
  }
  if (b->type == INTEGER) {
    switch (a->type) {
    case MPFR:
      mpfr_init2(cM, 12);
      /* Next call may change the precision of cM */
      mpfr_mul_exact_int(cM, a->value.mpfr, b->value.integer);
      if (mpfr_number_p(a->value.mpfr) &&
	  (!mpfr_number_p(cM))) {
	cExpr = addMemRef(makeMul(makeConstant(a->value.mpfr),
				  makeConstantInt(b->value.integer)));
	res = constantFromExpression(cExpr);
	freeThing(cExpr);
      } else {
	res = constantFromMpfr(cM);
      }
      mpfr_clear(cM);
      return res;      
      break;
    case SCALEDMPQ:
      mpq_init(cS);
      scaledMpqMulInt(&EC, cS, 
		      a->value.scaledMpq.expo, a->value.scaledMpq.significand, 
		      b->value.integer);
      res = constantFromScaledMpq(EC, cS);
      mpq_clear(cS);
      return res;      
      break;
    default:
      break;
    }
  }

  /* Here, none of the constants a and b is represented by an
     expression and they are of different types. 

     Convert both constants to scaled rational form, perform the
     operation and convert back.

     If one of the conversions to scaled rational form fails, convert
     to expression form and work with expression form.
  */
  mpq_init(aS);
  mpq_init(bS);
  if (tryConstantToScaledMpq(&EA, aS, a) &&
      tryConstantToScaledMpq(&EB, bS, b)) {
    mpq_init(cS);
    scaledMpqMul(&EC, cS, EA, aS, EB, bS);
    res = constantFromScaledMpq(EC, cS);
    mpq_clear(aS);
    mpq_clear(bS);
    mpq_clear(cS);
    return res;
  } 
  mpq_clear(aS);
  mpq_clear(bS);

  /* Could not convert to scaled rational form. Convert to
     expressions. 
  */
  aExpr = constantToExpression(a);
  bExpr = constantToExpression(b);
  cExpr = addMemRef(makeMul(aExpr, bExpr));
  res = constantFromExpression(cExpr);
  freeThing(cExpr);
  return res;  
}

static inline constant_t constantDiv(constant_t a, constant_t b) {
  node *aExpr, *bExpr, *cExpr;
  constant_t res;
  mpq_t aS, bS, cS, cQ;
  mp_exp_t EA, EB, EC;
  int cI;
  unsigned long int bAbs;

  /* Stupid inputs */
  if (a == NULL) return NULL;
  if (b == NULL) return NULL;

  /* Operation on neutral element */
  if (constantIsOne(b, 0)) return constantFromCopy(a);

  /* Handling expression representations */
  if ((a->type == EXPRESSION) ||
      (b->type == EXPRESSION)) {
    /* One of the two constants is an expression. 
       Get expression representations for both, 
       construct the expression for the operation
       and build the result constant.
    */
    aExpr = constantToExpression(a);
    bExpr = constantToExpression(b);
    cExpr = addMemRef(makeDiv(aExpr, bExpr));
    res = constantFromExpression(cExpr);
    freeThing(cExpr);
    return res;
  }

  /* Here, none of the constants a and b is represented by an
     expression

     Now handle the case when both constants are represented by the
     same type.

  */
  if (a->type == b->type) {
    switch (a->type) {
    case INTEGER:
      if (b->value.integer != 0) {
	if (tryExactIntDivision(&cI, a->value.integer, b->value.integer)) {
	  res = constantFromInt(cI);
	  return res;
	} 
	mpq_init(cQ);
	if (b->value.integer >= 0) {
	  bAbs = b->value.integer;
	  mpq_set_si(cQ,a->value.integer,bAbs);
	  mpq_canonicalize(cQ);
	} else {
	  bAbs = - b->value.integer;
	  mpq_set_si(cQ,a->value.integer,bAbs);
	  mpq_canonicalize(cQ);
	  mpq_neg(cQ, cQ);
	}
	res = constantFromMpq(cQ);
	mpq_clear(cQ);
	return res;
      }
    break;
    case MPFR:
      /* Do nothing for this case, let the common scaled MPQ logic handle it below. 
	 For cases when the MPFR values are Inf or NaN, the MPQ logic will handle the
	 case over to expression logic.
      */
      break;
    case SCALEDMPQ:
      mpq_init(cS);
      if (tryScaledMpqDiv(&EC, cS, 
			  a->value.scaledMpq.expo, a->value.scaledMpq.significand, 
			  b->value.scaledMpq.expo, b->value.scaledMpq.significand)) {
	res = constantFromScaledMpq(EC, cS);
	mpq_clear(cS);
	return res;
      } 
      mpq_clear(cS);
      break;
    default:
      break;
    }
  }

  /* Here, none of the constants a and b is represented by an
     expression and they are of different types or the division
     on the same type did not work.

     Convert both constants to scaled rational form, perform the
     operation and convert back.

     If one of the conversions to scaled rational form fails, convert
     to expression form and work with expression form.
  */
  mpq_init(aS);
  mpq_init(bS);
  if (tryConstantToScaledMpq(&EA, aS, a) &&
      tryConstantToScaledMpq(&EB, bS, b)) {
    mpq_init(cS);
    if (tryScaledMpqDiv(&EC, cS, EA, aS, EB, bS)) {
      res = constantFromScaledMpq(EC, cS);
      mpq_clear(aS);
      mpq_clear(bS);
      mpq_clear(cS);
      return res;
    }
    mpq_clear(cS);
  } 
  mpq_clear(aS);
  mpq_clear(bS);

  /* Could not convert to scaled rational form or the division did not work. 
     Convert to expressions. 
  */
  aExpr = constantToExpression(a);
  bExpr = constantToExpression(b);
  cExpr = addMemRef(makeDiv(aExpr, bExpr));
  res = constantFromExpression(cExpr);
  freeThing(cExpr);
  return res;  
}

static inline constant_t constantPow(constant_t a, constant_t b) {
  node *aExpr, *bExpr, *cExpr;
  constant_t res;
  mpq_t aS, bS, cS;
  mp_exp_t EA, EB, EC;
  int cI;
  mpfr_t cM;
  unsigned long int bAbs;
  int bb;
  mpz_t cZ;

  /* Stupid inputs */
  if (a == NULL) return NULL;
  if (b == NULL) return NULL;

  /* Operation on neutral element */
  if (constantIsOne(a, 0)) return constantFromCopy(a);
  if (constantIsOne(b, 0)) return constantFromCopy(a);
  if (constantIsZero(b, 0)) return constantFromInt(1);

  /* Handling expression representations */
  if ((a->type == EXPRESSION) ||
      (b->type == EXPRESSION)) {
    /* One of the two constants is an expression. 
       Get expression representations for both, 
       construct the expression for the operation
       and build the result constant.
    */
    aExpr = constantToExpression(a);
    bExpr = constantToExpression(b);
    cExpr = addMemRef(makePow(aExpr, bExpr));
    res = constantFromExpression(cExpr);
    freeThing(cExpr);
    return res;
  }

  /* Here, none of the constants a and b is represented by an
     expression

     Now handle the case when both constants are represented by the
     same type.

  */
  if (a->type == b->type) {
    switch (a->type) {
    case INTEGER:
      if (tryExactIntPower(&cI, a->value.integer, b->value.integer)) {
	res = constantFromInt(cI);
	return res;
      } 
      if (b->value.integer > 0) {
	bAbs = b->value.integer;
	mpz_init_set_si(cZ, a->value.integer);
	mpz_pow_ui(cZ, cZ, bAbs);
	res = constantFromMpz(cZ);
	mpz_clear(cZ);
	return res;
      }
    break;
    case MPFR:
      if (mpfr_number_p(a->value.mpfr) &&
	  mpfr_number_p(b->value.mpfr) &&
	  mpfr_integer_p(b->value.mpfr) && 
	  (mpfr_cmp_si(b->value.mpfr, 1) >= 0) &&
	  mpfr_fits_ulong_p(b->value.mpfr, GMP_RNDN)) {
	bAbs = mpfr_get_ui(b->value.mpfr, GMP_RNDN); /* exact */
	mpfr_init2(cM, 12);
	if (try_mpfr_pow_exact(cM, a->value.mpfr, bAbs)) {
	  if (mpfr_number_p(a->value.mpfr) && 
	      (!mpfr_number_p(cM))) {
	    cExpr = addMemRef(makePow(makeConstant(a->value.mpfr),
				      makeConstantUnsignedInt(bAbs)));
	    res = constantFromExpression(cExpr);
	    freeThing(cExpr);
	  } else {
	    res = constantFromMpfr(cM);
	  }
	  mpfr_clear(cM);
	  return res;
	}
	mpfr_clear(cM);
      }
      break;
    case SCALEDMPQ:
      mpq_init(cS);
      if (tryScaledMpqPow(&EC, cS, 
			  a->value.scaledMpq.expo, a->value.scaledMpq.significand, 
			  b->value.scaledMpq.expo, b->value.scaledMpq.significand)) {
	res = constantFromScaledMpq(EC, cS);
	mpq_clear(cS);
	return res;
      } 
      mpq_clear(cS);
      break;
    default:
      break;
    }
  }

  /* Handle the case when b is a machine integer and a is not */
  if (b->type == INTEGER) {
    switch (a->type) {
    case MPFR:
      if (mpfr_number_p(a->value.mpfr) && (b->value.integer >= 1)) {
	bAbs = b->value.integer;
	bb = bAbs;
	if ((bAbs >= 1) && (bb == b->value.integer)) {
	  mpfr_init2(cM, 12);
	  if (try_mpfr_pow_exact(cM, a->value.mpfr, bAbs)) {
	    if (mpfr_number_p(a->value.mpfr) &&
		(!mpfr_number_p(cM))) {
	      cExpr = addMemRef(makePow(makeConstant(a->value.mpfr),
					makeConstantInt(b->value.integer)));
	      res = constantFromExpression(cExpr);
	      freeThing(cExpr);
	    } else {
	      res = constantFromMpfr(cM);
	    }
	    mpfr_clear(cM);
	    return res;
	  }
	  mpfr_clear(cM);
	}
      }
      break;
    case SCALEDMPQ:
      mpq_init(cS);
      if (tryScaledMpqPowInt(&EC, cS, 
			     a->value.scaledMpq.expo, a->value.scaledMpq.significand, 
			     b->value.integer)) {
	res = constantFromScaledMpq(EC, cS);
	mpq_clear(cS);
	return res;
      } 
      mpq_clear(cS);
      break;
    default:
      break;
    }
  }

  /* Here, none of the constants a and b is represented by an
     expression and they are of different types or the powering
     on the same type did not work.

     Convert both constants to scaled rational form, perform the
     operation and convert back.

     If one of the conversions to scaled rational form fails, convert
     to expression form and work with expression form.
  */
  mpq_init(aS);
  mpq_init(bS);
  if (tryConstantToScaledMpq(&EA, aS, a) &&
      tryConstantToScaledMpq(&EB, bS, b)) {
    mpq_init(cS);
    if (tryScaledMpqPow(&EC, cS, EA, aS, EB, bS)) {
      res = constantFromScaledMpq(EC, cS);
      mpq_clear(aS);
      mpq_clear(bS);
      mpq_clear(cS);
      return res;
    }
    mpq_clear(cS);
  } 
  mpq_clear(aS);
  mpq_clear(bS);

  /* Could not convert to scaled rational form or the powering did not work. 
     Convert to expressions. 
  */
  aExpr = constantToExpression(a);
  bExpr = constantToExpression(b);
  cExpr = addMemRef(makePow(aExpr, bExpr));
  res = constantFromExpression(cExpr);
  freeThing(cExpr);
  return res;  
}

static inline constant_t constantNeg(constant_t a) {
  constant_t res;
  int cI;
  mpz_t cZ;
  node *cExpr;
  mpfr_t cM;
  mpq_t cQ;

  if (a == NULL) return NULL;
  switch (a->type) {
  case INTEGER:
    if (tryExactIntNegation(&cI, a->value.integer)) {
      res = constantFromInt(cI);
      return res;
    } 
    mpz_init_set_si(cZ, a->value.integer);
    mpz_neg(cZ,cZ);
    res = constantFromMpz(cZ);
    mpz_clear(cZ);
    return res;
    break;
  case EXPRESSION:
    cExpr = addMemRef(makeNeg(copyThing(a->value.expr)));
    res = constantFromExpression(cExpr);
    freeThing(cExpr);
    return res;
    break;
  case MPFR:
    mpfr_init2(cM, mpfr_get_prec(a->value.mpfr));
    mpfr_neg(cM, a->value.mpfr, GMP_RNDN); /* exact */
    res = constantFromMpfr(cM);
    mpfr_clear(cM);
    return res;
    break;
  case SCALEDMPQ:
    mpq_init(cQ);
    mpq_neg(cQ, a->value.scaledMpq.significand);
    res = constantFromScaledMpq(a->value.scaledMpq.expo, cQ);
    mpq_clear(cQ);
    return res;
    break;
  }

  return NULL;
}

/* Computes q and r such that 2^14 * q + r = a and q = floor(a *
   2^-14) 
*/
static inline void constantCutTwo14(constant_t *q, constant_t *r, constant_t a) {
  int qI, rI;
  node *qN, *rN;
  mpfr_t qM, rM;
  mp_exp_t qE, rE;
  mpq_t qQ, rQ;

  /* Handle stupid input */
  if (a == NULL) {
    *q = NULL;
    *r = NULL;
    return;
  }

  switch (a->type) {
  case INTEGER:
    /* The int type is guaranteed to have at least 16 bits. So a shift
       by 14 bits is guaranteed to work.
    */
    qI = a->value.integer >> 14;
    rI = a->value.integer - (qI << 14);
    *q = constantFromInt(qI);
    *r = constantFromInt(rI);
    break;
  case EXPRESSION:
    qN = addMemRef(makeFloor(makeMul(copyThing(a->value.expr),
				     makePow(makeConstantInt(2),
					     makeConstantInt(-14)))));
    rN = addMemRef(makeSub(copyThing(a->value.expr),
			   makeMul(makePow(makeConstantInt(2),
					   makeConstantInt(14)),
				   copyThing(qN))));
    *q = constantFromExpression(qN);
    *r = constantFromExpression(rN);
    freeThing(qN);
    freeThing(rN);
    break;
  case MPFR:
    mpfr_init2(qM, mpfr_get_prec(a->value.mpfr));
    mpfr_init2(rM, mpfr_get_prec(a->value.mpfr));
    mpfr_div_2si(qM, a->value.mpfr, 14, GMP_RNDN); /* exact */
    mpfr_floor(qM, qM); /* exact */
    mpfr_mul_2si(rM, qM, 14, GMP_RNDN); /* exact */
    mpfr_sub(rM, a->value.mpfr, rM, GMP_RNDN); /* exact, Sterbenz */
    *q = constantFromMpfr(qM);
    *r = constantFromMpfr(rM);
    mpfr_clear(rM);
    mpfr_clear(qM);
    break;
  case SCALEDMPQ: 
    mpq_init(qQ);
    mpq_init(rQ);
    scaledMpqFloor(&qE, qQ, a->value.scaledMpq.expo - 14, a->value.scaledMpq.significand);
    scaledMpqSub(&rE, rQ, a->value.scaledMpq.expo, a->value.scaledMpq.significand, qE + 14, qQ);
    *q = constantFromScaledMpq(qE, qQ);
    *r = constantFromScaledMpq(rE, rQ);
    mpq_clear(rQ);
    mpq_clear(qQ);
    break;
  }
}

static inline constant_t constantGcd(constant_t a, constant_t b) {
  constant_t c;
  mpq_t aq, bq, cq;

  /* Handle stupid input */
  if (a == NULL) return NULL;
  if (b == NULL) return NULL;

  /* General case: we define the gcd of constants to be 1, unless
     both constants are rationals, in which case we define

     gcd(a/b, c/d) = gcd(a*d, c*b)/(b*d).

  */
  mpq_init(aq);
  if (tryConstantToMpq(aq, a)) {
    mpq_init(bq);
    if (tryConstantToMpq(bq, b)) {
      mpq_init(cq);
      sollya_mpq_gcd(cq, aq, bq);
      c = constantFromMpq(cq);
      mpq_clear(cq);
    } else {
      c = constantFromInt(1);
    }
    mpq_clear(bq);
  } else {
    c = constantFromInt(1);
  }
  mpq_clear(aq);

  return c;
}

static inline void constantEvalMpfr(mpfr_t rop, constant_t c) {
  mp_prec_t p;
  mpfr_t cutoff;
  int res;
  sollya_mpfi_t y;

 /* Handle stupid input */
  if (c == NULL) {
    mpfr_set_nan(rop);
    return;
  }

  /* Evaluation depending on the representation type */
  switch (c->type) {
  case INTEGER:
    mpfr_set_si(rop, c->value.integer, GMP_RNDN);
    break;
  case EXPRESSION:
    mpfr_set_nan(rop);
    mpfr_init2(cutoff, 12);
    mpfr_set_si(cutoff, 0, GMP_RNDN);
    res = evaluateFaithfulWithCutOffFast(rop, c->value.expr, NULL, rop, cutoff, mpfr_get_prec(rop) + 10);
    if (res == 0) {
      sollya_mpfi_init2(y, mpfr_get_prec(rop) + 10);
      evaluateConstantExpressionToInterval(y, c->value.expr);
      if (!(sollya_mpfi_has_infinity(y) ||
	    sollya_mpfi_has_nan(y))) {
	mpfr_set_si(rop, 0, GMP_RNDN);
      }
      sollya_mpfi_clear(y);
    }
    mpfr_clear(cutoff);
    break;
  case MPFR:
    mpfr_set(rop, c->value.mpfr, GMP_RNDN);
    break;
  case SCALEDMPQ: 
    p = mpfr_get_prec(rop);
    mpfr_set_prec(rop, p + 10);
    mpfr_set_z_2exp(rop, mpq_numref(c->value.scaledMpq.significand), c->value.scaledMpq.expo, GMP_RNDN);   
    mpfr_div_z(rop, rop, mpq_denref(c->value.scaledMpq.significand), GMP_RNDN);
    mpfr_prec_round(rop, p, GMP_RNDN);
    break;
  }
}

static inline void constantEvalMpfi(sollya_mpfi_t rop, constant_t c) {
  mp_prec_t p;

 /* Handle stupid input */
  if (c == NULL) {
    sollya_mpfi_set_nan(rop);
    return;
  }

  /* Evaluation depending on the representation type */
  switch (c->type) {
  case INTEGER:
    sollya_mpfi_set_si(rop, c->value.integer);
    break;
  case EXPRESSION:
    evaluateConstantExpressionToInterval(rop, c->value.expr);
    break;
  case MPFR:
    sollya_mpfi_set_fr(rop, c->value.mpfr);
    break;
  case SCALEDMPQ: 
    p = sollya_mpfi_get_prec(rop);
    sollya_mpfi_set_prec(rop, p + 10);
    sollya_mpfi_set_z_2exp(rop, mpq_numref(c->value.scaledMpq.significand), c->value.scaledMpq.expo);   
    sollya_mpfi_div_z(rop, rop, mpq_denref(c->value.scaledMpq.significand));
    sollya_mpfi_prec_round(rop, p);
    break;
  }
}

static inline constant_t constantRoundDyadic(constant_t c, mp_prec_t prec) {
  mpfr_t v;
  constant_t res;

  /* Handle stupid inputs */
  if (c == NULL) return NULL;

  /* If c is already dyadic, return a copy of c */
  if (constantIsDyadic(c, 0)) return constantFromCopy(c);

  /* Here, c is not dyadic. Allocate an MPFR variable with precision
     prec and evaluate c to that precision.
  */
  mpfr_init2(v, prec);
  constantEvalMpfr(v, c);
  if (mpfr_number_p(v)) {
    res = constantFromMpfr(v);
  } else {
    res = constantFromCopy(c);
  }
  mpfr_clear(v);

  return res;
}

static inline constant_t constantRoundRational(constant_t c, mp_prec_t prec) {
  mpfr_t v;
  constant_t res;

  /* Handle stupid inputs */
  if (c == NULL) return NULL;

  /* If c is already rational, return a copy of c */
  if (constantIsRational(c, 0)) return constantFromCopy(c);

  /* Here, c is not rational. Allocate an MPFR variable with precision
     prec and evaluate c to that precision. 
  */
  mpfr_init2(v, prec);
  constantEvalMpfr(v, c);
  if (mpfr_number_p(v)) {
    res = constantFromMpfr(v);
  } else {
    res = constantFromCopy(c);
  }
  mpfr_clear(v);

  return res;
}

static inline constant_t constantRound(constant_t c, mp_prec_t prec) {
  mpfr_t v;
  constant_t res;

  /* Handle stupid inputs */
  if (c == NULL) return NULL;

  /* If c holds on prec bits, return a copy of c */
  if (constantHoldsOnPrecBits(c, prec, 0)) return constantFromCopy(c);

  /* Here, c does not hold on prec bits. Allocate an MPFR variable
     with precision prec and evaluate c to that precision.
  */
  mpfr_init2(v, prec);
  constantEvalMpfr(v, c);
  if (mpfr_number_p(v)) {
    res = constantFromMpfr(v);
  } else {
    res = constantFromCopy(c);
  }
  mpfr_clear(v);

  return res;
}

static inline int constantReferencesExpression(constant_t c, node *expr) {
  /* Handle stupid input */
  if (c == NULL) return 0;

  /* Evaluation depending on the representation type */
  switch (c->type) {
  case EXPRESSION:
    if (c->value.expr == expr) return 1;
    if (accessThruMemRef(c->value.expr) == expr) return 1;
    return 0;
    break;
  case INTEGER:
  case MPFR:
  case SCALEDMPQ: 
    return 0;
    break;
  }

  return 0;
}

static inline int constantUsesExpressionConstant(constant_t c) {
  /* Handle stupid input */
  if (c == NULL) return 0;

  /* Evaluation depending on the representation type */
  return (c->type == EXPRESSION);
}


static inline uint64_t constantHash(constant_t c) {
  uint64_t hash;
  
  /* Handle stupid input */
  if (c == NULL) return hashPointer(NULL);

  /* Check if the hash has already been computed */
  if (c->hash.hasHash) {
    return c->hash.hash;
  }
  
  /* Compute the hash */
  hash = hashInt((int) (c->type));
  switch (c->type) {
  case INTEGER:
    hash = hashCombine(hash, hashInt(c->value.integer));
    break;
  case EXPRESSION:
    hash = hashCombine(hash, hashThingNoPolynomialHandling(c->value.expr));
    break;
  case MPFR:
    hash = hashCombine(hash, hashMpfr(c->value.mpfr));
    break;
  case SCALEDMPQ:
    hash = hashCombine(hash,
		       hashCombine(hashInt64((int64_t) (c->value.scaledMpq.expo)),
				   hashMpq(c->value.scaledMpq.significand)));
    break;
  }

  /* Cache the result */
  c->hash.hash = hash;
  c->hash.hasHash = 1;
  
  /* Return the hash */
  return hash;
}

/* End of part for constants */

/* Start of part for sparse polynomials */

static inline sparse_polynomial_t __polynomialGetSparsePolynomial(polynomial_t);
static inline int sparsePolynomialPowConstant(sparse_polynomial_t *, sparse_polynomial_t, constant_t);
static inline int sparsePolynomialGetDegreeAsInt(sparse_polynomial_t);
static inline int sparsePolynomialCoefficientsAreRational(sparse_polynomial_t, int);
static inline sparse_polynomial_t sparsePolynomialPowUnsignedInt(sparse_polynomial_t, unsigned int);
static inline constant_t sparsePolynomialGetIthCoefficientAsConstantIntIndex(sparse_polynomial_t, int);
static inline int sparsePolynomialEvalMpz(mpz_t, sparse_polynomial_t, mpz_t);

static inline sparse_polynomial_t __sparsePolynomialAllocate() {
  return (sparse_polynomial_t) safeMalloc(sizeof(struct __sparse_polynomial_struct_t));
}

static inline void __sparsePolynomialFreeMem(sparse_polynomial_t p) {
  safeFree(p);
}

static inline void sparsePolynomialInitializeCaches() {
  constantInitializeCaches();
}

static inline void sparsePolynomialFreeCaches() {
  constantFreeCaches();
}

static inline void __sparsePolynomialAdjustDegree(sparse_polynomial_t p) {
  unsigned int i, k;

  if (p == NULL) return;
  for (k=0u,i=p->monomialCount-1u;i>=1u;i--,k++) {
    if (!constantIsZero(p->coeffs[i],0)) break;
  }
  if (k == 0u) return;
  for (i=p->monomialCount-k;i<p->monomialCount;i++) {
    constantFree(p->coeffs[i]);
    constantFree(p->monomialDegrees[i]);
  }
  p->monomialCount -= k;
  constantFree(p->deg);
  p->deg = constantFromCopy(p->monomialDegrees[p->monomialCount-1u]);
  p->coeffs = (constant_t *) safeRealloc(p->coeffs, 
					 ((size_t) (p->monomialCount)) * sizeof(constant_t));
  p->monomialDegrees = (constant_t *) safeRealloc(p->monomialDegrees, 
						  ((size_t) (p->monomialCount)) * sizeof(constant_t));
}

static inline int __sparsePolynomialAscendingDegrees(sparse_polynomial_t p, int defVal) {
  unsigned int i;
  int ord;
 
  if (p == NULL) return defVal;
  if (p->monomialCount <= 1u) return 1;
  for (i=1u;i<p->monomialCount;i++) {
    ord = constantIsGreater(p->monomialDegrees[i], p->monomialDegrees[i-1u], 19);
    if (ord == 19) return defVal;
    if (!ord) return 0;
  }
  return 1;
}

/* Given a array a of n constants, sorted in ascending order, find
   the lowest index of an element greater than or equal to d. 

   If n is zero, return n.
   Otherwise:

      - if the 0-th (first) element of the array is greater than or
        equal to d, return 0.

      - if the (n-1)-th (last) element of the array is less than d,
        return n.

      - otherwise return the lowest index of an element greater than
        or equal to d. 

   The algorithm is just supposed to work; it can have a ridiculous
   complexity.

*/
static inline unsigned int __sparsePolynomialFindDegreeNaive(constant_t d, 
						       constant_t *a, 
						       unsigned int n) {
  unsigned int i;
  
  if (n == 0u) return n;
  for (i=0u; i<n; i++) {
    if (constantIsGreaterOrEqual(a[i],d,0)) return i;
  }
  return n;
}

/* Given a array a of n constants, sorted in ascending order, find
   the lowest index of an element greater than or equal to d. 

   If n is zero, return n.
   Otherwise:

      - if the 0-th (first) element of the array is greater than or
        equal to d, return 0.

      - if the (n-1)-th (last) element of the array is less than d,
        return n.

      - otherwise return the lowest index of an element greater than
        or equal to d. 

  The complexity of the algorithm is supposed to be O(log(n)). 

  The algorithm can take g into account, which is supposed to be "a
  good guess" for the sought index for the first element greater than
  or equal to d.

*/
static inline unsigned int __sparsePolynomialFindDegree(constant_t d, 
							constant_t *a, 
							unsigned int n,
							unsigned int g) {
  unsigned int l, r, m, i, j, p;
  unsigned long long int t;
  int tv;

  /* Easy case */
  if (n == 0u) return n;

  /* Check if 0-th (first) element is greater than or equal to d. */
  tv = constantIsGreaterOrEqual(a[0],d,42);
  if (tv == 42) return __sparsePolynomialFindDegreeNaive(d, a, n);
  if (tv) return 0;
  
  /* Check if (n-1)-th (last) element is less than, i.e. not greater
     than nor equal to, d.
  */
  tv = constantIsGreaterOrEqual(a[n-1u],d,1);
  if (tv == 42) return __sparsePolynomialFindDegreeNaive(d, a, n);
  if (!tv) return n;

  /* If g is in the valid range of indices (0 <= g <= n - 1), try to
     find the answer up to 64 steps around g. 
  */
  if (g <= n - 1u) { /* 0 <= g is always true because g is unsigned */
    for (p=n+1u, i=g, j=0; (j<64u) && (i <= n - 1u); j++) {
      tv = constantIsGreaterOrEqual(a[i],d,42);
      if (tv == 42) return __sparsePolynomialFindDegreeNaive(d, a, n);
      if (tv) {
	/* The i-th element of the array is greater than or equal to
	   d. If we just made a step to the right, we know the
	   answer. If i is zero, we stop the search loop without an
	   answer (the guess was incorrect). Otherwise we make a step
	   to the left.
	*/
	if (p == i-1u) return i;
	if (i == 0u) break;
	p = i;
	i--;
      } else {
	/* The i-th element of the array is less than d. If we just
	   made a step to the left, we know the answer. If i is n-1,
	   we stop the search loop without an answer (the guess was
	   incorrect). Otherwise we make a step to the right.
	*/
	if (p == i+1u) return i+1u;
	if (i == n-1u) break;
	p = i;
	i++;
      }
    }
  }

  /* We still didn't find the right index. Perform a bisection.

     We maintain to indices l and r such that we are sure that:

     *    a[l] < d 
     *    d <= a[r] 
     
     In order to be sure that m = floor((l+r)/2) satisfies

     l + 1 <= m <= r - 1, 

     i.e. that the interval becomes strictly smaller in each step, 
     we have to suppose that 

     l <= r - 2,

     i.e. that 

     r - l > 1.

  */
  for (l=0u, r=n-1u; r - l > 1u; /* nothing */) {
    /* Compute the middle of the interval between l and r */
    t = l; t += r; t >>= 1; m = t;
    
    /* Check value at the middle of the interval */
    tv = constantIsGreaterOrEqual(a[m],d,42);
    if (tv == 42) return __sparsePolynomialFindDegreeNaive(d, a, n);
    if (tv) {
      /* Here, d <= a[m], so the new interval becomes [l;m] */
      r = m;
    } else {
      /* Here, a[m] < d, so the new interval becomes [m;r] */
      l = m;
    }
  }

  /* Here, we know that a[l] < d and d <= a[r] */
  for (i=l; l<=r; i++) {
    tv = constantIsGreaterOrEqual(a[i],d,42);
    if (tv == 42) return __sparsePolynomialFindDegreeNaive(d, a, n);
    if (tv) return i;
  }
  
  /* Can only be reached if comparison not possible */
  return __sparsePolynomialFindDegreeNaive(d, a, n);
}

static inline sparse_polynomial_t sparsePolynomialFromMpfrConstant(mpfr_t c) {
  sparse_polynomial_t res;
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->monomialCount = 1;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->coeffs[0] = constantFromMpfr(c);
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees[0] = constantFromInt(0);
  res->deg = constantFromCopy(res->monomialDegrees[0]);
  res->hash.hasHash = 0;
  return res;
}

static inline sparse_polynomial_t sparsePolynomialFromMpzConstant(mpz_t c) {
  sparse_polynomial_t res;
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->monomialCount = 1;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->coeffs[0] = constantFromMpz(c);
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees[0] = constantFromInt(0);
  res->deg = constantFromCopy(res->monomialDegrees[0]);
  res->hash.hasHash = 0;
  return res;
}

static inline sparse_polynomial_t sparsePolynomialFromMpqConstant(mpq_t c) {
  sparse_polynomial_t res;
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->monomialCount = 1;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  mpq_canonicalize(c);
  res->coeffs[0] = constantFromMpq(c);
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees[0] = constantFromInt(0);
  res->deg = constantFromCopy(res->monomialDegrees[0]);
  res->hash.hasHash = 0;
  return res;
}

static inline sparse_polynomial_t sparsePolynomialFromConstant(constant_t c) {
  sparse_polynomial_t res;
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->monomialCount = 1;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->coeffs[0] = constantFromCopy(c);
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees[0] = constantFromInt(0);
  res->deg = constantFromCopy(res->monomialDegrees[0]);
  res->hash.hasHash = 0;
  return res;
}

static inline sparse_polynomial_t sparsePolynomialFromIntConstant(int c) {
  sparse_polynomial_t res;
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->monomialCount = 1;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->coeffs[0] = constantFromInt(c);
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees[0] = constantFromInt(0);
  res->deg = constantFromCopy(res->monomialDegrees[0]);
  res->hash.hasHash = 0;
  return res;
}

static inline int __sparsePolynomialFromConstantExpression(sparse_polynomial_t *r, node *c) {
  sparse_polynomial_t res;
  if (c == NULL) return 0;
  if (!isConstant(c)) return 0;
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->monomialCount = 1;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->coeffs[0] = constantFromExpression(c);
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees[0] = constantFromInt(0);
  res->deg = constantFromCopy(res->monomialDegrees[0]);
  res->hash.hasHash = 0;
  *r = res;
  return 1;
}

static inline int __sparsePolynomialFromConstantExpressionOnlyRealCoeffs(sparse_polynomial_t *r, node *c) {
  sparse_polynomial_t res;
  if (c == NULL) return 0;
  if (!isConstant(c)) return 0;
  if (!containsOnlyRealNumbers(c)) return 0;
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->monomialCount = 1;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->coeffs[0] = constantFromExpression(c);
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees[0] = constantFromInt(0);
  res->deg = constantFromCopy(res->monomialDegrees[0]);
  res->hash.hasHash = 0;
  *r = res;
  return 1;
}

static inline sparse_polynomial_t sparsePolynomialFromIdentity() {
  sparse_polynomial_t res;
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->monomialCount = 1;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->coeffs[0] = constantFromInt(1);
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees[0] = constantFromCopy(res->coeffs[0]);
  res->deg = constantFromCopy(res->monomialDegrees[0]);
  res->hash.hasHash = 0;
  return res;
}

static inline sparse_polynomial_t sparsePolynomialFromMpfrCoefficients(mpfr_t *coeffs, unsigned int deg) {
  unsigned int i, startSize;
  sparse_polynomial_t res;
  constant_t c;
  
  if (coeffs == NULL) return NULL;
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->hash.hasHash = 0;
  startSize = deg + 1u;
  if (startSize == 0u) startSize = UINT_MAX;
  res->monomialCount = 0u;
  res->coeffs = (constant_t *) safeCalloc(startSize, sizeof(constant_t));
  res->monomialDegrees = (constant_t *) safeCalloc(startSize, sizeof(constant_t));
  for (i=0u;i<=deg;i++) {
    c = constantFromMpfr(coeffs[i]);
    if (!constantIsZero(c,0)) {
      res->coeffs[res->monomialCount] = c;
      res->monomialDegrees[res->monomialCount] = constantFromUnsignedInt(i);
      (res->monomialCount)++;
    } else {
      constantFree(c);
    }
  }
  /* If res->monomialCount still is zero, we never added anything
     because the two polynomials multiplied gave zero. We add a zero
     coefficient in this case. 
  */
  if (res->monomialCount == 0u) {
    res->coeffs[res->monomialCount] = constantFromInt(0);
    res->monomialDegrees[res->monomialCount] = constantFromInt(0);
    (res->monomialCount)++;
  }
  /* Set the degree of the polynomial */
  res->deg = constantFromCopy(res->monomialDegrees[res->monomialCount - 1u]);
  /* Adjust memory to the real amount of monomials */
  if (res->monomialCount != startSize) {
    res->coeffs = (constant_t *) safeRealloc(res->coeffs, 
					     ((size_t) (res->monomialCount)) * sizeof(constant_t));
    res->monomialDegrees = (constant_t *) safeRealloc(res->monomialDegrees, 
						      ((size_t) (res->monomialCount)) * sizeof(constant_t));
  }
  /* Adjust the degree for zero leading coefficients */
  __sparsePolynomialAdjustDegree(res);

  /* Return the polynomial */
  return res;
}

static inline int sparsePolynomialFromConstantExpressionCoefficients(sparse_polynomial_t *r, node **coeffs, unsigned int deg) {
  unsigned int i, startSize;
  sparse_polynomial_t res;
  constant_t c;
  
  if (coeffs == NULL) return 0;
  for (i=0;i<=deg;i++) {
    if ((coeffs[i] != NULL) && (!isConstant(coeffs[i]))) return 0;
  }
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->hash.hasHash = 0;
  startSize = deg + 1u;
  if (startSize == 0u) startSize = UINT_MAX;
  res->monomialCount = 0u;
  res->coeffs = (constant_t *) safeCalloc(startSize, sizeof(constant_t));
  res->monomialDegrees = (constant_t *) safeCalloc(startSize, sizeof(constant_t));
  for (i=0u;i<=deg;i++) {
    if (coeffs[i] != NULL) {
      c = constantFromExpression(coeffs[i]);
    } else {
      c = constantFromInt(0);
    }
    if (!constantIsZero(c,0)) {
      res->coeffs[res->monomialCount] = c;
      res->monomialDegrees[res->monomialCount] = constantFromUnsignedInt(i);
      (res->monomialCount)++;
    } else {
      constantFree(c);
    }
  }
  /* If res->monomialCount still is zero, we never added anything
     because the two polynomials multiplied gave zero. We add a zero
     coefficient in this case. 
  */
  if (res->monomialCount == 0u) {
    res->coeffs[res->monomialCount] = constantFromInt(0);
    res->monomialDegrees[res->monomialCount] = constantFromInt(0);
    (res->monomialCount)++;
  }
  /* Set the degree of the polynomial */
  res->deg = constantFromCopy(res->monomialDegrees[res->monomialCount - 1u]);
  /* Adjust memory to the real amount of monomials */
  if (res->monomialCount != startSize) {
    res->coeffs = (constant_t *) safeRealloc(res->coeffs, 
					     ((size_t) (res->monomialCount)) * sizeof(constant_t));
    res->monomialDegrees = (constant_t *) safeRealloc(res->monomialDegrees, 
						      ((size_t) (res->monomialCount)) * sizeof(constant_t));
  }
  /* Adjust the degree for zero leading coefficients */
  __sparsePolynomialAdjustDegree(res);

  /* Return the polynomial */
  *r = res;
  return 1;
}

static inline sparse_polynomial_t sparsePolynomialFromCopy(sparse_polynomial_t p) {
  if (p == NULL) return NULL;
  p->refCount++;
  return p;
}

static inline void sparsePolynomialFree(sparse_polynomial_t p) {
  unsigned int i;

  if (p == NULL) return;
  p->refCount--;
  if (p->refCount > 0u) return;
  constantFree(p->deg);
  for (i=0;i<p->monomialCount;i++) {
    constantFree(p->coeffs[i]);
    constantFree(p->monomialDegrees[i]);
  }
  safeFree(p->coeffs);
  safeFree(p->monomialDegrees);
  __sparsePolynomialFreeMem(p);
}

static inline int sparsePolynomialIsConstant(sparse_polynomial_t p, int defVal) {
  int degZero;

  if (p == NULL) return defVal;
  if (p->monomialCount == 0u) return 1;
  degZero = constantIsZero(p->deg, 42);
  if (degZero == 42) return defVal;
  if (!degZero) return 0;
  return 1;
}

static inline int sparsePolynomialConstantGetConstant(constant_t *c, sparse_polynomial_t p) {
  unsigned int i;
  constant_t t;

  if (p == NULL) return 0;
  if (!sparsePolynomialIsConstant(p, 0)) return 0;
  if (p->monomialCount == 0u) {
    *c = constantFromInt(0);
    return 1;
  }
  if (p->monomialCount == 1u) {
    *c = constantFromCopy(p->coeffs[0]);
    return 1;
  }
  *c = constantFromCopy(p->coeffs[0]);
  for (i=1;i<p->monomialCount;i++) {
    t = constantAdd(*c, p->coeffs[i]);
    constantFree(*c);
    *c = t;
  }
  return 1;
}

static inline int sparsePolynomialIsConstantZero(sparse_polynomial_t p, int defVal) {
  int isConst, res;
  constant_t c;

  if (p == NULL) return defVal;
  isConst = sparsePolynomialIsConstant(p, 42);
  if (isConst == 42) return defVal;
  if (!isConst) return 0;
  if (sparsePolynomialConstantGetConstant(&c, p)) {
    res = constantIsZero(c, defVal);
    constantFree(c);
    return res;
  }
  return defVal;
}

static inline int sparsePolynomialIsConstantOne(sparse_polynomial_t p, int defVal) {
  int isConst, res;
  constant_t c;

  if (p == NULL) return defVal;
  isConst = sparsePolynomialIsConstant(p, 42);
  if (isConst == 42) return defVal;
  if (!isConst) return 0;
  if (sparsePolynomialConstantGetConstant(&c, p)) {
    res = constantIsOne(c, defVal);
    constantFree(c);
    return res;
  }
  return defVal;
}

static inline int sparsePolynomialIsConstantInteger(sparse_polynomial_t p, int k, int defVal) {
  int isConst, res;
  constant_t c, d;

  if (p == NULL) return defVal;
  isConst = sparsePolynomialIsConstant(p, 42);
  if (isConst == 42) return defVal;
  if (!isConst) return 0;
  if (sparsePolynomialConstantGetConstant(&c, p)) {
    d = constantFromInt(k);
    res = constantIsEqual(c, d, defVal);
    constantFree(d);
    constantFree(c);
    return res;
  }
  return defVal;
}

static inline sparse_polynomial_t sparsePolynomialAdd(sparse_polynomial_t p, sparse_polynomial_t q) {
  sparse_polynomial_t res;
  unsigned int i,j,startSize;
  constant_t newCoeff, newMonomial;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;
  if (q == NULL) return NULL;
  
  /* Handle addition of zero polynomial */
  if (sparsePolynomialIsConstantZero(p, 0)) {
    return sparsePolynomialFromCopy(q);
  }
  if (sparsePolynomialIsConstantZero(q, 0)) {
    return sparsePolynomialFromCopy(p);
  }

  /* General case */
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->hash.hasHash = 0;
  startSize = p->monomialCount + q->monomialCount;
  if (startSize < p->monomialCount) {
    startSize = UINT_MAX;
  }
  res->coeffs = (constant_t *) safeCalloc(startSize, 
					  sizeof(constant_t));
  res->monomialDegrees = (constant_t *) safeCalloc(startSize, 
						   sizeof(constant_t));
  res->monomialCount = 0;
  for (i=0u, j=0u; (i<p->monomialCount) && (j<q->monomialCount); /* nothing */) {
    /* We have 3 cases:

       (i)   The degree of the i-th monomial of p and the one of the
             j-th monomial of q are the same.

	     In this case we add the coefficients and increment both i and j.
	   
       (ii)  The degree of the i-th monomial of p is greater than the one 
             of the j-th monomial of q.

	     In this case we copy the j-th coefficient of q and increment j.

       (iii) The degree of the i-th monomial of p is less than the one 
             of the j-th monomial of q. 

	     In this case we copy the i-th coefficient of p and increment i.
	    
       In case (i) it is possible that after addition, the coefficient
       gets zero. In this case we not add the monomial to the
       polynomial.

    */
    if (constantIsEqual(p->monomialDegrees[i], q->monomialDegrees[j], 0)) {
      /* Case (i) */
      newCoeff = constantAdd(p->coeffs[i], q->coeffs[j]);
      newMonomial = constantFromCopy(p->monomialDegrees[i]);
      i++;
      j++;
    } else {
      if (constantIsGreater(p->monomialDegrees[i], q->monomialDegrees[j], 0)) {
	/* Case (iii) */
	newCoeff = constantFromCopy(q->coeffs[j]);
	newMonomial = constantFromCopy(q->monomialDegrees[j]);
	j++;
      } else {
	/* Case (ii) */
	newCoeff = constantFromCopy(p->coeffs[i]);
	newMonomial = constantFromCopy(p->monomialDegrees[i]);
	i++;
      }
    }
    /* Now check if we have to add the monomial */
    if (!constantIsZero(newCoeff, 0)) {
      /* Add the coefficient */
      res->coeffs[res->monomialCount] = newCoeff;
      res->monomialDegrees[res->monomialCount] = newMonomial;
      (res->monomialCount)++;
    } else {
      /* Just free the new coefficient */
      constantFree(newCoeff);
      constantFree(newMonomial);
    }
  }
  /* Now add the left-over monomials of p and q */
  for (/* nothing */; i<p->monomialCount; i++) {
    res->coeffs[res->monomialCount] = constantFromCopy(p->coeffs[i]);
    res->monomialDegrees[res->monomialCount] = constantFromCopy(p->monomialDegrees[i]);
    (res->monomialCount)++;
  }
  for (/* nothing */; j<q->monomialCount; j++) {
    res->coeffs[res->monomialCount] = constantFromCopy(q->coeffs[j]);
    res->monomialDegrees[res->monomialCount] = constantFromCopy(q->monomialDegrees[j]);
    (res->monomialCount)++;
  }
  /* If res->monomialCount still is zero, we never added anything
     because the two polynomials cancelled out. We add a zero
     coefficient in this case. 
  */
  if (res->monomialCount == 0u) {
    res->coeffs[res->monomialCount] = constantFromInt(0);
    res->monomialDegrees[res->monomialCount] = constantFromInt(0);
    (res->monomialCount)++;
  }
  /* Set the degree of the polynomial */
  res->deg = constantFromCopy(res->monomialDegrees[res->monomialCount - 1u]);
  /* Adjust memory to the real amount of monomials */
  if (res->monomialCount != startSize) {
    res->coeffs = (constant_t *) safeRealloc(res->coeffs, 
					     ((size_t) (res->monomialCount)) * sizeof(constant_t));
    res->monomialDegrees = (constant_t *) safeRealloc(res->monomialDegrees, 
						      ((size_t) (res->monomialCount)) * sizeof(constant_t));
  }
  /* Adjust the degree for zero leading coefficients */
  __sparsePolynomialAdjustDegree(res);

  return res;
}

static inline sparse_polynomial_t sparsePolynomialNeg(sparse_polynomial_t);

static inline sparse_polynomial_t sparsePolynomialSub(sparse_polynomial_t p, sparse_polynomial_t q) {
  sparse_polynomial_t res;
  unsigned int i,j,startSize;
  constant_t newCoeff, newMonomial;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;
  if (q == NULL) return NULL;
  
  /* Handle subtraction of zero polynomial */
  if (sparsePolynomialIsConstantZero(p, 0)) {
    return sparsePolynomialNeg(q);
  }
  if (sparsePolynomialIsConstantZero(q, 0)) {
    return sparsePolynomialFromCopy(p);
  }

  /* General case */
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->hash.hasHash = 0;
  startSize = p->monomialCount + q->monomialCount;
  if (startSize < p->monomialCount) {
    startSize = UINT_MAX;
  }
  res->coeffs = (constant_t *) safeCalloc(startSize, 
					  sizeof(constant_t));
  res->monomialDegrees = (constant_t *) safeCalloc(startSize, 
						   sizeof(constant_t));
  res->monomialCount = 0;
  for (i=0u, j=0u; (i<p->monomialCount) && (j<q->monomialCount); /* nothing */) {
    /* We have 3 cases:

       (i)   The degree of the i-th monomial of p and the one of the
             j-th monomial of q are the same.

	     In this case we subtract the coefficients and increment
	     both i and j.
	   
       (ii)  The degree of the i-th monomial of p is greater than the one 
             of the j-th monomial of q.

	     In this case we negate the j-th coefficient of q and increment j.

       (iii) The degree of the i-th monomial of p is less than the one 
             of the j-th monomial of q. 

	     In this case we copy the i-th coefficient of p and increment i.
	    
       In case (i) it is possible that after addition, the coefficient
       gets zero. In this case we not add the monomial to the
       polynomial.

    */
    if (constantIsEqual(p->monomialDegrees[i], q->monomialDegrees[j], 0)) {
      /* Case (i) */
      newCoeff = constantSub(p->coeffs[i], q->coeffs[j]);
      newMonomial = constantFromCopy(p->monomialDegrees[i]);
      i++;
      j++;
    } else {
      if (constantIsGreater(p->monomialDegrees[i], q->monomialDegrees[j], 0)) {
	/* Case (ii) */
	newCoeff = constantNeg(q->coeffs[j]);
	newMonomial = constantFromCopy(q->monomialDegrees[j]);
	j++;
      } else {
	/* Case (iii) */
	newCoeff = constantFromCopy(p->coeffs[i]);
	newMonomial = constantFromCopy(p->monomialDegrees[i]);
	i++;
      }
    }
    /* Now check if we have to add the monomial */
    if (!constantIsZero(newCoeff, 0)) {
      /* Add the coefficient */
      res->coeffs[res->monomialCount] = newCoeff;
      res->monomialDegrees[res->monomialCount] = newMonomial;
      (res->monomialCount)++;
    } else {
      /* Just free the new coefficient */
      constantFree(newCoeff);
      constantFree(newMonomial);
    }
  }
  /* Now add the left-over monomials of p and q */
  for (/* nothing */; i<p->monomialCount; i++) {
    res->coeffs[res->monomialCount] = constantFromCopy(p->coeffs[i]);
    res->monomialDegrees[res->monomialCount] = constantFromCopy(p->monomialDegrees[i]);
    (res->monomialCount)++;
  }
  for (/* nothing */; j<q->monomialCount; j++) {
    res->coeffs[res->monomialCount] = constantNeg(q->coeffs[j]);
    res->monomialDegrees[res->monomialCount] = constantFromCopy(q->monomialDegrees[j]);
    (res->monomialCount)++;
  }
  /* If res->monomialCount still is zero, we never added anything
     because the two polynomials cancelled out. We add a zero
     coefficient in this case. 
  */
  if (res->monomialCount == 0u) {
    res->coeffs[res->monomialCount] = constantFromInt(0);
    res->monomialDegrees[res->monomialCount] = constantFromInt(0);
    (res->monomialCount)++;
  }
  /* Set the degree of the polynomial */
  res->deg = constantFromCopy(res->monomialDegrees[res->monomialCount - 1u]);
  /* Adjust memory to the real amount of monomials */
  if (res->monomialCount != startSize) {
    res->coeffs = (constant_t *) safeRealloc(res->coeffs, 
					     ((size_t) (res->monomialCount)) * sizeof(constant_t));
    res->monomialDegrees = (constant_t *) safeRealloc(res->monomialDegrees, 
						      ((size_t) (res->monomialCount)) * sizeof(constant_t));
  }
  /* Adjust the degree for zero leading coefficients */
  __sparsePolynomialAdjustDegree(res);

  return res;
}

static inline int sparsePolynomialEqual(sparse_polynomial_t p, sparse_polynomial_t q, int defVal) {
  unsigned int i;
  int monDegEqual, coeffEqual, degreeEqual;
  int res;
  sparse_polynomial_t d;

  if (p == NULL) return defVal;
  if (q == NULL) return defVal;
  if (p == q) return 1;
  if (p->hash.hasHash && q->hash.hasHash) {
    if (p->hash.hash != q->hash.hash) return 0;
  }
  if (p->monomialCount == 0u) {
    return sparsePolynomialIsConstantZero(q, defVal);
  }
  if (q->monomialCount == 0u) {
    return sparsePolynomialIsConstantZero(p, defVal);
  }
  if (p->monomialCount != q->monomialCount) {
    if (__sparsePolynomialAscendingDegrees(p, 0) &&
	__sparsePolynomialAscendingDegrees(q, 0))
      return 0;
    d = sparsePolynomialSub(p, q);
    res = sparsePolynomialIsConstantZero(d, defVal);
    sparsePolynomialFree(d);
    return res;
  }
  degreeEqual = constantIsEqual(p->deg, q->deg, 42);
  if (degreeEqual == 42) return defVal;
  if (!degreeEqual) return 0;  
  for (i=0;i<p->monomialCount;i++) {
    monDegEqual = constantIsEqual(p->monomialDegrees[i], q->monomialDegrees[i], 42);
    if (monDegEqual == 42) return defVal;
    if (!monDegEqual) return 0;
    coeffEqual = constantIsEqual(p->coeffs[i], q->coeffs[i], 42);
    if (coeffEqual == 42) return defVal;
    if (!coeffEqual) return 0;
  }

  return 1;
}

static inline int sparsePolynomialIsIdentity(sparse_polynomial_t p, int defVal) {
  sparse_polynomial_t q;
  int res;

  if (p == NULL) return defVal;
  if (sparsePolynomialGetDegreeAsInt(p) != 1) return 0;

  /* Can be optimized */
  q = sparsePolynomialFromIdentity();
  res = sparsePolynomialEqual(p, q, defVal);
  sparsePolynomialFree(q);

  return res;
}

static inline void __sparsePolynomialAddMonomialToArrays(constant_t *coeffs, 
						   constant_t *degrees,
						   unsigned int *n,
						   unsigned int *g,
						   constant_t c,
						   constant_t d) {
  unsigned int i, j;
  constant_t t;

  /* Find the index where to add the coefficient */
  j = __sparsePolynomialFindDegree(d, degrees, *n, *g);

  /* Keep the value of the index as a guess for the next time */
  *g = j;

  /* If the coefficient's degree is greater than all the others, it
     goes at the end.
  */
  if (j == *n) {
    /* Adding the coefficient at the end */
    coeffs[j] = c;
    degrees[j] = d;
    (*n)++;
    return;
  }

  /* Here, we know that 0 <= j <= n - 1 

     There are two cases: 

     a) degree[j] == d. In this case we add c to degree[j]. 
        
        We have two subcases: 

	(i)  After addition the new coefficient is non-zero. We leave
	     it untouched.

        (ii) After addition the new coefficient is zero. We remove it by
	     moving all following coefficients one to the left.

     b) degree[j] > d. In this case we deplace all existing
        coefficients and degrees of indices j thru n-1 and put the new
	coefficient and degree at index j.

  */
  if (constantIsEqual(degrees[j],d,0)) {
    /* Case a) */
    t = constantAdd(coeffs[j], c);
    constantFree(coeffs[j]);
    if (constantIsZero(t, 0)) {
      constantFree(t);
      constantFree(degrees[j]);
      for (i=j;i<*n;i++) {
	coeffs[i] = coeffs[i+1u];
	degrees[i] = degrees[i+1u];
      }
      (*n)--;
    } else {
      coeffs[j] = t;
    }
    constantFree(c);
    constantFree(d);
    return;
  }

  /* Case b) 

     Deplace the existing coefficients and degrees.
  */
  for (i=*n;i>j;i--) {
    coeffs[i] = coeffs[i-1u];
    degrees[i] = degrees[i-1u];
  }

  /* Account for the new entry */
  (*n)++;

  /* Put the new entry */
  coeffs[j] = c;
  degrees[j] = d;
}

static inline sparse_polynomial_t sparsePolynomialAddConstant(sparse_polynomial_t p, constant_t c) {
  sparse_polynomial_t t, res;

  /* Handle the stupid cases */
  if (p == NULL) return NULL;
  if (c == NULL) return NULL;

  /* Can be optimized */
  t = sparsePolynomialFromConstant(c);
  res = sparsePolynomialAdd(p, t);
  sparsePolynomialFree(t);

  return res;
}

static inline sparse_polynomial_t sparsePolynomialMul(sparse_polynomial_t p, sparse_polynomial_t q) {
  sparse_polynomial_t res;
  unsigned int i,j,startSize,g,d;
  constant_t pp, md, t;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;
  if (q == NULL) return NULL;
  
  /* Handle multiplication of zero polynomial */
  if (sparsePolynomialIsConstantZero(p, 0)) {
    return sparsePolynomialFromCopy(p);
  }
  if (sparsePolynomialIsConstantZero(q, 0)) {
    return sparsePolynomialFromCopy(q);
  }

  /* Handle multiplication of polynomial one */
  if (sparsePolynomialIsConstantOne(p, 0)) {
    return sparsePolynomialFromCopy(q);
  }
  if (sparsePolynomialIsConstantOne(q, 0)) {
    return sparsePolynomialFromCopy(p);
  }

  /* General case 

     Start by determining the maximum number of coefficients. This number is 
     upper bounded by

     min(p->monomialCount * q->monomialCount, deg(p) + deg(q) + 1)

  */
  if (!tryExactUnsignedIntMultiplication(&startSize, 
					 p->monomialCount, q->monomialCount)) {
    startSize = UINT_MAX;
  }
  t = constantAdd(p->deg, q->deg);
  if (tryConstantToUnsignedInt(&d, t)) {
    d++;
    if (d != 0u) {
      if (d < startSize) startSize = d;
    }
  }
  constantFree(t);

  /* Allocate the result polynomial */
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->hash.hasHash = 0;
  res->coeffs = (constant_t *) safeCalloc(startSize, 
					  sizeof(constant_t));
  res->monomialDegrees = (constant_t *) safeCalloc(startSize, 
						   sizeof(constant_t));
  res->monomialCount = 0;
  
  /* Produce all partial products p_i * q_i * x^(m_i + n_i) */
  g = 0;
  for (i=0u;i<p->monomialCount;i++) {
    for (j=0u;j<q->monomialCount;j++) {
      pp = constantMul(p->coeffs[i], q->coeffs[j]);
      if (!constantIsZero(pp, 0)) {
	md = constantAdd(p->monomialDegrees[i],q->monomialDegrees[j]);
	__sparsePolynomialAddMonomialToArrays(res->coeffs, res->monomialDegrees, 
					&(res->monomialCount), &g, pp, md);
      } else {
	constantFree(pp);
      }
    }
  }

  /* If res->monomialCount still is zero, we never added anything
     because the two polynomials multiplied gave zero. We add a zero
     coefficient in this case. 
  */
  if (res->monomialCount == 0u) {
    res->coeffs[res->monomialCount] = constantFromInt(0);
    res->monomialDegrees[res->monomialCount] = constantFromInt(0);
    (res->monomialCount)++;
  }
  /* Set the degree of the polynomial */
  res->deg = constantFromCopy(res->monomialDegrees[res->monomialCount - 1u]);
  /* Adjust memory to the real amount of monomials */
  if (res->monomialCount != startSize) {
    res->coeffs = (constant_t *) safeRealloc(res->coeffs, 
					     ((size_t) (res->monomialCount)) * sizeof(constant_t));
    res->monomialDegrees = (constant_t *) safeRealloc(res->monomialDegrees, 
						      ((size_t) (res->monomialCount)) * sizeof(constant_t));
  }
  /* Adjust the degree for zero leading coefficients */
  __sparsePolynomialAdjustDegree(res);

  return res;
}

static inline sparse_polynomial_t sparsePolynomialNeg(sparse_polynomial_t p) {
  unsigned int i;
  sparse_polynomial_t res;

  if (p == NULL) return NULL;
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->hash.hasHash = 0;
  res->deg = constantFromCopy(p->deg);
  res->monomialCount = p->monomialCount;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  for (i=0u;i<res->monomialCount;i++) {
    res->coeffs[i] = constantNeg(p->coeffs[i]);
    res->monomialDegrees[i] = constantFromCopy(p->monomialDegrees[i]);
  }
  __sparsePolynomialAdjustDegree(res);
  return res;
}

static inline sparse_polynomial_t sparsePolynomialCompose(sparse_polynomial_t p, sparse_polynomial_t q) {
  sparse_polynomial_t res, t, m;
  constant_t d;
  unsigned int i;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;
  if (q == NULL) return NULL;

  /* Handle the strange case when p has no monomials */
  if (p->monomialCount == 0u) return sparsePolynomialFromIntConstant(0);

  /* Perform Horner evaluation of p in q */
  res = sparsePolynomialFromConstant(p->coeffs[p->monomialCount-1u]);
  for (i=p->monomialCount-1u;i>=1u;i--) {
    d = constantSub(p->monomialDegrees[i], p->monomialDegrees[i-1u]);
    if (!sparsePolynomialPowConstant(&m, q, d)) {
      sollyaFprintf(stderr,"Error: sparsePolynomialCompose: monomial degrees not appropriately ordered\n");
      exit(1);
    }
    t = sparsePolynomialMul(res, m);
    sparsePolynomialFree(res);
    res = t;
    sparsePolynomialFree(m);
    constantFree(d);
    t = sparsePolynomialAddConstant(res, p->coeffs[i-1u]);
    sparsePolynomialFree(res);
    res = t;
  }
  if (!sparsePolynomialPowConstant(&m, q, p->monomialDegrees[0u])) {
    sollyaFprintf(stderr,"Error: sparsePolynomialCompose: monomial degrees not appropriately ordered\n");
    exit(1);
  }
  t = sparsePolynomialMul(res, m);
  sparsePolynomialFree(res);
  res = t;
  sparsePolynomialFree(m);

  return res;
}

 static inline void __sparsePolynomialGetLeadingCoefficient(constant_t *c, constant_t *d, sparse_polynomial_t *r, sparse_polynomial_t p) {
  sparse_polynomial_t res;
  unsigned int i;

  /* Handle stupid cases */
  if (p == NULL) {
    *c = NULL;
    *d = NULL;
    *r = NULL;
    return;
  }

  /* Handle a strange case that should never happen */
  if (p->monomialCount == 0u) {
    *c = constantFromInt(0);
    *d = constantFromInt(0);
    *r = sparsePolynomialFromIntConstant(0);
    return;
  }

  /* Handle the case where the polynomial has just one monomial */
  if (p->monomialCount == 1u) {
    *c = constantFromCopy(p->coeffs[0]);
    *d = constantFromCopy(p->monomialDegrees[0]);
    *r = sparsePolynomialFromIntConstant(0);
    return;
  }
  
  /* Here, the polynomial has at least two monomials */
  *c = constantFromCopy(p->coeffs[p->monomialCount - 1u]);
  *d = constantFromCopy(p->monomialDegrees[p->monomialCount - 1u]);
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->hash.hasHash = 0;
  res->monomialCount = p->monomialCount - 1u;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  for (i=0u;i<res->monomialCount;i++) {
    res->coeffs[i] = constantFromCopy(p->coeffs[i]);
    res->monomialDegrees[i] = constantFromCopy(p->monomialDegrees[i]);
  }
  res->deg = constantFromCopy(res->monomialDegrees[res->monomialCount-1u]);
  __sparsePolynomialAdjustDegree(res);
  *r = res;
}

static inline sparse_polynomial_t __sparsePolynomialFromMonomial(constant_t c, constant_t d) {
  sparse_polynomial_t res;
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->hash.hasHash = 0;
  res->monomialCount = 1;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->coeffs[0] = constantFromCopy(c);
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees[0] = constantFromCopy(d);
  res->deg = constantFromCopy(d);
  return res;
}

static inline void sparsePolynomialDiv(sparse_polynomial_t *quot, sparse_polynomial_t *rest, sparse_polynomial_t a, sparse_polynomial_t b) {
  sparse_polynomial_t q, r, rt, bt, t1, t2, qp, recprPB;
  constant_t rc, rd, bc, bd, qc, qd, recprB, one;
  

  /* Handle stupid cases */
  if ((a == NULL) || (b == NULL)) {
    *quot = NULL;
    *rest = NULL;
    return;
  }

  /* If both polynomials are equal, the quotient is one, the rest is
     zero.
  */
  if (sparsePolynomialEqual(a, b, 0)) {
    *quot = sparsePolynomialFromIntConstant(1);
    *rest = sparsePolynomialFromIntConstant(0);
    return;
  }
  
  /* If b is constant, do some special handling */
  if (sparsePolynomialConstantGetConstant(&bc, b)) {
    if (constantIsOne(bc, 0)) {
      *quot = sparsePolynomialFromCopy(a);
      *rest = sparsePolynomialFromIntConstant(0);
    } else {
      if (!constantIsZero(bc, 1)) {
	one = constantFromInt(1);
	recprB = constantDiv(one, bc);
	constantFree(one);
	recprPB = sparsePolynomialFromConstant(recprB);
	constantFree(recprB);
	*quot = sparsePolynomialMul(a, recprPB);
	*rest = sparsePolynomialFromIntConstant(0);
	sparsePolynomialFree(recprPB);
      } else {
	*quot = sparsePolynomialFromIntConstant(0);
	*rest = sparsePolynomialFromCopy(a);
      }
    }
    constantFree(bc);
    return;
  }

  /* Initialize q and r to

     q = 0
     r = a

     such that 

     a = q * b + r

  */
  q = sparsePolynomialFromIntConstant(0);
  r = sparsePolynomialFromCopy(a);
    
  /* Euclidean division */
  __sparsePolynomialGetLeadingCoefficient(&bc,&bd,&bt,b);
  while ((!constantIsGreater(b->deg,r->deg,1)) && 
	 (!sparsePolynomialIsConstantZero(r,1))) {
    __sparsePolynomialGetLeadingCoefficient(&rc,&rd,&rt,r);
    qc = constantDiv(rc,bc);
    constantFree(rc);
    qd = constantSub(rd,bd);
    constantFree(rd);
    qp = __sparsePolynomialFromMonomial(qc,qd);
    constantFree(qc);
    constantFree(qd);
    t1 = sparsePolynomialAdd(q, qp);
    sparsePolynomialFree(q);
    q = t1;
    t1 = sparsePolynomialMul(qp,bt);
    t2 = sparsePolynomialSub(rt,t1);
    sparsePolynomialFree(t1);
    sparsePolynomialFree(rt);
    sparsePolynomialFree(qp);
    sparsePolynomialFree(r);
    r = t2;
  }
  constantFree(bc);
  constantFree(bd);
  sparsePolynomialFree(bt);

  /* Set the result polynomials */
  *quot = q;
  *rest = r;
}

static inline void __sparsePolynomialCutIntoHalves(sparse_polynomial_t *r, sparse_polynomial_t *q, sparse_polynomial_t p) {
  unsigned int i, monomialCount1, monomialCount2;
  sparse_polynomial_t res;

  /* Handle stupid input */
  if (p == NULL) {
    *r = NULL;
    *q = NULL;
    return;
  }

  /* Handle the case when p has no monomials */
  if (p->monomialCount == 0u) {
    *r = sparsePolynomialFromIntConstant(0);
    *q = sparsePolynomialFromIntConstant(0);    
    return;
  }

  /* Handle the case when p has only one monomial */
  if (p->monomialCount == 1u) {
    *r = sparsePolynomialFromCopy(p);
    *q = sparsePolynomialFromIntConstant(0);
    return;
  }

  /* Now, we know that we have at least 2 monomials. This implies that
     both polynomials have at least one monomial. 
  */
  monomialCount1 = p->monomialCount >> 1;
  monomialCount2 = p->monomialCount - monomialCount1;

  /* Form polynomial p */
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->hash.hasHash = 0;
  res->monomialCount = monomialCount1;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  for (i=0u;i<monomialCount1;i++) {
    res->coeffs[i] = constantFromCopy(p->coeffs[i]);
    res->monomialDegrees[i] = constantFromCopy(p->monomialDegrees[i]);
  }
  res->deg = constantFromCopy(res->monomialDegrees[res->monomialCount-1u]);
  __sparsePolynomialAdjustDegree(res);
  *r = res;

  /* Form polynomial q */
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->hash.hasHash = 0;
  res->monomialCount = monomialCount2;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  for (i=0u;i<monomialCount2;i++) {
    res->coeffs[i] = constantFromCopy(p->coeffs[i + monomialCount1]);
    res->monomialDegrees[i] = constantFromCopy(p->monomialDegrees[i + monomialCount1]);
  }
  res->deg = constantFromCopy(res->monomialDegrees[res->monomialCount-1u]);
  __sparsePolynomialAdjustDegree(res);
  *q = res;
}

/* If p is of the form p(x) = x^m * (c + d * x^n) then decompose p into
   c, m, d and n and return a non-zero value.
   
   Otherwise return 0 and don't touch the pointers c, m, d and n.
*/
static inline int sparsePolynomialDecomposeTwoMonomials(constant_t *c, constant_t *m, constant_t *d, constant_t *n,
							sparse_polynomial_t p) {

  /* Handle stupid input */
  if (p == NULL) return 0;

  /* Check if there are really just two monomials */
  if (p->monomialCount != 2u) return 0;

  /* Here we are sure that we only have two monomials.

     Do the decomposition work.

  */
  *c = constantFromCopy(p->coeffs[0]);
  *d = constantFromCopy(p->coeffs[1]);
  *m = constantFromCopy(p->monomialDegrees[0]);
  *n = constantSub(p->monomialDegrees[1],p->monomialDegrees[0]);

  /* Signal success */
  return 1;
}

static inline sparse_polynomial_t __sparsePolynomialPowUnsignedIntAlternate(sparse_polynomial_t p, unsigned int n) {
  constant_t nC, bin;
  unsigned int i;
  sparse_polynomial_t r, q, rToTheI, qToTheNMinusI, binPoly, tempPoly1, tempPoly2, tempPoly3, res;
  
  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* Handle easy cases */
  if (n == 0u) return sparsePolynomialFromIntConstant(1);
  if (n == 1u) return sparsePolynomialFromCopy(p);
  if (n == 2u) return sparsePolynomialMul(p,p);
  
  /* Handle the case when the polynomial only has one monomial */
  if (p->monomialCount == 1u) {
    nC = constantFromUnsignedInt(n);
    res = __sparsePolynomialAllocate();
    res->refCount = 1;
    res->hash.hasHash = 0;
    res->monomialCount = 1;
    res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
    res->coeffs[0] = constantPow(p->coeffs[0],nC);
    res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
    res->monomialDegrees[0] = constantMul(p->monomialDegrees[0],nC);
    res->deg = constantFromCopy(res->monomialDegrees[0]);
    constantFree(nC);
    return res;
  }

  /* The polynomial has at least two monomials 

     Write p(x) = r(x) + q(x), with r and q having approximately half
     of the monomials,

     and obtain

     (p(x))^n = sum_i=0^n bin(n,i) * (r(x))^i * (q(x))^(n - i)
              = r(x)^n + sum_i=0^(n-1) bin(n,i) * (r(x))^i * (q(x))^(n - i)
  */
  __sparsePolynomialCutIntoHalves(&r, &q, p);
  res = sparsePolynomialPowUnsignedInt(r, n);
  for (i=0u;i<n;i++) {
    bin = constantFromBinomialUnsignedInt(n, i);
    binPoly = sparsePolynomialFromConstant(bin);
    constantFree(bin);
    rToTheI = sparsePolynomialPowUnsignedInt(r, i);
    qToTheNMinusI = sparsePolynomialPowUnsignedInt(q, n - i);
    tempPoly1 = sparsePolynomialMul(binPoly, rToTheI);
    tempPoly2 = sparsePolynomialMul(tempPoly1, qToTheNMinusI);
    tempPoly3 = sparsePolynomialAdd(res, tempPoly2);
    sparsePolynomialFree(binPoly);
    sparsePolynomialFree(rToTheI);
    sparsePolynomialFree(qToTheNMinusI);
    sparsePolynomialFree(tempPoly1);
    sparsePolynomialFree(tempPoly2);
    sparsePolynomialFree(res);
    res = tempPoly3;
  }
  sparsePolynomialFree(r);
  sparsePolynomialFree(q);

  return res;
}

static inline sparse_polynomial_t sparsePolynomialPowUnsignedInt(sparse_polynomial_t p, unsigned int n) {
  sparse_polynomial_t res, t, tmp;
  constant_t nC;
  
  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* Handle easy cases */
  if (n == 0u) return sparsePolynomialFromIntConstant(1);
  if (n == 1u) return sparsePolynomialFromCopy(p);
  if (n == 2u) return sparsePolynomialMul(p,p);
  
  /* Handle the case when the polynomial only has one monomial */
  if (p->monomialCount == 1u) {
    nC = constantFromUnsignedInt(n);
    res = __sparsePolynomialAllocate();
    res->refCount = 1;
    res->hash.hasHash = 0;
    res->monomialCount = 1;
    res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
    res->coeffs[0] = constantPow(p->coeffs[0],nC);
    res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
    res->monomialDegrees[0] = constantMul(p->monomialDegrees[0],nC);
    res->deg = constantFromCopy(res->monomialDegrees[0]);
    constantFree(nC);
    return res;
  }

  /* Use alternate algorithm for polynomials that might have
     non-rational coefficients 
  */
  if (!sparsePolynomialCoefficientsAreRational(p, 0)) {
    return __sparsePolynomialPowUnsignedIntAlternate(p, n);
  }

  /* Square and multiply loop */
  t = sparsePolynomialFromCopy(p);
  res = sparsePolynomialFromIntConstant(1);
  while (n > 0u) {
    if (n & 1u) {
      tmp = sparsePolynomialMul(res,t);
      sparsePolynomialFree(res);
      res = tmp;
    }
    n >>= 1;
    if (n > 0u) {
      tmp = sparsePolynomialMul(t,t);
      sparsePolynomialFree(t);
      t = tmp;
    }
  }
  sparsePolynomialFree(t);

  return res;
}

static inline int sparsePolynomialPowConstant(sparse_polynomial_t *r, sparse_polynomial_t p, constant_t n) {
  unsigned int nI;
  constant_t k, l, m;
  sparse_polynomial_t res, a, b, t;

  /* Make compiler happy */
  k = NULL;
  l = NULL;
  /* End of compiler happiness */

  if (p == NULL) return 0;
  if (n == NULL) return 0;
  if (!constantIsNonNegativeInteger(n,0)) return 0;
 
  /* Handle the cases n = 0 and n = 1 */
  if (constantIsZero(n,0)) {
    *r = sparsePolynomialFromIntConstant(1);
    return 1;
  }
  if (constantIsOne(n,0)) {
    *r = sparsePolynomialFromCopy(p);
    return 1;
  }

  /* Handle the case when p only has one monomial */
  if (p->monomialCount == 1u) {
    res = __sparsePolynomialAllocate();
    res->refCount = 1;
    res->hash.hasHash = 0;
    res->monomialCount = 1;
    res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
    res->coeffs[0] = constantPow(p->coeffs[0],n);
    res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
    res->monomialDegrees[0] = constantMul(p->monomialDegrees[0],n);
    res->deg = constantFromCopy(res->monomialDegrees[0]);
    *r = res;
    return 1;
  }
 
  /* Try to use the powering function with an unsigned int argument */
  if (tryConstantToUnsignedInt(&nI, n)) {
    *r = sparsePolynomialPowUnsignedInt(p, nI);
    return 1;
  }

  /* If we are here, n does not hold on an unsigned int 

     Chances are low that we actually have enough memory but we still
     try to handle the case. 

     We rewrite p^n as

     p^n = p^(2^14 * k + l) = (p^k)^(2^14) * p^l

     In order to compute t = p^k, t^(2^14) and p^l, we call ourselves
     recursively. As we know that l is integer and l <= 2^14-1 <
     2^16-1 and 2^14 < 2^16 - 1, we know that for these exponents, we
     will be able to use the branch where the exponent holds on an
     unsigned int.

  */
  constantCutTwo14(&k, &l, n);
  m = constantFromInt(1 << 14);
  
  if (!sparsePolynomialPowConstant(&t,p,k)) {
    constantFree(k);
    constantFree(l);
    constantFree(m);
    return 0;
  }
  constantFree(k);  

  if (!sparsePolynomialPowConstant(&a,t,m)) {
    sparsePolynomialFree(t);
    constantFree(l);
    constantFree(m);
    return 0;
  }
  sparsePolynomialFree(t);
  constantFree(m);

  if (!sparsePolynomialPowConstant(&b,p,l)) {
    sparsePolynomialFree(a);
    constantFree(l);
    return 0;
  }
  constantFree(l);

  *r = sparsePolynomialMul(a,b);

  sparsePolynomialFree(b);
  sparsePolynomialFree(a);
  return 1;
}

static inline int sparsePolynomialPow(sparse_polynomial_t *r, sparse_polynomial_t p, sparse_polynomial_t q) {
  constant_t n;
  int res;

  if (p == NULL) return 0;
  if (q == NULL) return 0;
  if (!sparsePolynomialConstantGetConstant(&n, q)) return 0;
  if (!constantIsNonNegativeInteger(n,0)) {
    constantFree(n);
    return 0;
  }
  res = sparsePolynomialPowConstant(r, p, n);
  constantFree(n);
  return res;
}

static inline sparse_polynomial_t sparsePolynomialDeriv(sparse_polynomial_t p) {
  unsigned int i;
  sparse_polynomial_t res;
  constant_t one;

  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* Handle the strange case when p has no monomials */
  if (p->monomialCount == 0u) {
    return sparsePolynomialFromIntConstant(0);
  }

  /* Handle the case when p is a constant polynomial */
  if (sparsePolynomialIsConstant(p, 0)) {
    return sparsePolynomialFromIntConstant(0);
  }

  /* Handle the general case */
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->hash.hasHash = 0;
  res->coeffs = (constant_t *) safeCalloc(p->monomialCount, sizeof(constant_t));
  res->monomialDegrees = (constant_t *) safeCalloc(p->monomialCount, sizeof(constant_t));
  one = constantFromInt(1);
  res->monomialCount = 0u;
  for (i=0u;i<p->monomialCount;i++) {
    if (constantIsGreaterOrEqual(p->monomialDegrees[i], one, 1)) {
      res->coeffs[res->monomialCount] = constantMul(p->coeffs[i],p->monomialDegrees[i]);
      res->monomialDegrees[res->monomialCount] = constantSub(p->monomialDegrees[i],one);
      (res->monomialCount)++;
    }
  }
  constantFree(one);

  /* If res->monomialCount still is zero, we never added anything. We
     add a zero coefficient in this case.
  */
  if (res->monomialCount == 0u) {
    res->coeffs[res->monomialCount] = constantFromInt(0);
    res->monomialDegrees[res->monomialCount] = constantFromInt(0);
    (res->monomialCount)++;
  }
  /* Set the degree of the polynomial */
  res->deg = constantFromCopy(res->monomialDegrees[res->monomialCount - 1u]);
  /* Adjust memory to the real amount of monomials */
  if (res->monomialCount != p->monomialCount) {
    res->coeffs = (constant_t *) safeRealloc(res->coeffs, 
					     ((size_t) (res->monomialCount)) * sizeof(constant_t));
    res->monomialDegrees = (constant_t *) safeRealloc(res->monomialDegrees, 
						      ((size_t) (res->monomialCount)) * sizeof(constant_t));
  }
  /* Adjust the degree for zero leading coefficients */
  __sparsePolynomialAdjustDegree(res);

  return res;
}

static inline sparse_polynomial_t __sparsePolynomialGcdBaseCase(sparse_polynomial_t p, sparse_polynomial_t q) {
  sparse_polynomial_t u, v, t, z;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;
  if (q == NULL) return NULL;

  /* Euclidian algorithm */
  u = sparsePolynomialFromCopy(p);
  v = sparsePolynomialFromCopy(q);
  while (!sparsePolynomialIsConstantZero(v, 1)) {
    sparsePolynomialDiv(&z, &t, u, v);
    sparsePolynomialFree(z);
    sparsePolynomialFree(u);
    u = v;
    v = t;
  }
  if (!sparsePolynomialIsConstantZero(v, 0)) {
    sparsePolynomialFree(u);
    u = sparsePolynomialFromIntConstant(1);
  }
  sparsePolynomialFree(v);

  /* Return a gcd */
  return u;
}

static inline int __sparsePolynomialCoefficientDenominatorLcm(mpz_t lcm, sparse_polynomial_t p) {
  unsigned int i;
  mpz_t t;

  /* Handle stupid input */
  if (p == NULL) return 0;

  /* Handle the case when the polynomial has no monomials */
  if (p->monomialCount == 0u) return 1;

  /* Go over the coefficients and compute the denominator lcm */
  mpz_init(t);
  mpz_set(t, lcm);
  for (i=0u;i<p->monomialCount;i++) {
    if (!tryConstantDenominatorLcm(t, p->coeffs[i])) {
      mpz_clear(t);
      return 0;
    }
  }
  mpz_set(lcm, t);
  mpz_clear(t);
  return 1;
}

static inline int __sparsePolynomialCoefficientIntegerBounds(mp_bitcnt_t *b, sparse_polynomial_t p) {
  unsigned int i;
  mp_bitcnt_t bound, bc;

  /* Handle stupid input */
  if (p == NULL) return 0;

  /* Handle the case when the polynomial has no monomials */
  if (p->monomialCount == 0u) {
    *b = (mp_bitcnt_t) 1u;
    return 1;
  }

  /* Go over the coefficients and get their integer size */
  bound = (mp_bitcnt_t) 1u;
  for (i=0u;i<p->monomialCount;i++) {
    if (!tryConstantIntegerBound(&bc, p->coeffs[i])) return 0;
    if (bc > bound) bound = bc;
  }
  *b = bound;
  return 1;
}

static inline void __sparsePolynomialGcdHeuristicGenPolyModsDiv(mpz_t c, mpz_t g, mpz_t zeta, mpz_t scratch, mpz_t scratch2) {
  mpz_fdiv_qr(scratch, c, g, zeta);
  mpz_mul_ui(scratch2, c, 2u);
  if (mpz_cmpabs(scratch2, zeta) >= 0) {
    mpz_sub(c, c, zeta);
    mpz_add_ui(scratch, scratch, 1u);
  }
  mpz_set(g, scratch);
}

static inline sparse_polynomial_t __sparsePolynomialGcdHeuristicGenPoly(mpz_t gamma, mpz_t zeta) {
  sparse_polynomial_t p, t, q;
  int i;
  mpz_t g, c, scratch, scratch2;
  constant_t cc, ci;

  p = sparsePolynomialFromIntConstant(0);
  i = 0;
  mpz_init(g);
  mpz_set(g, gamma);
  mpz_init(c);
  mpz_init(scratch);
  mpz_init(scratch2);
  while (mpz_sgn(g) != 0) {
    __sparsePolynomialGcdHeuristicGenPolyModsDiv(c, g, zeta, scratch, scratch2);
    cc = constantFromMpz(c);
    ci = constantFromInt(i);
    t = __sparsePolynomialFromMonomial(cc, ci);
    constantFree(cc);
    constantFree(ci);
    q = sparsePolynomialAdd(p, t);
    sparsePolynomialFree(p);
    sparsePolynomialFree(t);
    p = q;
    i++;
  }
  mpz_clear(g);
  mpz_clear(c);
  mpz_clear(scratch);
  mpz_clear(scratch2);
    
  return p;
}

static inline int __sparsePolynomialGcdHeuristicDivides(sparse_polynomial_t a, sparse_polynomial_t b) {
  sparse_polynomial_t q, r;
  int res;

  sparsePolynomialDiv(&q, &r, b, a);
  res = sparsePolynomialIsConstantZero(r, 0);
  sparsePolynomialFree(q);
  sparsePolynomialFree(r);
  return res;
}

static inline int __sparsePolynomialGcdHeuristic(sparse_polynomial_t *u, sparse_polynomial_t p, sparse_polynomial_t q) {
  int n;
  mpz_t lcm;
  sparse_polynomial_t g, h, lcmp, G;
  mp_bitcnt_t b, b1, b2;
  int i;
  mpz_t zeta, alpha, beta, gamma;
  
  /* Handle stupid inputs */
  if (p == NULL) return 0;
  if (q == NULL) return 0;

  /* Both polynomials must not have too many coefficients */
  if (p->monomialCount > SPARSE_POLYNOMIAL_GCD_HEURISTIC_MAX_MONOMIAL_COUNT) return 0;
  if (q->monomialCount > SPARSE_POLYNOMIAL_GCD_HEURISTIC_MAX_MONOMIAL_COUNT) return 0;

  /* The degrees of both polynomials must not be greater than some bound */
  n = sparsePolynomialGetDegreeAsInt(p);
  if (n < 0) return 0;
  if (n > SPARSE_POLYNOMIAL_GCD_HEURISTIC_MAX_DEGREE) return 0;
  n = sparsePolynomialGetDegreeAsInt(q);
  if (n < 0) return 0;
  if (n > SPARSE_POLYNOMIAL_GCD_HEURISTIC_MAX_DEGREE) return 0;
  
  /* Both polynomials must have rational coefficients */
  if (!sparsePolynomialCoefficientsAreRational(p, 0)) return 0;
  if (!sparsePolynomialCoefficientsAreRational(q, 0)) return 0;

  /* Compute the lcm of the denominators of all coefficients of both polynomials */
  mpz_init(lcm);
  mpz_set_ui(lcm, 1u);
  if (!__sparsePolynomialCoefficientDenominatorLcm(lcm, p)) {
    /* Could not compute lcm */
    mpz_clear(lcm);
    return 0;
  }
  if (!__sparsePolynomialCoefficientDenominatorLcm(lcm, q)) {
    /* Could not compute lcm */
    mpz_clear(lcm);
    return 0;
  }

  /* Multiply both polynomials by the lcm, to make them polynomials
     with integer coefficients 
  */
  lcmp = sparsePolynomialFromMpzConstant(lcm);
  mpz_clear(lcm);
  g = sparsePolynomialMul(p, lcmp);
  h = sparsePolynomialMul(q, lcmp);  
  sparsePolynomialFree(lcmp);

  /* Now get the maximum number of bits needed to store the
     coefficients of both polynomials 
  */
  if (!__sparsePolynomialCoefficientIntegerBounds(&b1, g)) {
    /* Could not compute the bound */
    sparsePolynomialFree(g);
    sparsePolynomialFree(h);
    return 0;
  }
  if (!__sparsePolynomialCoefficientIntegerBounds(&b2, h)) {
    /* Could not compute the bound */
    sparsePolynomialFree(g);
    sparsePolynomialFree(h);
    return 0;
  }
  b = b1;
  if (b2 > b1) b = b2;

  /* The bound on the integer coefficients must not be too large */
  if (b > SPARSE_POLYNOMIAL_GCD_HEURISTIC_MAX_INTEGER_BOUND) {
    sparsePolynomialFree(g);
    sparsePolynomialFree(h);
    return 0;
  }

  /* Now take zeta = 2^(b + 5) + 1. Clearly zeta > 2 * abs(c_i) for
     all coefficients c_i of both polynomials. 
     
     See Theorem 1 of 

     Bruce W. Char, Keith O. Geddes, Gaston H. Gonnet, GCDHEU:
     Heuristic polynomial GCD algorithm based on integer GCD
     computation, Journal of Symbolic Computation, Volume 7, Issue 1,
     1989, Pages 31-48.

     to understand why this is important.

  */
  mpz_init(zeta);
  mpz_ui_pow_ui(zeta, 2u, (unsigned int) (b + 5u));
  mpz_add_ui(zeta, zeta, 1u);

  /* Initialize alpha, beta and gamma */
  mpz_init(alpha);
  mpz_init(beta);
  mpz_init(gamma);
  
  /* Loop for the heuristic */
  for (i=0;i<SPARSE_POLYNOMIAL_GCD_HEURISTIC_TRIALS;i++) {
    /* Evaluate both polynomials at zeta, yielding alpha and beta 
       
       alpha = g(zeta)
       beta  = h(zeta)

    */
    if (!sparsePolynomialEvalMpz(alpha, g, zeta)) {
      /* Evaluation did not work */
      mpz_clear(zeta);
      mpz_clear(alpha);
      mpz_clear(beta);
      mpz_clear(gamma);
      sparsePolynomialFree(g);
      sparsePolynomialFree(h);
      return 0;
    }
    if (!sparsePolynomialEvalMpz(beta, h, zeta)) {
      /* Evaluation did not work */
      mpz_clear(zeta);
      mpz_clear(alpha);
      mpz_clear(beta);
      mpz_clear(gamma);
      sparsePolynomialFree(g);
      sparsePolynomialFree(h);
      return 0;
    }
    
    /* Compute integer gcd: gamma = gcd(alpha, beta) */
    mpz_gcd(gamma, alpha, beta);

    /* If gamma is zero, continue directly with the next zeta */
    if (mpz_sgn(gamma) != 0) {
      /* Lift the integer gcd onto the polynomials */
      G = __sparsePolynomialGcdHeuristicGenPoly(gamma, zeta);

      /* Check if G divides g and h */
      if (__sparsePolynomialGcdHeuristicDivides(G, g) &&
	  __sparsePolynomialGcdHeuristicDivides(G, h)) {
	/* G divides both g and h and hence is the gcd of g and h 

	   Set the result.
	*/
	*u = G;

	/* Clear the variables */
	mpz_clear(zeta);
	mpz_clear(alpha);
	mpz_clear(beta);
	mpz_clear(gamma);
	sparsePolynomialFree(g);
	sparsePolynomialFree(h);
	
	/* Return success */
	return 1;
      }
      /* The heuristic did not work, clear G */
      sparsePolynomialFree(G);
    }

    /* Prepare the next step */
    mpz_mul_ui(zeta, zeta, 4u);
    mpz_add_ui(zeta, zeta, 1u);
  }

  /* Clear variables */
  mpz_clear(zeta);
  mpz_clear(alpha);
  mpz_clear(beta);
  mpz_clear(gamma);
  sparsePolynomialFree(g);
  sparsePolynomialFree(h);

  /* The heuristic did not work */
  return 0;
}

static inline sparse_polynomial_t sparsePolynomialGcd(sparse_polynomial_t p, sparse_polynomial_t q) {
  sparse_polynomial_t u, v, t, z;
  constant_t a, b, c, d;
  
  /* Handle stupid inputs */
  if (p == NULL) return NULL;
  if (q == NULL) return NULL;
  
  /* Handle the case when both polynomials are constants.  

     In this case return the constant polynomial built from the gcd of
     these two constants. See the definition for gcd of two constants
     above.
  */
  if (sparsePolynomialIsConstant(p, 0) &&
      sparsePolynomialIsConstant(q, 0)) {
    a = sparsePolynomialGetIthCoefficientAsConstantIntIndex(p, 0);
    b = sparsePolynomialGetIthCoefficientAsConstantIntIndex(q, 0);
    c = constantGcd(a, b);
    u = sparsePolynomialFromConstant(c);
    constantFree(c);
    constantFree(b);
    constantFree(a);
    return u;
  }

  /* Handle the case when the input polynomials are equal. */
  if (sparsePolynomialEqual(p, q, 0)) {
    return sparsePolynomialFromCopy(p);
  }

  /* Try a heuristic algorithm. If it does not work, use the general
     Euclidian algorithm.
  */
  if (!__sparsePolynomialGcdHeuristic(&u, p, q)) {
    /* General case: Euclidian algorithm */
    u = __sparsePolynomialGcdBaseCase(p, q);
  }

  /* Make the polynomial u = gcd(a,b) unitary */
  if (!sparsePolynomialIsConstantZero(u, 1)) {
    __sparsePolynomialGetLeadingCoefficient(&c, &d, &z, u);
    sparsePolynomialFree(z);
    constantFree(d);
    v = sparsePolynomialFromConstant(c);
    constantFree(c);
    sparsePolynomialDiv(&z, &t, u, v);
    sparsePolynomialFree(t);
    sparsePolynomialFree(v);
    sparsePolynomialFree(u);
    u = z;
  }

  /* Normalize the polynomial u = gcd(p,q) in such a way that its
     leading coefficient is the gcd of the leading coefficients of the
     input polynomials p and q.
  */
  __sparsePolynomialGetLeadingCoefficient(&a, &d, &z, p);
  constantFree(d);
  sparsePolynomialFree(z);
  __sparsePolynomialGetLeadingCoefficient(&b, &d, &z, q);
  constantFree(d);
  sparsePolynomialFree(z);
  c = constantGcd(a, b);
  v = sparsePolynomialFromConstant(c);
  constantFree(a);
  constantFree(b);
  constantFree(c);
  z = sparsePolynomialMul(u, v);
  sparsePolynomialFree(v);
  sparsePolynomialFree(u);
  u = z;
  
  /* Return the normalized u = gcd(p,q) */
  return u;
}

int sparsePolynomialFromExpression(sparse_polynomial_t *r, node *p) {
  sparse_polynomial_t a, b, quot, rest;
  int res;

  /* Handle stupid inputs */
  if (p == NULL) return 0;  

  /* Try to decompose expressions built on c, x, +, -, *, /, -, ^ */
  switch (p->nodeType) {
  case MEMREF:
    if (p->polynomialRepresentation != NULL) {
      *r = __polynomialGetSparsePolynomial(p->polynomialRepresentation);
      return 1;
    }
    return sparsePolynomialFromExpression(r, getMemRefChild(p));
    break;
  case VARIABLE:
    *r = sparsePolynomialFromIdentity();
    return 1;
    break;
  case CONSTANT:
    *r = sparsePolynomialFromMpfrConstant(*(p->value));
    return 1;
    break;
  case NEG:
    if (!sparsePolynomialFromExpression(&a, p->child1)) 
      break;
    *r = sparsePolynomialNeg(a);
    sparsePolynomialFree(a);
    return 1;
    break;
  case ADD:
  case SUB:
  case MUL:
  case DIV:
  case POW:
    if (!sparsePolynomialFromExpression(&a, p->child1)) 
      break;
    if (!sparsePolynomialFromExpression(&b, p->child2)) {
      sparsePolynomialFree(a);
      break;
    }
    res = 0;
    switch (p->nodeType) {
    case ADD:
      *r = sparsePolynomialAdd(a, b);
      res = 1;
      break;
    case SUB:
      *r = sparsePolynomialSub(a, b);
      res = 1;
      break;
    case MUL:
      *r = sparsePolynomialMul(a, b);
      res = 1;
      break;
    case DIV:
      if (sparsePolynomialIsConstantZero(b, 1)) {
	res = 0;
      } else {
	sparsePolynomialDiv(&quot, &rest, a, b);
	if (sparsePolynomialIsConstantZero(rest, 0)) {
	  *r = quot;
	  res = 1;
	} else {
	  sparsePolynomialFree(quot);
	  res = 0;
	}
	sparsePolynomialFree(rest);
      }
      break;
    case POW:
      res = sparsePolynomialPow(r, a, b);
      break;
    }
    sparsePolynomialFree(a);
    sparsePolynomialFree(b);
    if (res) return 1;
    break;
  }

  /* Here, the expression contains some other basic function. It can
     be a polynomial only if it is a constant expression.
  */
  if (isConstant(p)) {
    if (!__sparsePolynomialFromConstantExpression(r, p)) return 0;
    return 1;
  }

  return 0;
}

static inline void sparsePolynomialGetDegree(mpz_t deg, sparse_polynomial_t p) {
  if (p == NULL) {
    mpz_set_si(deg, -1);
    return;
  }
  if (!tryConstantToMpz(deg, p->deg)) {
    mpz_set_si(deg, -1);
  }
}

static inline int sparsePolynomialGetDegreeAsInt(sparse_polynomial_t p) {
  int deg;

  /* Make compiler happy */
  deg = -1;
  /* End of compiler happiness */

  if (p == NULL) return -1;
  if (!tryConstantToInt(&deg, p->deg)) return -1;
  return deg;
}

static inline constant_t sparsePolynomialGetIthCoefficientAsConstant(sparse_polynomial_t p, mpz_t i) {
  constant_t ic, coeffsum, t;
  unsigned int j, k;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;

  /* The index must be non-negative */
  if (mpz_sgn(i) < 0) {
    return constantFromInt(0);
  }

  /* Handle the strange case when p has no monomials */
  if (p->monomialCount == 0u) {
    return constantFromInt(0);
  }

  /* Construct a constant from i */
  ic = constantFromMpz(i);

  /* The index must be no greater than the degree of p */
  if (constantIsGreater(ic, p->deg, 0)) {
    constantFree(ic);
    return constantFromInt(0);
  }

  /* Get the index j where the coefficients for degree i might
     start 
  */
  j = __sparsePolynomialFindDegree(ic, p->monomialDegrees, p->monomialCount, 0u); 

  /* If j is greater than the last index of the monomials of p, there
     exists no monomial of degree i.
  */
  if (j >= p->monomialCount) {
    constantFree(ic);
    return constantFromInt(0);
  }

  /* Here, we know that 0 <= j <= n - 1.

     We add up all coefficients associated with a degree equal to i.

  */
  coeffsum = constantFromInt(0);
  for (k=j;k<p->monomialCount;k++) {
    if (!constantIsEqual(ic, p->monomialDegrees[k], 0)) 
      break;
    t = constantAdd(coeffsum, p->coeffs[k]);
    constantFree(coeffsum);
    coeffsum = t;
  }
  constantFree(ic);
  return coeffsum;
}

node *sparsePolynomialGetIthCoefficient(sparse_polynomial_t p, mpz_t i) {
  constant_t c;
  node *res;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;

  /* Get coefficient as a constant */
  c = sparsePolynomialGetIthCoefficientAsConstant(p, i);

  /* Convert to an expression */
  res = addMemRef(constantToExpression(c));
  constantFree(c);

  /* Return the result */
  return res;
}

node *sparsePolynomialGetIthCoefficientIntIndex(sparse_polynomial_t p, int i) {
  constant_t ic, coeffsum, t;
  unsigned int j, k;
  node *res;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;

  /* The index must be non-negative */
  if (i < 0) {
    return addMemRef(makeConstantInt(0));
  }

  /* Handle the strange case when p has no monomials */
  if (p->monomialCount == 0u) {
    return addMemRef(makeConstantInt(0));
  }

  /* Construct a constant from i */
  ic = constantFromInt(i);

  /* The index must be no greater than the degree of p */
  if (constantIsGreater(ic, p->deg, 0)) {
    constantFree(ic);
    return addMemRef(makeConstantInt(0));
  }

  /* Get the index j where the coefficients for degree i might
     start 
  */
  j = __sparsePolynomialFindDegree(ic, p->monomialDegrees, p->monomialCount, 0u); 

  /* If j is greater than the last index of the monomials of p, there
     exists no monomial of degree i.
  */
  if (j >= p->monomialCount) {
    constantFree(ic);
    return addMemRef(makeConstantInt(0));
  }

  /* Here, we know that 0 <= j <= n - 1.

     We add up all coefficients associated with a degree equal to i.

  */
  coeffsum = constantFromInt(0);
  for (k=j;k<p->monomialCount;k++) {
    if (!constantIsEqual(ic, p->monomialDegrees[k], 0)) 
      break;
    t = constantAdd(coeffsum, p->coeffs[k]);
    constantFree(coeffsum);
    coeffsum = t;
  }
  constantFree(ic);
  res = constantToExpression(coeffsum);
  constantFree(coeffsum);
  
  return res;
}

static inline constant_t sparsePolynomialGetIthCoefficientAsConstantIntIndex(sparse_polynomial_t p, int i) {
  constant_t ic, coeffsum, t;
  unsigned int j, k;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;

  /* The index must be non-negative */
  if (i < 0) {
    return constantFromInt(0);
  }

  /* Handle the strange case when p has no monomials */
  if (p->monomialCount == 0u) {
    return constantFromInt(0);
  }

  /* Construct a constant from i */
  ic = constantFromInt(i);

  /* The index must be no greater than the degree of p */
  if (constantIsGreater(ic, p->deg, 0)) {
    constantFree(ic);
    return constantFromInt(0);
  }

  /* Get the index j where the coefficients for degree i might
     start 
  */
  j = __sparsePolynomialFindDegree(ic, p->monomialDegrees, p->monomialCount, 0u); 

  /* If j is greater than the last index of the monomials of p, there
     exists no monomial of degree i.
  */
  if (j >= p->monomialCount) {
    constantFree(ic);
    return constantFromInt(0);
  }

  /* Here, we know that 0 <= j <= n - 1.

     We add up all coefficients associated with a degree equal to i.

  */
  coeffsum = constantFromInt(0);
  for (k=j;k<p->monomialCount;k++) {
    if (!constantIsEqual(ic, p->monomialDegrees[k], 0)) 
      break;
    t = constantAdd(coeffsum, p->coeffs[k]);
    constantFree(coeffsum);
    coeffsum = t;
  }
  constantFree(ic);
  return coeffsum;
}

static inline int sparsePolynomialGetCoefficients(node ***coeffs, unsigned int *deg, sparse_polynomial_t p) {
  unsigned int degree, size, i, j, d, k;
  node **coefficients;
  constant_t c, t;

  /* Handle stupid inputs */
  if (p == NULL) {
    return 0;
  }

  /* Handle the strange case when p has no monomials */
  if (p->monomialCount == 0u) {
    *deg = 0u;
    *coeffs = (node **) safeCalloc(1, sizeof(node *));
    (*coeffs)[0] = addMemRef(makeConstantInt(0));
    return 1;
  }

  /* Check if the degree of the polynomial holds on an unsigned
     int. 
  */
  if (!tryConstantToUnsignedInt(&degree, p->deg)) {
    return 0;
  }
  
  /* Check if the size of the output (degree + 1) holds on an unsigned
     int. 
  */
  size = degree + 1u;
  if (size == 0u) return 0;

  /* Allocate the output array and initialize all its elements to
     NULL. 
  */
  coefficients = (node **) safeCalloc(size, sizeof(node *));
  for (i=0u;i<=degree;i++) {
    coefficients[i] = NULL;
  }

  /* Go over all monomials, add up equal degrees and set the output
     coefficients. 

     If ever there is a degree we can't place, desallocate everything
     and signal failure.
  */
  for (i=0u;i<p->monomialCount;i++) {
    if (!tryConstantToUnsignedInt(&d, p->monomialDegrees[i])) {
      for (j=0u;j<=degree;j++) {
	if (coefficients[j] != NULL) freeThing(coefficients[j]);
      }
      safeFree(coefficients);
      return 0;
    }
    if (d > degree) {
      for (j=0u;j<=degree;j++) {
	if (coefficients[j] != NULL) freeThing(coefficients[j]);
      }
      safeFree(coefficients);
      return 0;
    }
    for (j=i;j<p->monomialCount-1u;j++) {
      if (!constantIsEqual(p->monomialDegrees[i], p->monomialDegrees[j+1u], 0)) 
	break;
    }
    if (i == j) {
      c = constantFromCopy(p->coeffs[i]);
    } else {
      c = constantFromInt(0);
      for (k=i;k<=j;k++) {
	t = constantAdd(c, p->coeffs[k]);
	constantFree(c);
	c = t;
      }
      i = j;
    }
    if (coefficients[d] == NULL) {
      coefficients[d] = constantToExpression(c);
    } else {
      coefficients[d] = addMemRef(makeAdd(coefficients[d], constantToExpression(c)));
    }
    constantFree(c);
  }

  /* Fill in the holes in the output coefficient array and set the
     output 
  */
  for (i=0u;i<=degree;i++) {
    if (coefficients[i] == NULL) {
      coefficients[i] = addMemRef(makeConstantInt(0));
    }
  }
  *coeffs = coefficients;
  *deg = degree;
  return 1;
}

static inline unsigned int sparsePolynomialGetNumberOfNonZeroCoefficients(sparse_polynomial_t p) {

  /* Handle stupid inputs */
  if (p == NULL) {
    return 0u;
  }

  /* General case */
  return p->monomialCount;
}

static inline void __sparsePolynomialEvalMpfr(mpfr_t y, sparse_polynomial_t p, mpfr_t x, mpfr_t scratch) {
  unsigned int i, a, b, d;
  mp_prec_t prec;
  constant_t dc;

  /* Handle stupid inputs */
  if (p == NULL) {
    mpfr_set_nan(y);
    return;
  }

  /* Handle the strange case when p has no monomials */
  if (p->monomialCount == 0u) {
    mpfr_set_si(y, 0, GMP_RNDN);
    return;
  }

  /* Save precision of y and increase it for a moment */
  prec = mpfr_get_prec(y);
  mpfr_set_prec(y, prec + 25);
  mpfr_set_prec(y, prec + 10);

  /* Set precision of scratch variable */
  mpfr_set_prec(scratch, prec + 25);
  mpfr_set_prec(scratch, prec + 10);
  
  /* Perform Horner evaluation of p in x */
  constantEvalMpfr(y, p->coeffs[p->monomialCount-1u]);
  for (i=p->monomialCount-1u;i>=1u;i--) {
    /* y <- x^(p->monomialDegrees[i] - p->monomialDegrees[i-1u]) * y */
    if (tryConstantToUnsignedInt(&a, p->monomialDegrees[i]) &&
	tryConstantToUnsignedInt(&b, p->monomialDegrees[i-1u])) {
      if (a < b) {
	sollyaFprintf(stderr,"Error: __sparsePolynomialEvalMpfr: monomial degrees not appropriately ordered\n");
	exit(1);
      }
      d = a - b;
      if (d != 0u) {
	if (d == 1u) {
	  mpfr_mul(y, x, y, GMP_RNDN);
	} else {
	  mpfr_pow_ui(scratch, x, d, GMP_RNDN);
	  mpfr_mul(y, scratch, y, GMP_RNDN);
	}
      }
    } else {
      dc = constantSub(p->monomialDegrees[i], p->monomialDegrees[i-1u]);
      constantEvalMpfr(scratch, dc);
      constantFree(dc);
      mpfr_pow(scratch, x, scratch, GMP_RNDN);
      mpfr_mul(y, scratch, y, GMP_RNDN);
    }
    /* y <- p->coeffs[i-1u] + y */
    constantEvalMpfr(scratch, p->coeffs[i-1u]);
    mpfr_add(y, y, scratch, GMP_RNDN);
  }
  /* y <- x^(p->monomialDegrees[0u]) * y */
  if (tryConstantToUnsignedInt(&a, p->monomialDegrees[0u])) {
    if (a != 0u) {
      if (a == 1u) {
	mpfr_mul(y, x, y, GMP_RNDN);
      } else {
	mpfr_pow_ui(scratch, x, a, GMP_RNDN);
	mpfr_mul(y, scratch, y, GMP_RNDN);
      }
    }
  } else {
    constantEvalMpfr(scratch, p->monomialDegrees[0u]);
    mpfr_pow(scratch, x, scratch, GMP_RNDN);
    mpfr_mul(y, scratch, y, GMP_RNDN);
  }

  /* Round result to original precision of y */
  mpfr_prec_round(y, prec, GMP_RNDN);
}

static inline void sparsePolynomialEvalMpfr(mpfr_t y, sparse_polynomial_t p, mpfr_t x) {
  mpfr_t scratch, Y;

  /* Handle stupid inputs */
  if (p == NULL) {
    mpfr_set_nan(y);
    return;
  }

  /* Handle the strange case when p has no monomials */
  if (p->monomialCount == 0u) {
    mpfr_set_si(y, 0, GMP_RNDN);
    return;
  }

  /* Cover the case when x and y are the same MPFR variable */
  if (x == y) {
    mpfr_init2(Y, mpfr_get_prec(y));
    sparsePolynomialEvalMpfr(Y, p, x);
    mpfr_set(y, Y, GMP_RNDN); /* exact */
    mpfr_clear(Y);
    return;
  }

  /* Try to use a global variable as scratch space */
  if  (__sparsePolynomialEvalMpfr_var_used) {
    /* The global scratch space is taken; use new scratch space on the
       stack 
    */
    mpfr_init2(scratch, mpfr_get_prec(y) + 25);

    /* Use inner evaluation function */
    __sparsePolynomialEvalMpfr(y, p, x, scratch);
    
    /* Clear the scratch variable */
    mpfr_clear(scratch);

    return;
  }

  /* Here the global scratch space variable still is free; so start
     using it. 
  */
  __sparsePolynomialEvalMpfi_var_used = 1;

  /* Initialize the global scratch variable */
  if (__sparsePolynomialEvalMpfr_scratch_initialized) {
    mpfr_set_prec(__sparsePolynomialEvalMpfr_scratch, mpfr_get_prec(y) + 25);
  } else {
    mpfr_init2(__sparsePolynomialEvalMpfr_scratch, mpfr_get_prec(y) + 25);
    __sparsePolynomialEvalMpfr_scratch_initialized = 1;
  }
  
  /* Use inner evaluation function */
  __sparsePolynomialEvalMpfr(y, p, x, __sparsePolynomialEvalMpfr_scratch);

  /* Free the access to the global scratch space variable again */
  __sparsePolynomialEvalMpfi_var_used = 0;
}

static inline void __sparsePolynomialEvalMpfi(sollya_mpfi_t y, sparse_polynomial_t p, sollya_mpfi_t x, sollya_mpfi_t scratch) {
  unsigned int i, a, b, d;
  mp_prec_t prec;
  constant_t dc;

  /* Handle stupid inputs */
  if (p == NULL) {
    sollya_mpfi_set_nan(y);
    return;
  }

  /* Handle the strange case when p has no monomials */
  if (p->monomialCount == 0u) {
    sollya_mpfi_set_si(y, 0);
    return;
  }

  /* Save precision of y and increase it for a moment */
  prec = sollya_mpfi_get_prec(y);
  sollya_mpfi_set_prec(y, prec + 25);
  sollya_mpfi_set_prec(y, prec + 10);

  /* Set precision of scratch variable */
  sollya_mpfi_set_prec(scratch, prec + 25);
  sollya_mpfi_set_prec(scratch, prec + 10);
  
  /* Perform Horner evaluation of p in x */
  constantEvalMpfi(y, p->coeffs[p->monomialCount-1u]);
  for (i=p->monomialCount-1u;i>=1u;i--) {
    /* y <- x^(p->monomialDegrees[i] - p->monomialDegrees[i-1u]) * y */
    if (tryConstantToUnsignedInt(&a, p->monomialDegrees[i]) &&
	tryConstantToUnsignedInt(&b, p->monomialDegrees[i-1u])) {
      if (a < b) {
	sollyaFprintf(stderr,"Error: __sparsePolynomialEvalMpfi: monomial degrees not appropriately ordered\n");
	exit(1);
      }
      d = a - b;
      if (d != 0u) {
	if (d == 1u) {
	  sollya_mpfi_mul(y, x, y);
	} else {
	  sollya_mpfi_pow_ulong(scratch, x, d);
	  sollya_mpfi_mul(y, scratch, y);
	}
      }
    } else {
      dc = constantSub(p->monomialDegrees[i], p->monomialDegrees[i-1u]);
      constantEvalMpfi(scratch, dc);
      constantFree(dc);
      sollya_mpfi_pow(scratch, x, scratch);
      sollya_mpfi_mul(y, scratch, y);
    }
    /* y <- p->coeffs[i-1u] + y */
    constantEvalMpfi(scratch, p->coeffs[i-1u]);
    sollya_mpfi_add(y, y, scratch);
  }
  /* y <- x^(p->monomialDegrees[0u]) * y */
  if (tryConstantToUnsignedInt(&a, p->monomialDegrees[0u])) {
    if (a != 0u) {
      if (a == 1u) {
	sollya_mpfi_mul(y, x, y);
      } else {
	sollya_mpfi_pow_ulong(scratch, x, a);
	sollya_mpfi_mul(y, scratch, y);
      }
    }
  } else {
    constantEvalMpfi(scratch, p->monomialDegrees[0u]);
    sollya_mpfi_pow(scratch, x, scratch);
    sollya_mpfi_mul(y, scratch, y);
  }

  /* Round result to original precision of y */
  sollya_mpfi_prec_round(y, prec);
}

static inline void sparsePolynomialEvalMpfi(sollya_mpfi_t y, sparse_polynomial_t p, sollya_mpfi_t x) {
  sollya_mpfi_t scratch, Y;

  /* Handle stupid inputs */
  if (p == NULL) {
    sollya_mpfi_set_nan(y);
    return;
  }

  /* Handle the strange case when p has no monomials */
  if (p->monomialCount == 0u) {
    sollya_mpfi_set_si(y, 0);
    return;
  }

  /* Cover the case when x and y are the same MPFI variable */
  if (x == y) {
    sollya_mpfi_init2(Y, sollya_mpfi_get_prec(y));
    sparsePolynomialEvalMpfi(Y, p, x);
    sollya_mpfi_set(y, Y);
    sollya_mpfi_clear(Y);
    return;
  }

  /* Try to use a global variable as scratch space */
  if  (__sparsePolynomialEvalMpfi_var_used) {
    /* The global scratch space is taken; use new scratch space on the
       stack 
    */
    sollya_mpfi_init2(scratch, sollya_mpfi_get_prec(y) + 25);

    /* Use inner evaluation function */
    __sparsePolynomialEvalMpfi(y, p, x, scratch);
    
    /* Clear the scratch variable */
    sollya_mpfi_clear(scratch);

    return;
  }

  /* Here the global scratch space variable still is free; so start
     using it. 
  */
  __sparsePolynomialEvalMpfi_var_used = 1;

  /* Initialize the global scratch variable */
  if (__sparsePolynomialEvalMpfi_scratch_initialized) {
    sollya_mpfi_set_prec(__sparsePolynomialEvalMpfi_scratch, sollya_mpfi_get_prec(y) + 25);
  } else {
    sollya_mpfi_init2(__sparsePolynomialEvalMpfi_scratch, sollya_mpfi_get_prec(y) + 25);
    __sparsePolynomialEvalMpfi_scratch_initialized = 1;
  }
  
  /* Use inner evaluation function */
  __sparsePolynomialEvalMpfi(y, p, x, __sparsePolynomialEvalMpfi_scratch);

  /* Free the access to the global scratch space variable again */
  __sparsePolynomialEvalMpfi_var_used = 0;
}

static inline int __sparsePolynomialEvalMpz(mpz_t y, sparse_polynomial_t p, mpz_t x, mpz_t scratch) {
  unsigned int i, a, b, d;
  constant_t dc;

  /* Handle stupid inputs */
  if (p == NULL) {
    return 0;
  }

  /* Handle the strange case when p has no monomials */
  if (p->monomialCount == 0u) {
    mpz_set_si(y, 0);
    return 1;
  }
  
  /* Perform Horner evaluation of p in x */
  if (!tryConstantToMpz(y, p->coeffs[p->monomialCount-1u])) return 0;
  for (i=p->monomialCount-1u;i>=1u;i--) {
    /* y <- x^(p->monomialDegrees[i] - p->monomialDegrees[i-1u]) * y */
    if (tryConstantToUnsignedInt(&a, p->monomialDegrees[i]) &&
	tryConstantToUnsignedInt(&b, p->monomialDegrees[i-1u])) {
      if (a < b) {
	sollyaFprintf(stderr,"Error: __sparsePolynomialEvalMpz: monomial degrees not appropriately ordered\n");
	exit(1);
      }
      d = a - b;
      if (d != 0u) {
	if (d == 1u) {
	  mpz_mul(y, x, y);
	} else {
	  mpz_pow_ui(scratch, x, d);
	  mpz_mul(y, scratch, y);
	}
      }
    } else {
      dc = constantSub(p->monomialDegrees[i], p->monomialDegrees[i-1u]);
      if (!tryConstantToMpz(scratch, dc)) {
	constantFree(dc);
	return 0;
      }
      constantFree(dc);
      if (!sollya_mpz_pow(scratch, x, scratch)) return 0;
      mpz_mul(y, scratch, y);
    }
    /* y <- p->coeffs[i-1u] + y */
    if (!tryConstantToMpz(scratch, p->coeffs[i-1u])) return 0;
    mpz_add(y, y, scratch);
  }
  /* y <- x^(p->monomialDegrees[0u]) * y */
  if (tryConstantToUnsignedInt(&a, p->monomialDegrees[0u])) {
    if (a != 0u) {
      if (a == 1u) {
	mpz_mul(y, x, y);
      } else {
	mpz_pow_ui(scratch, x, a);
	mpz_mul(y, scratch, y);
      }
    }
  } else {
    if (!tryConstantToMpz(scratch, p->monomialDegrees[0u])) return 0;
    if (!sollya_mpz_pow(scratch, x, scratch)) return 0;
    mpz_mul(y, scratch, y);
  }

  /* Indicate success */
  return 1;
}

static inline int sparsePolynomialEvalMpz(mpz_t y, sparse_polynomial_t p, mpz_t x) {
  mpz_t scratch, Y;
  int res;

  /* Handle stupid inputs */
  if (p == NULL) {
    return 0;
  }

  /* Handle the strange case when p has no monomials */
  if (p->monomialCount == 0u) {
    mpz_set_si(y, 0);
    return 1;
  }

  /* Cover the case when x and y are the same MPFI variable */
  if (x == y) {
    mpz_init(Y);
    res = sparsePolynomialEvalMpz(Y, p, x);
    mpz_set(y, Y);
    mpz_clear(Y);
    return res;
  }

  /* Initialize a scratch variable */
  mpz_init(scratch);

  /* Use inner evaluation function */
  res = __sparsePolynomialEvalMpz(y, p, x, scratch);
    
  /* Clear the scratch variable */
  mpz_clear(scratch);

  /* Return success flag */
  return res;
}

static inline node *__sparsePolynomialGetExpressionCanonical(sparse_polynomial_t p) {
  node *res, *temp;
  unsigned int i;

  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* Handle the case when p has no monomials */
  if (p->monomialCount == 0u) {
    return addMemRef(makeConstantInt(0));
  }

  /* Handle the case when p only has one monomial */
  if (p->monomialCount == 1u) {
    if (constantIsZero(p->monomialDegrees[0],0)) {
      return addMemRef(constantToExpression(p->coeffs[0]));
    } 
    if (constantIsOne(p->monomialDegrees[0],0)) {
      if (constantIsOne(p->coeffs[0],0)) {
	return addMemRef(makeVariable());
      }
      return addMemRef(makeMul(constantToExpression(p->coeffs[0]),
			       makeVariable()));
    }
    if (constantIsOne(p->coeffs[0],0)) {
      return addMemRef(makePow(makeVariable(),
			       constantToExpression(p->monomialDegrees[0])));
    }
    return addMemRef(makeMul(constantToExpression(p->coeffs[0]),
			     makePow(makeVariable(),
				     constantToExpression(p->monomialDegrees[0]))));
  }

  /* Here, p has at least two monomials */
  if (constantIsZero(p->monomialDegrees[0],0)) {
    res = constantToExpression(p->coeffs[0]);
  } else {
    if (constantIsOne(p->monomialDegrees[0],0)) {
      if (constantIsOne(p->coeffs[0], 0)) {
	res = makeVariable();
      } else {
	res = makeMul(constantToExpression(p->coeffs[0]),
		      makeVariable());
      }
    } else {
      if (constantIsOne(p->coeffs[0], 0)) {
	res = makePow(makeVariable(),
		      constantToExpression(p->monomialDegrees[0]));
      } else {
	res = makeMul(constantToExpression(p->coeffs[0]),
		      makePow(makeVariable(),
			      constantToExpression(p->monomialDegrees[0])));
      }
    }
  }
  for (i=1u;i<p->monomialCount;i++) {
    if (constantIsZero(p->monomialDegrees[i],0)) {
      temp = constantToExpression(p->coeffs[i]);
    } else {
      if (constantIsOne(p->monomialDegrees[i],0)) {
	if (constantIsOne(p->coeffs[i], 0)) {
	  temp = makeVariable();
	} else {
	  temp = makeMul(constantToExpression(p->coeffs[i]),
			 makeVariable());
	}
      } else {
	if (constantIsOne(p->coeffs[i], 0)) {
	  temp = makePow(makeVariable(),
			 constantToExpression(p->monomialDegrees[i]));	  
	} else {
	  temp = makeMul(constantToExpression(p->coeffs[i]),
			 makePow(makeVariable(),
				 constantToExpression(p->monomialDegrees[i])));
	}
      }
    }
    res = makeAdd(res, temp);
  }

  return addMemRef(res);
}

static inline node *__sparsePolynomialGetExpressionHorner(sparse_polynomial_t p) {
  unsigned int i, a, b, d;
  constant_t dc;
  node *res;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;

  /* Handle the strange case when p has no monomials */
  if (p->monomialCount == 0u) return addMemRef(makeConstantInt(0)); 
  
  /* Handle the case when the polynomial only has one monomial */
  if (p->monomialCount == 1u) {
    if (constantIsZero(p->monomialDegrees[0u], 0)) {
      return addMemRef(constantToExpression(p->coeffs[0u]));
    } else {
      if (constantIsOne(p->monomialDegrees[0u], 0)) {
	if (constantIsOne(p->coeffs[0u], 0)) {
	  return addMemRef(makeVariable());
	} else {
	  return addMemRef(makeMul(makeVariable(), 
				   constantToExpression(p->coeffs[0u])));
	}
      } else {
	if (constantIsOne(p->coeffs[0u], 0)) {
	  return addMemRef(makePow(makeVariable(),
				   constantToExpression(p->monomialDegrees[0u])));
	} else {
	  return addMemRef(makeMul(makePow(makeVariable(),
					   constantToExpression(p->monomialDegrees[0u])),
				   constantToExpression(p->coeffs[0u])));
	}
      }
    }
  }

  /* Perform Horner "evaluation" of p 

     Here, the polynomial has at least 2 monomials. 
     
  */
  if (constantIsOne(p->coeffs[p->monomialCount-1u], 0)) {
    /* Use NULL as a marker for a coefficient equal to one */
    res = NULL;
  } else {
    res = constantToExpression(p->coeffs[p->monomialCount-1u]);
  }
  for (i=p->monomialCount-1u;i>=1u;i--) {
    /* y <- x^(p->monomialDegrees[i] - p->monomialDegrees[i-1u]) * y */
    if (tryConstantToUnsignedInt(&a, p->monomialDegrees[i]) &&
	tryConstantToUnsignedInt(&b, p->monomialDegrees[i-1u])) {
      if (a < b) {
	sollyaFprintf(stderr,"Error: __sparsePolynomialGetExpressionHorner: monomial degrees not appropriately ordered\n");
	exit(1);
      }
      d = a - b;
      if (d != 0u) {
	if (d == 1u) {
	  if (res == NULL) {
	    res = makeVariable();
	  } else {
	    res = makeMul(makeVariable(), res);
	  }
	} else {
	  if (res == NULL) {
	    res = makePow(makeVariable(),makeConstantInt(d));
	  } else {
	    res = makeMul(makePow(makeVariable(),makeConstantInt(d)),res);
	  }
	}
      }
    } else {
      dc = constantSub(p->monomialDegrees[i], p->monomialDegrees[i-1u]);
      if (!constantIsZero(dc,0)) {
	if (constantIsOne(dc,0)) {
	  if (res == NULL) {
	    res = makeVariable();
	  } else {
	    res = makeMul(makeVariable(), res);
	  }
	} else {
	  if (res == NULL) {
	    res = makePow(makeVariable(),constantToExpression(dc));
	  } else {
	    res = makeMul(makePow(makeVariable(),constantToExpression(dc)),res);
	  }
	}
      }
      constantFree(dc);
    }
    /* y <- p->coeffs[i-1u] + y */
    if (res == NULL) {
      res = makeAdd(constantToExpression(p->coeffs[i-1u]), makeConstantInt(1));
    } else {
      res = makeAdd(constantToExpression(p->coeffs[i-1u]), res);
    }
  }
  /* y <- x^(p->monomialDegrees[0u]) * y */
  if (tryConstantToUnsignedInt(&a, p->monomialDegrees[0u])) {
    if (a != 0u) {
      if (a == 1u) {
	if (res == NULL) {
	  res = makeVariable();
	} else {
	  res = makeMul(makeVariable(), res);
	}
      } else {
	if (res == NULL) {
	  res = makePow(makeVariable(),
			makeConstantInt(a));
	} else {
	  res = makeMul(makePow(makeVariable(),
				makeConstantInt(a)),
			res);
	}
      }
    }
  } else {
    if (!constantIsZero(p->monomialDegrees[0u],0)) {
      if (constantIsOne(p->monomialDegrees[0u],0)) {
	if (res == NULL) {
	  res = makeVariable();
	} else {
	  res = makeMul(makeVariable(), res);
	}
      } else {
	if (res == NULL) {
	  res = makePow(makeVariable(),
			constantToExpression(p->monomialDegrees[0u]));
	} else {
	  res = makeMul(makePow(makeVariable(),
				constantToExpression(p->monomialDegrees[0u])),
			res);
	}
      }
    }
  }
  
  if (res == NULL) {
    res = makeConstantInt(1);
  }

  /* Return the result */
  return addMemRef(res);
}

static inline node *sparsePolynomialGetExpression(sparse_polynomial_t p, int canonical) {
  if (canonical) return __sparsePolynomialGetExpressionCanonical(p);
  return __sparsePolynomialGetExpressionHorner(p);
}

void sparsePolynomialFPrintf(FILE *fd, sparse_polynomial_t p, int canonical) {
  node *t;

  /* Handle stupid cases */
  if (p == NULL) {
    sollyaFprintf(fd, "(null)");
  }

  /* Handle the general case.

     Can be optimized.

  */
  t = sparsePolynomialGetExpression(p,canonical);
  sollyaFprintf(fd, "%b",t);
  freeThing(t);
}

char *sparsePolynomialToString(sparse_polynomial_t p, int canonical) {
  node *t;
  char *str;
  char staticStr[8];
  int size, r;

  /* Handle stupid cases */
  if (p == NULL) return NULL;

  /* Handle the general case.

     Can be optimized.

  */
  t = sparsePolynomialGetExpression(p,canonical);
  size = sollya_snprintf(staticStr,8,"%b",t);
  if (size < 0) {
    freeThing(t);
    return NULL;
  }
  str = (char *) safeCalloc(size + 2, sizeof(char));
  r = sollya_snprintf(str,size,"%b",t);
  if (r < 0) {
    freeThing(t);
    safeFree(str);
    return NULL;
  }

  /* Return the string */
  return str;
}

static inline unsigned int sparsePolynomialGetMonomialCount(sparse_polynomial_t p) {
  
  /* Handle stupid input */
  if (p == NULL) return 0u;

  /* Return the number of monomials */
  return p->monomialCount;
}

static inline int sparsePolynomialCoefficientsAreDyadic(sparse_polynomial_t p, int defVal) {
  unsigned int i;
  int t;

  /* Handle stupid input */
  if (p == NULL) return 0;

  /* Handle the case when the polynomial has no monomials */
  if (p->monomialCount == 0u) return 1;

  /* Go over the coefficients and check if they are dyadic */
  for (i=0u;i<p->monomialCount;i++) {
    t = constantIsDyadic(p->coeffs[i], 99);
    if (t == 99) return defVal;
    if (!t) return 0;
  }
  return 1;
}

static inline int sparsePolynomialCoefficientsAreRational(sparse_polynomial_t p, int defVal) {
  unsigned int i;
  int t;

  /* Handle stupid input */
  if (p == NULL) return 0;

  /* Handle the case when the polynomial has no monomials */
  if (p->monomialCount == 0u) return 1;

  /* Go over the coefficients and check if they are rational */
  for (i=0u;i<p->monomialCount;i++) {
    t = constantIsRational(p->coeffs[i], 99);
    if (t == 99) return defVal;
    if (!t) return 0;
  }
  return 1;
}

static inline int sparsePolynomialCoefficientsHoldOnPrecBits(sparse_polynomial_t p, mp_prec_t prec, int defVal) {
  unsigned int i;
  int t;

  /* Handle stupid input */
  if (p == NULL) return 0;

  /* Handle the case when the polynomial has no monomials */
  if (p->monomialCount == 0u) return 1;

  /* Go over the coefficients and check if they hold on prec bits */
  for (i=0u;i<p->monomialCount;i++) {
    t = constantHoldsOnPrecBits(p->coeffs[i], prec, 99);
    if (t == 99) return defVal;
    if (!t) return 0;
  }
  return 1;
}

static inline sparse_polynomial_t sparsePolynomialRoundDyadic(sparse_polynomial_t p, mp_prec_t prec) { 
  sparse_polynomial_t res;
  unsigned int startSize, i;
  constant_t c;

  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* If the polynomial only has dyadic coefficients, just return a
     copy of it. 
  */
  if (sparsePolynomialCoefficientsAreDyadic(p, 0)) 
    return sparsePolynomialFromCopy(p);

  /* Here, at least one of the coefficients is not dyadic. Perform a
     "deep" copy while rounding the coefficients that are not
     dyadic. 
  */ 
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->hash.hasHash = 0;
  res->monomialCount = 0u;
  startSize = p->monomialCount;
  res->coeffs = (constant_t *) safeCalloc(startSize, sizeof(constant_t));
  res->monomialDegrees = (constant_t *) safeCalloc(startSize, sizeof(constant_t));
  for (i=0u;i<p->monomialCount;i++) {
    c = constantRoundDyadic(p->coeffs[i], prec);
    if (!constantIsZero(c, 0)) {
      res->coeffs[res->monomialCount] = c;
      res->monomialDegrees[res->monomialCount] = constantFromCopy(p->monomialDegrees[i]);
      (res->monomialCount)++;  
    } else {
      constantFree(c);
    }
  }

  /* If res->monomialCount still is zero, we never added anything
     because the two polynomials multiplied gave zero. We add a zero
     coefficient in this case. 
  */
  if (res->monomialCount == 0u) {
    res->coeffs[res->monomialCount] = constantFromInt(0);
    res->monomialDegrees[res->monomialCount] = constantFromInt(0);
    (res->monomialCount)++;
  }
  /* Set the degree of the polynomial */
  res->deg = constantFromCopy(res->monomialDegrees[res->monomialCount - 1u]);
  /* Adjust memory to the real amount of monomials */
  if (res->monomialCount != startSize) {
    res->coeffs = (constant_t *) safeRealloc(res->coeffs, 
					     ((size_t) (res->monomialCount)) * sizeof(constant_t));
    res->monomialDegrees = (constant_t *) safeRealloc(res->monomialDegrees, 
						      ((size_t) (res->monomialCount)) * sizeof(constant_t));
  }
  /* Adjust the degree for zero leading coefficients */
  __sparsePolynomialAdjustDegree(res);

  /* Return the polynomial */
  return res;
}

static inline sparse_polynomial_t sparsePolynomialRoundRational(sparse_polynomial_t p, mp_prec_t prec) {
  sparse_polynomial_t res;
  unsigned int startSize, i;
  constant_t c;

  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* If the polynomial only has rational coefficients, just return a
     copy of it. 
  */
  if (sparsePolynomialCoefficientsAreRational(p, 0)) 
    return sparsePolynomialFromCopy(p);

  /* Here, at least one of the coefficients is not rational. Perform a
     "deep" copy while rounding the coefficients that are not
     rational. 
  */ 
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->hash.hasHash = 0;
  res->monomialCount = 0u;
  startSize = p->monomialCount;
  res->coeffs = (constant_t *) safeCalloc(startSize, sizeof(constant_t));
  res->monomialDegrees = (constant_t *) safeCalloc(startSize, sizeof(constant_t));
  for (i=0u;i<p->monomialCount;i++) {
    c = constantRoundRational(p->coeffs[i], prec);
    if (!constantIsZero(c, 0)) {
      res->coeffs[res->monomialCount] = c;
      res->monomialDegrees[res->monomialCount] = constantFromCopy(p->monomialDegrees[i]);
      (res->monomialCount)++;  
    } else {
      constantFree(c);
    }
  }

  /* If res->monomialCount still is zero, we never added anything
     because the two polynomials multiplied gave zero. We add a zero
     coefficient in this case. 
  */
  if (res->monomialCount == 0u) {
    res->coeffs[res->monomialCount] = constantFromInt(0);
    res->monomialDegrees[res->monomialCount] = constantFromInt(0);
    (res->monomialCount)++;
  }
  /* Set the degree of the polynomial */
  res->deg = constantFromCopy(res->monomialDegrees[res->monomialCount - 1u]);
  /* Adjust memory to the real amount of monomials */
  if (res->monomialCount != startSize) {
    res->coeffs = (constant_t *) safeRealloc(res->coeffs, 
					     ((size_t) (res->monomialCount)) * sizeof(constant_t));
    res->monomialDegrees = (constant_t *) safeRealloc(res->monomialDegrees, 
						      ((size_t) (res->monomialCount)) * sizeof(constant_t));
  }
  /* Adjust the degree for zero leading coefficients */
  __sparsePolynomialAdjustDegree(res);

  /* Return the polynomial */
  return res;
}

static inline sparse_polynomial_t sparsePolynomialRound(sparse_polynomial_t p, mp_prec_t prec) { 
  sparse_polynomial_t res;
  unsigned int startSize, i;
  constant_t c;

  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* If the polynomial only has coefficients that hold on prec bits,
     just return a copy of it.
  */
  if (sparsePolynomialCoefficientsHoldOnPrecBits(p, prec, 0)) 
    return sparsePolynomialFromCopy(p);

  /* Here, at least one of the coefficients does not hold on prec
     bits. Perform a "deep" copy while rounding the coefficients that
     are not dyadic.
  */ 
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->hash.hasHash = 0;
  res->monomialCount = 0u;
  startSize = p->monomialCount;
  res->coeffs = (constant_t *) safeCalloc(startSize, sizeof(constant_t));
  res->monomialDegrees = (constant_t *) safeCalloc(startSize, sizeof(constant_t));
  for (i=0u;i<p->monomialCount;i++) {
    c = constantRound(p->coeffs[i], prec);
    if (!constantIsZero(c, 0)) {
      res->coeffs[res->monomialCount] = c;
      res->monomialDegrees[res->monomialCount] = constantFromCopy(p->monomialDegrees[i]);
      (res->monomialCount)++;  
    } else {
      constantFree(c);
    }
  }

  /* If res->monomialCount still is zero, we never added anything
     because the two polynomials multiplied gave zero. We add a zero
     coefficient in this case. 
  */
  if (res->monomialCount == 0u) {
    res->coeffs[res->monomialCount] = constantFromInt(0);
    res->monomialDegrees[res->monomialCount] = constantFromInt(0);
    (res->monomialCount)++;
  }
  /* Set the degree of the polynomial */
  res->deg = constantFromCopy(res->monomialDegrees[res->monomialCount - 1u]);
  /* Adjust memory to the real amount of monomials */
  if (res->monomialCount != startSize) {
    res->coeffs = (constant_t *) safeRealloc(res->coeffs, 
					     ((size_t) (res->monomialCount)) * sizeof(constant_t));
    res->monomialDegrees = (constant_t *) safeRealloc(res->monomialDegrees, 
						      ((size_t) (res->monomialCount)) * sizeof(constant_t));
  }
  /* Adjust the degree for zero leading coefficients */
  __sparsePolynomialAdjustDegree(res);

  /* Return the polynomial */
  return res;
}

static inline unsigned int sparsePolynomialGetReferenceCount(sparse_polynomial_t p) {
  if (p == NULL) return 0u;
  return p->refCount;
}

static inline int sparsePolynomialUsesExpressionConstant(sparse_polynomial_t p) {
  unsigned int i;
  
  /* Handle stupid input */
  if (p == NULL) return 0;

  /* Handle the case when the polynomial has no monomials */
  if (p->monomialCount == 0u) return 0;

  /* Go over all monomials and check the constants */
  for (i=0u;i<p->monomialCount;i++) {
    if (constantUsesExpressionConstant(p->coeffs[i])) return 1;
    if (constantUsesExpressionConstant(p->monomialDegrees[i])) return 1;
  }
  if (constantUsesExpressionConstant(p->deg)) return 1;

  /* The expression has nowhere been found */
  return 0;
}

static inline int sparsePolynomialReferencesExpression(sparse_polynomial_t p, node *expr) {
  unsigned int i;
  
  /* Handle stupid input */
  if (p == NULL) return 0;

  /* If the polynomial does not use any expression, it cannot
     reference the given expression 
  */
  if (!sparsePolynomialUsesExpressionConstant(p)) {
    return 0;
  }

  /* Handle the case when the polynomial has no monomials */
  if (p->monomialCount == 0u) return 0;

  /* Go over all monomials and check the constants */
  for (i=0u;i<p->monomialCount;i++) {
    if (constantReferencesExpression(p->coeffs[i], expr)) return 1;
    if (constantReferencesExpression(p->monomialDegrees[i], expr)) return 1;
  }
  if (constantReferencesExpression(p->deg, expr)) return 1;

  /* The expression has nowhere been found */
  return 0;
}

static inline uint64_t sparsePolynomialHash(sparse_polynomial_t p) {
  uint64_t hash;
  unsigned int i;
  
  /* Handle stupid input */
  if (p == NULL) return hashPointer(NULL);

  /* Check if the hash has already been computed */
  if (p->hash.hasHash) {
    return p->hash.hash;
  }
  
  /* Compute the hash */
  hash = hashUnsignedInt(p->monomialCount);
  hash = hashCombine(hash, constantHash(p->deg));
  for (i=0u;i<p->monomialCount;i++) {
    hash = hashCombine(hash, constantHash(p->coeffs[i]));
    hash = hashCombine(hash, constantHash(p->monomialDegrees[i]));
  }

  /* Cache the result */
  p->hash.hash = hash;
  p->hash.hasHash = 1;
  
  /* Return the hash */
  return hash;
}

/* End of part for sparse polynomials */

/* Start of part for general (composed) polynomials */

static inline void __polynomialSparsify(polynomial_t);
int polynomialUsesExpressionConstant(polynomial_t p);

static inline polynomial_t __polynomialAllocate() {
  return (polynomial_t) safeMalloc(sizeof(struct __polynomial_struct_t));
}

static inline void __polynomialFreeMem(polynomial_t p) {
  safeFree(p);
}

void polynomialInitializeCaches() {
  sparsePolynomialInitializeCaches();
}

void polynomialFreeCaches() {
  sparsePolynomialFreeCaches();
}

static inline polynomial_t __polynomialBuildFromSparse(sparse_polynomial_t p) {
  polynomial_t res;

  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* Handle general case */
  res = __polynomialAllocate();
  res->refCount = 1u;
  res->type = SPARSE;
  res->outputType = ANY_FORM;
  res->value.sparse = p;
  res->hash.hasHash = 0;
  res->usesExpressionConstant.cached = 0;

  /* Return the newly built polynomial */
  return res;
}

static inline sparse_polynomial_t __polynomialSparsifyComposition(polynomial_t p, 
								  sparse_polynomial_t q) {
  sparse_polynomial_t g, h, r;

  /* Handle stupid input */
  if (p == NULL) return NULL;
  if (q == NULL) return NULL;

  /* Consider the structure of p */
  switch (p->type) {
  case SPARSE:
    return sparsePolynomialCompose(p->value.sparse, q);
    break;
  case ADDITION:
    g = __polynomialSparsifyComposition(p->value.pair.g, q);
    h = __polynomialSparsifyComposition(p->value.pair.h, q);
    r = sparsePolynomialAdd(g, h);
    sparsePolynomialFree(g);
    sparsePolynomialFree(h);
    return r;
    break;
  case SUBTRACTION:
    g = __polynomialSparsifyComposition(p->value.pair.g, q);
    h = __polynomialSparsifyComposition(p->value.pair.h, q);
    r = sparsePolynomialSub(g, h);
    sparsePolynomialFree(g);
    sparsePolynomialFree(h);
    return r;
    break;
  case MULTIPLICATION:
    g = __polynomialSparsifyComposition(p->value.pair.g, q);
    h = __polynomialSparsifyComposition(p->value.pair.h, q);
    r = sparsePolynomialMul(g, h);
    sparsePolynomialFree(g);
    sparsePolynomialFree(h);
    return r;
    break;
  case COMPOSITION:
    h = __polynomialSparsifyComposition(p->value.pair.h, 
					q);
    r = __polynomialSparsifyComposition(p->value.pair.g, h);
    sparsePolynomialFree(h);
    return r;
    break;
  case NEGATE:
    g = __polynomialSparsifyComposition(p->value.g, q);
    r = sparsePolynomialNeg(g);
    sparsePolynomialFree(g);
    return r;
    break;
  case POWER:
    g = __polynomialSparsifyComposition(p->value.powering.g, q);
    if (!sparsePolynomialPowConstant(&r, g, p->value.powering.c)) {
      sollyaFprintf(stderr,"Error: __polynomialSparsifyComposition: could not compute power of sparse polynomial\n");
      exit(1);
    }
    sparsePolynomialFree(g);
    return r;
    break;
  default:
    break;
  }
  __polynomialSparsify(p);
  return sparsePolynomialCompose(p->value.sparse, q);
}

static inline void __polynomialSparsify(polynomial_t p) {
  sparse_polynomial_t sp;

  /* Handle stupid input */
  if (p == NULL) return;
  
  /* General case */
  switch (p->type) {
  case SPARSE:
    /* Nothing to do */
    return;
    break;
  case ADDITION:
    __polynomialSparsify(p->value.pair.g);
    __polynomialSparsify(p->value.pair.h);
    sp = sparsePolynomialAdd(p->value.pair.g->value.sparse,
			     p->value.pair.h->value.sparse);
    polynomialFree(p->value.pair.g);
    polynomialFree(p->value.pair.h);
    p->value.sparse = sp;
    break;
  case SUBTRACTION:
    __polynomialSparsify(p->value.pair.g);
    __polynomialSparsify(p->value.pair.h);
    sp = sparsePolynomialSub(p->value.pair.g->value.sparse,
			     p->value.pair.h->value.sparse);
    polynomialFree(p->value.pair.g);
    polynomialFree(p->value.pair.h);
    p->value.sparse = sp;
    break;
  case MULTIPLICATION:
    __polynomialSparsify(p->value.pair.g);
    __polynomialSparsify(p->value.pair.h);
    sp = sparsePolynomialMul(p->value.pair.g->value.sparse,
			     p->value.pair.h->value.sparse);
    polynomialFree(p->value.pair.g);
    polynomialFree(p->value.pair.h);
    p->value.sparse = sp;
    break;
  case COMPOSITION:
    __polynomialSparsify(p->value.pair.h);
    sp = __polynomialSparsifyComposition(p->value.pair.g, 
					 p->value.pair.h->value.sparse);
    polynomialFree(p->value.pair.g);
    polynomialFree(p->value.pair.h);
    p->value.sparse = sp;
    break;
  case NEGATE:
    __polynomialSparsify(p->value.g);
    sp = sparsePolynomialNeg(p->value.g->value.sparse);
    polynomialFree(p->value.g);
    p->value.sparse = sp;
    break;
  case POWER:
    __polynomialSparsify(p->value.powering.g);
    if (!sparsePolynomialPowConstant(&sp, 
				     p->value.powering.g->value.sparse,
				     p->value.powering.c)) {
      sollyaFprintf(stderr,"Error: __polynomialSparsify: could not compute power of sparse polynomial\n");
      exit(1);
    }
    polynomialFree(p->value.powering.g);
    constantFree(p->value.powering.c);
    p->value.sparse = sp;
    break;
  }
  if (p->usesExpressionConstant.cached) {
    if (p->usesExpressionConstant.res) {
      p->usesExpressionConstant.cached = 0;
    }
  }
  p->type = SPARSE;
}

static inline polynomial_t __polynomialExecuteCompositionCompose(polynomial_t p, polynomial_t q) {
  polynomial_t res, t;

  /* Make compiler happy */
  res = NULL;
  /* End of compiler happiness */

  /* Handle stupid input */
  if (p == NULL) return NULL;
  if (q == NULL) return NULL;

  /* General case */
  switch (p->type) {
  case SPARSE:
    __polynomialSparsify(q);
    res = __polynomialAllocate();
    res->refCount = 1u;
    res->hash.hasHash = 0;
    res->usesExpressionConstant.cached = 0;
    res->outputType = ANY_FORM;
    res->value.sparse = sparsePolynomialCompose(p->value.sparse,  
						q->value.sparse);
    res->type = SPARSE;
    break;
  case ADDITION:
  case SUBTRACTION:
  case MULTIPLICATION:
    res = __polynomialAllocate();
    res->refCount = 1u;
    res->hash.hasHash = 0;
    res->usesExpressionConstant.cached = 0;
    res->outputType = ANY_FORM;
    res->type = p->type;
    res->value.pair.g = __polynomialExecuteCompositionCompose(p->value.pair.g, q);
    res->value.pair.h = __polynomialExecuteCompositionCompose(p->value.pair.h, q);
    break;
  case COMPOSITION:
    t = __polynomialExecuteCompositionCompose(p->value.pair.h, q);
    res = __polynomialExecuteCompositionCompose(p->value.pair.g, t);
    polynomialFree(t);
    break;
  case NEGATE:
    res = __polynomialAllocate();
    res->refCount = 1u;
    res->hash.hasHash = 0;
    res->usesExpressionConstant.cached = 0;
    res->outputType = ANY_FORM;
    res->type = NEGATE;
    res->value.g = __polynomialExecuteCompositionCompose(p->value.g, q);
    break;
  case POWER:
    res = __polynomialAllocate();
    res->refCount = 1u;
    res->hash.hasHash = 0;
    res->usesExpressionConstant.cached = 0;
    res->outputType = ANY_FORM;
    res->type = POWER;
    res->value.powering.g = __polynomialExecuteCompositionCompose(p->value.powering.g, q);
    res->value.powering.c = constantFromCopy(p->value.powering.c);
    break;
  }  

  return res;
}

static inline void __polynomialExecuteComposition(polynomial_t p) {
  sparse_polynomial_t sp;
  constant_t c;
  polynomial_type_t t;
  polynomial_t g, h;

  /* Handle stupid input */
  if (p == NULL) return;
  
  /* General case */
  switch (p->type) {
  case SPARSE:
    /* Nothing to do */
    return;
    break;
  case ADDITION:
  case SUBTRACTION:
  case MULTIPLICATION:
    __polynomialExecuteComposition(p->value.pair.g);
    __polynomialExecuteComposition(p->value.pair.h);
    break;
  case COMPOSITION:
    __polynomialExecuteComposition(p->value.pair.g);
    __polynomialExecuteComposition(p->value.pair.h);
    t = p->value.pair.g->type;
    switch (t) {
    case SPARSE:
      __polynomialSparsify(p->value.pair.h);
      sp = sparsePolynomialCompose(p->value.pair.g->value.sparse, 
				   p->value.pair.h->value.sparse);
      polynomialFree(p->value.pair.g);
      polynomialFree(p->value.pair.h);
      p->value.sparse = sp;
      p->type = SPARSE;
      return;
      break;
    case ADDITION:
    case SUBTRACTION:
    case MULTIPLICATION:
      g = __polynomialExecuteCompositionCompose(p->value.pair.g->value.pair.g, 
						p->value.pair.h);
      h = __polynomialExecuteCompositionCompose(p->value.pair.g->value.pair.h, 
						p->value.pair.h);
      polynomialFree(p->value.pair.g);
      polynomialFree(p->value.pair.h);
      p->value.pair.g = g;
      p->value.pair.h = h;
      p->type = t;
      return;
      break;
    case COMPOSITION:
      /* Should never happen */
      __polynomialSparsify(p);
      return;
      break;
    case NEGATE:
      g = __polynomialExecuteCompositionCompose(p->value.pair.g->value.g, 
						p->value.pair.h);
      polynomialFree(p->value.pair.g);
      polynomialFree(p->value.pair.h);
      p->value.g = g;
      p->type = NEGATE;
      return;
      break;
    case POWER:
      g = __polynomialExecuteCompositionCompose(p->value.pair.g->value.powering.g, 
						p->value.pair.h);
      c = constantFromCopy(p->value.pair.g->value.powering.c);
      polynomialFree(p->value.pair.g);
      polynomialFree(p->value.pair.h);
      p->value.powering.g = g;
      p->value.powering.c = c;
      p->type = POWER;
      return;
      break;
    }
    break;
  case NEGATE:
    __polynomialExecuteComposition(p->value.g);
    break;
  case POWER:
    __polynomialExecuteComposition(p->value.powering.g);
    break;
  }
}

polynomial_t polynomialFromMpfrConstant(mpfr_t c) {
  return __polynomialBuildFromSparse(sparsePolynomialFromMpfrConstant(c));
}

polynomial_t polynomialFromMpzConstant(mpz_t c) {
  return __polynomialBuildFromSparse(sparsePolynomialFromMpzConstant(c));
}

polynomial_t polynomialFromMpqConstant(mpq_t c) {
  return __polynomialBuildFromSparse(sparsePolynomialFromMpqConstant(c));
}

polynomial_t polynomialFromIntConstant(int c) {
  return __polynomialBuildFromSparse(sparsePolynomialFromIntConstant(c));
}

polynomial_t polynomialFromIdentity() {
  return __polynomialBuildFromSparse(sparsePolynomialFromIdentity());
}

polynomial_t polynomialFromMpfrCoefficients(mpfr_t *coeffs, unsigned int deg) {
  return __polynomialBuildFromSparse(sparsePolynomialFromMpfrCoefficients(coeffs, deg));
}

static inline polynomial_t __polynomialFromConstant(constant_t c) {
  return __polynomialBuildFromSparse(sparsePolynomialFromConstant(c));
}

int polynomialFromConstantExpressionCoefficients(polynomial_t *r, node **coeffs, unsigned int deg) {
  sparse_polynomial_t sp;

  if (!sparsePolynomialFromConstantExpressionCoefficients(&sp, coeffs, deg)) return 0;
  *r = __polynomialBuildFromSparse(sp);
  return 1;
}

polynomial_t polynomialFromCopy(polynomial_t p) {
  if (p == NULL) return NULL;
  p->refCount++;
  return p;
}

void polynomialFree(polynomial_t p) {
  if (p == NULL) return;
  p->refCount--;
  if (p->refCount > 0u) return;
  switch (p->type) {
  case SPARSE:
    sparsePolynomialFree(p->value.sparse);
    break;
  case ADDITION:
  case SUBTRACTION:
  case MULTIPLICATION:
  case COMPOSITION:
    polynomialFree(p->value.pair.g);
    polynomialFree(p->value.pair.h);
    break;
  case NEGATE:
    polynomialFree(p->value.g);
    break;
  case POWER:
    polynomialFree(p->value.powering.g);
    constantFree(p->value.powering.c);
    break;
  }
  __polynomialFreeMem(p);
}

static inline int __polynomialGetDegreeAsIntCheap(polynomial_t p) {
  int gd, hd, res;
  unsigned int u, gdu, tu, ttu;

  /* Handle stupid input */
  if (p == NULL) return -1;

  /* Try to handle the case without sparsifying the polynomial. */
  switch (p->type) {
  case SPARSE:
    return sparsePolynomialGetDegreeAsInt(p->value.sparse);
    break;
  case ADDITION:
  case SUBTRACTION:
    gd = __polynomialGetDegreeAsIntCheap(p->value.pair.g);
    hd = __polynomialGetDegreeAsIntCheap(p->value.pair.h);
    if ((gd == 0) && (hd == 0)) return 0;
    if ((gd >= 0) && (hd >= 0) && (gd != hd)) 
      return (hd > gd ? hd : gd);
    break;
  case MULTIPLICATION:
    gd = __polynomialGetDegreeAsIntCheap(p->value.pair.g);
    hd = __polynomialGetDegreeAsIntCheap(p->value.pair.h);
    if ((gd != 0) && (hd != 0)) {
      if ((gd < 0) || (hd < 0)) return -1;
      if (tryExactIntAddition(&res, gd, hd)) {
	return res;
      } else {
	return -1;
      }
    }
    break;
  case COMPOSITION:
    gd = __polynomialGetDegreeAsIntCheap(p->value.pair.g);
    hd = __polynomialGetDegreeAsIntCheap(p->value.pair.h);
    if ((gd == 0) || (hd == 0)) return 0;
    if ((gd < 0) || (hd < 0)) return -1;
    if (tryExactIntMultiplication(&res, gd, hd)) {
      return res;
    } else {
      return -1;
    }    
    break;
  case NEGATE:
    return __polynomialGetDegreeAsIntCheap(p->value.g);
    break;
  case POWER:
    gd = __polynomialGetDegreeAsIntCheap(p->value.powering.g);
    if (gd == 0) return 0;
    if (tryConstantToUnsignedInt(&u, p->value.powering.c)) {
      if (u == 0u) {
	return 0;
      } else {
	if (gd < 0) return -1;
	gdu = gd;
	if (tryExactUnsignedIntMultiplication(&tu,gdu,u)) {
	  res = tu;
	  if (res < 0) return -1;
	  ttu = res;
	  if (ttu != tu) return -1;
	  return res;
	} else {
	  return -1;
	}
      }
    } else {
      return -1;
    }
    break;
  }

  /* The fallback would be too expensive */
  return -1;
}

static inline int __polynomialEqualCheap(polynomial_t p, polynomial_t q) {
  if (p == NULL) return 0;
  if (q == NULL) return 0;
  if (p == q) return 1;
  if (p->type != q->type) return 0;
  switch (p->type) {
  case SPARSE:
    return sparsePolynomialEqual(p->value.sparse, q->value.sparse, 0);
    break;
  case ADDITION:
  case SUBTRACTION:
  case MULTIPLICATION:
  case COMPOSITION:
    return (__polynomialEqualCheap(p->value.pair.g, q->value.pair.g) &&
	    __polynomialEqualCheap(p->value.pair.h, q->value.pair.h));
    break;
  case NEGATE:
    return __polynomialEqualCheap(p->value.g, q->value.g);
    break;
  case POWER:
    return (__polynomialEqualCheap(p->value.powering.g, q->value.powering.g) &&
	    constantIsEqual(p->value.powering.c, q->value.powering.c, 0));    
    break;
  }
  return 0;
}

static inline int __polynomialStructurallyEqualCheap(polynomial_t p, polynomial_t q) {
  if (p == NULL) return 0;
  if (q == NULL) return 0;
  if (p == q) return 1;
  if (p->type != q->type) return 0;
  if (p->outputType != q->outputType) return 0;
  switch (p->type) {
  case SPARSE:
    return sparsePolynomialEqual(p->value.sparse, q->value.sparse, 0);
    break;
  case ADDITION:
  case SUBTRACTION:
  case MULTIPLICATION:
  case COMPOSITION:
    return (__polynomialStructurallyEqualCheap(p->value.pair.g, q->value.pair.g) &&
	    __polynomialStructurallyEqualCheap(p->value.pair.h, q->value.pair.h));
    break;
  case NEGATE:
    return __polynomialStructurallyEqualCheap(p->value.g, q->value.g);
    break;
  case POWER:
    return (__polynomialStructurallyEqualCheap(p->value.powering.g, q->value.powering.g) &&
	    constantIsEqual(p->value.powering.c, q->value.powering.c, 0));    
    break;
  }
  return 0;
}

static inline int __polynomialIsConstantCheap(polynomial_t p) {
  if (p == NULL) return 0;
  switch (p->type) {
  case SPARSE:
    return sparsePolynomialIsConstant(p->value.sparse, 0);
    break;
  case ADDITION:
  case SUBTRACTION:
  case MULTIPLICATION:
    return (__polynomialIsConstantCheap(p->value.pair.g) &&
	    __polynomialIsConstantCheap(p->value.pair.h));
    break;
  case COMPOSITION:
    return (__polynomialIsConstantCheap(p->value.pair.g) ||
	    __polynomialIsConstantCheap(p->value.pair.h));
    break;
  case NEGATE:
    return __polynomialIsConstantCheap(p->value.g);
    break;
  case POWER:
    return (__polynomialIsConstantCheap(p->value.powering.g) ||
	    constantIsZero(p->value.powering.c, 0));
    break;
  }
  return 0;
}

static inline void __polynomialUnifyEqual(polynomial_t p, polynomial_t q) {
  if (p == NULL) return;
  if (q == NULL) return;
  if (p == q) return;
  if ((p->type == SPARSE) && (q->type == SPARSE)) {
    if (p->value.sparse == q->value.sparse) return;
    if (sparsePolynomialGetReferenceCount(p->value.sparse) > 
	sparsePolynomialGetReferenceCount(q->value.sparse)) {
      sparsePolynomialFree(q->value.sparse);
      q->value.sparse = sparsePolynomialFromCopy(p->value.sparse);
    } else {
      sparsePolynomialFree(p->value.sparse);
      p->value.sparse = sparsePolynomialFromCopy(q->value.sparse);
    }
    return;
  }
  if ((p->type == SPARSE) && (q->type != SPARSE)) {
    switch (q->type) {
    case ADDITION:
    case SUBTRACTION:
    case MULTIPLICATION:
    case COMPOSITION:
      polynomialFree(q->value.pair.g);
      polynomialFree(q->value.pair.h);
      break;
    case NEGATE:
      polynomialFree(q->value.g);
      break;
    case POWER:
      polynomialFree(q->value.powering.g);
      constantFree(q->value.powering.c);
      break;
    default:
      return;
    }
    q->type = SPARSE;
    q->value.sparse = sparsePolynomialFromCopy(p->value.sparse);
    return;
  }
  if ((q->type == SPARSE) && (p->type != SPARSE)) {
    switch (p->type) {
    case ADDITION:
    case SUBTRACTION:
    case MULTIPLICATION:
    case COMPOSITION:
      polynomialFree(p->value.pair.g);
      polynomialFree(p->value.pair.h);
      break;
    case NEGATE:
      polynomialFree(p->value.g);
      break;
    case POWER:
      polynomialFree(p->value.powering.g);
      constantFree(p->value.powering.c);
      break;
    default:
      return;
    }
    p->type = SPARSE;
    p->value.sparse = sparsePolynomialFromCopy(q->value.sparse);
    return;
  }
}

int polynomialEqual(polynomial_t p, polynomial_t q, int defVal) {
  int dp, dq, res;

  /* Handle stupid inputs */
  if (p == NULL) return defVal;
  if (q == NULL) return defVal;

  /* Pointer equality */
  if (p == q) return 1;

  /* If both polynomials are in sparse form, just use this */
  if ((p->type == SPARSE) && (q->type == SPARSE)) {
    res = sparsePolynomialEqual(p->value.sparse, q->value.sparse, 51);
    if (res == 51) return defVal;
    if (!res) return 0;
    __polynomialUnifyEqual(p,q);
    return 1;
  }

  /* Try to cheaply compute the degrees. If both degrees can be
     computed but are different, we know that the polynomials are
     different.
  */
  dp = __polynomialGetDegreeAsIntCheap(p);
  dq = __polynomialGetDegreeAsIntCheap(q);
  if ((dp >= 0) && (dq >= 0) && (dp != dq)) return 0;

  /* If p and q are written the same way and each sub-polynomial is
     the same, they are the same. 
  */
  if (__polynomialEqualCheap(p, q)) {
    __polynomialUnifyEqual(p, q);
    return 1;
  }

  /* General case
     
     Can perhaps be optimized still a little bit.
  
  */
  __polynomialSparsify(p);
  __polynomialSparsify(q);
  res = sparsePolynomialEqual(p->value.sparse, q->value.sparse, 51);
  if (res == 51) return defVal;
  if (!res) return 0;
  __polynomialUnifyEqual(p, q);
  return 1;
}

unsigned int polynomialGetNumberOfNonZeroCoefficients(polynomial_t p) {

  /* Handle stupid inputs */
  if (p == NULL) return 0u;

  /* General case */
  __polynomialSparsify(p);
  return sparsePolynomialGetNumberOfNonZeroCoefficients(p->value.sparse);
}

int polynomialStructurallyEqual(polynomial_t p, polynomial_t q, int defVal) {
  constant_t c0, cd;
  mpz_t deg;
  unsigned int nzc;
  
  /* Handle stupid inputs */
  if (p == NULL) return defVal;
  if (q == NULL) return defVal;

  /* Pointer equality */
  if (p == q) return 1;

  /* Use hashes if we have them already computed */
  if (p->hash.hasHash && q->hash.hasHash) {
    /* Both polynomials have a hash.
       
       If the hashes are different, they are structurally different.
       
    */
    if (p->hash.hash != q->hash.hash) return 0;
  }

  /* Do a first check using a cheap function */
  if (__polynomialStructurallyEqualCheap(p, q)) return 1;
  
  /* If the polynomials are not mathematically equal, they cannot be
     structurally equal.
  */
  if (!polynomialEqual(p, q, 1)) return 0;

  /* Try with the other default value */
  if (!polynomialEqual(p, q, 0)) {
	__polynomialSparsify(p);
	__polynomialSparsify(q);    
  }

  /* Retry */
  if (!polynomialEqual(p, q, 1)) return 0;
  
  /* Here the polynomials are mathematically equal or their
     mathematical equality cannot be decided.

     Start by assuring us that we know that they are mathematically
     equal.
     
  */
  if (polynomialEqual(p, q, 0)) {
    /* Here, the polynomials are mathematically equal. 

       If they have the same output form and that output form 
       is canonical or hornerized, they are structurally equal.
       
       If they have the same output form and that output form is
       ANY_FORM, we sparsify both polynomials and can then answer that
       they are structurally equal.

    */
    if (p->outputType == q->outputType) {
      if (p->outputType == ANY_FORM) {
	__polynomialSparsify(p);
	__polynomialSparsify(q);
      }
      return 1;
    }

    /* Here, the polynomials have different output forms.

       Start by checking if one of them is constant. As the other one
       is mathematically equal, it is constant, too.

       If the polynomials are constant, they are structurally the
       same.

    */
    if (polynomialIsConstant(p, 0)) return 1;

    /* Now check if the polynomials (testing one is sufficient as they
       are mathematically equal) have a non-zero constant coefficient, 
       a leading coefficient equal to 1 and only two non-zero coefficients.
       
       In this particular case the hornerized form c_0 + x^d is
       structurally the same as the canonical form c_0 + x^d.
    */
    c0 = __polynomialGetIthCoefficientAsConstantIntIndex(p, 0);
    if (!constantIsZero(c0, 1)) {
      mpz_init(deg);
      polynomialGetDegree(deg, p);
      cd = __polynomialGetIthCoefficientAsConstant(p, deg);
      mpz_clear(deg);
      if (constantIsOne(cd, 0)) {
	nzc = polynomialGetNumberOfNonZeroCoefficients(p);
	if (nzc == 2u) {
	  constantFree(cd);
	  constantFree(c0);
	  return 1;
	}
      }
      constantFree(cd);
    }
    constantFree(c0);

    /* Here, the polynomials have different output forms and are not
       constant.

       If none of these output forms is ANY_FORM, they are
       structurally different.
    */
    if ((p->outputType != ANY_FORM) &&
	(q->outputType != ANY_FORM)) {
      return 0;
    }

    /* Here at least one of the output forms is ANY_FORM and
       the other output form is not ANY_FORM.

       We sparsify both polynomials. Then they are structurally equal
       iff none of them is in CANONICAL_FORM as the output for
       sparsified ANY_FORMS is HORNER_FORM.

    */
    __polynomialSparsify(p);
    __polynomialSparsify(q);
    if ((p->outputType != CANONICAL_FORM) &&
	(q->outputType != CANONICAL_FORM)) {
      return 1;
    }
    
    /* Here the polynomials are structurally different */
    return 0;
  }

  /* Here, we don't know if the polynomials are 
     mathematically the same and we already retried. 

     We return the default answer.

  */
  return defVal;
}

int polynomialIsIdentity(polynomial_t p, int defVal) {
  int dp;

  /* Handle stupid inputs */
  if (p == NULL) return defVal;

  /* If the polynomial is sparse, just use the test function on sparse
     polynomials 
  */
  if (p->type == SPARSE) 
    return sparsePolynomialIsIdentity(p->value.sparse, defVal);
  
  /* If the polynomial is constant, it is not the identity function */
  if (__polynomialIsConstantCheap(p)) return 0;

  /* If the degree of the polynomial can be computed easily and it is
     not 1, the polynomial can't be the identity function. 
  */
  dp = __polynomialGetDegreeAsIntCheap(p);
  if ((dp >= 0) && (dp != 1)) return 0;

  /* Some optimisations are still possible */

  /* Sparsify the polynomial and use the test function on sparse
     polynomials 
  */
  __polynomialSparsify(p);
  return sparsePolynomialIsIdentity(p->value.sparse, defVal);  
}

int polynomialIsConstant(polynomial_t p, int defVal) {
  int deg;
  constant_t c;
  
  /* Handle stupid inputs */
  if (p == NULL) return defVal;

  /* If the polynomial is sparse, just use the test function on sparse
     polynomials.
  */
  if (p->type == SPARSE) 
    return sparsePolynomialIsConstant(p->value.sparse, defVal);
  
  /* Try a cheap function to determine if the polynomial is
     constant 
  */
  if (__polynomialIsConstantCheap(p)) return 1;

  /* Here, we still do not know if the polynomial is constant.

     We compute its degree.

     If the degree is larger than the largest machine integer, the
     function will return a negative value instead of the degree. 

     In this case, we sparsify the polynomial and call the test 
     function on the sparsified polynomial.

     Otherwise, 

     * if the degree is zero, we are sure the polynomial is 
       constant,

     * otherwise, we get the coefficient corresponding to the 
       degree and check if it is non-zero. If we cannot show it 
       is non-zero, we return the default answer. Otherwise, we
       are sure that the polynomial is not constant.

  */
  deg = polynomialGetDegreeAsInt(p);
  if (deg < 0) {
    __polynomialSparsify(p);    
    return sparsePolynomialIsConstant(p->value.sparse, defVal);
  }

  /* Here the degree of the polynomial holds on a machine integer. 

     If it is zero, the polynomial is constant.
  */
  if (deg == 0) return 1;

  /* Get the coefficient of the degree of the polynomial */
  c = __polynomialGetIthCoefficientAsConstantIntIndex(p, deg);

  /* Check if coefficient of degree is zero. */
  if (constantIsZero(c, 1)) {
    /* Here the coefficient of the polynomial corresponding 
       to its degree is said to be zero. This means it cannot
       be shown to be non-zero.

       We return the default value, as we do not know if the 
       polynomial is constant or not.

    */
    constantFree(c);
    return defVal;
  }

  /* Free the constant corresponding to coefficient */
  constantFree(c);
  
  /* Here the coefficient of degree deg != 0 is shown to be
     non-zero. Hence, the polynomial is not constant.
  */
  return 0;
}

void polynomialDiv(polynomial_t *quot, polynomial_t *rest, polynomial_t a, polynomial_t b) {
  sparse_polynomial_t sq, sr;
  polynomial_t q, r, qg, qh, rg, rh, t, rcp;
  constant_t one, tc, c, recprC;

  /* Handle stupid inputs */
  if ((a == NULL) || (b == NULL)) {
    *quot = NULL;
    *rest = NULL;
    return;
  }

  /* If a and b are both sparse, just use the sparse division */
  if ((a->type == SPARSE) && (b->type == SPARSE)) {
    sparsePolynomialDiv(&sq, &sr, a->value.sparse, b->value.sparse);

    *quot = __polynomialBuildFromSparse(sq);
    *rest = __polynomialBuildFromSparse(sr);
    return;
  }

  /* If b is constant, do some special handling */
  if (polynomialGetDegreeAsInt(b) == 0) {
    /* Get the constant corresponding to b */
    c = __polynomialGetIthCoefficientAsConstantIntIndex(b, 0);
    /* If we are sure that b = c is non-zero, we are sure we can
       compute 1/c and replace the division by a multiplication 
    */
    if (!constantIsZero(c, 1)) {
      /* Compute 1/c and represent it as a polynomial */
      one = constantFromInt(1);
      recprC = constantDiv(one, c);
      constantFree(one);
      rcp = __polynomialFromConstant(recprC);
      constantFree(recprC);
      *quot = polynomialMul(rcp, a);
      *rest = polynomialFromIntConstant(0);
      polynomialFree(rcp);
    } else {
      /* Here b = c = 0. We take q = 0 and r = a */
      *quot = polynomialFromIntConstant(0);
      *rest = polynomialFromCopy(a);
    }
    constantFree(c);
    return;
  }

  /* If a and b are equal, the answer is easy */
  if (__polynomialEqualCheap(a, b)) {
    *quot = polynomialFromIntConstant(1);
    *rest = polynomialFromIntConstant(0);
    return;
  }

  /* If a is a negation, do the division recursively and negate the
     quotient and rest. 
  */
  if (a->type == NEGATE) {
    polynomialDiv(&q, &r, a->value.g, b);
    *quot = polynomialNeg(q);
    *rest = polynomialNeg(r);
    polynomialFree(q);
    polynomialFree(r);
    return;
  }

  /* If a is an addition, just do the division
     recursively.
  */
  if (a->type == ADDITION) {
    polynomialDiv(&qg, &rg, a->value.pair.g, b);
    polynomialDiv(&qh, &rh, a->value.pair.h, b);
    *quot = polynomialAdd(qg, qh);
    *rest = polynomialAdd(rg, rh);
    polynomialFree(qg);
    polynomialFree(rg);
    return;
  }

  /* If a is a subtraction, just do the division
     recursively.
  */
  if (a->type == SUBTRACTION) {
    polynomialDiv(&qg, &rg, a->value.pair.g, b);
    polynomialDiv(&qh, &rh, a->value.pair.h, b);
    *quot = polynomialSub(qg, qh);
    *rest = polynomialSub(rg, rh);
    polynomialFree(qg);
    polynomialFree(rg);
    return;
  }

  /* If a is of the form a = p^c with c >= 2 and b = p, 
     the result is q = p^(c - 1) and r = 0.
  */
  if (a->type == POWER) {
    if ((!constantIsZero(a->value.powering.c,1)) &&
	(!constantIsOne(a->value.powering.c,1))) {
      if (polynomialEqual(a->value.powering.g, b, 0)) {
	one = constantFromInt(1);
	tc = constantSub(a->value.powering.c, one);
	constantFree(one);
	t = __polynomialFromConstant(tc);
	constantFree(tc);
	if (polynomialPow(quot, b, t)) {
	  *rest = polynomialFromIntConstant(0);
	  polynomialFree(t);
	  return;
	}
	polynomialFree(t);
      }
    }
  }

  /* General case 

     Some cases can still be optimized.

  */
  __polynomialSparsify(a);
  __polynomialSparsify(b);
  sparsePolynomialDiv(&sq, &sr, a->value.sparse, b->value.sparse);

  *quot = __polynomialBuildFromSparse(sq);
  *rest = __polynomialBuildFromSparse(sr);
}

polynomial_t polynomialHornerize(polynomial_t p) {
  polynomial_t res;

  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* If p is already hornerized, return a plain copy. */
  if (p->outputType == HORNER_FORM) {
    return polynomialFromCopy(p);
  }

  /* Here we have work to do. We perform a "deep" copy and set the
     output form to "hornerized". 
  */
  res = __polynomialAllocate();
  res->refCount = 1u;
  res->hash.hasHash = 0;
  res->usesExpressionConstant.cached = 0;
  res->type = p->type;
  res->outputType = HORNER_FORM;
  switch (p->type) {
  case SPARSE:
    res->value.sparse = sparsePolynomialFromCopy(p->value.sparse);
    break;
  case ADDITION:
  case SUBTRACTION:
  case MULTIPLICATION:
  case COMPOSITION:
    res->value.pair.g = polynomialFromCopy(p->value.pair.g);
    res->value.pair.h = polynomialFromCopy(p->value.pair.h);
    break;
  case NEGATE:
    res->value.g = polynomialFromCopy(p->value.g);
    break;
  case POWER:
    res->value.powering.g = polynomialFromCopy(p->value.powering.g);
    res->value.powering.c = constantFromCopy(p->value.powering.c);
    break;
  }
  return res;
}

polynomial_t polynomialCanonicalize(polynomial_t p) {
  polynomial_t res;

  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* If p is already canonical, return a plain copy. */
  if (p->outputType == CANONICAL_FORM) {
    return polynomialFromCopy(p);
  }

  /* Here we have work to do. We perform a "deep" copy and set the
     output form to "canonical". 
  */
  res = __polynomialAllocate();
  res->refCount = 1u;
  res->hash.hasHash = 0;
  res->usesExpressionConstant.cached = 0;
  res->type = p->type;
  res->outputType = CANONICAL_FORM;
  switch (p->type) {
  case SPARSE:
    res->value.sparse = sparsePolynomialFromCopy(p->value.sparse);
    break;
  case ADDITION:
  case SUBTRACTION:
  case MULTIPLICATION:
  case COMPOSITION:
    res->value.pair.g = polynomialFromCopy(p->value.pair.g);
    res->value.pair.h = polynomialFromCopy(p->value.pair.h);
    break;
  case NEGATE:
    res->value.g = polynomialFromCopy(p->value.g);
    break;
  case POWER:
    res->value.powering.g = polynomialFromCopy(p->value.powering.g);
    res->value.powering.c = constantFromCopy(p->value.powering.c);
    break;
  }
  return res;
}

void polynomialGetDegree(mpz_t deg, polynomial_t p) {
  mpz_t gd, hd;

  /* Handle stupid input */
  if (p == NULL) {
    mpz_set_si(deg,-1);
    return;
  }

  /* Try to handle the case without sparsifying the polynomial. */
  switch (p->type) {
  case SPARSE:
    sparsePolynomialGetDegree(deg, p->value.sparse);
    return;
    break;
  case ADDITION:
  case SUBTRACTION:
    mpz_init(gd);
    mpz_init(hd);
    polynomialGetDegree(gd, p->value.pair.g);
    polynomialGetDegree(hd, p->value.pair.h);
    if ((mpz_sgn(gd) < 0) || 
	(mpz_sgn(hd) < 0)) {
      mpz_clear(gd);
      mpz_clear(hd);
      mpz_set_si(deg,-1);
      return;
    }
    if (mpz_cmp(gd, hd) != 0) {
      if (mpz_cmp(gd, hd) > 0) {
	mpz_set(deg, gd);
      } else {
	mpz_set(deg, hd);      
      }
      mpz_clear(gd);
      mpz_clear(hd);
      return;
    }
    mpz_clear(gd);
    mpz_clear(hd);
    break;
  case MULTIPLICATION:
    mpz_init(gd);
    mpz_init(hd);
    polynomialGetDegree(gd, p->value.pair.g);
    polynomialGetDegree(hd, p->value.pair.h);
    if ((mpz_sgn(gd) < 0) || 
	(mpz_sgn(hd) < 0)) {
      mpz_clear(gd);
      mpz_clear(hd);
      mpz_set_si(deg,-1);
      return;
    }
    if ((mpz_sgn(gd) > 0) &&
	(mpz_sgn(hd) > 0)) {
      mpz_add(deg, gd, hd);
      mpz_clear(gd);
      mpz_clear(hd);
      return;
    }
    mpz_clear(gd);
    mpz_clear(hd);
    break;
  case COMPOSITION:
    mpz_init(gd);
    mpz_init(hd);
    polynomialGetDegree(gd, p->value.pair.g);
    polynomialGetDegree(hd, p->value.pair.h);
    if ((mpz_sgn(gd) < 0) || 
	(mpz_sgn(hd) < 0)) {
      mpz_clear(gd);
      mpz_clear(hd);
      mpz_set_si(deg,-1);
      return;
    }
    mpz_mul(deg, hd, gd);
    mpz_clear(gd);
    mpz_clear(hd);    
    return;
    break;
  case NEGATE:
    polynomialGetDegree(deg,p->value.g);
    return;
    break;
  case POWER:
    mpz_init(gd);
    polynomialGetDegree(gd, p->value.powering.g);
    if (mpz_sgn(gd) < 0) {
      mpz_clear(gd);
      mpz_set_si(deg,-1);
      return;
    }
    if (mpz_sgn(gd) == 0) {
      mpz_clear(gd);
      mpz_set_si(deg,0);
      return;
    }
    mpz_init(hd);
    if (!tryConstantToMpz(hd, p->value.powering.c)) {
      mpz_clear(gd);
      mpz_clear(hd);
      mpz_set_si(deg,-1);
      return;
    }
    if (mpz_sgn(hd) < 0) {
      mpz_clear(hd);
      mpz_clear(gd);
      mpz_set_si(deg,-1);
      return;
    }    
    mpz_mul(deg, gd, hd);
    mpz_clear(gd);
    mpz_clear(hd);    
    return;
    break;
  }

  /* Fall-back case */
  __polynomialSparsify(p);
  sparsePolynomialGetDegree(deg, p->value.sparse);
}

int polynomialGetDegreeAsInt(polynomial_t p) {
  int gd, hd, res;
  mpz_t degz;
  signed long t, tt;
  unsigned int u, gdu, tu, ttu;

  /* Handle stupid input */
  if (p == NULL) return -1;

  /* Try to handle the case without sparsifying the polynomial. */
  switch (p->type) {
  case SPARSE:
    return sparsePolynomialGetDegreeAsInt(p->value.sparse);
    break;
  case ADDITION:
  case SUBTRACTION:
    gd = polynomialGetDegreeAsInt(p->value.pair.g);
    hd = polynomialGetDegreeAsInt(p->value.pair.h);
    if ((gd == 0) && (hd == 0)) return 0;
    if ((gd >= 0) && (hd >= 0) && (gd != hd)) 
      return (hd > gd ? hd : gd);
    break;
  case MULTIPLICATION:
    gd = polynomialGetDegreeAsInt(p->value.pair.g);
    hd = polynomialGetDegreeAsInt(p->value.pair.h);
    if ((gd != 0) && (hd != 0)) {
      if ((gd < 0) || (hd < 0)) return -1;
      if (tryExactIntAddition(&res, gd, hd)) {
	return res;
      } else {
	return -1;
      }
    }
    break;
  case COMPOSITION:
    gd = polynomialGetDegreeAsInt(p->value.pair.g);
    hd = polynomialGetDegreeAsInt(p->value.pair.h);
    if ((gd == 0) || (hd == 0)) return 0;
    if ((gd < 0) || (hd < 0)) return -1;
    if (tryExactIntMultiplication(&res, gd, hd)) {
      return res;
    } else {
      return -1;
    }    
    break;
  case NEGATE:
    return polynomialGetDegreeAsInt(p->value.g);
    break;
  case POWER:
    gd = polynomialGetDegreeAsInt(p->value.powering.g);
    if (gd == 0) return 0;
    if (tryConstantToUnsignedInt(&u, p->value.powering.c)) {
      if (u == 0u) {
	return 0;
      } else {
	if (gd < 0) return -1;
	gdu = gd;
	if (tryExactUnsignedIntMultiplication(&tu,gdu,u)) {
	  res = tu;
	  if (res < 0) return -1;
	  ttu = res;
	  if (ttu != tu) return -1;
	  return res;
	} else {
	  return -1;
	}
      }
    } else {
      return -1;
    }
    break;
  }

  /* Fall-back case */
  mpz_init(degz);
  polynomialGetDegree(degz, p);
  if (!mpz_fits_slong_p(degz)) {
    res = -1;
  } else {
    t = mpz_get_si(degz);    
    res = t;
    tt = res;
    if (t != tt) {
      res = -1;
    }
  }
  mpz_clear(degz);
  return res;
}

/* Try to compute the i-th coefficient c of p^k, fail if this
   computation is not easy to perform.
*/
static inline int __polynomialGetIthCoefficientAsConstantIntIndexPowerCheap(constant_t *c, polynomial_t p, constant_t k, int i) {
  constant_t ic, mk, t, j, a, b, m, n, r, bin, aPowj, bPowr, prod;

  /* Make compiler happy */
  bin = NULL;
  
  /* Handle stupid inputs */
  if (p == NULL) return 0;
  if (k == NULL) return 0;
  if (i < 0) return 0;

  /* Sparsify the polynomial p */
  __polynomialSparsify(p);

  /* Try to decompose p into p = x^m * (a + b * x^n) */
  if (!sparsePolynomialDecomposeTwoMonomials(&a, &m, &b, &n, p->value.sparse))
    return 0;

  /* Here we have to compute the coefficient of degree i of 

     x^(m * k) * (a + b * x^n)^k

     This means we have to get the coefficient of degree

     t = i - m * k

     of the polynomial 

     (a + b * x^n)^k.

     If t = i - m * k is less than 0, we can simply return 0.

  */
  ic = constantFromInt(i);
  mk = constantMul(m, k);
  t = constantSub(ic, mk);

  /* Check if t is less than 0 */
  if (!(constantIsPositive(t, 1) || constantIsZero(t, 1))) {
    /* Here t = i - m * k is less than 0. 
       
       Assign 0 to c and return success.
       
    */
    *c = constantFromInt(0);
    constantFree(ic);
    constantFree(mk);
    constantFree(t);
    constantFree(a);
    constantFree(m);
    constantFree(b);
    constantFree(n);
    return 1;
  }

  /* Here t = i - m * k is greater than or equal to 0 

     We have to compute the coefficient of degree t of the polynomial

     (a + b * x^n)^k.

     Since 

     (a + b * x^n)^k = sum_j=0^k bin(k,j) * a^j * b^(k - j) * x^(n * j)

     the coefficient of degree t is 

     * equal to bin(k,j) * a^j * b^(k - j) if t is divisible by n and j = t / n
     * equal to 0                          if t is not divisible by n.

     We start by computing j = t / n and then test if j is a
     non-negative integer (hence if t is divisible by n).

  */
  j = constantDiv(t, n);

  /* Check if j is integer */
  if (constantIsNonNegativeInteger(j, 0)) {
    /* The value j = t / n is a non-negative integer. Hence t is divisible by n.

       We have to return c = bin(k,j) * a^j * b^(k - j).

       We first start by computing r = k - j. 

       We then check if r is positive or zero. 
       In this case we compute and return

       c = bin(k,j) * a^j * b^(k - j).

       Otherwise we return c = 0.

    */
    r = constantSub(k, j);
    if (constantIsZero(r, 0) || constantIsPositive(r,0)) {
      /* Here r >= 0 and we have to return 

	 c = bin(k,j) * a^j * b^r.

	 If we can compute the binomial coefficient bin(k,j)
	 we use it, otherwise we simply signal failure.
      */
      if (constantFromBinomial(&bin, k, j)) {
	/* We have 

	   bin = bin(k, j).
	   
	   Compute 

	   c = bin * a^j * b^r.

	*/
	aPowj = constantPow(a, j);
	bPowr = constantPow(b, r);
	prod = constantMul(aPowj, bPowr);
	*c = constantMul(bin, prod);
	constantFree(bin);
	constantFree(aPowj);
	constantFree(bPowr);
	constantFree(prod);
      } else {
	/* Signal failure. */
	constantFree(ic);
	constantFree(mk);
	constantFree(t);
	constantFree(a);
	constantFree(m);
	constantFree(b);
	constantFree(n);
	constantFree(j);
	constantFree(r);
	return 0;
      }
    } else {
      /* Here, r < 0 and we have to return c = 0. */
      *c = constantFromInt(0);
    }
    constantFree(r);
  } else {
    /* The value j = t / n is no non-negative integer. Hence t is not
       divisible by n and the coefficient of degree t of (a + b * x^n)^k is zero. 
    */
    *c = constantFromInt(0);
  }

  /* Free the intermediate values */
  constantFree(ic);
  constantFree(mk);
  constantFree(t);
  constantFree(a);
  constantFree(m);
  constantFree(b);
  constantFree(n);
  constantFree(j);

  /* Return success */
  return 1;
}

static inline constant_t __polynomialGetIthCoefficientAsConstantIntIndex(polynomial_t p, int i) {
  int deg;
  constant_t a, b, res;

  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* Handle case when i < 0 */
  if (i < 0) {
    return constantFromInt(0);
  }

  /* Handle case when i > degree */
  deg = polynomialGetDegreeAsInt(p);
  if ((deg >= 0) && 
      (i > deg)) {
    return constantFromInt(0);
  }

  /* General case */
  switch (p->type) {
  case SPARSE:
    return sparsePolynomialGetIthCoefficientAsConstantIntIndex(p->value.sparse, i);
    break;
  case ADDITION:
  case SUBTRACTION:
    a = __polynomialGetIthCoefficientAsConstantIntIndex(p->value.pair.g,i);
    b = __polynomialGetIthCoefficientAsConstantIntIndex(p->value.pair.h,i);
    if (p->type == ADDITION) {
      res = constantAdd(a, b);
    } else {
      res = constantSub(a, b);
    }
    constantFree(a);
    constantFree(b);
    return res;
    break;
  case NEGATE:
    a = __polynomialGetIthCoefficientAsConstantIntIndex(p->value.g,i);
    res = constantNeg(a);
    constantFree(a);
    return res;
    break;
  case POWER:
    if (__polynomialGetIthCoefficientAsConstantIntIndexPowerCheap(&res,
								  p->value.powering.g,
								  p->value.powering.c,
								  i)) {
      return res;
    }
    break;
  default:
    break;
  }
    
  /* Fall-back case */
  __polynomialSparsify(p);
  return sparsePolynomialGetIthCoefficientAsConstantIntIndex(p->value.sparse, i);
}

static inline constant_t __polynomialGetIthCoefficientAsConstant(polynomial_t p, mpz_t i) {
  mpz_t deg;
  signed long int t, ttt;
  int tt;

  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* Handle case when i < 0 */
  if (mpz_sgn(i) < 0) {
    return constantFromInt(0);
  }

  /* Handle case when i > degree */
  mpz_init(deg);
  polynomialGetDegree(deg, p);
  if ((mpz_sgn(deg) >= 0) && 
      (mpz_cmp(i,deg) > 0)) {
    mpz_clear(deg);
    return constantFromInt(0);
  }

  /* Check if index i holds on a machine int */
  if (mpz_fits_slong_p(i)) {
    t = mpz_get_si(i);
    tt = (int) t;
    ttt = (signed long int) tt;
    if (ttt == t) {
      /* Here i holds on a signed long int and is equal t, which is
	 equal to tt. 
      */
      return __polynomialGetIthCoefficientAsConstantIntIndex(p, tt);
    }
  }

  /* General case */
  __polynomialSparsify(p);
  return sparsePolynomialGetIthCoefficientAsConstant(p->value.sparse, i);
}

node *polynomialGetIthCoefficient(polynomial_t p, mpz_t i) {
  constant_t c;
  node *res;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;

  /* Get coefficient as a constant */
  c = __polynomialGetIthCoefficientAsConstant(p, i);

  /* Convert to an expression */
  res = addMemRef(constantToExpression(c));
  constantFree(c);

  /* Return the result */
  return res;
}

node *polynomialGetIthCoefficientIntIndex(polynomial_t p, int i) {
  constant_t c;
  node *res;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;

  /* Get coefficient as a constant */
  c = __polynomialGetIthCoefficientAsConstantIntIndex(p, i);

  /* Convert to an expression */
  res = addMemRef(constantToExpression(c));
  constantFree(c);

  /* Return the result */
  return res;
}

int polynomialGetCoefficients(node ***coeffs, unsigned int *deg, polynomial_t p) {
  /* Handle stupid input */
  if (p == NULL) return 0;

  /* General case */
  __polynomialSparsify(p);
  return sparsePolynomialGetCoefficients(coeffs, deg, p->value.sparse);
}

static inline node *__polynomialGetExpressionAnyForm(polynomial_t p) {

  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* Handle composition */
  if (p->type == COMPOSITION) {
    __polynomialExecuteComposition(p);
    return __polynomialGetExpressionAnyForm(p);
  }
  
  /* General case */
  switch (p->type) {
  case SPARSE:
    return sparsePolynomialGetExpression(p->value.sparse, 0); /* This must be Horner form */
    break;
  case ADDITION:
    return addMemRef(makeAdd(__polynomialGetExpressionAnyForm(p->value.pair.g),
			     __polynomialGetExpressionAnyForm(p->value.pair.h)));
    break;
  case SUBTRACTION:
    return addMemRef(makeSub(__polynomialGetExpressionAnyForm(p->value.pair.g),
			     __polynomialGetExpressionAnyForm(p->value.pair.h)));    
    break;
  case MULTIPLICATION:
    return addMemRef(makeMul(__polynomialGetExpressionAnyForm(p->value.pair.g),
			     __polynomialGetExpressionAnyForm(p->value.pair.h)));        
    break;
  case NEGATE:
    return addMemRef(makeNeg(__polynomialGetExpressionAnyForm(p->value.g)));
    break;
  case POWER:
    return addMemRef(makePow(__polynomialGetExpressionAnyForm(p->value.powering.g),
			     constantToExpression(p->value.powering.c)));            
    break;
  default:
    break;
  }

  /* Cannot be reached */
  return NULL;
}

node *polynomialGetExpressionExplicit(polynomial_t p) {

  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* If we can output the polynomial in any form, we do so. */
  if ((p->outputType == ANY_FORM) ||
      __polynomialIsConstantCheap(p)) {
    return __polynomialGetExpressionAnyForm(p);
  }

  /* Here, we have to output a "hornerized" or "canonicalized"
     form. 
  */
  __polynomialSparsify(p);
  return sparsePolynomialGetExpression(p->value.sparse, 
				       (p->outputType == CANONICAL_FORM));
}

node *polynomialGetExpression(polynomial_t p) {
  node *res;

  /* Handle stupid input */
  if (p == NULL) return polynomialGetExpressionExplicit(p);

  /* Try to build a lazy memrefed expression */
  res = addMemRefEvenOnNull(NULL);
  if (res != NULL) {
    if (res->nodeType == MEMREF) {
      res->polynomialRepresentation = polynomialFromCopy(p);
      return res;
    } else {
      freeThing(res);
    }
  }

  /* Could not create a lazy memrefed expression */
  return polynomialGetExpressionExplicit(p);
}

void polynomialFPrintf(FILE *fd, polynomial_t p) {
  node *t;

  /* Handle stupid cases */
  if (p == NULL) {
    sollyaFprintf(fd, "(null)");
  }

  /* Handle the general case.

     Can be optimized.

  */
  t = polynomialGetExpressionExplicit(p);
  sollyaFprintf(fd, "%b", t);
  freeThing(t);
}

char *polynomialToString(polynomial_t p) {
  node *t;
  char *str;
  char staticStr[8];
  int size, r;

  /* Handle stupid cases */
  if (p == NULL) return NULL;

  /* Handle the general case.

     Can be optimized.

  */
  t = polynomialGetExpressionExplicit(p);
  size = sollya_snprintf(staticStr,8,"%b",t);
  if (size < 0) {
    freeThing(t);
    return NULL;
  }
  str = (char *) safeCalloc(size + 2, sizeof(char));
  r = sollya_snprintf(str,size,"%b",t);
  if (r < 0) {
    freeThing(t);
    safeFree(str);
    return NULL;
  }

  /* Return the string */
  return str;
}

void polynomialEvalMpfr(mpfr_t y, polynomial_t p, mpfr_t x) {
  mpfr_t scratch, Y;

  /* Handle stupid inputs */
  if (p == NULL) {
    mpfr_set_nan(y);
    return;
  }

  /* If the polynomial is in sparse form, use the evaluation routine
     for sparse forms. 
  */
  if (p->type == SPARSE) {
    sparsePolynomialEvalMpfr(y, p->value.sparse, x);
    return;
  }

  /* If the polynomial is in negation form, call ourselves recursively
     and negate y. 
  */
  if (p->type == NEGATE) {
    polynomialEvalMpfr(y, p->value.g, x);
    mpfr_neg(y, y, GMP_RNDN);
    return;
  }

  /* Check if we can use y as a scratch space for recursion */
  if (x == y) {
    /* We can't. Call ourselves recursively. */
    mpfr_init2(Y, mpfr_get_prec(y));
    polynomialEvalMpfr(Y, p, x);
    mpfr_set(y, Y, GMP_RNDN);
    mpfr_clear(Y);
    return;
  }

  /* Here, we have an addition, subtraction, multiplication,
     composition or power form to handle and we may use y as a scratch
     space.

     We need another scratch variable. 

  */
  mpfr_init2(scratch, mpfr_get_prec(y));
  /* Handle composition apart */
  if (p->type == COMPOSITION) {
    polynomialEvalMpfr(scratch, p->value.pair.h, x);
    polynomialEvalMpfr(y, p->value.pair.g, scratch);
  } else {
    /* Everything else but composition */
    switch (p->type) {
    case ADDITION:
    case SUBTRACTION:
    case MULTIPLICATION:
      polynomialEvalMpfr(y, p->value.pair.g, x);
      polynomialEvalMpfr(scratch, p->value.pair.h, x);
      break;
    case POWER:
      polynomialEvalMpfr(y, p->value.powering.g, x);
      constantEvalMpfr(scratch, p->value.powering.c);
      break;
    default:
      break;
    }
    switch (p->type) {
    case ADDITION:
      mpfr_add(y, y, scratch, GMP_RNDN);
      break;
    case SUBTRACTION:
      mpfr_sub(y, y, scratch, GMP_RNDN);
      break;
    case MULTIPLICATION:
      mpfr_mul(y, y, scratch, GMP_RNDN);
      break;
    case POWER:
      mpfr_pow(y, y, scratch, GMP_RNDN);
      break;
    default:
      break;
    }
  }
  mpfr_clear(scratch);
}

void polynomialEvalMpfi(sollya_mpfi_t y, polynomial_t p, sollya_mpfi_t x) {
  mp_prec_t prec;
  sollya_mpfi_t scr, Y;
  sollya_mpfi_t *scratch;
  sollya_mpfi_t *reusedVars;
  int usingReused;

  /* Handle stupid inputs */
  if (p == NULL) {
    sollya_mpfi_set_nan(y);
    return;
  }

  /* If the polynomial is in sparse form, use the evaluation routine
     for sparse forms. 
  */
  if (p->type == SPARSE) {
    sparsePolynomialEvalMpfi(y, p->value.sparse, x);
    return;
  }

  /* If the polynomial is in negation form, call ourselves recursively
     and negate y. 
  */
  if (p->type == NEGATE) {
    polynomialEvalMpfi(y, p->value.g, x);
    sollya_mpfi_neg(y, y);
    return;
  }

  /* Check if we can use y as a scratch space for recursion */
  if (x == y) {
    /* We can't. Call ourselves recursively. */
    sollya_mpfi_init2(Y, sollya_mpfi_get_prec(y));
    polynomialEvalMpfi(Y, p, x);
    sollya_mpfi_set(y, Y);
    sollya_mpfi_clear(Y);
    return;
  }

  /* Here, we have an addition, subtraction, multiplication,
     composition or power form to handle and we may use y as a scratch
     space.

     We need another scratch variable. We can try to get it by two
     means.

  */
  prec = sollya_mpfi_get_prec(y);
  reusedVars = getReusedGlobalMPFIVars(1, prec);
  if (reusedVars == NULL) {
    usingReused = 0;
    mpfi_init2(scr, prec);
    scratch = &scr;
  } else {
    usingReused = 1;
    scratch = reusedVars;
  }

  /* Handle composition apart */
  if (p->type == COMPOSITION) {
      polynomialEvalMpfi(*scratch, p->value.pair.h, x);
      polynomialEvalMpfi(y, p->value.pair.g, *scratch);
  } else {
    /* Everything else but composition */
    switch (p->type) {
    case ADDITION:
    case SUBTRACTION:
    case MULTIPLICATION:
      polynomialEvalMpfi(y, p->value.pair.g, x);
      polynomialEvalMpfi(*scratch, p->value.pair.h, x);
      break;
    case POWER:
      polynomialEvalMpfi(y, p->value.powering.g, x);
      constantEvalMpfi(*scratch, p->value.powering.c);
      break;
    default:
      break;
    }
    switch (p->type) {
    case ADDITION:
      sollya_mpfi_add(y, y, *scratch);
      break;
    case SUBTRACTION:
      sollya_mpfi_sub(y, y, *scratch);
      break;
    case MULTIPLICATION:
      sollya_mpfi_mul(y, y, *scratch);
      break;
    case POWER:
      sollya_mpfi_pow(y, y, *scratch);
      break;
    default:
      break;
    }
  }

  if (usingReused) {
    returnReusedGlobalMPIVars(1);
  } else {
    sollya_mpfi_clear(scr);
  }
}

static inline int __polynomialCheapCheckConstantInteger(polynomial_t p, int k) {
  int deg;
  sollya_mpfi_t X, Y;
  mpfr_t temp;
  constant_t c, d;
  int res;
  
  /* Handle stupid inputs */
  if (p == NULL) return 0;

  /* General case 

     If the polynomial is sparse, we can use the appropriate function
     for sparse polynomials.

  */
  if (p->type == SPARSE) {
    return sparsePolynomialIsConstantInteger(p->value.sparse, k, 0);
  }

  /* The polynomial is not sparse.

     Try to determine its degree.

     If we can't determine the degree, we just fail.

     If we can and it is zero, we do the following:
     
     (i)  We evaluate the polynomial at [1;1] with interval
          arithmetic. If k is not in the image interval,
	  the polynomial can't be equal to that constant k.
     (ii) We get the coefficient of degree zero and 
          compare it to k.
  */
  deg = __polynomialGetDegreeAsIntCheap(p);

  /* Fail if we couldn't get the degree. Any fall-back would seem to
     be too expensive. 
  */
  if (deg < 0) return 0;

  /* If the degree is not zero, the polynomial can't be constant */
  if (deg != 0) return 0;

  /* Here the degree is zero. 

     Start by evaluating the polynomial at [1;1].

  */
  sollya_mpfi_init2(X, 12);
  sollya_mpfi_init2(Y, 8 * sizeof(int) + 10);
  sollya_mpfi_set_si(X, 1);
  polynomialEvalMpfi(Y, p, X);
  if (!sollya_mpfi_has_nan(Y)) {
    mpfr_init2(temp, 8 * sizeof(int) + 10);
    mpfr_set_si(temp, k, GMP_RNDN); /* exact, enough precision */
    if (!sollya_mpfi_fr_in_interval(temp, Y)) {
      /* k is not in the non-NaN image interval */
      mpfr_clear(temp);
      sollya_mpfi_clear(Y);
      sollya_mpfi_clear(X);
      return 0;
    }
    mpfr_clear(temp);
  }
  sollya_mpfi_clear(Y);
  sollya_mpfi_clear(X);

  /* Get the coefficient of degree zero and compare it to k */
  c = __polynomialGetIthCoefficientAsConstantIntIndex(p, 0);
  d = constantFromInt(k);
  res = constantIsEqual(c, d, 0);
  constantFree(d);
  constantFree(c);
  
  /* Return result */
  return res;
}

static inline int __polynomialCheapCheckConstantZero(polynomial_t p) {
  
  /* Handle stupid inputs */
  if (p == NULL) return 0;
  
  /* General case */
  switch (p->type) {
  case SPARSE:
    return sparsePolynomialIsConstantZero(p->value.sparse, 0);
    break;
  case ADDITION:
  case SUBTRACTION:
    return __polynomialCheapCheckConstantInteger(p, 0);
    break;
  case MULTIPLICATION:
    return (__polynomialCheapCheckConstantZero(p->value.pair.g) ||
	    __polynomialCheapCheckConstantZero(p->value.pair.h));
    break;
  case COMPOSITION:
    return __polynomialCheapCheckConstantInteger(p, 0);
    break;
  case NEGATE:
    return __polynomialCheapCheckConstantZero(p->value.g);
    break;
  case POWER:
    return __polynomialCheapCheckConstantZero(p->value.powering.g);
    break;
  }
  return 0;
}

static inline int __polynomialCheapCheckConstantOne(polynomial_t p) {
  
  /* Handle stupid inputs */
  if (p == NULL) return 0;

  /* General case */
  switch (p->type) {
  case SPARSE:
    return sparsePolynomialIsConstantOne(p->value.sparse, 0);
    break;
  case ADDITION:
  case SUBTRACTION:
    return __polynomialCheapCheckConstantInteger(p, 1);
    break;
  case MULTIPLICATION:
    if (__polynomialCheapCheckConstantOne(p->value.pair.g) &&
	__polynomialCheapCheckConstantOne(p->value.pair.h))
      return 1;
    return __polynomialCheapCheckConstantInteger(p, 1);
    break;
  case COMPOSITION:
    return __polynomialCheapCheckConstantInteger(p, 1);
    break;
  case NEGATE:
    return __polynomialCheapCheckConstantInteger(p, 1);
    break;
  case POWER:
    return __polynomialCheapCheckConstantOne(p->value.powering.g);
    break;
  }
  return 0;
}

polynomial_t polynomialAdd(polynomial_t p, polynomial_t q) {
  polynomial_t res;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;
  if (q == NULL) return NULL;
  
  /* If an easy check shows that p or q are constant zero, return a
     copy of the other 
  */
  if (__polynomialCheapCheckConstantZero(p)) 
    return polynomialFromCopy(q);
  if (__polynomialCheapCheckConstantZero(q)) 
    return polynomialFromCopy(p);

  /* If both polynomials are sparse, perform the operation on sparse
     polynomials. 
  */
  if ((p->type == SPARSE) && (q->type == SPARSE)) {
    return __polynomialBuildFromSparse(sparsePolynomialAdd(p->value.sparse,
							   q->value.sparse));
  }

  /* General case: construct the addition polynomial */
  res = __polynomialAllocate();
  res->refCount = 1u;
  res->hash.hasHash = 0;
  res->usesExpressionConstant.cached = 0;
  res->type = ADDITION;
  res->outputType = ANY_FORM;
  res->value.pair.g = polynomialFromCopy(p);
  res->value.pair.h = polynomialFromCopy(q);
  
  /* Return the result polynomial */
  return res;
}

polynomial_t polynomialSub(polynomial_t p, polynomial_t q) {
  polynomial_t res;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;
  if (q == NULL) return NULL;
  
  /* If an easy check shows that p or q are constant zero, return a
     copy of the other or its negation
  */
  if (__polynomialCheapCheckConstantZero(p)) 
    return polynomialNeg(q);
  if (__polynomialCheapCheckConstantZero(q)) 
    return polynomialFromCopy(p);

  /* If both polynomials are sparse, perform the operation on sparse
     polynomials. 
  */
  if ((p->type == SPARSE) && (q->type == SPARSE)) {
    return __polynomialBuildFromSparse(sparsePolynomialSub(p->value.sparse,
							   q->value.sparse));
  }

  /* General case: construct the addition polynomial */
  res = __polynomialAllocate();
  res->refCount = 1u;
  res->hash.hasHash = 0;
  res->usesExpressionConstant.cached = 0;
  res->type = SUBTRACTION;
  res->outputType = ANY_FORM;
  res->value.pair.g = polynomialFromCopy(p);
  res->value.pair.h = polynomialFromCopy(q);
  
  /* Return the result polynomial */
  return res;
}

polynomial_t polynomialMul(polynomial_t p, polynomial_t q) {
  polynomial_t res;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;
  if (q == NULL) return NULL;
  
  /* If an easy check shows that p or q are constant zero, return a
     copy of them.
  */
  if (__polynomialCheapCheckConstantZero(p)) 
    return polynomialFromCopy(p);
  if (__polynomialCheapCheckConstantZero(q)) 
    return polynomialFromCopy(q);

  /* If an easy check shows that p or q are constant one, return a
     copy of the other one.
  */
  if (__polynomialCheapCheckConstantOne(p)) 
    return polynomialFromCopy(q);
  if (__polynomialCheapCheckConstantOne(q)) 
    return polynomialFromCopy(p);
  
  /* If both polynomials are sparse, and at least one of them has only
     one monomial (or less), perform the operation on sparse
     polynomials.
  */
  if ((p->type == SPARSE) && (q->type == SPARSE)) {
    if ((sparsePolynomialGetMonomialCount(p->value.sparse) <= 1u) ||
	(sparsePolynomialGetMonomialCount(q->value.sparse) <= 1u)) {
      return __polynomialBuildFromSparse(sparsePolynomialMul(p->value.sparse,
							     q->value.sparse));
    }
  }

  /* General case: construct the multiplication polynomial */
  res = __polynomialAllocate();
  res->refCount = 1u;
  res->hash.hasHash = 0;
  res->usesExpressionConstant.cached = 0;
  res->type = MULTIPLICATION;
  res->outputType = ANY_FORM;
  res->value.pair.g = polynomialFromCopy(p);
  res->value.pair.h = polynomialFromCopy(q);
  
  /* Return the result polynomial */
  return res;
}

polynomial_t polynomialNeg(polynomial_t p) {
  polynomial_t res;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;

  /* Handle certain cases in an ad-hoc way */
  switch (p->type) {
  case SPARSE:
    return __polynomialBuildFromSparse(sparsePolynomialNeg(p->value.sparse));
    break;
  case SUBTRACTION:
    res = __polynomialAllocate();
    res->refCount = 1u;
    res->hash.hasHash = 0;
    res->usesExpressionConstant.cached = 0;
    res->type = SUBTRACTION;
    res->outputType = ANY_FORM;
    res->value.pair.g = polynomialFromCopy(p->value.pair.h);
    res->value.pair.h = polynomialFromCopy(p->value.pair.g);
    return res;
    break;
  case NEGATE:
    return polynomialFromCopy(p->value.g);
    break;
  default:
    break;
  }

  /* General case: construct the negation polynomial */
  res = __polynomialAllocate();
  res->refCount = 1u;
  res->hash.hasHash = 0;
  res->usesExpressionConstant.cached = 0;
  res->type = NEGATE;
  res->outputType = ANY_FORM;
  res->value.g = polynomialFromCopy(p);
  
  /* Return the result polynomial */
  return res;
}

polynomial_t polynomialCompose(polynomial_t p, polynomial_t q) {
  polynomial_t res, g, h;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;
  if (q == NULL) return NULL;
    
  /* If q easily proves constant, just construct the composition
     polynomial. This takes almost no space and no time and allows for
     later approximate evaluation of the constructed constant.
  */
  if (__polynomialIsConstantCheap(q)) {
    res = __polynomialAllocate();
    res->refCount = 1u;
    res->hash.hasHash = 0;
    res->usesExpressionConstant.cached = 0;
    res->type = COMPOSITION;
    res->outputType = ANY_FORM;
    res->value.pair.g = polynomialFromCopy(p);
    res->value.pair.h = polynomialFromCopy(q);
    return res;
  }

  /* If both polynomials are sparse, and polynomial q only has one
     monomial, perform the operation on sparse polynomials.
  */
  if ((p->type == SPARSE) && (q->type == SPARSE)) {
    if (sparsePolynomialGetMonomialCount(q->value.sparse) <= 1u) {
      return __polynomialBuildFromSparse(sparsePolynomialCompose(p->value.sparse,
								 q->value.sparse));
    }
  }

  /* If q is sparse and the identity polynomial, just return a copy of
     p. 
  */
  if ((q->type == SPARSE) &&
      (sparsePolynomialIsIdentity(q->value.sparse, 0))) {
    return polynomialFromCopy(p);
  }

  /* If p is of the form p = r^c, we have p(q) = (r(q))^c */
  if (p->type == POWER) {
    res = __polynomialAllocate();
    res->refCount = 1u;
    res->hash.hasHash = 0;
    res->usesExpressionConstant.cached = 0;
    res->type = POWER;
    res->outputType = ANY_FORM;
    res->value.powering.g = polynomialCompose(p->value.powering.g, q);
    res->value.powering.c = constantFromCopy(p->value.powering.c);
    return res;
  }

  /* If q is sparse and only has one monomial but is not constant,
     perform the composition by recursion.
  */
  if ((q->type == SPARSE) &&
      (sparsePolynomialGetMonomialCount(q->value.sparse) <= 1u) &&
      (!sparsePolynomialIsConstant(q->value.sparse,0))) {
    switch (p->type) {
    case SPARSE:
      return __polynomialBuildFromSparse(sparsePolynomialCompose(p->value.sparse,
								 q->value.sparse));
      break;
    case ADDITION:
      g = polynomialCompose(p->value.pair.g, q);
      h = polynomialCompose(p->value.pair.h, q);
      res = polynomialAdd(g, h);
      polynomialFree(g);
      polynomialFree(h);
      break;
    case SUBTRACTION:
      g = polynomialCompose(p->value.pair.g, q);
      h = polynomialCompose(p->value.pair.h, q);
      res = polynomialSub(g, h);
      polynomialFree(g);
      polynomialFree(h);
      break;
    case MULTIPLICATION:
      g = polynomialCompose(p->value.pair.g, q);
      h = polynomialCompose(p->value.pair.h, q);
      res = polynomialMul(g, h);
      polynomialFree(g);
      polynomialFree(h);
      break;
    case COMPOSITION:
      h = polynomialCompose(p->value.pair.h, q);
      res = polynomialCompose(p->value.pair.g, h);
      polynomialFree(h);
      return res;
      break;
    case NEGATE:
      g = polynomialCompose(p->value.g, q);
      res = polynomialNeg(g);
      polynomialFree(g);
      break;
    case POWER:
      res = __polynomialAllocate();
      res->refCount = 1u;
      res->hash.hasHash = 0;
      res->usesExpressionConstant.cached = 0;
      res->type = POWER;
      res->outputType = ANY_FORM;
      res->value.powering.g = polynomialCompose(p->value.powering.g, q);
      res->value.powering.c = constantFromCopy(p->value.powering.c);
      return res;      
      break;
    default:
      break;
    }
  }

  /* General case: construct the composed polynomial */
  res = __polynomialAllocate();
  res->refCount = 1u;
  res->hash.hasHash = 0;
  res->usesExpressionConstant.cached = 0;
  res->type = COMPOSITION;
  res->outputType = ANY_FORM;
  res->value.pair.g = polynomialFromCopy(p);
  res->value.pair.h = polynomialFromCopy(q);
  
  /* Return the result polynomial */
  return res;
}

polynomial_t polynomialDeriv(polynomial_t p) {
  polynomial_t res;
  constant_t one;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;

  /* Easy check if the polynomial actually is constant */
  if (__polynomialIsConstantCheap(p)) {
    return polynomialFromIntConstant(0);
  }

  /* Handle certain cases in an ad-hoc way */
  switch (p->type) {
  case SPARSE:
    return __polynomialBuildFromSparse(sparsePolynomialDeriv(p->value.sparse));
    break;
  case ADDITION:
  case SUBTRACTION:
    res = __polynomialAllocate();
    res->refCount = 1u;
    res->hash.hasHash = 0;
    res->usesExpressionConstant.cached = 0;
    res->type = p->type;
    res->outputType = ANY_FORM;
    res->value.pair.g = polynomialDeriv(p->value.pair.h);
    res->value.pair.h = polynomialDeriv(p->value.pair.g);
    return res;
    break;
  case NEGATE:
    res = __polynomialAllocate();
    res->refCount = 1u;
    res->hash.hasHash = 0;
    res->usesExpressionConstant.cached = 0;
    res->type = NEGATE;
    res->outputType = ANY_FORM;
    res->value.g = polynomialDeriv(p->value.g);
    return res;
    break;
  case POWER:
    if ((!constantIsZero(p->value.powering.c, 1)) &&
	(!constantIsOne(p->value.powering.c, 1))) {
      one = constantFromInt(1);
      res = __polynomialAllocate();
      res->refCount = 1u;
      res->hash.hasHash = 0;
      res->usesExpressionConstant.cached = 0;
      res->type = MULTIPLICATION;
      res->outputType = ANY_FORM;
      res->value.pair.g = __polynomialFromConstant(p->value.powering.c);
      res->value.pair.h = __polynomialAllocate();
      res->value.pair.h->refCount = 1u;
      res->value.pair.h->hash.hasHash = 0;
      res->value.pair.h->usesExpressionConstant.cached = 0;
      res->value.pair.h->type = POWER;
      res->value.pair.h->outputType = ANY_FORM;
      res->value.pair.h->value.powering.g = polynomialFromCopy(p->value.powering.g);
      res->value.pair.h->value.powering.c = constantSub(p->value.powering.c, one);
      constantFree(one);
      return res;
    }
    break;
  default:
    break;
  }
  
  /* General case: sparsify the polynomial and return the derivative
     of the sparse polynomial.
  */
  __polynomialSparsify(p);  
  return __polynomialBuildFromSparse(sparsePolynomialDeriv(p->value.sparse));
}

polynomial_t polynomialGcd(polynomial_t p, polynomial_t q) {
  
  /* Handle stupid inputs */
  if (p == NULL) return NULL;
  if (q == NULL) return NULL;

  /* General case: sparsify the polynomials and return the gcd 
     of the sparse polynomials.
  */
  __polynomialSparsify(p);
  __polynomialSparsify(q);
  return __polynomialBuildFromSparse(sparsePolynomialGcd(p->value.sparse, q->value.sparse));
}

static inline int __polynomialFromExpressionInner(polynomial_t *r, node *p, node *t) {
  polynomial_t a, b, quot, rest, zero;
  sparse_polynomial_t sr;
  int res;
  constant_t c;
  node *pp;

  /* Handle stupid inputs */
  if (p == NULL) return 0;  

  /* Try to decompose expressions built on c, x, +, -, *, /, -, ^ */
  switch (p->nodeType) {
  case MEMREF:
    if (p->polynomialRepresentation != NULL) {
      *r = polynomialFromCopy(p->polynomialRepresentation);
      return 1;
    }
    res = __polynomialFromExpressionInner(r, getMemRefChild(p), p);
    if (res && (p->polynomialRepresentation == NULL) && (!polynomialReferencesExpression(*r, p))) {
      p->polynomialRepresentation = polynomialFromCopy(*r);
    }
    return res;
    break;
  case VARIABLE:
    *r = polynomialFromIdentity();
    return 1;
    break;
  case CONSTANT:
    *r = polynomialFromMpfrConstant(*(p->value));
    return 1;
    break;
  case NEG:
    if (!__polynomialFromExpressionInner(&a, p->child1, NULL)) 
      break;
    *r = polynomialNeg(a);
    polynomialFree(a);
    return 1;
    break;
  case ADD:
  case SUB:
  case MUL:
  case DIV:
  case POW:
    if (!__polynomialFromExpressionInner(&a, p->child1, NULL)) 
      break;
    if (!__polynomialFromExpressionInner(&b, p->child2, NULL)) {
      polynomialFree(a);
      break;
    }
    res = 0;
    switch (p->nodeType) {
    case ADD:
      *r = polynomialAdd(a, b);
      res = 1;
      break;
    case SUB:
      *r = polynomialSub(a, b);
      res = 1;
      break;
    case MUL:
      *r = polynomialMul(a, b);
      res = 1;
      break;
    case DIV:
      polynomialDiv(&quot, &rest, a, b);
      zero = polynomialFromIntConstant(0);
      if (polynomialEqual(rest, zero, 0)) {
	if (polynomialGetDegreeAsInt(b) == 0) {
	  c = __polynomialGetIthCoefficientAsConstantIntIndex(b, 0);
	  if (!constantIsZero(c, 1)) {
	    *r = quot;
	    res = 1;
	  } else {
	    polynomialFree(quot);
	    res = 0;
	  }
	  constantFree(c);
	} else {
	  *r = quot;
	  res = 1;
	}
      } else {
	polynomialFree(quot);
	res = 0;
      }
      polynomialFree(rest);
      polynomialFree(zero);
      break;
    case POW:
      res = polynomialPow(r, a, b);
      break;
    }
    polynomialFree(a);
    polynomialFree(b);
    if (res) return 1;
    break;
  }

  /* Here, the expression contains some other basic function. It can
     be a polynomial only if it is a constant expression.
  */
  if (t != NULL) {
    pp = t;
  } else {
    pp = p;
  }
  if (isConstant(pp)) {
    if (!__sparsePolynomialFromConstantExpression(&sr, pp)) return 0;
    *r = __polynomialBuildFromSparse(sr);
    return 1;
  }

  return 0;
}

int polynomialFromExpression(polynomial_t *r, node *p) {
  return __polynomialFromExpressionInner(r, p, NULL);
}

static inline int __polynomialFromExpressionOnlyRealCoeffsInner(polynomial_t *r, node *p, node *t) {
  polynomial_t a, b, quot, rest, zero;
  sparse_polynomial_t sr;
  int res;
  constant_t c;
  node *pp;

  /* Handle stupid inputs */
  if (p == NULL) return 0;  

  /* Try to decompose expressions built on c, x, +, -, *, /, -, ^ */
  switch (p->nodeType) {
  case MEMREF:
    if (p->polynomialRepresentation != NULL) {
      *r = polynomialFromCopy(p->polynomialRepresentation);
      return 1;
    }
    res = __polynomialFromExpressionOnlyRealCoeffsInner(r, getMemRefChild(p), p);
    if (res && (p->polynomialRepresentation == NULL) && (!polynomialReferencesExpression(*r, p))) {
      p->polynomialRepresentation = polynomialFromCopy(*r);
    }
    return res;
    break;
  case VARIABLE:
    *r = polynomialFromIdentity();
    return 1;
    break;
  case CONSTANT:
    if (mpfr_number_p(*(p->value))) {
      *r = polynomialFromMpfrConstant(*(p->value));
      return 1;
    } else {
      return 0;
    }
    break;
  case NEG:
    if (!__polynomialFromExpressionOnlyRealCoeffsInner(&a, p->child1, NULL)) 
      break;
    *r = polynomialNeg(a);
    polynomialFree(a);
    return 1;
    break;
  case ADD:
  case SUB:
  case MUL:
  case DIV:
  case POW:
    if (!__polynomialFromExpressionOnlyRealCoeffsInner(&a, p->child1, NULL))
      break;
    if (!__polynomialFromExpressionOnlyRealCoeffsInner(&b, p->child2, NULL)) {
      polynomialFree(a);
      break;
    }
    res = 0;
    switch (p->nodeType) {
    case ADD:
      *r = polynomialAdd(a, b);
      res = 1;
      break;
    case SUB:
      *r = polynomialSub(a, b);
      res = 1;
      break;
    case MUL:
      *r = polynomialMul(a, b);
      res = 1;
      break;
    case DIV:
      polynomialDiv(&quot, &rest, a, b);
      zero = polynomialFromIntConstant(0);
      if (polynomialEqual(rest, zero, 0)) {
	if (polynomialGetDegreeAsInt(b) == 0) {
	  c = __polynomialGetIthCoefficientAsConstantIntIndex(b, 0);
	  if (!constantIsZero(c, 1)) {
	    *r = quot;
	    res = 1;
	  } else {
	    polynomialFree(quot);
	    res = 0;
	  }
	  constantFree(c);
	} else {
	  *r = quot;
	  res = 1;
	}
      } else {
	polynomialFree(quot);
	res = 0;
      }
      polynomialFree(rest);
      polynomialFree(zero);
      break;
    case POW:      
      res = polynomialPow(r, a, b);
      break;
    }
    polynomialFree(a);
    polynomialFree(b);
    if (res) return 1;
    break;
  }

  /* Here, the expression contains some other basic function. It can
     be a polynomial only if it is a constant expression.
  */
  if (t != NULL) {
    pp = t;
  } else {
    pp = p;
  }
  if (isConstant(pp)) {
    if (!__sparsePolynomialFromConstantExpressionOnlyRealCoeffs(&sr, pp)) return 0;
    *r = __polynomialBuildFromSparse(sr);
    return 1;
  }

  return 0;
}

int polynomialFromExpressionOnlyRealCoeffs(polynomial_t *r, node *p) {
  return __polynomialFromExpressionOnlyRealCoeffsInner(r, p, NULL);
}

polynomial_t polynomialPowUnsignedInt(polynomial_t p, unsigned int n) {
  polynomial_t res;

  /* Handle stupid case */
  if (p == NULL) return NULL;
  
  /* If n is zero, return the polynomial 1. */
  if (n == 0u) return polynomialFromIntConstant(1);
  
  /* If n is one, return a copy of the polynomial p. */
  if (n == 1u) return polynomialFromCopy(p);

  /* If n is two, return p * p. */
  if (n == 2u) return polynomialMul(p, p);

  /* If p is sparse and only has one monomial (or less), immediately
     perform the operation 
  */
  if (p->type == SPARSE) {
    if (sparsePolynomialGetMonomialCount(p->value.sparse) <= 1u) {
      return __polynomialBuildFromSparse(sparsePolynomialPowUnsignedInt(p->value.sparse, 
									n));
    }
  }

  /* General case: construct the powering polynomial */  
  res = __polynomialAllocate();
  res->refCount = 1u;
  res->hash.hasHash = 0;
  res->usesExpressionConstant.cached = 0;
  res->type = POWER;
  res->outputType = ANY_FORM;
  res->value.powering.g = polynomialFromCopy(p);
  res->value.powering.c = constantFromUnsignedInt(n);
  
  /* Return the result polynomial */
  return res;
}

int polynomialPow(polynomial_t *r, polynomial_t p, polynomial_t q) {
  constant_t n, two14;
  polynomial_t res;
  sparse_polynomial_t sp;

  /* Handle stupid inputs */
  if (p == NULL) return 0;
  if (q == NULL) return 0;
  
  /* Get the degree of q 

     If the integer-output degree accessor function 
     indicates failure, this means that q has a degree
     so huge that it does not hold on an integer. In this case,
     p^q is not a polynomial.

     If the integer-output degress accessor function indicates success
     but a degree of at least one, p^q is not a polynomial.

     So the only way that q is a constant is that the integer-output
     degree accessor function returns 0.

  */
  if (polynomialGetDegreeAsInt(q) != 0) return 0;

  /* Here we know that the degree of q is 0, i.e. that it is a
     constant.

     We get its coefficient of degree 0 to get its constant value.

  */
  n = __polynomialGetIthCoefficientAsConstantIntIndex(q, 0);
  if (!constantIsNonNegativeInteger(n,0)) {
    constantFree(n);
    return 0;
  }

  /* Handle the cases n = 0, n = 1 */
  if (constantIsZero(n, 0)) {
    constantFree(n);
    *r = polynomialFromIntConstant(1);
    return 1;
  }
  if (constantIsOne(n, 0)) {
    constantFree(n);
    *r = polynomialFromCopy(p);
    return 1;
  }

  /* If p is sparse and only has one monomial (or less) and n is
     reasonably small (say n <= 2^14), immediately perform the
     operation.
  */
  if (p->type == SPARSE) {
    if (sparsePolynomialGetMonomialCount(p->value.sparse) <= 1u) {
      two14 = constantFromUnsignedInt(1u << 14);
      if (!constantIsGreater(n, two14, 1)) {
	if (sparsePolynomialPowConstant(&sp, p->value.sparse, n)) {
	  constantFree(n);
	  constantFree(two14);
	  *r = __polynomialBuildFromSparse(sp);
	  return 1;
	}
      }
      constantFree(two14);
    }
  }

  /* General case: construct the powering polynomial */
  res = __polynomialAllocate();
  res->refCount = 1u;
  res->hash.hasHash = 0;
  res->usesExpressionConstant.cached = 0;
  res->type = POWER;
  res->outputType = ANY_FORM;
  res->value.powering.g = polynomialFromCopy(p);
  res->value.powering.c = n;
  
  /* Set the output and return success */
  *r = res;
  return 1;
}

int polynomialCoefficientsAreDyadic(polynomial_t p, int defVal) {

  /* Handle stupid input */
  if (p == NULL) return defVal;

  /* Try to handle the case without sparsifying the polynomial */
  switch (p->type) {
  case SPARSE:
    return sparsePolynomialCoefficientsAreDyadic(p->value.sparse, defVal);
    break;
  case ADDITION:
  case SUBTRACTION:
  case MULTIPLICATION:
  case COMPOSITION:
    if (polynomialCoefficientsAreDyadic(p->value.pair.g, 0) &&
	polynomialCoefficientsAreDyadic(p->value.pair.h, 0)) return 1;
    break;
  case NEGATE:
    return polynomialCoefficientsAreDyadic(p->value.g, defVal);
    break;
  case POWER:
    if (constantIsZero(p->value.powering.c, 0)) return 1;
    if (polynomialCoefficientsAreDyadic(p->value.powering.g, 0)) return 1;
    break;
  }

  /* We couldn't decide the case. Sparsify the polynomial. */
  __polynomialSparsify(p);  
  return sparsePolynomialCoefficientsAreDyadic(p->value.sparse, defVal);
}

int polynomialCoefficientsAreRational(polynomial_t p, int defVal) {

  /* Handle stupid input */
  if (p == NULL) return defVal;

  /* Try to handle the case without sparsifying the polynomial */
  switch (p->type) {
  case SPARSE:
    return sparsePolynomialCoefficientsAreRational(p->value.sparse, defVal);
    break;
  case ADDITION:
  case SUBTRACTION:
  case MULTIPLICATION:
  case COMPOSITION:
    if (polynomialCoefficientsAreRational(p->value.pair.g, 0) &&
	polynomialCoefficientsAreRational(p->value.pair.h, 0)) return 1;
    break;
  case NEGATE:
    return polynomialCoefficientsAreRational(p->value.g, defVal);
    break;
  case POWER:
    if (constantIsZero(p->value.powering.c, 0)) return 1;
    if (polynomialCoefficientsAreRational(p->value.powering.g, 0)) return 1;
    break;
  }

  /* We couldn't decide the case. Sparsify the polynomial. */
  __polynomialSparsify(p);  
  return sparsePolynomialCoefficientsAreRational(p->value.sparse, defVal);
}

int polynomialCoefficientsHoldOnPrecBits(polynomial_t p, mp_prec_t prec, int defVal) {

  /* Handle stupid input */
  if (p == NULL) return defVal;
  if (prec < 2) return defVal;
  
  /* Try to handle the case without sparsifying the polynomial */
  switch (p->type) {
  case SPARSE:
    return sparsePolynomialCoefficientsHoldOnPrecBits(p->value.sparse, prec, defVal);
    break;
  case ADDITION:
  case SUBTRACTION:
  case MULTIPLICATION:
  case COMPOSITION:
    if (polynomialCoefficientsHoldOnPrecBits(p->value.pair.g, prec, 0) &&
	polynomialCoefficientsHoldOnPrecBits(p->value.pair.h, prec, 0)) return 1;
    break;
  case NEGATE:
    return polynomialCoefficientsHoldOnPrecBits(p->value.g, prec, defVal);
    break;
  case POWER:
    if (constantIsZero(p->value.powering.c, 0)) return 1;
    if (polynomialCoefficientsHoldOnPrecBits(p->value.powering.g, prec, 0)) return 1;
    break;
  }

  /* We couldn't decide the case. Sparsify the polynomial. */
  __polynomialSparsify(p);  
  return sparsePolynomialCoefficientsHoldOnPrecBits(p->value.sparse, prec, defVal);
}

static inline polynomial_t __polynomialRoundDyadicAnyForm(polynomial_t p, mp_prec_t prec) {
  polynomial_t res;

  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* If the polynomial only has dyadic coefficients, just return a
     copy of it. 
  */
  if (polynomialCoefficientsAreDyadic(p, 0)) 
    return polynomialFromCopy(p);

  /* Handle composition */
  if (p->type == COMPOSITION) {
    __polynomialExecuteComposition(p);
    return __polynomialRoundDyadicAnyForm(p, prec);
  }

  /* Handle the problem recursively, keeping the expression structure
     of the polynomial 
  */
  res = __polynomialAllocate();
  res->refCount = 1u;
  res->hash.hasHash = 0;
  res->usesExpressionConstant.cached = 0;
  res->type = p->type;
  res->outputType = ANY_FORM;
  switch (p->type) {
  case SPARSE:
    res->value.sparse = sparsePolynomialRoundDyadic(p->value.sparse, prec);
    break;
  case ADDITION:
  case SUBTRACTION:
  case MULTIPLICATION:
    res->value.pair.g = __polynomialRoundDyadicAnyForm(p->value.pair.g, prec);
    res->value.pair.h = __polynomialRoundDyadicAnyForm(p->value.pair.h, prec);
    break;
  case NEGATE:
    res->value.g = __polynomialRoundDyadicAnyForm(p->value.g, prec);
    break;
  case POWER:
    res->value.powering.g = __polynomialRoundDyadicAnyForm(p->value.powering.g, prec);
    res->value.powering.c = constantFromCopy(p->value.powering.c);
    break;
  default:
    break;
  }
  return res;
}

polynomial_t polynomialRoundDyadic(polynomial_t p, mp_prec_t prec) {
  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* If the output type is "any form", use a special algorithm */
  if ((p->outputType == ANY_FORM) ||
      __polynomialIsConstantCheap(p)) {
    return __polynomialRoundDyadicAnyForm(p, prec);
  }

  /* If the polynomial only has dyadic coefficients, just return a
     copy of it. 
  */
  if (polynomialCoefficientsAreDyadic(p, 0)) 
    return polynomialFromCopy(p);

  /* Here, at least one of the coefficients is not dyadic. Round it to
     dyadic, by sparsifying it because the output form is "hornerized" or "canonicalized". 
  */
  __polynomialSparsify(p);
  return __polynomialBuildFromSparse(sparsePolynomialRoundDyadic(p->value.sparse, prec));		       
}

static inline polynomial_t __polynomialRoundRationalAnyForm(polynomial_t p, mp_prec_t prec) {
  polynomial_t res;

  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* If the polynomial only has rational coefficients, just return a
     copy of it. 
  */
  if (polynomialCoefficientsAreRational(p, 0)) 
    return polynomialFromCopy(p);

  /* Handle composition */
  if (p->type == COMPOSITION) {
    __polynomialExecuteComposition(p);
    return __polynomialRoundRationalAnyForm(p, prec);
  }

  /* Handle the problem recursively, keeping the expression structure
     of the polynomial 
  */
  res = __polynomialAllocate();
  res->refCount = 1u;
  res->hash.hasHash = 0;
  res->usesExpressionConstant.cached = 0;
  res->type = p->type;
  res->outputType = ANY_FORM;
  switch (p->type) {
  case SPARSE:
    res->value.sparse = sparsePolynomialRoundRational(p->value.sparse, prec);
    break;
  case ADDITION:
  case SUBTRACTION:
  case MULTIPLICATION:
    res->value.pair.g = __polynomialRoundRationalAnyForm(p->value.pair.g, prec);
    res->value.pair.h = __polynomialRoundRationalAnyForm(p->value.pair.h, prec);
    break;
  case NEGATE:
    res->value.g = __polynomialRoundRationalAnyForm(p->value.g, prec);
    break;
  case POWER:
    res->value.powering.g = __polynomialRoundRationalAnyForm(p->value.powering.g, prec);
    res->value.powering.c = constantFromCopy(p->value.powering.c);
    break;
  default:
    break;
  }
  return res;
}

polynomial_t polynomialRoundRational(polynomial_t p, mp_prec_t prec) {
  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* If the output type is "any form", use a special algorithm */
  if ((p->outputType == ANY_FORM) ||
      __polynomialIsConstantCheap(p)) {
    return __polynomialRoundRationalAnyForm(p, prec);
  }

  /* If the polynomial only has rational coefficients, just return a
     copy of it. 
  */
  if (polynomialCoefficientsAreRational(p, 0)) 
    return polynomialFromCopy(p);

  /* Here, at least one of the coefficients is not rational. Round it to
     rational, by sparsifying it because the output form is "hornerized" or "canonicalized". 
  */
  __polynomialSparsify(p);
  return __polynomialBuildFromSparse(sparsePolynomialRoundRational(p->value.sparse, prec));		       
}

static inline polynomial_t __polynomialRoundAnyForm(polynomial_t p, mp_prec_t prec) {
  polynomial_t res;

  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* If the polynomial only has coefficients that hold on prec bits,
     just return a copy of it.
  */
  if (polynomialCoefficientsHoldOnPrecBits(p, prec, 0)) 
    return polynomialFromCopy(p);

  /* Handle composition */
  if (p->type == COMPOSITION) {
    __polynomialExecuteComposition(p);
    return __polynomialRoundAnyForm(p, prec);
  }

  /* Handle the problem recursively, keeping the expression structure
     of the polynomial 
  */
  res = __polynomialAllocate();
  res->refCount = 1u;
  res->hash.hasHash = 0;
  res->usesExpressionConstant.cached = 0;
  res->type = p->type;
  res->outputType = ANY_FORM;
  switch (p->type) {
  case SPARSE:
    res->value.sparse = sparsePolynomialRound(p->value.sparse, prec);
    break;
  case ADDITION:
  case SUBTRACTION:
  case MULTIPLICATION:
    res->value.pair.g = __polynomialRoundAnyForm(p->value.pair.g, prec);
    res->value.pair.h = __polynomialRoundAnyForm(p->value.pair.h, prec);
    break;
  case NEGATE:
    res->value.g = __polynomialRoundAnyForm(p->value.g, prec);
    break;
  case POWER:
    res->value.powering.g = __polynomialRoundAnyForm(p->value.powering.g, prec);
    res->value.powering.c = constantFromCopy(p->value.powering.c);
    break;
  default:
    break;
  }
  return res;
}

polynomial_t polynomialRound(polynomial_t p, mp_prec_t prec) {
  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* If the output type is "any form", use a special algorithm */
  if ((p->outputType == ANY_FORM) ||
      __polynomialIsConstantCheap(p)) {
    return __polynomialRoundAnyForm(p, prec);
  }

  /* If the polynomial only has coefficients that hold on prec bits,
     just return a copy of it.
  */
  if (polynomialCoefficientsHoldOnPrecBits(p, prec, 0)) 
    return polynomialFromCopy(p);

  /* Here, at least one of the coefficients is not rational. Round it to
     rational, by sparsifying it because the output form is "hornerized" or "canonicalized". 
  */
  __polynomialSparsify(p);
  return __polynomialBuildFromSparse(sparsePolynomialRound(p->value.sparse, prec));		       
}

int polynomialIsHornerized(polynomial_t p) {
  if (p == NULL) return 0;
  return (p->outputType == HORNER_FORM);
}

int polynomialIsCanonicalized(polynomial_t p) {
  if (p == NULL) return 0;
  return (p->outputType == CANONICAL_FORM);
}

static inline sparse_polynomial_t __polynomialGetSparsePolynomial(polynomial_t p) {

  /* Handle stupid input */
  if (p == NULL) return NULL;
  
  /* Sparsify the polynomial and return a copy */
  __polynomialSparsify(p);
  return sparsePolynomialFromCopy(p->value.sparse);
}

static inline int __polynomialUsesExpressionConstantInner(polynomial_t p) {

  /* Handle stupid input */
  if (p == NULL) return 0;
  
  /* Handle the problem recursively */
  switch (p->type) {
  case SPARSE:
    return sparsePolynomialUsesExpressionConstant(p->value.sparse);
    break;
  case ADDITION:
  case SUBTRACTION:
  case MULTIPLICATION:
  case COMPOSITION:
    if (polynomialUsesExpressionConstant(p->value.pair.g)) return 1;
    if (polynomialUsesExpressionConstant(p->value.pair.h)) return 1;
    return 0;
    break;
  case NEGATE:
    return polynomialUsesExpressionConstant(p->value.g);
    break;
  case POWER:
    if (constantUsesExpressionConstant(p->value.powering.c)) return 1;
    if (polynomialUsesExpressionConstant(p->value.powering.g)) return 1;
    return 0;
    break;
  }
  return 0;
}

int polynomialUsesExpressionConstant(polynomial_t p) {
  int res;
  
  /* Handle stupid input */
  if (p == NULL) return 0;

  /* Check if the result is in the cache */
  if (p->usesExpressionConstant.cached) {
    return p->usesExpressionConstant.res;
  }

  /* Call inner function */
  res = __polynomialUsesExpressionConstantInner(p);

  /* Put the result into the cache */
  p->usesExpressionConstant.res = res;
  p->usesExpressionConstant.cached = 1;
  
  /* Return the result */
  return res;
}

int polynomialReferencesExpression(polynomial_t p, node *expr) {

  /* Handle stupid input */
  if (p == NULL) return 0;

  /* If the polynomial does not use any expression, it cannot
     reference the given expression 
  */
  if (!polynomialUsesExpressionConstant(p)) {
    return 0;
  }

  /* Handle the problem recursively */
  switch (p->type) {
  case SPARSE:
    return sparsePolynomialReferencesExpression(p->value.sparse, expr);
    break;
  case ADDITION:
  case SUBTRACTION:
  case MULTIPLICATION:
  case COMPOSITION:
    if (polynomialReferencesExpression(p->value.pair.g, expr)) return 1;
    if (polynomialReferencesExpression(p->value.pair.h, expr)) return 1;
    return 0;
    break;
  case NEGATE:
    return polynomialReferencesExpression(p->value.g, expr);
    break;
  case POWER:
    if (constantReferencesExpression(p->value.powering.c, expr)) return 1;
    if (polynomialReferencesExpression(p->value.powering.g, expr)) return 1;
    return 0;
    break;
  }
  return 0;
}

uint64_t polynomialHash(polynomial_t p) {
  uint64_t hash;
  
  /* Handle stupid input */
  if (p == NULL) return hashPointer(NULL);

  /* Check if the hash has already been computed */
  if (p->hash.hasHash) {
    return p->hash.hash;
  }
  
  /* Sparsify the polynomial and return a copy */
  __polynomialSparsify(p);
  hash = sparsePolynomialHash(p->value.sparse);

  /* Cache the result */
  p->hash.hash = hash;
  p->hash.hasHash = 1;
  
  /* Return the hash */
  return hash;
}
