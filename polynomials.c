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

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include "general.h"
#include "execute.h"
#include "infnorm.h"
#include "expression.h"
#include "polynomials.h"

/* Helper types */


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
  constant_t *coeffs;
  constant_t *monomialDegrees;
};

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

static inline mp_exp_t mpq_remove_powers_of_two(mpq_t op) {
  mp_bitcnt_t dyadNum, dyadDen;
  mp_exp_t expo;

  if (mpq_sgn(op) == 0) return ((mp_exp_t) 0);

  dyadNum = mpz_scan1(mpq_numref(op), 0);
  dyadDen = mpz_scan1(mpq_denref(op), 0);
  mpz_tdiv_q_2exp(mpq_numref(op), mpq_numref(op), dyadNum);
  mpz_tdiv_q_2exp(mpq_denref(op), mpq_denref(op), dyadDen);
  mpq_canonicalize(op);
  expo = dyadNum - dyadDen;

  return expo;
}

static inline void scaledMpqAdd(mp_exp_t *EC, mpq_t c, 
				mp_exp_t EA, mpq_t a, 
				mp_exp_t EB, mpq_t b) {
  if (EB <= EA) {
    *EC = EB;
    mpq_mul_2exp(c, a, EA - EB);
    mpq_add(c, c, b);
  } else {
    *EC = EA;
    mpq_mul_2exp(c, b, EB - EA);
    mpq_add(c, c, a);
  }
  mpq_canonicalize(c);
  *EC += mpq_remove_powers_of_two(c);
}

static inline void scaledMpqAddInt(mp_exp_t *EC, mpq_t c, 
				   mp_exp_t EA, mpq_t a, 
				   int b) {
  mpq_t t;

  /* Can be optimized */
  mpq_init(t);
  mpq_set_si(t, b, 1ui);
  mpq_canonicalize(t);
  scaledMpqAdd(EC, c, EA, a, 0, t);
  mpq_clear(t);
}

static inline void scaledMpqSub(mp_exp_t *EC, mpq_t c, 
				mp_exp_t EA, mpq_t a, 
				mp_exp_t EB, mpq_t b) {
  if (EB <= EA) {
    *EC = EB;
    mpq_mul_2exp(c, a, EA - EB);
    mpq_sub(c, c, b);
  } else {
    *EC = EA;
    mpq_mul_2exp(c, b, EB - EA);
    mpq_sub(c, c, a);
  }
  mpq_canonicalize(c);
  *EC += mpq_remove_powers_of_two(c);
}

static inline void scaledMpqSubInt(mp_exp_t *EC, mpq_t c, 
				   mp_exp_t EA, mpq_t a, 
				   int b) {
  mpq_t t;

  /* Can be optimized */
  mpq_init(t);
  mpq_set_si(t, b, 1ui);
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
  mpq_set_si(t, a, 1ui);
  mpq_canonicalize(t);
  scaledMpqSub(EC, c, 0, t, EB, b);
  mpq_clear(t);
}

static inline void scaledMpqMul(mp_exp_t *EC, mpq_t c, 
				mp_exp_t EA, mpq_t a, 
				mp_exp_t EB, mpq_t b) {
  *EC = EA + EB;
  mpq_mul(c, a, b);
  mpq_canonicalize(c);
  *EC += mpq_remove_powers_of_two(c);
}

static inline void scaledMpqMulInt(mp_exp_t *EC, mpq_t c, 
				   mp_exp_t EA, mpq_t a, 
				   int b) {
  mpq_t t;

  /* Can be optimized */
  mpq_init(t);
  mpq_set_si(t, b, 1ui);
  mpq_canonicalize(t);
  scaledMpqMul(EC, c, EA, a, 0, t);
  mpq_clear(t);
}

static inline int tryScaledMpqDiv(mp_exp_t *EC, mpq_t c, 
				  mp_exp_t EA, mpq_t a, 
				  mp_exp_t EB, mpq_t b) {
  if (mpq_sgn(b) == 0) return 0;
  *EC = EA - EB;
  mpq_div(c, a, b);
  mpq_canonicalize(c);
  *EC += mpq_remove_powers_of_two(c);
  return 1;
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

  /* Now canonicalize the fraction, i.e. reduce common factors */
  mpq_canonicalize(a);

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

  EE += mpq_remove_powers_of_two(aa);

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
		   mpz_sizeinbase(mpq_numref(a),2)));
  Emin = E - 2;
  Emax = E + 2;

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
    if (EA - EB >= Emax) return 1;
    if (EA - EB <= Emin) return 0;
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
    if (EA - EB >= Emax) return 0;
    if (EA - EB <= Emin) return 1;
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
  D = EA - EB - F;

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
  sizeOut = EA + sizeNum + 5 - sizeDen;
  if (sizeOut < 12) sizeOut = 12;
 
  /* This upper bound on the size of the output gives us 
     an idea of the precision we need to compute a first guess 
     on the output.
  */
  precOut = sizeOut;
  if (precOut < 12) precOut = 12;

  /* Compute a first guess on the output */
  mpfr_init2(t, precOut);
  mpq_canonicalize(a);
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
  ER += mpq_remove_powers_of_two(r);

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
    mpz_mul_2exp(den, den, (mp_bitcnt_t) (-Erest));    
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
  mpq_canonicalize(a);
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
    mpz_mul_2exp(r, r, (mp_bitcnt_t) (-E));
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
  if (n == 0) {
    *EC = 0;
    mpq_set_ui(c,1,1u);
    mpq_canonicalize(c);
    *EC += mpq_remove_powers_of_two(c);
    return 1;
  }
  if (n == 1) {
    *EC = EA;
    mpq_set(c, a);
    mpq_canonicalize(c);
    *EC += mpq_remove_powers_of_two(c);
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
  *EC += mpq_remove_powers_of_two(c);

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

  if (n == 0) {
    *EC = 0;
    mpq_set_ui(c,1,1u);
    mpq_canonicalize(c);
    *EC += mpq_remove_powers_of_two(c);
    return 1;
  }
  if (n == 1) {
    *EC = EA;
    mpq_set(c, a);
    mpq_canonicalize(c);
    *EC += mpq_remove_powers_of_two(c);
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
  *EC += mpq_remove_powers_of_two(c);

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
      cAbs = a - b;
    } else {
      cNeg = bNeg;
      cAbs = b - a;
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
      cAbs = a - b;
    } else {
      cNeg = bNeg;
      cAbs = b - a;
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

static inline int exactUint64Mul(uint64_t *ch, 
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
  pc = (Ea > Eb? Ea : Eb) - (Ea - pa < Eb - pb ? Ea - pa : Eb - pb) + 1;
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
  pc = (Ea > Eb? Ea : Eb) - (Ea - pa < Eb - pb ? Ea - pa : Eb - pb) + 1;
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
  pc = (Ea > Eb? Ea : Eb) - (Ea - pa < Eb - pb ? Ea - pa : Eb - pb) + 1;
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
  pc = (Ea > Eb? Ea : Eb) - (Ea - pa < Eb - pb ? Ea - pa : Eb - pb) + 1;
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
  pc = pa + pb;
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
  pc = pa + pb;
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

constant_t constantAdd(constant_t, constant_t);
constant_t constantSub(constant_t, constant_t);
constant_t constantMul(constant_t, constant_t);
constant_t constantDiv(constant_t, constant_t);
constant_t constantPow(constant_t, constant_t);
constant_t constantNeg(constant_t);

static inline constant_t __constantAllocate() {
  return (constant_t) safeMalloc(sizeof(struct __constant_struct_t));
}

static inline void __constantFreeMem(constant_t c) {
  safeFree(c);
}

constant_t constantFromMpfr(mpfr_t c) {
  constant_t res;
  int intval;

  res = __constantAllocate();
  res->refCount = 1;
  res->isZero.cached = 0;
  res->isOne.cached = 0;
  res->isNonNegativeInteger.cached = 0;
  res->isPositive.cached = 0;
  if (mpfr_is_machine_integer(&intval, c)) {
    res->type = INTEGER;
    res->value.integer = intval;
    return res;
  } 
  res->type = MPFR;
  mpfr_init2(res->value.mpfr, mpfr_get_prec(c));
  mpfr_set(res->value.mpfr, c, GMP_RNDN); /* exact, same precision */
  
  return res;
}

constant_t constantFromInt(int c) {
  constant_t res;

  res = __constantAllocate();
  res->refCount = 1;
  res->isZero.cached = 0;
  res->isOne.cached = 0;
  res->isNonNegativeInteger.cached = 0;
  res->isPositive.cached = 0;
  res->type = INTEGER;
  res->value.integer = c;
  
  return res;
}

constant_t constantFromUnsignedInt(unsigned int c) {
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


constant_t constantFromScaledMpq(mp_exp_t e, mpq_t c) {
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
  expo = e + dyadNum - dyadDen;

  if (mpz_cmp_si(den, 1) == 0) {
    /* The denominator is one, so we can actually 
       use a MPFR (or integer) representation.
    */
    p = mpz_sizeinbase(num, 2);
    if (p < 12) p = 12;
    mpfr_init2(mpfrval, p);
    mpfr_set_z_2exp(mpfrval, num, expo, GMP_RNDN); /* exact as enough precision */
    res = constantFromMpfr(mpfrval);
    mpfr_clear(mpfrval);
    mpz_clear(den);
    mpz_clear(num);
  
    return res;
  }

  /* Here, we are sure we must use a scaled MPQ representation */
  res = __constantAllocate();
  res->refCount = 1;
  res->isZero.cached = 0;
  res->isOne.cached = 0;
  res->isNonNegativeInteger.cached = 0;
  res->isPositive.cached = 0;
  res->type = SCALEDMPQ;
  res->value.scaledMpq.expo = expo;
  mpq_init(res->value.scaledMpq.significand);
  mpq_set_num(res->value.scaledMpq.significand, num);
  mpq_set_den(res->value.scaledMpq.significand, den);
  mpq_canonicalize(res->value.scaledMpq.significand);
  mpz_clear(den);
  mpz_clear(num);
  
  return res;
}

constant_t constantFromMpq(mpq_t c) {
  return constantFromScaledMpq((mp_exp_t) 0, c);
}

constant_t constantFromMpz(mpz_t c) {
  mpq_t q;
  constant_t res;
  
  mpq_init(q);
  mpq_set_z(q, c);
  res = constantFromMpq(q);
  mpq_clear(q);
  
  return res;
}

constant_t constantFromCopy(constant_t c) {
  if (c == NULL) return NULL;
  c->refCount++;
  return c;
}

constant_t constantFromExpression(node *c) {
  node *simplified;
  mpq_t rational;
  constant_t res;
  
  if (c == NULL) return NULL;
  if (!isConstant(c)) return NULL;
  if (accessThruMemRef(c)->nodeType == CONSTANT) {
    return constantFromMpfr(*(accessThruMemRef(c)->value));
  }
  mpq_init(rational);
  if (tryEvaluateConstantTermToMpq(rational, c)) {
    res = constantFromMpq(rational);
    mpq_clear(rational);
    return res;
  }
  simplified = simplifyRationalErrorfree(c);
  if (accessThruMemRef(simplified)->nodeType == CONSTANT) {
    res = constantFromMpfr(*(accessThruMemRef(simplified)->value));
    freeThing(simplified);
    return res;
  }
  
  res = __constantAllocate();
  res->refCount = 1;
  res->isZero.cached = 0;
  res->isOne.cached = 0;
  res->isNonNegativeInteger.cached = 0;
  res->isPositive.cached = 0;
  res->type = EXPRESSION;
  res->value.expr = simplified;
  
  return res;
}


void constantFree(constant_t c) {
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

int constantIsZero(constant_t a, int defVal) {
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

int constantIsOne(constant_t a, int defVal) {
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
    mpq_canonicalize(a->value.scaledMpq.significand);
    a->value.scaledMpq.expo += mpq_remove_powers_of_two(a->value.scaledMpq.significand);
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

int constantIsNonNegativeInteger(constant_t a, int defVal) {
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
    a->isNonNegativeInteger.res = !(s == 0);
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

int constantIsEqual(constant_t a, constant_t b, int defVal) {
  constant_t d;
  int res, sa, sb;
  mpq_t t;
  mp_exp_t G;
  
  if (a == NULL) return defVal;
  if (b == NULL) return defVal;
  if (a == b) return 1;
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
    mpq_canonicalize(t);
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
    if (a->value.scaledMpq.expo - b->value.scaledMpq.expo != G) {
      res = 0;
    } else {
      res = (mpq_cmp_si(t, 1, 1ul) == 0);
    }
    mpq_clear(t);
    return res;
    break;
  }
 
  return defVal;
}

int constantIsPositive(constant_t a, int defVal) {
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

int constantIsGreater(constant_t a, constant_t b, int defVal) {
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

int constantIsGreaterOrEqual(constant_t a, constant_t b, int defVal) {
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

node *constantToExpression(constant_t a) {
  mpfr_t num, den;
  mp_prec_t p;
  node *res;

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
      res = addMemRef(makeConstant(num));
      mpfr_clear(num);
      return res;
    }
    /* Here, we must take both the numerator and the denominator into
       account.
    */
    p = mpz_sizeinbase(mpq_numref(a->value.scaledMpq.significand), 2);
    if (p < 12) p = 12;
    mpfr_init2(num, p);
    mpfr_set_z_2exp(num, 
		    mpq_numref(a->value.scaledMpq.significand), 
		    a->value.scaledMpq.expo, 
		    GMP_RNDN); /* exact as enough precision */
    p = mpz_sizeinbase(mpq_denref(a->value.scaledMpq.significand), 2);
    if (p < 12) p = 12;
    mpfr_init2(den, p);
    mpfr_set_z(den, 
	       mpq_denref(a->value.scaledMpq.significand), 
	       GMP_RNDN); /* exact as enough precision */
    res = addMemRef(makeDiv(makeConstant(num),makeConstant(den)));
    mpfr_clear(num);
    mpfr_clear(den);
    return res;
    break;
  }
  return NULL;
}

int tryConstantToScaledMpq(mp_exp_t *E, mpq_t rop, constant_t a) {
  mpq_t t;
  mpz_t mant;
  mp_exp_t expo;

  if (a == NULL) return 0;
  switch (a->type) {
  case INTEGER:
    mpq_set_si(rop, a->value.integer, 1ul);
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
      mpq_set_si(rop, 0, 1ul);
      mpq_canonicalize(rop);
      return 1;
    }
    mpz_init(mant);
    expo = mpfr_get_z_2exp(mant, a->value.mpfr);
    mpq_set_z(rop, mant);
    mpq_canonicalize(rop);
    *E = expo + mpq_remove_powers_of_two(rop);      
    mpz_clear(mant);
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

int tryConstantToMpz(mpz_t r, constant_t a) {
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
      mpz_mul_2exp(den, den, (mp_bitcnt_t) (-a->value.scaledMpq.expo));
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
    mpz_mul_2exp(den, den, (mp_bitcnt_t) (-a->value.scaledMpq.expo));
  }
  mpz_fdiv_q(r, num, den);
  mpz_clear(den);
  mpz_clear(num);
  return 1;
}

int tryConstantToInt(int *r, constant_t a) {
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

int tryConstantToUnsignedInt(unsigned int *r, constant_t a) {
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
    sollyaFprintf(fd,"2^(%lld) * %r",
		  (long long int) c->value.scaledMpq.expo, 
		  c->value.scaledMpq.significand);
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
    size = sollya_snprintf(staticStr,8,"2^(%lld) * %r",
			   (long long int) c->value.scaledMpq.expo, 
			   c->value.scaledMpq.significand);
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
    r = sollya_snprintf(str,size,"2^(%lld) * %r",
			(long long int) c->value.scaledMpq.expo, 
			c->value.scaledMpq.significand);
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

constant_t constantAdd(constant_t a, constant_t b) {
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
      res = constantFromMpfr(cM);
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
    }
  }

  /* Handle the case when only one of the inputs is a machine integer */
  if (a->type == INTEGER) {
    switch (b->type) {
    case MPFR:
      mpfr_init2(cM, 12);
      /* Next call may change the precision of cM */
      mpfr_add_exact_int(cM, b->value.mpfr, a->value.integer);
      res = constantFromMpfr(cM);
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
    }
  }
  if (b->type == INTEGER) {
    switch (a->type) {
    case MPFR:
      mpfr_init2(cM, 12);
      /* Next call may change the precision of cM */
      mpfr_add_exact_int(cM, a->value.mpfr, b->value.integer);
      res = constantFromMpfr(cM);
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

constant_t constantSub(constant_t a, constant_t b) {
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
      res = constantFromMpfr(cM);
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
    }
  }

    /* Handle the case when only one of the inputs is a machine integer */
  if (a->type == INTEGER) {
    switch (b->type) {
    case MPFR:
      mpfr_init2(cM, 12);
      /* Next call may change the precision of cM */
      mpfr_int_sub_exact(cM, a->value.integer, b->value.mpfr);
      res = constantFromMpfr(cM);
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
    }
  }
  if (b->type == INTEGER) {
    switch (a->type) {
    case MPFR:
      mpfr_init2(cM, 12);
      /* Next call may change the precision of cM */
      mpfr_sub_exact_int(cM, a->value.mpfr, b->value.integer);
      res = constantFromMpfr(cM);
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

constant_t constantMul(constant_t a, constant_t b) {
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
      res = constantFromMpfr(cM);
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
    }
  }

  /* Handle the case when only one of the inputs is a machine integer */
  if (a->type == INTEGER) {
    switch (b->type) {
    case MPFR:
      mpfr_init2(cM, 12);
      /* Next call may change the precision of cM */
      mpfr_mul_exact_int(cM, b->value.mpfr, a->value.integer);
      res = constantFromMpfr(cM);
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
    }
  }
  if (b->type == INTEGER) {
    switch (a->type) {
    case MPFR:
      mpfr_init2(cM, 12);
      /* Next call may change the precision of cM */
      mpfr_mul_exact_int(cM, a->value.mpfr, b->value.integer);
      res = constantFromMpfr(cM);
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

constant_t constantDiv(constant_t a, constant_t b) {
  node *aExpr, *bExpr, *cExpr;
  constant_t res;
  mpq_t aS, bS, cS, cQ;
  mp_exp_t EA, EB, EC;
  int cI;
  mpfr_t cM;
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

constant_t constantPow(constant_t a, constant_t b) {
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
	  res = constantFromMpfr(cM);
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
	    res = constantFromMpfr(cM);
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

  /* Could not convert to scaled rational form or the division did not work. 
     Convert to expressions. 
  */
  aExpr = constantToExpression(a);
  bExpr = constantToExpression(b);
  cExpr = addMemRef(makePow(aExpr, bExpr));
  res = constantFromExpression(cExpr);
  freeThing(cExpr);
  return res;  
}

constant_t constantNeg(constant_t a) {
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
    cExpr = addMemRef(makeNeg(a->value.expr));
    res = constantFromExpression(cExpr);
    freeThing(cExpr);
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
    mpq_canonicalize(cQ);
    res = constantFromScaledMpq(a->value.scaledMpq.expo, cQ);
    mpq_clear(cQ);
    return res;
    break;
  }
}

/* Computes q and r such that 2^14 * q + r = a and q = floor(a *
   2^-14) 
*/
void constantCutTwo14(constant_t *q, constant_t *r, constant_t a) {
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

void constantEvalMpfr(mpfr_t rop, constant_t c) {
  mp_prec_t p;

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
    evaluateFaithful(rop, c->value.expr, rop, mpfr_get_prec(rop) + 10);
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

void constantEvalMpfi(sollya_mpfi_t rop, constant_t c) {
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
    evaluateConstantExpressionToSharpInterval(rop, c->value.expr);
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

/* End of part for constants */

/* Start of part for sparse polynomials */

static inline sparse_polynomial_t __sparsePolynomialAllocate() {
  return (sparse_polynomial_t) safeMalloc(sizeof(struct __sparse_polynomial_struct_t));
}

static inline void __sparsePolynomialFreeMem(sparse_polynomial_t p) {
  safeFree(p);
}

static inline void __sparsePolynomialAdjustDegree(sparse_polynomial_t p) {
  unsigned int i, k;

  if (p == NULL) return;
  for (k=0,i=p->monomialCount-1;i>=1;i--,k++) {
    if (!constantIsZero(p->coeffs[i],0)) break;
  }
  if (k == 0) return;
  for (i=p->monomialCount-k;i<p->monomialCount;i++) {
    constantFree(p->coeffs[i]);
    constantFree(p->monomialDegrees[i]);
  }
  p->monomialCount -= k;
  constantFree(p->deg);
  p->deg = constantFromCopy(p->monomialDegrees[p->monomialCount]);
  p->coeffs = (constant_t *) safeRealloc(p->coeffs, 
					 ((size_t) (p->monomialCount)) * sizeof(constant_t));
  p->monomialDegrees = (constant_t *) safeRealloc(p->monomialDegrees, 
						  ((size_t) (p->monomialCount)) * sizeof(constant_t));
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
  if ((0u <= g) && (g <= n - 1u)) {
    for (p=n+1u, i=g, j=0; (j<64) && (0u <= i) && (i <= n - 1u); j++) {
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

sparse_polynomial_t sparsePolynomialFromMpfrConstant(mpfr_t c) {
  sparse_polynomial_t res;
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->monomialCount = 1;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->coeffs[0] = constantFromMpfr(c);
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees[0] = constantFromInt(0);
  res->deg = constantFromCopy(res->monomialDegrees[0]);
  return res;
}

sparse_polynomial_t sparsePolynomialFromMpzConstant(mpz_t c) {
  sparse_polynomial_t res;
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->monomialCount = 1;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->coeffs[0] = constantFromMpz(c);
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees[0] = constantFromInt(0);
  res->deg = constantFromCopy(res->monomialDegrees[0]);
  return res;
}

sparse_polynomial_t sparsePolynomialFromMpqConstant(mpq_t c) {
  sparse_polynomial_t res;
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->monomialCount = 1;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->coeffs[0] = constantFromMpq(c);
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees[0] = constantFromInt(0);
  res->deg = constantFromCopy(res->monomialDegrees[0]);
  return res;
}

sparse_polynomial_t sparsePolynomialFromConstant(constant_t c) {
  sparse_polynomial_t res;
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->monomialCount = 1;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->coeffs[0] = constantFromCopy(c);
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees[0] = constantFromInt(0);
  res->deg = constantFromCopy(res->monomialDegrees[0]);
  return res;
}

sparse_polynomial_t sparsePolynomialFromIntConstant(int c) {
  sparse_polynomial_t res;
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->monomialCount = 1;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->coeffs[0] = constantFromInt(c);
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees[0] = constantFromInt(0);
  res->deg = constantFromCopy(res->monomialDegrees[0]);
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
  *r = res;
  return 1;
}

sparse_polynomial_t sparsePolynomialFromIdentity() {
  sparse_polynomial_t res;
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->monomialCount = 1;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->coeffs[0] = constantFromInt(1);
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees[0] = constantFromCopy(res->coeffs[0]);
  res->deg = constantFromCopy(res->monomialDegrees[0]);
  return res;
}

sparse_polynomial_t sparsePolynomialFromMpfrCoefficients(mpfr_t *coeffs, unsigned int deg) {
  unsigned int i, s;
  sparse_polynomial_t res;
  constant_t c;
  
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  res->deg = constantFromUnsignedInt(deg);
  s = deg + 1u;
  if (s == 0u) s = UINT_MAX;
  res->monomialCount = s;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  for (i=0;i<=deg;i++) {
    c = constantFromMpfr(coeffs[i]);
    if ((i==0) || (!constantIsZero(c,0))) {
      res->coeffs[i] = c;
      res->monomialDegrees[i] = constantFromUnsignedInt(i);
    } else {
      constantFree(c);
    }
  }
  res->deg = constantFromCopy(res->monomialDegrees[res->monomialCount-1u]);
  __sparsePolynomialAdjustDegree(res);
  return res;  
}

int sparsePolynomialFromConstantExpressionCoefficients(sparse_polynomial_t *r, node **coeffs, unsigned int deg) {
  unsigned int i, s;
  sparse_polynomial_t res;
  constant_t c;
  
  if (coeffs == NULL) return 0;
  for (i=0;i<=deg;i++) {
    if (!isConstant(coeffs[i])) return 0;
  }
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
  s = deg + 1u;
  if (s == 0u) s = UINT_MAX;
  res->monomialCount = s;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  for (i=0;i<=deg;i++) {
    c = constantFromExpression(coeffs[i]);
    if ((i==0) || (!constantIsZero(c,0))) {
      res->coeffs[i] = c;
      res->monomialDegrees[i] = constantFromUnsignedInt(i);
    } else {
      constantFree(c);
    }
  }
  res->deg = constantFromCopy(res->monomialDegrees[res->monomialCount-1u]);
  __sparsePolynomialAdjustDegree(res);
  *r = res;
  return 1;
}

sparse_polynomial_t sparsePolynomialFromCopy(sparse_polynomial_t p) {
  if (p == NULL) return NULL;
  p->refCount++;
  return p;
}

void sparsePolynomialFree(sparse_polynomial_t p) {
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

int sparsePolynomialEqual(sparse_polynomial_t p, sparse_polynomial_t q, int defVal) {
  unsigned int i;
  int monDegEqual, coeffEqual, degreeEqual;

  if (p == NULL) return defVal;
  if (q == NULL) return defVal;
  if (p == q) return 1;
  if (p->monomialCount != q->monomialCount) return 0;
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

int sparsePolynomialIsConstant(sparse_polynomial_t p, int defVal) {
  int degZero;

  if (p == NULL) return defVal;
  if (p->monomialCount == 0u) return 1;
  degZero = constantIsZero(p->deg, 42);
  if (degZero == 42) return defVal;
  if (!degZero) return 0;
  return 1;
}

int sparsePolynomialConstantGetConstant(constant_t *c, sparse_polynomial_t p) {
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

int sparsePolynomialIsConstantZero(sparse_polynomial_t p, int defVal) {
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

int sparsePolynomialIsConstantOne(sparse_polynomial_t p, int defVal) {
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

sparse_polynomial_t sparsePolynomialAdd(sparse_polynomial_t p, sparse_polynomial_t q) {
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
	newCoeff = constantFromCopy(p->coeffs[i]);
	newMonomial = constantFromCopy(p->monomialDegrees[i]);
	i++;
      } else {
	/* Case (ii) */
	newCoeff = constantFromCopy(q->coeffs[j]);
	newMonomial = constantFromCopy(q->monomialDegrees[j]);
	j++;
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
  if (res->monomialCount == 0) {
    res->coeffs[res->monomialCount] = constantFromInt(0);
    res->monomialDegrees[res->monomialCount] = constantFromInt(0);
    (res->monomialCount)++;
  }
  /* Set the degree of the polynomial */
  res->deg = constantFromCopy(res->monomialDegrees[res->monomialCount - 1]);
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

sparse_polynomial_t sparsePolynomialNeg(sparse_polynomial_t);

sparse_polynomial_t sparsePolynomialSub(sparse_polynomial_t p, sparse_polynomial_t q) {
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
      newCoeff = constantAdd(p->coeffs[i], q->coeffs[j]);
      newMonomial = constantFromCopy(p->monomialDegrees[i]);
      i++;
      j++;
    } else {
      if (constantIsGreater(p->monomialDegrees[i], q->monomialDegrees[j], 0)) {
	/* Case (iii) */
	newCoeff = constantFromCopy(p->coeffs[i]);
	newMonomial = constantFromCopy(p->monomialDegrees[i]);
	i++;
      } else {
	/* Case (ii) */
	newCoeff = constantNeg(q->coeffs[j]);
	newMonomial = constantFromCopy(q->monomialDegrees[j]);
	j++;
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
  if (res->monomialCount == 0) {
    res->coeffs[res->monomialCount] = constantFromInt(0);
    res->monomialDegrees[res->monomialCount] = constantFromInt(0);
    (res->monomialCount)++;
  }
  /* Set the degree of the polynomial */
  res->deg = constantFromCopy(res->monomialDegrees[res->monomialCount - 1]);
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

     b) degree[j] > d. In this case we deplace all existing
        coefficients and degrees of indices j thru n-1 and put the new
	coefficient and degree at index j.

  */
  if (constantIsEqual(degrees[j],d,0)) {
    /* Case a) */
    t = constantAdd(coeffs[j], c);
    constantFree(coeffs[j]);
    coeffs[j] = t;
    constantFree(c);
    constantFree(d);
    return;
  }

  /* Case b) 

     Deplace the existing coefficients and degrees.
  */
  for (i=*n;i>j;i--) {
    coeffs[i] = coeffs[i-1];
    degrees[i] = degrees[i-1];
  }

  /* Account for the new entry */
  (*n)++;

  /* Put the new entry */
  coeffs[j] = c;
  degrees[j] = d;
}

sparse_polynomial_t sparsePolynomialAddConstant(sparse_polynomial_t p, constant_t c) {
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

sparse_polynomial_t sparsePolynomialMul(sparse_polynomial_t p, sparse_polynomial_t q) {
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
  if (res->monomialCount == 0) {
    res->coeffs[res->monomialCount] = constantFromInt(0);
    res->monomialDegrees[res->monomialCount] = constantFromInt(0);
    (res->monomialCount)++;
  }
  /* Set the degree of the polynomial */
  res->deg = constantFromCopy(res->monomialDegrees[res->monomialCount - 1]);
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

sparse_polynomial_t sparsePolynomialNeg(sparse_polynomial_t p) {
  unsigned int i;
  sparse_polynomial_t res;

  if (p == NULL) return NULL;
  res = __sparsePolynomialAllocate();
  res->refCount = 1;
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

sparse_polynomial_t sparsePolynomialCompose(sparse_polynomial_t p, sparse_polynomial_t q) {
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
  res->monomialCount = 1;
  res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->coeffs[0] = constantFromCopy(c);
  res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
  res->monomialDegrees[0] = constantFromCopy(d);
  res->deg = constantFromCopy(d);
  return res;
}

void sparsePolynomialDiv(sparse_polynomial_t *quot, sparse_polynomial_t *rest, sparse_polynomial_t a, sparse_polynomial_t b) {
  sparse_polynomial_t q, r, rt, bt, t1, t2, qp;
  constant_t rc, rd, bc, bd, qc, qd;

  /* Handle stupid cases */
  if ((a == NULL) || (b == NULL)) {
    *quot = NULL;
    *rest = NULL;
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

sparse_polynomial_t sparsePolynomialPowUnsignedInt(sparse_polynomial_t p, unsigned int n) {
  sparse_polynomial_t res, t, tmp;
  constant_t nC;
  
  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* Handle easy cases */
  if (n == 0) return sparsePolynomialFromIntConstant(1);
  if (n == 1) return sparsePolynomialFromCopy(p);
  if (n == 2) return sparsePolynomialMul(p,p);
  
  /* Handle the case when the polynomial only has one monomial */
  if (p->monomialCount == 1) {
    nC = constantFromUnsignedInt(n);
    res = __sparsePolynomialAllocate();
    res->refCount = 1;
    res->monomialCount = 1;
    res->coeffs = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
    res->coeffs[0] = constantPow(p->coeffs[0],nC);
    res->monomialDegrees = (constant_t *) safeCalloc(res->monomialCount, sizeof(constant_t));
    res->monomialDegrees[0] = constantMul(p->monomialDegrees[0],nC);
    res->deg = constantFromCopy(res->monomialDegrees[0]);
    constantFree(nC);
    return res;
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

int sparsePolynomialPowConstant(sparse_polynomial_t *r, sparse_polynomial_t p, constant_t n) {
  unsigned int nI;
  constant_t k, l, m;
  sparse_polynomial_t res, a, b, t;

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
  if (p->monomialCount == 1) {
    res = __sparsePolynomialAllocate();
    res->refCount = 1;
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

int sparsePolynomialPow(sparse_polynomial_t *r, sparse_polynomial_t p, sparse_polynomial_t q) {
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

sparse_polynomial_t sparsePolynomialDeriv(sparse_polynomial_t p) {
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
  if (res->monomialCount == 0) {
    res->coeffs[res->monomialCount] = constantFromInt(0);
    res->monomialDegrees[res->monomialCount] = constantFromInt(0);
    (res->monomialCount)++;
  }
  /* Set the degree of the polynomial */
  res->deg = constantFromCopy(res->monomialDegrees[res->monomialCount - 1]);
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

int sparsePolynomialFromExpression(sparse_polynomial_t *r, node *p) {
  sparse_polynomial_t a, b, quot, rest;
  int res;

  /* Handle stupid inputs */
  if (p == NULL) return 0;  

  /* Try to decompose expressions built on c, x, +, -, *, /, -, ^ */
  switch (p->nodeType) {
  case MEMREF:
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
      sparsePolynomialDiv(&quot, &rest, a, b);
      if (sparsePolynomialIsConstantZero(rest, 0)) {
	*r = quot;
	res = 1;
      } else {
	sparsePolynomialFree(quot);
	res = 0;
      }
      sparsePolynomialFree(rest);
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

void sparsePolynomialGetDegree(mpz_t deg, sparse_polynomial_t p) {
  if (p == NULL) {
    mpz_set_si(deg, -1);
    return;
  }
  if (!tryConstantToMpz(deg, p->deg)) {
    mpz_set_si(deg, -1);
  }
}

int sparsePolynomialGetDegreeAsInt(sparse_polynomial_t p) {
  int deg;

  if (p == NULL) return -1;
  if (!tryConstantToInt(&deg, p->deg)) return -1;
  return deg;
}

node *sparsePolynomialGetIthCoefficient(sparse_polynomial_t p, mpz_t i) {
  constant_t ic, coeffsum, t;
  unsigned int j, k;
  node *res;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;

  /* The index must be non-negative */
  if (mpz_sgn(i) < 0) {
    return addMemRef(makeConstantInt(0));
  }

  /* Handle the strange case when p has no monomials */
  if (p->monomialCount == 0u) {
    return addMemRef(makeConstantInt(0));
  }

  /* Construct a constant from i */
  ic = constantFromMpz(i);

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

constant_t sparsePolynomialGetIthCoefficientAsConstantIntIndex(sparse_polynomial_t p, int i) {
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

int sparsePolynomialGetCoefficients(node ***coeffs, unsigned int *deg, sparse_polynomial_t p) {
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

void sparsePolynomialEvalMpfr(mpfr_t y, sparse_polynomial_t p, mpfr_t x) {
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

void sparsePolynomialEvalMpfi(sollya_mpfi_t y, sparse_polynomial_t p, sollya_mpfi_t x) {
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
      return addMemRef(makeMul(constantToExpression(p->coeffs[0]),
			       makeVariable()));
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
      res = makeMul(constantToExpression(p->coeffs[0]),
		    makeVariable());
    } else {
      res = makeMul(constantToExpression(p->coeffs[0]),
		    makePow(makeVariable(),
			    constantToExpression(p->monomialDegrees[0])));
    }
  }
  for (i=1u;i<p->monomialCount;i++) {
    if (constantIsZero(p->monomialDegrees[i],0)) {
      temp = constantToExpression(p->coeffs[i]);
    } else {
      if (constantIsOne(p->monomialDegrees[i],0)) {
	temp = makeMul(constantToExpression(p->coeffs[i]),
		      makeVariable());
      } else {
	temp = makeMul(constantToExpression(p->coeffs[i]),
		       makePow(makeVariable(),
			       constantToExpression(p->monomialDegrees[i])));
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
  
  /* Perform Horner "evaluation" of p */
  res = constantToExpression(p->coeffs[p->monomialCount-1u]);
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
	  res = makeMul(makeVariable(), res);
	} else {
	  res = makeMul(makePow(makeVariable(),makeConstantInt(d)),res);
	}
      }
    } else {
      dc = constantSub(p->monomialDegrees[i], p->monomialDegrees[i-1u]);
      if (!constantIsZero(dc,0)) {
	if (constantIsOne(dc,0)) {
	  res = makeMul(makeVariable(), res);
	} else {
	  res = makeMul(makePow(makeVariable(),constantToExpression(dc)),res);
	}
      }
      constantFree(dc);
    }
    /* y <- p->coeffs[i-1u] + y */
    res = makeAdd(constantToExpression(p->coeffs[i-1u]), res);
  }
  /* y <- x^(p->monomialDegrees[0u]) * y */
  if (tryConstantToUnsignedInt(&a, p->monomialDegrees[0u])) {
    if (a != 0u) {
      if (a == 1u) {
	res = makeMul(makeVariable(), res);
      } else {
	res = makeMul(makePow(makeVariable(),
			      makeConstantInt(a)),
		      res);
      }
    }
  } else {
    if (!constantIsZero(p->monomialDegrees[0u],0)) {
      if (constantIsOne(p->monomialDegrees[0u],0)) {
	res = makeMul(makeVariable(), res);
      } else {
	res = makeMul(makePow(makeVariable(),
			      constantToExpression(p->monomialDegrees[0u])),
		      res);
      }
    }
  }
  
  /* Return the result */
  return addMemRef(res);
}

node *sparsePolynomialGetExpression(sparse_polynomial_t p, int canonical) {
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

unsigned int sparsePolynomialGetMonomialCount(sparse_polynomial_t p) {
  
  /* Handle stupid input */
  if (p == NULL) return 0u;

  /* Return the number of monomials */
  return p->monomialCount;
}

/* End of part for sparse polynomials */

/* Start of part for general (composed) polynomials */

static inline polynomial_t __polynomialAllocate() {
  return (polynomial_t) safeMalloc(sizeof(struct __polynomial_struct_t));
}

static inline void __polynomialFreeMem(polynomial_t p) {
  safeFree(p);
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

  /* Return the newly built polynomial */
  return res;
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
    __polynomialSparsify(p->value.pair.g);
    __polynomialSparsify(p->value.pair.h);
    sp = sparsePolynomialCompose(p->value.pair.g->value.sparse,
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
    p->value.sparse = sp;
    break;
  }
  p->type = SPARSE;
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

int polynomialEqual(polynomial_t p, polynomial_t q, int defVal) {

  /* Handle stupid inputs */
  if (p == NULL) return defVal;
  if (q == NULL) return defVal;

  /* Pointer equality */
  if (p == q) return 1;

  /* General case
     
     Can be optimized.
  
  */
  __polynomialSparsify(p);
  __polynomialSparsify(q);
  return sparsePolynomialEqual(p->value.sparse, q->value.sparse, defVal);
}

void polynomialDiv(polynomial_t *quot, polynomial_t *rest, polynomial_t a, polynomial_t b) {
  sparse_polynomial_t sq, sr;

  /* Handle stupid inputs */
  if ((a == NULL) || (b == NULL)) {
    *quot = NULL;
    *rest = NULL;
    return;
  }

  /* General case 

     Some cases can be optimized.

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
    return;
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

node *polynomialGetIthCoefficient(polynomial_t p, mpz_t i) {
  mpz_t deg;

  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* Handle case when i < 0 */
  if (mpz_sgn(i) < 0) {
    return addMemRef(makeConstantInt(0));
  }

  /* Handle case when i > degree */
  mpz_init(deg);
  polynomialGetDegree(deg, p);
  if ((mpz_sgn(deg) >= 0) && 
      (mpz_cmp(i,deg) > 0)) {
    mpz_clear(deg);
    return addMemRef(makeConstantInt(0));
  }
  
  /* General case */
  __polynomialSparsify(p);
  return sparsePolynomialGetIthCoefficient(p->value.sparse, i);
}

node *polynomialGetIthCoefficientIntIndex(polynomial_t p, int i) {
  int deg;

  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* Handle case when i < 0 */
  if (i < 0) {
    return addMemRef(makeConstantInt(0));
  }

  /* Handle case when i > degree */
  deg = polynomialGetDegreeAsInt(p);
  if ((deg >= 0) && 
      (i > deg)) {
    return addMemRef(makeConstantInt(0));
  }
  
  /* General case */
  __polynomialSparsify(p);
  return sparsePolynomialGetIthCoefficientIntIndex(p->value.sparse, i);
}

static inline constant_t __polynomialGetIthCoefficientAsConstantIntIndex(polynomial_t p, int i) {
  int deg;

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
  __polynomialSparsify(p);
  return sparsePolynomialGetIthCoefficientAsConstantIntIndex(p->value.sparse, i);
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
    __polynomialSparsify(p);
    return sparsePolynomialGetExpression(p->value.sparse, 0);
  }
  
  /* General case */
  switch (p->type) {
  case SPARSE:
    return sparsePolynomialGetExpression(p->value.sparse, 0);
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
  }

  /* Cannot be reached */
  return NULL;
}

node *polynomialGetExpression(polynomial_t p) {

  /* Handle stupid input */
  if (p == NULL) return NULL;

  /* If we can output the polynomial in any form, we do so. */
  if (p->outputType == ANY_FORM) 
    return __polynomialGetExpressionAnyForm(p);

  /* Here, we have to output a "hornerized" or "canonicalized"
     form. 
  */
  __polynomialSparsify(p);
  return sparsePolynomialGetExpression(p->value.sparse, 
				       (p->outputType == CANONICAL_FORM));
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
  t = polynomialGetExpression(p);
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
  t = polynomialGetExpression(p);
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

  /* Here, we have an addition, subtraction, multiplication or power
     form to handle and we may use y as a scratch space.

     We need another scratch variable. 

  */
  mpfr_init2(scratch, mpfr_get_prec(y));
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

  /* Here, we have an addition, subtraction, multiplication or power
     form to handle and we may use y as a scratch space.

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
  }
  if (usingReused) {
    returnReusedGlobalMPIVars(1);
  } else {
    sollya_mpfi_clear(scr);
  }
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
    /* We don't know easily */
    return 0;
    break;
  case MULTIPLICATION:
    return (__polynomialCheapCheckConstantZero(p->value.pair.g) ||
	    __polynomialCheapCheckConstantZero(p->value.pair.h));
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
    /* We don't know easily */
    return 0;
    break;
  case MULTIPLICATION:
    return (__polynomialCheapCheckConstantOne(p->value.pair.g) &&
	    __polynomialCheapCheckConstantOne(p->value.pair.h));
    break;
  case NEGATE:
    /* We don't know easily */
    return 0;    
    break;
  case POWER:
    if (constantIsZero(p->value.powering.c,0)) return 1;
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
  if (__polynomialCheapCheckConstantZero(p)) 
    return polynomialFromCopy(q);
  if (__polynomialCheapCheckConstantZero(q)) 
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
    res->type = SUBTRACTION;
    res->outputType = ANY_FORM;
    res->value.pair.g = polynomialFromCopy(p->value.pair.h);
    res->value.pair.h = polynomialFromCopy(p->value.pair.g);
    return res;
    break;
  case NEGATE:
    return polynomialFromCopy(p->value.g);
    break;
  }

  /* General case: construct the negation polynomial */
  res = __polynomialAllocate();
  res->refCount = 1u;
  res->type = NEGATION;
  res->outputType = ANY_FORM;
  res->value.g = polynomialFromCopy(p);
  
  /* Return the result polynomial */
  return res;
}

polynomial_t polynomialCompose(polynomial_t p, polynomial_t q) {
  polynomial_t res;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;
  if (q == NULL) return NULL;
    
  /* If both polynomials are sparse, and polynomial q only has one
     monomial, perform the operation on sparse polynomials.
  */
  if ((p->type == SPARSE) && (q->type == SPARSE)) {
    if (sparsePolynomialGetMonomialCount(q->value.sparse) <= 1u) {
      return __polynomialBuildFromSparse(sparsePolynomialCompose(p->value.sparse,
								 q->value.sparse));
    }
  }

  /* General case: construct the composed polynomial */
  res = __polynomialAllocate();
  res->refCount = 1u;
  res->type = COMPOSITION;
  res->outputType = ANY_FORM;
  res->value.pair.g = polynomialFromCopy(p);
  res->value.pair.h = polynomialFromCopy(q);
  
  /* Return the result polynomial */
  return res;
}

polynomial_t polynomialDeriv(polynomial_t p) {
  polynomial_t res;

  /* Handle stupid inputs */
  if (p == NULL) return NULL;

  /* Handle certain cases in an ad-hoc way */
  switch (p->type) {
  case SPARSE:
    return __polynomialBuildFromSparse(sparsePolynomialDeriv(p->value.sparse));
    break;
  case ADDITION:
  case SUBTRACTION:
    res = __polynomialAllocate();
    res->refCount = 1u;
    res->type = p->type;
    res->outputType = ANY_FORM;
    res->value.pair.g = polynomialDeriv(p->value.pair.h);
    res->value.pair.h = polynomialDeriv(p->value.pair.g);
    return res;
    break;
  case NEGATE:
    res = __polynomialAllocate();
    res->refCount = 1u;
    res->type = NEGATE;
    res->outputType = ANY_FORM;
    res->value.g = polynomialDeriv(p->value.g);
    return res;
    break;
  }
  
  /* General case: sparsify the polynomial and return the derivative
     of the sparse polynomial.
  */
  __polynomialSparsify(p);  
  return __polynomialBuildFromSparse(sparsePolynomialDeriv(p->value.sparse));
}

int polynomialFromExpression(polynomial_t *r, node *p) {
  polynomial_t a, b, quot, rest, zero;
  sparse_polynomial_t sr;
  int res;

  /* Handle stupid inputs */
  if (p == NULL) return 0;  

  /* Try to decompose expressions built on c, x, +, -, *, /, -, ^ */
  switch (p->nodeType) {
  case MEMREF:
    return polynomialFromExpression(r, getMemRefChild(p));
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
    if (!polynomialFromExpression(&a, p->child1)) 
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
    if (!polynomialFromExpression(&a, p->child1)) 
      break;
    if (!polynomialFromExpression(&b, p->child2)) {
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
	*r = quot;
	res = 1;
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
  if (isConstant(p)) {
    if (!__sparsePolynomialFromConstantExpression(&sr, p)) return 0;
    *r = __polynomialBuildFromSparse(sr);
    return 1;
  }

  return 0;
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

  /* General case: construct the powering polynomial */  
  res = __polynomialAllocate();
  res->refCount = 1u;
  res->type = POWER;
  res->outputType = ANY_FORM;
  res->value.powering.g = polynomialFromCopy(p);
  res->value.powering.c = constantFromUnsignedInt(n);
  
  /* Return the result polynomial */
  return res;
}

int polynomialPow(polynomial_t *r, polynomial_t p, polynomial_t q) {
  constant_t n;
  int deg;
  polynomial_t res;

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

  /* General case: construct the powering polynomial */
  res = __polynomialAllocate();
  res->refCount = 1u;
  res->type = POWER;
  res->outputType = ANY_FORM;
  res->value.powering.g = polynomialFromCopy(p);
  res->value.powering.c = n;
  
  /* Set the output and return success */
  *r = res;
  return 1;
}

