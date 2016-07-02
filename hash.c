/*

  Copyright 2015-2016 by

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

#include <stdint.h>
#include <stdlib.h>
#include "hash.h"

static inline uint64_t __rotateLeft(uint64_t x, unsigned int k) {
  unsigned int s, r;
  uint64_t tmp1, tmp2, tmp3;

  s = k & ((unsigned int) 63);
  if (s == ((unsigned int) 0)) return x;
  r = 64 - s;
  tmp1 = x;
  tmp1 <<= s;
  tmp2 = x;
  tmp2 >>= r;
  tmp3 = tmp1 | tmp2;
  return tmp3;
}

static inline uint64_t __hashCombine_internal(uint64_t h1, uint64_t h2) {
  uint64_t tmp1, tmp2;
  tmp1 = h1 ^ h2;
  tmp2 = __rotateLeft(tmp1, 17);
  return tmp2;
}

static inline uint64_t __hashUint64_internal(uint64_t x) {
  uint64_t tmp1, tmp2;
  tmp1 = x ^ UINT64_C(0xcafebabedeadbeef);
  tmp2 = __rotateLeft(tmp1, 13);
  return tmp2;
}

static inline uint64_t __hashUint8_internal(uint8_t x) {
  uint64_t tmp1, tmp2;
  tmp1 = (uint64_t) x;
  tmp2 = tmp1; tmp1 <<= 8;
  tmp2 |= tmp1; tmp1 <<= 8;
  tmp2 |= tmp1; tmp1 <<= 8;
  tmp2 |= tmp1; tmp1 <<= 8;
  tmp2 |= tmp1; tmp1 <<= 8;
  tmp2 |= tmp1; tmp1 <<= 8;
  tmp2 |= tmp1; tmp1 <<= 8;
  tmp2 |= tmp1;
  return __hashUint64_internal(tmp2);
}

static inline uint64_t __hashUint16_internal(uint16_t x) {
  uint64_t tmp1, tmp2;
  tmp1 = (uint64_t) x;
  tmp2 = tmp1; tmp1 <<= 16;
  tmp2 |= tmp1; tmp1 <<= 16;
  tmp2 |= tmp1; tmp1 <<= 16;
  tmp2 |= tmp1;
  return __hashUint64_internal(tmp2);
}

static inline uint64_t __hashUint32_internal(uint32_t x) {
  uint64_t tmp1, tmp2;
  tmp1 = (uint64_t) x;
  tmp2 = tmp1; tmp1 <<= 32;
  tmp2 |= tmp1;
  return __hashUint64_internal(tmp2);
}

static inline uint64_t __hashUnsignedChar_internal(unsigned char x) {
  uint8_t tmp1;
  uint16_t tmp2;
  uint32_t tmp3;
  uint64_t tmp4;
  if (sizeof(uint8_t) >= sizeof(unsigned char)) {
    tmp1 = (uint8_t) x;
    return __hashUint8_internal(tmp1);
  }
  if (sizeof(uint16_t) >= sizeof(unsigned char)) {
    tmp2 = (uint16_t) x;
    return __hashUint16_internal(tmp2);
  }
  if (sizeof(uint32_t) >= sizeof(unsigned char)) {
    tmp3 = (uint32_t) x;
    return __hashUint32_internal(tmp3);
  }
  tmp4 = (uint64_t) x;
  return __hashUint64_internal(tmp4);
}

static inline uint64_t __hashUnsignedInt_internal(unsigned int x) {
  uint8_t tmp1;
  uint16_t tmp2;
  uint32_t tmp3;
  uint64_t tmp4;
  if (sizeof(uint8_t) >= sizeof(unsigned int)) {
    tmp1 = (uint8_t) x;
    return __hashUint8_internal(tmp1);
  }
  if (sizeof(uint16_t) >= sizeof(unsigned int)) {
    tmp2 = (uint16_t) x;
    return __hashUint16_internal(tmp2);
  }
  if (sizeof(uint32_t) >= sizeof(unsigned int)) {
    tmp3 = (uint32_t) x;
    return __hashUint32_internal(tmp3);
  }
  tmp4 = (uint64_t) x;
  return __hashUint64_internal(tmp4);
}

static inline uint64_t __hashUnsignedLong_internal(unsigned long x) {
  uint8_t tmp1;
  uint16_t tmp2;
  uint32_t tmp3;
  uint64_t tmp4;
  if (sizeof(uint8_t) >= sizeof(unsigned long)) {
    tmp1 = (uint8_t) x;
    return __hashUint8_internal(tmp1);
  }
  if (sizeof(uint16_t) >= sizeof(unsigned long)) {
    tmp2 = (uint16_t) x;
    return __hashUint16_internal(tmp2);
  }
  if (sizeof(uint32_t) >= sizeof(unsigned long)) {
    tmp3 = (uint32_t) x;
    return __hashUint32_internal(tmp3);
  }
  tmp4 = (uint64_t) x;
  return __hashUint64_internal(tmp4);
}

static inline uint64_t __hashPointer_internal(void *x) {
  uintptr_t X;
  uint8_t tmp1;
  uint16_t tmp2;
  uint32_t tmp3;
  uint64_t tmp4;
  X = (uintptr_t) x;
  if (sizeof(uint8_t) >= sizeof(uintptr_t)) {
    tmp1 = (uint8_t) X;
    return __hashUint8_internal(tmp1);
  }
  if (sizeof(uint16_t) >= sizeof(uintptr_t)) {
    tmp2 = (uint16_t) X;
    return __hashUint16_internal(tmp2);
  }
  if (sizeof(uint32_t) >= sizeof(uintptr_t)) {
    tmp3 = (uint32_t) X;
    return __hashUint32_internal(tmp3);
  }
  tmp4 = (uint64_t) X;
  return __hashUint64_internal(tmp4);
}

static inline uint64_t __hashChar_internal(char x) {
  char xx;
  unsigned char tmp;
  xx = x;
  tmp = *((unsigned char *) &xx);
  return __hashUnsignedChar_internal(tmp);
}

static inline uint64_t __hashInt_internal(int x) {
  int xx;
  unsigned int tmp;
  xx = x;
  tmp = *((unsigned int *) &xx);
  return __hashUnsignedInt_internal(tmp);
}

static inline uint64_t __hashInt64_internal(int64_t x) {
  int64_t xx;
  uint64_t tmp;
  xx = x;
  tmp = *((uint64_t *) &xx);
  return __hashUint64_internal(tmp);
}

static inline uint64_t __hashLong_internal(long x) {
  long xx;
  unsigned long tmp;
  xx = x;
  tmp = *((unsigned long *) &xx);
  return __hashUnsignedLong_internal(tmp);
}

typedef union {
  uint64_t l;
  double d;
} binary64_caster;

static inline uint64_t __hashDouble_internal(double x) {
  double xx;
  uint64_t tmp;
  binary64_caster caster;
  xx = x;
  caster.d = xx;
  tmp = caster.l;
  return __hashUint64_internal(tmp);
}

static inline uint64_t __hashMpfr_internal(mpfr_t x) {
  long exp;
  double mant;
  uint64_t tmp1, tmp2;

  exp = (long) 0;
  mant = mpfr_get_d_2exp(&exp, x, GMP_RNDN);

  tmp1 = __hashLong_internal(exp);
  tmp2 = __hashDouble_internal(mant);
  return __hashCombine_internal(tmp1, tmp2);
}

static inline uint64_t __hashMpfi_internal(sollya_mpfi_t x) {
  uint64_t tmp1, tmp2;
  
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  tmp1 = __hashMpfr_internal(&(x->left));
  tmp2 = __hashMpfr_internal(&(x->right));
  return __hashCombine_internal(tmp1, tmp2);
}

static inline uint64_t __hashMpz_internal(mpz_t x) {
  long exp;
  double mant;
  uint64_t tmp1, tmp2;

  exp = (long) 0;
  mant = mpz_get_d_2exp(&exp, x);
  tmp1 = __hashLong_internal(exp);
  tmp2 = __hashDouble_internal(mant);
  return __hashCombine_internal(tmp1, tmp2);
}

static inline uint64_t __hashMpq_internal(mpq_t x) {
  uint64_t tmp1, tmp2;
  
  mpq_canonicalize(x);
  tmp1 = __hashMpz_internal(mpq_numref(x));
  tmp2 = __hashMpz_internal(mpq_denref(x));
  return __hashCombine_internal(tmp1, tmp2);  
}

static inline uint64_t __hashString_internal(char *x) {
  uint64_t tmp;
  char *curr;
  
  tmp = UINT64_C(0);
  for (curr=x;*curr!='\0';curr++) {
    tmp = __hashCombine_internal(tmp, __hashChar_internal(*curr));
  }
  return tmp;
}

uint64_t hashCombine(uint64_t h1, uint64_t h2) {
  uint64_t res;
  res = __hashCombine_internal(h1, h2);
  return res;
}

uint64_t hashMpfr(mpfr_t x) {
  uint64_t res;
  res = __hashMpfr_internal(x);
  return res;
}

uint64_t hashMpfi(sollya_mpfi_t x) {
  uint64_t res;
  res = __hashMpfi_internal(x);
  return res;
}

uint64_t hashMpq(mpq_t x) {
  uint64_t res;
  res = __hashMpq_internal(x);
  return res;
}

uint64_t hashString(char *x) {
  uint64_t res;
  res = __hashString_internal(x);
  return res;
}

uint64_t hashChar(char x) {
  uint64_t res;
  res = __hashChar_internal(x);
  return res;
}

uint64_t hashInt(int x) {
  uint64_t res;
  res = __hashInt_internal(x);
  return res;
}

uint64_t hashUnsignedInt(unsigned int x) {
  uint64_t res;
  res = __hashUnsignedInt_internal(x);
  return res;
}

uint64_t hashInt64(int64_t x) {
  uint64_t res;
  res = __hashInt64_internal(x);
  return res;
}

uint64_t hashPointer(void *x) {
  uint64_t res;
  res = __hashPointer_internal(x);
  return res;
}


