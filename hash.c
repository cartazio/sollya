/*

  Copyright 2015 by

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

#include <stdint.h>
#include <stdlib.h>
#include "hash.h"

uint64_t hashChar(char x) {
  char xx;
  unsigned char X;
  uint64_t tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
  xx = x;
  X = *((unsigned char *) &xx);
  tmp1 = (uint64_t) X;
  tmp2 = tmp1 << 8;
  tmp3 = tmp2 << 8;
  tmp4 = tmp3 << 8;
  tmp5 = tmp4 << 8;
  tmp6 = tmp5 << 8;
  tmp7 = tmp6 << 8;
  tmp8 = tmp7 << 8;
  return hashCombine(hashCombine(hashCombine(hashCombine(hashCombine(hashCombine(hashCombine(tmp1, tmp2), tmp3), tmp4), tmp5), tmp6), tmp7), tmp8);
}

uint64_t hashInt(int x) {
  int xx;
  unsigned int X;
  uint64_t tmp1, tmp2, tmp3, tmp4;
  xx = x;
  X = *((unsigned int *) &xx);
  tmp1 = (uint64_t) X;
  tmp2 = tmp1 << 16;
  tmp3 = tmp2 << 16;
  tmp4 = tmp3 << 16;
  return hashCombine(hashCombine(hashCombine(tmp1, tmp2), tmp3), tmp4);
}

uint64_t hashUnsignedInt(unsigned int x) {
  unsigned int X;
  uint64_t tmp1, tmp2, tmp3, tmp4;
  X = x;
  tmp1 = (uint64_t) X;
  tmp2 = tmp1 << 16;
  tmp3 = tmp2 << 16;
  tmp4 = tmp3 << 16;
  return hashCombine(hashCombine(hashCombine(tmp1, tmp2), tmp3), tmp4);
}

uint64_t hashInt64(int64_t x) {
  int64_t xx;
  uint64_t X;
  uint64_t tmp1;
  xx = x;
  X = *((uint64_t *) &xx);
  tmp1 = X;
  tmp1 ^= UINT64_C(0x0f0f0f0f0f0f0f0f);
  return tmp1;
}

uint64_t hashLong(long x) {
  long xx;
  unsigned long X;
  uint64_t tmp1, tmp2, tmp3, tmp4;
  xx = x;
  X = *((unsigned long *) &xx);
  tmp1 = (uint64_t) X;
  tmp2 = tmp1 << 16;
  tmp3 = tmp2 << 16;
  tmp4 = tmp3 << 16;
  return hashCombine(hashCombine(hashCombine(tmp1, tmp2), tmp3), tmp4);
}

uint64_t hashDouble(double x) {
  double xx;
  uint64_t tmp1;

  xx = x;
  tmp1 = *((uint64_t *) &xx);
  tmp1 ^= UINT64_C(0x00ff00ff00ff00ff);
  return tmp1;
}

uint64_t hashMpfr(mpfr_t x) {
  long exp;
  double mant;
  uint64_t tmp1, tmp2;

  exp = (long) 0;
  mant = mpfr_get_d_2exp(&exp, x, GMP_RNDN);

  tmp1 = hashLong(exp);
  tmp2 = hashDouble(mant);
  return hashCombine(tmp1, tmp2);
}

uint64_t hashMpfi(sollya_mpfi_t x) {
  uint64_t tmp1, tmp2;
  
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  tmp1 = hashMpfr(&(x->left));
  tmp2 = hashMpfr(&(x->right));
  return hashCombine(tmp1, tmp2);
}

uint64_t hashMpz(mpz_t x) {
  long exp;
  double mant;
  uint64_t tmp1, tmp2;

  exp = (long) 0;
  mant = mpz_get_d_2exp(&exp, x);
  tmp1 = hashLong(exp);
  tmp2 = hashDouble(mant);
  return hashCombine(tmp1, tmp2);
}

uint64_t hashMpq(mpq_t x) {
  uint64_t tmp1, tmp2;
  
  mpq_canonicalize(x);
  tmp1 = hashMpz(mpq_numref(x));
  tmp2 = hashMpz(mpq_denref(x));
  return hashCombine(tmp1, tmp2);  
}

uint64_t hashString(char *x) {
  uint64_t tmp;
  char *curr;
  
  tmp = UINT64_C(0xcafebabedeadbeef);
  for (curr=x;*curr!='\0';curr++) {
    tmp = hashCombine(tmp, hashChar(*curr));
  }
  return tmp;
}

uint64_t hashPointer(void *x) {
  uint64_t tmp1, tmp2, tmp3, tmp4;
  tmp1 = (uint64_t) x;
  tmp2 = tmp1 << 16;
  tmp3 = tmp2 << 16;
  tmp4 = tmp3 << 16;
  return hashCombine(hashCombine(hashCombine(tmp1, tmp2), tmp3), tmp4);
}

uint64_t hashCombine(uint64_t h1, uint64_t h2) {
  uint64_t tmp;
  tmp = h1 ^ h2;
  tmp ^= UINT64_C(0xff00ff00ff00ff00);
  return tmp;
}

