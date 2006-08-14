#include <gmp.h>
#include <mpfr.h>
#include <stdio.h> /* fprintf, fopen, fclose, */
#include <errno.h>
#include <inttypes.h>
#include <stdlib.h>
#include "expression.h"
#include "double.h"
#include "main.h"


typedef union {
  int32_t i[2]; 
  double d;
} db_number;




int mpfr_round_to_double(mpfr_t rop, mpfr_t op) {
  double d;
  int res;

  if (mpfr_get_prec(op) < 53) {
    printMessage(1,"Warning: rounding a value computed on less than 53 bits to double precision.\n");
  }

  d = mpfr_get_d(op,GMP_RNDN);
  if (mpfr_set_d(rop,d,GMP_RNDN) != 0) {
    printMessage(1,"Warning: double rounding occured on invoking the double precision rounding operator.\n");
    printMessage(1,"Try to increase the working precision.\n");
  }
  
  res = mpfr_cmp(rop,op);

  return res;
}

int mpfr_round_to_doubledouble(mpfr_t rop, mpfr_t op) {
  double d;
  mpfr_t accu, temp, rest;
  mp_prec_t prec;
  int res;

  prec = mpfr_get_prec(op);
  if (prec < 106) {
    printMessage(1,"Warning: rounding a value computed on less than 106 bits to double-double precision.\n");
    prec = 106;
  }

  mpfr_init2(accu, prec);
  mpfr_init2(temp, prec);
  mpfr_init2(rest, prec);

  d = mpfr_get_d(op,GMP_RNDN);
  if (mpfr_set_d(accu,d,GMP_RNDN) != 0) {
    printMessage(1,"Warning: double rounding occured on invoking the double-double rounding operator.\n");
    printMessage(1,"The rounding occured on recasting to MPFR. This should not occur.\n");
  }
  if (mpfr_sub(rest,op,accu,GMP_RNDN) != 0) {
    printMessage(1,"Warning: double rounding occured on invoking the double-double rounding operator.\n");
    printMessage(1,"The rounding occured on substracting in MPFR. This should not occur.\n");
  }
  d = mpfr_get_d(rest,GMP_RNDN);
  if (mpfr_set_d(temp,d,GMP_RNDN) != 0) {
    printMessage(1,"Warning: double rounding occured on invoking the double-double rounding operator.\n");
    printMessage(1,"The rounding occured on recasting to MPFR. This should not occur.\n");
  }
  if (mpfr_add(accu,accu,temp,GMP_RNDN) != 0) {
    printMessage(1,"Warning: double rounding occured on invoking the double-double rounding operator.\n");
    printMessage(1,"The rounding occured on substracting in MPFR. This should not occur.\n");
  }
  if (mpfr_set(rop,accu,GMP_RNDN) != 0) {
    printMessage(1,"Warning: double rounding occured on invoking the double-double rounding operator.\n");
    printMessage(1,"Try to increase the working precision.\n");
  }

  res = mpfr_cmp(rop,op);

  mpfr_clear(accu);
  mpfr_clear(temp);
  mpfr_clear(rest);
  return res;
}

int mpfr_round_to_tripledouble(mpfr_t rop, mpfr_t op) {
  double d;
  mpfr_t accu, temp, rest;
  mp_prec_t prec;
  int res;

  prec = mpfr_get_prec(op);
  if (prec < 159) {
    printMessage(1,"Warning: rounding a value computed on less than 159 bits to triple-double precision\n");
    prec = 159;
  }

  mpfr_init2(accu, prec);
  mpfr_init2(temp, prec);
  mpfr_init2(rest, prec);

  d = mpfr_get_d(op,GMP_RNDN);
  if (mpfr_set_d(accu,d,GMP_RNDN) != 0) {
    printMessage(1,"Warning: double rounding occured on invoking the triple-double rounding operator.\n");
    printMessage(1,"The rounding occured on recasting to MPFR. This should not occur.\n");
  }
  if (mpfr_sub(rest,op,accu,GMP_RNDN) != 0) {
    printMessage(1,"Warning: double rounding occured on invoking the triple-double rounding operator.\n");
    printMessage(1,"The rounding occured on substracting in MPFR. This should not occur.\n");
  }
  d = mpfr_get_d(rest,GMP_RNDN);
  if (mpfr_set_d(temp,d,GMP_RNDN) != 0) {
    printMessage(1,"Warning: double rounding occured on invoking the triple-double rounding operator.\n");
    printMessage(1,"The rounding occured on recasting to MPFR. This should not occur.\n");
  }
  if (mpfr_add(accu,accu,temp,GMP_RNDN) != 0) {
    printMessage(1,"Warning: double rounding occured on invoking the triple-double rounding operator.\n");
    printMessage(1,"The rounding occured on substracting in MPFR. This should not occur.\n");
  }
  if (mpfr_sub(rest,op,accu,GMP_RNDN) != 0) {
    printMessage(1,"Warning: double rounding occured on invoking the triple-double rounding operator.\n");
    printMessage(1,"The rounding occured on substracting in MPFR. This should not occur.\n");
  }
  d = mpfr_get_d(rest,GMP_RNDN);
  if (mpfr_set_d(temp,d,GMP_RNDN) != 0) {
    printMessage(1,"Warning: double rounding occured on invoking the triple-double rounding operator.\n");
    printMessage(1,"The rounding occured on recasting to MPFR. This should not occur.\n");
  }
  if (mpfr_add(accu,accu,temp,GMP_RNDN) != 0) {
    printMessage(1,"Warning: double rounding occured on invoking the triple-double rounding operator.\n");
    printMessage(1,"The rounding occured on substracting in MPFR. This should not occur.\n");
  }
  if (mpfr_set(rop,accu,GMP_RNDN) != 0) {
    printMessage(1,"Warning: double rounding occured on invoking the triple-double rounding operator.\n");
    printMessage(1,"Try to increase the working precision.\n");
  }

  res = mpfr_cmp(rop,op);

  mpfr_clear(accu);
  mpfr_clear(temp);
  mpfr_clear(rest);

  return res;
}


int printDoubleInHexa(mpfr_t x) {
  int res;
  double d;
  mpfr_t temp;
  db_number xdb, endianessdb;

  mpfr_init2(temp,mpfr_get_prec(x));
  
  d = mpfr_get_d(x,GMP_RNDN);
  if (mpfr_set_d(temp,d,GMP_RNDN) != 0) {
    printMessage(1,"Warning: rounding occured unexpectedly on reconverting a double value.\n");
  }
  
  res = mpfr_cmp(temp,x);

  if (res) 
    printMessage(1,"Warning: rounding occured before printing a value as a double.\n");

  xdb.d = d;
  endianessdb.d = 1.0;
  if ((endianessdb.i[1] == 0x3ff00000) && (endianessdb.i[0] == 0)) {
    printf("0x%08x%08x\n",xdb.i[1],xdb.i[0]);
  } else {
    if ((endianessdb.i[0] == 0x3ff00000) && (endianessdb.i[1] == 0)) {
      printf("0x%08x%08x\n",xdb.i[0],xdb.i[1]);
    } else {
      printMessage(1,"Warning: could not figure out the endianess of the system. Will print 1.0 instead of the value.\n");
      printf("0x3ff0000000000000\n");
    }
  }

  mpfr_clear(temp);
  return res;
}


int readHexa(mpfr_t res, char *c) {
  int ret, msb, lsb, i;
  double x;
  char msbstr[9], lsbstr[9];
  db_number xdb, endianessdb;
  
  x = 1.0;
  c += 2; /* Skip over "0x" */
  for (i=0;i<9;i++) {
    msbstr[i] = '\0';
    lsbstr[i] = '\0';
  }
  for (i=0;i<8;i++) {
    msbstr[i] = *c;
    c++;
  }
  for (i=0;i<8;i++) {
    lsbstr[i] = *c;
    c++;
  }

  msb = strtol(msbstr,NULL,16);
  lsb = strtol(lsbstr,NULL,16);

  endianessdb.d = 1.0;
  if ((endianessdb.i[1] == 0x3ff00000) && (endianessdb.i[0] == 0)) {
    xdb.i[1] = msb;
    xdb.i[0] = lsb;
  } else {
    if ((endianessdb.i[0] == 0x3ff00000) && (endianessdb.i[1] == 0)) {
      xdb.i[0] = msb;
      xdb.i[1] = lsb;
    } else {
      printMessage(1,"Warning: could not figure out the endianess of the system. Will read 1.0 instead of the value.\n");
      xdb.d = 1.0;
    }
  }

  x = xdb.d;

  if (mpfr_set_d(res,x,GMP_RNDN) != 0) ret = 0; else ret = 1;
  return ret;
}
