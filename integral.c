#include <mpfr.h>
#include <mpfi.h>
#include "expression.h"
#include "infnorm.h"
#include "integral.h"

#include <stdio.h> /* fprintf, fopen, fclose, */
#include <stdlib.h> /* exit, free, mktemp */
#include <errno.h>

#define DEBUG 0
#define DEBUGMPFI 0



mpfi_t integral(node *func, rangetype interval, mp_prec_t prec, mpfr_t diam) {
  rangetype x,y;
  mpfr_t x1,x2,y1,y2;
  mpfi_t temp, val;

  rangetype sum;
  sum.a = (mpfr_t*) malloc(sizeof(mpfr_t));
  sum.b = (mpfr_t*) malloc(sizeof(mpfr_t));
  mpfr_init2(*(sum.a),prec);
  mpfr_init2(*(sum.b),prec);
  mpfr_set_d(*(sum.a),0.,GMP_RNDD);
  mpfr_set_d(*(sum.b),0.,GMP_RNDU);
  
  
  if (mpfr_equal_p(*(interval.a),*(interval.b))) {
    printf("Warning: the given interval is reduced to one point.\n");
    return sum;
  }

  if (mpfr_less_p(*(interval.b),*(interval.a))) {
    printf("Error: the interval is empty.\n");
    return sum;
  }

  mpfi_init2(temp,prec);
  mpfi_init2(val,prec);

  mpfr_init2(x1,prec);
  mpfr_init2(x2,prec);
  mpfr_set(x1, *(interval.a),GMP_RNDD);
  mpfr_add(x2, x1, diam, GMP_RNDN);
  x.a = &x1;
  x.b = &x2;
  
  mpfr_init2(y1,prec);
  mpfr_init2(y2,prec);
  y.a = &y1;
  y.b = &y2;


  while(mpfr_less_p(x2,*(interval.b))) {
    evaluateRangeFunction(y, func, x, prec);

    mpfi_set_fr(temp, x1);
    mpfi_set_fr(val, x2);
    mpfi_sub(temp, x2, x1);
    
    mpfi_interv_fr(val, *(y.a), *(y.b));
    mpfi_mul(temp, temp, val);
    
    mpfr_get_left(y1, temp);
    mpfr_get_right(y2, temp);
    mpfr_add(*(sum.a), *(sum.a), y1, GMP_RNDD);
    mpfr_add(*(sum.b), *(sum.b), y2, GMP_RNDU);
    
    mpfr_set(x1,x2,GMP_RNDD); // exact
    mpfr_add(x2, x1, diam, GMP_RNDN);
  }

  mpfr_set(x2,b,GMP_RNDU);
  evaluateRangeFunction(y, func, x, prec);

  mpfi_set_fr(temp, x1);
  mpfi_set_fr(val, x2);
  mpfi_sub(temp, x2, x1);
    
  mpfi_interv_fr(val, *(y.a), *(y.b));
  mpfi_mul(temp, temp, val);
    
  mpfr_get_left(y1, temp);
  mpfr_get_right(y2, temp);
  mpfr_add(*(sum.a), *(sum.a), y1, GMP_RNDD);
  mpfr_add(*(sum.b), *(sum.b), y2, GMP_RNDU);
  
 
  mpfi_clear(val); mpfi_clear(temp);
  mpfr_clear(x1); mpfi_clear(x2);  
  mpfr_clear(y1); mpfi_clear(y2);
  
  return sum;
}
