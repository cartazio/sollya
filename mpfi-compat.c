/*

  Copyright 2007-2014 by

  Laboratoire de l'Informatique du Parallelisme,
  UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668,

  Laboratoire d'Informatique de Paris 6, equipe PEQUAN,
  UPMC Universite Paris 06 - CNRS - UMR 7606 - LIP6, Paris, France,

  LORIA (CNRS, INPL, INRIA, UHP, U-Nancy 2)

  and by

  Centre de recherche INRIA Sophia-Antipolis Mediterranee, equipe APICS,
  Sophia Antipolis, France.

  Contributors Ch. Lauter, S. Chevillard, M. Joldes

  christoph.lauter@ens-lyon.org
  sylvain.chevillard@ens-lyon.org
  joldes@laas.fr

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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mpfi-compat.h"
#include "general.h"
#include "sollya-messaging.h"
#include <stdlib.h>
#include <stdio.h>

/* Functions that handle non regular intervals */

int sollya_mpfi_has_nan(sollya_mpfi_t op) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  return ( mpfr_nan_p(&(op->left))
           || mpfr_nan_p(&(op->right))
           );
}

static inline int sollya_mpfi_has_nan_opt(sollya_mpfi_t op) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  return ( mpfr_nan_p(&(op->left))
           || mpfr_nan_p(&(op->right))
           );
}

int sollya_mpfi_nan_p(sollya_mpfi_t op) {
  return sollya_mpfi_has_nan_opt(op);
}

int sollya_mpfi_set_nan(sollya_mpfi_t rop) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  mpfr_set_nan(&(rop->left));
  mpfr_set_nan(&(rop->right));
  return MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
}

static inline int sollya_mpfi_set_nan_opt(sollya_mpfi_t rop) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  mpfr_set_nan(&(rop->left));
  mpfr_set_nan(&(rop->right));
  return MPFI_FLAGS_BOTH_ENDPOINTS_INEXACT;
}

void sollya_mpfi_nan_normalize(sollya_mpfi_t rop) {
  if (sollya_mpfi_has_nan_opt(rop)) sollya_mpfi_set_nan_opt(rop);
}

static inline void sollya_mpfi_nan_normalize_opt(sollya_mpfi_t rop) {
  if (sollya_mpfi_has_nan_opt(rop)) sollya_mpfi_set_nan_opt(rop);
}

void sollya_mpfi_zero_sign_normalize(sollya_mpfi_t op) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  if (mpfr_zero_p(&(op->left))) {
    mpfr_mul(&(op->left), &(op->left), &(op->left), GMP_RNDN); /* (+/- 0)^2 = +0 */
  } 
  if (mpfr_zero_p(&(op->right))) { 
    mpfr_mul(&(op->right), &(op->right), &(op->right), GMP_RNDN); /* (+/- 0)^2 = +0 */   
    mpfr_neg(&(op->right), &(op->right), GMP_RNDN); /* - (+0) = -0 */
  }
}

int sollya_mpfi_is_empty(sollya_mpfi_t op) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  if (sollya_mpfi_has_nan_opt(op))
    return 0;
  else return mpfr_greater_p(&(op->left), &(op->right)); 

}

static inline int sollya_mpfi_is_empty_opt(sollya_mpfi_t op) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  if (sollya_mpfi_has_nan_opt(op))
    return 0;
  else return mpfr_greater_p(&(op->left), &(op->right)); 

}

int sollya_mpfi_set_empty(sollya_mpfi_t rop) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  mpfr_set_inf(&(rop->left),1);
  mpfr_set_inf(&(rop->right),-1);
  return MPFI_FLAGS_BOTH_ENDPOINTS_EXACT;
}

static inline int sollya_mpfi_set_empty_opt(sollya_mpfi_t rop) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  mpfr_set_inf(&(rop->left),1);
  mpfr_set_inf(&(rop->right),-1);
  return MPFI_FLAGS_BOTH_ENDPOINTS_EXACT;
}

void sollya_mpfi_empty_normalize(sollya_mpfi_t rop) {
  if (sollya_mpfi_is_empty_opt(rop)) sollya_mpfi_set_empty_opt(rop);
}

static inline void sollya_mpfi_empty_normalize_opt(sollya_mpfi_t rop) {
  if (sollya_mpfi_is_empty_opt(rop)) sollya_mpfi_set_empty_opt(rop);
}

int sollya_mpfi_set_full_range(sollya_mpfi_t rop) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  mpfr_set_inf(&(rop->left),-1);
  mpfr_set_inf(&(rop->right),1);
  return MPFI_FLAGS_BOTH_ENDPOINTS_EXACT;
}

int sollya_mpfi_set_negative_inf(sollya_mpfi_t rop) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  mpfr_set_inf(&(rop->left),-1);
  mpfr_set_inf(&(rop->right),-1);
  return MPFI_FLAGS_BOTH_ENDPOINTS_EXACT;
}

int sollya_mpfi_set_positive_inf(sollya_mpfi_t rop) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  mpfr_set_inf(&(rop->left),1);
  mpfr_set_inf(&(rop->right),1);
  return MPFI_FLAGS_BOTH_ENDPOINTS_EXACT;
}

/* Check for infinities and zeros */

int sollya_mpfi_has_zero(sollya_mpfi_t op) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  if ( sollya_mpfi_has_nan_opt(op) || sollya_mpfi_is_empty_opt(op) )
    return 0;
  else return ( mpfr_sgn(&(op->left)) * mpfr_sgn(&(op->right)) <= 0 );
}

int sollya_mpfi_has_zero_inside(sollya_mpfi_t op) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  if ( sollya_mpfi_has_nan_opt(op) || sollya_mpfi_is_empty_opt(op) )
    return 0;
  else return ( mpfr_sgn(&(op->left)) * mpfr_sgn(&(op->right)) < 0 );
}

int sollya_mpfi_is_zero(sollya_mpfi_t op) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  if ( sollya_mpfi_has_nan_opt(op) || sollya_mpfi_is_empty_opt(op) )
    return 0;
  else return ( (mpfr_sgn(&(op->left))==0) && (mpfr_sgn(&(op->right)) == 0) );
}

int mpfr_is_positive_infinity(mpfr_t op) {
  return (mpfr_inf_p(op) && (mpfr_sgn(op) > 0));
}

int mpfr_is_negative_infinity(mpfr_t op) {
  return (mpfr_inf_p(op) && (mpfr_sgn(op) < 0));
}

int sollya_mpfi_has_positive_numbers(sollya_mpfi_t op) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  if ( sollya_mpfi_has_nan_opt(op) || sollya_mpfi_is_empty_opt(op) )
    return 0;
  else return (mpfr_sgn(&(op->right)) > 0);
}

int sollya_mpfi_has_negative_numbers(sollya_mpfi_t op) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  if ( sollya_mpfi_has_nan_opt(op) || sollya_mpfi_is_empty_opt(op) )
    return 0;
  else return (mpfr_sgn(&(op->left)) < 0);
}

int sollya_mpfi_is_nonneg(sollya_mpfi_t op) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  if ( sollya_mpfi_has_nan_opt(op) || sollya_mpfi_is_empty_opt(op) )
    return 0;
  else return (mpfr_sgn(&(op->left)) >= 0);
}

int sollya_mpfi_is_nonpos(sollya_mpfi_t op) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  if ( sollya_mpfi_has_nan_opt(op) || sollya_mpfi_is_empty_opt(op) )
    return 0;
  else return (mpfr_sgn(&(op->right)) <= 0);
}

int sollya_mpfi_has_positive_infinity(sollya_mpfi_t op) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  if ( sollya_mpfi_has_nan_opt(op) || sollya_mpfi_is_empty_opt(op) )
    return 0;
  else return mpfr_is_positive_infinity(&(op->right));
}

int sollya_mpfi_is_positive_infinity(sollya_mpfi_t op) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  if ( sollya_mpfi_has_nan_opt(op) || sollya_mpfi_is_empty_opt(op) )
    return 0;
  else return mpfr_is_positive_infinity(&(op->left));
}

int sollya_mpfi_has_negative_infinity(sollya_mpfi_t op) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  if ( sollya_mpfi_has_nan_opt(op) || sollya_mpfi_is_empty_opt(op) )
    return 0;
  else return mpfr_is_negative_infinity(&(op->left));
}

int sollya_mpfi_is_negative_infinity(sollya_mpfi_t op) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  if ( sollya_mpfi_has_nan_opt(op) || sollya_mpfi_is_empty_opt(op) )
    return 0;
  else return mpfr_is_negative_infinity(&(op->right));
}

int sollya_mpfi_has_infinity(sollya_mpfi_t op) {
  return (sollya_mpfi_has_negative_infinity(op) || sollya_mpfi_has_positive_infinity(op));
}

int sollya_mpfi_is_infinity(sollya_mpfi_t op) {
  return (sollya_mpfi_is_negative_infinity(op) || sollya_mpfi_is_positive_infinity(op));
}

int sollya_mpfi_inf_p(sollya_mpfi_t op) {
  return sollya_mpfi_has_infinity(op);
}



/* Functions that create sollya_mpfi_t */

int sollya_mpfi_set(sollya_mpfi_t rop, sollya_mpfi_srcptr op) {
  int res;
  res = mpfi_set(rop,op);
  sollya_mpfi_nan_normalize_opt(rop);
  sollya_mpfi_empty_normalize_opt(rop);
  return res;
}

int sollya_mpfi_set_z_2exp(sollya_mpfi_t rop, mpz_t op1, mp_exp_t op2) {
  int res, ra, rb;
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  ra = mpfr_set_z_2exp(&(rop->left), op1, op2, GMP_RNDD);
  rb = mpfr_set_z_2exp(&(rop->right), op1, op2, GMP_RNDU);
  if ((ra == 0) && (rb == 0)) {
    res = MPFI_FLAGS_BOTH_ENDPOINTS_EXACT;
  } else {
    if ((ra != 0) && (rb != 0)) {
	res = MPFI_FLAGS_BOTH_ENDPOINTS_EXACT;
    } else {
      if (ra != 0) {
	res = MPFI_FLAGS_LEFT_ENDPOINT_INEXACT;
      } else {
	res = MPFI_FLAGS_RIGHT_ENDPOINT_INEXACT;
      }
    }
  }
  sollya_mpfi_nan_normalize_opt(rop);
  sollya_mpfi_empty_normalize_opt(rop);
  return res;
}

int mpfi_to_sollya_mpfi(sollya_mpfi_t rop, mpfi_t op) {
  int res;
  res = mpfi_set(rop,op);
  sollya_mpfi_nan_normalize_opt(rop);
  sollya_mpfi_empty_normalize_opt(rop);
  return res;
}

int sollya_init_and_convert_interval(sollya_mpfi_t rop, mpfi_t op) {
  int res;

  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */

  sollya_mpfi_init2(rop, mpfi_get_prec(op));
  if ((!sollya_mpfi_has_nan_opt(op)) &&
      (mpfr_cmp(&(op->left), &(op->right)) > 0)) {
    printMessage(1,SOLLYA_MSG_RANGE_BOUNDS_IN_INVERSE_ORDER,"Warning: the bounds of a given interval are given in inverse order. Will revert them.\n");
    res = sollya_mpfi_interv_fr(rop, &(op->right), &(op->left));
  } else {
    if (sollya_mpfi_has_nan_opt(op)) {
      if ((!!(mpfr_nan_p(&(op->left)))) ^ (!!(mpfr_nan_p(&(op->right))))) {
	printMessage(1,SOLLYA_MSG_ONLY_ONE_ENDPOINT_OF_RANGE_IS_NAN,"Warning: one bound of a given interval is NaN while the other is not. Will normalize the interval to have two NaN endpoints.\n");
      }
      res = sollya_mpfi_set_nan_opt(rop);
    } else {
      res = mpfi_to_sollya_mpfi(rop, op);
    }
  }

  return res;
}

int sollya_mpfi_to_mpfi(mpfi_t rop, sollya_mpfi_t op) {
  return mpfi_set(rop,op);
}

int sollya_mpfi_set_d(sollya_mpfi_t rop, double op) {
  return mpfi_set_d(rop,op);
}

int sollya_mpfi_set_fr(sollya_mpfi_t rop, mpfr_t op) {
  return mpfi_set_fr(rop,op);
}

int sollya_mpfi_set_q(sollya_mpfi_t rop, mpq_t op) {
  int res;
  res = mpfi_set_q(rop,op);
  sollya_mpfi_nan_normalize_opt(rop);
  sollya_mpfi_empty_normalize_opt(rop);
  return res;
}

int sollya_mpfi_set_si(sollya_mpfi_t rop, long op) {
  return mpfi_set_si(rop,op);
}

int sollya_mpfi_set_ui(sollya_mpfi_t rop, unsigned long op) {
  return mpfi_set_ui(rop,op);
}

int sollya_mpfi_interv_d(sollya_mpfi_t rop, double op1, double op2) {
  int res;
  res = mpfi_interv_d(rop,op1,op2);
  sollya_mpfi_nan_normalize_opt(rop);
  sollya_mpfi_empty_normalize_opt(rop);
  return res;
}

int sollya_mpfi_interv_fr(sollya_mpfi_t rop, mpfr_t op1, mpfr_t op2) {
  int res;
  res = mpfi_interv_fr(rop,op1,op2);
  sollya_mpfi_nan_normalize_opt(rop);
  sollya_mpfi_empty_normalize_opt(rop);
  return res;
}

int sollya_mpfi_interv_si(sollya_mpfi_t rop, long op1, long op2) {
  int res;
  res = mpfi_interv_si(rop,op1,op2);
  sollya_mpfi_nan_normalize_opt(rop);
  sollya_mpfi_empty_normalize_opt(rop);
  return res;
}

int sollya_mpfi_interv_d_safe(sollya_mpfi_t rop, double op1, double op2) {
  int res = sollya_mpfi_interv_d(rop, op1, op2);
  if (sollya_mpfi_is_empty_opt(rop)) {
    sollyaFprintf(stderr,"Error: trying to define an interval with reversed bounds.\nThis should never happen. Please report the bug to the developers.\n");
    exit(1);
  }
  return res;
}

int sollya_mpfi_interv_fr_safe(sollya_mpfi_t rop, mpfr_t op1, mpfr_t op2) {
  int res = sollya_mpfi_interv_fr(rop, op1, op2);
  if (sollya_mpfi_is_empty_opt(rop)) {
    sollyaFprintf(stderr,"Error: trying to define an interval with reversed bounds.\nThis should never happen. Please report the bug to the developers.\n");
    exit(1);
  }
  return res;
}

int sollya_mpfi_interv_si_safe(sollya_mpfi_t rop, long op1, long op2) {
  int res = sollya_mpfi_interv_si(rop, op1, op2);
  if (sollya_mpfi_is_empty_opt(rop)) {
    sollyaFprintf(stderr,"Error: trying to define an interval with reversed bounds.\nThis should never happen. Please report the bug to the developers.\n");
    exit(1);
  }
  return res;
}

int sollya_mpfi_interv_si_2exp(sollya_mpfi_t rop, long op1, mp_exp_t op2, long op3, mp_exp_t op4) {
  int res, ra, rb;
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  ra = mpfr_set_si_2exp(&(rop->left), op1, op2, GMP_RNDD);
  rb = mpfr_set_si_2exp(&(rop->right), op3, op4, GMP_RNDU);
  if ((ra == 0) && (rb == 0)) {
    res = MPFI_FLAGS_BOTH_ENDPOINTS_EXACT;
  } else {
    if ((ra != 0) && (rb != 0)) {
	res = MPFI_FLAGS_BOTH_ENDPOINTS_EXACT;
    } else {
      if (ra != 0) {
	res = MPFI_FLAGS_LEFT_ENDPOINT_INEXACT;
      } else {
	res = MPFI_FLAGS_RIGHT_ENDPOINT_INEXACT;
      }
    }
  }
  sollya_mpfi_nan_normalize_opt(rop);
  sollya_mpfi_empty_normalize_opt(rop);
  return res;
}

/* Elementary univariate functions */

#define define_simple_func(f)                                           \
  int sollya_mpfi_##f (sollya_mpfi_t rop, sollya_mpfi_t op) {           \
    int res;                                                            \
                                                                        \
    if (sollya_mpfi_is_empty_opt(op)) return sollya_mpfi_set_empty_opt(rop);    \
                                                                        \
    res = mpfi_##f (rop,op); sollya_mpfi_nan_normalize_opt(rop);            \
    return res;                                                         \
  }

define_simple_func(abs)
define_simple_func(acos)
define_simple_func(acosh)
define_simple_func(asin)
define_simple_func(asinh)
define_simple_func(atan)
define_simple_func(atanh)
define_simple_func(cosh)
define_simple_func(sinh)
define_simple_func(tanh)
define_simple_func(exp)
define_simple_func(expm1)
define_simple_func(sqr)
define_simple_func(sqrt)


#define define_trig_func(f, compute_range)                              \
  int sollya_mpfi_##f (sollya_mpfi_t rop, sollya_mpfi_t op) {           \
    int res;                                                            \
                                                                        \
    if (sollya_mpfi_is_empty_opt(op)) return sollya_mpfi_set_empty_opt(rop);    \
    if (sollya_mpfi_has_infinity(op)) {                                 \
      compute_range;                                                    \
      res = MPFI_FLAGS_BOTH_ENDPOINTS_EXACT;                            \
    }                                                                   \
    else res = mpfi_##f (rop,op);                                       \
                                                                        \
    sollya_mpfi_nan_normalize_opt(rop);                                     \
    return res;                                                         \
  }                                                                     \

define_trig_func(cos, sollya_mpfi_interv_si(rop,-1,1))
define_trig_func(sin, sollya_mpfi_interv_si(rop,-1,1))
define_trig_func(tan, sollya_mpfi_set_full_range(rop))


/* HACK ALERT: For performance reasons, we will access the internals
   of an mpfi_t !!!
*/
#define define_log_func(f, left_bound_of_the_domain)                    \
  int sollya_mpfi_##f (sollya_mpfi_t rop, sollya_mpfi_t op) {           \
    int res;                                                            \
                                                                        \
    if (sollya_mpfi_is_empty_opt(op)) return sollya_mpfi_set_empty_opt(rop);    \
    if (sollya_mpfi_has_nan_opt(op)) return sollya_mpfi_set_nan_opt(rop);       \
    if (mpfr_cmp_si(&(op->left), left_bound_of_the_domain) < 0)         \
      return sollya_mpfi_set_nan_opt(rop);                                  \
    if (mpfr_cmp_si(&(op->left), left_bound_of_the_domain) == 0) {      \
      if (mpfr_cmp_si(&(op->right), left_bound_of_the_domain) == 0) {   \
        sollya_mpfi_set_negative_inf(rop);                              \
        return MPFI_FLAGS_BOTH_ENDPOINTS_EXACT;                         \
      }                                                                 \
      else {                                                            \
        mpfr_set_inf(&(rop->left),-1);                                  \
        if (mpfr_##f (&(rop->right),&(op->right),GMP_RNDU) == 0)        \
          res = MPFI_FLAGS_BOTH_ENDPOINTS_EXACT;                        \
        else                                                            \
          res = MPFI_FLAGS_RIGHT_ENDPOINT_INEXACT;                      \
	sollya_mpfi_zero_sign_normalize(rop);                           \
      }                                                                 \
    }                                                                   \
    else res = mpfi_##f (rop,op);                                       \
                                                                        \
    sollya_mpfi_nan_normalize_opt(rop);                                     \
                                                                        \
    return res;                                                         \
  }

define_log_func(log, 0)
define_log_func(log2, 0)
define_log_func(log10, 0)
define_log_func(log1p, -1)


/* Tricky binary functions */


/* Cases for the division:
   # op1 = NaN or op2 = NaN -> NaN
   # op1 empty or op2 empty -> empty
   # If op1 = [0] :
   # If op2 = [0] -> NaN
   # Else -> [0]    the argument here is: even if op2 contains 0,
   the function op1/y is 0 everywhere
   so we define it at y=0 by
   continuity.
   # If op1 = [+Inf] or op1 = [-Inf] :
   # If op2 = [+Inf] or op2 = [-Inf] -> NaN
   # If 0 is inside op2 or op2=[0] -> [-Inf, Inf]    The argument when op2=[0] is:
   Inf/[0] is an Inf but we do not know its sign. So
   we return [-Inf, Inf] and we are sure.
   # If 0 is not inside op2 (this includes the case when 0 does not belong to op2)
   -> sgn(op2)*op1     the argument is: even
   if op2 contains an Inf,
   the function op1/y is constant
   to Inf everywhere, so we
   define it at Inf as this constant.


   Now, we know that op1 is neither [0], nor [-Inf] or [+Inf]
   #  If op2 = [0] -> [-Inf, Inf]    the argument here is: even if op1 contains 0,
   x/[0] is always -Inf or +Inf, hence we
   define it at 0 by continuity as -Inf or +Inf.
   Sadly, we do not know the sign of 0, so we have
   to return [-Inf, +Inf] as a result
   # If op2 = [-Inf] or [+Inf] -> [0]     (even if op1 contains an Inf: by continuity)

   Now, we know that neither op1 nor op2 are singular point intervals
   # If op2 does not contain 0 -> mpfi_div(op1, op2)
   # If op2 has 0 inside -> [-Inf, +Inf]
   # Else op2 = [0, b] or op2 = [a, 0]:
   # In both cases, if op1 has 0 inside -> [-Inf, +Inf]

   # If op2 = [0, b]
   # If op1 = [u, v] with 0<=u ->  [u/b, +Inf]
   # If op1 = [u, v] with v<=0 ->  [-Inf, v/b]
   # Else (op2 = [a, 0]) :
   # If op1 = [u, v] with 0<=u ->  [-Inf, u/a]
   # If op1 = [u, v] with v<=0 ->  [v/a, +Inf]
*/
int sollya_mpfi_div(sollya_mpfi_t rop, sollya_mpfi_t op1, sollya_mpfi_t op2) {
  int res;
  int sign;

  if (sollya_mpfi_is_empty_opt(op1) || sollya_mpfi_is_empty_opt(op2))
    return sollya_mpfi_set_empty_opt(rop);
  if (sollya_mpfi_has_nan_opt(op1) || sollya_mpfi_has_nan_opt(op2))
    return sollya_mpfi_set_nan_opt(rop);

  if (sollya_mpfi_is_zero(op1)) {
    if (sollya_mpfi_is_zero(op2)) return sollya_mpfi_set_nan_opt(rop);
    else return sollya_mpfi_set_si(rop, 0);
  }

  if (sollya_mpfi_is_infinity(op1)) {
    if (sollya_mpfi_is_infinity(op2)) return sollya_mpfi_set_nan_opt(rop);
    if (sollya_mpfi_is_zero(op2) || sollya_mpfi_has_zero_inside(op2))
      return sollya_mpfi_set_full_range(rop);

    /* Now the case op1=[+/-Inf], and op2 has constant sign */
    sign = (sollya_mpfi_has_positive_numbers(op2) ? +1 : -1);
    if (sollya_mpfi_is_negative_infinity(op1)) sign = -sign;
    if (sign > 0) return sollya_mpfi_set_positive_inf(rop);
    else return sollya_mpfi_set_negative_inf(rop);
  }

  /* Here op1 != [0], op1 != [+/-Inf] */
  if (sollya_mpfi_is_zero(op2)) return sollya_mpfi_set_full_range(rop);
  if (sollya_mpfi_is_infinity(op2)) return sollya_mpfi_set_si(rop, 0);

  /* Here op1 != [0], op1 != [+/-Inf], op2 != [0], op2 != [+/-Inf] */
  if (!sollya_mpfi_has_zero(op2)) return mpfi_div(rop, op1, op2);
  if (sollya_mpfi_has_zero_inside(op2)) return sollya_mpfi_set_full_range(rop);

  if (sollya_mpfi_has_zero_inside(op1)) return sollya_mpfi_set_full_range(rop);


  /* Now one of the bounds of op2 (and only one) is zero. */
  /* If op1 contains 0, it is one of its bounds (and only one) */

  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  if (sollya_mpfi_has_positive_numbers(op2)) {
    if (mpfr_sgn(&(op1->left)) >= 0) { /* Case 0<=u and op2=[0,b] */
      mpfr_set_inf(&(rop->right),1);
      if (mpfr_div(&(rop->left),&(op1->left),&(op2->right),GMP_RNDD) == 0)
        res = MPFI_FLAGS_BOTH_ENDPOINTS_EXACT;
      else
        res = MPFI_FLAGS_LEFT_ENDPOINT_INEXACT;
      sollya_mpfi_zero_sign_normalize(rop);
    }
    else { /* Case v<=0 and op2=[0,b] */
      mpfr_set_inf(&(rop->left),-1);
      if (mpfr_div(&(rop->right),&(op1->right),&(op2->right),GMP_RNDU) == 0)
        res = MPFI_FLAGS_BOTH_ENDPOINTS_EXACT;
      else
        res = MPFI_FLAGS_RIGHT_ENDPOINT_INEXACT;
      sollya_mpfi_zero_sign_normalize(rop);
    }
  }
  else { /* Case op2=[a,0]... */
    if (mpfr_sgn(&(op1->left)) >= 0) { /* ...and 0<=u */
      mpfr_set_inf(&(rop->left),-1);
      if (mpfr_div(&(rop->right),&(op1->left),&(op2->left),GMP_RNDU) == 0)
        res = MPFI_FLAGS_BOTH_ENDPOINTS_EXACT;
      else
        res = MPFI_FLAGS_RIGHT_ENDPOINT_INEXACT;
      sollya_mpfi_zero_sign_normalize(rop);
    }
    else { /* case v<=0 */
      mpfr_set_inf(&(rop->right),1);
      if (mpfr_div(&(rop->left),&(op1->right),&(op2->left),GMP_RNDD) == 0)
        res = MPFI_FLAGS_BOTH_ENDPOINTS_EXACT;
      else
        res = MPFI_FLAGS_LEFT_ENDPOINT_INEXACT;
      sollya_mpfi_zero_sign_normalize(rop);
    }
  }

  sollya_mpfi_nan_normalize_opt(rop);
  return res;
}


int sollya_mpfi_mul(sollya_mpfi_t rop, sollya_mpfi_t op1, sollya_mpfi_t op2) {
  /*
    # op1 = NaN or op2 = NaN -> NaN
    # op1 empty or op2 empty -> empty
    # If op1 = [+Inf] or op1 = [-Inf] :
    # If op2 = [0] -> NaN
    # If op2 has 0 inside -> [-Inf, Inf]    The argument here is: Inf*y is an +/-Inf for
    every y. Since y can have both signs, we
    can have bot +Inf and -Inf.
    # Else -> sgn(op2)*op1
    # If op2 = [-Inf] or [+Inf] :
    # If op1 = [0] -> NaN
    # If op1 has 0 inside -> [-Inf, +Inf]
    # Else -> sgn(op1)*op2

    Now, we know that neither op1 or op2 is [-Inf] or [Inf].
    # If op1 = [0] or op2 = [0] -> [0] the argument here is: even if the other contains an Inf,
    the function 0*y is 0 everywhere
    so we define it at y=Inf by
    continuity.

    Now, we know that neither op1 nor op2 are singular point intervals
    # Else -> mpfi_mul(op1, op2)
  */

  if (sollya_mpfi_is_empty_opt(op1) || sollya_mpfi_is_empty_opt(op2)) return sollya_mpfi_set_empty_opt(rop);
  if (sollya_mpfi_has_nan_opt(op1) || sollya_mpfi_has_nan_opt(op2))  return sollya_mpfi_set_nan_opt(rop);

  if (sollya_mpfi_is_infinity(op1)) {
    if (sollya_mpfi_is_zero(op2)) return sollya_mpfi_set_nan_opt(rop);
    if (sollya_mpfi_has_zero_inside(op2)) return sollya_mpfi_set_full_range(rop);
    if (sollya_mpfi_is_nonneg(op2)) return sollya_mpfi_set(rop, op1);
    /* else: op2<=0 */
    return sollya_mpfi_neg(rop, op1);
  }
  if (sollya_mpfi_is_infinity(op2)) {
    if (sollya_mpfi_is_zero(op1)) return sollya_mpfi_set_nan_opt(rop);
    if (sollya_mpfi_has_zero_inside(op1)) return sollya_mpfi_set_full_range(rop);
    if (sollya_mpfi_is_nonneg(op1)) return sollya_mpfi_set(rop, op2);
    /* else: op1<=0 */
    return sollya_mpfi_neg(rop, op2);
  }
  if (sollya_mpfi_is_zero(op1) || sollya_mpfi_is_zero(op2)) return sollya_mpfi_set_ui(rop, 0);

  /* else: default case: */
  return mpfi_mul(rop,op1,op2);
}



/* Easy binary operations */

int sollya_mpfi_add(sollya_mpfi_t rop, sollya_mpfi_t op1, sollya_mpfi_t op2) {
  int res;

  if (sollya_mpfi_is_empty_opt(op1)) return sollya_mpfi_set_empty_opt(rop);
  if (sollya_mpfi_is_empty_opt(op2)) return sollya_mpfi_set_empty_opt(rop);

  res = mpfi_add(rop,op1,op2); sollya_mpfi_nan_normalize_opt(rop);

  return res;
}

int sollya_mpfi_add_fr(sollya_mpfi_t rop, sollya_mpfi_t op1, mpfr_t op2) {
  int res;

  if (sollya_mpfi_is_empty_opt(op1)) return sollya_mpfi_set_empty_opt(rop);

  res = mpfi_add_fr(rop,op1,op2); sollya_mpfi_nan_normalize_opt(rop);

  return res;
}

int sollya_mpfi_add_ui(sollya_mpfi_t rop, sollya_mpfi_t op1, unsigned long op2) {
  int res;

  if (sollya_mpfi_is_empty_opt(op1)) return sollya_mpfi_set_empty_opt(rop);

  res = mpfi_add_ui(rop,op1,op2); sollya_mpfi_nan_normalize_opt(rop);

  return res;
}

int sollya_mpfi_sub(sollya_mpfi_t rop, sollya_mpfi_t op1, sollya_mpfi_t op2) {
  int res;

  if (sollya_mpfi_is_empty_opt(op1)) return sollya_mpfi_set_empty_opt(rop);
  if (sollya_mpfi_is_empty_opt(op2)) return sollya_mpfi_set_empty_opt(rop);

  res = mpfi_sub(rop,op1,op2); sollya_mpfi_nan_normalize_opt(rop);

  return res;
}

int sollya_mpfi_sub_fr(sollya_mpfi_t rop, sollya_mpfi_t op1, mpfr_t op2) {
  int res;
  mpfi_t tmp;

  if (sollya_mpfi_is_empty_opt(op1)) return sollya_mpfi_set_empty_opt(rop);

  mpfi_init2(tmp, mpfr_get_prec(op2));
  mpfi_set_fr(tmp, op2); /* exact */
  res = mpfi_sub(rop,op1,tmp);
  mpfi_clear(tmp);
  sollya_mpfi_nan_normalize_opt(rop);

  return res;
}

int sollya_mpfi_sub_ui(sollya_mpfi_t rop, sollya_mpfi_t op1, unsigned long op2) {
  int res;

  if (sollya_mpfi_is_empty_opt(op1)) return sollya_mpfi_set_empty_opt(rop);

  res = mpfi_sub_ui(rop,op1,op2); sollya_mpfi_nan_normalize_opt(rop);

  return res;
}

int sollya_mpfi_div_ui(sollya_mpfi_t rop, sollya_mpfi_t op1, unsigned long op2) {
  int res;
  mpfi_t temp;

  if (sollya_mpfi_is_empty_opt(op1)) return sollya_mpfi_set_empty_opt(rop);

  mpfi_init2(temp,8 * sizeof(op2));  mpfi_set_ui(temp,op2);
  res = sollya_mpfi_div(rop,op1,temp); sollya_mpfi_nan_normalize_opt(rop);

  mpfi_clear(temp);
  return res;
}

int sollya_mpfi_div_z(sollya_mpfi_t rop, sollya_mpfi_t op1, mpz_t op2) {
  int res, ra, rb;

  /* Empty input */
  if (sollya_mpfi_is_empty_opt(op1)) return sollya_mpfi_set_empty_opt(rop);

  /* Division by 0 */
  if (mpz_sgn(op2) == 0) return sollya_mpfi_div_ui(rop, op1, 0);

  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  
  /* General case */
  if (mpz_sgn(op2) > 0) {
    ra = mpfr_div_z(&(rop->left),&(op1->left),op2,GMP_RNDD);  
    rb = mpfr_div_z(&(rop->right),&(op1->right),op2,GMP_RNDU);  
  } else {
    ra = mpfr_div_z(&(rop->left),&(op1->right),op2,GMP_RNDD);  
    rb = mpfr_div_z(&(rop->right),&(op1->left),op2,GMP_RNDU);  
  }
  if ((ra == 0) && (rb == 0)) {
    res = MPFI_FLAGS_BOTH_ENDPOINTS_EXACT;
  } else {
    if ((ra != 0) && (rb != 0)) {
	res = MPFI_FLAGS_BOTH_ENDPOINTS_EXACT;
    } else {
      if (ra != 0) {
	res = MPFI_FLAGS_LEFT_ENDPOINT_INEXACT;
      } else {
	res = MPFI_FLAGS_RIGHT_ENDPOINT_INEXACT;
      }
    }
  }
  sollya_mpfi_nan_normalize_opt(rop);
  sollya_mpfi_empty_normalize_opt(rop);

  return res;
}


int sollya_mpfi_mul_ui(sollya_mpfi_t rop, sollya_mpfi_t op1, unsigned long op2) {
  int res;
  mpfi_t temp;

  if (sollya_mpfi_is_empty_opt(op1)) return sollya_mpfi_set_empty_opt(rop);

  mpfi_init2(temp,8 * sizeof(op2));  mpfi_set_ui(temp,op2);
  res = sollya_mpfi_mul(rop,op1,temp); sollya_mpfi_nan_normalize_opt(rop);

  mpfi_clear(temp);
  return res;
}

int sollya_mpfi_ui_div(sollya_mpfi_t rop, unsigned long op1, sollya_mpfi_t op2) {
  int res;
  mpfi_t temp;

  if (sollya_mpfi_is_empty_opt(op2)) return sollya_mpfi_set_empty_opt(rop);

  mpfi_init2(temp,8 * sizeof(op1));  mpfi_set_ui(temp,op1);
  res = sollya_mpfi_div(rop,temp, op2); sollya_mpfi_nan_normalize_opt(rop);

  mpfi_clear(temp);
  return res;
}



/* Other functions */

int sollya_mpfi_inv(sollya_mpfi_t rop, sollya_mpfi_t op) {
  return sollya_mpfi_ui_div(rop, 1, op);
}

int sollya_mpfi_blow(sollya_mpfi_t rop, sollya_mpfi_t op1, double op2) {
  int res;

  if (sollya_mpfi_is_empty_opt(op1)) return sollya_mpfi_set_empty_opt(rop);

  res = mpfi_blow(rop,op1,op2); sollya_mpfi_nan_normalize_opt(rop);
  return res;
}

void sollya_mpfi_blow_1ulp(sollya_mpfi_t rop) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  mpfr_nextbelow(&(rop->left));  
  mpfr_nextabove(&(rop->right));
}

int sollya_mpfi_bounded_p(sollya_mpfi_t op) {
  if (sollya_mpfi_has_nan_opt(op)) return 0;
  if (sollya_mpfi_is_empty_opt(op)) return 1;
  return mpfi_bounded_p(op);
}

void sollya_mpfi_clear(sollya_mpfi_t op) {
  mpfi_clear(op);
}

int sollya_mpfi_const_pi(sollya_mpfi_t rop) {
  return mpfi_const_pi(rop);
}

int sollya_mpfi_diam_abs(mpfr_t rop, sollya_mpfi_t op) {
  if (sollya_mpfi_has_nan_opt(op) || sollya_mpfi_is_empty_opt(op)) {
    mpfr_set_nan(rop);
    return 0;
  }
  if (sollya_mpfi_is_infinity(op)) return mpfr_set_ui(rop, 0, GMP_RNDN);

  /* else... */
  return mpfi_diam_abs(rop, op);
}

int sollya_mpfi_get_left(mpfr_t rop, sollya_mpfi_t op) {
  if (sollya_mpfi_has_nan_opt(op) || sollya_mpfi_is_empty_opt(op)) {
    mpfr_set_nan(rop);
    return 0;
  }
  return mpfi_get_left(rop,op);
}

int sollya_mpfi_get_right(mpfr_t rop, sollya_mpfi_t op) {
  if (sollya_mpfi_has_nan_opt(op) || sollya_mpfi_is_empty_opt(op)) {
    mpfr_set_nan(rop);
    return 0;
  }
  return mpfi_get_right(rop,op);
}

void sollya_mpfi_get_fr(mpfr_t rop, sollya_mpfi_t op) {
  if (sollya_mpfi_has_nan_opt(op) || sollya_mpfi_is_empty_opt(op))
    mpfr_set_nan(rop);
  else mpfi_get_fr(rop,op);
}

mp_prec_t sollya_mpfi_get_prec(sollya_mpfi_srcptr op) {
  return mpfi_get_prec(op);
}

const char *sollya_mpfi_get_version() {
  const char *res;
  res = (const char *) mpfi_get_version();
  return res;
}

void sollya_mpfi_init2(sollya_mpfi_t rop, mp_prec_t op) {
  mpfi_init2(rop,op);
}

int sollya_mpfi_intersect(sollya_mpfi_t rop, sollya_mpfi_t op1, sollya_mpfi_t op2) {
  int res;

  if (sollya_mpfi_is_empty_opt(op1) || sollya_mpfi_is_empty_opt(op2))
    return sollya_mpfi_set_empty_opt(rop);

  if (sollya_mpfi_has_nan_opt(op1) || sollya_mpfi_has_nan_opt(op2))
    return sollya_mpfi_set_nan_opt(rop);

  res = mpfi_intersect(rop,op1,op2);
  sollya_mpfi_empty_normalize_opt(rop);
  sollya_mpfi_nan_normalize_opt(rop);

  return res;
}

int sollya_mpfi_is_inside(sollya_mpfi_t op1, sollya_mpfi_t op2) {
  if (sollya_mpfi_is_empty_opt(op1)) return 0;
  if (sollya_mpfi_is_empty_opt(op2)) return 1;

  if (sollya_mpfi_has_nan_opt(op1) || sollya_mpfi_has_nan_opt(op2)) return 0;
  else return mpfi_is_inside(op1,op2);
}

int sollya_mpfi_mid(mpfr_t rop, sollya_mpfi_t op) {
  if (sollya_mpfi_has_nan_opt(op) || sollya_mpfi_is_empty_opt(op)) {
    mpfr_set_nan(rop);
    return 0;
  }
  return mpfi_mid(rop,op);
}

int sollya_mpfi_neg(sollya_mpfi_t rop, sollya_mpfi_t op) {
  int res;

  if (sollya_mpfi_is_empty_opt(op)) return sollya_mpfi_set_empty_opt(rop);

  res = mpfi_neg(rop,op);  sollya_mpfi_nan_normalize_opt(rop);
  return res;
}

void sollya_mpfi_set_prec(sollya_mpfi_t rop, mp_prec_t op) {
  mpfi_set_prec(rop,op);
}

int sollya_mpfi_prec_round(sollya_mpfi_t rop, mp_prec_t op) {
  int res, ra, rb;

  if (sollya_mpfi_is_empty_opt(rop)) {
    sollya_mpfi_set_prec(rop, op);
    return sollya_mpfi_set_empty_opt(rop);
  }
  
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  ra = mpfr_prec_round(&(rop->left), op, GMP_RNDD);
  rb = mpfr_prec_round(&(rop->right), op, GMP_RNDU);
  if ((ra == 0) && (rb == 0)) {
    res = MPFI_FLAGS_BOTH_ENDPOINTS_EXACT;
  } else {
    if ((ra != 0) && (rb != 0)) {
	res = MPFI_FLAGS_BOTH_ENDPOINTS_EXACT;
    } else {
      if (ra != 0) {
	res = MPFI_FLAGS_LEFT_ENDPOINT_INEXACT;
      } else {
	res = MPFI_FLAGS_RIGHT_ENDPOINT_INEXACT;
      }
    }
  }
  sollya_mpfi_nan_normalize_opt(rop);
  sollya_mpfi_empty_normalize_opt(rop);

  return res;
}


int sollya_mpfi_union(sollya_mpfi_t rop, sollya_mpfi_t op1, sollya_mpfi_t op2) {
  int res;

  if (sollya_mpfi_has_nan_opt(op1) || sollya_mpfi_has_nan_opt(op2)) return sollya_mpfi_set_nan_opt(rop);

  if (sollya_mpfi_is_empty_opt(op1)) res = sollya_mpfi_set(rop, op2);
  else {
    if (sollya_mpfi_is_empty_opt(op2)) res = sollya_mpfi_set(rop, op1);
    else res = mpfi_union(rop, op1, op2);
  }

  sollya_mpfi_nan_normalize_opt(rop);
  sollya_mpfi_empty_normalize_opt(rop);
  return res;
}

int sollya_mpfi_is_point_and_real(sollya_mpfi_t op) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  return ((!(mpfr_nan_p(&(op->left)) || mpfr_nan_p(&(op->right)))) &&
	  (!(mpfr_inf_p(&(op->left)) || mpfr_inf_p(&(op->right)))) &&
	  (mpfr_equal_p(&(op->left),&(op->right))));
}

int sollya_mpfi_is_quasi_point_and_real(sollya_mpfi_t op) {
  mp_exp_t el, er, ea, eb, d;

  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  if (!mpfr_number_p(&(op->left))) return 0;
  if (!mpfr_number_p(&(op->right))) return 0;
  if (mpfr_equal_p(&(op->left), &(op->right))) return 1;
  if (mpfr_get_prec(&(op->left)) != mpfr_get_prec(&(op->right))) return 0;
  if (mpfr_cmp(&(op->left), &(op->right)) > 0) return 0;
  if (mpfr_zero_p(&(op->left)) || mpfr_zero_p(&(op->right))) return 0;
  if (mpfr_sgn(&(op->left)) != mpfr_sgn(&(op->right))) return 0;
  el = mpfr_get_exp(&(op->left));
  er = mpfr_get_exp(&(op->right));
  ea = el; eb = er;
  if (eb > ea) {
    ea = er; eb = el;
  }
  d = ea - eb;
  if ((d < 0) || (d > 1)) return 0;
  mpfr_nextabove(&(op->left));
  mpfr_nextabove(&(op->left));
  if (mpfr_cmp(&(op->left), &(op->right)) >= 0) {
      mpfr_nextbelow(&(op->left));
      mpfr_nextbelow(&(op->left));
      return 1;
  }
  mpfr_nextbelow(&(op->left));
  mpfr_nextbelow(&(op->left));
  return 0;
}

int sollya_mpfi_equal_p(sollya_mpfi_t op1, sollya_mpfi_t op2) {
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  return (mpfr_equal_p(&(op1->left),&(op2->left)) && 
	  mpfr_equal_p(&(op1->right),&(op2->right)));
}

mp_exp_t sollya_mpfi_max_exp(sollya_mpfi_t op) {
  mp_exp_t rl, rr;

 /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  if (!mpfr_number_p(&(op->left))) return mpfr_get_emin_min();
  if (!mpfr_number_p(&(op->right))) return mpfr_get_emin_min();
  if (mpfr_zero_p(&(op->left))) {
    if (mpfr_zero_p(&(op->right))) {
      return mpfr_get_emin_min();
    } else {
      return mpfr_get_exp(&(op->right));
    }
  } else {
    if (mpfr_zero_p(&(op->right))) {
      return mpfr_get_exp(&(op->left));
    } else {
      rl = mpfr_get_exp(&(op->left));
      rr = mpfr_get_exp(&(op->right));
      return (rl > rr ? rl : rr);
    }
  }
  return mpfr_get_emin_min();
}

int sollya_mpfi_fr_in_interval(mpfr_t op1, sollya_mpfi_t op2) {
  if (!mpfr_number_p(op1)) return 0;
  if (sollya_mpfi_has_nan_opt(op2)) return 0;
  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  return (mpfr_lessequal_p(&(op2->left), op1) && 
	  mpfr_lessequal_p(op1, &(op2->right)));
}

int sollya_mpfi_enclosure_accurate_enough(sollya_mpfi_t op1, mp_prec_t op2) {
  mp_exp_t eA, eB;
  mp_prec_t p;
  mpfr_t diam, temp;
  int res;

  if (sollya_mpfi_has_nan_opt(op1)) return 0;
  if (sollya_mpfi_is_empty_opt(op1)) return 0;
  if (sollya_mpfi_has_infinity(op1)) return 0;
  if (sollya_mpfi_has_zero(op1)) return 0;
  if (op2 <= 2) return 0;

  /* HACK ALERT: For performance reasons, we will access the internals
     of an mpfi_t !!!
  */
  eA = mpfr_get_exp(&(op1->left));
  eB = mpfr_get_exp(&(op1->right));
  
  /* If there is a whole binade between the endpoints a and b, there
     is no accuracy 
  */
  if ((eB - eA) > 1) return 0;

  /* Compute the difference between the endpoints and compute 2^(-op2)
     * a, where a is the left endpoint 
  */
  p = mpfr_get_prec(&(op1->left));
  if (mpfr_get_prec(&(op1->right)) > p) p = mpfr_get_prec(&(op1->right));
  mpfr_init2(diam, p + 2);
  mpfr_init2(temp, p);
  mpfr_sub(diam, &(op1->right), &(op1->left), GMP_RNDN); /* exact, Sterbenz */
  mpfr_mul_2si(temp, &(op1->left), -op2, GMP_RNDN); /* exact, same precision */
  
  /* We got enough precision iff |b - a| <= 2^(-prec) * |a|, where a
     is the left endpoint, b the right one and prec = op2 
  */
  res = (mpfr_cmpabs(diam, temp) <= 0);
  mpfr_clear(temp);
  mpfr_clear(diam);

  return res;
}
