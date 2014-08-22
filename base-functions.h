/*

  Copyright 2014 by

  Centre de recherche INRIA Sophia-Antipolis Mediterranee, equipe APICS,
  Sophia Antipolis, France.

  Contributors S. Chevillard

  sylvain.chevillard@ens-lyon.org

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

#ifndef BASE_FUNCTIONS_H
#define BASE_FUNCTIONS_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gmp.h>
#include <mpfr.h>
#include "mpfi-compat.h"

#define SQRT 0
#define EXP 1
#define LOG 2
#define LOG_2 3
#define LOG_10 4
#define SIN 5
#define COS 6
#define TAN 7
#define ASIN 8
#define ACOS 9
#define ATAN 10
#define SINH 11
#define COSH 12
#define TANH 13
#define ASINH 14
#define ACOSH 15
#define ATANH 16
#define ABS 17
#define DOUBLE 18
#define DOUBLEDOUBLE 19
#define TRIPLEDOUBLE 20
#define ERF 21
#define ERFC 22
#define LOG_1P 23
#define EXP_M1 24
#define DOUBLEEXTENDED 25
#define CEIL 26
#define FLOOR 27
#define SINGLE 28
#define NEARESTINT 39
#define HALFPRECISION 30
#define QUAD 31

#define DECREASING 0         /* Indicates that for any x<y in the domain of f, f(x)>f(y) */
#define NONINCREASING 1      /* Indicates that for any x<y in the domain of f, f(x)>=f(y) */
#define MONOTONICITY_NONE 2  /* Indicates that f has no particular monotonic behavior */
#define NONDECREASING 3      /* Indicates that for any x<y in the domain of f, f(x)<=f(y) */
#define INCREASING 4         /* Indicates that for any x<y in the domain of f, f(x)<f(y) */

struct nodeStruct;
#ifndef NODE_TYPEDEF
#define NODE_TYPEDEF
typedef struct nodeStruct node;
#endif

typedef struct baseFunctionStruct baseFunction;
struct baseFunctionStruct
{
  int baseFunctionCode; /* The unique code defined above */
  char *functionName;   /* The name of the function as it is used inside Sollya */
  char *xmlString;      /* The xml code for the function */
  char *mpfrName;       /* The name of the function in mpfr. Used to generate code, e.g. as in implementconstant */
  int handledByImplementconst; /* A boolean. Functions that must not be handled by implementconstant are those functions having a discontinuity at a representable point */
  int isDefinedEverywhere; /* A boolean. True if the function is defined and finite on the whole real line */
  int onlyZeroIsZero; /* A boolean. True if the only real zero of the function is zero */
  int doesNotVanish; /* A boolean. True if the function never takes the value zero on the real line */
  int monotonicity; /* One of the monotonicity code, indicating the behavior of the function on its domain */
  void (*baseAutodiff)(sollya_mpfi_t *, sollya_mpfi_t, int, int *); /* Computes the vector of the f^(k)(x0)/k!, k=0..n. The last parameter is a silent parameter */
  int (*interval_eval)(sollya_mpfi_t, sollya_mpfi_t); /* Performs an interval evaluation ``à la'' mpfi */
  int (*point_eval)(mpfr_t, mpfr_t, mp_rnd_t); /* Performs an evaluation ``à la'' mpfr */
  node *(*diff_expr)(node *); /* If g if the argument, returns a tree representing diff(f o g) */
  node *(*simplify)(node *); /* If g is the argument (supposed already simplified as much as possible), returns a tree representing a simplification
                                (without introducing errors) of f(g).
                                Notice that g must either be used (eaten up) to construct the new tree, or be freed by the function  */
  int (*evalsign)(int *, node *); /* If s and g are the arguments, tries to determine the sign of f(g) assuming that g is a constant expression. In case of success, the sign is assigned to variable s and 1 is returned. Otherwise, s is left unchanged and 0 is returned */
};


extern static baseFunction *basefun_sqrt;
extern static baseFunction *basefun_exp;
extern static baseFunction *basefun_log;
extern static baseFunction *basefun_log2;
extern static baseFunction *basefun_log10;
extern static baseFunction *basefun_sin;
extern static baseFunction *basefun_cos;
extern static baseFunction *basefun_tan;
extern static baseFunction *basefun_asin;
extern static baseFunction *basefun_acos;
extern static baseFunction *basefun_atan;
extern static baseFunction *basefun_sinh;
extern static baseFunction *basefun_cosh;
extern static baseFunction *basefun_tanh;
extern static baseFunction *basefun_asinh;
extern static baseFunction *basefun_acosh;
extern static baseFunction *basefun_atanh;
extern static baseFunction *basefun_abs;
extern static baseFunction *basefun_double;
extern static baseFunction *basefun_single;
extern static baseFunction *basefun_halfprecision;
extern static baseFunction *basefun_quad;
extern static baseFunction *basefun_doubledouble;
extern static baseFunction *basefun_tripledouble;
extern static baseFunction *basefun_erf;
extern static baseFunction *basefun_erfc;
extern static baseFunction *basefun_log1p;
extern static baseFunction *basefun_expm1;
extern static baseFunction *basefun_doubleextended;
extern static baseFunction *basefun_ceil;
extern static baseFunction *basefun_floor;
extern static baseFunction *basefun_nearestint;

node *makeSqrt(node *op1);
node *makeExp(node *op1);
node *makeLog(node *op1);
node *makeLog2(node *op1);
node *makeLog10(node *op1);
node *makeSin(node *op1);
node *makeCos(node *op1);
node *makeTan(node *op1);
node *makeAsin(node *op1);
node *makeAcos(node *op1);
node *makeAtan(node *op1);
node *makeAbs(node *op1);
node *makeDouble(node *op1);
node *makeSingle(node *op1);
node *makeQuad(node *op1);
node *makeHalfPrecision(node *op1);
node *makeDoubledouble(node *op1);
node *makeTripledouble(node *op1);
node *makeErf(node *op1);
node *makeErfc(node *op1);
node *makeLog1p(node *op1);
node *makeExpm1(node *op1);
node *makeDoubleextended(node *op1);
node *makeCeil(node *op1);
node *makeFloor(node *op1);
node *makeNearestInt(node *op1);
node *makeSinh(node *op1);
node *makeCosh(node *op1);
node *makeTanh(node *op1);
node *makeAsinh(node *op1);
node *makeAcosh(node *op1);
node *makeAtanh(node *op1);

#endif /* ifdef BASE_FUNCTIONS */
