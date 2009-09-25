#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "sollya.h"
#define coeff(i,j,n) ((i)-1)*(n)+(j)-1
extern int mpfi_pow(mpfi_t res, mpfi_t x, mpfi_t y);
extern void fprintInterval(FILE *fd, mpfi_t interval);


void perturbPoints(mpfr_t *x, int p);
void system_solve(mpfr_t *res, mpfr_t *M, mpfr_t *b, int n, mp_prec_t prec);
void printMatrix(mpfr_t *M, int n);
node *constructPoly(mpfr_t *coeff, int n);
node *constructPolyFromRoots(mpfr_t *coeff, int n);
void findNewSetRoots(mpfr_t * newRoots, mpfr_t *roots, int nrRoots, mpfr_t *chebArray, int n);
void evaluateMpfiFunction(mpfi_t y, node *f, mpfi_t x, mp_prec_t prec);
int getIntervalAroundRoot(mpfi_t *res, node *dif, mpfr_t r, mp_prec_t prec);
/* This function performs the interpolation.
   See the commentaries below.
*/
node* interpolation(mpfr_t *coeffArray,  mpfi_t bound, node *f, mpfi_t x, int n);


#endif
