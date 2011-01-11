#include "sollya.h"

extern int polynomialDivide_mpfi(sollya_mpfi_t *quotient, int *quotient_degree, sollya_mpfi_t *rest, int *rest_degree, sollya_mpfi_t *p, int p_degree, sollya_mpfi_t *q, int q_degree, mp_prec_t prec);
/* Evaluates a symbolic polynomial at point x by Horner scheme */
extern void symbolic_poly_evaluation_horner(sollya_mpfi_t res, sollya_mpfi_t *coeffs_array, sollya_mpfi_t x, int degree);
extern void symbolic_poly_diff(sollya_mpfi_t *res, sollya_mpfi_t *coeff_array, int degree);

int mpfi_set_node( sollya_mpfi_t *r, node * c);

/*returns n chebyshev points in x*/
void getChebyshevPoints(sollya_mpfi_t *chebPoints, int n, sollya_mpfi_t x);

/*Returns n-1 chebyshev extrema in x*/ 
void getChebyshevExtrema(sollya_mpfi_t *chebPoints, int n, sollya_mpfi_t x);

/*The so called chebMatrix, contains 
the specific values of T_i(x_j) i=0..n-1, j=0..n-1
x_j \in [-1,1]. 
Note: It has to be computed only once for each n*/
void getChebMatrix(sollya_mpfi_t**chebMatrix, int n);

/*Returns the coeffs of the cheb polynomials, n*n*/
void getChebPolyCoeffs(mpz_t* chebMatrix, int n);

/*returns in chebCoeffs the coeffs of the finite Chebyshev basis expansion of the polynomial given by the monomial coefficients p, of degree n-1*/
void getPolyCoeffsChebBasis(sollya_mpfi_t *chebCoeffs, sollya_mpfi_t *p, int n);

/*returns in c the coeffs in the monomial basis for the polynomial p(a*x+b), where a and b are mpfi_s
 the polynomial p  of degree n-1 given by the monomial coefficients stored as mpfi_s in p*/
void getTranslatedPolyCoeffs(sollya_mpfi_t *c, sollya_mpfi_t *p, int n, sollya_mpfi_t a, sollya_mpfi_t b);

/*evaluate a function in n chebPoints*/
void getFunctionValues(sollya_mpfi_t* fValues, sollya_mpfi_t * chebPoints,node* f,int n);


/*this function computes the cheb coeffs for the interpolation polynomial of degree n-1 at chebpoints (we have n chebpoints and a n*n chebMatrix)*/
void getChebCoeffs(sollya_mpfi_t* coeffs, sollya_mpfi_t *chebMatrix, sollya_mpfi_t *fValues,int n);


/*wrapper to get directly the coeffs from the function*/
void getChebCoeffsFromFunction(sollya_mpfi_t* coeffs, sollya_mpfi_t * chebPoints, sollya_mpfi_t * chebMatrix,node *f,int n);


/*wrapper to get directly the coeffs in the chebyshev basis from a polynomial in the monomial basis(given a pointer to node, over a given interval x*/
void getChebCoeffsFromPolynomial(sollya_mpfi_t**coeffs, int*n, node *f, sollya_mpfi_t x);

/*wrapper to get directly the coeffs in the chebyshev basis up to degree n-1 (first n coeffs) and a bound for the remaining polynomial, from a polynomial in the monomial basis(given a pointer to node, over a given interval x*/
void getNChebCoeffsFromPolynomial(sollya_mpfi_t *coeffs, sollya_mpfi_t bound, node *f, sollya_mpfi_t x, int n);

/********************************Functions related to derivation and integration of polynomials in cheb basis*************************/
void getChebCoeffsDerivativePolynomial(sollya_mpfi_t*coeffs, sollya_mpfi_t *chebCoeffs, int n, sollya_mpfi_t x);

/*Computes the antiderivative of a polynomial in Chebyshev basis.
NOTE: the constant coefficient is set to zero, but it should be viewed as a constant*/
void getChebCoeffsIntegrationPolynomial(sollya_mpfi_t*coeffs, sollya_mpfi_t *chebCoeffs, int n, sollya_mpfi_t x);
/*************************************************************************************************************************************/

/*****************************************************************************/
/*************Functions related to bounding polynomials in ChebBasis**********/
/*****************************************************************************/

/* This function computes an interval bound for a polynomial in cheb basis. */
/*The coefficients are given in coeffs. We have n coeffs*/
/* by just adding the absolute values of the coeffs        */

void chebPolynomialBoundSimple(sollya_mpfi_t bound, int n, sollya_mpfi_t *coeffs);


/* This function computes an interval bound for a polynomial in cheb basis. */
/*The coefficients are given in coeffs. We have n coeffs*/
/* we will use a refined method for this        */
void chebPolynomialBoundRefined(sollya_mpfi_t bound, int n, sollya_mpfi_t *coeffs);

/* This function computes an interval bound for a polynomial in cheb basis. */
/*One day this function may become more complex*/
void chebPolynomialBoundDefault(sollya_mpfi_t bound, int n, sollya_mpfi_t *coeffs);

/* This function computes an interval bound for a polynomial in cheb basis. */
/*The coefficients are given in coeffs. We have n coeffs*/
/* by just adding the absolute values of the coeffs        */

void evaluateChebPolynomialClenshaw(sollya_mpfi_t bound, int n, sollya_mpfi_t *coeffs, mpfi_t x,mpfi_t x0);

/*wrapper to get directly the coeffs in the monomial basis from a polynomial in the chebBasis basis(given a pointer to node, over a given interval x*/
void getCoeffsFromChebPolynomial(sollya_mpfi_t**coeffs, sollya_mpfi_t *chebCoeffs, int n, sollya_mpfi_t x);



