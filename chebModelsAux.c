/*For compiling this file:
    gcc -fPIC -Wall -c chebModelsAux.c
    gcc -fPIC -shared -o chebModelsAux chebModelsAux.o
    
    From sollya:
    externalproc(CM, "./chebModelsAux", (function, range, integer,integer) ->list of range);

*/



//extern void fprintInterval(FILE *fd, sollya_mpfi_t interval);
//extern void pushTimeCounter(void);
//extern void popTimeCounter(char *s);

#include  "chebModelsAux.h"
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/*This function gets an mpfi from a node c:
if c constant after simplification-> ok
if not 0
*/
int mpfi_set_node( sollya_mpfi_t *r, node * c) {
  sollya_mpfi_t rr;
   sollya_mpfi_t *rrr;
  node *cc;
  sollya_mpfi_t dummy;
  sollya_mpfi_init2(rr,getToolPrecision());
  sollya_mpfi_init2(dummy,getToolPrecision());
  rrr=safeMalloc(sizeof(sollya_mpfi_t));
  sollya_mpfi_init2(*rrr,getToolPrecision());
  if (c!=NULL){
    cc=simplifyTreeErrorfree(c);
    switch (cc->nodeType){
      case PI_CONST: mpfi_const_pi(rr);
      break;
      case CONSTANT:mpfi_set_fr(rr,*(cc->value));
      break;
      default:  auto_diff(rrr,c,dummy,0); mpfi_set(rr, *rrr);
      break;
      }
    free(cc);
  }
  else mpfi_set_ui(rr,0);
  mpfi_set(*r,rr);
  //printf("\n in mpfi set node\n");
  //printInterval(rr);
  mpfi_clear(rr);
  mpfi_clear(dummy);
  mpfi_clear(*rrr);
  free(rrr);
  return 0;
}



/**************************************************************************************/
/**************************************************************************************/
/***************************Functions related to Chebyshev points**********************/
/**************************************************************************************/
/**************************************************************************************/


/*returns n chebyshev points in X*/

void getChebyshevPoints(sollya_mpfi_t *chebPoints, int n, sollya_mpfi_t x){
int i;
mpfr_t u, v;

sollya_mpfi_t ui, vi, temp1, temp2, mpfiPi, mpfiPiArg;

/*print("The chebyshev points are:");
  };
  for i from 1 to n do{  
    chebList[i-1]=(cos((2*i-1)*Pi/(2*n)))*(u-v)/2 + (u+v)/2;
  };
*/
sollya_mpfi_init2(ui,getToolPrecision());
sollya_mpfi_init2(vi, getToolPrecision());
sollya_mpfi_init2(temp1, getToolPrecision());
sollya_mpfi_init2(temp2, getToolPrecision());
sollya_mpfi_init2(mpfiPi, getToolPrecision());
sollya_mpfi_init2(mpfiPiArg, getToolPrecision());

mpfr_init2(u, getToolPrecision());
mpfr_init2(v, getToolPrecision());

sollya_mpfi_get_left(u,x);
sollya_mpfi_get_right(v,x);

sollya_mpfi_set_fr(ui,u);
sollya_mpfi_set_fr(vi,v);

sollya_mpfi_sub(temp1, ui, vi);
sollya_mpfi_div_ui(temp1, temp1,2);



sollya_mpfi_add(temp2, ui, vi);
sollya_mpfi_div_ui(temp2, temp2,2);


sollya_mpfi_const_pi(mpfiPi);
sollya_mpfi_div_ui(mpfiPi,mpfiPi,2*n);

for (i=1;i<=n;i++){
//if (2*i-1<n) 
sollya_mpfi_mul_ui(mpfiPiArg,mpfiPi,(2*i-1));
//else{
//sollya_mpfi_mul_ui(mpfiPiArg,mpfiPi,2*n-1-j);
//j=j+2;
//}

sollya_mpfi_cos(chebPoints[i-1],mpfiPiArg);

sollya_mpfi_mul( chebPoints[i-1], chebPoints[i-1], temp1);
sollya_mpfi_add( chebPoints[i-1], chebPoints[i-1], temp2);

}

sollya_mpfi_clear(ui);
sollya_mpfi_clear(vi);
sollya_mpfi_clear(temp1);
sollya_mpfi_clear(temp2);
sollya_mpfi_clear(mpfiPi);
sollya_mpfi_clear(mpfiPiArg);

mpfr_clear(u);
mpfr_clear(v);
}




void getChebyshevExtrema(sollya_mpfi_t *chebPoints, int n, sollya_mpfi_t x){
int i;
mpfr_t u, v;

sollya_mpfi_t ui, vi, temp1, temp2, mpfiPi, mpfiPiArg;

/*print("The chebyshev extremas are:");
  };
  for i from 1 to n-1 do{  
    chebExtremas[i-1]=(cos((i)*Pi/(n)))*(u-v)/2 + (u+v)/2;
  };
*/

sollya_mpfi_init2(ui,getToolPrecision());
sollya_mpfi_init2(vi, getToolPrecision());
sollya_mpfi_init2(temp1, getToolPrecision());
sollya_mpfi_init2(temp2, getToolPrecision());
sollya_mpfi_init2(mpfiPi, getToolPrecision());
sollya_mpfi_init2(mpfiPiArg, getToolPrecision());

mpfr_init2(u, getToolPrecision());
mpfr_init2(v, getToolPrecision());

sollya_mpfi_get_left(u,x);
sollya_mpfi_get_right(v,x);

sollya_mpfi_set_fr(ui,u);
sollya_mpfi_set_fr(vi,v);

sollya_mpfi_sub(temp1, ui, vi);
sollya_mpfi_div_ui(temp1, temp1,2);



sollya_mpfi_add(temp2, ui, vi);
sollya_mpfi_div_ui(temp2, temp2,2);


sollya_mpfi_const_pi(mpfiPi);
sollya_mpfi_div_ui(mpfiPi,mpfiPi,n);


for (i=1;i<=n-1;i++){
sollya_mpfi_mul_ui(mpfiPiArg,mpfiPi,i);
sollya_mpfi_cos(chebPoints[i-1],mpfiPiArg);

sollya_mpfi_mul( chebPoints[i-1], chebPoints[i-1], temp1);
sollya_mpfi_add( chebPoints[i-1], chebPoints[i-1], temp2);

}

sollya_mpfi_clear(ui);
sollya_mpfi_clear(vi);
sollya_mpfi_clear(temp1);
sollya_mpfi_clear(temp2);
sollya_mpfi_clear(mpfiPi);
sollya_mpfi_clear(mpfiPiArg);

mpfr_clear(u);
mpfr_clear(v);
}
/*The so called chebMatrix, contains 
the specific values of T_i(x_j) i=0..n-1, j=0..n-1
x_j \in [-1,1]. 
Note: It has to be computed only once for each n*/

void getChebMatrix(sollya_mpfi_t**chebMatrix, int n){
int i,j;
sollya_mpfi_t *chebPoints;
sollya_mpfi_t intrval;

sollya_mpfi_t temp;
sollya_mpfi_init2(temp, getToolPrecision());

chebPoints=safeMalloc((n)*sizeof(sollya_mpfi_t));
  for (i=0;i<n;i++){
  sollya_mpfi_init2(chebPoints[i],getToolPrecision());
  }

sollya_mpfi_init2(intrval,getToolPrecision());
sollya_mpfi_interv_si(intrval,-1,1);

getChebyshevPoints(chebPoints, n, intrval);

*chebMatrix= (sollya_mpfi_t *)safeMalloc((n*n)*sizeof(sollya_mpfi_t));  

for (i=0;i<n;i++){
  for (j=0;j<n;j++){
    sollya_mpfi_init2((*chebMatrix)[i*n+j],getToolPrecision());
  }
}

for (i=0;i<n;i++){
sollya_mpfi_set_ui((*chebMatrix)[i],1);
}

for (i=0;i<n;i++){
sollya_mpfi_set((*chebMatrix)[i+n],chebPoints[i]);
}

for (i=2;i<n;i++){
  for (j=0;j<n;j++){

  sollya_mpfi_mul(temp,(*chebMatrix)[(i-1)*n+j], (*chebMatrix)[n+j]);
  sollya_mpfi_mul_ui(temp,temp,2);
  sollya_mpfi_sub((*chebMatrix)[i*n+j],temp,(*chebMatrix)[(i-2)*n+j]);
  }
}


sollya_mpfi_clear(intrval);
sollya_mpfi_clear(temp);
for (i=0;i<n;i++){
  sollya_mpfi_clear(chebPoints[i]);
  }
free(chebPoints);
}

void getChebPolyCoeffs(mpz_t* chebMatrix, int n){
int i,j;
mpz_t temp;
mpz_init2(temp, getToolPrecision());
//printf("in getChebPolyCoeffs, %d\n",n);
for (i=0;i<n;i++)
  for (j=0;j<n;j++)
    mpz_set_ui(chebMatrix[i*n+j],0);

if (n>0) mpz_set_ui(chebMatrix[0],1);
if (n>1) mpz_set_ui(chebMatrix[n+1],1);

for (i=2;i<n;i++){
  mpz_neg(chebMatrix[i*n],chebMatrix[(i-2)*n]); //put constant coeff of T_{n-2}
  for (j=0;j<i;j++){
  mpz_mul_ui(temp,chebMatrix[(i-1)*n+j],2);
  mpz_sub(chebMatrix[i*n+j+1],temp,chebMatrix[(i-2)*n+j+1]);
  }
}

mpz_clear(temp);

}

/*returns in chebCoeffs the coeffs of the finite Chebyshev basis expansion of the polynomial given by the monomial coefficients p, of degree n-1*/
void getPolyCoeffsChebBasis(sollya_mpfi_t *chebCoeffs, sollya_mpfi_t *p, int n){
sollya_mpfi_t *pAux, temp;
mpz_t *chebMatrix;
int i,j;
  //printf("\nIn getPolyCoeffsChebBasis: %d \n",n);
  pAux=safeMalloc((n)*sizeof(sollya_mpfi_t));

  for (i=0;i<n;i++){
    //printInterval(p[i]);
    sollya_mpfi_init2(pAux[i],getToolPrecision());
    sollya_mpfi_set(pAux[i], p[i]);
  }
  
  chebMatrix=(mpz_t *)safeMalloc((n*n)*sizeof(mpz_t));  
  
  for (i=0;i<n*n;i++){
    mpz_init2(chebMatrix[i],getToolPrecision());
  }
  
  getChebPolyCoeffs(chebMatrix, n);
  sollya_mpfi_init2(temp,getToolPrecision());  
  
  for(i=n-1; i>=0;i--){
    
    mpfi_div_z(chebCoeffs[i],pAux[i],chebMatrix[i*n+i]);
    for(j=i-1;j>=0;j--){
      mpfi_mul_z(temp,chebCoeffs[i],chebMatrix[i*n+j]);
      mpfi_sub(pAux[j],pAux[j],temp);
    }
  }
  for (i=0;i<n;i++){
    sollya_mpfi_clear(pAux[i]);
  }
  free(pAux);
  
  for (i=0;i<n*n;i++){
    mpz_clear(chebMatrix[i]);
  }
  free(chebMatrix);
  
  sollya_mpfi_clear(temp);
  //printf("\nOut of getPolyCoeffsChebBasis: %d \n",n);
}

/*returns in c the coeffs in the monomial basis for the polynomial p(a*x+b), where a and b are mpfi_s
 the polynomial p  of degree n-1 given by the monomial coefficients stored as mpfi_s in p*/
void getTranslatedPolyCoeffs(sollya_mpfi_t *c, sollya_mpfi_t *p, int n, sollya_mpfi_t a, sollya_mpfi_t b){
sollya_mpfi_t *pAux, temp, pow, alpha, evalP;

int i;
 // printf("\nIn getTranslatedPolyCoeffs: %d \n",n);
  pAux=safeMalloc((n)*sizeof(sollya_mpfi_t));
  sollya_mpfi_init2(temp,getToolPrecision());  
  sollya_mpfi_init2(pow,getToolPrecision());  
  sollya_mpfi_init2(alpha,getToolPrecision());  
  sollya_mpfi_init2(evalP,getToolPrecision());  
  for (i=0;i<n;i++){
    //printInterval(p[i]);
    sollya_mpfi_init2(pAux[i],getToolPrecision());
    sollya_mpfi_set(pAux[i], p[i]);
  }
  
  //printInterval(a);
  //printInterval(b);  
  /*do the transformation : P(a*x+b) = \sum p_i*a^i *(x+b/a)^i*/
  for (i=0;i<n;i++){
    sollya_mpfi_set_ui(pow,i);
    sollya_mpfi_pow(temp, a,pow);
    sollya_mpfi_mul(pAux[i], pAux[i],temp);
  }
  sollya_mpfi_div(alpha, b,a);
  /*do the translation from P(x+alpha) --> Q(x)*/
  /*P(x+alpha)= \sum P^(i)(alpha)/i! *x^i */
  /*we have to compute the first n-1 derivatives in alpha: P(alpha), P'(alpha), ..., P^(n-1)(alpha)*/
  /*evaluate P in alpha, P is of degree n-1*/
  symbolic_poly_evaluation_horner(evalP, pAux,alpha,n-1);
  sollya_mpfi_set(c[0],evalP);
  mpfi_set_ui(temp, 1);
  for(i=1;i<n;i++){
   symbolic_poly_diff(pAux, pAux, n-i); //differentiate in place pAux
   symbolic_poly_evaluation_horner(evalP, pAux,alpha,n-i-1);//evaluate P^(i)(alpha)
   sollya_mpfi_mul_ui(temp,temp,i);   
   sollya_mpfi_div(c[i],evalP, temp);
  }
  

  for (i=0;i<n;i++){
    sollya_mpfi_clear(pAux[i]);
  }
  free(pAux);
  
  sollya_mpfi_clear(temp);
  sollya_mpfi_clear(pow);
  sollya_mpfi_clear(alpha);
  sollya_mpfi_clear(evalP);
 // printf("\nOut of getTranslatedPolyCoeffs: %d \n",n);
}

/***************************************************************/
/***************************************************************/
/**********Functions related to Obtaining Chebyshev Coeffs*********/
/***************************************************************/
/***************************************************************/


/*evaluate a function in n chebPoints*/
void getFunctionValues(sollya_mpfi_t* fValues, sollya_mpfi_t * chebPoints,node* f,int n){
int i;
for (i=0;i<n;i++){
  //evaluateMpfiFunction(fValues[i],f,chebPoints[i],getToolPrecision());
  
  auto_diff(&fValues[i],f,chebPoints[i],0);
  }

}


/*this function computes the cheb coeffs for the interpolation polynomial of degree n-1 at chebpoints (we have n chebpoints and a n*n chebMatrix)*/
void getChebCoeffs(sollya_mpfi_t* coeffs, sollya_mpfi_t *chebMatrix, sollya_mpfi_t *fValues,int n){

int i,j;
sollya_mpfi_t temp;
sollya_mpfi_init2(temp, getToolPrecision());

for (i=0;i<n;i++)
    sollya_mpfi_set_ui(coeffs[i], 0);

for (i=0;i<n;i++)
    sollya_mpfi_add(coeffs[0],coeffs[0],fValues[i]);
sollya_mpfi_div_ui(coeffs[0],coeffs[0], n);


for(i=1;i<n;i++){
  for(j=0;j<n;j++){
    sollya_mpfi_mul(temp, fValues[j], chebMatrix[i*n+j]);
    sollya_mpfi_add(coeffs[i],coeffs[i],temp);
  } 
sollya_mpfi_mul_ui(coeffs[i],coeffs[i], 2);
sollya_mpfi_div_ui(coeffs[i],coeffs[i], n);
}

sollya_mpfi_clear(temp);
}

/*wrapper to get directly the coeffs from the function*/
void getChebCoeffsFromFunction(sollya_mpfi_t* coeffs, sollya_mpfi_t * chebPoints, sollya_mpfi_t * chebMatrix,node *f,int n){

sollya_mpfi_t *fValues ;
int i;

fValues=safeMalloc((n+1)*sizeof(sollya_mpfi_t));
  for (i=0;i<=n;i++){
  sollya_mpfi_init2(fValues[i],getToolPrecision());
  }
  getFunctionValues(fValues, chebPoints,f,n);
  getChebCoeffs(coeffs,chebMatrix,fValues,n);

for (i=0;i<=n;i++){
  sollya_mpfi_clear(fValues[i]);
  }
free(fValues);

}

/*wrapper to get directly the coeffs in the chebyshev basis from a polynomial in the monomial basis(given a pointer to node, over a given interval x*/
void getChebCoeffsFromPolynomial(sollya_mpfi_t**coeffs, int*n, node *f, sollya_mpfi_t x){
  sollya_mpfi_t z1, z2, ui, vi;
  node **coefficients;
  int d,i;
  sollya_mpfi_t *p, *c;
  mpfr_t u,v;
  int verbosity =0;
  if (isPolynomial(f) ){
    //printf("We found a polynomial!!!\n");
    getCoefficients(&d, &coefficients, f);
    *n=d+1;
   
    *coeffs= (sollya_mpfi_t *)safeMalloc((d+1)*sizeof(sollya_mpfi_t));  
   
    p=safeMalloc((d+1)*sizeof(sollya_mpfi_t));
    c=safeMalloc((d+1)*sizeof(sollya_mpfi_t));
    for (i=0;i<d+1;i++){
      sollya_mpfi_init2((*coeffs)[i],getToolPrecision());
      sollya_mpfi_init2(p[i],getToolPrecision());
      sollya_mpfi_init2(c[i],getToolPrecision());
      mpfi_set_node(&(p[i]),coefficients[i]);   
    }
   
    /*Here we have the coeffs of the polynomial in p, over the interval x=[a,b]*/
    /*we need to compute the polynomial over [-1,1]*/
    /*we make the change of variable: x= y*(b-a)/2 + (b+a)/2, hence for y \in [-1,1] we have x\in [a,b]*/
    /* we compute P(x)=Q(y)*/
    sollya_mpfi_init2(ui, getToolPrecision());
    sollya_mpfi_init2(vi, getToolPrecision());
    
    
    mpfr_init2(u, getToolPrecision());
    mpfr_init2(v, getToolPrecision());
    
    sollya_mpfi_init2(z1, getToolPrecision());
    sollya_mpfi_init2(z2, getToolPrecision());
    
    sollya_mpfi_get_left(u,x);
    sollya_mpfi_get_right(v,x);

    sollya_mpfi_set_fr(ui,u);
    sollya_mpfi_set_fr(vi,v);
    
    sollya_mpfi_add(z2,ui,vi);
    sollya_mpfi_sub(z1,vi,ui);
    
    sollya_mpfi_div_ui(z1,z1,2);
    sollya_mpfi_div_ui(z2,z2,2);
    
    getTranslatedPolyCoeffs(c, p, d+1, z1,z2);
    if (verbosity>10) {
    printf("\nThe %d coefficients of the translated poly are :\n ",d+1);
      for (i=0;i<d+1;i++){
      printInterval(c[i]);
      }
    }
    
    getPolyCoeffsChebBasis(*coeffs, c, d+1);
    if (verbosity>10) {
    printf("\nThe %d coefficients in cheby basis are :\n ",d+1);
      for (i=0;i<d+1;i++){
      printInterval((*coeffs)[i]);
      }
    }
  }
  else{
     printf("The given function is not a polynomial, no modif is made \n");  
  }
}


/*wrapper to get directly the coeffs in the monomial basis from a polynomial in the chebBasis basis(given a pointer to node, over a given interval x*/
void getCoeffsFromChebPolynomial(sollya_mpfi_t**coeffs, sollya_mpfi_t *chebCoeffs, int n, sollya_mpfi_t x){
  sollya_mpfi_t z1, z2, ui, vi, temp;
 
  int j,i;
  sollya_mpfi_t *c;
  mpfr_t u,v;
  int verbosity =0;
  mpz_t *chebMatrix;
  
    sollya_mpfi_init2(temp, getToolPrecision());
   
   chebMatrix= (mpz_t *)safeMalloc((n*n)*sizeof(mpz_t));  
   
   for (i=0;i<n*n;i++){
    mpz_init2(chebMatrix[i], getToolPrecision());
   }
   getChebPolyCoeffs(chebMatrix, n);
   //printf("After computing cheby coeffs");
   
   *coeffs= (sollya_mpfi_t *)safeMalloc((n)*sizeof(sollya_mpfi_t));  
   c=(sollya_mpfi_t *)safeMalloc((n)*sizeof(sollya_mpfi_t));
   
   for (i=0;i<n;i++){
      sollya_mpfi_init2((*coeffs)[i],getToolPrecision());
      sollya_mpfi_init2(c[i],getToolPrecision());
      sollya_mpfi_set_ui(c[i],0);
   }
   
   for (j=0;j<n;j++){
    for (i=j;i<n;i++){
      mpfi_mul_z(temp, chebCoeffs[i], chebMatrix[i*n+j]);
      sollya_mpfi_add(c[j], c[j], temp);  
    }
   }
   
   /*we have in c_i the values of the coefs of P(2/(b-a)x- (b+a)/(b-a)) = \sum c_i (2/(b-a)x- (b+a)/(b-a))^i*/   
   /*we need to the translation*/
    
    /*we compute z1=2/(b-a); z2=-(b+a)/(b-a)*/
   
    sollya_mpfi_init2(ui, getToolPrecision());
    sollya_mpfi_init2(vi, getToolPrecision());
    
    
    mpfr_init2(u, getToolPrecision());
    mpfr_init2(v, getToolPrecision());
    
    sollya_mpfi_init2(z1, getToolPrecision());
    sollya_mpfi_init2(z2, getToolPrecision());
    
    sollya_mpfi_get_left(u,x);
    sollya_mpfi_get_right(v,x);

    sollya_mpfi_set_fr(ui,u);
    sollya_mpfi_set_fr(vi,v);
    
    sollya_mpfi_sub(z2,vi,ui);
        
    sollya_mpfi_ui_div(z1,2,z2);
    sollya_mpfi_add(temp, ui, vi);
    sollya_mpfi_div(z2, temp, z2);
    sollya_mpfi_neg(z2, z2);    
        getTranslatedPolyCoeffs((*coeffs), c, n, z1,z2);
    if (verbosity>10) {
    printf("\nThe %d coefficients of the translated poly are :\n ",n);
      for (i=0;i<n;i++){
      printInterval((*coeffs)[i]);
      }
    }
   
}
  



/*wrapper to get directly the coeffs in the chebyshev basis up to degree n-1 (first n coeffs) and a bound for the remaining polynomial, from a polynomial in the monomial basis(given a pointer to node, over a given interval x*/
void getNChebCoeffsFromPolynomial(sollya_mpfi_t *coeffs, sollya_mpfi_t bound, node *f, sollya_mpfi_t x, int n){
  
  sollya_mpfi_t **c;
  int d,i;
    
  c= (sollya_mpfi_t **)safeMalloc(sizeof(sollya_mpfi_t*));  
  getChebCoeffsFromPolynomial(c, &d, f, x);
 
  if (d<=n) {
    for(i=0;i<d;i++)
      sollya_mpfi_set(coeffs[i],(*c)[i]);
    for(i=d;i<n;i++)
    sollya_mpfi_set_ui(coeffs[i],0);
    
    sollya_mpfi_set_ui(bound,0);
  }
  else{
     for(i=0;i<n;i++)
       sollya_mpfi_set(coeffs[i],(*c)[i]);
    //printf("before bounding");
    
    chebPolynomialBoundSimple(bound, d-n, (&(*c)[n]));
  }
  
}


/********************************Functions related to derivation and integration of polynomials in cheb basis*************************/

void getChebCoeffsDerivativePolynomial(sollya_mpfi_t*coeffs, sollya_mpfi_t *chebCoeffs, int n, sollya_mpfi_t x){
  sollya_mpfi_t z1, z2, ui, vi, temp;
 
  int i;
  sollya_mpfi_t *c;
  mpfr_t u,v;
  int verbosity=0;
 
    c=(sollya_mpfi_t *)safeMalloc((n-1)*sizeof(sollya_mpfi_t));
    sollya_mpfi_init2(temp, getToolPrecision());
   
  
   for (i=0;i<n-1;i++){
      sollya_mpfi_init2(c[i],getToolPrecision());
      sollya_mpfi_set_ui(c[i],0);
   }
    
   if(n>1) {
    sollya_mpfi_mul_ui(c[n-2],chebCoeffs[n-1],2*(n-1));
   }
   if(n>2) {
    sollya_mpfi_mul_ui(c[n-3],chebCoeffs[n-2],2*(n-2));
   }
   for (i=n-3;i>0;i--){
     sollya_mpfi_mul_ui(c[i-1],chebCoeffs[i],2*i);
     sollya_mpfi_add(c[i-1],c[i-1],c[i+1]);
   }
   sollya_mpfi_div_ui(c[0],c[0],2);
   
   
   /*we have in c_i the values of the coefs of P'(y) = \sum c_i T_i(x)*/   
   /*we have to multiply by y'(x), which is z1=2/(b-a) */
    
    /*we compute z1=2/(b-a)*/
   
    sollya_mpfi_init2(ui, getToolPrecision());
    sollya_mpfi_init2(vi, getToolPrecision());
    
    
    mpfr_init2(u, getToolPrecision());
    mpfr_init2(v, getToolPrecision());
    
    sollya_mpfi_init2(z1, getToolPrecision());
    sollya_mpfi_init2(z2, getToolPrecision());
    
    sollya_mpfi_get_left(u,x);
    sollya_mpfi_get_right(v,x);

    sollya_mpfi_set_fr(ui,u);
    sollya_mpfi_set_fr(vi,v);
    
    sollya_mpfi_sub(z2,vi,ui);
        
    sollya_mpfi_ui_div(z1,2,z2);
   
     for (i=0;i<n-1;i++){
      sollya_mpfi_mul(c[i], c[i],z1);
     }
    
    if (verbosity>10) {
    printf("\nThe %d coefficients of the derivated polynomial are :\n ",n-1);
      for (i=0;i<n-1;i++){
      printInterval(c[i]);
      }
    }
    for (i=0;i<n-1;i++){
      sollya_mpfi_set(coeffs[i], c[i]);
    }
    for (i=0;i<n-1;i++){
      sollya_mpfi_clear(c[i]);
    }
    free(c);
    sollya_mpfi_clear(z1);
    sollya_mpfi_clear(z2);
    sollya_mpfi_clear(ui);
    sollya_mpfi_clear(vi);
    mpfr_clear(u);
    mpfr_clear(v);
}

/*Computes the antiderivative of a polynomial in Chebyshev basis.
NOTE: the constant coefficient is set to zero, but it should be viewed as a constant*/
void getChebCoeffsIntegrationPolynomial(sollya_mpfi_t*coeffs, sollya_mpfi_t *chebCoeffs, int n, sollya_mpfi_t x){
  sollya_mpfi_t z1, z2, ui, vi, temp;
 
  int i;
  sollya_mpfi_t *c;
  mpfr_t u,v;
  int verbosity=0;
 
    c=(sollya_mpfi_t *)safeMalloc((n+1)*sizeof(sollya_mpfi_t));
    sollya_mpfi_init2(temp, getToolPrecision());
   
  
   for (i=0;i<n+1;i++){
      sollya_mpfi_init2(c[i],getToolPrecision());
      sollya_mpfi_set_ui(c[i],0);
   }
   
   if(n>0){
   sollya_mpfi_div_ui(c[1],chebCoeffs[2],2);
   sollya_mpfi_sub(c[1],chebCoeffs[0], c[1]);
   }
   
    
   for (i=2;i<n-1;i++){
     sollya_mpfi_sub(c[i],chebCoeffs[i-1],chebCoeffs[i+1]);
     sollya_mpfi_div_ui(c[i],c[i],2*i);
    
   }
   
   
   if(n>1){
   sollya_mpfi_set(c[n-1],chebCoeffs[n-2]);
   sollya_mpfi_div_ui(c[n-1],c[n-1],2*(n-1));
   }
   
   if(n>0){
   sollya_mpfi_set(c[n],chebCoeffs[n-1]);
   sollya_mpfi_div_ui(c[n],c[n],2*(n));
   }
   
   
   
   /*we have in c_i the values of the coefs of \int P(y) = \sum c_i T_i(x) (the constant of integration in c_0 is not computed*/   
   /*we have to multiply by 1/y'(x), which is z1=(b-a)/2 */
    
    /*we compute z1=(b-a)/2*/
   
    sollya_mpfi_init2(ui, getToolPrecision());
    sollya_mpfi_init2(vi, getToolPrecision());
    
    
    mpfr_init2(u, getToolPrecision());
    mpfr_init2(v, getToolPrecision());
    
    sollya_mpfi_init2(z1, getToolPrecision());
    sollya_mpfi_init2(z2, getToolPrecision());
    
    sollya_mpfi_get_left(u,x);
    sollya_mpfi_get_right(v,x);

    sollya_mpfi_set_fr(ui,u);
    sollya_mpfi_set_fr(vi,v);
    
    sollya_mpfi_sub(z2,vi,ui);
        
    sollya_mpfi_div_ui(z1,z2,2);
   
     for (i=1;i<n+1;i++){
      sollya_mpfi_mul(c[i], c[i],z1);
     }
    
    if (verbosity>10) {
    printf("\nThe %d coefficients of the integrated polynomial are :\n ",n+1);
      for (i=0;i<n+1;i++){
      printInterval(c[i]);
      }
    }
    for (i=0;i<n+1;i++){
      sollya_mpfi_set(coeffs[i], c[i]);
    }
    for (i=0;i<n+1;i++){
      sollya_mpfi_clear(c[i]);
    }
    free(c);
    sollya_mpfi_clear(z1);
    sollya_mpfi_clear(z2);
    sollya_mpfi_clear(ui);
    sollya_mpfi_clear(vi);
    mpfr_clear(u);
    mpfr_clear(v);
}


 
/*****************************************************************************/
/*************Functions related to bounding polynomials in ChebBasis**********/
/*****************************************************************************/

/* This function computes an interval bound for a polynomial in cheb basis. */
/*The coefficients are given in coeffs. We have n coeffs*/
/* by just adding the absolute values of the coeffs        */

void chebPolynomialBoundSimple(sollya_mpfi_t bound, int n, sollya_mpfi_t *coeffs){
  sollya_mpfi_t temp, temp2, intrval;
  mp_prec_t prec;
  int i;
  //printf("\nin chebPolynomialBoundSimple\n");
  //for(i=0;i<n;i++) {
   // printInterval(coeffs[i]);
 // }
  prec = getToolPrecision();
  sollya_mpfi_init2(temp, prec);
  sollya_mpfi_init2(temp2, prec);
  sollya_mpfi_init2(intrval, prec);
  sollya_mpfi_set_ui(temp, 0);
  sollya_mpfi_interv_si(intrval, -1,1);
  if (n>0) sollya_mpfi_set(temp, coeffs[0]);
  for(i=1;i<n;i++) {
  sollya_mpfi_mul(temp2,coeffs[i], intrval);
  sollya_mpfi_add(temp, temp,temp2);
  }
  sollya_mpfi_set(bound, temp);
  sollya_mpfi_clear(temp);
  sollya_mpfi_clear(temp2);
  sollya_mpfi_clear(intrval);
  //printf("out of chebPolynomail bound simple");
}

/* This function computes an interval bound for a polynomial in cheb basis. */
/*The coefficients are given in coeffs. We have n coeffs*/
/* TODO      */

void chebPolynomialBoundRefined(sollya_mpfi_t bound, int n, sollya_mpfi_t *coeffs){
  sollya_mpfi_t temp, temp2, intrval;
  mp_prec_t prec;
  int i;
  //printf("\nin chebPolynomialBoundSimple\n");
  for(i=0;i<n;i++) {
    //printInterval(coeffs[i]);
  }
  prec = getToolPrecision();
  sollya_mpfi_init2(temp, prec);
  sollya_mpfi_init2(temp2, prec);
  sollya_mpfi_init2(intrval, prec);
  sollya_mpfi_set_ui(temp, 0);
  sollya_mpfi_interv_si(intrval, -1,1);
  if (n>0) sollya_mpfi_set(temp, coeffs[0]);
  for(i=1;i<n;i++) {
  sollya_mpfi_mul(temp2,coeffs[i], intrval);
  sollya_mpfi_add(temp, temp,temp2);
  }
  sollya_mpfi_set(bound, temp);
  sollya_mpfi_clear(temp);
  sollya_mpfi_clear(temp2);
  sollya_mpfi_clear(intrval);
}


/* This function computes an interval bound for a polynomial in cheb basis. */
/*One day this function may become more complex*/
void chebPolynomialBoundDefault(sollya_mpfi_t bound, int n, sollya_mpfi_t *coeffs){
chebPolynomialBoundSimple(bound, n,coeffs);
}


/* This function computes an interval bound for a polynomial in cheb basis. */
/*The coefficients are given in coeffs. We have n coeffs*/
/* by using Clenshaw's method --      */

void evaluateChebPolynomialClenshaw(sollya_mpfi_t bound, int n, sollya_mpfi_t *coeffs, mpfi_t x,mpfi_t x0 ){
  mp_prec_t prec;
  int i;
  sollya_mpfi_t z, zz, z1,b0,b1;
  mpfr_t a, b;
  prec = getToolPrecision();
  sollya_mpfi_init2(z, prec);
  sollya_mpfi_init2(zz, prec);
  sollya_mpfi_init2(z1, prec);
  sollya_mpfi_init2(b0, prec);
  sollya_mpfi_init2(b1, prec);
  mpfr_init2(a, prec);
  mpfr_init2(b, prec);
  
  sollya_mpfi_get_right(b,x);
  sollya_mpfi_get_left(a,x);
  sollya_mpfi_set_fr(z1,b);
  sollya_mpfi_sub_fr(z1,z1,a);
  sollya_mpfi_inv(z1,z1);
  sollya_mpfi_mul_ui(z,z1,2);
  /*z=2/(b-a)*/
  sollya_mpfi_set_fr(zz,b);
  mpfi_add_fr(zz,zz,a);
  sollya_mpfi_mul(zz,zz,z1);
  /*zz=(b+a)/(b-a)*/
  sollya_mpfi_mul(z,z,x0);
  sollya_mpfi_sub(z,z,zz);
  
  /*z=2/(b-a) * x0 - (b+a)/(b-a)*/
  
  /*Do the clenshaw algo*/
  /*b1:=[0.,0.];
    b0:=[0.,0.];
    for i from n to 2 by -1 do
      bb:=(((2&*z)&*b0) &-b1) &+L[i];
      b1:=b0;
      b0:=bb;
    end do;
    bb:=((z&*b0) &-b1) &+L[1];
    return bb;
  */
  sollya_mpfi_set_ui(b0,0);  
  sollya_mpfi_set_ui(b1,0);
  
  for(i=n-1;i>0; i--){
    sollya_mpfi_mul(zz,z,b0);
    sollya_mpfi_mul_ui(zz,zz,2);
    sollya_mpfi_sub(zz,zz,b1);
    sollya_mpfi_add(zz,zz,coeffs[i]);
    sollya_mpfi_set(b1,b0);
    sollya_mpfi_set(b0,zz);  
  }
  sollya_mpfi_mul(zz,z,b0);
  sollya_mpfi_sub(zz,zz,b1);
  sollya_mpfi_add(zz,zz,coeffs[0]);
  sollya_mpfi_set(bound, zz);  
  
  sollya_mpfi_clear(zz);
  sollya_mpfi_clear(z);
  sollya_mpfi_clear(z1);
  sollya_mpfi_clear(b0);
  sollya_mpfi_clear(b1);
  mpfr_clear(b); mpfr_clear(a);
}

/*
int CM(chain**resP, void **args) {
  node *f;
  sollya_mpfi_t x, bound;
  sollya_mpfi_t *coeffsF, *fValues, *chebMatrix, *chebCoeffs;
  mpz_t *chebCoeffMatrix;
  //sollya_mpfi_t *coeffsErrors, *rest;
  int n,i,j,verbosity;
  int d;

  //chain *ch;
  //node *resultPoly;
  //node *remainderL, *remainderR;
  //mpfr_t *ptr;
  //ch=NULL;
  sollya_mpfi_t *chebPoints;
  
  chain *ch;
 
  sollya_mpfi_t *zz;
  //sollya_mpfi_t *ptr;
  
  f = (node *)args[0];
  n = *( (int *)args[2] );
  verbosity=*( (int *)args[3] );
  sollya_mpfi_init2(x, getToolPrecision());
  sollya_mpfi_set(x, *( (sollya_mpfi_t *)args[1] ));
  
  chebCoeffs=safeMalloc((n)*sizeof(sollya_mpfi_t));
  
  for (i=0;i<n;i++){
  sollya_mpfi_init2(chebCoeffs[i],getToolPrecision());
  }
  sollya_mpfi_init2(bound, getToolPrecision());
  getNChebCoeffsFromPolynomial(chebCoeffs, bound, f, x, n);
  
    ch=NULL;
    
    zz = (sollya_mpfi_t*)safeMalloc(sizeof(sollya_mpfi_t));
    sollya_mpfi_init2(*zz, getToolPrecision());
    sollya_mpfi_set(*zz, bound);
    ch=addElement(ch,zz);
    
    for (i=n-1; i>=0;i--){
      zz = (sollya_mpfi_t*)safeMalloc(sizeof(sollya_mpfi_t));
      sollya_mpfi_init2(*zz, getToolPrecision());
      sollya_mpfi_set(*zz, chebCoeffs[i]);
      ch=addElement(ch,zz);
    }
    for (i=0;i<n;i++){
      sollya_mpfi_clear(chebCoeffs[i]);
    }
    free(chebCoeffs);
  
 
  
  *resP=ch;
    
  return 1;
}*/
