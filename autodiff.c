/*
  For compiling this file:
    gcc -fPIC -Wall -c autodiff.c
    gcc -fPIC -shared -o autodiff autodiff.o


  Within Sollya:
    > externalproc(AD, "./autodiff", (function, range, integer) -> list of range);

  And then, for instance:
    > AD(exp(x)*cos(x), [2.5; 2.6], 10);

*/


#include "sollya.h"

extern int mpfi_pow(mpfi_t res, mpfi_t x, mpfi_t y);
extern void fprintInterval(FILE *fd, mpfi_t interval);


/* This function performs the differentiation.
   See the commentaries below.
*/
void auto_diff(mpfi_t* res, node *f, mpfi_t x, int n);


int AD(chain **res, void **args) {
  node *f;
  mpfi_t x;
  int i,n;
  mpfi_t *res_array;
  mpfi_t *temp;

  f = (node *)args[0];
  n = *( (int *)args[2] );

  res_array = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
  for(i=0;i<=n;i++) mpfi_init2(res_array[i], getToolPrecision());

  
  mpfi_init2(x, getToolPrecision());
  mpfi_set(x, *( (mpfi_t *)args[1] ));

  auto_diff(res_array, f, x, n);

  *res = NULL;
  for(i=n;i>=0;i--) {
    temp = (mpfi_t *)safeMalloc(sizeof(mpfi_t));
    mpfi_init2(*temp, getToolPrecision());
    mpfi_set(*temp, res_array[i]);
    *res = addElement(*res, temp);
  }
  
  free(res_array);
  mpfi_clear(x);
	   
  return 1;
}


void computeBinomials(mpfi_t **res, int n) {
  int m,k;
  
  mpfi_set_ui(res[0][0], 1);
  for(m=1; m<=n; m++) {
    mpfi_set_ui(res[m][0], 1);
    for(k=1; k<=m-1; k++) {
      mpfi_add(res[m][k], res[m-1][k-1], res[m-1][k]);
    }
    mpfi_set_ui(res[m][m], 1);
  }

  return;
}


void multiplication_AD(mpfi_t *res, mpfi_t *f, mpfi_t *g, int n) {
  int i,j,p;
  mpfi_t temp;
  mpfi_t **binomial_array;

  binomial_array = (mpfi_t **)safeMalloc( (n+1)*sizeof(mpfi_t *));
  for(i=0;i<=n;i++) {
    binomial_array[i] = (mpfi_t *)safeMalloc( (n+1)*sizeof(mpfi_t) );
    for(j=0;j<=n;j++) {
      mpfi_init2(binomial_array[i][j], getToolPrecision());
    }
  }
  computeBinomials(binomial_array, n);

  mpfi_init2(temp, getToolPrecision());

  for(p=0;p<=n;p++) {
    i=0; j=p; mpfi_set_ui(res[p], 0);
    while(i<=p) {
      mpfi_mul(temp, f[i], g[j]);
      mpfi_mul(temp, temp, binomial_array[p][i]);
      mpfi_add(res[p], res[p], temp);

      i++;
      j--;
    }
  }

  mpfi_clear(temp);

  for(i=0;i<=n;i++) {
    for(j=0;j<=n;j++) {
      mpfi_clear(binomial_array[i][j]);
    }
    free(binomial_array[i]);
  }
  free(binomial_array);
  return;
}

void composition_AD(mpfi_t *res, mpfi_t *g, mpfi_t *f, int n) {
  mpfi_t *temp_array;
  int i;

  if(n==0) mpfi_set(res[0], g[0]);
  else {
    temp_array = (mpfi_t *)safeMalloc(n*sizeof(mpfi_t));
    for(i=0;i<=n-1;i++) {
      mpfi_init2(temp_array[i], getToolPrecision());
    }

    composition_AD(temp_array, g+1, f, n-1);
    multiplication_AD(res+1, temp_array, f+1, n-1);
    mpfi_set(res[0], g[0]);

    for(i=0;i<=n-1;i++) mpfi_clear(temp_array[i]);
    free(temp_array);
  }

  return ;
}


/* res is a reserved space for n+1 mpfi_t such that:
   res = [ f(x), f'(x), f''(x), ..., f^(n)(x) ]
*/
void auto_diff(mpfi_t* res, node *f, mpfi_t x, int n) {
  int i;
  mpfi_t *res1, *res2, *res3, *res4, *res5, *res6;
  mpfr_t minusOne;
  node *simplifiedChild1, *simplifiedChild2;
  mpfi_t temp1,temp2;
  switch (f->nodeType) {

  case VARIABLE:
    mpfi_set(res[0], x);
    if(n>=1) {
      mpfi_set_ui(res[1], 1);
      for(i=2;i<=n;i++) mpfi_set_ui(res[i], 0);
    }
    break;

  case PI_CONST:
    mpfi_const_pi(res[0]);
    for(i=1;i<=n;i++) mpfi_set_ui(res[i], 0);
    break;

  case CONSTANT:
    mpfi_set_fr(res[0], *(f->value));
    for(i=1;i<=n;i++) mpfi_set_ui(res[i], 0);
    break;

  case NEG:
    res1 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    for(i=0;i<=n;i++) {
      mpfi_init2(res1[i], getToolPrecision());
    }
    auto_diff(res1, f->child1, x, n);
    for(i=0;i<=n;i++)  mpfi_neg(res[i], res1[i]);

    for(i=0;i<=n;i++) mpfi_clear(res1[i]);
    free(res1);

    break;

  case ADD:
    res1 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    res2 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    for(i=0;i<=n;i++) {
      mpfi_init2(res1[i], getToolPrecision());
      mpfi_init2(res2[i], getToolPrecision());
    }
    auto_diff(res1, f->child1, x, n);
    auto_diff(res2, f->child2, x, n);
    for(i=0;i<=n;i++)  mpfi_add(res[i], res1[i], res2[i]);

    for(i=0;i<=n;i++) {
      mpfi_clear(res1[i]);
      mpfi_clear(res2[i]);
    }
    free(res1);
    free(res2);

    break;

  case SUB:
    res1 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    res2 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    for(i=0;i<=n;i++) {
      mpfi_init2(res1[i], getToolPrecision());
      mpfi_init2(res2[i], getToolPrecision());
    }
    auto_diff(res1, f->child1, x, n);
    auto_diff(res2, f->child2, x, n);
    for(i=0;i<=n;i++)  mpfi_sub(res[i], res1[i], res2[i]);

    for(i=0;i<=n;i++) {
      mpfi_clear(res1[i]);
      mpfi_clear(res2[i]);
    }
    free(res1);
    free(res2);

    break;

  case MUL:
    res1 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    res2 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    for(i=0;i<=n;i++) {
      mpfi_init2(res1[i], getToolPrecision());
      mpfi_init2(res2[i], getToolPrecision());
    }
    auto_diff(res1, f->child1, x, n);
    auto_diff(res2, f->child2, x, n);

    multiplication_AD(res, res1, res2, n);

    for(i=0;i<=n;i++) {
      mpfi_clear(res1[i]);
      mpfi_clear(res2[i]);
    }
    free(res1);
    free(res2);
    break;

  case DIV:
    res1 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    res2 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    res3 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    res4 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    for(i=0;i<=n;i++) {
      mpfi_init2(res1[i], getToolPrecision());
      mpfi_init2(res2[i], getToolPrecision());
      mpfi_init2(res3[i], getToolPrecision());
      mpfi_init2(res4[i], getToolPrecision());
    }

    auto_diff(res1, f->child2, x, n);

    mpfr_init2(minusOne, getToolPrecision());
    
    mpfr_set_si(minusOne, -1, GMP_RNDN);
    constantPower_diff(res2, minusOne, res1[0], n);
    composition_AD(res3, res2, res1, n);

    auto_diff(res4, f->child1, x, n);
    multiplication_AD(res, res3, res4, n);

    mpfr_clear(minusOne);
    for(i=0;i<=n;i++) {
      mpfi_clear(res1[i]);
      mpfi_clear(res2[i]);
      mpfi_clear(res3[i]);
      mpfi_clear(res4[i]);
    }
    free(res1);
    free(res2);
    free(res3);
    free(res4);

    break;

  case SQRT:
  case EXP:
  case LOG:
  case LOG_2:
  case LOG_10:
  case SIN:
  case COS:
  case TAN:
  case ASIN:
  case ACOS:
  case ATAN:
  case SINH:
  case COSH:
  case TANH:
  case ASINH:
  case ACOSH:
  case ATANH:
  case ABS:
  case DOUBLE:
  case DOUBLEDOUBLE:
  case TRIPLEDOUBLE:
  case ERF: 
  case ERFC:
  case LOG_1P:
  case EXP_M1:
  case DOUBLEEXTENDED:
  case CEIL:
  case FLOOR:
  case NEARESTINT:

    res1 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    res2 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    for(i=0;i<=n;i++) {
      mpfi_init2(res1[i], getToolPrecision());
      mpfi_init2(res2[i], getToolPrecision());
    }

    auto_diff(res1, f->child1, x, n);
    baseFunction_diff(res2, f->nodeType, res1[0], n);
    composition_AD(res, res2, res1, n);

    for(i=0;i<=n;i++) {
      mpfi_clear(res1[i]);
      mpfi_clear(res2[i]);
    }
    free(res1);
    free(res2);
    
    break;
  case POW:
    
    if (((f->child2)->nodeType==CONSTANT) && ((f->child1)->nodeType==VARIABLE)){
      constantPower_diff(res,*(f->child2->value) , x, n);
    }
    else{
      simplifiedChild2=simplifyTreeErrorfree(f->child2);
      simplifiedChild1=simplifyTreeErrorfree(f->child1);
      
      if ((simplifiedChild2->nodeType==CONSTANT) &&(simplifiedChild1->nodeType==CONSTANT)) { //we have the ct1^ct2 case
         // printf("We are in the  ct1^ct2 case");       
         mpfi_init2(temp1, getToolPrecision());
         mpfi_set_fr(temp1, *(simplifiedChild1->value));
         mpfi_init2(temp2, getToolPrecision());
         mpfi_set_fr(temp2, *(simplifiedChild2->value));
         mpfi_pow(res[0],temp1,temp2);
         for(i=1;i<=n;i++) mpfi_set_ui(res[i],0);

         mpfi_clear(temp1);
         mpfi_clear(temp2);
      }
      else if (simplifiedChild2->nodeType==CONSTANT) { //we have the f^p case
        //printf("We are in the  f^p case");        
        res1 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
        res2 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
        for(i=0;i<=n;i++) {
          mpfi_init2(res1[i], getToolPrecision());
          mpfi_init2(res2[i], getToolPrecision());
        }
        auto_diff(res1, f->child1, x, n);
        constantPower_diff(res2,*(simplifiedChild2->value) , res1[0], n);
        composition_AD(res, res2, res1, n); 
        for(i=0;i<=n;i++) {
          mpfi_clear(res1[i]);
          mpfi_clear(res2[i]); 
        }
        free(res1);
        free(res2);    
      } 
       else if (simplifiedChild1->nodeType==CONSTANT) { //we have the p^f case
        //printf("We are in the  p^f case");     
        res1 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
        res2 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
        for(i=0;i<=n;i++) {
          mpfi_init2(res1[i], getToolPrecision());
          mpfi_init2(res2[i], getToolPrecision());
        }
        auto_diff(res1, f->child2, x, n);
        powerFunction_diff(res2,*(simplifiedChild1->value) , res1[0], n);
        composition_AD(res, res2, res1, n); 
        for(i=0;i<=n;i++) {
          mpfi_clear(res1[i]);
          mpfi_clear(res2[i]); 
        }
        free(res1);
        free(res2);    
      } 
      else {
      //printf("We are in the  f^g case");     
      res1 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
      res2 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
      res3 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
      res4 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
      res5 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
      res6 = (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
      for(i=0;i<=n;i++) {
        mpfi_init2(res1[i], getToolPrecision());
        mpfi_init2(res2[i], getToolPrecision());
        mpfi_init2(res3[i], getToolPrecision());
        mpfi_init2(res4[i], getToolPrecision());
        mpfi_init2(res5[i], getToolPrecision());
        mpfi_init2(res6[i], getToolPrecision());
      }

      
      auto_diff(res1, f->child1, x, n);
      log_diff(res2,res1[0],n);
      composition_AD(res3, res2, res1, n);
      auto_diff(res4,f->child2,x,n);
      multiplication_AD(res5,res3,res4,n);
      exp_diff(res6,res5[0],n);
      composition_AD(res,res6,res5,n);
      
      for(i=0;i<=n;i++) {
        mpfi_clear(res1[i]);
        mpfi_clear(res2[i]);
        mpfi_clear(res3[i]);
        mpfi_clear(res4[i]);
        mpfi_clear(res5[i]);
        mpfi_clear(res6[i]);
      }
      free(res1);
      free(res2);
      free(res3);
      free(res4);
      free(res5);
      free(res6);

    }
    free_memory(simplifiedChild2);
    free_memory(simplifiedChild1);
  }
    break;

  case LIBRARYFUNCTION:
    break;

  default:
   fprintf(stderr,"Error: AD: unknown identifier (%d) in the tree\n",f->nodeType);
   exit(1);
  }

  return;
}
