/*

Copyright 2009 by 

Laboratoire de l'Informatique du Parallélisme, 
UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668

Contributors Ch. Lauter, S. Chevillard, N. Jourdan

christoph.lauter@ens-lyon.fr
sylvain.chevillard@ens-lyon.fr
nicolas.jourdan@ens-lyon.fr

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

*/


#include "taylorform.h"
#include "external.h"
#include <stdio.h>
#include <stdlib.h>


/***********************************************************/
/***********************************************************/
/***************functions related to autodiff!!!!!!!!!******/
/***********************************************************/
/***********************************************************/
void auto_diff(mpfi_t* res, node *f, mpfi_t x, int n);


/* Computes the successive derivatives of y -> y^p at point x */ 
/* [x^p    p*x^(p-1)   ...   p*(p-1)*...*(p-n+1)*x^(p-n) ]    */
void constantPower_diff(mpfi_t *res, mpfr_t p, mpfi_t x, int n) {
  mpfi_t expo, acc;
  int i;
  
  /* The precision of expo is set in such a way that expo will
     be a point interval during the algorithm */
  mpfi_init2(expo, mpfr_get_prec(p)+8*sizeof(int));
  mpfi_init2(acc, getToolPrecision());
  
  mpfi_set_fr(expo, p);
  mpfi_set_ui(acc, 1);
  
  for(i=0; i<=n; i++) {
    if (mpfi_is_zero(acc)) mpfi_set_ui(res[i],0);
    else {
      mpfi_pow(res[i], x, expo);
      mpfi_mul(res[i], res[i], acc);
      
      mpfi_mul(acc, acc, expo);
      mpfi_sub_ui(expo, expo, 1);
    }
  }

  mpfi_clear(expo);
  mpfi_clear(acc);

  return;
}


void exp_diff(mpfi_t *res, mpfi_t x, int n) {
  int i;
  mpfi_t temp;

  mpfi_init2(temp, getToolPrecision());

  mpfi_exp(temp, x);
  for(i=0;i<=n;i++) mpfi_set(res[i], temp);

  mpfi_clear(temp);
  return;
}

void expm1_diff(mpfi_t *res, mpfi_t x, int n) {
  int i;
  mpfi_t temp;

  mpfi_init2(temp, getToolPrecision());

  mpfi_exp(temp, x);
  for(i=0;i<=n;i++) mpfi_set(res[i], temp);

  mpfi_sub_ui(res[0], res[0], 1);

  mpfi_clear(temp);
  return;
}

void powerFunction_diff(mpfi_t *res, mpfr_t p, mpfi_t x, int n) { //the power function is: p^x, where p is a positive ct
  int i;
  mpfi_t temp1,temp2;

  mpfi_init2(temp1, getToolPrecision());
  mpfi_init2(temp2, getToolPrecision());
  mpfi_set_fr(temp1,p);
  mpfi_log(temp1,temp1);
  mpfi_mul(temp2,temp1,x);
  mpfi_exp(temp2, temp2);
  for(i=0;i<=n;i++) {
    mpfi_set(res[i], temp2);
    mpfi_mul(temp2,temp2,temp1);
  }
  mpfi_clear(temp1);
  mpfi_clear(temp2);
  return;
}



void log_diff(mpfi_t *res, mpfi_t x, int n) {
  mpfr_t minusOne;
  
  mpfi_log(res[0], x);

  if(n>=1) {
    mpfr_init2(minusOne, getToolPrecision());
    mpfr_set_si(minusOne, -1, GMP_RNDN);
    constantPower_diff(res+1, minusOne, x, n-1);
    mpfr_clear(minusOne);
  }
  return;
}

void log2_diff(mpfi_t *res, mpfi_t x, int n) {
  int i;
  mpfi_t log2;
  mpfi_init2(log2, getToolPrecision());

  mpfi_set_ui(log2, 2); mpfi_log(log2, log2);
  log_diff(res,x,n);
  for(i=0;i<=n;i++) mpfi_div(res[i], res[i], log2);

  mpfi_clear(log2);
  return;
}

void log10_diff(mpfi_t *res, mpfi_t x, int n) {
  int i;
  mpfi_t log10;
  mpfi_init2(log10, getToolPrecision());

  mpfi_set_ui(log10, 10); mpfi_log(log10, log10);
  log_diff(res,x,n);
  for(i=0;i<=n;i++) mpfi_div(res[i], res[i], log10);

  mpfi_clear(log10);
  return;
}

void sin_diff(mpfi_t *res, mpfi_t x, int n) {
  int i,index;
  mpfi_t vals[4];
  
  if(n==0) mpfi_sin(res[0], x);
  else {
    
    for(index=0;index<4;index++) mpfi_init2(vals[index], getToolPrecision());
    
    mpfi_sin(vals[0], x);  mpfi_cos(vals[1], x);
    mpfi_neg(vals[2],vals[0]);
    mpfi_neg(vals[3],vals[1]);

    for(i=0;i<=n;i++) mpfi_set(res[i], vals[i % 4]);

    for(index=0;index<4;index++) mpfi_clear(vals[index]);
  }

  return;
}

void cos_diff(mpfi_t *res, mpfi_t x, int n) {
  int i,index;
  mpfi_t vals[4];
  
  if(n==0) mpfi_cos(res[0], x);
  else {
    
    for(index=0;index<4;index++) mpfi_init2(vals[index], getToolPrecision());
    
    mpfi_cos(vals[0], x);  mpfi_sin(vals[3], x);
    mpfi_neg(vals[2],vals[0]);
    mpfi_neg(vals[1],vals[3]);

    for(i=0;i<=n;i++) mpfi_set(res[i], vals[i % 4]);

    for(index=0;index<4;index++) mpfi_clear(vals[index]);
  }

  return;
}

void sinh_diff(mpfi_t *res, mpfi_t x, int n) {
  int i,index;
  mpfi_t vals[2];
  
  if(n==0) mpfi_sinh(res[0], x);
  else {
    
    for(index=0;index<2;index++) mpfi_init2(vals[index], getToolPrecision());
    
    mpfi_sinh(vals[0], x);  mpfi_cosh(vals[1], x);

    for(i=0;i<=n;i++) mpfi_set(res[i], vals[i % 2]);

    for(index=0;index<2;index++) mpfi_clear(vals[index]);
  }

  return;
}

void cosh_diff(mpfi_t *res, mpfi_t x, int n) {
  int i,index;
  mpfi_t vals[2];
  
  if(n==0) mpfi_cosh(res[0], x);
  else {
    
    for(index=0;index<2;index++) mpfi_init2(vals[index], getToolPrecision());
    
    mpfi_cosh(vals[0], x);  mpfi_sinh(vals[1], x);

    for(i=0;i<=n;i++) mpfi_set(res[i], vals[i % 2]);

    for(index=0;index<2;index++) mpfi_clear(vals[index]);
  }

  return;
}

/*tan_diff : recurrence formula: p_(n+1)(u)=p'_n(u) * (1+u^2) , u=tan(x)
tan^(n)=p_(n)(u)
p_0=u;
p_1=1+u^2
*/
void tan_diff(mpfi_t *res, mpfi_t x, int n) {
  int i,index;
  mpfi_t *coeffs_array, *powers_array, u, partialSum,s1,oldCoeff,pow;
  
    
  coeffs_array = (mpfi_t *)safeMalloc( (n+2)*sizeof(mpfi_t));
  powers_array = (mpfi_t *)safeMalloc( (n+2)*sizeof(mpfi_t));

  for (index=0; index<=n+1; index++){
    mpfi_init2(coeffs_array[index], getToolPrecision());
    mpfi_init2(powers_array[index], getToolPrecision());
  }
 

  mpfi_init2(partialSum,getToolPrecision());
  mpfi_init2(oldCoeff,getToolPrecision());
  mpfi_init2(s1,getToolPrecision());
  mpfi_init2(pow,getToolPrecision());
  mpfi_init2(u,getToolPrecision());
  mpfi_tan(u, x);
  
  for (index=0; index<=n; index++){
    if (index==0){
      mpfi_set_ui(coeffs_array[0],0);
      mpfi_set_ui(coeffs_array[1],1);
      mpfi_set_ui(powers_array[0],1);
      mpfi_set(powers_array[1],u);
      
    }
    else if (index==1){
      mpfi_set_ui(coeffs_array[0],1);
      mpfi_set_ui(coeffs_array[1],0);
      mpfi_set_ui(coeffs_array[2],1);
      mpfi_set_ui(pow,2);
      mpfi_pow(powers_array[2],u,pow);
    }
    else{
     //newcoeff(k)=(k-1) coeff(k-1)+(k+1)coeff(k+1) except for the last two ones and first two ones
     //the new coefficent added at index+1
      mpfi_set(coeffs_array[index+1],coeffs_array[index]);
      mpfi_mul_ui(coeffs_array[index+1],coeffs_array[index+1],index);
      //keep the coeff at index      
      mpfi_set(oldCoeff,coeffs_array[index]);
      //modify at index with new value      
      mpfi_set(coeffs_array[index],coeffs_array[index-1]);
      mpfi_mul_ui(coeffs_array[index],coeffs_array[index],index-1);
      //inter-extremities case
      for (i=index-1; i>1;i--){
        mpfi_set(s1,oldCoeff);
        mpfi_mul_ui(s1,s1,i+1);
        mpfi_set(oldCoeff,coeffs_array[i]);
        mpfi_set(coeffs_array[i],coeffs_array[i-1]);
        mpfi_mul_ui(coeffs_array[i],coeffs_array[i],i-1);
        mpfi_add(coeffs_array[i],coeffs_array[i], s1);
      }
      //modify at coeff 0 with new value      
      mpfi_set(coeffs_array[0],coeffs_array[1]);
      //modify at 1 with new value      
      mpfi_set(coeffs_array[1],oldCoeff);
      mpfi_mul_ui(coeffs_array[1],coeffs_array[1],2);
      //add in powers array one more power
      mpfi_set_ui(pow,index+1);
      mpfi_pow(powers_array[index+1],u,pow);
    }
    //we have the coeffcients, we compute the result
    /*printf("index=%d",index);
    for (i=0;i<=index+1; i++){
      printInterval(coeffs_array[i]);
    }
    printf("\n");*/
    mpfi_set_ui(res[index],0);    
    for (i=0;i<=index+1; i++){
      mpfi_mul(partialSum,coeffs_array[i],powers_array[i]);
      mpfi_add(res[index],res[index],partialSum);
    }

  }
  for (index=0; index<=n+1; index++){
    mpfi_clear(coeffs_array[index]);
    mpfi_clear(powers_array[index]);
  }
  mpfi_clear(oldCoeff);
  mpfi_clear(partialSum);
  mpfi_clear(s1);
  mpfi_clear(pow);
  
  return;
}

/*tanh_diff : recurrence formula: p_(n+1)(u)=p'_n(u) * (1-u^2) , u=tanh(x)
tan^(n)=p_(n)(u)
p_0=u;
p_1=1-u^2
*/
void tanh_diff(mpfi_t *res, mpfi_t x, int n) {
  int i,index;
  mpfi_t *coeffs_array, *powers_array, u, partialSum,s1,oldCoeff,pow;
  
    
  coeffs_array = (mpfi_t *)safeMalloc( (n+2)*sizeof(mpfi_t));
  powers_array = (mpfi_t *)safeMalloc( (n+2)*sizeof(mpfi_t));

  for (index=0; index<=n+1; index++){
    mpfi_init2(coeffs_array[index], getToolPrecision());
    mpfi_init2(powers_array[index], getToolPrecision());
  }
 

  mpfi_init2(partialSum,getToolPrecision());
  mpfi_init2(oldCoeff,getToolPrecision());
  mpfi_init2(s1,getToolPrecision());
  mpfi_init2(pow,getToolPrecision());
  mpfi_init2(u,getToolPrecision());
  mpfi_tanh(u, x);
  
  for (index=0; index<=n; index++){
    if (index==0){
      mpfi_set_ui(coeffs_array[0],0);
      mpfi_set_ui(coeffs_array[1],1);
      mpfi_set_ui(powers_array[0],1);
      mpfi_set(powers_array[1],u);
      
    }
    else if (index==1){
      mpfi_set_ui(coeffs_array[0],1);
      mpfi_set_ui(coeffs_array[1],0);
      mpfi_set_si(coeffs_array[2],-1);
      mpfi_set_ui(pow,2);
      mpfi_pow(powers_array[2],u,pow);
    }
    else{
     //newcoeff(k)=-(k-1) coeff(k-1)+(k+1)coeff(k+1) except for the last two ones and first two ones
     //the new coefficent added at index+1
      mpfi_set(coeffs_array[index+1],coeffs_array[index]);
      mpfi_mul_si(coeffs_array[index+1],coeffs_array[index+1],-index);
      //keep the coeff at index      
      mpfi_set(oldCoeff,coeffs_array[index]);
      //modify at index with new value      
      mpfi_set(coeffs_array[index],coeffs_array[index-1]);
      mpfi_mul_si(coeffs_array[index],coeffs_array[index],-(index-1));
      //inter-extremities case
      for (i=index-1; i>1;i--){
        mpfi_set(s1,oldCoeff);
        mpfi_mul_ui(s1,s1,i+1);
        mpfi_set(oldCoeff,coeffs_array[i]);
        mpfi_set(coeffs_array[i],coeffs_array[i-1]);
        mpfi_mul_si(coeffs_array[i],coeffs_array[i],-(i-1));
        mpfi_add(coeffs_array[i],coeffs_array[i], s1);
      }
      //modify at coeff 0 with new value      
      mpfi_set(coeffs_array[0],coeffs_array[1]);
      //modify at 1 with new value      
      mpfi_set(coeffs_array[1],oldCoeff);
      mpfi_mul_ui(coeffs_array[1],coeffs_array[1],2);
      //add in powers array one more power
      mpfi_set_ui(pow,index+1);
      mpfi_pow(powers_array[index+1],u,pow);
    }
    //we have the coeffcients, we compute the result
    /*printf("index=%d",index);
    for (i=0;i<=index+1; i++){
      printInterval(coeffs_array[i]);
    }
    printf("\n");*/
    mpfi_set_ui(res[index],0);    
    for (i=0;i<=index+1; i++){
      mpfi_mul(partialSum,coeffs_array[i],powers_array[i]);
      mpfi_add(res[index],res[index],partialSum);
    }

  }
  for (index=0; index<=n+1; index++){
    mpfi_clear(coeffs_array[index]);
    mpfi_clear(powers_array[index]);
  }
  mpfi_clear(oldCoeff);
  mpfi_clear(partialSum);
  mpfi_clear(s1);
  mpfi_clear(pow);
  return;
}



/*atan_diff : reccurence formula: p_(n+1)=p'_n * (1+x^2) -2*n *x * p_n 
atan^(n)=p_(n)/((1+x^2)^n)
p_1=1;
*/
void atan_diff(mpfi_t *res, mpfi_t x, int n) {
  int i,index;
  mpfi_t *coeffs_array, *powers_array, u, partialSum,s1,oldCoeff,pow,a1,nominator;
  
    
  coeffs_array = (mpfi_t *)safeMalloc( (n)*sizeof(mpfi_t));
  powers_array = (mpfi_t *)safeMalloc( (n)*sizeof(mpfi_t));

  for (index=0; index<=n-1; index++){
    mpfi_init2(coeffs_array[index], getToolPrecision());
    mpfi_init2(powers_array[index], getToolPrecision());
  }
 

  mpfi_init2(partialSum,getToolPrecision());
  mpfi_init2(oldCoeff,getToolPrecision());
  mpfi_init2(s1,getToolPrecision());
  mpfi_init2(pow,getToolPrecision());
  mpfi_init2(u,getToolPrecision());
  mpfi_init2(a1,getToolPrecision());
  mpfi_atan(u, x);
  mpfi_init2(nominator,getToolPrecision());
  mpfi_sqr(nominator,x);
  mpfi_add_si(nominator, nominator,1);
  mpfi_set(res[0],u);
  //if (n==0) 
  //else
  //{
    for (index=0; index<=n-1; index++){
      if (index==0){
        mpfi_set_ui(coeffs_array[0],1);
        mpfi_set_ui(powers_array[0],1);
      }
      else if (index==1){
        mpfi_set_ui(coeffs_array[0],0);
        mpfi_set_si(coeffs_array[1],-2);
        mpfi_set_ui(pow,index);
        mpfi_pow(powers_array[1],x,pow);
      }
      else if (index==2){
        mpfi_set_si(coeffs_array[0],-2);
        mpfi_set_ui(coeffs_array[1],0);
        mpfi_set_ui(coeffs_array[2],6);
        mpfi_set_ui(pow,index);
        mpfi_pow(powers_array[2],x,pow);
      }
      else if (index==3){
        mpfi_set_si(coeffs_array[0],0);
        mpfi_set_ui(coeffs_array[1],24);
        mpfi_set_ui(coeffs_array[2],0);
        mpfi_set_si(coeffs_array[3],-24);
        mpfi_set_ui(pow,index);
        mpfi_pow(powers_array[3],x,pow);
      }
      else{
        //newcoeff(k)=(k-1 -2*n) coeff(k-1)+(k+1)coeff(k+1) except for the last two ones and first two ones
        //the new coefficent added at index
        mpfi_set(coeffs_array[index],coeffs_array[index-1]);
        mpfi_mul_si(coeffs_array[index],coeffs_array[index],index-1-2*(index));
        //keep the coeff at index-1     
        mpfi_set(oldCoeff,coeffs_array[index-1]);
        //modify at index with new value      
        mpfi_set(coeffs_array[index-1],coeffs_array[index-2]);
        mpfi_mul_si(coeffs_array[index-1],coeffs_array[index-1],index-2-2*(index));
        //inter-extremities case
        for (i=index-2; i>1;i--){
          mpfi_set(s1,oldCoeff);
          mpfi_mul_si(s1,s1,i+1);
          mpfi_set(oldCoeff,coeffs_array[i]);
          mpfi_set(coeffs_array[i],coeffs_array[i-1]);
          mpfi_mul_si(coeffs_array[i],coeffs_array[i],i-1-2*(index));
          mpfi_add(coeffs_array[i],coeffs_array[i], s1);
        } 
        //for a1 and ao special case
        //newa1=2*a2-2n*a0   newa0=a1      
        //compute 2*a2      
        mpfi_mul_si(oldCoeff, oldCoeff,2);      
        mpfi_set(a1, coeffs_array[1]);      
        //compute -2n*a0 
        mpfi_mul_si(coeffs_array[1],coeffs_array[0],-2*(index));
        //add
        mpfi_add(coeffs_array[1], coeffs_array[1],oldCoeff);
        //modify at coeff 0 with new value      
        mpfi_set(coeffs_array[0],a1);
        //add in powers array one more power
        mpfi_set_ui(pow,index);
        mpfi_pow(powers_array[index],x,pow);
      }
      //we have the coeffcients, we compute the result
     /* printf("index=%d",index);
      for (i=0;i<=index; i++){
        printInterval(coeffs_array[i]);
      }
      printf("\n");*/
      mpfi_set_ui(res[index+1],0);    
      for (i=0;i<=index; i++){
        mpfi_mul(partialSum,coeffs_array[i],powers_array[i]);
        mpfi_add(res[index+1],res[index+1],partialSum);
      } 
      mpfi_sqr(nominator,x);
      mpfi_add_si(nominator, nominator,1);
  
      mpfi_set_ui(pow,index+1);
      mpfi_pow(nominator,nominator,pow);
      mpfi_div(res[index+1],res[index+1],nominator);
    }
  //}
  for (index=0; index<=n-1; index++){
    mpfi_clear(coeffs_array[index]);
    mpfi_clear(powers_array[index]);
  }
  mpfi_clear(oldCoeff);
  mpfi_clear(partialSum);
  mpfi_clear(s1);
  mpfi_clear(pow);
  mpfi_clear(a1);
  mpfi_clear(nominator);
  return;
}


/*atanh_diff : reccurence formula: p_(n+1)=p'_n * (1-x^2) +2*n *x * p_n 
atan^(n)=p_(n)/((1-x^2)^n)
p_1=1;
*/
void atanh_diff(mpfi_t *res, mpfi_t x, int n) {
  int i,index;
  mpfi_t *coeffs_array, *powers_array, u, partialSum,s1,oldCoeff,pow,a1,nominator;
  
    
  coeffs_array = (mpfi_t *)safeMalloc( (n)*sizeof(mpfi_t));
  powers_array = (mpfi_t *)safeMalloc( (n)*sizeof(mpfi_t));

  for (index=0; index<=n-1; index++){
    mpfi_init2(coeffs_array[index], getToolPrecision());
    mpfi_init2(powers_array[index], getToolPrecision());
  }
 

  mpfi_init2(partialSum,getToolPrecision());
  mpfi_init2(oldCoeff,getToolPrecision());
  mpfi_init2(s1,getToolPrecision());
  mpfi_init2(pow,getToolPrecision());
  mpfi_init2(u,getToolPrecision());
  mpfi_init2(a1,getToolPrecision());
  mpfi_atanh(u, x);
  mpfi_init2(nominator,getToolPrecision());
  mpfi_sqr(nominator,x);
  mpfi_add_si(nominator, nominator,1);
  mpfi_set(res[0],u);
  //if (n==0) 
  //else
  //{
    for (index=0; index<=n-1; index++){
      if (index==0){
        mpfi_set_ui(coeffs_array[0],1);
        mpfi_set_ui(powers_array[0],1);
      }
      else if (index==1){
        mpfi_set_ui(coeffs_array[0],0);
        mpfi_set_si(coeffs_array[1],2);
        mpfi_set_ui(pow,index);
        mpfi_pow(powers_array[1],x,pow);
      }
      else if (index==2){
        mpfi_set_si(coeffs_array[0],2);
        mpfi_set_ui(coeffs_array[1],0);
        mpfi_set_ui(coeffs_array[2],6);
        mpfi_set_ui(pow,index);
        mpfi_pow(powers_array[2],x,pow);
      }
      else if (index==3){
        mpfi_set_si(coeffs_array[0],0);
        mpfi_set_ui(coeffs_array[1],24);
        mpfi_set_ui(coeffs_array[2],0);
        mpfi_set_si(coeffs_array[3],24);
        mpfi_set_ui(pow,index);
        mpfi_pow(powers_array[3],x,pow);
      }
      else{
        //newcoeff(k)=(-(k-1) +2*n) coeff(k-1)+(k+1)coeff(k+1) except for the last two ones and first two ones
        //the new coefficent added at index
        mpfi_set(coeffs_array[index],coeffs_array[index-1]);
        mpfi_mul_si(coeffs_array[index],coeffs_array[index],-(index-1)+2*(index));
        //keep the coeff at index-1     
        mpfi_set(oldCoeff,coeffs_array[index-1]);
        //modify at index with new value      
        mpfi_set(coeffs_array[index-1],coeffs_array[index-2]);
        mpfi_mul_si(coeffs_array[index-1],coeffs_array[index-1],-(index-2)+2*(index));
        //inter-extremities case
        for (i=index-2; i>1;i--){
          mpfi_set(s1,oldCoeff);
          mpfi_mul_si(s1,s1,i+1);
          mpfi_set(oldCoeff,coeffs_array[i]);
          mpfi_set(coeffs_array[i],coeffs_array[i-1]);
          mpfi_mul_si(coeffs_array[i],coeffs_array[i],-(i-1)+2*(index));
          mpfi_add(coeffs_array[i],coeffs_array[i], s1);
        } 
        //for a1 and ao special case
        //newa1=2*a2+2n*a0   newa0=a1      
        //compute 2*a2      
        mpfi_mul_si(oldCoeff, oldCoeff,2);      
        mpfi_set(a1, coeffs_array[1]);      
        //compute -2n*a0 
        mpfi_mul_si(coeffs_array[1],coeffs_array[0],2*(index));
        //add
        mpfi_add(coeffs_array[1], coeffs_array[1],oldCoeff);
        //modify at coeff 0 with new value      
        mpfi_set(coeffs_array[0],a1);
        //add in powers array one more power
        mpfi_set_ui(pow,index);
        mpfi_pow(powers_array[index],x,pow);
      }
      //we have the coeffcients, we compute the result
     /* printf("index=%d",index);
      for (i=0;i<=index; i++){
        printInterval(coeffs_array[i]);
      }
      printf("\n");*/
      mpfi_set_ui(res[index+1],0);    
      for (i=0;i<=index; i++){
        mpfi_mul(partialSum,coeffs_array[i],powers_array[i]);
        mpfi_add(res[index+1],res[index+1],partialSum);
      } 
      mpfi_sqr(nominator,x);
      mpfi_si_sub(nominator, 1,nominator);
  
      mpfi_set_ui(pow,index+1);
      mpfi_pow(nominator,nominator,pow);
      mpfi_div(res[index+1],res[index+1],nominator);
    }
  //}
  for (index=0; index<=n-1; index++){
    mpfi_clear(coeffs_array[index]);
    mpfi_clear(powers_array[index]);
  }
  mpfi_clear(oldCoeff);
  mpfi_clear(partialSum);
  mpfi_clear(s1);
  mpfi_clear(pow);
  mpfi_clear(a1);
  mpfi_clear(nominator);
  return;
}






/*asin_diff : reccurence formula: p_(n+1)=p'_n * (1-x^2) +(2*n-1) *x * p_n 
asin^(n)=p_(n)/((1-x^2)^((2*n-1)/2))
p_1=1;
p_2=x;
p_3=2x^2+1
*/
void asin_diff(mpfi_t *res, mpfi_t x, int n) {
  int i,index;
  mpfi_t *coeffs_array, *powers_array, u, partialSum,s1,oldCoeff,pow,a1,nominator;
  
  //the polynomial for the nth derivative has degree n-1, need n coeffs  
  coeffs_array = (mpfi_t *)safeMalloc( (n)*sizeof(mpfi_t));
  powers_array = (mpfi_t *)safeMalloc( (n)*sizeof(mpfi_t));

  for (index=0; index<=n-1; index++){
    mpfi_init2(coeffs_array[index], getToolPrecision());
    mpfi_init2(powers_array[index], getToolPrecision());
  }
 

  mpfi_init2(partialSum,getToolPrecision());
  mpfi_init2(oldCoeff,getToolPrecision());
  mpfi_init2(s1,getToolPrecision());
  mpfi_init2(pow,getToolPrecision());
  mpfi_init2(u,getToolPrecision());
  mpfi_init2(a1,getToolPrecision());
  mpfi_asin(u, x);
  //we need the nominator for dividing the polynominal by ((1-x^2)^((2*n-1)/2))
  mpfi_init2(nominator,getToolPrecision());
  mpfi_sqr(nominator,x);
  mpfi_si_sub(nominator, 1,nominator);
  //put the asin value in res[0]
  mpfi_set(res[0],u);
  //if (n==0) 
  //else
  //{
    for (index=0; index<=n-1; index++){
      if (index==0){
        mpfi_set_ui(coeffs_array[0],1);
        mpfi_set_ui(powers_array[0],1);
      }
      else if (index==1){
        mpfi_set_ui(coeffs_array[0],0);
        mpfi_set_si(coeffs_array[1],1);
        mpfi_set_ui(pow,index);
        mpfi_pow(powers_array[1],x,pow);
      }
      else if (index==2){
        mpfi_set_si(coeffs_array[0],1);
        mpfi_set_ui(coeffs_array[1],0);
        mpfi_set_ui(coeffs_array[2],2);
        mpfi_set_ui(pow,index);
        mpfi_pow(powers_array[2],x,pow);
      }
      else if (index==3){
        mpfi_set_si(coeffs_array[0],0);
        mpfi_set_ui(coeffs_array[1],9);
        mpfi_set_ui(coeffs_array[2],0);
        mpfi_set_si(coeffs_array[3],6);
        mpfi_set_ui(pow,index);
        mpfi_pow(powers_array[3],x,pow);
      }
      else{
        //newcoeff(k)=(-(k-1) +(2*n-1)) coeff(k-1)+(k+1)coeff(k+1) except for the last two ones and first two ones
        //the new coefficent added at index
        mpfi_set(coeffs_array[index],coeffs_array[index-1]);
        mpfi_mul_si(coeffs_array[index],coeffs_array[index],-(index-1)+(2*index-1));
        //keep the coeff at index-1     
        mpfi_set(oldCoeff,coeffs_array[index-1]);
        //modify at index with new value      
        mpfi_set(coeffs_array[index-1],coeffs_array[index-2]);
        mpfi_mul_si(coeffs_array[index-1],coeffs_array[index-1],-(index-2)+(2*index-1));
        //inter-extremities case
        for (i=index-2; i>1;i--){
          mpfi_set(s1,oldCoeff);
          mpfi_mul_si(s1,s1,i+1);
          mpfi_set(oldCoeff,coeffs_array[i]);
          mpfi_set(coeffs_array[i],coeffs_array[i-1]);
          mpfi_mul_si(coeffs_array[i],coeffs_array[i],-(i-1)+(2*index-1));
          mpfi_add(coeffs_array[i],coeffs_array[i], s1);
        } 
        //for a1 and ao special case
        //newa1=2*a2 +(2n-1)*a0   newa0=a1      
        //compute 2*a2      
        mpfi_mul_si(oldCoeff, oldCoeff,2);      
        mpfi_set(a1, coeffs_array[1]);      
        //compute (2n-1)*a0 
        mpfi_mul_si(coeffs_array[1],coeffs_array[0],2*index-1);
        //add
        mpfi_add(coeffs_array[1], coeffs_array[1],oldCoeff);
        //modify at coeff 0 with new value      
        mpfi_set(coeffs_array[0],a1);
        //add in powers array one more power
        mpfi_set_ui(pow,index);
        mpfi_pow(powers_array[index],x,pow);
      }
      //we have the coeffcients, we compute the result
      /*printf("index=%d",index);
      for (i=0;i<=index; i++){
        printInterval(coeffs_array[i]);
      }
      printf("\n");*/
      mpfi_set_ui(res[index+1],0);    
      for (i=0;i<=index; i++){
        mpfi_mul(partialSum,coeffs_array[i],powers_array[i]);
        mpfi_add(res[index+1],res[index+1],partialSum);
      } 
      mpfi_sqr(nominator,x);
      mpfi_si_sub(nominator, 1,nominator);
  
      mpfi_set_si(pow,(2*index+1));
      mpfi_div_si(pow,pow,2);
      mpfi_pow(nominator,nominator,pow);
      mpfi_div(res[index+1],res[index+1],nominator);
    }
  //}
  for (index=0; index<=n-1; index++){
    mpfi_clear(coeffs_array[index]);
    mpfi_clear(powers_array[index]);
  }
  mpfi_clear(oldCoeff);
  mpfi_clear(partialSum);
  mpfi_clear(s1);
  mpfi_clear(pow);
  mpfi_clear(a1);
  mpfi_clear(nominator);
  return;
}

/*acos_diff : except for the res[0], all the terms are equal to -asin^(n)(x)
*/
void acos_diff(mpfi_t *res, mpfi_t x, int n) {
  int i;
 
  
  asin_diff(res,x,n);
  mpfi_acos(res[0],x);
  for (i=1; i<=n;i++)
    mpfi_mul_si(res[i],res[i],-1);
  return;
}

/*asinh_diff : reccurence formula: p_(n+1)=p'_n * (1+x^2) -(2*n-1) *x * p_n 
asin^(n)=p_(n)/((1+x^2)^((2*n-1)/2))
p_1=1;
p_2=-x;
p_3=2x^2-1
*/
void asinh_diff(mpfi_t *res, mpfi_t x, int n) {
  int i,index;
  mpfi_t *coeffs_array, *powers_array, u, partialSum,s1,oldCoeff,pow,a1,nominator;
  
  //the polynomial for the nth derivative has degree n-1, need n coeffs  
  coeffs_array = (mpfi_t *)safeMalloc( (n)*sizeof(mpfi_t));
  powers_array = (mpfi_t *)safeMalloc( (n)*sizeof(mpfi_t));

  for (index=0; index<=n-1; index++){
    mpfi_init2(coeffs_array[index], getToolPrecision());
    mpfi_init2(powers_array[index], getToolPrecision());
  }
 

  mpfi_init2(partialSum,getToolPrecision());
  mpfi_init2(oldCoeff,getToolPrecision());
  mpfi_init2(s1,getToolPrecision());
  mpfi_init2(pow,getToolPrecision());
  mpfi_init2(u,getToolPrecision());
  mpfi_init2(a1,getToolPrecision());
  mpfi_asinh(u, x);
  //we need the nominator for dividing the polynominal by ((1+x^2)^((2*n-1)/2))
  mpfi_init2(nominator,getToolPrecision());
  mpfi_sqr(nominator,x);
  mpfi_add_si(nominator, nominator,1);
  //put the asin value in res[0]
  mpfi_set(res[0],u);
  //if (n==0) 
  //else
  //{
    for (index=0; index<=n-1; index++){
      if (index==0){
        mpfi_set_ui(coeffs_array[0],1);
        mpfi_set_ui(powers_array[0],1);
      }
      else if (index==1){
        mpfi_set_ui(coeffs_array[0],0);
        mpfi_set_si(coeffs_array[1],-1);
        mpfi_set_ui(pow,index);
        mpfi_pow(powers_array[1],x,pow);
      }
      else if (index==2){
        mpfi_set_si(coeffs_array[0],-1);
        mpfi_set_ui(coeffs_array[1],0);
        mpfi_set_ui(coeffs_array[2],2);
        mpfi_set_ui(pow,index);
        mpfi_pow(powers_array[2],x,pow);
      }
      else if (index==3){
        mpfi_set_si(coeffs_array[0],0);
        mpfi_set_ui(coeffs_array[1],9);
        mpfi_set_ui(coeffs_array[2],0);
        mpfi_set_si(coeffs_array[3],-6);
        mpfi_set_ui(pow,index);
        mpfi_pow(powers_array[3],x,pow);
      }
      else{
        //newcoeff(k)=(-(k-1) -(2*n-1)) coeff(k-1)+(k+1)coeff(k+1) except for the last two ones and first two ones
        //the new coefficent added at index
        mpfi_set(coeffs_array[index],coeffs_array[index-1]);
        mpfi_mul_si(coeffs_array[index],coeffs_array[index],(index-1)-(2*index-1));
        //keep the coeff at index-1     
        mpfi_set(oldCoeff,coeffs_array[index-1]);
        //modify at index with new value      
        mpfi_set(coeffs_array[index-1],coeffs_array[index-2]);
        mpfi_mul_si(coeffs_array[index-1],coeffs_array[index-1],(index-2)-(2*index-1));
        //inter-extremities case
        for (i=index-2; i>1;i--){
          mpfi_set(s1,oldCoeff);
          mpfi_mul_si(s1,s1,i+1);
          mpfi_set(oldCoeff,coeffs_array[i]);
          mpfi_set(coeffs_array[i],coeffs_array[i-1]);
          mpfi_mul_si(coeffs_array[i],coeffs_array[i],(i-1)-(2*index-1));
          mpfi_add(coeffs_array[i],coeffs_array[i], s1);
        } 
        //for a1 and ao special case
        //newa1=2*a2 -(2n-1)*a0   newa0=a1      
        //compute 2*a2      
        mpfi_mul_si(oldCoeff, oldCoeff,2);      
        mpfi_set(a1, coeffs_array[1]);      
        //compute -(2n-1)*a0 
        mpfi_mul_si(coeffs_array[1],coeffs_array[0],-(2*index-1));
        //add
        mpfi_add(coeffs_array[1], coeffs_array[1],oldCoeff);
        //modify at coeff 0 with new value      
        mpfi_set(coeffs_array[0],a1);
        //add in powers array one more power
        mpfi_set_ui(pow,index);
        mpfi_pow(powers_array[index],x,pow);
      }
      //we have the coeffcients, we compute the result
      /*printf("index=%d",index);
      for (i=0;i<=index; i++){
        printInterval(coeffs_array[i]);
      }
      printf("\n");*/
      mpfi_set_ui(res[index+1],0);    
      for (i=0;i<=index; i++){
        mpfi_mul(partialSum,coeffs_array[i],powers_array[i]);
        mpfi_add(res[index+1],res[index+1],partialSum);
      } 
      mpfi_sqr(nominator,x);
      mpfi_add_si(nominator, nominator,1);
  
      mpfi_set_si(pow,(2*index+1));
      mpfi_div_si(pow,pow,2);
      mpfi_pow(nominator,nominator,pow);
      mpfi_div(res[index+1],res[index+1],nominator);
    }
  //}
  for (index=0; index<=n-1; index++){
    mpfi_clear(coeffs_array[index]);
    mpfi_clear(powers_array[index]);
  }
  mpfi_clear(oldCoeff);
  mpfi_clear(partialSum);
  mpfi_clear(s1);
  mpfi_clear(pow);
  mpfi_clear(a1);
  mpfi_clear(nominator);
  return;
}


/*acosh_diff : reccurence formula: p_(n+1)=p'_n * (x^2-1) -(2*n-1) *x * p_n 
asin^(n)=p_(n)/((x^2-1)^((2*n-1)/2))
p_1=1;
p_2=-x;
p_3=2x^2+1
*/
void acosh_diff(mpfi_t *res, mpfi_t x, int n) {
  int i,index;
  mpfi_t *coeffs_array, *powers_array, u, partialSum,s1,oldCoeff,pow,a1,nominator;
  
  //the polynomial for the nth derivative has degree n-1, need n coeffs  
  coeffs_array = (mpfi_t *)safeMalloc( (n)*sizeof(mpfi_t));
  powers_array = (mpfi_t *)safeMalloc( (n)*sizeof(mpfi_t));

  for (index=0; index<=n-1; index++){
    mpfi_init2(coeffs_array[index], getToolPrecision());
    mpfi_init2(powers_array[index], getToolPrecision());
  }
 

  mpfi_init2(partialSum,getToolPrecision());
  mpfi_init2(oldCoeff,getToolPrecision());
  mpfi_init2(s1,getToolPrecision());
  mpfi_init2(pow,getToolPrecision());
  mpfi_init2(u,getToolPrecision());
  mpfi_init2(a1,getToolPrecision());
  mpfi_acosh(u, x);
  //we need the nominator for dividing the polynominal by ((1+x^2)^((2*n-1)/2))
  mpfi_init2(nominator,getToolPrecision());
  mpfi_sqr(nominator,x);
  mpfi_add_si(nominator, nominator,1);
  //put the asin value in res[0]
  mpfi_set(res[0],u);
  //if (n==0) 
  //else
  //{
    for (index=0; index<=n-1; index++){
      if (index==0){
        mpfi_set_ui(coeffs_array[0],1);
        mpfi_set_ui(powers_array[0],1);
      }
      else if (index==1){
        mpfi_set_ui(coeffs_array[0],0);
        mpfi_set_si(coeffs_array[1],-1);
        mpfi_set_ui(pow,index);
        mpfi_pow(powers_array[1],x,pow);
      }
      else if (index==2){
        mpfi_set_si(coeffs_array[0],1);
        mpfi_set_ui(coeffs_array[1],0);
        mpfi_set_ui(coeffs_array[2],2);
        mpfi_set_ui(pow,index);
        mpfi_pow(powers_array[2],x,pow);
      }
      else if (index==3){
        mpfi_set_si(coeffs_array[0],0);
        mpfi_set_si(coeffs_array[1],-9);
        mpfi_set_ui(coeffs_array[2],0);
        mpfi_set_si(coeffs_array[3],-6);
        mpfi_set_ui(pow,index);
        mpfi_pow(powers_array[3],x,pow);
      }
      else{
        //newcoeff(k)=((k-1) -(2*n-1)) coeff(k-1)-(k+1)coeff(k+1) except for the last two ones and first two ones
        //the new coefficent added at index
        mpfi_set(coeffs_array[index],coeffs_array[index-1]);
        mpfi_mul_si(coeffs_array[index],coeffs_array[index],(index-1)-(2*index-1));
        //keep the coeff at index-1     
        mpfi_set(oldCoeff,coeffs_array[index-1]);
        //modify at index with new value      
        mpfi_set(coeffs_array[index-1],coeffs_array[index-2]);
        mpfi_mul_si(coeffs_array[index-1],coeffs_array[index-1],(index-2)-(2*index-1));
        //inter-extremities case
        for (i=index-2; i>1;i--){
          mpfi_set(s1,oldCoeff);
          mpfi_mul_si(s1,s1,-(i+1));
          mpfi_set(oldCoeff,coeffs_array[i]);
          mpfi_set(coeffs_array[i],coeffs_array[i-1]);
          mpfi_mul_si(coeffs_array[i],coeffs_array[i],(i-1)-(2*index-1));
          mpfi_add(coeffs_array[i],coeffs_array[i], s1);
        } 
        //for a1 and ao special case
        //newa1=-2*a2 -(2n-1)*a0   newa0=-a1      
        //compute -2*a2      
        mpfi_mul_si(oldCoeff, oldCoeff,-2);      
        mpfi_set(a1, coeffs_array[1]);      
        //compute -(2n-1)*a0 
        mpfi_mul_si(coeffs_array[1],coeffs_array[0],-(2*index-1));
        //add
        mpfi_add(coeffs_array[1], coeffs_array[1],oldCoeff);
        //modify at coeff 0 with new value      
        mpfi_mul_si(coeffs_array[0],a1,-1);
        //add in powers array one more power
        mpfi_set_ui(pow,index);
        mpfi_pow(powers_array[index],x,pow);
      }
      //we have the coeffcients, we compute the result
      /*printf("index=%d",index);
      for (i=0;i<=index; i++){
        printInterval(coeffs_array[i]);
      }
      printf("\n");*/
      mpfi_set_ui(res[index+1],0);    
      for (i=0;i<=index; i++){
        mpfi_mul(partialSum,coeffs_array[i],powers_array[i]);
        mpfi_add(res[index+1],res[index+1],partialSum);
      } 
      mpfi_sqr(nominator,x);
      mpfi_sub_si(nominator, nominator,1);
  
      mpfi_set_si(pow,(2*index+1));
      mpfi_div_si(pow,pow,2);
      mpfi_pow(nominator,nominator,pow);
      mpfi_div(res[index+1],res[index+1],nominator);
    }
  //}
  for (index=0; index<=n-1; index++){
    mpfi_clear(coeffs_array[index]);
    mpfi_clear(powers_array[index]);
  }
  mpfi_clear(oldCoeff);
  mpfi_clear(partialSum);
  mpfi_clear(s1);
  mpfi_clear(pow);
  mpfi_clear(a1);
  mpfi_clear(nominator);
  return;
}

void baseFunction_diff(mpfi_t *res, int nodeType, mpfi_t x, int n) {
  mpfr_t oneHalf;
  int i;
  switch(nodeType) {
  
  case SQRT:
    mpfr_init2(oneHalf, getToolPrecision());
    mpfr_set_d(oneHalf, 0.5, GMP_RNDN);
    constantPower_diff(res, oneHalf, x, n);
    mpfr_clear(oneHalf);
    break;
  case EXP:
    exp_diff(res, x, n);
    break;
  case LOG:
    log_diff(res,x,n);
    break;
  case LOG_2:
    log2_diff(res,x,n);
    break;
  case LOG_10:
    log10_diff(res,x,n);
    break;
  case SIN:
    sin_diff(res,x,n);
    break;
  case COS:
    cos_diff(res,x,n);
    break;
  case TAN:
    tan_diff(res,x,n);
    break;
  case ASIN:
    asin_diff(res,x,n);
    break;
  case ACOS:
    acos_diff(res,x,n);
    break;
  case ATAN:
     atan_diff(res,x,n);
    break;
  case SINH:
    sinh_diff(res,x,n);
    break;
  case COSH:
    cosh_diff(res,x,n);
    break;
  case TANH:
    tanh_diff(res,x,n);
    break;
  case ASINH:
    asinh_diff(res,x,n);
    break;
  case ACOSH:
    acosh_diff(res,x,n);
    break;
  case ATANH:
    atanh_diff(res,x,n);
    break;
  case ABS:
    break;
  case DOUBLE:
    break;
  case DOUBLEDOUBLE:
    break;
  case TRIPLEDOUBLE:
    break;
  case ERF:
    break; 
  case ERFC:
    break;
  case LOG_1P:
    break;
  case EXP_M1:
    expm1_diff(res,x,n);
    break;
  case DOUBLEEXTENDED:
    break;
  case CEIL:
    break;
  case FLOOR:
    break;
  default:
    fprintf(stderr,"Error: AD: unknown unary function (%d) in the tree\n",nodeType);
    exit(1);
  }

  return;
}



/***********************************************************/
/***********************************************************/





/*this is the taylor model structure:
n- order: polynomial of degree n-1, remainder of order o(x^n)
rem_bound - bound for the remainder
poly_array - array of coeffs for the polynomial (mpfi's)
poly_bound - bound for the polynomial (helpful for computations)
x- interval on which the tm is computed
x0 - interval around the expansion point
*/
typedef struct tmdl {
int n; 
mpfi_t rem_bound;
mpfi_t *poly_array;
mpfi_t poly_bound;
mpfi_t x;
mpfi_t x0;

} tModel;
void ctMultiplication_TM(tModel*d,tModel*s, mpfi_t c);
void ctDivision_TM(tModel*d,tModel*s, mpfi_t c);
void polynomialBoundSharp(mpfi_t *bound,int n,mpfi_t *coeffs,mpfi_t x0Int,mpfi_t x);
  
/*This function creates an empty taylor model
*/
tModel* createEmptytModel(int n,  mpfi_t x0, mpfi_t x){
  tModel* t;
  int i;
 // printf("\nin createEmptyTM\n ");
  t= (tModel *)safeMalloc(sizeof(tModel));
  mpfi_init2(t->rem_bound, getToolPrecision());
  mpfi_init2(t->poly_bound,getToolPrecision());
  mpfi_init2(t->x,getToolPrecision());
  mpfi_set(t->x,x);
  mpfi_init2(t->x0, getToolPrecision());
  mpfi_set(t->x0,x0);
  t->n=n;
  t->poly_array= (mpfi_t *)safeMalloc(n*sizeof(mpfi_t));
  for(i=0;i<n;i++){
    mpfi_init2(t->poly_array[i], getToolPrecision());
    //mpfi_set_ui(t->poly_array[i],0);
  }
  return t;
}

/*This function gives a constant taylor model,
for exiting tModel t
*/

void consttModel(tModel*t, mpfi_t ct){ 
  int i,n;
  n=t->n;
  
  for(i=1;i<n;i++){
     mpfi_set_ui(t->poly_array[i],0);
  }
  
  mpfi_set(t->poly_array[0],ct);
  mpfi_set(t->poly_bound,ct);
  mpfi_set_ui(t->rem_bound,0); 

}

/*This function dealocates a taylor model
*/
void cleartModel(tModel *t){
  int i;
  for(i=0;i<t->n;i++) mpfi_clear(t->poly_array[i]);
  free(t->poly_array);
  mpfi_clear(t->rem_bound);
  mpfi_clear(t->poly_bound);  
  mpfi_clear(t->x);
  mpfi_clear(t->x0);
  free(t);
}

/*This function pretty prints a taylor model
*/
void printtModel(tModel *t){
  int i;
  printf("\nTaylor model of order, %d expanded in ", t->n);
  printInterval(t->x0);
  printf("\nCoeffs:");
  for(i=0;i<t->n;i++) {
    printInterval(t->poly_array[i]);
    printf(",");
  }  
  printf("r=");
  printInterval(t->rem_bound);
  printf(",b=");
  printInterval(t->poly_bound);  
  printf("\n");  
  }


/*This function sets a taylor model t 
with the values given by anoter tm tt
they implicitely have the same order expansion point
and interval

*/
void copytModel(tModel *t, tModel *tt){
  int i;
  for(i=0;i<tt->n;i++) {
    mpfi_set(t->poly_array[i],tt->poly_array[i]);
  }  
  mpfi_set(t->rem_bound,tt->rem_bound);
  mpfi_set(t->poly_bound,tt->poly_bound);  
  }

/*This function gets an mpfi from a node c:
if c constant after simplification-> ok
if not 0
*/
int mpfi_set_node( mpfi_t *r, node * c) {
  mpfi_t rr;
  node *cc;
  mpfi_init2(rr,getToolPrecision());
  if (c!=NULL){
    cc=simplifyTreeErrorfree(c);
    switch (cc->nodeType){
      case PI_CONST: mpfi_const_pi(rr);
      break;
      case CONSTANT:mpfi_set_fr(rr,*(cc->value));
      break;
      default: mpfi_set_ui(rr,0);
      break;
      }
    free(cc);
  }
  else mpfi_set_ui(rr,0);
  mpfi_set(*r,rr);
  mpfi_clear(rr);
  return 0;
}

void evaluateMpfiFunction(mpfi_t y, node *f, mpfi_t x, mp_prec_t prec){
rangetype yr, xr;

 //xr = (rangetype)safeMalloc(sizeof(rangetype));
 
 xr.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
 xr.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
 mpfr_init2(*(xr.a),getToolPrecision());
 mpfr_init2(*(xr.b),getToolPrecision());
 mpfi_get_left(*(xr.a), x);
 mpfi_get_right(*(xr.b), x); 
 
 
 yr.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
 yr.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
 mpfr_init2(*(yr.a),getToolPrecision());
 mpfr_init2(*(yr.b),getToolPrecision());
  
 evaluateRangeFunction(yr, f, xr, prec);
 
 mpfi_interv_fr(y, *(yr.a), *(yr.b));
 
 mpfr_clear(*(xr.a));
 mpfr_clear(*(xr.b));
 mpfr_clear(*(yr.a));
 mpfr_clear(*(yr.b));
 
  free(xr.a);
	free(xr.b);
	free(yr.a);
	free(yr.b);
 
}
/*This function computes an interval bound
for a polynomial given by coeffs and order n, 
on int x, using basic IA and Horner form
*/
void polynomialBound(mpfi_t *bound,int n,mpfi_t *coeffs,mpfi_t x){
  int i;
  mpfi_t r;
  mpfi_init2(r,getToolPrecision());
  mpfi_set(r,coeffs[n]);
  for (i=n-1;i>=0;i--){
  mpfi_mul(r,r,x);
  mpfi_add(r,r,coeffs[i]);
  }
  mpfi_set(*bound,r);
  mpfi_clear(r);
}
/*This function transforms a polynomial with interval coeffs
into a poly with mpfr coeffs and a small remainder
*/
void mpfr_get_poly(mpfr_t *rc, mpfi_t rest, int n, mpfi_t *gc, mpfi_t x){
  int i;
  mpfi_t *res;
  mpfi_t r;
  res= (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
  mpfi_init2(r,getToolPrecision());
  for (i=0; i<=n; i++){
    mpfi_mid(rc[i],gc[i]);
    mpfi_init2(res[i], getToolPrecision());
    mpfi_sub_fr(res[i],gc[i], rc[i]);
  }
  polynomialBound(&r,n,res,x);
  
  mpfi_set(rest,r);
  
  for (i=0; i<=n; i++){
    mpfi_clear(res[i]);  
  }
  free(res);
  mpfi_clear(r);
  return;

}

/*This function computes the tm for multiplication of two 
given tm's 
*/
void  multiplication_TM(tModel *t,tModel *c1, tModel *c2){
  //we will multiply two taylor models of order n; and obtain a new taylor model of order n;
  int n,i,j;
  mpfi_t *r;
  mpfi_t pow;
  tModel *tt;
  //printf("in multiplication TM");
  n=t->n;
  //aux tm for doing the multiplications
  
  tt=createEmptytModel(n,t->x0, t->x);
  
  r= (mpfi_t *)safeMalloc((n)*sizeof(mpfi_t));
  
  for(i=0;i<=n-1;i++){
    mpfi_init2(r[i], getToolPrecision());
    mpfi_set_ui(r[i],0);
    mpfi_set_ui(tt->poly_array[i],0);
    
  }
  mpfi_t temp1, temp2;
  mpfi_init2(temp1, getToolPrecision());
  mpfi_init2(temp2, getToolPrecision());
  
 /*We are multiplying taylor models, considering the relative error
 We are given:  (T1,delta1*(x-x0)^n
                (T2,delta2*(x-x0)^n
               
 The product is:(T1*T2|[0...n-1], {delta1*B(T2)+delta2*B(T1)+delta1*delta2*B((x-x0)^n)+B(T1*T2|[n...2*n-2])}  *(x-x0)^n */
 
  /*compute in temp1 delta1*delta2*B((x-x0)^n)*/
  mpfi_mul(temp1, c1->rem_bound, c2->rem_bound);
  mpfi_sub(temp2,t->x,t->x0);
  
  mpfi_init2(pow,getToolPrecision());
  mpfi_set_ui(pow,n);
  //printf("ajunge aici");
  //printInterval(temp2);
  //printInterval(pow);
  mpfi_pow(temp2,temp2,pow);
  //printInterval(temp2);
  mpfi_mul(temp1,temp2,temp1);
  
  
  
  //printf("compute in temp2 delta2*B(T1)");
  mpfi_mul(temp2, c1->poly_bound, c2->rem_bound);
  
  mpfi_add(tt->rem_bound,temp1,temp2);
  
  /*compute in temp1 delta1*B(T2)*/
  mpfi_mul(temp1, c1->rem_bound, c2->poly_bound);
  mpfi_add(tt->rem_bound, tt->rem_bound, temp1);
  
  /*compute the product of the two polynomials*/
  for(i=0; i<n;i++)
    for (j=0;j<n;j++){
      mpfi_mul(temp1,c1->poly_array[i], c2->poly_array[j]);
      if ((i+j)<n )
        mpfi_add(tt->poly_array[i+j],tt->poly_array[i+j],temp1);
      else
        mpfi_add(r[i+j-n],r[i+j-n],temp1);
    }
   /*compute bound for polynomial T1*T2|[n...2*n-2] scaled by (x-x0)^n*/
   polynomialBoundSharp(&temp1, n-2,r,t->x0,t->x);
   
   //printInterval(temp1);
   
   /*we add temp1 in the bound for the remainder*/
   mpfi_add(tt->rem_bound,tt->rem_bound,temp1);
   
  /*we compute the new polynomial bound for the new model*/
  polynomialBoundSharp(&temp1, n-1,tt->poly_array,t->x0,t->x);   
  mpfi_set(tt->poly_bound,temp1);
  
  mpfi_clear(temp1);
  mpfi_clear(temp2);
  mpfi_clear(pow);
  for(i=0;i<n-1;i++)
    mpfi_clear(r[i]);
  free(r); 
  
  /*set the result*/
  copytModel(t,tt);
  /*clear the aux tm*/
  cleartModel(tt);
 }

//-----------------------------------------------------------
//-----------------------------------------------------------

/*This function computes the tm for addition of two 
given tm's 
*/
void addition_TM(tModel *t,tModel *child1_tm, tModel *child2_tm){
  int i;
  int n;
  tModel *tt;
  n=t->n;
  //printf("in addition tm");
  tt=createEmptytModel(n,t->x0,t->x);
  for(i=0;i<n;i++)  mpfi_add(tt->poly_array[i], child1_tm->poly_array[i],child2_tm->poly_array[i]);
  mpfi_add(tt->rem_bound,child1_tm->rem_bound,child2_tm->rem_bound);
 // mpfi_add(tt->poly_bound,child1_tm->poly_bound,child2_tm->poly_bound);
  polynomialBoundSharp(&tt->poly_bound, n-1,tt->poly_array,t->x0,t->x);   
  copytModel(t,tt);
  cleartModel(tt);
}

/*base function*/
void base_TM(tModel *t,int nodeType, int n, mpfi_t x0, mpfi_t x){
  int i;
  tModel *tt;
  mpfi_t *nDeriv;
  mpfi_t fact;
  
 
  mpfi_init2(fact,getToolPrecision());
  mpfi_set_ui(fact,1);
  nDeriv= (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
  
  for(i=0;i<=n;i++){
    mpfi_init2(nDeriv[i], getToolPrecision());
  }
  
  tt=createEmptytModel(n,x0,x);
  
  baseFunction_diff( tt->poly_array,nodeType,x0, n-1);
  
  polynomialBoundSharp(&tt->poly_bound, n-1,tt->poly_array,t->x0,t->x);   
  
  baseFunction_diff(nDeriv,nodeType,x,n);
  //printInterval(x);
  for(i=1;i<n;i++){
  mpfi_mul_ui(fact,fact,i);
  mpfi_div(tt->poly_array[i],tt->poly_array[i],fact);
  }
  mpfi_mul_ui(fact,fact,n);
  mpfi_set(tt->rem_bound,nDeriv[n]);
  mpfi_div(tt->rem_bound,tt->rem_bound,fact);
  
  
  printf("basefunction interval:");
  printInterval(tt->rem_bound);
  copytModel(t,tt);
  cleartModel(tt);
}
/*polynomial_TM(tModel *t,node *f, int n, mpfr_t x0, mpfi_t x0Interv, mpfi_t x){
  int d;
  
  getDegree(&d,f);
  
  autodiff(& res,f,n-1,x0);
  for (i=0;i<n;i++)
  //setcoeffs tt
  
  polynomialSharpBound();
  AD(&res
  
  
  
  }
  
  
*/

/*composition*/
void composition_TM(tModel *t,tModel *g, tModel *f){
  int i;
  int n;
  tModel *tt, *partial_tmul,*tinterm,*tmul ;
  n=f->n;
  //mpfi_t zero;
  //mpfi_init2(zero,getToolPrecision());
  //mpfi_set_ui(zero,0);
  
  /*we have f taylor model in x0, over x*/
  /*we have g taylor model in f(x0) over bound of f(x)*/
  /* the resulting  tm is in x0 over x*/
  
  
  /*compute the */ 
  tt=createEmptytModel(n,f->x0,f->x);
  consttModel(tt,g->poly_array[0]);//initially, the tm contains g(f(x0));
  
  
  /*create the taylor polynomial for f(x)-f(x0)*/
  tinterm=createEmptytModel(n,f->x0,f->x);
  copytModel(tinterm,f);
  /*set the ct part of tinterm as 0*/
  mpfi_set_ui(tinterm->poly_array[0],0);
    
  
  tmul=createEmptytModel(n,f->x0,f->x);
  copytModel(tmul,tinterm);  
  
  
  partial_tmul=createEmptytModel(n,f->x0,f->x);
  copytModel(partial_tmul,tinterm);  
  
  
  for (i=1;i<n;i++){
  ctMultiplication_TM(tinterm,partial_tmul,g->poly_array[i]);
  addition_TM(tt,tt,tinterm);
  multiplication_TM(partial_tmul,partial_tmul,tmul);
 //  printtModel(partial_tmul);
  }
  
  //in tt we have now the good polynomial;
  //we have to compute the bounds.
  //(y-y0)^n
  
 
 

  //printtModel(tmul);  
  //printInterval(partial_tmul->rem_bound);
  mpfi_mul(partial_tmul->rem_bound,partial_tmul->rem_bound,g->rem_bound);
  mpfi_add(tt->rem_bound,tt->rem_bound,partial_tmul->rem_bound);
  polynomialBoundSharp(&tt->poly_bound, n-1,tt->poly_array,tt->x0,tt->x);   
  copytModel(t,tt);
  cleartModel(tt);
}




/*This function computes the tm for division
with a ct term of a given tm
*/
void ctDivision_TM(tModel*d,tModel*s, mpfi_t c){
  int i;
  int n;
  tModel *tt;
  n=s->n;
  tt=createEmptytModel(n,s->x0,s->x);
  copytModel(tt,s);
  for(i=0;i<n;i++)  mpfi_div(tt->poly_array[i],tt->poly_array[i], c);
  mpfi_div(tt->rem_bound,tt->rem_bound,c);
  mpfi_div(tt->poly_bound,tt->poly_bound,c);
  copytModel(d,tt);
  cleartModel(tt);
}


/*This function computes the tm for multiplication
with a ct term of a given tm
*/
void ctMultiplication_TM(tModel*d,tModel*s, mpfi_t c){
  int i;
  int n;
  n=s->n;
  tModel *tt;
  tt=createEmptytModel(n,s->x0,s->x);
  copytModel(tt,s);
  for(i=0;i<n;i++)  mpfi_mul(tt->poly_array[i], tt->poly_array[i],c);
  mpfi_mul(tt->rem_bound,tt->rem_bound,c);
  mpfi_mul(tt->poly_bound,tt->poly_bound,c);
  copytModel(d,tt);
  cleartModel(tt);
}



/*this function computes the taylor model for 1/x*/
void  varInv_TM(tModel *t,mpfi_t x0, mpfi_t x, int n){
    tModel *tt;
    mpfi_t *nDeriv;
    int i;
    printf("We are given the following parameters:");
    printInterval(x0);
    printInterval(x);
    printf("%d",n);    
    tt=createEmptytModel(n,x0,x); 
    mpfr_t minusOne;
    mpfi_t fact;    
    
    mpfr_init2(minusOne, getToolPrecision());
    mpfr_set_si(minusOne, -1, GMP_RNDN);

    constantPower_diff(tt->poly_array,minusOne, x0, n-1);
    
    nDeriv= (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    for(i=0;i<=n;i++){
      mpfi_init2(nDeriv[i], getToolPrecision());
    }
    constantPower_diff(nDeriv,minusOne, x, n);
    mpfi_init2(fact, getToolPrecision());
    mpfi_set_ui(fact,1);
    for(i=1;i<n;i++){
      mpfi_mul_ui(fact,fact,i);
      mpfi_div(tt->poly_array[i],tt->poly_array[i],fact);
    }
    mpfi_mul_ui(fact,fact,n);
    mpfi_set(tt->rem_bound,nDeriv[n]);
    mpfi_div(tt->rem_bound,tt->rem_bound,fact);
    printf("inversion function interval:");
    printInterval(tt->rem_bound);
    polynomialBoundSharp(&tt->poly_bound, n-1,tt->poly_array,tt->x0,tt->x);   
   
    copytModel(t,tt);
    cleartModel(tt);
    mpfr_clear(minusOne);    
    mpfi_clear(fact);

}


/*this function computes the taylor model for x^p*/
void  ctPowerVar_TM(tModel *t,mpfi_t x0, mpfi_t x, int n, mpfr_t p){
    tModel *tt;
    mpfi_t *nDeriv;
    int i;
    mpfi_t fact;
    printf("We are given the following parameters:");
    printInterval(x0);
    printInterval(x);
    printf("%d",n);    
    printMpfr(p);
    tt=createEmptytModel(n,x0,x); 
    
        
    
    constantPower_diff(tt->poly_array,p, x0, n-1);
    
    nDeriv= (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    for(i=0;i<=n;i++){
      mpfi_init2(nDeriv[i], getToolPrecision());
    }
    constantPower_diff(nDeriv,p, x, n);
    
    mpfi_init2(fact, getToolPrecision());
    mpfi_set_ui(fact,1);
    for(i=1;i<n;i++){
      mpfi_mul_ui(fact,fact,i);
      mpfi_div(tt->poly_array[i],tt->poly_array[i],fact);
    }
    mpfi_mul_ui(fact,fact,n);
    mpfi_set(tt->rem_bound,nDeriv[n]);
    mpfi_div(tt->rem_bound,tt->rem_bound,fact);
    printf("inversion function interval:");
    printInterval(tt->rem_bound);
    polynomialBoundSharp(&tt->poly_bound, n-1,tt->poly_array,tt->x0,tt->x);   
   
    copytModel(t,tt);
    cleartModel(tt);
    mpfi_clear(fact);

}

/*this function computes the taylor model for p^x*/
void  varCtPower_TM(tModel *t,mpfi_t x0, mpfi_t x, int n, mpfr_t p){
    tModel *tt;
    mpfi_t *nDeriv;
    int i;
    printf("We are given the following parameters:");
    printInterval(x0);
    printInterval(x);
    printf("%d",n);    
    tt=createEmptytModel(n,x0,x); 
    
    mpfi_t fact;    
    
    powerFunction_diff(tt->poly_array,p, x0, n-1);
    
    nDeriv= (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    for(i=0;i<=n;i++){
      mpfi_init2(nDeriv[i], getToolPrecision());
    }
    powerFunction_diff(nDeriv,p, x, n);
    mpfi_init2(fact, getToolPrecision());
    mpfi_set_ui(fact,1);
    for(i=1;i<n;i++){
      mpfi_mul_ui(fact,fact,i);
      mpfi_div(tt->poly_array[i],tt->poly_array[i],fact);
    }
    mpfi_mul_ui(fact,fact,n);
    mpfi_set(tt->rem_bound,nDeriv[n]);
    mpfi_div(tt->rem_bound,tt->rem_bound,fact);
    printf("p^x interval bound:");
    printInterval(tt->rem_bound);
    polynomialBoundSharp(&tt->poly_bound, n-1,tt->poly_array,tt->x0,tt->x);   
   
    copytModel(t,tt);
    cleartModel(tt);
    mpfi_clear(fact);

}





/*This function computes an interval bound
for a polynomial given by coeffs and order n, 
on int x, using basic IA and Horner form
*/
void polynomialBoundSharp(mpfi_t *bound,int n,mpfi_t *coeffs,mpfi_t x0,mpfi_t x){
  int i;
  mpfi_t r,xs;
  mpfi_init2(r,getToolPrecision());
  mpfi_init2(xs,getToolPrecision());
  mpfi_sub(xs,x,x0);
  mpfi_set(r,coeffs[n]);
  for (i=n-1;i>=0;i--){
  mpfi_mul(r,r,xs);
  mpfi_add(r,r,coeffs[i]);
  }
  mpfi_set(*bound,r);
  mpfi_clear(r);
  mpfi_clear(xs);
}

void reduceOrder_TM(tModel*d,tModel*s, int n){
  int i;
  int oldn;
  mpfi_t *remTerms;
  oldn=s->n;
  tModel *tt;
  mpfi_t pow;
  mpfi_t temp, temp2;
  tt=createEmptytModel(n,s->x0,s->x);
  for(i=0;i<n;i++){
    mpfi_set(tt->poly_array[i], s->poly_array[i]);
  }  
  //we are left with terms from n up to oldn-1
  remTerms= (mpfi_t *)safeMalloc((oldn-n)*sizeof(mpfi_t));
    for(i=n;i<oldn;i++){
      mpfi_init2(remTerms[i-n], getToolPrecision());
      mpfi_set(remTerms[i-n], s->poly_array[i]);
    }
  
  mpfi_init2(temp,getToolPrecision());  
  polynomialBoundSharp(&temp, oldn-n-1,remTerms,tt->x0,tt->x);     
      
  polynomialBoundSharp(&tt->poly_bound, n-1,tt->poly_array,tt->x0,tt->x);     
  
  mpfi_init2(pow,getToolPrecision());  
  mpfi_init2(temp2,getToolPrecision());  
  mpfi_set_ui(pow,oldn-n);
  mpfi_sub(temp2,tt->x,tt->x0);
  mpfi_pow(temp2,temp2,pow);
  
  mpfi_mul(temp2,temp2,s->rem_bound);
  mpfi_add(tt->rem_bound,temp2,temp);
  
  copytModel(d,tt);
  cleartModel(tt);
  mpfi_clear(temp);
  mpfi_clear(temp2);
  mpfi_clear(pow);
  for(i=0;i<oldn-n;i++){
     mpfi_clear(remTerms[i]);
  }
  free(remTerms);
  
}

void removeCoeffs_TM(tModel*d,tModel*s, int l){
  int i;
  int oldn,newn;
  
  
  tModel *tt;
  mpfi_t pow;
  mpfi_t temp;
  
  //we know that a0,...,al are 0;
  //we create a tm of order oldOrder - (l+1)
  oldn=s->n;
  newn=oldn-l-1;
  tt=createEmptytModel(newn,s->x0,s->x);
  
  for (i=l+1;i<oldn;i++){
    mpfi_set(tt->poly_array[i-l-1],s->poly_array[i]);
  }
  
  //remainder: it is the same, since this remove coeffs is equivalent to a formal simplification
  //by (x-x0)^(oldn-newn)
  //Delta * (x-x0)^oldn --> Delta*(x-x0)^newn 
  
  mpfi_set(tt->rem_bound,s->rem_bound);
  /*mpfi_init2(pow,getToolPrecision());  
  mpfi_init2(temp,getToolPrecision());  
  mpfi_set_ui(pow,oldn-newn);
  mpfi_sub(temp,tt->x,tt->x0);
  mpfi_pow(tt->rem_bound,temp,pow);
  */    
  polynomialBoundSharp(&tt->poly_bound, newn-1,tt->poly_array,tt->x0,tt->x);     
  
  copytModel(d,tt);
  cleartModel(tt);
  /*mpfi_clear(temp);
  mpfi_clear(pow);
  */
}



void taylor_model(tModel *t, node *f, int n, mpfi_t x0, mpfi_t x) {
  int i;
  
  node *simplifiedChild1, *simplifiedChild2;
  mpfi_t temp1,temp2;
  node **coefficients;
  mpfi_t *rpoly, *boundRpoly;
  tModel *tt,*tPoly, *child1_tm, *child2_tm, *ctPowVar_tm, *varCtPower_tm, *logx_tm, *expx_tm, *logf_tm;
  int d;
  mpfi_t powx,powy;
  mpfi_t ct;
  mpfr_t zero;
  
  
  
  /*
  //printf("We entered taylor model order: %d\n",n);
  if (isPolynomial(f) ){
  tt=createEmptytModel(n,x0,x);
  printf("We found a polynomial!!!\n");
  
  polynomial_tm(tt,f,n,x0,x);
  
  
  
  */
  
  
  
  /*we create a taylor model for this polynomial
  it is in fact in the canonical base;
  in order to get the polynomial in base (x-x0),
  we compose it with (x-x0, x0, x,x-[x0],0)*/
  /*
  mpfi_t translatedI;
  mpfi_init2(translatedI,getToolPrecision());
  mpfi_sub_fr(translatedI,x,x0);
  mpfr_init2(zero,getToolPrecision());
  mpfr_set_ui(zero,0,GMP_RNDN);
    
  tPoly=createEmptytModel(n,zero,translatedI);
  getCoefficients(&d, &coefficients, f);
  //for(i=0;i<=d;i++) if (coefficients[i]!=NULL) printTree(coefficients[i]);
    if (d<n){
      printf("The degree of polynomial smaller : %d than the taylor model degree\n",d);
      for(i=0;i<=d;i++) {
        printf("%d\n",i);
        
        mpfi_set_node(&(tPoly->poly_array[i]),coefficients[i]);
        printInterval(tPoly->poly_array[i]);
      }
      for(i=d+1;i<n;i++) {
        mpfi_set_ui(tPoly->poly_array[i],0);
      }
       printf("we've set the coefficients\n");
      //compute the bound for the polynomial, put the bound for tm= 0;
      polynomialBoundSharp(&tPoly->poly_bound,d,tPoly->poly_array,tPoly->x0, tPoly->x);
      mpfi_set_ui(tPoly->rem_bound,0);
      printf("passed the bond");
      printtModel(tPoly);
    }
    else {
      //printf("The degree of polynomial bigger : %d then the taylor model degree\n",d);
      for(i=0;i<n;i++) mpfi_set_node(&(tPoly->poly_array[i]),coefficients[i]);
      //we keep the bound on the n degree polynomial
        polynomialBoundSharp(&tPoly->poly_bound,n-1,tPoly->poly_array,tPoly->x0,tPoly->x);
        
      //for n pana la d, creaza polynomul, evalueaza-l pe interval si pe urma returneaza in rem_bound;
      rpoly=(mpfi_t *)safeMalloc((d+1-n)*sizeof(mpfi_t));
     
      for (i=n;i<=d;i++){
        //the polynomial has the form x^(n)*(coeff[n] + coeff[n+1]*x+...+ coeff[d]*x^(d-n) )
        //we create the polynomial of degree d-n; and evaluate it on the interval
     
        mpfi_init2(rpoly[i-n], getToolPrecision());
        mpfi_set_node( &rpoly[i-n], coefficients[i]);
      }
        polynomialBoundSharp(&tPoly->rem_bound,d-n,rpoly,tPoly->x0,tPoly->x);
        //printInterval(*boundRpoly);
        
    }
    child1_tm=createEmptytModel(n,x0,x);
    
    mpfi_set_fr(child1_tm->poly_array[0],x0);
    mpfi_set_ui(child1_tm->poly_array[1],1);
    mpfi_set_ui(child1_tm->rem_bound,0);
    mpfi_set(child1_tm->poly_bound,x);
    for (i=2;i<n;i++){
        mpfi_set_ui(child1_tm->poly_array[i],0);    
    }
    composition_TM(tt,tPoly,child1_tm);
    copytModel(t,tt);
    cleartModel(tPoly);
    cleartModel(tt);
    cleartModel(child1_tm);
  }
  else {*/
  
  switch (f->nodeType) {
  case VARIABLE:
  
  
  tt=createEmptytModel(n,x0,x); 
  mpfi_set(tt->poly_array[0],x0);
  if (n==1){
    mpfi_set_ui(tt->rem_bound,1);
    mpfi_set(tt->poly_bound,x0);
  }
  else{
    mpfi_set_ui(tt->rem_bound,0);
    mpfi_set_ui(tt->poly_array[1],1);
    for (i=2;i<n;i++){
      mpfi_set_ui(tt->poly_array[i],0);
    }
    mpfi_set(tt->poly_bound,x);
  }
  copytModel(t,tt);
  cleartModel(tt);
  break;
  case PI_CONST:
  case CONSTANT:
  //printf("in constant");
  mpfi_init2(ct, getToolPrecision());
  tt=createEmptytModel(n,x0,x); 
  mpfi_set_node( &ct, f);
  consttModel(tt,ct);
  mpfi_set_ui(tt->rem_bound,0);
  mpfi_set(tt->poly_bound,ct);
  copytModel(t,tt);
  cleartModel(tt);
  mpfi_clear(ct);
  //printf("done");
  break;
    
  
  
  case NEG:
    tt=createEmptytModel(n,x0,x);
    //create a new empty taylor model the child
    child1_tm=createEmptytModel(n, x0, x);
    //call taylor_model on the child
    taylor_model(child1_tm, f->child1,n,x0,x);
    //do the necessary chages from child to parent
    for(i=0;i<n;i++)  mpfi_neg(tt->poly_array[i], child1_tm->poly_array[i]);
    mpfi_neg(tt->rem_bound,child1_tm->rem_bound);
    mpfi_neg(tt->poly_bound,child1_tm->poly_bound);
    copytModel(t,tt);
    //clear old taylor models
    cleartModel(child1_tm);
    cleartModel(tt);
    break;

  case ADD:
  
  //create a new empty taylor model the children
    tt=createEmptytModel(n,x0,x); 
    child1_tm=createEmptytModel(n, x0,x);
    child2_tm=createEmptytModel(n, x0,x);
    //call taylor_model on the children
    taylor_model(child1_tm, f->child1,n,x0,x);
    taylor_model(child2_tm, f->child2,n,x0,x);
    //do the necessary chages from children to parent
    /*for(i=0;i<=n;i++)  mpfi_add(t->poly_array[i], child1_tm->poly_array[i],child2_tm->poly_array[i]);
    mpfi_add(t->rem_bound,child1_tm->rem_bound,child2_tm->rem_bound);
    mpfi_add(t->poly_bound,child1_tm->poly_bound,child2_tm->poly_bound);
    */
    addition_TM(tt,child1_tm, child2_tm);
    //clear old taylor model
    cleartModel(child1_tm);
     cleartModel(child2_tm);
     copytModel(t,tt);
     cleartModel(tt);
    break;

  case SUB:
    printf("in SUB");
    tt=createEmptytModel(n,x0,x); 
    //create a new empty taylor model the children
    child1_tm=createEmptytModel(n, x0,x);
    child2_tm=createEmptytModel(n, x0,x);
    //call taylor_model on the children
    taylor_model(child1_tm, f->child1,n,x0,x);
    taylor_model(child2_tm, f->child2,n,x0,x);
    //do the necessary chages from children to parent
    for(i=0;i<n;i++)  mpfi_sub(tt->poly_array[i], child1_tm->poly_array[i],child2_tm->poly_array[i]);
    mpfi_sub(tt->rem_bound,child1_tm->rem_bound,child2_tm->rem_bound);
    mpfi_sub(tt->poly_bound,child1_tm->poly_bound,child2_tm->poly_bound);
    
    //clear old taylor model
    cleartModel(child1_tm);
    cleartModel(child2_tm);
     copytModel(t,tt);
     cleartModel(tt);
    
    break;

  case MUL:
     tt=createEmptytModel(n,x0,x); 
   //create a new empty taylor model the children
    child1_tm=createEmptytModel(n,x0,x);
    child2_tm=createEmptytModel(n, x0,x);
    //call taylor_model on the children
    taylor_model(child1_tm, f->child1,n,x0,x);
    taylor_model(child2_tm, f->child2,n,x0,x);
    
    //do the necessary chages from children to parent
    multiplication_TM(tt,child1_tm, child2_tm);
    
    //clear old taylor model
    cleartModel(child1_tm);
    cleartModel(child2_tm);
     copytModel(t,tt);
     cleartModel(tt);
    
    break;

  case DIV:
  //child1 * inverse(child2)
  printf("division");
  tt=createEmptytModel(n,x0,x); 
    
  //check whether g(x0)is zero
  mpfi_t gx0, powx,rangeg;
  tModel *ttt, *inv_tm, *child1Extended_tm, *child2Extended_tm, *child1RemoveCoeffs_tm,*child2RemoveCoeffs_tm; 
  mpfi_init2(gx0, getToolPrecision());
  evaluateMpfiFunction(gx0,f->child2, x0, getToolPrecision())  ;
  ttt=createEmptytModel(n,x0,x);  
  if (mpfi_is_zero(gx0)){
  int k;  
    k=10;
    //create a new empty taylor model the child
    child1Extended_tm=createEmptytModel(n+k,x0,x);
    child2Extended_tm=createEmptytModel(n+k,x0,x);
    
    //call taylor_model for the children
    taylor_model(child1Extended_tm, f->child1,n+k,x0,x);
    taylor_model(child2Extended_tm, f->child2,n+k,x0,x);
    
    //reduce order taylor model
    int l;
    l=0;
    while((mpfi_is_zero(child1Extended_tm->poly_array[l])) && (mpfi_is_zero(child2Extended_tm->poly_array[l]))){
    l++;}
    
    printf("The order of the zero is: %d",l);
    
    
    //create a new empty taylor model the child
    child1RemoveCoeffs_tm=createEmptytModel(n+k-l,x0,x);
    child2RemoveCoeffs_tm=createEmptytModel(n+k-l,x0,x);
    
    
    removeCoeffs_TM(child1RemoveCoeffs_tm, child1Extended_tm, l-1);
    removeCoeffs_TM(child2RemoveCoeffs_tm, child2Extended_tm, l-1);
    
    
    
    child1_tm=createEmptytModel(n,x0,x);
    child2_tm=createEmptytModel(n,x0,x);
    
    reduceOrder_TM(child1_tm,child1RemoveCoeffs_tm,n);
    reduceOrder_TM(child2_tm,child2RemoveCoeffs_tm,n);
    printf("\nThe reduced tm are:\n");
    printtModel(child1_tm);
    printtModel(child2_tm);
    
    
    cleartModel(child1RemoveCoeffs_tm);
    cleartModel(child2RemoveCoeffs_tm);
    
    cleartModel(child1Extended_tm);
    cleartModel(child2Extended_tm);
    
    
  }//we have reduced the poles, we have taylor models of order n, we apply the basic division
  else{//just create tms or order n directly
  //create a new empty taylor model the child
    child1_tm=createEmptytModel(n,x0,x);
    child2_tm=createEmptytModel(n,x0,x);
    //call taylor_model for the children
    taylor_model(child1_tm, f->child1,n,x0,x);
    taylor_model(child2_tm, f->child2,n,x0,x);
  }
    mpfi_set(gx0,child2_tm->poly_array[0]);
    mpfi_init2(rangeg, getToolPrecision());
    
    mpfi_init2(powx, getToolPrecision());
    mpfi_set_ui(powx,n);
    
    mpfi_sub(rangeg, child2_tm->x,child2_tm->x0);
    mpfi_pow(rangeg,rangeg,powx);
    mpfi_mul(rangeg,rangeg,child2_tm->rem_bound);
    mpfi_add(rangeg,rangeg, child2_tm->poly_bound);
    
    inv_tm=createEmptytModel(n,gx0,rangeg);
    
    varInv_TM(inv_tm,gx0,rangeg,n);
    composition_TM(ttt,inv_tm,child2_tm);
    
    multiplication_TM(tt, ttt, child1_tm);
    //clear old child
    cleartModel(child1_tm);
    cleartModel(child2_tm);
    cleartModel(inv_tm);
    cleartModel(ttt);
    copytModel(t,tt);
    cleartModel(tt);
    mpfi_clear(rangeg);
    mpfi_clear(powx);
    
    mpfi_clear(gx0);
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
    tt=createEmptytModel(n,x0,x); 
   //create a new empty taylor model the child
    child1_tm=createEmptytModel(n,x0,x);
    child2_tm=createEmptytModel(n,x0,x);
    //call taylor_model on the child
    taylor_model(child1_tm, f->child1,n,x0,x);
    //compute tm for the basic case
    mpfi_t fx0,rangef,pow;
        
    mpfi_init2(fx0,getToolPrecision());
    mpfi_init2(rangef, getToolPrecision());
    
    mpfi_init2(pow, getToolPrecision());
    mpfi_set_ui(pow,n);
    
    evaluateMpfiFunction(fx0,f->child1,x0,getToolPrecision());
    
    mpfi_sub(rangef, x,x0);
    mpfi_pow(rangef,rangef,pow);
    mpfi_mul(rangef,rangef,child1_tm->rem_bound);
    mpfi_add(rangef,rangef, child1_tm->poly_bound);
    base_TM(child2_tm,f->nodeType,n,fx0, rangef);
    composition_TM(tt,child2_tm, child1_tm);
    //clear old child
    cleartModel(child1_tm);
    cleartModel(child2_tm);
    copytModel(t,tt);
    cleartModel(tt);
    mpfi_clear(fx0);    
    mpfi_clear(rangef);    
    mpfi_clear(pow);
        
  break;
  case POW:
  
  //create the result taylorModel
  tt=createEmptytModel(n,x0,x); 
  
      simplifiedChild2=simplifyTreeErrorfree(f->child2);
      simplifiedChild1=simplifyTreeErrorfree(f->child1);
      
      if ((simplifiedChild2->nodeType==CONSTANT) &&(simplifiedChild1->nodeType==CONSTANT)) { //we have the ct1^ct2 case
         // printf("We are in the  ct1^ct2 case");       
         mpfi_init2(temp1, getToolPrecision());
         mpfi_set_fr(temp1, *(simplifiedChild1->value));
         mpfi_init2(temp2, getToolPrecision());
         mpfi_set_fr(temp2, *(simplifiedChild2->value));
         mpfi_pow(temp1,temp1,temp2);
         consttModel(tt,temp1);
         mpfi_clear(temp1);
         mpfi_clear(temp2);
      }
      else if (simplifiedChild2->nodeType==CONSTANT) { //we have the f^p case
        printf("We are in the  f^p case");        
        mpfi_t powx,fx0,rangef;        
        //create a new empty taylor model the child
        child1_tm=createEmptytModel(n,x0,x);
        //call taylor_model for the child
        taylor_model(child1_tm, f->child1,n,x0,x);
        printf("\n\n-----------taylormodel child1: \n");
        printtModel(child1_tm);
        printf("-----------------------------\n");
        mpfi_init2(fx0,getToolPrecision());
        evaluateMpfiFunction(fx0,f->child1,x0,getToolPrecision());
        
        mpfi_init2(rangef, getToolPrecision());
        mpfi_init2(powx, getToolPrecision());
        mpfi_set_ui(powx,n);
        mpfi_sub(rangef, child1_tm->x,child1_tm->x0);
        mpfi_pow(rangef,rangef,powx);
        mpfi_mul(rangef,rangef,child1_tm->rem_bound);
        mpfi_add(rangef,rangef, child1_tm->poly_bound);
        //create tm for x^p over ragef in fx0
        ctPowVar_tm=createEmptytModel(n,fx0,rangef);
        
        ctPowerVar_TM(ctPowVar_tm,fx0,rangef,n,*(simplifiedChild2->value) );
        
        printf("\n\n-----------taylormodel child1: \n");
        printtModel(ctPowVar_tm);
        printf("-----------------------------\n");
        
        
        composition_TM(tt,ctPowVar_tm,child1_tm);
    
        //clear old child
        cleartModel(child1_tm);
        cleartModel(ctPowVar_tm);
        copytModel(t,tt);
        cleartModel(tt);
        mpfi_clear(rangef);
        mpfi_clear(powx);
        mpfi_clear(fx0);
      } 
       else if (simplifiedChild1->nodeType==CONSTANT) { //we have the p^f case
        
        mpfi_t powx,fx0,rangef;        
        //create a new empty taylor model the child
        child2_tm=createEmptytModel(n,x0,x);
        //call taylor_model for the child
        taylor_model(child2_tm, f->child2,n,x0,x);
        
        mpfi_init2(fx0,getToolPrecision());
        evaluateMpfiFunction(fx0,f->child1,x0,getToolPrecision());
        
        mpfi_init2(rangef, getToolPrecision());
        mpfi_init2(powx, getToolPrecision());
        mpfi_set_ui(powx,n);
        mpfi_sub(rangef, child2_tm->x,child2_tm->x0);
        mpfi_pow(rangef,rangef,powx);
        mpfi_mul(rangef,rangef,child2_tm->rem_bound);
        mpfi_add(rangef,rangef, child2_tm->poly_bound);
        //create tm for p^x over ragef in fx0
        varCtPower_tm=createEmptytModel(n,fx0,rangef);
        
        varCtPower_TM(varCtPower_tm,fx0,rangef,n,*(simplifiedChild1->value) );
        composition_TM(tt,varCtPower_tm,child2_tm);
    
        //clear old child
        cleartModel(child2_tm);
        cleartModel(varCtPower_tm);
        copytModel(t,tt);
        cleartModel(tt);
        mpfi_clear(rangef);
        mpfi_clear(powx);
        mpfi_clear(fx0);
        }
      else {
      //printf("We are in the  f^g case");  
      //exp(g log(f))
        mpfi_t powx,fx0,rangef;        
        //create a new empty taylor model the children
        child1_tm=createEmptytModel(n,x0,x);
        child2_tm=createEmptytModel(n,x0,x);
        //call taylor_model for child 2 = g
        taylor_model(child2_tm, f->child2,n,x0,x);
        
        //call taylor_model for child 1 = f
        taylor_model(child1_tm, f->child1,n,x0,x);
        
        //create  taylor_model for log (child 1) = log(f)
                
        mpfi_init2(fx0,getToolPrecision());
        evaluateMpfiFunction(fx0,f->child1,x0,getToolPrecision());
        
        mpfi_init2(rangef, getToolPrecision());
        mpfi_init2(powx, getToolPrecision());
        mpfi_set_ui(powx,n);
        mpfi_sub(rangef, child1_tm->x,child1_tm->x0);
        mpfi_pow(rangef,rangef,powx);
        mpfi_mul(rangef,rangef,child1_tm->rem_bound);
        mpfi_add(rangef,rangef, child1_tm->poly_bound);
        
        logx_tm=createEmptytModel(n,fx0,rangef);
        base_TM(logx_tm,LOG,n,fx0, rangef);
        logf_tm=createEmptytModel(n,x0,x);
        composition_TM(logf_tm,logx_tm,child1_tm);
        //-------------------------------------------
        
        
        //multiply g by log f
        ttt=createEmptytModel(n,x0,x);
        multiplication_TM(ttt,logf_tm,child2_tm);
        
        //------------------------------------------
        mpfi_init2(gx0,getToolPrecision());
        evaluateMpfiFunction(gx0,f->child2,x0,getToolPrecision());
        mpfi_log(fx0,fx0);
        mpfi_mul(gx0,gx0,fx0);
        
        mpfi_set_ui(powx,n);
        mpfi_sub(rangef, ttt->x,ttt->x0);
        mpfi_pow(rangef,rangef,powx);
        mpfi_mul(rangef,rangef,ttt->rem_bound);
        mpfi_add(rangef,rangef, ttt->poly_bound);
        
        expx_tm=createEmptytModel(n,gx0,rangef);
        base_TM(expx_tm,EXP,n,gx0, rangef);
        composition_TM(tt,expx_tm,ttt);
               
    
        //clear old child
        cleartModel(child2_tm);
        cleartModel(child1_tm);
        cleartModel(ttt);
        
        copytModel(t,tt);
        cleartModel(tt);
        mpfi_clear(rangef);
        mpfi_clear(powx);
        mpfi_clear(fx0);   
        mpfi_clear(gx0); 
    }
    free_memory(simplifiedChild2);
    free_memory(simplifiedChild1);
  break;
  
  
  
  
  
  case LIBRARYFUNCTION:
  break;

  default:
   fprintf(stderr,"Error: TM: unknown identifier (%d) in the tree\n",f->nodeType);
   exit(1);
//  }
}
  return;
}




void taylorform(node **T, chain **errors, mpfi_t **delta,
		node *f, int n,	mpfr_t x0, mpfi_t *d, int mode) {

  printf("taylorform: f  = ");
  printTree(f);
  printf("\n");
  printf("taylorform: n  = %d\n",n);
  printf("taylorform: x0 = ");
  printMpfr(x0);
  if (d != NULL) {
    printf("taylorform: d  = ");
    printInterval(d);
    printf("\n");
  } else {
    printf("taylorform: no domain d given\n");
  }
  if (mode == ABSOLUTE) {
    printf("taylorform: absolute mode\n");
  } else {
    if (mode == RELATIVE) {
      printf("taylorform: relative mode\n");
    } else {
      printf("taylorform: mode %d\n",mode);
    }
  }
tModel *t;
mpfi_t x0Int;
mpfi_init2(x0Int,getToolPrecision());
mpfi_set_fr(x0Int,x0);
t=createEmptytModel(n,x0Int,*d);
//printf("we have created an emptytm");  

taylor_model(t,f,n,x0Int,*d);

printtModel(t);
mpfi_clear(x0Int);
}
