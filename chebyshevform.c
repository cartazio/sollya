/*

Copyright 2009 by 

Laboratoire de l'Informatique du Parallélisme, 
UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668

Contributors M.Joldes

mioara.joldes@ens-lyon.fr


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

/*
  For compiling this file:
    gcc -fPIC -Wall -c chebyshevform.c
    gcc -fPIC -shared -o chebyshevform chebyshevform.o


  Within Sollya:
    > externalproc(CM, "./chebyshevform", (function, range,range, integer) ->list of function);

  And then, for instance:
    > CM(exp(x)*cos(x), [2.5; 2.6], 10);

*/

/*******************************************************************************/
/*********************Functions related to interpolation***********************/
/*******************************************************************************/
/*******************************************************************************/
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
node* interpolation(mpfr_t *coeffArray,  mpfi_t *bound, node *f, mpfi_t x, int n);

void findNewSetRoots(mpfr_t * newRoots, mpfr_t *roots, int nrRoots, mpfr_t *chebArray, int n){
  int i,j;
  mpfr_t dist,currentDist;
  mpfr_init2(dist, getToolPrecision());
  mpfr_init2(currentDist, getToolPrecision());
  /*for each cheby point, find the point that is closest to it*/
  
  /*printf("\nthe cheby points are\n");
     for(i=0;i<n;i++){
     printMpfr(chebArray[i]);
     }
     
     printf("\nthe roots found are\n");
     for(i=0;i<nrRoots;i++){
     printMpfr(roots[i]);
     }
   */  
    for(i=0;i<n;i++){
     mpfr_init2(newRoots[i], getToolPrecision());
     mpfr_sub(dist,roots[0],chebArray[i],GMP_RNDN);
     mpfr_abs(dist,dist,GMP_RNDN);
     mpfr_set(newRoots[i],roots[0],GMP_RNDN);
     /*consider the first root as being closer*/
     for(j=1;j<nrRoots;j++){
      mpfr_sub(currentDist,roots[j],chebArray[i],GMP_RNDN);
      mpfr_abs(currentDist,currentDist,GMP_RNDN);
     
     
      if ((mpfr_cmp(currentDist, dist))<=0) {
        mpfr_set(newRoots[i], roots[j],GMP_RNDN);
        mpfr_set(dist,currentDist,GMP_RNDN);
      }
     } 
     
    }
    /*printf("\nthe roots found are\n");
     for(i=0;i<n;i++){
     printMpfr(newRoots[i]);
     }
     */
     mpfr_clear(dist);
     mpfr_clear(currentDist);
     
   }
   
void evaluateMpfiFunction(mpfi_t y, node *f, mpfi_t x, mp_prec_t prec){
rangetype yr, xr;

 //xr = (rangetype)safeMalloc(sizeof(rangetype));
 
 xr.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
 xr.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
 mpfr_init2(*(xr.a),prec);
 mpfr_init2(*(xr.b),prec);
 mpfi_get_left(*(xr.a), x);
 mpfi_get_right(*(xr.b), x); 
 
 
 yr.a = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
 yr.b = (mpfr_t *) safeMalloc(sizeof(mpfr_t));
 mpfr_init2(*(yr.a),prec);
 mpfr_init2(*(yr.b),prec);
  
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
int getIntervalAroundRoot(mpfi_t *res, node *dif, mpfr_t r, mp_prec_t prec){

  mpfi_t roundedr, difr;
  
  int rounded=0;
  mp_prec_t precForRounding;
  
  precForRounding=prec;
  
  mpfi_init2(roundedr,precForRounding);
  mpfi_set_fr(roundedr, r);
  //printf("in getInterval around root");
  while ((rounded!=1)&&(precForRounding>20)){
    mpfi_init2(difr,precForRounding);
    evaluateMpfiFunction(difr, dif,roundedr,precForRounding);
    //printInterval(difr);
    if ((mpfi_has_zero(difr))>0){
     rounded=1; 
     mpfi_init2(*res, precForRounding);
     mpfi_set(*res,roundedr);
  
    }
    else{
      precForRounding=precForRounding-10;
      mpfi_clear(roundedr);
      mpfi_init2(roundedr,precForRounding);  
      mpfi_set_fr(roundedr,r);
   }   
 }
 mpfi_clear(difr);
 return rounded;
 }
        
void evaluateW(mpfi_t res, mpfi_t extrem, mpfi_t *rootIntervals, int n){
  mpfi_t prod, term;
  int i;
  mpfi_init2(prod, getToolPrecision());
  mpfi_init2(term, getToolPrecision());
  mpfi_set_ui(prod,1);
  for (i=0; i<n; i++){
    mpfi_sub(term, extrem, rootIntervals[i]);
    mpfi_mul(prod, term, prod);
  }
  mpfi_set(res, prod);
  mpfi_clear(prod);
  mpfi_clear(term);
}
      
void mpfi_get_mid(mpfr_t res, mpfi_t op1, mpfi_t op2){

  mpfr_t op1Fr, op2Fr;
  mpfi_t newInt;
  
  mpfi_init2(newInt, getToolPrecision());
  mpfr_init2(op1Fr, getToolPrecision());
  mpfr_init2(op2Fr, getToolPrecision());  
  
  mpfi_get_right(op1Fr, op1);
  mpfi_get_left(op2Fr, op2);
  mpfi_interv_fr(newInt, op1Fr, op2Fr);
  mpfi_mid(res,newInt);
  
  mpfi_clear(newInt);
  mpfr_clear(op1Fr);
  mpfr_clear(op2Fr);  

}



int findWExtremum(mpfi_t wExtrem,int i,mpfr_t chebExtrem, int freeDegrees, mpfi_t *rootsIntervals){
    mp_prec_t prec;
    
    
    mpfi_t leftRoot, rightRoot;
    mpfi_t extremPossible;
    mpfi_t extremInt;
    mpfi_t extremLeftInt, extremRightInt;
    mpfr_t extremPossibleFr, extremRight, extremLeft;
    
    mpfi_t wleft, wmiddle, wright;
    
    mpfr_t ecart;
    
    int hill,disjoint,notFound;
    
    mpfi_init2(leftRoot, getToolPrecision());
    mpfi_init2(rightRoot, getToolPrecision());
    mpfi_set(leftRoot, rootsIntervals[i]);
    mpfi_set(rightRoot, rootsIntervals[i+1]);
    
    /*mpfi_init2(leftBorder, getToolPrecision());
    mpfi_init2(rightBorder, getToolPrecision());
    mpfi_set(leftBorder, rootsIntervals[i]);
    mpfi_set(rightBorder, rootsIntervals[i+1]);
    */
        
    //between the two Roots lies the extremum
    //closer then ever to the chebExtrem
    
    mpfi_init2(extremPossible, getToolPrecision());
    mpfi_set_fr(extremPossible, chebExtrem); 
    
    //find the general direction: hill or valley
    mpfi_init2(wleft, getToolPrecision());
    mpfi_init2(wright, getToolPrecision());
    mpfi_init2(wmiddle, getToolPrecision());  
    
    
    evaluateW(wleft, leftRoot, rootsIntervals, freeDegrees);
    evaluateW(wright, rightRoot, rootsIntervals, freeDegrees);
    evaluateW(wmiddle, extremPossible, rootsIntervals, freeDegrees);
    
    printf("First evaluation:\n");
    
    printInterval(wleft);
    printInterval(wright);
    printInterval(wmiddle);
    
        
    if (((mpfi_cmp(wleft,wmiddle))>0) && ((mpfi_cmp(wmiddle,wright))<0)  ){
    hill=0;
    }
    else if (((mpfi_cmp(wleft,wmiddle))<0) && ((mpfi_cmp(wmiddle,wright))>0)  ){
    hill=1;
    }
    else{
    printf("Something really sucks in this algo....");
    hill=-1;
    return hill;
   }
   if (hill==1)
   printf("We are in a hill %d\n", hill);
   else
   printf("We are in a valley %d\n", hill);
   
   //find the direction around the current potential extremum
   notFound=1;
   mpfr_init2(extremPossibleFr, getToolPrecision());
   mpfr_set(extremPossibleFr, chebExtrem, GMP_RNDN);
    int iteration=0;
   while((notFound==1)&&(iteration<30)){
    iteration++;
    printf("We start searching...\n");
    disjoint=0;
    prec=getToolPrecision();
    
    mpfi_set_fr(extremPossible, extremPossibleFr);
    evaluateW(wmiddle, extremPossible, rootsIntervals, freeDegrees);
    //printf("We evaluated\n");
    //printInterval(wmiddle);
    mpfr_init2(ecart, getToolPrecision());
    
    while (disjoint!=1) {
    printf("The blowing factor is: %d", prec);
    
      mpfi_init2(extremInt, prec);
      mpfr_set_ui(ecart, 2, GMP_RNDN);  
      mpfr_pow_si(ecart, ecart, -(prec-10), GMP_RNDN);
      mpfi_set_fr(extremInt, extremPossibleFr);
      printf("Before blowing");
      printInterval(extremInt);
      mpfi_increase(extremInt, ecart);
      printf("After blowing");
      printInterval(extremInt);
     
      mpfr_init2(extremLeft, prec);
      mpfr_init2(extremRight,prec);
        
      mpfi_init2(extremLeftInt, prec);
      mpfi_init2(extremRightInt, prec);
   
      mpfi_get_left(extremLeft,extremInt);
      mpfi_get_right(extremRight, extremInt);
    
      mpfi_set_fr(extremLeftInt, extremLeft);
      mpfi_set_fr(extremRightInt, extremRight);
       
      evaluateW(wleft, extremLeftInt, rootsIntervals, freeDegrees);
      evaluateW(wright, extremLeftInt, rootsIntervals, freeDegrees);
    
      if (((mpfi_cmp(wleft,wmiddle))==0) || ((mpfi_cmp(wright,wmiddle))==0)  ){
        prec=prec-10;
      }
      else{
        disjoint=1;
      }
    }
    printf("First disjoint found!!");
    
    //compute the direction
    //if we are in a valley, three are 3 posibilities
    if (hill==0){
      //order: wleft>wmiddle, wmiddle<wright
      if (((mpfi_cmp(wleft,wmiddle))>0) && ((mpfi_cmp(wmiddle,wright))<0)  ){
      //we have found a min
      printf("We have found a min:");
      printInterval(extremInt);
      notFound=0;
      }
      //order: wleft<wmiddle, wmiddle<wright
      else if   ( ((mpfi_cmp(wleft,wmiddle))<0) && ((mpfi_cmp(wmiddle,wright))<0)  ){
        //we are in a valley and should search more to the left
        //newLeftRoot=oldLeftRoot;
        //newRightRoot=extrprintf("Before blowing");emPossible;
        mpfi_set(rightRoot,extremPossible);
        mpfi_get_mid(extremPossibleFr, leftRoot,rightRoot);
        
        printf("We are in a valley and should search more to the left\n");
        printf("Left:=");
        printInterval(leftRoot);
        printf("Middle:=");
        printInterval(extremPossible);
        printf("Right=");
        printInterval(rightRoot);
        
      }
      //order: wleft>wmiddle, wmiddle>wright
     else if   (((mpfi_cmp(wleft,wmiddle))>0) && ((mpfi_cmp(wmiddle,wright))>0)  ){
        //we are in a valley and should search more to the right
        //newLeftRoot=extremPossible;
        //newRightRoot=oldRightRoot;
        mpfi_set(leftRoot,extremPossible);
        mpfi_get_mid(extremPossibleFr, leftRoot,rightRoot);
        printf("We are in a valley and should search more to the right\n");
        printf("Left:=");
        printInterval(leftRoot);
        printf("Middle:=");
        printInterval(extremPossible);
        printf("Right=");
        printInterval(rightRoot);
      }
      else{ printf("We are in a valley but looks like a hill....\n");
      printf("Left:=");
        printInterval(extremLeftInt);
     
     printf("LeftValue:=");
        printInterval(wleft);
     
        printf("Middle:=");
        printInterval(extremPossible);
        
        printf("MiddleValue:=");
        printInterval(wmiddle);
        
        
        printf("Right=");
        printInterval(extremRightInt);
      
      printf("RightValue=");
        printInterval(wright);
      
      }
    }
    else{ //if we are in a hill, three are 3 posibilities
      //order: wleft<wmiddle, wmiddle>wright
      if (((mpfi_cmp(wleft,wmiddle))<0) && ((mpfi_cmp(wmiddle,wright))>0)  ){
      //we have found a max
      printf("We have found a max:");
      printInterval(extremInt);
      notFound=0;
      }
      //order: wleft<wmiddle, wmiddle<wright
      if   ( ((mpfi_cmp(wleft,wmiddle))<0) && ((mpfi_cmp(wmiddle,wright))<0)  ){
        //we are in a hill and should search more to the right
        //newLeftRoot=extremPossible;
        //newRightRoot=oldRightRoot;
         mpfi_set(leftRoot,extremPossible);
         mpfi_get_mid(extremPossibleFr, leftRoot,rightRoot);
         printf("We are in a hill and should search more to the right\n");
        printf("Left:=");
        printInterval(leftRoot);
        printf("Middle:=");
        printInterval(extremPossible);
        printf("Right=");
        printInterval(rightRoot);
       }
      //order: wleft>wmiddle, wmiddle>wright
      if   (((mpfi_cmp(wleft,wmiddle))>0) && ((mpfi_cmp(wmiddle,wright))>0)  ){
        //we are in a hill and should search more to the left        
        //newLeftRoot=oldLeftRoot;
        //newRightRoot=extremPossible;
        mpfi_set(rightRoot,extremPossible);
        mpfi_get_mid(extremPossibleFr, leftRoot,rightRoot);
        printf("We are in a hill and should search more to the left\n");
        printf("Left:=");
        printInterval(leftRoot);
        printf("Middle:=");
        printInterval(extremPossible);
        printf("Right=");
        printInterval(rightRoot);
      }
    }
      /*printf("\nThe extremum is between:");
      printInterval(leftRoot);
      printf("-----");
      printInterval(extremPossible);
      printf("-----");
      printInterval(rightRoot);*/
    safeMalloc((freeDegrees+1)*sizeof(mpfr_t));
    }
printf("End of While loop");    
/*printf("\nThe extremum is between:");
      printInterval(leftRoot);
      printf("-----");
      printInterval(extremPossible);
      printf("-----");
      printInterval(rightRoot); 
*/      
     mpfi_set(wExtrem,extremInt); 
     
     mpfi_clear(leftRoot);
     mpfi_clear(rightRoot);
     mpfi_clear(extremPossible);
    mpfi_clear(extremInt);
    mpfi_clear(extremLeftInt);
    mpfi_clear(extremRightInt);
    mpfr_clear(extremPossibleFr);
    mpfr_clear(extremRight);
    mpfr_clear(extremLeft);
    mpfi_clear(wleft);
    mpfi_clear(wmiddle);
    mpfi_clear(wright);  
return hill;        
}




 
   
node* interpolation (mpfr_t *coeffArray, mpfi_t *wbound, node *f, mpfi_t x, int n) {
  node *res, *dif, *w;
  chain *zeros;
  int nrRoots;
  mpfr_t *chebArray, *chebArrayExtrema, *M, *b, *lambdai_vect, *roots, *newRoots;
  mpfi_t *rootsIntervals;//, *wArrayExtrema;
  mpfr_t boundPos, boundNeg;
  mpfr_t u, v;
  mpfr_t  var1, var2, var3, zero_mpfr; 
  int freeDegrees=n;
  int j,i,r;
  int points=getToolPoints();
  int notEnoughRoots,iteration,oldPoints;
/*************************************************************/
  int verbosity=10;
  // Initialisations and precomputations
  
  mpfr_init2(u, getToolPrecision());
  mpfr_init2(v, getToolPrecision());
 
  mpfr_init2(zero_mpfr,53);
  mpfr_set_d(zero_mpfr,0.,GMP_RNDN);
 
  mpfi_get_left(u,x);
  mpfi_get_right(v,x);
  
  
  mpfr_init2(var1, getToolPrecision());
  mpfr_init2(var2, getToolPrecision());
  mpfr_init2(var3, getToolPrecision());
  
  chebArray = safeMalloc((freeDegrees+1)*sizeof(mpfr_t));
  for(j=1; j <=freeDegrees+1 ; j++) {
    mpfr_init2(chebArray[j-1], getToolPrecision());
  }
  
  
  chebArrayExtrema = safeMalloc((freeDegrees)*sizeof(mpfr_t));
  for(j=1; j <=freeDegrees ; j++) {
    mpfr_init2(chebArrayExtrema[j-1], getToolPrecision());
  }
  
  /*Chebyshev points where the extremas are likely to be found*/ 

  mpfr_const_pi(var1, GMP_RNDN);
  mpfr_div_si(var1, var1, (long)(freeDegrees+1), GMP_RNDN); // var1 = Pi/freeDegrees
  mpfr_sub(var2, u, v, GMP_RNDN);
  mpfr_div_2ui(var2, var2, 1, GMP_RNDN); // var2 = (u-v)/2
  mpfr_add(var3, u, v, GMP_RNDN);
  mpfr_div_2ui(var3, var3, 1, GMP_RNDN); // var3 = (u+v)/2
  
  
  
  for (i=1 ; i <= freeDegrees ; i++) {
    mpfr_mul_si(chebArrayExtrema[i-1], var1, i, GMP_RNDN);
    mpfr_cos(chebArrayExtrema[i-1], chebArrayExtrema[i-1], GMP_RNDN);
    mpfr_fma(chebArrayExtrema[i-1], chebArrayExtrema[i-1], var2, var3, GMP_RNDN); // x_i = [cos((i)*Pi/(freeDegrees+1))]*(u-v)/2 + (u+v)/2
  }
  
  
  
  
  /*************************************************************/
  /*                  Alternative Cheb points                  */
  mpfr_const_pi(var1, GMP_RNDN);
  mpfr_div_si(var1, var1, 2*((long)freeDegrees+1), GMP_RNDN); // var1 = Pi/(2*freeDegrees+2)
  mpfr_sub(var2, u, v, GMP_RNDN);
  mpfr_div_2ui(var2, var2, 1, GMP_RNDN); // var2 = (u-v)/2
  mpfr_add(var3, u, v, GMP_RNDN);
  mpfr_div_2ui(var3, var3, 1, GMP_RNDN); // var3 = (u+v)/2

  for (i=1 ; i <= freeDegrees+1 ; i++) {
    mpfr_mul_si(chebArray[i-1], var1, 2*i-1, GMP_RNDN);
    mpfr_cos(chebArray[i-1], chebArray[i-1], GMP_RNDN);
    mpfr_fma(chebArray[i-1], chebArray[i-1], var2, var3, GMP_RNDN); // x_i=[cos((2i-1)*Pi/(2freeDegrees+2))]*(u-v)/2 + (u+v)/2
  }
  /*************************************************************/
  

  if(verbosity>=8) {
    printf("Computed points set:\n");
    for(i=1;i<=freeDegrees+1;i++) printMpfr(chebArray[i-1]);
  }
  
  perturbPoints(chebArray, freeDegrees+1);

  if(verbosity>=8) {
    printf("Computed points set:\n");
    for(i=1;i<=freeDegrees+1;i++) printMpfr(chebArray[i-1]);
  }
/*************************************************************/

//prepare the data for solving the system
  M = safeMalloc((freeDegrees+1)*(freeDegrees+1)*sizeof(mpfr_t));
  b = safeMalloc((freeDegrees+1)*sizeof(mpfr_t));
  lambdai_vect = safeMalloc((freeDegrees+1)*sizeof(mpfr_t));

  for(i=1; i <= freeDegrees+1 ; i++) {
    for(j=1; j<= freeDegrees+1; j++) {
      mpfr_init2(M[coeff(i,j,freeDegrees+1)],getToolPrecision());
      mpfr_pow_ui(M[coeff(i,j,freeDegrees+1)],chebArray[i-1],j-1,GMP_RNDN);
      if(verbosity>=8) {
        printMpfr(M[coeff(i,j,freeDegrees+1)]);
      }
    }
  }
  
    
  for(i=1; i<= freeDegrees+1; i++) {
  mpfr_set_ui(M[coeff(i,1,freeDegrees+1)],1, GMP_RNDN);
  }
  
  
  for(i=1;i<=freeDegrees+1;i++){
    mpfr_init2(b[i-1], getToolPrecision());
    mpfr_init2(lambdai_vect[i-1], getToolPrecision());
  }

 for (i=1 ; i <= freeDegrees+1 ; i++) {
	r = evaluateFaithfulWithCutOffFast(var1, f, NULL, chebArray[i-1], zero_mpfr, getToolPrecision()); // var1=f(x_i)
	if(r==0) mpfr_set_d(var1, 0., GMP_RNDN);
	mpfr_set(b[i-1],var1,GMP_RNDN);
 }
  if(verbosity>=8) {
    printf("The computed matrix is "); printMatrix(M, freeDegrees+1);
  
    printf("The values for the function in the chosen points are:\n");
    for(i=1;i<=freeDegrees+1;i++) printMpfr(b[i-1]);
  }
  //solve the system M*Lambda=b
  system_solve( lambdai_vect , M, b, freeDegrees+1, getToolPrecision());


  //get coeffs of polynomial
  
  if(verbosity>=8) {
    printf("The coeffs of the polynomial are:\n");
    for(i=1;i<=freeDegrees+1;i++) printMpfr(lambdai_vect[i-1]);
  }
  res = constructPoly(lambdai_vect,freeDegrees);
  for(i=1;i<=freeDegrees+1;i++) mpfr_set(coeffArray[i-1],lambdai_vect[i-1],GMP_RNDN);
  
  
  
  
  dif = (node*)safeMalloc(sizeof(node));
  dif->nodeType = SUB;
  dif->child1=f;
  dif->child2=res;
  notEnoughRoots=1;
  iteration=0;
  oldPoints=points;
  while((notEnoughRoots)&&(iteration<10)){
    zeros =uncertifiedFindZeros(dif, u, v, points, getToolPrecision());
    nrRoots=lengthChain(zeros);
    roots= (mpfr_t *)safeMalloc((nrRoots)*sizeof(mpfr_t));
    newRoots= (mpfr_t *)safeMalloc((nrRoots)*sizeof(mpfr_t));
  
    for(i=0;i<nrRoots;i++){
      mpfr_init2(roots[i],getToolPrecision());
      mpfr_set(roots[i],*((mpfr_t*)zeros->value),GMP_RNDN);
      zeros=zeros->next;
    }
  
    rootsIntervals= (mpfi_t *)safeMalloc((nrRoots)*sizeof(mpfi_t));
   printf("We have found %d and we had %d cheby points:",nrRoots,freeDegrees+1); 
  
    if (nrRoots>(freeDegrees+1)) {
      
      printf("\n-------we have more cheby points than needed----------\n");
      notEnoughRoots=0;
      findNewSetRoots(newRoots, roots, nrRoots, chebArray, freeDegrees+1);
      w=constructPolyFromRoots(newRoots,freeDegrees+1);
    
    for(i=0;i<freeDegrees+1;i++){
    mpfi_init2(rootsIntervals[i],getToolPrecision());
   // printf("we are going for the interval root\n");
    if (0==(getIntervalAroundRoot(&rootsIntervals[i],dif, newRoots[i], getToolPrecision()))){printf("Huge Warning one root can not be made as an interval!!!!!!!!!!");}
    //printInterval(rootsIntervals[i]);
    }
   
  } 
  else if (nrRoots==(freeDegrees+1)){
  printf("\n-------we have the good nr of cheby points----------\n");
  notEnoughRoots=0;
  w=constructPolyFromRoots(roots,freeDegrees+1);
  for(i=0;i<freeDegrees+1;i++){
    mpfi_init2(rootsIntervals[i],getToolPrecision());
    if (0==(getIntervalAroundRoot(&rootsIntervals[i],dif, roots[i], getToolPrecision()))){printf("Huge Warning one root can not be made as an interval!!!!!!!!!!");}
    //printInterval(rootsIntervals[i]);
    }
  
  }
  else{
  printf("\n-------we don't have enough cheby points----------\n");
  points=2*points;
  iteration++;
  }
  
  }
  points=oldPoints;
  printf("Roots:\n");
  for(i=0;i<freeDegrees+1;i++){
    printInterval(rootsIntervals[i]);
    
  }
  
  printf("Extremas:\n");
  for(i=0;i<freeDegrees;i++){
    printMpfr(chebArrayExtrema[i]);
    
  }
  printf("W is");
  printTree(w);
  mpfr_init2(boundPos, getToolPrecision());
  mpfr_init2(boundNeg, getToolPrecision());
  uncertifiedInfnorm(boundPos, w, u,v, getToolPoints(), getToolPrecision());
  printMpfr(boundPos);
  mpfr_neg(boundNeg, boundPos, GMP_RNDN);
  mpfi_interv_fr(*wbound, boundNeg, boundPos);
  printInterval(*wbound);
  //mpfr_clear(boundNeg);
  //mpfr_clear(boundPos);
  /*
  wArrayExtrema= (mpfi_t *)safeMalloc((freeDegrees)*sizeof(mpfi_t));
  for(i=0;i<freeDegrees;i++){
    mpfi_init2(wArrayExtrema[i], getToolPrecision());
    printf("\n******************************Start of computations for the %d extremum*******************************************\n",i);
    findWExtremum(wArrayExtrema[i],i,chebArrayExtrema[i],freeDegrees,rootsIntervals);
    printf("\n******************************End of computations for the %d extremum*******************************************\n",i);    
    printf("We have successfully found %d extremum:\n", i);
    printInterval(wArrayExtrema[i]);
    
  } */
  
  
  for(i=1; i <= freeDegrees+1 ; i++) {
    for(j=1; j<= freeDegrees+1; j++) {
      mpfr_clear(M[coeff(i,j,freeDegrees+1)]);
    }
  }
  for(i=1;i<=freeDegrees+1;i++){
    mpfr_clear(b[i-1]);
    mpfr_clear(lambdai_vect[i-1]);
    mpfr_clear(chebArray[i-1]);
  }
  
  free(M);
  free(b);
  free(lambdai_vect);
  free(chebArray);
  
  mpfr_clears(u, v, var1, var2, var3, zero_mpfr,(mpfr_ptr) 0);
  printf("\n out of interpolation\n");
  return res; 

//:)


}

/***************************************************/


void system_solve(mpfr_t *res, mpfr_t *M, mpfr_t *b, int n, mp_prec_t prec) {
  chain *i_list=NULL;
  chain *j_list=NULL;
  chain *curri;
  chain *currj;
  int i0, j0, i, j, k;
  int *var;
  mpfr_t max,lambda;
  int *order_i = safeMalloc(n*sizeof(int));
  int *order_j = safeMalloc(n*sizeof(int));

  mpfr_init2(max, 53);
  mpfr_init2(lambda, prec);

  for(i=1;i<=n;i++) {
    var = safeMalloc(sizeof(int));
    *var = i;
    i_list = addElement(i_list, (void *)var);
  }
  for(j=1;j<=n;j++) {
    var = safeMalloc(sizeof(int));
    *var = j;
    j_list = addElement(j_list, (void *)var);
  }


  // Triangulation by Gaussian elimination
  i0 = j0 = -1;
  for(k=1;k<=n;k++) {
    mpfr_set_d(max, 0., GMP_RNDN); //exact

    // In this part, we search for the bigger element of the matrix
    curri = i_list;
    while(curri!=NULL) {
      currj = j_list;
      while(currj!=NULL) {
	i = *(int *)(curri->value);
	j = *(int *)(currj->value);
	if(mpfr_cmpabs(M[coeff(i,j,n)],max)>=0) {
	  i0 = i;
	  j0 = j;
	  mpfr_set(max, M[coeff(i,j,n)], GMP_RNDN);
	}
	currj = currj->next;
      }
      curri = curri->next;
    }
    
    i_list = removeInt(i_list, i0);
    j_list = removeInt(j_list, j0);

    order_i[k-1] = i0;
    order_j[k-1] = j0;

    // Here we update the matrix and the second member
    curri = i_list;
    while(curri!=NULL) {
      i = *(int *)(curri->value);
      mpfr_div(lambda, M[coeff(i,j0,n)], M[coeff(i0,j0,n)], GMP_RNDN);
      mpfr_neg(lambda, lambda, GMP_RNDN);

      currj = j_list;
      while(currj!=NULL) {
	j = *(int *)(currj->value);
	mpfr_fma(M[coeff(i,j,n)], lambda, M[coeff(i0,j,n)], M[coeff(i,j,n)], GMP_RNDN);
	currj = currj->next;
      }

      mpfr_fma(b[i-1], lambda, b[i0-1], b[i-1], GMP_RNDN);
      mpfr_set_d(M[coeff(i,j0,n)], 0., GMP_RNDN); // this line is not useful strictly speaking
      curri = curri->next;
    }
  }
  


  // Resolution of the system itself
  for(i=1;i<=n;i++) {
    var = safeMalloc(sizeof(int));
    *var = i;
    i_list = addElement(i_list, (void *)var);
  }

  for(k=n;k>=1;k--) {
    i0 = order_i[k-1];
    j0 = order_j[k-1];
    mpfr_div(res[j0-1], b[i0-1], M[coeff(i0,j0,n)], GMP_RNDN);
    i_list = removeInt(i_list, i0);

    curri = i_list;
    while(curri!=NULL) {
      i = *(int *)(curri->value);
      mpfr_neg(M[coeff(i,j0,n)], M[coeff(i,j0,n)], GMP_RNDN);
      mpfr_fma(b[i-1], M[coeff(i,j0,n)], res[j0-1], b[i-1], GMP_RNDN);
      curri=curri->next;
    }
  }

  free(order_i);
  free(order_j);
  freeChain(i_list, freeIntPtr);
  freeChain(j_list, freeIntPtr);
  mpfr_clear(max);
  mpfr_clear(lambda);
  return;
}


/*********************************************************************/

void myPrintValue(mpfr_t *x, mp_prec_t prec) {
  mpfr_t y;
  mpfr_init2(y,prec);
  mpfr_set(y,*x,GMP_RNDN);
  printValue(&y);
  mpfr_clear(y);
}

/*********************************************************************/

void printMatrix(mpfr_t *M, int n) {
  int i,j;
  printf("[");
  for(i=1;i<=n;i++) {
    for(j=1;j<=n;j++) {
      myPrintValue(&M[coeff(i,j,n)],53); if(j!=n)printf(", ");
    }
    if(i!=n) printf(";\n");
  }
  printf("]\n");
  return;
}
/*********************************************************************/

node *constructPoly(mpfr_t *coeff, int n) {
  int i=1;
  //chain *curr;
  node *temp1;
  node *temp2;
  node *temp3;
  node *temp4;
  node *temp5;
  node *temp6;

  node *poly;
  mpfr_t *ptr;

  poly = (node*)safeMalloc(sizeof(node));
  poly->nodeType = CONSTANT;
  ptr = (mpfr_t*)safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*ptr, getToolPrecision());
  mpfr_set(*ptr,coeff[0], GMP_RNDN);
  poly->value = ptr;

  
  for (i=1;i<=n;i++) {
    temp1 = (node*)safeMalloc(sizeof(node));
    temp1->nodeType = ADD;

    temp2 =(node*) safeMalloc(sizeof(node));
    temp2->nodeType = MUL;

    temp3 = (node*)safeMalloc(sizeof(node));
    temp3->nodeType = CONSTANT;
    ptr = (mpfr_t*)safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*ptr, getToolPrecision());
    mpfr_set(*ptr, coeff[i], GMP_RNDN);
    temp3->value = ptr;

    temp4 = (node*)safeMalloc(sizeof(node));
    temp4->nodeType = POW;

    temp5 = (node*)safeMalloc(sizeof(node));
    temp5->nodeType = VARIABLE;

    temp6 = (node*)safeMalloc(sizeof(node));
    temp6->nodeType = CONSTANT;
    ptr = (mpfr_t*)safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*ptr, getToolPrecision());
    mpfr_set_si(*ptr, i, GMP_RNDN);
    temp6->value = ptr;

    temp4->child1 = temp5;
    temp4->child2 = temp6;
    temp2->child1 = temp3;
    temp2->child2 = temp4;
    temp1->child1 = temp2;
    temp1->child2 = poly;
    poly = temp1;
   }

  return poly;
}
/*********************************************************************/



/*********************************************************************/

node *constructPolyFromRoots(mpfr_t *roots, int n) {
  int i=1;
  //chain *curr;
  node *temp1;
  node *temp2;
  node *temp3;
  node *temp4;
 
  node *poly;
  mpfr_t *ptr;

  
    poly= (node*)safeMalloc(sizeof(node));
    poly->nodeType = CONSTANT;
    ptr = (mpfr_t*)safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*ptr, getToolPrecision());
    mpfr_set_ui(*ptr, 1, GMP_RNDN);
    poly->value = ptr;

   for (i=0;i<n;i++) {
   
    temp1 = (node*)safeMalloc(sizeof(node));
    temp1->nodeType = SUB;

    temp2 =(node*) safeMalloc(sizeof(node));
    temp2->nodeType = VARIABLE;

    temp3 = (node*)safeMalloc(sizeof(node));
    temp3->nodeType = CONSTANT;
    ptr = (mpfr_t*)safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*ptr, getToolPrecision());
    mpfr_set(*ptr, roots[i], GMP_RNDN);
    temp3->value = ptr;

    temp4 = (node*)safeMalloc(sizeof(node));
    temp4->nodeType = MUL;
  
    temp1->child1 = temp2;
    temp1->child2 = temp3;

    temp4->child1=temp1;
    temp4->child2=poly;
    
    poly = temp4;
   }

  return poly;
}
/*********************************************************************/


/* Random pertubration of the points... */
/* There is p points : x0 ... x(p-1)    */
void perturbPoints(mpfr_t *x, int p) {
  mpfr_t perturb;
  mpfr_t var1, var2, var3;
  int i;
  gmp_randstate_t random_state;
   printf("ajunge aici");
  gmp_randinit_default(random_state);
  gmp_randseed_ui(random_state, 65845285);

  mpfr_init2(perturb, getToolPrecision());
  mpfr_init2(var1, getToolPrecision());
  mpfr_init2(var2, getToolPrecision());
  mpfr_init2(var3, getToolPrecision());
  printf("ajunge aici");
  for(i=1;i<=p-2;i++) {
    mpfr_urandomb(perturb, random_state);
    mpfr_mul_2ui(perturb, perturb, 1, GMP_RNDN);
    mpfr_sub_ui(perturb, perturb, 1, GMP_RNDN);
    mpfr_div_2ui(perturb, perturb, 2, GMP_RNDN); // perturb \in [-1/4; 1/4]

    mpfr_sub(var1, x[i], x[i-1], GMP_RNDN);
    mpfr_sub(var2, x[i+1], x[i], GMP_RNDN);
    if (mpfr_cmpabs(var1,var2)>0) mpfr_mul(var3, var2, perturb, GMP_RNDN);
    else mpfr_mul(var3, var1, perturb, GMP_RNDN);
    mpfr_add(x[i], x[i], var3, GMP_RNDN);
  }
   printf("ajunge aici");
  gmp_randclear(random_state);
  mpfr_clear(perturb);
  mpfr_clear(var1);
  mpfr_clear(var2);
  mpfr_clear(var3);
   printf("ajunge aici");
  return;
}
/*********************************************************************/



/*
void constructCoeffsPolyFromRoots(mpfi_t *coeffs, mpfi_t *roots, int n){

  int N, i,j,k;
  mpfi_t *coeffsW;
  int *indexes;
  
  
  
  
  
  indexes=safeMalloc(*sizeof(int));
  coeffsW =safeMalloc((n+1)*sizeof(mpfi_t));
  partialProd=safeMalloc((n+1)*sizeof(mpfi_t));
  for (i=0;i<n ; i++){  
    mpfi_init2(coeffsW[i],getToolPrecision());
    mpfi_set_ui( coeffsW[i], 0);
  }
  mpfi_init2(coeffsW[n],getToolPrecision());
  mpfi_set_ui( coeffsW[n], 1);
    
    
    
    
}



 

getCoeffsPolyFromRoots2 = proc(listRoots,n) {
  var N,i,j,k,comb,elem, newelem, partialProd, newProd,coeffW,pp;
  
  
  coeffW=[||];
  
  for i from 0 to n-1 do
  coeffW[i]=[0];
  coeffW[n]=[1];
  
    
  comb=[||];
  comb[0]=[||];
  comb[1]=[||];
  
  partialProd=[||];
  partialProd[0]=[||];
  partialProd[1]=[||];
  
  
  for i from 0 to n-1 do{
    comb[1]=[|i|].:comb[1];
    partialProd[1]=listRoots[i].:partialProd[1];
    coeffW[n-1]=coeffW[n-1]+listRoots[i];
  };
  coeffW[n-1]=coeffW[n-1]*(-1)^(n+n-1);
  //print(comb);
  
  for k from 1 to (n-1) do{ 
    comb[2]=[||];
    partialProd[2]=[||];
    
    index=0;
    for elem in comb[1] do {
          
      if (elem[0]<(n-1)) then {
        
        for j from (elem[0]+1) to (n-1) do {
          newelem=elem;
          newProd=partialProd[1][index];
          
          newelem=j.:newelem;
          newProd=newProd*listRoots[j];
          
          comb[2]=newelem.:comb[2];
          partialProd[2]=newProd.:partialProd[2];
          
          
        };
      };
      index=index+1;
    };
    for pp in partialProd[2] do{
    coeffW[n-k-1]=coeffW[n-k-1]+pp;
    };
    coeffW[n-k-1]=coeffW[n-k-1]*(-1)^(n+n-k-1);
    partialProd[1]=partialProd[2];
    comb[1]=comb[2];
  };
  //print(comb);
  //print(partialProd);
  return coeffW;
};  

*/


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/


/* This function performs the differentiation.
   See the commentaries below.
*/
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



/*Code name: ana are mere*/

/*********************************************************************/
/***********************************************************/
/***********************************************************/

/*cheby model structure:
n- order: polynomial of degree n-1, remainder of order o(x^n)
rem_bound - bound for the remainder
poly_array - array of coeffs for the polynomial (mpfi's)
poly_bound - bound for the polynomial (helpful for computations)
x- interval on which the tm is computed
x0 - interval around the translation point
*/
typedef struct cmdl {
int n; 
mpfi_t rem_bound;
mpfi_t *poly_array;
mpfi_t poly_bound;
mpfi_t x;
mpfi_t x0; //is used just for translation

} tModel;
void polynomialTranslateCM(mpfi_t *coeffsT,int n,mpfi_t *coeffs,mpfi_t x0);  
void ctMultiplication_CM(tModel*d,tModel*s, mpfi_t c);
void ctDivision_CM(tModel*d,tModel*s, mpfi_t c);
void polynomialBoundHorner(mpfi_t *bound,int n,mpfi_t *coeffs,mpfi_t x0,mpfi_t x);
void polynomialBoundSharp(mpfi_t *bound,int n,mpfi_t *coeffs,mpfi_t x0Int,mpfi_t x);
void polynomialBoundSharp(mpfi_t *bound,int n,mpfi_t *coeffs,mpfi_t x0,mpfi_t x);
void polynomial_CM(tModel *t,mpfi_t *r, int d, int n, mpfi_t x0, mpfi_t x);
/*This function creates an empty chebyshev model
*/
tModel* createEmptycModel(int n,  mpfi_t x0, mpfi_t x){
  tModel* t;
  int i;
 
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
  }
  return t;
}
/*the convention for all the following functions is:
the tmodel given as parameter has to be created previously 
*/

/*This function sets the cheby model t 
with constant ct;
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

/*This function dealocates a cheby model
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

/*This function pretty prints a cheby model
*/
void printcModel(tModel *t){
  int i;
  printf("\nCheby model of order, %d in base x - ", t->n);
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


/*This function sets a cheby model t 
with the values given by anoter tm tt
they implicitely have the same base
and interval

*/
void copycModel(tModel *t, tModel *tt){
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

/*void evaluateMpfiFunction(mpfi_t y, node *f, mpfi_t x, mp_prec_t prec){
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
 
}*/

/*This function transforms a polynomial with interval coeffs
into a poly with mpfr coeffs and a small remainder
Parameters:
--input: n  - degree of poly
         gc - array of given coeffs
         x0 - expasion point
         x  - interval
-- output: rc -mpfr coeffs
           errors - errors around the coeffs rc
           rest - remainder  
*/
void mpfr_get_poly(mpfr_t *rc, mpfi_t *errors, mpfi_t rest, int n, mpfi_t *gc, mpfi_t x0, mpfi_t x){
  int i;
  mpfi_t *res;
  mpfi_t r;
  res= (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
  mpfi_init2(r,getToolPrecision());
  for (i=0; i<=n; i++){
    mpfi_mid(rc[i],gc[i]);
    mpfi_init2(res[i], getToolPrecision());
    mpfi_sub_fr(res[i],gc[i], rc[i]);
    mpfi_set(errors[i],res[i]);
  }
  polynomialBoundHorner(&r,n,res,x0,x);
  
  mpfi_set(rest,r);
  
  for (i=0; i<=n; i++){
    mpfi_clear(res[i]);  
  }
  free(res);
  mpfi_clear(r);
  return;

}

/*This function transforms a polynomial with interval coeffs
into a poly with mpfr coeffs and a small remainder
Parameters:
--input: n  - degree of poly
         gc - array of given coeffs
         x0 - expasion point
         x  - interval
-- output: rc -mpfr coeffs
           errors - errors around the coeffs rc
           rest - remainder  
*/
void mpfr_get_poly2(mpfr_t *rc, mpfi_t *rest, int n, mpfi_t *gc, mpfi_t x0, mpfi_t x){
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
  polynomialBoundHorner(&r,n,res,x0,x);
  
  mpfi_set(*rest,r);
  
  for (i=0; i<=n; i++){
    mpfi_clear(res[i]);  
  }
  free(res);
  mpfi_clear(r);
  return;

}


/*This function computes the tm for multiplication of two 
given cm's 
*/
void  multiplication_CM(tModel *t,tModel *c1, tModel *c2){
  //we will multiply two cheby models of order n; and obtain a new cheby model of order n;
  int n,i,j;
  mpfi_t *r;
  tModel *tt, *ttt;
  mpfi_t temp1, temp2;
  
  n=t->n;
  
  //aux tm for doing the multiplications
  tt=createEmptycModel(n,t->x0, t->x);
  for(i=0;i<=n-1;i++){
   mpfi_set_ui(tt->poly_array[i],0);
  }
  
  mpfi_init2(temp1, getToolPrecision());
  mpfi_init2(temp2, getToolPrecision());
  
  
  /*absolute error*/
  /*We are multiplying cheby models, considering the absolute error
    We are given:  (C1,delta1)
                   (C2,delta2)
               
    The product is:(C1*C2-> intepolated to [x_0...x_n-1], {delta1*B(C2)+delta2*B(C1)+delta1*delta2+B(interp(C1,C2))}*/
 
   /*compute in temp1 delta1*delta2*/
   mpfi_mul(temp1, c1->rem_bound, c2->rem_bound);
   //printf("compute in temp2 delta2*B(T1)");
   mpfi_mul(temp2, c1->poly_bound, c2->rem_bound);
  
   mpfi_add(tt->rem_bound,temp1,temp2);
  
   /*compute in temp1 delta1*B(T2)*/
   mpfi_mul(temp1, c1->rem_bound, c2->poly_bound);
   mpfi_add(tt->rem_bound, tt->rem_bound, temp1);
  
   /*compute the product of the two polynomials*/
   
   r= (mpfi_t *)safeMalloc((2*n-1)*sizeof(mpfi_t));
   for(i=0;i<=2*n-2;i++){
    mpfi_init2(r[i],getToolPrecision());
    mpfi_set_ui(r[i],0);
   }
   
   for(i=0; i<n;i++)
     for (j=0;j<n;j++){
       mpfi_mul(temp1,c1->poly_array[i], c2->poly_array[j]);
       //if ((i+j)<n )
       //  mpfi_add(tt->poly_array[i+j],tt->poly_array[i+j],temp1);
       //else
        mpfi_add(r[i+j],r[i+j],temp1);
     }
     /*in r we have the product of the two polynomials*/
    /*interpolate r*/
    ttt=createEmptycModel(n,t->x0, t->x);
    polynomial_CM(ttt,r,2*n-1, n,t->x0,t->x);
    printf("\nAfter multiplication - interpolation we have:\n");
    printcModel(ttt);
    //polynomialBoundSharp(&temp1, 2*n-2,r,t->x0,t->x);
    for(i=0; i<n;i++){
    mpfi_set(tt->poly_array[i],ttt->poly_array[i]);
    }
    mpfi_set(tt->poly_bound, ttt->poly_bound);
    
    /*we add temp1 in the bound for the remainder*/
    mpfi_add(tt->rem_bound,tt->rem_bound,ttt->rem_bound);
   
    /*we compute the new polynomial bound for the new model*/
    //polynomialBoundSharp(&temp1, n-1,tt->poly_array,t->x0,t->x);   
    //mpfi_set(tt->poly_bound,temp1);
    /*for(i=0;i<2*n-1;i++)
      mpfi_clear(r[i]);
    free(r); 
  */
    
  mpfi_clear(temp1);
  mpfi_clear(temp2);
  
 
  
  /*set the result*/
  copycModel(t,tt);
  /*clear the aux tm*/
  cleartModel(tt);
  //cleartModel(ttt);
    
    printf("out of multiplication");
    
    
 }






//-----------------------------------------------------------
//-----------------------------------------------------------

/*This function computes the tm for addition of two 
given tm's 
The addition of two taylor models is the same, regardless the mode, 
absolute or relative - I put the parameter just to have some coherence
with the other functions
*/
void addition_CM(tModel *t,tModel *child1_tm, tModel *child2_tm){
  int i;
  int n;
  tModel *tt;
  
  n=t->n;
  tt=createEmptycModel(n,t->x0,t->x);
  for(i=0;i<n;i++)  
  mpfi_add(tt->poly_array[i], child1_tm->poly_array[i],child2_tm->poly_array[i]);
  
  mpfi_add(tt->rem_bound,child1_tm->rem_bound,child2_tm->rem_bound);
  polynomialBoundSharp(&tt->poly_bound, n-1,tt->poly_array,t->x0,t->x);   
  copycModel(t,tt);
  cleartModel(tt);
}

void polynomial_CM(tModel *t,mpfi_t *r, int d, int n, mpfi_t x0, mpfi_t x){
  int i, trueDegree;
  tModel *tt;
  mpfi_t *nDeriv;//, *nDeriv2;
  //mpfi_t fact, temp,pow;
  node *f, *p;
  mpfr_t *coeffArray;
  mpfr_t *coeffPoly;
  mpfi_t rest;
  mpfr_t fac;
  mpfi_t temp;
  
  tt=createEmptycModel(n,x0,x);
  printcModel(tt);
  /*Use INTERPOLATION to compute the coeffs*/
 
  trueDegree=d; 
  while (mpfi_is_zero(r[trueDegree-1])>0)  {
  trueDegree--;
  if (trueDegree==0) break;}
  printf("The true degree of the polynomial is: %d", trueDegree-1);
  
 
  if(trueDegree<=n){
  printf("\n\nsmaller\n\n");
    for (i=0;i<trueDegree;i++){
    mpfi_set(tt->poly_array[i], r[i]);
    }
    for (i=trueDegree;i<n;i++){
    mpfi_set_ui(tt->poly_array[i], 0);
    }
    
    mpfi_set_ui(tt->rem_bound,0);
    
    /*bound the polynomial obtained*/
    
     polynomialBoundSharp(&tt->poly_bound, n-1,tt->poly_array,tt->x0,tt->x);   
  }
  else{
  printf("trueDegree=  %d, n= %d ", trueDegree,n);
  
  coeffPoly=safeMalloc((trueDegree)*sizeof(mpfr_t));
  for (i=0;i<trueDegree;i++){
  mpfr_init2(coeffPoly[i], getToolPrecision());
  
  }
  mpfi_init2(rest, getToolPrecision());
  /*create a polynomial*/
  printf("before mpfr_get_poly");
  mpfr_get_poly2(coeffPoly, &rest, trueDegree-1, r, x0,x);
  printf("after mpfr_get_poly");  
  f=constructPoly(coeffPoly, trueDegree-1);
  
  printf("\n in cheby models for polynomial: ");
  printTree(f);
  printf("\n");
  
  coeffArray=safeMalloc((n)*sizeof(mpfr_t));
  for (i=0;i<n;i++){
  mpfr_init2(coeffArray[i], getToolPrecision());
  }
  printInterval(x);
  mpfi_init2(temp, getToolPrecision());
  p=interpolation(coeffArray,&temp, f,x,n-1);
  if (p==NULL) {printf("Interpolation did not work!!!");return;}
  /*Use AD for the base functions to bound derivatives up to nth derivative*/
  nDeriv= (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    for(i=0;i<=n;i++){
      mpfi_init2(nDeriv[i], getToolPrecision());
    }
  
  auto_diff(nDeriv,f,x,n);
  printf("zzzzzzzzzzzzzzzzzzz");
  mpfi_set(tt->rem_bound, nDeriv[n]);
  mpfr_init2(fac, getToolPrecision());
  mpfr_fac_ui(fac,n,GMP_RNDN);
  mpfi_div_fr(tt->rem_bound, tt->rem_bound, fac);
  mpfi_mul(tt->rem_bound,tt->rem_bound,temp);
  mpfr_clear(fac);
  
  for (i=0;i<n;i++){
  mpfi_set_fr(tt->poly_array[i], coeffArray[i]);
  printInterval(tt->poly_array[i]);
  }
  printf("zzzzzzzzzzzzzzzzzzz");
  
  //polynomialTranslateCM(tt->poly_array,n-1,tt->poly_array,tt->x0);
  
  /*bound the polynomial obtained*/
    
  polynomialBoundSharp(&tt->poly_bound, n-1,tt->poly_array,tt->x0,tt->x);   
  
  printf("zzzzzzzzzzzzzzzzzzz");
  printInterval(tt->poly_bound);
  /*
  for(i=0;i<=n;i++){
    mpfi_clear(nDeriv[i]);
  }
  free(nDeriv);  

  for(i=0;i<n;i++){
    mpfr_clear(coeffArray[i]);
  }
  free(coeffArray);  
  for(i=0;i<d;i++){
    mpfr_clear(coeffPoly[i]);
  }
  free(coeffPoly);  */
}

copycModel(t,tt);
  cleartModel(tt);
  
printf("\n\nout\n\n");
}


void polynomialNode_CM(tModel *t,node*f, int n, mpfi_t x0, mpfi_t x){
  int i;
  tModel *tt;
  mpfi_t *nDeriv;//, *nDeriv2;
  //mpfi_t fact, temp,pow;
  node  *p;
  mpfr_t *coeffArray;
  
  
  mpfr_t fac;
  tt=createEmptycModel(n,x0,x);
  
  /*Use INTERPOLATION to compute the coeffs*/
 
  
  
  
  printf("\n in cheby models for polynomial: ");
  printTree(f);
  printf("\n");
  
  coeffArray=safeMalloc((n)*sizeof(mpfr_t));
  for (i=0;i<n;i++){
  mpfr_init2(coeffArray[i], getToolPrecision());
  }
  p=interpolation(coeffArray,&(tt->rem_bound), f,x,n-1);
  if (p==NULL) return;
  /*Use AD for the base functions to bound derivatives up to nth derivative*/
  nDeriv= (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    for(i=0;i<=n;i++){
      mpfi_init2(nDeriv[i], getToolPrecision());
    }
  
  auto_diff(nDeriv,f,x,n);
  mpfi_mul(tt->rem_bound,tt->rem_bound, nDeriv[n]);
  mpfr_init2(fac, getToolPrecision());
  mpfr_fac_ui(fac,n,GMP_RNDN);
  mpfi_div_fr(tt->rem_bound, tt->rem_bound, fac);
  mpfr_clear(fac);
  
  for (i=0;i<n;i++){
  mpfi_set_fr(tt->poly_array[i], coeffArray[i]);
  }
  polynomialTranslateCM(tt->poly_array,n-1,tt->poly_array,tt->x0);
  
  /*bound the polynomial obtained*/
    
  polynomialBoundSharp(&tt->poly_bound, n-1,tt->poly_array,tt->x0,tt->x);   
  
  copycModel(t,tt);
  cleartModel(tt);
  
  for(i=0;i<=n;i++){
    mpfi_clear(nDeriv[i]);
  }
  free(nDeriv);  

  for(i=0;i<n;i++){
    mpfr_clear(coeffArray[i]);
  }
  free(coeffArray);  
}




/*This function computes a chebyshev model - there is no mode x0 is the point where we want the model to be translated
for a base function*/
void base_CM(tModel *t,int nodeType, int n, mpfi_t x0, mpfi_t x){
  int i;
  tModel *tt;
  mpfi_t *nDeriv;//, *nDeriv2;
  //mpfi_t fact, temp,pow;
  node *f, *v, *p;
  mpfr_t *coeffArray;
  mpfr_t fac;
  tt=createEmptycModel(n,x0,x);
  
  /*Use INTERPOLATION to compute the coeffs*/
  
  /*create a base function*/
  f = (node*)safeMalloc(sizeof(node));
  f->nodeType = nodeType;
  
  v = (node*)safeMalloc(sizeof(node));
  v->nodeType = VARIABLE;
  f->child1=v;
  printf("\n in cheby models for base functions: ");
  printTree(f);
  printf("\n");
  
  coeffArray=safeMalloc((n)*sizeof(mpfr_t));
  for (i=0;i<n;i++){
  mpfr_init2(coeffArray[i], getToolPrecision());
  }
  p=interpolation(coeffArray,&(tt->rem_bound), f,x,n-1);
  if (p==NULL) return;
  /*Use AD for the base functions to bound derivatives up to nth derivative*/
  nDeriv= (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    for(i=0;i<=n;i++){
      mpfi_init2(nDeriv[i], getToolPrecision());
    }
  
  baseFunction_diff(nDeriv,nodeType,x,n);
  mpfi_mul(tt->rem_bound,tt->rem_bound, nDeriv[n]);
  mpfr_init2(fac, getToolPrecision());
  mpfr_fac_ui(fac,n,GMP_RNDN);
  mpfi_div_fr(tt->rem_bound, tt->rem_bound, fac);
  mpfr_clear(fac);
  
  for (i=0;i<n;i++){
  mpfi_set_fr(tt->poly_array[i], coeffArray[i]);
  }
  polynomialTranslateCM(tt->poly_array,n-1,tt->poly_array,tt->x0);
  
  /*bound the polynomial obtained*/
    
  polynomialBoundSharp(&tt->poly_bound, n-1,tt->poly_array,tt->x0,tt->x);   
  
  copycModel(t,tt);
  cleartModel(tt);
  
  for(i=0;i<=n;i++){
    mpfi_clear(nDeriv[i]);
  }
  free(nDeriv);  
}


/*composition:
  Assumptions:
We are given a taylor model for the function f in x0, over x, order n
and a taylor model for basic function g in f(x0) over y \superset range(f,x), order n
We obtain a taylor model for g(f(x)) in x0 over x.

Note: this is just like in the proofs, but it is easier to create the tm for g 
outside this function, such that a parameter of type NodeType is not needed here
*/

/*This is not true anymore!!!
for taylor models of order 2, we are given exaclty the same parameters, except that the expansion point for g is not f(x0) but mid (range(f,x))*/
void composition_CM(tModel *t,tModel *g, tModel *f){
  int i;
  int n;
  tModel *tt, *partial_tmul,*tinterm,*tmul ;
  
  n=f->n;
  printf("in Composition CM");
  /*create the taylor model for f(x)-f(x0): M1 in the proofs*/
  tinterm=createEmptycModel(n,f->x0,f->x);
  copycModel(tinterm,f);
  /*if we create 2nd order tms*/
  //if (CM2==1){
  /*mpfi_sub(tinterm->poly_array[0],tinterm->poly_array[0],g->x0);*/
  /* second possibility, put this difference in the remainder*/
  /*mpfi_sub(tinterm->poly_array[0],tinterm->poly_array[0],g->x0);
  mpfi_add(tinterm->rem_bound,tinterm->rem_bound,tinterm->poly_array[0]);
  mpfi_set_ui(tinterm->poly_array[0],0);*/
  
  //}
  //else{
  /*set the ct part of tinterm as 0*/
  //mpfi_set_ui(tinterm->poly_array[0],0);
  //}
  tt=createEmptycModel(n,f->x0,f->x);
  consttModel(tt,g->poly_array[0]);
    
  tmul=createEmptycModel(n,f->x0,f->x);
  copycModel(tmul,tinterm);  
  
  partial_tmul=createEmptycModel(n,f->x0,f->x);
  copycModel(partial_tmul,tinterm);  
  
  
  for (i=1;i<n;i++){
  ctMultiplication_CM(tinterm,partial_tmul,g->poly_array[i]);
  addition_CM(tt,tt,tinterm);
  multiplication_CM(partial_tmul,partial_tmul,tmul);
  }
  
  //in tt we have now the good polynomial;
  //we have to compute the bounds.
  //printcModel(g);
  //printcModel(partial_tmul);
  //printcModel(tt);
   
  //mpfi_mul(partial_tmul->rem_bound,partial_tmul->rem_bound,g->rem_bound);
  mpfi_add(tt->rem_bound,tt->rem_bound,g->rem_bound);
  polynomialBoundSharp(&tt->poly_bound, n-1,tt->poly_array,tt->x0,tt->x);   
  copycModel(t,tt);
  //printcModel(tt);
  cleartModel(tt);
  cleartModel(partial_tmul);
  cleartModel(tmul);
  cleartModel(tinterm);
}


/*This function computes the tm for division
with a ct term of a given tm
*/
void ctDivision_CM(tModel*d,tModel*s, mpfi_t c){
  int i;
  int n;
  tModel *tt;
  n=s->n;
  tt=createEmptycModel(n,s->x0,s->x);
  copycModel(tt,s);
  for(i=0;i<n;i++)  mpfi_div(tt->poly_array[i],tt->poly_array[i], c);
  mpfi_div(tt->rem_bound,tt->rem_bound,c);
  mpfi_div(tt->poly_bound,tt->poly_bound,c);
  copycModel(d,tt);
  cleartModel(tt);
}


/*This function computes the tm for multiplication
with a ct term of a given tm
*/
void ctMultiplication_CM(tModel*d,tModel*s, mpfi_t c){
  int i;
  int n;
  n=s->n;
  tModel *tt;
  tt=createEmptycModel(n,s->x0,s->x);
  copycModel(tt,s);
  for(i=0;i<n;i++)  mpfi_mul(tt->poly_array[i], tt->poly_array[i],c);
  mpfi_mul(tt->rem_bound,tt->rem_bound,c);
  mpfi_mul(tt->poly_bound,tt->poly_bound,c);
  copycModel(d,tt);
  cleartModel(tt);
}



/*This function computes an interval bound
for a polynomial given by coeffs and order n, 
on int x, using basic IA and Horner form
*/
void polynomialBoundHorner(mpfi_t *bound,int n,mpfi_t *coeffs,mpfi_t x0,mpfi_t x){
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

/*This function computes a sharp interval bound
for a polynomial given by coeffs and order n, 
on int x, using dirtyfindzeros - for the moment
this function is not certified:
**Sturm's nr of roots test is not integrated**
**The mpfr's obtained for the roots are not expanded
into small intervals arround the true root**
*/
void polynomialBoundSharpUncertified(mpfi_t *bound,int n,mpfi_t *coeffs,mpfi_t x0,mpfi_t x){
  int i;
  mpfi_t r,partialBound;
 
  mpfr_t *polyCoeffs;
  node *poly, *diff_poly;
  mpfr_t a, b;
  chain *zeros;
  mpfi_t *rootsIntervals;
  mpfr_t inf, sup, pinf, psup;
  mpfi_t extr1;
  mpfi_t shiftx;
  mpfi_t zero;
  mpfi_t err;
  int points, nrRoots;
  points =100;
  
  mpfi_init2(r,getToolPrecision());
  mpfi_init2(partialBound, getToolPrecision());
  
  polyCoeffs= (mpfr_t *)safeMalloc((n+1)*sizeof(mpfr_t));
  for (i=0; i<=n; i++){
  mpfr_init2(polyCoeffs[i], getToolPrecision());  
  }
  
  mpfi_init2(shiftx,getToolPrecision());
  mpfi_sub(shiftx,x,x0);
  
  mpfi_init2(zero,getToolPrecision());
  mpfi_set_ui(zero,0);
  
  mpfi_init2(err,getToolPrecision());
  
  mpfr_init2(a, getToolPrecision());
  mpfr_init2(b, getToolPrecision());    
  mpfi_get_left(a,shiftx);
  mpfi_get_right(b,shiftx);
  
  //transform the interval coeffs into mpfr + I;
  //printf("\nWe transform into mpfrs the coeffs of poly\n");
  
  mpfr_get_poly2(polyCoeffs, &r, n,coeffs, x0,x);
  
 //transform the polynomial with mpfr_coeffs into a node * p
 poly=makePolynomial(polyCoeffs, n);
 
 //derivate the polynomial
 diff_poly = differentiate(poly);
 //find the zeros
 zeros =uncertifiedFindZeros(diff_poly, a, b, points, getToolPrecision());
 nrRoots=lengthChain(zeros);
 rootsIntervals= (mpfi_t *)safeMalloc((nrRoots)*sizeof(mpfi_t));
 for(i=0;i<nrRoots;i++){
  mpfi_init2(rootsIntervals[i],getToolPrecision());
  mpfi_set_fr(rootsIntervals[i],*((mpfr_t*)zeros->value));
 }
 
 
  //compute the the bound for each small interval for the initial polynomial
  mpfr_init2(inf,getToolPrecision());
  mpfr_init2(sup,getToolPrecision());
  mpfr_init2(pinf,getToolPrecision());
  mpfr_init2(psup,getToolPrecision());
  //take the values in the extremas of the interval
  mpfi_init2(extr1,getToolPrecision());
  mpfi_set_fr(extr1,a);
  polynomialBoundHorner(&partialBound,n, coeffs,zero,extr1);
  mpfi_get_left(inf, partialBound);
  mpfi_get_right(sup,partialBound);
  
    for (i=0;i<nrRoots;i++){
      polynomialBoundHorner(&partialBound,n, coeffs,zero,rootsIntervals[i]);
      mpfi_get_left(pinf, partialBound);
      mpfi_get_right(psup,partialBound);
      if (mpfr_less_p(pinf,inf)!=0) mpfr_set(inf, pinf,GMP_RNDN);
      if (mpfr_greater_p(psup,sup)!=0) mpfr_set(sup,psup,GMP_RNDU);
    }
    
  mpfi_init2(extr1,getToolPrecision());
  mpfi_set_fr(extr1,b);
  polynomialBoundHorner(&partialBound,n, coeffs,zero,extr1);  
  mpfi_get_left(pinf, partialBound);
  mpfi_get_right(psup,partialBound);
  if (mpfr_less_p(pinf,inf)!=0) mpfr_set(inf, pinf,GMP_RNDN);
  if (mpfr_greater_p(psup,sup)!=0) mpfr_set(sup,psup,GMP_RNDU);
  
  mpfi_interv_fr(partialBound, inf,sup);
  mpfi_add(r,r,partialBound);
  mpfi_set(*bound,r);
  
  
  
  mpfi_clear(r);
  mpfi_clear(partialBound);
  mpfi_clear(extr1);
  //freeChain(zeros,freeMpfrPtr);
  for(i=0;i<nrRoots;i++){
  mpfi_clear(rootsIntervals[i]);
  }
  free(rootsIntervals);
  
  for(i=0;i<=n;i++){
    mpfr_clear(polyCoeffs[i]);
  }
  free(polyCoeffs);
  
  mpfr_clear(a);
  mpfr_clear(b);  
  mpfr_clear(inf);
  mpfr_clear(sup);
  mpfr_clear(pinf);
  mpfr_clear(psup);  
  free_memory(poly);
  free_memory(diff_poly);  
  
}
/*This function computes an interval bound
for a polynomial given by coeffs and order n, 
on int x, using whatever the user wants
*/
void polynomialBoundSharp(mpfi_t *bound,int n,mpfi_t *coeffs,mpfi_t x0,mpfi_t x){
//polynomialBoundHorner(bound,n,coeffs,x0,x);
polynomialBoundSharpUncertified(bound,n,coeffs,x0,x);
}



/*this function computes a cheby model for 1/x*/
void  varInv_CM(tModel *t,mpfi_t x0, mpfi_t x, int n){
    tModel *tt;
    mpfi_t *nDeriv;
    int i;
    node *f, *v, *c, *p;
    mpfr_t *coeffArray;
    mpfr_t fac;
    mpfr_t *ptr;
     mpfr_t minusOne;
    tt=createEmptycModel(n,x0,x);
  
  /*Use INTERPOLATION to compute the coeffs*/
  
  /*create a base function*/
  f = (node*)safeMalloc(sizeof(node));
  f->nodeType = DIV;
  
  v = (node*)safeMalloc(sizeof(node));
  v->nodeType = VARIABLE;
  f->child2=v;
  
  c = (node*)safeMalloc(sizeof(node));
  c->nodeType = CONSTANT;
  ptr = (mpfr_t*)safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*ptr, getToolPrecision());
  mpfr_set_ui(*ptr,1, GMP_RNDN);
  c->value = ptr;
  f->child1=c;
  
  printf("\n in cheby models for base functions: ");
  printTree(f);
  printf("\n");
  
  coeffArray=safeMalloc((n)*sizeof(mpfr_t));
  for (i=0;i<n;i++){
  mpfr_init2(coeffArray[i], getToolPrecision());
  }
  p=interpolation(coeffArray,&(tt->rem_bound), f,x,n-1);
  if (p==NULL) return;
  /*Use AD for the base functions to bound derivatives up to nth derivative*/
  nDeriv= (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    for(i=0;i<=n;i++){
      mpfi_init2(nDeriv[i], getToolPrecision());
    }
  mpfr_init2(minusOne, getToolPrecision());
  mpfr_set_si(minusOne, -1, GMP_RNDN);
  
  nDeriv= (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
  for(i=0;i<=n;i++){
    mpfi_init2(nDeriv[i], getToolPrecision());
  }
  constantPower_diff(nDeriv,minusOne, x, n);
 
  mpfi_mul(tt->rem_bound,tt->rem_bound, nDeriv[n]);
  mpfr_init2(fac, getToolPrecision());
  mpfr_fac_ui(fac,n,GMP_RNDN);
  mpfi_div_fr(tt->rem_bound, tt->rem_bound, fac);
  mpfr_clear(fac);
  
  for (i=0;i<n;i++){
  mpfi_set_fr(tt->poly_array[i], coeffArray[i]);
  }
  polynomialTranslateCM(tt->poly_array,n-1,tt->poly_array,tt->x0);
  
  /*bound the polynomial obtained*/
    
  polynomialBoundSharp(&tt->poly_bound, n-1,tt->poly_array,tt->x0,tt->x);   
  
  copycModel(t,tt);
  cleartModel(tt);
  mpfr_clear(minusOne);    
  for(i=0;i<=n;i++){
      mpfi_clear(nDeriv[i]);
  }
  free(nDeriv); 

}

/*this function computes a cheby model for x^p*/
void  ctPowerVar_CM(tModel *t,mpfi_t x0, mpfi_t x, int n, mpfr_t p){
    tModel *tt;
    mpfi_t *nDeriv;
    int i;
    node *f, *v, *c, *poly;
    mpfr_t *coeffArray;
    mpfr_t fac;
    mpfr_t *ptr;
    tt=createEmptycModel(n,x0,x);
  
  /*Use INTERPOLATION to compute the coeffs*/
  
  /*create a base function*/
  f = (node*)safeMalloc(sizeof(node));
  f->nodeType = POW;
  
  v = (node*)safeMalloc(sizeof(node));
  v->nodeType = VARIABLE;
  f->child1=v;
  
  c = (node*)safeMalloc(sizeof(node));
  c->nodeType = CONSTANT;
  ptr = (mpfr_t*)safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*ptr, getToolPrecision());
  mpfr_set(*ptr,p, GMP_RNDN);
  c->value = ptr;
  f->child2=c;
  
  printf("\n in cheby models for base functions: ");
  printTree(f);
  printf("\n");
  
  coeffArray=safeMalloc((n)*sizeof(mpfr_t));
  for (i=0;i<n;i++){
  mpfr_init2(coeffArray[i], getToolPrecision());
  }
  poly=interpolation(coeffArray,&(tt->rem_bound), f,x,n-1);
  if (poly==NULL) return;
  /*Use AD for the base functions to bound derivatives up to nth derivative*/
  nDeriv= (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    for(i=0;i<=n;i++){
      mpfi_init2(nDeriv[i], getToolPrecision());
    }
  nDeriv= (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
  for(i=0;i<=n;i++){
    mpfi_init2(nDeriv[i], getToolPrecision());
  }
  constantPower_diff(nDeriv,p, x, n);
 
  mpfi_mul(tt->rem_bound,tt->rem_bound, nDeriv[n]);
  mpfr_init2(fac, getToolPrecision());
  mpfr_fac_ui(fac,n,GMP_RNDN);
  mpfi_div_fr(tt->rem_bound, tt->rem_bound, fac);
  mpfr_clear(fac);
  
  for (i=0;i<n;i++){
  mpfi_set_fr(tt->poly_array[i], coeffArray[i]);
  }
  polynomialTranslateCM(tt->poly_array,n-1,tt->poly_array,tt->x0);
  
  /*bound the polynomial obtained*/
    
  polynomialBoundSharp(&tt->poly_bound, n-1,tt->poly_array,tt->x0,tt->x);   
  
  copycModel(t,tt);
  cleartModel(tt);
    
  for(i=0;i<=n;i++){
      mpfi_clear(nDeriv[i]);
  }
  free(nDeriv); 
    

}

/*this function computes the cheby model for p^x*/
void  varCtPower_CM(tModel *t,mpfi_t x0, mpfi_t x, int n, mpfr_t p){
    tModel *tt;
    mpfi_t *nDeriv;
    int i;
    node *f, *v, *c, *poly;
    mpfr_t *coeffArray;
    mpfr_t fac;
    mpfr_t *ptr;
    tt=createEmptycModel(n,x0,x);
  
  /*Use INTERPOLATION to compute the coeffs*/
  
  /*create a base function*/
  f = (node*)safeMalloc(sizeof(node));
  f->nodeType = POW;
  
  v = (node*)safeMalloc(sizeof(node));
  v->nodeType = VARIABLE;
  f->child2=v;
  
  c = (node*)safeMalloc(sizeof(node));
  c->nodeType = CONSTANT;
  ptr = (mpfr_t*)safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*ptr, getToolPrecision());
  mpfr_set(*ptr,p, GMP_RNDN);
  c->value = ptr;
  f->child1=c;
  
  printf("\n in cheby models for base functions: ");
  printTree(f);
  printf("\n");
  
  coeffArray=safeMalloc((n)*sizeof(mpfr_t));
  for (i=0;i<n;i++){
  mpfr_init2(coeffArray[i], getToolPrecision());
  }
  poly=interpolation(coeffArray,&(tt->rem_bound), f,x,n-1);
  if (poly==NULL) return;
  /*Use AD for the base functions to bound derivatives up to nth derivative*/
  nDeriv= (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
    for(i=0;i<=n;i++){
      mpfi_init2(nDeriv[i], getToolPrecision());
    }
  nDeriv= (mpfi_t *)safeMalloc((n+1)*sizeof(mpfi_t));
  for(i=0;i<=n;i++){
    mpfi_init2(nDeriv[i], getToolPrecision());
  }
  powerFunction_diff(nDeriv,p, x, n);
 
 
  mpfi_mul(tt->rem_bound,tt->rem_bound, nDeriv[n]);
  mpfr_init2(fac, getToolPrecision());
  mpfr_fac_ui(fac,n,GMP_RNDN);
  mpfi_div_fr(tt->rem_bound, tt->rem_bound, fac);
  mpfr_clear(fac);
  
  for (i=0;i<n;i++){
  mpfi_set_fr(tt->poly_array[i], coeffArray[i]);
  }
  polynomialTranslateCM(tt->poly_array,n-1,tt->poly_array,tt->x0);
  
  /*bound the polynomial obtained*/
    
  polynomialBoundSharp(&tt->poly_bound, n-1,tt->poly_array,tt->x0,tt->x);   
  
  copycModel(t,tt);
  cleartModel(tt);
    
  for(i=0;i<=n;i++){
      mpfi_clear(nDeriv[i]);
  }
  free(nDeriv); 
}    
   


/**Working lineL ana are mere**/

void cheby_model(tModel *t, node *f, int n, mpfi_t x0, mpfi_t x) {
  int i;
  
  node *simplifiedChild1, *simplifiedChild2;
  
  mpfi_t temp1,temp2;
  node **coefficients;
  //mpfi_t *rpoly, *boundRpoly;
  tModel *tt, *child1_tm, *child2_tm,*ctPowVar_tm, *varCtPower_tm, *logx_tm, *expx_tm, *logf_tm;;//*tPoly, 
  node *p;
  /*used by division*/
  mpfi_t gx0,rangeg;
  tModel *ttt, *inv_tm;//, *child1Extended_tm, *child2Extended_tm, *child1RemoveCoeffs_tm,*child2RemoveCoeffs_tm; 
  //int orderUpperBound;  
  int d;
  //mpfi_t powx,powy;
  mpfi_t  minusOne; //ct,
  //mpfr_t zero;
  /*used by base functions*/
  mpfi_t rangef,fx0;//pow;
  /*used by tm2 in base functions*/
  mpfi_t fcoeff0;//,c;
  //mpfr_t fmid;
  //tModel *child2Mid_tm;
  
  mpfr_t *coeffArray;
  coeffArray=safeMalloc((n)*sizeof(mpfr_t));
  for (i=0;i<n;i++){
  mpfr_init2(coeffArray[i], getToolPrecision());
  }
  if (isPolynomial(f) ){
    //printf("We found a polynomial!!!\n");
    tt=createEmptycModel(n,x0,x); 
    getCoefficients(&d, &coefficients, f);
    //for(i=0;i<=d;i++) if (coefficients[i]!=NULL) printTree(coefficients[i]);
    if (d<n){
      //printf("The degree of polynomial smaller : %d than the requested degree\n",d);
      for(i=0;i<=d;i++) {
        //printf("%d\n",i);
        //
        mpfi_set_node(&(tt->poly_array[i]),coefficients[i]);
        //printInterval(tt->poly_array[i]);
      }
      
      for(i=d+1;i<n;i++) {
        //printf("%d\n",i);
        //printInterval(t->poly_array[i]);
        mpfi_set_ui(tt->poly_array[i],0);
      }
      //printcModel(tt);
      //printf("-------------Before poly translate-----------");
      polynomialTranslateCM(tt->poly_array,n-1,tt->poly_array,tt->x0);
      
      //printf("-------------After poly translate-----------");
      
      
      mpfi_set_ui(tt->rem_bound,0);
       polynomialBoundSharp(&tt->poly_bound, n-1,tt->poly_array,tt->x0,tt->x);   
      // printf("we set the coefficients\n");
      printf("\n\n*****************The polynomial model is:\n********* ");
      printcModel(tt);
    }
    else {
      printf("The degree of polynomial bigger : %d then the cheby model degree\n",d);
      
      polynomialNode_CM(tt, f, n, x0, x);
      printf("\nThe polynomial obtained is:");
      polynomialTranslateCM(tt->poly_array,n-1,tt->poly_array,x0);
    }
    copycModel(t,tt);   
    cleartModel(tt);
  }
  else {  
    switch (f->nodeType) {
    //case VARIABLE:
    //case PI_CONST:
    //case CONSTANT:
    case NEG:
  
      tt=createEmptycModel(n,x0,x);
      //create a new empty taylor model the child
      child1_tm=createEmptycModel(n, x0, x);
      //call cheby_model on the child
      cheby_model(child1_tm, f->child1,n,x0,x);
      //do the necessary chages from child to parent
      for(i=0;i<n;i++) 
      mpfi_neg(tt->poly_array[i], child1_tm->poly_array[i]);
    
      mpfi_neg(tt->rem_bound,child1_tm->rem_bound);
      mpfi_neg(tt->poly_bound,child1_tm->poly_bound);
      copycModel(t,tt);
      //clear old cheby models
      cleartModel(child1_tm);
      cleartModel(tt);
    
      break;

    case ADD:
    
      //create a new empty cheby model the children
      tt=createEmptycModel(n,x0,x); 
      child1_tm=createEmptycModel(n, x0,x);
      child2_tm=createEmptycModel(n, x0,x);
      //call cheby_model on the children
      cheby_model(child1_tm, f->child1,n,x0,x);
      cheby_model(child2_tm, f->child2,n,x0,x);
    
      addition_CM(tt,child1_tm, child2_tm);
      copycModel(t,tt);
      //clear old cheby model
      cleartModel(child1_tm);
      cleartModel(child2_tm);
      cleartModel(tt);
      break;

    case SUB:
    
      tt=createEmptycModel(n,x0,x); 
      //create a new empty cheby model the children
      child1_tm=createEmptycModel(n, x0,x);
      child2_tm=createEmptycModel(n, x0,x);
      //call cheby_model on the children
      cheby_model(child1_tm, f->child1,n,x0,x);
      cheby_model(child2_tm, f->child2,n,x0,x);
      
      //do the necessary chages from children to parent
      mpfi_init2(minusOne,getToolPrecision());
      mpfi_set_si(minusOne,-1);
      ctMultiplication_CM(child2_tm,child2_tm, minusOne);
      addition_CM(tt,child1_tm, child2_tm);
      
      copycModel(t,tt);
    
      //clear old cheby model
      cleartModel(child1_tm);
      cleartModel(child2_tm);
      cleartModel(tt);
      mpfi_clear(minusOne);
      break;

    case MUL:
     tt=createEmptycModel(n,x0,x); 
   //create a new empty taylor model the children
    child1_tm=createEmptycModel(n,x0,x);
    child2_tm=createEmptycModel(n, x0,x);
    //call cheby_model on the children
    cheby_model(child1_tm, f->child1,n,x0,x);
    cheby_model(child2_tm, f->child2,n,x0,x);
    printf("\n\n\n\n\n\n\n****************************************************\n");
    printf("****************************************************\n");
    printf("****************************************************\n");
    printf("****************************************************\n");    
    printcModel(child1_tm);
    printcModel(child2_tm);
    printf("****************************************************\n");
    printf("****************************************************\n");
    printf("****************************************************\n");
    printf("****************************************************\n\n\n\n\n");    
       
    //do the necessary chages from children to parent
    multiplication_CM(tt,child1_tm, child2_tm);
    copycModel(t,tt);
    //clear old cheby model
    cleartModel(child1_tm);
    cleartModel(child2_tm);     
    cleartModel(tt);
    break;

  case DIV:
    //child1 * inverse(child2)
  tt=createEmptycModel(n,x0,x); 
    
  
  ttt=createEmptycModel(n,x0,x); 
   
  //just create tms or order n directly
  //create a new empty cheby model the child
  child1_tm=createEmptycModel(n,x0,x);
  child2_tm=createEmptycModel(n,x0,x);
  //call cheby_model for the children
  cheby_model(child1_tm, f->child1,n,x0,x);
  cheby_model(child2_tm, f->child2,n,x0,x);
  mpfi_init2(gx0,getToolPrecision());
 // mpfi_set(gx0,child2_tm->poly_array[0]);
  mpfi_set_ui(gx0,0);
  mpfi_init2(rangeg, getToolPrecision());
  mpfi_add(rangeg,child2_tm->rem_bound,child2_tm->poly_bound);
  
  inv_tm=createEmptycModel(n,gx0,rangeg);
  varInv_CM(inv_tm,gx0,rangeg,n);
  composition_CM(ttt,inv_tm,child2_tm);
  multiplication_CM(tt, ttt, child1_tm);
  //clear old children
  cleartModel(child1_tm);
  cleartModel(child2_tm);
  cleartModel(inv_tm);
  cleartModel(ttt);
  copycModel(t,tt);
  cleartModel(tt);
  mpfi_clear(rangeg);
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
  
    tt=createEmptycModel(n,x0,x); 
   //create a new empty cheby model the child
    child1_tm=createEmptycModel(n,x0,x);
    child2_tm=createEmptycModel(n,x0,x);
    //call taylor_model on the child
    cheby_model(child1_tm, f->child1,n,x0,x);
    //compute tm for the basic case
    printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
    printf("\nThe cheb model for: ");
    printTree(f->child1);
    printf(" is\n");
    printcModel(child1_tm);
    printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");    
    
    mpfi_init2(rangef, getToolPrecision());
    mpfi_init2(fcoeff0,getToolPrecision());
    //mpfi_set(fcoeff0,child1_tm->poly_array[0]);
    mpfi_set_ui(fcoeff0,0);
    mpfi_add(rangef,child1_tm->rem_bound, child1_tm->poly_bound);
    printInterval(rangef);
    child2_tm=createEmptycModel(n,fcoeff0, rangef);
    printf("\n************* cheby model for elementary function************\n");
    base_CM(child2_tm,f->nodeType,n,fcoeff0, rangef);
    printf("\n************* result************\n");   
    
    printcModel(child2_tm);
    
    printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
    printf("\nThe cheb model for: ");
    printTree(f);
    printf(" is\n");
    printcModel(child2_tm);
    printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");    
    printf("We compose");
    composition_CM(tt,child2_tm, child1_tm);
    printf("\n\n*************after the composition************\n");   
    printcModel(tt);
    
    //clear old child
    cleartModel(child1_tm);
    cleartModel(child2_tm);
    copycModel(t,tt);
    cleartModel(tt);
    mpfi_clear(rangef);    
    mpfi_clear(fcoeff0);
    
        
  break;
  case POW:
  
    //create the result chebyshev model
    tt=createEmptycModel(n,x0,x); 
  
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
         copycModel(t,tt);
         cleartModel(tt);
         mpfi_clear(temp1);
         mpfi_clear(temp2);
      }
      else if (simplifiedChild2->nodeType==CONSTANT) { //we have the f^p case
        printf("We are in the  f^p case");        
        //mpfi_t powx,fx0,rangef;        
        //create a new empty cheby model for the child
        child1_tm=createEmptycModel(n,x0,x);
        //call taylor_model for the child
        cheby_model(child1_tm, f->child1,n,x0,x);
        //printf("\n\n-----------chebymodel child1: \n");
        //printcModel(child1_tm);
        //printf("-----------------------------\n");
        mpfi_init2(fx0,getToolPrecision());
        mpfi_set(fx0, child1_tm->poly_array[0]);
        //evaluateMpfiFunction(fx0,f->child1,x0,getToolPrecision());
        
        mpfi_init2(rangef, getToolPrecision());
        mpfi_add(rangef,child1_tm->rem_bound, child1_tm->poly_bound);
        
        //create tm for x^p over ragef in fx0
        ctPowVar_tm=createEmptycModel(n,fx0,rangef);
        
        ctPowerVar_CM(ctPowVar_tm,fx0,rangef,n,*(simplifiedChild2->value));
        
        //printf("\n\n-----------chebymodel child1: \n");
        //printcModel(ctPowVar_tm);
        //printf("-----------------------------\n");
        
        
        composition_CM(tt,ctPowVar_tm,child1_tm);
    
        //clear old child
        cleartModel(child1_tm);
        cleartModel(ctPowVar_tm);
        copycModel(t,tt);
        cleartModel(tt);
        mpfi_clear(rangef);
        mpfi_clear(fx0);
      } 
       else if (simplifiedChild1->nodeType==CONSTANT) { //we have the p^f case
         printf("We are in the  p^f case");  
        //create a new empty taylor model the child
        child2_tm=createEmptycModel(n,x0,x);
        //call cheby_model for the child
        cheby_model(child2_tm, f->child2,n,x0,x);
        
        mpfi_init2(fx0,getToolPrecision());
        mpfi_set(fx0, child2_tm->poly_array[0]);
        
        mpfi_init2(rangef, getToolPrecision());
        
        mpfi_add(rangef,child2_tm->rem_bound, child2_tm->poly_bound);
        
        //create tm for p^x over ragef in fx0
        varCtPower_tm=createEmptycModel(n,fx0,rangef);
        
        varCtPower_CM(varCtPower_tm,fx0,rangef,n,*(simplifiedChild1->value));
        printf("ajunge aici");
        composition_CM(tt,varCtPower_tm,child2_tm);
    
        //clear old child
        cleartModel(child2_tm);
        cleartModel(varCtPower_tm);
        copycModel(t,tt);
        cleartModel(tt);
        mpfi_clear(rangef);
        mpfi_clear(fx0);
        }
      else {
      //printf("We are in the  f^g case");  
      //exp(g log(f))
          
        //create a new empty cheby model the children
        child1_tm=createEmptycModel(n,x0,x);
        child2_tm=createEmptycModel(n,x0,x);
        //call cheby_model for child 2 = g
        cheby_model(child2_tm, f->child2,n,x0,x);
        
        //call cheby_model for child 1 = f
        cheby_model(child1_tm, f->child1,n,x0,x);
        
        //create  cheby_model for log (child 1) = log(f)
                
        mpfi_init2(fx0,getToolPrecision());
        mpfi_set(fx0, child1_tm->poly_array[0]);
        
        mpfi_init2(rangef, getToolPrecision());
        mpfi_add(rangef,child1_tm->rem_bound, child1_tm->poly_bound);
        
        
        logx_tm=createEmptycModel(n,fx0,rangef);
        base_CM(logx_tm,LOG,n,fx0, rangef);
        logf_tm=createEmptycModel(n,x0,x);
        composition_CM(logf_tm,logx_tm,child1_tm);
        //-------------------------------------------
        
        
        //multiply g by log f
        ttt=createEmptycModel(n,x0,x);
        multiplication_CM(ttt,logf_tm,child2_tm);
        
        //----exp()--------------------------------------
        mpfi_init2(gx0,getToolPrecision());
        mpfi_set(gx0,ttt->poly_array[0]);
        mpfi_add(rangef,ttt->rem_bound, ttt->poly_bound);
            
        expx_tm=createEmptycModel(n,gx0,rangef);
        base_CM(expx_tm,EXP,n,gx0, rangef);
        composition_CM(tt,expx_tm,ttt);
               
    
        //clear old child
        cleartModel(child2_tm);
        cleartModel(child1_tm);
        cleartModel(ttt);
        cleartModel(expx_tm);
        cleartModel(logx_tm);
        cleartModel(logf_tm);
        
        copycModel(t,tt);
        cleartModel(tt);
        mpfi_clear(rangef);
        
        mpfi_clear(fx0);   
        mpfi_clear(gx0); 
        
    }
    free_memory(simplifiedChild2);
    free_memory(simplifiedChild1);
  
  
  
  
  break;
  case LIBRARYFUNCTION:
  break;

  default:
   fprintf(stderr,"Error: CM: unknown identifier (%d) in the tree\n",f->nodeType);
   exit(1);

}
}
  return;
}
/*This function computes a translation for the polynomial P(x) to the polynomial P(x-x0)
*/

void polynomialTranslateCM(mpfi_t *coeffsT,int n,mpfi_t *coeffs,mpfi_t x0){
  mpfi_t temp;
  int i,k;
  mpfi_t *coeffArray;
  coeffArray=safeMalloc((n+1)*sizeof(mpfi_t));
  
  
  
  for (i=0;i<=n;i++){
    mpfi_init2(coeffArray[i], getToolPrecision());
    mpfi_set(coeffArray[i],coeffs[i]);
  }
  
  mpfi_init2(temp, getToolPrecision());
  
  for (i=1; i<=n; i++)
    for(k=n-i;k<=n-1;k++){
      mpfi_mul(temp, x0, coeffArray[k+1]);
      //for translation to x-x0, we have to multipliy by -1
      //mpfi_mul_si(temp,temp,-1);
      mpfi_add(coeffArray[k],coeffArray[k],temp);  
    }
    mpfi_clear(temp);
   
  for (i=0; i<=n;i++){
     printInterval(coeffArray[i]);
    mpfi_set(coeffsT[i], coeffArray[i]);
  }
  
  for (i=0; i<=n;i++){
    mpfi_clear(coeffArray[i]);
  }
  free(coeffArray);
  
  
  //printf("ana are mre");
}
node *constructPolyShifted(mpfr_t *coeff, int n, mpfi_t x0) {
  int i=1;
  //chain *curr;
  node *temp1;
  node *temp2;
  node *temp3;
  node *temp4;
  node *temp5;
  node *temp6;

  node *poly;
  mpfr_t *ptr;

  poly = (node*)safeMalloc(sizeof(node));
  poly->nodeType = CONSTANT;
  ptr = (mpfr_t*)safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*ptr, getToolPrecision());
  mpfr_set(*ptr,coeff[0], GMP_RNDN);
  poly->value = ptr;

  
  for (i=1;i<=n;i++) {
    temp1 = (node*)safeMalloc(sizeof(node));
    temp1->nodeType = ADD;

    temp2 =(node*) safeMalloc(sizeof(node));
    temp2->nodeType = MUL;

    temp3 = (node*)safeMalloc(sizeof(node));
    temp3->nodeType = CONSTANT;
    ptr = (mpfr_t*)safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*ptr, getToolPrecision());
    mpfr_set(*ptr, coeff[i], GMP_RNDN);
    temp3->value = ptr;

    temp4 = (node*)safeMalloc(sizeof(node));
    temp4->nodeType = POW;

    temp5 = (node*)safeMalloc(sizeof(node));
    temp5->nodeType = VARIABLE;

    temp6 = (node*)safeMalloc(sizeof(node));
    temp6->nodeType = CONSTANT;
    ptr = (mpfr_t*)safeMalloc(sizeof(mpfr_t));
    mpfr_init2(*ptr, getToolPrecision());
    mpfr_set_si(*ptr, i, GMP_RNDN);
    temp6->value = ptr;

    temp4->child1 = temp5;
    temp4->child2 = temp6;
    temp2->child1 = temp3;
    temp2->child2 = temp4;
    temp1->child1 = temp2;
    temp1->child2 = poly;
    poly = temp1;
   }

  return poly;
}
chain *constructChain(mpfi_t *err, int n){
  chain *l;
  mpfi_t *elem;
  int i;
  
  l = NULL;
  for(i=n;i>=0;i--){
    elem= (mpfi_t*)safeMalloc(sizeof(mpfi_t));
    mpfi_init2(*elem,getToolPrecision());
    mpfi_set(*elem,err[i]);
    l=addElement(l, elem);
  }
  return l;
}
void printMpfiChain(chain *c) {
  chain *curr=c;
  printf("[");
  while(curr!=NULL) {
    printInterval( *(mpfi_t *)(curr->value));
    curr=curr->next;
  }
  printf("]\n");
  return;
}


int CM(chain**resP, void **args) {
  node *f;
  mpfi_t x,x0 ;// , bound;
  mpfr_t *coeffsMpfr;
  mpfi_t *coeffsErrors, *rest;
  int n,i;
  chain *ch;
  node *resultPoly;
  node *remainderL, *remainderR;
  mpfr_t *ptr;
  ch=NULL;
  
  f = (node *)args[0];
  n = *( (int *)args[3] );

  mpfi_init2(x0, getToolPrecision());
  mpfi_set(x0, *( (mpfi_t *)args[1] ));
  
  mpfi_init2(x, getToolPrecision());
  mpfi_set(x, *( (mpfi_t *)args[2] ));
  
  tModel *t;
  
  t=createEmptycModel(n,x0,x);
  cheby_model(t,f,n,x0,x);
  
  
  printcModel(t);

  coeffsMpfr= (mpfr_t *)safeMalloc((n)*sizeof(mpfr_t));
  coeffsErrors= (mpfi_t *)safeMalloc((n)*sizeof(mpfi_t));

  rest= (mpfi_t*)safeMalloc(sizeof(mpfi_t));
  mpfi_init2(*rest,getToolPrecision());

  for(i=0;i<n;i++){
   mpfi_init2(coeffsErrors[i],getToolPrecision());
   mpfr_init2(coeffsMpfr[i],getToolPrecision());
  }
  //mpfr_get_poly(mpfr_t *rc, mpfi_t *errors, mpfi_t rest, int n, mpfi_t *gc, mpfi_t x0, mpfi_t x)
  mpfr_get_poly(coeffsMpfr, coeffsErrors, *rest, t->n -1,t->poly_array, t->x0,t->x);
 
 //create T; 
 resultPoly=constructPolyShifted(coeffsMpfr, t->n-1, t->x0);
  mpfi_set(*rest,t->rem_bound);  
  printInterval(*rest);
  
  
  
  remainderL = (node*)safeMalloc(sizeof(node));
  remainderL->nodeType = CONSTANT;
  ptr = (mpfr_t*)safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*ptr, getToolPrecision());
  mpfi_get_left(*ptr,t->rem_bound);
  remainderL->value = ptr;
  
  remainderR = (node*)safeMalloc(sizeof(node));
  remainderR->nodeType = CONSTANT;
  ptr = (mpfr_t*)safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*ptr, getToolPrecision());
  mpfi_get_right(*ptr,t->rem_bound);
  remainderR->value = ptr;
  
  
  
  
   ch=addElement(ch, remainderR);
   ch=addElement(ch, remainderL);
   ch=addElement(ch, resultPoly);
  *resP=ch;
  
  for(i=0;i<n;i++){
    mpfr_clear(coeffsMpfr[i]);
  }
  free(coeffsMpfr);
  
  cleartModel(t);

   
  mpfi_clear(x);
	mpfi_clear(x0);   
  return 1;
}


