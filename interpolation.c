/*
  For compiling this file:
    gcc -fPIC -Wall -c interpolation.c
    gcc -fPIC -shared -o interpolation interpolation.o


  Within Sollya:
    > externalproc(IP, "./interpolation", (function, range, integer) -> list of function);

  And then, for instance:
    > IP(exp(x)*cos(x), [2.5; 2.6], 10);

*/

#include "interpolation.h"

int IP(chain**resP, void **args) {
  node *f;
  mpfi_t x, bound;
  mpfr_t *coeffArray;
  node *resultPoly;
  node *remainderL, *remainderR;
  mpfr_t *ptr;
  chain *ch;
  ch=NULL;
  int n,i;
  
  f = (node *)args[0];
  n = *( (int *)args[2] );

  mpfi_init2(x, getToolPrecision());
  mpfi_set(x, *( (mpfi_t *)args[1] ));
  mpfi_init2(bound, getToolPrecision());
  
  coeffArray=safeMalloc((n+1)*sizeof(mpfr_t));
  for (i=0;i<=n;i++){
  mpfr_init2(coeffArray[i], getToolPrecision());
  
  }
  
  
  resultPoly = interpolation(coeffArray, bound, f, x, n);
  //printTree(*resP);
  remainderL = (node*)safeMalloc(sizeof(node));
  remainderL->nodeType = CONSTANT;
  ptr = (mpfr_t*)safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*ptr, getToolPrecision());
  mpfi_get_left(*ptr,bound);
  remainderL->value = ptr;
  
  remainderR = (node*)safeMalloc(sizeof(node));
  remainderR->nodeType = CONSTANT;
  ptr = (mpfr_t*)safeMalloc(sizeof(mpfr_t));
  mpfr_init2(*ptr, getToolPrecision());
  mpfi_get_right(*ptr,bound);
  remainderR->value = ptr;
  
  
  
  
   ch=addElement(ch, remainderR);
   ch=addElement(ch, remainderL);
   ch=addElement(ch, resultPoly);
  *resP=ch;
   
  mpfi_clear(x);
	   
  return 1;
}
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




 
   
node* interpolation (mpfr_t *coeffArray, mpfi_t wbound, node *f, mpfi_t x, int n) {
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
  int points=100;
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
  
  if (nrRoots>=(freeDegrees+1)) {
    printf("\n-------we have the more or enough cheby points----------\n");
    if (nrRoots>(freeDegrees+1)){
      findNewSetRoots(newRoots, roots, nrRoots, chebArray, freeDegrees+1);
      w=constructPolyFromRoots(newRoots,freeDegrees+1);
    
      for(i=0;i<freeDegrees+1;i++){
        mpfi_init2(rootsIntervals[i],getToolPrecision());
        // printf("we are going for the interval root\n");
        if (0==(getIntervalAroundRoot(&rootsIntervals[i],dif, newRoots[i], getToolPrecision()))){printf("Huge Warning one root can not be made as an     interval!!!!!!!!!!");}
       //printInterval(rootsIntervals[i]);
      }
   
    } 
  else{
  printf("\n-------we have the good nr of cheby points----------\n");
  w=constructPolyFromRoots(roots,freeDegrees+1);
  for(i=0;i<freeDegrees+1;i++){
    mpfi_init2(rootsIntervals[i],getToolPrecision());
    if (0==(getIntervalAroundRoot(&rootsIntervals[i],dif, roots[i], getToolPrecision()))){printf("Huge Warning one root can not be made as an interval!!!!!!!!!!");}
    //printInterval(rootsIntervals[i]);
    }
  
  }
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
  mpfi_interv_fr(wbound, boundNeg, boundPos);
  
  mpfr_clear(boundNeg);
  mpfr_clear(boundPos);
  }
  else{
  printf("We have the wrong nr of cheby points");
  }
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

