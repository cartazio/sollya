#include <stdio.h>
#include "expression.h"
#include "infnorm.h"

int constantImplementer(node *c, int gamma0, char *resName, int counter);

void implementCste(node *c) {
  int i, counter;

  printf("void\n");
  printf("mpfr_const_xxx (mpfr_ptr y, mp_prec_t prec)\n");
  printf("{\n");
  printf("  /* Declarations */\n");
  printf("  /* Initializations */\n");
  printf("\n");
  
  counter = 0;
  counter = constantImplementer(c, 0, "y", counter);

  printf("\n");
  printf("  /* Cleaning stuff */\n");
  for(i=1; i<=counter; i++) printf("  mpfr_t tmp%d;\n", i);
  for(i=1; i<=counter; i++) printf("  mpfr_init (tmp%d);\n", i);
  for(i=1; i<=counter; i++) printf("  mpfr_clear (tmp%d);\n", i);
  printf("}\n");

  return;
}

mp_exp_t mpfi_get_exp(mpfi_t x) {
  mpfr_t u,v;
  mp_exp_t Eu, Ev, E;
  mp_prec_t prec;

  prec = mpfi_get_prec(x);
  mpfr_init2(u, prec);
  mpfr_init2(v, prec);

  mpfi_get_left(u, x);
  mpfi_get_right(v, x);
  
  if (mpfr_zero_p(u)) E = mpfr_get_exp(v);
  else {
    if (mpfr_zero_p(v)) E = mpfr_get_exp(u);
    else {
      Eu = mpfr_get_exp(u);
      Ev = mpfr_get_exp(v);
      E = (Eu<=Ev) ? Ev : Eu;
    }
  }
  mpfr_clear(u);
  mpfr_clear(v);
  return E;
}
 
/* Let a be the constant given by the expression cste and f the function with */
/* node type nodeType. This functions generates code for the implementation   */
/* of f(a) in precision prec+gamma0, the result being stored in resName.      */
int unaryFunctionCase(int nodeType, node *cste, char *functionName, int gamma0, char *resName, char * tmpName, int counter) {
  mpfi_t a, b, u, v, tmp;
  mpfr_t alpha, beta;
  mp_prec_t prec = getToolPrecision();
  node *func, *deriv;
  int gamma;
  int toReturn;

  mpfi_init2(a, prec);
  mpfi_init2(b, prec);
  mpfi_init2(u, prec);
  mpfi_init2(v, prec);
  mpfi_init2(tmp, prec);
  mpfr_init2(alpha, prec);
  mpfr_init2(beta, prec);

  func = makeUnary(makeVariable(), nodeType);
  deriv = differentiate(func);
  
  evaluateInterval(a, cste, NULL, a);
  evaluateInterval(b, func, deriv, a);
  mpfi_div(u, a, b);
  evaluateInterval(tmp, deriv, NULL, a); 
  mpfi_mul(v, u, tmp);
  do {
    gamma = 2+mpfi_get_exp(v);

    mpfr_set_ui(beta, 1, GMP_RNDU);
    mpfr_div_2si(beta, beta, gamma+gamma0, GMP_RNDU);
    mpfr_ui_sub(alpha, 1, beta, GMP_RNDD);
    mpfr_add_ui(beta, beta, 1, GMP_RNDU);
    mpfi_interv_fr(tmp, alpha, beta);

    mpfi_mul(tmp, a, tmp);
    evaluateInterval(tmp, deriv, NULL, tmp);
    mpfi_mul(v, u, tmp);
  } while (gamma < 2+mpfi_get_exp(v));
  
  toReturn = constantImplementer(cste, gamma0+gamma, tmpName, counter);
  printf("  mpfr_set_prec (%s, prec+%d);\n", resName, gamma0+2);
  printf("  mpfr_%s (%s, %s, MPFR_RNDN);\n", functionName, resName, tmpName);

  mpfi_clear(a);
  mpfi_clear(b);
  mpfi_clear(u);
  mpfi_clear(v);
  mpfi_clear(tmp);
  mpfr_clear(alpha);
  mpfr_clear(beta);

  return toReturn;
}

void normalizeDivMul(node *c, chain **numerator, chain **denominator) {
  chain *num1 = NULL;
  chain *denom1 = NULL;
  chain *num2 = NULL;
  chain *denom2 = NULL;

  if (c->nodeType == MUL) {
    normalizeDivMul(c->child1, &num1, &denom1);
    normalizeDivMul(c->child2, &num2, &denom2);
    *numerator = concatChains(num1, num2);
    *denominator = concatChains(denom1, denom2);
  }
  else if (c->nodeType == DIV) {
    normalizeDivMul(c->child1, &num1, &denom1);
    normalizeDivMul(c->child2, &denom2, &num2);
    *numerator = concatChains(num1, num2);
    *denominator = concatChains(denom1, denom2);
  }
  else *numerator = addElement(*numerator, copyTree(c));
}


int implementDivMul(node *c, int gamma0, char *resName, int counter) {
  chain *numerator = NULL;
  chain *denominator = NULL;
  chain *curr;
  chain *bufferNum, *bufferDenom;
  int log2n, n, test;
  char tmpName[10] = "tmp";
  int toReturn = counter;
  int *tmp;

  normalizeDivMul(c, &numerator, &denominator);

  n = lengthChain(numerator) + lengthChain(denominator);
  log2n = 0;
  test = 1;
  /* Compute log2n such that 2^(log2n-1) <= n < 2^(log2n) */
  while (n>=1) { log2n++; if ((n%2)!=0) test=0; n = n/2;} 
  /* Ajust log2n in order to have 2^(log2n-1) < n <= 2^(log2n) */
  if(test) log2n--;

  curr = numerator;
  bufferNum = NULL;
  while(curr!=NULL) {
    sprintf(tmpName+3, "%d", toReturn+1);
    tmp = (int *)safeMalloc(sizeof(int));
    *tmp = toReturn+1;
    bufferNum = addElement(bufferNum, tmp);
    toReturn = constantImplementer(curr->value, gamma0+2+log2n, tmpName, toReturn+1);
    curr = curr->next;
  }
  curr = denominator;
  bufferDenom = NULL;
  while(curr!=NULL) {
    sprintf(tmpName+3, "%d", toReturn+1);
    tmp = safeMalloc(sizeof(int));
    *tmp = toReturn+1;
    bufferDenom = addElement(bufferDenom, tmp);
    toReturn = constantImplementer(curr->value, gamma0+2+log2n, tmpName, toReturn+1);
    curr = curr->next;
  }

  if ( (lengthChain(numerator)==1) && (lengthChain(denominator)==1) ) {
    printf("  mpfr_set_prec (%s, prec+%d);\n", resName, gamma0+2+log2n);
    printf("  mpfr_div (%s, tmp%d, tmp%d);\n", resName, *((int *)(bufferNum->value)), *((int *)(bufferDenom->value)));
  }
  else if (lengthChain(numerator)==1) {
    printf("  mpfr_set_prec (%s, prec+%d);\n", resName, gamma0+2+log2n);
    printf("  mpfr_mul (%s, tmp%d, tmp%d, MPFR_RNDN);\n", resName, *((int *)(bufferDenom->value)), *((int *)(bufferDenom->next->value)));
    curr = bufferDenom->next->next;
    while(curr!=NULL) {
      printf("  mpfr_mul (%s, %s, tmp%d, MPFR_RNDN);\n", resName, resName, *((int *)(curr->value)));
      curr = curr->next;
    }
    printf("  mpfr_div (%s, tmp%d, %s, MPFR_RNDN);\n", resName, *((int *)(bufferNum->value)), resName);
  }
  else if (lengthChain(denominator)<=1) {
    printf("  mpfr_set_prec (%s, prec+%d);\n", resName, gamma0+2+log2n);
    printf("  mpfr_mul (%s, tmp%d, tmp%d, MPFR_RNDN);\n", resName, *((int *)(bufferNum->value)), *((int *)(bufferNum->next->value)));
    curr = bufferNum->next->next;
    while(curr!=NULL) {
      printf("  mpfr_mul (%s, %s, tmp%d, MPFR_RNDN);\n", resName, resName, *((int *)(curr->value)));
      curr = curr->next;
    }
    if (lengthChain(denominator)==1)
      printf("  mpfr_div (%s, %s, tmp%d, MPFR_RNDN);\n", resName, resName, *((int *)(bufferDenom->value)));
  }
  else {
    printf("  mpfr_set_prec (%s, prec+%d);\n", resName, gamma0+2+log2n);
    printf("  mpfr_mul (%s, tmp%d, tmp%d, MPFR_RNDN);\n", resName, *((int *)(bufferNum->value)), *((int *)(bufferNum->next->value)));
    curr = bufferNum->next->next;
    while(curr!=NULL) {
      printf("  mpfr_mul (%s, %s, tmp%d, MPFR_RNDN);\n", resName, resName, *((int *)(curr->value)));
      curr = curr->next;
    }
    toReturn++;
    printf("  mpfr_set_prec (tmp%d, prec+%d);\n", toReturn, gamma0+2+log2n);
    printf("  mpfr_mul (tmp%d, tmp%d, tmp%d, MPFR_RNDN);\n", toReturn, *((int *)(bufferDenom->value)), *((int *)(bufferDenom->next->value)));
    curr = bufferDenom->next->next;
    while(curr!=NULL) {
      printf("  mpfr_mul (tmp%d, tmp%d, tmp%d, MPFR_RNDN);\n", toReturn, toReturn, *((int *)(curr->value)));
      curr = curr->next;
    }
    printf("  mpfr_div (%s, %s, tmp%d, MPFR_RNDN);\n", resName, resName, toReturn);
  }
 
  freeChain(bufferNum, freeIntPtr);
  freeChain(bufferDenom, freeIntPtr);
  freeChain(numerator, free_memory);
  freeChain(denominator, free_memory);  
  return toReturn;
}

int constantImplementer(node *c, int gamma0, char *resName,  int counter) {
  char tmpName[10] = "tmp";
  int toReturn;

  switch (c->nodeType) {
  case ADD:
  case SUB:
    /* toReturn = implementAddSub(c, gamma0, resName, counter); */
    break;
  case MUL:
  case DIV:
    toReturn = implementDivMul(c, gamma0, resName, counter);
    break;
  case POW:
    break;

  case CONSTANT:
    printf("  mpfr_set_prec (%s, prec+%d);\n", resName, gamma0);
    if (mpfr_integer_p(*(c->value)) && mpfr_fits_ulong_p(*(c->value), MPFR_RNDN)) {
      printf("  mpfr_set_ui (%s, %lu, MPFR_RNDN);\n", resName, mpfr_get_ui(*(c->value), GMP_RNDN));
    }
    else if (mpfr_integer_p(*(c->value)) && mpfr_fits_slong_p(*(c->value), MPFR_RNDN)) {
      printf("  mpfr_set_si (%s, %ld, MPFR_RNDN);\n", resName, mpfr_get_si(*(c->value), GMP_RNDN));
    }
    else {
      mpfr_printf("  mpfr_set_str (%s, \"%RNb\", 2, MPFR_RNDN);\n", resName, *(c->value));
    }
    toReturn = counter;
    break;
  case NEG:
    toReturn = constantImplementer(c->child1, gamma0, resName, counter);
    printf("  mpfr_neg (%s, %s, MPFR_RNDN);\n", resName, resName);
    break;
  case ABS:
    toReturn = constantImplementer(c->child1, gamma0, resName, counter);
    printf("  mpfr_abs (%s, %s, MPFR_RNDN);\n", resName, resName);
    break;
  case DOUBLE:        
    break;
  case DOUBLEDOUBLE:
    break;
  case TRIPLEDOUBLE:
    break;
  case DOUBLEEXTENDED:
    break;
  case CEIL:
    break;
  case FLOOR:
    break;
  case PI_CONST:
    printf("  mpfr_set_prec (%s, prec+%d);\n", resName, gamma0);
    printf("  mpfr_const_pi (%s, MPFR_RNDN);\n", resName);
    toReturn = counter;
    break;
  case SINGLE:
    break;
  case NEARESTINT:
    break;

  case SQRT:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "sqrt", gamma0, resName, tmpName, counter+1);
    break;
  case EXP:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "exp", gamma0, resName, tmpName, counter+1);
    break;
  case LOG:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "log", gamma0, resName, tmpName, counter+1);
    break;
  case LOG_2:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "log2", gamma0, resName, tmpName, counter+1);
    break;
  case LOG_10:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "log10", gamma0, resName, tmpName, counter+1);
    break;
  case SIN:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "sin", gamma0, resName, tmpName, counter+1);
    break;
  case COS:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "cos", gamma0, resName, tmpName, counter+1);
    break;
  case TAN:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "tan", gamma0, resName, tmpName, counter+1);
    break;
  case ASIN:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "asin", gamma0, resName, tmpName, counter+1);
    break;
  case ACOS:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "acos", gamma0, resName, tmpName, counter+1);
    break;
  case ATAN:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "atan", gamma0, resName, tmpName, counter+1);
    break;
  case SINH:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "sinh", gamma0, resName, tmpName, counter+1);
    break;
  case COSH:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "cosh", gamma0, resName, tmpName, counter+1);
    break;
  case TANH:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "tanh", gamma0, resName, tmpName, counter+1);
    break;
  case ASINH:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "asinh", gamma0, resName, tmpName, counter+1);
    break;
  case ACOSH:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "acosh", gamma0, resName, tmpName, counter+1);
    break;
  case ATANH:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "atanh", gamma0, resName, tmpName, counter+1);
    break;
  case ERF:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "erf", gamma0, resName, tmpName, counter+1);
    break;
  case ERFC:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "erfc", gamma0, resName, tmpName, counter+1);
    break;
  case LOG_1P:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "log1p", gamma0, resName, tmpName, counter+1);
    break;
  case EXP_M1:
    sprintf(tmpName+3, "%d", counter+1);
    toReturn = unaryFunctionCase(c->nodeType, c->child1, "expm1", gamma0, resName, tmpName, counter+1);
    break;
  case LIBRARYFUNCTION:
    break;

  default:
    printMessage(1, "Unknown identifier (%d) in the tree\n", c->nodeType);
  }

  return toReturn;
}
