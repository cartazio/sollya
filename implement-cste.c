#include <stdio.h>
#include "expression.h"
#include "infnorm.h"
#include "general.h"

#define BUFFERSIZE 64

/* Possible instructions:
   - mpfr_init2(var1, prec)
   - mpfr_set_prec(var1, prec)
   - 0ary function: name(var1, MPFR_RNDN)
   - unary function: name(var1, var2, MPFR_RNDN)
   - binary function: name(var1, var2, var3, MPFR_RNDN)
   - mpfr_set_ui(var1, valui, MPFR_RNDN)
   - mpfr_set_si(var1, valsi, MPFR_RNDN)
   - mpfr_set_str(var1, valstr, 2, MPFR_RNDN)
*/
#define INIT2 0
#define SETPREC 1
#define CONSTANTFUNC 2
#define UNARYFUNC 3
#define BINARYFUNC 4
#define SETUI 5
#define SETSI 6
#define SETSTR 7

struct implementCsteInstruction {
  int type;
  char var1[BUFFERSIZE];
  char var2[BUFFERSIZE];
  char var3[BUFFERSIZE];
  char name[BUFFERSIZE]
  long int prec;
  unsigned long int uival;
  long sival;
  char *strval
};

void fprintInstruction(FILE *output, struct implementCsteInstruction instr) {
  const char init_string[]="mpfr_init2";
  const char setprec_string[]="mpfr_set_prec";
  char *ptr;

  switch (instr.type) {
  case INIT2:
  case SETPREC:
    if (instr.type == INIT2) ptr=init_string;
    else ptr=setprec_string;

    if (instr.prec > 0)
      fprintf(output, "  %s (%s, prec+%d);\n", ptr, instr.var1, instr.prec);
    else if (instr.prec == 0)
      fprintf(output, "  %s (%s, prec);\n", ptr, instr.var1);
    else {
      fprintf(output, "  if (prec >= %d+MPFR_PREC_MIN)\n", -instr.prec);
      fprintf(output, "  {\n");
      fprintf(output, "    %s (%s, prec-%d);\n", ptr, instr.var1, -instr.prec);
      fprintf(output, "  }\n");
      fprintf(output, "  else\n");
      fprintf(output, "  {\n");
      fprintf(output, "    %s (%s, MPFR_PREC_MIN);\n", ptr, instr.var1);
      fprintf(output, "  }\n");
    }
    break;
  case CONSTANTFUNC:
    fprintf(output, "  %s (%s, MPFR_RNDN);\n", instr.name, instr.var1);
    break;
  case UNARYFUNC:
    fprintf(output, "  %s (%s, %s, MPFR_RNDN);\n", instr.name, instr.var1, instr.var2);
    break;
  case BINARYFUNC:
    fprintf(output, "  %s (%s, %s, %s, MPFR_RNDN);\n", instr.name, instr.var1, instr.var2, instr.var3);
    break;
  case SETUI:
    fprintf(output, "  mpfr_set_ui (%s, %lu, MPFR_RNDN);\n", instr.var1, instr.uival);
    break;
  case SETSI:
    fprintf(output, "  mpfr_set_si (%s, %ld, MPFR_RNDN);\n", instr.var1, instr.sival);
    break;
  case SETSTR:
    fprintf(output, "  mpfr_set_str (%s, %s, 2, MPFR_RNDN);\n", instr.var1, instr.strval);
    break;
  default: 
    fprintf(stderr, "Unknown instruction %d\n", instr.type);
  }
  return;
}

void constructName(char *res, int counter) {
  if (counter==0) sprintf(res, "y");
  else sprintf(res, "tmp%d", counter);
  return;
}

/* A program is given by a list of instructions, the index number of the last
   temporary variable used in the program, and a list giving, for each temporary
   variable, the maximum of the precision that it takes.      
*/
struct implementCsteProgram {
  chain *instructions;
  int counter;
  chain *precisions;
};

int constantImplementer(node *c, int gamma0, char *resName, int counter);

void implementCste(node *c) {
  int i;
  const char name[] = "something";
  FILE *output = stdout;
  struct implementCsteProgram program;

  program.instructions = NULL;
  program.counter = 0;
  program.precisions = NULL;

  counter = constantImplementer(c, 0, program);

  fprintf(output, "void\n");
  fprintf(output, "mpfr_const_%s (mpfr_ptr y, mp_prec_t prec)\n", name);
  fprintf(output, "{\n");
  fprintf(output, "  /* Declarations */\n");
  fprintf(output, "  /* Initializations */\n");
  fprintf(output, "\n");
  fprintf(output, "  /* Core */\n");
  fprintf(output, "\n");
  fprintf(output, "  /* Cleaning stuff */\n");
  fprintf(output, "}\n");

  return;
}

int ceil_log2n(int p) {
  int n,log2p, test;
  n = p;
  log2p = 0;
  test = 1;
  /* Compute log2p such that 2^(log2p-1) <= p < 2^(log2p) */
  while (n>=1) { log2p++; if ((n%2)!=0) test=0; n = n/2;} 
  /* Adjust log2p in order to have 2^(log2p-1) < n <= 2^(log2p) */
  if(test) log2p--;

  return log2p;
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
  int log2n, n;
  char tmpName[10] = "tmp";
  int toReturn = counter;
  int *tmp;

  normalizeDivMul(c, &numerator, &denominator);

  n = lengthChain(numerator) + lengthChain(denominator);
  log2n = ceil_log2n(n);
  curr = numerator;
  bufferNum = NULL;
  while(curr!=NULL) {
    sprintf(tmpName+3, "%d", toReturn+1);
    tmp = safeMalloc(sizeof(int));
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
  freeChain(numerator, (void (*)(void *))free_memory);
  freeChain(denominator, (void (*)(void *))free_memory);  
  return toReturn;
}

int implementAddSub(node *c, int gamma0, char *resName, int counter) {
  mpfi_t y, a, b, tmp, tmp2;
  int tmpa, tmpb;
  int toReturn = counter;
  char tmpName[10] = "tmp";
  mp_prec_t prec;

  prec = getToolPrecision();
  mpfi_init2(y, prec);
  mpfi_init2(a, prec);
  mpfi_init2(b, prec);
  mpfi_init2(tmp, prec);
  mpfi_init2(tmp2, prec);

  evaluateInterval(y, c, NULL, y);
  evaluateInterval(a, c->child1, NULL, a);
  evaluateInterval(b, c->child2, NULL, b);

  tmpa = toReturn+1;
  sprintf(tmpName+3, "%d", toReturn+1);
  mpfi_div(tmp, y, a); mpfi_div_ui(tmp, tmp, 3);
  toReturn = constantImplementer(c->child1, gamma0+1-mpfi_get_exp(tmp), tmpName, toReturn+1);

  tmpb = toReturn+1;
  sprintf(tmpName+3, "%d", toReturn+1);
  mpfi_div(tmp, y, b); mpfi_div_ui(tmp, tmp, 3);
  toReturn = constantImplementer(c->child2, gamma0+1-mpfi_get_exp(tmp), tmpName, toReturn+1);

  mpfi_abs(tmp, a);
  mpfi_abs(tmp2, b);
  mpfi_add(tmp, tmp, tmp2);
  mpfi_div(tmp, y, tmp);
  mpfi_div_ui(tmp, tmp, 3);
  printf("  mpfr_set_prec (%s, prec+%d);\n", resName, gamma0+2-mpfi_get_exp(tmp));
  if (c->nodeType==ADD)
    printf("  mpfr_add (%s, tmp%d, tmp%d, MPFR_RNDN);\n", resName, tmpa, tmpb);
  else
    printf("  mpfr_sub (%s, tmp%d, tmp%d, MPFR_RNDN);\n", resName, tmpa, tmpb);

  mpfi_clear(y);
  mpfi_clear(a);
  mpfi_clear(b);
  mpfi_clear(tmp);
  mpfi_clear(tmp2);
  return toReturn;
}

int implementPow(node *c, int gamma0, char *resName, int counter) {
  int toReturn = counter;
  char tmpName[10] = "tmp";
  int log2p, p, tmpNumber;
  node *tmpNode;

  if ( (c->child1->nodeType==CONSTANT) 
       && mpfr_integer_p(*(c->child1->value))
       && mpfr_fits_ulong_p(*(c->child1->value), MPFR_RNDN)
       && (c->child2->nodeType==CONSTANT) 
       && mpfr_integer_p(*(c->child2->value))
       && mpfr_fits_ulong_p(*(c->child2->value), MPFR_RNDN)) { /* Case n^p */
    printf("  mpfr_ui_pow_ui(%s, %lu, %lu, MPFR_RNDN);\n", resName, mpfr_get_ui(*(c->child1->value), MPFR_RNDN), mpfr_get_ui(*(c->child2->value), MPFR_RNDN));
    return counter;
  }

  if ( (c->child2->nodeType==CONSTANT) 
       && mpfr_integer_p(*(c->child2->value))
       && mpfr_fits_ulong_p(*(c->child2->value), GMP_RNDN) ) { /* Case x^p */
    p = mpfr_get_ui(*(c->child2->value), GMP_RNDN);
    log2p = ceil_log2n(p);
    sprintf(tmpName+3, "%d", toReturn+1);
    tmpNumber = toReturn + 1;
    toReturn = constantImplementer(c->child1, gamma0+log2p+3, tmpName, toReturn+1);
    printf("  mpfr_set_prec (%s, prec+%d);\n", resName, gamma0+2);
    printf("  mpfr_pow_ui (%s, tmp%d, %d, MPFR_RNDN);\n", resName, tmpNumber, p);
    return toReturn;
  }

  if ( (c->child2->nodeType==DIV)
       && (c->child2->child1->nodeType==CONSTANT)
       && (mpfr_cmp_ui(*(c->child2->child1->value), 1)==0)
       && (c->child2->child2->nodeType==CONSTANT)
       && mpfr_integer_p(*(c->child2->child2->value))
       && mpfr_fits_ulong_p(*(c->child2->child2->value), GMP_RNDN)
       ) { /* Case x^(1/p) */
    p = mpfr_get_ui(*(c->child2->child2->value), GMP_RNDN);
    log2p = ceil_log2n(p);
    sprintf(tmpName+3, "%d", toReturn+1);
    tmpNumber = toReturn + 1;
    toReturn = constantImplementer(c->child1, gamma0-log2p+3, tmpName, toReturn+1);
    printf("  mpfr_set_prec (%s, prec+%d);\n", resName, gamma0+2);
    printf("  mpfr_root (%s, tmp%d, %d, MPFR_RNDN);\n", resName, tmpNumber, p);
    return toReturn;    
  }

  /* else... case x^y with x possibly integer. Handled as exp(y*ln(x)) */
  tmpNode = makeExp(makeMul(copyTree(c->child2), makeLog(copyTree(c->child1))));
  toReturn = constantImplementer(tmpNode, gamma0, resName, counter);
  free_memory(tmpNode);
  return toReturn;
}

void implementCsteCase(node *c, int gamma0, struct implementCsteProgram program) {
  struct implementCsteInstruction *instr;

  instr = safeMalloc(sizeof(instr));
  instr->type = SETPREC;
  constructName(instr->var1, program.counter);
  instr->prec = gamma0;
  addElement(program.instructions, instr);

  if (mpfr_integer_p(*(c->value)) && mpfr_fits_ulong_p(*(c->value), GMP_RNDN)) {
    instr = safeMalloc(sizeof(instr));
    instr->type = SETUI;
    constructName(instr->var1, program.counter);
    instr->uival = mpfr_get_ui(*(c->value), GMP_RNDN);
    addElement(program.instructions, instr);
  }
  else if (mpfr_integer_p(*(c->value)) && mpfr_fits_slong_p(*(c->value), GMP_RNDN)) {
    instr = safeMalloc(sizeof(instr));
    instr->type = SETSI;
    constructName(instr->var1, program.counter);
    instr->sival = mpfr_get_si(*(c->value), GMP_RNDN);
    addElement(program.instructions, instr);
  }
  else {
    instr = safeMalloc(sizeof(instr));
    instr->type = SETSTR;
    constructName(instr->var1, program.counter);
    instr->strval = safeCalloc(mpfr_get_prec(*(c->value))+32); /* should be sufficient to store the string representing *(c->value) in binary */
    mpfr_sprintf(instr->strval, "%RNb", *(c->value));
    addElement(program.instructions, instr);
  }
}

constantImplementer(node *c, int gamma0, struct implementCsteProgram program) {
  char tmpName[10] = "tmp";

  switch (c->nodeType) {
  case ADD:
  case SUB:
    implementAddSub(c, gamma0, program);
    break;
  case MUL:
  case DIV:
    implementDivMul(c, gamma0, program);
    break;
  case POW:
    implementPow(c, gamma0, program);
    break;
  case CONSTANT:
    implementCsteCase(c, gamma0, program);
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
