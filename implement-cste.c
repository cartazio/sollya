#include <string.h>
#include <stdlib.h>
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
   - mpfr_ui_pow_ui(var1, valui, valui2, MPFR_RNDN);
   - mpfr_pow_ui(var1, var2, valui, MPFR_RNDN);
   - mpfr_root(var1, var2, valui, MPFR_RNDN);
   - valstr(var1, prec); for a library constant
*/
#define INIT2 0
#define SETPREC 1
#define CONSTANTFUNC 2
#define UNARYFUNC 3
#define BINARYFUNC 4
#define SETUI 5
#define SETSI 6
#define SETSTR 7
#define UIPOWUI 8
#define POWUI 9
#define ROOT 10
#define LIBRARYCONST 11

struct implementCsteInstruction {
  int type;
  char var1[BUFFERSIZE];
  char var2[BUFFERSIZE];
  char var3[BUFFERSIZE];
  char name[BUFFERSIZE];
  long int prec;
  unsigned long int uival;
  unsigned long int uival2;
  long sival;
  char *strval;
};

void free_implementCsteInstruction(void *instr) {
  if ( ((struct implementCsteInstruction *)instr)->strval!=NULL )
    free( ((struct implementCsteInstruction *)instr)->strval);
  free( (struct implementCsteInstruction *)instr);
  return;
}

void *copy_implementCsteInstructions(void *instr) {
  struct implementCsteInstruction *newInstr;
  newInstr = safeMalloc(sizeof(struct implementCsteInstruction));
  newInstr->type = ((struct implementCsteInstruction *)instr)->type;
  strcpy(newInstr->var1, ((struct implementCsteInstruction *)instr)->var1);
  strcpy(newInstr->var2, ((struct implementCsteInstruction *)instr)->var2);
  strcpy(newInstr->var3, ((struct implementCsteInstruction *)instr)->var3);
  strcpy(newInstr->name, ((struct implementCsteInstruction *)instr)->name);
  newInstr->prec = ((struct implementCsteInstruction *)instr)->prec;
  newInstr->uival = ((struct implementCsteInstruction *)instr)->uival;
  newInstr->uival2 = ((struct implementCsteInstruction *)instr)->uival2;
  newInstr->sival = ((struct implementCsteInstruction *)instr)->sival;
  if(((struct implementCsteInstruction *)instr)->strval != NULL) {
    newInstr->strval = safeCalloc(1+strlen(((struct implementCsteInstruction *)instr)->strval) , sizeof(char));
    strcpy(newInstr->strval, ((struct implementCsteInstruction *)instr)->strval);
  }
  else newInstr->strval = NULL;
  return (void *)newInstr;
}

void fprintInstruction(FILE *output, struct implementCsteInstruction instr) {
  const char init_string[]="mpfr_init2";
  const char setprec_string[]="mpfr_set_prec";
  const char *ptr;

  switch (instr.type) {
  case INIT2:
  case SETPREC:
  case LIBRARYCONST:
    if (instr.type == INIT2) ptr=init_string;
    else if (instr.type == SETPREC) ptr=setprec_string;
    else ptr = instr.strval;

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
    fprintf(output, "  mpfr_set_str (%s, \"%s\", 2, MPFR_RNDN);\n", instr.var1, instr.strval);
    break;
  case UIPOWUI:
    fprintf(output, "  mpfr_ui_pow_ui (%s, %lu, %lu, MPFR_RNDN);\n", instr.var1, instr.uival, instr.uival2);
    break;
  case POWUI:
    fprintf(output, "  mpfr_pow_ui (%s, %s, %lu, MPFR_RNDN);\n", instr.var1, instr.var2, instr.uival);
    break;
  case ROOT:
    fprintf(output, "  mpfr_root (%s, %s, %lu, MPFR_RNDN);\n", instr.var1, instr.var2, instr.uival);
    break;
  default: 
    fprintf(stderr, "Unknown instruction %d\n", instr.type);
  }
  return;
}

void constructName(char *res, int counter) {
  if (counter==0) strcpy(res, "y");
  else sprintf(res, "tmp%d", counter);
  return;
}

/* A program is given by a list of instructions, the index number of the first
   unused temporary variable in the program, the number of temporary variables
   used by the program and a list giving, for each temporary variable, the
   maximum of the precision that it takes.      
*/
struct implementCsteProgram {
  chain *instructions;
  int counter;
  int maxcounter;
  chain *precisions;
};

typedef struct implementCsteCouple {
  int var;
  long int prec;
} couple;

couple *makeCouple(int var, long int prec) {
  couple *res = safeMalloc(sizeof(couple));
  res->var=var;
  res->prec=prec;
  return res;
}

/* This functions looks for the variable var in program. If it is absent, it
   adds the couple (var, prec) to program->precisions. If there already is an
   occurence (var, prec2) in program->precisions, it compares prec and prec2
   and keeps the largest.
*/
void appendPrecisionProg(int var, long int prec, struct implementCsteProgram *program) {
  chain *curr;
  int test = 0;
  curr = program->precisions;
  while ( (curr != NULL) && (!test) ) {
    if ( ((couple *)(curr->value))->var == var) {
      test = 1;
      if (prec > ((couple *)(curr->value))->prec) ((couple *)(curr->value))->prec = prec;
    }
    else curr = curr->next;
  }
  if (!test) program->precisions = addElement(program->precisions, makeCouple(var, prec));
  return;
}

void incrementProgramCounter(struct implementCsteProgram *program) {
  program->counter++;
  if (program->counter >= program->maxcounter) program->maxcounter = program->counter;
  return;
}

/* These constructors allow one for adding an instruction to the end of a */
/* give program.                                                          */
void appendInit2Prog(int var1, long int prec, struct implementCsteProgram *program) {
  struct implementCsteInstruction *instr;
  instr = safeMalloc(sizeof(struct implementCsteInstruction));
  instr->type = INIT2;
  constructName(instr->var1, var1);
  strcpy(instr->var2, "");
  strcpy(instr->var3, "");
  strcpy(instr->name, "");
  instr->prec = prec;
  instr->strval = NULL;
  program->instructions = addElement(program->instructions, instr);
  return;
}

void appendSetprecProg(int var1, long int prec, struct implementCsteProgram *program) {
  struct implementCsteInstruction *instr;
  instr = safeMalloc(sizeof(struct implementCsteInstruction));
  instr->type = SETPREC;
  constructName(instr->var1, var1);
  strcpy(instr->var2, "");
  strcpy(instr->var3, "");
  strcpy(instr->name, "");
  instr->prec = prec;
  instr->strval = NULL;
  program->instructions = addElement(program->instructions, instr);
  appendPrecisionProg(var1, prec, program);
  return;
}

void appendConstantfuncProg(char *name, int var1, struct implementCsteProgram *program) {
  struct implementCsteInstruction *instr;
  instr = safeMalloc(sizeof(struct implementCsteInstruction));
  instr->type = CONSTANTFUNC;
  constructName(instr->var1, var1);
  strcpy(instr->var2, "");
  strcpy(instr->var3, "");
  strcpy(instr->name, name);
  instr->strval = NULL;
  program->instructions = addElement(program->instructions, instr);
  return;
}

void appendUnaryfuncProg(char *name, int var1, int var2, struct implementCsteProgram *program) {
  struct implementCsteInstruction *instr;
  instr = safeMalloc(sizeof(struct implementCsteInstruction));
  instr->type = UNARYFUNC;
  constructName(instr->var1, var1);
  constructName(instr->var2, var2);
  strcpy(instr->var3, "");
  strcpy(instr->name, name);
  instr->strval = NULL;
  program->instructions = addElement(program->instructions, instr);
  return;
}

void appendBinaryfuncProg(char *name, int var1, int var2, int var3, struct implementCsteProgram *program) {
  struct implementCsteInstruction *instr;
  instr = safeMalloc(sizeof(struct implementCsteInstruction));
  instr->type = BINARYFUNC;
  constructName(instr->var1, var1);
  constructName(instr->var2, var2);
  constructName(instr->var3, var3);
  strcpy(instr->name, name);
  instr->strval = NULL;
  program->instructions = addElement(program->instructions, instr);
  return;
}

void appendSetuiProg(int var1, unsigned long int val, struct implementCsteProgram *program) {
  struct implementCsteInstruction *instr;
  instr = safeMalloc(sizeof(struct implementCsteInstruction));
  instr->type = SETUI;
  constructName(instr->var1, var1);
  strcpy(instr->var2, "");
  strcpy(instr->var3, "");
  strcpy(instr->name, "");
  instr->uival = val;
  instr->strval = NULL;
  program->instructions = addElement(program->instructions, instr);
  return;
}

void appendSetsiProg(int var1, long int val, struct implementCsteProgram *program) {
  struct implementCsteInstruction *instr;
  instr = safeMalloc(sizeof(struct implementCsteInstruction));
  instr->type = SETSI;
  constructName(instr->var1, var1);
  strcpy(instr->var2, "");
  strcpy(instr->var3, "");
  strcpy(instr->name, "");
  instr->sival = val;
  instr->strval = NULL;
  program->instructions = addElement(program->instructions, instr);
  return;
}

void appendSetstrProg(int var1, mpfr_t val, struct implementCsteProgram *program) {
  struct implementCsteInstruction *instr;
  instr = safeMalloc(sizeof(struct implementCsteInstruction));
  instr->type = SETSTR;
  constructName(instr->var1, var1);
  strcpy(instr->var2, "");
  strcpy(instr->var3, "");
  strcpy(instr->name, "");
  instr->strval = safeCalloc(mpfr_get_prec(val)+32, sizeof(char)); /* should be sufficient to store the string representing val in binary */
  mpfr_sprintf(instr->strval, "%RNb", val);
  program->instructions = addElement(program->instructions, instr);
  return;
}

void appendUipowui(int var1, unsigned long int val1, unsigned long int val2, struct implementCsteProgram *program) {
  struct implementCsteInstruction *instr;
  instr = safeMalloc(sizeof(struct implementCsteInstruction));
  instr->type = UIPOWUI;
  constructName(instr->var1, var1);
  strcpy(instr->var2, "");
  strcpy(instr->var3, "");
  strcpy(instr->name, "");
  instr->uival = val1;
  instr->uival2 = val2;
  instr->strval = NULL;
  program->instructions = addElement(program->instructions, instr);
  return;  
}

void appendPowuiProg(int var1, int var2, unsigned long int val, struct implementCsteProgram *program) {
  struct implementCsteInstruction *instr;
  instr = safeMalloc(sizeof(struct implementCsteInstruction));
  instr->type = POWUI;
  constructName(instr->var1, var1);
  constructName(instr->var2, var2);
  strcpy(instr->var3, "");
  strcpy(instr->name, "");
  instr->uival = val;
  instr->strval = NULL;
  program->instructions = addElement(program->instructions, instr);
  return;  
}

void appendRootProg(int var1, int var2, unsigned long int val, struct implementCsteProgram *program) {
  struct implementCsteInstruction *instr;
  instr = safeMalloc(sizeof(struct implementCsteInstruction));
  instr->type = ROOT;
  constructName(instr->var1, var1);
  constructName(instr->var2, var2);
  strcpy(instr->var3, "");
  strcpy(instr->name, "");
  instr->uival = val;
  instr->strval = NULL;
  program->instructions = addElement(program->instructions, instr);
  return;  
}

void appendLibraryConstantProg(node *c, int gamma0, struct implementCsteProgram *program) {
  struct implementCsteInstruction *instr;
  instr = safeMalloc(sizeof(struct implementCsteInstruction));
  instr->type = LIBRARYCONST;
  constructName(instr->var1, program->counter);
  instr->prec = gamma0;
  strcpy(instr->var2, "");
  strcpy(instr->var3, "");
  strcpy(instr->name, "");
  instr->strval = safeCalloc(strlen(c->libFun->functionName)+1, sizeof(char));
  strcpy(instr->strval, c->libFun->functionName);

  program->instructions = addElement(program->instructions, instr);
  return;  
}

/* Prototype of the recursive function */
void constantImplementer(node *c, int gamma0, struct implementCsteProgram *program);

/* Main function */
void implementCste(node *c) {
  int i, test;
  const char name[] = "something";
  FILE *output = stdout;
  struct implementCsteProgram program;
  chain *curr;
  
  program.instructions = NULL;
  program.counter = 0;
  program.maxcounter = 0;
  program.precisions = NULL;

  constantImplementer(c, 0, &program);

  /* reverse the chain */
  curr = copyChain(program.instructions, copy_implementCsteInstructions);
  freeChain(program.instructions, free_implementCsteInstruction);
  program.instructions = curr;

  curr = program.precisions;
  while(curr != NULL) {
    appendInit2Prog( ((couple *)(curr->value))->var, ((couple *)(curr->value))->prec, &program);
    curr = curr->next;
  }

  fprintf(output, "void\n");
  fprintf(output, "mpfr_const_%s (mpfr_ptr y, mp_prec_t prec)\n", name);
  fprintf(output, "{\n");
  if(program.maxcounter>=2) fprintf(output, "  /* Declarations */\n");
  for(i=1; i<=program.maxcounter-1; i++) fprintf(output, "  mpfr_t tmp%d;\n", i);
  if(program.maxcounter>=2) fprintf(output, "\n");
  fprintf(output, "  /* Initializations */\n");

  test = 1;
  curr=program.instructions;
  while(curr!=NULL) {
    if (test 
        && ( ((struct implementCsteInstruction *)(curr->value))->type != INIT2 )
        ) {
      fprintf(output, "\n");
      fprintf(output, "  /* Core */\n");
      test = 0;
    }
    fprintInstruction(output, *(struct implementCsteInstruction *)(curr->value));
    curr = curr->next;
  }

  if(program.maxcounter>=2) {
    fprintf(output, "\n");
    fprintf(output, "  /* Cleaning stuff */\n");
  }
  for(i=1; i<=program.maxcounter-1; i++) fprintf(output, "  mpfr_clear(tmp%d);\n", i);
  fprintf(output, "}\n");

  freeChain(program.instructions, free_implementCsteInstruction);
  freeChain(program.precisions, free);
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

/* Returns the maximal exponent of a number among those contained in x */
mp_exp_t mpfi_max_exp(mpfi_t x) {
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

/* Return the smallest exponent among the exponents of the numbers contained
   in x. If 0 \in x, it returns NULL, else it returns a valid pointer E such
   that *E is the minimal exponent. */
mp_exp_t *mpfi_min_exp(mpfi_t x) {
  mpfr_t u,v;
  mp_exp_t Eu, Ev;
  mp_exp_t *E = NULL;
  mp_prec_t prec;

  prec = mpfi_get_prec(x);
  mpfr_init2(u, prec);
  mpfr_init2(v, prec);

  mpfi_get_left(u, x);
  mpfi_get_right(v, x);
  
  if (mpfr_sgn(u)*mpfr_sgn(v)>0) {
    E = safeMalloc(sizeof(mp_exp_t));
    Eu = mpfr_get_exp(u);
    Ev = mpfr_get_exp(v);
    *E = (Eu<=Ev) ? Eu : Ev;
  }

  mpfr_clear(u);
  mpfr_clear(v);
  return E;
}
 
/* Let a be the constant given by the expression cste and f the function with */
/* node type nodeType. This functions generates code for the implementation   */
/* of f(a) in precision prec+gamma0, the result being stored in resName.      */
int unaryFunctionCase(int nodeType, node *cste, char *functionName, int gamma0, struct implementCsteProgram *program) {
  mpfi_t a, b, u, v, tmp;
  mpfr_t alpha, beta;
  mp_prec_t prec = getToolPrecision();
  node *func, *deriv;
  int gamma;
  int counter;

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
  if (mpfi_has_zero(b)) {
    changeToWarningMode();
    fprintf(stderr, "Error in implementconstant: the following expression seems to be exactly zero:\n");
    fprintTree(stderr, makeUnary(copyTree(cste), nodeType));
    fprintf(stderr, "\nIf it is not exactly zero, increasing prec should solve the issue.\nAbort.\n");
    restoreMode();
    recoverFromError();
  }

  mpfi_div(u, a, b);
  evaluateInterval(tmp, deriv, NULL, a); 
  mpfi_mul(v, u, tmp);

  gamma = 2+mpfi_max_exp(v)-1;
  do {
    gamma++;
    mpfr_set_ui(beta, 1, GMP_RNDU);
    mpfr_div_2si(beta, beta, gamma+gamma0, GMP_RNDU);
    mpfr_ui_sub(alpha, 1, beta, GMP_RNDD);
    mpfr_add_ui(beta, beta, 1, GMP_RNDU);
    mpfi_interv_fr(tmp, alpha, beta);

    mpfi_mul(tmp, a, tmp);
    evaluateInterval(tmp, deriv, NULL, tmp);
    mpfi_mul(v, u, tmp);
  } while (gamma < 2+mpfi_max_exp(v));
  
  counter = program->counter;
  incrementProgramCounter(program);
  constantImplementer(cste, gamma0+gamma, program);
  program->counter = counter;
  appendSetprecProg(counter, gamma0+2, program);
  appendUnaryfuncProg(functionName, counter, counter+1, program);

  mpfi_clear(a);
  mpfi_clear(b);
  mpfi_clear(u);
  mpfi_clear(v);
  mpfi_clear(tmp);
  mpfr_clear(alpha);
  mpfr_clear(beta);
  free_memory(func);
  free_memory(deriv);
  return;
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

int implementDivMul(node *c, int gamma0, struct implementCsteProgram *program) {
  chain *numerator = NULL;
  chain *denominator = NULL;
  chain *curr;
  chain *bufferNum, *bufferDenom;
  int log2n, n;
  int *tmp;
  int counter;

  normalizeDivMul(c, &numerator, &denominator);

  n = lengthChain(numerator) + lengthChain(denominator);
  log2n = ceil_log2n(n);
  counter = program->counter;
  incrementProgramCounter(program);

  curr = numerator;
  bufferNum = NULL;
  while(curr!=NULL) {
    tmp = safeMalloc(sizeof(int));
    *tmp = program->counter;
    bufferNum = addElement(bufferNum, tmp);
    constantImplementer(curr->value, gamma0+2+log2n, program);
    curr = curr->next;
  }
  curr = denominator;
  bufferDenom = NULL;
  while(curr!=NULL) {
    tmp = safeMalloc(sizeof(int));
    *tmp = program->counter;
    bufferDenom = addElement(bufferDenom, tmp);
    constantImplementer(curr->value, gamma0+2+log2n, program);
    curr = curr->next;
  }

  program->counter = counter;
  if ( (lengthChain(numerator)==1) && (lengthChain(denominator)==1) ) {
    appendSetprecProg(counter, gamma0+2+log2n, program);
    appendBinaryfuncProg("mpfr_div", counter, *((int *)(bufferNum->value)), *((int *)(bufferDenom->value)), program);
  }
  else if (lengthChain(numerator)==1) {
    appendSetprecProg(counter, gamma0+2+log2n, program);
    appendBinaryfuncProg("mpfr_mul", counter, *((int *)(bufferDenom->value)), *((int *)(bufferDenom->next->value)), program);
    curr = bufferDenom->next->next;
    while(curr!=NULL) {
      appendBinaryfuncProg("mpfr_mul", counter, counter, *((int *)(curr->value)), program);
      curr = curr->next;
    }
    appendBinaryfuncProg("mpfr_div", counter, *((int *)(bufferNum->value)), counter, program);
  }
  else if (lengthChain(denominator)<=1) {
    appendSetprecProg(counter, gamma0+2+log2n, program);
    appendBinaryfuncProg("mpfr_mul", counter, *((int *)(bufferNum->value)), *((int *)(bufferNum->next->value)), program);
    curr = bufferNum->next->next;
    while(curr!=NULL) {
      appendBinaryfuncProg("mpfr_mul", counter, counter, *((int *)(curr->value)), program);
      curr = curr->next;
    }
    if (lengthChain(denominator)==1)
      appendBinaryfuncProg("mpfr_div", counter, counter, *((int *)(bufferDenom->value)), program);
  }
  else {
    incrementProgramCounter(program);
    appendSetprecProg(counter, gamma0+2+log2n, program);
    appendBinaryfuncProg("mpfr_mul", counter, *((int *)(bufferNum->value)), *((int *)(bufferNum->next->value)), program);
    curr = bufferNum->next->next;
    while(curr!=NULL) {
      appendBinaryfuncProg("mpfr_mul", counter, counter, *((int *)(curr->value)), program);
      curr = curr->next;
    }
    appendSetprecProg(program->counter, gamma0+2+log2n, program);
    appendBinaryfuncProg("mpfr_mul", program->counter,  *((int *)(bufferDenom->value)), *((int *)(bufferDenom->next->value)), program);
    curr = bufferDenom->next->next;
    while(curr!=NULL) {
      appendBinaryfuncProg("mpfr_mul", program->counter, program->counter, *((int *)(curr->value)), program);
      curr = curr->next;
    }
    appendBinaryfuncProg("mpfr_div", counter, counter, program->counter, program);
  }
 
  program->counter = counter;
  freeChain(bufferNum, freeIntPtr);
  freeChain(bufferDenom, freeIntPtr);
  freeChain(numerator, (void (*)(void *))free_memory);
  freeChain(denominator, (void (*)(void *))free_memory);  
  return;
}

int implementAddSub(node *c, int gamma0, struct implementCsteProgram *program) {
  mpfi_t y, a, b, tmp, tmp2;
  mp_exp_t *Ea, *Eb, *Ey;
  int tmpa, tmpb;
  mp_prec_t prec;
  int counter;

  prec = getToolPrecision();
  mpfi_init2(y, prec);
  mpfi_init2(a, prec);
  mpfi_init2(b, prec);
  mpfi_init2(tmp, prec);
  mpfi_init2(tmp2, prec);

  
  evaluateInterval(y, c, NULL, y);
  if (mpfi_has_zero(y)) {
    changeToWarningMode();
    fprintf(stderr, "Error in implementconstant: the following expression seems to be exactly zero:\n");
    fprintTree(stderr, c);
    fprintf(stderr, "\nIf it is not exactly zero, increasing prec should solve the issue.\nAbort.\n");
    restoreMode();
    recoverFromError();
  }
  evaluateInterval(a, c->child1, NULL, a);
  evaluateInterval(b, c->child2, NULL, b);

  counter = program->counter;
  incrementProgramCounter(program);

  tmpa = program->counter;
  mpfi_div(tmp, y, a); mpfi_div_ui(tmp, tmp, 3);
  Ea = mpfi_min_exp(tmp);
  if (Ea==NULL) {
    printMessage(0, "Unexpected error. Aborting\n");
    recoverFromError();
  }
  constantImplementer(c->child1, gamma0+1-*Ea, program);

  tmpb = program->counter;
  mpfi_div(tmp, y, b); mpfi_div_ui(tmp, tmp, 3);
  Eb = mpfi_min_exp(tmp);
  if (Eb==NULL) {
    printMessage(0, "Unexpected error. Aborting\n");
    recoverFromError();
  }
  constantImplementer(c->child2, gamma0+1-*Eb, program);

  mpfi_abs(tmp, a);
  mpfi_abs(tmp2, b);
  mpfi_add(tmp, tmp, tmp2);
  mpfi_div(tmp, y, tmp);
  mpfi_div_ui(tmp, tmp, 3);
  Ey = mpfi_min_exp(tmp);
  if (Ey==NULL) {
    printMessage(0, "Unexpected error. Aborting\n");
    recoverFromError();
  }
  appendSetprecProg(counter, gamma0+2-*Ey, program);
  if (c->nodeType==ADD)
    appendBinaryfuncProg("mpfr_add", counter, tmpa, tmpb, program);
  else
    appendBinaryfuncProg("mpfr_sub", counter, tmpa, tmpb, program);

  program->counter = counter;
  free(Ea);
  free(Eb);
  free(Ey);
  mpfi_clear(y);
  mpfi_clear(a);
  mpfi_clear(b);
  mpfi_clear(tmp);
  mpfi_clear(tmp2);
  return;
}

void implementPow(node *c, int gamma0, struct implementCsteProgram *program) {
  int log2p, p, tmpNumber, counter;
  node *tmpNode;
  mpfr_t tmp;

  counter = program->counter;
  if ( (c->child1->nodeType==CONSTANT) 
       && mpfr_integer_p(*(c->child1->value))
       && mpfr_fits_ulong_p(*(c->child1->value), GMP_RNDN)
       && (c->child2->nodeType==CONSTANT) 
       && mpfr_integer_p(*(c->child2->value))
       && mpfr_fits_ulong_p(*(c->child2->value), GMP_RNDN)) { /* Case n^p */
    appendSetprecProg(counter, gamma0, program);
    appendUipowui(counter, mpfr_get_ui(*(c->child1->value), GMP_RNDN), mpfr_get_ui(*(c->child2->value), GMP_RNDN), program);
    program->counter = counter;
    return;
  }

  if ( (c->child2->nodeType==CONSTANT) 
       && mpfr_integer_p(*(c->child2->value))
       && mpfr_fits_ulong_p(*(c->child2->value), GMP_RNDN) ) { /* Case x^p */
    p = mpfr_get_ui(*(c->child2->value), GMP_RNDN);
    log2p = ceil_log2n(p);
    incrementProgramCounter(program);
    constantImplementer(c->child1, gamma0+log2p+3, program);
    appendSetprecProg(counter, gamma0+2, program);
    appendPowuiProg(counter, counter+1, p, program);
    program->counter = counter;
    return;    
  }

  if ( (c->child2->nodeType==DIV)
       && (c->child2->child1->nodeType==CONSTANT)
       && (mpfr_cmp_ui(*(c->child2->child1->value), 1)==0)
       && (c->child2->child2->nodeType==CONSTANT)
       && mpfr_integer_p(*(c->child2->child2->value))
       && mpfr_fits_ulong_p(*(c->child2->child2->value), GMP_RNDN)
       ) { /* Case x^(1/p) (note that this does not handle the case when p=2^k */
    p = mpfr_get_ui(*(c->child2->child2->value), GMP_RNDN);
    log2p = ceil_log2n(p);
    incrementProgramCounter(program);
    constantImplementer(c->child1, gamma0-log2p+3, program);
    appendSetprecProg(counter, gamma0+2, program);
    appendRootProg(counter, counter+1, p, program);
    program->counter = counter;
    return; 
  }

  if (c->child2->nodeType==CONSTANT) {
    mpfr_init2(tmp, 64);
    if ( (mpfr_ui_div(tmp, 1, *(c->child2->value), GMP_RNDN) == 0)
         && mpfr_integer_p(tmp)
         && mpfr_fits_ulong_p(tmp, GMP_RNDN)
         ) { /* Case x^(1/p) where p is a power of 2 */
      p = mpfr_get_ui(tmp, GMP_RNDN);
      log2p = ceil_log2n(p);
      incrementProgramCounter(program);
      constantImplementer(c->child1, gamma0-log2p+3, program);
      appendSetprecProg(counter, gamma0+2, program);
      appendRootProg(counter, counter+1, p, program);
      program->counter = counter;
      mpfr_clear(tmp);
      return; 
    }
    mpfr_clear(tmp);
  }

  /* else... case x^y with x possibly integer. Handled as exp(y*ln(x)) */
  tmpNode = makeExp(makeMul(copyTree(c->child2), makeLog(copyTree(c->child1))));
  constantImplementer(tmpNode, gamma0, program);
  free_memory(tmpNode);
  program->counter = counter;
  return;
}

void implementCsteCase(node *c, int gamma0, struct implementCsteProgram *program) {
  appendSetprecProg(program->counter, gamma0, program);
  if (mpfr_integer_p(*(c->value)) && mpfr_fits_ulong_p(*(c->value), GMP_RNDN)) {
    appendSetuiProg(program->counter, mpfr_get_ui(*(c->value), GMP_RNDN), program);
  }
  else if (mpfr_integer_p(*(c->value)) && mpfr_fits_slong_p(*(c->value), GMP_RNDN)) {
    appendSetsiProg(program->counter, mpfr_get_si(*(c->value), GMP_RNDN), program);
  }
  else {
    appendSetstrProg(program->counter, *(c->value), program);
  }
}

void constantImplementer(node *c, int gamma0, struct implementCsteProgram *program) {
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
    constantImplementer(c->child1, gamma0, program);
    appendUnaryfuncProg("mpfr_neg", program->counter, program->counter, program);
    break;
  case ABS:
    constantImplementer(c->child1, gamma0, program);
    appendUnaryfuncProg("mpfr_abs", program->counter, program->counter, program);
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
    appendSetprecProg(program->counter, gamma0, program);
    appendConstantfuncProg("mpfr_const_pi", program->counter, program);
    break;
  case LIBRARYCONSTANT:
    appendPrecisionProg(program->counter, gamma0, program);
    appendLibraryConstantProg(c, gamma0, program);
    break;
  case SINGLE:
    break;
  case NEARESTINT:
    break;

  case SQRT:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_sqrt", gamma0, program);
    break;
  case EXP:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_exp", gamma0, program);
    break;
  case LOG:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_log", gamma0, program);
    break;
  case LOG_2:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_log2", gamma0, program);
    break;
  case LOG_10:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_log10", gamma0, program);
    break;
  case SIN:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_sin", gamma0, program);
    break;
  case COS:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_cos", gamma0, program);
    break;
  case TAN:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_tan", gamma0, program);
    break;
  case ASIN:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_asin", gamma0, program);
    break;
  case ACOS:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_acos", gamma0, program);
    break;
  case ATAN:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_atan", gamma0, program);
    break;
  case SINH:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_sinh", gamma0, program);
    break;
  case COSH:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_cosh", gamma0, program);
    break;
  case TANH:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_tanh", gamma0, program);
    break;
  case ASINH:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_asinh", gamma0, program);
    break;
  case ACOSH:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_acosh", gamma0, program);
    break;
  case ATANH:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_atanh", gamma0, program);
    break;
  case ERF:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_erf", gamma0, program);
    break;
  case ERFC:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_erfc", gamma0, program);
    break;
  case LOG_1P:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_log1p", gamma0, program);
    break;
  case EXP_M1:
    unaryFunctionCase(c->nodeType, c->child1, "mpfr_expm1", gamma0, program);
    break;
  case LIBRARYFUNCTION:
    break;

  default:
    printMessage(1, "Unknown identifier (%d) in the tree\n", c->nodeType);
  }

  incrementProgramCounter(program);
  return;
}
