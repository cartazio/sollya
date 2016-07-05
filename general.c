/*

  Copyright 2007-2016 by

  Laboratoire de l'Informatique du Parallelisme,
  UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668,

  Laboratoire d'Informatique de Paris 6 - Équipe PEQUAN
  Sorbonne Universités
  UPMC Univ Paris 06
  UMR 7606, LIP6
  Boîte Courrier 169
  4, place Jussieu
  F-75252 Paris Cedex 05
  France

  and by

  Centre de recherche INRIA Sophia-Antipolis Mediterranee, equipe APICS,
  Sophia Antipolis, France.

  Contributors Ch. Lauter, S. Chevillard

  christoph.lauter@ens-lyon.org
  sylvain.chevillard@ens-lyon.org

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

  This program is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

*/

#include "signalhandling.h"

#define _POSIX_SOURCE
#define _POSIX_C_SOURCE 199309L

#include <sys/types.h>
#include <gmp.h>
#include <mpfr.h>
#include <stdio.h> /* fprintf, fopen, fclose, */
#include <stdlib.h> /* exit, free, mktemp */
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <setjmp.h>
#include <stdarg.h>
#include <limits.h>
#include "general.h"
#include "plot.h"
#include "expression.h"
#include "remez.h"
#include "infnorm.h"
#include "assignment.h"

#define YYPARSE_PARAM scanner

#include "parser.h"
#include <termios.h>
#include <sys/resource.h>
#include <unistd.h>
#include "execute.h"
#include "sollya-messaging.h"
#include "bitfields.h"
#include "printf.h"
#include "assignment.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if TIME_WITH_SYS_TIME
#include <sys/time.h>
#include <time.h>
#else
#if HAVE_SYS_TIME_H
#include <sys/time.h>
#else
#include <time.h>
#endif
#endif

/* A constant for the global reused MPFI and MPFR variables */
#define GLOBAL_REUSED_VARS_MAX_ALLOC 2048

/* STATE OF THE TOOL */

int oldAutoPrint = 0;
int oldVoidPrint = 0;
int oldExternalProcedurePrint = 0;
int printMode = PRINT_MODE_LEGACY;
FILE *warnFile = NULL;
char *variablename = NULL;
bitfield suppressedMessages = NULL;
mp_prec_t defaultprecision = DEFAULTPRECISION;
mp_prec_t tools_precision = DEFAULTPRECISION;
char *temporyDirectory = NULL;
char *uniqueIdentifier = NULL;
unsigned int globalSeed = 0;
int defaultpoints = DEFAULTPOINTS;
int taylorrecursions = DEFAULTTAYLORRECURSIONS;
int dyadic = 0;
int verbosity = 1;
int activateMessageNumbers = 0;
int canonical = 0;
int fileNumber = 0;
int autosimplify = 1;
int timecounting = 0;
chain *timeStack=NULL;
int fullParentheses=0;
int midpointMode = 0;
int dieOnErrorMode = 0;
int rationalMode = 0;
int noRoundingWarnings = 0;
int hopitalrecursions = DEFAULTHOPITALRECURSIONS;
mpfr_t statediam;

int eliminatePromptBackup;
int libraryMode = 0;
chain *readStack = NULL;
chain *readStackTemp = NULL;
chain *readStack2 = NULL;
void *scanner = NULL;
int promptToBePrinted = 0;
int lastWasSyntaxError = 0;
int lastCorrectlyExecuted = 0;
int helpNotFinished = 0;
int noColor = 0;

node *parsedThing = NULL;
node *parsedThingIntern = NULL;

jmp_buf recoverEnvironment;
int handlingCtrlC = 0;
int lastHandledSignal = 0;
int recoverEnvironmentReady = 0;
int exitInsteadOfRecover = 1;
int numberBacktrace = 1;
int displayColor = -1;
int oldMode = 0;

int (*messageCallback)(sollya_msg_t, void *) = NULL;
void *messageCallbackData = NULL;
int lastMessageCallbackResult = 1;
int lastMessageSuppressedResult = -1;

chain *symbolTable = NULL;
chain *declaredSymbolTable = NULL;
chain *backtraceStack = NULL;

int oldrlwrapcompatible = 0;

FILE *inputFile = NULL;
int inputFileOpened = 0;
int flushOutput = 0;

node *memRefChainStart = NULL;

/* END OF STATE OF THE TOOL */

/* HELPER VARIABLES FOR THE LEXER/ PARSER */

char *newReadFilename = NULL;
char *newReadFilenameTemp = NULL;
FILE **tempFDPtr;
FILE *temp_fd;
char *constBuffer = NULL;
char *constBuffer2 = NULL;
char *tempString = NULL;
char *tempString2 = NULL;
node *tempNode = NULL;
int tempInteger = 0;
int *tempIntPtr = NULL;

/* END OF HELPER VARIABLES */

/* HELPER VARIABLES FOR LEGACY MINIPARSER */

node * temp_node;
mpfr_t *mpfr_temp;
mpfr_t *mpfr_temp2;
node *minitree;
int miniparserCharactersRead;
int miniparserSemicolonAtEnd;
int miniparserEofReached;


/* END OF HELPER VARIABLES */

/* GLOBAL VARIABLES FOR THE MEMORY ALLOCATION FUNCTIONS */

void *wrapSafeRealloc(void *ptr, size_t old_size, size_t new_size);
void wrapSafeFree(void *ptr, size_t size);

void *(*actualCalloc)(size_t, size_t) = calloc;
void *(*actualMalloc)(size_t) = malloc;
void (*actualFree)(void *) = free;
void *(*actualRealloc)(void *, size_t) = realloc;
void (*actualFreeWithSize)(void *, size_t) = wrapSafeFree;
void *(*actualReallocWithSize)(void *, size_t, size_t) = wrapSafeRealloc;

void *(*oldGMPMalloc)(size_t) = NULL;
void *(*oldGMPRealloc)(void *, size_t, size_t) = NULL;
void (*oldGMPFree)(void *, size_t) = NULL;

void (*replacement_mp_set_memory_functions)(void *(*)(size_t),
					    void *(*)(void *, size_t, size_t),
					    void (*)(void *, size_t)) = NULL;
void (*replacement_mp_get_memory_functions)(void *(**)(size_t),
					    void *(**)(void *, size_t, size_t),
					    void (**)(void *, size_t)) = NULL;


/* END OF GLOBAL VARIABLES FOR THE MEMORY ALLOCATION FUNCTIONS */


/* HELPER VARIABLES THAT WE NEED TO MAINTAIN GLOBALLY USED FOR CERTAIN
   FUNCTIONS 
*/

int __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_vars_used = 0;
int __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_x_initialized = 0;
int __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_y_initialized = 0;
int __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_temp_initialized = 0;
sollya_mpfi_t __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_x;
sollya_mpfi_t __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_y;
mpfr_t __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_temp;
int    __sparsePolynomialEvalMpfr_var_used = 0;
int    __sparsePolynomialEvalMpfr_scratch_initialized = 0;
mpfr_t __sparsePolynomialEvalMpfr_scratch;
int    __sparsePolynomialEvalMpfi_var_used = 0;
int    __sparsePolynomialEvalMpfi_scratch_initialized = 0;
sollya_mpfi_t __sparsePolynomialEvalMpfi_scratch;

/* END OF HELPER VARIABLES THAT WE NEED TO MAINTAIN GLOBALLY USED FOR CERTAIN
   FUNCTIONS 
*/

/* Globally reused MPFI and MPFR variables */
sollya_mpfi_t *globalReusedMPFIVars = NULL;
unsigned int globalReusedMPFIVarsAllocated = 0;
unsigned int globalReusedMPFIVarsUsed = 0;
unsigned int globalReusedMPFIVarsInitialized = 0;
unsigned int globalReusedMPFIVarsMaxAllocated = GLOBAL_REUSED_VARS_MAX_ALLOC;

mpfr_t *globalReusedMPFRVars = NULL;
unsigned int globalReusedMPFRVarsAllocated = 0;
unsigned int globalReusedMPFRVarsUsed = 0;
unsigned int globalReusedMPFRVarsInitialized = 0;
unsigned int globalReusedMPFRVarsMaxAllocated = GLOBAL_REUSED_VARS_MAX_ALLOC;

/* End of globally reused MPFI and MPFR variables */

/* A global variable to check if sollyaLibPrintmessage has ever been
   called. 
*/

int sollyaLibPrintmessageCalled = 0;

/* End of global variable to check if sollyaLibPrintmessage has ever been
   called. 
*/

/* Variables for the prompt-passed arguments */
char **argsArgv = NULL;
int argsArgc = 0;
/* End of variables for the prompt-passed arguments */

/* Global variables for sollya_getc and sollya_fread */
int inside_sollya_getc = 0;
int inside_sollya_fread = 0;
int inside_sollya_feof = 0;
int inside_sollya_ferror = 0;
int inside_sollya_fwrite = 0;
/* End of global variables for sollya_getc */

extern int yyparse(void *); 
extern void yylex_destroy(void *);
extern int yylex_init(void **);
extern int yylex(void *);
extern void yyset_in(FILE *, void *);
extern int parserCheckEof();

void freeGlobalReusedMPFIVars();
void freeGlobalReusedMPFRVars();
void freeTool();

int sollya_gettime(sollya_time_t *t) {
#if defined(HAVE_CLOCK_GETTIME) && HAVE_CLOCK_GETTIME
  struct timespec ts;
#if !defined(__CYGWIN__)
#if defined(HAVE_GETTIMEOFDAY) && HAVE_GETTIMEOFDAY
  struct timeval tv;
#endif
#endif
  int res;
  res = clock_gettime(CLOCK_MONOTONIC, &ts);
  t->seconds = (int64_t) (ts.tv_sec);
  t->nano_seconds = (int64_t) (ts.tv_nsec);
  if (res == 0) return 1;
#if !defined(__CYGWIN__)
#if defined(HAVE_GETTIMEOFDAY) && HAVE_GETTIMEOFDAY
  res = gettimeofday(&tv,NULL);
  t->seconds = (int64_t) (tv.tv_sec);
  t->nano_seconds = ((int64_t) (tv.tv_usec)) * ((int64_t) 1000);
  if (res == 0) return 1;
#endif
#endif
#if defined(HAVE_TIME) && HAVE_TIME
  t->seconds = (int64_t) time(NULL);
  t->nano_seconds = (int64_t) 0;  
  return 1;
#else
  return 0;
#endif
#else
#if defined(HAVE_GETTIMEOFDAY) && HAVE_GETTIMEOFDAY
  struct timeval tv;
  int res;
  res = gettimeofday(&tv,NULL);
  t->seconds = (int64_t) (tv.tv_sec);
  t->nano_seconds = ((int64_t) (tv.tv_usec)) * ((int64_t) 1000);
  if (res == 0) return 1;
#if defined(HAVE_TIME) && HAVE_TIME
  t->seconds = (int64_t) time(NULL);
  t->nano_seconds = (int64_t) 0;  
  return 1;
#else
  return 0;
#endif  
#else
#if defined(HAVE_TIME) && HAVE_TIME
  t->seconds = (int64_t) time(NULL);
  t->nano_seconds = (int64_t) 0;
  return 1;
#else
  return 0;
#endif    
#endif
#endif
}

int64_t sollya_timediff_ms(sollya_time_t *before, sollya_time_t *after) {
  int64_t secdiff, nanosecdiff, millisecdiff;
  int64_t ten_to_nine, ten_to_six;
  ten_to_six = (int64_t) 1000;
  ten_to_six *= (int64_t) 1000;
  ten_to_nine = ten_to_six * ((int64_t) 1000);
  secdiff = after->seconds - before->seconds;
  nanosecdiff = after->nano_seconds - before->nano_seconds;
  while (nanosecdiff < ((int64_t) 0)) {
    nanosecdiff += ten_to_nine;
    secdiff--;
  }
  millisecdiff = (nanosecdiff / ten_to_six) + (((int64_t) 1000) * secdiff);
  return millisecdiff;
}

sollya_time_t *sollya_gettime_var() {
  return (sollya_time_t *) safeMalloc(sizeof(sollya_time_t));
}

size_t sollya_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *fd) {
  size_t res;

  deferSignalHandling();

  if (inside_sollya_fwrite) {
    sollyaFprintf(stderr,"Error: sollya_fwrite called from within sollya_fwrite. Bad things are happening.\n");
    exit(1);
  }
  inside_sollya_fwrite = 1;

  res = fwrite(ptr, size, nmemb, fd);

  inside_sollya_fwrite = 0;

  resumeSignalHandling();

  return res;
}

size_t sollya_fread(void *ptr, size_t size, size_t nmemb, FILE *fd) {
  size_t res;

  deferSignalHandling();

  if (inside_sollya_fread) {
    sollyaFprintf(stderr,"Error: sollya_fread called from within sollya_fread. Bad things are happening.\n");
    exit(1);
  }
  inside_sollya_fread = 1;

  res = fread(ptr, size, nmemb, fd);

  inside_sollya_fread = 0;

  resumeSignalHandling();

  return res;
}

int sollya_getc(FILE *fd) {
  int res;

  deferSignalHandling();

  if (inside_sollya_getc) {
    sollyaFprintf(stderr,"Error: sollya_getc called from within sollya_getc. Bad things are happening.\n");
    exit(1);
  }
  inside_sollya_getc = 1;

  res = getc(fd);

  inside_sollya_getc = 0;

  resumeSignalHandling();

  return res;
}

int sollya_feof(FILE *fd) {
  int res;

  deferSignalHandling();

  if (inside_sollya_feof) {
    sollyaFprintf(stderr,"Error: sollya_feof called from within sollya_feof. Bad things are happening.\n");
    exit(1);
  }
  inside_sollya_feof = 1;

  res = feof(fd);

  inside_sollya_feof = 0;

  resumeSignalHandling();

  return res;
}

int sollya_ferror(FILE *fd) {
  int res;

  deferSignalHandling();

  if (inside_sollya_ferror) {
    sollyaFprintf(stderr,"Error: sollya_ferror called from within sollya_ferror. Bad things are happening.\n");
    exit(1);
  }
  inside_sollya_ferror = 1;

  res = ferror(fd);

  inside_sollya_ferror = 0;

  resumeSignalHandling();

  return res;
}

void makeToolDie() {
	freeTool();

	if (!eliminatePromptBackup) sollyaPrintf("\n");

	if (inputFileOpened) {
	  fclose(inputFile);
	  inputFileOpened = 0;
	}

	exit(2);
}

void considerDyingOnError() {
  if (!dieOnErrorMode) return;

  printMessage(1,SOLLYA_MSG_TOOL_DIES_ON_ERROR_AS_PER_DIE_ON_ERROR_MODE,"Warning: some syntax, typing or side-effect error has occurred.\nAs the die-on-error mode is activated, the tool will be exited.\n");

  makeToolDie();
}

void normalMode() {
  if (displayColor == 0) return;
  displayColor = 0;
  if (noColor) return;
  if (eliminatePromptBackup) return;
  deferSignalHandling();
  fflush(NULL);
  printf("\e[0m");
  fflush(NULL);
  resumeSignalHandling();
}

void redMode() {
  if (displayColor == 1) return;
  displayColor = 1;
  if (noColor) return;
  if (eliminatePromptBackup) return;
  deferSignalHandling();
  fflush(NULL);
  printf("\e[0m\e[31m");
  fflush(NULL);
  resumeSignalHandling();
}

void blueMode() {
  if (displayColor == 2) return;
  displayColor = 2;
  if (eliminatePromptBackup) return;
  if (noColor) return;
  deferSignalHandling();
  fflush(NULL);
  printf("\e[0m\e[34m");
  fflush(NULL);
  resumeSignalHandling();
}

void setDisplayColor(int i) {
  switch (i) {
  case 2:
    blueMode();
    break;
  case 1:
    redMode();
    break;
  default:
    normalMode();
    break;
  }
}

int getDisplayColor() {
  return displayColor;
}

void saveMode() {
  oldMode = displayColor;
}

void changeToWarningMode() {
  saveMode();
  warningMode();
}

void restoreMode() {
  setDisplayColor(oldMode);
}


void blinkMode() {
  if (eliminatePromptBackup) return;
  deferSignalHandling();
  fflush(NULL);
  printf("\e[5m");
  fflush(NULL);
  resumeSignalHandling();
}

void unblinkMode() {
  if (eliminatePromptBackup) return;
  deferSignalHandling();
  fflush(NULL);
  printf("\e[25m");
  fflush(NULL);
  resumeSignalHandling();
}


void parseMode() {
  normalMode();
}

void outputMode() {
  blueMode();
}

void warningMode() {
  redMode();
}


void *safeCalloc (size_t nmemb, size_t size) {
  void *ptr;
  size_t myNmemb, mySize;

  myNmemb = nmemb;
  if (myNmemb == 0) myNmemb = 1;
  mySize = size;
  if (mySize == 0) mySize = 1;

  deferSignalHandling();
  ptr = actualCalloc(myNmemb,mySize);
  resumeSignalHandling();
  if (ptr == NULL) {
    sollyaFprintf(stderr,"Error: calloc could not succeed. No more memory left.\n");
    exit(1);
  }
  return ptr;
}

void *safeMalloc (size_t size) {
  void *ptr;
  size_t mySize;
  if (size == 0) mySize = 1; else mySize = size;
  deferSignalHandling();
  ptr = actualMalloc(mySize);
  resumeSignalHandling();  
  if (ptr == NULL) {
    sollyaFprintf(stderr,"Error: malloc could not succeed. No more memory left.\n");
    exit(1);
  }
  return ptr;
}

void *safeRealloc (void *ptr, size_t size) {
  void *newPtr;
  deferSignalHandling();
  newPtr = actualRealloc(ptr,size);
  resumeSignalHandling();
  if ((size != 0) && (newPtr == NULL)) {
    sollyaFprintf(stderr,"Error: realloc could not succeed. No more memory left.\n");
    exit(1);
  }
  return newPtr;
}

void safeFree(void *ptr) {
  deferSignalHandling();
  actualFree(ptr);
  resumeSignalHandling();
}

/* The gmp signature for realloc is strange, we have to wrap our function */
void *wrapSafeRealloc(void *ptr, size_t old_size, size_t new_size) {
  UNUSED_PARAM(old_size);
  return (void *) safeRealloc(ptr,new_size);
}

/* Same for GMP's free */
void wrapSafeFree(void *ptr, size_t size) {
  UNUSED_PARAM(size);
  safeFree(ptr);
}

/* Cling-wrap the GMP functions again */
void *clingWrapSafeReallocWithSize(void *ptr, size_t old_size, size_t new_size) {
  void *newPtr;

  newPtr = actualReallocWithSize(ptr, old_size, new_size);
  if ((new_size != 0) && (newPtr == NULL)) {
    sollyaFprintf(stderr,"Error: realloc could not succeed. No more memory left.\n");
    exit(1);
  }

  return newPtr;
}

void clingWrapSafeFreeWithSize(void *ptr, size_t size) {
  actualFreeWithSize(ptr, size);
}


/* Wrap the GMP mp_set_memory_functions and mp_get_memory_functions
   once in order to enable opportunities for different memory
   management.
*/
void sollya_mp_set_memory_functions(void *(*alloc_func_ptr) (size_t),
				    void *(*realloc_func_ptr) (void *, size_t, size_t),
				    void (*free_func_ptr) (void *, size_t)) {
  if (replacement_mp_set_memory_functions != NULL) {
    replacement_mp_set_memory_functions(alloc_func_ptr, realloc_func_ptr, free_func_ptr);
    return;
  }
  mp_set_memory_functions(alloc_func_ptr, realloc_func_ptr, free_func_ptr);
}

void sollya_mp_get_memory_functions(void *(**alloc_func_ptr) (size_t),
				    void *(**realloc_func_ptr) (void *, size_t, size_t),
				    void (**free_func_ptr) (void *, size_t)) {
  if (replacement_mp_get_memory_functions != NULL) {
    replacement_mp_get_memory_functions(alloc_func_ptr, realloc_func_ptr, free_func_ptr);
    return;
  }
  mp_get_memory_functions(alloc_func_ptr, realloc_func_ptr, free_func_ptr);
}

/* Wrap the sollya_mp_set_memory_functions function into one that not
   only installs the functions but also stores the old values into the
   global backup variables. Do nothing if the backup function pointers
   are not NULL (i.e. if they have already been set.)
*/
void wrap_mp_set_memory_functions(void *(*alloc_func_ptr) (size_t),
				  void *(*realloc_func_ptr) (void *, size_t, size_t),
				  void (*free_func_ptr) (void *, size_t)) {
  if ((oldGMPMalloc != NULL) ||
      (oldGMPRealloc != NULL) ||
      (oldGMPFree != NULL)) {
    return;
  }

  sollya_mp_get_memory_functions(&oldGMPMalloc,&oldGMPRealloc,&oldGMPFree);
  sollya_mp_set_memory_functions(alloc_func_ptr, realloc_func_ptr, free_func_ptr);
}

/* Provide memory management functions for the parsers

   These functions are wrappers that are currently provided to cope
   with special behavior of memory allocation in the future.
*/
void *parserCalloc(size_t nmemb, size_t size) {
  return safeCalloc(nmemb, size);
}

void *parserMalloc(size_t size) {
  return safeMalloc(size);
}

void *parserRealloc(void *ptr, size_t size) {
  return safeRealloc(ptr, size);
}

void parserFree(void *ptr) {
  safeFree(ptr);
}


/* Other legacy stuff */
char *maskString(char *src) {
  char *buf;
  char *res;
  char *curr1, *curr2;

  if (src == NULL) return NULL;

  buf = (char *) safeCalloc(strlen(src)*2+1,sizeof(char));

  curr2 = buf;
  for (curr1=src;*curr1!='\0';curr1++) {
    switch (*curr1) {
    case '"':
      *curr2 = '\\';
      curr2++;
      *curr2 = '"';
      curr2++;
      break;
    case '\n':
      *curr2 = '\\';
      curr2++;
      *curr2 = 'n';
      curr2++;
      break;
    case '\t':
      *curr2 = '\\';
      curr2++;
      *curr2 = 't';
      curr2++;
      break;
    case '\a':
      *curr2 = '\\';
      curr2++;
      *curr2 = 'a';
      curr2++;
      break;
    case '\b':
      *curr2 = '\\';
      curr2++;
      *curr2 = 'b';
      curr2++;
      break;
    case '\f':
      *curr2 = '\\';
      curr2++;
      *curr2 = 'f';
      curr2++;
      break;
    case '\r':
      *curr2 = '\\';
      curr2++;
      *curr2 = 'r';
      curr2++;
      break;
    case '\v':
      *curr2 = '\\';
      curr2++;
      *curr2 = 'v';
      curr2++;
      break;
    default:
      *curr2 = *curr1;
      curr2++;
      break;
    }
  }

  res = (char *) safeCalloc(strlen(buf)+1,sizeof(char));
  strcpy(res,buf);
  safeFree(buf);

  return res;
}

void demaskString(char *dest, char *src) {
  char *curr, *curr2;
  char internalBuf[4];
  int i;

  for (i=0;i<4;i++) internalBuf[i] = '\0';
  curr2 = dest;
  for (curr=src;*curr != '\0';curr++) {
    if (*curr != '\\') {
      *curr2++ = *curr;
    } else {
      curr++;
      if (*curr == '\0') break;
      switch (*curr) {
      case '\\':
	*curr2++ = '\\';
	break;
      case '"':
	*curr2++ = '"';
	break;
      case '?':
	*curr2++ = '?';
	break;
      case '\'':
	*curr2++ = '\'';
	break;
      case 'n':
	*curr2++ = '\n';
	break;
      case 't':
	*curr2++ = '\t';
	break;
      case 'a':
	*curr2++ = '\a';
	break;
      case 'b':
	*curr2++ = '\b';
	break;
      case 'f':
	*curr2++ = '\f';
	break;
      case 'r':
	*curr2++ = '\r';
	break;
      case 'v':
	*curr2++ = '\v';
	break;
      case 'x':
	curr++;
        i = 0;
	while ((i < 2) && (*curr != '\0') &&
	       (((*curr >= '0') && (*curr <= '9')) ||
		((*curr >= 'A') && (*curr <= 'F')) ||
		((*curr >= 'a') && (*curr <= 'f')))) {
	  internalBuf[i] = *curr;
	  curr++; i++;
	}
	curr--;
	for (i=0;i<2;i++) {
	  if ((internalBuf[i] >= 'a') && (internalBuf[i] <= 'f')) {
	    internalBuf[i] = (internalBuf[i] - 'a') + 'A';
	  }
	}
	*curr2++ = (char) sollya_strtol(internalBuf,NULL,16);
	break;
      default:
	i = 0;
	while ((i < 3) && (*curr != '\0') && (*curr >= '0') && (*curr <= '7')) {
	  internalBuf[i] = *curr;
	  curr++; i++;
	}
	curr--;
	*curr2++ = (char) sollya_strtol(internalBuf,NULL,8);
      }
    }
  }
}

int removeSpaces(char *outbuf, char *inbuf) {
  char *temp, *temp2;
  int removed;

  removed = 0;
  temp = inbuf; temp2 = outbuf;
  while ((temp != NULL) && (*temp != '\0')) {
    if (*temp != ' ') {
      *temp2 = *temp;
      temp2++;
    } else {
      removed = 1;
    }
    temp++;
  }

  return removed;
}

int removeMidpointMode(char *outbuf, char *inbuf) {
  char *temp, *temp2;
  int removed;

  removed = 0;

  temp = inbuf;
  temp2 = outbuf;

  while (*temp != '\0') {
    if (*temp != '~') {
      *temp2 = *temp;
      temp2++;
      temp++;
    } else {
      temp++;
      while ((*temp != '\0') && (*temp != '~')) temp++;
      if (*temp != '\0') temp++;
      removed = 1;
    }
  }

  return removed;
}

void cutMidpointStringIntoTwo(char *str1, char *str2, char *str) {
  char *temp1, *temp2, *temp;

  temp1 = str1;
  temp2 = str2;
  temp = str;
  while (*temp != '\0') {
    if (*temp != '~') {
      *temp1 = *temp;
      *temp2 = *temp;
      temp++;
      temp1++;
      temp2++;
    } else {
      temp++;
      while (*temp != '\0') {
	if (*temp != '/') {
	  *temp1 = *temp;
	  temp++;
	  temp1++;
	} else {
	  temp++;
	  while (*temp != '\0') {
	    if (*temp != '~') {
	      *temp2 = *temp;
	      temp++;
	      temp2++;
	    } else {
	      temp++;
	      break;
	    }
	  }
	  break;
	}
      }
    }
  }
}

void newReadFileStarted() {

}

void carriageReturnLexed() {
  if (helpNotFinished) { outputMode(); sollyaPrintf("This is %s. Having typed 'help', you have got to a special prompt.\nType now a semicolon (';') for an introduction on the %s help system.\nType now 'help;' for getting a list of available commands.\nType now a command name followed by a semicolon (';') for help on this command.\n>> ",PACKAGE_NAME,PACKAGE_NAME); }
  if (promptToBePrinted) printPrompt();
}

void newTokenLexed() {
  helpNotFinished = 0;
  promptToBePrinted = 0;
  lastCorrectlyExecuted = 0;
}

/* Precondition: fd can only be one of stdout and stderr */
int sollyaVfprintfSpecial(FILE *fd, const char *format, va_list varlist) {
  int res;

  if (printMode == PRINT_MODE_LEGACY) return sollyaInternalVfprintf(fd,format,varlist);

  if ((printMode == PRINT_MODE_WARNING_TO_STDERR) &&
      (displayColor == 1))
    return sollyaInternalVfprintf(stderr,format,varlist);

  if ((printMode == PRINT_MODE_WARNING_TO_FILE) &&
      (displayColor == 1) &&
      (warnFile != NULL)) {
    res = sollyaInternalVfprintf(warnFile,format,varlist);
    return res;
  }

  return sollyaInternalVfprintf(fd,format,varlist);
}

int sollyaVfprintf(FILE *fd, const char *format, va_list varlist) {

  if ((fd == stdout) || (fd == stderr))
    return sollyaVfprintfSpecial(fd,format,varlist);

  return sollyaInternalVfprintf(fd,format,varlist);
}

int installMessageCallback(int (*msgHandler) (sollya_msg_t, void *), void *data) {
  messageCallback = msgHandler;
  messageCallbackData = data;
  lastMessageCallbackResult = 1;
  return 1;
}

void getMessageCallback(int (**fptr)(sollya_msg_t, void *), void **data) {
  if (fptr != NULL) {
    *fptr = messageCallback;
  }
  if (data != NULL) {
    if (messageCallback != NULL) {
      *data = messageCallbackData;
    } else {
      *data = NULL;
    }
  }
}

int uninstallMessageCallback() {
  messageCallback = NULL;
  messageCallbackData = NULL;
  lastMessageCallbackResult = 1;
  return 1;
}

int messageHasEnoughVerbosityAndIsNotSuppressed(int verb, int msgNum) {
  int suppressed; 
  
  if ((verb >= 0) && (verbosity < verb)) return 0;
  /* Check if message suppression is activated for that message */
  suppressed = 0;
  if ((suppressedMessages != NULL) &&
      (verb >= 0) &&
      (msgNum != SOLLYA_MSG_NO_MSG)) {
    if (msgNum != SOLLYA_MSG_CONTINUATION) {
      suppressed = getBitInBitfield(suppressedMessages, msgNum);
    } else {
      if (lastMessageSuppressedResult == -1) {
	suppressed = 0;
      } else {
	suppressed = lastMessageSuppressedResult;
      }
    }
  }
  if ((verb >= 0) && suppressed && (msgNum != SOLLYA_MSG_NO_MSG)) return 0;

  return 1;
}

int printMessageInner(int verb, int msgNum, const char *format, va_list varlist) {
  int oldColor;
  int res, suppressed;
  const char *myFormat;
  const char *tempStr;
  sollya_msg_decl_t myMsg;

  if ((verb >= 0) && (verbosity < verb)) return 0;

  /* Check if message suppression is activated for that message */
  suppressed = 0;
  if ((suppressedMessages != NULL) &&
      (verb >= 0) &&
      (msgNum != SOLLYA_MSG_NO_MSG)) {
    if (msgNum != SOLLYA_MSG_CONTINUATION) {
      suppressed = getBitInBitfield(suppressedMessages, msgNum);
    } else {
      if (lastMessageSuppressedResult == -1) {
	suppressed = 0;
      } else {
	suppressed = lastMessageSuppressedResult;
      }
    }
  }
  lastMessageSuppressedResult = suppressed;
  if ((verb >= 0) && suppressed && (msgNum != SOLLYA_MSG_NO_MSG)) return 0;

  /* If there is a message callback installed, call it.
     If it says no message is to be displayed, just bail out.

     Do not call the message callback handler for no messages and
     continuation messages. In the case of a continuation message,
     take the last "display/don't display" value instead of what the
     handler would return.
  */
  if ((msgNum != SOLLYA_MSG_CONTINUATION) || (messageCallback == NULL)) lastMessageCallbackResult = 1;
  if ((msgNum != SOLLYA_MSG_NO_MSG) &&
      (msgNum != SOLLYA_MSG_CONTINUATION) &&
      (messageCallback != NULL)) {
    myMsg.msg_id = msgNum;
    lastMessageCallbackResult = messageCallback(&myMsg, messageCallbackData);
  }
  if (!lastMessageCallbackResult) return 0;

  oldColor = displayColor;

  if ((verb >= 1) || (verb < 0)) warningMode(); else outputMode();

  if (activateMessageNumbers && (msgNum != SOLLYA_MSG_CONTINUATION) && (msgNum != SOLLYA_MSG_NO_MSG)) {
    myFormat = format;
    res = 0;
    if (((tempStr = sollya_strstr(format,"Warning")) != NULL) && (tempStr == format)) {
      if (verb >= 0) {
	res += sollyaFprintf(stdout,"Warning (%d)",msgNum);
      } else {
	res += sollyaFprintf(stderr,"Warning (%d)",msgNum);
      }
      myFormat = tempStr + strlen("Warning");
    } else {
      if (((tempStr = sollya_strstr(format,"Error")) != NULL) && (tempStr == format)) {
	if (verb >= 0) {
	  res += sollyaFprintf(stdout,"Error (%d)",msgNum);
	} else {
	  res += sollyaFprintf(stderr,"Error (%d)",msgNum);
	}
	myFormat = tempStr + strlen("Error");
      } else {
	if (((tempStr = sollya_strstr(format,"Information")) != NULL) && (tempStr == format)) {
	  if (verb >= 0) {
	    res += sollyaFprintf(stdout,"Information (%d)",msgNum);
	  } else {
	    res += sollyaFprintf(stderr,"Information (%d)",msgNum);
	  }
	  myFormat = tempStr + strlen("Information");
	} else {
	  if (verb >= 0) {
	    res += sollyaFprintf(stdout,"Message (%d): ",msgNum);
	  } else {
	    res += sollyaFprintf(stderr,"Message (%d): ",msgNum);
	  }
	}
      }
    }
    if (verb >= 0) {
      res += sollyaVfprintf(stdout,myFormat,varlist);
    } else {
      res += sollyaVfprintf(stderr,myFormat,varlist);
    }
  } else {
    if (verb >= 0) {
      res = sollyaVfprintf(stdout,format,varlist);
    } else {
      res = sollyaVfprintf(stderr,format,varlist);
    }
  }

  setDisplayColor(oldColor);

  return res;
}

int printMessage(int verb, int msgNum, const char *format, ...) {
  int res;
  va_list varlist;

  va_start(varlist,format);
  res = printMessageInner(verb, msgNum, format, varlist);
  va_end(varlist);

  return res;
}

int sollyaLibPrintmessage(int verb, int cont, const char *format, ...) {
  va_list varlist;
  int res;

  va_start(varlist,format);
  if (cont) {
    if (sollyaLibPrintmessageCalled) {
      sollyaLibPrintmessageCalled = 1;
      res = printMessageInner(verb, SOLLYA_MSG_CONTINUATION, format, varlist);
    } else {
      sollyaLibPrintmessageCalled = 1;
      res = printMessageInner(verb, SOLLYA_MSG_GENERIC_SOLLYA_LIBRARY_MSG, format, varlist);
    }
  } else {
    sollyaLibPrintmessageCalled = 1;
    res = printMessageInner(verb, SOLLYA_MSG_GENERIC_SOLLYA_LIBRARY_MSG, format, varlist);
  }
  va_end(varlist);
  
  return res;
}

int sollyaPrintf(const char *format, ...) {
  va_list varlist;
  int res;

  va_start(varlist,format);

  res = sollyaVfprintf(stdout,format,varlist);

  va_end(varlist);

  return res;
}

int sollyaFprintf(FILE *fd, const char *format, ...) {
  va_list varlist;
  int res;

  va_start(varlist,format);

  res = sollyaVfprintf(fd,format,varlist);

  va_end(varlist);

  return res;
}

/* Returns a string representing x written in binary, such
   that it can be read back by mpfr_set_str.
   x must be a valid number. If x is an Inf or NaN, NULL is
   returned.
   Note: the string must be freed afterwards, of course. */
char *mpfr_to_binary_str(mpfr_t x) {
  char *s;
  char *ptr;
  mp_exp_t e;
  mp_prec_t prec = mpfr_get_prec(x);

  if (!mpfr_number_p(x)) return NULL;

  if (mpfr_zero_p(x)) {
    s = calloc(5, sizeof(char));
    s[0] = '0';
    s[1] = 'p';
    s[2] = '+';
    s[3] = '0';
    return s;
  }

  /* The documentation of mpfr_get_str explains that it needs
     max ( mpfr_get_prec(x)+2, 7 ). This is achieved by
     mpfr_get_prec(x)+7. We must add 4 chars for "0." and "p+".
     Finally, 20 extra bits are sufficient to write a 64-bit
     exponent in base 10 */
  s = (char *)calloc(prec+7+4+20, sizeof(char));
  mpfr_get_str(s+2, &e, 2, 0, x, GMP_RNDN);
  if ( s[2] == '-' ) {
    s[0] = '-';
    ptr = s+1;
  }
  else ptr = s;

  ptr[0] = '0';
  ptr[1] = '.';
  ptr[prec+2] = 'p';
  if (e>=0) {
    ptr[prec+3] = '+';
    ptr = ptr + prec + 4;
  }
  else ptr = ptr + prec + 3;
  sollya_snprintf(ptr, prec+7+4+20, "%ld", (long)e);

  return s;
}

void freeCounter(void) {
  freeChain(timeStack, safeFree);
  timeStack=NULL;
  return;
}

void pushTimeCounter(void) {
  sollya_time_t *buf;
  if(timecounting==1) {
    buf = safeMalloc(sizeof(sollya_time_t));
    if(!sollya_gettime(buf))
      sollyaFprintf(stderr, "Error: unable to use the timer. Measures may be untrustable\n");
    timeStack = addElement(timeStack, buf);
  }
  return;
}

void popTimeCounter(char *s) {
  sollya_time_t *buf_init;
  sollya_time_t *buf_final;

  long long int days, hours, minutes, seconds, milliseconds, nanoseconds;

  chain *prev;
  if((timecounting==1)&&(timeStack!=NULL)) {
    buf_final = safeMalloc(sizeof(sollya_time_t));
    if(!sollya_gettime(buf_final))
      sollyaFprintf(stderr, "Error: unable to use the timer. Measures may be untrustable\n");
    buf_init = timeStack->value;

    seconds = (long long int)(buf_final->seconds) - (long long int)(buf_init->seconds);
    nanoseconds = (long long int)(buf_final->nano_seconds) - (long long int)(buf_init->nano_seconds);

    if (nanoseconds < 0) {
      nanoseconds += 1000000000ll;
      seconds--;
    }

    milliseconds = nanoseconds / 1000000ll;

    if((milliseconds>0)||(seconds>0)) {
      changeToWarningMode();
      sollyaPrintf("Information: %s spent ", s);
      if(seconds!=0) {
	minutes = seconds / 60;
	seconds = seconds % 60;
	hours = minutes / 60;
	minutes = minutes % 60;
	days = hours / 24;
	hours = hours % 24;

	if(days!=0) sollyaPrintf("%ld days, ", days);
	if(hours!=0) sollyaPrintf("%ld hours, ", hours);
	if(minutes!=0) sollyaPrintf("%ld minutes, ", minutes);
	if(seconds!=0) sollyaPrintf("%ld seconds, ", seconds);
      }
      sollyaPrintf("%ld ms\n", milliseconds);
      restoreMode();
    }

    prev = timeStack;
    timeStack = timeStack->next;
    safeFree(prev);
    safeFree(buf_init);
    safeFree(buf_final);
  }
  return;
}


void printPrompt(void) {
  if (eliminatePromptBackup) return;
  if (readStack != NULL) return;
  parseMode();
  if (oldrlwrapcompatible && (!noColor)) {
    deferSignalHandling();
    printf("\e[1A\n");
    fflush(NULL);
    resumeSignalHandling();
  }
  sollyaPrintf("> ");
}

void freeFunctionSpecialVariables() {
 __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_vars_used = 0;
 if (__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_x_initialized) {
   sollya_mpfi_clear(__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_x);
   __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_x_initialized = 0;
 }
 if (__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_y_initialized) {
   sollya_mpfi_clear(__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_y);
   __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_y_initialized = 0;
 }
 if (__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_temp_initialized) {
   mpfr_clear(__firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_temp);
   __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_temp_initialized = 0;
 }
 __sparsePolynomialEvalMpfr_var_used = 0;
 if (__sparsePolynomialEvalMpfr_scratch_initialized) {
   mpfr_clear(__sparsePolynomialEvalMpfr_scratch);
   __sparsePolynomialEvalMpfr_scratch_initialized = 0;
 }
 __sparsePolynomialEvalMpfi_var_used = 0;
 if (__sparsePolynomialEvalMpfi_scratch_initialized) {
   sollya_mpfi_clear(__sparsePolynomialEvalMpfi_scratch);
   __sparsePolynomialEvalMpfi_scratch_initialized = 0;
 }
}

void freeTool() {
  if(variablename != NULL) safeFree(variablename);
  if(newReadFilename != NULL) safeFree(newReadFilename);
  if(suppressedMessages != NULL) freeBitfield(suppressedMessages);

  if (!(eliminatePromptBackup == 1)) {
    removePlotFiles();
  }

  while ((readStack != NULL) && (readStack2 != NULL)) {
    temp_fd = *((FILE **) (readStack2->value));
    fclose(temp_fd);
    safeFree(readStack2->value);
    readStackTemp = readStack2->next;
    safeFree(readStack2);
    readStack2 = readStackTemp;
    safeFree(readStack->value);
    readStackTemp = readStack->next;
    safeFree(readStack);
    readStack = readStackTemp;
  }
  yylex_destroy(scanner);
  freeFunctionLibraries();
  freeConstantLibraries();
  freeProcLibraries();
  freeCounter();
  freeSymbolTable(symbolTable, freeThingOnVoid);
  symbolTable = NULL;
  freeDeclaredSymbolTable(declaredSymbolTable, freeThingOnVoid);
  freeBacktraceStack();
  freeFunctionSpecialVariables();
  freeGlobalReusedMPFIVars();
  freeGlobalReusedMPFRVars();
  polynomialFreeCaches();
  declaredSymbolTable = NULL;
  backtraceStack = NULL;
  mpfr_clear(statediam);
  safeFree(temporyDirectory); temporyDirectory = NULL;
  safeFree(uniqueIdentifier); uniqueIdentifier = NULL;
  mpfr_free_cache();
  normalMode();
}

void initToolDefaults() {
  int k;
  node *tempNode;
  chain *tempList;
  char *pidStr, *uniqueStr;
  unsigned int mySeed;
  char *c;

  if(variablename != NULL) safeFree(variablename);
  variablename = NULL;
  if(suppressedMessages != NULL) freeBitfield(suppressedMessages);
  suppressedMessages = NULL;
  defaultprecision = DEFAULTPRECISION;
  tools_precision = DEFAULTPRECISION;
  defaultpoints = DEFAULTPOINTS;
  taylorrecursions = DEFAULTTAYLORRECURSIONS;
  memRefChainStart = NULL;
  dyadic = 0;
  verbosity = 1;
  activateMessageNumbers = 0;
  canonical = 0;
  fileNumber = 0;
  autosimplify = 1;
  rationalMode = 0;
  timecounting = 0;
  if (timeStack != NULL) freeCounter();
  timeStack=NULL;
  fullParentheses=0;
  midpointMode = 0;
  noRoundingWarnings = eliminatePromptBackup;
  hopitalrecursions = DEFAULTHOPITALRECURSIONS;
  symbolTable = NULL;
  declaredSymbolTable = NULL;
  backtraceStack = NULL;
  mpfr_init2(statediam,10);
  mpfr_set_d(statediam,DEFAULTDIAM,GMP_RNDN);
  __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_vars_used = 0;
  __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_x_initialized = 0;
  __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_y_initialized = 0;
  __firstTryEvaluateFaithfulWithCutOffFastInternalImplementation_temp_initialized = 0;
  __sparsePolynomialEvalMpfr_var_used = 0;
  __sparsePolynomialEvalMpfr_scratch_initialized = 0;
  __sparsePolynomialEvalMpfi_var_used = 0;
  __sparsePolynomialEvalMpfi_scratch_initialized = 0;
  globalReusedMPFIVars = NULL;
  globalReusedMPFIVarsAllocated = 0;
  globalReusedMPFIVarsUsed = 0;
  globalReusedMPFIVarsInitialized = 0;
  globalReusedMPFIVarsMaxAllocated = GLOBAL_REUSED_VARS_MAX_ALLOC;
  globalReusedMPFRVars = NULL;
  globalReusedMPFRVarsAllocated = 0;
  globalReusedMPFRVarsUsed = 0;
  globalReusedMPFRVarsInitialized = 0;
  globalReusedMPFRVarsMaxAllocated = GLOBAL_REUSED_VARS_MAX_ALLOC;
  polynomialInitializeCaches();
  sollyaLibPrintmessageCalled = 0;

  if ((argsArgc < 1) || (argsArgv == NULL)) {
    tempNode = makeEmptyList();
  } else {
    tempList = NULL;
    for (k=argsArgc-1;k>=0;k--) {
      tempList = addElement(tempList, makeString(argsArgv[k]));      
    }
    tempNode = makeList(tempList);
  }
  tempNode = addMemRef(tempNode);
  symbolTable = addEntry(symbolTable, "__argv", tempNode, copyThingOnVoid);
  freeThing(tempNode);
  pidStr = getUniqueId();
  uniqueStr = (char *) safeCalloc(4 + strlen(PACKAGE_STRING) + 1 + strlen(pidStr) + 1 + 8 * sizeof(int) + 1, sizeof(char));
  mySeed = (unsigned int) ((unsigned int) time(NULL)) + ((unsigned int) globalSeed);
  globalSeed = rand_r(&mySeed);
  sollya_snprintf(uniqueStr, 4 + strlen(PACKAGE_STRING) + 1 + strlen(pidStr) + 1 + 8 * sizeof(int) + 1, "_id_%s_%s_%08d", PACKAGE_STRING, pidStr, rand_r(&mySeed));
  for (c=uniqueStr;*c!='\0';c++) {
    if ((*c == ' ') || 
	(*c == '\t') || 
	(*c == '\n') ||
	(*c == '.') ||
	(*c == '/')) {
      *c = '_';
    }
  }
  for (c=uniqueStr;*c!='\0';c++) {
    if (!((('A' <= *c) && (*c <= 'Z')) ||
	  (('a' <= *c) && (*c <= 'z')) ||
	  (('0' <= *c) && (*c <= '9')) ||
	  (*c == '_'))) {
      *c = '_';
    }
  }
  tempNode = makeString(uniqueStr);
  symbolTable = addEntry(symbolTable, "__unique_id", tempNode, copyThingOnVoid);
  freeThing(tempNode);
  safeFree(uniqueStr);
  parseMode();
}

void restartTool() {
  mpfr_clear(statediam);
  freeSymbolTable(symbolTable, freeThingOnVoid);
  symbolTable = NULL;
  freeDeclaredSymbolTable(declaredSymbolTable, freeThingOnVoid);
  freeBacktraceStack();
  freeFunctionSpecialVariables();
  freeGlobalReusedMPFIVars();
  freeGlobalReusedMPFRVars();
  polynomialFreeCaches();
  declaredSymbolTable = NULL;
  backtraceStack = NULL;
  freeFunctionLibraries();
  freeConstantLibraries();
  freeProcLibraries();
  mpfr_free_cache();
  initToolDefaults();
  parseMode();
}

char *getNameOfVariable() {
  char *str;
  if (variablename == NULL) return NULL;
  str = (char *) safeCalloc(strlen(variablename)+1,sizeof(char));
  strcpy(str,variablename);
  return str;
}

int setNameOfVariable(char *str) {
  if (variablename != NULL) return 0;
  variablename = (char *) safeCalloc(strlen(str)+1,sizeof(char));
  strcpy(variablename,str);
  return 1;
}

char *initTempDir() {
  char *res;
  char *staticRes;
  size_t len;

  staticRes = getenv("TMPDIR");
  if ((staticRes != NULL) && ((len = strlen(staticRes)) != 0)) {
    res = safeCalloc(len + 1, sizeof(char));
    strcpy(res, staticRes);
  } else {
    staticRes = getenv("%TEMP%");
    if ((staticRes != NULL) && ((len = strlen(staticRes)) != 0)) {
      res = safeCalloc(len + 1, sizeof(char));
      strcpy(res, staticRes);
    } else {
      staticRes = getenv("%TMP%");
      if ((staticRes != NULL) && ((len = strlen(staticRes)) != 0)) {
	res = safeCalloc(len + 1, sizeof(char));
	strcpy(res, staticRes);
      } else {
	staticRes = "/tmp";
	len = strlen(staticRes);
	res = safeCalloc(len + 1, sizeof(char));
	strcpy(res, staticRes);
      }
    }
  }

  return res;
}

char *initUniqueId() {
  char *res, *str;
  pid_t pid;
  size_t size;

  size = 8 * sizeof(signed long long int) + 10;

  pid = getpid();

  str = safeCalloc(size + 1, sizeof(char));
  sollya_snprintf(str, size, "%lld", (signed long long int) pid);

  res = safeCalloc(strlen(str) + 1, sizeof(char));
  strcpy(res, str);
  safeFree(str);

  return res;
}

char *getTempDir() {

  if (temporyDirectory == NULL) {
    temporyDirectory = initTempDir();
  }

  return temporyDirectory;
}

char *getUniqueId() {

  if (uniqueIdentifier == NULL) {
    uniqueIdentifier = initUniqueId();
  }

  return uniqueIdentifier;
}

mp_prec_t getToolPrecision() {
  return tools_precision;
}

void setToolPrecision(mp_prec_t prec) {
  tools_precision = prec;
  defaultprecision = prec;
}

int getToolPoints() {
  return defaultpoints;
}

void setToolPoints(int p) {
  if (p < 12) defaultpoints = 12; else defaultpoints = p;
}

int getToolTaylorRecursions() {
  return taylorrecursions;
}

void setToolTaylorRecursions(int i) {
  if (i < 0) taylorrecursions = 0; else taylorrecursions = i;
}

int getToolHopitalRecursions() {
  return hopitalrecursions;
}

void setToolHopitalRecursions(int i) {
  if (i < 0) hopitalrecursions = 0; else hopitalrecursions = i;
}


int getToolDiameter(mpfr_t rop) {
  return mpfr_set(rop,statediam,GMP_RNDN);
}

void setToolDiameter(mpfr_t op) {
  mpfr_set_prec(statediam,mpfr_get_prec(op));
  mpfr_set(statediam,op,GMP_RNDN);
}

int getDisplayMode() {
  return dyadic;
}

int setDisplayMode(int newMode) {
  if ((0 <= newMode) && (newMode <= 4)) {
    dyadic = newMode;
    return 1;
  } else {
    dyadic = 0;
    return 0;
  }
}

int getVerbosity() {
  return verbosity;
}

int setVerbosity(int newVerbosity) {
  if (newVerbosity >= 0) {
    verbosity = newVerbosity;
    return 1;
  } else {
    verbosity = 0;
    return 0;
  }
}

int getCanonical() {
  return canonical;
}

void setCanonical(int newCanonical) {
  canonical = (!(!newCanonical));
}

int getAutosimplify() {
  return autosimplify;
}

void setAutosimplify(int newAutosimplify) {
  autosimplify = (!(!newAutosimplify));
}

int getFullParentheses() {
  return fullParentheses;
}

void setFullParentheses(int newFullParentheses) {
  fullParentheses = (!(!newFullParentheses));
}

int getMidpointMode() {
  return midpointMode;
}

void setMidpointMode(int newMidpointMode) {
  midpointMode = (!(!newMidpointMode));
}

int getDieOnErrorMode() {
  return dieOnErrorMode;
}

void setDieOnErrorMode(int newDieOnErrorMode) {
  dieOnErrorMode = (!(!newDieOnErrorMode));
}

int getTimecounting() {
  return timecounting;
}

void setTimecounting(int newTimecounting) {
  timecounting = (!(!newTimecounting));
}

int getRoundingWarnings() {
  return (!noRoundingWarnings);
}

void setRoundingWarnings(int newRoundingWarnings) {
  noRoundingWarnings = (!newRoundingWarnings);
}

int getRationalMode() {
  return rationalMode;
}

void setRationalMode(int newRationalMode) {
  rationalMode = (!(!newRationalMode));
}

void setRecoverEnvironment(jmp_buf *env) {
  sollya_memmove(&recoverEnvironment,env,sizeof(recoverEnvironment));
  exitInsteadOfRecover = 0;
}

void invalidateRecoverEnvironment() {
  exitInsteadOfRecover = 1;
}

int initializeLibraryMode(void *(*myActualMalloc)(size_t),
			  void *(*myActualCalloc)(size_t, size_t),
			  void *(*myActualRealloc)(void *, size_t),
			  void (*myActualFree)(void*),
			  void *(*myActualReallocWithSize)(void *, size_t, size_t),
			  void (*myActualFreeWithSize)(void *, size_t),
			  int argc,
			  char **argv,
			  void (*my_mp_set_func)(void *(*)(size_t),
						 void *(*)(void *, size_t, size_t),
						 void (*)(void *, size_t)),
			  void (*my_mp_get_func)(void *(**)(size_t),
						 void *(**)(void *, size_t, size_t),
						 void (**)(void *, size_t))) {
  void *ptr;
  int k, allArgsOkay;
  libraryMode = 1;
  oldGMPMalloc = NULL;
  oldGMPRealloc = NULL;
  oldGMPFree = NULL;
  replacement_mp_set_memory_functions = NULL;
  replacement_mp_get_memory_functions = NULL;
  actualCalloc = calloc;
  actualMalloc = malloc;
  actualFree = free;
  actualRealloc = realloc;
  actualFreeWithSize = wrapSafeFree;
  actualReallocWithSize = wrapSafeRealloc;  
  if (my_mp_set_func != NULL) replacement_mp_set_memory_functions = my_mp_set_func;
  if (my_mp_get_func != NULL) replacement_mp_get_memory_functions = my_mp_get_func;
  if (myActualMalloc != NULL) actualMalloc = myActualMalloc;
  if (myActualCalloc != NULL) actualCalloc = myActualCalloc;
  if (myActualRealloc != NULL) actualRealloc = myActualRealloc;
  if (myActualFree != NULL) actualFree = myActualFree;
  if (myActualReallocWithSize != NULL) actualReallocWithSize = myActualReallocWithSize;
  if (myActualFreeWithSize != NULL) actualFreeWithSize = myActualFreeWithSize;
  messageCallback = NULL;
  lastMessageCallbackResult = 1;
  lastMessageSuppressedResult = -1;
  inputFileOpened = 0;
  flushOutput = 0;
  oldAutoPrint = 0;
  oldExternalProcedurePrint = 0;
  printMode = PRINT_MODE_LEGACY;
  warnFile = NULL;
  eliminatePromptBackup = 1;
  wrap_mp_set_memory_functions(safeMalloc,clingWrapSafeReallocWithSize,clingWrapSafeFreeWithSize);
  if ((argc >= 1) && (argv != NULL)) {
    for (k=0,allArgsOkay = 1;(k<argc)&&allArgsOkay;k++) {
      allArgsOkay = allArgsOkay && (argv[k] != NULL);
    }
    if (allArgsOkay) {
      argsArgc = argc;
      argsArgv = (char **) safeCalloc(argsArgc, sizeof(char *));
      for (k=0;k<argsArgc;k++) {
	argsArgv[k] = (char *) safeCalloc(strlen(argv[k])+1,sizeof(char));
	strcpy(argsArgv[k], argv[k]);
      }
    }
  }
  initToolDefaults();
  handlingCtrlC = 0;
  lastHandledSignal = 0;
  noRoundingWarnings = 0;
  ptr = safeMalloc(sizeof(char));
  safeFree(ptr);
  ptr = safeCalloc(1,sizeof(char));
  safeFree(ptr);
  ptr = safeMalloc(sizeof(char));
  ptr = safeRealloc(ptr, 2 * sizeof(char));
  safeFree(ptr);
  return 1;
}

int finalizeLibraryMode() {
  int k;
  if (argsArgv != NULL) {
    for (k=0;k<argsArgc;k++) {
      safeFree(argsArgv[k]);
    }
    safeFree(argsArgv);
    argsArgv = NULL;
  }
  if(variablename != NULL) safeFree(variablename);
  variablename = NULL;
  if(newReadFilename != NULL) safeFree(newReadFilename);
  newReadFilename = NULL;
  if(suppressedMessages != NULL) freeBitfield(suppressedMessages);
  suppressedMessages = NULL;

  if (!(eliminatePromptBackup == 1)) {
    removePlotFiles();
  }

  while ((readStack != NULL) && (readStack2 != NULL)) {
    temp_fd = *((FILE **) (readStack2->value));
    fclose(temp_fd);
    safeFree(readStack2->value);
    readStackTemp = readStack2->next;
    safeFree(readStack2);
    readStack2 = readStackTemp;
    safeFree(readStack->value);
    readStackTemp = readStack->next;
    safeFree(readStack);
    readStack = readStackTemp;
  }
  freeFunctionLibraries();
  freeConstantLibraries();
  freeProcLibraries();
  freeCounter();
  freeSymbolTable(symbolTable, freeThingOnVoid);
  symbolTable = NULL;
  freeDeclaredSymbolTable(declaredSymbolTable, freeThingOnVoid);
  freeBacktraceStack();
  freeFunctionSpecialVariables();
  freeGlobalReusedMPFIVars();
  freeGlobalReusedMPFRVars();
  polynomialFreeCaches();
  declaredSymbolTable = NULL;
  backtraceStack = NULL;
  mpfr_clear(statediam);
  safeFree(temporyDirectory); temporyDirectory = NULL;
  safeFree(uniqueIdentifier); uniqueIdentifier = NULL;
  mpfr_free_cache();
  uninstallMessageCallback();
  sollya_mp_set_memory_functions(oldGMPMalloc,oldGMPRealloc,oldGMPFree);
  actualCalloc = calloc;
  actualMalloc = malloc;
  actualFree = free;
  actualRealloc = realloc;
  actualFreeWithSize = wrapSafeFree;
  actualReallocWithSize = wrapSafeRealloc;
  replacement_mp_set_memory_functions = NULL;
  replacement_mp_get_memory_functions = NULL;
  oldGMPMalloc = NULL;
  oldGMPRealloc = NULL;
  oldGMPFree = NULL;
  libraryMode = 0;
  return 1;
}

double sollya_mpfr_get_d(mpfr_t op, mpfr_rnd_t rnd) {
  union {
    uint64_t i;
    double d;
  } xdb;

  /* Handle NaNs on our own */
  if (mpfr_nan_p(op)) {
    xdb.i = 0x7ff800007ff80000ull;
    return xdb.d;
  }

  return mpfr_get_d(op, rnd);
}

void allocateReusedGlobalMPFIVars() {
  if (globalReusedMPFIVars != NULL) return;
  if (globalReusedMPFIVarsAllocated != 0) return;
  globalReusedMPFIVars = (sollya_mpfi_t *) safeCalloc(globalReusedMPFIVarsMaxAllocated, sizeof(sollya_mpfi_t));
  globalReusedMPFIVarsAllocated = globalReusedMPFIVarsMaxAllocated;
  globalReusedMPFIVarsUsed = 0;
  globalReusedMPFIVarsInitialized = 0;
}

void allocateReusedGlobalMPFRVars() {
  if (globalReusedMPFRVars != NULL) return;
  if (globalReusedMPFRVarsAllocated != 0) return;
  globalReusedMPFRVars = (mpfr_t *) safeCalloc(globalReusedMPFRVarsMaxAllocated, sizeof(mpfr_t));
  globalReusedMPFRVarsAllocated = globalReusedMPFRVarsMaxAllocated;
  globalReusedMPFRVarsUsed = 0;
  globalReusedMPFRVarsInitialized = 0;
}

sollya_mpfi_t *getReusedGlobalMPFIVars(unsigned int n, mp_prec_t prec) {
  sollya_mpfi_t *ptr;
  unsigned int i;

  if (n == 0) return NULL;
  if ((globalReusedMPFIVars == NULL) || 
      (globalReusedMPFIVarsAllocated == 0)) allocateReusedGlobalMPFIVars();
  if ((globalReusedMPFIVarsAllocated - globalReusedMPFIVarsUsed) < n) return NULL;
  ptr = &(globalReusedMPFIVars[globalReusedMPFIVarsUsed]);
  globalReusedMPFIVarsUsed += n;
  for (i=globalReusedMPFIVarsInitialized;i<globalReusedMPFIVarsUsed;i++) {
    sollya_mpfi_init2(globalReusedMPFIVars[i],prec);
  }
  if (globalReusedMPFIVarsUsed > globalReusedMPFIVarsInitialized) {
    globalReusedMPFIVarsInitialized = globalReusedMPFIVarsUsed;
  }
  for (i=0;i<n;i++) {
    sollya_mpfi_set_prec(ptr[i],prec);
  }

  return ptr;
}

mpfr_t *getReusedGlobalMPFRVars(unsigned int n, mp_prec_t prec) {
  mpfr_t *ptr;
  unsigned int i;

  if (n == 0) return NULL;
  if ((globalReusedMPFRVars == NULL) || 
      (globalReusedMPFRVarsAllocated == 0)) allocateReusedGlobalMPFRVars();
  if ((globalReusedMPFRVarsAllocated - globalReusedMPFRVarsUsed) < n) return NULL;
  ptr = &(globalReusedMPFRVars[globalReusedMPFRVarsUsed]);
  globalReusedMPFRVarsUsed += n;
  for (i=globalReusedMPFRVarsInitialized;i<globalReusedMPFRVarsUsed;i++) {
    mpfr_init2(globalReusedMPFRVars[i],prec);
  }
  if (globalReusedMPFRVarsUsed > globalReusedMPFRVarsInitialized) {
    globalReusedMPFRVarsInitialized = globalReusedMPFRVarsUsed;
  }
  for (i=0;i<n;i++) {
    mpfr_set_prec(ptr[i],prec);
  }

  return ptr;
}

void returnReusedGlobalMPIVars(unsigned int n) {
  if (n == 0) return;
  if (n > globalReusedMPFIVarsUsed) {
    globalReusedMPFIVarsUsed = 0;
  } else {
    globalReusedMPFIVarsUsed -= n;
  }
}

void returnReusedGlobalMPFRVars(unsigned int n) {
  if (n == 0) return;
  if (n > globalReusedMPFRVarsUsed) {
    globalReusedMPFRVarsUsed = 0;
  } else {
    globalReusedMPFRVarsUsed -= n;
  }
}

void freeGlobalReusedMPFIVars() {
  unsigned int i;

  if (globalReusedMPFIVars == NULL) return;
  if (globalReusedMPFIVarsAllocated == 0) return;

  for (i=0;i<globalReusedMPFIVarsInitialized;i++) {
    sollya_mpfi_clear(globalReusedMPFIVars[i]);
  }
  safeFree(globalReusedMPFIVars);
  globalReusedMPFIVars = NULL;
  globalReusedMPFIVarsAllocated = 0;
  globalReusedMPFIVarsUsed = 0;
  globalReusedMPFIVarsInitialized = 0;
}

void freeGlobalReusedMPFRVars() {
  unsigned int i;

  if (globalReusedMPFRVars == NULL) return;
  if (globalReusedMPFRVarsAllocated == 0) return;

  for (i=0;i<globalReusedMPFRVarsInitialized;i++) {
    mpfr_clear(globalReusedMPFRVars[i]);
  }
  safeFree(globalReusedMPFRVars);
  globalReusedMPFRVars = NULL;
  globalReusedMPFRVarsAllocated = 0;
  globalReusedMPFRVarsUsed = 0;
  globalReusedMPFRVarsInitialized = 0;
}

/* Compatibility functions */

void *sollya_memmove_impl(void *dest, const void *src, size_t n) {
  void *tmp;
  tmp = safeMalloc(n);
  memcpy(tmp, src, n);
  memcpy(dest, tmp, n);
  safeFree(tmp);
  return dest;
}

void *sollya_memset_impl(void *s, int c, size_t n) {
  char b;
  char *a;
  size_t i;
  b = (char) c;
  a = (char *) s;
  for (i=((size_t) 0);i<n;i++) {
    a[i] = b;
  }
  return s;
}

char *sollya_strchr_impl(const char *s, int c) {
  char b;
  char *curr;
  b = (char) c;
  for (curr=((char *) s);*curr!='\0';curr++) {
    if (*curr == b) return curr;
  }
  if (b == '\0') return curr;
  return NULL;
}

char *sollya_strrchr_impl(const char *s, int c) {
  char b;
  char *curr;
  b = (char) c;
  for (curr=((char *) s);*curr!='\0';curr++);
  for (;curr>=s;curr--) {
    if (*curr == b) return curr;
  }
  return NULL;
}

static inline int __sollya_strtol_impl_aux(uint64_t *v, char c, int b) {
  uint64_t t;
  if (b < 2) return 0;
  if (b > 36) return 0;
  if ((c >= '0') && (c <= '9')) {
    t = (uint64_t) (c - '0');
    if (t >= ((uint64_t) b)) return 0;
    *v = t;
    return 1;
  }
  if ((c >= 'A') && (c <= 'Z')) {
    t = ((uint64_t) (c - 'A')) + ((uint64_t) 10);
    if (t >= ((uint64_t) b)) return 0;
    *v = t;
    return 1;
  }
  if ((c >= 'a') && (c <= 'z')) {
    t = ((uint64_t) (c - 'z')) + ((uint64_t) 10);
    if (t >= ((uint64_t) b)) return 0;
    *v = t;
    return 1;
  }
  return 0;
}

long int sollya_strtol_impl(const char *nptr, char **endptr, int base) {
  int b;
  char *curr;
  int sign;
  uint64_t v, c;
  long int value;
  b = base;
  sign = 1;
  for (curr = ((char *) nptr);*curr!='\0';curr++) {
    if (!((*curr == ' ') ||
	  (*curr == '\f') ||
	  (*curr == '\n') ||
	  (*curr == '\r') ||
	  (*curr == '\t') ||
	  (*curr == '\v'))) break;
  }
  if (*curr == '+') {
    curr++;
  } else {
    if (*curr == '-') {
      sign = -1;
      curr++;
    }
  }
  if (base == 0) {
    if ((*curr == '0') && (*(curr + 1) == 'x')) {
      b = 16;
      curr += 2;
    } else {
      if (*curr == '0') {
	b = 8;
	curr++;
      } else {
	b = 10;
      }
    }
  } else {
    if (base == 16) {
      if ((*curr == '0') && (*(curr + 1) == 'x')) {
	curr += 2;
      }
    }
  }
  if ((b < 2) || (b > 36)) {
    if (endptr != NULL) {
      *endptr = curr;
    }
    return 0;
  }
  v = (uint64_t) 0;
  value = (long int) 0;
  for (;*curr!='\0';curr++) {
    if (!__sollya_strtol_impl_aux(&c, *curr, b))
      break;
    v = v * ((uint64_t) b) + c;
    if (sign > 0) {
      if (v > ((uint64_t) LONG_MAX)) {
	value = LONG_MAX;
	break;
      } else {
	value = (long int) v;
      }
    } else {
      if (v > (((uint64_t) (-(LONG_MIN + ((long int) 16)))) + ((uint64_t) 16))) {
	value = LONG_MIN;
	break;
      } else {
	value = (long int) v;
      }
    }
  }
  if (endptr != NULL) {
    *endptr = curr;
  }
  return value;
}

#define __SOLLYA_DUP2_IMPL_BUF_SIZE 512

/* This is a lazy reimplementation of dup2 that will not work all the
   time
*/
int sollya_dup2_impl(int oldfd, int newfd) {
  int extraFDs[__SOLLYA_DUP2_IMPL_BUF_SIZE];
  int extraFDsCount;
  int tmp;
  int i;
  int gotIt;
  int closingExtraOkay;

  /* Handle stupid input */
  if (oldfd == newfd) return newfd;

  /* Initialize index extraFDsCount of extra file descriptors */
  extraFDsCount = 0;
  
  /* Get a copy of oldfd, to see if it is a valid 
     file descriptor 
  */
  tmp = dup(oldfd);
  if (tmp == -1) {
    /* The oldfd file descriptor does not seem to be valid */
    return -1;
  }
  extraFDs[extraFDsCount] = tmp;
  extraFDsCount++;

  /* Try to get a copy of the newfd file descriptor. Hence we can
     handle closing errors on newfd.

     The call may fail, which is not necessarily an issue as the newfd
     file descriptor may also not be open.

  */
  tmp = dup(newfd);
  if (tmp != -1) {
    extraFDs[extraFDsCount] = tmp;
    extraFDsCount++;
  }

  /* Close the newfd file descriptor. 

     The call may fail, which is not necessarily an issue as the newfd
     file descriptor may also not be open.

  */
  close(newfd);

  /* Try several times to get the right newfd copy of oldfd. */
  gotIt = 0;
  for (i=0;i<__SOLLYA_DUP2_IMPL_BUF_SIZE-2;i++) {
    tmp = dup(oldfd);
    if (tmp != -1) {
      if (tmp == newfd) {
	gotIt = 1;
	break;
      } else {
	extraFDs[extraFDsCount] = tmp;
	extraFDsCount++;
      }
    }
  }

  /* Close the extra file descriptors */
  closingExtraOkay = 1;
  for (i=0;i<extraFDsCount;i++) {
    if (close(extraFDs[i]) != 0) {
      closingExtraOkay = 0;
    }
  }
  
  /* Handle the case when one of the extra file descriptors could not
     be closed 
  */
  if (!closingExtraOkay) {
    return -1;
  }

  /* If we got the newfd file descriptor, return it, otherwise return
     error.
  */
  if (gotIt) {
    return newfd;
  }
  return -1;
}

static inline char *__sollya_strstr_impl_aux(char *haystack, char *needle) {
  char *currHaystack, *currNeedle;
  for (currHaystack=haystack, currNeedle=needle;
       ((*currHaystack != '\0') && (*currNeedle != '\0'));
       currHaystack++, currNeedle++) {
    if (*currHaystack != *currNeedle) return NULL;
  }
  if ((*currHaystack == '\0') && (*currNeedle != '\0')) return NULL;
  return haystack;
}

char *sollya_strstr_impl(const char *haystack, const char *needle) {
  char *curr, *tmp;
  for (curr=((char *) haystack);*curr!='\0';curr++) {
    tmp = __sollya_strstr_impl_aux(curr, ((char *) needle));
    if (tmp != NULL) return tmp;
  }
  return NULL;
}


/* End of compatibility functions */

static int general(int argc, char *argv[]) {
  struct termios termAttr;
  volatile int parseAbort, executeAbort;
  int i;
  FILE *fd = NULL;
  struct rlimit rlim;
  char *error;
  int doNotModifyStackSize;
#if !defined(__CYGWIN__)
  int repeatSetRLimit;
#endif
  volatile int lastWasError;
  int finishedBeforeParsing;
  int argsArgRead;
  char **temp;
  int k;
  int sollyaOptions;
  int frameCorruptionPrinted;

  oldGMPMalloc = NULL;
  oldGMPRealloc = NULL;
  oldGMPFree = NULL;
  replacement_mp_set_memory_functions = NULL;
  replacement_mp_get_memory_functions = NULL;
  actualCalloc = calloc;
  actualMalloc = malloc;
  actualFree = free;
  actualRealloc = realloc;
  actualFreeWithSize = wrapSafeFree;
  actualReallocWithSize = wrapSafeRealloc;    
  messageCallback = NULL;
  libraryMode = 0;
  lastMessageCallbackResult = 1;
  lastMessageSuppressedResult = -1;
  doNotModifyStackSize = 0;
  inputFileOpened = 0;
  flushOutput = 0;
  oldAutoPrint = 0;
  oldExternalProcedurePrint = 0;
  printMode = PRINT_MODE_LEGACY;
  warnFile = NULL;

  argsArgRead = 0;
  argsArgv = NULL;
  argsArgc = 0;
  for (i=1;i<argc;i++) {
    if (argsArgRead) {
      argsArgc++;
      temp = (char **) safeCalloc(argsArgc, sizeof(char *));
      if (argsArgc > 1) {
	for (k=0;k<argsArgc-1;k++) {
	  temp[k] = argsArgv[k];
	}
	safeFree(argsArgv);
      }
      argsArgv = temp;
      argsArgv[argsArgc-1] = (char *) safeCalloc(strlen(argv[i])+1, sizeof(char));
      strcpy(argsArgv[argsArgc-1], argv[i]);
    } else {
      if (strcmp(argv[i],"--help") == 0) {
	if ((!inputFileOpened) && (!eliminatePromptBackup)) {
	  if (tcgetattr(0,&termAttr) == -1) {
	    eliminatePromptBackup = 1;
	  } else {
	    if (tcgetattr(1,&termAttr) == -1) {
	      eliminatePromptBackup = 1;
	    }
	  }
	} 
	sollyaPrintf("This is %s connected to ",PACKAGE_STRING);
	if (eliminatePromptBackup)
	  sollyaPrintf("a regular file");
	else
	  sollyaPrintf("a terminal");
	sollyaPrintf(".\n\nUsage: %s [options]\n\nPossible options are:\n",PACKAGE_NAME);
	sollyaPrintf("--args : Transmit the following arguments verbatim to the interpreter\n");
	sollyaPrintf("--donotmodifystacksize : do not attempt to set the maximal stack size to the maximum size allowed on the system\n");
	sollyaPrintf("--flush : flush standard output and standard error after each command\n");
	sollyaPrintf("--help : print this help text\n");
	sollyaPrintf("--nocolor : do not color the output using ANSI escape sequences\n");
	sollyaPrintf("--noprompt : do not print a prompt symbol\n");
	sollyaPrintf("--oldautoprint : print commas between autoprinted elements separated by commas\n");
	sollyaPrintf("--oldexternalprocprint : print the signature of an external procedure when autoprinting\n");
	sollyaPrintf("--oldrlwrapcompatible : acheive some compatibilty with old rlwrap versions by emitting wrong ANSI sequences (deprecated)\n");
	sollyaPrintf("--warninfile[append] <file> : print warning messages into a file instead on the standard output\n");
	sollyaPrintf("--warnonstderr : print warning messages on error output instead on the standard output\n");
	sollyaPrintf("\nFor help on %s commands type \"help;\" on the %s prompt\n",PACKAGE_NAME,PACKAGE_NAME);
	sollyaPrintf("More documentation on %s is available on the %s website http://sollya.gforge.inria.fr/.\nFor bug reports send an email to %s.\n",PACKAGE_NAME,PACKAGE_NAME,PACKAGE_BUGREPORT);
	sollyaPrintf("\n%s is\n\nCopyright 2006-2016 by\n\n   Laboratoire de l'Informatique du Parallelisme,\n   UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668,\n\n   LORIA (CNRS, INPL, INRIA, UHP, U-Nancy 2)\n\n   Laboratoire d'Informatique de Paris 6 - Équipe PEQUAN\n   Sorbonne Universités\n   UPMC Univ Paris 06\n   UMR 7606, LIP6\n   Boîte Courrier 169\n   4, place Jussieu\n   F-75252 Paris Cedex 05\n   France\n\n   CNRS, LIP6, UPMC\n   Sorbonne Universités, UPMC Univ Paris 06,\n   CNRS, LIP6 UMR 7606, 4 place Jussieu 75005 Paris\n\nand by\n\n    INRIA Sophia-Antipolis Mediterranee, APICS Team,\n    Sophia-Antipolis, France.\n\nAll rights reserved.\n\nContributors:\n              Ch. Lauter\n	      S. Chevillard\n	      M. Joldes\n              N. Jourdan\n	      M. Mezzarobba.\n\nThis software is governed by the CeCILL-C license under French law and\nabiding by the rules of distribution of free software.  You can  use,\nmodify and/ or redistribute the software under the terms of the CeCILL-C\nlicense as circulated by CEA, CNRS and INRIA at the following URL\n\"http://www.cecill.info\".\n\nAs a counterpart to the access to the source code and  rights to copy,\nmodify and redistribute granted by the license, users are provided only\nwith a limited warranty  and the software's author,  the holder of the\neconomic rights,  and the successive licensors  have only  limited\nliability.\n\nIn this respect, the user's attention is drawn to the risks associated\nwith loading,  using,  modifying and/or developing or reproducing the\nsoftware by the user in light of its specific status of free software,\nthat may mean  that it is complicated to manipulate,  and  that  also\ntherefore means  that it is reserved for developers  and  experienced\nprofessionals having in-depth computer knowledge. Users are therefore\nencouraged to load and test the software's suitability as regards their\nrequirements in conditions enabling the security of their systems and/or\ndata to be ensured and,  more generally, to use and operate it in the\nsame conditions as regards security.\n\nThe fact that you are presently reading this means that you have had\nknowledge of the CeCILL-C license and that you accept its terms.\n\nThis program is distributed WITHOUT ANY WARRANTY; without even the\nimplied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n\nThis build of %s is based on GMP %s, MPFR %s and MPFI %s.\n",PACKAGE_STRING,PACKAGE_STRING,gmp_version,mpfr_get_version(),sollya_mpfi_get_version());
#if defined(HAVE_FPLLL_VERSION_STRING)
	sollyaPrintf("%s uses FPLLL as: \"%s\"\n",PACKAGE_STRING,HAVE_FPLLL_VERSION_STRING);
#endif
	sollyaPrintf("\n");
	return 1;
      } else
	if (strcmp(argv[i],"--args") == 0) argsArgRead = 1; else
	  if (strcmp(argv[i],"--nocolor") == 0) noColor = 1; else
	    if (strcmp(argv[i],"--noprompt") == 0) eliminatePromptBackup = 1; else
	      if (strcmp(argv[i],"--oldrlwrapcompatible") == 0) oldrlwrapcompatible = 1; else
		if (strcmp(argv[i],"--flush") == 0) flushOutput = 1; else
		  if (strcmp(argv[i],"--oldautoprint") == 0) oldAutoPrint = 1; else
		    if (strcmp(argv[i],"--oldexternalprocprint") == 0) oldExternalProcedurePrint = 1; else
		      if (strcmp(argv[i],"--warnonstderr") == 0) {
			if (printMode != PRINT_MODE_WARNING_TO_FILE) {
			  printMode = PRINT_MODE_WARNING_TO_STDERR;
			} else {
			  sollyaPrintf("Error: only one of the --warnonstderr or --warninfile options can be given.\n");
			  return 1;
			}
		      } else
			if ((strcmp(argv[i],"--warninfile") == 0) ||
			    (strcmp(argv[i],"--warninfileappend") == 0)){
			  if (printMode != PRINT_MODE_WARNING_TO_STDERR) {
			    if (i+1<argc) {
			      i++;
			      if (strcmp(argv[i-1],"--warninfileappend") == 0) {
				fd = fopen(argv[i],"a");
			      } else {
				fd = fopen(argv[i],"w");
			      }
			      if (fd != NULL) {
				warnFile = fd;
				fd = NULL;
				printMode = PRINT_MODE_WARNING_TO_FILE;
			      } else {
				sollyaPrintf("Error: the file \"%s\" could not be opened for warning display: %s\n",argv[i],strerror(errno));
				return 1;
			      }
			    } else {
			      sollyaPrintf("Error: no file argument is given for the --warninfile option.\n");
			      return 1;
			    }
			  } else {
			    sollyaPrintf("Error: only one of the --warnonstderr or --warninfile options can be given.\n");
			    return 1;
			  }
			} else {
			  if (strcmp(argv[i],"--donotmodifystacksize") == 0) doNotModifyStackSize = 1; else {
			    if (!inputFileOpened) {
			      fd = fopen(argv[i],"r");
			      if (fd != NULL) {
				inputFile = fd;
				inputFileOpened = 1;
				eliminatePromptBackup = 1;
				noColor = 1;
			      } else {
				sollyaPrintf("Error: the file \"%s\" could not be opened: %s\n",argv[i],strerror(errno));
				return 1;
			      }
			    } else {
			      sollyaOptions = 0;
			      for (k=i;(k<argc)&&(!sollyaOptions);k++) {
				if ((strcmp(argv[k], "--args") == 0) ||
				    (strcmp(argv[k], "--donotmodifystacksize") == 0) ||
				    (strcmp(argv[k], "--flush") == 0) ||
				    (strcmp(argv[k], "--help") == 0) ||
				    (strcmp(argv[k], "--nocolor") == 0) ||
				    (strcmp(argv[k], "--noprompt") == 0) ||
				    (strcmp(argv[k], "--oldautoprint") == 0) ||
				    (strcmp(argv[k], "--oldexternalprocprint") == 0) ||
				    (strcmp(argv[k], "--oldrlwrapcompatible") == 0) ||
				    (strcmp(argv[k], "--warninfile") == 0) ||
				    (strcmp(argv[k], "--warninfileappend") == 0) ||
				    (strcmp(argv[k], "--warnonstderr") == 0)) {
				  sollyaOptions = 1;
				}
			      }
			      if (sollyaOptions) {
				sollyaPrintf("Error: another input file besides \"%s\" has been indicated and opened and the order of options is ambiguous.\n",argv[i]);
				return 1;
			      } else {
				argsArgRead = 1;
				argsArgc++;
				temp = (char **) safeCalloc(argsArgc, sizeof(char *));
				if (argsArgc > 1) {
				  for (k=0;k<argsArgc-1;k++) {
				    temp[k] = argsArgv[k];
				  }
				  safeFree(argsArgv);
				}
				argsArgv = temp;
				argsArgv[argsArgc-1] = (char *) safeCalloc(strlen(argv[i])+1, sizeof(char));
				strcpy(argsArgv[argsArgc-1], argv[i]);
			      }
			    }
			  }
			}
    }
  }

  if ((!inputFileOpened) && (!eliminatePromptBackup)) {
    if (tcgetattr(0,&termAttr) == -1) {
      eliminatePromptBackup = 1;
    } else {
      if (tcgetattr(1,&termAttr) == -1) {
	eliminatePromptBackup = 1;
      }
    }
  } 

  yylex_init(&scanner);
  if (inputFileOpened) {
    yyset_in(inputFile, scanner);
  }

  if (!doNotModifyStackSize) {
#if defined(RLIMIT_STACK)
    if (getrlimit(RLIMIT_STACK,&rlim) == 0) {
      rlim.rlim_cur = rlim.rlim_max;
      setrlimit(RLIMIT_STACK,&rlim);
#if !defined(__CYGWIN__)
      repeatSetRLimit = 0;
      if (getrlimit(RLIMIT_STACK,&rlim) == 0) {
        repeatSetRLimit = (rlim.rlim_cur != rlim.rlim_max);
      } else {
        repeatSetRLimit = 1;
      }
      if (repeatSetRLimit) {
        if (setrlimit(RLIMIT_STACK,&rlim) != 0) {
          if ((error = strerror(errno)) != NULL) {
            printMessage(1,SOLLYA_MSG_ERROR_ON_INITIAL_SETUP,"Warning: during initial setup, the following error occurred: \"%s\"\nTry using --donotmodifystacksize when invoking the tool.\n",error);
          } else {
            printMessage(1,SOLLYA_MSG_ERROR_ON_INITIAL_SETUP,"Warning: during initial setup, an unknown error occurred.\nTry using --donotmodifystacksize when invoking the tool.\n");
          }
        }
      }
#endif
    } else {
      if ((error = strerror(errno)) != NULL) {
	printMessage(1,SOLLYA_MSG_ERROR_ON_INITIAL_SETUP,"Warning: during initial setup, the following error occurred: \"%s\"\nTry using --donotmodifystacksize when invoking the tool.\n",error);
      } else {
	printMessage(1,SOLLYA_MSG_ERROR_ON_INITIAL_SETUP,"Warning: during initial setup, an unknown error occurred.\nTry using --donotmodifystacksize when invoking the tool.\n");
      }
    }
#endif
  }
  initSignalHandler(eliminatePromptBackup);
  blockSignals(eliminatePromptBackup);
  wrap_mp_set_memory_functions(safeMalloc,clingWrapSafeReallocWithSize,clingWrapSafeFreeWithSize);
  initToolDefaults();

  exitInsteadOfRecover = 0;
  helpNotFinished = 0;
  printPrompt();
  lastWasError = 0;
  lastCorrectlyExecuted = 0;
  while (1) {
    if (lastHandledSignal != 0) {
      switch (lastHandledSignal) {
      case HANDLING_SIGINT:
	break;
      case HANDLING_SIGSEGV:
	printMessage(1,SOLLYA_MSG_HANDLED_SIGSEGV,"Warning: a SIGSEGV signal has been handled.\n");
	break;
      case HANDLING_SIGBUS:
	printMessage(1,SOLLYA_MSG_HANDLED_SIGBUS,"Warning: a SIGBUS signal has been handled.\n");
	break;
      case HANDLING_SIGFPE:
	printMessage(1,SOLLYA_MSG_HANDLED_SIGFPE,"Warning: a SIGFPE signal has been handled.\n");
	break;
      case HANDLING_SIGPIPE:
	printMessage(1,SOLLYA_MSG_HANDLED_SIGPIPE,"Warning: a SIGPIPE signal has been handled.\n");
	break;
      default:
	break;
      }
      lastHandledSignal = 0;
    }
    executeAbort = 0;
    parsedThing = NULL;
    lastWasSyntaxError = 0;
    finishedBeforeParsing = parserCheckEof();
    (void) finishedBeforeParsing;
    parseAbort = yyparse(scanner);
    lastWasError = 0;
    if (parsedThing != NULL) {

      handlingCtrlC = 0;
      lastHandledSignal = 0;
      if (!setjmp(recoverEnvironment)) {
	recoverEnvironmentReady = 1;
	frameCorruptionPrinted = 0;
	if (declaredSymbolTable != NULL) {
	  printMessage(1,SOLLYA_MSG_FRAME_STACK_HAS_BEEN_CORRUPTED,"Warning: a preceeding command interruption corrupted the variable frame stack.\n");
	  frameCorruptionPrinted = 1;
	  freeDeclaredSymbolTable(declaredSymbolTable, freeThingOnVoid);
	  declaredSymbolTable = NULL;
	}
	if (backtraceStack != NULL) {
	  if (!frameCorruptionPrinted) printMessage(1,SOLLYA_MSG_FRAME_STACK_HAS_BEEN_CORRUPTED,"Warning: a preceeding command interruption corrupted the variable frame stack.\n");
	  freeBacktraceStack();
	  backtraceStack = NULL;
	}
	initSignalHandler(eliminatePromptBackup);
	numberBacktrace = 1;
	if (timeStack != NULL) {
	  printMessage(4,SOLLYA_MSG_TIMING_STACK_HAS_BEEN_CORRUPTED,"Information: corrupted timing stack. Releasing the stack.\n");
	  freeCounter();
	}
        if (flushOutput) {
	  deferSignalHandling();
          fflush(stdout);
          fflush(stderr);
	  resumeSignalHandling();
        }
	libraryMode = 0;
	pushTimeCounter();
	executeAbort = executeCommand(parsedThing);
	lastCorrectlyExecuted = 1;
	popTimeCounter("full execution of the last parse chunk");
	if((!timecounting) && (timeStack!=NULL)) {
	  freeCounter();
	}
        if (flushOutput) {
	  deferSignalHandling();
          fflush(stdout);
          fflush(stderr);
	  resumeSignalHandling();
        }
	blockSignals(eliminatePromptBackup);
	recoverEnvironmentReady = 0;
      } else {
	displayColor = -1; normalMode();
	blockSignals(eliminatePromptBackup);
	lastWasError = 1;
	if (lastHandledSignal != 0) {
	  switch (lastHandledSignal) {
	  case HANDLING_SIGINT:
	    break;
	  case HANDLING_SIGSEGV:
	    printMessage(1,SOLLYA_MSG_HANDLED_SIGSEGV,"Warning: a SIGSEGV signal has been handled.\n");
	    break;
	  case HANDLING_SIGBUS:
	    printMessage(1,SOLLYA_MSG_HANDLED_SIGBUS,"Warning: a SIGBUS signal has been handled.\n");
	    break;
	  case HANDLING_SIGFPE:
	    printMessage(1,SOLLYA_MSG_HANDLED_SIGFPE,"Warning: a SIGFPE signal has been handled.\n");
	    break;
	  case HANDLING_SIGPIPE:
	    printMessage(1,SOLLYA_MSG_HANDLED_SIGPIPE,"Warning: a SIGPIPE signal has been handled.\n");
	    break;
	  default:
	    break;
	  }
	  lastHandledSignal = 0;
	}
	if (handlingCtrlC)
	  printMessage(1,SOLLYA_MSG_LAST_COMMAND_INTERRUPTED,"Warning: the last command has been interrupted. May leak memory.\n");
	else {
	  printMessage(1,SOLLYA_MSG_COMMAND_NOT_EXECUTABLE,"Warning: the last command could not be executed. May leak memory.\n");
          considerDyingOnError();
        }
	frameCorruptionPrinted = 0;
	if (declaredSymbolTable != NULL) {
	  if (!handlingCtrlC)
	    printMessage(1,SOLLYA_MSG_RELEASING_FRAME_STACK,"Warning: releasing the variable frame stack.\n");
	  else
	    printMessage(2,SOLLYA_MSG_RELEASING_FRAME_STACK,"Information: releasing the variable frame stack.\n");
	  frameCorruptionPrinted = 1;
	  freeDeclaredSymbolTable(declaredSymbolTable, freeThingOnVoid);
	}
	declaredSymbolTable = NULL;
	if (backtraceStack != NULL) {
	  if (!frameCorruptionPrinted) {
	    if (!handlingCtrlC)
	      printMessage(1,SOLLYA_MSG_RELEASING_FRAME_STACK,"Warning: releasing the variable frame stack.\n");
	    else
	      printMessage(2,SOLLYA_MSG_RELEASING_FRAME_STACK,"Information: releasing the variable frame stack.\n");
	  }
	  freeBacktraceStack();
	}
	backtraceStack = NULL;
	if (timeStack != NULL) {
	  printMessage(2,SOLLYA_MSG_TIMING_STACK_HAS_BEEN_CORRUPTED,"Information: corrupted timing stack. Releasing the stack.\n");
	  freeCounter();
	}
      }

      freeThing(parsedThing);
    } else {
      lastWasError = 1;
    }
    if (parseAbort || executeAbort) break;
    promptToBePrinted = 1;
    helpNotFinished = 0;
  }

  freeTool();

  if (argsArgv != NULL) {
    for (k=0;k<argsArgc;k++) {
      safeFree(argsArgv[k]);
    }
    safeFree(argsArgv);
    argsArgv = NULL;
  }
  
  if (!eliminatePromptBackup) sollyaPrintf("\n");

  if (inputFileOpened) {
    fclose(inputFile);
    inputFileOpened = 0;
  }

  if (warnFile != NULL) {
    fclose(warnFile);
    warnFile = NULL;
  }

  sollya_mp_set_memory_functions(oldGMPMalloc,oldGMPRealloc,oldGMPFree);

  if (lastWasError) {
    if (lastCorrectlyExecuted) {
      return 3;
    } else {
      return 4;
    }
  } else {
    return 0;
  }
}

int sollya_tool_main(int argc, char *argv[]) {
  return general(argc, argv);
}

#undef malloc
#undef realloc

void *rpl_malloc(size_t n) {
  if (n == 0)
    n = 1;
  return malloc (n);
}

void *rpl_realloc(void *ptr, size_t n) {
  if (n == 0)
    n = 1;
  return realloc (ptr, n);
}

