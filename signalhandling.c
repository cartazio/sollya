/*

  Copyright 2013-2016 by

  Laboratoire d'Informatique de Paris 6 - Équipe PEQUAN
  Sorbonne Universités
  UPMC Univ Paris 06
  UMR 7606, LIP6
  Boîte Courrier 169
  4, place Jussieu
  F-75252 Paris Cedex 05
  France

  Contributors Ch. Lauter

  christoph.lauter@ens-lyon.org

  This software is a computer program whose purpose is to provide an
  environment for safe floating-point code development. It is
  particularly targeted to the automated implementation of
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

#define _POSIX_SOURCE
#define _POSIX_C_SOURCE (200312L)

#include <stdlib.h>
#include <signal.h>
#include <setjmp.h>
#include <stdio.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "signalhandling.h"

#define UNUSED_PARAM(_unused_param_x) ((void)(_unused_param_x))

extern jmp_buf recoverEnvironment;
extern int handlingCtrlC;
extern int lastHandledSignal;
extern int recoverEnvironmentReady;
extern int exitInsteadOfRecover;
extern int libraryMode;

int blockedSignalCounter = 0;

int deferredMode = 0;
int deferredSignal = 0;
int deferredSignalIsDeferred = 0;

int deferredCount = 0;

int emulatedBlockedSignalCounter = 0;

int sollyaFprintf(FILE *fd, const char *format, ...);

void signalHandler(int i) {
  
  if (deferredMode) {
    if (deferredSignalIsDeferred) return;
    deferredSignal = i;
    deferredSignalIsDeferred = 1;
    return;
  }

  switch (i) {
  case SIGINT:
    handlingCtrlC = 1;
    lastHandledSignal = HANDLING_SIGINT;
    break;
  case SIGSEGV:
    lastHandledSignal = HANDLING_SIGSEGV;
    break;
  case SIGBUS:
    lastHandledSignal = HANDLING_SIGBUS;
    break;
  case SIGFPE:
    lastHandledSignal = HANDLING_SIGFPE;
    break;
  case SIGPIPE:
    lastHandledSignal = HANDLING_SIGPIPE;
    break;
  default:
    sollyaFprintf(stderr,"Error: must handle an unknown signal.\n");
    exit(1);
  }
  if (recoverEnvironmentReady) {
    if (exitInsteadOfRecover) {
      sollyaFprintf(stderr,"Error: the recover environment has not been initialized. Exiting.\n");
      exit(1);
    }
    longjmp(recoverEnvironment,1);
  }
}

RETSIGTYPE signalHandlerForSignal(int i) {
  if (!emulatedBlockedSignalCounter) {
    signalHandler(i);
  }
  return (RETSIGTYPE) 1; /* RETSIGTYPE is either int or void. When it
			    is int, nothing happens. When it is void,
			    C supports that cast and returns
			    nothing. 
			 */
}

static inline void deferSignalHandlingInner() {
  if (deferredMode) return;
  deferredMode = 1;
  deferredSignal = 0;
  deferredSignalIsDeferred = 0;
}

void deferSignalHandling() {
  deferredCount++;
  if (deferredCount > 0) 
    deferSignalHandlingInner();
}

static inline void resumeSignalHandlingInner() {
  if (!deferredMode) return;
  deferredMode = 0;
  if (!deferredSignalIsDeferred) return;
  deferredSignalIsDeferred = 0;
  signalHandler(deferredSignal);
  deferredSignal = 0;
}

void resumeSignalHandling() {
  deferredCount--;
  if (deferredCount <= 0) 
    resumeSignalHandlingInner();
}

void initSignalHandler(int nointeract) {
#if (defined(HAVE_SIGACTION) && HAVE_SIGACTION) && (defined(HAVE_SIGADDSET) && HAVE_SIGADDSET) && (defined(HAVE_SIGEMPTYSET) && HAVE_SIGEMPTYSET) && (defined(HAVE_SIGPROCMASK) && HAVE_SIGPROCMASK) && (!defined(__CYGWIN__)) 

  sigset_t mask;
  struct sigaction action;
  
#endif

  
  blockedSignalCounter = 0;

  if (libraryMode) return;

#if (defined(HAVE_SIGACTION) && HAVE_SIGACTION) && (defined(HAVE_SIGADDSET) && HAVE_SIGADDSET) && (defined(HAVE_SIGEMPTYSET) && HAVE_SIGEMPTYSET) && (defined(HAVE_SIGPROCMASK) && HAVE_SIGPROCMASK) && (!defined(__CYGWIN__)) 
  
  action.sa_handler = signalHandler;
  action.sa_flags = 0;
  sigemptyset(&(action.sa_mask));
  if (!nointeract) sigaddset(&(action.sa_mask),SIGINT);
  sigaddset(&(action.sa_mask),SIGSEGV);
  sigaddset(&(action.sa_mask),SIGBUS);
  sigaddset(&(action.sa_mask),SIGFPE);
  sigaddset(&(action.sa_mask),SIGPIPE);
  if (!nointeract) sigaction(SIGINT, &action, NULL);
  sigaction(SIGSEGV, &action, NULL);
  sigaction(SIGBUS, &action, NULL);
  sigaction(SIGFPE, &action, NULL);
  sigaction(SIGPIPE, &action, NULL);

  sigemptyset(&mask);
  if (!nointeract) sigaddset(&mask,SIGINT);
  sigaddset(&mask,SIGSEGV);
  sigaddset(&mask,SIGBUS);
  sigaddset(&mask,SIGFPE);
  sigaddset(&mask,SIGPIPE);
  sigprocmask(SIG_UNBLOCK, &mask, NULL);

#else
  emulatedBlockedSignalCounter = 1;
  if (!nointeract) {
    if (signal(SIGINT, signalHandlerForSignal) == SIG_IGN) {
      signal(SIGINT, SIG_IGN);
    }
  }
  if (signal(SIGSEGV, signalHandlerForSignal) == SIG_IGN) {
      signal(SIGSEGV, SIG_IGN);
  }
  if (signal(SIGBUS, signalHandlerForSignal) == SIG_IGN) {
      signal(SIGBUS, SIG_IGN);
  }
  if (signal(SIGFPE, signalHandlerForSignal) == SIG_IGN) {
      signal(SIGFPE, SIG_IGN);
  }
  if (signal(SIGPIPE, signalHandlerForSignal) == SIG_IGN) {
      signal(SIGPIPE, SIG_IGN);
  }
  emulatedBlockedSignalCounter = 0;
#endif
  
}

void blockSignals(int nointeract) {
#if (defined(HAVE_SIGACTION) && HAVE_SIGACTION) && (defined(HAVE_SIGADDSET) && HAVE_SIGADDSET) && (defined(HAVE_SIGEMPTYSET) && HAVE_SIGEMPTYSET) && (defined(HAVE_SIGPROCMASK) && HAVE_SIGPROCMASK) && (!defined(__CYGWIN__))
  
  sigset_t mask;

#endif
  
  blockedSignalCounter = 0;

  if (libraryMode) return;

#if (defined(HAVE_SIGACTION) && HAVE_SIGACTION) && (defined(HAVE_SIGADDSET) && HAVE_SIGADDSET) && (defined(HAVE_SIGEMPTYSET) && HAVE_SIGEMPTYSET) && (defined(HAVE_SIGPROCMASK) && HAVE_SIGPROCMASK) && (!defined(__CYGWIN__))
  
  sigemptyset(&mask);
  if (!nointeract) sigaddset(&mask,SIGINT);
  sigaddset(&mask,SIGSEGV);
  sigaddset(&mask,SIGBUS);
  sigaddset(&mask,SIGFPE);
  sigaddset(&mask,SIGPIPE);
  sigprocmask(SIG_BLOCK, &mask, NULL);

#else
  UNUSED_PARAM(nointeract);
  emulatedBlockedSignalCounter = 1;

#endif
  
}

void initSignalHandlerCounted(int nointeract) {
  blockedSignalCounter--;
  if (blockedSignalCounter < 0) blockedSignalCounter = 0;
  if (blockedSignalCounter == 0) initSignalHandler(nointeract);
}

void blockSignalsCounted(int nointeract) {
  if (blockedSignalCounter < 0) blockedSignalCounter = 0;
  if (blockedSignalCounter == 0) blockSignals(nointeract);
  blockedSignalCounter++;
}




