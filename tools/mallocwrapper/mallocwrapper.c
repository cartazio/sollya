/*

  Copyright 2012 by

  Laboratoire d'Informatique de Paris 6, equipe PEQUAN,
  UPMC Universite Paris 06 - CNRS - UMR 7606 - LIP6, Paris, France,

  Contributors Ch. Lauter

  christoph.lauter@ens-lyon.org

  This software allows for testing the correct use of the malloc,
  calloc, realloc and free system functions in the Sollya tool.

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
  herefore means  that it is reserved for developers  and  experienced
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

/* Compile this file like that:

   gcc -fPIC -Wall -c mallocwrapper.c 
   gcc -fPIC -shared -o mallocwrapper.so mallocwrapper.o -ldl

   You will need a GNU compatible compiler that accepts the definition
   of _GNU_SOURCE or one that supports RTLD_NEXT without _GNU_SOURCE
   being supported.

   Then use it doing stuff like

   export LD_LIBRARY_PATH=`pwd`
   export LD_PRELOAD=mallocwrapper.so 
   ls 

   You will see a line on stderr for any memory management function call your
   program does.

*/

#define _GNU_SOURCE 

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <dlfcn.h>


/* Function pointers to the system memory management functions */
void *(*systemMalloc)(size_t) = NULL;
void *(*systemCalloc)(size_t, size_t) = NULL;
void *(*systemRealloc)(void *, size_t) = NULL;
void (*systemFree)(void *) = NULL;

/* A global variable holding state if we are in dlsym or not */
int executingDlsym = 0;

/* A pretty huge memory block for the simplistic bootstrap
   implementations of the memory management functions below 
*/

#define BOOTSTRAP_MEMORY_SIZE 1048576

unsigned char bootstrapMemory[BOOTSTRAP_MEMORY_SIZE];
unsigned char *firstFreeInBootstrapMemory = bootstrapMemory;

/* Manual, very simplistic implementations of the functions needed to bootstrap the process */

void *manualMalloc(size_t size) {
  unsigned char *ptr;
  unsigned char *newFirstFreeInBootstrapMemory;

  newFirstFreeInBootstrapMemory = firstFreeInBootstrapMemory + size;
  if (newFirstFreeInBootstrapMemory >= bootstrapMemory + BOOTSTRAP_MEMORY_SIZE) {
    return NULL;
  }
  
  ptr = firstFreeInBootstrapMemory;
  firstFreeInBootstrapMemory = newFirstFreeInBootstrapMemory;

  return (void *) ptr;
}

void *manualCalloc(size_t nmemb, size_t size) {
  void *ptr;
  size_t allSize;

  allSize = nmemb * size;
  ptr = malloc(allSize);
  if (ptr != NULL) memset(ptr, 0, allSize);

  return ptr;
}

void *manualRealloc(void *ptr, size_t size) {
  void *newPtr;

  newPtr = malloc(size);
  if (newPtr != NULL) {
    memcpy(newPtr, ptr, size);
  } 
  free(ptr);

  return newPtr;
}

void manualFree(void *ptr) {
  return;
}

/* A function to check if a pointer comes out of the bootstrap memory or not */
int pointerComesFromBootstrapMemory(void *ptr) {
  unsigned char *myPtr = (unsigned char *) ptr;
  return ((bootstrapMemory <= myPtr) && (myPtr < bootstrapMemory + BOOTSTRAP_MEMORY_SIZE));
}

/* The actual replacement functions */

void *mallocInner(size_t size) {
  void *newPtr;

  /* Initialize the actual system malloc */
  if (systemMalloc == NULL) {
    if (executingDlsym > 0) return manualMalloc(size);
    executingDlsym++;
    systemMalloc = dlsym(RTLD_NEXT,"malloc");
    executingDlsym--;
    if (systemMalloc == NULL) {
      return NULL;
    }
  }

  /* Do the system malloc */
  newPtr = systemMalloc(size);

  return newPtr;
}

void *malloc(size_t size) {
  void *newPtr;

  newPtr = mallocInner(size);

  fprintf(stderr, "malloc(0x%zx) = %p\n", size, newPtr);

  return newPtr;
}

void *callocInner(size_t nmemb, size_t size) {
  void *newPtr;

  /* Initialize the actual system calloc */
  if (systemCalloc == NULL) {
    if (executingDlsym > 0) return manualCalloc(nmemb,size);
    executingDlsym++;
    systemCalloc = dlsym(RTLD_NEXT,"calloc");
    executingDlsym--;
    if (systemCalloc == NULL) {
      return NULL;
    }
  }

  /* Do the system calloc */
  newPtr = systemCalloc(nmemb,size);

  return newPtr;
}

void *calloc(size_t nmemb, size_t size) {
  void *newPtr;

  newPtr = callocInner(nmemb, size);

  fprintf(stderr, "calloc(0x%zx, 0x%zx) = %p\n", nmemb, size, newPtr);

  return newPtr;
}

void *reallocInner(void *ptr, size_t size) {
  void *newPtr;

  /* Initialize the actual system realloc */
  if (systemRealloc == NULL) {
    if (executingDlsym > 0) return manualRealloc(ptr,size);
    executingDlsym++;
    systemRealloc = dlsym(RTLD_NEXT,"realloc");
    executingDlsym--;
    if (systemRealloc == NULL) {
      free(ptr);
      return NULL;
    }
  }

  /* Check if the pointer to free comes from the bootstrap memory
     pool. If yes, replace it with memory from the system malloc.
  */
  if (pointerComesFromBootstrapMemory(ptr)) {
    newPtr = malloc(size);
    if (newPtr != NULL) {
      memcpy(newPtr, ptr, size);
    } 
    manualFree(ptr);
    return ptr;
  }

  /* Do the system realloc */
  newPtr = systemRealloc(ptr,size);

  return newPtr;
}

void *realloc(void *ptr, size_t size) {
  void *newPtr;

  newPtr = reallocInner(ptr, size);

  fprintf(stderr, "realloc(%p, 0x%zx) = %p\n", ptr, size, newPtr);

  return newPtr;
}

void freeInner(void *ptr) {

  /* Initialize the actual system free */
  if (systemFree == NULL) {
    if (executingDlsym > 0) return manualFree(ptr);
    executingDlsym++;
    systemFree = dlsym(RTLD_NEXT,"free");
    executingDlsym--;
    if (systemFree == NULL) {
      return;
    }
  }

  /* Check if the pointer to free comes from the bootstrap memory
     pool. If yes, do not actually free it with the system free.
  */
  if (pointerComesFromBootstrapMemory(ptr)) {
    manualFree(ptr);
  }

  /* Do the system free */
  systemFree(ptr);
}

void free(void *ptr) {

  freeInner(ptr);

  fprintf(stderr, "free(%p)\n", ptr);
}

