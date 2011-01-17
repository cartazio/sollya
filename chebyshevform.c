/*

Copyright 2011 by 

Laboratoire de l'Informatique du Parallélisme, 
UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668

and by

Laboratoire d'Informatique de Paris 6, equipe PEQUAN,
UPMC Universite Paris 06 - CNRS - UMR 7606 - LIP6, Paris, France.

Contributors Ch. Lauter, M. Joldes

christoph.lauter@ens-lyon.org
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


#include "chebyshevform.h"
#include "external.h"
#include "infnorm.h"
#include "autodiff.h"
#include "general.h"
#include <stdio.h>
#include <stdlib.h>


/* Compute a Chebyshev form of order n for f over dom 
   
   Make a base change such that Ch + sum errors_i * x^i + delta is an
   absolute-error Taylor form for f over dom, developed at a point in
   the middle of dom.
   
   Return the original Chebyshev coefficients (in the Chebyshev basis)
   as chebyshevCoefficients.

   Upon entry, delta has already been initialized with sollya_mpfi_init2.
   The (sollya_mpfi_t *) intervals in errors and chebyshevCoefficients 
   must be allocated (safeMalloc'ed) and initialized with sollya_mpfi_init2.

   The entries in errors and chebyshevCoefficients are in natural order,
   i.e. errors->value points to the 0-th (constant) radius interval
   around the 0-th (constant) coefficient in Ch.

   An error during computation can be signaled by setting Ch to NULL.
   In this case, no memory must be left allocated by the function 
   upon return.

   All computations should be performed using prec as working precision
   (try to refrain from using getToolsPrecision()).

*/
void chebyshevform(node **Ch, chain **errors, sollya_mpfi_t delta, 
		   chain **chebyshevCoefficients,
		   node *f, int n, sollya_mpfi_t dom, mp_prec_t prec) {

  /* Adjust n to the notion of degree in the taylor command */
  n++;

  /* Check if degree is at least 1, once it has been adjusted */
  if (n < 1) {
    printMessage(1,"Warning: the degree of a Chebyshev Model must be at least 0.\n");
    *Ch = NULL;
    return;
  } 

  /* TODO */
  *Ch = NULL;

}
