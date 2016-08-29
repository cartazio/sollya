/*

  Copyright 2008-2016 by

  Laboratoire de l'Informatique du Parallelisme,
  UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668,

  LORIA (CNRS, INPL, INRIA, UHP, U-Nancy 2),

  Centre de recherche INRIA Sophia-Antipolis Mediterranee, equipe APICS,
  Sophia Antipolis, France.

  and by

  Laboratoire d'Informatique de Paris 6 - Équipe PEQUAN
  Sorbonne Universités
  UPMC Univ Paris 06
  UMR 7606, LIP6
  Boîte Courrier 169
  4, place Jussieu
  F-75252 Paris Cedex 05
  France

  Contributors Ch. Lauter, S. Chevillard

  christoph.lauter@ens-lyon.org
  sylvain.chevillard@ens-lyon.org

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


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#undef malloc
#undef realloc

#ifdef HAVE_SPECIAL_FPLLL_INCLUDE
#if HAVE_SPECIAL_FPLLL_INCLUDE
#include <fplll/fplll.h>
#else
#include <fplll.h>
#endif
#else
#include <fplll.h>
#endif

#include <gmp.h>


#define coeff(i,j,n) ((i)-1)*(n)+(j)-1

/* Considering a matrix of nbpoints+1 x dim+1 rational numbers, performs the
   LLL reduction of the matrix (seen as a basis of row vectors) and returns the
   last reduced vector (which is hence a nbpoints+1 array of mpq_t).
   reducedVect is supposed to point to a large enough segment of allocated and
   already initialized mpq_t
*/
extern "C" void fplll_wrapper(mpq_t *reducedVect, mpq_t *exactMatrix, int dim, int nbpoints) {
  int i,j;
  ZZ_mat<mpz_t> * FPlllMat;
  Z_NR<mpz_t>  zval;
  wrapper *LLLwrapper;
  mpz_t mpzval;

  mpz_init(mpzval);

  FPlllMat = new ZZ_mat<mpz_t>(dim+1,nbpoints+1);

  for(j=1; j<=dim+1; j++) {
    for(i=1; i<=nbpoints+1; i++) {
      mpz_set_q(mpzval, exactMatrix[coeff(i,j,dim+1)]); /* Casts M[i,j] into a mpz_t */
      zval.set(mpzval);
      FPlllMat->Set(j-1,i-1,zval);
    }
  }

  // LLL reduction
  LLLwrapper = new wrapper(FPlllMat);
  LLLwrapper->LLL();

  // Converting all stuff into exact numbers
  for(i=1; i<=nbpoints+1; i++) {
    mpq_set_z(reducedVect[i-1], LLLwrapper->GetBase()->Get(dim, i-1).GetData());
  }

  // Cleaning
  delete FPlllMat;
  delete LLLwrapper;

  mpz_clear(mpzval);
  return;
}

