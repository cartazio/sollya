/*
Copyright 2011-2013 by 
  
  Centre de recherche INRIA Sophia-Antipolis Mediterranee, equipe APICS,
  Sophia Antipolis, France.
  
  Laboratoire de l'Informatique du Parallelisme, 
  UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668

  Laboratoire d'Informatique de Paris 6, equipe PEQUAN,
  UPMC Universite Paris 06 - CNRS - UMR 7606 - LIP6, Paris, France,

  Contributors S. Chevillard, M. Joldes, Ch. Lauter 

  sylvain.chevillard@ens-lyon.org
  mioara.joldes@ens-lyon.fr
  christoph.lauter@ens-lyon.org
  

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

#ifndef CHEBYSHEVFORM_H
#define CHEBYSHEVFORM_H

#include "expression.h"
#include "chain.h"
#include "chebModelsAux.h"
/*Cheb model (CM, cm) structure:
n - order: polynomial of degree n-1, remainder of order O(x^n)
x - interval on which the cm is computed
cheb_array - n chebyshev points (of first kind) kept; denoted x_i
             these are computed only once for each interval x
cheb_matix - T_j(x_i);  i, j in {0,..., n-1}; 
             This has to be computed only once on each interval x                       
poly_array - array of n coeffs (in Cheb Basis)
             corresponding to polynomial of degree n-1
rem_bound - bound for the remainder
poly_bound - bound for the polynomial (helpful for computations)
*/
typedef struct cmdl {
int n; 
sollya_mpfi_t x;
sollya_mpfi_t **cheb_array;
sollya_mpfi_t **cheb_matrix;
sollya_mpfi_t *poly_array;
sollya_mpfi_t rem_bound;
sollya_mpfi_t poly_bound;
} chebModel;




void chebyshevform(node **Ch, chain **errors, sollya_mpfi_t delta, 
		   chain **chebyshevCoefficients, node *f, int n, 
		   sollya_mpfi_t dom, mp_prec_t prec);

#endif /* CHEBYSHEVFORM_H */
