/*

Copyright 2008 by 

Laboratoire de l'Informatique du Parallélisme, 
UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668

Contributors Ch. Lauter, S. Chevillard, N. Jourdan

christoph.lauter@ens-lyon.org
sylvain.chevillard@ens-lyon.org
nicolas.jourdan@ens-lyon.fr

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

#ifndef REMEZ_H
#define REMEZ_H



#include <mpfr.h>
#include "expression.h"
#include "remez.h"
#include "chain.h"
#include "infnorm.h"

chain *uncertifiedFindZeros(node *tree, mpfr_t a, mpfr_t b, unsigned long int points, mp_prec_t prec);

node* remez(node *func, node *weight, chain* monom, mpfr_t a, mpfr_t b, mpfr_t *requestedQuality, mp_prec_t prec);
node* remezWithWeight(node *func, node *weight, chain *monomials, mpfr_t a, mpfr_t b, mp_prec_t prec);

rangetype guessDegree(node *func, node *weight, mpfr_t a, mpfr_t b, mpfr_t eps);
node *constructPolynomial(mpfr_t *coeff, chain *monomials, mp_prec_t prec);

int newtonFaithful(mpfr_t res, node *f, node *f_diff, mpfr_t a, mpfr_t b, mp_prec_t prec);

#endif /* ifdef REMEZ_H*/
