/*

  Copyright 2008-2010 by

  Laboratoire de l'Informatique du Parallelisme,
  UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668

  and by

  LORIA (CNRS, INPL, INRIA, UHP, U-Nancy 2)

  Contributors S. Chevillard, M. Joldes, Ch. Lauter

  sylvain.chevillard@ens-lyon.org
  joldes@laas.fr
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

#ifndef AUTODIFF_H
#define AUTODIFF_H

#include "mpfi-compat.h"
#include "expression.h"

void symbolic_poly_diff(sollya_mpfi_t *res, sollya_mpfi_t *coeff_array, int degree);
void symbolic_poly_evaluation_horner(sollya_mpfi_t res, sollya_mpfi_t *coeffs_array, sollya_mpfi_t x, int degree);
void symbolic_poly_evaluation_powers(sollya_mpfi_t res, sollya_mpfi_t *coeffs_array, sollya_mpfi_t *powers_array, sollya_mpfi_t x, int degree);

void libraryFunction_diff(sollya_mpfi_t *res, node *f, sollya_mpfi_t x, int n, int *silent);
void procedureFunction_diff(sollya_mpfi_t *res, node *f, sollya_mpfi_t x, int n, int *silent);

void powerFunction_diff(sollya_mpfi_t *res, mpfr_t p, sollya_mpfi_t x0, int n, int *silent);
void constantPower_diff(sollya_mpfi_t *res, sollya_mpfi_t x0, mpfr_t p, int n, int *silent);

void multiplication_AD(sollya_mpfi_t *res, sollya_mpfi_t *f, sollya_mpfi_t *g, int n);
void composition_AD(sollya_mpfi_t *res, sollya_mpfi_t *g, sollya_mpfi_t *f, int n);
void auto_diff_scaled(sollya_mpfi_t* res, node *f, sollya_mpfi_t x0, int n);
void auto_diff(sollya_mpfi_t* res, node *f, sollya_mpfi_t x0, int n);

#endif /* AUTODIFF_H */
