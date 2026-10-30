/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2026 Guido Dhondt                     */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include "CalculiX.h"

/**
 * densify
 * =======
 *
 * Convert a sparse matrix into a dense matrix.
 *
 * Usage
 * -----
 * double *A;
 * NNEW(A, double, *n **n);
 * densify(A, n, ad, au, irow, jq, nzs, inputformat);
  *
 * Authors
 * -------
 * 2025-09-25 Christoph Woelfle - initial commit
 *
 * Parameters
 * ----------
 * @param A output dense matrix
 * @param n row and column dimension of input and output matrix
 * @param ad diagonal elements of input matrix
 * @param au offdiagonal elements of input matrix
 * @param irow row indices of elements in au
 * @param jq index in irow of first nonzero of each column
 * @param nzs number of nonzeroes
 * @param inputformat flag for input matrix structure
 * ...
 */
void densify(double *A, const ITG *n, const double *ad, const double *au,
              const ITG *irow, const ITG *jq, const ITG *nzs,
              const ITG *inputformat) {
  /* convert a sparse matrix to a dense matrix */

  /* diagonal */
  for (ITG i = 0; i < *n; i++) {
    A[i * (*n + 1)] = ad[i];
  }

  /* offdiagonal */
  if (*inputformat == 0) {
    /* symmetric matrix */
    for (ITG j = 0; j < *n; j++) {
      for (ITG i = jq[j]; i < jq[j + 1]; i++) {
        A[j * *n + irow[i - 1] - 1] = au[i - 1];
        A[(irow[i - 1] - 1) * *n + j] = au[i - 1];
      }
    }
  } else if (*inputformat == 1) {
    /* asymmetric matrix with symmetric sparsity pattern*/
    for (ITG j = 0; j < *n; j++) {
      for (ITG i = jq[j]; i < jq[j + 1]; i++) {
        A[j * *n + irow[i - 1] - 1] = au[i - 1];
        A[(irow[i - 1] - 1) * *n + j] = au[nzs[2] + i - 1];
      }
    }
  } else if (*inputformat == 3) {
    /* general asymmetric matrix */
    for (ITG j = 0; j < *n; j++) {
      for (ITG i = jq[j]; i < jq[j + 1]; i++) {
        A[j * *n + irow[i - 1] - 1] = au[i - 1];
      }
    }
  } else {
    printf(" *ERROR in getcond.c (_densify): unsupported matrix format "
           "%" ITGFORMAT "\n\n",
           *inputformat);
    FORTRAN(stop, ());
  }
}
