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
#include <stdio.h>


/**
 * getcond
 * =======
 *
 * Computes the condition number of a small sparse square matrix.
 *
 * Note
 * ----
 *
 * The sparse matrix must be small enough to be stored temporarily as dense
 * matrix in memory. The condition number is computed via the singular value
 * decomposition provided by the LAPACK routine DGESDD.
 *
 * Usage
 * -----
 * cond = getcond(&neq[0], ad, au, irow, jq, nzs, &inputformat, &ierrcond);
 *
 * Authors
 * -------
 * 2025-09-26 Christoph Woelfle - initial commit
 *
 * Parameters
 * ----------
 * @param n row and column dimension of input matrix
 * @param ad diagonal elements of input matrix
 * @param au offdiagonal elements of input matrix
 * @param irow row indices of elements in au
 * @param jq index in irow of first nonzero of each column
 * @param nzs number of nonzeroes
 * @param inputformat flag for input matrix structure
 * @param ierr error return code
 * ...
 */
double getcond(const ITG *n, const double *ad, const double *au,
               const ITG *irow, const ITG *jq, const ITG *nzs,
               const ITG *inputformat, ITG *ierr) {
  double *A, *s, *work, cond;
  ITG *iwork;

  /* store as dense matrix */
  NNEW(A, double, *n **n);

  densify(A, n, ad, au, irow, jq, nzs, inputformat);

  /* compute condition number */
  char jobz[2] = "N";
  ITG lwork = 10 * *n;
  NNEW(s, double, *n);
  NNEW(iwork, ITG, 8 * *n);
  MNEW(work, double, lwork);
  ITG info = 0;

  /* query optimal size of lwork  */
  lwork = -1;
  FORTRAN(dgesdd,
          (jobz, n, n, A, n, s, NULL, n, NULL, n, work, &lwork, iwork, &info));
  lwork = work[0];
  RENEW(work, double, lwork);

  printf(" *INFO getcond.c: calling DGESDD\n");
  FORTRAN(dgesdd,
          (jobz, n, n, A, n, s, NULL, n, NULL, n, work, &lwork, iwork, &info));
  if (info != 0) {
    printf(" *ERROR in getcond.c: DGESDD returned code %" ITGFORMAT "\n\n",
           info);
    *ierr = info;
    FORTRAN(stop, ());
  }
  cond = s[0] / s[*n - 1];

  /* free memory */
  SFREE(A);
  SFREE(s);
  SFREE(iwork);
  SFREE(work);

  /* return condition number */
  *ierr = 0;
  return cond;
}
