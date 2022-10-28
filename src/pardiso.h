/*     CALCULIX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998 Guido Dhondt                          */
/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation; either version 2 of    */
/*     the License, or (at your option) any later version.               */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

void pardiso_main(double *ad,double *au,double *adb,double *aub,
         double *sigma,double *b,ITG *icol,ITG *irow,
	 ITG *neq,ITG *nzs,ITG *symmetryflag,ITG *inputformat,ITG *jq,
	 ITG *nzs3,ITG *nrhs);

void pardiso_factor(double *ad,double *au,double *adb,double *aub,
                double *sigma,ITG *icol,ITG *irow,
		ITG *neq,ITG *nzs,ITG *symmetryflag,ITG *inputformat,
		ITG *jq,ITG *nzs3);

void pardiso_solve(double *b,ITG *neq,ITG *symmetryflag,ITG *inputformat,
		   ITG *nrhs);

void pardiso_cleanup(ITG *neq,ITG *symmetryflag,ITG *inputformat);

void pardiso_main_cp(double *ad,double *au,double *adb,double *aub,
		     double *sigma,double *b,ITG *icol,ITG *irow,
		     ITG *neq,ITG *nzs,ITG *symmetryflag,ITG *inputformat,
		     ITG *jq,ITG *nzs3,ITG *nrhs);

void pardiso_factor_cp(double *ad,double *au,double *adb,double *aub,
		       double *sigma,ITG *icol,ITG *irow,
		       ITG *neq,ITG *nzs,ITG *symmetryflag,ITG *inputformat,
		       ITG *jq,ITG *nzs3,ITG *iexpl);

void pardiso_solve_cp(double *b,ITG *neq,ITG *symmetryflag,ITG *inputformat,
		   ITG *nrhs);

void pardiso_cleanup_cp(ITG *neq,ITG *symmetryflag,ITG *inputformat);

void FORTRAN(pardiso,(long long *pt,ITG *maxfct,ITG *mnum,ITG *mtype,ITG *phase,
                   ITG *neq,double *aupardiso,ITG *pointers,ITG *irowpardiso,
                   ITG *perm,ITG *nrhs,ITG *iparm,ITG *msglvl,double *b,
                   double *x,ITG *error));

//char envMKL[32];
