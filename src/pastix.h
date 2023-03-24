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

void pastix_main(double *ad, double *au, double *adb, double *aub, 
         double *sigma,double *b, ITG *icol, ITG *irow, 
	 ITG *neq, ITG *nzs,ITG *symmetryflag,ITG *inputformat,ITG *jq,
	 ITG *nzs3,ITG *nrhs);

void pastix_main_as(double *ad, double *au, double *adb, double *aub, 
     double *sigma,double *b, ITG *icol, ITG *irow, 
	 ITG *neq, ITG *nzs,ITG *symmetryflag,ITG *inputformat,
	 ITG *jq, ITG *nzs3,ITG *nrhs);
	 
void pastix_init(double *ad, double *au, double *adb, double *aub, 
                double *sigma,ITG *icol, ITG *irow, 
		ITG *neq, ITG *nzs,ITG *symmetryflag,ITG *inputformat,
		ITG *jq,ITG *nzs3);

void pastix_csc_conversion(double *ad, double *au, double *adb, double *aub, 
                double *sigma,ITG *icol, ITG *irow, 
		ITG *neq, ITG *nzs,ITG *symmetryflag,ITG *inputformat,
		ITG *jq,ITG *nzs3);

void pastix_analyze();

void pastix_factor_main(double *ad, double *au, double *adb, double *aub, 
        double *sigma,ITG *icol, ITG *irow, 
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3);

void pastix_factor_main_as(double *ad, double *au, double *adb, double *aub, 
        double *sigma,ITG *icol, ITG *irow, 
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3);

void pastix_factor_main_generic(double *ad, double *au, double *adb, double *aub, 
        double *sigma,ITG *icol, ITG *irow, 
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3);

void pastix_factor(double *ad, double *au, double *adb, double *aub, 
                double *sigma,ITG *icol, ITG *irow, 
		ITG *neq, ITG *nzs,ITG *symmetryflag,ITG *inputformat,
		ITG *jq,ITG *nzs3);


ITG pastix_solve(double *b,ITG *neq,ITG *symmetryflag,ITG *nrhs);

ITG pastix_solve_as(double *x, ITG *neq,ITG *symmetryflag,ITG *nrhs);

ITG pastix_solve_generic(double *x, ITG *neq,ITG *symmetryflag,ITG *nrhs);

void pastix_cleanup(ITG *neq,ITG *symmetryflag);

void pastix_set_globals(char mode);

