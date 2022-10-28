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

#ifndef __CCX_SPOOLES_H
#define __CCX_SPOOLES_H

/*
 * separated from CalculiX.h: otherwise everyone would have to include
 * the spooles header files
 */

#include <pthread.h>
#include <misc.h>
#include <FrontMtx.h>
#include <SymbFac.h>
#if USE_MT
#include <MT/spoolesMT.h>
#endif

/* increase this for debugging */
#define DEBUG_LVL	0

struct factorinfo 
{
	ITG size;
	double cpus[11];
	IV *newToOldIV, *oldToNewIV;
	SolveMap *solvemap;
	FrontMtx *frontmtx;
	SubMtxManager *mtxmanager;
	ETree *frontETree;
	ITG nthread;
	FILE *msgFile;

};

void spooles_factor(double *ad, double *au, double *adb, double *aub, 
                    double *sigma, ITG *icol, ITG *irow,
                    ITG *neq, ITG *nzs, ITG *symmetryflag,
                    ITG *inputformat, ITG *nzs3);

void spooles_solve(double *b, ITG *neq);

void spooles_cleanup();

void spooles_factor_rad(double *ad, double *au, double *adb, double *aub, 
			double *sigma, ITG *icol, ITG *irow,
			ITG *neq, ITG *nzs, ITG *symmetryflag,
			ITG *inputformat,ITG *iexpl);

void spooles_solve_rad(double *b, ITG *neq);

void spooles_cleanup_rad();

#endif
