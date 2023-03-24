/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2015 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                                                                       */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "CalculiX.h"

void insert(ITG *ipointer, ITG **mast1p, ITG **nextp, ITG *i1,
	    ITG *i2, ITG *ifree, ITG *nzs_){

  /* routine for the lower triangular matrix, excluding the diagonal */

  /*   inserts a new nonzero matrix position into the data structure 
       in FORTRAN notation: 
       - ipointer(i) points to a position in field mast1 containing
         the row number of a nonzero position in column i; 
         next(ipointer(i)) points a position in field mast1 containing
         the row number of another nonzero position in column i, and
         so on until no nonzero positions in column i are left; for 
         the position j in field mast1 containing the momentarily last
         nonzero number in column i we have next(j)=0 

       notice that in C the positions start at 0 and not at 1 as in 
       FORTRAN; the present routine is written in FORTRAN convention */

    ITG idof1,idof2,istart,*mast1=NULL,*next=NULL;

  mast1=*mast1p;
  next=*nextp;

  if(*i1==*i2) return;
  if(*i1<*i2){
    idof1=*i2;
    idof2=*i1-1;
  }
  else{
    idof1=*i1;
    idof2=*i2-1;
  }

  if(*ifree>=*nzs_){
      *nzs_=(ITG)(1.1**nzs_);
      RENEW(mast1,ITG,*nzs_);
      RENEW(next,ITG,*nzs_);
  }
  mast1[*ifree]=idof1;
  next[*ifree]=ipointer[idof2];
  ipointer[idof2]=++*ifree;

  *mast1p=mast1;
  *nextp=next;
  
  return;

}
