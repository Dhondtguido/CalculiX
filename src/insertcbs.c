/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2011 Guido Dhondt                          */

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

void insertcbs(ITG *ipointer, ITG **irowp, ITG **nextp, ITG *i1,
	    ITG *i2, ITG *ifree, ITG *nzs_){

  /*   inserts a new nonzero matrix position into the data structure 
       in FORTRAN notation: 
       - ipointer(i) points to a position in field irow containing
         the row number of a nonzero position in column i; 
         next(ipointer(i)) points a position in field irow containing
         the row number of another nonzero position in column i, and
         so on until no nonzero positions in column i are left; for 
         the position j in field irow containing the momentarily last
         nonzero number in column i we have next(j)=0 

       notice that in C the positions start at 0 and not at 1 as in 
       FORTRAN; the present routine is written in FORTRAN convention 

       routine specifically for the CBS-Method in CFD */
  
  ITG idof1,idof2,istart,*irow=NULL,*next=NULL;

  irow=*irowp;
  next=*nextp;

  if(*i1<*i2){
    idof1=*i1;
    idof2=*i2;
  }
  else{
    idof1=*i2;
    idof2=*i1;
  }

  if(ipointer[idof2-1]==0){
    ++*ifree;
    if(*ifree>*nzs_){
      *nzs_=(ITG)(1.1**nzs_);
      RENEW(irow,ITG,*nzs_);
      RENEW(next,ITG,*nzs_);
    }
    ipointer[idof2-1]=*ifree;
    irow[*ifree-1]=idof1;
    next[*ifree-1]=0;
  }
  else{
    istart=ipointer[idof2-1];
    while(1){
      if(irow[istart-1]==idof1) break;
      if(next[istart-1]==0){
	++*ifree;
	if(*ifree>*nzs_){
	  *nzs_=(ITG)(1.1**nzs_);
	  RENEW(irow,ITG,*nzs_);
	  RENEW(next,ITG,*nzs_);
	}
	next[istart-1]=*ifree;
	irow[*ifree-1]=idof1;
	next[*ifree-1]=0;
	break;
      }
      else{
	istart=next[istart-1];
      }
    }
  }

  *irowp=irow;
  *nextp=next;
  
  return;

}
