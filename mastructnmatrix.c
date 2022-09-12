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

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

void mastructnmatrix(ITG *icols,ITG *jqs,ITG **mast1p,ITG **irowsp,
                  ITG *ipointer,ITG *nzss,ITG *nactive,ITG *nnlconst){

  /* determines the structure of the (N^T*N)-Matrix;
     (i.e. the location of the nonzeros */

  ITG i,j,ii,jj,kk,index,jdof2,jdof1,nmast,ifree,kflag,indexe,isize,
      *mast1=NULL,*irows=NULL,*next=NULL,jstart;
      
  /* the indices in the comments follow FORTRAN convention, i.e. the
     fields start with 1 */

  mast1=*mast1p;
  irows=*irowsp;
  ifree=0;
  kflag=2;

  NNEW(next,ITG,*nzss);

  for(kk=0;kk<*nnlconst;kk++){
     jdof1=kk+1;
        
     for(jj=0;jj<*nactive;jj++){
	jdof2=jj+1;	
	if(jdof2>jdof1){
	   insert(ipointer,&mast1,&next,&jdof1,&jdof2,&ifree,nzss); 
	}
     }
  }
  
  /*   determination of the following fields:       
       
       - irows: row numbers, column per column
       - jqs(i)= location in field irows of the first SUBdiagonal
       nonzero in column i  */
  
  RENEW(irows,ITG,ifree);
  nmast=0;
  jqs[0]=1;
  for(i=0;i<*nactive;i++){
      index=ipointer[i];
      do{
	  if(index==0) break;
	  irows[nmast++]=mast1[index-1];
	  index=next[index-1];
      }while(1);
      jqs[i+1]=nmast+1;
  }
  
  /* sorting the row numbers within each column */
  
  for(i=0;i<*nactive;++i){
      if(jqs[i+1]-jqs[i]>0){
	  isize=jqs[i+1]-jqs[i];
	  FORTRAN(isortii,(&irows[jqs[i]-1],&mast1[jqs[i]-1],&isize,&kflag));
      }
  }
  
  /* removing duplicate entries */
  
  nmast=0;
  for(i=0;i<*nactive;i++){
      jstart=nmast+1;
      if(jqs[i+1]-jqs[i]>0){
	  irows[nmast++]=irows[jqs[i]-1];
	  for(j=jqs[i];j<jqs[i+1]-1;j++){
	      if(irows[j]==irows[nmast-1])continue;
	      irows[nmast++]=irows[j];
	  }
      }
      jqs[i]=jstart;
  }
  jqs[*nactive]=nmast+1;
  
  for(i=0;i<*nactive;i++){
      icols[i]=jqs[i+1]-jqs[i];
  }
  *nzss=jqs[*nactive]-1;
  
  SFREE(next);
  
  *mast1p=mast1;
  *irowsp=irows;
  
  return;
  
}
