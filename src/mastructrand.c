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

void mastructrand(ITG *icols,ITG *jqs,ITG **mast1p,ITG **irowsp,ITG *ipointer,
                  ITG *nzss,ITG *ndesi,double *corrlen,double *xo,double *yo,
                  double *zo,double *x,double *y,double *z,ITG *nx,ITG *ny,
                  ITG *nz){

  /* determines the structure of the covariance matrix;
     (i.e. the location of the nonzeros */

  ITG i,j,ii,index,jdof2,jdof1,nmast,ifree,kflag,isize,*mast1=NULL,
      *irows=NULL,*next=NULL,jstart,idesvar,*neighbor=NULL,nnodesinside;

  double *r=NULL,searchlen;

  /* the indices in the comments follow FORTRAN convention, i.e. the
     fields start with 1 */

  mast1=*mast1p;
  irows=*irowsp;
  ifree=0;
  kflag=2;
  searchlen=5**corrlen;

  NNEW(next,ITG,*nzss);
  NNEW(neighbor,ITG,*ndesi+6);
  NNEW(r,double,*ndesi+6);

  for(idesvar=0;idesvar<*ndesi;idesvar++){
      jdof1=idesvar+1;
      
      /* nodes within 3 times the correlation length */
      
      FORTRAN(near3d_se,(xo,yo,zo,x,y,z,nx,ny,nz,&xo[idesvar],
			 &yo[idesvar],&zo[idesvar],ndesi,neighbor,
			 r,&nnodesinside,&searchlen));
      
      for(ii=0;ii<nnodesinside;ii++){
	  jdof2=neighbor[ii];
	  insert(ipointer,&mast1,&next,&jdof1,&jdof2,&ifree,nzss);
      }
  }
  
  /*   determination of the following fields:       
       
       - irows: row numbers, column per column
       - jqs(i)= location in field irows of the first SUBdiagonal
       nonzero in column i  */
  
  RENEW(irows,ITG,ifree);
  nmast=0;
  jqs[0]=1;
  for(i=0;i<*ndesi;i++){
      index=ipointer[i];
      do{
	  if(index==0) break;
	  irows[nmast++]=mast1[index-1];
	  index=next[index-1];
      }while(1);
      jqs[i+1]=nmast+1;
  }
  
  /* sorting the row numbers within each column */
  
  for(i=0;i<*ndesi;++i){
      if(jqs[i+1]-jqs[i]>0){
	  isize=jqs[i+1]-jqs[i];
	  FORTRAN(isortii,(&irows[jqs[i]-1],&mast1[jqs[i]-1],&isize,&kflag));
      }
  }
  
  /* removing duplicate entries */
  
  nmast=0;
  for(i=0;i<*ndesi;i++){
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
  jqs[*ndesi]=nmast+1;
  
  for(i=0;i<*ndesi;i++){
      icols[i]=jqs[i+1]-jqs[i];
  }
  *nzss=jqs[*ndesi]-1;
  
  SFREE(next);SFREE(neighbor);SFREE(r);
  
  *mast1p=mast1;
  *irowsp=irows;
  
  return;
  
}
