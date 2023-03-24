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

void mastructfilter(ITG *icol,ITG *jq,ITG **mastp,ITG **irowp,
                  ITG *ipointer,ITG *nzs,ITG *ndesi,ITG *nodedesi,
		  double *xo,double *yo,double *zo,double *x,
		  double *y,double *z,ITG *nx,ITG *ny,ITG *nz,
		  char *objectset,double *filterrad){
		  

  /* determines the structure of the filter matrix */

  ITG i,j,jj,k,index,idof2,idof1,nmast,ifree,kflag,isize,
      *mast=NULL,*irow=NULL,*next=NULL,jstart,inode1,inode2,
      istart,ipos,*neighbor=NULL,nnodesinside;
  
  double *r=NULL;
  
  /* the indices in the comments follow FORTRAN convention, i.e. the
     fields start with 1 */

  mast=*mastp;
  irow=*irowp;
  ifree=0;
  kflag=2;

  NNEW(next,ITG,*nzs);
  NNEW(r,double,*ndesi);
  NNEW(neighbor,ITG,*ndesi);

  /* loop over all columns */
  
  for(k=0;k<*ndesi;k++){
     idof1=k+1;
     //inode1=nodedesi[k];  

     FORTRAN(near3d_se,(xo,yo,zo,x,y,z,nx,ny,nz,&xo[k],&yo[k],&zo[k],
                        ndesi,neighbor,r,&nnodesinside,filterrad));
     
     /* loop over all rows */
     
     for(j=0;j<nnodesinside;j++){
	//inode2=nodedesi[neighbor[j]-1];
	idof2=(neighbor[j]-1)+1;

	insert_cmatrix(ipointer,&mast,&next,&idof1,&idof2,&ifree,nzs); 
     }
  }
  
  /*   determination of the following fields:       
       
       - irow: row numbers, column per column
       - jq(i)= location in field irow of the first nonzero in column i*/
        
  RENEW(irow,ITG,ifree);
  nmast=0;
  jq[0]=1;
  for(i=0;i<*ndesi;i++){
      index=ipointer[i];
      do{
	  if(index==0) break;
	  irow[nmast++]=mast[index-1];
	  index=next[index-1];
      }while(1);
      jq[i+1]=nmast+1;
  }
  
  /* sorting the row numbers within each column */
  
  for(i=0;i<*ndesi;++i){
      if(jq[i+1]-jq[i]>0){
	  isize=jq[i+1]-jq[i];
	  FORTRAN(isortii,(&irow[jq[i]-1],&mast[jq[i]-1],&isize,&kflag));
      }
  }
  
  /* removing duplicate entries */
  
  nmast=0;
  for(i=0;i<*ndesi;i++){
      jstart=nmast+1;
      if(jq[i+1]-jq[i]>0){
	  irow[nmast++]=irow[jq[i]-1];
	  for(j=jq[i];j<jq[i+1]-1;j++){
	      if(irow[j]==irow[nmast-1])continue;
	      irow[nmast++]=irow[j];
	  }
      }
      jq[i]=jstart;
  }
  jq[*ndesi]=nmast+1;
  
  for(i=0;i<*ndesi;i++){
      icol[i]=jq[i+1]-jq[i];
  }
  *nzs=jq[*ndesi]-1;
  
  SFREE(next);SFREE(neighbor);SFREE(r);
  
  *mastp=mast;
  *irowp=irow;
  
  return;
  
}
