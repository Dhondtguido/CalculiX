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

void mastructdmatrix(ITG *icold,ITG *jqd,ITG **mast1p,ITG **irowdp,
                  ITG *ipointer,ITG *nzss,ITG *ndesibou,ITG *nodedesibou,
		  ITG *nodedesiinvbou,ITG *jqs,ITG *irows,ITG *icols,
		  ITG *ndesi,ITG *nodedesi){
		  

  /* determines the structure of the D-Matrix which contains the 
     covariance information of the boundary nodes;
     (i.e. the location of the nonzeros */

  ITG i,j,jj,kk,index,jdof2,jdof1,nmast,ifree,kflag,isize,
      *mast1=NULL,*irowd=NULL,*next=NULL,jstart,inode1,inode2,ipos1,
      ipos2,istart,idesvar;
      
  /* the indices in the comments follow FORTRAN convention, i.e. the
     fields start with 1 */

  mast1=*mast1p;
  irowd=*irowdp;
  ifree=0;
  kflag=2;

  NNEW(next,ITG,*nzss);

  for(kk=0;kk<*ndesibou;kk++){
     jdof1=kk+1;
     inode1=nodedesibou[kk];
     FORTRAN(nident,(nodedesi,&inode1,ndesi,&ipos1));     
     istart=jqs[ipos1-1]-1;   
     
     for(jj=0;jj<icols[ipos1-1];jj++){
     	idesvar=irows[istart+jj]-1;
	inode2=nodedesi[idesvar];
     	if(nodedesiinvbou[inode2-1]>0){
     	   FORTRAN(nident,(nodedesibou,&inode2,ndesibou,&ipos2));
	   jdof2=ipos2;    
     	   if(jdof2>jdof1){ 	   
     	      insert(ipointer,&mast1,&next,&jdof1,&jdof2,&ifree,nzss); 
     	   }
     	}
     }
  }
  
  /*   determination of the following fields:       
       
       - irowd: row numbers, column per column
       - jqd(i)= location in field irowd of the first SUBdiagonal
       nonzero in column i  */
  
  RENEW(irowd,ITG,ifree);
  nmast=0;
  jqd[0]=1;
  for(i=0;i<*ndesibou;i++){
      index=ipointer[i];
      do{
	  if(index==0) break;
	  irowd[nmast++]=mast1[index-1];
	  index=next[index-1];
      }while(1);
      jqd[i+1]=nmast+1;
  }
  
  /* sorting the row numbers within each column */
  
  for(i=0;i<*ndesibou;++i){
      if(jqd[i+1]-jqd[i]>0){
	  isize=jqd[i+1]-jqd[i];
	  FORTRAN(isortii,(&irowd[jqd[i]-1],&mast1[jqd[i]-1],&isize,&kflag));
      }
  }
  
  /* removing duplicate entries */
  
  nmast=0;
  for(i=0;i<*ndesibou;i++){
      jstart=nmast+1;
      if(jqd[i+1]-jqd[i]>0){
	  irowd[nmast++]=irowd[jqd[i]-1];
	  for(j=jqd[i];j<jqd[i+1]-1;j++){
	      if(irowd[j]==irowd[nmast-1])continue;
	      irowd[nmast++]=irowd[j];
	  }
      }
      jqd[i]=jstart;
  }
  jqd[*ndesibou]=nmast+1;
  
  for(i=0;i<*ndesibou;i++){
      icold[i]=jqd[i+1]-jqd[i];
  }
  *nzss=jqd[*ndesibou]-1;
  
  SFREE(next);
  
  *mast1p=mast1;
  *irowdp=irowd;
  
  return;
  
}
