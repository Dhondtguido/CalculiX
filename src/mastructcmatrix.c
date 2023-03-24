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

void mastructcmatrix(ITG *icolc,ITG *jqc,ITG **mast1p,ITG **irowcp,
                  ITG *ipointer,ITG *nzsc,ITG *ndesibou,ITG *nodedesibou,
		  ITG *nodedesiinvbou,ITG *jqs,ITG *irows,ITG *icols,
		  ITG *ndesi,ITG *nodedesi){
		  

  /* determines the structure of the D-Matrix which contains the 
     covariance information of the boundary nodes;
     (i.e. the location of the nonzeros */

  ITG i,j,jj,kk,index,jdof2,jdof1,nmast,ifree,kflag,isize,
      *mast1=NULL,*irowc=NULL,*next=NULL,jstart,inode1,inode2,
      inode3,istart,irow,ipos;
      
  /* the indices in the comments follow FORTRAN convention, i.e. the
     fields start with 1 */

  mast1=*mast1p;
  irowc=*irowcp;
  ifree=0;
  kflag=2;

  NNEW(next,ITG,*nzsc);

  for(kk=0;kk<*ndesi;kk++){
     jdof1=kk+1;
     inode1=nodedesi[kk];  
     
     for(jj=0;jj<*ndesi;jj++){
        inode2=nodedesi[jj];
        if(nodedesiinvbou[inode2-1]>0){
	   FORTRAN(nident,(nodedesibou,&inode2,ndesibou,&irow));
	   jdof2=irow;   
	   ipos=jj+1; 
           /* diagonal matrix entries */
	   if(inode1==inode2){
	      insert_cmatrix(ipointer,&mast1,&next,&jdof1,&jdof2,
	         &ifree,nzsc); 
           /* lower triangel matrix entries */
           }else if(inode1<inode2){
	      istart=jqs[kk]-1; 
	      FORTRAN(nident,(&irows[istart],&ipos,&icols[kk],&irow));
	      inode3=nodedesi[irows[istart+irow-1]-1];
	      if(inode3==inode2){
	         insert_cmatrix(ipointer,&mast1,&next,&jdof1,&jdof2,
	            &ifree,nzsc); 
	      }
 	   /* upper triangel matrix entries */
	   }else{
	      istart=jqs[jj]-1;  
	      FORTRAN(nident,(&irows[istart],&jdof1,&icols[jj],&irow));
	      inode3=nodedesi[irows[istart+irow-1]-1];
	      if(inode3==inode1){
	         insert_cmatrix(ipointer,&mast1,&next,&jdof1,&jdof2,
	            &ifree,nzsc); 
	      }
     	   }
	}
     }
  }
  
  /*   determination of the following fields:       
       
       - irowc: row numbers, column per column
       - jqc(i)= location in field irowc of the first nonzero in column i*/
       
  
  RENEW(irowc,ITG,ifree);
  nmast=0;
  jqc[0]=1;
  for(i=0;i<*ndesi;i++){
      index=ipointer[i];
      do{
	  if(index==0) break;
	  irowc[nmast++]=mast1[index-1];
	  index=next[index-1];
      }while(1);
      jqc[i+1]=nmast+1;
  }
  
  /* sorting the row numbers within each column */
  
  for(i=0;i<*ndesi;++i){
      if(jqc[i+1]-jqc[i]>0){
	  isize=jqc[i+1]-jqc[i];
	  FORTRAN(isortii,(&irowc[jqc[i]-1],&mast1[jqc[i]-1],&isize,&kflag));
      }
  }
  
  /* removing duplicate entries */
  
  nmast=0;
  for(i=0;i<*ndesi;i++){
      jstart=nmast+1;
      if(jqc[i+1]-jqc[i]>0){
	  irowc[nmast++]=irowc[jqc[i]-1];
	  for(j=jqc[i];j<jqc[i+1]-1;j++){
	      if(irowc[j]==irowc[nmast-1])continue;
	      irowc[nmast++]=irowc[j];
	  }
      }
      jqc[i]=jstart;
  }
  jqc[*ndesi]=nmast+1;
  
  for(i=0;i<*ndesi;i++){
      icolc[i]=jqc[i+1]-jqc[i];
  }
  *nzsc=jqc[*ndesi]-1;
  
  SFREE(next);
  
  *mast1p=mast1;
  *irowcp=irowc;
  
  return;
  
}
