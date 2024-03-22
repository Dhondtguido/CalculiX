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

void mastructmm(ITG *icol,ITG *jq,ITG **mastp,ITG **irowp,
                  ITG *ipointer,ITG *nzs,ITG *ndesi,ITG *nodedesi,
		  ITG *iponoelfa,ITG *inoelfa,ITG *nk,char *lakonfa,
		  ITG *konfa,ITG *ipkonfa,ITG *nodedesiinv,
		  ITG *nodedesipos){
		  
  /* determines the structure of the implicit filter matrix */

  ITG i,j,jj,k,index,idof2,idof1,nmast,ifree,kflag,isize,
      *mast=NULL,*irow=NULL,*next=NULL,jstart,inode1,inode2,
      istart,nope,ielem,nopesurf[9],indexe,inode2pos;

  /* the indices in the comments follow FORTRAN convention, i.e. the
     fields start with 1 */

  mast=*mastp;
  irow=*irowp;
  ifree=0;
  kflag=2;
  NNEW(next,ITG,*nzs);
  
  /* lists which external faces correspond to a given node i
     iponoelfa(i) points to an entry j in field inoelfa where:
     inoelfa(1,j): face number as catalogued in fields konfa, lakonfa
     inoelfa(2,j): local node number in the topology description
     inoelfa(3,j): pointer to the next face to which i belongs, or, if
		     none is left: zero */
  
  /* loop over all columns (number of design variables) */
  
  for(k=0;k<*ndesi;k++){
     
     inode1=nodedesi[k]; 
     index=iponoelfa[inode1-1];
     if(index==0)continue;
     
     idof1=k+1;
     
     /* loop over all surfaces belonging to node k */
     
     do{
       
       ielem=inoelfa[3*(index-1)];             
       index=inoelfa[3*(index-1)+3-1];
       
       if (strcmp1(&lakonfa[8*(ielem-1)+1],"4")==0){nope=4;}
       else if (strcmp1(&lakonfa[8*(ielem-1)+1],"8")==0){nope=8;}
       else if (strcmp1(&lakonfa[8*(ielem-1)+1],"3")==0){nope=3;}
       else if (strcmp1(&lakonfa[8*(ielem-1)+1],"6")==0){nope=6;}
       
       indexe=ipkonfa[ielem-1];
       for(jj=0;jj<nope;jj++){
          nopesurf[jj]=konfa[indexe+jj];
       }

       /* loop over all nodes on surface */
     
       for(j=0;j<nope;j++){
         
	 inode2=nopesurf[j];
	 inode2pos=nodedesipos[inode2-1]-1; 
         
	 if(nodedesiinv[inode2-1]==0)continue;
	 if(inode2pos>k){
	 
	   idof2=inode2pos+1;
           insert(ipointer,&mast,&next,&idof1,&idof2,&ifree,nzs); 
	 
         }
       }
       
     }while(index!=0);            
     
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
  
  SFREE(next);
  
  *mastp=mast;
  *irowp=irow;
  
  return;
  
}
