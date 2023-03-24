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

void mastructse(ITG *kon, ITG *ipkon, char *lakon, ITG *ne,
	      ITG *ipompc, ITG *nodempc, ITG *nmpc,
	      ITG *nactdof, ITG *jqs, ITG **mast1p, ITG **irowsp, 
              ITG *ipointer, ITG *nzss, ITG *mi, ITG *mortar,
	      ITG *nodedesi, ITG *ndesi,ITG *icoordinate,ITG *ielorien,
	      ITG *istartdesi,ITG *ialdesi){

  /* determines the structure of the sensitivity matrix;
     (i.e. the location of the nonzeros */

  char lakonl[2]=" \0";

  ITG i,j,k,id,index,jdof1,idof1,idof2,nmast,ifree,kdof1,
    node,ist,kflag,indexe,nope,isize,*mast1=NULL,ii,
    *irows=NULL,mt=mi[1]+1,*next=NULL,jstart,idesvar;

  mast1=*mast1p;
  irows=*irowsp;
  ifree=0;
  kflag=2;

  NNEW(next,ITG,*nzss);

  for(idesvar=0;idesvar<*ndesi;idesvar++){
      idof2=idesvar+1;

      for(ii=istartdesi[idesvar]-1;ii<istartdesi[idesvar+1]-1;ii++){
	  i=ialdesi[ii]-1;
      
	  if(ipkon[i]<0) continue;
	  if(strcmp1(&lakon[8*i],"F")==0)continue;
	  indexe=ipkon[i];
/* Bernhardi start */
	  if (strcmp1(&lakon[8*i+3],"8I")==0)nope=11;
	  else if(strcmp1(&lakon[8*i+3],"20")==0)nope=20;
/* Bernhardi end */
	  else if (strcmp1(&lakon[8*i+3],"8")==0)nope=8;
	  else if (strcmp1(&lakon[8*i+3],"10")==0)nope=10;
	  else if ((strcmp1(&lakon[8*i+3],"4")==0)||
		   (strcmp1(&lakon[8*i+2],"4")==0)) nope=4;
	  else if (strcmp1(&lakon[8*i+3],"15")==0)nope=15;
	  else if (strcmp1(&lakon[8*i+3],"6")==0)nope=6;
	  else if (strcmp1(&lakon[8*i],"E")==0){
	      if((strcmp1(&lakon[8*i+6],"C")==0)&&(*mortar==1)){
		  
		  /* face-to-face contact (all nodes already belong
		     to other elements */
		  
		  continue;
	      }else if(strcmp1(&lakon[8*i+6],"F")!=0){
		  
		  /* node-to-face contact */
		  
		  lakonl[0]=lakon[8*i+7];
		  nope=atoi(lakonl)+1;
	      }else{
		  
		  /* advection elements */
		  
		  continue;
	      }
	  }else continue;
	  
	  /* displacement degrees of freedom */
	  
	  for(j=0;j<nope;++j){
	      node=kon[indexe+j]-1;
	      for(k=1;k<4;k++){
		  idof1=nactdof[mt*node+k];
		  if(idof1>0){
		      insertfreq(ipointer,&mast1,&next,
				 &idof1,&idof2,&ifree,nzss);
		  }else if(*nmpc!=0){
		      if(idof1!=2*(idof1/2)){
			  id=(-idof1+1)/2-1;
			  ist=ipompc[id]-1;
			  index=nodempc[3*ist+2];
			  if(index==0) continue;
			  do{
			      jdof1=nactdof[mt*(nodempc[3*index-3]-1)
					    +nodempc[3*index-2]];
			      if(jdof1>0){
				  insertfreq(ipointer,&mast1,&next,
					     &jdof1,&idof2,&ifree,nzss);
			      }
			      index=nodempc[3*index-1];
			      if(index==0) break;
			  }while(1);
		      }
		  }
	      }
	  }
      }
  }
  
  /* determine irows and jqs */
  
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
  
  *nzss=jqs[*ndesi]-1;
  
  SFREE(next);
  
  *mast1p=mast1;
  *irowsp=irows;

  return;

}
