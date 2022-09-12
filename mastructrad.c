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

void mastructrad(ITG *ntr,ITG *nloadtr,char *sideload,ITG *ipointerrad,
              ITG **mast1radp,ITG **irowradp,ITG *nzsrad,
              ITG *jqrad,ITG *icolrad){

/* determining the structure of the viewfactor and the radiation
   matrix (both have the same structure). Only the structure of the
   lower half of the matrix is determined, since the structure of
   the upper half is identical */

  char crcav[4]="   \0";

  ITG three=3,i,j,k,l,ii,jj,icav,jcav,*mast1rad=NULL,*irowrad=NULL,
      ifree,nzsrad_,nmast,isubtract,isize,kflag,istart,istartold;

  mast1rad=*mast1radp;
  irowrad=*irowradp;

  kflag=2;
  nzsrad_=*nzsrad;

    /* determining the position of each nonzero matrix element in
       the SUPERdiagonal matrix */

  ifree=0;

  for(ii=1;ii<=*ntr;ii++){
      i=nloadtr[ii-1]-1;
      strcpy1(crcav,&sideload[20*i+17],three);
      icav=atoi(crcav);
      for(jj=1;jj<=ii;jj++){
	  j=nloadtr[jj-1]-1;
	  strcpy1(crcav,&sideload[20*j+17],three);
	  jcav=atoi(crcav);
	  if(icav==jcav){
	      insertrad(ipointerrad,&mast1rad,&irowrad,&ii,
                     &jj,&ifree,&nzsrad_);
	  }
      }
  }

  /*   storing the nonzero nodes in the SUPERdiagonal columns:
       mast1rad contains the row numbers (already done),
       irowrad the column numbers (done in the next lines) */
    
  for(i=0;i<*ntr;++i){
      if(ipointerrad[i]==0){
	  printf(" *ERROR in mastructrad: zero column\n");
	  printf("       DOF=%" ITGFORMAT "\n",i);
	  FORTRAN(stop,());
      }
      istart=ipointerrad[i];
      while(1){
	  istartold=istart;
	  istart=irowrad[istart-1];
	  irowrad[istartold-1]=i+1;
	  if(istart==0) break;
      }
  }
  
  nmast=ifree;
  
  /* summary */
  
  printf(" number of radiation equations\n");
  printf(" %" ITGFORMAT "\n",*ntr);
  printf(" number of nonzero radiation matrix elements\n");
  printf(" %" ITGFORMAT "\n",2*nmast-*ntr);
  printf(" \n");

    /* switching from a SUPERdiagonal inventory to a SUBdiagonal one:
       since the nonzeros are located in symmetric positions mast1rad
       can be considered to contain the column numbers and irowrad the
       row numbers; after sorting mast1rad the following results:

       - irowrad contains the row numbers of the SUBdiagonal
         nonzero's, column per column
       - mast1rad contains the column numbers

       Furthermore, the following fields are determined:       

       - icolrad(i)=# SUBdiagonal nonzero's in column i
       - jqrad(i)= location in field irow of the first SUBdiagonal
         nonzero in column i  */

    /* ordering the column numbers in mast1rad */
  
  FORTRAN(isortii,(mast1rad,irowrad,&nmast,&kflag));
  
  /* filtering out the diagonal elements and generating icolrad and jqrad */
  
  isubtract=0;
  for(i=0;i<*ntr;++i){icolrad[i]=0;}
  k=0;
  for(i=0;i<nmast;++i){
      if(mast1rad[i]==irowrad[i]){++isubtract;}
      else{
	  mast1rad[i-isubtract]=mast1rad[i];
	  irowrad[i-isubtract]=irowrad[i];
	  if(k!=mast1rad[i]){
	      for(l=k;l<mast1rad[i];++l){jqrad[l]=i+1-isubtract;}
	      k=mast1rad[i];
	  }
	  ++icolrad[k-1];
      }
  }
  nmast=nmast-isubtract;
  for(l=k;l<*ntr+1;++l){jqrad[l]=nmast+1;}

    /* sorting the row numbers within each column */
  
  for(i=0;i<*ntr;++i){
      if(jqrad[i+1]-jqrad[i]>0){
	  isize=jqrad[i+1]-jqrad[i];
	  FORTRAN(isortii,(&irowrad[jqrad[i]-1],&mast1rad[jqrad[i]-1],
                  &isize,&kflag));
      }
  }
  
  *nzsrad=jqrad[*ntr]-1;
  
  *mast1radp=mast1rad;
  *irowradp=irowrad;
  
  return;

}
	  
