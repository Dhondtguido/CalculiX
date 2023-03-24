/*     CalculiX - A 3-dimensional finite element program                   */
/*              Copyright (C) 1998-2005 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CalculiX.h"

void convert2rowbyrow(double *ad,double *au, ITG *icol,ITG *irow, 
		      ITG *jq,ITG *neq,ITG *nzs,double **aupardisop,
		      ITG **pointersp,ITG **icolpardisop){
  
  ITG i,j,k,l,*irowpardiso=NULL,kflag,kstart,n,ndim,*icolpardiso=NULL,
    *pointers=NULL;

  double *aupardiso=NULL;
	  
  /* lower triangular matrix is stored column by column in
     au; the diagonal terms are stored in ad; the row numbers of the
     subdiagonal terms are stored in irow, column by column (and
     ordered in ascending order for each column); the first entry
     in column i is stored in jq(i), the number of entries for
     column i in icol(i) */

  /* the matrix is stored as full matrix, row by row, in aupardiso
     (needed for shock smoothing), the corresponding row numbers
     are stored in pointers and the start of row i is stored in
     icolpardiso(i) (FORTRAN convention); */

  icolpardiso=*icolpardisop;aupardiso=*aupardisop;pointers=*pointersp;
	
  /* reordering lower triangular matrix */

  ndim=*nzs;
  NNEW(pointers,ITG,*neq+1);
  NNEW(irowpardiso,ITG,ndim);
  NNEW(icolpardiso,ITG,ndim);
  NNEW(aupardiso,double,ndim);
	  
  k=0;
  for(i=0;i<*neq;i++){
    for(j=0;j<icol[i];j++){
      icolpardiso[k]=i+1;
      irowpardiso[k]=irow[k];
      aupardiso[k]=au[k];
      k++;
    }
  }
	  
  /* pardiso needs the entries row per row; so sorting is
     needed */
	  
  kflag=2;
  FORTRAN(isortiid,(irowpardiso,icolpardiso,aupardiso,
		    &ndim,&kflag));
	  
  /* sorting each row */
	  
  k=0;
  pointers[0]=1;
  if(ndim>0){
    for(i=0;i<*neq;i++){
      j=i+1;
      kstart=k;
      do{

	/* end of row reached */

	if(irowpardiso[k]!=j){
	  n=k-kstart;
	  FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
			   &n,&kflag));
	  pointers[i+1]=k+1;
	  break;
	}else{

	  /* end of last row reached */

	  if(k+1==ndim){
	    n=k-kstart+1;
	    FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
			     &n,&kflag));
	    break;
	  }else{

	    /* end of row not yet reached */

	    k++;
	  }
	}
      }while(1);
    }
  }
  pointers[*neq]=ndim+1;
  SFREE(irowpardiso);

  /* composing the matrix: lower triangle + diagonal + upper triangle */

  ndim=*neq+2**nzs;
  RENEW(icolpardiso,ITG,ndim);
  RENEW(aupardiso,double,ndim);
  k=ndim;
  for(i=*neq-1;i>=0;i--){
    l=k+1;
    for(j=jq[i+1]-1;j>=jq[i];j--){
      icolpardiso[--k]=irow[j-1];
      aupardiso[k]=au[j-1];
    }
    icolpardiso[--k]=i+1;
    aupardiso[k]=ad[i];
    for(j=pointers[i+1]-1;j>=pointers[i];j--){
      icolpardiso[--k]=icolpardiso[j-1];
      aupardiso[k]=aupardiso[j-1];
    }
    pointers[i+1]=l;
  }
  pointers[0]=1;

  *icolpardisop=icolpardiso;*aupardisop=aupardiso;*pointersp=pointers;

  return;
}

