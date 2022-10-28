/*     CalculiX - A 3-dimensional finite element program                   */
/*              Copyright (C) 1998-2015 Guido Dhondt                          */

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

void preiter(double *ad, double **aup, double *b, ITG **icolp, ITG **irowp, 
	     ITG *neq, ITG *nzs, ITG *isolver, ITG *iperturb){
  
  ITG precFlg,niter=5000000,ndim,i,j,k,ier,*icol=NULL,*irow=NULL,
    *irow_save=NULL,*icol_save=NULL;
  double eps=1.e-4,*u=NULL,*au=NULL;

  if(*neq==0) return;

  /* icol(i) = # subdiagonal nonzeros in column i (i=1,neq)
     irow(i) = row number of entry i in au (i=1,nzs)
     ad(i) = diagonal term in column i of the matrix
     au(i) = subdiagonal nonzero term i; the terms are entered
             column per column */

  au=*aup;
  irow=*irowp;
  icol=*icolp;

  if(*iperturb>1){
    NNEW(irow_save,ITG,*nzs);
    NNEW(icol_save,ITG,*neq);
    for(i=0;i<*nzs;++i){
      irow_save[i]=irow[i];
    }
    for(i=0;i<*neq;++i){
      icol_save[i]=icol[i];
    }
  }

  if(*isolver==2) {precFlg=0;}
  else {precFlg=3;}

  ndim=*neq+*nzs;

  RENEW(au,double,ndim);
  RENEW(irow,ITG,ndim);
  RENEW(icol,ITG,ndim);

  k=*nzs;
  for(i=*neq-1;i>=0;--i){
    for(j=0;j<icol[i];++j){
      icol[--k]=i+1;
    }
  }

  k=*nzs;
  j=0;
  for(i=0;i<*neq;++i){
    au[k]=ad[i];
    irow[k]=++j;
    icol[k]=j;
    ++k;
  }

  /* rearranging the storage of the left hand side */

  FORTRAN(rearrange,(au,irow,icol,&ndim,neq));

  RENEW(irow,ITG,*neq);

  NNEW(u,double,*neq);

  ier=cgsolver(au,u,b,*neq,ndim,icol,irow,&eps,&niter,precFlg);

  printf("error condition (0=good, 1=bad) = %" ITGFORMAT "\n",ier);
  printf("# of iterations = %" ITGFORMAT "\n",niter);

  for(i=0;i<*neq;++i){
    b[i]=u[i];
  }

  SFREE(u);

  if(*iperturb>1){
    RENEW(irow,ITG,*nzs);
    RENEW(icol,ITG,*neq);
    for(i=0;i<*nzs;++i){
      irow[i]=irow_save[i];
    }
    for(i=0;i<*neq;++i){
      icol[i]=icol_save[i];
    }
    SFREE(irow_save);SFREE(icol_save);
  }

  *aup=au;
  *irowp=irow;
  *icolp=icol;

  return;
}
