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

#ifdef PARDISO

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CalculiX.h"
#include "pardiso.h"

/* next 3 lines are for the simulateous use of PARDISO and PaStiX */

#ifdef COMPANY
#include <mkl_service.h>
#endif

ITG *icolpardisocp=NULL,*pointerscp=NULL,iparmcp[64];
long long ptcp[64];
double *aupardisocp=NULL;
/* double dparm[64];  not used */
ITG mthread_mkl_cp=0;
char envMKL_cp[32];

void pardiso_factor_cp(double *ad, double *au, double *adb, double *aub, 
		       double *sigma,ITG *icol, ITG *irow, 
		       ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		       ITG *jq, ITG *nzs3,ITG *iexpl){

  char *env;
  /*  char env1[32]; */
  ITG i,j,k,l,maxfct=1,mnum=1,phase=12,nrhs=1,*perm=NULL,mtype,
    msglvl=0,error=0,*irowpardisocp=NULL,kflag,kstart,n,ifortran,
    lfortran,index,id,k2;
  ITG ndim,nthread,nthread_v;
  double *b=NULL,*x=NULL;

  /* messages only for implicit calculations (for explicit calculations
     the output should be kept to a minimum because of performance
     reasons*/
	
  if(*iexpl<=1){
    if(*symmetryflag==0){
      printf(" Factoring the system of equations using the symmetric pardiso solver\n");
    }else{
      printf(" Factoring the system of equations using the unsymmetric pardiso solver\n");
    }
  }

  iparmcp[0]=0;
  iparmcp[1]=3;
  /* set MKL_NUM_THREADS to min(CCX_NPROC_EQUATION_SOLVER,OMP_NUM_THREADS)
     must be done once  */
  if (mthread_mkl_cp == 0) {
    nthread=1;
    env=getenv("MKL_NUM_THREADS");
    if(env) {
      nthread=atoi(env);}
    else {
      env=getenv("OMP_NUM_THREADS");
      if(env) {nthread=atoi(env);}
    }
    env=getenv("CCX_NPROC_EQUATION_SOLVER");
    if(env) {
      nthread_v=atoi(env);
      if (nthread_v <= nthread) {nthread=nthread_v;}
    }
    if (nthread < 1) {nthread=1;}
    if(*iexpl<=1){
      sprintf(envMKL_cp,"MKL_NUM_THREADS=%" ITGFORMAT "",nthread);
    }
    putenv(envMKL_cp);
    mthread_mkl_cp=nthread;
  }
    
  if(*iexpl<=1){
    printf(" number of threads =% d\n\n",mthread_mkl_cp);
  }

  for(i=0;i<64;i++){ptcp[i]=0;}

  if(*symmetryflag==0){

    /* symmetric matrix; the subdiagonal entries are stored
       column by column in au, the diagonal entries in ad;
       pardiso needs the entries row per row */      

    mtype=-2;
      
    ndim=*neq+*nzs;
      
    NNEW(pointerscp,ITG,*neq+1);
    NNEW(icolpardisocp,ITG,ndim);
    NNEW(aupardisocp,double,ndim);
      
    k=ndim;
    l=*nzs;
      
    if(*sigma==0.){
      pointerscp[*neq]=ndim+1;
      for(i=*neq-1;i>=0;--i){
	for(j=0;j<icol[i];++j){
	  icolpardisocp[--k]=irow[--l];
	  aupardisocp[k]=au[l];
	}
	pointerscp[i]=k--;
	icolpardisocp[k]=i+1;
	aupardisocp[k]=ad[i];
      }
    }
    else{
      pointerscp[*neq]=ndim+1;
      for(i=*neq-1;i>=0;--i){
	for(j=0;j<icol[i];++j){
	  icolpardisocp[--k]=irow[--l];
	  aupardisocp[k]=au[l]-*sigma*aub[l];
	}
	pointerscp[i]=k--;
	icolpardisocp[k]=i+1;
	aupardisocp[k]=ad[i]-*sigma*adb[i];
      }
    }
  }else{

    if(*inputformat==3){

      /* off-diagonal terms  are stored column per
	 column from top to bottom in au;
	 diagonal terms are stored in ad  */

      /* structurally and numerically asymmetric */
	
      mtype=11;
	
      ndim=*neq+*nzs;
      NNEW(pointerscp,ITG,*neq+1);
      NNEW(irowpardisocp,ITG,ndim);	  
      NNEW(icolpardisocp,ITG,ndim);
      NNEW(aupardisocp,double,ndim);
	  
      k=0;
      k2=0;
      for(i=0;i<*neq;i++){
	for(j=0;j<icol[i];j++){
	  if(au[k]>1.e-12||au[k]<-1.e-12){
	    icolpardisocp[k2]=i+1;
	    irowpardisocp[k2]=irow[k];
	    aupardisocp[k2]=au[k];
	    k2++;		  
	  }
	  k++;	      
	}	  
      }  
      /* diagonal terms */  
      for(i=0;i<*neq;i++){
	icolpardisocp[k2]=i+1;
	irowpardisocp[k2]=i+1;
	aupardisocp[k2]=ad[i];
	k2++;	  
      }
      ndim=k2;
	  
      /* pardiso needs the entries row per row; so sorting is
	 needed */ 
	  
      kflag=2;
      FORTRAN(isortiid,(irowpardisocp,icolpardisocp,aupardisocp,
			&ndim,&kflag));
	  
      /* sorting each row */
	  
      k=0;
      pointerscp[0]=1;
      for(i=0;i<*neq;i++){
	j=i+1;
	kstart=k;
	do{
	  if(irowpardisocp[k]!=j ){
	    n=k-kstart;		  
	    FORTRAN(isortid,(&icolpardisocp[kstart],&aupardisocp[kstart],
			     &n,&kflag));
	    pointerscp[i+1]=k+1;
	    break;  
	  }else{
	    if(k+1==ndim){
	      n=k-kstart+1;	  
	      FORTRAN(isortid,(&icolpardisocp[kstart],
			       &aupardisocp[kstart],&n,&kflag));
	      break;	       
	    }else{
	      k++;	       
	    }  
	  }
	}while(1);
      }
      pointerscp[*neq]=ndim+1;
      SFREE(irowpardisocp);

    }else if(*inputformat==1){
	  
      /* lower triangular matrix is stored column by column in
	 au, followed by the upper triangular matrix row by row;
	 the diagonal terms are stored in ad */

      /* structurally symmetric, numerically asymmetric */
	
      mtype=1;
	
      /* reordering lower triangular matrix */

      ndim=*nzs;
      NNEW(pointerscp,ITG,*neq+1);
      NNEW(irowpardisocp,ITG,ndim);
      NNEW(icolpardisocp,ITG,ndim);
      NNEW(aupardisocp,double,ndim);
	  
      k=0;
      for(i=0;i<*neq;i++){
	for(j=0;j<icol[i];j++){
	  icolpardisocp[k]=i+1;
	  irowpardisocp[k]=irow[k];
	  aupardisocp[k]=au[k];
	  k++;
	}
      }
	  
      /* pardiso needs the entries row per row; so sorting is
	 needed */
	  
      kflag=2;
      FORTRAN(isortiid,(irowpardisocp,icolpardisocp,aupardisocp,
			&ndim,&kflag));
	  
      /* sorting each row */
	  
      k=0;
      pointerscp[0]=1;
      if(ndim>0){
	for(i=0;i<*neq;i++){
	  j=i+1;
	  kstart=k;
	  do{

	    /* end of row reached */

	    if(irowpardisocp[k]!=j){
	      n=k-kstart;
	      FORTRAN(isortid,(&icolpardisocp[kstart],&aupardisocp[kstart],
			       &n,&kflag));
	      pointerscp[i+1]=k+1;
	      break;
	    }else{

	      /* end of last row reached */

	      if(k+1==ndim){
		n=k-kstart+1;
		FORTRAN(isortid,(&icolpardisocp[kstart],&aupardisocp[kstart],
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
      pointerscp[*neq]=ndim+1;
      SFREE(irowpardisocp);

      /* composing the matrix: lower triangle + diagonal + upper triangle */

      ndim=*neq+2**nzs;
      RENEW(icolpardisocp,ITG,ndim);
      RENEW(aupardisocp,double,ndim);
      k=ndim;
      for(i=*neq-1;i>=0;i--){
	l=k+1;
	for(j=jq[i+1]-1;j>=jq[i];j--){
	  icolpardisocp[--k]=irow[j-1];
	  aupardisocp[k]=au[j+*nzs3-1];
	}
	icolpardisocp[--k]=i+1;
	aupardisocp[k]=ad[i];
	for(j=pointerscp[i+1]-1;j>=pointerscp[i];j--){
	  icolpardisocp[--k]=icolpardisocp[j-1];
	  aupardisocp[k]=aupardisocp[j-1];
	}
	pointerscp[i+1]=l;
      }
      pointerscp[0]=1;
    }
  }

/* next 3 lines are for the simulateous use of PARDISO and PaStiX */

#ifdef COMPANY
  mkl_domain_set_num_threads(mthread_mkl_cp,MKL_DOMAIN_PARDISO);
#endif

  FORTRAN(pardiso,(ptcp,&maxfct,&mnum,&mtype,&phase,neq,aupardisocp,
		   pointerscp,icolpardisocp,perm,&nrhs,iparmcp,&msglvl,
                   b,x,&error));

  return;
}

void pardiso_solve_cp(double *b, ITG *neq,ITG *symmetryflag,ITG *inputformat,
		   ITG *nrhs){

  ITG maxfct=1,mnum=1,phase=33,*perm=NULL,mtype,
    msglvl=0,i,error=0;
  double *x=NULL;

  /* if(*symmetryflag==0){
    printf(" Solving the system of equations using the symmetric pardiso solver\n");
  }else{
    printf(" Solving the system of equations using the unsymmetric pardiso solver\n");
    }*/

  if(*symmetryflag==0){
    mtype=-2;
  }else{
    if(*inputformat==3){
      mtype=11;
    }else{
      mtype=1;
    }
  }
  iparmcp[1]=3;
  
  /* pardiso_factor_cp has been called befor, MKL_NUM_THREADS=mthread_mkl_cp is set*/

  //printf(" number of threads =% d\n\n",mthread_mkl_cp);

  NNEW(x,double,*nrhs**neq);

  FORTRAN(pardiso,(ptcp,&maxfct,&mnum,&mtype,&phase,neq,aupardisocp,
		   pointerscp,icolpardisocp,perm,nrhs,iparmcp,&msglvl,
                   b,x,&error));

  for(i=0;i<*nrhs**neq;i++){b[i]=x[i];}
  SFREE(x);

  return;
}

void pardiso_cleanup_cp(ITG *neq,ITG *symmetryflag,ITG *inputformat){

  ITG maxfct=1,mnum=1,phase=-1,*perm=NULL,nrhs=1,mtype,
    msglvl=0,error=0;
  double *b=NULL,*x=NULL;

  if(*symmetryflag==0){
    mtype=-2;
  }else{
    if(*inputformat==3){
      mtype=11;
    }else{
      mtype=1;
    }
  }

  FORTRAN(pardiso,(ptcp,&maxfct,&mnum,&mtype,&phase,neq,aupardisocp,
		   pointerscp,icolpardisocp,perm,&nrhs,iparmcp,&msglvl,
                   b,x,&error));

  SFREE(icolpardisocp);
  SFREE(aupardisocp);
  SFREE(pointerscp);

  return;
}

void pardiso_main_cp(double *ad, double *au, double *adb, double *aub, 
		  double *sigma,double *b, ITG *icol, ITG *irow, 
		  ITG *neq, ITG *nzs,ITG *symmetryflag,ITG *inputformat,
		  ITG *jq, ITG *nzs3,ITG *nrhs){

  ITG iexpl=0;
  
  if(*neq==0) return;

  pardiso_factor_cp(ad,au,adb,aub,sigma,icol,irow, 
		    neq,nzs,symmetryflag,inputformat,jq,nzs3,&iexpl);

  pardiso_solve_cp(b,neq,symmetryflag,inputformat,nrhs);

  pardiso_cleanup_cp(neq,symmetryflag,inputformat);

  return;
}

#endif

