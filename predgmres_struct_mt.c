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
#include <string.h>
#include "CalculiX.h"

void predgmres_struct_mt(double *ad, double **aup, double *adb, double *aub, 
			 double *sigma,double *b,ITG *icol, ITG *irow, 
			 ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
			 ITG *jq, ITG *nzs3,ITG *nrhs){

  ITG i,j,k,l,*irowpardiso=NULL,kflag,kstart,n,
      k2,ndim,num_cpus,isym,itol,itmax,iter,ierr,iunit,lrgw,
      *igwk=NULL,ligw,*iwork=NULL,*icolpardiso=NULL,*ipointers=NULL,
      *nestart=NULL,*ja=NULL,nz_num;

  double *sb=NULL,*sx=NULL,*rgwk=NULL,tol,err,*rwork=NULL,*x=NULL,
    *aupardiso=NULL,*au=NULL,dgmrestol=1.e-12;

  au=*aup;

  /* variables for multithreading procedure */

  ITG sys_cpus;
  char *env,*envloc,*envsys;
      
  num_cpus=0;
  sys_cpus=0;
  
  /* explicit user declaration prevails */
  
  envsys=getenv("NUMBER_OF_CPUS");
  if(envsys){
      sys_cpus=atoi(envsys);
      if(sys_cpus<0) sys_cpus=0;
  }
  
  /* automatic detection of available number of processors */
  
  if(sys_cpus==0){
      sys_cpus = getSystemCPUs();
      if(sys_cpus<1) sys_cpus=1;
  }
  
  /* local declaration prevails, if strictly positive */
  
  envloc = getenv("CCX_NPROC_EQUATION_SOLVER");
  if(envloc){
      num_cpus=atoi(envloc);
      if(num_cpus<0){
	  num_cpus=0;
      }else if(num_cpus>sys_cpus){
	  num_cpus=sys_cpus;
      }
  }
  
  /* else global declaration, if any, applies */
  
  env = getenv("OMP_NUM_THREADS");
  if(num_cpus==0){
      if (env)
	  num_cpus = atoi(env);
      if (num_cpus < 1) {
	  num_cpus=1;
      }else if(num_cpus>sys_cpus){
	  num_cpus=sys_cpus;
      }
  }
  
  if(*neq<num_cpus) num_cpus=*neq;

  printf(" Factoring and solving the system of equations using the gmres solver\n");
  printf(" using up to %" ITGFORMAT " cpu(s).\n\n", num_cpus);

  /* rearranging the matrix format in a row-by-row format */

  if(*symmetryflag==0){

      /* symmetric matrix; the subdiagonal entries are stored
         column by column in au, the diagonal entries in ad;
         pardiso needs the entries row per row */ 

      /* adding the symmetric contributions */

      RENEW(au,double,2**nzs);
      memcpy(&au[*nzs],&au[0],sizeof(double)**nzs);
      
      /* that way this case is reduced to an asymmetric matrix
         with inputformat=1 (cf. further down; code was copied) */

      /* reordering lower triangular matrix */
      
      ndim=*nzs;
      NNEW(ipointers,ITG,*neq+1);
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
      ipointers[0]=1;
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
		      ipointers[i+1]=k+1;
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
      ipointers[*neq]=ndim+1;
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
	      aupardiso[k]=au[j+*nzs3-1];
	  }
	  icolpardiso[--k]=i+1;
	  aupardiso[k]=ad[i];
	  for(j=ipointers[i+1]-1;j>=ipointers[i];j--){
	      icolpardiso[--k]=icolpardiso[j-1];
	      aupardiso[k]=aupardiso[j-1];
	  }
	  ipointers[i+1]=l;
      }
      ipointers[0]=1;
      
      /* end copying */
      
      RENEW(au,double,*nzs);

  }else{

      if(*inputformat==3){

          /* off-diagonal terms  are stored column per
             column from top to bottom in au;
             diagonal terms are stored in ad  */

	  ndim=*neq+*nzs;
	  NNEW(ipointers,ITG,*neq+1);
	  NNEW(irowpardiso,ITG,ndim);	  
	  NNEW(icolpardiso,ITG,ndim);
	  NNEW(aupardiso,double,ndim);
	  
	  k=0;
	  k2=0;
	  for(i=0;i<*neq;i++){
	      for(j=0;j<icol[i];j++){
		  if(au[k]>1.e-12||au[k]<-1.e-12){
		      icolpardiso[k2]=i+1;
		      irowpardiso[k2]=irow[k];
		      aupardiso[k2]=au[k];
		      k2++;		  
		  }
		  k++;	      
	      }	  
	  }  
	  /* diagonal terms */  
	  for(i=0;i<*neq;i++){
	      icolpardiso[k2]=i+1;
	      irowpardiso[k2]=i+1;
	      aupardiso[k2]=ad[i];
	      k2++;	  
	  }
	  ndim=k2;
	  
	  /* pardiso needs the entries row per row; so sorting is
	     needed */ 
	  
	  kflag=2;
	  FORTRAN(isortiid,(irowpardiso,icolpardiso,aupardiso,
			    &ndim,&kflag));
	  
	  /* sorting each row */
	  
	  k=0;
	  ipointers[0]=1;
	  for(i=0;i<*neq;i++){
	      j=i+1;
	      kstart=k;
	      do{
		  if(irowpardiso[k]!=j ){
		      n=k-kstart;		  
		      FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
				       &n,&kflag));
		      ipointers[i+1]=k+1;
		      break;  
		  }else{
		      if(k+1==ndim){
			  n=k-kstart+1;	  
			  FORTRAN(isortid,(&icolpardiso[kstart],
                                  &aupardiso[kstart],&n,&kflag));
			  break;	       
		      }else{
			  k++;	       
		      }  
		  }
	      }while(1);
	  }
	  ipointers[*neq]=ndim+1;
	  SFREE(irowpardiso);

      }else if(*inputformat==1){
	  
          /* lower triangular matrix is stored column by column in
             au, followed by the upper triangular matrix row by row;
             the diagonal terms are stored in ad */

          /* reordering lower triangular matrix */

	  ndim=*nzs;
	  NNEW(ipointers,ITG,*neq+1);
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
	  ipointers[0]=1;
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
			  ipointers[i+1]=k+1;
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
	  ipointers[*neq]=ndim+1;
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
		  aupardiso[k]=au[j+*nzs3-1];
	      }
	      icolpardiso[--k]=i+1;
	      aupardiso[k]=ad[i];
	      for(j=ipointers[i+1]-1;j>=ipointers[i];j--){
		  icolpardiso[--k]=icolpardiso[j-1];
		  aupardiso[k]=aupardiso[j-1];
	      }
	      ipointers[i+1]=l;
	  }
	  ipointers[0]=1;
      }
  }

  NNEW(rwork,double,*neq);
  for(i=0;i<*neq;i++){rwork[i]=1./ad[i];}

  NNEW(nestart,ITG,num_cpus+1);
  NNEW(ja,ITG,*neq+1);

  /* splitting the matrix in disjunct blocks;
     entries outside these blocks are moved to the rhs */

  FORTRAN(createblock_struct,(neq,ipointers,icolpardiso,aupardiso,nestart,
			      &num_cpus,ja,&nz_num));

  RENEW(icolpardiso,ITG,nz_num);
  RENEW(aupardiso,double,nz_num);
  SFREE(ipointers);

  /* the start solution is zero since each iteration involves the calculation
     of an increment to the global solution */

  NNEW(x,double,*neq);

  dgmresmain(neq,b,x,&nz_num,icolpardiso,ja,aupardiso,&isym,&itol,&tol,&itmax,
	     &iter,&err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,&ligw,rwork,
	     iwork,nestart,&num_cpus,&dgmrestol);

  if(ierr>0){
      printf(" *WARNING in predgmres_struct_mt: error message from dgmresmain =%d\n\n",ierr);
  }

  memcpy(&b[0],&x[0],sizeof(double)**neq);
  SFREE(x);SFREE(nestart);SFREE(icolpardiso);SFREE(aupardiso);SFREE(ja);
  SFREE(rwork);

  *aup=au;

  return;
}
