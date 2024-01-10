/*     CalculiX - A 3-dimensional finite element program                 */
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

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

static ITG *jq1,*irow1,num_cpus,*n1,*ncol1,*itranspose1;

static double *x1,*yy=NULL,*au1,*yy1;

void mulmatvec_asymmain(double *au,ITG *jq,ITG *irow,ITG *ncol,double *x,
			double *y,ITG *itranspose,ITG *n){

  /* calculates y=y+au*x with
     ncol: number of columns in au and rows in x
     n: number of rows in y and au
     jq(i): location of first element in column j in au
     au(k): element in matrix au
     irow(k): corresponding row 

     if itranspose=1, y=y+au^T*x is calculated and now we have:
     ncol: number of columns in au and rows in y
     n: number of rows in au and x */
  
  ITG sys_cpus,*ithread=NULL,i;
  char *env,*envloc,*envsys;

  num_cpus = 0;
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

  envloc = getenv("CCX_NPROC_RESULTS");
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

  // next line is to be inserted in a similar way for all other paralell parts

  if(*ncol<num_cpus) num_cpus=*ncol;

  pthread_t tid[num_cpus];
  
  /* calculating the product of the matrix given by ad,au with vector x */

  if(*itranspose==0){
  
    /* allocating memory for the solution vector yy */

    NNEW(yy,double,num_cpus**n);

    au1=au;jq1=jq;irow1=irow;ncol1=ncol;x1=x;itranspose1=itranspose;
    n1=n;

    /* create threads and wait */
	
    NNEW(ithread,ITG,num_cpus);
    for(i=0; i<num_cpus; i++)  {
      ithread[i]=i;
      pthread_create(&tid[i], NULL, (void *)mulmatvec_asymmt1, (void *)&ithread[i]);
    }
    for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);

    SFREE(ithread);

    yy1=y;n1=n;
  
    /* collecting results: create threads and wait */
	
    NNEW(ithread,ITG,num_cpus);
    for(i=0; i<num_cpus; i++)  {
      ithread[i]=i;
      pthread_create(&tid[i], NULL, (void *)mulmatvec_asymct1, (void *)&ithread[i]);
    }
    for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);

    SFREE(ithread);SFREE(yy);
    
  }else{

    au1=au;jq1=jq;irow1=irow;ncol1=ncol;x1=x;yy1=y;itranspose1=itranspose;
    n1=n;

    /* create threads and wait */
	
    NNEW(ithread,ITG,num_cpus);
    for(i=0; i<num_cpus; i++)  {
      ithread[i]=i;
      pthread_create(&tid[i], NULL, (void *)mulmatvec_asymmt2, (void *)&ithread[i]);
    }
    for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);

    SFREE(ithread);

  }
	
}

/* subroutine for multithreading of mulmatvec_asym (itranspose=0) */

void *mulmatvec_asymmt1(ITG *i){

  ITG indexf,idelta,ncola,ncolb;

  indexf=*i**n1;

  idelta=(ITG)ceil(*ncol1/(double)num_cpus);
  ncola=*i*idelta+1;
  ncolb=(*i+1)*idelta;
  if(ncolb>*ncol1) ncolb=*ncol1;

  FORTRAN(mulmatvec_asym,(au1,jq1,irow1,x1,&yy[indexf],itranspose1,
			  &ncola,&ncolb));

  return NULL;
}

/* collecting the results (only needed for itranspose=0) */

void *mulmatvec_asymct1(ITG *i){

  ITG idelta,na,nb,j,k,index;

  idelta=(ITG)ceil(*n1/(double)num_cpus);
  na=*i*idelta;
  nb=(*i+1)*idelta;
  if(nb>*n1) nb=*n1;

  for(k=0;k<num_cpus;k++){
    index=k**n1;
    for(j=na;j<nb;j++){
      yy1[j]+=yy[j+index];
    }
  }

  return NULL;
}

/* subroutine for multithreading of mulmatvec_asym (itranspose=1) */

void *mulmatvec_asymmt2(ITG *i){

  ITG idelta,ncola,ncolb;

  idelta=(ITG)ceil(*ncol1/(double)num_cpus);
  ncola=*i*idelta+1;
  ncolb=(*i+1)*idelta;
  if(ncolb>*ncol1) ncolb=*ncol1;

  FORTRAN(mulmatvec_asym,(au1,jq1,irow1,x1,yy1,itranspose1,
			  &ncola,&ncolb));

  return NULL;
}

