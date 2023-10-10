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

static ITG *jq1,*irow1,num_cpus,*n1;

static double *x1,*yy1=NULL,*ad1,*au1;

void opmain(ITG *n,double *x,double *y,double *ad,double*au,ITG *jq,ITG *irow){

  ITG sys_cpus,*ithread=NULL,i,j,isum,idelta;
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

  if(*n<num_cpus) num_cpus=*n;

  pthread_t tid[num_cpus];

  /* allocating memory for the solution vector y */

  NNEW(yy1,double,num_cpus**n);

  x1=x;ad1=ad;au1=au;jq1=jq;irow1=irow;n1=n;
  
  /* calculating the product of the matrix given by ad,au with vector x */

    printf(" Using up to %" ITGFORMAT " cpu(s) for the matrix times vector multiplication.\n\n", num_cpus);

  /* create threads and wait */
  /* distributed force: body foce*/
	
  NNEW(ithread,ITG,num_cpus);
  for(i=0; i<num_cpus; i++)  {
    ithread[i]=i;
    pthread_create(&tid[i], NULL, (void *)opmt, (void *)&ithread[i]);
  }
  for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);

  for(i=0;i<*n;i++){
    y[i]=yy1[i];
  }
  for(i=0;i<*n;i++){
    for(j=1;j<num_cpus;j++){
      y[i]+=yy1[i+j**n];
    }
  }

  SFREE(ithread);SFREE(yy1);
	
}

/* subroutine for multithreading of op */

void *opmt(ITG *i){

  ITG indexf,idelta,na,nb;

  indexf=*i**n1;

  idelta=(ITG)ceil(*n1/(double)num_cpus);
  na=*i*idelta+1;
  nb=(*i+1)*idelta;
  if(nb>*n1) nb=*n1;

  FORTRAN(op,(x1,&yy1[indexf],ad1,au1,jq1,irow1,&na,&nb));

  return NULL;
}

