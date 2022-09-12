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

static ITG *neq1,*num_cpus1;

static double *sol1,*b1,*adl1;

void solveeq(double *adb,double *aub,double *adl,double *b,double *sol,
	     double *aux,ITG *irow,ITG *jq,ITG *neq,ITG *maxit,
	     ITG *num_cpus){

  /* 
     solving a system of equations by iteratively solving the
     lumped version
     The diagonal terms f the original system are stored in adb,
     the off-diagonal terms in aub
     Ref: The Finite Element Method for Fluid Dynamics,
          O.C. Zienkiewicz, R.L. Taylor & P. Nithiarasu
          6th edition (2006) ISBN 0 7506 6322 7
          p. 61
  */
  
  ITG i,k,*ithread=NULL;
  
  pthread_t tid[*num_cpus];

  /* first iteration */

  neq1=neq;sol1=sol;b1=b;adl1=adl;num_cpus1=num_cpus;
 
  NNEW(ithread,ITG,*num_cpus);
  for(i=0; i<*num_cpus; i++)  {
    ithread[i]=i;
    pthread_create(&tid[i], NULL, (void *)solveeqparmt, (void *)&ithread[i]);
  }
  for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
      
  SFREE(ithread);
  
  /*for(i=0;i<*neq;i++){
    sol[i]=b[i]*adl[i];
    }*/
  
  if(*maxit==1) return;

  /* iterating maxit times */
  
  for(k=1;k<*maxit;k++){

    /* multiplying the difference of the original matrix
       with the lumped matrix with the actual solution */

    FORTRAN(op,(neq,sol,aux,adb,aub,jq,irow));

    for(i=0;i<*neq;i++){
      sol[i]=(b[i]-aux[i])*adl[i];
    }
    
  }
  
  return;
}

/* subroutine for collecting results */

void *solveeqparmt(ITG *i){

  ITG neqa,neqb,neqdelta,j;
    
  neqdelta=(ITG)ceil(*neq1/(double)(*num_cpus1));
  neqa=*i*neqdelta;
  neqb=neqa+neqdelta;
  if(neqb>*neq1) neqb=*neq1;
  
  for(j=neqa;j<neqb;j++){
    sol1[j]=b1[j]*adl1[j];
  }

  return NULL;
}
