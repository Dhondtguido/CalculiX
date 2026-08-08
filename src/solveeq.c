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
#include "CalculiX.h"

void solveeq(double *adb,double *aub,double *adl,double *b,double *sol,
	     double *aux,ITG *irow,ITG *jq,ITG *neq,ITG *maxit,
	     ITG *num_cpus){

  /* 
     solving a system of equations by iteratively solving the
     lumped version
     The diagonal terms of the original system are stored in adb,
     the off-diagonal terms in aub
     Ref: The Finite Element Method for Fluid Dynamics,
          O.C. Zienkiewicz, R.L. Taylor & P. Nithiarasu
          6th edition (2006) ISBN 0 7506 6322 7
          p. 61
  */
  
  ITG i,k;
  
  /* first iteration */

  #pragma omp parallel for num_threads(*num_cpus)
  for(i=0;i<*neq;i++){
    sol[i]=b[i]*adl[i];
  }
  
  if(*maxit==1) return;

  /* iterating maxit times */
  
  for(k=1;k<*maxit;k++){

    /* multiplying the difference of the original matrix
       with the lumped matrix with the actual solution */

    FORTRAN(spmv,(neq,sol,aux,adb,aub,jq,irow,num_cpus));

    #pragma omp parallel for num_threads(*num_cpus)
    for(i=0;i<*neq;i++){
      sol[i]=(b[i]-aux[i])*adl[i];
    }
    
  }
  
  return;
}
