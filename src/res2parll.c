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

/*     Parallellization of subroutine of calcresidual.c          	 */


#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

static ITG *nqapar=NULL,*nqbpar=NULL;

static double *b1,*scal11,*fext1,*f1,*alpha1,*fextini1,*fini1,
    *adb1,*aux21;

void res2parll(double *b,double *scal1,double *fext,double *f,
		    double *alpha,double *fextini,double *fini,
		    double *adb,double *aux2,ITG *neq0,ITG *num_cpus){

    ITG i,idelta,isum;

    /* variables for multithreading procedure */

    ITG *ithread=NULL;

    pthread_t tid[*num_cpus];

    /* determining the equation bounds in each thread */

    NNEW(nqapar,ITG,*num_cpus);
    NNEW(nqbpar,ITG,*num_cpus);

    /* dividing the equation number range into num_cpus equal numbers of 
       active entries.  */

    idelta=(ITG)floor(*neq0/(double)(*num_cpus));
    isum=0;
    for(i=0;i<*num_cpus;i++){
	nqapar[i]=isum;
	if(i!=*num_cpus-1){
	    isum+=idelta;
	}else{
	    isum=*neq0;
	}
	nqbpar[i]=isum;
    }

    /* create threads and wait */
    
    b1=b;scal11=scal1;fext1=fext;f1=f;alpha1=alpha;fextini1=fextini;
    fini1=fini;adb1=adb;aux21=aux2;
    
    NNEW(ithread,ITG,*num_cpus);

    for(i=0; i<*num_cpus; i++)  {
      ithread[i]=i;
      pthread_create(&tid[i], NULL, (void *)res2parllmt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);

    SFREE(ithread);SFREE(nqapar);SFREE(nqbpar);

}

/* subroutine for multithreading of res2parll */

void *res2parllmt(ITG *i){

    ITG nqa,nqb,k;

    nqa=nqapar[*i];
    nqb=nqbpar[*i];
    
    for(k=nqa;k<nqb;++k){
	b1[k]=*scal11*(fext1[k]-f1[k])-alpha1[0]*(fextini1[k]-fini1[k])
	    -adb1[k]*aux21[k];
    }

    return NULL;
}
