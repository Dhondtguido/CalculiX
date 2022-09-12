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

static ITG *nkapar=NULL,*nkbpar=NULL,*mt1,*nactdof1;

static double *aux21,*vold1,*vini1,*dtime1,*accold1;

void res1parll(ITG *mt,ITG *nactdof,double *aux2,double *vold,
		    double *vini,double *dtime,double *accold,
		    ITG *nk,ITG *num_cpus){

    ITG i,idelta,isum;

    /* variables for multithreading procedure */

    ITG *ithread=NULL;

    pthread_t tid[*num_cpus];

    /* determining the node bounds in each thread */

    NNEW(nkapar,ITG,*num_cpus);
    NNEW(nkbpar,ITG,*num_cpus);

    /* dividing the node number range into num_cpus equal numbers of 
       active entries.  */

    idelta=(ITG)floor(*nk/(double)(*num_cpus));
    isum=0;
    for(i=0;i<*num_cpus;i++){
	nkapar[i]=isum;
	if(i!=*num_cpus-1){
	    isum+=idelta;
	}else{
	    isum=*nk;
	}
	nkbpar[i]=isum;
    }

    /* create threads and wait */
    
    mt1=mt;nactdof1=nactdof;aux21=aux2;vold1=vold;vini1=vini;dtime1=dtime;
    accold1=accold;
    
    NNEW(ithread,ITG,*num_cpus);

    for(i=0; i<*num_cpus; i++)  {
      ithread[i]=i;
      pthread_create(&tid[i], NULL, (void *)res1parllmt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);

    SFREE(ithread);SFREE(nkapar);SFREE(nkbpar);

}

/* subroutine for multithreading of res1parll */

void *res1parllmt(ITG *i){

    ITG nka,nkb,k,j;

    nka=nkapar[*i];
    nkb=nkbpar[*i];
    
    for(k=nka;k<nkb;++k){
	if(nactdof1[*mt1*k]>0){
	    aux21[nactdof1[*mt1*k]-1]=(vold1[*mt1*k]-vini1[*mt1*k])/(*dtime1);}
	for(j=1;j<*mt1;++j){
	    if(nactdof1[*mt1*k+j]>0){aux21[nactdof1[*mt1*k+j]-1]=accold1[*mt1*k+j];}
	}
    }

    return NULL;
}
