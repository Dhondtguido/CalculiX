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

/*     A parallel part in prediction.c	 */


#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

static ITG *nkapar=NULL,*nkbpar=NULL,*mt1,*nactdof1;

static double *dtime1,*veold1,*scal11,*accold1,*uam1,*v1,*vold1,*scal21;

void preparll(ITG *mt,double *dtime,double *veold,double *scal1,
		   double *accold,double *uam,ITG *nactdof,double *v,
		   double *vold,double *scal2,ITG *nk,ITG *num_cpus){

    ITG i,idelta,isum;

    /* variables for multithreading procedure */

    ITG *ithread=NULL;

    pthread_t tid[*num_cpus];

    /* determining the node bounds in each thread */

    NNEW(nkapar,ITG,*num_cpus);
    NNEW(nkbpar,ITG,*num_cpus);
    NNEW(uam1,double,*num_cpus);

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
    
    mt1=mt;nactdof1=nactdof;dtime1=dtime;veold1=veold;scal11=scal1;
    accold1=accold;v1=v;vold1=vold;scal21=scal2;
    
    NNEW(ithread,ITG,*num_cpus);

    for(i=0; i<*num_cpus; i++)  {
      ithread[i]=i;
      pthread_create(&tid[i], NULL, (void *)preparllmt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);

    for(i=0;i<*num_cpus;i++){
	if(uam1[i]>uam[0]){uam[0]=uam1[i];}
    }

    SFREE(ithread);SFREE(nkapar);SFREE(nkbpar);SFREE(uam1);

}

/* subroutine for multithreading of preparll */

void *preparllmt(ITG *i){

    ITG nka,nkb,k,j;

    double dextrapol;

    nka=nkapar[*i];
    nkb=nkbpar[*i];

    for(k=nka;k<nkb;k++){
	for(j=1;j<*mt1;j++){
//	for(j=0;j<*mt1;j++){
	    dextrapol=*dtime1*veold1[*mt1*k+j]+(*scal11)*accold1[*mt1*k+j];
	    if((fabs(dextrapol)>uam1[*i])&&(nactdof1[*mt1*k+j]>0)) {uam1[*i]=fabs(dextrapol);}
	    v1[*mt1*k+j]=vold1[*mt1*k+j]+dextrapol;
	    veold1[*mt1*k+j]=veold1[*mt1*k+j]+(*scal21)*accold1[*mt1*k+j];
	    accold1[*mt1*k+j]=0.;
	}
    }

    return NULL;
}
