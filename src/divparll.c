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

/*     Parallellization of the division of an array by another array,
       element per element	 */


#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

static ITG *neapar=NULL,*nebpar=NULL;

static double *var11=NULL,*var21=NULL;

void divparll(double *var1,double *var2,ITG *isize,ITG *num_cpus){

  
//    static ITG i,num_cpus,nepar,idelta,isum,*ipar=NULL;

    ITG i,idelta,isum;

    /* variables for multithreading procedure */

    ITG *ithread=NULL;

    pthread_t tid[*num_cpus];

    /* determining the element bounds in each thread */

    NNEW(neapar,ITG,*num_cpus);
    NNEW(nebpar,ITG,*num_cpus);

    /* dividing the element number range into num_cpus equal numbers of 
       active entries.  */

    idelta=(ITG)floor(*isize/(double)(*num_cpus));
    isum=0;
    for(i=0;i<*num_cpus;i++){
	neapar[i]=isum;
	if(i!=*num_cpus-1){
	    isum+=idelta;
	}else{
	    isum=*isize;
	}
	nebpar[i]=isum;
    }

    /* create threads and wait */
    var11=var1;var21=var2;
    NNEW(ithread,ITG,*num_cpus);

    for(i=0; i<*num_cpus; i++)  {
      ithread[i]=i;
      pthread_create(&tid[i], NULL, (void *)divparllmt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);

    SFREE(ithread);SFREE(neapar);SFREE(nebpar);

}

/* subroutine for multithreading of divparllmt */

void *divparllmt(ITG *i){

    ITG nea,neb,j;

    nea=neapar[*i];
    neb=nebpar[*i];

    for(j=nea;j<neb;j++){
    	var11[j]=var11[j]/var21[j];
    }

    return NULL;
}
