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

static char *lakon1;

static ITG *kon1,*ipkon1,*ne1,*mi1,num_cpus,*nkapar=NULL,*nkepar=NULL;

static double *co1,*qfx1,*h01;

void biosav(ITG *ipkon,ITG *kon,char *lakon,ITG *ne,double *co,
                double *qfx,double *h0,ITG *mi,ITG *inomat,ITG *nk){

    ITG i,j,*ithread=NULL,nkphi,idelta,isum;

    /* calculates the magnetic intensity due to currents in the phi-
       domain of an electromagnetic calculation */
      
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

    envloc = getenv("CCX_NPROC_BIOTSAVART");
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
	    num_cpus=atoi(env);
	if (num_cpus<1) {
	    num_cpus=1;
	}else if(num_cpus>sys_cpus){
	    num_cpus=sys_cpus;
	}
    }
    
    /* determining the nodal bounds in each thread */

    NNEW(nkapar,ITG,num_cpus);
    NNEW(nkepar,ITG,num_cpus);

    /* n1 is the number of nodes in the phi(magnetostatic)-domain in
       an electromagnetic calculation */

    nkphi=0;
    for(i=0;i<*nk;i++){
	if(inomat[i]==1) nkphi++;
    }
    if(nkphi<num_cpus) num_cpus=nkphi;

    idelta=nkphi/num_cpus;
    
    /* dividing the range from 1 to the number of phi-nodes */

    isum=0;
    for(i=0;i<num_cpus;i++){
	nkapar[i]=isum;
	if(i!=num_cpus-1){
	    isum+=idelta;
	}else{
	    isum=nkphi;
	}
	nkepar[i]=isum-1;
    }
    
    /* translating the bounds of the ranges to real node numbers */

//    i=0;
    i=-1;
    j=0;
    nkphi=-1;

    do{
	if(j==num_cpus) break;
	do{
	    if(nkapar[j]==nkphi){
		nkapar[j]=i;
		break;
	    }else{
		do{
		    i++;
		    if(inomat[i]==1){
			nkphi++;
			break;
		    }
		}while(1);
	    }
	}while(1);

	do{
	    if(nkepar[j]==nkphi){
		nkepar[j]=i;
		j++;
		break;
	    }else{
		do{
		    i++;
		    if(inomat[i]==1){
			nkphi++;
			break;
		    }
		}while(1);
	    }
	}while(1);
    }while(1);

    ipkon1=ipkon;kon1=kon;lakon1=lakon;ne1=ne;co1=co;qfx1=qfx;
    h01=h0;mi1=mi;
    
    printf(" Using up to %" ITGFORMAT " cpu(s) for the Biot-Savart calculation.\n\n", num_cpus);
    
    /* create threads and wait */
    
    pthread_t tid[num_cpus];
    
    NNEW(ithread,ITG,num_cpus);
    for(i=0;i<num_cpus;i++){
	ithread[i]=i;
	pthread_create(&tid[i],NULL,(void *)biotsavartmt,(void *)&ithread[i]);
    }
    for(i=0;i<num_cpus;i++)pthread_join(tid[i], NULL);
    
    SFREE(ithread);SFREE(nkapar);SFREE(nkepar);
    
    return;
    
}

/* subroutine for multithreading of biotsavart */

void *biotsavartmt(ITG *i){

    ITG nka,nkb;

    nka=nkapar[*i]+1;
    nkb=nkepar[*i]+1;

    FORTRAN(biotsavart,(ipkon1,kon1,lakon1,ne1,co1,qfx1,h01,mi1,&nka,
			&nkb));

    return NULL;
}

