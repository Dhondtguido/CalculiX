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

static ITG *nstate1_,*mi1,*islavsurf1,*ne01,*islavsurfold1,numstart,numfaces,
  num_cpus;

static double *xstate1,*xstateini1,*pslavsurf1,*pslavsurfold1;

void interpolatestatemain(ITG *ne,ITG *ipkon,ITG *kon,char *lakon,ITG *ne0,
			  ITG *mi,double *xstate,double *pslavsurf,
			  ITG *nstate_,double *xstateini,ITG *islavsurf,
			  ITG *islavsurfold,double *pslavsurfold,char *tieset,
			  ITG *ntie,ITG *itiefac){

  ITG sys_cpus,*ithread=NULL,i,j,num_cpus_ref;
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

  envloc = getenv("CCX_NPROC_INTERPOLSTATE");
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

  /* reference number of cpus; can be decrease in the following loops
     due to a too small amount of faces per contact tie */
  
  num_cpus_ref=num_cpus;

  /* loop over all ties */
  
  for(i=0;i<*ntie;i++){
    if(strcmp1(&tieset[i*243+80],"C")!=0) continue;

    /* parallellizing the loop over the contact slave faces for each tie
       - numstart is the location in islavsurf(old) befor the first
       face of the present contact tie
       - numfaces is the total number of contact faces in the present
       contact tie */
    
    numstart=itiefac[2*i]-1;
    numfaces=itiefac[2*i+1]-numstart;
  
    if(numfaces<num_cpus_ref){
      num_cpus=numfaces;
    }else{
      num_cpus=num_cpus_ref;
    }

    pthread_t tid[num_cpus];
	
    xstate1=xstate;xstateini1=xstateini;nstate1_=nstate_;mi1=mi;
    islavsurf1=islavsurf;pslavsurf1=pslavsurf;ne01=ne0;
    islavsurfold1=islavsurfold;pslavsurfold1=pslavsurfold;
 
    /* create threads and wait */
    /* distributed force: body foce*/
	
    NNEW(ithread,ITG,num_cpus);
    for(j=0; j<num_cpus; j++){
      ithread[j]=j;
      pthread_create(&tid[j], NULL, (void *)interpolatestatemainmt,
		     (void *)&ithread[j]);
    }
    for(j=0; j<num_cpus; j++)  pthread_join(tid[j], NULL);

    SFREE(ithread);
  }
	
}

/* subroutine for multithreading of rhs: distributed forces */

void *interpolatestatemainmt(ITG *i){

  ITG nfaa,nfab,nfadelta,kk,numpts;

  nfadelta=(ITG)floor(numfaces/(double)num_cpus);
  nfaa=numstart+(*i)*nfadelta+1;
  nfab=numstart+(*i+1)*nfadelta;
    
  if((*i==num_cpus-1)&&(nfab<numstart+numfaces))
    nfab=numstart+numfaces;

  for(kk=nfaa;kk<=nfab;kk++){
    numpts=islavsurfold1[2*kk+1]-islavsurfold1[2*kk-1];

    /* interpolating the state variables from the integration point 
       scheme from the last increment onto the integration point scheme
       in the present increment */

    if(numpts>2){
      FORTRAN(interpolateinface,(&kk,xstate1,xstateini1,&numpts,nstate1_,mi1,
				 islavsurf1,pslavsurf1,ne01,islavsurfold1,
				 pslavsurfold1));
    }
  }

  return NULL;
}

