/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2021 Guido Dhondt                     */

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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CalculiX.h"

static ITG *nk1,*nactdof1,num_cpus,mt1;

static double *vold1,*volddof1;

void resforccont(double *vold,ITG *nk,ITG *mi,double *aubi,ITG *irowbi,
		 ITG *jqbi,ITG *neqtot,ITG *ktot,double *fext,double *gapdisp,
		 double *auib,ITG *irowib,ITG *jqib,ITG *nactdof,
		 double *volddof,ITG *neq,double *qi_kbi){

  /* calculate the residual contact force fext-Kbi*Ui */

  ITG sys_cpus,*ithread=NULL,i,itranspose,mt=mi[1]+1;
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

  // next line is to be inserted in a similar way for all other parallel parts

  if(*nk<num_cpus) num_cpus=*nk;

  pthread_t tid[num_cpus];

  nk1=nk;nactdof1=nactdof;volddof1=volddof;vold1=vold;mt1=mt;

  /* create threads and wait */
	
  NNEW(ithread,ITG,num_cpus);
  for(i=0; i<num_cpus; i++)  {
    ithread[i]=i;
    pthread_create(&tid[i], NULL, (void *)resforccontmt, (void *)&ithread[i]);
  }
  for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);

  SFREE(ithread);

  /* We compute g as volddof=(Kbi*volddof)+(Kib*volddof) in qi_kbi
     to account for the missing terms due to the low triangle structure
     of the matrices */

  /* calculate Kbi*volddof */

  itranspose=0;
  mulmatvec_asymmain(aubi,jqbi,irowbi,&neq[0],volddof,qi_kbi,
			  &itranspose,neqtot);

  /* calculate Kib^T*volddof and add to g.
     transposed multiplication */

  itranspose=1;
  mulmatvec_asymmain(auib,jqib,irowib,neqtot,volddof,qi_kbi,
			  &itranspose,neqtot);

  /* add external force of BOUNDARY DOF */

  for(i=0;i<*neqtot;i++){
    gapdisp[i]=fext[ktot[i]-1]-qi_kbi[i];
  }
  
  return;
}

/* collecting the results */

void *resforccontmt(ITG *i){

  /* vold contains the displacements in the order of the nodes;
     volddof contains the displacements in the order of the
     degrees of freedom */

  ITG idelta,nka,nkb,j,k,index;
  
  idelta=(ITG)ceil(*nk1/(double)num_cpus);
  nka=*i*idelta;
  nkb=(*i+1)*idelta;
  if(nkb>*nk1) nkb=*nk1;

  for(j=nka;j<nkb;j++){
    for(k=1;k<4;k++){
      index=j*mt1+k;
      if(nactdof1[index]>0){
	volddof1[nactdof1[index]-1]=vold1[index];
      }
    }
  }
  
  return NULL;
}
