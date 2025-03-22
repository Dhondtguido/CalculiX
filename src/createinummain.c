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
#include "mortar.h"

static char *cflag1;

static ITG *ipkon1,*inum1,*kon1,*nk1,*ne1,*nelemload1,*nload1,*nodeboun1,
  *nboun1,*ndirboun1,*ithermal1,*mi1,*ielmat1,*ielprop1,num_cpus,*neapar=NULL,
  *nebpar=NULL;

static double *co1,*vold1,*prop1;

void FORTRAN(createinum,(ITG *ipkon,ITG *inum,ITG *kon,char *lakon,ITG *nk,
			 ITG *ne,char *cflag,ITG *nelemload,ITG *nload,
			 ITG *nodeboun,
			 ITG *nboun,ITG *ndirboun,ITG *ithermal,double *co,
			 double *vold,ITG *mi,ITG *ielmat,ITG *ielprop,
			 double *prop)){

  ITG intpointvarm,calcul_fn,calcul_f,calcul_qa,calcul_cauchy,ikin,
    intpointvart,mt=mi[1]+1,i,j;

  /* determine which nodes are used (for output purposes) */
      
  /* variables for multithreading procedure */
    
  ITG sys_cpus,*ithread=NULL;
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

  if(*ne<num_cpus) num_cpus=*ne;
    
  pthread_t tid[num_cpus];
    
    /* determining the element bounds in each thread */

    NNEW(neapar,ITG,num_cpus);
    NNEW(nebpar,ITG,num_cpus);
    elementcpuload(neapar,nebpar,ne,ipkon,&num_cpus);

    ipkon1=ipkon;inum1=inum;kon1=kon;lakon1=lakon;nk1=nk;ne1=ne;
    cflag1=cflag;nelemload1=nelemload;nload1=nload;nodeboun1=nodeboun;
    nboun1=nboun;ndirboun1=ndirboun;ithermal1=ithermal;co1=co;
    vold1=vold;mi1=mi;ielmat1=ielmat;ielprop1=ielprop;prop1=prop;

    /* calculating the active nodes */
	
    if(((*nmethod!=4)&&(*nmethod!=5))||((iperturb[0]>1)&&(*mscalmethod<0))){
      printf(" Using up to %" ITGFORMAT " cpu(s) for the active node calculation.\n\n", num_cpus);
    }
	
    /* create threads and wait */
	
    NNEW(ithread,ITG,num_cpus);
    for(i=0; i<num_cpus; i++)  {
      ithread[i]=i;
      pthread_create(&tid[i], NULL, (void *)createinummt, (void *)&ithread[i]);
    }
    for(i=0; i<num_cpus; i++)
      pthread_join(tid[i], NULL);
    
    SFREE(ithread);SFREE(neapar);SFREE(nebpar);

  /* calculating the matrix system internal force vector */
    
    resultsforc(nk,f,fn,nactdof,ipompc,nodempc,
		coefmpc,labmpc,nmpc,mi,fmpc,&calcul_fn,&calcul_f,
		&num_cpus);
  
  return;

}

/* subroutine for multithreading of resultsmech */

void *createinummt(ITG *i){

  ITG indexfn,indexqa,indexnal,nea,neb,list1,*ilist1=NULL;

  indexfn=*i*mt1**nk1;
  indexqa=*i*4;
  indexnal=*i;

  nea=neapar[*i]+1;
  neb=nebpar[*i]+1;

  FORTRAN(createinum,(ipkon1,inum1,kon1,lakon1,nk1,ne1,cflag1,nelemload1,
		      nload1,nodeboun1,nboun1,ndirboun1,ithermal1,co1,
		      vold1,mi1,
  
  list1=0;
  FORTRAN(resultsmech,(co1,kon1,ipkon1,lakon1,ne1,v1,stx1,elcon1,nelcon1,
		       rhcon1,nrhcon1,alcon1,nalcon1,alzero1,ielmat1,ielorien1,
		       norien1,orab1,ntmat1_,t01,t11,ithermal1,prestr1,
		       iprestr1,eme1,iperturb1,&fn1[indexfn],iout1,
		       &qa1[indexqa],vold1,nmethod1,veold1,dtime1,time1,ttime1,
		       plicon1,nplicon1,plkcon1,nplkcon1,xstateini1,xstiff1,
		       xstate1,npmat1_,matname1,mi1,ielas1,icmd1,ncmat1_,
		       nstate1_,stiini1,vini1,ener1,eei1,enerini1,istep1,iinc1,
		       springarea1,reltime1,&calcul_fn1,&calcul_qa1,
		       &calcul_cauchy1,nener1,&ikin1,&nal[indexnal],ne01,
		       thicke1,emeini1,pslavsurf1,pmastsurf1,mortar1,clearini1,
		       &nea,&neb,ielprop1,prop1,kscale1,&list1,ilist1,smscale1,
		       mscalmethod1,&energysms1[indexnal],t0g1,t1g1,
		       islavquadel1,aut1,irowt1,jqt1,mortartrafoflag1,
		       intscheme1,physcon1));

  return NULL;
}

/* subroutine for multithreading of resultstherm */

void *resultsthermmt(ITG *i){

  ITG indexfn,indexqa,indexnal,nea,neb;

  indexfn=*i*mt1**nk1;
  indexqa=*i*4;
  indexnal=*i;

  nea=neapar[*i]+1;
  neb=nebpar[*i]+1;

  FORTRAN(resultstherm,(co1,kon1,ipkon1,lakon1,v1,elcon1,nelcon1,rhcon1,
			nrhcon1,ielmat1,ielorien1,norien1,orab1,ntmat1_,t01,
			iperturb1,&fn1[indexfn],shcon1,nshcon1,iout1,
			&qa1[indexqa],vold1,ipompc1,nodempc1,coefmpc1,nmpc1,
			dtime1,time1,ttime1,plkcon1,nplkcon1,xstateini1,
			xstiff1,xstate1,npmat1_,matname1,mi1,ncmat1_,nstate1_,
			cocon1,ncocon1,qfx1,ikmpc1,ilmpc1,istep1,iinc1,
			springarea1,&calcul_fn1,&calcul_qa1,&nal[indexnal],
			&nea,&neb,ithermal1,nelemload1,nload1,nmethod1,
			reltime1,sideload1,xload1,xloadold1,pslavsurf1,
			pmastsurf1,mortar1,clearini1,plicon1,nplicon1,ielprop1,
			prop1,iponoel1,inoel1,network1,ipobody1,xbody1,ibody1));

  return NULL;
}

/* subroutine for multithreading of calcenergy */

void *calcenergymt(ITG *i){

  ITG indexenergy,nea,neb;

  indexenergy=*i*4;

  nea=neapar[*i]+1;
  neb=nebpar[*i]+1;

  FORTRAN(calcenergy,(ipkon1,lakon1,kon1,co1,ener1,mi1,ne1,
		      thicke1,ielmat1,&energy1[indexenergy],
		      ielprop1,prop1,&nea,&neb));

  return NULL;
}
