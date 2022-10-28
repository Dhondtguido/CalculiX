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

static char *lakon1,*matname1,*sideload1;

static ITG *nk1,*kon1,*ipkon1,*ne1,*ipompc1,*nodempc1,
  *nmpc1,*nodeforc1,*ndirforc1,*nforc1,*nelemload1,*nload1,
  *ipobody1,*nbody1,*nactdof1,*neq1,*nmethod1=NULL,*ikmpc1,*ilmpc1,
  *nelcon1,*nrhcon1,*nalcon1,*ielmat1,*ielorien1,*norien1,
  *ntmat_1,*ithermal1,*iprestr1,*iperturb1,*iexpl1,*nplicon1,
  *nplkcon1,*npmat_1,*istep1,*iinc1,*ibody1,*mi1,
  *ikactmech1,*nactmech1,*ielprop1,*nstate_1,*ntrans1,*inotr1,
  num_cpus,*neapar=NULL,*nebpar=NULL;

static double *co1,*coefmpc1,*xforcact1,*xloadact1,
  *xbodyact1,*xbodyact1,*cgr1,*elcon1,*rhcon1,*alcon1,
  *alzero1,*orab1,*t01,*t1act1,*vold1,*plicon1,*plkcon1,
  *ttime1,*time1,*dtime1,*physcon1,*xbodyold1,*reltime1,*veold1,
  *prop1,*sti1,*xstateini1,*xstate1,*trab1,*fext1=NULL,*fnext1=NULL;

void rhsmain(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
	     ITG *ipompc,ITG *nodempc,double *coefmpc,ITG *nmpc,
	     ITG *nodeforc,ITG *ndirforc,double *xforcact,ITG *nforc,
	     ITG *nelemload,char *sideload,double *xloadact,ITG *nload,
	     double *xbodyact,ITG *ipobody,ITG *nbody,double *cgr,double *fext,
	     ITG *nactdof,ITG *neq,ITG *nmethod,ITG *ikmpc,ITG *ilmpc,
	     double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,double *alcon,
	     ITG *nalcon,double *alzero,ITG *ielmat,ITG *ielorien,ITG *norien,
	     double *orab,ITG *ntmat_,double *t0,double *t1act,ITG *ithermal,
	     ITG *iprestr,double *vold,ITG *iperturb,ITG *iexpl,double *plicon,
	     ITG *nplicon,double *plkcon,ITG *nplkcon,ITG *npmat_,
	     double *ttime,double *time,ITG *istep,ITG *iinc,double *dtime,
	     double *physcon,ITG *ibody,double *xbodyold,double *reltime,
	     double *veold,char *matname,ITG *mi,ITG *ikactmech,ITG *nactmech,
	     ITG *ielprop,double *prop,double *sti,double *xstateini,
	     double *xstate,ITG *nstate_,ITG *ntrans,ITG *inotr,double *trab,
	     double *fnext){

  ITG sys_cpus,*ithread=NULL,i,j,isum,idelta,mt=mi[1]+1;
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

  // next line is to be inserted in a similar way for all other paralell parts

  if(*ne<num_cpus) num_cpus=*ne;

  pthread_t tid[num_cpus];



  /* determining the element bounds in each thread */

  NNEW(neapar,ITG,num_cpus);
  NNEW(nebpar,ITG,num_cpus);
  elementcpuload(neapar,nebpar,ne,ipkon,&num_cpus);

  /* allocating memory for the nodal external forces (only for explicit
     dynamic calculations) */

  if((*nmethod==4)&&(iperturb[0]>1)){
    NNEW(fnext1,double,num_cpus*mt**nk);
  }

  /* allocating memory for nmethod; if the Jacobian determinant
     in any of the elements is nonpositive, nmethod is set to
     zero */

  NNEW(nmethod1,ITG,num_cpus);
  for(j=0;j<num_cpus;j++){
    nmethod1[j]=*nmethod;
  }
	
  nk1=nk;co1=co;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;
  ipompc1=ipompc;nodempc1=nodempc;coefmpc1=coefmpc;nmpc1=nmpc;
  nodeforc1=nodeforc;
  ndirforc1=ndirforc;xforcact1=xforcact;nforc1=nforc;nelemload1=nelemload;
  sideload1=sideload;xloadact1=xloadact;nload1=nload;xbodyact1=xbodyact;
  ipobody1=ipobody;nbody1=nbody;cgr1=cgr;nactdof1=nactdof;neq1=neq;
  ikmpc1=ikmpc;ilmpc1=ilmpc;elcon1=elcon;nelcon1=nelcon;
  rhcon1=rhcon;nrhcon1=nrhcon;alcon1=alcon;nalcon1=nalcon;alzero1=alzero;
  ielmat1=ielmat;ielorien1=ielorien;norien1=norien;orab1=orab;ntmat_1=ntmat_;
  t01=t0;t1act1=t1act;ithermal1=ithermal;iprestr1=iprestr;vold1=vold;
  iperturb1=iperturb;
  iexpl1=iexpl;plicon1=plicon;nplicon1=nplicon;plkcon1=plkcon;nplkcon1=nplkcon;
  npmat_1=npmat_;ttime1=ttime;time1=time;istep1=istep;iinc1=iinc;dtime1=dtime;
  physcon1=physcon;ibody1=ibody;xbodyold1=xbodyold;reltime1=reltime;
  veold1=veold;matname1=matname;mi1=mi;ielprop1=ielprop;
  prop1=prop;sti1=sti;xstateini1=xstateini;xstate1=xstate;nstate_1=nstate_;
  ntrans1=ntrans;inotr1=inotr;trab1=trab;

  NNEW(fext1,double,num_cpus**neq);
  NNEW(ikactmech1,ITG,num_cpus**neq);
  NNEW(nactmech1,ITG,num_cpus);
  
  /* calculating the rhs */

  if(((*nmethod!=4)&&(*nmethod!=5))||((iperturb[0]>1)&&(*iexpl<=1))){
    printf(" Using up to %" ITGFORMAT " cpu(s) for the rhs calculation.\n\n", num_cpus);
  }

  /* create threads and wait */
  /* distributed force: body foce*/
	
  NNEW(ithread,ITG,num_cpus);
  for(i=0; i<num_cpus; i++)  {
    ithread[i]=i;
    pthread_create(&tid[i], NULL, (void *)rhsmt, (void *)&ithread[i]);
  }
  for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
	

  for(i=0;i<*neq;i++){
    fext[i]=fext1[i];
  }
  for(i=0;i<*neq;i++){
    for(j=1;j<num_cpus;j++){
      fext[i]+=fext1[i+j**neq];
    }
  }

  /* collecting the contributions of fnext; only for explicit
     dynamic calculations */
  
  if((*nmethod==4)&&(iperturb[0]>1)){
    for(i=0;i<mt**nk;i++){
      fnext[i]=fnext1[i];
    }
    for(i=0;i<mt**nk;i++){
      for(j=1;j<num_cpus;j++){
	fnext[i]+=fnext1[i+j*mt**nk];
      }
    }
    SFREE(fnext1);
  }

  for(j=0;j<num_cpus;j++){
    if(nmethod1[j]==0){
      *nmethod=0;
      break;
    }
  }

  SFREE(ithread);SFREE(neapar);SFREE(nebpar);
  SFREE(fext1);SFREE(nmethod1);

  /* merging ikactmech1 into ikactmech */

  FORTRAN(merge_ikactmech1,(ikactmech1,nactmech1,neq,ikactmech,nactmech,
			    &num_cpus));
  
  SFREE(ikactmech1);SFREE(nactmech1);
	
  FORTRAN(rhsnodef,(co,kon,ne,
		    ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
		    nforc,fext,nactdof,nmethod,ikmpc,ntmat_,iperturb,
		    mi,ikactmech,nactmech,ntrans,inotr,trab,fnext));
	
}



/* subroutine for multithreading of rhs: distributed forces */

void *rhsmt(ITG *i){

  ITG indexf,indexfnext,nea,neb;

  indexf=*i**neq1;
  indexfnext=0;
  if((nmethod1[0]==4)&&(iperturb1[0]>1)){
    indexfnext=*i*(mi1[1]+1)**nk1;
  }

  nea=neapar[*i]+1;
  neb=nebpar[*i]+1;

  FORTRAN(rhs,(co1,nk1,kon1,ipkon1,lakon1,ne1,ipompc1,nodempc1,
	       coefmpc1,nmpc1,nodeforc1,ndirforc1,xforcact1,
	       nforc1,nelemload1,sideload1,xloadact1,nload1,xbodyact1,ipobody1,
	       nbody1,cgr1,&fext1[indexf],nactdof1,neq1,
	       &nmethod1[*i],ikmpc1,ilmpc1,
	       elcon1,nelcon1,rhcon1,nrhcon1,alcon1,nalcon1,alzero1,
	       ielmat1,ielorien1,norien1,orab1,ntmat_1,
	       t01,t1act1,ithermal1,iprestr1,vold1,iperturb1,
	       iexpl1,plicon1,nplicon1,plkcon1,nplkcon1,
	       npmat_1,ttime1,time1,istep1,iinc1,dtime1,physcon1,ibody1,
	       xbodyold1,reltime1,veold1,matname1,mi1,&ikactmech1[indexf],
	       &nactmech1[*i],ielprop1,prop1,sti1,xstateini1,xstate1,nstate_1,
	       ntrans1,inotr1,trab1,&fnext1[indexfnext],&nea,&neb));

  return NULL;
}

