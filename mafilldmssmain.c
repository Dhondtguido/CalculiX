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

static char *lakon1,*sideload1,*matname1,*tieset1,*set1;

static ITG *nk1,*kon1,*ipkon1,*ne1,*nset1,
    *ipompc1,*nodempc1,*nmpc1,*nelemload1,
    *nload1,*ipobody1,*nbody1,*nactdof1,*jq1,*irow1,*neq1,
    *nmethod1=NULL,*ikmpc1,*ilmpc1,*nelcon1,
    *nrhcon1,*nalcon1,*ielmat1,*ielorien1,*norien1,*ntmat1_,*ithermal1,
    *iperturb1,*iexpl1,*nplicon1,*nplkcon1,*npmat1_,*mi1,*ncmat1_,
    *istep1,*iinc1,*ibody1,*nstate1_,
    *integerglob1,*istartset1,*iendset1,*ialset1,*ntie1,*nasym1,
    *mortar1,*ielprop1,*ne01,num_cpus,
    *neapar=NULL,*nebpar=NULL,*ndamp1,*nzs1;

static double *co1,*coefmpc1,*xload1,*xbody1,*cgr1,
    *ad1=NULL,*au1=NULL,*elcon1,*rhcon1,*alcon1,*alzero1,
    *orab1,*t01,*t11,*vold1,*sti1,*stx1,
    *plicon1,*plkcon1,*xstiff1,*dtime1,*physcon1,
    *ttime1,*time1,*xloadold1,*reltime1,*veold1,*springarea1,
    *xstateini1,*xstate1,*thicke1,*doubleglob1,*pslavsurf1,*pmastsurf1,
    *clearini1,*prop1,*freq1,*dacon1;

void mafilldmssmain(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	       ITG *ne,ITG *ipompc,ITG *nodempc,double *coefmpc, 
	       ITG *nmpc,ITG *nelemload,char *sideload,
	       double *xload,ITG *nload,double *xbody,ITG *ipobody,
	       ITG *nbody,double *cgr,
	       double *ad,double *au,ITG *nactdof, 
	       ITG *jq,ITG *irow,ITG *neq,
	       ITG *nmethod,ITG *ikmpc,ITG *ilmpc, 
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,ITG *ithermal,double *vold,
	       ITG *iperturb,double *sti,ITG *nzs,double *stx,
	       ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff, 
	       ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,double *physcon,double *ttime,double *time,
               ITG *istep,ITG *iinc,ITG *ibody,
	       double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
	       double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
	       ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
	       ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
	       double *clearini,ITG *ielprop,double *prop,ITG *ne0,
	       double *freq,ITG *ndamp,double *dacon,char *set,ITG *nset){

    ITG i,j,mt=mi[1]+1;
      
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

//    sys_cpus=1;

    /* automatic detection of available number of processors */

    if(sys_cpus==0){
	sys_cpus = getSystemCPUs();
	if(sys_cpus<1) sys_cpus=1;
    }

    /* local declaration prevails, if strictly positive */

    envloc = getenv("CCX_NPROC_STIFFNESS");
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
    
    /* determining the element bounds in each thread
       contact elements are included (because of contact
       damping) */
    
    NNEW(neapar,ITG,num_cpus);
    NNEW(nebpar,ITG,num_cpus);

    elementcpuload(neapar,nebpar,ne,ipkon,&num_cpus);

    /* allocating fields for mass and stiffness matrix */

    NNEW(ad1,double,num_cpus*neq[1]);
    NNEW(au1,double,(long long)num_cpus*nzs[1]);

    /* allocating memory for nmethod; if the Jacobian determinant
       in any of the elements is nonpositive, nmethod is set to
       zero */

    NNEW(nmethod1,ITG,num_cpus);
    for(j=0;j<num_cpus;j++){
	nmethod1[j]=*nmethod;
    }

    /* calculating the stiffness and/or mass matrix 
       (symmetric part) */

    co1=co;nk1=nk;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;
    ipompc1=ipompc;nodempc1=nodempc;coefmpc1=coefmpc;
    nmpc1=nmpc;
    nelemload1=nelemload;sideload1=sideload;xload1=xload;
    nload1=nload;xbody1=xbody;ipobody1=ipobody;nbody1=nbody;
    cgr1=cgr;nactdof1=nactdof;jq1=jq;irow1=irow;neq1=neq;
    ikmpc1=ikmpc;ilmpc1=ilmpc;
    elcon1=elcon;nelcon1=nelcon;rhcon1=rhcon;
    nrhcon1=nrhcon;alcon1=alcon;nalcon1=nalcon;alzero1=alzero;
    ielmat1=ielmat;ielorien1=ielorien;norien1=norien;orab1=orab;
    ntmat1_=ntmat_;t01=t0;t11=t1;ithermal1=ithermal;
    vold1=vold;iperturb1=iperturb;sti1=sti;
    stx1=stx;iexpl1=iexpl;plicon1=plicon;nplicon1=nplicon;
    plkcon1=plkcon;nplkcon1=nplkcon;xstiff1=xstiff;npmat1_=npmat_;
    dtime1=dtime;matname1=matname;mi1=mi;ncmat1_=ncmat_;
    physcon1=physcon;
    ttime1=ttime;time1=time;istep1=istep;iinc1=iinc;
    ibody1=ibody;xloadold1=xloadold;reltime1=reltime;
    veold1=veold;springarea1=springarea;nstate1_=nstate_;xstateini1=xstateini;
    xstate1=xstate;thicke1=thicke;integerglob1=integerglob;
    doubleglob1=doubleglob;tieset1=tieset;istartset1=istartset;
    iendset1=iendset;ialset1=ialset;ntie1=ntie;nasym1=nasym;
    pslavsurf1=pslavsurf;pmastsurf1=pmastsurf;mortar1=mortar;
    clearini1=clearini;ielprop1=ielprop;prop1=prop;ne01=ne0;
    freq1=freq;ndamp1=ndamp;dacon1=dacon;nzs1=nzs;set1=set;
    nset1=nset;

    /* calculating the stiffness/mass */
    
    printf(" Using up to %" ITGFORMAT " cpu(s) for the symmetric stiffness/mass contributions.\n\n", num_cpus);
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,num_cpus);
    for(i=0; i<num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)mafilldmssmt, (void *)&ithread[i]);
    }
    for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);SFREE(neapar);SFREE(nebpar);

    /* copying and accumulating the damping matrix */

    for(i=0;i<neq[1];i++){
	ad[i]=ad1[i];
    }

    for(i=0;i<neq[1];i++){
	for(j=1;j<num_cpus;j++){
	    ad[i]+=ad1[i+j*neq[1]];
	}
    }
    SFREE(ad1);

    for(i=0;i<nzs[1];i++){
	au[i]=au1[i];
    }

    for(i=0;i<nzs[1];i++){
	for(j=1;j<num_cpus;j++){
	    au[i]+=au1[i+(long long)j*nzs[1]];
	}
    }
    SFREE(au1);

    for(j=0;j<num_cpus;j++){
	if(nmethod1[j]==0){
	    *nmethod=0;
	    break;
	}
    }
    SFREE(nmethod1);
  
  return;

}

/* subroutine for multithreading of mafilldmss */

void *mafilldmssmt(ITG *i){

    ITG indexad,nea,neb;
    long long indexau;

    indexad=0;
    indexau=0;

    indexad=*i*neq1[1];
    indexau=(long long)*i*nzs1[1];

    nea=neapar[*i]+1;
    neb=nebpar[*i]+1;

    FORTRAN(mafilldmss,(co1,nk1,kon1,ipkon1,lakon1,ne1,
	    ipompc1,nodempc1,coefmpc1,nmpc1,
	    nelemload1,sideload1,xload1,nload1,xbody1,ipobody1,
	    nbody1,cgr1,&ad1[indexad],&au1[indexau],
	    nactdof1,jq1,irow1,neq1,&nmethod1[*i],
	    ikmpc1,ilmpc1,
	    elcon1,nelcon1,rhcon1,nrhcon1,alcon1,nalcon1,alzero1,ielmat1,
	    ielorien1,norien1,orab1,ntmat1_,
	    t01,t11,ithermal1,vold1,iperturb1,sti1,
	    stx1,iexpl1,plicon1,
            nplicon1,plkcon1,nplkcon1,
	    xstiff1,npmat1_,dtime1,matname1,mi1,
            ncmat1_,physcon1,
            ttime1,time1,istep1,iinc1,
	    ibody1,xloadold1,reltime1,veold1,springarea1,nstate1_,
            xstateini1,xstate1,thicke1,integerglob1,doubleglob1,
	    tieset1,istartset1,iendset1,ialset1,ntie1,nasym1,pslavsurf1,
	    pmastsurf1,mortar1,clearini1,ielprop1,prop1,ne01,
	    &nea,&neb,freq1,ndamp1,dacon1,set1,nset1));

    return NULL;
}
