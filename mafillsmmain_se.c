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

static char *lakon1,*sideload1,*matname1,*tieset1,*labmpc1,*set1;

static ITG *nk1,*kon1,*ipkon1,*ne1,*nodeboun1,*ndirboun1,*nboun1,
  *ipompc1,*nodempc1,*nmpc1,*nodeforc1,*ndirforc1,*nforc1,*nelemload1,
  *nload1,*ipobody1,*nbody1,*nactdof1,*neq1,
  *nmethod1=NULL,*ikmpc1,*ilmpc1,*ikboun1,*ilboun1,*nelcon1,
  *nrhcon1,*nalcon1,*ielmat1,*ielorien1,*norien1,*ntmat1_,*ithermal1,
  *iprestr1,*iperturb1,*iexpl1,*nplicon1,*nplkcon1,*npmat1_,
  *mi1,*ncmat1_,*mass1,*stiffness1,*buckling1,*rhsi1,*intscheme1,
  *nshcon1,*ncocon1,*istep1,*iinc1,*coriolis1,*ibody1,*nstate1_,
  *integerglob1,*istartset1,*iendset1,*ialset1,*ntie1,*nasym1,
  *mortar1,*ielprop1,*ne01,num_cpus,*ndesi1,*nodedesi1,
  *nzss1,*jqs1,*irows1,*icoordinate1,*istartelem1,*ialelem1,
  *cyclicsymmetry1,*ics1,*mcs1,*ieigenfrequency1,*neapar=NULL,
  *nebpar=NULL,*nset1;

static double *co1,*xboun1,*coefmpc1,*xforc1,*xload1,*xbody1,*cgr1,
    *elcon1,*rhcon1,*alcon1,*alzero1,*xdesi1,*v1,*sigma1,
    *orab1,*t01,*t11,*prestr1,*vold1,*sti1,*stx1,*cs1,
    *plicon1,*plkcon1,*xstiff1,*dtime1,*physcon1,*shcon1,*cocon1,
    *ttime1,*time1,*xloadold1,*reltime1,*veold1,*springarea1,
    *xstateini1,*xstate1,*thicke1,*doubleglob1,*pslavsurf1,*pmastsurf1,
  *clearini1,*prop1,*distmin1,*df1,*dfl1,*dxstiff1,*sigmak1;

void mafillsmmain_se(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	       ITG *ne,ITG *nodeboun,ITG *ndirboun,double *xboun, 
	       ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc, 
	       ITG *nmpc,ITG *nodeforc,ITG *ndirforc,
	       double *xforc,ITG *nforc,ITG *nelemload,char *sideload,
	       double *xload,ITG *nload,double *xbody,ITG *ipobody,
	       ITG *nbody,double *cgr,ITG *nactdof, 
	       ITG *neq,ITG *nmethod,ITG *ikmpc,ITG *ilmpc,
               ITG *ikboun, ITG *ilboun,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,ITG *ithermal,
	       double *prestr,ITG *iprestr,double *vold,
	       ITG *iperturb,double *sti,double *stx,
	       ITG *iexpl,
               double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
               double *xstiff, 
	       ITG *npmat_,double *dtime,char *matname,ITG *mi,
               ITG *ncmat_,ITG *mass,ITG *stiffness,ITG *buckling,ITG *rhsi,
               ITG *intscheme,double *physcon,double *shcon,ITG *nshcon,
               double *cocon,ITG *ncocon,double *ttime,double *time,
               ITG *istep,ITG *iinc,ITG *coriolis,ITG *ibody,
	       double *xloadold,double *reltime,double *veold,
               double *springarea,ITG *nstate_,double *xstateini,
	       double *xstate,double *thicke,
               ITG *integerglob,double *doubleglob,char *tieset,
	       ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
	       ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
	       double *clearini,ITG *ielprop,double *prop,ITG *ne0,
               double *fnext,double *distmin,ITG *ndesi,ITG *nodedesi,
	       double *df,ITG *nzss,ITG *jqs,ITG *irows,
	       ITG *icoordinate,double *dxstiff,double *xdesi,
	       ITG *istartelem,ITG *ialelem,double *v,double *sigma,
	       ITG *cyclicsymmetry,char *labmpc,ITG *ics,double *cs,
	       ITG *mcs,ITG *ieigenfrequency,char *set,ITG *nset,
	       double *sigmak){
	       
    ITG i,j;
     
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

    envloc = getenv("CCX_NPROC_SENS");
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

    /* allocating fields for mass and stiffness matrix and sensitivity matrix */
    
    /* maximum 20 design variables per element (20-node element) */

    if(!*cyclicsymmetry){
	NNEW(dfl1,double,num_cpus*60*20);
	NNEW(df1,double,num_cpus**nzss);
    }else{
	NNEW(dfl1,double,num_cpus*120*20);
	NNEW(df1,double,num_cpus*2**nzss);
    }

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
    nodeboun1=nodeboun;ndirboun1=ndirboun;xboun1=xboun;
    nboun1=nboun;ipompc1=ipompc;nodempc1=nodempc;coefmpc1=coefmpc;
    nmpc1=nmpc;nodeforc1=nodeforc;ndirforc1=ndirforc;xforc1=xforc;
    nforc1=nforc;nelemload1=nelemload;sideload1=sideload;xload1=xload;
    nload1=nload;xbody1=xbody;ipobody1=ipobody;nbody1=nbody;
    cgr1=cgr;nactdof1=nactdof;neq1=neq;
    ikmpc1=ikmpc;ilmpc1=ilmpc;ikboun1=ikboun;
    ilboun1=ilboun;elcon1=elcon;nelcon1=nelcon;rhcon1=rhcon;
    nrhcon1=nrhcon;alcon1=alcon;nalcon1=nalcon;alzero1=alzero;
    ielmat1=ielmat;ielorien1=ielorien;norien1=norien;orab1=orab;
    ntmat1_=ntmat_;t01=t0;t11=t1;ithermal1=ithermal;prestr1=prestr;
    iprestr1=iprestr;vold1=vold;iperturb1=iperturb;sti1=sti;
    stx1=stx;iexpl1=iexpl;plicon1=plicon;nplicon1=nplicon;
    plkcon1=plkcon;nplkcon1=nplkcon;xstiff1=xstiff;npmat1_=npmat_;
    dtime1=dtime;matname1=matname;mi1=mi;ncmat1_=ncmat_;mass1=mass;
    stiffness1=stiffness;buckling1=buckling;rhsi1=rhsi;intscheme1=intscheme;
    physcon1=physcon;shcon1=shcon;nshcon1=nshcon;cocon1=cocon;
    ncocon1=ncocon;ttime1=ttime;time1=time;istep1=istep;iinc1=iinc;
    coriolis1=coriolis;ibody1=ibody;xloadold1=xloadold;reltime1=reltime;
    veold1=veold;springarea1=springarea;nstate1_=nstate_;xstateini1=xstateini;
    xstate1=xstate;thicke1=thicke;integerglob1=integerglob;
    doubleglob1=doubleglob;tieset1=tieset;istartset1=istartset;
    iendset1=iendset;ialset1=ialset;ntie1=ntie;nasym1=nasym;
    pslavsurf1=pslavsurf;pmastsurf1=pmastsurf;mortar1=mortar;
    clearini1=clearini;ielprop1=ielprop;prop1=prop;ne01=ne0;
    distmin1=distmin;ndesi1=ndesi;nodedesi1=nodedesi;v1=v;
    nzss1=nzss;jqs1=jqs;irows1=irows;icoordinate1=icoordinate;
    dxstiff1=dxstiff;xdesi1=xdesi;istartelem1=istartelem;ialelem1=ialelem;
    sigma1=sigma;cyclicsymmetry1=cyclicsymmetry;labmpc1=labmpc;
    ics1=ics;cs1=cs;mcs1=mcs;ieigenfrequency1=ieigenfrequency;
    set1=set;nset1=nset;sigmak1=sigmak;

    /* calculating the stiffness/mass sensitivity */
    
    printf(" Using up to %" ITGFORMAT " cpu(s) for the calculation of the sensitivity of the external forces \n and/or the element stiffness matrices.\n\n", num_cpus);
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,num_cpus);
    for(i=0; i<num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)mafillsmsemt, (void *)&ithread[i]);
    }
    for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);SFREE(neapar);SFREE(nebpar);

    for(j=0;j<num_cpus;j++){
	if(nmethod1[j]==0){
	    *nmethod=0;
	    break;
	}
    }
    SFREE(nmethod1);
    
    /* passing of df */

    if(!*cyclicsymmetry){
//	if((iperturb[1]==1)&&(*ieigenfrequency!=1)){
	if(*ieigenfrequency!=1){

            /* nonlinear geometric: add df to df from results_se */

	    for(i=0;i<*nzss;i++){
		df[i]+=df1[i];
	    }
	}else{

	    for(i=0;i<*nzss;i++){
		df[i]=df1[i];
	    }
	}
	for(i=0;i<*nzss;i++){
	    for(j=1;j<num_cpus;j++){
		df[i]+=df1[i+j**nzss];
	    }
	}
    }else{
	for(i=0;i<2**nzss;i++){
	    df[i]=df1[i];
	}
	for(i=0;i<2**nzss;i++){
	    for(j=1;j<num_cpus;j++){
		df[i]+=df1[i+j*2**nzss];
	    }
	}
    }
    SFREE(df1);SFREE(dfl1);
        
        return;

}

/* subroutine for multithreading of mafillsm */

void *mafillsmsemt(ITG *i){

    ITG indexdf,indexdfl,nea,neb;

    if(!*cyclicsymmetry1){
	indexdf=*i**nzss1;
	indexdfl=*i*60*20;
    }else{
	indexdf=*i*2**nzss1;
	indexdfl=*i*120*20;
    }

    nea=neapar[*i]+1;
    neb=nebpar[*i]+1;

    if(!*cyclicsymmetry1){
	FORTRAN(mafillsmse,(co1,kon1,ipkon1,lakon1,ne1,ipompc1,nodempc1,
	    coefmpc1,nmpc1,nelemload1,sideload1,xload1,nload1,xbody1,
	    ipobody1,nbody1,cgr1,nactdof1,neq1,&nmethod1[*i],ikmpc1,
	    ilmpc1,elcon1,nelcon1,rhcon1,nrhcon1,alcon1,nalcon1,alzero1,
	    ielmat1,ielorien1,norien1,orab1,ntmat1_,t01,t11,ithermal1,
	    iprestr1,vold1,iperturb1,sti1,stx1,iexpl1,plicon1,nplicon1,
            plkcon1,nplkcon1,xstiff1,npmat1_,dtime1,matname1,mi1,
            ncmat1_,mass1,stiffness1,buckling1,rhsi1,intscheme1,physcon1,
            ttime1,time1,istep1,iinc1,coriolis1,ibody1,xloadold1,
	    reltime1,veold1,springarea1,nstate1_,xstateini1,xstate1,
            thicke1,integerglob1,doubleglob1,tieset1,istartset1,iendset1,
	    ialset1,ntie1,nasym1,pslavsurf1,pmastsurf1,mortar1,clearini1,
	    ielprop1,prop1,ne01,&nea,&neb,distmin1,ndesi1,nodedesi1,
	    &df1[indexdf],jqs1,irows1,&dfl1[indexdfl],
	    icoordinate1,dxstiff1,xdesi1,istartelem1,ialelem1,v1,sigma1,
	    ieigenfrequency1,set1,nset1,sigmak1));
    }else{
	FORTRAN(mafillsmcsse,(co1,kon1,ipkon1,lakon1,ne1,ipompc1,nodempc1,
	    coefmpc1,nmpc1,nelemload1,sideload1,xload1,nload1,xbody1,
	    ipobody1,nbody1,cgr1,nactdof1,neq1,&nmethod1[*i],ikmpc1,
	    ilmpc1,elcon1,nelcon1,rhcon1,nrhcon1,alcon1,nalcon1,alzero1,
	    ielmat1,ielorien1,norien1,orab1,ntmat1_,t01,t11,ithermal1,
	    iprestr1,vold1,iperturb1,sti1,stx1,iexpl1,plicon1,nplicon1,
            plkcon1,nplkcon1,xstiff1,npmat1_,dtime1,matname1,mi1,
            ncmat1_,mass1,stiffness1,buckling1,rhsi1,intscheme1,physcon1,
            ttime1,time1,istep1,iinc1,coriolis1,ibody1,xloadold1,
	    reltime1,veold1,springarea1,nstate1_,xstateini1,xstate1,
            thicke1,integerglob1,doubleglob1,tieset1,istartset1,iendset1,
	    ialset1,ntie1,nasym1,pslavsurf1,pmastsurf1,mortar1,clearini1,
	    ielprop1,prop1,ne01,&nea,&neb,distmin1,ndesi1,nodedesi1,
	    &df1[indexdf],jqs1,irows1,&dfl1[indexdfl],
	    icoordinate1,dxstiff1,xdesi1,istartelem1,ialelem1,v1,sigma1,
	    labmpc1,ics1,cs1,mcs1,nk1,nzss1,set1,nset1));
    }
	
	    
    return NULL;
}
