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

static ITG *nk1,*kon1,*ipkon1,*ne1,*nodeboun1,*ndirboun1,*nboun1,*ipompc1,
  *nodempc1,*nmpc1,*nodeforc1,*ndirforc1,*nforc1,*nelemload1,*nload1,*ipobody1,
  *nbody1,*nactdof1,*icol1,*jq1,*irow1,*neq1,*nzl1,*nmethod1=NULL,*ikmpc1,
  *ilmpc1,*ikboun1,*ilboun1,*nelcon1,*nrhcon1,*nalcon1,*ielmat1,*ielorien1,
  *norien1,*ntmat1_,*ithermal1,*iprestr1,*iperturb1,*nzs1,*iexpl1,*nplicon1,
  *nplkcon1,*npmat1_,*mi1,*ncmat1_,*mass1,*stiffness1,*buckling1,*rhsi1,
  *intscheme1,*nshcon1,*ncocon1,*istep1,*iinc1,*coriolis1,*ibody1,*nstate1_,
  *integerglob1,*istartset1,*iendset1,*ialset1,*ntie1,*nasym1,*mortar1,
  *ielprop1,*ne01,*kscale1,*iponoel1,*inoel1,*network1,*neaparm=NULL,
  *nebparm=NULL,*neapart=NULL,*nebpart=NULL,*nset1;

static double *co1,*xboun1,*coefmpc1,*xforc1,*xload1,*xbody1,*cgr1,*ad1=NULL,
  *au1=NULL,*bb1,*elcon1,*rhcon1,*alcon1,*alzero1,*orab1,*t01,*t11,*prestr1,
  *vold1,*sti1,*stx1,*adb1,*aub1,*plicon1,*plkcon1,*xstiff1,*dtime1,*physcon1,
  *shcon1,*cocon1,*ttime1,*time1,*xloadold1,*reltime1,*veold1,*springarea1,
  *xstateini1,*xstate1,*thicke1,*doubleglob1,*pslavsurf1,*pmastsurf1,
  *clearini1,*prop1;
  
void mafillsmasmain(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
		ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
		ITG *ipompc,ITG *nodempc,double *coefmpc,ITG *nmpc,
		ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
		ITG *nelemload,char *sideload,double *xload,ITG *nload,
		double *xbody,ITG *ipobody,ITG *nbody,double *cgr,double *ad,
		double *au,double *bb,ITG *nactdof,ITG *icol,ITG *jq,
		ITG *irow,ITG *neq,ITG *nzl,ITG *nmethod,ITG *ikmpc,
		ITG *ilmpc,ITG *ikboun,ITG *ilboun,double *elcon,ITG *nelcon,
		double *rhcon,ITG *nrhcon,double *alcon,ITG *nalcon,
		double *alzero,ITG *ielmat,ITG *ielorien,ITG *norien,
		double *orab,ITG *ntmat_,double *t0,double *t1,ITG *ithermal,
		double *prestr,ITG *iprestr,double *vold,ITG *iperturb,
		double *sti,ITG *nzs,double *stx,double *adb,double *aub,
		ITG *iexpl,double *plicon,ITG *nplicon,double *plkcon,
		ITG *nplkcon,double *xstiff,ITG *npmat_,double *dtime,
		char *matname,ITG *mi,ITG *ncmat_,ITG *mass,ITG *stiffness,
		ITG *buckling,ITG *rhsi,ITG *intscheme,double *physcon,
		double *shcon,ITG *nshcon,double *cocon,ITG *ncocon,
		double *ttime,double *time,ITG *istep,ITG *iinc,
		ITG *coriolis,ITG *ibody,double *xloadold,double *reltime,
		double *veold,double *springarea,ITG *nstate_,
		double *xstateini,double *xstate,double *thicke,
		ITG *integerglob,double *doubleglob,char *tieset,
		ITG *istartset,ITG *iendset,ITG *ialset,ITG *ntie,
		ITG *nasym,double *pslavsurf,double *pmastsurf,ITG *mortar,
		double *clearini,ITG *ielprop,double *prop,ITG *ne0,
		ITG *kscale,ITG *iponoel,ITG *inoel,ITG *network,
		char *set,ITG *nset){

  ITG i,j,isiz,nec;
      
  /* variables for multithreading procedure */
    
  ITG sys_cpus,*ithread=NULL,num_cpus;
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

  /* storing the symmetric matrix in asymmetric format */

  isiz=nzs[1];cpypardou(&au[nzs[2]],au,&isiz,&num_cpus);

  /* dynamic applications */

  if((*nmethod==2)||(*nmethod==4)){
    isiz=nzs[1];cpypardou(&aub[nzs[2]],aub,&isiz,&num_cpus);
  }
    
  /* number of contact elements */
    
  nec=*ne-*ne0;
    
  pthread_t tid[num_cpus];
    
  /* determining the element bounds in each thread 
     asymmetric mechanical contributions*/
    
  NNEW(neaparm,ITG,num_cpus);
  NNEW(nebparm,ITG,num_cpus);

  if(nec>0){

    /* contact elements */
    
    elementcpuload(neaparm,nebparm,&nec,&ipkon[*ne0],&num_cpus);

    for(i=0;i<num_cpus;i++){
      neaparm[i]+=(*ne0);
      nebparm[i]+=(*ne0);
    }
  }else{

    /* no contact elements */
    
    for(i=0;i<num_cpus;i++){
      neaparm[i]=(*ne0)+1;
      nebparm[i]=(*ne0);
    }
  }
    
  /* determining the element bounds in each thread 
     asymmetric thermal contributions*/
    
  NNEW(neapart,ITG,num_cpus);
  NNEW(nebpart,ITG,num_cpus);
    
  elementcpuload(neapart,nebpart,ne,ipkon,&num_cpus);

  /* allocating fields for mass and stiffness matrix */

  NNEW(ad1,double,num_cpus*neq[1]);
  NNEW(au1,double,(long long)num_cpus*(nzs[2]+nzs[1]));

  /* allocating memory for nmethod; if the Jacobian determinant
     in any of the elements is nonpositive, nmethod is set to
     zero */

  NNEW(nmethod1,ITG,num_cpus);
  for(j=0;j<num_cpus;j++){
    nmethod1[j]=*nmethod;
  }

  /* calculating the stiffness 
     (asymmetric part) */

  co1=co;nk1=nk;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;nodeboun1=nodeboun;
  ndirboun1=ndirboun;xboun1=xboun;nboun1=nboun;ipompc1=ipompc;nodempc1=nodempc;
  coefmpc1=coefmpc;nmpc1=nmpc;nodeforc1=nodeforc;ndirforc1=ndirforc;
  xforc1=xforc;nforc1=nforc;nelemload1=nelemload;sideload1=sideload;
  xload1=xload;nload1=nload;xbody1=xbody;ipobody1=ipobody;nbody1=nbody;
  cgr1=cgr;bb1=bb;nactdof1=nactdof;icol1=icol;jq1=jq;irow1=irow;
  neq1=neq;nzl1=nzl;ikmpc1=ikmpc;ilmpc1=ilmpc;ikboun1=ikboun;
  ilboun1=ilboun;elcon1=elcon;nelcon1=nelcon;rhcon1=rhcon;nrhcon1=nrhcon;
  alcon1=alcon;nalcon1=nalcon;alzero1=alzero;ielmat1=ielmat;ielorien1=ielorien;
  norien1=norien;orab1=orab;ntmat1_=ntmat_;t01=t0;t11=t1;ithermal1=ithermal;
  prestr1=prestr;iprestr1=iprestr;vold1=vold;iperturb1=iperturb;sti1=sti;
  nzs1=nzs;stx1=stx;adb1=adb;aub1=aub;iexpl1=iexpl;plicon1=plicon;
  nplicon1=nplicon;plkcon1=plkcon;nplkcon1=nplkcon;xstiff1=xstiff;
  npmat1_=npmat_;dtime1=dtime;matname1=matname;mi1=mi;ncmat1_=ncmat_;mass1=mass;
  stiffness1=stiffness;buckling1=buckling;rhsi1=rhsi;intscheme1=intscheme;
  physcon1=physcon;shcon1=shcon;nshcon1=nshcon;cocon1=cocon;ncocon1=ncocon;
  ttime1=ttime;time1=time;istep1=istep;iinc1=iinc;coriolis1=coriolis;
  ibody1=ibody;xloadold1=xloadold;reltime1=reltime;veold1=veold;
  springarea1=springarea;nstate1_=nstate_;xstateini1=xstateini;xstate1=xstate;
  thicke1=thicke;integerglob1=integerglob;doubleglob1=doubleglob;
  tieset1=tieset;istartset1=istartset;iendset1=iendset;ialset1=ialset;
  ntie1=ntie;nasym1=nasym;pslavsurf1=pslavsurf;pmastsurf1=pmastsurf;
  mortar1=mortar;clearini1=clearini;ielprop1=ielprop;prop1=prop;ne01=ne0;
  kscale1=kscale;iponoel1=iponoel;inoel1=inoel;network1=network;
  set1=set;nset1=nset;
  
  /* calculating the stiffness/mass */
    
  printf(" Using up to %" ITGFORMAT " cpu(s) for the asymmetric stiffness/mass contributions.\n\n", num_cpus);
    
  /* create threads and wait */
    
  NNEW(ithread,ITG,num_cpus);
  for(i=0; i<num_cpus; i++)  {
    ithread[i]=i;
    pthread_create(&tid[i], NULL, (void *)mafillsmasmt, (void *)&ithread[i]);
  }
  for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
    
  SFREE(ithread);SFREE(neaparm);SFREE(nebparm);
  SFREE(neapart);SFREE(nebpart);
    
  /* copying and accumulating the stiffnes matrices */
    
  for(i=0;i<neq[1];i++){
    for(j=0;j<num_cpus;j++){
      ad[i]+=ad1[i+j*neq[1]];
    }
  }
  SFREE(ad1);

  for(i=0;i<(nzs[2]+nzs[1]);i++){
    for(j=0;j<num_cpus;j++){
      au[i]+=au1[i+(long long)j*(nzs[2]+nzs[1])];
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

/* subroutine for multithreading of mafillsm */

void *mafillsmasmt(ITG *i){

  ITG indexad,neam,nebm,neat,nebt;
  long long indexau;

  indexad=0;
  indexau=0;

  indexad=*i*neq1[1];
  indexau=(long long)*i*(nzs1[2]+nzs1[1]);

  neam=neaparm[*i]+1;
  nebm=nebparm[*i]+1;

  neat=neapart[*i]+1;
  nebt=nebpart[*i]+1;

  FORTRAN(mafillsmas,(co1,nk1,kon1,ipkon1,lakon1,ne1,nodeboun1,ndirboun1,
		      xboun1,nboun1,ipompc1,nodempc1,coefmpc1,nmpc1,
		      nodeforc1,ndirforc1,xforc1,nforc1,nelemload1,sideload1,
		      xload1,nload1,xbody1,ipobody1,nbody1,cgr1,
		      &ad1[indexad],&au1[indexau],bb1,nactdof1,icol1,jq1,
		      irow1,neq1,nzl1,&nmethod1[*i],ikmpc1,ilmpc1,ikboun1,
		      ilboun1,elcon1,nelcon1,rhcon1,nrhcon1,alcon1,nalcon1,
		      alzero1,ielmat1,ielorien1,norien1,orab1,ntmat1_,t01,
		      t11,ithermal1,prestr1,iprestr1,vold1,iperturb1,sti1,
		      nzs1,stx1,adb1,aub1,iexpl1,plicon1,nplicon1,plkcon1,
		      nplkcon1,xstiff1,npmat1_,dtime1,matname1,mi1,ncmat1_,
		      mass1,stiffness1,buckling1,rhsi1,intscheme1,physcon1,
		      shcon1,nshcon1,cocon1,ncocon1,ttime1,time1,istep1,
		      iinc1,coriolis1,ibody1,xloadold1,reltime1,veold1,
		      springarea1,nstate1_,xstateini1,xstate1,thicke1,
		      integerglob1,doubleglob1,tieset1,istartset1,iendset1,
		      ialset1,ntie1,nasym1,pslavsurf1,pmastsurf1,mortar1,
		      clearini1,ielprop1,prop1,ne01,kscale1,iponoel1,inoel1,
		      network1,&neam,&nebm,&neat,&nebt,set1,nset1));

  return NULL;
}
