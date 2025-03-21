/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2011 Guido Dhondt                          */

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
#include <time.h>
#include <string.h>
#include "CalculiX.h"
#include "mortar.h"

/* function called before solver transforming the SPCs/MPCs, the matrix 
   and the right hand side for quadratic elements 
   (see Phd-thesis Sitzmann,Chapter 4)

   Author: Saskia Sitzmann */

void premortar(ITG *nzs,ITG *nzsc2,
	       double **auc2p,double **adc2p,ITG **irowc2p,ITG **icolc2p,
	       ITG **jqc2p,
	       double **aubdp,ITG **irowbdp,ITG **jqbdp,
	       double **aubdtilp,ITG **irowbdtilp,ITG **jqbdtilp,
	       double **aubdtil2p,ITG **irowbdtil2p,ITG **jqbdtil2p,
	       double **auddp,ITG **irowddp,ITG **jqddp,
	       double **auddtilp,ITG **irowddtilp,ITG **jqddtilp,
	       double **auddtil2p,ITG **irowddtil2p,ITG **jqddtil2p,
	       double **auddinvp,ITG **irowddinvp,ITG **jqddinvp,
	       ITG **jqtempp,ITG **irowtempp,ITG **icoltempp,ITG *nzstemp,
	       ITG *iit,
	       ITG *icol,ITG *irow,ITG *jq,
	       ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
	       ITG *imastnode,ITG *nmastnode,
	       double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,
	       ITG *ne,double *stn,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG *ielmat,
	       ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,ITG *ithermal,double *prestr,
	       ITG *iprestr,char *filab,double *eme,double *emn,
	       double *een,ITG *iperturb,ITG *nactdof,
	       ITG *iout,double *qa,
	       double *vold,double *b,ITG *nodeboun,ITG *ndirboun,
	       double *xbounact,double *xboun,ITG *nboun,ITG *ipompc,
	       ITG *nodempc,
	       double *coefmpc,char *labmpc,ITG *nmpc,ITG *nmethod,
	       ITG *neq,double *veold,double *accold,
	       double *dtime,double *time,
	       double *ttime,double *plicon,
	       ITG *nplicon,double *plkcon,ITG *nplkcon,
	       double *xstateini,double *xstiff,double *xstate,ITG *npmat_,
	       char *matname,ITG *mi,ITG *ielas,
	       ITG *icmd,ITG *ncmat_,ITG *nstate_,double *stiini,
	       double *vini,double *ener,
	       double *enern,double *emeini,double *xstaten,double *eei,
	       double *enerini,double *cocon,ITG *ncocon,char *set,
	       ITG *nset,ITG *istartset,
	       ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
	       char *prset,double *qfx,double *qfn,double *trab,
	       ITG *inotr,ITG *ntrans,ITG *nelemload,
	       ITG *nload,ITG *istep,ITG *iinc,
	       double *springarea,double *reltime,ITG *ne0,double *xforc,
	       ITG *nforc,double *thicke,
	       double *shcon,ITG *nshcon,char *sideload,double *xload,
	       double *xloadold,ITG *icfd,ITG *inomat,
	       ITG *islavquadel,ITG *islavsurf,
	       ITG *iponoels,ITG *inoels,                      
	       ITG *mortar,ITG *nslavnode,ITG *islavnode,ITG *nslavs,
	       ITG *ntie,
	       double *aut,ITG *irowt,ITG *jqt,
	       double *autinv,ITG *irowtinv,ITG *jqtinv,
	       char *tieset,ITG *itiefac  ,ITG *rhsi,
	       double *au,double *ad,double **f_cmp,double **f_csp,
	       double *t1act,double *cam,double *bet,double *gam,
	       double *epn,
	       double *xloadact,ITG *nodeforc,ITG *ndirforc,double *xforcact,
	       double *xbodyact,ITG *ipobody,ITG *nbody,double *cgr,
	       ITG *nzl,double *sti,ITG *iexpl,ITG *mass,ITG *buckling,
	       ITG *stiffness,
	       ITG *intscheme,double *physcon,ITG *coriolis,ITG *ibody,
	       ITG *integerglob,double *doubleglob,ITG *nasym,
	       double *alpham,double *betam,
	       double *pslavsurf,double *pmastsurf,
	       double *clearini,ITG *ielprop,double *prop,
	       ITG *islavact,double *cdn,ITG *memmpc_,
	       ITG *idamping,
	       ITG *iforbou,ITG *iperturb_sav,
	       ITG *itietri,double *cg,double *straight,ITG *koncont,
	       double *energyini,
	       double *energy,ITG *kscale,ITG *iponoel,ITG *inoel,ITG *nener,
	       char *orname,ITG *network,
	       char *typeboun,ITG *num_cpus,double *t0g,double *t1g,
	       double *smscale,ITG *mscalmethod,ITG *nslavquadel){
  
  ITG im,i,k,mt=mi[1]+1,*irowc2=NULL,*icolc2=NULL,*jqc2=NULL,*irowbd=NULL,
    *jqbd=NULL,*irowbdtil=NULL,*jqbdtil=NULL,*irowbdtil2=NULL,*jqbdtil2=NULL,
    *irowdd=NULL,*jqdd=NULL,*irowddtil=NULL,*jqddtil=NULL,
    *irowddinv=NULL,*jqddinv=NULL,*irowddtil2=NULL,*jqddtil2=NULL,
    *irowtemp=NULL,*icoltemp=NULL,*jqtemp=NULL,*inum=NULL,
    mortartrafoflag=1;
  
  double alpha,*auc2=NULL,*adc2=NULL,*aubd=NULL,*aux2=NULL,*cv=NULL,
    *aubdtil=NULL,*aubdtil2=NULL,*f=NULL,*fext=NULL,*cvini=NULL,
    *audd=NULL,*auddtil=NULL,*auddinv=NULL,*auddtil2=NULL,*v=NULL,
    *stx=NULL,*fn=NULL,*fini=NULL,*fextini=NULL,*dam=NULL,*damn=NULL,
    *fmpc2=NULL,*aub=NULL,*adb=NULL,
    *adc=NULL,*auc=NULL,*volddummy=NULL,
    *vectornull=NULL,*f_cs=NULL,*f_cm=NULL,*fnext=NULL;
    
  alpha=1-2*sqrt(*bet);
  
  auc2=*auc2p;adc2=*adc2p;irowc2=*irowc2p;icolc2=*icolc2p;jqc2=*jqc2p;
  aubd=*aubdp;irowbd=*irowbdp;jqbd=*jqbdp;
  aubdtil=*aubdtilp;irowbdtil=*irowbdtilp;jqbdtil=*jqbdtilp;
  aubdtil2=*aubdtil2p;irowbdtil2=*irowbdtil2p;jqbdtil2=*jqbdtil2p;
  audd=*auddp;irowdd=*irowddp;jqdd=*jqddp;
  auddtil=*auddtilp;irowddtil=*irowddtilp;jqddtil=*jqddtilp;
  auddtil2=*auddtil2p;irowddtil2=*irowddtil2p;jqddtil2=*jqddtil2p;
  auddinv=*auddinvp;irowddinv=*irowddinvp;jqddinv=*jqddinvp;
  irowtemp=*irowtempp;icoltemp=*icoltempp;jqtemp=*jqtempp;
  f_cs=*f_csp;f_cm=*f_cmp;
  
  NNEW(f_cs,double,neq[1]);
  NNEW(f_cm,double,neq[1]);
  
  // check for coupled thermo-mechanical calculation
  
  if(ithermal[0]>1){
    printf("\tpremortar: coupled thermo-mechanical calculations NOT");
    printf("supported yet!\n \tPlease use surface-to-surface penalty contact");
    printf("instead.\n\n STOP!\n");
    fflush(stdout);
    FORTRAN(stop,());
  }
  
  // fix for linear calculation in first iteration of first increment

  if(*nslavquadel>0){
    if(*iforbou==1 && *iit==1 && *iinc==1){  
      *ielas=1;  
      iperturb[0]=-1;  
      iperturb[1]=0;	  
    }
  }
  
  /* small sliding is automatically set active due to combined fix-point
     Newton approach 
     do NOT change this unless the additional derivates neglected here 
     have been implemented */
  
  *nzsc2=nzs[1];
  NNEW(auc2,double,*nzsc2);
  NNEW(adc2,double,neq[1]);
  NNEW(irowc2,ITG,*nzsc2);
  NNEW(icolc2,ITG,neq[1]);
  NNEW(jqc2,ITG,neq[1]+1); 
  if(*iit==1){
    NNEW(aubd,double,6**nslavs);
    NNEW(irowbd,ITG,6**nslavs);
    NNEW(jqbd,ITG,neq[1]+1);
    NNEW(aubdtil,double,6**nslavs);
    NNEW(irowbdtil,ITG,6**nslavs);
    NNEW(jqbdtil,ITG,neq[1]+1);
    NNEW(aubdtil2,double,6**nslavs);
    NNEW(irowbdtil2,ITG,6**nslavs);
    NNEW(jqbdtil2,ITG,neq[1]+1);
    NNEW(audd,double,3**nslavs);
    NNEW(irowdd,ITG,3**nslavs);
    NNEW(jqdd,ITG,neq[1]+1);
    NNEW(auddtil,double,3**nslavs);
    NNEW(irowddtil,ITG,3**nslavs);
    NNEW(jqddtil,ITG,neq[1]+1);
    NNEW(auddtil2,double,3**nslavs);
    NNEW(irowddtil2,ITG,3**nslavs);
    NNEW(jqddtil2,ITG,neq[1]+1);
    NNEW(auddinv,double,3**nslavs);
    NNEW(irowddinv,ITG,3**nslavs);
    NNEW(jqddinv,ITG,neq[1]+1);	    
  }
  
  /* storing the original stiffness matrix */
  
  NNEW(jqtemp,ITG,neq[1]+1);
  NNEW(irowtemp,ITG,nzs[1]);
  NNEW(icoltemp,ITG,neq[1]);
  for(i=0;i<3;i++){
    nzstemp[i]=nzs[i];}
  for (i=0;i<neq[1];i++){jqtemp[i]=jq[i];icoltemp[i]=icol[i];}
  jqtemp[neq[1]]=jq[neq[1]];
  for (i=0;i<nzs[1];i++){irowtemp[i]=irow[i];}
    
  /* fix for quadratic FE */
    
  if(*nslavquadel>0){
    NNEW(v,double,mt**nk);
    NNEW(stx,double,6*mi[0]**ne);
    NNEW(fn,double,mt**nk);
    NNEW(fmpc2,double,*nmpc);
    memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
    *iout=-1;
    NNEW(f,double,neq[1]);
    NNEW(fext,double,neq[1]);
    NNEW(inum,ITG,*nk);
  
    /* calculating the internal forces and tangent stiffness using
       modified shape functions for quadratic elements */
  
    results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	    elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	    ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
	    prestr,iprestr,filab,eme,emn,een,iperturb,
	    f,fn,nactdof,iout,qa,vold,b,nodeboun,
	    ndirboun,xbounact,nboun,ipompc,
	    nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
	    bet,gam,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
	    xstateini,xstiff,xstate,npmat_,epn,matname,mi,ielas,icmd,
	    ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
	    xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
	    ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc2,
	    nelemload,nload,ikmpc,ilmpc,istep,iinc,springarea,
	    reltime,ne0,thicke,shcon,nshcon,
	    sideload,xloadact,xloadold,icfd,inomat,pslavsurf,pmastsurf,
	    mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
	    islavsurf,ielprop,prop,energyini,energy,kscale,iponoel,
	    inoel,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
	    itiefac,tieset,smscale,mscalmethod,nbody,t0g,t1g,
	    islavquadel,aut,irowt,jqt,&mortartrafoflag,
	    intscheme,physcon,dam,damn);
  
    SFREE(v);SFREE(stx);SFREE(fn);SFREE(inum);SFREE(fmpc2);
    *iout=0;	    
  
    *rhsi=1;
    DMEMSET(ad,0,neq[1],0.);
    DMEMSET(au,0,nzs[1],0.);
  
    /* calculating the external forces fext and stiffness matrix au/ad using
       modified shape functions for quadratic elements */
  
    mafillsmmain(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,nboun,
		 ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
		 nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,
		 nbody,cgr,ad,au,fext,nactdof,icol,jq,irow,
		 neq,nzl,
		 nmethod,ikmpc,ilmpc,ikboun,ilboun,
		 elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		 ielmat,ielorien,norien,orab,ntmat_,
		 t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
		 nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
		 xstiff,npmat_,dtime,matname,mi,
		 ncmat_,mass,stiffness,buckling,rhsi,intscheme,
		 physcon,shcon,nshcon,cocon,ncocon,ttime,time,istep,iinc,
		 coriolis,ibody,xloadold,reltime,veold,springarea,nstate_,
		 xstateini,xstate,thicke,integerglob,doubleglob,
		 tieset,istartset,iendset,ialset,ntie,nasym,pslavsurf,
		 pmastsurf,mortar,clearini,ielprop,prop,ne0,fnext,kscale,
		 iponoel,inoel,network,ntrans,inotr,trab,smscale,
		 mscalmethod,set,nset,islavquadel,aut,irowt,jqt,
		 &mortartrafoflag);

    /* calculating the residual b */
  
    calcresidual(nmethod,neq,b,fext,f,iexpl,nactdof,aux2,vold,
		 vini,dtime,accold,nk,adb,aub,jq,irow,nzl,
		 &alpha,fextini,fini,islavnode,nslavnode,mortar,ntie,
		 mi,nzs,nasym,idamping,veold,adc,auc,cvini,cv,
		 alpham,num_cpus);
  
    /*  for(k=0;k<neq[1];++k){printf("f=%" ITGFORMAT ",%f\n",k,f[k]);}
	for(k=0;k<neq[1];++k){printf("fext=%" ITGFORMAT ",%f\n",k,fext[k]);}
	for(k=0;k<neq[1];++k){printf("b=%" ITGFORMAT ",%f\n",k,b[k]);}
	for(k=0;k<neq[1];++k){printf("ad=%" ITGFORMAT ",%f\n",k,ad[k]);}
	for(k=0;k<nzs[1];++k){printf("au=%" ITGFORMAT ",%f\n",k,au[k]);}*/
  
    SFREE(f);SFREE(fext);
  }
  
  /* update vold due to spcs to get gap right for rigid body movements */
  
  /* if(*iinc==1 && *iit==1 && *nmethod!=4){       
     NNEW(v,double,mt**nk);	       
     NNEW(volddummy,double,mt**nk);
     for(k=0;k<mt**nk;k++){
     volddummy[k]=0.0;}           
     memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);	       
     NNEW(vectornull,double,neq[1]);	       
     *iout=-1;
    
     FORTRAN(resultsini_mortar,(nk,v,ithermal,iperturb,
     nactdof,iout,volddummy,vectornull,nodeboun,
     ndirboun,
     xbounact,nboun,ipompc,nodempc,coefmpc,labmpc,
     nmpc,nmethod,cam,
     bet,gam,dtime,mi));
    
     memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);	     	
     SFREE(v);SFREE(vectornull);SFREE(volddummy);   
     }*/
  
  *ielas=0;
  *iout=0;
  iperturb[0]=iperturb_sav[0];
  iperturb[1]=iperturb_sav[1];
  
  *auc2p=auc2;*adc2p=adc2;*irowc2p=irowc2;*icolc2p=icolc2;*jqc2p=jqc2;
  *aubdp=aubd;*irowbdp=irowbd;*jqbdp=jqbd;
  *aubdtilp=aubdtil;*irowbdtilp=irowbdtil;*jqbdtilp=jqbdtil;
  *aubdtil2p=aubdtil2;*irowbdtil2p=irowbdtil2;*jqbdtil2p=jqbdtil2;
  *auddp=audd;*irowddp=irowdd;*jqddp=jqdd;
  *auddtilp=auddtil;*irowddtilp=irowddtil;*jqddtilp=jqddtil;
  *auddtil2p=auddtil2;*irowddtil2p=irowddtil2;*jqddtil2p=jqddtil2;
  *auddinvp=auddinv;*irowddinvp=irowddinv;*jqddinvp=jqddinv;
  *irowtempp=irowtemp;*icoltempp=icoltemp;*jqtempp=jqtemp;
  *f_csp=f_cs;*f_cmp=f_cm;
  
  return;
}
