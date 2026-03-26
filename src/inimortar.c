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
#include "CalculiX.h"
#include "mortar.h"

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

/*   - allocating fields
     - determining the transformation matrix for the shape functions
       for quadratic elements
     - determining the quadratic elements which contain slave nodes
     - create fields islavspc, islavmpc, imastspc... connecting slave
       and/or master nodes with the spc/mpc's applied in them
     - remove Lagrange multipliers in nodes common to two slave
       surfaces and in slave nodes with MPC's
 
     Author: Saskia Sitzmann
*/

void inimortar(double **enerp,ITG *mi,ITG *ne ,ITG *nslavs,ITG *nk,ITG *nener,
	       ITG **ipkonp,char **lakonp,ITG **konp,ITG *nkon,
	       ITG *maxprevcontel,double **xstatep,ITG *nstate_,
	       ITG **islavtiep,double **bpp,ITG **islavactp,
	       double **gapp,
	       double **cdispp,double **cstressp,double **cfsp,
	       double **bpinip,ITG **islavactinip,double **cstressinip,
	       ITG *ntie,char *tieset,ITG *nslavnode,ITG *islavnode,
	       ITG **islavnodeinvp,ITG **islavquadelp,double **pslavdualp,
	       double **autp,ITG **irowtp,ITG **jqtp,	
	       double **autinvp,ITG **irowtinvp,ITG **jqtinvp,
	       double **Bdp,ITG **irowbp,ITG **jqbp,
	       double **Bdhelpp,ITG **irowbhelpp,ITG **jqbhelpp,
	       double **Ddp,ITG **irowdp,ITG **jqdp,
	       double **Ddtilp,ITG **irowdtilp,ITG **jqdtilp,
	       double **Bdtilp,ITG **irowbtilp,ITG **jqbtilp,
	       ITG *itiefac,ITG *islavsurf,ITG *nboun,ITG *nmpc,
	       ITG **nslavspcp,ITG **islavspcp,ITG **nslavmpcp,ITG **islavmpcp,
	       ITG **nmastspcp,ITG **imastspcp,ITG **nmastmpcp,ITG **imastmpcp,
	       ITG *imastnode,ITG *nmastnode,ITG *nasym,ITG *mortar,
	       ITG **ielmatp,ITG **ielorienp,ITG *norien,ITG *ipompc,
	       ITG *nodempc,ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
	       char *jobnamef,char *set,double *co,double *vold,ITG *nset,
	       ITG *nslavquadel){

  char *lakon=NULL;
    
  ITG k,i,j,node,mt=mi[1]+1,*ipkon=NULL,*kon=NULL,*islavtie=NULL,
    *islavact=NULL,*islavactini=NULL,*islavnodeinv=NULL,*islavquadel=NULL,
    *irowt=NULL,*jqt=NULL,nmspc,nmmpc,nsspc,nsmpc,*irowtinv=NULL,*jqtinv=NULL,
    *irowbhelp=NULL,*jqbhelp=NULL,*irowb=NULL,*jqb=NULL,*irowd=NULL,*jqd=NULL,
    *irowdtil=NULL,*jqdtil=NULL,*irowbtil=NULL,*jqbtil=NULL,
    *nslavspc=NULL,*islavspc=NULL,*nslavmpc=NULL,*islavmpc=NULL,
    *nmastspc=NULL,*imastspc=NULL,*nmastmpc=NULL,*imastmpc=NULL,
    *ielmat=NULL,*ielorien=NULL;
  
  double *ener=NULL,*xstate=NULL,*bp=NULL,*gap=NULL,
    *cdisp=NULL,*cstress=NULL,*cfs=NULL,
    *bpini=NULL,*cstressini=NULL,*pslavdual=NULL,*aut=NULL,
    *autinv=NULL,*Bd=NULL,*Bdhelp=NULL,*Dd=NULL,*Ddtil=NULL,*Bdtil=NULL;
  
  ener=*enerp;ipkon=*ipkonp;lakon=*lakonp;kon=*konp;xstate=*xstatep;
  islavtie=*islavtiep;bp=*bpp;islavact=*islavactp;gap=*gapp;
  cdisp=*cdispp;cstress=*cstressp;
  cfs=*cfsp;
  bpini=*bpinip;islavactini=*islavactinip;cstressini=*cstressinip;
  islavnodeinv=*islavnodeinvp;islavquadel=*islavquadelp;pslavdual=*pslavdualp;
  aut=*autp;irowt=*irowtp;jqt=*jqtp;	
  autinv=*autinvp;irowtinv=*irowtinvp;jqtinv=*jqtinvp;
  Bd=*Bdp;irowb=*irowbp;jqb=*jqbp;
  Bdhelp=*Bdhelpp;irowbhelp=*irowbhelpp;jqbhelp=*jqbhelpp;
  Dd=*Ddp;irowd=*irowdp;jqd=*jqdp;
  Ddtil=*Ddtilp;irowdtil=*irowdtilp;jqdtil=*jqdtilp;
  Bdtil=*Bdtilp;irowbtil=*irowbtilp;jqbtil=*jqbtilp;
  nslavspc=*nslavspcp;islavspc=*islavspcp;nslavmpc=*nslavmpcp;
  islavmpc=*islavmpcp;
  nmastspc=*nmastspcp;imastspc=*imastspcp;nmastmpc=*nmastmpcp;
  imastmpc=*imastmpcp;
  ielmat=*ielmatp;ielorien=*ielorienp;
  
  RENEW(ielmat,ITG,mi[2]*(*ne+*nslavs));
  for(k=mi[2]**ne;k<mi[2]*(*ne+*nslavs);k++) ielmat[k]=1;
  if(*norien>0){      
    RENEW(ielorien,ITG,mi[2]*(*ne+*nslavs));
    for(k=mi[2]**ne;k<mi[2]*(*ne+*nslavs);k++) ielorien[k]=0;
  }  
  
  if(*nener==1){
    RENEW(ener,double,mi[0]*(*ne+*nslavs)*2);
  }
  RENEW(ipkon,ITG,*ne+*nslavs);
  RENEW(lakon,char,8*(*ne+*nslavs));
  
  /* adding one element per slave node, similar to 
     spring elements;
     needed for output in frd-format of CDISP and CSTRES */
  
  RENEW(kon,ITG,*nkon+*nslavs);
  if((*maxprevcontel==0)&&(*nslavs!=0)){
    RENEW(xstate,double,*nstate_*mi[0]*(*ne+*nslavs));
    for(k=*nstate_*mi[0]**ne;k<*nstate_*mi[0]*(*ne+*nslavs);k++){
      xstate[k]=0.;
    }
  }      
  
  for(k=0;k<*nslavs;k++){
    ipkon[*ne+k]=*nkon+k;
    kon[*nkon+k]=islavnode[k];
    strcpy1(&lakon[8*(*ne+k)]," S    C0",8);
  }
  NNEW(islavtie,ITG,*nslavs);
  NNEW(bp,double,*nslavs);
  NNEW(islavact,ITG,*nslavs);
  NNEW(gap,double,*nslavs);
  NNEW(cdisp,double,6**nslavs);
  
  /* allocation of temporary fields: stores the structure
     of the stiffness matrix without mortar contact */
  
  NNEW(cstress,double,mt**nslavs);
  NNEW(cfs,double,mt**nk);

  /* fields for cutback  */
  
  NNEW(bpini,double,*nslavs);
  NNEW(islavactini,ITG,*nslavs); 
  NNEW(cstressini,double,mt**nslavs);
  
  /* connect each slave node j with its tie i */
  
  for(i=0;i<*ntie;i++){
    if(tieset[i*(81*3)+80]=='C'){
      for(j=nslavnode[i];j<nslavnode[i+1];j++){
	islavtie[j]=i;
      }
    }
  }
  
  /* get results from last step */
  
  if(*maxprevcontel!=0){
    for(i=0;i<*ntie;i++){
      if(tieset[i*(81*3)+80]=='C'){
	if(*nstate_*mi[0]>0){
	  for(j=nslavnode[i];j<nslavnode[i+1];j++){
	    for(k=0;k<3;k++){
	      cstress[mt*j+k]=xstate[*nstate_*mi[0]*(*ne+j)+k];
	      xstate[*nstate_*mi[0]*(*ne+j)+k]=0.;		
	    }
	    if(xstate[*nstate_*mi[0]*(*ne+j)+3]>1.0 &&
	       xstate[*nstate_*mi[0]*(*ne+j)+3]<2.0){
	      islavact[j]=1;
	    }else if(xstate[*nstate_*mi[0]*(*ne+j)+3]>2.0 &&
		     xstate[*nstate_*mi[0]*(*ne+j)+3]<3.0){
	      islavact[j]=2;
	    }
	    xstate[*nstate_*mi[0]*(*ne+j)+3]=0.;
	  }
	}
      }
    }	
  }
  NNEW(islavnodeinv,ITG,*nk);
  NNEW(islavquadel,ITG,*ne);
  for( i=0;i<*ntie;i++){    
    if(tieset[i*(81*3)+80]=='C'){ 
      for(j=nmastnode[i];j<nmastnode[i+1];j++){
	node=imastnode[j];
	islavnodeinv[node-1]=-(j+1);
      }				
    }
  }    
  for( i=0;i<*ntie;i++){    
    if(tieset[i*(81*3)+80]=='C'){ 
      for(j=nslavnode[i];j<nslavnode[i+1];j++){
	node=islavnode[j];
	islavnodeinv[node-1]=j+1;
      }			
    }
  }
  
  /* coeffs for dual basis functions */
  
  NNEW(pslavdual,double,64*itiefac[2**ntie-1]);
  
  /*  T and T^-1 and coupling matrices in nodes  */
  
  NNEW(aut,double,3**nslavs);
  NNEW(irowt,ITG,3**nslavs);
  NNEW(jqt,ITG,*nk+1);	
  NNEW(autinv,double,3**nslavs);
  NNEW(irowtinv,ITG,3**nslavs);
  NNEW(jqtinv,ITG,*nk+1);
  
  buildtquad(ntie,ipkon,kon,nk,lakon,nslavnode,itiefac,tieset,
	     islavnode,islavsurf,&irowt,jqt,&aut,
	     &irowtinv,jqtinv,&autinv);
  
  FORTRAN(genislavquadel,(islavquadel,jqt,lakon,ipkon,kon,ne,nasym,
			  nslavquadel));
  
  /* checking for SPC's and MPC's on slave and master surface */

  NNEW(nslavspc,ITG,2**nslavs);
  NNEW(islavspc,ITG,*nboun);
  NNEW(nslavmpc,ITG,2**nslavs);
  NNEW(islavmpc,ITG,*nmpc);
  NNEW(nmastspc,ITG,2*nmastnode[*ntie]);
  NNEW(imastspc,ITG,*nboun);
  NNEW(nmastmpc,ITG,2*nmastnode[*ntie]);
  NNEW(imastmpc,ITG,*nmpc);
  
  FORTRAN(spcmpcmortar,(ntie,islavnode,imastnode,nslavnode,
			nmastnode,nboun,nmpc,
			ipompc,nodempc,ikboun,ilboun,ikmpc,ilmpc,
			nslavspc,islavspc,
			&nsspc,nslavmpc,islavmpc,&nsmpc,
			nmastspc,imastspc,&nmspc,nmastmpc,
			imastmpc,&nmmpc,jobnamef));
  
  RENEW(islavspc,ITG,nsspc);
  RENEW(islavmpc,ITG,nsmpc);
  RENEW(imastspc,ITG,nmspc);
  RENEW(imastmpc,ITG,nmmpc);
    
    /* remove the Lagrange Multiplier from nodes common to several
       slave surfaces and from node which belong to a MPC */
    
  FORTRAN(remlagrangemult,(ntie,tieset,islavnode,imastnode,nslavnode,nmastnode,
			   islavact,nodempc,nmpc,ipompc));
  
  NNEW(Bd,double,1);
  NNEW(irowb,ITG,1);
  NNEW(jqb,ITG,*nk+1);
  NNEW(Bdhelp,double,1);
  NNEW(irowbhelp,ITG,1);
  NNEW(jqbhelp,ITG,*nk+1);
  NNEW(Dd,double,1);
  NNEW(irowd,ITG,1);
  NNEW(jqd,ITG,*nk+1);
  NNEW(Ddtil,double,1);
  NNEW(irowdtil,ITG,1);
  NNEW(jqdtil,ITG,*nk+1);
  NNEW(Bdtil,double,1);
  NNEW(irowbtil,ITG,1);
  NNEW(jqbtil,ITG,*nk+1);
  
  *enerp=ener;*ipkonp=ipkon;*lakonp=lakon;*konp=kon;*xstatep=xstate;
  *islavtiep=islavtie;*bpp=bp;*islavactp=islavact;*gapp=gap;
  *cdispp=cdisp;*cstressp=cstress;
  *cfsp=cfs;
  *bpinip=bpini;*islavactinip=islavactini;*cstressinip=cstressini;
  *islavnodeinvp=islavnodeinv;*islavquadelp=islavquadel;*pslavdualp=pslavdual;
  *autp=aut;*irowtp=irowt;*jqtp=jqt;	
  *autinvp=autinv;*irowtinvp=irowtinv;*jqtinvp=jqtinv;
  *Bdp=Bd;*irowbp=irowb;*jqbp=jqb;
  *Bdhelpp=Bdhelp;*irowbhelpp=irowbhelp;*jqbhelpp=jqbhelp;
  *Ddp=Dd;*irowdp=irowd;*jqdp=jqd;
  *Ddtilp=Ddtil;*irowdtilp=irowdtil;*jqdtilp=jqdtil;
  *Bdtilp=Bdtil;*irowbtilp=irowbtil;*jqbtilp=jqbtil;
  *nslavspcp=nslavspc;*islavspcp=islavspc;*nslavmpcp=nslavmpc;
  *islavmpcp=islavmpc;
  *nmastspcp=nmastspc;*imastspcp=imastspc;*nmastmpcp=nmastmpc;
  *imastmpcp=imastmpc;
  *ielmatp=ielmat;*ielorienp=ielorien;
  
  return;
}
