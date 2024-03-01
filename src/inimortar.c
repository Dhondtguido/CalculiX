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
/**   function initializing mortar contact at the start of nonlingeo.c
 *   getting contact parameters, allocating needed fields, get results from last step,
 *   determine used mortar method, (transform and cataloque SPCs/MPCs)
 *   and generate islavelinv and  islavnodeinv
 *
 *   Author: Saskia Sitzmann
 *
 *  [in,out] enerp		?
 *  [out] islavactdoftiep   (i)=tie number for active dof i
 *  [out] bpp		friction bounds 
 *  [out] islavactp	(i) indicates, if slave node i is active (=-3 no-slave-node, =-2 no-LM-node, =-1 no-gap-node, =0 inactive node, =1 sticky node, =2 slipping/active node)  
 *  [out] gapp		(i) gap for node i on slave surface
 *  [out] slavnorp		slave normal
 *  [out] slavtanp		slave tangent 
 *  [out] cdispp		vector saving contact variables for frd-output 
 *  [out] cstressp	current Lagrange multiplier 
 *  [out] cfsp 		contact force 
 *  [out] cfmp 		not used any more 
 *  [out] bpinip		friction bounds at start of the increment
 *  [out] islavactinip	islavact at the start of the increment
 *  [out] cstressinip	Lagrange multiplier at start of the increment
 *  [out] islavnodeinvp     (i) slave node index for node i
 *  [out] islavelinvp       (i)==0 if there is no slave node in the element, >0 otherwise
 *  [out] pslavdualp	(:,i)coefficients \f$ \alpha_{ij}\f$, \f$ 1,j=1,..8\f$ for dual shape functions for face i
 *  [out] pslavdualpgp	(:,i)coefficients \f$ \alpha_{ij}\f$, \f$ 1,j=1,..8\f$ for Petrov-Galerkin shape functions for face i 
 *  [out] autlocp		transformation matrix \f$ T[p,q]\f$ for slave nodes \f$ p,q \f$
 *  [out] irowtlocp		field containing row numbers of autloc
 *  [out] jqtlocp	        pointer into field irowtloc
 *  [out] autlocinvp	transformation matrix \f$ T^{-1}[p,q]\f$ for slave nodes \f$ p,q \f$ 
 *  [out] irowtlocinvp	field containing row numbers of autlocinv
 *  [out] jqtlocinvp	pointer into field irowtlocinv
 *  [out] Bdp		coupling matrix \f$ B_d[p,q]=\int \psi_p \phi_q dS \f$, \f$ p \in S, q \in M \f$ 
 *  [out] irowbp		field containing row numbers of Bd
 *  [out] jqbp		pointer into field irowb
 *  [out] Bdhelpp		coupling matrix \f$ Bhelp_d[p,q]=\tilde{D}^{-1}\tilde{B}\f$, \f$ p \in S, q \in M \f$ 
 *  [out] irowbhelpp	field containing row numbers of Bdhelp
 *  [out] jqbhelpp		pointer into field irowbhelp
 *  [out] Ddp		coupling matrix \f$ D_d[p,q]=\int \psi_p \phi_q dS \f$, \f$ p,q \in S \f$ 
 *  [out] irowdp		field containing row numbers of Dd
 *  [out] jqdp		pointer into field irowd
 *  [out] Ddtilp		coupling matrix \f$ \tilde{D}_d[p,q]=\int \psi_p \tilde{\phi}_q dS \f$, \f$ p,q \in S \f$ 
 *  [out] irowdtilp	field containing row numbers of Ddtil
 *  [out] jqdtilp		pointer into field irowdtil 
 *  [out] Bdtilp		coupling matrix \f$ \tilde{B}_d[p,q]=\int \psi_p \tilde{\phi}_q dS \f$, \f$ p \in S, q \in M \f$ 
 *  [out] irowbtilp	field containing row numbers of Bdtil
 *  [out] jqbtilp		pointer into field irowbtil
 *  [out] Bpgdp		Petrov-Galerkin coupling matrix \f$ B_d^{PG}[p,q]=\int \tilde{\phi}_p \phi_q dS \f$, \f$ p \in S, q \in M \f$ 
 *  [out] irowbpgp		field containing row numbers of Bpgd
 *  [out] jqbpgp		pointer into field irowbpg
 *  [out] Dpgdp		Petrov-Galerkin coupling matrix \f$ D_d[p,q]=\int \tilde{\phi}_p \phi_q dS \f$, \f$ p,q \in S \f$ 
 *  [out] irowdpgp		field containing row numbers of Dpgd
 *  [out] jqdpgp		pointer into field irowdpg
 *  [out] Dpgdtilp		transformed Petrov-Galerkin coupling matrix \f$ D_d[p,q]=\int \tilde{\phi}_p \tilde{\phi}_q dS \f$, \f$ p,q \in S \f$ 
 *  [out] irowdpgtilp	field containing row numbers of Dpgdtil
 *  [out] jqdpgtilp	pointer into field irowdpgtil
 *  [out] Bpgdtilp		transformed Petrov-Galerkin coupling matrix \f$ B_d^{PG}[p,q]=\int \tilde{\phi}_p \tilde{\phi}_q dS \f$, \f$ p \in S, q \in M \f$ 
 *  [out] irowbpgtilp	field containing row numbers of Bpgdtil
 *  [out] jqbpgtilp	pointer into field irowbpgtil
 *  [in]  iflagdualquad   flag indicating what mortar contact is used (=1 quad-lin, =2 quad-quad, =3 PG quad-lin, =4 PG quad-quad)
 *  [in] itiefac 		pointer into field islavsurf: (1,i) beginning slave_i (2,i) end of slave_i
 *  [out] nslavspcp	(2*i) pointer to islavspc...
 *  [out] islavspcp         ... which stores SPCs for slave node i
 *  [out] nslavmpcp	(2*i) pointer to islavmpc...
 *  [out] islavmpcp	... which stores MPCs for slave node i
 *  [out] nmastspcp	(2*i) pointer to imastspc...
 *  [out] imastspcp         ... which stores SPCs for master node i
 *  [out] nmastmpcp	(2*i) pointer to imastmpc...
 *  [out] imastmpcp	... which stores MPCs for master node i 
 *  [out] nsspc            number of SPC for slave nodes 
 *  [out] nsmpc		number of MPC for slave nodes
 *  [out] nmspc            number of SPC for master nodes
 *  [out] nmmpc		number of MPC for master nodes 
 *  [in] iponoels         (i) pointer to inoels
 *  [in] inoels           (3,i) element number, local node number and pointer to another entry
**/
void inimortar(double **enerp,ITG *mi,ITG *ne ,ITG *nslavs,ITG *nk,ITG *nener,
	       ITG **ipkonp,char **lakonp,ITG **konp,ITG *nkon,
	       ITG *maxprevcontel,double **xstatep,ITG *nstate_,
	       ITG **islavactdoftiep,double **bpp,ITG **islavactp,
	       double **gapp,double **slavnorp,double **slavtanp,
	       double **cdispp,
	       double **cstressp,double **cfsp,double **cfmp,
	       double **cfsinip,double **cfsinitilp,double **cfstilp,
	       double **bpinip,ITG **islavactinip,double **cstressinip,
	       ITG *ntie,char *tieset,
	       ITG *nslavnode,ITG *islavnode,
	       ITG **islavnodeinvp,ITG **islavelinvp,double **pslavdualp,
	       double **pslavdualpgp,
	       double **autlocp,ITG **irowtlocp,ITG **jqtlocp,	
	       double **autlocinvp,ITG **irowtlocinvp,ITG **jqtlocinvp,
	       double **Bdp,ITG **irowbp,ITG **jqbp,
	       double **Bdhelpp,ITG **irowbhelpp,ITG **jqbhelpp,
	       double **Ddp,ITG **irowdp,ITG **jqdp,
	       double **Ddtilp,ITG **irowdtilp,ITG **jqdtilp,
	       double **Bdtilp,ITG **irowbtilp,ITG **jqbtilp,
	       double **Bpgdp,ITG **irowbpgp,ITG **jqbpgp,
	       double **Dpgdp,ITG **irowdpgp,ITG **jqdpgp,
	       double **Dpgdtilp,ITG **irowdpgtilp,ITG **jqdpgtilp,
	       double **Bpgdtilp,ITG **irowbpgtilp,ITG **jqbpgtilp,
	       ITG *iflagdualquad,ITG *itiefac,ITG *islavsurf,
	       ITG *nboun,ITG *ndirboun,ITG *nodeboun,double *xboun,
	       ITG *nmpc,ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
	       ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
	       ITG **nslavspcp,ITG **islavspcp,ITG **nslavmpcp,ITG **islavmpcp,
	       ITG **nmastspcp,ITG **imastspcp,ITG **nmastmpcp,ITG **imastmpcp,
	       ITG *nsspc,ITG *nsmpc,
	       ITG *imastnode,ITG *nmastnode,ITG *nmspc,ITG *nmmpc,
	       ITG *iponoels,ITG *inoels,
	       double *tietol,double *elcon,ITG *ncmat_,ITG *ntmat_,ITG *nasym,
	       double *vold,ITG *nset,char *set,ITG *mortar,
	       ITG *memmpc_,
	       ITG **ielmatp,ITG **ielorienp,ITG *norien,ITG *nmethod,
	       ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc){

  char *lakon=NULL;
    
  ITG k,i,j,node,mt=mi[1]+1,*ipkon=NULL,*kon=NULL,
    *islavactdoftie=NULL,*islavact=NULL,
    *islavactini=NULL,*islavnodeinv=NULL,*islavelinv=NULL,*irowtloc=NULL,
    *jqtloc=NULL,
    *irowtlocinv=NULL,*jqtlocinv=NULL,*irowbhelp=NULL,*jqbhelp=NULL,
    *irowb=NULL,*jqb=NULL,
    *irowd=NULL,*jqd=NULL,*irowdtil=NULL,*jqdtil=NULL,
    *irowbtil=NULL,*jqbtil=NULL,*irowbpg=NULL,*jqbpg=NULL,
    *irowdpg=NULL,*jqdpg=NULL,*irowdpgtil=NULL,*jqdpgtil=NULL,
    *irowbpgtil=NULL,*jqbpgtil=NULL,
    *nslavspc=NULL,*islavspc=NULL,*nslavmpc=NULL,*islavmpc=NULL,
    *nmastspc=NULL,*imastspc=NULL,*nmastmpc=NULL,*imastmpc=NULL,
    *ielmat=NULL,*ielorien=NULL;
  
  double *ener=NULL,*xstate=NULL,*bp=NULL,
    *gap=NULL,*slavnor=NULL,*slavtan=NULL,
    *cdisp=NULL,*cstress=NULL,*cfs=NULL,*cfsini=NULL,*cfsinitil=NULL,
    *cfstil=NULL,*cfm=NULL,*bpini=NULL,*cstressini=NULL,*pslavdual=NULL,
    *pslavdualpg=NULL,*autloc=NULL,
    *autlocinv=NULL,*Bd=NULL,*Bdhelp=NULL,*Dd=NULL,*Ddtil=NULL,*Bdtil=NULL,
    *Bpgd=NULL,*Dpgd=NULL,*Dpgdtil=NULL,*Bpgdtil=NULL;
  
  ener=*enerp;ipkon=*ipkonp;lakon=*lakonp;kon=*konp;xstate=*xstatep;
  islavactdoftie=*islavactdoftiep;bp=*bpp;islavact=*islavactp;gap=*gapp;
  slavnor=*slavnorp;slavtan=*slavtanp;cdisp=*cdispp;cstress=*cstressp;
  cfs=*cfsp;cfsini=*cfsinip;cfsinitil=*cfsinitilp;cfstil=*cfstilp;
  cfm=*cfmp;bpini=*bpinip;islavactini=*islavactinip;cstressini=*cstressinip;
  islavnodeinv=*islavnodeinvp;islavelinv=*islavelinvp;pslavdual=*pslavdualp;
  pslavdualpg=*pslavdualpgp;
  autloc=*autlocp;irowtloc=*irowtlocp;jqtloc=*jqtlocp;	
  autlocinv=*autlocinvp;irowtlocinv=*irowtlocinvp;jqtlocinv=*jqtlocinvp;
  Bd=*Bdp;irowb=*irowbp;jqb=*jqbp;
  Bdhelp=*Bdhelpp;irowbhelp=*irowbhelpp;jqbhelp=*jqbhelpp;
  Dd=*Ddp;irowd=*irowdp;jqd=*jqdp;
  Ddtil=*Ddtilp;irowdtil=*irowdtilp;jqdtil=*jqdtilp;
  Bdtil=*Bdtilp;irowbtil=*irowbtilp;jqbtil=*jqbtilp;
  Bpgd=*Bpgdp;irowbpg=*irowbpgp;jqbpg=*jqbpgp;
  Dpgd=*Dpgdp;irowdpg=*irowdpgp;jqdpg=*jqdpgp;
  Dpgdtil=*Dpgdtilp;irowdpgtil=*irowdpgtilp;jqdpgtil=*jqdpgtilp;
  Bpgdtil=*Bpgdtilp;irowbpgtil=*irowbpgtilp;jqbpgtil=*jqbpgtilp; 
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
  NNEW(islavactdoftie,ITG,*nslavs);
  NNEW(bp,double,*nslavs);
  NNEW(islavact,ITG,*nslavs);
  NNEW(gap,double,*nslavs);
  NNEW(slavnor,double,3**nslavs);
  NNEW(slavtan,double,6**nslavs);
  NNEW(cdisp,double,6**nslavs);
  
  /* allocation of temperary fields: stores the structure
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
	islavactdoftie[j]=i;
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
  NNEW(islavelinv,ITG,*ne);
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
  
  NNEW(autloc,double,3**nslavs);
  NNEW(irowtloc,ITG,3**nslavs);
  NNEW(jqtloc,ITG,*nk+1);	
  NNEW(autlocinv,double,3**nslavs);
  NNEW(irowtlocinv,ITG,3**nslavs);
  NNEW(jqtlocinv,ITG,*nk+1);
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
  NNEW(Bpgd,double,1);
  NNEW(irowbpg,ITG,1);
  NNEW(jqbpg,ITG,*nk+1);
  NNEW(Dpgd,double,1);
  NNEW(irowdpg,ITG,1);
  NNEW(jqdpg,ITG,*nk+1);
  NNEW(Dpgdtil,double,1);
  NNEW(irowdpgtil,ITG,1);
  NNEW(jqdpgtil,ITG,*nk+1);
  NNEW(Bpgdtil,double,1);
  NNEW(irowbpgtil,ITG,1);
  NNEW(jqbpgtil,ITG,*nk+1);
  
  /* iflagdualquad==1 : linear-quadratic (mortar==2: linmortar)
   * iflagdualquad==2 : quadratric-quadratic (mortar==3: mortar)
   * iflagdualquad==3 : petrov-galerkin lin-quad (mortar==4: pglinmortar)
   * iflagdualquad==4 : petrov-galerkin quad-quad (mortar==5: pgmortar) */
  
  if(*mortar==2){
    *iflagdualquad=1;
  }else if(*mortar==3){
    *iflagdualquad=2;
  }else if(*mortar==4){
    *iflagdualquad=3;
  }else if(*mortar==5){
    *iflagdualquad=4;
  }else{
    *iflagdualquad=2;
  }
  
  if((*iflagdualquad==3)||(*iflagdualquad==4)){
    NNEW(pslavdualpg,double,64*itiefac[2**ntie-1]);
  }else{
    NNEW(pslavdualpg,double,64);
  }
  
  buildtquad(ntie,ipkon,kon,nk,lakon,nslavnode,itiefac,tieset,
	     islavnode,islavsurf,&irowtloc,jqtloc,&autloc,
	     &irowtlocinv,jqtlocinv,&autlocinv,iflagdualquad);
  
  /* checking for SPC's and MPC's on slave and master surface */

  NNEW(nslavspc,ITG,2**nslavs);
  NNEW(islavspc,ITG,2**nboun);
  NNEW(nslavmpc,ITG,2**nslavs);
  NNEW(islavmpc,ITG,2**nmpc);
  NNEW(nmastspc,ITG,2*nmastnode[*ntie]);
  NNEW(imastspc,ITG,2**nboun);
  NNEW(nmastmpc,ITG,2*nmastnode[*ntie]);
  NNEW(imastmpc,ITG,2**nmpc);
  
  FORTRAN(genislavelinv,(islavelinv,jqtloc,lakon,ipkon,kon,ne,nasym));
  
  *enerp=ener;*ipkonp=ipkon;*lakonp=lakon;*konp=kon;*xstatep=xstate;
  *islavactdoftiep=islavactdoftie;*bpp=bp;*islavactp=islavact;*gapp=gap;
  *slavnorp=slavnor;*slavtanp=slavtan;*cdispp=cdisp;*cstressp=cstress;
  *cfsp=cfs;*cfsinip=cfsini;*cfsinitilp=cfsinitil;*cfstilp=cfstil;
  *cfmp=cfm;*bpinip=bpini;*islavactinip=islavactini;*cstressinip=cstressini;
  *islavnodeinvp=islavnodeinv;*islavelinvp=islavelinv;*pslavdualp=pslavdual;
  *pslavdualpgp=pslavdualpg;
  *autlocp=autloc;*irowtlocp=irowtloc;*jqtlocp=jqtloc;	
  *autlocinvp=autlocinv;*irowtlocinvp=irowtlocinv;*jqtlocinvp=jqtlocinv;
  *Bdp=Bd;*irowbp=irowb;*jqbp=jqb;
  *Bdhelpp=Bdhelp;*irowbhelpp=irowbhelp;*jqbhelpp=jqbhelp;
  *Ddp=Dd;*irowdp=irowd;*jqdp=jqd;
  *Ddtilp=Ddtil;*irowdtilp=irowdtil;*jqdtilp=jqdtil;
  *Bdtilp=Bdtil;*irowbtilp=irowbtil;*jqbtilp=jqbtil;
  *Bpgdp=Bpgd;*irowbpgp=irowbpg;*jqbpgp=jqbpg;
  *Dpgdp=Dpgd;*irowdpgp=irowdpg;*jqdpgp=jqdpg;
  *Dpgdtilp=Dpgdtil;*irowdpgtilp=irowdpgtil;*jqdpgtilp=jqdpgtil;
  *Bpgdtilp=Bpgdtil;*irowbpgtilp=irowbpgtil;*jqbpgtilp=jqbpgtil;
  *nslavspcp=nslavspc;*islavspcp=islavspc;*nslavmpcp=nslavmpc;
  *islavmpcp=islavmpc;
  *nmastspcp=nmastspc;*imastspcp=imastspc;*nmastmpcp=nmastmpc;
  *imastmpcp=imastmpc;
  *ielmatp=ielmat;*ielorienp=ielorien;
  
  return;
}
