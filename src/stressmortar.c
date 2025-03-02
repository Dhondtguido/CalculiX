/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2011 Guido Dhondt                     */

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
#define abs(a) (((a) < (0)) ? (-a) : (a))

/**  transforming the system back to the standard basis functions 
     \f$ \tilde{u}\rightarrow u \f$ and updating the active
     *       and inactive sets and the Langrange Multipliers (LM) 
     *       see Sitzmann Algorithm 2, p.71
     * 
     * Author: Saskia Sitzmann
     *
     *  [in] 	bhat		intermediate right hand side 
     *  [in] 	adc		intermediate system matrix, diagonal terms
     *  [in] 	auc		intermediate system matrix, nondiagonal terms
     *  [in] 	jqc             pointer to irowc
     *  [in] 	irowc		row numbers of auc
    *  [in] 	gap		(i) gap for slave node i
     *  [in,out] b		in: differenzial displacement out:real displacement
     *  [in,out] islavact	(i) indicates, if slave node i is active (=-3 no-slave-node, =-2 no-LM-node, =-1 no-gap-node, =0 inactive node, =1 sticky node, =2 slipping/active node) 
     *  [out] irowddinv	field containing row numbers of auddinv
     *  [out] jqddinv		pointer into field irowddinv
     *  [out] auddinv		coupling matrix \f$ \tilde{D}^{-1}_d[nactdof(i,p),nactdof(j,q)]\f$ for all active degrees od freedoms
     *  [in] irowt		field containing row numbers of aut
     *  [in] jqt	        pointer into field irowt
     *  [in] aut		transformation matrix \f$ T[p,q]\f$ for slave nodes \f$ p,q \f$ 
     *  [in] irowtinv	field containing row numbers of autinv
     *  [in] jqtinv	pointer into field irowtinv
     *  [in] autinv	transformation matrix \f$ T^{-1}[p,q]\f$ for slave nodes \f$ p,q \f$  
     *  [in] slavnor		slave normals
     *  [in] slavtan		slave tangents 
     *  [out]	iflagact	flag indicating if semi-smooth Newton has converged
     *  [in,out] cstress	current Lagrange multiplier 
     *  [in] cstressini	Lagrange multiplier at start of the increment
     *  [out] cdisp		vector saving contact variables for frd-output
     *  [out] f_cs            contact forces for active degrees of freedom
     *  [out] f_cm            not used any more
     *  [out] bp		current friction bound
     *  [out] cfs 		contact force 
     *  [out] cfm 		not used any more
     *  [in] islavnodeinv     (i) slave node index for node i
     *  [out] Bd		coupling matrix \f$ B_d[p,q]=\int \psi_p \phi_q dS \f$, \f$ p \in S, q \in M \f$ 
     *  [out] irowb		field containing row numbers of Bd
     *  [out] jqb		pointer into field irowb
     *  [out] Dd		coupling matrix \f$ D_d[p,q]=\int \psi_p \phi_q dS \f$, \f$ p,q \in S \f$ 
     *  [out] irowd		field containing row numbers of Dd
     *  [out] jqd		pointer into field irowd
     *  [out] Ddtil		coupling matrix \f$ \tilde{D}_d[p,q]=\int \psi_p \tilde{\phi}_q dS \f$, \f$ p,q \in S \f$ 
     *  [out] irowdtil	field containing row numbers of Ddtil
     *  [out] jqdtil		pointer into field irowdtil 
     */

void stressmortar(double *bhat,double *adc,double *auc,ITG *jqc,ITG *irowc,
		  ITG *neq,double *gap,double *b,ITG *islavact,
		  ITG *irowddinv,ITG *jqddinv,double *auddinv,ITG *irowt,
		  ITG *jqt,double *aut,ITG *irowtinv,ITG *jqtinv,
		  double *autinv,ITG *ntie,ITG *nslavnode,ITG *islavnode,
		  ITG *nmastnode,ITG *imastnode,double *slavnor,
		  double *slavtan,ITG *nactdof,ITG *iflagact,double *cstress,
		  double *cstressini,ITG *mi,double *cdisp,double *f_cs,
		  double *f_cm,ITG *iit,ITG *iinc,double *vold,double *vini,
		  double* bp,ITG *nk,ITG *nboun,ITG *ndirboun,
		  ITG *nodeboun,double *xboun,
		  ITG *nmpc,ITG *ipompc,ITG *nodempc,
		  double *coefmpc,
		  ITG *nslavmpc,ITG *islavmpc,
		  char *tieset,double  *elcon,double *tietol,
		  ITG *ncmat_,ITG *ntmat_,double *plicon,ITG *nplicon,
		  ITG *npmat_,ITG *nelcon,double *dtime,double *cfs,
		  double *cfm,ITG *islavnodeinv,double *Bd,ITG *irowb,
		  ITG *jqb,double *Dd,ITG *irowd,ITG *jqd,double *Ddtil,
		  ITG *irowdtil,ITG *jqdtil,double *Bdtil,ITG *irowbtil,
		  ITG *jqbtil,
		  ITG *nmethod,double *bet,
		  ITG *ithermal,ITG *iperturb,
		  char *labmpc,double *cam,double *veold,
		  double *accold,double *gam,double *cfsini,
		  double *cfstil,double *plkcon,ITG *nplkcon,char *filab,
		  double *f,double *fn,double *qa,ITG *nprint,char *prlab,
		  double *xforc,ITG *nforc){
  
  ITG i,j,l,jj,k,idof1,idof2,idof3,nodes,mt=mi[1]+1,nstick=0,nslip=0,ninacti=0,
    nnogap=0,nolm=0,ndiverg,nhelp,idof,node2,dirind,dirdep,index,ist,
    iout,num_cpus=1,keepset,derivmode,regmode,
    ndof,calcul_fn,calcul_f;
  
  double aux,stressnormal,ddispnormal,*unitmatrix=NULL,constant=1.E10,
    constantn=1.E10,constantt=1.E10,stresst[2],stressinit[2],stresstildet[2],
    disp_tildet[2],nw_t=0.0,*du=NULL,
    ch,coefdep,disp_t[2],scal,w_t[3],*bhat2=NULL,
    *fmpc=NULL,lm_t1_av,lm_t2_av,*rc=NULL,f_cs_tot[4],f_cm_tot[4],
    *vectornull=NULL,*u_oldt=NULL,
    *cstress2=NULL,*cstresstil=NULL,*cstressini2=NULL,*u_old=NULL,
    aninvloc,atauinvloc,gnc,gtc[2],*cold=NULL,*cold2=NULL,ln_old,
    ncf_n,max_ncf_n,ncf_t[2],max_ncf_t[2],mu,mumax,p0,beta,*b2=NULL,alpha;
  
  alpha=1-2*sqrt(*bet);
  keepset=0;
  mumax=-1.0;
  max_ncf_n=-1.0;max_ncf_t[0]=-1.0;max_ncf_t[1]=-1.0;
  lm_t1_av=0;lm_t2_av=0;
  *iflagact=0;
  
  NNEW(vectornull,double,neq[1]);
  NNEW(cstress2,double,neq[1]); 
  NNEW(bhat2,double,mt**nk);
  NNEW(b2,double,mt**nk);
  NNEW(cold,double,3*nslavnode[*ntie]);
  NNEW(cold2,double,3*nslavnode[*ntie]);
  NNEW(cstresstil,double,mt*nslavnode[*ntie]);
  NNEW(cstressini2,double,mt*nslavnode[*ntie]);
  NNEW(du,double,3*nslavnode[*ntie]);
  NNEW(u_old,double,3*nslavnode[*ntie]);
  NNEW(u_oldt,double,3*nslavnode[*ntie]);
  ndiverg=14;
  
  /* get cstress2 (Sitzmann,Equation (4.14)) */
  
  FORTRAN(opnonsym,(&neq[1],&aux,b,vectornull,adc,auc,jqc,irowc));
  for(i=0;i<neq[1];i++){f_cs[i]=bhat[i]-vectornull[i];}

  for(i=0;i<neq[1];i++){vectornull[i]=0.0;}	
  FORTRAN(opnonsymt,(&neq[1],&aux,f_cs,cstress2,vectornull,auddinv,jqddinv,
		     irowddinv));
  for(i=0;i<neq[1];i++){f_cs[i]=0.0;}
  
  NNEW(unitmatrix,double,neq[1]);
  
  // mechanical part
  
  for(i=0;i<neq[0];i++){
    unitmatrix[i]=1.;vectornull[i]=0.0;
    if(*nmethod==4){
      cstress2[i]*=1.0/(1.0+alpha);
    }
  }
  
  // thermal part
  
  for(i=neq[0];i<neq[1];i++){
    unitmatrix[i]=1.;vectornull[i]=0.0;
  }  
  for(i=0;i<mt**nk;i++){bhat2[i]=0.0;}
  iout=1;
  
  /* fill in missing results from SPCs/MPCs */
  
  FORTRAN(resultsini_mortar,(nk,b2,ithermal,iperturb,nactdof,&iout,vold,b,
			     nodeboun,ndirboun,xboun,nboun,ipompc,
			     nodempc,coefmpc,labmpc,nmpc,nmethod,cam,bet,
			     gam,dtime,mi));
  
  /* transformation delta tildeu^q->delta u:
     u=T tildeu^q*/
  
  for(i=0;i<*nk;i++){
    nodes=i+1;
    for(jj=jqt[nodes-1]-1;jj<jqt[nodes-1+1]-1;jj++){
      for(l=0;l<3;l++){
	idof1=mt*nodes-3+l;
	idof2=mt*irowt[jj]-3+l;
	bhat2[idof2]+=aut[jj]*b2[idof1];
      }
    }	
  }
  
  /* overwrite b with untransformed delta u*/
  
  for( i=0;i<*ntie;i++){      	
    if(tieset[i*(81*3)+80]=='C'){
      nhelp=0;
      for(j=nslavnode[i];j<nslavnode[i+1];j++){
	nodes=islavnode[j];
	for(l=0;l<3;l++){
	  b2[mt*nodes-3+l]=bhat2[mt*nodes-3+l];
	  idof1=nactdof[mt*nodes-3+l]-1;
	  if(idof1>-1){	      
	    if(*nmethod==4){
	      b[idof1]=(b2[mt*nodes-3+l])/(*bet**dtime**dtime);
	    }else{
	      b[idof1]=b2[mt*nodes-3+l];
	    }
	  }
	}
	if(islavact[j]>-1){nhelp++;}
      }
      ndiverg=max(ndiverg,(nhelp/100)+*ntie);
    }
  }
    
  /* normal formulation */  
  /* get du^hat */ 
  /* get uhat_k-1 */
    
  for(i=0;i<*nk;i++){
    nodes=i+1;
    for(jj=jqd[nodes-1]-1;jj<jqd[nodes-1+1]-1;jj++){
      for(k=0;k<3;k++){
	du[(islavnodeinv[irowd[jj]-1]-1)*3+k]+=Dd[jj]*b2[mt*nodes-3+k];
	u_oldt[(islavnodeinv[irowd[jj]-1]-1)*3+k]+=
	  Dd[jj]*(vold[mt*(nodes)-3+k]-vini[mt*(nodes)-3+k]);
	u_old[(islavnodeinv[irowd[jj]-1]-1)*3+k]+=
	  Dd[jj]*(vold[mt*(nodes)-3+k]);
      }
    }	    
  }
  for(i=0;i<*nk;i++){
    nodes=i+1;
    for(jj=jqb[nodes-1]-1;jj<jqb[nodes-1+1]-1;jj++){
      for(k=0;k<3;k++){
	du[(islavnodeinv[irowb[jj]-1]-1)*3+k]+=Bd[jj]*b2[mt*nodes-3+k];	
	u_oldt[(islavnodeinv[irowb[jj]-1]-1)*3+k]+=
	  Bd[jj]*(vold[mt*(nodes)-3+k]-vini[mt*(nodes)-3+k]);
	u_old[(islavnodeinv[irowb[jj]-1]-1)*3+k]+=
	  Bd[jj]*(vold[mt*(nodes)-3+k]);
      }
    }	  
  } 
  
  /* get lambda_scaled */
  
  for (i=0;i<*ntie;i++){
    if(tieset[i*(81*3)+80]=='C'){      
      for(j=nslavnode[i];j<nslavnode[i+1];j++){
	cold[3*j+0]=cstress[mt*j];
	cold[3*j+1]=cstress[mt*j+1];
	cold[3*j+2]=cstress[mt*j+2];      
	nodes=islavnode[j];	    	
	idof1=nactdof[mt*nodes-3]-1;	
	idof2=nactdof[mt*nodes-2]-1;	
	idof3=nactdof[mt*nodes-1]-1;
	ndof=0;
	if(idof1>-1){ndof++;}if(idof2>-1){ndof++;}if(idof3>-1){ndof++;}
	for(k=0;k<3;k++){	  
	  idof=nactdof[mt*nodes-3+k]-1;	  
	  if(idof>-1) {	  	  	    	  	    
	    cstress[mt*j+k]=cstress2[idof];	  	       	 	  
	  }else{	    	  	    
	    cstress[mt*j+k]=0.0;	     	  	    
	    for(jj=nslavmpc[2*j];jj<nslavmpc[2*(j)+1];jj++){
	      ist=islavmpc[jj];           	    	      
	      dirdep=nodempc[3*(ist-1)+1];           	    	      
	      coefdep=coefmpc[ist-1];           	    	      
	      index=nodempc[3*(ist-1)+2];	     	    	      
	      if(dirdep==k+1){               	      		
		while(index!=0){               				  
		  dirind=nodempc[3*(index-1)+1];
		  node2=nodempc[3*(index-1)];
		  ch=0.0;
		  if(node2>*nk){
		    idof=-1;
		  }else{
		    idof=nactdof[mt*node2-3+(dirind-1)]-1;
		  }
		  if(idof>-1){	       		 		    
		    ch=cstress2[idof];
		  } 	       
		  cstress[mt*j+dirdep-1]=cstress[mt*j+dirdep-1]-
		    coefmpc[index-1]*ch/coefdep;
		  index=nodempc[3*(index-1)+2];	       	      		
		}	      	    	      
	      }	     	  	    
	    }	     	     	  	   		  
	  }		
	}
      }
    }
  }
  
  /** generate cstressini2,cstresstil **/
  
  for (i=0;i<*ntie;i++){      	
    if(tieset[i*(81*3)+80]=='C'){      	
      for(j=nslavnode[i];j<nslavnode[i+1];j++){	    
	nodes=islavnode[j];
	for(jj=jqdtil[nodes-1]-1;jj<jqdtil[nodes-1+1]-1;jj++){
	  for(l=0;l<3;l++){
	    cstresstil[(islavnodeinv[irowdtil[jj]-1]-1)*mt+l]+=
	      Ddtil[jj]*cstress[(islavnodeinv[nodes-1]-1)*mt+l];
	    cstressini2[(islavnodeinv[irowdtil[jj]-1]-1)*mt+l]+=
	      Ddtil[jj]*cstressini[(islavnodeinv[nodes-1]-1)*mt+l];
	    cold2[(islavnodeinv[irowdtil[jj]-1]-1)*3+l]+=
	      Ddtil[jj]*cold[(islavnodeinv[nodes-1]-1)*3+l];
	  }
	}	  
      }	
    }    
  } 
  
  /* update slave nodes according to semi-smooth Newton */
  
  for (i=0;i<*ntie;i++){
    if(tieset[i*(81*3)+80]=='C'){      
      for(j=nslavnode[i];j<nslavnode[i+1];j++){
	nodes=islavnode[j];
	idof1=nactdof[mt*nodes-3]-1;	
	idof2=nactdof[mt*nodes-2]-1;	
	idof3=nactdof[mt*nodes-1]-1;
	ndof=0;
	if(idof1>-1){ndof++;}if(idof2>-1){ndof++;}if(idof3>-1){ndof++;}
	FORTRAN(getcontactparams,(&mu,&regmode,&aninvloc,&atauinvloc,
				  &p0,&beta,tietol,elcon,&i,ncmat_,ntmat_));
	mumax=max(mu,mumax);
	stressnormal=cstresstil[mt*j+0]*slavnor[3*j]+
	  cstresstil[mt*j+1]*slavnor[3*j+1]+cstresstil[mt*j+2]*slavnor[3*j+2];	
	stresst[0]=cstresstil[mt*j+0]*slavtan[6*j]+
	  cstresstil[mt*j+1]*slavtan[6*j+1]+cstresstil[mt*j+2]*slavtan[6*j+2];	
	stresst[1]=cstresstil[mt*j+0]*slavtan[6*j+3]+
	  cstresstil[mt*j+1]*slavtan[6*j+4]+cstresstil[mt*j+2]*slavtan[6*j+5];
	stressinit[0]=cstressini2[mt*j+0]*slavtan[6*j]+
	  cstressini2[mt*j+1]*slavtan[6*j+1]+
	  cstressini2[mt*j+2]*slavtan[6*j+2];	
	stressinit[1]=cstressini2[mt*j+0]*slavtan[6*j+3]+
	  cstressini2[mt*j+1]*slavtan[6*j+4]+
	  cstressini2[mt*j+2]*slavtan[6*j+5]; 
	stresstildet[0]=(stresst[0]-stressinit[0]);
	stresstildet[1]=(stresst[1]-stressinit[1]);
	ddispnormal=du[3*j+0]*slavnor[3*j]+du[3*j+1]*slavnor[3*j+1]+
	  du[3*j+2]*slavnor[3*j+2];	
	disp_tildet[0]=(du[3*j+0]+u_oldt[j*3])*slavtan[6*j]
	  +(du[3*j+1]+u_oldt[j*3+1])*slavtan[6*j+1]
	  +(du[3*j+2]+u_oldt[j*3+2])*slavtan[6*j+2];
	disp_tildet[1]=(du[3*j+0]+u_oldt[j*3])*slavtan[6*j+3]
	  +(du[3*j+1]+u_oldt[j*3+1])*slavtan[6*j+4]
	  +(du[3*j+2]+u_oldt[j*3+2])*slavtan[6*j+5];
	disp_t[0]=(du[3*j+0]+u_old[j*3])*slavtan[6*j]
	  +(du[3*j+1]+u_old[j*3+1])*slavtan[6*j+1]
	  +(du[3*j+2]+u_old[j*3+2])*slavtan[6*j+2];
	disp_t[1]=(du[3*j+0]+u_old[j*3])*slavtan[6*j+3]
	  +(du[3*j+1]+u_old[j*3+1])*slavtan[6*j+4]
	  +(du[3*j+2]+u_old[j*3+2])*slavtan[6*j+5];		
	ln_old=cold2[3*j+0]*slavnor[3*j]+cold2[3*j+1]*slavnor[3*j+1]+
	  cold2[3*j+2]*slavnor[3*j+2];
	
	double gnc_old,dgnc_old;
	derivmode=0;
	if(islavact[j]>-1){scal=Ddtil[jqdtil[nodes-1]-1];}else{scal=0.0;}
	FORTRAN(regularization_gn_c,(&ln_old,&derivmode,&regmode,&gnc_old,
				     &aninvloc,&p0,&beta,elcon,nelcon,&i,
				     ntmat_,plicon,nplicon,npmat_,ncmat_,
				     tietol,&scal));
	derivmode=1;
	FORTRAN(regularization_gn_c,(&ln_old,&derivmode,&regmode,&dgnc_old,
				     &aninvloc,&p0,&beta,elcon,nelcon,&i,ntmat_,
				     plicon,nplicon,npmat_,ncmat_,tietol,
				     &scal));
	derivmode=0;
	FORTRAN(regularization_gn_c,(&stressnormal,&derivmode,&regmode,&gnc,
				     &aninvloc,&p0,&beta,elcon,nelcon,&i,ntmat_,
				     plicon,nplicon,npmat_,ncmat_,tietol,
				     &scal));
	
	derivmode=0;
	FORTRAN(regularization_gt_c,(stresstildet,&derivmode,gtc,
				     &atauinvloc));
	
	constantt=min(constant,1.0/atauinvloc);	
	constantn=constant;
	
	if(mu>1.E-10){ 	  
	  bp[j]=mu*(stressnormal+constantn*(ddispnormal-gap[j]-gnc));	
	}else{	  
	  bp[j]=(stressnormal+constantn*(ddispnormal-gap[j]-gnc));	
	}
	
	w_t[0]=stresst[0]+constantt*(disp_tildet[0]-gtc[0]);	
	w_t[1]=stresst[1]+constantt*(disp_tildet[1]-gtc[1]);		
	nw_t=sqrt(w_t[0]*w_t[0]+w_t[1]*w_t[1]); 
	ncf_t[0]=0.0;ncf_t[1]=0.0;ncf_n=0.0;
	
	/* evaluate non-linear complementary functions 
	   (Sitzmann,Equation (3.54),(3.65)) */
	
	if(mu>1.E-10){
	  ncf_n=mu*stressnormal-max(0.0,bp[j]);
	  if(islavact[j]==0){
	    ncf_t[0]=stresst[0];
	    ncf_t[1]=stresst[1];
	  }else if(islavact[j]==1){
	    ncf_t[0]=disp_tildet[0]-gtc[0];
	    ncf_t[1]=disp_tildet[1]-gtc[1];
	  }else if(islavact[j]==2){
	    ncf_t[0]=stresst[0]-bp[j]*w_t[0]/nw_t;
	    ncf_t[1]=stresst[1]-bp[j]*w_t[1]/nw_t;
	  }else{
	    ncf_t[0]=0.0;ncf_t[1]=0.0;
	  }
	}else{
	  ncf_n=stressnormal-max(0.0,bp[j]);
	  ncf_t[0]=0.0;ncf_t[1]=0.0;
	}  
	if( ncf_n<0.0)ncf_n=-ncf_n;
	if( ncf_t[0]<0.0)ncf_t[0]=-ncf_t[0];
	if( ncf_t[1]<0.0)ncf_t[1]=-ncf_t[1];
	if(ncf_n>max_ncf_n && ndof>0 && islavact[j]>-1){max_ncf_n=ncf_n;}
	if(ncf_t[0]>max_ncf_t[0] && ndof>0 && islavact[j]>-1){
	  max_ncf_t[0]=ncf_t[0];}
	if(ncf_t[1]>max_ncf_t[1] && ndof>0 && islavact[j]>-1){
	  max_ncf_t[1]=ncf_t[1];}	
	
	if(keepset==0){
	  if(mu>1.E-10){ //contact tie with friction	
	    if (islavact[j]>-1){	
		
	      /*** Active/Inactive set condition cf: 
		   Sitzmann Equation (3.58),(3.70) ***/
		
	      if(bp[j]>-1E-10 && nw_t<(bp[j])){
		nstick++;			
		if (islavact[j]!=1) {*iflagact=1;}				
		islavact[j]=1;
		cdisp[6*j]=gap[j]-ddispnormal;		  
		cdisp[6*j+1]=disp_t[0];
		cdisp[6*j+2]=disp_t[1];
		cdisp[6*j+3]=stressnormal;				      
		cdisp[6*j+4]=stresst[0];				      
		cdisp[6*j+5]=stresst[1];
		lm_t1_av=lm_t1_av+abs(stresst[0]);
		lm_t2_av=lm_t2_av+abs(stresst[1]);
	      }else if(bp[j]>-1E-10 && nw_t>=(bp[j])){
		nslip++;						      
		if (islavact[j]!=2) {*iflagact=1;}	
		islavact[j]=2;
		cdisp[6*j]=gap[j]-ddispnormal;		  
		cdisp[6*j+1]=disp_t[0];
		cdisp[6*j+2]=disp_t[1];
		cdisp[6*j+3]=stressnormal;				      
		cdisp[6*j+4]=stresst[0];				      
		cdisp[6*j+5]=stresst[1];
		lm_t1_av=lm_t1_av+abs(stresst[0]);
		lm_t2_av=lm_t2_av+abs(stresst[1]);		  
	      }else{				      
		if (islavact[j]>0){ *iflagact=1;}	
		ninacti++;
		islavact[j]=0;				      
		cstress[mt*j+0]=0.;				      
		cstress[mt*j+1]=0.;				      
		cstress[mt*j+2]=0.;
		cdisp[6*j]=gap[j]-ddispnormal;
		cdisp[6*j+1]=0.;		        	      
		cdisp[6*j+2]=0.;		        	      
		cdisp[6*j+3]=0.;		        	      
		cdisp[6*j+4]=0.;		        	      
		cdisp[6*j+5]=0.;
		if(idof1>-1)cstress2[idof1]=0.0;
		if(idof2>-1)cstress2[idof2]=0.0;
		if(idof3>-1)cstress2[idof3]=0.0;
	      }
	    }else{
	      
	      /* nodes without LM contribution */
	      
	      if(islavact[j]==-1){	      
		nnogap++;		    
	      }else if(islavact[j]==-2){	      
		nolm++;	    
	      }	            
	      cstress[mt*j+0]=0.;		    
	      cstress[mt*j+1]=0.;		    
	      cstress[mt*j+2]=0.;
	      cdisp[6*j]=0.;		    	    
	      cdisp[6*j+1]=0.;		    	    
	      cdisp[6*j+2]=0.;		    	    
	      cdisp[6*j+3]=0.;		    	    
	      cdisp[6*j+4]=0.;		    	    
	      cdisp[6*j+5]=0.;
	      if(idof1>-1)cstress2[idof1]=0.0;
	      if(idof2>-1)cstress2[idof2]=0.0;
	      if(idof3>-1)cstress2[idof3]=0.0;
	    }   	
	  }else{ //no friction            	  
	    if (islavact[j]>-1){
	      
	      /*** Active/Inactive set condition cf: 
		   Sitzmann Equation (3.58) ***/
	      
	      if(bp[j]>-1E-10){	      
		nslip++;				      
		if (islavact[j]!=2) {*iflagact=1;}
		islavact[j]=2;
		cdisp[6*j]=gap[j]-ddispnormal;
		cdisp[6*j+1]=disp_t[0];
		cdisp[6*j+2]=disp_t[1];						
		cdisp[6*j+3]=stressnormal;			
		cdisp[6*j+4]=stresst[0];			
		cdisp[6*j+5]=stresst[1];			    	    
	      }else{				          
		if (islavact[j]!=0){ *iflagact=1;}
		ninacti++;				      
		islavact[j]=0;
		cstress[mt*j+0]=0.;
		cstress[mt*j+1]=0.;				      
		cstress[mt*j+2]=0.;
		cdisp[6*j]=ddispnormal-gap[j];		        	      
		cdisp[6*j+1]=0.;	      
		cdisp[6*j+2]=0;
		cdisp[6*j+3]=0.;		        	      
		cdisp[6*j+4]=0.;		        	      
		cdisp[6*j+5]=0.;
		if(idof1>-1)cstress2[idof1]=0.0;
		if(idof2>-1)cstress2[idof2]=0.0;
		if(idof3>-1)cstress2[idof3]=0.0;
	      }			  
	    }else{
	      
	      /* nodes without LM contribution */
	      
	      if(islavact[j]==-1){	      
		nnogap++;		    
	      }else{	      
		nolm++;	    
	      }	            	    
	      cstress[mt*j+0]=0.;		    	    
	      cstress[mt*j+1]=0.;		    	    
	      cstress[mt*j+2]=0.;
	      cdisp[6*j]=0.;			    	    
	      cdisp[6*j+1]=0.;           
	      cdisp[6*j+2]=0.0;			    		    	    
	      cdisp[6*j+3]=0.;		    	    
	      cdisp[6*j+4]=0.;		    	    
	      cdisp[6*j+5]=0.;
	      if(idof1>-1)cstress2[idof1]=0.0;
	      if(idof2>-1)cstress2[idof2]=0.0;
	      if(idof3>-1)cstress2[idof3]=0.0;
	    }			    	
	  }	
	}else{
	  
	  /* in case active set is fixed */
	  
	  if(islavact[j]>0){				     
	    cdisp[6*j]=gap[j]-ddispnormal;
	    cdisp[6*j+1]=disp_t[0];
	    cdisp[6*j+2]=disp_t[1];	      			
	    cdisp[6*j+3]=stressnormal;			
	    cdisp[6*j+4]=stresst[0];			
	    cdisp[6*j+5]=stresst[1];		   
	    
	  }else{
	    cstress[mt*j+0]=0.;			
	    cstress[mt*j+1]=0.;			
	    cstress[mt*j+2]=0.;
	    cdisp[6*j]=ddispnormal-gap[j];		        
	    cdisp[6*j+1]=0.;
	    cdisp[6*j+2]=0;						        
	    cdisp[6*j+3]=0.;		        
	    cdisp[6*j+4]=0.;		        
	    cdisp[6*j+5]=0.;
	    if(idof1>-1)cstress2[idof1]=0.0;
	    if(idof2>-1)cstress2[idof2]=0.0;
	    if(idof3>-1)cstress2[idof3]=0.0;
	  }	  	
	}
	
	/* update weighted dual gap */
	
	gap[j]=gap[j]-ddispnormal;     
      }       
    }   
  }
  
  if(keepset==0){
    
    /* relative convergence criteria for semi-smooth Newton */
    
    if(max_ncf_n>1.e-3){*iflagact=1;}
    if((max_ncf_t[0]/((lm_t1_av+0.001)/(nstick+nslip+0.001))>9.e-4 ||
	max_ncf_t[1]/((lm_t2_av+0.001)/(nstick+nslip+0.001))>9.e-4) &&
       mumax>1.E-10 ){*iflagact=1;} 
  }
  
  SFREE(unitmatrix);
  
  /* calculating the  contact forces 
     Ph.D. Thesis Stefan Hartmann eqn. (6.26) */
  // fill cstress for nogap nodes
  
  for(i=0;i<mt**nk;i++){cfs[i]=0.0;}
  if(*nmethod==4){for(i=0;i<mt**nk;i++){cfstil[i]=0.0;}}
  
  // get contact forces
  
  for(i=0;i<*nk;i++){
    for(j=jqd[i]-1;j<jqd[i+1]-1;j++){
      for(l=0;l<3;l++){
	cfs[mt*(i+1)-3+l]+=Dd[j]*cstress[mt*(islavnodeinv[irowd[j]-1]-1)+l];
      }
    }
  }
  for(i=0;i<*nk;i++){
    for(j=jqb[i]-1;j<jqb[i+1]-1;j++){
      for(l=0;l<3;l++){
	cfs[mt*(i+1)-3+l]+=Bd[j]*cstress[mt*(islavnodeinv[irowb[j]-1]-1)+l];
      }
    }
  }
  if(*nmethod==4){
    for(i=0;i<*nk;i++){
      for(j=jqdtil[i]-1;j<jqdtil[i+1]-1;j++){
	for(l=0;l<3;l++){
	  cfstil[mt*(i+1)-3+l]+=Ddtil[j]*
	    cstress[mt*(islavnodeinv[irowdtil[j]-1]-1)+l];
	}
      }
    }
    for(i=0;i<*nk;i++){
      for(j=jqbtil[i]-1;j<jqbtil[i+1]-1;j++){
	for(l=0;l<3;l++){
	  cfstil[mt*(i+1)-3+l]+=Bdtil[j]*
	    cstress[mt*(islavnodeinv[irowbtil[j]-1]-1)+l];
	}
	
      }
    }
    
    for(i=0;i<mt**nk;i++){cfsini[i]=0.0;}
    for(i=0;i<*nk;i++){
      for(j=jqd[i]-1;j<jqd[i+1]-1;j++){
	for(l=0;l<3;l++){
	  cfsini[mt*(i+1)-3+l]+=
	    Dd[j]*cstressini[mt*(islavnodeinv[irowd[j]-1]-1)+l];
	}
      }
    }
    for(i=0;i<*nk;i++){
      for(j=jqb[i]-1;j<jqb[i+1]-1;j++){
	for(l=0;l<3;l++){
	  cfsini[mt*(i+1)-3+l]+=
	    Bd[j]*cstressini[mt*(islavnodeinv[irowb[j]-1]-1)+l];
	}
      }
    } 
    
  }
  
  NNEW(fmpc,double,*nmpc);
  NNEW(rc,double,*nk*mt);
  calcul_fn=1;
  calcul_f=1;
  for(i=0;i<*nk;i++){
    for(l=0;l<3;l++){      
      if(*nmethod==4){
	
	/* modification for dynamic calculations */
	/* Sitzmann Equation (5.12) */
	
	rc[mt*(i+1)-3+l]=cfs[mt*(i+1)-3+l]*(1+alpha)-cfsini[mt*(i+1)-3+l]*
	  (alpha);
      }else{
	rc[mt*(i+1)-3+l]=cfs[mt*(i+1)-3+l];
      }
    }
  }
  
  resultsforc(nk,f_cs,rc,nactdof,ipompc,nodempc,coefmpc,labmpc,nmpc,mi,fmpc,
	      &calcul_fn,&calcul_f,&num_cpus);
  SFREE(fmpc);
  
  /* print total contact force per contact tie for debugging purpose */
  
  for( i=0;i<*ntie;i++){
    if(tieset[i*(81*3)+80]=='C'){
      for(jj=0;jj<3;jj++){	
	f_cs_tot[jj]=0.0;	
	f_cm_tot[jj]=0.0;    
      }
      for(j=nslavnode[i];j<nslavnode[i+1];j++){	  
	nodes=islavnode[j];	  
	for(l=0;l<3;l++){	    
	  idof1=nactdof[mt*nodes-3+l]-1;
	  if(idof1>-1){
	    f_cs_tot[l]+=f_cs[idof1];
	    f_cs_tot[l]+=f_cm[idof1];
	  }
	}
      }			
      for(j=nmastnode[i];j<nmastnode[i+1];j++){	  
	nodes=imastnode[j];	    
	for(l=0;l<3;l++){		
	  idof1=nactdof[mt*nodes-3+l]-1;
	  if(idof1>-1){		 
	    f_cm_tot[l]+=f_cs[idof1];
	    
	  }	    
	}
      }
    } 	      
  }
  
  /* transform contact traction for frd-output I(LM)=T*LM */
  
  for( i=0;i<*ntie;i++){
    if(tieset[i*(81*3)+80]=='C'){
      for(j=nslavnode[i];j<nslavnode[i+1];j++){
	for(jj=0;jj<mt;jj++){	
	  cstresstil[mt*j+jj]=0.0;	    
	}
      }
    }
  }
  for(i=0;i<*nk;i++){
    nodes=i+1;
    for(jj=jqt[nodes-1]-1;jj<jqt[nodes-1+1]-1;jj++){
      for(l=0;l<3;l++){
	idof1=mt*(islavnodeinv[nodes-1]-1)+l;
	idof2=mt*(islavnodeinv[irowt[jj]-1]-1)+l;
	cstresstil[idof2]+=aut[jj]*cstress[idof1];
      }
    }	
  }
  for( i=0;i<*ntie;i++){      	
    if(tieset[i*(81*3)+80]=='C'){	        
      for(j=nslavnode[i];j<nslavnode[i+1];j++){
	stressnormal=cstresstil[mt*j+0]*slavnor[3*j]+
	  cstresstil[mt*j+1]*slavnor[3*j+1]+cstresstil[mt*j+2]*slavnor[3*j+2];	
	stresst[0]=cstresstil[mt*j+0]*slavtan[6*j]+
	  cstresstil[mt*j+1]*slavtan[6*j+1]+cstresstil[mt*j+2]*slavtan[6*j+2];	
	stresst[1]=cstresstil[mt*j+0]*slavtan[6*j+3]+
	  cstresstil[mt*j+1]*slavtan[6*j+4]+cstresstil[mt*j+2]*slavtan[6*j+5];
	cdisp[6*j+3]=stressnormal;			
	cdisp[6*j+4]=stresst[0];			
	cdisp[6*j+5]=stresst[1];
      }
    }
  }
  
  /* needed for adaptive time stepping */
  
  if(*iit>ndiverg){*iflagact=0;}
  
  SFREE(bhat2);
  SFREE(b2);
  SFREE(cstress2);
  SFREE(u_old);
  SFREE(u_oldt);
  SFREE(du);
  SFREE(cstresstil);
  SFREE(cstressini2);
  SFREE(cold);
  SFREE(cold2);
  SFREE(vectornull);
  SFREE(rc);
  return;
}
