/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2011 Guido Dhondt                     */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                                                                       */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "CalculiX.h"
#include "mortar.h"

/*
  Calculates the coupling matrices \f$ B_d \f$ and \f$ D_d \f$, and insert them into the data structure

  Authors: Samoela Rakotonanahary, Saskia Sitzmann */

void bdfill(ITG **irowbdp,ITG *jqbd,double **aubdp,ITG *nzsbd,
	    ITG **irowbdtilp,ITG *jqbdtil,double **aubdtilp,ITG *nzsbdtil,
	    ITG **irowbdtil2p,ITG *jqbdtil2,double **aubdtil2p,ITG *nzsbdtil2,
	    ITG **irowddp,ITG *jqdd,double **auddp,
	    ITG **irowddtilp,ITG *jqddtil,double **auddtilp,
	    ITG **irowddtil2p,ITG *jqddtil2,double **auddtil2p,
	    ITG **irowddinvp,ITG *jqddinv,double **auddinvp,
	    ITG *irowt,ITG *jqt,double *aut,
	    ITG *irowtinv,ITG *jqtinv,double *autinv,
	    ITG *ntie,ITG *ipkon,ITG *kon,
	    char *lakon,ITG *nslavnode,ITG *nmastnode,ITG *imastnode,
	    ITG *islavnode,ITG *islavsurf,ITG *imastsurf,double *pmastsurf,
	    ITG *itiefac,char *tieset,ITG *neq,ITG *nactdof,double *co,
	    double *vold,
	    ITG *iponoels,ITG *inoels,ITG *mi,double *gapmints,double *gap,
	    double* pslavsurf,double* pslavdual,
	    ITG *nintpoint,double *slavnor,ITG *nk,
	    ITG *nmpc,ITG *ipompc,ITG *nodempc,double *coefmpc,
	    ITG *ikmpc,ITG *ilmpc,ITG *nslavmpc,ITG *islavmpc,
	    ITG *nmastmpc,ITG *imastmpc,
	    ITG *iit,ITG *iinc,ITG *islavactdof,ITG *islavact,ITG *islavnodeinv,
	    double **Bdp,ITG **irowbp,ITG *jqb,
	    double **Bdhelpp,ITG **irowbhelpp,ITG *jqbhelp,
	    double **Ddp,ITG **irowdp,ITG *jqd,
	    double **Ddtilp,ITG **irowdtilp,ITG *jqdtil,
	    double **Bdtilp,ITG **irowbtilp,ITG *jqbtil,
	    ITG *ithermal){
  
  ITG i,ii,j,jj,k,kk,l,ll,icounter,icounter2,idof1,idofs,idofm,nodesf,nodem,
    isn,imn,*mast1=NULL,nzsbhelp,idof0,*irowbd=NULL,*irowdd=NULL,
    *irowddinv=NULL,*irowddtil=NULL,*irowbdtil=NULL,*irowbdtil2=NULL,
    *irowddtil2=NULL,ifree,ifree1,ifree3,mt=mi[1]+1,istart,*iscontr=NULL,
    *imcontr=NULL,*igap=NULL,*idcontr1=NULL,*idcontr2=NULL,*igcontr=NULL,
    intpointl,*mast2=NULL,ifree2,ist,ist2,node1,dim,index,dirdep,dirind,index2,
    dirdep2,dirind2,node2,idof2,idof3,nzsddinv,nzsdd,*irowb=NULL,nzsb1,
    nzsbtil2,nzsdtil,nzsd,nzsb,nzsdinv,nzsbtil,nzsd1,nzsb2a,*mast3=NULL,
    *irowdtil=NULL,*irowb1=NULL,*jqb1=NULL,*irowb2=NULL,*jqb2=NULL,
    *irowbtil2=NULL,*jqbtil2=NULL,*irowd=NULL,*irowbtil=NULL,*irowdinv=NULL,
    *jqdinv=NULL,*jqb2a=NULL,*irowb2a=NULL,*jqd1=NULL,*irowd1=NULL,
    *irowbhelp=NULL,*irowbtil2t=NULL,*jqbtil2t=NULL,*irowdtil2=NULL,
    *jqdtil2=NULL,nzsddtil2,nzsdtil2,islavnodeentry,jrow,iadd;
  
  double contribution=0.0,*aubd=NULL,*audd=NULL,*auddtil=NULL,*auddtil2=NULL,
    *auddinv=NULL,*Bdhelp=NULL,*contr=NULL,*dcontr=NULL,*gcontr=NULL,*Bd=NULL,
    *aubdtil=NULL,*aubdtil2=NULL,c2,coefdep,coefdep2,c3,*Dd=NULL,*Ddinv=NULL,
    *Ddtil2=NULL,*dinvloc=NULL,*dloc=NULL,*Bd1=NULL,*Bd2=NULL,*Bdtil2=NULL,
    *Bdtil=NULL,*Ddtil=NULL,*Dd1=NULL,*Bd2a=NULL,*Bdtil2t=NULL,detdloc,e1,e2,
    e3;
  
  /*
   * calculate Dtilde_d,B_d,Btilde_d2 nodewise
   * perform trafo Ntilde to N
   * go from nodes to dofs considering SPC's and MPC's
   */
  
  irowbd=*irowbdp; aubd=*aubdp;
  irowdd=*irowddp; audd=*auddp;  
  irowddtil=*irowddtilp; auddtil=*auddtilp;
  irowddtil2=*irowddtil2p; auddtil2=*auddtil2p;
  irowddinv=*irowddinvp; auddinv=*auddinvp;  
  irowbdtil=*irowbdtilp; aubdtil=*aubdtilp;
  irowbdtil2=*irowbdtil2p; aubdtil2=*aubdtil2p;
  irowb =*irowbp; Bd=*Bdp;
  irowbhelp =*irowbhelpp; Bdhelp=*Bdhelpp;
  irowd =*irowdp; Dd=*Ddp;
  irowdtil =*irowdtilp; Ddtil=*Ddtilp;
  irowbtil =*irowbtilp; Bdtil=*Bdtilp;
  ifree=1;
  // position in the fields FORTRAN condition
  ifree1=1;ifree2=1;ifree3=1;
  
  NNEW(igap,ITG,nslavnode[*ntie]);
  NNEW(dinvloc,double,9*nslavnode[*ntie]);
  NNEW(dloc,double,9*nslavnode[*ntie]);
  
  nzsdd=3*nslavnode[*ntie];
  nzsddinv=3*nslavnode[*ntie];
  
  nzsb1=*nzsbd;
  nzsbtil2=*nzsbd;
  NNEW(Bd1,double,nzsb1);
  NNEW(Bdtil2,double,nzsb1);
  NNEW(mast2,ITG,nzsb1);
  NNEW(mast3,ITG,nzsbtil2);
  nzsdtil=*nk;
  RENEW(Ddtil,double,*nk);
  RENEW(irowdtil,ITG,nzsdtil);
  NNEW(mast1,ITG,nzsdtil);
  NNEW(irowb1,ITG,nzsb1);
  NNEW(jqb1,ITG,*nk+1);
  NNEW(irowbtil2,ITG,nzsbtil2);
  NNEW(jqbtil2,ITG,*nk+1);
  
  /* calculating the off-diagonal terms and storing them in aubd */
  
  /* meaning of the fields in FORTRAN notation:
     ipointer(i): points to an element in field aubd belonging to column i
     aubd(ipointer(i)): value of that element
     irowbd(ipointer(i)): row to which that element belongs
     mast1(ipointer(i)): points to another element in field aubd belonging
     to column i, unless zero.
  */
  
  /* build Ddtil ,Bdtil2 ,Bd in nodes **/
  
  for( i=0; i<*ntie; i++){    
    if(tieset[i*(81*3)+80]=='C'){ 
      for(j=nslavnode[i]; j<nslavnode[i+1]; j++){
	gap[j]=0.0;
      }
      for(l= itiefac[2*i]; l<=itiefac[2*i+1];l++){
	intpointl=islavsurf[2*l+1]-islavsurf[2*(l-1)+1];
	if(intpointl>0){
	  NNEW(contr,double,9*9*intpointl/7+1);
	  NNEW(iscontr,ITG,9*9*intpointl/7+1);
	  NNEW(imcontr,ITG,9*9*intpointl/7+1);
	  NNEW(dcontr,double,9*9*intpointl/7+1);
	  NNEW(idcontr1,ITG,9*9*intpointl/7+1);
	  NNEW(idcontr2,ITG,9*9*intpointl/7+1);
	  NNEW(gcontr,double,9*intpointl/7+1);
	  NNEW(igcontr,ITG,9*intpointl/7+1);
	  
	  /* build dual coupling matrices (see phd-thesis Sitzmann equation
	     (4.18) and (4.19)) */
	  
	  FORTRAN(createbd,(&i,&l,ipkon,kon,lakon,co,vold,gapmints,islavsurf,
			    imastsurf,pmastsurf,contr,iscontr,imcontr,
			    dcontr,idcontr1,idcontr2,gcontr,igcontr,mi,
			    pslavsurf,pslavdual,nslavnode,islavnode,nmastnode,
			    imastnode,&icounter,&icounter2,islavact));
	  
	  /* build B_d **/
	  
	  for(j=0; j<icounter;j++){				
	    contribution=-contr[j];				
	    nodesf=islavnode[iscontr[j]-1];				
	    nodem=imastnode[imcontr[j]-1];				
	    
	    if((contribution>1e-14 ||contribution<-1e-14)){ 
		
	      insertas(&irowb1,&mast2,&nodesf,&nodem, &ifree2,&nzsb1,
		       &contribution,&Bd1);				
	    }				               		  
	  }
	  
	  /* build D_tilde **/
	  
	  for(j=0; j<icounter2;j++){			
	    contribution=dcontr[j];			
	    nodesf=islavnode[idcontr1[j]-1];		        
	    nodem=islavnode[idcontr2[j]-1];
	    
	    if(nodesf==nodem){
		
	      // Slave nodes with LM contribution
		
	      if(contribution>1e-14 ||contribution<-1e-14){
		
		insertas(&irowdtil,&mast1,&nodesf,&nodem, &ifree1,&nzsdtil,
			 &contribution,&Ddtil);	
	      }
	    }else{
		
	      // nogap or noLM slave nodes...
		
	      if(contribution>1e-14 ||contribution<-1e-14){	
		
		// Bbtil2^T
		// ... are stored in Bdtil
		
		insertas(&irowbtil2,&mast3,&nodem,&nodesf, &ifree3,&nzsbtil2,
			 &contribution,&Bdtil2);
	      }						
	    }		        
	  }
	  
	  /* build dual gap **/
	  
	  for(j=0; j<sqrt(icounter2);j++){					
	    contribution=gcontr[j];			
	    ll=igcontr[j]-1;			
	    gap[ll]+=contribution;
	    igap[ll]=1;
	  }
	  SFREE(contr);
	  SFREE(iscontr);
	  SFREE(imcontr);
	  SFREE(dcontr);
	  SFREE(idcontr1);
	  SFREE(idcontr2);
	  SFREE(gcontr);
	  SFREE(igcontr);              
	}
      }
    }
  }
  
  // set slave faces without corresponding master faces to normal nodes N
  
  for(i=0;i<*ntie;i++){	  
    if(tieset[i*(81*3)+80]=='C'){	    
      for(j=nslavnode[i];j<nslavnode[i+1];j++){
	if(igap[j]<1){
	  islavact[j]=-3;}
      }
    }
  }
  
  /* sort Ddtil **/
  
  nzsdtil=ifree1-1;
  RENEW(irowdtil,ITG,nzsdtil);
  RENEW(Ddtil,double,nzsdtil);
  dim=*nk;
  matrixsort(Ddtil,mast1,irowdtil,jqdtil,&nzsdtil,&dim);  

  /* adding values for identical locations */
  
  icounter=0;
  for(i=0;i<*nk;i++){
    if(jqdtil[i]!=jqdtil[i+1]){
      irowdtil[icounter]=irowdtil[jqdtil[i]-1];
      Ddtil[icounter]=Ddtil[jqdtil[i]-1];
      icounter++;
      istart=icounter;
      for(j=jqdtil[i];j<jqdtil[i+1]-1;j++){
	if(irowdtil[j]==irowdtil[icounter-1]){
	  Ddtil[icounter-1]+=Ddtil[j];   
	}else{
	  irowdtil[icounter]=irowdtil[j];
	  Ddtil[icounter]=Ddtil[j];
	  icounter++;
	}
      }
    }else{ istart=icounter+1;}
    
    jqdtil[i]=istart;
  }
  jqdtil[*nk]=icounter+1; 
  nzsdtil=icounter;
  RENEW(irowdtil,ITG,nzsdtil);
  RENEW(Ddtil,double,nzsdtil); 
  SFREE(mast1);
  
  /* sort Bd1 **/
  
  nzsb1=ifree2-1;
  dim=*nk;
  matrixsort(Bd1,mast2,irowb1,jqb1,&nzsb1,&dim); 

  /* adding values for identical locations */
  
  icounter=0;
  for(i=0;i<*nk;i++){
    if(jqb1[i]!=jqb1[i+1]){
      irowb1[icounter]=irowb1[jqb1[i]-1];
      Bd1[icounter]=Bd1[jqb1[i]-1];
      icounter++;
      istart=icounter;
      for(j=jqb1[i];j<jqb1[i+1]-1;j++){
	if(irowb1[j]==irowb1[icounter-1]){
	  Bd1[icounter-1]+=Bd1[j];   
	}else{
	  irowb1[icounter]=irowb1[j];
	  Bd1[icounter]=Bd1[j];
	  icounter++;
	}
      }
    }else{ istart=icounter+1;}
    
    jqb1[i]=istart;
  }
  jqb1[*nk]=icounter+1; 

  // guido: next line was replaced by nzsb1=icounter;
  //nzsb=icounter;

  nzsb1=icounter;
  RENEW(irowb1,ITG,nzsb1);
  RENEW(Bd1,double,nzsb1); 
  SFREE(mast2);
  
  /* sort Bdtil2 **/
  
  nzsbtil2=ifree3-1;
  dim=*nk;
  matrixsort(Bdtil2,mast3,irowbtil2,jqbtil2,&nzsbtil2,&dim); 

  /* adding values for identical locations */
  
  icounter=0;
  for(i=0;i<*nk;i++){
    if(jqbtil2[i]!=jqbtil2[i+1]){
      irowbtil2[icounter]=irowbtil2[jqbtil2[i]-1];
      Bdtil2[icounter]=Bdtil2[jqbtil2[i]-1];
      icounter++;
      istart=icounter;
      for(j=jqbtil2[i];j<jqbtil2[i+1]-1;j++){
	if(irowbtil2[j]==irowbtil2[icounter-1]){
	  Bdtil2[icounter-1]+=Bdtil2[j];   
	}else{
	  irowbtil2[icounter]=irowbtil2[j];
	  Bdtil2[icounter]=Bdtil2[j];
	  icounter++;
	}
      }
    }else{ istart=icounter+1;}
    
    jqbtil2[i]=istart;
  }
  jqbtil2[*nk]=icounter+1; 
  nzsbtil2=icounter;
  RENEW(irowbtil2,ITG,nzsbtil2);
  RENEW(Bdtil2,double,nzsbtil2); 
  SFREE(mast3);
  
  /* build D_d=(D_d^til*T^-1)=(T^-T*D_d^til)^T **/
  /* need D columnwise stored **/
  
  dim=*nk;
  nzsd1=jqtinv[*nk];
  NNEW(Dd1,double,nzsd1);
  NNEW(irowd1,ITG,nzsd1);
  NNEW(jqd1,ITG,*nk+1);
  kk=1;	
  jqd1[0]=1;
  for(ii=0;ii<*nk;ii++){
    for(jj=jqtinv[ii]-1;jj<jqtinv[ii+1]-1;jj++){
      k=irowtinv[jj]-1;
      if(jqdtil[k+1]-jqdtil[k]==1){
	Dd1[kk-1]=Ddtil[jqdtil[k]-1]*autinv[jj];
	irowd1[kk-1]=irowtinv[jj];
	kk++;
      }else{
	// something went wrong
      }
    }
    jqd1[ii+1]=kk;
  }
  nzsd1=kk;
  RENEW(Dd1,double,nzsd1);
  RENEW(irowd1,ITG,nzsd1);
    
  nzsb2a=nslavnode[*ntie];
  NNEW(Bd2a,double,nzsb2a);
  NNEW(irowb2a,ITG,nzsb2a);
  NNEW(jqb2a,ITG,*nk+1);	     
  multi_rect(Bdtil2,irowbtil2,jqbtil2,dim,dim,
	     autinv,irowtinv,jqtinv,dim,dim,           
	     &Bd2a,&irowb2a,jqb2a,&nzsb2a);
  
  NNEW(Bdtil2t,double,nzsbtil2);	
  NNEW(irowbtil2t,ITG,nzsbtil2);	
  NNEW(jqbtil2t,ITG,*nk+1);
  dim=*nk;
  transpose(Bdtil2,jqbtil2,irowbtil2,&dim,
	    Bdtil2t,jqbtil2t,irowbtil2t);
    
  nzsd=nzsd1;
  RENEW(Dd,double,nzsd);
  RENEW(irowd,ITG,nzsd);
  add_rect(Dd1,irowd1,jqd1,dim,dim,
	   Bd2a,irowb2a,jqb2a,dim,dim,
	   &Dd,&irowd,jqd,&nzsd);     

  // guido: next lines up to an including the print statements
  // are supersfluous: Bd1 de factor = Bd sind Bd2 is zero
  //
  
  nzsb=nslavnode[*ntie];
  RENEW(Bd,double,nzsb);
  RENEW(irowb,ITG,nzsb);
    
  NNEW(Bd2,double,1);
  NNEW(irowb2,ITG,1);
  NNEW(jqb2,ITG,*nk+1);  
  add_rect(Bd1,irowb1,jqb1,dim,dim,
	   Bd2,irowb2,jqb2,dim,dim,
	   &Bd,&irowb,jqb,&nzsb);
    
  /* build D_d^inv=((T^T)^T*D_dtil^-1) **/
  /* need D_d^inv column and row-wise stored? **/
  /* 1. invert D_dtil **/

  nzsdinv=nzsdtil;
  NNEW(Ddinv,double,nzsdinv);
  NNEW(irowdinv,ITG,nzsdinv);
  NNEW(jqdinv,ITG,*nk+1);
  ifree=1;ifree2=1;
  jqdinv[0]=1;
  for(i=0 ;i<*nk;i++){
    if(jqdtil[i+1]-jqdtil[i]==1){
      contribution=1.0/Ddtil[jqdtil[i]-1];
      nodesf=irowdtil[jqdtil[i]-1];
      nodem=i+1;
      insertas_ws(&irowdinv,&nodesf,
		  &nodem,&ifree2,&nzsdinv,&contribution,&Ddinv);
    }else{     
      //something went wrong
    }
    jqdinv[i+1]=ifree2;
  }
  nzsbtil=nzsb1;
  RENEW(Bdtil,double,nzsbtil);
  RENEW(irowbtil,ITG,nzsbtil);
  
  /* Bd on master side plus contributions from slave nodes without LM
     contribution (no-gap and no-LM) */
  
  add_rect(Bd1,irowb1,jqb1,dim,dim,
	   Bdtil2t,irowbtil2t,jqbtil2t,dim,dim,
	   &Bdtil,&irowbtil,jqbtil,&nzsbtil);
  
  /* Bdhelp=Ddtil^-1*Bdtil **/
  
  nzsbhelp=nzsbtil;
  RENEW(Bdhelp,double,nzsbhelp);
  RENEW(irowbhelp,ITG,nzsbhelp);
  dim=*nk;
   
  jqbhelp[0]=jqbtil[0];
  for(ii=0;ii<*nk;ii++){
    for(jj=jqbtil[ii]-1;jj<jqbtil[ii+1]-1;jj++){
      k=irowbtil[jj]-1;
      if(jqdinv[k+1]-jqdinv[k]==1){
	Bdhelp[jj]=Ddinv[jqdinv[k]-1]*Bdtil[jj];
	irowbhelp[jj]=irowbtil[jj];
      }else{
	// something went wrong
      }
    }
    jqbhelp[ii+1]=jqbtil[ii+1];
  }
  
  /* Ddtil2=Dd^-1*D^T= ID **/
  /* this is only help matrix needed for changed procedure dealing 
     with quadratic elements */
  /* OLD,not needed any more */
  
  nzsdtil2=nzsdtil;
  NNEW(Ddtil2,double,nzsdtil2);
  NNEW(irowdtil2,ITG,nzsdtil2);
  NNEW(jqdtil2,ITG,*nk+1);
  dim=*nk;
  jqdtil2[0]=jqdtil[0];
  for(ii=0;ii<*nk;ii++){
    for(jj=jqdtil[ii]-1;jj<jqdtil[ii+1]-1;jj++){
      Ddtil2[jj]=1.0;
      irowdtil2[jj]=irowdtil[jj];
    }
    jqdtil2[ii+1]=jqdtil[ii+1];
  }  
  
  /* create matrices in active dofs */
  /* create ddinv aubd aubdtil **/
  
  /* audd <- Dd **/
  
  ifree=1;
  nzsdd=3*nzsd;
  RENEW(audd,double,nzsdd);
  RENEW(irowdd,ITG,nzsdd);
  NNEW(mast1,ITG,nzsdd);
  for(i=0;i<*nk;i++){
    for(j=jqd[i]-1;j<jqd[i+1]-1;j++){
      nodesf=irowd[j];
      nodem=i+1;
      contribution=Dd[j];
      isn=islavnodeinv[irowd[j]-1];
      imn=islavnodeinv[i];
      k=0;
      if(ithermal[0]==3){k=-1;}
      for(ll=k; ll<3; ll++){	  				
	idofs=nactdof[mt*(nodesf-1)+ll+1];	    				
	idofm=nactdof[mt*(nodem-1)+ll+1];					
	if((idofs>0)&&(idofm>0)&&(contribution>1e-18 ||contribution<-1e-18)){
	  //insertion for active dofs  	      				
	  insertas(&irowdd,&mast1,&idofs,&idofm,&ifree,&nzsdd,
		   &contribution,&audd);				        
	}else if((idofs>0)&&(idofm<=0)&&(contribution>1e-18 ||
					 contribution<-1e-18)){
	  // mpc on master node			                  
	  for(jj=nslavmpc[2*(imn-1)];jj<nslavmpc[2*(imn-1)+1];jj++){
	    ist=islavmpc[jj];                                           
	    dirdep=nodempc[3*(ist-1)+1];
	    coefdep=coefmpc[ist-1];                                           
	    index=nodempc[3*(ist-1)+2];					   
	    if(ll==(dirdep-1)){                                           
	      while(index!=0){
		node1=nodempc[3*(index-1)];
		dirind=nodempc[3*(index-1)+1];
		c2=-coefmpc[index-1]*contribution/coefdep;
		idofm=nactdof[mt*(node1-1)+(dirind-1)+1];
		if(idofm>0&&(c2>1e-14 ||c2<-1e-14)){
		  insertas(&irowdd,&mast1,&idofs,&idofm,&ifree,&nzsdd,
			   &c2,&audd);					    
		}                                            
		index=nodempc[3*(index-1)+2];
	      }	                                   
	    }					  
	  }					
	}else if((idofm>0)&&(idofs<=0)&&(contribution>1e-18 ||
					 contribution<-1e-18)){
	  
	  //mpc on slave node
	  
	  for(jj=nslavmpc[2*(isn-1)];jj<nslavmpc[2*(isn-1)+1];jj++){
	    ist=islavmpc[jj];                                           
	    dirdep=nodempc[3*(ist-1)+1];
	    coefdep=coefmpc[ist-1];                                           
	    index=nodempc[3*(ist-1)+2];					   
	    if(ll==(dirdep-1)){                                           
	      while(index!=0){
		node1=nodempc[3*(index-1)];
		dirind=nodempc[3*(index-1)+1];
		c2=-coefmpc[index-1]*contribution/coefdep;
		idofs=nactdof[mt*(node1-1)+(dirind-1)+1];
		if(idofs>0&&(c2>1e-14 ||c2<-1e-14)){
		  insertas(&irowdd,&mast1,&idofs,&idofm,&ifree,&nzsdd,
			   &c2,&audd);
		  
		}
		index=nodempc[3*(index-1)+2];
	      }	                                   
	    }					  
	  }
	}else if((idofs<=0)&&(idofm<=0)&&(contribution>1e-18 ||
					  contribution<-1e-18)){
	  
	  //mpc on master and slave node
	  
	  for(jj=nslavmpc[2*(imn-1)];jj<nslavmpc[2*(imn-1)+1];jj++){
	    ist=islavmpc[jj];                                           
	    dirdep=nodempc[3*(ist-1)+1];
	    coefdep=coefmpc[ist-1];                                           
	    index=nodempc[3*(ist-1)+2];					   
	    if(ll==(dirdep-1)){                                          
	      while(index!=0){
		node1=nodempc[3*(index-1)];
		dirind=nodempc[3*(index-1)+1];
		c2=-coefmpc[index-1]*contribution/coefdep;
		idofm=nactdof[mt*(node1-1)+(dirind-1)+1];		        
		for(kk=nslavmpc[2*(isn-1)];kk<nslavmpc[2*(isn-1)+1];kk++){
		  ist2=islavmpc[kk];
		  dirdep2=nodempc[3*(ist2-1)+1];
		  coefdep2=coefmpc[ist2-1];
		  index2=nodempc[3*(ist2-1)+2];
		  if(ll==(dirdep2-1)){
		    while(index2!=0){
		      node2=nodempc[3*(index2-1)];
		      dirind2=nodempc[3*(index2-1)+1];
		      c3=-coefmpc[index2-1]*c2/coefdep2;   
		      idofs=nactdof[mt*(node2-1)+(dirind2-1)+1];
		      if(idofs>0&&idofm>0&&(c3>1e-14 ||c3<-1e-14)){
			insertas(&irowdd,&mast1,&idofs,&idofm,&ifree,&nzsdd,
				 &c3,&audd);
		      }                                              
		      index2=nodempc[3*(index2-1)+2];
		    }	                                     
		  }					    
		}
		index=nodempc[3*(index-1)+2];
	      }	                                   
	    }					  
	  }
	}  				
      }        		       
    }
  }  
  nzsdd=ifree-1;
  RENEW(irowdd,ITG,nzsdd);
  RENEW(audd,double,nzsdd);
  dim=neq[1];
  matrixsort(audd,mast1,irowdd,jqdd,&nzsdd,&dim);
  SFREE(mast1);
  icounter=0;
  for(i=0;i<neq[1];i++){
    if(jqdd[i]!=jqdd[i+1]){
      irowdd[icounter]=irowdd[jqdd[i]-1];
      audd[icounter]=audd[jqdd[i]-1];
      icounter++;
      istart=icounter;
      for(j=jqdd[i];j<jqdd[i+1]-1;j++){
	if(irowdd[j]==irowdd[icounter-1]){
	  audd[icounter-1]+=audd[j];   
	}else{
	  irowdd[icounter]=irowdd[j];
	  audd[icounter]=audd[j];
	  icounter++;
	}
      }
    }else{ istart=icounter+1;}
    
    jqdd[i]=istart;
  }
  jqdd[neq[1]]=icounter+1; 
  nzsdd=icounter;
  RENEW(irowdd,ITG,nzsdd);
  RENEW(audd,double,nzsdd); 
  
  /* aubd <-Bd **/
  
  ifree=1;
  *nzsbd=3*nzsb;
  RENEW(aubd,double,*nzsbd);
  RENEW(irowbd,ITG,*nzsbd);
  NNEW(mast1,ITG,*nzsbd);
  for(i=0;i<*nk;i++){
    for(j=jqb[i]-1;j<jqb[i+1]-1;j++){
      nodesf=irowb[j];
      nodem=i+1;
      isn=islavnodeinv[irowb[j]-1];
      imn=islavnodeinv[i];
      contribution=Bd[j];
      k=0;
      if(ithermal[0]==3){k=-1;}
      for(ll=k; ll<3; ll++){	  				
	idofs=nactdof[mt*(nodesf-1)+ll+1];	    				
	idofm=nactdof[mt*(nodem-1)+ll+1];					
	if((idofs>0)&&(idofm>0)&&(contribution>1e-18 ||contribution<-1e-18)){
	  
	  //insertion for active dofs
	  
	  insertas(&irowbd,&mast1,&idofs,&idofm,&ifree,nzsbd,
		   &contribution,&aubd);				        
	}else if((idofs>0)&&(contribution>1e-18 ||contribution<-1e-18)){
	  
	  // mpc on master node
	  
	  if(imn<0){
	    for(jj=nmastmpc[2*(-imn-1)];jj<nmastmpc[2*(-imn-1)+1];jj++){
	      ist=imastmpc[jj];                                           
	      dirdep=nodempc[3*(ist-1)+1];
	      coefdep=coefmpc[ist-1];                                           
	      index=nodempc[3*(ist-1)+2];
	      if(ll==(dirdep-1)){                                           
		while(index!=0){
		  node1=nodempc[3*(index-1)];
		  dirind=nodempc[3*(index-1)+1];
		  c2=-coefmpc[index-1]*contribution/coefdep;
		  idofm=nactdof[mt*(node1-1)+(dirind-1)+1];
		  if(idofm>0&&(c2>1e-14 ||c2<-1e-14)){
		    insertas(&irowbd,&mast1,&idofs,&idofm,&ifree,nzsbd,
			     &c2,&aubd);
		  }                                            
		  index=nodempc[3*(index-1)+2];
		}	                                   
	      }					  
	    }
	  }else{
	    for(jj=nslavmpc[2*(imn-1)];jj<nslavmpc[2*(imn-1)+1];jj++){
	      ist=islavmpc[jj];                                           
	      dirdep=nodempc[3*(ist-1)+1];
	      coefdep=coefmpc[ist-1];                                           
	      index=nodempc[3*(ist-1)+2];
	      if(ll==(dirdep-1)){                                           
		while(index!=0){
		  node1=nodempc[3*(index-1)];
		  dirind=nodempc[3*(index-1)+1];
		  c2=-coefmpc[index-1]*contribution/coefdep;
		  idofm=nactdof[mt*(node1-1)+(dirind-1)+1];
		  if(idofm>0&&(c2>1e-14 ||c2<-1e-14)){
		    insertas(&irowbd,&mast1,&idofs,&idofm,&ifree,nzsbd,
			     &c2,&aubd);
		  }                                            
		  index=nodempc[3*(index-1)+2];
		}	                                   
	      }					  
	    }	    
	  }
	}else if((idofm>0)&&(contribution>1e-18 ||contribution<-1e-18)){
	  
	  //mpc on slave node
	  
	  for(jj=nslavmpc[2*(isn-1)];jj<nslavmpc[2*(isn-1)+1];jj++){
	    ist=islavmpc[jj];                                           
	    dirdep=nodempc[3*(ist-1)+1];
	    coefdep=coefmpc[ist-1];                                           
	    index=nodempc[3*(ist-1)+2];					   
	    if(ll==(dirdep-1)){                                           
	      while(index!=0){
		node1=nodempc[3*(index-1)];
		dirind=nodempc[3*(index-1)+1];
		c2=-coefmpc[index-1]*contribution/coefdep;
		idofs=nactdof[mt*(node1-1)+(dirind-1)+1];
		if(idofs>0&&(c2>1e-14 ||c2<-1e-14)){
		  insertas(&irowbd,&mast1,&idofs,&idofm,&ifree,nzsbd,
			   &c2,&aubd);
		  
		}
		index=nodempc[3*(index-1)+2];
	      }	                                   
	    }					  
	  }
	}else if((idofs<=0)&&(idofm<=0)&&(contribution>1e-18 ||
					  contribution<-1e-18)){

	  //mpc on master and slave node
	  
	  if(imn<0){
	    for(jj=nmastmpc[2*(-imn-1)];jj<nmastmpc[2*(-imn-1)+1];jj++){
	      ist=imastmpc[jj];                                           
	      dirdep=nodempc[3*(ist-1)+1];
	      coefdep=coefmpc[ist-1];                                           
	      index=nodempc[3*(ist-1)+2];
	      if(ll==(dirdep-1)){                                          
		while(index!=0){
		  node1=nodempc[3*(index-1)];
		  dirind=nodempc[3*(index-1)+1];
		  c2=-coefmpc[index-1]*contribution/coefdep;
		  idofm=nactdof[mt*(node1-1)+(dirind-1)+1];		        
		  for(kk=nslavmpc[2*(isn-1)];kk<nslavmpc[2*(isn-1)+1];kk++){
		    ist2=islavmpc[kk];
		    dirdep2=nodempc[3*(ist2-1)+1];
		    coefdep2=coefmpc[ist2-1];
		    index2=nodempc[3*(ist2-1)+2];  
		    if(ll==(dirdep2-1)){
		      while(index2!=0){
			node2=nodempc[3*(index2-1)];
			dirind2=nodempc[3*(index2-1)+1];
			c3=-coefmpc[index2-1]*c2/coefdep2;   
			idofs=nactdof[mt*(node2-1)+(dirind2-1)+1];
			if(idofs>0&&idofm>0&&(c3>1e-14 ||c3<-1e-14)){
			  insertas(&irowbd,&mast1,&idofs,&idofm,&ifree,nzsbd,
				   &c3,&aubd);    
			}                                             
			index2=nodempc[3*(index2-1)+2];
		      }	                                     
		    } 
		  } 			        
		  index=nodempc[3*(index-1)+2];
		}	                                   
	      }					  
	    }
	  }else{
	    for(jj=nslavmpc[2*(imn-1)];jj<nslavmpc[2*(imn-1)+1];jj++){
	      ist=islavmpc[jj];                                           
	      dirdep=nodempc[3*(ist-1)+1];
	      coefdep=coefmpc[ist-1];                                           
	      index=nodempc[3*(ist-1)+2];
	      if(ll==(dirdep-1)){                                          
		while(index!=0){
		  node1=nodempc[3*(index-1)];
		  dirind=nodempc[3*(index-1)+1];
		  c2=-coefmpc[index-1]*contribution/coefdep;
		  idofm=nactdof[mt*(node1-1)+(dirind-1)+1];		        
		  for(kk=nslavmpc[2*(isn-1)];kk<nslavmpc[2*(isn-1)+1];kk++){
		    ist2=islavmpc[kk];
		    dirdep2=nodempc[3*(ist2-1)+1];
		    coefdep2=coefmpc[ist2-1];
		    index2=nodempc[3*(ist2-1)+2];  
		    if(ll==(dirdep2-1)){
		      while(index2!=0){
			node2=nodempc[3*(index2-1)];
			dirind2=nodempc[3*(index2-1)+1]; 
			c3=-coefmpc[index2-1]*c2/coefdep2;   
			idofs=nactdof[mt*(node2-1)+(dirind2-1)+1];
			if(idofs>0&&idofm>0&&(c3>1e-14 ||c3<-1e-14)){
			  insertas(&irowbd,&mast1,&idofs,&idofm,&ifree,nzsbd,
				   &c3,&aubd);    
			}                                              
			index2=nodempc[3*(index2-1)+2];
		      }	                                     
		    } 
		  } 			        
		  index=nodempc[3*(index-1)+2];
		}	                                   
	      }					  
	    }	    
	  }
	}  				
      }        		       
    }
  }  
  *nzsbd=ifree-1;
  RENEW(irowbd,ITG,*nzsbd);
  RENEW(aubd,double,*nzsbd);
  dim=neq[1];
  matrixsort(aubd,mast1,irowbd,jqbd,nzsbd,&dim);
  SFREE(mast1);
  icounter=0;
  for(i=0;i<neq[1];i++){
    if(jqbd[i]!=jqbd[i+1]){
      irowbd[icounter]=irowbd[jqbd[i]-1];
      aubd[icounter]=aubd[jqbd[i]-1];
      icounter++;
      istart=icounter;
      for(j=jqbd[i];j<jqbd[i+1]-1;j++){
	if(irowbd[j]==irowbd[icounter-1]){
	  aubd[icounter-1]+=aubd[j];   
	}else{
	  irowbd[icounter]=irowbd[j];
	  aubd[icounter]=aubd[j];
	  icounter++;
	}
      }
    }else{ istart=icounter+1;}
    
    jqbd[i]=istart;
  }
  jqbd[neq[1]]=icounter+1; 
  *nzsbd=icounter;
  RENEW(irowbd,ITG,*nzsbd);
  RENEW(aubd,double,*nzsbd);   
  
  /* auddtil2 <- Id **/
  
  ifree=1;
  nzsddtil2=3*nzsd;
  RENEW(auddtil2,double,nzsddtil2);
  RENEW(irowddtil2,ITG,nzsddtil2);
  NNEW(mast1,ITG,nzsddtil2);
  for(i=0;i<*nk;i++){
    for(j=jqdtil2[i]-1;j<jqdtil2[i+1]-1;j++){
      nodesf=irowdtil2[j];
      nodem=i+1;
      contribution=Ddtil2[j];
      isn=islavnodeinv[irowdtil2[j]-1];
      imn=islavnodeinv[i];
      k=0;
      if(ithermal[0]==3){k=-1;}
      for(ll=k; ll<3; ll++){	  				
	idofs=nactdof[mt*(nodesf-1)+ll+1];	    				
	idofm=nactdof[mt*(nodem-1)+ll+1];					
	if((idofs>0)&&(idofm>0)&&(contribution>1e-18 ||contribution<-1e-18)){
	  
	  //insertion for active dofs
	  
	  insertas(&irowddtil2,&mast1,&idofs,&idofm,&ifree,&nzsddtil2,
		   &contribution,&auddtil2);				        
	}
      }        		       
    }
  }  
  nzsddtil2=ifree-1;
  RENEW(irowddtil2,ITG,nzsddtil2);
  RENEW(auddtil2,double,nzsddtil2);
  dim=neq[1];
  matrixsort(auddtil2,mast1,irowddtil2,jqddtil2,&nzsddtil2,&dim);
  SFREE(mast1); 
  icounter=0;
  for(i=0;i<neq[1];i++){
    if(jqddtil2[i]!=jqddtil2[i+1]){
      irowddtil2[icounter]=irowddtil2[jqddtil2[i]-1];
      auddtil2[icounter]=auddtil2[jqddtil2[i]-1];
      icounter++;
      istart=icounter;
      for(j=jqddtil2[i];j<jqddtil2[i+1]-1;j++){
	if(irowddtil2[j]==irowddtil2[icounter-1]){
	  auddtil2[icounter-1]+=auddtil2[j];   
	}else{
	  irowddtil2[icounter]=irowddtil2[j];
	  auddtil2[icounter]=auddtil2[j];
	  icounter++;
	}
      }
    }else{ istart=icounter+1;}
    
    jqddtil2[i]=istart;
  }
  jqddtil2[neq[1]]=icounter+1; 
  nzsddtil2=icounter;
  RENEW(irowddtil2,ITG,nzsddtil2);
  RENEW(auddtil2,double,nzsddtil2); 

  /* aubdtil2 <- Ddtil,Bdtil in one matrix due to MPCs **/
  
  ifree=1;
  *nzsbdtil2=3*nzsbhelp;
  RENEW(aubdtil2,double,*nzsbdtil2);
  RENEW(irowbdtil2,ITG,*nzsbdtil2);
  NNEW(mast1,ITG,*nzsbdtil2);
  
  for(i=0;i<*nk;i++){
    for(j=jqdtil[i]-1;j<jqdtil[i+1]-1;j++){
      nodesf=irowdtil[j];
      nodem=i+1;
      contribution=Ddtil[j];
      isn=islavnodeinv[irowdtil[j]-1];
      imn=islavnodeinv[i];
      k=0;
      if(ithermal[0]==3){k=-1;}
      for(ll=k; ll<3; ll++){	  				
	idofs=nactdof[mt*(nodesf-1)+ll+1];	    				
	idofm=nactdof[mt*(nodem-1)+ll+1];					
	if((idofs>0)&&(idofm>0)&&(contribution>1e-18 ||contribution<-1e-18)){
	  
	  //insertion for active dofs
	  
	  insertas(&irowbdtil2,&mast1,&idofs,&idofm,&ifree,nzsbdtil2,
		   &contribution,&aubdtil2);				        
	}else if((idofs>0)&&(contribution>1e-18 ||contribution<-1e-18)){
	  
	  // mpc on master node
	  
	  for(jj=nslavmpc[2*(imn-1)];jj<nslavmpc[2*(imn-1)+1];jj++){
	    ist=islavmpc[jj];
	    dirdep=nodempc[3*(ist-1)+1];
	    coefdep=coefmpc[ist-1];
	    index=nodempc[3*(ist-1)+2];
	    if(ll==(dirdep-1)){
	      while(index!=0){
		node1=nodempc[3*(index-1)];  
		dirind=nodempc[3*(index-1)+1];
		c2=-coefmpc[index-1]*contribution/coefdep;
		if(node1>*nk){
		  idofm=0;
		}else{
		  idofm=nactdof[mt*(node1-1)+(dirind-1)+1];
		}
		if(idofm>0&&(c2>1e-14 ||c2<-1e-14)){	
		  insertas(&irowbdtil2,&mast1,&idofs,&idofm,&ifree,nzsbdtil2,
			   &c2,&aubdtil2); 
		} 
		index=nodempc[3*(index-1)+2];
	      }	                                   
	    }					  
	  }					
	}else if((idofm>0)&&(contribution>1e-18 ||contribution<-1e-18)){
	  //mpc on slave node			                  
	  for(jj=nslavmpc[2*(isn-1)];jj<nslavmpc[2*(isn-1)+1];jj++){
	    ist=islavmpc[jj];
	    dirdep=nodempc[3*(ist-1)+1];
	    coefdep=coefmpc[ist-1];
	    index=nodempc[3*(ist-1)+2];
	    if(ll==(dirdep-1)){
	      while(index!=0){
		node1=nodempc[3*(index-1)];  
		dirind=nodempc[3*(index-1)+1];
		c2=-coefmpc[index-1]*contribution/coefdep;
		if(node1>*nk){
		  idofs=0;
		}else{
		  idofs=nactdof[mt*(node1-1)+(dirind-1)+1];
		}
		if(idofs>0&&(c2>1e-14 ||c2<-1e-14)){
		  insertas(&irowbdtil2,&mast1,&idofs,&idofm,&ifree,nzsbdtil2,
			   &c2,&aubdtil2);
		  
		}
		index=nodempc[3*(index-1)+2];
	      }	                                   
	    }					  
	  }										
	}else if((idofs<=0)&&(idofm<=0)&&(contribution>1e-18 ||contribution<-1e-18)){
	  //mpc on master and slave node			                  
	  for(jj=nslavmpc[2*(imn-1)];jj<nslavmpc[2*(imn-1)+1];jj++){
	    ist=islavmpc[jj];
	    dirdep=nodempc[3*(ist-1)+1];
	    coefdep=coefmpc[ist-1];
	    index=nodempc[3*(ist-1)+2];
	    if(ll==(dirdep-1)){                                          
	      while(index!=0){
		node1=nodempc[3*(index-1)]; 
		dirind=nodempc[3*(index-1)+1];
		c2=-coefmpc[index-1]*contribution/coefdep;
		if(node1>*nk){
		  idofm=0;
		}else{
		  idofm=nactdof[mt*(node1-1)+(dirind-1)+1];
		}
		for(kk=nslavmpc[2*(isn-1)];kk<nslavmpc[2*(isn-1)+1];kk++){  
		  ist2=islavmpc[kk];  
		  dirdep2=nodempc[3*(ist2-1)+1];  
		  coefdep2=coefmpc[ist2-1];  
		  index2=nodempc[3*(ist2-1)+2];  
		  if(ll==(dirdep2-1)){   
		    while(index2!=0){   
		      node2=nodempc[3*(index2-1)];    
		      dirind2=nodempc[3*(index2-1)+1];   
		      c3=-coefmpc[index2-1]*c2/coefdep2;
		      if(node2>*nk){
			idofs=0;
		      }else{
		        idofs=nactdof[mt*(node2-1)+(dirind2-1)+1]; 
		      }
		      if(idofs>0&&idofm>0&&(c3>1e-14 ||c3<-1e-14)){
			insertas(&irowbdtil2,&mast1,&idofs,&idofm,&ifree,
				 nzsbdtil2,&c3,&aubdtil2);    
		      }   
		      index2=nodempc[3*(index2-1)+2];  
		    } 
		  } 
		} 				   
		index=nodempc[3*(index-1)+2];
	      }	                                   
	    }					  
	  }
	}  				
      }        		       
    }
  }  
  
  for(i=0;i<*nk;i++){
    for(j=jqbtil[i]-1;j<jqbtil[i+1]-1;j++){
      nodesf=irowbtil[j];
      nodem=i+1;
      isn=islavnodeinv[irowbtil[j]-1];
      imn=islavnodeinv[i];
      contribution=Bdtil[j];
      k=0;
      if(ithermal[0]==3){k=-1;}
      for(ll=k; ll<3; ll++){	  				
	idofs=nactdof[mt*(nodesf-1)+ll+1];	    				
	idofm=nactdof[mt*(nodem-1)+ll+1];					
	if((idofs>0)&&(idofm>0)){
	  
	  //insertion for active dofs
	  
	  insertas(&irowbdtil2,&mast1,&idofs,&idofm,&ifree,nzsbdtil2,
		   &contribution,&aubdtil2);				        
	}else if((idofs>0)&&(contribution>1e-18 ||contribution<-1e-18)){
	  
	  // mpc on master node
	  
	  if(imn<0){

	    /* mpc2 to be used, since slave middle nodes may occur on
               the ind side */
	      
	    for(jj=nmastmpc[2*(-imn-1)];jj<nmastmpc[2*(-imn-1)+1];jj++){
	      ist=imastmpc[jj];
	      dirdep=nodempc[3*(ist-1)+1];
	      coefdep=coefmpc[ist-1];
	      index=nodempc[3*(ist-1)+2];
	      if(ll==(dirdep-1)){
		while(index!=0){
		  node1=nodempc[3*(index-1)];  
		  dirind=nodempc[3*(index-1)+1];
		  c2=-coefmpc[index-1]*contribution/coefdep;
		  if(node1>*nk){
		    idofm=0;
		  }else{
		    idofm=nactdof[mt*(node1-1)+(dirind-1)+1];
		  }
		  if(idofm>0){
		    insertas(&irowbdtil2,&mast1,&idofs,&idofm,&ifree,nzsbdtil2,
			     &c2,&aubdtil2); 
		  } 
		  index=nodempc[3*(index-1)+2];
		}	                                   
	      }					  
	    }
	  }else{
	    for(jj=nslavmpc[2*(imn-1)];jj<nslavmpc[2*(imn-1)+1];jj++){
	      ist=islavmpc[jj];
	      dirdep=nodempc[3*(ist-1)+1];
	      coefdep=coefmpc[ist-1];
	      index=nodempc[3*(ist-1)+2];
	      if(ll==(dirdep-1)){
		while(index!=0){
		  node1=nodempc[3*(index-1)];  
		  dirind=nodempc[3*(index-1)+1];
		  c2=-coefmpc[index-1]*contribution/coefdep;
		  if(node1>*nk){
		    idofm=0;
		  }else{
		    idofm=nactdof[mt*(node1-1)+(dirind-1)+1];
		  }
		  if(idofm>0){
		    insertas(&irowbdtil2,&mast1,&idofs,&idofm,&ifree,nzsbdtil2,
			     &c2,&aubdtil2); 
		  } 
		  index=nodempc[3*(index-1)+2];
		}	                                   
	      }					  
	    }	    
	  }
	}else if((idofm>0)&&(contribution>1e-18 ||contribution<-1e-18)){
	  
	  //mpc on slave node
	  
	  for(jj=nslavmpc[2*(isn-1)];jj<nslavmpc[2*(isn-1)+1];jj++){
	    ist=islavmpc[jj];
	    dirdep=nodempc[3*(ist-1)+1];
	    coefdep=coefmpc[ist-1];
	    index=nodempc[3*(ist-1)+2];
	    if(ll==(dirdep-1)){
	      while(index!=0){
		node1=nodempc[3*(index-1)];  
		dirind=nodempc[3*(index-1)+1];
		c2=-coefmpc[index-1]*contribution/coefdep;
		if(node1>*nk){
		  idofs=0;
		}else{		
		  idofs=nactdof[mt*(node1-1)+(dirind-1)+1];
		}
		if(idofs>0){
		  insertas(&irowbdtil2,&mast1,&idofs,&idofm,&ifree,nzsbdtil2,
			   &c2,&aubdtil2);
		  
		}
		index=nodempc[3*(index-1)+2];
	      }	                                   
	    }					  
	  }
	}else if((idofs<=0)&&(idofm<=0)&&(contribution>1e-18 ||
					  contribution<-1e-18)){
	  
	  //mpc on master and slave node
	  
	  if(imn<0){
	    for(jj=nmastmpc[2*((-imn)-1)];jj<nmastmpc[2*((-imn)-1)+1];jj++){
	      ist=imastmpc[jj];
	      dirdep=nodempc[3*(ist-1)+1];
	      coefdep=coefmpc[ist-1];
	      index=nodempc[3*(ist-1)+2];
	      if(ll==(dirdep-1)){                                          
		while(index!=0){
		  node1=nodempc[3*(index-1)]; 
		  dirind=nodempc[3*(index-1)+1];
		  c2=-coefmpc[index-1]*contribution/coefdep;
		  if(node1>*nk){
		    idofm=0;
		  }else{
		    idofm=nactdof[mt*(node1-1)+(dirind-1)+1];
		  }		  
		  for(kk=nslavmpc[2*(isn-1)];kk<nslavmpc[2*(isn-1)+1];kk++){  
		    ist2=islavmpc[kk];  
		    dirdep2=nodempc[3*(ist2-1)+1];  
		    coefdep2=coefmpc[ist2-1];  
		    index2=nodempc[3*(ist2-1)+2];  
		    if(ll==(dirdep2-1)){   
		      while(index2!=0){   
			node2=nodempc[3*(index2-1)];    
			dirind2=nodempc[3*(index2-1)+1];   
			c3=-coefmpc[index2-1]*c2/coefdep2;
			if(node2>*nk){
			  idofs=0;
			}else{
			  idofs=nactdof[mt*(node2-1)+(dirind2-1)+1];
			}
			if(idofs>0&&idofm>0){
			  insertas(&irowbdtil2,&mast1,&idofs,&idofm,&ifree,
				   nzsbdtil2,&c3,&aubdtil2);    
			}   
			index2=nodempc[3*(index2-1)+2];  
		      } 
		    } 
		  } 				   
		  index=nodempc[3*(index-1)+2];
		}	                                   
	      }					  
	    }
	  }else{
	    for(jj=nslavmpc[2*(imn-1)];jj<nslavmpc[2*(imn-1)+1];jj++){
	      ist=islavmpc[jj];
	      dirdep=nodempc[3*(ist-1)+1];
	      coefdep=coefmpc[ist-1];
	      index=nodempc[3*(ist-1)+2];
	      if(ll==(dirdep-1)){                                          
		while(index!=0){
		  node1=nodempc[3*(index-1)]; 
		  dirind=nodempc[3*(index-1)+1];
		  c2=-coefmpc[index-1]*contribution/coefdep;
		  if(node1>*nk){
		    idofm=0;
		  }else{
		    idofm=nactdof[mt*(node1-1)+(dirind-1)+1];
		  }
		  for(kk=nslavmpc[2*(isn-1)];kk<nslavmpc[2*(isn-1)+1];kk++){  
		    ist2=islavmpc[kk];  
		    dirdep2=nodempc[3*(ist2-1)+1];  
		    coefdep2=coefmpc[ist2-1];  
		    index2=nodempc[3*(ist2-1)+2];  
		    if(ll==(dirdep2-1)){   
		      while(index2!=0){   
			node2=nodempc[3*(index2-1)];    
			dirind2=nodempc[3*(index2-1)+1];   
			c3=-coefmpc[index2-1]*c2/coefdep2;
			if(node2>*nk){
			  idofs=0;
			}else{
			  idofs=nactdof[mt*(node2-1)+(dirind2-1)+1];
			}
			if(idofs>0&&idofm>0){
			  insertas(&irowbdtil2,&mast1,&idofs,&idofm,&ifree,
				   nzsbdtil2,&c3,&aubdtil2);    
			}   
			index2=nodempc[3*(index2-1)+2];  
		      } 
		    } 
		  } 				   
		  index=nodempc[3*(index-1)+2];
		}	                                   
	      }					  
	    }	    
	  }
	}  				
      }        		       
    }
  }  
  *nzsbdtil2=ifree-1;
  RENEW(irowbdtil2,ITG,*nzsbdtil2);
  RENEW(aubdtil2,double,*nzsbdtil2);
  dim=neq[1];
  matrixsort(aubdtil2,mast1,irowbdtil2,jqbdtil2,nzsbdtil2,&dim);
  SFREE(mast1);
  icounter=0;
  for(i=0;i<neq[1];i++){
    if(jqbdtil2[i]!=jqbdtil2[i+1]){
      irowbdtil2[icounter]=irowbdtil2[jqbdtil2[i]-1];
      aubdtil2[icounter]=aubdtil2[jqbdtil2[i]-1];
      icounter++;
      istart=icounter;
      for(j=jqbdtil2[i];j<jqbdtil2[i+1]-1;j++){
	if(irowbdtil2[j]==irowbdtil2[icounter-1]){
	  aubdtil2[icounter-1]+=aubdtil2[j];   
	}else{
	  irowbdtil2[icounter]=irowbdtil2[j];
	  aubdtil2[icounter]=aubdtil2[j];
	  icounter++;
	}
      }
    }else{ istart=icounter+1;}
    
    jqbdtil2[i]=istart;
  }
  jqbdtil2[neq[1]]=icounter+1; 
  *nzsbdtil2=icounter;
  RENEW(irowbdtil2,ITG,*nzsbdtil2);
  RENEW(aubdtil2,double,*nzsbdtil2);
  
  /* generate auddinv out of aubdtil2  inverting the 3x3 matricen **/
  
  ifree=1;
  nzsddinv=3*nzsd;
  RENEW(auddinv,double,nzsddinv);
  RENEW(irowddinv,ITG,nzsddinv);
  NNEW(mast1,ITG,nzsddinv);

  /* loop over all ties (ok, since there is no active slave nodes
     in 2 ties at the same time */
  
  for( i=0; i<*ntie; i++){    
    if(tieset[i*(81*3)+80]=='C'){ 
      for(j=nslavnode[i]; j<nslavnode[i+1]; j++){
	nodesf=islavnode[j];
	
	// get D_d 3x3
	
	if(islavact[j]>-1){
	  for(ll=0; ll<3; ll++){	  				
	    idofm=nactdof[mt*(nodesf-1)+ll+1]-1;
	    dloc[9*j+ll*3+0]=0.0;dloc[9*j+ll*3+1]=0.0;dloc[9*j+ll*3+2]=0.0;
	    dloc[9*j+ll*3+ll]=1.0;
	    // fix if ndof<3
	    if(idofm>-1){
	      for(k=jqbdtil2[idofm]-1;k<jqbdtil2[idofm+1]-1;k++){
		if(irowbdtil2[k]==nactdof[mt*(nodesf-1)+0+1]){
		  //dof1
		  dloc[9*j+ll*3+0]=aubdtil2[k];
		}else if(irowbdtil2[k]==nactdof[mt*(nodesf-1)+1+1]){
		  //dof2
		  dloc[9*j+ll*3+1]=aubdtil2[k];
		}else if(irowbdtil2[k]==nactdof[mt*(nodesf-1)+2+1]){
		  //dof3
		  dloc[9*j+ll*3+2]=aubdtil2[k];
		}
	      }
	    }
	  }
	  
	  // invert D_d 3x3
	  
	  detdloc=dloc[9*j+0]*dloc[9*j+4]*dloc[9*j+8]
	    +dloc[9*j+3]*dloc[9*j+7]*dloc[9*j+2]
	    +dloc[9*j+6]*dloc[9*j+1]*dloc[9*j+5]
	    -dloc[9*j+6]*dloc[9*j+4]*dloc[9*j+2]
	    -dloc[9*j+0]*dloc[9*j+7]*dloc[9*j+5]
	    -dloc[9*j+3]*dloc[9*j+1]*dloc[9*j+8];
	  if(detdloc<1.e-19&&detdloc>-1.e-19){
	    fflush(stdout); 
	    FORTRAN(stop,());
	  }
	  dinvloc[9*j+0]=+1.0/detdloc*(dloc[9*j+4]*dloc[9*j+8]-
				       dloc[9*j+7]*dloc[9*j+5]);
	  dinvloc[9*j+1]=-1.0/detdloc*(dloc[9*j+1]*dloc[9*j+8]-
				       dloc[9*j+7]*dloc[9*j+2]);
	  dinvloc[9*j+2]=+1.0/detdloc*(dloc[9*j+1]*dloc[9*j+5]-
				       dloc[9*j+4]*dloc[9*j+2]);
	  dinvloc[9*j+3]=-1.0/detdloc*(dloc[9*j+3]*dloc[9*j+8]-
				       dloc[9*j+6]*dloc[9*j+5]);
	  dinvloc[9*j+4]=+1.0/detdloc*(dloc[9*j+0]*dloc[9*j+8]-
				       dloc[9*j+6]*dloc[9*j+2]);
	  dinvloc[9*j+5]=-1.0/detdloc*(dloc[9*j+0]*dloc[9*j+5]-
				       dloc[9*j+3]*dloc[9*j+2]);
	  dinvloc[9*j+6]=+1.0/detdloc*(dloc[9*j+3]*dloc[9*j+7]-
				       dloc[9*j+6]*dloc[9*j+4]);
	  dinvloc[9*j+7]=-1.0/detdloc*(dloc[9*j+0]*dloc[9*j+7]-
				       dloc[9*j+6]*dloc[9*j+1]);
	  dinvloc[9*j+8]=+1.0/detdloc*(dloc[9*j+0]*dloc[9*j+4]-
				       dloc[9*j+3]*dloc[9*j+1]);
	  
	  //add to auddinv
	  
	  for(ll=0; ll<3; ll++){
	    idofm=nactdof[mt*(nodesf-1)+ll+1];
	    for(k=0;k<3;k++){
	      idofs=nactdof[mt*(nodesf-1)+k+1];
	      contribution=dinvloc[9*j+3*ll+k];
	      if(idofs>0 && idofm>0 && (contribution>1.e-18 ||
					contribution<-1.e-18)){	
		insertas(&irowddinv,&mast1,&idofs,&idofm,&ifree,&nzsddinv,
			 &contribution,&auddinv);		
	      }
	    }
	  }
	  if(ithermal[0]==3){
	    idofs=nactdof[mt*(nodesf)-4];
	    if(jqdinv[nodesf]-jqdinv[nodesf-1]==1){
	      contribution=Ddinv[jqdinv[nodesf-1]-1];
	    }else{
	      contribution=0;
	      // something went wrong!!!
	    }
	    if(idofs>0  && (contribution>1.e-18 || contribution<-1.e-18)){	
	      insertas(&irowddinv,&mast1,&idofs,&idofs,&ifree,&nzsddinv,
		       &contribution,&auddinv);		
	    }	    
	  }
	}
      }
    }
  }

  nzsddinv=ifree-1;
  RENEW(irowddinv,ITG,nzsddinv);
  RENEW(auddinv,double,nzsddinv);
  dim=neq[1];
  matrixsort(auddinv,mast1,irowddinv,jqddinv,&nzsddinv,&dim);
  SFREE(mast1); 
  icounter=0;
  for(i=0;i<neq[1];i++){
    if(jqddinv[i]!=jqddinv[i+1]){
      irowddinv[icounter]=irowddinv[jqddinv[i]-1];
      auddinv[icounter]=auddinv[jqddinv[i]-1];
      icounter++;
      istart=icounter;
      for(j=jqddinv[i];j<jqddinv[i+1]-1;j++){
	if(irowddinv[j]==irowddinv[icounter-1]){
	  auddinv[icounter-1]+=auddinv[j];   
	}else{
	  irowddinv[icounter]=irowddinv[j];
	  auddinv[icounter]=auddinv[j];
	  icounter++;
	}
      }
    }else{ istart=icounter+1;}
    
    jqddinv[i]=istart;
  }
  jqddinv[neq[1]]=icounter+1; 
  nzsddinv=icounter;
  RENEW(irowddinv,ITG,nzsddinv);
  RENEW(auddinv,double,nzsddinv); 
  
  /* aubdtil= Dtil^-1*Bdtil for active dofs **/
  /* gehe ueber aubdtil2 **/
  
  ifree=1;
  *nzsbdtil=*nzsbdtil2;
  RENEW(aubdtil,double,*nzsbdtil);
  RENEW(irowbdtil,ITG,*nzsbdtil);
  jqbdtil[0]=1;
  for(i=0;i<neq[0];i++){
    for(j=jqbdtil2[i]-1;j<jqbdtil2[i+1]-1;j++){
      idofs=irowbdtil2[j];
      idofm=i+1;
      if(islavactdof[idofm-1]<0 ||
	 (islavactdof[idofm-1]>0 &&
	  islavact[(ITG)floor(islavactdof[idofm-1]/10.)-1]<0) ){
	
	// only Mast Dofs and possibly N dofs due to mpcs
	
	if(islavactdof[idofs-1]<1 && islavactdof[idofm-1]>0 ){
	  printf("bdfill: something went wrong in aubdtil2! Stop!\n");
	  printf("bdfill: noder %" ITGFORMAT " %" ITGFORMAT " nodec %"
		 ITGFORMAT " %" ITGFORMAT " \n",
		 islavnode[(ITG) floor(islavactdof[idofs-1]/10.)-1],idofs,
		 islavnode[(ITG) floor(islavactdof[idofm-1]/10.)-1],idofm );
	  FORTRAN(stop,());
	}else if(islavactdof[idofs-1]<1 && islavactdof[idofm-1]<0 ){
	  printf("bdfill: something went wrong in aubdtil2! Stop!\n");
	  printf("bdfill: slavenode row   %" ITGFORMAT " %" ITGFORMAT
		 " master node column %" ITGFORMAT " %" ITGFORMAT " \n",
		 islavnode[(ITG) floor(islavactdof[idofs-1]/10.)-1],idofs,
		 imastnode[-((ITG) floor(islavactdof[idofm-1]/10.))-1],idofm );
	  FORTRAN(stop,()); 
	}else if(islavactdof[idofs-1]<1 && islavactdof[idofm-1]==0 ){
	  printf("bdfill: something went wrong in aubdtil2! Stop!\n");	
	  printf("bdfill: slavenode row   %" ITGFORMAT " %" ITGFORMAT
		 " master node column %" ITGFORMAT " \n",
		 islavnode[(ITG) floor(islavactdof[idofs-1]/10.)-1],idofs,
		 idofm );
	  FORTRAN(stop,());
	}
	islavnodeentry=floor(islavactdof[idofs-1]/10.);
	jrow= islavactdof[idofs-1]-10*islavnodeentry;
	nodesf=islavnode[islavnodeentry-1];	   
        idof1=nactdof[mt*nodesf-3];	   
        idof2=nactdof[mt*nodesf-2];	   
        idof3=nactdof[mt*nodesf-1];
	iadd=0;
        if(jrow==1){
	  // k=idof1
	  e1=aubdtil2[j];
	  if(j+1<jqbdtil2[i+1]-1 && irowbdtil2[j+1]==idof2){
	    e2=aubdtil2[j+1];++iadd;}
	  else{
	    e2=0.0;}
	  if(j+1<jqbdtil2[i+1]-1 && irowbdtil2[j+1]==idof3){
	    e3=aubdtil2[j+1];
	    ++iadd;}
	  else if(j+2<jqbdtil2[i+1]-1 && irowbdtil2[j+2]==idof3){
	    e3=aubdtil2[j+2];
	    ++iadd;}
	  else{e3=0.0;}
	  
	}else if(jrow==2){
	  // k=idof2
	  e1=0.0;
	  e2=aubdtil2[j];
	  if(j+1<jqbdtil2[i+1]-1 && irowbdtil2[j+1]==idof3){
	    e3=aubdtil2[j+1];
	    ++iadd;}
	  else{e3=0.0;}
	  
	}else{
	  // k=idof3
	  e1=0.0;
	  e2=0.0;
	  e3=aubdtil2[j];
	}
	contribution=dinvloc[9*(islavnodeentry-1)+0]*e1
	  +dinvloc[9*(islavnodeentry-1)+3]*e2
	  +dinvloc[9*(islavnodeentry-1)+6]*e3;
	if(idof1>0 && (contribution>1.e-18 || contribution<-1.e-18)){
          insertas_ws(&irowbdtil,&idof1,&idofm,&ifree,
		      nzsbdtil,&contribution,&aubdtil);
	}
	contribution=dinvloc[9*(islavnodeentry-1)+1]*e1
	  +dinvloc[9*(islavnodeentry-1)+4]*e2
	  +dinvloc[9*(islavnodeentry-1)+7]*e3;
	if(idof2>0 &&(contribution>1.e-18 || contribution<-1.e-18)){
          insertas_ws(&irowbdtil,&idof2,&idofm,&ifree,
		      nzsbdtil,&contribution,&aubdtil);
	}
	contribution=dinvloc[9*(islavnodeentry-1)+2]*e1
	  +dinvloc[9*(islavnodeentry-1)+5]*e2
	  +dinvloc[9*(islavnodeentry-1)+8]*e3;
	if(idof3>0 && (contribution>1.e-18 || contribution<-1.e-18)){
          insertas_ws(&irowbdtil,&idof3,&idofm,&ifree,
		      nzsbdtil,&contribution,&aubdtil);
	}
	j=j+iadd;
      }
      
    }
    jqbdtil[i+1]=ifree;
  }
  if(ithermal[0]==3){
    for(i=neq[0];i<neq[1];i++){
      for(j=jqbdtil2[i]-1;j<jqbdtil2[i+1]-1;j++){
	idofs=irowbdtil2[j];
	idofm=i+1;
	if(islavactdof[idofm-1]<0 ||
	   (islavactdof[idofm-1]>0 &&
	    islavact[(ITG) floor(islavactdof[idofm-1]/10.)-1]<0) ||
	   islavactdof[idofm-1]==0){
	  
	  // only Mast Dofs and possibly N dofs due to mpcs
	  
	  if(islavactdof[idofs-1]<1 && islavactdof[idofm-1]>0 ){
	    
	    printf("bdfill: something went wrong in aubdtil2! Stop!\n");
	    printf("bdfill: noder %" ITGFORMAT " %" ITGFORMAT " nodec %"
		   ITGFORMAT " %" ITGFORMAT " \n",
		   islavnode[(ITG) floor(islavactdof[idofs-1]/10.)-1],idofs,
		   islavnode[(ITG) floor(islavactdof[idofm-1]/10.)-1],idofm );
	    FORTRAN(stop,());
	    
	  }else if(islavactdof[idofs-1]<1 && islavactdof[idofm-1]<0 ){
	    
	    printf("bdfill: something went wrong in aubdtil2! Stop!\n");
	    printf("bdfill: slavenode row   %" ITGFORMAT " %" ITGFORMAT
		   " master node column %" ITGFORMAT " %" ITGFORMAT " \n",
		   islavnode[(ITG) floor(islavactdof[idofs-1]/10.)-1],idofs,
		   imastnode[-((ITG) floor(islavactdof[idofm-1]/10.))-1],idofm);
	    FORTRAN(stop,());
	    
	  }else if(islavactdof[idofs-1]<1 && islavactdof[idofm-1]==0 ){
	    
	    printf("bdfill: something went wrong in aubdtil2! Stop!\n");	
	    printf("bdfill: slavenode row   %" ITGFORMAT " %" ITGFORMAT
		   " master node column %" ITGFORMAT " \n",
		   islavnode[(ITG) floor(islavactdof[idofs-1]/10.)-1],idofs,
		   idofm );
	    FORTRAN(stop,());
	  }
	  islavnodeentry=floor(islavactdof[idofs-1]/10.);
	  jrow= islavactdof[idofs-1]-10*islavnodeentry;
	  nodesf=islavnode[islavnodeentry-1];	   
	  idof1=nactdof[mt*nodesf-3];	   
	  idof2=nactdof[mt*nodesf-2];	   
	  idof3=nactdof[mt*nodesf-1];
	  idof0=nactdof[mt*nodesf-4];
	  if(jqdinv[nodesf]-jqdinv[nodesf-1]==1){
            contribution=Ddinv[jqdinv[nodesf-1]-1]*aubdtil2[j];
	  }else{
	    contribution=0;
	    // something went wrong!!!
	  }
	  if(idof0>0 && (contribution>1.e-18 || contribution<-1.e-18)){
	    insertas_ws(&irowbdtil,&idof0,&idofm,&ifree,
			nzsbdtil,&contribution,&aubdtil);
	  }
	}
      
      }
      jqbdtil[i+1]=ifree;
    } 
  }
  
  *nzsbdtil=ifree-1;
  RENEW(irowbdtil,ITG,*nzsbdtil);
  RENEW(aubdtil,double,*nzsbdtil); 
  
  *irowbdp=irowbd; *aubdp=aubd; 
  *irowddp=irowdd; *auddp=audd; 
  *irowddtilp=irowddtil; *auddtilp=auddtil;
  *irowddtil2p=irowddtil2; *auddtil2p=auddtil2;
  *irowddinvp=irowddinv; *auddinvp=auddinv; 
  *irowddp=irowdd; *auddp=audd;
  *irowbdtilp=irowbdtil; *aubdtilp=aubdtil;
  *irowbdtil2p=irowbdtil2; *aubdtil2p=aubdtil2;  
  *irowbp=irowb; *Bdp=Bd;
  *irowbhelpp=irowbhelp; *Bdhelpp=Bdhelp;
  *irowdp=irowd; *Ddp=Dd;
  *irowdtilp=irowdtil; *Ddtilp=Ddtil;
  *irowbtilp=irowbtil; *Bdtilp=Bdtil;
  SFREE(Dd1);SFREE(irowd1);SFREE(jqd1);
  SFREE(Ddtil2);SFREE(irowdtil2);SFREE(jqdtil2);
  SFREE(Ddinv);SFREE(irowdinv);SFREE(jqdinv);
  SFREE(Bd1);SFREE(irowb1);SFREE(jqb1);
  SFREE(Bd2);SFREE(irowb2);SFREE(jqb2);
  SFREE(Bd2a);SFREE(irowb2a);SFREE(jqb2a);
  SFREE(Bdtil2);SFREE(irowbtil2);SFREE(jqbtil2);
  SFREE(Bdtil2t);SFREE(irowbtil2t);SFREE(jqbtil2t);
  SFREE(igap);
  SFREE(dinvloc);SFREE(dloc);

  return;
}
