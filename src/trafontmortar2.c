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

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

/*
  - Condense Lagrange Multiplier and embed contact conditions for K_{AX}
  - changing au due to N and T (normal and tangential
    direction at the slave surface) 
  - changing b due to N and T (normal and tangential
    direction at the slave surface) 
  phd-thesis Sitzmann, Chapter 3+4, equation (4.15) (Type=MORTAR)
  
  author: Saskia Sitzmann
 */

void trafontmortar2(ITG *neq,ITG *nzs,ITG *islavactdof,ITG *islavact,
		    ITG *nslavnode,double *f_da,double *f_atil,
		    double *au_dan,ITG *irow_dan,ITG *jq_dan,
		    double *au_dam,ITG *irow_dam,ITG *jq_dam,
		    double *au_dai,ITG *irow_dai,ITG *jq_dai,
		    double *au_daa,ITG *irow_daa,ITG *jq_daa,
		    double **au_antilp,ITG **irow_antilp,ITG *jq_antil,
		    double **au_amtilp,ITG **irow_amtilp,ITG *jq_amtil,
		    double **au_aitilp,ITG **irow_aitilp,ITG *jq_aitil,
		    double **au_aatilp,ITG **irow_aatilp,ITG *jq_aatil,
		    double *gap,
		    double *Bd,ITG *irowb,ITG *jqb,
		    double *Dd,ITG *irowd,ITG *jqd,
		    double *Ddtil,ITG *irowdtil,ITG *jqdtil,
		    double *au_bdtil2,ITG *irow_bdtil2,ITG *jq_bdtil2,
		    double *au_ddtil2i,ITG *irow_ddtil2i,ITG *jq_ddtil2i,
		    double *au_ddtil2a,ITG *irow_ddtil2a,ITG *jq_ddtil2a,
		    ITG *m_flagr,ITG *i_flagr,ITG *a_flagr,ITG *a_flag,
		    ITG *i_flag,ITG *m_flag,
		    ITG *row_ln,ITG *row_lm,ITG *row_li,ITG *row_la,
		    double *slavnor,double *slavtan,
		    double *vold,double *vini,double *cstress,
		    double *cstressini,
		    double *bp_old,ITG *nactdof,ITG *islavnode,
		    ITG *ntie,ITG *mi,ITG *nk,
		    ITG *nboun,ITG *ndirboun,ITG *nodeboun,double *xboun,
		    ITG *ipompc,ITG *nodempc,double *coefmpc,
		    ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
		    ITG *nslavspc,ITG *islavspc,ITG *nslavmpc,
		    ITG *islavmpc,char *tieset,
		    ITG *islavtie,ITG *nelcon,double  *elcon,
		    double *tietol,ITG *ncmat_,ITG *ntmat_,
		    double *plicon,ITG *nplicon,ITG *npmat_,double *dtime,
		    ITG *irowt,ITG *jqt,double *aut, 
		    ITG *irowtinv,ITG *jqtinv,double *autinv,
		    ITG *islavnodeinv,
		    ITG *iit,double *beta,ITG *ithermal,
		    double *plkcon,ITG *nplkcon){
  
  ITG i,j,jj,j2,k,l,idof1,idof2,idof3,iadd,jrow,islavnodeentry,
    mt=mi[1]+1,nodes,node,derivmode,regmode,
    *irow_antil=NULL,*irow_amtil=NULL,*irow_aitil=NULL,*irow_aatil=NULL,
    *irow_amtil1=NULL,*irow_amtil2=NULL,*irow_aitil1=NULL,*irow_aitil2=NULL,
    *irow_aatil1=NULL,*irow_aatil2=NULL,nzs_antil,nzs_amtil,nzs_aitil,
    nzs_aatil,*jq_amtil1=NULL,*jq_amtil2=NULL,nzs_amtil1,nzs_amtil2,nzs_aitil1,
    nzs_aitil2,nzs_aatil1,nzs_aatil2,*jq_aitil1=NULL,*jq_aitil2=NULL,
    *jq_aatil1=NULL,*jq_aatil2=NULL,ifree_antil,ifree_amtil1,ifree_amtil2,
    ifree_aitil1,ifree_aitil2,ifree_aatil1,ifree_aatil2;
  
  double t1,t2,e1,e2,e3,contribution,hpn,scal,bp,constant=1.E10,
    constantt=1.E10,lambda_n,*u_tilde=NULL,resreg[2],
    *cstress2=NULL,*cstressini2=NULL,that[6],n11,n22,aninvloc,gnc,dgnc,dgnc1,
    mu,p0,beta_e,atauinvloc,lambda_t[2],lambdaini_t[2],lambdatilde_t[2],ltu[2],
    ltslip[6],rphat[2],n[3],n2[3],t[6],utildep_t[2],rslip[6],*au_antil=NULL,
    *au_amtil=NULL,*au_aitil=NULL,*au_aatil=NULL,*au_amtil1=NULL,
    *au_amtil2=NULL,*au_aitil1=NULL,*au_aitil2=NULL,*au_aatil1=NULL,
    *au_aatil2=NULL;
  
  au_antil=*au_antilp;au_amtil=*au_amtilp;au_aitil=*au_aitilp;
  au_aatil=*au_aatilp;
  irow_antil= *irow_antilp;irow_amtil= *irow_amtilp;irow_aitil= *irow_aitilp;
  irow_aatil= *irow_aatilp;
  NNEW(u_tilde,double,3*nslavnode[*ntie]);
  NNEW(cstress2,double,mt*nslavnode[*ntie]);
  NNEW(cstressini2,double,mt*nslavnode[*ntie]);
  
  nzs_antil=jq_dan[*row_ln]-1;
  nzs_amtil=jq_dam[*row_lm]-1;
  nzs_aitil=jq_dai[*row_li]-1;
  nzs_aatil=jq_daa[*row_la]-1;
  
  /* get uhat_k-1=D_p*u^S+B_d*u^M and lambda_s,lambda_s^ini **/
  
  for (i=0;i<*ntie;i++){      	
    if(tieset[i*(81*3)+80]=='C'){	
      for(j=nslavnode[i];j<nslavnode[i+1];j++){	    
	nodes=islavnode[j];
	for(jj=jqdtil[nodes-1]-1;jj<jqdtil[nodes-1+1]-1;jj++){
	  if(islavnodeinv[irowdtil[jj]-1]>0 && islavnodeinv[nodes-1]>0){
	    for(l=0;l<3;l++){
	      cstress2[(islavnodeinv[irowdtil[jj]-1]-1)*mt+l]+=
		Ddtil[jj]*cstress[(islavnodeinv[nodes-1]-1)*mt+l];
	      cstressini2[(islavnodeinv[irowdtil[jj]-1]-1)*mt+l]+=
		Ddtil[jj]*cstressini[(islavnodeinv[nodes-1]-1)*mt+l];
	    }
	  }else{
	    printf("\ttrafoNTmortar: something went wrong in node %"
		   ITGFORMAT " or %" ITGFORMAT "\n",irowdtil[jj],nodes);
	    FORTRAN(stop,());	
	  }

	}
      }	
    }
  }  
  for(i=0;i<*nk;i++){
    nodes=i+1;
    for(jj=jqd[nodes-1]-1;jj<jqd[nodes-1+1]-1;jj++){
      if(islavnodeinv[irowd[jj]-1]>0 ){
	for(l=0;l<3;l++){
	  u_tilde[(islavnodeinv[irowd[jj]-1]-1)*3+l]+=
	    Dd[jj]*(vold[mt*(nodes)-3+l]-vini[mt*(nodes)-3+l]);
        
	}
      }else{
	printf("\ttrafoNTmortar: something went wrong in node %" ITGFORMAT
	       "\n",irowd[jj]);
	FORTRAN(stop,());
      }
    }
  }
  for(i=0;i<*nk;i++){
    nodes=i+1;
    for(jj=jqb[nodes-1]-1;jj<jqb[nodes-1+1]-1;jj++){
      if(islavnodeinv[irowb[jj]-1]>0 ){
	for(l=0;l<3;l++){
	  u_tilde[(islavnodeinv[irowb[jj]-1]-1)*3+l]+=
	    Bd[jj]*(vold[mt*(nodes)-3+l]-vini[mt*(nodes)-3+l]);
	}
      }else{
	printf("\ttrafoNTmortar: something went wrong in node %" ITGFORMAT
	       "\n",irowb[jj]);
	FORTRAN(stop,());
      }
    }
  }

  /* K_AN^til **/
  
  jq_antil[0]=1;
  ifree_antil=1;
  for(j=0;j<*row_ln;j++){
    //loop over columns  N   
    j2=j+1;
    for(i=jq_dan[j]-1;i<jq_dan[j+1]-1;i++){
      //loop over rows A  	  
      k=irow_dan[i]-1;         
      islavnodeentry=floor(islavactdof[a_flagr[k]-1]/10.);	  
      jrow= islavactdof[a_flagr[k]-1]-10*islavnodeentry; 
      node=islavnode[islavnodeentry-1];	   
      idof1=nactdof[mt*node-3]-1;	   
      idof2=nactdof[mt*node-2]-1;	   
      idof3=nactdof[mt*node-1]-1;
      for(l=0;l<3;l++){	     
	n[l]=slavnor[3*(islavnodeentry-1)+l];	     
	n2[l]=slavnor[3*(islavnodeentry-1)+l];	     	   
      }
      for(l=0;l<6;l++){	     
	t[l]=slavtan[6*(islavnodeentry-1)+l];	     
	that[l]=slavtan[6*(islavnodeentry-1)+l];	     	   
      }

      trafontspcmpc(n,t,n2,that,&islavnodeentry,nboun,ndirboun,nodeboun,xboun,
		    ipompc,nodempc,coefmpc,ikboun,ilboun,ikmpc,ilmpc,
		    nslavspc,islavspc,nslavmpc,islavmpc,
		    &node);
      
      /* calculate fields needed for Coulomb friction */
      
      utildep_t[0]=u_tilde[(islavnodeentry-1)*3]*that[0]+
	u_tilde[(islavnodeentry-1)*3+1]*that[1]+
	u_tilde[(islavnodeentry-1)*3+2]*that[2];
      utildep_t[1]=u_tilde[(islavnodeentry-1)*3]*that[3]+
	u_tilde[(islavnodeentry-1)*3+1]*that[4]+
	u_tilde[(islavnodeentry-1)*3+2]*that[5];
      lambda_n=n2[0]*cstress2[(islavnodeentry-1)*mt]+
	n2[1]*cstress2[(islavnodeentry-1)*mt+1]+
	n2[2]*cstress2[(islavnodeentry-1)*mt+2];
      lambda_t[0]=that[0]*cstress2[(islavnodeentry-1)*mt]+
	that[1]*cstress2[(islavnodeentry-1)*mt+1]+that[2]*
	cstress2[(islavnodeentry-1)*mt+2];
      lambda_t[1]=that[3]*cstress2[(islavnodeentry-1)*mt]+
	that[4]*cstress2[(islavnodeentry-1)*mt+1]+that[5]*
	cstress2[(islavnodeentry-1)*mt+2];
      lambdaini_t[0]=that[0]*cstressini2[(islavnodeentry-1)*mt]+
	that[1]*cstressini2[(islavnodeentry-1)*mt+1]+that[2]*
	cstressini2[(islavnodeentry-1)*mt+2];
      lambdaini_t[1]=that[3]*cstressini2[(islavnodeentry-1)*mt]+
	that[4]*cstressini2[(islavnodeentry-1)*mt+1]+that[5]*
	cstressini2[(islavnodeentry-1)*mt+2];
      lambdatilde_t[0]=lambda_t[0]-lambdaini_t[0];
      lambdatilde_t[1]=lambda_t[1]-lambdaini_t[1];	   
      bp=bp_old[islavnodeentry-1];
      
      /* perturbed lagrange,normal direction
	 (see phd-thesis Sitzmann,Chapter 3.4.1) */
      
      scal=Ddtil[jqdtil[node-1]-1];
      FORTRAN(getcontactparams,(&mu,&regmode,&aninvloc,&atauinvloc,
				&p0,&beta_e,tietol,elcon,
				&islavtie[islavnodeentry-1],ncmat_,
				ntmat_));
      constantt=min(constant,1.0/atauinvloc);
      derivmode=1;
      FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&dgnc,
				   &aninvloc,&p0,&beta_e,elcon,nelcon,
				   &islavtie[islavnodeentry-1],ntmat_,
				   plicon,nplicon,npmat_,ncmat_,tietol,&scal));
      
      if(islavact[islavnodeentry-1]==1 || islavact[islavnodeentry-1]==2){
	derivmode=1;
	  
	/* perturbed lagrange,tangential direction
	   (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	FORTRAN(regularization_slip_lin,(utildep_t,&bp,&atauinvloc,resreg,
					 &derivmode,islavact,lambda_t,
					 lambdatilde_t,&constantt,
					 &islavnodeentry,n2,t,that,&mu,rslip,
					 ltslip,ltu));
	rphat[0]=0.0;
	rphat[1]=0.0;
	
	if(jrow==4 && ithermal[0]<2){
	  printf(" *ERROR in trafontmortar2\n");
	  // something went wrong
	  FORTRAN(stop,());
	}else{
	  // mechanical part
	  if(idof1>-1 && idof2>-1 && idof3>-1){
	    // 3D	      	      
	    e1=au_dan[i];	      
	    e2=au_dan[i+1];	      
	    e3=au_dan[i+2];

	    contribution=(dgnc)*(n[0]*e1+n[1]*e2+n[2]*e3);
	    insertas_ws(&irow_antil,&(irow_dan[i]),&j2,&ifree_antil,
			&nzs_antil,&contribution,&au_antil);
	    
	    contribution=(rslip[0]*e1+rslip[1]*e2+rslip[2]*e3);
	    insertas_ws(&irow_antil,&(irow_dan[i+1]),&j2,&ifree_antil,
			&nzs_antil,&contribution,&au_antil);

	    contribution=(rslip[3]*e1+rslip[4]*e2+rslip[5]*e3);
	    insertas_ws(&irow_antil,&(irow_dan[i+2]),&j2,&ifree_antil,
			&nzs_antil,&contribution,&au_antil);
	    i=i+2;	     	       
	  }else if(idof2>-1 && idof3>-1){
	    //2D auf 3D
	    e1=au_dan[i];                    		
	    e2=au_dan[i+1];                   		
	    t1=rslip[4];				
	    t2=rslip[5];
	    n11=n2[1];
	    n22=n2[2];	     
	    contribution=(dgnc)*(n11*e1+n22*e2);
	    insertas_ws(&irow_antil,&(irow_dan[i]),&j2,&ifree_antil,
			&nzs_antil,&contribution,&au_antil);
	    contribution=(t1*e1+t2*e2);
	    insertas_ws(&irow_antil,&(irow_dan[i+1]),&j2,&ifree_antil,
			&nzs_antil,&contribution,&au_antil);
	    i=i+1;	       	      	     	       
	  }else if(idof1>-1 && idof3>-1){
	    //2D auf 3D
	    e1=au_dan[i];                    		
	    e2=au_dan[i+1];                    		
	    t1=rslip[3];				
	    t2=rslip[5];
	    n11=n2[0];
	    n22=n2[2];
	    contribution=(dgnc)*(n11*e1+n22*e2);
	    insertas_ws(&irow_antil,&(irow_dan[i]),&j2,&ifree_antil,
			&nzs_antil,&contribution,&au_antil);
	    contribution=(t1*e1+t2*e2);
	    insertas_ws(&irow_antil,&(irow_dan[i+1]),&j2,&ifree_antil,
			&nzs_antil,&contribution,&au_antil);
	    i=i+1;	       	      	       	     
	  }else if(idof1>-1 && idof2>-1){
	    //2D auf 3D
	    e1=au_dan[i];                    		
	    e2=au_dan[i+1];                    		
	    t1=rslip[3];				
	    t2=rslip[4];
	    n11=n2[0];
	    n22=n2[1];
	    contribution=(dgnc)*(n11*e1+n22*e2);
	    insertas_ws(&irow_antil,&(irow_dan[i]),&j2,&ifree_antil,
			&nzs_antil,&contribution,&au_antil);
	    contribution=(t1*e1+t2*e2);
	    insertas_ws(&irow_antil,&(irow_dan[i+1]),&j2,&ifree_antil,
			&nzs_antil,&contribution,&au_antil);
	    i=i+1;     	     
	  }else{
	    e1=au_dan[i];
	    if(idof1>-1){		 
	      n11=n2[0];
	    }else if(idof2>-1){
	      n11=n2[1];
	    }else{
	      n11=n2[2];
	    }
	    contribution=(dgnc)*(n11*e1+n22*e2);
	    insertas_ws(&irow_antil,&(irow_dan[i]),&j2,&ifree_antil,
			&nzs_antil,&contribution,&au_antil);
	  }
	}
      }else{
	printf("trafoNT2: something went wrong in K_dan!\n");
	FORTRAN(stop,());
      }	
    }
    jq_antil[j+1]=ifree_antil;
  }
  RENEW(irow_antil,ITG,ifree_antil-1);
  RENEW(au_antil,double,ifree_antil-1);
  
  /* K_AM^til **/
  
  ifree_amtil1=1;
  nzs_amtil1=jq_dam[*row_lm]-1;
  NNEW(irow_amtil1,ITG,nzs_amtil1);
  NNEW(jq_amtil1,ITG,*row_lm+1);
  NNEW(au_amtil1,double,nzs_amtil1);
  jq_amtil1[0]=1;
  for(j=0;j<*row_lm;j++){
    //loop over columns  M   
    j2=j+1;
    for(i=jq_dam[j]-1;i<jq_dam[j+1]-1;i++){
      //loop over rows A  	  
      k=irow_dam[i]-1;
      islavnodeentry=floor(islavactdof[a_flagr[k]-1]/10.);	  
      jrow= islavactdof[a_flagr[k]-1]-10*islavnodeentry; 
      node=islavnode[islavnodeentry-1];	   
      idof1=nactdof[mt*node-3]-1;	   
      idof2=nactdof[mt*node-2]-1;	   
      idof3=nactdof[mt*node-1]-1;
      for(l=0;l<3;l++){	     
	n[l]=slavnor[3*(islavnodeentry-1)+l];	     
	n2[l]=slavnor[3*(islavnodeentry-1)+l];	     	   
      }
      for(l=0;l<6;l++){	     
	t[l]=slavtan[6*(islavnodeentry-1)+l];	     
	that[l]=slavtan[6*(islavnodeentry-1)+l];	     	   
      }
      trafontspcmpc(n,t,n2,that,&islavnodeentry,
		    nboun,ndirboun,nodeboun,xboun,
		    ipompc,nodempc,coefmpc,
		    ikboun,ilboun,ikmpc,ilmpc,
		    nslavspc,islavspc,nslavmpc,islavmpc,
		    &node);
      
      /* calculate fields needed for Coulomb friction **/
      
      utildep_t[0]=u_tilde[(islavnodeentry-1)*3]*that[0]+
	u_tilde[(islavnodeentry-1)*3+1]*that[1]+
	u_tilde[(islavnodeentry-1)*3+2]*that[2];
      utildep_t[1]=u_tilde[(islavnodeentry-1)*3]*that[3]+
	u_tilde[(islavnodeentry-1)*3+1]*that[4]+
	u_tilde[(islavnodeentry-1)*3+2]*that[5];
      lambda_n=n2[0]*cstress2[(islavnodeentry-1)*mt]+
	n2[1]*cstress2[(islavnodeentry-1)*mt+1]+
	n2[2]*cstress2[(islavnodeentry-1)*mt+2];
      lambda_t[0]=that[0]*cstress2[(islavnodeentry-1)*mt]+
	that[1]*cstress2[(islavnodeentry-1)*mt+1]+
	that[2]*cstress2[(islavnodeentry-1)*mt+2];
      lambda_t[1]=that[3]*cstress2[(islavnodeentry-1)*mt]+
	that[4]*cstress2[(islavnodeentry-1)*mt+1]+
	that[5]*cstress2[(islavnodeentry-1)*mt+2];
      lambdaini_t[0]=that[0]*cstressini2[(islavnodeentry-1)*mt]+
	that[1]*cstressini2[(islavnodeentry-1)*mt+1]+
	that[2]*cstressini2[(islavnodeentry-1)*mt+2];
      lambdaini_t[1]=that[3]*cstressini2[(islavnodeentry-1)*mt]+
	that[4]*cstressini2[(islavnodeentry-1)*mt+1]+
	that[5]*cstressini2[(islavnodeentry-1)*mt+2];
      lambdatilde_t[0]=lambda_t[0]-lambdaini_t[0];
      lambdatilde_t[1]=lambda_t[1]-lambdaini_t[1];	   
      bp=bp_old[islavnodeentry-1];
      
      /* perturbed lagrange,normal direction  
	 (see phd-thesis Sitzmann,Chapter 3.4.1) */
      
      scal=Ddtil[jqdtil[node-1]-1];
      FORTRAN(getcontactparams,(&mu,&regmode,&aninvloc,&atauinvloc,
				&p0,&beta_e,tietol,elcon,
				&islavtie[islavnodeentry-1],ncmat_,
				ntmat_));
      
      constantt=min(constant,1.0/atauinvloc);
      derivmode=1;
      
      FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&dgnc,
				   &aninvloc,&p0,&beta_e,elcon,nelcon,
				   &islavtie[islavnodeentry-1],ntmat_,
				   plicon,nplicon,npmat_,ncmat_,tietol,&scal));
      
      if(islavact[islavnodeentry-1]==1 || islavact[islavnodeentry-1]==2){
	derivmode=1;
	  
	/* perturbed lagrange,tangential direction  
	   (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	FORTRAN(regularization_slip_lin,(utildep_t,&bp,&atauinvloc,resreg,
					 &derivmode,islavact,lambda_t,
					 lambdatilde_t,&constantt,
					 &islavnodeentry,n2,t,that,&mu,rslip,
					 ltslip,ltu));
	  
	rphat[0]=0.0;
	rphat[1]=0.0;
	if(jrow==4 && ithermal[0]<2){ 
	  printf(" *ERROR in trafontmortar2\n");
	  // something went wrong
	  FORTRAN(stop,());
	}else{
	
	  if(idof1>-1 && idof2>-1 && idof3>-1){
	    // 3D	      	      
	    e1=au_dam[i];	      
	    e2=au_dam[i+1];	      
	    e3=au_dam[i+2];	      		     
	    contribution=(dgnc)*(n[0]*e1+n[1]*e2+n[2]*e3);
	    insertas_ws(&irow_amtil1,&(irow_dam[i]),&j2,&ifree_amtil1,
			&nzs_amtil1,&contribution,&au_amtil1);
	    contribution=(rslip[0]*e1+rslip[1]*e2+rslip[2]*e3);
	    insertas_ws(&irow_amtil1,&(irow_dam[i+1]),&j2,&ifree_amtil1,
			&nzs_amtil1,&contribution,&au_amtil1);
	    contribution=(rslip[3]*e1+rslip[4]*e2+rslip[5]*e3);
	    insertas_ws(&irow_amtil1,&(irow_dam[i+2]),&j2,&ifree_amtil1,
			&nzs_amtil1,&contribution,&au_amtil1);
	    i=i+2;	     	       
	  }else if(idof2>-1 && idof3>-1){
	    //2D auf 3D
	    e1=au_dam[i];                    		
	    e2=au_dam[i+1];                   		
	    t1=rslip[4];				
	    t2=rslip[5];
	    n11=n2[1];
	    n22=n2[2];	     
	    contribution=(dgnc)*(n11*e1+n22*e2);
	    insertas_ws(&irow_amtil1,&(irow_dam[i]),&j2,&ifree_amtil1,
			&nzs_amtil1,&contribution,&au_amtil1);
	    contribution=(t1*e1+t2*e2);
	    insertas_ws(&irow_amtil1,&(irow_dam[i+1]),&j2,&ifree_amtil1,
			&nzs_amtil1,&contribution,&au_amtil1);
	    i=i+1;	       	      	     	       
	  }else if(idof1>-1 && idof3>-1){
	    //2D auf 3D
	    e1=au_dam[i];                    		
	    e2=au_dam[i+1];                    		
	    t1=rslip[3];				
	    t2=rslip[5];
	    n11=n2[0];
	    n22=n2[2];
	    contribution=(dgnc)*(n11*e1+n22*e2);
	    insertas_ws(&irow_amtil1,&(irow_dam[i]),&j2,&ifree_amtil1,
			&nzs_amtil1,&contribution,&au_amtil1);
	    contribution=(t1*e1+t2*e2);
	    insertas_ws(&irow_amtil1,&(irow_dam[i+1]),&j2,&ifree_amtil1,
			&nzs_amtil1,&contribution,&au_amtil1);
	    i=i+1;	       	      	       	     
	  }else if(idof1>-1 && idof2>-1){
	    //2D auf 3D
	    e1=au_dam[i];                    		
	    e2=au_dam[i+1];                    		
	    t1=rslip[3];				
	    t2=rslip[4];
	    n11=n2[0];
	    n22=n2[1];
	    contribution=(dgnc)*(n11*e1+n22*e2);
	    insertas_ws(&irow_amtil1,&(irow_dam[i]),&j2,&ifree_amtil1,
			&nzs_amtil1,&contribution,&au_amtil1);
	    contribution=(t1*e1+t2*e2);
	    insertas_ws(&irow_amtil1,&(irow_dam[i+1]),&j2,&ifree_amtil1,
			&nzs_amtil1,&contribution,&au_amtil1);
	    i=i+1;     	     
	  }else {
	    e1=au_dam[i];
	    if(idof1>-1){		 
	      n11=n2[0];
	    }else if(idof2>-1){
	      n11=n2[1];
	    }else{
	      n11=n2[2];
	    }
	    contribution=(dgnc)*(n11*e1);
	    insertas_ws(&irow_amtil1,&(irow_dam[i]),&j2,&ifree_amtil1,
			&nzs_amtil1,&contribution,&au_amtil1);
	  }
	}
      }else{
	printf("trafoNT2: something went wrong in K_dam!\n");
	FORTRAN(stop,());
      }	
    }
    jq_amtil1[j+1]=ifree_amtil1;
  }
  
  /* add diagonal terms **/
  
  nzs_amtil2=3*(jq_bdtil2[*row_lm]-1);
  NNEW(au_amtil2,double,3*(jq_bdtil2[*row_lm]-1));
  NNEW(irow_amtil2,ITG,3*(jq_bdtil2[*row_lm]-1));
  NNEW(jq_amtil2,ITG,*row_lm+1);
  ifree_amtil2=1;
  jq_amtil2[0]=1;
  for(j=0;j<*row_lm;j++){
    //loop over columns  M   
    j2=j+1;
    for(i=jq_bdtil2[j]-1;i<jq_bdtil2[j+1]-1;i++){
      //loop over rows A  
      k=irow_bdtil2[i]-1; 
      islavnodeentry=floor(islavactdof[a_flagr[k]-1]/10.);	  
      jrow= islavactdof[a_flagr[k]-1]-10*islavnodeentry; 
      node=islavnode[islavnodeentry-1];	   
      idof1=nactdof[mt*node-3]-1;	   
      idof2=nactdof[mt*node-2]-1;	   
      idof3=nactdof[mt*node-1]-1;
      for(l=0;l<3;l++){	     
	n[l]=slavnor[3*(islavnodeentry-1)+l];	     
	n2[l]=slavnor[3*(islavnodeentry-1)+l];	     	   
      }
      for(l=0;l<6;l++){	     
	t[l]=slavtan[6*(islavnodeentry-1)+l];	     
	that[l]=slavtan[6*(islavnodeentry-1)+l];	     	   
      }
      trafontspcmpc(n,t,n2,that,&islavnodeentry,
		    nboun,ndirboun,nodeboun,xboun,
		    ipompc,nodempc,coefmpc,
		    ikboun,ilboun,ikmpc,ilmpc,
		    nslavspc,islavspc,nslavmpc,islavmpc,
		    &node);
      
      /* calculate fields needed for Coulomb friction **/
      
      utildep_t[0]=u_tilde[(islavnodeentry-1)*3]*that[0]+
	u_tilde[(islavnodeentry-1)*3+1]*that[1]+
	u_tilde[(islavnodeentry-1)*3+2]*that[2];
      utildep_t[1]=u_tilde[(islavnodeentry-1)*3]*that[3]+
	u_tilde[(islavnodeentry-1)*3+1]*that[4]+
	u_tilde[(islavnodeentry-1)*3+2]*that[5];
      lambda_n=n2[0]*cstress2[(islavnodeentry-1)*mt]+
	n2[1]*cstress2[(islavnodeentry-1)*mt+1]+
	n2[2]*cstress2[(islavnodeentry-1)*mt+2];
      lambda_t[0]=that[0]*cstress2[(islavnodeentry-1)*mt]+
	that[1]*cstress2[(islavnodeentry-1)*mt+1]+
	that[2]*cstress2[(islavnodeentry-1)*mt+2];
      lambda_t[1]=that[3]*cstress2[(islavnodeentry-1)*mt]+
	that[4]*cstress2[(islavnodeentry-1)*mt+1]+
	that[5]*cstress2[(islavnodeentry-1)*mt+2];
      lambdaini_t[0]=that[0]*cstressini2[(islavnodeentry-1)*mt]+
	that[1]*cstressini2[(islavnodeentry-1)*mt+1]+
	that[2]*cstressini2[(islavnodeentry-1)*mt+2];
      lambdaini_t[1]=that[3]*cstressini2[(islavnodeentry-1)*mt]+
	that[4]*cstressini2[(islavnodeentry-1)*mt+1]+
	that[5]*cstressini2[(islavnodeentry-1)*mt+2];
      lambdatilde_t[0]=lambda_t[0]-lambdaini_t[0];
      lambdatilde_t[1]=lambda_t[1]-lambdaini_t[1];	   
      bp=bp_old[islavnodeentry-1];
      
      /* perturbed lagrange,normal direction 
	 (see phd-thesis Sitzmann,Chapter 3.4.1) */
      
      scal=Ddtil[jqdtil[node-1]-1];
      FORTRAN(getcontactparams,(&mu,&regmode,&aninvloc,&atauinvloc,
				&p0,&beta_e,tietol,elcon,
				&islavtie[islavnodeentry-1],ncmat_,
				ntmat_));
      
      constantt=min(constant,1.0/atauinvloc);
      derivmode=1;
      FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&dgnc,
				   &aninvloc,&p0,&beta_e,elcon,nelcon,
				   &islavtie[islavnodeentry-1],ntmat_,
				   plicon,nplicon,npmat_,ncmat_,tietol,&scal));
      
      if(islavact[islavnodeentry-1]==1 || islavact[islavnodeentry-1]==2){
	derivmode=1;
	  
	/* perturbed lagrange,tangential direction 
	   (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	FORTRAN(regularization_slip_lin,(utildep_t,&bp,&atauinvloc,resreg,
					 &derivmode,islavact,lambda_t,
					 lambdatilde_t,&constantt,
					 &islavnodeentry,n2,t,that,&mu,rslip,
					 ltslip,ltu));
	  
	rphat[0]=0.0;
	rphat[1]=0.0;
	if(jrow==4 && ithermal[0]<2){ 
	  printf(" *ERROR in trafontmortar2\n");
	  // something went wrong
	  FORTRAN(stop,());
	}else{
	  //mechanical part
	  if(idof1>-1 && idof2>-1 && idof3>-1){
	    // 3D	
	    if(jrow==1){
	      // k=idof1
	      iadd=0;
	      e1=au_bdtil2[i];
	      if(i+1<jq_bdtil2[j+1]-1 && a_flagr[irow_bdtil2[i+1]-1]-1==idof2){
		e2=au_bdtil2[i+1];++iadd;}else{e2=0.0;}
	      if(i+1<jq_bdtil2[j+1]-1 && a_flagr[irow_bdtil2[i+1]-1]-1==idof3){
		e3=au_bdtil2[i+1];++iadd;
	      }else if(i+2<jq_bdtil2[j+1]-1 &&
		       a_flagr[irow_bdtil2[i+2]-1]-1==idof3){
		e3=au_bdtil2[i+2];++iadd;
	      }else{e3=0.0;}
	    }else if(jrow==2){
	      // k=idof2
	      iadd=0;
	      e1=0.0;
	      e2=au_bdtil2[i];
	      if(i+1<jq_bdtil2[j+1]-1 && a_flagr[irow_bdtil2[i+1]-1]-1==idof3){
		e3=au_bdtil2[i+1];++iadd;}else{e3=0.0;}
	    }else{
	      // k=idof3
	      iadd=0;
	      e1=0.0;
	      e2=0.0;
	      e3=au_bdtil2[i];
	    }
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    insertas_ws(&irow_amtil2,&(a_flag[idof1]),&(j2),&ifree_amtil2,
			&nzs_amtil2,&contribution,&au_amtil2);
	    contribution=(ltslip[0]*e1+ltslip[1]*e2+ltslip[2]*e3);
	    insertas_ws(&irow_amtil2,&(a_flag[idof2]),&(j2),&ifree_amtil2,
			&nzs_amtil2,&contribution,&au_amtil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    insertas_ws(&irow_amtil2,&(a_flag[idof3]),&(j2),&ifree_amtil2,
			&nzs_amtil2,&contribution,&au_amtil2);
	    i=i+iadd;    
	  }else if(idof2>-1 && idof3>-1){
	    //2D auf 3D
	    if(jrow==2){
	      // k=idof2
	      iadd=0;
	      e1=0.0;
	      e2=au_bdtil2[i];
	      if(i+1<jq_bdtil2[j+1]-1 &&
		 a_flagr[irow_bdtil2[i+1]-1]-1==idof3){
		e3=au_bdtil2[i+1];++iadd;}else{e3=0.0;}
	    }else{
	      // k=idof3
	      iadd=0;
	      e1=0.0;
	      e2=0.0;
	      e3=au_bdtil2[i];
	    }	  
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    insertas_ws(&irow_amtil2,&(a_flag[idof2]),&(j2),&ifree_amtil2,
			&nzs_amtil2,&contribution,&au_amtil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    insertas_ws(&irow_amtil2,&(a_flag[idof3]),&(j2),&ifree_amtil2,
			&nzs_amtil2,&contribution,&au_amtil2);
	    i=i+iadd;    	       	     
	  }else if(idof1>-1 && idof3>-1){
	    //2D auf 3D	
	    if(jrow==1){
	      // k=idof1
	      iadd=0;
	      e1=au_bdtil2[i];
	      e2=0.0;
	      if(i+1<jq_bdtil2[j+1]-1 && a_flagr[irow_bdtil2[i+1]-1]-1==idof3){
		e3=au_bdtil2[i+1];++iadd;
	      }else{
		e3=0.0;}
	    }else{
	      // k=idof3
	      iadd=0;
	      e1=0.0;
	      e2=0.0;
	      e3=au_bdtil2[i];
	    }
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    insertas_ws(&irow_amtil2,&(a_flag[idof1]),&(j2),&ifree_amtil2,
			&nzs_amtil2,&contribution,&au_amtil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    insertas_ws(&irow_amtil2,&(a_flag[idof3]),&(j2),&ifree_amtil2,
			&nzs_amtil2,&contribution,&au_amtil2);
	    i=i+iadd;       	     
	  }else if(idof1>-1 && idof2>-1){
	    //2D auf 3D	
	    if(jrow==1){
	      // k=idof1
	      iadd=0;
	      e1=au_bdtil2[i];
	      if(i+1<jq_bdtil2[j+1]-1 && a_flagr[irow_bdtil2[i+1]-1]-1==idof2){
		e2=au_bdtil2[i+1];
		++iadd;
	      }else{
		e2=0.0;}
	      e3=0.0;
	    }else{
	      // k=idof3
	      iadd=0;
	      e1=0.0;
	      e2=au_bdtil2[i];
	      e3=0.0;
	    }
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    insertas_ws(&irow_amtil2,&(a_flag[idof1]),&(j2),&ifree_amtil2,
			&nzs_amtil2,&contribution,&au_amtil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    insertas_ws(&irow_amtil2,&(a_flag[idof2]),&(j2),&ifree_amtil2,
			&nzs_amtil2,&contribution,&au_amtil2);
	    i=i+iadd;		           
	  }else {
	    // 1D auf 3D                   
	    e1=0.0;e2=0.0;e3=0.0;
	    if(idof1>-1){
	      e1=au_bdtil2[i];
	      contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	      insertas_ws(&irow_amtil2,&(a_flag[idof1]),&(j2),&ifree_amtil2,
			  &nzs_amtil2,&contribution,&au_amtil2);	     
	    }else if(idof2>-1){
	      e2=au_bdtil2[i];
	      contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	      insertas_ws(&irow_amtil2,&(a_flag[idof2]),&(j2),&ifree_amtil2,
			  &nzs_amtil2,&contribution,&au_amtil2);		
	    }else{
	      e3=au_bdtil2[i];
	      contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	      insertas_ws(&irow_amtil2,&(a_flag[idof3]),&(j2),&ifree_amtil2,
			  &nzs_amtil2,&contribution,&au_amtil2);
	    }   
	  }
	}
      }
    }
    jq_amtil2[j+1]=ifree_amtil2;
  }
  
  nzs_amtil2=ifree_amtil2-1;

  add_rect(au_amtil1,irow_amtil1,jq_amtil1,*row_la,*row_lm,
	   au_amtil2,irow_amtil2,jq_amtil2,*row_la,*row_lm,
	   &au_amtil,&irow_amtil,jq_amtil,&nzs_amtil);
  
  SFREE(au_amtil1);SFREE(irow_amtil1);SFREE(jq_amtil1);
  SFREE(au_amtil2);SFREE(irow_amtil2);SFREE(jq_amtil2);
  
  /* K_AI **/
  
  nzs_aitil1=jq_dai[*row_li];
  NNEW(irow_aitil1,ITG,nzs_aitil1);
  NNEW(jq_aitil1,ITG,*row_li+1);
  NNEW(au_aitil1,double,nzs_aitil1);
  jq_aitil1[0]=1;
  ifree_aitil1=1;
  for(j=0;j<*row_li;j++){
    //loop over columns  N   
    for(i=jq_dai[j]-1;i<jq_dai[j+1]-1;i++){
      //loop over rows A  	  
      k=irow_dai[i]-1;         
      islavnodeentry=floor(islavactdof[a_flagr[k]-1]/10.);	  
      jrow= islavactdof[a_flagr[k]-1]-10*islavnodeentry; 
      node=islavnode[islavnodeentry-1];	   
      idof1=nactdof[mt*node-3]-1;	   
      idof2=nactdof[mt*node-2]-1;	   
      idof3=nactdof[mt*node-1]-1;
      for(l=0;l<3;l++){	     
	n[l]=slavnor[3*(islavnodeentry-1)+l];	     
	n2[l]=slavnor[3*(islavnodeentry-1)+l];	     	   
      }
      for(l=0;l<6;l++){	     
	t[l]=slavtan[6*(islavnodeentry-1)+l];	     
	that[l]=slavtan[6*(islavnodeentry-1)+l];	     	   
      }

      trafontspcmpc(n,t,n2,that,&islavnodeentry,nboun,ndirboun,nodeboun,xboun,
		    ipompc,nodempc,coefmpc,ikboun,ilboun,ikmpc,ilmpc,
		    nslavspc,islavspc,nslavmpc,islavmpc,
		    &node);
      
      /* calculate fields needed for Coulomb friction **/
      
      utildep_t[0]=u_tilde[(islavnodeentry-1)*3]*that[0]+
	u_tilde[(islavnodeentry-1)*3+1]*that[1]+
	u_tilde[(islavnodeentry-1)*3+2]*that[2];
      utildep_t[1]=u_tilde[(islavnodeentry-1)*3]*that[3]+
	u_tilde[(islavnodeentry-1)*3+1]*that[4]+
	u_tilde[(islavnodeentry-1)*3+2]*that[5];
      lambda_n=n2[0]*cstress2[(islavnodeentry-1)*mt]+
	n2[1]*cstress2[(islavnodeentry-1)*mt+1]+
	n2[2]*cstress2[(islavnodeentry-1)*mt+2];
      lambda_t[0]=that[0]*cstress2[(islavnodeentry-1)*mt]+
	that[1]*cstress2[(islavnodeentry-1)*mt+1]+
	that[2]*cstress2[(islavnodeentry-1)*mt+2];
      lambda_t[1]=that[3]*cstress2[(islavnodeentry-1)*mt]+
	that[4]*cstress2[(islavnodeentry-1)*mt+1]+
	that[5]*cstress2[(islavnodeentry-1)*mt+2];
      lambdaini_t[0]=that[0]*cstressini2[(islavnodeentry-1)*mt]+
	that[1]*cstressini2[(islavnodeentry-1)*mt+1]+
	that[2]*cstressini2[(islavnodeentry-1)*mt+2];
      lambdaini_t[1]=that[3]*cstressini2[(islavnodeentry-1)*mt]+
	that[4]*cstressini2[(islavnodeentry-1)*mt+1]+
	that[5]*cstressini2[(islavnodeentry-1)*mt+2];
      lambdatilde_t[0]=lambda_t[0]-lambdaini_t[0];
      lambdatilde_t[1]=lambda_t[1]-lambdaini_t[1];	   
      bp=bp_old[islavnodeentry-1];
      
      /* perturbed lagrange,normal direction 
	 (see phd-thesis Sitzmann,Chapter 3.4.1) */
      
      scal=Ddtil[jqdtil[node-1]-1];
      FORTRAN(getcontactparams,(&mu,&regmode,&aninvloc,&atauinvloc,
				&p0,&beta_e,tietol,elcon,
				&islavtie[islavnodeentry-1],ncmat_,
				ntmat_));
      
      constantt=min(constant,1.0/atauinvloc);
      derivmode=1;
      
      FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&dgnc,
				   &aninvloc,&p0,&beta_e,elcon,nelcon,
				   &islavtie[islavnodeentry-1],ntmat_,
				   plicon,nplicon,npmat_,ncmat_,tietol,&scal));

      if(islavact[islavnodeentry-1]==1 || islavact[islavnodeentry-1]==2){
	derivmode=1;
	  
	/* perturbed lagrange,tangential direction 
	   (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	FORTRAN(regularization_slip_lin,(utildep_t,&bp,&atauinvloc,resreg,
					 &derivmode,islavact,lambda_t,
					 lambdatilde_t,&constantt,
					 &islavnodeentry,n2,t,that,&mu,rslip,
					 ltslip,ltu));
	  
	rphat[0]=0.0;
	rphat[1]=0.0;
	if(jrow==4 && ithermal[0]<2){ 
	  printf(" *ERROR in trafontmortar2\n");
	  // something went wrong
	  FORTRAN(stop,());
	}else{
	  // mechanical part
	  if(idof1>-1 && idof2>-1 && idof3>-1){
	    // 3D	      	      
	    e1=au_dai[i];	      
	    e2=au_dai[i+1];	      
	    e3=au_dai[i+2];	      		     
	    contribution=(dgnc)*(n[0]*e1+n[1]*e2+n[2]*e3);
	    insertas_ws(&irow_aitil1,&(irow_dai[i]),&j2,&ifree_aitil1,
			&nzs_aitil1,&contribution,&au_aitil1);
	    contribution=(rslip[0]*e1+rslip[1]*e2+rslip[2]*e3);
	    insertas_ws(&irow_aitil1,&(irow_dai[i+1]),&j2,&ifree_aitil1,
			&nzs_aitil1,&contribution,&au_aitil1);
	    contribution=(rslip[3]*e1+rslip[4]*e2+rslip[5]*e3);
	    insertas_ws(&irow_aitil1,&(irow_dai[i+2]),&j2,&ifree_aitil1,
			&nzs_aitil1,&contribution,&au_aitil1);
	    i=i+2;	     	       
	  }else if(idof2>-1 && idof3>-1){
	    //2D auf 3D
	    e1=au_dai[i];                    		
	    e2=au_dai[i+1];                   		
	    t1=rslip[4];				
	    t2=rslip[5];
	    n11=n2[1];
	    n22=n2[2];	     
	    contribution=(dgnc)*(n11*e1+n22*e2);
	    insertas_ws(&irow_aitil1,&(irow_dai[i]),&j2,&ifree_aitil1,
			&nzs_aitil1,&contribution,&au_aitil1);
	    contribution=(t1*e1+t2*e2);
	    insertas_ws(&irow_aitil1,&(irow_dai[i+1]),&j2,&ifree_aitil1,
			&nzs_aitil1,&contribution,&au_aitil1);
	    i=i+1;	       	      	     	       
	  }else if(idof1>-1 && idof3>-1){
	    //2D auf 3D
	    e1=au_dai[i];                    		
	    e2=au_dai[i+1];                    		
	    t1=rslip[3];				
	    t2=rslip[5];
	    n11=n2[0];
	    n22=n2[2];
	    contribution=(dgnc)*(n11*e1+n22*e2);
	    insertas_ws(&irow_aitil1,&(irow_dai[i]),&j2,&ifree_aitil1,
			&nzs_aitil1,&contribution,&au_aitil1);
	    contribution=(t1*e1+t2*e2);
	    insertas_ws(&irow_aitil1,&(irow_dai[i+1]),&j2,&ifree_aitil1,
			&nzs_aitil1,&contribution,&au_aitil1);
	    i=i+1;	       	      	       	     
	  }else if(idof1>-1 && idof2>-1){
	    //2D auf 3D
	    e1=au_dai[i];                    		
	    e2=au_dai[i+1];                    		
	    t1=rslip[3];				
	    t2=rslip[4];
	    n11=n2[0];
	    n22=n2[1];
	    contribution=(dgnc)*(n11*e1+n22*e2);
	    insertas_ws(&irow_aitil1,&(irow_dai[i]),&j2,&ifree_aitil1,
			&nzs_aitil1,&contribution,&au_aitil1);
	    contribution=(t1*e1+t2*e2);
	    insertas_ws(&irow_aitil1,&(irow_dai[i+1]),&j2,&ifree_aitil1,
			&nzs_aitil1,&contribution,&au_aitil1);
	    i=i+1;     	     
	  }else {
	    e1=au_dai[i];
	    if(idof1>-1){		 
	      n11=n2[0];
	    }else if(idof2>-1){
	      n11=n2[1];
	    }else{
	      n11=n2[2];
	    }
	    contribution=(dgnc)*(n11*e1);
	    insertas_ws(&irow_aitil1,&(irow_dai[i]),&j2,&ifree_aitil1,
			&nzs_aitil1,&contribution,&au_aitil1);
	  }
	}
      }else{
	printf("\ttrafoNT2: something went wrong in K_dai!\n");
	FORTRAN(stop,());
      }	
    }
    jq_aitil1[j+1]=ifree_aitil1;
  } 
  RENEW(irow_aitil1,ITG,ifree_aitil1-1);
  RENEW(au_aitil1,double,ifree_aitil1-1);

  /* add diagonal terms **/

  nzs_aitil2=3*(jq_ddtil2i[*row_li]-1);
  NNEW(au_aitil2,double,3*(jq_ddtil2i[*row_li]-1));
  NNEW(irow_aitil2,ITG,3*(jq_ddtil2i[*row_li]-1));
  NNEW(jq_aitil2,ITG,*row_li+1);
  ifree_aitil2=1;
  jq_aitil2[0]=1;
  for(j=0;j<*row_li;j++){
    //loop over columns  I   
    j2=j+1;
    for(i=jq_ddtil2i[j]-1;i<jq_ddtil2i[j+1]-1;i++){
      //loop over rows A  
      k=irow_ddtil2i[i]-1; 
      islavnodeentry=floor(islavactdof[a_flagr[k]-1]/10.);	  
      jrow= islavactdof[a_flagr[k]-1]-10*islavnodeentry; 
      node=islavnode[islavnodeentry-1];	   
      idof1=nactdof[mt*node-3]-1;	   
      idof2=nactdof[mt*node-2]-1;	   
      idof3=nactdof[mt*node-1]-1;
      for(l=0;l<3;l++){	     
	n[l]=slavnor[3*(islavnodeentry-1)+l];	     
	n2[l]=slavnor[3*(islavnodeentry-1)+l];	     	   
      }
      for(l=0;l<6;l++){	     
	t[l]=slavtan[6*(islavnodeentry-1)+l];	     
	that[l]=slavtan[6*(islavnodeentry-1)+l];	     	   
      }
      
      trafontspcmpc(n,t,n2,that,&islavnodeentry,nboun,ndirboun,nodeboun,xboun,
		    ipompc,nodempc,coefmpc,ikboun,ilboun,ikmpc,ilmpc,
		    nslavspc,islavspc,nslavmpc,islavmpc,
		    &node);
      
      
      /* calculate fields needed for Coulomb friction **/
      
      utildep_t[0]=u_tilde[(islavnodeentry-1)*3]*
	that[0]+u_tilde[(islavnodeentry-1)*3+1]*
	that[1]+u_tilde[(islavnodeentry-1)*3+2]*that[2];
      utildep_t[1]=u_tilde[(islavnodeentry-1)*3]*that[3]+
	u_tilde[(islavnodeentry-1)*3+1]*that[4]+
	u_tilde[(islavnodeentry-1)*3+2]*that[5];
      lambda_n=n2[0]*cstress2[(islavnodeentry-1)*mt]+
	n2[1]*cstress2[(islavnodeentry-1)*mt+1]+
	n2[2]*cstress2[(islavnodeentry-1)*mt+2];
      lambda_t[0]=that[0]*cstress2[(islavnodeentry-1)*mt]+
	that[1]*cstress2[(islavnodeentry-1)*mt+1]+
	that[2]*cstress2[(islavnodeentry-1)*mt+2];
      lambda_t[1]=that[3]*cstress2[(islavnodeentry-1)*mt]+
	that[4]*cstress2[(islavnodeentry-1)*mt+1]+
	that[5]*cstress2[(islavnodeentry-1)*mt+2];
      lambdaini_t[0]=that[0]*cstressini2[(islavnodeentry-1)*mt]+
	that[1]*cstressini2[(islavnodeentry-1)*mt+1]+
	that[2]*cstressini2[(islavnodeentry-1)*mt+2];
      lambdaini_t[1]=that[3]*cstressini2[(islavnodeentry-1)*mt]+
	that[4]*cstressini2[(islavnodeentry-1)*mt+1]+
	that[5]*cstressini2[(islavnodeentry-1)*mt+2];
      lambdatilde_t[0]=lambda_t[0]-lambdaini_t[0];
      lambdatilde_t[1]=lambda_t[1]-lambdaini_t[1];	   
      bp=bp_old[islavnodeentry-1];
      
      /* perturbed lagrange,normal direction 
	 (see phd-thesis Sitzmann,Chapter 3.4.1) */
      
      scal=Ddtil[jqdtil[node-1]-1];
      FORTRAN(getcontactparams,(&mu,&regmode,&aninvloc,&atauinvloc,
				&p0,&beta_e,tietol,elcon,
				&islavtie[islavnodeentry-1],ncmat_,
				ntmat_));
      
      constantt=min(constant,1.0/atauinvloc);
      derivmode=1;
      FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&dgnc,
				   &aninvloc,&p0,&beta_e,elcon,nelcon,
				   &islavtie[islavnodeentry-1],ntmat_,
				   plicon,nplicon,npmat_,ncmat_,tietol,&scal));
      
      if(islavact[islavnodeentry-1]==1 || islavact[islavnodeentry-1]==2){ 	    
	derivmode=1;
	  
	/* perturbed lagrange,tangential direction 
	   (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	FORTRAN(regularization_slip_lin,(utildep_t,&bp,&atauinvloc,resreg,
					 &derivmode,islavact,lambda_t,
					 lambdatilde_t,&constantt,
					 &islavnodeentry,n2,t,that,&mu,rslip,
					 ltslip,ltu));
	  
	rphat[0]=0.0;
	rphat[1]=0.0;
	if(jrow==4 && ithermal[0]<2){ 
	  printf(" *ERROR in trafontmortar2\n");
	  // something went wrong
	  FORTRAN(stop,());
	}else{
	  // mechanical part	
	  if(idof1>-1 && idof2>-1 && idof3>-1){
	    // 3D	
	    if(jrow==1){
	      // k=idof1
	      iadd=0;
	      e1=au_ddtil2i[i];
	      if(i+1<jq_ddtil2i[j+1]-1 &&
		 a_flagr[irow_ddtil2i[i+1]-1]-1==idof2){
		e2=au_ddtil2i[i+1];++iadd;}else{e2=0.0;}
	      if(i+1<jq_ddtil2i[j+1]-1 &&
		 a_flagr[irow_ddtil2i[i+1]-1]-1==idof3){
		e3=au_ddtil2i[i+1];++iadd;
	      }else if(i+2<jq_ddtil2i[j+1]-1 &&
		       a_flagr[irow_ddtil2i[i+2]-1]-1==idof3){
		e3=au_ddtil2i[i+2];++iadd;}else{e3=0.0;}
	    }else if(jrow==2){
	      // k=idof2
	      iadd=0;
	      e1=0.0;
	      e2=au_ddtil2i[i];
	      if(i+1<jq_ddtil2i[j+1]-1 &&
		 a_flagr[irow_ddtil2i[i+1]-1]-1==idof3){
		e3=au_ddtil2i[i+1];++iadd;}else{e3=0.0;}
	    }else{
	      // k=idof3
	      iadd=0;
	      e1=0.0;
	      e2=0.0;
	      e3=au_ddtil2i[i];
	    }
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    insertas_ws(&irow_aitil2,&(a_flag[idof1]),&(j2),&ifree_aitil2,
			&nzs_aitil2,&contribution,&au_aitil2);
	    contribution=(ltslip[0]*e1+ltslip[1]*e2+ltslip[2]*e3);
	    insertas_ws(&irow_aitil2,&(a_flag[idof2]),&(j2),&ifree_aitil2,
			&nzs_aitil2,&contribution,&au_aitil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    insertas_ws(&irow_aitil2,&(a_flag[idof3]),&(j2),&ifree_aitil2,
			&nzs_aitil2,&contribution,&au_aitil2);
	    i=i+iadd;    
	  }else if(idof2>-1 && idof3>-1){
	    //2D auf 3D
	    if(jrow==2){
	      // k=idof2
	      iadd=0;
	      e1=0.0;
	      e2=au_ddtil2i[i];
	      if(i+1<jq_ddtil2i[j+1]-1 &&
		 a_flagr[irow_ddtil2i[i+1]-1]-1==idof3){
		e3=au_ddtil2i[i+1];++iadd;}else{e3=0.0;}
	    }else{
	      // k=idof3
	      iadd=0;
	      e1=0.0;
	      e2=0.0;
	      e3=au_ddtil2i[i];
	    }	  
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    insertas_ws(&irow_aitil2,&(a_flag[idof2]),&(j2),&ifree_aitil2,
			&nzs_aitil2,&contribution,&au_aitil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    insertas_ws(&irow_aitil2,&(a_flag[idof3]),&(j2),&ifree_aitil2,
			&nzs_aitil2,&contribution,&au_aitil2);
	    i=i+iadd;    	       	     
	  }else if(idof1>-1 && idof3>-1){
	    //2D auf 3D	
	    if(jrow==1){
	      // k=idof1
	      iadd=0;
	      e1=au_ddtil2i[i];
	      e2=0.0;
	      if(i+1<jq_ddtil2i[j+1]-1 &&
		 a_flagr[irow_ddtil2i[i+1]-1]-1==idof3){
		e3=au_ddtil2i[i+1];++iadd;}else{e3=0.0;}
	    }else{
	      // k=idof3
	      iadd=0;
	      e1=0.0;
	      e2=0.0;
	      e3=au_ddtil2i[i];
	    }
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    insertas_ws(&irow_aitil2,&(a_flag[idof1]),&(j2),&ifree_aitil2,
			&nzs_aitil2,&contribution,&au_aitil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    insertas_ws(&irow_aitil2,&(a_flag[idof3]),&(j2),&ifree_aitil2,
			&nzs_aitil2,&contribution,&au_aitil2);
	    i=i+iadd;       	     
	  }else if(idof1>-1 && idof2>-1){
	    //2D auf 3D	
	    if(jrow==1){
	      // k=idof1
	      iadd=0;
	      e1=au_ddtil2i[i];
	      if(i+1<jq_ddtil2i[j+1]-1 &&
		 a_flagr[irow_ddtil2i[i+1]-1]-1==idof2){
		e2=au_ddtil2i[i+1];++iadd;}else{e2=0.0;}
	      e3=0.0;
	    }else{
	      // k=idof3
	      iadd=0;
	      e1=0.0;
	      e2=au_ddtil2i[i];
	      e3=0.0;
	    }
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    insertas_ws(&irow_aitil2,&(a_flag[idof1]),&(j2),&ifree_aitil2,
			&nzs_aitil2,&contribution,&au_aitil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    insertas_ws(&irow_aitil2,&(a_flag[idof2]),&(j2),&ifree_aitil2,
			&nzs_aitil2,&contribution,&au_aitil2);
	    i=i+iadd;		           
	  }else {
	    // 1D auf 3D                   
	    e1=0.0;e2=0.0;e3=0.0;
	    if(idof1>-1){
	      e1=au_ddtil2i[i];
	      contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	      insertas_ws(&irow_aitil2,&(a_flag[idof1]),&(j2),&ifree_aitil2,
			  &nzs_aitil2,&contribution,&au_aitil2);	     
	    }else if(idof2>-1){
	      e2=au_ddtil2i[i];
	      contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	      insertas_ws(&irow_aitil2,&(a_flag[idof2]),&(j2),&ifree_aitil2,
			  &nzs_aitil2,&contribution,&au_aitil2);		
	    }else{
	      e3=au_ddtil2i[i];
	      contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	      insertas_ws(&irow_aitil2,&(a_flag[idof3]),&(j2),&ifree_aitil2,
			  &nzs_aitil2,&contribution,&au_aitil2);
	    }   
	  }
	}
      }
    }
    jq_aitil2[j+1]=ifree_aitil2;
  }
  
  nzs_aitil2=ifree_aitil2-1;
  add_rect(au_aitil1,irow_aitil1,jq_aitil1,*row_la,*row_li,
	   au_aitil2,irow_aitil2,jq_aitil2,*row_la,*row_li,
	   &au_aitil,&irow_aitil,jq_aitil,&nzs_aitil);
  
  SFREE(au_aitil1);SFREE(irow_aitil1);SFREE(jq_aitil1);
  SFREE(au_aitil2);SFREE(irow_aitil2);SFREE(jq_aitil2); 
  

  /* K_AA **/
  /* loop over columns **/
  
  nzs_aatil1=jq_daa[*row_la];
  NNEW(irow_aatil1,ITG,nzs_aatil1);
  NNEW(jq_aatil1,ITG,*row_la+1);
  NNEW(au_aatil1,double,nzs_aatil1);
  ifree_aatil1=1;
  jq_aatil1[0]=1;
  for(j=0;j<*row_la;j++){
    //loop over columns    
    j2=j+1;
    for(i=jq_daa[j]-1;i<jq_daa[j+1]-1;i++){
      //loop over rows	  
      k=irow_daa[i]-1;         
      islavnodeentry=floor(islavactdof[a_flagr[k]-1]/10.);	  
      jrow= islavactdof[a_flagr[k]-1]-10*islavnodeentry;	  
      node=islavnode[islavnodeentry-1];	   
      idof1=nactdof[mt*node-3]-1;	   
      idof2=nactdof[mt*node-2]-1;	   
      idof3=nactdof[mt*node-1]-1;
      
      /* get normal and tangetials **/
      
      for(l=0;l<3;l++){	     
	n[l]=slavnor[3*(islavnodeentry-1)+l];	     
	n2[l]=slavnor[3*(islavnodeentry-1)+l];	     	   
      }
      for(l=0;l<6;l++){	     
	t[l]=slavtan[6*(islavnodeentry-1)+l];	     
	that[l]=slavtan[6*(islavnodeentry-1)+l];	     	   
      }
      trafontspcmpc(n,t,n2,that,&islavnodeentry,nboun,ndirboun,nodeboun,xboun,
		    ipompc,nodempc,coefmpc,ikboun,ilboun,ikmpc,ilmpc,
		    nslavspc,islavspc,nslavmpc,islavmpc,
		    &node);
      
      /* calculate fields needed for Coulomb friction **/
      
      utildep_t[0]=u_tilde[(islavnodeentry-1)*3]*that[0]+
	u_tilde[(islavnodeentry-1)*3+1]*that[1]+
	u_tilde[(islavnodeentry-1)*3+2]*that[2];
      utildep_t[1]=u_tilde[(islavnodeentry-1)*3]*that[3]+
	u_tilde[(islavnodeentry-1)*3+1]*that[4]+
	u_tilde[(islavnodeentry-1)*3+2]*that[5];
      lambda_n=n2[0]*cstress2[(islavnodeentry-1)*mt]+
	n2[1]*cstress2[(islavnodeentry-1)*mt+1]+
	n2[2]*cstress2[(islavnodeentry-1)*mt+2];
      lambda_t[0]=that[0]*cstress2[(islavnodeentry-1)*mt]+
	that[1]*cstress2[(islavnodeentry-1)*mt+1]+
	that[2]*cstress2[(islavnodeentry-1)*mt+2];
      lambda_t[1]=that[3]*cstress2[(islavnodeentry-1)*mt]+
	that[4]*cstress2[(islavnodeentry-1)*mt+1]+
	that[5]*cstress2[(islavnodeentry-1)*mt+2];
      lambdaini_t[0]=that[0]*cstressini2[(islavnodeentry-1)*mt]+
	that[1]*cstressini2[(islavnodeentry-1)*mt+1]+
	that[2]*cstressini2[(islavnodeentry-1)*mt+2];
      lambdaini_t[1]=that[3]*cstressini2[(islavnodeentry-1)*mt]+
	that[4]*cstressini2[(islavnodeentry-1)*mt+1]+
	that[5]*cstressini2[(islavnodeentry-1)*mt+2];
      lambdatilde_t[0]=lambda_t[0]-lambdaini_t[0];
      lambdatilde_t[1]=lambda_t[1]-lambdaini_t[1];	   
      bp=bp_old[islavnodeentry-1];
      
      /* perturbed lagrange,normal direction 
	 (see phd-thesis Sitzmann,Chapter 3.4.1) */
      
      scal=Ddtil[jqdtil[node-1]-1];
      FORTRAN(getcontactparams,(&mu,&regmode,&aninvloc,&atauinvloc,
				&p0,&beta_e,tietol,elcon,
				&islavtie[islavnodeentry-1],ncmat_,
				ntmat_));
      
      constantt=min(constant,1.0/atauinvloc);
      derivmode=1;
      
      FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&dgnc,
				   &aninvloc,&p0,&beta_e,elcon,nelcon,
				   &islavtie[islavnodeentry-1],ntmat_,
				   plicon,nplicon,npmat_,ncmat_,tietol,&scal));
      
      if(islavact[islavnodeentry-1]==1 || islavact[islavnodeentry-1]==2){
	derivmode=1;
	  
	/* perturbed lagrange,tangential direction 
	   (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	FORTRAN(regularization_slip_lin,(utildep_t,&bp,&atauinvloc,resreg,
					 &derivmode,islavact,lambda_t,
					 lambdatilde_t,&constantt,
					 &islavnodeentry,n2,t,that,&mu,rslip,
					 ltslip,ltu));
	
	rphat[0]=0.0;
	rphat[1]=0.0;
	if(jrow==4 && ithermal[0]<2){ 
	  printf(" *ERROR in trafontmortar2\n");
	  // something went wrong
	  FORTRAN(stop,());
	}else{
	  // mechanical part
	  if(idof1>-1 && idof2>-1 && idof3>-1){	     	        
	    e1=au_daa[i];	        
	    e2=au_daa[i+1];	        
	    e3=au_daa[i+2];	
	    contribution=(dgnc)*(n[0]*e1+n[1]*e2+n[2]*e3);
	    insertas_ws(&irow_aatil1,&(irow_daa[i]),&j2,&ifree_aatil1,
			&nzs_aatil1,&contribution,&au_aatil1);
	    contribution=(rslip[0]*e1+rslip[1]*e2+rslip[2]*e3 );
	    insertas_ws(&irow_aatil1,&(irow_daa[i+1]),&j2,&ifree_aatil1,
			&nzs_aatil1,&contribution,&au_aatil1);
	    contribution=(rslip[3]*e1+rslip[4]*e2+rslip[5]*e3);
	    insertas_ws(&irow_aatil1,&(irow_daa[i+2]),&j2,&ifree_aatil1,
			&nzs_aatil1,&contribution,&au_aatil1);
	    i=i+2;	      
	  }else if(idof2>-1 && idof3>-1){
	    //2D auf 3D		              
	    e1=au_daa[i];	        
	    e2=au_daa[i+1];	    
	    t1=rslip[4];				
	    t2=rslip[5];
	    n11=n2[1];
	    n22=n2[2];
	    contribution=(dgnc)*(n11*e1+n22*e2);
	    insertas_ws(&irow_aatil1,&(irow_daa[i]),&j2,&ifree_aatil1,
			&nzs_aatil1,&contribution,&au_aatil1);		   
	    contribution=(t1*e1+t2*e2 );
	    insertas_ws(&irow_aatil1,&(irow_daa[i+1]),&j2,&ifree_aatil1,
			&nzs_aatil1,&contribution,&au_aatil1);
	    i=i+1;	     
	  }else if(idof1>-1 && idof3>-1){
	    //2D auf 3D		              
	    e1=au_daa[i];	        
	    e2=au_daa[i+1];	  		              
	    t1=rslip[3];				
	    t2=rslip[5];	
	    n11=n2[0];
	    n22=n2[2];		 
	    contribution=(dgnc)*(n11*e1+n22*e2);
	    insertas_ws(&irow_aatil1,&(irow_daa[i]),&j2,&ifree_aatil1,
			&nzs_aatil1,&contribution,&au_aatil1);		   
	    contribution=(t1*e1+t2*e2 );
	    insertas_ws(&irow_aatil1,&(irow_daa[i+1]),&j2,&ifree_aatil1,
			&nzs_aatil1,&contribution,&au_aatil1);
	    i=i+1;      	     	       
	  }else if(idof1>-1 && idof2>-1){
	    //2D auf 3D	
	    e1=au_daa[i];	        
	    e2=au_daa[i+1];	  		              
	    t1=rslip[3];				
	    t2=rslip[4];	
	    n11=n2[0];
	    n22=n2[1];	
	    contribution=(dgnc)*(n11*e1+n22*e2);
	    insertas_ws(&irow_aatil1,&(irow_daa[i]),&j2,&ifree_aatil1,
			&nzs_aatil1,&contribution,&au_aatil1);		   
	    contribution=(t1*e1+t2*e2 );
	    insertas_ws(&irow_aatil1,&(irow_daa[i+1]),&j2,&ifree_aatil1,
			&nzs_aatil1,&contribution,&au_aatil1);
	    i=i+1;  			       	       	     
	  }else {
	    // 1D auf 3D 
	    e1=au_daa[i];
	    if(idof1>-1){
	      n11=n2[0];
	    }else if(idof2>-1){
	      n11=n2[1];
	    }else{
	      n11=n2[2];
	    }
	    contribution=(dgnc)*(n11*e1);
	    insertas_ws(&irow_aatil1,&(irow_daa[i]),&j2,&ifree_aatil1,
			&nzs_aatil1,&contribution,&au_aatil1); 
	  } 
	}
      }  	  	
    }        
    jq_aatil1[j+1]=ifree_aatil1;
  }
  RENEW(irow_aatil1,ITG,ifree_aatil1-1);
  RENEW(au_aatil1,double,ifree_aatil1-1);

  /* add diagonal terms **/

  nzs_aatil2=3*(jq_ddtil2a[*row_la]-1);
  NNEW(au_aatil2,double,3*(jq_ddtil2a[*row_la]-1));
  NNEW(irow_aatil2,ITG,3*(jq_ddtil2a[*row_la]-1));
  NNEW(jq_aatil2,ITG,*row_la+1);
  ifree_aatil2=1;
  jq_aatil2[0]=1;
  for(j=0;j<*row_la;j++){
    //loop over columns  I   
    j2=j+1;
    for(i=jq_ddtil2a[j]-1;i<jq_ddtil2a[j+1]-1;i++){
      //loop over rows A  
      k=irow_ddtil2a[i]-1; 
      islavnodeentry=floor(islavactdof[a_flagr[k]-1]/10.);	  
      jrow= islavactdof[a_flagr[k]-1]-10*islavnodeentry; 
      node=islavnode[islavnodeentry-1];	   
      idof1=nactdof[mt*node-3]-1;	   
      idof2=nactdof[mt*node-2]-1;	   
      idof3=nactdof[mt*node-1]-1;
      for(l=0;l<3;l++){	     
	n[l]=slavnor[3*(islavnodeentry-1)+l];	     
	n2[l]=slavnor[3*(islavnodeentry-1)+l];	     	   
      }
      for(l=0;l<6;l++){	     
	t[l]=slavtan[6*(islavnodeentry-1)+l];	     
	that[l]=slavtan[6*(islavnodeentry-1)+l];	     	   
      }
      trafontspcmpc(n,t,n2,that,&islavnodeentry,
		    nboun,ndirboun,nodeboun,xboun,
		    ipompc,nodempc,coefmpc,
		    ikboun,ilboun,ikmpc,ilmpc,
		    nslavspc,islavspc,nslavmpc,islavmpc,
		    &node);
      
      /* calculate fields needed for Coulomb friction */
      
      utildep_t[0]=u_tilde[(islavnodeentry-1)*3]*that[0]+
	u_tilde[(islavnodeentry-1)*3+1]*that[1]+
	u_tilde[(islavnodeentry-1)*3+2]*that[2];
      utildep_t[1]=u_tilde[(islavnodeentry-1)*3]*that[3]+
	u_tilde[(islavnodeentry-1)*3+1]*that[4]+
	u_tilde[(islavnodeentry-1)*3+2]*that[5];
      lambda_n=n2[0]*cstress2[(islavnodeentry-1)*mt]+
	n2[1]*cstress2[(islavnodeentry-1)*mt+1]+
	n2[2]*cstress2[(islavnodeentry-1)*mt+2];
      lambda_t[0]=that[0]*cstress2[(islavnodeentry-1)*mt]+
	that[1]*cstress2[(islavnodeentry-1)*mt+1]+
	that[2]*cstress2[(islavnodeentry-1)*mt+2];
      lambda_t[1]=that[3]*cstress2[(islavnodeentry-1)*mt]+
	that[4]*cstress2[(islavnodeentry-1)*mt+1]+
	that[5]*cstress2[(islavnodeentry-1)*mt+2];
      lambdaini_t[0]=that[0]*cstressini2[(islavnodeentry-1)*mt]+
	that[1]*cstressini2[(islavnodeentry-1)*mt+1]+
	that[2]*cstressini2[(islavnodeentry-1)*mt+2];
      lambdaini_t[1]=that[3]*cstressini2[(islavnodeentry-1)*mt]+
	that[4]*cstressini2[(islavnodeentry-1)*mt+1]+
	that[5]*cstressini2[(islavnodeentry-1)*mt+2];
      lambdatilde_t[0]=lambda_t[0]-lambdaini_t[0];
      lambdatilde_t[1]=lambda_t[1]-lambdaini_t[1];	   
      bp=bp_old[islavnodeentry-1];
      
      /* perturbed lagrange,normal direction 
	 (see phd-thesis Sitzmann,Chapter 3.4.1) */
      
      scal=Ddtil[jqdtil[node-1]-1];
      FORTRAN(getcontactparams,(&mu,&regmode,&aninvloc,&atauinvloc,
				&p0,&beta_e,tietol,elcon,
				&islavtie[islavnodeentry-1],ncmat_,
				ntmat_));
      
      constantt=min(constant,1.0/atauinvloc);
      derivmode=1;
      
      FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&dgnc,
				   &aninvloc,&p0,&beta_e,elcon,nelcon,
				   &islavtie[islavnodeentry-1],ntmat_,
				   plicon,nplicon,npmat_,ncmat_,tietol,&scal));
      
      if(islavact[islavnodeentry-1]==1 || islavact[islavnodeentry-1]==2){ 	    
	derivmode=1;
	  
	/* perturbed lagrange,tangential direction 
	   (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	FORTRAN(regularization_slip_lin,(utildep_t,&bp,&atauinvloc,resreg,
					 &derivmode,islavact,lambda_t,
					 lambdatilde_t,&constantt,
					 &islavnodeentry,n2,t,that,&mu,rslip,
					 ltslip,ltu));
	  
	rphat[0]=0.0;
	rphat[1]=0.0;
	if(jrow==4 && ithermal[0]<2){ 
	  printf(" *ERROR in trafontmortar2\n");
	  // something went wrong
	  FORTRAN(stop,());
	}else{
	  // mechanical part
	  if(idof1>-1 && idof2>-1 && idof3>-1){
	    // 3D	
	    if(jrow==1){
	      // k=idof1
	      iadd=0;
	      e1=au_ddtil2a[i];
	      if(i+1<jq_ddtil2a[j+1]-1 &&
		 a_flagr[irow_ddtil2a[i+1]-1]-1==idof2){
		e2=au_ddtil2a[i+1];++iadd;}else{e2=0.0;}
	      if(i+1<jq_ddtil2a[j+1]-1 &&
		 a_flagr[irow_ddtil2a[i+1]-1]-1==idof3){
		e3=au_ddtil2a[i+1];++iadd;
	      }else if(i+2<jq_ddtil2a[j+1]-1 &&
		       a_flagr[irow_ddtil2a[i+2]-1]-1==idof3){
		e3=au_ddtil2a[i+2];++iadd;
	      }else{e3=0.0;}
	    }else if(jrow==2){
	      // k=idof2
	      iadd=0;
	      e1=0.0;
	      e2=au_ddtil2a[i];
	      if(i+1<jq_ddtil2a[j+1]-1 &&
		 a_flagr[irow_ddtil2a[i+1]-1]-1==idof3){
		e3=au_ddtil2a[i+1];++iadd;}else{e3=0.0;}
	    }else{
	      // k=idof3
	      iadd=0;
	      e1=0.0;
	      e2=0.0;
	      e3=au_ddtil2a[i];
	    }
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    insertas_ws(&irow_aatil2,&(a_flag[idof1]),&(j2),&ifree_aatil2,
			&nzs_aatil2,&contribution,&au_aatil2);
	    contribution=(ltslip[0]*e1+ltslip[1]*e2+ltslip[2]*e3);
	    insertas_ws(&irow_aatil2,&(a_flag[idof2]),&(j2),&ifree_aatil2,
			&nzs_aatil2,&contribution,&au_aatil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    insertas_ws(&irow_aatil2,&(a_flag[idof3]),&(j2),&ifree_aatil2,
			&nzs_aatil2,&contribution,&au_aatil2);
	    i=i+iadd;    
	  }else if(idof2>-1 && idof3>-1){
	    //2D auf 3D
	    if(jrow==2){
	      // k=idof2
	      iadd=0;
	      e1=0.0;
	      e2=au_ddtil2a[i];
	      if(i+1<jq_ddtil2a[j+1]-1 &&
		 a_flagr[irow_ddtil2a[i+1]-1]-1==idof3){
		e3=au_ddtil2a[i+1];++iadd;}else{e3=0.0;}
	    }else{
	      // k=idof3
	      iadd=0;
	      e1=0.0;
	      e2=0.0;
	      e3=au_ddtil2a[i];
	    }	  
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    insertas_ws(&irow_aatil2,&(a_flag[idof2]),&(j2),&ifree_aatil2,
			&nzs_aatil2,&contribution,&au_aatil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    insertas_ws(&irow_aatil2,&(a_flag[idof3]),&(j2),&ifree_aatil2,
			&nzs_aatil2,&contribution,&au_aatil2);
	    i=i+iadd;    	       	     
	  }else if(idof1>-1 && idof3>-1){
	    //2D auf 3D	
	    if(jrow==1){
	      // k=idof1
	      iadd=0;
	      e1=au_ddtil2a[i];
	      e2=0.0;
	      if(i+1<jq_ddtil2a[j+1]-1 &&
		 a_flagr[irow_ddtil2a[i+1]-1]-1==idof3){
		e3=au_ddtil2a[i+1];++iadd;}else{e3=0.0;}
	    }else{
	      // k=idof3
	      iadd=0;
	      e1=0.0;
	      e2=0.0;
	      e3=au_ddtil2a[i];
	    }
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    insertas_ws(&irow_aatil2,&(a_flag[idof1]),&(j2),&ifree_aatil2,
			&nzs_aatil2,&contribution,&au_aatil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    insertas_ws(&irow_aatil2,&(a_flag[idof3]),&(j2),&ifree_aatil2,
			&nzs_aatil2,&contribution,&au_aatil2);
	    i=i+iadd;       	     
	  }else if(idof1>-1 && idof2>-1){
	    //2D auf 3D	
	    if(jrow==1){
	      // k=idof1
	      iadd=0;
	      e1=au_ddtil2a[i];
	      if(i+1<jq_ddtil2a[j+1]-1 &&
		 a_flagr[irow_ddtil2a[i+1]-1]-1==idof2){
		e2=au_ddtil2a[i+1];++iadd;}else{e2=0.0;}
	      e3=0.0;
	    }else{
	      // k=idof3
	      iadd=0;
	      e1=0.0;
	      e2=au_ddtil2a[i];
	      e3=0.0;
	    }
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    insertas_ws(&irow_aatil2,&(a_flag[idof1]),&(j2),&ifree_aatil2,
			&nzs_aatil2,&contribution,&au_aatil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    insertas_ws(&irow_aatil2,&(a_flag[idof2]),&(j2),&ifree_aatil2,
			&nzs_aatil2,&contribution,&au_aatil2);
	    i=i+iadd;		           
	  }else {
	    // 1D auf 3D                   
	    e1=0.0;e2=0.0;e3=0.0;
	    if(idof1>-1){
	      e1=au_ddtil2a[i];
	      contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	      insertas_ws(&irow_aatil2,&(a_flag[idof1]),&(j2),&ifree_aatil2,
			  &nzs_aatil2,&contribution,&au_aatil2);	     
	    }else if(idof2>-1){
	      e2=au_ddtil2a[i];
	      contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	      insertas_ws(&irow_aatil2,&(a_flag[idof2]),&(j2),&ifree_aatil2,
			  &nzs_aatil2,&contribution,&au_aatil2);		
	    }else{
	      e3=au_ddtil2a[i];
	      contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	      insertas_ws(&irow_aatil2,&(a_flag[idof3]),&(j2),&ifree_aatil2,
			  &nzs_aatil2,&contribution,&au_aatil2);
	    }   
	  }
	}
      }
    }
    jq_aatil2[j+1]=ifree_aatil2;
  }
  
  nzs_aatil2=ifree_aatil2-1;
  add_rect(au_aatil1,irow_aatil1,jq_aatil1,*row_la,*row_la,
	   au_aatil2,irow_aatil2,jq_aatil2,*row_la,*row_la,
	   &au_aatil,&irow_aatil,jq_aatil,&nzs_aatil);
  
  SFREE(au_aatil1);SFREE(irow_aatil1);SFREE(jq_aatil1);
  SFREE(au_aatil2);SFREE(irow_aatil2);SFREE(jq_aatil2); 
  
  /* changing f_da due to N and T (normal and tangential
     direction at the slave surface **/

  for(k=0;k<*row_la;k++){      
    if(islavactdof[a_flagr[k]-1]>0){	
      islavnodeentry=floor(islavactdof[a_flagr[k]-1]/10.);
      jrow= islavactdof[a_flagr[k]-1]-10*islavnodeentry;
      node=islavnode[islavnodeentry-1];
      idof1=nactdof[mt*node-3]-1;
      idof2=nactdof[mt*node-2]-1;
      idof3=nactdof[mt*node-1]-1;
      
      /* get normal and tangetials **/
      
      for(l=0;l<3;l++){	     
	n[l]=slavnor[3*(islavnodeentry-1)+l];	     
	n2[l]=slavnor[3*(islavnodeentry-1)+l];	     
      }
      for(l=0;l<6;l++){	     
	t[l]=slavtan[6*(islavnodeentry-1)+l];	     
	that[l]=slavtan[6*(islavnodeentry-1)+l];	     
      }	
     
      trafontspcmpc(n,t,n2,that,&islavnodeentry,nboun,ndirboun,nodeboun,xboun,
		    ipompc,nodempc,coefmpc,ikboun,ilboun,ikmpc,ilmpc,
		    nslavspc,islavspc,nslavmpc,islavmpc,
		    &node);
      
      /* calculate needed fields for coulomb friction **/
      
      utildep_t[0]=u_tilde[(islavnodeentry-1)*3]*that[0]+
	u_tilde[(islavnodeentry-1)*3+1]*that[1]+
	u_tilde[(islavnodeentry-1)*3+2]*that[2];
      utildep_t[1]=u_tilde[(islavnodeentry-1)*3]*that[3]+
	u_tilde[(islavnodeentry-1)*3+1]*that[4]+
	u_tilde[(islavnodeentry-1)*3+2]*that[5];
      lambda_n=n2[0]*cstress2[(islavnodeentry-1)*mt]+
	n2[1]*cstress2[(islavnodeentry-1)*mt+1]+
	n2[2]*cstress2[(islavnodeentry-1)*mt+2];
      lambda_t[0]=that[0]*cstress2[(islavnodeentry-1)*mt]+
	that[1]*cstress2[(islavnodeentry-1)*mt+1]+
	that[2]*cstress2[(islavnodeentry-1)*mt+2];
      lambda_t[1]=that[3]*cstress2[(islavnodeentry-1)*mt]+
	that[4]*cstress2[(islavnodeentry-1)*mt+1]+
	that[5]*cstress2[(islavnodeentry-1)*mt+2];	   
      lambdaini_t[0]=that[0]*cstressini2[(islavnodeentry-1)*mt]+
	that[1]*cstressini2[(islavnodeentry-1)*mt+1]+
	that[2]*cstressini2[(islavnodeentry-1)*mt+2];	   
      lambdaini_t[1]=that[3]*cstressini2[(islavnodeentry-1)*mt]+
	that[4]*cstressini2[(islavnodeentry-1)*mt+1]+
	that[5]*cstressini2[(islavnodeentry-1)*mt+2];	   
      lambdatilde_t[0]=lambda_t[0]-lambdaini_t[0];	   
      lambdatilde_t[1]=lambda_t[1]-lambdaini_t[1];		
      bp=bp_old[islavnodeentry-1];
      scal=Ddtil[jqdtil[node-1]-1];
      
      FORTRAN(getcontactparams,(&mu,&regmode,&aninvloc,&atauinvloc,
				&p0,&beta_e,tietol,elcon,
				&islavtie[islavnodeentry-1],ncmat_,
				ntmat_));
      
      constantt=min(constant,1.0/atauinvloc);
      
      /* perturbed lagrange,normal direction 
	 (see phd-thesis Sitzmann,Chapter 3.4.1) */
      
      derivmode=1;
      FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&dgnc,
				   &aninvloc,&p0,&beta_e,elcon,nelcon,
				   &islavtie[islavnodeentry-1],ntmat_,
				   plicon,nplicon,npmat_,ncmat_,tietol,&scal));
      derivmode=0;
      
      FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&gnc,
				   &aninvloc,&p0,&beta_e,elcon,nelcon,
				   &islavtie[islavnodeentry-1],ntmat_,
				   plicon,nplicon,npmat_,ncmat_,tietol,&scal));

      if(islavact[islavnodeentry-1]==1 || islavact[islavnodeentry-1]==2){
	derivmode=1;
	  
	/* perturbed lagrange,tangential direction 
	   (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	FORTRAN(regularization_slip_lin,(utildep_t,&bp,&atauinvloc,resreg,
					 &derivmode,islavact,lambda_t,
					 lambdatilde_t,&constantt,
					 &islavnodeentry,n2,t,that,&mu,rslip,
					 ltslip,ltu));
	rphat[0]=0.0;
	rphat[1]=0.0;
	hpn=gap[islavnodeentry-1]+gnc-dgnc*lambda_n;
	dgnc1=dgnc;
	
	if(jrow==4 && ithermal[0]<2){ 
	  printf(" *ERROR in trafontmortar2\n");
	  // something went wrong
	  FORTRAN(stop,());
	}else{
	  // mechanical part
	  if(idof1>-1 && idof2>-1 && idof3>-1 &&jrow==1){
	    e1=f_da[k];
	    e2=f_da[k+1];
	    e3=f_da[k+2];
	    /* right side if solving K du=f*/ 
	    f_atil[k]=hpn+(dgnc1)*(n[0]*e1+n[1]*e2+n[2]*e3);
	    f_atil[k+1]=(rslip[0]*e1+rslip[1]*e2+rslip[2]*e3)+rphat[0]-ltu[0];
	    f_atil[k+2]=(rslip[3]*e1+rslip[4]*e2+rslip[5]*e3)+rphat[1]-ltu[1];
	    k=k+2;
	  }else if(idof2>-1 && idof3>-1){
	    //2D auf 3D
	    e1=f_da[k];
	    e2=f_da[k+1];
	    t1=rslip[4];
	    t2=rslip[5];
	    f_atil[k]=hpn+(dgnc1)*(n2[1]*e1+n2[2]*e2);
	    f_atil[k+1]=(t1*e1+t2*e2)+rphat[1]-ltu[1];
	    k=k+1;  
	  }else if(idof1>-1 && idof3>-1){
	    //2D auf 3D
	    e1=f_da[k];
	    e2=f_da[k+1];
	    t1=rslip[3];
	    t2=rslip[5];
	    f_atil[k]=hpn+(dgnc1)*(n2[0]*e1+n2[2]*e2);
	    f_atil[k+1]=(t1*e1+t2*e2)+rphat[1]-ltu[1];
	    k=k+1; 
	  }else if(idof1>-1 && idof2>-1){
	    //2D auf 3D
	    e1=f_da[k];
	    e2=f_da[k+1];
	    t1=rslip[3];
	    t2=rslip[4];
	    f_atil[k]=hpn+(dgnc1)*(n2[0]*e1+n2[1]*e2);
	    f_atil[k+1]=(t1*e1+t2*e2)+rphat[1]-ltu[1];
	    k=k+1; 				 
	  }else{
	    //1D auf 3D 
	    e2=f_da[k];
	    if(idof1>-1){
	      n11=n2[0];
	      f_atil[k]=hpn+(dgnc1)*(n11*e2);
	    }else if(idof2>-1) {
	      n11=n2[1];
	      f_atil[k]=hpn+(dgnc1)*(n11*e2);	    
	    }else{
	      n11=n2[2];
	      f_atil[k]=hpn+(dgnc1)*(n11*e2);	      
	    }
	  }
	}
      }
      
    }
  }
  
  SFREE(u_tilde);
  SFREE(cstress2);
  SFREE(cstressini2);
  
  *au_antilp=au_antil;*au_amtilp=au_amtil;*au_aitilp=au_aitil;
  *au_aatilp=au_aatil;
  *irow_antilp=irow_antil;*irow_amtilp=irow_amtil;*irow_aitilp=irow_aitil;
  *irow_aatilp=irow_aatil;
  
}
