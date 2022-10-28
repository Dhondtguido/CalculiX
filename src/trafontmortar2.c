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

/**
 *  Condense Lagrange Multiplier and embedd contact conditionds for \f$ K_{AX}\f$
 * changing au due to N and T (normal and tangential
 *    direction at the slave surface) 
 * 	changing b due to N and T (normal and tangential
 *	direction at the slave surface) 
 *  phd-thesis Sitzmann, Chapter 3+4, equation (4.15) (Tpye=MORTAR/LINMORTAR) or (4.24) (TYPE=PGLINMORTAR) 
 * 
 * author: Saskia Sitzmann
 *  [in] islavactdof      (i)=10*slavenodenumber+direction for active dof i
 *  [in] islavact		(i) indicates, if slave node i is active (=-3 no-slave-node, =-2 no-LM-node, =-1 no-gap-node, =0 inactive node, =1 sticky node, =2 slipping/active node)  
 *  [in] nslavnode	(i)pointer into field isalvnode for contact tie i 
 *  [in] nmastnode	(i)pointer into field imastnode for contact tie i 
 *  [in] f_da		\f$ r_A \f$ residual for active slave nodes
 *  [out] f_atil		\f$ r_A \f$ condensed residual for active slave nodes
 *  [in] au_dan		\f$ K_{AN}\f$	
 *  [in] irow_dan		rows of matrix au_dan
 *  [in] jq_dan		column pointer to irow_dan
 *  [in] au_dam		\f$ K_{AM}\f$
 *  [in] irow_dam		rows of matrix au_dam
 *  [in] jq_dam		column pointer to irow_dam
 *  [in] au_dai		\f$ K_{AI}\f$
 *  [in] irow_dai		rows of matrix au_dai
 *  [in] jq_dai		column pointer to irow_dai
 *  [in] au_daa		\f$ K_{AA}\f$
 *  [in] irow_daa		rows of matrix au_daa
 *  [in] jq_daa		column pointer to irow_daa
 *  [out] au_antilp	condensed \f$ K_{AN}\f$
 *  [out] irow_antilp	rows of matrix au_antil
 *  [out] jq_antil	column pointer to irow_antil
 *  [out] au_amtilp	condensed \f$ K_{AM}\f$
 *  [out] irow_amtilp	rows of matrix au_amtil	
 *  [out] jq_amtil	column pointer to irow_amtil
 *  [out] au_aitilp	condensed \f$ K_{AI}\f$
 *  [out] irow_aitilp	rows of matrix au_aitil
 *  [out] jq_aitil	column pointer to irow_aitil
 *  [out] au_aatilp	condensed \f$ K_{AA}\f$
 *  [out] irow_aatilp	rows of matrix au_aatil
 *  [out] jq_aatil	column pointer to irow_aatil
 *  [in] gap		(i) gap for node i on slave surface 
 *  [in] Bd		coupling matrix \f$ B_d[p,q]=\int \psi_p \phi_q dS \f$, \f$ p \in S, q \in M \f$ 
 *  [in] irowb		field containing row numbers of Bd
 *  [in] jqb		pointer into field irowb
 *  [out] Dd		coupling matrix \f$ D_d[p,q]=\int \psi_p \phi_q dS \f$, \f$ p,q \in S \f$ 
 *  [out] irowd		field containing row numbers of Dd
 *  [out] jqd		pointer into field irowd 
 *  [out] Ddtil		coupling matrix \f$ \tilde{D}_d[p,q]=\int \psi_p \tilde{\phi}_q dS \f$, \f$ p,q \in S \f$ 
 *  [out] irowdtil	field containing row numbers of Ddtil
 *  [out] jqdtil		pointer into field irowdtil  
 *  [in] au_bdtil2	\f$ B_d|_{AM} \f$ for active degrees of freedom
 *  [in] irow_bdtil2	rows of matrix au_bdtil2	
 *  [in] jq_bdtil2	pointer into field irow_bdtil2
 *  [in] au_ddtil2i	\f$ D_d|_{AI} \f$ for active degrees of freedom
 *  [in] irow_ddtil2i	rows of matrix au_bdtil2i
 *  [in] jq_ddtil2i 	pointer into field irow_bdtil2i
 *  [in] au_ddtil2a	\f$ D_d|_{AA} \f$ for active degrees of freedom
 *  [in] irow_ddtil2a	rows of matrix au_bdtil2a
 *  [in] jq_ddtil2a 	pointer into field irow_bdtil2a
 *  [in] m_flagr		field from local to global dof for master nodes
 *  [in] i_flagr		field from local to global dof for inactive nodes
 *  [in] a_flagr		field from local to global dof for active nodes
 *  [in] a_flag		field from global to local dof for active nodes
 *  [in] i_flag		field from global to local dof for inactive nodes
 *  [in] m_flag		field from global to local dof for master nodes
 *  [in] row_ln		number of \f$ N \f$ rows
 *  [in] row_lm		number of \f$ M \f$ rows (Master)
 *  [in] row_li		number of \f$ I \f$ rows (Inactive)
 *  [in] row_la		number of \f$ A \f$ rows (Slave)
 *  [in] slavnor		slave normal
 *  [in] slavtan		slave tangent 
 *  [in] cstress		current Lagrange multiplier 
 *  [in] cstressini	Lagrange multiplier at start of the increment 
 *  [in] bp_old		old friction bounds
 *  [in] nactdof 		(i,j) actual degree of freedom for direction i of node j 
 *  [in] islavnode	field storing the nodes of the slave surface 
 *  [in] imastnode	field storing the nodes of the master surfaces 
 *  [in] nslavspc		(2*i) pointer to islavspc...
 *  [in] islavspc         ... which stores SPCs for slave node i
 *  [in] nsspc            number of SPC for slave nodes
 *  [in] nslavmpc		(2*i) pointer to islavmpc...
 *  [in] islavmpc		... which stores MPCs for slave node i
 *  [in] nsmpc		number of MPC for slave nodes
 *  [in] nmastspc		(2*i) pointer to imastspc...
 *  [in] imastspc         ... which stores SPCs for master node i
 *  [in] nmspc            number of SPC for master nodes
 *  [in] nmastmpc		(2*i) pointer to imastmpc...
 *  [in] imastmpc		... which stores MPCs for master node i
 *  [in] nmmpc		number of MPC for master nodes 
 *  [in] islavactdoftie   (i)=tie number for active dof i
 *  [in] irowtloc		field containing row numbers of autloc
 *  [in] jqtloc	        pointer into field irowtloc
 *  [in] autloc		transformation matrix \f$ T[p,q]\f$ for slave nodes \f$ p,q \f$ 
 *  [in] irowtlocinv	field containing row numbers of autlocinv
 *  [in] jqtlocinv	pointer into field irowtlocinv
 *  [in] autlocinv	transformation matrix \f$ T^{-1}[p,q]\f$ for slave nodes \f$ p,q \f$  
 *  [in] islavnodeinv     (i) slave node index for node i
 *  [in] lambdaiwan       Lagrange multiplier splitted to Iwan elements
 *  [in] lambdaiwanini    Lagrange multiplier splitted to Iwan elements at start of increment 
 */
void trafontmortar2(ITG *neq,ITG *nzs,ITG *islavactdof,ITG *islavact,
		    ITG *nslavnode,ITG *nmastnode,double *f_da,double *f_atil,
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
		    ITG *imastnode,ITG *ntie,ITG *mi,ITG *nk,
		    ITG *nboun,ITG *ndirboun,ITG *nodeboun,double *xboun,
		    ITG *nmpc,ITG *ipompc,ITG *nodempc,double *coefmpc,
		    ITG *ikboun,ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
		    ITG *nslavspc,ITG *islavspc,ITG *nsspc,ITG *nslavmpc,
		    ITG *islavmpc,ITG *nsmpc,
		    ITG *nmastspc,ITG *imastspc,ITG *nmspc,ITG *nmastmpc,
		    ITG *imastmpc,ITG *nmmpc,
		    char *tieset,
		    ITG *islavactdoftie,ITG *nelcon,double  *elcon,
		    double *tietol,ITG *ncmat_,ITG *ntmat_,
		    double *plicon,ITG *nplicon,ITG *npmat_,double *dtime,
		    ITG *irowtloc,ITG *jqtloc,double *autloc, 
		    ITG *irowtlocinv,ITG *jqtlocinv,double *autlocinv,
		    ITG *islavnodeinv,double *lambdaiwan,double *lambdaiwanini,
		    ITG *iit,ITG *nmethod,double *beta,ITG *ithermal,
		    double *plkcon,ITG *nplkcon){
  
  ITG i,j,jj,j2,k,l,debug,idof1,idof2,idof3,iadd,iwan,jrow,jcol,islavnodeentry,
    jslavnodeentry,mt=mi[1]+1,nodes,node,derivmode,regmode,regmodet,yielded,
    nodem,*irow_antil=NULL,*irow_amtil=NULL,*irow_aitil=NULL,*irow_aatil=NULL,
    *irow_amtil1=NULL,*irow_amtil2=NULL,*irow_aitil1=NULL,*irow_aitil2=NULL,
    *irow_aatil1=NULL,*irow_aatil2=NULL,nzs_antil,nzs_amtil,nzs_aitil,
    nzs_aatil,*jq_amtil1=NULL,*jq_amtil2=NULL,nzs_amtil1,nzs_amtil2,nzs_aitil1,
    nzs_aitil2,nzs_aatil1,nzs_aatil2,*jq_aitil1=NULL,*jq_aitil2=NULL,
    *jq_aatil1=NULL,*jq_aatil2=NULL,ifree_antil,ifree_amtil1,ifree_amtil2,
    ifree_aitil1,ifree_aitil2,ifree_aatil1,ifree_aatil2;
  
  double t1,t2,e1,e2,e3,contribution,dut[2],hpn,scal,bp,up_n,constant=1.E10,
    constantt=1.E10,atau,lambda_n,nlambda_t,*u_tilde=NULL,resreg[2],
    *cstress2=NULL,*cstressini2=NULL,that[6],n11,n22,aninvloc,gnc,dgnc,dgnc1,
    mu,p0,beta_e,atauinvloc,lambda_t[2],lambdaini_t[2],lambdatilde_t[2],ltu[2],
    ltslip[6],rphat[2],n[3],n2[3],t[6],utildep_t[2],rslip[6],*au_antil=NULL,
    *au_amtil=NULL,*au_aitil=NULL,*au_aatil=NULL,*au_amtil1=NULL,
    *au_amtil2=NULL,*au_aitil1=NULL,*au_aitil2=NULL,*au_aatil1=NULL,
    *au_aatil2=NULL;
  
  debug=0;
  
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
	if(debug==1) printf("j %" ITGFORMAT " node %" ITGFORMAT
			    " vold %e %e %e\n",j,nodes,vold[mt*(nodes)-3],
			    vold[mt*(nodes)-2],vold[mt*(nodes)-1] );
	if(debug==1) printf("j %" ITGFORMAT " node %" ITGFORMAT
			    " vini %e %e %e\n",j,nodes,vini[mt*(nodes)-3],
			    vini[mt*(nodes)-2],vini[mt*(nodes)-1] );
	if(debug==1) printf("j %" ITGFORMAT " node %" ITGFORMAT
			    " uold %e %e %e\n",j,nodes,u_tilde[j*3],
			    u_tilde[j*3+1],u_tilde[j*3+2] );
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
	if(debug==1) printf("j %" ITGFORMAT " node %" ITGFORMAT
			    " vold %e %e %e\n",j,nodes,vold[mt*(nodes)-3],
			    vold[mt*(nodes)-2],vold[mt*(nodes)-1] );
	if(debug==1) printf("j %" ITGFORMAT " node %" ITGFORMAT
			    " vini %e %e %e\n",j,nodes,vini[mt*(nodes)-3],
			    vini[mt*(nodes)-2],vini[mt*(nodes)-1] );
	if(debug==1) printf("j %" ITGFORMAT " node %" ITGFORMAT
			    " uold %e %e %e\n",j,nodes,u_tilde[j*3],
			    u_tilde[j*3+1],u_tilde[j*3+2] );
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
		    nmpc,ipompc,nodempc,coefmpc,ikboun,ilboun,ikmpc,ilmpc,
		    nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
		    nmastspc,imastspc,nmspc,nmastmpc,imastmpc,nmmpc,
		    &debug,&node);
      
      /* calculate fields needed for Coulomb friction */
      
      up_n=   u_tilde[(islavnodeentry-1)*3]*n2[0]+
	u_tilde[(islavnodeentry-1)*3+1]*n2[1]+
	u_tilde[(islavnodeentry-1)*3+2]*n2[2];
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
      nlambda_t=sqrt(lambda_t[0]*lambda_t[0]+lambda_t[1]*lambda_t[1]);
      bp=bp_old[islavnodeentry-1];
      
      /* perturbed lagrange,normal direction
	 (see phd-thesis Sitzmann,Chapter 3.4.1) */
      
      scal=Ddtil[jqdtil[node-1]-1];
      FORTRAN(getcontactparams,(&mu,&regmode,&regmodet,&aninvloc,&atauinvloc,
				&p0,&beta_e,tietol,elcon,
				&islavactdoftie[islavnodeentry-1],ncmat_,
				ntmat_,&iwan));
      constantt=min(constant,1.0/atauinvloc);
      derivmode=1;
      FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&dgnc,
				   &aninvloc,&p0,&beta_e,elcon,nelcon,
				   &islavactdoftie[islavnodeentry-1],ntmat_,
				   plicon,nplicon,npmat_,ncmat_,tietol,&scal));
      
      if(islavact[islavnodeentry-1]==1 || islavact[islavnodeentry-1]==2){
	if(regmodet==1){
	  derivmode=1;
	  
	  /* perturbed lagrange,tangential direction
	     (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	  FORTRAN(regularization_slip_lin,(utildep_t,&bp,&atauinvloc,resreg,
					   &derivmode,islavact,lambda_t,
					   lambdatilde_t,&constantt,&debug,
					   &islavnodeentry,n2,t,that,&mu,rslip,
					   ltslip,ltu));
          rphat[0]=0.0;
	  rphat[1]=0.0;
	}else{
	  derivmode=1;
	  atau=1.0/atauinvloc;
	  
	  /* perturbed lagrange,tangential direction with Iwan elements  
	     (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	  FORTRAN(regularization_slip_iwan,(&lambda_n,utildep_t,&bp,&atau,
					    resreg,&derivmode,&regmodet,
					    lambdaiwan,lambdaiwanini,
					    &islavnodeentry,n,t,&mu,rslip,
					    ltslip,ltu,&yielded,iit,&debug,
					    &iwan,dut));
	  rphat[0]=0.0;
	  rphat[1]=0.0;
	}
	
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
	    if(debug==1){printf("\t au_dan %e %e %e \n \t au_antil",e1,e2,e3);}

	    contribution=(dgnc)*(n[0]*e1+n[1]*e2+n[2]*e3);
	    insertas_ws(&irow_antil,&(irow_dan[i]),&j2,&ifree_antil,
			&nzs_antil,&contribution,&au_antil);
	    if(debug==1){printf(" %e ",contribution/Ddtil[jqdtil[node-1]-1]);}
	    
	    contribution=(rslip[0]*e1+rslip[1]*e2+rslip[2]*e3);
	    insertas_ws(&irow_antil,&(irow_dan[i+1]),&j2,&ifree_antil,
			&nzs_antil,&contribution,&au_antil);
	    if(debug==1){printf(" %e ",contribution/Ddtil[jqdtil[node-1]-1]);}

	    contribution=(rslip[3]*e1+rslip[4]*e2+rslip[5]*e3);
	    insertas_ws(&irow_antil,&(irow_dan[i+2]),&j2,&ifree_antil,
			&nzs_antil,&contribution,&au_antil);
	    if(debug==1){printf(" %e \n",contribution/Ddtil[jqdtil[node-1]-1]);}
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
  if(debug==1)printf("\ttrafoNT2: au_dan %" ITGFORMAT " au_antil %" ITGFORMAT "\n",
	 jq_dan[*row_ln]-1,jq_antil[*row_ln]-1);
  
  /* K_AM^til **/
  
  debug=0;
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
      if(islavactdof[m_flagr[j]-1]<0){
	jslavnodeentry=floor(-islavactdof[m_flagr[j]-1]/10.);	  
	jcol= -islavactdof[m_flagr[j]-1]-10*jslavnodeentry;
      }else{
	jslavnodeentry=floor(islavactdof[m_flagr[j]-1]/10.);	  
	jcol= islavactdof[m_flagr[j]-1]-10*jslavnodeentry;	
      }
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
		    nmpc,ipompc,nodempc,coefmpc,
		    ikboun,ilboun,ikmpc,ilmpc,
		    nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
		    nmastspc,imastspc,nmspc,nmastmpc,imastmpc,nmmpc,
		    &debug,&node);
      
      /* calculate fields needed for Coulomb friction **/
      
      up_n=   u_tilde[(islavnodeentry-1)*3]*n2[0]+
	u_tilde[(islavnodeentry-1)*3+1]*n2[1]+
	u_tilde[(islavnodeentry-1)*3+2]*n2[2];
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
      nlambda_t=sqrt(lambda_t[0]*lambda_t[0]+lambda_t[1]*lambda_t[1]);
      bp=bp_old[islavnodeentry-1];
      
      /* perturbed lagrange,normal direction  
	 (see phd-thesis Sitzmann,Chapter 3.4.1) */
      
      scal=Ddtil[jqdtil[node-1]-1];
      FORTRAN(getcontactparams,(&mu,&regmode,&regmodet,&aninvloc,&atauinvloc,
				&p0,&beta_e,tietol,elcon,
				&islavactdoftie[islavnodeentry-1],ncmat_,
				ntmat_,&iwan));
      
      constantt=min(constant,1.0/atauinvloc);
      derivmode=1;
      
      FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&dgnc,
				   &aninvloc,&p0,&beta_e,elcon,nelcon,
				   &islavactdoftie[islavnodeentry-1],ntmat_,
				   plicon,nplicon,npmat_,ncmat_,tietol,&scal));
      
      if(islavact[islavnodeentry-1]==1 || islavact[islavnodeentry-1]==2){
	if(regmodet==1){
	  derivmode=1;
	  
	  /* perturbed lagrange,tangential direction  
	     (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	  FORTRAN(regularization_slip_lin,(utildep_t,&bp,&atauinvloc,resreg,
					   &derivmode,islavact,lambda_t,
					   lambdatilde_t,&constantt,&debug,
					   &islavnodeentry,n2,t,that,&mu,rslip,
					   ltslip,ltu));
	  
          rphat[0]=0.0;
	  rphat[1]=0.0;
        }else{
	  derivmode=1;
	  atau=1.0/atauinvloc;
	  
	  /* perturbed lagrange,tangential direction  
	     (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	  FORTRAN(regularization_slip_iwan,(&lambda_n,utildep_t,&bp,&atau,
					    resreg,&derivmode,&regmodet,
					    lambdaiwan,lambdaiwanini,
					    &islavnodeentry,n,t,&mu,rslip,
					    ltslip,ltu,&yielded,iit,&debug,
					    &iwan,dut));
	  
	  rphat[0]=0.0;
	  rphat[1]=0.0;
	}	
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
  if(debug==1)printf("\ttrafoNT2: size au_amtil1 %" ITGFORMAT " \n",ifree_amtil1-1);
  
  /* add diagonal terms **/
  
  debug=0;
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
      if(islavactdof[m_flagr[j]-1]<0){
	jslavnodeentry=floor(-islavactdof[m_flagr[j]-1]/10.);	  
	jcol= -islavactdof[m_flagr[j]-1]-10*jslavnodeentry;
	nodem=imastnode[jslavnodeentry-1];
      }else{
	jslavnodeentry=floor(islavactdof[m_flagr[j]-1]/10.);	  
	jcol= islavactdof[m_flagr[j]-1]-10*jslavnodeentry;
	nodem=islavnode[jslavnodeentry-1];	 
      }	     
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
      debug=0;
      trafontspcmpc(n,t,n2,that,&islavnodeentry,
		    nboun,ndirboun,nodeboun,xboun,
		    nmpc,ipompc,nodempc,coefmpc,
		    ikboun,ilboun,ikmpc,ilmpc,
		    nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
		    nmastspc,imastspc,nmspc,nmastmpc,imastmpc,nmmpc,
		    &debug,&node);
      if(debug==1){printf("ddtil %e %" ITGFORMAT " %" ITGFORMAT "\n",
			  Ddtil[jqdtil[node-1]-1],jqdtil[node-1],jqdtil[node]);}
      if(debug==1){printf("node %" ITGFORMAT " idof %" ITGFORMAT " %"
			  ITGFORMAT " %" ITGFORMAT " act %" ITGFORMAT "\n",
			  node,idof1,idof2,idof3,islavact[islavnodeentry-1]);}
      if(debug==1){printf("\t r %" ITGFORMAT " c %" ITGFORMAT " nodem %"
			  ITGFORMAT " idof %" ITGFORMAT " %" ITGFORMAT " %"
			  ITGFORMAT "  \n",jrow,jcol,nodem,
			  nactdof[mt*nodem-3]-1,nactdof[mt*nodem-2]-1,
			  nactdof[mt*nodem-1]-1);}
      
      /* calculate fields needed for Coulomb friction **/
      
      up_n=   u_tilde[(islavnodeentry-1)*3]*n2[0]+
	u_tilde[(islavnodeentry-1)*3+1]*n2[1]+
	u_tilde[(islavnodeentry-1)*3+2]*n2[2];
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
      nlambda_t=sqrt(lambda_t[0]*lambda_t[0]+lambda_t[1]*lambda_t[1]);
      bp=bp_old[islavnodeentry-1];
      
      /* perturbed lagrange,normal direction 
	 (see phd-thesis Sitzmann,Chapter 3.4.1) */
      
      scal=Ddtil[jqdtil[node-1]-1];
      FORTRAN(getcontactparams,(&mu,&regmode,&regmodet,&aninvloc,&atauinvloc,
				&p0,&beta_e,tietol,elcon,
				&islavactdoftie[islavnodeentry-1],ncmat_,
				ntmat_,&iwan));
      
      constantt=min(constant,1.0/atauinvloc);
      derivmode=1;
      FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&dgnc,
				   &aninvloc,&p0,&beta_e,elcon,nelcon,
				   &islavactdoftie[islavnodeentry-1],ntmat_,
				   plicon,nplicon,npmat_,ncmat_,tietol,&scal));
      
      if(islavact[islavnodeentry-1]==1 || islavact[islavnodeentry-1]==2){
	if(regmodet==1){
	  derivmode=1;
	  
	  /* perturbed lagrange,tangential direction 
	     (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	  FORTRAN(regularization_slip_lin,(utildep_t,&bp,&atauinvloc,resreg,
					   &derivmode,islavact,lambda_t,
					   lambdatilde_t,&constantt,&debug,
					   &islavnodeentry,n2,t,that,&mu,rslip,
					   ltslip,ltu));
	  
          rphat[0]=0.0;
	  rphat[1]=0.0;
	}else{    
	  derivmode=1;
	  atau=1.0/atauinvloc;
	  
	  /* perturbed lagrange,tangential direction 
	     (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	  FORTRAN(regularization_slip_iwan,(&lambda_n,utildep_t,&bp,&atau,
					    resreg,&derivmode,&regmodet,
					    lambdaiwan,lambdaiwanini,
					    &islavnodeentry,n,t,&mu,rslip,
					    ltslip,ltu,&yielded,iit,&debug,
					    &iwan,dut));
	  rphat[0]=0.0;
	  rphat[1]=0.0;
	}
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
	    if(debug==1){printf("\t n %e %e %e\n",n2[0],n2[1],n2[2]);}
	    if(debug==1){printf("\t au_bdtil %e %e %e cdof %" ITGFORMAT
				" rdof %" ITGFORMAT " %" ITGFORMAT " %"
				ITGFORMAT " \n \t au_amtil",e1,e2,e3,
				m_flagr[j]-1,idof1,idof2,idof3);}
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    insertas_ws(&irow_amtil2,&(a_flag[idof1]),&(j2),&ifree_amtil2,
			&nzs_amtil2,&contribution,&au_amtil2);
	    if(debug==1){printf(" %e ",contribution);}
	    contribution=(ltslip[0]*e1+ltslip[1]*e2+ltslip[2]*e3);
	    insertas_ws(&irow_amtil2,&(a_flag[idof2]),&(j2),&ifree_amtil2,
			&nzs_amtil2,&contribution,&au_amtil2);
	    if(debug==1){printf(" %e ",contribution);}
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    insertas_ws(&irow_amtil2,&(a_flag[idof3]),&(j2),&ifree_amtil2,
			&nzs_amtil2,&contribution,&au_amtil2);
	    if(debug==1){printf(" %e \n",contribution);}
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
	    if(debug==1){printf("\t ktslip %e %e %e \n",ltslip[3],ltslip[4],
				ltslip[5]);}
	    if(debug==1){printf("\t au_bdtil %e %e %e cdof %" ITGFORMAT
				" rdof %" ITGFORMAT " %" ITGFORMAT
				" \n \t au_amtil",e1,e2,e3,m_flagr[j],idof2,
				idof3);}
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    if(debug==1){printf(" %e ",contribution);}
	    insertas_ws(&irow_amtil2,&(a_flag[idof2]),&(j2),&ifree_amtil2,
			&nzs_amtil2,&contribution,&au_amtil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    if(debug==1){printf(" %e \n",contribution);}
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
	    if(debug==1){printf("\t ktslip %e %e %e \n",ltslip[3],ltslip[4],
				ltslip[5]);}
	    if(debug==1){printf("\t au_bdtil %e %e %e cdof %" ITGFORMAT
				" rdof %" ITGFORMAT " %" ITGFORMAT
				" \n \t au_amtil",e1,e2,e3,m_flagr[j],idof1,
				idof3);}
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    if(debug==1){printf(" %e ",contribution);}
	    insertas_ws(&irow_amtil2,&(a_flag[idof1]),&(j2),&ifree_amtil2,
			&nzs_amtil2,&contribution,&au_amtil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    if(debug==1){printf(" %e \n",contribution);}
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
	    if(debug==1){printf("\t ktslip %e %e %e \n",ltslip[3],ltslip[4],
				ltslip[5]);}
	    if(debug==1){printf("\t au_bdtil %e %e %e cdof %" ITGFORMAT
				" rdof %" ITGFORMAT " %" ITGFORMAT
				" \n \t au_amtil",e1,e2,e3,m_flagr[j],idof1,
				idof2);}
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    if(debug==1){printf(" %e ",contribution);}
	    insertas_ws(&irow_amtil2,&(a_flag[idof1]),&(j2),&ifree_amtil2,
			&nzs_amtil2,&contribution,&au_amtil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    if(debug==1){printf(" %e \n",contribution);}
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
  if(debug==1)printf("\ttrafoNT2: size au_amtil2 %" ITGFORMAT " \n",nzs_amtil2);

  add_rect(au_amtil1,irow_amtil1,jq_amtil1,*row_la,*row_lm,
	   au_amtil2,irow_amtil2,jq_amtil2,*row_la,*row_lm,
	   &au_amtil,&irow_amtil,jq_amtil,&nzs_amtil);
  
  SFREE(au_amtil1);SFREE(irow_amtil1);SFREE(jq_amtil1);
  SFREE(au_amtil2);SFREE(irow_amtil2);SFREE(jq_amtil2);
  if(debug==1)if(debug==1)printf("\ttrafoNT2: size au_amtil %" ITGFORMAT " \n",nzs_amtil);
  debug=0;
  
  /* K_AI **/
  
  nzs_aitil1=jq_dai[*row_li];
  NNEW(irow_aitil1,ITG,nzs_aitil1);
  NNEW(jq_aitil1,ITG,*row_li+1);
  //mastamtil1=NNEW(int,nzs_amtil);
  NNEW(au_aitil1,double,nzs_aitil1);
  jq_aitil1[0]=1;
  ifree_aitil1=1;
  for(j=0;j<*row_li;j++){
    //loop over columns  N   
    jcol=j+1;
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
		    nmpc,ipompc,nodempc,coefmpc,ikboun,ilboun,ikmpc,ilmpc,
		    nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
		    nmastspc,imastspc,nmspc,nmastmpc,imastmpc,nmmpc,
		    &debug,&node);
      
      /* calculate fields needed for Coulomb friction **/
      
      up_n=   u_tilde[(islavnodeentry-1)*3]*n2[0]+
	u_tilde[(islavnodeentry-1)*3+1]*n2[1]+
	u_tilde[(islavnodeentry-1)*3+2]*n2[2];
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
      nlambda_t=sqrt(lambda_t[0]*lambda_t[0]+lambda_t[1]*lambda_t[1]);
      bp=bp_old[islavnodeentry-1];
      
      /* perturbed lagrange,normal direction 
	 (see phd-thesis Sitzmann,Chapter 3.4.1) */
      
      scal=Ddtil[jqdtil[node-1]-1];
      FORTRAN(getcontactparams,(&mu,&regmode,&regmodet,&aninvloc,&atauinvloc,
				&p0,&beta_e,tietol,elcon,
				&islavactdoftie[islavnodeentry-1],ncmat_,
				ntmat_,&iwan));
      
      constantt=min(constant,1.0/atauinvloc);
      derivmode=1;
      
      FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&dgnc,
				   &aninvloc,&p0,&beta_e,elcon,nelcon,
				   &islavactdoftie[islavnodeentry-1],ntmat_,
				   plicon,nplicon,npmat_,ncmat_,tietol,&scal));

      if(islavact[islavnodeentry-1]==1 || islavact[islavnodeentry-1]==2){
	if(regmodet==1){
	  derivmode=1;
	  
	  /* perturbed lagrange,tangential direction 
	     (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	  FORTRAN(regularization_slip_lin,(utildep_t,&bp,&atauinvloc,resreg,
					   &derivmode,islavact,lambda_t,
					   lambdatilde_t,&constantt,&debug,
					   &islavnodeentry,n2,t,that,&mu,rslip,
					   ltslip,ltu));
	  
          rphat[0]=0.0;
	  rphat[1]=0.0;
	}else{
	  derivmode=1;
	  atau=1.0/atauinvloc;
	  
	  /* perturbed lagrange,tangential direction 
	     (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	  FORTRAN(regularization_slip_iwan,(&lambda_n,utildep_t,&bp,&atau,
					    resreg,&derivmode,&regmodet,
					    lambdaiwan,lambdaiwanini,
					    &islavnodeentry,n,t,&mu,rslip,
					    ltslip,ltu,&yielded,iit,&debug,
					    &iwan,dut));
	  
	  rphat[0]=0.0;
	  rphat[1]=0.0;
	}
	if(debug==1){printf("ddtil %e %" ITGFORMAT " %" ITGFORMAT "\n",
			    Ddtil[jqdtil[node-1]-1],jqdtil[node-1],
			    jqdtil[node]);}
	if(debug==1){printf("node %" ITGFORMAT " idof %" ITGFORMAT " %"
			    ITGFORMAT " %" ITGFORMAT " act %" ITGFORMAT "\n",
			    node,idof1,idof2,idof3,islavact[islavnodeentry-1]);}
	if(debug==1){printf("\t dgnc %e \n",dgnc);}
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
	    if(debug==1){printf("\t au_dai %e %e %e \n \t au_aatil",
				e1/Ddtil[jqdtil[node-1]-1],
				e2/Ddtil[jqdtil[node-1]-1],
				e3/Ddtil[jqdtil[node-1]-1]);}
	    if(debug==1){printf("%e ",contribution/Ddtil[jqdtil[node-1]-1]);}
	    insertas_ws(&irow_aitil1,&(irow_dai[i]),&j2,&ifree_aitil1,
			&nzs_aitil1,&contribution,&au_aitil1);
	    contribution=(rslip[0]*e1+rslip[1]*e2+rslip[2]*e3);
	    if(debug==1){printf("%e ",contribution/Ddtil[jqdtil[node-1]-1]);}
	    insertas_ws(&irow_aitil1,&(irow_dai[i+1]),&j2,&ifree_aitil1,
			&nzs_aitil1,&contribution,&au_aitil1);
	    contribution=(rslip[3]*e1+rslip[4]*e2+rslip[5]*e3);
	    if(debug==1){printf("%e \n",contribution/Ddtil[jqdtil[node-1]-1]);}
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
  if(debug==1)printf("\ttrafoNT2: size au_aitil1 %" ITGFORMAT " \n",ifree_aitil1-1); 

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
      if(islavactdof[i_flagr[j]-1]<0){
	jslavnodeentry=floor(-islavactdof[i_flagr[j]-1]/10.);	  
	jcol= -islavactdof[i_flagr[j]-1]-10*jslavnodeentry;
	nodem=imastnode[jslavnodeentry-1];
      }else{
	jslavnodeentry=floor(islavactdof[i_flagr[j]-1]/10.);	  
	jcol= islavactdof[i_flagr[j]-1]-10*jslavnodeentry;
	nodem=islavnode[jslavnodeentry-1];	 
      }	     
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
		    nmpc,ipompc,nodempc,coefmpc,ikboun,ilboun,ikmpc,ilmpc,
		    nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
		    nmastspc,imastspc,nmspc,nmastmpc,imastmpc,nmmpc,
		    &debug,&node);
      
      if(debug==1){printf("ddtil %e %" ITGFORMAT " %" ITGFORMAT "\n",
			  Ddtil[jqdtil[node-1]-1],jqdtil[node-1],jqdtil[node]);}
      if(debug==1){printf("node %" ITGFORMAT " idof %" ITGFORMAT " %"
			  ITGFORMAT " %" ITGFORMAT " act %" ITGFORMAT "\n",
			  node,idof1,idof2,idof3,islavact[islavnodeentry-1]);}
      if(debug==1){printf("\t r %" ITGFORMAT " c %" ITGFORMAT " nodem %"
			  ITGFORMAT " idof %" ITGFORMAT " %" ITGFORMAT " %"
			  ITGFORMAT "  \n",jrow,jcol,nodem,
			  nactdof[mt*nodem-3]-1,nactdof[mt*nodem-2]-1,
			  nactdof[mt*nodem-1]-1);}
      
      /* calculate fields needed for Coulomb friction **/
      
      up_n=   u_tilde[(islavnodeentry-1)*3]*n2[0]+
	u_tilde[(islavnodeentry-1)*3+1]*n2[1]+
	u_tilde[(islavnodeentry-1)*3+2]*n2[2];
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
      nlambda_t=sqrt(lambda_t[0]*lambda_t[0]+lambda_t[1]*lambda_t[1]);
      bp=bp_old[islavnodeentry-1];
      
      /* perturbed lagrange,normal direction 
	 (see phd-thesis Sitzmann,Chapter 3.4.1) */
      
      scal=Ddtil[jqdtil[node-1]-1];
      FORTRAN(getcontactparams,(&mu,&regmode,&regmodet,&aninvloc,&atauinvloc,
				&p0,&beta_e,tietol,elcon,
				&islavactdoftie[islavnodeentry-1],ncmat_,
				ntmat_,&iwan));
      
      constantt=min(constant,1.0/atauinvloc);
      derivmode=1;
      FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&dgnc,
				   &aninvloc,&p0,&beta_e,elcon,nelcon,
				   &islavactdoftie[islavnodeentry-1],ntmat_,
				   plicon,nplicon,npmat_,ncmat_,tietol,&scal));
      
      if(islavact[islavnodeentry-1]==1 || islavact[islavnodeentry-1]==2){ 	    
	if(regmodet==1){
	  derivmode=1;
	  
	  /* perturbed lagrange,tangential direction 
	     (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	  FORTRAN(regularization_slip_lin,(utildep_t,&bp,&atauinvloc,resreg,
					   &derivmode,islavact,lambda_t,
					   lambdatilde_t,&constantt,&debug,
					   &islavnodeentry,n2,t,that,&mu,rslip,
					   ltslip,ltu));
	  
          rphat[0]=0.0;
	  rphat[1]=0.0;
	}else{    
	  derivmode=1;
	  atau=1.0/atauinvloc;
	  
	  /* perturbed lagrange,tangential direction 
	     (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	  FORTRAN(regularization_slip_iwan,(&lambda_n,utildep_t,&bp,&atau,
					    resreg,&derivmode,&regmodet,
					    lambdaiwan,lambdaiwanini,
					    &islavnodeentry,n,t,&mu,rslip,
					    ltslip,ltu,&yielded,iit,&debug,
					    &iwan,dut));
	  rphat[0]=0.0;
	  rphat[1]=0.0;
	}
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
	      //if(debug==1){printf("");}
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
	    if(debug==1){printf("\t n %e %e %e\n",n2[0],n2[1],n2[2]);}
	    if(debug==1){printf("\t au_bdtil %e %e %e cdof %" ITGFORMAT
				" rdof %" ITGFORMAT " %" ITGFORMAT " %"
				ITGFORMAT " \n \t au_amtil",e1,e2,e3,
				i_flagr[j]-1,idof1,idof2,idof3);}
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    insertas_ws(&irow_aitil2,&(a_flag[idof1]),&(j2),&ifree_aitil2,
			&nzs_aitil2,&contribution,&au_aitil2);
	    if(debug==1){printf(" %e ",contribution);}
	    contribution=(ltslip[0]*e1+ltslip[1]*e2+ltslip[2]*e3);
	    insertas_ws(&irow_aitil2,&(a_flag[idof2]),&(j2),&ifree_aitil2,
			&nzs_aitil2,&contribution,&au_aitil2);
	    if(debug==1){printf(" %e ",contribution);}
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    insertas_ws(&irow_aitil2,&(a_flag[idof3]),&(j2),&ifree_aitil2,
			&nzs_aitil2,&contribution,&au_aitil2);
	    if(debug==1){printf(" %e \n",contribution);}
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
	    if(debug==1){printf("\t ktslip %e %e %e \n",ltslip[3],ltslip[4],
				ltslip[5]);}
	    if(debug==1){printf("\t au_bdtil %e %e %e cdof %" ITGFORMAT
				" rdof %" ITGFORMAT " %" ITGFORMAT
				" \n \t au_amtil",e1,e2,e3,i_flagr[j],idof2,
				idof3);}
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    if(debug==1){printf(" %e ",contribution);}
	    insertas_ws(&irow_aitil2,&(a_flag[idof2]),&(j2),&ifree_aitil2,
			&nzs_aitil2,&contribution,&au_aitil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    if(debug==1){printf(" %e \n",contribution);}
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
	    if(debug==1){printf("\t ktslip %e %e %e \n",ltslip[3],ltslip[4],
				ltslip[5]);}
	    if(debug==1){printf("\t au_bdtil %e %e %e cdof %" ITGFORMAT
				" rdof %" ITGFORMAT " %" ITGFORMAT
				" \n \t au_amtil",e1,e2,e3,i_flagr[j],idof1,
				idof3);}
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    if(debug==1){printf(" %e ",contribution);}
	    insertas_ws(&irow_aitil2,&(a_flag[idof1]),&(j2),&ifree_aitil2,
			&nzs_aitil2,&contribution,&au_aitil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    if(debug==1){printf(" %e \n",contribution);}
	    insertas_ws(&irow_aitil2,&(a_flag[idof3]),&(j2),&ifree_aitil2,
			&nzs_aitil2,&contribution,&au_aitil2);
	    i=i+iadd;       	     
	  }else if(idof1>-1 && idof2>-1){
	    //2D auf 3D	
	    if(jrow==1){
	      // k=idof1
	      iadd=0;
	      e1=au_ddtil2i[i];
	      //if(debug==1){printf("");}
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
	    if(debug==1){printf("\t ktslip %e %e %e \n",ltslip[3],ltslip[4],
				ltslip[5]);}
	    if(debug==1){printf("\t au_bdtil %e %e %e cdof %" ITGFORMAT
				" rdof %" ITGFORMAT " %" ITGFORMAT
				" \n \t au_amtil",e1,e2,e3,i_flagr[j],idof1,
				idof2);}
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    if(debug==1){printf(" %e ",contribution);}
	    insertas_ws(&irow_aitil2,&(a_flag[idof1]),&(j2),&ifree_aitil2,
			&nzs_aitil2,&contribution,&au_aitil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    if(debug==1){printf(" %e \n",contribution);}
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
  if(debug==1)printf("\ttrafoNT2: size au_aitil2 %" ITGFORMAT " \n",nzs_aitil2);
  add_rect(au_aitil1,irow_aitil1,jq_aitil1,*row_la,*row_li,
	   au_aitil2,irow_aitil2,jq_aitil2,*row_la,*row_li,
	   &au_aitil,&irow_aitil,jq_aitil,&nzs_aitil);
  
  SFREE(au_aitil1);SFREE(irow_aitil1);SFREE(jq_aitil1);
  SFREE(au_aitil2);SFREE(irow_aitil2);SFREE(jq_aitil2); 
  
  if(debug==1)printf("\ttrafoNT2: au_dai %" ITGFORMAT " au_aitil %" ITGFORMAT "\n",
	 jq_dai[*row_li]-1,jq_aitil[*row_li]-1);

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
      jslavnodeentry=floor(islavactdof[a_flagr[j]-1]/10.);	  
      jcol= islavactdof[a_flagr[j]-1]-10*jslavnodeentry;	  
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
		    nmpc,ipompc,nodempc,coefmpc,ikboun,ilboun,ikmpc,ilmpc,
		    nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
		    nmastspc,imastspc,nmspc,nmastmpc,imastmpc,nmmpc,
		    &debug,&node);
      
      if(debug==1){printf("node %" ITGFORMAT " jcol %" ITGFORMAT " jrow %"
			  ITGFORMAT " idof %" ITGFORMAT " %" ITGFORMAT " %"
			  ITGFORMAT " act %" ITGFORMAT "\n",node,jcol,jrow,
			  idof1,idof2,idof3,islavact[islavnodeentry-1]);}
      
      /* calculate fields needed for Coulomb friction **/
      
      up_n=   u_tilde[(islavnodeentry-1)*3]*n2[0]+
	u_tilde[(islavnodeentry-1)*3+1]*n2[1]+
	u_tilde[(islavnodeentry-1)*3+2]*n2[2];
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
      nlambda_t=sqrt(lambda_t[0]*lambda_t[0]+lambda_t[1]*lambda_t[1]);
      bp=bp_old[islavnodeentry-1];
      
      /* perturbed lagrange,normal direction 
	 (see phd-thesis Sitzmann,Chapter 3.4.1) */
      
      scal=Ddtil[jqdtil[node-1]-1];
      FORTRAN(getcontactparams,(&mu,&regmode,&regmodet,&aninvloc,&atauinvloc,
				&p0,&beta_e,tietol,elcon,
				&islavactdoftie[islavnodeentry-1],ncmat_,
				ntmat_,&iwan));
      
      constantt=min(constant,1.0/atauinvloc);
      derivmode=1;
      
      FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&dgnc,
				   &aninvloc,&p0,&beta_e,elcon,nelcon,
				   &islavactdoftie[islavnodeentry-1],ntmat_,
				   plicon,nplicon,npmat_,ncmat_,tietol,&scal));
      
      if(islavact[islavnodeentry-1]==1 || islavact[islavnodeentry-1]==2){
	if(regmodet==1){
	  derivmode=1;
	  
	  /* perturbed lagrange,tangential direction 
	     (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	  FORTRAN(regularization_slip_lin,(utildep_t,&bp,&atauinvloc,resreg,
					   &derivmode,islavact,lambda_t,lambdatilde_t,&constantt,&debug,
					   &islavnodeentry,n2,t,that,&mu,rslip,ltslip,ltu));
	
          rphat[0]=0.0;
	  rphat[1]=0.0;
	}else{
	  derivmode=1;
	  atau=1.0/atauinvloc;
	  
	  /* perturbed lagrange,tangential direction 
	     (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	  FORTRAN(regularization_slip_iwan,(&lambda_n,utildep_t,&bp,&atau,
					    resreg,&derivmode,&regmodet,
					    lambdaiwan,lambdaiwanini,
					    &islavnodeentry,n,t,&mu,rslip,
					    ltslip,ltu,&yielded,iit,&debug,
					    &iwan,dut));
	  
	  rphat[0]=0.0;
	  rphat[1]=0.0;
	}
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
	    if(debug==1){printf("\t au_daa %e %e %e \n \t au_aatil1",e1,e2,e3);}
	    contribution=(dgnc)*(n[0]*e1+n[1]*e2+n[2]*e3);
	    insertas_ws(&irow_aatil1,&(irow_daa[i]),&j2,&ifree_aatil1,
			&nzs_aatil1,&contribution,&au_aatil1);
	    if(debug==1){printf("%e ",contribution/Ddtil[jqdtil[node-1]-1]);}
	    contribution=(rslip[0]*e1+rslip[1]*e2+rslip[2]*e3 );
	    insertas_ws(&irow_aatil1,&(irow_daa[i+1]),&j2,&ifree_aatil1,
			&nzs_aatil1,&contribution,&au_aatil1);
	    if(debug==1){printf("%e ",contribution/Ddtil[jqdtil[node-1]-1]);}
	    contribution=(rslip[3]*e1+rslip[4]*e2+rslip[5]*e3);
	    insertas_ws(&irow_aatil1,&(irow_daa[i+2]),&j2,&ifree_aatil1,
			&nzs_aatil1,&contribution,&au_aatil1);
	    if(debug==1){printf("%e \n",contribution/Ddtil[jqdtil[node-1]-1]);}
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
  if(debug==1)printf("\ttrafoNT2: size au_aatil1 %" ITGFORMAT "\n",jq_aatil1[*row_la]-1);

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
      if(islavactdof[a_flagr[j]-1]<0){
	jslavnodeentry=floor(-islavactdof[a_flagr[j]-1]/10.);	  
	jcol= -islavactdof[a_flagr[j]-1]-10*jslavnodeentry;
	nodem=imastnode[jslavnodeentry-1];
      }else{
	jslavnodeentry=floor(islavactdof[a_flagr[j]-1]/10.);	  
	jcol= islavactdof[a_flagr[j]-1]-10*jslavnodeentry;
	nodem=islavnode[jslavnodeentry-1];	 
      }	     
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
		    nmpc,ipompc,nodempc,coefmpc,
		    ikboun,ilboun,ikmpc,ilmpc,
		    nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
		    nmastspc,imastspc,nmspc,nmastmpc,imastmpc,nmmpc,
		    &debug,&node);
      if(debug==1){printf("ddtil %e %" ITGFORMAT " %" ITGFORMAT "\n",
			  Ddtil[jqdtil[node-1]-1],jqdtil[node-1],jqdtil[node]);}
      if(debug==1){printf("node %" ITGFORMAT " idof %" ITGFORMAT " %"
			  ITGFORMAT " %" ITGFORMAT " act %" ITGFORMAT "\n",
			  node,idof1,idof2,idof3,islavact[islavnodeentry-1]);}
      if(debug==1){printf("\t r %" ITGFORMAT " c %" ITGFORMAT " nodem %"
			  ITGFORMAT " idof %" ITGFORMAT " %" ITGFORMAT " %"
			  ITGFORMAT "  \n",jrow,jcol,nodem,
			  nactdof[mt*nodem-3]-1,nactdof[mt*nodem-2]-1,
			  nactdof[mt*nodem-1]-1);}
      
      /* calculate fields needed for Coulomb friction **/ 
      up_n=   u_tilde[(islavnodeentry-1)*3]*n2[0]+
	u_tilde[(islavnodeentry-1)*3+1]*n2[1]+
	u_tilde[(islavnodeentry-1)*3+2]*n2[2];
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
      nlambda_t=sqrt(lambda_t[0]*lambda_t[0]+lambda_t[1]*lambda_t[1]);
      bp=bp_old[islavnodeentry-1];
      
      /* perturbed lagrange,normal direction 
	 (see phd-thesis Sitzmann,Chapter 3.4.1) */
      
      scal=Ddtil[jqdtil[node-1]-1];
      FORTRAN(getcontactparams,(&mu,&regmode,&regmodet,&aninvloc,&atauinvloc,
				&p0,&beta_e,tietol,elcon,
				&islavactdoftie[islavnodeentry-1],ncmat_,
				ntmat_,&iwan));
      
      constantt=min(constant,1.0/atauinvloc);
      derivmode=1;
      
      FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&dgnc,
				   &aninvloc,&p0,&beta_e,elcon,nelcon,
				   &islavactdoftie[islavnodeentry-1],ntmat_,
				   plicon,nplicon,npmat_,ncmat_,tietol,&scal));
      
      if(islavact[islavnodeentry-1]==1 || islavact[islavnodeentry-1]==2){ 	    
	if(regmodet==1){
	  derivmode=1;
	  
	  /* perturbed lagrange,tangential direction 
	     (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	  FORTRAN(regularization_slip_lin,(utildep_t,&bp,&atauinvloc,resreg,
					   &derivmode,islavact,lambda_t,
					   lambdatilde_t,&constantt,&debug,
					   &islavnodeentry,n2,t,that,&mu,rslip,
					   ltslip,ltu));
	  
          rphat[0]=0.0;
	  rphat[1]=0.0;
	}else{    
	  derivmode=1;
	  atau=1.0/atauinvloc;
	  
	  /* perturbed lagrange,tangential direction 
	     (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	  FORTRAN(regularization_slip_iwan,(&lambda_n,utildep_t,&bp,&atau,
					    resreg,&derivmode,&regmodet,
					    lambdaiwan,lambdaiwanini,
					    &islavnodeentry,n,t,&mu,rslip,
					    ltslip,ltu,&yielded,iit,&debug,
					    &iwan,dut));
	  
	  rphat[0]=0.0;
	  rphat[1]=0.0;
	}
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
	    if(debug==1){printf("\t n %e %e %e\n",n2[0],n2[1],n2[2]);}
	    if(debug==1){printf("\t au_bdtil %e %e %e cdof %" ITGFORMAT
				" rdof %" ITGFORMAT " %" ITGFORMAT " %"
				ITGFORMAT " \n \t au_aatil",e1,e2,e3,
				a_flagr[j]-1,idof1,idof2,idof3);}
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    insertas_ws(&irow_aatil2,&(a_flag[idof1]),&(j2),&ifree_aatil2,
			&nzs_aatil2,&contribution,&au_aatil2);
	    if(debug==1){printf(" %e ",contribution/Ddtil[jqdtil[node-1]-1]);}
	    contribution=(ltslip[0]*e1+ltslip[1]*e2+ltslip[2]*e3);
	    insertas_ws(&irow_aatil2,&(a_flag[idof2]),&(j2),&ifree_aatil2,
			&nzs_aatil2,&contribution,&au_aatil2);
	    if(debug==1){printf(" %e ",contribution/Ddtil[jqdtil[node-1]-1]);}
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    insertas_ws(&irow_aatil2,&(a_flag[idof3]),&(j2),&ifree_aatil2,
			&nzs_aatil2,&contribution,&au_aatil2);
	    if(debug==1){printf(" %e \n",contribution/Ddtil[jqdtil[node-1]-1]);}
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
	    if(debug==1){printf("\t ktslip %e %e %e \n",ltslip[3],ltslip[4],
				ltslip[5]);}
	    if(debug==1){printf("\t au_bdtil %e %e %e cdof %" ITGFORMAT
				" rdof %" ITGFORMAT " %" ITGFORMAT
				" \n \t au_aatil",e1,e2,e3,a_flagr[j],idof2,
				idof3);}
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    if(debug==1){printf(" %e ",contribution);}
	    insertas_ws(&irow_aatil2,&(a_flag[idof2]),&(j2),&ifree_aatil2,
			&nzs_aatil2,&contribution,&au_aatil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    if(debug==1){printf(" %e \n",contribution);}
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
	      //if(debug==1){printf("");}
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
	    if(debug==1){printf("\t ktslip %e %e %e \n",ltslip[3],ltslip[4],
				ltslip[5]);}
	    if(debug==1){printf("\t au_bdtil %e %e %e cdof %" ITGFORMAT
				" rdof %" ITGFORMAT " %" ITGFORMAT
				" \n \t au_aatil",e1,e2,e3,a_flagr[j],idof1,
				idof3);}
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    if(debug==1){printf(" %e ",contribution);}
	    insertas_ws(&irow_aatil2,&(a_flag[idof1]),&(j2),&ifree_aatil2,
			&nzs_aatil2,&contribution,&au_aatil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    if(debug==1){printf(" %e \n",contribution);}
	    insertas_ws(&irow_aatil2,&(a_flag[idof3]),&(j2),&ifree_aatil2,
			&nzs_aatil2,&contribution,&au_aatil2);
	    i=i+iadd;       	     
	  }else if(idof1>-1 && idof2>-1){
	    //2D auf 3D	
	    if(jrow==1){
	      // k=idof1
	      iadd=0;
	      e1=au_ddtil2a[i];
	      //if(debug==1){printf("");}
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
	    if(debug==1){printf("\t ktslip %e %e %e \n",ltslip[3],ltslip[4],
				ltslip[5]);}
	    if(debug==1){printf("\t au_bdtil %e %e %e cdof %" ITGFORMAT
				" rdof %" ITGFORMAT " %" ITGFORMAT
				" \n \t au_aatil",e1,e2,e3,a_flagr[j],idof1,
				idof2);}
	    contribution=(n2[0]*e1+n2[1]*e2+n2[2]*e3);
	    if(debug==1){printf(" %e ",contribution);}
	    insertas_ws(&irow_aatil2,&(a_flag[idof1]),&(j2),&ifree_aatil2,
			&nzs_aatil2,&contribution,&au_aatil2);
	    contribution=(ltslip[3]*e1+ltslip[4]*e2+ltslip[5]*e3);
	    if(debug==1){printf(" %e \n",contribution);}
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
  if(debug==1)printf("\ttrafoNT2: size au_aatil2 %" ITGFORMAT " \n",nzs_aatil2);
  add_rect(au_aatil1,irow_aatil1,jq_aatil1,*row_la,*row_la,
	   au_aatil2,irow_aatil2,jq_aatil2,*row_la,*row_la,
	   &au_aatil,&irow_aatil,jq_aatil,&nzs_aatil);
  
  SFREE(au_aatil1);SFREE(irow_aatil1);SFREE(jq_aatil1);
  SFREE(au_aatil2);SFREE(irow_aatil2);SFREE(jq_aatil2); 
  
  if(debug==1)printf("\ttrafoNT2: au_daa %" ITGFORMAT " au_aatil %" ITGFORMAT "\n",
	 jq_daa[*row_la]-1,jq_aatil[*row_la]-1);
  
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
      if(debug==1){printf("node %" ITGFORMAT " idof %" ITGFORMAT " %"
			  ITGFORMAT " %" ITGFORMAT " act %" ITGFORMAT "\n",
			  node,idof1,idof2,idof3,islavact[islavnodeentry-1]);}
      
      trafontspcmpc(n,t,n2,that,&islavnodeentry,nboun,ndirboun,nodeboun,xboun,
		    nmpc,ipompc,nodempc,coefmpc,ikboun,ilboun,ikmpc,ilmpc,
		    nslavspc,islavspc,nsspc,nslavmpc,islavmpc,nsmpc,
		    nmastspc,imastspc,nmspc,nmastmpc,imastmpc,nmmpc,
		    &debug,&node);
      
      /* calculate needed fields for coulomb friction **/
      
      up_n=u_tilde[(islavnodeentry-1)*3]*n2[0]+
	u_tilde[(islavnodeentry-1)*3+1]*n2[1]+
	u_tilde[(islavnodeentry-1)*3+2]*n2[2];
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
      nlambda_t=sqrt(lambda_t[0]*lambda_t[0]+lambda_t[1]*lambda_t[1]);
      bp=bp_old[islavnodeentry-1];
      scal=Ddtil[jqdtil[node-1]-1];
      
      FORTRAN(getcontactparams,(&mu,&regmode,&regmodet,&aninvloc,&atauinvloc,
				&p0,&beta_e,tietol,elcon,
				&islavactdoftie[islavnodeentry-1],ncmat_,
				ntmat_,&iwan));
      
      constantt=min(constant,1.0/atauinvloc);
      
      /* perturbed lagrange,normal direction 
	 (see phd-thesis Sitzmann,Chapter 3.4.1) */
      
      derivmode=1;
      FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&dgnc,
				   &aninvloc,&p0,&beta_e,elcon,nelcon,
				   &islavactdoftie[islavnodeentry-1],ntmat_,
				   plicon,nplicon,npmat_,ncmat_,tietol,&scal));
      derivmode=0;
      
      FORTRAN(regularization_gn_c,(&lambda_n,&derivmode,&regmode,&gnc,
				   &aninvloc,&p0,&beta_e,elcon,nelcon,
				   &islavactdoftie[islavnodeentry-1],ntmat_,
				   plicon,nplicon,npmat_,ncmat_,tietol,&scal));

      if(islavact[islavnodeentry-1]==1 || islavact[islavnodeentry-1]==2){
	yielded=0;
	if(regmodet==1){
	  derivmode=1;
	  
	  /* perturbed lagrange,tangential direction 
	     (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	  FORTRAN(regularization_slip_lin,(utildep_t,&bp,&atauinvloc,resreg,
					   &derivmode,islavact,lambda_t,
					   lambdatilde_t,&constantt,&debug,
					   &islavnodeentry,n2,t,that,&mu,rslip,
					   ltslip,ltu));
          rphat[0]=0.0;
	  rphat[1]=0.0;
	}else{
	  derivmode=1;
	  atau=1.0/atauinvloc;
	  
	  /* perturbed lagrange,tangential direction 
	     (see phd-thesis Sitzmann,Chapter 3.4.2) */
	  
	  FORTRAN(regularization_slip_iwan,(&lambda_n,utildep_t,&bp,&atau,
					    resreg,&derivmode,&regmodet,
					    lambdaiwan,lambdaiwanini,
					    &islavnodeentry,n,t,&mu,rslip,
					    ltslip,ltu,&yielded,iit,&debug,
					    &iwan,dut));
	  
	  rphat[0]=0.0;
	  rphat[1]=0.0;
	}
	hpn=gap[islavnodeentry-1]+gnc-dgnc*lambda_n;
	dgnc1=dgnc;
	
	if(debug==1){
	  printf("\t lm= %e %e %e nlm=%e bp=%e\n",lambda_n,lambda_t[0],
		 lambda_t[1],nlambda_t,bp);
	  printf("\t lminit= %e %e utilt= %e %e\n",lambdaini_t[0],
		 lambdaini_t[1],utildep_t[0],utildep_t[1] );
	  if(regmodet==2){
	    printf("\t d= %e %e\n",resreg[0],resreg[1]);
	  }
	  printf("\t uold(%" ITGFORMAT ")= %e %e %e u_t= %e %e %e \n",k,
		 u_tilde[(islavnodeentry-1)*3],u_tilde[(islavnodeentry-1)*3+1],
		 u_tilde[(islavnodeentry-1)*3+2],up_n,
		 utildep_t[0],utildep_t[1]);
	  if(mu>1.E-10){
	    printf("\t rslip1: %e %e %e\n",rslip[0],rslip[1],rslip[2]);
	    printf("\t rslip2: %e %e %e\n",rslip[3],rslip[4],rslip[5]);
	    printf("\t ltslip1: %e %e %e\n",ltslip[0],ltslip[1],ltslip[2]);
	    printf("\t ltslip2: %e %e %e\n",ltslip[3],ltslip[4],ltslip[5]);
	    printf("\t rphat= %e %e \n",rphat[0],rphat[1]);
	    printf("\t ltu= %e %e hpn %e  gap %e\n",ltu[0],ltu[1],hpn,
		   gap[islavnodeentry-1]);
	  }else{
	    printf("\t no friction\n");
	    printf("\t rslip1: %e %e %e\n",rslip[0],rslip[1],rslip[2]);
	    printf("\t rslip2: %e %e %e\n",rslip[3],rslip[4],rslip[5]);	
	    printf("\t ltu= %e %e hpn %e  gap %e\n",ltu[0],ltu[1],hpn,
		   gap[islavnodeentry-1]);	    
	  }
	  
	}
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
	    if(debug==1){printf("\t f_d %e %e %e f_til %e %e %e \n",f_da[k],
				f_da[k+1],f_da[k+2],
				f_atil[k]/Ddtil[jqdtil[node-1]-1],
				f_atil[k+1]/Ddtil[jqdtil[node-1]-1],
				f_atil[k+2]/Ddtil[jqdtil[node-1]-1]);}
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
