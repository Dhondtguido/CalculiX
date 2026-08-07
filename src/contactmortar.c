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

/*
  function to include contact conditions with the dual mortar method in the 
  transformed system
 
  see phd-thesis Sitzmann, Algorithm 2, p.71
 
  Author: Saskia Sitzmann   */

void contactmortar(ITG *ncont,ITG *ntie,char *tieset,ITG *nset,char *set,
		   ITG *istartset,ITG *iendset,ITG *ialset,ITG *itietri,
		   char *lakon,ITG *ipkon,ITG *kon,ITG *koncont,ITG *ne,
		   double *cg,double *straight,double *co,double *vold,
		   ITG *ielmat,double *elcon,ITG *istep,ITG *iinc,ITG *iit,
		   ITG *ncmat_,ITG *ntmat_,ITG *ne0,double *vini,ITG *nmethod,
		   ITG *neq,ITG *nzs,ITG *nactdof,ITG *itiefac,ITG *islavsurf,
		   ITG *islavnode,ITG *imastnode,ITG *nslavnode,ITG *nmastnode,
		   double *ad,double **aup,double *b,ITG **irowp,ITG *icol,
		   ITG *jq,ITG *imastop,ITG *iponoels,ITG *inoels,ITG *nzsc,
		   double **aucp,double *adc,ITG **irowcp,ITG *jqc,
		   ITG *islavact,double *gap,double *slavnor,double *slavtan,
		   double *bhat,ITG **irowbdp,ITG *jqbd,double **aubdp,
		   ITG **irowbdtilp,ITG *jqbdtil,double **aubdtilp,
		   ITG **irowbdtil2p,ITG *jqbdtil2,double **aubdtil2p,
		   ITG **irowddp,ITG *jqdd,double **auddp,ITG **irowddtilp,
		   ITG *jqddtil,double **auddtilp,ITG **irowddtil2p,
		   ITG *jqddtil2,double **auddtil2p,ITG **irowddinvp,
		   ITG *jqddinv,double **auddinvp,ITG *irowt,ITG *jqt,
		   double *aut,ITG *irowtinv,ITG *jqtinv,
		   double *autinv,ITG *mi,ITG *ipe,ITG *ime,double *tietol,
		   double *cstress,double *cstressini,
		   double *bp_old,ITG *nk,ITG *nboun,
		   ITG *ndirboun,ITG *nodeboun,double *xboun,ITG *nmpc,
		   ITG *ipompc,ITG *nodempc,double *coefmpc,ITG *ikboun,
		   ITG *ilboun,ITG *ikmpc,ITG *ilmpc,
		   ITG *nslavspc,ITG *islavspc,ITG *nslavmpc,
		   ITG *islavmpc,ITG *nmastmpc,ITG *imastmpc,
		   double *pslavdual,ITG *islavactdof,
		   ITG *islavtie,double *plicon,ITG *nplicon,ITG *npmat_,
		   ITG *nelcon,double *dtime,ITG *islavnodeinv,double **Bdp,
		   ITG **irowbp,ITG *jqb,double **Bdhelpp,ITG **irowbhelpp,
		   ITG *jqbhelp,double **Ddp,ITG **irowdp,ITG *jqd,
		   double **Ddtilp,ITG **irowdtilp,ITG *jqdtil,double **Bdtilp,
		   ITG **irowbtilp,ITG *jqbtil,
		   double *bet,double *cfsinitil,
		   double *reltime,ITG *ithermal,double *plkcon,ITG *nplkcon){
  
  ITG i,j,k,ntrimax,*nx=NULL,*ny=NULL,*nz=NULL,nintpoint=0,
    nzsbd,*irowbd=NULL,*irowdd=NULL,*irowddinv=NULL,*irowddtil=NULL,
    *irowbdtil=NULL,nzs2,*irowddtil2=NULL,*irowbdtil2=NULL,l,nstart,kflag,
    ntri,ii,regmode,derivmode,*irowc=NULL,*imastsurf=NULL,
    *irow=NULL,*irowb=NULL,*irowbhelp=NULL,*irowd=NULL,
    *irowdtil=NULL,*irowbtil=NULL,nacti,ninacti,nnogap,nstick,nnolm,nnoslav,nzsbdtil,
    nzsbdtil2;
    
  double *xo=NULL,*yo=NULL,*zo=NULL,*x=NULL,*y=NULL,*z=NULL,*aubd=NULL,
    *cstresstil=NULL,scal,*audd=NULL,*auddtil=NULL,*auddtil2=NULL,
    *auddinv=NULL,*auc=NULL,*pmastsurf=NULL,*gapmints=NULL,*au=NULL,
    *pslavsurf=NULL,*aubdtil=NULL,*aubdtil2=NULL,*Bd=NULL,*Bdhelp=NULL,
    *Dd=NULL,*Ddtil=NULL,*Bdtil=NULL,
    mu,fkninv,fktauinv,p0,beta,*rs=NULL,*rsb=NULL;
  
  double aninvloc,gnc,xlnold,lt[2],ltold;
  
  double *u_old=NULL,*u_oldt=NULL;
  ITG mt=mi[1]+1,nodes,jj;
  
  irow=*irowp;au=*aup;auc=*aucp;irowc=*irowcp;
  aubd=*aubdp;irowbd=*irowbdp;
  aubdtil=*aubdtilp;irowbdtil=*irowbdtilp;
  aubdtil2=*aubdtil2p;irowbdtil2=*irowbdtil2p;
  irowdd=*irowddp;audd=*auddp;
  irowddtil=*irowddtilp;auddtil=*auddtilp;
  irowddtil2=*irowddtil2p;auddtil2=*auddtil2p;
  irowddinv=*irowddinvp;auddinv=*auddinvp;
  irowb=*irowbp;Bd =*Bdp;irowbhelp=*irowbhelpp;Bdhelp =*Bdhelpp;
  irowd=*irowdp;Dd =*Ddp;
  irowdtil=*irowdtilp;Ddtil =*Ddtilp;
  irowbtil=*irowbtilp;Bdtil =*Bdtilp;
  
  NNEW(rs,double,(mt)**nk);
  NNEW(rsb,double,neq[1]);
  
  /* create field islavactdof;
     coupling the equation degrees of freedom with the corresponding
     slave nodes and doing the same for the master nodes */
  
  FORTRAN(genislavactdof,(ntie,tieset,nactdof,nslavnode,
			  nmastnode,imastnode,islavactdof,
			  islavnode,mi,ithermal));
  
  /* The update of the normals and tangentials as well as 
     the segmentation of the contact surface needed
     for the calculation of the coupling matrices is done only once 
     per increment, since a combined fix-point
     Newton approach in implemented, see phd-thesis Saskia Sitzmann, 
     Chapter 3 introduction  */
  
  if(*iit==1){
      
    /* update the location of the center of gravity of 
       the master triangles and the coefficients of their
       bounding planes needed for the search algorithm in 
       slavintmortar->neartriangle */
      
    FORTRAN(updatecont,(koncont,ncont,co,vold,
			cg,straight,mi));
      
    /* determining the size of the auxiliary fields 
       (needed for the master triangle search for any
       given location on the slave faces) */
      
    ntrimax=0;	
    for(i=0;i<*ntie;i++){
      if(itietri[2*i+1]-itietri[2*i]+1>ntrimax){		
	ntrimax=itietri[2*i+1]-itietri[2*i]+1;} 	
    }
      
    /* For the first step, first increment, first iteration an initial guess
       for the active set is generated analogous to node-to-surface penalty */
      
    if ((*iinc==1)&&(*iit==1)&&(*istep==1)){	    
      NNEW(xo,double,ntrimax);	    
      NNEW(yo,double,ntrimax);	    
      NNEW(zo,double,ntrimax);	    
      NNEW(x,double,ntrimax);	    
      NNEW(y,double,ntrimax);	    
      NNEW(z,double,ntrimax);	   
      NNEW(nx,ITG,ntrimax);	   
      NNEW(ny,ITG,ntrimax);	    
      NNEW(nz,ITG,ntrimax);	    
      FORTRAN(genfirstactif,(tieset,ntie,itietri,cg,
			     straight,co,vold,xo,yo,zo,x,y,z,nx,ny,nz,
			     istep,iinc,iit,mi,imastop,
			     nslavnode,islavnode,
			     set,nset,istartset,iendset,ialset,
			     islavact,tietol));
	  
      SFREE(xo);SFREE(yo);SFREE(zo);SFREE(x);SFREE(y);SFREE(z);
      SFREE(nx);SFREE(ny);SFREE(nz);
    }
    fflush(stdout);

    /* counting the active/inactive.... nodes */
    
    nacti=0;ninacti=0;nnogap=0;nstick=0;nnolm=0;nnoslav=0;
    for (i=0;i<*ntie;i++){	    	    
      if(tieset[i*(81*3)+80]=='C'){	   	      
	for(j=nslavnode[i];j<nslavnode[i+1];j++){	
	  if(islavact[j]<0){
	    islavact[j]=-3;}
	  if(islavact[j]==2){nacti++;}
	  if(islavact[j]==0){ninacti++;}
	  if(islavact[j]==-3){nnoslav++;}	    
	}	    
      }	
    }
      
    /* calculating the normals,tangents in the nodes of the slave
       surface */
      
    FORTRAN(nortanslav,(tieset,ntie,ipkon,kon,lakon,set,
			co,vold,nset,islavsurf,itiefac,
			islavnode,nslavnode,slavnor,slavtan,mi));
    
    /* Calculating the location of the matched slave/master
       integration points, see phd-thesis Saskia Sitzmann, Appendix A */
    
    NNEW(xo,double,ntrimax);	
    NNEW(yo,double,ntrimax);	
    NNEW(zo,double,ntrimax);	
    NNEW(x,double,ntrimax);	
    NNEW(y,double,ntrimax);	
    NNEW(z,double,ntrimax);	
    NNEW(nx,ITG,ntrimax);	
    NNEW(ny,ITG,ntrimax);	
    NNEW(nz,ITG,ntrimax);
      
    NNEW(imastsurf,ITG,66);	
    NNEW(gapmints,double,66);	
    NNEW(pmastsurf,double,132);	
    NNEW(pslavsurf,double,198);	
    islavsurf[1]=0;	
    for(i=0;i<*ntie;i++){	    
      ii=i+1;	    
      if(tieset[i*(81*3)+80]=='C'){		
	nstart=itietri[2*i]-1;		
	ntri=(itietri[2*i+1]-nstart);		
	for(j=0;j<ntri;j++){		    
	  xo[j]=cg[(nstart+j)*3];		    
	  x[j]=xo[j];		   
	  nx[j]=j+1;		    
	  yo[j]=cg[(nstart+j)*3+1];		    
	  y[j]=yo[j];		    
	  ny[j]=j+1;		    
	  zo[j]=cg[(nstart+j)*3+2];		    
	  z[j]=zo[j];		    
	  nz[j]=j+1;		
	}
	kflag=2;		
	FORTRAN(dsort,(x,nx,&ntri,&kflag));		
	FORTRAN(dsort,(y,ny,&ntri,&kflag));		
	FORTRAN(dsort,(z,nz,&ntri,&kflag));		
	for(l=itiefac[2*i];l<=itiefac[2*i+1];l++){		    
	  RENEW(imastsurf,ITG,nintpoint+ntri*8*7);		    
	  RENEW(gapmints,double,nintpoint+ntri*8*7);		    
	  RENEW(pmastsurf,double,2*(nintpoint+ntri*8*7));
	  RENEW(pslavsurf,double,3*(nintpoint+ntri*8*7));
	  FORTRAN(slavintmortar,(ntie,itietri,ipkon,kon,lakon,
				 straight,&nintpoint,koncont,co,vold,
				 xo,yo,zo,x,y,z,nx,ny,nz,iinc,
				 islavsurf,imastsurf,pmastsurf,
				 islavnode,nslavnode,imastop,gapmints,
				 islavact,mi,ncont,ipe,ime,pslavsurf,
				 &ii,&l,&ntri,tietol,reltime,nmethod));
	}	    
      }	
    }
    if (nintpoint!=0){	    
      RENEW(imastsurf,ITG,nintpoint);	
    }else{	    
      RENEW(imastsurf,ITG,1);	
    }	
    if (nintpoint!=0){	    
      RENEW(gapmints,double,nintpoint);	
    }else{	    
      RENEW(gapmints,double,1);	
    }
    if (nintpoint!=0){    
      RENEW(pmastsurf,double,2*nintpoint);	
    }else{	    
      RENEW(pmastsurf,double,2);	
    }		
    if (nintpoint!=0){	    
      RENEW(pslavsurf,double,3*nintpoint);	
    }else{	    
      RENEW(pslavsurf,double,3);
    }
    SFREE(xo);SFREE(yo);SFREE(zo);SFREE(x);SFREE(y);SFREE(z);SFREE(nx);    
    SFREE(ny);SFREE(nz);
      
    /* check SPC's and MPC's on slave nodes for compatibility and set all 
       slave nodes involved in SPCs/MPCs to no-LM nodes */
      
    FORTRAN(remlagrangemult,(ntie,tieset,islavnode,imastnode,nslavnode,
			     nmastnode,islavact,nodempc,nmpc,ipompc));
      
    nacti=0;ninacti=0;nnogap=0;nstick=0;nnolm=0;nnoslav=0;	
    for (i=0;i<*ntie;i++){	
      if(tieset[i*(81*3)+80]=='C'){	    
	for(j=nslavnode[i];j<nslavnode[i+1];j++){		
	  if(islavact[j]==2){nacti++;}
	  if(islavact[j]==1){nstick++;}		
	  if(islavact[j]==0){ninacti++;}
	  if(islavact[j]==-1){nnogap++;}
	  if(islavact[j]==-2){nnolm++;}
	  if(islavact[j]==-3){nnoslav++;}
	}	 
      }	 	
    }
      
    /* calculating the coeffs of dual basis functions (Sitzmann, Chapter 3.3.) 
       and redistribute contributions of nogap
       and noLM nodes other slave nodes (Sitzmann, Chapter 4.3.)*/
      
    FORTRAN(gendualcoeffs,(tieset,ntie,ipkon,kon,lakon,co,vold,islavact,
			   islavsurf,itiefac,islavnode,nslavnode,
			   mi,pslavsurf,pslavdual));
    
    /* calculate all mortar coupling matrices as well as the dual gap 
       via the segmentation of the slave surface */
      
    nzsbd=6*nslavnode[*ntie];
    nzsbdtil=6*nslavnode[*ntie];
    
    bdfill(&irowbd,jqbd,&aubd,&nzsbd,&irowbdtil,jqbdtil,&aubdtil,&nzsbdtil,
	   &irowbdtil2,jqbdtil2,&aubdtil2,&nzsbdtil2,&irowdd,jqdd,&audd,
	   &irowddtil,jqddtil,&auddtil,&irowddtil2,jqddtil2,&auddtil2,
	   &irowddinv,jqddinv,&auddinv,irowt,jqt,aut,irowtinv,
	   jqtinv,autinv,ntie,ipkon,kon,lakon,nslavnode,nmastnode,
	   imastnode,islavnode,islavsurf,imastsurf,pmastsurf,itiefac,tieset,
	   neq,nactdof,co,vold,iponoels,inoels,mi,gapmints,gap,pslavsurf,
	   pslavdual,&nintpoint,slavnor,nk,nmpc,ipompc,nodempc,coefmpc,ikmpc,
	   ilmpc,nslavmpc,islavmpc,nmastmpc,imastmpc,
	   iit,iinc,islavactdof,islavact,islavnodeinv,&Bd,&irowb,jqb,
	   &Bdhelp,&irowbhelp,jqbhelp,&Dd,&irowd,jqd,&Ddtil,&irowdtil,jqdtil,
	   &Bdtil,&irowbtil,jqbtil,ithermal);
    
    SFREE(imastsurf);SFREE(pmastsurf);SFREE(gapmints);SFREE(pslavsurf);      
    fflush(stdout);  
  }

  /* end of loop preparing fields at the start of each increment */
  
  /* get uhat_k-1 for first increment and first iteration**/
  
  if(*iit==1){     
    NNEW(u_old,double,3*nslavnode[*ntie]);
    NNEW(u_oldt,double,3*nslavnode[*ntie]);
    NNEW(cstresstil,double,mt*nslavnode[*ntie]);
    for (i=0;i<*ntie;i++){      	
      if(tieset[i*(81*3)+80]=='C'){      	
	for(j=nslavnode[i];j<nslavnode[i+1];j++){	    
	  nodes=islavnode[j];
	  for(jj=jqd[nodes-1]-1;jj<jqd[nodes-1+1]-1;jj++){      
	    u_oldt[(islavnodeinv[irowd[jj]-1]-1)*3]+=
	      Dd[jj]*(vold[mt*(nodes)-3]-vini[mt*(nodes)-3]);
	    u_oldt[(islavnodeinv[irowd[jj]-1]-1)*3+1]+=
	      Dd[jj]*(vold[mt*(nodes)-2]-vini[mt*(nodes)-2]);
	    u_oldt[(islavnodeinv[irowd[jj]-1]-1)*3+2]+=
	      Dd[jj]*(vold[mt*(nodes)-1]-vini[mt*(nodes)-1]);
	    u_old[(islavnodeinv[irowd[jj]-1]-1)*3]+=
	      Dd[jj]*(vold[mt*(nodes)-3]);
	    u_old[(islavnodeinv[irowd[jj]-1]-1)*3+1]+=
	      Dd[jj]*(vold[mt*(nodes)-2]);
	    u_old[(islavnodeinv[irowd[jj]-1]-1)*3+2]+=
	      Dd[jj]*(vold[mt*(nodes)-1]);
	  }	  
	}	
      }    
    }
    for (i=0;i<*ntie;i++){      	
      if(tieset[i*(81*3)+80]=='C'){      	
	for(j=nmastnode[i];j<nmastnode[i+1];j++){	    
	  nodes=imastnode[j];
	  for(jj=jqb[nodes-1]-1;jj<jqb[nodes-1+1]-1;jj++){      
	    u_oldt[(islavnodeinv[irowb[jj]-1]-1)*3]+=
	      Bd[jj]*(vold[mt*(nodes)-3]-vini[mt*(nodes)-3]);
	    u_oldt[(islavnodeinv[irowb[jj]-1]-1)*3+1]+=
	      Bd[jj]*(vold[mt*(nodes)-2]-vini[mt*(nodes)-2]);
	    u_oldt[(islavnodeinv[irowb[jj]-1]-1)*3+2]+=
	      Bd[jj]*(vold[mt*(nodes)-1]-vini[mt*(nodes)-1]);
	    u_old[(islavnodeinv[irowb[jj]-1]-1)*3]+=
	      Bd[jj]*(vold[mt*(nodes)-3]);
	    u_old[(islavnodeinv[irowb[jj]-1]-1)*3+1]+=
	      Bd[jj]*(vold[mt*(nodes)-2]);
	    u_old[(islavnodeinv[irowb[jj]-1]-1)*3+2]+=
	      Bd[jj]*(vold[mt*(nodes)-1]);
	  }	  
	}	
      }    
    }
    for (i=0;i<*ntie;i++){      	
      if(tieset[i*(81*3)+80]=='C'){      	
	for(j=nslavnode[i];j<nslavnode[i+1];j++){	    
	  nodes=islavnode[j];
	  for(jj=jqdtil[nodes-1]-1;jj<jqdtil[nodes-1+1]-1;jj++){     
	    cstresstil[(islavnodeinv[irowdtil[jj]-1]-1)*(mt)]+=
	      Ddtil[jj]*cstress[(islavnodeinv[nodes-1]-1)*(mt)+0];
	    cstresstil[(islavnodeinv[irowdtil[jj]-1]-1)*(mt)+1]+=
	      Ddtil[jj]*cstress[(islavnodeinv[nodes-1]-1)*(mt)+1];
	    cstresstil[(islavnodeinv[irowdtil[jj]-1]-1)*(mt)+2]+=
	      Ddtil[jj]*cstress[(islavnodeinv[nodes-1]-1)*(mt)+2];
	  }	  
	}	
      }    
    }     
  }
  
  nacti=0;ninacti=0;nnogap=0;nstick=0;nnolm=0;nnoslav=0;		
  for (i=0;i<*ntie;i++){	  
    if(tieset[i*(81*3)+80]=='C'){	    
      for(j=nslavnode[i];j<nslavnode[i+1];j++){
	
	/* adjust active set for first iteration of first increment **/
	
	if((*iit==1)){		 
	  xlnold=cstresstil[(j)*(mt)+0]*slavnor[(j*3)+0]+
	    cstresstil[(j)*(mt)+1]*slavnor[(j*3)+1]+
	    cstresstil[(j)*(mt)+2]*slavnor[(j*3)+2];
	  lt[0]=cstresstil[(j)*(mt)+0]*slavtan[(j*6)+0]+
	    cstresstil[(j)*(mt)+1]*slavtan[(j*6)+1]+
	    cstresstil[(j)*(mt)+2]*slavtan[(j*6)+2];	
	  lt[1]=cstresstil[(j)*(mt)+0]*slavtan[(j*6)+3]+
	    cstresstil[(j)*(mt)+1]*slavtan[(j*6)+4]+
	    cstresstil[(j)*(mt)+2]*slavtan[(j*6)+5];		 
	  ltold=sqrt(lt[0]*lt[0]+lt[1]*lt[1]);
	  xlnold=cstresstil[(j)*(mt)+0]*slavnor[(j*3)+0]+
	    cstresstil[(j)*(mt)+1]*slavnor[(j*3)+1]+
	    cstresstil[(j)*(mt)+2]*slavnor[(j*3)+2];
	  FORTRAN(getcontactparams,(&mu,&regmode,&fkninv,&fktauinv,
				    &p0,&beta,tietol,elcon,&i,ncmat_,ntmat_));
	  derivmode=0;
	  if(islavact[j]>-1){
	    scal=Ddtil[jqdtil[islavnode[j]-1]-1];
	  }else{
	    scal=0.0;}
	  FORTRAN(regularization_gn_c,(&xlnold,&derivmode,&regmode,&gnc,
				       &aninvloc,&p0,&beta,elcon,nelcon,&i,
				       ntmat_,plicon,nplicon,npmat_,ncmat_,
				       tietol,&scal));
	  if(mu>1.E-10){     
	    bp_old[j]=mu*(xlnold);
	    
	  }else{
	    bp_old[j]=(xlnold);
	  }
	  jj=j+1;
	  ltold=sqrt((lt[0])*(lt[0])+(lt[1])*(lt[1]));
	  
	  // in case of NO friction node must set "slip"
	  
	  if(mu>1.E-10 ){
	    if(*iinc==1){
	      if(islavact[j]==0 && gap[j]<1.e-9 ) {
		islavact[j]=1;}	    
	      if(islavact[j]>0 && bp_old[j]<1.e-14){
		bp_old[j]=1;}
	      if(islavact[j]>0 && ltold <1.e-5){
		islavact[j]=1;}// first step
	    }
	  }else{
	    if(*iinc==1){
	      if(gap[j]>1E-10 && islavact[j]>0 ){
		islavact[j]=0;
		bp_old[j]=0.0;}
	      if(gap[j]<1E-10 && islavact[j]==0 ){
		islavact[j]=2;}
	      if(islavact[j]==1){
		islavact[j]=2;}
	    }
	  }	 
	}	
	if(islavact[j]==2){nacti++;}
	if(islavact[j]==1){nstick++;}		
	if(islavact[j]==0){ninacti++;}
	if(islavact[j]==-1){nnogap++;}
	if(islavact[j]==-2){nnolm++;}	
	if(islavact[j]==-3){nnoslav++;}
      }	  
    }	
  }
  
  if(*iinc==1 && *iit==1 && *istep==1 ){
    
    /* set initial value for bp_old in first iteration of first
       increment of first step */
    
    for (i=0;i<*ntie;i++){	 
      if(tieset[i*(81*3)+80]=='C'){	 
	for(j=nslavnode[i];j<nslavnode[i+1];j++){	  
	  bp_old[j]=1.0; 	 
	}	 
      }       
    }   
  }
  if(*iit==1){
    SFREE(u_old);SFREE(u_oldt);SFREE(cstresstil);}
  
  /* modifying the stiffnes matrix K with the coupling matrices;the
     expanded (symmetric) matrix is described in asymmetric form by
     the fields auc, adc, irowc, jqc and nzsc, bhat */
  
  nzsbd=jqbd[neq[1]]-1;		
  *nzsc=nzs[1];	
  nzs2=nzs[1];

  /* modifying the stiffnes matrix K with the coupling matrices;
     embedding of the contact conditions
     the expanded (symmetric) matrix is described in asymmetric form by
     the fields auc, adc, irowc, jqc and nzsc, bhat */
  
  /* k needed in semi-smooth Newton in tangential direction:
     k=1 stick is assumed in all nodes, k=2 stick or slip is assumed 
     according to active set entry*/
  
  if(*iit==1 && *iinc==1){k=1;}else{k=2;}
    
    /* normal formulation (Sitzmann, Chapter 4.1.)*/
    
  multimortar(&au,ad,&irow,jq,&nzs2,&auc,adc,&irowc,jqc,nzsc,aubd,irowbd,
	      jqbd,aubdtil,irowbdtil,jqbdtil,aubdtil2,irowbdtil2,jqbdtil2,
	      irowdd,jqdd,audd,irowddtil2,jqddtil2,auddtil2,irowddinv,
	      jqddinv,auddinv,Bd,irowb,jqb,Dd,irowd,jqd,Ddtil,irowdtil,
	      jqdtil,neq,b,bhat,islavnode,imastnode,nslavnode,nmastnode,
	      islavact,islavactdof,gap,slavnor,slavtan,vold,vini,cstress,
	      cstressini,bp_old,nactdof,ntie,mi,nk,nboun,ndirboun,nodeboun,
	      xboun,nmpc,ipompc,nodempc,coefmpc,ikboun,ilboun,ikmpc,ilmpc,
	      nslavspc,islavspc,nslavmpc,islavmpc,
	      tieset,islavtie,
	      nelcon,elcon,tietol,ncmat_,ntmat_,plicon,nplicon,npmat_,dtime,
	      irowt,jqt,aut,irowtinv,jqtinv,autinv,
	      islavnodeinv,&k,nmethod,bet,ithermal,
	      plkcon,nplkcon);

  nzs[0]=jq[neq[1]]-1; 
  nzs[1]=jq[neq[1]]-1;
  
  fflush(stdout);
  
  /* calculating icol and icolc (needed for SPOOLES) */
  
  for(i=0;i<neq[1];i++){	
    icol[i]=jq[i+1]-jq[i];
  } 
  
  /* nzlc is the number of the rightmost column with 
     nonzero off-diagonal terms */
  
  *irowp=irow;*aup=au;
  *aucp=auc;*irowcp=irowc;
  *aubdp=aubd;*irowbdp=irowbd;
  *aubdtilp=aubdtil;*irowbdtilp=irowbdtil;
  *aubdtil2p=aubdtil2;*irowbdtil2p=irowbdtil2;
  *auddp=audd;*irowddp=irowdd;
  *auddinvp=auddinv;*irowddinvp=irowddinv;   
  *auddtilp=auddtil;*irowddtilp=irowddtil;
  *auddtil2p=auddtil2;*irowddtil2p=irowddtil2;
  *Bdp=Bd;*irowbp= irowb;
  *Bdhelpp=Bdhelp;*irowbhelpp= irowbhelp;
  *Ddp=Dd;*irowdp= irowd;
  *Ddtilp=Ddtil;*irowdtilp= irowdtil;
  *Bdtilp=Bdtil;*irowbtilp= irowbtil;
  
  SFREE(rs);SFREE(rsb);
  fflush(stdout);
  return;
}
