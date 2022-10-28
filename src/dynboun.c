/*     CalculiX - A 3-dimensional finite element program                   */
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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

#ifdef SPOOLES
#include "spooles.h"
#endif
#ifdef SGI
#include "sgi.h"
#endif
#ifdef TAUCS
#include "tau.h"
#endif
#ifdef PARDISO
#include "pardiso.h"
#endif
#ifdef PASTIX
#include "pastix.h"
#endif

void dynboun(double *amta,ITG *namta,ITG *nam,double *ampli, double *time,
             double *ttime,double *dtime,double *xbounold,double *xboun,
             double *xbounact,ITG *iamboun,ITG *nboun,ITG *nodeboun,
             ITG *ndirboun, double *ad, double *au, double *adb,
             double *aub, ITG *icol, ITG *irow, ITG *neq, ITG *nzs,
             double *sigma, double *b, ITG *isolver,
             double *alpham, double *betam, ITG *nzl,
             ITG *init,double *bact, double *bmin, ITG *jq, 
             char *amname,double *bv, double *bprev, double *bdiff,
             ITG *nactmech, ITG *icorrect, ITG *iprev,double *reltime){

  ITG idiff[3],i,j,ic,ir,im,symmetryflag=0,nrhs=1,inputformat=0;

  double *xbounmin=NULL,*xbounplus=NULL,*bplus=NULL,
    *ba=NULL,deltatime,deltatime2,deltatimesq,timemin,ttimemin,
    timeplus,ttimeplus,*aux=NULL,*b1=NULL,*b2=NULL,*bnew=NULL;

#ifdef SGI
  ITG token=1;
#endif
    
  NNEW(xbounmin,double,*nboun);
  NNEW(xbounplus,double,*nboun);

  /* time increment for the calculation of the change of the
     particular solution (needed to account for nonzero
     SPC's) */

  deltatime=*dtime;
  deltatime2=2.*deltatime;
  deltatimesq=deltatime*deltatime;

  /* the SPC value at timemin is stored in xbounmin */

  if(*init==1){

    /* at the start of a new step it is assumed that the previous step
       has reached steady state (at least for the SPC conditions) */

    for(i=0;i<*nboun;i++){
      xbounmin[i]=xbounold[i];
      xbounact[i]=xbounold[i];
    }
  }
  else{
    timemin=*time-deltatime;
    ttimemin=*ttime-deltatime;
    FORTRAN(temploadmodal,(amta,namta,nam,ampli,&timemin,&ttimemin,dtime,
			   xbounold,xboun,xbounmin,iamboun,nboun,nodeboun,ndirboun,
			   amname,reltime));
  }

  /* the SPC value at timeplus is stored in xbounplus */

  timeplus=*time+deltatime;
  ttimeplus=*ttime+deltatime;
  FORTRAN(temploadmodal,(amta,namta,nam,ampli,&timeplus,&ttimeplus,dtime,
			 xbounold,xboun,xbounplus,iamboun,nboun,nodeboun,ndirboun,
			 amname,reltime));

  NNEW(bplus,double,neq[1]);
  NNEW(ba,double,neq[1]);
  NNEW(b1,double,neq[1]);
  NNEW(b2,double,neq[1]);

  /* check whether boundary conditions changed 
     comparision of min with prev */

  if(*init==1){
    for(i=0;i<*nboun;i++){
      ic=neq[1]+i;
      for(j=jq[ic]-1;j<jq[ic+1]-1;j++){
	ir=irow[j]-1;
	bmin[ir]=bmin[ir]-au[j]*xbounmin[i];
      }
    }
    if(*isolver==0){
#ifdef SPOOLES
      spooles_solve(bmin,&neq[1]);
#endif
    }
    else if(*isolver==4){
#ifdef SGI
      sgi_solve(bmin,token);
#endif
    }
    else if(*isolver==5){
#ifdef TAUCS
      tau_solve(bmin,&neq[1]);
#endif
    }
    if(*isolver==7){
#ifdef PARDISO
      pardiso_solve(bmin,&neq[1],&symmetryflag,&inputformat,&nrhs);
#endif
    }
    if(*isolver==8){
#ifdef PASTIX
      pastix_solve(bmin,&neq[1],&symmetryflag,&nrhs);
#endif
    }
  }

  /* check whether boundary conditions changed 
     comparision of act with min */
  
  idiff[1]=0;
  for(i=0;i<*nboun;i++){
    if(fabs(xbounact[i]-xbounmin[i])>1.e-10){
      idiff[1]=1;
      break;
    }
  }
  if(*init==1){
    for(i=0;i<*nboun;i++){
      ic=neq[1]+i;
      for(j=jq[ic]-1;j<jq[ic+1]-1;j++){
	ir=irow[j]-1;
	bact[ir]=bact[ir]-au[j]*xbounact[i];
      }
    }
    if(*isolver==0){
#ifdef SPOOLES
      spooles_solve(bact,&neq[1]);
#endif
    }
    else if(*isolver==4){
#ifdef SGI
      sgi_solve(bact,token);
#endif
    }
    else if(*isolver==5){
#ifdef TAUCS
      tau_solve(bact,&neq[1]);
#endif
    }
    if(*isolver==7){
#ifdef PARDISO
      pardiso_solve(bact,&neq[1],&symmetryflag,&inputformat,&nrhs);
#endif
    }
    if(*isolver==8){
#ifdef PASTIX
      pastix_solve(bact,&neq[1],&symmetryflag,&nrhs);
#endif
    }
  }

  /* check whether boundary conditions changed 
     comparision of plus with act */
  
  idiff[2]=0;
  for(i=0;i<*nboun;i++){
    if(fabs(xbounplus[i]-xbounact[i])>1.e-10){
      idiff[2]=1;
      break;
    }
  }
  if(idiff[2]==1){
    for(i=0;i<*nboun;i++){
      ic=neq[1]+i;
      for(j=jq[ic]-1;j<jq[ic+1]-1;j++){
	ir=irow[j]-1;
	bplus[ir]=bplus[ir]-au[j]*xbounplus[i];
      }
    }
    if(*isolver==0){
#ifdef SPOOLES
      spooles_solve(bplus,&neq[1]);
#endif
    }
    else if(*isolver==4){
#ifdef SGI
      sgi_solve(bplus,token);
#endif
    }
    else if(*isolver==5){
#ifdef TAUCS
      tau_solve(bplus,&neq[1]);
#endif
    }
    if(*isolver==7){
#ifdef PARDISO
      pardiso_solve(bplus,&neq[1],&symmetryflag,&inputformat,&nrhs);
#endif
    }
    if(*isolver==8){
#ifdef PASTIX
      pastix_solve(bplus,&neq[1],&symmetryflag,&nrhs);
#endif
    }
  }
  
  if((idiff[1]!=0)||(idiff[2]!=0)){

    /* present value is not zero */

    if(idiff[2]==0){
      for(i=0;i<neq[1];i++){bplus[i]=bact[i];}
    }
    for(i=0;i<neq[1];i++){
	  
      /* bv is the velocity */
	  
      bv[i]=(bplus[i]-bmin[i])/deltatime2;
	  
      /* ba is the acceleration */
	  
      ba[i]=(bmin[i]-2.*bact[i]+bplus[i])/deltatimesq;
	  
      b1[i]=ba[i]+*alpham*bv[i];
      b2[i]=*betam*bv[i];

      bmin[i]=bact[i];
      bact[i]=bplus[i];
    }
    NNEW(bnew,double,neq[1]);
    FORTRAN(op,(&neq[1],b1,bplus,adb,aub,jq,irow));
    for(i=0;i<neq[1];i++){bnew[i]=-bplus[i];}
    FORTRAN(op,(&neq[1],b2,bplus,ad,au,jq,irow));
    if(*icorrect==2){
      for(i=0;i<neq[1];i++){
	bnew[i]-=bplus[i];
	b[i]+=bnew[i];
      }
    }else if(*icorrect==0){
      for(i=0;i<neq[1];i++){
	bnew[i]-=bplus[i];
	bdiff[i]=bnew[i]-bprev[i];
	b[i]+=bdiff[i];
	//	      printf("dynboun %e,%e,%e,%e\n",bprev[i],bnew[i],bdiff[i],b[i]);
      }
      memcpy(&bprev[0],&bnew[0],sizeof(double)*neq[1]);
    }else{
      for(i=0;i<neq[1];i++){
	bnew[i]-=bplus[i];
	bdiff[i]+=bnew[i]-bprev[i];
	b[i]+=bdiff[i];
      }
      memcpy(&bprev[0],&bnew[0],sizeof(double)*neq[1]);
    }
    SFREE(bnew);
    *nactmech=neq[1];
    *iprev=1;
  }else if((*iprev!=0)&&(*icorrect!=2)){

    /* present value of b is zero, previous value was not zero */

    if(*icorrect==0){
      for(i=0;i<neq[1];i++){
	bdiff[i]=-bprev[i];
	b[i]+=bdiff[i];
	//	      printf("dynboun %e,%e,%e,%e\n",bprev[i],bdiff[i],b[i]);
      }
      //	  memset(&bprev[0],0.,sizeof(double)*neq[1]);
      DMEMSET(bprev,0,neq[1],0.);
    }else{
      for(i=0;i<neq[1];i++){
	bdiff[i]+=-bprev[i];
	b[i]+=bdiff[i];
      }
      //	  memset(&bprev[0],0.,sizeof(double)*neq[1]);
      DMEMSET(bprev,0,neq[1],0.);
    }
    *nactmech=neq[1];
    *iprev=0;
  }
  
  SFREE(xbounmin);SFREE(xbounplus);
  SFREE(bplus);SFREE(ba);SFREE(b1);SFREE(b2);
  
  return;
}
