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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
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


void prediction(double *uam, ITG *nmethod, double *bet, double *gam, 
		double *dtime,
		ITG *ithermal, ITG *nk, double *veold, double *accold, double *v,
		ITG *iinc, ITG *idiscon, double *vold, ITG *nactdof, ITG *mi,
		ITG *num_cpus){

  ITG j,k,mt=mi[1]+1;
  double dextrapol,scal1,scal2;

  uam[0]=0.;
  uam[1]=0.;
  if(*nmethod==4){
      
    scal1=0.5*(1.-2.**bet)**dtime**dtime;
    scal2=(1.-*gam)**dtime;
	
    if(*ithermal<2){
      preparll(&mt,dtime,veold,&scal1,accold,uam,nactdof,v,vold,&scal2,
	       nk,num_cpus);
    }else if(*ithermal==2){
      for(k=0;k<*nk;++k){
	for(j=0;j<mt;j++){
	  v[mt*k+j]=vold[mt*k+j];
	}
      }
      for(k=0;k<*nk;++k){
	dextrapol=*dtime*veold[mt*k];
	if(fabs(dextrapol)>100.) dextrapol=100.*dextrapol/fabs(dextrapol);
	if((fabs(dextrapol)>uam[1])&&(nactdof[mt*k]>0)) {
	  uam[1]=fabs(dextrapol);}
	v[mt*k]+=dextrapol;
      }
    }else{
      for(k=0;k<*nk;++k){
	for(j=0;j<mt;++j){
	  dextrapol=*dtime*veold[mt*k+j]+scal1*accold[mt*k+j];
	  if((j==0)&&fabs(dextrapol)>100.) dextrapol=100.*dextrapol/fabs(dextrapol);
	  if(j==0){
	    if((fabs(dextrapol)>uam[1])&&(nactdof[mt*k]>0)) {
	      uam[1]=fabs(dextrapol);}
	  }else{
	    if((fabs(dextrapol)>uam[0])&&(nactdof[mt*k+j]>0)) {
	      uam[0]=fabs(dextrapol);}
	  }
	  v[mt*k+j]=vold[mt*k+j]+dextrapol;
	  veold[mt*k+j]=veold[mt*k+j]+scal2*accold[mt*k+j];
	  accold[mt*k+j]=0.;
	}
      }
    }
  }
    
  /* for the static case: extrapolation of the previous increment
     (if any within the same step) */
    
  else{
    if(*iinc>1){
      if(*ithermal<2){
	for(k=0;k<*nk;++k){
	  for(j=1;j<mt;++j){
	    if(*idiscon==0){
	      dextrapol=*dtime*veold[mt*k+j];
	      if((fabs(dextrapol)>uam[0])&&(nactdof[mt*k+j]>0)) {
		uam[0]=fabs(dextrapol);}	
	      v[mt*k+j]=vold[mt*k+j]+dextrapol;
	    }else{
	      v[mt*k+j]=vold[mt*k+j];
	    }
	  }
	}
      }else if(*ithermal==2){
	for(k=0;k<*nk;++k){
	  for(j=0;j<mt;++j){
	    v[mt*k+j]=vold[mt*k+j];
	  }
	}
	for(k=0;k<*nk;++k){
	  if(*idiscon==0){
	    dextrapol=*dtime*veold[mt*k];
	    if(fabs(dextrapol)>100.) dextrapol=100.*dextrapol/fabs(dextrapol);
	    if((fabs(dextrapol)>uam[1])&&(nactdof[mt*k]>0)) {
	      uam[1]=fabs(dextrapol);}	
	    v[mt*k]+=dextrapol;
	  }
	}
      }else{
	for(k=0;k<*nk;++k){
	  for(j=0;j<mt;++j){
	    if(*idiscon==0){
	      dextrapol=*dtime*veold[mt*k+j];
	      if((j==0)&&fabs(dextrapol)>100.) dextrapol=100.*dextrapol/fabs(dextrapol);
	      if(j==0){
		if((fabs(dextrapol)>uam[1])&&(nactdof[mt*k+j]>0)) {
		  uam[1]=fabs(dextrapol);}
	      }else{
		if((fabs(dextrapol)>uam[0])&&(nactdof[mt*k+j]>0)) {
		  uam[0]=fabs(dextrapol);}
	      }	
	      v[mt*k+j]=vold[mt*k+j]+dextrapol;
	    }else{
	      v[mt*k+j]=vold[mt*k+j];
	    }
	  }
	}
      }
    }
    else{
      for(k=0;k<*nk;++k){
	for(j=0;j<mt;++j){
	  v[mt*k+j]=vold[mt*k+j];
	}
      }
    }
  }
  *idiscon=0;

  return;
}
