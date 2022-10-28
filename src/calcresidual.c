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


void calcresidual(ITG *nmethod, ITG *neq, double *b, double *fext, double *f,
		  ITG *iexpl, ITG *nactdof, double *aux2, double *vold,
		  double *vini, double *dtime, double *accold, ITG *nk, double *adb,
		  double *aub, ITG *jq, ITG *irow, ITG *nzl, double *alpha,
		  double *fextini, double *fini, ITG *islavnode, ITG *nslavnode,
		  ITG *mortar, ITG *ntie,double *f_cm,
		  double* f_cs, ITG *mi,ITG *nzs,ITG *nasym,ITG *idamping,
		  double *veold,double *adc,double *auc,double *cvini,double *cv,
		  double *alpham,ITG *num_cpus){

  ITG j,k,mt=mi[1]+1;
  double scal1;
      
  /* residual for a static analysis */
      
  if(*nmethod!=4){
    for(k=0;k<neq[1];++k){
      b[k]=fext[k]-f[k];
    }
  }
      
  /* residual for implicit dynamics */
      
  else if(*iexpl<=1){
    for(k=0;k<*nk;++k){
      if(nactdof[mt*k]>0){
	aux2[nactdof[mt*k]-1]=(vold[mt*k]-vini[mt*k])/(*dtime);}
      for(j=1;j<mt;++j){
	if(nactdof[mt*k+j]>0){
	  aux2[nactdof[mt*k+j]-1]=accold[mt*k+j];}
      }
    }
    if(*nasym==0){
      FORTRAN(op,(&neq[1],aux2,b,adb,aub,jq,irow)); 
    }else{
      FORTRAN(opas,(&neq[1],aux2,b,adb,aub,jq,irow,nzs)); 
    }
    scal1=1.+alpha[0];
    for(k=0;k<neq[0];++k){
      b[k]=scal1*(fext[k]-f[k])-alpha[0]*(fextini[k]-fini[k])-b[k];
    } 
    for(k=neq[0];k<neq[1];++k){
      b[k]=fext[k]-f[k]-b[k];
    } 

    /* correction for damping */

    if(*idamping==1){
      for(k=0;k<*nk;++k){
	if(nactdof[mt*k]>0){
	  aux2[nactdof[mt*k]-1]=0.;}
	for(j=1;j<mt;++j){
	  if(nactdof[mt*k+j]>0){
	    aux2[nactdof[mt*k+j]-1]=veold[mt*k+j];}
	}
      }
      if(*nasym==0){
	FORTRAN(op,(&neq[1],aux2,cv,adc,auc,jq,irow));
      }else{
	FORTRAN(opas,(&neq[1],aux2,cv,adc,auc,jq,irow,nzs)); 
      }
      for(k=0;k<neq[0];++k){
	b[k]-=scal1*cv[k]-alpha[0]*cvini[k];
      }
    }
  }

  /* residual for explicit dynamics */
    
  else{
    res1parll(&mt,nactdof,aux2,vold,vini,dtime,accold,nk,num_cpus);
    scal1=1.+alpha[0];
    res2parll(b,&scal1,fext,f,alpha,fextini,fini,adb,
	      aux2,&neq[0],num_cpus);
    for(k=neq[0];k<neq[1];++k){
      b[k]=fext[k]-f[k]-adb[k]*aux2[k];
    } 

    /* correction for damping */

    if(*idamping==1){
      res3parll(&mt,nactdof,aux2,veold,nk,num_cpus);
      res4parll(cv,alpham,adb,aux2,b,&scal1,alpha,
		cvini,&neq[0],num_cpus);
    }
  }

  return;
}
