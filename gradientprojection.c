/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2015 Guido Dhondt                     */

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

#ifdef ARPACK

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
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
#ifdef MATRIXSTORAGE
#include "matrixstorage.h"
#endif
#ifdef PARDISO
#include "pardiso.h"
#endif
#ifdef PASTIX
#include "pastix.h"
#endif

void gradientprojection(ITG *nobject,char *objectset,double *dgdxglob,
			double *g0,ITG *ndesi,ITG *nodedesi,ITG *nk,
			ITG *isolver,char *set,ITG *nset,ITG *istartset,
			ITG *iendset,ITG *ialset,
			double *gradproj,char *gradprojname,ITG *nactive,
			double *objnorm,ITG *ipoacti,ITG *iconstacti,
			ITG *inameacti,ITG *nnlconst){
               
  /* finding a feasible direction based on the sensitivity information */
   
  ITG nzss,*mast1=NULL,*irows=NULL,*icols=NULL,*jqs=NULL,*ipointer=NULL,
    symmetryflag=0,inputformat=0,i,iconst,iter=0,nactiveold,*ipoactiold=NULL,
    *iconstactiold=NULL,iscaleflag,nrhs=1,*inameactiold=NULL,istart;
             
  double *au=NULL,*ad=NULL,*adb=NULL,*aub=NULL,sigma=0,*rhs=NULL,
    *vector=NULL,*xlambd=NULL,*xtf=NULL,*dfdx=NULL;  
  
  /* scale all sensitivity vectors to unit length 1 */
  
  iscaleflag=1;   
  for(i=0;i<*nobject;i++){
    istart=i+1;
    FORTRAN(scalesen,(dgdxglob,dfdx,nk,nodedesi,ndesi,objectset,&iscaleflag,
                      &istart));
  }
  
  nactiveold=*nactive+1;
  	       
  while((*nactive<nactiveold)&&(*nactive>0)){
  
    /* Initialization of final gradient vector */
    if(iter>0){
      for(i=0;i<3**nk;i++){
        gradproj[i]=0.;
      }
    }
  
    nactiveold=*nactive;
    iter=iter+1;
   
    /* determining the structure of the N-Matrix */

    nzss=20000000;
    NNEW(mast1,ITG,nzss);
    NNEW(irows,ITG,1);
    NNEW(icols,ITG,*nactive);
    NNEW(jqs,ITG,*nactive+1);
    NNEW(ipointer,ITG,*nactive);
    NNEW(rhs,double,*nactive);
    NNEW(xlambd,double,*nactive);
    NNEW(xtf,double,*nactive);
    NNEW(vector,double,*ndesi);
   
    mastructnmatrix(icols,jqs,&mast1,&irows,ipointer,&nzss,nactive,nnlconst);

    RENEW(irows,ITG,nzss);
    SFREE(ipointer);
   
    /* determining the entries of the N-Matrix */
  
    NNEW(ad,double,*nactive);
    NNEW(au,double,nzss);    
  
    FORTRAN(nmatrix,(ad,au,jqs,irows,ndesi,nodedesi,dgdxglob,nactive,
  		     nobject,nnlconst,ipoacti,nk));

    /* Calculate inverse of the N-matrix */

    NNEW(adb,double,*nactive);
    NNEW(aub,double,nzss); 

    for(i=0;i<*nactive;i++){
      adb[i]=1.0;
    }
     
    /* LU decomposition of the left hand matrix */

    if(*isolver==0){
#ifdef SPOOLES
      spooles_factor(ad,au,adb,aub,&sigma,icols,irows,nactive,&nzss,
                     &symmetryflag,&inputformat,&nzss);
#else
      printf("*ERROR in projectgrad: the SPOOLES library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==4){
#ifdef SGI
      token=1;
      sgi_factor(ad,au,adb,aub,&sigma,icols,irows,nactive,&nzss,token);
#else
      printf("*ERROR in projectgrad: the SGI library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==5){
#ifdef TAUCS
      tau_factor(ad,&au,adb,aub,&sigma,icols,&irows,nactive,&nzss);
#else
      printf("*ERROR in projectgrad: the TAUCS library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==7){
#ifdef PARDISO
      pardiso_factor(ad,au,adb,aub,&sigma,icols,irows,nactive,&nzss,
		     &symmetryflag,&inputformat,jqs,&nzss);
#else
      printf("*ERROR in projectgrad: the PARDISO library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==8){
#ifdef PASTIX
      pastix_factor_main(ad,au,adb,aub,&sigma,icols,irows,nactive,&nzss,
			 &symmetryflag,&inputformat,jqs,&nzss);
#else
      printf("*ERROR in projectgrad: the PASTIX library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
        
    /* solve the system nactive-times */

    for(i=0;i<*nactive;i++){
      xtf[i]=0.00;
    }
        
    FORTRAN(preprojectgrad,(vector,ndesi,nodedesi,dgdxglob,nactive,nobject,
			    nnlconst,ipoacti,nk,rhs,objectset,xtf));
        
    /* Calculate the projected gradient and the lagrange multipliers */
        
    for(iconst=1;iconst<=*nactive;iconst++){
            
      for(i=0;i<*nactive;i++){
        rhs[i]=0.00;
      }
            
      rhs[iconst-1]=1.0;
            
      /* solve the system */
            
      if(*isolver==0){
#ifdef SPOOLES
        spooles_solve(rhs,nactive);
#endif
      }
      else if(*isolver==4){
#ifdef SGI
        sgi_solve(rhs,token);
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
        tau_solve(rhs,nactive);
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
        pardiso_solve(rhs,nactive,&symmetryflag,&inputformat,&nrhs);
#endif
      }
      else if(*isolver==8){
#ifdef PASTIX
        pastix_solve(rhs,nactive,&symmetryflag,&nrhs);
#endif
      }
            
      for(i=0;i<*ndesi;i++){
        vector[i]=0.00;
      }
            
      /* carry out matrix multiplications */
           
      FORTRAN(projectgrad,(vector,ndesi,nodedesi,dgdxglob,nactive,nobject,
                           nnlconst,ipoacti,nk,rhs,&iconst,objectset,xlambd,
			   xtf,objnorm,gradproj,g0));
    }
        
    /* clean the system */
        
    if(*isolver==0){
#ifdef SPOOLES
      spooles_cleanup();
#endif
    }
    else if(*isolver==4){
#ifdef SGI
      sgi_cleanup(token);
#endif
    }
    else if(*isolver==5){
#ifdef TAUCS
      tau_cleanup();
#endif
    }
    else if(*isolver==7){
#ifdef PARDISO
      pardiso_cleanup(nactive,&symmetryflag,&inputformat);
#endif
    }
    else if(*isolver==8){
#ifdef PASTIX
#endif
    }

    /* write the results of the gradient projection in the dat-file */
        
    FORTRAN(writelm,(&iter,xlambd,nactive,nnlconst,objectset,nobject,
                     ipoacti,iconstacti,inameacti,nodedesi,dgdxglob,nk));
        
    /* check if the langrange multipliers of the active constraints
       have the "correct" sign ("correct" in the sense of active constraints) */
                 
    NNEW(ipoactiold,ITG,*nactive);
    NNEW(iconstactiold,ITG,*nactive);
    NNEW(inameactiold,ITG,*nactive);
      
    FORTRAN(checkprojectgrad,(&nactiveold,nactive,ipoacti,ipoactiold,
			      objectset,xlambd,nnlconst,iconstacti,
			      iconstactiold,inameacti,inameactiold,g0,
			      nobject,ndesi,nodedesi,dgdxglob,nk));

    SFREE(mast1);SFREE(irows);SFREE(icols);SFREE(jqs);SFREE(ad);SFREE(au);
    SFREE(adb);SFREE(aub);SFREE(rhs);SFREE(xlambd);SFREE(xtf);SFREE(vector);
    SFREE(ipoactiold);SFREE(iconstactiold);SFREE(inameactiold);
        
  }
  
  /* calucaltion of final feasable direction */
       
  FORTRAN(calcfeasibledirection_gp,(ndesi,nodedesi,dgdxglob,nactive,nobject,nk,
				    gradproj,gradprojname));
        
  return;
  
} 

#endif
