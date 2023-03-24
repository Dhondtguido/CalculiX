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


void dudsmain(ITG *isolver,double *au,double *ad,double *aub,double*adb,
	 ITG *icol,ITG *irow,ITG *jq,ITG *neq,ITG *nzs,double *df,ITG *jqs,
	 ITG *irows,ITG *ndesi,double *duds){
	               
  /* calculating the matrix duds */

  ITG symmetryflag=0,inputformat=0,idesvar,idof,j,nrhs=1,num_cpus;
  
  double sigma=0,*dudsvec=NULL;	         
		 
  /* variables for multithreading procedure */
  
  ITG sys_cpus,*ithread=NULL;
  char *env,*envloc,*envsys;
  
  num_cpus=0;
  sys_cpus=0;
  
  /* explicit user declaration prevails */
  
  envsys=getenv("NUMBER_OF_CPUS");
  if(envsys){
      sys_cpus=atoi(envsys);
      if(sys_cpus<0) sys_cpus=0;
  }
  
  /* automatic detection of available number of processors */
  
  if(sys_cpus==0){
      sys_cpus = getSystemCPUs();
      if(sys_cpus<1) sys_cpus=1;
  }
  
  /* local declaration prevails, if strictly positive */
  
  envloc = getenv("CCX_NPROC_SENS");
  if(envloc){
      num_cpus=atoi(envloc);
      if(num_cpus<0){
  	  num_cpus=0;
      }else if(num_cpus>sys_cpus){
  	  num_cpus=sys_cpus;
      }
      
  }
  
  /* else global declaration, if any, applies */
  
  env = getenv("OMP_NUM_THREADS");
  if(num_cpus==0){
      if (env)
  	  num_cpus = atoi(env);
      if (num_cpus < 1) {
  	  num_cpus=1;
      }else if(num_cpus>sys_cpus){
  	  num_cpus=sys_cpus;
      }
  }
  
  pthread_t tid[num_cpus];
 
  /* factorize the system */

	if(*isolver==0){
#ifdef SPOOLES
	    spooles_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
			   &symmetryflag,&inputformat,&nzs[2]);
#else
	    printf(" *ERROR in dudsmain: the SPOOLES library is not linked\n\n");
	    FORTRAN(stop,());
#endif
	}
	else if(*isolver==4){
#ifdef SGI
	    token=1;
	    sgi_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],token);
#else
	    printf(" *ERROR in dudsmain: the SGI library is not linked\n\n");
	    FORTRAN(stop,());
#endif
	}
	else if(*isolver==5){
#ifdef TAUCS
	    tau_factor(ad,&au,adb,aub,&sigma,icol,&irow,&neq[1],&nzs[1]);
#else
	    printf(" *ERROR in dudsmain: the TAUCS library is not linked\n\n");
	    FORTRAN(stop,());
#endif
	}
	else if(*isolver==7){
#ifdef PARDISO
	    pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
			   &symmetryflag,&inputformat,jq,&nzs[2]);
#else
	    printf(" *ERROR in dudsmain: the PARDISO library is not linked\n\n");
	    FORTRAN(stop,());
#endif
    }
	else if(*isolver==8){
#ifdef PASTIX
	    pastix_factor_main(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
			   &symmetryflag,&inputformat,jq,&nzs[2]);
#else
	    printf(" *ERROR in dudsmain: the PASTIX library is not linked\n\n");
	    FORTRAN(stop,());
#endif
	}

  /* Computation of the matrix duds */
  /* duds = K^-1 * ( dF/ds + dK/ds * u )
  /*      = K^-1 * df */
  /* ndesi system of equations have to be solved to determine duds */
  
  NNEW(dudsvec,double,neq[1]);
  
  for(idesvar=0;idesvar<*ndesi;idesvar++){
  
     /* initialize vector dudsvec */
        
     for(j=0;j<neq[1];j++){
     
       dudsvec[j]=0.;
     
     }
     
     /* assembly of the vector taken from df */
     
     for(j=jqs[idesvar];j<jqs[idesvar+1]-1;j++){

	idof=irows[j-1]-1;
	dudsvec[idof]=df[j-1];

     }

     /*Alternative way to assemble the vector dudsvec /*
     /*FORTRAN(dudsassembly,(jqs,irows,neq,df,dudsvec,&idesvar));*/
          
     /* solve the system */
		    
		    if(*isolver==0){
#ifdef SPOOLES
			spooles_solve(dudsvec,&neq[1]);
#endif
		    }
		    else if(*isolver==4){
#ifdef SGI
			sgi_solve(dudsvec,token);
#endif
		    }
		    else if(*isolver==5){
#ifdef TAUCS
			tau_solve(dudsvec,&neq[1]);
#endif
		    }
		    else if(*isolver==7){
#ifdef PARDISO
			pardiso_solve(dudsvec,&neq[1],&symmetryflag,&inputformat,&nrhs);
#endif
                    }	      
		    else if(*isolver==8){
#ifdef PASTIX
			pastix_solve(dudsvec,&neq[1],&symmetryflag,&nrhs);
#endif
                    }	      

     /* copy results vector in duds */

     for(j=0;j<neq[1];j++){
       
        duds[(neq[1]*idesvar)+j]=dudsvec[j];
     
     }
  
  }
  
  SFREE(dudsvec);        
	    
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
	    pardiso_cleanup(&neq[1],&symmetryflag,&inputformat);
#endif
	}
	else if(*isolver==8){
#ifdef PASTIX
#endif
	}
       
  return;
  
} 

#endif
