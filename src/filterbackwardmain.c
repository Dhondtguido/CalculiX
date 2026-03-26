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

static char *lakonfa1;

static ITG num_cpus;

static ITG *ndesi1,*nodedesi1,*nodedesiinv1,*iregion1,*irow1,*icol1,*jq1,
  nzs1,*nsurfs1,*ipkonfa1,*konfa1,*nodedesipos1,ndesifaces1,*idesiface1;
   
static double *au1=NULL,*ad1=NULL,*aub1=NULL,*adb1=NULL,*area1=NULL,*co1;

void filterbackwardmain(double *co, double *dgdxglob, ITG *nobject,
                         ITG *nk,ITG *nodedesi, ITG *ndesi,char *objectset, 
		         double *xdesi,ITG *nobjectstart,ITG *iponoelfa,
			 ITG *inoelfa,char *lakonfa,ITG *konfa,
			 ITG *ipkonfa,ITG *nodedesiinv,ITG *istartdesi,
			 ITG *ialdesi,ITG *ipkon,char *lakon,ITG *ipoface,
			 ITG *nodface,ITG *kon,ITG *iregion,ITG *isolver,
			 double *dgdx,ITG *ne,ITG *nsurfs){

  /* backward filtering of the sensitivities */

  ITG *nx=NULL,*ny=NULL,*nz=NULL,i,j,*mast=NULL,*irow=NULL,*icol=NULL,*jq=NULL,
    *ipointer=NULL,nzs,symmetryflag=0,inputformat=0,nrhs=1,iobject,inode,icopy,
    *nodedesipos=NULL,*irowf=NULL,*icolf=NULL,*jqf=NULL,ndesifaces,
    *idesiface=NULL;
     
  double *xo=NULL,*yo=NULL,*zo=NULL,*x=NULL,*y=NULL,*z=NULL,filterrad,
    *weighting=NULL,*ad=NULL,*adb=NULL,*au=NULL,*aub=NULL,sigma=0,
    *rhs=NULL,*adf=NULL,*auf=NULL,*temparray=NULL,
    *area=NULL;
    
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
    
  /* check that the number of cpus does not supercede the number
     of design variables */
    
  if(*ndesi<num_cpus) num_cpus=*ndesi;
    
  pthread_t tid[num_cpus];
  
  /*--------------------------------------------------------------------------*/
  /* Determination of mass matrix M and it's derivative K                     */
  /*--------------------------------------------------------------------------*/

  
  /* nodedesipos: indicates for those nodes, which are design variables, their 
     location in nodedesi */
  
  NNEW(nodedesipos,ITG,*nk);
  for(i=0;i<*ndesi;i++){
    nodedesipos[nodedesi[i]-1]=i+1;
  }

  /* Catalogue all external design faces: must contain enough design nodes */
  
  NNEW(idesiface,ITG,*nsurfs+1);
  
  FORTRAN(identdesifaces,(iregion,nsurfs,ipkonfa,lakonfa,konfa,&ndesifaces,
  			  idesiface,nodedesiinv));
  
  RENEW(idesiface,ITG,ndesifaces+1);
  
  /* define the non-zero entries (both matrices have the same structure) */ 

  nzs=20000000;
  NNEW(mast,ITG,nzs);
  NNEW(irow,ITG,1);
  NNEW(icol,ITG,*ndesi);
  NNEW(jq,ITG,*ndesi+1);
  NNEW(ipointer,ITG,*ndesi);

  mastructmm(icol,jq,&mast,&irow,ipointer,&nzs,ndesi,nodedesi,iponoelfa,
             inoelfa,nk,lakonfa,konfa,ipkonfa,nodedesiinv,nodedesipos);
  					      
  SFREE(mast);SFREE(ipointer);    
  RENEW(irow,ITG,nzs);
  NNEW(aub,double,nzs);  
  NNEW(adb,double,*ndesi);
  NNEW(au,double,nzs);
  NNEW(ad,double,*ndesi);
  NNEW(aub1,double,(long long)num_cpus*nzs);  
  NNEW(adb1,double,num_cpus**ndesi);
  NNEW(au1,double,(long long)num_cpus*nzs);
  NNEW(ad1,double,num_cpus**ndesi);
  NNEW(area,double,*ndesi);
  NNEW(area1,double,num_cpus**ndesi);
  

  /* calculate the entries of the filter matrix */
  			
  ndesi1=ndesi;nodedesi1=nodedesi;co1=co;iregion1=iregion;irow1=irow;
  nodedesiinv1=nodedesiinv;icol1=icol;jq1=jq;nzs1=nzs;nsurfs1=nsurfs;
  ipkonfa1=ipkonfa;konfa1=konfa;lakonfa1=lakonfa;nodedesipos1=nodedesipos;
  idesiface1=idesiface;ndesifaces1=ndesifaces;
  
  printf(" Using up to %" ITGFORMAT " cpu(s) for the mass matrix entries.\n\n", num_cpus);
  
  /* create threads and wait */
  
  NNEW(ithread,ITG,num_cpus);
  for(i=0; i<num_cpus; i++)  {
    ithread[i]=i;     
    pthread_create(&tid[i], NULL, (void *)mafillmmmt, (void *)&ithread[i]);
  }
  for(i=0; i<num_cpus; i++)
    pthread_join(tid[i], NULL);
  
  SFREE(ithread);
  
  /* copying accumulating the filter matrices from the different threads */
  
  for(i=0;i<*ndesi;i++){
    for(j=0;j<num_cpus;j++){
      ad[i]+=ad1[i+j**ndesi];
    }
  }
  for(i=0;i<nzs;i++){
    for(j=0;j<num_cpus;j++){
      au[i]+=au1[i+(long long)j*nzs];
    }
  }
  for(i=0;i<*ndesi;i++){
    for(j=0;j<num_cpus;j++){
      adb[i]+=adb1[i+j**ndesi];
    }
  }
  for(i=0;i<nzs;i++){
    for(j=0;j<num_cpus;j++){
      aub[i]+=aub1[i+(long long)j*nzs];
    }
  }
  for(i=0;i<*ndesi;i++){
    for(j=0;j<num_cpus;j++){
      area[i]+=area1[i+j**ndesi];
    }
  }
  
  SFREE(ad1);SFREE(au1);SFREE(aub1);SFREE(adb1);SFREE(nodedesipos);
  SFREE(idesiface);SFREE(area1);    

  /* copying of unfiltered sensitivities in dgdxglob */
  
  icopy=0;        
  for(iobject=*nobjectstart;iobject<*nobject;iobject++){
    FORTRAN(copysens,(&dgdx[iobject**ndesi],dgdxglob,&iobject,&icopy,nk,
		      ndesi,nodedesi));
  }

  /*--------------------------------------------------------------------------*/
  /* Explicit Filtering                                                       */
  /*--------------------------------------------------------------------------*/
  
  if(strcmp1(&objectset[89],"E")==0){
  
    printf(" Calculation of backward filtered sensitivities with explicit method\n\n");

    /* prepare for near3d_se */
    
    NNEW(xo,double,*ndesi);
    NNEW(yo,double,*ndesi);
    NNEW(zo,double,*ndesi);
    NNEW(x,double,*ndesi);
    NNEW(y,double,*ndesi);
    NNEW(z,double,*ndesi);
    NNEW(nx,ITG,*ndesi);
    NNEW(ny,ITG,*ndesi);
    NNEW(nz,ITG,*ndesi);
    
    FORTRAN(prefilter,(co,nodedesi,ndesi,xo,yo,zo,x,y,z,nx,ny,nz,
                       objectset,&filterrad));
    
    /* define structure of the filter matrix */  

    nzs=20000000;
    NNEW(mast,ITG,nzs);
    NNEW(irowf,ITG,1);
    NNEW(icolf,ITG,*ndesi);
    NNEW(jqf,ITG,*ndesi+1);
    NNEW(ipointer,ITG,*ndesi);
    
    mastructfilter(icolf,jqf,&mast,&irowf,ipointer,&nzs,ndesi,
                   nodedesi,xo,yo,zo,x,y,z,nx,ny,nz,&filterrad);
		  
    SFREE(mast);SFREE(ipointer);    
    RENEW(irowf,ITG,nzs);
    NNEW(adf,double,*ndesi);
    NNEW(auf,double,nzs);   
    NNEW(weighting,double,*ndesi);
    NNEW(temparray,double,*ndesi**nobject);
 
    FORTRAN(mafillfilter,(adf,auf,jqf,irowf,ndesi,nodedesi,&filterrad,
                          co,weighting,objectset,xdesi,area));
			       
    FORTRAN(filterbackward_exp,(adf,auf,jqf,irowf,ndesi,nodedesi,dgdxglob,
                                dgdx,nobject,nk,nobjectstart,weighting));
    
    SFREE(weighting);SFREE(irowf);SFREE(jqf);SFREE(icolf);SFREE(adf);
    SFREE(auf);SFREE(xo);SFREE(yo);SFREE(zo);SFREE(x);SFREE(y);SFREE(z);
    SFREE(nx);SFREE(ny);SFREE(nz);SFREE(temparray);

  /*--------------------------------------------------------------------------*/
  /* Implicit Filtering                                                       */
  /*--------------------------------------------------------------------------*/
  
  }else if(strcmp1(&objectset[89],"I")==0){

    printf(" Calculation of backward filtered sensitivities with implicit method\n\n");

    /* assembly of implicit filter matrix */
    
    FORTRAN(filterbackward_imp,(ndesi,au,ad,aub,adb,jq,objectset));
         
    /* Solve the system of equations */
     
    /* LU decomposition of the left hand matrix */

    if(*isolver==0){
#ifdef SPOOLES
      spooles_factor(ad,au,adb,aub,&sigma,icol,irow,ndesi,&nzs,
                     &symmetryflag,&inputformat,&nzs);
#else
      printf("*ERROR in filterbackward: the SPOOLES library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==4){
#ifdef SGI
      token=1;
      sgi_factor(ad,au,adb,aub,&sigma,icol,irow,ndesi,&nzs,token);
#else
      printf("*ERROR in filterbackward: the SGI library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==5){
#ifdef TAUCS
      tau_factor(ad,&au,adb,aub,&sigma,icol,&irow,ndesi,&nzs);
#else
      printf("*ERROR in filterbackward: the TAUCS library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==7){
#ifdef PARDISO
      pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,ndesi,&nzs,
                     &symmetryflag,&inputformat,jq,&nzs);
#else
      printf("*ERROR in filterbackward: the PARDISO library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==8){
#ifdef PASTIX
      pastix_factor_main(ad,au,adb,aub,&sigma,icol,irow,ndesi,&nzs,
                         &symmetryflag,&inputformat,jq,&nzs);
#else
      printf("*ERROR in filterbackward: the PASTIX library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }

    /* solve the system nobject-times */

    NNEW(rhs,double,*ndesi);
    
    for(iobject=*nobjectstart;iobject<*nobject;iobject++){
      for(inode=0;inode<*ndesi;inode++){
         rhs[inode]=dgdx[iobject**ndesi+inode];
      }  
            
      if(*isolver==0){
#ifdef SPOOLES
        spooles_solve(rhs,ndesi);
#endif
      }
      else if(*isolver==4){
#ifdef SGI
        sgi_solve(rhs,token);
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
        tau_solve(rhs,ndesi);
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
        pardiso_solve(rhs,ndesi,&symmetryflag,&inputformat,&nrhs);
#endif
      }
      else if(*isolver==8){
#ifdef PASTIX
        pastix_solve(rhs,ndesi,&symmetryflag,&nrhs);
#endif
      }

      /* copying of normalized sensitivities in dgdxglob */
  
      icopy=1;
      FORTRAN(copysens,(rhs,dgdxglob,&iobject,&icopy,nk,ndesi,nodedesi));

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
      pardiso_cleanup(ndesi,&symmetryflag,&inputformat);
#endif
    }
    else if(*isolver==8){
#ifdef PASTIX
#endif
    }

  /*--------------------------------------------------------------------------*/
  /* No filtering, just normalization with mass matrix                        */
  /*--------------------------------------------------------------------------*/
  
  }else{

    printf(" Scaling of sensitivities with mass matrix\n\n");
    
    /* Solve the system of equations */
     
    /* LU decomposition of the left hand matrix */

    if(*isolver==0){
#ifdef SPOOLES
      spooles_factor(adb,aub,ad,au,&sigma,icol,irow,ndesi,&nzs,
                     &symmetryflag,&inputformat,&nzs);
#else
      printf("*ERROR in filterbackward: the SPOOLES library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==4){
#ifdef SGI
      token=1;
      sgi_factor(adb,aub,ad,au,&sigma,icol,irow,ndesi,&nzs,token);
#else
      printf("*ERROR in filterbackward: the SGI library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==5){
#ifdef TAUCS
      tau_factor(adb,&aub,ad,au,&sigma,icol,&irow,ndesi,&nzs);
#else
      printf("*ERROR in filterbackward: the TAUCS library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==7){
#ifdef PARDISO
      pardiso_factor(adb,aub,ad,au,&sigma,icol,irow,ndesi,&nzs,
                     &symmetryflag,&inputformat,jq,&nzs);
#else
      printf("*ERROR in filterbackward: the PARDISO library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==8){
#ifdef PASTIX
      pastix_factor_main(adb,aub,ad,au,&sigma,icol,irow,ndesi,&nzs,
                         &symmetryflag,&inputformat,jq,&nzs);
#else
      printf("*ERROR in filterbackward: the PASTIX library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }

    /* solve the system nobject-times */
     
    NNEW(rhs,double,*ndesi);
    
    for(iobject=*nobjectstart;iobject<*nobject;iobject++){
      for(inode=0;inode<*ndesi;inode++){
         rhs[inode]=dgdx[iobject**ndesi+inode];
      }  
            
      if(*isolver==0){
#ifdef SPOOLES
        spooles_solve(rhs,ndesi);
#endif
      }
      else if(*isolver==4){
#ifdef SGI
        sgi_solve(rhs,token);
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
        tau_solve(rhs,ndesi);
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
        pardiso_solve(rhs,ndesi,&symmetryflag,&inputformat,&nrhs);
#endif
      }
      else if(*isolver==8){
#ifdef PASTIX
        pastix_solve(rhs,ndesi,&symmetryflag,&nrhs);
#endif
      }

      /* copying of normalized sensitivities in dgdxglob */
  
      icopy=1;
      FORTRAN(copysens,(rhs,dgdxglob,&iobject,&icopy,nk,ndesi,nodedesi));

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
      pardiso_cleanup(ndesi,&symmetryflag,&inputformat);
#endif
    }
    else if(*isolver==8){
#ifdef PASTIX
#endif
    }

  }

  /* free the fields */
  
  SFREE(au),SFREE(ad);SFREE(aub);SFREE(adb);SFREE(irow);SFREE(jq);SFREE(icol);
  SFREE(rhs);SFREE(area);
   
  return;
    
} 

/* subroutine for multithreading of mafillmm */

void *mafillmmmt(ITG *i){

  ITG indexau,indexad,nsurfa,nsurfb,nsurfdelta;

  indexad=*i**ndesi1;
  indexau=*i*nzs1;
    
  nsurfdelta=(ITG)ceil(ndesifaces1/(double)num_cpus);
  nsurfa=*i*nsurfdelta+1;
  nsurfb=(*i+1)*nsurfdelta;
  if(nsurfb>ndesifaces1) nsurfb=ndesifaces1;

  FORTRAN(mafillmm,(co1,nodedesiinv1,iregion1,&au1[indexau],
		    &ad1[indexad],&aub1[indexau],&adb1[indexad],irow1,
		    jq1,ipkonfa1,konfa1,lakonfa1,nodedesipos1,idesiface1,
		    &nsurfa,&nsurfb,&area1[indexad]));	    

  return NULL;
}

#endif
