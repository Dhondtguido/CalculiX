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


void randomfieldmain(ITG *kon,ITG *ipkon,char *lakon,ITG *ne,ITG *nmpc,
		     ITG *nactdof,ITG *mi,ITG *nodedesi,ITG *ndesi,
		     ITG *istartdesi,
		     ITG *ialdesi,double *co,double *physcon,ITG *isolver,
		     ITG *ntrans,
		     ITG *nk,ITG *inotr,double *trab,char *jobnamec,ITG *nboun,
		     double *cs,	 
		     ITG *mcs,ITG *inum,ITG *nmethod,ITG *kode,char *filab,
		     ITG *nstate_,
		     ITG *istep,char *description,char *set,ITG *nset,
		     ITG *iendset,
		     char *output,ITG *istartset,ITG *ialset,double *extnor,
		     ITG *irandomtype,double *randomval,ITG *irobustdesign,
		     ITG *ndesibou,
		     ITG *nodedesibou,ITG *nodedesiinvbou){
	               
  /* generating the random field for robust assessment */

  char bmat[2]="I", which[3]="SA", howmny[2]="A",*objectset=NULL,*orname=NULL;
  
  ITG nzss,nzsd,*mast1=NULL,*irows=NULL,*icols=NULL,*jqs=NULL,*ipointer=NULL,
    *nx=NULL,*ny=NULL,*nz=NULL,symmetryflag=0,inputformat=0,mei[4],ido,ldz,
    iparam[11],mxiter,info,ncv,jrow,lworkl,nev,ipntr[14],m,
    *select=NULL,rvec=1,i,k,kflag,idesvar,node,inorm=0,irand=1,
    iinc=1,noddiam=-1,ngraph,iobject,icoordinate,nrhs=1,ithermal=0,
    ishape=0,inode,*irowd=NULL,*icold=NULL,*jqd=NULL,nzsc,*irowc=NULL,
    *icolc=NULL,*jqc=NULL,nevold=0,imodes=0,niter=0,ifeasd=0;
  
  double *xo=NULL,*yo=NULL,*zo=NULL,*x=NULL,*y=NULL,*z=NULL,*au=NULL,*ad=NULL,
    *adb=NULL,*aub=NULL,*resid=NULL,*workd=NULL,*workl=NULL,tol,
    *temp_array=NULL,*randval=NULL,fmin,fmax,pi,trace=0.0,dsum=0.,
    *d=NULL,*zz=NULL,*acvector=NULL,*randfield=NULL,abserr,
    actreliability,*stn=NULL,reliability,corrlen,*acscalar=NULL,
    *add=NULL,*aud=NULL,sigma=0,*rhs=NULL,*vector=NULL,*adbd=NULL,
    *aubd=NULL,*auc=NULL,delta,time=0.;
	         
  reliability=physcon[10];
  corrlen=physcon[11];

  /* Assigning values to eigenvalue parameters mei */
  
  mei[1]=*ndesi;  //number of requested Lanczos vectors
  mei[2]=1000;    //number of max iterations
  mei[3]=0;       //matrix storage
  
  /* copying the frequency parameters */

  pi=4.*atan(1.);
  mxiter=1000;
  tol=0.01;
  fmin=2*pi*-1;
  fmax=2*pi*-1;
  
  /* prepare for near3d_se; routine "prefilter" can w/o any changes be taken 
     from the routine "filtermain.c" */

  NNEW(xo,double,*ndesi);
  NNEW(yo,double,*ndesi);
  NNEW(zo,double,*ndesi);
  NNEW(x,double,*ndesi);
  NNEW(y,double,*ndesi);
  NNEW(z,double,*ndesi);
  NNEW(nx,ITG,*ndesi);
  NNEW(ny,ITG,*ndesi);
  NNEW(nz,ITG,*ndesi);
  
  printf("\n Sorting the design variables for the computation of the covariance matrix\n\n");
  
  for(m=0;m<*ndesi;m++){
    xo[m]=co[3*(nodedesi[m]-1)];
    x[m]=xo[m];
    nx[m]=m+1;
    yo[m]=co[3*(nodedesi[m]-1)+1];
    y[m]=yo[m];
    ny[m]=m+1;
    zo[m]=co[3*(nodedesi[m]-1)+2];
    z[m]=zo[m];
    nz[m]=m+1;
  }
  kflag=2;
  FORTRAN(dsort,(x,nx,ndesi,&kflag));
  FORTRAN(dsort,(y,ny,ndesi,&kflag));
  FORTRAN(dsort,(z,nz,ndesi,&kflag));

  /* determining the structure of the covariance matrix */

  printf(" Determining the structure of the covariance matrix\n\n");

  nzss=20000000;
  NNEW(mast1,ITG,nzss);
  NNEW(irows,ITG,1);
  NNEW(icols,ITG,*ndesi);
  NNEW(jqs,ITG,*ndesi+1);
  NNEW(ipointer,ITG,*ndesi);

  mastructrand(icols,jqs,&mast1,&irows,ipointer,&nzss,ndesi,&corrlen,
               xo,yo,zo,x,y,z,nx,ny,nz);

  SFREE(mast1);SFREE(ipointer);    
  RENEW(irows,ITG,nzss);

  SFREE(xo);SFREE(yo);SFREE(zo);SFREE(x);SFREE(y);SFREE(z);SFREE(nx);
  SFREE(ny);SFREE(nz);
	      
  /* determining the entries of the autocovariance matrix */

  printf(" Calculating the entries of the covariance matrix\n\n");
  
  NNEW(ad,double,*ndesi);
  NNEW(au,double,nzss);    

  FORTRAN(autocovmatrix,(co,ad,au,jqs,irows,ndesi,nodedesi,&corrlen,
			 randomval,irobustdesign));


  /* ***************************************************************** */
  /* start of conditional random field computation                     */
  /* ***************************************************************** */

  /* in case of a conditional random field, assemble the D-matrix and 
     the C-matrix */
  /* have to be assembled:  */
  /* Cov^c = Cov - C * D^(-1) * C^(T) */
  /* D-matrix: contains the covariance information of the boundary nodes */
  /* C-matrix: contains the linking information between the design variables and
     the boundary nodes subset) */
  
  if(irobustdesign[2]==1){
    nzsc=20000000;
    NNEW(mast1,ITG,nzsc);
    NNEW(irowc,ITG,1);
    NNEW(icolc,ITG,*ndesi);
    NNEW(jqc,ITG,*ndesi+1);
    NNEW(ipointer,ITG,*ndesi);

    mastructcmatrix(icolc,jqc,&mast1,&irowc,ipointer,&nzsc,ndesibou,
		    nodedesibou,nodedesiinvbou,jqs,irows,icols,ndesi,nodedesi);

    SFREE(mast1);SFREE(ipointer);    
    RENEW(irowc,ITG,nzsc);

    NNEW(auc,double,nzsc);    

    FORTRAN(cmatrix,(ad,au,jqs,irows,icols,ndesi,nodedesi,auc,jqc,irowc,
		     nodedesibou));
  
    nzsd=20000000;
    NNEW(mast1,ITG,nzsd);
    NNEW(irowd,ITG,1);
    NNEW(icold,ITG,*ndesi);
    NNEW(jqd,ITG,*ndesi+1);
    NNEW(ipointer,ITG,*ndesi);
  
    mastructdmatrix(icold,jqd,&mast1,&irowd,ipointer,&nzsd,ndesibou,
		    nodedesibou,nodedesiinvbou,jqs,irows,icols,ndesi,nodedesi);

    SFREE(mast1);SFREE(ipointer);    
    RENEW(irowd,ITG,nzsd);

    printf(" Calculating the entries of the conditional covariance matrix\n\n");
  
    NNEW(add,double,*ndesibou);
    NNEW(aud,double,nzsd);    
     
    FORTRAN(dmatrix,(ad,au,jqs,irows,icols,ndesi,nodedesi,add,aud,
		     jqd,irowd,ndesibou,nodedesibou));

    NNEW(adbd,double,*ndesibou);
    NNEW(aubd,double,nzsd); 
  
    for(i=0;i<*ndesibou;i++){
      adbd[i]=1.0;
    }

    /* LU decomposition of the left hand matrix */

    if(*isolver==0){
#ifdef SPOOLES
      spooles_factor(add,aud,adbd,aubd,&sigma,icold,irowd,ndesibou,&nzsd,
		     &symmetryflag,&inputformat,&nzsd);
#else
      printf(" *ERROR in randomfieldmain: the SPOOLES library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==4){
#ifdef SGI
      token=1;
      sgi_factor(add,aud,adbd,aubd,&sigma,icold,irowd,ndesibou,&nzsd,token);
#else
      printf(" *ERROR in randomfieldmain: the SGI library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==5){
#ifdef TAUCS
      tau_factor(add,&aud,adbd,aubd,&sigma,icold,&irowd,ndesibou,&nzsd);
#else
      printf(" *ERROR in randomfieldmain: the TAUCS library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==7){
#ifdef PARDISO
      pardiso_factor(add,aud,adbd,aubd,&sigma,icold,irowd,ndesibou,&nzsd,
		     &symmetryflag,&inputformat,jqd,&nzsd);
#else
      printf(" *ERROR in randomfieldmain: the PARDISO library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==8){
#ifdef PASTIX
      pastix_factor_main(add,aud,adbd,aubd,&sigma,icold,irowd,ndesibou,&nzsd,
		     &symmetryflag,&inputformat,jqd,&nzsd);
#else
      printf(" *ERROR in randomfieldmain: the PASTIX library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }

    NNEW(rhs,double,*ndesibou);
    NNEW(vector,double,*ndesi);

    for(idesvar=1;idesvar<=*ndesi;idesvar++){
	
      /* initialize rhs */
     
      for(i=0;i<*ndesibou;i++){
	rhs[i]=0.0;
      }
      for(i=0;i<*ndesi;i++){
	vector[i]=0.0;
      }

      /* assembly of rhs */
        
      FORTRAN(precondrandomfield,(auc,jqc,irowc,rhs,&idesvar));

      /* solve the system */
	    
      if(*isolver==0){
#ifdef SPOOLES
	spooles_solve(rhs,ndesibou);
#endif
      }
      else if(*isolver==4){
#ifdef SGI
	sgi_solve(rhs,token);
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
	tau_solve(rhs,ndesibou);
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	pardiso_solve(rhs,ndesibou,&symmetryflag,&inputformat,&nrhs);
#endif
      }
      else if(*isolver==8){
#ifdef PASTIX
	pastix_solve(rhs,ndesibou,&symmetryflag,&nrhs);
#endif
      }
	
      /* carry out matrix multiplications */
	    
      FORTRAN(condrandomfield,(ad,au,jqs,irows,ndesi,rhs,vector,&idesvar,jqc,
			       auc,irowc));
       	    
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
      pardiso_cleanup(ndesibou,&symmetryflag,&inputformat);
#endif
    }
    else if(*isolver==8){
#ifdef PASTIX
#endif
    }

    SFREE(aubd);SFREE(adbd);SFREE(irowd);SFREE(jqd);SFREE(rhs);
    SFREE(vector);SFREE(icold);SFREE(aud);SFREE(add);
    SFREE(irowc);SFREE(jqc);SFREE(auc);SFREE(icolc);

  }
   
  /* ************************************************************ */
  /* end of conditional random field computation                  */
  /* ************************************************************ */

  
  /* assign a starting value sigma for the calculation of the eigenvalues
     sigma is set to trace(autocovmatrix) = maximum eigenvalue */
  
  for(i=0;i<*ndesi;i++){
    trace+=ad[i];
  }

  printf(" Trace of the covariance matrix: %e\n\n", trace);
  /* Calculate the eigenmodes and eigenvalues of the (conditional) 
     autocovariance matrix */

  NNEW(adb,double,*ndesi);
  NNEW(aub,double,nzss); 
  
  for(i=0;i<*ndesi;i++){
    adb[i]=1.0;
  }
    
  /* LU decomposition of the left hand matrix */

  printf(" LU decomposition of the covariance matrix\n\n");

  if(*isolver==0){
#ifdef SPOOLES
    spooles_factor(ad,au,adb,aub,&trace,icols,irows,ndesi,&nzss,
		   &symmetryflag,&inputformat,&nzss);
#else
    printf(" *ERROR in randomfieldmain: the SPOOLES library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==4){
#ifdef SGI
    token=1;
    sgi_factor(ad,au,adb,aub,&trace,icols,irows,ndesi,&nzss,token);
#else
    printf(" *ERROR in randomfieldmain: the SGI library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==5){
#ifdef TAUCS
    tau_factor(ad,&au,adb,aub,&trace,icols,&irows,ndesi,&nzss);
#else
    printf(" *ERROR in randomfieldmain: the TAUCS library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  /*  else if(*isolver==6){
#ifdef MATRIXSTORAGE
    matrixstorage(ad,&au,adb,aub,&trace,icols,&irows,ndesi,&nzss,
		  ntrans,inotr,trab,co,nk,nactdof,jobnamec,mi,ipkon,
		  lakon,kon,ne,mei,nboun,nmpc,cs,mcs,&ithermal,nmethod);
#else
    printf(" *ERROR in randomfieldmain: the MATRIXSTORAGE library is not linked\n\n");
    FORTRAN(stop,());
#endif
}*/
  else if(*isolver==7){
#ifdef PARDISO
    pardiso_factor(ad,au,adb,aub,&trace,icols,irows,ndesi,&nzss,
		   &symmetryflag,&inputformat,jqs,&nzss);
#else
    printf(" *ERROR in randomfieldmain: the PARDISO library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==8){
#ifdef PASTIX
#ifdef PARDISO
    pardiso_factor(ad,au,adb,aub,&trace,icols,irows,ndesi,&nzss,
		   &symmetryflag,&inputformat,jqs,&nzss);
#else
    pastix_factor_main(ad,au,adb,aub,&trace,icols,irows,ndesi,&nzss,
		  &symmetryflag,&inputformat,jqs,&nzss);
#endif
#else
    printf(" *ERROR in randomfieldmain: the PASTIX library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  
  /* calculating the eigenvalues and eigenmodes */
  
  /* initialization of relative error measure and number of eigenvalues*/
    
  actreliability=0.;
  nev=1;
  ncv=4*nev;
  
  NNEW(resid,double,*ndesi);
  NNEW(zz,double,(long long)ncv**ndesi);
  NNEW(workd,double,3**ndesi);
  NNEW(workl,double,48);
  NNEW(temp_array,double,*ndesi);
  NNEW(select,ITG,4); 
  NNEW(d,double,nev);
  
  while(actreliability<=reliability){

    niter+=1;
    
    if(nev>=*ndesi){
      printf(" *WARNING in randomfieldmain: number of requested\n");
      printf("          eigenvectors must be smaller than the \n");
      printf("          number of design variables; \n");
      printf("          the number of eigenvectors is set to \n");
      printf("          the number of design variables minus one. \n");
      nev=*ndesi-1;
      break;
    }

    printf(" Calculating the eigenvalues and the eigenmodes\n\n");
     
    ncv=4*nev;
    ido=0;
    ldz=*ndesi;
    iparam[0]=1;
    iparam[2]=mxiter;
    iparam[3]=1;
    iparam[6]=3;
    info=0;

    RENEW(zz,double,(long long)ncv**ndesi);
     
    lworkl=ncv*(8+ncv);
    RENEW(workl,double,lworkl);
    FORTRAN(dsaupd,(&ido,bmat,ndesi,which,&nev,&tol,resid,&ncv,zz,&ldz,
		    iparam,ipntr,workd,workl,&lworkl,&info));

    while((ido==-1)||(ido==1)||(ido==2)){
      if(ido==-1){
	FORTRAN(op,(ndesi,&workd[ipntr[0]-1],temp_array,adb,aub,jqs,irows));
      }
      if((ido==-1)||(ido==1)){

	/* solve the linear equation system  */

	if(ido==-1){
	  if(*isolver==0){
#ifdef SPOOLES
	    spooles_solve(temp_array,ndesi);
#endif
	  }
	  else if(*isolver==4){
#ifdef SGI
	    sgi_solve(temp_array,token);
#endif
	  }
	  else if(*isolver==5){
#ifdef TAUCS
	    tau_solve(temp_array,ndesi);
#endif
	  }
	  else if(*isolver==7){
#ifdef PARDISO
	    pardiso_solve(temp_array,ndesi,&symmetryflag,&inputformat,&nrhs);
#endif
	  }
	  else if(*isolver==8){
#ifdef PASTIX
#ifdef PARDISO
	    pardiso_solve(temp_array,ndesi,&symmetryflag,&inputformat,&nrhs);
#else
	    pastix_solve(temp_array,ndesi,&symmetryflag,&nrhs);
#endif
#endif
	  }
	   
	  for(jrow=0;jrow<*ndesi;jrow++){
	    workd[ipntr[1]-1+jrow]=temp_array[jrow];
	  }
	}
	else if(ido==1){
	  if(*isolver==0){
#ifdef SPOOLES
	    spooles_solve(&workd[ipntr[2]-1],ndesi);
#endif
	  }
	  else if(*isolver==4){
#ifdef SGI
	    sgi_solve(&workd[ipntr[2]-1],token);
#endif
	  }
	  else if(*isolver==5){
#ifdef TAUCS
	    tau_solve(&workd[ipntr[2]-1],ndesi);
#endif
	  }
	  else if(*isolver==7){
#ifdef PARDISO
	    pardiso_solve(&workd[ipntr[2]-1],ndesi,&symmetryflag,&inputformat,&nrhs);
#endif
	  }
	  else if(*isolver==8){
#ifdef PASTIX
#ifdef PARDISO
	    pardiso_solve(&workd[ipntr[2]-1],ndesi,&symmetryflag,&inputformat,&nrhs);
#else
	    pastix_solve(&workd[ipntr[2]-1],ndesi,&symmetryflag,&nrhs);
#endif
#endif
	  }
	   
	  for(jrow=0;jrow<*ndesi;jrow++){
	    workd[ipntr[1]-1+jrow]=workd[ipntr[2]-1+jrow];
	  }
	}

      }

      if(ido==2){
	FORTRAN(op,(ndesi,&workd[ipntr[0]-1],&workd[ipntr[1]-1],
		    adb,aub,jqs,irows));
      }

      FORTRAN(dsaupd,(&ido,bmat,ndesi,which,&nev,&tol,resid,&ncv,zz,&ldz,
		      iparam,ipntr,workd,workl,&lworkl,&info));
    }
    
    RENEW(select,ITG,ncv); 
    RENEW(d,double,nev);
     
    FORTRAN(dseupd,(&rvec,howmny,select,d,zz,&ldz,&trace,bmat,ndesi,which,
		    &nev,&tol,resid,&ncv,zz,&ldz,iparam,ipntr,workd,workl,
		    &lworkl,&info));
  
    /* calculating the error measure and estimation of number of 
       eigenvalues for next iteration */

    for(i=nev-nevold-1;i>=0;i--){
      imodes+=1;
      dsum+=d[i];
      actreliability=dsum/trace;
      FORTRAN(writerandomfield,(&d[i],&actreliability,&imodes));
    }

    printf("\n Reliability iteration: %" ITGFORMAT "\n", niter);
    printf(" Number of eigenvalues calculated: %" ITGFORMAT "\n", nev);
    printf(" Corresponding reliability: %e\n\n", actreliability);

    nevold=nev;
    delta=reliability-actreliability;
    nev=(ITG)(nev+delta/(d[0]/trace)+1);

  } /*while loop end*/

  /* storing the mean in the frd file */

  NNEW(acvector,double,3**nk);
  NNEW(acscalar,double,*nk);
  
  ++*kode;
  for(idesvar=0;idesvar<*ndesi;idesvar++){
     node=nodedesi[idesvar]-1;
     if(irobustdesign[1]==1){
  	acscalar[node]=randomval[0];
	for(i=0;i<3;i++){
	   acvector[3*node+i]=randomval[0]*extnor[3*node+i];
	}
     }else{
        acscalar[node]=randomval[2*node];
	for(i=0;i<3;i++){
	   acvector[3*node+i]=randomval[2*node]*extnor[3*node+i];
	}
     }
  }

  irand=3;
  frd_sen(co,nk,stn,inum,nmethod,kode,filab,&time,nstate_,istep,
    	  &iinc,&k,&noddiam,description,mi,&ngraph,ne,cs,set,nset,
    	  istartset,iendset,ialset,jobnamec,output,
    	  acscalar,&iobject,objectset,ntrans,inotr,trab,&idesvar,orname,
    	  &icoordinate,&inorm,&irand,&ishape,&ifeasd);
  irand=4;
  frd_sen(co,nk,stn,inum,nmethod,kode,filab,&time,nstate_,istep,
    	  &iinc,&k,&noddiam,description,mi,&ngraph,ne,cs,set,nset,
    	  istartset,iendset,ialset,jobnamec,output,
    	  acvector,&iobject,objectset,ntrans,inotr,trab,&idesvar,orname,
    	  &icoordinate,&inorm,&irand,&ishape,&ifeasd); 
  irand=1;
      
  /* storing the eigenvectors of the autocovariance matrix in the frd file 
     acvector = autocovariance vector */
  
  for(k=nev;k>=1;k--){
      
    /* generating and writing the autocorrelation vector */
      
    ++*kode;
    for(idesvar=0;idesvar<*ndesi;idesvar++){
      node=nodedesi[idesvar]-1;
      acscalar[node]=zz[(k-1)**ndesi+idesvar];
      for(i=0;i<3;i++){
	acvector[3*node+i]=zz[(k-1)**ndesi+idesvar]*extnor[3*node+i];
      }
    }
      
    frd_sen(co,nk,stn,inum,nmethod,kode,filab,&d[k-1],nstate_,
	    istep,
	    &iinc,&k,&noddiam,description,mi,&ngraph,ne,cs,set,nset,
	    istartset,iendset,ialset,jobnamec,output,
	    acscalar,&iobject,objectset,ntrans,inotr,trab,&idesvar,orname,
	    &icoordinate,&inorm,&irand,&ishape,&ifeasd);
    irand=2;
    frd_sen(co,nk,stn,inum,nmethod,kode,filab,&d[k-1],nstate_,
	    istep,
	    &iinc,&k,&noddiam,description,mi,&ngraph,ne,cs,set,nset,
	    istartset,iendset,ialset,jobnamec,output,
	    acvector,&iobject,objectset,ntrans,inotr,trab,&idesvar,orname,
	    &icoordinate,&inorm,&irand,&ishape,&ifeasd); 
    irand=1;
     
  }

  /* generating and writing the local reliability */
    
  ++*kode;
  for(idesvar=0;idesvar<*ndesi;idesvar++){
    node=nodedesi[idesvar]-1;
    acscalar[node]=0.0;
  }
  for(k=nev;k>=1;k--){
    for(idesvar=0;idesvar<*ndesi;idesvar++){
      node=nodedesi[idesvar]-1;
      acscalar[node]+=zz[(k-1)**ndesi+idesvar]*zz[(k-1)**ndesi+idesvar]*
	d[k-1];
    }
  }
  
  for(idesvar=0;idesvar<*ndesi;idesvar++){
    node=nodedesi[idesvar]-1;
    if(irobustdesign[1]==1){
      acscalar[node]=acscalar[node]/(randomval[1]*randomval[1]);
    }else{
	acscalar[node]=acscalar[node]/(randomval[2*(node)+1]*
				       randomval[2*(node)+1]);
    }
  }

  irand=5;
  frd_sen(co,nk,stn,inum,nmethod,kode,filab,&time,nstate_,
	  istep,
	  &iinc,&k,&noddiam,description,mi,&ngraph,ne,cs,set,nset,
	  istartset,iendset,ialset,jobnamec,output,
	  acscalar,&iobject,objectset,ntrans,inotr,trab,&idesvar,orname,
	  &icoordinate,&inorm,&irand,&ishape,&ifeasd);
  /*
    -----------
    free memory
    -----------
  */
    
  SFREE(temp_array);
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
#ifdef PARDISO
    pardiso_cleanup(ndesi,&symmetryflag,&inputformat);
#endif
#endif
  }

  if(info!=0){
    printf(" *ERROR in randomfieldmain: info=%" ITGFORMAT "\n",info);
    printf("        # of converged eigenvalues=%" ITGFORMAT "\n\n",iparam[4]);
  }         
      
  SFREE(select);SFREE(workd);SFREE(workl);SFREE(resid);
  SFREE(irows);SFREE(icols);SFREE(jqs);SFREE(au);SFREE(ad);
  SFREE(adb);SFREE(aub);SFREE(d);SFREE(zz);SFREE(randval);
  SFREE(randfield);SFREE(acvector);SFREE(acscalar);
       
  return;
  
}

#endif
