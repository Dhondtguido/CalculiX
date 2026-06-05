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

#include <unistd.h>
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

void compfluidfem(double **cop,ITG *nk,ITG **ipkonp,ITG **konp,char **lakonp,
		  ITG *ne,char **sidefacep,ITG *ifreestream,
		  ITG *nfreestream,ITG *isolidsurf,ITG *neighsolidsurf,
		  ITG *nsolidsurf,ITG *iponoel,ITG *inoel,ITG *nshcon,
		  double *shcon,
		  ITG *nrhcon,double *rhcon,double **voldp,ITG *ntmat_,
		  ITG *nodeboun,
		  ITG *ndirboun,ITG *nboun,ITG **ipompcp,ITG **nodempcp,
		  ITG *nmpc,
		  ITG **ikmpcp,ITG **ilmpcp,ITG *ithermal,ITG *ikboun,
		  ITG *ilboun,
		  ITG *iturbulent,ITG *isolver,ITG *iexpl,
		  double *ttime,
		  double *time,double *dtime,ITG *nodeforc,ITG *ndirforc,
		  double *xforc,
		  ITG *nforc,ITG *nelemload,char *sideload,double *xload,
		  ITG *nload,
		  double *xbody,ITG *ipobody,ITG *nbody,ITG **ielmatp,
		  char *matname,
		  ITG *mi,ITG *ncmat_,double *physcon,ITG *istep,ITG *iinc,
		  ITG *ibody,double *xloadold,double *xboun,
		  double **coefmpcp,ITG *nmethod,double *xforcold,
		  double *xforcact,
		  ITG *iamforc,ITG *iamload,double *xbodyold,double *xbodyact,
		  double *t1old,double *t1,double *t1act,ITG *iamt1,
		  double *amta,
		  ITG *namta,ITG *nam,double *ampli,double *xbounold,
		  double *xbounact,
		  ITG *iamboun,ITG *itg,ITG *ntg,char *amname,double *t0,
		  ITG **nelemfacep,
		  ITG *nface,double *cocon,ITG *ncocon,double *xloadact,
		  double *tper,
		  ITG *jmax,ITG *jout,char *set,ITG *nset,ITG *istartset,
		  ITG *iendset,ITG *ialset,char *prset,char *prlab,ITG *nprint,
		  double *trab,ITG *inotr,ITG *ntrans,char *filab,
		  char **labmpcp,
		  double *sti,ITG *norien,double *orab,char *jobnamef,
		  char *tieset,
		  ITG *ntie,ITG *mcs,ITG *ics,double *cs,ITG *nkon,ITG *mpcfree,
		  ITG *memmpc_,double **fmpcp,ITG *nef,ITG **inomatp,
		  double *qfx,ITG *kode,ITG *ipface,ITG *ielprop,double *prop,
		  char *orname,double *tincf,ITG *ifreesurface,ITG *nktot,
		  ITG *ielorien,ITG *nelold,ITG *nkold,ITG *nknew,ITG *nelnew){

  /* main computational fluid dynamics routine */

  /* References:

     Zienkiewicz, O.C., Taylor, R.L. and Nithiarasu, P., "The Finite
     Element Method for Fluid Dynamics", 6th Edition, Elsevier (2006)

     Menter, F.R., "Two-Equation Eddy-Viscosity Turbulence Models
     for Engineering Applications", AIAA Journal(1994), 32(8), 
     1598-1605                                                       */
  
  char cflag[1],*labmpc=NULL,*lakon=NULL,*sideface=NULL;

  ITG *ipointer=NULL,*mast1=NULL,*irowv=NULL,*irowp=NULL,
    *icolv=NULL,*icolp=NULL,*jqv=NULL,*jqp=NULL,*nactdoh=NULL,i,j,
    *nx=NULL,*ny=NULL,*nz=NULL,nzs,neqp,nzsv,nzsp,iexplicit,nzlp,nnstep,
    iconvergence=0,iit,symmetryflag=0,inputformat=0,compressible,
    *inum=NULL,iqfx=0,isti=0,nrhs=1,
    mt=mi[1]+1,*ipvar=NULL,*ipvarf=NULL,nvar_,nvarf_,
    nfield,ndim,iorienloc,maxit=1,
    *integerglob=NULL,*nodempc=NULL,*ipompc=NULL,*ikmpc=NULL,*ilmpc=NULL,
    *ipkon=NULL,*kon=NULL,*ielmat=NULL,*nelemface=NULL,
    *inomat=NULL,ithermalref,ipower=1,*jyy=NULL,*jqvr=NULL,
    *irowvr=NULL,ierr=0,im,*nodfreesurf=NULL,istart,num_cpus;

  double *yy=NULL,*xsolidsurf=NULL,*dt=NULL,*vcon=NULL,*x=NULL,*y=NULL,*z=NULL,
    *xo=NULL,*yo=NULL,*zo=NULL,*adbv=NULL,*aubv=NULL,
    *adbp=NULL,*aubp=NULL,timef,ttimef,ttimf,xg[3],
    dtimef,*sol=NULL,*aux=NULL,*stn=NULL,theta1,theta2,*adb=NULL,
    *qfn=NULL,*aub=NULL,sigma=0.,*dh=NULL,reltimef,*thicke=NULL,
    shockcoef,*sa=NULL,*varf=NULL,dgravity,
    *adlv=NULL,*adlp=NULL,*var=NULL,
    *vconini=NULL,*doubleglob=NULL,*coefmpc=NULL,*fmpc=NULL,*co=NULL,
    *vold=NULL,reltime,*dhel=NULL,*voldo=NULL,dtimefo,temp,ratio,
    *aubvr=NULL,*coefmodmpc=NULL,*voldini=NULL,*depth=NULL,*vcono=NULL,
    vconmax[7],vmax[7],*v=NULL,*b1=NULL,*b2=NULL;

  nodempc=*nodempcp;ipompc=*ipompcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
  coefmpc=*coefmpcp;labmpc=*labmpcp;fmpc=*fmpcp;co=*cop;
  ipkon=*ipkonp;lakon=*lakonp;kon=*konp;ielmat=*ielmatp;
  nelemface=*nelemfacep;sideface=*sidefacep;
  vold=*voldp;inomat=*inomatp;

  /* standard: shockcoef=0 */

#ifdef SGI
  ITG token;
#endif

  /* open frd-file for fluids */

  FORTRAN(openfilefluidfem,(jobnamef));

  /* variables for multithreading procedure */

  ITG sys_cpus;
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
  
  envloc = getenv("CCX_NPROC_CFD");
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
  
  // next line is to be inserted in a similar way for all other parallel parts
  
  if(*ne<num_cpus) num_cpus=*ne;
  
  printf(" Using up to %" ITGFORMAT " cpu(s) for CFD.\n", num_cpus);
  
  *kode=0;
  
  /*  *iexpl==0:  structure:implicit, fluid:semi-implicit
   *iexpl==1:  structure:implicit, fluid:compressible explicit
   *iexpl==2:  structure:explicit, fluid:incompressible semi-implicit
   *iexpl==3:  structure:explicit, fluid:compressible semi-implicit  */

  if(*iexpl==1){

    /* compressible explicit */
    
    iexplicit=1;
    theta1=0.5;theta2=0.;compressible=1;
    if(*tincf<=1.25){*tincf=1.25;}
  }else{

    /* incompressible implicit */
    
    iexplicit=0;
    theta1=1.0;theta2=1.0;compressible=0;
    if(*tincf<=1.){*tincf=1.;}
  }

  /* if initial conditions are specified for the temperature, 
     it is assumed that the temperature is an unknown */

  ithermalref=ithermal[0];
  if(ithermal[0]==1){
    ithermal[0]=2;
  }

  /* check whether energy equation has to be included in the
     calculations */
  
  if(ithermal[0]>1){
    istart=0;
  }else{
    istart=1;
  }

  if(*mcs>1){
    printf(" *ERROR in compfluid: for CFD only one cyclic symmetry\n");
    printf("        conditions is allowed\n");
    FORTRAN(stop,());
  }

  /* determining the matrix structure */

  nzs=1000000;
  
  NNEW(ipointer,ITG,3**nk);
  NNEW(mast1,ITG,nzs);
  NNEW(irowv,ITG,nzs);
  NNEW(irowp,ITG,nzs);
  NNEW(icolv,ITG,*nk);
  NNEW(icolp,ITG,*nk);
  NNEW(jqv,ITG,*nk+1);
  NNEW(jqp,ITG,*nk+1);
  NNEW(nactdoh,ITG,*nk);

  mastructffem(nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,nboun,ipompc,
	       nodempc,nmpc,nactdoh,icolv,icolp,jqv,jqp,
	       &mast1,&irowv,&irowp,&neqp,ipointer,&nzsv,&nzsp,
	       &nzs,&compressible,inomat);

  SFREE(ipointer);SFREE(mast1);

  /* initialization */

  NNEW(yy,double,*nk);
  NNEW(jyy,ITG,*nk);
  NNEW(xsolidsurf,double,*nsolidsurf);
  NNEW(dh,double,*nk);
  NNEW(vcon,double,mt**nk);
  NNEW(vcono,double,mt**nk);
  NNEW(x,double,*nsolidsurf);
  NNEW(y,double,*nsolidsurf);
  NNEW(z,double,*nsolidsurf);
  NNEW(xo,double,*nsolidsurf);
  NNEW(yo,double,*nsolidsurf);
  NNEW(zo,double,*nsolidsurf);
  NNEW(nx,ITG,*nsolidsurf);
  NNEW(ny,ITG,*nsolidsurf);
  NNEW(nz,ITG,*nsolidsurf);
  NNEW(dhel,double,*ne);
  if(*ifreesurface==1){
    NNEW(depth,double,*nk);
    NNEW(nodfreesurf,ITG,*nk);
  }

  /* perform initial calculations such as distances needed for
     turbulence calculations and initial values for the conservative
     variables (including the density) */
  
  FORTRAN(initialcfdfem,(yy,nk,co,ne,ipkon,kon,lakon,x,y,z,xo,yo,zo,nx,ny,nz,
			 isolidsurf,neighsolidsurf,xsolidsurf,dh,nshcon,shcon,
			 nrhcon,rhcon,vold,vcon,ntmat_,iponoel,inoel,
			 nsolidsurf,iturbulent,physcon,
			 &compressible,matname,inomat,mi,
			 ithermal,dhel,jyy,ifreesurface,nbody,ipobody,
			 ibody,xbody,depth,nodfreesurf,&dgravity,xg));
  
  SFREE(x);SFREE(y);SFREE(z);SFREE(xo);SFREE(yo);SFREE(zo);SFREE(nx);SFREE(ny);
  SFREE(nz);SFREE(dhel);

  if(compressible){
    NNEW(coefmodmpc,double,*memmpc_);
    FORTRAN(normmpc,(nmpc,ipompc,nodempc,coefmpc,inomat,coefmodmpc,ikboun,
		     nboun));
  }
  
  /* calculating the shape functions, their derivatives and the
     Jacobian determinant in the integration points of the elements */

  nvar_=35**ne;
  NNEW(ipvar,ITG,*ne);
  NNEW(var,double,nvar_);

  nvarf_=8**ne;
  NNEW(ipvarf,ITG,*ne);
  NNEW(varf,double,nvarf_);

  calcshapef(&nvar_,ipvar,&var,ne,lakon,co,ipkon,kon,
             nelemface,sideface,nface,&nvarf_,ipvarf,
             &varf,iturbulent,yy);

  /* composing those left hand sides which do not depend on the increment */

  /* lhs for the energy, velocity, pressure for compressible fluids
     and turbulence */
  
  NNEW(adbv,double,*nk);
  NNEW(aubv,double,nzsv);

  FORTRAN(mafillvlhs,(nk,kon,ipkon,lakon,ne,icolv,jqv,irowv,&nzsv,adbv,aubv,
		      ipvar,var));

  NNEW(adlv,double,*nk);
  FORTRAN(lump,(adbv,aubv,adlv,irowv,jqv,nk));

  if(compressible){
    convert2rowbyrow(adbv,aubv,icolv,irowv,jqv,nk,&nzsv,&aubvr,&jqvr,
		     &irowvr);
  }

  /* lhs for the pressure for incompressible fluids  */

  if(compressible==0){

    NNEW(adbp,double,neqp);
    NNEW(aubp,double,nzsp);

    FORTRAN(mafillplhs,(kon,ipkon,lakon,ne,ipompc,nodempc,coefmpc,nmpc,nactdoh,
			icolp,jqp,irowp,&neqp,&nzlp,&nzsp,adbp,aubp,ipvar,var));

    if(neqp>0){

      /* LU decomposition of the left hand matrix */

      if(*isolver==0){
#ifdef SPOOLES
	spooles_factor(adbp,aubp,adb,aub,&sigma,icolp,irowp,&neqp,&nzsp,
		       &symmetryflag,&inputformat,&nzsp);
#else
	printf(" *ERROR in compfluid: the SPOOLES library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==4){
#ifdef SGI
	token=1;
	sgi_factor(adbp,aubp,adb,aub,&sigma,icolp,irowp,&neqp,&nzsp,token);
#else
	printf(" *ERROR in compfluid: the SGI library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
	tau_factor(adbp,&aubp,adb,aub,&sigma,icolp,&irowp,&neqp,&nzsp);
#else
	printf(" *ERROR in compfluid: the TAUCS library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	pardiso_factor(adbp,aubp,adb,aub,&sigma,icolp,irowp,&neqp,&nzsp,
		       &symmetryflag,&inputformat,jqp,&nzsp);
#else
	printf(" *ERROR in compfluid: the PARDISO library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==8){
#ifdef PASTIX
	pastix_factor_main(adbp,aubp,adb,aub,&sigma,icolp,irowp,&neqp,&nzsp,
			   &symmetryflag,&inputformat,jqp,&nzsp);
#else
	printf(" *ERROR in linstatic: the PASTIX library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      
    }
  }

  /* starting the main loop */

  NNEW(v,double,mt**nk);
      
  /* inserting the velocity and temperature conditions 
     for incompressible materials*/
    
  if(compressible==0){

      /* inserting SPC's/MPC's to the physical variables */
      FORTRAN(applybounfem,(nodeboun,ndirboun,xbounact,nk,
                            vold,isolidsurf,xsolidsurf,ifreestream,iturbulent,vcon,
                            shcon,nshcon,ntmat_,physcon,v,
                            &compressible,nodempc,ipompc,coefmpc,
                            inomat,mi,ilboun,ilmpc,labmpc,
                            coefmodmpc,&iexplicit,nboun,
                            nmpc,nfreestream,nsolidsurf));

      /* calculating the conservative variables from the physical
	 variables */

      FORTRAN(phys2con,(inomat,vold,ntmat_,shcon,nshcon,
		    physcon,&compressible,vcon,rhcon,nrhcon,ithermal,
		    mi,ifreesurface,&ierr,&dgravity,depth,nk,&num_cpus));

  }

  /* initialization of voldo: solution from last fluid increment
     and dtimefo: last time step */
  
  NNEW(voldo,double,mt**nk);
  memcpy(voldo,vold,sizeof(double)*mt**nk);
  dtimefo=1.;
  
  /* ttime is the total time up to the start of the mechanical increment
     time is the step time up to the end of the mechanical increment 
     dtime is the present increment size */

  reltime=*time/(*tper);

  ttimef=*ttime;
  timef=*time-*dtime;
  NNEW(dt,double,*nk);

  if(iexplicit){
    NNEW(sa,double,*nk);
    shockcoef=physcon[13];
  }

  /* for compressible fluids: store the initial values of vold and
     vcon in case the calculation has to be restarted with a higher
     shock smoothing */
  
  if(iexplicit){
    NNEW(vconini,double,mt**nk);
    memcpy(vconini,vcon,sizeof(double)*mt**nk);
    NNEW(voldini,double,mt**nk);
    memcpy(voldini,vold,sizeof(double)*mt**nk);
  }
  
  do{

    /* this loop is needed in case the shock smoothing coefficient
       has to be increased for compressible flow */

    iit=0;

    do{

      iit++;

      /* determining a new time increment */

      if((*nmethod==4)||((*nmethod==1)&&((iit/ipower)*ipower==iit))){
	FORTRAN(compdt,(nk,dt,nshcon,shcon,vold,ntmat_,iponoel,
			inoel,&dtimef,ielmat,dh,cocon,ncocon,
			ithermal,mi,
			vcon,&compressible,tincf,&ierr,ifreesurface,
			&dgravity,&iit));

	/* check whether an error occurred: if so, increase the shock
	   coefficient (only for compressible fluids) */

	if(ierr==1){
	  if(iexplicit==0) FORTRAN(stop,());
	  ierr=0;
	  ipower=1;
	  dtimefo=1.;
	  ttimef=*ttime;
	  timef=*time-*dtime;
	  if(shockcoef==0.){shockcoef=0.001;}else{
	    shockcoef*=2;
	    if(shockcoef>2.){
	      printf("shock coefficient > 2.; stop\n");
	      FORTRAN(stop,());
	    }
	  }
	  memcpy(vcon,vconini,sizeof(double)*mt**nk);
	  memcpy(vold,voldini,sizeof(double)*mt**nk);
	  memcpy(voldo,voldini,sizeof(double)*mt**nk);
	  DMEMSET(sa,0,*nk,0.);
	  break;
	}
	
	if(*nmethod==1) ipower*=2;
	if(iexplicit){
	  printf("iteration %" ITGFORMAT " shock coefficient %e\n",iit,shockcoef);
	}else{
	  printf("iteration %" ITGFORMAT "\n",iit);
	}
      }

      timef+=dtimef;
      if((*time<timef)&&(*nmethod==4)){
	dtimef-=timef-*time;
	if(dtimef<1.e-6**dtime){
	  iconvergence=1;
	  break;
	}
	timef=*time;
	iconvergence=1;
      }
      reltimef=timef/(*tper);
      if(reltimef>1.) reltimef=1.;
      
      /* predicting the new physical variables by linear 
         extrapolation (only for explicit compressible calculations) */

      if (compressible==1) {
        memcpy(vcono,vcon,sizeof(double)*mt**nk);
        ratio = dtimef / (2.*dtimefo);
        #pragma omp parallel for private(j, temp)
        for (i=0;i<*nk;i++) {
          for (j=0;j<mt;j++) {
            temp = vold[i*mt+j];
            vold[i*mt+j] += (vold[i*mt+j]-voldo[i*mt+j])*ratio;
            voldo[i*mt+j] = temp;
          }
        }
        dtimefo = dtimef;
      }

      /* determining the instantaneous load */

      if(*nmethod==1){

	/* TO DO: check whether ialset for submodels is correct, or
           whether the renumbering of nodes and elements in CFD
           has consequences here */
	
	/* boundary conditions at end of mechanical increment */

	if(iit==1){
	  FORTRAN(temploadfem,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,
			       xload,xloadact,iamload,nload,ibody,xbody,nbody,
			       xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
			       namta,nam,ampli,time,&reltime,ttime,dtime,
			       ithermal,
			       nmethod,xbounold,xboun,xbounact,iamboun,nboun,
			       nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
			       co,
			       vold,itg,ntg,amname,ikboun,ilboun,nelemload,
			       sideload,mi,ntrans,trab,inotr,vold,integerglob,
			       doubleglob,tieset,istartset,iendset,ialset,ntie,
			       nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc,set,
			       nset,cocon,ncocon,rhcon,nrhcon,shcon,nshcon,
			       ielmat,ielprop,prop,iponoel,inoel,ipkon,kon,
			       lakon,ipobody,ntmat_));
	}
	  
      }else if(*nmethod==4){

	/* boundary conditions at end of fluid increment */

	FORTRAN(temploadfem,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,
			     xload,xloadact,iamload,nload,ibody,xbody,nbody,
			     xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
			     namta,nam,ampli,&timef,&reltimef,&ttimef,&dtimef,
			     ithermal,nmethod,xbounold,xboun,xbounact,iamboun,
			     nboun,nodeboun,ndirboun,nodeforc,ndirforc,istep,
			     iinc,co,vold,itg,ntg,amname,ikboun,ilboun,
			     nelemload,
			     sideload,mi,ntrans,trab,inotr,vold,integerglob,
			     doubleglob,tieset,istartset,iendset,ialset,ntie,
			     nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc,set,nset,
			     cocon,ncocon,rhcon,nrhcon,shcon,nshcon,
			     ielmat,ielprop,prop,iponoel,inoel,ipkon,kon,
			     lakon,ipobody,ntmat_));
      }

      /* STEP 1: velocity correction */

      /* b1 contains the rhs of all steps except step 3;
         first the energy rhs is stored for all nodes, then the
         momentum rhs in x for all nodes etc.... */
      
      NNEW(b1,double,mt**nk);
      NNEW(b2,double,3**nk);

      FORTRAN(mafillv1rhs,(co,nk,kon,ipkon,lakon,ne,
                           ipompc,nodempc,coefmpc,nmpc,
                           nelemload,sideload,xloadact,nload,xbodyact,
                           ipobody,nbody,b2,nactdoh,nmethod,ikmpc,ilmpc,
                           rhcon,nrhcon,ielmat,ntmat_,ithermal,
                           vold,vcon,mi,physcon,shcon,nshcon,ttime,
                           &timef,istep,ibody,xloadold,iturbulent,
                           nelemface,sideface,nface,&compressible,
                           &dtimef,ipvar,var,ipvarf,varf,ipface,ifreesurface,
                           depth,&dgravity,cocon,ncocon,iinc,&theta1,
                           &reltimef,b1,&num_cpus));

      /* solving the equations */
      
      NNEW(aux,double,*nk);
      for(i=istart;i<4;i++){
	solveeq(adbv,aubv,adlv,&b1[i**nk],&v[i**nk],aux,
		irowv,jqv,nk,&maxit,&num_cpus);
      }

      for(i=5;i<mt;i++){
	solveeq(adbv,aubv,adlv,&b1[i**nk],&v[i**nk],aux,
			 irowv,jqv,nk,&maxit,&num_cpus);
      }
      SFREE(aux);

      /* STEP 2: pressure correction */

      FORTRAN(mafillprhs,(nk,kon,ipkon,lakon,ipompc,nodempc,coefmpc,
                          nmpc,b1,nactdoh,mi,v,&theta1,ne,
                          &dtimef,ipvar,var,&compressible,&num_cpus));

      if(compressible){
	      
	/* explicit compressible calculations: solve */

	NNEW(aux,double,neqp);
	solveeq(adbv,aubv,adlv,&b1[4**nk],&v[4**nk],aux,
			 irowv,jqv,nk,&maxit,&num_cpus);
	SFREE(aux);
	
      }else if(neqp>0){

	/* solving the system of equations (only for liquids) */

	if(*isolver==0){
#ifdef SPOOLES
	  spooles_solve(&b1[4**nk],&neqp);
#endif
	}
	else if(*isolver==4){
#ifdef SGI
	  sgi_solve(&b1[4**nk],token);
#endif
	}
	else if(*isolver==5){
#ifdef TAUCS
	  tau_solve(&b1[4**nk],&neqp);
#endif
	}
	else if(*isolver==7){
#ifdef PARDISO
	  pardiso_solve(&b1[4**nk],&neqp,&symmetryflag,&inputformat,&nrhs);
#endif
	}
	else if(*isolver==8){
#ifdef PASTIX
	  pastix_solve(&b1[4**nk],&neqp,&symmetryflag,&nrhs);
#endif
	}

	/* copying the solution into field sol */

	#pragma omp parallel for num_threads(num_cpus)
	for(j=0;j<*nk;j++){
	  b1[4**nk+j]=b1[4**nk+j]/(theta1*theta2*dtimef*dtimef);
	}

	/* storing the pressure (incompressible) or density (compressible)
	   correction in v */

	FORTRAN(resultsp,(nk,nactdoh,v,&b1[4**nk],mi));
	
      }

      SFREE(b1);

      if(compressible==0){
      
	/* inserting the pressure boundary conditions for liquids */

	FORTRAN(applybounp,(nodeboun,ndirboun,nboun,xbounact,nk,vold,v,
			    ipompc,nodempc,coefmpc,nmpc,inomat,mi));
      }
      
      /* STEP 3: velocity correction */


      if(compressible==0){
        FORTRAN(mafillv2rhs,(kon,ipkon,lakon,b2,v,ne,mi,
                             &dtimef,ipvar,var,nk,&num_cpus));
      }

      /* solving the equations for V** */
      
      NNEW(sol,double,*nk);
      NNEW(aux,double,*nk);
      for(j=0;j<3;j++){
        solveeq(adbv,aubv,adlv,&b2[j**nk],sol,aux,irowv,
                jqv,nk,&maxit,&num_cpus);
        #pragma omp parallel for num_threads(num_cpus)
        for(i=0;i<*nk;i++){
          v[(j+1)**nk+i]+=sol[i];
        }
      }
      SFREE(b2);SFREE(aux);SFREE(sol);
      
      /* update the conservative variables
	 (for incompressible fluids: pressure instead of density  */

      FORTRAN(updatecon,(vold,vcon,v,nk,ithermal,iturbulent,mi,
		     &compressible,&num_cpus));

      /* shock smoothing the solution (only for compressible fluids in the
         transsonic and supersonic range) */

      if((compressible)&&(shockcoef>0.)){

        /* shocksmoothing rho * total energy density, rho * velocity and rho */

        NNEW(aux,double,*nk);
        #pragma omp parallel private(j) shared(aux) num_threads(num_cpus)
        {
          for(j=istart;j<5;j++){
            #pragma omp for
            for(i=0;i<*nk;i++){
              aux[i]=0;
              for(ITG k=jqvr[i]-1;k<jqvr[i+1]-1;k++){
                aux[i]+=aubvr[k]*vcon[j**nk+irowvr[k]-1];
              }
            }
            #pragma omp for
            for(i=0;i<*nk;i++){
              vcon[j**nk+i]+=sa[i]*aux[i]*adlv[i];
            }
          }
        }
        SFREE(aux);

      }

      /* deriving the physical variables from the conservative
	 variables */

      FORTRAN(con2phys,(vold,vcon,nk,ntmat_,shcon,nshcon,rhcon,nrhcon,
                        physcon,ithermal,&compressible,iturbulent,inomat,mi,
                        &ierr,ifreesurface,&dgravity,depth,&num_cpus));

      /* check whether an error occurred: if so, increase the shock
	 coefficient (only for compressible fluids) */
      
      if(ierr==1){
	if(iexplicit==0) FORTRAN(stop,());
	ierr=0;
	ipower=1;
	dtimefo=1.;
	ttimef=*ttime;
	timef=*time-*dtime;
	if(shockcoef==0.){shockcoef=0.001;}else{
	  shockcoef*=2;
	  if(shockcoef>2.){
	    printf("shock coefficient > 2.; stop\n");
	    FORTRAN(stop,());
	  }
	}
	memcpy(vcon,vconini,sizeof(double)*mt**nk);
	memcpy(vold,voldini,sizeof(double)*mt**nk);
	memcpy(voldo,voldini,sizeof(double)*mt**nk);
	DMEMSET(sa,0,neqp,0.);
	break;
      }
      
      /* inserting SPC's/MPC's to the physical variables */
      FORTRAN(applybounfem,(nodeboun,ndirboun,xbounact,nk,
                            vold,isolidsurf,xsolidsurf,ifreestream,iturbulent,vcon,
                            shcon,nshcon,ntmat_,physcon,v,
                            &compressible,nodempc,ipompc,coefmpc,
                            inomat,mi,ilboun,ilmpc,labmpc,
                            coefmodmpc,&iexplicit,nboun,
                            nmpc,nfreestream,nsolidsurf));

      /* calculating the conservative variables from the physical
	 variables */

      FORTRAN(phys2con,(inomat,vold,ntmat_,shcon,nshcon,
		    physcon,&compressible,vcon,rhcon,nrhcon,ithermal,
		    mi,ifreesurface,&ierr,&dgravity,depth,nk,&num_cpus));

      /* check whether an error occurred: if so, increase the shock
	 coefficient (only for compressible fluids) */
      
      if(ierr==1){
	if(iexplicit==0) FORTRAN(stop,());
	ierr=0;
	ipower=1;
	dtimefo=1.;
	ttimef=*ttime;
	timef=*time-*dtime;
	if(shockcoef==0.){shockcoef=0.001;}else{
	  shockcoef*=2;
	  if(shockcoef>2.){
	    printf("shock coefficient > 2.; stop\n");
	    FORTRAN(stop,());
	  }
	}
	memcpy(vcon,vconini,sizeof(double)*mt**nk);
	memcpy(vold,voldini,sizeof(double)*mt**nk);
	memcpy(voldo,voldini,sizeof(double)*mt**nk);
	DMEMSET(sa,0,neqp,0.);
	break;
      }
      

      /* calculating the pressure gradient for the shock smoothing in
	 the next iteration */

      if((iexplicit)&&(shockcoef>0.)){

        FORTRAN(presgradient,(iponoel,inoel,sa,&shockcoef,
                              &dtimef,ipkon,kon,lakon,vold,mi,
                              nactdoh,nk,&num_cpus));

      }

      /* calculate the sum of the square of the conservative
         variables and their change in this iteration */

      FORTRAN(calcdev,(vold,vcon,v,nk,iturbulent,mi,
                       vconmax,vmax,&iexplicit,&num_cpus));

      /* check convergence */
      
      FORTRAN(cfdconv,(nmethod,&iconvergence,ithermal,&iit,
		       iturbulent,&dtimef,vconmax,vmax));

      // why next line? Look at couettecylsegcomfem.inp and fort.12
      // problems if all nodes are in spc/mpc: v=0 for all nodes:
      // immediate convergence; cf. example couettecylsegcompfem
      if((compressible)&&(*nmethod==1)) iconvergence=0;
    
      if(((iit/jout[1])*jout[1]==iit)||(iconvergence==1)||
	 (iit==jmax[1])){

	/* calculating the stress and the heat flow at the
	   integration points, if requested */

	if((strcmp1(&filab[3306],"SF  ")==0)||
	   (strcmp1(&filab[3480],"SVF ")==0)||
	   (strcmp1(&filab[2088],"TURB")==0))isti=1;
	if(strcmp1(&filab[3393],"HFLF")==0)iqfx=1;
	for(i=0;i<*nprint;i++){
	  if(strcmp1(&prlab[6*i],"SVF")==0) isti=1;
	  if(strcmp1(&prlab[6*i],"HFLF")==0)iqfx=1;
	}
	if((isti==1)||(iqfx==1)){
	  FORTRAN(calcstressheatfluxfem,(kon,lakon,ipkon,ielmat,ntmat_,vold,
					 matname,mi,shcon,nshcon,iturbulent,
					 &compressible,ipvar,var,sti,qfx,cocon,
					 ncocon,ne,&isti,&iqfx,ithermal,rhcon,
					 nrhcon,vcon,nk));
	}
 
	/* extrapolating the stresses */

	if((strcmp1(&filab[3306],"SF  ")==0)||
	   (strcmp1(&filab[3480],"SVF ")==0)||
	   (strcmp1(&filab[2088],"TURB")==0)){
	  nfield=6;
	  ndim=6;
	  if((*norien>0)&&(strcmp1(&filab[2962],"L")==0)){
	    iorienloc=1;
	  }else{
	    iorienloc=0;
	  }
	  strcpy1(&cflag[0],&filab[2962],1);
	  NNEW(stn,double,6**nk);
	  NNEW(inum,ITG,*nk);
	  FORTRAN(extrapolatefem,(sti,stn,ipkon,inum,kon,lakon,
				  &nfield,nk,ne,mi,&ndim,orab,ielorien,co,
				  &iorienloc));
	  SFREE(inum);
	}

	/* extrapolating the heat flow */

	  
	if(strcmp1(&filab[3393],"HFLF")==0){
	  nfield=3;
	  ndim=3;
	  if((*norien>0)&&(strcmp1(&filab[3049],"L")==0)){
	    iorienloc=1;
	  }else{
	    iorienloc=0;
	  }
	  strcpy1(&cflag[0],&filab[3049],1);
	  NNEW(qfn,double,3**nk);
	  NNEW(inum,ITG,*nk);
	  FORTRAN(extrapolatefem,(qfx,qfn,ipkon,inum,kon,lakon,
				  &nfield,nk,ne,mi,&ndim,orab,ielorien,co,
				  &iorienloc));
	  SFREE(inum);
	}
	 
	/* check whether the Mach number is requested */

	if((strcmp1(&filab[1914],"MACH")==0)|| 
	   (strcmp1(&filab[3219],"TTF")==0)||
	   (strcmp1(&filab[3132],"PTF")==0)){
	  FORTRAN(calcmach,(vold,v,nk,ntmat_,shcon,nshcon,physcon,inomat,mi));
	}

	/* print output */

	if(iconvergence==1){
	  ttimf=(*ttime)+(*time);
	}else{
	  ttimf=ttimef+dtimef;
	}
      
	FORTRAN(printoutfluidfem,(set,nset,istartset,iendset,ialset,nprint,
				  prlab,prset,v,ipkon,lakon,sti,
				  mi,ithermal,co,kon,qfx,
				  &ttimf,trab,inotr,ntrans,orab,ielorien,norien,
				  vold,ielmat,thicke,
				  physcon,ielprop,prop,orname,vcon,nk,nknew,
				  nelnew));

	/* lift and drag force */

	FORTRAN(printoutfacefem,(co,ntmat_,vold,shcon,nshcon,
				 &compressible,istartset,iendset,
				 ipkon,lakon,kon,ialset,prset,&ttimf,nset,set,
				 nprint,prlab,ielmat,mi,nelnew));

	/* frd output */

	nnstep=6;
	FORTRAN(frdfluidfem,(co,nk,kon,ipkon,lakon,ne,v,vold,
			     kode,&ttimf,ielmat,matname,&nnstep,
			     physcon,filab,inomat,ntrans,inotr,trab,mi,stn,qfn,
			     istep,sa,&compressible,nactdoh,yy,jyy,ithermal,
			     shcon,nshcon,rhcon,nrhcon,ntmat_,vcon,depth,xg,
			     nodfreesurf,dt,&shockcoef,&iexplicit,nkold,
			     nelold));

	if((strcmp1(&filab[3306],"SF  ")==0)||
	   (strcmp1(&filab[3480],"SVF ")==0)||
	   (strcmp1(&filab[2088],"TURB")==0)){SFREE(stn);}
	if(strcmp1(&filab[3393],"HFLF")==0){SFREE(qfn);}

      }

      if((iit==jmax[1])||(iconvergence==1)) break;
      
      ttimef+=dtimef;

    }while(1);

    if((iit==jmax[1])||(iconvergence==1)) break;

  }while(1);
  
  if((compressible==0)&&(neqp>0)){
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
      pardiso_cleanup(&neqp,&symmetryflag,&inputformat);
#endif
    }
  }

  if(iexplicit){SFREE(sa);SFREE(voldini);SFREE(vconini);}

  if(compressible){SFREE(coefmodmpc);}

  SFREE(yy);SFREE(xsolidsurf);SFREE(dt);SFREE(dh);SFREE(vcon);
  SFREE(jyy);SFREE(vcono);

  SFREE(irowv);SFREE(irowp);
  SFREE(icolv);SFREE(icolp);
  SFREE(jqv);SFREE(jqp);
  SFREE(nactdoh);

  SFREE(adbv);SFREE(aubv);
  if((compressible==0)||(iexplicit)){SFREE(adbp);SFREE(aubp);}
  SFREE(adlv);if(iexplicit) SFREE(adlp);

  if(iexplicit){SFREE(irowvr);SFREE(jqvr);SFREE(aubvr);}

  SFREE(v);SFREE(var);SFREE(ipvar);SFREE(varf);SFREE(ipvarf);
  SFREE(voldo);

  ithermal[0]=ithermalref;
  
  if(*ifreesurface==1){SFREE(depth);SFREE(nodfreesurf);}

  *nodempcp=nodempc;*ipompcp=ipompc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
  *coefmpcp=coefmpc;*labmpcp=labmpc;*fmpcp=fmpc;*cop=co;
  *ipkonp=ipkon;*lakonp=lakon;*konp=kon;*ielmatp=ielmat;
  *nelemfacep=nelemface;*sidefacep=sideface;
  *voldp=vold;*inomatp=inomat;
  
  FORTRAN(closefilefluid,());

  *iexpl=1;
  
  return;
  
}
