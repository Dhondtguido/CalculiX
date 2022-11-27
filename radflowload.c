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

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
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
#ifdef PARDISO
#include "pardiso.h"
#endif
#ifdef PASTIX
#include "pastix.h"
#endif

static char *sideload1,*covered1=NULL;

static ITG *kontri1,*nloadtr1,*idist=NULL,*ntrit1,*mi1,*jqrad1,
  *irowrad1,*nzsrad1,num_cpus,*ntri1,*ntr1,ng1;

static double *vold1,*co1,*pmid1,*e11,*e21,*e31,*adview=NULL,*auview=NULL,*dist=NULL,
  *area1,sidemean1;

void radflowload(ITG *itg,ITG *ieg,ITG *ntg,ITG *ntr,double *adrad,
                 double *aurad,double *bcr,ITG *ipivr,
                 double *ac,double *bc,ITG *nload,char *sideload,
                 ITG *nelemload,double *xloadact,char *lakon,ITG *ipiv,
                 ITG *ntmat_,double *vold,double *shcon,
                 ITG *nshcon,ITG *ipkon,ITG *kon,double *co,
                 ITG *kontri,
                 ITG *ntri,ITG *nloadtr,double *tarea,double *tenv,
                 double *physcon,double *erad,double **adviewp, 
                 double **auviewp,
                 ITG *nflow,ITG *ikboun,
                 double *xbounact,ITG *nboun,ITG *ithermal,
                 ITG *iinc,ITG *iit,double *cs, ITG *mcs, ITG *inocs, 
                 ITG *ntrit,ITG *nk, double *fenv,ITG *istep,double *dtime,
                 double *ttime,double *time,ITG *ilboun,ITG *ikforc,
                 ITG *ilforc,double *xforcact,ITG *nforc,double *cam,
                 ITG *ielmat,ITG *nteq,double *prop,ITG *ielprop,ITG *nactdog,
                 ITG *nacteq,ITG *nodeboun,ITG *ndirboun,
                 ITG *network, double *rhcon, ITG *nrhcon, ITG *ipobody,
                 ITG *ibody, double *xbodyact, ITG *nbody,ITG *iviewfile,
                 char *jobnamef, double *ctrl, double *xloadold,
                 double *reltime, ITG *nmethod, char *set, ITG *mi,
		 ITG * istartset,ITG* iendset,ITG *ialset,ITG *nset,
                 ITG *ineighe, ITG *nmpc, ITG *nodempc,ITG *ipompc,
                 double *coefmpc,char *labmpc, ITG *iemchange,ITG *nam, 
                 ITG *iamload,ITG *jqrad,ITG *irowrad,ITG *nzsrad,
                 ITG *icolrad,ITG *ne,ITG *iaxial,double *qa,
                 double *cocon,ITG *ncocon,ITG *iponoel,ITG *inoel,
                 ITG *nprop,char *amname,ITG *namta,double *amta,ITG *iexpl){
  
  /* network=0: no network
     network=1: purely thermal (presence of "Dx"- and/or of "D " network 
                elements; declared by the user to be purely thermal 
                (on the *STEP card); simultaneous solution)
     network=2: purely thermal (alternating solution; this becomes a 
                simultaneous solution in the absence of "Dx"-elements)
     network=3: general case (temperatures, fluxes and pressures unknown)
     network=4: purely aerodynamic, i.e. only fluxes and pressures unknown

     "D "-elements (D followed by a blank) alone do not trigger the 
     alternating solution (are not counted in envtemp.f as true network 
     elements) */

  ITG nhrs=1,info=0,i,j,iin=0,icntrl,icutb=0,iin_abs=0,mt=mi[1]+1,im,
    symmetryflag=2,inputformat=1,node,ichannel,*ithread=NULL,iplausi,
    *itreated=NULL,*istack=NULL,ndata,*jumpup=NULL,*jumpdo=NULL,
    *istackb=NULL,*ixnode=NULL,*iyload=NULL;

  static ITG ifactorization=0;

  double camt[2],camf[2],camp[2],qat,qaf,ramt,ramf,ramp,
    cam1t=0.,cam1f=0.,cam1p=0.,sidemean,qa0,qau,ea,*prop_store=NULL,
    cam2t=0.,cam2f=0.,cam2p=0.,dtheta=1.,*v=NULL,cama[2],cam1a=0.,
    cam2a=0.,vamt=0.,vamf=0.,vamp=0.,vama=0.,cam0t=0.,cam0f=0.,
    cam0p=0.,cam0a=0.,sigma=0.,*adbrad=NULL,*aubrad=NULL,*q=NULL,
    *area=NULL,*pmid=NULL,*e1=NULL,*e2=NULL,*e3=NULL,
    qamt,qamf,qamtold,qamfold,*sfr=NULL,*hfr=NULL,*sba=NULL,*hba=NULL;
  
  adview=*adviewp;auview=*auviewp;

  qa0=ctrl[20];qau=ctrl[21];ea=ctrl[23];

  /* check whether there are any gas temperature nodes; this check should
     NOT be done on nteq, since also for zero equations the temperature
     of the gas nodes with boundary conditions must be stored in v
     (in initialgas) */ 

  NNEW(v,double,mt**nk);

  /* gas networks */

  if(*ntg!=0) {
#ifdef COMPANY
    NNEW(prop_store,double,*nprop);
    memcpy(&prop_store[0],&prop[0],sizeof(double)**nprop);
    FORTRAN(propertynet,(ieg,nflow,prop,ielprop,lakon,&iin,
			 prop_store,ttime,time,nam,amname,namta,amta));
    FORTRAN(checkinputvaluesnet,(ieg,nflow,prop,ielprop,lakon));
#else
    if(*iit==-1) FORTRAN(checkinputvaluesnet,(ieg,nflow,prop,ielprop,lakon));
#endif	  
    icntrl=0;
    while(icntrl==0) {
	  
      if(iin==0){
	     
	/* since in nonlingeo.c radflowload.c is called before
	   results.c at the start of a new increment vold does 
	   not contain the spc's at the end of the increment yet;
	   for this purpose voldwithspc is locally introduced
	   in radflowload.c  */

	memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);

	/* resetting ineighe to 0 for renewed call of
	   radflowload (not for cut-backs without
	   leaving radflowload */

	/* for a cut-back iin is reset to 0, iin_abs is not */

	if(iin_abs==0) DMEMSET(ineighe,0,*ntg,0);

	/* initialization pressurized flow 
	   (no free surface: gas networks or
	   water networks with fully wetted perimeter*/

	FORTRAN(initialnet,(itg,ieg,ntg,ac,bc,lakon,v,
			    ipkon,kon,nflow,
			    ikboun,nboun,prop,ielprop,nactdog,ndirboun,
			    nodeboun,xbounact,ielmat,ntmat_,shcon,nshcon,
			    physcon,ipiv,nteq,rhcon,nrhcon,ipobody,ibody,
			    xbodyact,co,nbody,network,&iin_abs,vold,set,
			    istep,iit,mi,ineighe,ilboun,&ichannel,iaxial,
			    nmpc,labmpc,ipompc,nodempc,coefmpc,ttime,time,
			    iponoel,inoel));
      
	/* initialization for channels with free surface */

	if(ichannel==1){
	  ndata=100;
	  NNEW(itreated,ITG,*nflow);
	  NNEW(istack,ITG,2**ntg);
	  NNEW(istackb,ITG,2**ntg);
	  NNEW(sfr,double,ndata**nflow);
	  NNEW(hfr,double,ndata**nflow);
	  NNEW(sba,double,ndata**nflow);
	  NNEW(hba,double,ndata**nflow);
	  NNEW(jumpup,ITG,*nflow);
	  NNEW(jumpdo,ITG,*nflow);
	  NNEW(ixnode,ITG,*nload);
	  NNEW(iyload,ITG,*nload);
	  FORTRAN(initialchannel,(itg,ieg,ntg,lakon,v,ipkon,kon,nflow,
				  ikboun,nboun,prop,ielprop,ndirboun,
				  nodeboun,xbounact,ielmat,ntmat_,shcon,nshcon,
				  physcon,rhcon,nrhcon,ipobody,ibody,
				  xbodyact,co,nbody,network,vold,set,
				  istep,iit,mi,ineighe,ilboun,ttime,time,
				  itreated,iponoel,inoel,istack,sfr,hfr,
				  sba,hba,&ndata,jumpup,jumpdo,istackb,
				  nelemload,ixnode,iyload,nload,sideload,
				  xloadact,cocon,ncocon,iinc,nforc,ikforc,
				  ilforc,xforcact));
	  SFREE(itreated);SFREE(istack);SFREE(hfr);SFREE(hba);SFREE(istackb);
	  SFREE(ixnode);SFREE(iyload);
	  break;
	}

	/* storing the residual in the rhs vector */

	FORTRAN(resultnet,(itg,ieg,ntg,bc,nload,sideload,
			   nelemload,xloadact,
			   lakon,ntmat_,v,shcon,nshcon,ipkon,kon,co,nflow,
			   iinc,istep,dtime,ttime,time,
			   ikforc,ilforc,xforcact,
			   nforc,ielmat,nteq,prop,ielprop,nactdog,nacteq,&iin,
			   physcon,camt,camf,camp,rhcon,nrhcon,ipobody,
			   ibody,xbodyact,nbody,&dtheta,vold,xloadold,
			   reltime,nmethod,set,mi,ineighe,cama,&vamt,
			   &vamf,&vamp,&vama,nmpc,nodempc,ipompc,coefmpc,
			   labmpc,iaxial,&qat,&qaf,&ramt,&ramf,&ramp,
			   cocon,ncocon,iponoel,inoel,&iplausi));

	/* iniializing qamt and qamf (mean typical energy flow
	   and mass flow */

	if(qau>1.e-10){qamt=qau;}
	else if(qa0>1.e-10){qamt=qa0;}
	else if(qat>1.e-10){qamt=qat;}
	else {qamt=1.e-2;}

	if(qau>1.e-10){qamf=qau;}
	else if(qa0>1.e-10){qamf=qa0;}
	else if(qaf>1.e-10){qamf=qaf;}
	else {qamf=1.e-2;}
      }
	  
      iin++;
      iin_abs++;
      printf("      gas iteration %" ITGFORMAT " \n \n",iin);
	  
      /* store actual values of typical energy flow and 
	 mass flow */

      qamtold=qamt;
      qamfold=qamf;

      /* filling the lhs matrix */
         
      FORTRAN(mafillnet,(itg,ieg,ntg,ac,nload,sideload,
			 nelemload,xloadact,lakon,ntmat_,v,
			 shcon,nshcon,ipkon,kon,co,nflow,iinc,
			 istep,dtime,ttime,time,
			 ielmat,nteq,prop,ielprop,nactdog,nacteq,
			 physcon,rhcon,nrhcon,ipobody,ibody,xbodyact,
			 nbody,vold,xloadold,reltime,nmethod,set,mi,
			 nmpc,nodempc,ipompc,coefmpc,labmpc,iaxial,
			 cocon,ncocon,iponoel,inoel));
	  
      /* solving the system of equations */

      if(*nteq>0){
	FORTRAN(dgesv,(nteq,&nhrs,ac,nteq,ipiv,bc,nteq,&info)); 
      }

      /*spooles(ac,au,adb,aub,&sigma,bc,icol,irow,nteq,nteq,
	&symmetryflag,&inputformat);*/
	  
      if (info!=0) {
	printf(" *WARNING in radflowload: singular matrix\n");
	    
	FORTRAN(mafillnet,(itg,ieg,ntg,ac,nload,sideload,
			   nelemload,xloadact,lakon,ntmat_,v,
			   shcon,nshcon,ipkon,kon,co,nflow,iinc,
			   istep,dtime,ttime,time,
			   ielmat,nteq,prop,ielprop,nactdog,nacteq,
			   physcon,rhcon,nrhcon,ipobody,ibody,xbodyact,
			   nbody,vold,xloadold,reltime,nmethod,set,mi,
			   nmpc,nodempc,ipompc,coefmpc,labmpc,iaxial,
			   cocon,ncocon,iponoel,inoel));
	    
	FORTRAN(equationcheck,(ac,nteq,nactdog,itg,ntg,nacteq,network));
	    
	iin=0;

      }
      else {

	/* storing the residual in the rhs vector */

	FORTRAN(resultnet,(itg,ieg,ntg,bc,nload,sideload,nelemload,
			   xloadact,lakon,ntmat_,v,shcon,nshcon,ipkon,kon,co,
			   nflow,iinc,istep,dtime,ttime,time,ikforc,ilforc,xforcact,
			   nforc,ielmat,nteq,prop,ielprop,nactdog,nacteq,
			   &iin,physcon,camt,camf,camp,rhcon,nrhcon,ipobody,
			   ibody,xbodyact,nbody,&dtheta,vold,xloadold,
			   reltime,nmethod,set,mi,ineighe,cama,&vamt,
			   &vamf,&vamp,&vama,nmpc,nodempc,ipompc,coefmpc,labmpc,
			   iaxial,&qat,&qaf,&ramt,&ramf,&ramp,cocon,ncocon,iponoel,
			   inoel,&iplausi));

	/* updating the mean typical energy flow and mass flow */

	if(qau<1.e-10){
	  if(qat>ea*qamt){qamt=(qamtold*iin+qat)/(iin+1);}
	  else {qamt=qamtold;}
	  if(qaf>ea*qamf){qamf=(qamfold*iin+qaf)/(iin+1);}
	  else {qamf=qamfold;}
	}

	/* printing the largest corrections */
	    
	if(*network!=4){ 
	  cam2t=cam1t;
	  cam1t=cam0t;
	  cam0t=camt[0];
	  printf
	    ("      mean typical energy flow since start of network iterations= %e\n",qamt);
	  printf
	    ("      largest energy flow residual in present network iteration= %e\n",ramt);
	  printf
	    ("      largest change of gas temperature since start of network iterations= %e\n",vamt);
	  if((ITG)camt[1]==0){
	    printf
	      ("      largest correction to gas temperature in present network iteration= %e\n\n",
	       fabs(camt[0]));
	  }else{
	    printf
	      ("      largest correction to gas temperature in present network iteration= %e in node %" ITGFORMAT "\n\n",
	       fabs(camt[0]),(ITG)camt[1]);
	  }
	}
	      
	if(*network>2){
	  cam2f=cam1f;
	  cam1f=cam0f;
	  cam0f=camf[0];
	  printf
	    ("      mean typical mass flow since start of network iterations= %e\n",qamf);
	  printf
	    ("      largest mass flow residual in present network iteration= %e\n",ramf);
	  printf("      largest change of gas massflow since start of network iterations= %e\n",vamf);
	  if((ITG)camf[1]==0){
	    printf("      largest correction to gas massflow in present network iteration= %e\n\n",
		   fabs(camf[0]));
	  }else{
	    printf("      largest correction to gas massflow in present network iteration= %e in node %" ITGFORMAT "\n\n",
		   fabs(camf[0]),(ITG)camf[1]);
	  }
		  
	  cam2p=cam1p;
	  cam1p=cam0p;
	  cam0p=camp[0];
	  printf
	    ("      largest element equation residual in present network iteration= %e\n",ramp);
	  printf("      largest change of gas pressure since start of network iterations= %e\n",vamp);
	  if((ITG)camp[1]==0){
	    printf("      largest correction to gas pressure in present network iteration= %e\n\n",
		   fabs(camp[0]));
	  }else{
	    printf("      largest correction to gas pressure in present network iteration= %e in node %" ITGFORMAT "\n\n",
		   fabs(camp[0]),(ITG)camp[1]);
	  }
		  
	  cam2a=cam1a;
	  cam1a=cam0a;
	  cam0a=cama[0];
	  printf("      largest change of geometry since start of network iterations= %e\n",vama);
	  if((ITG)cama[1]==0){
	    printf("      largest correction to geometry in present network iteration= %e\n",
		   fabs(cama[0]));
	  }else{
	    printf("      largest correction to geometry in present network iteration= %e in node %" ITGFORMAT "\n",
		   fabs(cama[0]),(ITG)cama[1]);
	  }
	}	      
      }
	  
      printf("\n");
	  
      /* for purely thermal calculations no iterations are
	 deemed necessary */
	  
      if(*network<=2) {icntrl=1;}
      else {

	/* check the convergence */

	checkconvnet(&icutb,&iin,
		     &cam1t,&cam1f,&cam1p,&cam2t,&cam2f,&cam2p,&cam0t,&cam0f,
		     &cam0p,&icntrl,&dtheta,ctrl,&cam1a,&cam2a,&cam0a,
		     &vamt,&vamf,&vamp,&vama,qa,&qamt,&qamf,&ramt,&ramf,&ramp,
		     &iplausi,iaxial);
      }
    }

    /* storing network output as boundary conditions for
       the structure */

    FORTRAN(flowbc,(ntg,itg,cam,vold,v,nload,sideload,
			nelemload,xloadact,nactdog,network,mi,ne,ipkon,lakon,kon));

    /* extra output in file jobname.net (only for gas networks and 
       liquid channels) */

#ifdef NETWORKOUT
    if(*network>2){
      FORTRAN(flowoutput,(itg,ieg,ntg,nteq,bc,lakon,ntmat_,
			  v,shcon,nshcon,ipkon,kon,co,nflow, dtime,ttime,time,
			  ielmat,prop,ielprop,nactdog,nacteq,&iin,physcon,
			  camt,camf,camp,rhcon,nrhcon,
			  vold,jobnamef,set,istartset,iendset,ialset,nset,
			  mi,iaxial,istep,iit,ipobody,ibody,xbodyact,nbody,
			  &ndata,sfr,sba,jumpup,jumpdo));
      if(ichannel==1){
	SFREE(sfr);SFREE(sba);SFREE(jumpup);SFREE(jumpdo);
      }
      
    }
#endif
#ifdef COMPANY
    memcpy(&prop[0],&prop_store[0],sizeof(double)**nprop);
    SFREE(prop_store);
#endif	  
  }
      
  /* radiation */

  if(*ntr>0){
      
    /* variables for multithreading procedure */
      
    ITG sys_cpus;
    char *env,*envloc,*envsys;
      
    num_cpus = 0;
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
      
    envloc = getenv("CCX_NPROC_VIEWFACTOR");
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
      
    // next line is to be inserted in a similar way for all other paralell parts
      
    if(*ntr<num_cpus) num_cpus=*ntr;
      
    pthread_t tid[num_cpus];

    /* update of vold for the boundary conditions (only in the first
       iteration of each increment and if not already done in 
       initialnet.f) */

    /* in principle ok, does lead for thermomech.inp to much more
       iterations, results the same; */
      
    /*     if(*iit<1){
	   if(*ntg==0){
	   for(i=0;i<*nboun;i++){
	   vold[mt*(nodeboun[i]-1)+ndirboun[i]]=xbounact[i];
	   }
	   }
	   }*/
      
    /*the default sink temperature is updated at the start of each
      increment */
      
    for(i=0;i<*ntr;i++){
      node=nelemload[2*nloadtr[i]-1];
      if(node!=0){
	tenv[i]=vold[mt*(node-1)]-physcon[0];
      }else if(*iit<=0){
	tenv[i]=xloadact[2*nloadtr[i]-1]-physcon[0];
      }
    }
     
    /*     for pure thermal steps the viewfactors have to be
	   calculated only once, for thermo-mechanical steps
	   (ithermal=3) they are recalculated in each iteration
	   unless they are read from file */

    if(((*ithermal==3)&&(*iviewfile>=0))||(*iit==-1)){
      if(*iviewfile<0){

	/* reading viewfactors from file */

	FORTRAN(readview,(ntr,adview,auview,fenv,nzsrad,ithermal,
			  jobnamef));

      }else{

	/* determining geometric data to calculate the viewfactors */

	NNEW(area,double,*ntrit);
	NNEW(pmid,double,3**ntrit);
	NNEW(e1,double,3**ntrit);
	NNEW(e2,double,3**ntrit);
	NNEW(e3,double,4**ntrit);

	FORTRAN(geomview,(vold,co,pmid,e1,e2,e3,kontri,area,
			  cs,mcs,inocs,ntrit,nk,mi,&sidemean));

	RENEW(adview,double,num_cpus**ntr);
	RENEW(auview,double,num_cpus*2**nzsrad);
	      
	NNEW(dist,double,num_cpus**ntrit);
	NNEW(idist,ITG,num_cpus**ntrit);

	DMEMSET(adview,0,num_cpus**ntr,0.);
	DMEMSET(auview,0,num_cpus*2**nzsrad,0.);

	sideload1=sideload;vold1=vold;co1=co;pmid1=pmid;
	e11=e1;e21=e2;e31=e3;kontri1=kontri;ntr1=ntr;
	nloadtr1=nloadtr;area1=area;ntri1=ntri;
	ntrit1=ntrit;mi1=mi;jqrad1=jqrad;irowrad1=irowrad;
	nzsrad1=nzsrad;sidemean1=sidemean;

	/* size of the square mesh used to detect
	   the visibility of a triangle; the denser
	   the mesh,the more accurate the results */

	ng1=1280;
	//	      ng1=2560;
	NNEW(covered1,char,num_cpus*ng1*ng1);

	/* calculating the viewfactors */
	      
	printf(" Using up to %" ITGFORMAT " cpu(s) for the viewfactor calculation.\n\n", num_cpus);
  
	/* create threads and wait */
	      
	NNEW(ithread,ITG,num_cpus);
	for(i=0; i<num_cpus; i++)  {
	  ithread[i]=i;
	  pthread_create(&tid[i], NULL, (void *)calcviewmt, (void *)&ithread[i]);
	}
	for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
	      
	for(i=0;i<*ntr;i++){
	  for(j=1;j<num_cpus;j++){
	    adview[i]+=adview[i+j**ntr];
	  }
	}
	RENEW(adview,double,*ntr);
	      
	for(i=0;i<2**nzsrad;i++){
	  for(j=1;j<num_cpus;j++){
	    auview[i]+=auview[i+j*2**nzsrad];
	  }
	}
	RENEW(auview,double,2**nzsrad);

	SFREE(dist);SFREE(idist);SFREE(e1);SFREE(e2);SFREE(e3);
	SFREE(pmid);SFREE(ithread);SFREE(covered1);

	/* postprocessing the viewfactors */

	FORTRAN(postview,(ntr,sideload,nelemload,kontri,ntri,nloadtr,
			  tenv,adview,auview,area,fenv,jqrad,irowrad,
			  nzsrad));

	SFREE(area);

	if(*iviewfile>=2){
		  
	  /* writing viewfactors to file */
		  
	  FORTRAN(writeview,(ntr,adview,auview,fenv,nzsrad,
			     jobnamef));
	}
	      
	if(*iviewfile==3){
		  
	  /* calculation of viewfactors only */
		  
	  FORTRAN(stop,());
	}

      }
    }

    /* assembling the radiation matrix */
      
    FORTRAN(radmatrix,(ntr,adrad,aurad,bcr,sideload,nelemload,
		       xloadact,lakon,vold,ipkon,kon,co,nloadtr,tarea,tenv,physcon,
		       erad,adview,auview,ithermal,iinc,iit,fenv,istep,dtime,ttime,
		       time,iviewfile,xloadold,reltime,nmethod,mi,iemchange,nam,
		       iamload,jqrad,irowrad,nzsrad));

    /* factoring the system of equations */

    /* the left hand side of the radiation matrix has probably
       changed if
       - the viewfactors were updated
       - a new step was started and NO CHANGE is not active
       - the emissivity coefficients were changed
       - a new increment was started in a stationary calculation
       (since the emissivity coefficients are ramped)
       in that case the LU decomposition has to be repeated
       (i.e. call of dgesv) */

    if(((*ithermal==3)&&(*iviewfile>=0))||
       ((*iit==-1)&&(*iviewfile!=-2))||
       (*iemchange==1)||
       ((*iit==0)&&(abs(*nmethod)==1))){

#if defined(PASTIX)
      // not symmetric
      ITG symmetry = 1;
      pastix_factor_main_cp(adrad,aurad,adbrad,aubrad,&sigma,icolrad,
			    irowrad,ntr,nzsrad,&symmetry,&inputformat,jqrad,
			    nzsrad);
      ifactorization=1;
#elif defined(PARDISO)
      if(ifactorization==1) pardiso_cleanup_cp(ntr,&symmetryflag,&inputformat);
      pardiso_factor_cp(adrad,aurad,adbrad,aubrad,&sigma,icolrad,
			irowrad,ntr,nzsrad,&symmetryflag,&inputformat,jqrad,
			nzsrad,iexpl);
      ifactorization=1;
#elif defined(SPOOLES)
      if(ifactorization==1) spooles_cleanup_rad();
      spooles_factor_rad(adrad,aurad,adbrad,aubrad,&sigma,
			 icolrad,irowrad,ntr,nzsrad,
			 &symmetryflag,&inputformat,iexpl);
      ifactorization=1;
#else
      printf(" *ERROR in radflowload: the SPOOLES library is not linked\n\n");
      FORTRAN(stop,());
#endif

    }

    /* solving the system of equations */

#if defined(PASTIX)
    ITG symmetry = 1;
    ITG nrhs = 1;
    pastix_solve_cp(bcr,ntr,&symmetry,&nrhs);
#elif defined(PARDISO)
    ITG nrhs = 1;
    pardiso_solve_cp(bcr,ntr,&symmetryflag,&inputformat,&nrhs);

#elif defined(SPOOLES)
    spooles_solve_rad(bcr,ntr);
#endif
	
    if (info!=0){
      printf(" *ERROR IN RADFLOWLOAD: SINGULAR MATRIX*\n");}   
      
    else{ 
      NNEW(q,double,*ntr);
      FORTRAN(radresult,(ntr,xloadact,bcr,nloadtr,tarea,
			 tenv,physcon,erad,auview,fenv,
			 irowrad,jqrad,nzsrad,q));
      SFREE(q);
    }
      
  }

  SFREE(v);

  *adviewp=adview;*auviewp=auview;

  return;

} 

/* subroutine for multithreading of calcview */

void *calcviewmt(ITG *i){

  ITG indexad,indexau,indexdi,ntria,ntrib,nedelta,indexcovered;

  indexad=*i**ntr1;
  indexau=*i*2**nzsrad1;
  indexdi=*i**ntrit1;
  indexcovered=*i*ng1*ng1;
    
  nedelta=(ITG)ceil(*ntri1/(double)num_cpus);
  ntria=*i*nedelta+1;
  ntrib=(*i+1)*nedelta;
  if(ntrib>*ntri1) ntrib=*ntri1;

  //    printf("i=%" ITGFORMAT ",ntria=%" ITGFORMAT ",ntrib=%" ITGFORMAT "\n",i,ntria,ntrib);
  //    printf("indexad=%" ITGFORMAT ",indexau=%" ITGFORMAT ",indexdi=%" ITGFORMAT "\n",indexad,indexau,indexdi);

  FORTRAN(calcview,(sideload1,vold1,co1,pmid1,e11,e21,e31,
		    kontri1,nloadtr1,&adview[indexad],
		    &auview[indexau],&dist[indexdi],&idist[indexdi],area1,
		    ntrit1,mi1,jqrad1,irowrad1,nzsrad1,&sidemean1,
		    &ntria,&ntrib,&covered1[indexcovered],&ng1));

  return NULL;
}
    
