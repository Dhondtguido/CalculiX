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

void arpack(double *co, ITG *nk, ITG **konp, ITG **ipkonp, char **lakonp,
	    ITG *ne, 
	    ITG *nodeboun, ITG *ndirboun, double *xboun, ITG *nboun, 
	    ITG *ipompc, ITG *nodempc, double *coefmpc, char *labmpc, 
	    ITG *nmpc, 
	    ITG *nodeforc, ITG *ndirforc,double *xforc, ITG *nforc, 
	    ITG *nelemload, char *sideload, double *xload,
	    ITG *nload, ITG *nactdof, 
	    ITG *icol, ITG *jq, ITG **irowp, ITG *neq, ITG *nzl, 
	    ITG *nmethod, ITG *ikmpc, ITG *ilmpc, ITG *ikboun, 
	    ITG *ilboun,
	    double *elcon, ITG *nelcon, double *rhcon, ITG *nrhcon,
	    double *shcon, ITG *nshcon, double *cocon, ITG *ncocon,
	    double *alcon, ITG *nalcon, double *alzero, ITG **ielmatp,
	    ITG **ielorienp, ITG *norien, double *orab, ITG *ntmat_,
	    double *t0, double *t1, double *t1old, 
	    ITG *ithermal,double *prestr, ITG *iprestr, 
	    double *vold,ITG *iperturb, double *sti, ITG *nzs,  
	    ITG *kode, ITG *mei, double *fei,char *filab,
	    ITG *iexpl, double *plicon, ITG *nplicon, double *plkcon,
	    ITG *nplkcon,
	    double **xstatep, ITG *npmat_, char *matname, ITG *mi,
	    ITG *ncmat_, ITG *nstate_, double **enerp, char *jobnamec,
	    char *output, char *set, ITG *nset, ITG *istartset,
	    ITG *iendset, ITG *ialset, ITG *nprint, char *prlab,
	    char *prset, ITG *nener, ITG *isolver, double *trab, 
	    ITG *inotr, ITG *ntrans, double *ttime, double *fmpc,
	    ITG *ipobody, ITG *ibody,double *xbody, ITG *nbody,
	    double *thicke, ITG *nslavs, double *tietol, ITG *nkon,
	    ITG *mpcinfo,ITG *ntie,ITG *istep,ITG *mcs,ITG *ics,
	    char *tieset,double *cs,ITG *nintpoint,ITG *mortar,
	    ITG *ifacecount,ITG **islavsurfp,double **pslavsurfp,
	    double **clearinip,ITG *nmat,char *typeboun,
	    ITG *ielprop,double *prop,char *orname,ITG *inewton,
	    double *t0g,double *t1g){

  /* calls the Arnoldi Package (ARPACK) */
  
  char bmat[2]="G", which[3]="LM", howmny[2]="A", fneig[132]="",
    description[13]="            ",*lakon=NULL,jobnamef[396]="",*labmpc2=NULL;

  ITG *inum=NULL,k,ido,ldz,iparam[11],ipntr[14],lworkl,ngraph=1,im,
    info,rvec=1,*select=NULL,lfin,j,lint,iout,ielas=1,icmd=0,mt=mi[1]+1,
    iinc=1,nev,ncv,mxiter,jrow,ifreebody,kmax,
    mass[2]={1,1}, stiffness=1, buckling=0, rhsi=0, intscheme=0,noddiam=-1,
    coriolis=0,symmetryflag=0,inputformat=0,*ipneigh=NULL,*neigh=NULL,ne0,
    *integerglob=NULL,nasym=0,zero=0,ncont=0,*itietri=NULL,kref,
    *koncont=NULL,ismallsliding=0,*itiefac=NULL,*islavsurf=NULL,
    *islavnode=NULL,*imastnode=NULL,*nslavnode=NULL,*nmastnode=NULL,
    *imastop=NULL,*iponoels=NULL,*inoels=NULL,*ipe=NULL,*ime=NULL,
    mpcfree,memmpc_,icascade,maxlenmpc,iit=-1,*irow=NULL,nherm=1,
    icfd=0,*inomat=NULL,*ipkon=NULL,*kon=NULL,*ielmat=NULL,*ielorien=NULL,
    *islavact=NULL,maxprevcontel,iex,nslavs_prev_step,icutb=0,
    iflagact=0,*islavsurfold=NULL,ialeatoric=0,kscale=1,network=0,
    *iponoel=NULL,*inoel=NULL,nrhs=1,igreen=0,node,idir,mscalmethod=0,
    *jqw=NULL,*iroww=NULL,nzsw,
    *islavelinv=NULL,*irowtloc=NULL,*jqtloc=NULL,nboun2,
    *ndirboun2=NULL,*nodeboun2=NULL,nmpc2,*ipompc2=NULL,*nodempc2=NULL,
    *ikboun2=NULL,*ilboun2=NULL,*ikmpc2=NULL,*ilmpc2=NULL,mortartrafoflag=0;

  double *stn=NULL,*v=NULL,*resid=NULL,*z=NULL,*workd=NULL,
    *workl=NULL,*d=NULL,sigma=1,*temp_array=NULL,*pslavsurfold=NULL,
    *een=NULL,sum,cam[5],*f=NULL,*fn=NULL,qa[4],*fext=NULL,
    *epn=NULL,*fnr=NULL,*fni=NULL,*emn=NULL,*emeini=NULL,
    *xstateini=NULL,*xstiff=NULL,*stiini=NULL,*vini=NULL,freq,*stx=NULL,
    *enern=NULL,*xstaten=NULL,*eei=NULL,*enerini=NULL,*fnext=NULL,
    *physcon=NULL,*qfx=NULL,*qfn=NULL,tol,fmin,fmax,pi,*cgr=NULL,
    *xloadold=NULL,reltime=1.,*vr=NULL,*vi=NULL,*stnr=NULL,*stni=NULL,
    *vmax=NULL,*stnmax=NULL,*springarea=NULL,*eenmax=NULL,
    *doubleglob=NULL,*cg=NULL,*straight=NULL,*clearini=NULL,
    *xmastnor=NULL,*areaslav=NULL,*xnoels=NULL,theta=0.,
    *di=NULL,sigmai=0,*workev=NULL,*ener=NULL,*xstate=NULL,*dc=NULL,
    *au=NULL,*ad=NULL,*b=NULL,*aub=NULL,*adb=NULL,*pslavsurf=NULL,
    *pmastsurf=NULL,*cdn=NULL,*energyini=NULL,*energy=NULL,
    *cdnr=NULL,*cdni=NULL,*eme=NULL,alea=0.1,*smscale=NULL,zmax,*auw=NULL,
    *autloc=NULL,*xboun2=NULL,*coefmpc2=NULL;

  FILE *f1;

  /* dummy arguments for the results call */

  double *veold=NULL,*accold=NULL,bet,gam,dtime=1.,time=1.;

#ifdef SGI
  ITG token;
#endif

  irow=*irowp;ener=*enerp;xstate=*xstatep;ipkon=*ipkonp;lakon=*lakonp;
  kon=*konp;ielmat=*ielmatp;ielorien=*ielorienp;

  islavsurf=*islavsurfp;pslavsurf=*pslavsurfp;clearini=*clearinip;

  if(*nmethod==13){
    *nmethod=2;
    igreen=1;
  }
  
  for(k=0;k<3;k++){
    strcpy1(&jobnamef[k*132],&jobnamec[k*132],132);
  }
  
  if(*mortar!=1){
    maxprevcontel=*nslavs;
  }else if(*mortar==1){
    maxprevcontel=*nintpoint;
    if(*nstate_!=0){
      if(maxprevcontel!=0){
	NNEW(islavsurfold,ITG,2**ifacecount+2);
	NNEW(pslavsurfold,double,3**nintpoint);
	memcpy(&islavsurfold[0],&islavsurf[0],
	       sizeof(ITG)*(2**ifacecount+2));
	memcpy(&pslavsurfold[0],&pslavsurf[0],
	       sizeof(double)*(3**nintpoint));
      }
    }
  }
  nslavs_prev_step=*nslavs;

  if((strcmp1(&filab[870],"PU  ")==0)||
     (strcmp1(&filab[1479],"PHS ")==0)||
     (strcmp1(&filab[1566],"MAXU")==0)||
     (strcmp1(&filab[1653],"MAXS")==0)){
    printf(" *WARNING in arpack: PU, PHS, MAXU or MAX was selected in a frequency calculation without cyclic symmetry;\n this is not correct; output request is removed;\n");
    strcpy1(&filab[870],"    ",4);
    strcpy1(&filab[1479],"    ",4);
    strcpy1(&filab[1566],"    ",4);
    strcpy1(&filab[1653],"    ",4);
  }

  /* copying the frequency parameters */

  pi=4.*atan(1.);

  if(igreen==0){
    nev=mei[0];
    ncv=mei[1];
    mxiter=mei[2];
    tol=fei[0];
    fmin=2*pi*fei[1];
    fmax=2*pi*fei[2];
  }else{

    /* Green functions: number of "eigenmodes" is the number of
       unit forces */
    
    nev=*nforc;
    ncv=*nforc;
  }
  
  /* assigning the body forces to the elements */ 

  if(*nbody>0){
    if(*inewton==1){
      printf(" *ERROR in arpack: generalized gravity loading is not allowed in frequency calculations");
      FORTRAN(stop,());
    }
  }

  ne0=*ne;

  /* contact conditions */
  
  if(*iperturb!=0){

    memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];icascade=mpcinfo[2];
    maxlenmpc=mpcinfo[3];

    inicont(nk,&ncont,ntie,tieset,nset,set,istartset,iendset,ialset,&itietri,
	    lakon,ipkon,kon,&koncont,nslavs,tietol,&ismallsliding,&itiefac,
	    &islavsurf,&islavnode,&imastnode,&nslavnode,&nmastnode,
	    mortar,&imastop,nkon,&iponoels,&inoels,&ipe,&ime,ne,ifacecount,
	    iperturb,ikboun,nboun,co,istep,&xnoels);

    if(ncont!=0){

      NNEW(cg,double,3*ncont);
      NNEW(straight,double,16*ncont);
	  
      /* 11 instead of 10: last position is reserved for the
	 local contact spring element number; needed as
	 pointer into springarea */
	  
      if(*mortar==0){
	RENEW(kon,ITG,*nkon+11**nslavs);
	NNEW(springarea,double,2**nslavs);
	if(*nener==1){
	  RENEW(ener,double,2*mi[0]*(*ne+*nslavs));
	  DMEMSET(ener,0,2*mi[0]*(*ne+*nslavs),0.);
	}
	RENEW(ipkon,ITG,*ne+*nslavs);
	RENEW(lakon,char,8*(*ne+*nslavs));
	      
	if(*norien>0){
	  RENEW(ielorien,ITG,mi[2]*(*ne+*nslavs));
	  for(k=mi[2]**ne;k<mi[2]*(*ne+*nslavs);k++) ielorien[k]=0;
	}
	      
	RENEW(ielmat,ITG,mi[2]*(*ne+*nslavs));
	for(k=mi[2]**ne;k<mi[2]*(*ne+*nslavs);k++) ielmat[k]=1;
	      
	if((maxprevcontel==0)&&(*nslavs!=0)){
	  RENEW(xstate,double,*nstate_*mi[0]*(*ne+*nslavs));
	  for(k=*nstate_*mi[0]**ne;k<*nstate_*mi[0]*(*ne+*nslavs);k++){
	    xstate[k]=0.;
	  }
	}
	maxprevcontel=*nslavs;
	      
	NNEW(areaslav,double,*ifacecount);
	NNEW(xmastnor,double,3*nmastnode[*ntie]);
      }else if(*mortar==1){
	NNEW(islavact,ITG,nslavnode[*ntie]);
	DMEMSET(islavact,0,nslavnode[*ntie],1);
	if((*istep==1)||(nslavs_prev_step==0)) NNEW(clearini,double,3*9**ifacecount);
	NNEW(xmastnor,double,3*nmastnode[*ntie]);

	*nintpoint=0;
	      
	precontact(&ncont,ntie,tieset,nset,set,istartset,
		   iendset,ialset,itietri,lakon,ipkon,kon,koncont,ne,
		   cg,straight,co,vold,istep,&iinc,&iit,itiefac,
		   islavsurf,islavnode,imastnode,nslavnode,nmastnode,
		   imastop,mi,ipe,ime,tietol,&iflagact,
		   nintpoint,&pslavsurf,xmastnor,cs,mcs,ics,clearini,
		   nslavs);
	      
	/* changing the dimension of element-related fields */
	      
	RENEW(kon,ITG,*nkon+22**nintpoint);
	RENEW(springarea,double,2**nintpoint);
	RENEW(pmastsurf,double,6**nintpoint);
	      
	if(*nener==1){
	  RENEW(ener,double,2*mi[0]*(*ne+*nintpoint));
	  DMEMSET(ener,0,2*mi[0]*(*ne+*nslavs),0.);
	}
	RENEW(ipkon,ITG,*ne+*nintpoint);
	RENEW(lakon,char,8*(*ne+*nintpoint));
	      
	if(*norien>0){
	  RENEW(ielorien,ITG,mi[2]*(*ne+*nintpoint));
	  for(k=mi[2]**ne;k<mi[2]*(*ne+*nintpoint);k++) ielorien[k]=0;
	}
	RENEW(ielmat,ITG,mi[2]*(*ne+*nintpoint));
	for(k=mi[2]**ne;k<mi[2]*(*ne+*nintpoint);k++) ielmat[k]=1;

	/* interpolating the state variables */

	if(*nstate_!=0){
	  if(maxprevcontel!=0){
	    RENEW(xstateini,double,
		  *nstate_*mi[0]*(ne0+maxprevcontel));
	    for(k=*nstate_*mi[0]*ne0;
		k<*nstate_*mi[0]*(ne0+maxprevcontel);++k){
	      xstateini[k]=xstate[k];
	    }
	  }
		  
	  RENEW(xstate,double,*nstate_*mi[0]*(ne0+*nintpoint));
	  for(k=*nstate_*mi[0]*ne0;k<*nstate_*mi[0]*(ne0+*nintpoint);k++){
	    xstate[k]=0.;
	  }
		  
	  if((*nintpoint>0)&&(maxprevcontel>0)){
	    iex=2;
		      
	    /* interpolation of xstate */
		      
	    FORTRAN(interpolatestate,(ne,ipkon,kon,lakon,
				      &ne0,mi,xstate,pslavsurf,nstate_,
				      xstateini,islavsurf,islavsurfold,
				      pslavsurfold,tieset,ntie,itiefac));
		      
	  }

	  if(maxprevcontel!=0){
	    SFREE(islavsurfold);SFREE(pslavsurfold);
	  }

	  maxprevcontel=*nintpoint;
		  
	  RENEW(xstateini,double,*nstate_*mi[0]*(ne0+*nintpoint));
	  for(k=0;k<*nstate_*mi[0]*(ne0+*nintpoint);++k){
	    xstateini[k]=xstate[k];
	  }
	}
      }
	  
      /* generating contact spring elements */

      contact(&ncont,ntie,tieset,nset,set,istartset,iendset,
	      ialset,itietri,lakon,ipkon,kon,koncont,ne,cg,straight,nkon,
	      co,vold,ielmat,cs,elcon,istep,&iinc,&iit,ncmat_,ntmat_,
	      &ne0,nmethod,
	      iperturb,ikboun,nboun,mi,imastop,nslavnode,islavnode,islavsurf,
	      itiefac,areaslav,iponoels,inoels,springarea,tietol,&reltime,
	      imastnode,nmastnode,xmastnor,filab,mcs,ics,&nasym,
	      xnoels,mortar,pslavsurf,pmastsurf,clearini,&theta,
	      xstateini,xstate,nstate_,&icutb,&ialeatoric,jobnamef,&alea,
	      auw,jqw,iroww,&nzsw);
	  
      printf("number of contact spring elements=%" ITGFORMAT "\n\n",*ne-ne0);
	  
      /* determining the structure of the stiffness/mass matrix */
	  
      remastructar(ipompc,&coefmpc,&nodempc,nmpc,
		   &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,ikboun,ilboun,
		   labmpc,nk,&memmpc_,&icascade,&maxlenmpc,
		   kon,ipkon,lakon,ne,nactdof,icol,jq,&irow,isolver,
		   neq,nzs,nmethod,ithermal,iperturb,mass,mi,ics,cs,
		   mcs,mortar,typeboun,&iit,&network,iexpl);
    }
  }

  /* field for initial values of state variables (needed if
     previous static step was viscoplastic and for contact */

  if((*nstate_!=0)&&((*mortar==0)||(ncont==0))){
    NNEW(xstateini,double,*nstate_*mi[0]*(ne0+*nslavs));
    for(k=0;k<*nstate_*mi[0]*(ne0+*nslavs);++k){
      xstateini[k]=xstate[k];
    }
  }

  /* determining the internal forces and the stiffness coefficients */

  NNEW(f,double,neq[1]);

  /* allocating a field for the stiffness matrix */

  NNEW(xstiff,double,(long long)27*mi[0]**ne);

  iout=-1;
  NNEW(v,double,mt**nk);
  memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
  NNEW(fn,double,mt**nk);
  NNEW(stx,double,6*mi[0]**ne);
  NNEW(eme,double,6*mi[0]**ne);
  if(*ithermal>1){
    NNEW(qfx,double,3*mi[0]**ne);
  }
  NNEW(inum,ITG,*nk);

  /* following fields are needed if *ener==1 or kode<-100 */
    
  NNEW(vini,double,mt**nk);
  NNEW(eei,double,6*mi[0]**ne);
  NNEW(stiini,double,6*mi[0]*ne0);
  NNEW(emeini,double,6*mi[0]*ne0);
  if(*nener==1) NNEW(enerini,double,2*mi[0]*ne0);
    
  if(*iperturb==0){
    results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	    elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	    ielorien,norien,orab,ntmat_,t0,t0,ithermal,
	    prestr,iprestr,filab,eme,emn,een,iperturb,
	    f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	    ndirboun,xboun,nboun,ipompc,
	    nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
	    &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	    xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	    &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
	    emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
	    iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
	    fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
	    &reltime,&ne0,thicke,shcon,nshcon,
	    sideload,xload,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
	    mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
	    islavsurf,ielprop,prop,energyini,energy,&kscale,iponoel,
	    inoel,nener,orname,&network,ipobody,xbody,ibody,typeboun,
	    itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
	    islavelinv,autloc,irowtloc,jqtloc,&nboun2,
	    ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
	    labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
	    &intscheme,physcon);
  }else{
    results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	    elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	    ielorien,norien,orab,ntmat_,t1old,t1old,ithermal,
	    prestr,iprestr,filab,eme,emn,een,iperturb,
	    f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	    ndirboun,xboun,nboun,ipompc,
	    nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
	    &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	    xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	    &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
	    emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
	    iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
	    fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
	    &reltime,&ne0,thicke,shcon,nshcon,
	    sideload,xload,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
	    mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
	    islavsurf,ielprop,prop,energyini,energy,&kscale,iponoel,
	    inoel,nener,orname,&network,ipobody,xbody,ibody,typeboun,
	    itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
	    islavelinv,autloc,irowtloc,jqtloc,&nboun2,
	    ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
	    labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
	    &intscheme,physcon);
  }
  
  SFREE(eei);SFREE(stiini);SFREE(emeini);SFREE(vini);
  if(*nener==1) SFREE(enerini);
    
  SFREE(f);SFREE(v);SFREE(fn);SFREE(stx);SFREE(eme);
  if(*ithermal>1)SFREE(qfx);SFREE(inum);
  iout=1;

  /* for the frequency analysis linear strain and elastic properties
     are used */

  iperturb[1]=0;

  /* filling in the matrix */

  NNEW(ad,double,neq[1]);
  NNEW(au,double,nzs[2]);

  NNEW(adb,double,neq[1]);
  NNEW(aub,double,nzs[1]);

  NNEW(fext,double,neq[1]);

  if(*iperturb==0){
    mafillsmmain(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,nboun,
		 ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
		 nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
		 ad,au,fext,nactdof,icol,jq,irow,neq,nzl,nmethod,
		 ikmpc,ilmpc,ikboun,ilboun,
		 elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		 ielorien,norien,orab,ntmat_,
		 t0,t0,ithermal,prestr,iprestr,vold,iperturb,sti,
		 nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
		 xstiff,npmat_,&dtime,matname,mi,
		 ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
		 physcon,shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,
		 &coriolis,ibody,xloadold,&reltime,veold,springarea,nstate_,
		 xstateini,xstate,thicke,integerglob,doubleglob,
		 tieset,istartset,iendset,ialset,ntie,&nasym,pslavsurf,
		 pmastsurf,mortar,clearini,ielprop,prop,&ne0,fnext,&kscale,
		 iponoel,inoel,&network,ntrans,inotr,trab,smscale,&mscalmethod,
		 set,nset,islavelinv,autloc,irowtloc,jqtloc,&mortartrafoflag);
  }
  else{
    mafillsmmain(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,nboun,
		 ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
		 nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
		 ad,au,fext,nactdof,icol,jq,irow,neq,nzl,nmethod,
		 ikmpc,ilmpc,ikboun,ilboun,
		 elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		 ielorien,norien,orab,ntmat_,
		 t0,t1old,ithermal,prestr,iprestr,vold,iperturb,sti,
		 nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
		 xstiff,npmat_,&dtime,matname,mi,
		 ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
		 physcon,shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,
		 &coriolis,ibody,xloadold,&reltime,veold,springarea,nstate_,
		 xstateini,xstate,thicke,integerglob,doubleglob,
		 tieset,istartset,iendset,ialset,ntie,&nasym,pslavsurf,
		 pmastsurf,mortar,clearini,ielprop,prop,&ne0,fnext,&kscale,
		 iponoel,inoel,&network,ntrans,inotr,trab,smscale,&mscalmethod,
		 set,nset,islavelinv,autloc,irowtloc,jqtloc,&mortartrafoflag);

    if(nasym==1){
      RENEW(au,double,nzs[2]+nzs[1]);
      RENEW(aub,double,nzs[2]+nzs[1]);
      symmetryflag=2;
      inputformat=1;
	  
      mafillsmasmain(co,nk,kon,ipkon,lakon,ne,nodeboun,
		     ndirboun,xboun,nboun,
		     ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
		     nforc,nelemload,sideload,xload,nload,xbody,ipobody,
		     nbody,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,
		     nmethod,ikmpc,ilmpc,ikboun,ilboun,
		     elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		     ielmat,ielorien,norien,orab,ntmat_,
		     t0,t1old,ithermal,prestr,iprestr,vold,iperturb,sti,
		     nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
		     xstiff,npmat_,&dtime,matname,mi,
		     ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
		     physcon,shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,
		     &coriolis,ibody,xloadold,&reltime,veold,springarea,nstate_,
		     xstateini,xstate,thicke,
		     integerglob,doubleglob,tieset,istartset,iendset,
		     ialset,ntie,&nasym,pslavsurf,pmastsurf,mortar,clearini,
		     ielprop,prop,&ne0,&kscale,iponoel,inoel,&network,set,nset);
    }
  }

  SFREE(fext);

  if(*nmethod==0){

    /* error occurred in mafill: storing the geometry in frd format */

    ++*kode;
    NNEW(inum,ITG,*nk);for(k=0;k<*nk;k++) inum[k]=1;
    if(strcmp1(&filab[1044],"ZZS")==0){
      NNEW(neigh,ITG,40**ne);
      NNEW(ipneigh,ITG,*nk);
    }

    frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	kode,filab,een,t1,fn,&time,epn,ielmat,matname,enern,xstaten,
	nstate_,istep,&iinc,ithermal,qfn,&j,&noddiam,trab,inotr,
	ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni,nmat,
	ielprop,prop,sti);
    
    if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);}
    SFREE(inum);FORTRAN(stop,());

  }

  /* LU decomposition of the left hand matrix */

  if(igreen==0){
    if(nasym==1){sigma=0.;}else{sigma=1.;}
  }else{
    if(*nforc>0){sigma=xforc[0];}
  }

  if(*isolver==0){
#ifdef SPOOLES
    spooles_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
                   &symmetryflag,&inputformat,&nzs[2]);
#else
    printf(" *ERROR in arpack: the SPOOLES library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==4){
#ifdef SGI
    token=1;
    sgi_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],token);
#else
    printf(" *ERROR in arpack: the SGI library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==5){
#ifdef TAUCS
    tau_factor(ad,&au,adb,aub,&sigma,icol,&irow,&neq[1],&nzs[1]);
#else
    printf(" *ERROR in arpack: the TAUCS library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==6){
#ifdef MATRIXSTORAGE
    matrixstorage(ad,&au,adb,aub,&sigma,icol,&irow,&neq[1],&nzs[1],
		  ntrans,inotr,trab,co,nk,nactdof,jobnamec,mi,ipkon,
                  lakon,kon,ne,mei,nboun,nmpc,cs,mcs,ithermal,nmethod);

    strcpy2(fneig,jobnamec,132);
    strcat(fneig,".frd");
    if((f1=fopen(fneig,"ab"))==NULL){
      printf(" *ERROR in frd: cannot open frd file for writing...");
      exit(0);
    }
    fprintf(f1," 9999\n");
    fclose(f1);
    FORTRAN(stopwithout201,());
#else
    printf(" *ERROR in arpack: the MATRIXSTORAGE library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==7){
#ifdef PARDISO
    pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
		   &symmetryflag,&inputformat,jq,&nzs[2]);
#else
    printf(" *ERROR in arpack: the PARDISO library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==8){
#ifdef PASTIX
#ifdef PARDISO
    pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
		   &symmetryflag,&inputformat,jq,&nzs[2]);
#else
    pastix_factor_main(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
		       &symmetryflag,&inputformat,jq,&nzs[2]);
#endif
#else
    printf(" *ERROR in arpack: the PASTIX library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }

  if(igreen==0){
  
    /* calculating the eigenvalues and eigenmodes */
  
    printf(" Calculating the eigenvalues and the eigenmodes\n\n");

    ido=0;
    ldz=neq[1];
    iparam[0]=1;
    iparam[2]=mxiter;
    iparam[3]=1;
    iparam[6]=3;
    
    info=0;

    NNEW(resid,double,neq[1]);
    NNEW(z,double,(long long)ncv*neq[1]);
    NNEW(workd,double,3*neq[1]);

    /* for the limits: look for "ierr = -3" in dsaupd.f and dnaupd.f */
    
    if(nasym==1){

      /* nev+2 <= ncv <= n */
      
      if(nev+2>neq[1]){
	printf(" *ERROR in arpack: too many eigenvalues requested\n");
	printf("                   number is reduced to: %d\n\n",neq[1]-2);
	nev=neq[1]-2;
      }
      if(ncv<nev+2){
	ncv=nev+2;
      }else{
	if(ncv>neq[1]){ncv=neq[1];}
      }
      lworkl=3*ncv*(2+ncv);
      NNEW(workl,double,lworkl);
      FORTRAN(dnaupd,(&ido,bmat,&neq[1],which,&nev,&tol,resid,&ncv,z,&ldz,iparam,ipntr,workd,
		      workl,&lworkl,&info));
    }else{

      /* nev+1 <= ncv <= n */
      
      if(nev+1>neq[1]){
	printf(" *ERROR in arpack: too many eigenvalues requested\n");
	printf("                   number is reduced to: %d\n\n",neq[1]-1);
	nev=neq[1]-1;
      }
      if(ncv<nev+1){
	ncv=nev+1;
      }else{
	if(ncv>neq[1]){ncv=neq[1];}
      }
      lworkl=ncv*(8+ncv);
      NNEW(workl,double,lworkl);
      FORTRAN(dsaupd,(&ido,bmat,&neq[1],which,&nev,&tol,resid,&ncv,z,&ldz,iparam,ipntr,workd,
		      workl,&lworkl,&info));
    }
    
    NNEW(temp_array,double,neq[1]);
    
    while((ido==-1)||(ido==1)||(ido==2)){
      if(ido==-1){
	if(nasym==1){
	  FORTRAN(opas,(&neq[1],&workd[ipntr[0]-1],temp_array,adb,aub,jq,irow,nzs));
	}else{
	  FORTRAN(op,(&neq[1],&workd[ipntr[0]-1],temp_array,adb,aub,jq,irow));
	}
      }
      if((ido==-1)||(ido==1)){
	
	/* solve the linear equation system  */
	
	if(ido==-1){
	  if(*isolver==0){
#ifdef SPOOLES
	    spooles_solve(temp_array,&neq[1]);
#endif
	  }
	  else if(*isolver==4){
#ifdef SGI
	    sgi_solve(temp_array,token);
#endif
	  }
	  else if(*isolver==5){
#ifdef TAUCS
	    tau_solve(temp_array,&neq[1]);
#endif
	  }
	  else if(*isolver==7){
#ifdef PARDISO
	    pardiso_solve(temp_array,&neq[1],&symmetryflag,&inputformat,&nrhs);
#endif
	  }
	  else if(*isolver==8){
#ifdef PASTIX
#ifdef PARDISO
	    pardiso_solve(temp_array,&neq[1],&symmetryflag,&inputformat,&nrhs);
#else	    
	    if( pastix_solve(temp_array,&neq[1],&symmetryflag,&nrhs)==-1 )
	      printf(" *WARNING in arpack: solving step didn't converge! Continuing anyway\n");
#endif	    
#endif
	  }
	  for(jrow=0;jrow<neq[1];jrow++){
	    workd[ipntr[1]-1+jrow]=temp_array[jrow];
	  }
	}
	else if(ido==1){
	  if(*isolver==0){
#ifdef SPOOLES
	    spooles_solve(&workd[ipntr[2]-1],&neq[1]);
#endif
	  }
	  else if(*isolver==4){
#ifdef SGI
	    sgi_solve(&workd[ipntr[2]-1],token);
#endif
	  }
	  else if(*isolver==5){
#ifdef TAUCS
	    tau_solve(&workd[ipntr[2]-1],&neq[1]);
#endif
	  }
	  else if(*isolver==7){
#ifdef PARDISO
	    pardiso_solve(&workd[ipntr[2]-1],&neq[1],&symmetryflag,&inputformat,&nrhs);
#endif
	  }
	  else if(*isolver==8){
#ifdef PASTIX
#ifdef PARDISO
	    pardiso_solve(&workd[ipntr[2]-1],&neq[1],&symmetryflag,&inputformat,&nrhs);
#else
	    if( pastix_solve(&workd[ipntr[2]-1],&neq[1],&symmetryflag,&nrhs)==-1 )
	      printf(" *WARNING in arpack: solving step didn't converge! Continuing anyway\n");
#endif
#endif
	  }
	  for(jrow=0;jrow<neq[1];jrow++){
	    workd[ipntr[1]-1+jrow]=workd[ipntr[2]-1+jrow];
	  }
	}
	
      }
      
      if(ido==2){
	if(nasym==1){
	  FORTRAN(opas,(&neq[1],&workd[ipntr[0]-1],&workd[ipntr[1]-1],
			adb,aub,jq,irow,nzs));
	}else{
	  FORTRAN(op,(&neq[1],&workd[ipntr[0]-1],&workd[ipntr[1]-1],
		      adb,aub,jq,irow));
	}
      }
      
      if(nasym==1){
	FORTRAN(dnaupd,(&ido,bmat,&neq[1],which,&nev,&tol,resid,&ncv,z,&ldz,
                        iparam,ipntr,workd,workl,&lworkl,&info));
      }else{
	FORTRAN(dsaupd,(&ido,bmat,&neq[1],which,&nev,&tol,resid,&ncv,z,&ldz,
			iparam,ipntr,workd,workl,&lworkl,&info));
      }
    }
    
    if(info!=0){
      printf(" *ERROR in d[n,s]aupd: info=%" ITGFORMAT "\n",info);
      printf("       # of converged eigenvalues=%" ITGFORMAT "\n\n",iparam[4]);
    }         
    
    NNEW(select,ITG,ncv);
    
    if(nasym==1){
      NNEW(d,double,nev+1);
      NNEW(di,double,nev+1);
      NNEW(workev,double,3*ncv);
      FORTRAN(dneupd,(&rvec,howmny,select,d,di,z,&ldz,&sigma,&sigmai,
		      workev,bmat,&neq[1],which,&nev,&tol,resid,
		      &ncv,z,&ldz,iparam,ipntr,workd,workl,&lworkl,&info));
      SFREE(workev);
      NNEW(dc,double,2*nev);
      
      /* storing as complex number and taking the square root */
      
      for(j=0;j<nev;j++){
	dc[2*j]=sqrt(sqrt(d[j]*d[j]+di[j]*di[j])+d[j])/sqrt(2.);
	dc[2*j+1]=sqrt(sqrt(d[j]*d[j]+di[j]*di[j])-d[j])/sqrt(2.);
	if(di[j]<0.) dc[2*j+1]=-dc[2*j+1];
      }
      FORTRAN(writeevcomplex,(dc,&nev,&fmin,&fmax));
      SFREE(di);SFREE(dc);
    }else{
      NNEW(d,double,nev);
      FORTRAN(dseupd,(&rvec,howmny,select,d,z,&ldz,&sigma,bmat,&neq[1],which,
		      &nev,&tol,resid,&ncv,z,&ldz,iparam,ipntr,workd,workl,&lworkl,&info));
      FORTRAN(writeev,(d,&nev,&fmin,&fmax));
    }
    SFREE(select);SFREE(workd);SFREE(workl);SFREE(resid);
	
    if(info!=0){
      printf(" *ERROR in d[n,s]eupd: info=%" ITGFORMAT "\n",info);
    }         

    /* check the normalization of the eigenmodes */
    
    for(j=0;j<nev;j++){
      zmax=0.;
      kmax=-1;
      kref=j*neq[1];
      if(nasym==1){
	FORTRAN(opas,(&neq[1],&z[kref],temp_array,
		      adb,aub,jq,irow,nzs));
      }else{
	FORTRAN(op,(neq,&z[kref],temp_array,
		    adb,aub,jq,irow));
      }
      sum=0;
      for(k=0;k<neq[1];k++){sum+=z[kref+k]*temp_array[k];}
      printf("U^T*M*U=%f for eigenmode %d\n",j+1,sum);

      /* check for the maximum */

      /*      for(k=0;k<neq[1];k++){
	if(fabs(z[kref+k])>zmax){
	  zmax=fabs(z[kref+k]);
	  kmax=k;
	}
      }
      zmax=z[kref+k];*/
      
      /* normalizing the eigenmode (if not yet normalized by ARPACK */
      
      /*      if((fabs(sum-1.)>1.e-5)&&(sum>0.)){
	printf("sum=%e\n");
	printf("normalizing the eigenmode \n");
	if(zmax>0.){
	  for(k=0;k<neq[1];k++){z[kref+k]/=sqrt(sum);}
	}else{
	  for(k=0;k<neq[1];k++){z[kref+k]/=-sqrt(sum);}
	}
      }else if(zmax<0.){
	for(k=0;k<neq[1];k++){z[kref+k]*=-1.;}
	}*/
      
      if(fabs(sum-1.)>1.e-5){
	printf("sum=%e\n",sum);
	printf("normalizing the eigenmode \n");
	for(k=0;k<neq[1];k++){z[kref+k]/sqrt(sum);}
      }
    }
    
    SFREE(temp_array);

  }else{

    /* Green functions */

    /* storing omega_0^2 into d */
    
    NNEW(d,double,*nforc);
    for(j=0;j<*nforc;j++){d[j]=xforc[0];}
    NNEW(z,double,neq[1]);

  }
    
  
  /* writing the eigenvalues and mass matrix to a binary file */
  
  if(mei[3]==1){
    
    strcpy2(fneig,jobnamec,132);
    strcat(fneig,".eig");
      
    if((f1=fopen(fneig,"wb"))==NULL){
      printf(" *ERROR in arpack: cannot open eigenvalue file for writing...");

      exit(0);
    }

    /* storing a zero as indication that this was not a
       cyclic symmetry calculation */

    if(fwrite(&zero,sizeof(ITG),1,f1)!=1){
      printf(" *ERROR saving the cyclic symmetry flag to the eigenvalue file...");
      exit(0);
    }

    /* Hermitian */

    if(fwrite(&nherm,sizeof(ITG),1,f1)!=1){
      printf(" *ERROR saving the Hermitian flag to the eigenvalue file...");
      exit(0);
    }

    /* perturbation parameter iperturb[0] */

    if(fwrite(&iperturb[0],sizeof(ITG),1,f1)!=1){
      printf(" *ERROR saving the perturbation flag to the eigenvalue file...");
      exit(0);
    }

    /* reference displacements */

    if(iperturb[0]==1){
      if(fwrite(vold,sizeof(double),mt**nk,f1)!=mt**nk){
	printf(" *ERROR saving the reference displacements to the eigenvalue file...");
	exit(0);
      }
    }

    /* storing the number of eigenvalues */

    if(fwrite(&nev,sizeof(ITG),1,f1)!=1){
      printf(" *ERROR saving the number of eigenvalues to the eigenvalue file...");
      exit(0);
    }

    /* the eigenfrequencies are stored as radians/time */

    if(fwrite(d,sizeof(double),nev,f1)!=nev){
      printf(" *ERROR saving the eigenfrequencies to the eigenvalue file...");
      exit(0);
    }

    /* storing the stiffness matrix */

    if(fwrite(ad,sizeof(double),neq[1],f1)!=neq[1]){
      printf(" *ERROR saving the diagonal of the stiffness matrix to the eigenvalue file...");
      exit(0);
    }
    if(fwrite(au,sizeof(double),nzs[2],f1)!=nzs[2]){
      printf(" *ERROR saving the off-diagonal terms of the stiffness matrix to the eigenvalue file...");
      exit(0);
    }

    /* storing the mass matrix */

    if(fwrite(adb,sizeof(double),neq[1],f1)!=neq[1]){
      printf(" *ERROR saving the diagonal of the mass matrix to the eigenvalue file...");
      exit(0);
    }
    if(fwrite(aub,sizeof(double),nzs[1],f1)!=nzs[1]){
      printf(" *ERROR saving the off-diagonal terms of the mass matrix to the eigenvalue file...");
      exit(0);
    }
  }

  if(igreen==0){
    
    /* calculating the participation factors and the relative effective
       modal mass */
    
    if(*ithermal!=2){
      FORTRAN(effectivemodalmass,(neq,nactdof,mi,adb,aub,jq,irow,&nev,z,co,nk));
    }
  }

  SFREE(au);SFREE(ad);

  /* calculating the displacements and the stresses and storing */
  /* the results in frd format for each valid eigenmode */

  NNEW(v,double,mt**nk);
  NNEW(fn,double,mt**nk);
  NNEW(stn,double,6**nk);
  NNEW(inum,ITG,*nk);
  NNEW(stx,double,6*mi[0]**ne);
  NNEW(eme,double,6*mi[0]**ne);
  if(*ithermal>1){
    NNEW(qfn,double,3**nk);
    NNEW(qfx,double,3*mi[0]**ne);
  }
  
  if(strcmp1(&filab[261],"E   ")==0) NNEW(een,double,6**nk);
  if(strcmp1(&filab[2697],"ME  ")==0) NNEW(emn,double,6**nk);
  if(strcmp1(&filab[522],"ENER")==0) NNEW(enern,double,*nk);
  if(strcmp1(&filab[2175],"CONT")==0) NNEW(cdn,double,6**nk);

  //  NNEW(temp_array,double,neq[1]);

  lfin=0;
  for(j=0;j<nev;++j){
    lint=lfin;

    if(igreen==0){
      lfin=lfin+neq[1];
      
      /*      sum=0.;
	      for(k=0;k<neq[1];++k)
	      temp_array[k]=0.;
	      if(nasym==1){
	      FORTRAN(opas,(&neq[1],&z[lint],temp_array,adb,aub,jq,irow,nzs));
	      }else{
	      FORTRAN(op,(&neq[1],&z[lint],temp_array,adb,aub,jq,irow));
	      }
	      for(k=0;k<neq[1];++k)
	      sum+=z[lint+k]*temp_array[k];
	      for(k=0;k<neq[1];++k)
	      z[lint+k]=z[lint+k]/sqrt(sum);*/
    }else{

      /* Green function: solve the system ([K]-w**2*[M]){X}={Y} */
      
      node=nodeforc[2*j];
      idir=ndirforc[j];
      
      /* check whether degree of freedom is active */
      
      if(nactdof[mt*(node-1)+idir]==0){
	printf(" *ERROR in linstatic: degree of freedom corresponding to node %d \n and direction %d is not active: no unit force can be applied\n",node,idir);
	FORTRAN(stop,());
      }
      
      /* defining a unit force on the rhs */
      
      DMEMSET(z,0,neq[1],0.);
      z[nactdof[mt*(node-1)+idir]-1]=1.;
      
      /* solving the system */
      
      if(neq[1]>0){
	if(*isolver==0){
#ifdef SPOOLES
	  spooles_solve(z,&neq[1]);
#endif
	}
	else if(*isolver==7){
#ifdef PARDISO
	  pardiso_solve(z,&neq[1],&symmetryflag,&inputformat,&nrhs);
#endif
	}
	else if(*isolver==8){
#ifdef PASTIX
#ifdef PARDISO
	  pardiso_solve(z,&neq[1],&symmetryflag,&inputformat,&nrhs);
#else
	  if( pastix_solve(z,&neq[1],&symmetryflag,&nrhs)==-1 )
          printf(" *WARNING in arpack: solving step didn't converge! Continuing anyway\n");
#endif
#endif
	  
	}
      }
    }

    if(mei[3]==1){
      if(fwrite(&z[lint],sizeof(double),neq[1],f1)!=neq[1]){
	printf(" *ERROR saving data to the eigenvalue file...");
	exit(0);
      }
    }

    if(igreen==0){
      
      /* check whether the frequency belongs to the requested
	 interval */
      
      if(fmin>-0.5){
	if(fmin*fmin>d[j]) continue;
      }
      if(fmax>-0.5){
	if(fmax*fmax<d[j]) continue;
      }
    }

    if(*nprint>0) FORTRAN(writehe,(&j));

    NNEW(eei,double,6*mi[0]*ne0);

    /* following fields are needed if *ener==1 or kode<-100 */
    
    NNEW(vini,double,mt**nk);
    NNEW(stiini,double,6*mi[0]*ne0);
    NNEW(emeini,double,6*mi[0]*ne0);
    if(*nener==1) NNEW(enerini,double,2*mi[0]*ne0);

    DMEMSET(v,0,mt**nk,0.);
    if(*iperturb==0){
      results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
	      stx,elcon,
	      nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,ielorien,
	      norien,orab,ntmat_,t0,t0,ithermal,
	      prestr,iprestr,filab,eme,emn,een,iperturb,
	      f,fn,nactdof,&iout,qa,vold,&z[lint],
	      nodeboun,ndirboun,xboun,nboun,ipompc,
	      nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
	      &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	      xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	      &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
	      emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
	      iendset,
	      ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	      nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
	      &reltime,&ne0,thicke,shcon,nshcon,
	      sideload,xload,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
	      mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
	      islavsurf,ielprop,prop,energyini,energy,&kscale,iponoel,
              inoel,nener,orname,&network,ipobody,xbody,ibody,typeboun,
	      itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
	      islavelinv,autloc,irowtloc,jqtloc,&nboun2,
	      ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
	      labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
	      &intscheme,physcon);}
    else{
      results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
	      stx,elcon,
	      nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,ielorien,
	      norien,orab,ntmat_,t1old,t1old,ithermal,
	      prestr,iprestr,filab,eme,emn,een,iperturb,
	      f,fn,nactdof,&iout,qa,vold,&z[lint],
	      nodeboun,ndirboun,xboun,nboun,ipompc,
	      nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
	      &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	      xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	      &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
	      xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
	      ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	      nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,
	      &ne0,thicke,shcon,nshcon,
	      sideload,xload,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
	      mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
	      islavsurf,ielprop,prop,energyini,energy,&kscale,iponoel,
              inoel,nener,orname,&network,ipobody,xbody,ibody,typeboun,
	      itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
	      islavelinv,autloc,irowtloc,jqtloc,&nboun2,
	      ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
	      labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
	      &intscheme,physcon);
    }
    SFREE(eei);SFREE(stiini);SFREE(emeini);SFREE(vini);
    if(*nener==1) SFREE(enerini);

    ++*kode;
    if(igreen==0){
      if(d[j]>=0.){
	freq=sqrt(d[j])/6.283185308;
      }else{
	freq=0.;
      }
    }else{
      freq=*ttime+j;
    }

    if(strcmp1(&filab[1044],"ZZS")==0){
      NNEW(neigh,ITG,40**ne);
      NNEW(ipneigh,ITG,*nk);
    }

    if(*iperturb==0){
      frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	  kode,filab,een,t1,fn,&freq,epn,ielmat,matname,enern,xstaten,
	  nstate_,istep,&iinc,ithermal,qfn,&j,&noddiam,trab,inotr,
	  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	  mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	  thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni,nmat,
	  ielprop,prop,sti);
    }else{
      frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	  kode,filab,een,t1old,fn,&freq,epn,ielmat,matname,enern,xstaten,
	  nstate_,istep,&iinc,ithermal,qfn,&j,&noddiam,trab,inotr,
	  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	  mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	  thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni,nmat,
	  ielprop,prop,sti);
    }

    if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);}

    if(*iperturb!=0){
      for(k=0;k<*nstate_*mi[0]*(ne0+maxprevcontel);++k){
	xstate[k]=xstateini[k];
      }	  
    }
  }
  
  /*--------------------------------------------------------------------*/
  /*
    -----------
    free memory
    -----------
  */
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
#ifdef PARDISO
    pardiso_cleanup(&neq[1],&symmetryflag,&inputformat);
#endif
#endif
  }

  if(igreen==0){
    if((fmax>-0.5)&&(fmax*fmax>d[nev-1])){
      printf("\n*WARNING: not all frequencies in the requested interval might be found;\nincrease the number of requested frequencies\n");
    }
  }

  if(mei[3]==1){
    fclose(f1);
  }

  if(*iperturb!=0){
    if(ncont!=0){
      *ne=ne0;
      if(*nener==1){
	if(*mortar==1)RENEW(ener,double,mi[0]**ne*2);
      }
      RENEW(ipkon,ITG,*ne);
      RENEW(lakon,char,8**ne);
      RENEW(kon,ITG,*nkon);
      if(*norien>0){
	RENEW(ielorien,ITG,mi[2]**ne);
      }
      RENEW(ielmat,ITG,mi[2]**ne);
      SFREE(cg);SFREE(straight);
      SFREE(imastop);SFREE(itiefac);SFREE(islavnode);
      SFREE(nslavnode);SFREE(iponoels);SFREE(inoels);SFREE(imastnode);
      SFREE(nmastnode);SFREE(itietri);SFREE(koncont);SFREE(xnoels);
      SFREE(springarea);SFREE(xmastnor);

      if(*mortar==0){
	SFREE(areaslav);
      }else if(*mortar==1){
	SFREE(pmastsurf);SFREE(ipe);SFREE(ime);
	SFREE(islavact);
      }
    }
  }
  
  SFREE(adb);SFREE(aub);
  //SFREE(temp_array);

  SFREE(v);SFREE(fn);SFREE(stn);SFREE(inum);SFREE(stx);SFREE(eme);
  SFREE(z);SFREE(d);SFREE(xstiff);

  //SFREE(ipobody);

  if(*ithermal>1){SFREE(qfn);SFREE(qfx);}

  if(*nstate_!=0){SFREE(xstateini);}

  if(strcmp1(&filab[261],"E   ")==0) SFREE(een);
  if(strcmp1(&filab[2697],"ME  ")==0) SFREE(emn);
  if(strcmp1(&filab[522],"ENER")==0) SFREE(enern);
  if(strcmp1(&filab[2175],"CONT")==0) SFREE(cdn);

  if(*iperturb!=0){
    mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
    mpcinfo[3]=maxlenmpc;
  }

  if(igreen==1){
    *nmethod=13;
  }

  *irowp=irow;*enerp=ener;*xstatep=xstate;*ipkonp=ipkon;*lakonp=lakon;
  *konp=kon;*ielmatp=ielmat;*ielorienp=ielorien;

  *islavsurfp=islavsurf;*pslavsurfp=pslavsurf;*clearinip=clearini;

  return;
}

#endif
