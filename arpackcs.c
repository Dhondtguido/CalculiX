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

void arpackcs(double *co, ITG *nk, ITG **konp, ITG **ipkonp, char **lakonp,
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
	      double *alcon, ITG *nalcon, double *alzero, ITG **ielmatp,
	      ITG **ielorienp, ITG *norien, double *orab, ITG *ntmat_,
	      double *t0, double *t1, double *t1old, 
	      ITG *ithermal,double *prestr, ITG *iprestr, 
	      double *vold,ITG *iperturb, double *sti, ITG *nzs,  
	      ITG *kode, ITG *mei, double *fei,
	      char *filab,
	      ITG *iexpl, double *plicon, ITG *nplicon, double *plkcon,
	      ITG *nplkcon,
	      double **xstatep, ITG *npmat_, char *matname, ITG *mi,
	      ITG *ics, double *cs, ITG *mpcend, ITG *ncmat_,
	      ITG *nstate_, ITG *mcs, ITG *nkon,
	      char *jobnamec, char *output, char *set, ITG *nset, 
	      ITG *istartset,
	      ITG *iendset, ITG *ialset, ITG *nprint, char *prlab,
	      char *prset, ITG *nener, ITG *isolver, double *trab, 
	      ITG *inotr, ITG *ntrans, double *ttime, double *fmpc,
	      ITG *ipobody, ITG *ibody, double *xbody, ITG *nbody,
	      ITG *nevtot, double *thicke, ITG *nslavs, double *tietol, 
	      ITG *mpcinfo,ITG *ntie,ITG *istep,
	      char *tieset,ITG *nintpoint,ITG *mortar,ITG *ifacecount,
	      ITG **islavsurfp,double **pslavsurfp,double **clearinip,
	      ITG *nmat,char *typeboun,ITG *ielprop,double *prop,
	      char *orname,ITG *inewton,double *t0g,double *t1g){

  /* calls the Arnoldi Package (ARPACK) for cyclic symmetry calculations */
  
  char bmat[2]="G", which[3]="LM", howmny[2]="A",*lakont=NULL,
    description[13]="            ",fneig[132]="",filabcp[9]="        ",
    lakonl[2]=" \0",*lakon=NULL,jobnamef[396]="",*turdir=NULL,
    *labmpc2=NULL;

  ITG *inum=NULL,k,ido,ldz,iparam[11],ipntr[14],lworkl,idir,nherm=1,
    info,rvec=1,*select=NULL,lfin,j,lint,iout=1,nm,index,inode,id,i,idof,
    ielas=1,icmd=0,kk,l,nkt,icntrl,*kont=NULL,*ipkont=NULL,*inumt=NULL,
    *ielmatt=NULL,net,imag,icomplex,kkv,kk6,iinc=1,nev,ncv,kscale=1,
    mxiter,lprev,ilength,ij,i1,i2,iel,ielset,node,indexe,nope,ml1,
    *inocs=NULL,*ielcs=NULL,jj,l1,l2,ngraph,is,jrow,kflag=2,
    *inotrt=NULL,symmetryflag=0,inputformat=0,ifreebody,*iy=NULL,
    mass=1, stiffness=1, buckling=0, rhsi=0, intscheme=0,*ncocon=NULL,
    coriolis=0,iworsttime,l3,iray,mt,kkx,im,ne0,*integerglob=NULL,
    *nshcon=NULL,one=1,ncont=0,*itietri=NULL,neq2,
    *koncont=NULL,ismallsliding=0,*itiefac=NULL,*islavsurf=NULL,
    *islavnode=NULL,*imastnode=NULL,*nslavnode=NULL,*nmastnode=NULL,
    *imastop=NULL,*iponoels=NULL,*inoels=NULL,*ipe=NULL,*ime=NULL,
    mpcfree,memmpc_,icascade,maxlenmpc,nkon0,iit=-1,*irow=NULL,nasym=0,
    kmax1,kmax2,icfd=0,*inomat=NULL,*ipkon=NULL,*kon=NULL,*ielmat=NULL,
    *ielorien=NULL,*islavact=NULL,*islavsurfold=NULL,nslavs_prev_step,
    maxprevcontel,iflagact=0,*nmc=NULL,icutb=0,ialeatoric=0,
    *iponoel=NULL,*inoel=NULL,network=0,ioffr,nrhs=1,
    ioffrl,igreen=0,mscalmethod=0,kref,*jqw=NULL,*iroww=NULL,nzsw,
    *islavelinv=NULL,*irowtloc=NULL,*jqtloc=NULL,nboun2,
    *ndirboun2=NULL,*nodeboun2=NULL,nmpc2,*ipompc2=NULL,*nodempc2=NULL,
    *ikboun2=NULL,*ilboun2=NULL,*ikmpc2=NULL,*ilmpc2=NULL,mortartrafoflag=0;

  double *stn=NULL,*v=NULL,*resid=NULL,*z=NULL,*workd=NULL,*vr=NULL,
    *workl=NULL,*d=NULL,sigma,*temp_array=NULL,*vini=NULL,
    *een=NULL,cam[5],*f=NULL,*fn=NULL,qa[4],*fext=NULL,*emn=NULL,
    *epn=NULL,*stiini=NULL,*fnr=NULL,*fni=NULL,fnreal,fnimag,*emeini=NULL,
    *xstateini=NULL,theta=0,pi,*coefmpcnew=NULL,*xstiff=NULL,*vi=NULL,
    *vt=NULL,*fnt=NULL,*stnt=NULL,*eent=NULL,*cot=NULL,t[3],ctl,stl,
    *t1t=NULL,freq,*stx=NULL,*enern=NULL,*enernt=NULL,*xstaten=NULL,
    *eei=NULL,*enerini=NULL,*cocon=NULL,*qfx=NULL,*qfn=NULL,*qfnt=NULL,
    tol,fmin,fmax,xreal,ximag,*cgr=NULL,*xloadold=NULL,reltime=1.,constant,
    vreal,vimag,*stnr=NULL,*stni=NULL,stnreal,stnimag,*vmax=NULL,
    *stnmax=NULL,vl[4],stnl[6],dd,v1,v2,v3,bb,cc,al[3],cm,cn,tt,
    worstpsmax,vray[3],worstumax,p1[3],p2[3],q[3],tan[3],*springarea=NULL,
    *stxt=NULL,*eenmax=NULL,eenl[6],*emnt=NULL,*clearini=NULL,
    *doubleglob=NULL,*shcon=NULL,*cg=NULL,*straight=NULL,*cdn=NULL,
    *xmastnor=NULL,*areaslav=NULL,*dc=NULL,*di=NULL,*xnoels=NULL,
    *workev=NULL,*ener=NULL,*xstate=NULL,cdnimag,*temp_array2=NULL,
    sigmai=0,amp,ampmax,*zstorage=NULL,*au=NULL,*ad=NULL,cdnreal,
    *b=NULL,*aub=NULL,*adb=NULL,*pslavsurf=NULL,*pmastsurf=NULL,
    *cdnt=NULL,*cdnr=NULL,*cdni=NULL,*eme=NULL,alea=0.1,sum,
    *pslavsurfold=NULL,*energyini=NULL,*energy=NULL,xn[3],e1[3],e2[3],
    *smscale=NULL,*auw=NULL,*autloc=NULL,*xboun2=NULL,*coefmpc2=NULL,
    *dstorage=NULL,*distorage=NULL,*physcon=NULL;

  FILE *f1;

  /* dummy arguments for the results call */

  double *veold=NULL,*accold=NULL,bet,gam,dtime,time=1.;

  ITG *ipneigh=NULL,*neigh=NULL;
  
#ifdef SGI
  ITG token;
#endif
  
  irow=*irowp;xstate=*xstatep;ipkon=*ipkonp;lakon=*lakonp;
  kon=*konp;ielmat=*ielmatp;ielorien=*ielorienp;

  islavsurf=*islavsurfp;pslavsurf=*pslavsurfp;clearini=*clearinip;

  if(*nmethod==13){
    *nmethod=2;
    igreen=1;
  }

  if(*nener==1){NNEW(ener,double,2*mi[0]**ne);}
  
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
  
  mt=mi[1]+1;
  pi=4.*atan(1.);
  constant=180./pi;
  
  /* copying the frequency parameters */
  
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
      printf(" *ERROR in arpackcs: generalized gravity loading is not allowed in frequency calculations");
      FORTRAN(stop,());
    }
  }
  
  ne0=*ne;nkon0=*nkon;
  
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
	  DMEMSET(ener,2*mi[0]**ne,2*mi[0]*(*ne+*nslavs),0.);
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
	  DMEMSET(ener,2*mi[0]**ne,2*mi[0]*(*ne+*nintpoint),0.);
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
	      xstateini,xstate,nstate_,&icutb,&ialeatoric,jobnamef,
	      &alea,auw,jqw,iroww,&nzsw);
	  
      printf("number of contact spring elements=%" ITGFORMAT "\n\n",*ne-ne0);
	  
      /* determining the structure of the stiffness/mass matrix */
	  
      remastructar(ipompc,&coefmpc,&nodempc,nmpc,
		   &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,ikboun,ilboun,
		   labmpc,nk,&memmpc_,&icascade,&maxlenmpc,
		   kon,ipkon,lakon,ne,nactdof,icol,jq,&irow,isolver,
		   neq,nzs,nmethod,ithermal,iperturb,&mass,mi,ics,cs,
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
      
  SFREE(f);SFREE(v);SFREE(fn);SFREE(stx);SFREE(eme);SFREE(inum);
  iout=1;
  
  /* for the frequency analysis linear strain and elastic properties
     are used */
  
  iperturb[1]=0;
  
  /* determining the maximum number of sectors to be plotted */
  
  ngraph=1;
  for(j=0;j<*mcs;j++){
    if(cs[17*j+4]>ngraph) ngraph=cs[17*j+4];
  }
  
  /* assigning nodes and elements to sectors */
  
  NNEW(inocs,ITG,*nk);
  NNEW(ielcs,ITG,*ne);
  ielset=cs[12];
  if((*mcs!=1)||(ielset!=0)){
    for(i=0;i<*nk;i++) inocs[i]=-1;
    for(i=0;i<*ne;i++) ielcs[i]=-1;
  }
  
  for(i=0;i<*mcs;i++){
    is=cs[17*i+4];
    if((is==1)&&(*mcs==1)) continue;
    ielset=cs[17*i+12];
    if(ielset==0) continue;
    for(i1=istartset[ielset-1]-1;i1<iendset[ielset-1];i1++){
      if(ialset[i1]>0){
	iel=ialset[i1]-1;
	if(ipkon[iel]<0) continue;
	ielcs[iel]=i;
	indexe=ipkon[iel];
	if(strcmp1(&lakon[8*iel+3],"2")==0)nope=20;
	else if (strcmp1(&lakon[8*iel+3],"8")==0)nope=8;
	else if (strcmp1(&lakon[8*iel+3],"10")==0)nope=10;
	else if (strcmp1(&lakon[8*iel+3],"4")==0)nope=4;
	else if (strcmp1(&lakon[8*iel+3],"15")==0)nope=15;
	else if (strcmp1(&lakon[8*iel+3],"6")==0)nope=6;
	else if (strcmp1(&lakon[8*iel],"ES")==0){
	  lakonl[0]=lakon[8*iel+7];
	  nope=atoi(lakonl)+1;}
	else continue;
	      
	for(i2=0;i2<nope;++i2){
	  node=kon[indexe+i2]-1;
	  inocs[node]=i;
	}
      }
      else{
	iel=ialset[i1-2]-1;
	do{
	  iel=iel-ialset[i1];
	  if(iel>=ialset[i1-1]-1) break;
	  if(ipkon[iel]<0) continue;
	  ielcs[iel]=i;
	  indexe=ipkon[iel];
	  if(strcmp1(&lakon[8*iel+3],"2")==0)nope=20;
	  else if (strcmp1(&lakon[8*iel+3],"8")==0)nope=8;
	  else if (strcmp1(&lakon[8*iel+3],"10")==0)nope=10;
	  else if (strcmp1(&lakon[8*iel+3],"4")==0)nope=4;
	  else if (strcmp1(&lakon[8*iel+3],"15")==0)nope=15;
	  else {nope=6;}
	  for(i2=0;i2<nope;++i2){
	    node=kon[indexe+i2]-1;
	    inocs[node]=i;
	  }
	}while(1);
      }
    } 
  }
  
  /* loop over the nodal diameters */
  
  for(nm=cs[1];nm<=cs[2];++nm){
      
    NNEW(ad,double,neq[1]);
    NNEW(au,double,nzs[1]);
      
    NNEW(adb,double,neq[1]);
    NNEW(aub,double,nzs[1]);
      
    NNEW(fext,double,neq[1]);
    if(*iperturb==0){
      FORTRAN(mafillsmcs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun
			  ,nboun,
			  ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
			  nforc,nelemload,sideload,xload,nload,xbody,ipobody,
			  nbody,cgr,
			  ad,au,fext,nactdof,icol,jq,irow,&neq[1],nzl,nmethod,
			  ikmpc,ilmpc,ikboun,ilboun,
			  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
			  ielorien,norien,orab,ntmat_,
			  t0,t0,ithermal,prestr,iprestr,vold,iperturb,sti,
			  nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
			  xstiff,npmat_,&dtime,matname,mi,
			  ics,cs,&nm,ncmat_,labmpc,&mass,&stiffness,&buckling,
			  &rhsi,
			  &intscheme,mcs,&coriolis,ibody,xloadold,&reltime,
			  ielcs,veold,
			  springarea,thicke,integerglob,doubleglob,
			  tieset,istartset,iendset,ialset,ntie,&nasym,pslavsurf,
			  pmastsurf,mortar,clearini,ielprop,prop,&ne0,&kscale,
			  xstateini,xstate,nstate_,set,nset));
    }
    else{
      FORTRAN(mafillsmcs,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,
			  nboun,
			  ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
			  nforc,nelemload,sideload,xload,nload,xbody,ipobody,
			  nbody,cgr,
			  ad,au,fext,nactdof,icol,jq,irow,&neq[1],nzl,nmethod,
			  ikmpc,ilmpc,ikboun,ilboun,
			  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
			  ielorien,norien,orab,ntmat_,
			  t0,t1old,ithermal,prestr,iprestr,vold,iperturb,sti,
			  nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
			  xstiff,npmat_,&dtime,matname,mi,
			  ics,cs,&nm,ncmat_,labmpc,&mass,&stiffness,&buckling,
			  &rhsi,
			  &intscheme,mcs,&coriolis,ibody,xloadold,&reltime,
			  ielcs,veold,
			  springarea,thicke,integerglob,doubleglob,
			  tieset,istartset,iendset,ialset,ntie,&nasym,pslavsurf,
			  pmastsurf,mortar,clearini,ielprop,prop,&ne0,&kscale,
			  xstateini,xstate,nstate_,set,nset));
	  
      if(nasym==1){
	RENEW(au,double,nzs[2]+nzs[1]);
	RENEW(aub,double,nzs[2]+nzs[1]);
	symmetryflag=2;
	inputformat=1;

	FORTRAN(mafillsmcsas,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,
			      nboun,ipompc,nodempc,coefmpc,nmpc,nodeforc,
			      ndirforc,xforc,nforc,nelemload,sideload,xload,
			      nload,xbody,ipobody,nbody,cgr,ad,au,fext,nactdof,
			      icol,jq,irow,&neq[1],nzl,nmethod,ikmpc,ilmpc,
			      ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,alcon,
			      nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
			      t0,t1,ithermal,prestr,iprestr,vold,iperturb,sti,
			      nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,
			      nplkcon,xstiff,npmat_,&dtime,matname,mi,ics,cs,
			      &nm,ncmat_,labmpc,&mass,&stiffness,&buckling,
			      &rhsi,&intscheme,mcs,&coriolis,ibody,xloadold,
			      &reltime,ielcs,veold,springarea,thicke,
			      integerglob,doubleglob,tieset,istartset,iendset,
			      ialset,ntie,&nasym,nstate_,xstateini,xstate,
			      pslavsurf,pmastsurf,mortar,clearini,ielprop,prop,
			      &ne0,&kscale,set,nset));
	      
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
	  nstate_,istep,&iinc,ithermal,qfn,&j,&nm,trab,inotr,
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
      spooles_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],nzs,&symmetryflag,
		     &inputformat,&nzs[2]);
#else
      printf(" *ERROR in arpackcs: the SPOOLES library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==4){
#ifdef SGI
      token=1;
      sgi_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],nzs,token);
#else
      printf(" *ERROR in arpackcs: the SGI library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==5){
#ifdef TAUCS
      tau_factor(ad,&au,adb,aub,&sigma,icol,&irow,&neq[1],nzs);
#else
      printf(" *ERROR in arpackcs: the TAUCS library is not linked\n\n");
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
      pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],nzs,
		     &symmetryflag,&inputformat,jq,&nzs[2]);
#else
      printf(" *ERROR in arpack: the PARDISO library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==8){
#ifdef PASTIX
#ifdef PARDISO
      pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],nzs,
		     &symmetryflag,&inputformat,jq,&nzs[2]);
#else
      pastix_factor_main(ad,au,adb,aub,&sigma,icol,irow,&neq[1],nzs,
		    &symmetryflag,&inputformat,jq,&nzs[2]);
#endif
#else
      printf(" *ERROR in arpack: the PASTIX library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
      
    //      SFREE(au);SFREE(ad);
      
    /* calculating the eigenvalues and eigenmodes */

    if(igreen==0){
      
      printf(" Calculating the eigenvalues and the eigenmodes\n\n");
	
      ido=0;
      ldz=neq[1];
      for(k=0;k<11;k++)iparam[k]=0;
      iparam[0]=1;
      iparam[2]=mxiter;
      iparam[3]=1;
      iparam[6]=3;
	
      info=0;
	
      NNEW(resid,double,neq[1]);
      NNEW(z,double,(long long)ncv*neq[1]);
      NNEW(workd,double,3*neq[1]);
	
      if(nasym==1){
	lworkl=3*ncv*(2+ncv);
	NNEW(workl,double,lworkl);
	FORTRAN(dnaupd,(&ido,bmat,&neq[1],which,&nev,&tol,resid,&ncv,z,&ldz,iparam,ipntr,workd,
			workl,&lworkl,&info));
      }else{
	lworkl=ncv*(8+ncv);
	NNEW(workl,double,lworkl);
	FORTRAN(dsaupd,(&ido,bmat,&neq[1],which,&nev,&tol,resid,&ncv,z,&ldz,
			iparam,ipntr,workd,workl,&lworkl,&info));
      }
	
      NNEW(temp_array,double,neq[1]);
	
      while((ido==-1)||(ido==1)||(ido==2)){
	if(ido==-1){
	  if(nasym==1){
	    FORTRAN(opas,(&neq[1],&workd[ipntr[0]-1],temp_array,adb,aub,jq,irow,nzs));
	  }else{
	    FORTRAN(op,(&neq[1],&workd[ipntr[0]-1],temp_array,adb,aub,
			jq,irow));
	  }
	}
	  
	/* solve the linear equation system  */
	  
	if((ido==-1)||(ido==1)){
	    
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
	      if(pastix_solve(temp_array,&neq[1],&symmetryflag,&nrhs)==-1)
		printf(" *WARNING in arpackcs: solving step did not converge! Continuing anyway!\n");
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
	      if(pastix_solve(&workd[ipntr[2]-1],&neq[1],&symmetryflag,&nrhs)==-1)
		printf(" *WARNING in arpackcs: solving step did not converge! Continuing anyway!\n");
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
	    FORTRAN(op,(neq,&workd[ipntr[0]-1],&workd[ipntr[1]-1],
			adb,aub,jq,irow));
	  }
	}
	  
	if(nasym==1){
	  FORTRAN(dnaupd,(&ido,bmat,&neq[1],which,&nev,&tol,resid,&ncv,z,&ldz,
			  iparam,ipntr,workd,workl,&lworkl,&info));
	}else{
	  FORTRAN(dsaupd,(&ido,bmat,&neq[1],which,&nev,&tol,resid,&ncv,
			  z,&ldz,iparam,ipntr,workd,workl,&lworkl,&info));
	}
      }
	
      if(info!=0){
	printf(" *ERROR in d[n,s]aupd: info=%" ITGFORMAT "\n",info);
	printf("       # of converged eigenvalues=%" ITGFORMAT "\n\n",iparam[4]);
      }         
	
      NNEW(temp_array2,double,nev);
      NNEW(iy,ITG,nev);
	
      NNEW(select,ITG,ncv);
	
      if(nasym==1){
	NNEW(d,double,nev+1);
	NNEW(di,double,nev+1);
	NNEW(workev,double,3*ncv);
	FORTRAN(dneupd,(&rvec,howmny,select,d,di,z,&ldz,&sigma,&sigmai,
			workev,bmat,&neq[1],which,&nev,&tol,resid,
			&ncv,z,&ldz,iparam,ipntr,workd,workl,&lworkl,&info));
	SFREE(workev);

	/* storing the size of the eigenfrequencies */

	for(j=0;j<nev;j++){temp_array2[j]=sqrt(sqrt(d[j]*d[j]+di[j]*di[j]));iy[j]=j;}
      }else{
	NNEW(d,double,nev);
	FORTRAN(dseupd,(&rvec,howmny,select,d,z,&ldz,&sigma,bmat,&neq[1],which,&nev,
			&tol,resid,&ncv,z,&ldz,iparam,ipntr,workd,workl,&lworkl,&info));
	for(j=0;j<nev;j++){temp_array2[j]=d[j];iy[j]=j;}
      }
      SFREE(select);SFREE(resid);SFREE(workd);SFREE(workl);
	
      if(info!=0){
	printf(" *ERROR in d[n,s]eupd: info=%" ITGFORMAT "\n",info);
      }         

      /* checking the ascending order of the size of the eigenfrequencies */

      FORTRAN(dsort,(temp_array2,iy,&nev,&kflag));
      for(j=0;j<nev;j++){
	if(iy[j]!=j){

	  /* if two identical eigenfrequencies are swapped, do nothing */
	  
	  if((iy[iy[j]]==j)&&(fabs((temp_array2[j]-temp_array2[iy[j]])/
				   temp_array2[j])<1.e-5)) continue;
	  printf(" Eigenfrequencies are not in ascending order; sorting;\n\n");
	  NNEW(zstorage,double,nev*neq[1]);
	  if(nasym==1){
	    NNEW(dstorage,double,nev);NNEW(distorage,double,nev);
	  }
	  for(k=0;k<nev;k++){
	    memcpy(&zstorage[k*neq[1]],&z[iy[k]*neq[1]],sizeof(double)*neq[1]);
	    if(nasym==1){
	      dstorage[k]=d[iy[k]];
	      distorage[k]=di[iy[k]];
	    }else{d[k]=temp_array2[k];}
	  }
	  if(nasym==1){
	    memcpy(d,dstorage,sizeof(double)*nev);
	    memcpy(di,distorage,sizeof(double)*nev);
	    SFREE(dstorage);SFREE(distorage);
	  }
	  memcpy(z,zstorage,sizeof(double)*nev*neq[1]);SFREE(zstorage);
	  break;
	}
      }
      SFREE(temp_array2);SFREE(iy);
      
      /* writing the eigenvalues and eigenfrequencies to file */
      
      if(nasym==1){
	  
	/* storing the eigenfrequency (= sqrt(eigenvalue)) as complex number */
	  
	NNEW(dc,double,2*nev);
	NNEW(nmc,ITG,nev);
	for(j=0;j<nev;j++){
	  dc[2*j]=sqrt(sqrt(d[j]*d[j]+di[j]*di[j])+d[j])/sqrt(2.);
	  dc[2*j+1]=sqrt(sqrt(d[j]*d[j]+di[j]*di[j])-d[j])/sqrt(2.);
	  if(di[j]<0.) dc[2*j+1]=-dc[2*j+1];
	  nmc[j]=nm;
	}
	FORTRAN(writeevcscomplex,(dc,&nev,nmc,&fmin,&fmax));
	SFREE(di);SFREE(dc);SFREE(nmc);
      }else{
	FORTRAN(writeevcs,(d,&nev,&nm,&fmin,&fmax));
      }

      /* check the normalization of the eigenmodes */

      for(j=0;j<nev;j++){
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
	printf("U^T*M*U=%f for eigenmode %d\n",sum,j+1);

	/* normalizing the eigenmode (if not yet normalized by ARPACK */

	if(fabs(sum-1.)>1.e-5){
	  printf("normalizing the eigenmode \n");
	  for(k=0;k<neq[1];k++){z[kref+k]/=sqrt(sum);}
	}
	  
      }
	
      SFREE(temp_array);
	
      /* for double eigenmodes: the eigenmode for which the largest
	 amplitude has the lowest dof comes first */
	
      neq2=neq[1]/2;
      for (j=0;j<nev;j+=2){
	  
	ampmax=0.;kmax1=0;
	for(k=0;k<neq2;k++){
	  amp=z[j*neq[1]+k]*z[j*neq[1]+k]+
	    z[j*neq[1]+neq2+k]*z[j*neq[1]+neq2+k];
	  if(amp>ampmax){ampmax=amp;kmax1=k;}
	}
	  
	ampmax=0.;kmax2=0;
	for(k=0;k<neq2;k++){
	  amp=z[(j+1)*neq[1]+k]*z[(j+1)*neq[1]+k]+
	    z[(j+1)*neq[1]+neq2+k]*z[(j+1)*neq[1]+neq2+k];
	  if(amp>ampmax){ampmax=amp;kmax2=k;}
	}
	  
	if(kmax2<kmax1){
	  printf("exchange mode %d and mode %d!\n",j+1,j+2);
	  NNEW(zstorage,double,neq[1]);
	  memcpy(zstorage,&z[j*neq[1]],sizeof(double)*neq[1]);
	  memcpy(&z[j*neq[1]],&z[(j+1)*neq[1]],sizeof(double)*neq[1]);
	  memcpy(&z[(j+1)*neq[1]],zstorage,sizeof(double)*neq[1]);
	  SFREE(zstorage);
	}
      }
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
	  
      /* the first time the file is erased before writing, all subsequent
	 times the data is appended */
	  
      if(*nevtot==0){
	if((f1=fopen(fneig,"wb"))==NULL){
	  printf(" *ERROR in arpack: cannot open eigenvalue file for writing...");
	  exit(0);
	}
	      
	/* storing a one as indication that this was a
	   cyclic symmetry calculation */
	      
	if(fwrite(&one,sizeof(ITG),1,f1)!=1){
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
	      
      }else{
	if((f1=fopen(fneig,"ab"))==NULL){
	  printf(" *ERROR in arpack: cannot open eigenvalue file for writing...");
	  exit(0);
	}
      }
	  
      /* nodal diameter */
	  
      if(fwrite(&nm,sizeof(ITG),1,f1)!=1){
	printf(" *ERROR saving the nodal diameter to the eigenvalue file...");
	exit(0);
      }
	  
      /* number of eigenfrequencies */
	  
      if(fwrite(&nev,sizeof(ITG),1,f1)!=1){
	printf(" *ERROR saving the number of eigenfrequencies to the eigenvalue file...");
	exit(0);
      }
	  
      /* the eigenfrequencies are stored as radians/time */
	  
      if(fwrite(d,sizeof(double),nev,f1)!=nev){
	printf(" *ERROR saving the eigenfrequencies to the eigenvalue file...");
	exit(0);
      }
      if(*nevtot==0){
	      
	/* stiffness matrix */
	      
	if(fwrite(ad,sizeof(double),neq[1],f1)!=neq[1]){
	  printf(" *ERROR saving the diagonal of the stiffness matrix to the eigenvalue file...");
	  exit(0);
	}
	if(fwrite(au,sizeof(double),nzs[1],f1)!=nzs[1]){
	  printf(" *ERROR saving the off-diagonal terms of the stiffness matrix to the eigenvalue file...");
	  exit(0);
	}
	      
	/* mass matrix */
	      
	if(fwrite(adb,sizeof(double),neq[1],f1)!=neq[1]){
	  printf(" *ERROR saving the diagonal of the mass matrix to the eigenvalue file...");
	  exit(0);
	}
	if(fwrite(aub,sizeof(double),nzs[1],f1)!=nzs[1]){
	  printf(" *ERROR saving the off-diagonal terms of the mass matrix to the eigenvalue file...");
	  exit(0);
	}
      }

      //	  SFREE(ad);SFREE(au);
    }

    if(igreen==0){
	
      /* calculating the participation factors and the relative effective
	 modal mass */
	
      if(*ithermal!=2){
	FORTRAN(effectivemodalmass,(neq,nactdof,mi,adb,aub,jq,irow,&nev,z,co,nk));
      }
    }
      
    /* calculating the displacements and the stresses and storing */
    /* the results in frd format for each valid eigenmode */

    /* for energy calculations in other sectors the stress and the
       mechanical strain have to be calculated */

    if((ngraph>1)&&(strcmp1(&filab[522],"ENER")==0)){
      strcpy1(&filabcp[0],&filab[174],4);
      strcpy1(&filabcp[4],&filab[2697],4);
      strcpy1(&filab[174],"S   ",4);
      strcpy1(&filab[2697],"ME  ",4);
    }

      
    NNEW(v,double,2*mt**nk);
    NNEW(fn,double,2*mt**nk);
    NNEW(inum,ITG,*nk);
    NNEW(stx,double,2*6*mi[0]**ne);
    NNEW(eme,double,2*6*mi[0]**ne);

    if((strcmp1(&filab[174],"S   ")==0)||(strcmp1(&filab[1653],"MAXS")==0)|| 
       (strcmp1(&filab[1479],"PHS ")==0)||(strcmp1(&filab[1044],"ZZS ")==0)||
       (strcmp1(&filab[1044],"ERR ")==0)) 
      NNEW(stn,double,12**nk);
      
    if((strcmp1(&filab[261],"E   ")==0)||(strcmp1(&filab[2523],"MAXE")==0)) 
      NNEW(een,double,12**nk);
    if(strcmp1(&filab[522],"ENER")==0) NNEW(enern,double,2**nk);
    if(strcmp1(&filab[2697],"ME  ")==0) NNEW(emn,double,12**nk);
    if(((strcmp1(&filab[2175],"CONT")==0)||(strcmp1(&filab[3915],"PCON")==0))
       &&(*mortar==1)) NNEW(cdn,double,12**nk);
      
    NNEW(coefmpcnew,double,*mpcend);
      
    /* creating total fields for ngraph segments */

    NNEW(cot,double,3**nk*ngraph);
    if(*ntrans>0){NNEW(inotrt,ITG,2**nk*ngraph);}
    if((strcmp1(&filab[0],"U  ")==0)||(strcmp1(&filab[870],"PU  ")==0))
	  
      // real and imaginary part of the displacements
	  
      NNEW(vt,double,2*mt**nk*ngraph);
    if(strcmp1(&filab[87],"NT  ")==0)
      NNEW(t1t,double,*nk*ngraph);
    if((strcmp1(&filab[174],"S   ")==0)||(strcmp1(&filab[1479],"PHS ")==0)||
       (strcmp1(&filab[1044],"ZZS ")==0)||(strcmp1(&filab[1044],"ERR ")==0))
	  
      // real and imaginary part of the stresses
	  
      NNEW(stnt,double,2*6**nk*ngraph);
    if(strcmp1(&filab[261],"E   ")==0) NNEW(eent,double,2*6**nk*ngraph);
    if((strcmp1(&filab[348],"RF  ")==0)||(strcmp1(&filab[2610],"PRF ")==0))
	  
      // real and imaginary part of the forces
	  
      NNEW(fnt,double,2*mt**nk*ngraph);
    if(strcmp1(&filab[2697],"ME  ")==0) NNEW(emnt,double,2*6**nk*ngraph);
    if(((strcmp1(&filab[2175],"CONT")==0)||(strcmp1(&filab[3915],"PCON")==0))
       &&(*mortar==1)) NNEW(cdnt,double,2*6**nk*ngraph);
    if(strcmp1(&filab[522],"ENER")==0)
	  
      // real and imaginary part of the internal energy
	  
      NNEW(enernt,double,*nk*ngraph);

    // stresses at the integration points for the error estimator and contact conditions

    if((strcmp1(&filab[1044],"ZZS ")==0)||(strcmp1(&filab[1044],"ERR ")==0)||
       ((strcmp1(&filab[2175],"CONT")==0)&&(*mortar==0)))
      NNEW(stxt,double,2*6*mi[0]**ne*ngraph);
      
    NNEW(kont,ITG,*nkon*ngraph);
    NNEW(ipkont,ITG,*ne*ngraph);
    for(l=0;l<*ne*ngraph;l++)ipkont[l]=-1;
    NNEW(lakont,char,8**ne*ngraph);
    NNEW(inumt,ITG,*nk*ngraph);
    NNEW(ielmatt,ITG,mi[2]**ne*ngraph);
      
    nkt=ngraph**nk;
    net=ngraph**ne;
      
    /* copying the coordinates of the first sector */
      
    for(l=0;l<3**nk;l++){cot[l]=co[l];}
    if(*ntrans>0){for(l=0;l<*nk;l++){inotrt[2*l]=inotr[2*l];}}
    for(l=0;l<*nkon;l++){kont[l]=kon[l];}
    for(l=0;l<*ne;l++){ipkont[l]=ipkon[l];}
    for(l=0;l<8**ne;l++){lakont[l]=lakon[l];}
    for(l=0;l<*ne;l++){ielmatt[mi[2]*l]=ielmat[mi[2]*l];}
      
    /* generating the coordinates for the other sectors */
      
    icntrl=1;
      
    FORTRAN(rectcyl,(cot,v,fn,stn,qfn,een,cs,nk,&icntrl,t,filab,&imag,mi,emn));
      
    for(jj=0;jj<*mcs;jj++){
      is=cs[17*jj+4];
      for(i=1;i<is;i++){
	      
	theta=i*2.*pi/cs[17*jj];
	      
	for(l=0;l<*nk;l++){
	  if(inocs[l]==jj){
	    cot[3*l+i*3**nk]=cot[3*l];
	    cot[1+3*l+i*3**nk]=cot[1+3*l]+theta;
	    cot[2+3*l+i*3**nk]=cot[2+3*l];
	    if(*ntrans>0){inotrt[2*l+i*2**nk]=inotrt[2*l];}
	  }
	}
	for(l=0;l<*nkon;l++){kont[l+i**nkon]=kon[l]+i**nk;}
	for(l=0;l<*ne;l++){
	  if(ielcs[l]==jj){
	    if(ipkon[l]>=0){
	      ipkont[l+i**ne]=ipkon[l]+i**nkon;
	      ielmatt[mi[2]*(l+i**ne)]=ielmat[mi[2]*l];
	      for(l1=0;l1<8;l1++){
		l2=8*l+l1;
		lakont[l2+i*8**ne]=lakon[l2];
	      }
	    }
	  }
	}
      }
    }
      
    icntrl=-1;
      
    FORTRAN(rectcyl,(cot,vt,fnt,stnt,qfnt,eent,cs,&nkt,&icntrl,t,filab,
		     &imag,mi,emnt));
      
    /* check that the tensor fields which are extrapolated from the
       integration points are requested in global coordinates */
      
    if(strcmp1(&filab[174],"S   ")==0){
      if((strcmp1(&filab[179],"L")==0)&&(*norien>0)){
	printf("\n*WARNING in arpackcs: element fields in cyclic symmetry calculations\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
	strcpy1(&filab[179],"G",1);
      }
    }
      
    if(strcmp1(&filab[261],"E   ")==0){
      if((strcmp1(&filab[266],"L")==0)&&(*norien>0)){
	printf("\n*WARNING in arpackcs: element fields in cyclic symmetry calculation\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
	strcpy1(&filab[266],"G",1);
      }
    }
      
    if(strcmp1(&filab[1479],"PHS ")==0){
      if((strcmp1(&filab[1484],"L")==0)&&(*norien>0)){
	printf("\n*WARNING in arpackcs: element fields in cyclic symmetry calculation\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
	strcpy1(&filab[1484],"G",1);
      }
    }
      
    if(strcmp1(&filab[1653],"MAXS")==0){
      if((strcmp1(&filab[1658],"L")==0)&&(*norien>0)){
	printf("\n*WARNING in arpackcs: element fields in cyclic symmetry calculation\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
	strcpy1(&filab[1658],"G",1);
      }
    }   
      
    if(strcmp1(&filab[2523],"MAXE")==0){
      if((strcmp1(&filab[2528],"L")==0)&&(*norien>0)){
	printf("\n*WARNING in arpackcs: element fields in cyclic symmetry calculation\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
	strcpy1(&filab[1658],"G",1);
      }
    }   
      
    /* allocating fields for magnitude and phase information of
       displacements and stresses */
      
    if(strcmp1(&filab[870],"PU")==0){
      NNEW(vr,double,mt*nkt);
      NNEW(vi,double,mt*nkt);
    }
      
    if(strcmp1(&filab[1479],"PHS")==0){
      NNEW(stnr,double,6*nkt);
      NNEW(stni,double,6*nkt);
    }
      
    if(strcmp1(&filab[2610],"PRF")==0){
      NNEW(fnr,double,mt*nkt);
      NNEW(fni,double,mt*nkt);
    }
      
    if(strcmp1(&filab[3915],"PCON")==0){
      NNEW(cdnr,double,6*nkt);
      NNEW(cdni,double,6*nkt);
    }
      
    if(strcmp1(&filab[1566],"MAXU")==0){
      NNEW(vmax,double,4*nkt);
    }
      
    if(strcmp1(&filab[1653],"MAXS")==0){
      NNEW(stnmax,double,7*nkt);
    }
      
    if(strcmp1(&filab[2523],"MAXE")==0){
      NNEW(eenmax,double,7*nkt);
    }

    /* determining a local coordinate system based on the rotation axis */
      
    FORTRAN(localaxescs,(cs,mcs,e1,e2,xn));
    NNEW(turdir,char,nev);
      
    /* start of output calculations */
      
    lfin=0;
    for(j=0;j<nev;++j){
      lint=lfin;
      if(igreen==0){
	lfin=lfin+neq[1];
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
	    if(pastix_solve(z,&neq[1],&symmetryflag,&nrhs)==-1)
		printf(" *WARNING in arpackcs: solving step did not converge! Continuing anyway!\n");
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
	  
      if(*nprint>0)FORTRAN(writehe,(&j));
	  
      NNEW(eei,double,6*mi[0]*ne0);

      /* following fields are needed if *ener==1 or kode<-100 */
    
      NNEW(vini,double,mt**nk);
      NNEW(stiini,double,6*mi[0]*ne0);
      NNEW(emeini,double,6*mi[0]*ne0);
      if(*nener==1) NNEW(enerini,double,2*mi[0]*ne0);
	  
      DMEMSET(v,0,2*mt**nk,0.);
	  
      /* calculating the strains, stresses... (real and imaginary part) for
	 one specific eigenvalue */

      for(k=0;k<neq[1];k+=neq[1]/2){
	      
	if(k==0) {kk=0;kkv=0;kk6=0;kkx=0;if(*nprint>0)FORTRAN(writere,());}
	else {kk=*nk;kkv=mt**nk;kk6=6**nk;kkx=6*mi[0]**ne;
	  if(*nprint>0)FORTRAN(writeim,());}
	      
	/* generating the cyclic MPC's (needed for nodal diameters
	   different from 0 */
	      
	for(i=0;i<*nmpc;i++){
	  index=ipompc[i]-1;
	  /* check whether thermal mpc */
	  if(nodempc[3*index+1]==0) continue;
	  coefmpcnew[index]=coefmpc[index];
	  while(1){
	    index=nodempc[3*index+2];
	    if(index==0) break;
	    index--;
		      
	    icomplex=0;
	    inode=nodempc[3*index];
	    if(strcmp1(&labmpc[20*i],"CYCLIC")==0){
	      icomplex=atoi(&labmpc[20*i+6]);}
	    else if(strcmp1(&labmpc[20*i],"SUBCYCLIC")==0){
	      for(ij=0;ij<*mcs;ij++){
		lprev=cs[ij*17+13];
		ilength=cs[ij*17+3];
		FORTRAN(nident,(&ics[lprev],&inode,&ilength,&id));
		if(id!=0){
		  if(ics[lprev+id-1]==inode){icomplex=ij+1;break;}
		}
	      }
	    }
		      
	    if(icomplex!=0){
	      idir=nodempc[3*index+1];
	      idof=nactdof[mt*(inode-1)+idir]-1;
	      if(idof<=-1){xreal=1.;ximag=1.;}
	      else{xreal=z[lint+idof];ximag=z[lint+idof+neq[1]/2];}
	      if(k==0) {
		if(fabs(xreal)<1.e-30)xreal=1.e-30;
		coefmpcnew[index]=coefmpc[index]*
		  (cs[17*(icomplex-1)+14]+ximag/xreal*cs[17*(icomplex-1)+15]);}
	      else {
		if(fabs(ximag)<1.e-30)ximag=1.e-30;
		coefmpcnew[index]=coefmpc[index]*
		  (cs[17*(icomplex-1)+14]-xreal/ximag*cs[17*(icomplex-1)+15]);}
	    }
	    else{coefmpcnew[index]=coefmpc[index];}
	  }
	}
	      
	if(*iperturb==0){
	  results(co,nk,kon,ipkon,lakon,ne,&v[kkv],&stn[kk6],inum,
		  &stx[kkx],elcon,
		  nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,ielorien,
		  norien,orab,ntmat_,t0,t0,ithermal,
		  prestr,iprestr,filab,&eme[kkx],&emn[kk6],&een[kk6],iperturb,
		  f,&fn[kkv],nactdof,&iout,qa,vold,&z[lint+k],
		  nodeboun,ndirboun,xboun,nboun,ipompc,
		  nodempc,coefmpcnew,labmpc,nmpc,nmethod,cam,&neq[1],veold,
		  accold,
		  &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		  xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
		  ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,&enern[kk],
		  emeini,
		  xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
		  ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
		  nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,
		  &ne0,thicke,shcon,nshcon,
		  sideload,xload,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		  mortar,islavact,&cdn[kk6],islavnode,nslavnode,ntie,clearini,
		  islavsurf,ielprop,prop,energyini,energy,&kscale,iponoel,
		  inoel,nener,orname,&network,ipobody,xbody,ibody,typeboun,
		  itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
		  islavelinv,autloc,irowtloc,jqtloc,&nboun2,
		  ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
		  labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
		  &intscheme,physcon);}
	else{
	  results(co,nk,kon,ipkon,lakon,ne,&v[kkv],&stn[kk6],inum,
		  &stx[kkx],elcon,
		  nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,ielorien,
		  norien,orab,ntmat_,t1old,t1old,ithermal,
		  prestr,iprestr,filab,&eme[kkx],&emn[kk6],&een[kk6],iperturb,
		  f,&fn[kkv],nactdof,&iout,qa,vold,&z[lint+k],
		  nodeboun,ndirboun,xboun,nboun,ipompc,
		  nodempc,coefmpcnew,labmpc,nmpc,nmethod,cam,&neq[1],veold,
		  accold,
		  &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		  xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
		  ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,&enern[kk],
		  emeini,
		  xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
		  ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
		  nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,
		  &ne0,thicke,shcon,nshcon,
		  sideload,xload,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		  mortar,islavact,&cdn[kk6],islavnode,nslavnode,ntie,clearini,
		  islavsurf,ielprop,prop,energyini,energy,&kscale,iponoel,
		  inoel,nener,orname,&network,ipobody,xbody,ibody,typeboun,
		  itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
		  islavelinv,autloc,irowtloc,jqtloc,&nboun2,
		  ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
		  labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
		  &intscheme,physcon);
	}
	      
      }
      SFREE(eei);SFREE(stiini);SFREE(emeini);SFREE(vini);
      if(*nener==1) SFREE(enerini);

      /* determining the turning direction */
	  
      FORTRAN(turningdirection,(v,e1,e2,xn,mi,nk,&turdir[j],lakon,ipkon,
				kon,ne,co));
	  
      if(strcmp1(&filab[1566],"MAXU")==0){
	      
	/* determining the ray vector; the components of the
	   ray vector are the coordinates of the node in node set
	   RAY */
	      
	iray=0;
	for(i=0;i<*nset;i++){
	  if(strcmp1(&set[81*i],"RAYN")==0){
	    iray=ialset[istartset[i]-1];
	    vray[0]=co[3*iray-3];
	    vray[1]=co[3*iray-2];
	    vray[2]=co[3*iray-1];
	    break;
	  }
	}
	if(iray==0){
	  printf("/n *ERROR in arpackcs: no light ray vector/n/n");
	  FORTRAN(stop,());
	}
	      
	/* initialization */
	      
	for(l1=0;l1<4**nk;l1++){vmax[l1]=0.;}
	      
	/* vector p1 is a point on the rotation axis
	   vector p2 is a unit vector along the axis */
	      
	for(l2=0;l2<3;l2++){p1[l2]=cs[5+l2];}
	for(l2=0;l2<3;l2++){p2[l2]=cs[8+l2]-p1[l2];}
	dd=sqrt(p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]);
	for(l2=0;l2<3;l2++){p2[l2]/=dd;}
	      
	/* determine the time for the worst displacement
	   orthogonal to a give light ray vector ; */
	      
	for(l1=0;l1<*nk;l1++){
		  
	  /*  determining a vector through node (l1+1) and
	      orthogonal to the rotation axis */
		  
	  for(l2=0;l2<3;l2++){q[l2]=co[3*l1+l2]-p1[l2];}
	  dd=q[0]*p2[0]+q[1]*p2[1]+q[2]*p2[2];
	  for(l2=0;l2<3;l2++){q[l2]-=dd*p2[l2];}
		  
	  /* determining a vector tan orthogonal to vector q
	     and the ray vector */
		  
	  tan[0]=q[1]*vray[2]-q[2]*vray[1];
	  tan[1]=q[2]*vray[0]-q[0]*vray[2];
	  tan[2]=q[0]*vray[1]-q[1]*vray[0];
		  
	  printf("tangent= %" ITGFORMAT ",%e,%e,%e\n",l1,tan[0],tan[1],tan[2]);
		  
	  worstumax=0.;
	  iworsttime=0;
	  for(l3=0;l3<360;l3++){
	    ctl=cos(l3/constant);
	    stl=sin(l3/constant);
	    for(l2=1;l2<4;l2++){
	      l=mt*l1+l2;
	      vl[l2]=ctl*v[l]-stl*v[l+mt**nk];
	    }
		      
	    /* displacement component along the tangent vector
	       (no absolute value!) */
		      
	    dd=vl[1]*tan[0]+vl[2]*tan[1]+vl[3]*tan[2];
	    if(dd>worstumax){
	      worstumax=dd;
	      iworsttime=l3;
	    }
	  }
	  ctl=cos(iworsttime/constant);
	  stl=sin(iworsttime/constant);
	  for(l2=1;l2<4;l2++){
	    l=mt*l1+l2;
	    vl[l2]=ctl*v[l]-stl*v[l+mt**nk];
	  }
	  vmax[4*l1]=1.*iworsttime;
	  vmax[4*l1+1]=vl[1];
	  vmax[4*l1+2]=vl[2];
	  vmax[4*l1+3]=vl[3];
		  
	}
      }
	  
      /* determine the worst principal stress anywhere
	 in the structure as a function of time; 
	 the worst principal stress is the maximum
	 of the absolute value of the principal stresses */
	  
      if(strcmp1(&filab[1653],"MAXS")==0){
	      
	/* determining the set of nodes for the 
	   worst principal stress calculation */
	      
	ielset=0;
	for(i=0;i<*nset;i++){
	  if(strcmp1(&set[81*i],"STRESSDOMAINN")==0){
	    ielset=i+1;
	    break;
	  }
	}
	if(ielset==0){
	  printf("\n *ERROR in arpackcs: no node set for MAXS\n");
	  printf("       (must have the name STRESSDOMAIN)\n\n");
	  FORTRAN(stop,());
	}
	      
	for(i1=istartset[ielset-1]-1;i1<iendset[ielset-1];i1++){
	  if(ialset[i1]>0){
	    l1=ialset[i1]-1;
		      
	    worstpsmax=0.;
	    for(l3=0;l3<360;l3++){
	      ctl=cos(l3/constant);
	      stl=sin(l3/constant);
	      for(l2=0;l2<6;l2++){
		l=6*l1+l2;
		stnl[l2]=ctl*stn[l]-stl*stn[l+6**nk];
	      }
			  
	      /* determining the eigenvalues */
			  
	      v1=stnl[0]+stnl[1]+stnl[2];
	      v2=stnl[1]*stnl[2]+stnl[0]*stnl[2]+stnl[0]*stnl[1]-
		(stnl[5]*stnl[5]+stnl[4]*stnl[4]+stnl[3]*stnl[3]);
	      v3=stnl[0]*(stnl[1]*stnl[2]-stnl[5]*stnl[5])
		-stnl[3]*(stnl[3]*stnl[2]-stnl[4]*stnl[5])
		+stnl[4]*(stnl[3]*stnl[5]-stnl[4]*stnl[1]);
	      bb=v2-v1*v1/3.;
	      cc=-2.*v1*v1*v1/27.+v1*v2/3.-v3;
	      if(fabs(bb)<=1.e-10){
		if(fabs(cc)>1.e-10){
		  al[0]=-pow(cc,(1./3.));
		}else{
		  al[0]=0.;
		}
		al[1]=al[0];
		al[2]=al[0];
	      }else{
		cm=2.*sqrt(-bb/3.);
		cn=3.*cc/(cm*bb);
		if(fabs(cn)>1.){
		  if(cn>1.){
		    cn=1.;
		  }else{
		    cn=-1.;
		  }
		}
		tt=(atan2(sqrt(1.-cn*cn),cn))/3.;
		al[0]=cm*cos(tt);
		al[1]=cm*cos(tt+2.*pi/3.);
		al[2]=cm*cos(tt+4.*pi/3.);
	      }
	      for(l2=0;l2<3;l2++){
		al[l2]+=v1/3.;
	      }
	      dd=fabs(al[0]);
	      if(fabs(al[1])>dd) dd=fabs(al[1]);
	      if(fabs(al[2])>dd) dd=fabs(al[2]);
	      if(dd>worstpsmax){
		worstpsmax=dd;
		stnmax[7*l1]=dd;
		for(l2=1;l2<7;l2++){
		  stnmax[7*l1+l2]=stnl[l2-1];
		}
	      }
	    }
		      
	  }else{
	    l1=ialset[i1-2]-1;
	    do{
	      l1=l1-ialset[i1];
	      if(l1>=ialset[i1-1]-1) break;
			  
	      worstpsmax=0.;
	      for(l3=0;l3<360;l3++){
		ctl=cos(l3/constant);
		stl=sin(l3/constant);
		for(l2=0;l2<6;l2++){
		  l=6*l1+l2;
		  stnl[l2]=ctl*stn[l]-stl*stn[l+6**nk];
		}
			      
		/* determining the eigenvalues */
			      
		v1=stnl[0]+stnl[1]+stnl[2];
		v2=stnl[1]*stnl[2]+stnl[0]*stnl[2]+stnl[0]*stnl[1]-
		  (stnl[5]*stnl[5]+stnl[4]*stnl[4]+stnl[3]*stnl[3]);
		v3=stnl[0]*(stnl[1]*stnl[2]-stnl[5]*stnl[5])
		  -stnl[3]*(stnl[3]*stnl[2]-stnl[4]*stnl[5])
		  +stnl[4]*(stnl[3]*stnl[5]-stnl[4]*stnl[1]);
		bb=v2-v1*v1/3.;
		cc=-2.*v1*v1*v1/27.+v1*v2/3.-v3;
		if(fabs(bb)<=1.e-10){
		  if(fabs(cc)>1.e-10){
		    al[0]=-pow(cc,(1./3.));
		  }else{
		    al[0]=0.;
		  }
		  al[1]=al[0];
		  al[2]=al[0];
		}else{
		  cm=2.*sqrt(-bb/3.);
		  cn=3.*cc/(cm*bb);
		  if(fabs(cn)>1.){
		    if(cn>1.){
		      cn=1.;
		    }else{
		      cn=-1.;
		    }
		  }
		  tt=(atan2(sqrt(1.-cn*cn),cn))/3.;
		  al[0]=cm*cos(tt);
		  al[1]=cm*cos(tt+2.*pi/3.);
		  al[2]=cm*cos(tt+4.*pi/3.);
		}
		for(l2=0;l2<3;l2++){
		  al[l2]+=v1/3.;
		}
		dd=fabs(al[0]);
		if(fabs(al[1])>dd) dd=fabs(al[1]);
		if(fabs(al[2])>dd) dd=fabs(al[2]);
		if(dd>worstpsmax){
		  worstpsmax=dd;
		  stnmax[7*l1]=dd;
		  for(l2=1;l2<7;l2++){
		    stnmax[7*l1+l2]=stnl[l2-1];
		  }
		}
	      }
			  
	    }while(1);
	  }
	}
      }
	  
      /* determine the worst principal strain anywhere
	 in the structure as a function of time; 
	 the worst principal strain is the maximum
	 of the absolute value of the principal strains,
	 times its original sign */
	  
      if(strcmp1(&filab[2523],"MAXE")==0){

	/* determining the set of nodes for the 
	   worst principal strain calculation */

	ielset=0;
	for(i=0;i<*nset;i++){
	  if(strcmp1(&set[81*i],"STRAINDOMAINN")==0){
	    ielset=i+1;
	    break;
	  }
	}
	if(ielset==0){
	  printf("\n *ERROR in arpackcs: no node set for MAXE\n");
	  printf("       (must have the name STRAINDOMAIN)\n\n");
	  FORTRAN(stop,());
	}

	for(i1=istartset[ielset-1]-1;i1<iendset[ielset-1];i1++){
	  if(ialset[i1]>0){
	    l1=ialset[i1]-1;

	    worstpsmax=0.;
	    for(l3=0;l3<360;l3++){
	      ctl=cos(l3/constant);
	      stl=sin(l3/constant);
	      for(l2=0;l2<6;l2++){
		l=6*l1+l2;
		eenl[l2]=ctl*een[l]-stl*een[l+6**nk];
	      }
	    
	      /* determining the eigenvalues */
	    
	      v1=eenl[0]+eenl[1]+eenl[2];
	      v2=eenl[1]*eenl[2]+eenl[0]*eenl[2]+eenl[0]*eenl[1]-
		(eenl[5]*eenl[5]+eenl[4]*eenl[4]+eenl[3]*eenl[3]);
	      v3=eenl[0]*(eenl[1]*eenl[2]-eenl[5]*eenl[5])
		-eenl[3]*(eenl[3]*eenl[2]-eenl[4]*eenl[5])
		+eenl[4]*(eenl[3]*eenl[5]-eenl[4]*eenl[1]);
	      bb=v2-v1*v1/3.;
	      cc=-2.*v1*v1*v1/27.+v1*v2/3.-v3;
	      if(fabs(bb)<=1.e-10){
		if(fabs(cc)>1.e-10){
		  al[0]=-pow(cc,(1./3.));
		}else{
		  al[0]=0.;
		}
		al[1]=al[0];
		al[2]=al[0];
	      }else{
		cm=2.*sqrt(-bb/3.);
		cn=3.*cc/(cm*bb);
		if(fabs(cn)>1.){
		  if(cn>1.){
		    cn=1.;
		  }else{
		    cn=-1.;
		  }
		}
		tt=(atan2(sqrt(1.-cn*cn),cn))/3.;
		al[0]=cm*cos(tt);
		al[1]=cm*cos(tt+2.*pi/3.);
		al[2]=cm*cos(tt+4.*pi/3.);
	      }
	      for(l2=0;l2<3;l2++){
		al[l2]+=v1/3.;
	      }
	      dd=fabs(al[0]);
	      if(fabs(al[1])>dd) dd=fabs(al[1]);
	      if(fabs(al[2])>dd) dd=fabs(al[2]);
	      if(dd>worstpsmax){
		worstpsmax=dd;
		eenmax[7*l1]=dd;
		for(l2=1;l2<7;l2++){
		  eenmax[7*l1+l2]=eenl[l2-1];
		}
	      }
	    }
	  
	  }else{
	    l1=ialset[i1-2]-1;
	    do{
	      l1=l1-ialset[i1];
	      if(l1>=ialset[i1-1]-1) break;

	      worstpsmax=0.;
	      for(l3=0;l3<360;l3++){
		ctl=cos(l3/constant);
		stl=sin(l3/constant);
		for(l2=0;l2<6;l2++){
		  l=6*l1+l2;
		  eenl[l2]=ctl*een[l]-stl*een[l+6**nk];
		}
	      
		/* determining the eigenvalues */
	      
		v1=eenl[0]+eenl[1]+eenl[2];
		v2=eenl[1]*eenl[2]+eenl[0]*eenl[2]+eenl[0]*eenl[1]-
		  (eenl[5]*eenl[5]+eenl[4]*eenl[4]+eenl[3]*eenl[3]);
		v3=eenl[0]*(eenl[1]*eenl[2]-eenl[5]*eenl[5])
		  -eenl[3]*(eenl[3]*eenl[2]-eenl[4]*eenl[5])
		  +eenl[4]*(eenl[3]*eenl[5]-eenl[4]*eenl[1]);
		bb=v2-v1*v1/3.;
		cc=-2.*v1*v1*v1/27.+v1*v2/3.-v3;
		if(fabs(bb)<=1.e-10){
		  if(fabs(cc)>1.e-10){
		    al[0]=-pow(cc,(1./3.));
		  }else{
		    al[0]=0.;
		  }
		  al[1]=al[0];
		  al[2]=al[0];
		}else{
		  cm=2.*sqrt(-bb/3.);
		  cn=3.*cc/(cm*bb);
		  if(fabs(cn)>1.){
		    if(cn>1.){
		      cn=1.;
		    }else{
		      cn=-1.;
		    }
		  }
		  tt=(atan2(sqrt(1.-cn*cn),cn))/3.;
		  al[0]=cm*cos(tt);
		  al[1]=cm*cos(tt+2.*pi/3.);
		  al[2]=cm*cos(tt+4.*pi/3.);
		}
		for(l2=0;l2<3;l2++){
		  al[l2]+=v1/3.;
		}
		dd=fabs(al[0]);
		if(fabs(al[1])>dd) dd=fabs(al[1]);
		if(fabs(al[2])>dd) dd=fabs(al[2]);
		if(dd>worstpsmax){
		  worstpsmax=dd;
		  eenmax[7*l1]=dd;
		  for(l2=1;l2<7;l2++){
		    eenmax[7*l1+l2]=eenl[l2-1];
		  }
		}
	      }
	    
	    }while(1);
	  }
	}
      }
    
      /* mapping the results to the other sectors */

      for(l=0;l<*nk;l++){inumt[l]=inum[l];}

      icntrl=2;imag=1;

      FORTRAN(rectcyl,(co,v,fn,stn,qfn,een,cs,nk,&icntrl,t,filab,&imag,mi,emn));

      if((strcmp1(&filab[0],"U  ")==0)||(strcmp1(&filab[870],"PU  ")==0)){
	for(l=0;l<mt**nk;l++){vt[l]=v[l];}
	for(l=0;l<mt**nk;l++){vt[l+mt**nk*ngraph]=v[l+mt**nk];}}
      if(strcmp1(&filab[87],"NT  ")==0){
	if(*iperturb==0){
	  for(l=0;l<*nk;l++){t1t[l]=t1[l];};
	}else{
	  for(l=0;l<*nk;l++){t1t[l]=t1old[l];};
	}
      }
      if((strcmp1(&filab[174],"S   ")==0)||(strcmp1(&filab[1479],"PHS ")==0)){
	for(l=0;l<6**nk;l++){stnt[l]=stn[l];}
	for(l=0;l<6**nk;l++){stnt[l+6**nk*ngraph]=stn[l+6**nk];}}
      if(strcmp1(&filab[261],"E   ")==0){
	for(l=0;l<6**nk;l++){eent[l]=een[l];};
	for(l=0;l<6**nk;l++){eent[l+6**nk*ngraph]=een[l+6**nk];}}
      if((strcmp1(&filab[348],"RF  ")==0)||(strcmp1(&filab[2610],"PRF ")==0)){
	for(l=0;l<mt**nk;l++){fnt[l]=fn[l];}
	for(l=0;l<mt**nk;l++){fnt[l+mt**nk*ngraph]=fn[l+mt**nk];}}
      if(strcmp1(&filab[522],"ENER")==0){
	for(l=0;l<*nk;l++){enernt[l]=enern[l];}}
      //      for(l=0;l<*nk;l++){enernt[l+*nk*ngraph]=enern[l+*nk];}}
      if((strcmp1(&filab[1044],"ZZS ")==0)||(strcmp1(&filab[1044],"ERR ")==0)||
	 ((strcmp1(&filab[2175],"CONT")==0)&&(*mortar==0))){
	for(l=0;l<6*mi[0]**ne;l++){stxt[l]=stx[l];}
	for(l=0;l<6*mi[0]**ne;l++){stxt[l+6*mi[0]**ne*ngraph]=stx[l+6*mi[0]**ne];}}
      if(strcmp1(&filab[2697],"ME  ")==0){
	for(l=0;l<6**nk;l++){emnt[l]=emn[l];};
	for(l=0;l<6**nk;l++){emnt[l+6**nk*ngraph]=emn[l+6**nk];}}
      if(((strcmp1(&filab[2175],"CONT")==0)||(strcmp1(&filab[3915],"PCON")==0))
	 &&(*mortar==1)){
	for(l=0;l<6**nk;l++){cdnt[l]=cdn[l];}
	for(l=0;l<6**nk;l++){cdnt[l+6**nk*ngraph]=cdn[l+6**nk];}}

      for(jj=0;jj<*mcs;jj++){
	ilength=cs[17*jj+3];
	is=cs[17*jj+4];
	lprev=cs[17*jj+13];
	for(i=1;i<is;i++){
        
	  for(l=0;l<*nk;l++){inumt[l+i**nk]=inum[l];}
        
	  theta=i*nm*2.*pi/cs[17*jj];
	  ctl=cos(theta);
	  stl=sin(theta);
        
	  if((strcmp1(&filab[0],"U  ")==0)||(strcmp1(&filab[870],"PU  ")==0)){
	    for(l1=0;l1<*nk;l1++){
	      if(inocs[l1]==jj){

		/* check whether node lies on axis */

		ml1=-l1-1;
		FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		    for(l2=0;l2<4;l2++){
		      l=mt*l1+l2;
		      vt[l+mt**nk*i]=v[l];
		    }
		    continue;
		  }
		}
		for(l2=0;l2<4;l2++){
		  l=mt*l1+l2;
		  vt[l+mt**nk*i]=ctl*v[l]-stl*v[l+mt**nk];
		}
	      }
	    }
	  }
        
	  /* imaginary part of the displacements in cylindrical
	     coordinates */

	  if((strcmp1(&filab[0],"U  ")==0)||(strcmp1(&filab[870],"PU  ")==0)){
	    for(l1=0;l1<*nk;l1++){
	      if(inocs[l1]==jj){

		/* check whether node lies on axis */

		ml1=-l1-1;
		FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		    for(l2=0;l2<4;l2++){
		      l=mt*l1+l2;
		      vt[l+mt**nk*(i+ngraph)]=v[l+mt**nk];
		    }
		    continue;
		  }
		}
		for(l2=0;l2<4;l2++){
		  l=mt*l1+l2;
		  vt[l+mt**nk*(i+ngraph)]=stl*v[l]+ctl*v[l+mt**nk];
		}
	      }
	    }
	  }
        
	  if(strcmp1(&filab[87],"NT  ")==0){
	    if(*iperturb==0){
	      for(l=0;l<*nk;l++){
		if(inocs[l]==jj) t1t[l+*nk*i]=t1[l];
	      }
	    }else{
	      for(l=0;l<*nk;l++){
		if(inocs[l]==jj) t1t[l+*nk*i]=t1old[l];
	      }
	    }
	  }
        
	  if((strcmp1(&filab[174],"S   ")==0)||(strcmp1(&filab[1479],"PHS ")==0)){
	    for(l1=0;l1<*nk;l1++){
	      if(inocs[l1]==jj){

		/* check whether node lies on axis */

		ml1=-l1-1;
		FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		    for(l2=0;l2<6;l2++){
		      l=6*l1+l2;
		      stnt[l+6**nk*i]=stn[l];
		    }
		    continue;
		  }
		}
		for(l2=0;l2<6;l2++){
		  l=6*l1+l2;
		  stnt[l+6**nk*i]=ctl*stn[l]-stl*stn[l+6**nk];
		}
	      }
	    }
	  }
        
	  /* imaginary part of the stresses in cylindrical
	     coordinates */
        
	  if((strcmp1(&filab[174],"S   ")==0)||(strcmp1(&filab[1479],"PHS ")==0)){
	    for(l1=0;l1<*nk;l1++){
	      if(inocs[l1]==jj){

		/* check whether node lies on axis */

		ml1=-l1-1;
		FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		    for(l2=0;l2<6;l2++){
		      l=6*l1+l2;
		      stnt[l+6**nk*(i+ngraph)]=stn[l+6**nk];
		    }
		    continue;
		  }
		}
		for(l2=0;l2<6;l2++){
		  l=6*l1+l2;
		  stnt[l+6**nk*(i+ngraph)]=stl*stn[l]+ctl*stn[l+6**nk];
		}
	      }
	    }
	  }
        
	  if(strcmp1(&filab[261],"E   ")==0){
	    for(l1=0;l1<*nk;l1++){
	      if(inocs[l1]==jj){

		/* check whether node lies on axis */

		ml1=-l1-1;
		FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		    for(l2=0;l2<6;l2++){
		      l=6*l1+l2;
		      eent[l+6**nk*i]=een[l];
		    }
		    continue;
		  }
		}
		for(l2=0;l2<6;l2++){
		  l=6*l1+l2;
		  eent[l+6**nk*i]=ctl*een[l]-stl*een[l+6**nk];
		}
	      }
	    }
	  }
        
	  /* imaginary part of the strains in cylindrical
	     coordinates */
        
	  if(strcmp1(&filab[261],"E   ")==0){
	    for(l1=0;l1<*nk;l1++){
	      if(inocs[l1]==jj){

		/* check whether node lies on axis */

		ml1=-l1-1;
		FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		    for(l2=0;l2<6;l2++){
		      l=6*l1+l2;
		      eent[l+6**nk*(i+ngraph)]=een[l+6**nk];
		    }
		    continue;
		  }
		}
		for(l2=0;l2<6;l2++){
		  l=6*l1+l2;
		  eent[l+6**nk*(i+ngraph)]=stl*een[l]+ctl*een[l+6**nk];
		}
	      }
	    }
	  }
        
	  /* real part of the mechanical strains */

	  if(strcmp1(&filab[2697],"ME  ")==0){
	    for(l1=0;l1<*nk;l1++){
	      if(inocs[l1]==jj){

		/* check whether node lies on axis */

		ml1=-l1-1;
		FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		    for(l2=0;l2<6;l2++){
		      l=6*l1+l2;
		      emnt[l+6**nk*i]=emn[l];
		    }
		    continue;
		  }
		}
		for(l2=0;l2<6;l2++){
		  l=6*l1+l2;
		  emnt[l+6**nk*i]=ctl*emn[l]-stl*emn[l+6**nk];
		}
	      }
	    }
	  }
        
	  /* imaginary part of the mechanical strains in cylindrical
	     coordinates */
        
	  if(strcmp1(&filab[2697],"ME  ")==0){
	    for(l1=0;l1<*nk;l1++){
	      if(inocs[l1]==jj){

		/* check whether node lies on axis */

		ml1=-l1-1;
		FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		    for(l2=0;l2<6;l2++){
		      l=6*l1+l2;
		      emnt[l+6**nk*(i+ngraph)]=emn[l+6**nk];
		    }
		    continue;
		  }
		}
		for(l2=0;l2<6;l2++){
		  l=6*l1+l2;
		  emnt[l+6**nk*(i+ngraph)]=stl*emn[l]+ctl*emn[l+6**nk];
		}
	      }
	    }
	  }
        
	  /* real part of the contact stresses */

	  if(((strcmp1(&filab[2175],"CONT")==0)||(strcmp1(&filab[3915],"PCON")==0))
	     &&(*mortar==1)){
	    for(l1=0;l1<*nk;l1++){
	      if(inocs[l1]==jj){

		/* check whether node lies on axis */

		ml1=-l1-1;
		FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		    for(l2=0;l2<6;l2++){
		      l=6*l1+l2;
		      cdnt[l+6**nk*i]=cdn[l];
		    }
		    continue;
		  }
		}
		for(l2=0;l2<6;l2++){
		  l=6*l1+l2;
		  cdnt[l+6**nk*i]=ctl*cdn[l]-stl*cdn[l+6**nk];
		}
	      }
	    }
	  }
        
	  /* imaginary part of the contact stresses */
        
	  if(((strcmp1(&filab[2175],"CONT")==0)||(strcmp1(&filab[3915],"PCON")==0))
	     &&(*mortar==1)){
	    for(l1=0;l1<*nk;l1++){
	      if(inocs[l1]==jj){

		/* check whether node lies on axis */

		ml1=-l1-1;
		FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		    for(l2=0;l2<6;l2++){
		      l=6*l1+l2;
		      cdnt[l+6**nk*(i+ngraph)]=cdn[l+6**nk];
		    }
		    continue;
		  }
		}
		for(l2=0;l2<6;l2++){
		  l=6*l1+l2;
		  cdnt[l+6**nk*(i+ngraph)]=stl*cdn[l]+ctl*cdn[l+6**nk];
		}
	      }
	    }
	  }
        
	  if((strcmp1(&filab[348],"RF  ")==0)||(strcmp1(&filab[2610],"PRF ")==0)){
	    for(l1=0;l1<*nk;l1++){
	      if(inocs[l1]==jj){

		/* check whether node lies on axis */

		ml1=-l1-1;
		FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		    for(l2=0;l2<4;l2++){
		      l=mt*l1+l2;
		      fnt[l+mt**nk*i]=fn[l];
		    }
		    continue;
		  }
		}
		for(l2=0;l2<4;l2++){
		  l=mt*l1+l2;
		  fnt[l+mt**nk*i]=ctl*fn[l]-stl*fn[l+mt**nk];
		}
	      }
	    }
	  }
        
	  /* imaginary part of the forces in cylindrical
	     coordinates */

	  if((strcmp1(&filab[348],"RF  ")==0)||(strcmp1(&filab[2610],"PRF ")==0)){
	    for(l1=0;l1<*nk;l1++){
	      if(inocs[l1]==jj){
	      
		/* check whether node lies on axis */
	      
		ml1=-l1-1;
		FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
		if(id!=0){
		  if(ics[lprev+id-1]==ml1){
		    for(l2=0;l2<4;l2++){
		      l=mt*l1+l2;
		      fnt[l+mt**nk*(i+ngraph)]=fn[l+mt**nk];
		    }
		    continue;
		  }
		}
		for(l2=0;l2<4;l2++){
		  l=mt*l1+l2;
		  fnt[l+mt**nk*(i+ngraph)]=stl*fn[l]+ctl*fn[l+mt**nk];
		}
	      }
	    }
	  }

	}
      }
	
      icntrl=-2;imag=0;

      FORTRAN(rectcyl,(cot,vt,fnt,stnt,qfnt,eent,cs,&nkt,&icntrl,t,filab,
		       &imag,mi,emnt));

      FORTRAN(rectcylvi,(cot,&vt[mt**nk*ngraph],&fnt[mt**nk*ngraph],
			 &stnt[6**nk*ngraph],qfnt,&eent[6**nk*ngraph],cs,&nkt,&icntrl,
			 t,filab,&imag,mi,&emnt[6**nk*ngraph]));
        
      /* internal energy calculation */

      for(jj=0;jj<*mcs;jj++){
	ilength=cs[17*jj+3];
	is=cs[17*jj+4];
	lprev=cs[17*jj+13];
	for(i=1;i<is;i++){
        
	  for(l=0;l<*nk;l++){inumt[l+i**nk]=inum[l];}
        
	  theta=i*nm*2.*pi/cs[17*jj];
	  ctl=cos(theta);
	  stl=sin(theta);

	  if(strcmp1(&filab[522],"ENER")==0){
	    ioffr=6**nk*i;
	    for(l1=0;l1<*nk;l1++){
	      ioffrl=6*l1+ioffr;
	      enernt[l1+*nk*i]=(stnt[ioffrl]*emnt[ioffrl]+
	                        stnt[1+ioffrl]*emnt[1+ioffrl]+
                                stnt[2+ioffrl]*emnt[2+ioffrl])/2.+
		stnt[3+ioffrl]*emnt[3+ioffrl]+
		stnt[4+ioffrl]*emnt[4+ioffrl]+
		stnt[5+ioffrl]*emnt[5+ioffrl];
	    }
	  }

	}
      }

      /* determining magnitude and phase angle for the displacements */

      if(strcmp1(&filab[870],"PU")==0){
	for(l1=0;l1<nkt;l1++){
	  for(l2=0;l2<4;l2++){
	    l=mt*l1+l2;
	    vreal=vt[l];
	    vimag=vt[l+mt**nk*ngraph];
	    vr[l]=sqrt(vreal*vreal+vimag*vimag);
	    if(vr[l]<1.e-10){
	      vi[l]=0.;
	    }else if(fabs(vreal)<1.e-10){
	      if(vimag>0){vi[l]=90.;}
	      else{vi[l]=-90.;}
	    }
	    else{
	      vi[l]=atan(vimag/vreal)*constant;
	      if(vreal<0) vi[l]+=180.;
	    }
	  }
	}
      }

      /* determining magnitude and phase for the stress */

      if(strcmp1(&filab[1479],"PHS")==0){
	for(l1=0;l1<nkt;l1++){
	  for(l2=0;l2<6;l2++){
	    l=6*l1+l2;
	    stnreal=stnt[l];
	    stnimag=stnt[l+6**nk*ngraph];
	    stnr[l]=sqrt(stnreal*stnreal+stnimag*stnimag);
	    if(stnr[l]<1.e-10){
	      stni[l]=0.;
	    }else if(fabs(stnreal)<1.e-10){
	      if(stnimag>0){stni[l]=90.;}
	      else{stni[l]=-90.;}
	    }
	    else{
	      stni[l]=atan(stnimag/stnreal)*constant;
	      if(stnreal<0) stni[l]+=180.;
	    }
	  }
	}
      }

      /* determining magnitude and phase angle for the forces */

      if(strcmp1(&filab[2610],"PRF")==0){
	for(l1=0;l1<nkt;l1++){
	  for(l2=0;l2<4;l2++){
	    l=mt*l1+l2;
	    fnreal=fnt[l];
	    fnimag=fnt[l+mt**nk*ngraph];
	    fnr[l]=sqrt(fnreal*fnreal+fnimag*fnimag);
	    if(fnr[l]<1.e-10){
	      fni[l]=0.;
	    }else if(fabs(fnreal)<1.e-10){
	      if(fnimag>0){fni[l]=90.;}
	      else{fni[l]=-90.;}
	    }
	    else{
	      fni[l]=atan(fnimag/fnreal)*constant;
	      if(fnreal<0) fni[l]+=180.;
	    }
	  }
	}
      }

      /* determining magnitude and phase for the contact stress */

      if((strcmp1(&filab[3915],"PCON")==0)&&(*mortar==1)){
	for(l1=0;l1<nkt;l1++){
	  for(l2=0;l2<6;l2++){
	    l=6*l1+l2;
	    cdnreal=cdnt[l];
	    cdnimag=cdnt[l+6**nk*ngraph];
	    cdnr[l]=sqrt(cdnreal*cdnreal+cdnimag*cdnimag);
	    if(cdnr[l]<1.e-10){
	      cdni[l]=0.;
	    }else if(fabs(cdnreal)<1.e-10){
	      if(cdnimag>0){cdni[l]=90.;}
	      else{cdni[l]=-90.;}
	    }
	    else{
	      cdni[l]=atan(cdnimag/cdnreal)*constant;
	      if(cdnreal<0) cdni[l]+=180.;
	    }
	  }
	}
      }

      ++*kode;
      if(d[j]>=0.){
	freq=sqrt(d[j])/6.283185308;
      }else{
	freq=0.;
      }

      if(strcmp1(&filab[1044],"ZZS")==0){
	NNEW(neigh,ITG,40*net);
	NNEW(ipneigh,ITG,nkt);
      }

      /* deactivating S and ME request if only needed for ENER */

      if((ngraph>1)&&(strcmp1(&filab[522],"ENER")==0)){
	strcpy1(&filab[174],&filabcp[0],4);
	strcpy1(&filab[2697],&filabcp[4],4);
      }

      if(igreen==1) *nmethod=13;

      /* change on 20210510: if the nodal diameter exeeds half the
         number of segments change the sign of the axis */
      
      //      if(nm>cs[0]/2){
      if(nm-floor(nm/cs[0])*cs[0]>cs[0]/2){
	for(k=5;k<11;k++){cs[k]=-cs[k];}
      }
      frd(cot,&nkt,kont,ipkont,lakont,&net,vt,stnt,inumt,nmethod,
	  kode,filab,eent,t1t,fnt,&freq,epn,ielmatt,matname,enernt,xstaten,
	  nstate_,istep,&iinc,ithermal,qfn,&j,&nm,trab,inotrt,
	  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	  mi,stxt,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,&net,
	  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emnt,
	  thicke,jobnamec,output,qfx,cdnt,mortar,cdnr,cdni,nmat,
	  ielprop,prop,sti);
      if(nm>cs[0]/2){
	for(k=5;k<11;k++){cs[k]=-cs[k];}
      }

      if(igreen==1) *nmethod=2;

      /* reactivating S and ME request if only needed for ENER */

      if((ngraph>1)&&(strcmp1(&filab[522],"ENER")==0)){
	strcpy1(&filabcp[0],&filab[174],4);
	strcpy1(&filabcp[4],&filab[2697],4);
	strcpy1(&filab[174],"S   ",4);
	strcpy1(&filab[2697],"ME  ",4);
      }

      if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);}

      if(*iperturb!=0){
	for(k=0;k<*nstate_*mi[0]*(ne0+maxprevcontel);++k){
	  xstate[k]=xstateini[k];
	}	  
      }

      /* end loop over the eigenvalues */

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

    /* output of the turning direction */
      
    FORTRAN(writeturdircs,(xn,turdir,&nev,&nm));
    SFREE(turdir);

    if(igreen==0){
      if((fmax>-0.5)&&(fmax*fmax>d[nev-1])){
	printf("\n*WARNING: not all frequencies in the requested interval might be found;\nincrease the number of requested frequencies\n");
      }
    }

    SFREE(adb);SFREE(aub);SFREE(coefmpcnew);SFREE(ad);SFREE(au);

    /* deactivating S and ME request if only needed for ENER */

    if((ngraph>1)&&(strcmp1(&filab[522],"ENER")==0)){
      strcpy1(&filab[174],&filabcp[0],4);
      strcpy1(&filab[2697],&filabcp[4],4);
      if(strcmp1(&filab[174],"S   ")!=0){SFREE(stn);SFREE(stnt);}
      if(strcmp1(&filab[2697],"ME  ")!=0){SFREE(emn);SFREE(emnt);}
    }

    if((strcmp1(&filab[174],"S   ")==0)||(strcmp1(&filab[1653],"MAXS")==0)|| 
       (strcmp1(&filab[1479],"PHS ")==0)||(strcmp1(&filab[1044],"ZZS ")==0)||
       (strcmp1(&filab[1044],"ERR ")==0)) 
      SFREE(stn);

    SFREE(v);SFREE(fn);SFREE(inum);SFREE(stx);SFREE(eme);SFREE(z);SFREE(d);

    if((strcmp1(&filab[261],"E   ")==0)||(strcmp1(&filab[2523],"MAXE")==0)) SFREE(een);
    if(strcmp1(&filab[522],"ENER")==0) SFREE(enern);
    if(strcmp1(&filab[2697],"ME  ")==0) SFREE(emn);

    if((strcmp1(&filab[0],"U  ")==0)||(strcmp1(&filab[870],"PU  ")==0)) SFREE(vt);
    if(strcmp1(&filab[87],"NT  ")==0) SFREE(t1t);
    if((strcmp1(&filab[174],"S   ")==0)||(strcmp1(&filab[1479],"PHS ")==0)||
       (strcmp1(&filab[1044],"ZZS ")==0)||(strcmp1(&filab[1044],"ERR ")==0)) SFREE(stnt);
    if(strcmp1(&filab[261],"E   ")==0) SFREE(eent);
    if((strcmp1(&filab[348],"RF  ")==0)||(strcmp1(&filab[2610],"PRF ")==0)) SFREE(fnt);
    if(strcmp1(&filab[522],"ENER")==0) SFREE(enernt);
    if((strcmp1(&filab[1044],"ZZS ")==0)||(strcmp1(&filab[1044],"ERR ")==0)||
       ((strcmp1(&filab[2175],"CONT")==0)&&(*mortar==0))) SFREE(stxt);
    if(strcmp1(&filab[2697],"ME  ")==0) SFREE(emnt);
    if(((strcmp1(&filab[2175],"CONT")==0)||(strcmp1(&filab[3915],"PCON")==0))
       &&(*mortar==1)) SFREE(cdnt);

    SFREE(cot);SFREE(kont);SFREE(ipkont);SFREE(lakont);SFREE(inumt);SFREE(ielmatt);
    if(*ntrans>0){SFREE(inotrt);}

    if(mei[3]==1){
      (*nevtot)+=nev;
      fclose(f1);
    }

    /* end loop over the nodal diameters */

  }

  if(*iperturb!=0){
    if(ncont!=0){
      *ne=ne0;*nkon=nkon0;
      //	  if(*nener==1){SFREE(ener);}
      RENEW(ipkon,ITG,*ne);
      RENEW(lakon,char,8**ne);
      RENEW(kon,ITG,*nkon);
      if(*norien>0){
	RENEW(ielorien,ITG,mi[2]**ne);
      }
      RENEW(ielmat,ITG,mi[2]**ne);
      SFREE(cg);
      SFREE(straight);
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

  if(*nener==1){SFREE(ener);}
  SFREE(inocs);SFREE(ielcs);SFREE(xstiff);
  
  //  if(*nbody>0) SFREE(ipobody);

  if(*nstate_!=0) SFREE(xstateini);

  if(strcmp1(&filab[870],"PU")==0){SFREE(vr);SFREE(vi);}
  if(strcmp1(&filab[1479],"PHS")==0){SFREE(stnr);SFREE(stni);}
  if(strcmp1(&filab[3915],"PCON")==0){SFREE(cdnr);SFREE(cdni);}
  if(strcmp1(&filab[1566],"MAXU")==0){SFREE(vmax);}
  if(strcmp1(&filab[1653],"MAXS")==0){SFREE(stnmax);}
  if(strcmp1(&filab[2523],"MAXE")==0){SFREE(eenmax);}
  if(strcmp1(&filab[2610],"PRF")==0){SFREE(fnr);SFREE(fni);}

  if(*iperturb!=0){
    mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
    mpcinfo[3]=maxlenmpc;
  }

  *irowp=irow;*xstatep=xstate;*ipkonp=ipkon;*lakonp=lakon;
  *konp=kon;*ielmatp=ielmat;*ielorienp=ielorien;

  *islavsurfp=islavsurf;*pslavsurfp=pslavsurf;*clearinip=clearini;

  return;
}


#endif
