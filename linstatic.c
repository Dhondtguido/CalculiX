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
#ifdef PARDISO
#include "pardiso.h"
#endif
#ifdef PASTIX
#include "pastix.h"
#endif

void linstatic(double *co,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,
	       ITG *ne,
	       ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
	       ITG *ipompc,ITG *nodempc,double *coefmpc,char *labmpc,
	       ITG *nmpc,
	       ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
	       ITG *nelemload,char *sideload,double *xload,
	       ITG *nload,ITG *nactdof,
	       ITG **icolp,ITG *jq,ITG **irowp,ITG *neq,ITG *nzl,
	       ITG *nmethod,ITG *ikmpc,ITG *ilmpc,ITG *ikboun,
	       ITG *ilboun,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG **ielmatp,
	       ITG **ielorienp,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,double *t1old,
	       ITG *ithermal,double *prestr,ITG *iprestr,
	       double *vold,ITG *iperturb,double *sti,ITG *nzs, 
	       ITG *kode,char *filab,double *eme,
	       ITG *iexpl,double *plicon,ITG *nplicon,double *plkcon,
	       ITG *nplkcon,
	       double **xstatep,ITG *npmat_,char *matname,ITG *isolver,
	       ITG *mi,ITG *ncmat_,ITG *nstate_,double *cs,ITG *mcs,
	       ITG *nkon,double **enerp,double *xbounold,
	       double *xforcold,double *xloadold,
	       char *amname,double *amta,ITG *namta,
	       ITG *nam,ITG *iamforc,ITG *iamload,
	       ITG *iamt1,ITG *iamboun,double *ttime,char *output,
	       char *set,ITG *nset,ITG *istartset,
	       ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
	       char *prset,ITG *nener,double *trab,
	       ITG *inotr,ITG *ntrans,double *fmpc,ITG *ipobody,ITG *ibody,
	       double *xbody,ITG *nbody,double *xbodyold,double *timepar,
	       double *thicke,char *jobnamec,char *tieset,ITG *ntie,
	       ITG *istep,ITG *nmat,ITG *ielprop,double *prop,char *typeboun,
	       ITG *mortar,ITG *mpcinfo,double *tietol,ITG *ics,
	       char *orname,ITG *itempuser,double *t0g,double *t1g){
  
  char description[13]="            ",*lakon=NULL,stiffmatrix[132]="",
    fneig[132]="",jobnamef[396]="",*labmpc2=NULL;

  ITG *inum=NULL,k,*icol=NULL,*irow=NULL,ielas=0,icmd=0,iinc=1,nasym=0,i,j,ic,ir,
    mass[2]={0,0},stiffness=1,buckling=0,rhsi=1,intscheme=0,*ncocon=NULL,
    *nshcon=NULL,mode=-1,noddiam=-1,coriolis=0,iout,
    *itg=NULL,ntg=0,symmetryflag=0,inputformat=0,ngraph=1,im,
    mt=mi[1]+1,ne0,*integerglob=NULL,iglob=0,*ipneigh=NULL,*neigh=NULL,
    icfd=0,*inomat=NULL,*islavact=NULL,*islavnode=NULL,*nslavnode=NULL,
    *islavsurf=NULL,nretain,*iretain=NULL,*noderetain=NULL,*ndirretain=NULL,
    nmethodl,nintpoint,ifacecount,memmpc_,mpcfree,icascade,maxlenmpc,
    ncont=0,*itietri=NULL,*koncont=NULL,nslavs=0,ismallsliding=0,
    *itiefac=NULL,*imastnode=NULL,*nmastnode=NULL,*imastop=NULL,iitsta,
    *iponoels=NULL,*inoels=NULL,*ipe=NULL,*ime=NULL,iit=-1,iflagact=0,
    icutb=0,*kon=NULL,*ipkon=NULL,*ielmat=NULL,ialeatoric=0,kscale=1,
    *iponoel=NULL,*inoel=NULL,zero=0,nherm=1,nev=*nforc,node,idir,
    *ielorien=NULL,network=0,nrhs=1,iperturbsav,mscalmethod=0,*jqw=NULL,
    *iroww=NULL,nzsw,*islavelinv=NULL,*irowtloc=NULL,*jqtloc=NULL,nboun2,
    *ndirboun2=NULL,*nodeboun2=NULL,nmpc2,*ipompc2=NULL,*nodempc2=NULL,
    *ikboun2=NULL,*ilboun2=NULL,*ikmpc2=NULL,*ilmpc2=NULL,mortartrafoflag=0;

  double *stn=NULL,*v=NULL,*een=NULL,cam[5],*xstiff=NULL,*stiini=NULL,*tper,
    *f=NULL,*fn=NULL,qa[4],*fext=NULL,*epn=NULL,*xstateini=NULL,
    *vini=NULL,*stx=NULL,*enern=NULL,*xbounact=NULL,*xforcact=NULL,
    *xloadact=NULL,*t1act=NULL,*ampli=NULL,*xstaten=NULL,*eei=NULL,
    *enerini=NULL,*cocon=NULL,*shcon=NULL,*physcon=NULL,*qfx=NULL,
    *qfn=NULL,sigma=0.,*cgr=NULL,*xbodyact=NULL,*vr=NULL,*vi=NULL,
    *stnr=NULL,*stni=NULL,*vmax=NULL,*stnmax=NULL,*springarea=NULL,
    *eenmax=NULL,*fnr=NULL,*fni=NULL,*emn=NULL,*clearini=NULL,ptime,
    *emeini=NULL,*doubleglob=NULL,*au=NULL,*ad=NULL,*b=NULL,*aub=NULL,
    *adb=NULL,*pslavsurf=NULL,*pmastsurf=NULL,*cdn=NULL,*cdnr=NULL,
    *cdni=NULL,*submatrix=NULL,*xnoels=NULL,*cg=NULL,*straight=NULL,
    *areaslav=NULL,*xmastnor=NULL,theta=0.,*ener=NULL,*xstate=NULL,
    *fnext=NULL,*energyini=NULL,*energy=NULL,*d=NULL,alea=0.1,*smscale=NULL,
    *auw=NULL,*autloc=NULL,*xboun2=NULL,*coefmpc2=NULL;

  FILE *f1,*f2;
  
#ifdef SGI
  ITG token;
#endif

  /* dummy arguments for the results call */

  double *veold=NULL,*accold=NULL,bet,gam,dtime,time,reltime=1.;

  irow=*irowp;ener=*enerp;xstate=*xstatep;ipkon=*ipkonp;lakon=*lakonp;
  kon=*konp;ielmat=*ielmatp;ielorien=*ielorienp;icol=*icolp;
  
  for(k=0;k<3;k++){
    strcpy1(&jobnamef[k*132],&jobnamec[k*132],132);
  }

  tper=&timepar[1];

  time=*tper;
  dtime=*tper;

  ne0=*ne;

  /* determining the global values to be used as boundary conditions
     for a submodel */

  /* iglob=-1 if global results are from a *FREQUENCY calculation
     iglob=0 if no global results are used by boundary conditions
     iglob=1 if global results are from a *STATIC calculation */

  ITG irefine=0;
  getglobalresults(&jobnamec[396],&integerglob,&doubleglob,nboun,iamboun,xboun,
		   nload,sideload,iamload,&iglob,nforc,iamforc,xforc,
                   ithermal,nk,t1,iamt1,&sigma,&irefine);

  /* reading temperatures from frd-file */
  
  if((itempuser[0]==2)&&(itempuser[1]!=itempuser[2])) {
    utempread(t1,&itempuser[2],jobnamec);
  }      

  /* allocating fields for the actual external loading */

  NNEW(xbounact,double,*nboun);
  for(k=0;k<*nboun;++k){xbounact[k]=xbounold[k];}
  NNEW(xforcact,double,*nforc);
  NNEW(xloadact,double,2**nload);
  NNEW(xbodyact,double,7**nbody);
  /* copying the rotation axis and/or acceleration vector */
  for(k=0;k<7**nbody;k++){xbodyact[k]=xbody[k];}
  if(*ithermal==1){
    NNEW(t1act,double,*nk);
    for(k=0;k<*nk;++k){t1act[k]=t1old[k];}
  }
  
  /* assigning the body forces to the elements */ 

  /*  if(*nbody>0){
    ifreebody=*ne+1;
    NNEW(ipobody,ITG,2*ifreebody**nbody);
    for(k=1;k<=*nbody;k++){
      FORTRAN(bodyforce,(cbody,ibody,ipobody,nbody,set,istartset,
			 iendset,ialset,&inewton,nset,&ifreebody,&k));
      RENEW(ipobody,ITG,2*(*ne+ifreebody));
    }
    RENEW(ipobody,ITG,2*(ifreebody-1));
    }*/

  /* contact conditions */
  
  //  if(*icontact==1){
  if(*mortar>-2){

    memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];icascade=mpcinfo[2];
    maxlenmpc=mpcinfo[3];

    inicont(nk,&ncont,ntie,tieset,nset,set,istartset,iendset,ialset,&itietri,
	    lakon,ipkon,kon,&koncont,&nslavs,tietol,&ismallsliding,&itiefac,
	    &islavsurf,&islavnode,&imastnode,&nslavnode,&nmastnode,
	    mortar,&imastop,nkon,&iponoels,&inoels,&ipe,&ime,ne,&ifacecount,
	    iperturb,ikboun,nboun,co,istep,&xnoels);

    if(ncont!=0){

      NNEW(cg,double,3*ncont);
      NNEW(straight,double,16*ncont);
	  
      /* 11 instead of 10: last position is reserved for the
	 local contact spring element number; needed as
	 pointer into springarea */
	  
      if(*mortar==0){
	RENEW(kon,ITG,*nkon+11*nslavs);
	NNEW(springarea,double,2*nslavs);
	if(*nener==1){
	  RENEW(ener,double,2*mi[0]*(*ne+nslavs));
	  DMEMSET(ener,2*mi[0]**ne,2*mi[0]*(*ne+nslavs),0.);
	}
	RENEW(ipkon,ITG,*ne+nslavs);
	RENEW(lakon,char,8*(*ne+nslavs));
	      
	if(*norien>0){
	  RENEW(ielorien,ITG,mi[2]*(*ne+nslavs));
	  for(k=mi[2]**ne;k<mi[2]*(*ne+nslavs);k++) ielorien[k]=0;
	}
	      
	RENEW(ielmat,ITG,mi[2]*(*ne+nslavs));
	for(k=mi[2]**ne;k<mi[2]*(*ne+nslavs);k++) ielmat[k]=1;
	      
	if(nslavs!=0){
	  RENEW(xstate,double,*nstate_*mi[0]*(*ne+nslavs));
	  for(k=*nstate_*mi[0]**ne;k<*nstate_*mi[0]*(*ne+nslavs);k++){
	    xstate[k]=0.;
	  }
	}
	      
	NNEW(areaslav,double,ifacecount);
	NNEW(xmastnor,double,3*nmastnode[*ntie]);
      }else if(*mortar==1){
	NNEW(islavact,ITG,nslavnode[*ntie]);
	DMEMSET(islavact,0,nslavnode[*ntie],1);
	NNEW(clearini,double,3*9*ifacecount);
	NNEW(xmastnor,double,3*nmastnode[*ntie]);


	nintpoint=0;
	      
	precontact(&ncont,ntie,tieset,nset,set,istartset,
		   iendset,ialset,itietri,lakon,ipkon,kon,koncont,ne,
		   cg,straight,co,vold,istep,&iinc,&iit,itiefac,
		   islavsurf,islavnode,imastnode,nslavnode,nmastnode,
		   imastop,mi,ipe,ime,tietol,&iflagact,
		   &nintpoint,&pslavsurf,xmastnor,cs,mcs,ics,clearini,
		   &nslavs);
	      
	/* changing the dimension of element-related fields */
	      
	RENEW(kon,ITG,*nkon+22*nintpoint);
	RENEW(springarea,double,2*nintpoint);
	RENEW(pmastsurf,double,6*nintpoint);
	      
	if(*nener==1){
	  RENEW(ener,double,2*mi[0]*(*ne+nintpoint));
	  DMEMSET(ener,2*mi[0]**ne,2*mi[0]*(*ne+nintpoint),0.);
	}
	RENEW(ipkon,ITG,*ne+nintpoint);
	RENEW(lakon,char,8*(*ne+nintpoint));
	      
	if(*norien>0){
	  RENEW(ielorien,ITG,mi[2]*(*ne+nintpoint));
	  for(k=mi[2]**ne;k<mi[2]*(*ne+nintpoint);k++) ielorien[k]=0;
	}
	RENEW(ielmat,ITG,mi[2]*(*ne+nintpoint));
	for(k=mi[2]**ne;k<mi[2]*(*ne+nintpoint);k++) ielmat[k]=1;

	/* interpolating the state variables */

	if(*nstate_!=0){
		  
	  RENEW(xstate,double,*nstate_*mi[0]*(ne0+nintpoint));
	  for(k=*nstate_*mi[0]*ne0;k<*nstate_*mi[0]*(ne0+nintpoint);k++){
	    xstate[k]=0.;
	  }
		  
	  RENEW(xstateini,double,*nstate_*mi[0]*(ne0+nintpoint));
	  for(k=0;k<*nstate_*mi[0]*(ne0+nintpoint);++k){
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
		   neq,nzs,nmethod,ithermal,iperturb,mass,mi,ics,cs,
		   mcs,mortar,typeboun,&iit,&network,iexpl);
    }

    /* field for initial values of state variables (needed for contact */

    if((*nstate_!=0)&&((*mortar==0)||(ncont==0))){
      NNEW(xstateini,double,*nstate_*mi[0]*(ne0+nslavs));
      for(k=0;k<*nstate_*mi[0]*(ne0+nslavs);++k){
	xstateini[k]=xstate[k];
      }
    }
  }

  /* allocating a field for the instantaneous amplitude */

  NNEW(ampli,double,*nam);

  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,xload,
		    xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,xbodyact,
		    t1old,t1,t1act,iamt1,nk,amta,
		    namta,nam,ampli,&time,&reltime,ttime,&dtime,ithermal,nmethod,
		    xbounold,xboun,xbounact,iamboun,nboun,
		    nodeboun,ndirboun,nodeforc,ndirforc,istep,&iinc,
		    co,vold,itg,&ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
		    ntrans,trab,inotr,veold,integerglob,doubleglob,tieset,istartset,
		    iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc,
		    ipobody,iponoel,inoel,ipkon,kon,ielprop,prop,ielmat,
		    shcon,nshcon,rhcon,nrhcon,cocon,ncocon,ntmat_,lakon,
		    set,nset));

  /* determining the internal forces and the stiffness coefficients */

  NNEW(f,double,*neq);

  /* allocating a field for the stiffness matrix */

  NNEW(xstiff,double,(long long)27*mi[0]**ne);

  /* for a *STATIC,PERTURBATION analysis with submodel boundary
     conditions from a *FREQUENCY analysis iperturb[0]=1 has to be
     temporarily set to iperturb[0]=0 in order for f to be calculated in
     resultsini and subsequent results* routines */

  if((*nmethod==1)&&(iglob<0)&&(iperturb[0]>0)){
    iperturbsav=iperturb[0];
    iperturb[0]=0;
  }

  iout=-1;
  NNEW(v,double,mt**nk);
  NNEW(fn,double,mt**nk);
  NNEW(stx,double,6*mi[0]**ne);
  NNEW(inum,ITG,*nk);
  results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	  ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
	  prestr,iprestr,filab,eme,emn,een,iperturb,
	  f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	  ndirboun,xbounact,nboun,ipompc,
	  nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,veold,accold,
	  &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	  xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	  &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
	  emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
	  iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
	  fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
	  &reltime,&ne0,thicke,shcon,nshcon,
	  sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
	  mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
	  islavsurf,ielprop,prop,energyini,energy,&kscale,iponoel,
          inoel,nener,orname,&network,ipobody,xbodyact,ibody,typeboun,
	  itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
	  islavelinv,autloc,irowtloc,jqtloc,&nboun2,
	  ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
	  labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
	  &intscheme,physcon);
  SFREE(v);SFREE(fn);SFREE(stx);SFREE(inum);
  iout=1;

  if((*nmethod==1)&&(iglob<0)&&(iperturb[0]>0)){
    iperturb[0]=iperturbsav;
  }
  
  /* determining the system matrix and the external forces */

  NNEW(ad,double,*neq);
  NNEW(fext,double,*neq);

  if(*nmethod==11){
      
    /* determining the nodes and the degrees of freedom in those nodes
       belonging to the substructure */
      
    NNEW(iretain,ITG,*nk);
    NNEW(noderetain,ITG,*nk);
    NNEW(ndirretain,ITG,*nk);
    nretain=0;
      
    for(i=0;i<*nboun;i++){
      if(strcmp1(&typeboun[i],"C")==0){
	iretain[nretain]=i+1;
	noderetain[nretain]=nodeboun[i];
	ndirretain[nretain]=ndirboun[i];
	nretain++;
      }
    }
 
    /* nretain!=0: submatrix application
       nretain==0: Green function application */
      
    if(nretain>0){
      RENEW(iretain,ITG,nretain);
      RENEW(noderetain,ITG,nretain);
      RENEW(ndirretain,ITG,nretain);
    }else{
      SFREE(iretain);SFREE(noderetain);SFREE(ndirretain);
    }
      
    /* creating the right size au */

    NNEW(au,double,nzs[2]);
    rhsi=0;
    //      nmethodl=2;
    nmethodl=*nmethod;

    /* providing for the mass matrix in case of Green functions */

    if(nretain==0){
      mass[0]=1.;
      NNEW(adb,double,*neq);
      NNEW(aub,double,nzs[1]);
    }

  }else{

    /* linear static calculation */

    NNEW(au,double,*nzs);
    nmethodl=*nmethod;

    /* if submodel calculation with a global model obtained by
       a *FREQUENCY calculation: replace stiffness matrix K by
       K-sigma*M */

    if(iglob<0){
      mass[0]=1;
      NNEW(adb,double,*neq);
      NNEW(aub,double,nzs[1]);
    }
	  
  }

  mafillsmmain(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xbounact,nboun,
	       ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
	       nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,
	       nbody,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,&nmethodl,
	       ikmpc,ilmpc,ikboun,ilboun,
	       elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	       ielorien,norien,orab,ntmat_,
	       t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
	       nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
	       xstiff,npmat_,&dtime,matname,mi,
	       ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,physcon,
	       shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,&coriolis,
	       ibody,xloadold,&reltime,veold,springarea,nstate_,
	       xstateini,xstate,thicke,integerglob,doubleglob,
	       tieset,istartset,iendset,ialset,ntie,&nasym,pslavsurf,
	       pmastsurf,mortar,clearini,ielprop,prop,&ne0,fnext,&kscale,
	       iponoel,inoel,&network,ntrans,inotr,trab,smscale,&mscalmethod,
	       set,nset,islavelinv,autloc,irowtloc,jqtloc,&mortartrafoflag);

  /* check for negative Jacobians */

  if(nmethodl==0) *nmethod=0;

  if(nasym==1){
    RENEW(au,double,2*nzs[1]);
    symmetryflag=2;
    inputformat=1;
      
    mafillsmasmain(co,nk,kon,ipkon,lakon,ne,nodeboun,
		   ndirboun,xbounact,nboun,
		   ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
		   nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,
		   nbody,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,
		   nmethod,ikmpc,ilmpc,ikboun,ilboun,
		   elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		   ielmat,ielorien,norien,orab,ntmat_,
		   t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
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

  /* determining the right hand side */

  NNEW(b,double,*neq);
  for(k=0;k<*neq;++k){
    b[k]=fext[k]-f[k];
  }
  SFREE(fext);SFREE(f);

  /* generation of a substructure stiffness matrix (nretain>0) or treating
     Green functions (nretain=0) */

  if(*nmethod==11){

    /* recovering omega_0^2 for Green applications */

    if(nretain==0){
      if(*nforc>0){sigma=xforc[0];}
    }

    /* factorizing the matrix */

    if(*neq>0){
      if(*isolver==0){
#ifdef SPOOLES
	spooles_factor(ad,au,adb,aub,&sigma,icol,irow,neq,nzs,&symmetryflag,
		       &inputformat,&nzs[2]);
#else
	printf(" *ERROR in linstatic: the SPOOLES library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,neq,nzs,
		       &symmetryflag,&inputformat,jq,&nzs[2]);
#else
	printf(" *ERROR in linstatic: the PARDISO library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==8){
#ifdef PASTIX
	pastix_factor_main(ad,au,adb,aub,&sigma,icol,irow,neq,nzs,
		      &symmetryflag,&inputformat,jq,&nzs[2]);
#else
	printf(" *ERROR in linstatic: the PASTIX library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
    }

    /* solving the system of equations with appropriate rhs */

    if(nretain>0){

      /* substructure calculations */

      NNEW(submatrix,double,nretain*nretain);
	  
      for(i=0;i<nretain;i++){
	DMEMSET(b,0,*neq,0.);
	ic=*neq+iretain[i]-1;
	for(j=jq[ic]-1;j<jq[ic+1]-1;j++){
	  ir=irow[j]-1;
	  b[ir]-=au[j];
	}
	      
	/* solving the system */
	      
	if(*neq>0){
	  if(*isolver==0){
#ifdef SPOOLES
	    spooles_solve(b,neq);
#endif
	  }
	  else if(*isolver==7){
#ifdef PARDISO
	    pardiso_solve(b,neq,&symmetryflag,&inputformat,&nrhs);
#endif
		      
	  }
	  else if(*isolver==8){
#ifdef PASTIX
	    pastix_solve(b,neq,&symmetryflag,&nrhs);
#endif
		      
	  }
	}
	      
	/* calculating the internal forces */
	      
	NNEW(v,double,mt**nk);
	NNEW(fn,double,mt**nk);
	NNEW(stn,double,6**nk);
	NNEW(inum,ITG,*nk);
	NNEW(stx,double,6*mi[0]**ne);
	      
	if(strcmp1(&filab[261],"E   ")==0) NNEW(een,double,6**nk);
	if(strcmp1(&filab[2697],"ME  ")==0) NNEW(emn,double,6**nk);
	if(strcmp1(&filab[522],"ENER")==0) NNEW(enern,double,*nk);
	      
	NNEW(eei,double,6*mi[0]**ne);
	if(*nener==1){
	  NNEW(stiini,double,6*mi[0]**ne);
	  NNEW(emeini,double,6*mi[0]**ne);
	  NNEW(enerini,double,2*mi[0]**ne);}
	      
	/* replacing the appropriate boundary value by unity */
	      
	xbounact[iretain[i]-1]=1.;
	      
	results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
		elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
		prestr,iprestr,filab,eme,emn,een,iperturb,
		f,fn,nactdof,&iout,qa,vold,b,nodeboun,ndirboun,
		xbounact,nboun,ipompc,
		nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,veold,
		accold,&bet,
		&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
		ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
		xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
		ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
		nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,
		&ne0,thicke,shcon,nshcon,
		sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
		islavsurf,ielprop,prop,energyini,energy,&kscale,iponoel,
		inoel,nener,orname,&network,ipobody,xbodyact,ibody,typeboun,
		itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
		islavelinv,autloc,irowtloc,jqtloc,&nboun2,
		ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
		labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
		&intscheme,physcon);
	      
	xbounact[iretain[i]-1]=0.;
	      
	SFREE(v);SFREE(stn);SFREE(inum);SFREE(stx);
	      
	if(strcmp1(&filab[261],"E   ")==0) SFREE(een);
	if(strcmp1(&filab[2697],"ME  ")==0) SFREE(emn);
	if(strcmp1(&filab[522],"ENER")==0) SFREE(enern);
	      
	SFREE(eei);if(*nener==1){SFREE(stiini);SFREE(emeini);SFREE(enerini);}
	      
	/* storing the internal forces in the substructure
	   stiffness matrix */
	      
	for(j=0;j<nretain;j++){
	  submatrix[i*nretain+j]=fn[mt*(noderetain[j]-1)+ndirretain[j]];
	}
	      
	SFREE(fn);
	      
      }

      /* cleanup */
      
      if(*isolver==0){
#ifdef SPOOLES
	spooles_cleanup();
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	pardiso_cleanup(&neq[0],&symmetryflag,&inputformat);
#endif
      }
      else if(*isolver==8){
#ifdef PASTIX
#endif
      }
      
      SFREE(iretain);
	  
      FORTRAN(writesubmatrix,(submatrix,noderetain,ndirretain,&nretain,jobnamec));
	  
      SFREE(submatrix);SFREE(noderetain);SFREE(ndirretain);

    }else{

      /* Green function applications (nretain=0) */

      /* storing omega_0^2 into d */

      NNEW(d,double,*nforc);
      for(i=0;i<*nforc;i++){d[i]=xforc[0];}

      strcpy2(fneig,jobnamec,132);
      strcat(fneig,".eig");
      
      if((f2=fopen(fneig,"wb"))==NULL){
	printf(" *ERROR in linstatic: cannot open eigenvalue file for writing...");
	      
	exit(0);
      }
	  
      /* storing a zero as indication that this was not a
	 cyclic symmetry calculation */
	  
      if(fwrite(&zero,sizeof(ITG),1,f2)!=1){
	printf(" *ERROR saving the cyclic symmetry flag to the eigenvalue file...");
	exit(0);
      }
	  
      /* Hermitian */
	  
      if(fwrite(&nherm,sizeof(ITG),1,f2)!=1){
	printf(" *ERROR saving the Hermitian flag to the eigenvalue file...");
	exit(0);
      }

      /* perturbation parameter iperturb[0] */
	  
      if(fwrite(&iperturb[0],sizeof(ITG),1,f2)!=1){
	printf(" *ERROR saving the perturbation flag to the eigenvalue file...");
	exit(0);
      }
	      
      /* reference displacements */

      if(iperturb[0]==1){
	if(fwrite(vold,sizeof(double),mt**nk,f2)!=mt**nk){
	  printf(" *ERROR saving the reference displacements to the eigenvalue file...");
	  exit(0);
	}
      }
	  
      /* storing the number of eigenvalues */
	  
      if(fwrite(&nev,sizeof(ITG),1,f2)!=1){
	printf(" *ERROR saving the number of eigenvalues to the eigenvalue file...");
	exit(0);
      }
	  
      /* the eigenfrequencies are stored as radians/time */
	  
      if(fwrite(d,sizeof(double),nev,f2)!=nev){
	printf(" *ERROR saving the eigenfrequencies to the eigenvalue file...");
	exit(0);
      }
	  
      /* storing the stiffness matrix */
	  
      if(fwrite(ad,sizeof(double),neq[1],f2)!=neq[1]){
	printf(" *ERROR saving the diagonal of the stiffness matrix to the eigenvalue file...");
	exit(0);
      }
      if(fwrite(au,sizeof(double),nzs[2],f2)!=nzs[2]){
	printf(" *ERROR saving the off-diagonal terms of the stiffness matrix to the eigenvalue file...");
	exit(0);
      }
	  
      /* storing the mass matrix */
	  
      if(fwrite(adb,sizeof(double),neq[1],f2)!=neq[1]){
	printf(" *ERROR saving the diagonal of the mass matrix to the eigenvalue file...");
	exit(0);
      }
      if(fwrite(aub,sizeof(double),nzs[1],f2)!=nzs[1]){
	printf(" *ERROR saving the off-diagonal terms of the mass matrix to the eigenvalue file...");
	exit(0);
      }

      SFREE(d);
	  
      /* calculating each Green function */

      for(i=0;i<*nforc;i++){
	node=nodeforc[2*i];
	idir=ndirforc[i];

	/* check whether degree of freedom is active */

	if(nactdof[mt*(node-1)+idir]==0){
	  printf(" *ERROR in linstatic: degree of freedom corresponding to node %d \n and direction %d is not active: no unit force can be applied\n",node,idir);
	  FORTRAN(stop,());
	}

	/* defining a unit force on the rhs */

	DMEMSET(b,0,*neq,0.);
	b[nactdof[mt*(node-1)+idir]-1]=1.;
	      
	/* solving the system */
	      
	if(*neq>0){
	  if(*isolver==0){
#ifdef SPOOLES
	    spooles_solve(b,neq);
#endif
	  }
	  else if(*isolver==7){
#ifdef PARDISO
	    pardiso_solve(b,neq,&symmetryflag,&inputformat,&nrhs);
#endif
		      
	  }
	  else if(*isolver==8){
#ifdef PASTIX
	    pastix_solve(b,neq,&symmetryflag,&nrhs);
#endif
		      
	  }
	}
	      
	/* storing the Green function */

	if(fwrite(b,sizeof(double),*neq,f2)!=*neq){
	  printf(" *ERROR saving data to the eigenvalue file...");
	  exit(0);
	}

	/* calculating the displacements and the stresses and storing */
	/* the results in frd format for each valid eigenmode */
	      
	NNEW(v,double,mt**nk);
	NNEW(fn,double,mt**nk);
	NNEW(stn,double,6**nk);
	NNEW(inum,ITG,*nk);
	NNEW(stx,double,6*mi[0]**ne);
	      
	if(strcmp1(&filab[261],"E   ")==0) NNEW(een,double,6**nk);
	if(strcmp1(&filab[2697],"ME  ")==0) NNEW(emn,double,6**nk);
	if(strcmp1(&filab[522],"ENER")==0) NNEW(enern,double,*nk);
	if(strcmp1(&filab[2175],"CONT")==0) NNEW(cdn,double,6**nk);
	      
	NNEW(eei,double,6*mi[0]**ne);
	if(*nener==1){
	  NNEW(stiini,double,6*mi[0]**ne);
	  NNEW(emeini,double,6*mi[0]**ne);
	  NNEW(enerini,double,2*mi[0]**ne);}
	      
	results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
		elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
		prestr,iprestr,filab,eme,emn,een,iperturb,
		f,fn,nactdof,&iout,qa,vold,b,nodeboun,ndirboun,
		xbounact,nboun,ipompc,
		nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,veold,
		accold,&bet,
		&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
		ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
		xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
		ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
		nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,
		&ne0,thicke,shcon,nshcon,
		sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
		islavsurf,ielprop,prop,energyini,energy,&kscale,iponoel,
		inoel,nener,orname,&network,ipobody,xbodyact,ibody,typeboun,
		itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
		islavelinv,autloc,irowtloc,jqtloc,&nboun2,
		ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
		labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
		&intscheme,physcon);
	      
	SFREE(eei);
	if(*nener==1){
	  SFREE(stiini);SFREE(emeini);SFREE(enerini);}
	      
	/*	      memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);
		      memcpy(&sti[0],&stx[0],sizeof(double)*6*mi[0]*ne0);*/
	      
	++*kode;
	time=1.*i;
	      
	/* for cyclic symmetric sectors: duplicating the results */
	      
	if(*mcs>0){
	  ptime=*ttime+time;
	  frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,t1act,
		 fn,&ptime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
		 nstate_,istep,&iinc,iperturb,ener,mi,output,ithermal,
		 qfn,ialset,istartset,iendset,trab,inotr,ntrans,orab,
		 ielorien,norien,sti,veold,&noddiam,set,nset,emn,thicke,
		 jobnamec,&ne0,cdn,mortar,nmat,qfx,ielprop,prop);
	}
	else{
	  if(strcmp1(&filab[1044],"ZZS")==0){
	    NNEW(neigh,ITG,40**ne);
	    NNEW(ipneigh,ITG,*nk);
	  }
	  ptime=*ttime+time;
	  frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	      kode,filab,een,t1act,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	      nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	      ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	      mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	      cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	      thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni,nmat,
	      ielprop,prop,sti);
	  if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);}
	}
	      
	SFREE(v);SFREE(stn);SFREE(inum);
	SFREE(stx);SFREE(fn);
	      
	if(strcmp1(&filab[261],"E   ")==0) SFREE(een);
	if(strcmp1(&filab[2697],"ME  ")==0) SFREE(emn);
	if(strcmp1(&filab[522],"ENER")==0) SFREE(enern);
	if(strcmp1(&filab[2175],"CONT")==0) SFREE(cdn);
	      
      }


      fclose(f2);

    }
	  
    SFREE(au);SFREE(ad);SFREE(b);
      
    SFREE(xbounact);SFREE(xforcact);SFREE(xloadact);SFREE(t1act);SFREE(ampli);
    SFREE(xbodyact);

    //    if(*nbody>0) SFREE(ipobody);

    SFREE(xstiff);
      
    if(iglob!=0){SFREE(integerglob);SFREE(doubleglob);}
      
    return;


  }else if(*nmethod!=0){

    /* linear static applications */

    if(*isolver==0){
#ifdef SPOOLES
      spooles(ad,au,adb,aub,&sigma,b,icol,irow,neq,nzs,&symmetryflag,
              &inputformat,&nzs[2]);
#else
      printf(" *ERROR in linstatic: the SPOOLES library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if((*isolver==2)||(*isolver==3)){
      if(nasym>0){
	printf(" *ERROR in nonlingeo: the iterative solver cannot be used for asymmetric matrices\n\n");
	FORTRAN(stop,());
      }
      preiter(ad,&au,b,&icol,&irow,neq,nzs,isolver,iperturb);
    }
    else if(*isolver==4){
#ifdef SGI
      if(nasym>0){
	printf(" *ERROR in nonlingeo: the SGI solver cannot be used for asymmetric matrices\n\n");
	FORTRAN(stop,());
      }
      token=1;
      sgi_main(ad,au,adb,aub,&sigma,b,icol,irow,neq,nzs,token);
#else
      printf(" *ERROR in linstatic: the SGI library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==5){
#ifdef TAUCS
      if(nasym>0){
	printf(" *ERROR in nonlingeo: the TAUCS solver cannot be used for asymmetric matrices\n\n");
	FORTRAN(stop,());
      }
      tau(ad,&au,adb,aub,&sigma,b,icol,&irow,neq,nzs);
#else
      printf(" *ERROR in linstatic: the TAUCS library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==7){
#ifdef PARDISO
      pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,neq,nzs,
		   &symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
#else
      printf(" *ERROR in linstatic: the PARDISO library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==8){
#ifdef PASTIX
      pastix_main(ad,au,adb,aub,&sigma,b,icol,irow,neq,nzs,
		  &symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
#else
      printf(" *ERROR in linstatic: the PASTIX library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }

    /* saving of ad and au for sensitivity analysis */

    for(i=0;i<*ntie;i++){
      if(strcmp1(&tieset[i*243+80],"D")==0){
	    
	strcpy2(stiffmatrix,jobnamec,132);
	strcat(stiffmatrix,".stm");
	    
	if((f1=fopen(stiffmatrix,"wb"))==NULL){
	  printf(" *ERROR in linstatic: cannot open stiffness matrix file for writing...");
	  exit(0);
	}
	    
	/* storing the stiffness matrix */

	/* nzs,irow,jq and icol have to be stored too, since the static analysis
	   can involve contact, whereas in the sensitivity analysis contact is not
	   taken into account while determining the structure of the stiffness
	   matrix (in mastruct.c)
	*/
	    
	if(fwrite(&nasym,sizeof(ITG),1,f1)!=1){
	  printf(" *ERROR saving the symmetry flag to the stiffness matrix file...");
	  exit(0);
	}
	if(fwrite(nzs,sizeof(ITG),3,f1)!=3){
	  printf(" *ERROR saving the number of subdiagonal nonzeros to the stiffness matrix file...");
	  exit(0);
	}
	if(fwrite(irow,sizeof(ITG),nzs[2],f1)!=nzs[2]){
	  printf(" *ERROR saving irow to the stiffness matrix file...");
	  exit(0);
	}
	if(fwrite(jq,sizeof(ITG),neq[1]+1,f1)!=neq[1]+1){
	  printf(" *ERROR saving jq to the stiffness matrix file...");
	  exit(0);
	}
	if(fwrite(icol,sizeof(ITG),neq[1],f1)!=neq[1]){
	  printf(" *ERROR saving icol to the stiffness matrix file...");
	  exit(0);
	}
	if(fwrite(ad,sizeof(double),neq[1],f1)!=neq[1]){
	  printf(" *ERROR saving the diagonal of the stiffness matrix to the stiffness matrix file...");
	  exit(0);
	}
	if(fwrite(au,sizeof(double),nzs[2],f1)!=nzs[2]){
	  printf(" *ERROR saving the off-diagonal terms of the stiffness matrix to the tiffness matrix file...");
	  exit(0);
	}
	fclose(f1);

	break;
      }
    }
    
    SFREE(ad);SFREE(au);
    if(iglob<0){SFREE(adb);SFREE(aub);}

    /* calculating the displacements and the stresses and storing */
    /* the results in frd format for each valid eigenmode */

    NNEW(v,double,mt**nk);
    NNEW(fn,double,mt**nk);
    NNEW(stn,double,6**nk);
    NNEW(inum,ITG,*nk);
    NNEW(stx,double,6*mi[0]**ne);
  
    if(strcmp1(&filab[261],"E   ")==0) NNEW(een,double,6**nk);
    if(strcmp1(&filab[2697],"ME  ")==0) NNEW(emn,double,6**nk);
    if(strcmp1(&filab[522],"ENER")==0) NNEW(enern,double,*nk);
    if(strcmp1(&filab[2175],"CONT")==0) NNEW(cdn,double,6**nk);

    NNEW(eei,double,6*mi[0]**ne);
    if(*nener==1){
      NNEW(stiini,double,6*mi[0]**ne);
      NNEW(emeini,double,6*mi[0]**ne);
      NNEW(enerini,double,2*mi[0]**ne);}

    results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	    elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	    ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
	    prestr,iprestr,filab,eme,emn,een,iperturb,
            f,fn,nactdof,&iout,qa,vold,b,nodeboun,ndirboun,xbounact,nboun,
	    ipompc,
	    nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,veold,accold,&bet,
            &gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	    xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
            ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
            xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
            ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	    nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,
            &ne0,thicke,shcon,nshcon,
            sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
            mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
	    islavsurf,ielprop,prop,energyini,energy,&kscale,iponoel,
            inoel,nener,orname,&network,ipobody,xbodyact,ibody,typeboun,
	    itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
	    islavelinv,autloc,irowtloc,jqtloc,&nboun2,
	    ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
	    labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
	    &intscheme,physcon);

    SFREE(eei);
    if(*nener==1){
      SFREE(stiini);SFREE(emeini);SFREE(enerini);}

    memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);
    memcpy(&sti[0],&stx[0],sizeof(double)*6*mi[0]*ne0);

    ++*kode;

    /* for cyclic symmetric sectors: duplicating the results */

    if(*mcs>0){
      ptime=*ttime+time;
      frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,t1act,
	     fn,&ptime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
	     nstate_,istep,&iinc,iperturb,ener,mi,output,ithermal,
	     qfn,ialset,istartset,iendset,trab,inotr,ntrans,orab,
	     ielorien,norien,sti,veold,&noddiam,set,nset,emn,thicke,
	     jobnamec,&ne0,cdn,mortar,nmat,qfx,ielprop,prop);
    }
    else{
      if(strcmp1(&filab[1044],"ZZS")==0){
	NNEW(neigh,ITG,40**ne);
	NNEW(ipneigh,ITG,*nk);
      }
      ptime=*ttime+time;
      frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	  kode,filab,een,t1act,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	  nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	  mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	  thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni,nmat,ielprop,
	  prop,sti);
      if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);}
    }

    /* updating the .sta file */

    iitsta=1;
    FORTRAN(writesta,(istep,&iinc,&icutb,&iitsta,ttime,&time,&dtime));

    SFREE(v);SFREE(stn);SFREE(inum);
    SFREE(b);SFREE(stx);SFREE(fn);

    if(strcmp1(&filab[261],"E   ")==0) SFREE(een);
    if(strcmp1(&filab[2697],"ME  ")==0) SFREE(emn);
    if(strcmp1(&filab[522],"ENER")==0) SFREE(enern);
    if(strcmp1(&filab[2175],"CONT")==0) SFREE(cdn);

  }
  else {

    /* error occurred in mafill: storing the geometry in frd format */

    ++*kode;
    NNEW(inum,ITG,*nk);for(k=0;k<*nk;k++) inum[k]=1;
    if(strcmp1(&filab[1044],"ZZS")==0){
      NNEW(neigh,ITG,40**ne);
      NNEW(ipneigh,ITG,*nk);
    }
    ptime=*ttime+time;
    frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	kode,filab,een,t1,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni,nmat,ielprop,
	prop,sti);
    if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);}
    SFREE(inum);FORTRAN(stop,());

  }

  if(*mortar>-2){
    if(ncont!=0){
      *ne=ne0;
      if(*nener==1) RENEW(ener,double,2*mi[0]**ne);
      RENEW(ipkon,ITG,*ne);
      RENEW(lakon,char,8**ne);
      RENEW(kon,ITG,*nkon);
      if(*norien>0){
	RENEW(ielorien,ITG,mi[2]**ne);
      }
      RENEW(ielmat,ITG,mi[2]**ne);
      SFREE(cg);SFREE(straight);
      SFREE(imastop);SFREE(itiefac);SFREE(islavnode);SFREE(islavsurf);
      SFREE(nslavnode);SFREE(iponoels);SFREE(inoels);SFREE(imastnode);
      SFREE(nmastnode);SFREE(itietri);SFREE(koncont);SFREE(xnoels);
      SFREE(springarea);SFREE(xmastnor);

      if(*mortar==0){
	SFREE(areaslav);
      }else if(*mortar==1){
	SFREE(pmastsurf);SFREE(ipe);SFREE(ime);SFREE(pslavsurf);
	SFREE(islavact);SFREE(clearini);
      }
    }
    mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
    mpcinfo[3]=maxlenmpc;
  }

  /* updating the loading at the end of the step; 
     important in case the amplitude at the end of the step
     is not equal to one */

  for(k=0;k<*nboun;++k){xbounold[k]=xbounact[k];}
  for(k=0;k<*nforc;++k){xforcold[k]=xforcact[k];}
  for(k=0;k<2**nload;++k){xloadold[k]=xloadact[k];}
  for(k=0;k<7**nbody;k=k+7){xbodyold[k]=xbodyact[k];}
  if(*ithermal==1){
    for(k=0;k<*nk;++k){t1old[k]=t1act[k];}
    for(k=0;k<*nk;++k){vold[mt*k]=t1act[k];}
  }

  SFREE(xbounact);SFREE(xforcact);SFREE(xloadact);SFREE(t1act);SFREE(ampli);
  SFREE(xbodyact);

  //  if(*nbody>0) SFREE(ipobody);

  SFREE(xstiff);

  if(iglob!=0){SFREE(integerglob);SFREE(doubleglob);}

  *irowp=irow;*enerp=ener;*xstatep=xstate;*ipkonp=ipkon;*lakonp=lakon;
  *konp=kon;*ielmatp=ielmat;*ielorienp=ielorien;*icolp=icol;

  (*ttime)+=(*tper);
 
  return;
}
