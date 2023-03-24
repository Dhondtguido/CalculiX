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
#ifdef PARDISO
#include "pardiso.h"
#endif
#ifdef PASTIX
#include "pastix.h"
#endif


void electromagnetics(double **cop,ITG *nk,ITG **konp,ITG **ipkonp,
		      char **lakonp,ITG *ne,
		      ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
		      ITG **ipompcp,ITG **nodempcp,double **coefmpcp,
		      char **labmpcp,
		      ITG *nmpc,
		      ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
		      ITG **nelemloadp,char **sideloadp,double *xload,
		      ITG *nload,
		      ITG *nactdof,
		      ITG **icolp,ITG **jqp,ITG **irowp,ITG *neq,ITG *nzl,
		      ITG *nmethod,ITG **ikmpcp,ITG **ilmpcp,ITG *ikboun,
		      ITG *ilboun,
		      double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
		      double *alcon,ITG *nalcon,double *alzero,ITG **ielmatp,
		      ITG **ielorienp,ITG *norien,double *orab,ITG *ntmat_,
		      double *t0,double *t1,double *t1old,
		      ITG *ithermal,double *prestr,ITG *iprestr,
		      double **voldp,ITG *iperturb,double *sti,ITG *nzs, 
		      ITG *kode,char *filab,
		      ITG *idrct,ITG *jmax,ITG *jout,double *timepar,
		      double *eme,
		      double *xbounold,double *xforcold,double *xloadold,
		      double *veold,double *accold,
		      char *amname,double *amta,ITG *namta,ITG *nam,
		      ITG *iamforc,ITG **iamloadp,
		      ITG *iamt1,double *alpha,ITG *iexpl,
		      ITG *iamboun,double *plicon,ITG *nplicon,double *plkcon,
		      ITG *nplkcon,
		      double **xstatep,ITG *npmat_,ITG *istep,double *ttime,
		      char *matname,double *qaold,ITG *mi,
		      ITG *isolver,ITG *ncmat_,ITG *nstate_,
		      double *cs,ITG *mcs,ITG *nkon,double **enerp,ITG *mpcinfo,
		      char *output,
		      double *shcon,ITG *nshcon,double *cocon,ITG *ncocon,
		      double *physcon,ITG *nflow,double *ctrl,
		      char **setp,ITG *nset,ITG **istartsetp,
		      ITG **iendsetp,ITG **ialsetp,ITG *nprint,char *prlab,
		      char *prset,ITG *nener,ITG *ikforc,ITG *ilforc,
		      double *trab,
		      ITG *inotr,ITG *ntrans,double **fmpcp,ITG *ipobody,
		      ITG *ibody,double *xbody,ITG *nbody,double *xbodyold,
		      ITG *ielprop,double *prop,ITG *ntie,char **tiesetp,
		      ITG *itpamp,ITG *iviewfile,char *jobnamec,
		      double **tietolp,
		      ITG *nslavs,double *thicke,ITG *ics,ITG *nalset,
		      ITG *nmpc_,
		      ITG *nmat,char *typeboun,ITG *iaxial,ITG *nload_,
		      ITG *nprop,
		      ITG *network,char *orname,double *t0g,double *t1g){

  char description[13]="            ",*lakon=NULL,jobnamef[396]="",
    *labmpc=NULL,kind1[2]="E",kind2[2]="E",*set=NULL,*tieset=NULL,
    cflag[1]=" ",*sideloadref=NULL,*sideload=NULL,*env,*envsys,
    *labmpc2=NULL; 
 
  ITG *inum=NULL,k,iout=0,icntrl,iinc=0,jprint=0,iit=-1,jnz=0,
    icutb=0,istab=0,ifreebody,uncoupled=0,maxfaces,indexe,nope,
    iperturb_sav[2],*icol=NULL,*irow=NULL,ielas=0,icmd=0,j,
    memmpc_,mpcfree,icascade,maxlenmpc,*nodempc=NULL,*iaux=NULL,
    *itg=NULL,*ineighe=NULL,null=0,iactive[3],neqterms,ntflag,
    *ieg=NULL,ntg=0,ntr,*kontri=NULL,*nloadtr=NULL,index,
    *ipiv=NULL,ntri,mode=-1,noddiam=-1,nasym=0,
    ntrit,*inocs=NULL,*nacteq=NULL,kkv,kkx,kk6,
    *nactdog=NULL,nteq,nmastnode,imast,massact[2],
    *ipkon=NULL,*kon=NULL,*ielorien=NULL,nmethodact,ne2=0,
    *ielmat=NULL,inext,itp=0,symmetryflag=0,inputformat=0,
    iitterm=0,ngraph=1,ithermalact=2,*islavact=NULL,
    *ipompc=NULL,*ikmpc=NULL,*ilmpc=NULL,i0ref,irref,icref,
    *islavnode=NULL,*imastnode=NULL,*nslavnode=NULL,mortar=0,
    mt=mi[1]+1,*nactdofinv=NULL,inode,idir,*islavsurf=NULL,
    iemchange=0,nzsrad,*mast1rad=NULL,*irowrad=NULL,*icolrad=NULL,
    *jqrad=NULL,*ipointerrad=NULL,*integerglob=NULL,im,ne0,
    mass[2]={0,0},stiffness=1,buckling=0,rhsi=1,intscheme=0,idiscon=0,
    coriolis=0,*ipneigh=NULL,*neigh=NULL,i,icfd=0,id,node,networknode,
    iflagact=0,*nodorig=NULL,*ipivr=NULL,*inomat=NULL,*nodface=NULL,
    *ipoface=NULL,*istartset=NULL,*iendset=NULL,*ialset=NULL,
    *nelemloadref=NULL,*iamloadref=NULL,nloadref,kscale=1,
    *nelemload=NULL,*iamload=NULL,*idefload=NULL,ialeatoric=0,
    *iponoel=NULL,*inoel=NULL,inoelsize,nrhs=1,neqfreq,nzsfreq,
    *irowfreq=NULL,*icolfreq=NULL,*jqfreq=NULL,*jq=NULL,
    *itiefac=NULL,mscalmethod=0,nkon0,*nintpoint=0,num_cpus,sys_cpus,
    *islavelinv=NULL,*irowtloc=NULL,*jqtloc=NULL,nboun2,
    *ndirboun2=NULL,*nodeboun2=NULL,nmpc2,*ipompc2=NULL,*nodempc2=NULL,
    *ikboun2=NULL,*ilboun2=NULL,*ikmpc2=NULL,*ilmpc2=NULL,mortartrafoflag=0;

  double *stn=NULL,*v=NULL,*een=NULL,cam[5],*epn=NULL,*cdn=NULL,
    *f=NULL,*fn=NULL,qa[4]={0.,0.,-1.,0.},qam[2]={0.,0.},dtheta,theta,
    err,ram[4]={0.,0.,0.,0.},*springarea=NULL,*h0=NULL,
    ram1[2]={0.,0.},ram2[2]={0.,0.},deltmx,*clearini=NULL,
    uam[2]={0.,0.},*vini=NULL,*ac=NULL,qa0,qau,ea,ptime,
    *t1act=NULL,qamold[2],*xbounact=NULL,*bc=NULL,
    *xforcact=NULL,*xloadact=NULL,*fext=NULL,h0scale=1.,
    reltime,time,bet=0.,gam=0.,*aux1=NULL,*aux2=NULL,dtime,*fini=NULL,
    *fextini=NULL,*veini=NULL,*xstateini=NULL,*h0ref=NULL,
    *ampli=NULL,*eei=NULL,*t1ini=NULL,*tinc,*tper,*tmin,*tmax,
    *xbounini=NULL,*xstiff=NULL,*stx=NULL,*cv=NULL,*cvini=NULL,
    *enern=NULL,*coefmpc=NULL,*xstaten=NULL,
    *enerini=NULL,*emn=NULL,*xmastnor=NULL,*fnext=NULL,
    *tarea=NULL,*tenv=NULL,*erad=NULL,*fnr=NULL,*fni=NULL,
    *adview=NULL,*auview=NULL,*qfx=NULL,*adaux=NULL,
    *qfn=NULL,*co=NULL,*vold=NULL,*fenv=NULL,sigma=0.,
    *xbodyact=NULL,*cgr=NULL,dthetaref,*vr=NULL,*vi=NULL,
    *stnr=NULL,*stni=NULL,*vmax=NULL,*stnmax=NULL,*fmpc=NULL,*ener=NULL,
    *f_cm=NULL,*f_cs=NULL,*tietol=NULL,
    *xstate=NULL,*eenmax=NULL,*adrad=NULL,*aurad=NULL,*bcr=NULL,
    *emeini=NULL,*doubleglob=NULL,*au=NULL,
    *ad=NULL,*b=NULL,*aub=NULL,*adb=NULL,*pslavsurf=NULL,*pmastsurf=NULL,
    *cdnr=NULL,*cdni=NULL,*energyini=NULL,*energy=NULL,*adfreq=NULL,
    *aufreq=NULL,*bfreq=NULL,om,*smscale=NULL,*autloc=NULL,*xboun2=NULL,
    *coefmpc2=NULL;

#ifdef SGI
  ITG token;
#endif

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
    sys_cpus=getSystemCPUs();
    if(sys_cpus<1) sys_cpus=1;
  }
  
  /* else global declaration, if any, applies */
  
  env = getenv("OMP_NUM_THREADS");
  if(num_cpus==0){
    if(env)
      num_cpus=atoi(env);
    if(num_cpus<1) {
      num_cpus=1;
    }else if(num_cpus>sys_cpus){
      num_cpus=sys_cpus;
    }
  }
 
  ne0=*ne;

  /* next line is needed to avoid that elements with negative ipkon
     are taken into account in extrapolate.f */

  strcpy1(&filab[2],"C",1);

  if(*nmethod==8){
    *nmethod=1;
  }else if(*nmethod==9){
    *nmethod=4;
  }else if(*nmethod==10){
    *nmethod=2;
  }
 
  for(k=0;k<3;k++){
    strcpy1(&jobnamef[k*132],&jobnamec[k*132],132);
  }
  
  qa0=ctrl[20];qau=ctrl[21];ea=ctrl[23];deltmx=ctrl[26];
  i0ref=ctrl[0];irref=ctrl[1];icref=ctrl[3];
  
  memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];icascade=mpcinfo[2];
  maxlenmpc=mpcinfo[3];
  
  icol=*icolp;irow=*irowp;jq=*jqp;co=*cop;vold=*voldp;
  ipkon=*ipkonp;lakon=*lakonp;kon=*konp;ielorien=*ielorienp;
  ielmat=*ielmatp;ener=*enerp;xstate=*xstatep;
  
  ipompc=*ipompcp;labmpc=*labmpcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
  fmpc=*fmpcp;nodempc=*nodempcp;coefmpc=*coefmpcp;

  set=*setp;istartset=*istartsetp;iendset=*iendsetp;ialset=*ialsetp;
  tieset=*tiesetp;tietol=*tietolp;

  nelemload=*nelemloadp;iamload=*iamloadp;
  sideload=*sideloadp;

  tinc=&timepar[0];
  tper=&timepar[1];
  tmin=&timepar[2];
  tmax=&timepar[3];
  
  /* invert nactdof */
  
  NNEW(nactdofinv,ITG,mt**nk);
  NNEW(nodorig,ITG,*nk);
  FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
			 ipkon,lakon,kon,ne));
  SFREE(nodorig);
  
  /* allocating a field for the stiffness matrix */
  
  NNEW(xstiff,double,(long long)27*mi[0]**ne);
  
  /* allocating force fields */
  
  NNEW(f,double,neq[1]);
  NNEW(fext,double,neq[1]);
  
  NNEW(b,double,neq[1]);
  NNEW(vini,double,mt**nk);
  
  /* allocating fields for the actual external loading */
  
  NNEW(xbounact,double,*nboun);
  NNEW(xbounini,double,*nboun);
  for(k=0;k<*nboun;++k){xbounact[k]=xbounold[k];}
  NNEW(xforcact,double,*nforc);
  NNEW(xloadact,double,2**nload);
  NNEW(xbodyact,double,7**nbody);
  /* copying the rotation axis and/or acceleration vector */
  for(k=0;k<7**nbody;k++){xbodyact[k]=xbody[k];}
  
  /* assigning the body forces to the elements */ 
  
  /* if(*nbody>0){
    ifreebody=*ne+1;
    NNEW(ipobody,ITG,2*ifreebody**nbody);
    for(k=1;k<=*nbody;k++){
      FORTRAN(bodyforce,(cbody,ibody,ipobody,nbody,set,istartset,
			 iendset,ialset,&inewton,nset,&ifreebody,&k));
      RENEW(ipobody,ITG,2*(*ne+ifreebody));
    }
    RENEW(ipobody,ITG,2*(ifreebody-1));
    }*/
  
  /* for thermal calculations: forced convection and cavity
     radiation*/
  
  if(*ithermal>1){
    NNEW(itg,ITG,*nload+3**nflow);
    NNEW(ieg,ITG,*nflow);
    /* max 6 triangles per face, 4 entries per triangle */
    NNEW(kontri,ITG,24**nload);
    NNEW(nloadtr,ITG,*nload);
    NNEW(nacteq,ITG,4**nk);
    NNEW(nactdog,ITG,4**nk);
    NNEW(v,double,mt**nk);
    FORTRAN(envtemp,(itg,ieg,&ntg,&ntr,sideload,nelemload,
		     ipkon,kon,lakon,ielmat,ne,nload,
		     kontri,&ntri,nloadtr,nflow,ndirboun,nactdog,
		     nodeboun,nacteq,nboun,ielprop,prop,&nteq,
		     v,network,physcon,shcon,ntmat_,co,
		     vold,set,nshcon,rhcon,nrhcon,mi,nmpc,nodempc,
		     ipompc,labmpc,ikboun,&nasym,ttime,&time,
		     iaxial));
    SFREE(v);
      
    if((*mcs>0)&&(ntr>0)){
      NNEW(inocs,ITG,*nk);
      radcyc(nk,kon,ipkon,lakon,ne,cs,mcs,nkon,ialset,istartset,
	     iendset,&kontri,&ntri,&co,&vold,&ntrit,inocs,mi);
    }
    else{ntrit=ntri;}
      
    nzsrad=100*ntr;
    NNEW(mast1rad,ITG,nzsrad);
    NNEW(irowrad,ITG,nzsrad);
    NNEW(icolrad,ITG,ntr);
    NNEW(jqrad,ITG,ntr+1);
    NNEW(ipointerrad,ITG,ntr);
      
    if(ntr>0){
      mastructrad(&ntr,nloadtr,sideload,ipointerrad,
		  &mast1rad,&irowrad,&nzsrad,
		  jqrad,icolrad);
    }
      
    /* determine the network elements belonging to a given node (for usage
       in user subroutine film */

    //      if(ntg>0){
    if((*network>0)||(ntg>0)){
      NNEW(iponoel,ITG,*nk);
      NNEW(inoel,ITG,2**nkon);
      if(*network>0){
	FORTRAN(networkelementpernode,(iponoel,inoel,lakon,ipkon,kon,
				       &inoelsize,nflow,ieg,ne,network));
      }
      RENEW(inoel,ITG,2*inoelsize);
    }
      
    SFREE(ipointerrad);SFREE(mast1rad);
    RENEW(irowrad,ITG,nzsrad);
      
    RENEW(itg,ITG,ntg);
    NNEW(ineighe,ITG,ntg);
    RENEW(kontri,ITG,4*ntrit);
    RENEW(nloadtr,ITG,ntr);
      
    NNEW(adview,double,ntr);
    NNEW(auview,double,2*nzsrad);
    NNEW(tarea,double,ntr);
    NNEW(tenv,double,ntr);
    NNEW(fenv,double,ntr);
    NNEW(erad,double,ntr);
      
    NNEW(ac,double,nteq*nteq);
    NNEW(bc,double,nteq);
    NNEW(ipiv,ITG,nteq);
    NNEW(adrad,double,ntr);
    NNEW(aurad,double,2*nzsrad);
    NNEW(bcr,double,ntr);
    NNEW(ipivr,ITG,ntr);
  }
  if(*ithermal>1){NNEW(qfx,double,3*mi[0]**ne);}
  
  /* allocating a field for the instantaneous amplitude */
  
  NNEW(ampli,double,*nam);
  
  NNEW(fini,double,neq[1]);
  
  /* allocating fields for nonlinear dynamics */
  
  if((*nmethod==2)||(*nmethod==4)){
    mass[0]=1;
    mass[1]=1;
    NNEW(aux2,double,neq[1]);
    NNEW(fextini,double,neq[1]);
    NNEW(cv,double,neq[1]);
    NNEW(cvini,double,neq[1]);
    NNEW(veini,double,mt**nk);
    NNEW(adb,double,neq[1]);
    NNEW(aub,double,nzs[1]);
  }

  qa[0]=qaold[0];
  qa[1]=qaold[1];
  
  /* normalizing the time */
  
  FORTRAN(checktime,(itpamp,namta,tinc,ttime,amta,tmin,&inext,&itp,istep,tper));
  dtheta=(*tinc)/(*tper);
  dthetaref=dtheta;
  if(dtheta<=1.e-6){
    printf("\n *ERROR in electromagnetics\n");
    printf(" increment size smaller than one millionth of step size\n");
    printf(" increase increment size\n\n");
  }
  *tmin=*tmin/(*tper);
  *tmax=*tmax/(*tper);
  theta=0.;
  
  /* calculating an initial flux norm */
  
  if(*ithermal!=2){
    if(qau>1.e-10){qam[0]=qau;}
    else if(qa0>1.e-10){qam[0]=qa0;}
    else if(qa[0]>1.e-10){qam[0]=qa[0];}
    else {qam[0]=1.e-2;}
  }
  if(*ithermal>1){
    if(qau>1.e-10){qam[1]=qau;}
    else if(qa0>1.e-10){qam[1]=qa0;}
    else if(qa[1]>1.e-10){qam[1]=qa[1];}
    else {qam[1]=1.e-2;}
  }
  
  /*********************************************************************/
  
  /* calculating the initial quasi-static magnetic intensity due to 
     the coil current */
  
  /*********************************************************************/
  
  /* calculate the current density in the coils

     in this section nload, nforc, nbody and nam are set to zero; the
     electrical potential is supposed to be given (in the form of a
     temperature), the current is calculated (in the form of heat
     flux) by thermal analogy  */
  
  reltime=1.;
  time=0.;
  dtime=0.;
  ithermalact=2;

  nmethodact=1;
  massact[0]=0;
  massact[1]=0;

  if(*ithermal<=1){
    NNEW(qfx,double,3*mi[0]**ne);
    NNEW(t0,double,*nk);
  }
  if(strcmp1(&filab[3567],"ECD ")==0){NNEW(qfn,double,3**nk);}
  
  /* the coil current is assumed to be applied at once, i.e. as 
     step loading; the calculation, however, is a quasi-static
     calculation */

  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,&null,xloadold,xload,
		    xloadact,iamload,&null,ibody,xbody,&null,xbodyold,xbodyact,
		    t1old,t1,t1act,iamt1,nk,amta,
		    namta,&null,ampli,&time,&reltime,ttime,&dtime,&ithermalact,
		    nmethod,
		    xbounold,xboun,xbounact,iamboun,nboun,
		    nodeboun,ndirboun,nodeforc,ndirforc,istep,&iinc,
		    co,vold,itg,&ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
		    ntrans,trab,inotr,veold,integerglob,doubleglob,tieset,
		    istartset,
		    iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc,
		    ipobody,iponoel,inoel,ipkon,kon,ielprop,prop,ielmat,
		    shcon,nshcon,rhcon,nrhcon,cocon,ncocon,ntmat_,lakon,set,
		    nset));
  
  cam[0]=0.;cam[1]=0.;
  
  /* deactivating all elements except the shells */
  
  for(i=0;i<*ne;i++){
    if(strcmp1(&lakon[8*i+6],"L")!=0){
      ipkon[i]=-ipkon[i]-2;
    }
  }
  
  remastruct(ipompc,&coefmpc,&nodempc,nmpc,
	     &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,ikboun,ilboun,
	     labmpc,nk,&memmpc_,&icascade,&maxlenmpc,
	     kon,ipkon,lakon,ne,nactdof,icol,jq,&irow,isolver,
	     neq,nzs,&nmethodact,&f,&fext,&b,&aux2,&fini,&fextini,
	     &adb,&aub,&ithermalact,iperturb,mass,mi,iexpl,&mortar,
             typeboun,&cv,&cvini,&iit,network,itiefac,&ne0,&nkon0,
	     nintpoint,islavsurf,pmastsurf,tieset,ntie,&num_cpus);
  
  /* invert nactdof */
  
  SFREE(nactdofinv);
  NNEW(nactdofinv,ITG,mt**nk);
  NNEW(nodorig,ITG,*nk);
  FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
			 ipkon,lakon,kon,ne));
  SFREE(nodorig);
  
  iout=-1;
  
  NNEW(fn,double,mt**nk);
  NNEW(inum,ITG,*nk);
  NNEW(v,double,mt**nk);

  results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	  ielorien,norien,orab,ntmat_,t0,t1act,&ithermalact,
	  prestr,iprestr,filab,eme,emn,een,iperturb,f,fn,nactdof,&iout,
	  qa,vold,b,nodeboun,ndirboun,xbounact,nboun,ipompc,nodempc,coefmpc,
	  labmpc,nmpc,&nmethodact,cam,&neq[1],veold,accold,&bet,
          &gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	  xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
          ncmat_,nstate_,sti,vini,ikboun,ilboun,ener,enern,emeini,xstaten,
          eei,enerini,alcon,nalcon,set,nset,istartset,iendset,
          ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	  nelemload,&null,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,
	  ne,thicke,shcon,nshcon,
	  sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
          &mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
	  islavsurf,ielprop,prop,energyini,energy,&kscale,iponoel,
          inoel,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
	  itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
	  islavelinv,autloc,irowtloc,jqtloc,&nboun2,
	  ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
	  labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
	  &intscheme,physcon);
  
  SFREE(fn);SFREE(inum);SFREE(v);
  
  iout=1;
  
  NNEW(ad,double,neq[1]);
  NNEW(au,double,nzs[1]);
  
  mafillsmmain(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xbounold,nboun,
	       ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
	       &null,nelemload,sideload,xloadact,&null,xbodyact,ipobody,
	       &null,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,
	       &nmethodact,ikmpc,ilmpc,ikboun,ilboun,
	       elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
	       ielmat,ielorien,norien,orab,ntmat_,
	       t0,t1act,&ithermalact,prestr,iprestr,vold,iperturb,sti,
	       nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
	       xstiff,npmat_,&dtime,matname,mi,
	       ncmat_,massact,&stiffness,&buckling,&rhsi,&intscheme,
	       physcon,shcon,nshcon,alcon,nalcon,ttime,&time,istep,&iinc,
	       &coriolis,ibody,xloadold,&reltime,veold,springarea,nstate_,
	       xstateini,xstate,thicke,integerglob,doubleglob,
	       tieset,istartset,iendset,ialset,ntie,&nasym,pslavsurf,
	       pmastsurf,&mortar,clearini,ielprop,prop,&ne0,fnext,&kscale,
	       iponoel,inoel,network,ntrans,inotr,trab,smscale,&mscalmethod,
	       set,nset,islavelinv,autloc,irowtloc,jqtloc,&mortartrafoflag);
  
  if(nmethodact==0){
      
    /* error occurred in mafill: storing the geometry in frd format */
      
    ++*kode;
      
    ptime=*ttime+time;
    frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,&nmethodact,
	kode,filab,een,t1,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	nstate_,istep,&iinc,&ithermalact,qfn,&mode,&noddiam,trab,inotr,
	ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	thicke,jobnamec,output,qfx,cdn,&mortar,cdnr,cdni,nmat,ielprop,
	prop,sti);
      
    FORTRAN(stop,());
      
  }
  
  for(k=0;k<neq[1];++k){
    b[k]=fext[k]-f[k];
  }
  
  if(*isolver==0){
#ifdef SPOOLES
    spooles(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],
	    &symmetryflag,&inputformat,&nzs[2]);
#else
    printf(" *ERROR in electromagnetics: the SPOOLES library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if((*isolver==2)||(*isolver==3)){
    preiter(ad,&au,b,&icol,&irow,&neq[1],&nzs[1],isolver,iperturb);
  }
  else if(*isolver==4){
#ifdef SGI
    token=1;
    sgi_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],token);
#else
    printf(" *ERROR in electromagnetics: the SGI library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==5){
#ifdef TAUCS
    tau(ad,&au,adb,aub,&sigma,b,icol,&irow,&neq[1],&nzs[1]);
#else
    printf(" *ERROR in electromagnetics: the TAUCS library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==7){
#ifdef PARDISO
    pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],
		 &symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
#else
    printf(" *ERROR in electromagnetics: the PARDISO library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  else if(*isolver==8){
#ifdef PASTIX
    pastix_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],
		&symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
#else
    printf(" *ERROR in electromagnetics: the PASTIX library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }

  SFREE(au);SFREE(ad);

  NNEW(v,double,mt**nk);
  //  memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
  
  NNEW(fn,double,mt**nk);
  
  NNEW(inum,ITG,*nk);
  results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	  ielorien,norien,orab,ntmat_,t0,t1act,&ithermalact,
	  prestr,iprestr,filab,eme,emn,een,iperturb,
	  f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	  ndirboun,xbounact,nboun,ipompc,
	  nodempc,coefmpc,labmpc,nmpc,&nmethodact,cam,&neq[1],veold,accold,
	  &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	  xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	  &icmd,ncmat_,nstate_,sti,vini,ikboun,ilboun,ener,enern,
	  emeini,xstaten,eei,enerini,alcon,nalcon,set,nset,istartset,
	  iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
	  fmpc,nelemload,&null,ikmpc,ilmpc,istep,&iinc,springarea,
	  &reltime,ne,thicke,shcon,nshcon,
	  sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
          &mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
	  islavsurf,ielprop,prop,energyini,energy,&kscale,iponoel,
          inoel,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
	  itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
	  islavelinv,autloc,irowtloc,jqtloc,&nboun2,
	  ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
	  labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
	  &intscheme,physcon);
  
  //  memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);
  
  /* reactivating the non-shell elements (for mesh output purposes) 
     deactivating the initial temperature for the non-shell nodes */
  
  for(i=0;i<*ne;i++){
    if(strcmp1(&lakon[8*i+6],"L")!=0){
      ipkon[i]=-ipkon[i]-2;
    }else if(ipkon[i]!=-1){

      /* copy shell results */

      indexe=ipkon[i];
      if(strcmp1(&lakon[8*i+3],"6")==0){nope=6;}
      else if(strcmp1(&lakon[8*i+3],"8")==0){nope=8;}
      else if(strcmp1(&lakon[8*i+3],"1")==0){nope=15;}
      else{nope=20;}
      for(j=0;j<nope;j++){
	node=kon[indexe+j];
	vold[mt*(node-1)]=v[mt*(node-1)];
      }
    }
  }

  /* deactivating the output of temperatures */

  if(strcmp1(&filab[87],"NT  ")==0){
    ntflag=1;
    strcpy1(&filab[87],"    ",4);
  }else{ntflag=0;}

  /* createinum is called in order to store the nodes and elements
     of the complete structure, not only of the coil */

  FORTRAN(createinum,(ipkon,inum,kon,lakon,nk,ne,&cflag[0],nelemload,
		      nload,nodeboun,nboun,ndirboun,ithermal,co,vold,mi,ielmat,
                      ielprop,prop));

  ++*kode;
  if(*mcs!=0){
    ptime=*ttime+time;
    frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,&nmethodact,kode,filab,een,
	   t1act,fn,&ptime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
	   nstate_,istep,&iinc,iperturb,ener,mi,output,&ithermalact,qfn,
	   ialset,istartset,iendset,trab,inotr,ntrans,orab,ielorien,
	   norien,stx,veold,&noddiam,set,nset,emn,thicke,jobnamec,ne,
	   cdn,&mortar,nmat,qfx,ielprop,prop);
  }else{
      
    ptime=*ttime+time;
    frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,&nmethodact,
	kode,filab,een,t1act,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	nstate_,istep,&iinc,&ithermalact,qfn,&mode,&noddiam,trab,inotr,
	ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	thicke,jobnamec,output,qfx,cdn,&mortar,cdnr,cdni,nmat,ielprop,
	prop,sti);
      
  }
  SFREE(inum);SFREE(v);SFREE(fn);

  /* reactivating the temperature output, if previously deactivated */

  if(ntflag==1){
    strcpy1(&filab[87],"NT  ",4);
  }

  NNEW(inomat,ITG,*nk);

  /* calculating the magnetic intensity caused by the current */

  FORTRAN(assigndomtonodes,(ne,lakon,ipkon,kon,ielmat,inomat,
			    elcon,ncmat_,ntmat_,mi,&ne2));

  NNEW(h0ref,double,3**nk);
  NNEW(h0,double,3**nk);

  biosav(ipkon,kon,lakon,ne,co,qfx,h0ref,mi,inomat,nk);

  if(*ithermal<=1){SFREE(qfx);SFREE(t0);}
  if(strcmp1(&filab[3567],"ECD ")==0)SFREE(qfn);
  
  /* deactivating the shell elements */
  
  for(i=0;i<*ne;i++){
    if(strcmp1(&lakon[8*i+6],"L")==0){
      ipkon[i]=-ipkon[i]-2;
    }
  }
  
  /**************************************************************/
  /* creating connecting MPC's between the domains              */
  /**************************************************************/

  /* creating contact ties between the domains */

  if(*istep==1){
      
    NNEW(nodface,ITG,5*6**ne);
    NNEW(ipoface,ITG,*nk);
      
    RENEW(set,char,81*(*nset+3));
    RENEW(istartset,ITG,*nset+3);
    RENEW(iendset,ITG,*nset+3);
    RENEW(ialset,ITG,*nalset+6**ne);
    RENEW(tieset,char,243*(*ntie+5));
    RENEW(tietol,double,4*(*ntie+5));
      
    FORTRAN(createtiedsurfs,(nodface,ipoface,set,istartset,
			     iendset,ialset,tieset,inomat,ne,ipkon,lakon,kon,
			     ntie,tietol,nalset,nk,nset,iactive));
      
    SFREE(nodface);SFREE(ipoface);
    RENEW(set,char,81**nset);
    RENEW(istartset,ITG,*nset);
    RENEW(iendset,ITG,*nset);
    RENEW(ialset,ITG,*nalset);
    RENEW(tieset,char,243**ntie);
    RENEW(tietol,double,4**ntie);
      
    /* tied contact constraints: generate appropriate MPC's */
      
    tiedcontact(ntie,tieset,nset,set,istartset,iendset,ialset,
		lakon,ipkon,kon,tietol,nmpc, &mpcfree,&memmpc_,
		&ipompc,&labmpc,&ikmpc,&ilmpc,&fmpc,&nodempc,&coefmpc,
		ithermal,co,vold,&icfd,nmpc_,mi,nk,istep,ikboun,nboun,
		kind1,kind2,jobnamef);

    /* mapping h0ref from the phi domain onto the border of
       the A and A-V domains */

    FORTRAN(calch0interface,(nmpc,ipompc,nodempc,coefmpc,h0ref));

  }
  
  /**************************************************************/
  /* creating the A.n MPC                                       */
  /**************************************************************/

  /* identifying the interfaces between the A and A-V domains
     and the phi-domain */

  FORTRAN(generateeminterfaces,(istartset,iendset,
				ialset,iactive,ipkon,lakon,kon,ikmpc,nmpc,
				&maxfaces));
  
  for(i=1;i<3;i++){
    imast=iactive[i];
    if(imast==0) continue;

    /* determining the normals on the face */

    NNEW(imastnode,ITG,8*maxfaces);
    NNEW(xmastnor,double,3*8*maxfaces);

    FORTRAN(normalsoninterface,(istartset,iendset,
				ialset,&imast,ipkon,kon,lakon,imastnode,
				&nmastnode,xmastnor,co));

    /* enlarging the fields for MPC's */

    *nmpc_=*nmpc_+nmastnode;
    RENEW(ipompc,ITG,*nmpc_);
    RENEW(labmpc,char,20**nmpc_+1);
    RENEW(ikmpc,ITG,*nmpc_);
    RENEW(ilmpc,ITG,*nmpc_);
    RENEW(fmpc,double,*nmpc_);
      
    /* determining the maximum number of terms;
       expanding nodempc and coefmpc to accommodate
       those terms */
      
    neqterms=3*nmastnode;
    index=memmpc_;
    (memmpc_)+=neqterms;
    RENEW(nodempc,ITG,3*memmpc_);
    RENEW(coefmpc,double,memmpc_);
    for(k=index;k<memmpc_;k++){
      nodempc[3*k-1]=k+1;
    }
    nodempc[3*memmpc_-1]=0;

    /* creating the A.n MPC's */

    FORTRAN(createinterfacempcs,(imastnode,xmastnor,&nmastnode,
				 ikmpc,ilmpc,nmpc,ipompc,nodempc,coefmpc,labmpc,&mpcfree,
				 ikboun,nboun));

    SFREE(imastnode);SFREE(xmastnor);
  }
  
  /* determining the new matrix structure */
      
  remastructem(ipompc,&coefmpc,&nodempc,nmpc,
	       &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,ikboun,ilboun,
	       labmpc,nk,&memmpc_,&icascade,&maxlenmpc,
	       kon,ipkon,lakon,ne,nactdof,icol,jq,&irow,isolver,
	       neq,nzs,nmethod,&f,&fext,&b,&aux2,&fini,&fextini,
	       &adb,&aub,ithermal,iperturb,mass,mi,ielmat,elcon,
	       ncmat_,ntmat_,inomat,network);

  /**************************************************************/
  /* starting the loop over the increments                      */
  /**************************************************************/
  
  /* saving the distributed loads (volume heating will be
     added because of Joule heating) */

  //  if((*ithermal==3)||(*nmethod==2)){
  if((*nmethod==2)||(*nmethod==4)){
    nloadref=*nload;
    NNEW(nelemloadref,ITG,2**nload);
    if(*nam>0) NNEW(iamloadref,ITG,2**nload);
    NNEW(sideloadref,char,20**nload);
      
    memcpy(&nelemloadref[0],&nelemload[0],sizeof(ITG)*2**nload);
    if(*nam>0) memcpy(&iamloadref[0],&iamload[0],sizeof(ITG)*2**nload);
    memcpy(&sideloadref[0],&sideload[0],sizeof(char)*20**nload);
      
    /* generating new fields; ne2 is the number of elements
       in domain 2 = A,V-domain (the only domain with heating) */
      
    (*nload_)+=ne2;
    RENEW(nelemload,ITG,2**nload_);
    if(*nam>0) RENEW(iamload,ITG,2**nload_);
    RENEW(xloadact,double,2**nload_);
    RENEW(sideload,char,20**nload_);
  }

  if((*ithermal==1)||(*ithermal>=3)){
    NNEW(t1ini,double,*nk);
    NNEW(t1act,double,*nk);
    for(k=0;k<*nk;++k){t1act[k]=t1old[k];}
  }
  
  while(1.-theta>1.e-6){
      
    if(icutb==0){
	  
      /* previous increment converged: update the initial values */
	  
      iinc++;
      jprint++;
	  
      /* vold is copied into vini */
	  
      memcpy(&vini[0],&vold[0],sizeof(double)*mt**nk);
	  
      for(k=0;k<*nboun;++k){xbounini[k]=xbounact[k];}
      if((*ithermal==1)||(*ithermal>=3)){
	for(k=0;k<*nk;++k){t1ini[k]=t1act[k];}
      }
      for(k=0;k<neq[1];++k){
	fini[k]=f[k];
      }
      /*	  if((*nmethod==2)||(*nmethod==4)){
		  for(k=0;k<mt**nk;++k){
		  veini[k]=veold[k];
		  }
		  for(k=0;k<neq[1];++k){
		  fextini[k]=fext[k];
		  }
		  }*/
    }
      
    /* check for max. # of increments */
      
    if(iinc>*jmax){
      printf(" *ERROR: max. # of increments reached\n\n");
      FORTRAN(stop,());
    }

    if(*nmethod!=2){
      printf(" increment %" ITGFORMAT " attempt %" ITGFORMAT " \n",iinc,
	     icutb+1);
      printf(" increment size= %e\n",dtheta**tper);
      printf(" sum of previous increments=%e\n",theta**tper);
      printf(" actual step time=%e\n",(theta+dtheta)**tper);
      printf(" actual total time=%e\n\n",*ttime+(theta+dtheta)**tper);
      
      printf(" iteration 1\n\n");
    }
      
    qamold[0]=qam[0];
    qamold[1]=qam[1];
      
    /* determining the actual loads at the end of the new increment*/
      
    reltime=theta+dtheta;
    time=reltime**tper;
    dtime=dtheta**tper;

    /* restoring the distributed loading before adding the
       Joule heating */

    //    if((*ithermal==3)||(*nmethod==2)){
    if((*nmethod==2)||(*nmethod==4)){
      *nload=nloadref;
      DMEMSET(nelemload,0,2**nload_,0);
      memcpy(&nelemload[0],&nelemloadref[0],sizeof(ITG)*2**nload);
      if(*nam>0){
	DMEMSET(iamload,0,2**nload_,0);
	memcpy(&iamload[0],&iamloadref[0],sizeof(ITG)*2**nload);
      }
      DMEMSET(xloadact,0,2**nload_,0.);
      DMEMSET(sideload,0,'\0',0.);memcpy(&sideload[0],&sideloadref[0],sizeof(char)*20**nload);
    }
      
    /* determining the actual loading */

    //      for(i=0;i<3**nk;i++){h0[i]/=h0scale;}
    FORTRAN(tempload_em,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,xload,
			 xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,
			 xbodyact,
			 t1old,t1,t1act,iamt1,nk,amta,
			 namta,nam,ampli,&time,&reltime,ttime,&dtime,ithermal,
			 nmethod,
			 xbounold,xboun,xbounact,iamboun,nboun,
			 nodeboun,ndirboun,nodeforc,ndirforc,istep,&iinc,
			 co,vold,itg,&ntg,amname,ikboun,ilboun,nelemload,
			 sideload,mi,
			 ntrans,trab,inotr,veold,integerglob,doubleglob,tieset,
			 istartset,
			 iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,
			 coefmpc,
			 &h0scale,inomat,ipobody,iponoel,inoel,ipkon,kon,lakon,
			 ielprop,prop,ielmat,shcon,nshcon,rhcon,nrhcon,ntmat_,
			 cocon,
			 ncocon,set,nset));
    for(i=0;i<3**nk;i++){h0[i]=h0ref[i]*h0scale;}

    for(i=0;i<3;i++){cam[i]=0.;}for(i=3;i<5;i++){cam[i]=0.5;}
    if(*ithermal>1){
      radflowload(itg,ieg,&ntg,&ntr,adrad,aurad,bcr,ipivr,
		  ac,bc,nload,sideload,nelemload,xloadact,lakon,ipiv,ntmat_,
		  vold,shcon,nshcon,ipkon,kon,co,
		  kontri,&ntri,nloadtr,tarea,tenv,physcon,erad,&adview,&auview,
		  nflow,ikboun,xbounact,nboun,ithermal,&iinc,&iit,
		  cs,mcs,inocs,&ntrit,nk,fenv,istep,&dtime,ttime,&time,ilboun,
		  ikforc,ilforc,xforcact,nforc,cam,ielmat,&nteq,prop,ielprop,
		  nactdog,nacteq,nodeboun,ndirboun,network,
		  rhcon,nrhcon,ipobody,ibody,xbodyact,nbody,iviewfile,jobnamef,
		  ctrl,xloadold,&reltime,nmethod,set,mi,istartset,iendset,
		  ialset,nset,
		  ineighe,nmpc,nodempc,ipompc,coefmpc,labmpc,&iemchange,nam,
		  iamload,
		  jqrad,irowrad,&nzsrad,icolrad,ne,iaxial,qa,cocon,ncocon,
		  iponoel,
		  inoel,nprop,amname,namta,amta,iexpl);
    }
      
    /* prediction of the next solution (only for temperature)
       for the potentials the solution at the end of the previous
       increment is not changed (veold=0) */
      
    NNEW(v,double,mt**nk);
    
    prediction_em(uam,nmethod,&bet,&gam,&dtime,ithermal,nk,veold,v,
		  &iinc,&idiscon,vold,nactdof,mi);
      
    NNEW(fn,double,mt**nk);
      
    iout=-1;
    iperturb_sav[0]=iperturb[0];
    iperturb_sav[1]=iperturb[1];
      
    /* first iteration in first increment: heat tangent 
       a dummy field stx is created in order not to change 
       the values of sti from the previous increment
       (if this is not done the electric field would be
        zero due to v=vold; this would lead to a zero
        heat power in jouleheating.f ) */
      
    NNEW(inum,ITG,*nk);
    NNEW(stx,double,6**ne*mi[0]);
    resultsinduction(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
		     elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		     ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
		     prestr,iprestr,filab,eme,emn,een,iperturb,
		     f,fn,nactdof,&iout,qa,vold,b,nodeboun,
		     ndirboun,xbounact,nboun,ipompc,
		     nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,
		     accold,
		     &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		     xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
		     &icmd,ncmat_,nstate_,stx,vini,ikboun,ilboun,ener,enern,
		     emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
		     iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,
		     ntrans,
		     fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
		     &reltime,ne,xforc,nforc,thicke,shcon,nshcon,
		     sideload,xloadact,xloadold,&icfd,inomat,h0,islavnode,
		     nslavnode,ntie,ielprop,prop,iactive,energyini,energy,
		     iponoel,inoel,orname,network,ipobody,xbodyact,ibody,nbody);
    SFREE(inum);SFREE(stx);
     
    /* the calculation of the electromagnetic fields is (quasi)linear,
       i.e. the solution of the equations is the fields;
       only the temperature calculation is nonlinear,
       i.e. the solution of the equations is a differential temperature */
 
    memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);
      
    iout=0;
      
    SFREE(fn);SFREE(v);
      
    /***************************************************************/
    /* iteration counter and start of the loop over the iterations */
    /***************************************************************/
      
    iit=1;
    icntrl=0;
    ctrl[0]=i0ref;ctrl[1]=irref;ctrl[3]=icref;
      
    while(icntrl==0){
	  
      if(iit!=1){
	      
	printf(" iteration %" ITGFORMAT "\n\n",iit);

	/* restoring the distributed loading before adding the
	   Joule heating */

	//	if((*ithermal==3)||(*nmethod==2)){
	if((*nmethod==2)||(*nmethod==4)){
	  *nload=nloadref;
	  DMEMSET(nelemload,0,2**nload_,0);
	  memcpy(&nelemload[0],&nelemloadref[0],sizeof(ITG)*2**nload);
	  if(*nam>0){
	    DMEMSET(iamload,0,2**nload_,0);
	    memcpy(&iamload[0],&iamloadref[0],sizeof(ITG)*2**nload);
	  }
	  DMEMSET(xloadact,0,2**nload_,0.);
	  DMEMSET(sideload,0,20**nload_,'\0');memcpy(&sideload[0],&sideloadref[0],sizeof(char)*20**nload);
	}
      
	/* determining the actual loading */
	      
	//	      for(i=0;i<3**nk;i++){h0[i]/=h0scale;}
	FORTRAN(tempload_em,(xforcold,xforc,xforcact,iamforc,nforc,
			     xloadold,xload,
			     xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,
			     xbodyact,
			     t1old,t1,t1act,iamt1,nk,amta,
			     namta,nam,ampli,&time,&reltime,ttime,&dtime,
			     ithermal,nmethod,
			     xbounold,xboun,xbounact,iamboun,nboun,
			     nodeboun,ndirboun,nodeforc,ndirforc,istep,&iinc,
			     co,vold,itg,&ntg,amname,ikboun,ilboun,nelemload,
			     sideload,mi,
			     ntrans,trab,inotr,veold,integerglob,doubleglob,
			     tieset,istartset,
			     iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,
			     nodempc,coefmpc,
			     &h0scale,inomat,ipobody,iponoel,inoel,ipkon,kon,
			     lakon,
			     ielprop,prop,ielmat,shcon,nshcon,rhcon,nrhcon,
			     ntmat_,cocon,
			     ncocon,set,nset));
	for(i=0;i<3**nk;i++){h0[i]=h0ref[i]*h0scale;}

	for(i=0;i<3;i++){cam[i]=0.;}for(i=3;i<5;i++){cam[i]=0.5;}
	if(*ithermal>1){
	  radflowload(itg,ieg,&ntg,&ntr,adrad,aurad,
		      bcr,ipivr,ac,bc,nload,sideload,nelemload,xloadact,lakon,
		      ipiv,
		      ntmat_,vold,shcon,nshcon,ipkon,kon,co,
		      kontri,&ntri,nloadtr,tarea,tenv,physcon,erad,&adview,
		      &auview,
		      nflow,ikboun,xbounact,nboun,ithermal,&iinc,&iit,
		      cs,mcs,inocs,&ntrit,nk,fenv,istep,&dtime,ttime,&time,
		      ilboun,
		      ikforc,ilforc,xforcact,nforc,cam,ielmat,&nteq,prop,
		      ielprop,
		      nactdog,nacteq,nodeboun,ndirboun,network,
		      rhcon,nrhcon,ipobody,ibody,xbodyact,nbody,iviewfile,
		      jobnamef,
		      ctrl,xloadold,&reltime,nmethod,set,mi,istartset,iendset,
		      ialset,
		      nset,ineighe,nmpc,nodempc,ipompc,coefmpc,labmpc,
		      &iemchange,nam,
		      iamload,jqrad,irowrad,&nzsrad,icolrad,ne,iaxial,qa,cocon,
		      ncocon,iponoel,inoel,nprop,amname,namta,amta,iexpl);
	}
	      
      }
	  
      /* add Joule heating  */

      if(*nmethod==4){
	NNEW(idefload,ITG,*nload_);
	DMEMSET(idefload,0,*nload_,1);
	FORTRAN(jouleheating,(ipkon,lakon,kon,co,elcon,nelcon,
			      mi,ne,sti,ielmat,nelemload,sideload,xloadact,
			      nload,nload_,
			      iamload,nam,idefload,ncmat_,ntmat_,
			      alcon,nalcon,ithermal,vold,t1,nmethod));
	SFREE(idefload);
      }

      if(*ithermal==3){
	for(k=0;k<*nk;++k){t1act[k]=vold[mt*k];}
      }
	      
      /* calculating the local stiffness matrix and external loading */
	  
      NNEW(ad,double,neq[1]);
      NNEW(au,double,nzs[1]);
	  
      FORTRAN(mafillem,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
			xbounact,nboun,
			ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
			nforc,nelemload,sideload,xloadact,nload,xbodyact,
			ipobody,
			nbody,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,
			nmethod,ikmpc,ilmpc,ikboun,ilboun,
			elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
			ielmat,ielorien,norien,orab,ntmat_,
			t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
			nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
			xstiff,npmat_,&dtime,matname,mi,
			ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
			physcon,shcon,nshcon,cocon,ncocon,ttime,&time,istep,
			&iinc,
			&coriolis,ibody,xloadold,&reltime,veold,springarea,
			nstate_,
			xstateini,xstate,thicke,integerglob,doubleglob,
			tieset,istartset,iendset,ialset,ntie,&nasym,iactive,h0,
			pslavsurf,pmastsurf,&mortar,clearini,ielprop,prop,
			iponoel,inoel,network));
	      
      iperturb[0]=iperturb_sav[0];
      iperturb[1]=iperturb_sav[1];

      /* calculating the residual (f is only for the temperature
	 nonzero) */

      calcresidual_em(nmethod,neq,b,fext,f,iexpl,nactdof,aux1,aux2,vold,
		      vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,fextini,
		      fini,
		      islavnode,nslavnode,&mortar,ntie,f_cm,f_cs,mi,
		      nzs,&nasym,ithermal);
	  
      if(*nmethod==0){
	      
	/* error occurred in mafill: storing the geometry in frd format */
	      
	*nmethod=0;
	++*kode;
	NNEW(inum,ITG,*nk);for(k=0;k<*nk;k++) inum[k]=1;
	      
	ptime=*ttime+time;
	frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	    kode,filab,een,t1,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	    nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	    thicke,jobnamec,output,qfx,cdn,&mortar,cdnr,cdni,nmat,
	    ielprop,prop,sti);
	      
      }
	  
      if(*nmethod==2){

	/* frequency calculation 
	   expanding the matrices K and M */

	neqfreq=2*neq[0];
	nzsfreq=8*nzs[0]+neqfreq;
	NNEW(adfreq,double,neqfreq);
	NNEW(aufreq,double,nzsfreq);
	NNEW(irowfreq,ITG,nzsfreq);
	NNEW(iaux,ITG,nzsfreq);
	NNEW(icolfreq,ITG,neqfreq);
	NNEW(jqfreq,ITG,neqfreq+1);
	NNEW(bfreq,double,neqfreq);

	om=alpha[0];
	FORTRAN(mafillfreq_em,(ad,au,adb,aub,irow,jq,&neq[0],
			       adfreq,aufreq,irowfreq,iaux,jqfreq,
			       icolfreq,&neqfreq,&nzsfreq,&om,&symmetryflag,
			       &inputformat,b,bfreq));

	SFREE(iaux);

	RENEW(ad,double,neqfreq);
	memcpy(ad,adfreq,sizeof(double)*neqfreq);
	SFREE(adfreq);

	RENEW(au,double,nzsfreq);
	memcpy(au,aufreq,sizeof(double)*nzsfreq);
	SFREE(aufreq);

	RENEW(irow,ITG,nzsfreq);
	memcpy(irow,irowfreq,sizeof(ITG)*nzsfreq);
	SFREE(irowfreq);

	RENEW(icol,ITG,neqfreq);
	memcpy(icol,icolfreq,sizeof(ITG)*neqfreq);
	SFREE(icolfreq);

	RENEW(jq,ITG,neqfreq+1);
	memcpy(jq,jqfreq,sizeof(ITG)*(neqfreq+1));
	SFREE(jqfreq);

	RENEW(b,double,neqfreq);
	memcpy(b,bfreq,sizeof(double)*neqfreq);
	SFREE(bfreq);
	      
	neq[0]=neqfreq;neq[1]=neqfreq;neq[2]=neqfreq;
	nzs[0]=nzsfreq;nzs[1]=nzsfreq;nzs[2]=nzsfreq;

	/* summary */

	printf(" Expanding the system of equations for harmonic analysis:\n");
	printf(" number of equations\n");
	printf(" %" ITGFORMAT "\n",neqfreq);
	printf(" number of nonzero off-diagonal matrix elements\n");
	printf(" %" ITGFORMAT "\n",nzsfreq);
	printf("\n");

      }else if(*nmethod==4){
	//	  }else if((*nmethod==2)||(*nmethod==4)){
	      
	/* transient calculation */

	/* electromagnetic part */
	      
	if(*ithermal!=2){
	  for(k=0;k<neq[0];++k){
	    ad[k]=adb[k]/dtime+ad[k];
	  }
	  for(k=0;k<nzs[0];++k){
	    au[k]=aub[k]/dtime+au[k];
	  }
		  
	  /* upper triangle of asymmetric matrix */
		  
	  if(nasym>0){
	    for(k=nzs[2];k<nzs[2]+nzs[0];++k){
	      au[k]=aub[k]/dtime+au[k];
	    }
	  }
	}
	      
	/* thermal part */
	      
	if(*ithermal>1){
	  for(k=neq[0];k<neq[1];++k){
	    ad[k]=adb[k]/dtime+ad[k];
	  }
	  for(k=nzs[0];k<nzs[1];++k){
	    au[k]=aub[k]/dtime+au[k];
	  }
		  
	  /* upper triangle of asymmetric matrix */
		  
	  if(nasym>0){
	    for(k=nzs[2]+nzs[0];k<nzs[2]+nzs[1];++k){
	      au[k]=aub[k]/dtime+au[k];
	    }
	  }
	}
      }

      NNEW(adaux,double,neq[1]);
      FORTRAN(preconditioning,(ad,au,b,&neq[1],irow,jq,adaux));

	  
      if(*isolver==0){
#ifdef SPOOLES
	if(*ithermal<2){
	  spooles(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
		  &symmetryflag,&inputformat,&nzs[2]);
	}else{
	  spooles(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],
		  &symmetryflag,&inputformat,&nzs[2]);
	}
#else
	printf(" *ERROR in electromagnetics: the SPOOLES library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if((*isolver==2)||(*isolver==3)){
	if(nasym>0){
	  printf(" *ERROR in electromagnetics: the iterative solver cannot be used for asymmetric matrices\n\n");
	  FORTRAN(stop,());
	}
	preiter(ad,&au,b,&icol,&irow,&neq[1],&nzs[1],isolver,iperturb);
      }
      else if(*isolver==4){
#ifdef SGI
	if(nasym>0){
	  printf(" *ERROR in electromagnetics: the SGI solver cannot be used for asymmetric matrices\n\n");
	  FORTRAN(stop,());
	}
	token=1;
	if(*ithermal<2){
	  sgi_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],token);
	}else{
	  sgi_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],token);
	}
#else
	printf(" *ERROR in electromagnetics: the SGI library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==5){
	if(nasym>0){
	  printf(" *ERROR in electromagnetics: the TAUCS solver cannot be used for asymmetric matrices\n\n");
	  FORTRAN(stop,());
	}
#ifdef TAUCS
	tau(ad,&au,adb,aub,&sigma,b,icol,&irow,&neq[1],&nzs[1]);
#else
	printf(" *ERROR in electromagnetics: the TAUCS library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	if(*ithermal<2){
	  pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
		       &symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
	}else{
	  pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],
		       &symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
	}
#else
	printf(" *ERROR in electromagnetics: the PARDISO library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==8){
#ifdef PASTIX
	if(*ithermal<2){
	  pastix_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
		      &symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
	}else{
	  pastix_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],
		      &symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
	}
#else
	printf(" *ERROR in electromagnetics: the PASTIX library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }

      SFREE(au);SFREE(ad);

      for(i=0;i<neq[1];i++){b[i]*=adaux[i];}
      SFREE(adaux);
	  
      if(*nmethod==2){
	neq[0]=neqfreq/2;neq[1]=neq[0];
	nzs[0]=(nzsfreq-neqfreq)/8;nzs[1]=nzs[0];
	RENEW(irow,ITG,nzs[0]);
	RENEW(icol,ITG,neq[0]);
	RENEW(jq,ITG,neq[0]+1);
	break;
      }

      /* calculating the electromagnetic fields and temperatures 
	 only the temperature calculation is differential */
	  
      NNEW(v,double,mt**nk);
      memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
	  
      NNEW(fn,double,mt**nk);
	  
      NNEW(inum,ITG,*nk);
      resultsinduction(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
		       elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		       ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
		       prestr,iprestr,filab,eme,emn,een,iperturb,
		       f,fn,nactdof,&iout,qa,vold,b,nodeboun,
		       ndirboun,xbounact,nboun,ipompc,
		       nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,
		       accold,
		       &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,
		       nplkcon,
		       xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
		       &icmd,ncmat_,nstate_,sti,vini,ikboun,ilboun,ener,enern,
		       emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,
		       istartset,
		       iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,
		       ntrans,
		       fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
		       &reltime,ne,xforc,nforc,thicke,shcon,nshcon,
		       sideload,xloadact,xloadold,&icfd,inomat,h0,islavnode,
		       nslavnode,ntie,ielprop,prop,iactive,energyini,energy,
		       iponoel,inoel,orname,network,ipobody,xbodyact,ibody,
		       nbody);
      SFREE(inum);
	  
      if(*ithermal!=2){
	if(cam[0]>uam[0]){uam[0]=cam[0];}      
	if(qau<1.e-10){
	  if(qa[0]>ea*qam[0]){qam[0]=(qamold[0]*jnz+qa[0])/(jnz+1);}
	  else {qam[0]=qamold[0];}
	}
      }
      if(*ithermal>1){
	if(cam[1]>uam[1]){uam[1]=cam[1];}      
	if(qau<1.e-10){
	  if(qa[1]>ea*qam[1]){qam[1]=(qamold[1]*jnz+qa[1])/(jnz+1);}
	  else {qam[1]=qamold[1];}
	}
      }
	  
      memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);
	  
      SFREE(v);SFREE(fn);
	  
      /* calculating the residual */
	  
      calcresidual_em(nmethod,neq,b,fext,f,iexpl,nactdof,aux1,aux2,vold,
		      vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,fextini,fini,
		      islavnode,nslavnode,&mortar,ntie,f_cm,f_cs,mi,
		      nzs,&nasym,ithermal);
	  
      /* calculating the maximum residual (only thermal part)*/
	  
      for(k=0;k<2;++k){
	ram2[k]=ram1[k];
	ram1[k]=ram[k];
	ram[k]=0.;
      }

      if(*ithermal!=2) ram[0]=0.;

      if(*ithermal>1){
	for(k=neq[0];k<neq[1];++k){
	  err=fabs(b[k]);
	  if(err>ram[1]){ram[1]=err;ram[3]=k+0.5;}
	}
      }
	  
      /* printing residuals */
	  
      if(*ithermal>1){
	if(ram[1]<1.e-6) ram[1]=0.;      
	printf(" average flux= %f\n",qa[1]);
	printf(" time avg. flux= %f\n",qam[1]);
	if((ITG)((double)nactdofinv[(ITG)ram[3]]/mt)+1==0){
	  printf(" largest residual flux= %f\n",
		 ram[1]);
	}else{
	  inode=(ITG)((double)nactdofinv[(ITG)ram[3]]/mt)+1;
	  idir=nactdofinv[(ITG)ram[3]]-mt*(inode-1);
	  printf(" largest residual flux= %f in node %" ITGFORMAT " and dof %" ITGFORMAT "\n",ram[1],inode,idir);
	}
	printf(" largest increment of temp= %e\n",uam[1]);
	if((ITG)cam[4]==0){
	  printf(" largest correction to temp= %e\n\n",
		 cam[1]);
	}else{
	  inode=(ITG)((double)nactdofinv[(ITG)cam[4]]/mt)+1;
	  idir=nactdofinv[(ITG)cam[4]]-mt*(inode-1);
	  printf(" largest correction to temp= %e in node %" ITGFORMAT " and dof %" ITGFORMAT "\n\n",cam[1],inode,idir);
	}
      }
	  
      /* athermal electromagnetic calculations are linear:
	 set iit=2 to force convergence */

      if(*ithermal<=1) iit=2;
	  
      // MPADD: need for fake energy values!
      double energy[4] = {0, 0, 0, 0};
      double allwk     = 0.0;
      double energyref = 0.0;
      double emax, enres,enetoll, reswk, dampwk, allwkini;

      checkconvergence(co,nk,kon,ipkon,lakon,ne,stn,nmethod, 
		       kode,filab,een,t1act,&time,epn,ielmat,matname,enern, 
		       xstaten,nstate_,istep,&iinc,iperturb,ener,mi,output,
		       ithermal,qfn,&mode,&noddiam,trab,inotr,ntrans,orab,
		       ielorien,norien,description,sti,&icutb,&iit,&dtime,qa,
		       vold,qam,ram1,ram2,ram,cam,uam,&ntg,ttime,&icntrl,
		       &theta,&dtheta,veold,vini,idrct,tper,&istab,tmax, 
		       nactdof,b,tmin,ctrl,amta,namta,itpamp,&inext,&dthetaref,
		       &itp,&jprint,jout,&uncoupled,t1,&iitterm,nelemload,
		       nload,nodeboun,nboun,itg,ndirboun,&deltmx,&iflagact,
		       set,nset,istartset,iendset,ialset,emn,thicke,jobnamec,
		       &mortar,nmat,ielprop,prop,&ialeatoric,&kscale,
		       energy, &allwk, &energyref,&emax, &enres, &enetoll,        //MPADD
		       energyini, &allwkini ,&allwk, &reswk, &ne0, &ne0, &dampwk, //MPADD
		       &dampwk, energy);                                          //MPADD
    }
      
    /*********************************************************/
    /*   end of the iteration loop                          */
    /*********************************************************/
	  
    if(*nmethod==2) break;
      
    /* icutb=0 means that the iterations in the increment converged,
       icutb!=0 indicates that the increment has to be reiterated with
       another increment size (dtheta) */
      
    if(((qa[0]>ea*qam[0])||(qa[1]>ea*qam[1]))&&(icutb==0)){jnz++;}
    iit=0;
      
    if(icutb!=0){
      memcpy(&vold[0],&vini[0],sizeof(double)*mt**nk);
	  
      for(k=0;k<*nboun;++k){xbounact[k]=xbounini[k];}
      if((*ithermal==1)||(*ithermal>=3)){
	for(k=0;k<*nk;++k){t1act[k]=t1ini[k];}
      }
      for(k=0;k<neq[1];++k){
	f[k]=fini[k];
      }
	  
      qam[0]=qamold[0];
      qam[1]=qamold[1];
    }
      
    if((jout[0]==jprint)&&(icutb==0)){
	  
      jprint=0;
	  
      /* calculating the displacements and the stresses and storing */
      /* the results in frd format  */
	  
      NNEW(v,double,mt**nk);
      NNEW(fn,double,mt**nk);
      if(*ithermal>1) NNEW(qfn,double,3**nk);
      if((strcmp1(&filab[3741],"EMFE")==0)||
	 (strcmp1(&filab[3828],"EMFB")==0)) NNEW(stn,double,6**nk);
      NNEW(inum,ITG,*nk);
	  
      memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);

      iout=2;
      icmd=3;
	  
      resultsinduction(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
		       elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		       ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
		       prestr,iprestr,filab,eme,emn,een,iperturb,
		       f,fn,nactdof,&iout,qa,vold,b,nodeboun,
		       ndirboun,xbounact,nboun,ipompc,
		       nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,
		       accold,
		       &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,
		       nplkcon,
		       xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
		       &icmd,
		       ncmat_,nstate_,sti,vini,ikboun,ilboun,ener,enern,emeini,
		       xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
		       iendset,
		       ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
		       nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
		       &reltime,ne,xforc,nforc,thicke,shcon,nshcon,
		       sideload,xloadact,xloadold,&icfd,inomat,h0,islavnode,
		       nslavnode,ntie,ielprop,prop,iactive,energyini,energy,
		       iponoel,inoel,orname,network,ipobody,xbodyact,ibody,
		       nbody);
	  
      memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);
	  
      /* add Joule heating  */

      /*   NNEW(idefload,ITG,*nload_);
      DMEMSET(idefload,0,*nload_,1);
      FORTRAN(jouleheating,(ipkon,lakon,kon,co,elcon,nelcon,
			    mi,ne,sti,ielmat,nelemload,sideload,xloadact,
			    nload,nload_,
			    iamload,nam,idefload,ncmat_,ntmat_,
			    alcon,nalcon,ithermal,vold,t1,nmethod));
			    SFREE(idefload);*/

      /* printing the joule heating */
      
      /*  FORTRAN(printoutebhe,(set,nset,istartset,iendset,ialset,nprint,prlab,
			    prset,t1,ipkon,lakon,sti,ener,mi,ithermal,co,kon,
			    ttime,ne,vold,ielmat,thicke,&mortar,&time,ielprop,
			    prop,nelemload,nload,sideload,xloadact,rhcon,nrhcon,
			    ntmat_,ipobody,ibody,xbody,nbody,nmethod));*/
	  
      iout=0;
      icmd=0;
      //	  FORTRAN(networkinum,(ipkon,inum,kon,lakon,ne,itg,&ntg));
	  
      ++*kode;
      if(*mcs!=0){
	ptime=*ttime+time;
	frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,
	       t1act,fn,&ptime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
	       nstate_,istep,&iinc,iperturb,ener,mi,output,ithermal,qfn,
	       ialset,istartset,iendset,trab,inotr,ntrans,orab,ielorien,
	       norien,stx,veold,&noddiam,set,nset,emn,thicke,jobnamec,ne,
	       cdn,&mortar,nmat,qfx,ielprop,prop);
      }else{
	      
	ptime=*ttime+time;
	frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	    kode,filab,een,t1act,fn,&ptime,epn,ielmat,matname,
	    enern,xstaten,
	    nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	    thicke,jobnamec,output,qfx,cdn,&mortar,cdnr,cdni,nmat,
	    ielprop,prop,sti);
	      
      }
	  
      SFREE(v);SFREE(fn);SFREE(inum);
      if(*ithermal>1){SFREE(qfn);}
      if((strcmp1(&filab[3741],"EMFE")==0)||
	 (strcmp1(&filab[3828],"EMFB")==0)) SFREE(stn);
	  
    }
      
  }
  
  /*********************************************************/
  /*   end of the increment loop                          */
  /*********************************************************/
  
  if(jprint!=0){

    if(*nmethod==2){

      /* frequency step */

      NNEW(v,double,2*mt**nk);
      
      NNEW(fn,double,2*mt**nk);

      /* stx is used for the real and imaginary part
	 of -d(A+grad v)/dt at the integration points (first
         3 components) and for the real and imaginary part of
         B=grad x A= rot A (last 3 components) */
      
      NNEW(stx,double,2*6**ne*mi[0]);
      
      if((strcmp1(&filab[3741],"EMFE")==0)||
	 (strcmp1(&filab[3828],"EMFB")==0)){
	NNEW(stn,double,2*6**nk);
      }
      NNEW(inum,ITG,*nk);
      
      iout=1;
      icmd=3;

      for(k=0;k<2*neq[1];k+=neq[1]){

	if(k==0){
	  kkv=0;kk6=0;kkx=0;
	}else{
	  kkv=mt**nk;kk6=6**nk;kkx=6*mi[0]**ne;
	}
	  
	resultsinduction(co,nk,kon,ipkon,lakon,ne,&v[kkv],&stn[kk6],inum,elcon,
			 nelcon,
			 rhcon,nrhcon,alcon,nalcon,alzero,ielmat,ielorien,
			 norien,
			 orab,ntmat_,t0,t1act,ithermal,prestr,iprestr,filab,eme,
			 emn,een,iperturb,
			 f,&fn[kkv],nactdof,&iout,qa,vold,&b[k],nodeboun,
			 ndirboun,xbounact,nboun,ipompc,
			 nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,
			 accold,
			 &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,
			 nplkcon,
			 xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
			 &icmd,
			 ncmat_,nstate_,&stx[kkx],vini,ikboun,ilboun,ener,
			 enern,emeini,
			 xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
			 iendset,
			 ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
			 fmpc,
			 nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
			 &reltime,ne,xforc,nforc,thicke,shcon,nshcon,
			 sideload,xloadact,xloadold,&icfd,inomat,h0,islavnode,
			 nslavnode,ntie,ielprop,prop,iactive,energyini,energy,
			 iponoel,inoel,orname,network,ipobody,xbodyact,ibody,
			 nbody);
      }

      /* correct the electric field with contributions coming from 
         exp(i*omega*t):
      d(X*exp(i*omega*t))/dt=(dX/dt+X*i*omega)*exp(i*omega*t),
      where X=-A-grad(v) is complex */
      
      FORTRAN(correctem,(stx,stn,prlab,
			 nprint,ne,ipkon,lakon,elcon,ncmat_,ntmat_,nk,&om,
			 filab,mi,ielmat));
	  
      /* add Joule heating  */

      NNEW(idefload,ITG,*nload_);
      DMEMSET(idefload,0,*nload_,1);
      FORTRAN(jouleheating,(ipkon,lakon,kon,co,elcon,nelcon,
			    mi,ne,stx,ielmat,nelemload,sideload,xloadact,
			    nload,nload_,
			    iamload,nam,idefload,ncmat_,ntmat_,
			    alcon,nalcon,ithermal,vold,t1,nmethod));
      SFREE(idefload);

      SFREE(stx);SFREE(fn);

      /* printing the joule heating */
      
      FORTRAN(printoutebhe,(set,nset,istartset,iendset,ialset,nprint,prlab,
			    prset,t1,ipkon,lakon,stx,ener,mi,ithermal,co,kon,
			    ttime,ne,vold,ielmat,thicke,&mortar,&time,ielprop,
			    prop,nelemload,nload,sideload,xloadact,rhcon,nrhcon,
			    ntmat_,ipobody,ibody,xbody,nbody,nmethod));

      /* storing the real and imagninary part of the solution */
      
      ++*kode;
      ptime=*ttime+time;
      frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	  kode,filab,een,t1act,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	  nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	  mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	  thicke,jobnamec,output,qfx,cdn,&mortar,cdnr,cdni,nmat,
	  ielprop,prop,sti);

      SFREE(v);
      
      if((strcmp1(&filab[3741],"EMFE")==0)||
	 (strcmp1(&filab[3828],"EMFB")==0)){
	SFREE(stn);
      }
      
    }else{
      
      mode=-1;
     
      /* calculating the displacements and the stresses and storing  
	 the results in frd format */
      
      NNEW(v,double,mt**nk);
      NNEW(fn,double,mt**nk);
      if(*ithermal>1) NNEW(qfn,double,3**nk);
      if((strcmp1(&filab[3741],"EMFE")==0)||
	 (strcmp1(&filab[3828],"EMFB")==0)) NNEW(stn,double,6**nk);
      NNEW(inum,ITG,*nk);
	  
      memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
      iout=2;
      icmd=3;
	  
      resultsinduction(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
		       elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		       ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
		       prestr,iprestr,filab,eme,emn,een,iperturb,
		       f,fn,nactdof,&iout,qa,vold,&b[(mode+1)*neq[1]],nodeboun,
		       ndirboun,xbounact,nboun,ipompc,
		       nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,
		       accold,
		       &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,
		       nplkcon,
		       xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
		       &icmd,
		       ncmat_,nstate_,sti,vini,ikboun,ilboun,ener,enern,emeini,
		       xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
		       iendset,
		       ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
		       nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
		       &reltime,ne,xforc,nforc,thicke,shcon,nshcon,
		       sideload,xloadact,xloadold,&icfd,inomat,h0,islavnode,
		       nslavnode,ntie,ielprop,prop,iactive,energyini,energy,
		       iponoel,inoel,orname,network,ipobody,xbodyact,ibody,
		       nbody);
	  
      memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);
	  
      iout=0;
      icmd=0;
	  
      ++*kode;
      if(*mcs>0){
	ptime=*ttime+time;
	frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,
	       t1act,fn,&ptime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
	       nstate_,istep,&iinc,iperturb,ener,mi,output,ithermal,qfn,
	       ialset,istartset,iendset,trab,inotr,ntrans,orab,ielorien,
	       norien,stx,veold,&noddiam,set,nset,emn,thicke,jobnamec,ne,
	       cdn,&mortar,nmat,qfx,ielprop,prop);
      }else{
	      
	ptime=*ttime+time;
	frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	    kode,filab,een,t1act,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	    nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	    thicke,jobnamec,output,qfx,cdn,&mortar,cdnr,cdni,nmat,
	    ielprop,prop,sti);
	      
      }
	  
      SFREE(v);SFREE(fn);SFREE(inum);
      if(*ithermal>1){SFREE(qfn);}
      if((strcmp1(&filab[3741],"EMFE")==0)||
	 (strcmp1(&filab[3828],"EMFB")==0)) SFREE(stn);
    }
      
  }

  /* restoring the distributed loading  */

  //  if((*ithermal==3)||(*nmethod==2)){
  if((*nmethod==2)||(*nmethod==4)){
    *nload=nloadref;
    (*nload_)-=ne2;
    RENEW(nelemload,ITG,2**nload);memcpy(&nelemload[0],&nelemloadref[0],
					 sizeof(ITG)*2**nload);
    if(*nam>0){
      RENEW(iamload,ITG,2**nload);
      memcpy(&iamload[0],&iamloadref[0],sizeof(ITG)*2**nload);
    }
    RENEW(sideload,char,20**nload);memcpy(&sideload[0],&sideloadref[0],
					  sizeof(char)*20**nload);
      
    /* freeing the temporary fields */
      
    SFREE(nelemloadref);if(*nam>0){SFREE(iamloadref);};
    SFREE(sideloadref);
  }
  
  /* setting the velocity to zero at the end of a quasistatic or stationary
     step */
  
  if(abs(*nmethod)==1){
    for(k=0;k<mt**nk;++k){veold[k]=0.;}
  }
  
  /* updating the loading at the end of the step; 
     important in case the amplitude at the end of the step
     is not equal to one */
  
  for(k=0;k<*nboun;++k){
      
    /* thermal boundary conditions are updated only if the
       step was thermal or thermomechanical */
      
    if(ndirboun[k]==0){
      if(*ithermal<2) continue;
	  
      /* mechanical boundary conditions are updated only
	 if the step was not thermal or the node is a
	 network node */
	  
    }else if((ndirboun[k]>0)&&(ndirboun[k]<4)){
      node=nodeboun[k];
      FORTRAN(nident,(itg,&node,&ntg,&id));
      networknode=0;
      if(id>0){
	if(itg[id-1]==node) networknode=1;
      }
      if((*ithermal==2)&&(networknode==0)) continue;
    }
    xbounold[k]=xbounact[k];
  }
  for(k=0;k<*nforc;++k){xforcold[k]=xforcact[k];}
  for(k=0;k<2**nload;++k){xloadold[k]=xloadact[k];}
  for(k=0;k<7**nbody;k=k+7){xbodyold[k]=xbodyact[k];}
  if(*ithermal==1){
    for(k=0;k<*nk;++k){t1old[k]=t1act[k];}
    for(k=0;k<*nk;++k){vold[mt*k]=t1act[k];}
  }
  else if(*ithermal>1){
    for(k=0;k<*nk;++k){t1[k]=vold[mt*k];}
    if(*ithermal>=3){
      for(k=0;k<*nk;++k){t1old[k]=t1act[k];}
    }
  }
  
  qaold[0]=qa[0];
  qaold[1]=qa[1];
  
  SFREE(f);
  SFREE(b);
  SFREE(xbounact);SFREE(xforcact);SFREE(xloadact);SFREE(xbodyact);
  //  if(*nbody>0) SFREE(ipobody);
  SFREE(fext);SFREE(ampli);SFREE(xbounini);SFREE(xstiff);
  if((*ithermal==1)||(*ithermal>=3)){SFREE(t1act);SFREE(t1ini);}
  
  if(*ithermal>1){
    SFREE(itg);SFREE(ieg);SFREE(kontri);SFREE(nloadtr);
    SFREE(nactdog);SFREE(nacteq);SFREE(ineighe);
    SFREE(tarea);SFREE(tenv);SFREE(fenv);SFREE(qfx);
    SFREE(erad);SFREE(ac);SFREE(bc);SFREE(ipiv);
    SFREE(bcr);SFREE(ipivr);SFREE(adview);SFREE(auview);SFREE(adrad);
    SFREE(aurad);SFREE(irowrad);SFREE(jqrad);SFREE(icolrad);
    if((*mcs>0)&&(ntr>0)){SFREE(inocs);}
    if((*network>0)||(ntg>0)){SFREE(iponoel);SFREE(inoel);}
  }
  
  SFREE(fini);
  if((*nmethod==2)||(*nmethod==4)){
    SFREE(aux2);SFREE(fextini);SFREE(veini);
    SFREE(adb);SFREE(aub);SFREE(cv);SFREE(cvini);
  }
  
  SFREE(vini);SFREE(h0ref);SFREE(h0);SFREE(inomat);
  
  /* reset icascade */
  
  if(icascade==1){icascade=0;}
  
  mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
  mpcinfo[3]=maxlenmpc;
  
  *icolp=icol;*irowp=irow;*jqp=jq;*cop=co;*voldp=vold;
  
  *ipompcp=ipompc;*labmpcp=labmpc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
  *fmpcp=fmpc;*nodempcp=nodempc;*coefmpcp=coefmpc;
  
  *ipkonp=ipkon;*lakonp=lakon;*konp=kon;*ielorienp=ielorien;
  *ielmatp=ielmat;*enerp=ener;*xstatep=xstate;

  *setp=set;*istartsetp=istartset;*iendsetp=iendset;*ialsetp=ialset;
  *tiesetp=tieset;*tietolp=tietol;
  
  *nelemloadp=nelemload;*iamloadp=iamload;
  *sideloadp=sideload;
  
  (*tmin)*=(*tper);
  (*tmax)*=(*tper);
  
  SFREE(nactdofinv);
 
  if(*nmethod==1){
    *nmethod=8;
  }else if(*nmethod==4){
    *nmethod=9;
  }else if(*nmethod==2){
    *nmethod=10;
  }

  (*ttime)+=(*tper);
  
  return;
}
