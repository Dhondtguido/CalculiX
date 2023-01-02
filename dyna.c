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

void dyna(double **cop,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,ITG *ne,
	  ITG **nodebounp,ITG **ndirbounp,double **xbounp,ITG *nboun,
	  ITG **ipompcp,ITG **nodempcp,double **coefmpcp,char **labmpcp,
	  ITG *nmpc,ITG *nodeforc,ITG *ndirforc,double *xforc,
	  ITG *nforc,ITG *nelemload,char *sideload,double *xload,
	  ITG *nload,
	  ITG **nactdofp,ITG *neq,ITG *nzl,ITG *icol,ITG *irow,
	  ITG *nmethod,ITG **ikmpcp,ITG **ilmpcp,ITG **ikbounp,
	  ITG **ilbounp,double *elcon,ITG *nelcon,double *rhcon,
	  ITG *nrhcon,double *cocon,ITG *ncocon,
	  double *alcon,ITG *nalcon,double *alzero,
	  ITG **ielmatp,ITG **ielorienp,ITG *norien,double *orab,
	  ITG *ntmat_,double **t0p,
	  double **t1p,ITG *ithermal,double *prestr,ITG *iprestr,
	  double **voldp,ITG *iperturb,double **stip,ITG *nzs,
	  double *timepar,double *xmodal,
	  double **veoldp,char *amname,double *amta,
	  ITG *namta,ITG *nam,ITG *iamforc,ITG *iamload,
	  ITG **iamt1p,ITG *jout,
	  ITG *kode,char *filab,double **emep,double *xforcold,
	  double *xloadold,
	  double **t1oldp,ITG **iambounp,double **xbounoldp,ITG *iexpl,
	  double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
	  double **xstatep,ITG *npmat_,char *matname,ITG *mi,
	  ITG *ncmat_,ITG *nstate_,double **enerp,char *jobnamec,
	  double *ttime,char *set,ITG *nset,ITG *istartset,
	  ITG *iendset,ITG **ialsetp,ITG *nprint,char *prlab,
	  char *prset,ITG *nener,double *trab,
	  ITG **inotrp,ITG *ntrans,double **fmpcp,ITG *ipobody,ITG *ibody,
	  double *xbody,ITG *nbody,double *xbodyold,ITG *istep,
	  ITG *isolver,ITG *jq,char *output,ITG *mcs,ITG *nkon,
	  ITG *mpcend,ITG *ics,double *cs,ITG *ntie,char *tieset,
	  ITG *idrct,ITG *jmax,
	  double *ctrl,ITG *itpamp,double *tietol,ITG *nalset,
	  ITG *ikforc,ITG *ilforc,double *thicke,ITG *nslavs,ITG *nmat,
	  char *typeboun,ITG *ielprop,double *prop,char *orname,
	  double *t0g,double *t1g){

  char fneig[132]="",description[13]="            ",*lakon=NULL,*labmpc=NULL,
    *labmpcold=NULL,lakonl[9]="        \0",*tchar1=NULL,*tchar2=NULL,
    *tchar3=NULL,cflag[1]=" ",jobnamef[396]="",*labmpc2=NULL;

  ITG nev,i,j,k,idof,*inum=NULL,id,intscheme=0,
    iinc=0,jprint=0,l,iout,ielas=0,icmd=3,iprescribedboundary,init,ifreebody,
    mode=-1,noddiam=-1,*kon=NULL,*ipkon=NULL,*ielmat=NULL,*ielorien=NULL,
    *inotr=NULL,*nodeboun=NULL,*ndirboun=NULL,*iamboun=NULL,*ikboun=NULL,
    *ilboun=NULL,*nactdof=NULL,*ipompc=NULL,*nodempc=NULL,*ikmpc=NULL,
    *ilmpc=NULL,nsectors,nmpcold,mpcendold,*ipompcold=NULL,*nodempcold=NULL,
    *ikmpcold=NULL,*ilmpcold=NULL,kflag=2,nmd,nevd,*nm=NULL,*iamt1=NULL,
    *itg=NULL,ntg=0,symmetryflag=0,inputformat=0,dashpot,lrw,liw,iddebdf=0,
    *iwork=NULL,ngraph=1,nkg,neg,ncont=0,ne0,nkon0,*itietri=NULL,
    *koncont=NULL,konl[20],imat,nope,kodem,indexe,j1,jdof,icutb=0,
    *ipneigh=NULL,*neigh=NULL,inext,itp=0,*islavact=NULL,
    ismallsliding=0,isteadystate,mpcfree,im,cyclicsymmetry,
    memmpc_,imax,*icole=NULL,*irowe=NULL,*jqe=NULL,nzse[3],
    nalset_=*nalset,*ialset=*ialsetp,*istartset_=NULL,*iendset_=NULL,
    *itiefac=NULL,*islavsurf=NULL,*islavnode=NULL,mt=mi[1]+1,
    *imastnode=NULL,*nslavnode=NULL,*nmastnode=NULL,mortar=0,*imastop=NULL,
    *iponoels=NULL,*inoels=NULL,*imddof=NULL,nmddof,kchdep=0,
    *ikactcont=NULL,nactcont,nactcont_=100,*ikactmech=NULL,nactmech,
    iabsload=0,*ipe=NULL,*ime=NULL,iprev=1,inonlinmpc=0,ielem,
    *imdnode=NULL,nmdnode,*imdboun=NULL,nmdboun,*imdmpc=NULL,
    nmdmpc,intpointvar,kmin,kmax,i1,ifacecount,*izdof=NULL,
    nzdof,iload,iforc,*iponoel=NULL,*inoel=NULL,*imdelem=NULL,nmdelem,
    nasym=0,*nshcon=NULL,nherm,icfd=0,*inomat=NULL,
    ialeatoric=0,network=0,iperturbsav,mscalmethod=0,
    *islavelinv=NULL,*irowtloc=NULL,*jqtloc=NULL,nboun2,
    *ndirboun2=NULL,*nodeboun2=NULL,nmpc2,*ipompc2=NULL,*nodempc2=NULL,
    *ikboun2=NULL,*ilboun2=NULL,*ikmpc2=NULL,*ilmpc2=NULL,mortartrafoflag=0;

  long long i2;

  double *d=NULL,*z=NULL,*b=NULL,*zeta=NULL,*stiini=NULL,
    *cd=NULL,*cv=NULL,*xforcact=NULL,*xloadact=NULL,*cc=NULL,
    *t1act=NULL,*ampli=NULL,*aa=NULL,*bb=NULL,*aanew=NULL,*bj=NULL,
    *v=NULL,*aamech=NULL,*emn=NULL,*cdn=NULL,ptime,
    *stn=NULL,*stx=NULL,*een=NULL,*adb=NULL,*xstiff=NULL,*bjp=NULL,
    *aub=NULL,*temp_array1=NULL,*temp_array2=NULL,*aux=NULL,
    *f=NULL,*fn=NULL,*xbounact=NULL,*epn=NULL,*xstateini=NULL,
    *enern=NULL,*xstaten=NULL,*eei=NULL,*enerini=NULL,*qfn=NULL,
    *qfx=NULL,*xbodyact=NULL,*cgr=NULL,*au=NULL,*vbounact=NULL,
    *abounact=NULL,dtime,reltime,*t0=NULL,*t1=NULL,*t1old=NULL,
    physcon[1],zetaj,dj,ddj,h1,h2,h3,h4,h5,h6,sum,aai,bbi,tstart,tend,
    qa[4],cam[5],*accold=NULL,bet,gam,*ad=NULL,sigma=0.,alpham,betam,
    *bact=NULL,*bmin=NULL,*co=NULL,*xboun=NULL,*xbounold=NULL,*vold=NULL,
    *eme=NULL,*ener=NULL,*coefmpc=NULL,*fmpc=NULL,*coefmpcold,*veold=NULL,
    *xini=NULL,*rwork=NULL,*adc=NULL,*auc=NULL,*zc=NULL,*rpar=NULL,
    *cg=NULL,*straight=NULL,xl[27],voldl[mt*9],elas[21],fnl[27],t0l,t1l,
    elconloc[21],veoldl[mt*9],setnull,deltmx,fextmax,dd,dtheta,dthetaref,
    theta,*vini=NULL,dthetaold,*bcont=NULL,*vr=NULL,*vi=NULL,*bcontini=NULL,
    *stnr=NULL,*stni=NULL,*vmax=NULL,*stnmax=NULL,precision,resultmaxprev,
    resultmax,func,funcp,fexp,fexm,fcos,fsin,sump,*bp=NULL,h14,senergy=0.0,
    *bv=NULL,*cstr=NULL,*aube=NULL,*adbe=NULL,*sti=*stip,time0=0.0,
    time=0.0,*xforcdiff=NULL,*xloaddiff=NULL,*xbodydiff=NULL,*t1diff=NULL,
    *xboundiff=NULL,*bprev=NULL,*bdiff=NULL,*areaslav=NULL,venergy=0.0,
    *springarea=NULL,*fextold=NULL,*eenmax=NULL,*fnr=NULL,*fni=NULL,
    *xmastnor=NULL,*emeini=NULL,*xstate=NULL,*clearini=NULL,
    *shcon=NULL,*xmr=NULL,*xmi=NULL,*xnoels=NULL,*pslavsurf=NULL,
    *pmastsurf=NULL,*cdnr=NULL,*cdni=NULL,*tinc,*tper,*tmin,*tmax,
    *energyini=NULL,*energy=NULL,alea,*fext=NULL,*smscale=NULL,
    *autloc=NULL,*xboun2=NULL,*coefmpc2=NULL,*fnext=NULL;

  FILE *f1;

  /* dummy variables for nonlinmpc */

  ITG *iaux=NULL,maxlenmpc,icascade=0,newstep=0,iit=-1,idiscon;

#ifdef SGI
  ITG token;
#endif

  /*    if iabsload=0: aamech is modified by the present incremental 
	contribution of fext
	iabsload=1: the last incremental contribution is
	subtracted before the new one is added to fext;
	this latter incremental contribution is used
	to update aamech
	iabsload=2: aamech is determined by the absolute
	contribution of fext (no incremental procedure
	for the load; this is necessary if
	- nonlinear MPC's are applied or
	- user dloads are applied */

  co=*cop;kon=*konp;ipkon=*ipkonp;lakon=*lakonp;ielmat=*ielmatp;
  ielorien=*ielorienp;inotr=*inotrp;nodeboun=*nodebounp;
  ndirboun=*ndirbounp;iamboun=*iambounp;xboun=*xbounp;xstate=*xstatep;
  xbounold=*xbounoldp;ikboun=*ikbounp;ilboun=*ilbounp;nactdof=*nactdofp;
  vold=*voldp;eme=*emep;ener=*enerp;ipompc=*ipompcp;nodempc=*nodempcp;
  coefmpc=*coefmpcp;labmpc=*labmpcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
  fmpc=*fmpcp;veold=*veoldp;iamt1=*iamt1p;t0=*t0p;t1=*t1p;t1old=*t1oldp;
  
  for(k=0;k<3;k++){
    strcpy1(&jobnamef[k*132],&jobnamec[k*132],132);
  }

  tinc=&timepar[0];
  tper=&timepar[1];
  tmin=&timepar[2];
  tmax=&timepar[3];

  if(ithermal[0]<=1){
    kmin=1;kmax=3;
  }else if(ithermal[0]==2){
    kmin=0;kmax=mi[1];if(kmax>2)kmax=2;
  }else{
    kmin=0;kmax=3;
  }

  dtime=*tinc;

  alpham=xmodal[0];
  betam=xmodal[1];

  dd=ctrl[16];deltmx=ctrl[26];alea=ctrl[54];

  ne0=*ne;
  
  /* determining nzl */

  *nzl=0;
  for(i=neq[1];i>0;i--){
    if(icol[i-1]>0){
      *nzl=i;
      break;
    }
  }


  /* opening the eigenvalue file and checking for cyclic symmetry */

  strcpy2(fneig,jobnamec,132);
  strcat(fneig,".eig");

  if((f1=fopen(fneig,"rb"))==NULL){
    printf(" *ERROR in dyna: cannot open eigenvalue file for reading");
    printf(" *INFO  in dyna: if there are problems reading the .eig file this may be due to:\n");
    printf("        1) the nonexistence of the .eig file\n");
    printf("        2) other boundary conditions than in the input deck\n");
    printf("           which created the .eig file\n\n");
    exit(0);
  }

  if(fread(&cyclicsymmetry,sizeof(ITG),1,f1)!=1){
    printf(" *ERROR in dyna reading the cyclic symmetry flag in the eigenvalue file");
    printf(" *INFO  in dyna: if there are problems reading the .eig file this may be due to:\n");
    printf("        1) the nonexistence of the .eig file\n");
    printf("        2) other boundary conditions than in the input deck\n");
    printf("           which created the .eig file\n\n");
    exit(0);
  }

  if(fread(&nherm,sizeof(ITG),1,f1)!=1){
    printf(" *ERROR in dyna reading the Hermitian flag in the eigenvalue file");
    printf(" *INFO  in dyna: if there are problems reading the .eig file this may be due to:\n");
    printf("        1) the nonexistence of the .eig file\n");
    printf("        2) other boundary conditions than in the input deck\n");
    printf("           which created the .eig file\n\n");
    exit(0);
  }

  if(nherm!=1){
    printf(" *ERROR in dyna: the eigenvectors in the .eig-file result\n");
    printf("       from a non-Hermitian eigenvalue problem. The modal\n");
    printf("       dynamic procedure cannot handle that yet\n\n");
    FORTRAN(stop,());
  }

  iperturbsav=iperturb[0];
  if(fread(&iperturb[0],sizeof(ITG),1,f1)!=1){
    printf(" *ERROR in dyna reading the perturbation flag in the eigenvalue file");
    printf(" *INFO  in dyna: if there are problems reading the .eig file this may be due to:\n");
    printf("        1) the nonexistence of the .eig file\n");
    printf("        2) other boundary conditions than in the input deck\n");
    printf("           which created the .eig file\n\n");
    exit(0);
  }
  if(iperturbsav!=iperturb[0]){
    if(iperturb[0]==0){
      printf(" *ERROR in dyna: the .eig-file was created without perturbation info\n");
      printf("        the present *MODAL DYNAMIC step, however, contains the\n");
      printf("        perturbation parameter on the *STEP card.\n");
    }else if(iperturb[0]==1){
      printf(" *ERROR in dyna: the .eig-file was created with perturbation info\n");
      printf("        the present *MODAL DYNAMIC step, however, does not contain the\n");
      printf("        perturbation parameter on the *STEP card.\n");
    }
    FORTRAN(stop,());
  }

  if(iperturb[0]==1){
    NNEW(vold,double,mt**nk);
    if(fread(vold,sizeof(double),mt**nk,f1)!=mt**nk){
      printf(" *ERROR in dyna reading the reference displacements in the eigenvalue file...");
      printf(" *INFO  in dyna: if there are problems reading the .eig file this may be due to:\n");
      printf("        1) the nonexistence of the .eig file\n");
      printf("        2) other boundary conditions than in the input deck\n");
      printf("           which created the .eig file\n\n");
      exit(0);
    }
  }

  /* copying the displacements from the previous step as initial conditions in
     vini (only if 
     - not cyclic symmetric; initial conditions for cyclic symmetric
       calculations are assumed to be zero
     - AND not perturbed: if this is a perturbation step the initial
       conditions are zero */

  if(!cyclicsymmetry){
    NNEW(vini,double,mt**nk);
    if(iperturb[0]==0){
      memcpy(&vini[0],&vold[0],sizeof(double)*mt**nk);
    }
  }

  /* creating imddof containing the degrees of freedom
     retained by the user and imdnode containing the nodes */

  nmddof=0;nmdnode=0;nmdboun=0;nmdmpc=0;nmdelem=0;

  NNEW(imddof,ITG,*nk*3);
  NNEW(imdnode,ITG,*nk);
  NNEW(imdboun,ITG,*nboun);
  NNEW(imdmpc,ITG,*nmpc);
  FORTRAN(createmddof,(imddof,&nmddof,istartset,iendset,
		       ialset,nactdof,ithermal,mi,imdnode,&nmdnode,
		       ikmpc,ilmpc,ipompc,nodempc,nmpc,
		       imdmpc,&nmdmpc,imdboun,&nmdboun,ikboun,
		       nboun,nset,ntie,tieset,set,lakon,kon,ipkon,labmpc,
		       ilboun,filab,prlab,prset,nprint,ne,&cyclicsymmetry));

  /* if results are requested in too many nodes, it is faster to 
     calculate the results in all nodes */

  if((nmdnode>*nk/2)&&(!cyclicsymmetry)){
    nmdnode=0;nmddof=0;nmdboun=0;nmdmpc=0;
  }

  if(nmdnode!=0){
    if(!cyclicsymmetry){
      for(i=0;i<*nload;i++){
	iload=i+1;
	FORTRAN(addimdnodedload,(nelemload,sideload,ipkon,kon,lakon,
				 &iload,imdnode,&nmdnode,ikmpc,ilmpc,ipompc,nodempc,nmpc,
				 imddof,&nmddof,nactdof,mi,imdmpc,&nmdmpc,imdboun,&nmdboun,
				 ikboun,nboun,ilboun,ithermal));
      }
      for(i=0;i<*nforc;i++){
	iforc=i+1;
	FORTRAN(addimdnodecload,(nodeforc,&iforc,imdnode,&nmdnode,xforc,
				 ikmpc,ilmpc,ipompc,nodempc,nmpc,imddof,&nmddof,
				 nactdof,mi,imdmpc,&nmdmpc,imdboun,&nmdboun,
				 ikboun,nboun,ilboun,ithermal));
      }
    }
	  
    /* determining the elements belonging to a given node */
      
    NNEW(iponoel,ITG,*nk);
    NNEW(inoel,ITG,2**nkon);
    FORTRAN(elementpernode,(iponoel,inoel,lakon,ipkon,kon,ne));
    NNEW(imdelem,ITG,*ne);

    /* storing the elements in which integration point results
       are needed; storing the nodes which are needed to
       calculate these results */

    FORTRAN(createmdelem,(imdnode,&nmdnode,
			  ikmpc,ilmpc,ipompc,nodempc,nmpc,imddof,&nmddof,
			  nactdof,mi,imdmpc,&nmdmpc,imdboun,&nmdboun,
			  ikboun,nboun,ilboun,ithermal,imdelem,&nmdelem,
			  iponoel,inoel,prlab,prset,nprint,lakon,set,nset,
			  ialset,ipkon,kon,istartset,iendset,nforc,
			  ikforc,ilforc));

    RENEW(imdelem,ITG,nmdelem);
    SFREE(iponoel);SFREE(inoel);
  }

  /* if results are requested in too many nodes, it is faster to 
     calculate the results in all nodes */

  if((nmdnode>*nk/2)&&(!cyclicsymmetry)){
    nmdnode=0;nmddof=0;nmdboun=0;nmdmpc=0;
  }
  
  /* subtracting 1 to comply with the C-convention */

  for(i=0;i<nmddof;i++){imddof[i]-=1;}
  RENEW(imddof,ITG,nmddof);
  RENEW(imdnode,ITG,nmdnode);
  RENEW(imdboun,ITG,nmdboun);
  RENEW(imdmpc,ITG,nmdmpc);

  nsectors=1;

  /* reading the eigenvalues / eigenmodes */

  if(!cyclicsymmetry){

    nkg=*nk;
    neg=*ne;

    if(fread(&nev,sizeof(ITG),1,f1)!=1){
      printf(" *ERROR in dyna reading the number of eigenvalues in the eigenvalue file");
      printf(" *INFO  in dyna: if there are problems reading the .eig file this may be due to:\n");
      printf("        1) the nonexistence of the .eig file\n");
      printf("        2) other boundary conditions than in the input deck\n");
      printf("           which created the .eig file\n\n");
      exit(0);
    }
      
    NNEW(d,double,nev);
      
    if(fread(d,sizeof(double),nev,f1)!=nev){
      printf(" *ERROR in dyna reading the eigenvalues in the eigenvalue file");
      printf(" *INFO  in dyna: if there are problems reading the .eig file this may be due to:\n");
      printf("        1) the nonexistence of the .eig file\n");
      printf("        2) other boundary conditions than in the input deck\n");
      printf("           which created the .eig file\n\n");
      exit(0);
    }

    for(i=0;i<nev;i++){
      if(d[i]>0){d[i]=sqrt(d[i]);}else{d[i]=0.;}
    }
      
    NNEW(ad,double,neq[1]);
    NNEW(adb,double,neq[1]);
    NNEW(au,double,nzs[2]);
    NNEW(aub,double,nzs[1]);
      
    if(fread(ad,sizeof(double),neq[1],f1)!=neq[1]){
      printf(" *ERROR in dyna reading the diagonal of the stiffness matrix in the eigenvalue file");
      printf(" *INFO  in dyna: if there are problems reading the .eig file this may be due to:\n");
      printf("        1) the nonexistence of the .eig file\n");
      printf("        2) other boundary conditions than in the input deck\n");
      printf("           which created the .eig file\n\n");
      exit(0);
    }
      
    if(fread(au,sizeof(double),nzs[2],f1)!=nzs[2]){
      printf(" *ERROR in dyna reading the off-diagonals of the stiffness matrix in the eigenvalue file");
      printf(" *INFO  in dyna: if there are problems reading the .eig file this may be due to:\n");
      printf("        1) the nonexistence of the .eig file\n");
      printf("        2) other boundary conditions than in the input deck\n");
      printf("           which created the .eig file\n\n");
      exit(0);
    }
      
    if(fread(adb,sizeof(double),neq[1],f1)!=neq[1]){
      printf(" *ERROR in dyna reading the diagonal of the mass matrix in the eigenvalue file");
      printf(" *INFO  in dyna: if there are problems reading the .eig file this may be due to:\n");
      printf("        1) the nonexistence of the .eig file\n");
      printf("        2) other boundary conditions than in the input deck\n");
      printf("           which created the .eig file\n\n");
      exit(0);
    }
      
    if(fread(aub,sizeof(double),nzs[1],f1)!=nzs[1]){
      printf(" *ERROR in dyna reading the off-diagonals of the mass matrix in the  eigenvalue file");
      printf(" *INFO  in dyna: if there are problems reading the .eig file this may be due to:\n");
      printf("        1) the nonexistence of the .eig file\n");
      printf("        2) other boundary conditions than in the input deck\n");
      printf("           which created the .eig file\n\n");
      exit(0);
    }
      
    NNEW(z,double,neq[1]*nev);
      
    if(fread(z,sizeof(double),neq[1]*nev,f1)!=neq[1]*nev){
      printf(" *ERROR in dyna reading the eigenvectors in the eigenvalue file");
      printf(" *INFO  in dyna: if there are problems reading the .eig file this may be due to:\n");
      printf("        1) the nonexistence of the .eig file\n");
      printf("        2) other boundary conditions than in the input deck\n");
      printf("           which created the .eig file\n\n");
      exit(0);
    }
  }
  else{
    nev=0;
    do{
      if(fread(&nmd,sizeof(ITG),1,f1)!=1){
	break;
      }
      if(fread(&nevd,sizeof(ITG),1,f1)!=1){
	printf(" *ERROR in dyna reading the number of eigenvalues for nodal diameter %" ITGFORMAT " in the eigenvalue file",nmd);
	printf(" *INFO  in dyna: if there are problems reading the .eig file this may be due to:\n");
	printf("        1) the nonexistence of the .eig file\n");
	printf("        2) other boundary conditions than in the input deck\n");
	printf("           which created the .eig file\n\n");
	exit(0);
      }
      if(nev==0){
	NNEW(d,double,nevd);
	NNEW(nm,ITG,nevd);
      }else{
	RENEW(d,double,nev+nevd);
	RENEW(nm,ITG,nev+nevd);
      }
	  
      if(fread(&d[nev],sizeof(double),nevd,f1)!=nevd){
	printf(" *ERROR in dyna reading the eigenvalues for nodal diameter %" ITGFORMAT " in the eigenvalue file",nmd);
	printf(" *INFO  in dyna: if there are problems reading the .eig file this may be due to:\n");
	printf("        1) the nonexistence of the .eig file\n");
	printf("        2) other boundary conditions than in the input deck\n");
	printf("           which created the .eig file\n\n");
	exit(0);
      }

      for(i=nev;i<nev+nevd;i++){
	if(d[i]>0){d[i]=sqrt(d[i]);}else{d[i]=0.;}
      }

      for(i=nev;i<nev+nevd;i++){nm[i]=nmd;}
	  
      if(nev==0){
	NNEW(ad,double,neq[1]);
	NNEW(au,double,nzs[1]);

	if(fread(ad,sizeof(double),neq[1],f1)!=neq[1]){
	  printf(" *ERROR in dyna reading the diagonal of the stiffness matrix in the eigenvalue file");
	  printf(" *INFO  in dyna: if there are problems reading the .eig file this may be due to:\n");
	  printf("        1) the nonexistence of the .eig file\n");
	  printf("        2) other boundary conditions than in the input deck\n");
	  printf("           which created the .eig file\n\n");
	  exit(0);
	}
	      
	if(fread(au,sizeof(double),nzs[1],f1)!=nzs[1]){
	  printf(" *ERROR in dyna reading the off-diagonals of the stiffness matrix in the eigenvalue file");
	  printf(" *INFO  in dyna: if there are problems reading the .eig file this may be due to:\n");
	  printf("        1) the nonexistence of the .eig file\n");
	  printf("        2) other boundary conditions than in the input deck\n");
	  printf("           which created the .eig file\n\n");
	  exit(0);
	}

	SFREE(ad);SFREE(au);

	NNEW(adb,double,neq[1]);
	NNEW(aub,double,nzs[1]);

	if(fread(adb,sizeof(double),neq[1],f1)!=neq[1]){
	  printf(" *ERROR in dyna reading the diagonal of the mass matrix in the eigenvalue file");
	  printf(" *INFO  in dyna: if there are problems reading the .eig file this may be due to:\n");
	  printf("        1) the nonexistence of the .eig file\n");
	  printf("        2) other boundary conditions than in the input deck\n");
	  printf("           which created the .eig file\n\n");
	  exit(0);
	}
	      
	if(fread(aub,sizeof(double),nzs[1],f1)!=nzs[1]){
	  printf(" *ERROR in dyna reading the off-diagonals of the mass matrix in the eigenvalue file");
	  printf(" *INFO  in dyna: if there are problems reading the .eig file this may be due to:\n");
	  printf("        1) the nonexistence of the .eig file\n");
	  printf("        2) other boundary conditions than in the input deck\n");
	  printf("           which created the .eig file\n\n");
	  exit(0);
	}
      }
	  
      if(nev==0){
	NNEW(z,double,neq[1]*nevd);
      }else{
	RENEW(z,double,(long long)neq[1]*(nev+nevd));
      }
	  
      if(fread(&z[(long long)neq[1]*nev],sizeof(double),neq[1]*nevd,f1)!=neq[1]*nevd){
	printf(" *ERROR in dyna reading the eigenvectors for nodal diameter %" ITGFORMAT " in the eigenvalue file",nmd);
	printf(" *INFO  in dyna: if there are problems reading the .eig file this may be due to:\n");
	printf("        1) the nonexistence of the .eig file\n");
	printf("        2) other boundary conditions than in the input deck\n");
	printf("           which created the .eig file\n\n");
	exit(0);
      }
      nev+=nevd;
    }while(1);

    /* determining the maximum amount of segments */

    for(i=0;i<*mcs;i++){
      //	  if(cs[17*i]>nsectors) nsectors=cs[17*i];
      if(cs[17*i]>nsectors) nsectors=(ITG)(cs[17*i]+0.5);
    }

    /* determining the maximum number of sectors to be plotted */

    for(j=0;j<*mcs;j++){
      if(cs[17*j+4]>ngraph) ngraph=(ITG)cs[17*j+4];
    }
    nkg=*nk*ngraph;
    neg=*ne*ngraph;

    /* allocating field for the expanded structure */

    RENEW(co,double,3**nk*nsectors);

    /* next line is necessary for multiple cyclic symmetry
       conditions */

    for(i=3**nk;i<3**nk*nsectors;i++){co[i]=0.;}
    if(*ithermal!=0){
      RENEW(t0,double,*nk*nsectors);
      RENEW(t1old,double,*nk*nsectors);
      RENEW(t1,double,*nk*nsectors);
      if(*nam>0) RENEW(iamt1,ITG,*nk*nsectors);
    }
    RENEW(nactdof,ITG,mt**nk*nsectors);
    if(*ntrans>0) RENEW(inotr,ITG,2**nk*nsectors);
    RENEW(kon,ITG,*nkon*nsectors);
    RENEW(ipkon,ITG,*ne*nsectors);
    for(i=*ne;i<*ne*nsectors;i++){ipkon[i]=-1;}
    RENEW(lakon,char,8**ne*nsectors);
    RENEW(ielmat,ITG,mi[2]**ne*nsectors);
    if(*norien>0) RENEW(ielorien,ITG,mi[2]**ne*nsectors);
    //      RENEW(z,double,(long long)neq[1]*nev*nsectors/2);

    RENEW(nodeboun,ITG,*nboun*nsectors);
    RENEW(ndirboun,ITG,*nboun*nsectors);
    if(*nam>0) RENEW(iamboun,ITG,*nboun*nsectors);
    RENEW(xboun,double,*nboun*nsectors);
    RENEW(xbounold,double,*nboun*nsectors);
    RENEW(ikboun,ITG,*nboun*nsectors);
    RENEW(ilboun,ITG,*nboun*nsectors);

    NNEW(ipompcold,ITG,*nmpc);
    NNEW(nodempcold,ITG,3**mpcend);
    NNEW(coefmpcold,double,*mpcend);
    NNEW(labmpcold,char,20**nmpc);
    NNEW(ikmpcold,ITG,*nmpc);
    NNEW(ilmpcold,ITG,*nmpc);

    for(i=0;i<*nmpc;i++){ipompcold[i]=ipompc[i];}
    for(i=0;i<3**mpcend;i++){nodempcold[i]=nodempc[i];}
    for(i=0;i<*mpcend;i++){coefmpcold[i]=coefmpc[i];}
    for(i=0;i<20**nmpc;i++){labmpcold[i]=labmpc[i];}
    for(i=0;i<*nmpc;i++){ikmpcold[i]=ikmpc[i];}
    for(i=0;i<*nmpc;i++){ilmpcold[i]=ilmpc[i];}
    nmpcold=*nmpc;
    mpcendold=*mpcend;

    RENEW(ipompc,ITG,*nmpc*nsectors);
    RENEW(nodempc,ITG,3**mpcend*nsectors);
    RENEW(coefmpc,double,*mpcend*nsectors);
    RENEW(labmpc,char,20**nmpc*nsectors+1);
    RENEW(ikmpc,ITG,*nmpc*nsectors);
    RENEW(ilmpc,ITG,*nmpc*nsectors);
    RENEW(fmpc,double,*nmpc*nsectors);

    /* determining the space needed to expand the
       contact surfaces */

    NNEW(tchar1,char,81);
    NNEW(tchar2,char,81);
    NNEW(tchar3,char,81);
    for(i=0; i<*ntie; i++){
      if(tieset[i*(81*3)+80]=='C'){

	//a contact constraint was found, so increase nalset

	memcpy(tchar2,&tieset[i*(81*3)+81],81);
	tchar2[80]='\0';
	memcpy(tchar3,&tieset[i*(81*3)+81+81],81);
	tchar3[80]='\0';
	for(j=0; j<*nset; j++){
	  memcpy(tchar1,&set[j*81],81);
	  tchar1[80]='\0';
	  if(strcmp(tchar1,tchar2)==0){

	    //dependent nodal surface was found

	    (*nalset)+=(iendset[j]-istartset[j]+1)*(nsectors);
	  }
	  else if(strcmp(tchar1,tchar3)==0){

	    //independent element face surface was found

	    (*nalset)+=(iendset[j]-istartset[j]+1)*(nsectors);
	  }
	}
      }
    }
    SFREE(tchar1);
    SFREE(tchar2);
    SFREE(tchar3);

    RENEW(ialset,ITG,*nalset);

    /* save the information in istarset and isendset */

    NNEW(istartset_,ITG,*nset);
    NNEW(iendset_,ITG,*nset);
    for(j=0; j<*nset; j++){
      istartset_[j]=istartset[j];
      iendset_[j]=iendset[j];
    }

    /* reallocating the fields for the nodes in which the
       solution has to be calculated */

    RENEW(imddof,ITG,neq[1]/2*nsectors);
    RENEW(imdnode,ITG,*nk*nsectors);
    RENEW(imdboun,ITG,*nboun*nsectors);
    RENEW(imdmpc,ITG,*nmpc*nsectors);

    if(iperturb[0]==1){
      RENEW(vold,double,mt**nk*nsectors);
    }else{
      SFREE(vold);NNEW(vold,double,mt**nk*nsectors);
    }

    expand(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,nboun,
	   ipompc,nodempc,coefmpc,labmpc,nmpc,nodeforc,ndirforc,xforc,
	   nforc,nelemload,sideload,xload,nload,nactdof,neq,
	   nmethod,ikmpc,ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
	   alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
	   t0,ithermal,prestr,iprestr,vold,iperturb,sti,nzs,
	   adb,aub,filab,eme,plicon,nplicon,plkcon,nplkcon,
	   xstate,npmat_,matname,mi,ics,cs,mpcend,ncmat_,
	   nstate_,mcs,nkon,ener,jobnamec,output,set,nset,istartset,
	   iendset,ialset,nprint,prlab,prset,nener,trab,
	   inotr,ntrans,ttime,fmpc,&nev,&z,iamboun,xbounold,
	   &nsectors,nm,icol,irow,nzl,nam,ipompcold,nodempcold,coefmpcold,
	   labmpcold,&nmpcold,xloadold,iamload,t1old,t1,iamt1,xstiff,
	   &icole,&jqe,&irowe,isolver,nzse,&adbe,&aube,iexpl,
	   ibody,xbody,nbody,cocon,ncocon,tieset,ntie,imddof,&nmddof,
	   imdnode,&nmdnode,imdboun,&nmdboun,imdmpc,&nmdmpc,&izdof,&nzdof,
	   &nherm,xmr,xmi,typeboun,ielprop,prop,orname,itiefac,t0g,t1g);

    RENEW(imddof,ITG,nmddof);
    RENEW(imdnode,ITG,nmdnode);
    RENEW(imdboun,ITG,nmdboun);
    RENEW(imdmpc,ITG,nmdmpc);

    SFREE(veold);
    NNEW(veold,double,mt**nk);
    RENEW(eme,double,6*mi[0]**ne);
    RENEW(sti,double,6*mi[0]**ne);

    //      RENEW(xstiff,double,(long long)27*mi[0]**ne);
    if(*nener==1) RENEW(ener,double,2*mi[0]**ne);
  }

  fclose(f1);
	  
  /* checking for steadystate calculations */

  if(*tper<0){
    precision=-*tper;
    *tper=1.e10;
    isteadystate=1;
  }else{
    isteadystate=0;
  }

  /* checking for nonlinear MPC's */

  for(i=0;i<*nmpc;i++){
    if((strcmp1(&labmpc[20*i]," ")!=0)&&
       (strcmp1(&labmpc[20*i],"CYCLIC")!=0)&&
       (strcmp1(&labmpc[20*i],"SUBCYCLIC")!=0)){
      inonlinmpc=1;
      iabsload=2;
      break;
    }
  }

  /* in case this is a perturbation step (which cannot occur because
     of routine modaldynamics.f unless called by TRACE) the (d)loads are
     calculated in a non-incremental way in order to take the change
     of facial areas into account */
  
  if(iperturb[1]==1) iabsload=2;

  /* normalizing the time */

  FORTRAN(checktime,(itpamp,namta,tinc,ttime,amta,tmin,&inext,&itp,istep,tper));
  dtheta=(*tinc)/(*tper);
  dthetaref=dtheta;
  dthetaold=dtheta;

  *tmin=*tmin/(*tper);
  *tmax=*tmax/(*tper);
  theta=0.;

  /* check for rigid body modes 
     if there is a jump of 1.e4 in two subsequent eigenvalues
     all eigenvalues preceding the jump are considered to
     be rigid body modes and their frequency is set to zero */

  setnull=1.;
  for(i=nev-2;i>-1;i--){
    if(fabs(d[i])<0.0001*fabs(d[i+1])) setnull=0.;
    d[i]*=setnull;
  }

  /* check whether there are dashpot elements */

  dashpot=0;
  for(i=0;i<*ne;i++){
    if(ipkon[i]<0) continue;
    if(strcmp1(&lakon[i*8],"ED")==0){
      dashpot=1;break;}
  }

  if(dashpot){

    if(cyclicsymmetry){
      printf(" *ERROR in dyna: dashpots are not allowed in combination with cyclic symmetry\n");
      FORTRAN(stop,());
    }

    liw=51;
    NNEW(iwork,ITG,liw);
    lrw=130+42*nev;
    NNEW(rwork,double,lrw);
    NNEW(xini,double,2*nev);
    NNEW(adc,double,neq[1]);
    NNEW(auc,double,nzs[1]);
    FORTRAN(mafilldm,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xboun,nboun,
		      ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
		      nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
		      adc,auc,nactdof,icol,jq,irow,neq,nzl,nmethod,
		      ikmpc,ilmpc,ikboun,ilboun,
		      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		      ielorien,norien,orab,ntmat_,
		      t0,t0,ithermal,prestr,iprestr,vold,iperturb,sti,
		      nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
		      xstiff,npmat_,&dtime,matname,mi,ncmat_,
		      ttime,&time0,istep,&iinc,ibody,clearini,&mortar,springarea,
		      pslavsurf,pmastsurf,&reltime,&nasym));

    /*  zc = damping matrix * eigenmodes */

    NNEW(zc,double,neq[1]*nev);
    for(i=0;i<nev;i++){
      FORTRAN(op,(&neq[1],&z[i*neq[1]],&zc[i*neq[1]],adc,auc,
		  jq,irow));
    }

    /* cc is the reduced damping matrix (damping matrix mapped onto
       space spanned by eigenmodes) */

    NNEW(cc,double,nev*nev);
    for(i=0;i<nev;i++){
      for(j=0;j<=i;j++){
	for(k=0;k<neq[1];k++){
	  cc[i*nev+j]+=z[j*neq[1]+k]*zc[i*neq[1]+k];
	}
      }
    }

    /* symmetric part of cc matrix */

    for(i=0;i<nev;i++){
      for(j=i;j<nev;j++){
	cc[i*nev+j]=cc[j*nev+i];
      }
    }
    SFREE(zc);
  }

  /* contact conditions */

  /* storing the element and topology information before introducing 
     contact elements */

  ne0=*ne;nkon0=*nkon;

  NNEW(zeta,double,nev);
  NNEW(cstr,double,6);

  /* calculating the damping coefficients*/

  if(xmodal[10]<0){
    for(i=0;i<nev;i++){
      if(fabs(d[i])>(1.e-10)){
	zeta[i]=(alpham+betam*d[i]*d[i])/(2.*d[i]);
      }
      else {
	printf(" *WARNING in dyna: one of the frequencies is zero\n");
	printf("         no Rayleigh mass damping allowed\n");
	zeta[i]=0.;
      }

      /* if the nodal diameter exceeds half the number of sectors
	 the sign of the damping has to be reversed (omega is negative) */

      /*	if(cyclicsymmetry){
		if(nm[i]>nsectors/2) zeta[i]*=-1.;
		}*/
    }
  }
  //  else{
  else if((ITG)xmodal[10]>0){

    /*copy the damping coefficients for every eigenfrequency from xmodal[11....] */

    if(nev<(ITG)xmodal[10]){
      imax=nev;
      printf(" *WARNING in dyna: too many modal damping coefficients applied\n");
      printf("         damping coefficients corresponding to nonexisting eigenvalues are ignored\n");
    }
    else{
      imax=(ITG)xmodal[10];
    }
    for(i=0; i<imax; i++){
      zeta[i]=xmodal[11+i];     

      /* if the nodal diameter exceeds half the number of sectors
	 the sign of the damping has to be reversed (omega is negative) */

      /*   if(cyclicsymmetry){
	   if(nm[i]>nsectors/2) zeta[i]*=-1.;
	   }*/
    }
    
  }

  /* modal decomposition of the initial conditions */
  /* for cyclic symmetric structures the initial conditions
     are assumed to be zero */
  
  NNEW(cd,double,nev);
  NNEW(cv,double,nev);

  if(!cyclicsymmetry){
    NNEW(temp_array1,double,neq[1]);
    NNEW(temp_array2,double,neq[1]);
    for(i=0;i<neq[1];i++){temp_array1[i]=0;temp_array2[i]=0;}
    
    /* displacement initial conditions */
    
    for(i=0;i<*nk;i++){
      for(j=0;j<mt;j++){
	if(nactdof[mt*i+j]>0){
	  idof=nactdof[mt*i+j]-1;
	  //	  temp_array1[idof]=vold[mt*i+j];
	  temp_array1[idof]=vini[mt*i+j];
	}
      }
    }

    SFREE(vini);
    
    FORTRAN(op,(&neq[1],temp_array1,temp_array2,adb,aub,jq,irow));
    
    for(i=0;i<neq[1];i++){
      for(k=0;k<nev;k++){
	cd[k]+=z[k*neq[1]+i]*temp_array2[i];
      }
    }
    
    /* velocity initial conditions */
    
    for(i=0;i<neq[1];i++){temp_array1[i]=0;temp_array2[i]=0;}
    for(i=0;i<*nk;i++){
      for(j=0;j<mt;j++){
	if(nactdof[mt*i+j]>0){
	  idof=nactdof[mt*i+j]-1;
	  temp_array1[idof]=veold[mt*i+j];
	}
      }
    }
    
    FORTRAN(op,(&neq[1],temp_array1,temp_array2,adb,aub,jq,irow));
    
    for(i=0;i<neq[1];i++){
      for(k=0;k<nev;k++){
	cv[k]+=z[k*neq[1]+i]*temp_array2[i];
      }
    }
    
    SFREE(temp_array1);SFREE(temp_array2);
		      
  }
  NNEW(xforcact,double,*nforc);
  NNEW(xforcdiff,double,*nforc);
  NNEW(xloadact,double,2**nload);
  NNEW(xloaddiff,double,2**nload);
  NNEW(xbodyact,double,7**nbody);
  NNEW(xbodydiff,double,7**nbody);

  /* copying the rotation axis and/or acceleration vector */

  for(k=0;k<7**nbody;k++){xbodyact[k]=xbody[k];}
  for(k=0;k<7**nbody;k++){xbodydiff[k]=xbody[k];}
  NNEW(xbounact,double,*nboun);
  NNEW(xboundiff,double,*nboun);
  if(*ithermal==1) {NNEW(t1act,double,*nk);
    NNEW(t1diff,double,*nk);}
  
  /* assigning the body forces to the elements */ 

  /*  if(*nbody>0){
    ifreebody=*ne+1;
    NNEW(ipobody,ITG,2**ne);
    for(k=1;k<=*nbody;k++){
      FORTRAN(bodyforce,(cbody,ibody,ipobody,nbody,set,istartset,
			 iendset,ialset,&inewton,nset,&ifreebody,&k));
      RENEW(ipobody,ITG,2*(*ne+ifreebody));
    }
    RENEW(ipobody,ITG,2*(ifreebody-1));
    }*/
	       
  NNEW(fext,double,neq[1]); /* load rhs vector */
  NNEW(b,double,neq[1]); /* displacement solution vector */
  NNEW(bp,double,neq[1]); /* velocity solution vector */
  NNEW(bj,double,nev); /* response modal decomposition */
  NNEW(bjp,double,nev); /* derivative of the response modal decomposition */
  NNEW(ampli,double,*nam); /* instantaneous amplitude */

  /* constant coefficient of the linear amplitude function */

  NNEW(aa,double,nev); 
  NNEW(aanew,double,nev);
  NNEW(aamech,double,nev);

  /* linear coefficient of the linear amplitude function */

  NNEW(bb,double,nev);
  
  NNEW(v,double,mt**nk);
  NNEW(fn,double,mt**nk);
  NNEW(stn,double,6**nk);
  NNEW(inum,ITG,*nk);
  strcpy1(&cflag[0],&filab[4],1);
  FORTRAN(createinum,(ipkon,inum,kon,lakon,nk,ne,&cflag[0],nelemload,
		      nload,nodeboun,nboun,ndirboun,ithermal,co,vold,mi,ielmat,
                      ielprop,prop));

  if(*ithermal>1) {NNEW(qfn,double,3**nk);
    NNEW(qfx,double,3*mi[0]**ne);}

  if(strcmp1(&filab[261],"E   ")==0) NNEW(een,double,6**nk);
  if(strcmp1(&filab[522],"ENER")==0) NNEW(enern,double,*nk);
  if(strcmp1(&filab[609],"SDV ")==0) NNEW(xstaten,double,*nstate_**nk);
  if(strcmp1(&filab[2697],"ME  ")==0) NNEW(emn,double,6**nk);

  NNEW(eei,double,6*mi[0]**ne);
  if(*nener==1){
    NNEW(stiini,double,6*mi[0]**ne);
    NNEW(emeini,double,6*mi[0]**ne);
    NNEW(enerini,double,2*mi[0]**ne);}

  /*  check for nonzero SPC's */

  iprescribedboundary=0;
  for(i=0;i<*nboun;i++){
    if(fabs(xboun[i])>1.e-10){
      iprescribedboundary=1;
      if((*isolver==2)||(*isolver==3)){
	printf(" *ERROR in dyna: the iterative\n");
	printf("        solver cannot be used for modal dynamic\n");
	printf("        calculations with prescribed boundary\n");
	printf("        conditions\n");
	FORTRAN(stop,());
      }
      nmdnode=0;nmddof=0;nmdboun=0;nmdmpc=0;
      break;
    }
  }

  /*     calculating the instantaneous loads (forces, surface loading, 
	 centrifugal and gravity loading or temperature) at time 0 
	 setting iabsload to 2 if user subroutine dload is used */

  /*  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,
      xload,xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,
      xbodyact,t1old,t1,t1act,iamt1,nk,
      amta,namta,nam,ampli,&time0,&reltime,ttime,&dtime,ithermal,nmethod,
      xbounold,xboun,xbounact,iamboun,nboun,
      nodeboun,ndirboun,nodeforc,ndirforc,istep,&iinc,
      co,vold,itg,&ntg,amname,ikboun,ilboun,nelemload,sideload,mi));*/
      
  FORTRAN(temploaddiff,(xforcold,xforc,xforcact,iamforc,nforc,
			xloadold,xload,xloadact,iamload,nload,ibody,xbody,
			nbody,xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
			namta,nam,ampli,&time,&reltime,ttime,&dtime,ithermal,
			nmethod,xbounold,xboun,xbounact,iamboun,nboun,nodeboun,
			ndirboun,nodeforc,
			ndirforc,istep,&iinc,co,vold,itg,&ntg,amname,ikboun,ilboun,
			nelemload,sideload,mi,
			xforcdiff,xloaddiff,xbodydiff,t1diff,xboundiff,&iabsload,
			&iprescribedboundary,ntrans,trab,inotr,veold,nactdof,bcont,
			fn,ipobody,iponoel,inoel,ipkon,kon,lakon,ielprop,prop,ielmat,
			shcon,nshcon,rhcon,nrhcon,ntmat_,cocon,ncocon));

  if(iabsload==2) NNEW(fextold,double,neq[1]);

  /*  calculating the instantaneous loading vector at time 0 */

  NNEW(ikactmech,ITG,neq[1]);
  nactmech=0;
  rhsmain(co,nk,kon,ipkon,lakon,ne,
	  ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
	  nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,nbody,
	  cgr,fext,nactdof,&neq[1],nmethod,
	  ikmpc,ilmpc,elcon,nelcon,rhcon,nrhcon,alcon,
	  nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,t0,t1act,
	  ithermal,iprestr,vold,iperturb,iexpl,plicon,
	  nplicon,plkcon,nplkcon,npmat_,ttime,&time0,istep,&iinc,&dtime,
	  physcon,ibody,xbodyold,&reltime,veold,matname,mi,ikactmech,
	  &nactmech,ielprop,prop,sti,xstateini,xstate,nstate_,ntrans,
	  inotr,trab,fnext);
  
  /*  correction for nonzero SPC's */
  
  if(iprescribedboundary){

    if(cyclicsymmetry){
      printf(" *ERROR in dyna: prescribed boundaries are not allowed in combination with cyclic symmetry\n");
      FORTRAN(stop,());
    }

    if(*idrct!=1){
      printf(" *ERROR in dyna: variable increment length is not allwed in combination with prescribed boundaries\n");
      FORTRAN(stop,());
    }
      
    /* LU decomposition of the stiffness matrix */
      
    if(*isolver==0){
#ifdef SPOOLES
      spooles_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
		     &symmetryflag,&inputformat,&nzs[2]);
#else
      printf(" *ERROR in dyna: the SPOOLES library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==4){
#ifdef SGI
      token=1;
      sgi_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],token);
#else
      printf(" *ERROR in dyna: the SGI library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==5){
#ifdef TAUCS
      tau_factor(ad,&au,adb,aub,&sigma,icol,&irow,&neq[1],&nzs[1]);
#else
      printf(" *ERROR in dyna: the TAUCS library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==7){
#ifdef PARDISO
      pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
		     &symmetryflag,&inputformat,jq,&nzs[2]);
#else
      printf(" *ERROR in dyna: the PARDISO library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }
    else if(*isolver==8){
#ifdef PASTIX
      pastix_factor_main(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
		    &symmetryflag,&inputformat,jq,&nzs[2]);
#else
      printf(" *ERROR in arpack: the PASTIX library is not linked\n\n");
      FORTRAN(stop,());
#endif
    }

    NNEW(bact,double,neq[1]);
    NNEW(bmin,double,neq[1]);
    NNEW(bv,double,neq[1]);
    NNEW(bprev,double,neq[1]);
    NNEW(bdiff,double,neq[1]);

    init=1;
    dynboun(amta,namta,nam,ampli,&time0,ttime,&dtime,xbounold,xboun,
	    xbounact,iamboun,nboun,nodeboun,ndirboun,ad,au,adb,
	    aub,icol,irow,neq,nzs,&sigma,fext,isolver,
	    &alpham,&betam,nzl,&init,bact,bmin,jq,amname,bv,
	    bprev,bdiff,&nactmech,&iabsload,&iprev,&reltime);
    init=0;
  }

  iit=1;

  /* load at the start of a new step:
     mechanical loading without contact */

  if(!cyclicsymmetry){
    for(i=0;i<nev;i++){
      i2=(long long)i*neq[1];
      aamech[i]=0.;
      if(nactmech<neq[1]/2){
	for(j=0;j<nactmech;j++){
	  aamech[i]+=z[i2+ikactmech[j]]*fext[ikactmech[j]];
	}
      }else{
	for(j=0;j<neq[1];j++){
	  aamech[i]+=z[i2+j]*fext[j];
	}
      }
      aanew[i]=aamech[i];
    }
  }else{
    for(i=0;i<nev;i++){aamech[i]=0.;}
    for(j=0;j<nactmech;j++){
      FORTRAN(nident,(izdof,&ikactmech[j],&nzdof,&id));
      if(id!=0){
	if(izdof[id-1]==ikactmech[j]){
	  for(i=0;i<nev;i++){
	    aamech[i]+=z[(long long)i*nzdof+id-1]*fext[ikactmech[j]];
	  }
	}else{printf(" *ERROR in dyna\n");FORTRAN(stop,());}
      }else{printf(" *ERROR in dyna\n");FORTRAN(stop,());}
    }
    memcpy(&aanew[0],&aamech[0],sizeof(double)*nev);
  }

  /* check whether integration point values are requested; if not,
     the stress fields do not have to be allocated */

  intpointvar=0;
  if(*ithermal<=1){

    /* mechanical */

    if((strcmp1(&filab[174],"S")==0)||
       (strcmp1(&filab[261],"E")==0)||
       (strcmp1(&filab[348],"RF")==0)||
       (strcmp1(&filab[435],"PEEQ")==0)||
       (strcmp1(&filab[522],"ENER")==0)||
       (strcmp1(&filab[609],"SDV")==0)||
       (strcmp1(&filab[1044],"ZZS")==0)||
       (strcmp1(&filab[1044],"ERR")==0)||
       (strcmp1(&filab[1479],"PHS")==0)||
       (strcmp1(&filab[1653],"MAXS")==0)||
       (strcmp1(&filab[2175],"CONT")==0)||
       (strcmp1(&filab[2262],"CELS")==0)) intpointvar=1;
    for(i=0;i<*nprint;i++){
      if((strcmp1(&prlab[6*i],"S")==0)||
	 (strcmp1(&prlab[6*i],"E")==0)||
	 (strcmp1(&prlab[6*i],"PEEQ")==0)||
	 (strcmp1(&prlab[6*i],"ENER")==0)||
	 (strcmp1(&prlab[6*i],"SDV")==0)||
	 (strcmp1(&prlab[6*i],"CDIS")==0)||
	 (strcmp1(&prlab[6*i],"CSTR")==0)||
	 (strcmp1(&prlab[6*i],"CELS")==0)||
	 (strcmp1(&prlab[6*i],"RF")==0)) {intpointvar=1;break;}
    }
  }else{

    /* thermal */

    if((strcmp1(&filab[696],"HFL")==0)||
       (strcmp1(&filab[783],"RFL")==0)) intpointvar=1;
    for(i=0;i<*nprint;i++){
      if((strcmp1(&prlab[6*i],"HFL")==0)||
	 (strcmp1(&prlab[6*i],"RFL")==0)) {intpointvar=1;break;}
    }
  }

  if((intpointvar==1)) NNEW(stx,double,6*mi[0]**ne);

  /* major loop */

  resultmaxprev=0.;
  resultmax=0.;

  while(1.-theta>1.e-6){
    
    time0=time;
    
    //    printf("\nnew increment\n");
    
    if(*nener==1){
      memcpy(&enerini[0],&ener[0],sizeof(double)*2*mi[0]*ne0);
      if(*ithermal!=2){
	memcpy(&stiini[0],&sti[0],sizeof(double)*6*mi[0]*ne0);
	memcpy(&emeini[0],&eme[0],sizeof(double)*6*mi[0]*ne0);
      }
    }

    iinc++;
    jprint++;

    if(dashpot)RENEW(rpar,double,4+nev*(3+nev));
      
    /* check for max. # of increments */
    
    if(iinc>*jmax){
      printf(" *ERROR in dyna: max. # of increments reached\n\n");
      FORTRAN(stop,());
    }
    
    if(iinc>1){
      memcpy(&cd[0],&bj[0],sizeof(double)*nev);
      memcpy(&cv[0],&bjp[0],sizeof(double)*nev);
    }

      
    if((*idrct!=1)&&(iinc!=1)){
	
      /* increasing the increment size */
	
      dthetaold=dtheta;
      dtheta=dthetaref*dd;
	
      /* check increment length whether
	 - it does not exceed tmax
	 - the step length is not exceeded
	 - a time point is not exceeded  */
	
      dthetaref=dtheta;
      checkinclength(&time0,ttime,&theta,&dtheta,idrct,tper,tmax,
		     tmin,ctrl, amta,namta,itpamp,&inext,&dthetaref,&itp,
		     &jprint,jout);
    }
      
    reltime=theta+dtheta;
    time=reltime**tper;
    dtime=dtheta**tper;
      
    /* calculating the instantaneous loads (forces, surface loading, 
       centrifugal and gravity loading or temperature) */
      
    FORTRAN(temploaddiff,(xforcold,xforc,xforcact,iamforc,nforc,
			  xloadold,xload,xloadact,iamload,nload,ibody,xbody,
			  nbody,xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
			  namta,nam,ampli,&time,&reltime,ttime,&dtime,ithermal,
			  nmethod,xbounold,xboun,xbounact,iamboun,nboun,nodeboun,
			  ndirboun,nodeforc,
			  ndirforc,istep,&iinc,co,vold,itg,&ntg,amname,ikboun,ilboun,
			  nelemload,sideload,mi,
			  xforcdiff,xloaddiff,xbodydiff,t1diff,xboundiff,&iabsload,
			  &iprescribedboundary,ntrans,trab,inotr,veold,nactdof,bcont,
			  fn,ipobody,iponoel,inoel,ipkon,kon,lakon,ielprop,prop,ielmat,
			  shcon,nshcon,rhcon,nrhcon,ntmat_,cocon,ncocon));
	      
    /* calculating the instantaneous loading vector */
	      
    if(iabsload!=2){
      rhsmain(co,nk,kon,ipkon,lakon,ne,
	      ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcdiff,
	      nforc,nelemload,sideload,xloaddiff,nload,xbodydiff,
	      ipobody,nbody,cgr,fext,nactdof,&neq[1],nmethod,
	      ikmpc,ilmpc,elcon,nelcon,rhcon,nrhcon,
	      alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
	      t0,t1diff,ithermal,iprestr,vold,iperturb,iexpl,plicon,
	      nplicon,plkcon,nplkcon,
	      npmat_,ttime,&time,istep,&iinc,&dtime,physcon,ibody,
	      xbodyold,&reltime,veold,matname,mi,ikactmech,&nactmech,
	      ielprop,prop,sti,xstateini,xstate,nstate_,ntrans,inotr,
	      trab,fnext);
    }else{
      rhsmain(co,nk,kon,ipkon,lakon,ne,
	      ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
	      nforc,nelemload,sideload,xloadact,nload,xbodyact,
	      ipobody,nbody,cgr,fext,nactdof,&neq[1],nmethod,
	      ikmpc,ilmpc,elcon,nelcon,rhcon,nrhcon,
	      alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
	      t0,t1act,ithermal,iprestr,vold,iperturb,iexpl,plicon,
	      nplicon,plkcon,nplkcon,
	      npmat_,ttime,&time,istep,&iinc,&dtime,physcon,ibody,
	      xbodyold,&reltime,veold,matname,mi,ikactmech,&nactmech,
	      ielprop,prop,sti,xstateini,xstate,nstate_,ntrans,inotr,
	      trab,fnext);
    }
	      
    /* correction for nonzero SPC's */
      
    if(iprescribedboundary){
      dynboun(amta,namta,nam,ampli,&time,ttime,&dtime,
	      xbounold,xboun,
	      xbounact,iamboun,nboun,nodeboun,ndirboun,ad,au,adb,
	      aub,icol,irow,neq,nzs,&sigma,fext,isolver,
	      &alpham,&betam,nzl,&init,bact,bmin,jq,amname,bv,
	      bprev,bdiff,&nactmech,&iabsload,&iprev,&reltime);
    }

    if(*idrct==0){
      fextmax=0.;
      if(iabsload!=2){
	if(nactmech<neq[1]/2){
	  for(i=0;i<nactmech;i++){
	    if(fabs(fext[ikactmech[i]])>fextmax) fextmax=fabs(fext[ikactmech[i]]);
	  }
	}else{
	  for(i=0;i<neq[1];i++){
	    if(fabs(fext[i])>fextmax) fextmax=fabs(fext[i]);
	  }
	}
      }else{

	/* fextmax is to be calculated from the difference of fext and fextold */

	if(nactmech<neq[1]/2){
	  for(i=0;i<nactmech;i++){
	    if(fabs(fext[ikactmech[i]]-fextold[ikactmech[i]])>fextmax) 
	      fextmax=fabs(fext[ikactmech[i]]-fextold[ikactmech[i]]);
	  }
	}else{
	  for(i=0;i<neq[1];i++){
	    if(fabs(fext[i])>fextmax) fextmax=fabs(fext[i]-fextold[i]);
	  }
	}
	      
	/* copy fext into fextold */

	if(nactmech<neq[1]/2){
	  for(i=0;i<nactmech;i++){
	    fextold[ikactmech[i]]=fext[ikactmech[i]];
	  }
	}else{
	  memcpy(&fextold[0],&fext[0],sizeof(double)*neq[1]);
	}
      }

      /* check for size of mechanical force */

      if((fextmax>deltmx)&&(((itp==1)&&(dtheta>*tmin))||(itp==0))){

	/* force increase too big: increment size is decreased */

	if(iabsload==0) iabsload=1;
	dtheta=dtheta*deltmx/fextmax;
	dthetaref=dtheta;
	if(itp==1){
	  inext--;
	  itp=0;
	}
	      
	/* check whether the new increment size is not too small */
	      
	if(dtheta<*tmin){
	  dtheta=*tmin;
	}

	reltime=theta+dtheta;
	time=reltime**tper;
	dtime=dtheta**tper;
	      
	/* calculating the instantaneous loads (forces, surface loading, 
	   centrifugal and gravity loading or temperature) */
	      
	FORTRAN(temploaddiff,(xforcold,xforc,xforcact,iamforc,nforc,
			      xloadold,xload,xloadact,iamload,nload,ibody,xbody,
			      nbody,xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
			      namta,nam,ampli,&time,&reltime,ttime,&dtime,ithermal,
			      nmethod,xbounold,xboun,xbounact,iamboun,nboun,nodeboun,
			      ndirboun,nodeforc,
			      ndirforc,istep,&iinc,co,vold,itg,&ntg,amname,ikboun,ilboun,
			      nelemload,sideload,mi,
			      xforcdiff,xloaddiff,xbodydiff,t1diff,xboundiff,&iabsload,
			      &iprescribedboundary,ntrans,trab,inotr,veold,nactdof,bcont,
			      fn,ipobody,iponoel,inoel,ipkon,kon,lakon,ielprop,prop,ielmat,
			      shcon,nshcon,rhcon,nrhcon,ntmat_,cocon,ncocon));
	      
	/* calculating the instantaneous loading vector */
	      
	if(iabsload!=2){
	  rhsmain(co,nk,kon,ipkon,lakon,ne,
		  ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcdiff,
		  nforc,nelemload,sideload,xloaddiff,nload,xbodydiff,
		  ipobody,nbody,cgr,fext,nactdof,&neq[1],nmethod,
		  ikmpc,ilmpc,elcon,nelcon,rhcon,nrhcon,
		  alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
		  t0,t1diff,ithermal,iprestr,vold,iperturb,iexpl,plicon,
		  nplicon,plkcon,nplkcon,
		  npmat_,ttime,&time,istep,&iinc,&dtime,physcon,ibody,
		  xbodyold,&reltime,veold,matname,mi,ikactmech,&nactmech,
		  ielprop,prop,sti,xstateini,xstate,nstate_,ntrans,inotr,
		  trab,fnext);
	}else{
	  rhsmain(co,nk,kon,ipkon,lakon,ne,
		  ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
		  nforc,nelemload,sideload,xloadact,nload,xbodyact,
		  ipobody,nbody,cgr,fext,nactdof,&neq[1],nmethod,
		  ikmpc,ilmpc,elcon,nelcon,rhcon,nrhcon,
		  alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
		  t0,t1act,ithermal,iprestr,vold,iperturb,iexpl,plicon,
		  nplicon,plkcon,nplkcon,
		  npmat_,ttime,&time,istep,&iinc,&dtime,physcon,ibody,
		  xbodyold,&reltime,veold,matname,mi,ikactmech,&nactmech,
		  ielprop,prop,sti,xstateini,xstate,nstate_,ntrans,inotr,
		  trab,fnext);
	}
	      
	/* correction for nonzero SPC's */
	      
	if(iprescribedboundary){
	  dynboun(amta,namta,nam,ampli,&time,ttime,&dtime,
		  xbounold,xboun,
		  xbounact,iamboun,nboun,nodeboun,ndirboun,ad,au,adb,
		  aub,icol,irow,neq,nzs,&sigma,fext,isolver,
		  &alpham,&betam,nzl,&init,bact,bmin,jq,amname,bv,
		  bprev,bdiff,&nactmech,&iabsload,&iprev,&reltime);
	}
	if(iabsload==1) iabsload=0;
      }

    }

    /* step length is OK for mechanical load
       calculating equation for linearized loading */

    /* load: actual mechanical load +
       contact from last increment */

    if(!cyclicsymmetry){
      for(i=0;i<nev;i++){
	i2=(long long)i*neq[1];
	aa[i]=aanew[i];
	if(iabsload==2){aamech[i]=0.;}
	if(nactmech<neq[1]/2){
	  for(j=0;j<nactmech;j++){
	    aamech[i]+=z[i2+ikactmech[j]]*fext[ikactmech[j]];
	  }
	}else{
	  for(j=0;j<neq[1];j++){
	    aamech[i]+=z[i2+j]*fext[j];
	  }
	}
	      
	aanew[i]=aamech[i];
	      
	bb[i]=(aanew[i]-aa[i])/dtime;
	aa[i]=aanew[i]-bb[i]*time;
      }
    }else{
      for(i=0;i<nev;i++){
	memcpy(&aa[0],&aanew[0],sizeof(double)*nev);
	if(iabsload==2){aamech[i]=0.;}
      }
      for(j=0;j<nactmech;j++){
	FORTRAN(nident,(izdof,&ikactmech[j],&nzdof,&id));
	if(id!=0){
	  if(izdof[id-1]==ikactmech[j]){
	    for(i=0;i<nev;i++){
	      aamech[i]+=z[(long long)i*nzdof+id-1]*fext[ikactmech[j]];
	    }
	  }else{printf(" *ERROR in dyna\n");FORTRAN(stop,());}
	}else{printf(" *ERROR in dyna\n");FORTRAN(stop,());}
      }
      memcpy(&aanew[0],&aamech[0],sizeof(double)*nev);
      for(i=0;i<nev;i++){
	bb[i]=(aanew[i]-aa[i])/(dtime);
	aa[i]=aanew[i]-bb[i]*time;
      }	
    }
      
    /* calculating the response due to unchanged contact force during
       the increment */

    if(dashpot){
      FORTRAN(subspace,(d,aa,bb,cc,&alpham,&betam,&nev,xini,
			cd,cv,&time,rwork,&lrw,&iinc,jout,rpar,bj,
			iwork,&liw,&iddebdf,bjp));
      if(iddebdf==2){
	liw=56+2*nev;
	RENEW(iwork,ITG,liw);
	for(i=0;i<liw;i++){iwork[i]=0;}
	lrw=250+20*nev+4*nev*nev;
	RENEW(rwork,double,lrw);
	for(i=0;i<lrw;i++){rwork[i]=0.;}
	iddebdf=1;
	FORTRAN(subspace,(d,aa,bb,cc,&alpham,&betam,&nev,xini,
			  cd,cv,&time,rwork,&lrw,&iinc,jout,rpar,bj,
			  iwork,&liw,&iddebdf,bjp));
      }
    }
    else{
      for(l=0;l<nev;l++){
	zetaj=zeta[l];
	dj=d[l];
	      
	/* zero eigenfrequency: rigid body mode */
	      
	if(fabs(d[l])<=1.e-10){
	  aai=aa[l];
	  bbi=bb[l];
	  tstart=time0;
	  tend=time;
	  sum=tend*(aai*time+
		    tend*((bbi*time-aai)/2.-bbi*tend/3.))-
	    tstart*(aai*time+
		    tstart*((bbi*time-aai)/2.-bbi*tstart/3.));
	  sump=tend*(aai+bbi*tend/2.)-tstart*(aai+bbi*tstart/2.);
	  bj[l]=sum+cd[l]+dtime*cv[l];
	  bjp[l]=sump+cv[l];
	}
	      
	/*   subcritical damping */
	      
	else if(zetaj<1.-1.e-6){
	  ddj=dj*sqrt(1.-zetaj*zetaj);
	  h1=zetaj*dj;
	  h2=h1*h1+ddj*ddj;
	  h3=h1*h1-ddj*ddj;
	  h4=2.*h1*ddj/h2;
	  h14=h1/ddj;
	  tstart=0.;		      
	  FORTRAN(fsub,(&time,&dtime,&aa[l],&bb[l],&ddj,
			&h1,&h2,&h3,&h4,&func,&funcp));
	  sum=func;sump=funcp;
	  FORTRAN(fsub,(&time,&tstart,&aa[l],&bb[l],&ddj,
			&h1,&h2,&h3,&h4,&func,&funcp));
	  sum-=func;sump-=funcp;
	  fexp=exp(-h1*dtime);
	  fsin=sin(ddj*dtime);
	  fcos=cos(ddj*dtime);

	  bj[l]=sum/ddj+fexp*(fcos+zetaj/sqrt(1.-zetaj*zetaj)*fsin)*cd[l]+
	    fexp*fsin*cv[l]/ddj;
	  bjp[l]=sump/ddj+fexp*((-h1+ddj*h14)*fcos+(-ddj-h1*h14)*fsin)*cd[l]
	    +fexp*(-h1*fsin+ddj*fcos)*cv[l]/ddj;

	}
	      
	/*      supercritical damping */
	      
	else if(zetaj>1.+1.e-6){
	  ddj=dj*sqrt(zetaj*zetaj-1.);
	  h1=ddj-zetaj*dj;
	  h2=ddj+zetaj*dj;
	  h3=1./h1;
	  h4=1./h2;
	  h5=h3*h3;
	  h6=h4*h4;
	  tstart=0.;		      
	  FORTRAN(fsuper,(&time,&dtime,&aa[l],&bb[l],
			  &h1,&h2,&h3,&h4,&h5,&h6,&func,&funcp));
	  sum=func;sump=funcp;
	  FORTRAN(fsuper,(&time,&tstart,&aa[l],&bb[l],
			  &h1,&h2,&h3,&h4,&h5,&h6,&func,&funcp));
	  sum-=func;sump-=funcp;
		  
	  fexm=exp(h1*dtime);
	  fexp=exp(-h2*dtime);
	  h14=zetaj*dj/ddj;
	  bj[l]=sum/(2.*ddj)+(fexm+fexp)*cd[l]/2.+zetaj*(fexm-fexp)/(2.*
								     sqrt(zetaj*zetaj-1.))*cd[l]+(fexm-fexp)*cv[l]/(2.*ddj);
	  bjp[l]=sump/(2.*ddj)+(h1*fexm-h2*fexp)*cd[l]/2.
	    +(h14*cd[l]+cv[l]/ddj)*(h1*fexm+h2*fexp)/2.;
	}
	      
	/* critical damping */
	      
	else{
	  h1=zetaj*dj;
	  h2=1./h1;
	  h3=h2*h2;
	  h4=h2*h3;
	  tstart=0.;
	  FORTRAN(fcrit,(&time,&dtime,&aa[l],&bb[l],&zetaj,&dj,
			 &ddj,&h1,&h2,&h3,&h4,&func,&funcp));
	  sum=func;sump=funcp;
	  FORTRAN(fcrit,(&time,&tstart,&aa[l],&bb[l],&zetaj,&dj,
			 &ddj,&h1,&h2,&h3,&h4,&func,&funcp));
	  sum-=func;sump-=funcp;
	  fexp=exp(-h1*dtime);
	  bj[l]=sum+fexp*((1.+h1*dtime)*cd[l]+dtime*cv[l]);
	  bjp[l]=sump+fexp*(-h1*h1*dtime*cd[l]+
			    (1.-h1*dtime)*cv[l]);
	}
      }
    }
      
    /* composing the response */
      
    if(iprescribedboundary){
      if(nmdnode==0){
	memcpy(&b[0],&bmin[0],sizeof(double)*neq[1]);
	memcpy(&bp[0],&bv[0],sizeof(double)*neq[1]);
      }else{
	for(i=0;i<nmddof;i++){
	  b[imddof[i]]=bmin[imddof[i]];
	  bp[imddof[i]]=bv[imddof[i]];
	}
      }
    }
    else{
      if(nmdnode==0){
	DMEMSET(b,0,neq[1],0.);
	DMEMSET(bp,0,neq[1],0.);
      }else{
	for(i=0;i<nmddof;i++){
	  b[imddof[i]]=0.;
	  bp[imddof[i]]=0.;
	}
      }
    }
      
    if(!cyclicsymmetry){
      if(nmdnode==0){
	for(i=0;i<neq[1];i++){
	  for(j=0;j<nev;j++){
	    b[i]+=bj[j]*z[(long long)j*neq[1]+i];
	    bp[i]+=bjp[j]*z[(long long)j*neq[1]+i];
	  }
	}
      }else{
	for(i=0;i<nmddof;i++){
	  for(j=0;j<nev;j++){
	    b[imddof[i]]+=bj[j]*z[(long long)j*neq[1]+imddof[i]];
	    bp[imddof[i]]+=bjp[j]*z[(long long)j*neq[1]+imddof[i]];
	  }
	}
      }
    }else{
      for(i=0;i<nmddof;i++){
	FORTRAN(nident,(izdof,&imddof[i],&nzdof,&id));
	if(id!=0){
	  if(izdof[id-1]==imddof[i]){
	    for(j=0;j<nev;j++){
	      b[imddof[i]]+=bj[j]*z[(long long)j*nzdof+id-1];
	      bp[imddof[i]]+=bjp[j]*z[(long long)j*nzdof+id-1];
	    }
	  }else{printf(" *ERROR in dyna\n");FORTRAN(stop,());}
	}else{printf(" *ERROR in dyna\n");FORTRAN(stop,());}
      }
    }
      
    /* update nonlinear MPC-coefficients (e.g. for rigid
       body MPC's */

    if(inonlinmpc==1){
      FORTRAN(nonlinmpc,(co,vold,ipompc,nodempc,coefmpc,labmpc,
			 nmpc,ikboun,ilboun,nboun,xbounact,aux,iaux,
			 &maxlenmpc,ikmpc,ilmpc,&icascade,
			 kon,ipkon,lakon,ne,&reltime,&newstep,xboun,fmpc,
			 &iit,&idiscon,&ncont,trab,ntrans,ithermal,mi,
			 &kchdep));
    }
      
    /* calculating displacements/temperatures */
      
    FORTRAN(dynresults,(nk,v,ithermal,nactdof,vold,nodeboun,
			ndirboun,xbounact,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,
			b,bp,veold,&dtime,mi,imdnode,&nmdnode,imdboun,&nmdboun,
			imdmpc,&nmdmpc,nmethod,&time));
	  
    theta+=dtheta;
    //      (*ttime)+=dtime;
      
    /* check whether a time point was reached */
      
    if((*itpamp>0)&&(*idrct==0)){
      if(itp==1){
	jprint=*jout;
      }else{
	jprint=*jout+1;
      }
    }
      
    /* check whether output is needed */
      
    if((*jout==jprint)||(1.-theta<=1.e-6)){
      iout=2;
      jprint=0;
    }else if((*nener==1)){
      iout=-2;
    }else{
      iout=0;
    }
      
    if((iout==2)||(iout==-2)){

      /* deactivating the elements for which the stresses are not
	 needed */

      if(nmdnode>0){
	if((intpointvar==1)){
	  for(k=0;k<ne0;k++){
	    if(ipkon[k]<-1){
	      printf(" *ERROR in dyna: contact remeshing of quadratic elements is not allowed\n\n");
	      FORTRAN(stop,());
	    }else if(ipkon[k]!=-1){
	      ipkon[k]=-ipkon[k]-2;
	    }
	  }
	  for(k=0;k<nmdelem;k++){
	    ielem=imdelem[k]-1;
	    ipkon[ielem]=-2-ipkon[ielem];
	  }
	}
      }

      results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
	      stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
	      ielmat,ielorien,norien,orab,ntmat_,t0,t1,
	      ithermal,prestr,iprestr,filab,eme,emn,een,
	      iperturb,f,fn,nactdof,&iout,qa,
	      vold,b,nodeboun,ndirboun,xbounact,nboun,
	      ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],
	      veold,accold,&bet,&gam,&dtime,&time,ttime,
	      plicon,nplicon,plkcon,nplkcon,
	      xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	      &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,
	      enern,emeini,xstaten,eei,enerini,cocon,ncocon,
	      set,nset,istartset,iendset,ialset,nprint,prlab,prset,
	      qfx,qfn,trab,inotr,ntrans,fmpc,nelemload,nload,ikmpc,
	      ilmpc,istep,&iinc,springarea,&reltime,&ne0,
	      thicke,shcon,nshcon,
	      sideload,xload,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
	      &mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
	      islavsurf,ielprop,prop,energyini,energy,&iit,iponoel,
	      inoel,nener,orname,&network,ipobody,xbodyact,ibody,typeboun,
	      itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
	      islavelinv,autloc,irowtloc,jqtloc,&nboun2,
	      ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
	      labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
	      &intscheme,physcon);

      /* restoring */

      if(nmdnode>0){
	if((intpointvar==1)){
	  for(k=0;k<ne0;k++){
	    if(ipkon[k]<-1){ipkon[k]=-2-ipkon[k];}
	  }
	}
      }
	
      if((*ithermal!=2)&&(intpointvar==1)){
	for(k=0;k<6*mi[0]*ne0;++k){
	  sti[k]=stx[k];
	}
      }
    }
    if(iout==2){
      (*kode)++;
      if(strcmp1(&filab[1044],"ZZS")==0){
	NNEW(neigh,ITG,40**ne);
	NNEW(ipneigh,ITG,*nk);
      }

      ptime=*ttime+time;
      frd(co,&nkg,kon,ipkon,lakon,&neg,v,stn,inum,nmethod,
	  kode,filab,een,t1,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	  nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	  mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	  thicke,jobnamec,output,qfx,cdn,&mortar,cdnr,cdni,nmat,ielprop,
	  prop,sti);
	
      if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);}
    }
      
    if(isteadystate==1){
	  
      /* calculate maximum displacement/temperature */
	  
      resultmax=0.;
      if(*ithermal<2){
	for(i=1;i<mt**nk;i=i+mt){
	  if(fabs(v[i])>resultmax) resultmax=fabs(v[i]);}
	for(i=2;i<mt**nk;i=i+mt){
	  if(fabs(v[i])>resultmax) resultmax=fabs(v[i]);}
	for(i=3;i<mt**nk;i=i+mt){
	  if(fabs(v[i])>resultmax) resultmax=fabs(v[i]);}
      }else if(*ithermal==2){
	for(i=0;i<mt**nk;i=i+mt){
	  if(fabs(v[i])>resultmax) resultmax=fabs(v[i]);}
      }else{
	printf(" *ERROR in dyna: coupled temperature-displacement calculations are not allowed\n");
      }
      if(fabs((resultmax-resultmaxprev)/resultmax)<precision){
	break;
      }else{resultmaxprev=resultmax;}
    }
     
  }
      
  if((intpointvar==1)) SFREE(stx);

  /* calculating the displacements and velocities in all nodes as 
     initial condition for the next step; only needed if
     - nonzero initial conditions are allowed (-> no cyclic symmetry)
     - the output was restricted (-> nmdnode nonzero) */

  if((nmdnode!=0)&&(!cyclicsymmetry)){

    /* determining the solution in the independent nodes */

    DMEMSET(b,0,neq[1],0.);
    DMEMSET(bp,0,neq[1],0.);

    for(i=0;i<neq[1];i++){
      for(j=0;j<nev;j++){
	b[i]+=bj[j]*z[(long long)j*neq[1]+i];
	bp[i]+=bjp[j]*z[(long long)j*neq[1]+i];
      }
    }
      
    /* update nonlinear MPC-coefficients (e.g. for rigid
       body MPC's */

    if(inonlinmpc==1){
      FORTRAN(nonlinmpc,(co,vold,ipompc,nodempc,coefmpc,labmpc,
			 nmpc,ikboun,ilboun,nboun,xbounact,aux,iaux,
			 &maxlenmpc,ikmpc,ilmpc,&icascade,
			 kon,ipkon,lakon,ne,&reltime,&newstep,xboun,fmpc,
			 &iit,&idiscon,&ncont,trab,ntrans,ithermal,mi,
			 &kchdep));
    }
      
    /* calculating displacements/temperatures */
      
    nmdnode=0;
    FORTRAN(dynresults,(nk,v,ithermal,nactdof,vold,nodeboun,
			ndirboun,xbounact,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,
			b,bp,veold,&dtime,mi,imdnode,&nmdnode,imdboun,&nmdboun,
			imdmpc,&nmdmpc,nmethod,&time));
  }
  
  SFREE(eei);
  SFREE(vbounact);
  SFREE(abounact);

  if(*nener==1){SFREE(stiini);SFREE(emeini);SFREE(enerini);}

  if(strcmp1(&filab[261],"E   ")==0) SFREE(een);
  if(strcmp1(&filab[522],"ENER")==0) SFREE(enern);
  if(strcmp1(&filab[609],"SDV ")==0) SFREE(xstaten);
  if(strcmp1(&filab[2697],"ME  ")==0) SFREE(emn);
  if(*ithermal>1) {SFREE(qfn);SFREE(qfx);}

  /* updating the loading at the end of the step; 
     important in case the amplitude at the end of the step
     is not equal to one */

  /* this should not be done if the loading is user
     defined: user-defined loading is kept from one step
     to the next */

  for(k=0;k<*nboun;++k){
    if((xboun[k]<1.2357111316)||(xboun[k]>1.2357111318)){
      xboun[k]=xbounact[k];
    }
  }
  for(k=0;k<*nforc;++k){
    if((xforc[k]<1.2357111316)||(xforc[k]>1.2357111318)){
      xforc[k]=xforcact[k];
    }
  }
  for(k=0;k<2**nload;++k){xload[k]=xloadact[k];}
  for(k=0;k<7**nbody;k=k+7){xbody[k]=xbodyact[k];}
  if(*ithermal==1){
    for(k=0;k<*nk;++k){t1[k]=t1act[k];}
  }

  if(!cyclicsymmetry){

    /* copying the actual displacements into field vold
       serves as initial conditions in the next step */

    memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);
  }
  
  SFREE(v);SFREE(fn);SFREE(stn);SFREE(inum);SFREE(adb);SFREE(d);SFREE(fext);
  SFREE(aub);SFREE(z);SFREE(b);SFREE(zeta);SFREE(bj);SFREE(cd);SFREE(cv);
  SFREE(xforcact);SFREE(xloadact);SFREE(xbounact);SFREE(aa);SFREE(bb);SFREE(aanew);
  SFREE(ampli);SFREE(xbodyact);SFREE(bjp);SFREE(bp);SFREE(aamech);SFREE(ikactmech);
  SFREE(xforcdiff);SFREE(xloaddiff);SFREE(xboundiff),SFREE(xbodydiff);

  if(*ithermal==1) {SFREE(t1act);SFREE(t1diff);}

  if(iprescribedboundary){
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
    SFREE(bact);SFREE(bmin);SFREE(bv);SFREE(bprev);SFREE(bdiff);
  }

  if(!cyclicsymmetry){
    SFREE(ad);SFREE(au);
  }else{
    SFREE(adbe); SFREE(aube);SFREE(icole); SFREE(irowe); SFREE(jqe);SFREE(izdof);
    SFREE(nm);

    *nk/=nsectors;
    *ne/=nsectors;
    *nkon/=nsectors;
    *nboun/=nsectors;
    neq[1]=neq[1]*2/nsectors;

    RENEW(ialset,ITG,nalset_);

    /* restore the infomration in istartset and iendset */

    for(j=0; j<*nset; j++){
      istartset[j]=istartset_[j];
      iendset[j]=iendset_[j];
    }
    SFREE(istartset_);
    SFREE(iendset_);

    RENEW(co,double,3**nk);
    if((*ithermal!=0)&&(*nam>0)) RENEW(iamt1,ITG,*nk);
    RENEW(nactdof,ITG,mt**nk);
    if(*ntrans>0) RENEW(inotr,ITG,2**nk);
    RENEW(kon,ITG,*nkon);
    RENEW(ipkon,ITG,*ne);
    RENEW(lakon,char,8**ne);
    RENEW(ielmat,ITG,mi[2]**ne);
    if(*norien>0) RENEW(ielorien,ITG,mi[2]**ne);
    RENEW(nodeboun,ITG,*nboun);
    RENEW(ndirboun,ITG,*nboun);
    if(*nam>0) RENEW(iamboun,ITG,*nboun);
    RENEW(xboun,double,*nboun);
    RENEW(xbounold,double,*nboun);
    RENEW(ikboun,ITG,*nboun);
    RENEW(ilboun,ITG,*nboun);

    /* recovering the original multiple point constraints */

    RENEW(ipompc,ITG,*nmpc);
    RENEW(nodempc,ITG,3**mpcend);
    RENEW(coefmpc,double,*mpcend);
    RENEW(labmpc,char,20**nmpc+1);
    RENEW(ikmpc,ITG,*nmpc);
    RENEW(ilmpc,ITG,*nmpc);
    RENEW(fmpc,double,*nmpc);

    *nmpc=nmpcold;
    *mpcend=mpcendold;
    for(i=0;i<*nmpc;i++){ipompc[i]=ipompcold[i];}
    for(i=0;i<3**mpcend;i++){nodempc[i]=nodempcold[i];}
    for(i=0;i<*mpcend;i++){coefmpc[i]=coefmpcold[i];}
    for(i=0;i<20**nmpc;i++){labmpc[i]=labmpcold[i];}
    for(i=0;i<*nmpc;i++){ikmpc[i]=ikmpcold[i];}
    for(i=0;i<*nmpc;i++){ilmpc[i]=ilmpcold[i];}
    SFREE(ipompcold);SFREE(nodempcold);SFREE(coefmpcold);
    SFREE(labmpcold);SFREE(ikmpcold);SFREE(ilmpcold);

    RENEW(vold,double,mt**nk);
    RENEW(veold,double,mt**nk);
    RENEW(eme,double,6*mi[0]**ne);
    RENEW(sti,double,6*mi[0]**ne);
    if(*nener==1)RENEW(ener,double,2*mi[0]**ne);

    /* distributed loads */

    for(i=0;i<*nload;i++){
      if(nelemload[2*i+1]<nsectors){
	nelemload[2*i]-=*ne*nelemload[2*i+1];
      }else{
	nelemload[2*i]-=*ne*(nelemload[2*i+1]-nsectors);
      }
    }

    /*  sorting the elements with distributed loads */

    if(*nload>0){
      if(*nam>0){
	FORTRAN(isortiiddc,(nelemload,iamload,xload,xloadold,sideload,nload,&kflag));
      }else{
	FORTRAN(isortiddc,(nelemload,xload,xloadold,sideload,nload,&kflag));
      }
    }
      
    /* point loads */
      
    for(i=0;i<*nforc;i++){
      if(nodeforc[2*i+1]<nsectors){
	nodeforc[2*i]-=*nk*nodeforc[2*i+1];
      }else{
	nodeforc[2*i]-=*nk*(nodeforc[2*i+1]-nsectors);
      }
    }
  }

  // if(*nbody>0) SFREE(ipobody);

  if(dashpot){
    SFREE(xini);SFREE(rwork);SFREE(adc);SFREE(auc);SFREE(cc);
    SFREE(rpar);SFREE(iwork);}

  SFREE(cstr);

  SFREE(imddof);SFREE(imdnode);SFREE(imdboun);SFREE(imdmpc);SFREE(imdelem);

  if(iabsload==2) SFREE(fextold);

  *ialsetp=ialset;
  *cop=co;*konp=kon;*ipkonp=ipkon;*lakonp=lakon;*ielmatp=ielmat;
  *ielorienp=ielorien;*inotrp=inotr;*nodebounp=nodeboun;
  *ndirbounp=ndirboun;*iambounp=iamboun;*xbounp=xboun;
  *xbounoldp=xbounold;*ikbounp=ikboun;*ilbounp=ilboun;*nactdofp=nactdof;
  *voldp=vold;*emep=eme;*enerp=ener;*ipompcp=ipompc;*nodempcp=nodempc;
  *coefmpcp=coefmpc;*labmpcp=labmpc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
  *fmpcp=fmpc;*veoldp=veold;*iamt1p=iamt1;*t0p=t0;*t1oldp=t1old;*t1p=t1;
  *stip=sti;*xstatep=xstate;

  (*ttime)+=(*tper);

  return;
}

