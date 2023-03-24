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

void steadystate(double **cop,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,ITG *ne,
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
		 double **voldp,ITG *iperturb,double *sti,ITG *nzs,
		 double *timepar,double *xmodal,
		 double **veoldp,char *amname,double *amta,
		 ITG *namta,ITG *nam,ITG *iamforc,ITG *iamload,
		 ITG **iamt1p,ITG *jout,
		 ITG *kode,char *filab,double **emep,double *xforcold,
		 double *xloadold,
		 double **t1oldp,ITG **iambounp,double **xbounoldp,ITG *iexpl,
		 double *plicon,ITG *nplicon,double *plkcon,ITG *nplkcon,
		 double *xstate,ITG *npmat_,char *matname,ITG *mi,
		 ITG *ncmat_,ITG *nstate_,double **enerp,char *jobnamec,
		 double *ttime,char *set,ITG *nset,ITG *istartset,
		 ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
		 char *prset,ITG *nener,double *trab,
		 ITG **inotrp,ITG *ntrans,double **fmpcp,ITG *ipobody,
		 ITG *ibody,
		 double *xbody,ITG *nbody,double *xbodyold,ITG *istep,
		 ITG *isolver,ITG *jq,char *output,ITG *mcs,ITG *nkon,
		 ITG *ics,double *cs,ITG *mpcend,double *ctrl,
		 ITG *ikforc,ITG *ilforc,double *thicke,ITG *nmat,
		 char *typeboun,ITG *ielprop,double *prop,char *orname,
		 ITG *ndamp,double *dacon,double *t0g,double *t1g){

  char fneig[132]="",description[13]="            ",*lakon=NULL,*labmpc=NULL,
    *labmpcold=NULL,cflag[1]=" ",*labmpc2=NULL;

  ITG nev,i,j,k,*inum=NULL,nsectors,im,intscheme=0,
    iinc=0,l,iout,ielas,icmd=0,iprescribedboundary,ndata,nmd,nevd,
    ndatatot,*iphaseforc=NULL,*iphaseload=NULL,*iphaseboun=NULL,
    *isave=NULL,nfour,ii,ir,ic,mode,noddiam=-1,*nm=NULL,*islavact=NULL,
    *kon=NULL,*ipkon=NULL,*ielmat=NULL,*ielorien=NULL,*inotr=NULL,
    *nodeboun=NULL,*ndirboun=NULL,*iamboun=NULL,*ikboun=NULL,jj,
    *ilboun=NULL,*nactdof=NULL,*ipompc=NULL,*nodempc=NULL,*ikmpc=NULL,
    *ilmpc=NULL,*ipompcold=NULL,*nodempcold=NULL,*ikmpcold=NULL,
    *ilmpcold=NULL,nmpcold,mpcendold,kflag=2,*iamt1=NULL,ifreebody,
    *itg=NULL,ntg=0,symmetryflag=0,inputformat=0,dampmat=0,nrhs=1,
    *ipiv=NULL,info,nev2,ngraph=1,nkg,neg,iflag=1,idummy=1,imax,
    nzse[3],mt=mi[1]+1,*ikactmech=NULL,nactmech,id,nasym=0,
    *imddof=NULL,nmddof,*imdnode=NULL,nmdnode,*imdboun=NULL,nmdboun,
    *imdmpc=NULL,nmdmpc,*izdof=NULL,nzdof,cyclicsymmetry,mortar=0,
    *ikactmechr=NULL,*ikactmechi=NULL,nactmechr,nactmechi,intpointvar,
    iforc,iload,ne0,*iponoel=NULL,*inoel=NULL,*imdelem=NULL,
    nmdelem,*integerglob=NULL,*nshcon=NULL,nherm,icfd=0,*inomat=NULL,
    *islavnode=NULL,*nslavnode=NULL,*islavsurf=NULL,iit=-1,
    network=0,kscale=0,nmethodact=1,iperturbsav,coriolis,*itiefac=NULL,
    mscalmethod=0,*islavelinv=NULL,*irowtloc=NULL,*jqtloc=NULL,nboun2,
    *ndirboun2=NULL,*nodeboun2=NULL,nmpc2,*ipompc2=NULL,*nodempc2=NULL,
    *ikboun2=NULL,*ilboun2=NULL,*ikmpc2=NULL,*ilmpc2=NULL,mortartrafoflag=0;

  long long i2;

  double *d=NULL,*z=NULL,*stiini=NULL,*vini=NULL,*freqnh=NULL,
    *xforcact=NULL,*xloadact=NULL,y,*fr=NULL,*fi=NULL,*cc=NULL,
    *t1act=NULL,*ampli=NULL,*aa=NULL,*bb=NULL,*vr=NULL,*vi=NULL,
    *stn=NULL,*stx=NULL,*een=NULL,*adb=NULL,*xstiff=NULL,ptime,
    *aub=NULL,*bjr=NULL,*bji=NULL,*xbodyr=NULL,*cdn=NULL,
    *f=NULL,*fn=NULL,*xbounact=NULL,*epn=NULL,*xstateini=NULL,
    *enern=NULL,*xstaten=NULL,*eei=NULL,*enerini=NULL,*qfn=NULL,
    *qfx=NULL,*xbodyact=NULL,*cgr=NULL,*au=NULL,*xbodyi=NULL,
    time,dtime,reltime=1.,*co=NULL,*xboun=NULL,*xbounold=NULL,
    physcon[1],qa[4],cam[5],accold[1],bet,gam,*emn=NULL,timem,
    *ad=NULL,sigma=0.,alpham,betam,*fnr=NULL,*fni=NULL,*emeini=NULL,
    fmin,fmax,bias,*freq=NULL,*xforcr=NULL,dd,pi,vreal,constant,
    *xforci=NULL,*xloadr=NULL,*xloadi=NULL,*xbounr=NULL,*xbouni=NULL,
    *br=NULL,*bi=NULL,*ubr=NULL,*ubi=NULL,*mubr=NULL,*mubi=NULL,
    *wsave=NULL,*r=NULL,*xbounacttime=NULL,*btot=NULL,breal,tmin,tmax,
    *vold=NULL,*eme=NULL,*ener=NULL,*coefmpc=NULL,*fmpc=NULL,
    *coefmpcold=NULL,*t0=NULL,*t1=NULL,*t1old=NULL,*adc=NULL,*auc=NULL,
    *am=NULL,*bm=NULL,*zc=NULL,*e=NULL,*stnr=NULL,*stni=NULL,
    *vmax=NULL,*stnmax=NULL,*va=NULL,*vp=NULL,*fric=NULL,*springarea=NULL,
    *stna=NULL,*stnp=NULL,*bp=NULL,*eenmax=NULL,*clearini=NULL,
    *doubleglob=NULL,*shcon=NULL,*veold=NULL,*xmr=NULL,*xmi=NULL,*eig=NULL,
    *ax=NULL,*bx=NULL,*pslavsurf=NULL,*pmastsurf=NULL,xnull=0.,
    *cdnr=NULL,*cdni=NULL,*energyini=NULL,*energy=NULL,*v=NULL,*b=NULL,
    *cco=NULL,*smscale=NULL,*autloc=NULL,*xboun2=NULL,*coefmpc2=NULL,
    *umubr=NULL,*umubi=NULL,*fnext=NULL;

  /* dummy arguments for the call of expand*/

  char* tieset=NULL;
  ITG *jqe=NULL,*icole=NULL,*irowe=NULL,ntie=0;
  double *adbe=NULL,*aube=NULL;

  FILE *f1;

  ITG *ipneigh=NULL,*neigh=NULL;

#ifdef SGI
  ITG token;
#endif

  co=*cop;kon=*konp;ipkon=*ipkonp;lakon=*lakonp;ielmat=*ielmatp;
  ielorien=*ielorienp;inotr=*inotrp;nodeboun=*nodebounp;
  ndirboun=*ndirbounp;iamboun=*iambounp;xboun=*xbounp;veold=*veoldp;
  xbounold=*xbounoldp;ikboun=*ikbounp;ilboun=*ilbounp;nactdof=*nactdofp;
  vold=*voldp;eme=*emep;ener=*enerp;ipompc=*ipompcp;nodempc=*nodempcp;
  coefmpc=*coefmpcp;labmpc=*labmpcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
  fmpc=*fmpcp;iamt1=*iamt1p;t0=*t0p;t1=*t1p;t1old=*t1oldp;

  pi=4.*atan(1.);
  iout=2;

  alpham=xmodal[0];
  betam=xmodal[1];

  fmin=2.*pi*xmodal[2];
  fmax=2.*pi*xmodal[3];
  ndata=floor(xmodal[4]);
  bias=xmodal[5];
  nfour=floor(xmodal[6]);
  if(nfour>0){
    tmin=xmodal[7];
    tmax=xmodal[8];
  }

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
    printf(" *ERROR in steadystate: cannot open eigenvalue file for reading");
    exit(0);
  }

  if(fread(&cyclicsymmetry,sizeof(ITG),1,f1)!=1){
    printf(" *ERROR in steadystate reading the cyclic symmetry flag in the eigenvalue file");
    printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
    printf("        1) the nonexistence of the .eig file\n");
    printf("        2) other boundary conditions than in the input deck\n");
    printf("           which created the .eig file\n\n");
    exit(0);
  }

  if(fread(&nherm,sizeof(ITG),1,f1)!=1){
    printf(" *ERROR in steadystate reading the Hermitian flag in the eigenvalue file");
    printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
    printf("        1) the nonexistence of the .eig file\n");
    printf("        2) other boundary conditions than in the input deck\n");
    printf("           which created the .eig file\n\n");
    exit(0);
  }

  iperturbsav=iperturb[0];
  if(fread(&iperturb[0],sizeof(ITG),1,f1)!=1){
    printf(" *ERROR in steadystate reading the perturbation flag in the eigenvalue file");
    printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
    printf("        1) the nonexistence of the .eig file\n");
    printf("        2) other boundary conditions than in the input deck\n");
    printf("           which created the .eig file\n\n");
    exit(0);
  }
  if(iperturbsav!=iperturb[0]){
    if(iperturb[0]==0){
      printf(" *ERROR in steadystate: the .eig-file was created without perturbation info\n");
      printf("        the present *STEADY STATE DYNAMICS step, however, contains the\n");
      printf("        perturbation parameter on the *STEP card.\n\n");
    }else if(iperturb[0]==1){
      printf(" *ERROR in steadystate: the .eig-file was created with perturbation info\n\n");
      printf("        the present *STEADY STATE DYNAMICS step, however, does not contain the\n");
      printf("        perturbation parameter on the *STEP card.\n");
    }
    FORTRAN(stop,());
  }

  if(iperturb[0]==1){
    if(fread(vold,sizeof(double),mt**nk,f1)!=mt**nk){
      printf(" *ERROR in steadystate reading the reference displacements in the eigenvalue file...");
      printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
      printf("        1) the nonexistence of the .eig file\n");
      printf("        2) other boundary conditions than in the input deck\n");
      printf("           which created the .eig file\n\n");
      exit(0);
    }
  }

  /* determining retained nodes, dofs, spcs and mpcs based on:
     1) the .dat output
     2) the .frd output
     3) contact definitions (master/slave)
     4) nonlinear mpc's  */

  nmddof=0;nmdnode=0;nmdboun=0;nmdmpc=0;nmdelem=0;

  NNEW(imddof,ITG,*nk*3);
  NNEW(imdnode,ITG,*nk);
  NNEW(imdboun,ITG,*nboun);
  NNEW(imdmpc,ITG,*nmpc);
  FORTRAN(createmddof,(imddof,&nmddof,istartset,iendset,
		       ialset,nactdof,ithermal,mi,imdnode,&nmdnode,
		       ikmpc,ilmpc,ipompc,nodempc,nmpc,
		       imdmpc,&nmdmpc,imdboun,&nmdboun,ikboun,
		       nboun,nset,&ntie,tieset,set,lakon,kon,ipkon,labmpc,
		       ilboun,filab,prlab,prset,nprint,ne,&cyclicsymmetry));

  /* if results are requested in too many nodes, it is faster to 
     calculate the results in all nodes */

  if((nmdnode>*nk/2)&&(!cyclicsymmetry)){
    nmdnode=0;nmddof=0;nmdboun=0;nmdmpc=0;
  }
  
  /* determining retained nodes, dofs, spcs and mpcs based
     on USER defined cloads and dloads */

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

  if(!cyclicsymmetry){

    nkg=*nk;
    neg=*ne;

    if(fread(&nev,sizeof(ITG),1,f1)!=1){
      printf(" *ERROR in steadystate reading the number of eigenvalues in the eigenvalue file");
      printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
      printf("        1) the nonexistence of the .eig file\n");
      printf("        2) other boundary conditions than in the input deck\n");
      printf("           which created the .eig file\n\n");
      exit(0);
    }

    if(nherm==1){
      NNEW(d,double,nev);
      if(fread(d,sizeof(double),nev,f1)!=nev){
	printf(" *ERROR in steadystate reading the eigenvalues in the eigenvalue file...");
	printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
	printf("        1) the nonexistence of the .eig file\n");
	printf("        2) other boundary conditions than in the input deck\n");
	printf("           which created the .eig file\n\n");
	exit(0);
      }

      for(i=0;i<nev;i++){
	if(d[i]<0){d[i]=0.;}
      }
    }else{
      NNEW(d,double,2*nev);
      if(fread(d,sizeof(double),2*nev,f1)!=2*nev){
	printf(" *ERROR in steadystate reading the eigenvalues in the eigenvalue file...");
	printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
	printf("        1) the nonexistence of the .eig file\n");
	printf("        2) other boundary conditions than in the input deck\n");
	printf("           which created the .eig file\n\n");
	exit(0);
      }
    }
      
    NNEW(ad,double,neq[1]);
    NNEW(adb,double,neq[1]);
    NNEW(au,double,nzs[2]);
    NNEW(aub,double,nzs[1]);
      
    /* reading the stiffness matrix */
      
    if(fread(ad,sizeof(double),neq[1],f1)!=neq[1]){
      printf(" *ERROR in steadystate reading the diagonal of the stiffness matrix in the eigenvalue file");
      printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
      printf("        1) the nonexistence of the .eig file\n");
      printf("        2) other boundary conditions than in the input deck\n");
      printf("           which created the .eig file\n\n");
      exit(0);
    }
      
    if(fread(au,sizeof(double),nzs[2],f1)!=nzs[2]){
      printf(" *ERROR in steadystate reading the off-diagonals of the stiffness matrix in the eigenvalue file");
      printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
      printf("        1) the nonexistence of the .eig file\n");
      printf("        2) other boundary conditions than in the input deck\n");
      printf("           which created the .eig file\n\n");
      exit(0);
    }
      
    /* reading the mass matrix */
      
    if(fread(adb,sizeof(double),neq[1],f1)!=neq[1]){
      printf(" *ERROR in steadystate reading the diagonal of the mass matrix in the eigenvalue file");
      printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
      printf("        1) the nonexistence of the .eig file\n");
      printf("        2) other boundary conditions than in the input deck\n");
      printf("           which created the .eig file\n\n");
      exit(0);
    }
      
    if(fread(aub,sizeof(double),nzs[1],f1)!=nzs[1]){
      printf(" *ERROR in steadystate reading the off-diagonals of the mass matrix in the eigenvalue file");
      printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
      printf("        1) the nonexistence of the .eig file\n");
      printf("        2) other boundary conditions than in the input deck\n");
      printf("           which created the .eig file\n\n");
      exit(0);
    }
      
    /* reading the eigenvectors */

    if(nherm==1){
      NNEW(z,double,neq[1]*nev);
      if(fread(z,sizeof(double),neq[1]*nev,f1)!=neq[1]*nev){
	printf(" *ERROR in complexfreq reading the eigenvectors in the eigenvalue file...");
	printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
	printf("        1) the nonexistence of the .eig file\n");
	printf("        2) other boundary conditions than in the input deck\n");
	printf("           which created the .eig file\n\n");
	exit(0);
      }
    }else{
      NNEW(z,double,2*neq[1]*nev);
      if(fread(z,sizeof(double),2*neq[1]*nev,f1)!=2*neq[1]*nev){
	printf(" *ERROR in complexfreq reading the eigenvectors in the eigenvalue file...");
	printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
	printf("        1) the nonexistence of the .eig file\n");
	printf("        2) other boundary conditions than in the input deck\n");
	printf("           which created the .eig file\n\n");
	exit(0);
      }
    }

    /* reading the orthogonality matrices */

    if(nherm!=1){
      NNEW(xmr,double,nev*nev);
      NNEW(xmi,double,nev*nev);
      if(fread(xmr,sizeof(double),nev*nev,f1)!=nev*nev){
	printf(" *ERROR in steadystate reading the real orthogonality matrix to the eigenvalue file...");
	printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
	printf("        1) the nonexistence of the .eig file\n");
	printf("        2) other boundary conditions than in the input deck\n");
	printf("           which created the .eig file\n\n");
	exit(0);
      }
	  
      if(fread(xmi,sizeof(double),nev*nev,f1)!=nev*nev){
	printf(" *ERROR in steadystate reading the imaginary orthogonality matrix to the eigenvalue file...");
	printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
	printf("        1) the nonexistence of the .eig file\n");
	printf("        2) other boundary conditions than in the input deck\n");
	printf("           which created the .eig file\n\n");
	exit(0);
      }
    }
  }
  else{
    nev=0;
    do{
      if(fread(&nmd,sizeof(ITG),1,f1)!=1){
	break;
      }

      if(fread(&nevd,sizeof(ITG),1,f1)!=1){
	printf(" *ERROR in steadystate reading the number of eigenvalues for nodal diameter %" ITGFORMAT " in the eigenvalue file",nmd);
	printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
	printf("        1) the nonexistence of the .eig file\n");
	printf("        2) other boundary conditions than in the input deck\n");
	printf("           which created the .eig file\n\n");
	exit(0);
      }

      /* reading the eigenvalues (complex for non-Hermitian systems) */

      if(nev==0){
	if(nherm==1){NNEW(d,double,nevd);
	}else{NNEW(d,double,2*nevd);}
	NNEW(nm,ITG,nevd);
      }else{
	if(nherm!=1){
	  printf(" *ERROR in steadystate: non-Hermitian systems cannot\n");
	  printf("        be combined with multiple modal diameters\n");
	  printf("        in cyclic symmetry calculations\n\n");
	  FORTRAN(stop,());
	}
	RENEW(d,double,nev+nevd);
	RENEW(nm,ITG,nev+nevd);
      }
	  
      if(nherm==1){
	if(fread(&d[nev],sizeof(double),nevd,f1)!=nevd){
	  printf(" *ERROR in steadystate reading the eigenvalues for nodal diameter %" ITGFORMAT " in the eigenvalue file",nmd);
	  printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
	  printf("        1) the nonexistence of the .eig file\n");
	  printf("        2) other boundary conditions than in the input deck\n");
	  printf("           which created the .eig file\n\n");
	  exit(0);
	}
	for(i=nev;i<nev+nevd;i++){
	  if(d[i]<0){d[i]=0.;}
	}
      }else{
	if(fread(&d[nev],sizeof(double),2*nevd,f1)!=2*nevd){
	  printf(" *ERROR in steadystate reading the eigenvalues in the eigenvalue file...");
	  printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
	  printf("        1) the nonexistence of the .eig file\n");
	  printf("        2) other boundary conditions than in the input deck\n");
	  printf("           which created the .eig file\n\n");
	  exit(0);
	}
      }

      for(i=nev;i<nev+nevd;i++){nm[i]=nmd;}
	  
      if(nev==0){
	NNEW(ad,double,neq[1]);
	NNEW(au,double,nzs[1]);

	if(fread(ad,sizeof(double),neq[1],f1)!=neq[1]){
	  printf(" *ERROR in steadystate reading the diagonal of the stiffness matrix in the eigenvalue file");
	  printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
	  printf("        1) the nonexistence of the .eig file\n");
	  printf("        2) other boundary conditions than in the input deck\n");
	  printf("           which created the .eig file\n\n");
	  exit(0);
	}
	      
	if(fread(au,sizeof(double),nzs[1],f1)!=nzs[1]){
	  printf(" *ERROR in steadystate reading the off-diagonals of the stiffness matrix in the eigenvalue file");
	  printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
	  printf("        1) the nonexistence of the .eig file\n");
	  printf("        2) other boundary conditions than in the input deck\n");
	  printf("           which created the .eig file\n\n");
	  exit(0);
	}

	SFREE(ad);SFREE(au);
	      
	NNEW(adb,double,neq[1]);
	NNEW(aub,double,nzs[1]);

	if(fread(adb,sizeof(double),neq[1],f1)!=neq[1]){
	  printf(" *ERROR in steadystate reading the diagonal of the mass matrix in the eigenvalue file");
	  printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
	  printf("        1) the nonexistence of the .eig file\n");
	  printf("        2) other boundary conditions than in the input deck\n");
	  printf("           which created the .eig file\n\n");
	  exit(0);
	}
	      
	if(fread(aub,sizeof(double),nzs[1],f1)!=nzs[1]){
	  printf(" *ERROR in steadystate reading the off-diagonals of the mass matrix in the eigenvalue file");
	  printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
	  printf("        1) the nonexistence of the .eig file\n");
	  printf("        2) other boundary conditions than in the input deck\n");
	  printf("           which created the .eig file\n\n");
	  exit(0);
	}
      }
	  
      /* reading the eigenvectors */

      if(nev==0){
	NNEW(z,double,neq[1]*nevd);
      }else{
	RENEW(z,double,(long long)neq[1]*(nev+nevd));
      }
	  
      if(fread(&z[(long long)neq[1]*nev],sizeof(double),neq[1]*nevd,f1)!=neq[1]*nevd){
	printf(" *ERROR in steadystate reading the eigenvectors for nodal diameter %" ITGFORMAT " in the eigenvalue file",nmd);
	printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
	printf("        1) the nonexistence of the .eig file\n");
	printf("        2) other boundary conditions than in the input deck\n");
	printf("           which created the .eig file\n\n");
	exit(0);
      }

      /* reading the orthogonality matrices */
	  
      if(nherm!=1){
	NNEW(xmr,double,nev*nev);
	NNEW(xmi,double,nev*nev);
	if(fread(xmr,sizeof(double),nev*nev,f1)!=nev*nev){
	  printf(" *ERROR in steadystate reading the real orthogonality matrix to the eigenvalue file...");
	  printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
	  printf("        1) the nonexistence of the .eig file\n");
	  printf("        2) other boundary conditions than in the input deck\n");
	  printf("           which created the .eig file\n\n");
	  exit(0);
	}
	      
	if(fread(xmi,sizeof(double),nev*nev,f1)!=nev*nev){
	  printf(" *ERROR in steadystate reading the imaginary orthogonality matrix to the eigenvalue file...");
	  printf(" *INFO  in steadystate: if there are problems reading the .eig file this may be due to:\n");
	  printf("        1) the nonexistence of the .eig file\n");
	  printf("        2) other boundary conditions than in the input deck\n");
	  printf("           which created the .eig file\n\n");
	  exit(0);
	}
      }

      nev+=nevd;
    }while(1);

    /* determining the maximum amount of sectors */

    for(i=0;i<*mcs;i++){
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
	   labmpcold,&nmpcold,xloadold,iamload,t1old,t1,iamt1,xstiff,&icole,&jqe,
	   &irowe,isolver,nzse,&adbe,&aube,iexpl,
	   ibody,xbody,nbody,cocon,ncocon,tieset,&ntie,imddof,&nmddof,
	   imdnode,&nmdnode,imdboun,&nmdboun,imdmpc,&nmdmpc,&izdof,&nzdof,
	   &nherm,xmr,xmi,typeboun,ielprop,prop,orname,itiefac,t0g,t1g);

    RENEW(imddof,ITG,nmddof);
    RENEW(imdnode,ITG,nmdnode);
    RENEW(imdboun,ITG,nmdboun);
    RENEW(imdmpc,ITG,nmdmpc);

    SFREE(veold);
    NNEW(veold,double,mt**nk);
    RENEW(eme,double,6*mi[0]**ne);

    if(*nener==1) RENEW(ener,double,2*mi[0]**ne);
  }

  fclose(f1);

  /* allocating space for the friction coefficients */

  if(nherm==1){
    NNEW(fric,double,nev);
  }else{
    NNEW(fric,double,2*nev);
  }

  /* check whether there is structural damping or whether there are dashpot elements
     both lead to the creation of a damping matrix */

  if(*ndamp>0){

    /* structural damping */

    dampmat=1;
  }else{

    /* check for dashpot elements */

    for(i=0;i<*ne;i++){
      if(ipkon[i]==-1) continue;
      if(strcmp1(&lakon[i*8],"ED")==0){
	dampmat=1;break;}
    }
  }

  if(dampmat){
    if(cyclicsymmetry){
      printf(" *ERROR in steadystate: structural damping and/or dashpot elements\n");
      printf("        are not allowed in combination with cyclic symmetry\n");
      FORTRAN(stop,());
    }

    if(nherm!=1){
      printf("ERROR in steadystate: structural damping and/or dashpot elements \n");
      printf("      cannot be combined with non-Hermitian systems (in the present version of CalculiX)\n");
      FORTRAN(stop,());
    }

    /* cc is the reduced damping matrix (damping matrix mapped onto
       space spanned by eigenmodes) */

    NNEW(cc,double,nev*nev);
  }

  /* check whether there are centrifugal forces; if so, they are converted
     into Coriolis forces */

  coriolis=0;
  for(i=0;i<*nbody;i++){
    if(ibody[3*i]==1){
      ibody[3*i]=4;
      coriolis=1;
    }
  }

  if(coriolis){
    if(nherm!=1){
      printf(" *WARNING in steadystate: centrifugal loading cannot be combined \n");
      printf("          with non-Hermitian systems (in the present version of CalculiX)\n");
      printf("          it is deactivated\n");
      for(i=0;i<*nbody;i++){
	if(ibody[3*i]==4){
	  ibody[3*i]=0;
	}
	coriolis=0;
      }
    }else{

      /* cc is the reduced damping matrix (damping matrix mapped onto
         space spanned by eigenmodes) */

      NNEW(cc,double,nev*nev);
    }
  }

  NNEW(inum,ITG,*nk);
  strcpy1(&cflag[0],&filab[4],1);
  FORTRAN(createinum,(ipkon,inum,kon,lakon,nk,ne,&cflag[0],nelemload,
		      nload,nodeboun,nboun,ndirboun,ithermal,co,vold,mi,ielmat,
                      ielprop,prop));

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
 
  /**********************************************************************

  HARMONIC EXCITATION

  **********************************************************************/
  
  if(nfour<=0){

    /* harmonic excitation */

    NNEW(ikactmechr,ITG,neq[1]);
    NNEW(ikactmechi,ITG,neq[1]);
    nactmechr=0;nactmechi=0;
      
    /* result fields */
      
    if(intpointvar==1){
      NNEW(fn,double,mt**nk);
      NNEW(stnr,double,6**nk);
      NNEW(stni,double,6**nk);
      NNEW(stx,double,6*mi[0]**ne);
      NNEW(eei,double,6*mi[0]**ne);

      if(*ithermal>1) {NNEW(qfn,double,3**nk);
	NNEW(qfx,double,3*mi[0]**ne);}
	  
      if(strcmp1(&filab[261],"E   ")==0) NNEW(een,double,6**nk);
      if(strcmp1(&filab[522],"ENER")==0) NNEW(enern,double,*nk);
      if(strcmp1(&filab[2697],"ME  ")==0) NNEW(emn,double,6**nk);
	  
      if(*nener==1){
	NNEW(stiini,double,6*mi[0]**ne);
	NNEW(emeini,double,6*mi[0]**ne);
	NNEW(enerini,double,2*mi[0]**ne);}
    }
      
    /* determining the frequency data points */
      
    NNEW(freq,double,ndata*(nev+1));
      
    ndatatot=0.;
    freq[0]=fmin;
    if(fabs(fmax-fmin)<1.e-10){
      ndatatot=1;
    }else{

      /* copy the eigenvalues and sort them in ascending order
         (important for values from distinct nodal diameters */

      NNEW(e,double,nev);
      if(nherm==1){
	for(i=0;i<nev;i++){e[i]=sqrt(d[i]);}
      }else{

	/* for complex eigenvalues: sorting the real part */

	for(i=0;i<nev;i++){
	  e[i]=sqrt(sqrt(d[2*i]*d[2*i]+d[2*i+1]*d[2*i+1])+d[2*i])/sqrt(2.);
	}
      }
      FORTRAN(dsort,(e,&idummy,&nev,&iflag));
	  
      for(i=0;i<nev;i++){
	if(i!=0){
	  if(fabs(e[i]-e[i-1])<1.e-5){continue;}
	}
	if(e[i]>=fmin){
	  if(e[i]<=fmax){
	    for(j=1;j<ndata;j++){
	      y=-1.+2.*j/((double)(ndata-1));
	      if(fabs(y)<1.e-10){freq[ndatatot+j]=
		  (freq[ndatatot]+e[i])/2.;}
	      else{
		freq[ndatatot+j]=(freq[ndatatot]+e[i])/2.+
		  (e[i]-freq[ndatatot])*pow(fabs(y),1./bias)*
		  y/(2.*fabs(y));
	      }
	    }
	    ndatatot+=ndata-1;
	  }
	  else{break;}
	}
      }

      SFREE(e);

      for(j=1;j<ndata;j++){
	y=-1.+2.*j/((double)(ndata-1));
	if(fabs(y)<1.e-10){freq[ndatatot+j]=(freq[ndatatot]+fmax)/2.;}
	else{
	  freq[ndatatot+j]=(freq[ndatatot]+fmax)/2.+
	    (fmax-freq[ndatatot])*pow(fabs(y),1./bias)*
	    y/(2.*fabs(y));
	}
      }
      ndatatot+=ndata;
    }
    RENEW(freq,double,ndatatot);
      
    /*  check for nonzero SPC's */
      
    iprescribedboundary=0;
    for(i=0;i<*nboun;i++){
      if(fabs(xboun[i])>1.e-10){
	iprescribedboundary=1;
	if((*isolver==2)||(*isolver==3)){
	  printf(" *ERROR in steadystate: the iterative solver\n");
	  printf("        cannot be used for steady state dynamic\n");
	  printf("        calculations with prescribed boundary\n");
	  printf("        conditions\n");
	  FORTRAN(stop,());
	}
	nmdnode=0;nmddof=0;nmdboun=0;nmdmpc=0;
	break;
      }
    }

    if((iprescribedboundary)&&(cyclicsymmetry)){
      printf(" *ERROR in steadystate: prescribed boundaries are not allowed in combination with cyclic symmetry\n");
      FORTRAN(stop,());
    }

    /* calculating the damping coefficients = friction coefficient*2*eigenvalue 
       this is done for Rayleigh damping and direct damping
       (both activated with the *MODAL DAMPING card) */
      
    if(xmodal[10]<0){

      /* Rayleigh damping */

      for(i=0;i<nev;i++){
	if(nherm==1){
	  if(sqrt(d[i])>(1.e-10)){
	    fric[i]=(alpham+betam*d[i]);
	  }
	  else {
	    printf(" *WARNING in steadystate: one of the frequencies is zero\n");
	    printf("         no Rayleigh mass damping allowed\n");
	    fric[i]=0.;
	  }
	}else{
	  fric[2*i]=sqrt(sqrt(d[2*i]*d[2*i]+d[2*i+1]*d[2*i+1])+d[2*i])/
	    sqrt(2.);
	  fric[2*i+1]=sqrt(sqrt(d[2*i]*d[2*i]+d[2*i+1]*d[2*i+1])-d[2*i])/
	    sqrt(2.);
	  if(d[2*i+1]<0.) fric[2*i+1]=-fric[2*i+1];
	  fric[2*i]=alpham+betam*fric[2*i];
	  fric[2*i+1]=betam*fric[2*i+1];
	}
      }
    }
    else if((ITG)xmodal[10]>0){

      /* direct damping */

      /*      if(iprescribedboundary){
	printf(" *ERROR in steadystate: prescribed boundaries are not allowed in combination with direct modal damping\n");
	FORTRAN(stop,());
	}*/
	      
      /*copy the damping coefficients for every eigenfrequencie from xmodal[11....] */
      if(nev<(ITG)xmodal[10]){
	imax=nev;
	printf(" *WARNING in steadystate: too many modal damping coefficients applied\n");
	printf("         damping coefficients corresponding to nonexisting eigenvalues are ignored\n");
      }
      else{
	imax=(ITG)xmodal[10];
      }
      for(i=0; i<imax; i++){
	if(nherm==1){
	  fric[i]=2.*sqrt(d[i])*xmodal[11+i];    
	}else{
	  fric[2*i]=sqrt(sqrt(d[2*i]*d[2*i]+d[2*i+1]*d[2*i+1])+d[2*i])/
	    sqrt(2.);
	  fric[2*i+1]=sqrt(sqrt(d[2*i]*d[2*i]+d[2*i+1]*d[2*i+1])-d[2*i])/
	    sqrt(2.);
	  if(d[2*i+1]<0.) fric[2*i+1]=-fric[2*i+1];
	  fric[2*i]=2.*fric[2*i]*xmodal[11+i];
	  fric[2*i+1]=2.*fric[2*i+1]*xmodal[11+i];
	} 
      }
	  
    }
      
    /* check whether the loading is real or imaginary */
      
    NNEW(iphaseforc,ITG,*nforc);
    for (i=0;i<*nforc;i++){
      if(nodeforc[2*i+1]>=nsectors){
	iphaseforc[i]=1;
      }
    }
      
    NNEW(iphaseload,ITG,*nload);
    for (i=0;i<*nload;i++){
      if(nelemload[2*i+1]>=nsectors){
	iphaseload[i]=1;
      }
    }
      
    if(iprescribedboundary){
      NNEW(iphaseboun,ITG,*nboun);
      for (i=0;i<*nboun;i++){
	if(nodeboun[i]>*nk){
	  iphaseboun[i]=1;
	  nodeboun[i]=nodeboun[i]-*nk;
	}
      }
    }
      
    /* allocating actual loading fields */
      
    NNEW(xforcact,double,*nforc);
    NNEW(xforcr,double,*nforc);
    NNEW(xforci,double,*nforc);
      
    NNEW(xloadact,double,2**nload);
    NNEW(xloadr,double,2**nload);
    NNEW(xloadi,double,2**nload);
      
    NNEW(xbodyact,double,7**nbody);
    NNEW(xbodyr,double,7**nbody);
    NNEW(xbodyi,double,7**nbody);
    /* copying the rotation axis and/or acceleration vector */
    for(k=0;k<7**nbody;k++){xbodyact[k]=xbody[k];}
      
    NNEW(xbounact,double,*nboun);
      
    if(*ithermal==1) NNEW(t1act,double,*nk);
      
    /* assigning the body forces to the elements */ 

    /*   if(*nbody>0){
      ifreebody=*ne+1;
      NNEW(ipobody,ITG,2*ifreebody**nbody);
      for(k=1;k<=*nbody;k++){
	FORTRAN(bodyforce,(cbody,ibody,ipobody,nbody,set,istartset,
			   iendset,ialset,&inewton,nset,&ifreebody,&k));
	RENEW(ipobody,ITG,2*(*ne+ifreebody));
      }
      RENEW(ipobody,ITG,2*(ifreebody-1));
      }*/
      
    NNEW(br,double,neq[1]); /* load rhs vector */
    NNEW(bi,double,neq[1]); /* load rhs vector */
      
    if(iprescribedboundary){
      NNEW(xbounr,double,*nboun);
      NNEW(xbouni,double,*nboun);
	  
      NNEW(fr,double,neq[1]); /* force corresponding to real particular solution */
      NNEW(fi,double,neq[1]); /* force corresponding to imaginary particular solution */
	  
      NNEW(ubr,double,neq[1]); /* real particular solution */
      NNEW(ubi,double,neq[1]); /* imaginary particular solution */
	  
      NNEW(mubr,double,neq[1]); /* mass times real particular solution */
      NNEW(mubi,double,neq[1]); /* mass times imaginary particular solution */
    }
      
    NNEW(bjr,double,nev); /* real response modal decomposition */
    NNEW(bji,double,nev); /* imaginary response modal decomposition */
      
    NNEW(ampli,double,*nam); /* instantaneous amplitude */
      
    if(nherm==1){
      NNEW(aa,double,nev); /* modal coefficients of the real loading */
      NNEW(bb,double,nev); /* modal coefficients of the imaginary loading */
    }else{
      NNEW(aa,double,2*nev); /* modal coefficients of the real loading */
      NNEW(bb,double,2*nev); /* modal coefficients of the imaginary loading */
    }
    
    /* result fields */
      
    NNEW(vr,double,mt**nk);
    NNEW(vi,double,mt**nk);
      
    if(iprescribedboundary){
	  
      /* LU decomposition of the stiffness matrix */
	  
      if(*isolver==0){
#ifdef SPOOLES
	spooles_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
		       &symmetryflag,&inputformat,&nzs[2]);
#else
	printf(" *ERROR in steadystate: the SPOOLES library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==4){
#ifdef SGI
	token=1;
	sgi_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],token);
#else
	printf(" *ERROR in steadystate: the SGI library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
	tau_factor(ad,&au,adb,aub,&sigma,icol,&irow,&neq[1],&nzs[1]);
#else
	printf(" *ERROR in steadystate: the TAUCS library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
		       &symmetryflag,&inputformat,jq,&nzs[2]);
#else
	printf(" *ERROR in steadystate: the PARDISO library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==8){
#ifdef PASTIX
	pastix_factor_main(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
		      &symmetryflag,&inputformat,jq,&nzs[2]);
#else
	printf(" *ERROR in steadystate: the PASTIX library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
    }
      
    for(l=0;l<ndatatot;l=l+*jout){
      for(i=0;i<6*mi[0]**ne;i++){eme[i]=0.;}
      time=freq[l]/(2.*pi);
      timem=-time;
      ptime=time;

      /* calculating cc */

      if(dampmat||coriolis){

	/* determining the elastic constants in xstiff (needed in mafilldmss) */

	iout=-1;
	NNEW(xstiff,double,(long long)27*mi[0]**ne);
	NNEW(v,double,mt**nk);
	NNEW(f,double,*neq);
	if(intpointvar!=1){
	  NNEW(fn,double,mt**nk);
	  NNEW(stx,double,6*mi[0]**ne);
	}
	results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
		elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
		prestr,iprestr,filab,eme,emn,een,iperturb,
		f,fn,nactdof,&iout,qa,vold,b,nodeboun,
		ndirboun,xbounact,nboun,ipompc,
		nodempc,coefmpc,labmpc,nmpc,&nmethodact,cam,neq,veold,accold,
		&bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
		&icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
		emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
		iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
		fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
		&reltime,&ne0,thicke,shcon,nshcon,
		sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		&mortar,islavact,cdn,islavnode,nslavnode,&ntie,clearini,
		islavsurf,ielprop,prop,energyini,energy,&kscale,iponoel,
		inoel,nener,orname,&network,ipobody,xbodyact,ibody,typeboun,
		itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
		islavelinv,autloc,irowtloc,jqtloc,&nboun2,
		ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
		labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
		&intscheme,physcon);
	SFREE(v);SFREE(f);
	if(intpointvar!=1){SFREE(fn);SFREE(stx);}
	iout=2;

	for(i=0;i<nev*nev;i++){cc[i]=0.;}

	if(dampmat){
	      
	  NNEW(adc,double,neq[1]);
	  NNEW(auc,double,nzs[1]);
	  mafilldmssmain(co,nk,kon,ipkon,lakon,ne,
			 ipompc,nodempc,coefmpc,nmpc,
			 nelemload,sideload,xload,nload,xbody,ipobody,
			 nbody,cgr,adc,auc,nactdof,jq,irow,neq,
			 nmethod,ikmpc,ilmpc,
			 elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
			 ielmat,ielorien,norien,orab,ntmat_,
			 t0,t0,ithermal,vold,iperturb,sti,
			 nzs,stx,iexpl,plicon,nplicon,plkcon,nplkcon,
			 xstiff,npmat_,&dtime,matname,mi,
			 ncmat_,physcon,ttime,&time,istep,&iinc,
			 ibody,xloadold,&reltime,veold,springarea,nstate_,
			 xstateini,xstate,thicke,integerglob,doubleglob,
			 tieset,istartset,iendset,ialset,&ntie,&nasym,pslavsurf,
			 pmastsurf,&mortar,clearini,ielprop,prop,&ne0,
			 &freq[l],ndamp,dacon,set,nset);

	  /*  zc = damping matrix * eigenmodes */
		  
	  NNEW(zc,double,neq[1]*nev);
	  for(i=0;i<nev;i++){
	    FORTRAN(op,(&neq[1],&z[(long long)i*neq[1]],&zc[i*neq[1]],adc,auc,
			jq,irow));
	  }
	  SFREE(adc);SFREE(auc);
	      
	  /* cc is the reduced damping matrix (damping matrix mapped onto
	     space spanned by eigenmodes) */
		  
	  for(i=0;i<nev;i++){
	    for(j=0;j<=i;j++){
	      for(k=0;k<neq[1];k++){
		cc[i*nev+j]+=z[(long long)j*neq[1]+k]*zc[i*neq[1]+k];
	      }
	    }
	  }
	  SFREE(zc);
	      
	  /* symmetric part of cc matrix */
	      
	  for(i=0;i<nev;i++){
	    for(j=i;j<nev;j++){
	      cc[i*nev+j]=cc[j*nev+i];
	    }
	  }

	}else if(coriolis){
	      
	  NNEW(adc,double,neq[1]);
	  NNEW(auc,double,nzs[1]);
	  FORTRAN(mafillcorio,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
			       xboun,nboun,
			       ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
			       nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
			       adc,auc,nactdof,icol,jq,irow,neq,nzl,nmethod,
			       ikmpc,ilmpc,ikboun,ilboun,
			       elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
			       ielorien,norien,orab,ntmat_,
			       t0,t0,ithermal,prestr,iprestr,vold,iperturb,sti,
			       nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
			       xstiff,npmat_,&dtime,matname,mi,ncmat_,
			       ttime,&time,istep,&iinc,ibody,ielprop,prop));

	  /*  zc = damping matrix * eigenmodes */
	      
	  NNEW(zc,double,neq[1]*nev);
	  for(i=0;i<nev;i++){
	    FORTRAN(op,(&neq[1],&z[(long long)i*neq[1]],&zc[i*neq[1]],adc,auc,
			jq,irow));
	  }
	  SFREE(adc);SFREE(auc);
	      
	  /* cco is the reduced coriolis matrix (coriolis matrix mapped onto
	     space spanned by eigenmodes) */
		  
	  NNEW(cco,double,nev*nev);
	  for(i=0;i<nev;i++){
	    for(j=0;j<=i;j++){
	      for(k=0;k<neq[1];k++){
		cco[i*nev+j]+=z[(long long)j*neq[1]+k]*zc[i*neq[1]+k];
	      }
	    }
	  }
	  SFREE(zc);

	  /* multiplying with the driving frequency */

	  for(i=0;i<nev;i++){
	    for(j=0;j<=i;j++){
	      cco[i*nev+j]*=freq[l];
	    }
	  }

	  /* transpose is the negative of the matrix */

	  for(i=0;i<nev;i++){
	    for(j=i+1;j<nev;j++){
	      cco[i*nev+j]=-cco[j*nev+i];
	    }
	  }

	  /* adding cco to cc */

	  for(i=0;i<nev;i++){
	    for(j=0;j<nev;j++){
	      cc[i*nev+j]+=cco[j*nev+i];
	    }
	  }
	  SFREE(cco);
	}
	SFREE(xstiff);
      }
	  
      /* calculating the instantaneous loads (forces, surface loading, 
	 centrifugal and gravity loading or temperature) */
	  
      FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
			xloadold,xload,xloadact,iamload,nload,ibody,xbody,
			nbody,xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,namta,
			nam,ampli,&time,&reltime,ttime,&dtime,ithermal,nmethod,
			xbounold,xboun,xbounact,iamboun,nboun,nodeboun,ndirboun,
			nodeforc,ndirforc,istep,&iinc,co,vold,itg,&ntg,amname,
			ikboun,ilboun,nelemload,sideload,mi,
			ntrans,trab,inotr,veold,integerglob,doubleglob,tieset,istartset,
			iendset,ialset,&ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc,
			ipobody,iponoel,inoel,ipkon,kon,ielprop,prop,ielmat,
			shcon,nshcon,rhcon,nrhcon,cocon,ncocon,ntmat_,lakon,set,nset));
	  
      /* real part of forces */
	  
      for (i=0;i<*nforc;i++){
	xforcr[i]=xforcact[i]*(1-iphaseforc[i]);
      }
	  
      for (i=0;i<*nload;i++){
	for(j=0;j<2;j++){
	  xloadr[2*i+j]=xloadact[2*i+j]*(1-iphaseload[i]);
	}
      }

      for(i=0;i<*nbody;i++){
	for(j=0;j<7;j++){
	  xbodyr[7*i+j]=xbodyact[7*i+j];
	}
	if(ibody[3*i+2]==2){
	  xbodyr[7*i]=0.;
	}
      }
	  
      /* calculating the instantaneous loading vector */
	  
      rhsmain(co,nk,kon,ipkon,lakon,ne,
	      ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcr,
	      nforc,nelemload,sideload,xloadr,nload,xbodyr,
	      ipobody,nbody,cgr,br,nactdof,&neq[1],nmethod,
	      ikmpc,ilmpc,elcon,nelcon,rhcon,nrhcon,
	      alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
	      t0,t1act,ithermal,iprestr,vold,iperturb,iexpl,plicon,
	      nplicon,plkcon,nplkcon,
	      npmat_,ttime,&time,istep,&iinc,&dtime,physcon,ibody,
	      xbodyold,&reltime,veold,matname,mi,ikactmechr,
	      &nactmechr,ielprop,prop,sti,xstateini,xstate,nstate_,
	      ntrans,inotr,trab,fnext);
	  
      /* real modal coefficients */
	  
      if(!iprescribedboundary){

	if(nherm==1){
	  if(!cyclicsymmetry){
	    for(i=0;i<nev;i++){
	      i2=(long long)i*neq[1];
	      aa[i]=0.;
	      if(nactmechr<neq[1]/2){
		for(j=0;j<nactmechr;j++){
		  aa[i]+=z[i2+ikactmechr[j]]*br[ikactmechr[j]];
		}
	      }else{
		for(j=0;j<neq[1];j++){
		  aa[i]+=z[i2+j]*br[j];
		}
	      }
	    }
	  }else{
	    for(i=0;i<nev;i++){aa[i]=0.;}
	    for(j=0;j<nactmechr;j++){
	      for(i=0;i<nev;i++){
		FORTRAN(nident,(izdof,&ikactmechr[j],&nzdof,&id));
		if(id!=0){
		  if(izdof[id-1]==ikactmechr[j]){
		    aa[i]+=z[(long long)i*nzdof+id-1]*br[ikactmechr[j]];
		  }else{printf(" *ERROR in steadystate\n");FORTRAN(stop,());}
		}else{printf(" *ERROR in steadystate\n");FORTRAN(stop,());}
	      }
	    }
	  }
	}else{
	  if(!cyclicsymmetry){
	    for(i=0;i<nev;i++){
	      aa[2*i]=0.;aa[2*i+1]=0.;
	      if(nactmechr<neq[1]/2){
		for(j=0;j<nactmechr;j++){
		  aa[2*i]+=z[(long long)2*i*neq[1]+ikactmechr[j]]*br[ikactmechr[j]];
		  aa[2*i+1]+=z[(long long)(2*i+1)*neq[1]+ikactmechr[j]]*br[ikactmechr[j]];
		}
	      }else{
		for(j=0;j<neq[1];j++){
		  aa[2*i]+=z[(long long)2*i*neq[1]+j]*br[j];
		  aa[2*i+1]+=z[(long long)(2*i+1)*neq[1]+j]*br[j];
		}
	      }
	    }
	  }else{
	    for(i=0;i<nev;i++){aa[2*i]=0.;aa[2*i+1]=0.;}
	    for(j=0;j<nactmechr;j++){
	      for(i=0;i<nev;i++){
		FORTRAN(nident,(izdof,&ikactmechr[j],&nzdof,&id));
		if(id!=0){
		  if(izdof[id-1]==ikactmechr[j]){
		    aa[i]+=z[(long long)2*i*nzdof+id-1]*br[ikactmechr[j]];
		    aa[i]+=z[(long long)(2*i+1)*nzdof+id-1]*br[ikactmechr[j]];
		  }else{printf(" *ERROR in steadystate\n");FORTRAN(stop,());}
		}else{printf(" *ERROR in steadystate\n");FORTRAN(stop,());}
	      }
	    }
	  }
	}
	  
      }else{
	  
	/* correction for nonzero SPC's */
	      
	/* next statement makes sure that br is reset to zero at the
	   start of rhs.f */
	nactmechr=neq[1];

	/* real part of boundary conditions */
	      
	for (i=0;i<*nboun;i++){
	  xbounr[i]=xbounact[i]*(1-iphaseboun[i]);
	  if(typeboun[i]=='A') xbounr[i]/=(-freq[l]*freq[l]);
	}
	      
	for(j=0;j<neq[1];j++){fr[j]=0.;ubr[j]=0.;}
	for(i=0;i<*nboun;i++){
	  ic=neq[1]+i;
	  for(j=jq[ic]-1;j<jq[ic+1]-1;j++){
	    ir=irow[j]-1;
	    fr[ir]=fr[ir]-au[j]*xbounr[i];
	    ubr[ir]=fr[ir];
	  }
	}

	/* check for zero force (points to rigid body motion) */

	dd=0.;
	for(j=0;j<neq[1];j++){
	  dd=dd+fr[j]*fr[j];
	}
	if(dd!=0.){

	  /* solve equation system */
		
	  if(*isolver==0){
#ifdef SPOOLES
	    spooles_solve(ubr,&neq[1]);
#endif
	  }
	  else if(*isolver==4){
#ifdef SGI
	    sgi_solve(ubr,token);
#endif
	  }
	  else if(*isolver==5){
#ifdef TAUCS
	    tau_solve(ubr,&neq[1]);
#endif
	  }
	  else if(*isolver==7){
#ifdef PARDISO
	    pardiso_solve(ubr,&neq[1],&symmetryflag,&inputformat,&nrhs);
#endif
	  }
	  else if(*isolver==8){
#ifdef PASTIX
	    pastix_solve(ubr,&neq[1],&symmetryflag,&nrhs);
#endif
	  }
	}else{

	  /* rigid body mode: all displacements are equal */
		
	  for(j=0;j<neq[1];j++){
	    ubr[j]=xbounr[0];
	  }
	}
		
	FORTRAN(op,(&neq[1],ubr,mubr,adb,aub,jq,irow));
      }
	  
      /* imaginary part of forces */
	  
      for (i=0;i<*nforc;i++){
	xforci[i]=xforcact[i]*iphaseforc[i];
      }
	  
      for (i=0;i<*nload;i++){
	for(j=0;j<2;j++){
	  xloadi[2*i+j]=xloadact[2*i+j]*iphaseload[i];
	}
      }
	  
      for(i=0;i<*nbody;i++){
	for(j=0;j<7;j++){
	  xbodyi[7*i+j]=xbodyact[7*i+j];
	}
	if(ibody[3*i+2]==1){
	  xbodyi[7*i]=0.;
	}
      }
	  
      /* calculating the instantaneous loading vector */
	  
      rhsmain(co,nk,kon,ipkon,lakon,ne,
	      ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforci,
	      nforc,nelemload,sideload,xloadi,nload,xbodyi,
	      ipobody,nbody,cgr,bi,nactdof,&neq[1],nmethod,
	      ikmpc,ilmpc,elcon,nelcon,rhcon,nrhcon,
	      alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
	      t0,t1act,ithermal,iprestr,vold,iperturb,iexpl,plicon,
	      nplicon,plkcon,nplkcon,
	      npmat_,ttime,&time,istep,&iinc,&dtime,physcon,ibody,
	      xbodyold,&reltime,veold,matname,mi,ikactmechi,
	      &nactmechi,ielprop,prop,sti,xstateini,xstate,nstate_,
	      ntrans,inotr,trab,fnext);
	  
      /* imaginary modal coefficients */
	  
      if(!iprescribedboundary){

	if(nherm==1){
	  if(!cyclicsymmetry){
	    for(i=0;i<nev;i++){
	      i2=(long long)i*neq[1];
	      bb[i]=0.;
	      if(nactmechi<neq[1]/2){
		for(j=0;j<nactmechi;j++){
		  bb[i]+=z[i2+ikactmechi[j]]*bi[ikactmechi[j]];
		}
	      }else{
		for(j=0;j<neq[1];j++){
		  bb[i]+=z[i2+j]*bi[j];
		}
	      }
	    }
	  }else{
	    for(i=0;i<nev;i++){bb[i]=0.;}
	    for(j=0;j<nactmechi;j++){
	      for(i=0;i<nev;i++){
		FORTRAN(nident,(izdof,&ikactmechi[j],&nzdof,&id));
		if(id!=0){
		  if(izdof[id-1]==ikactmechi[j]){
		    bb[i]+=z[(long long)i*nzdof+id-1]*bi[ikactmechi[j]];
		  }else{printf(" *ERROR in steadystate\n");FORTRAN(stop,());}
		}else{printf(" *ERROR in steadystate\n");FORTRAN(stop,());}
	      }
	    }
	  }
	}else{
	  if(!cyclicsymmetry){
	    for(i=0;i<nev;i++){
	      bb[2*i]=0.;bb[2*i+1]=0.;
	      if(nactmechi<neq[1]/2){
		for(j=0;j<nactmechi;j++){
		  bb[2*i]+=z[(long long)2*i*neq[1]+ikactmechi[j]]*bi[ikactmechi[j]];
		  bb[2*i+1]+=z[(long long)(2*i+1)*neq[1]+ikactmechi[j]]*bi[ikactmechi[j]];
		}
	      }else{
		for(j=0;j<neq[1];j++){
		  bb[2*i]+=z[(long long)2*i*neq[1]+j]*bi[j];
		  bb[2*i+1]+=z[(long long)(2*i+1)*neq[1]+j]*bi[j];
		}
	      }
	    }
	  }else{
	    for(i=0;i<nev;i++){bb[2*i]=0.;bb[2*i+1]=0.;}
	    for(j=0;j<nactmechi;j++){
	      for(i=0;i<nev;i++){
		FORTRAN(nident,(izdof,&ikactmechi[j],&nzdof,&id));
		if(id!=0){
		  if(izdof[id-1]==ikactmechi[j]){
		    bb[2*i]+=z[(long long)2*i*nzdof+id-1]*bi[ikactmechi[j]];
		    bb[2*i+1]+=z[(long long)(2*i+1)*nzdof+id-1]*bi[ikactmechi[j]];
		  }else{printf(" *ERROR in steadystate\n");FORTRAN(stop,());}
		}else{printf(" *ERROR in steadystate\n");FORTRAN(stop,());}
	      }
	    }
	  }
	}

      }else{
	  
	/* correction for nonzero SPC's */
	      
	/* next statement makes sure that bi is reset to zero at the
	   start of rhs.f */

	nactmechi=neq[1];
	      
	/* imaginary part of boundary conditions */
	      
	for (i=0;i<*nboun;i++){
	  xbouni[i]=xbounact[i]*iphaseboun[i];
	  if(typeboun[i]=='A') xbouni[i]/=(-freq[l]*freq[l]);
	}
	      
	for(j=0;j<neq[1];j++){fi[j]=0.;ubi[j]=0.;}
	for(i=0;i<*nboun;i++){
	  ic=neq[1]+i;
	  for(j=jq[ic]-1;j<jq[ic+1]-1;j++){
	    ir=irow[j]-1;
	    fi[ir]=fi[ir]-au[j]*xbouni[i];
	    ubi[ir]=fi[ir];
	  }
	}

	/* check for zero force (points to rigid body motion) */

	dd=0.;
	for(j=0;j<neq[1];j++){
	  dd=dd+fi[j]*fi[j];
	}
	if(dd!=0.){

	  /* solve equation system */
		
	  if(*isolver==0){
#ifdef SPOOLES
	    spooles_solve(ubi,&neq[1]);
#endif
	  }
	  else if(*isolver==4){
#ifdef SGI
	    sgi_solve(ubi,token);
#endif
	  }
	  else if(*isolver==5){
#ifdef TAUCS
	    tau_solve(ubi,&neq[1]);
#endif
	  }
	  else if(*isolver==7){
#ifdef PARDISO
	    pardiso_solve(ubi,&neq[1],&symmetryflag,&inputformat,&nrhs);
#endif
	  }
	  else if(*isolver==8){
#ifdef PASTIX
	    pastix_solve(ubi,&neq[1],&symmetryflag,&nrhs);
#endif
	  }
	}else{

	  /* rigid body mode: all displacements are equal */
		
	  for(j=0;j<neq[1];j++){
	    ubi[j]=xbouni[0];
	  }
	}
	      
	FORTRAN(op,(&neq[1],ubi,mubi,adb,aub,jq,irow));
	  
	/* correction for prescribed boundary conditions */
	  
	/*	for(i=0;i<neq[1];i++){
	  br[i]+=freq[l]*(freq[l]*mubr[i]+alpham*mubi[i]+betam*fi[i]);
	  bi[i]+=freq[l]*(freq[l]*mubi[i]-alpham*mubr[i]-betam*fr[i]);
	  }*/
	  
	/* real and imaginary modal coefficients */

	/*	for(i=0;i<nev;i++){
	  aa[i]=0.;
	  for(j=0;j<neq[1];j++){
	    aa[i]+=z[(long long)i*neq[1]+j]*br[j];
	  }
	}
	      
	for(i=0;i<nev;i++){
	  bb[i]=0.;
	  for(j=0;j<neq[1];j++){
	    bb[i]+=z[(long long)i*neq[1]+j]*bi[j];
	  }
	  }*/

	/* calculating the eigenmodes x mass x particular solution */

	NNEW(umubr,double,nev);
	NNEW(umubi,double,nev);

	for(i=0;i<nev;i++){
	  for(j=0;j<neq[1];j++){
	    umubr[i]+=z[(long long)i*neq[1]+j]*mubr[j];
	  }
	}
	      
	for(i=0;i<nev;i++){
	  for(j=0;j<neq[1];j++){
	    umubi[i]+=z[(long long)i*neq[1]+j]*mubi[j];
	  }
	}
	  
	/* real and imaginary modal coefficients */

	for(i=0;i<nev;i++){
	  aa[i]=0.;
	  for(j=0;j<neq[1];j++){
	    aa[i]+=z[(long long)i*neq[1]+j]*br[j];
	  }
	  aa[i]+=freq[l]*(freq[l]*umubr[i]+fric[i]*umubi[i]);
	}
	      
	for(i=0;i<nev;i++){
	  bb[i]=0.;
	  for(j=0;j<neq[1];j++){
	    bb[i]+=z[(long long)i*neq[1]+j]*bi[j];
	  }
	  bb[i]+=freq[l]*(freq[l]*umubi[i]-fric[i]*umubr[i]);
	}
	SFREE(umubr);SFREE(umubi);
      }
	  
      /* calculating the modal coefficients */
	  
      if(nherm==1){
	if((dampmat==0)&&(coriolis==0)){
	  for(i=0;i<nev;i++){
	    dd=pow(d[i]-pow(freq[l],2),2)+
	      pow(fric[i],2)*pow(freq[l],2);
	    bjr[i]=(aa[i]*(d[i]-freq[l]*freq[l])+
		    bb[i]*fric[i]*freq[l])/dd;
	    bji[i]=(bb[i]*(d[i]-freq[l]*freq[l])-
		    aa[i]*fric[i]*freq[l])/dd;
	  }
	  /*    printf("old l=%" ITGFORMAT ",bjr=%f,bji=%f\n",l,bjr[0],bji[0]);*/
	}else{
	  nev2=2*nev;
	  NNEW(am,double,nev2*nev2);
	  NNEW(bm,double,nev2);
	  NNEW(ipiv,ITG,nev2);
		  
	  for(i=0;i<nev2;i++){
	    for(j=0;j<nev2;j++){
	      am[i*nev2+j]=0.;
	    }
	    bm[i]=0.;
	  }
	  for(i=0;i<nev;i++){
	    am[i*nev2+i]=d[i]-freq[l]*freq[l];
	    am[(i+nev)*nev2+i]=-fric[i]*freq[l];
	    bm[i]=aa[i];
	    am[i*nev2+nev+i]=-am[(i+nev)*nev2+i];
	    am[(i+nev)*nev2+nev+i]=am[i*nev2+i];
	    bm[nev+i]=bb[i];
	    for(j=0;j<nev;j++){
	      am[(j+nev)*nev2+i]=am[(j+nev)*nev2+i]
		-cc[i*nev+j];//*freq[l];
	      am[j*nev2+nev+i]=am[j*nev2+nev+i]
		+cc[i*nev+j];//*freq[l];
	    }
	  }
		  
	  /* solving the system of equations */
		  
	  FORTRAN(dgesv,(&nev2,&nrhs,am,&nev2,ipiv,bm,&nev2,&info));
	  if(info!=0){
	    printf(" *ERROR in steadystate: fatal termination of dgesv\n");
	    printf("       info=%" ITGFORMAT "\n",info);
	    FORTRAN(stop,());
	  }
		  
	  /* storing the solution in bjr and bji */
		  
	  for(i=0;i<nev;i++){
	    bjr[i]=bm[i];
	    bji[i]=bm[nev+i];
	  }

	  SFREE(am);SFREE(bm);SFREE(ipiv);
	}
      }else{
	nev2=2*nev;
	NNEW(am,double,nev2*nev2);
	NNEW(bm,double,nev2);
	NNEW(ipiv,ITG,nev2);
		  
	NNEW(ax,double,nev);
	NNEW(bx,double,nev);
	for(i=0;i<nev;i++){
	  ax[i]=-pow(freq[l],2)-freq[l]*fric[2*i+1]+d[2*i];
	  bx[i]=-freq[l]*fric[2*i]-d[2*i+1];
	}
	for(i=0;i<nev;i++){
	  for(j=0;j<nev;j++){
	    am[j*nev2+i]=xmr[j*nev+i]*ax[j]+xmi[j*nev+i]*bx[j];
	    am[(j+nev)*nev2+i]=xmr[j*nev+i]*bx[j]-xmi[j*nev+i]*bx[j];
	    am[j*nev2+nev+i]=xmi[j*nev+i]*ax[j]-xmr[j*nev+i]*bx[j];
	    am[(j+nev)*nev2+nev+i]=xmi[j*nev+i]*bx[j]+xmr[j*nev+i]*ax[j];
	  }
	  bm[i]=aa[i]-bb[nev+i];
	  bm[nev+i]=bb[i]+aa[nev+i];
	}
	SFREE(ax);SFREE(bx);
		  
	/* solving the system of equations */
	      
	FORTRAN(dgesv,(&nev2,&nrhs,am,&nev2,ipiv,bm,&nev2,&info));
	if(info!=0){
	  printf(" *ERROR in steadystate: fatal termination of dgesv\n");
	  printf("       info=%" ITGFORMAT "\n",info);
	  FORTRAN(stop,());
	}
	      
	/* storing the solution in bjr and bji */
	      
	for(i=0;i<nev;i++){
	  bjr[i]=bm[i];
	  bji[i]=bm[nev+i];
	}
	      
	SFREE(am);SFREE(bm);SFREE(ipiv);
      }

      /* storing the participation factors */

      if(nherm==1){
	NNEW(eig,double,nev);
	for(i=0;i<nev;i++){
	  eig[i]=sqrt(d[i]);
	}
      }else{
	      
	NNEW(eig,double,2*nev);
	for(i=0;i<nev;i++){
	  eig[2*i]=sqrt(sqrt(d[2*i]*d[2*i]+d[2*i+1]*d[2*i+1])+d[2*i])
	    /sqrt(2.);
	  eig[2*i+1]=sqrt(sqrt(d[2*i]*d[2*i]+d[2*i+1]*d[2*i+1])-d[2*i])
	    /sqrt(2.);
	  if(d[2*i+1]<0.) eig[2*i+1]=-eig[2*i+1];
	}
      }
	  
      mode=0;
      FORTRAN(writepf,(eig,bjr,bji,&time,&nev,&mode,&nherm));
      SFREE(eig);
      
      /* calculating the real response */
      
      if(iprescribedboundary){
	if(nmdnode==0){
	  memcpy(&br[0],&ubr[0],sizeof(double)*neq[1]);
	}else{
	  for(i=0;i<nmddof;i++){
	    br[imddof[i]]=ubr[imddof[i]];
	  }
	}
      }
      else{
	if(nmdnode==0){
	  DMEMSET(br,0,neq[1],0.);
	}else{
	  for(i=0;i<nmddof;i++){
	    br[imddof[i]]=0.;
	  }
	}
      }
	  
      if(!cyclicsymmetry){
	if(nmdnode==0){
	  for(i=0;i<neq[1];i++){
	    for(j=0;j<nev;j++){
	      br[i]+=bjr[j]*z[(long long)j*neq[1]+i];
	    }
	  }
	}else{
	  for(i=0;i<nmddof;i++){
	    for(j=0;j<nev;j++){
	      br[imddof[i]]+=bjr[j]*z[(long long)j*neq[1]+imddof[i]];
	    }
	  }
	}
      }else{
	for(i=0;i<nmddof;i++){
	  FORTRAN(nident,(izdof,&imddof[i],&nzdof,&id));
	  if(id!=0){
	    if(izdof[id-1]==imddof[i]){
	      for(j=0;j<nev;j++){
		br[imddof[i]]+=bjr[j]*z[(long long)j*nzdof+id-1];
	      }
	    }else{printf(" *ERROR in steadystate\n");FORTRAN(stop,());}
	  }else{printf(" *ERROR in steadystate\n");FORTRAN(stop,());}
	}
      }
	  
      if(nmdnode==0){
	DMEMSET(vr,0,mt**nk,0.);
      }else{
	for(jj=0;jj<nmdnode;jj++){
	  i=imdnode[jj]-1;
	  for(j=1;j<4;j++){
	    vr[mt*i+j]=0.;
	  }
	}
      }
	  
      /* calculating the imaginary response */
      
      if(iprescribedboundary){
	if(nmdnode==0){
	  memcpy(&bi[0],&ubi[0],sizeof(double)*neq[1]);
	}else{
	  for(i=0;i<nmddof;i++){
	    bi[imddof[i]]=ubi[imddof[i]];
	  }
	}
      }
      else{
	if(nmdnode==0){
	  DMEMSET(bi,0,neq[1],0.);
	}else{
	  for(i=0;i<nmddof;i++){
	    bi[imddof[i]]=0.;
	  }
	}
      }
	  
      if(!cyclicsymmetry){
	if(nmdnode==0){
	  for(i=0;i<neq[1];i++){
	    for(j=0;j<nev;j++){
	      bi[i]+=bji[j]*z[(long long)j*neq[1]+i];
	    }
	  }
	}else{
	  for(i=0;i<nmddof;i++){
	    for(j=0;j<nev;j++){
	      bi[imddof[i]]+=bji[j]*z[(long long)j*neq[1]+imddof[i]];
	    }
	  }
	}
      }else{
	for(i=0;i<nmddof;i++){
	  FORTRAN(nident,(izdof,&imddof[i],&nzdof,&id));
	  if(id!=0){
	    if(izdof[id-1]==imddof[i]){
	      for(j=0;j<nev;j++){
		bi[imddof[i]]+=bji[j]*z[(long long)j*nzdof+id-1];
	      }
	    }else{printf(" *ERROR in steadystate\n");FORTRAN(stop,());}
	  }else{printf(" *ERROR in steadystate\n");FORTRAN(stop,());}
	}
      }
	  
	  
      if(nmdnode==0){
	DMEMSET(vi,0,mt**nk,0.);
      }else{
	for(jj=0;jj<nmdnode;jj++){
	  i=imdnode[jj]-1;
	  for(j=1;j<4;j++){
	    vi[mt*i+j]=0.;
	  }
	}
      }

      /* real response */

      if(iprescribedboundary){
      
	/* calculating displacements/temperatures */
      
	FORTRAN(dynresults,(nk,vr,ithermal,nactdof,vold,nodeboun,
			    ndirboun,xbounr,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,
			    br,bi,veold,&dtime,mi,imdnode,&nmdnode,imdboun,&nmdboun,
			    imdmpc,&nmdmpc,nmethod,&timem));

	results(co,nk,kon,ipkon,lakon,ne,vr,stnr,inum,
		stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		ielmat,ielorien,norien,orab,ntmat_,t0,t1,
		ithermal,prestr,iprestr,filab,eme,emn,een,
		iperturb,f,fn,nactdof,&iout,qa,
		vold,br,nodeboun,ndirboun,xbounr,nboun,
		ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],
		veold,accold,&bet,&gam,&dtime,&time,&xnull,
		plicon,nplicon,plkcon,nplkcon,
		xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
		&icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,
		enern,emeini,xstaten,eei,enerini,cocon,ncocon,
		set,nset,istartset,iendset,ialset,nprint,prlab,prset,
		qfx,qfn,trab,inotr,ntrans,fmpc,nelemload,nload,
		ikmpc,ilmpc,istep,&iinc,springarea,&reltime,&ne0,
		thicke,shcon,nshcon,
		sideload,xload,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		&mortar,islavact,cdn,islavnode,nslavnode,&ntie,clearini,
		islavsurf,ielprop,prop,energyini,energy,&iit,iponoel,
		inoel,nener,orname,&network,ipobody,xbodyact,ibody,
		typeboun,itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
		islavelinv,autloc,irowtloc,jqtloc,&nboun2,
		ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
		labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
		&intscheme,physcon);}
      else{
      
	/* calculating displacements/temperatures */
      
	FORTRAN(dynresults,(nk,vr,ithermal,nactdof,vold,nodeboun,
			    ndirboun,xbounact,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,
			    br,bi,veold,&dtime,mi,imdnode,&nmdnode,imdboun,&nmdboun,
			    imdmpc,&nmdmpc,nmethod,&timem));

	results(co,nk,kon,ipkon,lakon,ne,vr,stnr,inum,
		stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		ielmat,ielorien,norien,orab,ntmat_,t0,t1,
		ithermal,prestr,iprestr,filab,eme,emn,een,
		iperturb,f,fn,nactdof,&iout,qa,
		vold,br,nodeboun,ndirboun,xbounact,nboun,
		ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],
		veold,accold,&bet,&gam,&dtime,&time,&xnull,
		plicon,nplicon,plkcon,nplkcon,
		xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
		&icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,
		enern,emeini,xstaten,eei,enerini,cocon,ncocon,
		set,nset,istartset,iendset,ialset,nprint,prlab,prset,
		qfx,qfn,trab,inotr,ntrans,fmpc,nelemload,nload,
		ikmpc,ilmpc,istep,&iinc,springarea,&reltime,&ne0,
		thicke,shcon,nshcon,
		sideload,xload,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		&mortar,islavact,cdn,islavnode,nslavnode,&ntie,
		clearini,islavsurf,ielprop,prop,energyini,energy,&iit,
		iponoel,inoel,nener,orname,&network,ipobody,xbodyact,
		ibody,typeboun,itiefac,tieset,smscale,&mscalmethod,nbody,
		t0g,t1g,islavelinv,autloc,irowtloc,jqtloc,&nboun2,
		ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
		labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
		&intscheme,physcon);
	      
	if(nmdnode==0){
	  DMEMSET(br,0,neq[1],0.);
	}else{
	  for(i=0;i<nmddof;i++){
	    br[imddof[i]]=0.;
	  }
	}
      }
	  
      (*kode)++;
	  
      mode=-1;
      if(strcmp1(&filab[1044],"ZZS")==0){
	NNEW(neigh,ITG,40**ne);
	NNEW(ipneigh,ITG,*nk);
      }

      frd(co,&nkg,kon,ipkon,lakon,&neg,vr,stnr,inum,nmethod,
	  kode,filab,een,t1,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	  nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	  mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,&neg,
	  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	  thicke,jobnamec,output,qfx,cdn,&mortar,cdnr,cdni,nmat,
	  ielprop,prop,sti);

      if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);}

      /* imaginary response */

      if(iprescribedboundary){
      
	/* calculating displacements/temperatures */
      
	FORTRAN(dynresults,(nk,vi,ithermal,nactdof,vold,nodeboun,
			    ndirboun,xbouni,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,
			    bi,br,veold,&dtime,mi,imdnode,&nmdnode,imdboun,&nmdboun,
			    imdmpc,&nmdmpc,nmethod,&time));

	results(co,nk,kon,ipkon,lakon,ne,vi,stni,inum,
		stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		ielmat,ielorien,norien,orab,ntmat_,t0,t1,
		ithermal,prestr,iprestr,filab,eme,emn,een,
		iperturb,f,fn,nactdof,&iout,qa,
		vold,bi,nodeboun,ndirboun,xbouni,nboun,
		ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],
		veold,accold,&bet,&gam,&dtime,&time,&xnull,
		plicon,nplicon,plkcon,nplkcon,
		xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
		&icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,
		enern,emeini,xstaten,eei,enerini,cocon,ncocon,
		set,nset,istartset,iendset,ialset,nprint,prlab,prset,
		qfx,qfn,trab,inotr,ntrans,fmpc,nelemload,nload,
		ikmpc,ilmpc,istep,&iinc,springarea,&reltime,&ne0,
		thicke,shcon,nshcon,
		sideload,xload,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		&mortar,islavact,cdn,islavnode,nslavnode,&ntie,clearini,
		islavsurf,ielprop,prop,energyini,energy,&iit,iponoel,
		inoel,nener,orname,&network,ipobody,xbodyact,ibody,
		typeboun,itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
		islavelinv,autloc,irowtloc,jqtloc,&nboun2,
		ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
		labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
		&intscheme,physcon);}
      else{ 
      
	/* calculating displacements/temperatures */
      
	FORTRAN(dynresults,(nk,vi,ithermal,nactdof,vold,nodeboun,
			    ndirboun,xbounact,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,
			    bi,br,veold,&dtime,mi,imdnode,&nmdnode,imdboun,&nmdboun,
			    imdmpc,&nmdmpc,nmethod,&time));

	results(co,nk,kon,ipkon,lakon,ne,vi,stni,inum,
		stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		ielmat,ielorien,norien,orab,ntmat_,t0,t1,
		ithermal,prestr,iprestr,filab,eme,emn,een,
		iperturb,f,fn,nactdof,&iout,qa,
		vold,bi,nodeboun,ndirboun,xbounact,nboun,
		ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],
		veold,accold,&bet,&gam,&dtime,&time,&xnull,
		plicon,nplicon,plkcon,nplkcon,
		xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
		&icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,
		enern,emeini,xstaten,eei,enerini,cocon,ncocon,
		set,nset,istartset,iendset,ialset,nprint,prlab,prset,
		qfx,qfn,trab,inotr,ntrans,fmpc,nelemload,nload,
		ikmpc,ilmpc,istep,&iinc,springarea,&reltime,&ne0,
		thicke,shcon,nshcon,
		sideload,xload,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		&mortar,islavact,cdn,islavnode,nslavnode,&ntie,
		clearini,islavsurf,ielprop,prop,energyini,energy,&iit,
		iponoel,inoel,nener,orname,&network,ipobody,xbodyact,
		ibody,typeboun,itiefac,tieset,smscale,&mscalmethod,nbody,
		t0g,t1g,islavelinv,autloc,irowtloc,jqtloc,&nboun2,
		ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
		labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
		&intscheme,physcon);

	if(nmdnode==0){
	  DMEMSET(bi,0,neq[1],0.);
	}else{
	  for(i=0;i<nmddof;i++){
	    bi[imddof[i]]=0.;
	  }
	}
      }
	  
      /* calculating the magnitude and phase */
	  
      if(strcmp1(&filab[870],"PU")==0){
	      
	constant=180./pi;
	NNEW(va,double,mt**nk);
	NNEW(vp,double,mt**nk);
	      
	if(*ithermal<=1){
	  if(nmdnode==0){
	    for(i=0;i<*nk;i++){
	      for(j=1;j<4;j++){
		vreal=vr[mt*i+j];
		va[mt*i+j]=sqrt(vr[mt*i+j]*vr[mt*i+j]+vi[mt*i+j]*vi[mt*i+j]);
		if(va[mt*i+j]<1.e-10){
		  vp[mt*i+j]=0.;
		}else if(fabs(vreal)<1.e-10){
		  if(vi[mt*i+j]>0.){vp[mt*i+j]=90.;}
		  else{vp[mt*i+j]=-90.;}
		}
		else{
		  vp[mt*i+j]=atan(vi[mt*i+j]/vreal)*constant;
		  if(vreal<0.) vp[mt*i+j]+=180.;
		}
	      }
	    }
	  }else{
	    for(jj=0;jj<nmdnode;jj++){
	      i=imdnode[jj]-1;
	      for(j=1;j<4;j++){
		vreal=vr[mt*i+j];
		va[mt*i+j]=sqrt(vr[mt*i+j]*vr[mt*i+j]+vi[mt*i+j]*vi[mt*i+j]);
		if(va[mt*i+j]<1.e-10){
		  vp[mt*i+j]=0.;
		}else if(fabs(vreal)<1.e-10){
		  if(vi[mt*i+j]>0.){vp[mt*i+j]=90.;}
		  else{vp[mt*i+j]=-90.;}
		}
		else{
		  vp[mt*i+j]=atan(vi[mt*i+j]/vreal)*constant;
		  if(vreal<0.) vp[mt*i+j]+=180.;
		}
	      }
	    }
	  }
	}
	else{
	  if(nmdnode==0){
	    for(i=0;i<*nk;i++){
	      vreal=vr[mt*i];
	      va[mt*i]=sqrt(vr[mt*i]*vr[mt*i]+vi[mt*i]*vi[mt*i]);
		if(va[mt*i+j]<1.e-10){
		  vp[mt*i+j]=0.;
		}else if(fabs(vreal)<1.e-10){
		if(vi[mt*i]>0){vp[mt*i]=90.;}
		else{vp[mt*i]=-90.;}
	      }
	      else{
		vp[mt*i]=atan(vi[mt*i]/vreal)*constant;
		if(vreal<0.) vp[mt*i]+=180.;
	      }
	    }
	  }else{
	    for(jj=0;jj<nmdnode;jj++){
	      i=imdnode[jj]-1;
	      vreal=vr[mt*i];
	      va[mt*i]=sqrt(vr[mt*i]*vr[mt*i]+vi[mt*i]*vi[mt*i]);
		if(va[mt*i+j]<1.e-10){
		  vp[mt*i+j]=0.;
		}else if(fabs(vreal)<1.e-10){
		if(vi[mt*i]>0){vp[mt*i]=90.;}
		else{vp[mt*i]=-90.;}
	      }
	      else{
		vp[mt*i]=atan(vi[mt*i]/vreal)*constant;
		if(vreal<0.) vp[mt*i]+=180.;
	      }
	    }
	  }
	}
      }
	  
      if(strcmp1(&filab[1479],"PHS")==0){
	      
	constant=180./pi;
	NNEW(stna,double,6**nk);
	NNEW(stnp,double,6**nk);
	      
	if(*ithermal<=1){
	  if(nmdnode==0){
	    for(i=0;i<*nk;i++){
	      for(j=0;j<6;j++){
		vreal=stnr[6*i+j];
		stna[6*i+j]=sqrt(stnr[6*i+j]*stnr[6*i+j]+stni[6*i+j]*stni[6*i+j]);
		if(fabs(vreal)<1.e-10){
		  if(stni[6*i+j]>0.){stnp[6*i+j]=90.;}
		  else{stnp[6*i+j]=-90.;}
		}
		else{
		  stnp[6*i+j]=atan(stni[6*i+j]/vreal)*constant;
		  if(vreal<0.) stnp[6*i+j]+=180.;
		}
	      }
	    }
	  }else{
	    for(jj=0;jj<nmdnode;jj++){
	      i=imdnode[jj]-1;
	      for(j=0;j<6;j++){
		vreal=stnr[6*i+j];
		stna[6*i+j]=sqrt(stnr[6*i+j]*stnr[6*i+j]+stni[6*i+j]*stni[6*i+j]);
		if(fabs(vreal)<1.e-10){
		  if(stni[6*i+j]>0.){stnp[6*i+j]=90.;}
		  else{stnp[6*i+j]=-90.;}
		}
		else{
		  stnp[6*i+j]=atan(stni[6*i+j]/vreal)*constant;
		  if(vreal<0.) stnp[6*i+j]+=180.;
		}
	      }
	    }
	  }
	}
      }
 
      //	  (*kode)++;
      mode=0;
	  
      if(strcmp1(&filab[1044],"ZZS")==0){
	NNEW(neigh,ITG,40**ne);
	NNEW(ipneigh,ITG,*nk);
      }

      frd(co,&nkg,kon,ipkon,lakon,&neg,vi,stni,inum,nmethod,
	  kode,filab,een,t1,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	  nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	  mi,stx,va,vp,stna,stnp,vmax,stnmax,&ngraph,veold,ener,&neg,
	  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	  thicke,jobnamec,output,qfx,cdn,&mortar,cdnr,cdni,nmat,
	  ielprop,prop,sti);

      if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);}

      SFREE(va);SFREE(vp);SFREE(stna);SFREE(stnp);
	  
    }
      
    /* restoring the imaginary loading */
      
    SFREE(iphaseforc);SFREE(xforcr);SFREE(xforci);

    SFREE(iphaseload);SFREE(xloadr);SFREE(xloadi);
      
    SFREE(xbodyr);SFREE(xbodyi);
      
    if(iprescribedboundary){
      for (i=0;i<*nboun;i++){
	if(iphaseboun[i]==1){
	  nodeboun[i]=nodeboun[i]+*nk;
	}
      }
      SFREE(iphaseboun);
    }
      
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
      
    SFREE(br);SFREE(bi);SFREE(bjr);SFREE(bji),SFREE(freq);
    SFREE(xforcact);SFREE(xloadact);SFREE(xbounact);SFREE(aa);SFREE(bb);
    SFREE(ampli);SFREE(xbodyact);SFREE(vr);SFREE(vi);

    //if(*nbody>0) SFREE(ipobody);
      
    if(*ithermal==1) SFREE(t1act);
      
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
      SFREE(xbounr);SFREE(xbouni);SFREE(fr);SFREE(fi);SFREE(ubr);SFREE(ubi);
      SFREE(mubr);SFREE(mubi);
    }

    SFREE(ikactmechr);SFREE(ikactmechi);
      
    if(intpointvar==1){
      SFREE(fn);
      SFREE(stnr);SFREE(stni);SFREE(stx);SFREE(eei);

      if(*ithermal>1) {SFREE(qfn);SFREE(qfx);}
	  
      if(strcmp1(&filab[261],"E   ")==0) SFREE(een);
      if(strcmp1(&filab[522],"ENER")==0) SFREE(enern);
      if(strcmp1(&filab[2697],"ME  ")==0) SFREE(emn);
	  
      if(*nener==1){SFREE(stiini);SFREE(emeini);SFREE(enerini);}
    }

  }else{
 
    /**********************************************************************

      NONHARMONIC PERIODIC EXCITATION

    **********************************************************************/

    NNEW(ikactmech,ITG,neq[1]);
    nactmech=0;
      
    NNEW(xforcact,double,nfour**nforc);
    NNEW(xloadact,double,nfour*2**nload);
    NNEW(xbodyact,double,nfour*7**nbody);
    NNEW(xbounact,double,nfour**nboun);
    NNEW(xbounacttime,double,nfour**nboun);
    if(*ithermal==1) NNEW(t1act,double,*nk);

    NNEW(r,double,nfour);
    NNEW(wsave,double,2*nfour);
    NNEW(isave,ITG,15);
      
    /*  check for nonzero SPC's */
      
    iprescribedboundary=0;
    for(i=0;i<*nboun;i++){
      if(fabs(xboun[i])>1.e-10){
	iprescribedboundary=1;
	break;
      }
    }

    if((iprescribedboundary)&&(cyclicsymmetry)){
      printf(" *ERROR in steadystate: prescribed boundaries are not allowed in combination with cyclic symmetry\n");
      FORTRAN(stop,());
    }

    /* calculating the damping coefficients = friction coefficient*2*eigenvalue */
      
    if(xmodal[10]<0){
      for(i=0;i<nev;i++){
	if(sqrt(d[i])>(1.e-10)){
	  fric[i]=(alpham+betam*d[i]);
	}
	else {
	  printf(" *WARNING in steadystate: one of the frequencies is zero\n");
	  printf("         no Rayleigh mass damping allowed\n");
	  fric[i]=0.;
	}
      }
    }
    else{
      if(iprescribedboundary){
	printf(" *ERROR in steadystate: prescribed boundaries are not allowed in combination with direct modal damping\n");
	FORTRAN(stop,());
      }
	      
      /*copy the damping coefficients for every eigenfrequencie from xmodal[11....] */
      if(nev<(ITG)xmodal[10]){
	imax=nev;
	printf(" *WARNING in steadystate: too many modal damping coefficients applied\n");
	printf("         damping coefficients corresponding to nonexisting eigenvalues are ignored\n");
      }
      else{
	imax=(ITG)xmodal[10];
      }
      for(i=0; i<imax; i++){
	fric[i]=2.*sqrt(d[i])*xmodal[11+i];     
      }
	  
    }

    /* determining the load time history */
      
    NNEW(ampli,double,*nam); /* instantaneous amplitude */

    for(l=0;l<nfour;l++){

      time=tmin+(tmax-tmin)*(double)l/(double)nfour;
	      
      FORTRAN(tempload,(xforcold,xforc,&xforcact[l**nforc],iamforc,nforc,
			xloadold,xload,&xloadact[l*2**nload],iamload,nload,ibody,xbody,
			nbody,xbodyold,&xbodyact[l*7**nbody],t1old,t1,t1act,
			iamt1,nk,amta,namta,nam,ampli,&time,&reltime,ttime,&dtime,
			ithermal,nmethod,xbounold,xboun,&xbounact[l**nboun],iamboun,nboun,
			nodeboun,ndirboun,nodeforc,ndirforc,istep,&iinc,co,vold,itg,&ntg,
			amname,ikboun,ilboun,nelemload,sideload,mi,
			ntrans,trab,inotr,veold,integerglob,doubleglob,tieset,istartset,
			iendset,ialset,&ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc,
			ipobody,iponoel,inoel,ipkon,kon,ielprop,prop,ielmat,
			shcon,nshcon,rhcon,nrhcon,cocon,ncocon,ntmat_,lakon,set,nset));
	  
    }

    SFREE(ampli);

    for(i=0;i<l**nboun;i++){xbounacttime[i]=xbounact[i];}

    /* determining the load frequency history:
       frequency transform of the load time history */

    for(i=0;i<*nforc;i++){
      for(l=0;l<nfour;l++){
	r[l]=xforcact[l**nforc+i];
      }
      FORTRAN(drffti,(&nfour,wsave,isave));
      FORTRAN(drfftf,(&nfour,r,wsave,isave));
      for(l=0;l<nfour;l++){
	xforcact[l**nforc+i]=r[l]/nfour*2.;
      }
      xforcact[i]=xforcact[i]/2.;
    }

    for(i=0;i<*nload;i++){
      for(l=0;l<nfour;l++){
	r[l]=xloadact[l*2**nload+2*i];
      }
      FORTRAN(drffti,(&nfour,wsave,isave));
      FORTRAN(drfftf,(&nfour,r,wsave,isave));
      for(l=0;l<nfour;l++){
	xloadact[l*2**nload+2*i]=r[l]/nfour*2.;
      }
      xloadact[2*i]=xloadact[2*i]/2.;
    }

    for(i=0;i<*nbody;i++){
      for(l=0;l<nfour;l++){
	r[l]=xbodyact[l*7**nbody+7*i];
      }
      FORTRAN(drffti,(&nfour,wsave,isave));
      FORTRAN(drfftf,(&nfour,r,wsave,isave));
      for(l=0;l<nfour;l++){
	xbodyact[l*7**nbody+7*i]=r[l]/nfour*2.;
      }
      xbodyact[7*i]=xbodyact[7*i]/2.;
    }

    if(iprescribedboundary){
      for(i=0;i<*nboun;i++){
	for(l=0;l<nfour;l++){
	  r[l]=xbounact[l**nboun+i];
	}
	FORTRAN(drffti,(&nfour,wsave,isave));
	FORTRAN(drfftf,(&nfour,r,wsave,isave));
	for(l=0;l<nfour;l++){
	  xbounact[l**nboun+i]=r[l]/nfour*2.;
	}
	xbounact[i]=xbounact[i]/2.;
      }
	  
      /* LU decomposition of the stiffness matrix */
	  
      if(*isolver==0){
#ifdef SPOOLES
	spooles_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
		       &symmetryflag,&inputformat,&nzs[2]);
#else
	printf(" *ERROR in steadystate: the SPOOLES library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==4){
#ifdef SGI
	token=1;
	sgi_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],token);
#else
	printf(" *ERROR in steadystate: the SGI library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
	tau_factor(ad,&au,adb,aub,&sigma,icol,&irow,&neq[1],&nzs[1]);
#else
	printf(" *ERROR in steadystate: the TAUCS library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
		       &symmetryflag,&inputformat,jq,&nzs[2]);
#else
	printf(" *ERROR in steadystate: the PARDISO library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==8){
#ifdef PASTIX
	pastix_factor_main(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
		      &symmetryflag,&inputformat,jq,&nzs[2]);
#else
	printf(" *ERROR in steadystate: the PASTIX library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }


    }

    SFREE(r);SFREE(wsave);SFREE(isave);

    /* determining the frequency data points */
      
    NNEW(freqnh,double,ndata*(nev+1));
      
    ndatatot=0.;
    freqnh[0]=fmin;
    if(fabs(fmax-fmin)<1.e-10){
      ndatatot=1;
    }else{

      /* copy the eigenvalues and sort them in ascending order
         (important for values from distinct nodal diameters */

      NNEW(e,double,nev);
      for(i=0;i<nev;i++){e[i]=sqrt(d[i]);}
      FORTRAN(dsort,(e,&idummy,&nev,&iflag));

      for(i=0;i<nev;i++){
	if(i!=0){
	  if(fabs(e[i]-e[i-1])<1.e-5){continue;}
	}
	if(e[i]>=fmin){
	  if(e[i]<=fmax){
	    for(j=1;j<ndata;j++){
	      y=-1.+2.*j/((double)(ndata-1));
	      if(fabs(y)<1.e-10){freqnh[ndatatot+j]=
		  (freqnh[ndatatot]+e[i])/2.;}
	      else{
		freqnh[ndatatot+j]=(freqnh[ndatatot]+e[i])/2.+
		  (e[i]-freqnh[ndatatot])*pow(fabs(y),1./bias)
		  *y/(2.*fabs(y));
	      }
	    }
	    ndatatot+=ndata-1;
	  }
	  else{break;}
	}
      }
      SFREE(e);
      for(j=1;j<ndata;j++){
	y=-1.+2.*j/((double)(ndata-1));
	if(fabs(y)<1.e-10){freqnh[ndatatot+j]=
	    (freqnh[ndatatot]+fmax)/2.;}
	else{
	  freqnh[ndatatot+j]=(freqnh[ndatatot]+fmax)/2.+
	    (fmax-freqnh[ndatatot])*pow(fabs(y),1./bias)*
	    y/(2.*fabs(y));
	}
      }
      ndatatot+=ndata;
    }
    RENEW(freqnh,double,ndatatot);

    for(ii=0;ii<ndatatot;ii++){
      for(i=0;i<6*mi[0]**ne;i++){eme[i]=0.;}

      sprintf(description,"%12f",freqnh[ii]/(2.*pi));
	  
      NNEW(xforcr,double,*nforc);
      NNEW(xloadr,double,2**nload);
      NNEW(xbodyr,double,7**nbody);
      for(k=0;k<7**nbody;k++){xbodyr[k]=xbody[k];}
      if(iprescribedboundary){
	NNEW(xbounr,double,*nboun);
	NNEW(fr,double,neq[1]); /* force corresponding to real particular solution */
	NNEW(ubr,double,neq[1]); /* real particular solution */
	NNEW(mubr,double,neq[1]); /* mass times real particular solution */
      }
	  
      /* assigning the body forces to the elements */ 

      /*     if(*nbody>0){
	ifreebody=*ne+1;
	NNEW(ipobody,ITG,2*ifreebody**nbody);
	for(k=1;k<=*nbody;k++){
	  FORTRAN(bodyforce,(cbody,ibody,ipobody,nbody,set,istartset,
			     iendset,ialset,&inewton,nset,&ifreebody,&k));
	  RENEW(ipobody,ITG,2*(*ne+ifreebody));
	}
	RENEW(ipobody,ITG,2*(ifreebody-1));
	}*/
	  
      NNEW(br,double,neq[1]); /* load rhs vector (real part) */
      NNEW(bi,double,neq[1]); /* load rhs vector (imaginary part) */
      NNEW(btot,double,nfour*neq[1]);
      NNEW(bp,double,nfour*neq[1]);
	  
      NNEW(bjr,double,nev); /* real response modal decomposition */
      NNEW(bji,double,nev); /* imaginary response modal decomposition */
	  
      NNEW(aa,double,nev); /* modal coefficients of the real loading */
      NNEW(bb,double,nev); /* modal coefficients of the imaginary loading */
	  
      /* loop over all Fourier frequencies */
	  
      NNEW(freq,double,nfour);
	  
      for(l=0;l<nfour;l++){
	      
	/* frequency */
	      
	freq[l]=freqnh[ii]*floor((l+1.)/2.+0.1);

	/* calculating cc */

	if((dampmat||coriolis)){

	  /* determining the elastic constants in xstiff (needed in mafilldmss) */

	  iout=-1;
	  NNEW(xstiff,double,(long long)27*mi[0]**ne);
	  NNEW(v,double,mt**nk);
	  NNEW(f,double,*neq);
	  if(intpointvar!=1){
	    NNEW(fn,double,mt**nk);
	    NNEW(stx,double,6*mi[0]**ne);
	  }
	  results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
		  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		  ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
		  prestr,iprestr,filab,eme,emn,een,iperturb,
		  f,fn,nactdof,&iout,qa,vold,b,nodeboun,
		  ndirboun,xbounact,nboun,ipompc,
		  nodempc,coefmpc,labmpc,nmpc,&nmethodact,cam,neq,veold,accold,
		  &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		  xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
		  &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
		  emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
		  iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
		  fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
		  &reltime,&ne0,thicke,shcon,nshcon,
		  sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		  &mortar,islavact,cdn,islavnode,nslavnode,&ntie,clearini,
		  islavsurf,ielprop,prop,energyini,energy,&kscale,iponoel,
		  inoel,nener,orname,&network,ipobody,xbodyact,ibody,typeboun,
		  itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
		  islavelinv,autloc,irowtloc,jqtloc,&nboun2,
		  ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
		  labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
		  &intscheme,physcon);
	  SFREE(v);SFREE(f);
	  if(intpointvar!=1){SFREE(fn);SFREE(stx);}
	  iout=2;

	  for(i=0;i<nev*nev;i++){cc[i]=0.;}
		  
	  if(dampmat){

	    NNEW(adc,double,neq[1]);
	    NNEW(auc,double,nzs[1]);
	    mafilldmssmain(co,nk,kon,ipkon,lakon,ne,
			   ipompc,nodempc,coefmpc,nmpc,
			   nelemload,sideload,xload,nload,xbody,ipobody,
			   nbody,cgr,adc,auc,nactdof,jq,irow,neq,
			   nmethod,ikmpc,ilmpc,
			   elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
			   ielmat,ielorien,norien,orab,ntmat_,
			   t0,t0,ithermal,vold,iperturb,sti,
			   nzs,stx,iexpl,plicon,nplicon,plkcon,nplkcon,
			   xstiff,npmat_,&dtime,matname,mi,
			   ncmat_,physcon,ttime,&time,istep,&iinc,
			   ibody,xloadold,&reltime,veold,springarea,nstate_,
			   xstateini,xstate,thicke,integerglob,doubleglob,
			   tieset,istartset,iendset,ialset,&ntie,&nasym,pslavsurf,
			   pmastsurf,&mortar,clearini,ielprop,prop,&ne0,
			   &freq[l],ndamp,dacon,set,nset);
		      
	    SFREE(xstiff);
		      
	    /*  zc = damping matrix * eigenmodes */
		      
	    NNEW(zc,double,neq[1]*nev);
	    for(i=0;i<nev;i++){
	      FORTRAN(op,(&neq[1],&z[(long long)i*neq[1]],&zc[i*neq[1]],
			  adc,auc,jq,irow));
	    }
	    SFREE(adc);SFREE(auc);
		      
	    /* cc is the reduced damping matrix (damping matrix mapped onto
	       space spanned by eigenmodes) */
		      
	    for(i=0;i<nev*nev;i++){cc[i]=0.;}
	    for(i=0;i<nev;i++){
	      for(j=0;j<=i;j++){
		for(k=0;k<neq[1];k++){
		  cc[i*nev+j]+=z[(long long)j*neq[1]+k]*zc[i*neq[1]+k];
		}
	      }
	    }
	    SFREE(zc);
		      
	    /* symmetric part of cc matrix */
		      
	    for(i=0;i<nev;i++){
	      for(j=i;j<nev;j++){
		cc[i*nev+j]=cc[j*nev+i];
	      }
	    }
		      
	  }else if(coriolis){
		      
	    NNEW(adc,double,neq[1]);
	    NNEW(auc,double,nzs[1]);
	    FORTRAN(mafillcorio,(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
				 xboun,nboun,
				 ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforc,
				 nforc,nelemload,sideload,xload,nload,xbody,ipobody,nbody,cgr,
				 adc,auc,nactdof,icol,jq,irow,neq,nzl,nmethod,
				 ikmpc,ilmpc,ikboun,ilboun,
				 elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
				 ielorien,norien,orab,ntmat_,
				 t0,t0,ithermal,prestr,iprestr,vold,iperturb,sti,
				 nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
				 xstiff,npmat_,&dtime,matname,mi,ncmat_,
				 ttime,&time,istep,&iinc,ibody,ielprop,prop));
		      
	    /*  zc = damping matrix * eigenmodes */
		      
	    NNEW(zc,double,neq[1]*nev);
	    for(i=0;i<nev;i++){
	      FORTRAN(op,(&neq[1],&z[(long long)i*neq[1]],&zc[i*neq[1]],adc,auc,
			  jq,irow));
	    }
	    SFREE(adc);SFREE(auc);
		      
	    /* cco is the reduced coriolis matrix (coriolis matrix mapped onto
	       space spanned by eigenmodes) */
		      
	    NNEW(cco,double,nev*nev);
	    for(i=0;i<nev;i++){
	      for(j=0;j<=i;j++){
		for(k=0;k<neq[1];k++){
		  cco[i*nev+j]+=z[(long long)j*neq[1]+k]*zc[i*neq[1]+k];
		}
	      }
	    }
	    SFREE(zc);
		      
	    /* multiplying with the driving frequency */
		      
	    for(i=0;i<nev;i++){
	      for(j=0;j<=i;j++){
		cco[i*nev+j]*=freq[l];
	      }
	    }
		      
	    /* transpose is the negative of the matrix */
		      
	    for(i=0;i<nev;i++){
	      for(j=i+1;j<nev;j++){
		cco[i*nev+j]=-cco[j*nev+i];
	      }
	    }
		      
	    /* adding cco to cc */
		      
	    for(i=0;i<nev;i++){
	      for(j=0;j<nev;j++){
		cc[i*nev+j]+=cco[j*nev+i];
	      }
	    }
	    SFREE(cco);
	  }
	  SFREE(xstiff);
	}
	      
	/* loading for this frequency */
	      
	for(i=0;i<*nforc;i++){
	  xforcr[i]=xforcact[l**nforc+i];
	}
	      
	for(i=0;i<*nload;i++){
	  xloadr[2*i]=xloadact[l*2**nload+2*i];
	}
	      
	for(i=0;i<*nbody;i++){
	  xbodyr[7*i]=xbodyact[l*7**nbody+7*i];
	}
	      
	/* calculating the instantaneous loading vector */
	      
	rhsmain(co,nk,kon,ipkon,lakon,ne,
		ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcr,
		nforc,nelemload,sideload,xloadr,nload,xbodyr,
		ipobody,nbody,cgr,br,nactdof,&neq[1],nmethod,
		ikmpc,ilmpc,elcon,nelcon,rhcon,nrhcon,
		alcon,nalcon,alzero,ielmat,ielorien,norien,orab,ntmat_,
		t0,t1act,ithermal,iprestr,vold,iperturb,iexpl,plicon,
		nplicon,plkcon,nplkcon,
		npmat_,ttime,&time,istep,&iinc,&dtime,physcon,ibody,
		xbodyold,&reltime,veold,matname,mi,ikactmech,&nactmech,
                ielprop,prop,sti,xstateini,xstate,nstate_,ntrans,inotr,
                trab,fnext);
	  
	/* real modal coefficients */
	  
	if(!iprescribedboundary){
	  if(!cyclicsymmetry){
	    for(i=0;i<nev;i++){
	      i2=(long long)i*neq[1];
	      aa[i]=0.;
	      if(nactmech<neq[1]/2){
		for(j=0;j<nactmech;j++){
		  aa[i]+=z[i2+ikactmech[j]]*br[ikactmech[j]];
		}
	      }else{
		for(j=0;j<neq[1];j++){
		  aa[i]+=z[i2+j]*br[j];
		}
	      }
	    }
	  }else{
	    for(i=0;i<nev;i++){aa[i]=0.;}
	    for(j=0;j<nactmech;j++){
	      for(i=0;i<nev;i++){
		FORTRAN(nident,(izdof,&ikactmech[j],&nzdof,&id));
		if(id!=0){
		  if(izdof[id-1]==ikactmech[j]){
		    aa[i]+=z[(long long)i*nzdof+id-1]*br[ikactmech[j]];
		  }else{printf(" *ERROR in steadystate\n");FORTRAN(stop,());}
		}else{printf(" *ERROR in steadystate\n");FORTRAN(stop,());}
	      }
	    }
	  }
	  
	  /* imaginary modal coefficients */
		  
	  for(i=0;i<nev;i++){
	    bb[i]=0.;
	  }

	}else{

	  /* prescribed boundary conditions */
	      
	  /* next statement makes sure that br is reset to zero at the
	     start of rhs.f */

	  nactmech=neq[1];
	      
	  for(i=0;i<neq[1];i++){bi[i]=0.;}
		  
	  for(i=0;i<*nboun;i++){
	    xbounr[i]=xbounact[l**nboun+i];
	  }
		  
	  for(j=0;j<neq[1];j++){fr[j]=0.;ubr[j]=0.;}
	  for(i=0;i<*nboun;i++){
	    ic=neq[1]+i;
	    for(j=jq[ic]-1;j<jq[ic+1]-1;j++){
	      ir=irow[j]-1;
	      fr[ir]=fr[ir]-au[j]*xbounr[i];
	      ubr[ir]=fr[ir];
	    }
	  }

	  /* check for zero force (points to rigid body motion) */
		  
	  dd=0.;
	  for(j=0;j<neq[1];j++){
	    dd=dd+fr[j]*fr[j];
	  }
	  if(dd!=0.){

	    /* solve equation system */
		    
	    if(*isolver==0){
#ifdef SPOOLES
	      spooles_solve(ubr,&neq[1]);
#endif
	    }
	    else if(*isolver==4){
#ifdef SGI
	      sgi_solve(ubr,token);
#endif
	    }
	    else if(*isolver==5){
#ifdef TAUCS
	      tau_solve(ubr,&neq[1]);
#endif
	    }
	    else if(*isolver==7){
#ifdef PARDISO
	      pardiso_solve(ubr,&neq[1],&symmetryflag,&inputformat,&nrhs);
#endif
	    }
	    else if(*isolver==8){
#ifdef PASTIX
	      pastix_solve(ubr,&neq[1],&symmetryflag,&nrhs);
#endif
	    }
	  }else{
		    
	    /* rigid body mode: all displacements are equal */
		    
	    for(j=0;j<neq[1];j++){
	      ubr[j]=xbounr[0];
	    }
	  }
		  
	  FORTRAN(op,(&neq[1],ubr,mubr,adb,aub,jq,irow));
		  
	  /*	  for(i=0;i<neq[1];i++){
	    br[i]+=freq[l]*(freq[l]*mubr[i]);
	    bi[i]+=freq[l]*(-alpham*mubr[i]-betam*fr[i]);
	    }*/
	      
	  /* real and imaginary modal coefficients */
		  
	  /*	  for(i=0;i<nev;i++){
	    aa[i]=0.;
	    for(j=0;j<neq[1];j++){
	      aa[i]+=z[(long long)i*neq[1]+j]*br[j];
	    }
	  }
		  
	  for(i=0;i<nev;i++){
	    bb[i]=0.;
	    for(j=0;j<neq[1];j++){
	      bb[i]+=z[(long long)i*neq[1]+j]*bi[j];
	    }
	    }*/

	  /* calculating the eigenmodes x mass x particular solution */

	  NNEW(umubr,double,nev);

	  for(i=0;i<nev;i++){
	    for(j=0;j<neq[1];j++){
	      umubr[i]+=z[(long long)i*neq[1]+j]*mubr[j];
	    }
	  }
	  
	  /* real and imaginary modal coefficients */

	  for(i=0;i<nev;i++){
	    aa[i]=0.;
	    for(j=0;j<neq[1];j++){
	      aa[i]+=z[(long long)i*neq[1]+j]*br[j];
	    }
	    aa[i]+=freq[l]*(freq[l]*umubr[i]);
	  }
	      
	  for(i=0;i<nev;i++){
	    bb[i]=0.;
	    for(j=0;j<neq[1];j++){
	      bb[i]+=z[(long long)i*neq[1]+j]*bi[j];
	    }
	    bb[i]+=freq[l]*(-fric[i]*umubr[i]);
	  }
	  SFREE(umubr);
	}
	      
	/* calculating the modal coefficients */
	      
	if(dampmat==0){
	  for(i=0;i<nev;i++){
	    dd=pow(d[i]-pow(freq[l],2),2)+
	      pow(fric[i],2)*pow(freq[l],2);
	    bjr[i]=(aa[i]*(d[i]-freq[l]*freq[l])+
		    bb[i]*fric[i]*freq[l])/dd;
	    bji[i]=(bb[i]*(d[i]-freq[l]*freq[l])-
		    aa[i]*fric[i]*freq[l])/dd;
	  }
	}else{
	  nev2=2*nev;
	  NNEW(am,double,nev2*nev2);
	  NNEW(bm,double,nev2);
	  NNEW(ipiv,ITG,nev2);
		  
	  for(i=0;i<nev2;i++){
	    for(j=0;j<nev2;j++){
	      am[i*nev2+j]=0.;
	    }
	    bm[i]=0.;
	  }
	  for(i=0;i<nev;i++){
	    am[i*nev2+i]=d[i]-freq[l]*freq[l];
	    am[(i+nev)*nev2+i]=-fric[i]*freq[l];
	    bm[i]=aa[i];
	    am[i*nev2+nev+i]=-am[(i+nev)*nev2+i];
	    am[(i+nev)*nev2+nev+i]=am[i*nev2+i];
	    bm[nev+i]=bb[i];
	    for(j=0;j<nev;j++){
	      am[(j+nev)*nev2+i]=am[(j+nev)*nev2+i]
		-cc[i*nev+j];//*freq[l];
	      am[j*nev2+nev+i]=am[j*nev2+nev+i]
		+cc[i*nev+j];//*freq[l];
	    }
	  }
		  
	  /* solving the system of equations */
		  
	  FORTRAN(dgesv,(&nev2,&nrhs,am,&nev2,ipiv,bm,&nev2,&info));
	  if(info!=0){
	    printf(" *ERROR in steadystate: fatal termination of dgesv\n");
	    printf("       info=%" ITGFORMAT "\n",info);
	    /*		  FORTRAN(stop,());*/
	  }
		  
	  /* storing the solution in bjr and bji */
		  
	  for(i=0;i<nev;i++){
	    bjr[i]=bm[i];
	    bji[i]=bm[nev+i];
	  }
		  
	  SFREE(am);SFREE(bm);SFREE(ipiv);
	}
      
	/* calculating the real response */
	      
	if(iprescribedboundary){
	  if(nmdnode==0){
	    memcpy(&br[0],&ubr[0],sizeof(double)*neq[1]);
	  }else{
	    for(i=0;i<nmddof;i++){
	      br[imddof[i]]=ubr[imddof[i]];
	    }
	  }
	}
	else{
	  if(nmdnode==0){
	    DMEMSET(br,0,neq[1],0.);
	  }else{
	    for(i=0;i<nmddof;i++){
	      br[imddof[i]]=0.;
	    }
	  }
	}
	      
	if(!cyclicsymmetry){
	  if(nmdnode==0){
	    for(i=0;i<neq[1];i++){
	      for(j=0;j<nev;j++){
		br[i]+=bjr[j]*z[(long long)j*neq[1]+i];
	      }
	    }
	  }else{
	    for(i=0;i<nmddof;i++){
	      for(j=0;j<nev;j++){
		br[imddof[i]]+=bjr[j]*z[(long long)j*neq[1]+imddof[i]];
	      }
	    }
	  }
	}else{
	  for(i=0;i<nmddof;i++){
	    FORTRAN(nident,(izdof,&imddof[i],&nzdof,&id));
	    if(id!=0){
	      if(izdof[id-1]==imddof[i]){
		for(j=0;j<nev;j++){
		  br[imddof[i]]+=bjr[j]*z[(long long)j*nzdof+id-1];
		}
	      }else{printf(" *ERROR in steadystate\n");FORTRAN(stop,());}
	    }else{printf(" *ERROR in steadystate\n");FORTRAN(stop,());}
	  }
	}
	  
	/* calculating the imaginary response */
	      
	if(nmdnode==0){
	  DMEMSET(bi,0,neq[1],0.);
	}else{
	  for(i=0;i<nmddof;i++){
	    bi[imddof[i]]=0.;
	  }
	}
	  
	if(!cyclicsymmetry){
	  if(nmdnode==0){
	    for(i=0;i<neq[1];i++){
	      for(j=0;j<nev;j++){
		bi[i]+=bji[j]*z[(long long)j*neq[1]+i];
	      }
	    }
	  }else{
	    for(i=0;i<nmddof;i++){
	      for(j=0;j<nev;j++){
		bi[imddof[i]]+=bji[j]*z[(long long)j*neq[1]+imddof[i]];
	      }
	    }
	  }
	}else{
	  for(i=0;i<nmddof;i++){
	    FORTRAN(nident,(izdof,&imddof[i],&nzdof,&id));
	    if(id!=0){
	      if(izdof[id-1]==imddof[i]){
		for(j=0;j<nev;j++){
		  bi[imddof[i]]+=bji[j]*z[(long long)j*nzdof+id-1];
		}
	      }else{printf(" *ERROR in steadystate\n");FORTRAN(stop,());}
	    }else{printf(" *ERROR in steadystate\n");FORTRAN(stop,());}
	  }
	}
	      
	if(nmdnode==0){
	      
	  /* magnitude and phase of the response */

	  for(i=0;i<neq[1];i++){
	    breal=br[i];
	    br[i]=sqrt(br[i]*br[i]+bi[i]*bi[i]);
	    if(fabs(breal)<1.e-10){
	      if(bi[i]>0.){bi[i]=pi/2.;}
	      else{bi[i]=-pi/2.;}
	    }
	    else{
	      bi[i]=atan(bi[i]/breal);
	      if(breal<0.){bi[i]+=pi;}
	    }
	  }
		  
	  /* correction for the sinus terms */
		  
	  if((l!=0)&&(2*(ITG)floor(l/2.+0.1)==l)){
	    for(i=0;i<neq[1];i++){
	      //			  bi[i]-=pi/2.;}
	      bi[i]+=pi/2.;}
	  }
		  
	  /* contribution to the time response */
		  
	  for(j=0;j<nfour;j++){
	    time=tmin+2.*pi/freqnh[ii]*(double)j/(double)nfour;
	    for(i=0;i<neq[1];i++){
	      btot[j*neq[1]+i]+=br[i]*cos(freq[l]*time+bi[i]);
	      bp[j*neq[1]+i]-=freq[l]*br[i]*sin(freq[l]*time+bi[i]);
	    }
	  }
	}else{
	      
	  /* magnitude and phase of the response */

	  for(jj=0;jj<nmddof;jj++){
	    i=imddof[jj];
	    breal=br[i];
	    br[i]=sqrt(br[i]*br[i]+bi[i]*bi[i]);
	    if(fabs(breal)<1.e-10){
	      if(bi[i]>0.){bi[i]=pi/2.;}
	      else{bi[i]=-pi/2.;}
	    }
	    else{
	      bi[i]=atan(bi[i]/breal);
	      if(breal<0.){bi[i]+=pi;}
	    }
	  }
		  
	  /* correction for the sinus terms */
		  
	  if((l!=0)&&(2*(ITG)floor(l/2.+0.1)==l)){
	    for(jj=0;jj<nmddof;jj++){
	      i=imddof[jj];
	      bi[i]+=pi/2.;}
	  }
		  
	  /* contribution to the time response */
		  
	  for(j=0;j<nfour;j++){
	    time=tmin+2.*pi/freqnh[ii]*(double)j/(double)nfour;
	    for(jj=0;jj<nmddof;jj++){
	      i=imddof[jj];
	      btot[j*neq[1]+i]+=br[i]*cos(freq[l]*time+bi[i]);
	      bp[j*neq[1]+i]-=freq[l]*br[i]*sin(freq[l]*time+bi[i]);
	    }
	  }
	}

	/* resetting the part of br occupied by the variables to be printed
	   to zero */

	if(!iprescribedboundary){
	  if(nmdnode==0){
	    DMEMSET(br,0,neq[1],0.);
	  }else{
	    for(i=0;i<nmddof;i++){
	      br[imddof[i]]=0.;
	    }
	  }
	}
	      
      }
	  
      SFREE(xforcr);SFREE(xloadr);SFREE(xbodyr);SFREE(br);SFREE(bi);SFREE(freq);
      SFREE(bjr);SFREE(bji);SFREE(aa);SFREE(bb);

      if(iprescribedboundary) {SFREE(xbounr);SFREE(fr);SFREE(ubr);SFREE(mubr);}
	  
      /* result fields */
	  
      NNEW(vr,double,mt**nk);

      if(intpointvar==1){
	NNEW(fn,double,mt**nk);
	NNEW(stn,double,6**nk);
	NNEW(stx,double,6*mi[0]**ne);

	if(*ithermal>1) {
	  NNEW(qfn,double,3**nk);
	  NNEW(qfx,double,3*mi[0]**ne);}
	      
	if(strcmp1(&filab[261],"E   ")==0) NNEW(een,double,6**nk);
	if(strcmp1(&filab[522],"ENER")==0) NNEW(enern,double,*nk);
	if(strcmp1(&filab[2697],"ME  ")==0) NNEW(emn,double,6**nk);
	      
	NNEW(eei,double,6*mi[0]**ne);
	if(*nener==1){
	  NNEW(stiini,double,6*mi[0]**ne);
	  NNEW(emeini,double,6*mi[0]**ne);
	  NNEW(enerini,double,2*mi[0]**ne);}
      }
	  
      /* storing the results */
	  
      for(l=0;l<nfour;l++){
	time=tmin+2.*pi/freqnh[ii]*(double)l/(double)nfour;
	ptime=time;
	      
	if(nmdnode==0){
	  DMEMSET(vr,0,mt**nk,0.);
	}else{
	  for(jj=0;jj<nmdnode;jj++){
	    i=imdnode[jj]-1;
	    for(j=1;j<4;j++){
	      vr[mt*i+j]=0.;
	    }
	  }
	}
      
	/* calculating displacements/temperatures */

	*nmethod=4;
	FORTRAN(dynresults,(nk,vr,ithermal,nactdof,vold,nodeboun,
			    ndirboun,&xbounacttime[l**nboun],nboun,ipompc,nodempc,
			    coefmpc,labmpc,nmpc,&btot[l*neq[1]],&bp[l*neq[1]],veold,&dtime,mi,
			    imdnode,&nmdnode,imdboun,&nmdboun,imdmpc,&nmdmpc,nmethod,&time));
	*nmethod=5;

	results(co,nk,kon,ipkon,lakon,ne,vr,stn,inum,
		stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		ielmat,ielorien,norien,orab,ntmat_,t0,t1,
		ithermal,prestr,iprestr,filab,eme,emn,een,
		iperturb,f,fn,nactdof,&iout,qa,
		vold,&btot[l*neq[1]],nodeboun,ndirboun,
		&xbounacttime[l**nboun],nboun,
		ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],
		veold,accold,&bet,&gam,&dtime,&time,&xnull,
		plicon,nplicon,plkcon,nplkcon,
		xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
		&icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,
		enern,emeini,xstaten,eei,enerini,cocon,ncocon,
		set,nset,istartset,iendset,ialset,nprint,prlab,prset,
		qfx,qfn,trab,inotr,ntrans,fmpc,nelemload,nload,ikmpc,
		ilmpc,istep,&iinc,springarea,&reltime,&ne0,
		thicke,shcon,nshcon,
		sideload,xload,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		&mortar,islavact,cdn,islavnode,nslavnode,&ntie,clearini,
		islavsurf,ielprop,prop,energyini,energy,&iit,iponoel,
		inoel,nener,orname,&network,ipobody,xbodyact,ibody,
		typeboun,itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
		islavelinv,autloc,irowtloc,jqtloc,&nboun2,
		ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
		labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
		&intscheme,physcon);
	  
	(*kode)++;
	mode=-1;
	      
	if(strcmp1(&filab[1044],"ZZS")==0){
	  NNEW(neigh,ITG,40**ne);
	  NNEW(ipneigh,ITG,*nk);
	}

	frd(co,&nkg,kon,ipkon,lakon,&neg,vr,stn,inum,nmethod,
	    kode,filab,een,t1,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	    nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,&neg,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	    thicke,jobnamec,output,qfx,cdn,&mortar,cdnr,cdni,nmat,
	    ielprop,prop,sti);

	if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);}

      }
	  
      SFREE(vr);SFREE(btot);SFREE(bp);

      if(intpointvar==1){
	SFREE(fn);SFREE(stn);SFREE(stx);SFREE(eei);
	if(*ithermal>1) {SFREE(qfn);SFREE(qfx);}
	      
	if(strcmp1(&filab[261],"E   ")==0) SFREE(een);
	if(strcmp1(&filab[522],"ENER")==0) SFREE(enern);
	if(strcmp1(&filab[2697],"ME  ")==0) SFREE(emn);
	      
	if(*nener==1){SFREE(stiini);SFREE(emeini);SFREE(enerini);}
      }
	  
    }
    SFREE(xforcact);SFREE(xloadact);SFREE(xbodyact);SFREE(xbounact);
    SFREE(xbounacttime);SFREE(freqnh);
    if(*ithermal==1) SFREE(t1act);
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
    }

    SFREE(ikactmech);

  }

  SFREE(adb);SFREE(aub);SFREE(z);SFREE(d);SFREE(inum);

  if(!cyclicsymmetry){
    SFREE(ad);SFREE(au);
  }else{
    SFREE(izdof);SFREE(nm);

    *nk/=nsectors;
    *ne/=nsectors;
    *nboun/=nsectors;
    neq[1]=neq[1]*2/nsectors;

    RENEW(co,double,3**nk);
    if(*ithermal!=0){
      RENEW(t0,double,*nk);
      RENEW(t1old,double,*nk);
      RENEW(t1,double,*nk);
      if(*nam>0) RENEW(iamt1,ITG,*nk);
    }
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
    if(*nener==1)RENEW(ener,double,2*mi[0]**ne);

    /* distributed loads */

    for(i=0;i<*nload;i++){
      if(nelemload[2*i]<=*ne*nsectors){
	nelemload[2*i]-=*ne*nelemload[2*i+1];
      }else{
	nelemload[2*i]-=*ne*(nsectors+nelemload[2*i+1]-1);
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
      if(nodeforc[2*i]<=*nk*nsectors){
	nodeforc[2*i]-=*nk*nodeforc[2*i+1];
      }else{
	nodeforc[2*i]-=*nk*(nsectors+nodeforc[2*i+1]-1);
      }
    }
  }

  SFREE(fric);

  if(dampmat||coriolis){SFREE(cc);}

  /* reactivating centrigual forces */

  for(i=0;i<*nbody;i++){
    if((ibody[3*i]==0)||(ibody[3*i]==4)){
      ibody[3*i]=1;
    }
  }

  if(nherm!=1){SFREE(xmr);SFREE(xmi);}

  SFREE(imddof);SFREE(imdnode);SFREE(imdboun);SFREE(imdmpc);SFREE(imdelem);

  *cop=co;*konp=kon;*ipkonp=ipkon;*lakonp=lakon;*ielmatp=ielmat;
  *ielorienp=ielorien;*inotrp=inotr;*nodebounp=nodeboun;
  *ndirbounp=ndirboun;*iambounp=iamboun;*xbounp=xboun;*veoldp=veold;
  *xbounoldp=xbounold;*ikbounp=ikboun;*ilbounp=ilboun;*nactdofp=nactdof;
  *voldp=vold;*emep=eme;*enerp=ener;*ipompcp=ipompc;*nodempcp=nodempc;
  *coefmpcp=coefmpc;*labmpcp=labmpc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
  *fmpcp=fmpc;*iamt1p=iamt1;*t0p=t0;*t1oldp=t1old;*t1p=t1;

  //  (*ttime)+=(*tper);

  return;
}
