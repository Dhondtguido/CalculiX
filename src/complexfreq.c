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

void complexfreq(double **cop,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,ITG *ne,
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
		 double *xstate,ITG *npmat_,char *matname,ITG *mi,
		 ITG *ncmat_,ITG *nstate_,double **enerp,char *jobnamec,
		 double *ttime,char *set,ITG *nset,ITG *istartset,
		 ITG *iendset,ITG **ialsetp,ITG *nprint,char *prlab,
		 char *prset,ITG *nener,double *trab,
		 ITG **inotrp,ITG *ntrans,double **fmpcp,ITG *ipobody,
		 ITG *ibody,
		 double *xbody,ITG *nbody,double *xbodyold,ITG *istep,
		 ITG *isolver,ITG *jq,char *output,ITG *mcs,ITG *nkon,
		 ITG *mpcend,ITG *ics,double *cs,ITG *ntie,char *tieset,
		 ITG *idrct,ITG *jmax,
		 double *ctrl,ITG *itpamp,double *tietol,ITG *nalset,
		 ITG *ikforc,ITG *ilforc,double *thicke,
		 char *jobnamef,ITG *mei,ITG *nmat,ITG *ielprop,double *prop,
		 char *orname,char *typeboun,double *t0g,double *t1g){

  char fneig[132]="",description[13]="            ",*lakon=NULL,*labmpc=NULL,
    *lakont=NULL,*turdir=NULL,*labmpc2=NULL;

  ITG nev,i,j,k,idof,*inum=NULL,id,
    iinc=0,l,iout=1,ielas,icmd=3,ifreebody,mode,m,nherm,
    *kon=NULL,*ipkon=NULL,*ielmat=NULL,*ielorien=NULL,*islavact=NULL,
    *inotr=NULL,*nodeboun=NULL,*ndirboun=NULL,*iamboun=NULL,*ikboun=NULL,
    *ilboun=NULL,*nactdof=NULL,*ipompc=NULL,*nodempc=NULL,*ikmpc=NULL,
    *ilmpc=NULL,nsectors,nmd,nevd,*nm=NULL,*iamt1=NULL,*islavnode=NULL,
    ngraph=1,nkg,neg,ne0,ij,lprev,nope,indexe,ilength,*nslavnode=NULL,
    *ipneigh=NULL,*neigh=NULL,index,im,cyclicsymmetry,inode,
    *ialset=*ialsetp,mt=mi[1]+1,kmin,kmax,i1,iit=-1,network=0,
    *iter=NULL,lint,lfin,kk,kkv,kk6,kkx,icomplex,igeneralizedforce,
    idir,*inumt=NULL,icntrl,imag,jj,is,l1,*inocs=NULL,ml1,l2,nkt,net,
    *ipkont=NULL,*ielmatt=NULL,*inotrt=NULL,*kont=NULL,node,iel,*ielcs=NULL,
    ielset,*istartnmd=NULL,*iendnmd=NULL,inmd,neqact,*nshcon=NULL,
    *ipev=NULL,icfd=0,*inomat=NULL,mortar=0,*islavsurf=NULL,
    *iponoel=NULL,*inoel=NULL,iperturbsav,nevcomplex,*itiefac=NULL,
    mscalmethod=0,*islavelinv=NULL,*irowtloc=NULL,*jqtloc=NULL,nboun2,
    *ndirboun2=NULL,*nodeboun2=NULL,nmpc2,*ipompc2=NULL,*nodempc2=NULL,
    *ikboun2=NULL,*ilboun2=NULL,*ikmpc2=NULL,*ilmpc2=NULL,mortartrafoflag=0,
    intscheme=0;

  long long i2;

  double *d=NULL,*z=NULL,*stiini=NULL,*cc=NULL,*v=NULL,*zz=NULL,*emn=NULL,
    *stn=NULL,*stx=NULL,*een=NULL,*adb=NULL,*xstiff=NULL,*cdn=NULL,
    *aub=NULL,*f=NULL,*fn=NULL,*epn=NULL,*xstateini=NULL,
    *enern=NULL,*xstaten=NULL,*eei=NULL,*enerini=NULL,*qfn=NULL,
    *qfx=NULL,*cgr=NULL,*au=NULL,dtime,reltime,*t0=NULL,*t1=NULL,*t1old=NULL,
    sum,qa[4],cam[5],accold[1],bet,gam,*ad=NULL,alpham,betam,
    *co=NULL,*xboun=NULL,*xbounold=NULL,*vold=NULL,*emeini=NULL,
    *eme=NULL,*ener=NULL,*coefmpc=NULL,*fmpc=NULL,*veold=NULL,
    *adc=NULL,*auc=NULL,*zc=NULL,*fnr=NULL,*fni=NULL,setnull,deltmx,dd,
    theta,*vini=NULL,*vr=NULL,*vi=NULL,*stnr=NULL,*stni=NULL,*vmax=NULL,
    *stnmax=NULL,*cstr=NULL,*sti=*stip,time0=0.0,time=0.0,zero=0.0,
    *springarea=NULL,*eenmax=NULL,*aa=NULL,*bb=NULL,*xx=NULL,
    *eiga=NULL,*eigb=NULL,*eigxx=NULL,*temp=NULL,*coefmpcnew=NULL,xreal,
    ximag,t[3],*vt=NULL,*t1t=NULL,*stnt=NULL,*eent=NULL,*fnt=NULL,*enernt=NULL,
    *stxt=NULL,pi,ctl,stl,*cot=NULL,*qfnt=NULL,vreal,vimag,constant,stnreal,
    stnimag,freq,*emnt=NULL,*shcon=NULL,*eig=NULL,*clearini=NULL,
    *eigxr=NULL,*eigxi=NULL,*xmac=NULL,*bett=NULL,*betm=NULL,*xmaccpx=NULL,
    fmin=0.,fmax=1.e30,*xmr=NULL,*xmi=NULL,*zi=NULL,*eigx=NULL,
    *pslavsurf=NULL,*pmastsurf=NULL,*cdnr=NULL,*cdni=NULL,*tinc,*tper,
    *tmin,*tmax,*energyini=NULL,*energy=NULL,e1[3],e2[3],xn[3],*smscale=NULL,
    *autloc=NULL,*xboun2=NULL,*coefmpc2=NULL,*physcon=NULL;

  FILE *f1;

#ifdef SGI
  ITG token;
#endif

  pi=4.*atan(1.);
  constant=180./pi;

  co=*cop;kon=*konp;ipkon=*ipkonp;lakon=*lakonp;ielmat=*ielmatp;
  ielorien=*ielorienp;inotr=*inotrp;nodeboun=*nodebounp;
  ndirboun=*ndirbounp;iamboun=*iambounp;xboun=*xbounp;
  xbounold=*xbounoldp;ikboun=*ikbounp;ilboun=*ilbounp;nactdof=*nactdofp;
  vold=*voldp;eme=*emep;ener=*enerp;ipompc=*ipompcp;nodempc=*nodempcp;
  coefmpc=*coefmpcp;labmpc=*labmpcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
  fmpc=*fmpcp;veold=*veoldp;iamt1=*iamt1p;t0=*t0p;t1=*t1p;t1old=*t1oldp;

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

  NNEW(xstiff,double,(long long)27*mi[0]**ne);

  dtime=*tinc;

  alpham=xmodal[0];
  betam=xmodal[1];

  dd=ctrl[16];deltmx=ctrl[26];

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
    printf(" *ERROR in complexfreq: cannot open eigenvalue file for reading");
    exit(0);
  }

  printf(" *INFO in complexfreq: if there are problems reading the .eig file this may be due to:\n");
  printf("       1) the nonexistence of the .eig file\n");
  printf("       2) other boundary conditions than in the input deck\n");
  printf("          which created the .eig file\n\n");

  if(fread(&cyclicsymmetry,sizeof(ITG),1,f1)!=1){
    printf(" *ERROR in complexfreq reading the cyclic symmetry flag in the eigenvalue file");
    exit(0);
  }

  if(fread(&nherm,sizeof(ITG),1,f1)!=1){
    printf(" *ERROR in complexfreq reading the Hermitian flag in the eigenvalue file");
    exit(0);
  }

  if(nherm!=1){
    printf(" *ERROR in complexfreq: the eigenvectors in the .eig-file result\n");
    printf("        from a non-Hermitian eigenvalue problem. The complex\n");
    printf("        frequency procedure cannot handle that yet\n\n");
    FORTRAN(stop,());
  }

  iperturbsav=iperturb[0];
  if(fread(&iperturb[0],sizeof(ITG),1,f1)!=1){
    printf(" *ERROR in complexfreq reading the perturbation flag in the eigenvalue file");
    exit(0);
  }
  if(iperturbsav!=iperturb[0]){
    if(iperturb[0]==0){
      printf(" *ERROR in complexfreq: the .eig-file was created without perturbation info\n");
      printf("        the present *COMPLEX FREQUENCY step, however, contains the\n");
      printf("        perturbation parameter on the *STEP card.\n");
    }else if(iperturb[0]==1){
      printf(" *ERROR in complexfreq: the .eig-file was created with perturbation info\n");
      printf("        the present *COMPLEX FREQUENCY step, however, does not contain the\n");
      printf("        perturbation parameter on the *STEP card.\n");
    }
    FORTRAN(stop,());
  }

  if(iperturb[0]==1){
    if(fread(vold,sizeof(double),mt**nk,f1)!=mt**nk){
      printf(" *ERROR in complexfreq reading the reference displacements in the eigenvalue file...");
      exit(0);
    }
  }

  nsectors=1;

  if(!cyclicsymmetry){

    nkg=*nk;
    neg=*ne;

    if(fread(&nev,sizeof(ITG),1,f1)!=1){
      printf(" *ERROR in complexfreq reading the number of eigenvalues in the eigenvalue file...");
      exit(0);
    }
      

    if(nherm==1){
      NNEW(d,double,nev);
      if(fread(d,sizeof(double),nev,f1)!=nev){
	printf(" *ERROR in complexfreq reading the eigenvalues in the eigenvalue file...");
	exit(0);
      }
      //	  for(i=0;i<nev;i++){printf("eigenvalue %d %e\n",i,d[i]);}
    }else{
      NNEW(d,double,2*nev);
      if(fread(d,sizeof(double),2*nev,f1)!=2*nev){
	printf(" *ERROR in complexfreq reading the eigenvalues in the eigenvalue file...");
	exit(0);
      }
    }
      
    NNEW(ad,double,neq[1]);
    NNEW(adb,double,neq[1]);
    NNEW(au,double,nzs[2]);
    NNEW(aub,double,nzs[1]);
      
    /* reading the stiffness matrix */

    if(fread(ad,sizeof(double),neq[1],f1)!=neq[1]){
      printf(" *ERROR in complexfreq reading the diagonal of the stiffness matrix in the eigenvalue file...");
      exit(0);
    }
      
    if(fread(au,sizeof(double),nzs[2],f1)!=nzs[2]){
      printf(" *ERROR in complexfreq reading the off-diagonal terms of the stiffness matrix in the eigenvalue file...");
      exit(0);
    }
      
    /* reading the mass matrix */

    if(fread(adb,sizeof(double),neq[1],f1)!=neq[1]){
      printf(" *ERROR in complexfreq reading the diagonal of the mass matrix in eigenvalue file...");
      exit(0);
    }
      
    if(fread(aub,sizeof(double),nzs[1],f1)!=nzs[1]){
      printf(" *ERROR in complexfreq reading the off-diagonals of the mass matrix in the eigenvalue file...");
      exit(0);
    }
      
    /* reading the eigenvectors */

    if(nherm==1){
      NNEW(z,double,neq[1]*nev);
      if(fread(z,sizeof(double),neq[1]*nev,f1)!=neq[1]*nev){
	printf(" *ERROR in complexfreq reading the eigenvectors in the eigenvalue file...");
	exit(0);
      }
    }else{
      NNEW(z,double,2*neq[1]*nev);
      if(fread(z,sizeof(double),2*neq[1]*nev,f1)!=2*neq[1]*nev){
	printf(" *ERROR in complexfreq reading the eigenvectors in the eigenvalue file...");
	exit(0);
      }
    }

    NNEW(nm,ITG,nev);
    for(i=0;i<nev;i++){nm[i]=-1;}
  }
  else{

    if(*nmethod==6){
      printf(" *ERROR in complexfreq: Coriolis forces cannot\n");
      printf("       be combined with cyclic symmetry\n\n");
      FORTRAN(stop,());
    }

    nev=0;
    do{
      if(fread(&nmd,sizeof(ITG),1,f1)!=1){
	break;
      }
      if(fread(&nevd,sizeof(ITG),1,f1)!=1){
	printf(" *ERROR in complexfreq reading the number of eigenvalues in the eigenvalue file...");
	exit(0);
      }

      /* reading the eigenvalues */

      if(nev==0){
	if(nherm==1){NNEW(d,double,nevd);
	}else{NNEW(d,double,2*nevd);}
	NNEW(nm,ITG,nevd);
      }else{
	printf(" *ERROR in complexfreq: flutter forces cannot\n");
	printf("       be combined with multiple modal diameters\n");
	printf("       in cyclic symmetry calculations\n\n");
	FORTRAN(stop,());
      }
	  
      if(nherm==1){
	if(fread(&d[nev],sizeof(double),nevd,f1)!=nevd){
	  printf(" *ERROR in complexfreq reading the eigenvalues in the eigenvalue file...");
	  exit(0);
	}
      }else{
	if(fread(&d[nev],sizeof(double),2*nevd,f1)!=2*nevd){
	  printf(" *ERROR in complexfreq reading the eigenvalues in the eigenvalue file...");
	  exit(0);
	}
      }
      for(i=nev;i<nev+nevd;i++){nm[i]=nmd;}
	  
      /* reading the stiffness and the mass matrix */

      if(nev==0){
	    
	NNEW(ad,double,neq[1]);
	NNEW(au,double,nzs[1]);

	if(fread(ad,sizeof(double),neq[1],f1)!=neq[1]){
	  printf(" *ERROR in complexfreq reading the diagonal of the stiffness matrix in the eigenvalue file...");
	  exit(0);
	}
	      
	if(fread(au,sizeof(double),nzs[1],f1)!=nzs[1]){
	  printf(" *ERROR in complexfreq reading the off-diagonals of the stiffness matrix in the eigenvalue file...");
	  exit(0);
	}

	SFREE(ad);SFREE(au);
	      
	NNEW(adb,double,neq[1]);
	NNEW(aub,double,nzs[1]);

	if(fread(adb,sizeof(double),neq[1],f1)!=neq[1]){
	  printf(" *ERROR in complexfreq reading the diagonal of the mass matrix in the eigenvalue file...");
	  exit(0);
	}
	      
	if(fread(aub,sizeof(double),nzs[1],f1)!=nzs[1]){
	  printf(" *ERROR in complexfreq reading the off-diagonals of the mass matrix in the eigenvalue file...");
	  exit(0);
	}
      }
	  
      /* reading the eigenmodes */
 
      if(nev==0){
	NNEW(z,double,neq[1]*nevd);
      }else{
	RENEW(z,double,(long long)neq[1]*(nev+nevd));
      }
	  
      if(fread(&z[neq[1]*nev],sizeof(double),neq[1]*nevd,f1)!=neq[1]*nevd){
	printf(" *ERROR in complexfreq reading eigenvectors in the eigenvalue file...");
	exit(0);
      }
      nev+=nevd;
    }while(1);

    /* removing double eigenmodes */

    j=-1;
    for(i=0;i<nev;i++){
      if((i/2)*2==i){
	j++;
	if(nherm==1){
	  d[j]=d[i];
	}else{
	  d[2*j]=d[2*i];d[2*j+1]=d[2*i+1];
	}
	nm[j]=nm[i];
	for(k=0;k<neq[1];k++){
	  z[j*neq[1]+k]=z[i*neq[1]+k];
	}
      }
    }
    nev=j+1;
    if(nherm==1){RENEW(d,double,nev);}else{RENEW(d,double,2*nev);}
    RENEW(nm,ITG,nev);
    RENEW(z,double,neq[1]*nev);

    /* determining the maximum amount of segments */

    for(i=0;i<*mcs;i++){
      if(cs[17*i]>nsectors) nsectors=(ITG)(cs[17*i]+0.5);
    }

    /* determining the maximum number of sectors to be plotted */

    for(j=0;j<*mcs;j++){
      if(cs[17*j+4]>ngraph) ngraph=(ITG)cs[17*j+4];
    }
    nkg=*nk*ngraph;
    neg=*ne*ngraph;

  }

  fclose(f1);

  /* assigning nodes and elements to sectors */

  if(cyclicsymmetry){
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
	  else {nope=6;}
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
  }
  
  /* check for rigid body modes 
     if there is a jump of 1.e4 in two subsequent eigenvalues
     all eigenvalues preceding the jump are considered to
     be rigid body modes and their frequency is set to zero
     deactivated for Hermitian matrices on 20190215 */

  if(nherm==1){
    //      setnull=1.;
    //      for(i=0;i<nev;i++){d[i]=fabs(d[i]);}
    //      for(i=nev-2;i>-1;i--){
    //	  if(fabs(d[i])<10.) d[i]=0.;
    //	  if(fabs(d[i])<0.0001*fabs(d[i+1])) setnull=0.;
    //	  d[i]*=setnull;
    //	  }
    //	  for(i=0;i<nev;i++){printf("eigenvalue %d %e\n",i,d[i]);}
  }else{
    setnull=1.;
    for(i=nev-2;i>-1;i--){
      if(sqrt(d[2*i]*d[2*i]+d[2*i+1]*d[2*i+1])<
	 0.0001*sqrt(d[2*i+2]*d[2*i+2]+d[2*i+3]*d[2*i+3])) setnull=0.;
      d[2*i]*=setnull;
      d[2*i+1]*=setnull;
    }
  }

  /* determining the frequency ranges corresponding to one
     and the same nodal diameter */

  if(cyclicsymmetry){
    NNEW(istartnmd,ITG,nev);
    NNEW(iendnmd,ITG,nev);
    nmd=0;
    inmd=nm[0];
    istartnmd[0]=1;
    for(i=1;i<nev;i++){
      if(nm[i]==inmd) continue;
      iendnmd[nmd]=i;
      nmd++;
      istartnmd[nmd]=i+1;
      inmd=nm[i];
    }
    iendnmd[nmd]=nev;
    nmd++;
    RENEW(istartnmd,ITG,nmd);
    RENEW(iendnmd,ITG,nmd);
  }

  if(*nmethod==6){

    /* Coriolis */

    neqact=neq[1];

    /* converting centrifugal force in Coriolis force */

    for(k=0;k<*nbody;k++){
      if(ibody[3*k]==1) ibody[3*k]=4;
    }

    /* assigning the body forces to the elements */ 

    /*    ifreebody=*ne+1;
    NNEW(ipobody,ITG,2**ne);
    for(k=1;k<=*nbody;k++){
      FORTRAN(bodyforce,(cbody,ibody,ipobody,nbody,set,istartset,
			 iendset,ialset,&inewton,nset,&ifreebody,&k));
      RENEW(ipobody,ITG,2*(*ne+ifreebody));
    }
    RENEW(ipobody,ITG,2*(ifreebody-1));*/
      
    if(cyclicsymmetry){
      printf(" *ERROR in complexfreq: Coriolis is not allowed in combination with cyclic symmetry\n");
      FORTRAN(stop,());
    }

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
			 ttime,&time0,istep,&iinc,ibody,ielprop,prop));

    /*  zc = damping matrix * eigenmodes */

    NNEW(zc,double,neq[1]*nev);
    for(i=0;i<nev;i++){
      FORTRAN(op_corio,(&neq[1],&z[(long long)i*neq[1]],&zc[i*neq[1]],adc,auc,
			jq,irow));
    }
    SFREE(adc);SFREE(auc);

    /* cc is the reduced damping matrix (damping matrix mapped onto
       space spanned by eigenmodes) */

    NNEW(cc,double,nev*nev);
    for(i=0;i<nev;i++){
      for(j=0;j<=i;j++){
	for(k=0;k<neq[1];k++){
	  cc[i*nev+j]+=z[(long long)j*neq[1]+k]*zc[i*neq[1]+k];
	}
      }
    }

    /* symmetric part of cc matrix */

    for(i=0;i<nev;i++){
      for(j=i+1;j<nev;j++){
	cc[i*nev+j]=-cc[j*nev+i];
      }
    }
    SFREE(zc);

    /* solving for the complex eigenvalues */

    nevcomplex=mei[0];

    NNEW(aa,double,4*nev*nev);
    NNEW(bb,double,4*nev*nev);
    NNEW(xx,double,4*nev*nevcomplex);
    NNEW(temp,double,4*nev*nev);
    NNEW(eiga,double,2*nev);
    NNEW(eigb,double,2*nev);
    NNEW(eigxx,double,2*nevcomplex);
    NNEW(iter,ITG,nev);

    FORTRAN(coriolissolve,(cc,&nev,aa,bb,xx,eiga,eigb,eigxx,
			   iter,d,temp,&nevcomplex));
      
    SFREE(aa);SFREE(bb);SFREE(temp);SFREE(eiga);SFREE(eigb);SFREE(iter);SFREE(cc);

  }else{

    /* flutter */

    /* complex force is being read (e.g. due to fluid flow) */

    if(cyclicsymmetry){
      neqact=neq[1]/2;
    }else{
      neqact=neq[1];
    }
    NNEW(zc,double,2*neqact*nev);
    NNEW(aa,double,4*nev*nev);

    FORTRAN(readforce,(zc,&neqact,&nev,nactdof,ikmpc,nmpc,
		       ipompc,nodempc,mi,coefmpc,jobnamef,
		       aa,&igeneralizedforce));

    NNEW(bb,double,4*nev*nev);
    NNEW(xx,double,4*nev*nev);
    NNEW(eiga,double,2*nev);
    NNEW(eigb,double,2*nev);
    NNEW(eigxx,double,2*nev);
    NNEW(iter,ITG,nev);
    FORTRAN(forcesolve,(zc,&nev,aa,bb,xx,eiga,eigb,eigxx,iter,d,
			&neq[1],z,istartnmd,iendnmd,&nmd,&cyclicsymmetry,
			&neqact,&igeneralizedforce));
    SFREE(aa);SFREE(bb);SFREE(eiga);SFREE(eigb);SFREE(iter);SFREE(zc);

    nevcomplex=nev;
      
  }

  /* sorting the eigenvalues and eigenmodes according to the size of the
     eigenvalues */
      
  NNEW(ipev,ITG,nev);
  NNEW(eigxr,double,nev);
  NNEW(aa,double,2*nev);
  NNEW(bb,double,4*nev*nev);
  
  FORTRAN(sortev,(&nev,&nmd,eigxx,&cyclicsymmetry,xx,
		  eigxr,ipev,istartnmd,iendnmd,aa,bb,&nevcomplex));
  
  SFREE(ipev);SFREE(eigxr);SFREE(aa);SFREE(bb);

  /* storing the eigenvalues in the .dat file */

  if(cyclicsymmetry){
    FORTRAN(writeevcscomplex,(eigxx,&nev,nm,&fmin,&fmax));
  }else{
    FORTRAN(writeevcomplex,(eigxx,&nevcomplex,&fmin,&fmax));
  }

  /* storing the participation factors */

  NNEW(eigxr,double,nev);
  NNEW(eigxi,double,nev);
  if(nherm==1){
    NNEW(eig,double,nev);
    for(l=0;l<nev;l++){
      if(d[l]<0.){
	eig[l]=0.;
      }else{
	eig[l]=sqrt(d[l]);
      }
    }
  }else{

    NNEW(eig,double,2*nev);
    for(l=0;l<nev;l++){
      eig[2*l]=sqrt(sqrt(d[2*l]*d[2*l]+d[2*l+1]*d[2*l+1])+d[2*l])/sqrt(2.);
      eig[2*l+1]=sqrt(sqrt(d[2*l]*d[2*l]+d[2*l+1]*d[2*l+1])-d[2*l])/sqrt(2.);
      if(d[2*l+1]<0.) eig[2*l+1]=-eig[2*l+1];
    }
  }
  for(l=0;l<nevcomplex;l++){
    mode=l+1;
    for(k=0;k<nev;k++){
      eigxr[k]=xx[2*l*nev+2*k];
      eigxi[k]=xx[2*l*nev+2*k+1];
    }
    FORTRAN(writepf,(eig,eigxr,eigxi,&zero,&nev,&mode,&nherm));
  }
  SFREE(eigxr);SFREE(eigxi);SFREE(eig);SFREE(d);

  if(cyclicsymmetry){
      
    /* assembling the new eigenmodes */

    /* storage in zz: per eigenmode first the complete real part of
       the eigenvector, then the complete imaginary part */
      
    NNEW(zz,double,(long long)2*nev*neqact);
    for(l=0;l<nev;l++){
      for(i=0;i<neqact;i++){
	for(k=0;k<nev;k++){

	  /* real part */

	  zz[(long long)2*l*neqact+i]+=
	    (xx[2*l*nev+2*k]*z[(long long)2*k*neqact+i]
	     -xx[2*l*nev+2*k+1]*z[(long long)(2*k+1)*neqact+i]);

	  /* imaginary part */

	  zz[(long long)(2*l+1)*neqact+i]+=
	    (xx[2*l*nev+2*k]*z[(long long)(2*k+1)*neqact+i]
	     +xx[2*l*nev+2*k+1]*z[(long long)2*k*neqact+i]);
	}
      }
    }

    /* calculating the scalar product of all old eigenmodes with
       all new eigenmodes => nev x nev matrix */

    NNEW(xmac,double,nev*nev);
    NNEW(xmaccpx,double,4*nev*nev);
    NNEW(bett,double,nev);
    NNEW(betm,double,nev);
    FORTRAN(calcmac,(&neq[1],z,zz,&nev,xmac,xmaccpx,istartnmd,
		     iendnmd,&nmd,&cyclicsymmetry,&neqact,bett,betm,
		     &nevcomplex));
    FORTRAN(writemaccs,(xmac,&nev,nm));

    SFREE(xmac);SFREE(bett);SFREE(betm);SFREE(xmaccpx);
    SFREE(z);
      
    /* normalizing the eigenmodes */
      
    NNEW(z,double,neq[1]);
    for(l=0;l<nev;l++){
      sum=0.;
      DMEMSET(z,0,neq[1],0.);
      FORTRAN(op,(&neq[1],&zz[l*neq[1]],z,adb,aub,jq,irow));
      for(k=0;k<neq[1];k++){
	sum+=zz[l*neq[1]+k]*z[k];
      }
	  
      sum=sqrt(sum);
      for(k=0;k<neq[1];k++){
	zz[l*neq[1]+k]/=sum;
      }
    }
    SFREE(z);

    /* calculating the mass-weighted internal products (eigenvectors are not 
       necessarily orthogonal, since the matrix of the eigenvalue problem is
       not necessarily Hermitian)
       = orthogonality matrices */

    if(mei[3]==1){
	  
      NNEW(xmr,double,nev*nev);
      NNEW(xmi,double,nev*nev);
      NNEW(z,double,neq[1]);
      NNEW(zi,double,neq[1]);
	  
      for(l=0;l<nev;l++){
	DMEMSET(z,0,neq[1],0.);
	FORTRAN(op,(&neq[1],&zz[l*neq[1]],z,adb,aub,jq,irow));
	for(m=l;m<nev;m++){
	  for(k=0;k<neq[1];k++){
	    xmr[l*nev+m]+=zz[m*neq[1]+k]*z[k];
	  }
	}
	      
	memcpy(&zi[0],&zz[(2*l+1)*neqact],sizeof(double)*neqact);
	for(k=0;k<neqact;k++){zi[neqact+k]=-zz[2*l*neqact+k];}
	DMEMSET(z,0,neq[1],0.);
	FORTRAN(op,(&neq[1],zi,z,adb,aub,jq,irow));
	for(m=l;m<nev;m++){
	  for(k=0;k<neq[1];k++){
	    xmi[l*nev+m]+=zz[m*neq[1]+k]*z[k];
	  }
	}
      }
	  
      /* Hermitian part of the matrix */
	  
      for(l=0;l<nev;l++){
	for(m=0;m<l;m++){
	  xmr[l*nev+m]=xmr[m*nev+l];
	  xmi[l*nev+m]=-xmi[m*nev+l];
	}
      }
	  
      for(l=0;l<nev;l++){
	for(m=0;m<nev;m++){
	  printf(" %f",xmr[m*nev+l]);
	}
	printf("\n");
      }
      printf("\n");
      for(l=0;l<nev;l++){
	for(m=0;m<nev;m++){
	  printf(" %f",xmi[m*nev+l]);
	}
	printf("\n");
      }
      SFREE(z);SFREE(zi);
    }

  }else{
      
    /* no cyclic symmmetry */
      
    /* assembling the new eigenmodes */
      
    NNEW(zz,double,2*nevcomplex*neq[1]);
    for(l=0;l<nevcomplex;l++){
      for(j=0;j<2;j++){
	for(i=0;i<neq[1];i++){
	  for(k=0;k<nev;k++){
	    zz[(2*l+j)*neq[1]+i]+=xx[2*l*nev+2*k+j]*
	      z[(long long)k*neq[1]+i];
	  }
	}
      }
    }

    /* calculating the scalar product of all old eigenmodes with
       all new eigenmodes => nev x nevcomplex matrix */

    NNEW(xmac,double,nev*nevcomplex);
    NNEW(xmaccpx,double,4*nev*nevcomplex);
    NNEW(bett,double,nev);
    NNEW(betm,double,nevcomplex);
    FORTRAN(calcmac,(&neq[1],z,zz,&nev,xmac,xmaccpx,istartnmd,
		     iendnmd,&nmd,&cyclicsymmetry,&neqact,bett,betm,
		     &nevcomplex));
    FORTRAN(writemac,(xmac,&nev,&nevcomplex));
    SFREE(xmac);SFREE(bett);SFREE(betm);SFREE(xmaccpx);

    SFREE(z);
      
    /* normalizing the eigenmodes */
      
    NNEW(z,double,neq[1]);
    for(l=0;l<nevcomplex;l++){
      sum=0.;
	  
      /* Ureal^T*M*Ureal */
	  
      DMEMSET(z,0,neq[1],0.);
      FORTRAN(op,(&neq[1],&zz[2*l*neq[1]],z,adb,aub,jq,irow));
      for(k=0;k<neq[1];k++){
	sum+=zz[2*l*neq[1]+k]*z[k];
      }
	  
      /* Uimag^T*M*Uimag */
	  
      DMEMSET(z,0,neq[1],0.);
      FORTRAN(op,(&neq[1],&zz[(2*l+1)*neq[1]],z,adb,aub,jq,irow));
      for(k=0;k<neq[1];k++){
	sum+=zz[(2*l+1)*neq[1]+k]*z[k];
      }
	  
      sum=sqrt(sum);
      for(k=0;k<2*neq[1];k++){
	zz[2*l*neq[1]+k]/=sum;
      }
    }
    SFREE(z);

    /* calculating the mass-weighted internal products (eigenvectors are not 
       necessarily orthogonal, since the matrix of the eigenvalue problem is
       not necessarily symmetric)
       = orthogonality matrices */

    /* mei[3]==1 means that STORAGE=YES was selected on the 
     *COMPLEX FREQUENCY card */

    if(mei[3]==1){
	  
      NNEW(xmr,double,nevcomplex*nevcomplex);
      NNEW(xmi,double,nevcomplex*nevcomplex);
      NNEW(z,double,neq[1]);
	  
      for(l=0;l<nevcomplex;l++){
	sum=0.;
	      
	/* M*Ureal */
	      
	DMEMSET(z,0,neq[1],0.);
	FORTRAN(op,(&neq[1],&zz[2*l*neq[1]],z,adb,aub,jq,irow));
	      
	/* Ureal^T*M*Ureal and Uimag^T*M*Ureal */
	      
	for(m=l;m<nevcomplex;m++){
	  for(k=0;k<neq[1];k++){
	    xmr[l*nevcomplex+m]+=zz[2*m*neq[1]+k]*z[k];
	  }
	  for(k=0;k<neq[1];k++){
	    xmi[l*nevcomplex+m]-=zz[(2*m+1)*neq[1]+k]*z[k];
	  }
	}
	      
	/* M*Uimag */
	      
	DMEMSET(z,0,neq[1],0.);
	FORTRAN(op,(&neq[1],&zz[(2*l+1)*neq[1]],z,adb,aub,jq,irow));
	      
	/* Ureal^T*M*Uimag and Uimag^T*M*Uimag */
	      
	for(m=l;m<nevcomplex;m++){
	  for(k=0;k<neq[1];k++){
	    xmr[l*nevcomplex+m]+=zz[(2*m+1)*neq[1]+k]*z[k];
	  }
	  for(k=0;k<neq[1];k++){
	    xmi[l*nevcomplex+m]+=zz[2*m*neq[1]+k]*z[k];
	  }
	}
      }
	  
      /* Hermitian part of the matrix */
	  
      for(l=0;l<nevcomplex;l++){
	for(m=0;m<l;m++){
	  xmr[l*nevcomplex+m]=xmr[m*nevcomplex+l];
	  xmi[l*nevcomplex+m]=-xmi[m*nevcomplex+l];
	}
      }
	  
      for(l=0;l<nevcomplex;l++){
	for(m=0;m<nevcomplex;m++){
	  printf(" %f",xmr[m*nevcomplex+l]);
	}
	printf("\n");
      }
      printf("\n");
      for(l=0;l<nevcomplex;l++){
	for(m=0;m<nevcomplex;m++){
	  printf(" %f",xmi[m*nevcomplex+l]);
	}
	printf("\n");
      }
      SFREE(z);
    }

  }

  /*storing new eigenmodes and eigenvalues to *.eig-file for later use in 
    steady states dynamic analysis*/

  if(mei[3]==1){

    nherm=0;
      
    if(!cyclicsymmetry){
      if((f1=fopen(fneig,"wb"))==NULL){
	printf(" *ERROR in complexfreq: cannot open eigenvalue file for writing...");
	exit(0);
      }
	  
      /* storing a zero as indication that this was not a
	 cyclic symmetry calculation */
	  
      if(fwrite(&cyclicsymmetry,sizeof(ITG),1,f1)!=1){
	printf(" *ERROR in complexfreq saving the cyclic symmetry flag to the eigenvalue file...");
	exit(0);
      }
	  
      /* not Hermitian */
	  
      if(fwrite(&nherm,sizeof(ITG),1,f1)!=1){
	printf(" *ERROR in complexfreq saving the Hermitian flag to the eigenvalue file...");
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
	  
      if(fwrite(&nevcomplex,sizeof(ITG),1,f1)!=1){
	printf(" *ERROR in complexfreq saving the number of eigenfrequencies to the eigenvalue file...");
	exit(0);
      }
	  
      /* the eigenfrequencies are stored as (radians/time)**2
	 squaring the complexe eigenvalues first */
	  
      NNEW(eigx,double,2*nevcomplex);
      for(i=0;i<nevcomplex;i++){
	eigx[2*i]=eigxx[2*i]*eigxx[2*i]-eigxx[2*i+1]*eigxx[2*i+1];
	eigx[2*i+1]=2.*eigxx[2*i]*eigxx[2*i+1];
      }
	  
      if(fwrite(eigx,sizeof(double),2*nevcomplex,f1)!=2*nevcomplex){
	printf(" *ERROR in complexfreq saving the eigenfrequencies to the eigenvalue file...");
	exit(0);
      }
	  
      SFREE(eigx);
	  
      /* storing the stiffness matrix */
	  
      if(fwrite(ad,sizeof(double),neq[1],f1)!=neq[1]){
	printf(" *ERROR in complexfreq saving the diagonal of the stiffness matrix to the eigenvalue file...");
	exit(0);
      }
      if(fwrite(au,sizeof(double),nzs[2],f1)!=nzs[2]){
	printf(" *ERROR in complexfreq saving the off-diagonal entries of the stiffness matrix to the eigenvalue file...");
	exit(0);
      }
	  
      /* storing the mass matrix */
	  
      if(fwrite(adb,sizeof(double),neq[1],f1)!=neq[1]){
	printf(" *ERROR in complexfreq saving the diagonal of the mass matrix to the eigenvalue file...");
	exit(0);
      }
      if(fwrite(aub,sizeof(double),nzs[1],f1)!=nzs[1]){
	printf(" *ERROR in complexfreq saving the off-diagonal entries of the mass matrix to the eigenvalue file...");
	exit(0);
      }
	  
      /* storing the eigenvectors */
	  
      lfin=0;
      lint=0;
      for(j=0;j<nevcomplex;++j){
	lint=lfin;
	lfin=lfin+2*neq[1];
	if(fwrite(&zz[lint],sizeof(double),2*neq[1],f1)!=2*neq[1]){
	  printf(" *ERROR in complexfreq saving the eigenvectors to the eigenvalue file...");
	  exit(0);
	}
      }
	  
      /* storing the orthogonality matrices */
	  
      if(fwrite(xmr,sizeof(double),nevcomplex*nevcomplex,f1)!=nevcomplex*nevcomplex){
	printf(" *ERROR in complexfreq saving the real orthogonality matrix to the eigenvalue file...");
	exit(0);
      }
	  
      if(fwrite(xmi,sizeof(double),nevcomplex*nevcomplex,f1)!=nevcomplex*nevcomplex){
	printf(" *ERROR in complexfreq saving the imaginary orthogonality matrix to the eigenvalue file...");
	exit(0);
      }
	  
    }else{
	  
      /* opening .eig file for writing */
	  
      if((f1=fopen(fneig,"wb"))==NULL){
	printf(" *ERROR in complexfreq: cannot open eigenvalue file for writing...");
	exit(0);
      }
      /* storing a one as indication that this was a
	 cyclic symmetry calculation */
	  
      if(fwrite(&cyclicsymmetry,sizeof(ITG),1,f1)!=1){
	printf(" *ERROR in complexfreq saving the cyclic symmetry flag to the eigenvalue file...");
	exit(0);
      }
	  
      /* not Hermitian */
	  
      if(fwrite(&nherm,sizeof(ITG),1,f1)!=1){
	printf(" *ERROR in complexfreq saving the Hermitian flag to the eigenvalue file...");
	exit(0);
      }

      /* perturbation parameter iperturb[0] */
	  
      if(fwrite(&iperturb[0],sizeof(ITG),1,f1)!=1){
	printf(" *ERROR saving the perturbation flag to the eigenvalue file...");
	exit(0);
      }
	  
      if(fwrite(&nm[0],sizeof(ITG),1,f1)!=1){
	printf(" *ERROR in complexfreq saving the nodal diameter to the eigqenvalue file...");
	exit(0);
      }
      if(fwrite(&nev,sizeof(ITG),1,f1)!=1){
	printf(" *ERROR in complexfreq saving the number of eigenvalues to the eigenvalue file...");
	exit(0);
      }
	  
      /* the eigenfrequencies are stored as (radians/time)**2
	 squaring the complexe eigenvalues first */
	  
      NNEW(eigx,double,2*nev);
      for(i=0;i<nev;i++){
	eigx[2*i]=eigxx[2*i]*eigxx[2*i]-eigxx[2*i+1]*eigxx[2*i+1];
	eigx[2*i+1]=2.*eigxx[2*i]*eigxx[2*i+1];
      }
	  
      if(fwrite(eigx,sizeof(double),2*nev,f1)!=2*nev){
	printf(" *ERROR in complexfreq saving the eigenfrequencies to the eigenvalue file...");
	exit(0);
      }
	  
      SFREE(eigx);
	  
      /* storing the stiffness matrix */
	  
      if(fwrite(ad,sizeof(double),*neq,f1)!=*neq){
	printf(" *ERROR in complexfreq saving the diagonal of the stiffness matrix to the eigenvalue file...");
	exit(0);
      }
      if(fwrite(au,sizeof(double),*nzs,f1)!=*nzs){
	printf(" *ERROR in complexfreq saving the off-diagonal terms of the stiffness matrix to the eigenvalue file...");
	exit(0);
      }
	  
      /* storing the mass matrix */
	  
      if(fwrite(adb,sizeof(double),*neq,f1)!=*neq){
	printf(" *ERROR in complexfreq saving the diagonal of the mass matrix to the eigenvalue file...");
	exit(0);
      }
      if(fwrite(aub,sizeof(double),*nzs,f1)!=*nzs){
	printf(" *ERROR in complexfreq saving the off-diagonal terms of the mass matrix to the eigenvalue file...");
	exit(0);
      }
	  
      /* storing the eigenvectors */
	  
      lfin=0;
      for(j=0;j<nev;++j){
	lint=lfin;
	lfin=lfin+neq[1];
	if(fwrite(&zz[lint],sizeof(double),neq[1],f1)!=neq[1]){
	  printf(" *ERROR in complexfreq saving the eigenvectors to the eigenvalue file...");
	  exit(0);
	}
      }
	  
      /* storing the orthogonality matrices */
	  
      if(fwrite(xmr,sizeof(double),nev*nev,f1)!=nev*nev){
	printf(" *ERROR in complexfreq saving the real orthogonality matrix to the eigenvalue file...");
	exit(0);
      }
	  
      if(fwrite(xmi,sizeof(double),nev*nev,f1)!=nev*nev){
	printf(" *ERROR in complexfreq saving the imaginary orthogonality matrix to the eigenvalue file...");
	exit(0);
      }
    }
    SFREE(xmr);SFREE(xmi);
      
    fclose(f1);
  }

  SFREE(adb);SFREE(aub);
  if(!cyclicsymmetry){SFREE(ad);SFREE(au);}

  /* calculating the displacements and the stresses and storing */
  /* the results in frd format for each valid eigenmode */

  NNEW(v,double,2*mt**nk);
  NNEW(fn,double,2*mt**nk);
  if((strcmp1(&filab[174],"S   ")==0)||(strcmp1(&filab[1653],"MAXS")==0)|| 
     (strcmp1(&filab[1479],"PHS ")==0)||(strcmp1(&filab[1044],"ZZS ")==0)||
     (strcmp1(&filab[1044],"ERR ")==0)) 
    NNEW(stn,double,12**nk);

  if((strcmp1(&filab[261],"E   ")==0)||(strcmp1(&filab[2523],"MAXE")==0)) 
    NNEW(een,double,12**nk);
  if(strcmp1(&filab[522],"ENER")==0) NNEW(enern,double,2**nk);
  if(strcmp1(&filab[2697],"ME  ")==0) NNEW(emn,double,12**nk);

  NNEW(inum,ITG,*nk);
  NNEW(stx,double,2*6*mi[0]**ne);
  
  NNEW(coefmpcnew,double,*mpcend);

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
  if(strcmp1(&filab[261],"E   ")==0) 
    NNEW(eent,double,2*6**nk*ngraph);
  if((strcmp1(&filab[348],"RF  ")==0)||(strcmp1(&filab[2610],"PRF ")==0))

    // real and imaginary part of the forces

    NNEW(fnt,double,2*mt**nk*ngraph);
  if(strcmp1(&filab[522],"ENER")==0)
    NNEW(enernt,double,*nk*ngraph);
  if((strcmp1(&filab[1044],"ZZS ")==0)||(strcmp1(&filab[1044],"ERR ")==0))
    NNEW(stxt,double,2*6*mi[0]**ne*ngraph);
  if(strcmp1(&filab[2697],"ME  ")==0) 
    NNEW(emnt,double,2*6**nk*ngraph);

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

  if(cyclicsymmetry){

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
  }

  /* check that the tensor fields which are extrapolated from the
     integration points are requested in global coordinates */

  if(strcmp1(&filab[174],"S   ")==0){
    if((strcmp1(&filab[179],"L")==0)&&(*norien>0)){
      printf("\n *WARNING in complexfreq: element fields in cyclic symmetry calculations\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
      strcpy1(&filab[179],"G",1);
    }
  }

  if(strcmp1(&filab[261],"E   ")==0){
    if((strcmp1(&filab[266],"L")==0)&&(*norien>0)){
      printf("\n *WARNING in complexfreq: element fields in cyclic symmetry calculation\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
      strcpy1(&filab[266],"G",1);
    }
  }

  if(strcmp1(&filab[1479],"PHS ")==0){
    if((strcmp1(&filab[1484],"L")==0)&&(*norien>0)){
      printf("\n *WARNING in complexfreq: element fields in cyclic symmetry calculation\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
      strcpy1(&filab[1484],"G",1);
    }
  }

  if(strcmp1(&filab[1653],"MAXS")==0){
    if((strcmp1(&filab[1658],"L")==0)&&(*norien>0)){
      printf("\n *WARNING in complexfreq: element fields in cyclic symmetry calculation\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
      strcpy1(&filab[1658],"G",1);
    }
  }   

  if(strcmp1(&filab[2523],"MAXE")==0){
    if((strcmp1(&filab[2528],"L")==0)&&(*norien>0)){
      printf("\n *WARNING in complexfreq: element fields in cyclic symmetry calculation\n cannot be requested in local orientations;\n the global orientation will be used \n\n");
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

  if(strcmp1(&filab[1566],"MAXU")==0){
    NNEW(vmax,double,4*nkt);
  }

  if(strcmp1(&filab[1653],"MAXS")==0){
    NNEW(stnmax,double,7*nkt);
  }

  if(strcmp1(&filab[2523],"MAXE")==0){
    NNEW(eenmax,double,7*nkt);
  }

  /* storing the results */

  if(!cyclicsymmetry) (neq[1])*=2;

  /* determining a local coordinate system based on the rotation axis */

  if(*nmethod==6){
    FORTRAN(localaxes,(ibody,nbody,xbody,e1,e2,xn));
    NNEW(turdir,char,nev);
  }
  
  lfin=0;
  //  for(j=0;j<nev;++j){
  for(j=0;j<nevcomplex;++j){
    lint=lfin;
    lfin=lfin+neq[1];

    /* calculating the cosine and sine */

    for(k=0;k<*mcs;k++){
      theta=nm[j]*2.*pi/cs[17*k];
      cs[17*k+14]=cos(theta);
      cs[17*k+15]=sin(theta);
    }

    if(*nprint>0)FORTRAN(writehe,(&j));

    NNEW(eei,double,6*mi[0]**ne);
    if(*nener==1){
      NNEW(stiini,double,6*mi[0]**ne);
      NNEW(emeini,double,6*mi[0]**ne);
      NNEW(enerini,double,2*mi[0]**ne);}

    DMEMSET(v,0,2*mt**nk,0.);

    for(k=0;k<neq[1];k+=neq[1]/2){

      for(i=0;i<6*mi[0]**ne;i++){eme[i]=0.;}

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
            else{xreal=zz[lint+idof];ximag=zz[lint+idof+neq[1]/2];}
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
		prestr,iprestr,filab,eme,&emn[kk6],&een[kk6],iperturb,
		f,&fn[kkv],nactdof,&iout,qa,vold,&zz[lint+k],
		nodeboun,ndirboun,xboun,nboun,ipompc,
		nodempc,coefmpcnew,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
		&bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
		ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,&enern[kk],emeini,
		xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
		ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
		nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,
		&ne0,thicke,shcon,nshcon,
		sideload,xload,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		&mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
		islavsurf,ielprop,prop,energyini,energy,&iit,iponoel,
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
		norien,orab,ntmat_,t0,t1old,ithermal,
		prestr,iprestr,filab,eme,&emn[kk6],&een[kk6],iperturb,
		f,&fn[kkv],nactdof,&iout,qa,vold,&zz[lint+k],
		nodeboun,ndirboun,xboun,nboun,ipompc,
		nodempc,coefmpcnew,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
		&bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
		ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,&enern[kk],emeini,
		xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
		ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
		nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,
		&ne0,thicke,shcon,nshcon,
		sideload,xload,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		&mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
		islavsurf,ielprop,prop,energyini,energy,&iit,iponoel,
		inoel,nener,orname,&network,ipobody,xbody,ibody,typeboun,
		itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
		islavelinv,autloc,irowtloc,jqtloc,&nboun2,
		ndirboun2,nodeboun2,xboun2,&nmpc2,ipompc2,nodempc2,coefmpc2,
		labmpc2,ikboun2,ilboun2,ikmpc2,ilmpc2,&mortartrafoflag,
		&intscheme,physcon);
      }

    }
    SFREE(eei);
    if(*nener==1){SFREE(stiini);SFREE(emeini);SFREE(enerini);}

    /* determining the turning direction */

    if(*nmethod==6){
      FORTRAN(turningdirection,(v,e1,e2,xn,mi,nk,&turdir[j],lakon,ipkon,
				kon,ne,co));
    }
   
    /* changing the basic results into cylindrical coordinates
       (only for cyclic symmetry structures */

    for(l=0;l<*nk;l++){inumt[l]=inum[l];}

    if(cyclicsymmetry){
      icntrl=2;imag=1;      
      FORTRAN(rectcyl,(co,v,fn,stn,qfn,een,cs,nk,&icntrl,t,filab,&imag,mi,emn));
    }

    /* copying the basis results (real and imaginary) into
       a new field */

    if((strcmp1(&filab[0],"U  ")==0)||(strcmp1(&filab[870],"PU  ")==0)){
      for(l=0;l<mt**nk;l++){vt[l]=v[l];}
      for(l=0;l<mt**nk;l++){vt[l+mt**nk*ngraph]=v[l+mt**nk];}}
    if(strcmp1(&filab[87],"NT  ")==0)
      for(l=0;l<*nk;l++){t1t[l]=t1[l];};
    if((strcmp1(&filab[174],"S   ")==0)||(strcmp1(&filab[1479],"PHS ")==0)){
      for(l=0;l<6**nk;l++){stnt[l]=stn[l];}
      for(l=0;l<6**nk;l++){stnt[l+6**nk*ngraph]=stn[l+6**nk];}}
    if(strcmp1(&filab[261],"E   ")==0){
      for(l=0;l<6**nk;l++){eent[l]=een[l];};
      for(l=0;l<6**nk;l++){eent[l+6**nk*ngraph]=een[l+6**nk];}}
    if((strcmp1(&filab[348],"RF  ")==0)||(strcmp1(&filab[2610],"PRF ")==0)){
      for(l=0;l<mt**nk;l++){fnt[l]=fn[l];}
      for(l=0;l<mt**nk;l++){fnt[l+mt**nk*ngraph]=fn[l+mt**nk];}}
    if(strcmp1(&filab[522],"ENER")==0)
      for(l=0;l<*nk;l++){enernt[l]=enern[l];};
    if((strcmp1(&filab[1044],"ZZS ")==0)||(strcmp1(&filab[1044],"ERR ")==0)){
      for(l=0;l<6*mi[0]**ne;l++){stxt[l]=stx[l];}
      for(l=0;l<6*mi[0]**ne;l++){stxt[l+6*mi[0]**ne*ngraph]=stx[l+6*mi[0]**ne];}}
    if(strcmp1(&filab[2697],"ME  ")==0){
      for(l=0;l<6**nk;l++){emnt[l]=emn[l];};
      for(l=0;l<6**nk;l++){emnt[l+6**nk*ngraph]=emn[l+6**nk];}}
    
    /* mapping the results to the other sectors
       (only for cyclic symmetric structures */

    if(cyclicsymmetry){

      for(jj=0;jj<*mcs;jj++){
	ilength=cs[17*jj+3];
	is=cs[17*jj+4];
	lprev=cs[17*jj+13];
	for(i=1;i<is;i++){
	  
	  for(l=0;l<*nk;l++){inumt[l+i**nk]=inum[l];}
	  
	  theta=i*nm[j]*2.*pi/cs[17*jj];
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
	    for(l=0;l<*nk;l++){
	      if(inocs[l]==jj) t1t[l+*nk*i]=t1[l];
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
	  
	  if(strcmp1(&filab[522],"ENER")==0){
	    for(l=0;l<*nk;l++){
	      if(inocs[l]==jj) enernt[l+*nk*i]=0.;
	    }
	  }
	}
      }
      
      icntrl=-2;imag=0;
      
      FORTRAN(rectcyl,(cot,vt,fnt,stnt,qfnt,eent,cs,&nkt,&icntrl,t,filab,
		       &imag,mi,emnt));
      
      FORTRAN(rectcylvi,(cot,&vt[mt**nk*ngraph],&fnt[mt**nk*ngraph],
                         &stnt[6**nk*ngraph],qfnt,&eent[6**nk*ngraph],
                         cs,&nkt,&icntrl,t,filab,&imag,mi,&emnt[6**nk*ngraph]));
      
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

    ++*kode;

    /* storing the real part of the eigenfrequencies in freq */

    freq=eigxx[2*j]/6.283185308;
    if(strcmp1(&filab[1044],"ZZS")==0){
      NNEW(neigh,ITG,40*net);
      NNEW(ipneigh,ITG,nkt);
    }
    frd(cot,&nkt,kont,ipkont,lakont,&net,vt,stnt,inumt,nmethod,
	kode,filab,eent,t1t,fnt,&freq,epn,ielmatt,matname,enernt,xstaten,
	nstate_,istep,&iinc,ithermal,qfn,&j,&nm[j],trab,inotrt,
	ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	mi,stxt,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,&net,
	cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emnt,
	thicke,jobnamec,output,qfx,cdn,&mortar,cdnr,cdni,nmat,ielprop,prop,sti);
    if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);}
    
  }   // end loop over the eigenfrequencies
  
  SFREE(xstiff);
  //if(*nbody>0) SFREE(ipobody);
  SFREE(cstr);SFREE(zz);SFREE(eigxx);SFREE(xx);

  if(cyclicsymmetry){
    SFREE(istartnmd);SFREE(iendnmd);
  }else{
    (neq[1])/=2;
  }

  if(*nmethod==6){
  
    FORTRAN(writeturdir,(xn,turdir,&nevcomplex));
    SFREE(turdir);
      
    /* reconverting Coriolis force into centrifugal force */
      
    for(k=0;k<*nbody;k++){
      if(ibody[3*k]==4) ibody[3*k]=1;
    }
  }

  SFREE(nm);SFREE(coefmpcnew);

  if((strcmp1(&filab[174],"S   ")==0)||(strcmp1(&filab[1653],"MAXS")==0)|| 
     (strcmp1(&filab[1479],"PHS ")==0)||(strcmp1(&filab[1044],"ZZS ")==0)||
     (strcmp1(&filab[1044],"ERR ")==0)) 
    SFREE(stn);

  SFREE(v);SFREE(fn);SFREE(inum);SFREE(stx);

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
  if((strcmp1(&filab[1044],"ZZS ")==0)||(strcmp1(&filab[1044],"ERR ")==0)) SFREE(stxt);
  if(strcmp1(&filab[2697],"ME  ")==0) SFREE(emnt);

  SFREE(cot);SFREE(kont);SFREE(ipkont);SFREE(lakont);SFREE(inumt);SFREE(ielmatt);
  if(*ntrans>0){SFREE(inotrt);}

  *ialsetp=ialset;
  *cop=co;*konp=kon;*ipkonp=ipkon;*lakonp=lakon;*ielmatp=ielmat;
  *ielorienp=ielorien;*inotrp=inotr;*nodebounp=nodeboun;
  *ndirbounp=ndirboun;*iambounp=iamboun;*xbounp=xboun;
  *xbounoldp=xbounold;*ikbounp=ikboun;*ilbounp=ilboun;*nactdofp=nactdof;
  *voldp=vold;*emep=eme;*enerp=ener;*ipompcp=ipompc;*nodempcp=nodempc;
  *coefmpcp=coefmpc;*labmpcp=labmpc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
  *fmpcp=fmpc;*veoldp=veold;*iamt1p=iamt1;*t0p=t0;*t1oldp=t1old;*t1p=t1;
  *stip=sti;

  return;
}

