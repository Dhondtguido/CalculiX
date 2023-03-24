/*     CalculiX - A 3-Dimensional finite element program                   */
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

void expand(double *co, ITG *nk, ITG *kon, ITG *ipkon, char *lakon,
	     ITG *ne, ITG *nodeboun, ITG *ndirboun, double *xboun, ITG *nboun, 
	     ITG *ipompc, ITG *nodempc, double *coefmpc, char *labmpc,
             ITG *nmpc, ITG *nodeforc, ITG *ndirforc,double *xforc, 
             ITG *nforc, ITG *nelemload, char *sideload, double *xload,
             ITG *nload, ITG *nactdof, ITG *neq, 
	     ITG *nmethod, ITG *ikmpc, ITG *ilmpc, ITG *ikboun, ITG *ilboun,
	     double *elcon, ITG *nelcon, double *rhcon, ITG *nrhcon,
	     double *alcon, ITG *nalcon, double *alzero, ITG *ielmat,
	     ITG *ielorien, ITG *norien, double *orab, ITG *ntmat_,
	     double *t0,ITG *ithermal,double *prestr, ITG *iprestr, 
	     double *vold,ITG *iperturb, double *sti, ITG *nzs,  
	     double *adb, double *aub,char *filab, double *eme,
             double *plicon, ITG *nplicon, double *plkcon,ITG *nplkcon,
             double *xstate, ITG *npmat_, char *matname, ITG *mi,
	     ITG *ics, double *cs, ITG *mpcend, ITG *ncmat_,
             ITG *nstate_, ITG *mcs, ITG *nkon, double *ener,
             char *jobnamec, char *output, char *set, ITG *nset,ITG *istartset,
             ITG *iendset, ITG *ialset, ITG *nprint, char *prlab,
             char *prset, ITG *nener, double *trab, 
             ITG *inotr, ITG *ntrans, double *ttime, double *fmpc,
	     ITG *nev, double **zp, ITG *iamboun, double *xbounold,
             ITG *nsectors, ITG *nm,ITG *icol,ITG *irow,ITG *nzl, ITG *nam,
             ITG *ipompcold, ITG *nodempcold, double *coefmpcold,
             char *labmpcold, ITG *nmpcold, double *xloadold, ITG *iamload,
             double *t1old,double *t1,ITG *iamt1, double *xstiff,ITG **icolep,
	     ITG **jqep,ITG **irowep,ITG *isolver,
	     ITG *nzse,double **adbep,double **aubep,ITG *iexpl,ITG *ibody,
	     double *xbody,ITG *nbody,double *cocon,ITG *ncocon,
	     char* tieset,ITG* ntie,ITG *imddof,ITG *nmddof,
	     ITG *imdnode,ITG *nmdnode,ITG *imdboun,ITG *nmdboun,
  	     ITG *imdmpc,ITG *nmdmpc, ITG **izdofp, ITG *nzdof,ITG *nherm,
	     double *xmr,double *xmi,char *typeboun,ITG *ielprop,double *prop,
	    char *orname,ITG *itiefac,double *t0g,double *t1g){

  /* calls the Arnoldi Package (ARPACK) for cyclic symmetry calculations */
  
    char *filabt,lakonl[2]=" \0",*labmpc2=NULL;

    ITG *inum=NULL,k,idir,j,iout=0,index,inode,id,i,idof,im,
      ielas,icmd,kk,l,nkt,icntrl,imag=1,icomplex,kkv,kk6,iterm,
      lprev,ilength,ij,i1,i2,iel,ielset,node,indexe,nope,ml1,nelem,
      *inocs=NULL,*ielcs=NULL,jj,l1,l2,is,nlabel,*nshcon=NULL,
      nodeleft,*noderight=NULL,numnodes,ileft,kflag=2,itr,locdir,
      neqh,j1,nodenew,mt=mi[1]+1,istep=1,iinc=1,iit=-1,ne0,
      network=0,noderight_,*izdof=*izdofp,iload,iforc,*iznode=NULL,
      nznode,ll,icfd=0,*inomat=NULL,mortar=0,*islavact=NULL,*ipobody=NULL,
      *islavnode=NULL,*nslavnode=NULL,*islavsurf=NULL,idirnew,
      *iponoel=NULL,*inoel=NULL,mscalmethod=0,intscheme=0,
      *islavelinv=NULL,*irowtloc=NULL,*jqtloc=NULL,nboun2,
      *ndirboun2=NULL,*nodeboun2=NULL,nmpc2,*ipompc2=NULL,*nodempc2=NULL,
      *ikboun2=NULL,*ilboun2=NULL,*ikmpc2=NULL,*ilmpc2=NULL,mortartrafoflag=0;

    long long lint;

    double *stn=NULL,*v=NULL,*temp_array=NULL,*vini=NULL,*csmass=NULL,
      *een=NULL,cam[5],*f=NULL,*fn=NULL,qa[4],*epn=NULL,summass,
      *stiini=NULL,*emn=NULL,*emeini=NULL,*clearini=NULL,
      *xstateini=NULL,theta,pi,*coefmpcnew=NULL,t[3],ctl,stl,
      *stx=NULL,*enern=NULL,*xstaten=NULL,*eei=NULL,*enerini=NULL,
      *qfx=NULL,*qfn=NULL,xreal,ximag,*vt=NULL,sum,*voldt=NULL,
      *coefright=NULL,coef,a[9],ratio,reltime,*physcon=NULL,
      *shcon=NULL,*springarea=NULL,*z=*zp, *zdof=NULL, *thicke=NULL,
      *sumi=NULL,
      *vti=NULL,*pslavsurf=NULL,*pmastsurf=NULL,*cdn=NULL,
      *energyini=NULL,*energy=NULL,*smscale=NULL,
      *autloc=NULL,*xboun2=NULL,*coefmpc2=NULL;
    
    /* dummy arguments for the results call */
    
    double *veold=NULL,*accold=NULL,bet,gam,dtime,time;

    pi=4.*atan(1.);
    neqh=neq[1]/2;

    ne0=*ne;

    noderight_=10;
    NNEW(noderight,ITG,noderight_);
    NNEW(coefright,double,noderight_);

    NNEW(v,double,2*mt**nk);
    NNEW(vt,double,mt**nk**nsectors);
    
    NNEW(fn,double,2*mt**nk);
    NNEW(stn,double,12**nk);
    NNEW(inum,ITG,*nk);
    NNEW(stx,double,6*mi[0]**ne);
    
    nlabel=55;
    NNEW(filabt,char,87*nlabel);
    for(i=1;i<87*nlabel;i++) filabt[i]=' ';
    filabt[0]='U';
    
    NNEW(temp_array,double,neq[1]);
    NNEW(coefmpcnew,double,*mpcend);
    
    nkt=*nsectors**nk;
 
    /* assigning nodes and elements to sectors */
    
    NNEW(inocs,ITG,*nk);
    NNEW(ielcs,ITG,*ne);
    ielset=cs[12];
    if((*mcs!=1)||(ielset!=0)){
	for(i=0;i<*nk;i++) inocs[i]=-1;
	for(i=0;i<*ne;i++) ielcs[i]=-1;
    }
    NNEW(csmass,double,*mcs);
    if(*mcs==1) csmass[0]=1.;
    
    for(i=0;i<*mcs;i++){
	is=cs[17*i];
	//	if(is==1) continue;
	ielset=cs[17*i+12];
	if(ielset==0) continue;
	for(i1=istartset[ielset-1]-1;i1<iendset[ielset-1];i1++){
	    if(ialset[i1]>0){
		iel=ialset[i1]-1;
		if(ipkon[iel]<0) continue;
		ielcs[iel]=i;
		indexe=ipkon[iel];
		if(*mcs==1){
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
		}else{
		  nelem=iel+1;
		  FORTRAN(calcmass,(ipkon,lakon,kon,co,mi,&nelem,ne,thicke,
                        ielmat,&nope,t0,rhcon,nrhcon,ntmat_,
			ithermal,&csmass[i],ielprop,prop));
		}
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
		    if(*mcs==1){
		      if(strcmp1(&lakon[8*iel+3],"2")==0)nope=20;
		      else if (strcmp1(&lakon[8*iel+3],"8")==0)nope=8;
		      else if (strcmp1(&lakon[8*iel+3],"10")==0)nope=10;
		      else if (strcmp1(&lakon[8*iel+3],"4")==0)nope=4;
		      else if (strcmp1(&lakon[8*iel+3],"15")==0)nope=15;
		      else {nope=6;}
		    }else{
		      nelem=iel+1;
		      FORTRAN(calcmass,(ipkon,lakon,kon,co,mi,&nelem,ne,thicke,
                        ielmat,&nope,t0,rhcon,nrhcon,ntmat_,
			ithermal,&csmass[i],ielprop,prop));
		    }
		    for(i2=0;i2<nope;++i2){
			node=kon[indexe+i2]-1;
			inocs[node]=i;
		    }
		}while(1);
	    }
	} 
//	printf("expand.c mass = %" ITGFORMAT ",%e\n",i,csmass[i]);
    }

    /* copying imdnode into iznode 
       iznode contains the nodes in which output is requested and
       the nodes in which loading is applied */

    NNEW(iznode,ITG,*nk);
    for(j=0;j<*nmdnode;j++){iznode[j]=imdnode[j];}
    nznode=*nmdnode;

/* expanding imddof, imdnode, imdboun and imdmpc */

    for(i=1;i<*nsectors;i++){
	for(j=0;j<*nmddof;j++){
	    imddof[i**nmddof+j]=imddof[j]+i*neqh;
	}
	for(j=0;j<*nmdnode;j++){
	    imdnode[i**nmdnode+j]=imdnode[j]+i**nk;
	}
	for(j=0;j<*nmdboun;j++){
	    imdboun[i**nmdboun+j]=imdboun[j]+i**nboun;
	}
	for(j=0;j<*nmdmpc;j++){
	    imdmpc[i**nmdmpc+j]=imdmpc[j]+i**nmpc;
	}
    }
    (*nmddof)*=(*nsectors);
    (*nmdnode)*=(*nsectors);
    (*nmdboun)*=(*nsectors);
    (*nmdmpc)*=(*nsectors);

/* creating a field with the degrees of freedom in which the eigenmodes
   are needed:
   1. all dofs in which the solution is needed (=imddof)
   2. all dofs in which loading was applied
 */	

    NNEW(izdof,ITG,neqh**nsectors);
    for(j=0;j<*nmddof;j++){izdof[j]=imddof[j];}
    *nzdof=*nmddof;
    
    /* generating the coordinates for the other sectors */
    
    icntrl=1;
    
    FORTRAN(rectcyl,(co,v,fn,stn,qfn,een,cs,nk,&icntrl,t,filabt,&imag,mi,emn));
    
    for(jj=0;jj<*mcs;jj++){
	is=(ITG)(cs[17*jj]+0.5);
	for(i=1;i<is;i++){
	    
	    theta=i*2.*pi/cs[17*jj];
	    
	    for(l=0;l<*nk;l++){
		if(inocs[l]==jj){
		    co[3*l+i*3**nk]=co[3*l];
		    co[1+3*l+i*3**nk]=co[1+3*l]+theta;
		    co[2+3*l+i*3**nk]=co[2+3*l];
		    if(*ntrans>0) inotr[2*l+i*2**nk]=inotr[2*l];
		}
	    }
	    for(l=0;l<*nkon;l++){kon[l+i**nkon]=kon[l]+i**nk;}
	    for(l=0;l<*ne;l++){
		if(ielcs[l]==jj){
		    if(ipkon[l]>=0){
			ipkon[l+i**ne]=ipkon[l]+i**nkon;
			ielmat[mi[2]*(l+i**ne)]=ielmat[mi[2]*l];
			if(*norien>0) ielorien[l+i**ne]=ielorien[l];
			for(l1=0;l1<8;l1++){
			    l2=8*l+l1;
			    lakon[l2+i*8**ne]=lakon[l2];
			}
		    }else{
			ipkon[l+i**ne]=ipkon[l];
		    }	
		}
	    }
	}
    }
    
    icntrl=-1;
    
    FORTRAN(rectcyl,(co,vt,fn,stn,qfn,een,cs,&nkt,&icntrl,t,filabt,&imag,mi,emn));

/* expand nactdof */

    for(i=1;i<*nsectors;i++){
	lint=i*mt**nk;
	for(j=0;j<mt**nk;j++){
	    if(nactdof[j]>0){
		nactdof[lint+j]=nactdof[j]+i*neqh;
	    }else{
		nactdof[lint+j]=0;
	    }
	}
    }
    
/* copying the boundary conditions
   (SPC's must be defined in cylindrical coordinates) */
    
    for(i=1;i<*nsectors;i++){
	for(j=0;j<*nboun;j++){
	    nodeboun[i**nboun+j]=nodeboun[j]+i**nk;
	    ndirboun[i**nboun+j]=ndirboun[j];
	    xboun[i**nboun+j]=xboun[j];
	    xbounold[i**nboun+j]=xbounold[j];
	    if(*nam>0) iamboun[i**nboun+j]=iamboun[j];
	    ikboun[i**nboun+j]=ikboun[j]+8*i**nk;
	    ilboun[i**nboun+j]=ilboun[j]+i**nboun;
	}
    }
    
    /* distributed loads */
    
    for(i=0;i<*nload;i++){
	if(nelemload[2*i+1]<*nsectors){
	    nelemload[2*i]+=*ne*nelemload[2*i+1];
	}else{
	    nelemload[2*i]+=*ne*(nelemload[2*i+1]-(*nsectors));
	}
	iload=i+1;
	FORTRAN(addizdofdload,(nelemload,sideload,ipkon,kon,lakon,
		nactdof,izdof,nzdof,mi,&iload,iznode,&nznode,nk,
		imdnode,nmdnode));
    }

    /* body loads */

    if(*nbody>0){
	printf("\n *WARNING in expand: body loads are not allowed for modal dynamics\n          and steady state dynamics calculations in cyclic symmetric structures;\n          the load is removed\n\n");
	*nbody=0;
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

    i=0;
    while(i<*nforc){
	node=nodeforc[2*i];
	
	/* checking for a cylindrical transformation;
	   comparison with the cyclic symmetry system */

            /* carthesian coordinate system */

	    if(nodeforc[2*i+1]<*nsectors){
		nodeforc[2*i]+=*nk*nodeforc[2*i+1];
	    }else{
		nodeforc[2*i]+=*nk*(nodeforc[2*i+1]-(*nsectors));
	    }
	    i++;iforc=i;
	    FORTRAN(addizdofcload,(nodeforc,ndirforc,nactdof,mi,izdof,
		    nzdof,&iforc,iznode,&nznode,nk,imdnode,nmdnode,xforc,
		    ntrans,inotr));
    }
    
    /* loop over all eigenvalues; the loop starts from the highest eigenvalue
       so that the reuse of z is not a problem
       z before: real and imaginary part for a segment for all eigenvalues
       z after: real part for all segments for all eigenvalues */

    if(*nherm==1){
	NNEW(zdof,double,(long long)*nev**nzdof);
    }else{
	NNEW(zdof,double,(long long)2**nev**nzdof);
	NNEW(sumi,double,*nev);
    }

    for(j=*nev-1;j>-1;--j){
	lint=(long long)2*j*neqh;

	/* calculating the cosine and sine of the phase angle */

	for(jj=0;jj<*mcs;jj++){
	    theta=nm[j]*2.*pi/cs[17*jj];
	    cs[17*jj+14]=cos(theta);
	    cs[17*jj+15]=sin(theta);
	}
	
	/* generating the cyclic MPC's (needed for nodal diameters
	   different from 0 */
	
	//	NNEW(eei,double,6*mi[0]**ne);

	DMEMSET(v,0,2*mt**nk,0.);
	
	for(k=0;k<2*neqh;k+=neqh){
	    
	    for(i=0;i<6*mi[0]**ne;i++){eme[i]=0.;}
	    
	    if(k==0) {kk=0;kkv=0;kk6=0;}
	    else {kk=*nk;kkv=mt**nk;kk6=6**nk;}
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
			else{xreal=z[lint+idof];ximag=z[lint+idof+neqh];}
			if(k==0) {
			    if(fabs(xreal)<1.e-30)xreal=1.e-30;
			    coefmpcnew[index]=coefmpc[index]*
				(cs[17*(icomplex-1)+14]+
                                 ximag/xreal*cs[17*(icomplex-1)+15]);}
			else {
			    if(fabs(ximag)<1.e-30)ximag=1.e-30;
			    coefmpcnew[index]=coefmpc[index]*
				(cs[17*(icomplex-1)+14]-
                                 xreal/ximag*cs[17*(icomplex-1)+15]);}
		    }
		    else{coefmpcnew[index]=coefmpc[index];}
		}
	    }
	    
	    results(co,nk,kon,ipkon,lakon,ne,&v[kkv],&stn[kk6],inum,
	      stx,elcon,
	      nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,ielorien,
	      norien,orab,ntmat_,t0,t0,ithermal,
	      prestr,iprestr,filab,eme,&emn[kk6],&een[kk6],iperturb,
	      f,&fn[kkv],nactdof,&iout,qa,vold,&z[lint+k],
	      nodeboun,ndirboun,xboun,nboun,ipompc,
	      nodempc,coefmpcnew,labmpc,nmpc,nmethod,cam,&neqh,veold,accold,
	      &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	      xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
	      ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,&enern[kk],emeini,
	      xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
	      ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	      nelemload,nload,ikmpc,ilmpc,&istep,&iinc,springarea,&reltime,
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
	//	SFREE(eei);

	/* mapping the results to the other sectors */

	if(*nherm!=1)NNEW(vti,double,mt**nk**nsectors);
	
	icntrl=2;imag=1;
	
	FORTRAN(rectcylexp,(co,v,fn,stn,qfn,een,cs,nk,&icntrl,t,filabt,&imag,mi,
			    iznode,&nznode,nsectors,nk,emn));
	
	/* basis sector */

	for(ll=0;ll<nznode;ll++){
	    l1=iznode[ll]-1;
	    for(l2=0;l2<mt;l2++){
		l=mt*l1+l2;
		vt[l]=v[l];
		if(*nherm!=1)vti[l]=v[l+mt**nk];
	    }
	}

	/* other sectors */

	for(jj=0;jj<*mcs;jj++){
	    ilength=cs[17*jj+3];
	    lprev=cs[17*jj+13];
	    for(i=1;i<*nsectors;i++){
		
		theta=i*nm[j]*2.*pi/cs[17*jj];
		ctl=cos(theta);
		stl=sin(theta);
		
		for(ll=0;ll<nznode;ll++){
		    l1=iznode[ll]-1;
		    if(inocs[l1]==jj){
			
			/* check whether node lies on axis */
			
			ml1=-l1-1;
			FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
			if(id!=0){
			    if(ics[lprev+id-1]==ml1){
				for(l2=0;l2<mt;l2++){
				    l=mt*l1+l2;
				    vt[l+mt**nk*i]=v[l];
				    if(*nherm!=1)vti[l+mt**nk*i]=v[l+mt**nk];
				}
				continue;
			    }
			}
			for(l2=0;l2<mt;l2++){
			    l=mt*l1+l2;
			    vt[l+mt**nk*i]=ctl*v[l]-stl*v[l+mt**nk];
			    if(*nherm!=1)vti[l+mt**nk*i]=stl*v[l]+ctl*v[l+mt**nk];
			}
		    }
		}
	    }
	}
	
	icntrl=-2;imag=0;
	
	FORTRAN(rectcylexp,(co,vt,fn,stn,qfn,een,cs,&nkt,&icntrl,t,filabt,
			    &imag,mi,iznode,&nznode,nsectors,nk,emn));
	
/* storing the displacements into the expanded eigenvectors */

	for(ll=0;ll<nznode;ll++){
	    i=iznode[ll]-1;
	    for(j1=0;j1<mt;j1++){

		for(k=0;k<*nsectors;k++){
		    /* C-convention */
		    idof=nactdof[mt*(k**nk+i)+j1]-1;
		    if(idof>-1){
			FORTRAN(nident,(izdof,&idof,nzdof,&id));
			if(id!=0){
			    if(izdof[id-1]==idof){
				if(*nherm==1){
				    zdof[(long long)j**nzdof+id-1]=vt[k*mt**nk+mt*i+j1];
				}else{
				    zdof[(long long)2*j**nzdof+id-1]=vt[k*mt**nk+mt*i+j1];
				    zdof[(long long)(2*j+1)**nzdof+id-1]=vti[k*mt**nk+mt*i+j1];
				}
			    }
			}
		    }
		}
	    }	    
	}

	if(*nherm!=1) SFREE(vti);
	
/* normalizing the eigenvectors with the mass */

	sum=0.;
	summass=0.;
	for(i=0;i<*mcs;i++){
	  if (nm[j]==0||(nm[j]==(ITG)((cs[17*i]/2))&&(fmod(cs[17*i],2.)==0.))){
	    sum+=cs[17*i]*csmass[i];
	  }else{
	    sum+=cs[17*i]*csmass[i]/2.;
	  }
	  summass+=csmass[i];
	}
	if(fabs(summass)>1.e-20){
	  sum=sqrt(sum/summass);
	}else{
	  printf(" *ERROR in expand.c: total mass of structure is zero\n");
	  printf("        maybe no element sets were specified for the\n");
	  printf("        cyclic symmetry ties\n");
	  FORTRAN(stop,());
	}

	if(*nherm==1){
	    for(i=0;i<*nzdof;i++){zdof[(long long)j**nzdof+i]/=sum;}
	}else{
	    for(i=0;i<*nzdof;i++){zdof[(long long)(2*j)**nzdof+i]/=sum;}
	    for(i=0;i<*nzdof;i++){zdof[(long long)(2*j+1)**nzdof+i]/=sum;}
	    sumi[j]=sqrt(sum);
	}
    }

    /* mapping the reference displacements to the other sectors 
       expanding vold */

    if(iperturb[0]==1){
	NNEW(voldt,double,mt**nk**nsectors);
	
	icntrl=2;
	
	FORTRAN(rectcylvold,(co,vold,cs,&icntrl,mi,
			     iznode,&nznode,nsectors,nk));
	
	/* basis sector */
	
	for(ll=0;ll<nznode;ll++){
	    l1=iznode[ll]-1;
	    for(l2=0;l2<mt;l2++){
		l=mt*l1+l2;
		voldt[l]=vold[l];
	    }
	}
	
	/* other sectors */
	
	for(jj=0;jj<*mcs;jj++){
	    ilength=cs[17*jj+3];
	    lprev=cs[17*jj+13];
	    for(i=1;i<*nsectors;i++){
		
		for(ll=0;ll<nznode;ll++){
		    l1=iznode[ll]-1;
		    if(inocs[l1]==jj){
			
			/* check whether node lies on axis */
			
			ml1=-l1-1;
			FORTRAN(nident,(&ics[lprev],&ml1,&ilength,&id));
			if(id!=0){
			    if(ics[lprev+id-1]==ml1){
				for(l2=0;l2<mt;l2++){
				    l=mt*l1+l2;
				    voldt[l+mt**nk*i]=vold[l];
				}
				continue;
			    }
			}
			for(l2=0;l2<mt;l2++){
			    l=mt*l1+l2;
			    voldt[l+mt**nk*i]=vold[l];
			}
		    }
		}
	    }
	}
	
	icntrl=-2;
	
	FORTRAN(rectcylvold,(co,voldt,cs,&icntrl,mi,iznode,&nznode,nsectors,nk));
	
	memcpy(&vold[0],&voldt[0],sizeof(double)*mt*nkt);
	SFREE(voldt);
	
    }

/* copying zdof into z */
  
    if(*nherm==1){
	RENEW(z,double,(long long)*nev**nzdof);
	memcpy(&z[0],&zdof[0],(long long)sizeof(double)**nev**nzdof);
    }else{
	RENEW(z,double,(long long)2**nev**nzdof);
	memcpy(&z[0],&zdof[0],(long long)sizeof(double)*2**nev**nzdof);
	for(i=0;i<*nev;i++){
	    for(j=0;j<*nev;j++){
		xmr[i**nev+j]/=(sumi[i]*sumi[j]);
		xmi[i**nev+j]/=(sumi[i]*sumi[j]);
	    }
	}
	SFREE(sumi);
    }
    SFREE(zdof);

/* copying the multiple point constraints */

    *nmpc=0;
    *mpcend=0;
    for(i=0;i<*nsectors;i++){
	if(i==0){
	    ileft=*nsectors-1;
	}else{
	    ileft=i-1;
	}
	for(j=0;j<*nmpcold;j++){
	    if(noderight_>10){
		noderight_=10;
		RENEW(noderight,ITG,noderight_);
		RENEW(coefright,double,noderight_);
	    }
	    ipompc[*nmpc]=*mpcend+1;
	    ikmpc[*nmpc]=ikmpc[j]+8*i**nk;
	    ilmpc[*nmpc]=ilmpc[j]+i**nmpcold;
	    strcpy1(&labmpc[20**nmpc],&labmpcold[20*j],20);
	    if(strcmp1(&labmpcold[20*j],"CYCLIC")==0){

                /* identifying the nodes on the master side
                   corresponding to a given slave node */

		index=ipompcold[j]-1;
		nodeleft=nodempcold[3*index];

                /* slave node term direction */

		idir=nodempcold[3*index+1];
		index=nodempcold[3*index+2]-1;
		numnodes=0;
		do{
		    node=nodempcold[3*index];
		    if(nodempcold[3*index+1]==idir){
			noderight[numnodes]=node;
			coefright[numnodes]=coefmpcold[index];
			numnodes++;
			if(numnodes>=noderight_){
			    noderight_=(ITG)(1.5*noderight_);
			    RENEW(noderight,ITG,noderight_);
			    RENEW(coefright,double,noderight_);
			}
		    }else{

                        /* master node term direction */

			idirnew=nodempcold[3*index+1];
		    }
		    index=nodempcold[3*index+2]-1;
		    if(index==-1) break;
		}while(1);

                /* if not successful:
                   The cyclic MPC consists of one term corresponding
                   to the slave node and one or more corresponding to
                   the master node; the direction of the first master
                   node term can correspond to the direction of the
                   slave node term, but does not have to. Any master
                   node term direction will actually do */

		if(numnodes==0){

                    /* identifying the nodes on the master side
                       corresponding to a given slave node */

		    index=ipompcold[j]-1;
		    nodeleft=nodempcold[3*index];
		    idir=idirnew;
		    index=nodempcold[3*index+2]-1;
		    numnodes=0;
		    do{
			node=nodempcold[3*index];
			if(nodempcold[3*index+1]==idir){
			    noderight[numnodes]=node;
			    coefright[numnodes]=coefmpcold[index];
			    numnodes++;
			    if(numnodes>=noderight_){
				noderight_=(ITG)(1.5*noderight_);
				RENEW(noderight,ITG,noderight_);
				RENEW(coefright,double,noderight_);
			    }
			}
			index=nodempcold[3*index+2]-1;
			if(index==-1) break;
		    }while(1);

		}

		if(numnodes>0){
		    sum=0.;
		    for(k=0;k<numnodes;k++){
			sum+=coefright[k];
		    }
		    for(k=0;k<numnodes;k++){
			coefright[k]/=sum;
		    }
		}else{coefright[0]=1.;}
		nodempc[3**mpcend]=nodeleft+i**nk;
		nodempc[3**mpcend+1]=idir;
		nodempc[3**mpcend+2]=*mpcend+2;
		coefmpc[*mpcend]=1.;
		for(k=0;k<numnodes;k++){
		    (*mpcend)++;
		    nodempc[3**mpcend]=noderight[k]+ileft**nk;
		    nodempc[3**mpcend+1]=idir;
		    nodempc[3**mpcend+2]=*mpcend+2;
		    coefmpc[*mpcend]=-coefright[k];
		}
		nodempc[3**mpcend+2]=0;
		(*mpcend)++;
	    }else{
		index=ipompcold[j]-1;
		iterm=0;
		do{
		    iterm++;
		    node=nodempcold[3*index];
		    idir=nodempcold[3*index+1];
		    coef=coefmpcold[index];

		    /* check whether SUBCYCLIC MPC: if the current node
                       is an independent node of a CYCLIC MPC, the
                       node in the new MPC should be the cylic previous
                       one */

		    nodenew=node+i**nk;
		    if(strcmp1(&labmpcold[20*j],"SUBCYCLIC")==0){
			for(ij=0;ij<*mcs;ij++){
			    lprev=cs[ij*17+13];
			    ilength=cs[ij*17+3];
			    FORTRAN(nident,(&ics[lprev],&node,&ilength,&id));
			    if(id!=0){
				if(ics[lprev+id-1]==node){
				    nodenew=node+ileft**nk;
				    break;
				}
			    }
			}
		    }
		    
		    /* modification of the MPC coefficients if
                       cylindrical coordinate system is active
                       it is assumed that all terms in the MPC are
                       either in the radial, the circumferential
                       or axial direction  */
		    
		    if(*ntrans<=0){itr=0;}
		    else if(inotr[2*node-2]==0){itr=0;}
		    else{itr=inotr[2*node-2];}
		    
		    if(iterm==1) locdir=-1;
		    
		    if((itr!=0)&&(idir!=0)){
			if(trab[7*itr-1]<0){
			    FORTRAN(transformatrix,(&trab[7*itr-7],
						    &co[3*node-3],a));
			    if(iterm==1){
				for(k=0;k<3;k++){
				    if(fabs(a[3*k+idir-1]-coef)<1.e-10){
					FORTRAN(transformatrix,(&trab[7*itr-7],
								&co[3*nodenew-3],a));
					coef=a[3*k+idir-1];
					locdir=k;
					break;
				    }
				    if(fabs(a[3*k+idir-1]+coef)<1.e-10){
					FORTRAN(transformatrix,(&trab[7*itr-7],
								&co[3*nodenew-3],a));
					coef=-a[3*k+idir-1];
					locdir=k;
					break;
				    }
				}
			    }else{
				if(locdir!=-1){
				    if(fabs(a[3*locdir+idir-1])>1.e-10){
					ratio=coef/a[3*locdir+idir-1];
				    }else{ratio=0.;}
				    FORTRAN(transformatrix,(&trab[7*itr-7],
							    &co[3*nodenew-3],a));
				    coef=ratio*a[3*locdir+idir-1];
				}
			    }		
			}
		    }
		    
		    nodempc[3**mpcend]=nodenew;
		    nodempc[3**mpcend+1]=idir;
		    coefmpc[*mpcend]=coef;
		    index=nodempcold[3*index+2]-1;
		    if(index==-1) break;
		    nodempc[3**mpcend+2]=*mpcend+2;
		    (*mpcend)++;
		}while(1);
		nodempc[3**mpcend+2]=0;
		(*mpcend)++;
	    }
	    (*nmpc)++;
	}
    }

    /* copying the temperatures */

    if(*ithermal!=0){
	for(i=1;i<*nsectors;i++){
	    lint=i**nk;
	    for(j=0;j<*nk;j++){
		t0[lint+j]=t0[j];
		t1old[lint+j]=t1old[j];
		t1[lint+j]=t1[j];
	    }
	}
	if(*nam>0){
	    for(i=1;i<*nsectors;i++){
		lint=i**nk;
		for(j=0;j<*nk;j++){
		    iamt1[lint+j]=iamt1[j];
		}
	    }
	}
    }

    /* copying the contact definition */

//    if(*nmethod==4){
      
      /* first find the startposition to append the expanded contact fields*/
      
/*      for(j=0; j<*nset; j++){
	if(iendset[j]>tint){
	  tint=iendset[j];
	}
      }
      tint++;*/
      /* now append and expand the contact definitons*/
/*      NNEW(tchar1,char,81);
      NNEW(tchar2,char,81);
      NNEW(tchar3,char,81);
      for(i=0; i<*ntie; i++){
	if(tieset[i*(81*3)+80]=='C'){
	  memcpy(tchar2,&tieset[i*(81*3)+81],81);
	  tchar2[80]='\0';
	  memcpy(tchar3,&tieset[i*(81*3)+81+81],81);
	  tchar3[80]='\0';*/
	  //a contact constraint was found, so append and expand the information
/*	  for(j=0; j<*nset; j++){
	    memcpy(tchar1,&set[j*81],81);
	    tchar1[80]='\0';
	    if(strcmp(tchar1,tchar2)==0){*/
	      /* dependent nodal surface was found,copy the original information first */
/*	      tnstart=tint;
	      for(k=0; k<iendset[j]-istartset[j]+1; k++){
		ialset[tint-1]=ialset[istartset[j]-1+k];
		tint++;
	      }*/
	      /* now append the expanded information */
/*	      for(l=1; l<*nsectors; l++){
		for(k=0; k<iendset[j]-istartset[j]+1; k++){
		  ialset[tint-1]=(ialset[istartset[j]-1+k]!=-1)?ialset[istartset[j]-1+k]+*nk*l:-1;
		  tint++;
		}
	      }
	      tnend=tint-1;*/
	      /* now replace the information in istartset and iendset*/
/*	      istartset[j]=tnstart;
	      iendset[j]=tnend;
	    }
	    else if(strcmp(tchar1,tchar3)==0){*/
	      /* independent element face surface was found */
/*	      tnstart=tint;
	      for(k=0; k<iendset[j]-istartset[j]+1; k++){
		ialset[tint-1]=ialset[istartset[j]-1+k];
		tint++;
	      }*/
	      /* now append the expanded information*/
/*	      for(l=1; l<*nsectors; l++){
		for(k=0; k<iendset[j]-istartset[j]+1; k++){
		  tint2=((ITG)(ialset[istartset[j]-1+k]))/10;
		  ialset[tint-1]=(ialset[istartset[j]-1+k]!=-1)?(tint2+*ne*l)*10+(ialset[istartset[j]-1+k]-(tint2*10)):-1;
		  tint++;
		}
	      }
	      tnend=tint-1;*/
	      /* now replace the information in istartset and iendset*/
/*	      istartset[j]=tnstart;
	      iendset[j]=tnend;
	    }
	  }
	}
      }
      SFREE(tchar1);
      SFREE(tchar2);
      SFREE(tchar3);
    }    */
    
    *nk=nkt;
    (*ne)*=(*nsectors);
    (*nkon)*=(*nsectors);
    (*nboun)*=(*nsectors);
    neq[1]=neqh**nsectors;

    *zp=z;*izdofp=izdof;
    
    SFREE(temp_array);SFREE(coefmpcnew);SFREE(noderight);SFREE(coefright);
    SFREE(v);SFREE(vt);SFREE(fn);SFREE(stn);SFREE(inum);SFREE(stx);
    SFREE(inocs);SFREE(ielcs);SFREE(filabt);SFREE(iznode);SFREE(csmass);

    return;
}
