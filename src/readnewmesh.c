
/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                                                                       */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <ctype.h>

#include "CalculiX.h"
#include "readfrd.h"

ITG nkold,mpcrfna,mpcrfnb;

void readnewmesh(char *jobnamec,ITG *nboun,ITG *nodeboun,ITG *iamboun,
		 double *xboun,ITG *nload,char *sideload,ITG *iamload,
		 ITG *nforc,ITG *nodeforc,
		 ITG *iamforc,double *xforc,ITG *ithermal,ITG *nk,
		 double **t1p,ITG **iamt1p,ITG *ne,char **lakonp,ITG **ipkonp,
		 ITG **konp,ITG *istartset,ITG *iendset,ITG *ialset,
		 char *set,ITG *nset,char *filab,double **cop,ITG **ipompcp,
		 ITG **nodempcp,double **coefmpcp,ITG *nmpc,ITG *nmpc_,
		 char **labmpcp,ITG *mpcfree,ITG *memmpc_,ITG **ikmpcp,
		 ITG **ilmpcp,ITG *nk_,ITG *ne_,ITG *nkon_,ITG *istep,
		 ITG *nprop_,ITG **ielpropp,ITG *ne1d,ITG *ne2d,ITG **iponorp,
		 double **thicknp,double **thickep,ITG *mi,double **offsetp,
		 ITG **iponoelp,ITG **rigp,ITG **ne2bounp,ITG **ielorienp,
		 ITG **inotrp,double **t0p,double **t0gp,double **t1gp,
		 double **prestrp,double **voldp,double **veoldp,ITG **ielmatp,
		 ITG *irobustdesign,ITG **irandomtypep,double **randomvalp,
		 ITG *nalset,ITG *nalset_,ITG *nkon,double *xnor,
		 ITG *iaxial,ITG *network,ITG *nlabel,ITG *iuel,ITG *iperturb,
		 ITG *iprestr,ITG *ntie,char *tieset,ITG **iparentelp,
		 ITG *ikboun,ITG *ifreebody,ITG **ipobodyp,ITG *nbody,
		 ITG **iprfnp,ITG **konrfnp,double **ratiorfnp,ITG *nodempcref,
		 double *coefmpcref,ITG *memmpcref_,ITG *mpcfreeref,
		 ITG *maxlenmpcref,ITG *maxlenmpc,ITG *norien,double *tietol)
{

  char masterrfnfile[132]="",fnrfn[132]="",*inpc=NULL,*labmpc=NULL,*lakon=NULL,
    masterurffile[132]="";

  ITG *integerglob=NULL,iglob=1,irefine=1,*inodestet=NULL,nnodestet=0,i,
    istart,j,iquadratic,nenew,nline,nset_=0,*ipoinp=NULL,
    *inp=NULL,*ipoinpc=NULL,idummy[2]={0,0},nuel_=0,inp_size,nentries=18,
    nkold_,neold_,nkonold_,*ielprop=NULL,*iponor=NULL,*iponoel=NULL,
    *rig=NULL,*ne2boun=NULL,*ielorien=NULL,*inotr=NULL,im,mt=mi[1]+1,
    *ielmat=NULL,*irandomtype=NULL,*iparentel=NULL,*ipompc=NULL,*ikmpc=NULL,
    *ilmpc=NULL,*nodempc=NULL,*kon=NULL,*ipkon=NULL,*iamt1=NULL,*ipobody=NULL,
    limit,*ipobody2=NULL,ifreebody2,index,index2,*iprfn=NULL,*konrfn=NULL,
    *jq=NULL,*irow=NULL,*icol=NULL,*loc=NULL,*irowt=NULL,*jqt=NULL,
    *itemp=NULL,*ixcol=NULL,*ipoface=NULL,*nodface=NULL;

  double *doubleglob=NULL,sigma=0.,*thickn=NULL,*thicke=NULL,*offset=NULL,
    *t0=NULL,*t0g=NULL,*t1g=NULL,*prestr=NULL,*vold=NULL,*veold=NULL,
    *randomval=NULL,*t1=NULL,*coefmpc=NULL,*co=NULL,*ratiorfn=NULL,
    *au=NULL;

  ielprop=*ielpropp;iponor=*iponorp;thickn=*thicknp;thicke=*thickep;
  offset=*offsetp;iponoel=*iponoelp;rig=*rigp;ne2boun=*ne2bounp;
  ielorien=*ielorienp;inotr=*inotrp;t0=*t0p;t0g=*t0gp;t1g=*t1gp;
  prestr=*prestrp;vold=*voldp;veold=*veoldp;ielmat=*ielmatp;
  irandomtype=*irandomtypep;randomval=*randomvalp;t1=*t1p;
  iparentel=*iparentelp;ipompc=*ipompcp;labmpc=*labmpcp;ikmpc=*ikmpcp;
  ilmpc=*ilmpcp;nodempc=*nodempcp;coefmpc=*coefmpcp;co=*cop;kon=*konp;
  ipkon=*ipkonp;lakon=*lakonp;iamt1=*iamt1p;ipobody=*ipobodyp;
  iprfn=*iprfnp;konrfn=*konrfnp;ratiorfn=*ratiorfnp;
  
  /* adding the refined mesh (only in the first step) */
  
  if(*istep==1){

    nkold=*nk;
  
    /* reading the input file with information on the refined mesh */

    NNEW(ipoinp,ITG,2*nentries);

    strcpy2(fnrfn,jobnamec,132);
    strcat(fnrfn,".rfn");

    readinput(fnrfn,&inpc,&nline,&nset_,ipoinp,&inp,&ipoinpc,idummy,&nuel_,
	      &inp_size); 

    /* determine new values for nk_, ne_ and nkon_ */

    nkold_=*nk_;
    neold_=*ne_;
    nkonold_=*nkon_;
    FORTRAN(allocation_rfn,(nk_,ne_,nkon_,ipoinp,ipoinpc,inpc,inp));

    /* reallocating the fields depending on nk_, ne_ and nkon_ */

    RENEW(co,double,3**nk_);
    RENEW(kon,ITG,*nkon_);
    RENEW(ipkon,ITG,*ne_);
    for(i=neold_;i<*ne_;i++) ipkon[i]=-1;
    RENEW(lakon,char,8**ne_);
    if(*nprop_>0){
      RENEW(ielprop,ITG,*ne_);
      for(i=neold_;i<*ne_;i++) ielprop[i]=-1;
    }
    if((*ne1d!=0)||(*ne2d!=0)){
      RENEW(iponor,ITG,2**nkon_);
      for(i=2*nkonold_;i<2**nkon_;i++) iponor[i]=-1;
      RENEW(thickn,double,2**nk_);
      RENEW(thicke,double,mi[2]**nkon_);
      RENEW(offset,double,2**ne_);
      RENEW(iponoel,ITG,*nk_);
      RENEW(rig,ITG,*nk_);
      RENEW(ne2boun,ITG,2**nk_);
    }
    RENEW(ielorien,ITG,mi[2]**ne_);
    RENEW(inotr,ITG,2**nk_);
    RENEW(t0,double,*nk_);
    RENEW(t1,double,*nk_);
    if((*ne1d!=0)||(*ne2d!=0)){
      RENEW(t0g,double,2**nk_);
      RENEW(t1g,double,2**nk_);
    }
    DMEMSET(t0,nkold_,*nk_,1.2357111319);
    DMEMSET(t1,nkold_,*nk_,1.2357111319);
    RENEW(iamt1,ITG,*nk_);
    RENEW(prestr,double,6*mi[0]**ne_);
    RENEW(vold,double,mt**nk_);
    RENEW(veold,double,mt**nk_);
    RENEW(ielmat,ITG,mi[2]**ne_);
    if(irobustdesign[0]>0){
      RENEW(irandomtype,ITG,*nk_);
      RENEW(randomval,double,2**nk_);
    }

    NNEW(iparentel,ITG,*ne_);

    FORTRAN(calinput_rfn,(co,filab,set,istartset,iendset,ialset,nset,&nset_,
			  nalset,nalset_,mi,kon,ipkon,lakon,nkon,ne,ne_,
			  iponor,xnor,istep,ipoinp,inp,iaxial,ipoinpc,
			  network,nlabel,iuel,&nuel_,ielmat,inpc,iperturb,
			  iprestr,nk,nk_,ntie,tieset,iparentel,tietol));

    RENEW(iparentel,ITG,*ne);

    /* transferring the material and orientation information from the
       parent elements */
    
    for(i=0;i<*ne_;i++){
      if(iparentel[i]>0){
	ielmat[i]=ielmat[iparentel[i]-1];
	//	ielorien[i]=ielorien[iparentel[i]-1];
      }
    }

    if(*norien>0){
      for(i=0;i<*ne_;i++){
	if(iparentel[i]>0){
	  //	  ielmat[i]=ielmat[iparentel[i]-1];
	  ielorien[i]=ielorien[iparentel[i]-1];
	}
      }
    }
  
    /* get the nodes and topology of the refined mesh of the part of the 
       mesh which was refined */
  
    strcpy2(masterrfnfile,jobnamec,132);
    strcat(masterrfnfile,".rfn.frd");
  
    getglobalresults(masterrfnfile,&integerglob,&doubleglob,nboun,iamboun,xboun,
		     nload,sideload,iamload,&iglob,nforc,iamforc,xforc,
		     ithermal,nk,t1,iamt1,&sigma,&irefine);

    /* determine all nodes at the free surface of the part of the old
       mesh which is refined */
  
    NNEW(inodestet,ITG,*nk);
    NNEW(ipoface,ITG,*nk);
    NNEW(nodface,ITG,20**ne);

    FORTRAN(getnodesinitetmesh,(ne,lakon,ipkon,kon,istartset,iendset,ialset,set,
				nset,filab,inodestet,&nnodestet,nodface,
				ipoface,nk));

    SFREE(ipoface);SFREE(nodface);

    RENEW(inodestet,ITG,nnodestet);

    //  for(i=0;i<nnodestet;i++) {printf("%d\n",inodestet[i]);}

    /* create MPC's for nodes in the old mesh in which SPC's or
       point forces were defined */
  
    *nmpc_+=3*nnodestet;
    RENEW(ipompc,ITG,*nmpc_);
    RENEW(labmpc,char,20**nmpc_+1);
    RENEW(ikmpc,ITG,*nmpc_);
    RENEW(ilmpc,ITG,*nmpc_);
  
    istart=*memmpc_+1;
    nodempc[3**memmpc_-1]=istart;
  
    nenew=integerglob[1];

    iquadratic=0;
    for(i=0;i<nenew;i++){
      if(ipkon[i]>=0){
	if(strcmp1(&lakon[8*i],"C3D10   ")==0){
	  iquadratic=1;
	  break;
	}
      }
    }
    if(iquadratic==0) {
      *memmpc_+=3*5*nnodestet;
    }else {
      *memmpc_+=3*11*nnodestet;
    }

    RENEW(nodempc,ITG,3**memmpc_);
    RENEW(coefmpc,double,*memmpc_);
    for(j=istart;j<*memmpc_;j++){
      nodempc[3*j-1]=j+1;
    }
    nodempc[3**memmpc_-1]=0;

    mpcrfna=*nmpc+1;

    FORTRAN(genmpc,(inodestet,&nnodestet,co,doubleglob,integerglob,ipompc,
		    nodempc,coefmpc,nmpc,nmpc_,labmpc,mpcfree,ikmpc,ilmpc));

    SFREE(inodestet);mpcrfnb=*nmpc;
    
    /*   for(i=0;i<*nmpc;i++){
      j=i+1;
      FORTRAN(writempc,(ipompc,nodempc,coefmpc,labmpc,&j));
      }*/

    SFREE(integerglob);SFREE(doubleglob);

    if(ithermal[0]>0){
    
      /* create MPC's to determine the temperatures t0 and t1 in
	 the new mesh based on the values in the old mesh */

      strcpy2(masterurffile,jobnamec,132);
      strcat(masterurffile,".urf.frd");

      /* reading the old mesh */

      getglobalresults(masterurffile,&integerglob,&doubleglob,nboun,iamboun,
		       xboun,nload,sideload,iamload,&iglob,nforc,iamforc,xforc,
		       ithermal,nk,t1,iamt1,&sigma,&irefine);

      NNEW(iprfn,ITG,*nk-nkold+1);
      NNEW(konrfn,ITG,20*(*nk-nkold));
      NNEW(ratiorfn,double,20*(*nk-nkold));

      FORTRAN(genratio,(co,doubleglob,integerglob,&nkold,nk,iprfn,konrfn,
			ratiorfn));

      RENEW(konrfn,ITG,iprfn[*nk-nkold]);
      RENEW(ratiorfn,double,iprfn[*nk-nkold]);

      /* interpolating the initial temperatures t0 */

      for(i=0;i<*nk-nkold;i++){
	t0[nkold+i]=0.;
	for(j=0;j<iprfn[i+1]-iprfn[i];j++){
	  t0[nkold+i]+=ratiorfn[iprfn[i]+j]*t0[konrfn[iprfn[i]+j]-1];
	}
      }

    }
    
  }

  /* the remaining lines are executed in each step */

  /* rearranging the equations connecting the nodes at the surface
     of the unrefined mesh to the nodes at the surface of the refined
     mesh in case some of the former are subject to SPC's or MPC's */

  nnodestet=mpcrfnb-mpcrfna+1;
  NNEW(jq,ITG,nnodestet+1);
  NNEW(irow,ITG,10*nnodestet);
  NNEW(au,double,10*nnodestet);
  NNEW(icol,ITG,nnodestet);
  NNEW(ixcol,ITG,nnodestet);
  NNEW(loc,ITG,10*nnodestet);
  NNEW(itemp,ITG,10*nnodestet);
  NNEW(irowt,ITG,10*nnodestet);
  NNEW(jqt,ITG,*nk+1);
  NNEW(inodestet,ITG,nnodestet);

  FORTRAN(modifympc,(inodestet,&nnodestet,co,doubleglob,integerglob,ipompc,
		     nodempc,coefmpc,nmpc,nmpc_,labmpc,mpcfree,ikmpc,ilmpc,
		     jq,irow,icol,loc,irowt,jqt,itemp,au,ixcol,ikboun,
		     nboun,nodeboun,&mpcrfna,&mpcrfnb,nodempcref,coefmpcref,
		     memmpcref_,mpcfreeref,maxlenmpcref,memmpc_,maxlenmpc,
		     istep));

  SFREE(jq);SFREE(irow);SFREE(icol);SFREE(loc);SFREE(irowt);SFREE(jqt);
  SFREE(itemp);SFREE(au);SFREE(ixcol);SFREE(inodestet);
    
  /*for(i=0;i<*nmpc;i++){
    j=i+1;
    FORTRAN(writempc,(ipompc,nodempc,coefmpc,labmpc,&j));
    }*/

  /* transferring the body load information from the
     parent elements */

   if(*nbody>0){
    limit=(ITG)(1.1**ne);
    if(limit<100) limit=100;
    NNEW(ipobody2,ITG,2*limit);
    ifreebody2=*ne+1;
  
    for(i=0;i<*ne;i++){
      if(ipkon[i]<0) continue;
      
      if(iparentel[i]==0){
	index=i+1;
      }else{
	index=iparentel[i];
      }
      
      index2=i+1;
      ipobody2[2*index2-2]=ipobody[2*index-2];
      index=ipobody[2*index-1];

      do{
	if(index!=0){
	  ipobody2[2*index2-1]=ifreebody2;
	  index2=ifreebody2;
	  
	  ifreebody2++;
	  if(ifreebody2>limit){
	    limit=(ITG)(1.1*limit);
	    RENEW(ipobody2,ITG,2*limit);
	  }

	  ipobody2[2*index2-2]=ipobody[2*index-2];
	  index=ipobody[2*index-1];
	}else{
	  ipobody2[2*index2-1]=0;
	  break;
	}
      }while(1);
    }
  
    RENEW(ipobody,ITG,2*(ifreebody2-1));
    memcpy(ipobody,ipobody2,sizeof(ITG)*2*(ifreebody2-1));
    *ifreebody=ifreebody2;
    SFREE(ipobody2);
    }

  /* remove the refine request, if any */
  
  if(strcmp1(&filab[4089],"RM")==0){
    strcpy1(&filab[4089],"  ",2);
  }

  if(ithermal[0]>0){

    /* interpolating the initial temperatures t1 */

    for(i=0;i<*nk-nkold;i++){
      t1[nkold+i]=0.;
      for(j=0;j<iprfn[i+1]-iprfn[i];j++){
	t1[nkold+i]+=ratiorfn[iprfn[i]+j]*t1[konrfn[iprfn[i]+j]-1];
      }
    }

  }

  *ielpropp=ielprop;*iponorp=iponor;*thicknp=thickn;*thickep=thicke;
  *offsetp=offset;*iponoelp=iponoel;*rigp=rig;*ne2bounp=ne2boun;
  *ielorienp=ielorien;*inotrp=inotr;*t0p=t0;*t0gp=t0g;*t1gp=t1g;
  *prestrp=prestr;*voldp=vold;*veoldp=veold;*ielmatp=ielmat;
  *irandomtypep=irandomtype;*randomvalp=randomval;*t1p=t1;
  *iparentelp=iparentel;*ipompcp=ipompc;*labmpcp=labmpc;*ikmpcp=ikmpc;
  *ilmpcp=ilmpc;*nodempcp=nodempc;*coefmpcp=coefmpc;*cop=co;*konp=kon;
  *ipkonp=ipkon;*lakonp=lakon;*iamt1p=iamt1;*ipobodyp=ipobody;
  *iprfnp=iprfn;*konrfnp=konrfn;*ratiorfnp=ratiorfn;
  
  return;

}
