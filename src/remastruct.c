/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2015 Guido Dhondt                          */

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
#include "CalculiX.h"

void remastruct(ITG *ipompc, double **coefmpcp, ITG **nodempcp, ITG *nmpc,
		ITG *mpcfree, ITG *nodeboun, ITG *ndirboun, ITG *nboun,
		ITG *ikmpc, ITG *ilmpc, ITG *ikboun, ITG *ilboun,
		char *labmpc, ITG *nk,
		ITG *memmpc_, ITG *icascade, ITG *maxlenmpc,
		ITG *kon, ITG *ipkon, char *lakon, ITG *ne,
		ITG *nactdof, ITG *icol, ITG *jq, ITG **irowp, ITG *isolver,
		ITG *neq, ITG *nzs,ITG *nmethod, double **fp,
		double **fextp, double **bp, double **aux2p, double **finip,
		double **fextinip,double **adbp, double **aubp, ITG *ithermal,
		ITG *iperturb, ITG *mass, ITG *mi,ITG *iexpl,ITG *mortar,
		char *typeboun,double **cvp,double **cvinip,ITG *iit,
		ITG *network,ITG *itiefac,ITG *ne0,ITG *nkon0,ITG *nintpoint,
		ITG *islavsurf,double *pmastsurf,char*tieset,ITG *ntie,
		ITG *num_cpus){

  /* reconstructs the nonzero locations in the stiffness and mass
     matrix after a change in MPC's */

  char *lakontot=NULL;
  
  ITG *nodempc=NULL,*mast1=NULL,*ipointer=NULL,mpcend,isiz,
    callfrommain,i,*irow=NULL,mt,im,*ipkontot=NULL,*kontot=NULL,
    netot;

  double *coefmpc=NULL,*f=NULL,*fext=NULL,*b=NULL,*aux2=NULL,
    *fini=NULL,*fextini=NULL,*adb=NULL,*aub=NULL,*cv=NULL,*cvini=NULL;
    
  nodempc=*nodempcp;coefmpc=*coefmpcp;irow=*irowp;
  f=*fp;fext=*fextp;b=*bp;aux2=*aux2p;fini=*finip;
  fextini=*fextinip;adb=*adbp;aub=*aubp;cv=*cvp;cvini=*cvinip;

  mt=mi[1]+1;

  /* decascading the MPC's */

  if(*icascade>0){
    if(*iexpl<=1) printf(" Decascading the MPC's\n\n");
	
    callfrommain=0;
    cascade(ipompc,&coefmpc,&nodempc,nmpc,
	    mpcfree,nodeboun,ndirboun,nboun,ikmpc,
	    ilmpc,ikboun,ilboun,&mpcend,
	    labmpc,nk,memmpc_,icascade,maxlenmpc,
	    &callfrommain,iperturb,ithermal);
  }

  /* determining the matrix structure */
    
  if(*iexpl<=1) printf(" Determining the structure of the matrix:\n\n");
 
  if(nzs[1]<10) nzs[1]=10;   
  NNEW(mast1,ITG,nzs[1]);
  NNEW(ipointer,ITG,mt**nk);

  if(*mortar==1){

    /* for face-to-face penalty: take all possible contact
       elements into account: this number is constant during an
       increment. Therefore, the location of nonzero's does not change 
       and the renumbering during the equation solution (PaStiX) has only 
       to be done once*/
      
    NNEW(ipkontot,ITG,*ne0+*nintpoint);
    NNEW(kontot,ITG,*nkon0+22**nintpoint);
    NNEW(lakontot,char,8*(*ne0+*nintpoint));
    
    isiz=*ne0;cpyparitg(ipkontot,ipkon,&isiz,num_cpus);
    isiz=*nkon0;cpyparitg(kontot,kon,&isiz,num_cpus);
    memcpy(&lakontot[0],&lakon[0],sizeof(char)*8**ne0);
    
    FORTRAN(totalcontact,(tieset,ntie,&netot,ipkontot,kontot,lakontot,
			  islavsurf,itiefac,pmastsurf,ne0,nkon0));

    //   printf(" Maximal possible contact elements = %d\n\n",netot-*ne0);
    
    mastruct(nk,kontot,ipkontot,lakontot,&netot,nodeboun,ndirboun,nboun,
	     ipompc,nodempc,nmpc,nactdof,icol,jq,&mast1,&irow,isolver,neq,
	     ikmpc,ilmpc,ipointer,nzs,nmethod,ithermal,
	     ikboun,ilboun,iperturb,mi,mortar,typeboun,labmpc,
	     iit,icascade,network,iexpl);
    
    SFREE(ipkontot);SFREE(kontot);SFREE(lakontot);
    
  }else{

    /* take actual contact elements into account */
      
    mastruct(nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,nboun,ipompc,
	     nodempc,nmpc,nactdof,icol,jq,&mast1,&irow,isolver,neq,
	     ikmpc,ilmpc,ipointer,nzs,nmethod,ithermal,
	     ikboun,ilboun,iperturb,mi,mortar,typeboun,labmpc,
	     iit,icascade,network,iexpl);
  }

  SFREE(ipointer);SFREE(mast1);
  RENEW(irow,ITG,nzs[2]);
    
  *nodempcp=nodempc;*coefmpcp=coefmpc;*irowp=irow;

  /* reallocating fields the size of which depends on neq[1] or *nzs */

  RENEW(f,double,neq[1]);DMEMSET(f,0,neq[1],0.);
  RENEW(fext,double,neq[1]);DMEMSET(fext,0,neq[1],0.);
  RENEW(b,double,neq[1]);DMEMSET(b,0,neq[1],0.);
  RENEW(fini,double,neq[1]);

  /* for static calculations fini has to be set to f at the
     start of the calculation; in dynamic calculations this is
     not needed, since the initial accelerations has already
     been calculated */

  if((*nmethod!=4)&&(*iit==-1)) DMEMSET(fini,0,neq[1],0.);

  if(*nmethod==4){
    RENEW(aux2,double,neq[1]);DMEMSET(aux2,0,neq[1],0.);
    RENEW(cv,double,neq[1]);
    RENEW(cvini,double,neq[1]);
    RENEW(fextini,double,neq[1]);

    /* the mass matrix is diagonal in an explicit dynamic
       calculation and is not changed by contact; this
       assumes that the number of degrees of freedom does
       not change  */

    if(*iexpl<=1){
      RENEW(adb,double,neq[1]);for(i=0;i<neq[1];i++) adb[i]=0.;
      RENEW(aub,double,nzs[1]);for(i=0;i<nzs[1];i++) aub[i]=0.;
      mass[0]=1;
    }
  }

  *fp=f;*fextp=fext;*bp=b;*aux2p=aux2;*finip=fini;
  *fextinip=fextini;*adbp=adb;*aubp=aub;*cvp=cv;*cvinip=cvini;

  return;
}


