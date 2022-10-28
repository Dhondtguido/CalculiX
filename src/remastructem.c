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

void remastructem(ITG *ipompc, double **coefmpcp, ITG **nodempcp, ITG *nmpc,
              ITG *mpcfree, ITG *nodeboun, ITG *ndirboun, ITG *nboun,
              ITG *ikmpc, ITG *ilmpc, ITG *ikboun, ITG *ilboun,
              char *labmpc, ITG *nk,
              ITG *memmpc_, ITG *icascade, ITG *maxlenmpc,
              ITG *kon, ITG *ipkon, char *lakon, ITG *ne,
              ITG *nactdof, ITG *icol, ITG *jq, ITG **irowp, ITG *isolver,
              ITG *neq, ITG *nzs,ITG *nmethod, double **fp,
              double **fextp, double **bp, double **aux2p, double **finip,
              double **fextinip,double **adbp, double **aubp, ITG *ithermal,
	      ITG *iperturb, ITG *mass, ITG *mi, ITG *ielmat, double *elcon,
	      ITG *ncmat_, ITG *ntmat_,ITG *inomat,ITG *network){

    /* reconstructs the nonzero locations in the stiffness and mass
       matrix after a change in MPC's */

    ITG *nodempc=NULL,*mast1=NULL,*ipointer=NULL,mpcend,
        callfrommain,i,*irow=NULL,mt;

    double *coefmpc=NULL,*f=NULL,*fext=NULL,*b=NULL,*aux2=NULL,
        *fini=NULL,*fextini=NULL,*adb=NULL,*aub=NULL;
    
    nodempc=*nodempcp;coefmpc=*coefmpcp;irow=*irowp;
    f=*fp;fext=*fextp;b=*bp;aux2=*aux2p;fini=*finip;
    fextini=*fextinip;adb=*adbp;aub=*aubp;

    mt=mi[1]+1;

    /* decascading the MPC's */

    printf(" Decascading the MPC's\n\n");
   
    callfrommain=0;
    cascade(ipompc,&coefmpc,&nodempc,nmpc,
	    mpcfree,nodeboun,ndirboun,nboun,ikmpc,
	    ilmpc,ikboun,ilboun,&mpcend,
	    labmpc,nk,memmpc_,icascade,maxlenmpc,
            &callfrommain,iperturb,ithermal);

    /* determining the matrix structure */
    
    printf(" Determining the structure of the matrix:\n");
 
    if(nzs[1]<10) nzs[1]=10;   
    NNEW(mast1,ITG,nzs[1]);
    NNEW(ipointer,ITG,mt**nk);
    RENEW(irow,ITG,nzs[1]);
    
    mastructem(nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,nboun,ipompc,
	     nodempc,nmpc,nactdof,icol,jq,&mast1,&irow,isolver,neq,
	     ikmpc,ilmpc,ipointer,nzs,ithermal,mi,ielmat,elcon,
	     ncmat_,ntmat_,inomat,network);

    SFREE(ipointer);SFREE(mast1);
    RENEW(irow,ITG,nzs[2]);
    
    *nodempcp=nodempc;*coefmpcp=coefmpc;*irowp=irow;

    /* reallocating fields the size of which depends on neq[1] or *nzs */

    RENEW(f,double,neq[1]);for(i=0;i<neq[1];i++) f[i]=0.;
    RENEW(fext,double,neq[1]);for(i=0;i<neq[1];i++) fext[i]=0.;
    RENEW(b,double,neq[1]);for(i=0;i<neq[1];i++) b[i]=0.;
    RENEW(fini,double,neq[1]);for(i=0;i<neq[1];i++) fini[i]=0.;

    if((*nmethod==2)||(*nmethod==4)){
	RENEW(aux2,double,neq[1]);for(i=0;i<neq[1];i++) aux2[i]=0.;
	RENEW(fextini,double,neq[1]);for(i=0;i<neq[1];i++) fextini[i]=0.;
	RENEW(adb,double,neq[1]);for(i=0;i<neq[1];i++) adb[i]=0.;
	RENEW(aub,double,nzs[1]);for(i=0;i<nzs[1];i++) aub[i]=0.;
	mass[0]=1;
    }

    *fp=f;*fextp=fext;*bp=b;*aux2p=aux2;*finip=fini;
    *fextinip=fextini;*adbp=adb;*aubp=aub;

    return;
}


