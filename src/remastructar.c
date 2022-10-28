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

void remastructar(ITG *ipompc, double **coefmpcp, ITG **nodempcp, ITG *nmpc,
              ITG *mpcfree, ITG *nodeboun, ITG *ndirboun, ITG *nboun,
              ITG *ikmpc, ITG *ilmpc, ITG *ikboun, ITG *ilboun,
              char *labmpc, ITG *nk,
              ITG *memmpc_, ITG *icascade, ITG *maxlenmpc,
              ITG *kon, ITG *ipkon, char *lakon, ITG *ne,
              ITG *nactdof, ITG *icol, ITG *jq, ITG **irowp, ITG *isolver,
              ITG *neq, ITG *nzs,ITG *nmethod, ITG *ithermal,
	      ITG *iperturb, ITG *mass, ITG *mi, ITG *ics, double *cs,
	      ITG *mcs,ITG *mortar,char *typeboun,ITG *iit,
	      ITG *network,ITG *iexpl){

    /* reconstructs the nonzero locations in the stiffness and mass
       matrix after a change in MPC's or the generation of contact
       spring elements: version for frequency calculations (called
       by arpack, arpackcs and linstatic)  */

    ITG *nodempc=NULL,*mast1=NULL,*ipointer=NULL,mpcend,
        callfrommain,i,*irow=NULL,mt;

    double *coefmpc=NULL;
    
    nodempc=*nodempcp;coefmpc=*coefmpcp;irow=*irowp;

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
    RENEW(irow,ITG,nzs[1]);
  
    if((*mcs==0)||(cs[1]<0)){

	NNEW(ipointer,ITG,mt**nk);
    
	mastruct(nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,nboun,ipompc,
	     nodempc,nmpc,nactdof,icol,jq,&mast1,&irow,isolver,neq,
	     ikmpc,ilmpc,ipointer,nzs,nmethod,ithermal,
	     ikboun,ilboun,iperturb,mi,mortar,typeboun,labmpc,
	     iit,icascade,network,iexpl);

    }else{
      
      NNEW(ipointer,ITG,8**nk);
      
      mastructcs(nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,nboun,
		 ipompc,nodempc,nmpc,nactdof,icol,jq,&mast1,&irow,isolver,
		 neq,ikmpc,ilmpc,ipointer,nzs,nmethod,
		 ics,cs,labmpc,mcs,mi,mortar);
    }

    SFREE(ipointer);SFREE(mast1);
    RENEW(irow,ITG,nzs[2]);
    
    *nodempcp=nodempc;*coefmpcp=coefmpc;*irowp=irow;

    return;
}


