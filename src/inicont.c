/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2007 Guido Dhondt                          */

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

void inicont(ITG * nk,ITG *ncont, ITG *ntie, char *tieset, ITG *nset, char *set,
               ITG *istartset, ITG *iendset, ITG *ialset, ITG **itietrip,
               char *lakon, ITG *ipkon, ITG *kon, ITG **koncontp,
               ITG *nslavs, double *tietol, ITG *ismallsliding, ITG **itiefacp,
               ITG **islavsurfp, ITG **islavnodep, ITG **imastnodep,
               ITG **nslavnodep, ITG **nmastnodep, ITG *mortar,
               ITG **imastopp,ITG *nkon,ITG **iponoelsp,ITG **inoelsp,
	       ITG **ipep, ITG **imep, ITG *ne, ITG *ifacecount,
	       ITG *iperturb, ITG *ikboun, ITG *nboun, double *co,
	       ITG *istep,double **xnoelsp){
    
  char kind1[2]="C",kind2[2]="-";
    
  ITG *itietri=NULL,*koncont=NULL,*itiefac=NULL, *islavsurf=NULL,im,
      *islavnode=NULL,*imastnode=NULL,*nslavnode=NULL,*nmastnode=NULL,
      nmasts,*ipe=NULL,*ime=NULL,*imastop=NULL,
      *iponoels=NULL,*inoels=NULL,ifreenoels,ifreeme,*ipoface=NULL,
      *nodface=NULL,ncone;
    
  double *xnoels=NULL;
    
  itietri=*itietrip;koncont=*koncontp;itiefac=*itiefacp;islavsurf=*islavsurfp;
  islavnode=*islavnodep;imastnode=*imastnodep;nslavnode=*nslavnodep;
  nmastnode=*nmastnodep;imastop=*imastopp,iponoels=*iponoelsp;
  inoels=*inoelsp;ipe=*ipep;ime=*imep;xnoels=*xnoelsp;

  /* determining the number of slave entities (nodes or faces, ncone),
     and the number of master triangles (ncont) */

  FORTRAN(allocont,(ncont,ntie,tieset,nset,set,istartset,iendset,
	  ialset,lakon,&ncone,tietol,ismallsliding,kind1,kind2,mortar,
          istep));
  if(*ncont==0) return;

  NNEW(itietri,ITG,2**ntie);
  NNEW(koncont,ITG,4**ncont);
  
  /* triangulation of the master side */
  
  FORTRAN(triangucont,(ncont,ntie,tieset,nset,set,istartset,iendset,ialset,
		       itietri,lakon,ipkon,kon,koncont,kind1,kind2,co,nk,
		       mortar));

  NNEW(ipe,ITG,*nk);
  NNEW(ime,ITG,12**ncont);
  DMEMSET(ipe,0,*nk,0.);
  DMEMSET(ime,0,12**ncont,0.);
  NNEW(imastop,ITG,3**ncont);

  FORTRAN(trianeighbor,(ipe,ime,imastop,ncont,koncont,
		        &ifreeme));

  if((*mortar==-1)||(*mortar==0)){SFREE(ipe);SFREE(ime);}
  else{RENEW(ime,ITG,4*ifreeme);}

  /* catalogueing the external faces (only for node-to-face
     contact with a nodal slave surface */

  NNEW(ipoface,ITG,*nk);
  NNEW(nodface,ITG,5*6**ne);
  FORTRAN(findsurface,(ipoface,nodface,ne,ipkon,kon,lakon,ntie,
		 tieset));
    
  NNEW(itiefac,ITG,2**ntie);
  RENEW(islavsurf,ITG,2*6**ne);DMEMSET(islavsurf,0,12**ne,0);
  NNEW(islavnode,ITG,8*ncone);
  NNEW(nslavnode,ITG,*ntie+1);
  NNEW(iponoels,ITG,*nk);
  NNEW(inoels,ITG,2**nkon);
  NNEW(xnoels,double,*nkon);
  
  NNEW(imastnode,ITG,3**ncont);
  NNEW(nmastnode,ITG,*ntie+1);
  
  /* catalogueing the slave faces and slave nodes 
     catalogueing the master nodes (only for Mortar contact) */

  FORTRAN(tiefaccont,(lakon,ipkon,kon,ntie,tieset,nset,set,
       istartset,iendset,ialset,itiefac,islavsurf,islavnode,
       imastnode,nslavnode,nmastnode,nslavs,&nmasts,ifacecount,
       iponoels,inoels,&ifreenoels,mortar,ipoface,nodface,nk,
       xnoels));

  RENEW(islavsurf,ITG,2**ifacecount+2);
  RENEW(islavnode,ITG,*nslavs);
  RENEW(inoels,ITG,2*ifreenoels);
  RENEW(xnoels,double,ifreenoels);
  SFREE(ipoface);SFREE(nodface);
  
  RENEW(imastnode,ITG,nmasts);

  *itietrip=itietri;*koncontp=koncont;
  *itiefacp=itiefac;*islavsurfp=islavsurf;
  *islavnodep=islavnode;*imastnodep=imastnode;
  *nslavnodep=nslavnode;*nmastnodep=nmastnode;
  *imastopp=imastop;*iponoelsp=iponoels;*inoelsp=inoels;
  *ipep=ipe;*imep=ime;*xnoelsp=xnoels;
  
  return;
}
