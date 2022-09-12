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

void tiedcontact(ITG *ntie, char *tieset, ITG *nset, char *set,
		 ITG *istartset, ITG *iendset, ITG *ialset,
		 char *lakon, ITG *ipkon, ITG *kon,
		 double *tietol,
		 ITG *nmpc, ITG *mpcfree, ITG *memmpc_,
		 ITG **ipompcp, char **labmpcp, ITG **ikmpcp, ITG **ilmpcp,
		 double **fmpcp, ITG **nodempcp, double **coefmpcp,
		 ITG *ithermal, double *co, double *vold, ITG *nef,
		 ITG *nmpc_, ITG *mi, ITG *nk,ITG *istep,ITG *ikboun,
		 ITG *nboun,char *kind1,char *kind2,char *jobnamef){

  char *labmpc=NULL;

  ITG *itietri=NULL,*koncont=NULL,nconf,i,k,*nx=NULL,im,
      *ny=NULL,*nz=NULL,*ifaceslave=NULL,*istartfield=NULL,
      *iendfield=NULL,*ifield=NULL,ntrimax,index,
      ncont,ncone,*ipompc=NULL,*ikmpc=NULL,
      *ilmpc=NULL,*nodempc=NULL,ismallsliding=0,neq,neqterms,
      nmpctied,mortar=0,*ipe=NULL,*ime=NULL,*imastop=NULL,ifreeme;

  double *xo=NULL,*yo=NULL,*zo=NULL,*x=NULL,*y=NULL,*z=NULL,
    *cg=NULL,*straight=NULL,*fmpc=NULL,*coefmpc=NULL;

  ipompc=*ipompcp;labmpc=*labmpcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
  fmpc=*fmpcp;nodempc=*nodempcp;coefmpc=*coefmpcp;

  /* identifying the slave surfaces as nodal or facial surfaces */

  NNEW(ifaceslave,ITG,*ntie);

  FORTRAN(identifytiedface,(tieset,ntie,set,nset,ifaceslave,kind1));

  /* determining the number of triangles of the triangulation
     of the master surface and the number of entities on the
     slave side */

  FORTRAN(allocont,(&ncont,ntie,tieset,nset,set,istartset,iendset,
	  ialset,lakon,&ncone,tietol,&ismallsliding,kind1,
	  kind2,&mortar,istep));

  if(ncont==0){
      SFREE(ifaceslave);return;
  }

  /* allocation of space for the triangulation; 
     koncont(1..3,i): nodes belonging to triangle i
     koncont(4,i): face label to which the triangle belongs =
     10*element+side number */

  NNEW(itietri,ITG,2**ntie);
  NNEW(koncont,ITG,4*ncont);

  /* triangulation of the master surface */

  FORTRAN(triangucont,(&ncont,ntie,tieset,nset,set,istartset,iendset,ialset,
		       itietri,lakon,ipkon,kon,koncont,kind1,kind2,co,nk,
		       &mortar));
  
  /* catalogueing the neighbors of the master triangles */
  
  RENEW(ipe,ITG,*nk);
  RENEW(ime,ITG,12*ncont);
  DMEMSET(ipe,0,*nk,0.);
  DMEMSET(ime,0,12*ncont,0.);
  NNEW(imastop,ITG,3*ncont);

  FORTRAN(trianeighbor,(ipe,ime,imastop,&ncont,koncont,
		        &ifreeme));

  SFREE(ipe);SFREE(ime);

  /* allocation of space for the center of gravity of the triangles
     and the 4 describing planes */

  NNEW(cg,double,3*ncont);
  NNEW(straight,double,16*ncont);
  
  FORTRAN(updatecont,(koncont,&ncont,co,vold,cg,straight,mi));
  
  /* determining the nodes belonging to the slave face surfaces */

  NNEW(istartfield,ITG,*ntie);
  NNEW(iendfield,ITG,*ntie);
  NNEW(ifield,ITG,8*ncone);

  FORTRAN(nodestiedface,(tieset,ntie,ipkon,kon,lakon,set,istartset,
       iendset,ialset,nset,ifaceslave,istartfield,iendfield,ifield,
       &nconf,&ncone,kind1));

  /* determining the maximum number of equations neq */

  if(*nef>0){
    if(ithermal[1]<=1){
      neq=4;
    }else{
      neq=5;
    }
  }else{
    if(ithermal[1]<=1){
      neq=3;
    }else if(ithermal[1]==2){
      neq=1;
    }else{
      neq=4;
    }
  }
  neq*=(ncone+nconf);

  /* reallocating the MPC fields for the new MPC's
     ncone: number of MPC'S due to nodal slave surfaces
     nconf: number of MPC's due to facal slave surfaces */  

  RENEW(ipompc,ITG,*nmpc_+neq);
  RENEW(labmpc,char,20*(*nmpc_+neq)+1);
  RENEW(ikmpc,ITG,*nmpc_+neq);
  RENEW(ilmpc,ITG,*nmpc_+neq);
  RENEW(fmpc,double,*nmpc_+neq);

  /* determining the maximum number of terms;
     expanding nodempc and coefmpc to accommodate
     those terms */
  
  neqterms=9*neq;
  index=*memmpc_;
  (*memmpc_)+=neqterms;
  RENEW(nodempc,ITG,3**memmpc_);
  RENEW(coefmpc,double,*memmpc_);
  for(k=index;k<*memmpc_;k++){
      nodempc[3*k-1]=k+1;
  }
  nodempc[3**memmpc_-1]=0;

  /* determining the size of the auxiliary fields */
  
  ntrimax=0;
  for(i=0;i<*ntie;i++){
    if(itietri[2*i+1]-itietri[2*i]+1>ntrimax)
      ntrimax=itietri[2*i+1]-itietri[2*i]+1;
  }
  NNEW(xo,double,ntrimax);
  NNEW(yo,double,ntrimax);
  NNEW(zo,double,ntrimax);
  NNEW(x,double,ntrimax);
  NNEW(y,double,ntrimax);
  NNEW(z,double,ntrimax);
  NNEW(nx,ITG,ntrimax);
  NNEW(ny,ITG,ntrimax);
  NNEW(nz,ITG,ntrimax);
  
  /* generating the tie MPC's */

  FORTRAN(gentiedmpc,(tieset,ntie,itietri,ipkon,kon,
		      lakon,set,istartset,iendset,ialset,cg,straight,
		      koncont,co,xo,yo,zo,x,y,z,nx,ny,nz,nset,
		      ifaceslave,istartfield,iendfield,ifield,
		      ipompc,nodempc,coefmpc,nmpc,&nmpctied,mpcfree,ikmpc,ilmpc,
		      labmpc,ithermal,tietol,nef,&ncont,imastop,ikboun,nboun,
		      kind1,jobnamef));

  (*nmpc_)+=nmpctied;
  
  SFREE(xo);SFREE(yo);SFREE(zo);SFREE(x);SFREE(y);SFREE(z);SFREE(nx);
  SFREE(ny);SFREE(nz);SFREE(imastop);

  SFREE(ifaceslave);SFREE(istartfield);SFREE(iendfield);SFREE(ifield);
  SFREE(itietri);SFREE(koncont);SFREE(cg);SFREE(straight);

  /* reallocating the MPC fields */

  /*  RENEW(ipompc,ITG,nmpc_);
  RENEW(labmpc,char,20*nmpc_+1);
  RENEW(ikmpc,ITG,nmpc_);
  RENEW(ilmpc,ITG,nmpc_);
  RENEW(fmpc,double,nmpc_);*/

  *ipompcp=ipompc;*labmpcp=labmpc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
  *fmpcp=fmpc;*nodempcp=nodempc;*coefmpcp=coefmpc;

  /*  for(i=0;i<*nmpc;i++){
    j=i+1;
    FORTRAN(writempc,(ipompc,nodempc,coefmpc,labmpc,&j));
    }*/

  return;
}
