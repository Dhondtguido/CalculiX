/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2007 Guido Dhondt                     */

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

void contact(ITG *ncont, ITG *ntie, char *tieset,ITG *nset,char *set,
	     ITG *istartset, ITG *iendset, ITG *ialset,ITG *itietri,
	     char *lakon, ITG *ipkon, ITG *kon, ITG *koncont, ITG *ne,
	     double *cg, double *straight, ITG *ifree, double *co,
	     double *vold, ITG *ielmat, double *cs, double *elcon,
             ITG *istep,ITG *iinc,ITG *iit,ITG *ncmat_,ITG *ntmat_,
             ITG *ne0,ITG *nmethod,
             ITG *iperturb, ITG *ikboun, ITG *nboun, ITG *mi,
             ITG *imastop,ITG *nslavnode,ITG *islavnode,ITG *islavsurf,
             ITG *itiefac,double *areaslav,ITG *iponoels,ITG *inoels,
             double *springarea, double *tietol, double *reltime,
	     ITG *imastnode, ITG *nmastnode, double *xmastnor,
	     char *filab, ITG *mcs, ITG *ics,
             ITG *nasym,double *xnoels,ITG *mortar,double *pslavsurf,
             double *pmastsurf,double *clearini,double *theta,
             double *xstateini,double *xstate,ITG *nstate_,ITG *icutb,
             ITG *ialeatoric,char *jobnamef,double *alea,double *auw,
	     ITG *jqw,ITG *iroww,ITG *nzsw){

    ITG i,ntrimax,*nx=NULL,*ny=NULL,*nz=NULL,im;
    
    double *xo=NULL,*yo=NULL,*zo=NULL,*x=NULL,*y=NULL,*z=NULL;

    /* next call is only for node-to-face penalty contact
       setting up bordering planes for the master triangles;
       these planes are common between neighboring traingles */

    if(*mortar<=0){

	DMEMSET(xmastnor,0,3*nmastnode[*ntie],0.);
    
	FORTRAN(updatecontpen,(koncont,ncont,co,vold,
			cg,straight,mi,imastnode,nmastnode,xmastnor,
			ntie,tieset,nset,set,istartset,
			iendset,ialset,ipkon,lakon,kon,cs,mcs,ics));
    }
    
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
    
    if(*mortar<=0){
    
	FORTRAN(gencontelem_n2f,(tieset,ntie,itietri,ne,ipkon,kon,lakon,
	  cg,straight,ifree,koncont,
          co,vold,xo,yo,zo,x,y,z,nx,ny,nz,ielmat,elcon,istep,
          iinc,iit,ncmat_,ntmat_,nmethod,mi,
          imastop,nslavnode,islavnode,islavsurf,itiefac,areaslav,iponoels,
          inoels,springarea,
          set,nset,istartset,iendset,ialset,tietol,reltime,
	  filab,nasym,xnoels,icutb,ne0,jobnamef,mortar,auw,jqw,iroww,nzsw,
	  imastnode,nmastnode));

    }else if(*mortar==1){

	FORTRAN(gencontelem_f2f,(tieset,ntie,itietri,ne,ipkon,kon,
	  lakon,cg,straight,ifree,koncont,co,vold,xo,yo,zo,x,y,z,nx,ny,nz,
          ielmat,elcon,istep,iinc,iit,ncmat_,ntmat_,mi,imastop,islavsurf,
	  itiefac,springarea,tietol,reltime,filab,nasym,pslavsurf,pmastsurf,
	  clearini,theta,xstateini,xstate,nstate_,ne0,icutb,ialeatoric,
	  nmethod,jobnamef,alea));

    }

    SFREE(xo);SFREE(yo);SFREE(zo);SFREE(x);SFREE(y);SFREE(z);SFREE(nx);
    SFREE(ny);SFREE(nz);
  
    return;
}
