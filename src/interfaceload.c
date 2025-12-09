
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

#include "CalculiX.h"

void interfaceload(ITG *ne,ITG *ipkon,ITG *kon,char *lakon,ITG *nk,char*set,
		   ITG *istartset,ITG *iendset,ITG *ialset,ITG *nset,
		   ITG *nset_,ITG *nalset,double *co,double *vold,ITG *mi,
		   double *cs,ITG *mcs,ITG *ics)
{

  ITG imastset,nmastface,*koncont=NULL,ncont,*ipe=NULL,*ime=NULL,
    *imastop=NULL,*imastnode=NULL,nmasts;

  double *cg=NULL,*straight=NULL,*xmastnor=NULL;

  FORTRAN(calcglobmastsurf,(ne,ipkon,kon,lakon,nk,set,istartset,iendset,
			    ialset,nset,nset_,nalset,&imastset,&nmastface));

  NNEW(koncont,ITG,4*(6*nmastface));

  FORTRAN(triangucont_load,(&ncont,istartset,iendset,ialset,lakon,ipkon,
			    kon,koncont,imastset));

  RENEW(koncont,ITG,4*ncont);

  /* establishing the neighborhood relationships within the
     triangulation */
  
  NNEW(ipe,ITG,*nk);
  NNEW(ime,ITG,12**ncont);
  DMEMSET(ipe,0,*nk,0.);
  DMEMSET(ime,0,12**ncont,0.);
  NNEW(imastop,ITG,3**ncont);

  FORTRAN(trianeighbor,(ipe,ime,imastop,ncont,koncont,
			&ifreeme));

  RENEW(ime,ITG,4*ifreeme);

  NNEW(imastnode,ITG,3**ncont);

  /* catalogue the master nodes */
  
  FORTRAN(catmastnodes,(lakon,ipkon,kon,istartset,iendset,ialset,imastnode,
			&nmasts,&imastset));

  /* calculate geometric data of the master surface triangulation 
     (normals, bounding planes..) */
      
  NNEW(cg,double,3*ncont);
  NNEW(straight,double,16*ncont);
  NNEW(xmastnor,double,3*nmasts);

  FORTRAN(updatecontpen_load(koncont,co,vold,cg,straight,mi,imastnode,
			     &nmasts,xmastnor,istartset,iendset,ialset,
			     ipkon,lakon,kon,cs,mcs,ics,&imastset));
  
  return;
}
