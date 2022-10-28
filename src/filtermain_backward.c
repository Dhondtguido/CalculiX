/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2015 Guido Dhondt                     */

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

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

void filtermain_backward(double *co, double *dgdxglob, ITG *nobject,
                         ITG *nk,ITG *nodedesi, ITG *ndesi, 
		         char *objectset,double *xdesi,double *distmin,
			 ITG *nobjectstart){

  /* filtering the sensitivities */

  ITG *nx=NULL,*ny=NULL,*nz=NULL,i,*ithread=NULL,*mast=NULL,*irow=NULL,
    *icol=NULL,*jq=NULL,*ipointer=NULL,nzs;
    
  double *xo=NULL,*yo=NULL,*zo=NULL,*x=NULL,*y=NULL,*z=NULL,*au=NULL,
    filterrad,*denominator=NULL,*filterval=NULL;

  /* if no radius is defined no filtering is performed
     the radius applies to all objective functions */
    
  if(*nobject==0){return;}
  if(strcmp1(&objectset[81],"     ")==0){
    for(i=2**nk**nobjectstart+1;i<2**nk**nobject;i=i+2){
      dgdxglob[i]=dgdxglob[i-1];
    }
    return;
  }
    
  /* prepare for near3d_se */
    
  NNEW(xo,double,*ndesi);
  NNEW(yo,double,*ndesi);
  NNEW(zo,double,*ndesi);
  NNEW(x,double,*ndesi);
  NNEW(y,double,*ndesi);
  NNEW(z,double,*ndesi);
  NNEW(nx,ITG,*ndesi);
  NNEW(ny,ITG,*ndesi);
  NNEW(nz,ITG,*ndesi);
    
  FORTRAN(prefilter,(co,nodedesi,ndesi,xo,yo,zo,x,y,z,nx,ny,nz,
                     objectset,&filterrad));
    
  /* define structure and the input of the filter matrix */  

  nzs=20000000;
  NNEW(mast,ITG,nzs);
  NNEW(irow,ITG,1);
  NNEW(icol,ITG,*ndesi);
  NNEW(jq,ITG,*ndesi+1);
  NNEW(ipointer,ITG,*ndesi);
    
  mastructfilter(icol,jq,&mast,&irow,ipointer,&nzs,ndesi,nodedesi,
                 xo,yo,zo,x,y,z,nx,ny,nz,objectset,&filterrad);
		  
  SFREE(mast);SFREE(ipointer);    
  RENEW(irow,ITG,nzs);
  NNEW(au,double,nzs);   
  NNEW(denominator,double,*ndesi);
  NNEW(filterval,double,*ndesi);
  
  FORTRAN(filtermatrix,(au,jq,irow,icol,ndesi,nodedesi,&filterrad,co,nk,
                        denominator,objectset,filterval,xdesi,distmin));

			       
  FORTRAN(filter_backward,(au,jq,irow,icol,ndesi,nodedesi,dgdxglob,nobject,
                          nk,nobjectstart,objectset));  
    
  SFREE(denominator);SFREE(filterval);SFREE(au),SFREE(irow);SFREE(jq);
  SFREE(icol);SFREE(xo);SFREE(yo);SFREE(zo);SFREE(x);SFREE(y);
  SFREE(z);SFREE(nx);SFREE(ny);SFREE(nz);
    
  return;
    
} 
