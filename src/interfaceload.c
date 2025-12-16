
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
		   double *cs,ITG *mcs,ITG *ics,ITG **nelemloadp,
		   char **sideloadp,double **xloadp,double **xloadoldp,
		   ITG **iamloadp,ITG *nam,ITG *nload,ITG *nload_,
		   ITG **imastloadp,double **pmastloadp){

  char *sideload=NULL;
  
  ITG imastset,nmastface,*koncont=NULL,ncont,*ipe=NULL,*ime=NULL,
    *imastop=NULL,*imastnode=NULL,nmasts,*nelemload=NULL,*iamload=NULL,
    *nx=NULL,*ny=NULL,*nz=NULL,kflag,nintpoint,*imastload=NULL;

  double *cg=NULL,*straight=NULL,*xmastnor=NULL,*xload=NULL,*xloadold=NULL,
    *x=NULL,*y=NULL,*z=NULL,*xo=NULL,*yo=NULL,*zo=NULL,*pmastload=NULL;

  nelemload=*nelemloadp;sideload=*sideloadp;xload=*xloadp;xloadold=*xloadoldp;
  iamload=*iamloadp;imastload=*imastloadp;pmastload=*pmastloadp;
  
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

  NNEW(imastnode,ITG,nmasts);

  /* calculate geometric data of the master surface triangulation 
     (normals, bounding planes..) */
      
  NNEW(cg,double,3*ncont);
  NNEW(straight,double,16*ncont);
  NNEW(xmastnor,double,3*nmasts);

  FORTRAN(updatecontpen_load(koncont,co,vold,cg,straight,mi,imastnode,
			     &nmasts,xmastnor,istartset,iendset,ialset,
			     ipkon,lakon,kon,cs,mcs,ics,&imastset));
  
  SFREE(imastnode);SFREE(xmastnor);
 
  NNEW(nelemloadcpy,ITG,2**nload);
  NNEW(sideloadcpy,char,20**nload);
  memcpy(nelemloadcpy,nelemload,sizeof(ITG)*2**nload);
  memcpy(sideloadcpy,sideload,sizeof(char)*20**nload);

  /* preparing the triangles for nearest neighbor search */
  
  NNEW(xo,double,ncont);	
  NNEW(yo,double,ncont);	
  NNEW(zo,double,ncont);	
  NNEW(x,double,ncont);	
  NNEW(y,double,ncont);	
  NNEW(z,double,ncont);	
  NNEW(nx,ITG,ncont);	
  NNEW(ny,ITG,ncont);	
  NNEW(nz,ITG,ncont);
  
  for(j=0;j<ncont;j++){		    
    xo[j]=cg[j*3];		    
    x[j]=xo[j];		   
    nx[j]=j+1;		    
    yo[j]=cg[j*3+1];		    
    y[j]=yo[j];		    
    ny[j]=j+1;		    
    zo[j]=cg[j*3+2];		    
    z[j]=zo[j];		    
    nz[j]=j+1;		
  }
  kflag=2;		
  FORTRAN(dsort,(x,nx,&ncont,&kflag));		
  FORTRAN(dsort,(y,ny,&ncont,&kflag));		
  FORTRAN(dsort,(z,nz,&ncont,&kflag));
  
  SFREE(cg);

  /* check for interface loads (on shell elements) */
  
  RENEW(nelemload,ITG,2*(*nload+nmastface));
  RENEW(sideload,ITG,20*(*nload+nmastface));
  RENEW(xload,double,2*(*nload+nmastface));
  if(*nam>0){
    RENEW(iamload,ITG,2*(*nload+nmastface));
  }

  nintpoint=0;
  for(i=0;i<*nload;i++){
    
    if(strcmp1(&sideloadcpy[20*i],"I")==0){

      RENEW(imastload,ITG,2*(nintpoint+10000));
      RENEW(pmastload,double,3*(nintpoint+10000));
      FORTRAN(mastintpoints,(ipkon,kon,lakon,straight,&nintpoint,koncont,
			     co,vold,xo,yo,zo,x,y,z,nx,ny,nz,imastop,mi,
			     &ncont,ipe,ime,nelemload,sideload,nload,
			     nload_,imastload,pmastload,&nelemloadcpy[2*i],
			     &sideloadcpy[20*i],xload,
			     iamload,nam));
    }
  }
    
  RENEW(imastload,ITG,2*nintpoint);
  RENEW(pmastload,double,3*nintpoint);
  RENEW(nelemload,ITG,2*(*nload));
  RENEW(sideload,ITG,20*(*nload));
  RENEW(xload,double,2*(*nload));
  if(*nam>0){
    RENEW(iamload,ITG,2*nload);
  }

  SFREE(nelemloadcpy);SFREE(sideloadcpy);
  SFREE(x);SFREE(y);SFREE(z);SFREE(xo);SFREE(yo);SFREE(zo);
  SFREE(nx);SFREE(ny);SFREE(nz);
  SFREE(koncont);SFREE(ipe);SFREE(ime);SFREE(imastop);
  SFREE(straight);

  *nelemloadp=nelemload;*sideloadp=sideload;*xloadp=xload;*xloadoldp=xloadold;
  *iamloadp=iamload;*imastloadp=imastload;*pmastloadp=pmastload;
  
  return;
}
