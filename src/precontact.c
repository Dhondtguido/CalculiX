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

#include <stdio.h> 
#include <math.h> 
#include <stdlib.h>
#include <time.h>
#include "CalculiX.h"

void precontact(ITG *ncont, ITG *ntie, char *tieset, ITG *nset, char *set,
        ITG *istartset, ITG *iendset, ITG *ialset, ITG *itietri,
        char *lakon, ITG *ipkon, ITG *kon, ITG *koncont, ITG *ne,
        double *cg, double *straight, double *co,double *vold,
        ITG *istep,ITG *iinc,ITG *iit,ITG *itiefac,
        ITG *islavsurf, ITG *islavnode, ITG *imastnode,
        ITG *nslavnode, ITG *nmastnode,ITG *imastop,ITG *mi,
	ITG *ipe, ITG *ime,double *tietol,ITG *iflagact,
	ITG *nintpoint,double **pslavsurfp,double *xmastnor,double *cs,
	ITG *mcs,ITG *ics,double *clearini,ITG *nslavs){

    /* authors: S. Rakotonanahary, S. Sitzmann and J. Hokkanen */

    ITG i,j,ntrimax,*nx=NULL,*ny=NULL,*nz=NULL,im,
        l,nstart,kflag,ntri,ii;
    
    double *xo=NULL,*yo=NULL,*zo=NULL,*x=NULL,*y=NULL,*z=NULL,
        *pslavsurf=NULL,*clearslavnode=NULL;
    
    pslavsurf=*pslavsurfp;
    
    /* update the location of the center of gravity of 
       the master triangles and the coefficients of their
       bounding planes */
    
    DMEMSET(xmastnor,0,3*nmastnode[*ntie],0.);
    
    FORTRAN(updatecontpen,(koncont,ncont,co,vold,
			   cg,straight,mi,imastnode,nmastnode,xmastnor,
			   ntie,tieset,nset,set,istartset,
			   iendset,ialset,ipkon,lakon,kon,cs,mcs,ics));
    
    /* determining the size of the auxiliary fields 
       (needed for the master triangle search for any
	   given location on the slave faces */	
    
    ntrimax=0;	
    for(i=0;i<*ntie;i++){	    
	if(itietri[2*i+1]-itietri[2*i]+1>ntrimax)		
	    ntrimax=itietri[2*i+1]-itietri[2*i]+1;  	
    }
    
    /* only at the start of a new step */
    
    if ((*istep==1)&&(*iinc==1)&&(*iit<=0)){	    
	NNEW(xo,double,ntrimax);	    
	NNEW(yo,double,ntrimax);	    
	NNEW(zo,double,ntrimax);	    
	NNEW(x,double,ntrimax);	    
	NNEW(y,double,ntrimax);	    
	NNEW(z,double,ntrimax);	   
	NNEW(nx,ITG,ntrimax);	   
	NNEW(ny,ITG,ntrimax);	    
	NNEW(nz,ITG,ntrimax);

	NNEW(clearslavnode,double,3**nslavs);
	    
	FORTRAN(adjustcontactnodes,(tieset,ntie,itietri,cg,straight,
		co,vold,xo,yo,zo,x,y,z,nx,ny,nz,istep,iinc,iit,
		mi,imastop,nslavnode,islavnode,set,nset,istartset,
		iendset,ialset,tietol,clearini,clearslavnode,itiefac,
                ipkon,kon,lakon,islavsurf));
	    
	SFREE(clearslavnode);
	SFREE(xo);SFREE(yo);SFREE(zo);SFREE(x);SFREE(y);SFREE(z);SFREE(nx);	    
	SFREE(ny);SFREE(nz);	    
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
    
    /* Calculating the location of the matched slave/master
       integration points */
    
    RENEW(pslavsurf,double,198);
    
    /* pointer of islavsurf into field pslavsurf and
       pmastsurf */
    
    islavsurf[1]=0;	
    
    /* loop over all ties */
    
    for(i=0;i<*ntie;i++){
	ii=i+1;	   
	
	/* only active contact ties are treated */
	
	if(tieset[i*(81*3)+80]=='C'){		
	    nstart=itietri[2*i]-1;		
	    ntri=itietri[2*i+1]-nstart;		
	    for(j=0;j<ntri;j++){		    
		xo[j]=cg[(nstart+j)*3];		    
		x[j]=xo[j];		   
		nx[j]=j+1;		    
		yo[j]=cg[(nstart+j)*3+1];		    
		y[j]=yo[j];		    
		ny[j]=j+1;		    
		zo[j]=cg[(nstart+j)*3+2];		    
		z[j]=zo[j];		    
		nz[j]=j+1;		
	    }
	    kflag=2;		
	    FORTRAN(dsort,(x,nx,&ntri,&kflag));		
	    FORTRAN(dsort,(y,ny,&ntri,&kflag));		
	    FORTRAN(dsort,(z,nz,&ntri,&kflag));	
	    
	    /* loop over all slave faces belonging to the tie */
	    
	    for(l=itiefac[2*i];l<=itiefac[2*i+1];l++){
		RENEW(pslavsurf,double,3*(*nintpoint+ntri*66));		    
		FORTRAN(slavintpoints,(ntie,itietri,ipkon,kon,
			lakon,straight,nintpoint,koncont,co,vold,
                        xo,yo,zo,x,y,z,nx,ny,nz,islavsurf,
	        	islavnode,nslavnode,imastop,
	        	mi,ncont,ipe,ime,pslavsurf,&ii,&l,&ntri));
	    }	    
	}	
    }
    SFREE(xo);SFREE(yo);SFREE(zo);SFREE(x);SFREE(y);SFREE(z);SFREE(nx);	    
    SFREE(ny);SFREE(nz);
    
    *pslavsurfp=pslavsurf;
    
    return;
}
