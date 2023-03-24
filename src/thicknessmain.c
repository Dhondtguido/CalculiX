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

static char *objectset1;

static ITG *nodedesiboun1,ndesiboun1,*nx1,*ny1,*nz1,num_cpus,ifree1,
    *iobject1,*nk1,*nobject1;

/* y1 had to be replaced by yy1, else the following compiler error
   popped up: 

   thicknessmain.c:42: error: ‘y1’ redeclared as different kind of symbol */

static double *xo1,*yo1,*zo1,*x1,*yy1,*z1,*co1,*dgdxglob1,*coini1,*g01,
    *extnor1;
    

void thicknessmain(double *co,ITG *nobject,ITG *nk,ITG *nodedesi,ITG *ndesi,
                char *objectset,char *set,ITG *nset,ITG *istartset,
		ITG *iendset,ITG *ialset,ITG *iobject,ITG *nodedesiinv,
		double *dgdxglob,double *extnor,double *coini,double *g0){
		
    /* calculation of distance between nodes */

    ITG *nx=NULL,*ny=NULL,*nz=NULL,ifree,i,j,*ithread=NULL,ndesiboun,
        *nodedesiboun=NULL;
    
    double *xo=NULL,*yo=NULL,*zo=NULL,*x=NULL,*y=NULL,*z=NULL;

    /* prepare for near3d */
    
    NNEW(xo,double,*nk);
    NNEW(yo,double,*nk);
    NNEW(zo,double,*nk);
    NNEW(x,double,*nk);
    NNEW(y,double,*nk);
    NNEW(z,double,*nk);
    NNEW(nx,ITG,*nk);
    NNEW(ny,ITG,*nk);
    NNEW(nz,ITG,*nk);
    NNEW(nodedesiboun,ITG,*ndesi); 
    
    FORTRAN(prethickness,(co,xo,yo,zo,x,y,z,nx,ny,nz,&ifree,nodedesiinv,
     			&ndesiboun,nodedesiboun,set,nset,objectset,
     			iobject,istartset,iendset,ialset));
     			 
    RENEW(nodedesiboun,ITG,ndesiboun);
    
    /* variables for multithreading procedure */
    
    ITG sys_cpus;
    char *env,*envloc,*envsys;
    
    num_cpus = 0;
    sys_cpus=0;
    
    /* explicit user declaration prevails */
    
    envsys=getenv("NUMBER_OF_CPUS");
    if(envsys){
       sys_cpus=atoi(envsys);
       if(sys_cpus<0) sys_cpus=0;
    }
    
    /* automatic detection of available number of processors */
    
    if(sys_cpus==0){
     	sys_cpus = getSystemCPUs();
     	if(sys_cpus<1) sys_cpus=1;
    }
    
    /* local declaration prevails, if strictly positive */
    
    envloc = getenv("CCX_NPROC_SENS");
    if(envloc){
       num_cpus=atoi(envloc);
       if(num_cpus<0){
     	  num_cpus=0;
       }else if(num_cpus>sys_cpus){
     	  num_cpus=sys_cpus;
       }

    }
    
    /* else global declaration, if any, applies */
    
    env = getenv("OMP_NUM_THREADS");
    if(num_cpus==0){
       if (env)
     	  num_cpus = atoi(env);
       if (num_cpus < 1) {
     	  num_cpus=1;
       }else if(num_cpus>sys_cpus){
     	 num_cpus=sys_cpus;
       }
    }
    
    /* check that the number of cpus does not supercede the number
       of design variables in nodedesiboun */
    
    if(ndesiboun<num_cpus) num_cpus=ndesiboun;
    
    pthread_t tid[num_cpus];

    NNEW(g01,double,num_cpus**nobject);
    
    nodedesiboun1=nodedesiboun;coini1=coini;nk1=nk;nobject1=nobject;
    ndesiboun1=ndesiboun;objectset1=objectset;xo1=xo;yo1=yo;zo1=zo;
    x1=x;yy1=y;z1=z;nx1=nx;ny1=ny;nz1=nz;ifree1=ifree;co1=co;  
    iobject1=iobject;dgdxglob1=dgdxglob;extnor1=extnor;
    
    /* assessment of actual wallthickness */
    /* create threads and wait */
  
    NNEW(ithread,ITG,num_cpus);
    for(i=0; i<num_cpus; i++)  {
       ithread[i]=i;
       pthread_create(&tid[i], NULL, (void *)thicknessmt, (void *)&ithread[i]);
    }
    for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);

    /* Assembling g0 */
     	 
    g0[*iobject-1]=g01[*iobject-1];
    for(j=1;j<num_cpus;j++){
     g0[*iobject-1]+=g01[*iobject-1+j*(*nobject-1)];
    }
    
    SFREE(xo);SFREE(yo);SFREE(zo);SFREE(g01);
    SFREE(x);SFREE(y);SFREE(z);SFREE(nx);SFREE(ny);SFREE(nz);
    SFREE(ithread);SFREE(nodedesiboun);   
                                     
    return;
    
} 

/* subroutine for multithreading of wallthickness assessment */

void *thicknessmt(ITG *i){

    ITG indexr,indexg0,ndesia,ndesib,ndesidelta;

    indexr=*i*ifree1;
    indexg0=*i**nobject1;
    
    ndesidelta=(ITG)ceil(ndesiboun1/(double)num_cpus);
    ndesia=*i*ndesidelta+1;
    ndesib=(*i+1)*ndesidelta;
    if(ndesib>ndesiboun1) ndesib=ndesiboun1;
    
    //printf("indexr=%" ITGFORMAT","ndesia=%" ITGFORMAT",ndesib=%" ITGFORMAT"\n",indexr,ndesia,ndesib);

    FORTRAN(thickness,(nodedesiboun1,&ndesiboun1,objectset1,xo1,yo1,zo1,
		       x1,yy1,z1,nx1,ny1,nz1,co1,&ifree1,&ndesia,&ndesib,
		       iobject1,dgdxglob1,nk1,extnor1,&g01[indexg0],coini1));
       
    return NULL;
}
