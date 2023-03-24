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

static ITG *nobject1,*nk1,*nodedesi1,*ndesi1,*nx1,*ny1,*nz1,
    num_cpus,ifree1;

/* y1 had to be replaced by yy1, else the following compiler error
   popped up: 

   transitionmain.c:42: error: ‘y1’ redeclared as different kind of symbol */

static double *feasdir1,*xo1,*yo1,*zo1,*x1,*yy1,*z1,*co1;

void transitionmain(double *co, double *feasdir, ITG *nobject, ITG *nk,
                ITG *nodedesi, ITG *ndesi, char *objectset,ITG *ipkon,
		ITG *kon,char *lakon,ITG *ipoface,ITG *nodface,
		ITG *nodedesiinv){

    /* reduction of the sensitivities in the transition from the design
       space to the non-design space */

    ITG *nx=NULL,*ny=NULL,*nz=NULL,ifree,i,*ithread=NULL;
    
    double *xo=NULL,*yo=NULL,*zo=NULL,*x=NULL,*y=NULL,*z=NULL;

    if(*nobject==0){return;}
    if(strcmp1(&objectset[86],"BOU")==0){
    
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
    
       FORTRAN(pretransition,(ipkon,kon,lakon,co,nk,ipoface,nodface,
                           nodedesiinv,xo,yo,zo,x,y,z,nx,ny,nz,&ifree));

       RENEW(xo,double,ifree);
       RENEW(yo,double,ifree);
       RENEW(zo,double,ifree);
       RENEW(x,double,ifree);
       RENEW(y,double,ifree);
       RENEW(z,double,ifree);
       RENEW(nx,ITG,ifree);
       RENEW(ny,ITG,ifree);
       RENEW(nz,ITG,ifree);

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
          of design variables */
    
       if(*ndesi<num_cpus) num_cpus=*ndesi;
    
       pthread_t tid[num_cpus];

       feasdir1=feasdir;nobject1=nobject;nk1=nk;nodedesi1=nodedesi;
       ndesi1=ndesi;objectset1=objectset;xo1=xo;yo1=yo;zo1=zo;
       x1=x;yy1=y;z1=z;nx1=nx;ny1=ny;nz1=nz;ifree1=ifree;co1=co;
       
       /* transition */
    
       printf(" Using up to %" ITGFORMAT " cpu(s) for transition to sensitivities.\n\n", num_cpus);
    
       /* create threads and wait */
  
       NNEW(ithread,ITG,num_cpus);
       for(i=0; i<num_cpus; i++)  {
	  ithread[i]=i;
	  pthread_create(&tid[i], NULL, (void *)transitionmt, (void *)&ithread[i]);
       }
       for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);

       SFREE(xo);SFREE(yo);SFREE(zo);SFREE(ithread);
       SFREE(x);SFREE(y);SFREE(z);SFREE(nx);SFREE(ny);SFREE(nz);
                   
    }
       
    return;
    
} 

/* subroutine for multithreading of transition */

void *transitionmt(ITG *i){

    ITG ndesia,ndesib,ndesidelta;
    
    ndesidelta=(ITG)ceil(*ndesi1/(double)num_cpus);
    ndesia=*i*ndesidelta+1;
    ndesib=(*i+1)*ndesidelta;
    if(ndesib>*ndesi1) ndesib=*ndesi1;

    FORTRAN(transition,(feasdir1,nobject1,nk1,nodedesi1,ndesi1,objectset1,
                        xo1,yo1,zo1,x1,yy1,z1,nx1,ny1,nz1,co1,&ifree1,
                        &ndesia,&ndesib));

    return NULL;
}
