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

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

static ITG num_cpus,*nx1,*ny1,*nz1,*ifatet1,*netet1,*kontet1,*konl1,*nkf1,
  *iparent1,*ikf1;

/* y1 had to be replaced by yy1, else the following compiler error
   popped up: 

interpoltetmain.c:28:20: error: ‘y1’ redeclared as different kind of symbol
 static double *x1,*y1,*z1,*xo1,*yo1,*zo1,*planfa1,*cotet1,*ratio1,*co1;
                    ^~
In file included from /usr/include/features.h:423:0,
                 from /usr/include/unistd.h:25,
                 from interpoltetmain.c:18:
/usr/include/bits/mathcalls.h:227:1: note: previous declaration of ‘y1’ was here
__MATHCALL (y1,, (_Mdouble_)); */

static double *x1,*yy1,*z1,*xo1,*yo1,*zo1,*planfa1,*cotet1,*ratio1,*co1;

void interpoltetmain(double *planfa,ITG *ifatet,ITG *netet_,ITG *kontet,
		     double *cotet,ITG *iparent,double *bc,double *co,
		     ITG *nkf,ITG *ikf,ITG *konl,double *ratio){

  /* the master mesh is described by:
     netet_: upper limit of tetrahedral numbers
     kontet(4,*): topology of the tetrahedral elements
     cotet(3,*): coordinates of the nodes of the tets = 
     coordinates of the fluid element centers
     ifatet(4,*): faces of the tets
     planfa(16,*): coefficients of the facial planes
     bc(4,*): coordinates and radius of the inscribed sphere
     iparent(*): parent elements of the tets numbered successively
     the parent number is used in kontet, ifatet 
     and bc

     nkf: number of fluid nodes numbered successively
     ikf(i): node number of the fluid nodes according to the input
     deck (used in co)
     co(3,*): coordinates of the fluid nodes */
     
  
  ITG i,netet,isiz,kflag,*nx=NULL,*ny=NULL,*nz=NULL;

  double *x=NULL,*y=NULL,*z=NULL,*xo=NULL,*yo=NULL,*zo=NULL;
      
  /* variables for multithreading procedure */
    
  ITG sys_cpus,*ithread=NULL;
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

  envloc = getenv("CCX_NPROC_CFD");
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

  /* initialization of additional fields */
  
  NNEW(x,double,*netet_);
  NNEW(y,double,*netet_);
  NNEW(z,double,*netet_);
  NNEW(nx,ITG,*netet_);
  NNEW(ny,ITG,*netet_);
  NNEW(nz,ITG,*netet_);
  
  netet=0;
  for(i=0;i<*netet_;i++){
    if(kontet[4*i]==0) continue;
    nx[netet]=netet+1;
    ny[netet]=netet+1;
    nz[netet]=netet+1;
    x[netet]=bc[4*i];
    y[netet]=bc[4*i+1];
    z[netet]=bc[4*i+2];
    iparent[netet]=i+1;
    netet++;
  }

  /* reallocate/allocate fields to the correct size
     (netet instead of netet_) */
  
  RENEW(x,double,netet);
  RENEW(y,double,netet);
  RENEW(z,double,netet);
  NNEW(xo,double,netet);
  NNEW(yo,double,netet);
  NNEW(zo,double,netet);
  RENEW(nx,ITG,netet);
  RENEW(ny,ITG,netet);
  RENEW(nz,ITG,netet);
  
  isiz=netet;
  cpypardou(xo,x,&isiz,&num_cpus);
  cpypardou(yo,y,&isiz,&num_cpus);
  cpypardou(zo,z,&isiz,&num_cpus);

  kflag=2;
  FORTRAN(dsort,(x,nx,&netet,&kflag));
  FORTRAN(dsort,(y,ny,&netet,&kflag));
  FORTRAN(dsort,(z,nz,&netet,&kflag));

  if(*nkf<num_cpus) num_cpus=*nkf;
    
  pthread_t tid[num_cpus];

  /* calculating the stiffness and/or mass matrix 
     (symmetric part) */

  x1=x;yy1=y;z1=z;xo1=xo;yo1=yo;zo1=zo;nx1=nx;ny1=ny;nz1=nz;planfa1=planfa;
  ifatet1=ifatet;netet1=&netet;kontet1=kontet;cotet1=cotet;iparent1=iparent;
  co1=co;konl1=konl;ratio1=ratio;ikf1=ikf;nkf1=nkf;
    
  /* create threads and wait */
    
  NNEW(ithread,ITG,num_cpus);
  for(i=0; i<num_cpus; i++)  {
    ithread[i]=i;
    pthread_create(&tid[i], NULL, (void *)interpoltetmt, (void *)&ithread[i]);
  }
  for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
    
  SFREE(ithread);
  
  return;

}

/* subroutine for multithreading of interpoltet */

void *interpoltetmt(ITG *i){

  ITG nkfa,nkfb,nkfdelta;
    
  // ceil -> floor

  nkfdelta=(ITG)floor(*nkf1/(double)num_cpus);
  nkfa=*i*nkfdelta+1;
  nkfb=(*i+1)*nkfdelta;
  // next line! -> all parallel sections
  if((*i==num_cpus-1)&&(nkfb<*nkf1)) nkfb=*nkf1;

  FORTRAN(interpoltet,(x1,yy1,z1,xo1,yo1,zo1,nx1,ny1,nz1,planfa1,ifatet1,
		       netet1,kontet1,cotet1,iparent1,co1,&nkfa,&nkfb,
		       konl1,ratio1,ikf1));

  return NULL;
}
