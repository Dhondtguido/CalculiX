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

static ITG *nelt1,*ia1,*ja1,*isym1,*itol1,*itmax1,*iter1=NULL,*igwk1=NULL,
    *lrgw1,*ligw1,*iwork1,*nestart1,*iunit1,*ierr1=NULL;

static double *b1,*x1,*a1,*tol1=NULL,*err1=NULL,*sb1,*sx1,*rgwk1=NULL,
    *rwork1;

void dgmresmain(ITG *nef,double *b,double *x,ITG *nelt,ITG *ia,ITG *ja,
		double *a,ITG *isym,ITG *itol,double *tol,ITG *itmax,
		ITG *iter,double *err,ITG *ierr,ITG *iunit,double *sb,
		double *sx,double *rgwk,ITG *lrgw,ITG *igwk,ITG *ligw,
		double *rwork,ITG *iwork,ITG *nestart,ITG *num_cpus,
		double *dgmrestol){

    ITG i;
      
    /* variables for multithreading procedure */
    
    ITG *ithread=NULL;
    
    pthread_t tid[*num_cpus];

    *itol=0;
    *isym=0;
    *itmax=0;
    *iunit=0;

    *ligw=20;
    
    *lrgw=131+16*nestart[1];

    NNEW(tol1,double,*num_cpus);
    for(i=0;i<*num_cpus;i++){tol1[i]=*dgmrestol;}
    NNEW(iter1,ITG,*num_cpus);
    NNEW(err1,double,*num_cpus);
    NNEW(ierr1,ITG,*num_cpus);
    NNEW(rgwk1,double,*lrgw**num_cpus);

    NNEW(igwk1,ITG,*ligw**num_cpus);
    for(i=0;i<*num_cpus;i++){
	igwk1[i**ligw]=10;
	igwk1[i**ligw+1]=10;
	igwk1[i**ligw+2]=0;
	igwk1[i**ligw+3]=1;
	igwk1[i**ligw+4]=10;
    }

    b1=b;x1=x;nelt1=nelt;ia1=ia;ja1=ja;a1=a;isym1=isym;itol1=itol;
    itmax1=itmax;iunit1=iunit;sb1=sb;sx1=sx;lrgw1=lrgw;ligw1=ligw;
    rwork1=rwork;iwork1=iwork;nestart1=nestart;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)dgmres1mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    

    *ierr=0;
    for(i=0;i<*num_cpus;i++){
	if(ierr1[i]>*ierr){*ierr=ierr1[i];}
    }
    SFREE(ithread);SFREE(tol1);SFREE(iter1);SFREE(err1);SFREE(ierr1);
    SFREE(rgwk1);SFREE(igwk1);
  
    return;

}

/* subroutine for multithreading of dgmres1 */

void *dgmres1mt(ITG *i){

    ITG n,indexrgwk,indexigwk;

    /* number of equations for this thread */

    n=nestart1[*i+1]-nestart1[*i];
    indexrgwk=*i**lrgw1;
    indexigwk=*i**ligw1;

    FORTRAN(dgmres1,(&n,&b1[nestart1[*i]],&x1[nestart1[*i]],nelt1,
		     ia1,&ja1[nestart1[*i]],a1,isym1,itol1,&tol1[*i],
		     itmax1,&iter1[*i],&err1[*i],&ierr1[*i],iunit1,sb1,sx1,
		     &rgwk1[indexrgwk],lrgw1,&igwk1[indexigwk],
		     ligw1,&rwork1[nestart1[*i]],iwork1));

    return NULL;
}
