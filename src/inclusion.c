/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2021 Guido Dhondt                     */

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

static ITG *nacti1,num_cpus_loc;

static double *al1,*alnew1,*gmatrix1;

void inclusion(double *gmatrix,double *cvec,ITG *iacti,ITG *nacti,double *fric,
	      double *atol,double *rtol,double *alglob,ITG *kitermax,
	      double *auw,ITG *jqw,ITG *iroww,ITG *nslavs,double *al,
	      double *alnew,double *r,double *omega,ITG *masslesslinear,
	       double *fullr,ITG *num_cpus){

  /* determining the RHS of the global system for massless contact */
  
  ITG iscvg=0,icont=0,i,j,in,it1,it2,irow,*ithread=NULL;

  double err,alsize,altan,altanmax,ratio,value;

  pthread_t tid[*num_cpus];

  /* determine the relaxation vector */

  if(*masslesslinear>0){
    for(i=0;i<3**nslavs;i++){
      if(iacti[i]!=0){
	r[iacti[i]-1]=fullr[i];
      }
    }
  }else{
    FORTRAN(relaxval_al,(r,gmatrix,nacti));
  }

  /* check that num_cpus does not exceed nacti */

  if(*nacti<*num_cpus){
    num_cpus_loc=*nacti;
  }else{
    num_cpus_loc=*num_cpus;
  }

  while((icont<=*kitermax)&&(iscvg==0)){
  
    /* alnew=G*al */

    nacti1=nacti;alnew1=alnew;al1=al;gmatrix1=gmatrix;

    /* create threads and wait */
	
    NNEW(ithread,ITG,num_cpus_loc);
    for(i=0; i<num_cpus_loc; i++)  {
      ithread[i]=i;
      pthread_create(&tid[i], NULL, (void *)gmatrixtimesalmt, (void *)&ithread[i]);
    }
    for(i=0; i<num_cpus_loc; i++)  pthread_join(tid[i], NULL);

    SFREE(ithread);
    
    /* alnew=al-omega*r*(alnew+c) */

    for(i=0;i<*nacti;i++){
      alnew[i]=al[i]-(*omega)*r[i]*(alnew[i]+cvec[i]);
    }

    /* projection for normal and tangential contact */

    err=0.;
    alsize=0.;

    for(i=0;i<*nslavs;i++){
    
      if(iacti[3*i]==0) continue;
    
      in=iacti[3*i]-1;
      it1=in+1;
      it2=in+2;

      /* F_normal */

      if(alnew[in]<0.) alnew[in]=0.;

      /* F_tangential */

      altan=sqrt(alnew[it1]*alnew[it1]+alnew[it2]*alnew[it2]);

      /* F_normal * fric_coefficient */

      altanmax=fric[i]*alnew[in];

      /* evaluate stick or slip */

      if(altan>altanmax){

	/* if slip, adjust velocities according to Coulomb model */

	ratio=altanmax/altan;
	alnew[it1]=ratio*alnew[it1];
	alnew[it2]=ratio*alnew[it2];
      }

      /* evaluate change of solution */

      err+=(alnew[in]-al[in])*(alnew[in]-al[in])+
	(alnew[it1]-al[it1])*(alnew[it1]-al[it1])+
	(alnew[it2]-al[it2])*(alnew[it2]-al[it2]);
      alsize+=al[in]*al[in]+al[it1]*al[it1]+al[it2]*al[it2];

      /* update al */

      al[in]=alnew[in];
      al[it1]=alnew[it1];
      al[it2]=alnew[it2];
    }

    err=sqrt(err);
    alsize=sqrt(alsize);

    /* check for convergence */

    if(err<=(alsize**rtol+(*atol))){
      iscvg=1;
    }

    icont++;
  }

  if(icont>*kitermax){
    printf(" *WARNING in inclusion: maximum number of iterations\n");
    printf("          for massless contact reached: %d\n",*kitermax);
    printf("          with error norm: %e\n\n",err);
  }

  /* calculate the contact force in the global rectangular system:
     alglob=Wb*al */

  for(i=0;i<3**nslavs;i++){
    if(iacti[i]!=0){
      for(j=jqw[i]-1;j<jqw[i+1]-1;j++){
	value=auw[j];
	irow=iroww[j]-1;
	alglob[irow]+=value*al[iacti[i]-1];
      }
    }
  }
  
  return;
}

void *gmatrixtimesalmt(ITG *i){

  /* full matrix times vector multiplication */

  ITG idelta,ia,ib,j,k;
  
  idelta=(ITG)ceil(*nacti1/(double)num_cpus_loc);
  ia=*i*idelta;
  ib=(*i+1)*idelta;
  if(ib>*nacti1) ib=*nacti1;
    
  for(j=ia;j<ib;j++){
    alnew1[j]=0.;
    for(k=0;k<*nacti1;k++){
      alnew1[j]+=gmatrix1[k**nacti1+j]*al1[k];
    }
  }
  
  return NULL;
}
