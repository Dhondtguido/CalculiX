/*     CalculiX - A 3-dimensional finite element program                   */
/*              Copyright (C) 1998-2015 Guido Dhondt                       */
/*	        Copyright (c) 2016 Jean-Marie Verdun 			 */

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

#ifdef TAUCS

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CalculiX.h"
#include "tau.h"
#include <taucs.h>

taucs_ccs_matrix aa[1];
void* F=NULL;
char* taufactor[]={ "taucs.factor.LLT=true","taucs.factor.mf=true",
  "taucs.factor.ordering=amd",NULL };
char* taufactorooc[]={ "taucs.factor.LLT=true","taucs.factor.mf=true","taucs.factor.ordering=amd","taucs.ooc=true",
                  "taucs.ooc.basename=/tmp/scratch",
                  "taucs.ooc.memory=500000000.0",NULL };
char* tausolve[]={ "taucs.factor=false",NULL };
char* tausolveooc[]={"taucs.factor=false","taucs.ooc.basename=/tmp/scratch","taucs.ooc.memory=500000000.0",NULL };
ITG *irowtau=NULL,*pointtau=NULL;
double *autau=NULL;
ITG* perm;


void tau_factor(double *ad, double **aup, double *adb, double *aub, 
                double *sigma,ITG *icol, ITG **irowp, 
                ITG *neq, ITG *nzs){

  ITG i,j,k,l,*irow=NULL;
  long long ndim;
  double *au=NULL;
  double memory_mb = -1.0;
  ITG    mb = -1;
  ITG ret;

  printf(" Factoring the system of equations using TAUCS\n\n");

  taucs_logfile("stdout");

  au=*aup;
  irow=*irowp;

  ndim=*neq+*nzs;

  NNEW(autau,double,ndim);
  NNEW(irowtau,ITG,ndim);
  NNEW(pointtau,ITG,*neq+1);

  k=ndim;
  l=*nzs;

  if(*sigma==0.){
    pointtau[*neq]=ndim;
    for(i=*neq-1;i>=0;--i){
      for(j=0;j<icol[i];++j){
	irowtau[--k]=irow[--l]-1;
	autau[k]=au[l];
      }
      pointtau[i]=--k;
      irowtau[k]=i;
      autau[k]=ad[i];
    }
  }
  else{
    pointtau[*neq]=ndim;
    for(i=*neq-1;i>=0;--i){
      for(j=0;j<icol[i];++j){
	irowtau[--k]=irow[--l]-1;
	autau[k]=au[l]-*sigma*aub[l];
      }
      pointtau[i]=--k;
      irowtau[k]=i;
      autau[k]=ad[i]-*sigma*adb[i];
    }
  }

  /* defining the matrix */

  aa->n = *neq;
  aa->m = *neq;
  aa->flags  = TAUCS_SYMMETRIC | TAUCS_LOWER | TAUCS_DOUBLE;
  aa->colptr = pointtau;
  aa->rowind = irowtau;
  aa->values.d = autau;

  if(*neq<5000){
      taucs_linsolve(aa,&F,0,NULL,NULL,taufactor,NULL);
  }
  else{
      ret = taucs_linsolve(aa,&F,0,NULL,NULL,taufactorooc,NULL);
           
      printf(" Return Code from Factoring %" ITGFORMAT "\n\n",ret);

  }
  
  *aup=au;
  *irowp=irow;

  return;
}

void tau_solve(double *b,ITG *neq){

  ITG i;
  /*static ITG j;*/
  double *x=NULL;
  ITG ret;

  NNEW(x,double,*neq);
  
  if(*neq<5000){
      taucs_linsolve(aa,&F,1,x,b,tausolve,NULL);
  }
  else{
      ret =  taucs_linsolve(aa,&F,1,x,b,tausolveooc,NULL);
      
//      ret = taucs_ooc_solve_llt(F, x, b);
       
      printf(" Return Code from Solving %" ITGFORMAT "\n\n",ret);
      
//      taucs_io_delete(F);

  }

  for(i=0;i<=*neq-1;++i){
    b[i]=x[i];
  }
  SFREE(x);/*
  if (mb > 0)
    memory_mb = (double) mb;
  else
    memory_mb = ((double) (-mb)) * taucs_available_memory_size()/1048576.0;  
  */
  /*j++;printf("%" ITGFORMAT "\n",j);*/

  return;
}

void tau_cleanup(){

//  taucs_linsolve(NULL,&F,0,NULL,NULL,NULL,NULL);
  SFREE(pointtau);
  SFREE(irowtau);
  SFREE(autau);

  return;
}

void tau(double *ad, double **aup, double *adb, double *aub, double *sigma,
         double *b, ITG *icol, ITG **irowp, 
         ITG *neq, ITG *nzs){

  if(*neq==0) return;

     
  tau_factor(ad,aup,adb,aub,sigma,icol,irowp, 
             neq,nzs);

  tau_solve(b,neq);
  
  tau_cleanup();
  

  return;
}

#endif
