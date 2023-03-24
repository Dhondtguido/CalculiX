/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2015 Guido Dhondt                          */

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
#include "CalculiX.h"

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

void crackfrd(ITG *nk,ITG *ngraph,ITG *noddiam,double *cs,ITG *kode,ITG *inum,
	      ITG *nmethod,double *time,ITG *istep,ITG *iinc,ITG *mode,
	      char *description,char *set,ITG *nset,ITG *istartset,
	      ITG *iendset,ITG *ialset,char *jobnamec,char *output,
	      double *dkeqglob,double *wk1glob,double *wk2glob,double *wk3glob,
	      double *phiglob,double *dadnglob,double *dnglob,
	      double *acrackglob,double *xkeqminglob,double *xkeqmaxglob,
	      ITG *iincglob,double *domstepglob,double *rglob){

  /* stores the results in frd format

     iselect selects which nodes are to be stored:
     iselect=-1 means only those nodes for which inum negative
     ist, i.e. network nodes
     iselect=+1 means only those nodes for which inum positive
     ist, i.e. structural nodes
     iselect=0  means both of the above */
  
  FILE *f1;
  
  char m1[4]=" -1",m2[4]=" -2",m3[4]=" -3",fneig[132]="",filab[87]="";

  ITG icounter=0,nkcoords;

  ITG null,i,j,noutloc,iset,iselect,nfield[2],ifieldstate[15],icompstate[15],
    ioutall=0,nout,noutplus,noutmin,nstate_;

  double pi,oner,nullr,*xstaten=NULL;

  strcpy2(fneig,jobnamec,132);
  strcat(fneig,".frd");

  if((f1=fopen(fneig,"ab"))==NULL){
    printf("*EOR in frd: cannot open frd file for writing...");
    exit(0);
  }

  /* check whether all results have to be stored (also those
     corresponding to inactive nodes or elements) */
  
  if(strcmp1(&output[3],"a")==0) ioutall=1;
  
  pi=4.*atan(1.);
  null=0;
  nullr=0;
  oner=1.;

  /* nkcoords is the number of nodes at the time when 
     the nodal coordinates are stored in the frd file */

  nkcoords=*nk;

  /* determining nout, noutplus and noutmin 
     nout: number of structural and network nodes
     noutplus: number of structural nodes
     noutmin: number of network nodes */

  if(*nmethod!=0){
    nout=0;
    noutplus=0;
    noutmin=0;
    if(ioutall==0){
      for(i=0;i<*nk;i++){
	if(inum[i]==0) continue;
	nout++;
	if(inum[i]>0) noutplus++;
	if(inum[i]<0) noutmin++;
      }
    }else{
      for(i=0;i<*nk;i++){
	nout++;
	if(inum[i]>0) noutplus++;
	if(inum[i]<0) noutmin++;
      }
    }
  }else{
    nout=*nk;
  }
  
  /* storing the internal state variables in the nodes */
  
  iselect=1;
  nstate_=15;
    
  frdset(filab,set,&iset,istartset,iendset,ialset,
	 inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	 ngraph);
    
  frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	    &noutloc,description,kode,nmethod,f1,output,istep,iinc);

  fprintf(f1," -4  CT3D-MIS  %3" ITGFORMAT "    1\n",nstate_);
  fprintf(f1," -5  DOM_STEP    1    1    1    1\n");
  fprintf(f1," -5  DeltaKEQ    1    1    1    1\n");
  fprintf(f1," -5  KEQMIN      1    1    1    1\n");
  fprintf(f1," -5  KEQMAX      1    1    1    1\n");
  fprintf(f1," -5  K1WORST     1    1    1    1\n");
  fprintf(f1," -5  K2WORST     1    1    1    1\n");
  fprintf(f1," -5  K3WORST     1    1    1    1\n");
  fprintf(f1," -5  PHI(DEG)    1    1    1    1\n");
  fprintf(f1," -5  R           1    1    1    1\n");
  fprintf(f1," -5  DADN        1    1    1    1\n");
  fprintf(f1," -5  KTH         1    1    1    1\n");
  fprintf(f1," -5  INC         1    1    1    1\n");
  fprintf(f1," -5  CYCLES      1    1    1    1\n");
  fprintf(f1," -5  CRLENGTH    1    1    1    1\n");
  fprintf(f1," -5  DOM_SLIP    1    1    1    1\n");

  for(i=0;i<nstate_;i++){
    ifieldstate[i]=1;icompstate[i]=i;
  }
  nfield[0]=nstate_;

  NNEW(xstaten,double,nstate_**nk);
  
  for(i=0;i<*nk;i++){
    xstaten[nstate_*i]=domstepglob[i];
    xstaten[nstate_*i+1]=dkeqglob[i];
    xstaten[nstate_*i+2]=xkeqminglob[i];
    xstaten[nstate_*i+3]=xkeqmaxglob[i];
    xstaten[nstate_*i+4]=wk1glob[i];
    xstaten[nstate_*i+5]=wk2glob[i];
    xstaten[nstate_*i+6]=wk3glob[i];
    xstaten[nstate_*i+7]=phiglob[i]*180./pi;
    xstaten[nstate_*i+8]=rglob[i];
    xstaten[nstate_*i+9]=dadnglob[i];
    xstaten[nstate_*i+10]=nullr;
    xstaten[nstate_*i+11]=1.0*iincglob[i];
    xstaten[nstate_*i+12]=dnglob[i];
    xstaten[nstate_*i+13]=acrackglob[i];
    xstaten[nstate_*i+14]=nullr;
  }

  frdselect(xstaten,xstaten,&iset,&nkcoords,inum,m1,istartset,iendset,
	    ialset,ngraph,&nstate_,ifieldstate,icompstate,
	    nfield,&iselect,m2,f1,output,m3);

  SFREE(xstaten);
  
  fclose(f1);
  return;
  
}
