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

void writeoldmesh(ITG *nk,ITG *ne,double *co,ITG *ipkon,
		  ITG *kon,char *lakon,ITG *mi,
		  char *matname,ITG *ithermal,char *jobnamec,
		  char *output,ITG *nmat){

  /* writing the old mesh in a mesh refinement calculation */

  FILE *f1;
  
  char filabtmp[5]="    ",*description=NULL,*set=NULL,
    foldmesh[132]="",fneig[132]="";

  ITG *inum=NULL,
    nmethod=0,kode=1,*ielmatold=NULL,nstate_=0,istep,iinc,mode=-1,
    noddiam=-1,*inotr=NULL,ntrans,*ielorien=NULL,norien,*ipneigh=NULL,
    *neigh=NULL,ngraph,nset,*istartset=NULL,*iendset=NULL,*ialset=NULL,
    mortar=0,*ielprop=NULL;

  double *v=NULL,*stn=NULL,*een=NULL,*t1=NULL,*fn=NULL,
    time=0.,*epn=NULL,*enern=NULL,*xstaten=NULL,*qfn=NULL,*trab=NULL,
    *orab=NULL,*stx=NULL,*vr=NULL,*vi=NULL,*stnr=NULL,*stni=NULL,
    *vmax=NULL,*stnmax=NULL,*veold=NULL,*ener=NULL,*cs=NULL,*eenmax=NULL,
    *fnr=NULL,*fni=NULL,*emn=NULL,*thicke=NULL,*qfx=NULL,*cdn=NULL,
    *cdnr=NULL,*cdni=NULL,*prop=NULL,*sti=NULL;

  strcpy2(foldmesh,jobnamec,132);
  strcat(foldmesh,".urf");

  strcpy2(fneig,foldmesh,132);
  strcat(fneig,".frd");

  if((f1=fopen(fneig,"wb"))==NULL){
    printf(" *ERROR in frd: cannot open frd file for writing...");
    exit(0);
  }
  
  fclose(f1);

  NNEW(ielmatold,ITG,mi[2]**ne);
  
  /* creating the tetrahedral mesh in frd format*/
   
  frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,&nmethod,
      &kode,filabtmp,een,t1,fn,&time,epn,ielmatold,matname,enern,xstaten,
      &nstate_,&istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
      &ntrans,orab,ielorien,&norien,description,ipneigh,neigh,
      mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
      cs,set,&nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
      thicke,foldmesh,output,qfx,cdn,&mortar,cdnr,cdni,nmat,ielprop,
      prop,sti);
  
  strcat(foldmesh,".frd");
  if((f1=fopen(foldmesh,"ab"))==NULL){
    printf(" *ERROR in frd: cannot open frd file for writing...");
    exit(0);
  }
  fprintf(f1," 9999\n");
  fclose(f1);

  SFREE(ielmatold);
  
  return;
  
}
