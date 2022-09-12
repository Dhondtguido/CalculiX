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

void writenewmesh(ITG *nktet,ITG *netet_,double *cotet,ITG *iquad,
		  ITG *kontet,ITG *iedgmid,ITG *iedtet,ITG *mi,
		  char *matname,ITG *ithermal,char *jobnamec,
		  char *output,ITG *nmat){

  /* writing the new mesh in a mesh refinement calculation */

  FILE *f1;
  
  char *lakonnew=NULL,filabnew[5]="    ",*description=NULL,*set=NULL,
    fnewmesh[132]="",fneig[132]="",outputnew[5]="bin ";

  ITG nknew,nenew,*ipkonnew=NULL,*konnew=NULL,i,j,netet,*inum=NULL,
    nmethod=0,kode=1,*ielmatnew=NULL,nstate_=0,istep,iinc,mode=-1,
    noddiam=-1,*inotr=NULL,ntrans,*ielorien=NULL,norien,*ipneigh=NULL,
    *neigh=NULL,ngraph,nset,*istartset=NULL,*iendset=NULL,*ialset=NULL,
    mortar=0,*ielprop=NULL;

  double *conew=NULL,*v=NULL,*stn=NULL,*een=NULL,*t1=NULL,*fn=NULL,
    time=0.,*epn=NULL,*enern=NULL,*xstaten=NULL,*qfn=NULL,*trab=NULL,
    *orab=NULL,*stx=NULL,*vr=NULL,*vi=NULL,*stnr=NULL,*stni=NULL,
    *vmax=NULL,*stnmax=NULL,*veold=NULL,*ener=NULL,*cs=NULL,*eenmax=NULL,
    *fnr=NULL,*fni=NULL,*emn=NULL,*thicke=NULL,*qfx=NULL,*cdn=NULL,
    *cdnr=NULL,*cdni=NULL,*prop=NULL,*sti=NULL;

  strcpy2(fnewmesh,jobnamec,132);
  strcat(fnewmesh,".rfn");

  nknew=*nktet;
  nenew=*netet_;

  NNEW(conew,double,3*nknew);
  memcpy(&conew[0],&cotet[0],sizeof(double)*3*nknew);
   
  NNEW(ipkonnew,ITG,nenew);
  NNEW(konnew,ITG,10*nenew);
  NNEW(lakonnew,char,8*nenew);
  NNEW(ielmatnew,ITG,mi[2]*nenew);

  for(i=0;i<nenew;i++){
    ipkonnew[i]=-1;
  }
  
  if(*iquad==0) {
    netet=0;
    for(i=0;i<nenew;i++){
      if(kontet[i*4]==0) continue;
      ipkonnew[i]=4*netet;
      for(j=0;j<4;j++) {
	konnew[4*netet+j]=kontet[4*i+j];
      }
      strcpy1(&lakonnew[8*i],"C3D4    ",8);
      netet++;
    }
    RENEW(konnew,ITG,4*netet);
  }
  else if(*iquad==1) {
    netet=0;
    for(i=0;i<nenew;i++){
      if(kontet[i*4]==0) continue;
      ipkonnew[i]=10*netet;
      for(j=0;j<4;j++) {
	konnew[10*netet+j]=kontet[4*i+j];
      }
      for(j=4;j<10;j++){
	konnew[10*netet+j]=iedgmid[iedtet[6*i+j-4]-1];
      } 
      
      strcpy1(&lakonnew[8*i],"C3D10   ",8);
      netet++;
    }
    RENEW(konnew,ITG,10*netet);
  }

  strcpy2(fneig,fnewmesh,132);
  strcat(fneig,".frd");

  if((f1=fopen(fneig,"wb"))==NULL){
    printf(" *ERROR in frd: cannot open frd file for writing...");
    exit(0);
  }
  
  fclose(f1);
  
  /* creating the tetrahedral mesh in frd format*/
   
  frd(conew,&nknew,konnew,ipkonnew,lakonnew,&nenew,v,stn,inum,&nmethod,
      &kode,filabnew,een,t1,fn,&time,epn,ielmatnew,matname,enern,xstaten,
      &nstate_,&istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
      &ntrans,orab,ielorien,&norien,description,ipneigh,neigh,
      mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,&nenew,
      cs,set,&nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
      thicke,fnewmesh,outputnew,qfx,cdn,&mortar,cdnr,cdni,nmat,ielprop,
      prop,sti);
  
  strcat(fnewmesh,".frd");
  if((f1=fopen(fnewmesh,"ab"))==NULL){
    printf(" *ERROR in frd: cannot open frd file for writing...");
    exit(0);
  }
  fprintf(f1," 9999\n");
  fclose(f1);

  SFREE(conew);SFREE(ipkonnew);SFREE(konnew);SFREE(lakonnew);SFREE(ielmatnew);
  
  return;
  
}
