/*     CalculiX - A 3-dimensional finite element program                   */
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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

void radcyc(ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
	    double *cs, ITG *mcs, ITG *nkon,ITG *ialset, ITG *istartset,
            ITG *iendset,ITG **kontrip,ITG *ntri,
            double **cop, double **voldp,ITG *ntrit, ITG *inocs,
            ITG *mi){

  /* duplicates triangular faces for cyclic radiation conditions */

  char *filab=NULL;

  ITG i,is,nsegments,idtie,nkt,icntrl,imag=0,*kontri=NULL,mt=mi[1]+1,
     node,i1,i2,nope,iel,indexe,j,k,ielset,node1,node2,node3,l,jj;

  double *vt=NULL,*fnt=NULL,*stnt=NULL,*eent=NULL,*qfnt=NULL,t[3],theta,
     pi,*v=NULL,*fn=NULL,*stn=NULL,*een=NULL,*qfn=NULL,*co=NULL,
     *vold=NULL,*emnt=NULL,*emn=NULL;

  pi=4.*atan(1.);
  
  kontri=*kontrip;co=*cop;vold=*voldp;

  /* determining the maximum number of sectors */

  nsegments=1;
  for(j=0;j<*mcs;j++){
      if(cs[17*j]>nsegments) nsegments=(ITG)(cs[17*j]);
  }

  /* assigning nodes and elements to sectors */

  ielset=cs[12];
  if((*mcs!=1)||(ielset!=0)){
    for(i=0;i<*nk;i++) inocs[i]=-1;
  }

  for(i=0;i<*mcs;i++){
    is=cs[17*i+4];
    if(is==1) continue;
    ielset=cs[17*i+12];
    if(ielset==0) continue;
    for(i1=istartset[ielset-1]-1;i1<iendset[ielset-1];i1++){
      if(ialset[i1]>0){
        iel=ialset[i1]-1;
        if(ipkon[iel]<0) continue;
        indexe=ipkon[iel];
        if(strcmp1(&lakon[8*iel+3],"2")==0)nope=20;
        else if (strcmp1(&lakon[8*iel+3],"8")==0)nope=8;
        else if (strcmp1(&lakon[8*iel+3],"10")==0)nope=10;
        else if (strcmp1(&lakon[8*iel+3],"4")==0)nope=4;
        else if (strcmp1(&lakon[8*iel+3],"15")==0)nope=15;
        else {nope=6;}
        for(i2=0;i2<nope;++i2){
          node=kon[indexe+i2]-1;
          inocs[node]=i;
        }
      }
      else{
        iel=ialset[i1-2]-1;
        do{
          iel=iel-ialset[i1];
          if(iel>=ialset[i1-1]-1) break;
          if(ipkon[iel]<0) continue;
          indexe=ipkon[iel];
          if(strcmp1(&lakon[8*iel+3],"2")==0)nope=20;
          else if (strcmp1(&lakon[8*iel+3],"8")==0)nope=8;
          else if (strcmp1(&lakon[8*iel+3],"10")==0)nope=10;
          else if (strcmp1(&lakon[8*iel+3],"4")==0)nope=4;
          else if (strcmp1(&lakon[8*iel+3],"15")==0)nope=15;
          else {nope=6;}
          for(i2=0;i2<nope;++i2){
            node=kon[indexe+i2]-1;
            inocs[node]=i;
          }
        }while(1);
      }
    } 
  }

  /* duplicating triangular faces 
     only those faces are duplicated the nodes of which belong to
     the same cyclic symmetry. non-integer cyclic symmety numbers are
     reduced to the next lower integer. */

  *ntrit=nsegments**ntri;
  RENEW(kontri,ITG,4**ntrit);
  for(i=4**ntri;i<4**ntrit;i++) kontri[i]=0;

  for(i=0;i<*ntri;i++){
    node1=kontri[4*i];
    if(inocs[node1-1]<0) continue;
    idtie=inocs[node1-1];
    node2=kontri[4*i+1];
    if((inocs[node2-1]<0)||(inocs[node2-1]!=idtie)) continue;
    node3=kontri[4*i+2];
    if((inocs[node3-1]<0)||(inocs[node3-1]!=idtie)) continue;
    idtie=cs[17*idtie];
    for(k=1;k<idtie;k++){
      j=i+k**ntri;
      kontri[4*j]=node1+k**nk;
      kontri[4*j+1]=node2+k**nk;
      kontri[4*j+2]=node3+k**nk;
      kontri[4*j+3]=kontri[4*i+3];
    }
  }

  RENEW(co,double,3**nk*nsegments);
  RENEW(vold,double,mt**nk*nsegments);
  nkt=*nk*nsegments;
      
  /* generating the coordinates for the other sectors */
  
  icntrl=1;
  
  FORTRAN(rectcyl,(co,v,fn,stn,qfn,een,cs,nk,&icntrl,t,filab,&imag,mi,emn));
  
  for(jj=0;jj<*mcs;jj++){
    is=(ITG)(cs[17*jj]);
    for(i=1;i<is;i++){
      
      theta=i*2.*pi/cs[17*jj];
      
      for(l=0;l<*nk;l++){
        if(inocs[l]==jj){
	  co[3*l+i*3**nk]=co[3*l];
	  co[1+3*l+i*3**nk]=co[1+3*l]-theta;
	  co[2+3*l+i*3**nk]=co[2+3*l];
        }
      }
    }
  }

  icntrl=-1;
    
  FORTRAN(rectcyl,(co,vt,fnt,stnt,qfnt,eent,cs,&nkt,&icntrl,t,filab,
		   &imag,mi,emnt));

  *kontrip=kontri;*cop=co;*voldp=vold;

  return;
}

