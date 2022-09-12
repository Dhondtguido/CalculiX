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

void storecontactdof(ITG *nope,ITG *nactdof, ITG *mt, ITG *konl, ITG **ikactcontp, 
		     ITG *nactcont,ITG *nactcont_, double *bcont, double *fnl, 
		     ITG *ikmpc, ITG *nmpc, ITG *ilmpc,ITG *ipompc, ITG *nodempc, 
		     double *coefmpc){

  ITG j,j1,jdof,id,k,l,ist,index,node,ndir,*ikactcont=*ikactcontp;

  for(j=0;j<*nope;j++){
    for(j1=0;j1<3;j1++){
      jdof=nactdof[*mt*(konl[j]-1)+j1+1];
      if(jdof>0){
	      
	jdof--;
	FORTRAN(nident,(ikactcont,&jdof,nactcont,&id));
	do{
	  if(id>0){
	    if(ikactcont[id-1]==jdof){
	      break;
	    }
	  }
	  (*nactcont)++;
	  if(*nactcont>*nactcont_){
	    *nactcont_=(ITG)(1.1**nactcont_);
	    RENEW(ikactcont,ITG,*nactcont_);
	  }
	  k=*nactcont-1;
	  l=k-1;
	  while(k>id){
	    ikactcont[k--]=ikactcont[l--];
	  }
	  ikactcont[id]=jdof;
	  break;
	}while(1);
	      
	bcont[jdof]-=fnl[3*j+j1];
      }else{
	jdof=8*(konl[j]-1)+j1+1;
	FORTRAN(nident,(ikmpc,&jdof,nmpc,&id));
	if(id>0){
	  if(ikmpc[id-1]==jdof){
	    id=ilmpc[id-1];
	    ist=ipompc[id-1];
	    index=nodempc[3*ist-1];
	    if(index==0) continue;
	    do{
	      node=nodempc[3*index-3];
	      ndir=nodempc[3*index-2];
	      jdof=nactdof[*mt*(node-1)+ndir];
	      if(jdof>0){
			      
		jdof--;
		FORTRAN(nident,(ikactcont,&jdof,nactcont,&id));
		do{
		  if(id>0){
		    if(ikactcont[id-1]==jdof){
		      break;
		    }
		  }
		  (*nactcont)++;
		  if(*nactcont>*nactcont_){
		    *nactcont_=(ITG)(1.1**nactcont_);
		    RENEW(ikactcont,ITG,*nactcont_);
		  }
		  k=*nactcont-1;
		  l=k-1;
		  while(k>id){
		    ikactcont[k--]=ikactcont[l--];
		  }
		  ikactcont[id]=jdof;
		  break;
		}while(1);
			      
		/*			      bcont[jdof]+=coefmpc[index-1]*
					      fnl[3*j+j1]/coefmpc[ist-1];*/
		bcont[jdof]-=coefmpc[index-1]*
		  fnl[3*j+j1]/coefmpc[ist-1];
	      }
	      index=nodempc[3*index-1];
	      if(index==0) break;
	    }while(1);
	  }
	}
      }
    }
  }

  *ikactcontp=ikactcont;
  
  return;
}

