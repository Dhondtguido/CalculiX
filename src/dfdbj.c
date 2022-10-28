/*
  CalculiX - A 3-dimensional finite element program
  Copyright (C) 1998-2015 Guido Dhondt

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation(version 2);
     

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
/*                                                                       */
/*     author: Reinhold Fischer                                          */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

void dfdbj(double *bcont,double **dbcontp,ITG *neq,ITG *nope,ITG *konl,
	   ITG* nactdof,double *s,double *z,ITG *ikmpc,ITG *ilmpc,
	   ITG *ipompc,ITG *nodempc,ITG *nmpc,double *coefmpc,
	   double *fnl,ITG *nev,ITG **ikactcontp,ITG **ilactcontp,
           ITG *nactcont,ITG *nactcont_,ITG *mi, ITG *cyclicsymmetry,
           ITG *izdof, ITG *nzdof){

  ITG j,j1,jdof,kdof,k,k1,l,id,index,ist,id1,ist1,index1,id2,ist2,index2,
    jdbcontcol,i1,i3,i4,mt=mi[1]+1,im,*ikactcont=*ikactcontp,
    *ilactcont=*ilactcontp,kdofm1;

  double d1,sl,*dbcont=*dbcontp;
  
  for(j=0; j<*nope; j++){
    i1=mt*(konl[j]-1)+1;
    for(j1=0; j1<3; j1++){
      jdof=nactdof[i1+j1];
      if(jdof>0){
	jdof--;
	FORTRAN(nident,(ikactcont,&jdof,nactcont,&id));
	do{
	  if(id>0){
	    if(ikactcont[id-1]==jdof){
	      jdbcontcol=ilactcont[id-1];
	      break;
	    }
	  }
	  (*nactcont)++;
	  if(*nactcont>*nactcont_){
	    *nactcont_=(ITG)(1.1**nactcont_);
	    RENEW(ikactcont,ITG,*nactcont_);
	    RENEW(ilactcont,ITG,*nactcont_);
	    RENEW(dbcont,double,*nev**nactcont_);
	  }
	  k=*nactcont-1;
	  l=k-1;
	  while(k>id){
	    ikactcont[k]=ikactcont[l];
	    ilactcont[k--]=ilactcont[l--];
	  }
	  jdbcontcol=*nactcont;
	  ikactcont[id]=jdof;
	  ilactcont[id]=*nactcont;
	  //		  memset(&dbcont[(*nactcont-1)**nev],0,sizeof(double)**nev);
	  DMEMSET(dbcont,(*nactcont-1)**nev,*nactcont**nev,0.);
	  break;
	}while(1);
	bcont[jdof]-=fnl[j*3+j1];
	i4=(jdbcontcol-1)**nev;
	i3=(3*j+j1);
	for(k=0; k<*nope; k++){
	  for(k1=0; k1<3; k1++){
	    sl=s[(3*k+k1)*60+i3];
	    kdof=nactdof[mt*(konl[k]-1)+k1+1];
	    if(kdof>0){
	      if(!(*cyclicsymmetry)){
		for(l=0; l<*nev; l++){
		  dbcont[i4+l]-=sl*z[(long long)l**neq+kdof-1];
		}
	      }else{
		kdofm1=kdof-1;
		FORTRAN(nident,(izdof,&kdofm1,nzdof,&id));
		if(id!=0){
		  if(izdof[id-1]==kdofm1){
		    for(l=0; l<*nev; l++){
		      dbcont[i4+l]-=sl*z[l**nzdof+id-1];
		    }
		  }else{printf(" *ERROR in dfdbj\n");FORTRAN(stop,());}
		}else{printf(" *ERROR in dfdbj\n");FORTRAN(stop,());}
	      }
	    }
	    else{
	      kdof=8*(konl[k]-1)+k1+1;
	      FORTRAN(nident,(ikmpc,&kdof,nmpc,&id));
	      if(id>0){
		id--;
		if(ikmpc[id]==kdof){
		  id=ilmpc[id];
		  ist=ipompc[id-1];
		  ist--;
		  index=nodempc[ist*3+2];
		  if(index==0) continue;
		  index--;
		  do{
		    kdof=nactdof[mt*(nodempc[index*3]-1)+nodempc[index*3+1]];
		    d1=sl*coefmpc[index]/coefmpc[ist];
		    if(kdof>0){
		      if(!(*cyclicsymmetry)){
			for(l=0; l<*nev; l++){
			  dbcont[i4+l]+=d1*z[(long long)l**neq+kdof-1];
			}
		      }
		    }else{
		      kdofm1=kdof-1;
		      FORTRAN(nident,(izdof,&kdofm1,nzdof,&id));
		      if(id!=0){
			if(izdof[id-1]==kdofm1){
			  for(l=0; l<*nev; l++){
			    dbcont[i4+l]+=d1*z[l**nzdof+id-1];
			  }
			}else{printf(" *ERROR in dfdbj\n");FORTRAN(stop,());}
		      }else{printf(" *ERROR in dfdbj\n");FORTRAN(stop,());}
		    }
		    index=nodempc[index*3+2];
		    if(index==0) break;
		    index--;
		  }while(1);
		}
	      }
	    }
	  }
	}
      }
      else{
	jdof=8*(konl[j]-1)+j1+1;
	FORTRAN(nident,(ikmpc,&jdof,nmpc,&id1));
	if(id1>0){
	  id1--;
	  if(ikmpc[id1]==jdof){
	    id1=ilmpc[id1];
	    ist1=ipompc[id1-1];
	    ist1--;
	    index1=nodempc[ist1*3+2];
	    if(index1==0) continue;
	    index1--;
	    do{
	      jdof=nactdof[mt*(nodempc[index1*3]-1)+nodempc[index1*3+1]];
	      if(jdof>0){
		jdof--;
		FORTRAN(nident,(ikactcont,&jdof,nactcont,&id));
		do{
		  if(id>0){
		    if(ikactcont[id-1]==jdof){
		      jdbcontcol=ilactcont[id-1];
		    }
		  }
		  (*nactcont)++;
		  if(*nactcont>*nactcont_){
		    *nactcont_=(ITG)(1.1**nactcont_);
		    RENEW(ikactcont,ITG,*nactcont_);
		    RENEW(ilactcont,ITG,*nactcont_);
		    RENEW(dbcont,double,*nev**nactcont_);
		  }
		  k=*nactcont-1;
		  l=k-1;
		  do{
		    ikactcont[k]=ikactcont[l];
		    ilactcont[k--]=ilactcont[l--];
		  }while(k>id);
		  jdbcontcol=*nactcont;
		  ikactcont[id]=jdof;
		  ilactcont[id]=*nactcont;
		  //				  memset(&dbcont[(*nactcont-1)**nev],0,sizeof(double)**nev);		  
		  DMEMSET(dbcont,(*nactcont-1)**nev,*nactcont**nev,0.);
		  break;
		}while(1);
		bcont[jdof]+=coefmpc[index1]*fnl[j*3+j1]/coefmpc[ist1];
		i4=(jdbcontcol-1)**nev;
		i3=(3*j+j1);
		for(k=0; k<*nope; k++){
		  for(k1=0; k1<3; k1++){
		    sl=s[(3*k+k1)*60+i3];
		    kdof=nactdof[mt*(konl[k]-1)+k1+1];
		    if(kdof>0){
		      d1=sl*coefmpc[index1]/coefmpc[ist1];
		      if(!(*cyclicsymmetry)){
			for(l=0; l<*nev; l++){
			  dbcont[i4+l]+=d1*z[(long long)l**neq+kdof-1];
			}
		      }else{
			kdofm1=kdof-1;
			FORTRAN(nident,(izdof,&kdofm1,nzdof,&id));
			if(id!=0){
			  if(izdof[id-1]==kdofm1){
			    for(l=0; l<*nev; l++){
			      dbcont[i4+l]+=d1*z[l**nzdof+id-1];
			    }
			  }else{printf(" *ERROR in dfdbj\n");FORTRAN(stop,());}
			}else{printf(" *ERROR in dfdbj\n");FORTRAN(stop,());}
		      }
		    }
		    else{
		      kdof=8*(konl[k]-1)+k1+1;
		      FORTRAN(nident,(ikmpc,&kdof,nmpc,&id2));
		      if(id2>0){
			id2--;
			if(ikmpc[id2]==kdof){
			  id2=ilmpc[id2];
			  ist2=ipompc[id2-1];
			  ist2--;
			  index2=nodempc[ist2*3+2];
			  if(index2==0) continue;
			  index2--;
			  do{
			    kdof=nactdof[mt*(nodempc[index2*3]-1)+nodempc[index2*3+1]];
			    if(kdof>0){
			      d1=sl*coefmpc[index1]*coefmpc[index2]/(coefmpc[ist1]*coefmpc[ist2]);
			      if(!(*cyclicsymmetry)){
				for(l=0; l<*nev; l++){
				  dbcont[i4+l]-=d1*z[(long long)l**neq+kdof-1];
				}
			      }else{
				kdofm1=kdof-1;
				FORTRAN(nident,(izdof,&kdofm1,nzdof,&id));
				if(id!=0){
				  if(izdof[id-1]==kdofm1){
				    for(l=0; l<*nev; l++){
				      dbcont[i4+l]-=d1*z[l**nzdof+id-1];
				    }
				  }else{printf(" *ERROR in dfdbj\n");FORTRAN(stop,());}
				}else{printf(" *ERROR in dfdbj\n");FORTRAN(stop,());}
			      }
			    }
			    index2=nodempc[index2*3+2];
			    if(index2==0) break;
			    index2--;
			  }while(1);
			}
		      }
		    }
		  }
		}
	      }
	      index1=nodempc[index1*3+2];
	      if(index1==0) break;
	      index1--;
	    }while(1);
	  }
	}
      }
    }
  }
  *dbcontp=dbcont;
  *ikactcontp=ikactcont;
  *ilactcontp=ilactcont;
}

/*!
  !     CalculiX - A 3-dimensional finite element program
  !              Copyright (C) 1998-2015 Guido Dhondt
  !
  !     This program is free software; you can redistribute it and/or
  !     modify it under the terms of the GNU General Public License as
  !     published by the Free Software Foundation(version 2);
  !     
  !
  !     This program is distributed in the hope that it will be useful,
  !     but WITHOUT ANY WARRANTY; without even the implied warranty of 
  !     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
  !     GNU General Public License for more details.
  !
  !     You should have received a copy of the GNU General Public License
  !     along with this program; if not, write to the Free Software
  !     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
  !
  subroutine dfdbj(bcont,dbcont,neq,nope,konl,nactdof,s,z,
  &  ikmpc,ilmpc,ipompc,nodempc,nmpc,coefmpc,fnl,nev,iactcont,
  &  nactcont)
  !
  !     calculates the derivative of the contact forces with respect
  !     to the modal variables
  !
  implicit none
  !
  integer j,j1,neq,nope,konl(*),nactdof(0:3,*),jdof,kdof,
  &  k,k1,l,id,ikmpc(*),ilmpc(*),ipompc(*),nodempc(3,*),nmpc,
  &  index,ist,id1,ist1,index1,id2,ist2,index2,nev ,iactcont(*),
  &  nactcont,jdofcont
  !
  real*8 bcont(*),dbcont(nev,*),s(60,60),z(neq,*),coefmpc(*),
  &  fnl(3,9)
  !
  do j=1,nope
  do j1=1,3
  jdof=nactdof(j1,konl(j))
  if(jdof.ne.0) then
  call nident(iactcont,jdof,nactcont,id)
  jdofcont=0
  if(id.gt.0)then
  if(iactcont(id).eq.jdof) then
  jdofcont=id
  endif
  endif
  if(jdofcont.eq.0) then
  nactcont=nactcont+1
  do k=nactcont,id+2,-1
  iactcont(k)=iactcont(k-1)
  do l=1,nev
  dbcont(l,k)=dbcont(l,k-1)
  enddo
  enddo
  jdofcont=id+1
  iactcont(jdofcont)=jdof
  do l=1,nev
  dbcont(l,jdofcont)=0.d0
  enddo
  endif
  bcont(jdof)=bcont(jdof)-fnl(j1,j)
  do k=1,nope
  do k1=1,3
  kdof=nactdof(k1,konl(k))
  if(kdof.ne.0) then
  do l=1,nev
  dbcont(l,jdofcont)=dbcont(l,jdofcont)-
  &                       s(3*(j-1)+j1,3*(k-1)+k1)*z(kdof,l)
  enddo
  else
  kdof=8*(konl(k)-1)+k1
  call nident(ikmpc,kdof,nmpc,id)
  if(id.gt.0) then
  if(ikmpc(id).eq.kdof) then
  id=ilmpc(id)
  ist=ipompc(id)
  index=nodempc(3,ist)
  if(index.eq.0) cycle
  do
  kdof=nactdof(nodempc(2,index),
  &                                        nodempc(1,index))
  if(kdof.ne.0) then
  do l=1,nev
  dbcont(l,jdofcont)=
  &                                      dbcont(l,jdofcont)+
  &                                   s(3*(j-1)+j1,3*(k-1)+k1)*
  &                                   coefmpc(index)*z(kdof,l)/
  &                                   coefmpc(ist)
  enddo
  endif
  index=nodempc(3,index)
  if(index.eq.0) exit
  enddo
  endif
  endif
  endif
  enddo
  enddo
  else
  jdof=8*(konl(j)-1)+j1
  call nident(ikmpc,jdof,nmpc,id1)
  if(id1.gt.0) then
  if(ikmpc(id1).eq.jdof) then
  id1=ilmpc(id1)
  ist1=ipompc(id1)
  index1=nodempc(3,ist1)
  if(index1.eq.0) cycle
  do
  jdof=nactdof(nodempc(2,index1),
  &                               nodempc(1,index1))
  if(jdof.ne.0) then
  call nident(iactcont,jdof,nactcont,id)
  jdofcont=0
  if(id.gt.0)then
  if(iactcont(id).eq.jdof) then
  jdofcont=id
  endif
  endif
  if(jdofcont.eq.0) then
  nactcont=nactcont+1
  do k=nactcont,id+2,-1
  iactcont(k)=iactcont(k-1)
  do l=1,nev
  dbcont(l,k)=dbcont(l,k-1)
  enddo
  enddo
  jdofcont=id+1
  iactcont(jdofcont)=jdof
  do l=1,nev
  dbcont(l,jdofcont)=0.d0
  enddo
  endif
  bcont(jdofcont)=bcont(jdofcont)+
  &                          coefmpc(index1)*
  &                          fnl(j1,j)/coefmpc(ist1)
  do k=1,nope
  do k1=1,3
  kdof=nactdof(k1,konl(k))
  if(kdof.ne.0) then
  do l=1,nev
  dbcont(l,jdofcont)=
  &                                      dbcont(l,jdofcont)
  &                                   +s(3*(j-1)+j1,3*(k-1)+k1)
  &                                   *coefmpc(index1)*z(kdof,l)/
  &                                   coefmpc(ist1)
  enddo
  else
  kdof=8*(konl(k)-1)+k1
  call nident(ikmpc,kdof,nmpc,id2)
  if(id2.gt.0) then
  if(ikmpc(id2).eq.kdof) then
  id2=ilmpc(id2)
  ist2=ipompc(id2)
  index2=nodempc(3,ist2)
  if(index2.eq.0) cycle
  do
  !
  !                   translated to the left to avoid exceedance
  !                   of 72 columns
  !
  kdof=nactdof(nodempc(2,index2),
  &                            nodempc(1,index2))
  if(kdof.ne.0) then
  do l=1,nev
  dbcont(l,jdofcont)=dbcont(l,jdofcont)
  &                                  -s(3*(j-1)+j1,3*(k-1)+k1)
  &                                  *coefmpc(index1)
  &                                  *coefmpc(index2)*z(kdof,l)/
  &                                  (coefmpc(ist1)*coefmpc(ist2))
  enddo
  endif
  index2=nodempc(3,index2)
  if(index2.eq.0) exit
  !
  !                   end of translation
  !
  enddo
  endif
  endif
  endif
  enddo
  enddo
  endif
  index1=nodempc(3,index1)
  if(index1.eq.0) exit
  enddo
  endif
  endif
  endif
  enddo
  enddo
  !
  return
  end
*/
