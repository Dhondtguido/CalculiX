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

void resultsforc(ITG *nk,double *f,double *fn,ITG *nactdof,ITG *ipompc,
		 ITG *nodempc,double *coefmpc,char *labmpc,ITG *nmpc,
		 ITG *mi,double *fmpc,ITG *calcul_fn,ITG *calcul_f,
                 ITG *num_cpus){

    ITG i,j,ist,node,ndir,mt=mi[1]+1,index,index2;

    double forcempc;
    
/*     subtracting the mpc force (for each linear mpc there is one
       force; the actual force in a node belonging to the mpc is
       obtained by multiplying this force with the nodal coefficient.
       The force has to be subtracted from f, since it does not
       appear on the rhs of the equation system */
    
    if(*calcul_fn==1){
	for(i=0;i<*nmpc;i++){
	    if(strcmp1(&labmpc[20*i],"FLUID")==0) continue;
	    ist=ipompc[i]-1;
	    node=nodempc[3*ist]-1;
	    ndir=nodempc[3*ist+1];
	    if(ndir>3) continue;
	    forcempc=fn[mt*node+ndir]/coefmpc[ist];
	    fmpc[i]=forcempc;
	    fn[mt*node+ndir]=0.;
	    index=nodempc[3*ist+2]-1;
	    if(index==-1) continue;
	    do{
		node=nodempc[3*index]-1;
		ndir=nodempc[3*index+1];
		fn[mt*node+ndir]-=coefmpc[index]*forcempc;
		index=nodempc[3*index+2]-1;
		if(index==-1) break;
	    }while(1);
	}
    }
    
//     calculating the system force vector

    if(*calcul_f==1){
	forparll(&mt,nactdof,f,fn,nk,num_cpus);
/*	for(i=0;i<*nk;i++){
	    for(j=0;j<mt;j++){
		if(nactdof[mt*i+j]>0){
		    f[nactdof[mt*i+j]-1]=fn[mt*i+j];
		}
	    }
	    }*/
    }
    
    /* adding the mpc force again to fn */

    if(*calcul_fn==1){
	for(i=0;i<*nmpc;i++){
	    if(strcmp1(&labmpc[20*i],"FLUID")==0) continue;
	    ist=ipompc[i]-1;
	    node=nodempc[3*ist]-1;
	    ndir=nodempc[3*ist+1];
	    if(ndir>3) continue;
	    forcempc=fmpc[i];
	    fn[mt*node+ndir]=forcempc*coefmpc[ist];
	    index=nodempc[3*ist+2]-1;

            /* nodes not belonging to the structure have to be taken out */
	    
	    if(strcmp1(&labmpc[20*i],"MEANROT")==0){
		index2=nodempc[3*index+2]-1;
		index2=nodempc[3*index2+2];
		if(index2==0) continue;
	    }else if(strcmp1(&labmpc[20*i],"PRETENSION")==0){
		if(nodempc[3*index+2]==0) continue;
	    }else if(strcmp1(&labmpc[20*i],"RIGID")==0){
		index2=nodempc[3*index+2]-1;
		index2=nodempc[3*index2+2]-1;
		index2=nodempc[3*index2+2]-1;
		index2=nodempc[3*index2+2]-1;
		index2=nodempc[3*index2+2];
		if(index2==0) continue;
	    }else{
		if(index==-1) continue;
	    }

	    do{
		node=nodempc[3*index]-1;
		ndir=nodempc[3*index+1];
		fn[mt*node+ndir]+=coefmpc[index]*forcempc;
		index=nodempc[3*index+2]-1;

                /* nodes not belonging to the structure have to be taken out */

		if(strcmp1(&labmpc[20*i],"MEANROT")==0){
		    index2=nodempc[3*index+2]-1;
		    index2=nodempc[3*index2+2];
		    if(index2==0) break;
		}else if(strcmp1(&labmpc[20*i],"PRETENSION")==0){
		    if(nodempc[3*index+2]==0) break;
		}else if(strcmp1(&labmpc[20*i],"RIGID")==0){
		    index2=nodempc[3*index+2]-1;
		    index2=nodempc[3*index2+2]-1;
		    index2=nodempc[3*index2+2]-1;
		    index2=nodempc[3*index2+2]-1;
		    index2=nodempc[3*index2+2];
		    if(index2==0) break;
		}else{
		    if(index==-1) break;
		}

	    }while(1);
	}
    }
    
    return;

}
/*c!
c!     CalculiX - A 3-dimensional finite element program
c!              Copyright (C) 1998-2015 Guido Dhondt
c!
c!     This program is free software; you can redistribute it and/or
c!     modify it under the terms of the GNU General Public License as
c!     published by the Free Software Foundation(version 2);
c!     
c!
c!     This program is distributed in the hope that it will be useful,
c!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
c!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
c!     GNU General Public License for more details.
c!
c!     You should have received a copy of the GNU General Public License
c!     along with this program; if not, write to the Free Software
c!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
c!
c      subroutine resultsforc(nk,f,fn,nactdof,ipompc,nodempc,
c     &  coefmpc,labmpc,nmpc,mi,fmpc,calcul_fn,calcul_f)
c!
c!     calculating the equation system internal force vector
c!     (one entry for each active degree of freedom)
c!
c      implicit none
c!
c      character*20 labmpc(*)
c!
c      integer mi(*),nactdof(0:mi(2),*),ipompc(*),nodempc(3,*),nk,i,j,
c     &  nmpc,ist,ndir,node,index,calcul_fn,calcul_f
c!
c      real*8 f(*),fn(0:mi(2),*),coefmpc(*),fmpc(*),forcempc
c!
c!     subtracting the mpc force (for each linear mpc there is one
c!     force; the actual force in a node belonging to the mpc is
c!     obtained by multiplying this force with the nodal coefficient.
c!     The force has to be subtracted from f, since it does not
c!     appear on the rhs of the equation system
c!
c      if(calcul_fn.eq.1)then
c        do i=1,nmpc
c            if(labmpc(i)(1:5).eq.'FLUID') cycle
c            ist=ipompc(i)
c            node=nodempc(1,ist)
c            ndir=nodempc(2,ist)
c            if(ndir.gt.3) cycle
c            forcempc=fn(ndir,node)/coefmpc(ist)
c            fmpc(i)=forcempc
c            fn(ndir,node)=0.d0
c            index=nodempc(3,ist)
c            if(index.eq.0) cycle
c            do
c               node=nodempc(1,index)
c               ndir=nodempc(2,index)
c               fn(ndir,node)=fn(ndir,node)-coefmpc(index)*forcempc
c               index=nodempc(3,index)
c               if(index.eq.0) exit
c            enddo
c         enddo
c      endif
c!
c!     calculating the system force vector
c!
c      if(calcul_f.eq.1) then
c         do i=1,nk
c            do j=0,mi(2)
c               if(nactdof(j,i).gt.0) then
c                  f(nactdof(j,i))=fn(j,i)
c               endif
c            enddo
c         enddo
c      endif
c!
c!     adding the mpc force again to fn
c!
c      if(calcul_fn.eq.1)then
c         do i=1,nmpc
c            if(labmpc(i)(1:5).eq.'FLUID') cycle
c            ist=ipompc(i)
c            node=nodempc(1,ist)
c            ndir=nodempc(2,ist)
c            if(ndir.gt.3) cycle
c            forcempc=fmpc(i)
c            fn(ndir,node)=forcempc*coefmpc(ist)
c            index=nodempc(3,ist)
c!
c!           nodes not belonging to the structure have to be
c!           taken out
c!
c            if(labmpc(i)(1:7).eq.'MEANROT') then
c               if(nodempc(3,nodempc(3,index)).eq.0) cycle
c            elseif(labmpc(i)(1:10).eq.'PRETENSION') then
c               if(nodempc(3,index).eq.0) cycle
c            elseif(labmpc(i)(1:5).eq.'RIGID') then
c               if(nodempc(3,nodempc(3,nodempc(3,nodempc(3,nodempc(3,inde
c     &x))))).eq.0) cycle
c            else
c               if(index.eq.0) cycle
c            endif
c            do
c               node=nodempc(1,index)
c               ndir=nodempc(2,index)
c               fn(ndir,node)=fn(ndir,node)+coefmpc(index)*forcempc
c               index=nodempc(3,index)
c               if(labmpc(i)(1:7).eq.'MEANROT') then
c                  if(nodempc(3,nodempc(3,index)).eq.0) exit
c               elseif(labmpc(i)(1:10).eq.'PRETENSION') then
c                  if(nodempc(3,index).eq.0) exit
c               elseif(labmpc(i)(1:5).eq.'RIGID') then
c                  if(nodempc(3,nodempc(3,nodempc(3,nodempc(3,nodempc(3,i
c     &ndex))))).eq.0) exit
c               else
c                  if(index.eq.0) exit
c               endif
c            enddo
c         enddo
c      endif
c!
c      return
c      end*/
