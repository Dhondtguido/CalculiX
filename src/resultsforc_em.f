!
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
      subroutine resultsforc_em(nk,f,fn,nactdof,ipompc,nodempc,
     &  coefmpc,labmpc,nmpc,mi,fmpc,calcul_fn,calcul_f,inomat)
!
!     calculating the equation system internal force vector
!     (one entry for each active degree of freedom)
!
!     This routine is based on resultsforc. Whereas resultsforc
!     applies to mechanical and thermal degrees of freedom,
!     resultsforcem is limited to thermal degrees of freedom.
!     (application for electromagnetic calculations)
!
      implicit none
!
      character*20 labmpc(*)
!
      integer mi(*),nactdof(0:mi(2),*),ipompc(*),nodempc(3,*),nk,i,j,
     &  nmpc,ist,ndir,node,index,calcul_fn,calcul_f,inomat(*)
!
      real*8 f(*),fn(0:mi(2),*),coefmpc(*),fmpc(*),forcempc
!
!     subtracting the mpc force (for each linear mpc there is one
!     force; the actual force in a node belonging to the mpc is
!     obtained by multiplying this force with the nodal coefficient.
!     The force has to be subtracted from f, since it does not
!     appear on the rhs of the equation system
!
      if(calcul_fn.eq.1)then
        do i=1,nmpc
            ist=ipompc(i)
            node=nodempc(1,ist)
            if(inomat(node).eq.0) cycle
            ndir=nodempc(2,ist)
c            if(ndir.ne.0) cycle
            forcempc=fn(ndir,node)/coefmpc(ist)
            fmpc(i)=forcempc
            fn(ndir,node)=0.d0
            index=nodempc(3,ist)
            if(index.eq.0) cycle
            do
               node=nodempc(1,index)
               ndir=nodempc(2,index)
               fn(ndir,node)=fn(ndir,node)-coefmpc(index)*forcempc
               index=nodempc(3,index)
               if(index.eq.0) exit
            enddo
         enddo
      endif
!
!     calculating the system force vector
!
      if(calcul_f.eq.1) then
         do i=1,nk
            if(inomat(i).eq.0) cycle
            if(nactdof(0,i).gt.0) then
               f(nactdof(0,i))=fn(0,i)
            endif
c            do j=0,mi(2)
c               if(nactdof(j,i).ne.0) then
c                  f(nactdof(j,i))=fn(j,i)
c               endif
c            enddo
         enddo
      endif
!
!     adding the mpc force again to fn
!
      if(calcul_fn.eq.1)then
         do i=1,nmpc
            ist=ipompc(i)
            node=nodempc(1,ist)
            if(inomat(node).eq.0) cycle
            ndir=nodempc(2,ist)
c            if(ndir.ne.0) cycle
            forcempc=fmpc(i)
            fn(ndir,node)=forcempc*coefmpc(ist)
            index=nodempc(3,ist)
!
!           nodes not belonging to the structure have to be
!           taken out
!
            if(labmpc(i)(1:7).eq.'MEANROT') then
               if(nodempc(3,nodempc(3,index)).eq.0) cycle
            elseif(labmpc(i)(1:10).eq.'PRETENSION') then
               if(nodempc(3,index).eq.0) cycle
            elseif(labmpc(i)(1:5).eq.'RIGID') then
               if(nodempc(3,nodempc(3,nodempc(3,nodempc(3,nodempc(3,inde
     &x))))).eq.0) cycle
            else
               if(index.eq.0) cycle
            endif
            do
               node=nodempc(1,index)
               ndir=nodempc(2,index)
               fn(ndir,node)=fn(ndir,node)+coefmpc(index)*forcempc
               index=nodempc(3,index)
               if(labmpc(i)(1:7).eq.'MEANROT') then
                  if(nodempc(3,nodempc(3,index)).eq.0) exit
               elseif(labmpc(i)(1:10).eq.'PRETENSION') then
                  if(nodempc(3,index).eq.0) exit
               elseif(labmpc(i)(1:5).eq.'RIGID') then
                  if(nodempc(3,nodempc(3,nodempc(3,nodempc(3,nodempc(3,i
     &ndex))))).eq.0) exit
               else
                  if(index.eq.0) exit
               endif
            enddo
         enddo
      endif
!
      return
      end
