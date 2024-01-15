!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2023 Guido Dhondt
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
      subroutine scalesen(dgdxglob,feasdir,nk,nodedesi,ndesi,objectset,
     &     iscaleflag,iobject,ne2d)
!
!     Scaling the sensitivities      
!
!     iscaleflag=1: length of the vector is scaled to 1 
!     iscaleflag=2: greatest vector value is scaled to 1
!     iscaleglag=3: sensitivities of the objective function are multiplied
!                   with -1 in case of a minimization task
!     iscaleflag=4: greatest vector value of the final feasible direction 
!                   sensitivity feasdir is scaled to 1
!
      implicit none
!
      character*81 objectset(5,*)
!     
      integer nk,nodedesi(*),i,j,k,m,ndesi,iobject,iscaleflag,node,ne2d
!     
      real*8 dgdxglob(2,nk,*),feasdir(2,*),dd,dd2
!     
      if(iscaleflag.eq.1) then
!
!       normalization over all design nodes of a filtered design
!       response
!
        if(objectset(5,iobject)(81:81).ne.'G') then
          dd=0.d0
          do i=1,ndesi
            node=nodedesi(i)
            dd=dd+dgdxglob(2,node,iobject)**2
          enddo
          dd=dsqrt(dd)
          do i=1,ndesi
            node=nodedesi(i)
            dgdxglob(2,node,iobject)=dgdxglob(2,node,iobject)/dd
            if(ne2d.ne.0) then
              dgdxglob(2,node+1,iobject)=dgdxglob(2,node,iobject)
              dgdxglob(2,node+2,iobject)=dgdxglob(2,node,iobject)
            endif
          enddo
        endif
      elseif(iscaleflag.eq.2) then
        if(objectset(5,iobject)(81:81).ne.'G') then
          dd=0.d0
          do i=1,ndesi
            node=nodedesi(i)
            dd=max(dd,abs(dgdxglob(2,node,iobject)))
          enddo
          do i=1,ndesi
            node=nodedesi(i)
            dgdxglob(2,node,iobject)=dgdxglob(2,node,iobject)/dd
          enddo
        endif
      elseif(iscaleflag.eq.3) then
        if(objectset(2,1)(17:19).eq.'MIN') then
          do i=1,ndesi
            node=nodedesi(i)
!           NEXT LINE IS NOT NEEDED ???????????????
c            dgdxglob(1,node,1)=-dgdxglob(1,node,1)
            dgdxglob(2,node,1)=-dgdxglob(2,node,1)
          enddo
        endif
      elseif(iscaleflag.eq.4) then
!
!       scaling such that the maximum over all design variables is one
!
        dd=0.d0
        dd2=0.d0
        do i=1,ndesi
          node=nodedesi(i)
          dd2=max(dd2,abs(feasdir(2,node)))
        enddo
        if(dd2.le.0.d0) then
          dd2=1.0
        endif
        if(ne2d.eq.0) then
          do i=1,ndesi
            node=nodedesi(i)
            feasdir(2,node)=feasdir(2,node)/dd2
          enddo
        else
!
!         for 2d-calculations: copy the results to the
!         other expanded nodes as well
!
          do i=1,ndesi
            node=nodedesi(i)
            feasdir(2,node)=feasdir(2,node)/dd2
            feasdir(1,node+1)=feasdir(1,node)
            feasdir(1,node+2)=feasdir(1,node)
            feasdir(2,node+1)=feasdir(2,node)
            feasdir(2,node+2)=feasdir(2,node)
          enddo
        endif
      endif
!     
      return        
      end
