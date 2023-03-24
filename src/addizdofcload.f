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
      subroutine addizdofcload(nodeforc,ndirforc,nactdof,mi,izdof,
     &  nzdof,iforc,iznode,nznode,nk,imdnode,nmdnode,xforc,ntrans,
     &  inotr)
!
!     adds the dof in which a point force was applied to iznode, izdof
!     and to ** imdnode if user-defined load **
!     (needed in dyna.c and steadystate.c)
!
      implicit none
!
      integer nodeforc(2,*),ndirforc(*),iforc,node,j,jdof,mi(*),nk,
     &  nactdof(0:mi(2),*),izdof(*),nzdof,iznode(*),nznode,nodebasis,
     &  imdnode(*),nmdnode,ntrans,itr,inotr(2,*)
!
      real*8 xforc(*)
!
      node=nodeforc(1,iforc)
!
!     adding the nodes in the basis sector to iznode
!
      nodebasis=mod(node,nk)
      call addimd(iznode,nznode,nodebasis)
c!
c!     user-defined load
c!
c      if((xforc(iforc).lt.1.2357111318d0).and.
c     &     (xforc(iforc).gt.1.2357111316d0)) then
c         call addimd(imdnode,nmdnode,node)
c      endif
!
      if(ntrans.eq.0) then
         itr=0
      else
         itr=inotr(1,node)
      endif
!
      if(itr.eq.0) then
!        
!        no local transformation
!
         j=ndirforc(iforc)
!
!        C-convention!
!
         jdof=nactdof(j,node)-1
         if(jdof.gt.0) call addimd(izdof,nzdof,jdof)
      else
!
!        local transformation: loop over all dofs
!
         do j=1,3
            jdof=nactdof(j,node)-1
            if(jdof.gt.0) call addimd(izdof,nzdof,jdof)
         enddo
      endif
!
      return
      end

