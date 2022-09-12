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
      subroutine addimdnodecload(nodeforc,iforc,imdnode,nmdnode,xforc,
     &              ikmpc,ilmpc,ipompc,
     &              nodempc,nmpc,imddof,nmddof,
     &              nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &              ikboun,nboun,ilboun,ithermal)
!
!     adds the dof in which a user-defined point force was applied to imdnode
!     (needed in dyna.c and steadystate.c)
!
      implicit none
!
      integer nodeforc(2,*),iforc,node,imdnode(*),nmdnode,ikmpc(*),
     &  ilmpc(*),ipompc(*),nodempc(3,*),nmpc,imddof(*),nmddof,
     &  mi(*),nactdof(0:mi(2),*),imdmpc(*),nmdmpc,imdboun(*),nmdboun,
     &  ikboun(*),nboun,ilboun(*),ithermal(*),k
!
      real*8 xforc(*)
!
      node=nodeforc(1,iforc)
!
!     user-defined load
!
      if((xforc(iforc).lt.1.2357111318d0).and.
     &     (xforc(iforc).gt.1.2357111316d0)) then
!
         call addimd(imdnode,nmdnode,node)
!
!        add the degrees of freedom corresponding to the node
!
         if(ithermal(1).ne.2) then
            do k=1,3
               call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &              nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &              nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &              ikboun,nboun,ilboun)
            enddo
         else
            k=0
            call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &           nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &           nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,ikboun,
     &           nboun,ilboun)
         endif
      endif
!
      return
      end

