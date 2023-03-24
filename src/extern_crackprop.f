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
      subroutine extern_crackprop(ieled,nedg,ibounedg,nbounedg,
     &     ibounnod,nbounnod,iedg,iedno,ier)
!
!     determines which edges and nodes are external
!     
!     an edge is external if it belongs to only one element   
!     a node is external if it belongs to an external edge
!
      implicit none
!
      integer ieled(2,*),nedg,ibounedg(*),i,j,k,nbounedg,ier,
     &     ibounnod(*),nbounnod,id,node,iedg(3,*),iedno(2,*)
!
!     store the "external" edges of the crack shapes
!
      nbounedg=0
      do i=1,nedg
         if(ieled(2,i).eq.0) then
            nbounedg=nbounedg+1
            ibounedg(nbounedg)=i
         endif
      enddo
!
!     store the "external" nodes of the crack shapes
!
      nbounnod=0
      do i=1,nbounedg
         do j=1,2
            node=iedg(j,ibounedg(i))
            call nident(ibounnod,node,nbounnod,id)
            if(id.gt.0) then
               if(ibounnod(id).eq.node) cycle
            endif
            nbounnod=nbounnod+1
            do k=nbounnod,id+2,-1
               ibounnod(k)=ibounnod(k-1)
            enddo
            ibounnod(id+1)=node
         enddo
      enddo
!
!     store the external edges to which an external node belongs
!     in iedno
!
      do i=1,nbounedg
         do j=1,2
            node=iedg(j,ibounedg(i))
            call nident(ibounnod,node,nbounnod,id)
            if(iedno(1,id).eq.0) then
               iedno(1,id)=i
            elseif(iedno(2,id).eq.0) then
              iedno(2,id)=i
            else
              write(*,*) '*ERROR in extern_crackprop: a node on the'
              write(*,*) '       boundary of a crack belongs to more'
              write(*,*) '       than two external edges: crack mesh'
              write(*,*) '       seems to be corrupt.'
              ier=1
              return
            endif
         enddo
      enddo
!
      return
      end
