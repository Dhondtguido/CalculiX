!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine networkneighbor(nelem,node,nelemnei,nodenei,ibranch,
     &  iponoel,inoel,ipkon,kon)
!
!     looks for the neighboring element and neighboring end node of
!     node "node" of element "nelem". If the neighboring end node
!     belongs to more than 2 elements ibranch=1, else ibranch=0
!     
      implicit none
!
      integer nelem,node,nelemnei,nodenei,ibranch,index,indexe,
     &  iponoel(*),inoel(2,*),ipkon(*),kon(*),iel
!
      nelemnei=0
      ibranch=0
!
      index=iponoel(node)
      if(index.eq.0) then
         write(*,*) '*ERROR in networkneighbor:node',node
         write(*,*) '       does not belong to network element',nelem
         call exit(201)
      endif
!
      do
         iel=inoel(1,index)
!
         if(iel.eq.nelem) then
            index=inoel(2,index)
            if(index.eq.0) exit
            cycle
         endif
!
!        neighboring element; check whether a neighboring element
!        was already found
!
         if(nelemnei.ne.0) then
            ibranch=1
            exit
         endif
!
         nelemnei=iel
         indexe=ipkon(iel)
         if(kon(indexe+1).eq.node) then
            nodenei=kon(indexe+3)
         else
            nodenei=kon(indexe+1)
         endif
         index=inoel(2,index)
         if(index.eq.0) exit
      enddo
!
      return
      end
