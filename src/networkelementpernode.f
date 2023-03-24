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
      subroutine networkelementpernode(iponoel,inoel,lakon,ipkon,kon,
     &       inoelsize,nflow,ieg,ne,network)
!
      implicit none
!
      character*8 lakon(*)
!
      integer iponoel(*),inoel(2,*),ipkon(*),kon(*),i,j,k,
     &  inoelfree,nope,indexe,node,inoelsize,nflow,ieg(*),ne,
     &  network
!
!     determining the elements belonging to the nodes of
!     the elements
!
!     network<=1: simultaneous procedure
!     network>1: alternating procedure
!
      inoelfree=1
!
      if(network.gt.1) then
         do k=1,nflow
            i=ieg(k)
            indexe=ipkon(i)
            do j=1,3
               node=kon(indexe+j)
               if(node.eq.0) cycle
               inoel(1,inoelfree)=i
               inoel(2,inoelfree)=iponoel(node)
               iponoel(node)=inoelfree
               inoelfree=inoelfree+1
            enddo
         enddo
      else
         do i=1,ne
            if(lakon(i)(1:1).eq.'D') then
               indexe=ipkon(i)
               do j=1,3
                  node=kon(indexe+j)
                  if(node.eq.0) cycle
                  inoel(1,inoelfree)=i
                  inoel(2,inoelfree)=iponoel(node)
                  iponoel(node)=inoelfree
                  inoelfree=inoelfree+1
               enddo
            endif
         enddo
      endif
!
!     size of field inoel
!
      inoelsize=inoelfree-1
!
      return
      end
