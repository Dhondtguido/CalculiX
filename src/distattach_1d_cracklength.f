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
      subroutine distattach_1d_cracklength(xig,pneigh,a,p,ratio,
     &     xo,yo,zo,x,y,z,nx,ny,nz,n)
!
!     calculates the minimum distance from a node with local
!     coordinate xig on a line described by 2 nodes with coordinates
!     in pneigh from the set of all n front nodes belonging to a
!     given crack (stored in fields xo,yo,zo,x,y,z,nx,ny,nz).
!
      implicit none
!
      integer i,j,nx(*),ny(*),nz(*),n,neigh(1),kneigh
!
      real*8 ratio(*),pneigh(3,*),a,xig,p(3),xo(*),yo(*),zo(*),
     &     x(*),y(*),z(*)
!
      ratio(1)=(1.d0-xig)/2.d0
      ratio(2)=(1.d0+xig)/2.d0
!
!     calculating the position in the face
!      
      do i=1,3
        p(i)=ratio(1)*pneigh(i,1)+ratio(2)*pneigh(i,2)
      enddo
!
!     calculating the smallest distance from the front nodes
!     belonging to the crack
!
      kneigh=1
      call near3d(xo,yo,zo,x,y,z,nx,ny,nz,p(1),p(2),p(3),
     &     n,neigh,kneigh)
!
!     calculating the negative distance (maximum is sought, not
!     minimum)
!
      a=-((xo(neigh(1))-p(1))**2+(yo(neigh(1))-p(2))**2
     &     +(zo(neigh(1))-p(3))**2)
!
      return
      end
      
