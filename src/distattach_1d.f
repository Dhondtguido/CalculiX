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
      subroutine distattach_1d(xig,pneigh,pnode,a,p,ratio,nterms)
!
!     calculates the distance between the node with coordinates
!     in "pnode" and the node with local coordinate xig 
!     on a line described by "nterms" nodes with coordinates
!     in pneigh
!
      implicit none
!
      integer nterms,i,j
!
      real*8 ratio(3),pneigh(3,*),pnode(3),a,xig,p(3)
!
!
!
      if(nterms.eq.2) then
         ratio(1)=(1.d0-xig)/2.d0
         ratio(2)=(1.d0+xig)/2.d0
      elseif(nterms.eq.3) then
         ratio(1)=xig*(xig-1.d0)/2.d0
         ratio(2)=(1.d0-xig)*(1.d0+xig)
         ratio(3)=xig*(xig+1.d0)/2.d0
      else
         write(*,*) '*ERROR in distattach1d: case with ',nterms
         write(*,*) '       terms is not covered'
         call exit(201)
      endif
!
!     calculating the position in the face
!      
      do i=1,3
         p(i)=0.d0
         do j=1,nterms
            p(i)=p(i)+ratio(j)*pneigh(i,j)
         enddo
      enddo
!
!     calculating the distance
!
      a=(pnode(1)-p(1))**2+(pnode(2)-p(2))**2+(pnode(3)-p(3))**2
!
      return
      end
      
