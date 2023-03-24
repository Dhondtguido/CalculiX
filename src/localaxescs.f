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
      subroutine localaxescs(cs,mcs,e1,e2,xn)
!
!     determines a local axis system based on the rotation axis
!     defined on a CYCLIC SYMMETRY MODEL card 
!
      implicit none
!
      integer mcs,i,imax
!
      real*8 cs(17,*),xn(3),e1(3),e2(3),xmax,dd
!
!     xn: axis direction; first cyclic symmetry definition is taken
!
      do i=1,3
         xn(1)=cs(9,1)-cs(6,1)
         xn(2)=cs(10,1)-cs(7,1)
         xn(3)=cs(11,1)-cs(8,1)
      enddo
      dd=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
      do i=1,3
         xn(i)=xn(i)/dd
      enddo
!
!     e1: unit vector orthogonal to xn
!
      if(xn(1).eq.0.d0) then
         e1(1)=1.d0
         e1(2)=0.d0
         e1(3)=0.d0
      elseif(xn(2).eq.0.d0) then
         e1(1)=0.d0
         e1(2)=1.d0
         e1(3)=0.d0
      elseif(xn(3).eq.0.d0) then
         e1(1)=0.d0
         e1(2)=0.d0
         e1(3)=1.d0
      else
!
!        determining the maximum entry in xn in absolute value
!
         xmax=0.d0
         if(dabs(xn(1)).gt.xmax) then
            xmax=dabs(xn(1))
            imax=1
         endif
         if(dabs(xn(2)).gt.xmax) then
            xmax=dabs(xn(2))
            imax=2
         endif
         if(dabs(xn(3)).gt.xmax) then
            xmax=dabs(xn(3))
            imax=3
         endif
!
!        creating a vector orthogonal to xn using the maximum
!        component value of xn
!
         e1(1)=1.d0
         e1(2)=1.d0
         e1(3)=1.d0
!
         e1(imax)=-(xn(1)+xn(2)+xn(3)-xn(imax))/xn(imax)
!
!        normalizing e1
!
         dd=dsqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
         do i=1,3
            e1(i)=e1(i)/dd
         enddo
      endif
!
!     e2 = n x e1
!
      e2(1)=xn(2)*e1(3)-xn(3)*e1(2)
      e2(2)=xn(3)*e1(1)-xn(1)*e1(3)
      e2(3)=xn(1)*e1(2)-xn(2)*e1(1)
!
c      write(*,*) 'localaxes',e1(1),e1(2),e1(3)
c      write(*,*) 'localaxes',e2(1),e2(2),e2(3)
c      write(*,*) 'localaxes',xn(1),xn(2),xn(3)
!         
      return
      end

