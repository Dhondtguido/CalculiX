!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998 Guido Dhondt
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
      subroutine straighteq3d(col,straight)
!
!     calculate the equation of the planes through the
!     edges of a triangle and perpendicular to the triangle together
!     with the plane of the triangle itself with 
!     (col(1,1),col(2,1),col(3,1)),(col(1,2),col(2,2),col(3,2)),
!     (col(1,3),col(2,3),col(3,3))
!     as vertices. The equation of the plane through the edge 
!     opposite nodet(1) is of the form
!     straight(1)*x+straight(2)*y+straight(3)*z+straight(4)=0, such that the
!     vector (straight(1),straight(2),straight(3)) points outwards; 
!     for the edge opposite of nodet(2) the equation is 
!     straight(5)*x+straight(6)*y+straight(7)*z+straight(8)=0 and for the edge
!     oppositie of nodet(3) it is
!     straight(9)*x+straight(10)*y+straight(11)*z+straight(12)=0. 
!     Here too, the normals
!     (straight(5),straight(6),straight(7)) and 
!     (straight(9),straight(10),straight(11)) point
!     outwards of the triangle. The equation of the triangle plane is
!     straight(13)*x+straight(14)*y+straight(15)*z+straight(16)=0 such
!     that the triangle is numbered clockwise when looking in the 
!     direction of vector (straight(13),straight(14),straight(15)).
!
      implicit none
!
      integer i
!
      real*8 col(3,3),straight(16),p12(3),p23(3),p31(3),dd
!
!     sides of the triangle
!
      do i=1,3
         p12(i)=col(i,2)-col(i,1)
         p23(i)=col(i,3)-col(i,2)
         p31(i)=col(i,1)-col(i,3)
      enddo
!
!     normalized vector normal to the triangle: xn = p12 x p23
!
      straight(13)=p12(2)*p23(3)-p12(3)*p23(2)
      straight(14)=p12(3)*p23(1)-p12(1)*p23(3)
      straight(15)=p12(1)*p23(2)-p12(2)*p23(1)
      dd=dsqrt(straight(13)*straight(13)+straight(14)*straight(14)+
     &         straight(15)*straight(15))
      do i=13,15
         straight(i)=straight(i)/dd
      enddo
!
!     p12 x xn
!
      straight(9)=p12(2)*straight(15)-p12(3)*straight(14)
      straight(10)=p12(3)*straight(13)-p12(1)*straight(15)
      straight(11)=p12(1)*straight(14)-p12(2)*straight(13)
      dd=dsqrt(straight(9)*straight(9)+straight(10)*straight(10)+
     &         straight(11)*straight(11))
      do i=9,11
         straight(i)=straight(i)/dd
      enddo
!
!     p23 x xn
!
      straight(1)=p23(2)*straight(15)-p23(3)*straight(14)
      straight(2)=p23(3)*straight(13)-p23(1)*straight(15)
      straight(3)=p23(1)*straight(14)-p23(2)*straight(13)
      dd=dsqrt(straight(1)*straight(1)+straight(2)*straight(2)+
     &         straight(3)*straight(3))
      do i=1,3
         straight(i)=straight(i)/dd
      enddo
!
!     p31 x xn
!
      straight(5)=p31(2)*straight(15)-p31(3)*straight(14)
      straight(6)=p31(3)*straight(13)-p31(1)*straight(15)
      straight(7)=p31(1)*straight(14)-p31(2)*straight(13)
      dd=dsqrt(straight(5)*straight(5)+straight(6)*straight(6)+
     &         straight(7)*straight(7))
      do i=5,7
         straight(i)=straight(i)/dd
      enddo
!
!     determining the inhomogeneous terms
!
      straight(12)=-straight(9)*col(1,1)-straight(10)*col(2,1)-
     &             straight(11)*col(3,1)
      straight(4)=-straight(1)*col(1,2)-straight(2)*col(2,2)-
     &             straight(3)*col(3,2)
      straight(8)=-straight(5)*col(1,3)-straight(6)*col(2,3)-
     &             straight(7)*col(3,3)
      straight(16)=-straight(13)*col(1,1)-straight(14)*col(2,1)-
     &             straight(15)*col(3,1)
!
      return
      end

