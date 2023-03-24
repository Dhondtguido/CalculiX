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
      subroutine straighteq2d(col,straight)
!
!     calculate the equation of the edges of a triangle with 
!     (col(1,1),col(2,1)),(col(1,2),col(2,2)),(col(1,3),col(2,3))
!     as vertices. The equation of the edge 
!     opposite notet(1) is of the form
!     straight(1)*x+straight(2)*y+straight(3)=0, such that the
!     vector (straight(1),straight(2)) points outwards; for the edge
!     opposite of nodet(2) the equation is 
!     straight(4)*x+straight(5)*y+straight(6)=0 and for the edge
!     oppositie of nodet(3) it is
!     straight(7)*x+straight(8)*y+straight(9)=0. Here too, the normals
!     (straight(4),straight(5)) and (straight(7),straight(8)) point
!     outwards of the triangle.
!
      implicit none
!
      real*8 col(2,3),straight(9),x1,y1,dd
!
!     edge opposite of 1
!
      x1=col(1,3)-col(1,2)
      y1=col(2,3)-col(2,2)
      dd=dsqrt(x1*x1+y1*y1)
!
      straight(1)=y1/dd
      straight(2)=-x1/dd
!
      straight(3)=-(straight(1)*col(1,3)+
     &             straight(2)*col(2,3))
!
      if(straight(1)*col(1,1)+straight(2)*col(2,1)+
     &   straight(3).gt.0.d0) then
         straight(1)=-straight(1)
         straight(2)=-straight(2)
         straight(3)=-straight(3)
      endif
!
!     edge opposite of 2
!
      x1=col(1,1)-col(1,3)
      y1=col(2,1)-col(2,3)
      dd=dsqrt(x1*x1+y1*y1)
!
      straight(4)=y1/dd
      straight(5)=-x1/dd
!
      straight(6)=-(straight(4)*col(1,1)+
     &             straight(5)*col(2,1))
!
      if(straight(4)*col(1,2)+straight(5)*col(2,2)+
     &   straight(6).gt.0.d0) then
         straight(4)=-straight(4)
         straight(5)=-straight(5)
         straight(6)=-straight(6)
      endif
!
!     edge opposite of 3
!
      x1=col(1,2)-col(1,1)
      y1=col(2,2)-col(2,1)
      dd=dsqrt(x1*x1+y1*y1)
!
      straight(7)=y1/dd
      straight(8)=-x1/dd
!
      straight(9)=-(straight(7)*col(1,2)+
     &             straight(8)*col(2,2))
!
      if(straight(7)*col(1,3)+straight(8)*col(2,3)+
     &   straight(9).gt.0.d0) then
         straight(7)=-straight(7)
         straight(8)=-straight(8)
         straight(9)=-straight(9)
      endif
!
      return
      end

