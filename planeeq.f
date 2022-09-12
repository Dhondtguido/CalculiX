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
      subroutine planeeq(cotet,nodef,planfal)
!
!     calculate the equation of a plane through the points nodef(1),
!     nodef(2) and nodef(3). The equation of the plane is
!     planfal(1)*x+planfal(2)*y+planfal(3)*z+planfal(4)=0; the 
!     coordinates of the nodes are stored in cotet
!
!     the vector (planefal(1),planfal(2),planfal(3)) is orthogonal
!     to the plane and is normalized
!
      implicit none
!
      integer nodef(3)
      real*8 cotet(3,*),planfal(4),x1,y1,z1,x2,y2,z2,dd
!
      x1=cotet(1,nodef(1))-cotet(1,nodef(2))
      y1=cotet(2,nodef(1))-cotet(2,nodef(2))
      z1=cotet(3,nodef(1))-cotet(3,nodef(2))
!
      x2=cotet(1,nodef(1))-cotet(1,nodef(3))
      y2=cotet(2,nodef(1))-cotet(2,nodef(3))
      z2=cotet(3,nodef(1))-cotet(3,nodef(3))
!
      planfal(1)=y1*z2-y2*z1
      planfal(2)=x2*z1-x1*z2
      planfal(3)=x1*y2-x2*y1
!
      dd=dsqrt(planfal(1)*planfal(1)+planfal(2)*planfal(2)+
     &         planfal(3)*planfal(3))
!
      if(dd.lt.1.d-10) then
         planfal(1)=0.d0
         planfal(2)=0.d0
         planfal(3)=0.d0
         planfal(4)=0.d0
         return
      endif
!
      planfal(1)=planfal(1)/dd
      planfal(2)=planfal(2)/dd
      planfal(3)=planfal(3)/dd
!
      planfal(4)=-(planfal(1)*cotet(1,nodef(1))+
     &             planfal(2)*cotet(2,nodef(1))+
     &             planfal(3)*cotet(3,nodef(1)))
!
      return
      end

