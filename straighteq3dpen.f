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
      subroutine straighteq3dpen(col,straight,xnor,noeq)
!
!     calculate the equation of the planes through the
!     edges of a triangle and perpendicular to the representative
!     normal of each triangle edge together
!     with the plane of the triangle itself with 
!     (col(1,1),col(2,1),col(3,1)),(col(1,2),col(2,2),col(3,2)),
!     (col(1,3),col(2,3),col(3,3))
!     as vertices (nodet(1),nodet(2),nodet(3)). 
!     The equation of the plane through the edge 
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
!     edgenor(*,1) is the representative normal for the triangle edge going
!     through the nodes 2-3, edgenor(*,2) for the edge 3-1 and
!     edgenor(*,3) for the edge 1-2
!
!     noeq is the triangle nodes for which the opposite side should not get
!     an equation (all coefficients zero; this is needed for plane strain,
!     plane stress and axisymmetric elements
!
      implicit none
!
      integer i,noeq
!
      real*8 col(3,3),straight(16),p12(3),p23(3),p31(3),dd,scal1,
     &       scal2,xp12(3),xp23(3),xp31(3),xq12(3),xq23(3),xq31(3),
     &       xnor(3,3),dxp,dxq,edgenor(3,3),dd12,dd23,dd31
!
!     sides of the triangle
!
      do i=1,3
         p12(i)=col(i,2)-col(i,1)
         p23(i)=col(i,3)-col(i,2)
         p31(i)=col(i,1)-col(i,3)
      enddo
      dd12=p12(1)*p12(1)+p12(2)*p12(2)+p12(3)*p12(3)
      dd23=p23(1)*p23(1)+p23(2)*p23(2)+p23(3)*p23(3)
      dd31=p31(1)*p31(1)+p31(2)*p31(2)+p31(3)*p31(3)
!
!     calculating the representative normal for each triangle edge
!     edgenor(*,1) for p23 | edgenor(*,2) for p31 | edgenor(*,3) for p12
!  
         scal1=(xnor(1,2)*p23(1)+xnor(2,2)*p23(2)+xnor(3,2)*p23(3))/dd23
         scal2=(xnor(1,3)*p23(1)+xnor(2,3)*p23(2)+xnor(3,3)*p23(3))/dd23
         do i=1,3
            xp23(i)=xnor(i,2)-scal1*p23(i)
            xq23(i)=xnor(i,3)-scal2*p23(i)
         enddo
         dxp=dsqrt(xp23(1)*xp23(1)+xp23(2)*xp23(2)+xp23(3)*xp23(3))
         dxq=dsqrt(xq23(1)*xq23(1)+xq23(2)*xq23(2)+xq23(3)*xq23(3))
         do i=1,3
            xp23(i)=xp23(i)/dxp
            xq23(i)=xq23(i)/dxq
         enddo
!        
       scal1=(xnor(1,3)*p31(1)+xnor(2,3)*p31(2)+xnor(3,3)*p31(3))/dd31
       scal2=(xnor(1,1)*p31(1)+xnor(2,1)*p31(2)+xnor(3,1)*p31(3))/dd31
       do i=1,3
          xp31(i)=xnor(i,3)-scal1*p31(i)
          xq31(i)=xnor(i,1)-scal2*p31(i)
       enddo
       dxp=dsqrt(xp31(1)*xp31(1)+xp31(2)*xp31(2)+xp31(3)*xp31(3))
       dxq=dsqrt(xq31(1)*xq31(1)+xq31(2)*xq31(2)+xq31(3)*xq31(3))
       do i=1,3
          xp31(i)=xp31(i)/dxp
          xq31(i)=xq31(i)/dxq
       enddo
!
       scal1=(xnor(1,1)*p12(1)+xnor(2,1)*p12(2)+xnor(3,1)*p12(3))/dd12
       scal2=(xnor(1,2)*p12(1)+xnor(2,2)*p12(2)+xnor(3,2)*p12(3))/dd12
       do i=1,3
          xp12(i)=xnor(i,1)-scal1*p12(i)
          xq12(i)=xnor(i,2)-scal2*p12(i)
       enddo
       dxp=dsqrt(xp12(1)*xp12(1)+xp12(2)*xp12(2)+xp12(3)*xp12(3))
       dxq=dsqrt(xq12(1)*xq12(1)+xq12(2)*xq12(2)+xq12(3)*xq12(3))
       do i=1,3
          xp12(i)=xp12(i)/dxp
          xq12(i)=xq12(i)/dxq
       enddo 
!    
       do i=1,3
          edgenor(i,1)=xp23(i)+xq23(i)
          edgenor(i,2)=xp31(i)+xq31(i)
          edgenor(i,3)=xp12(i)+xq12(i)
       enddo
!     
!     normalized vector normal to each side of the triangle
! 
       dd=dsqrt(edgenor(1,1)*edgenor(1,1)+edgenor(2,1)*edgenor(2,1)
     &      +edgenor(3,1)*edgenor(3,1))
       do i=1,3
          edgenor(i,1)=edgenor(i,1)/dd
       enddo
       dd=dsqrt(edgenor(1,2)*edgenor(1,2)+edgenor(2,2)*edgenor(2,2)
     &      +edgenor(3,2)*edgenor(3,2))
       do i=1,3
          edgenor(i,2)=edgenor(i,2)/dd
       enddo
       dd=dsqrt(edgenor(1,3)*edgenor(1,3)+edgenor(2,3)*edgenor(2,3)
     &      +edgenor(3,3)*edgenor(3,3))
       do i=1,3
          edgenor(i,3)=edgenor(i,3)/dd
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
!     p12 x edgenor(*,3)
!
      if(noeq.ne.3) then
         straight(9)=p12(2)*edgenor(3,3)-p12(3)*edgenor(2,3)
         straight(10)=p12(3)*edgenor(1,3)-p12(1)*edgenor(3,3)
         straight(11)=p12(1)*edgenor(2,3)-p12(2)*edgenor(1,3)
         dd=dsqrt(straight(9)*straight(9)+straight(10)*straight(10)+
     &        straight(11)*straight(11))
         do i=9,11
            straight(i)=straight(i)/dd
         enddo
      else
         do i=9,11
            straight(i)=0.d0
         enddo
      endif
!
!     p23 x edgenor(*,1)
!
      if(noeq.ne.1) then
         straight(1)=p23(2)*edgenor(3,1)-p23(3)*edgenor(2,1)
         straight(2)=p23(3)*edgenor(1,1)-p23(1)*edgenor(3,1)
         straight(3)=p23(1)*edgenor(2,1)-p23(2)*edgenor(1,1)
         dd=dsqrt(straight(1)*straight(1)+straight(2)*straight(2)+
     &        straight(3)*straight(3))
         do i=1,3
            straight(i)=straight(i)/dd
         enddo
      else
         do i=1,3
            straight(i)=0.d0
         enddo
      endif
!
!     p31 x edgenor(*,2)
!
      if(noeq.ne.2) then
         straight(5)=p31(2)*edgenor(3,2)-p31(3)*edgenor(2,2)
         straight(6)=p31(3)*edgenor(1,2)-p31(1)*edgenor(3,2)
         straight(7)=p31(1)*edgenor(2,2)-p31(2)*edgenor(1,2)
         dd=dsqrt(straight(5)*straight(5)+straight(6)*straight(6)+
     &        straight(7)*straight(7))
         do i=5,7
            straight(i)=straight(i)/dd
         enddo
      else
         do i=5,7
            straight(i)=0.d0
         enddo
      endif
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

