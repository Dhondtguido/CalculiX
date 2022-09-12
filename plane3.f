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
      subroutine plane3(co,nodep,a,b,c,d)
!
!     calculate the equation of the plane through the
!     nodes nodep(1),nodep(2) and nodep(3) in the form
!     a*x+b*y+c*z+d=0 such that the triangle through the
!     nodes nodep(1),nodep(2),nopep(3) is numbered clockwise 
!     when looking in the direction of vector (a,b,c)
!
      implicit none
!
      integer nodep(3),i
!
      real*8 co(3,*),a,b,c,d,dd,p12(3),p23(3),p31(3)
!
!     sides of the triangle
!
      do i=1,3
         p12(i)=co(i,nodep(2))-co(i,nodep(1))
         p23(i)=co(i,nodep(3))-co(i,nodep(2))
         p31(i)=co(i,nodep(1))-co(i,nodep(3))
      enddo
!
!     normalized vector normal to the triangle: xn = p12 x p23
!
      a=p12(2)*p23(3)-p12(3)*p23(2)
      b=p12(3)*p23(1)-p12(1)*p23(3)
      c=p12(1)*p23(2)-p12(2)*p23(1)
      dd=dsqrt(a*a+b*b+c*c)
      a=a/dd
      b=b/dd
      c=c/dd
!
!     determining the inhomogeneous term
!
      d=-a*co(1,nodep(1))-b*co(2,nodep(1))-c*co(3,nodep(1))
!
      return
      end

