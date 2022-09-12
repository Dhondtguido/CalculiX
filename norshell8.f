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
      subroutine norshell8(xi,et,xl,xnor)
!
!     calculates the normal on a quadratic shell element in a point
!     with local coordinates xi and et. The coordinates of the nodes
!     belonging to the element are stored in xl
!
      implicit none
!
      integer i,j,k
!
      real*8 shp(4,8),xs(3,2),xl(3,8),xnor(3)
!
      real*8 xi,et
!
!     shape functions and their glocal derivatives for an element
!     described with two local parameters and three global ones.
!
!     local derivatives of the shape functions: xi-derivative
!
      shp(1,1)=(1.d0-et)*(2.d0*xi+et)/4.d0
      shp(1,2)=(1.d0-et)*(2.d0*xi-et)/4.d0
      shp(1,3)=(1.d0+et)*(2.d0*xi+et)/4.d0
      shp(1,4)=(1.d0+et)*(2.d0*xi-et)/4.d0
      shp(1,5)=-xi*(1.d0-et)
      shp(1,6)=(1.d0-et*et)/2.d0
      shp(1,7)=-xi*(1.d0+et)
      shp(1,8)=-(1.d0-et*et)/2.d0
!
!     local derivatives of the shape functions: eta-derivative
!
      shp(2,1)=(1.d0-xi)*(2.d0*et+xi)/4.d0
      shp(2,2)=(1.d0+xi)*(2.d0*et-xi)/4.d0
      shp(2,3)=(1.d0+xi)*(2.d0*et+xi)/4.d0
      shp(2,4)=(1.d0-xi)*(2.d0*et-xi)/4.d0
      shp(2,5)=-(1.d0-xi*xi)/2.d0
      shp(2,6)=-et*(1.d0+xi)
      shp(2,7)=(1.d0-xi*xi)/2.d0
      shp(2,8)=-et*(1.d0-xi)
!
!     computation of the local derivative of the global coordinates
!     (xs)
!
      do i=1,3
        do j=1,2
          xs(i,j)=0.d0
          do k=1,8
            xs(i,j)=xs(i,j)+xl(i,k)*shp(j,k)
          enddo
        enddo
      enddo
!
!     computation of the jacobian determinant
!
      xnor(1)=xs(2,1)*xs(3,2)-xs(3,1)*xs(2,2)
      xnor(2)=xs(1,2)*xs(3,1)-xs(3,2)*xs(1,1)
      xnor(3)=xs(1,1)*xs(2,2)-xs(2,1)*xs(1,2)
!
      return
      end
