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
      subroutine shape2l(xi,xl,xsj,xs,shp,iflag)
!
!     shape functions and derivatives for a 2-node linear
!     isoparametric 1-D element. -1<=xi<=1
!
!     iflag=2: calculate the value of the shape functions,
!              their derivatives w.r.t. the local coordinates
!              and the Jacobian (size of tangent vector to the
!              curved line)
!
      implicit none
!
      integer i,k,iflag
!
      real*8 shp(7,3),xs(3,7),xsi(2,3),xl(3,3),sh(3),xsj(3),xi
!
!
!
!     shape functions and their glocal derivatives for an element
!     described with one local parameter and three global ones.
!
!     local derivatives of the shape functions: xi-derivative
!
      shp(1,1)=-0.5d0
      shp(1,2)=0.5d0
!
!     shape functions
!
      shp(4,1)=(1.d0-xi)/2.d0
      shp(4,2)=(1.d0+xi)/2.d0
!
!     computation of the local derivative of the global coordinates
!     (xs)
!
      do i=1,3
         xs(i,1)=0.d0
         do k=1,2
            xs(i,1)=xs(i,1)+xl(i,k)*shp(1,k)
         enddo
      enddo
!
!     computation of the jacobian vector
!
c      xsj(1)=dsqrt(xs(1,1)**2+xs(2,1)**2+xs(3,1)**2)
!
      xsj(1)=xs(1,1)
      xsj(2)=xs(2,1)
      xsj(3)=xs(3,1)
!     
      return
      end



