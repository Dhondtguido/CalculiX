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
      subroutine shape8h(xi,et,ze,xl,xsj,shp,iflag)
!
!     shape functions and derivatives for a 8-node linear isoparametric
!     solid element
!
!     iflag=1: calculate only the value of the shape functions
!     iflag=2: calculate the value of the shape functions and
!              the Jacobian determinant
!     iflag=3: calculate the value of the shape functions, the
!              value of their derivatives w.r.t. the global
!              coordinates and the Jacobian determinant
!
      implicit none
!
      integer i,j,k,iflag
!
      real*8 shp(4,20),xs(3,3),xsi(3,3),xl(3,20),sh(3)
!
      real*8 xi,et,ze,xsj,omg,omh,omr,opg,oph,opr
!
!
!
!     shape functions and their glocal derivatives
!
      omg=1.d0-xi
      omh=1.d0-et
      omr=1.d0-ze
      opg=1.d0+xi
      oph=1.d0+et
      opr=1.d0+ze
!
!     shape functions
!
      shp(4, 1)=omg*omh*omr/8.d0
      shp(4, 2)=opg*omh*omr/8.d0
      shp(4, 3)=opg*oph*omr/8.d0
      shp(4, 4)=omg*oph*omr/8.d0
      shp(4, 5)=omg*omh*opr/8.d0
      shp(4, 6)=opg*omh*opr/8.d0
      shp(4, 7)=opg*oph*opr/8.d0
      shp(4, 8)=omg*oph*opr/8.d0
!
      if(iflag.eq.1) return
!
!     local derivatives of the shape functions: xi-derivative
!
      shp(1, 1)=-omh*omr
      shp(1, 2)=omh*omr
      shp(1, 3)=oph*omr
      shp(1, 4)=-oph*omr
      shp(1, 5)=-omh*opr
      shp(1, 6)=omh*opr
      shp(1, 7)=oph*opr
      shp(1, 8)=-oph*opr
!
!     local derivatives of the shape functions: eta-derivative
!
      shp(2, 1)=-omg*omr
      shp(2, 2)=-opg*omr
      shp(2, 3)=opg*omr
      shp(2, 4)=omg*omr
      shp(2, 5)=-omg*opr
      shp(2, 6)=-opg*opr
      shp(2, 7)=opg*opr
      shp(2, 8)=omg*opr
!
!     local derivatives of the shape functions: zeta-derivative
!
      shp(3, 1)=-omg*omh
      shp(3, 2)=-opg*omh
      shp(3, 3)=-opg*oph
      shp(3, 4)=-omg*oph
      shp(3, 5)=omg*omh
      shp(3, 6)=opg*omh
      shp(3, 7)=opg*oph
      shp(3, 8)=omg*oph
!
!     computation of the local derivative of the global coordinates
!     (xs)
!
      do i=1,3
        do j=1,3
          xs(i,j)=0.d0
          do k=1,8
            xs(i,j)=xs(i,j)+xl(i,k)*shp(j,k)
          enddo
        enddo
      enddo
!
!     computation of the jacobian determinant
!
      xsj=xs(1,1)*(xs(2,2)*xs(3,3)-xs(2,3)*xs(3,2))
     &   -xs(1,2)*(xs(2,1)*xs(3,3)-xs(2,3)*xs(3,1))
     &   +xs(1,3)*(xs(2,1)*xs(3,2)-xs(2,2)*xs(3,1))
!
      if(iflag.eq.2) then
         xsj=xsj/512.d0
         return
      endif
!
!     computation of the global derivative of the local coordinates
!     (xsi) (inversion of xs)
!
      xsi(1,1)=(xs(2,2)*xs(3,3)-xs(3,2)*xs(2,3))/xsj
      xsi(1,2)=(xs(1,3)*xs(3,2)-xs(1,2)*xs(3,3))/xsj
      xsi(1,3)=(xs(1,2)*xs(2,3)-xs(2,2)*xs(1,3))/xsj
      xsi(2,1)=(xs(2,3)*xs(3,1)-xs(2,1)*xs(3,3))/xsj
      xsi(2,2)=(xs(1,1)*xs(3,3)-xs(3,1)*xs(1,3))/xsj
      xsi(2,3)=(xs(1,3)*xs(2,1)-xs(1,1)*xs(2,3))/xsj
      xsi(3,1)=(xs(2,1)*xs(3,2)-xs(3,1)*xs(2,2))/xsj
      xsi(3,2)=(xs(1,2)*xs(3,1)-xs(1,1)*xs(3,2))/xsj
      xsi(3,3)=(xs(1,1)*xs(2,2)-xs(2,1)*xs(1,2))/xsj
!
!     computation of the global derivatives of the shape functions
!
      do k=1,8
        do j=1,3
          sh(j)=shp(1,k)*xsi(1,j)+shp(2,k)*xsi(2,j)+shp(3,k)*xsi(3,j)
        enddo
        do j=1,3
          shp(j,k)=sh(j)
        enddo
      enddo
!
      xsj=xsj/512.d0
!
      return
      end
