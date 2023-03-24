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
      subroutine shape4q(xi,et,xl,xsj,xs,shp,iflag)
!
!     iflag=1: calculate only the value of the shape functions
!     iflag=2: calculate the value of the shape functions,
!              their derivatives w.r.t. the local coordinates
!              and the Jacobian vector (local normal to the
!              surface)
!     iflag=3: calculate the value of the shape functions, the
!              value of their derivatives w.r.t. the global
!              coordinates and the Jacobian vector (local normal
!              to the surface)
!     iflag=4: calculate the value of the shape functions, the
!              value of their 1st and 2nd order derivatives 
!              w.r.t. the local coordinates, the Jacobian vector 
!              (local normal to the surface)
!     iflag=5: calculate the value of the shape functions and
!              their derivatives w.r.t. the local coordinates
!
      implicit none
!
      integer i,j,k,iflag
!
      real*8 shp(7,4),xs(3,7),xsi(2,3),xl(3,8),sh(3),xsj(3),xi,et,
     &  xip,xim,etp,etm
!
!
!
!     shape functions and their glocal derivatives for an element
!     described with two local parameters and three global ones.
!
      xip=1.d0+xi
      xim=1.d0-xi
      etp=1.d0+et
      etm=1.d0-et
!
!     shape functions
!
      shp(4,1)=xim*etm/4.d0
      shp(4,2)=xip*etm/4.d0
      shp(4,3)=xip*(etp)/4.d0
      shp(4,4)=xim*(etp)/4.d0
!
      if(iflag.eq.1) return
!
!     local derivatives of the shape functions: xi-derivative
!      
      shp(1,1)=-etm/4.d0
      shp(1,2)=etm/4.d0
      shp(1,3)=(etp)/4.d0
      shp(1,4)=-(etp)/4.d0
!
!     local derivatives of the shape functions: eta-derivative
!
      shp(2,1)=-xim/4.d0
      shp(2,2)=-xip/4.d0
      shp(2,3)=xip/4.d0
      shp(2,4)=xim/4.d0
!
      if(iflag.eq.5) return
!
!     computation of the local derivative of the global coordinates
!     (xs)
!
      do i=1,3
        do j=1,2
          xs(i,j)=0.d0
          do k=1,4
            xs(i,j)=xs(i,j)+xl(i,k)*shp(j,k)
          enddo
        enddo
      enddo
!
!     computation of the jacobian vector
!
      xsj(1)=xs(2,1)*xs(3,2)-xs(3,1)*xs(2,2)
      xsj(2)=xs(1,2)*xs(3,1)-xs(3,2)*xs(1,1)
      xsj(3)=xs(1,1)*xs(2,2)-xs(2,1)*xs(1,2)
!
      if(iflag.eq.3) then
!     
!     computation of the global derivative of the local coordinates
!     (xsi) (inversion of xs)
!     
c         xsi(1,1)=xs(2,2)/xsj(3)
c         xsi(2,1)=-xs(2,1)/xsj(3)
c         xsi(1,2)=-xs(1,2)/xsj(3)
c         xsi(2,2)=xs(1,1)/xsj(3)
c         xsi(1,3)=-xs(2,2)/xsj(1)
c         xsi(2,3)=xs(2,1)/xsj(1)
         if(dabs(xsj(3)).gt.1.d-10) then
            xsi(1,1)=xs(2,2)/xsj(3)
            xsi(2,2)=xs(1,1)/xsj(3)
            xsi(1,2)=-xs(1,2)/xsj(3)
            xsi(2,1)=-xs(2,1)/xsj(3)
            if(dabs(xsj(2)).gt.1.d-10) then
               xsi(2,3)=xs(1,1)/(-xsj(2))
               xsi(1,3)=-xs(1,2)/(-xsj(2))
            elseif(dabs(xsj(1)).gt.1.d-10) then
               xsi(2,3)=xs(2,1)/xsj(1)
               xsi(1,3)=-xs(2,2)/xsj(1)
            else
               xsi(2,3)=0.d0
               xsi(1,3)=0.d0
            endif
         elseif(dabs(xsj(2)).gt.1.d-10) then
            xsi(1,1)=xs(3,2)/(-xsj(2))
            xsi(2,3)=xs(1,1)/(-xsj(2))
            xsi(1,3)=-xs(1,2)/(-xsj(2))
            xsi(2,1)=-xs(3,1)/(-xsj(2))
            if(dabs(xsj(1)).gt.1.d-10) then
               xsi(1,2)=xs(3,2)/xsj(1)
               xsi(2,2)=-xs(3,1)/xsj(1)
            else
               xsi(1,2)=0.d0
               xsi(2,2)=0.d0
            endif
         else
            xsi(1,2)=xs(3,2)/xsj(1)
            xsi(2,3)=xs(2,1)/xsj(1)
            xsi(1,3)=-xs(2,2)/xsj(1)
            xsi(2,2)=-xs(3,1)/xsj(1)
            xsi(1,1)=0.d0
            xsi(2,1)=0.d0
         endif
!     
!     computation of the global derivatives of the shape functions
!     
         do k=1,4
            do j=1,3
               sh(j)=shp(1,k)*xsi(1,j)+shp(2,k)*xsi(2,j)
            enddo
            do j=1,3
               shp(j,k)=sh(j)
            enddo
         enddo
!
      elseif(iflag.eq.4) then
!
!     local 2nd order derivatives of the shape functions: xi,xi-derivative
!     
         shp(5,1)=0.d0
         shp(5,2)=0.d0
         shp(5,3)=0.d0
         shp(5,4)=0.d0
!
!     local 2nd order derivatives of the shape functions: xi,eta-derivative
!     
         shp(6,1)=0.25d0
         shp(6,2)=-0.25d0
         shp(6,3)=0.25d0
         shp(6,4)=-0.25d0
!     
!     local 2nd order derivatives of the shape functions: eta,eta-derivative
!     
         shp(7,1)=0.d0
         shp(7,2)=0.d0
         shp(7,3)=0.d0
         shp(7,4)=0.d0
!
!     computation of the local 2nd derivatives of the global coordinates
!     (xs)
!
         do i=1,3
            xs(i,5)=0.d0
            xs(i,7)=0.d0
         enddo
         do i=1,3
            xs(i,6)=0.d0
            do k=1,4
               xs(i,6)=xs(i,6)+xl(i,k)*shp(6,k)
            enddo
         enddo
      endif
!     
      return
      end
