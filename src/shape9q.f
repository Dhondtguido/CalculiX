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
      subroutine shape9q(xi,et,xl,xsj,xs,shp,iflag)
!
!     shape functions and derivatives for a 9-node quadratic
!     isoparametric quadrilateral element. -1<=xi,et<=1 
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
!
      implicit none
!
      integer i,j,k,iflag
!
      real*8 shp(7,9),xs(3,7),xsi(2,3),xl(3,9),sh(3),xsj(3),xi,et,
     &  fxi1,fxi2,fxi3,fet1,fet2,fet3,dfxi1,dfxi2,dfxi3,dfet1,dfet2,
     &  dfet3,ddfxi1,ddfxi2,ddfxi3,ddfet1,ddfet2,ddfet3
!
!
!
!     shape functions and their glocal derivatives for an element
!     described with two local parameters and three global ones.
!
!     shape functions in one dimension
!
      fxi1=xi*(xi-1.d0)/2.d0
      fxi2=(1.d0-xi)*(1.d0+xi)
      fxi3=xi*(xi+1.d0)/2.d0
!
      fet1=et*(et-1.d0)/2.d0
      fet2=(1.d0-et)*(1.d0+et)
      fet3=et*(et+1.d0)/2.d0
!
!     shape functions
!
      shp(4,1)=fxi1*fet1
      shp(4,2)=fxi3*fet1
      shp(4,3)=fxi3*fet3
      shp(4,4)=fxi1*fet3
      shp(4,5)=fxi2*fet1
      shp(4,6)=fxi3*fet2
      shp(4,7)=fxi2*fet3
      shp(4,8)=fxi1*fet2
      shp(4,9)=fxi2*fet2
!
      if(iflag.eq.1) return
!
!     derivative of the shape functions in one dimension
!
      dfxi1=(2.d0*xi-1.d0)/2.d0
      dfxi2=-2.d0*xi
      dfxi3=(2.d0*xi+1.d0)/2.d0
!
      dfet1=(2.d0*et-1.d0)/2.d0
      dfet2=-2.d0*et
      dfet3=(2.d0*et+1.d0)/2.d0
!
!     local derivatives of the shape functions: xi-derivative
!
      shp(1,1)=dfxi1*fet1
      shp(1,2)=dfxi3*fet1
      shp(1,3)=dfxi3*fet3
      shp(1,4)=dfxi1*fet3
      shp(1,5)=dfxi2*fet1
      shp(1,6)=dfxi3*fet2
      shp(1,7)=dfxi2*fet3
      shp(1,8)=dfxi1*fet2
      shp(1,9)=dfxi2*fet2
!
!     local derivatives of the shape functions: eta-derivative
!
      shp(2,1)=fxi1*dfet1
      shp(2,2)=fxi3*dfet1
      shp(2,3)=fxi3*dfet3
      shp(2,4)=fxi1*dfet3
      shp(2,5)=fxi2*dfet1
      shp(2,6)=fxi3*dfet2
      shp(2,7)=fxi2*dfet3
      shp(2,8)=fxi1*dfet2
      shp(2,9)=fxi2*dfet2
!
!     computation of the local derivative of the global coordinates
!     (xs)
!
      do i=1,3
        do j=1,2
          xs(i,j)=0.d0
          do k=1,9
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
         do k=1,9
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
!     second derivative of the shape functions in one dimension
!
         ddfxi1=1.d0
         ddfxi2=-2.d0
         ddfxi3=1.d0
!     
         ddfet1=1.d0
         ddfet2=-2.d0
         ddfet3=1.d0
!
!     local 2nd order derivatives of the shape functions: xi,xi-derivative
!     
         shp(5,1)=ddfxi1*fet1
         shp(5,2)=ddfxi3*fet1
         shp(5,3)=ddfxi3*fet3
         shp(5,4)=ddfxi1*fet3
         shp(5,5)=ddfxi2*fet1
         shp(5,6)=ddfxi3*fet2
         shp(5,7)=ddfxi2*fet3
         shp(5,8)=ddfxi1*fet2
         shp(5,9)=ddfxi2*fet2
!     
!     local 2nd order derivatives of the shape functions: xi,eta-derivative
!     
         shp(6,1)=dfxi1*dfet1
         shp(6,2)=dfxi3*dfet1
         shp(6,3)=dfxi3*dfet3
         shp(6,4)=dfxi1*dfet3
         shp(6,5)=dfxi2*dfet1
         shp(6,6)=dfxi3*dfet2
         shp(6,7)=dfxi2*dfet3
         shp(6,8)=dfxi1*dfet2
         shp(6,9)=dfxi2*dfet2
!     
!     local 2nd order derivatives of the shape functions: eta,eta-derivative
!     
         shp(7,1)=fxi1*ddfet1
         shp(7,2)=fxi3*ddfet1
         shp(7,3)=fxi3*ddfet3
         shp(7,4)=fxi1*ddfet3
         shp(7,5)=fxi2*ddfet1
         shp(7,6)=fxi3*ddfet2
         shp(7,7)=fxi2*ddfet3
         shp(7,8)=fxi1*ddfet2
         shp(7,9)=fxi2*ddfet2
!     
!     computation of the local 2nd derivatives of the global coordinates
!     (xs)
!
         do i=1,3
            do j=5,7
               xs(i,j)=0.d0
               do k=1,9
                  xs(i,j)=xs(i,j)+xl(i,k)*shp(j,k)
               enddo
            enddo
         enddo
      endif
!     
      return
      end
