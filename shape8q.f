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
      subroutine shape8q(xi,et,xl,xsj,xs,shp,iflag)
!
!     shape functions and derivatives for a 8-node quadratic
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
!     iflag=5: calculate the value of the shape functions and
!              their derivatives w.r.t. the local coordinates
!
      implicit none
!
      integer i,j,k,iflag
!
      real*8 shp(7,8),xs(3,7),xsi(2,3),xl(3,8),sh(3),xsj(3),xi,et,
     &  xip,xim,xim2,etp,etm,etm2,xipet,ximet
!
!     shape functions and their glocal derivatives for an element
!     described with two local parameters and three global ones.
!
      xip=1.d0+xi
      xim=1.d0-xi
      xim2=xip*xim
!
      etp=1.d0+et
      etm=1.d0-et
      etm2=etp*etm
!
      xipet=xi+et
      ximet=xi-et
!
!     shape functions
!
      shp(4,1)=xim*etm*(-xipet-1.d0)/4.d0
      shp(4,2)=xip*etm*(ximet-1.d0)/4.d0
      shp(4,3)=xip*etp*(xipet-1.d0)/4.d0
      shp(4,4)=xim*etp*(-ximet-1.d0)/4.d0
      shp(4,5)=xim2*etm/2.d0
      shp(4,6)=xip*etm2/2.d0
      shp(4,7)=xim2*etp/2.d0
      shp(4,8)=xim*etm2/2.d0
!
      if(iflag.eq.1) return
!
!     local derivatives of the shape functions: xi-derivative
!
      shp(1,1)=etm*(xi+xipet)/4.d0
      shp(1,2)=etm*(xi+ximet)/4.d0
      shp(1,3)=etp*(xi+xipet)/4.d0
      shp(1,4)=etp*(xi+ximet)/4.d0
      shp(1,5)=-xi*etm
      shp(1,6)=etm2/2.d0
      shp(1,7)=-xi*etp
      shp(1,8)=-etm2/2.d0
!
!     local derivatives of the shape functions: eta-derivative
!
      shp(2,1)=xim*(et+xipet)/4.d0
      shp(2,2)=xip*(et-ximet)/4.d0
      shp(2,3)=xip*(et+xipet)/4.d0
      shp(2,4)=xim*(et-ximet)/4.d0
      shp(2,5)=-xim2/2.d0
      shp(2,6)=-et*xip
      shp(2,7)=xim2/2.d0
      shp(2,8)=-et*xim
!
      if(iflag.eq.5) return
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
         do k=1,8
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
         shp(5,1)=etm/2.d0
         shp(5,2)=etm/2.d0
         shp(5,3)=etp/2.d0
         shp(5,4)=etp/2.d0
         shp(5,5)=-etm
         shp(5,6)=0.d0
         shp(5,7)=-etp
         shp(5,8)=0.d0
!
!     local 2nd order derivatives of the shape functions: xi,eta-derivative
!     
         shp(6,1)=(1.d0-2.d0*xipet)/4.d0
         shp(6,2)=(-1.d0-2.d0*ximet)/4.d0
         shp(6,3)=(1.d0+2.d0*xipet)/4.d0
         shp(6,4)=(-1.d0+2.d0*ximet)/4.d0
         shp(6,5)=xi
         shp(6,6)=-et
         shp(6,7)=-xi
         shp(6,8)=et
!     
!     local 2nd order derivatives of the shape functions: eta,eta-derivative
!     
         shp(7,1)=xim/2.d0
         shp(7,2)=xip/2.d0
         shp(7,3)=xip/2.d0
         shp(7,4)=xim/2.d0
         shp(7,5)=0.d0
         shp(7,6)=-xip
         shp(7,7)=0.d0
         shp(7,8)=-xim
!
!     computation of the local 2nd derivatives of the global coordinates
!     (xs)
!
         do i=1,3
            do j=5,7
               xs(i,j)=0.d0
               do k=1,8
                  xs(i,j)=xs(i,j)+xl(i,k)*shp(j,k)
               enddo
            enddo
         enddo
      endif
!     
      return
      end
