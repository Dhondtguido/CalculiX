!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2023 Guido Dhondt
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
      subroutine shape13p(xi,et,ze,xl,xsj,shp,iflag)
!
!     shape functions and derivatives for a 13-node quadratic
!     isoparametric paramid element. 0<=xi,et<=1,-1<=ze<=1,xi+et<=1.
!
!     iflag=1: calculate only the value of the shape functions
!     iflag=2: calculate the value of the shape functions and
!              the Jacobian determinant
!     iflag=3: calculate the value of the shape functions, the
!              value of their derivatives w.r.t. the global
!              coordinates and the Jacobian determinant
!
!
!     Written September 2023 fgr
!     The shape functions is based on the 13-node Bedrosian Pyramid
!     described by 
!       Robert S. Browning
!       A SECOND-ORDER 19-NODE PYRAMID FINITE ELEMENT SUITABLE FOR
!               LUMPED MASS EXPLICIT DYNAMIC METHODS IN
!                  NONLINEAR SOLID MECHANICS
!     http://dx.doi.org/10.13140/RG.2.2.26801.20322
!
!     and G. Bedrosian
!      SHAPE FUNCTIONS AND INTEGRATION FORMULAS
!      FOR THREE-DIMENSIONAL FINITE ELEMENT ANALYSIS 
!     https://onlinelibrary.wiley.com/doi/10.1002/nme.1620350106
!
      implicit none
!
      integer i,j,k,iflag
!
      real*8 shp(4,13),xs(3,3),xsi(3,3),xl(3,13),sh(3)
!
      real*8 xi,et,ze,xsj, a, b, c
!
!     shape functions and their glocal derivatives
!
!     shape functions
!
!
      a = xi*et*ze/(1.d0-ze);
      b = 0.5d0/(1.d0-ze);
      c = ze/(1.d0-ze);
!
      shp(4,1)= 0.25d0*(-xi-et-1.d0)*
     &          ((1.d0-xi)*(1.d0-et)-ze+a)
      shp(4,2)= 0.25d0*(xi-et-1.d0)*
     &          ((1.d0+xi)*(1.d0-et)-ze-a)
      shp(4,3)= 0.25d0*(xi+et-1.d0)*
     &          ((1.d0+xi)*(1.d0+et)-ze+a)
      shp(4,4)= 0.25d0*(-xi+et-1.d0)*
     &          ((1.d0-xi)*(1.d0+et)-ze-a)
      shp(4,5)= ze*(2.d0*ze-1.d0)
      shp(4,6)= (1.d0+xi-ze)*
     &          (1.d0-xi-ze)*(1.d0-et-ze)*b
      shp(4,7)= (1.d0+et-ze)*
     &          (1.d0-et-ze)*(1.d0+xi-ze)*b
      shp(4,8)= (1.d0+xi-ze)*
     &          (1.d0-xi-ze)*(1.d0+et-ze)*b
      shp(4,9)= (1.d0+et-ze)*
     &          (1.d0-et-ze)*(1.d0-xi-ze)*b
      shp(4,10)= (1.d0-xi-ze)*(1.d0-et-ze)*c
      shp(4,11)= (1.d0+xi-ze)*(1.d0-et-ze)*c
      shp(4,12)= (1.d0+xi-ze)*(1.d0+et-ze)*c
      shp(4,13)= (1.d0-xi-ze)*(1.d0+et-ze)*c
!
      if(iflag.eq.1) return
!
!     local derivatives of the shape functions: xi-derivative
!
      shp(1,1)= (ze*ze+(2.d0*xi+2.d0*et-1.d0)*ze
     &         +(2.d0*et-2.d0)*xi+et*et-et)/(4.d0*ze-4.d0)
      shp(1,2)=-(ze*ze+(-2.d0*xi+2.d0*et-1.d0)*ze
     &         +(2.d0-2.d0*et)*xi+et*et-et)/(4.d0*ze-4.d0)
      shp(1,3)=-(ze*ze+(-2.d0*xi-2.d0*et-1.d0)*ze
     &         +(2.d0*et+2.d0)*xi+et*et+et)/(4.d0*ze-4.d0)
      shp(1,4)= (ze*ze+(2.d0*xi-2.d0*et-1.d0)*ze
     &         +(-2.d0*et-2.d0)*xi+et*et+et)/(4.d0*ze-4.d0)
!
      shp(1,5)=0.d0
!	  
      shp(1,6)=-(xi*ze+(et-1.d0)*xi)/(ze-1.d0)
      shp(1,7)=-(ze*ze-2.d0*ze-et*et+1.d0)/(2.d0*ze-2.d0)
      shp(1,8)=-(xi*ze+(-et-1.d0)*xi)/(ze-1.d0)
      shp(1,9)= (ze*ze-2.d0*ze-et*et+1.d0)/(2.d0*ze-2.d0)
!
      shp(1,10)=-(ze*ze+(et-1.d0)*ze)/(ze-1.d0)
      shp(1,11)= (ze*ze+(et-1.d0)*ze)/(ze-1.d0)
      shp(1,12)= (ze*ze+(-et-1.d0)*ze)/(ze-1.d0)
      shp(1,13)=-(ze*ze+(-et-1.d0)*ze)/(ze-1.d0)
!
!     local derivatives of the shape functions: eta-derivative
!
      shp(2,1)= (ze*ze+(2.d0*xi+2.d0*et-1.d0)*ze+xi*xi
     &         +(2.d0*et-1.d0)*xi-2.d0*et)/(4.d0*ze-4.d0)
      shp(2,2)= (ze*ze+(-2.d0*xi+2.d0*et-1.d0)*ze+xi*xi
     &         +(1.d0-2.d0*et)*xi-2.d0*et)/(4.d0*ze-4.d0)
      shp(2,3)=-(ze*ze+(-2.d0*xi-2.d0*et-1.d0)*ze+xi*xi
     &         +(2.d0*et+1.d0)*xi+2.d0*et)/(4.d0*ze-4.d0)
      shp(2,4)=-(ze*ze+(2.d0*xi-2.d0*et-1.d0)*ze+xi*xi
     &  +(-2.d0*et-1.d0)*xi+2.d0*et)/(4.d0*ze-4.d0)
!
      shp(2,5)=0.d0
!
      shp(2,6)= (ze*ze-2.d0*ze-xi*xi+1.d0)/(2.d0*ze-2.d0)
      shp(2,7)=-(et*ze-et*xi-et)/(ze-1.d0)
      shp(2,8)=-(ze*ze-2.d0*ze-xi*xi+1.d0)/(2.d0*ze-2.d0)
      shp(2,9)=-(et*ze+et*xi-et)/(ze-1.d0)
!
      shp(2,10)=-(ze*ze+(xi-1.d0)*ze)/(ze-1.d0)
      shp(2,11)=-(ze*ze+(-xi-1.d0)*ze)/(ze-1.d0)
      shp(2,12)= (ze*ze+(-xi-1.d0)*ze)/(ze-1.d0)
      shp(2,13)= (ze*ze+(xi-1.d0)*ze)/(ze-1.d0)
!
!     local derivatives of the shape functions: mu/zeta-derivative
!
      shp(3,1)= ((xi+et+1.d0)*ze*ze+(-2.d0*xi-2.d0*et-2.d0)*ze
     &          -et*xi*xi+(-et*et-et+1.d0)*xi+et+1.d0)
     &         /(4.d0*ze*ze-8.d0*ze+4.d0)
      shp(3,2)=-((xi-et-1.d0)*ze*ze+(-2.d0*xi+2.d0*et+2.d0)*ze
     &         +et*xi*xi+(-et*et-et+1.d0)*xi-et-1.d0)
     &         /(4.d0*ze*ze-8.d0*ze+4.d0)
      shp(3,3)=-((xi+et-1.d0)*ze*ze+(-2.d0*xi-2.d0*et+2.d0)*ze
     &         -et*xi*xi+(-et*et+et+1.d0)*xi+et-1.d0)
     &         /(4.d0*ze*ze-8.d0*ze+4.d0)
      shp(3,4)= ((xi-et+1.d0)*ze*ze+(-2.d0*xi+2.d0*et-2.d0)*ze
     &         +et*xi*xi+(-et*et+et+1.d0)*xi-et+1.d0)
     &         /(4.d0*ze*ze-8.d0*ze+4.d0)
!
      shp(3,5)=4.d0*ze-1.d0
!
      shp(3,6)= (2.d0*ze*ze*ze+(et-6.d0)*ze*ze+(6.d0-2.d0*et)*ze
     &         +et*xi*xi+et-2.d0)/(2.d0*ze*ze-4.d0*ze+2.d0)
      shp(3,7)= (2.d0*ze*ze*ze+(-xi-6.d0)*ze*ze+(2.d0*xi+6.d0)*ze
     &         +(-et*et-1.d0)*xi-2.d0)/(2.d0*ze*ze-4.d0*ze+2.d0)
      shp(3,8)= (2.d0*ze*ze*ze+(-et-6.d0)*ze*ze+(2.d0*et+6.d0)*ze
     &         -et*xi*xi-et-2.d0)/(2.d0*ze*ze-4.d0*ze+2.d0)
      shp(3,9)= (2.d0*ze*ze*ze+(xi-6.d0)*ze*ze+(6.d0-2.d0*xi)*ze
     &         +(et*et+1.d0)*xi-2.d0)/(2.d0*ze*ze-4.d0*ze+2.d0)
!
      shp(3,10)=-(2.d0*ze*ze*ze+(xi+et-5.d0)*ze*ze
     &          +(-2.d0*xi-2.d0*et+4.d0)*ze+(1.d0-et)*xi
     &          +et-1.d0)/(ze*ze-2.d0*ze+1.d0)
      shp(3,11)=-(2.d0*ze*ze*ze+(-xi+et-5.d0)*ze*ze
     &          +(2.d0*xi-2.d0*et+4.d0)*ze+(et-1.d0)*xi
     &          +et-1.d0)/(ze*ze-2.d0*ze+1.d0)
      shp(3,12)=-(2.d0*ze*ze*ze+(-xi-et-5.d0)*ze*ze
     &          +(2.d0*xi+2.d0*et+4.d0)*ze+(-et-1.d0)*xi
     &          -et-1.d0)/(ze*ze-2.d0*ze+1.d0)
      shp(3,13)=-(2.d0*ze*ze*ze+(xi-et-5.d0)*ze*ze
     &          +(-2.d0*xi+2.d0*et+4.d0)*ze+(et+1.d0)*xi
     &          -et-1.d0)/(ze*ze-2.d0*ze+1.d0)
!
!     computation of the local derivative of the global coordinates
!     (xs)
!
      do i=1,3
        do j=1,3
          xs(i,j)=0.d0
          do k=1,13
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
      if(iflag.eq.2) return
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
      do k=1,13
        do j=1,3
          sh(j)=shp(1,k)*xsi(1,j)+shp(2,k)*xsi(2,j)+shp(3,k)*xsi(3,j)
        enddo
        do j=1,3
          shp(j,k)=sh(j)
        enddo
      enddo
!
      return
      end
