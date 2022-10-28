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
!     d{K(X)}/dX
!
!     author: Yannick Muller
!
      subroutine dKdX(x,u,uprime,rpar,ipar)
!
      implicit none
      integer ipar
      real*8 x,u(1),uprime(1),rpar(*),zk0,phi
!
!     defining the parameters
      phi=rpar(1)
      zk0=rpar(3)

      uprime(1)=datan(1.d0)*0.315d0/(phi)*x**1.6d0*
     &    ((zk0*u(1))**1.75d0-
     &    (dabs(1.d0-u(1)))**1.75d0*(1.d0-u(1))/dabs(1.d0-u(1)))
     &    -2.d0*u(1)/x
!
      return
!     
      end
!
