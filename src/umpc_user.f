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
      subroutine umpc_user(x,u,f,a,jdof,n,force,iit,idiscon)
!
!     updates the coefficients in a user mpc
!
!     INPUT:
!
!     x(3,n)             Carthesian coordinates of the nodes in the
!                        user mpc.
!     u(3,n)             Actual displacements of the nodes in the
!                        user mpc.     
!     jdof               Actual degrees of freedom of the mpc terms
!     n                  number of terms in the user mpc
!     force              Actual value of the mpc force
!     iit                iteration number
!
!     OUTPUT:
!
!     f                  Actual value of the mpc. If the mpc is
!                        exactly satisfied, this value is zero
!     a(n)               coefficients of the linearized mpc
!     jdof               Corrected degrees of freedom of the mpc terms
!     idiscon            0: no discontinuity
!                        1: discontinuity
!                        If a discontinuity arises the previous
!                        results are not extrapolated at the start of
!                        a new increment
!
      implicit none
!
      integer jdof(*),n,iit,idiscon
!
      real*8 x(3,*),u(3,*),f,a(*),force
!
!
      return
      end













