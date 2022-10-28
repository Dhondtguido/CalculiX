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
      subroutine sdvini(statev,coords,nstatv,ncrds,noel,npt,
     &  layer,kspt)
!
!     user subroutine sdvini
!
!
!     INPUT:
!
!     coords(1..3)       global coordinates of the integration point
!     nstatv             number of internal variables (must be
!                        defined by the user with the *DEPVAR card)
!     ncrds              number of coordinates
!     noel               element number
!     npt                integration point number
!     layer              not used
!     kspt               not used
!
!     OUTPUT:
!
!     statev(1..nstatv)  initial value of the internal state
!                        variables
!       
      implicit none
!
      integer nstatv,ncrds,noel,npt,layer,kspt,i
!
      real*8 statev(nstatv),coords(ncrds)
!
!     code for retrieving the internal state variables
!
c      do i=1,13
      do i=1,nstatv
         statev(i)=1.d0
      enddo
      return
      end

