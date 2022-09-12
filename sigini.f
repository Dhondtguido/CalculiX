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
      subroutine sigini(sigma,coords,ntens,ncrds,noel,npt,layer,
     &  kspt,lrebar,rebarn)
!
!     user subroutine sigini
!
!     INPUT:
!
!     coords             coordinates of the integration point
!     ntens              number of stresses to be defined
!     ncrds              number of coordinates
!     noel               element number
!     npt                integration point number
!     layer              currently not used
!     kspt               currently not used 
!     lrebar             currently not used (value: 0)
!     rebarn             currently not used
!
!     OUTPUT:
!
!     sigma(1..ntens)    residual stress values in the integration
!                        point. If ntens=6 the order of the 
!                        components is 11,22,33,12,13,23
!           
      implicit none
!
      character*80 rebarn
      integer ntens,ncrds,noel,npt,layer,kspt,lrebar
      real*8 sigma(*),coords(*)
!
      sigma(1)=-100.d0*coords(2)
      sigma(2)=-100.d0*coords(2)
      sigma(3)=-100.d0*coords(2)
      sigma(4)=0.d0
      sigma(5)=0.d0
      sigma(6)=0.d0
!
      return
      end

