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
      subroutine fricheat(et,f,vnorm,time,ciname,noel,nelems,jfaces,
     &    nelemm,jfacem,um,kstep,kinc,area,pressure,coords)
!
!     user subroutine fricheat (only for surface-to-surface contact)
!
!     INPUT:
!
!     time(1)            step time at the end of the increment
!     time(2)            total time at the end of the increment
!     ciname             surface interaction name
!     noel               element number of the contact spring element
!     nelems             slave element number
!     jfaces             local slave face number
!     nelemm             master element number
!     jfacem             local master face number
!     um                 friction coefficient
!     kstep              step number
!     kinc               increment number
!     area               slave area corresponding to the contact
!                        spring element
!     pressure           actual pressure
!     coords(1..3)       coordinates of the slave integration point
!
!     OUTPUT:
!
!     et                 portion of the work converted into heat
!                        (0 <= et <= 1)
!     f                  portion of the heat going into the slave
!                        surface; 0 <= f <= 1; the portion of the 
!                        heat going into the master surface is 1-f
!     vnorm              differential velocity between the surfaces
!                        in friction (>0)
!           
      implicit none
!
      character*80 ciname
!
      integer noel,nelems,jfaces,nelemm,jfacem,kstep,kinc
!
      real*8 vnorm,et,f,time(*),um,area,pressure,coords(3)
!
!
!
!     insert code here
!
      et=0.9
      f=0.3
      vnorm=200.
!
      return
      end

