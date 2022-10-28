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
      subroutine uamplitude(time,name,amplitude)
!
!     user subroutine uamplitude: user defined amplitude definition
!
!     INPUT:
!
!     name               amplitude name
!     time               time at which the amplitude is to be
!                        evaluated
!
!     OUTPUT:
!
!     amplitude          value of the amplitude at time
!           
      implicit none
!
      character*80 name
!
      real*8 time,amplitude
!
      if(name(1:9).eq.'QUADRATIC') then
         amplitude=time**2
      else
         write(*,*) '*ERROR in uamplitude: unknown amplitude'
         call exit(201)
      endif
!
      return
      end

