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
      subroutine writestadiv(istep,j,icutb,l,ttime,time,dtime)
!
      implicit none
!
!     writes increment statistics in the .sta file
!     the close and open guarantees that the computer buffer is
!     emptied each time a new line is written. That way the file
!     is always up to data (also during the calculation)
!
!     this version of writesummary is meant for increments which did
!     not converge
!
      integer istep,j,icutb,l
      real*8 ttime,time,dtime
!
      write(8,100) istep,j,icutb+1,l,ttime+time-dtime,time-dtime,dtime
      flush(8)
!
 100  format(1x,i5,1x,i10,1x,i5,'U',1x,i4,3(1x,e13.6))
!
      return
      end
