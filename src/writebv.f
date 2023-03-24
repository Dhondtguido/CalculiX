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
      subroutine writebv(x,nx)
!
!     writes the buckling force factor to unit 5
!
      implicit none
!
      integer j,nx
      real*8 x(nx),pi
!
      pi=4.d0*datan(1.d0)
!
      write(5,*)
      write(5,*) '    B U C K L I N G   F A C T O R   O U T P U T'
      write(5,*)
      write(5,*) 'MODE NO       BUCKLING'
      write(5,*) '               FACTOR'
      write(5,*)
      do j=1,nx
         write(5,'(i7,2x,e14.7)') j,x(j)
      enddo
!
      return
      end

