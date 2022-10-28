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
      subroutine getversion(version)
!
!     this routine is called if an inconsistency is noticed between
!     the element count and the number of elements stored in the frd-
!     file. Such an inconsistency will lead to a crash while reading
!     a binary frd-file
!
      implicit none
!
      character*80 version
!
      integer i
!
!
      do i=1,80
         version(i:i)=' '
      enddo
!
      version='Version DEVELOPMENT'
!
      return
      end
