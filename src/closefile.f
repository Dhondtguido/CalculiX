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
      subroutine closefile()
      implicit none
!
!     closes files at the end of the calculation
!
      logical rout
!
!     closing the .inp file
!
      close(1)
!
!     closing the .dat file
!
c      flush(5)
      close(5)
!
!     closing the .sta file
!
      close(8)
!
!     closing the .cvg file
!
      close(11)
!
!     closing the .rout file
!
      inquire(15,opened=rout)
      if(rout) close(15)
!     
      return
      end
