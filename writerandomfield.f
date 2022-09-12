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
      subroutine writerandomfield(d,relerr,imodes)
!
!     writes the error measures of the randomfield in the .dat file
!
      implicit none
!
      integer imodes
      real*8 d,relerr
!
      if(imodes.eq.1) then
         write(5,*)
         write(5,*) 'SPECTRAL DECOMPOSITION OF RANDOMFIELD'
         write(5,*) 'MODESHAPE   EIGENVALUE   GLOBAL RELIABILITY'  
         write(5,*)
         write(5,'(1x,i3.3,6x,e13.4,e13.4)') imodes,d,relerr
      else
         write(5,'(1x,i3.3,6x,e13.4,e13.4)') imodes,d,relerr
      endif
!
      return
      end

