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
      subroutine writeevcomplex(x,nx,fmin,fmax)
!
!     writes the eigenvalues in the .dat file
!
      implicit none
!
      integer j,nx
      real*8 pi,fmin,fmax
      complex*16 x(nx)
!     
      pi=4.d0*datan(1.d0)
!
      write(5,*)
      write(5,*) '    E I G E N V A L U E   O U T P U T'
      write(5,*)
      write(5,*) 'MODE NO                     FREQUENCY          
     &  '
      write(5,*) '                     REAL PART         IMAGINARY PART'
      write(5,*) '            (RAD/TIME)   (CYCLES/TIME)   (RAD/TIME)'
      write(5,*)
      do j=1,nx
         write(5,'(i7,4(2x,e14.7))') j,dreal(x(j)),
     &        dreal(x(j))/(2.d0*pi),dimag(x(j))
      enddo
!
      return
      end

