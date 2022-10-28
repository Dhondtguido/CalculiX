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
      subroutine writeev(x,nx,xmin,xmax)
!
!     writes the eigenvalues in the .dat file and replaces the 
!     eigenvalue by its square root = frequency (in rad/time)
!
      implicit none
!
      integer j,nx
      real*8 x(nx),pi,xmin,xmax,xnull
!     
!     for real symmetric matrices the eigenvalue is real;
!     the frequency, which is the square root of the eigenvalue
!     can be real or complex (in the latter case buckling occurs)
!
      pi=4.d0*datan(1.d0)
      xnull=0.d0
!
      write(5,*)
      write(5,*) '    E I G E N V A L U E   O U T P U T'
      write(5,*)
      write(5,*) 'MODE NO    EIGENVALUE                       FREQUENCY          
     &  '
      write(5,*) '                                    REAL PART          
     &   IMAGINARY PART'
      write(5,*) '                          (RAD/TIME)      (CYCLES/TIME
     &     (RAD/TIME)'
      write(5,*)
!
      do j=1,nx
!
!        user-defined minimum frequency
!
         if(xmin.gt.-0.5d0) then
            if(xmin*xmin.gt.x(j)) cycle
         endif
!
!        user-defined maximum frequency
!
         if(xmax.gt.-0.5d0) then
            if(xmax*xmax.lt.x(j)) exit
         endif
!
         if(x(j).lt.0.d0) then
            write(5,'(i7,4(2x,e14.7))') j,x(j),xnull,
     &         xnull,dsqrt(-x(j))
         else
            write(5,'(i7,4(2x,e14.7))') j,x(j),dsqrt(x(j)),
     &     dsqrt(x(j))/(2.d0*pi),xnull
         endif
      enddo
!
      return
      end

