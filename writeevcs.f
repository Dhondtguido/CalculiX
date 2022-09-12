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
      subroutine writeevcs(x,nx,nm,xmin,xmax)
!
!     writes the eigenvalues to unit 3 and replaces the 
!     eigenvalue by its square root = frequency (in rad/time)
!
!     nm is the nodal diameter
!
      implicit none
!
      integer j,nx,nm
      real*8 x(nx),pi,xmin,xmax,xnull
!
      pi=4.d0*datan(1.d0)
      xnull=0.d0
!
      write(5,*)
      write(5,*) '    E I G E N V A L U E   O U T P U T'
      write(5,*)
      write(5,*) ' NODAL   MODE NO    EIGENVALUE                      FR
     &EQUENCY'
      write(5,*) 'DIAMETER                                    REAL PART
     &           IMAGINARY PART'
      write(5,*) '                                   (RAD/TIME)      (CY
     &CLES/TIME)   (RAD/TIME)'
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
!
         if(x(j).lt.0.d0) then
            write(5,'(i5,4x,i7,4(2x,e14.7))') nm,j,x(j),xnull,
     &         xnull,dsqrt(-x(j))
         else
            write(5,'(i5,4x,i7,4(2x,e14.7))') nm,j,x(j),dsqrt(x(j)),
     &     dsqrt(x(j))/(2.d0*pi),xnull
         endif
      enddo
!
      return
      end

