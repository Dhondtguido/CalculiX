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
      subroutine insertsorti(ix,n)
!
!     simple insertion sort routine for very small n
!
!     https://en.wikipedia.org/wiki/Insertion_sort
!
!     Author: Lukas Mayrhofer
!      
      implicit none
!
      integer n,ix(*)
!
      integer i,j,itmp
!
      do i=2,n
         itmp=ix(i)
         do j=i-1,1,-1
            if(itmp.lt.ix(j)) then
               ix(j+1)=ix(j)
            else
               exit
            endif
         enddo
         ix(j+1)=itmp
      enddo
!
      return
      end
