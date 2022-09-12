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
      subroutine mult(matrix,trans,n)
!
      implicit none
!
      integer i,j,k,n
      real*8 matrix(3,3),trans(3,3),a(3,3)
!
!     3x3 matrix multiplication. If n=1 then
!        matrix=trans^T*matrix,
!     if n=2 then
!        matrix=matrix*trans.
!
      if(n.eq.1) then
         do i=1,3
            do j=1,3
               a(i,j)=0.d0
               do k=1,3
                  a(i,j)=a(i,j)+trans(k,i)*matrix(k,j)
               enddo
            enddo
         enddo
      elseif(n.eq.2) then
         do i=1,3
            do j=1,3
               a(i,j)=0.d0
               do k=1,3
                  a(i,j)=a(i,j)+matrix(i,k)*trans(k,j)
               enddo
            enddo
         enddo
      endif
!
      do i=1,3
         do j=1,3
            matrix(i,j)=a(i,j)
         enddo
      enddo
!
      return
      end
         
