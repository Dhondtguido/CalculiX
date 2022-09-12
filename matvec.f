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
!     y=A*x for real sparse matrices (compressed row format)
!
      subroutine matvec(n,x,y,nelt,ia,ja,a,isym)
!
      implicit none
!
      integer ia(*),ja(*),i,j,n,nelt,isym,nflnei
      real*8 y(*),x(*),a(*)
!
      do i=1,n
         y(i)=a(ja(i)+1)*x(ia(ja(i)+1))
         do j=ja(i)+2,ja(i+1)
               y(i)=y(i)+a(j)*x(ia(j))
         enddo
      enddo
!
      return
      end
