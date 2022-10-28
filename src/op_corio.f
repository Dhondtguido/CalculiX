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
!     y=A*x for real sparse antisymmetric matrices,
!     i.e. the transpose is the negative matrix: A^T=-A
!
!     storage of the matrix:
!        au: first lower triangle
!        ad: diagonal terms
!
      subroutine op_corio(n,x,y,ad,au,jq,irow)
!
      implicit none
!
      integer irow(*),n,j,l,i,jq(*)
      real*8 y(*),x(*),au(*),ad(*)
!
!     diagonal terms
!
      do i=1,n
         y(i)=ad(i)*x(i)
      enddo
!
!     off-diagonal terms
!
      do j=1,n
         do l=jq(j),jq(j+1)-1
            i=irow(l)
            y(i)=y(i)+au(l)*x(j)
            y(j)=y(j)-au(l)*x(i)
         enddo
      enddo
!
      return
      end
