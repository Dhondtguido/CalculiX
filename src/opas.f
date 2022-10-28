!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
!     y=A*x for real sparse asymmetric matrices
!     
!     storage of the matrix:
!     au: first lower triangle, then upper triangle
!     ad: diagonal terms
!     
      subroutine opas(n,x,y,ad,au,jq,irow,nzs)
!     
      implicit none
!     
      integer irow(*),n,nzs(3),j,l,i,jq(*)
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
          y(j)=y(j)+au(l+nzs(3))*x(i)
        enddo
      enddo
!     
      return
      end
