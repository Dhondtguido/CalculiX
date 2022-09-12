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
!     y=A*x for real sparse matrices (symmetric and non-symmetric)
!
!     storage of the matrix in a:
!        - first the lower triangular terms
!        - then, if the matrix is non-symmetric, the upper triangular terms
!        - finally the diagonal terms
!
      subroutine matvec_struct(n,x,y,nelt,ia,ja,a,isym)
      use omp_lib
!
      implicit none
!
      integer ia(*),ja(*),i,j,l,n,nelt,isym,nd,na
      real*8 y(*),x(*),a(*)
!
!     number of off-diagonal terms
!
      nd=nelt-n
!
      if(isym.eq.0) then
         na=nd/2
!
!        non-symmetric
!
!        diagonal terms
!
c$omp parallel default(none)
c$omp& shared(n,x,a,y,nd,ja,ia,na)
c$omp& private(i,l,j)
c$omp do
         do i=1,n
            y(i)=a(nd+i)*x(i)
         enddo
c$omp end do
!     
!        off-diagonal terms
!     
!        number of upper triangular terms
!
c$omp do
         do j=1,n
            do l=ja(j),ja(j+1)-1
               i=ia(l)
c$omp atomic
               y(i)=y(i)+a(l)*x(j)
            enddo
         enddo
c$omp end do
!
c$omp do
         do j=1,n
            do l=ja(j),ja(j+1)-1
               i=ia(l)
               y(j)=y(j)+a(l+na)*x(i)
            enddo
         enddo
c$omp end do
c$omp end parallel
!
      else
!
!        symmetric
!     
!        diagonal terms
!
c$omp parallel default(none)
c$omp& shared(n,x,a,y,nd,ja,ia)
c$omp& private(i,l,j)
c$omp do
         do i=1,n
            y(i)=a(nd+i)*x(i)
         enddo
c$omp end do
!
!        off-diagonal terms
!     
c$omp do
         do j=1,n
            do l=ja(j),ja(j+1)-1
               i=ia(l)
c$omp atomic
               y(i)=y(i)+a(l)*x(j)
            enddo
         enddo
c$omp end do
!
c$omp do
         do j=1,n
            do l=ja(j),ja(j+1)-1
               i=ia(l)
               y(j)=y(j)+a(l)*x(i)
            enddo
         enddo
c$omp end do
c$omp end parallel
      endif
!
      return
      end
