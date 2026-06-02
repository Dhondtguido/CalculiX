!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2026 Guido Dhondt
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
C> spmv
C> ====
C>
C> Parallelized computation of sparse matrix vector products (SpMV).
C>
C>    A = (ad,au,jq,irow,nzs,n) is an n-by-n sparse matrix
C>    x                         is a dense vector of length n
C>    y = A*x                   is a dense vector of length n
C>
C> Note
C> ----
C> inputformat==0 is assumed for A (symmetric matrix).
C>
C> Authors
C> -------
C> 2026-04-08 Christoph Woelfle - initial commit
C>
C> Parameters
C> ----------
C> @param n dimensions of matrix and vector
C> @param x input vector
C> @param y output vector
C> @param ad diagonal elements of input matrix
C> @param au offdiagonal elements of input matrix
C> @param jq base-1 start indices of the columns of A in au and irow
C> @param irow row/column indices of the offdiagonal elements
C> @param num_cpus number of threads to use
C>
      subroutine spmv(n,x,y,ad,au,jq,irow,num_cpus)
!
      use iso_c_binding
!
      implicit none
!
      integer irow(*),n,j,l,i,jq(n+1),num_cpus
      real*8 y(n),x(n),au(*),ad(n)
!
!$omp parallel num_threads(num_cpus)
!
!$omp do private(i,l), reduction(+:y)
      do j=1,n
        do l=jq(j),jq(j+1)-1
          i=irow(l)
          y(i)=y(i)+au(l)*x(j)
          y(j)=y(j)+au(l)*x(i)
        enddo
      enddo
!
!$omp do
      do i=1,n
         y(i)=y(i)+ad(i)*x(i)
      end do
!$omp end do
!
!$omp end parallel
      end subroutine
