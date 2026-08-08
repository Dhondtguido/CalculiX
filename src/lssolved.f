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
C> lssolved
C> ========
C>
C> Solve overdetermined or underdetermined linear equation systems in
C> the least squares sense using the singular value decomposition and
C> double precision.
C>
C> This wraps DGELSS from LAPACK.
C>
C> Notes
C> =====
C>
C> Singular values s(i) are treated as zero if s(i)/s(1) <= cutoff.
C> For cutoff < 0, the cutoff is set to machine precision by DGELSS.
C>
C> a is overwritten with the first min(m,n) right singular
C> vectors, rowwise. b is overwritten with nrhs least squares solutions
C> in its nrhs columns, where the row length of the solutions is n.
C> Refer to the documentation of LAPACK for DGELSS for more information.
C>
C> Authors
C> -------
C> 2026-03-24 Christoph Woelfle - initial commit
C>
C> Parameters
C> ----------
C> @param m number of rows of a
C> @param n number of columns of a
C> @param nrhs number of right hand sides
C> @param a input matrix and on output right singular vectors
C> @param b right hand sides and on output least squares solutions
C> @param cutoff singular value cutoff factor
C> @param ier error code
C>
      subroutine lssolved(m,n,nrhs,a,b,cutoff,ier)
!
      implicit none
!
      integer, intent(in) :: m,n,nrhs
      integer, intent(out) :: ier
      real*8, intent(inout) :: a(m,n),b(m,nrhs),cutoff
!
      integer mx,mn,lwork,rank
      real*8, dimension(:), allocatable :: s,work
!
      ier = 0
      mx = max(m,n)
      mn = min(m,n)
      lwork = 3*mn + max(2*mn,mx,nrhs)
      allocate(s(mn))
      allocate(work(lwork))
!
!     query optimal workspace size
!
      lwork = -1
      call dgelss(m,n,nrhs,a,m,b,m,s,cutoff,rank,work,lwork,ier)
      if (ier .ne. 0) then
         write(*,*)
     &        "ERROR in lssolved.f when computing optimal workspace: ",
     &        ier
         return
      endif
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))
!
!     compute least squares solution
!
      call dgelss(m,n,nrhs,a,m,b,m,s,cutoff,rank,work,lwork,ier)
      if (ier .ne. 0) then
         write(*,*)
     &        "ERROR in lssolved.f when computing solution: ",
     &        ier
         return
      endif
!
      deallocate(s)
      deallocate(work)
!
      end
