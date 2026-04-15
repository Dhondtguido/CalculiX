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
C> eigsymd
C> =======
C>
C> Compute the eigenvalues and (if requested) the eigenvectors of a
C> double precision symmetric matrix.
C>
C> This wraps DSYEV from LAPACK.
C>
C> Notes
C> -----
C> Does NOT support strided a.
C>
C> Authors
C> -------
C> 2026-03-12 Christoph Woelfle - initial commit
C>
C> Parameters
C> ----------
C> @param n row and column dimension of a
C> @param a input matrix
C> @param w eigenvalues in ascending order
C> @param matz (0 if only eigenvalues requested, else 1)
C> @param z eigenvectors of a if requested
C> @param ier error code
C>
      subroutine eigsymd(n,a,w,matz,z,ier)
!
      implicit none
      integer n, matz, ier, lwork
      real*8 a(n,n), w(n), z(n,n)
      real*8, dimension(:), allocatable :: work
!
      ier = 0
!
      lwork=max(1,3*n-1)
      allocate(work(lwork))
!
!     DSYEV might overwrite a with eigenvalues, work with copy of a in z
!
      if (matz .ne. 0) then
         z = a
      end if
!
!     query optimal workspace size for larger matrices
!
      if (n .gt. 3) then
         lwork = -1
         if (matz .eq. 0) then
            call dsyev('V','U',n,a,n,w,work,lwork,ier)
         else
            call dsyev('V','U',n,z,n,w,work,lwork,ier)
         endif
         if (ier .ne. 0) then
            write(*,*)
     &           "ERROR in eigsymd.f when computing",
     &           "optimal work array: ",
     &           ier
            return
         endif
         lwork = int(work(1))
         deallocate(work)
         allocate(work(lwork))
      endif
!
!     solve the eigenproblem
!
      if (matz .eq. 0) then
         call dsyev('V','U',n,a,n,w,work,lwork,ier)
      else
         call dsyev('V','U',n,z,n,w,work,lwork,ier)
      endif
      deallocate(work)
      if (ier .ne. 0) then
         write(*,*)
     &        "ERROR in eigsymd.f when solving eigenproblem: ",
     &        ier
         return
      endif
!
      end
