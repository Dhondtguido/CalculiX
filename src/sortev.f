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
!     ADDITIONAL INPUT Parameters:
!     ipev:     Position of Eigenvalues, saves original Position of 
!            Eigenvalues before sorting
!     eigxr:   Real Part of Eigenvalues out of eigxx, used for
!            sorting Eigenvalues in increasing order
!
      subroutine sortev(nev,nmd,eigxx,cyclicsymmetry,x,eigxr,ipev,
     &     istartnmd,iendnmd,a,b,nevcomplex)
!     
!     sorts the eigenvalues and eigenvectors of complex frequency
!     
      implicit none
!     
      integer nev,i,j,k,l,m,istartnmd(*),iendnmd(*),nmd,
     &     cyclicsymmetry,ipev(*),nevcomplex
!     
      real*8 eigxr(*)
!     
      complex*16 eigxx(*),a(nev),b(nev,*),x(nev,*)
!     
      if (cyclicsymmetry.eq.0)then
!     
!     sorting the eigenvalues according to their size
!     
         do i=1,nevcomplex
            ipev(i)=i
            eigxr(i)=cdabs(eigxx(i))
         enddo
         call dsort(eigxr,ipev,nevcomplex,2)
!     
!     sorting the eigenvectors
!     
         do i=1,nevcomplex
            a(i)=eigxx(ipev(i))
            do j=1,nev
               b(j,i)=x(j,ipev(i))
            enddo
         enddo
!     
!     copying in the original fields
!     
         do i=1,nevcomplex
            eigxx(i)=a(i)
            do j=1,nev
               x(j,i)=b(j,i)
               enddo
            enddo
         else
!     
!     Cyclic Symmetry
!     
            do l=1,nmd
!     
!     sorting the eigenvalues according to their size
!     
!     
            do i=istartnmd(l),iendnmd(l)
               if (l.eq.1) then
                  ipev(i)=i
                  eigxr(i)=cdabs(eigxx(i))
                  k=i
               else
                  k=i-istartnmd(l)+1
                  ipev(k)=i
                  eigxr(i)=cdabs(eigxx(i))
               endif
            enddo
            call dsort (eigxr,ipev,k,2)
!     
!     sorting the eigenvectors
!     
            do i=istartnmd(l),iendnmd(l)
               if (l.eq.1) then
               m=ipev(i)
               a(i)=eigxx(m)
               do j=istartnmd(l),iendnmd(l)
                  b(j,i)=x(j,m)
               enddo
            else
               k=i-istartnmd(l)+1
               a(i)=eigxx(ipev(k))
               do j=istartnmd(l),iendnmd(l)
                  b(j,i)=x(j,m)
               enddo
            endif
            enddo
         enddo
!     
!     copying in the original fields
!     
         do l=1,nmd
            if((a(istartnmd(l)).ne.0).and.
     &           (b(istartnmd(l),istartnmd(l)).ne.0))then
               do i=istartnmd(l),iendnmd(l)
               eigxx(i)=a(i)
               do j=istartnmd(l),iendnmd(l)
                  x(i,j)=b(i,j)
               enddo
            enddo
         endif
      enddo   
      endif
!     
      return
      end
      
