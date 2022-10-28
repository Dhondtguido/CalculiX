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
!      You should have received a copy of the GNU General Public License
!      along with this program; if not, write to the Free Software
!      Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!     
!      INPUT Parameters:
!      neq:      number of equations, length of eigenvectors 
!                (non-cyclic-symmetry)
!      nev:      number of eigenvectors
!      neqact:   length of eigenvectors (cyclic symmetry)
!      z:        Matrix of eigenvectors out of frequency analysis
!                storage: for each eigenvalue: first complete real part,
!                         then complete imaginary part of eigenvector
!                         (imaginary part only for cyclic symmetry)
!      zz:       Matrix of eigenvectors out of complex frequency analysis
!                storage: for each eigenvalue: first complete real part,
!                         then complete imaginary part of eigenvector
!      istartnmd:first number of nev of nodal diameter l
!      iendnmd:  last number ov nev of nodal diameter l
!      nmd:      Number of nodal diameters (cyclic symmetry)
!      
!      OUTPUT Parameters
!      xmac:     Matrix that contains absolute MAC-values of all vectors 
!                combined
!                with all vectors; in case of cyclic symmetry: only eigenvectors
!                of the same nodal diameter are taken in account
!      xmaccpx:  Matrix that contains the complex MAC-values of all vectors 
!                combined with all vectors; in case of cyclic symmetry: 
!                only eigenvectors of the same nodal diameter are taken 
!                in account
!
      subroutine calcmac(neq,z,zz,nev,xmac,xmaccpx,istartnmd,
     &     iendnmd,nmd,cyclicsymmetry,neqact,bett,betm,nevcomplex)
!     
!     calculates the Modal Assurance Criterium MAC=<z,zz>/(||z||*||zz||)
!     
      implicit none
!     
      integer neq,nev,i,j,k,l,istartnmd(*),iendnmd(*),nmd,
     &     cyclicsymmetry,neqact,nevcomplex
!     
      real*8 bett(*),betm(*),xmac(nev,*),z(neq,*),zz(2*neqact,*)
!     
      complex*16 xmaccpx(nev,*)
!     
      if(cyclicsymmetry.eq.0)then
!     size of vectors
         do i=1,nev
            do k=1,neq
               bett(i)=bett(i)+z(k,i)**2
            enddo
         enddo
         do i=1,nevcomplex
            do k=1,neq
               betm(i)=betm(i)+zz(k,i)**2+zz(k+neq,i)**2
            enddo
         enddo
!
         do i=1,nev
            bett(i)=dsqrt(bett(i))
         enddo
         do i=1,nevcomplex
            betm(i)=dsqrt(betm(i)) 
         enddo
!     Calculation of MAC
         do i=1,nev
            do j=1,nevcomplex
               do k=1,neq  
                  xmaccpx(i,j)=xmaccpx(i,j)+z(k,i)
     &                 *(zz(k,j)+zz(k+neq,j)*(0.d0,1.d0))
               enddo
               xmac(i,j)=cdabs(xmaccpx(i,j))
               xmac(i,j)=xmac(i,j)/bett(i)/betm(j)
            enddo
         enddo
!     
!     Cyclic Symmetry
!     size of vectors
!     
      else
         do i=1,nev
            do k=1,neqact
               bett(i)=bett(i)+z(k,i)**2+z(k+neqact,i)**2
               betm(i)=betm(i)+zz(k,i)**2+zz(k+neqact,i)**2
            enddo
            bett(i)=dsqrt(bett(i))
            betm(i)=dsqrt(betm(i))
         enddo
!     Calculation of MAC
         do l=1,nmd
            do i=istartnmd(l),iendnmd(l)
               do j=istartnmd(l),iendnmd(l)
                  xmac(i,j)=0
                  do k=1,neqact,2
                     xmaccpx(i,j)=xmaccpx(i,j)+
     &                    (z(k,i)-z(k+neqact,i)*(0.d0,1.d0))*
     &                    (zz(k,j)+zz(k+neqact,j)*(0.d0,1.d0))
                  enddo
                  xmac(i,j)=cdabs(xmaccpx(i,j))
                  xmac(i,j)=(xmac(i,j))/bett(i)/betm(j)
               enddo
            enddo
         enddo
      endif
!     
      return
      end
      
