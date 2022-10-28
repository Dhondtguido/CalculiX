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
      subroutine forcesolve(zc,nev,a,b,x,eiga,eigb,eigxx,iter,d,
     &  neq,z,istartnmd,iendnmd,nmd,cyclicsymmetry,neqact,
     &  igeneralizedforce)
!
!     solves for the complex eigenfrequencies due to Coriolis 
!     forces
!
      implicit none
!
      logical wantx
!
      integer nev,neq,iter(*),i,j,k,l,istartnmd(*),iendnmd(*),nmd,
     &  cyclicsymmetry,neqact,igeneralizedforce
!
      real*8 z(neq,*),d(*)
!
      complex*16 a(nev,*),b(nev,*),x(nev,*),eiga(*),eigb(*),eigxx(*),
     &  zc(neqact,*)
!
      if(igeneralizedforce.eq.0) then
!
!        no generalized force: multiplication with the eigenmodes
!        is necessary
!
         if(cyclicsymmetry.eq.0) then
            do i=1,nev
               do j=1,nev
                  do k=1,neq
                     a(i,j)=a(i,j)+z(k,i)*zc(k,j)
                  enddo
               enddo
               write(*,*) 
     &              'aerodynamic stiffness/structural stiffness = ',
     &              a(i,i)/d(i)
               a(i,i)=a(i,i)+d(i)
               b(i,i)=(1.d0,0.d0)
            enddo
         else
!     
!     cyclic symmetry
!     
            do l=1,nmd
               do i=istartnmd(l),iendnmd(l)
                  do j=istartnmd(l),iendnmd(l)
                     do k=1,neqact
                        a(i,j)=a(i,j)+z(k,i)*zc(k,j)-
     &                       z(k+neqact,i)*zc(k,j)*(0.d0,1.d0)
                     enddo
                  enddo
                  write(*,*) 
     &                 'aerodynamic stiffness/structural stiffness = ',
     &                 a(i,i)/d(i)
                  a(i,i)=a(i,i)+d(i)
                  b(i,i)=(1.d0,0.d0)
               enddo
            enddo
         endif
      else
!
!        generalized force: the a-matrix is (apart from the diagonal)
!        known
!
         if(cyclicsymmetry.eq.0) then
            do i=1,nev
               write(*,*) 
     &              'aerodynamic stiffness/structural stiffness = ',
     &              a(i,i)/d(i)
               a(i,i)=a(i,i)+d(i)
               b(i,i)=(1.d0,0.d0)
            enddo
         else
!     
!     cyclic symmetry
!     
            do l=1,nmd
               do i=istartnmd(l),iendnmd(l)
                  write(*,*) 
     &                 'aerodynamic stiffness/structural stiffness = ',
     &                 a(i,i)/d(i)
                  a(i,i)=a(i,i)+d(i)
                  b(i,i)=(1.d0,0.d0)
               enddo
            enddo
         endif
      endif
!     
      wantx=.true.
!
!     solving for the complex eigenvalues
!     
      call dlzhes(nev,a,nev,b,nev,x,nev,wantx)
      call dlzit(nev,a,nev,b,nev,x,nev,wantx,iter,eiga,eigb)
!     
      do i=1,nev
         if(iter(i).eq.-1) then
            write(*,*) '*ERROR in forcesolve: fatal error'
            write(*,*) '       in dlzit'
            call exit(201)
         elseif(cdabs(eigb(i)).lt.1.d-10) then
            write(*,*) '*ERROR in forcesolve: eigenvalue'
            write(*,*) '       out of bounds'
            call exit(201)
         else
            eigxx(i)=cdsqrt(eiga(i)/eigb(i))
         endif
      enddo
!     
      return
      end
      
      
