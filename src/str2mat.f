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
      subroutine str2mat(str,ckl,vj,cauchy)
!
!     converts the stress in spatial coordinates into material coordinates 
!     or the strain in material coordinates into spatial coordinates. 
!
!     INPUT:
!
!     str(6):     Cauchy stress, Kirchhoff stress or Lagrange strain
!                 component order: 11,22,33,12,13,23
!     ckl(3,3):   the inverse deformation gradient
!     vj:         Jakobian determinant
!     cauchy:     integer variable
!                 if 1: str contains the Cauchy stress
!                 if 0: str contains the Kirchhoff stress or
!                           Lagrange strain
!
!     OUTPUT:
!
!     str(6):     Piola-Kirchhoff stress of the second kind (PK2) or
!                 Euler strain
!
      implicit none
!
      integer cauchy
!
      integer i,m1,m2
!
      real*8 str(6),s(6),ckl(3,3),vj
!
      do i=1,6
         if(i.eq.1) then
            m1=1
            m2=1
         elseif(i.eq.2) then
            m1=2
            m2=2
         elseif(i.eq.3) then
            m1=3
            m2=3
         elseif(i.eq.4) then
            m1=2
            m2=1
         elseif(i.eq.5) then
            m1=3
            m2=1
         else
            m1=3
            m2=2
         endif
!
         s(i)=(str(1)*ckl(m1,1)*ckl(m2,1)+
     &        str(2)*ckl(m1,2)*ckl(m2,2)+
     &        str(3)*ckl(m1,3)*ckl(m2,3)+
     &        str(4)*(ckl(m1,1)*ckl(m2,2)+ckl(m1,2)*ckl(m2,1))+
     &        str(5)*(ckl(m1,1)*ckl(m2,3)+ckl(m1,3)*ckl(m2,1))+
     &        str(6)*(ckl(m1,2)*ckl(m2,3)+ckl(m1,3)*ckl(m2,2)))
!
      enddo
!
      if(cauchy.eq.1) then
         do i=1,6
            str(i)=s(i)*vj
         enddo
      else
         do i=1,6
            str(i)=s(i)
         enddo
      endif
!
      return
      end
