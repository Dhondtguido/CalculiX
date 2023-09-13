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
      subroutine stiff2mat(stiff,ckl,vj,cauchy)
!
!     stiff(21):   stiffness constants in the spatial description, i.e.
!                 the derivative of the Cauchy stress or the Kirchhoff
!                 stress with respect to the Eulerian strain
!     ckl(3,3):   inverse deformation gradient
!     vj:         Jacobian determinant
!     cauchy:     if 1: stiff is written in terms of Cauchy stress
!                 if 0: stiff is written in terms of Kirchhoff stress
!
!     OUTPUT:
!
!     stiff(21):   stiffness constants in the material description,i.e.
!                 the derivative of the second Piola-Kirchhoff stress (PK2)
!                 with respect to the Lagrangian strain
!
      implicit none
!
      integer cauchy
!
      integer kk(84),i,nt,k,l,m,n
!
      real*8 stiff(21),e(21),ckl(3,3),vj
!
      data kk /1,1,1,1,1,1,2,2,2,2,2,2,1,1,3,3,2,2,3,3,3,3,3,3,
     &  1,1,1,2,2,2,1,2,3,3,1,2,1,2,1,2,1,1,1,3,2,2,1,3,3,3,1,3,
     &  1,2,1,3,1,3,1,3,1,1,2,3,2,2,2,3,3,3,2,3,1,2,2,3,1,3,2,3,
     &  2,3,2,3/
!
      nt=0
      do i=1,21
         k=kk(nt+1)
         l=kk(nt+2)
         m=kk(nt+3)
         n=kk(nt+4)
         nt=nt+4
         e(i)=stiff(1)*ckl(k,1)*ckl(l,1)*ckl(m,1)*ckl(n,1)
     &       +stiff(2)*(ckl(k,2)*ckl(l,2)*ckl(m,1)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,1)*ckl(m,2)*ckl(n,2))
     &       +stiff(3)*ckl(k,2)*ckl(l,2)*ckl(m,2)*ckl(n,2)
     &       +stiff(4)*(ckl(k,3)*ckl(l,3)*ckl(m,1)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,1)*ckl(m,3)*ckl(n,3))
     &       +stiff(5)*(ckl(k,3)*ckl(l,3)*ckl(m,2)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,2)*ckl(m,3)*ckl(n,3))
     &       +stiff(6)*ckl(k,3)*ckl(l,3)*ckl(m,3)*ckl(n,3)
     &       +stiff(7)*(ckl(k,2)*ckl(l,1)*ckl(m,1)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,2)*ckl(m,1)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,1)*ckl(m,2)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,1)*ckl(m,1)*ckl(n,2))
     &       +stiff(8)*(ckl(k,2)*ckl(l,2)*ckl(m,2)*ckl(n,1)+
     &                 ckl(k,2)*ckl(l,2)*ckl(m,1)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,1)*ckl(m,2)*ckl(n,2)+
     &                 ckl(k,1)*ckl(l,2)*ckl(m,2)*ckl(n,2))
     &       +stiff(9)*(ckl(k,3)*ckl(l,3)*ckl(m,2)*ckl(n,1)+
     &                 ckl(k,3)*ckl(l,3)*ckl(m,1)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,1)*ckl(m,3)*ckl(n,3)+
     &                 ckl(k,1)*ckl(l,2)*ckl(m,3)*ckl(n,3))
     &       +stiff(10)*(ckl(k,2)*ckl(l,1)*ckl(m,2)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,2)*ckl(m,2)*ckl(n,1)+
     &                 ckl(k,2)*ckl(l,1)*ckl(m,1)*ckl(n,2)+
     &                 ckl(k,1)*ckl(l,2)*ckl(m,1)*ckl(n,2))
     &       +stiff(11)*(ckl(k,3)*ckl(l,1)*ckl(m,1)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,3)*ckl(m,1)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,1)*ckl(m,3)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,1)*ckl(m,1)*ckl(n,3))
     &       +stiff(12)*(ckl(k,2)*ckl(l,2)*ckl(m,3)*ckl(n,1)+
     &                 ckl(k,3)*ckl(l,1)*ckl(m,2)*ckl(n,2)+
     &                 ckl(k,1)*ckl(l,3)*ckl(m,2)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,2)*ckl(m,3)*ckl(n,1))
     &       +stiff(13)*(ckl(k,3)*ckl(l,3)*ckl(m,3)*ckl(n,1)+
     &                 ckl(k,3)*ckl(l,3)*ckl(m,1)*ckl(n,3)+
     &                 ckl(k,3)*ckl(l,1)*ckl(m,3)*ckl(n,3)+
     &                 ckl(k,1)*ckl(l,3)*ckl(m,3)*ckl(n,3))
     &       +stiff(14)*(ckl(k,3)*ckl(l,1)*ckl(m,2)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,3)*ckl(m,2)*ckl(n,1)+
     &                 ckl(k,2)*ckl(l,1)*ckl(m,3)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,2)*ckl(m,3)*ckl(n,1)+
     &                 ckl(k,3)*ckl(l,1)*ckl(m,1)*ckl(n,2)+
     &                 ckl(k,1)*ckl(l,3)*ckl(m,1)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,1)*ckl(m,1)*ckl(n,3)+
     &                 ckl(k,1)*ckl(l,2)*ckl(m,1)*ckl(n,3))
     &       +stiff(15)*(ckl(k,3)*ckl(l,1)*ckl(m,3)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,3)*ckl(m,3)*ckl(n,1)+
     &                 ckl(k,3)*ckl(l,1)*ckl(m,1)*ckl(n,3)+
     &                 ckl(k,1)*ckl(l,3)*ckl(m,1)*ckl(n,3))
     &       +stiff(16)*(ckl(k,3)*ckl(l,2)*ckl(m,1)*ckl(n,1)+
     &                 ckl(k,2)*ckl(l,3)*ckl(m,1)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,1)*ckl(m,3)*ckl(n,2)+
     &                 ckl(k,1)*ckl(l,1)*ckl(m,2)*ckl(n,3))
     &       +stiff(17)*(ckl(k,3)*ckl(l,2)*ckl(m,2)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,3)*ckl(m,2)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,2)*ckl(m,3)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,2)*ckl(m,2)*ckl(n,3))
     &       +stiff(18)*(ckl(k,3)*ckl(l,3)*ckl(m,3)*ckl(n,2)+
     &                 ckl(k,3)*ckl(l,3)*ckl(m,2)*ckl(n,3)+
     &                 ckl(k,3)*ckl(l,2)*ckl(m,3)*ckl(n,3)+
     &                 ckl(k,2)*ckl(l,3)*ckl(m,3)*ckl(n,3))
     &       +stiff(19)*(ckl(k,3)*ckl(l,2)*ckl(m,2)*ckl(n,1)+
     &                 ckl(k,2)*ckl(l,3)*ckl(m,2)*ckl(n,1)+
     &                 ckl(k,3)*ckl(l,2)*ckl(m,1)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,3)*ckl(m,1)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,1)*ckl(m,3)*ckl(n,2)+
     &                 ckl(k,1)*ckl(l,2)*ckl(m,3)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,1)*ckl(m,2)*ckl(n,3)+
     &                 ckl(k,1)*ckl(l,2)*ckl(m,2)*ckl(n,3))
     &       +stiff(20)*(ckl(k,3)*ckl(l,2)*ckl(m,3)*ckl(n,1)+
     &                 ckl(k,2)*ckl(l,3)*ckl(m,3)*ckl(n,1)+
     &                 ckl(k,3)*ckl(l,1)*ckl(m,3)*ckl(n,2)+
     &                 ckl(k,1)*ckl(l,3)*ckl(m,3)*ckl(n,2)+
     &                 ckl(k,3)*ckl(l,2)*ckl(m,1)*ckl(n,3)+
     &                 ckl(k,2)*ckl(l,3)*ckl(m,1)*ckl(n,3)+
     &                 ckl(k,3)*ckl(l,1)*ckl(m,2)*ckl(n,3)+
     &                 ckl(k,1)*ckl(l,3)*ckl(m,2)*ckl(n,3))
     &       +stiff(21)*(ckl(k,3)*ckl(l,2)*ckl(m,3)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,3)*ckl(m,3)*ckl(n,2)+
     &                 ckl(k,3)*ckl(l,2)*ckl(m,2)*ckl(n,3)+
     &                 ckl(k,2)*ckl(l,3)*ckl(m,2)*ckl(n,3))
      enddo
!
      if(cauchy.eq.1) then
         do i=1,21
            stiff(i)=e(i)*vj
         enddo
      else
         do i=1,21
            stiff(i)=e(i)
         enddo
      endif
!
      return
      end
