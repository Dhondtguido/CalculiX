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
      subroutine stiff2mat(elas,ckl,vj,cauchy)
!
!     elas(21):   stiffness constants in the spatial description, i.e.
!                 the derivative of the Cauchy stress or the Kirchhoff
!                 stress with respect to the Eulerian strain
!     ckl(3,3):   inverse deformation gradient
!     vj:         Jacobian determinant
!     cauchy:     if 1: elas is written in terms of Cauchy stress
!                 if 0: elas is written in terms of Kirchhoff stress
!
!     OUTPUT:
!
!     elas(21):   stiffness constants in the material description,i.e.
!                 the derivative of the second Piola-Kirchhoff stress (PK2)
!                 with respect to the Lagrangian strain
!
      implicit none
!
      integer cauchy
!
      integer kk(84),i,nt,k,l,m,n
!
      real*8 elas(21),e(21),ckl(3,3),vj
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
         e(i)=elas(1)*ckl(k,1)*ckl(l,1)*ckl(m,1)*ckl(n,1)
     &       +elas(2)*(ckl(k,2)*ckl(l,2)*ckl(m,1)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,1)*ckl(m,2)*ckl(n,2))
     &       +elas(3)*ckl(k,2)*ckl(l,2)*ckl(m,2)*ckl(n,2)
     &       +elas(4)*(ckl(k,3)*ckl(l,3)*ckl(m,1)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,1)*ckl(m,3)*ckl(n,3))
     &       +elas(5)*(ckl(k,3)*ckl(l,3)*ckl(m,2)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,2)*ckl(m,3)*ckl(n,3))
     &       +elas(6)*ckl(k,3)*ckl(l,3)*ckl(m,3)*ckl(n,3)
     &       +elas(7)*(ckl(k,2)*ckl(l,1)*ckl(m,1)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,2)*ckl(m,1)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,1)*ckl(m,2)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,1)*ckl(m,1)*ckl(n,2))
     &       +elas(8)*(ckl(k,2)*ckl(l,2)*ckl(m,2)*ckl(n,1)+
     &                 ckl(k,2)*ckl(l,2)*ckl(m,1)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,1)*ckl(m,2)*ckl(n,2)+
     &                 ckl(k,1)*ckl(l,2)*ckl(m,2)*ckl(n,2))
     &       +elas(9)*(ckl(k,3)*ckl(l,3)*ckl(m,2)*ckl(n,1)+
     &                 ckl(k,3)*ckl(l,3)*ckl(m,1)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,1)*ckl(m,3)*ckl(n,3)+
     &                 ckl(k,1)*ckl(l,2)*ckl(m,3)*ckl(n,3))
     &       +elas(10)*(ckl(k,2)*ckl(l,1)*ckl(m,2)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,2)*ckl(m,2)*ckl(n,1)+
     &                 ckl(k,2)*ckl(l,1)*ckl(m,1)*ckl(n,2)+
     &                 ckl(k,1)*ckl(l,2)*ckl(m,1)*ckl(n,2))
     &       +elas(11)*(ckl(k,3)*ckl(l,1)*ckl(m,1)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,3)*ckl(m,1)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,1)*ckl(m,3)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,1)*ckl(m,1)*ckl(n,3))
     &       +elas(12)*(ckl(k,2)*ckl(l,2)*ckl(m,3)*ckl(n,1)+
     &                 ckl(k,3)*ckl(l,1)*ckl(m,2)*ckl(n,2)+
     &                 ckl(k,1)*ckl(l,3)*ckl(m,2)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,2)*ckl(m,3)*ckl(n,1))
     &       +elas(13)*(ckl(k,3)*ckl(l,3)*ckl(m,3)*ckl(n,1)+
     &                 ckl(k,3)*ckl(l,3)*ckl(m,1)*ckl(n,3)+
     &                 ckl(k,3)*ckl(l,1)*ckl(m,3)*ckl(n,3)+
     &                 ckl(k,1)*ckl(l,3)*ckl(m,3)*ckl(n,3))
     &       +elas(14)*(ckl(k,3)*ckl(l,1)*ckl(m,2)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,3)*ckl(m,2)*ckl(n,1)+
     &                 ckl(k,2)*ckl(l,1)*ckl(m,3)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,2)*ckl(m,3)*ckl(n,1)+
     &                 ckl(k,3)*ckl(l,1)*ckl(m,1)*ckl(n,2)+
     &                 ckl(k,1)*ckl(l,3)*ckl(m,1)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,1)*ckl(m,1)*ckl(n,3)+
     &                 ckl(k,1)*ckl(l,2)*ckl(m,1)*ckl(n,3))
     &       +elas(15)*(ckl(k,3)*ckl(l,1)*ckl(m,3)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,3)*ckl(m,3)*ckl(n,1)+
     &                 ckl(k,3)*ckl(l,1)*ckl(m,1)*ckl(n,3)+
     &                 ckl(k,1)*ckl(l,3)*ckl(m,1)*ckl(n,3))
     &       +elas(16)*(ckl(k,3)*ckl(l,2)*ckl(m,1)*ckl(n,1)+
     &                 ckl(k,2)*ckl(l,3)*ckl(m,1)*ckl(n,1)+
     &                 ckl(k,1)*ckl(l,1)*ckl(m,3)*ckl(n,2)+
     &                 ckl(k,1)*ckl(l,1)*ckl(m,2)*ckl(n,3))
     &       +elas(17)*(ckl(k,3)*ckl(l,2)*ckl(m,2)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,3)*ckl(m,2)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,2)*ckl(m,3)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,2)*ckl(m,2)*ckl(n,3))
     &       +elas(18)*(ckl(k,3)*ckl(l,3)*ckl(m,3)*ckl(n,2)+
     &                 ckl(k,3)*ckl(l,3)*ckl(m,2)*ckl(n,3)+
     &                 ckl(k,3)*ckl(l,2)*ckl(m,3)*ckl(n,3)+
     &                 ckl(k,2)*ckl(l,3)*ckl(m,3)*ckl(n,3))
     &       +elas(19)*(ckl(k,3)*ckl(l,2)*ckl(m,2)*ckl(n,1)+
     &                 ckl(k,2)*ckl(l,3)*ckl(m,2)*ckl(n,1)+
     &                 ckl(k,3)*ckl(l,2)*ckl(m,1)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,3)*ckl(m,1)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,1)*ckl(m,3)*ckl(n,2)+
     &                 ckl(k,1)*ckl(l,2)*ckl(m,3)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,1)*ckl(m,2)*ckl(n,3)+
     &                 ckl(k,1)*ckl(l,2)*ckl(m,2)*ckl(n,3))
     &       +elas(20)*(ckl(k,3)*ckl(l,2)*ckl(m,3)*ckl(n,1)+
     &                 ckl(k,2)*ckl(l,3)*ckl(m,3)*ckl(n,1)+
     &                 ckl(k,3)*ckl(l,1)*ckl(m,3)*ckl(n,2)+
     &                 ckl(k,1)*ckl(l,3)*ckl(m,3)*ckl(n,2)+
     &                 ckl(k,3)*ckl(l,2)*ckl(m,1)*ckl(n,3)+
     &                 ckl(k,2)*ckl(l,3)*ckl(m,1)*ckl(n,3)+
     &                 ckl(k,3)*ckl(l,1)*ckl(m,2)*ckl(n,3)+
     &                 ckl(k,1)*ckl(l,3)*ckl(m,2)*ckl(n,3))
     &       +elas(21)*(ckl(k,3)*ckl(l,2)*ckl(m,3)*ckl(n,2)+
     &                 ckl(k,2)*ckl(l,3)*ckl(m,3)*ckl(n,2)+
     &                 ckl(k,3)*ckl(l,2)*ckl(m,2)*ckl(n,3)+
     &                 ckl(k,2)*ckl(l,3)*ckl(m,2)*ckl(n,3))
      enddo
!
      if(cauchy.eq.1) then
         do i=1,21
            elas(i)=e(i)*vj
         enddo
      else
         do i=1,21
            elas(i)=e(i)
         enddo
      endif
!
      return
      end
