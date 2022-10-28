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
      subroutine hgforce (fn,elas,a,gs,vl,mi,konl)
!
!     hourglass control forces for 8-node solid mean strain element
!
!     Reference: Flanagan, D.P., Belytschko, T.; "Uniform  strain hexahedron
!     and quadrilateral with orthogonal Hourglass control". Int. J. Num.
!     Meth. Engg., Vol. 17, 679-706, 1981. 
!
!     author: Otto-Ernst Bernhardi
!
      implicit none
      integer i,j,k,mi(*),konl(20)
      real*8 gs(8,4),a,elas(1),ahr
      real*8 vl(0:mi(2),20), fn(0:mi(2),*)
      real*8 hglf(3,4)
!
      ahr=elas(1)*a
c     write(6,*) "force:", ahr
!
      do i=1,3
         do k=1,4    
            hglf(i,k)=0.0d0
            do j=1,8
               hglf(i,k)=hglf(i,k)+gs(j,k)*vl(i,j)
            enddo
            hglf(i,k)=hglf(i,k)*ahr
         enddo
      enddo
      do i=1,3
         do j=1,8
            do k=1,4
               fn(i,konl(j))=fn(i,konl(j))+hglf(i,k)*gs(j,k)
            enddo
         enddo
      enddo
      return
      end
