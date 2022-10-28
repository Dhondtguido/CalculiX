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
      subroutine hgstiffness(s,elas,a,gs)
!
!     hourglass control stiffness for 8-node solid mean strain element
!
!     Reference: Flanagan, D.P., Belytschko, T.; "Uniform  strain hexahedron
!     and quadrilateral with orthogonal Hourglass control". Int. J. Num.
!     Meth. Engg., Vol. 17, 679-706, 1981. 
!
!     author: Otto-Ernst Bernhardi
!
      implicit none
!
      integer ii1,jj1,ii,jj,m1
!
      real*8 s(60,60),gs(8,4),a,elas(1),hgls,ahr
!
!
!
      ahr=elas(1)*a
c     write(6,*) "stiffness:", ahr
!
      jj1=1
      do jj=1,8
         ii1=1
         do ii=1,jj
           hgls=0.0d0
           do m1=1,4
              hgls=hgls+gs(jj,m1)*gs(ii,m1)
           enddo
           hgls=hgls*ahr
           s(ii1,jj1)=s(ii1,jj1)+hgls
           s(ii1+1,jj1+1)=s(ii1+1,jj1+1)+hgls
           s(ii1+2,jj1+2)=s(ii1+2,jj1+2)+hgls
           ii1=ii1+3
         enddo
         jj1=jj1+3
      enddo
      return
      end
