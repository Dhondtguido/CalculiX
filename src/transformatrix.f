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
      subroutine transformatrix(xab,p,a)
!
!     determines the transformation matrix a in a point p for a carthesian 
!     (xab(7)>0) or cylindrical transformation (xab(7)<0)
!
!     the matrix a(i,j) corresponds to matrix T, p. 100 in 
!     Dhondt, G., The Finite Element Method for Three-Dimensional
!                 Thermomechanical Applications, Wiley (2004).
!
      implicit none
!
      integer j
!
      real*8 xab(7),p(3),a(3,3),e1(3),e2(3),e3(3),dd
!
      if(xab(7).gt.0) then
!
!        carthesian transformation
!
         e1(1)=xab(1)
         e1(2)=xab(2)
         e1(3)=xab(3)
!
         e2(1)=xab(4)
         e2(2)=xab(5)
         e2(3)=xab(6)
!
         dd=dsqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
         do j=1,3
            e1(j)=e1(j)/dd
         enddo
!
         dd=e1(1)*e2(1)+e1(2)*e2(2)+e1(3)*e2(3)
         do j=1,3
            e2(j)=e2(j)-dd*e1(j)
         enddo
!
         dd=dsqrt(e2(1)*e2(1)+e2(2)*e2(2)+e2(3)*e2(3))
         do j=1,3
            e2(j)=e2(j)/dd
         enddo
!
         e3(1)=e1(2)*e2(3)-e2(2)*e1(3)
         e3(2)=e1(3)*e2(1)-e1(1)*e2(3)
         e3(3)=e1(1)*e2(2)-e2(1)*e1(2)
!
      else
!
!        cylindrical coordinate system in point p
!
         e1(1)=p(1)-xab(1)
         e1(2)=p(2)-xab(2)
         e1(3)=p(3)-xab(3)
!
         e3(1)=xab(4)-xab(1)
         e3(2)=xab(5)-xab(2)
         e3(3)=xab(6)-xab(3)
!
         dd=dsqrt(e3(1)*e3(1)+e3(2)*e3(2)+e3(3)*e3(3))
!
         do j=1,3
            e3(j)=e3(j)/dd
         enddo
!
         dd=e1(1)*e3(1)+e1(2)*e3(2)+e1(3)*e3(3)
!
         do j=1,3
            e1(j)=e1(j)-dd*e3(j)
         enddo
!
         dd=dsqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
!
!        check whether p belongs to the cylindrical coordinate axis
!        if so, an arbitrary vector perpendicular to the axis can
!        be taken
!
         if(dd.lt.1.d-10) then
c            write(*,*) '*WARNING in transformatrix: point on axis'
            if(dabs(e3(1)).gt.1.d-10) then
               e1(2)=1.d0
               e1(3)=0.d0
               e1(1)=-e3(2)/e3(1)
            elseif(dabs(e3(2)).gt.1.d-10) then
               e1(3)=1.d0
               e1(1)=0.d0
               e1(2)=-e3(3)/e3(2)
            else
               e1(1)=1.d0
               e1(2)=0.d0
               e1(3)=-e3(1)/e3(3)
            endif
            dd=dsqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
         endif
!
         do j=1,3
            e1(j)=e1(j)/dd
         enddo
!
         e2(1)=e3(2)*e1(3)-e1(2)*e3(3)
         e2(2)=e3(3)*e1(1)-e1(3)*e3(1)
         e2(3)=e3(1)*e1(2)-e1(1)*e3(2)
!
      endif
!
!     finding the transformation matrix
!
      do j=1,3
         a(j,1)=e1(j)
         a(j,2)=e2(j)
         a(j,3)=e3(j)
      enddo
!
      return
      end













