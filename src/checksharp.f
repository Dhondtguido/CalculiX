!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine checksharp(nexternedg,iedgextfa,cotet,ifacext,isharp)
!     
!     check which edges of the unrefined mesh are sharp
!     
!     the check is done on the angle between the normals on the
!     adjacent faces
!     
!     To this end middle nodes are NOT taken into account, i.e. quadratic
!     faces are reduced to linear faces      
!     
      implicit none
!     
      integer nexternedg,iedgextfa(2,*),ifacext(6,*),isharp(*),iflag,
     &     imastfa,i,j,k
!     
      real*8 cotet(3,*),xi,et,xn1(3),xn2(3),dd,xl(3,3),xsj(3),xs(3,7),
     &     shp(7,3)
!     
!
!     
      iflag=2
      xi=1.d0/3.d0
      et=1.d0/3.d0
!     
      do i=1,nexternedg
!     
!     first neighboring face
!     
        imastfa=iedgextfa(1,i)
        do j=1,3
          do k=1,3
            xl(k,j)=cotet(k,ifacext(j,imastfa))
          enddo
        enddo
        call shape3tri(xi,et,xl,xsj,xs,shp,iflag)
        dd=dsqrt(xsj(1)*xsj(1)+xsj(2)*xsj(2)+xsj(3)*xsj(3))
        do j=1,3
          xn1(j)=xsj(j)/dd
        enddo
!     
!     second neighboring face
!     
        imastfa=iedgextfa(2,i)
        do j=1,3
          do k=1,3
            xl(k,j)=cotet(k,ifacext(j,imastfa))
          enddo
        enddo
        call shape3tri(xi,et,xl,xsj,xs,shp,iflag)
        dd=dsqrt(xsj(1)*xsj(1)+xsj(2)*xsj(2)+xsj(3)*xsj(3))
        do j=1,3
          xn2(j)=xsj(j)/dd
        enddo
!     
!     if the normals are nearly parallel, the edge is no sharp edge
!     "nearly parallel" means that the angle between the vectors
!     is smaller than 0.0464 degrees.
!     
        if(dabs(xn1(1)*xn2(1)+xn1(2)*xn2(2)+xn1(3)*xn2(3)-1.d0)
     &       .gt.1.d-10) isharp(i)=1
      enddo
!     
      return
      end
