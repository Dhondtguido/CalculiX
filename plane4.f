!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998 Guido Dhondt
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
      subroutine plane4(co,node,nodep,a,b,c,d)
!
!     given are 5 nodes: node and nodep(1) up to nodep(4)
!
!     first, the node among the nodep's with the largest
!     distance from node is eliminated. The remaining nodes 
!     are stored into noden(1) up to noden(3).
!
!     Then, the equation of the plane through the
!     nodes noden(1),noden(2) and noden(3) in the form
!     a*x+b*y+c*z+d=0 such that the triangle through the
!     nodes noden(1),noden(2),nopen(3) is numbered clockwise 
!     when looking in the direction of vector (a,b,c)
!
      implicit none
!
      integer nodep(4),i,j,noden(3),node,kflag,idist(4),n
!
      real*8 co(3,*),a,b,c,d,dd,p12(3),p23(3),p31(3),dist(4)
!
      kflag=2
      n=4
!
!     determining the distance of the nodep's to node
!
      do i=1,4
         dist(i)=((co(1,nodep(i))-co(1,node))**2+
     &            (co(2,nodep(i))-co(2,node))**2+
     &            (co(3,nodep(i))-co(3,node))**2)
         idist(i)=nodep(i)
c         write(*,*) nodep(i),dist(i)
      enddo
!
!     sorting the distances
!      
      call dsort(dist,idist,n,kflag)
!
!     storing the 3 closest nodes in noden
!
      j=0
      do i=1,4
         if(nodep(i).eq.idist(4)) cycle
         j=j+1
         noden(j)=nodep(i)
      enddo
c      write(*,*) 'noden ',(noden(i),i=1,3)
!
!     sides of the triangle
!
      do i=1,3
         p12(i)=co(i,noden(2))-co(i,noden(1))
         p23(i)=co(i,noden(3))-co(i,noden(2))
         p31(i)=co(i,noden(1))-co(i,noden(3))
      enddo
!
!     normalized vector normal to the triangle: xn = p12 x p23
!
      a=p12(2)*p23(3)-p12(3)*p23(2)
      b=p12(3)*p23(1)-p12(1)*p23(3)
      c=p12(1)*p23(2)-p12(2)*p23(1)
      dd=dsqrt(a*a+b*b+c*c)
      a=a/dd
      b=b/dd
      c=c/dd
!
!     determining the inhomogeneous term
!
      d=-a*co(1,noden(1))-b*co(2,noden(1))-c*co(3,noden(1))
!
      return
      end

