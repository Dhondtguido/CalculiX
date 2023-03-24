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
      real*8 function fuvertex(n,x,cotet,kontet,ipoeln,ieln,node,iedge,
     &     ipoeled,ieled,iedgmid,iedtet)
!     
!     determination of the quality of the ball of linear tetrahedron
!     elements around a node (P.L. George)
!     
      implicit none
!     
      integer n,kontet(4,*),ipoeln(*),ieln(2,*),node,ine(2,6),indexe,
     &     ielem,nodes(4),n1,n2,i,j,iedge,ipoeled(*),ieled(2,*),
     &     iedgmid(*),iedtet(6,*)
!     
      real*8 cotet(3,*),alpha,x(n),volume,surface(4),totsurface,
     &     edgelength(6),radius,hmax,quality
!
      data ine /1,2,2,3,1,3,1,4,2,4,3,4/
!     
      fuvertex=0.d0
!     
!     alpha is the proporionality factor
!     
      alpha=dsqrt(6.d0)/12.d0
!
      indexe=ipoeln(node)
!
      do
        if(indexe.eq.0) exit
! 
!       an element belonging to the ball of node
!
        ielem=ieln(1,indexe)
!     
        do i=1,4
          nodes(i)=kontet(i,ielem)
          if(nodes(i).eq.node) then
            do j=1,3
              cotet(j,node)=x(j)
            enddo
          endif
        enddo
!     
!     calculating the volume of the element
!     
        call calcvol(nodes(1),nodes(2),nodes(3),nodes(4),cotet,
     &       volume)
        if(volume.le.0.d0) volume=1.d-30
!     
!     calculating area of each face in the element
!     
        call calcsurf(nodes(1),nodes(2),nodes(3),cotet,
     &       surface(1))
        call calcsurf(nodes(2),nodes(3),nodes(4),cotet,
     &       surface(2))
        call calcsurf(nodes(3),nodes(4),nodes(1),cotet,
     &       surface(3))
        call calcsurf(nodes(4),nodes(1),nodes(2),cotet,
     &       surface(4))
!     
!     calculating the total surface
!     
        totsurface=surface(1)+surface(2)+surface(3)+surface(4)
!     
!     radius of the inscribed sphere
!     
        radius=3.d0*volume/totsurface
!     
!     length of each edge
!     
        do j=1,6
          n1=nodes(ine(1,j))
          n2=nodes(ine(2,j))
          edgelength(j)=dsqrt((cotet(1,n1)-cotet(1,n2))**2+
     &         (cotet(2,n1)-cotet(2,n2))**2+
     &         (cotet(3,n1)-cotet(3,n2))**2)
        enddo
!     
!     maximum edge length
!     
        hmax=maxval(edgelength)
!     
!     quality of the element
!     
        quality=alpha*hmax/radius
!
!     worst quality of elements treated so far
!
        fuvertex=max(fuvertex,quality)
!
        indexe=ieln(2,indexe)
      enddo
      write(*,100) x(1),x(2),x(3),fuvertex
 100  format('fuvertex ',4(1x,f15.8))
!     
      return
      end
