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
      real*8 function fumid(n,x,cotet,kontet,ipoeln,ieln,node,iedge,
     &     ipoeled,ieled,iedgmid,iedtet)
!     
!     determination of the quality of the shell of quadratic tetrahedron
!     elements around a node; consists of a combination of:
!     - the quality of the linear tetrahedrons in which it can be
!       subdivided (P.L. George)
!     - a penalty function for negative Jacobians at the integration
!       points
!     
      implicit none
!     
      integer n,kontet(4,*),ipoeln(*),ieln(2,*),node,ine(2,6),indexe,
     &     ielem,nodes(10),n1,n2,i,j,k,m,iedge,ipoeled(*),ieled(2,*),
     &     iedgmid(*),iedtet(6,*),i6(4,8),nodesedge(4),iflag
!     
      real*8 cotet(3,*),alpha,x(n),volume,surface(4),totsurface,
     &     edgelength(6),radius,hmax,quality,xi,et,ze,shp(4,10),xsj,
     &     xl(3,10)
!
!     nodes belonging to the linear sub-tetrahedra
!
      data i6 /8,9,10,4,1,5,7,8,7,6,3,10,9,8,10,7,
     &     8,9,5,7,9,10,6,7,5,6,7,9,5,2,6,9/
!
!     edges of a tetrahedron
!
      data ine /1,2,2,3,1,3,1,4,2,4,3,4/
!
      data iflag /2/
!
      include "gauss.f"
!      
      fumid=0.d0
!     
!     alpha is the proporionality factor
!     
      alpha=dsqrt(6.d0)/12.d0
!
      indexe=ipoeled(iedge)
!
      do
        if(indexe.eq.0) exit
        ielem=ieled(1,indexe)
!     
!     determine the coordinates of the nodes belonging to the element
!     
!     vertex nodes
!     
        do i=1,4
          nodes(i)=kontet(i,ielem)
          do k=1,3
            xl(k,i)=cotet(k,nodes(i))
          enddo
        enddo
!     
!     middle nodes
!     
        do i=1,6
          nodes(i+4)=iedgmid(iedtet(i,ielem))
          do k=1,3
            xl(k,i+4)=cotet(k,nodes(i+4))
          enddo
!     
!     replace the coordinates by the optimization variables
!     
          if(nodes(i+4).eq.node) then
            do m=1,3
              cotet(m,node)=x(m)
              xl(m,i+4)=x(m)
            enddo
          endif
        enddo
!     
!     loop over the linear sub-tetrahedrons
!     
        do i=1,8
!     
          call calcvol(nodes(i6(1,i)),nodes(i6(2,i)),nodes(i6(3,i)),
     &         nodes(i6(4,i)),cotet,volume)
          if(volume.le.1.d-15) volume=1.d-15
!     
!     calculating the area of each face of the element
!     
          call calcsurf(nodes(i6(1,i)),nodes(i6(2,i)),nodes(i6(3,i)),
     &         cotet,surface(1))
          call calcsurf(nodes(i6(2,i)),nodes(i6(3,i)),nodes(i6(4,i)),
     &         cotet,surface(2))
          call calcsurf(nodes(i6(3,i)),nodes(i6(4,i)),nodes(i6(1,i)),
     &         cotet,surface(3))
          call calcsurf(nodes(i6(4,i)),nodes(i6(1,i)),nodes(i6(2,i)),
     &         cotet,surface(4))
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
          nodesedge(1)=nodes(i6(1,i))
          nodesedge(2)=nodes(i6(2,i))
          nodesedge(3)=nodes(i6(3,i))
          nodesedge(4)=nodes(i6(4,i))
!     
          do j=1,6
            n1=nodesedge(ine(1,j))
            n2=nodesedge(ine(2,j))
            edgelength(j)=dsqrt((cotet(1,n1)-cotet(1,n2))**2+
     &           (cotet(2,n1)-cotet(2,n2))**2+
     &           (cotet(3,n1)-cotet(3,n2))**2)
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
          fumid=max(fumid,quality)
        enddo
!     
!     determine the Jacobian in each integration point
!     notice: the penalty due to a negative Jacobian must be
!     bigger than the size of quality (volume >= 1.d-15)        
!     
        do j=1,4
          xi=gauss3d5(1,j)
          et=gauss3d5(2,j)
          ze=gauss3d5(3,j)
          call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
          if(xsj.le.0.d0) then
            fumid=fumid+1.d30*(1.d0-xsj)
            exit
          endif
        enddo
!
        indexe=ieled(2,indexe)
      enddo
!     
      return
      end
