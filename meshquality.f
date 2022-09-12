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
      subroutine meshquality(netet_,kontet,cotet,quality,ielem)
!
!     calculate the element quality. The measure used is proportional
!     to the ratio of the longest edge divided by the radius of the
!     inscribed sphere. The proportionality constant is such that
!     the quality is 1 for an equilateral tetrahedron. For all other
!     elements it exceeds 1. The bigger this number, the worse the
!     quality
!
!     if ielem>0 the quality for this element is calculated, if
!     ielem=0 the quality for all elements in the mesh is determined
!
      implicit none
!
      integer netet_,kontet(4,*),ielem,i,j,nodes(4),n1,n2,ine(2,6)
!
      real*8 cotet(3,*),quality(*),volume,surface(4),totsurface,alpha,
     &     edgelength(6),hmax,radius
!
      data ine /1,2,2,3,1,3,1,4,2,4,3,4/
!
!     alpha is the proporionality factor
!            
      alpha=dsqrt(6.d0)/12.d0
!     
      if(ielem.eq.0) then
!
!     calculating the quality for all elements
!
        do i=1,netet_
          if(kontet(1,i).ne.0) then
            do j=1,4
              nodes(j)=kontet(j,i)
            enddo
!
!     calculating the volume of the element
!
            call calcvol(nodes(1),nodes(2),nodes(3),nodes(4),cotet,
     &           volume)
            if(volume.le.0.d0) volume=1.d-30
!
!     calculating area of each face in the element
!
            call calcsurf(nodes(1),nodes(2),nodes(3),cotet,
     &           surface(1))
            call calcsurf(nodes(2),nodes(3),nodes(4),cotet,
     &           surface(2))
            call calcsurf(nodes(3),nodes(4),nodes(1),cotet,
     &           surface(3))
            call calcsurf(nodes(4),nodes(1),nodes(2),cotet,
     &           surface(4))
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
     &             (cotet(2,n1)-cotet(2,n2))**2+
     &             (cotet(3,n1)-cotet(3,n2))**2)
            enddo
!
!     maximum edge length
!
            hmax=maxval(edgelength)
!
!     quality
!
            quality(i)=alpha*hmax/radius
!
!     are the next lines really needed?
!
c            if(quality(i).lt.1.d0) then
c              quality(i)=1.d0/quality(i)
c           endif
          endif
        enddo
      else
!
!     calculating the quality for just one element
!
        i=ielem
!        
        if(kontet(1,i).ne.0) then
          do j=1,4
            nodes(j)=kontet(j,i)
          enddo
!     
!     calculating the volume of the element
!     
          call calcvol(nodes(1),nodes(2),nodes(3),nodes(4),cotet,
     &         volume)
          if(volume.le.0.d0) volume=1.d-30
!     
!     calculating area of each face in the element
!     
          call calcsurf(nodes(1),nodes(2),nodes(3),cotet,
     &         surface(1))
          call calcsurf(nodes(2),nodes(3),nodes(4),cotet,
     &         surface(2))
          call calcsurf(nodes(3),nodes(4),nodes(1),cotet,
     &         surface(3))
          call calcsurf(nodes(4),nodes(1),nodes(2),cotet,
     &         surface(4))
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
     &           (cotet(2,n1)-cotet(2,n2))**2+
     &           (cotet(3,n1)-cotet(3,n2))**2)
          enddo
!     
!     maximum edge length
!     
          hmax=maxval(edgelength)
!     
!     quality
!     
          quality(i)=alpha*hmax/radius
!     
!     are the next lines really needed?
!     
c          if(quality(i).lt.1.d0) then
c            quality(i)=1.d0/quality(i)
c          endif
        endif
      endif
!            
      return
      end
