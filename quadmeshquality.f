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
      subroutine quadmeshquality(netet_,cotet,kontet,iedtet,
     &     iedgmid,qualityjac,ielem)
!     
!     calculate the element quality. At first the jacobian in all 4
!     integration points of the 10-node tet are calculated. If the minimum
!     value is negative the quality is set to a very big number. If not,
!     the quality is the difference of the extremes divided by the sum
!     of the extremes. This quality measure is good if the vertex nodes
!     of the tetrahedron are fixed.
!     
!     if ielem>0 the quality for this element is calculated, if
!     ielem=0 the quality for all elements in the mesh is determined
!     
      implicit none
!     
      integer netet_,iedtet(6,*),kontet(4,*),k,iedgmid(*),i,j,iflag,
     &     ielem
!     
      real*8 cotet(3,*),qualityjac(*),xjac(4),xsj,xi,et,ze,xl(3,10),
     &     shp(4,10),maxjac,minjac
!     
      include "gauss.f"
!     
      iflag=2
!     
      if(ielem.eq.0) then
!     
!     calculating the quality for all elements
!     
        do i=1,netet_
          if(kontet(1,i).eq.0) cycle
!     
!     coordinates of the vertex nodes
!     
          do j=1,4
            do k=1,3
              xl(k,j)=cotet(k,kontet(j,i))
            enddo
          enddo
!     
!     coordinates of the midnodes
!     
          do j=1,6
            do k=1,3
              xl(k,j+4)=cotet(k,iedgmid(iedtet(j,i)))
            enddo
          enddo
!     
!     determine the jacobian at the four integration points
!     
          do j=1,4
            xi=gauss3d5(1,j)
            et=gauss3d5(2,j)
            ze=gauss3d5(3,j)
            call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
            xjac(j)=xsj
          enddo
!     
!     find the maximum and minimum jacobian
!     
          maxjac=max(xjac(1),xjac(2),xjac(3),xjac(4))
          minjac=min(xjac(1),xjac(2),xjac(3),xjac(4))
!     
!     calculate the quality
!     
          if(minjac.le.0.d0) then
            qualityjac(i)=1.d10-1.d10*minjac
          else
            qualityjac(i)=(maxjac-minjac)/(maxjac+minjac)
          endif
        enddo
      else
!     
!     calculating the quality for just one element
!     
        i=ielem
!     
!     coordinates of the vertex nodes
!     
        do j=1,4
          do k=1,3
            xl(k,j)=cotet(k,kontet(j,i))
          enddo
        enddo
!     
!     coordinates of the midnodes
!     
        do j=1,6
          do k=1,3
            xl(k,j+4)=cotet(k,iedgmid(iedtet(j,i)))
          enddo
        enddo
!     
!     determine the jacobian at the four integration points
!     
        do j=1,4
          xi=gauss3d5(1,j)
          et=gauss3d5(2,j)
          ze=gauss3d5(3,j)
          call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
          xjac(j)=xsj
        enddo
!     
!     find the maximum and minimum jacobian
!     
        maxjac=max(xjac(1),xjac(2),xjac(3),xjac(4))
        minjac=min(xjac(1),xjac(2),xjac(3),xjac(4))
!     
!     calculate the quality
!     
        if(minjac.le.0.d0) then
          qualityjac(i)=1.d10-1.d10*minjac
        else
          qualityjac(i)=(maxjac-minjac)/(maxjac+minjac)
        endif
      endif
!     
      return
      end
