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
      subroutine calcnormal(nelem,jface,lakon,co,xn,indexe,kon)
!
!     determine the local normal on face "jface" of element "nelem".
!
      implicit none
!
      logical line,quad
!
      character*8 lakon(*)
!
      integer i,j,ifacequad(3,4),ifacetria(3,3),nelem,jface,nopes,
     &  nface,nodef(8),ifaceq(8,6),ifacet(6,4),ifacew1(4,5),
     &  ifacew2(8,5),nope,kon(*),indexe,iflag
!
      real*8 xn(3),xt(3),xd(3),dd,co(3,*),xl(3,8),xi,et,
     &  xi3(3),et3(3),xi4(4),et4(4),xs(3,7),shp(7,8)
!
!     nodes per face for hex elements
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
!
!     nodes per face for tet elements
!
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
!
!     nodes per face for linear wedge elements
!
      data ifacew1 /1,3,2,0,
     &             4,5,6,0,
     &             1,2,5,4,
     &             2,3,6,5,
     &             3,1,4,6/
!
!     nodes per face for quadratic wedge elements
!
      data ifacew2 /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             3,1,4,6,9,13,12,15/
!
!     nodes per face for quad elements
!
      data ifacequad /1,2,5,
     &                2,3,6,
     &                3,4,7,
     &                4,1,8/
!
!     nodes per face for tria elements
!
      data ifacetria /1,2,4,
     &                2,3,5,
     &                3,1,6/
!
!     flag for shape functions
!
      data iflag /2/
!
      data xi3 /0.5d0,0.5d0,0.d0/
      data et3 /0.d0,0.5d0,0.5d0/
      data xi4 /0.d0,1.d0,0.d0,-1.d0/
      data et4 /-1.d0,0.d0,1.d0,0.d0/
!
!     line=.true. means that the surface is reduced to a line,
!     i.e. it is a face of a plane stress, plane strain,
!     axisymmetric or shell element
!     initialization:
!
      line=.false.
!
!     nodes: #nodes in the face
!     the nodes are stored in nodef(*)
!
      if(lakon(nelem)(4:4).eq.'2') then
         nopes=8
         nface=6
      elseif(lakon(nelem)(3:4).eq.'D8') then
         nopes=4
         nface=6
      elseif(lakon(nelem)(4:5).eq.'10') then
         nopes=6
         nface=4
      elseif(lakon(nelem)(3:4).eq.'D4') then
         nopes=3
         nface=4
      elseif(lakon(nelem)(4:5).eq.'15') then
         if(jface.le.2) then
            nopes=6
         else
            nopes=8
         endif
         nface=5
         nope=15
      elseif(lakon(nelem)(3:4).eq.'D6') then
         if(jface.le.2) then
            nopes=3
         else
            nopes=4
         endif
         nface=5
         nope=6
      elseif((lakon(nelem)(2:2).eq.'8').or.
     &        (lakon(nelem)(4:4).eq.'8')) then
!     
!     8-node 2-D elements
!     
         nopes=3
         nface=4
         quad=.true.
c         if(lakon(nelem)(4:4).eq.'8') then
            line=.true.
            jface=jface-2
c         endif
      elseif((lakon(nelem)(2:2).eq.'6').or.
     &        (lakon(nelem)(4:4).eq.'6')) then
!     
!     6-node 2-D elements
!     
         nopes=3
         nface=3
c         if(lakon(nelem)(4:4).eq.'6') then
            line=.true.
            jface=jface-2
c         endif
      elseif((lakon(nelem)(2:2).eq.'4').or.
     &        (lakon(nelem)(4:4).eq.'4')) then
!     
!     4-node 2-D elements
!     
         nopes=2
         nface=4
         quad=.true.
c         if(lakon(nelem)(4:4).eq.'4') then
            line=.true.
            jface=jface-2
c         endif
      elseif((lakon(nelem)(2:2).eq.'3').or.
     &        (lakon(nelem)(4:4).eq.'3')) then
!     
!     3-node 2-D elements
!     
         nopes=2
         nface=3
c         if(lakon(nelem)(4:4).eq.'3') then
            line=.true.
            jface=jface-2
c         endif
      endif
!     
!     determining the nodes of the face
!     
      if(nface.eq.3) then
         do i=1,nopes
            nodef(i)=kon(indexe+ifacetria(i,jface))
         enddo
      elseif(nface.eq.4) then
         if(quad) then
            do i=1,nopes
               nodef(i)=kon(indexe+ifacequad(i,jface))
            enddo
         else
            do i=1,nopes
               nodef(i)=kon(indexe+ifacet(i,jface))
            enddo
         endif
      elseif(nface.eq.5) then
         if(nope.eq.6) then
            do i=1,nopes
               nodef(i)=kon(indexe+ifacew1(i,jface))
            enddo
         elseif(nope.eq.15) then
            do i=1,nopes
               nodef(i)=kon(indexe+ifacew2(i,jface))
            enddo
         endif
      elseif(nface.eq.6) then
         do i=1,nopes
            nodef(i)=kon(indexe+ifaceq(i,jface))
         enddo
      endif
!
!     storing the nodes in the face
!
      do i=1,nopes
         do j=1,3
            xl(j,i)=co(j,nodef(i))
         enddo
      enddo
!
      iflag=2
      if(nopes.eq.2) then
!
!        side of a 3-node triangular or a 4-node
!        quadrilateral element
!
         do j=1,3
            xt(j)=xl(j,2)-xl(j,1)
         enddo
         if(nface.eq.3) then
            nope=3
            do i=1,nope
               do j=1,3
                  xl(j,i)=co(j,kon(indexe+i))
               enddo
            enddo
            xi=xi3(jface)
            et=et3(jface)
            call shape3tri(xi,et,xl,xd,xs,shp,iflag)
         elseif(nface.eq.4) then
            nope=4
            do i=1,nope
               do j=1,3
                  xl(j,i)=co(j,kon(indexe+i))
               enddo
            enddo
            xi=xi4(jface)
            et=et4(jface)
            call shape4q(xi,et,xl,xd,xs,shp,iflag)
         endif
         xn(1)=xt(2)*xd(3)-xt(3)*xd(2)
         xn(2)=xt(3)*xd(1)-xt(1)*xd(3)
         xn(3)=xt(1)*xd(2)-xt(2)*xd(1)
      elseif(nopes.eq.3) then
!
!        side of a 6-node triangular element (2D) or a
!        8-node quadrilateral element (2D) or a
!        4-node tetrahedral element (3D)
!         
         if(line) then
            xi=0.d0
            call shape3l(xi,xl,xt,xs,shp,iflag)
            if(nface.eq.3) then
               nope=6
               do i=1,nope
                  do j=1,3
                     xl(j,i)=co(j,kon(indexe+i))
                  enddo
               enddo
               xi=xi3(jface)
               et=et3(jface)
               call shape6tri(xi,et,xl,xd,xs,shp,iflag)
            elseif(nface.eq.4) then
               nope=8
               do i=1,nope
                  do j=1,3
                     xl(j,i)=co(j,kon(indexe+i))
                  enddo
               enddo
               xi=xi4(jface)
               et=et4(jface)
               call shape8q(xi,et,xl,xd,xs,shp,iflag)
            endif
            xn(1)=xt(2)*xd(3)-xt(3)*xd(2)
            xn(2)=xt(3)*xd(1)-xt(1)*xd(3)
            xn(3)=xt(1)*xd(2)-xt(2)*xd(1)
         else
            xi=0.d0
            et=0.d0
            call shape3tri(xi,et,xl,xn,xs,shp,iflag)
         endif
      elseif(nopes.eq.4) then
!
!        side of a 8-node hex element
!
         xi=0.d0
         et=0.d0
         call shape4q(xi,et,xl,xn,xs,shp,iflag)
      elseif(nopes.eq.6) then
!
!        side of a 10-node tet element
!
         xi=1.d0/3.d0
         et=1.d0/3.d0
         call shape6tri(xi,et,xl,xn,xs,shp,iflag)
      elseif(nopes.eq.8) then
         xi=0.d0
         et=0.d0
         call shape8q(xi,et,xl,xn,xs,shp,iflag)
      endif
!
!     normalizing
!
      dd=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
      do i=1,3
         xn(i)=xn(i)/dd
      enddo
!     
      return
      end

