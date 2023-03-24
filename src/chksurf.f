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
      subroutine chksurf(lakon,kon,ipkon,neigh,ipneigh,co,itypflag,node,
     & icont,iscount,angmax)
!
!     icont=1: element surfaces adjacent to a surface node have normal
!                vectors which have an angle of less than 10 degree
!               -> free surface assumed
!     icont=0: -> edge assumed
!
!     also counts the free surfaces adjacent to a node
!
!     author: Sascha Merz
!
      implicit none
!
      integer kon(*),ipkon(*),ielem,i,j,k,indexe,
     & neigh(2,*),ipneigh(*),index,m,nvertex,itypflag,isurf,node,index1,
     & ielem1,ncount,ntos8h(3,8),ntos4tet(3,4),iston8h(4,6),
     & isnode, isidx,ifreesur(3),iston20h(8,6),iston10tet(6,4),
     & iscount,lnod,icont,iston4tet(3,4)
!
      real*8 co(3,*),angle,shpder8q(2,4,8),shpder6tri(2,3,6),
     & vectors(3,3),vlen(2),lastvec(3),angtmp,xl(3,8),
     & shpder4q(2,4,4),shpder3tri(2,3,3),angmax
!     
!     ntosX(j,k) returns the three surface id's j for the corner node k
!     for the element surfaces adjacent to the node
!
      data ntos8h /1,3,6,1,3,4,1,4,5,1,5,6,2,3,6,2,3,4,2,4,5,2,5,6/
!
      data ntos4tet /1,2,4,1,2,3,1,3,4,2,3,4/
!
!     istonX(j,k) returns the nodes j of the element surface k
!
      data iston8h /1,2,3,4,5,8,7,6,1,5,6,2,2,6,7,3,3,7,8,4,4,8,5,1/
!
      data iston20h /1,2,3,4,9,10,11,12,5,8,7,6,16,15,14,13,
     &              1,5,6,2,17,13,18,9,2,6,7,3,18,14,19,10,
     &              3,7,8,4,19,15,20,11,4,8,5,1,20,16,17,12/
!
      data iston4tet /1,2,3,1,4,2,2,4,3,3,4,1/
!
      data iston10tet /1,2,3,5,6 ,7,1,4,2,8 ,9,5,
     &                 2,4,3,9,10,6,3,4,1,10,8,7/
!
!     shpder8q contains the first derivative of the shape functions 
!     of a 8 node quadrilateral element with shpder8q(i,j,k) where i
!     can be 1 for the xi-derivative or 2 for the eta-derivative j
!     can be 1-4 for the location in the corner nodes to be evaluated
!     for the element nodes k
!
      data shpder8q /
     & -1.5d0,-1.5d0, 0.5d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.5d0,
     & -0.5d0, 0.0d0, 1.5d0,-1.5d0, 0.0d0, 0.5d0, 0.0d0, 0.0d0,
     &  0.0d0, 0.0d0, 0.0d0,-0.5d0, 1.5d0, 1.5d0,-0.5d0, 0.0d0,
     &  0.0d0,-0.5d0, 0.0d0, 0.0d0, 0.5d0, 0.0d0,-1.5d0, 1.5d0,
     &  2.0d0, 0.0d0,-2.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,
     &  0.0d0, 0.0d0, 0.0d0, 2.0d0, 0.0d0,-2.0d0, 0.0d0, 0.0d0,
     &  0.0d0, 0.0d0, 0.0d0, 0.0d0,-2.0d0, 0.0d0, 2.0d0, 0.0d0,
     &  0.0d0, 2.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,-2.0d0/
!
!     same as above for a 4 node linear quadrilateral element
!
      data shpder4q /
     & -0.5d0,-0.5d0,-0.5d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,-0.5d0,
     &  0.5d0, 0.0d0, 0.5d0,-0.5d0, 0.0d0,-0.5d0, 0.0d0, 0.0d0,
     &  0.0d0, 0.0d0, 0.0d0, 0.5d0, 0.5d0, 0.5d0, 0.5d0, 0.0d0,
     &  0.0d0, 0.5d0, 0.0d0, 0.0d0,-0.5d0, 0.0d0,-0.5d0, 0.5d0/
!
!     same as above for a 6 node quadratic triangular element
!
      data shpder6tri /
     & -3.0d0,-3.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0,
     & -1.0d0, 0.0d0, 3.0d0, 0.0d0,-1.0d0, 0.0d0,
     &  0.0d0,-1.0d0, 0.0d0,-1.0d0, 0.0d0, 3.0d0,
     &  4.0d0, 0.0d0,-4.0d0,-4.0d0, 0.0d0, 0.0d0,
     &  0.0d0, 0.0d0, 0.0d0, 4.0d0, 4.0d0, 0.0d0,
     &  0.0d0, 4.0d0, 0.0d0, 0.0d0,-4.0d0,-4.0d0/
!
!     same as above for a 3 node linear triangular element
!
      data shpder3tri /
     & -1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,
     &  1.0d0, 0.0d0, 1.0d0, 0.0d0, 1.0d0, 0.0d0,
     &  0.0d0, 1.0d0, 0.0d0, 1.0d0, 0.0d0, 1.0d0/
!
      character*8 lakon(*)
!
      index=ipneigh(node)
      icont=1
      iscount=0
      angmax=0.d0
!
      do
         if(index.eq.0) exit
         ielem=neigh(1,index)
!
         if(lakon(ielem)(1:5).eq.'C3D20'.and.itypflag.eq.1) then
            nvertex=8
         elseif(lakon(ielem)(1:5).eq.'C3D10'.and.itypflag.eq.2) then
            nvertex=4
         elseif(lakon(ielem)(1:4).eq.'C3D8'.and.itypflag.eq.3) then
            nvertex=8
         elseif(lakon(ielem)(1:4).eq.'C3D4'.and.itypflag.eq.4) then
            nvertex=4
         else
            index=neigh(2,index)
            cycle
         endif
!
!        find the index of the node in the element
!
         indexe=ipkon(ielem)
         do m=1,nvertex
            if(kon(indexe+m).eq.node) exit
         enddo
!
!        the local node number is m
!
!        now every surface has to be checked
!
         do j=1,3
            ifreesur(j)=0
         enddo
!
         do isurf=1,3
!
!           finding the global node numbers of the
!           nodes of the surface
!            
            if(nvertex.eq.4) then
!
!              isidx: index of the surface neighbouring the node
!               
               isidx=ntos4tet(isurf,m)
            elseif(nvertex.eq.8) then
               isidx=ntos8h(isurf,m)
            endif
!
!           find out, if there is any element neighbouring 'node',
!           which has also those nodes (-> surface is within volume)
!            
            index1=ipneigh(node)
            do
               if(index1.eq.0) exit
               ielem1=neigh(1,index1)
               if(
     &            .not.(
     &                  lakon(ielem1)(1:5).eq.'C3D20'.and.itypflag.eq.1
     &                  .or.
     &                  lakon(ielem1)(1:5).eq.'C3D10'.and.itypflag.eq.2
     &                  .or.
     &                  lakon(ielem1)(1:4).eq.'C3D8'.and.itypflag.eq.3
     &                  .or.
     &                  lakon(ielem1)(1:4).eq.'C3D4'.and.itypflag.eq.4
     &                 )
     &            .or.ielem.eq.ielem1
     &           ) then
                  index1=neigh(2,index1)
                  cycle
               endif
!
!              check every corner node in the element
!               
               ncount=0
               do k=1,3
                  if(nvertex.eq.4) then
                     isnode=kon(indexe+iston4tet(k,isidx))
                  elseif(nvertex.eq.8) then
                     isnode=kon(indexe+iston8h(k,isidx))
                  endif
                  do j=1,nvertex
                     if(kon(ipkon(ielem1)+j).eq.isnode)
     &                    ncount=ncount+1
                  enddo
               enddo
!
               if(ncount.eq.3) then
!
!                 surface isurf is not a free surface
!
                  ifreesur(isurf)=1
               endif
!
               index1=neigh(2,index1)
            enddo
         enddo
!
         do isurf=1,3
            if(ifreesur(isurf).eq.0) then
               iscount=iscount+1
               do i=1,3
                  do j=1,3
                     vectors(j,i)=0.d0
                  enddo
               enddo
!
!              free surface: find out local node number
!              of the 'surface element' neighbouring the
!              node to be evaluated
!              
               if(nvertex.eq.8) then
                  isidx=ntos8h(isurf,m)
                  do j=1,4
                     if(      (isidx.eq.1.and.m.eq.1)
     &                    .or.(isidx.eq.2.and.m.eq.5)
     &                    .or.(isidx.eq.3.and.m.eq.1)
     &                    .or.(isidx.eq.4.and.m.eq.2)
     &                    .or.(isidx.eq.5.and.m.eq.3)
     &                    .or.(isidx.eq.6.and.m.eq.4)) then
                        lnod=1
                     elseif(  (isidx.eq.1.and.m.eq.2)
     &                    .or.(isidx.eq.2.and.m.eq.8)
     &                    .or.(isidx.eq.3.and.m.eq.5)
     &                    .or.(isidx.eq.4.and.m.eq.6)
     &                    .or.(isidx.eq.5.and.m.eq.7)
     &                    .or.(isidx.eq.6.and.m.eq.8)) then
                        lnod=2
                     elseif(  (isidx.eq.1.and.m.eq.3)
     &                    .or.(isidx.eq.2.and.m.eq.7)
     &                    .or.(isidx.eq.3.and.m.eq.6)
     &                    .or.(isidx.eq.4.and.m.eq.7)
     &                    .or.(isidx.eq.5.and.m.eq.8)
     &                    .or.(isidx.eq.6.and.m.eq.5)) then
                        lnod=3
                     elseif(  (isidx.eq.1.and.m.eq.4)
     &                    .or.(isidx.eq.2.and.m.eq.6)
     &                    .or.(isidx.eq.3.and.m.eq.2)
     &                    .or.(isidx.eq.4.and.m.eq.3)
     &                    .or.(isidx.eq.5.and.m.eq.4)
     &                    .or.(isidx.eq.6.and.m.eq.1)) then
                        lnod=4
                     endif
                  enddo
!     
c                  do k=1,8
c                     write(*,*) 'node',node,' nodes',
c     &                    kon(indexe+iston20h(k,isidx)),'lnod',lnod
c                  enddo
!
                  if(itypflag.eq.1) then
!
!                    get coordinates of the 2d-element nodes
!
                     do k=1,8
                        do j=1,3
                           xl(j,k)=co(j,kon(indexe+iston20h(k,isidx)))
                        enddo
                     enddo
!
!                    vectors(j,i) (i=1,2) is the j-derivative for the
!                    coordinates i.
!
                     do k=1,8
                        do i=1,3
                           do j=1,2
                              vectors(j,i)=vectors(j,i)+
     &                             xl(i,k)*shpder8q(j,lnod,k)
                           enddo
                        enddo
                     enddo
!
                  elseif(itypflag.eq.3) then
                     do k=1,4
                        do j=1,3
                           xl(j,k)=co(j,kon(indexe+iston8h(k,isidx)))
                        enddo
                     enddo
                     do k=1,4
                        do i=1,3
                           do j=1,2
                              vectors(j,i)=vectors(j,i)+
     &                             xl(i,k)*shpder4q(j,lnod,k)
                           enddo
                        enddo
                     enddo
                  endif
!
               elseif(nvertex.eq.4) then
                  isidx=ntos4tet(isurf,m)
                  do j=1,3
                     if(      (isidx.eq.1.and.m.eq.1)
     &                    .or.(isidx.eq.2.and.m.eq.1)
     &                    .or.(isidx.eq.3.and.m.eq.2)
     &                    .or.(isidx.eq.4.and.m.eq.3)) then
                        lnod=1
                     elseif(  (isidx.eq.1.and.m.eq.2)
     &                    .or.(isidx.eq.2.and.m.eq.4)
     &                    .or.(isidx.eq.3.and.m.eq.4)
     &                    .or.(isidx.eq.4.and.m.eq.4)) then
                        lnod=2
                     elseif(  (isidx.eq.1.and.m.eq.3)
     &                    .or.(isidx.eq.2.and.m.eq.2)
     &                    .or.(isidx.eq.3.and.m.eq.3)
     &                    .or.(isidx.eq.4.and.m.eq.1)) then
                        lnod=3
                     endif
                  enddo
!
                  if(itypflag.eq.2) then
                     do k=1,6
                        do j=1,3
                           xl(j,k)=co(j,kon(indexe+iston10tet(k,isidx)))
                        enddo
                     enddo
                     do k=1,6
                        do i=1,3
                           do j=1,2
                              vectors(j,i)=vectors(j,i)+
     &                             xl(i,k)*shpder6tri(j,lnod,k)
                           enddo
                        enddo
                     enddo
                  elseif(itypflag.eq.4) then
                     do k=1,3
                        do j=1,3
                           xl(j,k)=co(j,kon(indexe+iston4tet(k,isidx)))
                        enddo
                     enddo
                     do k=1,3
                        do i=1,3
                           do j=1,2
                              vectors(j,i)=vectors(j,i)+
     &                             xl(i,k)*shpder3tri(j,lnod,k)
                           enddo
                        enddo
                     enddo
                  endif
!
               endif
!
!              vectors(3,i) is the normal vector of the surface in the
!              evaluated node 'node'
!
               vectors(3,1)=vectors(1,2)*vectors(2,3)
     &              -vectors(1,3)*vectors(2,2)
               vectors(3,2)=vectors(1,3)*vectors(2,1)
     &              -vectors(1,1)*vectors(2,3)
               vectors(3,3)=vectors(1,1)*vectors(2,2)
     &              -vectors(1,2)*vectors(2,1)
               vlen(2)=dsqrt(vectors(3,1)*vectors(3,1)
     &              +vectors(3,2)*vectors(3,2)
     &              +vectors(3,3)*vectors(3,3))
!     
               if(iscount.gt.1) then
                  angtmp=dabs((vectors(3,1)*lastvec(1)
     &                 +vectors(3,2)*lastvec(2)
     &                 +vectors(3,3)*lastvec(3))
     &                 /(vlen(1)*vlen(2)))
                  if(angtmp.lt.1.d0) then
                     angle=57.29577951d0*dacos(angtmp)
                     if(angle.gt.angmax) angmax=angle
                  else
                     angle=0.d0
                  endif
!
!                 if the angle between the normal vectors of two
!                 surfaces is greater than 10 degree, than it is
!                 assumed that an edge (dicontinuity) is present
!
                  if(angle.ge.10.d0) then
                     icont=0
                  endif
               endif
!
               do j=1,3
                  lastvec(j)=vectors(3,j)
               enddo
               vlen(1)=vlen(2)
!                  
            endif
         enddo
         index=neigh(2,index)
      enddo
      return
      end
