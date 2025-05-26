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
      subroutine extsurface(nodface,ipoface,ipkon,lakon,
     &  kon,nea,neb,ielem,iflag,ifree)
!
      implicit none
!
!     catalogues the external surface of a set of elements:
!     if iflag=0: from element nea up to element neb
!     if iflag=1: for all elements in field ielem(i),i=nea,neb
!
!     for the format of the data base: see comments further down
!
      character*8 lakon(*)
!     
      integer ipkon(*),kon(*),nodface(5,*),ipoface(*),
     &     ithree,ifour,ifaceq(8,6),
     &     ifacet(6,4),ifacew2(8,5),ifree,ifreenew,index,indexold,
     &     i,j,k,nodes(4),indexe,ii,
     &     ifacew1(4,5),iflag,ielem(*),nea,neb
!     
!     nodes belonging to the element faces
!     
      data ifaceq /4,3,2,1,11,10,9,12,
     &     5,6,7,8,13,14,15,16,
     &     1,2,6,5,9,18,13,17,
     &     2,3,7,6,10,19,14,18,
     &     3,4,8,7,11,20,15,19,
     &     4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,
     &     1,2,4,5,9,8,
     &     2,3,4,6,10,9,
     &     1,4,3,8,10,7/
      data ifacew1 /1,3,2,0,
     &     4,5,6,0,
     &     1,2,5,4,
     &     2,3,6,5,
     &     3,1,4,6/
      data ifacew2 /1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     3,1,4,6,9,13,12,15/
!     
      ithree=3
      ifour=4
!     
!     determining the external element faces of the mesh; 
!     the faces are catalogued by the three lowest nodes numbers
!     in ascending order. ipoface(i) points to a face for which
!     node i is the lowest node and nodface(1,ipoface(i)) and
!     nodface(2,ipoface(i)) are the next lower ones. 
!     nodface(3,ipoface(i)) contains the element number,
!     nodface(4,ipoface(i)) the face number and nodface(5,ipoface(i))
!     is a pointer to the next surface for which node i is the
!     lowest node; if there are no more such surfaces the pointer
!     has the value zero
!     An external element face is one which belongs to one element
!     only
!     
      ifree=1
      do i=1,6*(neb-nea+1)-1
        nodface(5,i)=i+1
      enddo
!     
      do ii=1,nea,neb
!     
        if(iflag.eq.1) then
          i=ielem(ii)
        else
          i=ii
        endif
!     
!     only active elements
!     
        if(ipkon(i).lt.0) cycle
!     
!     only 3D elements
!     
        if(lakon(i)(1:3).ne.'C3D') cycle
!     
        indexe=ipkon(i)
!     
        if((lakon(i)(4:4).eq.'2').or.(lakon(i)(4:4).eq.'8')) then
          do j=1,6
            do k=1,4
              nodes(k)=kon(indexe+ifaceq(k,j))
            enddo
            call insertsorti(nodes,ifour)
            indexold=0
            index=ipoface(nodes(1))
            do
!     
!     adding a surface which has not been 
!     catalogued so far
!     
              if(index.eq.0) then
                ifreenew=nodface(5,ifree)
                nodface(1,ifree)=nodes(2)
                nodface(2,ifree)=nodes(3)
                nodface(3,ifree)=i
                nodface(4,ifree)=j
                nodface(5,ifree)=ipoface(nodes(1))
                ipoface(nodes(1))=ifree
                ifree=ifreenew
                exit
              endif
!     
!     removing a surface which has already
!     been catalogued
!     
              if((nodface(1,index).eq.nodes(2)).and.
     &             (nodface(2,index).eq.nodes(3))) then
                if(indexold.eq.0) then
                  ipoface(nodes(1))=nodface(5,index)
                else
                  nodface(5,indexold)=nodface(5,index)
                endif
                nodface(5,index)=ifree
                ifree=index
                exit
              endif
              indexold=index
              index=nodface(5,index)
            enddo
          enddo
        elseif((lakon(i)(4:4).eq.'4').or.(lakon(i)(4:5).eq.'10')) then
          do j=1,4
            do k=1,3
              nodes(k)=kon(indexe+ifacet(k,j))
            enddo
            call insertsorti(nodes,ithree)
            indexold=0
            index=ipoface(nodes(1))
            do
!     
!     adding a surface which has not been 
!     catalogued so far
!     
              if(index.eq.0) then
                ifreenew=nodface(5,ifree)
                nodface(1,ifree)=nodes(2)
                nodface(2,ifree)=nodes(3)
                nodface(3,ifree)=i
                nodface(4,ifree)=j
                nodface(5,ifree)=ipoface(nodes(1))
                ipoface(nodes(1))=ifree
                ifree=ifreenew
                exit
              endif
!     
!     removing a surface which has already
!     been catalogued
!     
              if((nodface(1,index).eq.nodes(2)).and.
     &             (nodface(2,index).eq.nodes(3))) then
                if(indexold.eq.0) then
                  ipoface(nodes(1))=nodface(5,index)
                else
                  nodface(5,indexold)=nodface(5,index)
                endif
                nodface(5,index)=ifree
                ifree=index
                exit
              endif
              indexold=index
              index=nodface(5,index)
            enddo
          enddo
        else
          do j=1,5
            if(j.le.2) then
              do k=1,3
                nodes(k)=kon(indexe+ifacew2(k,j))
              enddo
              call insertsorti(nodes,ithree)
            else
              do k=1,4
                nodes(k)=kon(indexe+ifacew2(k,j))
              enddo
              call insertsorti(nodes,ifour)
            endif
            indexold=0
            index=ipoface(nodes(1))
            do
!     
!     adding a surface which has not been 
!     catalogued so far
!     
              if(index.eq.0) then
                ifreenew=nodface(5,ifree)
                nodface(1,ifree)=nodes(2)
                nodface(2,ifree)=nodes(3)
                nodface(3,ifree)=i
                nodface(4,ifree)=j
                nodface(5,ifree)=ipoface(nodes(1))
                ipoface(nodes(1))=ifree
                ifree=ifreenew
                exit
              endif
!     
!     removing a surface which has already
!     been catalogued
!     
              if((nodface(1,index).eq.nodes(2)).and.
     &             (nodface(2,index).eq.nodes(3))) then
                if(indexold.eq.0) then
                  ipoface(nodes(1))=nodface(5,index)
                else
                  nodface(5,indexold)=nodface(5,index)
                endif
                nodface(5,index)=ifree
                ifree=index
                exit
              endif
              indexold=index
              index=nodface(5,index)
            enddo
          enddo
        endif
      enddo
!     
      return
      end
