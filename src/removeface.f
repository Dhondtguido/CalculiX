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
      subroutine removeface(nodface,ipoface,ipkon,lakon,
     &  kon,i,j,ifree)
!
      implicit none
!
!     removes a face (local face j of element i)
!     from the database (ipoface,nodface). For the
!     structure of the database cf. extsurface.f
!
      character*8 lakon(*)
!     
      integer ipkon(*),kon(*),nodface(5,*),ipoface(*),ithree,ifour,
     &     ifaceq(8,6),ifacet(6,4),ifacew2(8,5),ifree,index,indexold,
     &     i,j,k,nodes(4),indexe,ifacew1(4,5)
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
!     check whether active element
!     
      if(ipkon(i).lt.0) return
!     
!     check whether the element is 3-dimensional
!     
      if(lakon(i)(1:3).ne.'C3D') return
!     
      indexe=ipkon(i)
!     
      if((lakon(i)(4:4).eq.'2').or.(lakon(i)(4:4).eq.'8')) then
        do k=1,4
          nodes(k)=kon(indexe+ifaceq(k,j))
        enddo
        call insertsorti(nodes,ifour)
        indexold=0
        index=ipoface(nodes(1))
        do
          if(index.eq.0) exit
!     
!     removing a surface which has already
!     been catalogued
!     
          if((nodface(1,index).eq.nodes(2)).and.
     &         (nodface(2,index).eq.nodes(3))) then
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
      elseif((lakon(i)(4:4).eq.'4').or.(lakon(i)(4:5).eq.'10')) then
        do k=1,3
          nodes(k)=kon(indexe+ifacet(k,j))
        enddo
        call insertsorti(nodes,ithree)
        indexold=0
        index=ipoface(nodes(1))
        do
          if(index.eq.0) exit
!     
!     removing a surface which has already
!     been catalogued
!     
          if((nodface(1,index).eq.nodes(2)).and.
     &         (nodface(2,index).eq.nodes(3))) then
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
      else
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
          if(index.eq.0) exit
!     
!     removing a surface which has already
!     been catalogued
!     
          if((nodface(1,index).eq.nodes(2)).and.
     &         (nodface(2,index).eq.nodes(3))) then
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
      endif
!     
      return
      end
