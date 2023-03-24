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
      subroutine getnodesinitetmesh(ne,lakon,ipkon,kon,istartset,
     &     iendset,ialset,set,nset,filab,inodestet,nnodestet,
     &     nodface,ipoface,nk)
!     
!     reading the initial tet mesh which will be refined 
!     
      implicit none
!     
      character*8 lakon(*)
      character*81 set(*),elset
      character*87 filab(*)
!     
      integer ipkon(*),kon(*),istartset(*),iendset(*),ialset(*),
     &     inodestet(*),nnodestet,i,j,k,m,node,ne,nset,indexe,id,
     &     nodface(5,*),ipoface(*),nodes(3),ifree,ithree,kflag,indexold,
     &     index,ifreenew,iel,iset,j1,nk,iaux,ifacet(6,4)
!     
!     nodes belonging to the element faces
!     
      data ifacet /1,3,2,7,6,5,
     &     1,2,4,5,9,8,
     &     2,3,4,6,10,9,
     &     1,4,3,8,10,7/
!     
      kflag=1
      ithree=3
!     
!     identify the set name/s if they exist
!     
      read(filab(48)(27:87),'(a61)') elset(1:61)
!     
      do i=62,81
        elset(i:i)=' '
      enddo
!     
      call cident81(set,elset,nset,id)
      iset=nset+1
      if(id.gt.0) then
        if(elset.eq.set(id)) then
          iset=id
        endif
      endif
!     
!     determining the external element faces of the unrefined mesh 
!     the faces are catalogued by the three lowes nodes numbers
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
      do i=1,4*ne-1
        nodface(5,i)=i+1
      enddo
!     
      if(iset.le.nset) then     
        do j1=istartset(iset),iendset(iset)
          if(ialset(j1).gt.0) then
            i=ialset(j1)
!     
!     the elements belonging to the unrefined mesh have been deactivated,
!     i.e. the transformation ipkon(k)=-2-ipkon(k) was performed in
!     modelchanges.f
!     
            indexe=-2-ipkon(i)
            if(indexe.lt.0) cycle
            if((lakon(i)(1:4).eq.'C3D4').or.
     &           (lakon(i)(1:5).eq.'C3D10')) then
              do j=1,4
                do k=1,3
                  nodes(k)=kon(indexe+ifacet(k,j))
                enddo
                call isortii(nodes,iaux,ithree,kflag)
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
     &                 (nodface(2,index).eq.nodes(3))) then
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
          endif
        enddo
      else
        do i=1,ne
!     
!     the elements belonging to the unrefined mesh have been deactivated,
!     i.e. the transformation ipkon(k)=-2-ipkon(k) was performed in
!     modelchanges.f
!     
          indexe=-2-ipkon(i)
          if(indexe.lt.0) cycle
          if((lakon(i)(1:4).eq.'C3D4').or.
     &         (lakon(i)(1:5).eq.'C3D10')) then
            do j=1,4
              do k=1,3
                nodes(k)=kon(indexe+ifacet(k,j))
              enddo
              call isortii(nodes,iaux,ithree,kflag)
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
     &               (nodface(2,index).eq.nodes(3))) then
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
      endif
!     
!     storing the nodes belonging to the external faces in inodestet
!     
      do i=1,nk
        index=ipoface(i)
        do
          if(index.eq.0) exit
          iel=nodface(3,index)
          j=nodface(4,index)
          indexe=-2-ipkon(iel)
          if(lakon(iel)(4:4).eq.'4') then
            do k=1,3
              node=kon(indexe+ifacet(k,j))
              call nident(inodestet,node,nnodestet,id)
              if(id.gt.0) then
                if(inodestet(id).eq.node) cycle
              endif
              nnodestet=nnodestet+1
              do m=nnodestet,id+2,-1
                inodestet(m)=inodestet(m-1)
              enddo
              inodestet(id+1)=node
            enddo
          else
            do k=1,6
              node=kon(indexe+ifacet(k,j))
              call nident(inodestet,node,nnodestet,id)
              if(id.gt.0) then
                if(inodestet(id).eq.node) cycle
              endif
              nnodestet=nnodestet+1
              do m=nnodestet,id+2,-1
                inodestet(m)=inodestet(m-1)
              enddo
              inodestet(id+1)=node
            enddo
          endif
          index=nodface(5,index)
        enddo
      enddo
!     
      return
      end

