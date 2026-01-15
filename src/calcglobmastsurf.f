!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine calcglobmastsurf(ne,ipkon,kon,lakon,nk,
     &     set,istartset,iendset,ialset,nset,nalset,
     &     imastset,nmastface,interfaceload)
!     
!     preliminary calculations for I-type loading applications:
!     determining the external faces of the mesh and storing
!     them in a master surface with the name EXTERNALFACES1235711.
!     
      implicit none
!     
      character*8 lakon(*)
      character*81 set(*),noelset
!     
      integer nodes(4),ne,ipkon(*),kon(*),
     &     indexe,ifaceq(8,6),ifacet(7,4),index,ifacew(8,5),ithree,
     &     ifour,iaux,kflag,nalset,id,nk,i,j,k,m,ifree,indexold,
     &     ifreenew,istartset(*),iendset(*),ialset(*),nset,imastset,
     &     nmastface,interfaceload
!
      integer,dimension(:),allocatable::ipoface
      integer,dimension(:,:),allocatable::nodface
!     
!     determine the external faces of the volumetric elements and
!     store them in a global master surface
!     
      data ifaceq /4,3,2,1,11,10,9,12,
     &     5,6,7,8,13,14,15,16,
     &     1,2,6,5,9,18,13,17,
     &     2,3,7,6,10,19,14,18,
     &     3,4,8,7,11,20,15,19,
     &     4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,11,
     &     1,2,4,5,9,8,12,
     &     2,3,4,6,10,9,13,
     &     1,4,3,8,10,7,14/
      data ifacew /1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     4,6,3,1,12,15,9,13/
!
!     an interface load has already been defined in a prevous step   
!
      if(interfaceload.eq.2) then
        noelset='EXTERNALFACES1235711T'
        do i=22,81
          noelset(i:i)=' '
        enddo
        call cident81(set,noelset,nset,id)
        if(id.eq.0) then
          write(*,*) '*ERROR in calcglobmastsurf: a facial surface'
          write(*,*) '       with the name EXTERNALFACES1235711'
          write(*,*) '       does not exist.'
          call exit(201)
        elseif(set(id).ne.noelset) then
          write(*,*) '*ERROR in calcglobmastsurf: a facial surface'
          write(*,*) '       with the name EXTERNALFACES1235711'
          write(*,*) '       does not exist.'
          call exit(201)
        endif
        imastset=id
        nmastface=iendset(imastset)-istartset(imastset)+1
        return
      endif
!     
!     first occurrence of an interface load in the input deck:    
!     external faces have to be determined    
!     
      kflag=1
      ithree=3
      ifour=4
!
      allocate(ipoface(nk))
      allocate(nodface(5,6*ne))
!     
!     determining the external element faces of the volumetric mesh 
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
      do i=1,6*ne-1
        nodface(5,i)=i+1
      enddo
      do i=1,ne
        if(ipkon(i).lt.0) cycle
        if(lakon(i)(1:3).ne.'C3D') cycle
        indexe=ipkon(i)
        if((lakon(i)(4:4).eq.'2').or.(lakon(i)(4:4).eq.'8')) then
          do j=1,6
            do k=1,4
              nodes(k)=kon(indexe+ifaceq(k,j))
            enddo
            call isortii(nodes,iaux,ifour,kflag)
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
                nodes(k)=kon(indexe+ifacew(k,j))
              enddo
              call isortii(nodes,iaux,ithree,kflag)
            else
              do k=1,4
                nodes(k)=kon(indexe+ifacew(k,j))
              enddo
              call isortii(nodes,iaux,ifour,kflag)
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
!     storing the external faces in a global master surface
!
      noelset='EXTERNALFACES1235711T'
      do i=22,81
        noelset(i:i)=' '
      enddo
      call cident81(set,noelset,nset,id)
      if(id.gt.0) then
        if(set(id).eq.noelset) then
          write(*,*) '*ERROR in calcglobmastsurf: a facial surface'
          write(*,*) '       with the name EXTERNALFACES1235711'
          write(*,*) '       already exists.'
          call exit(201)
        endif
      endif
      nset=nset+1
      do j=nset,id+2,-1
        istartset(j)=istartset(j-1)
        iendset(j)=iendset(j-1)
        set(j)=set(j-1)
      enddo
      set(id+1)=noelset
      istartset(id+1)=nalset+1
!
      do m=1,nk
        index=ipoface(m)
        do
          if(index.eq.0) exit
          nalset=nalset+1     
          i=nodface(3,index)
          j=nodface(4,index)
          ialset(nalset)=10*i+j
          index=nodface(5,index)
        enddo
      enddo
!
      iendset(id+1)=nalset
!
      imastset=id+1
      nmastface=iendset(imastset)-istartset(imastset)+1
!
      deallocate(ipoface)
      deallocate(nodface)
!     
      return
      end
