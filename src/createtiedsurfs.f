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
      subroutine createtiedsurfs(nodface,ipoface,set,istartset,
     &     iendset,ialset,tieset,inomat,ne,ipkon,lakon,kon,ntie,
     &     tietol,nalset,nk,nset,iactive)
!     
      implicit none
!     
!     creates ties for the surfaces in between domains of an
!     electromagnetic calculation
!     
      character*8 lakon(*)
      character*81 set(*),tieset(3,*),surfset
!     
      integer ipkon(*),kon(*),ne,nodface(5,*),ipoface(*),istartset(*),
     &     iendset(*),ialset(*),inomat(*),ithree,ifour,ifaceq(8,6),iset,
     &     ifacet(6,4),ifacew(8,5),ifree,ifreenew,index,indexold,id,
     &     i,j,k,iactive(3),ntie,nodes(4),nalset,nk,nset,indexe
!     
      real*8 tietol(4,*)
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
      data ifacew /1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     4,6,3,1,12,15,9,13/
!     
      ithree=3
      ifour=4
!     
!     determining the external element faces of the electromagnetic mesh 
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
            call insertsorti(nodes,ifour)
c     call isortii(nodes,iaux,ifour,kflag)
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
c     call isortii(nodes,iaux,ithree,kflag)
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
              call insertsorti(nodes,ithree)
c     call isortii(nodes,iaux,ithree,kflag)
            else
              do k=1,4
                nodes(k)=kon(indexe+ifacew(k,j))
              enddo
              call insertsorti(nodes,ifour)
c     call isortii(nodes,iaux,ifour,kflag)
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
!     boundary of domain 1 (phi-domain)
!     
      surfset(1:22)='ELECTROMAGNETICSZONE1T'
      do i=23,81
        surfset(i:i)=' '
      enddo
      call cident81(set,surfset,nset,id)
      nset=nset+1
      do j=nset,id+2,-1
        istartset(j)=istartset(j-1)
        iendset(j)=iendset(j-1)
        set(j)=set(j-1)
      enddo
      iset=id+1
ccc   to remove start          
c     iset=nset
ccc   to remove end
      istartset(iset)=nalset+1
      do i=1,nk
        if(ipoface(i).eq.0) cycle
        if(inomat(i).ne.1) cycle
        index=ipoface(i)
        do
          if(index.eq.0) exit
          nalset=nalset+1
          ialset(nalset)=10*nodface(3,index)+nodface(4,index)
          index=nodface(5,index)
        enddo
      enddo
      if(istartset(iset).gt.nalset) then
        iactive(1)=0
        do j=iset+1,nset
          istartset(j-1)=istartset(j)
          iendset(j-1)=iendset(j)
          set(j-1)=set(j)
        enddo
        nset=nset-1
      else
        iactive(1)=iset
ccc   to remove start          
c     set(iset)(1:22)='ELECTROMAGNETICSZONE1T'
c     do i=23,81
c     set(iset)(i:i)=' '
c     enddo
ccc   to remove end
        set(iset)=surfset
        iendset(iset)=nalset
      endif
!     
!     boundary of domain 2 (A-V-domain)
!     
      surfset(1:22)='ELECTROMAGNETICSZONE2T'
      do i=23,81
        surfset(i:i)=' '
      enddo
      call cident81(set,surfset,nset,id)
      nset=nset+1
      do j=nset,id+2,-1
        istartset(j)=istartset(j-1)
        iendset(j)=iendset(j-1)
        set(j)=set(j-1)
      enddo
      iset=id+1
ccc   to remove start          
c     iset=nset
ccc   to remove end
      istartset(iset)=nalset+1
      do i=1,nk
        if(ipoface(i).eq.0) cycle
        if(inomat(i).ne.2) cycle
        index=ipoface(i)
        do
          if(index.eq.0) exit
          nalset=nalset+1
          ialset(nalset)=10*nodface(3,index)+nodface(4,index)
          index=nodface(5,index)
        enddo
      enddo
      if(istartset(iset).gt.nalset) then
        iactive(2)=0
        do j=iset+1,nset
          istartset(j-1)=istartset(j)
          iendset(j-1)=iendset(j)
          set(j-1)=set(j)
        enddo
        nset=nset-1
      else
        iactive(2)=iset
ccc   to remove start          
c     set(iset)(1:22)='ELECTROMAGNETICSZONE2T'
c     do i=23,81
c     set(iset)(i:i)=' '
c     enddo
ccc   to remove end
        set(iset)=surfset
        iendset(iset)=nalset
      endif
!     
!     boundary of domain 3 (A-domain)
!     
      surfset(1:22)='ELECTROMAGNETICSZONE3T'
      do i=23,81
        surfset(i:i)=' '
      enddo
      call cident81(set,surfset,nset,id)
      nset=nset+1
      do j=nset,id+2,-1
        istartset(j)=istartset(j-1)
        iendset(j)=iendset(j-1)
        set(j)=set(j-1)
      enddo
      iset=id+1
ccc   to remove start          
c     iset=nset
ccc   to remove end
      istartset(iset)=nalset+1
      do i=1,nk
        if(ipoface(i).eq.0) cycle
        if(inomat(i).ne.3) cycle
        index=ipoface(i)
        do
          if(index.eq.0) exit
          nalset=nalset+1
          ialset(nalset)=10*nodface(3,index)+nodface(4,index)
          index=nodface(5,index)
        enddo
      enddo
      if(istartset(iset).gt.nalset) then
        iactive(3)=0
        do j=iset+1,nset
          istartset(j-1)=istartset(j)
          iendset(j-1)=iendset(j)
          set(j-1)=set(j)
        enddo
        nset=nset-1
      else
        iactive(3)=iset
ccc   to remove start          
c     set(iset)(1:22)='ELECTROMAGNETICSZONE3T'
c     do i=23,81
c     set(iset)(i:i)=' '
c     enddo
ccc   to remove end
        set(iset)=surfset
        iendset(iset)=nalset
      endif
!     
!     create contact ties between the domains
!     
      do i=1,3
        do j=1,3
          if((i.eq.j).or.((i.eq.3).and.(j.eq.2))) cycle
          if(iactive(i)*iactive(j).gt.0) then
            ntie=ntie+1
            tietol(1,ntie)=0.d0
            tietol(2,ntie)=0.d0
            write(tieset(1,ntie)(1:1),'(i1)')i
            write(tieset(1,ntie)(2:2),'(i1)')j
            do k=3,80
              tieset(1,ntie)(k:k)=' '
            enddo
            tieset(1,ntie)(81:81)='E'
!     
!     slave set 
!     this set does not have to be identified as nodal
!     or facial at this stage
!     
            tieset(2,ntie)(1:20)='ELECTROMAGNETICSZONE'
            write(tieset(2,ntie)(21:21),'(i1)')i
            do k=22,81
              tieset(2,ntie)(k:k)=' '
            enddo
!     
!     master set
!     
            tieset(3,ntie)(1:22)='ELECTROMAGNETICSZONE T'
            write(tieset(3,ntie)(21:21),'(i1)')j
            do k=23,81
              tieset(3,ntie)(k:k)=' '
            enddo
          endif
        enddo
      enddo
!     
      return
      end
