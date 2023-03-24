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
      subroutine adjacentbounodes(ifront,ifrontrel,nfront,iedno,iedg,
     &     nnfront,ifront2,ifrontrel2,ibounnod,nbounnod,istartfront,
     &     iendfront,ibounedg,istartcrackfro,iendcrackfro,ncrack,
     &     ibounnod2,istartcrackbou,iendcrackbou,isubsurffront,iedno2,
     &     stress,stress2,iresort,ieled,kontri,costruc,costruc2,temp,
     &     temp2,nstep,ier)
!     
!     sorting the crack boundary nodes according to adjacency
!     
      implicit none
!     
      integer ifront(*),ifrontrel(*),nfront,iedno(2,*),iedg(3,*),
     &     ichange,nnfront,nfront2,node,noderel,iedge,ibounnod(*),
     &     nbounnod,i,j,k,istartfront(*),iendfront(*),ifront2(*),
     &     ifrontrel2(*),id,ibounedg(*),nodestart,istartcrackfro(*),
     &     iendcrackfro(*),ncrack,iact,nodeprev,noderelprev,kflag,
     &     ibounnod2(*),istartcrackbou(*),iendcrackbou(*),nbounnod2,
     &     isubsurffront(*),isubsurf,iedno2(2,*),iresort(*),itri,
     &     nodenext,ieled(2,*),kontri(3,*),nstep,ier
!
      real*8 stress(6,nstep,*),stress2(6,nstep,*),costruc(3,*),
     &     costruc2(3,*),temp(nstep,*),temp2(nstep,*)
!     
!     nnfront is the number of distinct fronts
!     
      nnfront=0
      nfront2=0
      nbounnod2=0
      ncrack=0
      kflag=2
c      do i=1,nfront
c        write(*,*) 'adjacentbounodes ',i,ifront(i)
c      enddo
!     
      do
!     
!     for all recorded front nodes i the entry in ifrontrel(i) is
!     set to 0. If no zero is left, the job is done.
!     
        ichange=0
!     
        loop: do i=1,nfront
        if(ifrontrel(i).gt.0) then
          ichange=1
!     
!     find the end of a front: start with an arbitrary node
!     
c          write(*,*) 'look for new front '
          node=ifront(i)
c          write(*,*) node
          nodestart=node
          noderel=ifrontrel(i)
!
!     take one of the edges to which the node belongs
!     check whether the nodes belonging to the edge are in the
!     same order as in the element topology of the crack
!
          iedge=ibounedg(iedno(1,noderel))
!
!     crack triangle belonging to the edge
!
          itri=ieled(1,iedge)
!
!     finding the other node of the edge
!
          if(iedg(1,iedge).ne.node) then
            nodenext=iedg(1,iedge)
          else
            nodenext=iedg(2,iedge)
          endif
!
!     node and nodenext also belong to itri; check whether the
!     order in the element topology is nodenext -> node
!     (because the order is reverted as soon as an external node
!     is found (subsurface crack) or the start node is retreated
!     (surface crack)
!     
          do j=1,3
            if(kontri(j,itri).eq.node) exit
          enddo
          k=j+1
          if(k.gt.3) k=1
c          k=mod(j+1,3)
!
!     take the other edge is the order is not the same
!
          if(kontri(k,itri).eq.nodenext) then
            iedge=ibounedg(iedno(2,noderel))
          endif
!
          do
!     
!     find the other node on the edge
!     
            if(iedg(1,iedge).ne.node) then
              node=iedg(1,iedge)
            else
              node=iedg(2,iedge)
            endif
c            write(*,*) node
            if(node.eq.nodestart) then
!     
!     subsurface crack
!     
              call nident(ifront,node,nfront,id)
              ifrontrel(id)=0
              isubsurf=1
              exit
            endif
!     
!     check whether this is a front node
!     
            call nident(ifront,node,nfront,id)
            if(id.gt.0) then
              if(ifront(id).eq.node) then
!     
!     front node! search for other edge adjacent to the node
!     
                noderel=ifrontrel(id)
                if(ibounedg(iedno(1,noderel)).eq.iedge) then
                  iedge=ibounedg(iedno(2,noderel))
                else
                  iedge=ibounedg(iedno(1,noderel))
                endif
                cycle
              endif
            endif
!
            isubsurf=0
            exit
          enddo
!
!     if surface crack:          
!     node is a boundary node which is external to the structure
!     but adjacent to a front node (i.e. a boundary node internal
!     to the structure)
!
!     if subsurface crack:
!     node is an arbitary node belonging to the crack front
!     
!     finding the relative node number
!     
          call nident(ibounnod,node,nbounnod,noderel)
c          write(*,*) 'adjacentbounodes start of new front'
!     
!     defining a new crack
!     
          ncrack=ncrack+1
          nodestart=node
!     
!     defining a new front segment
!     
          nnfront=nnfront+1
          iact=1
!     
!     new node belonging to a segment (a segment contains the internal
!     nodes plus the immediately adjacent external nodes
!     
          nfront2=nfront2+1
          ifront2(nfront2)=node
!
          nbounnod2=nbounnod2+1
          ibounnod2(nbounnod2)=node
!          
c          write(*,*) node
          ifrontrel2(nfront2)=noderel
          istartfront(nnfront)=nfront2
          isubsurffront(nnfront)=isubsurf
          istartcrackfro(ncrack)=nfront2
          istartcrackbou(ncrack)=nbounnod2
!     
!     looking for other nodes belonging to the segment
!     
          do
!     
!     find the other node on the edge
!     
            if(iedg(1,iedge).ne.node) then
              nodeprev=node
              noderelprev=noderel
              node=iedg(1,iedge)
            else
              nodeprev=node
              noderelprev=noderel
              node=iedg(2,iedge)
            endif
!     
!     reached starting node
!     
            if(node.eq.nodestart) then
              iendcrackbou(ncrack)=nbounnod2
              iendfront(nnfront)=nfront2
              exit
            endif
!     
!     check whether this is a front node
!     
            call nident(ifront,node,nfront,id)
            if(id.gt.0) then
              if(ifront(id).eq.node) then
!
!     front node !
!     
!     check whether new front is to be activated
!     
                if(iact.eq.0) then
!     
!     new front
!     
                  iact=1
                  nnfront=nnfront+1
                  isubsurffront(nnfront)=isubsurf
!     
!     first node of new front
!     
                  nfront2=nfront2+1
                  ifront2(nfront2)=nodeprev
                  ifrontrel2(nfront2)=noderelprev
                  istartfront(nnfront)=nfront2
                endif
!     
                noderel=ifrontrel(id)
                if(ibounedg(iedno(1,noderel)).eq.iedge) then
                  iedge=ibounedg(iedno(2,noderel))
                else
                  iedge=ibounedg(iedno(1,noderel))
                endif
!     
                nfront2=nfront2+1
                ifront2(nfront2)=node
!
                nbounnod2=nbounnod2+1
                ibounnod2(nbounnod2)=node
!                
c                write(*,*) node
                ifrontrel2(nfront2)=noderel
!     
!     setting the entry in frontrel for a node which was treated
!     to zero                  
!     
                ifrontrel(id)=0
                cycle
              endif
            endif
!
!     no front node!
!     
            call nident(ibounnod,node,nbounnod,noderel)
            if(ibounedg(iedno(1,noderel)).eq.iedge) then
              iedge=ibounedg(iedno(2,noderel))
            else
              iedge=ibounedg(iedno(1,noderel))
            endif
!
            nbounnod2=nbounnod2+1
            ibounnod2(nbounnod2)=node
!     
            if(iact.eq.1) then
!              
!     last node of the segment (external boundary node immediately
!     adjacent to an internal node)
!     
              nfront2=nfront2+1
              ifront2(nfront2)=node
c              write(*,*) node
              ifrontrel2(nfront2)=noderel
!     
!     finishing the actual segment
!     
              iendfront(nnfront)=nfront2
              iact=0
            endif
          enddo
!     
!     end of crack
!     
          iendcrackfro(ncrack)=iendfront(nnfront)
        endif
      enddo loop                ! end loop over the front nodes  
!     
!     if no front node is left untreated, leave
!     
      if(ichange.eq.0) exit
      enddo
!
!     adjust the dependent fields of ibounnod
!
      do i=1,nbounnod2
        node=ibounnod2(i)
        call nident(ibounnod,node,nbounnod,id)
        do j=1,2
          iedno2(j,i)=iedno(j,id)
        enddo
        do k=1,nstep
          temp2(k,i)=temp(k,id)
        enddo
        do k=1,nstep
          do j=1,6
            stress2(j,k,i)=stress(j,k,id)
          enddo
        enddo
        do j=1,3
          costruc2(j,i)=costruc(j,id)
        enddo
      enddo
!
!     copy the ibounnod2 field and its dependent fields
!
      do i=1,nbounnod2
        ibounnod(i)=ibounnod2(i)
        do j=1,2
          iedno(j,i)=iedno2(j,i)
        enddo
        do k=1,nstep
          temp(k,i)=temp2(k,i)
        enddo
        do k=1,nstep
          do j=1,6
            stress(j,k,i)=stress2(j,k,i)
          enddo
        enddo
        do j=1,3
          costruc(j,i)=costruc2(j,i)
        enddo
      enddo
!
!     calculating in which order ibounnod was resorted
!      
      do i=1,nbounnod2
        iresort(i)=i
      enddo
!
      nfront=nfront
      call isortii(ibounnod2,iresort,nbounnod2,kflag)
!     
!     copy the temporary fields (ifront2 and ifrontrel2)
!     
      do i=1,nfront2
        ifront(i)=ifront2(i)
        ifrontrel(i)=iresort(ifrontrel2(i))
      enddo
      nfront=nfront2
!
c      write(*,*)
c      write(*,*) 'adjacentbounodes: front nodes'
c      write(*,*)
c      do i=1,nfront
c        write(*,*) i,ifront(i)
c      enddo
!     
      if(nfront.eq.0) then
        write(*,*) '*ERROR in adjacentbounodes: no front node found'
        ier=1
c         call exit(201)
       endif
!       
      return
      end

