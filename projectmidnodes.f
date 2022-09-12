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
      subroutine projectmidnodes(nktet_,ipoed,iedgmid,iexternedg,
     &     iedgext,cotet,nktet,iedg,iexternfa,ifacext,itreated,
     &     ilist,isharp,ipofa,ifac,iedgextfa,ifacexted,jfix,co,idimsh,
     &     ipoeled,ieled,kontet,c1,jflag,iedtet,ibadnodes,nbadnodes,
     &     iwrite)
!     
!     1. projects all midnodes on external edges of the unrefined mesh
!        on the parent external edges if they are sharp
!     
!     2. projects all other midnodes belonging to external faces on
!        neighboring external faces of the unrefined mesh      
!     
      implicit none
!     
      integer nktet_,ipoed(*),index,i,j,k,iedgmid(*),iexternedg(*),
     &     imasted,iedgext(3,*),nterms,node1,node2,iedg(3,*),
     &     nktet,imastfa,ifacext(6,*),iexternfa(*),itreated(*),ii,node,
     &     isharp(*),ilist(*),nlist,ifac(4,*),iedgextfa(2,*),
     &     ifacexted(3,*),id,indexe,isol,ipofa(*),jfix(*),idimsh(*),
     &     ipoeled(*),ieled(2,*),kontet(4,*),jflag,iedtet(6,*),
     &     ibadnodes(*),nbadnodes,iwrite
!     
      real*8 pneigh(3,9),cotet(3,*),pnode(3),ratio(9),dist,xi,et,
     &     pnodeproj(3),co(3,*),c1
!
      nbadnodes=0
!
!     set all nodes on "non-treated"
!      
      do i=1,nktet
        itreated(i)=0
      enddo
!     
!     loop over all edges: projection on sharp edges from the unrefined mesh 
!     
      loop1: do i=1,nktet_
        index=ipoed(i)
!     
        do
          if(index.eq.0) cycle loop1
!     
          if(iexternedg(index).gt.0) then
!     
!     recovering master edge
!     
            imasted=iexternedg(index)
!
!     if the master edge is not sharp: loop
!
            if(isharp(imasted).ne.1) then
              index=iedg(3,index)
              cycle
            endif
!     
!     edge is a subset of an external sharp edge
!     of the unrefined mesh
!     
!     end nodes belonging to the edge
!     
            node1=iedg(1,index)
            node2=iedg(2,index)
!            
            if(iedgext(2,imasted).ne.0) then
              nterms=3
              do j=1,3
                do k=1,3
                  pneigh(k,j)=cotet(k,iedgext(j,imasted))
                enddo
              enddo
            else
              nterms=2
            endif
!     
!     edge is a subset of an external sharp edge
!     of the unrefined mesh
!     
            if((jfix(node1).ne.1).or.(jfix(node2).ne.1)) then
              node=iedgmid(index)
!     
              if(nterms.eq.3) then
                do k=1,3
                  pnode(k)=cotet(k,node)
                enddo
!     
!     projection for quadratic master edges
!     
                call attach_1d(pneigh,pnode,nterms,ratio,dist,xi)
                call checkjac(cotet,node,pnode,kontet,c1,jflag,
     &               iedtet,iedgmid,ipoeled,ieled,index,ibadnodes,
     &               nbadnodes,iwrite)
!     
                itreated(node)=1
              endif
            endif
          endif
!     
          index=iedg(3,index)
        enddo
      enddo loop1
!     
!     projection of the external nodes not treated yet onto the
!     faces of the unrefined mesh
!     
      loop2: do i=1,nktet_
!     
!     loop over all faces of the refined mesh
!     
        index=ipofa(i)
        do
          if(index.eq.0) cycle loop2
!     
!     if no external face: loop
!     
          if(iexternfa(index).le.0) then
            index=ifac(4,index)
            cycle
          endif
!     
!     external face; treat the midnodes of the face
!     
          do ii=1,3
            if(ii.eq.1) then
              node1=ifac(1,index)
              node2=ifac(2,index)
            elseif(ii.eq.2) then
              node1=ifac(2,index)
              node2=ifac(3,index)
            else
              node1=ifac(1,index)
              node2=ifac(3,index)
            endif
!     
!     determine the number of the edge
!     
            if(node1.gt.node2) then
              node=node2
              node2=node1
              node1=node
            endif
            indexe=ipoed(node1)
            do
              if(iedg(2,indexe).eq.node2) exit
              indexe=iedg(3,indexe)
            enddo
!     
!     check whether a middle node was treated
!     if so, nothing to do -> check the next edge
!     
            if(itreated(iedgmid(indexe)).ne.0) cycle
!     
!     project the middle node
!
            if((jfix(node1).eq.1).and.(jfix(node2).eq.1)) then
              cycle
            else
              node=iedgmid(indexe)
              do k=1,3
                pnode(k)=cotet(k,node)
              enddo
            endif
!     
!     parent face
!     
            imastfa=iexternfa(index)
            ilist(1)=imastfa
            nlist=1
!     
!     start the loop looking for the correct face; 
!     starting with the parent face 
!     
            do
              if(ifacext(4,imastfa).ne.0) then
                nterms=6
                do j=1,6
                  do k=1,3
                    pneigh(k,j)=co(k,ifacext(j,imastfa))
                  enddo
                enddo
              else
                nterms=3
                do j=1,3
                  do k=1,3
                    pneigh(k,j)=co(k,ifacext(j,imastfa))
                  enddo
                enddo
              endif
!     
              do k=1,3
                pnodeproj(k)=pnode(k)
              enddo
              call attach_2d(pneigh,pnodeproj,nterms,ratio,dist,
     &             xi,et)
!     
!     check whether this face is the correct one;
!     if dabs(xi)=1 or dabs(et)=1 or xi+et=0 this
!     may not be the case
!     
!     the solution is found (isol=1) unless proved
!     otherwise
!     
              isol=1
!     
              if((dabs(et).lt.1.d-10).and.
     &           (dabs(xi).gt.1.d-10)) then
!     
!     take neighboring face across edge 1-2 unless sharp
!     
                imasted=ifacexted(1,imastfa)
                if(isharp(imasted).eq.0) then
                  if(iedgextfa(1,imasted).eq.imastfa) then
                    imastfa=iedgextfa(2,imasted)
                  else
                    imastfa=iedgextfa(1,imasted)
                  endif
                  isol=0
                endif
!     
              elseif((dabs(xi+et-1.d0).lt.1.d-10).and.
     &               (dabs(et).gt.1.d-10)) then
!     
!     take neighboring face across edge 2-3 unless sharp
!     
                imasted=ifacexted(2,imastfa)
                if(isharp(imasted).eq.0) then
                  if(iedgextfa(1,imasted).eq.imastfa) then
                    imastfa=iedgextfa(2,imasted)
                  else
                    imastfa=iedgextfa(1,imasted)
                  endif
                  isol=0
                endif
!     
              elseif((dabs(xi).lt.1.d-10).and.
     &               (dabs(xi+et-1.d0).gt.1.d-10)) then
!     
!     take neighboring face across edge 3-1 unless sharp
!     
                imasted=ifacexted(3,imastfa)
                if(isharp(imasted).eq.0) then
                  if(iedgextfa(1,imasted).eq.imastfa) then
                    imastfa=iedgextfa(2,imasted)
                  else
                    imastfa=iedgextfa(1,imasted)
                  endif
                  isol=0
                endif
              endif
!     
!     if solution is found: copy projected coordinates
!     else continue with a neighbor
!     
              if(isol.eq.1) then
                call checkjac(cotet,node,pnodeproj,kontet,c1,jflag,
     &               iedtet,iedgmid,ipoeled,ieled,indexe,ibadnodes,
     &               nbadnodes,iwrite)
                itreated(node)=2
                exit
              else
!     
!     update list; exit if an element in the list is
!     revisited
!     
                call nident(ilist,imastfa,nlist,id)
                if(id.gt.0) then
                  if(ilist(id).eq.imastfa) then
                    if(nlist.eq.2) then
                      isol=1
                      call checkjac(cotet,node,pnodeproj,kontet,c1,
     &                     jflag,iedtet,iedgmid,ipoeled,ieled,indexe,
     &                     ibadnodes,nbadnodes,iwrite)
                      itreated(node)=2
                      exit
                    endif
!                      
                    write(*,*) '*WARNING in projectmidnodes:'
                    write(*,*) '         face in list revisited'
                    write(*,*) '         no projection for node',
     &                   node
                    write(*,*)
                    exit
                  endif
                endif
                nlist=nlist+1
                do k=nlist,id+2,-1
                  ilist(k)=ilist(k-1)
                enddo
                ilist(id+1)=imastfa
!     
              endif
            enddo
!     
          enddo
!     
!     treat the next face
!     
          index=ifac(4,index)
        enddo
      enddo loop2
!     
      return
      end
