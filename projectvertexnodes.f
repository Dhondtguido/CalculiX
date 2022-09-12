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
      subroutine projectvertexnodes(ipoed,iexternedg,iedgext,cotet,
     &     nktet,iedg,iexternfa,ifacext,itreated,ilist,isharp,ipofa,
     &     ifac,iedgextfa,ifacexted,co,idimsh,ipoeln,ieln,kontet,c1,
     &     jflag,ibadnodes,nbadnodes,iwrite)
!     
!     1. projects vertex nodes lying on external edges of the unrefined mesh
!        on the parent external edges if these are sharp
!
!     2. projects all other vertex nodes belonging to external faces
!        on neighboring external faces of the unrefined mesh      
!     
      implicit none
!     
      integer ipoed(*),index,i,j,k,iexternedg(*),imasted,iedgext(3,*),
     &     nterms,node1,node2,iedg(3,*),nktet,imastfa,ifacext(6,*),
     &     iexternfa(*),itreated(*),ii,node,isharp(*),idimsh(*),
     &     ilist(*),nlist,ifac(4,*),iedgextfa(2,*),ifacexted(3,*),
     &     id,isol,ipofa(*),ipoeln(*),ieln(2,*),kontet(4,*),jflag,
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
!     loop over all edges: projection on edges from the unrefined mesh 
!     
      loop1: do i=1,nktet
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
              if(idimsh(node1).ne.2) itreated(node1)=1
              if(idimsh(node2).ne.2) itreated(node2)=1
            else
              nterms=2
            endif
!     
!     projection is only needed if the master edge is
!     quadratic
!     
            if(nterms.eq.3) then
!     
!     attach the first end node of the edge
!     
              if(itreated(node1).eq.0) then
                do k=1,3
                  pnode(k)=cotet(k,node1)
                enddo
                call attach_1d(pneigh,pnode,nterms,ratio,dist,xi)
                call checkvol(cotet,node1,pnode,ipoeln,ieln,kontet,c1,
     &               jflag,ibadnodes,nbadnodes,iwrite)
                itreated(node1)=1
              endif
!     
!     attach the other end node of the edge
!     
              if(itreated(node2).eq.0) then
                do k=1,3
                  pnode(k)=cotet(k,node2)
                enddo
                call attach_1d(pneigh,pnode,nterms,ratio,dist,xi)
                call checkvol(cotet,node2,pnode,ipoeln,ieln,kontet,c1,
     &               jflag,ibadnodes,nbadnodes,iwrite)
                itreated(node2)=1
              endif
            else
              itreated(node1)=1
              itreated(node2)=1
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
      loop2: do i=1,nktet
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
!     external face; treat the vertex nodes of the face
!     
          do ii=1,3
            node=ifac(ii,index)
            if(itreated(node).ne.0) cycle
!     
!     parent face
!     
            imastfa=iexternfa(index)
            ilist(1)=imastfa
            nlist=1
            do k=1,3
              pnode(k)=cotet(k,node)
            enddo
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
                call checkvol(cotet,node,pnodeproj,ipoeln,ieln,kontet,
     &               c1,jflag,ibadnodes,nbadnodes,iwrite)
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
                      call checkvol(cotet,node,pnodeproj,ipoeln,ieln,
     &                     kontet,c1,jflag,ibadnodes,nbadnodes,iwrite)
                      itreated(node)=2
                      exit
                    endif
!                      
                    write(*,*) '*WARNING in projectvertexnodes:'
                    write(*,*) '         face in list revisited'
                    write(*,*) '         no projection for node',
     &                   node
                    write(*,*)
                    itreated(node)=2
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
          enddo
!     
          index=ifac(4,index)
        enddo
      enddo loop2
!     
      return
      end
