!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2021 Guido Dhondt
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
      subroutine boundarymesh(nbounedg,ibounedg,ieled,
     &     ibounel,nbounel,iedg,ne,kontri,ipoed,
     &     ipkon,lakon,ncenter,nkon,kon,mastelnr,ntri)
!     
!     simplifying the part of the crack mesh not connected to the
!     crack boundary
!     
      implicit none
!
      character*8 lakon(*)
!     
      integer iedge,ieled(2,*),nbounedg,mastelnr(*),ntri,
     &     ibounel(*),ielem,nbounel,i,j,k,id,iedg(3,*),
     &     ibounedg(*),ne,iel,node,node1,node2,
     &     kontri(3,*),noperel(2,3),ipoed(*),index,ipkon(*),node1or,
     &     node2or,ncenter,nkon,kon(*),neafterprop
!     
      data noperel /1,2,2,3,3,1/
!
!     present number of elements
!
      neafterprop=ne
!     
!     finding all elements connected to the boundary of the crack
!     before the propagation in the actual increment =: boundary elements
!
!     nbounel: number of boundary elements
!     ibounel(*): element numbers in ascending order
!
!     crack elements which are no boundary elements will be called
!     internal elements
!
      nbounel=0
      do i=1,nbounedg
        iedge=ibounedg(i)
        ielem=ieled(1,iedge)
!
!       check whether the element is already catalogued
!
        call nident(ibounel,ielem,nbounel,id)
        if(id.gt.0) then
          if(ibounel(id).eq.ielem) cycle
        endif
        nbounel=nbounel+1
        do k=nbounel,id+2,-1
          ibounel(k)=ibounel(k-1)
        enddo
        ibounel(id+1)=ielem
      enddo
!
!     look for edges of boundary elements which
!     - are no boundary edges AND
!     - do not belong to any other boundary element
!     these edges form the boundary between the boundary elements
!     and internal elements. 
!
!     creation of new elements by connecting these edges with a
!     fictitious center node
!
      do i=1,nbounel
        iel=ibounel(i)
!
!       loop over the edges belonging to the element
!
        do j=1,3
!
!         nodes belonging to the edge
!
          node1=kontri(noperel(1,j),iel)
          node2=kontri(noperel(2,j),iel)
!
!         sorting in ascending order
!
          node1or=node1
          node2or=node2
          if(node2.lt.node1) then
            node=node1
            node1=node2
            node2=node
          endif
!
!         determining the number of the edge in iedg(3,*)          
!
          index=ipoed(node1)
          do
            if(iedg(2,index).eq.node2) exit
            index=iedg(3,index)
            if(index.eq.0) then
              write(*,*) '*ERROR in boundarymesh: edge database corrupt'
              call exit(201)
            endif
          enddo
          iedge=index
!
!         check whether boundary edge         
!
          if(ieled(2,iedge).eq.0) cycle
!
!         check whether the edge belongs to a boundary element         
!
          do k=1,2
            ielem=ieled(k,iedge)
            if(ielem.eq.iel) cycle
            call nident(ibounel,ielem,nbounel,id)
            if(id.gt.0) then
              if(ibounel(id).eq.ielem) cycle
            endif
!
!           create a new element
!
            ne=ne+1
            lakon(ne)='C3D6  L '
            ipkon(ne)=nkon
            kon(nkon+1)=node2or
            kon(nkon+2)=node1or
            kon(nkon+3)=ncenter
            kon(nkon+4)=node2or
            kon(nkon+5)=node1or
            kon(nkon+6)=ncenter
            kon(nkon+7)=node2or
            kon(nkon+8)=node1or
            kon(nkon+9)=ncenter
            nkon=nkon+9
          enddo
        enddo
      enddo
!
!     deactivating the inner trangles of the crack mesh
!
      do j=1,ntri
!
!       global element number of local S3-element j
!
        i=mastelnr(j)
!
!       element must be active
!
        if(ipkon(i).lt.0) cycle
!
!       element must be a three-noded shell element
!
        if(lakon(i).ne.'C3D6  L ') cycle
!
!       element must not be a boundary element
!
        call nident(ibounel,j,nbounel,id)
        if(id.gt.0) then
          if(ibounel(id).eq.j) cycle
        endif
!
!       deactivate the element
!
        if(ipkon(i).ge.0) ipkon(i)=-ipkon(i)-2
      enddo
!      
      return
      end

