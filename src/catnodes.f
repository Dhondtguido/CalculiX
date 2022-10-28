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
      subroutine catnodes(ifreenn,inn,iponn,iedg,ipoed,nktet_,
     &  iexternnode,idimsh,isharp,iexternedg)
!
!     catalogues the nodes according to the number of sharp edges
!     in the unrefined mesh they belong to. This number is stored
!     in idimsh(node)
!
!     stores the neighboring nodes of node which are to be used 
!     for smoothing in fields iponn and inn (cf. Documentation)
!     
      implicit none
!
      integer ifreenn,inn(2,*),iponn(*),iedg(3,*),ipoed(*),nktet_,
     &     iexternnode(*),idimsh(*),isharp(*),iexternedg(*),node1,
     &     node2,index,imasted,i
!
!     calculating idimsh
!
      loop1: do i=1,nktet_
        index=ipoed(i)
        do
          if(index.eq.0) cycle loop1
          if(iexternedg(index).gt.0) then
!
!           edge is part of an external edge of the unrefined mesh
!
            imasted=iexternedg(index)
            if(isharp(imasted).eq.1) then
              node1=iedg(1,index)
              node2=iedg(2,index)
              idimsh(node1)=idimsh(node1)+1
              idimsh(node2)=idimsh(node2)+1
            endif
          endif
          index=iedg(3,index)
        enddo
      enddo loop1
!
!     determining iponn and inn
!
      loop2: do i=1,nktet_
        index=ipoed(i)
        do
          if(index.eq.0) cycle loop2
          node1=iedg(1,index)
          node2=iedg(2,index)
!
!         node1
!
!         case 1: internal node
!
          if(iexternnode(node1).eq.0) then
            inn(1,ifreenn)=node2
            inn(2,ifreenn)=iponn(node1)
            iponn(node1)=ifreenn
            ifreenn=ifreenn+1
!
!     case 2: external node, no sharp edges attached and
!             actual edge is external            
!
          elseif(idimsh(node1).eq.0) then
            if(iexternedg(index).ne.0) then
              inn(1,ifreenn)=node2
              inn(2,ifreenn)=iponn(node1)
              iponn(node1)=ifreenn
              ifreenn=ifreenn+1
            endif
!
!     case 3: external node, 2 sharp edges attached,
!             actual edge is external and sharp           
!
          elseif(idimsh(node1).eq.2) then
            if((iexternedg(index).gt.0).and.
     &           (isharp(iexternedg(index)).eq.1)) then
              inn(1,ifreenn)=node2
              inn(2,ifreenn)=iponn(node1)
              iponn(node1)=ifreenn
              ifreenn=ifreenn+1
            endif
          endif
!
!         node2
!
!         case 1: internal node
!
          if(iexternnode(node2).eq.0) then
            inn(1,ifreenn)=node1
            inn(2,ifreenn)=iponn(node2)
            iponn(node2)=ifreenn
            ifreenn=ifreenn+1
!
!     case 2: external node, no sharp edges attached and
!             actual edge is external            
!
          elseif(idimsh(node2).eq.0) then
            if(iexternedg(index).ne.0) then
              inn(1,ifreenn)=node1
              inn(2,ifreenn)=iponn(node2)
              iponn(node2)=ifreenn
              ifreenn=ifreenn+1
            endif
!
!     case 3: external node, 2 sharp edges attached,
!             actual edge is external and sharp           
!
          elseif(idimsh(node2).eq.2) then
            if((iexternedg(index).gt.0).and.
     &           (isharp(iexternedg(index)).eq.1)) then
              inn(1,ifreenn)=node1
              inn(2,ifreenn)=iponn(node2)
              iponn(node2)=ifreenn
              ifreenn=ifreenn+1
            endif
          endif
          index=iedg(3,index)
        enddo
      enddo loop2
!
      return
      end
