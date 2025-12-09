!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine catmastnodes(lakon,ipkon,kon,istartset,iendset,ialset,
     &     imastnode,nmasts,imastset)
!     
!     catalogueing the master nodes
!     
      implicit none
!     
      character*8 lakon(*)
!     
      logical exist
!     
      integer j,k,l,istartset(*),iendset(*),ialset(*),
     &     ifacem,nelemm,nmasts,jfacem,indexe,nopem,ipkon(*),kon(*),id,
     &     ifaceq(8,6),ifacet(6,4),ifacew1(4,5),ifacew2(8,5),node,
     &     imastnode(*),imastset
!     
!     nodes per face for hex elements
!     
      data ifaceq /4,3,2,1,11,10,9,12,
     &     5,6,7,8,13,14,15,16,
     &     1,2,6,5,9,18,13,17,
     &     2,3,7,6,10,19,14,18,
     &     3,4,8,7,11,20,15,19,
     &     4,1,5,8,12,17,16,20/
!     
!     nodes per face for tet elements
!     
      data ifacet /1,3,2,7,6,5,
     &     1,2,4,5,9,8,
     &     2,3,4,6,10,9,
     &     1,4,3,8,10,7/
!     
!     nodes per face for linear wedge elements
!     
      data ifacew1 /1,3,2,0,
     &     4,5,6,0,
     &     1,2,5,4,
     &     2,3,6,5,
     &     3,1,4,6/
!     
!     nodes per face for quadratic wedge elements
!     
      data ifacew2 /1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     3,1,4,6,9,13,12,15/
!     
      nmasts=0
!     
      do j=istartset(imastset),iendset(imastset)
!     
!     create imastnode and determine nmasts
!     
        ifacem=ialset(j)
        nelemm=int(ifacem/10)
        if(ipkon(nelemm).lt.0) cycle
        jfacem=ifacem-nelemm*10
        indexe=ipkon(nelemm)
!     
        if(lakon(nelemm)(4:5).eq.'20') then
          nopem=8
        elseif(lakon(nelemm)(4:4).eq.'8') then
          nopem=4
        elseif(lakon(nelemm)(4:5).eq.'10') then
          nopem=6
        elseif(lakon(nelemm)(4:4).eq.'4') then
          nopem=3
        endif
!     
        if(lakon(nelemm)(4:4).eq.'6') then
          if(jfacem.le.2) then
            nopem=3
          else
            nopem=4
          endif
        endif
        if(lakon(nelemm)(4:5).eq.'15') then
          if(jfacem.le.2) then
            nopem=6
          else
            nopem=8
          endif
        endif   
!     
        do l=1,nopem
          if((lakon(nelemm)(4:4).eq.'2').or.
     &         (lakon(nelemm)(4:4).eq.'8')) then
            node=kon(indexe+ifaceq(l,jfacem))
          elseif((lakon(nelemm)(4:4).eq.'4').or.
     &           (lakon(nelemm)(4:5).eq.'10')) then
            node=kon(indexe+ifacet(l,jfacem))
          elseif(lakon(nelemm)(4:4).eq.'6') then
            node=kon(indexe+ifacew1(l,jfacem))
          elseif(lakon(nelemm)(4:5).eq.'15') then
            node=kon(indexe+ifacew2(l,jfacem))
          endif
          call nident(imastnode,node,
     &         nmasts,id)
          exist=.false.
          if(id.gt.0) then
            if(imastnode(id).eq.node) then
              exist=.true.
            endif
          endif
          if(exist) cycle
          nmasts=nmasts+1
          do k=nmasts,id+2,-1
            imastnode(k)=imastnode(k-1)
          enddo
          imastnode(id+1)=node
        enddo
!     
      enddo
!     
      return
      end
