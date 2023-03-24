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
      subroutine addizdofdload(nelemload,sideload,ipkon,kon,lakon,
     &  nactdof,izdof,nzdof,mi,iload,iznode,nznode,nk,imdnode,nmdnode)
!
!     adds the nodes belonging to a facial load to iznode, izdof
!     and to imdnode if user-defined load
!     (needed in dyna.c and steadystate.c)
!
      implicit none
!
      character*8 lakon(*),lakonl
      character*20 sideload(*)
!
      integer mi(*),nelemload(2,*),ipkon(*),kon(*),nactdof(0:mi(2),*),
     &  izdof(*),nzdof,iload,j,ii,nopes,node,indexe,jdof,ifaceq(8,6),
     &  ifacew(8,5),ifacet(7,4),ig,ielem,nope,iznode(*),nznode,
     &  nodebasis,nk,imdnode(*),nmdnode
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,11,
     &             1,2,4,5,9,8,12,
     &             2,3,4,6,10,9,13,
     &             1,4,3,8,10,7,14/
      data ifacew /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
!
      ielem=nelemload(1,iload)
      lakonl=lakon(ielem)
      indexe=ipkon(ielem)
!
      if(sideload(iload)(1:1).eq.'P') then
         read(sideload(iload)(2:2),'(i1)') ig
!
!        surface pressure: number of nodes belonging to the face
!
         if(lakonl(4:4).eq.'2') then
            nopes=8
         elseif(lakonl(4:4).eq.'8') then
            nopes=4
         elseif(lakonl(4:5).eq.'10') then
            nopes=6
         elseif(lakonl(4:4).eq.'4') then
            nopes=3
         elseif(lakonl(4:5).eq.'15') then
            if(ig.le.2) then
               nopes=6
            else
               nopes=8
            endif
         elseif(lakonl(4:4).eq.'6') then
            if(ig.le.2) then
               nopes=3
            else
               nopes=4
            endif
         endif
!     
         do ii=1,nopes
            if((lakonl(4:4).eq.'2').or.(lakonl(4:4).eq.'8')) then
               node=kon(indexe+ifaceq(ii,ig))
            elseif((lakonl(4:5).eq.'10').or.(lakonl(4:4).eq.'4')) then
               node=kon(indexe+ifacet(ii,ig))
            elseif((lakonl(4:5).eq.'15').or.(lakonl(4:4).eq.'6')) then
               node=kon(indexe+ifacew(ii,ig))
            endif
!
!           adding the nodes in the basis sector to iznode
!
            nodebasis=mod(node,nk)
            call addimd(iznode,nznode,nodebasis)
!
!           user-defined load
!
            if(sideload(iload)(3:4).eq.'NU') then
               call addimd(imdnode,nmdnode,node)
            endif
!
            do j=1,3
!     
!              C-convention!
!     
               jdof=nactdof(j,node)-1
               if(jdof.gt.0) call addimd(izdof,nzdof,jdof)
c               if(jdof.gt.-1) call addimd(izdof,nzdof,jdof)
            enddo
         enddo
      elseif(sideload(iload)(1:1).eq.'B') then
!
!        volumetric load; number of nodes in the element
!   
         if(lakonl(4:4).eq.'2') then
            nope=20
         elseif(lakonl(4:4).eq.'8') then
            nope=8
         elseif(lakonl(4:5).eq.'10') then
            nope=10
         elseif(lakonl(4:4).eq.'4') then
            nope=4
         elseif(lakonl(4:5).eq.'15') then
            nope=15
         elseif(lakonl(4:4).eq.'6') then
            nope=6
         endif
!
         do ii=1,nope
            node=kon(indexe+ii)
!     
!     adding the nodes in the basis sector to iznode
!     
            nodebasis=mod(node,nk)
            call addimd(iznode,nznode,nodebasis)
!     
!     user-defined load
!     
            if(sideload(iload)(3:4).eq.'NU') then
               call addimd(imdnode,nmdnode,node)
            endif
!     
            do j=1,3
!     
!     C-convention!
!     
               jdof=nactdof(j,node)-1
               if(jdof.gt.0) call addimd(izdof,nzdof,jdof)
c               if(jdof.gt.-1) call addimd(izdof,nzdof,jdof)
            enddo
         enddo
      endif
!     
      return
      end

