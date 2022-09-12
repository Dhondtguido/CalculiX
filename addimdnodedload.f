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
      subroutine addimdnodedload(nelemload,sideload,ipkon,kon,lakon,
     &  iload,imdnode,nmdnode,ikmpc,ilmpc,ipompc,nodempc,nmpc,imddof,
     &  nmddof,nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,ikboun,nboun,
     &  ilboun,ithermal)
!
!     adds the nodes belonging to a user-defined facial load to imdnode 
!     (needed in dyna.c and steadystate.c)
!
      implicit none
!
      character*8 lakon(*),lakonl
      character*20 sideload(*)
!
      integer nelemload(2,*),ipkon(*),kon(*),iload,ii,nopes,node,
     &  indexe,
     &  ifaceq(8,6),ifacew(8,5),ifacet(6,4),ig,ielem,nope,imdnode(*),
     &  nmdnode,ikmpc(*),
     &  ilmpc(*),ipompc(*),nodempc(3,*),nmpc,imddof(*),nmddof,
     &  mi(*),nactdof(0:mi(2),*),imdmpc(*),nmdmpc,imdboun(*),nmdboun,
     &  ikboun(*),nboun,ilboun(*),ithermal,k
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
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
      if((sideload(iload)(1:1).eq.'P').and.
     &   (sideload(iload)(3:4).eq.'NU')) then
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
!           user-defined load
!
            if(sideload(iload)(3:4).eq.'NU') then
               call addimd(imdnode,nmdnode,node)
!
!        add the degrees of freedom corresponding to the node
!
               if(ithermal.ne.2) then
                  do k=1,3
                     call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                    nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                    nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &                    ikboun,nboun,ilboun)
                  enddo
               else
                  k=0
                  call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                 nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                 nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,ikboun,
     &                 nboun,ilboun)
               endif
            endif
!
         enddo
      elseif((sideload(iload)(1:1).eq.'B').and.
     &       (sideload(iload)(3:4).eq.'NU')) then
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
!     user-defined load
!     
c            if(sideload(iload)(3:4).eq.'NU') then
               call addimd(imdnode,nmdnode,node)
!     
!     add the degrees of freedom corresponding to the node
!     
               if(ithermal.ne.2) then
                  do k=1,3
                     call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                    nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                    nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &                    ikboun,nboun,ilboun)
                  enddo
               else
                  k=0
                  call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &                 nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &                 nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,ikboun,
     &                 nboun,ilboun)
               endif
c            endif
         enddo
!
      endif
!     
      return
      end

