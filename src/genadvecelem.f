!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine genadvecelem(inodesd,ipkon,ne,lakon,kon,nload,
     &     sideload,nelemload,nkon,network)
!
!     generates elements simulating advection between network
!     elements and structural faces
!
      implicit none
!
      character*8 lakon(*)
      character*20 sideload(*)
!
      integer inodesd(*),nnodesd,ipkon(*),ne,i,j,indexe,node,id,kon(*),
     &  nload,ifaceq(8,6),ifacet(6,4),ifacew1(4,5),ifacew2(8,5),
     &  nodef(8),nelemload(2,*),nope,jface,k,nopes,nkon,nelem,
     &  nface,network
!
!     nodes per face for hex elements
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
!
!     nodes per face for tet elements
!
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
!
!     nodes per face for linear wedge elements
!
      data ifacew1 /1,3,2,0,
     &             4,5,6,0,
     &             1,2,5,4,
     &             2,3,6,5,
     &             3,1,4,6/
!
!     nodes per face for quadratic wedge elements
!
      data ifacew2 /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             3,1,4,6,9,13,12,15/
!
!     catalogueing the nodes belonging to "Dx"-elements (specific
!     network elements, for which "D" is followed by some
!     specification such as restrictor or vortex) unless the
!     network is declared by the user as a thermal network
!     (parameter THERMAL NETWORK on the *STEP card; network=1)
!
      nnodesd=0
      do i=1,ne
         if(ipkon(i).lt.0) cycle
         if(lakon(i)(1:1).ne.'D') cycle
         if((lakon(i)(2:2).eq.' ').or.(network.eq.1)) cycle
         indexe=ipkon(i)
         do j=1,3,2
            node=kon(indexe+j)
            if(node.eq.0) cycle
            call nident(inodesd,node,nnodesd,id)
            if(id.gt.0) then
               if(inodesd(id).eq.node) cycle
            endif
            nnodesd=nnodesd+1
            do k=nnodesd,id+2,-1
               inodesd(k)=inodesd(k-1)
            enddo
            inodesd(id+1)=node
         enddo
      enddo
!
!     check whether forced convection film condition has generic
!     network nodes as sink nodes
!
      do i=1,nload
         if((sideload(i)(1:1).eq.'F').and.
     &      (sideload(i)(3:4).eq.'FC')) then
            node=nelemload(2,i)
            call nident(inodesd,node,nnodesd,id)
            if(id.gt.0) then
               if(inodesd(id).eq.node) cycle
            endif
            nelem=nelemload(1,i)
            indexe=ipkon(nelem)
            if(indexe.lt.0) cycle
!
!           new advection element is generated
!
            ne=ne+1
            ipkon(ne)=nkon
            lakon(ne)(1:7)='ESPRNGF'
            read(sideload(i)(2:2),'(i1)') jface
!     
            if(lakon(nelem)(4:4).eq.'2') then
               nopes=8
               nface=6
            elseif(lakon(nelem)(4:4).eq.'8') then
               nopes=4
               nface=6
            elseif(lakon(nelem)(4:5).eq.'10') then
               nopes=6
               nface=4
            elseif(lakon(nelem)(4:4).eq.'4') then
               nopes=3
               nface=4
            elseif(lakon(nelem)(4:5).eq.'15') then
               if(jface.le.2) then
                  nopes=6
               else
                  nopes=8
               endif
               nface=5
               nope=15
            elseif(lakon(nelem)(4:4).eq.'6') then
               if(jface.le.2) then
                  nopes=3
               else
                  nopes=4
               endif
               nface=5
               nope=6
            else
               cycle
            endif
!     
!           determining the nodes of the face
!     
            if(nface.eq.4) then
               do k=1,nopes
                  nodef(k)=kon(indexe+ifacet(k,jface))
               enddo
            elseif(nface.eq.5) then
               if(nope.eq.6) then
                  do k=1,nopes
                     nodef(k)=kon(indexe+ifacew1(k,jface))
                  enddo
               elseif(nope.eq.15) then
                  do k=1,nopes
                     nodef(k)=kon(indexe+ifacew2(k,jface))
                  enddo
               endif
            elseif(nface.eq.6) then
               do k=1,nopes
                  nodef(k)=kon(indexe+ifaceq(k,jface))
               enddo
            endif
!     
            do k=1,nopes
               kon(nkon+k)=nodef(k)
            enddo
            nkon=nkon+nopes+1
            kon(nkon)=node
!     
            write(lakon(ne)(8:8),'(i1)') nopes
!     
!           copying the loading
!     
            nload=nload+1
            nelemload(1,nload)=ne
!     
!           pointer to the original load
!     
            nelemload(2,nload)=i
            sideload(nload)='                    '
!     
!           deactivating the original load
!     
            sideload(i)(1:1)=' '
!     
         endif
      enddo
!
      return
      end

