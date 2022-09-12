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
      subroutine nodestiedface(tieset,ntie,ipkon,kon,
     &  lakon,set,istartset,iendset,ialset,nset,ifaceslave,
     &  istartfield,iendfield,ifield,nconf,ncone,kind)
!
!     identifies slave nodes in tied slave faces
!
      implicit none
!
      character*1 kind
      character*8 lakon(*)
      character*81 tieset(3,*),slavset,set(*)
!
      integer ntie,nset,istartset(*),iendset(*),ialset(*),ifree,
     &  ipkon(*),kon(*),node,ifaceslave(*),i,j,k,l,
     &  ifaceq(8,6),ifacet(6,4),ilength,id,ncone,
     &  ifacew1(4,5),ifacew2(8,5),nelem,jface,indexe,
     &  nnodelem,nface,nope,nodef(8),
     &  ifield(*),istartfield(*),iendfield(*),nconf
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
      ifree=1
!
      do i=1,ntie
         ilength=0
         if(tieset(1,i)(81:81).ne.kind) cycle
         if(ifaceslave(i).eq.0) cycle
         slavset=tieset(2,i)
c         do j=1,nset
c            if(set(j).eq.slavset) exit
c         enddo
         call cident81(set,slavset,nset,id)
         j=nset+1
         if(id.gt.0) then
           if(slavset.eq.set(id)) then
             j=id
           endif
         endif
!
         istartfield(i)=ifree
         do j=istartset(j),iendset(j)
            ncone=ncone-1
            nelem=int(ialset(j)/10.)
            jface=ialset(j)-10*nelem
!
            indexe=ipkon(nelem)
            if(lakon(nelem)(4:4).eq.'2') then
               nnodelem=8
               nface=6
            elseif(lakon(nelem)(4:4).eq.'8') then
               nnodelem=4
               nface=6
            elseif(lakon(nelem)(4:5).eq.'10') then
               nnodelem=6
               nface=4
            elseif(lakon(nelem)(4:4).eq.'4') then
               nnodelem=3
               nface=4
            elseif(lakon(nelem)(4:5).eq.'15') then
               if(jface.le.2) then
                  nnodelem=6
               else
                  nnodelem=8
               endif
               nface=5
               nope=15
            elseif(lakon(nelem)(4:4).eq.'6') then
               if(jface.le.2) then
                  nnodelem=3
               else
                  nnodelem=4
               endif
               nface=5
               nope=6
            else
               cycle
            endif
!     
!     determining the nodes of the face
!     
            if(nface.eq.4) then
               do k=1,nnodelem
                  nodef(k)=kon(indexe+ifacet(k,jface))
               enddo
            elseif(nface.eq.5) then
               if(nope.eq.6) then
                  do k=1,nnodelem
                     nodef(k)=kon(indexe+ifacew1(k,jface))
                  enddo
               elseif(nope.eq.15) then
                  do k=1,nnodelem
                     nodef(k)=kon(indexe+ifacew2(k,jface))
                  enddo
               endif
            elseif(nface.eq.6) then
               do k=1,nnodelem
                  nodef(k)=kon(indexe+ifaceq(k,jface))
               enddo
            endif
!
!           inserting the nodes in ifield
!
            do k=1,nnodelem
               node=nodef(k)
               call nident(ifield(istartfield(i)),node,ilength,id)
               id=istartfield(i)+id-1
               if(id.gt.istartfield(i)-1) then
                  if(ifield(id).eq.node) cycle
               endif
               do l=ifree,id+2,-1
                  ifield(l)=ifield(l-1)
               enddo
               ifield(id+1)=node
               ifree=ifree+1
               ilength=ilength+1
            enddo
         enddo
         iendfield(i)=ifree-1
      enddo
!
      nconf=ifree-1
!     
      return
      end

