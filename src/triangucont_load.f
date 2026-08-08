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
      subroutine triangucont_load(ncont,istartset,
     &     iendset,ialset,lakon,ipkon,kon,koncont,
     &     imastset)
!
!     generate a triangulation of the external faces of the volumetric
!     mesh
!
      implicit none
!
      character*8 lakon(*)
!
      integer ncont,j,k,l,istartset(*),iendset(*),ialset(*),
     &     nelem,jface,indexe,ipkon(*),nope,
     &     ifaceq(8,6),ifacet(6,4),ifacew1(4,5),ifacew2(8,5),node,
     &     ntrifac,itrifac3(3,1),itrifac4(3,2),itrifac6(3,4),
     &     itrifac8(3,6),itrifac(3,8),nnodelem,nface,nodef(9),kon(*),
     &     koncont(4,*),imastset
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
!     triangulation for three-node face
!
      data itrifac3 /1,2,3/
!
!     triangulation for four-node face
!
      data itrifac4 /1,2,4,2,3,4/
!
!     triangulation for six-node face
!
      data itrifac6 /1,4,6,4,2,5,6,5,3,4,5,6/
!
!     triangulation for eight-node face
!
      data itrifac8 /1,5,8,5,2,6,7,6,3,8,7,4,8,5,7,5,6,7/
!     
      ncont=0
!     
      do j=istartset(imastset),iendset(imastset)
!     
        nelem=int(ialset(j)/10.d0)
        indexe=ipkon(nelem)
        if(indexe.lt.0) cycle
!     
        jface=ialset(j)-10*nelem
!     
!     for dynamic massless contact: store contact face number
!     in element label
!     
        if(lakon(nelem)(4:5).eq.'20') then
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
!     number of triangles
!     
        if(nnodelem.eq.3) then
          ntrifac=1
          do l=1,ntrifac
            do k=1,3
              itrifac(k,l)=itrifac3(k,l)
            enddo
          enddo
        elseif(nnodelem.eq.4) then
          ntrifac=2
          do l=1,ntrifac
            do k=1,3
              itrifac(k,l)=itrifac4(k,l)
            enddo
          enddo
        elseif(nnodelem.eq.6) then
          ntrifac=4
          do l=1,ntrifac
            do k=1,3
              itrifac(k,l)=itrifac6(k,l)
            enddo
          enddo
        elseif(nnodelem.eq.8) then
          ntrifac=6
          do l=1,ntrifac
            do k=1,3
              itrifac(k,l)=itrifac8(k,l)
            enddo
          enddo
        endif
!     
!     storing the topology of the triangles
!     
        do l=1,ntrifac
!     
          ncont=ncont+1
          do k=1,3
            node=nodef(itrifac(k,l))
            koncont(k,ncont)=node
          enddo
          koncont(4,ncont)=ialset(j)
!     
        enddo
      enddo
!     
      return
      end
