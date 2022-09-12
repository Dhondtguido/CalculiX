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
      subroutine triangulate(ics,rcs0,zcs0,ncsnodes,
     &  rcscg,rcs0cg,zcscg,zcs0cg,nrcg,nzcg,jcs,kontri,straight,
     &  ne,ipkon,kon,lakon,lcs,netri,ifacetet,inodface)
!
!     generate a triangulation of the independent side (= right side)
!
!     the element faces of the independent side are identified and
!     triangulated. The nodes belonging to the faces are stored in
!     field inodface, face after face. For a triangle i the value
!     ifacetet(i) points to the last node in field inodface of the
!     face the triangle belongs to.
!
      implicit none
!
      character*8 lakon(*)
!
      integer jcs(*),l,j,ics(*),nodef(8),ifacetet(*),
     &  nrcg(*),node,ncsnodes,id,ifaceq(8,6),ifacet(6,4),
     &  ifacew1(4,5),iface(8,6),nodelem(20),nnodelem,nzcg(*),
     &  itrifac3(3,1),itrifac4(3,2),itrifac6(3,4),itrifac8(3,6),
     &  itrifac(3,6),ifacew2(8,5),lcs(*),inodface(*),nnodface,
     &  k,kflag,i,ne,ipkon(*),kon(*),indexe,nope,nface,nodface,jface,
     &  netri,ntrifac,kontri(3,*),m
!
      real*8 straight(9,*),zcscg(*),rcscg(*),zcs0cg(*),
     &  rcs0cg(*),cgl(2),col(2,3),rcs0(*),zcs0(*)
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
!     pointer into field inodface
!
      nnodface=0
!
!     sort the nodes on the independent side
!
      do j=1,ncsnodes
         jcs(j)=abs(ics(j))
         lcs(j)=j
      enddo
!
      kflag=2
      call isortii(jcs,lcs,ncsnodes,kflag)
!
      netri=0
!
!     check the elements adjacent to the independent nodes
!
      do i=1,ne
         indexe=ipkon(i)
         if(lakon(i)(4:4).eq.'2') then
            nope=20
            nface=6
            nodface=8
         elseif(lakon(i)(4:4).eq.'8') then
            nope=8
            nface=6
            nodface=4
         elseif(lakon(i)(4:5).eq.'10') then
            nope=10
            nface=4
            nodface=6
         elseif(lakon(i)(4:4).eq.'4') then
            nope=4
            nface=4
            nodface=3
         elseif(lakon(i)(4:5).eq.'15') then
            nope=15
            nface=5
            nodface=8
         elseif(lakon(i)(4:4).eq.'6') then
            nope=6
            nface=5
            nodface=4
         else
            cycle
         endif
!
!        check which nodes of the element belong to the independent set
!
         nnodelem=0
         do j=1,nope
            nodelem(j)=0
            node=kon(indexe+j)
            call nident(jcs,node,ncsnodes,id)
            if(id.le.0) cycle
            if(jcs(id).ne.node) cycle
            nodelem(j)=node
            nnodelem=nnodelem+1
         enddo
         if(nnodelem.eq.0) cycle
!
         if(nface.eq.4) then
            do j=1,nface
               do k=1,nodface
                  iface(k,j)=ifacet(k,j)
               enddo
            enddo
         elseif(nface.eq.5) then
            if(nope.eq.6) then
               do j=1,nface
                  do k=1,nodface
                     iface(k,j)=ifacew1(k,j)
                  enddo
               enddo
            elseif(nope.eq.15) then
               do j=1,nface
                  do k=1,nodface
                     iface(k,j)=ifacew2(k,j)
                  enddo
               enddo
            endif
         elseif(nface.eq.6) then
            do j=1,nface
               do k=1,nodface
                  iface(k,j)=ifaceq(k,j)
               enddo
            enddo
         endif
!
!        check which face of the element belongs to the independent side 
!
         jface=0
         loop: do m=1,nface
            nnodelem=nodface
!
!           several faces of one and the same element may belong
!           to the master surface
!
            do k=1,nodface
               if(iface(k,m).eq.0) then
                  nnodelem=k-1
                  exit
               endif
               if(nodelem(iface(k,m)).eq.0) cycle loop
            enddo
            jface=m
c            exit
c         enddo loop
c         if(jface.eq.0) cycle
!
!        store the node numbers in a local face field
!
            do k=1,nnodelem
               nodef(k)=nodelem(iface(k,jface))
               inodface(nnodface+k)=nodef(k)
            enddo
            nnodface=nnodface+nnodelem
!
!        number of triangles
!
            if(nnodelem.eq.3) then
               ntrifac=1
               do j=1,ntrifac
                  do k=1,3
                     itrifac(k,j)=itrifac3(k,j)
                  enddo
               enddo
            elseif(nnodelem.eq.4) then
               ntrifac=2
               do j=1,ntrifac
                  do k=1,3
                     itrifac(k,j)=itrifac4(k,j)
                  enddo
               enddo
            elseif(nnodelem.eq.6) then
               ntrifac=4
               do j=1,ntrifac
                  do k=1,3
                     itrifac(k,j)=itrifac6(k,j)
                  enddo
               enddo
            elseif(nnodelem.eq.8) then
               ntrifac=6
               do j=1,ntrifac
                  do k=1,3
                     itrifac(k,j)=itrifac8(k,j)
                  enddo
               enddo
            endif
!
            do j=1,ntrifac
!
!           new triangle
!
               netri=netri+1
               do l=1,2
                  cgl(l)=0.d0
               enddo
               do k=1,3
                  node=nodef(itrifac(k,j))
                  kontri(k,netri)=node
                  call nident(jcs,node,ncsnodes,id)
                  col(1,k)=rcs0(lcs(id))
                  col(2,k)=zcs0(lcs(id))
                  do l=1,2
                     cgl(l)=cgl(l)+col(l,k)
                  enddo
               enddo
!
!           center of gravity of the triangle
!
               rcscg(netri)=cgl(1)/3.d0
               zcscg(netri)=cgl(2)/3.d0
!
!           determining the equations of the straight lines bordering
!           the triangle
!            
               call straighteq2d(col,straight(1,netri))
!
               ifacetet(netri)=nnodface
!
            enddo
         enddo loop
      enddo
!
      if(netri.eq.0) then
         write(*,*) '*ERROR in triangulate: no faces found on the'
         write(*,*) '       independent side'
         call exit(201)
      endif
!
!     initialization of near2d
!
      do i=1,netri
         nrcg(i)=i
         nzcg(i)=i
         rcs0cg(i)=rcscg(i)
         zcs0cg(i)=zcscg(i)
      enddo
!
      kflag=2
      call dsort(rcscg,nrcg,netri,kflag)
      call dsort(zcscg,nzcg,netri,kflag)
!
      return
      end

