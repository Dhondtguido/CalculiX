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
      subroutine calcexternalwork(co,vold,istartset,iendset,ipkon,
     &  lakon,kon,ialset,nset,set,mi,externalfaces,nelemload,
     &  nload,xload,sideload,delexternalwork,voldprev)
!
!     calculation and printout of the lift and drag forces
!
      implicit none
!
      character*8 lakonl,lakon(*)
      character*20 sideload(*)
      character*81 set(*),externalfaces
!
      integer konl(20),ifaceq(8,6),nelem,i,j,i1,j1,
     &  jj,ig,nope,nopes,igpres,id,
     &  mint2d,ifacet(6,4),ifacew(8,5),iflag,indexe,jface,istartset(*),
     &  iendset(*),ipkon(*),kon(*),iset,ialset(*),nset,ipos,
     &  mi(*),nelemload(2,*),nload,ipres
!
      real*8 co(3,*),xs2(3,7),pres,voldl2(3,8),xl2(3,8),xsj2(3),
     &  shp2(7,8),vold(0:mi(2),*),xi,et,weight,voldprev(0:mi(2),*),
     &  xlocal20(3,9,6),xlocal4(3,1,4),xlocal10(3,3,4),xlocal6(3,1,5),
     &  xlocal15(3,4,5),xlocal8(3,4,6),xlocal8r(3,1,6),
     &  disp(3),xload(2,*),delexternalwork
!
      include "gauss.f"
      include "xlocal.f"
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
      data iflag /3/
!     
!           determining the set with the external faces
!
      ipos=index(externalfaces,' ')
      externalfaces(ipos:ipos)='T'
      do iset=1,nset
         if(set(iset)(1:ipos).eq.externalfaces(1:ipos)) exit
      enddo
      externalfaces(ipos:ipos)=' '
!
      delexternalwork=0.d0
!     
      do jj=istartset(iset),iendset(iset)
!     
         jface=ialset(jj)
!     
         nelem=int(jface/10.d0)
c** start change 20191219
         indexe=ipkon(nelem)
         if(indexe.lt.0) cycle
c** end chenge 20191219
         ig=jface-10*nelem
!     
!        determining the distributed load on the face
!     
         ipres=0
         call nident2(nelemload,nelem,nload,id)
         do
            if((id.eq.0).or.(nelemload(1,id).ne.nelem)) exit
            if(sideload(id)(1:1).ne.'P') then
               id=id-1
               cycle
            endif
            igpres=ichar(sideload(id)(2:2))-48
            if(ig.ne.igpres) then
               id=id-1
               cycle
            else
               ipres=1
               pres=xload(1,id)
               exit
            endif
         enddo
c         if(ipres.eq.0) then
c            write(*,*) '*ERROR in calcexternalwork: external face',
c     &           ig,' of element ',nelem,' has no external load'
c            call exit(201)
c         endif
!     
c         indexe=ipkon(nelem)
         lakonl=lakon(nelem)
!     
!        determining the number of nodes in the face
!
         if(lakonl(4:4).eq.'2') then
            nope=20
            nopes=8
         elseif(lakonl(4:4).eq.'8') then
            nope=8
            nopes=4
         elseif(lakonl(4:5).eq.'10') then
            nope=10
            nopes=6
         elseif(lakonl(4:4).eq.'4') then
            nope=4
            nopes=3
         elseif(lakonl(4:5).eq.'15') then
            nope=15
         elseif(lakonl(4:4).eq.'6') then
            nope=6
         endif
!
!        determining the number of integration points in the face
!     
         if(lakonl(4:5).eq.'8R') then
            mint2d=1
         elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) 
     &           then
            if((lakonl(7:7).eq.'A').or.(lakonl(7:7).eq.'S').or.
     &           (lakonl(7:7).eq.'E')) then
               mint2d=2
            else
               mint2d=4
            endif
         elseif(lakonl(4:4).eq.'2') then
            mint2d=9
         elseif(lakonl(4:5).eq.'10') then
            mint2d=3
         elseif(lakonl(4:4).eq.'4') then
            mint2d=1
         endif
!     
!        local topology
!     
         do i=1,nope
            konl(i)=kon(indexe+i)
         enddo
!     
!        treatment of wedge faces
!     
         if(lakonl(4:4).eq.'6') then
            mint2d=1
            if(ig.le.2) then
               nopes=3
            else
               nopes=4
            endif
         endif
         if(lakonl(4:5).eq.'15') then
            if(ig.le.2) then
               mint2d=3
               nopes=6
            else
               mint2d=4
               nopes=8
            endif
         endif

!
!        coordinates and differential displacements of the nodes belonging
!        to the face face
!
         if((nope.eq.20).or.(nope.eq.8)) then
            do i=1,nopes
               do j=1,3
                  voldl2(j,i)=vold(j,konl(ifaceq(i,ig)))
                  xl2(j,i)=co(j,konl(ifaceq(i,ig)))
     &                 +voldl2(j,i)
                  voldl2(j,i)=voldl2(j,i)-voldprev(j,konl(ifaceq(i,ig)))
               enddo
            enddo
         elseif((nope.eq.10).or.(nope.eq.4)) then
            do i=1,nopes
               do j=1,3
                  voldl2(j,i)=vold(j,konl(ifacet(i,ig)))
                  xl2(j,i)=co(j,konl(ifacet(i,ig)))
     &                 +voldl2(j,i)
                 voldl2(j,i)=voldl2(j,i)-voldprev(j,konl(ifacet(i,ig)))
               enddo
            enddo
         else
            do i=1,nopes
               do j=1,3
                  voldl2(j,i)=vold(j,konl(ifacew(i,ig)))
                  xl2(j,i)=co(j,konl(ifacew(i,ig)))
     &                 +voldl2(j,i)
                  voldl2(j,i)=voldl2(j,i)-voldprev(j,konl(ifacew(i,ig)))
               enddo
            enddo
         endif
!     
         do i=1,mint2d
!     
!          local coordinates of the surface integration
!           point within the surface local coordinate system
!     
            if((lakonl(4:5).eq.'8R').or.
     &           ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
               xi=gauss2d1(1,i)
               et=gauss2d1(2,i)
               weight=weight2d1(i)
            elseif((lakonl(4:4).eq.'8').or.
     &              (lakonl(4:6).eq.'20R').or.
     &              ((lakonl(4:5).eq.'15').and.(nopes.eq.8))) then
               xi=gauss2d2(1,i)
               et=gauss2d2(2,i)
               weight=weight2d2(i)
            elseif(lakonl(4:4).eq.'2') then
               xi=gauss2d3(1,i)
               et=gauss2d3(2,i)
               weight=weight2d3(i)
            elseif((lakonl(4:5).eq.'10').or.
     &              ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
               xi=gauss2d5(1,i)
               et=gauss2d5(2,i)
               weight=weight2d5(i)
            elseif((lakonl(4:4).eq.'4').or.
     &              ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
               xi=gauss2d4(1,i)
               et=gauss2d4(2,i)
               weight=weight2d4(i)
            endif
!     
!           local surface normal
!     
            if(nopes.eq.8) then
               call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
            elseif(nopes.eq.4) then
               call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
            elseif(nopes.eq.6) then
               call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
            else
               call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
            endif
!
!           differential displacements at the integration point
!
            do j1=1,3
               disp(j1)=0.d0
               do i1=1,nopes
                  disp(j1)=disp(j1)+shp2(4,i1)*voldl2(j1,i1)
               enddo
            enddo
!
!           total external work (pressure is a negative load)
!
            delexternalwork=delexternalwork-pres*weight*
     &        (disp(1)*xsj2(1)+disp(2)*xsj2(2)+disp(3)*xsj2(3))
         enddo
      enddo
!
      return
      end
