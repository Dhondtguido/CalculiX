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
      subroutine printoutintfluidfem(prlab,ipkon,lakon,stx,
     &  mi,ii,nelem,qfx,orab,ielorien,norien,co,kon,
     &  ielmat,thicke,ielprop,prop,nelel,ithermal,orname,nelnew)
!
!     stores integration point results for element "nelem" in the .dat file
!
!     nelem is the element number from the input deck
!     nelel is a local number, created e.g. by renumbering; thus
!     far only applicable for FVM-CFD; for all other applications
!     nelem=nelel
!
      implicit none
!
      character*6 prlab(*)
      character*8 lakon(*)
      character*80 orname(*)
!
      integer ipkon(*),mi(*),nelem,l,ii,mint3d,j,k,nope,
     &  ielorien(mi(3),*),norien,kon(*),konl,indexe,m,iorien,iflag,
     &  ielmat(mi(3),*),nopes,mint2d,kk,ki,kl,nlayer,ilayer,
     &  null,ielprop(*),nelel,ithermal(*),nelef,nelnew(*)
!
      real*8 stx(6,mi(1),*),
     &  qfx(3,mi(1),*),xi,et,ze,xl(3,20),xsj,shp(4,20),
     &  coords(3,27),weight,orab(7,*),co(3,*),a(3,3),b(3,3),c(3,3),
     &  qfxl(3),thicke(mi(3),*),xsj2(3),shp2(7,8),xl2(3,8),xs2(3,7),
     &  thickness,tlayer(4),dlayer(4),xlayer(mi(3),4),
     &  prop(*)
!
      include "gauss.f"
!
      data iflag /1/
      data null /0/
!
      write(*,*) 'printoutintfluidfem ',nelel
      write(*,*) 'printoutintfluidfem new ',nelnew(nelel)
      nelef=nelnew(nelel)
!
      if(ipkon(nelef).lt.0) then
         return
      else
         indexe=ipkon(nelef)
      endif
!
!     check whether transformation is necessary (if orientation
!     is applied and output in local system is requested)
!
      if(lakon(nelef)(7:8).ne.'LC') then
         if((norien.eq.0).or.(prlab(ii)(6:6).eq.'G')) then
            iorien=0
         else
            iorien=max(0,ielorien(1,nelef))
         endif
      elseif(lakon(nelef)(4:5).eq.'20') then
!     
!     composite materials
!     
         if((norien.eq.0).or.(prlab(ii)(6:6).eq.'G')) then
            iorien=0
         else
!
!           check whether at least one layer has a transformation
!
            iorien=0
            do k=1,mi(3)
               if(ielorien(k,nelef).ne.0) then
                  iorien=max(0,ielorien(k,nelef))
                  exit
               endif
            enddo
         endif
!
!     determining the number of layers
!     
         mint2d=4
         nopes=8
!
         nlayer=0
         do k=1,mi(3)
            if(ielmat(k,nelef).ne.0) then
               nlayer=nlayer+1
            endif
         enddo
!     
!     determining the layer thickness and global thickness
!     at the shell integration points
!     
         iflag=1
         do kk=1,mint2d
            xi=gauss3d2(1,kk)
            et=gauss3d2(2,kk)
            call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
            tlayer(kk)=0.d0
            do k=1,nlayer
               thickness=0.d0
               do j=1,nopes
                  thickness=thickness+thicke(k,indexe+j)*shp2(4,j)
               enddo
               tlayer(kk)=tlayer(kk)+thickness
               xlayer(k,kk)=thickness
            enddo
         enddo
         iflag=3
!     
         ilayer=0
         do k=1,4
            dlayer(k)=0.d0
         enddo
!
!
      elseif(lakon(nelef)(4:5).eq.'15') then
!     
!     composite materials
!     
         if((norien.eq.0).or.(prlab(ii)(6:6).eq.'G')) then
            iorien=0
         else
!
!           check whether at least one layer has a transformation
!
            iorien=0
            do k=1,mi(3)
               if(ielorien(k,nelef).ne.0) then
                  iorien=max(0,ielorien(k,nelef))
                  exit
               endif
            enddo
         endif
!
!     determining the number of layers
!     
         mint2d=3
         nopes=6
!
         nlayer=0
         do k=1,mi(3)
            if(ielmat(k,nelef).ne.0) then
               nlayer=nlayer+1
            endif
         enddo
!     
!     determining the layer thickness and global thickness
!     at the shell integration points
!     
         iflag=1
         do kk=1,mint2d
            xi=gauss3d10(1,kk)
            et=gauss3d10(2,kk)
            call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
            tlayer(kk)=0.d0
            do k=1,nlayer
               thickness=0.d0
               do j=1,nopes
                  thickness=thickness+thicke(k,indexe+j)*shp2(4,j)
               enddo
               tlayer(kk)=tlayer(kk)+thickness
               xlayer(k,kk)=thickness
            enddo
         enddo
         iflag=3
!     
         ilayer=0
         do k=1,3
            dlayer(k)=0.d0
         enddo
!
      endif
!
!     number of integration points
!
      if(lakon(nelef)(4:5).eq.'8R') then
         mint3d=1
      elseif(lakon(nelef)(4:7).eq.'20RB') then
         if((lakon(nelef)(8:8).eq.'R').or.
     &      (lakon(nelef)(8:8).eq.'C')) then
            mint3d=50
         else
            call beamintscheme(lakon(nelef),mint3d,ielprop(nelel),prop,
     &           null,xi,et,ze,weight)
         endif
      elseif((lakon(nelef)(4:4).eq.'8').or.
     &       (lakon(nelef)(4:6).eq.'20R')) then
         if(lakon(nelef)(7:8).eq.'LC') then
            mint3d=8*nlayer
         else
            mint3d=8
         endif
      elseif(lakon(nelef)(4:4).eq.'2') then
         mint3d=27
      elseif(lakon(nelef)(4:5).eq.'10') then
         mint3d=4
      elseif(lakon(nelef)(4:4).eq.'4') then
         mint3d=1
      elseif(lakon(nelef)(4:5).eq.'15') then
         if(lakon(nelef)(7:8).eq.'LC') then
            mint3d=6*nlayer
         else
            mint3d=9
         endif
      elseif(lakon(nelef)(4:4).eq.'6') then
         mint3d=2
      elseif(lakon(nelef)(1:1).eq.'U') then
         mint3d=ichar(lakon(nelef)(6:6))
      else
         return
      endif
!
!     calculation of the integration point coordinates for
!     output in the local system (if needed)
!
      if((iorien.ne.0).or.(prlab(ii)(1:4).eq.'COOR')) then
         if(lakon(nelef)(4:4).eq.'2') then
            nope=20
         elseif(lakon(nelef)(4:4).eq.'8') then
            nope=8
         elseif(lakon(nelef)(4:5).eq.'10') then
            nope=10
         elseif(lakon(nelef)(4:4).eq.'4') then
            nope=4
         elseif(lakon(nelef)(4:5).eq.'15') then
            nope=15
         elseif(lakon(nelef)(4:4).eq.'6') then
            nope=6
         endif
!
         do j=1,nope
            konl=kon(indexe+j)
            do k=1,3
               xl(k,j)=co(k,konl)
            enddo
         enddo
!
         do j=1,mint3d
            if((lakon(nelef)(4:5).eq.'8R').or.
     &         (lakon(nelef)(1:4).eq.'F3D8')) then
               xi=gauss3d1(1,j)
               et=gauss3d1(2,j)
               ze=gauss3d1(3,j)
               weight=weight3d1(j)
            elseif(lakon(nelef)(4:8).eq.'20RB') then
               if((lakon(nelef)(8:8).eq.'B').or.
     &            (lakon(nelef)(8:8).eq.'C')) then
                  xi=gauss3d13(1,j)
                  et=gauss3d13(2,j)
                  ze=gauss3d13(3,j)
                  weight=weight3d13(j)
               else
                  call beamintscheme(lakon(nelef),mint3d,ielprop(nelel),
     &                  prop,j,xi,et,ze,weight)
               endif
            elseif((lakon(nelef)(4:4).eq.'8').or.
     &             (lakon(nelef)(4:6).eq.'20R'))
     &              then
               if(lakon(nelef)(7:8).ne.'LC') then
                  xi=gauss3d2(1,j)
                  et=gauss3d2(2,j)
                  ze=gauss3d2(3,j)
                  weight=weight3d2(j)
               else
                  kl=mod(j,8)
                  if(kl.eq.0) kl=8
!
                  xi=gauss3d2(1,kl)
                  et=gauss3d2(2,kl)
                  ze=gauss3d2(3,kl)
                  weight=weight3d2(kl)
!
                  ki=mod(j,4)
                  if(ki.eq.0) ki=4
!
                  if(kl.eq.1) then
                     ilayer=ilayer+1
                     if(ilayer.gt.1) then
                        do k=1,4
                           dlayer(k)=dlayer(k)+xlayer(ilayer-1,k)
                        enddo
                     endif
                  endif
                  ze=2.d0*(dlayer(ki)+(ze+1.d0)/2.d0*xlayer(ilayer,ki))/
     &                 tlayer(ki)-1.d0
                  weight=weight*xlayer(ilayer,ki)/tlayer(ki)
               endif
            elseif(lakon(nelef)(4:4).eq.'2') then
               xi=gauss3d3(1,j)
               et=gauss3d3(2,j)
               ze=gauss3d3(3,j)
               weight=weight3d3(j)
            elseif(lakon(nelef)(4:5).eq.'10') then
               xi=gauss3d5(1,j)
               et=gauss3d5(2,j)
               ze=gauss3d5(3,j)
               weight=weight3d5(j)
            elseif(lakon(nelef)(4:4).eq.'4') then
               xi=gauss3d4(1,j)
               et=gauss3d4(2,j)
               ze=gauss3d4(3,j)
               weight=weight3d4(j)
            elseif(lakon(nelef)(4:5).eq.'15') then
               if(lakon(nelef)(7:8).ne.'LC') then
                  xi=gauss3d8(1,j)
                  et=gauss3d8(2,j)
                  ze=gauss3d8(3,j)
                  weight=weight3d8(j)
               else
                  kl=mod(j,6)
                  if(kl.eq.0) kl=6
!
                  xi=gauss3d10(1,kl)
                  et=gauss3d10(2,kl)
                  ze=gauss3d10(3,kl)
                  weight=weight3d10(kl)
!
                  ki=mod(j,3)
                  if(ki.eq.0) ki=3
!
                  if(kl.eq.1) then
                     ilayer=ilayer+1
                     if(ilayer.gt.1) then
                        do k=1,3
                           dlayer(k)=dlayer(k)+xlayer(ilayer-1,k)
                        enddo
                     endif
                  endif
                  ze=2.d0*(dlayer(ki)+(ze+1.d0)/2.d0*xlayer(ilayer,ki))/
     &                 tlayer(ki)-1.d0
                  weight=weight*xlayer(ilayer,ki)/tlayer(ki)
               endif
            elseif(lakon(nelef)(1:4).eq.'C3D6') then
               xi=gauss3d7(1,j)
               et=gauss3d7(2,j)
               ze=gauss3d7(3,j)
               weight=weight3d7(j)
            elseif(lakon(nelef)(1:4).eq.'F3D6') then
               xi=gauss3d14(1,j)
               et=gauss3d14(2,j)
               ze=gauss3d14(3,j)
               weight=weight3d14(j)
            endif
!
            if(nope.eq.20) then
               call shape20h(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.8) then
               call shape8h(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.10) then
               call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.4) then
               call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.15) then
               call shape15w(xi,et,ze,xl,xsj,shp,iflag)
            else
               call shape6w(xi,et,ze,xl,xsj,shp,iflag)
            endif
!
            do k=1,3
               coords(k,j)=0.d0
               do l=1,nope
                  coords(k,j)=coords(k,j)+xl(k,l)*shp(4,l)
               enddo
            enddo
         enddo
      endif
!
      if((prlab(ii)(1:4).eq.'S   ').or.(prlab(ii)(1:4).eq.'SVF ')) then
         do j=1,mint3d
!
!           composite materials
!
            if(lakon(nelef)(7:8).eq.'LC') then
               if((norien.eq.0).or.(prlab(ii)(6:6).eq.'G')) then
                  iorien=0
               elseif(lakon(nelef)(4:5).eq.'20') then
                  kl=mod(j,8)
                  if(kl.eq.0) kl=8
                  ilayer=(j-kl)/8+1
                  iorien=max(0,ielorien(ilayer,nelef))
               elseif(lakon(nelef)(4:5).eq.'15') then
                  kl=mod(j,6)
                  if(kl.eq.0) kl=6
                  ilayer=(j-kl)/6+1
                  iorien=max(0,ielorien(ilayer,nelef))
               endif
            endif
!
            if(iorien.eq.0) then
               write(5,'(i10,1x,i3,1p,6(1x,e13.6))') nelem,j,
     &              (stx(k,j,nelef),k=1,6)
            else
               call transformatrix(orab(1,iorien),coords(1,j),a)
               b(1,1)=stx(1,j,nelef)
               b(2,2)=stx(2,j,nelef)
               b(3,3)=stx(3,j,nelef)
               b(1,2)=stx(4,j,nelef)
               b(1,3)=stx(5,j,nelef)
               b(2,3)=stx(6,j,nelef)
               b(2,1)=b(1,2)
               b(3,1)=b(1,3)
               b(3,2)=b(2,3)
               do k=1,3
                  do l=1,3
                     c(k,l)=0.d0
                     do m=1,3
                        c(k,l)=c(k,l)+b(k,m)*a(m,l)
                     enddo
                  enddo
               enddo
               do k=1,3
                  do l=k,3
                     b(k,l)=0.d0
                     do m=1,3
                        b(k,l)=b(k,l)+a(m,k)*c(m,l)
                     enddo
                  enddo
               enddo
               write(5,'(i10,1x,i3,1p,6(1x,e13.6),1x,a20)') nelem,j,
     &              b(1,1),b(2,2),b(3,3),b(1,2),b(1,3),b(2,3),
     &              orname(iorien)(1:20)
            endif
         enddo
      elseif(((prlab(ii)(1:4).eq.'HFL ').or.(prlab(ii)(1:4).eq.'HFLF'))
     &        .and.(ithermal(1).gt.1)) then
         do j=1,mint3d
!
!           composite materials
!
            if(lakon(nelef)(7:8).eq.'LC') then
               if((norien.eq.0).or.(prlab(ii)(6:6).eq.'G')) then
                  iorien=0
               elseif(lakon(nelef)(4:5).eq.'20') then
                  kl=mod(j,8)
                  if(kl.eq.0) kl=8
                  ilayer=(j-kl)/8+1
                  iorien=max(0,ielorien(ilayer,nelef))
               elseif(lakon(nelef)(4:5).eq.'15') then
                  kl=mod(j,6)
                  if(kl.eq.0) kl=6
                  ilayer=(j-kl)/6+1
                  iorien=max(0,ielorien(ilayer,nelef))
               endif
            endif
!
            if(iorien.eq.0) then
               write(5,'(i10,1x,i3,1p,3(1x,e13.6))') nelem,j,
     &              (qfx(k,j,nelef),k=1,3)
            else
               do k=1,3
                  qfxl(k)=qfx(k,j,nelef)
               enddo
               call transformatrix(orab(1,iorien),coords(1,j),a)
               write(5,'(i10,1x,i3,1p,3(1x,e13.6),1x,a20)') nelem,j,
     &              qfxl(1)*a(1,1)+qfxl(2)*a(2,1)+qfxl(3)*a(3,1),
     &              qfxl(1)*a(1,2)+qfxl(2)*a(2,2)+qfxl(3)*a(3,2),
     &              qfxl(1)*a(1,3)+qfxl(2)*a(2,3)+qfxl(3)*a(3,3),
     &              orname(iorien)(1:20)
            endif
         enddo
      elseif(prlab(ii)(1:4).eq.'COOR') then
         do j=1,mint3d
            write(5,'(i10,1x,i3,1p,3(1x,e13.6))') nelem,j,
     &           (coords(k,j),k=1,3)
         enddo
      endif
!
      return
      end
