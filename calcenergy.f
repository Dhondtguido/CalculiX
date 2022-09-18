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
      subroutine calcenergy(ipkon,lakon,kon,co,ener,mi,ne,
     &     thicke,ielmat,energy,ielprop,prop,nea,neb)
!
!     calculates the energy in a *DYNAMIC calculation
!
      implicit none
!
      character*8 lakon(*),lakonl
!
      integer ipkon(*),nelem,kon(*),mi(*),nope,indexe,i,j,k,
     &  konl(20),mint3d,jj,iflag,ne,ki,kl,ilayer,nlayer,kk,
     &  nopes,ielmat(mi(3),*),mint2d,null,ielprop(*),nea,neb
!
      real*8 ener(2,mi(1),*),enerinttot,enerint,co(3,*),prop(*),
     &  xl(3,20),xi,et,ze,xsj,shp(4,20),weight,enerkintot,enerkin,
     &  a,gs(8,4),dlayer(4),tlayer(4),thickness,
     &  thicke(mi(3),*),xlayer(mi(3),4),shp2(7,8),xs2(3,7),xsj2(3),
     &  xl2(3,8),energy(*),enerelctot,enervictot,enerelc,enervic
!
      include "gauss.f"
!
      iflag=2
      null=0
!     
      enerinttot=0.d0
      enerkintot=0.d0
      enerelctot=0.d0
      enervictot=0.d0
!
      do nelem=nea,neb
         if(ipkon(nelem).lt.0) cycle
!
         lakonl=lakon(nelem)
         indexe=ipkon(nelem)
!     
         if(lakonl(1:5).eq.'C3D8I') then
            nope=11
         elseif(lakonl(4:4).eq.'2') then
            nope=20
         elseif(lakonl(4:4).eq.'8') then
            nope=8
         elseif(lakonl(4:5).eq.'10') then
            nope=10
         elseif(lakonl(4:4).eq.'4') then
            nope=4
         elseif(lakonl(4:5).eq.'15') then
            nope=15
         elseif(lakonl(4:5).eq.'6') then
            nope=6
         else
            nope=0
         endif
!     
!     composite materials
!     
         if(lakonl(7:8).eq.'LC') then
!     
!        determining the number of layers
!
            nlayer=0
            do k=1,mi(3)
               if(ielmat(k,nelem).ne.0) then
                  nlayer=nlayer+1
               endif
            enddo
!
            if(lakonl(4:4).eq.'2') then
               mint2d=4
               nopes=8
!     
!     determining the layer thickness and global thickness
!     at the shell integration points
!     
               iflag=1
               indexe=ipkon(nelem)
               do kk=1,mint2d
                  xi=gauss3d2(1,kk)
                  et=gauss3d2(2,kk)
                  call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
                  tlayer(kk)=0.d0
                  do i=1,nlayer
                     thickness=0.d0
                     do j=1,nopes
                        thickness=thickness+thicke(i,indexe+j)*shp2(4,j)
                     enddo
                     tlayer(kk)=tlayer(kk)+thickness
                     xlayer(i,kk)=thickness
                  enddo
               enddo
               iflag=2
!     
               ilayer=0
               do i=1,4
                  dlayer(i)=0.d0
               enddo
            elseif(lakonl(4:5).eq.'15') then
               mint2d=3
               nopes=6
!     
!     determining the layer thickness and global thickness
!     at the shell integration points
!     
               iflag=1
               indexe=ipkon(nelem)
               do kk=1,mint2d
                  xi=gauss3d10(1,kk)
                  et=gauss3d10(2,kk)
                  call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
                  tlayer(kk)=0.d0
                  do i=1,nlayer
                     thickness=0.d0
                     do j=1,nopes
                        thickness=thickness+thicke(i,indexe+j)*shp2(4,j)
                     enddo
                     tlayer(kk)=tlayer(kk)+thickness
                     xlayer(i,kk)=thickness
                  enddo
               enddo
               iflag=2
!     
               ilayer=0
               do i=1,3
                  dlayer(i)=0.d0
               enddo
            endif
!     
         endif
!     
         do j=1,nope
            konl(j)=kon(indexe+j)
            do k=1,3
               xl(k,j)=co(k,konl(j))
            enddo
         enddo
!     
         enerint=0.d0
         enerkin=0.d0
         enerelc=0.d0
         enervic=0.d0
!
         if(lakonl(4:5).eq.'8R') then
            mint3d=1
         elseif(lakonl(4:7).eq.'20RB') then
            if((lakonl(8:8).eq.'R').or.(lakonl(8:8).eq.'C')) then
               mint3d=50
            else
               call beamintscheme(lakonl,mint3d,ielprop(nelem),prop,
     &              null,xi,et,ze,weight)
            endif
         elseif((lakonl(4:4).eq.'8').or.
     &           (lakonl(4:6).eq.'20R')) then
            if(lakonl(7:8).eq.'LC') then
               mint3d=8*nlayer
            else
               mint3d=8
            endif
         elseif(lakonl(4:4).eq.'2') then
            mint3d=27
         elseif(lakonl(4:5).eq.'10') then
            mint3d=4
         elseif(lakonl(4:4).eq.'4') then
            mint3d=1
         elseif(lakonl(4:5).eq.'15') then
            if(lakonl(7:8).eq.'LC') then
               mint3d=6*nlayer
            else
               mint3d=9
            endif
         elseif(lakonl(4:5).eq.'6') then
            mint3d=2
         elseif(lakonl(1:2).eq.'ES') then
            if(lakonl(7:7).eq.'C') then
               enerelc=ener(1,1,nelem)
               enervic=ener(2,1,nelem)
            else
               enerint=ener(1,1,nelem)
            endif
            mint3d=0
         else
            cycle
         endif
!     
         do jj=1,mint3d
            if(lakonl(4:5).eq.'8R') then
               xi=gauss3d1(1,jj)
               et=gauss3d1(2,jj)
               ze=gauss3d1(3,jj)
               weight=weight3d1(jj)
            elseif(lakonl(4:7).eq.'20RB') then
               if((lakonl(8:8).eq.'R').or.(lakonl(8:8).eq.'C')) then
                  xi=gauss3d13(1,jj)
                  et=gauss3d13(2,jj)
                  ze=gauss3d13(3,jj)
                  weight=weight3d13(jj)
               else
                  call beamintscheme(lakonl,mint3d,ielprop(nelem),prop,
     &                 kk,xi,et,ze,weight)
               endif
            elseif((lakonl(4:4).eq.'8').or.
     &              (lakonl(4:6).eq.'20R'))
     &              then
               if(lakonl(7:8).ne.'LC') then
                  xi=gauss3d2(1,jj)
                  et=gauss3d2(2,jj)
                  ze=gauss3d2(3,jj)
                  weight=weight3d2(jj)
               else
                  kl=mod(jj,8)
                  if(kl.eq.0) kl=8
!     
                  xi=gauss3d2(1,kl)
                  et=gauss3d2(2,kl)
                  ze=gauss3d2(3,kl)
                  weight=weight3d2(kl)
!     
                  ki=mod(jj,4)
                  if(ki.eq.0) ki=4
!     
                  if(kl.eq.1) then
                     ilayer=ilayer+1
                     if(ilayer.gt.1) then
                        do i=1,4
                           dlayer(i)=dlayer(i)+xlayer(ilayer-1,i)
                        enddo
                     endif
                  endif
                  ze=2.d0*(dlayer(ki)+(ze+1.d0)/2.d0*xlayer(ilayer,ki))/
     &                 tlayer(ki)-1.d0
                  weight=weight*xlayer(ilayer,ki)/tlayer(ki)
               endif
            elseif(lakonl(4:4).eq.'2') then
               xi=gauss3d3(1,jj)
               et=gauss3d3(2,jj)
               ze=gauss3d3(3,jj)
               weight=weight3d3(jj)
            elseif(lakonl(4:5).eq.'10') then
               xi=gauss3d5(1,jj)
               et=gauss3d5(2,jj)
               ze=gauss3d5(3,jj)
               weight=weight3d5(jj)
            elseif(lakonl(4:4).eq.'4') then
               xi=gauss3d4(1,jj)
               et=gauss3d4(2,jj)
               ze=gauss3d4(3,jj)
               weight=weight3d4(jj)
            elseif(lakonl(4:5).eq.'15') then
               if(lakonl(7:8).ne.'LC') then
                  xi=gauss3d8(1,jj)
                  et=gauss3d8(2,jj)
                  ze=gauss3d8(3,jj)
                  weight=weight3d8(jj)
               else
                  kl=mod(jj,6)
                  if(kl.eq.0) kl=6
!     
                  xi=gauss3d10(1,kl)
                  et=gauss3d10(2,kl)
                  ze=gauss3d10(3,kl)
                  weight=weight3d10(kl)
!     
                  ki=mod(jj,3)
                  if(ki.eq.0) ki=3
!     
                  if(kl.eq.1) then
                     ilayer=ilayer+1
                     if(ilayer.gt.1) then
                        do i=1,3
                           dlayer(i)=dlayer(i)+xlayer(ilayer-1,i)
                        enddo
                     endif
                  endif
                  ze=2.d0*(dlayer(ki)+(ze+1.d0)/2.d0*xlayer(ilayer,ki))/
     &                 tlayer(ki)-1.d0
                  weight=weight*xlayer(ilayer,ki)/tlayer(ki)
               endif
            else
               xi=gauss3d7(1,jj)
               et=gauss3d7(2,jj)
               ze=gauss3d7(3,jj)
               weight=weight3d7(jj)
            endif
!
            if(lakonl(1:5).eq.'C3D8R') then
               call shape8hr(xl,xsj,shp,gs,a)
            elseif(lakonl(1:5).eq.'C3D8I') then
               call shape8hu(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.20) then
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
            enerint=enerint+weight*xsj*ener(1,jj,nelem)
            enerkin=enerkin+weight*xsj*ener(2,jj,nelem)
         enddo
!     
         enerinttot=enerinttot+enerint
         enerkintot=enerkintot+enerkin
         enerelctot=enerelctot+enerelc
         enervictot=enervictot+enervic
!
      enddo
!
      energy(1)=enerinttot
      energy(2)=enerkintot
      energy(3)=enerelctot
      energy(4)=enervictot
!
      return
      end
      
