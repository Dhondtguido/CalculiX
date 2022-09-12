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
      subroutine biotsavart(ipkon,kon,lakon,ne,co,qfx,h0,mi,nka,nkb)
!
      implicit none
!
!     calculates the magnetic intensity due to currents in the phi-
!     domain of an electromagnetic calculation
!
      character*8 lakon(*)
!
      integer ipkon(*),kon(*),ne,i,nka,nkb,mint3d,konl(26),
     &  j,k,indexe,kk,iflag,mi(*),nope,l
!
      real*8 co(3,*),qfx(3,mi(1),*),h0(3,*),xl(3,26),r(3),c2,
     &  con(3),pgauss(3),c1,xi,et,ze,xsj,shp(4,20),weight
!
      include "gauss.f"
!
      c1=1.d0/(16.d0*datan(1.d0))
      iflag=2
!
      do j=nka,nkb
!
         do k=1,3
            con(k)=co(k,j)
         enddo
!
         do i=1,ne
            if(ipkon(i).lt.0) cycle
!
!           currents are supposed to be modeled by shell elements
!           only
!
            if(lakon(i)(7:7).ne.'L') cycle
!
            if(lakon(i)(4:5).eq.'8R') then
               mint3d=1
               nope=8
            elseif(lakon(i)(4:4).eq.'8') then
               mint3d=8
               nope=8
            elseif(lakon(i)(4:6).eq.'20R') then
               mint3d=8
               nope=20
            elseif(lakon(i)(4:4).eq.'2') then
               mint3d=27
               nope=20
            elseif(lakon(i)(4:5).eq.'15') then
               mint3d=9
               nope=15
            elseif(lakon(i)(4:4).eq.'6') then
               mint3d=2
               nope=6
            endif
!
            indexe=ipkon(i)
!
            do l=1,nope
               konl(l)=kon(indexe+l)
               do k=1,3
                  xl(k,l)=co(k,konl(l))
               enddo
            enddo
!
            do kk=1,mint3d
!     
               if(lakon(i)(4:5).eq.'8R') then
                  xi=gauss3d1(1,kk)
                  et=gauss3d1(2,kk)
                  ze=gauss3d1(3,kk)
                  weight=weight3d1(kk)
               elseif((lakon(i)(4:4).eq.'8').or.
     &                 (lakon(i)(4:6).eq.'20R'))
     &                 then
                  xi=gauss3d2(1,kk)
                  et=gauss3d2(2,kk)
                  ze=gauss3d2(3,kk)
                  weight=weight3d2(kk)
               elseif(lakon(i)(4:4).eq.'2') then
                  xi=gauss3d3(1,kk)
                  et=gauss3d3(2,kk)
                  ze=gauss3d3(3,kk)
                  weight=weight3d3(kk)
               elseif(lakon(i)(4:5).eq.'15') then
                  xi=gauss3d8(1,kk)
                  et=gauss3d8(2,kk)
                  ze=gauss3d8(3,kk)
                  weight=weight3d8(kk)
               elseif(lakon(i)(4:4).eq.'6') then
                  xi=gauss3d7(1,kk)
                  et=gauss3d7(2,kk)
                  ze=gauss3d7(3,kk)
                  weight=weight3d7(kk)
               endif
!     
!              shape functions
!
               if(nope.eq.20) then
                  call shape20h(xi,et,ze,xl,xsj,shp,iflag)
               elseif(nope.eq.8) then
                  call shape8h(xi,et,ze,xl,xsj,shp,iflag)
               elseif(nope.eq.15) then
                  call shape15w(xi,et,ze,xl,xsj,shp,iflag)
               elseif(nope.eq.6) then
                  call shape6w(xi,et,ze,xl,xsj,shp,iflag)
               endif
!     
!              coordinates of the gauss point
!
               do k=1,3
                  pgauss(k)=0.d0
                  do l=1,nope
                     pgauss(k)=pgauss(k)+shp(4,l)*xl(k,l)
                  enddo
               enddo
!
!              distance from node to gauss point
!
               do k=1,3
                  r(k)=con(k)-pgauss(k)
               enddo

               c2=weight*xsj/((r(1)*r(1)+r(2)*r(2)+r(3)*r(3))**(1.5d0))
!
               h0(1,j)=h0(1,j)+c2*
     &                   (qfx(2,kk,i)*r(3)-qfx(3,kk,i)*r(2))
               h0(2,j)=h0(2,j)+c2*
     &                   (qfx(3,kk,i)*r(1)-qfx(1,kk,i)*r(3))
               h0(3,j)=h0(3,j)+c2*
     &                   (qfx(1,kk,i)*r(2)-qfx(2,kk,i)*r(1))
            enddo
         enddo
!
         do k=1,3
            h0(k,j)=h0(k,j)*c1
         enddo
      enddo
!     
      return
      end
