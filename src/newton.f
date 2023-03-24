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
      subroutine newton(icalccg,ne,ipkon,lakon,kon,t0,co,rhcon,
     &       nrhcon,ntmat_,physcon,nelem,cgr,bodyf,ielmat,ithermal,
     &       vold,mi)
!
!     assigns the body forces to the elements by use of field ipobody
!
      implicit none
!
      character*8 lakon(*)
!
      integer i,j,k,ne,icalccg,ipkon(*),nope,konl(20),kon(*),two,id,
     &  nrhcon(*),ntmat_,nelem,indexe,imat,mi(*),ielmat(mi(3),*),
     &  ithermal(*),iflag
!
      real*8 xi,et,ze,weight,xl(3,20),shp(4,20),xsj,rho,cgr(4,*),
     &  t0l,t0(*),rhcon(0:1,ntmat_,*),physcon(*),co(3,*),dd,bodyf(3),
     &  vold(0:mi(2),*)
!
!
!
c      data two /2/
      two=2
!
      if(icalccg.eq.0) then
!
!        first call: calculate the center of gravity of all elements
!
         icalccg=1
         do i=1,ne
            if(ipkon(i).lt.0) cycle
            if(lakon(i)(4:4).eq.'2') then
               nope=20
               xi=0.d0
               et=0.d0
               ze=0.d0
               weight=8.d0
            elseif(lakon(i)(4:4).eq.'8') then
               nope=8
               xi=0.d0
               et=0.d0
               ze=0.d0
               weight=8.d0
            elseif(lakon(i)(4:5).eq.'10') then
               nope=10
               xi=0.25d0
               et=0.25d0
               ze=0.25d0
               weight=1.d0/6.d0
            elseif(lakon(i)(4:4).eq.'4') then
               nope=4
               xi=0.25d0
               et=0.25d0
               ze=0.25d0
               weight=1.d0/6.d0
            elseif(lakon(i)(4:5).eq.'15') then
               nope=15
               xi=1.d0/3.d0
               et=1.d0/3.d0
               ze=0.d0
               weight=1.d0
            elseif(lakon(i)(4:4).eq.'6') then
               nope=6
               xi=1.d0/3.d0
               et=1.d0/3.d0
               ze=0.d0
               weight=1.d0
            else
               cycle
            endif
!
!           coordinates of the nodes in the deformed configuration
!
            indexe=ipkon(i)
            do j=1,nope
               konl(j)=kon(indexe+j)
               do k=1,3
                  xl(k,j)=co(k,konl(j))+vold(k,konl(j))
               enddo
            enddo
!
!           calculation of the shape functions
!           in the gauss point
!
            iflag=1
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
!           calculation of the center of gravity
!
            do k=1,3
               cgr(k,i)=0.d0
               do j=1,nope
                  cgr(k,i)=cgr(k,i)+shp(4,j)*xl(k,j)
               enddo
            enddo
!
!           determining the density
!
            imat=ielmat(1,nelem)
            if(ithermal(1).eq.0) then
               rho=rhcon(1,1,imat)
            else
!
!              calculation of the initial temperature
!
               t0l=0.d0
               do j=1,nope
                  t0l=t0l+shp(4,j)*t0(konl(j))
               enddo
               call ident2(rhcon(0,1,imat),t0l,nrhcon(imat),two,id)
               if(nrhcon(imat).eq.0) then
                  continue
               elseif(nrhcon(imat).eq.1) then
                  rho=rhcon(1,1,imat)
               elseif(id.eq.0) then
                  rho=rhcon(1,1,imat)
               elseif(id.eq.nrhcon(imat)) then
                  rho=rhcon(1,id,imat)
               else
                  rho=rhcon(1,id,imat)+
     &                 (rhcon(1,id+1,imat)-rhcon(1,id,imat))*
     &                 (t0l-rhcon(0,id,imat))/
     &                 (rhcon(0,id+1,imat)-rhcon(0,id,imat))
               endif
            endif
!
!           coordinates of the nodes in the undeformed configuration
!
            indexe=ipkon(i)
            do j=1,nope
               konl(j)=kon(indexe+j)
               do k=1,3
                  xl(k,j)=co(k,konl(j))+vold(k,konl(j))
               enddo
            enddo
!
!           calculation of the Jacobian determinant
!           in the gauss point
!
            iflag=2
            if(nope.eq.20) then
               if(lakon(i)(7:7).eq.'A') then
                  call shape20h_ax(xi,et,ze,xl,xsj,shp,iflag)
               elseif((lakon(i)(7:7).eq.'E').or.
     &                (lakon(i)(7:7).eq.'S')) then
                  call shape20h_pl(xi,et,ze,xl,xsj,shp,iflag)
               else
                  call shape20h(xi,et,ze,xl,xsj,shp,iflag)
               endif
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
!           calculating the "constant" term in the gravity force
!
            cgr(4,i)=physcon(3)*rho*xsj*weight
!
         enddo
      endif
!
!     calculating the force per unit mass
!
      do j=1,3
         bodyf(j)=0.d0
      enddo
!
      do i=1,ne
         if(ipkon(i).lt.0) cycle
         if(i.eq.nelem) cycle
         dd=(cgr(1,i)-cgr(1,nelem))**2+(cgr(2,i)-cgr(2,nelem))**2+
     &      (cgr(3,i)-cgr(3,nelem))**2
         do j=1,3
            bodyf(j)=bodyf(j)+(cgr(j,i)-cgr(j,nelem))*cgr(4,i)/
     &               (dd*dsqrt(dd))
         enddo
      enddo
c      write(*,*) 'newton',nelem,bodyf(1),bodyf(2),bodyf(3)
!
      return
      end

