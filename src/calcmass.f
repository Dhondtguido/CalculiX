!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine calcmass(ipkon,lakon,kon,co,
     &     mi,nelem,ne,thicke,ielmat,nope,t0,rhcon,nrhcon,ntmat_,
     &     ithermal,csmasstot,ielprop,prop)
!     
!     calculates the mass of element "nelem"
!     
      implicit none
!     
      character*8 lakon(*)
!     
      integer ipkon(*),nelem,kon(*),mi(*),nope,indexe,i,j,k,konl(20),
     &     mint3d,jj,iflag,ne,ki,kl,ilayer,nlayer,kk,ntmat_,i1,imat,
     &     nopes,ielmat(mi(3),*),mint2d,nrhcon(*),ithermal(*),null,
     &     ielprop(*)
!     
      real*8 csmasstot,csmass,co(3,*),prop(*),
     &     xl(3,20),xi,et,ze,xsj,shp(4,20),weight,
     &     a,gs(8,4),dlayer(4),tlayer(4),thickness,t0l,rho,
     &     thicke(mi(3),*),xlayer(mi(3),4),shp2(7,8),xs2(3,7),xsj2(3),
     &     xl2(3,8),t0(*),rhcon(0:1,ntmat_,*)
!     
      include "gauss.f"
!     
      data iflag /2/
      null=0
!     
      if(ipkon(nelem).lt.0) return
!     
      indexe=ipkon(nelem)
!     
      if(lakon(nelem)(1:5).eq.'C3D8I') then
        nope=11
      elseif(lakon(nelem)(4:4).eq.'2') then
        nope=20
      elseif(lakon(nelem)(4:4).eq.'8') then
        nope=8
      elseif(lakon(nelem)(4:5).eq.'10') then
        nope=10
      elseif(lakon(nelem)(4:4).eq.'4') then
        nope=4
      elseif(lakon(nelem)(4:5).eq.'15') then
        nope=15
      elseif(lakon(nelem)(4:5).eq.'6') then
        nope=6
      else
        nope=0
      endif
!     
!     composite materials
!     
      if(lakon(nelem)(7:8).ne.'LC') then
        imat=ielmat(1,nelem)
      else
!     
!     determining the number of layers
!     
        nlayer=0
        do k=1,mi(3)
          if(ielmat(k,nelem).ne.0) then
            nlayer=nlayer+1
          endif
        enddo
!     
        if(lakon(nelem)(4:4).eq.'2') then
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
        elseif(lakon(nelem)(4:5).eq.'15') then
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
      csmass=0.d0
!     
      if(lakon(nelem)(4:5).eq.'8R') then
        mint3d=1
      elseif(lakon(nelem)(4:7).eq.'20RB') then
        if((lakon(nelem)(8:8).eq.'R').or.
     &       (lakon(nelem)(8:8).eq.'C')) then
          mint3d=50
        else
          call beamintscheme(lakon(nelem),mint3d,ielprop(nelem),prop,
     &         null,xi,et,ze,weight)
        endif
      elseif((lakon(nelem)(4:4).eq.'8').or.
     &       (lakon(nelem)(4:6).eq.'20R')) then
        if(lakon(nelem)(7:8).eq.'LC') then
          mint3d=8*nlayer
        else
          mint3d=8
        endif
      elseif(lakon(nelem)(4:4).eq.'2') then
        mint3d=27
      elseif(lakon(nelem)(4:5).eq.'10') then
        mint3d=4
      elseif(lakon(nelem)(4:4).eq.'4') then
        mint3d=1
      elseif(lakon(nelem)(4:5).eq.'15') then
        if(lakon(nelem)(7:8).eq.'LC') then
          mint3d=6*nlayer
        else
          mint3d=9
        endif
      elseif(lakon(nelem)(4:5).eq.'6') then
        mint3d=2
      else
        mint3d=0
      endif
!     
      do jj=1,mint3d
        if(lakon(nelem)(4:5).eq.'8R') then
          xi=gauss3d1(1,jj)
          et=gauss3d1(2,jj)
          ze=gauss3d1(3,jj)
          weight=weight3d1(jj)
        elseif(lakon(nelem)(4:7).eq.'20RB') then
          if((lakon(nelem)(8:8).eq.'R').or.
     &         (lakon(nelem)(8:8).eq.'C')) then
            xi=gauss3d13(1,kk)
            et=gauss3d13(2,kk)
            ze=gauss3d13(3,kk)
            weight=weight3d13(kk)
          else
            call beamintscheme(lakon(nelem),mint3d,ielprop(nelem),
     &           prop,kk,xi,et,ze,weight)
          endif
        elseif((lakon(nelem)(4:4).eq.'8').or.
     &         (lakon(nelem)(4:6).eq.'20R'))
     &         then
          if(lakon(nelem)(7:8).ne.'LC') then
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
     &           tlayer(ki)-1.d0
            weight=weight*xlayer(ilayer,ki)/tlayer(ki)
            imat=ielmat(ilayer,nelem)
          endif
        elseif(lakon(nelem)(4:4).eq.'2') then
          xi=gauss3d3(1,jj)
          et=gauss3d3(2,jj)
          ze=gauss3d3(3,jj)
          weight=weight3d3(jj)
        elseif(lakon(nelem)(4:5).eq.'10') then
          xi=gauss3d5(1,jj)
          et=gauss3d5(2,jj)
          ze=gauss3d5(3,jj)
          weight=weight3d5(jj)
        elseif(lakon(nelem)(4:4).eq.'4') then
          xi=gauss3d4(1,jj)
          et=gauss3d4(2,jj)
          ze=gauss3d4(3,jj)
          weight=weight3d4(jj)
        elseif(lakon(nelem)(4:5).eq.'15') then
          if(lakon(nelem)(7:8).ne.'LC') then
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
     &           tlayer(ki)-1.d0
            weight=weight*xlayer(ilayer,ki)/tlayer(ki)
            imat=ielmat(ilayer,nelem)
          endif
        else
          xi=gauss3d7(1,jj)
          et=gauss3d7(2,jj)
          ze=gauss3d7(3,jj)
          weight=weight3d7(jj)
        endif
!     
        if(lakon(nelem)(1:5).eq.'C3D8R') then
          call shape8hr(xl,xsj,shp,gs,a)
        elseif(lakon(nelem)(1:5).eq.'C3D8I') then
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
!     calculating the initial temperature in the integration
!     point
!     
        t0l=0.d0
        if(ithermal(1).eq.1) then
          if(lakon(nelem)(4:5).eq.'8 ') then
            do i1=1,nope
              t0l=t0l+t0(konl(i1))/8.d0
            enddo
          elseif(lakon(nelem)(4:6).eq.'20 ') then
            call lintemp(t0,konl,nope,jj,t0l)
          elseif(lakon(nelem)(4:6).eq.'10T') then
            call linscal10(t0,konl,t0l,null,shp)
          else
            do i1=1,nope
              t0l=t0l+shp(4,i1)*t0(konl(i1))
            enddo
          endif
        endif
!     
!     calculating the density
!     
        call materialdata_rho(rhcon,nrhcon,imat,rho,
     &       t0l,ntmat_,ithermal)
!     
        csmass=csmass+weight*xsj*rho
      enddo
!     
      csmasstot=csmasstot+csmass
c     write(*,*) 'calcmass ',nelem,csmass,csmasstot
!     
      return
      end
      
