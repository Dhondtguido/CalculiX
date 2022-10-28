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
      subroutine jouleheating(ipkon,lakon,kon,co,elcon,nelcon,
     &     mi,ne,sti,ielmat,nelemload,sideload,xload,nload,nload_,
     &     iamload,nam,idefload,ncmat_,ntmat_,
     &     alcon,nalcon,ithermal,vold,t1,nmethod)
!     
!     determines the effect of Joule heating
!     
      implicit none
!     
      character*8 lakon(*)
      character*20 label,sideload(*)
!     
      integer ipkon(*),nelem,kon(*),mi(*),nope,indexe,j,k,null,
     &     mint3d,jj,iflag,ne,nelemload(2,*),iamload(2,*),nload,nload_,
     &     ielmat(mi(3),*),konl(20),idefload(*),iamplitude,isector,nam,
     &     one,nelcon(2,*),nalcon(2,*),ithermal(*),i1,ncmat_,ntmat_,
     &     imat,nmethod
!     
      real*8 co(3,*),xl(3,20),xi,et,ze,xsj,shp(4,20),weight,xload(2,*),
     &     sti(6,mi(1),*),alpha(6),heat,elcon(0:ncmat_,ntmat_,*),volume,
     &     elconloc(ncmat_),t1l,alcon(0:6,ntmat_,*),vold(0:mi(2),*),
     &     t1(*)
!     
      include "gauss.f"
!     
      data iflag /2/
!     
      null=0
      one=1
!     
      do nelem=1,ne
        if(ipkon(nelem).lt.0) cycle
!     
!     only elements belonging to the A,V-domain experience
!     Joule heating
!     
        if(int(elcon(2,1,ielmat(1,nelem))).ne.2) cycle
!     
        imat=ielmat(1,nelem)
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
        endif
!     
        do j=1,nope
          konl(j)=kon(indexe+j)
          do k=1,3
            xl(k,j)=co(k,konl(j))
          enddo
        enddo
!     
        heat=0.d0
        volume=0.d0
!     
        if(lakon(nelem)(4:5).eq.'8R') then
          mint3d=1
        elseif((lakon(nelem)(4:4).eq.'8').or.
     &         (lakon(nelem)(4:6).eq.'20R')) then
          mint3d=8
        elseif(lakon(nelem)(4:4).eq.'2') then
          mint3d=27
        elseif(lakon(nelem)(4:5).eq.'10') then
          mint3d=4
        elseif(lakon(nelem)(4:4).eq.'4') then
          mint3d=1
        elseif(lakon(nelem)(4:5).eq.'15') then
          mint3d=9
        elseif(lakon(nelem)(4:5).eq.'6') then
          mint3d=2
        endif
!     
!     loop over the integration points
!     
        do jj=1,mint3d
          if(lakon(nelem)(4:5).eq.'8R') then
            xi=gauss3d1(1,jj)
            et=gauss3d1(2,jj)
            ze=gauss3d1(3,jj)
            weight=weight3d1(jj)
          elseif((lakon(nelem)(4:4).eq.'8').or.
     &           (lakon(nelem)(4:6).eq.'20R'))
     &           then
            xi=gauss3d2(1,jj)
            et=gauss3d2(2,jj)
            ze=gauss3d2(3,jj)
            weight=weight3d2(jj)
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
            xi=gauss3d8(1,jj)
            et=gauss3d8(2,jj)
            ze=gauss3d8(3,jj)
            weight=weight3d8(jj)
          else
            xi=gauss3d7(1,jj)
            et=gauss3d7(2,jj)
            ze=gauss3d7(3,jj)
            weight=weight3d7(jj)
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
!     calculating the temperature
!     
          t1l=0.d0
          if(ithermal(1).eq.1) then
            if(lakon(nelem)(4:5).eq.'8 ') then
              do i1=1,nope
                t1l=t1l+t1(konl(i1))/8.d0
              enddo
            elseif(lakon(nelem)(4:6).eq.'20 ')then
              call linscal(t1,konl,nope,jj,t1l,one)
            elseif(lakon(nelem)(4:6).eq.'10T') then
              call linscal10(t1,konl,t1l,null,shp)
            else
              do i1=1,nope
                t1l=t1l+shp(4,i1)*t1(konl(i1))
              enddo
            endif
          elseif(ithermal(1).ge.2) then
            if(lakon(nelem)(4:5).eq.'8 ') then
              do i1=1,nope
                t1l=t1l+vold(0,konl(i1))/8.d0
              enddo
            elseif(lakon(nelem)(4:6).eq.'20 ')then
              call linscal(vold,konl,nope,jj,t1l,mi(2))
            elseif(lakon(nelem)(4:6).eq.'10T') then
              call linscal10(vold,konl,t1l,mi(2),shp)
            else
              do i1=1,nope
                t1l=t1l+shp(4,i1)*vold(0,konl(i1))
              enddo
            endif
          endif
!     
!     material data (electric conductivity and
!     magnetic permeability)
!     
          call materialdata_em(elcon,nelcon,alcon,nalcon,
     &         imat,ntmat_,t1l,elconloc,ncmat_,alpha)
!
          if(nmethod.ne.2) then
!
!           time dependent current: 
!           heat=sigma*E**2
!
            heat=heat+weight*xsj*alpha(1)*
     &           (sti(1,jj,nelem)**2+
     &            sti(2,jj,nelem)**2+
     &            sti(3,jj,nelem)**2)
          else
!
!           alernating current: 
!           heat=sigma*(E_real**2+E_imaginary**2)
!
            heat=heat+weight*xsj*alpha(1)*
     &           (sti(1,jj,nelem)**2+
     &            sti(2,jj,nelem)**2+
     &            sti(3,jj,nelem)**2+
     &            sti(1,jj,nelem+ne)**2+
     &            sti(2,jj,nelem+ne)**2+
     &            sti(3,jj,nelem+ne)**2)
          endif
          volume=volume+weight*xsj
        enddo
!     
        heat=heat/volume
!     
!     adding the Joule heating to the distributed loading
!     
        label='BF                  '    
        iamplitude=0
        isector=0
        call loadadd(nelem,label,heat,nelemload,sideload,
     &       xload,nload,nload_,iamload,iamplitude,nam,isector,
     &       idefload)
      enddo
!     
      return
      end
      
