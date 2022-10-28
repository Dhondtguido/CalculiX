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
      subroutine printoutelem(prlab,ipkon,lakon,kon,co,
     &     ener,mi,ii,nelem,energytot,volumetot,enerkintot,ne,
     &     stx,nodes,thicke,ielmat,ielem,iface,mortar,ielprop,prop,
     &     sideload,nload,nelemload,xload,bhetot,xmasstot,xinertot,
     &     cg,ithermal,rhcon,nrhcon,ntmat_,t1,vold,ipobody,ibody,
     &     xbody,nbody)
!
!     stores whole element results for element "nelem" in the .dat file
!
      implicit none
!
      character*6 prlab(*)
      character*8 lakon(*)
      character*20 sideload(*)
!
      integer ipkon(*),nelem,ii,kon(*),mi(*),nope,indexe,i,j,k,
     &  konl(20),iface,mortar,ielem,ielprop(*),nvol,nmas,nbhe,
     &  mint3d,jj,nener,iflag,nkin,ne,nodes,ki,kl,ilayer,nlayer,kk,
     &  nopes,ielmat(mi(3),*),mint2d,null,id,nload,nelemload(2,*),
     &  ithermal(*),nrhcon(*),ntmat_,imat,i1,ipobody(2,*),ibody(3,*),
     &  index,nbody
!
      real*8 ener(2,mi(1),*),energytot,volumetot,energy,volume,co(3,*),
     &  xl(3,20),xi,et,ze,xsj,shp(4,20),weight,enerkintot,enerkin,
     &  stx(6,mi(1),*),a,gs(8,4),dlayer(4),tlayer(4),thickness,
     &  thicke(mi(3),*),xlayer(mi(3),4),shp2(7,8),xs2(3,7),xsj2(3),
     &  xl2(3,8),prop(*),dflux,xload(2,*),bhe,bhetot,xmass,xmasstot,
     &  xiner(6),xinertot(6),cg(3),t1l,rho,rhcon(0:1,ntmat_,*),
     &  t1(*),vold(0:mi(2),*),dxmass,xlint(3),xbody(7,*),om
!
      include "gauss.f"
!
      data iflag /2/
!
      if(ipkon(nelem).lt.0) return
      indexe=ipkon(nelem)
      null=0
!
      nener=0
      nkin=0
      nvol=0
      nmas=0
      nbhe=0
!
      if((prlab(ii)(1:4).eq.'ELSE').or.(prlab(ii)(1:4).eq.'CELS')) then
         nener=1
      elseif(prlab(ii)(1:4).eq.'ELKE') then
         nkin=1
      elseif(prlab(ii)(1:4).eq.'EVOL') then
         nvol=1
      elseif(prlab(ii)(1:4).eq.'EMAS') then
         nmas=1
      elseif(prlab(ii)(1:4).eq.'EBHE') then
         nbhe=1
      endif
!
!     for contact displacements and stresses no integration has
!     to be performed
!
      if((prlab(ii)(1:5).eq.'CDIS ').or.
     &        (prlab(ii)(1:5).eq.'CDIST')) then
!
!        contact displacements
!
         if(mortar.eq.0) then
            write(5,'(i10,1p,1x,e13.6,1p,1x,e13.6,1p,1x,e13.6)') nodes,
     &           stx(1,1,nelem),stx(2,1,nelem),stx(3,1,nelem)
         elseif(mortar.eq.1) then
            write(5,'(i10,1x,i10,1p,1x,e13.6,1p,1x,e13.6,1p,1x,e13.6)') 
     &           ielem,iface,
     &           stx(1,1,nelem),stx(2,1,nelem),stx(3,1,nelem)
         endif
         return
      elseif((prlab(ii)(1:5).eq.'CSTR ').or.
     &        (prlab(ii)(1:5).eq.'CSTRT')) then
!
!        contact stresses
!
         if(mortar.eq.0) then
            write(5,'(i10,1p,1x,e13.6,1p,1x,e13.6,1p,1x,e13.6)') nodes,
     &           stx(4,1,nelem),stx(5,1,nelem),stx(6,1,nelem)
         elseif(mortar.eq.1) then
            write(5,'(i10,1x,i10,1p,1x,e13.6,1p,1x,e13.6,1p,1x,e13.6)') 
     &           ielem,iface,
     &           stx(4,1,nelem),stx(5,1,nelem),stx(6,1,nelem)
         endif
         return
      elseif(prlab(ii)(1:4).eq.'EBHE') then
!
!        body heating: check whether there is any body heating in
!        this element
!
         dflux=0.d0
         if(nload.gt.0) then
            call nident2(nelemload,nelem,nload,id)
            do
               if((id.eq.0).or.(nelemload(1,id).ne.nelem)) exit
               if(sideload(id)(1:2).ne.'BF') then
                  id=id-1
                  cycle
               endif
               dflux=xload(1,id)
               exit
            enddo
         endif
!
!        if no body heat flux: print and leave
!
         if(dflux.eq.0.d0) then
            if((prlab(ii)(1:5).eq.'EBHE ').or.
     &           (prlab(ii)(1:5).eq.'EBHET')) then
               write(5,'(i10,1p,1x,e13.6)') nelem,dflux
            endif
            return
         endif
      elseif(prlab(ii)(1:4).eq.'CENT') then
!
!        centrifugal loading
!
         om=0.d0
         if(nbody.gt.0) then
            index=nelem
            do
               j=ipobody(1,index)
               if(j.eq.0) exit
               if(ibody(1,j).eq.1) then
                  om=xbody(1,j)
                  exit
               endif
               index=ipobody(2,index)
               if(index.eq.0) exit
            enddo
         endif
         write(5,'(i10,1p,1x,e13.6)') nelem,om
      endif
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
!        composite materials
!
      if(lakon(nelem)(7:8).ne.'LC') then
         imat=ielmat(1,nelem)
      else
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
         if(lakon(nelem)(4:4).eq.'2') then
!
            mint2d=4
            nopes=8
!
!           determining the layer thickness and global thickness
!           at the shell integration points
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
   !
            mint2d=3
            nopes=6
   !
   !        determining the layer thickness and global thickness
   !        at the shell integration points
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
      energy=0.d0
      volume=0.d0
      bhe=0.d0
      enerkin=0.d0
      xmass=0.d0
      do j=1,6
         xiner(j)=0.d0
      enddo
!
      if(lakon(nelem)(4:5).eq.'8R') then
         mint3d=1
      elseif(lakon(nelem)(4:7).eq.'20RB') then
         if((lakon(nelem)(8:8).eq.'R').or.
     &      (lakon(nelem)(8:8).eq.'C')) then
            mint3d=50
         else
            call beamintscheme(lakon(nelem),mint3d,ielprop(nelem),prop,
     &           null,xi,et,ze,weight)
         endif
      elseif((lakon(nelem)(4:4).eq.'8').or.
     &        (lakon(nelem)(4:6).eq.'20R')) then
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
         if(nener.eq.1)then
            energy=ener(1,1,nelem)
         endif
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
               xi=gauss3d13(1,jj)
               et=gauss3d13(2,jj)
               ze=gauss3d13(3,jj)
               weight=weight3d13(jj)
            else
               call beamintscheme(lakon(nelem),mint3d,ielprop(nelem),
     &              prop,jj,xi,et,ze,weight)
            endif
         elseif((lakon(nelem)(4:4).eq.'8').or.
     &           (lakon(nelem)(4:6).eq.'20R'))
     &           then
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
     &              tlayer(ki)-1.d0
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
     &              tlayer(ki)-1.d0
               weight=weight*xlayer(ilayer,ki)/tlayer(ki)
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
         if(nener.eq.1) then
            energy=energy+weight*xsj*ener(1,jj,nelem)
         elseif(nkin.eq.1) then
            enerkin=enerkin+weight*xsj*ener(2,jj,nelem)
         elseif(nvol.eq.1) then
            volume=volume+weight*xsj
         elseif(nmas.eq.1) then
!
!           coordinates of the integration point
!
            do k=1,3
               xlint(k)=0.d0
            enddo
            do j=1,nope
               do k=1,3
                  xlint(k)=xlint(k)+xl(k,j)*shp(4,j)
               enddo
            enddo
!
!           temperature of the integration point
!
            t1l=0.d0
            if(ithermal(1).eq.1) then
               if((lakon(nelem)(4:5).eq.'8 ').or.
     &              (lakon(nelem)(4:5).eq.'8I')) then
                  do i1=1,8
                     t1l=t1l+t1(konl(i1))/8.d0
                  enddo
               elseif(lakon(nelem)(4:6).eq.'20 ') then
                  call lintemp(t1,konl,nope,jj,t1l)
               elseif(lakon(nelem)(4:6).eq.'10T') then
                  call linscal10(t1,konl,t1l,null,shp)
               else
                  do i1=1,nope
                     t1l=t1l+shp(4,i1)*t1(konl(i1))
                  enddo
               endif
            elseif(ithermal(1).ge.2) then
               if((lakon(nelem)(4:5).eq.'8 ').or.
     &              (lakon(nelem)(4:5).eq.'8I')) then
                  do i1=1,8
                     t1l=t1l+vold(0,konl(i1))/8.d0
                  enddo
               elseif(lakon(nelem)(4:6).eq.'20 ') then
                  call lintemp_th1(vold,konl,nope,jj,t1l,mi)
               elseif(lakon(nelem)(4:6).eq.'10T') then
                  call linscal10(vold,konl,t1l,mi(2),shp)
               else
                  do i1=1,nope
                     t1l=t1l+shp(4,i1)*vold(0,konl(i1))
                  enddo
               endif
            endif
!
!           calculating the density
!
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           t1l,ntmat_,ithermal(1))
!
            dxmass=weight*xsj*rho
            xmass=xmass+dxmass
            xiner(1)=xiner(1)+xlint(1)*xlint(1)*dxmass
            xiner(2)=xiner(2)+xlint(2)*xlint(2)*dxmass
            xiner(3)=xiner(3)+xlint(3)*xlint(3)*dxmass
            xiner(4)=xiner(4)+xlint(1)*xlint(2)*dxmass
            xiner(5)=xiner(5)+xlint(1)*xlint(3)*dxmass
            xiner(6)=xiner(6)+xlint(2)*xlint(3)*dxmass
            do i1=1,3
               cg(i1)=cg(i1)+xlint(i1)*dxmass
            enddo
         elseif(nbhe.eq.1) then
            bhe=bhe+dflux*weight*xsj
         endif
      enddo
!
      if(nener.eq.1) then
         energytot=energytot+energy
      elseif(nkin.eq.1) then
         enerkintot=enerkintot+enerkin
      elseif(nvol.eq.1) then
         volumetot=volumetot+volume
      elseif(nmas.eq.1) then
         xmasstot=xmasstot+xmass
         do j=1,6
            xinertot(j)=xinertot(j)+xiner(j)
         enddo
      elseif(nbhe.eq.1) then
         bhetot=bhetot+bhe
      endif
!     
!     writing to file
!     
      if((prlab(ii)(1:5).eq.'ELSE ').or.
     &     (prlab(ii)(1:5).eq.'ELSET')) then
         write(5,'(i10,1p,1x,e13.6)') nelem,energy
      elseif((prlab(ii)(1:5).eq.'CELS ').or.
     &        (prlab(ii)(1:5).eq.'CELST')) then
         if(mortar.eq.0) then
            write(5,'(i10,1p,1x,e13.6)') nodes,energy
         elseif(mortar.eq.1) then
            write(5,'(i10,1x,i10,1p,1x,e13.6)') ielem,iface,energy
         endif
      elseif((prlab(ii)(1:5).eq.'EVOL ').or.
     &        (prlab(ii)(1:5).eq.'EVOLT')) then
         write(5,'(i10,1p,1x,e13.6)') nelem,volume
      elseif((prlab(ii)(1:5).eq.'EMAS ').or.
     &        (prlab(ii)(1:5).eq.'EMAST')) then
         write(5,'(i10,1p,7(1x,e13.6))') nelem,xmass,
     &     (xiner(i),i=1,6)
      elseif((prlab(ii)(1:5).eq.'ELKE ').or.
     &        (prlab(ii)(1:5).eq.'ELKET')) then
         write(5,'(i10,1p,1x,e13.6)') nelem,enerkin
      elseif((prlab(ii)(1:5).eq.'EBHE ').or.
     &        (prlab(ii)(1:5).eq.'EBHET')) then
         write(5,'(i10,1p,1x,e13.6)') nelem,bhe
      endif
!
      return
      end
      
